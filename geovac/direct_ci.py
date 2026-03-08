"""
DirectCISolver — Excitation-driven FCI Hamiltonian construction
================================================================

Replaces the O(N^2_SD) pairwise loop in LatticeIndex._assemble_hamiltonian_full
with excitation-driven sparse matrix construction in O(N_SD x N_connected).

Algorithm overview:
  For each determinant I, generate all singly and doubly excited determinants J
  connected by the Hamiltonian. Build the sparse H matrix in COO format, then
  convert to CSR for eigsh.

  Singles: O(N_SD x n_el x n_sp) — iterate occupied p, all virtual r with same
    spin, compute Slater-Condon single-excitation matrix element.

  Doubles: O(N_SD x n_occ_pairs x |eri_targets|) — for each occupied pair (p,q),
    iterate over spatial orbital pairs (r,s) that have non-zero ERI entries.
    Precomputed spatial_targets dict avoids the O(n_virt^2) brute-force search.

Performance optimisation (v0.9.7):
  All matrix element computation uses dense NumPy arrays instead of sparse matrix
  element access and dict lookups. The one-electron Hamiltonian H1_spatial is
  densified (n_spatial x n_spatial), and ERIs are stored in a dense 4D array
  (n_spatial^4). This eliminates scipy's slow _validate_indices path which
  dominated >80% of singles wall time.

Reference: Knowles & Handy, Chem. Phys. Lett. 111, 315 (1984).

Author: GeoVac Development Team
Date: March 2026
"""

import time
from typing import Dict, List, Set, Tuple

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


class DirectCISolver:
    """
    Excitation-driven FCI Hamiltonian construction and diagonalization.

    Borrows all data structures from a pre-built LatticeIndex: sd_basis,
    _sd_index, _eri, _h1_diag, _H1_spatial, _slater_condon_* methods.

    Parameters
    ----------
    lattice_index : LatticeIndex
        Fully constructed LatticeIndex with ERI and SD basis already built.
    """

    def __init__(self, lattice_index: 'LatticeIndex') -> None:
        self._idx = lattice_index
        self.n_sd: int = lattice_index.n_sd
        self.n_sp: int = lattice_index.n_sp
        self.n_el: int = lattice_index.n_electrons
        self.n_spatial: int = lattice_index.n_sp // 2

        t0 = time.perf_counter()

        # Precompute occupied sets for each SD (O(N_SD x n_el))
        self._occ_sets: List[frozenset] = [
            frozenset(sd) for sd in lattice_index.sd_basis
        ]

        # Precompute spatial ERI targets for sparse double excitations
        self._spatial_targets: Dict[
            Tuple[int, int], Set[Tuple[int, int]]
        ] = self._build_spatial_targets()

        # Dense arrays for fast element access (replaces sparse/dict lookups)
        self._H1_dense: np.ndarray = self._precompute_H1_dense()
        self._eri_4d: np.ndarray = self._precompute_eri_dense()
        self._h1_diag_arr: np.ndarray = np.array(lattice_index._h1_diag)

        n_targets = sum(len(v) for v in self._spatial_targets.values())
        dt = time.perf_counter() - t0
        print(
            f"[DirectCI] Precomputed: {dt:.3f}s "
            f"(N_SD={self.n_sd:,}, spatial_target_pairs={n_targets})"
        )

    # ------------------------------------------------------------------
    # Precomputation methods
    # ------------------------------------------------------------------

    def _precompute_H1_dense(self) -> np.ndarray:
        """
        Convert sparse H1_spatial to dense array for fast element access.

        The sparse matrix element access path (_validate_indices) in scipy
        costs ~24 us per access. Dense array access costs ~100 ns.
        For 118k calls this saves ~2.8s.

        Returns
        -------
        H1 : np.ndarray, shape (n_spatial, n_spatial)
        """
        return np.asarray(self._idx._H1_spatial.todense())

    def _precompute_eri_dense(self) -> np.ndarray:
        """
        Build dense 4D ERI array from sparse dict.

        For nmax=4 (30 spatial orbitals): 30^4 = 810k entries = 6.5 MB.
        For nmax=5 (55 spatial orbitals): 55^4 = 9.15M entries = 73 MB.
        Both are acceptable for dense storage.

        Returns
        -------
        eri : np.ndarray, shape (n_spatial, n_spatial, n_spatial, n_spatial)
        """
        n = self.n_spatial
        eri = np.zeros((n, n, n, n))
        for (a, b, c, d), val in self._idx._eri.items():
            eri[a, b, c, d] = val
        return eri

    def _build_spatial_targets(
        self,
    ) -> Dict[Tuple[int, int], Set[Tuple[int, int]]]:
        """
        Build sparse ERI target lookup for double excitations.

        For each spatial pair (a, b) appearing as the first two indices of a
        non-zero ERI entry, collects all canonical spatial pairs (lo, hi) that
        appear as the last two indices. Also ensures the reverse symmetry
        <ab|cd> = <cd|ab> so excitations can be discovered from either direction.

        Returns
        -------
        targets : dict
            Maps (sp_a, sp_b) -> set of (sp_lo, sp_hi) with sp_lo <= sp_hi.
        """
        targets: Dict[Tuple[int, int], Set[Tuple[int, int]]] = {}

        if not hasattr(self._idx, '_eri'):
            return targets

        threshold = self._idx.threshold
        for (a, b, c, d), val in self._idx._eri.items():
            if abs(val) < threshold:
                continue
            pair_cd = (min(c, d), max(c, d))
            pair_ab = (min(a, b), max(a, b))

            # Forward: occupied (a,b) → virtual (c,d)
            targets.setdefault((a, b), set()).add(pair_cd)
            # Reverse symmetry: occupied (c,d) → virtual (a,b)
            targets.setdefault((c, d), set()).add(pair_ab)
            targets.setdefault((d, c), set()).add(pair_ab)

        return targets

    # ------------------------------------------------------------------
    # Hamiltonian assembly
    # ------------------------------------------------------------------

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Build the FCI Hamiltonian via excitation-driven iteration.

        Complexity: O(N_SD x n_el x n_sp) for singles
                  + O(N_SD x C(n_el,2) x |eri_targets|) for doubles
        versus O(N^2_SD) for the pairwise loop.

        Returns
        -------
        H : csr_matrix, shape (n_sd, n_sd)
        """
        t0 = time.perf_counter()

        idx = self._idx
        sd_basis = idx.sd_basis
        sd_index = idx._sd_index
        n_el = self.n_el
        n_sp = self.n_sp
        n_sd = self.n_sd
        occ_sets = self._occ_sets
        threshold = idx.threshold
        compute_phase = idx._compute_phase
        compute_double_phase = idx._compute_double_phase
        spatial_targets = self._spatial_targets

        # Dense arrays — eliminates sparse matrix and dict lookup overhead
        H1 = self._H1_dense
        eri = self._eri_4d
        h1_diag = self._h1_diag_arr

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        # --- Phase 1: Diagonal --- O(N_SD x n_el^2)
        # Inlined _slater_condon_diagonal using dense arrays
        for I, sd_I in enumerate(sd_basis):
            h_diag = 0.0
            for p in sd_I:
                h_diag += h1_diag[p >> 1]

            n = len(sd_I)
            for i in range(n):
                pi = sd_I[i]
                sp_i = pi >> 1
                sig_i = pi & 1
                for j in range(i + 1, n):
                    pj = sd_I[j]
                    sp_j = pj >> 1
                    h_diag += eri[sp_i, sp_j, sp_i, sp_j]
                    if (pj & 1) == sig_i:
                        h_diag -= eri[sp_i, sp_j, sp_j, sp_i]

            if abs(h_diag) >= threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # --- Phase 2: Singles --- O(N_SD x n_el x n_sp)
        # Inlined Slater-Condon single excitation using dense arrays.
        # Matrix element: h1[p,r] + Σ_{q∈occ, q≠p} [eri(p,q,r,q) - δ(σ)·eri(p,q,q,r)]
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = occ_sets[I]

            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for r in range(n_sp):
                    if r in occ_I:
                        continue
                    if (r & 1) != sig_p:
                        continue  # spin conservation

                    sp_r = r >> 1

                    # Inline Slater-Condon single-excitation matrix element
                    me = H1[sp_p, sp_r]
                    for q in sd_I:
                        if q == p:
                            continue
                        sp_q = q >> 1
                        me += eri[sp_p, sp_q, sp_r, sp_q]
                        if (q & 1) == sig_p:
                            me -= eri[sp_p, sp_q, sp_q, sp_r]

                    if abs(me) < threshold:
                        continue

                    # Build target SD
                    new_sd = list(sd_I)
                    new_sd[kp] = r
                    new_sd_t = tuple(sorted(new_sd))
                    J = sd_index.get(new_sd_t)
                    if J is None or J <= I:
                        continue  # upper triangle only

                    phase = compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        t_singles = time.perf_counter() - t0 - t_diag

        # --- Phase 3: Doubles --- O(N_SD x C(n_el,2) x |targets|)
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = occ_sets[I]

            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for kq in range(kp + 1, n_el):
                    q = sd_I[kq]
                    sp_q = q >> 1
                    sig_q = q & 1

                    # Collect candidate spatial target pairs from ERI
                    cands = spatial_targets.get((sp_p, sp_q), set())
                    cands_rev = spatial_targets.get((sp_q, sp_p), set())
                    if cands_rev:
                        cands = cands | cands_rev

                    if not cands:
                        continue

                    for sp_lo, sp_hi in cands:
                        # Enumerate spin-orbital pairs (r, s) with r < s
                        if sp_lo < sp_hi:
                            spin_pairs = [
                                ((sp_lo << 1) | sr, (sp_hi << 1) | ss)
                                for sr in (0, 1) for ss in (0, 1)
                            ]
                        else:
                            # sp_lo == sp_hi: only up/down pair
                            spin_pairs = [(sp_lo << 1, (sp_lo << 1) | 1)]

                        for r, s in spin_pairs:
                            if r in occ_I or s in occ_I:
                                continue
                            sig_r = r & 1
                            sig_s = s & 1

                            # Spin conservation
                            if sig_p + sig_q != sig_r + sig_s:
                                continue

                            sp_r = r >> 1
                            sp_s = s >> 1

                            # Matrix element (dense array access)
                            me = 0.0
                            if sig_p == sig_r and sig_q == sig_s:
                                me += eri[sp_p, sp_q, sp_r, sp_s]
                            if sig_p == sig_s and sig_q == sig_r:
                                me -= eri[sp_p, sp_q, sp_s, sp_r]

                            if abs(me) < threshold:
                                continue

                            # Build target SD
                            new_sd = list(sd_I)
                            new_sd[kp] = r
                            new_sd[kq] = s
                            new_sd_t = tuple(sorted(new_sd))
                            J = sd_index.get(new_sd_t)
                            if J is None or J <= I:
                                continue

                            phase = compute_double_phase(
                                sd_I, kp, kq, r, s
                            )
                            off_rows.append(I)
                            off_cols.append(J)
                            off_vals.append(phase * me)

        t_doubles = time.perf_counter() - t0 - t_diag - t_singles

        # --- Assemble sparse matrix ---
        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)),
            shape=(n_sd, n_sd),
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)),
            shape=(n_sd, n_sd),
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[DirectCI] Hamiltonian assembled: shape={H.shape}, "
            f"nnz={H.nnz:,}, time={elapsed:.3f}s "
            f"(diag={t_diag:.1f}s, singles={t_singles:.1f}s, "
            f"doubles={t_doubles:.1f}s)"
        )
        return H.tocsr()

    def solve(
        self,
        n_states: int = 1,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Assemble Hamiltonian and compute lowest eigenvalues.

        Parameters
        ----------
        n_states : int
            Number of lowest eigenvalues to compute.

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        H = self.assemble_hamiltonian()
        k = min(n_states, self.n_sd - 2)
        rng = np.random.RandomState(42)
        v0 = rng.randn(H.shape[0])
        eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
        order = np.argsort(eigvals)
        return eigvals[order], eigvecs[:, order]
