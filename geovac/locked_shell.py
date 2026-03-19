"""
LockedShellMolecule — Closed-shell-aware molecular solver
==========================================================

Closed shells and subshells are single quantum states (Paper 16 Type B).
Don't expand them into determinants. Lock them as single states with
precomputed energies and couple to the active space via J/K integrals.

Molecular wavefunction = tensor product of locked units, NOT determinant
expansion over all orbitals.

Example: LiH
  - Li 1s^2: LOCKED (1 state, E_core precomputed)
  - Bond pair in {Li 2s, H 1s}: ACTIVE (C(4,2) = 6 SDs)
  - Total: 6 determinants, not 367,290

Scaling: O(n_active^4) where n_active ~ 4-20 spin-orbitals.

Author: GeoVac Development Team
Date: March 2026
"""

import time
from itertools import combinations
from math import comb
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


class LockedShellMolecule:
    """
    Molecular solver that locks closed shells as single states.

    Parameters
    ----------
    Z_A, Z_B : int
        Nuclear charges.
    nmax_A, nmax_B : int
        Basis truncation per atom.
    R : float
        Internuclear distance (bohr).
    n_electrons : int
        Total electron count.
    locked_config : dict, optional
        Explicit locked configuration. Maps atom index (0 or 1) to list
        of (n, l) subshells to lock. Each locked subshell must be fully
        occupied (2*(2l+1) electrons).
        Example: {0: [(1, 0)]} locks Li 1s on atom A.
        If None, auto-detects closed shells.
    active_nmax : int, optional
        Maximum n for active orbitals (default: 2). Higher values include
        more virtual orbitals for correlation recovery.
    vee_method : str
        V_ee method (default: 'slater_full').

    Example
    -------
    >>> mol = LockedShellMolecule(Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
    ...     R=3.015, n_electrons=4, locked_config={0: [(1, 0)]})
    >>> E, psi = mol.solve()
    """

    def __init__(
        self,
        Z_A: int,
        Z_B: int,
        nmax_A: int,
        nmax_B: int,
        R: float,
        n_electrons: int,
        locked_config: Optional[Dict[int, List[Tuple[int, int]]]] = None,
        active_nmax: int = 2,
        vee_method: str = 'slater_full',
    ) -> None:
        from .lattice_index import MolecularLatticeIndex

        self.Z_A = Z_A
        self.Z_B = Z_B
        self.R = R
        self.n_electrons = n_electrons
        self.active_nmax = active_nmax

        t0 = time.perf_counter()

        # Build molecular index (skip SD enumeration — we build our own)
        self._parent = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons,
            vee_method=vee_method,
            enumerate_sds=False,
        )

        self.V_NN = self._parent.V_NN
        self.n_spatial = self._parent._n_spatial
        self.n_spatial_A = self._parent._n_spatial_A

        # Identify locked and active orbitals
        if locked_config is None:
            locked_config = self._auto_detect_locked()
        self._locked_config = locked_config

        self._classify_orbitals(locked_config)

        # Densify integrals
        self._h1_diag = np.array(self._parent._h1_diag)
        self._H1_dense = np.asarray(self._parent._H1_spatial.todense())
        self._eri_4d = self._build_eri_dense()

        # Compute locked-shell energy and effective H1
        self.E_locked = self._compute_locked_energy()
        self._h1_eff = self._build_effective_h1()

        # Enumerate active-space SDs
        self._enumerate_active_sds()

        elapsed = time.perf_counter() - t0

        full_sd_est = comb(self._parent.n_sp, n_electrons)
        print(
            f"[LockedShell] {len(self._locked_spatial)} locked spatial orbitals, "
            f"{len(self._active_spatial)} active spatial ({self.n_active_sp} spin-orbs), "
            f"{self.n_sd} active SDs "
            f"(vs ~{full_sd_est:,} full, ~{full_sd_est // max(self.n_sd, 1):,}x reduction), "
            f"E_locked = {self.E_locked:.6f} Ha, "
            f"setup = {elapsed:.3f}s"
        )

    # ------------------------------------------------------------------
    # Orbital classification
    # ------------------------------------------------------------------

    def _auto_detect_locked(self) -> Dict[int, List[Tuple[int, int]]]:
        """
        Auto-detect closed shells based on electron configuration.

        For each atom, identify (n, l) subshells that would be fully
        occupied in the ground-state Aufbau configuration and lock them
        if they are deep core orbitals (n < valence n).
        """
        config: Dict[int, List[Tuple[int, int]]] = {}

        for atom_idx, Z in enumerate([self.Z_A, self.Z_B]):
            if Z <= 1:
                continue  # H has no core

            # Aufbau filling: find the highest occupied shell
            electrons_remaining = Z
            shells_to_lock = []
            for n in range(1, 10):
                for l in range(n):
                    capacity = 2 * (2 * l + 1)
                    if electrons_remaining >= capacity:
                        electrons_remaining -= capacity
                        shells_to_lock.append((n, l))
                    else:
                        break
                if electrons_remaining == 0:
                    break

            # Lock all but the outermost occupied shell
            if len(shells_to_lock) > 1:
                # Keep the last shell unlocked (it's the valence)
                config[atom_idx] = shells_to_lock[:-1]
            elif len(shells_to_lock) == 1 and electrons_remaining > 0:
                # Single complete inner shell + partial outer: lock the inner
                config[atom_idx] = shells_to_lock

        return config

    def _classify_orbitals(
        self, locked_config: Dict[int, List[Tuple[int, int]]]
    ) -> None:
        """
        Partition spatial orbitals into locked and active sets.

        Locked: all spatial orbitals matching (n, l) in locked_config.
        Active: spatial orbitals with n <= active_nmax that aren't locked.
        """
        sp_states = self._parent.sp_states
        spatial_atom = self._parent._spatial_atom

        self._locked_spatial: List[int] = []
        self._active_spatial: List[int] = []
        self._locked_spinorb: List[int] = []
        self._active_spinorb: List[int] = []

        n_locked_electrons = 0

        for sp_idx in range(self.n_spatial):
            # sp_states has pairs: (n,l,m,0), (n,l,m,1) for each spatial
            n, l, m, _ = sp_states[sp_idx * 2]
            atom = spatial_atom[sp_idx]

            # Check if this orbital is locked
            is_locked = False
            if atom in locked_config:
                for locked_n, locked_l in locked_config[atom]:
                    if n == locked_n and l == locked_l:
                        is_locked = True
                        break

            if is_locked:
                self._locked_spatial.append(sp_idx)
                self._locked_spinorb.extend([sp_idx * 2, sp_idx * 2 + 1])
                n_locked_electrons += 2  # doubly occupied
            elif n <= self.active_nmax:
                self._active_spatial.append(sp_idx)
                self._active_spinorb.extend([sp_idx * 2, sp_idx * 2 + 1])

        self._locked_spinorb_set = frozenset(self._locked_spinorb)
        self._active_spinorb_set = frozenset(self._active_spinorb)
        self.n_locked_el = n_locked_electrons
        self.n_active_el = self.n_electrons - n_locked_electrons
        self.n_active_sp = len(self._active_spinorb)

        if self.n_active_el < 0:
            raise ValueError(
                f"Locked {n_locked_electrons} electrons but molecule has "
                f"only {self.n_electrons}"
            )

        locked_desc = []
        for atom_idx, shells in sorted(locked_config.items()):
            atom_name = "A" if atom_idx == 0 else "B"
            for n, l in shells:
                l_name = "spdfg"[l]
                locked_desc.append(f"{atom_name}:{n}{l_name}")
        print(
            f"[LockedShell] Locked: [{', '.join(locked_desc)}] "
            f"({n_locked_electrons}e), "
            f"Active: {self.n_active_el}e in {len(self._active_spatial)} "
            f"spatial orbs (n <= {self.active_nmax})"
        )

    # ------------------------------------------------------------------
    # Integral construction
    # ------------------------------------------------------------------

    def _build_eri_dense(self) -> np.ndarray:
        """Build dense 4D ERI array from parent's sparse dict."""
        n = self.n_spatial
        eri = np.zeros((n, n, n, n))
        for (a, b, c, d), val in self._parent._eri.items():
            eri[a, b, c, d] = val
        return eri

    def _compute_locked_energy(self) -> float:
        """
        Compute energy of locked shells: h1 + J - K.

        E_locked = sum_c h1[c,c] + sum_{c<c'} [J(c,c') - K(c,c')]
        """
        eri = self._eri_4d
        e_h1 = 0.0
        for c in self._locked_spinorb:
            e_h1 += self._h1_diag[c >> 1]

        e_jk = 0.0
        locked = self._locked_spinorb
        for i in range(len(locked)):
            ci = locked[i]
            sp_i = ci >> 1
            sig_i = ci & 1
            for j in range(i + 1, len(locked)):
                cj = locked[j]
                sp_j = cj >> 1
                e_jk += eri[sp_i, sp_j, sp_i, sp_j]
                if (cj & 1) == sig_i:
                    e_jk -= eri[sp_i, sp_j, sp_j, sp_i]

        return e_h1 + e_jk

    def _build_effective_h1(self) -> np.ndarray:
        """
        Build effective one-electron integrals: h1 + locked-shell J/K.

        h_eff[a,b] = h1[a,b] + sum_{gamma in locked_spatial}
                     [2*eri(a,gamma,b,gamma) - eri(a,gamma,gamma,b)]
        """
        eri = self._eri_4d
        h_eff = self._H1_dense.copy()

        for gamma in self._locked_spatial:
            for a in range(self.n_spatial):
                for b in range(self.n_spatial):
                    h_eff[a, b] += (
                        2.0 * eri[a, gamma, b, gamma]
                        - eri[a, gamma, gamma, b]
                    )

        return h_eff

    # ------------------------------------------------------------------
    # Active-space SD enumeration
    # ------------------------------------------------------------------

    def _enumerate_active_sds(self) -> None:
        """Enumerate SDs over active spin-orbitals only."""
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(self._active_spinorb, self.n_active_el)
        )
        self.n_sd = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }

    # ------------------------------------------------------------------
    # Hamiltonian assembly
    # ------------------------------------------------------------------

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Build active-space Hamiltonian using effective integrals.

        Uses h_eff (h1 + locked J/K) for one-electron terms and
        original ERIs for active-active two-electron terms.
        """
        t0 = time.perf_counter()

        sd_basis = self.sd_basis
        sd_index = self._sd_index
        n_el = self.n_active_el
        n_sd = self.n_sd
        threshold = self._parent.threshold
        h_eff = self._h1_eff
        eri = self._eri_4d
        active_sp = self._active_spinorb

        h_eff_diag = np.array([h_eff[i, i] for i in range(self.n_spatial)])

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        occ_sets: List[frozenset] = [frozenset(sd) for sd in sd_basis]

        # --- Diagonal ---
        for I, sd_I in enumerate(sd_basis):
            h_diag = 0.0
            for p in sd_I:
                h_diag += h_eff_diag[p >> 1]

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

        # --- Singles ---
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = occ_sets[I]

            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for r in active_sp:
                    if r in occ_I:
                        continue
                    if (r & 1) != sig_p:
                        continue

                    sp_r = r >> 1
                    me = h_eff[sp_p, sp_r]
                    for q in sd_I:
                        if q == p:
                            continue
                        sp_q = q >> 1
                        me += eri[sp_p, sp_q, sp_r, sp_q]
                        if (q & 1) == sig_p:
                            me -= eri[sp_p, sp_q, sp_q, sp_r]

                    if abs(me) < threshold:
                        continue

                    new_sd = list(sd_I)
                    new_sd[kp] = r
                    new_sd_t = tuple(sorted(new_sd))
                    J = sd_index.get(new_sd_t)
                    if J is None or J <= I:
                        continue

                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        # --- Doubles ---
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

                    for ir in range(len(active_sp)):
                        r = active_sp[ir]
                        if r in occ_I:
                            continue
                        sig_r = r & 1

                        for js in range(ir + 1, len(active_sp)):
                            s = active_sp[js]
                            if s in occ_I:
                                continue
                            sig_s = s & 1

                            if sig_p + sig_q != sig_r + sig_s:
                                continue

                            sp_r = r >> 1
                            sp_s = s >> 1

                            me = 0.0
                            if sig_p == sig_r and sig_q == sig_s:
                                me += eri[sp_p, sp_q, sp_r, sp_s]
                            if sig_p == sig_s and sig_q == sig_r:
                                me -= eri[sp_p, sp_q, sp_s, sp_r]

                            if abs(me) < threshold:
                                continue

                            new_sd = list(sd_I)
                            new_sd[kp] = r
                            new_sd[kq] = s
                            new_sd_t = tuple(sorted(new_sd))
                            J_idx = sd_index.get(new_sd_t)
                            if J_idx is None or J_idx <= I:
                                continue

                            phase = self._compute_double_phase(
                                sd_I, kp, kq, r, s
                            )
                            off_rows.append(I)
                            off_cols.append(J_idx)
                            off_vals.append(phase * me)

        # Assemble
        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)), shape=(n_sd, n_sd)
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)), shape=(n_sd, n_sd)
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[LockedShell] H assembled: shape={H.shape}, nnz={H.nnz:,}, "
            f"time={elapsed:.3f}s"
        )
        return H.tocsr()

    # ------------------------------------------------------------------
    # Phase computation
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_phase(sd: Tuple[int, ...], kp: int, r: int) -> float:
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _compute_double_phase(
        sd: Tuple[int, ...], kp: int, kq: int, r: int, s: int
    ) -> float:
        p = sd[kp]
        q = sd[kq]
        n_swap_p = kp
        remaining_1 = [a for a in sd if a != p]
        kr = sum(1 for a in remaining_1 if a < r)
        phase1 = (-1) ** (n_swap_p + kr)
        intermediate = sorted(remaining_1[:kr] + [r] + remaining_1[kr:])
        kq_new = intermediate.index(q)
        remaining_2 = [a for a in intermediate if a != q]
        ks = sum(1 for a in remaining_2 if a < s)
        phase2 = (-1) ** (kq_new + ks)
        return float(phase1 * phase2)

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def solve(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve active-space problem with locked shells.

        Returns
        -------
        eigvals : np.ndarray
            Total energies: E_locked + E_active + V_NN.
        eigvecs : np.ndarray
            Active-space CI vectors.
        """
        H = self.assemble_hamiltonian()

        if self.n_sd <= 2:
            # Dense diagonalization for tiny matrices
            H_dense = H.toarray()
            eigvals, eigvecs = np.linalg.eigh(H_dense)
            eigvals = eigvals[:n_states]
            eigvecs = eigvecs[:, :n_states]
        else:
            k = min(n_states, self.n_sd - 2)
            rng = np.random.RandomState(42)
            v0 = rng.randn(H.shape[0])
            eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
            order = np.argsort(eigvals)
            eigvals, eigvecs = eigvals[order], eigvecs[:, order]

        # Add locked energy + nuclear repulsion
        eigvals = eigvals + self.E_locked + self.V_NN

        E_active = eigvals[0] - self.E_locked - self.V_NN
        print(
            f"[LockedShell] E_locked = {self.E_locked:.6f}, "
            f"E_active = {E_active:.6f}, "
            f"V_NN = {self.V_NN:.6f}, "
            f"E_total = {eigvals[0]:.6f} Ha"
        )
        return eigvals, eigvecs
