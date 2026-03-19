"""
FrozenCoreLatticeIndex — Hierarchical molecular solver with frozen cores
========================================================================

Key insight from Paper 16: Type C [N-1,1] structure factorizes.
For molecules: S_N on S^(3N-1) → S_{N_core} × S_{N_active} on
S^(3N_core-1) × S^(3N_active-1).

This class:
1. Partitions spin-orbitals into core (frozen, always doubly occupied)
   and active (variational)
2. Precomputes core energy E_core (one-time cost)
3. Builds effective Hamiltonian with core-active J/K integrals folded
   into the one-electron operator
4. Solves only the active-space FCI (dramatically smaller)

Scaling: C(n_active_sp, n_active_el) instead of C(n_sp, n_el)
  LiH nmax=3: 367,290 SDs → 1,431 SDs (257× reduction)
  LiH nmax=5: 64,684,950 SDs → 19,503 SDs (3,317× reduction)

Reference: frozen-core CI is standard quantum chemistry
(Szabo & Ostlund, Chapter 6).

Author: GeoVac Development Team
Date: March 2026
"""

import time
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


class FrozenCoreLatticeIndex:
    """
    Frozen-core FCI solver for hierarchical molecules.

    Wraps a pre-built LatticeIndex or MolecularLatticeIndex, partitions
    the spin-orbital basis into frozen core and active space, and solves
    the active-space FCI with core-active Coulomb (J) and exchange (K)
    integrals folded into an effective one-electron Hamiltonian.

    Parameters
    ----------
    parent : LatticeIndex or MolecularLatticeIndex
        Fully constructed parent index with ERI and H1 already built.
        Must have vee_method='slater_full' for accurate integrals.
    frozen_orbitals : list of int
        Spatial orbital indices to freeze (each doubly occupied).
        Example: [0] freezes the 1s orbital (spin-orbitals 0 and 1).
    n_active_electrons : int
        Number of electrons in the active space.

    Attributes
    ----------
    E_core : float
        Total energy of the frozen core (h1 + J - K).
    n_active_sp : int
        Number of active spin-orbitals.
    n_sd : int
        Number of active-space Slater determinants.
    sd_basis : list of tuple
        Active-space SDs (using original spin-orbital indices).

    Example
    -------
    >>> from geovac import MolecularLatticeIndex
    >>> mol = MolecularLatticeIndex(Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
    ...                              R=3.015, n_electrons=4,
    ...                              vee_method='slater_full')
    >>> fc = FrozenCoreLatticeIndex(mol, frozen_orbitals=[0], n_active_electrons=2)
    >>> E, psi = fc.solve()
    """

    def __init__(
        self,
        parent: 'LatticeIndex',
        frozen_orbitals: List[int],
        n_active_electrons: int,
    ) -> None:
        self.parent = parent
        self.n_active_electrons = n_active_electrons
        self.threshold = parent.threshold

        # Frozen spatial orbitals → frozen spin-orbitals (alpha + beta)
        self.frozen_spatial: List[int] = sorted(frozen_orbitals)
        self.frozen_spinorb: List[int] = []
        for sp in self.frozen_spatial:
            self.frozen_spinorb.extend([sp << 1, (sp << 1) | 1])
        self.frozen_spinorb_set: frozenset = frozenset(self.frozen_spinorb)
        self.n_frozen_el: int = len(self.frozen_spinorb)

        # Active spin-orbitals (everything not frozen)
        self.active_spinorb: List[int] = [
            i for i in range(parent.n_sp) if i not in self.frozen_spinorb_set
        ]
        self.n_active_sp: int = len(self.active_spinorb)
        self.active_spinorb_set: frozenset = frozenset(self.active_spinorb)

        # Validate electron count
        total_electrons = self.n_frozen_el + n_active_electrons
        if total_electrons != parent.n_electrons:
            raise ValueError(
                f"Frozen ({self.n_frozen_el}) + active ({n_active_electrons}) = "
                f"{total_electrons} ≠ parent n_electrons ({parent.n_electrons})"
            )

        t0 = time.perf_counter()

        # 1. Densify parent integrals for fast access
        self._h1_diag = np.array(parent._h1_diag)
        self._H1_dense = np.asarray(parent._H1_spatial.todense())
        self._eri_4d = self._build_eri_dense()

        # 2. Compute frozen-core energy
        self.E_core = self._compute_core_energy()

        # 3. Build effective one-electron integrals (h1 + core J/K)
        self._h1_eff = self._build_effective_h1()

        # 4. Enumerate active-space SDs
        self._enumerate_active_sd_basis()

        elapsed = time.perf_counter() - t0
        if parent.n_sd > 0:
            reduction_str = (
                f"(vs {parent.n_sd:,} full, "
                f"{parent.n_sd / max(self.n_sd, 1):.0f}x reduction), "
            )
        else:
            # Parent skipped SD enumeration (too large)
            from math import comb as _comb
            n_sd_full_est = _comb(parent.n_sp, parent.n_electrons)
            reduction_str = (
                f"(vs ~{n_sd_full_est:,.0f} full, "
                f"~{n_sd_full_est / max(self.n_sd, 1):,.0f}x reduction), "
            )
        print(
            f"[FrozenCore] {len(self.frozen_spatial)} frozen spatial orbitals, "
            f"{self.n_active_sp} active spin-orbitals, "
            f"{self.n_sd:,} active SDs "
            f"{reduction_str}"
            f"E_core = {self.E_core:.6f} Ha, "
            f"setup = {elapsed:.3f}s"
        )

    # ------------------------------------------------------------------
    # Construction methods
    # ------------------------------------------------------------------

    def _build_eri_dense(self) -> np.ndarray:
        """Build dense 4D ERI array from parent's sparse dict."""
        n = self.parent.n_sp // 2
        eri = np.zeros((n, n, n, n))
        for (a, b, c, d), val in self.parent._eri.items():
            eri[a, b, c, d] = val
        return eri

    def _compute_core_energy(self) -> float:
        """
        Compute frozen-core energy via Slater-Condon diagonal formula.

        E_core = Σ_{c∈core} h1(c,c) + Σ_{c<c'∈core} [J(c,c') - K(c,c')]

        where J and K are computed from the ERI using:
          J(c,c') = <cc'|cc'> = eri[sp_c, sp_c', sp_c, sp_c']
          K(c,c') = <cc'|c'c> = eri[sp_c, sp_c', sp_c', sp_c]
                    (only for same-spin pairs)
        """
        eri = self._eri_4d

        # One-electron: sum of h1 diagonal for all core spin-orbitals
        e_h1 = 0.0
        for c in self.frozen_spinorb:
            e_h1 += self._h1_diag[c >> 1]

        # Two-electron: Coulomb J - Exchange K for all core pairs
        e_jk = 0.0
        core = self.frozen_spinorb
        for i in range(len(core)):
            ci = core[i]
            sp_i = ci >> 1
            sig_i = ci & 1
            for j in range(i + 1, len(core)):
                cj = core[j]
                sp_j = cj >> 1
                sig_j = cj & 1
                # Coulomb: always present
                e_jk += eri[sp_i, sp_j, sp_i, sp_j]
                # Exchange: only for same-spin pairs
                if sig_i == sig_j:
                    e_jk -= eri[sp_i, sp_j, sp_j, sp_i]

        return e_h1 + e_jk

    def _build_effective_h1(self) -> np.ndarray:
        """
        Build effective one-electron integrals for active space.

        For each pair of active spin-orbitals (a, b):
          h_eff(a, b) = h1(a, b) + Σ_{c∈core} [<ac|bc> - δ(σ_a,σ_c)·<ac|cb>]

        In spatial orbital notation (spin-diagonal):
          h_eff_spatial(α, β) = h1_spatial(α, β)
              + Σ_{γ∈frozen_spatial} [2·eri(α,γ,β,γ) - eri(α,γ,γ,β)]

        The factor 2 comes from summing over both alpha and beta spins
        of each frozen spatial orbital. Exchange only acts on same-spin,
        so each frozen spatial contributes 1× exchange (not 2×).

        Returns
        -------
        h_eff : np.ndarray, shape (n_spatial_total, n_spatial_total)
            Effective one-electron integrals in the full spatial index.
        """
        n_spatial = self.parent.n_sp // 2
        eri = self._eri_4d
        h_eff = self._H1_dense.copy()

        for gamma in self.frozen_spatial:
            for a in range(n_spatial):
                for b in range(n_spatial):
                    # 2J - K from each doubly-occupied frozen orbital
                    h_eff[a, b] += (
                        2.0 * eri[a, gamma, b, gamma]
                        - eri[a, gamma, gamma, b]
                    )

        return h_eff

    def _enumerate_active_sd_basis(self) -> None:
        """
        Enumerate active-space Slater determinants.

        Each active SD is a sorted tuple of n_active_electrons spin-orbital
        indices chosen from the active (non-frozen) spin-orbitals. The full
        SD is the frozen core + active electrons.

        We store active SDs using ORIGINAL spin-orbital indices (not
        re-indexed) so that ERI and H1 lookups remain valid.
        """
        t0 = time.perf_counter()

        # Generate all combinations of active spin-orbitals
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(self.active_spinorb, self.n_active_electrons)
        )
        self.n_sd: int = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }

        elapsed = time.perf_counter() - t0
        print(
            f"[FrozenCore] {self.n_sd:,} active SDs enumerated in {elapsed:.3f}s"
        )

    # ------------------------------------------------------------------
    # Hamiltonian assembly
    # ------------------------------------------------------------------

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Build the active-space FCI Hamiltonian via Slater-Condon rules.

        Uses the effective one-electron integrals (h_eff = h1 + core J/K)
        and the original ERIs restricted to active orbitals.

        The core energy E_core is NOT added here — it is a constant shift
        applied to eigenvalues after diagonalization.

        Returns
        -------
        H : csr_matrix, shape (n_sd, n_sd)
        """
        t0 = time.perf_counter()

        sd_basis = self.sd_basis
        sd_index = self._sd_index
        n_el = self.n_active_electrons
        n_sd = self.n_sd
        threshold = self.threshold
        h_eff = self._h1_eff
        eri = self._eri_4d

        # Precompute h_eff diagonal
        n_spatial = self.parent.n_sp // 2
        h_eff_diag = np.array([h_eff[i, i] for i in range(n_spatial)])

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        # Precompute occupied sets for fast membership tests
        occ_sets: List[frozenset] = [frozenset(sd) for sd in sd_basis]

        # --- Phase 1: Diagonal --- O(N_SD × n_el²)
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
                    # Coulomb
                    h_diag += eri[sp_i, sp_j, sp_i, sp_j]
                    # Exchange (same spin only)
                    if (pj & 1) == sig_i:
                        h_diag -= eri[sp_i, sp_j, sp_j, sp_i]

            if abs(h_diag) >= threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # --- Phase 2: Singles --- O(N_SD × n_el × n_active_sp)
        active_sp_arr = self.active_spinorb
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = occ_sets[I]

            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for r in active_sp_arr:
                    if r in occ_I:
                        continue
                    if (r & 1) != sig_p:
                        continue  # spin conservation

                    sp_r = r >> 1

                    # Slater-Condon single-excitation matrix element
                    # Using h_eff (already includes core J/K)
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

                    # Build target SD
                    new_sd = list(sd_I)
                    new_sd[kp] = r
                    new_sd_t = tuple(sorted(new_sd))
                    J = sd_index.get(new_sd_t)
                    if J is None or J <= I:
                        continue  # upper triangle only

                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        t_singles = time.perf_counter() - t0 - t_diag

        # --- Phase 3: Doubles --- O(N_SD × C(n_el,2) × n_active²)
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

                    # Iterate over virtual pairs from active space
                    for ir in range(len(active_sp_arr)):
                        r = active_sp_arr[ir]
                        if r in occ_I:
                            continue
                        sig_r = r & 1

                        for js in range(ir + 1, len(active_sp_arr)):
                            s = active_sp_arr[js]
                            if s in occ_I:
                                continue
                            sig_s = s & 1

                            # Spin conservation
                            if sig_p + sig_q != sig_r + sig_s:
                                continue

                            sp_r = r >> 1
                            sp_s = s >> 1

                            # Matrix element
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
                            J_idx = sd_index.get(new_sd_t)
                            if J_idx is None or J_idx <= I:
                                continue

                            phase = self._compute_double_phase(
                                sd_I, kp, kq, r, s
                            )
                            off_rows.append(I)
                            off_cols.append(J_idx)
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
            f"[FrozenCore] H assembled: shape={H.shape}, nnz={H.nnz:,}, "
            f"time={elapsed:.3f}s "
            f"(diag={t_diag:.1f}s, singles={t_singles:.1f}s, "
            f"doubles={t_doubles:.1f}s)"
        )
        return H.tocsr()

    # ------------------------------------------------------------------
    # Phase computation (Slater-Condon sign rules)
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_phase(sd: Tuple[int, ...], kp: int, r: int) -> float:
        """
        Fermionic sign for single excitation sd[kp] → r.

        Phase = (-1)^{kp + kr} where kr = #{a ∈ SD : a < r, a ≠ sd[kp]}.
        """
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _compute_double_phase(
        sd: Tuple[int, ...], kp: int, kq: int, r: int, s: int
    ) -> float:
        """
        Fermionic sign for double excitation sd[kp],sd[kq] → r,s.

        Computed by applying two sequential single excitations.
        """
        p = sd[kp]
        q = sd[kq]

        # First excitation: remove p, insert r
        n_swap_p = kp
        remaining_1 = [a for a in sd if a != p]
        kr = sum(1 for a in remaining_1 if a < r)
        phase1 = (-1) ** (n_swap_p + kr)

        # Second excitation: remove q, insert s
        intermediate = sorted(remaining_1[:kr] + [r] + remaining_1[kr:])
        kq_new = intermediate.index(q)
        remaining_2 = [a for a in intermediate if a != q]
        ks = sum(1 for a in remaining_2 if a < s)
        phase2 = (-1) ** (kq_new + ks)

        return float(phase1 * phase2)

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def solve(
        self,
        n_states: int = 1,
        add_nuclear_repulsion: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve active-space FCI with frozen core.

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
            Total energies: E_core + E_active [+ V_NN if molecular].
        eigvecs : np.ndarray, shape (n_sd, n_states)
            Active-space CI vectors.
        """
        H = self.assemble_hamiltonian()
        k = min(n_states, self.n_sd - 2)
        rng = np.random.RandomState(42)
        v0 = rng.randn(H.shape[0])
        eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
        order = np.argsort(eigvals)
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]

        # Add frozen-core energy
        eigvals = eigvals + self.E_core

        # Add nuclear repulsion if parent is a molecular index
        V_NN = 0.0
        if add_nuclear_repulsion and hasattr(self.parent, 'V_NN'):
            V_NN = self.parent.V_NN
            eigvals = eigvals + V_NN

        E_active = eigvals[0] - self.E_core - V_NN
        print(
            f"[FrozenCore] E_core = {self.E_core:.6f}, "
            f"E_active = {E_active:.6f}, "
            f"V_NN = {V_NN:.6f}, "
            f"E_total = {eigvals[0]:.6f} Ha"
        )
        return eigvals, eigvecs

    # ------------------------------------------------------------------
    # Convenience constructor
    # ------------------------------------------------------------------

    @classmethod
    def from_molecular(
        cls,
        Z_A: int,
        Z_B: int,
        nmax_A: int,
        nmax_B: int,
        R: float,
        n_core_A: int = 0,
        n_core_B: int = 0,
        vee_method: str = 'slater_full',
        fci_method: str = 'auto',
    ) -> 'FrozenCoreLatticeIndex':
        """
        Build a frozen-core molecular solver from scratch.

        Parameters
        ----------
        Z_A, Z_B : int
            Nuclear charges.
        nmax_A, nmax_B : int
            Basis truncation per atom.
        R : float
            Internuclear distance (bohr).
        n_core_A, n_core_B : int
            Number of core ELECTRONS on each atom (must be even).
            Example: Li has 2 core electrons (1s²), so n_core_A=2.
        vee_method : str
            V_ee method (default: 'slater_full').
        fci_method : str
            FCI method for parent construction.

        Returns
        -------
        FrozenCoreLatticeIndex
        """
        from .lattice_index import MolecularLatticeIndex

        if n_core_A % 2 != 0 or n_core_B % 2 != 0:
            raise ValueError("Core electron counts must be even (doubly-occupied orbitals)")

        n_electrons = Z_A + Z_B  # neutral molecule
        n_active = n_electrons - n_core_A - n_core_B

        # Check if full SD enumeration would OOM
        # Skip it when the full basis exceeds ~10M SDs
        from math import comb as _comb
        n_spatial_est = nmax_A**2 + nmax_B**2  # rough upper bound
        n_sp_est = 2 * n_spatial_est
        n_sd_est = _comb(n_sp_est, n_electrons)
        skip_sds = n_sd_est > 10_000_000

        # Build molecular index (skip SD enumeration if too large)
        parent = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=nmax_A, nmax_B=nmax_B,
            R=R,
            n_electrons=n_electrons,
            vee_method=vee_method,
            fci_method=fci_method,
            enumerate_sds=not skip_sds,
        )

        # Determine frozen spatial orbitals
        # Atom A orbitals: indices 0 to n_spatial_A - 1
        # Atom B orbitals: indices n_spatial_A to n_spatial - 1
        # Core orbitals are the lowest-n orbitals on each atom
        frozen_spatial = []

        # Atom A core: first n_core_A/2 spatial orbitals
        n_core_spatial_A = n_core_A // 2
        for i in range(n_core_spatial_A):
            frozen_spatial.append(i)

        # Atom B core: first n_core_B/2 spatial orbitals of atom B
        n_core_spatial_B = n_core_B // 2
        n_spatial_A = parent._n_spatial_A
        for i in range(n_core_spatial_B):
            frozen_spatial.append(n_spatial_A + i)

        return cls(parent, frozen_spatial, n_active)
