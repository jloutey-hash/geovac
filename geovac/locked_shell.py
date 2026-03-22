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


# ---------------------------------------------------------------------------
# Clementi-Raimondi effective nuclear charges
# Source: Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963)
# Key: (Z, n, l) -> Z_eff for that orbital
# ---------------------------------------------------------------------------
CLEMENTI_RAIMONDI_ZEFF = {
    (1, 1, 0): 1.000,   # H 1s
    (2, 1, 0): 1.688,   # He 1s
    (3, 1, 0): 2.691,   # Li 1s (core)
    (3, 2, 0): 1.279,   # Li 2s (valence)
    (3, 2, 1): 1.279,   # Li 2p
    (4, 1, 0): 3.685,   # Be 1s
    (4, 2, 0): 1.912,   # Be 2s
    (5, 1, 0): 4.680,   # B 1s
    (5, 2, 0): 2.576,   # B 2s
    (5, 2, 1): 2.421,   # B 2p
    (6, 1, 0): 5.673,   # C 1s
    (6, 2, 0): 3.217,   # C 2s
    (6, 2, 1): 3.136,   # C 2p
    (7, 1, 0): 6.665,   # N 1s
    (7, 2, 0): 3.847,   # N 2s
    (7, 2, 1): 3.834,   # N 2p
    (8, 1, 0): 7.658,   # O 1s
    (8, 2, 0): 4.492,   # O 2s
    (8, 2, 1): 4.453,   # O 2p
    (9, 1, 0): 8.650,   # F 1s
    (9, 2, 0): 5.128,   # F 2s
    (9, 2, 1): 5.100,   # F 2p
    (10, 1, 0): 9.642,  # Ne 1s
    (10, 2, 0): 5.758,  # Ne 2s
    (10, 2, 1): 5.758,  # Ne 2p
    (11, 3, 0): 2.507,  # Na 3s
    (17, 3, 1): 6.116,  # Cl 3p
}


def get_valence_zeff(Z: int) -> float:
    """Get Clementi-Raimondi Z_eff for the outermost valence orbital.

    Parameters
    ----------
    Z : int
        Atomic number.

    Returns
    -------
    float
        Z_eff for the valence orbital.  Falls back to 1.0 if not tabulated.
    """
    _valence = {
        1: 1.000,   # H  1s
        2: 1.688,   # He 1s
        3: 1.279,   # Li 2s
        4: 1.912,   # Be 2s
        5: 2.421,   # B  2p
        6: 3.136,   # C  2p
        7: 3.834,   # N  2p
        8: 4.453,   # O  2p
        9: 5.100,   # F  2p
        10: 5.758,  # Ne 2p
        11: 2.507,  # Na 3s
        17: 6.116,  # Cl 3p
    }
    return _valence.get(Z, 1.0)


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
        active_method: str = 'lcao',
        prolate_kwargs: Optional[Dict] = None,
        hyperspherical_kwargs: Optional[Dict] = None,
    ) -> None:
        from .lattice_index import MolecularLatticeIndex

        self.Z_A = Z_A
        self.Z_B = Z_B
        self.R = R
        self.n_electrons = n_electrons
        self.active_nmax = active_nmax
        self._active_method = active_method
        self._prolate_kwargs = prolate_kwargs or {}
        self._hyperspherical_kwargs = hyperspherical_kwargs or {}

        t0 = time.perf_counter()

        if active_method == 'hyperspherical':
            self._init_hyperspherical(locked_config, t0)
            return

        if active_method == 'prolate':
            self._init_prolate(locked_config, t0)
            return

        # --- LCAO path (original) ---
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

        # One-electron integrals (dense — n_spatial × n_spatial, always small)
        self._h1_diag = np.array(self._parent._h1_diag)
        self._H1_dense = np.asarray(self._parent._H1_spatial.todense())

        # Two-electron integrals: keep parent's sparse _eri dict as primary.
        # Dense _eri_4d is built lazily only when the matrix method needs it.
        # Selection rules (Gaunt c^k for same-atom, monopole F⁰ for cross-atom)
        # give O(M²) non-zero entries — densifying wastes O(M⁴) memory.
        self._eri_4d: Optional[np.ndarray] = None

        # Compute locked-shell energy and effective H1 from sparse dict
        self.E_locked = self._compute_locked_energy()
        self._h1_eff = self._build_effective_h1()

        # Enumerate active-space SDs
        self._enumerate_active_sds()

        elapsed = time.perf_counter() - t0

        full_sd_est = comb(self._parent.n_sp, n_electrons)
        n_eri = len(self._parent._eri)
        m4 = self.n_spatial ** 4
        eri_density = n_eri / m4 if m4 > 0 else 0.0
        print(
            f"[LockedShell] {len(self._locked_spatial)} locked spatial orbitals, "
            f"{len(self._active_spatial)} active spatial ({self.n_active_sp} spin-orbs), "
            f"{self.n_sd} active SDs "
            f"(vs ~{full_sd_est:,} full, ~{full_sd_est // max(self.n_sd, 1):,}x reduction), "
            f"E_locked = {self.E_locked:.6f} Ha, "
            f"ERI: {n_eri}/{m4} = {eri_density:.2%} density, "
            f"setup = {elapsed:.3f}s"
        )

    # ------------------------------------------------------------------
    # Prolate spheroidal active-space initialization
    # ------------------------------------------------------------------

    def _init_prolate(
        self,
        locked_config: Optional[Dict[int, List[Tuple[int, int]]]],
        t0: float,
    ) -> None:
        """Initialize the prolate spheroidal active-space path.

        Architecture:
          1. Compute Z_eff = Z - n_core for each atom (screened charges).
          2. Solve active MOs at Z_eff (correct orbital shapes for valence).
          3. Correct h_eff with:
             (a) Nuclear attraction difference: -(Z_bare - Z_eff)/r per nucleus
             (b) Core Coulomb: +2*J(core, active) per core orbital
             (c) Core exchange: -K(core, active) per core orbital
             Terms (a) and (b) nearly cancel at long range; the residual is
             the short-range penetration correction that breaks gerade/ungerade
             symmetry for heteronuclear molecules.
          4. E_locked = hydrogenic eigenvalues + intra-core V_ee +
             cross-nuclear attraction (analytical).
          5. V_NN uses the REAL (unscreened) nuclear charges.

        For LiH (locked Li 1s^2):
          Z_eff_A = 1, Z_eff_B = 1 (orbital solve — screened)
          Z_bare_A = 3, Z_bare_B = 1 (h_eff correction — real charges)
          Active MOs polarize toward H via the penetration correction.
        """
        from .prolate_active_space import ProlateActiveSpace

        # Use REAL charges for V_NN
        self.V_NN = self.Z_A * self.Z_B / self.R

        # Determine locked configuration
        if locked_config is None:
            locked_config = self._auto_detect_locked_simple()
        self._locked_config = locked_config

        # Build locked_orbitals list and compute screened charges
        locked_orbitals: List[Dict] = []
        n_locked = 0
        n_core_A = 0
        n_core_B = 0
        for atom_idx, shells in locked_config.items():
            Z = self.Z_A if atom_idx == 0 else self.Z_B
            atom_label = 'A' if atom_idx == 0 else 'B'
            for n_q, l_q in shells:
                n_elec = 2 * (2 * l_q + 1)
                n_locked += n_elec
                if atom_idx == 0:
                    n_core_A += n_elec
                else:
                    n_core_B += n_elec
                locked_orbitals.append({
                    'Z': float(Z), 'n': n_q, 'l': l_q,
                    'atom': atom_label, 'n_elec': n_elec,
                })

        self.n_locked_el = n_locked
        self.n_active_el = self.n_electrons - n_locked

        # Screened charges for the orbital solve
        Z_eff_A = max(float(self.Z_A - n_core_A), 0.01)
        Z_eff_B = max(float(self.Z_B - n_core_B), 0.01)

        # Active spec
        active_spec = self._prolate_kwargs.get('active_spec', None)
        if active_spec is None:
            active_spec = [(0, 0, 0), (0, 1, 0)]  # 1sigma, 2sigma

        # Filter out prolate-specific keys from kwargs
        pas_kwargs = {
            k: v for k, v in self._prolate_kwargs.items()
            if k not in ('core_spec', 'active_spec')
        }

        # Build ProlateActiveSpace:
        #   Z_A, Z_B = Z_eff (screened charges for orbital solve)
        #   Z_bare_A, Z_bare_B = real charges (for energy correction)
        #   locked_orbitals = core specifications (for core-active J/K + E_locked)
        self._prolate_active = ProlateActiveSpace(
            Z_A=Z_eff_A,
            Z_B=Z_eff_B,
            R=self.R,
            n_electrons=self.n_active_el,
            core_spec=[],
            active_spec=active_spec,
            locked_orbitals=locked_orbitals if locked_orbitals else None,
            Z_bare_A=float(self.Z_A),
            Z_bare_B=float(self.Z_B),
            **pas_kwargs,
        )
        self._prolate_active.build()

        # E_locked is computed by ProlateActiveSpace (analytical path)
        self.E_locked = self._prolate_active.E_locked

        elapsed = time.perf_counter() - t0
        print(
            f"[LockedShell/prolate] {n_locked}e locked, "
            f"{self.n_active_el}e active "
            f"(Z_eff={Z_eff_A:.2f},{Z_eff_B:.2f} / "
            f"Z_bare={self.Z_A},{self.Z_B}), "
            f"E_locked = {self.E_locked:.6f} Ha, "
            f"setup = {elapsed:.3f}s"
        )

    # ------------------------------------------------------------------
    # Hyperspherical (Level 4) active-space initialization
    # ------------------------------------------------------------------

    def _init_hyperspherical(
        self,
        locked_config: Optional[Dict[int, List[Tuple[int, int]]]],
        t0: float,
    ) -> None:
        """Initialize the Level 4 molecule-frame hyperspherical active-space path.

        Architecture:
          1. Determine locked shells and count core electrons per atom.
          2. Compute screened charges Z_eff = Z - n_core for each atom.
          3. Compute E_locked analytically (hydrogenic eigenvalues + intra-core
             V_ee + cross-nuclear attraction).
          4. Store Z_eff for the Level 4 solver, which treats the 2 active
             electrons in the screened nuclear potential.
          5. V_NN uses bare (unscreened) nuclear charges.

        The Level 4 solver (level4_multichannel.py) handles V_ee between
        the active electrons as an angular eigenvalue in hyperspherical
        coordinates — no integral approximation needed.

        For LiH (locked Li 1s^2):
          Z_eff_A = 1, Z_eff_B = 1 (screened)
          Active: 2 electrons in hyperspherical coordinates
          V_ee: angular eigenvalue (exact in the partial-wave basis)
        """
        # V_NN uses bare charges
        self.V_NN = self.Z_A * self.Z_B / self.R

        # Determine locked configuration
        if locked_config is None:
            locked_config = self._auto_detect_locked_simple()
        self._locked_config = locked_config

        # Count core electrons per atom and compute screened charges
        n_locked = 0
        n_core_A = 0
        n_core_B = 0
        locked_orbitals: List[Dict] = []

        for atom_idx, shells in locked_config.items():
            Z = self.Z_A if atom_idx == 0 else self.Z_B
            atom_label = 'A' if atom_idx == 0 else 'B'
            for n_q, l_q in shells:
                n_elec = 2 * (2 * l_q + 1)
                n_locked += n_elec
                if atom_idx == 0:
                    n_core_A += n_elec
                else:
                    n_core_B += n_elec
                locked_orbitals.append({
                    'Z': float(Z), 'n': n_q, 'l': l_q,
                    'atom': atom_label, 'n_elec': n_elec,
                })

        self.n_locked_el = n_locked
        self.n_active_el = self.n_electrons - n_locked

        if self.n_active_el != 2:
            raise ValueError(
                f"Hyperspherical active space requires exactly 2 active "
                f"electrons, got {self.n_active_el} "
                f"({self.n_electrons} total - {n_locked} locked)"
            )

        # Screened charges for the Level 4 solver.
        # zeff_mode selects how Z_eff is computed:
        #   'screened'  — Z_eff = Z - n_core (integer screening, default)
        #   'clementi'  — Clementi-Raimondi atomic HF Z_eff (asymmetric)
        #   'bare'      — Z_eff = Z (no screening, for testing)
        zeff_mode = self._hyperspherical_kwargs.get('zeff_mode', 'screened')
        self._zeff_mode = zeff_mode

        if zeff_mode == 'clementi':
            # Apply Clementi-Raimondi Z_eff ONLY to atoms with locked cores.
            # Atoms with no locked electrons keep bare Z (no screening).
            self._Z_eff_A = (get_valence_zeff(self.Z_A) if n_core_A > 0
                             else float(self.Z_A))
            self._Z_eff_B = (get_valence_zeff(self.Z_B) if n_core_B > 0
                             else float(self.Z_B))
        elif zeff_mode == 'bare':
            self._Z_eff_A = float(self.Z_A)
            self._Z_eff_B = float(self.Z_B)
        else:  # 'screened'
            self._Z_eff_A = max(float(self.Z_A - n_core_A), 0.01)
            self._Z_eff_B = max(float(self.Z_B - n_core_B), 0.01)

        # Store core counts for penetration correction
        self._n_core_A = n_core_A
        self._n_core_B = n_core_B

        # Compute E_locked analytically
        self.E_locked = self._compute_locked_energy_analytic_hyper(
            locked_orbitals
        )

        elapsed = time.perf_counter() - t0
        print(
            f"[LockedShell/hyperspherical] {n_locked}e locked, "
            f"{self.n_active_el}e active "
            f"(Z_eff={self._Z_eff_A:.2f},{self._Z_eff_B:.2f} "
            f"[{zeff_mode}] / Z_bare={self.Z_A},{self.Z_B}), "
            f"E_core = {self._E_core_isolated:.6f}, "
            f"V_cross_nuc = {self._V_cross_nuc:.6f}, "
            f"setup = {elapsed:.3f}s"
        )

    def _compute_core_energy_isolated(
        self,
        locked_orbitals: List[Dict],
    ) -> float:
        """Compute R-independent isolated core energy.

        E_core = sum_c n_c * (-Z_c^2/2n_c^2)     [hydrogenic eigenvalues]
               + sum_{c<c'} V_ee(c,c')              [intra-core repulsion]

        No cross-nuclear terms — those go in V_cross_nuc (R-dependent).

        For Li 1s^2:
          E_core = 2*(-9/2) + (5*3/8) = -9.0 + 1.875 = -7.125 Ha
        """
        # One-electron energies
        e_h1 = 0.0
        for spec in locked_orbitals:
            Z, n_q = spec['Z'], spec['n']
            n_elec = spec['n_elec']
            e_h1 += n_elec * (-Z**2 / (2.0 * n_q**2))

        # Intra-core V_ee via exact Slater F^0
        e_jk = 0.0
        for i in range(len(locked_orbitals)):
            si = locked_orbitals[i]
            ni_elec = si['n_elec']
            # Intra-orbital J (alpha-beta pair)
            if ni_elec >= 2:
                F0_ii = self._compute_slater_f0_pair(
                    si['Z'], si['n'], si['l'],
                    si['Z'], si['n'], si['l'],
                )
                e_jk += F0_ii

            # Inter-orbital J-K
            for j in range(i + 1, len(locked_orbitals)):
                sj = locked_orbitals[j]
                nj_elec = sj['n_elec']
                if si['Z'] != sj['Z']:
                    continue  # cross-atom core-core V_ee small at typical R
                F0_ij = self._compute_slater_f0_pair(
                    si['Z'], si['n'], si['l'],
                    sj['Z'], sj['n'], sj['l'],
                )
                e_jk += ni_elec * nj_elec * F0_ij
                if si['l'] == 0 and sj['l'] == 0:
                    e_jk -= (ni_elec // 2) * (nj_elec // 2) * 2 * F0_ij

        return e_h1 + e_jk

    def _compute_core_cross_nuclear(
        self,
        locked_orbitals: List[Dict],
        R: float,
    ) -> float:
        """Compute R-dependent core-to-other-nucleus attraction.

        V_cross = -sum_c n_c * Z_other * <c|1/r_other|c>

        For Li 1s^2 at distance R from H:
          <1s_Z|1/r_B|1s_Z> = (1/R)[1 - (1+ZR)exp(-2ZR)]

        This is an attractive (negative) R-dependent term.  At large R,
        V_cross -> -n_core * Z_other / R (screens V_NN).
        """
        e_cross = 0.0
        for spec in locked_orbitals:
            Z_core = spec['Z']
            n_q, l_q = spec['n'], spec['l']
            n_elec = spec['n_elec']
            atom = spec['atom']
            Z_other = self.Z_B if atom == 'A' else self.Z_A

            if Z_other == 0 or l_q != 0:
                continue

            # <n s_Z | 1/r_other | n s_Z>
            if n_q == 1:
                zr = Z_core * R
                v_cross = (1.0 / R) * (
                    1.0 - (1.0 + zr) * np.exp(-2.0 * zr)
                )
            else:
                # Higher-n s-orbitals: numerical quadrature
                from scipy.integrate import quad

                def _integrand(r: float, n: int, Z: float) -> float:
                    if n == 2:
                        R_nl = ((1.0 / (2.0 * np.sqrt(2.0))) * Z**1.5
                                * (2.0 - Z * r) * np.exp(-Z * r / 2.0))
                    elif n == 3:
                        rho = 2.0 * Z * r / 3.0
                        R_nl = ((2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5
                                * (27.0 - 18.0 * rho + 2.0 * rho**2)
                                * np.exp(-rho / 2.0))
                    else:
                        return 0.0
                    return R_nl**2 * r**2 / np.sqrt(r**2 + R**2)

                v_cross, _ = quad(
                    _integrand, 0, 50.0 / Z_core,
                    args=(n_q, Z_core), limit=200,
                )

            e_cross -= n_elec * Z_other * v_cross

        return e_cross

    def _compute_locked_energy_analytic_hyper(
        self,
        locked_orbitals: List[Dict],
    ) -> float:
        """Compute E_locked = E_core(isolated) + V_cross_nuc(R).

        Split into:
          E_core:      R-independent (hydrogenic + intra-core V_ee)
          V_cross_nuc: R-dependent (core attracted to other nucleus)

        Both components stored separately for energy decomposition.
        """
        self._E_core_isolated = self._compute_core_energy_isolated(
            locked_orbitals
        )
        self._V_cross_nuc = self._compute_core_cross_nuclear(
            locked_orbitals, self.R
        )
        return self._E_core_isolated + self._V_cross_nuc

    def _auto_detect_locked_simple(self) -> Dict[int, List[Tuple[int, int]]]:
        """Auto-detect locked shells without building MolecularLatticeIndex.

        Same logic as _auto_detect_locked but doesn't need _parent.
        """
        config: Dict[int, List[Tuple[int, int]]] = {}
        for atom_idx, Z in enumerate([self.Z_A, self.Z_B]):
            if Z <= 1:
                continue
            electrons_remaining = Z
            shells_to_lock = []
            for n in range(1, Z + 1):
                for l in range(n):
                    capacity = 2 * (2 * l + 1)
                    if electrons_remaining >= capacity:
                        electrons_remaining -= capacity
                        shells_to_lock.append((n, l))
                    else:
                        break
                if electrons_remaining == 0:
                    break
            if electrons_remaining > 0 and len(shells_to_lock) >= 1:
                config[atom_idx] = list(shells_to_lock)
            elif len(shells_to_lock) > 1:
                config[atom_idx] = shells_to_lock[:-1]
        return config

    def _compute_locked_energy_analytic(
        self,
        locked_config: Dict[int, List[Tuple[int, int]]],
    ) -> float:
        """Compute locked core energy from exact hydrogenic eigenvalues.

        E_core = sum_c [-Z^2/(2n^2)] * n_elec_c
                + sum_{c<c'} [J(c,c') - K(c,c')]

        For Li 1s^2: E = 2 * (-9/2) + F0(1s,1s) = -9 + 5*3/8 = -7.125 Ha
        """
        from scipy.integrate import quad

        e_h1 = 0.0
        all_core_orbs: List[Tuple[int, int, int, int]] = []  # (Z, n, l, n_elec)

        for atom_idx, shells in locked_config.items():
            Z = self.Z_A if atom_idx == 0 else self.Z_B
            for n_q, l_q in shells:
                n_elec = 2 * (2 * l_q + 1)
                e_h1 += n_elec * (-Z**2 / (2.0 * n_q**2))
                all_core_orbs.append((Z, n_q, l_q, n_elec))

        # For multiple core orbitals, add intra-core V_ee
        # For Li 1s^2 (single orbital, 2 electrons): J(1s,1s) via Slater F0
        e_jk = 0.0
        for i in range(len(all_core_orbs)):
            Zi, ni, li, ni_elec = all_core_orbs[i]
            # Intra-orbital J (alpha-beta pair within same spatial orbital)
            if ni_elec >= 2:
                F0_ii = self._compute_slater_f0_pair(Zi, ni, li, Zi, ni, li)
                # C(n_elec, 2) pairs, each contributing J; minus K for same-spin
                # For 2 electrons in one orbital: 1 pair, J only (opposite spin)
                e_jk += F0_ii

            # Inter-orbital J-K
            for j in range(i + 1, len(all_core_orbs)):
                Zj, nj, lj, nj_elec = all_core_orbs[j]
                if Zi != Zj:
                    continue  # cross-atom core-core V_ee is small at typical R
                F0_ij = self._compute_slater_f0_pair(Zi, ni, li, Zi, nj, lj)
                # ni_elec * nj_elec J pairs, minus same-spin K
                e_jk += ni_elec * nj_elec * F0_ij
                # Exchange: only s-s contributes in monopole approx
                if li == 0 and lj == 0:
                    e_jk -= (ni_elec // 2) * (nj_elec // 2) * 2 * F0_ij

        return e_h1 + e_jk

    @staticmethod
    def _compute_slater_f0_pair(
        Z1: float, n1: int, l1: int,
        Z2: float, n2: int, l2: int,
    ) -> float:
        """Compute Slater F^0 integral between two hydrogenic orbitals.

        F0 = integral_0^inf integral_0^inf |R_{n1l1}(r1)|^2 |R_{n2l2}(r2)|^2
             * (1/r_>) * r1^2 * r2^2 dr1 dr2

        Only supports l=0 currently. Uses exact analytical results.
        """
        from scipy.integrate import quad

        if l1 != 0 or l2 != 0:
            return 0.0  # Monopole approximation: l>0 → 0

        # Hydrogenic radial wavefunction squared * r^2
        def rho_r2(n: int, Z: float, r: float) -> float:
            if n == 1:
                R_nl = 2.0 * Z**1.5 * np.exp(-Z * r)
            elif n == 2:
                R_nl = (1.0 / (2.0 * np.sqrt(2.0))) * Z**1.5 * (2.0 - Z * r) * np.exp(-Z * r / 2.0)
            elif n == 3:
                rho_val = 2.0 * Z * r / 3.0
                R_nl = (2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5 * (27.0 - 18.0 * rho_val + 2.0 * rho_val**2) * np.exp(-rho_val / 2.0)
            else:
                return 0.0
            return R_nl**2 * r**2

        def inner(r1: float) -> float:
            rho1 = rho_r2(n1, Z1, r1)
            if rho1 < 1e-30:
                return 0.0

            def integrand_lt(r2: float) -> float:
                return rho_r2(n2, Z2, r2) / r1

            def integrand_gt(r2: float) -> float:
                return rho_r2(n2, Z2, r2) / r2

            val_lt, _ = quad(integrand_lt, 0, r1, limit=100)
            val_gt, _ = quad(integrand_gt, r1, 50.0 / Z2, limit=100)
            return rho1 * (val_lt + val_gt)

        result, _ = quad(inner, 0, 50.0 / Z1, limit=200)
        return result

    # ------------------------------------------------------------------
    # Orbital classification
    # ------------------------------------------------------------------

    def _auto_detect_locked(self) -> Dict[int, List[Tuple[int, int]]]:
        """
        Auto-detect closed shells based on electron configuration.

        For each atom, identify (n, l) subshells that would be fully
        occupied in the ground-state Aufbau configuration and lock them
        if they are deep core orbitals (n < valence n).

        Only considers orbitals within the atom's basis (n <= nmax).
        """
        config: Dict[int, List[Tuple[int, int]]] = {}
        nmax_list = [self._parent.nmax_A, self._parent.nmax_B]

        for atom_idx, Z in enumerate([self.Z_A, self.Z_B]):
            if Z <= 1:
                continue  # H has no core

            nmax = nmax_list[atom_idx]

            # Aufbau filling within basis: n up to nmax only
            electrons_remaining = Z
            shells_to_lock = []
            for n in range(1, nmax + 1):
                for l in range(n):
                    capacity = 2 * (2 * l + 1)
                    if electrons_remaining >= capacity:
                        electrons_remaining -= capacity
                        shells_to_lock.append((n, l))
                    else:
                        break
                if electrons_remaining == 0:
                    break

            # Lock core subshells, keep valence unlocked.
            # If electrons_remaining > 0, there's a partially occupied shell
            # above all filled shells — all filled shells are core, lock them.
            # If electrons_remaining == 0, the last filled shell is valence —
            # lock all but that one.
            if electrons_remaining > 0 and len(shells_to_lock) >= 1:
                # All filled subshells are core (valence is partially filled above)
                config[atom_idx] = list(shells_to_lock)
            elif len(shells_to_lock) > 1:
                # Last filled subshell is valence; lock all below it
                config[atom_idx] = shells_to_lock[:-1]

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

    @staticmethod
    def eri_selection_rules() -> str:
        """
        Document the geometric selection rules that enforce O(M^2) ERI
        density in the molecular basis (Paper 0 principle).

        Same-atom ERIs (Gaunt c^k, Paper 14 Sec III.B):
          - Triangle inequality: |l_a - l_c| <= k <= l_a + l_c
          - Parity: (l_a + l_c + k) must be even
          - m-conservation: m_a + m_b = m_c + m_d
          - k-matching: multipole rank k must be same for bra and ket
          Result: ~4% ERI density at nmax=3

        Cross-atom ERIs (monopole F^0, Paper 12 Sec III):
          - Only (a,b,a,b) Coulomb and (a,b,b,a) exchange patterns
          - J(n_a,l_a; n_b,l_b) depends on radial parts only (no m)
          - Exchange K restricted to s-s pairs (Mulliken approximation)
          - Higher multipole cross-atom integrals neglected
          Result: O(M_A * M_B) entries, not O(M_A^2 * M_B^2)

        Combined: O(M^2) total ERI entries, giving ~0.5% density at
        nmax=3 (28 spatial orbitals). This sparsity propagates through
        the locked-shell adapter to DirectCISolver without densification.
        """
        return "O(M^2) ERI density from Gaunt + monopole selection rules"

    def _get_eri(self, a: int, b: int, c: int, d: int) -> float:
        """Look up ERI from parent's sparse dict. O(1) per call."""
        return self._parent._eri.get((a, b, c, d), 0.0)

    def _ensure_eri_dense(self) -> np.ndarray:
        """
        Lazily build dense 4D ERI array from parent's sparse dict.

        Only needed by the matrix-method assemble_hamiltonian() which
        requires fast random access in the Slater-Condon inner loop.
        The Direct CI path bypasses this entirely.
        """
        if self._eri_4d is None:
            n = self.n_spatial
            eri = np.zeros((n, n, n, n))
            for (a, b, c, d), val in self._parent._eri.items():
                eri[a, b, c, d] = val
            self._eri_4d = eri
        return self._eri_4d

    def _compute_locked_energy(self) -> float:
        """
        Compute energy of locked shells: h1 + J - K.

        Uses parent's sparse _eri dict — O(n_locked²) lookups.
        E_locked = sum_c h1[c,c] + sum_{c<c'} [J(c,c') - K(c,c')]
        """
        get = self._get_eri
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
                e_jk += get(sp_i, sp_j, sp_i, sp_j)
                if (cj & 1) == sig_i:
                    e_jk -= get(sp_i, sp_j, sp_j, sp_i)

        return e_h1 + e_jk

    def _build_effective_h1(self) -> np.ndarray:
        """
        Build effective one-electron integrals: h1 + locked-shell J/K.

        Uses parent's sparse _eri dict — O(n_locked × n_spatial²) lookups.
        h_eff[a,b] = h1[a,b] + sum_{gamma in locked_spatial}
                     [2*eri(a,gamma,b,gamma) - eri(a,gamma,gamma,b)]
        """
        get = self._get_eri
        h_eff = self._H1_dense.copy()

        for gamma in self._locked_spatial:
            for a in range(self.n_spatial):
                for b in range(self.n_spatial):
                    h_eff[a, b] += (
                        2.0 * get(a, gamma, b, gamma)
                        - get(a, gamma, gamma, b)
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
        eri = self._ensure_eri_dense()
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
    # Direct CI adapter
    # ------------------------------------------------------------------

    def _build_direct_ci_adapter(self) -> '_ActiveSpaceAdapter':
        """
        Build a lightweight adapter that remaps active orbitals to
        contiguous indices and duck-types as a LatticeIndex for
        DirectCISolver.
        """
        return _ActiveSpaceAdapter(self)

    # ------------------------------------------------------------------
    # Hyperspherical (Level 4) solver
    # ------------------------------------------------------------------

    def _solve_hyperspherical(
        self,
        n_states: int = 1,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Solve the 2-electron active space using Level 4 hyperspherical coords.

        Uses Z_eff charges in the multipole nuclear expansion (long-range
        screening) plus a core penetration correction V_pen(r) for
        short-range core-active Coulomb (the attractive extra potential
        felt by an electron penetrating the locked core).

        V_pen(r_A) = -(n_core/r_A)(1 + Z_core r_A) exp(-2 Z_core r_A)

        Combined with -Z_eff/r_A from multipole, this gives the full
        effective potential: -Z_bare/r_A + V_core(r_A).

        Returns total energies = E_locked + E_elec + V_NN(Z_bare).

        Supported hyperspherical_kwargs:
          l_max : int (default 4) — max angular momentum per electron
          m_max : int (default 0) — max |m| (0=sigma-only)
          n_alpha : int (default 200) — FD grid for correlation angle
          n_Re : int (default 400) — grid for hyperradial equation
          R_e_min, R_e_max : float — hyperradial boundaries
          n_coupled : int (default 1) — radial solver mode
              1: single-channel adiabatic
             -1: direct 2D (variational, slower)
          origin : str (default 'charge_center' for heteronuclear)
              'midpoint', 'charge_center', or explicit float z0
          zeff_mode : str (default 'screened')
              'screened': Z_eff = Z - n_core (integer, symmetric for LiH)
              'clementi': Clementi-Raimondi atomic HF Z_eff (asymmetric)
              'bare': Z_eff = Z (no screening)
          verbose : bool (default True)
        """
        from .level4_multichannel import solve_level4_h2_multichannel

        kw = dict(self._hyperspherical_kwargs)

        # Defaults
        kw.setdefault('l_max', 4)
        kw.setdefault('m_max', 0)
        kw.setdefault('n_alpha', 200)
        kw.setdefault('n_Re', 400)
        kw.setdefault('n_coupled', 1)
        kw.setdefault('verbose', True)

        # Auto-scale R_e_max so the dissociation wavefunction fits on grid.
        # At dissociation, electrons sit near nuclei at ±R/2 from midpoint,
        # giving R_e ~ R/sqrt(2).  Need R_e_max > R/sqrt(2) with margin.
        R_e_max_min = max(15.0, self.R * 0.8 + 5.0)
        kw.setdefault('R_e_max', R_e_max_min)
        # Scale n_Re if R_e_max grew, to maintain grid density
        if kw['R_e_max'] > 15.0 and 'n_Re' not in self._hyperspherical_kwargs:
            kw['n_Re'] = max(400, int(kw['R_e_max'] / 15.0 * 400))

        # Effective charges for the active-electron solver (set by zeff_mode
        # in _init_hyperspherical: 'screened', 'clementi', or 'bare').
        Z_A_solver = float(self._Z_eff_A)
        Z_B_solver = float(self._Z_eff_B)
        homonuclear = abs(Z_A_solver - Z_B_solver) < 1e-10

        # Use charge-center origin for heteronuclear by default
        if not homonuclear:
            kw.setdefault('origin', 'charge_center')
        else:
            kw.setdefault('origin', 'midpoint')

        # Don't pass internal config keys to the solver
        kw.pop('E_exact', None)
        kw.pop('D_e_exact', None)
        kw.pop('pk_E_val', None)
        kw.pop('pk_C_core', None)
        kw.pop('pk_beta_core', None)
        kw.pop('zeff_mode', None)

        # Core penetration correction (off by default).
        # V_pen(r) = -(n_core/r)(1 + Z_core r) exp(-2 Z_core r)
        # adds short-range attraction where active electrons penetrate the
        # locked core.  In practice this causes the hyperspherical solver
        # to collapse into the core well (no PES minimum).  Enable with
        # hyperspherical_kwargs={'penetration': True} for diagnostics.
        core_pots: Optional[List[dict]] = None
        if kw.pop('penetration', False):
            core_pots = []
            for atom_idx, shells in self._locked_config.items():
                Z = self.Z_A if atom_idx == 0 else self.Z_B
                atom_label = 'A' if atom_idx == 0 else 'B'
                n_core_on_atom = sum(
                    2 * (2 * l_q + 1) for (n_q, l_q) in shells
                )
                for n_q, l_q in shells:
                    if l_q != 0:
                        continue
                    core_pots.append({
                        'Z_core': float(Z),
                        'n_core': n_core_on_atom,
                        'n_q': n_q,
                        'l_q': l_q,
                        'atom': atom_label,
                    })

        result = solve_level4_h2_multichannel(
            R=self.R,
            Z_A=Z_A_solver,
            Z_B=Z_B_solver,
            E_exact=None,
            D_e_exact=None,
            core_potentials=core_pots,
            pk_potentials=None,
            **kw,
        )

        E_elec = result['E_elec']
        # Decomposed total energy:
        #   E_core_isolated : R-independent (hydrogenic + intra-core V_ee)
        #   V_cross_nuc     : R-dependent core-to-other-nucleus attraction
        #   E_elec          : active 2e energy from Level 4 solver
        #   V_NN            : bare nuclear repulsion
        E_total = (self._E_core_isolated + self._V_cross_nuc
                   + E_elec + self.V_NN)

        # Store decomposition for diagnostics
        self._E_elec = E_elec

        print(
            f"[LockedShell/hyperspherical] "
            f"E_core = {self._E_core_isolated:.6f}, "
            f"V_cross_nuc = {self._V_cross_nuc:.6f}, "
            f"E_elec(active) = {E_elec:.6f}, "
            f"V_NN = {self.V_NN:.6f}, "
            f"E_total = {E_total:.6f} Ha"
        )

        eigvals = np.array([E_total])
        # Store the Level 4 result for diagnostics
        self._level4_result = result
        # Hyperspherical wavefunction is in (R_e, alpha) space, not SD space
        eigvecs = result.get('wavefunction', np.array([1.0]))

        return eigvals, eigvecs

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def solve(
        self,
        n_states: int = 1,
        fci_method: str = 'auto',
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve active-space problem with locked shells.

        Parameters
        ----------
        n_states : int
            Number of lowest eigenvalues to compute.
        fci_method : str
            Assembly method: 'auto' (direct for N_SD >= 5000),
            'direct' (excitation-driven), or 'matrix' (dense loops).
            Ignored when active_method='prolate'.

        Returns
        -------
        eigvals : np.ndarray
            Total energies: E_locked + E_active + V_NN.
        eigvecs : np.ndarray
            Active-space CI vectors.
        """
        # --- Hyperspherical (Level 4) path ---
        if self._active_method == 'hyperspherical':
            return self._solve_hyperspherical(n_states=n_states)

        # --- Prolate path ---
        if self._active_method == 'prolate':
            return self._prolate_active.solve(
                E_locked=self.E_locked,
                V_NN=self.V_NN,
                n_states=n_states,
            )

        # --- LCAO path ---
        method = fci_method
        if method == 'auto':
            method = 'direct' if self.n_sd >= 5000 else 'matrix'

        if method == 'direct':
            from .direct_ci import DirectCISolver
            adapter = self._build_direct_ci_adapter()
            solver = DirectCISolver(adapter)
            eigvals, eigvecs = solver.solve(n_states=n_states)
        else:
            H = self.assemble_hamiltonian()

            if self.n_sd <= 2:
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


class _ActiveSpaceAdapter:
    """
    Adapter that remaps LockedShellMolecule's active orbitals to
    contiguous indices, duck-typing as a LatticeIndex for DirectCISolver.

    Sparsity preservation (Paper 0 principle):
      Parent's _eri dict has O(M²) entries due to geometric selection rules
      (Gaunt c^k for same-atom, monopole F⁰ for cross-atom). This adapter
      remaps the sparse dict directly — O(|_eri|) — instead of iterating
      O(M_active⁴) over a dense array and filtering by threshold.

    DirectCISolver accesses:
      .sd_basis, .n_sd, .n_sp, .n_electrons, ._sd_index,
      ._H1_spatial (sparse), ._eri (dict), ._h1_diag (list),
      .threshold, ._compute_phase(), ._compute_double_phase()
    """

    def __init__(self, locked_mol: 'LockedShellMolecule') -> None:
        active_spatial = locked_mol._active_spatial
        n_active_spatial = len(active_spatial)

        # Spatial orbital remap: original index -> contiguous 0..N-1
        sp_remap: Dict[int, int] = {
            old: new for new, old in enumerate(active_spatial)
        }

        # Spin-orbital remap: old_sp*2+spin -> new_sp*2+spin
        spinorb_remap: Dict[int, int] = {}
        for old_sp, new_sp in sp_remap.items():
            spinorb_remap[old_sp * 2] = new_sp * 2
            spinorb_remap[old_sp * 2 + 1] = new_sp * 2 + 1

        self.n_sp = n_active_spatial * 2
        self.n_electrons = locked_mol.n_active_el
        self.n_sd = locked_mol.n_sd
        self.threshold = locked_mol._parent.threshold

        # Remap SD basis to contiguous indices
        self.sd_basis: List[Tuple[int, ...]] = [
            tuple(sorted(spinorb_remap[s] for s in sd))
            for sd in locked_mol.sd_basis
        ]
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }

        # Remap effective H1 to active-only block
        h1_full = locked_mol._h1_eff
        h1_active = np.zeros((n_active_spatial, n_active_spatial))
        for i, old_i in enumerate(active_spatial):
            for j, old_j in enumerate(active_spatial):
                h1_active[i, j] = h1_full[old_i, old_j]
        self._H1_spatial = csr_matrix(h1_active)

        # h1_diag uses effective (dressed) diagonal
        self._h1_diag: List[float] = [
            h1_full[old, old] for old in active_spatial
        ]

        # Remap ERI from parent's sparse dict — O(|_eri|) not O(M⁴).
        # Filter to active-active entries and remap indices.
        active_sp_set = frozenset(active_spatial)
        self._eri: Dict[Tuple[int, int, int, int], float] = {}
        for (a, b, c, d), val in locked_mol._parent._eri.items():
            if (a in active_sp_set and b in active_sp_set
                    and c in active_sp_set and d in active_sp_set):
                self._eri[
                    (sp_remap[a], sp_remap[b], sp_remap[c], sp_remap[d])
                ] = val

        n_parent = len(locked_mol._parent._eri)
        n_active = len(self._eri)
        m4 = n_active_spatial ** 4
        print(
            f"[Adapter] ERI remap: {n_active}/{n_parent} parent entries "
            f"are active-active (density {n_active}/{m4} = "
            f"{n_active / m4:.2%} vs dense M^4)"
        )

    @staticmethod
    def _compute_phase(
        sd: Tuple[int, ...], kp: int, r: int
    ) -> float:
        """Fermionic sign for single excitation."""
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _compute_double_phase(
        sd: Tuple[int, ...], kp: int, kq: int, r: int, s: int
    ) -> float:
        """Fermionic sign for double excitation."""
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
