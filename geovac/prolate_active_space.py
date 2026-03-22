"""
ProlateActiveSpace — Molecular active space in prolate spheroidal MOs
=====================================================================

Replaces the LCAO graph active space with eigenstates of the prolate
spheroidal one-electron problem.  Bonding/antibonding splitting emerges
from the eigenvalue spectrum, not from off-diagonal matrix elements.

Architecture:
  1. Solve 1-electron prolate spheroidal eigenproblems for ALL states
     (core + active) at the physical Z_A, Z_B, R.
  2. Lock core orbitals (computed on the prolate grid), compute E_locked
     from their eigenvalues + intra-core V_ee.
  3. Build effective h1 for active orbitals: bare eigenvalues + core-active
     Coulomb/exchange integrals evaluated on the prolate grid.
  4. Build a small CI Hamiltonian in the active MO basis and diagonalize.

For LiH:
  Core MO:    1sigma (n_ang=0, n_rad=0) = Li 1s, 2 electrons
  Active MOs: 2sigma (n_ang=1, n_rad=0) = bonding (Li 2s + H 1s)
              3sigma (n_ang=0, n_rad=1) = antibonding (Li 2s - H 1s)

Reference: Paper 11 (prolate spheroidal lattice).

Author: GeoVac Development Team
Date: March 2026
"""

import time
from itertools import combinations
from math import comb
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import eigh

from geovac.prolate_scf import (
    get_orbital_on_grid,
    compute_vee_integral,
)


# ============================================================
# Atomic orbital on prolate grid (for core-active coupling)
# ============================================================

def atomic_orbital_on_prolate_grid(
    Z: float,
    n: int,
    l: int,
    atom_label: str,
    R: float,
    xi: np.ndarray,
    eta: np.ndarray,
    w_xi: np.ndarray,
    w_eta: np.ndarray,
) -> Dict:
    """Evaluate a hydrogenic orbital on a prolate spheroidal grid.

    Places a hydrogenic (n, l, m=0) orbital centered on atom A or B
    onto the (xi, eta) grid. Returns a dict compatible with
    compute_vee_integral.

    Parameters
    ----------
    Z : float
        Nuclear charge for the orbital (e.g. 3 for Li 1s).
    n : int
        Principal quantum number.
    l : int
        Angular momentum (only l=0 supported currently).
    atom_label : str
        'A' or 'B' — which nucleus the orbital is centered on.
    R : float
        Internuclear distance (bohr).
    xi, eta, w_xi, w_eta : np.ndarray
        Prolate spheroidal grid coordinates and weights.
    """
    if l != 0:
        raise NotImplementedError("Only s-orbitals (l=0) supported for core")

    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    if atom_label == 'A':
        r = (R / 2.0) * (XI - ETA)
    else:
        r = (R / 2.0) * (XI + ETA)

    # Hydrogenic radial wavefunction R_{n,l=0}(r)
    if n == 1:
        psi_r = 2.0 * Z**1.5 * np.exp(-Z * r)
    elif n == 2:
        psi_r = (1.0 / (2.0 * np.sqrt(2.0))) * Z**1.5 * (2.0 - Z * r) * np.exp(-Z * r / 2.0)
    elif n == 3:
        rho = 2.0 * Z * r / 3.0
        psi_r = (2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5 * (27.0 - 18.0 * rho + 2.0 * rho**2) * np.exp(-rho / 2.0)
    else:
        from scipy.special import genlaguerre
        from math import factorial
        rho = 2.0 * Z * r / n
        L = genlaguerre(n - 1, 1)(rho)
        norm = np.sqrt((2.0 * Z / n)**3 * factorial(n - 1) / (2.0 * n * factorial(n)**2))
        psi_r = norm * np.exp(-rho / 2.0) * rho * L

    psi_2d = psi_r / np.sqrt(4.0 * np.pi)

    # Normalize on the prolate grid
    J = (R / 2.0)**3 * (XI**2 - ETA**2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
    norm_sq = np.sum(psi_2d**2 * J * W_XI * W_ETA) * 2.0 * np.pi
    if norm_sq > 1e-30:
        psi_2d /= np.sqrt(norm_sq)

    return {
        'psi': psi_2d,
        'xi': xi, 'eta': eta,
        'w_xi': w_xi, 'w_eta': w_eta,
        'R': R,
    }


# ============================================================
# Symmetry labels
# ============================================================

def _make_label(m: int, n_angular: int, n_radial: int) -> str:
    """Generate orbital symmetry label from quantum numbers."""
    # Count the ordinal: n_radial+1 for the n-th state of that symmetry
    ordinal = n_radial + 1
    if m == 0:
        sym = 'sigma_g' if n_angular % 2 == 0 else 'sigma_u'
    elif m == 1:
        sym = 'pi_u' if n_angular % 2 == 0 else 'pi_g'
    else:
        sym = f'delta(m={m})'
    return f'{ordinal}{sym}'


# ============================================================
# Slater F^0 integral (analytical, for locked core energy)
# ============================================================

def _slater_f0(
    Z1: float, n1: int, l1: int,
    Z2: float, n2: int, l2: int,
) -> float:
    """Compute Slater F^0 integral between two hydrogenic s-orbitals.

    Uses known analytical results for efficiency.
    F^0(1s,1s;Z) = 5Z/8, F^0(1s,2s;Z) = 17Z/81, F^0(2s,2s;Z) = 77Z/512.
    """
    if l1 != 0 or l2 != 0:
        return 0.0

    # Use exact analytical values when possible (same Z, same atom)
    if abs(Z1 - Z2) < 1e-10:
        Z = Z1
        key = (min(n1, n2), max(n1, n2))
        if key == (1, 1):
            return 5.0 * Z / 8.0
        elif key == (1, 2):
            return 17.0 * Z / 81.0
        elif key == (2, 2):
            return 77.0 * Z / 512.0
        elif key == (1, 3):
            return 815.0 * Z / 19683.0  # exact
        elif key == (2, 3):
            return 3.0 * Z * 2363.0 / 131072.0  # exact
        elif key == (3, 3):
            return 7861.0 * Z / 131072.0  # exact

    # Fall back to numerical integration
    from scipy.integrate import quad

    def _R_nl_sq_r2(n: int, Z: float, r: float) -> float:
        if n == 1:
            R_nl = 2.0 * Z**1.5 * np.exp(-Z * r)
        elif n == 2:
            R_nl = (1.0 / (2.0 * np.sqrt(2.0))) * Z**1.5 * (2.0 - Z * r) * np.exp(-Z * r / 2.0)
        elif n == 3:
            rho = 2.0 * Z * r / 3.0
            R_nl = (2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5 * (27.0 - 18.0 * rho + 2.0 * rho**2) * np.exp(-rho / 2.0)
        else:
            return 0.0
        return R_nl**2 * r**2

    def inner(r1: float) -> float:
        rho1 = _R_nl_sq_r2(n1, Z1, r1)
        if rho1 < 1e-30:
            return 0.0
        val_lt, _ = quad(lambda r2: _R_nl_sq_r2(n2, Z2, r2) / r1, 0, r1, limit=100)
        val_gt, _ = quad(lambda r2: _R_nl_sq_r2(n2, Z2, r2) / r2, r1, 50.0 / Z2, limit=100)
        return rho1 * (val_lt + val_gt)

    result, _ = quad(inner, 0, 50.0 / Z1, limit=200)
    return result


# ============================================================
# Cross-nuclear attraction integral
# ============================================================

def _cross_nuclear_attraction(
    Z: float, n: int, l: int, R: float,
) -> float:
    """Compute <n,l,m=0 on A | 1/r_B | n,l,m=0 on A> for hydrogenic orbital.

    For 1s(Z) at distance R from point charge:
      <1s|1/r_B|1s> = (1/R) * [1 - (1 + Z*R)*exp(-2*Z*R)]

    For 2s(Z):
      <2s|1/r_B|2s> = (1/R) * [1 - (1 + Z*R + Z²R²/3 + Z³R³/12 + ...)*exp(-Z*R)]

    Parameters
    ----------
    Z : nuclear charge of the orbital
    n : principal quantum number
    l : angular momentum
    R : distance to the other nucleus

    Returns
    -------
    float : <nl|1/r_B|nl> (always positive)
    """
    if R < 1e-10:
        return 0.0

    if n == 1 and l == 0:
        # Exact for 1s: Coulomb integral of 1s density with point charge at R
        x = 2.0 * Z * R
        return (1.0 / R) * (1.0 - (1.0 + x) * np.exp(-x))

    elif n == 2 and l == 0:
        # Exact for 2s(Z): <2s|1/r|2s> via integration
        # |R_{2s}|^2 = Z^3/8 * (2-Zr)^2 * exp(-Zr) ; rho(r) = |R_{2s}|^2 * r^2/(4pi)
        x = Z * R
        return (1.0 / R) * (1.0 - np.exp(-x) * (
            1.0 + x + x**2 / 3.0 + x**3 / 12.0 + x**4 / 48.0
        ))

    else:
        # Numerical integration for general (n, l)
        from scipy.integrate import quad

        def rho_r2(r: float) -> float:
            if n == 3 and l == 0:
                rho_val = 2.0 * Z * r / 3.0
                R_nl = (2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5 * (
                    27.0 - 18.0 * rho_val + 2.0 * rho_val**2
                ) * np.exp(-rho_val / 2.0)
            else:
                return 0.0
            return R_nl**2 * r**2

        # V(R) = integral_0^R rho(r)*r^2/R dr + integral_R^inf rho(r)*r dr
        # where rho includes the 4pi angular part already in R_nl^2 * r^2
        # Actually: V(R) = 4pi * [1/R * int_0^R rho(r) r^2 dr + int_R^inf rho(r) r dr]
        # With rho(r) = |R_nl(r)|^2 / (4pi), the 4pi cancels:
        # V(R) = 1/R * int_0^R |R_nl|^2 r^2 dr + int_R^inf |R_nl|^2 r dr
        def inner_part(r: float) -> float:
            return rho_r2(r)

        def outer_part(r: float) -> float:
            if n == 3 and l == 0:
                rho_val = 2.0 * Z * r / 3.0
                R_nl = (2.0 / (81.0 * np.sqrt(3.0))) * Z**1.5 * (
                    27.0 - 18.0 * rho_val + 2.0 * rho_val**2
                ) * np.exp(-rho_val / 2.0)
            else:
                return 0.0
            return R_nl**2 * r

        v_inner, _ = quad(inner_part, 0, R, limit=200)
        v_outer, _ = quad(outer_part, R, 100.0 / Z, limit=200)
        return v_inner / R + v_outer


# ============================================================
# ProlateActiveSpace
# ============================================================

class ProlateActiveSpace:
    """
    Active space solver using prolate spheroidal molecular orbitals.

    Solves 1-electron prolate spheroidal eigenproblems for core and active
    states, then builds and diagonalizes a multi-electron CI matrix.

    Parameters
    ----------
    Z_A, Z_B : float
        Nuclear charges (full, unscreened).
    R : float
        Internuclear distance (bohr).
    n_electrons : int
        Number of active electrons.
    core_spec : list of (m, n_angular, n_radial), optional
        Core MOs to lock. Each is doubly occupied (2 electrons).
        Default: [] (no core).
    active_spec : list of (m, n_angular, n_radial), optional
        Active MOs for CI. Default: [(0,0,0), (0,1,0)] for sigma_g, sigma_u.
    N_xi_solve : int
        Radial grid for eigenvalue solve.
    N_grid : int
        Quadrature grid size for V_ee integrals.
    xi_max_grid : float
        Maximum xi for quadrature.
    verbose : bool
        Print progress info.
    """

    def __init__(
        self,
        Z_A: float,
        Z_B: float,
        R: float,
        n_electrons: int = 2,
        core_spec: Optional[List[Tuple[int, int, int]]] = None,
        active_spec: Optional[List[Tuple[int, int, int]]] = None,
        N_xi_solve: int = 5000,
        N_grid: int = 80,
        xi_max_grid: float = 15.0,
        verbose: bool = True,
        # Legacy API compatibility
        orbital_spec: Optional[List[Tuple[int, int]]] = None,
        # Hybrid mode: atomic core orbitals instead of prolate core
        locked_orbitals: Optional[List[Dict]] = None,
        # Bare nuclear charges for heteronuclear correction
        # When set, Z_A/Z_B are treated as Z_eff for the eigenvalue solve,
        # and the difference (Z_bare - Z_eff) is added as a correction to h_eff.
        Z_bare_A: Optional[float] = None,
        Z_bare_B: Optional[float] = None,
    ) -> None:
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.R = R
        self.n_electrons = n_electrons
        self.N_xi_solve = N_xi_solve
        self.N_grid = N_grid
        self.xi_max_grid = xi_max_grid
        self.verbose = verbose

        # Bare nuclear charges for heteronuclear correction.
        # Z_A, Z_B = screened charges for orbital solve.
        # Z_bare_A, Z_bare_B = actual nuclear charges for energy correction.
        self.Z_bare_A = Z_bare_A if Z_bare_A is not None else Z_A
        self.Z_bare_B = Z_bare_B if Z_bare_B is not None else Z_B

        # Handle legacy orbital_spec (2-tuples)
        if orbital_spec is not None and active_spec is None:
            active_spec = [(m, n_ang, 0) for m, n_ang in orbital_spec]
        if core_spec is None:
            core_spec = []
        if active_spec is None:
            active_spec = [(0, 0, 0), (0, 1, 0)]

        self.core_spec = core_spec
        self.active_spec = active_spec
        self.n_core_orbitals = len(core_spec)
        self.n_orbitals = len(active_spec)
        self.n_spinorb = 2 * self.n_orbitals

        # Hybrid mode: use atomic orbitals for core instead of prolate MOs
        # Each entry: {'Z': float, 'n': int, 'l': int, 'atom': str, 'n_elec': int}
        self._locked_orbitals = locked_orbitals or []
        self._use_atomic_core = len(self._locked_orbitals) > 0

        # Internal state
        self._core_orbitals: List[Dict] = []
        self._active_orbitals: List[Dict] = []
        self._orbital_energies: np.ndarray = np.array([])
        self._orbital_labels: List[str] = []
        self._core_labels: List[str] = []
        self._eri: Optional[np.ndarray] = None
        self._h_eff: Optional[np.ndarray] = None
        self.E_locked: float = 0.0

    def build(self) -> 'ProlateActiveSpace':
        """Solve orbitals, compute integrals, build effective Hamiltonian."""
        t0 = time.perf_counter()

        self._solve_all_orbitals()
        self._compute_locked_energy()
        self._compute_eri()
        self._compute_effective_h1()

        dt = time.perf_counter() - t0
        if self.verbose:
            labels = ', '.join(self._orbital_labels)
            energies_str = ', '.join(f'{e:.6f}' for e in self._orbital_energies)
            core_str = ', '.join(self._core_labels) if self._core_labels else 'none'
            print(
                f"[ProlateActive] core=[{core_str}], "
                f"active=[{labels}], "
                f"eps = [{energies_str}], "
                f"{self.n_electrons}e active, "
                f"E_locked = {self.E_locked:.6f}, "
                f"setup = {dt:.3f}s"
            )
        return self

    # ------------------------------------------------------------------
    # Orbital computation
    # ------------------------------------------------------------------

    def _solve_all_orbitals(self) -> None:
        """Solve for core and active 1-electron MOs.

        In hybrid mode (_use_atomic_core=True): core orbitals are hydrogenic
        wavefunctions evaluated on the prolate grid. Active orbitals are
        prolate spheroidal eigenstates.

        In standard mode: both core and active are prolate eigenstates.
        """
        # First solve active orbitals (we need their grid for atomic core)
        self._active_orbitals = []
        energies = []
        labels = []
        for m, n_ang, n_rad in self.active_spec:
            orb = get_orbital_on_grid(
                R=self.R, Z_A=self.Z_A, Z_B=self.Z_B,
                n_angular=n_ang, n_radial=n_rad,
                N_xi_solve=self.N_xi_solve,
                N_xi_grid=self.N_grid, N_eta_grid=self.N_grid,
                xi_max_grid=self.xi_max_grid, m=m,
            )
            self._active_orbitals.append(orb)
            energies.append(orb['E_elec'])
            labels.append(_make_label(m, n_ang, n_rad))

        self._orbital_energies = np.array(energies)
        self._orbital_labels = labels

        # Core orbitals
        self._core_orbitals = []
        self._core_labels = []

        if self._use_atomic_core:
            # Hybrid mode: use atomic orbitals on prolate grid
            ref = self._active_orbitals[0]  # use first active orbital's grid
            for spec in self._locked_orbitals:
                orb = atomic_orbital_on_prolate_grid(
                    Z=spec['Z'], n=spec['n'], l=spec['l'],
                    atom_label=spec['atom'], R=self.R,
                    xi=ref['xi'], eta=ref['eta'],
                    w_xi=ref['w_xi'], w_eta=ref['w_eta'],
                )
                self._core_orbitals.append(orb)
                self._core_labels.append(
                    f"{spec['atom']}:{spec['n']}{'spdf'[spec['l']]}"
                )
        else:
            # Standard mode: prolate eigenstates as core
            for m, n_ang, n_rad in self.core_spec:
                orb = get_orbital_on_grid(
                    R=self.R, Z_A=self.Z_A, Z_B=self.Z_B,
                    n_angular=n_ang, n_radial=n_rad,
                    N_xi_solve=self.N_xi_solve,
                    N_xi_grid=self.N_grid, N_eta_grid=self.N_grid,
                    xi_max_grid=self.xi_max_grid, m=m,
                )
                self._core_orbitals.append(orb)
                self._core_labels.append(_make_label(m, n_ang, n_rad))

    # ------------------------------------------------------------------
    # Locked core energy
    # ------------------------------------------------------------------

    def _compute_locked_energy(self) -> None:
        """Compute E_locked from core orbitals.

        In hybrid mode (atomic core): uses exact hydrogenic eigenvalues
        + analytical Slater F^0 integrals. This is R-independent and gives
        the correct dissociation limit.

        In standard mode (prolate core): uses prolate eigenvalues + numerical
        J/K on the prolate grid.

        For Li 1s^2: E_locked = 2*(-9/2) + F0(1s,1s) = -9 + 15/8 = -7.125 Ha
        """
        if not self._core_orbitals and not self._locked_orbitals:
            self.E_locked = 0.0
            return

        if self._use_atomic_core:
            self._compute_locked_energy_analytic()
        else:
            self._compute_locked_energy_prolate()

    def _compute_locked_energy_analytic(self) -> None:
        """Compute E_locked from exact hydrogenic eigenvalues + Slater F^0
        + cross-nuclear attraction.

        E_locked = sum_c n_c * (-Z_c^2/2n_c^2)     [hydrogenic eigenvalues]
                 + sum_{c<c'} V_ee(c,c')              [intra-core repulsion]
                 + sum_c n_c * V_cross_nuc(c)          [core-OTHER-nucleus attraction]

        The cross-nuclear term is critical: it provides the R-dependent
        attraction that creates the PES minimum. A Li 1s electron on atom A
        is attracted to the H nucleus on atom B with:
          V_cross = -Z_B * <1s_A|1/r_B|1s_A>
                  = -Z_B/R * [1 - (1 + Z_A*R)*exp(-2*Z_A*R)]

        This screens the nuclear repulsion: V_NN + V_cross ≈ (Z_A-n_core)*Z_B/R.
        """
        e_h1 = 0.0
        for spec in self._locked_orbitals:
            Z, n, l = spec['Z'], spec['n'], spec['l']
            n_elec = spec.get('n_elec', 2)
            e_h1 += n_elec * (-Z**2 / (2.0 * n**2))

        # Intra-core V_ee via exact Slater F^0
        e_jk = 0.0
        n_core = len(self._locked_orbitals)
        for i in range(n_core):
            si = self._locked_orbitals[i]
            ni_elec = si.get('n_elec', 2)
            if ni_elec >= 2:
                F0_ii = _slater_f0(si['Z'], si['n'], si['l'],
                                   si['Z'], si['n'], si['l'])
                e_jk += F0_ii

            for j in range(i + 1, n_core):
                sj = self._locked_orbitals[j]
                nj_elec = sj.get('n_elec', 2)
                F0_ij = _slater_f0(si['Z'], si['n'], si['l'],
                                   sj['Z'], sj['n'], sj['l'])
                e_jk += ni_elec * nj_elec * F0_ij
                if si['l'] == 0 and sj['l'] == 0:
                    e_jk -= (ni_elec // 2) * (nj_elec // 2) * 2 * F0_ij

        # Cross-nuclear attraction: core electrons on atom X attracted to
        # the OTHER nucleus Y. For 1s(Z) on atom X at distance R from Y:
        #   <1s|Z_Y/r_Y|1s> = Z_Y/R * [1 - (1 + Z*R)*exp(-2*Z*R)]
        e_cross = 0.0
        for spec in self._locked_orbitals:
            Z_core = spec['Z']
            n_q = spec['n']
            l_q = spec['l']
            n_elec = spec.get('n_elec', 2)
            atom = spec['atom']

            # Which nucleus is the OTHER one?
            if atom == 'A':
                Z_other = self.Z_B
            else:
                Z_other = self.Z_A

            if Z_other == 0 or self.R < 1e-10:
                continue

            # Nuclear attraction integral for hydrogenic ns orbital
            v_cross = _cross_nuclear_attraction(Z_core, n_q, l_q, self.R)
            e_cross -= n_elec * Z_other * v_cross

        self.E_locked = e_h1 + e_jk + e_cross
        if self.verbose:
            print(
                f"[ProlateActive] Analytic E_locked = {self.E_locked:.6f} Ha "
                f"(h1={e_h1:.4f}, Vee={e_jk:.4f}, Vcross={e_cross:.4f})"
            )

    def _compute_locked_energy_prolate(self) -> None:
        """Compute E_locked from prolate core MO eigenvalues + numerical J/K."""
        # One-electron: 2 electrons per core MO
        e_h1 = sum(2.0 * orb['E_elec'] for orb in self._core_orbitals)

        # Intra-core V_ee
        e_jk = 0.0
        n_core = len(self._core_orbitals)
        for i in range(n_core):
            J_ii = compute_vee_integral(
                self._core_orbitals[i], self._core_orbitals[i],
                self._core_orbitals[i], self._core_orbitals[i],
            )
            e_jk += J_ii

            for j in range(i + 1, n_core):
                J_ij = compute_vee_integral(
                    self._core_orbitals[i], self._core_orbitals[j],
                    self._core_orbitals[i], self._core_orbitals[j],
                )
                K_ij = compute_vee_integral(
                    self._core_orbitals[i], self._core_orbitals[j],
                    self._core_orbitals[j], self._core_orbitals[i],
                )
                e_jk += 4 * J_ij - 2 * K_ij

        self.E_locked = e_h1 + e_jk

    # ------------------------------------------------------------------
    # Two-electron integrals (active-active)
    # ------------------------------------------------------------------

    def _compute_eri(self) -> None:
        """Compute the (ij|kl) ERI tensor in the active MO basis."""
        n = self.n_orbitals
        eri = np.zeros((n, n, n, n))

        for i in range(n):
            for j in range(i, n):
                for k in range(n):
                    for l_idx in range(k, n):
                        if (i, j) > (k, l_idx):
                            continue
                        val = compute_vee_integral(
                            self._active_orbitals[i],
                            self._active_orbitals[j],
                            self._active_orbitals[k],
                            self._active_orbitals[l_idx],
                        )
                        # 8-fold symmetry for real orbitals
                        for a, b, c, d in [
                            (i, j, k, l_idx), (j, i, l_idx, k),
                            (k, l_idx, i, j), (l_idx, k, j, i),
                            (i, j, l_idx, k), (j, i, k, l_idx),
                            (k, l_idx, j, i), (l_idx, k, i, j),
                        ]:
                            eri[a, b, c, d] = val

        self._eri = eri

    # ------------------------------------------------------------------
    # Effective one-electron Hamiltonian
    # ------------------------------------------------------------------

    def _compute_effective_h1(self) -> None:
        """Build h_eff for the active space.

        Three modes, selected automatically:

        1. **No core** (homonuclear H₂, HeH⁺, etc.):
           h_eff = diag(eps), no corrections.

        2. **Atomic core, Z_bare == Z_A** (bare-charge active MOs):
           h_eff = diag(eps) + analytical core J − numerical core K.

        3. **Atomic core, Z_bare != Z_A** (Z_eff active MOs + penetration):
           h_eff = diag(eps) + <p|V_pen|q>
           where V_pen combines nuclear correction and core Coulomb into
           a single short-range local potential:

             V_pen(r_A) = −(n_core/r_A)(1 + Z_c·r_A) exp(−2Z_c·r_A)

           This eliminates the large numerical cancellation between
           −dZ/r_A (nuclear correction) and +V_J(r_A) (core Coulomb).
           Exchange is added numerically from the core orbital on the grid.
        """
        n = self.n_orbitals
        h_eff = np.diag(self._orbital_energies.copy())

        has_core = self._use_atomic_core and len(self._locked_orbitals) > 0
        dZ_A = self.Z_bare_A - self.Z_A
        dZ_B = self.Z_bare_B - self.Z_B
        has_nuclear_corr = abs(dZ_A) > 1e-10 or abs(dZ_B) > 1e-10

        if has_core and has_nuclear_corr:
            # Mode 3: penetration potential (Z_eff + V_pen)
            self._apply_penetration_potential(h_eff, n)
        elif has_core:
            # Mode 2: analytical core J + numerical K (bare charges)
            self._compute_core_active_analytic(h_eff, n)
        elif has_nuclear_corr:
            # Nuclear correction without core (unusual but handle it)
            self._apply_nuclear_correction(h_eff, n)
        elif self._core_orbitals and not self._use_atomic_core:
            # Prolate core: numerical J/K
            for core_orb in self._core_orbitals:
                for p in range(n):
                    J_cp = compute_vee_integral(
                        core_orb, self._active_orbitals[p],
                        core_orb, self._active_orbitals[p],
                    )
                    h_eff[p, p] += 2.0 * J_cp
                    for q in range(p, n):
                        K_cpq = compute_vee_integral(
                            core_orb, self._active_orbitals[p],
                            self._active_orbitals[q], core_orb,
                        )
                        h_eff[p, q] -= K_cpq
                        if p != q:
                            h_eff[q, p] -= K_cpq

        if self.verbose and (self._core_orbitals or self._locked_orbitals):
            h_eff_diag = [h_eff[i, i] for i in range(n)]
            print(
                f"[ProlateActive] h_eff diagonal (with corrections): "
                f"[{', '.join(f'{v:.4f}' for v in h_eff_diag)}]"
            )

        # --- Orbital rotation to natural heteronuclear MOs ---
        # If h_eff has significant off-diagonal elements (from the nuclear
        # correction breaking g/u symmetry), the Z_eff MOs are not the natural
        # basis.  Diagonalize h_eff to get natural MOs and rotate the ERI.
        # This absorbs the symmetry-breaking into orbital shapes, so CI handles
        # only correlation, not orbital polarization.
        off_diag_norm = np.linalg.norm(h_eff - np.diag(np.diag(h_eff)))
        if off_diag_norm > 1e-8:
            eps_nat, U = eigh(h_eff)
            h_eff = np.diag(eps_nat)

            # Rotate ERI: (pq|rs)' = sum_{abcd} U[a,p]*U[b,q]*U[c,r]*U[d,s]*(ab|cd)
            if self._eri is not None:
                eri_old = self._eri
                eri_new = np.einsum('ap,bq,cr,ds,abcd->pqrs', U, U, U, U, eri_old)
                self._eri = eri_new

            if self.verbose:
                print(
                    f"[ProlateActive] Orbital rotation: "
                    f"off-diag norm = {off_diag_norm:.4f}, "
                    f"natural eps = [{', '.join(f'{e:.4f}' for e in eps_nat)}]"
                )

        self._h_eff = h_eff

    # ------------------------------------------------------------------
    # Analytical core-active Coulomb
    # ------------------------------------------------------------------

    def _apply_nuclear_correction(
        self, h_eff: np.ndarray, n: int,
    ) -> None:
        """Add nuclear attraction correction when Z_bare != Z_eff (no core)."""
        dZ_A = self.Z_bare_A - self.Z_A
        dZ_B = self.Z_bare_B - self.Z_B
        R = self.R
        ref = self._active_orbitals[0]
        xi, eta = ref['xi'], ref['eta']
        w_xi, w_eta = ref['w_xi'], ref['w_eta']

        XI, ETA = np.meshgrid(xi, eta, indexing='ij')
        W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')

        weight_A = 2.0 * np.pi * (R**2 / 4.0) * (XI + ETA) * W_XI * W_ETA
        weight_B = 2.0 * np.pi * (R**2 / 4.0) * (XI - ETA) * W_XI * W_ETA

        for p in range(n):
            psi_p = self._active_orbitals[p]['psi']
            for q in range(p, n):
                psi_q = self._active_orbitals[q]['psi']
                psi_pq = psi_p * psi_q
                vnuc_corr = 0.0
                if abs(dZ_A) > 1e-10:
                    vnuc_corr -= dZ_A * np.sum(psi_pq * weight_A)
                if abs(dZ_B) > 1e-10:
                    vnuc_corr -= dZ_B * np.sum(psi_pq * weight_B)
                h_eff[p, q] += vnuc_corr
                if p != q:
                    h_eff[q, p] += vnuc_corr

    def _apply_penetration_potential(
        self, h_eff: np.ndarray, n: int,
    ) -> None:
        """Add combined penetration correction for Z_eff MOs with atomic core.

        Instead of computing the nuclear correction (−dZ/r) and core Coulomb
        (V_J) separately — both ~2/r at long range, nearly canceling — compute
        the NET local potential:

          V_pen(r_A) = V_J(r_A) − n_core/r_A
                     = −(n_core/r_A)(1 + Z_c·r_A) exp(−2Z_c·r_A)

        This is attractive (negative), short-range (decays as exp(−2Z_c·r)),
        and well-conditioned numerically.

        Exchange is added from the numerical core orbital on the grid.
        """
        R = self.R
        ref = self._active_orbitals[0]
        xi, eta = ref['xi'], ref['eta']
        w_xi, w_eta = ref['w_xi'], ref['w_eta']

        XI, ETA = np.meshgrid(xi, eta, indexing='ij')
        W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
        J_vol = (R / 2.0)**3 * (XI**2 - ETA**2)

        for spec in self._locked_orbitals:
            Z_c = spec['Z']
            n_q = spec['n']
            l_q = spec['l']
            n_elec = spec.get('n_elec', 2)
            atom = spec['atom']

            if l_q != 0:
                # Non-s core: fall back to separate nuclear + core J/K
                self._apply_nuclear_correction(h_eff, n)
                core_idx = self._locked_orbitals.index(spec)
                if core_idx < len(self._core_orbitals):
                    core_orb = self._core_orbitals[core_idx]
                    for p in range(n):
                        J_cp = compute_vee_integral(
                            core_orb, self._active_orbitals[p],
                            core_orb, self._active_orbitals[p])
                        h_eff[p, p] += 2.0 * J_cp
                        for q in range(p, n):
                            K_cpq = compute_vee_integral(
                                core_orb, self._active_orbitals[p],
                                self._active_orbitals[q], core_orb)
                            h_eff[p, q] -= K_cpq
                            if p != q:
                                h_eff[q, p] -= K_cpq
                continue

            # Distance to the atom hosting this core shell
            if atom == 'A':
                r_core = (R / 2.0) * (XI - ETA)
            else:
                r_core = (R / 2.0) * (XI + ETA)
            r_safe = np.maximum(r_core, 1e-12)

            # Combined penetration potential:
            # V_pen = -(n_elec/r) * (1 + zeta*r) * exp(-2*zeta*r)
            # where zeta = Z_c for 1s (the orbital exponent)
            if n_q == 1:
                zeta = Z_c
                V_pen = -(n_elec / r_safe) * (
                    1.0 + zeta * r_safe
                ) * np.exp(-2.0 * zeta * r_safe)
            elif n_q == 2:
                # For 2s: V_pen = V_J(2s) - n_elec/r
                # V_J(2s) = (n_elec/r)[1 - exp(-Zr)(1+Zr+Z²r²/3+Z³r³/12+Z⁴r⁴/48)]
                x = Z_c * r_safe
                V_pen = -(n_elec / r_safe) * np.exp(-x) * (
                    1.0 + x + x**2 / 3.0 + x**3 / 12.0 + x**4 / 48.0
                )
            else:
                # Higher shells: just use -n_elec/r * exp(-2Z_c*r/n) as approx
                zeta = Z_c / n_q
                V_pen = -(n_elec / r_safe) * np.exp(-2.0 * zeta * r_safe)

            # Integrate V_pen against active MO densities
            weight_pen = 2.0 * np.pi * V_pen * J_vol * W_XI * W_ETA

            for p in range(n):
                psi_p = self._active_orbitals[p]['psi']
                for q in range(p, n):
                    psi_q = self._active_orbitals[q]['psi']
                    corr = np.sum(psi_p * psi_q * weight_pen)
                    h_eff[p, q] += corr
                    if p != q:
                        h_eff[q, p] += corr

            if self.verbose:
                pen_diag = []
                for p in range(n):
                    psi_p = self._active_orbitals[p]['psi']
                    pen_diag.append(np.sum(psi_p**2 * weight_pen))
                print(
                    f"[ProlateActive] Penetration V_pen ({atom}:"
                    f"{n_q}{'spdf'[l_q]}, {n_elec}e): "
                    f"diag = [{', '.join(f'{v:.4f}' for v in pen_diag)}]"
                )

            # Exchange: numerical (on prolate grid)
            core_idx = self._locked_orbitals.index(spec)
            if core_idx < len(self._core_orbitals):
                core_orb = self._core_orbitals[core_idx]
                for p in range(n):
                    for q in range(p, n):
                        K_cpq = compute_vee_integral(
                            core_orb, self._active_orbitals[p],
                            self._active_orbitals[q], core_orb,
                        )
                        h_eff[p, q] -= K_cpq
                        if p != q:
                            h_eff[q, p] -= K_cpq

    def _compute_core_active_analytic(
        self, h_eff: np.ndarray, n: int,
    ) -> None:
        """Add core-active Coulomb (analytical) and exchange (numerical).

        For a closed ns² core on atom X with nuclear charge Z_core:

          V_J(r) = (n_elec / r_X) * [1 - (1 + Z_c*r_X)*exp(-2*Z_c*r_X)]

        where Z_c = Z_core/n (the orbital exponent, Z/n for the ns orbital).
        This Coulomb potential is smooth and exact; we integrate it against
        active MO densities on the prolate grid.

        Exchange uses the numerical compute_vee_integral with the core orbital
        placed on the grid.  Exchange is ~50% of Coulomb and less sensitive to
        grid resolution.
        """
        ref = self._active_orbitals[0]
        xi = ref['xi']
        eta = ref['eta']
        w_xi = ref['w_xi']
        w_eta = ref['w_eta']
        R = self.R

        XI, ETA = np.meshgrid(xi, eta, indexing='ij')
        W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
        J_vol = (R / 2.0) ** 3 * (XI ** 2 - ETA ** 2)

        for spec in self._locked_orbitals:
            Z_c = spec['Z']
            n_q = spec['n']
            l_q = spec['l']
            n_elec = spec.get('n_elec', 2)
            atom = spec['atom']

            if l_q != 0:
                # For non-s core, fall back to numerical (rare case)
                core_idx = self._locked_orbitals.index(spec)
                if core_idx < len(self._core_orbitals):
                    core_orb = self._core_orbitals[core_idx]
                    for p in range(n):
                        J_cp = compute_vee_integral(
                            core_orb, self._active_orbitals[p],
                            core_orb, self._active_orbitals[p],
                        )
                        h_eff[p, p] += 2.0 * J_cp
                        for q in range(p, n):
                            K_cpq = compute_vee_integral(
                                core_orb, self._active_orbitals[p],
                                self._active_orbitals[q], core_orb,
                            )
                            h_eff[p, q] -= K_cpq
                            if p != q:
                                h_eff[q, p] -= K_cpq
                continue

            # --- Analytical Coulomb potential of ns² shell ---
            # Orbital exponent: zeta = Z/n for hydrogenic ns
            zeta = Z_c / n_q
            if atom == 'A':
                r_core = (R / 2.0) * (XI - ETA)  # distance to A
            else:
                r_core = (R / 2.0) * (XI + ETA)  # distance to B

            # V_J(r) = n_elec * (1/r) * [1 - (1 + 2*zeta*r)*exp(-2*2*zeta*r)]
            # Wait — for 1s(Z): V_J(r) = n_elec * (1/r)*[1-(1+Z*r)*exp(-2Z*r)]
            # For ns(Z): the density is more complex. Use the exact 1/r_A form
            # from _cross_nuclear_attraction but as a function of r on the grid.
            # For 1s: V(r) = (1/r)[1 - (1 + Z*r)*exp(-2Z*r)]
            # For 2s: V(r) = (1/r)[1 - exp(-Z*r)(1+Z*r+Z²r²/3+Z³r³/12+Z⁴r⁴/48)]
            r_safe = np.maximum(r_core, 1e-12)

            if n_q == 1:
                x = 2.0 * Z_c * r_safe
                V_J = n_elec * (1.0 / r_safe) * (
                    1.0 - (1.0 + x) * np.exp(-x)
                )
            elif n_q == 2:
                x = Z_c * r_safe  # Note: Z_c, not Z_c/2
                V_J = n_elec * (1.0 / r_safe) * (
                    1.0 - np.exp(-x) * (
                        1.0 + x + x**2 / 3.0 + x**3 / 12.0 + x**4 / 48.0
                    )
                )
            else:
                # Fallback: use 1/R point-charge approximation
                V_J = n_elec / r_safe

            # Integrate V_J against active MO densities
            weight_J = 2.0 * np.pi * V_J * J_vol * W_XI * W_ETA
            for p in range(n):
                psi_p = self._active_orbitals[p]['psi']
                for q in range(p, n):
                    psi_q = self._active_orbitals[q]['psi']
                    J_pq = np.sum(psi_p * psi_q * weight_J)
                    h_eff[p, q] += J_pq
                    if p != q:
                        h_eff[q, p] += J_pq

            if self.verbose:
                J_diag = []
                for p in range(n):
                    psi_p = self._active_orbitals[p]['psi']
                    J_diag.append(np.sum(psi_p**2 * weight_J))
                print(
                    f"[ProlateActive] Analytic core J ({spec['atom']}:"
                    f"{n_q}{'spdf'[l_q]}, {n_elec}e): "
                    f"diag = [{', '.join(f'{v:.4f}' for v in J_diag)}]"
                )

            # --- Exchange: numerical (requires core on grid) ---
            core_idx = self._locked_orbitals.index(spec)
            if core_idx < len(self._core_orbitals):
                core_orb = self._core_orbitals[core_idx]
                for p in range(n):
                    for q in range(p, n):
                        K_cpq = compute_vee_integral(
                            core_orb, self._active_orbitals[p],
                            self._active_orbitals[q], core_orb,
                        )
                        h_eff[p, q] -= K_cpq
                        if p != q:
                            h_eff[q, p] -= K_cpq

    # ------------------------------------------------------------------
    # CI solver
    # ------------------------------------------------------------------

    def solve(
        self,
        E_locked: Optional[float] = None,
        V_NN: float = 0.0,
        n_states: int = 1,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build and diagonalize the active-space CI Hamiltonian.

        Parameters
        ----------
        E_locked : float, optional
            Override for locked core energy. If None, uses self.E_locked
            computed from core_spec prolate MOs.
        V_NN : float
            Nuclear repulsion energy.
        n_states : int
            Number of eigenvalues to return.

        Returns
        -------
        eigvals : np.ndarray
            Total energies: E_active + E_locked + V_NN.
        eigvecs : np.ndarray
            CI coefficient vectors.
        """
        if self._h_eff is None:
            self.build()

        if E_locked is None:
            E_locked = self.E_locked

        t0 = time.perf_counter()

        # Enumerate Slater determinants
        all_spinorbs = list(range(self.n_spinorb))
        sd_basis = list(combinations(all_spinorbs, self.n_electrons))
        n_sd = len(sd_basis)
        sd_index = {sd: i for i, sd in enumerate(sd_basis)}

        h_eff = self._h_eff
        eri = self._eri

        H_CI = np.zeros((n_sd, n_sd))

        # --- Diagonal elements ---
        for I, sd_I in enumerate(sd_basis):
            val = 0.0
            for p in sd_I:
                val += h_eff[p >> 1, p >> 1]
            n_occ = len(sd_I)
            for i in range(n_occ):
                pi = sd_I[i]
                for j in range(i + 1, n_occ):
                    pj = sd_I[j]
                    val += eri[pi >> 1, pj >> 1, pi >> 1, pj >> 1]
                    if (pi & 1) == (pj & 1):
                        val -= eri[pi >> 1, pj >> 1, pj >> 1, pi >> 1]
            H_CI[I, I] = val

        # --- Off-diagonal: singles ---
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = frozenset(sd_I)
            n_occ = len(sd_I)

            for kp in range(n_occ):
                p = sd_I[kp]
                sig_p = p & 1

                for r in all_spinorbs:
                    if r in occ_I or (r & 1) != sig_p:
                        continue

                    me = h_eff[p >> 1, r >> 1]
                    for q in sd_I:
                        if q == p:
                            continue
                        me += eri[p >> 1, q >> 1, r >> 1, q >> 1]
                        if (q & 1) == sig_p:
                            me -= eri[p >> 1, q >> 1, q >> 1, r >> 1]

                    if abs(me) < 1e-12:
                        continue

                    new_sd = list(sd_I)
                    new_sd[kp] = r
                    new_sd_t = tuple(sorted(new_sd))
                    J_idx = sd_index.get(new_sd_t)
                    if J_idx is None or J_idx <= I:
                        continue

                    phase = self._phase_single(sd_I, kp, r)
                    H_CI[I, J_idx] += phase * me
                    H_CI[J_idx, I] += phase * me

        # --- Off-diagonal: doubles ---
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = frozenset(sd_I)
            n_occ = len(sd_I)

            for kp in range(n_occ):
                p = sd_I[kp]
                sig_p = p & 1
                for kq in range(kp + 1, n_occ):
                    q = sd_I[kq]
                    sig_q = q & 1
                    for ir in range(len(all_spinorbs)):
                        r = all_spinorbs[ir]
                        if r in occ_I:
                            continue
                        for js in range(ir + 1, len(all_spinorbs)):
                            s = all_spinorbs[js]
                            if s in occ_I:
                                continue
                            if sig_p + sig_q != (r & 1) + (s & 1):
                                continue

                            me = 0.0
                            if sig_p == (r & 1) and sig_q == (s & 1):
                                me += eri[p >> 1, q >> 1, r >> 1, s >> 1]
                            if sig_p == (s & 1) and sig_q == (r & 1):
                                me -= eri[p >> 1, q >> 1, s >> 1, r >> 1]

                            if abs(me) < 1e-12:
                                continue

                            new_sd = list(sd_I)
                            new_sd[kp] = r
                            new_sd[kq] = s
                            new_sd_t = tuple(sorted(new_sd))
                            J_idx = sd_index.get(new_sd_t)
                            if J_idx is None or J_idx <= I:
                                continue

                            phase = self._phase_double(sd_I, kp, kq, r, s)
                            H_CI[I, J_idx] += phase * me
                            H_CI[J_idx, I] += phase * me

        # Diagonalize
        eigvals, eigvecs = eigh(H_CI)
        eigvals = eigvals[:n_states]
        eigvecs = eigvecs[:, :n_states]

        # Add locked energy and nuclear repulsion
        eigvals = eigvals + E_locked + V_NN

        dt = time.perf_counter() - t0
        if self.verbose:
            E_active = eigvals[0] - E_locked - V_NN
            print(
                f"[ProlateActive] CI: {n_sd} SDs, "
                f"E_active = {E_active:.6f}, "
                f"E_locked = {E_locked:.6f}, "
                f"V_NN = {V_NN:.6f}, "
                f"E_total = {eigvals[0]:.6f} Ha, "
                f"time = {dt:.3f}s"
            )

        self._sd_basis = sd_basis
        self._H_CI = H_CI
        self._n_sd = n_sd

        return eigvals, eigvecs

    # ------------------------------------------------------------------
    # Phase computation (Slater-Condon)
    # ------------------------------------------------------------------

    @staticmethod
    def _phase_single(sd: Tuple[int, ...], kp: int, r: int) -> float:
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _phase_double(
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
