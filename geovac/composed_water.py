"""
Composed natural geometry solver for H₂O (symmetric stretch).

Extends the composed geometry framework to a 10-electron triatomic with:
  - 1 core block:  O 1s² (2 electrons, Level 3 hyperspherical, Z=8)
  - 2 bond pairs:  O–H_A and O–H_B (2 electrons each, Level 4, Z_A=Z_eff, Z_B=1)
  - 2 lone pairs:  on O (2 electrons each, Level 3, Z_eff)

Block decomposition (7 blocks, matching build_composed_h2o in composed_qubit.py):
  Block 1: Core (Z=8, 2e)
  Block 2: Bond pair 1 O-side (Z_eff=6)
  Block 3: Bond pair 1 H-side (Z=1)
  Block 4: Bond pair 2 O-side (Z_eff=6)
  Block 5: Bond pair 2 H-side (Z=1)
  Block 6: Lone pair 1 (Z_eff=6)
  Block 7: Lone pair 2 (Z_eff=6)

Total energy at each R (symmetric stretch, R_OH1 = R_OH2 = R):
  E_total(R) = E_core
              + 2 × E_bond(R)
              + 2 × E_lone_pair          [R-independent]
              + V_NN(R)                  [nuclear repulsion]
              + V_core_H(R)             [core–ligand, 2 copies]
              + E_coupling(R)           [inter-fiber exchange, 6 pairs]

Inter-fiber coupling (6 pairs):
  Bond 1 ↔ Bond 2:       angle = θ_HOH
  Bond 1 ↔ Lone pair 1:  angle = θ_bond_lone  (×2 by symmetry with Bond 2)
  Bond 1 ↔ Lone pair 2:  angle = θ_bond_lone  (×2 by symmetry with Bond 2)
  Lone pair 1 ↔ Lone pair 2: angle = θ_lone_lone

Nuclear repulsion for H₂O symmetric stretch with angle θ:
  V_NN = 2 × Z_O × Z_H / R + Z_H² / d_HH
  where d_HH = 2R sin(θ/2)
  But bond-pair solver already includes Z_A–Z_B nuclear attraction, so
  V_NN here is only the H–H repulsion that is NOT in the bond-pair Hamiltonian.

References:
  - Paper 17: Composed natural geometries
  - Paper 14 Sec IV: Composed qubit Hamiltonians (H₂O 7-block structure)
  - Paper 13: Hyperspherical lattice (core & lone pair solver)
"""

import time
import numpy as np
from scipy.optimize import curve_fit
from typing import Optional, Dict, Any, List

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.nuclear_lattice import HARTREE_TO_CM, AMU_TO_ME
from geovac.lone_pair import (
    solve_lone_pair,
    extract_channel_data_level3,
    compute_inter_fiber_overlap,
)
from geovac.inter_fiber_coupling import (
    full_exchange_inter_fiber_energy,
    extract_channel_data,
)


# Reference data for H₂O
REFERENCE_DATA: Dict[str, float] = {
    'R_eq': 1.809,            # bohr (0.957 Å)
    'theta_HOH': 104.5,       # degrees
    'omega_e_sym': 3657.0,    # cm-1 (ν₁ symmetric stretch)
    'D_e': 0.370,             # Ha (~10.1 eV, H₂O → O + 2H)
    'E_total_expt': -76.438,  # Ha (experimental total energy)
}


def _morse_potential(R: np.ndarray, E_min: float, D_e: float,
                     a: float, R_eq: float) -> np.ndarray:
    """Morse potential for curve fitting."""
    return E_min + D_e * (1.0 - np.exp(-a * (R - R_eq)))**2


def _v_cross_nuc_1s(Z_core: float, n_core: int, Z_other: float,
                    R: float) -> float:
    """
    Core-to-other-nucleus attraction for 1s^n_core core.

    V_cross = -n_core * Z_other * <1s_Z|1/r_other|1s_Z>
    where <1s_Z|1/r_B|1s_Z> = (1/R) * [1 - (1 + Z*R) * exp(-2*Z*R)]
    """
    zr = Z_core * R
    expectation = (1.0 / R) * (1.0 - (1.0 + zr) * np.exp(-2.0 * zr))
    return -n_core * Z_other * expectation


class ComposedWaterSolver:
    """
    Composed natural geometry solver for H₂O (symmetric stretch).

    Parameters
    ----------
    Z_O : int
        Oxygen nuclear charge (default 8).
    Z_H : int
        Hydrogen nuclear charge (default 1).
    n_core : int
        Number of core electrons on O (default 2).
    l_max : int
        Maximum angular momentum in Level 4 solver.
    n_alpha : int
        Number of alpha grid points.
    theta_HOH : float
        H-O-H bond angle in radians (default 1.8238 ≈ 104.5°).
    theta_bond_lone : float
        Angle between bond axis and lone pair axis (radians).
        Default 2.0944 ≈ 120° (tetrahedral-ish approximation).
    theta_lone_lone : float
        Angle between the two lone pair axes (radians).
        Default 1.9897 ≈ 114°.
    pk_mode : str
        'ab_initio', 'manual', or 'none'.
    pk_A, pk_B : float or None
        Manual PK parameters.
    include_coupling : bool
        If True, add inter-fiber exchange coupling to the PES.
    coupling_pairs : str
        Which inter-fiber pairs to couple:
          'all'       — all 6 pairs (1 bond-bond, 4 bond-lone, 1 lone-lone)
          'bond_bond' — only bond 1 ↔ bond 2
          'none'      — no coupling (equivalent to include_coupling=False)
    coupling_scale : float
        Scaling factor for inter-fiber coupling (default 1.0).
    verbose : bool
        Print progress.
    """

    def __init__(
        self,
        Z_O: int = 8,
        Z_H: int = 1,
        n_core: int = 2,
        l_max: int = 2,
        n_alpha: int = 100,
        theta_HOH: float = 1.8238,
        theta_bond_lone: float = 2.0944,
        theta_lone_lone: float = 1.9897,
        pk_mode: str = 'ab_initio',
        pk_A: Optional[float] = None,
        pk_B: Optional[float] = None,
        include_coupling: bool = True,
        coupling_pairs: str = 'all',
        coupling_scale: float = 1.0,
        verbose: bool = True,
    ) -> None:
        if n_core != 2:
            raise ValueError("Only n_core=2 (1s² core) is supported.")
        if coupling_pairs not in ('all', 'bond_bond', 'none'):
            raise ValueError(
                f"coupling_pairs must be 'all', 'bond_bond', or 'none', "
                f"got '{coupling_pairs}'")

        self.Z_O = float(Z_O)
        self.Z_H = float(Z_H)
        self.n_core = n_core
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.verbose = verbose

        # Geometry angles
        self.theta_HOH = theta_HOH
        self.theta_bond_lone = theta_bond_lone
        self.theta_lone_lone = theta_lone_lone

        # Effective charge
        self.Z_eff = self.Z_O - self.n_core  # 6

        # Coupling settings
        self.include_coupling = include_coupling
        self.coupling_pairs = coupling_pairs
        if coupling_pairs == 'none':
            self.include_coupling = False
        self.coupling_scale = coupling_scale

        # PK mode
        if pk_mode not in ('ab_initio', 'manual', 'none'):
            raise ValueError(f"pk_mode must be 'ab_initio', 'manual', or 'none'")
        self.pk_mode = pk_mode
        self.pk_potentials: Optional[List[dict]] = None
        self.ab_initio_pk: Optional[AbInitioPK] = None

        if pk_mode == 'manual':
            if pk_A is None:
                pk_A = 5.0 * (self.Z_O / 3.0)
            if pk_B is None:
                pk_B = 7.0 * (self.Z_O / 3.0) ** 2
            self.pk_A = pk_A
            self.pk_B = pk_B
            self.pk_potentials = [{
                'C_core': pk_A,
                'beta_core': pk_B,
                'atom': 'A',
            }]
        else:
            self.pk_A = 0.0
            self.pk_B = 0.0

        # Pipeline state
        self.core: Optional[CoreScreening] = None
        self.E_core: Optional[float] = None
        self.E_lone_pair: Optional[float] = None
        self._lone_pair_result: Optional[dict] = None
        self._lone_pair_channel_data: Optional[dict] = None
        self.pes_result: Optional[dict] = None
        self.spectro: Optional[Dict[str, float]] = None
        self.timings: Dict[str, float] = {}

        # Reference data
        self.ref = REFERENCE_DATA

    def solve_core(self) -> float:
        """
        Step 1: Solve the O 1s² core (R-independent).

        Returns
        -------
        E_core : float
            Core ground-state energy (Ha).
        """
        t0 = time.time()

        self.core = CoreScreening(
            Z=int(self.Z_O), l_max=self.l_max, n_alpha=200,
        )
        self.core.solve(verbose=self.verbose)
        self.E_core = self.core.energy

        # Derive ab initio PK from the core
        if self.pk_mode == 'ab_initio':
            self.ab_initio_pk = AbInitioPK(
                self.core, n_core=self.n_core,
            )
            self.pk_A = self.ab_initio_pk.A
            self.pk_B = self.ab_initio_pk.B
            self.pk_potentials = [{
                'C_core': self.pk_A,
                'beta_core': self.pk_B,
                'atom': 'A',
            }]

        self.timings['core'] = time.time() - t0

        if self.verbose:
            print(f"\n  [H2O] Core: E_core = {self.E_core:.6f} Ha")
            print(f"  Z_eff = {self.Z_eff:.3f}")
            if self.pk_mode != 'none':
                print(f"  PK ({self.pk_mode}): A={self.pk_A:.4f},"
                      f" B={self.pk_B:.4f}")
                if self.ab_initio_pk is not None:
                    print(f"  r_core = {self.ab_initio_pk.r_core:.4f} bohr")
            print(f"  Time: {self.timings['core']:.1f}s")

        return self.E_core

    def solve_lone_pair(self) -> float:
        """
        Step 2: Solve one lone pair (R-independent, reused for both).

        Returns
        -------
        E_lone : float
            Lone pair ground-state energy (Ha).
        """
        t0 = time.time()

        self._lone_pair_result = solve_lone_pair(
            Z_eff=self.Z_eff,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            verbose=self.verbose,
        )
        self.E_lone_pair = self._lone_pair_result['energy']

        # Extract channel data once for coupling calculations
        self._lone_pair_channel_data = extract_channel_data_level3(
            self._lone_pair_result,
            Z=self.Z_eff,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
        )

        self.timings['lone_pair'] = time.time() - t0

        if self.verbose:
            print(f"\n  [H2O] Lone pair: E = {self.E_lone_pair:.6f} Ha")
            print(f"  Time: {self.timings['lone_pair']:.1f}s")

        return self.E_lone_pair

    def _solve_bond_at_R(self, R: float, n_Re: int = 300,
                         return_full: bool = False) -> Any:
        """
        Solve a single O–H bond pair at distance R.

        Uses Level 4 multichannel with Z_A=Z_eff(O), Z_B=Z_H.
        """
        result = solve_level4_h2_multichannel(
            R=R,
            Z_A=self.Z_eff,
            Z_B=self.Z_H,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            n_Re=n_Re,
            verbose=False,
            pk_potentials=self.pk_potentials,
            origin='charge_center',
        )
        if return_full:
            return result
        return result['E_elec']

    def _nuclear_repulsion(self, R: float) -> float:
        """
        Full nuclear-nuclear repulsion for H₂O.

        The Level 4 solver returns E_elec (electronic energy only, no V_NN).
        We must add all nuclear repulsion terms:
          V_NN = 2 × Z_O × Z_H / R  (O–H repulsion, 2 bonds)
                + Z_H² / d_HH        (H–H repulsion)

        where d_HH = 2R sin(θ/2).

        The Level 4 solver uses Z_A=Z_eff (not Z_O) for the electron-nuclear
        attraction. The difference (Z_O-Z_eff)*Z_H/R = n_core*Z_H/R per bond
        is compensated by V_cross_nuc (core electron attraction to H nuclei).
        """
        d_HH = 2.0 * R * np.sin(self.theta_HOH / 2.0)
        V_OH = 2.0 * self.Z_O * self.Z_H / R
        V_NN = self.Z_H ** 2 / d_HH
        return V_OH + V_NN

    def _cross_nuclear(self, R: float) -> float:
        """
        Core-to-ligand nuclear attraction (both H atoms).

        Each H nucleus attracts the O core electrons.
        Returns total V_cross for both H atoms.
        """
        V_cross_one = _v_cross_nuc_1s(
            self.Z_O, self.n_core, self.Z_H, R)
        return 2.0 * V_cross_one

    def _compute_bond_bond_exchange(
        self, R: float, bond_result: Dict[str, Any],
    ) -> float:
        """Bond 1 ↔ Bond 2 exchange at θ_HOH."""
        exch = full_exchange_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_H,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            bond_angle=self.theta_HOH,
        )
        return exch['E_exchange']

    def _compute_bond_lone_exchange(
        self, R: float, bond_result: Dict[str, Any],
    ) -> float:
        """
        Bond ↔ Lone pair exchange at θ_bond_lone.

        Uses the generic inter-fiber overlap between Level 4 (bond)
        and Level 3 (lone pair) channel data, with S*F⁰ factorization.

        Returns exchange energy for ONE bond-lone pair.
        """
        from geovac.inter_fiber_coupling import (
            extract_origin_density_algebraic,
            slater_f0_integral,
        )

        # Extract bond channel data
        bond_channel_data = extract_channel_data(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_H,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
        )

        # Compute overlap S between bond and lone pair
        overlap = compute_inter_fiber_overlap(
            bond_channel_data,
            self._lone_pair_channel_data,
            bond_angle=self.theta_bond_lone,
        )
        S_avg = overlap['S_avg']

        # Estimate F⁰ from densities of the two fibers
        r_grid_bond, P_bond = extract_origin_density_algebraic(
            bond_channel_data, n_r=300, r_max=10.0,
        )
        r_grid_lone, P_lone = extract_origin_density_algebraic(
            self._lone_pair_channel_data, n_r=300, r_max=10.0,
        )
        F0 = slater_f0_integral(r_grid_bond, P_bond, P_lone)

        # Exchange = -S * F⁰ (attractive)
        return -abs(S_avg) * F0

    def _compute_lone_lone_exchange(self) -> float:
        """
        Lone pair 1 ↔ Lone pair 2 exchange at θ_lone_lone.

        Both fibers are identical Level 3 systems; R-independent.
        """
        from geovac.inter_fiber_coupling import (
            extract_origin_density_algebraic,
            slater_f0_integral,
        )

        overlap = compute_inter_fiber_overlap(
            self._lone_pair_channel_data,
            self._lone_pair_channel_data,
            bond_angle=self.theta_lone_lone,
        )
        S_avg = overlap['S_avg']

        # F⁰ from lone-pair density (self-interaction monopole)
        r_grid, P_r = extract_origin_density_algebraic(
            self._lone_pair_channel_data, n_r=300, r_max=10.0,
        )
        F0 = slater_f0_integral(r_grid, P_r, P_r)

        return -abs(S_avg) * F0

    def _compute_all_coupling(
        self, R: float, bond_result: Dict[str, Any],
    ) -> Dict[str, float]:
        """
        Compute inter-fiber exchange couplings.

        Which pairs are included depends on self.coupling_pairs:
          'all'       — all 6 pairs
          'bond_bond' — only bond 1 ↔ bond 2

        Returns dict with individual contributions and total.
        """
        # 1. Bond–Bond: 1 pair at θ_HOH (always included when coupling is on)
        E_bb = self._compute_bond_bond_exchange(R, bond_result)

        # 2. Bond–Lone: 4 pairs (2 bonds × 2 lone pairs)
        if self.coupling_pairs == 'all':
            E_bl_one = self._compute_bond_lone_exchange(R, bond_result)
            E_bl_total = 4.0 * E_bl_one
        else:
            E_bl_one = 0.0
            E_bl_total = 0.0

        # 3. Lone–Lone: 1 pair at θ_lone_lone (R-independent)
        if self.coupling_pairs == 'all':
            E_ll = self._compute_lone_lone_exchange()
        else:
            E_ll = 0.0

        E_total = self.coupling_scale * (E_bb + E_bl_total + E_ll)

        return {
            'E_bond_bond': E_bb,
            'E_bond_lone_one': E_bl_one,
            'E_bond_lone_total': E_bl_total,
            'E_lone_lone': E_ll,
            'E_coupling_total': E_total,
        }

    def scan_pes(self, R_grid: Optional[np.ndarray] = None,
                 n_Re: int = 300) -> Dict[str, Any]:
        """
        Step 3: Scan PES for symmetric stretch.

        E_total(R) = E_core + 2*E_bond(R) + 2*E_lone + V_NN(R)
                    + V_cross(R) + E_coupling(R)
        """
        if self.E_core is None:
            raise RuntimeError("Call solve_core() first.")
        if self.E_lone_pair is None:
            raise RuntimeError("Call solve_lone_pair() first.")

        if R_grid is None:
            R_grid = np.arange(1.4, 4.5, 0.2)

        n_R = len(R_grid)
        t0 = time.time()

        if self.verbose:
            print(f"\n  [H2O] PES scan: {n_R} R-points, n_Re={n_Re}")
            if self.include_coupling:
                print(f"  {'R':>6s}  {'E_bond':>10s}  {'V_NN':>8s}"
                      f"  {'V_cross':>8s}  {'E_coupl':>8s}"
                      f"  {'E_total':>12s}  {'time':>6s}")
                print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*8}"
                      f"  {'-'*8}  {'-'*12}  {'-'*6}")
            else:
                print(f"  {'R':>6s}  {'E_bond':>10s}  {'V_NN':>8s}"
                      f"  {'V_cross':>8s}  {'E_total':>12s}  {'time':>6s}")
                print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*8}"
                      f"  {'-'*12}  {'-'*6}")

        E_bond_arr = np.zeros(n_R)
        V_NN_arr = np.zeros(n_R)
        V_cross_arr = np.zeros(n_R)
        E_coupling_arr = np.zeros(n_R)
        E_total = np.zeros(n_R)
        wall_times = np.zeros(n_R)

        for i, R in enumerate(R_grid):
            ti = time.time()
            try:
                # Solve ONE bond pair (both are identical by symmetry)
                if self.include_coupling:
                    bond_result = self._solve_bond_at_R(
                        R, n_Re=n_Re, return_full=True)
                    E_bond_one = bond_result['E_elec']
                else:
                    E_bond_one = self._solve_bond_at_R(R, n_Re=n_Re)
                    bond_result = None
                E_bond_total = 2.0 * E_bond_one

                V_NN = self._nuclear_repulsion(R)
                V_cross = self._cross_nuclear(R)

                if self.include_coupling and bond_result is not None:
                    coupling = self._compute_all_coupling(R, bond_result)
                    E_coupl = coupling['E_coupling_total']
                else:
                    E_coupl = 0.0

                E_bond_arr[i] = E_bond_total
                V_NN_arr[i] = V_NN
                V_cross_arr[i] = V_cross
                E_coupling_arr[i] = E_coupl
                E_total[i] = (self.E_core
                              + E_bond_total
                              + 2.0 * self.E_lone_pair
                              + V_NN + V_cross + E_coupl)

            except Exception as e:
                if self.verbose:
                    print(f"  {R:6.3f}  FAILED: {e}")
                E_total[i] = np.nan

            wall_times[i] = time.time() - ti

            if self.verbose and not np.isnan(E_total[i]):
                if self.include_coupling:
                    print(f"  {R:6.3f}  {E_bond_total:10.4f}  {V_NN:8.4f}"
                          f"  {V_cross:8.4f}  {E_coupl:8.4f}"
                          f"  {E_total[i]:12.6f}  {wall_times[i]:6.1f}s")
                else:
                    print(f"  {R:6.3f}  {E_bond_total:10.4f}  {V_NN:8.4f}"
                          f"  {V_cross:8.4f}  {E_total[i]:12.6f}"
                          f"  {wall_times[i]:6.1f}s")

        self.timings['pes_scan'] = time.time() - t0

        # Remove failed points
        valid = ~np.isnan(E_total)
        R_valid = R_grid[valid]
        E_valid = E_total[valid]

        if len(R_valid) < 3:
            raise RuntimeError("Too few valid PES points for analysis.")

        # Find minimum
        i_min = np.argmin(E_valid)
        R_eq = R_valid[i_min]
        E_min = E_valid[i_min]

        # Dissociation limit: largest-R point
        E_dissoc = E_valid[-1]
        D_e = E_dissoc - E_min

        self.pes_result = {
            'R': R_grid.tolist(),
            'E_bond': E_bond_arr.tolist(),
            'V_NN': V_NN_arr.tolist(),
            'V_cross_nuc': V_cross_arr.tolist(),
            'E_coupling': E_coupling_arr.tolist(),
            'E_total': E_total.tolist(),
            'wall_times': wall_times.tolist(),
            'R_eq': float(R_eq),
            'E_min': float(E_min),
            'E_dissoc': float(E_dissoc),
            'D_e': float(D_e),
            'R_valid': R_valid.tolist(),
            'E_valid': E_valid.tolist(),
        }

        if self.verbose:
            ref_R = self.ref.get('R_eq', '?')
            ref_D = self.ref.get('D_e', '?')
            print(f"\n  R_eq = {R_eq:.3f} bohr  (ref: {ref_R})")
            print(f"  E_min = {E_min:.6f} Ha")
            print(f"  E_dissoc = {E_dissoc:.6f} Ha")
            print(f"  D_e = {D_e:.6f} Ha  (ref: {ref_D})")
            print(f"  PES time: {self.timings['pes_scan']:.1f}s"
                  f"  (avg {wall_times[valid].mean():.1f}s/pt)")

        return self.pes_result

    def fit_spectroscopic_constants(
        self, fit_window: float = 1.0,
    ) -> Dict[str, float]:
        """
        Step 4: Fit Morse potential near minimum.

        Extracts R_eq, D_e, ω_e (symmetric stretch frequency).
        """
        if self.pes_result is None:
            raise RuntimeError("Call scan_pes() first.")

        t0 = time.time()

        R_valid = np.array(self.pes_result['R_valid'])
        E_valid = np.array(self.pes_result['E_valid'])
        R_eq_grid = self.pes_result['R_eq']
        E_min_grid = self.pes_result['E_min']

        # Select points near minimum
        mask = np.abs(R_valid - R_eq_grid) < fit_window
        R_fit = R_valid[mask]
        E_fit = E_valid[mask]

        if len(R_fit) < 4:
            R_fit = R_valid
            E_fit = E_valid

        D_e_guess = max(self.pes_result['D_e'], 0.01)

        try:
            popt, _ = curve_fit(
                _morse_potential, R_fit, E_fit,
                p0=[E_min_grid, D_e_guess, 1.0, R_eq_grid],
                maxfev=10000,
            )
            E_min_fit, D_e_fit, a_fit, R_eq_fit = popt
        except RuntimeError:
            E_min_fit = E_min_grid
            D_e_fit = D_e_guess
            a_fit = 1.0
            R_eq_fit = R_eq_grid
            if self.verbose:
                print("  WARNING: Morse fit failed, using grid values.")

        D_e_fit = abs(D_e_fit)
        a_fit = abs(a_fit)

        # For symmetric stretch of bent triatomic:
        # mu_eff = M_H for symmetric stretch when O is approximately stationary
        M_H = 1.00782503
        M_O = 15.999
        mu_sym_amu = M_H
        mu_sym_au = mu_sym_amu * AMU_TO_ME

        omega_e_au = a_fit * np.sqrt(2.0 * D_e_fit / mu_sym_au)
        omega_e_cm = omega_e_au * HARTREE_TO_CM

        self.spectro = {
            'R_eq': float(R_eq_fit),
            'D_e': float(D_e_fit),
            'a': float(a_fit),
            'E_min': float(E_min_fit),
            'omega_e_sym': float(omega_e_cm),
        }

        self.timings['morse_fit'] = time.time() - t0

        if self.verbose:
            ref_R = self.ref.get('R_eq', '---')
            ref_w = self.ref.get('omega_e_sym', '---')
            print(f"\n  [H2O] Morse fit (symmetric stretch):")
            print(f"    R_eq    = {R_eq_fit:.4f} bohr  (ref: {ref_R})")
            print(f"    D_e     = {D_e_fit:.6f} Ha")
            print(f"    a       = {a_fit:.4f} bohr^-1")
            print(f"    omega_e = {omega_e_cm:.1f} cm-1  (ref: {ref_w})")

            if isinstance(ref_R, float):
                err_R = abs(R_eq_fit - ref_R) / ref_R * 100
                print(f"    R_eq error: {err_R:.1f}%")
            if isinstance(ref_w, float):
                err_w = abs(omega_e_cm - ref_w) / ref_w * 100
                print(f"    omega_e error: {err_w:.1f}%")

        return self.spectro

    def run_all(self, R_grid: Optional[np.ndarray] = None,
                n_Re: int = 300,
                fit_window: float = 1.0) -> Dict[str, Any]:
        """Run the complete H₂O composed pipeline."""
        t_total = time.time()

        if self.verbose:
            print("=" * 64)
            print("H₂O Composed Water Pipeline")
            print(f"  Z_O={self.Z_O:.0f}, Z_H={self.Z_H:.0f},"
                  f" n_core={self.n_core}, l_max={self.l_max}")
            print(f"  Z_eff={self.Z_eff:.3f}, pk_mode={self.pk_mode}")
            print(f"  θ_HOH={np.degrees(self.theta_HOH):.1f}°,"
                  f" θ_bond-lone={np.degrees(self.theta_bond_lone):.1f}°,"
                  f" θ_lone-lone={np.degrees(self.theta_lone_lone):.1f}°")
            coupl_str = "ON" if self.include_coupling else "OFF"
            pairs_str = self.coupling_pairs if self.include_coupling else "none"
            print(f"  Inter-fiber coupling: {coupl_str}"
                  f" (pairs={pairs_str}, scale={self.coupling_scale:.2f})")
            print("=" * 64)

        self.solve_core()
        self.solve_lone_pair()
        self.scan_pes(R_grid=R_grid, n_Re=n_Re)
        self.fit_spectroscopic_constants(fit_window=fit_window)

        self.timings['total'] = time.time() - t_total

        if self.verbose:
            self._print_summary()

        return {
            'E_core': self.E_core,
            'E_lone_pair': self.E_lone_pair,
            'pes': self.pes_result,
            'spectro': self.spectro,
            'timings': self.timings,
        }

    def _print_summary(self) -> None:
        """Print comprehensive results table."""
        print("\n" + "=" * 64)
        print("H₂O Composed Water Results")
        print("=" * 64)

        print("\nPipeline timing:")
        for step in ['core', 'lone_pair', 'pes_scan', 'morse_fit']:
            t = self.timings.get(step, 0.0)
            print(f"  {step:20s}  {t:8.1f} sec")
        print(f"  {'TOTAL':20s}  {self.timings.get('total', 0.0):8.1f} sec")

        print(f"\nCore:")
        print(f"  E_core = {self.E_core:.6f} Ha")
        print(f"  E_lone_pair = {self.E_lone_pair:.6f} Ha (×2)")

        ref_R = self.ref.get('R_eq', '---')
        ref_w = self.ref.get('omega_e_sym', '---')
        ref_D = self.ref.get('D_e', '---')
        print(f"\nPES:")
        print(f"  R_eq  = {self.pes_result['R_eq']:.3f} bohr  (ref: {ref_R})")
        print(f"  E_min = {self.pes_result['E_min']:.6f} Ha")
        print(f"  D_e   = {self.pes_result['D_e']:.6f} Ha  (ref: {ref_D})")

        s = self.spectro
        if s:
            print(f"\nSpectroscopic constants:")
            print(f"  {'':14s} {'Computed':>10s} {'Reference':>10s} {'Unit':>6s}")
            print(f"  {'R_eq':14s} {s['R_eq']:10.3f} {str(ref_R):>10s} {'bohr':>6s}")
            print(f"  {'omega_e (sym)':14s} {s['omega_e_sym']:10.1f}"
                  f" {str(ref_w):>10s} {'cm-1':>6s}")
            print(f"  {'D_e':14s} {s['D_e']:10.4f} {str(ref_D):>10s} {'Ha':>6s}")

            if isinstance(ref_R, float):
                err_R = abs(s['R_eq'] - ref_R) / ref_R * 100
                print(f"\n  R_eq error: {err_R:.1f}%")
