"""
LiH via the composed graph: core screening + Level 4 valence + nuclear lattice.

Architecture:
  Layer 0 (nuclear): Morse SU(2) x SO(3) -- from fitted PES
  Layer 1 (core):    Li 1s^2 via hyperspherical (Z=3, 2e) -- solved once
  Layer 2 (valence): 2 valence e- in Level 4 mol-frame hyperspherical
                     with Z_A_eff = 1.0 (screened) + Phillips-Kleinman
                     pseudopotential for core-valence Pauli repulsion

The total energy at each R is:
  E_total(R) = E_core + V_cross_nuc(R) + E_elec(R) + V_NN_bare(R)
where:
  E_core ~ -7.28 Ha               (Li+ 1s^2 energy, computed once)
  V_cross_nuc(R)                   (core electrons attracted to H nucleus)
  E_elec(R)                        (Level 4 two-electron energy, Z_A=Z_eff)
  V_NN_bare(R) = Z_A_bare * Z_B / R    (bare nuclear repulsion = 3/R)

The Phillips-Kleinman pseudopotential V_PK = C * exp(-beta * r^2) / r^2
adds a repulsive barrier at short range, preventing valence electrons from
collapsing into the core orbital.  Without it, the Level 4 solver produces
a monotonically decreasing PES (R_eq -> 0) because the constant Z_eff
approach lacks the Pauli exclusion short-range repulsion.
"""

import time
import numpy as np
from scipy.optimize import curve_fit
from typing import Optional, Dict, Any, List

from geovac.core_screening import CoreScreening
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.nuclear_lattice import NuclearLattice, HARTREE_TO_CM, AMU_TO_ME


# LiH physical constants
M_LI_AMU = 7.016003   # 7-Li atomic mass in amu
M_H_AMU = 1.00782503  # 1-H atomic mass in amu
MU_LIH_AMU = M_LI_AMU * M_H_AMU / (M_LI_AMU + M_H_AMU)
MU_LIH_AU = MU_LIH_AMU * AMU_TO_ME  # reduced mass in atomic units (m_e)

# Experimental values for comparison
EXPT = {
    'R_eq': 3.015,       # bohr
    'D_e': 0.0920,       # Ha
    'omega_e': 1405.65,  # cm-1
    'omega_e_xe': 23.20, # cm-1
    'B_e': 7.5131,       # cm-1
    'alpha_e': 0.2132,   # cm-1
    'E_Li_plus': -7.2799, # Ha (exact Li+ 1s^2 energy)
}

# Phillips-Kleinman pseudopotential defaults for Li 1s^2 core.
# C_core: amplitude of Gaussian/r^2 repulsive barrier (Ha * bohr^2).
# beta_core: Gaussian exponent (1/bohr^2), ~Z_core^2 matches 1s orbital width.
# Calibrated to reproduce R_eq ~ 3.0 bohr for LiH (expt: 3.015).
PK_C_CORE_DEFAULT = 5.0
PK_BETA_CORE_DEFAULT = 7.0


def _morse_potential(R: np.ndarray, E_min: float, D_e: float,
                     a: float, R_eq: float) -> np.ndarray:
    """Morse potential for curve fitting."""
    return E_min + D_e * (1.0 - np.exp(-a * (R - R_eq)))**2


def _v_cross_nuc_1s(Z_core: float, n_core: int, Z_other: float,
                    R: float) -> float:
    """
    Analytical core-to-other-nucleus attraction for 1s^n_core core.

    V_cross = -n_core * Z_other * <1s_Z|1/r_other|1s_Z>

    where <1s_Z|1/r_B|1s_Z> = (1/R) * [1 - (1 + Z*R) * exp(-2*Z*R)]

    At large R: V_cross -> -n_core * Z_other / R
    At R=0: V_cross -> -n_core * Z_other * 2*Z_core  (finite)
    """
    zr = Z_core * R
    expectation = (1.0 / R) * (1.0 - (1.0 + zr) * np.exp(-2.0 * zr))
    return -n_core * Z_other * expectation


class LiHComposedSolver:
    """
    LiH via the composed graph: core screening + Level 4 valence + nuclear lattice.

    Parameters
    ----------
    l_max : int
        Maximum partial wave for angular expansion (default 2 for speed).
    n_alpha : int
        Number of alpha grid points for Level 4 solver.
    zeff_mode : str or float
        'screened' (Z_eff=1.0, default), 'clementi' (Z_eff=1.279), or float.
    pk_C_core : float
        Phillips-Kleinman pseudopotential amplitude (Ha * bohr^2).
    pk_beta_core : float
        Phillips-Kleinman Gaussian exponent (1/bohr^2).
    use_pk : bool
        Whether to use the PK pseudopotential (default True).
    verbose : bool
        Print progress information.
    """

    def __init__(
        self,
        l_max: int = 2,
        n_alpha: int = 100,
        zeff_mode: str = 'screened',
        pk_C_core: float = PK_C_CORE_DEFAULT,
        pk_beta_core: float = PK_BETA_CORE_DEFAULT,
        use_pk: bool = True,
        verbose: bool = True,
    ) -> None:
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.verbose = verbose

        self.Z_A_bare = 3.0   # Li bare nuclear charge
        self.Z_B = 1.0        # H nuclear charge
        self.n_core = 2       # Li 1s^2 core electrons

        # Set Z_A_eff for valence solver
        if zeff_mode == 'screened':
            self.Z_A_eff = self.Z_A_bare - self.n_core  # 1.0
        elif zeff_mode == 'clementi':
            self.Z_A_eff = 1.279
        else:
            self.Z_A_eff = float(zeff_mode)
        self.zeff_mode = zeff_mode

        # Phillips-Kleinman pseudopotential
        self.use_pk = use_pk
        self.pk_potentials: Optional[List[dict]] = None
        if use_pk:
            self.pk_potentials = [{
                'C_core': pk_C_core,
                'beta_core': pk_beta_core,
                'atom': 'A',
            }]

        # Will be populated by pipeline steps
        self.core: Optional[CoreScreening] = None
        self.E_core: Optional[float] = None
        self.pes_result: Optional[dict] = None
        self.spectro: Optional[Dict[str, float]] = None
        self.nuclear: Optional[NuclearLattice] = None
        self.timings: Dict[str, float] = {}

    def solve_core(self) -> float:
        """
        Step 1: Solve Li+ core (Z=3, 2 electrons) for E_core.

        Returns
        -------
        E_core : float
            Li+ ground-state energy (Ha).
        """
        t0 = time.time()

        self.core = CoreScreening(Z=3, l_max=self.l_max, n_alpha=200)
        self.core.solve(verbose=self.verbose)
        self.E_core = self.core.energy

        self.timings['core'] = time.time() - t0

        if self.verbose:
            err = (abs(self.E_core - EXPT['E_Li_plus'])
                   / abs(EXPT['E_Li_plus']) * 100)
            print(f"\n  Core solve: E_core = {self.E_core:.6f} Ha"
                  f"  (exact: {EXPT['E_Li_plus']:.4f}, error: {err:.2f}%)")
            print(f"  Z_A_eff = {self.Z_A_eff:.3f} ({self.zeff_mode})")
            if self.use_pk:
                pk = self.pk_potentials[0]
                print(f"  PK: C_core={pk['C_core']:.1f},"
                      f" beta_core={pk['beta_core']:.1f}")
            print(f"  Time: {self.timings['core']:.1f}s")

        return self.E_core

    def _solve_valence_at_R(self, R: float, n_Re: int = 300) -> float:
        """
        Solve Level 4 valence problem at a single internuclear distance.

        Uses solve_level4_h2_multichannel directly, passing pk_potentials
        for core-valence Pauli repulsion.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        n_Re : int
            Number of hyperradial grid points.

        Returns
        -------
        E_elec : float
            Electronic energy of 2 valence electrons (Ha).
        """
        result = solve_level4_h2_multichannel(
            R=R,
            Z_A=self.Z_A_eff,
            Z_B=self.Z_B,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            n_Re=n_Re,
            verbose=False,
            pk_potentials=self.pk_potentials,
        )
        return result['E_elec']

    def scan_pes(self, R_grid: Optional[np.ndarray] = None,
                 n_Re: int = 300) -> Dict[str, Any]:
        """
        Step 2: Scan PES using direct Level 4 solver at each R-point.

        E_composed(R) = E_core + V_cross_nuc(R) + E_elec(R) + V_NN_bare(R)

        Parameters
        ----------
        R_grid : ndarray or None
            Internuclear distances (bohr).
        n_Re : int
            Number of hyperradial grid points per R-point.
        """
        if self.E_core is None:
            raise RuntimeError("Call solve_core() first.")

        if R_grid is None:
            R_grid = np.concatenate([
                np.linspace(2.0, 2.5, 3),
                np.linspace(2.7, 4.0, 10),
                np.linspace(4.5, 7.0, 5),
            ])

        n_R = len(R_grid)
        t0 = time.time()

        if self.verbose:
            print(f"\n  PES scan: {n_R} R-points, n_Re={n_Re}")
            print(f"  {'R':>6s}  {'E_elec':>10s}  {'V_NN_b':>8s}  {'V_cross':>8s}"
                  f"  {'E_comp':>12s}  {'time':>6s}")
            print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*12}  {'-'*6}")

        E_elec_arr = np.zeros(n_R)
        V_NN_bare_arr = np.zeros(n_R)
        V_cross_arr = np.zeros(n_R)
        E_composed = np.zeros(n_R)
        wall_times = np.zeros(n_R)

        for i, R in enumerate(R_grid):
            ti = time.time()
            try:
                E_elec = self._solve_valence_at_R(R, n_Re=n_Re)

                V_NN_bare = self.Z_A_bare * self.Z_B / R
                V_cross = _v_cross_nuc_1s(
                    self.Z_A_bare, self.n_core, self.Z_B, R)

                E_elec_arr[i] = E_elec
                V_NN_bare_arr[i] = V_NN_bare
                V_cross_arr[i] = V_cross
                E_composed[i] = (self.E_core + V_cross
                                 + E_elec + V_NN_bare)
            except Exception as e:
                if self.verbose:
                    print(f"  {R:6.3f}  FAILED: {e}")
                E_composed[i] = np.nan
            wall_times[i] = time.time() - ti

            if self.verbose and not np.isnan(E_composed[i]):
                print(f"  {R:6.3f}  {E_elec:10.4f}  {V_NN_bare:8.4f}"
                      f"  {V_cross:8.4f}  {E_composed[i]:12.6f}"
                      f"  {wall_times[i]:6.1f}s")

        self.timings['pes_scan'] = time.time() - t0

        # Remove failed points
        valid = ~np.isnan(E_composed)
        R_valid = R_grid[valid]
        E_valid = E_composed[valid]

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
            'R': R_grid,
            'E_elec': E_elec_arr,
            'V_NN_bare': V_NN_bare_arr,
            'V_cross_nuc': V_cross_arr,
            'E_composed': E_composed,
            'wall_times': wall_times,
            'R_eq': R_eq,
            'E_min': E_min,
            'E_dissoc': E_dissoc,
            'D_e': D_e,
            'R_valid': R_valid,
            'E_valid': E_valid,
        }

        if self.verbose:
            print(f"\n  R_eq = {R_eq:.3f} bohr"
                  f"  (expt: {EXPT['R_eq']:.3f}, LCAO: 2.5)")
            print(f"  E_min = {E_min:.6f} Ha")
            print(f"  E_dissoc = {E_dissoc:.6f} Ha")
            print(f"  D_e = {D_e:.6f} Ha  (expt: {EXPT['D_e']:.4f})")
            print(f"  PES time: {self.timings['pes_scan']:.1f}s"
                  f"  (avg {wall_times[valid].mean():.1f}s/pt)")

        return self.pes_result

    def fit_spectroscopic_constants(self) -> Dict[str, float]:
        """
        Step 3: Fit Morse potential to PES near minimum.

        Returns
        -------
        dict with R_eq, D_e, a, omega_e, omega_e_xe, B_e, alpha_e
        """
        if self.pes_result is None:
            raise RuntimeError("Call scan_pes() first.")

        t0 = time.time()

        R_valid = self.pes_result['R_valid']
        E_valid = self.pes_result['E_valid']
        R_eq_grid = self.pes_result['R_eq']
        E_min_grid = self.pes_result['E_min']

        # Select points near minimum for Morse fit
        mask = np.abs(R_valid - R_eq_grid) < 2.5
        R_fit = R_valid[mask]
        E_fit = E_valid[mask]

        if len(R_fit) < 4:
            R_fit = R_valid
            E_fit = E_valid

        D_e_guess = max(self.pes_result['D_e'], 0.01)
        a_guess = 1.0

        try:
            popt, _ = curve_fit(
                _morse_potential, R_fit, E_fit,
                p0=[E_min_grid, D_e_guess, a_guess, R_eq_grid],
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

        mu = MU_LIH_AU
        omega_e_au = a_fit * np.sqrt(2.0 * D_e_fit / mu)
        omega_e_cm = omega_e_au * HARTREE_TO_CM
        omega_e_xe_au = omega_e_au**2 / (4.0 * D_e_fit)
        omega_e_xe_cm = omega_e_xe_au * HARTREE_TO_CM
        B_e_au = 1.0 / (2.0 * mu * R_eq_fit**2)
        B_e_cm = B_e_au * HARTREE_TO_CM

        if omega_e_au > 0:
            alpha_e_au = (6.0 * B_e_au**2 / omega_e_au
                          * (a_fit * R_eq_fit - 1.0))
            alpha_e_cm = alpha_e_au * HARTREE_TO_CM
        else:
            alpha_e_cm = 0.0

        self.spectro = {
            'R_eq': R_eq_fit,
            'D_e': D_e_fit,
            'a': a_fit,
            'E_min': E_min_fit,
            'omega_e': omega_e_cm,
            'omega_e_xe': omega_e_xe_cm,
            'B_e': B_e_cm,
            'alpha_e': abs(alpha_e_cm),
            'mu_au': mu,
        }

        self.timings['morse_fit'] = time.time() - t0

        if self.verbose:
            print(f"\n  Morse fit:")
            print(f"    R_eq   = {R_eq_fit:.4f} bohr  (expt: {EXPT['R_eq']})")
            print(f"    D_e    = {D_e_fit:.6f} Ha  (expt: {EXPT['D_e']})")
            print(f"    a      = {a_fit:.4f} bohr^-1")
            print(f"    w_e    = {omega_e_cm:.1f} cm-1  (expt: {EXPT['omega_e']})")
            print(f"    w_exe  = {omega_e_xe_cm:.1f} cm-1"
                  f"  (expt: {EXPT['omega_e_xe']})")
            print(f"    B_e    = {B_e_cm:.2f} cm-1  (expt: {EXPT['B_e']})")

        return self.spectro

    def build_nuclear_lattice(self, J_max: int = 10) -> NuclearLattice:
        """Step 4: Build rovibrational spectrum from spectroscopic constants."""
        if self.spectro is None:
            raise RuntimeError("Call fit_spectroscopic_constants() first.")

        t0 = time.time()

        self.nuclear = NuclearLattice(
            D_e=self.spectro['D_e'],
            omega_e=self.spectro['omega_e'],
            B_e=self.spectro['B_e'],
            alpha_e=self.spectro['alpha_e'],
            omega_e_xe=self.spectro['omega_e_xe'],
            J_max=J_max,
        )

        self.timings['nuclear'] = time.time() - t0

        if self.verbose:
            nv = self.nuclear.vib.n_states
            nu01_ha = (self.nuclear.vib.morse_energy(1)
                       - self.nuclear.vib.morse_energy(0))
            print(f"\n  Nuclear lattice: v_max={nv-1}, J_max={J_max}")
            print(f"    Total states: {self.nuclear.n_states}")
            print(f"    Fundamental v01 = {nu01_ha * HARTREE_TO_CM:.1f} cm-1")

        return self.nuclear

    def run_all(self, R_grid: Optional[np.ndarray] = None,
                n_Re: int = 300, J_max: int = 10) -> Dict[str, Any]:
        """Run the complete LiH composed graph pipeline."""
        t_total = time.time()

        if self.verbose:
            print("=" * 64)
            print("LiH Composed Graph Pipeline")
            print("=" * 64)

        self.solve_core()
        self.scan_pes(R_grid=R_grid, n_Re=n_Re)
        self.fit_spectroscopic_constants()
        self.build_nuclear_lattice(J_max=J_max)

        self.timings['total'] = time.time() - t_total

        if self.verbose:
            self._print_summary()

        return {
            'E_core': self.E_core,
            'pes': self.pes_result,
            'spectro': self.spectro,
            'nuclear': self.nuclear,
            'timings': self.timings,
        }

    def _print_summary(self) -> None:
        """Print comprehensive results table."""
        print("\n" + "=" * 64)
        print("LiH Composed Graph Results")
        print("=" * 64)

        print("\nPipeline timing:")
        for step in ['core', 'pes_scan', 'morse_fit', 'nuclear']:
            t = self.timings.get(step, 0.0)
            print(f"  {step:20s}  {t:8.1f} sec")
        print(f"  {'TOTAL':20s}  {self.timings.get('total', 0.0):8.1f} sec")

        err_core = (abs(self.E_core - EXPT['E_Li_plus'])
                    / abs(EXPT['E_Li_plus']) * 100)
        print(f"\nCore:")
        print(f"  E_core = {self.E_core:.6f} Ha"
              f"  (exact: {EXPT['E_Li_plus']:.4f}, error: {err_core:.2f}%)")

        print(f"\nPES:")
        print(f"  R_eq  = {self.pes_result['R_eq']:.3f} bohr"
              f"  (expt: {EXPT['R_eq']}, LCAO: 2.5)")
        print(f"  E_min = {self.pes_result['E_min']:.6f} Ha")
        print(f"  D_e   = {self.pes_result['D_e']:.6f} Ha"
              f"  (expt: {EXPT['D_e']})")

        print(f"\nSpectroscopic constants:")
        hdr = (f"  {'':14s} {'Composed':>10s} {'Expt':>10s}"
               f" {'LCAO':>10s} {'Unit':>6s}")
        print(hdr)
        s = self.spectro
        print(f"  {'R_eq':14s} {s['R_eq']:10.3f}"
              f" {EXPT['R_eq']:10.3f} {'2.5':>10s} {'bohr':>6s}")
        print(f"  {'omega_e':14s} {s['omega_e']:10.1f}"
              f" {EXPT['omega_e']:10.1f} {'---':>10s} {'cm-1':>6s}")
        print(f"  {'omega_e_xe':14s} {s['omega_e_xe']:10.1f}"
              f" {EXPT['omega_e_xe']:10.1f} {'---':>10s} {'cm-1':>6s}")
        print(f"  {'B_e':14s} {s['B_e']:10.2f}"
              f" {EXPT['B_e']:10.2f} {'---':>10s} {'cm-1':>6s}")
        print(f"  {'D_e':14s} {s['D_e']:10.4f}"
              f" {EXPT['D_e']:10.4f} {'0.093':>10s} {'Ha':>6s}")

        if self.nuclear is not None:
            nv = self.nuclear.vib.n_states
            print(f"\nRovibrational spectrum (v_max={nv-1}):")
            for v in range(min(3, nv)):
                for J in [0, 1]:
                    E = self.nuclear.rovibrational_energy(v, J)
                    E_cm = E * HARTREE_TO_CM
                    print(f"  v={v}, J={J}:  {E_cm:10.1f} cm-1")

        print("=" * 64)
