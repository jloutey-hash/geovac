"""
General 2-core + 2-valence diatomic solver via composed graph.

Generalizes the LiH pipeline (lih_composed.py) to any diatomic with:
  - 2 core electrons on atom A (solved via hyperspherical)
  - 2 valence electrons in Level 4 molecule-frame hyperspherical
  - Nuclear motion via Morse/rotational lattice from fitted PES

Architecture:
  Layer 0 (nuclear): Morse SU(2) x SO(3) -- from fitted PES
  Layer 1 (core):    Atom A^{n_core+} 1s^2 via hyperspherical (Z=Z_A)
  Layer 2 (valence): 2 valence e- in Level 4 with Z_A_eff + PK pseudopotential

The total energy at each R is:
  E_total(R) = E_core + V_cross_nuc(R) + E_elec(R) + V_NN_bare(R)

Phillips-Kleinman pseudopotential parameters scale with Z_A:
  B_core ~ Z_A^2 (core orbital width scales as 1/Z_A)
  A_core ~ Z_A (barrier height scales linearly)
"""

import time
import numpy as np
from scipy.optimize import curve_fit
from typing import Optional, Dict, Any, List

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.nuclear_lattice import NuclearLattice, HARTREE_TO_CM, AMU_TO_ME


# Reference data for known molecules (used for comparison, not computation)
REFERENCE_DATA: Dict[str, Dict[str, float]] = {
    'LiH': {
        'R_eq': 3.015,
        'D_e': 0.0920,
        'omega_e': 1405.65,
        'omega_e_xe': 23.20,
        'B_e': 7.5131,
        'alpha_e': 0.2132,
        'E_core': -7.2799,  # exact Li+ 1s^2 energy
    },
    'BeH+': {
        'R_eq': 2.479,      # ~1.312 Å
        'omega_e': 2222.0,   # cm-1 (NIST/theoretical)
        'B_e': 10.31,        # cm-1 (estimate)
        'E_core': -13.6556,  # exact Be2+ 1s^2 energy
    },
}


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
    """
    zr = Z_core * R
    expectation = (1.0 / R) * (1.0 - (1.0 + zr) * np.exp(-2.0 * zr))
    return -n_core * Z_other * expectation


class ComposedDiatomicSolver:
    """
    General 2-core + 2-valence diatomic solver via composed graph.

    Parameters
    ----------
    Z_A : int
        Nuclear charge of the heavy atom (the one with the core).
    Z_B : int
        Nuclear charge of the light atom.
    n_core : int
        Number of core electrons (must be 2 for current implementation).
    M_A : float
        Atomic mass of atom A in amu.
    M_B : float
        Atomic mass of atom B in amu.
    l_max : int
        Maximum angular momentum in Level 4 solver.
    m_max : int
        Maximum |m| per electron. 0 = sigma-only, 1 = sigma+pi.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit. E.g., {0: 4, 1: 2} uses l_max=4
        for sigma and l_max=2 for pi. If None, uses l_max for all m values.
    n_alpha : int
        Number of alpha grid points.
    use_pauli : bool
        Whether to apply Phillips-Kleinman pseudopotential.
    pk_mode : str
        'ab_initio' - derive PK from core wavefunction (default).
        'manual' - use user-specified (pk_A, pk_B).
        'none' - disable PK entirely.
    pk_A : float or None
        PK barrier height parameter (used when pk_mode='manual').
    pk_B : float or None
        PK barrier width parameter (used when pk_mode='manual').
    pk_channel_mode : str
        'channel_blind' (default) - PK applied uniformly to all channels.
        'l_dependent' - PK weight per electron is δ_{l_i, 0}, so only
        channels with l=0 character see the PK barrier.  This prevents
        unphysical repulsion in higher-l channels that are automatically
        orthogonal to the 1s² core by angular momentum.
    zeff_mode : str or float
        'screened' (Z_A_eff = Z_A - n_core), or a float value.
    verbose : bool
        Print progress information.
    label : str
        Molecule label for display (e.g. 'LiH', 'BeH+').
    """

    def __init__(
        self,
        Z_A: int,
        Z_B: int,
        n_core: int = 2,
        M_A: float = 7.016,
        M_B: float = 1.008,
        l_max: int = 2,
        m_max: int = 0,
        l_max_per_m: Optional[Dict[int, int]] = None,
        n_alpha: int = 100,
        use_pauli: bool = True,
        pk_mode: str = 'ab_initio',
        pk_A: Optional[float] = None,
        pk_B: Optional[float] = None,
        pk_channel_mode: str = 'channel_blind',
        zeff_mode: str = 'screened',
        level4_method: str = 'adiabatic',
        verbose: bool = True,
        label: str = '',
    ) -> None:
        if n_core != 2:
            raise ValueError("Only n_core=2 (1s^2 core) is supported.")

        self.Z_A_bare = float(Z_A)
        self.Z_B = float(Z_B)
        self.n_core = n_core
        self.M_A = M_A
        self.M_B = M_B
        self.l_max = l_max
        self.m_max = m_max
        self.l_max_per_m = l_max_per_m
        self.n_alpha = n_alpha
        self.verbose = verbose
        self.label = label or f"Z_A={Z_A},Z_B={Z_B}"

        # Reduced mass
        self.mu_amu = M_A * M_B / (M_A + M_B)
        self.mu_au = self.mu_amu * AMU_TO_ME

        # Z_A_eff for valence solver
        if zeff_mode == 'screened':
            self.Z_A_eff = self.Z_A_bare - self.n_core
        else:
            self.Z_A_eff = float(zeff_mode)
        self.zeff_mode = zeff_mode

        # PK pseudopotential mode
        if pk_mode not in ('ab_initio', 'manual', 'none', 'algebraic'):
            raise ValueError(
                f"pk_mode must be 'ab_initio', 'manual', 'none', or "
                f"'algebraic', got '{pk_mode}'"
            )
        if pk_channel_mode not in ('channel_blind', 'l_dependent'):
            raise ValueError(
                f"pk_channel_mode must be 'channel_blind' or 'l_dependent',"
                f" got '{pk_channel_mode}'"
            )
        self.pk_mode = pk_mode
        self.pk_channel_mode = pk_channel_mode

        # Level 4 solver method
        if level4_method not in ('adiabatic', 'variational_2d'):
            raise ValueError(
                f"level4_method must be 'adiabatic' or 'variational_2d', "
                f"got '{level4_method}'"
            )
        self.level4_method = level4_method

        # Handle backward compatibility: use_pauli=False overrides pk_mode
        if not use_pauli:
            self.pk_mode = 'none'
        self.use_pauli = (self.pk_mode != 'none')

        self.pk_potentials: Optional[List[dict]] = None
        self.pk_projector: Optional[dict] = None
        self.ab_initio_pk: Optional[AbInitioPK] = None

        if self.pk_mode == 'manual':
            # Manual PK: use provided or Z-scaled defaults
            if pk_A is None:
                pk_A = 5.0 * (self.Z_A_bare / 3.0)
            if pk_B is None:
                pk_B = 7.0 * (self.Z_A_bare / 3.0) ** 2
            self.pk_A = pk_A
            self.pk_B = pk_B
            self.pk_potentials = [{
                'C_core': pk_A,
                'beta_core': pk_B,
                'atom': 'A',
                'channel_mode': pk_channel_mode,
            }]
        elif self.pk_mode in ('ab_initio', 'algebraic'):
            # Will be set after solve_core()
            self.pk_A = 0.0
            self.pk_B = 0.0
        else:  # 'none'
            self.pk_A = 0.0
            self.pk_B = 0.0

        # Pipeline state
        self.core: Optional[CoreScreening] = None
        self.E_core: Optional[float] = None
        self.pes_result: Optional[dict] = None
        self.spectro: Optional[Dict[str, float]] = None
        self.nuclear: Optional[NuclearLattice] = None
        self.timings: Dict[str, float] = {}

        # Reference data (if available)
        self.ref = REFERENCE_DATA.get(label, {})

    @classmethod
    def LiH(cls, l_max: int = 2, **kwargs) -> "ComposedDiatomicSolver":
        """LiH: Li (Z=3) + H (Z=1), 2 core electrons on Li.

        Uses fitted PK parameters by default for backward compatibility.
        Use pk_mode='ab_initio' for zero-parameter derivation.
        """
        defaults = dict(
            Z_A=3, Z_B=1, n_core=2,
            M_A=7.016003, M_B=1.00782503,
            label='LiH',
            pk_mode='manual', pk_A=5.0, pk_B=7.0,
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    @classmethod
    def LiH_ab_initio(cls, l_max: int = 2, **kwargs) -> "ComposedDiatomicSolver":
        """LiH with ab initio PK (zero fitted parameters)."""
        defaults = dict(
            Z_A=3, Z_B=1, n_core=2,
            M_A=7.016003, M_B=1.00782503,
            label='LiH',
            pk_mode='ab_initio',
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    @classmethod
    def LiH_algebraic_pk(cls, l_max: int = 2,
                          **kwargs) -> "ComposedDiatomicSolver":
        """LiH with algebraic PK projector (rank-1, channel-space).

        Replaces the Gaussian barrier with the exact Phillips-Kleinman
        projector V_PK = E_shift * |core⟩⟨core| expressed in the Level 4
        channel basis.  Zero fitted parameters.
        """
        defaults = dict(
            Z_A=3, Z_B=1, n_core=2,
            M_A=7.016003, M_B=1.00782503,
            label='LiH',
            pk_mode='algebraic',
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    @classmethod
    def BeH_plus(cls, l_max: int = 2, **kwargs) -> "ComposedDiatomicSolver":
        """BeH+: Be (Z=4) + H (Z=1), 2 core electrons on Be.

        Uses Z-scaled manual PK by default for backward compatibility.
        Use pk_mode='ab_initio' for zero-parameter derivation.
        """
        defaults = dict(
            Z_A=4, Z_B=1, n_core=2,
            M_A=9.01218, M_B=1.00782503,
            label='BeH+',
            pk_mode='manual',  # Z-scaled defaults from __init__
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    @classmethod
    def BeH_plus_ab_initio(cls, l_max: int = 2,
                           **kwargs) -> "ComposedDiatomicSolver":
        """BeH+ with ab initio PK (zero fitted parameters)."""
        defaults = dict(
            Z_A=4, Z_B=1, n_core=2,
            M_A=9.01218, M_B=1.00782503,
            label='BeH+',
            pk_mode='ab_initio',
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    def solve_core(self) -> float:
        """
        Step 1: Solve the core (Z=Z_A, 2 electrons) for E_core.

        If pk_mode='ab_initio', derives PK parameters from the core solution.

        Returns
        -------
        E_core : float
            Core ground-state energy (Ha).
        """
        t0 = time.time()

        self.core = CoreScreening(
            Z=int(self.Z_A_bare), l_max=self.l_max, n_alpha=200,
        )
        self.core.solve(verbose=self.verbose)
        self.E_core = self.core.energy

        # Derive ab initio PK parameters from the core solution
        if self.pk_mode == 'ab_initio':
            self.ab_initio_pk = AbInitioPK(
                self.core, n_core=self.n_core,
            )
            self.pk_A = self.ab_initio_pk.A
            self.pk_B = self.ab_initio_pk.B
            pk_d = self.ab_initio_pk.pk_dict(atom='A')
            pk_d['channel_mode'] = self.pk_channel_mode
            self.pk_potentials = [pk_d]
        elif self.pk_mode == 'algebraic':
            self.ab_initio_pk = AbInitioPK(
                self.core, n_core=self.n_core,
            )
            self.pk_A = self.ab_initio_pk.A
            self.pk_B = self.ab_initio_pk.B
            self.pk_projector = self.ab_initio_pk.algebraic_projector(atom='A')
            # No Gaussian pk_potentials — the projector replaces them

        self.timings['core'] = time.time() - t0

        if self.verbose:
            ref_E = self.ref.get('E_core')
            if ref_E:
                err = abs(self.E_core - ref_E) / abs(ref_E) * 100
                print(f"\n  [{self.label}] Core: E_core = {self.E_core:.6f} Ha"
                      f"  (ref: {ref_E:.4f}, err: {err:.2f}%)")
            else:
                print(f"\n  [{self.label}] Core: E_core = {self.E_core:.6f} Ha")
            print(f"  Z_A_eff = {self.Z_A_eff:.3f} ({self.zeff_mode})")
            if self.use_pauli:
                mode_str = f" ({self.pk_mode})"
                ch_str = (f", channel_mode={self.pk_channel_mode}"
                          if self.pk_channel_mode != 'channel_blind' else "")
                if self.pk_mode == 'algebraic':
                    print(f"  PK: algebraic projector,"
                          f" E_shift={self.pk_projector['energy_shift']:.4f} Ha")
                else:
                    print(f"  PK: A={self.pk_A:.4f}, B={self.pk_B:.4f}"
                          f"{mode_str}{ch_str}")
                if self.ab_initio_pk is not None:
                    print(f"  r_core = {self.ab_initio_pk.r_core:.4f} bohr")
            else:
                print(f"  PK: disabled (pk_mode='none')")
            print(f"  Time: {self.timings['core']:.1f}s")

        return self.E_core

    def _solve_valence_at_R(self, R: float, n_Re: int = 300) -> float:
        """Solve Level 4 valence problem at a single R."""
        # n_coupled=-1 activates the direct 2D variational solver,
        # which bypasses the adiabatic approximation entirely.
        n_coupled = -1 if self.level4_method == 'variational_2d' else 1
        result = solve_level4_h2_multichannel(
            R=R,
            Z_A=self.Z_A_eff,
            Z_B=self.Z_B,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            n_Re=n_Re,
            verbose=False,
            m_max=self.m_max,
            l_max_per_m=self.l_max_per_m,
            pk_potentials=self.pk_potentials,
            pk_projector=self.pk_projector,
            n_coupled=n_coupled,
        )
        return result['E_elec']

    def scan_pes(self, R_grid: Optional[np.ndarray] = None,
                 n_Re: int = 300) -> Dict[str, Any]:
        """
        Step 2: Scan PES.

        E_composed(R) = E_core + V_cross_nuc(R) + E_elec(R) + V_NN_bare(R)
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
            print(f"\n  [{self.label}] PES scan: {n_R} R-points, n_Re={n_Re}")
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
        self, fit_window: float = 1.5,
    ) -> Dict[str, float]:
        """
        Step 3: Fit Morse potential to PES near minimum.

        Parameters
        ----------
        fit_window : float
            Half-width in bohr around R_eq for fitting.
        """
        if self.pes_result is None:
            raise RuntimeError("Call scan_pes() first.")

        t0 = time.time()

        R_valid = self.pes_result['R_valid']
        E_valid = self.pes_result['E_valid']
        R_eq_grid = self.pes_result['R_eq']
        E_min_grid = self.pes_result['E_min']

        # Select points near minimum for Morse fit
        mask = np.abs(R_valid - R_eq_grid) < fit_window
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

        mu = self.mu_au
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
            ref_R = self.ref.get('R_eq', '---')
            ref_w = self.ref.get('omega_e', '---')
            ref_B = self.ref.get('B_e', '---')
            print(f"\n  [{self.label}] Morse fit:")
            print(f"    R_eq   = {R_eq_fit:.4f} bohr  (ref: {ref_R})")
            print(f"    D_e    = {D_e_fit:.6f} Ha")
            print(f"    a      = {a_fit:.4f} bohr^-1")
            print(f"    w_e    = {omega_e_cm:.1f} cm-1  (ref: {ref_w})")
            print(f"    w_exe  = {omega_e_xe_cm:.1f} cm-1")
            print(f"    B_e    = {B_e_cm:.2f} cm-1  (ref: {ref_B})")

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
            print(f"\n  [{self.label}] Nuclear lattice:"
                  f" v_max={nv-1}, J_max={J_max}")
            print(f"    Total states: {self.nuclear.n_states}")
            print(f"    Fundamental v01 = {nu01_ha * HARTREE_TO_CM:.1f} cm-1")

        return self.nuclear

    def run_all(self, R_grid: Optional[np.ndarray] = None,
                n_Re: int = 300, J_max: int = 10,
                fit_window: float = 1.5) -> Dict[str, Any]:
        """Run the complete composed graph pipeline."""
        t_total = time.time()

        if self.verbose:
            print("=" * 64)
            print(f"{self.label} Composed Graph Pipeline")
            m_str = f", m_max={self.m_max}" if self.m_max > 0 else ""
            print(f"  Z_A={self.Z_A_bare:.0f}, Z_B={self.Z_B:.0f},"
                  f" n_core={self.n_core}, l_max={self.l_max}{m_str}")
            print("=" * 64)

        self.solve_core()
        self.scan_pes(R_grid=R_grid, n_Re=n_Re)
        self.fit_spectroscopic_constants(fit_window=fit_window)
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
        print(f"{self.label} Composed Graph Results")
        m_str = f", m_max={self.m_max}" if self.m_max > 0 else ""
        print(f"  Z_A={self.Z_A_bare:.0f}, Z_B={self.Z_B:.0f},"
              f" l_max={self.l_max}{m_str}")
        print("=" * 64)

        print("\nPipeline timing:")
        for step in ['core', 'pes_scan', 'morse_fit', 'nuclear']:
            t = self.timings.get(step, 0.0)
            print(f"  {step:20s}  {t:8.1f} sec")
        print(f"  {'TOTAL':20s}  {self.timings.get('total', 0.0):8.1f} sec")

        ref_E = self.ref.get('E_core')
        if ref_E:
            err = abs(self.E_core - ref_E) / abs(ref_E) * 100
            print(f"\nCore:")
            print(f"  E_core = {self.E_core:.6f} Ha"
                  f"  (ref: {ref_E:.4f}, err: {err:.2f}%)")
        else:
            print(f"\nCore:")
            print(f"  E_core = {self.E_core:.6f} Ha")

        ref_R = self.ref.get('R_eq', '---')
        ref_D = self.ref.get('D_e', '---')
        print(f"\nPES:")
        print(f"  R_eq  = {self.pes_result['R_eq']:.3f} bohr  (ref: {ref_R})")
        print(f"  E_min = {self.pes_result['E_min']:.6f} Ha")
        print(f"  D_e   = {self.pes_result['D_e']:.6f} Ha  (ref: {ref_D})")

        print(f"\nSpectroscopic constants:")
        ref_w = self.ref.get('omega_e', '---')
        ref_B = self.ref.get('B_e', '---')
        s = self.spectro
        hdr = f"  {'':14s} {'Computed':>10s} {'Reference':>10s} {'Unit':>6s}"
        print(hdr)
        print(f"  {'R_eq':14s} {s['R_eq']:10.3f} {str(ref_R):>10s} {'bohr':>6s}")
        print(f"  {'omega_e':14s} {s['omega_e']:10.1f} {str(ref_w):>10s} {'cm-1':>6s}")
        print(f"  {'omega_e_xe':14s} {s['omega_e_xe']:10.1f} {'---':>10s} {'cm-1':>6s}")
        print(f"  {'B_e':14s} {s['B_e']:10.2f} {str(ref_B):>10s} {'cm-1':>6s}")
        print(f"  {'D_e':14s} {s['D_e']:10.4f} {str(ref_D):>10s} {'Ha':>6s}")

        if self.use_pauli:
            ch_str = (f", channel_mode={self.pk_channel_mode}"
                      if self.pk_channel_mode != 'channel_blind' else "")
            print(f"\nPK pseudopotential ({self.pk_mode}{ch_str}):"
                  f" A={self.pk_A:.4f}, B={self.pk_B:.4f}")
            if self.ab_initio_pk is not None:
                print(f"  r_core = {self.ab_initio_pk.r_core:.4f} bohr")
        else:
            print(f"\nPK pseudopotential: disabled")

        if self.nuclear is not None:
            nv = self.nuclear.vib.n_states
            print(f"\nRovibrational spectrum (v_max={nv-1}):")
            for v in range(min(3, nv)):
                for J in [0, 1]:
                    E = self.nuclear.rovibrational_energy(v, J)
                    E_cm = E * HARTREE_TO_CM
                    print(f"  v={v}, J={J}:  {E_cm:10.1f} cm-1")

        print("=" * 64)
