"""
Composed natural geometry solver for linear triatomic molecules (BeH2).

Extends ComposedDiatomicSolver to linear A-B-A or B-A-B geometries where:
  - 2 core electrons on the central atom (solved via hyperspherical)
  - 2 bond pairs, each with 2 valence electrons in Level 4 mol-frame hyperspherical
  - For symmetric stretch, both bonds are identical (D_inf_h symmetry)

Architecture (5 blocks, matching Paper 14 Sec IV):
  Block 1: Core (Z=Z_center, 2 electrons, Level 3 hyperspherical)
  Block 2: Bond pair 1, center-side (Z_eff = Z_center - n_core)
  Block 3: Bond pair 1, ligand-side (Z = Z_ligand)
  Block 4: Bond pair 2, center-side (Z_eff, identical to Block 2 by symmetry)
  Block 5: Bond pair 2, ligand-side (Z = Z_ligand, identical to Block 3)

Total energy at each R (symmetric stretch, R1 = R2 = R):
  E_total(R) = E_core + 2*V_cross_nuc(R) + 2*E_bond(R) + V_NN(R) + V_inter(R)

V_inter(R) is the mean-field inter-bond electron-electron repulsion,
modeled as point charges at hydrogenic expectation values. This is a
scalar correction at each R and does not affect the qubit Hamiltonian.

Nuclear repulsion for linear L-C-L geometry:
  V_NN = 2*Z_C*Z_L/R + Z_L^2/(2R)

References:
  - Paper 17: Composed natural geometries (LiH, BeH+)
  - Paper 14 Sec IV: Composed qubit Hamiltonians (BeH2 block structure)
  - Paper 13: Hyperspherical lattice (core solver)
"""

import json
import time
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
from typing import Optional, Dict, Any, List

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.nuclear_lattice import HARTREE_TO_CM, AMU_TO_ME
from geovac.inter_fiber_coupling import (
    monopole_inter_fiber_energy,
    exchange_inter_fiber_energy,
    compute_overlap_diagnostic,
    direct_exchange_inter_fiber_energy,
    multipole_exchange_inter_fiber_energy,
    full_exchange_inter_fiber_energy,
)


# Reference data for BeH2
REFERENCE_DATA: Dict[str, Dict[str, float]] = {
    'BeH2': {
        'R_eq': 2.507,          # bohr (1.326 A)
        'omega_e_sym': 2345.0,  # cm-1 (nu_1 symmetric stretch)
        'omega_e_asym': 2435.0, # cm-1 (nu_3 asymmetric stretch)
        'D_e': 0.147,           # Ha (~4.0 eV, BeH2 -> Be + 2H)
        'E_core': -13.6556,     # Be2+ 1s^2 exact energy
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


class ComposedTriatomicSolver:
    """
    Composed natural geometry solver for linear triatomic molecules.

    Handles symmetric L-C-L geometry (e.g., H-Be-H) with:
    - Central atom C with n_core core electrons
    - Two identical ligand atoms L
    - Two bond pairs, each solved with Level 4 multichannel

    Parameters
    ----------
    Z_center : int
        Nuclear charge of the central atom.
    Z_ligand : int
        Nuclear charge of each ligand atom.
    n_core : int
        Number of core electrons on the central atom.
    M_center : float
        Atomic mass of center atom in amu.
    M_ligand : float
        Atomic mass of ligand atoms in amu.
    l_max : int
        Maximum angular momentum in Level 4 solver.
    n_alpha : int
        Number of alpha grid points.
    pk_mode : str
        'ab_initio', 'manual', or 'none'.
    pk_A, pk_B : float or None
        Manual PK parameters.
    include_interbond : bool
        If True, add inter-bond electron-electron repulsion to the PES.
        Only affects triatomics with n_bonds > 1. Default True.
    interbond_mode : str
        'classical' - point-charge model (original, default).
        'monopole' - k=0 Slater F^0 from Level 4 wavefunction densities.
        'exchange' - S_avg * F^0 factorized exchange.
        'direct_exchange' - channel-resolved exchange without S*F^0 factorization.
        'multipole_exchange' - generalized multipole exchange up to k_max.
        'full_exchange' - off-diagonal 1-RDM exchange contraction.
        'full_exchange_kinetic' - full exchange + Löwdin kinetic orthogonalization.
    interbond_scale : float
        Scaling factor for the inter-bond repulsion. Default 1.0 (unscaled).
    k_max : int
        Maximum multipole order for interbond_mode='multipole_exchange'.
        Default 2 (monopole + dipole + quadrupole).
    bond_angle : float
        Angle between bond axes (radians). Default π (linear L-C-L).
        For bent molecules like H₂O: 1.824 rad (104.5°).
        Affects the inter-fiber exchange coupling via Legendre polynomial
        rotation phases P_l(cos θ) in place of the linear (-1)^l phases.
    verbose : bool
        Print progress.
    label : str
        Molecule label for display.
    """

    def __init__(
        self,
        Z_center: int,
        Z_ligand: int,
        n_core: int = 2,
        M_center: float = 9.01218,
        M_ligand: float = 1.00782503,
        l_max: int = 2,
        n_alpha: int = 100,
        pk_mode: str = 'ab_initio',
        pk_A: Optional[float] = None,
        pk_B: Optional[float] = None,
        include_interbond: bool = True,
        interbond_mode: str = 'classical',
        interbond_scale: float = 1.0,
        k_max: int = 2,
        bond_angle: float = np.pi,
        verbose: bool = True,
        label: str = '',
    ) -> None:
        if n_core != 2:
            raise ValueError("Only n_core=2 (1s^2 core) is supported.")

        self.Z_center = float(Z_center)
        self.Z_ligand = float(Z_ligand)
        self.n_core = n_core
        self.M_center = M_center
        self.M_ligand = M_ligand
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.verbose = verbose
        self.label = label or f"Z_C={Z_center},Z_L={Z_ligand}"
        self.n_bonds = 2  # linear L-C-L has 2 bond pairs

        # Effective charge: full (unpartitioned)
        self.Z_eff = self.Z_center - self.n_core

        # Inter-bond repulsion settings
        self.include_interbond = include_interbond
        _valid_modes = ('classical', 'monopole', 'exchange',
                        'direct_exchange', 'multipole_exchange',
                        'full_exchange', 'full_exchange_kinetic')
        if interbond_mode not in _valid_modes:
            raise ValueError(
                f"interbond_mode must be one of {_valid_modes}, "
                f"got '{interbond_mode}'"
            )
        self.interbond_mode = interbond_mode
        self.interbond_scale = interbond_scale
        self.k_max = k_max
        self.bond_angle = bond_angle

        # PK mode
        if pk_mode not in ('ab_initio', 'manual', 'none'):
            raise ValueError(f"pk_mode must be 'ab_initio', 'manual', or 'none'")
        self.pk_mode = pk_mode
        self.pk_potentials: Optional[List[dict]] = None
        self.ab_initio_pk: Optional[AbInitioPK] = None

        if pk_mode == 'manual':
            if pk_A is None:
                pk_A = 5.0 * (self.Z_center / 3.0)
            if pk_B is None:
                pk_B = 7.0 * (self.Z_center / 3.0) ** 2
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
        self.pes_result: Optional[dict] = None
        self.spectro: Optional[Dict[str, float]] = None
        self.timings: Dict[str, float] = {}

        # Reference data
        self.ref = REFERENCE_DATA.get(label, {})

    @classmethod
    def BeH2(cls, l_max: int = 2, **kwargs) -> "ComposedTriatomicSolver":
        """BeH2: Be (Z=4) center + 2x H (Z=1) ligands, 2 core electrons."""
        defaults = dict(
            Z_center=4, Z_ligand=1, n_core=2,
            M_center=9.01218, M_ligand=1.00782503,
            label='BeH2',
            pk_mode='ab_initio',
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    @classmethod
    def BeH2_manual_pk(cls, l_max: int = 2, **kwargs) -> "ComposedTriatomicSolver":
        """BeH2 with manual Z-scaled PK parameters."""
        defaults = dict(
            Z_center=4, Z_ligand=1, n_core=2,
            M_center=9.01218, M_ligand=1.00782503,
            label='BeH2',
            pk_mode='manual',
        )
        defaults.update(kwargs)
        return cls(l_max=l_max, **defaults)

    def solve_core(self) -> float:
        """
        Step 1: Solve the core (Z=Z_center, 2 electrons).

        The core is R-independent — compute once and cache.

        Returns
        -------
        E_core : float
            Core ground-state energy (Ha).
        """
        t0 = time.time()

        self.core = CoreScreening(
            Z=int(self.Z_center), l_max=self.l_max, n_alpha=200,
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
            ref_E = self.ref.get('E_core')
            if ref_E:
                err = abs(self.E_core - ref_E) / abs(ref_E) * 100
                print(f"\n  [{self.label}] Core: E_core = {self.E_core:.6f} Ha"
                      f"  (ref: {ref_E:.4f}, err: {err:.2f}%)")
            else:
                print(f"\n  [{self.label}] Core: E_core = {self.E_core:.6f} Ha")
            print(f"  Z_eff = {self.Z_eff:.3f}")
            if self.pk_mode != 'none':
                print(f"  PK ({self.pk_mode}): A={self.pk_A:.4f},"
                      f" B={self.pk_B:.4f}")
                if self.ab_initio_pk is not None:
                    print(f"  r_core = {self.ab_initio_pk.r_core:.4f} bohr")
            if self.include_interbond:
                print(f"  Inter-bond repulsion: {self.interbond_mode}"
                      f" (scale={self.interbond_scale:.2f})")
            print(f"  Time: {self.timings['core']:.1f}s")

        return self.E_core

    def _solve_bond_at_R(self, R: float, n_Re: int = 300,
                         return_full: bool = False) -> Any:
        """
        Solve a single bond pair at internuclear distance R.

        Uses Level 4 multichannel solver with Z_A=Z_eff, Z_B=Z_ligand.
        Heteronuclear (Z_eff != Z_ligand in general).

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        n_Re : int
            Grid points for hyperradial equation.
        return_full : bool
            If True, return the full Level 4 result dict (needed for
            monopole density extraction). Default False.

        Returns
        -------
        E_bond : float (if return_full=False)
            Electronic energy of one bond pair (Ha).
        result : dict (if return_full=True)
            Full Level 4 result dict.
        """
        result = solve_level4_h2_multichannel(
            R=R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            n_Re=n_Re,
            verbose=False,
            pk_potentials=self.pk_potentials,
        )
        if return_full:
            return result
        return result['E_elec']

    def _nuclear_repulsion(self, R: float) -> float:
        """
        Nuclear repulsion for linear L-C-L geometry.

        V_NN = 2*Z_C*Z_L/R + Z_L^2/(2R)

        For BeH2: V_NN = 2*4*1/R + 1*1/(2R) = 8.5/R
        """
        return (2.0 * self.Z_center * self.Z_ligand / R
                + self.Z_ligand ** 2 / (2.0 * R))

    def _cross_nuclear(self, R: float) -> float:
        """
        Core-to-ligand nuclear attraction (both ligands).

        Each H nucleus attracts the Be core electrons.
        By symmetry, both contributions are identical.

        Returns total V_cross for both ligands.
        """
        V_cross_one = _v_cross_nuc_1s(
            self.Z_center, self.n_core, self.Z_ligand, R)
        return 2.0 * V_cross_one

    def _inter_bond_repulsion(self, R: float) -> float:
        """
        Mean-field inter-bond electron-electron repulsion (point charge model).

        In linear H-Be-H, bond pair 1 points along +z and bond pair 2 along -z.
        Each bond pair has 2 valence electrons modeled as point charges:
          - Be-side electron: at distance r_Be from Be along the bond axis
          - H-side electron: at distance R from Be (at the H nucleus)

        Geometry (z-axis):
          Bond 1: e_Be at +r_Be, e_H at +R
          Bond 2: e_Be at -r_Be, e_H at -R

        4 inter-bond e-e pairs:
          1. Be(+r_Be) <-> Be(-r_Be):  distance = 2*r_Be
          2. H(+R)     <-> H(-R):      distance = 2*R
          3. Be(+r_Be) <-> H(-R):      distance = r_Be + R
          4. H(+R)     <-> Be(-r_Be):  distance = R + r_Be

        r_Be = <r>_{1s, Z_eff} = 3/(2*Z_eff) for the hydrogenic Be-side valence.

        Returns
        -------
        V_inter : float
            Inter-bond repulsion energy (Ha), always positive.
        """
        r_Be = 1.5 / self.Z_eff  # 3/(2*Z_eff), hydrogenic <r> for 1s

        V_BeBe = 1.0 / (2.0 * r_Be)         # Be-side <-> Be-side
        V_HH = 1.0 / (2.0 * R)               # H-side <-> H-side
        V_cross = 2.0 / (R + r_Be)           # 2 cross terms (symmetric)

        return self.interbond_scale * (V_BeBe + V_HH + V_cross)

    def _monopole_coupling(self, R: float,
                           bond_result: Dict[str, Any]) -> float:
        """
        Monopole (k=0) inter-fiber coupling from Level 4 wavefunction.

        Extracts the one-electron density from the Level 4 result,
        transforms to shared-center (Be) coordinates, and computes the
        Slater F^0 integral between the two identical fiber densities.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict from _solve_bond_at_R(return_full=True).

        Returns
        -------
        V_monopole : float
            Monopole inter-fiber repulsion energy (Ha), always positive.
        """
        mono = monopole_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
        )
        return self.interbond_scale * mono['E_monopole']

    def _exchange_coupling(self, R: float,
                           bond_result: Dict[str, Any]) -> float:
        """
        Exchange inter-fiber coupling via E_exch = -S_avg(R) * J_BeSide(R).

        Uses the Be-side electron density and inter-fiber channel overlap
        to compute an exchange energy that is more negative at short R
        (where overlap is larger), pulling R_eq inward toward experiment.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict.

        Returns
        -------
        E_exchange : float
            Exchange energy (Ha), negative (attractive).
        """
        exch = exchange_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            bond_angle=self.bond_angle,
        )
        return self.interbond_scale * exch['E_exchange']

    def _direct_exchange_coupling(self, R: float,
                                  bond_result: Dict[str, Any]) -> float:
        """
        Direct exchange: E = -sum_ch (-1)^{l1+l2} * F^0_ch.

        Avoids S*F^0 factorization. Each channel's monopole F^0 is weighted
        by its own exchange parity.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict.

        Returns
        -------
        E_exchange : float
            Direct exchange energy (Ha), negative (attractive).
        """
        exch = direct_exchange_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            bond_angle=self.bond_angle,
        )
        return self.interbond_scale * exch['E_exchange']

    def _multipole_exchange_coupling(self, R: float,
                                     bond_result: Dict[str, Any]) -> float:
        """
        Multipole exchange summing E_k from k=0 to k_max.

        Includes Gaunt angular coupling between channels at each
        multipole order, with selection rules |l-l'| <= k <= l+l'.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict.

        Returns
        -------
        E_exchange : float
            Multipole exchange energy (Ha), typically negative.
        """
        exch = multipole_exchange_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            k_max=self.k_max,
            bond_angle=self.bond_angle,
        )
        return self.interbond_scale * exch['E_exchange']

    def _full_exchange_coupling(self, R: float,
                                bond_result: Dict[str, Any]) -> float:
        """
        Full exchange using off-diagonal 1-RDM contraction.

        Contracts the channel-resolved 1-RDM (including off-diagonal
        elements) with the parity transform and l1-indexed F^0 matrix.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict.

        Returns
        -------
        E_exchange : float
            Full exchange energy (Ha), negative (attractive).
        """
        exch = full_exchange_inter_fiber_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            bond_angle=self.bond_angle,
        )
        return self.interbond_scale * exch['E_exchange']

    def _full_exchange_kinetic_coupling(self, R: float,
                                        bond_result: Dict[str, Any]) -> float:
        """
        Full exchange + Löwdin kinetic orthogonalization correction.

        Combines the off-diagonal 1-RDM exchange (attractive) with the
        kinetic energy penalty from orthogonalizing overlapping fiber
        wavefunctions (repulsive).

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        bond_result : dict
            Full Level 4 result dict.

        Returns
        -------
        V_inter : float
            Combined exchange + kinetic correction (Ha).
        """
        from geovac.inter_fiber_coupling import kinetic_orthogonalization_energy

        # Full exchange (attractive, negative)
        E_exch = self._full_exchange_coupling(R, bond_result)

        # Kinetic orthogonalization (repulsive, positive)
        kin = kinetic_orthogonalization_energy(
            bond_result, R,
            Z_A=self.Z_eff,
            Z_B=self.Z_ligand,
            l_max=self.l_max,
            n_alpha=self.n_alpha,
            pk_potentials=self.pk_potentials,
            n_sample_Re=10,
            bond_angle=self.bond_angle,
        )
        delta_T = kin['delta_T']

        return E_exch + self.interbond_scale * delta_T

    def scan_pes(self, R_grid: Optional[np.ndarray] = None,
                 n_Re: int = 300) -> Dict[str, Any]:
        """
        Step 2: Scan PES for symmetric stretch.

        E_total(R) = E_core + V_cross_nuc(R) + 2*E_bond(R) + V_NN(R)

        Uses symmetry: both bond pairs are identical, so we solve
        one and double the result.
        """
        if self.E_core is None:
            raise RuntimeError("Call solve_core() first.")

        if R_grid is None:
            R_grid = np.arange(1.5, 5.1, 0.1)

        n_R = len(R_grid)
        t0 = time.time()

        if self.verbose:
            hdr_extra = "  {'V_inter':>8s}" if self.include_interbond else ""
            print(f"\n  [{self.label}] PES scan: {n_R} R-points, n_Re={n_Re}")
            if self.include_interbond:
                print(f"  {'R':>6s}  {'E_bond':>10s}  {'V_NN':>8s}  {'V_cross':>8s}"
                      f"  {'V_inter':>8s}  {'E_total':>12s}  {'time':>6s}")
                print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*8}"
                      f"  {'-'*8}  {'-'*12}  {'-'*6}")
            else:
                print(f"  {'R':>6s}  {'E_bond':>10s}  {'V_NN':>8s}  {'V_cross':>8s}"
                      f"  {'E_total':>12s}  {'time':>6s}")
                print(f"  {'-'*6}  {'-'*10}  {'-'*8}  {'-'*8}"
                      f"  {'-'*12}  {'-'*6}")

        E_bond_arr = np.zeros(n_R)
        V_NN_arr = np.zeros(n_R)
        V_cross_arr = np.zeros(n_R)
        V_inter_arr = np.zeros(n_R)
        E_total = np.zeros(n_R)
        wall_times = np.zeros(n_R)

        for i, R in enumerate(R_grid):
            ti = time.time()
            try:
                # Solve ONE bond pair, then double (symmetry)
                need_full = (self.include_interbond
                             and self.interbond_mode in (
                                 'monopole', 'exchange',
                                 'direct_exchange', 'multipole_exchange',
                                 'full_exchange', 'full_exchange_kinetic'))
                if need_full:
                    bond_result = self._solve_bond_at_R(
                        R, n_Re=n_Re, return_full=True)
                    E_bond_one = bond_result['E_elec']
                else:
                    E_bond_one = self._solve_bond_at_R(R, n_Re=n_Re)
                    bond_result = None
                E_bond_total = 2.0 * E_bond_one

                V_NN = self._nuclear_repulsion(R)
                V_cross = self._cross_nuclear(R)

                if self.include_interbond:
                    if self.interbond_mode == 'monopole':
                        V_inter = self._monopole_coupling(R, bond_result)
                    elif self.interbond_mode == 'exchange':
                        V_inter = self._exchange_coupling(R, bond_result)
                    elif self.interbond_mode == 'direct_exchange':
                        V_inter = self._direct_exchange_coupling(
                            R, bond_result)
                    elif self.interbond_mode == 'multipole_exchange':
                        V_inter = self._multipole_exchange_coupling(
                            R, bond_result)
                    elif self.interbond_mode == 'full_exchange':
                        V_inter = self._full_exchange_coupling(
                            R, bond_result)
                    elif self.interbond_mode == 'full_exchange_kinetic':
                        V_inter = self._full_exchange_kinetic_coupling(
                            R, bond_result)
                    else:
                        V_inter = self._inter_bond_repulsion(R)
                else:
                    V_inter = 0.0

                E_bond_arr[i] = E_bond_total
                V_NN_arr[i] = V_NN
                V_cross_arr[i] = V_cross
                V_inter_arr[i] = V_inter
                E_total[i] = (self.E_core + V_cross + E_bond_total
                              + V_NN + V_inter)

            except Exception as e:
                if self.verbose:
                    print(f"  {R:6.3f}  FAILED: {e}")
                E_total[i] = np.nan

            wall_times[i] = time.time() - ti

            if self.verbose and not np.isnan(E_total[i]):
                if self.include_interbond:
                    print(f"  {R:6.3f}  {E_bond_total:10.4f}  {V_NN:8.4f}"
                          f"  {V_cross:8.4f}  {V_inter:8.4f}"
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
            'V_inter': V_inter_arr.tolist(),
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
        self, fit_window: float = 1.5,
    ) -> Dict[str, float]:
        """
        Step 3: Fit Morse potential to PES near minimum.

        Extracts R_eq, D_e, omega_e (symmetric stretch frequency).
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

        # For symmetric stretch of linear triatomic, the reduced mass
        # involves both ligand atoms moving symmetrically:
        # mu_sym = M_L * (M_C + M_L) / (M_C + 2*M_L)
        # But for the Be-H bond Morse fit, each bond has its own
        # reduced mass: mu_bond = M_C * M_L / (M_C + M_L)
        # The symmetric stretch frequency involves coordinated motion
        # of both H atoms: mu_sym = M_H for symmetric stretch of
        # linear symmetric triatomic (standard result).
        # Actually: mu_sym = M_L * M_C / (M_C + 2*M_L) for A1 mode
        # Wait - for symmetric stretch of X-Y-X linear molecule,
        # only Y is stationary and both X move, so mu = M_X.
        # More precisely: mu_eff = M_L for symmetric stretch when center
        # atom is stationary by symmetry.
        # Standard result: omega_1 = sqrt(k / M_L) where k is the
        # bond force constant.
        # From Morse: k = 2 * D_e * a^2, so omega = a * sqrt(2*D_e/M_L)
        mu_sym_amu = self.M_ligand  # symmetric stretch: center stationary
        mu_sym_au = mu_sym_amu * AMU_TO_ME

        omega_e_au = a_fit * np.sqrt(2.0 * D_e_fit / mu_sym_au)
        omega_e_cm = omega_e_au * HARTREE_TO_CM

        # Also compute per-bond reduced mass spectroscopic constants
        mu_bond_amu = self.M_center * self.M_ligand / (self.M_center + self.M_ligand)
        mu_bond_au = mu_bond_amu * AMU_TO_ME

        omega_e_xe_au = omega_e_au**2 / (4.0 * D_e_fit)
        omega_e_xe_cm = omega_e_xe_au * HARTREE_TO_CM

        B_e_au = 1.0 / (2.0 * mu_bond_au * R_eq_fit**2)
        B_e_cm = B_e_au * HARTREE_TO_CM

        self.spectro = {
            'R_eq': float(R_eq_fit),
            'D_e': float(D_e_fit),
            'a': float(a_fit),
            'E_min': float(E_min_fit),
            'omega_e_sym': float(omega_e_cm),
            'omega_e_xe': float(omega_e_xe_cm),
            'B_e': float(B_e_cm),
            'mu_sym_au': float(mu_sym_au),
            'mu_bond_au': float(mu_bond_au),
        }

        self.timings['morse_fit'] = time.time() - t0

        if self.verbose:
            ref_R = self.ref.get('R_eq', '---')
            ref_w = self.ref.get('omega_e_sym', '---')
            print(f"\n  [{self.label}] Morse fit (symmetric stretch):")
            print(f"    R_eq     = {R_eq_fit:.4f} bohr  (ref: {ref_R})")
            print(f"    D_e      = {D_e_fit:.6f} Ha")
            print(f"    a        = {a_fit:.4f} bohr^-1")
            print(f"    omega_e  = {omega_e_cm:.1f} cm-1  (ref: {ref_w})")
            print(f"    omega_exe = {omega_e_xe_cm:.1f} cm-1")

            if isinstance(ref_R, float):
                err_R = abs(R_eq_fit - ref_R) / ref_R * 100
                print(f"    R_eq error: {err_R:.1f}%")
            if isinstance(ref_w, float):
                err_w = abs(omega_e_cm - ref_w) / ref_w * 100
                print(f"    omega_e error: {err_w:.1f}%")

        return self.spectro

    def run_all(self, R_grid: Optional[np.ndarray] = None,
                n_Re: int = 300,
                fit_window: float = 1.5) -> Dict[str, Any]:
        """Run the complete composed triatomic pipeline."""
        t_total = time.time()

        if self.verbose:
            print("=" * 64)
            print(f"{self.label} Composed Triatomic Pipeline")
            print(f"  Z_center={self.Z_center:.0f}, Z_ligand={self.Z_ligand:.0f},"
                  f" n_core={self.n_core}, l_max={self.l_max}")
            print(f"  Z_eff={self.Z_eff:.3f}, pk_mode={self.pk_mode}")
            if self.include_interbond:
                ib_str = (f"ON ({self.interbond_mode},"
                          f" scale={self.interbond_scale:.2f})")
            else:
                ib_str = "OFF"
            print(f"  Inter-bond repulsion: {ib_str}")
            print("=" * 64)

        self.solve_core()
        self.scan_pes(R_grid=R_grid, n_Re=n_Re)
        self.fit_spectroscopic_constants(fit_window=fit_window)

        self.timings['total'] = time.time() - t_total

        if self.verbose:
            self._print_summary()

        return {
            'E_core': self.E_core,
            'pes': self.pes_result,
            'spectro': self.spectro,
            'timings': self.timings,
        }

    def _print_summary(self) -> None:
        """Print comprehensive results table."""
        print("\n" + "=" * 64)
        print(f"{self.label} Composed Triatomic Results")
        print("=" * 64)

        print("\nPipeline timing:")
        for step in ['core', 'pes_scan', 'morse_fit']:
            t = self.timings.get(step, 0.0)
            print(f"  {step:20s}  {t:8.1f} sec")
        print(f"  {'TOTAL':20s}  {self.timings.get('total', 0.0):8.1f} sec")

        print(f"\nCore:")
        print(f"  E_core = {self.E_core:.6f} Ha")

        s = self.spectro
        ref_R = self.ref.get('R_eq', '---')
        ref_w = self.ref.get('omega_e_sym', '---')
        ref_D = self.ref.get('D_e', '---')
        print(f"\nPES:")
        print(f"  R_eq  = {self.pes_result['R_eq']:.3f} bohr  (ref: {ref_R})")
        print(f"  E_min = {self.pes_result['E_min']:.6f} Ha")
        print(f"  D_e   = {self.pes_result['D_e']:.6f} Ha  (ref: {ref_D})")

        if s:
            print(f"\nSpectroscopic constants:")
            print(f"  {'':14s} {'Computed':>10s} {'Reference':>10s} {'Unit':>6s}")
            print(f"  {'R_eq':14s} {s['R_eq']:10.3f} {str(ref_R):>10s} {'bohr':>6s}")
            print(f"  {'omega_e (sym)':14s} {s['omega_e_sym']:10.1f} {str(ref_w):>10s} {'cm-1':>6s}")
            print(f"  {'D_e':14s} {s['D_e']:10.4f} {str(ref_D):>10s} {'Ha':>6s}")

            if isinstance(ref_R, float):
                err_R = abs(s['R_eq'] - ref_R) / ref_R * 100
                print(f"\n  R_eq error: {err_R:.1f}%")
            if isinstance(ref_w, float):
                err_w = abs(s['omega_e_sym'] - ref_w) / ref_w * 100
                print(f"  omega_e error: {err_w:.1f}%")
            if isinstance(ref_D, float):
                err_D = abs(s['D_e'] - ref_D) / ref_D * 100
                print(f"  D_e error: {err_D:.1f}%")

        print(f"\nZ_eff = {self.Z_eff:.3f}")
        print(f"PK pseudopotential ({self.pk_mode}):"
              f" A={self.pk_A:.4f}, B={self.pk_B:.4f}")
        if self.include_interbond:
            if self.interbond_mode == 'monopole':
                print(f"Inter-bond coupling: monopole F^0"
                      f" (scale={self.interbond_scale:.2f})")
            elif self.interbond_mode == 'exchange':
                print(f"Inter-bond coupling: exchange S*J"
                      f" (scale={self.interbond_scale:.2f})")
            elif self.interbond_mode == 'direct_exchange':
                print(f"Inter-bond coupling: direct exchange"
                      f" (scale={self.interbond_scale:.2f})")
            elif self.interbond_mode == 'multipole_exchange':
                print(f"Inter-bond coupling: multipole exchange"
                      f" k_max={self.k_max}"
                      f" (scale={self.interbond_scale:.2f})")
            elif self.interbond_mode == 'full_exchange':
                print(f"Inter-bond coupling: full 1-RDM exchange"
                      f" (scale={self.interbond_scale:.2f})")
            elif self.interbond_mode == 'full_exchange_kinetic':
                print(f"Inter-bond coupling: full exchange + kinetic orthog"
                      f" (scale={self.interbond_scale:.2f})")
            else:
                r_Be = 1.5 / self.Z_eff
                print(f"Inter-bond coupling: classical"
                      f" (scale={self.interbond_scale:.2f},"
                      f" r_Be={r_Be:.4f} bohr)")
        else:
            print(f"Inter-bond repulsion: OFF")

        # PES shape diagnostics
        R_valid = np.array(self.pes_result['R_valid'])
        E_valid = np.array(self.pes_result['E_valid'])
        i_min = np.argmin(E_valid)

        is_bound = self.pes_result['D_e'] > 0
        has_repulsive_wall = i_min > 0
        approaches_dissoc = E_valid[-1] > E_valid[i_min]

        print(f"\nPES shape diagnostics:")
        print(f"  Bound state:       {'YES' if is_bound else 'NO'}")
        print(f"  Repulsive wall:    {'YES' if has_repulsive_wall else 'NO'}")
        print(f"  Approaches dissoc: {'YES' if approaches_dissoc else 'NO'}")
        print("=" * 64)

    def save_results(self) -> Dict[str, str]:
        """Save PES data and spectroscopic constants to debug/data/."""
        base = Path(__file__).parent.parent / 'debug'

        paths: Dict[str, str] = {}

        # PES data
        if self.pes_result is not None:
            pes_path = base / 'data' / 'beh2_composed_pes.json'
            pes_path.parent.mkdir(parents=True, exist_ok=True)
            pes_data = {
                'molecule': self.label,
                'Z_center': self.Z_center,
                'Z_ligand': self.Z_ligand,
                'Z_eff': self.Z_eff,
                'n_bonds': self.n_bonds,
                'n_core': self.n_core,
                'l_max': self.l_max,
                'pk_mode': self.pk_mode,
                'pk_A': self.pk_A,
                'pk_B': self.pk_B,
                'include_interbond': self.include_interbond,
                'interbond_scale': self.interbond_scale,
                'E_core': self.E_core,
                'pes': self.pes_result,
            }
            with open(pes_path, 'w') as f:
                json.dump(pes_data, f, indent=2)
            paths['pes'] = str(pes_path)

        # Spectroscopic constants
        if self.spectro is not None:
            spec_path = base / 'data' / 'beh2_spectroscopic_constants.json'
            spec_data = {
                'molecule': self.label,
                'computed': self.spectro,
                'reference': {
                    'R_eq': 2.507,
                    'omega_e_sym': 2345.0,
                    'D_e': 0.147,
                    'source': 'Bernath et al. 2002; Shayesteh et al. 2003',
                },
                'errors': {},
            }
            ref = self.ref
            if 'R_eq' in ref:
                spec_data['errors']['R_eq_pct'] = abs(
                    self.spectro['R_eq'] - ref['R_eq']) / ref['R_eq'] * 100
            if 'omega_e_sym' in ref:
                spec_data['errors']['omega_e_pct'] = abs(
                    self.spectro['omega_e_sym'] - ref['omega_e_sym']
                ) / ref['omega_e_sym'] * 100
            if 'D_e' in ref:
                spec_data['errors']['D_e_pct'] = abs(
                    self.spectro['D_e'] - ref['D_e']) / ref['D_e'] * 100

            with open(spec_path, 'w') as f:
                json.dump(spec_data, f, indent=2)
            paths['spectro'] = str(spec_path)

        if self.verbose:
            for name, path in paths.items():
                print(f"  Saved {name}: {path}")

        return paths

    def plot_pes(self, save_path: Optional[str] = None) -> None:
        """Generate PES plot and save to file."""
        if self.pes_result is None:
            raise RuntimeError("Call scan_pes() first.")

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            if self.verbose:
                print("  WARNING: matplotlib not available, skipping plot.")
            return

        R_valid = np.array(self.pes_result['R_valid'])
        E_valid = np.array(self.pes_result['E_valid'])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Left: full PES
        ax1.plot(R_valid, E_valid, 'b-o', markersize=4, label='E_total(R)')
        ax1.axhline(y=self.pes_result['E_dissoc'], color='r', linestyle='--',
                     alpha=0.5, label='Dissociation limit')
        ax1.axvline(x=self.pes_result['R_eq'], color='g', linestyle='--',
                     alpha=0.5, label=f'R_eq = {self.pes_result["R_eq"]:.3f}')
        ax1.axvline(x=2.507, color='orange', linestyle=':', alpha=0.5,
                     label='R_eq (expt) = 2.507')
        ax1.set_xlabel('R (bohr)')
        ax1.set_ylabel('E (Ha)')
        ax1.set_title(f'{self.label} Composed PES (symmetric stretch)')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Right: relative to dissociation (binding curve)
        E_rel = E_valid - self.pes_result['E_dissoc']
        ax2.plot(R_valid, E_rel * HARTREE_TO_CM, 'b-o', markersize=4)
        ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
        ax2.set_xlabel('R (bohr)')
        ax2.set_ylabel('E - E_dissoc (cm-1)')
        ax2.set_title(f'{self.label} Binding curve')
        ax2.grid(True, alpha=0.3)

        # Add Morse fit if available
        if self.spectro is not None:
            R_fine = np.linspace(R_valid[0], R_valid[-1], 200)
            E_morse = _morse_potential(
                R_fine, self.spectro['E_min'], self.spectro['D_e'],
                self.spectro['a'], self.spectro['R_eq'])
            ax1.plot(R_fine, E_morse, 'r-', alpha=0.5, label='Morse fit')
            ax1.legend(fontsize=8)

        plt.tight_layout()

        if save_path is None:
            save_path = str(Path(__file__).parent.parent
                           / 'debug' / 'plots' / 'beh2_composed_pes.png')

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=150)
        plt.close(fig)

        if self.verbose:
            print(f"  Plot saved: {save_path}")
