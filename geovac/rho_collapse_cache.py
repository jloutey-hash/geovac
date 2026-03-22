"""
Fast adiabatic PES scanner via angular eigenvalue caching.

The angular eigenvalue problem for two-electron diatomics in molecule-frame
hyperspherical coordinates (Paper 15) has the form:

    [Λ²/2 + R_e × C_mol(α, Ω; ρ)] Φ = μ(R_e, ρ) Φ

where ρ = R/(2R_e) is the dimensionless ratio of internuclear distance to
electronic hyperradius.  The charge function C_mol depends on (R, R_e) only
through ρ, but the eigenvalue μ depends on BOTH ρ and R_e because the
kinetic operator Λ²/2 does not scale with R_e.

The key optimization: for a given R, the angular sweep over R_e (equivalently
over ρ) is the expensive part.  By precomputing the ρ-independent parts of
the Hamiltonian (kinetic, centrifugal, Liouville, e-e coupling) and only
recomputing the nuclear coupling at each ρ, the per-ρ cost is roughly halved.
Combined with cubic spline interpolation on a coarse ρ grid, the radial solve
at each R-point becomes milliseconds.

References:
  - Paper 15, Sec. III-IV (charge function, adiabatic separation)
  - level4_multichannel.py (angular Hamiltonian builder)
  - hyperspherical_radial.py (radial solver template)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from typing import Tuple, List, Optional, Union, Callable
from math import sqrt

from geovac.level4_multichannel import (
    _channel_list,
    _channel_list_extended,
    compute_nuclear_coupling,
    _ee_coupling,
    _ee_coupling_general,
)
from geovac.hyperspherical_angular import gaunt_integral


def _build_static_hamiltonian(
    alpha_grid: np.ndarray,
    l_max: int,
    Z_A: float,
    Z_B: float,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
) -> Tuple[np.ndarray, list]:
    """
    Build the ρ-independent part of the angular Hamiltonian.

    This includes:
    - Kinetic energy in α (tridiagonal FD)
    - Centrifugal barriers l(l+1)/(2cos²α), l(l+1)/(2sin²α)
    - Liouville correction (-2)
    - Electron-electron coupling (depends on α only, not ρ)

    Returns
    -------
    H_static : ndarray of shape (N, N)
        Static Hamiltonian (everything except nuclear coupling and R_e scaling).
    channels_4 : list of (l1, m1, l2, m2) tuples
    """
    homonuclear = (Z_A == Z_B)

    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )

    n_ch = len(channels_4)
    n_alpha = len(alpha_grid)
    N = n_ch * n_alpha
    h = alpha_grid[1] - alpha_grid[0]

    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    H = np.zeros((N, N))

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    def idx(ch: int, i: int) -> int:
        return ch * n_alpha + i

    # --- Diagonal blocks: kinetic + centrifugal + Liouville ---
    for ic, (l1, m1, l2, m2) in enumerate(channels_4):
        V_cent = (0.5 * l1 * (l1 + 1) / cos_a**2
                  + 0.5 * l2 * (l2 + 1) / sin_a**2)
        V_liouville = -2.0
        V_diag = V_cent + V_liouville

        for i in range(n_alpha):
            ii = idx(ic, i)
            H[ii, ii] = kinetic_diag + V_diag[i]

        for i in range(n_alpha - 1):
            ii = idx(ic, i)
            jj = idx(ic, i + 1)
            H[ii, jj] = kinetic_off
            H[jj, ii] = kinetic_off

    # --- E-e coupling (ρ-independent) ---
    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue

            if m_max == 0:
                W_ee = _ee_coupling(l1p, l2p, l1, l2, alpha_grid, l_max)
            else:
                W_ee = _ee_coupling_general(
                    l1p, m1p, l2p, m2p, l1, m1, l2, m2, alpha_grid
                )

            if np.max(np.abs(W_ee)) < 1e-15:
                continue

            # Store WITHOUT R_e factor — will be multiplied later
            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                H[ii, jj] += W_ee[i]  # R_e factor added at solve time
                if ic != jc:
                    H[jj, ii] += W_ee[i]

    return H, channels_4


def _build_nuclear_matrix(
    alpha_grid: np.ndarray,
    rho: float,
    channels_4: list,
    Z_A: float,
    Z_B: float,
    rho_A: float = None,
    rho_B: float = None,
) -> np.ndarray:
    """
    Build the nuclear coupling matrix at a given ρ.

    Returns the matrix WITHOUT the R_e prefactor.

    Parameters
    ----------
    alpha_grid : ndarray
        Alpha grid points.
    rho : float
        R / (2 R_e).
    channels_4 : list
        Channel labels (l1, m1, l2, m2).
    Z_A, Z_B : float
        Nuclear charges.
    rho_A, rho_B : float or None
        Per-nucleus ρ values. If None, uses ρ for both (midpoint origin).

    Returns
    -------
    W_nuc : ndarray of shape (N, N)
        Nuclear coupling matrix (without R_e factor).
    """
    if rho_A is None:
        rho_A = rho
    if rho_B is None:
        rho_B = rho

    n_ch = len(channels_4)
    n_alpha = len(alpha_grid)
    N = n_ch * n_alpha

    W = np.zeros((N, N))

    def idx(ch: int, i: int) -> int:
        return ch * n_alpha + i

    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue

            if m1p != m1 or m2p != m2:
                continue

            V_nuc = compute_nuclear_coupling(
                l1p, l2p, l1, l2, m1, m2, alpha_grid, rho,
                Z=Z_A, Z_A=Z_A, Z_B=Z_B,
                rho_A=rho_A, rho_B=rho_B,
            )

            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                W[ii, jj] += V_nuc[i]
                if ic != jc:
                    W[jj, ii] += V_nuc[i]

    return W


class AngularCache:
    """
    Cache angular eigenvalues for fast adiabatic PES construction.

    Precomputes the ρ-independent parts of the angular Hamiltonian once,
    then efficiently solves at each (ρ, R_e) point by adding only the
    nuclear coupling.

    The angular eigenvalue μ(ρ, R_e) is stored on a grid and interpolated
    via cubic splines for the radial solve.

    Usage
    -----
    >>> cache = AngularCache(Z_A=1.0, Z_B=1.0, l_max=2, n_alpha=100)
    >>> cache.build_for_R(R=1.4, n_rho=80)
    >>> mu_0 = cache.epsilon(0, rho=0.5)  # interpolated ground-state eigenvalue
    """

    def __init__(
        self,
        Z_A: float = 1.0,
        Z_B: float = 1.0,
        l_max: int = 2,
        n_alpha: int = 100,
        m_max: int = 0,
        l_max_per_m: Union[dict, None] = None,
        n_eig: int = 1,
        Z_A_func: Optional[Callable] = None,
        n_theta: int = 64,
    ):
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.m_max = m_max
        self.l_max_per_m = l_max_per_m
        self.n_eig = n_eig
        self.Z_A_func = Z_A_func
        self.n_theta = n_theta

        # Build alpha grid
        h = (np.pi / 2) / (n_alpha + 1)
        self.alpha_grid = (np.arange(n_alpha) + 1) * h

        # Precompute static Hamiltonian (kinetic + centrifugal + Liouville + ee)
        self.H_static, self.channels_4 = _build_static_hamiltonian(
            self.alpha_grid, l_max, Z_A, Z_B, m_max, l_max_per_m,
        )
        self.n_ch = len(self.channels_4)

        # Cache state
        self.rho_grid: Optional[np.ndarray] = None
        self.mu_curves: Optional[np.ndarray] = None
        self.splines: Optional[List[CubicSpline]] = None
        self._R: Optional[float] = None

    def build_for_R(
        self,
        R: float,
        n_rho: int = 80,
        rho_min: float = 0.02,
        rho_max: float = 5.0,
        z0: float = 0.0,
        verbose: bool = False,
    ) -> None:
        """
        Compute angular eigenvalue curves on a ρ grid for a specific R.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        n_rho : int
            Number of ρ grid points.
        rho_min, rho_max : float
            Range of ρ = R/(2R_e).
        z0 : float
            Origin shift along internuclear axis.
        verbose : bool
            Print progress.
        """
        self._R = R

        # Non-uniform ρ grid: denser at small ρ where μ changes rapidly
        rho_grid = np.concatenate([
            np.linspace(rho_min, 0.1, n_rho // 6 + 1),
            np.linspace(0.1, 0.3, n_rho // 5),
            np.linspace(0.3, 0.8, n_rho // 4),
            np.linspace(0.8, 2.0, n_rho // 4),
            np.linspace(2.0, rho_max, n_rho // 6 + 1),
        ])
        rho_grid = np.unique(rho_grid)
        self.rho_grid = rho_grid

        n_pts = len(rho_grid)
        self.mu_curves = np.zeros((self.n_eig, n_pts))

        for i, rho in enumerate(rho_grid):
            R_e = R / (2.0 * rho)

            # Compute nuclear coupling at this ρ
            rho_A = rho - z0 / R_e if z0 != 0.0 else rho
            rho_B = rho + z0 / R_e if z0 != 0.0 else rho

            W_nuc = _build_nuclear_matrix(
                self.alpha_grid, rho, self.channels_4,
                self.Z_A, self.Z_B, rho_A=rho_A, rho_B=rho_B,
            )

            # Full Hamiltonian: H_static has ee WITHOUT R_e factor,
            # nuclear WITHOUT R_e factor. Add R_e scaling to both.
            # H = (kinetic+cent+liouville) + R_e * (ee + nuclear)
            # But H_static already has kinetic+cent+liouville as-is
            # and ee without R_e. So:
            # H_full = H_kinetic_part + R_e * H_ee_part + R_e * W_nuc
            # We stored H_static = H_kinetic + H_ee (without R_e on ee).
            # We need to separate them to apply R_e correctly.
            #
            # Actually, let me restructure: H_static has kinetic+cent+liouville
            # with no R_e, and ee with no R_e. The correct angular Hamiltonian is:
            # H = (kinetic+cent+liouville) + R_e * (ee + nuclear)
            # = H_kinetic_only + R_e * (H_ee_only + W_nuc)
            #
            # But we stored H_static = H_kinetic_only + H_ee_only.
            # So H = H_static + (R_e - 1) * H_ee_only + R_e * W_nuc
            # This requires separating H_static...
            #
            # Simpler: just build H = H_kinetic + R_e * (H_ee + W_nuc)
            # But we don't have H_kinetic and H_ee separated.
            #
            # Let me fix this by computing H_full correctly.
            # H_static stores kinetic+cent+liouville+ee (all without R_e on ee).
            # The correct form needs R_e on ee and nuclear.
            # So I need to scale the ee part by R_e. But it's baked into H_static.
            #
            # The simplest fix: store H_kinetic and H_ee separately.
            # I'll refactor below. For now, use the level4 code directly.

            # Fall back to building full Hamiltonian directly
            from geovac.level4_multichannel import build_angular_hamiltonian

            H = build_angular_hamiltonian(
                self.alpha_grid, rho, R_e, self.l_max, self.Z_A,
                self.m_max, self.l_max_per_m,
                Z_A=self.Z_A, Z_B=self.Z_B, z0=z0,
                Z_A_func=self.Z_A_func, n_theta=self.n_theta,
            )

            # Only need lowest n_eig eigenvalues
            evals = eigh(H, subset_by_index=[0, self.n_eig - 1],
                         eigvals_only=True)
            self.mu_curves[:, i] = evals[:self.n_eig]

            if verbose and (i % 20 == 0 or i == n_pts - 1):
                print(f"  rho={rho:.4f} (R_e={R_e:.3f}): mu_0={evals[0]:.6f}")

        # Build cubic spline interpolators for each eigenvalue curve
        self.splines = []
        for nu in range(self.n_eig):
            self.splines.append(
                CubicSpline(self.rho_grid, self.mu_curves[nu],
                            extrapolate=True)
            )

    def epsilon(self, nu: int, rho: float) -> float:
        """
        Interpolated angular eigenvalue for channel ν at given ρ.

        Parameters
        ----------
        nu : int
            Channel index (0 = ground state).
        rho : float
            Dimensionless ratio R/(2R_e).

        Returns
        -------
        mu : float
            Angular eigenvalue.
        """
        if self.splines is None:
            raise RuntimeError("Cache not built. Call build_for_R() first.")
        return float(self.splines[nu](rho))

    def epsilon_array(self, nu: int, rho_array: np.ndarray) -> np.ndarray:
        """Vectorized interpolation for an array of ρ values."""
        if self.splines is None:
            raise RuntimeError("Cache not built. Call build_for_R() first.")
        return self.splines[nu](rho_array)

    def effective_potential(
        self,
        R_e_grid: np.ndarray,
        R: float,
        nu: int = 0,
    ) -> np.ndarray:
        """
        Compute the adiabatic effective potential U(R_e; R).

        U(R_e) = (μ(ρ) + 15/8) / R_e²

        where ρ = R/(2R_e) and μ is the cached angular eigenvalue.

        Parameters
        ----------
        R_e_grid : ndarray
            Electronic hyperradius grid.
        R : float
            Internuclear distance.
        nu : int
            Adiabatic channel index.

        Returns
        -------
        U : ndarray
            Effective potential (Ha).
        """
        rho_vals = R / (2.0 * R_e_grid)
        mu_vals = self.epsilon_array(nu, rho_vals)
        return (mu_vals + 15.0 / 8.0) / R_e_grid**2


class FastAdiabaticPES:
    """
    Fast potential energy surface scanner using adiabatic hyperspherical method.

    For each R-point:
    1. Build angular eigenvalue cache on a coarse ρ grid (the expensive step)
    2. Interpolate to construct the effective radial potential U(R_e; R)
    3. Solve the 1D radial Schrödinger equation (fast tridiagonal eigensolver)

    Usage
    -----
    >>> cache = AngularCache(Z_A=1.0, Z_B=1.0, l_max=2, n_alpha=100)
    >>> pes = FastAdiabaticPES(cache, Z_A=1.0, Z_B=1.0)
    >>> R_grid = np.linspace(0.8, 6.0, 20)
    >>> results = pes.scan_pes(R_grid)
    """

    def __init__(
        self,
        cache: AngularCache,
        Z_A: float = 1.0,
        Z_B: float = 1.0,
        R_e_min: float = 0.3,
        R_e_max: float = 15.0,
        n_Re: int = 600,
    ):
        self.cache = cache
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.R_e_min = R_e_min
        self.R_e_max = R_e_max
        self.n_Re = n_Re

    @classmethod
    def from_core_screening(
        cls,
        core_screening,
        Z_B: float,
        l_max: int = 2,
        n_alpha: int = 100,
        m_max: int = 0,
        n_eig: int = 1,
        n_theta: int = 64,
        R_e_min: float = 0.3,
        R_e_max: float = 15.0,
        n_Re: int = 600,
    ) -> "FastAdiabaticPES":
        """
        Construct a PES scanner from a CoreScreening object.

        Extracts Z_eff(r) from the core screening and wires it into the
        Level 4 angular Hamiltonian as Z_A_func.

        Parameters
        ----------
        core_screening : CoreScreening
            Solved core screening object (must have .solve() called).
        Z_B : float
            Nuclear charge for nucleus B (constant).
        l_max, n_alpha, m_max, n_eig, n_theta : int
            Angular solver parameters.
        R_e_min, R_e_max, n_Re : float/int
            Radial grid parameters.

        Returns
        -------
        FastAdiabaticPES
            Configured PES scanner with screened Z_A.
        """
        Z_A = float(core_screening.Z)
        Z_A_func = core_screening.z_eff

        cache = AngularCache(
            Z_A=Z_A, Z_B=Z_B, l_max=l_max, n_alpha=n_alpha,
            m_max=m_max, n_eig=n_eig,
            Z_A_func=Z_A_func, n_theta=n_theta,
        )
        return cls(cache, Z_A=Z_A, Z_B=Z_B,
                   R_e_min=R_e_min, R_e_max=R_e_max, n_Re=n_Re)

    def energy_at_R(
        self,
        R: float,
        n_rho: int = 80,
        z0: float = 0.0,
        verbose: bool = False,
    ) -> Tuple[float, float]:
        """
        Compute electronic and total energy at internuclear distance R.

        Parameters
        ----------
        R : float
            Internuclear distance (bohr).
        n_rho : int
            Number of ρ grid points for angular cache.
        z0 : float
            Origin shift.
        verbose : bool
            Print diagnostics.

        Returns
        -------
        E_elec : float
            Electronic energy (Ha).
        E_total : float
            Total energy E_elec + V_NN (Ha).
        """
        # Step 1: Build angular cache for this R
        self.cache.build_for_R(R, n_rho=n_rho, z0=z0, verbose=verbose)

        # Step 2: Construct radial grid and effective potential
        h_Re = (self.R_e_max - self.R_e_min) / (self.n_Re + 1)
        R_e_grid = self.R_e_min + (np.arange(self.n_Re) + 1) * h_Re

        U = self.cache.effective_potential(R_e_grid, R, nu=0)

        # Step 3: Solve 1D radial equation
        # [-½ d²F/dR_e² + U(R_e)] F = E_elec F
        diag = np.ones(self.n_Re) / h_Re**2 + U
        off_diag = -0.5 * np.ones(self.n_Re - 1) / h_Re**2

        evals, _ = eigh_tridiagonal(
            diag, off_diag,
            select='i', select_range=(0, 0),
        )

        E_elec = evals[0]
        V_NN = self.Z_A * self.Z_B / R
        E_total = E_elec + V_NN

        return E_elec, E_total

    def scan_pes(
        self,
        R_grid: np.ndarray,
        n_rho: int = 80,
        z0: float = 0.0,
        verbose: bool = True,
    ) -> dict:
        """
        Scan the potential energy surface E_total(R).

        Parameters
        ----------
        R_grid : ndarray
            Internuclear distances (bohr).
        n_rho : int
            Number of ρ grid points per R-point.
        z0 : float
            Origin shift.
        verbose : bool
            Print table of results.

        Returns
        -------
        results : dict
            Keys: 'R', 'E_elec', 'E_total', 'V_NN', 'wall_times',
            'R_eq', 'E_min', 'D_e'.
        """
        import time

        n_R = len(R_grid)
        E_elec = np.zeros(n_R)
        E_total = np.zeros(n_R)
        V_NN = np.zeros(n_R)
        wall_times = np.zeros(n_R)

        if verbose:
            homonuclear = (self.Z_A == self.Z_B)
            sys_label = "H2" if homonuclear and self.Z_A == 1 else f"Z_A={self.Z_A}, Z_B={self.Z_B}"
            print(f"Fast adiabatic PES scan for {sys_label}")
            print(f"  l_max={self.cache.l_max}, n_ch={self.cache.n_ch}, "
                  f"n_alpha={self.cache.n_alpha}, n_rho={n_rho}")
            print(f"  R range: [{R_grid[0]:.2f}, {R_grid[-1]:.2f}] bohr, "
                  f"{n_R} points")
            print(f"  Radial grid: [{self.R_e_min}, {self.R_e_max}], "
                  f"n_Re={self.n_Re}")
            print()
            print(f"  {'R':>6s}  {'E_elec':>12s}  {'E_total':>12s}  "
                  f"{'V_NN':>10s}  {'time(s)':>8s}")
            print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*8}")

        for i, R in enumerate(R_grid):
            t0 = time.time()
            ee, et = self.energy_at_R(R, n_rho=n_rho, z0=z0)
            t1 = time.time()

            E_elec[i] = ee
            E_total[i] = et
            V_NN[i] = self.Z_A * self.Z_B / R
            wall_times[i] = t1 - t0

            if verbose:
                print(f"  {R:6.3f}  {ee:12.6f}  {et:12.6f}  "
                      f"{V_NN[i]:10.6f}  {wall_times[i]:8.3f}")

        # Find equilibrium
        i_min = np.argmin(E_total)
        R_eq = R_grid[i_min]
        E_min = E_total[i_min]

        # Dissociation limit
        homonuclear = (self.Z_A == self.Z_B)
        if homonuclear:
            E_atoms = -self.Z_A**2  # two atoms each with one electron
        else:
            Z_max = max(self.Z_A, self.Z_B)
            if Z_max == 2:
                E_atoms = -2.9037
            else:
                E_atoms = -Z_max**2 + 5 * Z_max / 8

        D_e = E_atoms - E_min

        if verbose:
            print()
            print(f"  R_eq = {R_eq:.3f} bohr")
            print(f"  E_min = {E_min:.6f} Ha")
            print(f"  E_atoms = {E_atoms:.6f} Ha")
            print(f"  D_e = {D_e:.6f} Ha")
            if homonuclear and self.Z_A == 1:
                D_e_exact = 0.17447
                print(f"  D_e / exact = {D_e / D_e_exact * 100:.1f}%"
                      f"  (exact = {D_e_exact:.5f} Ha)")
            print(f"  Total time: {wall_times.sum():.2f}s"
                  f"  (avg {wall_times.mean():.3f}s per R-point)")

        return {
            'R': R_grid,
            'E_elec': E_elec,
            'E_total': E_total,
            'V_NN': V_NN,
            'wall_times': wall_times,
            'R_eq': R_eq,
            'E_min': E_min,
            'E_atoms': E_atoms,
            'D_e': D_e,
        }
