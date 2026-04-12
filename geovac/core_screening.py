"""
Core electron screening function from hyperspherical two-electron wavefunctions.

Extracts Z_eff(r) -- the effective nuclear charge seen by an external electron
at distance r from the nucleus -- from a solved two-electron hyperspherical
system.

Two methods are available:

  histogram (original): Bins |psi(R,alpha)|^2 onto a radial grid by mapping
      (R, alpha) -> (r1, r2) = (R cos alpha, R sin alpha). Simple but limited
      by histogram resolution.

  algebraic (new): Uses the exact coordinate transform identity
      n(r) = (2/h_alpha) * integral |F(R)|^2 D(arccos(r/R); R) / sqrt(R^2-r^2) dR
      where D(alpha; R) is the angular density interpolated to the exact alpha
      via cubic spline. Eliminates binning artifacts entirely.

Key relations:
    r1 = R * cos(alpha),  r2 = R * sin(alpha)
    Z_eff(r) = Z - N_core(r)
    N_core(r) = integral_0^r n(r') dr'

where n(r) is the radial number density normalized to 2 (two core electrons).
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
from typing import Optional, Tuple

from geovac.hyperspherical_radial import solve_helium
from geovac.hyperspherical_angular import solve_angular


def compute_core_density(
    Z: float,
    l_max: int = 2,
    n_alpha: int = 200,
    n_radial: int = 600,
    r_max: float = 20.0,
    N_R_angular: int = 200,
    N_R_radial: int = 2000,
    R_max: float = 30.0,
    verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the spherically-averaged one-electron density of a two-electron core.

    Solves the two-electron hyperspherical problem at nuclear charge Z, then
    computes the radial number density n(r) by binning the squared wavefunction
    over hyperradius R and hyperangle alpha.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        Number of FD grid points for hyperangle alpha.
    n_radial : int
        Number of output r grid points.
    r_max : float
        Maximum r for the output grid (bohr).
    N_R_angular : int
        Number of R points for computing adiabatic curves.
    N_R_radial : int
        Number of grid points for the radial solve.
    R_max : float
        Maximum hyperradius (bohr).
    verbose : bool
        Print progress.

    Returns
    -------
    r_grid : ndarray of shape (n_radial,)
        Radial distance grid (bohr), bin centers.
    n_r : ndarray of shape (n_radial,)
        Radial number density n(r) normalized so that integral n(r) dr = 2.
    """
    # Step 1: Solve the two-electron problem
    result = solve_helium(
        Z=Z, l_max=l_max, n_alpha=n_alpha,
        N_R_angular=N_R_angular, N_R_radial=N_R_radial,
        R_max=R_max, verbose=verbose,
    )

    F_R = result['wavefunction']  # shape (N_R_radial,) for single-channel
    R_rad = result['R_grid_radial']
    R_ang = result['R_grid_angular']

    # Step 2: Interpolate F(R) onto the angular R grid
    F_spline = CubicSpline(R_rad, F_R)
    F_on_ang = F_spline(R_ang)

    # Step 3: Set up output r grid and alpha grid
    dr = r_max / n_radial
    r_grid = np.linspace(dr / 2, r_max - dr / 2, n_radial)

    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha
    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    n_r = np.zeros(n_radial)

    # Step 4: Trapezoid weights for the non-uniform R_ang grid
    dR = np.zeros(len(R_ang))
    dR[0] = 0.5 * (R_ang[1] - R_ang[0])
    dR[-1] = 0.5 * (R_ang[-1] - R_ang[-2])
    for i in range(1, len(R_ang) - 1):
        dR[i] = 0.5 * (R_ang[i + 1] - R_ang[i - 1])

    # Step 5: Loop over R grid, compute angular eigenvectors, bin density
    n_l = l_max + 1

    for i_R, R in enumerate(R_ang):
        F_sq = F_on_ang[i_R] ** 2
        if F_sq < 1e-30:
            continue

        # Compute angular eigenvector at this R
        _, vecs = solve_angular(R, Z, l_max, n_alpha, n_channels=1)
        u = vecs[0]  # shape (n_l * n_alpha,)

        # Angular density at each alpha point: sum over partial waves
        ang_density = np.zeros(n_alpha)
        for l_idx in range(n_l):
            ang_density += u[l_idx * n_alpha:(l_idx + 1) * n_alpha] ** 2

        # Weight = |F(R)|^2 * dR * angular_density
        weights = F_sq * dR[i_R] * ang_density

        # Map to individual electron distances
        r1 = R * cos_a
        r2 = R * sin_a

        # Bin contributions from both electrons (vectorized)
        bins1 = (r1 / dr).astype(int)
        bins2 = (r2 / dr).astype(int)

        mask1 = (bins1 >= 0) & (bins1 < n_radial)
        mask2 = (bins2 >= 0) & (bins2 < n_radial)

        np.add.at(n_r, bins1[mask1], weights[mask1])
        np.add.at(n_r, bins2[mask2], weights[mask2])

    # Convert from accumulated probability to density: n(r) = accumulated / dr
    n_r_density = n_r / dr

    # Renormalize to exactly 2 (compensate for discretization errors)
    total = np.trapezoid(n_r_density, r_grid)
    if total > 0:
        n_r_density *= 2.0 / total

    return r_grid, n_r_density


def compute_core_density_algebraic(
    Z: float,
    l_max: int = 2,
    n_alpha: int = 200,
    n_radial: int = 600,
    r_max: float = 20.0,
    N_R_angular: int = 200,
    N_R_radial: int = 2000,
    R_max: float = 30.0,
    verbose: bool = False,
    helium_result: Optional[dict] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the one-electron density algebraically from channel coefficients.

    Instead of histogram binning, evaluates the exact coordinate-transform
    identity:

        n(r) = (2/h_alpha) * sum_R |F(R)|^2 dR * D(arccos(r/R); R) / sqrt(R^2 - r^2)

    where D(alpha; R) = sum_l |u_l(alpha)|^2 is the angular density in the
    Liouville-substituted FD basis, interpolated to the exact hyperangle
    alpha = arccos(r/R) via cubic spline. The factor 2 uses singlet exchange
    symmetry. The Jacobian 1/sqrt(R^2 - r^2) comes from the delta function
    delta(R cos alpha - r).

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        Number of FD grid points for hyperangle alpha.
    n_radial : int
        Number of output r grid points.
    r_max : float
        Maximum r for the output grid (bohr).
    N_R_angular : int
        Number of R points for computing adiabatic curves.
    N_R_radial : int
        Number of grid points for the radial solve.
    R_max : float
        Maximum hyperradius (bohr).
    verbose : bool
        Print progress.
    helium_result : dict, optional
        Pre-solved helium result from solve_helium(). If None, solves
        internally (avoids double-solve when called from CoreScreening).

    Returns
    -------
    r_grid : ndarray of shape (n_radial,)
        Radial distance grid (bohr).
    n_r : ndarray of shape (n_radial,)
        Radial number density n(r) normalized so that integral n(r) dr = 2.
    """
    # Step 1: Solve the two-electron problem (or use pre-solved)
    if helium_result is None:
        result = solve_helium(
            Z=Z, l_max=l_max, n_alpha=n_alpha,
            N_R_angular=N_R_angular, N_R_radial=N_R_radial,
            R_max=R_max, verbose=verbose,
        )
    else:
        result = helium_result

    F_R = result['wavefunction']
    R_rad = result['R_grid_radial']
    R_ang = result['R_grid_angular']

    # Step 2: Interpolate F(R) onto the angular R grid
    F_spline = CubicSpline(R_rad, F_R)
    F_on_ang = F_spline(R_ang)

    # Step 3: Set up output r grid
    dr = r_max / n_radial
    r_grid = np.linspace(dr / 2, r_max - dr / 2, n_radial)

    # Alpha FD grid (must match solve_angular)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    n_l = l_max + 1

    # Step 4: Trapezoid weights for R_ang
    dR = np.zeros(len(R_ang))
    dR[0] = 0.5 * (R_ang[1] - R_ang[0])
    dR[-1] = 0.5 * (R_ang[-1] - R_ang[-2])
    for i in range(1, len(R_ang) - 1):
        dR[i] = 0.5 * (R_ang[i + 1] - R_ang[i - 1])

    # Step 5: Build angular density splines at each R
    # D(alpha; R) = sum_l |u_l(alpha)|^2 interpolated as cubic spline
    # Boundary conditions: u(0) = u(pi/2) = 0 (Dirichlet), so D(0) = D(pi/2) = 0
    ang_splines = []
    for i_R, R in enumerate(R_ang):
        F_sq = F_on_ang[i_R] ** 2
        if F_sq < 1e-30:
            ang_splines.append(None)
            continue

        _, vecs = solve_angular(R, Z, l_max, n_alpha, n_channels=1)
        u = vecs[0]

        ang_density = np.zeros(n_alpha)
        for l_idx in range(n_l):
            ang_density += u[l_idx * n_alpha:(l_idx + 1) * n_alpha] ** 2

        # Extend with boundary zeros for clean spline interpolation
        alpha_ext = np.concatenate([[0.0], alpha_grid, [np.pi / 2]])
        dens_ext = np.concatenate([[0.0], ang_density, [0.0]])
        ang_splines.append(CubicSpline(alpha_ext, dens_ext))

    # Step 6: Evaluate algebraic density at each r
    # n(r) = (2/h_alpha) * sum_R |F(R)|^2 dR * D(arccos(r/R); R) / sqrt(R^2 - r^2)
    n_r = np.zeros(n_radial)

    for i_R, R in enumerate(R_ang):
        F_sq = F_on_ang[i_R] ** 2
        if F_sq < 1e-30 or ang_splines[i_R] is None:
            continue

        # r values that can be reached from this R: r < R
        # The singularity at r=R is integrable (D ~ alpha^2 -> 0 while
        # Jacobian ~ 1/alpha -> inf, product ~ alpha -> 0). Use absolute
        # tolerance based on alpha grid spacing to avoid numerical issues.
        mask = r_grid < R * np.cos(h_alpha * 0.5)
        if not np.any(mask):
            continue

        r_valid = r_grid[mask]

        # Exact hyperangle for each r at this R
        alpha_vals = np.arccos(r_valid / R)

        # Evaluate angular density at exact alpha values
        ang_vals = ang_splines[i_R](alpha_vals)
        ang_vals = np.maximum(ang_vals, 0.0)

        # Jacobian from delta(R cos alpha - r)
        jacobian = 1.0 / np.sqrt(R**2 - r_valid**2)

        # Accumulate
        n_r[mask] += F_sq * dR[i_R] * ang_vals * jacobian

    # Apply symmetry factor and FD normalization
    n_r *= 2.0 / h_alpha

    # Renormalize to exactly 2
    total = np.trapezoid(n_r, r_grid)
    if verbose:
        print(f"Algebraic density raw integral: {total:.6f} (target: 2.0)")
    if total > 0:
        n_r *= 2.0 / total

    return r_grid, n_r


def compute_z_eff(
    Z: float,
    r_grid: np.ndarray,
    l_max: int = 2,
    n_alpha: int = 200,
    n_radial: int = 600,
    r_max: float = 20.0,
    N_R_angular: int = 200,
    N_R_radial: int = 2000,
    R_max: float = 30.0,
    verbose: bool = False,
) -> np.ndarray:
    """
    Compute Z_eff(r) = Z - N_core(r) on a user-supplied r grid.

    N_core(r) is the cumulative number of core electrons enclosed within
    radius r, obtained by integrating the one-electron density.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    r_grid : ndarray
        Radial distances at which to evaluate Z_eff (bohr).
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        Number of FD grid points for hyperangle alpha.
    n_radial : int
        Number of internal r grid points for density computation.
    r_max : float
        Maximum r for the internal density grid (bohr).
    N_R_angular : int
        Number of R points for adiabatic curves.
    N_R_radial : int
        Number of grid points for the radial solve.
    R_max : float
        Maximum hyperradius (bohr).
    verbose : bool
        Print progress.

    Returns
    -------
    z_eff : ndarray, same shape as r_grid
        Effective nuclear charge at each r.
    """
    r_internal, n_r = compute_core_density(
        Z, l_max, n_alpha, n_radial, r_max,
        N_R_angular, N_R_radial, R_max, verbose,
    )

    # Cumulative integral: N_core(r) = integral_0^r n(r') dr'
    N_core = np.zeros_like(r_internal)
    N_core[1:] = cumulative_trapezoid(n_r, r_internal)

    # Interpolate N_core to user grid, clamping to [0, 2]
    N_core_spline = CubicSpline(r_internal, N_core, extrapolate=True)
    r_arr = np.asarray(r_grid, dtype=float)
    N_at_r = N_core_spline(r_arr)
    N_at_r = np.clip(N_at_r, 0.0, 2.0)

    # Handle boundaries: r <= 0 -> Z, r very large -> Z - 2
    z_eff = Z - N_at_r
    z_eff = np.where(r_arr <= 0, Z, z_eff)

    return z_eff


class CoreScreening:
    """
    Reusable core screening object for a two-electron system.

    Usage
    -----
    >>> cs = CoreScreening(Z=3)
    >>> cs.solve()
    >>> z_eff_at_1bohr = cs.z_eff(1.0)
    >>> z_eff_array = cs.z_eff(np.array([0.1, 0.5, 1.0, 2.0]))

    Parameters
    ----------
    method : str
        Density extraction method: 'algebraic' (exact coordinate transform
        with spline interpolation) or 'histogram' (original binning method).
    """

    def __init__(
        self,
        Z: float,
        l_max: int = 2,
        n_alpha: int = 200,
        n_radial: int = 600,
        r_max: float = 20.0,
        N_R_angular: int = 200,
        N_R_radial: int = 2000,
        R_max: float = 30.0,
        method: str = 'algebraic',
    ) -> None:
        if method not in ('algebraic', 'histogram'):
            raise ValueError(f"method must be 'algebraic' or 'histogram', got '{method}'")
        self.Z = Z
        self.l_max = l_max
        self.n_alpha = n_alpha
        self.n_radial = n_radial
        self.r_max = r_max
        self.N_R_angular = N_R_angular
        self.N_R_radial = N_R_radial
        self.R_max = R_max
        self.method = method

        self._solved = False
        self._r_grid: Optional[np.ndarray] = None
        self._n_r: Optional[np.ndarray] = None
        self._N_core: Optional[np.ndarray] = None
        self._N_core_spline: Optional[CubicSpline] = None
        self._n_r_spline: Optional[CubicSpline] = None
        self._helium_result: Optional[dict] = None

    def solve(self, verbose: bool = False) -> None:
        """Solve the two-electron problem and compute the screening function."""
        self._helium_result = solve_helium(
            Z=self.Z, l_max=self.l_max, n_alpha=self.n_alpha,
            N_R_angular=self.N_R_angular, N_R_radial=self.N_R_radial,
            R_max=self.R_max, verbose=verbose,
        )

        if self.method == 'algebraic':
            self._r_grid, self._n_r = compute_core_density_algebraic(
                self.Z, self.l_max, self.n_alpha, self.n_radial,
                self.r_max, self.N_R_angular, self.N_R_radial,
                self.R_max, verbose=verbose,
                helium_result=self._helium_result,
            )
        else:
            self._r_grid, self._n_r = compute_core_density(
                self.Z, self.l_max, self.n_alpha, self.n_radial,
                self.r_max, self.N_R_angular, self.N_R_radial,
                self.R_max, verbose=verbose,
            )

        # Cumulative integral for N_core
        self._N_core = np.zeros_like(self._r_grid)
        self._N_core[1:] = cumulative_trapezoid(self._n_r, self._r_grid)

        # Build interpolation splines
        self._N_core_spline = CubicSpline(
            self._r_grid, self._N_core, extrapolate=True
        )
        self._n_r_spline = CubicSpline(
            self._r_grid, self._n_r, extrapolate=True
        )

        self._solved = True

    def _check_solved(self) -> None:
        if not self._solved:
            raise RuntimeError("Call .solve() before querying screening data.")

    def z_eff(self, r: float | np.ndarray) -> float | np.ndarray:
        """
        Effective nuclear charge at distance r from the nucleus.

        Parameters
        ----------
        r : float or ndarray
            Radial distance(s) in bohr.

        Returns
        -------
        float or ndarray
            Z_eff(r) = Z - N_core(r).
        """
        self._check_solved()
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        N_at_r = np.clip(self._N_core_spline(r_arr), 0.0, 2.0)
        result = self.Z - N_at_r
        result = np.where(r_arr <= 0, self.Z, result)

        if np.ndim(r) == 0:
            return float(result[0])
        return result

    def density(self, r: float | np.ndarray) -> float | np.ndarray:
        """
        One-electron radial number density n(r) at distance r.

        Normalized so that integral n(r) dr = 2.

        Parameters
        ----------
        r : float or ndarray
            Radial distance(s) in bohr.

        Returns
        -------
        float or ndarray
            n(r) at the requested points.
        """
        self._check_solved()
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        result = np.clip(self._n_r_spline(r_arr), 0.0, None)

        if np.ndim(r) == 0:
            return float(result[0])
        return result

    @property
    def r_grid(self) -> np.ndarray:
        """The internal r grid used for the density computation."""
        self._check_solved()
        return self._r_grid

    @property
    def n_r(self) -> np.ndarray:
        """The radial number density on the internal grid."""
        self._check_solved()
        return self._n_r

    @property
    def N_core(self) -> np.ndarray:
        """The cumulative core electron count on the internal grid."""
        self._check_solved()
        return self._N_core

    @property
    def energy(self) -> float:
        """Ground-state energy from the hyperspherical solve."""
        self._check_solved()
        return self._helium_result['energy']
