"""
Hyperradial solver for the adiabatic hyperspherical method.

Solves the 1D Schrodinger equation in the hyperradius R:
    [-1/2 d^2F/dR^2 + V_eff(R)] F(R) = E F(R)

where V_eff(R) = mu(R)/R^2 + 15/(8R^2), with mu(R) from the angular solve.

Uses self-adjoint finite differences (same strategy as Paper 11's
prolate spheroidal radial solver).

References:
  - Macek, J. Phys. B 1, 831 (1968)
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import CubicSpline
from typing import Tuple

from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
)


def solve_radial(
    V_eff_func,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R: int = 2000,
    n_states: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the hyperradial Schrodinger equation.

    [-1/2 d^2F/dR^2 + V_eff(R) F(R)] = E F(R)

    with Dirichlet boundary conditions F(R_min) = F(R_max) = 0.

    Parameters
    ----------
    V_eff_func : callable or ndarray
        Effective potential V_eff(R). If callable, evaluated on grid.
    R_min : float
        Left boundary (bohr). Must be > 0.
    R_max : float
        Right boundary (bohr).
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenstates to return.

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha).
    F : ndarray of shape (n_states, N_R)
        Radial wavefunctions on the grid.
    R_grid : ndarray of shape (N_R,)
        Grid points.
    """
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    if callable(V_eff_func):
        V = V_eff_func(R_grid)
    else:
        V = np.asarray(V_eff_func)

    # FD for -1/2 d^2/dR^2 + V(R)
    diag = np.ones(N_R) / h**2 + V
    off_diag = -0.5 * np.ones(N_R - 1) / h**2

    evals, evecs = eigh_tridiagonal(
        diag, off_diag,
        select='i', select_range=(0, n_states - 1)
    )

    for i in range(n_states):
        norm = np.sqrt(h * np.sum(evecs[:, i]**2))
        if norm > 0:
            evecs[:, i] /= norm

    return evals, evecs.T, R_grid


def solve_helium(
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    N_R_angular: int = 200,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_channels: int = 1,
    verbose: bool = True,
) -> dict:
    """
    Full single-channel adiabatic hyperspherical solver for He.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        FD grid points for alpha.
    N_R_angular : int
        Number of R points for computing mu(R).
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Number of grid points for the radial solve.
    n_channels : int
        Number of adiabatic channels (1 for single-channel).
    verbose : bool
        Print progress information.

    Returns
    -------
    result : dict
        Keys: 'energy', 'R_grid_angular', 'mu_curves', 'V_eff',
        'R_grid_radial', 'wavefunction', 'Z', 'l_max', 'n_channels'.
    """
    import time

    t0 = time.time()

    # --- Step 1: Compute adiabatic curves ---
    # Denser grid near small R where mu(R) varies rapidly
    R_grid_ang = np.concatenate([
        np.linspace(0.1, 1.0, N_R_angular // 3),
        np.linspace(1.0, 5.0, N_R_angular // 3),
        np.linspace(5.0, R_max, N_R_angular // 3 + 1),
    ])
    R_grid_ang = np.unique(R_grid_ang)

    if verbose:
        print(f"Computing adiabatic curves on {len(R_grid_ang)} R points...")
        print(f"  Z={Z}, l_max={l_max}, n_alpha={n_alpha}, "
              f"n_channels={n_channels}")

    mu = compute_adiabatic_curve(
        R_grid_ang, Z, l_max, n_alpha, n_channels
    )

    t1 = time.time()
    if verbose:
        print(f"  Angular solve: {t1 - t0:.2f}s")
        V_eff_last = effective_potential(
            R_grid_ang[-1:], mu[0, -1:]
        )[0]
        print(f"  V_eff(R={R_grid_ang[-1]:.1f}) = {V_eff_last:.4f} Ha "
              f"(should approach {-Z**2/2:.1f})")

    # --- Step 2: Interpolate V_eff for radial solver ---
    V_eff_angular = effective_potential(R_grid_ang, mu[0])
    V_eff_spline = CubicSpline(R_grid_ang, V_eff_angular, extrapolate=True)

    # --- Step 3: Solve hyperradial equation ---
    if verbose:
        print(f"Solving hyperradial equation (N_R={N_R_radial})...")

    E, F, R_grid_rad = solve_radial(
        V_eff_spline, R_min, R_max, N_R_radial, n_states=1
    )

    t2 = time.time()
    E_exact = -2.903724
    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")
        print(f"\n  Ground state energy: {E[0]:.6f} Ha")
        print(f"  Exact He:            {E_exact:.6f} Ha")
        err = abs(E[0] - E_exact) / abs(E_exact) * 100
        print(f"  Error:               {err:.2f}%")
        print(f"  Total time:          {t2 - t0:.2f}s")

    return {
        'energy': E[0],
        'R_grid_angular': R_grid_ang,
        'mu_curves': mu,
        'V_eff': V_eff_angular,
        'R_grid_radial': R_grid_rad,
        'wavefunction': F[0],
        'Z': Z,
        'l_max': l_max,
        'n_channels': n_channels,
    }
