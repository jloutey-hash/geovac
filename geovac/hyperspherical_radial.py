"""
Hyperradial solver for the adiabatic hyperspherical method.

Solves the 1D Schrodinger equation in the hyperradius R:
    [-1/2 d^2F/dR^2 + V_eff(R)] F(R) = E F(R)

where V_eff(R) = mu(R)/R^2 + 15/(8R^2), with mu(R) from the angular solve.

For coupled channels, solves the matrix equation:
    [-1/2 d^2/dR^2 δ_μν + V_eff_μ δ_μν - P_μν d/dR - 1/2 Q_μν] F_ν = E F_μ

Uses self-adjoint finite differences (same strategy as Paper 11's
prolate spheroidal radial solver).

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import CubicSpline
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from typing import Tuple, Optional

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


def solve_coupled_radial(
    V_eff_splines: list,
    P_splines: np.ndarray,
    Q_splines: Optional[np.ndarray] = None,
    n_channels: int = 2,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R: int = 3000,
    n_states: int = 5,
    sigma: float = -3.0,
    include_Q: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the coupled-channel hyperradial equation.

    [-1/2 d²/dR² δ_μν + V_eff_μ(R) δ_μν - P_μν(R) d/dR (- 1/2 Q_μν)] F_ν = E F_μ

    Builds a block-tridiagonal sparse matrix and uses shift-invert eigsh.

    The Q term (second-derivative coupling) is optional. When n_channels is
    small, the closure approximation Q ≈ P² can be unreliable and is omitted
    by default. The P-only formulation captures the dominant non-adiabatic
    coupling.

    Parameters
    ----------
    V_eff_splines : list of callables
        V_eff_μ(R) for each channel μ, length n_channels.
    P_splines : ndarray of callables, shape (n_channels, n_channels)
        P_μν(R) interpolating functions.
    Q_splines : ndarray of callables, optional
        Q_μν(R) interpolating functions. Only used if include_Q=True.
    n_channels : int
        Number of coupled channels.
    R_min, R_max : float
        Radial boundaries.
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenvalues to find.
    sigma : float
        Shift for shift-invert mode (target energy region).
    include_Q : bool
        Whether to include the Q (second-derivative) coupling term.

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha), sorted ascending.
    F : ndarray of shape (n_states, n_channels, N_R)
        Channel wavefunctions on the grid.
    R_grid : ndarray of shape (N_R,)
        Grid points.
    """
    n_ch = n_channels
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    # Evaluate potentials and couplings on grid
    V_eff = np.zeros((n_ch, N_R))
    for mu in range(n_ch):
        V_eff[mu] = V_eff_splines[mu](R_grid)

    P_grid = np.zeros((n_ch, n_ch, N_R))
    for mu in range(n_ch):
        for nu in range(n_ch):
            P_grid[mu, nu] = P_splines[mu][nu](R_grid)

    Q_grid = np.zeros((n_ch, n_ch, N_R))
    if include_Q and Q_splines is not None:
        for mu in range(n_ch):
            for nu in range(n_ch):
                Q_grid[mu, nu] = Q_splines[mu][nu](R_grid)

    # Build block-tridiagonal sparse matrix of dimension (N_R * n_ch) x (N_R * n_ch)
    # Indexing: global index = i * n_ch + mu  (grid point i, channel mu)
    dim = N_R * n_ch
    H = lil_matrix((dim, dim))

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    for i in range(N_R):
        for mu in range(n_ch):
            row = i * n_ch + mu

            # Diagonal block: kinetic + V_eff (+ Q if included)
            H[row, row] = kinetic_diag + V_eff[mu, i] - 0.5 * Q_grid[mu, mu, i]

            # Off-diagonal channels at same grid point (from Q)
            if include_Q:
                for nu in range(n_ch):
                    if nu != mu:
                        col = i * n_ch + nu
                        H[row, col] += -0.5 * Q_grid[mu, nu, i]

        # Off-diagonal in R: coupling to i+1
        if i < N_R - 1:
            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i + 1) * n_ch + nu
                    # Kinetic: delta_mu_nu * kinetic_off
                    # P coupling: -P_mu_nu * 1/(2h)
                    val = 0.0
                    if mu == nu:
                        val += kinetic_off
                    val -= P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

        # Off-diagonal in R: coupling to i-1
        if i > 0:
            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i - 1) * n_ch + nu
                    # Kinetic: delta_mu_nu * kinetic_off
                    # P coupling: +P_mu_nu * 1/(2h)
                    val = 0.0
                    if mu == nu:
                        val += kinetic_off
                    val += P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

    H_csr = csr_matrix(H)

    # Shift-invert near sigma to find states near ground state
    evals, evecs = eigsh(H_csr, k=n_states, sigma=sigma, which='LM')

    # Sort by energy
    idx = np.argsort(evals)
    evals = evals[idx]
    evecs = evecs[:, idx]

    # Reshape eigenvectors: (N_R * n_ch, n_states) -> (n_states, n_ch, N_R)
    F = np.zeros((n_states, n_ch, N_R))
    for s in range(n_states):
        for mu in range(n_ch):
            F[s, mu, :] = evecs[mu::n_ch, s]

        # Normalize each state
        norm_sq = h * np.sum(F[s]**2)
        if norm_sq > 0:
            F[s] /= np.sqrt(norm_sq)

    return evals, F, R_grid


def solve_helium(
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    N_R_angular: int = 200,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_channels: int = 1,
    coupled: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Full hyperspherical solver for He.

    Single-channel adiabatic (coupled=False) or coupled-channel (coupled=True).

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
        Number of adiabatic channels.
    coupled : bool
        If True, include non-adiabatic coupling between channels.
    verbose : bool
        Print progress information.

    Returns
    -------
    result : dict
        Keys: 'energy', 'R_grid_angular', 'mu_curves', 'V_eff',
        'R_grid_radial', 'wavefunction', 'Z', 'l_max', 'n_channels',
        'coupled'.
        If coupled: also 'P_coupling', 'channel_weights'.
    """
    import time

    t0 = time.time()

    # --- Step 1: Compute adiabatic curves ---
    R_grid_ang = np.concatenate([
        np.linspace(0.1, 1.0, N_R_angular // 3),
        np.linspace(1.0, 5.0, N_R_angular // 3),
        np.linspace(5.0, R_max, N_R_angular // 3 + 1),
    ])
    R_grid_ang = np.unique(R_grid_ang)

    if verbose:
        mode = "coupled" if coupled else "single-channel"
        print(f"Computing adiabatic curves ({mode}, {n_channels} ch) "
              f"on {len(R_grid_ang)} R points...")
        print(f"  Z={Z}, l_max={l_max}, n_alpha={n_alpha}")

    if coupled and n_channels > 1:
        # Use coupling module with Hellmann-Feynman P and DBOC
        from geovac.hyperspherical_coupling import compute_coupling_matrices

        coupling = compute_coupling_matrices(
            R_grid_ang, Z, l_max, n_alpha, n_channels
        )
        mu_curves = coupling['mu']
        P = coupling['P']
        DBOC = coupling['DBOC']

        t1 = time.time()
        if verbose:
            print(f"  Angular + coupling solve: {t1 - t0:.2f}s")
            print(f"  DBOC peak (ch 0): {np.max(DBOC[0]):.6f} Ha")

        # --- Step 2: Build V_eff splines ---
        # V_eff_μ(R) = μ(R)/R² + 15/(8R²)
        # Note: DBOC is NOT added to diagonal. The off-diagonal P coupling
        # implicitly includes the DBOC through the coupled-channel equations.
        # Adding DBOC explicitly would double-count because the DBOC and
        # off-diagonal coupling nearly cancel (~97% cancellation for He).
        V_eff_splines = []
        for ch in range(n_channels):
            V_eff_ch = effective_potential(R_grid_ang, mu_curves[ch])
            V_eff_splines.append(
                CubicSpline(R_grid_ang, V_eff_ch, extrapolate=True)
            )

        P_splines = [[None] * n_channels for _ in range(n_channels)]
        for mu_idx in range(n_channels):
            for nu_idx in range(n_channels):
                P_splines[mu_idx][nu_idx] = CubicSpline(
                    R_grid_ang, P[mu_idx, nu_idx], extrapolate=True
                )

        # --- Step 3: Solve coupled radial equation ---
        if verbose:
            print(f"Solving coupled-channel radial equation "
                  f"(N_R={N_R_radial}, {n_channels} ch)...")

        E_all, F_all, R_grid_rad = solve_coupled_radial(
            V_eff_splines, P_splines,
            n_channels=n_channels,
            R_min=R_min, R_max=R_max, N_R=N_R_radial,
            n_states=min(5, n_channels * 3),
            sigma=-3.0,
        )

        t2 = time.time()
        E_exact = -2.903724

        # Channel weights for ground state
        h_rad = R_grid_rad[1] - R_grid_rad[0]
        weights = np.array([
            h_rad * np.sum(F_all[0, ch]**2)
            for ch in range(n_channels)
        ])
        weights /= weights.sum()

        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")
            print(f"\n  Ground state energy: {E_all[0]:.6f} Ha")
            print(f"  Exact He:            {E_exact:.6f} Ha")
            err = abs(E_all[0] - E_exact) / abs(E_exact) * 100
            print(f"  Error:               {err:.2f}%")
            print(f"  Channel weights:     {weights}")
            print(f"  Total time:          {t2 - t0:.2f}s")

        V_eff_angular = effective_potential(R_grid_ang, mu_curves[0])

        return {
            'energy': E_all[0],
            'energies': E_all,
            'R_grid_angular': R_grid_ang,
            'mu_curves': mu_curves,
            'V_eff': V_eff_angular,
            'DBOC': DBOC,
            'R_grid_radial': R_grid_rad,
            'wavefunction': F_all[0],
            'channel_weights': weights,
            'P_coupling': P,
            'Z': Z,
            'l_max': l_max,
            'n_channels': n_channels,
            'coupled': True,
        }

    else:
        # Single-channel adiabatic (original code path)
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
        V_eff_spline = CubicSpline(
            R_grid_ang, V_eff_angular, extrapolate=True
        )

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
            'coupled': False,
        }
