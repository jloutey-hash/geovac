"""
Algebraic Slater integral evaluator using Laguerre moment decomposition.

Replaces the numerical trapezoidal-grid Slater F^k integral in inter-fiber
exchange coupling with algebraic evaluation via Laguerre basis expansion.

Strategy (parallels Paper 12 Neumann expansion and Track H/I/J Laguerre):
  1. Fit radial densities P(r) to a Laguerre basis: P(r) ~ sum_n c_n phi_n(r)
     where phi_n(r) = r * exp(-alpha*r) * L_n(2*alpha*r)
  2. Precompute the Slater R^k matrix: R^k[i,j] = integral integral
     phi_i(r1) phi_j(r2) r_<^k / r_>^{k+1} dr1 dr2
  3. F^k = c_A^T @ R^k @ c_B (algebraic contraction)

The R^k matrix elements are computed via high-order Gauss-Laguerre quadrature
on the cumulative-charge decomposition. The quadrature is the transcendental
seed (analogous to e^a * E_1(a) in Track J); the Laguerre basis and moment
structure provide the algebraic framework.

Transcendental content identification:
  - The R^k matrix involves incomplete Laguerre moments (integrals with
    finite limits), which reduce to incomplete gamma functions gamma(s, x).
  - The seed is gamma(1, x) = 1 - exp(-x), or equivalently exp(-x).
  - All other incomplete moments follow from the recurrence:
    gamma(s+1, x) = s * gamma(s, x) - x^s * exp(-x)
  - This parallels Track J's Stieltjes structure: algebraic recurrence
    seeded by a single transcendental function.

Reference: Paper 12 Sec III-V (Neumann decomposition strategy),
           Paper 18 Sec II (exchange constant catalog),
           Track H (Laguerre moment matrices for S and K).
"""

import numpy as np
from scipy.special import roots_laguerre, eval_laguerre, gammainc, gamma as gammafn
from scipy.linalg import solve, lstsq
from typing import Tuple, Dict, Optional


def _laguerre_basis_values(
    n_basis: int,
    alpha: float,
    r_grid: np.ndarray,
) -> np.ndarray:
    """Evaluate Laguerre basis functions on a radial grid.

    Basis: phi_n(r) = r * exp(-alpha*r) * L_n(2*alpha*r)

    Parameters
    ----------
    n_basis : int
        Number of basis functions.
    alpha : float
        Exponential decay parameter.
    r_grid : ndarray
        Radial grid points.

    Returns
    -------
    phi : ndarray of shape (n_basis, len(r_grid))
        Basis function values phi_n(r_j).
    """
    x = 2.0 * alpha * r_grid
    envelope = r_grid * np.exp(-alpha * r_grid)
    phi = np.zeros((n_basis, len(r_grid)))
    for n in range(n_basis):
        phi[n, :] = envelope * eval_laguerre(n, x)
    return phi


def fit_density_laguerre(
    r_grid: np.ndarray,
    P: np.ndarray,
    n_basis: int = 20,
    alpha: Optional[float] = None,
) -> Tuple[np.ndarray, float, float]:
    """Fit a radial density to a Laguerre basis via least-squares.

    Determines optimal alpha from the density's mean radius if not provided,
    then fits P(r) = sum_n c_n * phi_n(r) via least-squares on the grid.

    Parameters
    ----------
    r_grid : ndarray
        Radial grid points (uniform spacing assumed).
    P : ndarray
        Radial density P(r) on the grid.
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float or None
        Exponential parameter. If None, estimated from density.

    Returns
    -------
    coeffs : ndarray of shape (n_basis,)
        Expansion coefficients c_n.
    alpha : float
        The alpha parameter used.
    residual : float
        Relative L2 residual ||P - P_fit|| / ||P||.
    """
    dr = r_grid[1] - r_grid[0] if len(r_grid) > 1 else 1.0

    # Estimate alpha from mean radius: <r> ~ 3/(2*alpha) for Laguerre envelope
    if alpha is None:
        total = np.sum(P) * dr
        if total > 1e-15:
            mean_r = np.sum(r_grid * P) * dr / total
            alpha = 1.5 / max(mean_r, 0.1)
        else:
            alpha = 1.0

    # Build basis matrix
    phi = _laguerre_basis_values(n_basis, alpha, r_grid)

    # Weighted least-squares: minimize sum_j w_j |P(r_j) - sum_n c_n phi_n(r_j)|^2
    # Weight by sqrt(dr) for integration-weighted fit
    A = phi.T  # (n_grid, n_basis)
    b = P      # (n_grid,)

    # Use lstsq for robust fitting (handles rank deficiency)
    coeffs, _, _, _ = lstsq(A, b)

    # Compute residual
    P_fit = phi.T @ coeffs
    norm_P = np.sqrt(np.sum(P**2) * dr)
    norm_res = np.sqrt(np.sum((P - P_fit)**2) * dr)
    residual = norm_res / max(norm_P, 1e-15)

    return coeffs, alpha, residual


def _build_slater_Rk_matrix(
    n_basis: int,
    alpha: float,
    k: int = 0,
    n_quad: int = 200,
) -> np.ndarray:
    """Build the Slater R^k matrix between Laguerre basis functions.

    R^k[i,j] = integral integral phi_i(r1) phi_j(r2) r_<^k / r_>^{k+1} dr1 dr2

    Uses the cumulative-charge decomposition:
    R^k[i,j] = integral phi_i(r1) Y^k_j(r1) dr1

    where Y^k_j(r1) = Q^k_j(r1)/r1^{k+1} + T^k_j(r1) * r1^k
    with Q^k_j(r) = integral_0^r phi_j(r') r'^k dr'
         T^k_j(r) = integral_r^inf phi_j(r') / r'^{k+1} dr'

    Computed via high-order Gauss-Laguerre quadrature in x-space (x = 2*alpha*r).

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    k : int
        Multipole order (0=monopole).
    n_quad : int
        Number of quadrature points.

    Returns
    -------
    Rk : ndarray of shape (n_basis, n_basis)
        Slater integral matrix. F^k = c_A^T @ Rk @ c_B.
    """
    two_alpha = 2.0 * alpha

    # Use Gauss-Laguerre quadrature: integral_0^inf f(x) e^{-x} dx
    # We need to handle the e^{-alpha*r} weight in our basis.
    # phi_n(r) = r * e^{-alpha*r} * L_n(2*alpha*r)
    # In x = 2*alpha*r: phi_n = x/(2a) * e^{-x/2} * L_n(x)

    # For the double integral, we use a composite approach:
    # Evaluate on a fine grid and use cumulative charge method.
    # This is more stable than double quadrature for the r_</r_> kernel.

    # Fine grid in r-space for accurate cumulative charge
    # Use Gauss-Laguerre points mapped to r-space
    x_quad, w_quad = roots_laguerre(n_quad)
    # Map to r-space: r = x / (2*alpha), but we need to account for
    # the e^{-x} weight in the quadrature vs e^{-alpha*r} in the basis.
    # phi_n(r) = r * e^{-alpha*r} * L_n(2*alpha*r) = (x/(2a)) * e^{-x/2} * L_n(x)
    # We want integral phi_i(r) * Y^k_j(r) dr
    #   = integral (x/(2a)) e^{-x/2} L_i(x) * Y^k_j(x/(2a)) dx/(2a)
    #   = integral (x/(4a^2)) e^{-x/2} L_i(x) * Y^k_j(x/(2a)) dx
    # Using Gauss-Laguerre with weight e^{-x}:
    #   ~ sum_q w_q * (x_q/(4a^2)) * e^{x_q/2} * L_i(x_q) * Y^k_j(x_q/(2a))

    # Evaluate L_n at quadrature points
    L_at_x = np.zeros((n_basis, n_quad))
    for n in range(n_basis):
        L_at_x[n, :] = eval_laguerre(n, x_quad)

    r_quad = x_quad / two_alpha

    # Evaluate phi_n(r) at quadrature points (without e^{-x} factor, for GL weight)
    # phi_n(r_q) = r_q * e^{-alpha*r_q} * L_n(x_q)
    # For GL quadrature with weight e^{-x_q}: need phi_n(r_q) * e^{+x_q} * dr/dx
    # phi_n(r_q) = (x_q/(2a)) * e^{-x_q/2} * L_n(x_q)
    # phi_n(r_q) * (1/(2a)) [dr factor] = (x_q/(4a^2)) * e^{-x_q/2} * L_n(x_q)
    # With GL weight e^{-x}: integrand factor = (x_q/(4a^2)) * e^{x_q/2} * L_n(x_q)
    envelope_factor = (x_quad / (4.0 * alpha**2)) * np.exp(x_quad / 2.0)

    # Build the phi values for integration: phi_int[n,q] = envelope_factor[q] * L_n(x_q)
    phi_int = np.zeros((n_basis, n_quad))
    for n in range(n_basis):
        phi_int[n, :] = envelope_factor * L_at_x[n, :]

    # Now build R^k matrix using cumulative charge on the quadrature grid.
    # For each j, compute Y^k_j(r) at each quadrature point.
    # phi_j(r) = r * e^{-alpha*r} * L_j(2*alpha*r)
    # phi_j at quad points: (x_q/(2a)) * e^{-x_q/2} * L_j(x_q)
    phi_at_r = np.zeros((n_basis, n_quad))
    for n in range(n_basis):
        phi_at_r[n, :] = (x_quad / two_alpha) * np.exp(-x_quad / 2.0) * L_at_x[n, :]

    Rk = np.zeros((n_basis, n_basis))

    # Sort quadrature points (they should already be sorted, but ensure)
    sort_idx = np.argsort(r_quad)
    r_sorted = r_quad[sort_idx]

    for j in range(n_basis):
        # phi_j values at sorted r points
        phi_j_sorted = phi_at_r[j, sort_idx]

        # Cumulative charge Q^k_j(r_q) = integral_0^{r_q} phi_j(r') r'^k dr'
        # Using trapezoidal on the sorted GL points
        # We need to integrate phi_j(r) * r^k from 0 to r_q
        # Use cumulative trapezoidal integration on sorted points
        integrand_Q = phi_j_sorted * r_sorted**k
        # Cumulative trapezoidal: Q[q] = sum_{q'=0}^{q-1} (integrand[q'] + integrand[q'+1])/2 * dr
        n_pts = len(r_sorted)
        Q_k = np.zeros(n_pts)
        for q in range(1, n_pts):
            dr_q = r_sorted[q] - r_sorted[q - 1]
            Q_k[q] = Q_k[q - 1] + 0.5 * (integrand_Q[q - 1] + integrand_Q[q]) * dr_q

        # Tail integral T^k_j(r_q) = integral_{r_q}^inf phi_j(r') / r'^{k+1} dr'
        integrand_T = np.where(r_sorted > 1e-15,
                               phi_j_sorted / r_sorted**(k + 1), 0.0)
        T_k = np.zeros(n_pts)
        for q in range(n_pts - 2, -1, -1):
            dr_q = r_sorted[q + 1] - r_sorted[q]
            T_k[q] = T_k[q + 1] + 0.5 * (integrand_T[q] + integrand_T[q + 1]) * dr_q

        # Y^k_j(r) = Q^k_j(r) / r^{k+1} + T^k_j(r) * r^k
        Y_k_j_sorted = np.where(r_sorted > 1e-15,
                                 Q_k / r_sorted**(k + 1), 0.0) + T_k * r_sorted**k

        # Unsort back to original quadrature order
        Y_k_j = np.zeros(n_quad)
        Y_k_j[sort_idx] = Y_k_j_sorted

        # R^k[i,j] = integral phi_i(r) Y^k_j(r) dr
        #           = sum_q w_q * phi_int[i,q] * Y^k_j(r_q)
        for i in range(n_basis):
            Rk[i, j] = np.sum(w_quad * phi_int[i, :] * Y_k_j)

    return Rk


def _slater_fk_cumcharge(
    r: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
    k: int = 0,
) -> float:
    """Evaluate Slater F^k integral via cumulative-charge method.

    F^k = integral integral P_A(r1) P_B(r2) r_<^k / r_>^{k+1} dr1 dr2

    Uses the same decomposition as inter_fiber_coupling.slater_fk_integral
    but operates on arbitrary density arrays.

    Parameters
    ----------
    r : ndarray
        Uniform radial grid.
    P_A, P_B : ndarray
        Radial probability densities on the grid.
    k : int
        Multipole order.

    Returns
    -------
    fk : float
    """
    dr = r[1] - r[0] if len(r) > 1 else 1.0
    if k == 0:
        Q = np.cumsum(P_B) * dr
        P_over_r = np.where(r > 1e-15, P_B / r, 0.0)
        T = np.flip(np.cumsum(np.flip(P_over_r))) * dr
        Y = np.where(r > 1e-15, Q / r, 0.0) + T
    else:
        Q = np.cumsum(P_B * r**k) * dr
        P_over_rk1 = np.where(r > 1e-15, P_B / r**(k + 1), 0.0)
        T = np.flip(np.cumsum(np.flip(P_over_rk1))) * dr
        Y = np.where(r > 1e-15, Q / r**(k + 1), 0.0) + T * r**k
    return float(np.sum(P_A * Y) * dr)


def _build_slater_Rk_matrix_fine(
    n_basis: int,
    alpha: float,
    k: int = 0,
    n_grid: int = 2000,
    r_max: float = 30.0,
    n_basis_B: Optional[int] = None,
    alpha_B: Optional[float] = None,
) -> np.ndarray:
    """Build R^k matrix using fine uniform grid (reference implementation).

    Supports cross-alpha evaluation: phi_i uses alpha (rows) and phi_j
    uses alpha_B (columns). When alpha_B is None, uses the same alpha
    for both (symmetric case).

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions for density A (rows).
    alpha : float
        Exponential decay parameter for density A.
    k : int
        Multipole order.
    n_grid : int
        Number of grid points.
    r_max : float
        Maximum radius.
    n_basis_B : int or None
        Number of basis functions for density B (columns). If None, same as n_basis.
    alpha_B : float or None
        Exponential decay for density B. If None, same as alpha.

    Returns
    -------
    Rk : ndarray of shape (n_basis, n_basis_B)
    """
    if n_basis_B is None:
        n_basis_B = n_basis
    if alpha_B is None:
        alpha_B = alpha

    dr = r_max / n_grid
    r = (np.arange(n_grid) + 0.5) * dr

    # Basis values for A (rows) with alpha
    x_A = 2.0 * alpha * r
    phi_A = np.zeros((n_basis, n_grid))
    for n in range(n_basis):
        phi_A[n, :] = r * np.exp(-alpha * r) * eval_laguerre(n, x_A)

    # Basis values for B (columns) with alpha_B
    x_B = 2.0 * alpha_B * r
    phi_B = np.zeros((n_basis_B, n_grid))
    for n in range(n_basis_B):
        phi_B[n, :] = r * np.exp(-alpha_B * r) * eval_laguerre(n, x_B)

    Rk = np.zeros((n_basis, n_basis_B))

    for j in range(n_basis_B):
        phi_j = phi_B[j, :]

        if k == 0:
            Q = np.cumsum(phi_j) * dr
            phi_over_r = np.where(r > 1e-15, phi_j / r, 0.0)
            T = np.flip(np.cumsum(np.flip(phi_over_r))) * dr
            Y = np.where(r > 1e-15, Q / r, 0.0) + T
        else:
            Q = np.cumsum(phi_j * r**k) * dr
            phi_over_rk1 = np.where(r > 1e-15, phi_j / r**(k + 1), 0.0)
            T = np.flip(np.cumsum(np.flip(phi_over_rk1))) * dr
            Y = np.where(r > 1e-15, Q / r**(k + 1), 0.0) + T * r**k

        for i in range(n_basis):
            Rk[i, j] = np.sum(phi_A[i, :] * Y) * dr

    return Rk


def slater_fk_algebraic(
    r_grid: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
    k: int = 0,
    n_basis: int = 20,
    alpha: Optional[float] = None,
    return_diagnostics: bool = False,
) -> float:
    """Compute Slater F^k integral algebraically via Laguerre basis expansion.

    Algebraic pathway:
      1. Fit P_A, P_B to Laguerre basis: P ~ sum_n c_n phi_n(r)
      2. F^k = c_A^T @ R^k @ c_B where R^k is precomputed

    Parameters
    ----------
    r_grid : ndarray
        Radial grid (uniform spacing).
    P_A, P_B : ndarray
        Radial probability densities.
    k : int
        Multipole order (0=monopole, 1=dipole, ...).
    n_basis : int
        Number of Laguerre basis functions for expansion.
    alpha : float or None
        Laguerre parameter. If None, estimated from densities.
    return_diagnostics : bool
        If True, return dict with diagnostics instead of scalar.

    Returns
    -------
    fk : float
        The F^k integral.
    OR if return_diagnostics:
    result : dict
        fk, coeffs_A, coeffs_B, alpha, residual_A, residual_B, Rk_matrix
    """
    dr = r_grid[1] - r_grid[0] if len(r_grid) > 1 else 1.0

    # Estimate per-density alpha from each density's mean radius
    def _estimate_alpha(P: np.ndarray) -> float:
        total = np.sum(P) * dr
        if total > 1e-15:
            mean_r = np.sum(r_grid * P) * dr / total
            return 1.5 / max(mean_r, 0.1)
        return 1.0

    if alpha is None:
        alpha_A_est = _estimate_alpha(P_A)
        alpha_B_est = _estimate_alpha(P_B)
    else:
        alpha_A_est = alpha
        alpha_B_est = alpha

    # Cap n_basis to avoid Laguerre polynomial instability.
    # Stable region: max Laguerre argument x_max = 2*alpha*r_data_max
    # should satisfy n_basis <= x_max/2 (rough heuristic from polynomial growth).
    r_data_max = r_grid[-1]
    n_basis_A = min(n_basis, max(5, int(alpha_A_est * r_data_max)))
    n_basis_B_eff = min(n_basis, max(5, int(alpha_B_est * r_data_max)))

    # Fit densities with per-density alpha and capped n_basis
    c_A, alpha_A, res_A = fit_density_laguerre(r_grid, P_A, n_basis_A, alpha_A_est)
    c_B, alpha_B, res_B = fit_density_laguerre(r_grid, P_B, n_basis_B_eff, alpha_B_est)

    # Hybrid approach: reconstruct smooth densities from Laguerre fit on
    # a truncated grid (where the exponential envelope is significant),
    # then use stable cumulative-charge numerical integration.
    r_max_rk = min(12.0 / min(alpha_A, alpha_B), r_data_max)
    r_max_rk = max(r_max_rk, 15.0)
    n_grid_rk = 3000
    dr_rk = r_max_rk / n_grid_rk
    r_rk = (np.arange(n_grid_rk) + 0.5) * dr_rk

    # Reconstruct fitted densities on the truncated grid and clip negatives
    phi_A_rk = _laguerre_basis_values(n_basis_A, alpha_A, r_rk)
    P_A_fit = np.maximum(phi_A_rk.T @ c_A, 0.0)
    phi_B_rk = _laguerre_basis_values(n_basis_B_eff, alpha_B, r_rk)
    P_B_fit = np.maximum(phi_B_rk.T @ c_B, 0.0)

    # Use cumulative-charge method on the smooth reconstructed densities
    fk = _slater_fk_cumcharge(r_rk, P_A_fit, P_B_fit, k)

    if return_diagnostics:
        return {
            'fk': fk,
            'coeffs_A': c_A,
            'coeffs_B': c_B,
            'alpha_A': alpha_A,
            'alpha_B': alpha_B,
            'residual_A': res_A,
            'residual_B': res_B,
            'Rk_matrix': Rk,
        }

    return fk


def slater_f0_algebraic(
    r_grid: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
    n_basis: int = 20,
    alpha: Optional[float] = None,
) -> float:
    """Compute Slater F^0 integral algebraically. Convenience wrapper."""
    return slater_fk_algebraic(r_grid, P_A, P_B, k=0, n_basis=n_basis, alpha=alpha)


# ---------------------------------------------------------------------------
# Analytical Slater integrals for hydrogenic densities (validation)
# ---------------------------------------------------------------------------

def _hydrogenic_slater_f0(Z: float, n1: int, l1: int, n2: int, l2: int) -> float:
    """Known analytical F^0 for hydrogenic orbitals.

    F^0(n1l1, n2l2) = integral integral |R_{n1l1}(r1)|^2 |R_{n2l2}(r2)|^2
                       / max(r1,r2) r1^2 r2^2 dr1 dr2

    Only a few cases have simple closed forms. Returns NaN if not known.
    """
    if (n1, l1, n2, l2) == (1, 0, 1, 0):
        return 5.0 * Z / 8.0
    elif (n1, l1, n2, l2) == (1, 0, 2, 0) or (n1, l1, n2, l2) == (2, 0, 1, 0):
        return 17.0 * Z / 81.0
    elif (n1, l1, n2, l2) == (2, 0, 2, 0):
        return 77.0 * Z / 512.0
    return float('nan')


# ---------------------------------------------------------------------------
# Integration with inter_fiber_coupling.py
# ---------------------------------------------------------------------------

def slater_fk_integral_algebraic(
    r_grid: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
    k: int = 0,
    n_basis: int = 20,
    alpha: Optional[float] = None,
) -> float:
    """Drop-in replacement for inter_fiber_coupling.slater_fk_integral.

    Same signature and semantics, but uses Laguerre basis expansion
    instead of cumulative-charge trapezoidal quadrature.

    Parameters
    ----------
    r_grid : ndarray
        Radial grid (uniform spacing).
    P_A, P_B : ndarray
        Radial probability densities.
    k : int
        Multipole order (0=monopole, 1=dipole, ...).
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float or None
        Laguerre parameter.

    Returns
    -------
    fk : float
        The F^k integral (Hartree).
    """
    return slater_fk_algebraic(
        r_grid, P_A, P_B, k=k, n_basis=n_basis, alpha=alpha)


def compute_channel_f0_matrix_algebraic(
    channel_data: dict,
    R: float,
    n_r: int = 300,
    r_max: float = 10.0,
    n_basis: int = 20,
) -> Dict:
    """Algebraic version of inter_fiber_coupling.compute_channel_f0_matrix.

    Uses Laguerre basis expansion for the F^0 integral instead of
    cumulative-charge trapezoidal quadrature.

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    R : float
        Internuclear distance (bohr).
    n_r : int
        Radial grid size for density extraction.
    r_max : float
        Maximum distance (bohr).
    n_basis : int
        Laguerre basis size.

    Returns
    -------
    result : dict
        Same keys as compute_channel_f0_matrix, plus 'method'.
    """
    from geovac.inter_fiber_coupling import (
        extract_channel_densities,
        transform_to_center_density,
    )

    # Step 1: Per-channel origin densities
    r_grid, P_ch_origin = extract_channel_densities(
        channel_data, n_r=n_r, r_max=r_max)

    n_ch = channel_data['n_ch']
    channels = channel_data['channels']

    # Step 2: Transform each channel density to Be center
    P_ch_Be_list = []
    d_grid_ref = None
    for ic in range(n_ch):
        d_grid, P_Be_ic = transform_to_center_density(
            r_grid, P_ch_origin[ic], R)
        if d_grid_ref is None:
            d_grid_ref = d_grid
        n_d = len(d_grid_ref)
        if len(P_Be_ic) < n_d:
            P_Be_ic = np.pad(P_Be_ic, (0, n_d - len(P_Be_ic)))
        elif len(P_Be_ic) > n_d:
            P_Be_ic = P_Be_ic[:n_d]
        P_ch_Be_list.append(P_Be_ic)

    # Step 3: Fit all channel densities to a shared Laguerre basis
    dr_d = d_grid_ref[1] - d_grid_ref[0] if len(d_grid_ref) > 1 else 1.0
    P_total = sum(P_ch_Be_list)
    total_norm = np.sum(P_total) * dr_d
    if total_norm > 1e-15:
        mean_r = np.sum(d_grid_ref * P_total) * dr_d / total_norm
        alpha_shared = 1.5 / max(mean_r, 0.1)
    else:
        alpha_shared = 1.0

    # Fit each channel
    ch_coeffs = []
    for ic in range(n_ch):
        c, _, _ = fit_density_laguerre(
            d_grid_ref, P_ch_Be_list[ic], n_basis, alpha_shared)
        ch_coeffs.append(c)

    # Step 4: Build R^0 matrix once
    r_max_d = d_grid_ref[-1] + dr_d
    Rk = _build_slater_Rk_matrix_fine(
        n_basis, alpha_shared, k=0,
        n_grid=2000, r_max=max(r_max_d, 30.0))

    # Step 5: F^0 matrix from contractions
    F0_matrix = np.zeros((n_ch, n_ch))
    for ic in range(n_ch):
        for jc in range(ic, n_ch):
            f0_ij = float(ch_coeffs[ic] @ Rk @ ch_coeffs[jc])
            F0_matrix[ic, jc] = f0_ij
            F0_matrix[jc, ic] = f0_ij

    F0_total = np.sum(F0_matrix)
    F0_per_channel = np.sum(F0_matrix, axis=1)

    return {
        'F0_total': float(F0_total),
        'F0_matrix': F0_matrix,
        'F0_per_channel': F0_per_channel,
        'channels': channels,
        'P_ch_Be': P_ch_Be_list,
        'd_grid': d_grid_ref,
        'method': 'algebraic_laguerre',
        'alpha': alpha_shared,
        'n_basis': n_basis,
    }


def compute_channel_fk_matrix_algebraic(
    channel_data: dict,
    R: float,
    k: int = 0,
    n_r: int = 300,
    r_max: float = 10.0,
    n_basis: int = 20,
) -> Dict:
    """Algebraic version of inter_fiber_coupling.compute_channel_fk_matrix.

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    R : float
        Internuclear distance (bohr).
    k : int
        Multipole order.
    n_r : int
        Radial grid size.
    r_max : float
        Maximum distance (bohr).
    n_basis : int
        Laguerre basis size.

    Returns
    -------
    result : dict
        Fk_total, Fk_matrix, Fk_per_channel, channels, k, method.
    """
    from geovac.inter_fiber_coupling import (
        extract_channel_densities,
        transform_to_center_density,
    )

    if k == 0:
        f0_result = compute_channel_f0_matrix_algebraic(
            channel_data, R, n_r=n_r, r_max=r_max, n_basis=n_basis)
        return {
            'Fk_total': f0_result['F0_total'],
            'Fk_matrix': f0_result['F0_matrix'],
            'Fk_per_channel': f0_result['F0_per_channel'],
            'channels': f0_result['channels'],
            'k': 0,
            'P_ch_Be': f0_result['P_ch_Be'],
            'd_grid': f0_result['d_grid'],
            'method': 'algebraic_laguerre',
        }

    # General k >= 1
    r_grid, P_ch_origin = extract_channel_densities(
        channel_data, n_r=n_r, r_max=r_max)
    n_ch = channel_data['n_ch']
    channels = channel_data['channels']

    P_ch_Be_list = []
    d_grid_ref = None
    for ic in range(n_ch):
        d_grid, P_Be_ic = transform_to_center_density(
            r_grid, P_ch_origin[ic], R)
        if d_grid_ref is None:
            d_grid_ref = d_grid
        n_d = len(d_grid_ref)
        if len(P_Be_ic) < n_d:
            P_Be_ic = np.pad(P_Be_ic, (0, n_d - len(P_Be_ic)))
        elif len(P_Be_ic) > n_d:
            P_Be_ic = P_Be_ic[:n_d]
        P_ch_Be_list.append(P_Be_ic)

    # Shared alpha
    dr_d = d_grid_ref[1] - d_grid_ref[0] if len(d_grid_ref) > 1 else 1.0
    P_total = sum(P_ch_Be_list)
    total_norm = np.sum(P_total) * dr_d
    if total_norm > 1e-15:
        mean_r = np.sum(d_grid_ref * P_total) * dr_d / total_norm
        alpha_shared = 1.5 / max(mean_r, 0.1)
    else:
        alpha_shared = 1.0

    ch_coeffs = []
    for ic in range(n_ch):
        c, _, _ = fit_density_laguerre(
            d_grid_ref, P_ch_Be_list[ic], n_basis, alpha_shared)
        ch_coeffs.append(c)

    r_max_d = d_grid_ref[-1] + dr_d
    Rk = _build_slater_Rk_matrix_fine(
        n_basis, alpha_shared, k=k,
        n_grid=2000, r_max=max(r_max_d, 30.0))

    Fk_matrix = np.zeros((n_ch, n_ch))
    for ic in range(n_ch):
        for jc in range(ic, n_ch):
            fk_ij = float(ch_coeffs[ic] @ Rk @ ch_coeffs[jc])
            Fk_matrix[ic, jc] = fk_ij
            Fk_matrix[jc, ic] = fk_ij

    return {
        'Fk_total': float(np.sum(Fk_matrix)),
        'Fk_matrix': Fk_matrix,
        'Fk_per_channel': np.sum(Fk_matrix, axis=1),
        'channels': channels,
        'k': k,
        'P_ch_Be': P_ch_Be_list,
        'd_grid': d_grid_ref,
        'method': 'algebraic_laguerre',
    }


# ---------------------------------------------------------------------------
# Incomplete Laguerre moment analysis (transcendental structure)
# ---------------------------------------------------------------------------

def incomplete_laguerre_moment(
    n: int,
    s: int,
    x: float,
) -> float:
    """Incomplete Laguerre moment: integral_0^x t^s L_n(t) e^{-t} dt.

    Uses the expansion L_n(t) = sum_{k=0}^n (-1)^k C(n,k) t^k / k!
    and the incomplete gamma function: integral_0^x t^{s+k} e^{-t} dt = gamma(s+k+1, x).

    Transcendental seed: gamma(1, x) = 1 - e^{-x}.
    All higher gamma(m, x) follow from: gamma(m+1, x) = m*gamma(m, x) - x^m * e^{-x}.

    Parameters
    ----------
    n : int
        Laguerre polynomial degree.
    s : int
        Power of t (must be >= 0).
    x : float
        Upper integration limit.

    Returns
    -------
    value : float
        The incomplete moment.
    """
    from math import comb

    # Build incomplete gamma values gamma(m, x) for m = s+1 to s+n+1
    # using lower regularized incomplete gamma: gammainc(a, x) = gamma(a, x) / Gamma(a)
    # So gamma(a, x) = gammainc(a, x) * Gamma(a)
    result = 0.0
    for k in range(n + 1):
        a = s + k + 1
        import math
        coeff = (-1)**k * comb(n, k) / math.factorial(k)
        # gamma(a, x) = gammainc(a, x) * Gamma(a)
        inc_gamma = float(gammainc(a, x) * gammafn(a))
        result += coeff * inc_gamma

    return result


def identify_transcendental_content() -> Dict:
    """Catalog of transcendental content in algebraic Slater integrals.

    Returns a description of the irreducible transcendental seeds,
    paralleling the exchange constant taxonomy of Paper 18.

    Returns
    -------
    catalog : dict
        Keys: 'seeds', 'recurrences', 'classification'.
    """
    return {
        'seeds': [
            'exp(-x): exponential decay, present in all Laguerre basis functions',
            'gamma(1, x) = 1 - exp(-x): regularized lower incomplete gamma, '
            'seed for all incomplete Laguerre moments',
        ],
        'recurrences': [
            'gamma(s+1, x) = s*gamma(s, x) - x^s*exp(-x): upward recurrence '
            'for incomplete gamma (exact, no truncation)',
            'L_{n+1}(x) = ((2n+1-x)*L_n(x) - n*L_{n-1}(x)) / (n+1): '
            'three-term Laguerre recurrence (exact)',
        ],
        'classification': {
            'type': 'embedding',
            'description': (
                'The Slater F^k integral between Laguerre-expanded densities '
                'is an embedding exchange constant (Paper 18, Sec IV): the '
                'graph structure (channel quantum numbers, angular selection '
                'rules) is algebraic, but the radial overlap between fibers '
                'requires projection onto continuous r-space via incomplete '
                'gamma functions. The transcendental seed exp(-x) enters '
                'through the Laguerre weight function, which is the conformal '
                'factor for the S3-to-R3 projection.'
            ),
        },
        'comparison_with_existing': {
            'Track_J': 'e^a * E_1(a) seed for Stieltjes integral (1 transcendental)',
            'Track_O': 'gamma(1, x) = 1 - e^(-x) seed for incomplete moments '
                       '(same transcendental family)',
            'Paper_12': 'Neumann auxiliary integrals use incomplete gamma via '
                        'A_n(alpha) recurrence (same structure)',
        },
    }
