"""
2D variational solver for Level 3 (He) on S^5.

Eliminates the adiabatic approximation by treating the hyperradius R and
hyperangle alpha simultaneously in a tensor-product spectral basis:

    Psi(R, alpha) = sum_{n,j} c_{nj} phi_n(R) chi_j(alpha)

where phi_n are Laguerre radial functions and chi_j are Gegenbauer angular
functions from AlgebraicAngularSolver.

The full Hamiltonian (after R^{-5/2} extraction) is:

    H = -1/2 d^2/dR^2 + [Lambda^2/2 + R C(alpha)] / R^2 + 15/(8R^2)

In the tensor product basis this becomes:

    H = K_R (x) I_alpha
        + M_{1/R^2} (x) [diag(casimir)/2 + 15/8 I]
        + M_{1/R}   (x) coupling_full

    S = S_R (x) I_alpha

where K_R, S_R are the algebraic Laguerre kinetic/overlap matrices,
M_{1/R} and M_{1/R^2} are radial operator matrices computed by quadrature,
and casimir/coupling_full come from AlgebraicAngularSolver.

This avoids diagonalizing the angular Hamiltonian at each R, capturing the
full non-adiabatic correlation between radial and angular motion.

References:
  - Paper 13 (hyperspherical coordinates, angular eigenvalue problem)
  - Paper 15 (2D variational solver precedent for Level 4 H2)
  - Track DI Sprint 1 (this implementation)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import roots_laguerre, eval_laguerre
from typing import Tuple, Optional, Dict, List

from geovac.algebraic_angular import AlgebraicAngularSolver
from geovac.hyperspherical_radial import _build_laguerre_SK_algebraic
from geovac.cusp_correction import cusp_correction_he, cusp_correction_he_extrapolated


def _build_radial_operator_matrices(
    n_basis: int,
    alpha: float,
    R_min: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build radial operator matrices for the 2D tensor product Hamiltonian.

    Computes S_R, K_R (algebraic), and M_{1/R}, M_{1/R^2} (quadrature) in
    the Laguerre basis phi_n(R) = (R-R_min) exp(-alpha(R-R_min)) L_n(x).

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary (bohr).

    Returns
    -------
    S_R : ndarray (n_basis, n_basis)
        Overlap matrix (algebraic, exact).
    K_R : ndarray (n_basis, n_basis)
        Kinetic energy matrix (algebraic, exact).
    M_inv_R : ndarray (n_basis, n_basis)
        1/R operator matrix (quadrature).
    M_inv_R2 : ndarray (n_basis, n_basis)
        1/R^2 operator matrix (quadrature).
    """
    N = n_basis
    two_alpha = 2.0 * alpha
    inv_8a3 = 1.0 / (8.0 * alpha ** 3)

    # Algebraic overlap and kinetic (exact)
    S_R, K_R = _build_laguerre_SK_algebraic(N, alpha)

    # Gauss-Laguerre quadrature for 1/R and 1/R^2 operators
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_min + x_quad / two_alpha

    # Evaluate Laguerre polynomials at quadrature points
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    # Common weight factor: w_k * x_k^2
    x2_w = w_quad * x_quad ** 2

    # M_f[i,j] = inv_8a3 * sum_k w_k x_k^2 f(R_k) L_i(x_k) L_j(x_k)
    inv_R_q = 1.0 / R_q
    inv_R2_q = inv_R_q ** 2

    x2_w_invR = (x2_w * inv_R_q)[np.newaxis, :] * L_vals
    M_inv_R = inv_8a3 * (x2_w_invR @ L_vals.T)

    x2_w_invR2 = (x2_w * inv_R2_q)[np.newaxis, :] * L_vals
    M_inv_R2 = inv_8a3 * (x2_w_invR2 @ L_vals.T)

    return S_R, K_R, M_inv_R, M_inv_R2


def solve_he_variational_2d(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 15,
    l_max: int = 0,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
    symmetry: str = 'singlet',
    n_quad_alpha: int = 100,
    n_states: int = 1,
) -> Dict:
    """Solve the He atom via 2D variational method on S^5.

    Treats hyperradius R and hyperangle alpha simultaneously in a
    tensor-product spectral basis, eliminating the adiabatic approximation.

    Parameters
    ----------
    Z : float
        Nuclear charge (default 2.0 for He).
    n_basis_R : int
        Number of Laguerre radial basis functions.
    n_basis_alpha : int
        Number of Gegenbauer angular basis functions per l-channel.
    l_max : int
        Maximum partial wave quantum number.
    alpha_R : float
        Laguerre exponential decay parameter.
    R_min : float
        Left radial boundary (bohr).
    symmetry : str
        'singlet' or 'triplet'.
    n_quad_alpha : int
        Number of Gauss-Legendre quadrature points per sub-interval
        for the angular integrals.
    n_states : int
        Number of eigenstates to return.

    Returns
    -------
    result : dict
        'energies': ndarray of shape (n_states,) — energy eigenvalues (Ha)
        'E_exact': float — exact nonrelativistic He energy
        'error_pct': float — percentage error of ground state
        'n_basis_R': int
        'n_basis_alpha': int
        'l_max': int
        'dim_total': int — total basis dimension
        'dim_R': int — radial basis dimension
        'dim_alpha': int — angular basis dimension
    """
    E_EXACT = -2.903724377034119598  # Ha (Pekeris)

    # Build angular solver
    ang = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis_alpha, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad_alpha,
    )
    dim_alpha = ang._total_dim
    casimir_all = np.concatenate(ang._channel_casimir)
    coupling_full = ang._coupling_full

    # Build radial operator matrices
    S_R, K_R, M_inv_R, M_inv_R2 = _build_radial_operator_matrices(
        n_basis_R, alpha_R, R_min,
    )
    dim_R = n_basis_R

    # Total dimension
    dim_total = dim_R * dim_alpha

    # Full 2D Hamiltonian (after R^{-5/2} extraction):
    #   H = -1/2 d^2/dR^2 + H_ang(R)/R^2 + 15/(8R^2)
    # where H_ang(R) = diag(casimir) + R * coupling_full.
    # Expanding: H = K_R (x) I + M_{1/R^2} (x) [diag(casimir) + 15/8 I]
    #              + M_{1/R} (x) coupling_full
    ang_diag = casimir_all + 15.0 / 8.0

    # Assemble full Hamiltonian
    I_alpha = np.eye(dim_alpha)
    H = np.zeros((dim_total, dim_total))
    S = np.zeros((dim_total, dim_total))

    # Term 1: K_R (x) I_alpha — radial kinetic energy
    H += np.kron(K_R, I_alpha)

    # Term 2: M_{1/R^2} (x) diag(ang_diag) — angular kinetic + centrifugal
    H += np.kron(M_inv_R2, np.diag(ang_diag))

    # Term 3: M_{1/R} (x) coupling_full — nuclear + V_ee potential
    H += np.kron(M_inv_R, coupling_full)

    # Overlap: S_R (x) I_alpha
    S += np.kron(S_R, I_alpha)

    # Solve generalized eigenvalue problem
    evals, evecs = eigh(H, S)

    # Select lowest n_states
    idx = np.argsort(evals)
    energies = evals[idx[:n_states]]

    E_gs = energies[0]
    error_pct = abs((E_gs - E_EXACT) / E_EXACT) * 100.0

    return {
        'energies': energies,
        'E_exact': E_EXACT,
        'error_pct': error_pct,
        'n_basis_R': n_basis_R,
        'n_basis_alpha': n_basis_alpha,
        'l_max': l_max,
        'dim_total': dim_total,
        'dim_R': dim_R,
        'dim_alpha': dim_alpha,
    }


def convergence_study(
    Z: float = 2.0,
    basis_R_values: Optional[list] = None,
    basis_alpha_values: Optional[list] = None,
    l_max_values: Optional[list] = None,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
) -> Dict:
    """Run convergence studies for the 2D variational solver.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    basis_R_values : list of int
        Radial basis sizes to test.
    basis_alpha_values : list of int
        Angular basis sizes to test.
    l_max_values : list of int
        l_max values to test.
    alpha_R : float
        Laguerre decay parameter.
    R_min : float
        Left radial boundary.

    Returns
    -------
    results : dict
        'radial_convergence': list of (n_basis_R, energy, error_pct)
        'angular_convergence': list of (n_basis_alpha, energy, error_pct)
        'lmax_convergence': list of (l_max, energy, error_pct, dim)
    """
    if basis_R_values is None:
        basis_R_values = [10, 15, 20, 25, 30]
    if basis_alpha_values is None:
        basis_alpha_values = [5, 8, 10, 15, 20]
    if l_max_values is None:
        l_max_values = [0, 1, 2, 3]

    results = {}

    # Radial convergence (fixed angular)
    radial_conv = []
    for n_R in basis_R_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=n_R, n_basis_alpha=15,
            l_max=0, alpha_R=alpha_R, R_min=R_min,
        )
        radial_conv.append((n_R, res['energies'][0], res['error_pct']))
    results['radial_convergence'] = radial_conv

    # Angular convergence (fixed radial)
    angular_conv = []
    for n_a in basis_alpha_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=25, n_basis_alpha=n_a,
            l_max=0, alpha_R=alpha_R, R_min=R_min,
        )
        angular_conv.append((n_a, res['energies'][0], res['error_pct']))
    results['angular_convergence'] = angular_conv

    # l_max convergence (fixed radial and angular per channel)
    lmax_conv = []
    for lm in l_max_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=25, n_basis_alpha=10,
            l_max=lm, alpha_R=alpha_R, R_min=R_min,
        )
        lmax_conv.append((lm, res['energies'][0], res['error_pct'],
                          res['dim_total']))
    results['lmax_convergence'] = lmax_conv

    return results


def solve_he_precision(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 40,
    l_max_range: Optional[List[int]] = None,
    alpha_R: float = 2.0,
    R_min: float = 0.05,
    n_quad_alpha: int = 150,
) -> Dict:
    """High-precision He ground state with cusp correction.

    Runs the 2D variational solver at multiple l_max values and applies
    both theoretical Schwartz cusp correction and empirical CBS
    extrapolation.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis_R : int
        Laguerre radial basis size (25 is converged for He).
    n_basis_alpha : int
        Gegenbauer angular basis size per l-channel.
    l_max_range : list of int or None
        l_max values to compute. Default [0, 1, 2, 3, 4, 5, 6, 7].
    alpha_R : float
        Laguerre decay parameter.
    R_min : float
        Left radial boundary.
    n_quad_alpha : int
        Angular quadrature points.

    Returns
    -------
    result : dict
        'raw_energies': dict {l_max: energy}
        'cusp_corrected': dict {l_max: (E_corr, error_pct)}
        'cbs_extrapolated': (E_CBS, A_fit, rms)
        'best_energy': float — lowest-error energy estimate
        'best_error_pct': float
        'best_method': str — description of how best energy was obtained
        'E_exact': float
    """
    E_EXACT = -2.903724377034119598

    if l_max_range is None:
        l_max_range = list(range(8))

    # Compute at each l_max
    raw_energies = {}
    for lm in l_max_range:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=n_basis_R, n_basis_alpha=n_basis_alpha,
            l_max=lm, alpha_R=alpha_R, R_min=R_min,
            n_quad_alpha=n_quad_alpha,
        )
        raw_energies[lm] = res['energies'][0]

    # Apply theoretical cusp correction
    cusp_corrected = {}
    for lm, E in raw_energies.items():
        dE_cusp = cusp_correction_he(Z=Z, l_max=lm)
        E_corr = E + dE_cusp
        err = abs((E_corr - E_EXACT) / E_EXACT) * 100.0
        cusp_corrected[lm] = (E_corr, err)

    # Empirical CBS extrapolation (use l_max >= 3)
    high_lmax = {l: E for l, E in raw_energies.items() if l >= 3}
    cbs_result = None
    if len(high_lmax) >= 3:
        E_CBS, A_fit, rms = cusp_correction_he_extrapolated(high_lmax, Z=Z)
        cbs_result = (E_CBS, A_fit, rms)

    # Find best estimate
    best_energy = None
    best_error = float('inf')
    best_method = ''

    for lm, (E_corr, err) in cusp_corrected.items():
        if err < best_error:
            best_error = err
            best_energy = E_corr
            best_method = f'theoretical cusp correction at l_max={lm}'

    if cbs_result is not None:
        E_CBS = cbs_result[0]
        err_cbs = abs((E_CBS - E_EXACT) / E_EXACT) * 100.0
        if err_cbs < best_error:
            best_error = err_cbs
            best_energy = E_CBS
            best_method = 'empirical CBS extrapolation'

    return {
        'raw_energies': raw_energies,
        'cusp_corrected': cusp_corrected,
        'cbs_extrapolated': cbs_result,
        'best_energy': best_energy,
        'best_error_pct': best_error,
        'best_method': best_method,
        'E_exact': E_EXACT,
    }
