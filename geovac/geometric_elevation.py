"""
Geometric elevation analysis for Level 4 piecewise angular Hamiltonian.

Track Z investigation: can the piecewise structure of V_nuc (from
min(s,rho)/max(s,rho) split-region Legendre expansion) be resolved by
embedding in a higher-dimensional geometry?

Three avenues investigated, all NEGATIVE:
  1. Algebraic blow-up: resolves kink locally but integration limits are
     transcendental in rho -> no global polynomial P(rho, mu) = 0.
  2. Lie algebra elevation: min/max is C^0 not C^1, incompatible with
     matrix elements of any finite-dimensional Lie group representation.
  3. S3 x S3 hybrid: two-center nuclear potential cannot be smooth on any
     single compact manifold; equivalent to existing approach up to relabeling.

Master obstruction: the min/max boundary s = rho is PHYSICAL (electron
equidistant from origin and nucleus). It moves with rho, making the angular
domain partition intrinsically rho-dependent. No geometric elevation can
make this partition rho-independent.

This module provides computational verification of the obstructions:
  - Sheet decomposition of V_nuc on the blown-up space
  - Smoothness verification within each sheet
  - Transcendence of matrix elements as functions of rho
  - Lie representation incompatibility test (C^0 not C^1 verification)

References:
  - Paper 15, Section V.A (split-region Legendre expansion)
  - Paper 18 (exchange constant classification)
  - Track S (algebraic_curve_level4.py, linearity residuals)
  - Track W (cusp_graph.py, S^5 dimensionality obstruction)
"""

import numpy as np
from typing import Tuple, Dict, Any, List, Optional


def split_region_sheets(
    alpha_grid: np.ndarray,
    rho: float,
) -> Dict[str, Any]:
    """Decompose the alpha domain into sheets based on min/max boundaries.

    For electron 1 (s1 = cos(alpha)): boundary at alpha_1 = arccos(rho)
    For electron 2 (s2 = sin(alpha)): boundary at alpha_2 = arcsin(rho)

    Creates 4 sheets (when 0 < rho < 1):
      Sheet (++): s1 >= rho AND s2 >= rho  (alpha near pi/4 for rho < 1/sqrt(2))
      Sheet (+-): s1 >= rho AND s2 < rho
      Sheet (-+): s1 < rho AND s2 >= rho
      Sheet (--): s1 < rho AND s2 < rho

    When rho >= 1, s1 = cos(alpha) <= 1 < rho always, so electron 1 is
    always in the s < rho regime. Similarly for rho <= 0.

    Parameters
    ----------
    alpha_grid : ndarray
        Alpha grid points in [0, pi/2].
    rho : float
        Nuclear distance parameter R/(2 R_e).

    Returns
    -------
    result : dict
        Sheet decomposition data:
        - 'alpha_1': boundary for electron 1 (arccos(rho) or None)
        - 'alpha_2': boundary for electron 2 (arcsin(rho) or None)
        - 'n_sheets': number of distinct sheets (1, 2, or 4)
        - 'sheet_masks': dict of boolean masks for each sheet
        - 'sheet_names': list of sheet names
    """
    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    # Boundary for electron 1: cos(alpha) = rho
    if 0 < rho < 1:
        alpha_1 = np.arccos(rho)
    else:
        alpha_1 = None

    # Boundary for electron 2: sin(alpha) = rho
    if 0 < rho < 1:
        alpha_2 = np.arcsin(rho)
    else:
        alpha_2 = None

    # Build masks
    if alpha_1 is not None and alpha_2 is not None:
        mask_s1_ge = cos_a >= rho  # s1 >= rho
        mask_s2_ge = sin_a >= rho  # s2 >= rho
        sheets = {
            '++': mask_s1_ge & mask_s2_ge,
            '+-': mask_s1_ge & ~mask_s2_ge,
            '-+': ~mask_s1_ge & mask_s2_ge,
            '--': ~mask_s1_ge & ~mask_s2_ge,
        }
        # Remove empty sheets
        sheets = {k: v for k, v in sheets.items() if np.any(v)}
        n_sheets = len(sheets)
    elif alpha_1 is not None:
        mask_s1_ge = cos_a >= rho
        sheets = {'+': mask_s1_ge, '-': ~mask_s1_ge}
        sheets = {k: v for k, v in sheets.items() if np.any(v)}
        n_sheets = len(sheets)
    elif alpha_2 is not None:
        mask_s2_ge = sin_a >= rho
        sheets = {'+': mask_s2_ge, '-': ~mask_s2_ge}
        sheets = {k: v for k, v in sheets.items() if np.any(v)}
        n_sheets = len(sheets)
    else:
        sheets = {'all': np.ones(len(alpha_grid), dtype=bool)}
        n_sheets = 1

    return {
        'alpha_1': alpha_1,
        'alpha_2': alpha_2,
        'n_sheets': n_sheets,
        'sheet_masks': sheets,
        'sheet_names': list(sheets.keys()),
    }


def verify_within_sheet_smoothness(
    rho: float,
    R_e: float = 1.5,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    n_rho_test: int = 20,
    delta_rho: float = 0.01,
) -> Dict[str, Any]:
    """Verify that V_nuc matrix elements are smooth WITHIN each sheet.

    For a fixed sheet (fixed sign pattern), V_nuc should be a smooth
    (actually polynomial) function of rho. This is verified by checking
    that finite-difference derivatives up to order 3 are smooth.

    The test uses a narrow rho range around the given rho value,
    staying within a single sheet.

    Parameters
    ----------
    rho : float
        Central rho value.
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum.
    n_basis : int
        Basis functions per channel.
    n_quad : int
        Quadrature points.
    n_rho_test : int
        Number of rho test points.
    delta_rho : float
        Half-width of rho range.

    Returns
    -------
    result : dict
        Smoothness analysis.
    """
    from geovac.algebraic_curve_level4 import extract_level4_matrices

    # Stay within a single sheet by using small rho range
    rho_test = np.linspace(rho - delta_rho, rho + delta_rho, n_rho_test)
    rho_test = rho_test[rho_test > 0.01]  # avoid rho=0

    V_nuc_list = []
    for r in rho_test:
        _, V_nuc, _, _, channels = extract_level4_matrices(
            r, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        )
        V_nuc_list.append(V_nuc.copy())

    dim = V_nuc_list[0].shape[0]

    # Check smoothness via polynomial fit residual
    max_poly_residuals = {}
    for deg in [1, 2, 3]:
        max_res = 0.0
        for i in range(dim):
            for j in range(i, dim):
                vals = np.array([V[i, j] for V in V_nuc_list])
                if np.max(np.abs(vals)) < 1e-15:
                    continue
                # Polynomial fit
                coeffs = np.polyfit(rho_test, vals, deg)
                fit = np.polyval(coeffs, rho_test)
                res = np.max(np.abs(vals - fit)) / (np.max(np.abs(vals)) + 1e-30)
                max_res = max(max_res, res)
        max_poly_residuals[deg] = max_res

    return {
        'rho_center': rho,
        'rho_range': (rho_test[0], rho_test[-1]),
        'n_rho': len(rho_test),
        'poly_residuals': max_poly_residuals,
        'is_smooth_within_sheet': max_poly_residuals.get(3, 1.0) < 1e-4,
    }


def verify_kink_at_boundary(
    R_e: float = 1.5,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    n_rho: int = 100,
) -> Dict[str, Any]:
    """Verify derivative discontinuity at the sheet boundary.

    Scans V_nuc matrix elements across rho values that cross the
    boundary cos(alpha*) = rho (i.e., alpha* = pi/4 corresponds to
    rho* = cos(pi/4) = 1/sqrt(2) ~= 0.707). The kink appears in the
    first derivative dV/drho.

    Parameters
    ----------
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum.
    n_basis : int
        Basis functions per channel.
    n_quad : int
        Quadrature points.
    n_rho : int
        Number of rho scan points.

    Returns
    -------
    result : dict
        Kink analysis data.
    """
    from geovac.algebraic_curve_level4 import extract_level4_matrices

    # Scan through key boundary region
    # For alpha in [0, pi/2], the boundaries are:
    #   cos(alpha) = rho at alpha = arccos(rho)
    #   sin(alpha) = rho at alpha = arcsin(rho)
    # These cross at alpha = pi/4, rho = 1/sqrt(2)
    rho_scan = np.linspace(0.1, 1.2, n_rho)

    V_nuc_list = []
    for rho in rho_scan:
        _, V_nuc, _, _, _ = extract_level4_matrices(
            rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        )
        V_nuc_list.append(V_nuc.copy())

    dim = V_nuc_list[0].shape[0]

    # Compute first and second derivatives for the (0,0) element
    element_derivatives = {}
    for i in range(min(dim, 3)):
        for j in range(i, min(dim, 3)):
            vals = np.array([V[i, j] for V in V_nuc_list])
            if np.max(np.abs(vals)) < 1e-15:
                continue
            dv = np.gradient(vals, rho_scan)
            d2v = np.gradient(dv, rho_scan)

            # Detect kinks: jumps in d2v
            d2v_jumps = np.abs(np.diff(d2v))
            median_jump = np.median(d2v_jumps)
            max_jump = np.max(d2v_jumps)
            kink_ratio = max_jump / (median_jump + 1e-30)

            element_derivatives[(i, j)] = {
                'values': vals,
                'first_deriv': dv,
                'second_deriv': d2v,
                'max_d2v_jump': max_jump,
                'median_d2v_jump': median_jump,
                'kink_ratio': kink_ratio,
                'has_kink': kink_ratio > 5.0,
            }

    return {
        'rho_scan': rho_scan,
        'element_derivatives': element_derivatives,
        'rho_boundary': 1.0 / np.sqrt(2.0),
    }


def verify_lie_representation_obstruction(
    n_alpha: int = 1000,
    rho_values: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Verify that min/max is C^0 but not C^1, ruling out Lie group matrix elements.

    Matrix elements of finite-dimensional Lie group representations are
    real-analytic. We verify that f_k(s, rho) = min(s,rho)^k / max(s,rho)^{k+1}
    is NOT C^1 at s = rho by computing left and right derivatives.

    Parameters
    ----------
    n_alpha : int
        Number of alpha grid points for evaluation.
    rho_values : list of float, optional
        Rho values to test. Default: [0.3, 0.5, 0.707, 0.9].

    Returns
    -------
    result : dict
        Derivative discontinuity data for each k and rho.
    """
    if rho_values is None:
        rho_values = [0.3, 0.5, 1.0 / np.sqrt(2.0), 0.9]

    results = {}
    for rho in rho_values:
        for k in [0, 1, 2]:
            # f_k(s, rho) as function of s
            s_left = np.linspace(max(rho - 0.2, 0.01), rho - 1e-8, 500)
            s_right = np.linspace(rho + 1e-8, min(rho + 0.2, 0.99), 500)

            # s < rho: f_k = s^k / rho^{k+1}
            f_left = s_left**k / rho**(k + 1)
            # s > rho: f_k = rho^k / s^{k+1}
            f_right = rho**k / s_right**(k + 1)

            # Left derivative: df/ds at s -> rho^-
            # d/ds [s^k / rho^{k+1}] = k * s^{k-1} / rho^{k+1}
            if k > 0:
                deriv_left = k * rho**(k - 1) / rho**(k + 1)  # = k / rho^2
            else:
                deriv_left = 0.0

            # Right derivative: df/ds at s -> rho^+
            # d/ds [rho^k / s^{k+1}] = -(k+1) * rho^k / s^{k+2}
            deriv_right = -(k + 1) * rho**k / rho**(k + 2)  # = -(k+1) / rho^2

            # Derivative jump
            deriv_jump = abs(deriv_right - deriv_left)

            # Continuity check: both limits -> f(rho, rho) = 1/rho
            f_at_boundary = 1.0 / rho
            f_left_limit = rho**k / rho**(k + 1)  # = 1/rho
            f_right_limit = rho**k / rho**(k + 1)  # = 1/rho

            results[(rho, k)] = {
                'f_continuous': abs(f_left_limit - f_right_limit) < 1e-14,
                'deriv_left': deriv_left,
                'deriv_right': deriv_right,
                'deriv_jump': deriv_jump,
                'is_C1': deriv_jump < 1e-10,
                'is_C0_not_C1': abs(f_left_limit - f_right_limit) < 1e-14
                               and deriv_jump > 1e-10,
            }

    return {
        'rho_values': rho_values,
        'k_values': [0, 1, 2],
        'results': results,
        'all_C0_not_C1': all(
            v['is_C0_not_C1'] for v in results.values()
        ),
    }


def verify_transcendental_matrix_elements(
    R_e: float = 1.5,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
) -> Dict[str, Any]:
    """Verify that blown-up sheet integrals produce transcendental rho-dependence.

    On the blown-up space (Avenue 1), matrix elements within a sheet are:
      M_{kj}(rho) = integral_{0}^{arccos(rho)} phi_k(a) * [rho^p / cos(a)^q] * phi_j(a) da

    The upper limit arccos(rho) makes M_{kj} a transcendental function of rho.
    This is verified by checking that polynomial fits of increasing degree
    do NOT converge to machine precision.

    Parameters
    ----------
    R_e, l_max, n_basis, n_quad : as for extract_level4_matrices.

    Returns
    -------
    result : dict
        Polynomial fit residuals showing transcendental behavior.
    """
    from geovac.algebraic_curve_level4 import extract_level4_matrices

    # Scan V_nuc across a range of rho values within a single sheet
    # Use rho < 0.5 to stay in the regime where most of [0, pi/2]
    # has cos(alpha) > rho (i.e., electron 1 farther than nucleus)
    rho_scan = np.linspace(0.05, 0.45, 40)

    V_nuc_list = []
    for rho in rho_scan:
        _, V_nuc, _, _, _ = extract_level4_matrices(
            rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        )
        V_nuc_list.append(V_nuc.copy())

    dim = V_nuc_list[0].shape[0]

    # For the largest non-zero element, fit polynomials of increasing degree
    # and check residual convergence
    best_element = None
    best_max_val = 0.0
    for i in range(dim):
        for j in range(i, dim):
            vals = np.array([V[i, j] for V in V_nuc_list])
            mv = np.max(np.abs(vals))
            if mv > best_max_val:
                best_max_val = mv
                best_element = (i, j)
                best_vals = vals.copy()

    if best_element is None:
        return {'error': 'No nonzero matrix elements found'}

    # Polynomial fit at increasing degree
    poly_residuals = {}
    for deg in range(1, 16):
        coeffs = np.polyfit(rho_scan, best_vals, deg)
        fit = np.polyval(coeffs, rho_scan)
        res = np.max(np.abs(best_vals - fit)) / best_max_val
        poly_residuals[deg] = res

    # Transcendental: residuals decay slower than exponential in degree
    # Polynomial: residuals hit machine precision at degree = true degree
    residual_values = [poly_residuals[d] for d in sorted(poly_residuals.keys())]
    hits_machine_eps = any(r < 1e-12 for r in residual_values)

    return {
        'element': best_element,
        'max_value': best_max_val,
        'poly_residuals': poly_residuals,
        'is_polynomial': hits_machine_eps,
        'is_transcendental': not hits_machine_eps,
    }


def assess_blowup_cost(
    l_max: int = 2,
    n_basis: int = 10,
    n_rho: int = 130,
) -> Dict[str, Any]:
    """Estimate computational cost of algebraic blow-up vs current approach.

    Current approach: n_rho eigensolves of dim x dim matrices
    Blow-up approach: n_rho evaluations on 4 sheets with rho-dependent limits,
    then eigensolves of the same dimension.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum.
    n_basis : int
        Basis functions per channel.
    n_rho : int
        Number of rho grid points.

    Returns
    -------
    result : dict
        Cost comparison.
    """
    from geovac.level4_multichannel import _channel_list

    channels = _channel_list(l_max, homonuclear=True)
    n_ch = len(channels)
    dim = n_ch * n_basis

    # Current cost: n_rho * O(dim^3) eigensolves + n_rho * O(dim^2) V_nuc builds
    current_eigensolve_cost = n_rho * dim**3
    current_vnuc_cost = n_rho * dim**2

    # Blow-up cost: same eigensolves (cannot be avoided) +
    # 4 * n_rho * O(dim^2) V_nuc builds (one per sheet) +
    # n_rho * O(dim) boundary calculations
    blowup_vnuc_cost = 4 * n_rho * dim**2
    blowup_boundary_cost = n_rho * dim

    return {
        'l_max': l_max,
        'n_basis': n_basis,
        'n_channels': n_ch,
        'matrix_dim': dim,
        'n_rho': n_rho,
        'current_eigensolve_flops': current_eigensolve_cost,
        'current_vnuc_flops': current_vnuc_cost,
        'blowup_vnuc_flops': blowup_vnuc_cost,
        'blowup_boundary_flops': blowup_boundary_cost,
        'blowup_overhead_factor': blowup_vnuc_cost / current_vnuc_cost,
        'eigensolve_dominates': current_eigensolve_cost > 10 * current_vnuc_cost,
        'blowup_speedup': 'none',
        'reason': (
            'Eigensolves dominate cost and cannot be eliminated. '
            'Blow-up adds 4x V_nuc overhead without reducing eigensolve count. '
            'Net effect: 0-4x SLOWER, not faster.'
        ),
    }
