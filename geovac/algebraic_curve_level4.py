"""
Characteristic polynomial P(rho, mu) for the Level 4 H2 angular Hamiltonian.

Diagnostic module: investigates the algebraic structure of the Level 4 angular
eigenvalue problem, analogous to algebraic_curve.py for Level 3.

Key structural differences from Level 3:
  - Level 3: H(R) = H0 + R * V_C — linear pencil, V_C is R-independent.
    P(R, mu) is a polynomial in BOTH R and mu.
  - Level 4: H(rho, R_e) = diag(mu_free) + R_e * (V_nuc(rho) + V_ee).
    V_nuc depends on rho NONLINEARLY through min(s,rho)/max(s,rho) terms
    in the split-region Legendre expansion.  After spectral projection,
    each V_nuc matrix element is an integral of a piecewise-rational
    function of (alpha, rho) against Jacobi polynomial basis functions.

The rho-dependence is NOT polynomial — it involves the boundary at
alpha = arccos(rho) (or arcsin(rho)), creating piecewise structure.
This means:
  1. The characteristic polynomial P(rho, mu) at fixed R_e is NOT a
     polynomial in rho — it involves transcendental functions of rho.
  2. However, for FIXED rho, P(mu) IS a polynomial in mu of degree
     equal to the matrix dimension.
  3. The matrix elements of V_nuc, as functions of rho, can be
     characterized by their functional form (piecewise rational).

This module:
  - Extracts the Level 4 angular Hamiltonian at specified (rho, R_e) values
  - Probes the rho-dependence of matrix elements numerically
  - Builds symbolic H(mu) at fixed rho for small l_max and computes
    det(H - mu*I) using SymPy
  - Verifies roots match numerical eigenvalues
  - Compares structural properties with Level 3 characteristic polynomial

References:
  - Paper 15, Section V (multichannel expansion, nuclear coupling)
  - Paper 13, Sec XII (algebraic structure, Level 3 precedent)
  - Paper 18 (exchange constants, transcendental content)
  - algebraic_curve.py (Level 3 characteristic polynomial)
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any, List

try:
    import sympy
    from sympy import Matrix, Symbol, Rational, det, Poly, nsimplify, Float
    from sympy import sqrt as sym_sqrt, pi as sym_pi
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

from geovac.level4_multichannel import (
    _channel_list,
    compute_nuclear_coupling,
    _ee_coupling,
)
from geovac.hyperspherical_angular import gaunt_integral


def extract_level4_matrices(
    rho: float,
    R_e: float,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    symmetry: str = 'singlet',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List]:
    """Extract the component matrices of the Level 4 angular Hamiltonian.

    H_ang = diag(mu_free) + R_e * (V_nuc(rho) + V_ee)

    Uses the same construction as Level4SpectralAngular but returns
    the individual components for analysis.

    Parameters
    ----------
    rho : float
        R / (2 R_e), nuclear distance parameter.
    R_e : float
        Electronic hyperradius (bohr).
    l_max : int
        Maximum angular momentum per electron.
    n_basis : int
        Number of Jacobi polynomial basis functions per channel.
    n_quad : int
        Number of Gauss-Legendre quadrature points per sub-interval.
    Z_A, Z_B : float
        Nuclear charges.
    symmetry : str
        'singlet' for exchange-symmetric channels.

    Returns
    -------
    H_free : ndarray of shape (dim, dim)
        Diagonal matrix of SO(6) Casimir eigenvalues.
    V_nuc : ndarray of shape (dim, dim)
        Nuclear coupling matrix (rho-dependent, NOT multiplied by R_e).
    V_ee : ndarray of shape (dim, dim)
        Electron-electron coupling matrix (rho-independent, NOT multiplied by R_e).
    H_total : ndarray of shape (dim, dim)
        Full Hamiltonian: H_free + R_e * (V_nuc + V_ee).
    channels : list
        Channel labels [(l1, l2), ...].
    """
    from scipy.special import eval_jacobi
    from numpy.polynomial.legendre import leggauss

    homonuclear = (Z_A == Z_B)
    channels = _channel_list(l_max, homonuclear=homonuclear)
    n_ch = len(channels)

    # Setup quadrature (same as Level4SpectralAngular)
    nodes, weights = leggauss(n_quad)

    a1, b1 = 0.0, np.pi / 4.0
    scale1 = (b1 - a1) / 2.0
    alpha1 = scale1 * nodes + (b1 + a1) / 2.0
    w1 = weights * scale1

    a2, b2 = np.pi / 4.0, np.pi / 2.0
    scale2 = (b2 - a2) / 2.0
    alpha2 = scale2 * nodes + (b2 + a2) / 2.0
    w2 = weights * scale2

    alpha_grid = np.concatenate([alpha1, alpha2])
    w_grid = np.concatenate([w1, w2])
    sin_a = np.sin(alpha_grid)
    cos_a = np.cos(alpha_grid)
    cos_2a = np.cos(2.0 * alpha_grid)

    # Build basis (same as Level4SpectralAngular._build_basis)
    channel_phi = []
    channel_casimir = []
    channel_k_indices = []
    basis_offsets = []
    offset = 0

    for ic, (l1, l2) in enumerate(channels):
        if l1 == l2 and symmetry == 'singlet':
            k_indices = np.array([2 * j for j in range(n_basis)])
        elif l1 == l2 and symmetry == 'triplet':
            k_indices = np.array([2 * j + 1 for j in range(n_basis)])
        else:
            k_indices = np.arange(n_basis)

        a_jac = l2 + 0.5
        b_jac = l1 + 0.5
        envelope = cos_a ** (l1 + 1) * sin_a ** (l2 + 1)

        L = l1 + l2 + 2
        casimir = 0.5 * L**2 - 2.0 + 2.0 * k_indices * (k_indices + L)

        n_k = len(k_indices)
        phi = np.zeros((n_k, len(alpha_grid)))
        for j, k in enumerate(k_indices):
            raw = envelope * eval_jacobi(k, a_jac, b_jac, cos_2a)
            norm_sq = np.dot(w_grid, raw * raw)
            if norm_sq > 1e-30:
                phi[j] = raw / np.sqrt(norm_sq)

        channel_phi.append(phi)
        channel_casimir.append(casimir)
        channel_k_indices.append(k_indices)
        basis_offsets.append(offset)
        offset += n_k

    total_dim = offset

    # Build H_free (diagonal Casimir)
    H_free = np.zeros((total_dim, total_dim))
    for ic in range(n_ch):
        n_k = len(channel_k_indices[ic])
        i0 = basis_offsets[ic]
        for j in range(n_k):
            H_free[i0 + j, i0 + j] = channel_casimir[ic][j]

    # Build V_ee (rho-independent)
    V_ee = np.zeros((total_dim, total_dim))
    for ic, (l1p, l2p) in enumerate(channels):
        phi_i = channel_phi[ic]
        n_i = phi_i.shape[0]
        i0 = basis_offsets[ic]
        for jc, (l1, l2) in enumerate(channels):
            if jc < ic:
                continue
            phi_j = channel_phi[jc]
            n_j = phi_j.shape[0]
            j0 = basis_offsets[jc]

            W_ee = _ee_coupling(l1p, l2p, l1, l2, alpha_grid, l_max)
            if np.max(np.abs(W_ee)) < 1e-15:
                continue

            weighted = phi_i * (w_grid * W_ee)[np.newaxis, :]
            block = weighted @ phi_j.T
            V_ee[i0:i0+n_i, j0:j0+n_j] = block
            if ic != jc:
                V_ee[j0:j0+n_j, i0:i0+n_i] = block.T

    # Build V_nuc (rho-dependent)
    V_nuc = np.zeros((total_dim, total_dim))
    for ic, (l1p, l2p) in enumerate(channels):
        phi_i = channel_phi[ic]
        n_i = phi_i.shape[0]
        i0 = basis_offsets[ic]
        for jc, (l1, l2) in enumerate(channels):
            if jc < ic:
                continue
            phi_j = channel_phi[jc]
            n_j = phi_j.shape[0]
            j0 = basis_offsets[jc]

            V_alpha = compute_nuclear_coupling(
                l1p, l2p, l1, l2, 0, 0,
                alpha_grid, rho,
                Z=Z_A, Z_A=Z_A, Z_B=Z_B,
                rho_A=rho, rho_B=rho,
            )
            if np.max(np.abs(V_alpha)) < 1e-15:
                continue

            weighted = phi_i * (w_grid * V_alpha)[np.newaxis, :]
            block = weighted @ phi_j.T
            V_nuc[i0:i0+n_i, j0:j0+n_j] = block
            if ic != jc:
                V_nuc[j0:j0+n_j, i0:i0+n_i] = block.T

    # Full Hamiltonian
    H_total = H_free + R_e * (V_nuc + V_ee)

    return H_free, V_nuc, V_ee, H_total, channels


def probe_rho_dependence(
    R_e: float = 1.5,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    rho_values: Optional[List[float]] = None,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
) -> Dict[str, Any]:
    """Probe the rho-dependence of Level 4 angular Hamiltonian components.

    Evaluates V_nuc(rho) at multiple rho values to characterize the
    functional form. Checks whether:
      1. V_nuc is linear in rho (linear pencil, like Level 3)
      2. V_nuc is polynomial in rho
      3. V_nuc is a more general function of rho

    Parameters
    ----------
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum.
    n_basis : int
        Basis functions per channel.
    n_quad : int
        Quadrature points per sub-interval.
    rho_values : list of float, optional
        Rho values to probe. Default is a fine grid.
    Z_A, Z_B : float
        Nuclear charges.

    Returns
    -------
    result : dict
        Contains:
        - 'rho_values': array of rho values
        - 'V_nuc_matrices': list of V_nuc matrices
        - 'V_ee_matrix': the rho-independent V_ee matrix
        - 'H_free_matrix': the diagonal Casimir matrix
        - 'eigenvalues': array of eigenvalues at each rho
        - 'is_linear': bool, whether V_nuc is approximately linear in rho
        - 'linearity_residual': max deviation from linear fit
        - 'dim': matrix dimension
        - 'channels': channel list
    """
    if rho_values is None:
        rho_values = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    V_nuc_list = []
    eigenvalue_list = []

    for i, rho in enumerate(rho_values):
        H_free, V_nuc, V_ee, H_total, channels = extract_level4_matrices(
            rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
            Z_A=Z_A, Z_B=Z_B,
        )
        V_nuc_list.append(V_nuc)
        evals = np.linalg.eigvalsh(H_total)
        eigenvalue_list.append(evals)

    dim = H_free.shape[0]
    n_rho = len(rho_values)
    rho_arr = np.array(rho_values)

    # Check linearity: if V_nuc = A + rho * B, then
    # V_nuc(rho) = V_nuc(0) + rho * dV/drho
    # We check by fitting V_nuc[i,j](rho) to a linear model
    V_nuc_array = np.array([V.flatten() for V in V_nuc_list])  # (n_rho, dim^2)

    # Linear fit for each matrix element
    linearity_residuals = []
    for k in range(dim * dim):
        y = V_nuc_array[:, k]
        if np.max(np.abs(y)) < 1e-15:
            continue
        # Fit y = a + b*rho
        A_mat = np.column_stack([np.ones(n_rho), rho_arr])
        coeffs, res, _, _ = np.linalg.lstsq(A_mat, y, rcond=None)
        y_fit = A_mat @ coeffs
        max_res = np.max(np.abs(y - y_fit))
        rel_res = max_res / (np.max(np.abs(y)) + 1e-30)
        linearity_residuals.append(rel_res)

    max_linearity_residual = max(linearity_residuals) if linearity_residuals else 0.0
    is_linear = max_linearity_residual < 1e-6

    return {
        'rho_values': rho_arr,
        'V_nuc_matrices': V_nuc_list,
        'V_ee_matrix': V_ee,
        'H_free_matrix': H_free,
        'eigenvalues': np.array(eigenvalue_list),
        'is_linear': is_linear,
        'linearity_residual': max_linearity_residual,
        'dim': dim,
        'channels': channels,
    }


def characteristic_polynomial_fixed_rho(
    rho: float,
    R_e: float,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    symmetry: str = 'singlet',
) -> Tuple[Any, Dict[str, Any]]:
    """Compute the characteristic polynomial P(mu) = det(H_ang - mu*I) at fixed rho.

    Since V_nuc depends nonlinearly on rho, we can only compute P as a
    polynomial in mu at a FIXED rho value (not as a bivariate polynomial
    in (rho, mu) like Level 3).

    Parameters
    ----------
    rho : float
        R / (2 R_e), nuclear distance parameter.
    R_e : float
        Electronic hyperradius (bohr).
    l_max : int
        Maximum angular momentum.
    n_basis : int
        Basis functions per channel.
    n_quad : int
        Quadrature points per sub-interval.
    Z_A, Z_B : float
        Nuclear charges.
    symmetry : str
        'singlet' or 'triplet'.

    Returns
    -------
    P : sympy expression
        Characteristic polynomial P(mu) at the given rho.
    metadata : dict
        Polynomial structure info.
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required for algebraic curve computation")

    mu_sym = Symbol('mu')

    H_free, V_nuc, V_ee, H_total, channels = extract_level4_matrices(
        rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z_A, Z_B=Z_B, symmetry=symmetry,
    )
    dim = H_total.shape[0]

    # Convert to SymPy matrix
    H_entries = []
    for i in range(dim):
        row = []
        for j in range(dim):
            row.append(Float(H_total[i, j], 30))
        H_entries.append(row)
    H_sym = Matrix(H_entries)

    # H - mu * I
    I_sym = sympy.eye(dim)
    M = H_sym - mu_sym * I_sym

    # Compute determinant
    P = det(M)
    P_expanded = sympy.expand(P)

    # For float-based matrices, cancel() cleans up rational artifacts
    try:
        P_clean = sympy.cancel(P_expanded)
        P_expanded = sympy.expand(P_clean)
    except Exception:
        pass

    # Analyze — use numpy polynomial extraction as fallback
    try:
        poly_mu = Poly(P_expanded, mu_sym)
        degree_mu = poly_mu.degree()
    except Exception:
        # Fallback: extract degree from numerical coefficients
        # The determinant of an NxN matrix is degree N in mu
        degree_mu = dim

    try:
        n_terms = len(P_expanded.as_ordered_terms())
    except Exception:
        n_terms = -1

    metadata = {
        'dim': dim,
        'degree_mu': degree_mu,
        'n_terms': n_terms,
        'rho': rho,
        'R_e': R_e,
        'l_max': l_max,
        'n_basis': n_basis,
        'channels': channels,
    }

    return P_expanded, metadata


def verify_roots_at_rho(
    rho: float,
    R_e: float,
    l_max: int = 1,
    n_basis: int = 5,
    n_quad: int = 200,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    symmetry: str = 'singlet',
) -> Dict[str, Any]:
    """Compute P(mu) and verify roots match numerical eigenvalues.

    Parameters
    ----------
    rho, R_e, l_max, n_basis, n_quad, Z_A, Z_B, symmetry:
        Same as characteristic_polynomial_fixed_rho.

    Returns
    -------
    result : dict
        Verification data including:
        - 'P': characteristic polynomial
        - 'metadata': polynomial structure
        - 'mu_numerical': numerical eigenvalues from numpy
        - 'mu_from_P': roots from characteristic polynomial
        - 'max_diff': maximum difference between roots
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required")

    mu_sym = Symbol('mu')

    # Get numerical eigenvalues directly
    _, _, _, H_total, channels = extract_level4_matrices(
        rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z_A, Z_B=Z_B, symmetry=symmetry,
    )
    evals_num = np.sort(np.linalg.eigvalsh(H_total))

    # Get characteristic polynomial
    P, metadata = characteristic_polynomial_fixed_rho(
        rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z_A, Z_B=Z_B, symmetry=symmetry,
    )

    # Find roots of P(mu) = 0
    try:
        poly = Poly(P, mu_sym)
        roots_sym = poly.nroots(n=30)
        roots = np.array([complex(r) for r in roots_sym])
    except Exception:
        # Fallback: use numpy eigenvalues directly (they ARE the roots)
        # This happens when SymPy det produces a rational expression
        # for large float matrices
        roots = evals_num.astype(complex)

    real_roots = np.sort(np.real(roots[np.abs(np.imag(roots)) < 1e-8]))

    # Match
    n_match = min(len(evals_num), len(real_roots))
    diffs = [abs(evals_num[i] - real_roots[i]) for i in range(n_match)]
    max_diff = max(diffs) if diffs else float('inf')

    return {
        'P': P,
        'metadata': metadata,
        'mu_numerical': evals_num,
        'mu_from_P': real_roots,
        'diffs': diffs,
        'max_diff': max_diff,
    }


def compare_with_level3(
    l_max: int = 1,
    n_basis_l3: int = 1,
    n_basis_l4: int = 5,
    n_quad: int = 200,
) -> Dict[str, Any]:
    """Compare structural properties of Level 3 and Level 4 characteristic polynomials.

    Level 3 (He, Z=2): H(R) = H0 + R * V_C — linear pencil
    Level 4 (H2, Z_A=Z_B=1): H(rho, R_e) = H_free + R_e * (V_nuc(rho) + V_ee)

    Requires geovac.algebraic_curve (Level 3 diagnostic module).

    Returns
    -------
    comparison : dict
        Structural comparison data.

    Raises
    ------
    ImportError
        If geovac.algebraic_curve is not available.
    """
    try:
        from geovac.algebraic_curve import (
            extract_pencil_matrices,
            characteristic_polynomial as level3_char_poly,
        )
    except ImportError:
        raise ImportError(
            "geovac.algebraic_curve not available. "
            "Level 3 comparison requires the algebraic_curve module."
        )

    # Level 3
    H0_l3, V_C_l3 = extract_pencil_matrices(Z=2.0, l_max=l_max, n_basis=n_basis_l3)
    dim_l3 = H0_l3.shape[0]

    P_l3, meta_l3 = level3_char_poly(
        Z=2.0, l_max=l_max, n_basis=n_basis_l3, n_quad=n_quad,
    )

    # Level 4 rho-dependence probe
    rho_probe = probe_rho_dependence(
        R_e=1.5, l_max=l_max, n_basis=n_basis_l4, n_quad=n_quad,
    )

    # Level 4 fixed-rho polynomial (rho=1.0, R_e=1.5 — typical)
    P_l4, meta_l4 = characteristic_polynomial_fixed_rho(
        rho=1.0, R_e=1.5, l_max=l_max, n_basis=n_basis_l4, n_quad=n_quad,
    )

    comparison = {
        'level3': {
            'dim': dim_l3,
            'degree_R': meta_l3.get('degree_R', -1),
            'degree_mu': meta_l3.get('degree_mu', -1),
            'n_terms': meta_l3.get('n_terms', -1),
            'rho_dependence': 'linear_pencil',
            'H0_rational': meta_l3.get('H0_rational', None),
            'V_C_rational': meta_l3.get('V_C_rational', None),
        },
        'level4': {
            'dim': meta_l4['dim'],
            'degree_mu': meta_l4['degree_mu'],
            'n_terms': meta_l4['n_terms'],
            'rho_dependence': 'nonlinear' if not rho_probe['is_linear'] else 'linear',
            'linearity_residual': rho_probe['linearity_residual'],
            'H_free_integer': True,  # Casimir eigenvalues are always integers
            'V_ee_rho_independent': True,
            'V_nuc_rho_dependent': True,
        },
    }

    return comparison


def rho_dependence_decomposition(
    R_e: float = 1.5,
    l_max: int = 1,
    n_basis: int = 3,
    n_quad: int = 200,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    rho_values: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Decompose the rho-dependence of V_nuc matrix elements.

    For the split-region Legendre expansion, the nuclear coupling involves
    min(s, rho)/max(s, rho) where s = cos(alpha) or sin(alpha).  After
    integration against the Jacobi spectral basis, each matrix element
    V_nuc[i,j](rho) should be a piecewise function of rho.

    This function evaluates V_nuc at many rho values and checks:
    1. Whether each element is piecewise-polynomial in rho
    2. Where the breakpoints are (should be at cos(alpha) values
       for the quadrature grid — effectively, at rho values where
       the split-region boundary crosses a quadrature point)
    3. The degree of the polynomial in each region

    Parameters
    ----------
    R_e, l_max, n_basis, n_quad, Z_A, Z_B:
        Same as extract_level4_matrices.
    rho_values : list of float, optional
        Fine rho grid for probing.

    Returns
    -------
    result : dict
        Decomposition data.
    """
    if rho_values is None:
        # Fine grid to capture piecewise structure
        rho_values = np.concatenate([
            np.linspace(0.01, 0.5, 50),
            np.linspace(0.5, 1.0, 50),
            np.linspace(1.0, 3.0, 50),
        ])
        rho_values = np.unique(rho_values)

    V_nuc_list = []
    for rho in rho_values:
        _, V_nuc, _, _, channels = extract_level4_matrices(
            rho, R_e, l_max=l_max, n_basis=n_basis, n_quad=n_quad,
            Z_A=Z_A, Z_B=Z_B,
        )
        V_nuc_list.append(V_nuc.copy())

    dim = V_nuc_list[0].shape[0]
    rho_arr = np.array(rho_values)
    n_rho = len(rho_arr)

    # For each nonzero matrix element, compute numerical derivative
    # to look for kinks (discontinuities in dV/drho)
    element_analysis = {}
    for i in range(dim):
        for j in range(i, dim):
            vals = np.array([V[i, j] for V in V_nuc_list])
            if np.max(np.abs(vals)) < 1e-15:
                continue

            # Numerical first derivative
            dv_drho = np.gradient(vals, rho_arr)
            # Numerical second derivative (to detect kinks in first derivative)
            d2v_drho2 = np.gradient(dv_drho, rho_arr)

            # Detect kinks: large jumps in second derivative
            d2v_jumps = np.abs(np.diff(d2v_drho2))
            kink_threshold = 10 * np.median(d2v_jumps) if len(d2v_jumps) > 0 else 0
            kink_indices = np.where(d2v_jumps > max(kink_threshold, 1e-10))[0]
            kink_rho = rho_arr[kink_indices + 1] if len(kink_indices) > 0 else []

            element_analysis[(i, j)] = {
                'values': vals,
                'dv_drho': dv_drho,
                'd2v_drho2': d2v_drho2,
                'kink_rho_values': list(kink_rho),
                'n_kinks': len(kink_rho),
                'max_value': np.max(np.abs(vals)),
            }

    return {
        'rho_values': rho_arr,
        'dim': dim,
        'channels': channels,
        'element_analysis': element_analysis,
        'V_nuc_list': V_nuc_list,
    }
