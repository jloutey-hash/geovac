"""
Characteristic polynomial P(R, mu) for the Level 3 He angular Hamiltonian.

Diagnostic module: computes det(H_ang(R) - mu*I) symbolically using SymPy,
where H_ang(R) = H0 + R*V_C is the linear matrix pencil from the algebraic
angular solver.

The characteristic polynomial P(R, mu) = 0 defines an algebraic curve in
(R, mu) space.  Unlike the Taylor series of mu(R) (Track G, which diverges
for R > ~5 bohr), the implicit relation P(R, mu) = 0 is exact for all R.

Key findings:
  - H0 has INTEGER entries: SO(6) Casimir eigenvalues 2(l+k+1)^2 - 2.
  - V_C entries are NOT rational: they involve pi and sqrt(2).
    Nuclear coupling (diagonal in l): nuc[0,0] = -Z*32/(3*pi),
    nuc[1,1] = -Z*1024/(105*pi).
    V_ee coupling: vee[0,0] = 8*sqrt(2)/(3*pi),
    vee[0,1] = 16*sqrt(2)/(15*pi), vee[0,2] = 32/(35*pi).
  - Therefore P(R,mu) does NOT have rational coefficients; its coefficients
    lie in Q(pi, sqrt(2)) — a transcendental extension of Q.
  - Despite this, P(R,mu) = 0 exactly determines mu(R) as a function of R
    for ALL R, including R > 5 where Track G's Taylor series diverges.
  - The transcendental content enters through the coupling matrix V_C (the
    "exchange constants" of Paper 18), not through the algebraic structure
    of the characteristic polynomial itself.

Polynomial structure (n_basis=1 per channel):
  l_max=0 (dim 1): P = V_C[0,0]*R - mu, degree 1 in both R and mu.
  l_max=1 (dim 2): degree 2 in R and mu, 5 terms.
  l_max=2 (dim 3): degree 3 in R and mu, 9 terms.
  l_max=L (dim L+1): degree L+1 in both R and mu.
  l_max=2, n_basis=2 (dim 6): degree 6 in both, 27 terms.

References:
  - Paper 13, Sec III (angular Hamiltonian), Sec XII (algebraic structure)
  - Paper 13, Sec XII.B (transcendence of mu(R) as a function)
  - Paper 18, Sec II.C (mu(R) as flow constant)
  - Track G (perturbation series divergence, v2.0.9)
"""

import numpy as np
from scipy.linalg import eigh
from typing import Tuple, Optional, Dict, Any, List

try:
    import sympy
    from sympy import Matrix, Symbol, Rational, det, Poly, nsimplify, Float
    from sympy import sqrt as sym_sqrt
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

from geovac.algebraic_angular import AlgebraicAngularSolver


def extract_pencil_matrices(
    Z: float = 2.0,
    l_max: int = 0,
    n_basis: int = 1,
    symmetry: str = 'singlet',
    n_quad: int = 200,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract H0 and V_C matrices from the AlgebraicAngularSolver.

    The angular Hamiltonian is H(R) = H0 + R * V_C where:
      - H0 = diag(casimir_all) — SO(6) Casimir eigenvalues (integers)
      - V_C = coupling_full — nuclear + V_ee coupling (from quadrature)

    Parameters
    ----------
    Z : float
        Nuclear charge (default 2.0 for He).
    l_max : int
        Maximum partial wave.
    n_basis : int
        Number of basis functions per l-channel.
    symmetry : str
        'singlet' or 'triplet'.
    n_quad : int
        Quadrature points (higher = more accurate coupling matrix).

    Returns
    -------
    H0 : ndarray of shape (dim, dim)
        Diagonal matrix of SO(6) Casimir eigenvalues.
    V_C : ndarray of shape (dim, dim)
        R-independent coupling matrix (nuclear + V_ee).
    """
    solver = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad,
    )
    casimir_all = np.concatenate(solver._channel_casimir)
    H0 = np.diag(casimir_all)
    V_C = solver._coupling_full.copy()
    return H0, V_C


def rationalize_matrix(
    M: np.ndarray,
    tol: float = 1e-10,
    max_denom: int = 10000,
) -> 'Matrix':
    """Convert a numerical matrix to a SymPy Matrix with rational entries.

    Uses nsimplify to find exact rational representations of floating-point
    matrix elements.  Verifies that all entries are successfully rationalized.

    Parameters
    ----------
    M : ndarray
        Numerical matrix.
    tol : float
        Tolerance for rationalization.
    max_denom : int
        Maximum denominator for rational approximation.

    Returns
    -------
    M_sym : sympy.Matrix
        Matrix with Rational entries.

    Raises
    ------
    ValueError
        If any entry cannot be rationalized within tolerance.
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required for algebraic curve computation")

    n, m = M.shape
    entries = []
    for i in range(n):
        row = []
        for j in range(m):
            val = float(M[i, j])
            if abs(val) < tol:
                row.append(Rational(0))
            else:
                # Try to rationalize
                rat = nsimplify(val, rational=True, tolerance=tol)
                # Verify
                diff = abs(float(rat) - val)
                if diff > tol:
                    # Try harder with explicit rational_conversion
                    rat = Rational(val).limit_denominator(max_denom)
                    diff = abs(float(rat) - val)
                    if diff > tol:
                        raise ValueError(
                            f"Cannot rationalize M[{i},{j}] = {val:.15e} "
                            f"(best rational: {rat}, diff: {diff:.2e})"
                        )
                row.append(rat)
        entries.append(row)
    return Matrix(entries)


def characteristic_polynomial(
    Z: float = 2.0,
    l_max: int = 0,
    n_basis: int = 1,
    symmetry: str = 'singlet',
    n_quad: int = 200,
    rationalize: bool = True,
    tol: float = 1e-10,
) -> Tuple[Any, Dict[str, Any]]:
    """Compute the characteristic polynomial P(R, mu) = det(H0 + R*V_C - mu*I).

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_basis : int
        Number of basis functions per l-channel.
    symmetry : str
        'singlet' or 'triplet'.
    n_quad : int
        Quadrature points for coupling matrix.
    rationalize : bool
        If True, attempt to rationalize matrix entries.
    tol : float
        Tolerance for rationalization.

    Returns
    -------
    P : sympy expression
        The characteristic polynomial P(R, mu) as a SymPy expression.
    metadata : dict
        Dictionary with polynomial structure information:
        - 'dim': matrix dimension
        - 'degree_R': degree in R
        - 'degree_mu': degree in mu
        - 'n_terms': number of terms in expanded polynomial
        - 'H0_rational': whether H0 has rational entries
        - 'V_C_rational': whether V_C has rational entries
        - 'P_rational': whether P has rational coefficients
        - 'H0_sym': SymPy matrix for H0
        - 'V_C_sym': SymPy matrix for V_C
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required for algebraic curve computation")

    R_sym = Symbol('R')
    mu_sym = Symbol('mu')

    H0_np, V_C_np = extract_pencil_matrices(
        Z=Z, l_max=l_max, n_basis=n_basis,
        symmetry=symmetry, n_quad=n_quad,
    )

    dim = H0_np.shape[0]

    # H0 should have integer entries (SO(6) Casimir)
    H0_rational = True
    H0_entries = []
    for i in range(dim):
        row = []
        for j in range(dim):
            val = H0_np[i, j]
            int_val = int(round(val))
            if abs(val - int_val) < 1e-12:
                row.append(Rational(int_val))
            else:
                H0_rational = False
                row.append(nsimplify(val, rational=True, tolerance=tol))
        H0_entries.append(row)
    H0_sym = Matrix(H0_entries)

    # V_C: try to rationalize
    V_C_rational = True
    try:
        V_C_sym = rationalize_matrix(V_C_np, tol=tol)
    except ValueError:
        V_C_rational = False
        # Fall back to high-precision floats
        V_C_entries = []
        for i in range(dim):
            row = []
            for j in range(dim):
                row.append(Float(V_C_np[i, j], 30))
            V_C_entries.append(row)
        V_C_sym = Matrix(V_C_entries)

    # Build H(R) - mu*I
    I_sym = sympy.eye(dim)
    H_sym = H0_sym + R_sym * V_C_sym - mu_sym * I_sym

    # Compute determinant
    P = det(H_sym)
    P_expanded = sympy.expand(P)

    # Analyze polynomial structure
    # Degree in mu: should be dim
    poly_mu = Poly(P_expanded, mu_sym) if V_C_rational else None
    degree_mu = dim  # Always equals matrix dimension

    # Degree in R: at most dim (each row contributes at most one R factor)
    poly_R = None
    try:
        poly_R = Poly(P_expanded, R_sym)
        degree_R = poly_R.degree()
    except Exception:
        degree_R = dim  # Fallback

    # Number of terms
    try:
        n_terms = len(P_expanded.as_ordered_terms())
    except Exception:
        n_terms = -1

    # Check if all coefficients of the bivariate polynomial are rational
    P_rational = V_C_rational and H0_rational

    metadata = {
        'dim': dim,
        'degree_R': degree_R,
        'degree_mu': degree_mu,
        'n_terms': n_terms,
        'H0_rational': H0_rational,
        'V_C_rational': V_C_rational,
        'P_rational': P_rational,
        'H0_sym': H0_sym,
        'V_C_sym': V_C_sym,
    }

    return P_expanded, metadata


def solve_characteristic_polynomial(
    P: Any,
    R_val: float,
) -> np.ndarray:
    """Solve P(R_val, mu) = 0 for mu at a given R value.

    Substitutes R = R_val into the characteristic polynomial and finds
    all roots (eigenvalues mu).

    Parameters
    ----------
    P : sympy expression
        Characteristic polynomial P(R, mu).
    R_val : float
        Value of hyperradius R (bohr).

    Returns
    -------
    roots : ndarray
        Sorted real roots of P(R_val, mu) = 0.
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required")

    R_sym = Symbol('R')
    mu_sym = Symbol('mu')

    # Substitute R value
    P_at_R = P.subs(R_sym, Rational(R_val) if isinstance(R_val, int)
                    else Float(R_val, 30))

    # Convert to polynomial in mu and find roots
    try:
        poly = Poly(P_at_R, mu_sym)
        roots_sym = poly.nroots(n=30)
        roots = np.array([complex(r) for r in roots_sym])
    except Exception:
        # Fallback: use numpy polynomial root finding
        # Extract coefficients
        coeffs = Poly(P_at_R, mu_sym).all_coeffs()
        coeffs_float = [float(c) for c in coeffs]
        roots = np.roots(coeffs_float)

    # Filter real roots and sort
    real_roots = np.sort(np.real(roots[np.abs(np.imag(roots)) < 1e-10]))
    return real_roots


def verify_against_numerical(
    Z: float = 2.0,
    l_max: int = 0,
    n_basis: int = 1,
    symmetry: str = 'singlet',
    R_values: Optional[List[float]] = None,
    n_quad: int = 200,
    tol: float = 1e-10,
) -> Tuple[Any, Dict[str, Any], List[Dict[str, Any]]]:
    """Compute P(R,mu) and verify roots match numerical eigenvalues.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_basis : int
        Basis functions per channel.
    symmetry : str
        'singlet' or 'triplet'.
    R_values : list of float, optional
        R-points for verification. Default includes points beyond
        Taylor divergence radius.
    n_quad : int
        Quadrature points.
    tol : float
        Rationalization tolerance.

    Returns
    -------
    P : sympy expression
        Characteristic polynomial.
    metadata : dict
        Polynomial structure data.
    verification : list of dict
        Verification results at each R-point.
    """
    if R_values is None:
        R_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]

    P, metadata = characteristic_polynomial(
        Z=Z, l_max=l_max, n_basis=n_basis,
        symmetry=symmetry, n_quad=n_quad, tol=tol,
    )

    solver = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad,
    )
    dim = metadata['dim']

    verification = []
    for R_val in R_values:
        # Numerical eigenvalues
        evals_num, _ = solver.solve(R_val, n_channels=dim)

        # Roots from characteristic polynomial
        roots_P = solve_characteristic_polynomial(P, R_val)

        # Match eigenvalues
        max_diff = 0.0
        diffs = []
        for i in range(min(len(evals_num), len(roots_P))):
            d = abs(evals_num[i] - roots_P[i])
            diffs.append(d)
            max_diff = max(max_diff, d)

        verification.append({
            'R': R_val,
            'mu_numerical': evals_num.tolist(),
            'mu_from_P': roots_P.tolist(),
            'diffs': diffs,
            'max_diff': max_diff,
        })

    return P, metadata, verification


def explicit_mu_lmax0(
    Z: float = 2.0,
    n_basis: int = 1,
    n_quad: int = 200,
) -> Any:
    """For l_max=0, n_basis=1: mu(R) is an explicit rational function of R.

    The 1x1 matrix H = [casimir] + R*[V_C_00] gives:
        mu(R) = casimir + R * V_C_00

    This is LINEAR in R — the simplest possible case.

    Returns
    -------
    mu_expr : sympy expression
        mu as a function of R.
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy is required")

    R_sym = Symbol('R')

    H0, V_C = extract_pencil_matrices(Z=Z, l_max=0, n_basis=1, n_quad=n_quad)

    casimir_val = Rational(int(round(H0[0, 0])))
    vc_val = nsimplify(V_C[0, 0], rational=True, tolerance=1e-10)

    mu_expr = casimir_val + R_sym * vc_val
    return mu_expr


def polynomial_structure_report(
    Z: float = 2.0,
    l_max_values: Optional[List[int]] = None,
    n_basis: int = 1,
    n_quad: int = 200,
) -> List[Dict[str, Any]]:
    """Generate a structural report of P(R,mu) at each l_max.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max_values : list of int, optional
        l_max values to analyze. Default [0, 1, 2].
    n_basis : int
        Basis functions per channel.
    n_quad : int
        Quadrature points.

    Returns
    -------
    reports : list of dict
        One report per l_max with polynomial structure data.
    """
    if l_max_values is None:
        l_max_values = [0, 1, 2]

    reports = []
    for lm in l_max_values:
        P, meta = characteristic_polynomial(
            Z=Z, l_max=lm, n_basis=n_basis, n_quad=n_quad,
        )

        report = {
            'l_max': lm,
            'dim': meta['dim'],
            'degree_R': meta['degree_R'],
            'degree_mu': meta['degree_mu'],
            'n_terms': meta['n_terms'],
            'H0_rational': meta['H0_rational'],
            'V_C_rational': meta['V_C_rational'],
            'P_rational': meta['P_rational'],
            'P_expression': str(P),
        }
        reports.append(report)

    return reports


# ─────────────────────────────────────────────────────────────────────
# Algebraic V_eff matrix elements (Track P2)
# ─────────────────────────────────────────────────────────────────────

def _stieltjes_ordinary_laguerre(
    N: int, a: float
) -> np.ndarray:
    """Stieltjes integral matrix for ordinary Laguerre polynomials (alpha=0).

    J0[i,j] = int_0^inf L_i(x) L_j(x) e^{-x} / (x + a) dx

    where a > 0.  Uses the same Miller backward recurrence strategy as
    prolate_spheroidal_lattice._stieltjes_matrix, specialised to alpha=0.

    The single-polynomial Stieltjes transforms are:
        S_n(a) = int_0^inf L_n(x) e^{-x} / (x + a) dx

    Base case:  S_0(a) = e^a E_1(a).
    Recurrence: (n+1) S_{n+1} = (2n+1+a) S_n - n S_{n-1}
    (derived from multiplying x L_j three-term recurrence by 1/(x+a)).

    Parameters
    ----------
    N : int
        Matrix dimension.
    a : float
        Shift parameter (must be > 0).

    Returns
    -------
    J0 : ndarray of shape (N, N)
        Symmetric Stieltjes integral matrix.
    """
    from scipy.special import exp1
    import math

    if a <= 0:
        raise ValueError(f"Shift parameter a must be > 0, got {a}")

    # h_n = Gamma(n+1)/n! = 1 for ordinary Laguerre (alpha=0)
    # (orthonormality: int L_i L_j e^{-x} dx = delta_{ij})

    # Step 1: S_0(a) = e^a * E_1(a)
    if a > 500:
        S0_val = 0.0
        term = 1.0 / a
        for k in range(40):
            S0_val += term
            term *= -(k + 1) / a
            if abs(term) < 1e-16 * abs(S0_val):
                break
    else:
        S0_val = float(np.exp(a) * exp1(a))

    # S_1 = (1 + a) S_0 - 1  (from n=0 inhomogeneous equation, h_0 = 1)
    S1_val = (1.0 + a) * S0_val - 1.0

    # Step 2: Miller backward recurrence for S_n
    # For small a, more padding is needed for stable convergence.
    M_pad = N + max(80, int(40 / max(a, 0.01)))
    f = np.zeros(M_pad + 1)
    f[M_pad] = 0.0
    f[M_pad - 1] = 1.0
    for n in range(M_pad - 1, 0, -1):
        # (n+1) S_{n+1} = (2n+1+a) S_n - n S_{n-1}
        # => S_{n-1} = [(2n+1+a) f[n] - (n+1) f[n+1]] / n
        f[n - 1] = ((2 * n + 1 + a) * f[n] - (n + 1) * f[n + 1]) / max(n, 1)

    S_arr = np.zeros(N)
    S_arr[0] = S0_val
    if N > 1 and abs(f[1]) > 1e-300:
        scale = S1_val / f[1]
        for n in range(1, N):
            S_arr[n] = f[n] * scale

    # Step 3: Hybrid forward/backward recurrence for J0[i,j]
    # Inhomogeneous equation (from substituting x L_j = -(j+1)L_{j+1} + (2j+1)L_j - j L_{j-1}):
    # (j+1) J0[i,j+1] = (2j+1+a) J0[i,j] - j J0[i,j-1] - delta_{ij}
    J0 = np.zeros((N, N))
    M_J = N + max(80, int(40 / max(a, 0.01)))

    for i in range(N):
        row_fwd = np.zeros(max(N, 2))
        row_fwd[0] = S_arr[i]
        if N > 1 or M_J > 1:
            # j=0 step: (1) J0[i,1] = (1+a) J0[i,0] - delta_{i,0}
            row_fwd[1] = (1.0 + a) * S_arr[i] - (1.0 if i == 0 else 0.0)

        for j in range(1, min(i + 1, N - 1)):
            # (j+1) J0[i,j+1] = (2j+1+a) J0[i,j] - j J0[i,j-1] - delta_{ij}
            row_fwd[j + 1] = ((2 * j + 1 + a) * row_fwd[j]
                              - j * row_fwd[j - 1]
                              - (1.0 if i == j else 0.0)) / (j + 1)

        for j in range(min(i + 2, N)):
            J0[i, j] = row_fwd[j]

        # Backward for j > i
        if i + 2 < N:
            bridge = row_fwd[i + 1] if i + 1 < N else 0.0

            g = np.zeros(M_J + 1)
            g[M_J] = 0.0
            g[M_J - 1] = 1.0
            for j in range(M_J - 1, i + 1, -1):
                g[j - 1] = ((2 * j + 1 + a) * g[j]
                            - (j + 1) * g[j + 1]) / max(j, 1)

            if abs(g[i + 1]) > 1e-300:
                sc = bridge / g[i + 1]
                for j in range(i + 2, N):
                    J0[i, j] = g[j] * sc

    J0 = 0.5 * (J0 + J0.T)
    return J0


def _stieltjes_x2_weighted(
    N: int, a: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute weighted Stieltjes integrals for ordinary Laguerre polynomials.

    Returns two matrices:

    J1[i,j] = int_0^inf x^2 L_i(x) L_j(x) e^{-x} / (x + a) dx
    J2[i,j] = int_0^inf x^2 L_i(x) L_j(x) e^{-x} / (x + a)^2 dx

    These are computed algebraically from the moment matrices M0, M1 and the
    basic Stieltjes matrix J0, using the identities:

        x^2 / (x+a) = (x+a) - 2a + a^2/(x+a)
        => J1 = M1 + a*M0 - 2a*M0 + a^2*J0 = M1 - a*M0 + a^2*J0

        x^2 / (x+a)^2 = 1 - 2a/(x+a) + a^2/(x+a)^2
        => J2 = M0 - 2a*J0 + a^2 * dJ0

    where dJ0 = -dJ0/da. The derivative dJ0/da is computed by differentiating
    the Stieltjes recurrence with respect to a.

    Parameters
    ----------
    N : int
        Matrix dimension.
    a : float
        Shift parameter (must be > 0).

    Returns
    -------
    J1 : ndarray of shape (N, N)
        x^2-weighted Stieltjes matrix (for 1/R integrals).
    J2 : ndarray of shape (N, N)
        x^2-weighted squared-Stieltjes matrix (for 1/R^2 integrals).
    """
    from geovac.hyperspherical_radial import _laguerre_moment_matrices

    M0, M1, _ = _laguerre_moment_matrices(N)
    J0 = _stieltjes_ordinary_laguerre(N, a)

    # J1 = int x^2 L_i L_j e^{-x} / (x+a) dx
    # Identity: x^2/(x+a) = x + a - 2a + a^2/(x+a) = x - a + a^2/(x+a)
    # => J1 = M1 - a*M0 + a^2*J0
    J1 = M1 - a * M0 + a**2 * J0

    # J2 = int x^2 L_i L_j e^{-x} / (x+a)^2 dx
    # Identity: x^2/(x+a)^2 = 1 - 2a/(x+a) + a^2/(x+a)^2
    # J_sq = int L_i L_j e^{-x} / (x+a)^2 dx  (analytical recurrence)
    J_sq = _stieltjes_x2_squared_analytical(N, a)

    J2 = M0 - 2 * a * J0 + a**2 * J_sq

    return J1, J2


def _stieltjes_x2_squared_analytical(
    N: int, a: float
) -> np.ndarray:
    """Compute int x^2 L_i L_j e^{-x} / (x+a)^2 dx analytically.

    Uses the identity:
        x^2/(x+a)^2 = 1 - 2a/(x+a) + a^2/(x+a)^2

    and computes int L_i L_j e^{-x}/(x+a)^2 dx via the backward
    recurrence for the squared Stieltjes transform.

    The single-polynomial squared Stieltjes transforms:
        T_n(a) = int_0^inf L_n(x) e^{-x} / (x+a)^2 dx = -d/da S_n(a)

    satisfy the same three-term recurrence as S_n with modified
    inhomogeneous terms (obtained by differentiating the S_n recurrence).

    Parameters
    ----------
    N : int
        Matrix dimension.
    a : float
        Shift parameter (must be > 0).

    Returns
    -------
    J_sq : ndarray of shape (N, N)
        int L_i L_j e^{-x} / (x+a)^2 dx
    """
    from scipy.special import exp1

    if a <= 0:
        raise ValueError(f"Shift parameter a must be > 0, got {a}")

    # T_n(a) = -d/da S_n(a) = int L_n(x) e^{-x} / (x+a)^2 dx
    # From S_0 = e^a E_1(a), we get:
    # T_0 = -d/da [e^a E_1(a)] = -[e^a E_1(a) + e^a(-1/a e^{-a})]
    #      = -e^a E_1(a) + 1/a = -S_0 + 1/a
    # Recurrence: differentiate (n+1) S_{n+1} = (2n+1+a) S_n - n S_{n-1}
    # => (n+1) T_{n+1} = (2n+1+a) T_n - n T_{n-1} + S_n

    if a > 500:
        S0_val = 0.0
        term = 1.0 / a
        for k in range(40):
            S0_val += term
            term *= -(k + 1) / a
            if abs(term) < 1e-16 * abs(S0_val):
                break
    else:
        S0_val = float(np.exp(a) * exp1(a))

    T0_val = -S0_val + 1.0 / a

    # Need S_n values for the inhomogeneous term
    # First compute S_n via the same backward recurrence
    S1_val = (1.0 + a) * S0_val - 1.0
    M_pad = N + 80
    f = np.zeros(M_pad + 1)
    f[M_pad] = 0.0
    f[M_pad - 1] = 1.0
    for n in range(M_pad - 1, 0, -1):
        f[n - 1] = ((2 * n + 1 + a) * f[n] - (n + 1) * f[n + 1]) / max(n, 1)

    S_arr = np.zeros(N)
    S_arr[0] = S0_val
    if N > 1 and abs(f[1]) > 1e-300:
        scale = S1_val / f[1]
        for n in range(1, N):
            S_arr[n] = f[n] * scale

    # T_1 from the differentiated n=0 equation:
    # Differentiate S_1 = (1+a) S_0 - 1 w.r.t. a:
    #   dS_1/da = S_0 + (1+a) dS_0/da
    #   -T_1 = S_0 - (1+a) T_0
    #   T_1 = (1+a) T_0 - S_0
    T1_val = (1.0 + a) * T0_val - S_arr[0]

    # Forward recurrence for T_n (inhomogeneous):
    # Differentiate (n+1) S_{n+1} = (2n+1+a) S_n - n S_{n-1} w.r.t. a:
    #   (n+1) dS_{n+1}/da = S_n + (2n+1+a) dS_n/da - n dS_{n-1}/da
    #   -(n+1) T_{n+1} = S_n - (2n+1+a) T_n + n T_{n-1}
    #   (n+1) T_{n+1} = (2n+1+a) T_n - n T_{n-1} - S_n
    T_arr = np.zeros(N)
    T_arr[0] = T0_val
    if N > 1:
        T_arr[1] = T1_val
    for n in range(1, N - 1):
        T_arr[n + 1] = ((2 * n + 1 + a) * T_arr[n]
                        - n * T_arr[n - 1]
                        - S_arr[n]) / (n + 1)

    # Build J_sq[i,j] = int L_i L_j e^{-x}/(x+a)^2 dx
    #
    # Bilinear recurrence (derived from (x+a) L_j three-term recurrence):
    #   (j+1) J_sq[i,j+1] = (2j+1+a) J_sq[i,j] - j J_sq[i,j-1] - J0[i,j]
    # Seeds: J_sq[i,0] = T_arr[i], J_sq[i,1] = (1+a)*T_arr[i] - J0[i,0]

    J0 = _stieltjes_ordinary_laguerre(N, a)

    # Build J_sq[i,j] via forward recurrence in j for each i.
    # The driving term J0[i,j] is nonzero for all j (unlike delta_{ij}
    # in the J0 recurrence), so the hybrid forward/backward approach
    # is not applicable.  Pure forward recurrence is used instead.
    # For N <= 30 (typical spectral basis size), forward recurrence
    # is sufficiently stable.
    J_sq = np.zeros((N, N))

    for i in range(N):
        row = np.zeros(max(N, 2))
        row[0] = T_arr[i]
        if N > 1:
            row[1] = (1.0 + a) * T_arr[i] - J0[i, 0]

        for j in range(1, N - 1):
            row[j + 1] = ((2 * j + 1 + a) * row[j]
                          - j * row[j - 1]
                          - J0[i, j]) / (j + 1)

        J_sq[i, :N] = row[:N]

    J_sq = 0.5 * (J_sq + J_sq.T)
    return J_sq


def algebraic_veff_matrix_lmax0(
    Z: float = 2.0,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_quad: int = 200,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Compute the V_eff matrix algebraically for l_max=0 (Level 3 He).

    At l_max=0, mu(R) = c0 + c1*R (explicit linear function) where:
      c0 = H0[0,0] = 0 (SO(6) Casimir for l=0, k=0 singlet)
      c1 = V_C[0,0] = nuclear + V_ee coupling

    So V_eff(R) = (c0 + 15/8)/R^2 + c1/R.

    The matrix element in the Laguerre spectral basis is:
      V_nm = (1/(8*alpha^3)) * int_0^inf x^2 L_n(x) L_m(x) e^{-x} * V_eff(R) dx

    where x = 2*alpha*(R - R_min), R = R_min + x/(2*alpha) = (a+x)/(2*alpha),
    and a = 2*alpha*R_min.

    Substituting V_eff = A/R^2 + B/R = A*4*alpha^2/(a+x)^2 + B*2*alpha/(a+x):
      V_nm = (1/(8*alpha^3)) * [A*4*alpha^2 * J2[n,m] + B*2*alpha * J1[n,m]]
           = A/(2*alpha) * J2[n,m] + B/(4*alpha^2) * J1[n,m]

    where:
      J1[n,m] = int x^2 L_n L_m e^{-x} / (x+a) dx
      J2[n,m] = int x^2 L_n L_m e^{-x} / (x+a)^2 dx
      A = c0 + 15/8
      B = c1

    Parameters
    ----------
    Z : float
        Nuclear charge (default 2.0 for He).
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary (must be > 0).
    n_quad : int
        Quadrature points for angular solver.

    Returns
    -------
    V_mat : ndarray of shape (n_basis, n_basis)
        Algebraic V_eff matrix.
    info : dict
        Diagnostic information including c0, c1, A, B, a.
    """
    from geovac.hyperspherical_radial import _laguerre_moment_matrices

    H0, V_C = extract_pencil_matrices(Z=Z, l_max=0, n_basis=1, n_quad=n_quad)

    c0 = float(H0[0, 0])  # Should be 0 for l=0 singlet
    c1 = float(V_C[0, 0])  # Nuclear + V_ee coupling

    A = c0 + 15.0 / 8.0  # Coefficient of 1/R^2
    B = c1                # Coefficient of 1/R

    a = 2.0 * alpha * R_min  # Shift parameter

    N = n_basis

    # Compute weighted Stieltjes matrices
    M0, M1, _ = _laguerre_moment_matrices(N)
    J0 = _stieltjes_ordinary_laguerre(N, a)

    # J1 = int x^2 L_n L_m e^{-x} / (x+a) dx
    # = M1 - a*M0 + a^2*J0
    J1 = M1 - a * M0 + a**2 * J0

    # J2 = int x^2 L_n L_m e^{-x} / (x+a)^2 dx
    # = M0 - 2a*J0 + a^2 * J_sq
    J_sq = _stieltjes_x2_squared_analytical(N, a)
    J2 = M0 - 2 * a * J0 + a**2 * J_sq

    # V_nm = A/(2*alpha) * J2 + B/(4*alpha^2) * J1
    V_mat = A / (2.0 * alpha) * J2 + B / (4.0 * alpha**2) * J1

    info = {
        'c0': c0,
        'c1': c1,
        'A': A,
        'B': B,
        'a': a,
        'J0_cond': np.linalg.cond(J0),
        'J1_trace': np.trace(J1),
        'J2_trace': np.trace(J2),
    }

    return V_mat, info


def algebraic_veff_matrix_quadrature(
    Z: float = 2.0,
    l_max: int = 0,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_quad_angular: int = 200,
) -> np.ndarray:
    """Compute V_eff matrix via Gauss-Laguerre quadrature (reference method).

    This is the existing quadrature approach, extracted as a standalone function
    for comparison with the algebraic method.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave for angular solver.
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary.
    n_quad_angular : int
        Quadrature points for angular solver.

    Returns
    -------
    V_mat : ndarray of shape (n_basis, n_basis)
        V_eff matrix from quadrature.
    """
    from scipy.special import roots_laguerre, eval_laguerre
    from geovac.algebraic_angular import AlgebraicAngularSolver

    N = n_basis
    two_alpha = 2.0 * alpha

    # Build angular solver to get V_eff(R)
    solver = AlgebraicAngularSolver(
        Z=Z, n_basis=1, l_max=l_max,
        symmetry='singlet', n_quad=n_quad_angular,
    )

    def V_eff_func(R_arr: np.ndarray) -> np.ndarray:
        V = np.zeros_like(R_arr)
        for k, R in enumerate(R_arr):
            if R > 1e-10:
                evals, _ = solver.solve(R, n_channels=l_max + 1)
                mu_val = evals[0]  # Ground state adiabatic potential
                V[k] = mu_val / R**2 + 15.0 / (8.0 * R**2)
            else:
                V[k] = 1e10  # Regularise
        return V

    # Gauss-Laguerre quadrature
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_min + x_quad / two_alpha

    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    inv_8a3 = 1.0 / (8.0 * alpha**3)
    x2_w = w_quad * x_quad**2
    V_q = V_eff_func(R_q)
    x2_wVL = (x2_w * V_q)[np.newaxis, :] * L_vals
    V_mat = inv_8a3 * (x2_wVL @ L_vals.T)

    return V_mat


def algebraic_he_energy_lmax0(
    Z: float = 2.0,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_quad: int = 200,
    matrix_method_SK: str = 'algebraic',
) -> Tuple[float, Dict[str, Any]]:
    """Compute He ground state energy using algebraic V_eff at l_max=0.

    Combines algebraic S, K (from Track H) with algebraic V_eff (this track)
    to form the full Hamiltonian. The generalized eigenvalue problem
    H c = E S c is solved, with H = K + V_eff.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary.
    n_quad : int
        Quadrature points for angular solver.
    matrix_method_SK : str
        'algebraic' or 'quadrature' for S and K matrices.

    Returns
    -------
    E : float
        Ground state energy (Ha).
    info : dict
        Diagnostic information.
    """
    from geovac.hyperspherical_radial import _build_laguerre_SK_algebraic

    # S and K from Track H (algebraic)
    if matrix_method_SK == 'algebraic':
        S, K = _build_laguerre_SK_algebraic(n_basis, alpha)
    else:
        # Fall back to quadrature for comparison
        from geovac.hyperspherical_radial import _build_laguerre_matrices_dirichlet
        from geovac.algebraic_angular import AlgebraicAngularSolver
        solver = AlgebraicAngularSolver(
            Z=Z, n_basis=1, l_max=0, symmetry='singlet', n_quad=n_quad,
        )

        def dummy_V(R):
            return np.zeros_like(R)

        S, K, _, _, _, _, _, _ = _build_laguerre_matrices_dirichlet(
            dummy_V, n_basis, alpha, R_min, matrix_method='quadrature'
        )

    # V_eff from algebraic computation
    V_mat, veff_info = algebraic_veff_matrix_lmax0(
        Z=Z, n_basis=n_basis, alpha=alpha, R_min=R_min, n_quad=n_quad,
    )

    H = K + V_mat
    evals, evecs = eigh(H, S)

    E_ground = evals[0]

    info = {
        'E_ground': E_ground,
        'veff_info': veff_info,
        'S_cond': np.linalg.cond(S),
        'H_dim': n_basis,
    }

    return E_ground, info


def analyze_discriminant_lmax1(
    Z: float = 2.0,
    n_quad: int = 200,
) -> Dict[str, Any]:
    """Analyze the discriminant of the l_max=1 characteristic polynomial.

    At l_max=1, n_basis=1, the characteristic polynomial is quadratic in mu:
        P(R, mu) = mu^2 + b(R)*mu + c(R) = 0

    where b(R) and c(R) are polynomials in R (degree 1 and 2 respectively).
    The discriminant Delta(R) = b^2 - 4c determines whether mu(R) can be
    expressed via radicals.

    This function characterizes Delta(R) to determine:
    1. Its degree as a polynomial in R
    2. Whether the resulting integral has known closed forms
    3. The nature of the radical sqrt(Delta) as a function of R

    Returns
    -------
    analysis : dict
        Discriminant analysis results.
    """
    if not HAS_SYMPY:
        raise ImportError("SymPy required for discriminant analysis")

    R_sym = Symbol('R')
    mu_sym = Symbol('mu')

    P, meta = characteristic_polynomial(
        Z=Z, l_max=1, n_basis=1, n_quad=n_quad,
    )

    # Extract coefficients of mu^2, mu, 1
    poly_mu = sympy.Poly(P, mu_sym)
    coeffs = poly_mu.all_coeffs()  # [coeff of mu^2, coeff of mu, constant]

    a2 = coeffs[0]  # Should be 1 (monic)
    a1 = coeffs[1]  # Linear in R
    a0 = coeffs[2]  # Quadratic in R

    # Discriminant
    disc = sympy.expand(a1**2 - 4 * a2 * a0)

    # Analyze degree in R
    try:
        disc_poly = sympy.Poly(disc, R_sym)
        disc_degree = disc_poly.degree()
    except Exception:
        disc_degree = None

    # Check if discriminant is a perfect square (would give rational mu)
    # For degree 2, disc(R) = pR^2 + qR + r, which is generally not a perfect square

    # Evaluate numerically at several R values
    disc_values = {}
    for R_val in [0.5, 1.0, 2.0, 5.0, 10.0]:
        disc_val = float(disc.subs(R_sym, R_val))
        disc_values[R_val] = disc_val

    # mu(R) = [-a1 +/- sqrt(disc)] / (2*a2)
    # V_eff(R) = mu(R)/R^2 + 15/(8R^2)
    # The integral involves sqrt(polynomial in R) * Laguerre functions,
    # which is generally an elliptic integral for degree >= 3.

    analysis = {
        'P_expression': str(P),
        'a2': str(a2),
        'a1': str(a1),
        'a0': str(a0),
        'discriminant': str(disc),
        'disc_degree_R': disc_degree,
        'disc_values': disc_values,
        'disc_sign_at_R1': disc_values.get(1.0, None),
        'integral_type': (
            'rational' if disc_degree is not None and disc_degree <= 0
            else 'square_root_linear' if disc_degree == 1
            else 'square_root_quadratic' if disc_degree == 2
            else 'unknown'
        ),
        'note': (
            'For disc_degree=2, V_eff involves sqrt(aR^2+bR+c)/R^2, '
            'the Laguerre integral is generally non-elementary (elliptic type).'
        ),
    }

    return analysis


def detect_recurrence_lmax1(
    Z: float = 2.0,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_quad_angular: int = 200,
    max_order: int = 8,
) -> Dict[str, Any]:
    """Detect linear recurrence in V_nm matrix elements for l_max=1 (Tier 2c).

    Computes V_nm via quadrature for n,m = 0..n_basis-1 and searches for a
    finite-order linear recurrence in n at fixed m, with coefficients that
    are polynomial in m.

    Uses the Hankel matrix null-vector method: if V_{n,m} satisfies a
    recurrence of order p in n, then the (p+1)x(p+1) Hankel matrix
    H[i,j] = V_{i+j, m} has a null vector.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Number of basis functions (determines matrix size).
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary.
    n_quad_angular : int
        Quadrature points for angular solver.
    max_order : int
        Maximum recurrence order to search for.

    Returns
    -------
    result : dict
        Recurrence detection results.
    """
    V_quad = algebraic_veff_matrix_quadrature(
        Z=Z, l_max=1, n_basis=n_basis, alpha=alpha,
        R_min=R_min, n_quad_angular=n_quad_angular,
    )

    results_by_m = {}
    for m in range(min(5, n_basis)):
        col = V_quad[:, m]

        best_order = None
        best_sv_ratio = 0.0

        for p in range(2, min(max_order + 1, n_basis // 2)):
            # Build Hankel matrix
            n_rows = n_basis - 2 * p
            if n_rows < 2:
                break
            H = np.zeros((n_rows, p + 1))
            for i in range(n_rows):
                for j in range(p + 1):
                    H[i, j] = col[i + j]

            # Check for null vector via SVD
            sv = np.linalg.svd(H, compute_uv=False)
            if len(sv) > 0 and sv[0] > 0:
                ratio = sv[-1] / sv[0]
                if ratio < 1e-10:
                    best_order = p
                    best_sv_ratio = ratio
                    break

        results_by_m[m] = {
            'order_found': best_order,
            'sv_ratio': best_sv_ratio,
        }

    has_recurrence = any(r['order_found'] is not None for r in results_by_m.values())

    return {
        'n_basis': n_basis,
        'results_by_m': results_by_m,
        'has_recurrence': has_recurrence,
        'note': (
            'If recurrence found, V_eff can be computed from boundary values '
            'without pointwise quadrature.'
        ),
    }
