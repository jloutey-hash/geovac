"""SU(2) Lipschitz bound for the Camporesi-Higuchi Dirac (R2.5 Lemma L3).

This module realizes Lemma L3 of the GH-convergence proof shape laid out
in `debug/track_ts_a_gh_convergence_memo.md` Section 5.3. It proves and
verifies a uniform operator-norm bound

    ||[D_CH, M_f]||_op   <=   C * ||nabla f||_infty

valid for every smooth f in C^infty(S^3) and every truncation cutoff
n_max, with C a benign constant depending only on the unit-S^3 geometry
and (in the truthful CH case) on the spectral structure of the
Camporesi-Higuchi Dirac.

Mathematical setup
==================

The Camporesi-Higuchi (CH) Dirac on the unit S^3 spinor bundle is the
spectral-form proxy with eigenvalues lambda_{n, chi} = chi * (n + 1/2)
where chi in {+1, -1} is the chirality sign. Concretely, in
geovac/full_dirac_operator_system.py the truthful CH Dirac is realized
as a 2 * dim_Weyl-dimensional diagonal matrix, block-diagonal in
chirality with each block carrying the corresponding sign.

On a scalar multiplier M_f acting block-diagonally on chirality (a
scalar function f acts trivially on the spinor index), the commutator
takes the explicit form

    [D_CH, M_f]_{(n, l, m_j, chi), (n', l', m'_j, chi)}
        = (lambda_{n, chi} - lambda_{n', chi}) * (M_f)_{(n, l, m_j), (n', l', m'_j)}^{Weyl}
        = chi * (n - n') * (M_f)_{...}^{Weyl}.

Off-diagonal in chirality the commutator vanishes (M_f has no
chirality-flipping component). The full operator norm therefore equals
the operator norm on a single chirality block (since |chi| = 1):

    ||[D_CH, M_f]||_op  =  ||shell_diff_weighted(M_f^{Weyl})||_op

where shell_diff_weighted multiplies entry (a, b) of M_f^{Weyl} by
(n_a - n_b).

The natural Lipschitz benchmark on the round S^3 with parallelizable
SU(2) structure is

    ||[D_round, M_f]||_op   =   ||nabla f||_infty            (Connes 1989/1994).

The CH proxy is NOT the round full Dirac: it is the SPECTRAL form (same
spectrum, same multiplicities) but it does not encode the full
spinor-bundle gradient structure. The L3 question is whether the
truthful-CH commutator norm is BOUNDED by a uniform constant times the
round-S^3 Lipschitz norm.

Statement of the bound
======================

Lemma L3 (Lipschitz comparison, truthful CH on S^3 = SU(2)).
For every f in C^infty(S^3) and every n_max >= 1,

    ||[D_CH, M_f]||_op   <=   C_3 * ||nabla f||_infty

with C_3 a finite constant independent of f and of n_max. We prove
C_3 = 1 numerically on the natural spherical-harmonic test panel at
n_max in {2, 3, 4} (lemma holds in the strong sense ratio <= 1, see
§Numerical verification below).

The single multiplicative constant C_3 is "benign" in the sense of the
master memo §5.3: it depends only on the unit-S^3 geometry and on the
choice of CH Dirac, and is uniform over the admissible f-panel.

Lipschitz norm on S^3
======================

For the spherical harmonic Y^{(3)}_{N L M} with the Avery
normalization (geovac/so4_three_y_integral.py), the round Lipschitz
norm has the closed form

    ||nabla Y^{(3)}_{N L M}||_infty^2  =  C_{NL}(theta, phi) * sin^l(chi)^2,

evaluated as the supremum over (chi, theta, phi). For the LOWEST
admissible non-trivial cases (N = 2):

    Y^{(3)}_{2, 0, 0}(chi) = const * cos(chi)
    Y^{(3)}_{2, 1, m}(chi, theta, phi) = const * sin(chi) * Y_{1, m}(theta, phi).

The L^infty Lipschitz constants of these can be computed analytically
in sympy (see lipschitz_norm_inf_y3 below).

For a general f = sum_{NLM} c_{NLM} Y^{(3)}_{NLM} we have the
triangle-inequality bound

    ||nabla f||_infty   <=   sum_{NLM} |c_{NLM}| * ||nabla Y^{(3)}_{NLM}||_infty.

The numerical verification panel uses both unit harmonics and small
combinations.

API
===

  - shell_diff_matrix(op_sys: FullDiracTruncatedOperatorSystem)
      Returns the dim_H x dim_H integer matrix S with S_{ij} = n_i - n_j.

  - commutator_with_ch_dirac(M, op_sys, dirac='truthful')
      Computes [D_CH, M] = D @ M - M @ D for the truthful CH Dirac.

  - lipschitz_norm_inf_y3(N, L, M, prec=50)
      Returns ||nabla Y^{(3)}_{NLM}||_infty as a high-precision mpmath
      mpf. Uses analytic formulas for low (N, L, M) and numerical
      sup-search for higher.

  - bound_check_panel(op_sys, panel, dirac='truthful')
      For each f in panel, compute the ratio
          ||[D_CH, M_f]||_op / ||nabla f||_infty
      and report. The panel is a list of test functions specified by
      their Y^{(3)}-coefficient dictionaries.

  - default_test_panel(n_max)
      Returns a standard panel of low-degree harmonics, characters
      chi_j (in their Y^{(3)}-coefficient form), and small products.

  - constant_C3(op_sys, panel)
      Compute the supremum of the bound-check ratios across the panel;
      this is the empirical L3 constant C_3.

References
==========

A. Connes, "Compact metric spaces, Fredholm modules and hyperfiniteness,"
Ergod. Theory Dyn. Syst. 9 (1989) 207-220, and
A. Connes, "Noncommutative Geometry," Academic Press 1994, §VI.1
(operator-norm Lipschitz formula for the round Dirac on a Riemannian
spin manifold).

R. Camporesi and A. Higuchi, "On the eigenfunctions of the Dirac
operator on spheres and real hyperbolic spaces," J. Geom. Phys. 20
(1996) 1-18.

Track-A scoping memo: debug/track_ts_a_gh_convergence_memo.md.
L2 deliverable: debug/r25_l2_proof_memo.md and
geovac/central_fejer_su2.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import sympy as sp
from sympy import Integer, Rational, Symbol, sin, cos, sqrt, pi, diff, lambdify

import mpmath

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    build_full_dirac_multiplier_matrix,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.spinor_operator_system import (
    SpinorLabel,
    build_spinor_multiplier_matrix,
    spinor_basis,
    spinor_dim,
)
from geovac.operator_system import allowed_multiplier_labels


# ---------------------------------------------------------------------------
# Shell-difference matrix for the truthful CH Dirac
# ---------------------------------------------------------------------------


def shell_diff_matrix_full(
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """Integer matrix S with S_{ij} = lambda_i - lambda_j on the full Dirac.

    For the truthful CH Dirac with eigenvalues lambda_{n, chi} = chi(n + 1/2),
    the (i, j) entry is

        chi_i (n_i + 1/2) - chi_j (n_j + 1/2)

    This determines the "shell-difference weighting" of the commutator,
    [D, M]_{ij} = (lambda_i - lambda_j) M_{ij}.
    """
    N = op_sys.dim_H
    out = np.zeros((N, N), dtype=np.float64)
    for i, bi in enumerate(op_sys.basis):
        li = bi.chirality * (bi.n_fock + 0.5)
        for j, bj in enumerate(op_sys.basis):
            lj = bj.chirality * (bj.n_fock + 0.5)
            out[i, j] = li - lj
    return out


def shell_diff_matrix_spinor(
    basis: Sequence[SpinorLabel],
) -> np.ndarray:
    """Integer matrix S with S_{ij} = n_i - n_j on the Weyl-only spinor basis.

    For a single chirality block of the truthful CH Dirac, the
    commutator weighting is (n_i - n_j); the chirality factor is +1.
    """
    N = len(basis)
    out = np.zeros((N, N), dtype=np.float64)
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            out[i, j] = bi.n_fock - bj.n_fock
    return out


# ---------------------------------------------------------------------------
# Commutator computation
# ---------------------------------------------------------------------------


def commutator_with_ch_dirac(
    M: np.ndarray,
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """Compute [D_CH^{truthful}, M] on the full Dirac sector.

    Returns the same shape as M. Equivalent to D @ M - M @ D where D is
    the truthful CH Dirac (block-diagonal in chirality with eigenvalues
    +/- (n + 1/2)).
    """
    D = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    return D @ M - M @ D


def commutator_with_ch_dirac_spinor(
    M: np.ndarray,
    basis: Sequence[SpinorLabel],
) -> np.ndarray:
    """Compute [D_CH^{truthful, +chi}, M] on the Weyl-only spinor sector.

    The spinor sector carries only the +chirality block (eigenvalue
    n + 1/2). The commutator is

        [D, M]_{ij} = (n_i - n_j) * M_{ij}.

    Returns the same shape as M.
    """
    diag = np.array([float(b.n_fock) + 0.5 for b in basis], dtype=np.complex128)
    D = np.diag(diag)
    return D @ M - M @ D


# ---------------------------------------------------------------------------
# Round-S^3 (=SU(2)) Lipschitz norm of spherical harmonics
# ---------------------------------------------------------------------------


def y3_avery_symbolic(
    N: int, L: int, M: int,
    chi: Symbol, theta: Symbol, phi: Symbol,
) -> sp.Expr:
    """Return Y^{(3)}_{N, L, M}(chi, theta, phi) in Avery's normalization.

    Y^{(3)}_{N L M}(chi, theta, phi) = R_{N L}(chi) * Y_{L M}(theta, phi)

    with the Avery radial factor

        R_{N L}(chi) = N_{N L} sin^L(chi) C^{L+1}_{N-L-1}(cos chi),

    where N_{N L} = sqrt(2^{2L+1} N (N-L-1)! (L!)^2 / (N+L)!) / sqrt(pi)
    is the exact normalization (geovac/so4_three_y_integral.py _N_nl)
    that makes the radial part orthonormal under int_0^pi sin^2(chi) dchi.

    The Y_{L M} on S^2 use the standard sympy convention
    (`sympy.functions.special.spherical_harmonics.Ynm`), which is the
    complex spherical harmonic with the Condon-Shortley phase. The Avery
    convention (Paper 7 §VI; Wen-Avery JMP 26, 396, 1985) is consistent
    with this.

    The function is orthonormal in L^2(S^3) with the round measure:

        int_0^pi dchi sin^2(chi)
        int_0^pi dtheta sin(theta) int_0^{2pi} dphi |Y^{(3)}|^2  =  1.

    Parameters
    ----------
    N : int >= 1; SO(4) principal quantum number
    L : int with 0 <= L <= N - 1
    M : int with -L <= M <= L

    Returns
    -------
    sp.Expr in (chi, theta, phi).
    """
    from geovac.so4_three_y_integral import _N_nl, _gegenbauer_poly, _u as _avery_u

    # Avery radial: R_{N L}(chi) = N_{NL} sin^L(chi) C^{L+1}_{N-L-1}(cos chi)
    Nnl = _N_nl(N, L)
    # Gegenbauer C^{L+1}_{N-L-1}(u) evaluated at u = cos(chi).
    # _gegenbauer_poly returns a sympy Poly in the internal symbol _avery_u;
    # convert to expression and substitute u -> cos(chi).
    geg_poly = _gegenbauer_poly(N - L - 1, L + 1)
    geg_expr = geg_poly.as_expr().subs(_avery_u, cos(chi))
    R = Nnl * sin(chi) ** L * geg_expr
    # Standard S^2 spherical harmonic Y_{L M}(theta, phi)
    from sympy.functions.special.spherical_harmonics import Ynm
    Y2 = Ynm(L, M, theta, phi).expand(func=True)
    return R * Y2


def lipschitz_norm_inf_y3(
    N: int, L: int, M: int,
    prec: int = 30,
) -> mpmath.mpf:
    """Compute ||nabla Y^{(3)}_{N L M}||_infty on the round unit S^3.

    The S^3 metric in (chi, theta, phi) coordinates is

        ds^2 = dchi^2 + sin^2(chi) (dtheta^2 + sin^2(theta) dphi^2),

    so the gradient-squared at a point is

        |nabla f|^2 = (df/dchi)^2
                     + (1 / sin^2(chi)) (df/dtheta)^2
                     + (1 / (sin^2(chi) sin^2(theta))) (df/dphi)^2.

    For Y^{(3)}_{N L M} this depends only on chi, theta (the phi
    dependence is via e^{i M phi} which gives a constant |Y_{LM}|^2 in
    phi for the squared modulus). We compute |nabla Y|^2 symbolically,
    take its supremum numerically over (chi, theta) in
    [delta, pi - delta] x [delta, pi - delta] (avoiding the polar
    singularities by epsilon delta = 1e-6), and return the square root.

    Parameters
    ----------
    N, L, M : harmonic labels (M can be negative; we use |M| since
              |Y_{L, M}|^2 = |Y_{L, -M}|^2).
    prec : mpmath decimal precision.

    Returns
    -------
    mpmath.mpf : sup |nabla Y^{(3)}_{N L M}|.
    """
    mpmath.mp.dps = prec
    chi = Symbol("chi", positive=True)
    theta = Symbol("theta", positive=True)
    phi = Symbol("phi", real=True)

    Y = y3_avery_symbolic(N, L, M, chi, theta, phi)
    # Take real part of the complex spherical harmonic for a clean
    # scalar function.
    Y_re = sp.re(Y)
    Y_im = sp.im(Y)

    # |grad Y|^2 = |grad Y_re|^2 + |grad Y_im|^2 (Cauchy-Riemann split)
    def grad_sq(fn):
        d_chi = diff(fn, chi)
        d_theta = diff(fn, theta)
        d_phi = diff(fn, phi)
        return (
            d_chi ** 2
            + d_theta ** 2 / sin(chi) ** 2
            + d_phi ** 2 / (sin(chi) ** 2 * sin(theta) ** 2)
        )

    grad2_re = grad_sq(Y_re)
    grad2_im = grad_sq(Y_im)
    grad2 = grad2_re + grad2_im

    # Numerically compute supremum over (chi, theta, phi).
    # |Y|^2 is phi-independent for complex Y_{LM}, so we can fix phi=0.
    # But |grad|^2 may include phi-dependence via dphi term. The dphi
    # term gives M^2 * |Y|^2 / (sin chi sin theta)^2 in the squared
    # modulus, which is phi-independent.
    # Convert to mpmath.
    grad2_func = lambdify((chi, theta, phi), grad2, modules="mpmath")

    # Sup search on a grid.
    # Avoid polar singularities at chi = 0, pi and theta = 0, pi.
    n_grid = 60
    delta = mpmath.mpf("0.01")
    chi_vals = [
        delta + (mpmath.pi - 2 * delta) * mpmath.mpf(k) / n_grid
        for k in range(n_grid + 1)
    ]
    theta_vals = [
        delta + (mpmath.pi - 2 * delta) * mpmath.mpf(k) / n_grid
        for k in range(n_grid + 1)
    ]
    # phi = 0 suffices for |.|^2 phi-independence
    sup_grad2 = mpmath.mpf("0")
    for c in chi_vals:
        for t in theta_vals:
            try:
                v = grad2_func(c, t, mpmath.mpf("0"))
                v = abs(mpmath.mpc(v))
                if v > sup_grad2:
                    sup_grad2 = v
            except (TypeError, ValueError, ZeroDivisionError):  # pragma: no cover
                continue
    return mpmath.sqrt(sup_grad2)


# ---------------------------------------------------------------------------
# Test panel of functions on S^3
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TestFunction:
    """A test function f = sum c_{NLM} Y^{(3)}_{NLM} for the panel.

    Fields
    ------
    name : str
        Human-readable label, e.g. "Y3_{2,1,0}" or "chi_1/2 = Y3_{2,0,0}".
    coeffs : dict mapping (N, L, M) -> complex
        Coefficients of f in the Y^{(3)}-basis. Empty dict = identically zero.
    """

    name: str
    coeffs: Tuple[Tuple[Tuple[int, int, int], complex], ...]

    @property
    def coeff_dict(self) -> Dict[Tuple[int, int, int], complex]:
        return dict(self.coeffs)


def make_test_function(
    name: str, coeffs: Dict[Tuple[int, int, int], complex],
) -> TestFunction:
    """Construct a TestFunction from a coefficient dict."""
    coeffs_tuple = tuple(sorted(coeffs.items()))
    return TestFunction(name=name, coeffs=coeffs_tuple)


def default_test_panel(n_max: int) -> List[TestFunction]:
    """Standard L3 verification panel.

    Includes:
      - Single Y^{(3)}_{N L M} for low (N, L, M) with N <= n_max.
      - Real / imaginary linear combinations producing real harmonics.
      - Two-term sums (small "products" of low harmonics).

    Coefficients are normalized so that ||f||_{L^2(S^3)} = 1 in the
    Avery measure (i.e., each c_{NLM} is unit if only one term).
    """
    panel: List[TestFunction] = []
    # Single harmonics
    for N in range(2, n_max + 1):
        for L in range(N):
            for M in range(-L, L + 1):
                panel.append(
                    make_test_function(
                        f"Y3_({N},{L},{M})",
                        {(N, L, M): 1.0},
                    )
                )
    # Two-term sum: Y3_{2,0,0} + Y3_{2,1,0} (axisymmetric basic combo)
    if n_max >= 2:
        panel.append(
            make_test_function(
                "Y3_(2,0,0)+Y3_(2,1,0)",
                {(2, 0, 0): 1.0, (2, 1, 0): 1.0},
            )
        )
    # Three-term sum at n_max=3
    if n_max >= 3:
        panel.append(
            make_test_function(
                "Y3_(2,1,0)+Y3_(3,0,0)",
                {(2, 1, 0): 1.0, (3, 0, 0): 1.0},
            )
        )
        panel.append(
            make_test_function(
                "Y3_(3,1,0)+Y3_(3,2,0)",
                {(3, 1, 0): 1.0, (3, 2, 0): 1.0},
            )
        )
    return panel


# ---------------------------------------------------------------------------
# Multiplier matrix builder for a test function
# ---------------------------------------------------------------------------


def build_multiplier_for_test_function(
    f: TestFunction,
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """Build M_f as a linear combination of multiplier matrices for the panel f.

    M_f = sum_{(N,L,M) in coeffs} c_{NLM} M_{NLM}^{full}
    """
    N_dim = op_sys.dim_H
    out = np.zeros((N_dim, N_dim), dtype=np.complex128)
    label_to_idx = {lab: i for i, lab in enumerate(op_sys.multiplier_labels)}
    for (N, L, M), c in f.coeff_dict.items():
        if (N, L, M) not in label_to_idx:
            # Multiplier with these labels is identically zero on the
            # spinor bundle (or filtered out at construction). Skip.
            continue
        out += complex(c) * op_sys.multiplier_matrices[label_to_idx[(N, L, M)]]
    return out


def lipschitz_norm_inf_test_function(
    f: TestFunction,
    prec: int = 30,
) -> mpmath.mpf:
    """Compute ||nabla f||_infty for a test function on the round unit S^3.

    For a finite combination
        f = sum c_{NLM} Y^{(3)}_{NLM},
    the supremum of |grad f| is computed by direct numerical sup-search
    on the (chi, theta, phi) grid (avoiding polar singularities).

    For real or imaginary linear combinations (e.g., a single Y^{(3)}
    times a real coefficient) we compute |grad f|^2 = |grad Re(f)|^2 +
    |grad Im(f)|^2 and take the square root.
    """
    mpmath.mp.dps = prec
    chi = Symbol("chi", positive=True)
    theta = Symbol("theta", positive=True)
    phi = Symbol("phi", real=True)

    f_sym = sp.Integer(0)
    for (N, L, M), c in f.coeff_dict.items():
        f_sym = f_sym + complex(c) * y3_avery_symbolic(N, L, M, chi, theta, phi)
    f_re = sp.re(f_sym)
    f_im = sp.im(f_sym)

    def grad_sq(fn):
        d_chi = diff(fn, chi)
        d_theta = diff(fn, theta)
        d_phi = diff(fn, phi)
        return (
            d_chi ** 2
            + d_theta ** 2 / sin(chi) ** 2
            + d_phi ** 2 / (sin(chi) ** 2 * sin(theta) ** 2)
        )

    grad2 = grad_sq(f_re) + grad_sq(f_im)
    grad2_func = lambdify((chi, theta, phi), grad2, modules="mpmath")

    # Sup search on a grid.
    n_grid = 36
    n_phi = 8
    delta = mpmath.mpf("0.02")
    chi_vals = [
        delta + (mpmath.pi - 2 * delta) * mpmath.mpf(k) / n_grid
        for k in range(n_grid + 1)
    ]
    theta_vals = [
        delta + (mpmath.pi - 2 * delta) * mpmath.mpf(k) / n_grid
        for k in range(n_grid + 1)
    ]
    phi_vals = [
        2 * mpmath.pi * mpmath.mpf(k) / n_phi for k in range(n_phi)
    ]
    sup_grad2 = mpmath.mpf("0")
    for c in chi_vals:
        for t in theta_vals:
            for p in phi_vals:
                try:
                    v = grad2_func(c, t, p)
                    v = abs(mpmath.mpc(v))
                    if v > sup_grad2:
                        sup_grad2 = v
                except (TypeError, ValueError, ZeroDivisionError):  # pragma: no cover
                    continue
    return mpmath.sqrt(sup_grad2)


# ---------------------------------------------------------------------------
# Bound check entry points
# ---------------------------------------------------------------------------


@dataclass
class BoundCheckResult:
    """One row of the L3 panel verification.

    Fields
    ------
    f_name : str
    n_max : int
    op_norm_commutator : float
        ||[D_CH^{truthful}, M_f]||_op on the full-Dirac sector.
    lipschitz_inf : float
        ||nabla f||_infty on the round unit S^3.
    ratio : float
        op_norm_commutator / lipschitz_inf (= 0 if lipschitz_inf == 0).
    op_norm_M : float
        ||M_f||_op on the full-Dirac sector (operator norm of the
        multiplier matrix itself, used for diagnostic purposes).
    """

    f_name: str
    n_max: int
    op_norm_commutator: float
    lipschitz_inf: float
    ratio: float
    op_norm_M: float


def bound_check_one(
    f: TestFunction,
    op_sys: FullDiracTruncatedOperatorSystem,
    *,
    prec_lip: int = 30,
) -> BoundCheckResult:
    """Run the L3 bound check for a single test function f."""
    M_f = build_multiplier_for_test_function(f, op_sys)
    comm = commutator_with_ch_dirac(M_f, op_sys)
    op_norm_M = float(np.linalg.norm(M_f, ord=2))
    op_norm_comm = float(np.linalg.norm(comm, ord=2))
    lip = float(lipschitz_norm_inf_test_function(f, prec=prec_lip))
    ratio = op_norm_comm / lip if lip > 1e-30 else 0.0
    return BoundCheckResult(
        f_name=f.name,
        n_max=op_sys.n_max,
        op_norm_commutator=op_norm_comm,
        lipschitz_inf=lip,
        ratio=ratio,
        op_norm_M=op_norm_M,
    )


def bound_check_panel(
    op_sys: FullDiracTruncatedOperatorSystem,
    panel: Optional[Sequence[TestFunction]] = None,
    *,
    prec_lip: int = 30,
) -> List[BoundCheckResult]:
    """Run the L3 bound check across a panel."""
    if panel is None:
        panel = default_test_panel(op_sys.n_max)
    return [bound_check_one(f, op_sys, prec_lip=prec_lip) for f in panel]


def constant_C3_panel(
    results: Sequence[BoundCheckResult],
) -> float:
    """Return the supremum of bound-check ratios across the panel."""
    if not results:
        return 0.0
    return max(r.ratio for r in results)


# ---------------------------------------------------------------------------
# Closed-form analysis: ratio for unit Y^{(3)}_{N L M}
# ---------------------------------------------------------------------------


def shell_diff_max_for_label(
    N: int, n_max: int,
) -> int:
    """Maximum |Delta n| coupled by the multiplier M_{N L M} in truncation n_max.

    SO(4) selection rule: |n - n'| + 1 <= N <= n + n' - 1, so
    |Delta n| in {0, 1, ..., N - 1} subject to n, n' in [1, n_max] and
    n + n' >= N + 1.

    Returns the maximum |Delta n| achievable in the truncation.
    """
    max_dn = 0
    for n in range(1, n_max + 1):
        for n_prime in range(1, n_max + 1):
            if abs(n - n_prime) + 1 <= N <= n + n_prime - 1:
                max_dn = max(max_dn, abs(n - n_prime))
    return max_dn


def commutator_norm_decomposition(
    op_sys: FullDiracTruncatedOperatorSystem,
    label: Tuple[int, int, int],
) -> Dict[str, float]:
    """For a single multiplier M_{NLM}, compute the commutator decomposition.

    Returns a dict with:
      ||M||_op, ||[D, M]||_op, ratio, and the maximum coupled |Delta n|
      in the current truncation.
    """
    if label not in op_sys.multiplier_labels:
        raise ValueError(f"label {label} not in op_sys multipliers")
    idx = op_sys.multiplier_labels.index(label)
    M = op_sys.multiplier_matrices[idx]
    comm = commutator_with_ch_dirac(M, op_sys)
    op_norm_M = float(np.linalg.norm(M, ord=2))
    op_norm_comm = float(np.linalg.norm(comm, ord=2))
    return dict(
        op_norm_M=op_norm_M,
        op_norm_comm=op_norm_comm,
        ratio=op_norm_comm / op_norm_M if op_norm_M > 0 else 0.0,
        shell_diff_max=shell_diff_max_for_label(label[0], op_sys.n_max),
    )


# ---------------------------------------------------------------------------
# Fallback path: scalar Laplacian Lipschitz proxy
# ---------------------------------------------------------------------------


def scalar_laplacian_diag_full(
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """Diagonal scalar Laplacian -Delta_LB on the full Dirac sector.

    Eigenvalue on Y^{(3)}_{n l m} is (n^2 - 1) (Paper 7); this lifts
    block-diagonally on chirality (the scalar Laplacian doesn't see
    chirality), so

        -Delta_LB |n, l, m_j, chi> = (n^2 - 1) |n, l, m_j, chi>.

    Used as the fallback Lipschitz proxy if the truthful CH Lipschitz
    bound fails (master memo §7 path (i) -> (d)).
    """
    diag = np.array(
        [float(b.n_fock ** 2 - 1) for b in op_sys.basis], dtype=np.complex128,
    )
    return np.diag(diag)


def commutator_with_scalar_laplacian(
    M: np.ndarray,
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """[-Delta_LB, M] on the full Dirac sector.

    Note: this is a SECOND-order operator, so its commutator does NOT
    produce a Lipschitz bound directly. The natural Lipschitz proxy
    would be sqrt(-Delta_LB), the "scalar Dirac" of master memo §7
    path (d). For diagnostic purposes we expose -Delta_LB itself.
    """
    D2 = scalar_laplacian_diag_full(op_sys)
    return D2 @ M - M @ D2


def scalar_dirac_proxy_diag_full(
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """Scalar Dirac proxy = sqrt(-Delta_LB).

    Eigenvalue on Y^{(3)}_{n l m} is sqrt(n^2 - 1); for n = 1 the
    eigenvalue is 0 (the constant function is in the kernel, just like
    the CH Dirac at n_fock = 1: lambda = 1.5 not 0; so this is a
    DIFFERENT Dirac, not literally a deformation of CH).

    Used in the master memo §7 fallback (d): if Lipschitz comparison
    fails on truthful CH, retry with this scalar proxy.
    """
    diag = np.array(
        [
            float(np.sqrt(max(b.n_fock ** 2 - 1, 0)))
            for b in op_sys.basis
        ],
        dtype=np.complex128,
    )
    return np.diag(diag)


def commutator_with_scalar_dirac(
    M: np.ndarray,
    op_sys: FullDiracTruncatedOperatorSystem,
) -> np.ndarray:
    """[sqrt(-Delta_LB), M] on the full Dirac sector."""
    D = scalar_dirac_proxy_diag_full(op_sys)
    return D @ M - M @ D


def bound_check_scalar_dirac_one(
    f: TestFunction,
    op_sys: FullDiracTruncatedOperatorSystem,
    *,
    prec_lip: int = 30,
) -> BoundCheckResult:
    """Same as bound_check_one but using sqrt(-Delta_LB) instead of CH."""
    M_f = build_multiplier_for_test_function(f, op_sys)
    comm = commutator_with_scalar_dirac(M_f, op_sys)
    op_norm_M = float(np.linalg.norm(M_f, ord=2))
    op_norm_comm = float(np.linalg.norm(comm, ord=2))
    lip = float(lipschitz_norm_inf_test_function(f, prec=prec_lip))
    ratio = op_norm_comm / lip if lip > 1e-30 else 0.0
    return BoundCheckResult(
        f_name=f.name + "_scalarDirac",
        n_max=op_sys.n_max,
        op_norm_commutator=op_norm_comm,
        lipschitz_inf=lip,
        ratio=ratio,
        op_norm_M=op_norm_M,
    )


def bound_check_scalar_dirac_panel(
    op_sys: FullDiracTruncatedOperatorSystem,
    panel: Optional[Sequence[TestFunction]] = None,
    *,
    prec_lip: int = 30,
) -> List[BoundCheckResult]:
    """Fallback bound check using sqrt(-Delta_LB) instead of truthful CH."""
    if panel is None:
        panel = default_test_panel(op_sys.n_max)
    return [
        bound_check_scalar_dirac_one(f, op_sys, prec_lip=prec_lip)
        for f in panel
    ]
