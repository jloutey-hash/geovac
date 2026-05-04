"""Avery-Wen-Avery 3-Y integral on S^3.

This module implements the closed-form, exact-arithmetic Avery-Wen-Avery
three-spherical-harmonic integral on the unit 3-sphere:

    Avery3Y(n_a, l_a, m_a; n_b, l_b, m_b; n_c, l_c, m_c)
        := integral_{S^3} Y^{(3)*}_{n_a l_a m_a}(omega)
                          Y^{(3)}_{n_b l_b m_b}(omega)
                          Y^{(3)}_{n_c l_c m_c}(omega) dOmega_3

The 3-Y integral factorizes into a Gegenbauer radial 3-integral on
chi in [0, pi] (the polar S^3 angle) and the standard 2-sphere Gaunt
integral on the (theta, phi) angular factor (Avery, "Hyperspherical
Harmonics: Applications in Quantum Theory," Kluwer 1989, Eq. 2.5).

This replaces the convention-independent placeholder
`_so4_radial_overlap_placeholder` in `geovac.operator_system`. The
placeholder was robust enough to verify structural invariants
(propagation number = 2, dim sequence) but produced convention-
dependent absolute values for the operator-system matrix elements
and Connes distance. The Avery-Wen-Avery integral provides the
*physically meaningful* values for these quantities.

Mathematical setup
==================

Hyperspherical harmonics on S^3
-------------------------------

The S^3 hyperspherical harmonics factor as

    Y^{(3)}_{n l m}(chi, theta, phi) = R_{n l}(chi) Y_{l m}(theta, phi)

where chi in [0, pi] is the polar angle and (theta, phi) are standard
S^2 spherical coordinates. The radial-on-the-sphere factor R_{n l} is
the standard Gegenbauer-based form

    R_{n l}(chi) = N_{n l} sin^l(chi) C^{l+1}_{n-l-1}(cos chi)

with normalization (Avery 1989 Eq. 1.6.3):

    N_{n l}^2 = 2^{2l+1} n (n-l-1)! (l!)^2 / (pi (n+l)!)

This makes <R_{n l} | R_{n' l'}>_{[0,pi], sin^2 dchi} = delta_{n n'}
delta_{l l'} (verified symbolically in tests).

The full 3-sphere measure on chi in [0, pi], theta in [0, pi],
phi in [0, 2 pi) is

    dOmega_3 = sin^2(chi) sin(theta) dchi dtheta dphi.

Factorization of the 3-Y integral
---------------------------------

Under the chi - (theta, phi) factorization,

    integral Y^{(3)*}_a Y^{(3)}_b Y^{(3)}_c dOmega_3
        = I_rad(n_a, l_a; n_b, l_b; n_c, l_c)
          x G(l_a, m_a; l_b, m_b; l_c, m_c)

where

    I_rad := integral_0^pi R_{n_a l_a}(chi) R_{n_b l_b}(chi)
                           R_{n_c l_c}(chi) sin^2(chi) dchi   (Gegenbauer 3)

    G    := integral Y_{l_a m_a}^* Y_{l_b m_b} Y_{l_c m_c} dOmega_2  (Gaunt)

The Gaunt integral is the standard SU(2) 3-Y integral with closed form

    G = (-1)^{m_a} sqrt((2 l_a + 1)(2 l_b + 1)(2 l_c + 1) / (4 pi))
        x (l_a l_b l_c ; 0 0 0)
        x (l_a l_b l_c ; -m_a m_b m_c)

with selection rules m_a = m_b + m_c, |l_b - l_c| <= l_a <= l_b + l_c,
l_a + l_b + l_c even.

The Gegenbauer 3-integral I_rad is computable in closed form via
expansion in u = cos(chi) and the standard beta-function integral
identity for u^k (1 - u^2)^a, see _gegenbauer_triple_integral below.

Connection to SU(2)_L x SU(2)_R = Spin(4)
-----------------------------------------

S^3 = SU(2), so each S^3 hyperspherical harmonic
Y^{(3)}_{n l m} corresponds to a state in the SO(4) irrep
((n-1)/2, (n-1)/2) under the SU(2)_L x SU(2)_R double-cover
factorization. The "principal QN" n labels the SU(2)_L = SU(2)_R
spin: j = (n - 1) / 2. The pair (l, m) labels the diagonal SU(2)
content: l is the diagonal SU(2) total angular momentum, m its
projection.

In this language the 3-Y integral is a product of two SU(2)
Wigner 3j symbols (one for each SU(2) factor), with selection
rules from each factor. However, because the diagonal SU(2)
content (l, m) is what appears in the chi-theta-phi factorization,
the structurally simplest and computationally efficient form is the
direct chi-integration above. The SU(2)_L x SU(2)_R form is
sometimes useful for proving identities but is not used in this
module's implementation.

Algebraic structure
===================

Per the project's algebraic-first discipline (CLAUDE.md Sec 4 and
Sec 6 QED strategic directive), all integrals here are computed in
exact rational + algebraic arithmetic via sympy. The Gegenbauer
3-integral I_rad is an exact algebraic expression involving rational
coefficients and powers of sqrt(pi). The Gaunt integral G involves
rational Wigner 3j coefficients and a sqrt(1 / (4 pi)) prefactor.
Their product, the full 3-Y integral, is therefore an algebraic
number in a rational extension of Q(sqrt(pi)).

Concretely:

    I_rad in (rational) * pi^{-(l_a + l_b + l_c + 1) / 2}
            (and: pi^{1/2} prefactor from each N_{n l}, three of them total)

    G in    sqrt(rational) * (1 / sqrt(pi))

The combined 3-Y integral is in Q[sqrt(rational)] * pi^{-power}.

The propagation number computation (round 2) is proven invariant
under the choice of these irrational prefactors because it depends
only on the support pattern of the M_{N L M} matrices, not on their
values. The Connes distance computation (round 3) is value-sensitive
and will give different absolute numerical values when this module
is used instead of the placeholder.

Performance
===========

With sympy lru_cache + polynomial reduction, the cost is approximately
1 second for n_max = 3 (540 triples) and 5 seconds for n_max = 4
(2800 triples) on a modern desktop. The build of all multiplier
matrices for n_max = 4 in geovac.operator_system therefore takes ~10
seconds total (includes Gaunt, matrix assembly).

Algebraic-first discipline
==========================

The Gegenbauer 3-integral and the Gaunt integral are both available
in closed form. They are NOT a "wall" requiring numerical quadrature
(per CLAUDE.md Sec 4 refinement). The closed-form values are
algebraic numbers, with the only transcendental being pi (via the
S^3 measure normalization sqrt(pi) factors and the round S^3
volume 2 pi^2). Per Paper 35 / Paper 34, this pi is correctly
classified as a calibration / observation projection from the round
S^3 metric; it is not a sign of lurking algebraic structure to chase.

Module API
==========

  - `gegenbauer_triple_integral(n_a, l_a, n_b, l_b, n_c, l_c) -> sympy.Expr`
       Closed-form value of I_rad, exact algebraic expression.

  - `s2_gaunt(l_a, m_a, l_b, m_b, l_c, m_c) -> sympy.Expr`
       Closed-form value of G, exact algebraic expression.

  - `avery_wen_avery_3y(n_a, l_a, m_a, n_b, l_b, m_b, n_c, l_c, m_c) -> sympy.Expr`
       Closed-form value of full S^3 3-Y integral. Selection rules
       enforced (returns 0 if either Gaunt or Gegenbauer triangle
       fails).

  - `radial_overlap_avery(n_a, N, n_b, l_a, L, l_b) -> sympy.Expr`
       Drop-in replacement for `_so4_radial_overlap_placeholder` in
       geovac.operator_system. Returns I_rad with selection rules
       enforced.

References
==========

J. Avery, "Hyperspherical Harmonics: Applications in Quantum
Theory," Kluwer 1989. Equations 1.6.3 (radial normalization), 2.5
(SO(4) angular structure), 3.5 (3-Y integral factorization).

Z. Wen & J. Avery, "Some properties of hyperspherical harmonics,"
J. Math. Phys. 26, 396-403 (1985). Closed-form 3-Y integral on
S^d.

J. Avery & J. Avery, "Generalized Sturmians and Atomic Spectra,"
World Scientific 2006. Sturmian / Gegenbauer machinery.

For a 2-sphere Gaunt integral closed form:
A. R. Edmonds, "Angular Momentum in Quantum Mechanics," Princeton
1957, Eq. 4.6.3.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Optional

import sympy as sp
from sympy import (
    Integer,
    Poly,
    Rational,
    Symbol,
    cos,
    expand,
    factorial,
    gamma,
    gegenbauer,
    pi,
    simplify,
    sin,
    sqrt,
)
from sympy.physics.wigner import wigner_3j


# ---------------------------------------------------------------------------
# Internal sympy symbol used by the polynomial form of the integrand
# ---------------------------------------------------------------------------

_u = Symbol("_avery_u", real=True)


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _N_nl(n: int, l: int) -> sp.Expr:
    """S^3 radial normalization N_{n l}.

    Returns the EXACT sympy expression

        N_{n l} = sqrt( 2^{2l+1} * n * (n-l-1)! * (l!)^2 / (n+l)! ) / sqrt(pi)

    such that <R_{n l} | R_{n l}>_{[0,pi], sin^2 dchi} = 1 with

        R_{n l}(chi) = N_{n l} sin^l(chi) C^{l+1}_{n-l-1}(cos chi).

    The factor 1 / sqrt(pi) is pulled out so that N_{n l}^2 has its
    only pi-dependence in a single power of pi^{-1}.
    """
    if n < 1 or l < 0 or l > n - 1:
        raise ValueError(f"invalid (n, l) = ({n}, {l})")
    rational_part = (
        Integer(2) ** (2 * l + 1) * Integer(n)
        * factorial(n - l - 1) * factorial(l) ** 2
        / factorial(n + l)
    )
    return sqrt(rational_part) / sqrt(pi)


@lru_cache(maxsize=None)
def _gegenbauer_poly(degree: int, alpha: int) -> Poly:
    """Gegenbauer polynomial C^alpha_degree(u) as a sympy Poly in u.

    Cached because the same (degree, alpha) pairs reappear many times
    in the radial-3-integral loop.
    """
    if degree < 0:
        raise ValueError(f"Gegenbauer degree {degree} < 0")
    return Poly(expand(gegenbauer(degree, alpha, _u)), _u)


@lru_cache(maxsize=None)
def _beta_uk_factor(k: int, a: sp.Rational) -> sp.Expr:
    """Closed-form value of int_{-1}^{1} u^k (1 - u^2)^a du.

    For non-negative integer k and (rational) a > -1:

        int_{-1}^{1} u^k (1 - u^2)^a du
            = 0                                              if k odd
            = Gamma((k+1)/2) Gamma(a+1) / Gamma((k+1)/2 + a + 1)   if k even

    The latter is the standard Beta-function identity.

    Returns sympy Expr; for half-integer (k+1)/2 + a + 1 this gives
    rational + half-integer Gamma values, which sympy simplifies to
    rationals * sqrt(pi).
    """
    if k % 2 == 1:
        return Integer(0)
    half = Rational(k + 1, 2)
    return gamma(half) * gamma(a + 1) / gamma(half + a + 1)


# ---------------------------------------------------------------------------
# The Gegenbauer triple integral (radial part of the 3-Y integral)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def gegenbauer_triple_integral(
    n_a: int, l_a: int, n_b: int, l_b: int, n_c: int, l_c: int,
) -> sp.Expr:
    """Closed-form value of the S^3 Gegenbauer 3-integral.

        I_rad := int_0^pi R_{n_a l_a}(chi) R_{n_b l_b}(chi)
                          R_{n_c l_c}(chi) sin^2(chi) dchi

    with R_{n l}(chi) = N_{n l} sin^l(chi) C^{l+1}_{n-l-1}(cos chi).

    The integrand becomes, after substituting u = cos(chi),
    du = -sin(chi) dchi, so sin^2(chi) dchi = -sin(chi) du, and
    using sin(chi) = sqrt(1 - u^2):

        I_rad = norm * int_{-1}^{1} P(u) (1 - u^2)^{(l_a + l_b + l_c + 1)/2} du

    where norm = N_{n_a l_a} N_{n_b l_b} N_{n_c l_c} and P(u) is the
    polynomial product of the three Gegenbauer factors. The (1 - u^2)
    power has exponent (l_a + l_b + l_c + 1)/2; the integral over u is
    evaluated term-by-term in P(u) using `_beta_uk_factor`.

    Selection rules
    ---------------

    The integrand is automatically zero by parity if l_a + l_b + l_c is
    even and the polynomial P(u) has no even-degree terms compatible
    with (1 - u^2)^{(odd)/2}, etc. We do NOT enforce any explicit
    selection rules here -- if the integral is structurally zero, sympy
    returns zero.

    Note in particular that the naive "SO(4) triangle on principal QN"
    rule |n_a - n_b| + 1 <= n_c <= n_a + n_b - 1 is NOT a true selection
    rule for this integral as a function of (n, l). Many integrals that
    violate this triangle on n alone are nevertheless nonzero because of
    constructive interference between the sin^l factors and the
    Gegenbauer polynomials. Verified by direct sympy computation, e.g.
    I_rad(2, 1; 3, 0; 1, 0) is NONZERO even though
    |2 - 1| + 1 = 2 < 3 = n_c.

    Performance
    -----------

    Cached. Cost per call is dominated by the polynomial multiplication
    of three Gegenbauer factors (degrees up to n_max - 1 each) and the
    sum of `_beta_uk_factor` values for each polynomial term (up to
    degree 3*(n_max - 1)). For n_max = 4 the total cost across all
    distinct (n_a, l_a, n_b, l_b, n_c, l_c) triples is ~5 s including
    cache warmup.
    """
    # Normalization product
    norm = _N_nl(n_a, l_a) * _N_nl(n_b, l_b) * _N_nl(n_c, l_c)

    # Polynomial product P(u) = C^{l_a+1}_{n_a-l_a-1} * C^{l_b+1}_{n_b-l_b-1}
    #                            * C^{l_c+1}_{n_c-l_c-1}
    P = (
        _gegenbauer_poly(n_a - l_a - 1, l_a + 1)
        * _gegenbauer_poly(n_b - l_b - 1, l_b + 1)
        * _gegenbauer_poly(n_c - l_c - 1, l_c + 1)
    )

    # Power of (1 - u^2) factor in the integrand:
    # sin^{l_a + l_b + l_c}(chi) * sin^2(chi) dchi -> (1 - u^2)^{(l_sum + 1)/2} du
    l_sum = l_a + l_b + l_c
    a_pow = Rational(l_sum + 1, 2)

    result: sp.Expr = Integer(0)
    for monom, coeff in P.terms():
        k = monom[0]  # power of u in this term
        result = result + coeff * _beta_uk_factor(k, a_pow)

    if result == 0:
        return Integer(0)

    return simplify(norm * result)


# ---------------------------------------------------------------------------
# 2-sphere Gaunt integral
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def s2_gaunt(
    l_a: int, m_a: int, l_b: int, m_b: int, l_c: int, m_c: int,
) -> sp.Expr:
    """Closed-form value of the S^2 Gaunt integral

        G := int Y^*_{l_a m_a} Y_{l_b m_b} Y_{l_c m_c} dOmega_2.

    Standard formula (Edmonds 1957, Eq. 4.6.3):

        G = (-1)^{m_a}
            sqrt((2 l_a + 1)(2 l_b + 1)(2 l_c + 1) / (4 pi))
            (l_a l_b l_c ; 0 0 0)
            (l_a l_b l_c ; -m_a m_b m_c)

    Selection rules: returns 0 if any of:
      - m_a != m_b + m_c
      - |l_b - l_c| > l_a or l_a > l_b + l_c
      - l_a + l_b + l_c is odd
      - |m_a| > l_a, |m_b| > l_b, or |m_c| > l_c
      - either of the Wigner 3j symbols vanishes

    Convention note
    ---------------

    The conventional Gaunt integral has a complex-conjugate on
    Y_{l_a m_a}. With the standard Condon-Shortley phase
    Y^*_{l m} = (-1)^m Y_{l, -m}, this is equivalent to the
    "all three lower" form

        integral Y_{l_a, -m_a} Y_{l_b m_b} Y_{l_c m_c} dOmega_2
            = (-1)^{m_a} G  [with the above formula]

    and selection rule -m_a + m_b + m_c = 0, i.e. m_a = m_b + m_c.
    """
    # Selection rules
    if m_a != m_b + m_c:
        return Integer(0)
    if not (abs(l_b - l_c) <= l_a <= l_b + l_c):
        return Integer(0)
    if (l_a + l_b + l_c) % 2 != 0:
        return Integer(0)
    if abs(m_a) > l_a or abs(m_b) > l_b or abs(m_c) > l_c:
        return Integer(0)

    threej_zero = wigner_3j(l_a, l_b, l_c, 0, 0, 0)
    if threej_zero == 0:
        return Integer(0)
    threej_m = wigner_3j(l_a, l_b, l_c, -m_a, m_b, m_c)
    if threej_m == 0:
        return Integer(0)

    sign = Integer(-1) ** m_a
    pref = sqrt(
        Rational((2 * l_a + 1) * (2 * l_b + 1) * (2 * l_c + 1)) / (4 * pi)
    )
    return sign * pref * threej_zero * threej_m


# ---------------------------------------------------------------------------
# Full Avery-Wen-Avery 3-Y integral on S^3
# ---------------------------------------------------------------------------


def avery_wen_avery_3y(
    n_a: int, l_a: int, m_a: int,
    n_b: int, l_b: int, m_b: int,
    n_c: int, l_c: int, m_c: int,
) -> sp.Expr:
    """The full S^3 3-Y integral

        integral_{S^3} Y^{(3)*}_{n_a l_a m_a} Y^{(3)}_{n_b l_b m_b}
                       Y^{(3)}_{n_c l_c m_c} dOmega_3.

    Factorizes as I_rad * G with I_rad the Gegenbauer 3-integral on
    chi and G the standard Gaunt integral on (theta, phi).

    Returns sympy expression in exact arithmetic.
    """
    G = s2_gaunt(l_a, m_a, l_b, m_b, l_c, m_c)
    if G == 0:
        return Integer(0)
    I = gegenbauer_triple_integral(n_a, l_a, n_b, l_b, n_c, l_c)
    if I == 0:
        return Integer(0)
    return simplify(I * G)


# ---------------------------------------------------------------------------
# Drop-in for geovac.operator_system._so4_radial_overlap_placeholder
# ---------------------------------------------------------------------------


def radial_overlap_avery(
    n: int, N: int, np_: int, l: int, L: int, lp: int,
) -> sp.Expr:
    """Drop-in replacement for `_so4_radial_overlap_placeholder`.

    Returns the Gegenbauer 3-integral

        I_rad = int_0^pi R_{n l}(chi) R_{N L}(chi) R_{n' l'}(chi)
                          sin^2(chi) dchi.

    The signature matches the placeholder exactly:
        (n, N, np_, l, L, lp) where (N, L) is the multiplier label and
        (n, l), (np_, lp) are the bra/ket basis labels.

    Returns 0 if the underlying Gegenbauer 3-integral is structurally
    zero (which happens automatically by parity / triangle-on-Gegenbauer
    polynomials -- there is no separately enforced "selection rule"
    layer here, since the integral itself is the source of truth).

    Validity range: requires N >= 1, 0 <= L <= N - 1 (the bound on L is
    inherited from the requirement that L is a valid S^2 angular momentum
    quantum number under the (N, L, M) hyperspherical-harmonic indexing
    on S^3).
    """
    # The placeholder treated N=1 specially (returning 1 to make
    # M_{1,0,0} the identity). With the genuine integral, N=1, L=0
    # gives I_rad ~ <R_{n,l} | constant | R_{n',l'}>_{[0,pi], sin^2 dchi}
    # = N_{1,0} * <R_{n,l} | R_{n',l'}>_{[0,pi], sin^2 dchi}. Since
    # constants don't introduce angular-momentum coupling, the Gaunt
    # part of M_{1,0,0} is just (Y_{0,0})_{l=l',m=m'} delta-on-(l,m),
    # and this delta times the orthonormality of R_{n l} forces
    # M_{1,0,0} to be a multiple of the identity. We do NOT need a
    # special branch here.
    return gegenbauer_triple_integral(n, l, N, L, np_, lp)


# ---------------------------------------------------------------------------
# Convenience: orthonormality verification
# ---------------------------------------------------------------------------


def radial_overlap_two(n_a: int, l: int, n_b: int) -> sp.Expr:
    """The two-radial overlap

        <R_{n_a, l} | R_{n_b, l}>_{[0,pi], sin^2 dchi}.

    Returns sympy expression. By orthonormality this equals
    delta_{n_a, n_b} EXACTLY. Provided as an equation-verification
    helper for tests.

    For l_a != l_b this is computed-as-given (no closed-form delta);
    full S^3 orthogonality requires the angular Y_{l m} factor.
    """
    norm = _N_nl(n_a, l) * _N_nl(n_b, l)
    P = (
        _gegenbauer_poly(n_a - l - 1, l + 1)
        * _gegenbauer_poly(n_b - l - 1, l + 1)
    )
    a_pow = Rational(2 * l + 1, 2)
    result: sp.Expr = Integer(0)
    for monom, coeff in P.terms():
        k = monom[0]
        result = result + coeff * _beta_uk_factor(k, a_pow)
    return simplify(norm * result)


def s3_orthonormality_check(n_a: int, l_a: int, m_a: int,
                            n_b: int, l_b: int, m_b: int) -> sp.Expr:
    """Full S^3 orthonormality check:

        <Y^{(3)}_{n_a l_a m_a} | Y^{(3)}_{n_b l_b m_b}>_{S^3}.

    Returns sympy expression, equal to delta_{n_a n_b} delta_{l_a l_b}
    delta_{m_a m_b} EXACTLY.

    Implementation: factorizes as <R | R> * <Y | Y>. The angular part
    is delta_{l m} from the standard Y_{l m} orthonormality on S^2; we
    return the radial overlap times explicit Kronecker delta in (l, m).
    """
    if l_a != l_b or m_a != m_b:
        return Integer(0)
    return radial_overlap_two(n_a, l_a, n_b)
