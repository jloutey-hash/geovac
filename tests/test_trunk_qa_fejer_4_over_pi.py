"""Trunk QA — Claim A: the 4/pi asymptotic constant of Paper 38 (L2 rate).

CONTEXT
=======
The existing test ``tests/test_central_fejer_su2.py::test_asymptotic_constant_value``
is CIRCULAR: ``asymptotic_rate_constant()`` returns a hardcoded ``Rational(4)/pi``
and the test asserts equality with ``Rational(4)/pi``.  It cannot fail and proves
nothing about the actual theorem

    lim_{n->oo}  n * gamma_n / log(n)  =  4/pi                              (T1i)

where gamma_n = pi - 4 * T_n / (pi * Z_n) is the SU(2) central-Fejer mass-
concentration moment with the *independently* test-backed closed-form sum rule

    T_n = sum_{1<=k1,k2<=n, k1+k2 odd} sqrt(k1 k2) [1/(k1-k2)^2 - 1/(k1+k2)^2],
    Z_n = n(n+1)/2.

This file ATTEMPTS a genuine, falsifiable derivation of (T1i) from the gamma_n
closed form via the Stein-Weiss / Euler-Maclaurin reduction, deriving the
constant from independent series (sum_{d odd} 1/d^2 = pi^2/8 and the odd-harmonic
sum sum_{d odd<=D} 1/d ~ (1/2) log D) rather than hardcoding 4/pi.

DERIVATION (the thing under test)
=================================
gamma_n = pi - (4/pi) * (T_n / Z_n) = (4/pi) * defect_n,
    defect_n := pi^2/4 - T_n / Z_n.

The claim n*gamma_n/log n -> 4/pi is EQUIVALENT to: defect_n * n / log n -> 1,
i.e. the *defect log-coefficient* equals exactly 1.  We derive that 1 as a
difference of two structurally distinct, independently pinned sub-coefficients:

  (A) Triangle truncation of the diagonal-dominant series sum_{d odd} 1/d^2.
      tail_odd(D) := pi^2/8 - sum_{d odd<=D} 1/d^2 = 1/(2D) + O(1/D^2)  (Euler-Maclaurin).
      This produces defect log-coefficient EXACTLY 2.

  (B) The sqrt(a(a+d)) ~ a + d/2 correction, summed against 1/d^2, giving the
      odd-harmonic sum sum_{d odd<=D} 1/d ~ (1/2) log D.
      This produces defect log-coefficient EXACTLY 1.

  Net defect log-coefficient = 2 - 1 = 1, hence gamma_n ~ (4/pi) log n / n.

Each sub-coefficient is established by an independent series identity, and the
composition could have produced a DIFFERENT constant (e.g. if (B) gave 3/2 the
limit would be 2/pi, not 4/pi).  The tests below would FAIL if any of these
ingredients were wrong.  Nothing here references the package's hardcoded
``asymptotic_rate_constant()``.

VERDICT (see module docstring summary at bottom): the leading constant 4/pi is
DERIVABLE.  The two sub-coefficients (2 and 1) are pinned exactly by Euler-
Maclaurin; their net (=1) is corroborated to high precision by an independent
doubling (Richardson) estimator on the truthful closed-form gamma_n.
"""

from __future__ import annotations

import math

import mpmath
import pytest

# We deliberately use ONLY the closed-form sum rule for gamma_n / T_n, which is
# itself independently test-backed (test_central_fejer_su2.py
# ::test_T_n_sum_rule_matches_gamma_quadrature checks it against Gaussian
# quadrature of the actual kernel moment).  We do NOT import
# asymptotic_rate_constant() — that is the circular object under review.
from geovac.central_fejer_su2 import T_n_via_sum_rule, normalization_constant


PREC = 50


def _Z(n: int) -> mpmath.mpf:
    return mpmath.mpf(normalization_constant(n))


def _gamma(n: int) -> mpmath.mpf:
    """gamma_n via the (independently test-backed) closed-form sum rule."""
    mpmath.mp.dps = PREC
    T = T_n_via_sum_rule(n, prec=PREC)
    return mpmath.pi - 4 * T / (mpmath.pi * _Z(n))


def _defect(n: int) -> mpmath.mpf:
    """defect_n = pi^2/4 - T_n/Z_n  (== gamma_n * pi / 4)."""
    mpmath.mp.dps = PREC
    T = T_n_via_sum_rule(n, prec=PREC)
    return mpmath.pi ** 2 / 4 - T / _Z(n)


# ---------------------------------------------------------------------------
# Consistency: defect_n = (pi/4) * gamma_n exactly (closed-form bookkeeping).
# This is the bridge identity gamma_n = (4/pi) * defect_n.  If it failed, the
# whole derivation route would be invalid.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [3, 5, 10, 50, 100])
def test_gamma_equals_4_over_pi_times_defect(n):
    mpmath.mp.dps = PREC
    lhs = _gamma(n)
    rhs = (4 / mpmath.pi) * _defect(n)
    assert abs(lhs - rhs) < mpmath.mpf("1e-40"), (
        f"gamma_n != (4/pi) defect_n at n={n}: {lhs} vs {rhs}"
    )


# ---------------------------------------------------------------------------
# Sub-coefficient (A): triangle truncation of sum_{d odd} 1/d^2 gives EXACTLY 2.
# Independent series identity: tail_odd(D) = pi^2/8 - sum_{d odd<=D} 1/d^2
#                                          = 1/(2D) + O(1/D^2).
# This is pure Euler-Maclaurin and does not reference 4/pi at all.
# ---------------------------------------------------------------------------


def _tail_odd(D: int) -> mpmath.mpf:
    mpmath.mp.dps = PREC
    s = sum(mpmath.mpf(1) / d ** 2 for d in range(1, D + 1, 2))
    return mpmath.pi ** 2 / 8 - s


@pytest.mark.parametrize("D", [500, 1000, 2000, 4000])
def test_odd_zeta2_tail_is_half_over_D(D):
    """tail_odd(D) * 2D -> 1, i.e. tail_odd(D) ~ 1/(2D) (Euler-Maclaurin leading)."""
    val = _tail_odd(D) * 2 * D
    # Converges to 1 from below as 1 - O(1/D); at D=500 already within 2e-6.
    assert abs(val - 1) < 2e-6 * (4000 / D + 1), f"tail_odd({D})*2D = {val}"


def test_subcoeff_A_triangle_truncation_is_two():
    """The a-only triangle truncation contributes defect log-coefficient 2.

    defect_aonly(n) := pi^2/4 - (2/Z_n) sum_{a=1}^{n-1} a * S2(n-a),
        S2(D) = sum_{d odd<=D} 1/d^2.
    Claim: defect_aonly(n) * n / log n -> 2.  Pinned by a doubling estimator
    on the genuinely computed partial sums (NOT by assuming the answer).
    """
    mpmath.mp.dps = PREC

    def S2(D):
        return sum(mpmath.mpf(1) / d ** 2 for d in range(1, D + 1, 2))

    def defect_aonly(n):
        s = sum(a * S2(n - a) for a in range(1, n))
        return mpmath.pi ** 2 / 4 - 2 * s / _Z(n)

    def doubling(n):
        # n*defect = c log n + b + O(1/n); (2n f(2n) - n f(n))/log2 -> c.
        return (2 * n * defect_aonly(2 * n) - n * defect_aonly(n)) / mpmath.log(2)

    c200 = float(doubling(200))
    c400 = float(doubling(400))
    # Converging to 2 from above; deviation roughly halving.
    assert abs(c400 - 2.0) < abs(c200 - 2.0), "not converging toward 2"
    assert abs(c400 - 2.0) < 0.02, f"a-only doubling coeff = {c400}, expected ->2"


def test_subcoeff_B_sqrt_correction_is_one():
    """The sqrt(a(a+d)) ~ a + d/2 correction contributes defect log-coefficient 1.

    contrib_dhalf(n) := (2/Z_n) sum_{a=1}^{n-1} (1/2) S1(n-a),
        S1(D) = sum_{d odd<=D} 1/d ~ (1/2) log D.
    This piece ADDS to T_n/Z_n (reduces the defect), with log-coefficient 1.
    Pinned by the same independent doubling estimator.
    """
    mpmath.mp.dps = PREC

    def S1(D):
        return sum(mpmath.mpf(1) / d for d in range(1, D + 1, 2))

    def contrib(n):
        s = sum(mpmath.mpf(1) / 2 * S1(n - a) for a in range(1, n))
        return 2 * s / _Z(n)

    def doubling(n):
        return (2 * n * contrib(2 * n) - n * contrib(n)) / mpmath.log(2)

    c200 = float(doubling(200))
    c400 = float(doubling(400))
    assert abs(c400 - 1.0) < abs(c200 - 1.0), "not converging toward 1"
    assert abs(c400 - 1.0) < 0.02, f"sqrt-correction doubling coeff = {c400}, expected ->1"


def test_odd_harmonic_sum_log_coefficient_is_half():
    """sum_{d odd<=D} 1/d - (1/2) log D -> const (so the log coefficient is 1/2).

    Independent identity underpinning sub-coefficient (B).
    """
    mpmath.mp.dps = PREC

    def S1(D):
        return sum(mpmath.mpf(1) / d for d in range(1, D + 1, 2))

    vals = [S1(D) - mpmath.mpf(1) / 2 * mpmath.log(D) for D in [1000, 2000, 4000]]
    # The residual should be (essentially) constant: differences -> 0.
    assert abs(float(vals[1] - vals[0])) < 1e-6
    assert abs(float(vals[2] - vals[1])) < 1e-6
    # And the residual is the known constant ln2 + gamma_E/2 + (3/?) ... we only
    # need that the LOG coefficient is 1/2, which the constancy above proves.


# ---------------------------------------------------------------------------
# The headline: NET defect log-coefficient = 2 - 1 = 1, established directly on
# the TRUTHFUL closed-form gamma_n (not on the approximations).  A doubling
# (Richardson) estimator removes the additive constant b, exposing the leading
# log coefficient.  If the true coefficient were anything other than 1, this
# would converge to that other value and the assertion would fail.
# ---------------------------------------------------------------------------


def _defect_doubling(n: int) -> mpmath.mpf:
    """n*defect = c log n + b + O(1/n); Richardson removes b, returns c."""
    mpmath.mp.dps = PREC
    return (2 * n * _defect(2 * n) - n * _defect(n)) / mpmath.log(2)


def test_net_defect_coefficient_converges_to_one():
    """defect log-coefficient -> 1, derived (= 2 - 1), on the true gamma_n.

    This is the operational content of lim n*gamma_n/log n = 4/pi:
    gamma_n = (4/pi) defect_n and defect coeff -> 1  =>  4/pi.
    """
    cs = {n: float(_defect_doubling(n)) for n in [50, 100, 200]}
    # Monotone convergence toward 1 with deviation roughly halving per doubling.
    assert abs(cs[100] - 1.0) < abs(cs[50] - 1.0)
    assert abs(cs[200] - 1.0) < abs(cs[100] - 1.0)
    # At n=200 the doubling estimator is within ~1.5% of 1.0.
    assert abs(cs[200] - 1.0) < 0.02, f"net defect coeff (n=200) = {cs[200]}"


def test_gamma_doubling_converges_to_4_over_pi():
    """n*gamma_n/log n -> 4/pi: doubling estimator on the TRUE gamma_n.

    a_n := (2n gamma_{2n} - n gamma_n)/log 2 -> 4/pi.  This is exactly
    (4/pi) times the defect doubling estimator; we check the *value* lands on
    4/pi and NOT on a nearby false constant (2/pi, the circle-Fejer value, is
    explicitly excluded).
    """
    mpmath.mp.dps = PREC
    target = 4.0 / math.pi
    decoy = 2.0 / math.pi  # the circle-Fejer (unweighted) constant

    def a(n):
        return float((2 * n * _gamma(2 * n) - n * _gamma(n)) / mpmath.log(2))

    a50, a100, a200 = a(50), a(100), a(200)
    # Monotone approach to 4/pi from above.
    assert a200 < a100 < a50
    assert abs(a200 - target) < abs(a100 - target) < abs(a50 - target)
    assert abs(a200 - target) < 0.02, f"a_200 = {a200}, target 4/pi = {target}"
    # Decisively closer to 4/pi than to the circle 2/pi decoy.
    assert abs(a200 - target) < 0.25 * abs(a200 - decoy), (
        f"a_200={a200} not clearly resolving 4/pi vs 2/pi decoy"
    )


def test_constant_is_not_circle_fejer_2_over_pi():
    """Sanity guard: the SU(2) constant is 4/pi, NOT the circle-Fejer 2/pi.

    The SU(2) constant is twice the circle constant (an OBSERVATION -- which
    of the sin^2(chi/2) conjugacy-class weight and the sqrt(2j+1) Plancherel
    weight supplies the factor 2 is not isolated here).  This test documents
    that the derivation discriminates between the two candidate constants --
    i.e. it COULD have failed had the analysis been wrong.
    """
    mpmath.mp.dps = PREC
    a200 = float((2 * 200 * _gamma(400) - 200 * _gamma(200)) / mpmath.log(2))
    four_over_pi = 4.0 / math.pi
    two_over_pi = 2.0 / math.pi
    assert abs(a200 - four_over_pi) < abs(a200 - two_over_pi)


# ---------------------------------------------------------------------------
# The CIRCLE constant, computed (Paper 38 Remark "Connection to the circle
# Fejer estimate", eq. circle_fejer_moment; added at the 2026-09-02 trunk
# certification run).  Paper 38 formerly stated the circle Fejer first-moment
# constant as 4/pi, "the same on both sides"; the probability-normalised
# circle constant is 2/pi.  Two independent routes:
#   (i)  quadrature of the kernel moment  m_n = (1/pi) int_0^pi theta F_n dtheta,
#        F_n = (1/n) sin^2(n theta/2)/sin^2(theta/2);
#   (ii) the exact closed form  m_n = pi/2 - (4/pi) sum_{k odd < n} (1 - k/n)/k^2
#        (from int_0^pi theta cos(k theta) dtheta = ((-1)^k - 1)/k^2).
# The doubling estimator (2n m_{2n} - n m_n)/log 2 -> 2/pi on both, and the
# ratio to the SU(2) estimator -> 2.
# ---------------------------------------------------------------------------


def _circle_moment_closed_form(n: int) -> mpmath.mpf:
    mpmath.mp.dps = PREC
    s = sum((1 - mpmath.mpf(k) / n) / k ** 2 for k in range(1, n, 2))
    return mpmath.pi / 2 - (4 / mpmath.pi) * s


def _circle_moment_quadrature(n: int) -> mpmath.mpf:
    mpmath.mp.dps = 30

    def F(theta):
        if theta == 0:
            return mpmath.mpf(n)
        return mpmath.sin(n * theta / 2) ** 2 / (n * mpmath.sin(theta / 2) ** 2)

    # Subdivide at the kernel's zeros (theta = 2 pi m / n) so the oscillatory
    # integrand is resolved; the integrand is smooth on each panel.
    pts = [2 * mpmath.pi * m / n for m in range(0, n // 2 + 1)]
    if pts[-1] < mpmath.pi:
        pts.append(mpmath.pi)
    return mpmath.quad(lambda t: t * F(t), pts) / mpmath.pi


@pytest.mark.parametrize("n", [2, 3, 8, 25, 64])
def test_circle_fejer_closed_form_matches_quadrature(n):
    """Route (ii) closed form == route (i) quadrature of the actual kernel."""
    cf = _circle_moment_closed_form(n)
    q = _circle_moment_quadrature(n)
    assert abs(cf - q) < mpmath.mpf("1e-18"), f"n={n}: closed {cf} vs quad {q}"


def test_circle_fejer_kernel_is_probability_normalised():
    """(1/2pi) int F_n = 1 -- the normalisation under which the constant is 2/pi."""
    mpmath.mp.dps = 30
    for n in [3, 10, 40]:
        pts = [2 * mpmath.pi * m / n for m in range(0, n // 2 + 1)]
        if pts[-1] < mpmath.pi:
            pts.append(mpmath.pi)
        F = lambda t: (mpmath.mpf(n) if t == 0 else
                       mpmath.sin(n * t / 2) ** 2 / (n * mpmath.sin(t / 2) ** 2))
        total = 2 * mpmath.quad(F, pts) / (2 * mpmath.pi)
        assert abs(total - 1) < mpmath.mpf("1e-18"), f"n={n}: mass {total}"


def test_circle_fejer_constant_is_2_over_pi():
    """Doubling estimator on the circle moment -> 2/pi (and NOT 4/pi).

    n m_n = c log n + b + O(1/n); (2n m_{2n} - n m_n)/log 2 -> c.  Computed
    on the exact closed form; the value must land on 2/pi and be decisively
    closer to it than to the SU(2) 4/pi.
    """
    mpmath.mp.dps = PREC
    target = 2.0 / math.pi
    decoy = 4.0 / math.pi

    def a(n):
        return float((2 * n * _circle_moment_closed_form(2 * n)
                      - n * _circle_moment_closed_form(n)) / mpmath.log(2))

    a100, a200, a400 = a(100), a(200), a(400)
    assert abs(a400 - target) < abs(a200 - target) < abs(a100 - target)
    # Measured |a_400 - 2/pi| = 7.2e-7 (2026-09-02); pinned to 1e-5 so the
    # printed "0.63662 at n = 400" (Paper 38 rem:circle_fejer) is backed to
    # all five digits (3e-6 keeps the rounding at 0.63662; 1e-5 did not).
    assert abs(a400 - target) < 3e-6, f"a_400 = {a400}, target 2/pi = {target}"
    assert abs(a400 - target) < 0.05 * abs(a400 - decoy)


def test_su2_constant_is_twice_circle():
    """SU(2) doubling estimator / circle doubling estimator -> 2 (OBSERVATION).

    Both estimators are computed on their exact closed forms; neither 4/pi
    nor 2/pi is assumed.
    """
    mpmath.mp.dps = PREC

    def a_su2(n):
        return (2 * n * _gamma(2 * n) - n * _gamma(n)) / mpmath.log(2)

    def a_circle(n):
        return (2 * n * _circle_moment_closed_form(2 * n)
                - n * _circle_moment_closed_form(n)) / mpmath.log(2)

    r100 = float(a_su2(100) / a_circle(100))
    r200 = float(a_su2(200) / a_circle(200))
    r800 = float(a_su2(800) / a_circle(800))
    assert abs(r200 - 2.0) < abs(r100 - 2.0), "ratio not converging toward 2"
    assert abs(r200 - 2.0) < 0.03, f"SU(2)/circle estimator ratio = {r200}"
    # Paper 38 rem:circle_fejer prints "2.008 at n = 800" (measured 2.008248,
    # 2026-09-02); pin it here so the printed value has a frozen backing.
    assert abs(r800 - 2.008248) < 5e-4, f"SU(2)/circle estimator ratio at n=800 = {r800} (printed 2.008)"
    assert abs(r800 - 2.0) < abs(r200 - 2.0), "ratio not converging toward 2 at n=800"


def test_su2_doubling_estimator_at_n800():
    """Pin the printed SU(2) doubling-estimator value (Paper 38 rem:circle_fejer
    and App. A): a_800 = (1600 gamma_1600 - 800 gamma_800)/log 2 = 1.278491,
    three digits of 4/pi = 1.273240, approached from above (residual 5.25e-3,
    shrinking ~1.85x per doubling).  Added 2026-09-02 (trunk DELTA #2): the
    value had been printed as "a_1600" under an index-convention slip.
    """
    mpmath.mp.dps = PREC
    a800 = float((2 * 800 * _gamma(1600) - 800 * _gamma(800)) / mpmath.log(2))
    assert abs(a800 - 1.278491) < 2e-4, f"a_800 = {a800}"
    assert 0.0 < a800 - 4.0 / math.pi < 6e-3, f"a_800 residual {a800 - 4/math.pi}"


@pytest.mark.slow
def test_su2_doubling_estimator_at_n1600():
    """a_1600 = (3200 gamma_3200 - 1600 gamma_1600)/log 2 = 1.276067 (residual
    2.83e-3; ratio to the a_800 residual ~1.86).  Slow: gamma_3200 ~25 s."""
    mpmath.mp.dps = PREC
    g1600, g3200 = _gamma(1600), _gamma(3200)
    a800 = float((2 * 800 * g1600 - 800 * _gamma(800)) / mpmath.log(2))
    a1600 = float((2 * 1600 * g3200 - 1600 * g1600) / mpmath.log(2))
    assert abs(a1600 - 1.276067) < 2e-4, f"a_1600 = {a1600}"
    shrink = (a800 - 4.0 / math.pi) / (a1600 - 4.0 / math.pi)
    assert 1.7 < shrink < 2.0, f"residual shrink per doubling {shrink}"
