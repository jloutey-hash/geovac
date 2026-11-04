"""Paper 53, Remark rem:height_constant -- the Bochner-Riesz height constant.

The paper asserted the plane Berezin map is *gradient-non-expansive*,
||grad B f||_inf <= ||grad f||_inf, and used that as the height of its
pointed/proper propinquity assembly. It is false over the unit-Lipschitz ball.

B_Lambda is a radial convolution, so grad(B f) = B(grad f) and Young's
inequality makes the SHARP constant the kernel's L1 norm (its Lebesgue
constant), which equals 1 only when the kernel is non-negative a.e. The 2D
Bochner-Riesz kernel K_s(r) ~ J_{s+1}(r)/r^{s+1} oscillates, so it is negative
somewhere and its L1 norm exceeds 1 strictly at every finite Cesaro order s.

The wrong answers this test rejects:
  (a) "the map is gradient-non-expansive" -- i.e. Lebesgue constant <= 1;
  (b) "the ratio rises to 1 as Lambda -> infinity" -- the constant is
      Lambda-INDEPENDENT by scaling, so nothing about it moves with Lambda;
  (c) a Lebesgue constant that decreases to 1 at FINITE s -- it approaches 1
      only as s -> infinity.

Deliberately NOT asserted: what this costs the assembly. That is flagged
[OPEN -- PI adjudication] in the paper and is not a numerical question.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")
from scipy import integrate, optimize, special  # noqa: E402


def _bessel_zeros(v: float, n: int, step: float = 0.05) -> np.ndarray:
    """First `n` positive zeros of J_v (scipy's jn_zeros is integer-order only)."""
    rmax = (n + v) * np.pi + 20.0
    r = np.arange(1e-6, rmax, step)
    f = special.jv(v, r)
    idx = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    return np.array([optimize.brentq(lambda x: special.jv(v, x), r[i], r[i + 1])
                     for i in idx[:n]])


def lebesgue_constant(s: float, n_waves: int = 200) -> float:
    """||K_s||_1 / (int K_s) for the 2D Bochner-Riesz kernel at Cesaro order s.

    In the radial measure the integrand is J_{s+1}(x) x^{-s}.  It is summed
    half-wave by half-wave between consecutive zeros of J_{s+1}, so the sign
    changes are resolved rather than averaged away, and the |.| tail past the
    last zero A is added analytically: |J_nu(x)| averages (2/pi) sqrt(2/(pi x))
    over a period, giving (2/pi) sqrt(2/pi) A^{1/2-s} / (s - 1/2).  The
    denominator is the closed form int_0^inf J_{s+1}(x) x^{-s} dx
    = 1 / (2^s Gamma(s+1)).

    Integrating the |.| tail numerically to infinity is what produced the
    retracted table (5.72 at s = 0.6): quad truncates a tail that decays only
    like A^{1/2-s}.
    """
    zeros = _bessel_zeros(s + 1.0, n_waves)
    edges = np.concatenate(([1e-12], zeros))

    def f(x: float) -> float:
        return special.jv(s + 1.0, x) * x ** (-s)

    absol = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        v, _ = integrate.quad(f, a, b, limit=200)
        absol += abs(v)
    A = zeros[-1]
    tail = (2.0 / np.pi) * np.sqrt(2.0 / np.pi) * A ** (0.5 - s) / (s - 0.5)
    signed = 1.0 / (2.0 ** s * special.gamma(s + 1.0))
    return (absol + tail) / signed


# Paper 53, Sec. not_br table: (order s, printed value, half a unit in its last digit).
PRINTED = ((0.6, 7.10, 0.005), (0.75, 3.29, 0.005), (1.0, 2.03, 0.005),
           (1.5, 1.42, 0.005), (2.0, 1.23, 0.005), (3.0, 1.09, 0.005),
           (5.0, 1.02, 0.005), (8.0, 1.0045, 0.00005))


@pytest.mark.parametrize("s,printed,half_ulp", PRINTED)
def test_paper53_lebesgue_table_entry(s, printed, half_ulp):
    """Every printed L_s rounds to the recomputed value (rejects the old
    numerically-truncated tail: 5.72 at s = 0.6, 3.23 at s = 0.75)."""
    L = lebesgue_constant(s)
    assert abs(L - printed) <= half_ulp, f"s={s}: recomputed {L:.5f} vs printed {printed}"


def test_paper53_lebesgue_table_stable_in_zero_count():
    """The 'stable to the digits shown from 200 to 3200 zeros' clause, at the
    two orders where the tail is largest."""
    for s in (0.6, 0.75):
        assert abs(lebesgue_constant(s, 200) - lebesgue_constant(s, 1600)) < 1e-3


def test_paper53_bochner_riesz_is_not_gradient_non_expansive():
    """Rejects (a), (c): L_s > 1 at every tabulated finite order, strictly
    decreasing in s, and approaching 1 (not plateauing above it)."""
    Ls = [lebesgue_constant(s) for s, _, _ in PRINTED]
    assert all(x > 1.0 + 1e-3 for x in Ls), Ls
    assert all(x > y for x, y in zip(Ls, Ls[1:])), f"not monotone decreasing: {Ls}"
    assert Ls[-1] < 1.01, f"s=8 constant {Ls[-1]} is not approaching 1"


def test_paper53_height_constant_does_not_depend_on_lambda():
    """Rejects (b): the paper reported a ratio 'rising to 1 as Lambda -> inf'.

    The kernel at cutoff Lambda is K_s(Lambda r) Lambda^2, a pure rescaling, so
    its L1 norm is invariant. Verified directly on the scaled integrand rather
    than argued: any Lambda-dependence in a measured ratio is a property of the
    test function, not of the operator.
    """
    s = 1.0
    zeros = _bessel_zeros(s + 1.0, 60)

    def l1_at(lam: float) -> float:
        edges = np.concatenate(([1e-9], zeros / lam))

        def f(r: float) -> float:
            x = lam * r
            return special.jv(s + 1.0, x) / x ** (s + 1.0) * lam ** 2 * r

        signed = absol = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            v, _ = integrate.quad(f, a, b, limit=200)
            signed += v
            absol += abs(v)
        return absol / signed

    ref = l1_at(1.0)
    for lam in (2.0, 8.0, 40.0):
        assert abs(l1_at(lam) - ref) < 1e-6, (
            f"Lebesgue constant moved with Lambda={lam}: {l1_at(lam)} vs {ref}")
