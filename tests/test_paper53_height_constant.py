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

    Normalisation-independent by construction, so it is exactly the sharp
    gradient gain. Integrated half-wave by half-wave between consecutive zeros
    of J_{s+1}, so the sign changes are resolved rather than averaged away --
    integrating straight through would cancel them and manufacture the very
    answer the paper claimed.
    """
    zeros = _bessel_zeros(s + 1.0, n_waves)
    edges = np.concatenate(([1e-9], zeros))

    def f(r: float) -> float:                       # radial measure folded in
        return special.jv(s + 1.0, r) / r ** (s + 1.0) * r

    signed = absol = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        v, _ = integrate.quad(f, a, b, limit=200)
        signed += v
        absol += abs(v)
    tail, _ = integrate.quad(lambda r: abs(f(r)), edges[-1], np.inf, limit=200)
    return (absol + tail) / signed


@pytest.mark.slow
def test_paper53_bochner_riesz_is_not_gradient_non_expansive():
    """Rejects (a): the Lebesgue constant exceeds 1 at every finite order."""
    for s, expected in ((0.75, 3.23), (1.0, 2.01), (2.0, 1.23)):
        L = lebesgue_constant(s)
        assert L > 1.0 + 1e-3, f"s={s}: Lebesgue constant {L} claims non-expansive"
        assert abs(L - expected) < 0.05, f"s={s}: {L} drifted from {expected}"


@pytest.mark.slow
def test_paper53_height_constant_is_monotone_toward_one_only_at_infinity():
    """Rejects (c): it decreases in s but never reaches 1 at finite order."""
    orders = [0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    Ls = [lebesgue_constant(s) for s in orders]
    assert all(x > y for x, y in zip(Ls, Ls[1:])), f"not monotone decreasing: {Ls}"
    assert Ls[-1] > 1.0, f"reached non-expansive at finite s: {Ls[-1]}"
    # ... and the approach is genuine, not a plateau above some other value
    assert Ls[-1] < 1.10, f"s=5 constant {Ls[-1]} is not approaching 1"


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
