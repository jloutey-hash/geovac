"""Paper 38, "Other amplitudes" (thm:main_unconditional): a tapered amplitude gives O(1/n).

The lifted-state proof uses the Fejer amplitude only as a real central amplitude
in the spinor window, so d_GH <= gamma(a) for every real a_j.  This test backs the
[MEASURED] statement that the Cesaro taper a_j = sqrt(2j+1) (1 - 2j/(n+1)) gives
n * gamma(a) bounded (about 6.35 .. 6.86 for n = 25 .. 400), below the Fejer
moment, whose n * gamma grows like (4/pi) log n.  So 4/pi is the Fejer kernel's
constant, not the convergence rate.

Independent route: the moment is integrated here from the characters with numpy /
scipy quadrature, not via geovac.central_fejer_su2 (one cross-check against
gamma_rate(use_cesaro=True) at small n ties the two routes together).

Rejected wrong answers:
  * "the rate is (4/pi) log n/n for every admissible amplitude" -- then
    n*gamma(taper) would grow by ~ (4/pi) ln 4 ~ 1.77 between n = 25 and n = 100;
  * "the taper is no better than Fejer" -- then gamma(taper) >= gamma(Fejer);
  * "the taper still carries a smaller log rate, c log n/n with 0 < c < 4/pi" --
    then the doubling differences n*gamma(2n) - n*gamma(n) plateau at c ln 2
    instead of halving (fire-tested 2026-10-03 with a = 0.75 taper + 0.25,
    which plateaus at 0.196).
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad


@pytest.fixture(autouse=True)
def _restore_mp_dps():
    """gamma_rate sets mpmath's global precision; restore it (xdist hygiene)."""
    import mpmath
    saved = mpmath.mp.dps
    yield
    mpmath.mp.dps = saved


def _gamma(n: int, taper: bool) -> float:
    """First moment of |sum_j a_j chi_j|^2 / sum a_j^2 against the rotation angle
    chi on SU(2) (dual-Coxeter metric), Haar measure sin^2(chi/2)/pi dchi."""
    js = np.arange(n) / 2.0
    a = np.sqrt(2 * js + 1)
    if taper:
        a = a * np.clip(1 - 2 * js / (n + 1), 0.0, None)
    z = float(np.sum(a ** 2))
    w = 2 * js + 1

    def f(c: float) -> float:
        s = math.sin(c / 2)
        h = float(np.sum(a * np.sin(w * c / 2))) / s if s > 1e-14 else float(np.sum(a * w))
        return h * h / z * c * s * s / math.pi

    pts = np.linspace(0.0, 2 * math.pi, 4 * n + 1)
    return sum(quad(f, pts[i], pts[i + 1], limit=200)[0] for i in range(len(pts) - 1))


def test_independent_route_matches_module_at_small_n() -> None:
    from geovac.central_fejer_su2 import gamma_rate
    for n in (4, 7):
        assert _gamma(n, True) == pytest.approx(float(gamma_rate(n, prec=30, use_cesaro=True)), rel=1e-8)
        assert _gamma(n, False) == pytest.approx(float(gamma_rate(n, prec=30)), rel=1e-8)


@pytest.mark.parametrize("n", [3, 5, 10, 25, 50])
def test_taper_beats_fejer(n: int) -> None:
    assert _gamma(n, True) < _gamma(n, False)


def _n_gamma_sum_rule(n: int, taper: bool) -> float:
    """n * gamma(a) in closed form.  With b_k = a_j the coefficient of sin(k chi/2),
    k = 2j + 1 = 1..n:  gamma = sum_{k,l} b_k b_l [J(k-l) - J(k+l)] / (2 pi sum b^2),
    J(m) = int_0^{2pi} chi cos(m chi/2) d chi = 2 pi^2 (m = 0), -8/m^2 (m odd), 0 (m even)."""
    w = np.arange(1, n + 1, dtype=float)
    a = np.sqrt(w)
    if taper:
        a = a * np.clip(1 - (w - 1) / (n + 1), 0.0, None)
    k, l = np.meshgrid(w, w, indexing="ij")

    def jm(m: np.ndarray) -> np.ndarray:
        m = np.abs(m)
        safe = np.where(m == 0, 1.0, m)
        return np.where(m == 0, 2 * math.pi ** 2, np.where(m % 2 == 1, -8.0 / safe ** 2, 0.0))

    g = float((a[:, None] * a[None, :] * (jm(k - l) - jm(k + l))).sum()) / (2 * math.pi * float((a ** 2).sum()))
    return n * g


def test_sum_rule_matches_quadrature() -> None:
    for taper in (True, False):
        assert _n_gamma_sum_rule(25, taper) == pytest.approx(25 * _gamma(25, taper), rel=1e-10)


def test_taper_doubling_differences_vanish() -> None:
    """O(1/n), not c log n / n for c >~ 0.03 (over n <= 1600).  For gamma ~ c log n / n the doubling
    difference T(2n) - T(n), T = n * gamma, tends to c ln 2; for O(1/n) it tends to 0.
    The taper's differences halve at each doubling (0.144, 0.078, 0.041, 0.021, 0.011
    for n = 50..800); an amplitude with a residual log rate plateaus instead."""
    t = {n: _n_gamma_sum_rule(n, True) for n in (50, 100, 200, 400, 800, 1600)}
    d = [t[2 * n] - t[n] for n in (50, 100, 200, 400, 800)]
    for d0, d1 in zip(d, d[1:]):
        assert 0 < d1 < 0.6 * d0
    assert d[-1] < 0.02
    # the printed values (Paper 38, "Other amplitudes")
    for n, v in ((50, 6.60), (200, 6.82), (400, 6.86)):
        assert t[n] == pytest.approx(v, abs=0.006)


def test_tapered_moment_is_order_one_over_n() -> None:
    vals = {n: n * _gamma(n, True) for n in (25, 50, 100)}
    assert 6.0 < vals[25] < 7.0
    assert 6.0 < vals[100] < 7.0
    # bounded, not logarithmic: a (4/pi) log n law would add ~1.77 from 25 to 100
    assert vals[100] - vals[25] < 0.6
    # whereas the Fejer moment does grow like log n over the same range
    fej = {n: n * _gamma(n, False) for n in (25, 100)}
    assert fej[100] - fej[25] > 1.0
