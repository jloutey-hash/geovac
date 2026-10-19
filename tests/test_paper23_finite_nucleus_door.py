"""Paper 23 sec:finite_nucleus -- a finite nucleus releases n from being a whole number.

Record: docs/predictions/2026-10-18_finite_nucleus.md (FN0-FN4, committed before the computation in
bd2435d1, with the FN3 amendment 5ba129aa). Written as its own step after the result; it does not
import the computation script, and it takes every derivative by numerical differentiation where the
computation used the Kummer and Whittaker recurrences. Each test names the wrong answer it rejects.

Model: a uniformly charged sphere of radius R; atomic units; x = ZR. Inside, the electron's radial
equation is an isotropic oscillator, omega^2 = Z/R^3, energy shifted by 3Z/(2R); outside, Coulomb.
"""
from __future__ import annotations

import mpmath as mp
import pytest


def _u_in(r, E, l, Z, R):
    w = mp.sqrt(Z / R**3)
    eps = E + 3 * Z / (2 * R)
    b = l + mp.mpf(3) / 2
    return r ** (l + 1) * mp.e ** (-w * r * r / 2) * mp.hyp1f1(b / 2 - eps / (2 * w), b, w * r * r)


def _logder(f, r):
    return mp.diff(f, r) / f(r)


def _nu(n: int, l: int, Z, R, guess):
    def g(nu):
        E = -Z**2 / (2 * nu**2)
        return (_logder(lambda r: _u_in(r, E, l, Z, R), R)
                - _logder(lambda r: mp.whitw(nu, l + mp.mpf(1) / 2, 2 * Z * r / nu), R))
    return mp.findroot(g, (n - guess, n - mp.mpf("1.001") * guess))


def _defect(n, l, x, Z=1):
    lead = {0: -mp.mpf(2) / 5 * x**2, 1: -x**4 / 105 * (1 - mp.mpf(1) / n**2)}.get(l)
    return n - _nu(n, l, mp.mpf(Z), x / Z, lead)


def _phase_over_pi(k, l, Z, R):
    E = k**2 / 2
    Li = _logder(lambda r: _u_in(r, E, l, Z, R), R)
    eta = -Z / k
    F = mp.coulombf(l, eta, k * R)
    G = mp.coulombg(l, eta, k * R)
    dF = k * mp.diff(lambda p: mp.coulombf(l, eta, p), k * R)
    dG = k * mp.diff(lambda p: mp.coulombg(l, eta, p), k * R)
    return mp.atan((dF - Li * F) / (Li * G - dG)) / mp.pi


def test_fn2_s_wave_defect_is_two_fifths_and_the_same_for_every_n():
    """Rejects a leading coefficient other than 2/5 (for example 2/3, the coefficient of <r^2> applied
    to R^2 without the uniform sphere's 3/5), an n-dependent s-wave defect, and a next-order term
    other than -(80/63) ZR."""
    mp.mp.dps = 40
    x = mp.mpf(1) / 1000
    d1, d5 = _defect(1, 0, x), _defect(5, 0, x)
    assert abs(d1 - d5) < mp.mpf(10) ** -7 * abs(d1)
    ratio = d1 / (-mp.mpf(2) / 5 * x**2)
    assert abs(ratio - (1 - mp.mpf(80) / 63 * x)) < mp.mpf(10) ** -5


def test_fn2_p_wave_defect_carries_the_product():
    """Rejects a p-wave defect without the factor (1 - 1/n^2), or with c_1 other than 1/105."""
    mp.mp.dps = 40
    x = mp.mpf(1) / 1000
    d2, d5 = _defect(2, 1, x), _defect(5, 1, x)
    assert abs(d2 / d5 - (mp.mpf(3) / 4) / (mp.mpf(24) / 25)) < mp.mpf(10) ** -5
    assert abs(d2 / (-x**4 / 105 * mp.mpf(3) / 4) - 1) < mp.mpf(10) ** -3


def test_fn1_defect_depends_on_ZR_only():
    """Rejects a defect that depends on Z and R separately: Z = 2 at R/2 equals Z = 1 at R."""
    mp.mp.dps = 40
    x = mp.mpf(1) / 500
    assert abs(_defect(3, 0, x, Z=1) - _defect(3, 0, x, Z=2)) < mp.mpf(10) ** -25


def test_fn1_node_label_survives():
    """Rejects a node count other than n - l - 1 at ZR = 0.01 for 2s and 3p (the exterior function's
    sign is matched to the interior's at r = R, so a convention flip is not counted as a node)."""
    mp.mp.dps = 25
    x = mp.mpf(1) / 100
    for n, l in ((2, 0), (3, 1)):
        nu = n - _defect(n, l, x)
        E = -1 / (2 * nu**2)
        out = lambda r: mp.whitw(nu, l + mp.mpf(1) / 2, 2 * r / nu)
        sgn = mp.sign(_u_in(x, E, l, 1, x)) * mp.sign(out(x))
        prev, cnt = None, 0
        grid = [x / 2] + [mp.mpf(i) * 4 * n * n / 1500 for i in range(1, 1500)]  # one point inside
        for r in grid:
            s = mp.sign(_u_in(r, E, l, 1, x) if r < x else sgn * out(r))
            if prev is not None and s != prev and s != 0:
                cnt += 1
            prev = s
        assert cnt == n - l - 1


def test_fn3_seaton_limit_across_threshold():
    """Rejects a discontinuity of the defect at threshold: the bound defects extrapolated to
    n -> infinity equal the phase shift over pi extrapolated to k -> 0."""
    mp.mp.dps = 40
    x = mp.mpf(1) / 100
    ns = [10, 20, 40]
    ds = [_defect(n, 0, x) for n in ns]
    # delta(n) = a + b/n^2 + c/n^4
    A = mp.matrix([[1, 1 / mp.mpf(n) ** 2, 1 / mp.mpf(n) ** 4] for n in ns])
    bound_inf = mp.lu_solve(A, mp.matrix(ds))[0]
    ks = [mp.mpf(1) / 20, mp.mpf(1) / 40, mp.mpf(1) / 60]
    ph = [_phase_over_pi(k, 0, mp.mpf(1), x) for k in ks]
    B = mp.matrix([[1, k**2, k**4] for k in ks])
    cont_0 = mp.lu_solve(B, mp.matrix(ph))[0]
    assert abs(bound_inf - cont_0) < mp.mpf(10) ** -10 * abs(bound_inf)


def test_fn3_continuum_carries_the_J13_factor():
    """Rejects the registered FN3 form (phase over pi = the continued defect, no factor), which the
    record scores as falsified at k = 2, and accepts only the amended form with 1/(1 - e^(-2 pi eta))."""
    mp.mp.dps = 40
    x = mp.mpf(1) / 1000
    bound_ratio = _defect(1, 0, x) / (-mp.mpf(2) / 5 * x**2)
    for k in (mp.mpf(1), mp.mpf(2)):
        eta = 1 / k
        fac = 1 / (1 - mp.e ** (-2 * mp.pi * eta))
        ratio_amended = _phase_over_pi(k, 0, mp.mpf(1), x) / (-mp.mpf(2) / 5 * x**2 * fac)
        assert abs(ratio_amended - bound_ratio) < mp.mpf(10) ** -6
    ratio_registered = _phase_over_pi(mp.mpf(2), 0, mp.mpf(1), x) / (-mp.mpf(2) / 5 * x**2)
    assert abs(ratio_registered - bound_ratio) > mp.mpf("0.03")


@pytest.mark.slow
def test_fn4_s_wave_needs_a_log_with_coefficient_two_fifths():
    """Rejects a pure power series in ZR for the s-wave defect (it stalls at a residual near 1e-8),
    and a (ZR)^2 log(ZR) coefficient other than 2/5 or a next-order coefficient other than -80/63."""
    mp.mp.dps = 60
    xs = [mp.mpf(1) / 4000 * 2 ** (mp.mpf(j) / 2) for j in range(14)]
    red = [_defect(1, 0, x) / (-mp.mpf(2) / 5 * x**2) for x in xs]

    def fit(cols):
        A = mp.matrix([[c(x) for c in cols] for x in xs])
        sol = mp.lu_solve(A.T * A, A.T * mp.matrix(red))
        return sol, mp.norm(A * sol - mp.matrix(red))

    P = lambda k: (lambda t: t**k)
    PL = lambda k: (lambda t: t**k * mp.log(t))
    _, res_pow = fit([P(k) for k in range(8)])
    sol, res_log = fit([P(k) for k in range(7)] + [PL(2), PL(3), PL(4)])
    assert res_pow > mp.mpf(10) ** -10 and res_log < mp.mpf(10) ** -13
    assert abs(sol[0] - 1) < mp.mpf(10) ** -9
    assert abs(sol[1] + mp.mpf(80) / 63) < mp.mpf(10) ** -7
    assert abs(sol[7] - mp.mpf(2) / 5) < mp.mpf(10) ** -4
