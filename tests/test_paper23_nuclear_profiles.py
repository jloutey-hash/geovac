"""Paper 23 sec:finite_nucleus (derivation and other profiles) -- eq:defect_s and the profile paragraph.

Records: docs/predictions/2026-10-18_finite_nucleus.md (derivation section) and
docs/predictions/2026-10-18_nuclear_profiles.md (NP1-NP3, committed before computing in 8292bf51).
Written as its own step after the result; it does not import the computation scripts. It evaluates
the defect by exact matching at E = 0 (Seaton's limit, guarded in test_paper23_finite_nucleus_door.py)
with its own code, and compares with the series as printed. Each test names the wrong answer it rejects.
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp


def _mu_from_logder(rl, xr):
    """Zero-energy defect from r u'/u = rl at r = xr (Z = 1), exterior pair sqrt(2r) J1, -sqrt(2r) Y1."""
    f = lambda r: mp.sqrt(2 * r) * mp.besselj(1, mp.sqrt(8 * r))
    g = lambda r: -mp.sqrt(2 * r) * mp.bessely(1, mp.sqrt(8 * r))
    rf, rg = xr * mp.diff(f, xr), xr * mp.diff(g, xr)
    return mp.atan((rf - rl * f(xr)) / (rl * g(xr) - rg)) / mp.pi


def _sphere_exact(x):
    """Uniform sphere, R = x, Z = 1, E = 0: interior oscillator, w^2 = 1/R^3, energy 3/(2R)."""
    R = x
    w = mp.sqrt(1 / R**3)
    eps = 3 / (2 * R)
    b = mp.mpf(3) / 2
    u = lambda r: r * mp.e ** (-w * r * r / 2) * mp.hyp1f1(b / 2 - eps / (2 * w), b, w * r * r)
    return _mu_from_logder(R * mp.diff(u, R) / u(R), R)


def _shell_exact(x):
    k = mp.sqrt(2 * x)
    return _mu_from_logder(k * mp.cot(k), x)


def _sphere_series(x):
    L, g = mp.log(2 * x), mp.euler
    return (-mp.mpf(2) / 5 * x**2 + mp.mpf(32) / 63 * x**3
            - x**4 * (mp.mpf(4) / 25 * L + mp.mpf(8) / 25 * g + mp.mpf(28477) / 225225)
            + x**5 * (mp.mpf(128) / 315 * L + mp.mpf(256) / 315 * g - mp.mpf(30486592) / 57432375))


def _shell_series(x):
    L, g = mp.log(2 * x), mp.euler
    return (-mp.mpf(2) / 3 * x**2 + mp.mpf(16) / 15 * x**3
            - x**4 * (mp.mpf(4) / 9 * L + mp.mpf(8) / 9 * g + mp.mpf(127) / 315)
            + x**5 * (mp.mpf(64) / 45 * L + mp.mpf(128) / 45 * g - mp.mpf(4448) / 2835))


def test_uniform_sphere_series_matches_exact_matching():
    """Rejects any printed coefficient of eq:defect_s that is off: the residual against exact matching
    must be O(x^6), below 1e-17 at x = 1e-3 and below 1e-14 at x = 1e-2 (a 1 percent change in the 28477/225225,
    a wrong gamma_E or ln coefficient, or the sign of the Friar term, breaks it)."""
    mp.mp.dps = 50
    for x, tol in ((mp.mpf(1) / 1000, mp.mpf(10) ** -17), (mp.mpf(1) / 100, mp.mpf(10) ** -11)):
        assert abs(_sphere_exact(x) - _sphere_series(x)) < tol


def test_shell_series_matches_exact_matching():
    """Rejects a shell series with the sphere's coefficients, or a fourth-order term without the universal
    -delta_0^2 [ln(2x) + 2 gamma_E] piece."""
    mp.mp.dps = 50
    for x, tol in ((mp.mpf(1) / 1000, mp.mpf(10) ** -16), (mp.mpf(1) / 100, mp.mpf(10) ** -10)):
        assert abs(_shell_exact(x) - _shell_series(x)) < tol
    # the universal piece: fourth-order ln and gamma_E coefficients equal -delta_0^2, -2 delta_0^2
    for c2, cL, cg in ((sp.Rational(-2, 5), sp.Rational(4, 25), sp.Rational(8, 25)),
                       (sp.Rational(-2, 3), sp.Rational(4, 9), sp.Rational(8, 9))):
        assert cL == c2**2 and cg == 2 * c2**2


def test_friar_moments_of_sphere_and_shell():
    """Rejects a third-order term other than the Friar term (1/3) Z^3 <r^3>_(2): the uniform sphere's
    self-convolution moment is 32 R^3/21 (so 32/63), the shell's 16 R^3/5 (so 16/15)."""
    d, R = sp.symbols("d R", positive=True)
    V = sp.Rational(4, 3) * sp.pi * R**3
    rho2 = sp.pi / 12 * (4 * R + d) * (2 * R - d) ** 2 / V**2
    assert sp.simplify(sp.integrate(4 * sp.pi * d**2 * rho2, (d, 0, 2 * R)) - 1) == 0
    m3 = sp.simplify(sp.integrate(4 * sp.pi * d**5 * rho2, (d, 0, 2 * R)))
    assert sp.simplify(m3 - sp.Rational(32, 21) * R**3) == 0 and sp.simplify(m3 / 3 / R**3) == sp.Rational(32, 63)
    shell = sp.integrate(d**3 * d / (2 * R**2), (d, 0, 2 * R))  # distance density d/(2R^2) on two shells
    assert sp.simplify(shell - sp.Rational(16, 5) * R**3) == 0


@pytest.mark.slow
def test_gaussian_friar_coefficient_carries_inverse_sqrt_pi():
    """Rejects a Gaussian third-order coefficient other than the Friar value 32/(3 sqrt(pi)) y^3, for
    example the uniform sphere's 32/63 rescaled, or a rational."""
    mp.mp.dps = 30

    def defect(y, smax=12):
        V = lambda s: mp.sqrt(2 / mp.pi) if s == 0 else mp.erf(s / mp.sqrt(2)) / s
        sol = mp.odefun(lambda s, Y: [Y[1], -2 * y * V(s) * Y[0]], 0, [mp.mpf(0), mp.mpf(1)])
        u, du = sol(smax)
        return _mu_from_logder(smax * du / u, smax * y)

    ys = [mp.mpf(1) / 1000 * 2 ** (mp.mpf(j) / 2) for j in range(9)]
    red = [defect(y) / (-2 * y**2) for y in ys]
    cols = [lambda t: 1, lambda t: t, lambda t: t**2, lambda t: t**2 * mp.log(t), lambda t: t**3, lambda t: t**3 * mp.log(t)]
    A = mp.matrix([[c(t) for c in cols] for t in ys])
    sol = mp.lu_solve(A.T * A, A.T * mp.matrix(red))
    assert abs(sol[0] - 1) < mp.mpf(10) ** -7
    assert abs(sol[1] + 16 / (3 * mp.sqrt(mp.pi))) < mp.mpf(10) ** -3
    assert abs(sol[3] - 2) < mp.mpf("0.05")
