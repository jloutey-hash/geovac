"""Guard for Paper 23 eq:defect_n4: the s-wave quantum defect of a uniformly charged nucleus
depends on n from fourth order in x = ZR,

    delta_n0 - delta_m0 = -(4/25) [f(n) - f(m)] x^4 - (<r^4>/15) (1/n^2 - 1/m^2) + O(x^5 ln x),
    f(n) = psi(n) - ln n + 1/(2n),  <r^4> = 3 x^4 / 7 (uniform sphere, Z = 1),

the Pachucki-Patkos-Yerokhin (2025) Eq. (30) finite-size energy read as a defect.

Independent route from the scoring script debug/eft_checks.py (which uses an interior power series,
the Whittaker function and a hand-written secant): here the interior is integrated by mpmath's Taylor
ODE solver, the exterior is written with Tricomi's U, and the root is found by mpmath.findroot.

Wrong answers this rejects: no n-dependence at x^4 (the old "the same at every n"); the <r^4> term
dropped; the PPY bracket's -1/n left unconverted; a sign error on the psi term.
"""
import mpmath as mp
import pytest

DPS = 40
NS = (1, 2, 3)
XS = ("4e-4", "6e-4", "9e-4", "1.35e-3", "2e-3", "3e-3", "4.5e-3")


def _exterior_logderiv(nu, r):
    # u = W_{nu,1/2}(z), z = 2r/nu, W = e^{-z/2} z U(1-nu, 2, z)
    z = 2 * r / nu
    a = 1 - nu
    dlnW_dz = -mp.mpf(1) / 2 + 1 / z - a * mp.hyperu(a + 1, 3, z) / mp.hyperu(a, 2, z)
    return (2 / nu) * dlnW_dz


def _interior_logderiv(nu, R):
    E = -1 / (2 * nu**2)
    # u'' = (-(3 R^2 - r^2)/R^3 - 2E) u, u(0) = 0, u'(0) = 1
    f = mp.odefun(lambda r, y: [y[1], (-(3 * R**2 - r**2) / R**3 - 2 * E) * y[0]], 0, [0, 1])
    u, du = f(R)
    return du / u


def _mu(n, x):
    g = lambda nu: _interior_logderiv(nu, x) - _exterior_logderiv(nu, x)
    nu = mp.findroot(g, mp.mpf(n) + mp.mpf(2) / 5 * x**2, tol=mp.mpf(10)**(-(DPS - 8)))
    return n - nu


@pytest.fixture(scope="module")
def mus():
    mp.mp.dps = DPS
    xs = [mp.mpf(s) for s in XS]
    return xs, {n: [_mu(n, x) for x in xs] for n in NS}


def _f(n):
    n = mp.mpf(n)
    return mp.digamma(n) - mp.log(n) + 1 / (2 * n)


def _predicted(n, m):
    return -mp.mpf(4) / 25 * (_f(n) - _f(m)) - mp.mpf(3) / 7 / 15 * (1 / mp.mpf(n)**2 - 1 / mp.mpf(m)**2)


def _fit(xs, ys, cols):
    A = mp.matrix([[c(x) for c in cols] for x in xs])
    sol, _ = mp.qr_solve(A, mp.matrix(ys))
    return [sol[i] for i in range(len(cols))]


COLS = [lambda x: 1, lambda x: x * mp.log(x), lambda x: x, lambda x: x**2 * mp.log(x)]


@pytest.mark.parametrize("m", [2, 3])
def test_defect_depends_on_n_at_fourth_order(mus, m):
    xs, mu = mus
    ys = [(mu[1][i] - mu[m][i]) / xs[i]**4 for i in range(len(xs))]
    c4 = _fit(xs, ys, COLS)[0]
    pred = _predicted(1, m)
    assert abs(pred) > mp.mpf("0.01")
    assert abs(c4 / pred - 1) < mp.mpf("1e-5"), (c4, pred)


def test_no_n_dependence_below_fourth_order(mus):
    xs, mu = mus
    # an x^2 or x^3 n-dependence would make (mu1 - mu2)/x^4 grow like 1/x^2 or 1/x at small x
    for i, x in enumerate(xs):
        assert abs((mu[1][i] - mu[2][i]) / x**4) < mp.mpf("0.02")


def test_no_log_at_fourth_order(mus):
    xs, mu = mus
    ys = [(mu[1][i] - mu[2][i]) / xs[i]**4 for i in range(len(xs))]
    sol = _fit(xs, ys, [lambda x: mp.log(x)] + COLS)
    assert abs(sol[0] / sol[1]) < mp.mpf("1e-4"), sol[:2]
