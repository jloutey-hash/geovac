"""Paper 23 sec:finite_nucleus, "A finite magnetization": the hyperfine factor of an s level.

Backs (i) the leading change -2 Z r_Z with the Zemach radius, the same for every n; (ii) with a point
charge, an entire function of x = ZR whose coefficients are polynomials in 1/n^2, the x^2 one being
1 + 1/(5 n^2) for a uniform-ball magnetization; (iii) with a uniform charge, the (4/5) x^2 ln x term,
-2 times the leading quantum defect -(2/5) x^2, for every magnetization profile and n.
Record: docs/predictions/2026-10-19_zemach_door.md.

Independent of debug/zemach_door.py, which used 40-digit mpmath and Whittaker functions: here the
Zemach radii come from Zemach's Fourier formula r_Z = -(4/pi) int (G_E G_M - 1) dq / q^2 instead of
the convolution, the point-charge series from the Laguerre form of the s density, and the bound
states from float64 inward integration (scipy DOP853) matched to the interior series.
Z = 1, non-relativistic, infinite nuclear mass, spin magnetization only.
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# ---------------------------------------------------------------- (i) Zemach radii, Fourier route
def _G(kind, q):
    if kind == "point":
        return mp.mpf(1)
    if kind == "ball":
        return 3 * (mp.sin(q) - q * mp.cos(q)) / q**3 if q > mp.mpf("1e-3") else 1 - q**2 / 10
    if kind == "shell":
        return mp.sin(q) / q if q > mp.mpf("1e-6") else 1 - q**2 / 6


def zemach_radius(charge, magn, L=400):
    mp.mp.dps = 20
    f = lambda q: (_G(charge, q) * _G(magn, q) - 1) / q**2
    pts = [mp.pi * k for k in range(0, int(L / mp.pi) + 1)]
    body = mp.quad(f, pts)
    return -4 / mp.pi * (body - 1 / pts[-1])     # tail of -1/q^2; the oscillating tail is O(L^-3)


@pytest.mark.parametrize("magn,expected", [("point", 0.75), ("ball", 36 / 35), ("shell", 1.2)])
def test_zemach_radius_uniform_charge(magn, expected):
    assert abs(float(zemach_radius("ball", magn)) - expected) < 1e-6


# ---------------------------------------------------------------- (ii) point charge, exact
r, x = sp.symbols("r x", positive=True)


def point_charge_ball_series(n, order=4):
    dens = sp.exp(-2 * r / n) * sp.assoc_laguerre(n - 1, 1, 2 * r / n)**2 / n**2   # |psi|^2 / |psi(0)|^2
    H = 3 / x**3 * sp.integrate(sp.expand(r**2 * dens), (r, 0, x))
    ser = sp.series(H, x, 0, order + 1).removeO()
    assert not ser.has(sp.log)
    return [sp.nsimplify(ser.coeff(x, k)) for k in range(order + 1)]


def test_point_charge_series_polynomial_in_inverse_n_squared():
    c = {n: point_charge_ball_series(n) for n in range(1, 5)}
    for n in range(1, 5):
        assert c[n][0] == 1
        assert c[n][1] == sp.Rational(-3, 2)
        assert c[n][2] == 1 + sp.Rational(1, 5) / n**2
    t = sp.symbols("t")
    for k, deg in ((3, 1), (4, 2)):
        P = sp.interpolate([(sp.Rational(1, n**2), c[n][k]) for n in range(1, deg + 2)], t)
        for n in range(deg + 2, 5):
            assert sp.simplify(P.subs(t, sp.Rational(1, n**2)) - c[n][k]) == 0, (k, n)


# ---------------------------------------------------------------- (iii) uniform charge, float64
def _interior(R, E, nterms=40):
    c0 = -3.0 / R - 2.0 * E
    a = np.zeros(nterms)
    a[1] = 1.0
    a[3] = c0 * a[1] / 6.0
    for k in range(4, nterms):
        a[k] = (c0 * a[k - 2] + a[k - 4] / R**3) / (k * (k - 1))
    ks = np.arange(nterms)
    uR = np.sum(a * R**ks)
    dR = np.sum(ks[1:] * a[1:] * R**(ks[1:] - 1))
    odd = ks[1::2]
    sq = sum(a[i] * a[j] * R**(i + j + 1) / (i + j + 1) for i in odd for j in odd)
    return uR, dR, sq


def _exterior(R, E, n):
    kappa = math.sqrt(-2.0 * E)
    rmax = 45.0 * n * n

    def rhs(rr, y):
        return [y[1], (-2.0 / rr - 2.0 * E) * y[0], y[0] ** 2]

    sol = solve_ivp(rhs, (rmax, R), [1e-12, -kappa * 1e-12, 0.0], method="DOP853",
                    rtol=1e-13, atol=1e-30)
    u, du, I = sol.y[:, -1]
    return u, du, -I


def _mismatch(E, R, n):
    uR, dR, _ = _interior(R, E)
    u, du, _ = _exterior(R, E, n)
    return dR / uR - du / u


def hyperfine(n, R):
    Ept = -0.5 / n**2
    lo, hi = Ept * (1 + 1e-10), Ept * (1 - 3 * R * R)
    E = brentq(_mismatch, lo, hi, args=(R, n), xtol=1e-18, rtol=1e-15, maxiter=200)
    uR, _, sq_in = _interior(R, E)
    u, _, I = _exterior(R, E, n)
    N = sq_in + (uR / u) ** 2 * I
    n3 = n**3
    return {"point": n3 / (4 * N) - 1,
            "ball": 3 * n3 * sq_in / (4 * R**3 * N) - 1,
            "shell": n3 * uR**2 / (4 * R**2 * N) - 1}


XS = np.logspace(-3.0, -1.6, 12)


@pytest.fixture(scope="module")
def fits():
    out = {}
    for n in (1, 2):
        rows = [hyperfine(n, xx) for xx in XS]
        A = np.column_stack([XS, XS**2 * np.log(XS), XS**2, XS**3 * np.log(XS), XS**3, XS**4 * np.log(XS), XS**4])
        for prof in ("point", "ball", "shell"):
            y = np.array([rw[prof] for rw in rows])
            out[(n, prof)] = np.linalg.lstsq(A, y, rcond=None)[0]
    return out


@pytest.mark.parametrize("prof,zemach", [("point", 0.75), ("ball", 36 / 35), ("shell", 1.2)])
@pytest.mark.parametrize("n", [1, 2])
def test_leading_term_is_minus_two_zemach_for_every_n(fits, n, prof, zemach):
    assert abs(fits[(n, prof)][0] - (-2 * zemach)) < 2e-5


@pytest.mark.parametrize("prof", ["point", "ball", "shell"])
@pytest.mark.parametrize("n", [1, 2])
def test_log_term_is_minus_two_leading_defect(fits, n, prof):
    leading_defect = -2 / 5          # Paper 23 eq:defect_s, coefficient of x^2
    assert abs(fits[(n, prof)][1] - (-2 * leading_defect)) < 1e-3
