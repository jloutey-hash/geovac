"""Guard for Paper 59's derived bond-angle singular term (records docs/predictions/2026-10-27_water_cusp_singularity.md,
docs/predictions/2026-10-27_water_cusp_next_term.md; WA6 GO).

Claims guarded, as the paper states them:
  * the local cusp-cusp computation reproduces the two-centre 1s Coulomb integral's first odd terms,
    (2/315) R^7 + (16/14175) R^9, exactly, with no R, R^3, R^5 term;
  * at water's equal bond lengths its first two orders alpha7 (1-x)^(7/2) + alpha9 (1-x)^(9/2) match the
    channel-sum Legendre coefficients with nothing fitted, the remainder falling like l^-4.
The derivation is copied from debug/wa6_local.py; the channel formula from tests/test_paper59_water_bond_angle.py.
"""
import mpmath as mp
import numpy as np
import pytest

from test_paper59_water_bond_angle import F_coeffs, WATER

import itertools

import sympy as sp

kx, ky, kz, Q = sp.symbols('kx ky kz Q', positive=True)
wx, wy, wz, Rw = sp.symbols('wx wy wz Rw', positive=True)
K = (kx, ky, kz)
W = (wx, wy, wz)
ax, ay, az, tt = sp.symbols('ax ay az tt')
A_ = (ax, ay, az)


def dk(f, j):
    return sp.diff(f, K[j]) + 2 * K[j] * sp.diff(f, Q)


def dw(f, j):
    return sp.diff(f, W[j]) + W[j] / Rw * sp.diff(f, Rw)


def taylor_poly(expr, deg):
    """Taylor polynomial in (ax, ay, az) to total degree deg."""
    e = sp.series(expr.subs({ax: tt * ax, ay: tt * ay, az: tt * az}), tt, 0, deg + 1).removeO()
    return sp.expand(e.subs(tt, 1))


def ft_cusp_poly(P, sign):
    """FT of |a| P(a), evaluated at sign*k, as a polynomial in k over powers of Q."""
    out = 0
    base = -8 * sp.pi * Q ** -2
    for term in sp.Add.make_args(sp.expand(P)):
        c, mon = term.as_coeff_Mul()
        powers = [sp.degree(mon, v) if mon.has(v) else 0 for v in A_]
        coeff = sp.simplify(term / (A_[0] ** powers[0] * A_[1] ** powers[1] * A_[2] ** powers[2]))
        f = base
        for j, pw in enumerate(powers):
            for _ in range(pw):
                f = sp.I * dk(f, j)
        f = f.subs({kx: sign * kx, ky: sign * ky, kz: sign * kz}, simultaneous=True)
        out += coeff * f
    return sp.expand(out)


def inv_ft(Shat):
    """Inverse FT of a sum of k-monomials times Q^-m."""
    out = 0
    for term in sp.Add.make_args(sp.expand(Shat)):
        m = sp.degree(term.subs(Q, 1 / sp.Symbol('q_')), sp.Symbol('q_'))  # power of Q^-1
        mon = sp.simplify(term * Q ** m)
        powers = [sp.degree(mon, v) if mon.has(v) else 0 for v in K]
        coeff = sp.simplify(mon / (kx ** powers[0] * ky ** powers[1] * kz ** powers[2]))
        g = sp.gamma(sp.Rational(3, 2) - m) / (2 ** (2 * m) * sp.pi ** sp.Rational(3, 2) * sp.gamma(m)) * Rw ** (2 * m - 3)
        for j, pw in enumerate(powers):
            for _ in range(pw):
                g = -sp.I * dw(g, j)
        out += coeff * g
    return sp.expand(out)


def singular_part(G1, G2, maxdeg=2):
    """S(w) from cusp polynomials: G1 Taylor in a about Y, G2 Taylor in b about Z (both in A_ variables)."""
    S = 0
    P1 = G1 if G1.is_polynomial(*A_) else taylor_poly(G1, maxdeg)
    P2 = G2 if G2.is_polynomial(*A_) else taylor_poly(G2, maxdeg)
    def by_degree(P):
        d = {}
        for term in sp.Add.make_args(P):
            deg = sum(sp.degree(term, v) if term.has(v) else 0 for v in A_)
            d[deg] = d.get(deg, 0) + term
        return d
    D1, D2 = by_degree(P1), by_degree(P2)
    for d1, d2 in itertools.product(D1, D2):
        if d1 + d2 > maxdeg:
            continue
        # f1~(a) = f1(-a): FT at -k
        Shat = ft_cusp_poly(D1[d1], -1) * ft_cusp_poly(D2[d2], +1) * 4 * sp.pi / Q
        S += inv_ft(Shat)
    return sp.expand(S)


def water_symmetric(zX, zH, D):
    """alpha7, alpha9 with S = alpha7 (1-x)^(7/2) + alpha9 (1-x)^(9/2) + ..., D1 = D2 = D."""
    zX, zH, D = sp.nsimplify(zX), sp.nsimplify(zH), sp.nsimplify(D)
    NX, NH = sp.sqrt(zX ** 3 / sp.pi), sp.sqrt(zH ** 3 / sp.pi)
    u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
    c, s = sp.symbols('c s')
    Y = sp.Matrix([0, 0, D])
    Z = sp.Matrix([D * s, 0, D * c])
    rr = sp.sqrt(ax ** 2 + ay ** 2 + (D + az) ** 2)
    GY = -NH * NX * sp.exp(-zX * rr) * sp.sinh(zH * u) / u
    PY = taylor_poly(GY, 2)
    # about Z = Rot Y, Rot about the y axis taking z to (s, 0, c): G_Z(b) = G_Y(Rot^T b)
    bx, bz = c * ax - s * az, s * ax + c * az
    PZ = sp.expand(PY.subs({ax: bx, az: bz}, simultaneous=True))
    S = singular_part(PY, PZ)
    w = Y - Z
    S = S.subs({wx: w[0], wy: w[1], wz: w[2]})
    S = sp.expand(S)
    S = S.subs(s ** 2, 1 - c ** 2)
    S = sp.expand(S)
    assert not S.has(s), "odd powers of sin(gamma) left"
    v = sp.symbols('v', positive=True)               # v^2 = 1 - c
    R = D * sp.sqrt(2) * v
    S = sp.expand(S.subs(c, 1 - v ** 2).subs(Rw, R))
    ser = sp.expand(S)
    a7 = ser.coeff(v, 7)
    a9 = ser.coeff(v, 9)
    lower = [ser.coeff(v, k) for k in (1, 3, 5)]
    return sp.simplify(a7), sp.simplify(a9), lower




def _leg(a, l):
    mp.mp.dps = 30
    return float((2 * l + 1) * mp.mpf(2) ** a * mp.gamma(a + 1) ** 2
                 / (mp.gamma(a + l + 2) * mp.gamma(a + 1 - l)) * (-1) ** l)


def test_two_centre_odd_terms():
    u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
    G = -(1 / sp.pi) * sp.sinh(2 * u) / u
    R = sp.symbols('R', positive=True)
    val = sp.expand(sp.simplify(singular_part(G, G).subs({wx: 0, wy: 0, wz: R, Rw: R})))
    assert sp.simplify(val - (sp.Rational(2, 315) * R ** 7 + sp.Rational(16, 14175) * R ** 9)) == 0
    Jx = 1 / R - sp.exp(-2 * R) * (1 / R + sp.Rational(11, 8) + sp.Rational(3, 4) * R + R ** 2 / 6)
    ser = sp.series(Jx, R, 0, 10).removeO()
    assert [ser.coeff(R, k) for k in (1, 3, 5)] == [0, 0, 0]
    assert ser.coeff(R, 7) == sp.Rational(2, 315) and ser.coeff(R, 9) == sp.Rational(16, 14175)


@pytest.mark.slow
def test_symmetric_water_series_from_derived_terms():
    for zX, zH in ((1, 1), ('2.25', '1.24')):
        a7, a9, lower = water_symmetric(zX, zH, str(WATER))
        assert lower == [0, 0, 0]
        a7, a9 = float(a7), float(a9)
        zXf, zHf = float(sp.nsimplify(zX)), float(sp.nsimplify(zH))
        assert abs(a7 - zHf ** 5 * zXf ** 3 * np.exp(-2 * zXf * WATER) / 630 * (2 * WATER ** 2) ** 3.5) < 1e-14
        F = F_coeffs(WATER, WATER, zXf, zHf, 40, level=1)
        dev2 = {l: abs(F[l] / (a7 * _leg(3.5, l) + a9 * _leg(4.5, l)) - 1) for l in (30, 40)}
        dev0 = abs(F[40] / (a7 * _leg(3.5, 40)) - 1)
        assert dev2[40] < 1.5e-3, (zX, dev2)
        assert 2.5 < dev2[30] / dev2[40] < 4.5, (zX, dev2)       # l^-4 gives (4/3)^4 = 3.16
        assert dev0 > 10 * dev2[40], (zX, dev0, dev2)


# Contact form and the zero-contact case (record docs/predictions/2026-10-27_water_zero_contact.md, WA7; the
# zero-contact numbers are reported only in the paper).  Derivation copied from debug/wa7_local.py.

FREE = [(1.0, 1.0), (2.0, -2.0 ** -2.5)]     # phi_1 - 2^(-5/2) phi_2: first-order cusp cancels


def _F_combo(zX, H1, H2, lmax, level):
    from numpy.polynomial.legendre import leggauss
    from scipy.special import spherical_jn
    from test_paper59_water_bond_angle import _G, _tanh_sinh
    t, wt = _tanh_sinh(60 * 2 ** level, 0.07 / 2 ** level)
    u, wu = leggauss(400 * 2 ** level)
    u, wu = 0.5 * (u + 1), 0.5 * wu
    k = 3.0 * u / (1 - u)
    wk = wu * 3.0 / (1 - u) ** 2
    G1 = sum(c * _G(k[:, None], t[None, :], WATER, zX, z) for z, c in H1)
    G2 = sum(c * _G(k[:, None], t[None, :], WATER, zX, z) for z, c in H2)
    out = []
    for l in range(lmax + 1):
        J = spherical_jn(l, np.outer(k, t) * WATER)
        out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * ((J * G1) @ wt) * ((J * G2) @ wt)))
    return np.array(out)


def _free_G_poly(zX, D, deg):
    NX = sp.sqrt(zX ** 3 / sp.pi)
    u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
    rr = sp.sqrt(ax ** 2 + ay ** 2 + (D + az) ** 2)
    odd = 0
    for z, c in FREE:
        z, c = sp.nsimplify(z), (1 if c == 1.0 else -sp.Rational(1, 2) ** sp.Rational(5, 2))
        odd += c * sp.sqrt(z ** 3 / sp.pi) * sp.sinh(z * u)
    return taylor_poly(-NX * sp.exp(-zX * rr) * odd / u, deg)


def test_contact_factorization():
    """Unit contacts give pi^2 |w|^7 / 630; 1s contacts g = -zeta_H chi_X chi_Y at Y rebuild C."""
    S = singular_part(sp.Integer(1), sp.Integer(1), maxdeg=0)
    assert sp.simplify(S - sp.pi ** 2 * Rw ** 7 / 630) == 0
    zX, zH, D = sp.symbols('zX zH D', positive=True)
    g = -zH * sp.sqrt(zH ** 3 / sp.pi) * sp.sqrt(zX ** 3 / sp.pi) * sp.exp(-zX * D)
    assert sp.simplify(sp.pi ** 2 / 630 * g * g - zH ** 5 * zX ** 3 * sp.exp(-2 * zX * D) / 630) == 0
    assert sp.simplify(_free_G_poly(sp.Integer(1), sp.nsimplify(WATER), 1)) == 0    # zero contact


@pytest.mark.slow
def test_zero_contact_drops_the_l8_term():
    zX, D = sp.Integer(1), sp.nsimplify(str(WATER))
    c, s = sp.symbols('c s')
    PY = _free_G_poly(zX, D, 4)
    PZ = sp.expand(PY.subs({ax: c * ax - s * az, az: s * ax + c * az}, simultaneous=True))
    S = singular_part(PY, PZ, maxdeg=6)
    w = sp.Matrix([0, 0, D]) - sp.Matrix([D * s, 0, D * c])
    S = sp.expand(sp.expand(S.subs({wx: w[0], wy: w[1], wz: w[2]})).subs(s ** 2, 1 - c ** 2))
    v = sp.symbols('v', positive=True)
    S = sp.expand(S.subs(c, 1 - v ** 2).subs(Rw, D * sp.sqrt(2) * v))
    assert [S.coeff(v, k) for k in (1, 3, 5, 7, 9)] == [0, 0, 0, 0, 0]
    b11, b13 = float(S.coeff(v, 11)), float(S.coeff(v, 13))
    F = _F_combo(1.0, FREE, FREE, 30, level=1)
    L1 = b11 * _leg(5.5, 30)
    L2 = L1 + b13 * _leg(6.5, 30)
    dev2, dev1 = abs(F[30] / L2 - 1), abs(F[30] / L1 - 1)
    assert dev2 < 2e-2, dev2
    assert dev1 > 5 * dev2, (dev1, dev2)
    assert F[30] > 0
