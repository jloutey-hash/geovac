"""WA6 derivation: the non-analytic part of (XY|XZ) near the coalescence Y = Z, through degree 9 in w = Y - Z.

rho1 = chi_X chi_Y.  Its non-smooth part at Y is |a| G1(a), a = r - Y, with
G1(a) = -N_H N_X exp(-zX |Y + a|) sinh(zH |a|)/|a| (analytic in a).  Likewise at Z.
Singular part S(w) = FT^-1[ F1hat(-k) F2hat(k) 4 pi / k^2 ], F hat[|a| a^alpha] = (i d_k)^alpha (-8 pi / k^4),
FT^-1[k^-2m] = Gamma(3/2 - m) / (2^2m pi^3/2 Gamma(m)) |w|^(2m-3), k^alpha -> (-i d_w)^alpha.
Only cusp-cusp pairs depend on w non-analytically; degrees dA + dB <= 2 fix the (1-x)^(7/2) and (1-x)^(9/2)
coefficients at D1 = D2 (record docs/predictions/2026-10-27_water_cusp_next_term.md).
Usage: python debug/wa6_local.py   (prints alpha7, alpha9 for E1 and E2; checks the two-centre R^7, R^9 terms)
"""
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


def two_centre_check():
    """rho = e^-2r/pi at both centres (no third centre): the R^7 and R^9 coefficients of J(R)."""
    u = sp.sqrt(ax ** 2 + ay ** 2 + az ** 2)
    G = -(1 / sp.pi) * sp.sinh(2 * u) / u
    S = singular_part(G, G)
    R = sp.symbols('R', positive=True)
    val = sp.simplify(S.subs({wx: 0, wy: 0, wz: R, Rw: R}))
    print("two-centre singular part:", sp.expand(val), " expected 2R^7/315 + 16R^9/14175")


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


if __name__ == "__main__":
    two_centre_check()
    for name, (zX, zH) in {'E1': (1, 1), 'E2': ('2.25', '1.24')}.items():
        a7, a9, lower = water_symmetric(zX, zH, '1.809')
        print(name, "alpha7 =", sp.N(a7, 15), " alpha9 =", sp.N(a9, 15), " alpha9/alpha7 =", sp.N(a9 / a7, 12),
              " lower odd terms:", [sp.N(v, 5) for v in lower])
        C = sp.nsimplify(zH) ** 5 * sp.nsimplify(zX) ** 3 * sp.exp(-sp.nsimplify(zX) * 2 * sp.nsimplify('1.809')) / 630
        print("   leading check C (2 D^2)^(7/2):", sp.N(C * (2 * sp.nsimplify('1.809') ** 2) ** sp.Rational(7, 2), 15))
