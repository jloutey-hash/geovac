"""Exact engine for the matched-exponent slice of the four-centre 1s integral (XY|ZW) (design pass for record FX;
computes no scored quantity).

Geometry as record FC (debug/fc_engine.py): X at the origin, Y = D1 z-hat, Z = a n-hat, W = b n-hat on one ray,
x = n-hat . z-hat; rho1 = e^{-zX|r-X|} e^{-zY|r-Y|}, rho2 = e^{-zZ|r-Z|} e^{-zW|r-W|} (unnormalised).

Why the matched condition.  With s = |r - Z|, t = |r - W|, mu = cos(angle(r, n)):
    s^2 = r^2 + a^2 - 2 r a mu,  t^2 = r^2 + b^2 - 2 r b mu   =>   b s^2 - a t^2 = (b - a)(r^2 - a b),
independent of mu.  At fixed r the pair (s, t) runs on the conic S^2 - T^2 = c(r), S = sqrt(b) s, T = sqrt(a) t.  The
exponent zZ s + zW t = (zZ/sqrt b) S + (zW/sqrt a) T is a light-cone coordinate S + T of that conic (so the angular
integral becomes a Laplace transform of a Laurent polynomial) iff zZ/sqrt(b) = zW/sqrt(a), i.e. zW/zZ = sqrt(a/b)
(the other light-cone direction needs zW < 0).  Off it, exponent = alpha (S+T) + beta c/(S+T): a Bessel-K kernel.
Rationality: kappa = sqrt(a/b) must be rational (a/b a ratio of squares) so that the endpoints are rational.

Closed form (kappa = sqrt(a/b), gamma(r) = s^2 - kappa^2 t^2 = (b - a)(r^2 - a b)/b, u = s + kappa t, gamma/u = s - kappa t):
    s = (u + gamma/u)/2,   s ds = (u - gamma^2/u^3)/4 du,   mu = (r^2 + a^2 - s^2)/(2 r a),
    rho2_L(r) = (2L+1)/(2 r a) int_{u_-}^{u_+} e^{-zZ u} P_L(mu(u)) (u - gamma^2 u^-3)/4 du,
    u_- = |r - a| + kappa |r - b|  (mu = 1),   u_+ = (r + a) + kappa (r + b)  (mu = -1).
The integrand is e^{-zZ u} times a Laurent polynomial with ODD powers u^-(2L+3) .. u^(2L+1), so
    rho2_L(r) = q+ e^{-zZ u+} + q- e^{-zZ u-} + R_L(r) [E1(zZ u-) - E1(zZ u+)],
    R_L(r) = (2L+1)/(2 r a) Res_{u=0} e^{-zZ u} P_L(mu(u)) (u - gamma^2 u^-3)/4      (nonzero: L = 0 gives
    -(2L+1)/(2ra) zZ^2 gamma^2 / 8).
So rho2_L is exp + E1 in r on the slice; the potential multipole psi_L of rho1 carries E1((zX + zY) r) (record NA's
constant set).  F_L = int r^2 rho2_L psi_L dr therefore contains integrals of E1 x E1 x exp x rational, outside the
{e, E1, ln, gamma} class that na_exact/aabc_exact reduce to.  See debug/fx_design_memo.md.

API:
    MATCHED: example cases;  rho2_split(case, r, L) -> {('e', u): q, ('E1', u): q} exact Fractions (rational r);
    rho2_value(split, mp);  residue_poly(case, L) -> sympy expression R_L(r) (symbolic r).
"""
from __future__ import annotations

from fractions import Fraction

MATCHED = {
    # kappa = 1/2: a = 1, b = 4, zW = zZ/2.  X, Y as FC's G2.
    'M1': dict(a='1', b='4', D1='6/5', zX='13/10', zY='1', zZ='6/5', zW='3/5'),
    # kappa = 2/3: a = 4/5, b = 9/5, zW = 2 zZ/3; D1 between a and b.
    'M2': dict(a='4/5', b='9/5', D1='6/5', zX='13/10', zY='1', zZ='3/2', zW='1'),
}


def params(case):
    d = MATCHED[case] if isinstance(case, str) else case
    p = {k: Fraction(v) for k, v in d.items()}
    k2 = p['a'] / p['b']
    kap = _qsqrt(k2)
    if kap is None or p['zW'] != kap * p['zZ']:
        raise ValueError('not a rational matched slice: need sqrt(a/b) rational and zW = sqrt(a/b) zZ')
    p['kappa'] = kap
    return p


def _qsqrt(q):
    q = Fraction(q)
    from math import isqrt
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    return Fraction(n, d) if n * n == q.numerator and d * d == q.denominator else None


def _poly_mul(A, B):
    out = {}
    for i, x in A.items():
        for j, y in B.items():
            out[i + j] = out.get(i + j, 0) + x * y
    return {k: v for k, v in out.items() if v != 0}


def laurent_h(p, r, L):
    """{n: coef} of h(u) = P_L(mu(u)) (u - gamma^2 u^-3)/4 at rational r (exact)."""
    a, b = p['a'], p['b']
    g = (b - a) * (r * r - a * b) / b
    # s^2 = (u^2 + 2 g + g^2 u^-2)/4 ; mu = (r^2 + a^2 - s^2)/(2 r a)
    c = 1 / (2 * r * a)
    mu = {2: -c / 4, 0: c * (r * r + a * a - g / 2), -2: -c * g * g / 4}
    mu = {k: v for k, v in mu.items() if v != 0}
    P0, P1 = {0: Fraction(1)}, mu
    if L == 0:
        PL = P0
    else:
        for l in range(1, L):
            t = _poly_mul(mu, P1)
            nxt = {k: Fraction(2 * l + 1, l + 1) * v for k, v in t.items()}
            for k, v in P0.items():
                nxt[k] = nxt.get(k, 0) - Fraction(l, l + 1) * v
            P0, P1 = P1, {k: v for k, v in nxt.items() if v != 0}
        PL = P1
    return _poly_mul(PL, {1: Fraction(1, 4), -3: -g * g / 4})


def rho2_split(case, r, L):
    """Exact split of rho2_L(r) at rational r: {('e', u): q} for q e^{-zZ u}, {('E1', u): q} for q E1(zZ u)."""
    p = params(case)
    r = Fraction(r)
    a, b, kap, z = p['a'], p['b'], p['kappa'], p['zZ']
    um = abs(r - a) + kap * abs(r - b)
    up = (r + a) + kap * (r + b)
    h = laurent_h(p, r, L)
    pref = Fraction(2 * L + 1) / (2 * r * a)
    out = {}

    def put(k, v):
        if v != 0:
            out[k] = out.get(k, 0) + v

    for n, c in h.items():
        c = c * pref
        if n >= 0:
            # int e^{-zu} u^n = G(u2) - G(u1), G = -e^{-zu} sum_k n!/(n-k)! u^{n-k} / z^{k+1}
            for (u, sgn) in ((up, 1), (um, -1)):
                fac, s = Fraction(1), Fraction(0)
                for k in range(n + 1):
                    s += fac * u ** (n - k) / z ** (k + 1)
                    fac *= (n - k)
                put(('e', u), -sgn * c * s)
        else:
            m = -n
            # int_{u1}^{u2} e^{-zu} u^-m = u1^{1-m} E_m(z u1) - u2^{1-m} E_m(z u2)
            for (u, sgn) in ((um, 1), (up, -1)):
                x = z * u
                al, be = Fraction(1), Fraction(0)          # E_1 = 1*E1 + 0*e^{-x}
                for j in range(1, m):
                    al, be = -x * al / j, (1 - x * be) / j
                w = sgn * c * u ** (1 - m)
                put(('E1', u), w * al)
                put(('e', u), w * be)
    return {k: v for k, v in out.items() if v != 0}


def rho2_value(split, case, mp):
    z = params(case)['zZ']
    zz = mp.mpf(z.numerator) / z.denominator
    s = mp.mpf(0)
    for (kind, u), q in split.items():
        uu = zz * mp.mpf(u.numerator) / u.denominator
        s += (mp.mpf(q.numerator) / q.denominator) * (mp.exp(-uu) if kind == 'e' else mp.e1(uu))
    return s


def residue_poly(case, L):
    """R_L(r), the E1(zZ u-) coefficient of rho2_L(r), as a sympy rational function of r (exact)."""
    import sympy as sp
    p = params(case)
    r, u = sp.symbols('r u', positive=True)
    a, b, z = (sp.Rational(p[k].numerator, p[k].denominator) for k in ('a', 'b', 'zZ'))
    g = (b - a) * (r ** 2 - a * b) / b
    s2 = (u ** 2 + 2 * g + g ** 2 / u ** 2) / 4
    mu = (r ** 2 + a ** 2 - s2) / (2 * r * a)
    h = sp.expand(sp.legendre(L, mu) * (u - g ** 2 / u ** 3) / 4)
    nmin = 2 * L + 3
    ser = sp.series(sp.exp(-z * u), u, 0, nmin + 1).removeO()
    res = sp.expand(h * ser).coeff(u, -1)
    return sp.factor(sp.Rational(2 * L + 1) / (2 * r * a) * res), r
