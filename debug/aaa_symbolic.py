"""Symbolic layer for debug/aaa_exact.py (p_pi pair PY and 3d-sigma pair DA on A, against 1s_B 1s_C).

Kernel decomposition (Fourier route).  For an A-pair density P(w) e^{-mu s} (P a polynomial, w = r - A, s = |w|,
mu = 2 zA) the transform is P(i grad_k) [8 pi mu / (k^2 + mu^2)^2].  For P invariant about an axis e it is a
polynomial in (k.e)^2 with coefficients R_j(k^2); the potential V-hat = 4 pi rho-hat / k^2 then reads

    V(w) = sum_j (e.grad_w)^{2j} U_{2j}(s),     U_{2j}-hat = (-1)^j 4 pi R_j(k^2) / k^2,

each U_{2j} radial, decaying, and in the family {1/s} + {s^m e^{-mu s}} (partial fractions in k^2: 4 pi/k^2 <-> 1/s,
4 pi/(k^2 + M)^{b+1} <-> (-1)^b/b! d^b/dM^b e^{-sqrt(M) s}/s).  The identity is exact as tempered distributions.
Family <-> engine names: 1/s = 'NA', s^m e^{-mu s} = (-1)^{m+1} d^{m+1}Y/dmu^{m+1}  (key m + 1).

Also the multipole potentials V_l(s) of the same densities by Gauss's formula (for the direct-quadrature checks).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache

import sympy as sp

zA = sp.symbols('zeta_A', positive=True)
s = sp.symbols('s', positive=True)
MU = 2 * zA


def density_poly(kind: str):
    """(normalization, polynomial in (wx, wy, wz)) with e = z axis.  kind 'p': (zA^5/pi) wz^2;
    'd': N_d^2 h2^2, h2 = (3 wz^2 - w^2)/2, N_d^2 = 5 (2 zA)^7 / (4 pi 720)."""
    wx, wy, wz = sp.symbols('wx wy wz', real=True)
    if kind == 'p':
        return zA ** 5 / sp.pi, wz ** 2, (wx, wy, wz)
    if kind == 'd':
        h2 = (3 * wz ** 2 - (wx ** 2 + wy ** 2 + wz ** 2)) / 2
        return 5 * (2 * zA) ** 7 / (4 * sp.pi * 720), sp.expand(h2 ** 2), (wx, wy, wz)
    raise ValueError(kind)


def _inv_ft(expr_K, K, M):
    """inverse Fourier transform of a radial function given as a rational function of K = k^2 with poles at
    K = 0 (simple) and K = -M; returns a sympy expression in s (mu = sqrt(M) substituted later)."""
    pf = sp.apart(sp.together(expr_K), K)
    out = 0
    for term in sp.Add.make_args(pf):
        num, den = sp.fraction(sp.factor(term))
        den = sp.factor(den)
        # identify den = c * K^a * (K + M)^b
        pd = sp.Poly(den, K)
        a = 0
        while pd.eval(0) == 0:
            pd = sp.Poly(sp.cancel(pd.as_expr() / K), K)
            a += 1
        b = sp.degree(pd, K)
        c = sp.cancel(pd.as_expr() / (K + M) ** b)
        assert not c.has(K), (term, c)
        assert not sp.Poly(num, K).degree() > 0, term
        coef = num / c
        if a == 1 and b == 0:
            out += coef / (4 * sp.pi * s)
        elif a == 0 and b >= 1:
            Mm = sp.symbols('Mm', positive=True)
            base = sp.exp(-sp.sqrt(Mm) * s) / s
            g = (-1) ** (b - 1) / sp.factorial(b - 1) * sp.diff(base, Mm, b - 1)
            out += coef / (4 * sp.pi) * g.subs(Mm, M)
        else:
            raise ValueError(f'unexpected term {term} (a={a}, b={b}): delta or growing piece')
    return out


@lru_cache(maxsize=None)
def kernels(kind: str):
    """{2j: U_{2j}(s)} sympy expressions in s, zA with V = sum_j (e.grad)^{2j} U_{2j}."""
    norm, poly, (wx, wy, wz) = density_poly(kind)
    kx, ky, kz = sp.symbols('kx ky kz', real=True)
    K, M = sp.symbols('K M', positive=True)
    fhat = 8 * sp.pi * MU / (kx ** 2 + ky ** 2 + kz ** 2 + MU ** 2) ** 2
    # P(i grad_k) fhat
    rho = 0
    for mon, co in sp.Poly(poly, wx, wy, wz).terms():
        d = fhat
        for var, pw in zip((kx, ky, kz), mon):
            if pw:
                d = sp.diff(d, var, pw)
        rho += co * sp.I ** sum(mon) * d
    rho = sp.simplify(norm * rho)
    rho = rho.subs(ky, 0)
    rho = sp.together(sp.expand(rho.subs(kx, sp.sqrt(K - kz ** 2))))
    rho = sp.together(rho.subs(MU ** 2, M)) if False else rho
    num, den = sp.fraction(sp.together(rho))
    den = sp.factor(den)
    assert not den.has(kz), den
    pn = sp.Poly(sp.expand(num), kz)
    out = {}
    for (pw,), co in pn.terms():
        assert pw % 2 == 0
        j = pw // 2
        Rj = sp.together(co / den)
        Uhat = (-1) ** j * 4 * sp.pi * Rj / K
        Uhat = sp.together(Uhat.subs(MU ** 2, M))   # write mu^2 as M (zA^2 = M/4)
        Uhat = sp.together(Uhat.subs(zA, sp.sqrt(M) / 2))
        U = _inv_ft(Uhat, K, M)
        U = sp.expand(sp.simplify(U.subs(M, MU ** 2)))
        out[2 * j] = U
    return out


def to_family(expr) -> dict:
    """sympy expression alpha/s + sum c_m s^m e^{-mu s} -> {'NA': alpha, m+1: (-1)^{m+1} c_m} (sympy in zA)."""
    e = sp.expand(sp.powsimp(sp.expand(expr)))
    out: dict = {}
    for term in sp.Add.make_args(e):
        ex = [f for f in sp.Mul.make_args(term) if f.func == sp.exp]
        rest = sp.Mul(*[f for f in sp.Mul.make_args(term) if f.func != sp.exp])
        if not ex:
            m = sp.degree(rest * s, s) - 1 if (rest * s).is_polynomial(s) else None
            co = sp.simplify(rest * s)
            assert not co.has(s), term
            out['NA'] = out.get('NA', 0) + co
            continue
        arg = sp.simplify(sum(f.args[0] for f in ex))
        assert sp.simplify(arg + MU * s) == 0, (term, arg)
        pw = sp.Poly(sp.expand(rest * s), s)   # rest = co * s^m, m >= -1
        for (q,), co in pw.terms():
            m = q - 1
            key = m + 1
            out[key] = out.get(key, 0) + (-1) ** (m + 1) * co
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def family_at(fam: dict, zval: Fraction) -> dict:
    """evaluate the family coefficients at a rational zA -> {key: Fraction}."""
    out = {}
    for k, v in fam.items():
        q = sp.nsimplify(v.subs(zA, sp.Rational(zval.numerator, zval.denominator)))
        assert q.is_Rational, (k, v, q)
        out[k] = Fraction(int(q.p), int(q.q))
    return out


@lru_cache(maxsize=None)
def multipole_potentials(kind: str):
    """{l: V_l(s)} for the A density written as sum_l f_l(s) P_l(u), u = (w.e)/s, by Gauss's formula
    V_l = 4 pi/(2l+1) [s^{-l-1} int_0^s t^{l+2} f_l + s^l int_s^oo t^{1-l} f_l]."""
    t = sp.symbols('t', positive=True)
    e = sp.exp(-MU * t)
    if kind == 'p':
        n = zA ** 5 / sp.pi
        fl = {0: n * t ** 2 * e / 3, 2: n * t ** 2 * e * sp.Rational(2, 3)}
    else:
        n = 5 * (2 * zA) ** 7 / (4 * sp.pi * 720)
        fl = {0: n * t ** 4 * e / 5, 2: n * t ** 4 * e * sp.Rational(2, 7), 4: n * t ** 4 * e * sp.Rational(18, 35)}
    V = {}
    for l, f in fl.items():
        inner = sp.integrate(t ** (l + 2) * f, (t, 0, s))
        outer = sp.integrate(t ** (1 - l) * f, (t, s, sp.oo))
        V[l] = sp.expand(sp.simplify(4 * sp.pi / (2 * l + 1) * (s ** (-l - 1) * inner + s ** l * outer)))
    return V, fl


def check_identity(kind: str, npts: int = 4) -> float:
    """max |sum_j d_z^{2j} U_{2j} - sum_l V_l P_l(u)| at random points (Cartesian derivatives, exact sympy, then
    evaluated at 40 digits)."""
    x, y, z = sp.symbols('x y z', real=True)
    Rr = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    Us = kernels(kind)
    V, _ = multipole_potentials(kind)
    lhs = 0
    for j2, U in Us.items():
        lhs += sp.diff(U.subs(s, Rr), z, j2)
    u = z / Rr
    rhs = sum(Vl.subs(s, Rr) * sp.legendre(l, u) for l, Vl in V.items())
    worst = 0.0
    import random
    random.seed(7)
    for _ in range(npts):
        pt = {x: sp.Rational(random.randint(-30, 30), 17), y: sp.Rational(random.randint(-30, 30), 19),
              z: sp.Rational(random.randint(-30, 30), 13), zA: sp.Rational(random.choice((9, 4, 7)), 10)}
        d = sp.N((lhs - rhs).subs(pt), 40)
        mag = sp.N(rhs.subs(pt), 40)
        worst = max(worst, abs(float(d / mag)))
    return worst


if __name__ == '__main__':
    import time
    for kind in ('p', 'd'):
        t0 = time.time()
        Us = kernels(kind)
        for j2, U in Us.items():
            print(kind, 'U', j2, '=', U)
            print('   family:', to_family(U))
        print(kind, 'identity residual (relative, 4 random points):', check_identity(kind), round(time.time() - t0, 1), 's')
