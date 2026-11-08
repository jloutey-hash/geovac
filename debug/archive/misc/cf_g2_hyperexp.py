"""CF5-CF7: exact hyperexponential solutions of an operator piece; hyperexp_exact and order1_from_s copied
verbatim from debug/correlated_triangle_g32_7_refine.py (that script runs on import, so it is not imported)."""
from __future__ import annotations
from fractions import Fraction
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib4 as S


def hyperexp_exact(L, f0, combo, N, maxprimes=200):
    """exact P (in t) for s = t^f0 prod (1 - t/a)^{e_a} P(t), via mod p nullspaces + CRT"""
    got, last = [], None
    for p in H.PRIMES[:maxprimes]:
        NT = 2 * (N + 1) + 60 + 2 * G.op_degree(L)
        pre = S.prefactor(combo, NT, p)
        Lt = S.t_series_ops(L, p)
        ep = G.frac_mod(Fraction(f0), p)
        cols = []
        for j in range(N + 1):
            s = [0] * j + pre[: NT - j]
            cols.append(S.apply_t(Lt, ep, s, p, NT))
        valid = max(NT - 2 * G.op_degree(L) - 8, N + 31)
        rows = [[col[n] for col in cols] for n in range(min(valid, NT))]
        vecs = G.nullspace_mod(rows, N + 1, p)
        if len(vecs) != 1:
            continue
        v = vecs[0]
        idx = max(i for i, x in enumerate(v) if x)
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % 4:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is not None and fr == last:
            return fr
        last = fr
    return None


def order1_from_s(f0, combo, Pt):
    """rho = Dn theta - Dn (theta s / s), as a polynomial operator in x (s even/odd in t handled: P must be even)."""
    t, x = sp.symbols("t x")
    Ps = sum(sp.Rational(c.numerator, c.denominator) * t**i for i, c in enumerate(Pt))
    assert sp.expand(Ps - Ps.subs(t, -t)) == 0, "P not even in t"
    s_log = sp.Integer(f0) * sp.log(t) + sum(sp.Rational(e.numerator, e.denominator) * sp.log(1 - t / sp.Rational(a.numerator, a.denominator))
                                             for a, e in combo.items() if e != 0) + sp.log(Ps)
    r = sp.together(sp.expand(t * sp.diff(s_log, t)))
    num, den = sp.fraction(sp.cancel(r))
    num, den = sp.expand(num.subs(t, sp.sqrt(x))), sp.expand(den.subs(t, sp.sqrt(x)))
    pn = sp.Poly(num, x).all_coeffs()[::-1]
    pd = sp.Poly(den, x).all_coeffs()[::-1]
    lc = sp.ilcm(*[sp.Rational(c).q for c in pn + pd])
    op = [flint.fmpz_poly([-int(c * lc) for c in pn]), flint.fmpz_poly([int(c * lc) for c in pd])]
    return G.content_free(op), sp.factor(Ps.subs(t, sp.sqrt(x)))


