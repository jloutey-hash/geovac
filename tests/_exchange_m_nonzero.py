"""Helper for tests/test_paper58_exchange_m_nonzero.py -- the m != 0 two-centre exchange integrals through the
r12-weighted masters f_n of Pachucki (arXiv:0905.1635), derived from his label 39 without the 1/p^2 of labels 45/46.

Equation (derived from label 39, first order in the r12 exponent w1, by the inverse Laplace transform t -> d/dr):

    q f_{m+1}'' = c f_{m+1} + m^2 (d^4 - Sig d^2 + e) f_{m-1} + m (m-1)^2 (m-2) f_{m-3}'' + S_m(r)          (**)

    q = (u2^2-u3^2)(w2^2-w3^2),  c = (u2^2-u3^2+w2^2-w3^2)(u2^2 w2^2-u3^2 w3^2),  Sig = u2^2+u3^2+w2^2+w3^2,
    e = (u2^2-w3^2)(u3^2-w2^2),  S_m = (-1)^{m+1} L^{-1}[d^m P/dw1^m at w1 = 0],  P = label 13 with u1 = t.

Derivatives in the four exponents (u3, u2, w2, w3) at the exchange point u3 = w2 = a, u2 = w3 = b are taken by
Leibniz, lowest multi-index first; each step is u'' = known, solved in the decaying span (the homogeneous solutions
A + B r do not decay, so the solution is unique; the solver REPORTS how many unknowns it left free).

Densities (per electron N_A N_B P e^{-a r_A - b r_B}; X = x1 x2 + y1 y2, rho^2 = x^2 + y^2):
    2ppi  = (2p+1_A 1s_B | 2p-1_A 1s_B)  = -(a^5 b^3/pi^2)(1/2) int int X psi psi / r12
    3ddel = (3d+2_A 1s_B | 3d-2_A 1s_B)  =  (a^7 b^3/(12 pi^2)) int int (2 X^2 - rho1^2 rho2^2) psi psi / r12
with X = (r1A^2 + r2A^2 - r12^2)/2 - z1 z2 and rho^2 = r_A^2 - z^2.  A monomial r12^{2k} prod r^{n_i} / r12 maps to
(4 pi)^2 / R d^alpha f_{2k}, alpha = n + 1.

Copied from the derivation's scripts (debug/pachucki_m_core.py, debug/pachucki_m_cases.py, frozen 2026-11-18) onto
the reviewed engine tests/_pachucki_extension.py (imported as X).  One change: `solve_fn` returns the number of
unknowns each `solve_D2` left free (the derivation set them to 0 silently).  The python-flint fallback for sympy's
heuristic gcd is copied too and installed only by `install_gcd_fallback()` (Q(t) runs).
Nothing here reads or imports debug/.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from typing import Dict, List, Sequence, Tuple

import sympy as sp

import _pachucki_extension as X

Func = X.Func
fadd, fscale, clean = X.fadd, X.fscale, X.clean
CASES_M = ("2ppi", "3ddel")

# sigma = 2 part of real 2p_x 2p_x exchange: (xx|xx) = pppi0 + ppdel / 2, ppdel = SIGMA2_RATIO(a, b) * 3ddel
SIGMA2_RATIO = lambda a, b: 3 * (b / a) ** 2          # noqa: E731

# ------------------------------------------------------------------ sympy heuristic gcd: exact fallback
_GCD_INSTALLED = [False]


def _gcd_ZZ_safe(f, g):
    from sympy.polys.heuristicgcd import heugcd
    from sympy.polys.polyerrors import HeuristicGCDFailed
    try:
        return heugcd(f, g)
    except HeuristicGCDFailed:
        return _gcd_ZZ_flint(f, g)


def _gcd_ZZ_flint(f, g):
    """(h, f/h, g/h) for univariate integer PolyElements, by python-flint's exact fmpz_poly gcd."""
    import flint
    ring = f.ring
    assert ring.ngens == 1
    fa = flint.fmpz_poly([int(c) for c in reversed(f.to_dense())])
    ga = flint.fmpz_poly([int(c) for c in reversed(g.to_dense())])
    h = fa.gcd(ga)
    cf, cg = fa // h, ga // h
    assert cf * h == fa and cg * h == ga

    def back(q):
        return ring.from_dense([ring.domain(int(c)) for c in reversed(q.coeffs())])
    return back(h), back(cf), back(cg)


def install_gcd_fallback() -> None:
    """sympy 1.14's PolyElement._gcd_ZZ is heugcd with no fallback and raises HeuristicGCDFailed on large univariate
    polynomials met in the Q(t) runs; route that failure to the exact flint gcd."""
    import sympy.polys.rings as rings
    if not _GCD_INSTALLED[0]:
        rings.PolyElement._gcd_ZZ = _gcd_ZZ_safe
        _GCD_INSTALLED[0] = True


# ------------------------------------------------------------------ label 13 (P with u1 = t)
u2s, u3s, w2s, w3s, w1s, ts = sp.symbols("u2 u3 w2 w3 w1 t")
VARIDX = {u3s: 0, u2s: 1, w2s: 2, w3s: 3}


def P_logs():
    """t-dependent logarithms of P(w1, t; w2, u2; w3, u3): list of (C(w1, t, theta), c(theta), eps), the term being
    C ln(t + c + eps w1).  The remaining terms of label 13 are polynomial in t (no r > 0 content)."""
    u2, u3, w2, w3, w1, t = u2s, u3s, w2s, w3s, w1s, ts
    A3 = t**2 * w1**2 + u2**2 * w2**2 - u3**2 * w3**2 + w1 * w2 * (t**2 + u2**2 - w3**2)
    A4 = t**2 * w1**2 - u2**2 * w2**2 + u3**2 * w3**2 + w1 * w3 * (t**2 + u3**2 - w2**2)
    A5 = u2 * (u2 + w1) * (t**2 + u3**2 - w2**2) - u3**2 * (t**2 + u2**2 - w3**2)
    A6 = u3 * (u3 + w1) * (t**2 + u2**2 - w3**2) - u2**2 * (t**2 + u3**2 - w2**2)
    D3 = w3**2 - (w1 + w2)**2
    D4 = w2**2 - (w1 + w3)**2
    D5 = u3**2 - (u2 + w1)**2
    D6 = u2**2 - (u3 + w1)**2
    return [
        (-A3 / D3 + A6 / D6, u2 + w3, 0),
        (-A4 / D4 + A5 / D5, u3 + w2, 0),
        (t * w1 + A3 / D3 - A5 / D5, u2 + w2, 1),
        (t * w1 + A4 / D4 - A6 / D6, u3 + w3, 1),
    ]


def P_full():
    """P(w1, u1 = t; w2, u2; w3, u3), label 13 verbatim (checks only)."""
    u2, u3, w2, w3, w1, u1 = u2s, u3s, w2s, w3s, w1s, ts

    def G(a1, a2, a3):
        return (sp.log(a1 + a2) - sp.log(a1 + a3)) / ((a2 - a3) * (a2 + a3))
    return (u1 * w1 * ((u1 + w2)**2 - u3**2) * G(u2 + w1, u3, u1 + w2)
            + u1 * w1 * ((u1 + u3)**2 - w2**2) * G(w1 + w3, w2, u1 + u3)
            - (u1**2 * w1**2 + u2**2 * w2**2 - u3**2 * w3**2 + w1 * w2 * (u1**2 + u2**2 - w3**2)) * G(u1 + u2, w3, w1 + w2)
            - (u1**2 * w1**2 - u2**2 * w2**2 + u3**2 * w3**2 + w1 * w3 * (u1**2 + u3**2 - w2**2)) * G(u1 + u3, w2, w1 + w3)
            + (u2 * (u2 + w1) * (u1**2 + u3**2 - w2**2) - u3**2 * (u1**2 + u2**2 - w3**2)) * G(u1 + w2, u3, u2 + w1)
            + (u3 * (u3 + w1) * (u1**2 + u2**2 - w3**2) - u2**2 * (u1**2 + u3**2 - w2**2)) * G(u1 + w3, u2, u3 + w1)
            - w1 * (w2 * (u1**2 - u2**2 + w3**2) + w3 * (u1**2 + w2**2 - u3**2)) * G(u2 + u3, w1, w2 + w3)
            - w1 * (u2 * (u1**2 - w2**2 + u3**2) + u3 * (u1**2 + u2**2 - w3**2)) * G(w2 + w3, w1, u2 + u3))


def jet_eval(J: X.Jets, expr, base):
    """A rational sympy expression in (u3, u2, w2, w3) as a scalar jet; base = the four variable jets."""
    ctx = J.ctx
    if expr.is_Rational:
        return J.const(ctx.c(Fraction(int(expr.p), int(expr.q)))) if expr != 0 else {}
    if expr.is_Symbol:
        return base[VARIDX[expr]]
    if expr.is_Add:
        out = {}
        for a in expr.args:
            out = J.add(out, jet_eval(J, a, base))
        return out
    if expr.is_Mul:
        out = J.const(ctx.one)
        for a in expr.args:
            out = J.mul(out, jet_eval(J, a, base))
        return out
    if expr.is_Pow:
        n = int(expr.exp)
        b = jet_eval(J, expr.base, base)
        if n < 0:
            b, n = J.inv(b), -n
        out = J.const(ctx.one)
        for _ in range(n):
            out = J.mul(out, b)
        return out
    raise ValueError(expr)


def _fact(al) -> int:
    f = 1
    for x in al:
        f *= factorial(x)
    return f


# ------------------------------------------------------------------ the solver of u'' = src
class MCalc(X.Calc):
    """X.Calc plus the solver of u'' = src in the decaying span (records self.last_free like solve_L0)."""

    def D2(self, A: Func) -> Func:
        return self.diff(self.diff(A))

    def solve_D2(self, src: Func) -> Func:
        ctx = self.ctx
        zero = ctx.zero
        rates = set()
        for k in list(src) + list(self.FM1) + list(self.F0):
            if k[0] == "E" and k[1] != 0:
                rates.add(k[1])
            if k[0] == "Ei":
                rates.add(k[1] + k[2])
                if k[1] == 0:
                    rates.add(k[2])
        mus = sorted([m for m in rates if m < 0] if ctx.kind == "Q" else [m for m in rates if self._neg(m)], key=str)
        n = 2 + 2 * len(mus)
        Z = [zero] * (n + 1)

        def unit(i):
            u = [zero] * (n + 1)
            u[i] = ctx.one
            return u

        def aadd(x, y, s=1):
            return [xi + s * yi for xi, yi in zip(x, y)]

        def ascale(x, c):
            return [xi * c for xi in x]

        target = {k: {p: [c] + [zero] * n for p, c in d.items()} for k, d in src.items()}
        sol, conds = {}, []

        def sub_image(img: Func, coeff):
            for k, d in img.items():
                o = target.setdefault(k, {})
                for p, c in d.items():
                    o[p] = aadd(o.get(p, Z), ascale(coeff, c), -1)

        # 1. f key:  (r^p f0)'' = p(p-1) r^{p-2} f0 + (1-2p) r^{p-2} FM1 - r^{p-1} F0
        tf = target.pop(("f",), {})
        xf = {}
        for q_, v in tf.items():
            p = q_ + 2
            if p in (0, 1):
                conds.append(v)
                continue
            xf[p] = ascale(v, 1 / ctx.c(p * (p - 1)))
        xf[0], xf[1] = unit(1), unit(2)
        for p, x in xf.items():
            sub_image(fadd(fscale(self.FM1, ctx.c(1 - 2 * p), p - 2), fscale(self.F0, ctx.c(-1), p - 1)), x)
        sol[("f",)] = xf

        def solve_key_rec(lam, T):
            # lam^2 x_q + 2 lam (q+1) x_{q+1} + (q+1)(q+2) x_{q+2} = T_q
            if not T:
                return {}
            qmax, qmin = max(T), min(T)
            lo = min(qmin, 0) - 4
            x = {}
            for q_ in range(qmax, lo - 1, -1):
                v = aadd(T.get(q_, Z), ascale(x.get(q_ + 1, Z), lam * (2 * (q_ + 1))), -1)
                v = aadd(v, ascale(x.get(q_ + 2, Z), ctx.c((q_ + 1) * (q_ + 2))), -1)
                x[q_] = ascale(v, 1 / (lam * lam))
            conds.append(x[lo])
            conds.append(x[lo + 1])
            return x

        # 2. Ei keys, lam != 0: other-key image of r^p e^{lam r} Ei(mu r):  E(lam+mu): (2p-1) r^{p-2} + (2lam+mu) r^{p-1}
        for k in sorted([k for k in list(target) if k[0] == "Ei" and k[1] != 0], key=str):
            T = target.pop(k)
            x = solve_key_rec(k[1], T)
            sol[k] = x
            for p, v in x.items():
                sub_image(clean({("E", k[1] + k[2], 0): {p - 2: ctx.c(2 * p - 1), p - 1: 2 * k[1] + k[2]}}), v)
        # 3. Ei keys, lam = 0
        for i, mu in enumerate(mus):
            k = ("Ei", zero, mu)
            T = target.pop(k, {})
            x = {}
            for q_, v in T.items():
                p = q_ + 2
                if p in (0, 1):
                    conds.append(v)
                    continue
                x[p] = ascale(v, 1 / ctx.c(p * (p - 1)))
            x[0], x[1] = unit(3 + 2 * i), unit(4 + 2 * i)
            sol[k] = x
            for p, v in x.items():
                sub_image(clean({("E", mu, 0): {p - 2: ctx.c(2 * p - 1), p - 1: mu}}), v)
        if any(k[0] == "Ei" for k in target):
            raise RuntimeError("Ei key with lam = 0 and non-candidate rate: %s" % [k for k in target if k[0] == "Ei"])
        # 4. E keys
        for k in sorted(list(target), key=str):
            T = target.pop(k)
            if k[1] == 0:
                conds.extend(T.values())
                continue
            sol[k] = solve_key_rec(k[1], T)
        vals = self._solve_conds(conds, n)
        out: Func = {}
        for k, d in sol.items():
            for p, v in d.items():
                val = v[0] + sum(v[i] * vals[i] for i in range(1, n + 1))
                if val != 0:
                    out.setdefault(k, {})[p] = val
        return clean(out)


# ------------------------------------------------------------------ scalar jets and the source S_m
SCALAR_EXPRS = {
    "q": (u2s**2 - u3s**2) * (w2s**2 - w3s**2),
    "c": (u2s**2 - u3s**2 + w2s**2 - w3s**2) * (u2s**2 * w2s**2 - u3s**2 * w3s**2),
    "Sig": u2s**2 + u3s**2 + w2s**2 + w3s**2,
    "e": (u2s**2 - w3s**2) * (u3s**2 - w2s**2),
}


class RJets:
    """Scalar jets (q, c, Sig, e) and the elementary source S_m at the exchange point, as DERIVATIVE values."""

    def __init__(self, ctx: X.Ctx, D: Sequence[tuple], calc: MCalc):
        self.ctx, self.D, self.calc = ctx, list(D), calc
        self.J = X.Jets(ctx, D)
        J = self.J
        a, b = ctx.a, ctx.b
        self.base = [J.var(0, a), J.var(1, b), J.var(2, a), J.var(3, b)]       # u3, u2, w2, w3
        for nm, ex in SCALAR_EXPRS.items():
            setattr(self, nm, self.deriv(jet_eval(J, ex, self.base)))

    def deriv(self, S):
        return {al: self.ctx.c(_fact(al)) * v for al, v in S.items()}

    def source(self, m: int) -> Dict[tuple, Func]:
        """S_m = (-1)^{m+1} L^{-1}[d^m P / dw1^m at w1 = 0] (r > 0), theta-derivative values on D."""
        J, ctx, calc = self.J, self.ctx, self.calc
        tot: Dict[tuple, Func] = {}
        sgn = ctx.c((-1) ** (m + 1))
        for (C, cexpr, eps) in P_logs():
            lin = J.scale(jet_eval(J, cexpr, self.base), -1)
            Ej = J.exp_jet(lin)                                   # e^{-c(theta) r}
            for k in range(0, (m if eps else 0) + 1):
                Cj = sp.diff(C, w1s, m - k).subs(w1s, 0) * comb(m, k)
                if k == 0:
                    basefn = {al: fscale(Fn, ctx.c(-1), -1) for al, Fn in Ej.items()}          # -e^{-cr}/r
                else:
                    basefn = {al: fscale(Fn, ctx.c((-1) ** (k - 1)), k - 1) for al, Fn in Ej.items()}
                for i in range(0, 3):
                    ci = sp.diff(Cj, ts, i).subs(ts, 0) / factorial(i)
                    if ci == 0:
                        continue
                    cj = jet_eval(J, ci, self.base)
                    if not cj:
                        continue
                    fn = basefn
                    for _ in range(i):
                        fn = {al: calc.diff(Fn) for al, Fn in fn.items()}
                    tot = J.fj_add(tot, J.fj_scal(cj, fn))
                assert sp.diff(Cj, ts, 3) == 0
        return {al: fscale(tot.get(al, {}), sgn * ctx.c(_fact(al))) for al in self.D}


def binom(al, U) -> int:
    b = 1
    for s_, u_ in zip(al, U):
        b *= comb(s_, u_)
    return b


def subs_of(al):
    return itertools.product(*[range(x + 1) for x in al])


def coef_m2(m: int) -> int:
    return m * m


def coef_m3(m: int) -> int:
    return m * (m - 1) ** 2 * (m - 2)


def solve_fn(calc: MCalc, RJ: RJets, D: Sequence[tuple], m: int, flow: dict, flow3: dict = None):
    """f_{m+1} jets on D from (**); flow = jets of f_{m-1} on D, flow3 = jets of f_{m-3} (m >= 3).
    Returns (jets, [unknowns left free by each solve_D2])."""
    ctx = calc.ctx
    S = RJ.source(m)
    out: Dict[tuple, Func] = {}
    free: List[int] = []
    q0 = RJ.q[(0, 0, 0, 0)]
    assert RJ.c.get((0, 0, 0, 0), 0) == 0
    m2, m3 = ctx.c(coef_m2(m)), ctx.c(coef_m3(m))
    for al in D:
        rhs = dict(S[al])
        for U in subs_of(al):
            T = tuple(s - u for s, u in zip(al, U))
            bn = ctx.c(binom(al, U))
            if sum(U) and RJ.q.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(calc.D2(out[T]), bn * RJ.q[U]), -1)
            if sum(U) and RJ.c.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(out[T], bn * RJ.c[U]))
            lo = flow[T]
            if sum(U) == 0:
                rhs = fadd(rhs, fscale(calc.D2(calc.D2(lo)), m2))
            if RJ.Sig.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(calc.D2(lo), -m2 * bn * RJ.Sig[U]))
            if RJ.e.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(lo, m2 * bn * RJ.e[U]))
            if sum(U) == 0 and m >= 3:
                rhs = fadd(rhs, fscale(calc.D2(flow3[T]), m3))
        out[al] = calc.solve_D2(fscale(rhs, 1 / q0))
        free.append(calc.last_free)
    return out, free


# ------------------------------------------------------------------ the two cases
_s1A, _s1B, _s2A, _s2B, _v, _R = sp.symbols("s1A s1B s2A s2B v R", positive=True)


def density_m(case: str):
    """(polynomial in r1A^2, r1B^2, r2A^2, r2B^2, r12^2 and R; prefactor N2pi2(a, b) = pi^2 N_A^2 N_B^2 x sign)."""
    z1 = (_s1A + _R**2 - _s1B) / (2 * _R)
    z2 = (_s2A + _R**2 - _s2B) / (2 * _R)
    Xd = (_s1A + _s2A - _v) / 2 - z1 * z2
    return {
        "2ppi": (Xd, lambda a, b: -a**5 * b**3 / 2),
        "3ddel": (2 * Xd**2 - (_s1A - z1**2) * (_s2A - z2**2), lambda a, b: a**7 * b**3 / 12),
    }[case]


@lru_cache(maxsize=None)
def case_terms_m(case: str):
    """{k: [(alpha, {power of R: Fraction})]}: the r12^{2k} part, alpha = (2 i1 + 1, 2 j1 + 1, 2 i2 + 1, 2 j2 + 1)."""
    poly = sp.Poly(sp.expand(density_m(case)[0]), _s1A, _s1B, _s2A, _s2B, _v)
    terms: Dict[int, Dict[tuple, Dict[int, Fraction]]] = {}
    for (i1, j1, i2, j2, k), co in poly.terms():
        al = (2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1)
        d = terms.setdefault(int(k), {}).setdefault(al, {})
        for term in sp.Add.make_args(sp.expand(co)):
            cc, e = term.as_coeff_exponent(_R)
            cc = sp.Rational(cc)
            d[int(e)] = d.get(int(e), Fraction(0)) + Fraction(int(cc.p), int(cc.q))
    return {k: sorted(v.items()) for k, v in terms.items()}


@lru_cache(maxsize=None)
def downsets(case: str) -> Dict[int, Tuple[tuple, ...]]:
    """f_{2k} is needed on D_k = downset of the multi-indices of every r12^{2k'} part with k' >= k."""
    terms = case_terms_m(case)
    out, gens = {}, []
    for k in range(max(terms), -1, -1):
        gens = gens + [al for al, _ in terms.get(k, [])]
        out[k] = tuple(X.downset(gens))
    return out


def solve_case(ctx: X.Ctx, case: str):
    """(calc, {k: jets of f_{2k} on D_k}, free unknowns of the f_0 solves, [free unknowns of every solve_D2])."""
    Dk = downsets(case)
    calc0, f0, Fs, P2, nfree0 = X.solve_all(ctx, list(Dk[0]))
    calc = MCalc(ctx)
    fl = {0: f0}
    free: List[int] = []
    for k in range(1, max(Dk) + 1):
        RJ = RJets(ctx, Dk[k], calc)
        fl[k], fr = solve_fn(calc, RJ, Dk[k], 2 * k - 1, fl[k - 1], fl[k - 2] if k >= 2 else None)
        free += fr
    return calc, fl, nfree0, free


def assemble_K(ctx: X.Ctx, case: str, fl) -> Func:
    """K = 16 N2pi2 / R sum_k sum_alpha c_alpha(R) d^alpha f_{2k}."""
    pref = 16 * density_m(case)[1](ctx.a, ctx.b)
    K: Func = {}
    for k, lst in case_terms_m(case).items():
        for al, d in lst:
            for e, cc in d.items():
                K = fadd(K, fscale(fl[k][al], pref * ctx.c(cc), e - 1))
    return K
