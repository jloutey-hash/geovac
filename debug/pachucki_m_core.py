"""r12-weighted two-centre masters f_{m+1} (Pachucki's f_n, arXiv:0905.1635 label 05) at the exchange point,
as exact jets in the four exponents -- WITHOUT the 1/p^2 of labels 45/46.

Route.  Label 39 (eq_2, the first-order equation in the r12 exponent w1):
      sigma dg/dw1 + (1/2) dsigma/dw1 g + P(w1,u1; w2,u2; w3,u3) = 0,        sigma = label 12, P = label 13.
With u1 = t:  sigma = sigma0 + sigma2 w1^2 + sigma4 w1^4,
      sigma0 = c - q t^2,   c = (u2^2-u3^2+w2^2-w3^2)(u2^2 w2^2-u3^2 w3^2),   q = (u2^2-u3^2)(w2^2-w3^2),
      sigma2 = t^4 - Sig t^2 + e,   Sig = u2^2+u3^2+w2^2+w3^2,   e = (u2^2-w3^2)(u3^2-w2^2),      sigma4 = t^2.
g = sum_n (-w1)^n g_n / n!  (g_n = Laplace transform of f_n, label 40).  The w1^m coefficient of label 39 is
      sigma0 g_{m+1} + m^2 sigma2 g_{m-1} + m (m-1)^2 (m-2) sigma4 g_{m-3} + (-1)^{m+1} P_m = 0,  P_m = d^m P/dw1^m at 0.
Inverse Laplace transform for r > 0 (t -> d/dr; delta terms at r = 0 drop):
      q f_{m+1}'' = c f_{m+1} + m^2 (d^4 - Sig d^2 + e) f_{m-1} + m (m-1)^2 (m-2) f_{m-3}'' + S_m(r),          (**)
      S_m = (-1)^{m+1} L^{-1}[P_m].
(**) holds for ALL exponents; q = -(a^2-b^2)^2 != 0 at the exchange point; no division by p^2 anywhere.  (At the
exchange locus it is label 61 differentiated twice.)  P_m(t) is  sum poly(t) ln(t + const) + sum poly(t)/(t+const)^k
+ polynomial, so S_m is elementary: derivatives of e^{-c r}/r and r^{k-1} e^{-c r}, c in {u2+w3, u3+w2, u2+w2, u3+w3}.

theta-derivatives of (**) at the exchange point (Leibniz), solved lowest multi-index first; each is u'' = known, and
the solution is unique in the exponentially decaying span.
"""
from __future__ import annotations

import itertools
import sys
from fractions import Fraction
from math import comb, factorial
from typing import Dict, Sequence

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402

# sympy 1.14: PolyElement._gcd_ZZ is heugcd with NO fallback; it raises HeuristicGCDFailed on the large univariate
# polynomials met in the Q(t) runs.  Fall back to python-flint's exact fmpz_poly gcd (univariate rings only).
import flint as _flint  # noqa: E402
import sympy.polys.rings as _rings  # noqa: E402
from sympy.polys.heuristicgcd import heugcd as _heugcd  # noqa: E402
from sympy.polys.polyerrors import HeuristicGCDFailed as _HGF  # noqa: E402


def _gcd_ZZ_safe(f, g):
    try:
        return _heugcd(f, g)
    except _HGF:
        ring = f.ring
        assert ring.ngens == 1
        fa = _flint.fmpz_poly([int(c) for c in reversed(f.to_dense())])
        ga = _flint.fmpz_poly([int(c) for c in reversed(g.to_dense())])
        h = fa.gcd(ga)
        cf, cg = fa // h, ga // h
        assert cf * h == fa and cg * h == ga

        def back(q):
            return ring.from_dense([ring.domain(int(c)) for c in reversed(q.coeffs())])
        return back(h), back(cf), back(cg)


_rings.PolyElement._gcd_ZZ = _gcd_ZZ_safe

Func = E.Func
fadd, fscale, clean = E.fadd, E.fscale, E.clean

u2s, u3s, w2s, w3s, w1s, ts = sp.symbols("u2 u3 w2 w3 w1 t")
VARIDX = {u3s: 0, u2s: 1, w2s: 2, w3s: 3}


def P_logs():
    """t-dependent logarithms of P(w1,t; w2,u2; w3,u3) (label 13 with u1 = t): list of (C(w1,t,theta), c(theta), eps),
    the term being C * ln(t + c + eps w1).  Terms with t-independent logs are polynomial in t (no r > 0 content)."""
    u2, u3, w2, w3, w1, t = u2s, u3s, w2s, w3s, w1s, ts
    A3 = t**2 * w1**2 + u2**2 * w2**2 - u3**2 * w3**2 + w1 * w2 * (t**2 + u2**2 - w3**2)
    A4 = t**2 * w1**2 - u2**2 * w2**2 + u3**2 * w3**2 + w1 * w3 * (t**2 + u3**2 - w2**2)
    A5 = u2 * (u2 + w1) * (t**2 + u3**2 - w2**2) - u3**2 * (t**2 + u2**2 - w3**2)
    A6 = u3 * (u3 + w1) * (t**2 + u2**2 - w3**2) - u2**2 * (t**2 + u3**2 - w2**2)
    D3 = w3**2 - (w1 + w2)**2          # Gamma(t+u2, w3, w1+w2)
    D4 = w2**2 - (w1 + w3)**2          # Gamma(t+u3, w2, w1+w3)
    D5 = u3**2 - (u2 + w1)**2          # Gamma(t+w2, u3, u2+w1)
    D6 = u2**2 - (u3 + w1)**2          # Gamma(t+w3, u2, u3+w1)
    return [
        (-A3 / D3 + A6 / D6, u2 + w3, 0),
        (-A4 / D4 + A5 / D5, u3 + w2, 0),
        (t * w1 + A3 / D3 - A5 / D5, u2 + w2, 1),
        (t * w1 + A4 / D4 - A6 / D6, u3 + w3, 1),
    ]


def P_full():
    """P(w1,u1=t; w2,u2; w3,u3), label 13 verbatim (used only for checks)."""
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


def jet_eval(J: E.Jets, expr, base):
    """Evaluate a rational sympy expression in (u3,u2,w2,w3) as a scalar jet.  base = the four variable jets."""
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


def _fact(al):
    f = 1
    for x in al:
        f *= factorial(x)
    return f


class MCalc(E.Calc):
    """Adds the solver of u'' = src in the decaying span."""

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


class RJets:
    """Scalar jets (q, c, Sig, e) and the elementary source S_m, at the exchange point, as DERIVATIVE values."""

    def __init__(self, ctx: E.Ctx, D: Sequence[tuple], calc: MCalc):
        self.ctx, self.D, self.calc = ctx, list(D), calc
        self.J = E.Jets(ctx, D)
        J = self.J
        a, b = ctx.a, ctx.b
        self.base = [J.var(0, a), J.var(1, b), J.var(2, a), J.var(3, b)]       # u3, u2, w2, w3
        u2, u3, w2, w3 = u2s, u3s, w2s, w3s
        self.q = self.deriv(jet_eval(J, (u2**2 - u3**2) * (w2**2 - w3**2), self.base))
        self.c = self.deriv(jet_eval(J, (u2**2 - u3**2 + w2**2 - w3**2) * (u2**2 * w2**2 - u3**2 * w3**2), self.base))
        self.Sig = self.deriv(jet_eval(J, u2**2 + u3**2 + w2**2 + w3**2, self.base))
        self.e = self.deriv(jet_eval(J, (u2**2 - w3**2) * (u3**2 - w2**2), self.base))

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


def binom(al, U):
    b = 1
    for s_, u_ in zip(al, U):
        b *= comb(s_, u_)
    return b


def subs_of(al):
    return itertools.product(*[range(x + 1) for x in al])


def solve_fn(calc: MCalc, RJ: RJets, D: Sequence[tuple], m: int, flow: dict, flow3: dict = None, check=True,
             verbose=False, out: dict = None, stop=None, progress=None):
    """f_{m+1} jets on D from (**).  flow = jets of f_{m-1} (on D), flow3 = jets of f_{m-3} (m >= 3).
    out: already solved jets (resume); stop(): return True to interrupt; progress(out): called after each solve."""
    ctx = calc.ctx
    S = RJ.source(m)
    out = {} if out is None else out
    q0 = RJ.q[(0, 0, 0, 0)]
    for al in D:
        if al in out:
            continue
        if stop is not None and stop():
            break
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
                rhs = fadd(rhs, fscale(calc.D2(calc.D2(lo)), ctx.c(m * m)))
            if RJ.Sig.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(calc.D2(lo), -ctx.c(m * m) * bn * RJ.Sig[U]))
            if RJ.e.get(U, 0) != 0:
                rhs = fadd(rhs, fscale(lo, ctx.c(m * m) * bn * RJ.e[U]))
            if sum(U) == 0 and m >= 3:
                rhs = fadd(rhs, fscale(calc.D2(flow3[T]), ctx.c(m * (m - 1) ** 2 * (m - 2))))
        assert RJ.c.get((0, 0, 0, 0), 0) == 0
        out[al] = calc.solve_D2(fscale(rhs, 1 / q0))
        if check:
            assert fadd(fscale(calc.D2(out[al]), q0), rhs, -1) == {}, al
        if verbose:
            print("  f_%d" % (m + 1), al, len(out[al]), flush=True)
        if progress is not None:
            progress(out)
    return out


# ---------------------------------------------------------------- nine-function span
def nine_span(ctx):
    a, b, one, zero = ctx.a, ctx.b, ctx.one, ctx.zero

    def psi(x, y):
        return {("Ei", -2 * x, 2 * (x - y)): {0: one}, ("Ei", -2 * x, (x - y)): {0: -2 * one}, ("E", -2 * x, 1): {0: -one}}
    nine = [{("f",): {0: one}}, {("Ei", zero, -(a + b)): {0: one}}, {("E", -(a + b), 0): {0: one}},
            {("E", -2 * a, 0): {0: one}}, {("E", -2 * b, 0): {0: one}},
            {("Ei", 2 * a, -2 * (a + b)): {0: one}}, {("Ei", 2 * b, -2 * (a + b)): {0: one}}, psi(a, b), psi(b, a)]
    lead = [("f",), ("Ei", zero, -(a + b)), ("E", -(a + b), 0), ("E", -2 * a, 0), ("E", -2 * b, 0),
            ("Ei", 2 * a, -2 * (a + b)), ("Ei", 2 * b, -2 * (a + b)), ("Ei", -2 * a, 2 * (a - b)),
            ("Ei", -2 * b, 2 * (b - a))]
    return nine, lead


def in_nine(ctx, A):
    nine, lead = nine_span(ctx)
    rest = {k: dict(v) for k, v in A.items()}
    co = []
    for g, l in zip(nine, lead):
        c = dict(rest.get(l, {}))
        co.append(c)
        for p, x in c.items():
            rest = fadd(rest, fscale(g, x, p), -1)
    return co if rest == {} else None
