"""Pachucki exchange-point engine, generic coefficient field (extension of tests/_pachucki_exchange.py).

Differences from the test helper (which stays the reviewed reference; this file is checked against it):
  * sources d^alpha F and d^alpha p^2 are built as truncated Taylor JETS in the four exponents around the exchange
    point (transcription of arXiv:0905.1635 labels 20, 22-23, 30-34), for any downward-closed set of multi-indices
    alpha (not only subsets) -> needed for l > 0 densities (powers of r_1A, r_1B, r_2A, r_2B);
  * the equation (r f_alpha')' = -F_alpha + r sum_{0 != U <= alpha} binom(alpha, U) (d_U p^2) f_{alpha-U}
    is solved by a direct triangular solver (no generic rref), so the coefficient field can be Q (Fraction) or
    Q(t) (sympy FracField elements; a = 1, b = t by scaling).
Function algebra keys:  ('E', lam, j) r^k LC^j e^{lam r};  ('Ei', lam, mu) r^k e^{lam r} Ei(mu r);  ('f',) r^k f0.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Sequence, Tuple

Func = Dict[tuple, Dict[int, object]]


class Ctx:
    """Coefficient field. kind 'Q' (Fraction) or 'Qt' (sympy QQ(t)); a, b are field elements."""

    def __init__(self, kind: str, a=None, b=None):
        self.kind = kind
        if kind == "Q":
            self.one, self.zero = Fraction(1), Fraction(0)
            self.a, self.b = Fraction(a), Fraction(b)
        else:
            import sympy as sp
            self.t = sp.Symbol("t")
            self.K = sp.QQ.frac_field(self.t)
            self.one, self.zero = self.K.one, self.K.zero
            self.a, self.b = self.one, self.K.from_sympy(self.t)

    def c(self, x):
        if self.kind == "Q":
            return Fraction(x)
        x = Fraction(x)
        return self.K(x.numerator) / self.K(x.denominator)


# ------------------------------------------------------------------ Func algebra
def clean(A: Func) -> Func:
    out = {}
    for k, d in A.items():
        dd = {p: c for p, c in d.items() if c != 0}
        if dd:
            out[k] = dd
    return out


def fadd(A: Func, B: Func, s=1) -> Func:
    out = {k: dict(v) for k, v in A.items()}
    for k, d in B.items():
        o = out.setdefault(k, {})
        for p, c in d.items():
            o[p] = o.get(p, 0) + s * c
    return clean(out)


def fscale(A: Func, c, shift: int = 0) -> Func:
    return clean({k: {p + shift: c * v for p, v in d.items()} for k, d in A.items()})


def _kmul(k1, k2):
    if k1[0] == "E" and k2[0] == "E":
        return ("E", k1[1] + k2[1], k1[2] + k2[2])
    if k1[0] == "E" and k2[0] == "Ei":
        assert k1[2] == 0
        return ("Ei", k1[1] + k2[1], k2[2])
    if k1[0] == "Ei" and k2[0] == "E":
        return _kmul(k2, k1)
    raise ValueError((k1, k2))


def fmul(A: Func, B: Func) -> Func:
    out: Func = {}
    for k1, d1 in A.items():
        for k2, d2 in B.items():
            k = _kmul(k1, k2)
            o = out.setdefault(k, {})
            for p1, c1 in d1.items():
                for p2, c2 in d2.items():
                    o[p1 + p2] = o.get(p1 + p2, 0) + c1 * c2
    return clean(out)


# ------------------------------------------------------------------ jets in the 4 exponents (u3, u2, w2, w3)
def downset(gens: Sequence[tuple]) -> List[tuple]:
    s = set()
    for g in gens:
        for al in itertools.product(*[range(x + 1) for x in g]):
            s.add(al)
    return sorted(s, key=lambda x: (sum(x), x))


class Jets:
    def __init__(self, ctx: Ctx, D: Sequence[tuple]):
        self.ctx, self.D = ctx, list(D)
        self.Dset = set(self.D)
        self.splits = {al: [(be, tuple(x - y for x, y in zip(al, be)))
                            for be in itertools.product(*[range(x + 1) for x in al])] for al in self.D}
        z = (0, 0, 0, 0)
        self.Z = z

    # scalar jets: dict alpha -> field element
    def const(self, c):
        return {self.Z: c}

    def var(self, i, c0):
        e = [0, 0, 0, 0]
        e[i] = 1
        j = {self.Z: c0}
        if tuple(e) in self.Dset:
            j[tuple(e)] = self.ctx.one
        return j

    def add(self, A, B, s=1):
        out = dict(A)
        for k, v in B.items():
            out[k] = out.get(k, 0) + s * v
        return {k: v for k, v in out.items() if v != 0}

    def mul(self, A, B):
        out = {}
        for al in self.D:
            acc = 0
            for be, ga in self.splits[al]:
                x, y = A.get(be), B.get(ga)
                if x is not None and y is not None:
                    acc = acc + x * y
            if acc != 0:
                out[al] = acc
        return out

    def scale(self, A, c):
        return {k: c * v for k, v in A.items() if c * v != 0}

    def inv(self, A):
        a0 = A[self.Z]
        N = self.add(A, self.const(a0), -1)            # nilpotent part
        Nn = self.scale(N, -1 / a0)                    # -N/a0
        out, term = self.const(self.ctx.one), self.const(self.ctx.one)
        for _ in range(max(sum(al) for al in self.D) + 1):
            term = self.mul(term, Nn)
            if not term:
                break
            out = self.add(out, term)
        return self.scale(out, 1 / a0)

    def log_nilpart(self, A):
        """log(A/A0) (zero constant term)."""
        a0 = A[self.Z]
        N = self.scale(self.add(A, self.const(a0), -1), 1 / a0)
        out, term, k = {}, self.const(self.ctx.one), 0
        while True:
            k += 1
            term = self.mul(term, N)
            if not term:
                break
            out = self.add(out, self.scale(term, self.ctx.c(Fraction((-1) ** (k + 1), k))))
        return out

    def pow_list(self, N, kmax=60):
        out = [self.const(self.ctx.one)]
        while len(out) <= kmax:
            nxt = self.mul(out[-1], N)
            if not nxt:
                break
            out.append(nxt)
        return out

    # function jets: dict alpha -> Func
    def fj_scal(self, S, Fj):
        out = {}
        for al in self.D:
            acc = {}
            for be, ga in self.splits[al]:
                x, y = S.get(be), Fj.get(ga)
                if x is not None and y is not None:
                    acc = fadd(acc, fscale(y, x))
            if acc:
                out[al] = acc
        return out

    def fj_mul(self, A, B):
        out = {}
        for al in self.D:
            acc = {}
            for be, ga in self.splits[al]:
                x, y = A.get(be), B.get(ga)
                if x is not None and y is not None:
                    acc = fadd(acc, fmul(x, y))
            if acc:
                out[al] = acc
        return out

    def fj_add(self, A, B, s=1):
        out = {k: v for k, v in A.items()}
        for k, v in B.items():
            out[k] = fadd(out.get(k, {}), v, s)
        return {k: v for k, v in out.items() if v}

    def exp_jet(self, lin):
        """e^{lam(theta) r}; lin = scalar jet of lam (linear in theta)."""
        l0 = lin.get(self.Z, 0)
        N = {k: v for k, v in lin.items() if k != self.Z}
        out = {}
        for k, Pk in enumerate(self.pow_list(N)):
            fk = Fraction(1, factorial(k))
            for al, c in Pk.items():
                out.setdefault(al, {})
                out[al] = fadd(out[al], {("E", l0, 0): {k: c * self.ctx.c(fk)}})
        return out

    def ei_jet(self, lin):
        """Ei(mu(theta) r); d^k/dmu^k Ei(mu r) = sum_j C(k-1,j) r^j e^{mu r} (-1)^{k-1-j} (k-1-j)! mu^{-(k-j)}."""
        m0 = lin[self.Z]
        N = {k: v for k, v in lin.items() if k != self.Z}
        out = {self.Z: {("Ei", self.ctx.zero, m0): {0: self.ctx.one}}}
        for k, Pk in enumerate(self.pow_list(N)):
            if k == 0:
                continue
            der: Func = {}
            for j in range(k):
                c = comb(k - 1, j) * (-1) ** (k - 1 - j) * factorial(k - 1 - j)
                der = fadd(der, {("E", m0, 0): {j: self.ctx.c(c) / m0 ** (k - j)}})
            der = fscale(der, self.ctx.c(Fraction(1, factorial(k))))
            for al, c in Pk.items():
                out[al] = fadd(out.get(al, {}), fscale(der, c))
        return out

    def scal_to_fj(self, S):
        return {al: {("E", self.ctx.zero, 0): {0: c}} for al, c in S.items()}


def sources(ctx: Ctx, D: Sequence[tuple]):
    """Jets of F (labels 30-34) and p^2 (label 20) at the exchange point u3 = w2 = a, u2 = w3 = b.
    Returns (Fs, P2): Fs[alpha] = d^alpha F (Func), P2[alpha] = d^alpha p^2 (field element)."""
    J = Jets(ctx, D)
    a, b = ctx.a, ctx.b
    u3, u2, w2, w3 = J.var(0, a), J.var(1, b), J.var(2, a), J.var(3, b)
    add, mul, inv = J.add, J.mul, J.inv
    sq = lambda x: mul(x, x)  # noqa: E731
    neg = lambda x: J.scale(x, -1)  # noqa: E731
    # p^2, label 20
    p2 = mul(mul(add(add(sq(u2), sq(u3), -1), add(sq(w2), sq(w3), -1)), add(mul(sq(u2), sq(w2)), mul(sq(u3), sq(w3)), -1)),
             inv(mul(add(sq(u2), sq(u3), -1), add(sq(w2), sq(w3), -1))))
    # ln|x|, labels 33-34: x = (u2+w3-u3-w2)(u2+u3)(w2+w3)/((u2+u3+w2+w3)(u2-u3)(w2-w3)); ln|x| = LC + log(x/x0)
    S = add(add(u2, u3), add(w2, w3))
    x = mul(mul(mul(add(add(u2, w3), add(u3, w2), -1), add(u2, u3)), add(w2, w3)),
            inv(mul(mul(S, add(u2, u3, -1)), add(w2, w3, -1))))
    lnN = J.log_nilpart(x)
    lnabs = J.fj_add({J.Z: {("E", ctx.zero, 1): {0: ctx.one}}}, J.scal_to_fj(lnN))
    E = J.exp_jet
    Ei = J.ei_jet
    one_fj = {J.Z: {("E", ctx.zero, 0): {0: ctx.one}}}
    M = J.fj_mul
    F1m = J.fj_add(M(E(add(u3, w2, -1)), Ei(neg(add(u2, u3)))), M(E(add(w2, u3, -1)), Ei(neg(add(w2, w3)))), -1)
    F2m = J.fj_add(M(E(add(w3, u2, -1)), Ei(neg(add(w2, w3)))), M(E(add(u2, w3, -1)), Ei(neg(add(u2, u3)))), -1)
    br3 = J.fj_add(J.fj_add(J.fj_add(lnabs, Ei(add(add(u3, w2), add(u2, w3), -1)), -1), Ei(add(w2, w3, -1))),
                   Ei(add(u3, u2, -1)))
    F3m = J.fj_add(M(E(neg(add(u3, w2))), br3), M(E(add(u3, w2)), Ei(neg(S))), -1)
    br4 = J.fj_add(J.fj_add(J.fj_add(lnabs, Ei(add(add(u2, w3), add(u3, w2), -1)), -1), Ei(add(w3, w2, -1))),
                   Ei(add(u2, u3, -1)))
    F4m = J.fj_add(M(E(neg(add(u2, w3))), br4), M(E(add(u2, w3)), Ei(neg(S))), -1)
    du = inv(add(sq(u2), sq(u3), -1))
    dw = inv(add(sq(w2), sq(w3), -1))
    ap = add(mul(w2, dw), mul(u3, du))
    am = add(mul(w2, dw), mul(u3, du), -1)
    bp = add(mul(u2, du), mul(w3, dw))
    bm = add(mul(u2, du), mul(w3, dw), -1)
    half = ctx.c(Fraction(1, 2))
    F = {}
    for co, Fi in ((ap, F1m), (bp, F2m), (am, F3m), (bm, F4m)):
        F = J.fj_add(F, J.fj_scal(J.scale(co, half), Fi))
    # strip zero-rate E keys produced by e^{0 r}: they are genuine (constants/powers times LC^j); keep as is
    Fs, P2 = {}, {}
    for al in D:
        fac = 1
        for x_ in al:
            fac *= factorial(x_)
        Fs[al] = fscale(F.get(al, {}), ctx.c(fac))
        P2[al] = ctx.c(fac) * p2.get(al, 0)
    return Fs, P2


# ------------------------------------------------------------------ calculus
class Calc:
    def __init__(self, ctx: Ctx):
        self.ctx = ctx
        a, b = ctx.a, ctx.b
        c = 1 / (2 * (a * a - b * b))
        self.FM1 = clean({   # label 60 at u = a, w = b; ln|(w-u)/(w+u)| = -LC
            ("Ei", 2 * b, -2 * (a + b)): {0: c},
            ("Ei", 2 * a, -2 * (a + b)): {0: -c},
            ("Ei", -2 * b, 2 * (b - a)): {0: -c},
            ("Ei", -2 * b, (b - a)): {0: 2 * c},
            ("E", -2 * b, 1): {0: c},
            ("Ei", -2 * a, 2 * (a - b)): {0: c},
            ("Ei", -2 * a, (a - b)): {0: -2 * c},
            ("E", -2 * a, 1): {0: -c},
        })
        self.F0 = self.diff(self.FM1)

    def diff(self, A: Func) -> Func:
        out: Func = {}

        def add(k, p, c):
            o = out.setdefault(k, {})
            o[p] = o.get(p, 0) + c

        for k, d in A.items():
            for p, c in d.items():
                if p != 0:
                    add(k, p - 1, p * c)
                if k[0] == "E":
                    if k[1] != 0:
                        add(k, p, k[1] * c)
                elif k[0] == "Ei":
                    lam, mu = k[1], k[2]
                    if lam != 0:
                        add(k, p, lam * c)
                    add(("E", lam + mu, 0), p - 1, c)
                elif k[0] == "f":
                    for kk, dd in self.FM1.items():
                        for pp, cc in dd.items():
                            add(kk, p + pp - 1, -c * cc)
        return clean(out)

    def L0(self, A: Func) -> Func:
        return self.diff(fscale(self.diff(A), 1, 1))

    # ---- direct solver of (r u')' = src in the decaying span
    def solve_L0(self, src: Func) -> Func:
        """Unique solution u of (r u')' = src in the decaying span
        {r^k f0, r^k e^{lam r} Ei(mu r) (lam+mu < 0), r^k e^{lam r} LC^j (lam < 0)}, or raise.
        The triangular recursion leaves free x_{f,0} and x_{Ei(0,mu),0} for every candidate E1 rate mu < 0;
        every coordinate is affine in these unknowns; termination conditions at the bottom fix them."""
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
        mus = sorted([m for m in rates if m < 0] if self.ctx.kind == "Q" else
                     [m for m in rates if self._neg(m)], key=str)
        n = 1 + len(mus)
        Z = [zero] * (n + 1)

        def aff(v):
            return [v] + [zero] * n

        def unit(i):
            u = [zero] * (n + 1)
            u[i] = ctx.one
            return u

        def aadd(x, y, s=1):
            return [xi + s * yi for xi, yi in zip(x, y)]

        def ascale(x, c):
            return [xi * c for xi in x]

        target = {k: {p: aff(c) for p, c in d.items()} for k, d in src.items()}
        sol = {}
        conds = []

        def sub_image(img: Func, coeff):
            for k, d in img.items():
                o = target.setdefault(k, {})
                for p, c in d.items():
                    o[p] = aadd(o.get(p, Z), ascale(coeff, c), -1)

        # 1. f coordinate: p^2 x_p = T_{p-1}; x_0 = unknown 1
        tf = target.pop(("f",), {})
        xf = {}
        for p1, v in tf.items():
            p = p1 + 1
            if p == 0:
                conds.append(v)
                continue
            xf[p] = ascale(v, 1 / ctx.c(p * p))
        xf[0] = unit(1)
        for p, x in xf.items():
            img = fadd(fscale(self.FM1, ctx.c(-2 * p), p - 1), fscale(self.F0, ctx.c(-1), p))
            sub_image(img, x)
        sol[("f",)] = xf

        def solve_key_rec(lam, T):
            if not T:
                return {}
            qmax, qmin = max(T), min(T)
            lo = min(qmin, 0) - 3
            x = {}
            for q in range(qmax, lo, -1):
                v = aadd(T.get(q, Z), ascale(x.get(q, Z), lam * (2 * q + 1)), -1)
                v = aadd(v, ascale(x.get(q + 1, Z), ctx.c((q + 1) ** 2)), -1)
                x[q - 1] = ascale(v, 1 / (lam * lam))
            conds.append(x[lo])
            conds.append(x[lo + 1])
            return x

        # 2. Ei keys with lam != 0
        for k in sorted([k for k in list(target) if k[0] == "Ei" and k[1] != 0], key=str):
            T = target.pop(k)
            x = solve_key_rec(k[1], T)
            sol[k] = x
            for p, v in x.items():
                sub_image(clean({("E", k[1] + k[2], 0): {p: 2 * k[1] + k[2], p - 1: ctx.c(2 * p)}}), v)
        # 3. Ei keys with lam == 0: p^2 x_p = T_{p-1}, x_0 free (one unknown per candidate rate)
        for i, mu in enumerate(mus):
            k = ("Ei", zero, mu)
            T = target.pop(k, {})
            x = {}
            for q, v in T.items():
                p = q + 1
                if p == 0:
                    conds.append(v)
                    continue
                x[p] = ascale(v, 1 / ctx.c(p * p))
            x[0] = unit(2 + i)
            sol[k] = x
            for p, v in x.items():
                sub_image(clean({("E", mu, 0): {p: mu, p - 1: ctx.c(2 * p)}}), v)
        if any(k[0] == "Ei" for k in target):
            raise RuntimeError("Ei key with lam = 0 and non-candidate rate: %s" % [k for k in target if k[0] == "Ei"])
        # 4. E keys
        for k in sorted(list(target), key=str):
            T = target.pop(k)
            if k[1] == 0:
                conds.extend(T.values())          # pure powers do not decay: must vanish
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

    def _neg(self, m):
        """m < 0 for an element of Q(t) linear in t, evaluated for t > 1 (both t = 2 and t -> oo)."""
        import sympy as sp
        e = self.ctx.K.to_sympy(m)
        return bool(e.subs(self.ctx.t, 2) < 0 and sp.limit(e, self.ctx.t, sp.oo) <= 0
                    and e.subs(self.ctx.t, 1) <= 0)

    def _solve_conds(self, conds, n):
        """Solve rows [c0, c1..cn] : c0 + sum ci xi = 0 for x; raise if inconsistent; record #free."""
        ctx = self.ctx
        M = [list(r) for r in conds if any(x != 0 for x in r)]
        piv = []
        row = 0
        for j in range(1, n + 1):
            pr = next((i for i in range(row, len(M)) if M[i][j] != 0), None)
            if pr is None:
                continue
            M[row], M[pr] = M[pr], M[row]
            pv = M[row][j]
            M[row] = [x / pv for x in M[row]]
            for i in range(len(M)):
                if i != row and M[i][j] != 0:
                    f = M[i][j]
                    M[i] = [x - f * y for x, y in zip(M[i], M[row])]
            piv.append(j)
            row += 1
        for i in range(row, len(M)):
            if M[i][0] != 0:
                raise RuntimeError("inconsistent: not in span")
        vals = [ctx.zero] * (n + 1)
        for i, j in enumerate(piv):
            vals[j] = -M[i][0]
        self.last_free = n - len(piv)
        return vals

def residual(calc: Calc, al, fS, Fs, P2):
    """(r f_al')' + F_al - r sum binom(al,U) P2[U] f_{al-U}  (should be identically 0)."""
    ctx = calc.ctx
    res = fadd(calc.L0(fS[al]), Fs[al])
    for U in itertools.product(*[range(x + 1) for x in al]):
        if sum(U) == 0:
            continue
        T = tuple(s - u for s, u in zip(al, U))
        bn = 1
        for s_, u_ in zip(al, U):
            bn *= comb(s_, u_)
        if P2[U] != 0:
            res = fadd(res, fscale(fS[T], -ctx.c(bn) * P2[U], 1))
    return res


def solve_all(ctx: Ctx, D: Sequence[tuple], verbose=False, check=True):
    calc = Calc(ctx)
    Fs, P2 = sources(ctx, D)
    assert P2[(0, 0, 0, 0)] == 0
    fS = {(0, 0, 0, 0): {("f",): {0: ctx.one}}}
    nfree = 0
    for al in D:
        if sum(al) == 0:
            continue
        src = fscale(Fs[al], -1)
        for U in itertools.product(*[range(x + 1) for x in al]):
            if sum(U) == 0 or P2[U] == 0:
                continue
            T = tuple(s - u for s, u in zip(al, U))
            bn = 1
            for s_, u_ in zip(al, U):
                bn *= comb(s_, u_)
            src = fadd(src, fscale(fS[T], ctx.c(bn) * P2[U], 1))
        fS[al] = calc.solve_L0(src)
        nfree += calc.last_free
        if check:
            assert residual(calc, al, fS, Fs, P2) == {}, al
        if verbose:
            print(al, len(fS[al]), flush=True)
    return calc, fS, Fs, P2, nfree


def keyset(A: Func):
    return sorted({k for k in A}, key=str)
