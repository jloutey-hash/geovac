"""Helper for tests/test_paper58_pachucki_extension.py -- Pachucki's ODE at the exchange point, any multi-index,
coefficient field Q or Q(t).

Extends the reviewed helper tests/_pachucki_exchange.py (imported as PX; its transcription of arXiv:0905.1635
labels 20, 22-23, 30-34, 59-60 stays the reference and the jet engine here is tested against it).  Two changes:

  * the sources d^alpha F and d^alpha p^2 are truncated Taylor JETS in the four exponents (u3, u2, w2, w3) around
    the exchange point u3 = w2 = a, u2 = w3 = b, for any downward-closed set of multi-indices alpha (PX covers
    subsets only: every order <= 1).  Needed for l > 0 densities, which are polynomials in r_A^2, r_B^2:
        r^n e^{-u r} = (-1)^{n+1} d_u^{n+1} [e^{-u r} / r];
  * the equation
        (r f_alpha')' = -F_alpha + r sum_{0 != U <= alpha} binom(alpha, U) (d_U p^2) f_{alpha-U}            (*)
    is solved by a direct triangular solver (no ansatz window), so the field can be Q (Fraction) or Q(t)
    (sympy FracField; a = 1, b = t, by the scaling K(R; a, b) = a K(aR; 1, b/a)).

Function algebra keys (as in PX):  ('E', lam, j) r^k LC^j e^{lam r};  ('Ei', lam, mu) r^k e^{lam r} Ei(mu r);
('f',) r^k f0.  LC = ln((a+b)/(b-a)) is an indeterminate.

Also here, and independent of all of the above: `neumann_integral`, a direct evaluation of the two-centre
integrals in prolate spheroidal coordinates (Neumann expansion of 1/r12); it never touches the closed form.

Nothing here reads or imports debug/.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from typing import Dict, List, Sequence

import sympy as sp

import _pachucki_exchange as PX

Func = Dict[tuple, Dict[int, object]]
CASES = ("1s", "2psig", "3dsig", "ppsig", "pppi0")
NAMES9 = ["f0", "-E1((a+b)R)", "e^{-(a+b)R}", "e^{-2aR}", "e^{-2bR}", "-e^{2aR}E1(2(a+b)R)", "-e^{2bR}E1(2(a+b)R)",
          "psi_a", "psi_b"]


class Ctx:
    """Coefficient field.  kind 'Q' (Fraction; a, b given) or 'Qt' (sympy QQ(t); a = 1, b = t)."""

    def __init__(self, kind: str, a=None, b=None):
        self.kind = kind
        if kind == "Q":
            self.one, self.zero = Fraction(1), Fraction(0)
            self.a, self.b = Fraction(a), Fraction(b)
        elif kind == "Qt":
            self.t = sp.Symbol("t")
            self.K = sp.QQ.frac_field(self.t)
            self.one, self.zero = self.K.one, self.K.zero
            self.a, self.b = self.one, self.K.from_sympy(self.t)
        else:
            raise ValueError(kind)

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
    """Truncated Taylor coefficients (NOT derivatives: coefficient of theta^alpha) on a downset D."""

    def __init__(self, ctx: Ctx, D: Sequence[tuple]):
        self.ctx, self.D = ctx, list(D)
        self.Dset = set(self.D)
        self.splits = {al: [(be, tuple(x - y for x, y in zip(al, be)))
                            for be in itertools.product(*[range(x + 1) for x in al])] for al in self.D}
        self.Z = (0, 0, 0, 0)

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
    """d^alpha F (labels 30-34, as a Func) and d^alpha p^2 (label 20, a field element) at the exchange point,
    for every alpha in the downset D."""
    J = Jets(ctx, D)
    a, b = ctx.a, ctx.b
    u3, u2, w2, w3 = J.var(0, a), J.var(1, b), J.var(2, a), J.var(3, b)
    add, mul, inv = J.add, J.mul, J.inv
    sq = lambda x: mul(x, x)  # noqa: E731
    neg = lambda x: J.scale(x, -1)  # noqa: E731
    # p^2, label 20
    p2 = mul(mul(add(add(sq(u2), sq(u3), -1), add(sq(w2), sq(w3), -1)),
                 add(mul(sq(u2), sq(w2)), mul(sq(u3), sq(w3)), -1)),
             inv(mul(add(sq(u2), sq(u3), -1), add(sq(w2), sq(w3), -1))))
    # ln|x|, labels 33-34: x = (u2+w3-u3-w2)(u2+u3)(w2+w3)/((u2+u3+w2+w3)(u2-u3)(w2-w3)); ln|x| = LC + log(x/x0)
    S = add(add(u2, u3), add(w2, w3))
    x = mul(mul(mul(add(add(u2, w3), add(u3, w2), -1), add(u2, u3)), add(w2, w3)),
            inv(mul(mul(S, add(u2, u3, -1)), add(w2, w3, -1))))
    lnN = J.log_nilpart(x)
    lnabs = J.fj_add({J.Z: {("E", ctx.zero, 1): {0: ctx.one}}}, J.scal_to_fj(lnN))
    E = J.exp_jet
    Ei = J.ei_jet
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
    Fs, P2 = {}, {}
    for al in D:
        fac = 1
        for x_ in al:
            fac *= factorial(x_)                       # Taylor coefficient -> derivative
        Fs[al] = fscale(F.get(al, {}), ctx.c(fac))
        P2[al] = ctx.c(fac) * p2.get(al, 0)
    return Fs, P2


# ------------------------------------------------------------------ calculus and the solver
class Calc:
    def __init__(self, ctx: Ctx):
        self.ctx = ctx
        a, b = ctx.a, ctx.b
        c = 1 / (2 * (a * a - b * b))
        self.FM1 = clean({   # label 60 at u = a, w = b; ln|(w-u)/(w+u)| = -LC   (same entries as PX.Algebra.fm1)
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
        self.last_free = 0

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
                    add(("E", lam + mu, 0), p - 1, c)          # d Ei(mu r) = e^{mu r} / r
                elif k[0] == "f":
                    for kk, dd in self.FM1.items():            # f0' = -F^(-1)/r, label 59
                        for pp, cc in dd.items():
                            add(kk, p + pp - 1, -c * cc)
                else:
                    raise KeyError(k)
        return clean(out)

    def L0(self, A: Func) -> Func:
        """(r A')'"""
        return self.diff(fscale(self.diff(A), 1, 1))

    def solve_L0(self, src: Func) -> Func:
        """The solution u of (r u')' = src in the decaying span
        {r^k f0, r^k e^{lam r} Ei(mu r) (lam + mu < 0), r^k e^{lam r} LC^j (lam < 0)}, or raise.
        The triangular recursion leaves free x_{f,0} and x_{Ei(0,mu),0} for every candidate E1 rate mu < 0; every
        coordinate is affine in these unknowns; the termination conditions fix them (self.last_free counts the
        unknowns left undetermined)."""
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
        """m < 0 for an element of Q(t) linear in t, for t > 1 (heuristic for the candidate list only: the
        solution is checked afterwards by the residual and by `audit_intermediates`)."""
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

    def derivatives(self, A: Func, n: int) -> List[Func]:
        out = [A]
        for _ in range(n):
            out.append(self.diff(out[-1]))
        return out


def multibinom(al, U) -> int:
    bn = 1
    for s_, u_ in zip(al, U):
        bn *= comb(s_, u_)
    return bn


def solve_all(ctx: Ctx, D: Sequence[tuple]):
    """Solve (*) for every alpha in the downset D, lowest first.  Returns (calc, fS, Fs, P2, number of free
    unknowns summed over alpha).  The residual is NOT checked here (the test does it with its own formula)."""
    calc = Calc(ctx)
    Fs, P2 = sources(ctx, D)
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
            src = fadd(src, fscale(fS[T], ctx.c(multibinom(al, U)) * P2[U], 1))
        fS[al] = calc.solve_L0(src)
        nfree += calc.last_free
    return calc, fS, Fs, P2, nfree


# ------------------------------------------------------------------ the l >= 0 map
_rA2, _rB2, _R = sp.symbols("rA2 rB2 R", positive=True)


def density(case: str):
    """(P, N2pi2) with P(r_A^2, r_B^2; R) the per-electron polynomial of the orbital product
    N_A N_B P e^{-a r_A - b r_B} and N2pi2(a, b) = pi^2 N_A^2 N_B^2 as a function of sympy/field a, b.
    A at the origin, B at z = R:  z_A = (r_A^2 + R^2 - r_B^2)/(2R),  z_B = z_A - R."""
    zA = (_rA2 + _R ** 2 - _rB2) / (2 * _R)
    zB = zA - _R
    return {
        "1s": (sp.Integer(1), lambda a, b: a ** 3 * b ** 3),
        "2psig": (zA, lambda a, b: a ** 5 * b ** 3),                           # 2p0 on A (rate a), 1s on B
        "3dsig": (3 * zA ** 2 - _rA2, lambda a, b: a ** 7 * b ** 3 / 18),      # 3d0 on A (rate a), 1s on B
        "ppsig": (zA * zB, lambda a, b: a ** 5 * b ** 5),                      # 2p0 on A, 2p0 on B
        "pppi0": (_rA2 - zA ** 2, lambda a, b: a ** 5 * b ** 5 / 4),           # conj(2p+1_A) 2p+1_B: (x^2+y^2)/2
    }[case]


@lru_cache(maxsize=None)
def case_terms(case: str):
    """[(alpha, {power of R: Fraction})]: K = 16 N2pi2 / R * sum_alpha c_alpha(R) d^alpha f, alpha in the order
    (u3, u2, w2, w3) = exponents of r1A, r1B, r2A, r2B.  A monomial r^n (n even here) costs n + 1 derivatives and
    the sign (-1)^{n+1}; the four signs multiply to +1."""
    P = sp.Poly(sp.expand(density(case)[0]), _rA2, _rB2)
    mono = list(zip(P.monoms(), P.coeffs()))
    out = {}
    for (i1, j1), c1 in mono:
        for (i2, j2), c2 in mono:
            al = (2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1)
            d = out.setdefault(al, {})
            for term in sp.Add.make_args(sp.expand(c1 * c2)):
                cc, k = term.as_coeff_exponent(_R)
                cc = sp.Rational(cc)
                d[int(k)] = d.get(int(k), Fraction(0)) + Fraction(int(cc.p), int(cc.q))
    return sorted(out.items())


def case_downset(case: str) -> List[tuple]:
    return downset([al for al, _ in case_terms(case)])


def assemble_K(ctx: Ctx, case: str, fS) -> Func:
    """K = N_A^2 N_B^2 (4 pi)^2 / R * sum c_alpha(R) f_alpha."""
    n2 = density(case)[1](ctx.a, ctx.b)
    K: Func = {}
    for al, d in case_terms(case):
        for k, cc in d.items():
            K = fadd(K, fscale(fS[al], 16 * n2 * ctx.c(cc), k - 1))
    return K


def nine_functions(ctx: Ctx):
    """[(function as a Func, its leading key)] for the nine functions, in the order NAMES9."""
    a, b, one = ctx.a, ctx.b, ctx.one

    def psi(x, y):   # e^{-2xR}[Ei(2(x-y)R) - 2 Ei((x-y)R) - LC]
        return {("Ei", -2 * x, 2 * (x - y)): {0: one}, ("Ei", -2 * x, (x - y)): {0: -2 * one},
                ("E", -2 * x, 1): {0: -one}}

    return [
        ({("f",): {0: one}}, ("f",)),
        ({("Ei", ctx.zero, -(a + b)): {0: one}}, ("Ei", ctx.zero, -(a + b))),
        ({("E", -(a + b), 0): {0: one}}, ("E", -(a + b), 0)),
        ({("E", -2 * a, 0): {0: one}}, ("E", -2 * a, 0)),
        ({("E", -2 * b, 0): {0: one}}, ("E", -2 * b, 0)),
        ({("Ei", 2 * a, -2 * (a + b)): {0: one}}, ("Ei", 2 * a, -2 * (a + b))),
        ({("Ei", 2 * b, -2 * (a + b)): {0: one}}, ("Ei", 2 * b, -2 * (a + b))),
        (psi(a, b), ("Ei", -2 * a, 2 * (a - b))),
        (psi(b, a), ("Ei", -2 * b, 2 * (b - a))),
    ]


def in_nine(ctx: Ctx, A: Func):
    """Laurent coefficient dicts of A on the nine functions, or None if A is not in their span."""
    rest = {k: dict(v) for k, v in A.items()}
    co = []
    for g, lead in nine_functions(ctx):
        c = dict(rest.get(lead, {}))
        co.append(c)
        for p, x in c.items():
            rest = fadd(rest, fscale(g, x, p), -1)
    return co if rest == {} else None


# ------------------------------------------------------------------ what "for every t > 1" needs
def audit_intermediates(ctx: Ctx, D, fS, Fs, P2) -> dict:
    """Over Q(t): irreducible factors (with the largest multiplicity) of every coefficient denominator of
    f_alpha, F_alpha, d_alpha p^2, alpha in D; and every function key with its decay rate."""
    KK, tS = ctx.K, ctx.t

    def den_factors(x, acc):
        d = sp.Poly(KK.to_sympy(x.denom), tS)
        if d.degree() <= 0:
            return
        for f, m in sp.factor_list(d.as_expr())[1]:
            acc[str(f)] = max(acc.get(str(f), 0), int(m))

    dens = {"f_alpha": {}, "F_alpha": {}, "d_alpha_p2": {}}
    keys = set()
    for al in D:
        for k, d in fS[al].items():
            keys.add(k)
            for c in d.values():
                den_factors(c, dens["f_alpha"])
        for k, d in Fs[al].items():
            keys.add(k)
            for c in d.values():
                den_factors(c, dens["F_alpha"])
        if P2[al] != 0:
            den_factors(P2[al], dens["d_alpha_p2"])
    rows = []
    for k in sorted(keys, key=str):
        if k[0] == "f":
            rows.append({"kind": "f"})
        elif k[0] == "E":
            rows.append({"kind": "E", "lam": str(sp.expand(KK.to_sympy(k[1]))), "log_power": int(k[2])})
        else:
            rows.append({"kind": "Ei", "lam": str(sp.expand(KK.to_sympy(k[1]))),
                         "mu": str(sp.expand(KK.to_sympy(k[2])))})
    return {"denominator_factors_max_multiplicity": dens, "keys": rows}


def key_verdict(row: dict) -> dict:
    """For a key row of `audit_intermediates` (lam, mu as strings, linear in t): the decay rate, whether it is
    < 0 for EVERY real t > 1, and the t >= 1 at which an Ei argument vanishes."""
    t = sp.Symbol("t")
    if row["kind"] == "f":
        return {"rate": None, "decays": True, "ei_zero": []}
    lam = sp.sympify(row["lam"], locals={"t": t})
    mu = sp.sympify(row["mu"], locals={"t": t}) if row["kind"] == "Ei" else sp.Integer(0)
    rate = sp.expand(lam + mu)
    for e in (lam, mu, rate):
        assert sp.Poly(e, t).degree() <= 1, e
    r1, slope = rate.subs(t, 1), sp.diff(rate, t)
    decays = bool(r1 <= 0 and slope <= 0 and not (r1 == 0 and slope == 0))   # linear: < 0 on (1, oo)
    ei_zero = []
    if row["kind"] == "Ei":
        if mu == 0:
            ei_zero = ["identically"]
        else:
            ei_zero = [str(x) for x in sp.solve(mu, t) if x.is_real and x > 1]
    return {"rate": str(rate), "decays": decays, "ei_zero": ei_zero}


# ------------------------------------------------------------------ conversions
def qt_to_fraction(ctx: Ctx, x, tv: Fraction) -> Fraction:
    if isinstance(x, (int, Fraction)):
        return Fraction(x)
    def ev(poly):
        return sum((Fraction(int(c)) * tv ** e for (e,), c in poly.terms()), Fraction(0))

    return ev(x.numer) / ev(x.denom)


def specialise(ctx: Ctx, A: Func, tv: Fraction) -> Func:
    """A Func over Q(t) at t = tv, as a Func over Q (keys and coefficients Fractions)."""
    out: Func = {}
    for k, d in A.items():
        kk = (k[0],) + tuple(x if (i == 1 and k[0] == "E") else qt_to_fraction(ctx, x, tv)
                             for i, x in enumerate(k[1:]))
        o = out.setdefault(kk, {})
        for p, c in d.items():
            o[p] = o.get(p, Fraction(0)) + qt_to_fraction(ctx, c, tv)
    return clean(out)


def coeff_strings_at(cft: Dict[str, Dict[str, str]], tv: Fraction) -> List[Dict[int, Fraction]]:
    """Stored closed-form coefficients {name: {power: rational function of t as a string}} at t = tv."""
    t = sp.Symbol("t")
    tq = sp.Rational(tv.numerator, tv.denominator)
    out = []
    for nm in NAMES9:
        d = {}
        for p, s in cft.get(nm, {}).items():
            v = sp.Rational(sp.sympify(s, locals={"t": t}).subs(t, tq))
            if v != 0:
                d[int(p)] = Fraction(int(v.p), int(v.q))
        out.append(d)
    return out


# ------------------------------------------------------------------ numerics of the closed form (mpmath)
def evaluate(a: Fraction, b: Fraction, A: Func, x, f0=None):
    """A(x) for a Func over Q at the ambient mpmath precision; LC = ln((a+b)/(b-a))."""
    import mpmath as mp

    def q(v):
        v = Fraction(v)
        return mp.mpf(v.numerator) / v.denominator

    lc = mp.log(q((a + b) / (b - a)))
    s = mp.mpf(0)
    for k, d in A.items():
        if k[0] == "f":
            base = f0
        elif k[0] == "E":
            base = mp.exp(q(k[1]) * x) * lc ** k[2]
        else:
            base = mp.exp(q(k[1]) * x) * mp.ei(q(k[2]) * x)
        for p, c in d.items():
            s += q(c) * x ** p * base
    return s


def f0_quadrature(a: Fraction, b: Fraction, x, fm1: Func = None):
    """f0(x) = int_x^oo F^(-1)(s)/s ds (label 59; F^(-1) = label 60), the closed form's only quadrature."""
    import mpmath as mp
    if fm1 is None:
        fm1 = Calc(Ctx("Q", a, b)).FM1
    return mp.quad(lambda s: evaluate(a, b, fm1, s) / s, [x, x + 1, x + 2, x + 4, x + 10, x + 30, mp.inf])


# ------------------------------------------------------------------ independent integrals (no closed form)
def neumann_integral(case: str, a, b, R, dps: int = 30, lmax: int = 40, guard: int = 25):
    """The two-centre integral by direct integration.  Prolate spheroidal coordinates xi = (rA + rB)/R,
    eta = (rA - rB)/R; Neumann expansion of 1/r12, m = 0 part (the densities are axially symmetric):

        1/r12 -> (2/R) sum_l (2l+1) P_l(xi<) Q_l(xi>) P_l(eta1) P_l(eta2),
        d^3r = (R/2)^3 (xi^2 - eta^2) dxi deta dphi.

    For a density W(xi, eta) e^{-al xi - be eta} (W polynomial, Jacobian included; al = R(a+b)/2, be = R(a-b)/2):

        I = pref * sum_l (2l+1) 2 int_1^oo Q_l g_l e^{-al x} [int_1^x P_l g_l e^{-al y} dy] dx,
        g_l(xi) = int_{-1}^{1} W(xi, eta) P_l(eta) e^{-be eta} deta.

    The inner integral is elementary; the eta moments and the outer integral are mpmath quadratures.
    case: one of CASES (exchange integral of normalised STO products), or 'f0' (Pachucki's master integral, label
    01, at the exchange point).  Shares no code with the engine above.  Working precision dps + guard digits; the
    l sum stops when a term is below 10^-(dps + min(8, guard/3)) of the total.  Returns (value, per-l terms)."""
    import mpmath as mp
    a, b, R = (sp.Rational(str(Fraction(v))) for v in (a, b, R))
    xi, eta = sp.symbols("xi eta")
    rA = R * (xi + eta) / 2
    zA = R * (1 + xi * eta) / 2          # = (rA^2 + R^2 - rB^2)/(2R) with rB = R(xi - eta)/2
    zB = zA - R
    jac = xi ** 2 - eta ** 2
    if case == "f0":
        W, N2pi2 = sp.Integer(1), None    # 1/(rA rB) = (2/R)^2 / jac cancels the Jacobian
    else:
        W, N2pi2 = {
            "1s": (jac, a ** 3 * b ** 3),
            "2psig": (jac * zA, a ** 5 * b ** 3),
            "3dsig": (jac * (3 * zA ** 2 - rA ** 2), a ** 7 * b ** 3 / 18),
            "ppsig": (jac * zA * zB, a ** 5 * b ** 5),
            "pppi0": (jac * (rA ** 2 - zA ** 2), a ** 5 * b ** 5 / 4),
        }[case]
    Wp = sp.Poly(sp.expand(W), xi, eta)
    with mp.workdps(dps + guard):
        al = mp.mpf(R * (a + b) / 2)
        be = mp.mpf(R * (a - b) / 2)
        Rm = mp.mpf(R)
        kmax = max(m[1] for m in Wp.monoms())
        total = mp.mpf(0)
        per = []
        for l in range(lmax + 1):
            Pl = [mp.mpf(c) for c in reversed(sp.Poly(sp.legendre(l, xi), xi).all_coeffs())]   # ascending
            mom = [mp.quad(lambda e, k=k: e ** k * mp.legendre(l, e) * mp.exp(-be * e), mp.linspace(-1, 1, 5))
                   for k in range(kmax + 1)]
            g = {}
            for (j, k), c in zip(Wp.monoms(), Wp.coeffs()):
                g[j] = g.get(j, mp.mpf(0)) + mp.mpf(c) * mom[k]
            gl = [g.get(j, mp.mpf(0)) for j in range(max(g) + 1)]
            q = [mp.mpf(0)] * (len(Pl) + len(gl) - 1)                    # q = P_l g_l
            for i, x in enumerate(Pl):
                for j, y in enumerate(gl):
                    q[i + j] += x * y
            h = [mp.mpf(0)] * len(q)                                     # d/dx[-h e^{-al x}] = q e^{-al x}
            d = list(q)
            j = 0
            while any(d):
                for i, x in enumerate(d):
                    h[i] += x / al ** (j + 1)
                d = [i * d[i] for i in range(1, len(d))]
                j += 1
            h1 = mp.polyval(h[::-1], 1) * mp.exp(-al)

            def outer(x):
                inner = h1 - mp.polyval(h[::-1], x) * mp.exp(-al * x)
                return mp.legenq(l, 0, x, type=3).real * mp.polyval(gl[::-1], x) * mp.exp(-al * x) * inner

            term = (2 * l + 1) * 2 * mp.quad(outer, [1, mp.mpf("1.25"), 2, 4, 8, 16, 40, mp.inf])
            per.append(term)
            total += term
            if l > kmax + 2 and abs(term) < mp.mpf(10) ** (-(dps + min(8, guard // 3))) * abs(total):
                break
        if case == "f0":
            pref = Rm * (2 / Rm) * (Rm / 2) ** 6 * (2 / Rm) ** 4 / 4   # R/(4pi)^2 (2/R)(2pi)^2 (R/2)^6 (2/R)^4
        else:
            pref = mp.mpf(sp.N(N2pi2, dps + 30)) / mp.pi ** 2 * (2 / Rm) * (Rm / 2) ** 6 * (2 * mp.pi) ** 2
        return +(pref * total), per


# ------------------------------------------------------------------ the b -> a limit (a = 1, b = t = 1 + eps)
LIM_SYMBOLS = sp.symbols("R e2 Lg X2 X4 J1 J2 J3")     # e2 = e^{-2R}, Lg = gamma + ln R, X2 = E1(2R),
#                                                        X4 = e^{2R} E1(4R); J1..J3 the integrals outside the class


def limit_b_to_a(cft: Dict[str, Dict[str, str]]) -> dict:
    """Expand K(R; 1, 1 + eps) = sum_i c_i(1 + eps, R) g_i(R; 1 + eps) in eps, c_i given as stored strings.

    The elementary and E1 pieces by Taylor's theorem; the psi brackets by Ei(y) = gamma + ln|y| + sum y^k/(k k!)
    (ln eps cancels inside each bracket); f0 = sum_k eps^k I_k, I_k = int_R^oo FM1_k(s)/s ds with FM1_k the eps^k
    coefficient of label 60, each I_k reduced by parts to the class
        {R^j e2, R^j e2 Lg, X2, R^j X4}  plus  J1 = int e^{-2s} Lg/s,  J2 = int e^{2s} E1(4s)/s,  J3 = int E1(2s)/s.
    Returns {"max_pole", "poles": {m: coefficient of eps^-m}, "limit": eps^0 term (sympy, in LIM_SYMBOLS),
             "J_in_Ik": J-content of every I_k}."""
    R, e2, Lg, X2, X4, J1, J2, J3 = LIM_SYMBOLS
    eps = sp.Symbol("epsilon", positive=True)
    T = sp.Symbol("t")
    cf_t = {nm: sum(sp.sympify(v, locals={"t": T}) * R ** int(p) for p, v in c.items()) for nm, c in cft.items()}
    maxpole = max(dict(sp.factor_list(sp.denom(sp.factor(sp.sympify(v, locals={"t": T}))))[1]).get(T - 1, 0)
                  for c in cft.values() for v in c.values())
    nord = maxpole + 2

    def ser(e, n=nord):
        return sp.series(e, eps, 0, n).removeO()

    x = sp.Symbol("x", positive=True)
    E1d = [None]
    g = -sp.exp(-x) / x
    for _ in range(1, nord + 2):
        E1d.append(g)
        g = sp.diff(g, x)

    def E1_shift(x0, h):                # E1(x0 + h) - E1(x0)
        return sum(h ** k / sp.factorial(k) * E1d[k].subs(x, x0) for k in range(1, nord + 1))

    def Ssum(y):
        return sum(y ** k / (k * sp.factorial(k)) for k in range(1, nord + 1))

    ER = sp.exp(-2 * R)
    fun = {
        "-E1((a+b)R)": -(X2 + E1_shift(2 * R, eps * R)),
        "e^{-(a+b)R}": ER * sp.exp(-eps * R),
        "e^{-2aR}": ER,
        "e^{-2bR}": ER * sp.exp(-2 * eps * R),
        "-e^{2aR}E1(2(a+b)R)": -(X4 + sp.exp(2 * R) * E1_shift(4 * R, 2 * eps * R)),
        "-e^{2bR}E1(2(a+b)R)": -sp.exp(2 * eps * R) * (X4 + sp.exp(2 * R) * E1_shift(4 * R, 2 * eps * R)),
        "psi_a": ER * (-Lg - sp.log(1 + eps / 2) + Ssum(-2 * eps * R) - 2 * Ssum(-eps * R)),
        "psi_b": ER * sp.exp(-2 * eps * R) * (-Lg - sp.log(1 + eps / 2) + Ssum(2 * eps * R) - 2 * Ssum(eps * R)),
    }
    fun = {k: ser(sp.expand(v)) for k, v in fun.items()}
    c = 1 / (2 * (1 - (1 + eps) ** 2))           # label 60 prefactor 1/(2(a^2 - b^2))
    FM1 = c * (fun["psi_a"] - fun["psi_b"] + fun["-e^{2bR}E1(2(a+b)R)"] - fun["-e^{2aR}E1(2(a+b)R)"])
    FM1 = ser(sp.expand(FM1), nord - 1)

    def to_class(e):
        e = sp.expand(e)
        e = sp.expand(e.subs(sp.exp(-2 * R), e2).subs(sp.exp(2 * R), 1 / e2))
        assert not e.has(sp.exp), e.atoms(sp.exp)
        return e

    @lru_cache(maxsize=None)
    def I_e(j):                         # int_R^oo s^j e^{-2s} ds
        if j >= 0:
            return sp.expand(sum(sp.factorial(j) / sp.factorial(i) * R ** i / sp.Integer(2) ** (j - i + 1)
                                 for i in range(j + 1)) * e2)
        if j == -1:
            return X2
        return sp.expand(-R ** (j + 1) * e2 / (j + 1) + sp.Rational(2, j + 1) * I_e(j + 1))

    @lru_cache(maxsize=None)
    def I_eL(j):                        # int_R^oo s^j e^{-2s} Lg(s) ds
        if j >= 0:
            P = sp.expand(I_e(j) / e2)
            Ps = sp.Poly(P, R)
            rest = sum(cf * I_e(m[0] - 1) for m, cf in zip(Ps.monoms(), Ps.coeffs()))
            return sp.expand(P * e2 * Lg + rest)
        if j == -1:
            return J1
        return sp.expand(-R ** (j + 1) * e2 * Lg / (j + 1) + sp.Rational(2, j + 1) * I_eL(j + 1)
                         - sp.Rational(1, j + 1) * I_e(j))

    @lru_cache(maxsize=None)
    def I_X4(j):                        # int_R^oo s^j X4(s) ds,  X4' = 2 X4 - 1/s
        if j >= 0:
            return sp.expand(-R ** j * X4 / 2 - (sp.Rational(j, 2) * I_X4(j - 1) if j > 0 else 0)
                             + sp.Rational(1, 2) * I_e(j - 1))
        if j == -1:
            return J2
        return sp.expand(-R ** (j + 1) * X4 / (j + 1) - sp.Rational(2, j + 1) * I_X4(j + 1)
                         + sp.Rational(1, j + 1) * I_e(j))

    @lru_cache(maxsize=None)
    def I_X2(j):                        # int_R^oo s^j X2(s) ds,  X2' = -e^{-2s}/s
        if j == -1:
            return J3
        return sp.expand(-R ** (j + 1) * X2 / (j + 1) + I_e(j) / (j + 1))

    table = {e2: I_e, e2 * Lg: I_eL, X4: I_X4, X2: I_X2}

    def integrate_over_s(G):            # int_R^oo G(s)/s ds
        out = 0
        for term in sp.Add.make_args(sp.expand(G / R)):
            cR, rest = term.as_independent(e2, Lg, X4, X2)
            cc, j = cR.as_coeff_exponent(R)
            assert not cc.has(R), term
            if rest not in table:
                raise ValueError(("unexpected", term))
            out += cc * table[rest](int(j))
        return sp.expand(out)

    # the reduction rules, by differentiation (d/dR of each = -integrand), j = -4..4
    rep = {e2: sp.exp(-2 * R), Lg: sp.EulerGamma + sp.log(R), X2: sp.expint(1, 2 * R),
           X4: sp.exp(2 * R) * sp.expint(1, 4 * R)}
    for j in range(-4, 5):
        for fn, ig in ((I_e, e2), (I_eL, e2 * Lg), (I_X4, X4), (I_X2, X2)):
            v = fn(j)
            if v.has(J1, J2, J3):
                continue
            assert sp.simplify(sp.diff(v.subs(rep), R) + (R ** j * ig).subs(rep)) == 0, (fn.__name__, j)

    FM1c = to_class(FM1)
    Ik = [integrate_over_s(sp.expand(FM1c.coeff(eps, k))) for k in range(nord - 1)]
    f0exp = sum(eps ** k * Ik[k] for k in range(nord - 1))
    total = 0
    for nm, cfn in cf_t.items():
        cser = sp.series(cfn.subs(T, 1 + eps), eps, 0, 1).removeO()
        gi = f0exp if nm == "f0" else to_class(fun[nm])
        total += sp.expand(cser * gi)
    total = sp.expand(total)
    return {"max_pole": int(maxpole),
            "poles": {m: sp.simplify(total.coeff(eps, -m)) for m in range(1, maxpole + 2)},
            "limit": sp.expand(total.coeff(eps, 0)),
            "J_in_Ik": [[sp.factor(Ik[k].coeff(J)) for J in (J1, J2, J3)] for k in range(nord - 1)]}


def four_coords(e):
    """Coordinates (rational functions of R) of e on (e2 Lg, e2, X2, X4), or raise if e is not in the span."""
    R, e2, Lg, X2, X4, J1, J2, J3 = LIM_SYMBOLS
    e = sp.expand(e)
    cL = e.coeff(e2 * Lg)
    rest = sp.expand(e - cL * e2 * Lg)
    out = [cL, rest.coeff(e2), rest.coeff(X2), rest.coeff(X4)]
    left = sp.expand(rest - out[1] * e2 - out[2] * X2 - out[3] * X4)
    if left != 0 or any(c.free_symbols - {R} for c in out):
        raise ValueError("not in the four-function span: %s" % left)
    return out


def four_dcoords(c4):
    """d/dR in the coordinates of `four_coords`:  e2' = -2 e2,  Lg' = 1/R,  X2' = -e2/R,  X4' = 2 X4 - e2/R."""
    R = LIM_SYMBOLS[0]
    cL, ce, c2, c4_ = c4
    return [sp.diff(cL, R) - 2 * cL, sp.diff(ce, R) - 2 * ce + cL / R - c2 / R - c4_ / R,
            sp.diff(c2, R), sp.diff(c4_, R) + 2 * c4_]


def four_operators(lim, order: int, degree: int) -> List[List[List[int]]]:
    """All operators sum_i c_i(R) D^i, order <= `order`, deg c_i <= `degree`, annihilating `lim` (an element of the
    four-function span), as primitive ascending integer coefficient lists."""
    R = LIM_SYMBOLS[0]
    ders = [four_coords(lim)]
    for _ in range(order):
        ders.append([sp.expand(x) for x in four_dcoords(ders[-1])])
    ncol = (order + 1) * (degree + 1)
    rows: Dict[tuple, List[Fraction]] = {}
    for i in range(order + 1):
        for bi in range(4):
            for term in sp.Add.make_args(sp.expand(ders[i][bi])):
                cc, p = term.as_coeff_exponent(R)
                cc = sp.Rational(cc)
                if cc == 0:
                    continue
                for j in range(degree + 1):
                    rows.setdefault((bi, int(p) + j), [Fraction(0)] * ncol)[i * (degree + 1) + j] += \
                        Fraction(int(cc.p), int(cc.q))
    out = []
    from math import gcd
    for v in PX.nullspace_q([rows[k] for k in sorted(rows)]):
        g = 0
        for x_ in v:
            g = gcd(g, x_)
        v = [x_ // g for x_ in v]
        out.append([PX.trim(v[i * (degree + 1):(i + 1) * (degree + 1)]) for i in range(order + 1)])
    return out


def parse_operator(strings: Sequence[str]) -> List[List[int]]:
    """Frozen operator (list of polynomial strings in R, p_0..p_r, integer coefficients) -> ascending int lists."""
    Rs = sp.Symbol("R")
    out = []
    for s in strings:
        P = sp.Poly(sp.sympify(s, locals={"R": Rs}), Rs)
        out.append([int(c) for c in reversed(P.all_coeffs())])
    return out


def primitive(op: Sequence[Sequence[int]]) -> List[List[int]]:
    """Content removed, leading coefficient of the top-order coefficient positive."""
    from math import gcd
    g = 0
    for c in op:
        for x in c:
            g = gcd(g, int(x))
    op = [PX.trim([int(x) // g for x in c]) for c in op]
    if op[-1][-1] < 0:
        op = [[-x for x in c] for c in op]
    return op
