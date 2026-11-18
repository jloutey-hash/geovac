"""Pachucki-ODE route to the two-centre 1s exchange integral (derivation instrument).

Source: K. Pachucki, PRA 80, 032520 (2009) [arXiv:0905.1635]:
  master  f(r) = int d3r1/4pi d3r2/4pi e^{-u3 r1A - u2 r1B - w2 r2A - w3 r2B} r/(r1A r1B r2A r2B r12)  (Eq. 1)
  ODE     r f'' + f' - p^2 r f + F(r) = 0                                                    (Eq. 29)
  p^2     Eq. (21);  F = (1/2)[a+ F1- + b+ F2- + a- F3- + b- F4-]   Eqs. (23,24,31-35).

Exchange point theta0: u3 = w2 = a, u2 = w3 = b  (Sec. VII; p^2 = 0 there).
Map:  K(R) = (1s_A^a 1s_B^b | 1s_A^a 1s_B^b) = 16 a^3 b^3 / R * d_u3 d_u2 d_w2 d_w3 f |_theta0.

Derivation used here (no guessing):
  differentiate Eq. (29) by a multi-index S of exponents and set theta = theta0 (p^2 = 0):
      (r f_S')' = -F_S + r * sum_{0 != U <= S} binom(S,U) (d_U p^2) f_{S-U}           (*)
  f_S decays at r -> oo (the integrand does), the kernel of (r d)' is {1, ln r} (non-decaying),
  so f_S is the UNIQUE decaying solution of (*).  We solve (*) exactly by an ansatz in a finite
  set of decaying functions with Laurent-polynomial coefficients; if the linear system has a
  solution, it IS f_S.  The only non-elementary function allowed is f0 = f(theta0) itself,
  with the exact rule  r f0' = -F^{(-1)}  (Pachucki Eq. 62-63; F^{(-1)} = int F, checked here).

Function representation: dict key -> {power k: coeff}, coeff a sympy expression in L3 = ln|.|
constant (polynomial).  Keys:
  ('E', lam)        r^k e^{lam r}
  ('Ei', lam, mu)   r^k e^{lam r} Ei(mu r)
  ('f',)            r^k f0(r)
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import sympy as sp

r = sp.Symbol("r", positive=True)
u2, u3, w2, w3 = sp.symbols("u2 u3 w2 w3", positive=True)
THETA = (u3, u2, w2, w3)          # order: r1A, r1B, r2A, r2B
LC = sp.Symbol("LC")              # the constant ln|...| of Eqs. (34,35) at theta0


def pachucki_F_and_p2():
    Ei = sp.Ei
    E = sp.exp
    S = u2 + u3 + w2 + w3
    # |.| argument of Eqs. (34),(35); its sign near theta0 (a<b) is negative -> use -(...)
    arg = (u2 + w3 - u3 - w2) * (u2 + u3) * (w2 + w3) / (S * (u2 - u3) * (w2 - w3))
    lnabs = sp.log(-arg)
    F1m = E(r * (u3 - w2)) * Ei(-r * (u2 + u3)) - E(r * (w2 - u3)) * Ei(-r * (w2 + w3))
    F2m = E(r * (w3 - u2)) * Ei(-r * (w2 + w3)) - E(r * (u2 - w3)) * Ei(-r * (u2 + u3))
    F3m = (E(-r * (u3 + w2)) * (lnabs - Ei(r * (u3 + w2 - u2 - w3)) + Ei(r * (w2 - w3)) + Ei(r * (u3 - u2)))
           - E(r * (u3 + w2)) * Ei(-r * S))
    F4m = (E(-r * (u2 + w3)) * (lnabs - Ei(r * (u2 + w3 - u3 - w2)) + Ei(r * (w3 - w2)) + Ei(r * (u2 - u3)))
           - E(r * (u2 + w3)) * Ei(-r * S))
    ap = w2 / (w2**2 - w3**2) + u3 / (u2**2 - u3**2)
    am = w2 / (w2**2 - w3**2) - u3 / (u2**2 - u3**2)
    bp = u2 / (u2**2 - u3**2) + w3 / (w2**2 - w3**2)
    bm = u2 / (u2**2 - u3**2) - w3 / (w2**2 - w3**2)
    F = sp.Rational(1, 2) * (ap * F1m + bp * F2m + am * F3m + bm * F4m)
    p2 = ((u2**2 - u3**2 + w2**2 - w3**2) * (u2**2 * w2**2 - u3**2 * w3**2)
          / ((u2**2 - u3**2) * (w2**2 - w3**2)))
    return F, p2, lnabs


def theta0(a, b):
    return {u3: a, w2: a, u2: b, w3: b}


# ------------------------------------------------------------------ function algebra
Func = Dict[tuple, Dict[int, sp.Expr]]


def fadd(A: Func, B: Func, s=1) -> Func:
    out = {k: dict(v) for k, v in A.items()}
    for k, d in B.items():
        o = out.setdefault(k, {})
        for p, c in d.items():
            o[p] = sp.expand(o.get(p, 0) + s * c)
    return clean(out)


def clean(A: Func) -> Func:
    out = {}
    for k, d in A.items():
        dd = {p: c for p, c in d.items() if sp.expand(c) != 0}
        if dd:
            out[k] = dd
    return out


def fscale(A: Func, c, shift=0) -> Func:
    return clean({k: {p + shift: sp.expand(c * v) for p, v in d.items()} for k, d in A.items()})


FM1: Func = {}      # filled with F^{(-1)} once known (derivative rule of f0)


def fdiff(A: Func) -> Func:
    out: Func = {}

    def add(k, p, c):
        o = out.setdefault(k, {})
        o[p] = sp.expand(o.get(p, 0) + c)

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
                add(("E", lam + mu), p - 1, c)
            elif k[0] == "f":
                for kk, dd in FM1.items():          # r^p * (-F^{(-1)}/r)
                    for pp, cc in dd.items():
                        add(kk, p + pp - 1, -c * cc)
    return clean(out)


def L0(A: Func) -> Func:
    """(r A')'"""
    return fdiff(fscale(fdiff(A), 1, 1))


def to_func(expr, lnconst=None) -> Func:
    """Parse a sympy expression in r (sums of coeff * r^k * exp(lam r) * Ei(mu r))."""
    expr = sp.expand(sp.powsimp(sp.expand(expr), combine="exp"))
    if lnconst is not None:
        expr = expr.subs(lnconst[0], LC)
    out: Func = {}
    for t in sp.Add.make_args(expr):
        coeff, k, lam, mu = sp.Integer(1), 0, sp.Integer(0), None
        for fac in sp.Mul.make_args(t):
            if not fac.has(r):
                coeff *= fac
            elif fac == r:
                k += 1
            elif fac.is_Pow and fac.base == r:
                k += int(fac.exp)
            elif isinstance(fac, sp.exp):
                lam += sp.expand(fac.args[0] / r)
            elif fac.is_Pow and isinstance(fac.base, sp.exp):
                lam += sp.expand(fac.base.args[0] * fac.exp / r)
            elif isinstance(fac, sp.Ei):
                assert mu is None
                mu = sp.nsimplify(sp.expand(fac.args[0] / r))
            else:
                raise ValueError(f"unparsed factor {fac} in {t}")
        lam = sp.nsimplify(lam)
        assert not lam.has(r)
        key = ("E", lam) if mu is None else ("Ei", lam, mu)
        o = out.setdefault(key, {})
        o[k] = sp.expand(o.get(k, 0) + coeff)
    return clean(out)


def fequal(A: Func, B: Func) -> bool:
    return not clean(fadd(A, B, -1))


# ------------------------------------------------------------------ ansatz solver
def solve_L0(src: Func, keys: List[tuple], pmin: int, pmax: int) -> Func:
    """Find u = sum_{key in keys, pmin<=k<=pmax} x r^k B_key with L0 u = src (exact)."""
    cols = [(k, p) for k in keys for p in range(pmin, pmax + 1)]
    imgs = [L0({k: {p: sp.Integer(1)}}) for (k, p) in cols]
    rows = set()
    for im in imgs + [src]:
        for k, d in im.items():
            for p in d:
                rows.add((k, p))
    rows = sorted(rows, key=str)
    ridx = {rw: i for i, rw in enumerate(rows)}
    M = sp.zeros(len(rows), len(cols))
    bvec = sp.zeros(len(rows), 1)
    for j, im in enumerate(imgs):
        for k, d in im.items():
            for p, c in d.items():
                M[ridx[(k, p)], j] = c
    for k, d in src.items():
        for p, c in d.items():
            bvec[ridx[(k, p)], 0] = c
    from sympy.polys.matrices import DomainMatrix
    dom = sp.QQ.frac_field(LC)
    Aug = DomainMatrix.from_Matrix(M.row_join(bvec)).convert_to(dom)
    R, piv = Aug.rref()
    ncol = len(cols)
    if ncol in piv:
        raise RuntimeError("no solution in ansatz")
    Rm = R.to_Matrix()
    sol: Func = {}
    for i, pc in enumerate(piv):
        val = sp.factor(Rm[i, ncol])
        if val != 0:
            k, p = cols[pc]
            sol.setdefault(k, {})[p] = sp.expand(val)
    # free columns -> 0 ; check
    if not fequal(L0(sol), src):
        raise RuntimeError("check failed")
    free = [c for j, c in enumerate(cols) if j not in piv]
    return clean(sol), free


def keys_closure(A: Func) -> set:
    ks = set(A.keys())
    while True:
        new = set(fdiff({k: {0: sp.Integer(1)} for k in ks}).keys()) | ks
        new |= set(fdiff({k: {5: sp.Integer(1)} for k in ks}).keys())
        if new == ks:
            return ks
        ks = new


def multi_derivs(F, p2, alphas, sub):
    """dict alpha -> (F_alpha at theta0 as Func, d_alpha p^2 at theta0)."""
    # incremental: d^alpha from a cached parent, one variable at a time (memoised)
    cache = {(0, 0, 0, 0): (F, p2)}

    def get(al):
        if al in cache:
            return cache[al]
        i = max(j for j in range(4) if al[j])
        par = tuple(x - (1 if j == i else 0) for j, x in enumerate(al))
        e, q = get(par)
        cache[al] = (sp.diff(e, THETA[i]), sp.diff(q, THETA[i]))
        return cache[al]

    out = {}
    for al in sorted(alphas, key=sum):
        e, q = get(al)
        out[al] = (e.subs(sub), sp.nsimplify(q.subs(sub)))
    return out


def sub_multi(al):
    return list(itertools.product(*[range(n + 1) for n in al]))


def binom_multi(al, be):
    out = 1
    for x, y in zip(al, be):
        out *= sp.binomial(x, y)
    return out
