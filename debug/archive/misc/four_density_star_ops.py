"""Exact operator toolkit for the FD1-FD3 star computation (shared by steps F..).

Operators are function-left lists [c_0, c_1, ...] meaning sum_k c_k(w) delta^k, with w = t^2,
delta = w d/dw = theta_t / 2, and c_k in the rational function field QQ(w) (sympy sparse fields,
much faster than generic cancel).  Q-form sum_i t^(2i) Q_i(theta_t) maps to sum_i w^i Q_i(2 delta)
(functions left: w^i Q_i(2 delta) is already function-left).
"""
from __future__ import annotations

from fractions import Fraction
from typing import List

import sympy as sp
from sympy import QQ
from sympy.polys.fields import field

K, W = field("w", QQ)
u = sp.Symbol("u")
w_sym = sp.Symbol("w")
t_sym = sp.Symbol("t")

Op = List


def d(c):
    """delta c = w dc/dw."""
    return W * c.diff(W)


def from_Q(Q) -> Op:
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    L = [K(0)] * n
    for i, q in enumerate(Q):
        p = sp.Poly(sp.expand(sp.sympify(q).subs(u, 2 * u)), u)
        for (k,), c in p.terms():
            L[k] += K(sp.Rational(c)) * W**i
    return L


def from_theta_t(Lt: List[sp.Expr]) -> Op:
    """function-left sum_k c_k(t) theta_t^k with c_k even rational in t -> delta form."""
    n = len(Lt)
    out = [K(0)] * n
    # theta_t^k = 2^k delta^k
    for k, c in enumerate(Lt):
        e = sp.together(sp.sympify(c))
        e = sp.cancel(e)
        num, den = sp.fraction(e)
        num = sp.Poly(num, t_sym)
        den = sp.Poly(den, t_sym)
        def to_w(p):
            acc = K(0)
            for (e_,), cf in p.terms():
                assert e_ % 2 == 0, "odd power of t"
                acc += K(sp.Rational(cf)) * W ** (e_ // 2)
            return acc
        out[k] = to_w(num) / to_w(den) * 2**k
    return out


def to_sympy(L: Op, var=w_sym):
    return [c.as_expr().subs(w_sym, var) for c in L]


def trim(L: Op) -> Op:
    L = list(L)
    while len(L) > 1 and L[-1] == 0:
        L.pop()
    return L


def compose(A: Op, B: Op) -> Op:
    out = [K(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        if a == 0:
            continue
        for j, b in enumerate(B):
            if b == 0:
                continue
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = d(bi)
    return trim(out)


def add(A: Op, B: Op) -> Op:
    n = max(len(A), len(B))
    return trim([(A[i] if i < len(A) else K(0)) + (B[i] if i < len(B) else K(0)) for i in range(n)])


def scale(c, A: Op) -> Op:
    return [c * a for a in A]


def rdiv(L: Op, R: Op):
    """L = Q o R + rem, ord rem < ord R."""
    L = list(L)
    r = len(R) - 1
    Q = [K(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        if L[k] == 0:
            continue
        c = L[k] / R[-1]
        Q[k - r] = c
        sub = compose([K(0)] * (k - r) + [c], R)
        for i in range(len(sub)):
            L[i] -= sub[i]
    return Q, L[:r]


def ldiv(L: Op, Lft: Op):
    """L = Lft o Q + rem? (left division) via adjoints."""
    Q, rem = rdiv(adjoint(L), adjoint(Lft))
    return adjoint(Q), rem


def adjoint(L: Op) -> Op:
    """formal adjoint w.r.t. dw/w: delta* = -delta, (c delta^k)* = (-delta)^k c."""
    out = [K(0)] * len(L)
    for k, c in enumerate(L):
        term = compose([K(0)] * k + [K((-1) ** k)], [c])
        for i, x in enumerate(term):
            out[i] += x
    return trim(out)


def to_Q(L: Op):
    """clear denominators on the left; return (Q, D): sum_i t^(2i) Q_i(theta_t) = D(w) * L."""
    D = None
    for c in L:
        if c != 0:
            D = c.denom if D is None else D.lcm(c.denom)
    Q = {}
    for k, c in enumerate(L):
        if c == 0:
            continue
        e = c * K(D)
        num, den = e.numer, e.denom   # den must be a constant
        assert den.is_ground, "denominator not cleared"
        dc = sp.Rational(den.LC)
        for (i,), cf in num.terms():
            Q[i] = Q.get(i, 0) + sp.Rational(cf) / dc * (u / 2) ** k
    return [sp.expand(Q.get(i, 0)) for i in range(max(Q) + 1)], D.as_expr()


def apply_fn(L: Op, f: sp.Expr, var=w_sym):
    """apply to an explicit sympy function of w (var)."""
    acc, g = 0, f
    for c in L:
        acc += c.as_expr().subs(w_sym, var) * g
        g = var * sp.diff(g, var)
    return sp.simplify(acc)
