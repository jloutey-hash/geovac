"""Operator utilities for test QG (sympy): theta-form <-> D-form, local exponents, composition."""
from __future__ import annotations
import json
import sympy as sp

u, z, y, a = sp.symbols("u z y a")


def load_Q(path, r):
    d = json.loads(open(path).read())[str(r)]
    P, dd = d["P"], d["d"]
    return [sp.expand(sum(sp.Integer(int(v)) * u ** int(k.split(",")[1]) for k, v in P.items()
                          if int(k.split(",")[0]) == i)) for i in range(dd + 1)]


def theta_to_D(Q, var=z, step=1):
    """sum_i var^(step i) Q_i(theta) -> [a_0..a_r] with L = sum a_k(var) d^k/dvar^k."""
    from sympy.functions.combinatorial.numbers import stirling
    r = max(sp.Poly(q, u).degree() for q in Q)
    A = [sp.Integer(0)] * (r + 1)
    for i, q in enumerate(Q):
        for (j,), c in sp.Poly(q, u).terms():
            for k in range(j + 1):
                A[k] += c * var ** (step * i) * stirling(j, k, kind=2) * var ** k
    return [sp.expand(x) for x in A]


def indicial(A, z0, var=z):
    """indicial polynomial in a at var = z0 (finite) for L = sum A_k d^k."""
    terms = []
    for k, ak in enumerate(A):
        if ak == 0:
            continue
        p = sp.Poly(sp.expand(ak.subs(var, z0 + y)), y)
        o = min(m[0] for m in p.monoms())
        terms.append((o - k, k, p.coeff_monomial(y ** o)))
    mn = min(t[0] for t in terms)
    ff = lambda kk: sp.prod([a - i for i in range(kk)])
    return sp.factor(sum(c * ff(k) for o, k, c in terms if o == mn))


def indicial_inf(A, var=z):
    """exponents rho with f ~ var^(-rho) at infinity."""
    terms = []
    for k, ak in enumerate(A):
        if ak == 0:
            continue
        p = sp.Poly(ak, var)
        o = p.degree()
        terms.append((o - k, k, p.LC()))
    mx = max(t[0] for t in terms)
    ff = lambda kk: sp.prod([(-a) - i for i in range(kk)])
    return sp.factor(sum(c * ff(k) for o, k, c in terms if o == mx))


def mellin_Y_from_phi(Qphi):
    """Y(t) = C t^3 int_0^inf z^(1/2) (z + t^2)^(-3/2) phi(z) dz.  From L_phi = sum_i z^i P_i(theta_z):
    L_Y = sum_i t^(2i) P_i(T/2 - 3/2) (T/2)_i prod_{m=0}^{d-i-1} (m + 3/2 - T/2),   T = theta_t.
    (Mellin: M[S](s) = B(s, 3/2 - s) M[phi](s), S(x) = x^(-3/2) Y / C, x = t^2.)  Returned as
    sum_i t^(2i) QY_i(u) with integer coefficients, u = theta_t."""
    d = len(Qphi) - 1
    out = []
    for i, P in enumerate(Qphi):
        e = P.subs(u, u / 2 - sp.Rational(3, 2))
        e *= sp.prod([u / 2 + m for m in range(i)])
        e *= sp.prod([m + sp.Rational(3, 2) - u / 2 for m in range(d - i)])
        out.append(sp.expand(e))
    den = sp.lcm([sp.denom(c) for q in out for c in sp.Poly(q, u).coeffs()])
    out = [sp.expand(q * den) for q in out]
    g = sp.gcd_list([c for q in out for c in sp.Poly(q, u).coeffs()])
    return [sp.expand(q / g) for q in out]


def compose_theta(A, B, step=2, var=None):
    """(sum_i var^(step i) A_i(theta)) o (sum_j var^(step j) B_j(theta)) = sum var^(step(i+j)) A_i(theta) B_j(theta + step i)."""
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for i, Ai in enumerate(A):
        for j, Bj in enumerate(B):
            out[i + j] += sp.expand(Ai * Bj.subs(u, u + step * i))
    return [sp.expand(x) for x in out]


# ---- operators in theta-form with function-left rational coefficients: L = sum_k c_k(t) theta^k
tt = sp.Symbol("t")


def th_fn(c):
    return sp.cancel(tt * sp.diff(c, tt))


def Q_to_coeffs(Q, step=2):
    """sum_i t^(step i) Q_i(theta) -> [c_0..c_r] function-left."""
    r = max(sp.Poly(q, u).degree() for q in Q)
    c = [sp.Integer(0)] * (r + 1)
    for i, q in enumerate(Q):
        for (j,), v in sp.Poly(q, u).terms():
            c[j] += v * tt ** (step * i)
    return [sp.expand(x) for x in c]


def comp(A, B):
    """(sum a_k theta^k) o (sum b_j theta^j), function-left."""
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        if a == 0:
            continue
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = th_fn(bi)
    return [sp.cancel(x) for x in out]


def rdiv(L, R):
    """L = Q o R + rem, ord rem < ord R."""
    L = [sp.cancel(c) for c in L]
    while len(L) > 1 and L[-1] == 0:
        L = L[:-1]
    r = len(R) - 1
    if len(L) - 1 < r:
        return [sp.Integer(0)], L
    Qo = [sp.Integer(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        if c == 0:
            continue
        Qo[k - r] = c
        sub = comp([0] * (k - r) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    rem = L[:r]
    while len(rem) > 1 and rem[-1] == 0:
        rem = rem[:-1]
    return Qo, rem


def monic(L):
    return [sp.cancel(c / L[-1]) for c in L]


def gcrd(A, B):
    A, B = monic(A), monic(B)
    while not (len(B) == 1 and B[0] == 0):
        _, R = rdiv(A, B)
        if all(c == 0 for c in R):
            return B
        A, B = B, monic(R)
    return A


def apply_fn(L, f):
    """L applied to an explicit sympy function f(t)."""
    acc, g = 0, f
    for c in L:
        acc += c * g
        g = tt * sp.diff(g, tt)
    return sp.simplify(acc)


def coeffs_to_Q(L, step=2):
    """function-left rational theta-coefficients -> sum_i t^(step i) Q_i(u) (cleared on the left)."""
    den = sp.lcm([sp.denom(sp.together(c)) for c in L])
    cs = [sp.expand(sp.cancel(c * den)) for c in L]
    Q = {}
    for k, c in enumerate(cs):
        pc = sp.Poly(c, tt)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            assert e % step == 0, (e, step)
            Q[e // step] = Q.get(e // step, 0) + v * u ** k
    Q = [sp.expand(Q.get(i, 0)) for i in range(max(Q) + 1)]
    g = sp.gcd_list([c for q in Q for c in sp.Poly(q, u).coeffs()])
    return [sp.expand(q / g) for q in Q]


def Q_to_P(Q):
    P = {}
    for i, q in enumerate(Q):
        for (j,), v in sp.Poly(q, u).terms():
            P[(i, j)] = int(v)
    return P
