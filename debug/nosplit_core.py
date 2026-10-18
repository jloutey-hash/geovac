"""nosplit core: exact operator algebra in x = t^2, delta = x d/dx (theta = t d/dt = 2 delta).
Operators: lists [c_0, c_1, ...] of sympy rational functions of x, meaning sum_k c_k(x) delta^k
(functions on the left).  The chain's operators from D10 (R3), D7 (R1), theta - 2 are rewritten here."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

x, r, eps = sp.symbols("x r eps")
t = sp.Symbol("t")
HERE = Path(__file__).parent


def dl(e):
    return sp.cancel(x * sp.diff(e, x))


def compose(A, B):
    """A o B."""
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = dl(bi)
    return [sp.cancel(c) for c in out]


def rdiv(L, R):
    """L = Q o R + rem, ord rem < ord R."""
    L = [sp.cancel(c) for c in L]
    m = len(R) - 1
    if len(L) - 1 < m:
        return [sp.Integer(0)], L + [sp.Integer(0)] * (m - len(L))
    Q = [sp.Integer(0)] * (len(L) - m)
    for k in range(len(L) - 1, m - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        Q[k - m] = c
        sub = compose([0] * (k - m) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return Q, [sp.cancel(c) for c in L[:m]]


def theta_to_delta(cs_t):
    """theta-form operator with coefficients even in t -> delta-form in x."""
    out = []
    for k, c in enumerate(cs_t):
        cx = sp.cancel(sp.sympify(c).subs(t, sp.sqrt(x)))
        out.append(sp.cancel(cx * 2**k))
    return out


def chain_ops():
    R3t = [sp.sympify(s, locals={"t": t}) for s in json.loads((HERE / "data" / "density_loop_D10.json").read_text())["R3"]]
    R3 = theta_to_delta(R3t)
    R1 = [-(3 + x), 2 * (1 - x)]            # (theta-3) - t^2 (theta+1)
    T2 = [sp.Integer(-2), sp.Integer(2)]     # theta - 2
    L3 = [sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(32)]
    L3 = [sp.cancel(a + b) for a, b in zip(L3, compose([-x], [1, 6, 12, 8]))]  # 4 th^3 - t^2 (th+1)^3 = 32 d^3 - x (2d+1)^3
    return {"R3": R3, "R1": R1, "T2": T2, "L3": L3}


def indicial(L, p):
    """Indicial polynomial I(r) of L at p (rational number or 'oo'), and whether L is regular
    singular (or ordinary) there: deg I == ord L.  Exponents in the local variable x - p,
    x (p = 0) or 1/x (p = oo)."""
    n = len(L) - 1
    terms = []
    if p == "oo":
        w = sp.Symbol("w")
        for k, c in enumerate(L):
            terms.append(sp.cancel(c.subs(x, 1 / w)) * (-r) ** k)
        var = w
    elif p == 0:
        for k, c in enumerate(L):
            terms.append(c * r**k)
        var = x
    else:
        psi = [sp.Integer(1)]
        for k in range(n):
            psi.append(sp.cancel((p + eps) * ((r - k) * psi[-1] + eps * sp.diff(psi[-1], eps))))
        for k, c in enumerate(L):
            terms.append(sp.cancel(c.subs(x, p + eps) * eps ** (-k) * psi[k]))
        var = eps
    tot = sp.cancel(sp.Add(*terms))
    num, den = sp.fraction(sp.factor(tot))
    m = min(sp.Poly(sp.expand(num), var).monoms()[-1][0], 10**9) - sp.Poly(den, var).monoms()[-1][0] \
        if False else None
    # leading coefficient of the Laurent expansion in var
    ser = sp.series(tot, var, 0, 1).removeO() if False else None
    nd = sp.Poly(sp.expand(num), var)
    dd = sp.Poly(sp.expand(den), var)
    lo_n = min(e[0] for e in nd.monoms())
    lo_d = min(e[0] for e in dd.monoms())
    I = sp.expand(nd.as_expr().coeff(var, lo_n) / dd.as_expr().coeff(var, lo_d))
    # nd's lowest coefficient may be a polynomial in r; ensure the ratio is polynomial in r
    I = sp.factor(sp.cancel(I))
    deg = sp.Poly(I, r).degree()
    return I, deg == n


def exponents(L, p):
    I, reg = indicial(L, p)
    roots = sp.roots(sp.Poly(I, r))
    if sum(roots.values()) != sp.Poly(I, r).degree():
        roots = {rt: 1 for rt in sp.Poly(I, r).all_roots()}
    return I, reg, roots


def singular_points(L):
    """Finite singular points of L (zeros of the monic coefficients' denominators), as
    irreducible factors over Q, plus 0 and oo."""
    Lm = [sp.cancel(c / L[-1]) for c in L]
    den = sp.lcm([sp.denom(sp.together(c)) for c in Lm])
    fac = sp.factor_list(sp.Poly(den, x))[1]
    return [f.as_expr() for f, _ in fac]
