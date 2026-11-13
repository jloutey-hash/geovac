"""Score record SA (docs/predictions/2026-11-13_star5_attainment.md): the five-density star reaches rung five.

SA1/SA4: local exponents of phi(L_4) at t = -4 and phi(L_6) at t = -6 (exact, sympy).
SA2: leading coefficient of phi(L_c(theta - lam)), c <= 3, lam in -6..12, nonzero at t = +-4 (exact).
SA3 + control: eps^(1/2) I5(eps) -> pi^2/(8 sqrt 2); eps I4(eps) -> (pi/2)^(3/2)/sqrt(6 pi) (float, scaled Bessels).
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.special import i0e, k1e, roots_laguerre

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import test_paper63_tree_rung as TR  # noqa: E402

t, th, s, rho = sp.symbols("t theta s rho")


def t_operator(P):
    """sum_i t^(2i) Q_i(theta) from the u^k theta^j dict (left-multiplied by t^(2K))."""
    return TR._t_lists(P)


def theta_to_d(Q):
    """sum_i t^(2i) Q_i(theta) -> coefficients a_k(t) of d^k/dt^k (theta^j = sum S(j,k) t^k d^k)."""
    a = {}
    for i, q in enumerate(Q):
        for (j,), cf in sp.Poly(sp.expand(q), th).terms():
            for k in range(j + 1):
                a[k] = a.get(k, 0) + cf * sp.functions.combinatorial.numbers.stirling(j, k) * t ** (2 * i + k)
    return {k: sp.expand(v) for k, v in a.items() if sp.expand(v) != 0}


def exponents_at(P, t0):
    a = theta_to_d(t_operator(P))
    vk = {}
    for k, ak in a.items():
        p = sp.Poly(sp.expand(ak.subs(t, s + t0)), s)
        v = min(m[0] for m in p.monoms())
        vk[k] = (v, p.coeff_monomial(s**v))
    low = min(v - k for k, (v, _) in vk.items())
    ind = sum(c * sp.ff(rho, k) for k, (v, c) in vk.items() if v - k == low)
    return sp.factor(ind), sp.solve(ind, rho)


def leading(P):
    Q = t_operator(P)
    J = max(sp.Poly(q, th).degree() for q in Q)
    return sp.factor(sum(sp.Poly(q, th).coeff_monomial(th**J) * t ** (2 * i) for i, q in enumerate(Q)))


out = {}
L4 = TR._sym_op(4)
P4 = TR._phi(L4)
lead4 = leading(P4)
ind4, ex4 = exponents_at(P4, -4)
out["SA1"] = {"leading": str(lead4), "indicial_at_-4": str(ind4), "exponents": [str(e) for e in ex4]}
print("SA1 leading", lead4, "| exponents at -4:", ex4)

bad = []
for c in range(0, 4):
    Lc = TR._sym_op(c) if c > 0 else {(0, 1): sp.Integer(1)}
    for lam in range(-6, 13):
        lc = leading(TR._phi(TR._shift(Lc, lam)))
        if lc.subs(t, 4) == 0 or lc.subs(t, -4) == 0:
            bad.append((c, lam, str(lc)))
out["SA2_bad"] = bad
print("SA2 zero leading coefficients at +-4:", bad)

P6 = TR._phi(TR._sym_op(6))
ind6, ex6 = exponents_at(P6, -6)
out["SA4"] = {"indicial_at_-6": str(ind6), "exponents": [str(e) for e in ex6]}
print("SA4 exponents at -6:", ex6)

# ---------------------------------------------------------------- SA3 + control (float)
ul, wl = roots_laguerre(80)


def Rs(x):
    """e^x R(x) = int_0^oo (x+u) k1e(x+u) e^-u du (Gauss-Laguerre)."""
    X = x[:, None] + ul[None, :]
    return (X * k1e(X)) @ wl


def branch(eps, kind):
    y = np.linspace(np.log(1.0), np.log(1.0 + 80.0 / eps), 400001)
    x = np.exp(y)
    if kind == 5:
        f = i0e((4 - eps) * x) * Rs(x) ** 4 * np.exp(-eps * x) * x**-2
    else:
        f = i0e((3 - eps) * x) * Rs(x) ** 3 * np.exp(-eps * x) * x**-1
    return np.trapezoid(f * x, y)


C5, C4 = np.pi**2 / (8 * np.sqrt(2)), (np.pi / 2) ** 1.5 / np.sqrt(6 * np.pi)
rows = []
for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
    I5, I4 = branch(eps, 5), branch(eps, 4)
    rows.append({"eps": eps, "sqrt_eps_I5": float(np.sqrt(eps) * I5), "eps_I4": float(eps * I4),
                 "sqrt_eps_I4": float(np.sqrt(eps) * I4)})
    print(f"eps={eps:g}: sqrt(eps) I5 = {np.sqrt(eps)*I5:.5f} (C5 {C5:.5f}) | eps I4 = {eps*I4:.5f} (C4 {C4:.5f})"
          f" | sqrt(eps) I4 = {np.sqrt(eps)*I4:.3f}", flush=True)
out["SA3"] = rows
out["C5"], out["C4"] = C5, C4
Path(__file__).resolve().parent.joinpath("data", "star5_attainment.json").write_text(json.dumps(out, indent=1))
