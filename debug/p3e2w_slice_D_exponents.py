"""P3e-3 step D: local exponents of an operator sum_i t^{2i} Q_i(theta) at t = 0, 1, 2 and infinity
(indicial polynomials), for L_G and M6 = L_G/(theta+2)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp

u, t, rho, s = sp.symbols("u t rho s")
HERE = Path(__file__).parent


def theta_coeffs(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    L = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            L[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in L]


def indicial_at(Lc, t0):
    """apply sum_k c_k(t) theta^k to (t - t0)^rho; leading power of (t - t0)."""
    f = sp.Integer(1)   # represent g(t) (t-t0)^rho by g; theta[(t-t0)^rho g] = (t-t0)^rho [t g' + rho t g/(t-t0)]
    acc = 0
    cur = f
    for c in Lc:
        acc += c * cur
        cur = sp.together(t * sp.diff(cur, t) + rho * t * cur / (t - t0))
    expr = sp.together(acc)
    num, den = sp.fraction(expr)
    num = sp.expand(num.subs(t, t0 + s))
    den = sp.expand(den.subs(t, t0 + s))
    pn = sp.Poly(num, s)
    pd = sp.Poly(den, s)
    ln = min(m[0] for m in pn.monoms())
    ld = min(m[0] for m in pd.monoms())
    return sp.factor(pn.coeff_monomial(s**ln) / pd.coeff_monomial(s**ld)), ln - ld


def main():
    A = json.loads((HERE / "data" / "p3e2w_slice_A.json").read_text())
    Bd = json.loads((HERE / "data" / "p3e2w_slice_B.json").read_text())
    out = {}
    for name, Q in [("L_G", [sp.sympify(q) for q in A["QG"]]), ("M6", [sp.sympify(q) for q in Bd["M6"]])]:
        Lc = theta_coeffs(Q)
        print(name, "leading:", sp.factor(Lc[-1]), flush=True)
        rec = {"t=0": str(sp.factor(Q[0].subs(u, rho)))}
        print("  t=0:", rec["t=0"], flush=True)
        for t0 in (1, -1, 2, -2):
            ip, _ = indicial_at(Lc, t0)
            rec[f"t={t0}"] = str(ip)
            print(f"  t={t0}:", ip, flush=True)
        rec["t=inf (from top Q, exponents mu of t^mu)"] = str(sp.factor(Q[-1].subs(u, rho)))
        print("  t=inf: Q_top(rho) =", rec["t=inf (from top Q, exponents mu of t^mu)"], " order", len(Lc) - 1, flush=True)
        out[name] = rec
    (HERE / "data" / "p3e2w_slice_D.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
