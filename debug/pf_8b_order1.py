"""PF1-PF4 (copy of debug/ga3_8b_order1.py): closed-form solutions of the order-one composition factors (residues of d log s; algebraic iff rational residues and zero polynomial part).
Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_8b_order1.py <pieces json> <out tag>"""
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
x, t = sp.symbols("x t")
d = json.loads(Path(sys.argv[1]).read_text())
tag = sys.argv[2]
out = {}
for name in sorted([k for k in d if k.startswith("factor_") and isinstance(d[k], list)], key=lambda k: int(k.split("_")[1])):
    if len(d[name]) != 2:
        continue
    c0, c1 = [sp.Poly(sp.sympify(s.replace("^", "**"), locals={"x": x}), x).as_expr() for s in d[name]]
    r = sp.cancel(-c0 / (2 * x * c1))          # d/dx log s
    num, den = sp.fraction(r)
    res = {}
    s_expr = 1
    algebraic = True
    for f, m in sp.factor_list(den)[1]:
        if m != 1:
            algebraic = False
            res[str(f)] = f"pole of order {m}: not a simple pole"
            continue
        q = sp.rem(sp.expand(num * sp.invert(sp.diff(den, x), f, x)), f, x)
        q = sp.simplify(q)
        res[str(f)] = str(q)
        if q.is_number and q.is_rational:
            s_expr *= f ** q
        else:
            algebraic = False
    pol_part = sp.cancel(r - sum(sp.Rational(sp.sympify(res[str(f)])) * sp.diff(f, x) / f for f, m in sp.factor_list(den)[1] if m == 1 and sp.sympify(res[str(f)]).is_rational))
    algebraic = algebraic and pol_part == 0
    s_t = sp.simplify(s_expr.subs(x, t**2))
    chk = sp.simplify(c1.subs(x, t**2) * t * sp.diff(s_t, t) + c0.subs(x, t**2) * s_t)
    out[name] = {"dlog_residues": res, "polynomial_part_of_dlog": str(pol_part), "solution": str(sp.factor(s_expr)),
                 "operator_kills_solution": bool(chk == 0), "algebraic": bool(algebraic)}
    print(name, out[name], flush=True)
(Path(__file__).parent / "data" / f"pf_8b_order1_{tag}.json").write_text(json.dumps(out, indent=1))
