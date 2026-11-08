"""MG1-MG3 (copy of cf_g2_8b_order1.py, parametrised): the order-one composition factors c1(x) theta + c0(x):
residues of d/dx log s = -c0/(2 x c1) at every root of x c1 (exact, sympy), the closed-form solution s and a
check that the operator kills it.
usage: python mg_8b_order1.py <pieces json> <out tag>"""
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
x, t = sp.symbols("x t")
d = json.loads(Path(sys.argv[1]).read_text())
out = {}
for name in [k for k, v in d.items() if k.startswith("factor_") and isinstance(v, list) and len(v) == 2]:
    c0, c1 = [sp.Poly(sp.sympify(s.replace("^", "**"), locals={"x": x}), x).as_expr() for s in d[name]]
    r = sp.cancel(-c0 / (2 * x * c1))          # d/dx log s
    num, den = sp.fraction(r)
    res = {}
    s_expr = 1
    for f, m in sp.factor_list(den)[1]:
        assert m == 1, (name, f, m)
        q = sp.rem(sp.expand(num * sp.invert(sp.diff(den, x), f, x)), f, x)
        q = sp.simplify(q)
        res[str(f)] = str(q)
        if q.is_number:
            s_expr *= f ** q
    pol_part = sp.cancel(r - sum(sp.Rational(sp.sympify(res[str(f)])) * sp.diff(f, x) / f for f, m in sp.factor_list(den)[1]))
    s_t = sp.simplify(s_expr.subs(x, t**2))
    chk = sp.simplify(c1.subs(x, t**2) * t * sp.diff(s_t, t) + c0.subs(x, t**2) * s_t)
    out[name] = {"dlog_residues": res, "polynomial_part_of_dlog": str(pol_part), "solution": str(sp.factor(s_expr)),
                 "operator_kills_solution": chk == 0}
    print(name, out[name], flush=True)
(Path(__file__).parent / "data" / f"mg_8b_order1_{sys.argv[2]}.json").write_text(json.dumps(out, indent=1))
