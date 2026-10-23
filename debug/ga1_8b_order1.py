"""GA1: the order-one composition factors c1(x) theta + c0(x): residues of d/dx log s = -c0/(2 x c1) at every
root of x c1 (exact, sympy; for each irreducible factor f the residue reduced mod f must be a constant);
then the closed-form solution s and a check that the operator kills it.  Copy of debug/cf_g2_8b_order1.py
made generic over the factor list.  A solution is algebraic iff every residue is rational and the
polynomial part of d log s vanishes.
usage: python ga1_8b_order1.py <pieces json> <out tag>"""
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
(Path(__file__).parent / "data" / f"ga1_8b_order1_{tag}.json").write_text(json.dumps(out, indent=1))
