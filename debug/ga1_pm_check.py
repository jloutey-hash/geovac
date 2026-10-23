"""PM recheck of GA1 (docs/predictions/2026-10-22_geminal_shapes.md) on a code path independent of the sub-agent's:
the Legendre operator at gamma = 3/2 rebuilt from c2b_feyn_A.leg_op, moved to x and twisted by x^(1/2) by hand;
the stored order-two factor of Y_lin's annihilator converted to d/dx by hand; the stored map applied with the
sympy residual of tests/test_paper59_c2b_elliptic.py.  Also: the stored factor is not C2b itself.
Usage: python debug/ga1_pm_check.py"""
import ast
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "debug"))
sys.path.insert(0, str(ROOT / "tests"))
from c2b_feyn_A import leg_op, t  # noqa: E402
from test_paper59_c2b_elliptic import _residual  # noqa: E402

x = sp.Symbol("x")
d = json.loads((ROOT / "debug" / "data" / "ga1_pieces.json").read_text())


def poly(s_):
    return sp.expand(sp.sympify(s_.replace("^", "**"), locals={"x": x}))


c0, c1, c2 = [poly(v) for v in d["factor_5"]]
target = [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]           # theta = 2x d/dx
l0, l1, l2 = leg_op(sp.Rational(3, 2))                           # H = t^(-1/2) w(z), d/dt
e2 = sp.expand(sp.cancel(l2 / t).subs(t, sp.sqrt(x)))
e1 = sp.expand(l1.subs(t, sp.sqrt(x)))
e0 = sp.expand(sp.cancel(l0 / t).subs(t, sp.sqrt(x)))
b0, b1, b2 = e0, 2 * e2 + 2 * e1, 4 * x * e2
L = -sp.Rational(1, 2) / x                                       # H = x^(-1/2) y
tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
src = [sp.expand(sp.cancel(v * den)) for v in tw]
m = d["map_Legendre_to_factor_5"]["map"]
m = ast.literal_eval(m) if isinstance(m, str) else m
N0, N1, E = poly(m["N0"]), poly(m["N1"]), sp.expand(sp.sympify(m["E"]))
U, V = _residual(target, src, N0, N1, E)
ok = U == 0 and V == 0
cs = sp.Poly(N0, x).all_coeffs()[::-1]
k = max(range(len(cs)), key=lambda i: abs(cs[i]))
Up, Vp = _residual(target, src, N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
brk = not (Up == 0 and Vp == 0)
C2b = json.loads((ROOT / "tests" / "data" / "paper59_c2b_2026_10_22.json").read_text())["C2b_theta"]
same = all(sp.cancel(poly(a) / poly(b)).is_number for a, b in zip(d["factor_5"], C2b)) and \
    len({sp.cancel(poly(a) / poly(b)) for a, b in zip(d["factor_5"], C2b)}) == 1
print("factor_5(A y) = 0 for y in Sol(twisted L_Leg):", ok)
print("1/1000 perturbation breaks it:", brk)
print("factor_5 is C2b itself up to scale:", same, " (x-degrees", sp.degree(c2, x), "vs", sp.degree(poly(C2b[2]), x), ")")
json.dump({"exact": bool(ok), "perturbed_breaks": bool(brk), "factor_is_c2b": bool(same)},
          open(ROOT / "debug" / "data" / "ga1_pm_check.json", "w"), indent=1)
