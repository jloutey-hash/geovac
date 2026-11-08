"""PM recheck of CF5 (docs/predictions/2026-10-22_c2b_feynman_g2.md) on a code path independent of the
sub-agent's: L_Leg at gamma = 2 rebuilt from c2b_feyn_A.leg_op (Legendre substitution in t), moved to x by
hand and twisted by x^(1/2); the theta-form C2b factor converted to d/dx by hand; the stored map applied
with the sympy residual of tests/test_paper59_c2b_elliptic.py (no flint).
Also: the local exponents of the stored C2b at x = 0 and the L9g exponent match are not rechecked here.
Usage: python debug/cf_g2_pm_check.py"""
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
d = json.loads((ROOT / "debug" / "data" / "cf_g2_pieces.json").read_text())


def poly(s_):
    return sp.expand(sp.sympify(s_.replace("^", "**"), locals={"x": x}))


c0, c1, c2 = [poly(v) for v in d["factor_C2b"]]
target = [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]

l0, l1, l2 = leg_op(sp.Integer(2))                 # in d/dt, for H = t^(-1/2) w(z)
e2 = sp.expand(sp.cancel(l2 / t).subs(t, sp.sqrt(x)))
e1 = sp.expand(l1.subs(t, sp.sqrt(x)))
e0 = sp.expand(sp.cancel(l0 / t).subs(t, sp.sqrt(x)))
src = [e0, 2 * e2 + 2 * e1, 4 * x * e2]          # H in d/dx
# twist y = x^(1/2) H: H = x^(-1/2) y
L = sp.Rational(-1, 2) / x
b0, b1, b2 = src
# H = phi y, phi'/phi = L: H' = phi (y' + L y), H'' = phi (y'' + 2 L y' + (L^2 + L') y)
tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
tw = [sp.expand(sp.cancel(v * den)) for v in tw]

m = d["map_CF5_Leg_to_C2b"]
N0, N1, E = poly(m["N0"]), poly(m["N1"]), sp.expand(sp.sympify(m["E"]))
U, V = _residual(target, tw, N0, N1, E)
print("C2b(A y) = 0 for y in Sol(twisted L_Leg):", U == 0 and V == 0)
cs = sp.Poly(N0, x).all_coeffs()[::-1]
k = max(range(len(cs)), key=lambda i: abs(cs[i]))
Up, Vp = _residual(target, tw, N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
print("1/1000 perturbation breaks it:", not (Up == 0 and Vp == 0))
z2 = (x - 3) ** 2 / (4 * x)
j = sp.factor(sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2))
print("j(gamma = 2):", j)
json.dump({"exact": bool(U == 0 and V == 0), "perturbed_breaks": bool(not (Up == 0 and Vp == 0)), "j": str(j)},
          open(ROOT / "debug" / "data" / "cf_g2_pm_check.json", "w"), indent=1)
