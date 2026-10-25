"""PM recheck of a Legendre map (MG1, MG2, GA3; docs/predictions/2026-10-23_*.md) on a code path independent of
the sub-agents': the Legendre operator at the given gamma rebuilt from c2b_feyn_A.leg_op, moved to x and twisted
by hand; the stored order-two factor (theta form) converted to d/dx by hand; the stored map applied with the
sympy residual of tests/test_paper59_c2b_elliptic.py (no flint).  Generalises debug/ga1_pm_check.py and
debug/cf_g2_pm_check.py to any gamma, any pieces file and either twist.
Usage: python debug/pm_map_check.py <gamma> <pieces json> <factor key> <map json> [<map key>] [<out tag>]
  <map json>: a ga3_map_*/mg_map_* file (map under key "map") or a pieces file with the map under <map key>."""
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


def poly(s_):
    return sp.expand(sp.sympify(str(s_).replace("^", "**"), locals={"x": x}))


def main(gamma, pieces, fkey, mapfile, mkey=None, tag=None):
    g = sp.Rational(gamma)
    d = json.loads(Path(pieces).read_text())
    c = [poly(v) for v in d[fkey]]
    assert len(c) == 3, "factor must have order two"
    c0, c1, c2 = c
    target = [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]           # theta = 2x d/dx
    md = json.loads(Path(mapfile).read_text())
    m = md[mkey] if mkey else md["map"]
    if isinstance(m, dict) and "map" in m and "N0" not in m:
        m = m["map"]
    m = ast.literal_eval(m) if isinstance(m, str) else m
    twist = m.get("twist", "x^(1/2)")
    l0, l1, l2 = leg_op(g)                                            # H = t^(-1/2) w(z), d/dt
    e2 = sp.expand(sp.cancel(l2 / t).subs(t, sp.sqrt(x)))
    e1 = sp.expand(l1.subs(t, sp.sqrt(x)))
    e0 = sp.expand(sp.cancel(l0 / t).subs(t, sp.sqrt(x)))
    for e in (e0, e1, e2):
        assert not e.has(sp.sqrt(x))
    b0, b1, b2 = e0, 2 * e2 + 2 * e1, 4 * x * e2                      # H in d/dx
    if twist == "x^(1/2)":
        L = -sp.Rational(1, 2) / x                                    # H = x^(-1/2) y
        tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
    else:
        tw = [b0, b1, b2]
    den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
    src = [sp.expand(sp.cancel(v * den)) for v in tw]
    N0, N1, E = poly(m["N0"]), poly(m["N1"]), sp.expand(sp.sympify(str(m["E"])))
    U, V = _residual(target, src, N0, N1, E)
    ok = U == 0 and V == 0
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    Up, Vp = _residual(target, src, N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
    brk = not (Up == 0 and Vp == 0)
    z = (x + 1 - g**2) ** 2 / (4 * x)
    j = sp.factor(sp.cancel(64 * (z + 3) ** 3 / (1 - z) ** 2))
    print(f"gamma = {g}, twist {twist}, E = {sp.factor(E)}, deg N0 = {sp.degree(N0, x)}, deg N1 = {sp.degree(N1, x)}")
    print(f"{fkey}(A y) = 0 for y in Sol(twisted L_Leg({g})):", ok)
    print("1/1000 perturbation breaks it:", brk)
    print("j(z):", j)
    out = {"gamma": str(g), "pieces": str(pieces), "factor": fkey, "map_file": str(mapfile), "twist": twist,
           "E": str(E), "exact": bool(ok), "perturbed_breaks": bool(brk), "j": str(j)}
    tag = tag or f"g{str(g).replace('/', 'over')}_{fkey}"
    json.dump(out, open(ROOT / "debug" / "data" / f"pm_map_check_{tag}.json", "w"), indent=1)
    return ok and brk


if __name__ == "__main__":
    a = sys.argv[1:]
    good = main(*a)
    sys.exit(0 if good else 1)
