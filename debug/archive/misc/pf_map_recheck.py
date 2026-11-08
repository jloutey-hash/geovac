"""PF1: re-check of a Legendre map on a code path independent of pf_map.py / cf_g2_map.search: generalises
debug/pm_map_check.py (the PM's recheck, B = 1 only) to the two-scale source L(t, B; i gamma), whose Legendre argument
z = (t^2 + B^2 - gamma^2)/(2 B t) is the B = 1 argument at gamma' = gamma/B in t' = t/B.  So the source is built from
c2b_feyn_A.leg_op(gamma/B) (the slice code, not pf_map.leg_op_general), moved to x' = x/B^2 by the chain rule, twisted by
hand, and the stored map is applied with the sympy residual of tests/test_paper59_c2b_elliptic.py (no flint); a 1/1000
change in the largest N0 coefficient must break it.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_map_recheck.py <B> <gamma> <pieces json> <factor key> <map json> <map key or -> <out tag>"""
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


def source_x(B, g, twist):
    """the twisted Legendre operator in d/dx for H(t) = t'^(-1/2) w(z), t' = t/B, from leg_op(g/B) in t'"""
    l0, l1, l2 = leg_op(g / B)                       # c2 H'' + c1 H' + c0 H = 0 in d/dt'
    tp = sp.Symbol("tp", positive=True)
    l0, l1, l2 = (sp.expand(c.subs(t, tp)) for c in (l0, l1, l2))
    # d/dt' = B d/dt: c2 B^2 H_tt + c1 B H_t + c0 H, with t' = t/B; then x = t^2 (d/dt = 2 t d/dx, d^2/dt^2 = 4x d^2/dx^2 + 2 d/dx)
    c2 = sp.expand(l2.subs(tp, t / B) * B**2)
    c1 = sp.expand(l1.subs(tp, t / B) * B)
    c0 = sp.expand(l0.subs(tp, t / B))
    # operator on H in t: c2 H'' + c1 H' + c0 H; in x: c2 (4x H_xx + 2 H_x) + c1 (2 t H_x) + c0 H
    b2 = 4 * x * c2
    b1 = 2 * c2 + 2 * t * c1
    b0 = c0
    # clear odd powers of t: the operator is even or odd in t overall; divide by t if needed
    ops = [b0, b1, b2]
    tst = [sp.expand(o.subs(t, -t)) for o in ops]
    if all(sp.expand(o - q) == 0 for o, q in zip(ops, tst)):
        pass
    elif all(sp.expand(o + q) == 0 for o, q in zip(ops, tst)):
        ops = [sp.expand(sp.cancel(o / t)) for o in ops]
    else:
        raise RuntimeError("source not of definite parity in t")
    ops = [sp.expand(o.subs(t, sp.sqrt(x))) for o in ops]
    for o in ops:
        assert not o.has(sp.sqrt(x))
    # H = t'^(-1/2) w = (x/B^2)^(-1/4) w; leg_op's H is the same function of t' (its gauge), so the x-operator above acts
    # on H; the stored maps are from the operator twisted by a further x^(1/2) (H = x^(-1/2) y) or not twisted
    b0, b1, b2 = ops
    if twist == "x^(1/2)":
        L = -sp.Rational(1, 2) / x
        tw = [b0 + L * b1 + (L**2 + sp.diff(L, x)) * b2, b1 + 2 * L * b2, b2]
    else:
        tw = [b0, b1, b2]
    den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in tw])
    return [sp.expand(sp.cancel(v * den)) for v in tw]


def main(B, g, pieces, fkey, mapfile, mkey, tag):
    B, g = sp.Rational(B), sp.Rational(g)
    d = json.loads(Path(pieces).read_text())
    c0, c1, c2 = [poly(v) for v in d[fkey]]
    target = [c0, 2 * x * c1 + 4 * x * c2, 4 * x**2 * c2]          # theta = 2x d/dx
    md = json.loads(Path(mapfile).read_text())
    m = md if mkey == "-" else md[mkey]
    m = ast.literal_eval(m) if isinstance(m, str) else m
    twist = m.get("twist", "x^(1/2)")
    src = source_x(B, g, twist)
    N0, N1, E = poly(m["N0"]), poly(m["N1"]), sp.expand(sp.sympify(str(m["E"])))
    U, V = _residual(target, src, N0, N1, E)
    ok = U == 0 and V == 0
    cs = sp.Poly(N0, x).all_coeffs()[::-1]
    k = max(range(len(cs)), key=lambda i: abs(cs[i]))
    Up, Vp = _residual(target, src, N0 + sp.Rational(1, 1000) * cs[k] * x**k, N1, E)
    brk = not (Up == 0 and Vp == 0)
    print(f"B = {B}, gamma = {g}, twist {twist}, E = {sp.factor(E)}: {fkey}(A y) = 0 for y in Sol(source): {ok}; 1/1000 breaks: {brk}")
    out = {"B": str(B), "gamma": str(g), "z": f"(t^2 + {B**2 - g**2})/({2 * B} t)", "pieces": str(pieces), "factor": fkey,
           "map_file": str(mapfile), "twist": twist, "E": str(E), "exact": bool(ok), "perturbed_breaks": bool(brk)}
    json.dump(out, open(ROOT / "debug" / "data" / f"pf_map_recheck_{tag}.json", "w"), indent=1)
    return ok and brk


if __name__ == "__main__":
    a = sys.argv[1:]
    sys.exit(0 if main(*a) else 1)
