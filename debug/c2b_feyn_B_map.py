"""CF2 (docs/predictions/2026-10-22_c2b_feynman.md, committed d68e6fe5): an exact order-<=1 map from the
gauged Legendre operator L_Leg (degree -1/2, z = (t^2 + 1 - g^2)/(2t), twist t^(-1/2), g = 3/2) to C2b.

L_Leg is even in t: c2 = t e2(x), c1 = e1(x), c0 = t e0(x), so in x = t^2 it is
4x e2 H_xx + (2 e2 + 2 e1) H_x + e0 H.  The search is over Q(x) (stronger than the registered Q(t)),
with the EC Stage E linear algebra: A = (N0 + N1 d/dx)/E, E = ((4x-1)(4x-25))^k (4x+5)^m x^n.
Control: the same call recovers the EC map from the EC source (c2b_ident_E_map source_op).
Usage: python debug/c2b_feyn_B_map.py"""
import json
import sys
from pathlib import Path

import flint
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import c2b_ident_ops as O  # noqa: E402
from c2b_ident_E_map import fq, int_nullspace, apply_c2b_to_image, gauge  # noqa: E402
from c2b_feyn_A import leg_op, t  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def leg_x():
    c0, c1, c2 = leg_op(sp.Rational(3, 2))
    e2 = sp.expand(sp.cancel(c2 / t).subs(t, sp.sqrt(x)))
    e1 = sp.expand(c1.subs(t, sp.sqrt(x)))
    e0 = sp.expand(sp.cancel(c0 / t).subs(t, sp.sqrt(x)))
    b = [e0, 2 * e2 + 2 * e1, 4 * x * e2]
    g = sp.gcd_list([sp.expand(v) for v in b])
    return [sp.expand(sp.cancel(v / g)) for v in b]


def search(src, Eexpr, dN):
    b = [fq(v) for v in src]
    ac, _ = O.c2b()
    a = [flint.fmpq_poly([flint.fmpq(int(c)) for c in aj]) for aj in ac]
    E = fq(Eexpr)
    X = flint.fmpq_poly([0, 1])
    zero = flint.fmpq_poly([])
    cols, labels = [], []
    for j in range(dN + 1):
        cols.append(apply_c2b_to_image(a, b, X**j, zero, E)); labels.append(("N0", j))
    for j in range(dN + 2):
        cols.append(apply_c2b_to_image(a, b, zero, X**j, E)); labels.append(("N1", j))
    deg = max(max(c[0].degree(), c[1].degree()) for c in cols) + 1
    rows = []
    for r in range(deg):
        rows.append([c[0].coeffs()[r] if r < len(c[0].coeffs()) else 0 for c in cols])
        rows.append([c[1].coeffs()[r] if r < len(c[1].coeffs()) else 0 for c in cols])
    NS, nul = int_nullspace(rows, len(cols))
    out = {"E": str(sp.factor(Eexpr)), "dN": dN, "nullity": nul}
    if nul >= 1:
        v = [int(NS[i, 0]) for i in range(len(cols))]
        N0 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N0")
        N1 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N1")
        U, V = apply_c2b_to_image(a, b, fq(N0), fq(N1), E)
        out.update({"N0": str(N0), "N1": str(N1), "exact": U.is_zero() and V.is_zero(),
                    "order_one": N1 != 0})
    return out


if __name__ == "__main__":
    src0 = leg_x()
    tw = gauge(src0, sp.Rational(1, 2) / x)
    gg = sp.gcd_list(tw)
    src = [sp.expand(sp.cancel(v / gg)) for v in tw] if "notwist" not in sys.argv else src0
    print("L_Leg in x:", [sp.factor(v) for v in src])
    res = {"source": [str(v) for v in src], "runs": []}
    q = (4 * x - 1) * (4 * x - 25)
    for k, m, dN in [(0, 1, 8), (1, 1, 12), (2, 1, 16), (3, 1, 20), (5, 1, 22), (5, 2, 24)]:
        r = search(src, q**k * (4 * x + 5) ** m, dN)
        print(f"k = {k}, m = {m}, dN = {dN}: nullity {r['nullity']}", {kk: r[kk] for kk in ('exact', 'order_one') if kk in r}, flush=True)
        res["runs"].append(r)
        if r["nullity"] >= 1:
            break
    json.dump(res, open(DATA / "c2b_feyn_B_map.json", "w"), indent=1)
