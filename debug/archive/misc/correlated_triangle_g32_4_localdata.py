"""GT stage 4: local data of the exact pieces (F = the order-3 left quotient, K6, and later the others)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc


def report(name, L):
    out = {"order": len(L) - 1, "x_degree": G.op_degree(L)}
    out["t0"] = [str(e) for e in Lc.exps_t0(L)[0]]
    out["inf_growth_mu"] = [str(e) for e in Lc.exps_inf_growth(L)[0]]
    for x0 in Lc.XPTS:
        out[f"x={x0}"] = [str(e) for e in Lc.exps_finite_t(L, x0)[0]]
    out["leading_factors"] = [(str(f), int(e)) for f, e in L[-1].factor()[1]]
    out["monodromy"] = {}
    for x0 in [Fraction(0)] + Lc.XPTS + ["inf"]:
        out["monodromy"][str(x0)] = Lc.monodromy_type(L, x0)
    print(f"== {name}: order {out['order']}, x-degree {out['x_degree']}")
    print("   t=0 exponents:", out["t0"], "  infinity growth mu:", out["inf_growth_mu"])
    for x0 in Lc.XPTS:
        print(f"   x={x0}:", out[f'x={x0}'])
    print("   leading coefficient factors:", [(f[:60], e) for f, e in out["leading_factors"]])
    for k, v in out["monodromy"].items():
        print(f"   local solutions at {k} (x-local exponents, multiplicity, log-free):", [c["exponents"] for c in v])
    return out


if __name__ == "__main__":
    d = json.loads((G.DATA / "correlated_triangle_g32_K6_F.json").read_text())
    F = G.parse_op(d["F"])
    K6 = G.parse_op(d["K6"])
    rec = {"F": report("F (order-3 left quotient)", F), "K6": report("K6", K6)}
    (G.DATA / "correlated_triangle_g32_4_localdata.json").write_text(json.dumps(rec, indent=1))
