"""CF5-CF7 stage 8 (gamma = 2): local data of every composition factor; irreducibility of the order-2 and
order-3 factors (exponent-bounded split searches of record g32, with the planted split V4 = B2q o rb o ra
and B2 as positive controls); the order-one factors' solutions (exponents at every point).
usage: python cf_g2_8_pieces.py"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C

C.set_points([Fraction(1), Fraction(2), Fraction(3), Fraction(4)])
d = json.loads((C.DATA / "cf_g2_pieces.json").read_text())
ops = {k: G.parse_op(v) for k, v in d.items() if isinstance(v, list)}
XP = [Fraction(1), Fraction(4), Fraction(9), Fraction(16)]
rec = {}
T0 = time.time()
for name in ("factor_ra", "factor_rb", "factor_B2q", "factor_rho1", "factor_C2b", "factor_F"):
    L = ops[name]
    rec[name] = C.local_report(L, XP, name=name)
    if len(L) - 1 == 1:
        # order-1 factor c1 theta + c0: solution exponents = local exponents (t at 0, (x - x0) at x0, t-growth at inf)
        rec[name]["solution_exponents"] = {"t=0 (t-power)": rec[name]["t0"],
                                            **{f"x={x0} ((x-x0)-power)": rec[name]["points"][str(x0)]["exps"] for x0 in XP},
                                            "inf (t-growth)": rec[name]["inf_mu"]}
    print(f"({time.time()-T0:.0f}s)", flush=True)
rec["split"] = {}
for name in ("chain_V4", "B2", "factor_C2b", "factor_F"):
    L = ops[name]
    print(f"== split search {name} (order {len(L)-1})", flush=True)
    r = {"sol": S.split_search(L, "sol")}
    if len(L) - 1 == 3:
        r["wronskian"] = S.split_search(L, "wronskian")
    rec["split"][name] = r
print(f"({time.time()-T0:.0f}s)", flush=True)
(C.DATA / "cf_g2_8_pieces.json").write_text(json.dumps(rec, indent=1, default=str))
