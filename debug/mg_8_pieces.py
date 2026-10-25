"""MG1-MG3 stage 8 (copy of cf_g2_8_pieces.py, parametrised): local data of every composition factor at every
candidate point x = t^2 (the four pair masses of the given gamma plus every rational root of L's leading
coefficient); irreducibility of the order-2 and order-3 factors (exponent-bounded split searches of record g32,
with the planted splits V = B2q o rb o ra and B2 as positive controls).
usage: python mg_8_pieces.py <pieces json> <gamma> <out tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C

PF, GAMMA, TAG = sys.argv[1], Fraction(sys.argv[2]), sys.argv[3]
TPOS = sorted({abs(1 - GAMMA), Fraction(2), 1 + GAMMA, 2 * GAMMA})
C.set_points(TPOS)
d = json.loads(Path(PF).read_text())
ops = {k: G.parse_op(v) for k, v in d.items() if isinstance(v, list)}
XP = sorted(set([v * v for v in TPOS] + [a for a, e in C.rational_roots(ops["L"]) if a != 0]))
rec = {"gamma": str(GAMMA), "x_points": [str(v) for v in XP]}
T0 = time.time()
for name in [k for k in ops if k.startswith("factor_")]:
    L = ops[name]
    rec[name] = C.local_report(L, XP, name=name)
    if len(L) - 1 == 1:
        rec[name]["solution_exponents"] = {"t=0 (t-power)": rec[name]["t0"],
                                            **{f"x={x0} ((x-x0)-power)": rec[name]["points"][str(x0)]["exps"] for x0 in XP},
                                            "inf (t-growth)": rec[name]["inf_mu"]}
    print(f"({time.time()-T0:.0f}s)", flush=True)
rec["split"] = {}
for name in [k for k in ("chain_V", "B2", "factor_C2b", "factor_C2", "factor_F") if k in ops]:
    L = ops[name]
    print(f"== split search {name} (order {len(L)-1})", flush=True)
    r = {"sol": S.split_search(L, "sol")}
    if len(L) - 1 == 3:
        r["wronskian"] = S.split_search(L, "wronskian")
    rec["split"][name] = r
print(f"({time.time()-T0:.0f}s)", flush=True)
(C.DATA / f"mg_8_pieces_{TAG}.json").write_text(json.dumps(rec, indent=1, default=str))
