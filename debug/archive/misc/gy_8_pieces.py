"""GY stage 8 (copy of ga3_8_pieces.py with gy_lib): local data of every composition factor (exponents at t = 0, at the
candidate points, growth / Newton polygon at infinity), apparent-point checks on the leading coefficient's factors of
degree >= 2 and on its rational roots, irreducibility of the order-2 and order-3 factors by the exponent-bounded split
searches WHERE THEY ARE VALID (Fuchsian at infinity), with the chain's planted splits as positive controls.
usage: python gy_8_pieces.py <pieces json> <out tag> [tpos comma list]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C
import gy_lib as A

d = json.loads(Path(sys.argv[1]).read_text())
tag = sys.argv[2]
if len(sys.argv) > 3:
    pts = sys.argv[3]
    if pts.startswith("x:"):
        A.set_xpoints([v for v in pts[2:].split(",")])
    else:
        A.set_points([v for v in pts.split(",")])
ops = {k: G.parse_op(v) for k, v in d.items() if isinstance(v, list) and v and isinstance(v[0], str)}
XP = A.XPTS
rec = {}
T0 = time.time()
factors = sorted([k for k in ops if k.startswith("factor_")], key=lambda k: int(k.split("_")[1]))
for name in factors + ["L"]:
    L = ops[name]
    rec[name] = A.local_report(L, name=name)
    if len(L) - 1 == 1:
        rec[name]["solution_exponents"] = {"t=0 (t-power)": rec[name]["t0"],
                                            **{f"x={x0} ((x-x0)-power)": rec[name]["points"][str(x0)]["exps"] for x0 in XP},
                                            "inf (t-growth, if Fuchsian)": rec[name]["inf_mu"]}
    rec[name]["rational_roots_of_lead"] = [(str(r_), e) for r_, e in C.rational_roots(L)]
    for f, e in C.lead_factors(L):
        if f.degree() >= 2:
            rec[name].setdefault("apparent", {})[f"degree {f.degree()} (mult {e})"] = C.apparent_check(L, f, nprimes=2)
    print(f"({time.time()-T0:.0f}s)", flush=True)
rec["split"] = {}
for name in factors:
    L = ops[name]
    if len(L) - 1 < 2:
        continue
    if not A.fuchsian_inf(L):
        rec["split"][name] = {"skipped": "not Fuchsian at infinity; the exponent-bounded split searches do not apply"}
        print(f"== split search {name}: skipped (irregular at infinity)", flush=True)
        continue
    print(f"== split search {name} (order {len(L)-1})", flush=True)
    r = {"sol": S.split_search(L, "sol")}
    if len(L) - 1 == 3:
        r["wronskian"] = S.split_search(L, "wronskian")
    rec["split"][name] = r
chains = sorted([k for k in ops if k.startswith("chain_R")], key=lambda k: int(k[7:]))
for name in chains:
    L = ops[name]
    if len(L) - 1 in (2, 3) and A.fuchsian_inf(L):
        print(f"== planted-split control {name} (order {len(L)-1})", flush=True)
        r = {"sol": S.split_search(L, "sol")}
        if len(L) - 1 == 3:
            r["wronskian"] = S.split_search(L, "wronskian")
        rec["split"][name + " (planted control)"] = r
print(f"({time.time()-T0:.0f}s)", flush=True)
(C.DATA / f"gy_8_pieces_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))
