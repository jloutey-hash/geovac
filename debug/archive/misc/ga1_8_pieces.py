"""GA1 stage 8: local data of every composition factor; irreducibility of the order-2 and order-3 factors
(exponent-bounded split searches of record g32: 'sol' = order-one right factors, 'wronskian' = order-two
right factors of an order-three piece), with planted splits (the chain's reducible partial compositions)
as positive controls; the order-one factors' solutions (exponents at every point).  Copy of
debug/cf_g2_8_pieces.py made generic over the factor list of a ga1 pieces json.
usage: python ga1_8_pieces.py <pieces json> <out tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C
import ga1_lib as A

d = json.loads(Path(sys.argv[1]).read_text())
tag = sys.argv[2]
ops = {k: G.parse_op(v) for k, v in d.items() if isinstance(v, list) and v and isinstance(v[0], str)}
XP = A.XPTS
rec = {}
T0 = time.time()
factors = sorted([k for k in ops if k.startswith("factor_")], key=lambda k: int(k.split("_")[1]))
for name in factors:
    L = ops[name]
    rec[name] = C.local_report(L, XP, name=name)
    if len(L) - 1 == 1:
        rec[name]["solution_exponents"] = {"t=0 (t-power)": rec[name]["t0"],
                                            **{f"x={x0} ((x-x0)-power)": rec[name]["points"][str(x0)]["exps"] for x0 in XP},
                                            "inf (t-growth)": rec[name]["inf_mu"]}
    for f, e in C.lead_factors(L):
        if f.degree() >= 2:
            rec[name].setdefault("apparent", {})[f"degree {f.degree()} (mult {e})"] = C.apparent_check(L, f, nprimes=2)
    print(f"({time.time()-T0:.0f}s)", flush=True)
rec["split"] = {}
for name in factors:
    L = ops[name]
    if len(L) - 1 < 2:
        continue
    print(f"== split search {name} (order {len(L)-1})", flush=True)
    r = {"sol": S.split_search(L, "sol")}
    if len(L) - 1 == 3:
        r["wronskian"] = S.split_search(L, "wronskian")
    rec["split"][name] = r
# planted-split controls: the chain's partial compositions of order 2 and 3 are reducible by construction
chains = sorted([k for k in ops if k.startswith("chain_R")], key=lambda k: int(k[7:]))
for name in chains:
    L = ops[name]
    if len(L) - 1 in (2, 3):
        print(f"== planted-split control {name} (order {len(L)-1})", flush=True)
        r = {"sol": S.split_search(L, "sol")}
        if len(L) - 1 == 3:
            r["wronskian"] = S.split_search(L, "wronskian")
        rec["split"][name + " (planted control)"] = r
print(f"({time.time()-T0:.0f}s)", flush=True)
(C.DATA / f"ga1_8_pieces_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))
