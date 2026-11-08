"""GA1 stage 7: an exact composition series of L (the minimal annihilator of Y_lin,log at gamma = 3/2), by the
generic recursive peeler of ga1_lib (record-g32 stages 3, 5, 7 made generic: the gamma = 2 chain script was
written for L9's particular exponents).  Every inclusion is checked by exact pseudo-division over Q.
usage: python ga1_7_chain.py <op json> <out tag> [dmax]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import ga1_lib as A

T0 = time.time()
L = C.load_Q(sys.argv[1])
tag = sys.argv[2]
dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 80
factors, chain = A.composition_series(L, dmax, name="L")
print(f"composition factors (right to left): {[(f['name'], A.order(f['op'])) for f in factors]}  ({time.time()-T0:.0f}s)", flush=True)
divs = {}
for i in range(len(chain)):
    ok, _ = A.divides(L, chain[i])
    divs[f"R{i+1} (order {A.order(chain[i])}) divides L"] = ok
    if i:
        ok2, _ = A.divides(chain[i], chain[i - 1])
        divs[f"R{i+1} = . o R{i}"] = ok2
    print(f"   R{i+1}: order {A.order(chain[i])}, x-degree {G.op_degree(chain[i])}; divides L: {ok}", flush=True)
ok_top = A.order(chain[-1]) == A.order(L) and A.divides(chain[-1], L)[0] and A.divides(L, chain[-1])[0]
print("top of the chain equals L (mutual exact division):", ok_top, flush=True)
ops = {f"factor_{i+1}": f["op"] for i, f in enumerate(factors)}
ops.update({f"chain_R{i+1}": c for i, c in enumerate(chain)})
ops["L"] = L
G.save_op(C.DATA / f"ga1_pieces_{tag}.json", ops,
          {"note": "gamma = 3/2, f = r exp(-gamma r). Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient "
                   "lists as flint polynomial strings. factor_i: composition factors, i = 1 the rightmost; chain_Ri = factor_i o ... o factor_1 "
                   "(exact right divisors of L).",
           "factors": [{"name": f["name"], "order": A.order(f["op"]), "x_degree": G.op_degree(f["op"]), "how": str(f["how"])} for f in factors],
           "divisions_exact": divs, "top_equals_L": ok_top, "seconds": time.time() - T0})
print("saved", flush=True)
