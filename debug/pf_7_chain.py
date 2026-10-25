"""PF1-PF4 stage 7 (copy of debug/ga3_7_chain.py with pf_lib): an exact composition series of L by the generic recursive peeler; every inclusion checked by exact pseudo-division over Q.
Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_7_chain.py <op json> <out tag> [dmax]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import pf_lib as A

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
G.save_op(C.DATA / f"pf_pieces_{tag}.json", ops,
          {"note": "PF off-slice exp geminal (a=1, b=2, gamma_1=3/2, gamma_2=3). Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient "
                   "lists as flint polynomial strings. factor_i: composition factors, i = 1 the rightmost; chain_Ri = factor_i o ... o factor_1 "
                   "(exact right divisors of L).",
           "factors": [{"name": f["name"], "order": A.order(f["op"]), "x_degree": G.op_degree(f["op"]), "how": str(f["how"])} for f in factors],
           "divisions_exact": divs, "top_equals_L": ok_top, "seconds": time.time() - T0})
print("saved", flush=True)
