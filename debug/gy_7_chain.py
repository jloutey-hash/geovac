"""GY stage 7 (copy of ga3_7_chain.py with gy_lib): an exact composition series of the operator L (Q-form json from
gy_5e, or a theta-form key of a pieces json) by the generic recursive peeler; every inclusion checked by exact
pseudo-division over Q.  usage: python gy_7_chain.py <op json> <key|Q> <out tag> [dmax] [tpos comma list]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import cf_g2_lib as C
import gy_lib as A

T0 = time.time()
d = json.loads(Path(sys.argv[1]).read_text())
key = sys.argv[2]
L = C.load_Q(sys.argv[1]) if key == "Q" else G.parse_op(d[key])
tag = sys.argv[3]
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 80
if len(sys.argv) > 5:
    pts = sys.argv[5]
    if pts.startswith("x:"):
        A.set_xpoints([v for v in pts[2:].split(",")])
    else:
        A.set_points([v for v in pts.split(",")])
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
G.save_op(C.DATA / f"gy_pieces_{tag}.json", ops,
          {"note": "gamma = 1, Gaussian geminal g = exp(-gamma r^2). Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient "
                   "lists as flint polynomial strings. factor_i: composition factors, i = 1 the rightmost; chain_Ri = factor_i o ... o factor_1 "
                   "(exact right divisors of L).",
           "factors": [{"name": f["name"], "order": A.order(f["op"]), "x_degree": G.op_degree(f["op"]), "how": str(f["how"]),
                        "fuchsian_at_infinity": A.fuchsian_inf(f["op"])} for f in factors],
           "divisions_exact": divs, "top_equals_L": ok_top, "points_used": [str(v) for v in A.TPOS], "seconds": time.time() - T0})
print("saved", flush=True)
