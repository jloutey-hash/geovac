"""GY stage 8b: deeper right/left-factor probes of the composition factors (gy_lib.probe with a larger degree bound):
for every log-free Frobenius exponent of F and of F* at t = 0 and the candidate points, the minimal annihilator
(order < ord F, x-degree <= dmax) of that local solution mod p.  A hit of order r < ord F is a proper factor the
chain missed; no hit to dmax is the stated bound.  usage: python gy_8b_probe.py <pieces json> <key> <dmax> <out tag> [x:pts|tpts]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import gy_lib as A

d = json.loads(Path(sys.argv[1]).read_text())
key, dmax, tag = sys.argv[2], int(sys.argv[3]), sys.argv[4]
if len(sys.argv) > 5:
    pts = sys.argv[5]
    A.set_xpoints(pts[2:].split(",")) if pts.startswith("x:") else A.set_points(pts.split(","))
F = G.parse_op(d[key])
t0 = time.time()
print(f"== probe {key}: order {len(F)-1}, x-degree {G.op_degree(F)}, dmax {dmax}, points {A.XPTS}", flush=True)
res = A.probe(F, dmax, rmax=len(F) - 2, verbose=True)
hits = [r for r in res if r["r"] < len(F) - 1]
print(f"proper-factor hits: {hits}  ({time.time()-t0:.0f}s)", flush=True)
(G.DATA / f"gy_8b_probe_{tag}.json").write_text(json.dumps({"key": key, "order": len(F) - 1, "dmax": dmax, "points": [str(x) for x in A.XPTS],
                                                              "hits": [{k: str(v) for k, v in h.items()} for h in hits], "seconds": time.time() - t0}, indent=1))
