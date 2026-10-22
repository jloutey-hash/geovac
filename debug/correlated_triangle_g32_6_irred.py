"""GT stage 6: split searches (order-1 right factors; for order 3 also order-2 right factors) on the
order-3 pieces F and C2, with the planted split V4 = B2q o rb o ra (order 3, right factors of order 1
and 2 by construction) and B2 (order 2, contains R2's algebraic solution) as positive controls."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S

d = json.loads((G.DATA / "correlated_triangle_g32_pieces.json").read_text())
ops = {k: G.parse_op(d[k]) for k in ("V4", "B2", "C2", "F")}
rec = {}
for name in ("V4", "B2", "C2", "F"):
    L = ops[name]
    t0 = time.time()
    print(f"== {name} (order {len(L)-1})", flush=True)
    r = {}
    r["sol"] = S.split_search(L, "sol")
    if len(L) - 1 == 3:
        r["wronskian"] = S.split_search(L, "wronskian")
    print(f"   ({time.time()-t0:.0f}s)", flush=True)
    rec[name] = r
(G.DATA / "correlated_triangle_g32_6_irred.json").write_text(json.dumps(rec, indent=1, default=str))
