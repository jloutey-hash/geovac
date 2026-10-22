"""GS step A2: L (step A, order 17) applied to G_b by the cached grid (generic_slice_num)."""
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_num as N  # noqa: E402
from generic_slice_B import coeff_list, Q  # noqa: E402

dps, hinv = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
g = N.Grid(2, dps, hinv)
print(f"grid {len(g.nodes)} nodes, f built in {g.build_time:.0f}s", flush=True)
L = coeff_list(Q)
rec = {}
for ts in ["0.5", "1.3", "2.5", "3.5"]:
    t1 = time.time()
    G = N.apply(g, [sp.Integer(1)], ts)
    v, big, err = N.apply(g, L, ts)
    rec[ts] = {"G": mp.nstr(G[0], 30), "G_quad_err": mp.nstr(G[2], 5), "LG": mp.nstr(v, 8),
               "largest": mp.nstr(big, 8), "quad_err": mp.nstr(err, 5)}
    print(ts, rec[ts], f"({time.time()-t1:.0f}s)", flush=True)
Path(__file__).with_name("data").joinpath(f"generic_slice_A2_{dps}_{hinv}.json").write_text(json.dumps(rec, indent=1))
