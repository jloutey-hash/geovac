"""Minimality audit, step 2b: numerical test of the ring's elliptic content.
O' = (theta-2) o M has the hyperexponential solution t^2/(t^2-4) with M y = t^2, so
M h = -3 t^2 is solved by h_p = 3 t^2/(4 - t^2) (= r).  G carries the elliptic (Legendre) piece iff
h = (theta - r) G differs from h_p.  Compute h on the momentum form at several t."""
import json, sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
from minimality_audit_ring import numeric_h

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
ts = ["0.5", "1.3", "2.5", "3.5"]
res = numeric_h(ts, dps)
out = {}
for tv, d in res.items():
    t = mp.mpf(tv); hp = 3 * t**2 / (4 - t**2)
    alg = (1 - t**2 / 4) ** mp.mpf(-1.5)
    out[tv] = {k: mp.nstr(v, dps - 2) for k, v in d.items()}
    out[tv]["h_minus_hp"] = mp.nstr(d["h"] - hp, 8)
    out[tv]["(G+1)/alg"] = mp.nstr((d["G"] + 1) / alg, dps - 2)
    print(tv, "G =", mp.nstr(d["G"], 15), " h =", mp.nstr(d["h"], 15), " h - h_p =", mp.nstr(d["h"] - hp, 8),
          " (G+1)/alg =", mp.nstr((d["G"] + 1) / alg, 15), flush=True)
Path(__file__).parent.joinpath("data", f"minimality_audit_ring_num_dps{dps}.json").write_text(json.dumps(out, indent=1))
