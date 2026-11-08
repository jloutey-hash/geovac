"""Unscored diagnostic: the ~1.5e-8 basis-vs-ODE r_m floor for P2 at N_r = 35, 40.
(a) Laguerre at N_r = 45, 50 (lambda_40 and a lambda search at 50); (b) ODE with r_max = 80 fm.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction as Fr

import mpmath as mp

import deuteron_second_laguerre as L
import deuteron_second_ode as O

out = {}
pot, hb = L.POTS["P2"]
for N, lam in ((45, Fr(936, 625)), (50, Fr(936, 625))):
    r = L.solve(N, lam, pot, hb)
    out[f"lag_N{N}_lam{lam}"] = {"B": mp.nstr(-r["E"], 18), "r_m": mp.nstr(mp.sqrt(r["r2"]) / 2, 18)}
    print(N, lam, out[f"lag_N{N}_lam{lam}"], flush=True)
lam50, E50, _, edge = L.lam_opt(50, pot, hb)
r = L.solve(50, lam50, pot, hb)
out["lag_N50_opt"] = {"lambda": str(lam50), "edge": edge, "B": mp.nstr(-r["E"], 18), "r_m": mp.nstr(mp.sqrt(r["r2"]) / 2, 18)}
print(out["lag_N50_opt"], flush=True)
fp, fhb = O.fpot("P2")
O.RMAX = 80.0
runs = [O.bound_state(fp, fhb, h) for h in O.HS]
r2 = O.richardson([x["r2"] for x in runs][-3:])
B = O.richardson([x["B"] for x in runs][-3:])
out["ode_rmax80"] = {"B": B[0], "r2_raw": [x["r2"] for x in runs], "r_m": 0.5 * math.sqrt(r2[0]),
                     "r2_last_change": r2[1]}
print(out["ode_rmax80"], flush=True)
json.dump(out, open("data/deuteron_second_floor.json", "w"), indent=1)
