"""Unscored sensitivity for DS1: P1 by the ODE route with hbar c = 197.33 instead of the registered 197.327
(potential strengths and kinetic coefficient both rescaled; ranges fixed from the thesis lambdas)."""
from __future__ import annotations

import json
import math

import deuteron_second_ode as O

out = {}
for hbc in (197.327, 197.33):
    pot = [("y", hbc * 7.40758, 1 / 0.321543), ("y", -hbc * 3.22749, 1 / 0.643087)]
    hb = hbc * hbc / 938.92
    runs = [O.bound_state(pot, hb, h) for h in O.HS]
    B = O.richardson([x["B"] for x in runs][-3:])[0]
    r2 = O.richardson([x["r2"] for x in runs][-3:])[0]
    at = O.richardson([O.scattering_length(pot, hb, h)[0] for h in O.HS][-3:])[0]
    out[str(hbc)] = {"B": B, "r_m": 0.5 * math.sqrt(r2), "a_t": at,
                     "rel_B_vs_2.27203": (B - 2.27203) / 2.27203, "rel_a_t_vs_5.4739": (at - 5.4739) / 5.4739,
                     "rel_r_m_vs_1.97625": (0.5 * math.sqrt(r2) - 1.97625) / 1.97625}
    print(hbc, out[str(hbc)], flush=True)
json.dump(out, open("data/deuteron_second_hbarc.json", "w"), indent=1)
