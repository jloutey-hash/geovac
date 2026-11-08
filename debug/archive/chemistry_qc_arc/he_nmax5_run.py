"""Run the exact-ingredient certification at n_max = 5 and dump the result."""
from __future__ import annotations

import json
import time
from pathlib import Path

import mpmath as mp

from debug.he_graph_native_certify import (build_matrix_exact, ground_state_mp,
                                           realise_mp)

OUT = Path("debug/data/he_gnci_nmax5.json")
OUT.parent.mkdir(parents=True, exist_ok=True)

t0 = time.time()
entries, cfg = build_matrix_exact(2, 5)
t_asm = time.time() - t0
print(f"dim={len(cfg)}  exact assembly {t_asm:.1f}s", flush=True)

res = {}
for dps in (55, 75):
    t1 = time.time()
    with mp.workdps(dps):
        H = realise_mp(entries, len(cfg))
        lam, resid, lam_f = ground_state_mp(H)
        res[dps] = dict(lam=mp.nstr(lam, 50, strip_zeros=False),
                        resid=mp.nstr(resid, 4),
                        float64=repr(lam_f),
                        seconds=round(time.time() - t1, 1))
    print(f"  dps={dps}  {res[dps]}", flush=True)

with mp.workdps(90):
    a = mp.mpf(res[55]["lam"])
    b = mp.mpf(res[75]["lam"])
    agree = int(mp.floor(-mp.log10(abs(a - b) / abs(b)))) if a != b else 10 ** 6
res["agree_digits"] = agree
res["dim"] = len(cfg)
res["assembly_seconds"] = round(t_asm, 1)
OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
print("agree_digits =", agree, "-> wrote", OUT, flush=True)
