"""Probe A -- leg (a2), cheap ladder version: the DEPLOYABLE spectrum-aware design.

Instead of bisecting for the exact minimal band degree (expensive: the bands near the top
of the spectrum merge and the LP grows), evaluate r*_band(d) on a ladder of degrees
anchored on the two known answers -- d_aware and d_generic -- and report where the band
design crosses eps.  Also report the resulting polynomial's error on the TRUE spectrum,
which is the only thing that matters operationally.

Bands: accuracy on [lam_k^law/(1+w), lam_k^law*(1+w)], sampled at spacing <= 1/(4d) in the
rescaled variable (a degree-d polynomial oscillates on scale ~1/d), and re-verified on a 4x
finer band grid.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from probeA_band import spectrum, law_spectrum, r_star_band                  # noqa: E402
import probeA_specaware_lp as core                                           # noqa: E402
from probeA_specaware_scan import HEADROOM                                   # noqa: E402

core.BND_PER_DEG = 6
OUT = os.path.join(HERE, "data", "probeA_band2.json")

# (s, n, d_aware, d_generic) from debug/data/probeA_scan.json, eps = 1e-3
CASES = [(1.4, 6, 24, 243), (2.0, 8, 27, 209), (3.0, 12, 33, 203)]

rows = []
for s, n, d_aw, d_gen in CASES:
    lt, ll = spectrum(s, n), law_spectrum(s, n)
    lo, hi = float(lt[0]), float(lt[-1])
    kap = hi / lo
    miss = float(np.max([np.min(np.abs(np.log(tv / ll))) for tv in lt]))
    print(f"--- s={s} n={n} kappa={kap:.1f}  d_aware={d_aw} d_generic={d_gen}  "
          f"law worst per-eigenvalue miss={100*(np.exp(miss)-1):.1f}%", flush=True)
    for w in (0.05, 0.30):
        for d in sorted({d_aw, 2 * d_aw, 4 * d_aw, d_gen}):
            r, coef, npts = r_star_band(ll, w, lo, hi, d)
            xt = lt / hi
            yt = HEADROOM * np.sqrt(lo / lt)
            pt = core.cheb_design(np.clip(xt, -1, 1), d) @ coef
            err = float(np.max(np.abs(pt - yt) / yt))
            rows.append(dict(s=s, n=n, kappa=kap, w=w, d=int(d), r_band=r,
                             err_on_true=err, npts=npts,
                             covered=bool(np.exp(miss) - 1.0 <= w)))
            print(f"    w={w:.2f} d={d:4d}  r*_band={r:.3e}  err_on_TRUE={err:.3e}  "
                  f"(acc pts {npts}, bands cover truth: "
                  f"{np.exp(miss)-1.0 <= w})", flush=True)
            json.dump(rows, open(OUT, "w"))
print("BAND2_DONE")
