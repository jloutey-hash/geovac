"""Probe A -- legs (b) and (c), split out of probeA_robustness.py.

(c) QSP head-room: does the aware degree depend on whether the polynomial is allowed to
    reach the QSP ceiling exactly (|p| <= 1) or is held to |p| <= 0.9?  Phase-factor
    synthesis is ill-conditioned when 1 - |p|^2 is tiny, so a usable degree must survive
    the stricter bound.

(b) Subnormalisation of the SHIFT arm.  Shifting to [-1,1] needs a block-encoding of
    B = (2S - (lam_max+lam_min) I)/(lam_max - lam_min); the natural LCU of S and I has
    1-norm ~3 (2*lam_max + lam_max + lam_min over lam_max - lam_min), so without
    re-amplification the spectrum lands in [-1/3, 1/3] -- interior, not at the endpoint.
    Arm "squeeze=3" repeats the shift arm with the spectrum compressed by 1/3 and asks
    whether the sqrt(kappa) advantage survives.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from geovac.sturmian_sigma_law import sw_cross_block, sigma_spectrum   # noqa: E402
import probeA_specaware_lp as core                                     # noqa: E402
from probeA_specaware_scan import bisect_degree, HEADROOM, N_GRID      # noqa: E402

core.BND_PER_DEG = 6
OUT = os.path.join(HERE, "data", "probeA_legs_bc.json")


def spectrum(s, n):
    sig = sigma_spectrum(sw_cross_block(s, n))
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def bisect_bound(x, y, eps, seed, dmax, bound):
    old = core.min_rel_error
    core.min_rel_error = lambda a, b, d, bound=bound: old(a, b, d, bound=bound)
    try:
        return bisect_degree(x, y, eps, seed, dmax)
    finally:
        core.min_rel_error = old


out = []
print("(c) QSP head-room -- aware/p60 minimal degree at |p| <= 1.0 vs |p| <= 0.9")
for s in (2.0,):
    for n in (6, 10, 14, 18):
        lt = spectrum(s, n)
        kap = float(lt[-1] / lt[0])
        x, y, _ = core.make_problem(lt, "aware", "p60", headroom=HEADROOM)
        res = {}
        for b in (1.0, 0.9):
            d, r, coef, _ = bisect_bound(x, y, 1e-3, max(4 * n, int(0.30 * kap)), 900, b)
            res[b] = int(d)
            out.append(dict(leg="qspbound", s=s, n=n, kappa=kap, bound=b, d=int(d)))
            json.dump(out, open(OUT, "w"))
        print(f"    s={s} n={n:2d} kappa={kap:7.1f}  d(|p|<=1.0)={res[1.0]:4d}  "
              f"d(|p|<=0.9)={res[0.9]:4d}  ratio={res[0.9]/max(res[1.0],1):.2f}", flush=True)

print()
print("(b) shift arm under LCU subnormalisation (spectrum squeezed into [-1/q, 1/q])")
for s in (2.0,):
    for n in (4, 6, 8, 10):
        lam = spectrum(s, n)
        lo, hi = float(lam[0]), float(lam[-1])
        kap = hi / lo
        for q in (1.0, 3.0):
            for mode in ("aware", "generic"):
                if mode == "aware":
                    lm = lam
                else:
                    lm = np.unique(np.concatenate([
                        lo * (hi / lo) ** np.linspace(0, 1, N_GRID),
                        0.5 * (lo + hi) - 0.5 * (hi - lo) *
                        np.cos(np.pi * np.arange(N_GRID) / (N_GRID - 1)), lam]))
                x = (2 * lm - (lo + hi)) / (hi - lo) / q
                y = HEADROOM * np.sqrt(lo / lm)
                seed = max(8, int((2.0 if q == 1.0 else 0.35 * np.sqrt(kap)) * np.sqrt(kap)))
                d, r, coef, _ = bisect_degree(x, y, 1e-3, seed, 900)
                out.append(dict(leg="squeeze", s=s, n=n, kappa=kap, q=q, mode=mode,
                                d=int(d)))
                json.dump(out, open(OUT, "w"))
                print(f"    s={s} n={n:2d} kappa={kap:7.2f} squeeze={q} {mode:7s} d={d}",
                      flush=True)
