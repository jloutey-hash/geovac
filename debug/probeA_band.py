"""Probe A -- leg (a2), done properly: the DEPLOYABLE spectrum-aware design.

The pure "aware" arm designs the polynomial on the exact eigenvalues.  The deployable
version would design it on the closed-form sigma law with a safety band around each
predicted eigenvalue,

    accuracy required on   [lam_k^law / (1+w),  lam_k^law * (1+w)]   for every k,

and never sees the true spectrum.  Question: does d_band stay near d_aware (exploit
survives on the law) or degrade to d_generic (exploit needs the exact spectrum)?

Each band is sampled uniformly in the RESCALED variable at spacing <= 1/(4d), because a
degree-d polynomial oscillates on scale ~1/d and a fixed handful of samples per band
silently under-resolves it (that error is what a first pass hit).  Every reported degree
is re-verified on a 4x finer band grid AND on the true spectrum.
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
from probeA_specaware_scan import HEADROOM                             # noqa: E402

core.BND_PER_DEG = 6
OUT = os.path.join(HERE, "data", "probeA_band.json")


def spectrum(s, n):
    sig = sigma_spectrum(sw_cross_block(s, n))
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def law_spectrum(s, n):
    k = np.arange(1, n + 1)
    gap = np.clip((s ** 2 / 24.0) * (k * np.pi) ** 2 / n ** 2, 1e-12, 1.0)
    sig = 1.0 - gap
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def band_grid(lam_law, w, lo, hi, d, dens=4.0):
    """Union of bands, each sampled at spacing <= 1/(dens*d) in x = lam/hi."""
    out = []
    for lk in lam_law:
        a, b = max(lk / (1.0 + w), 1e-12), min(lk * (1.0 + w), hi)
        m = int(np.ceil(dens * d * (b - a) / hi)) + 3
        out.append(np.linspace(a, b, min(m, 400)))
    return np.unique(np.concatenate(out))


def r_star_band(lam_law, w, lo, hi, d):
    lm = band_grid(lam_law, w, lo, hi, d, 4.0)
    x, y = lm / hi, HEADROOM * np.sqrt(lo / lm)
    r, coef = core.min_rel_error(x, y, d)
    fine = band_grid(lam_law, w, lo, hi, d, 16.0)
    xf, yf = fine / hi, HEADROOM * np.sqrt(lo / fine)
    pf = core.cheb_design(np.clip(xf, -1, 1), d) @ coef
    return float(np.max(np.abs(pf - yf) / yf)), coef, len(lm)


def min_degree_band(lam_law, w, lo, hi, eps, d_seed, d_max=700):
    d = max(4, d_seed)
    r, coef, _ = r_star_band(lam_law, w, lo, hi, d)
    if r <= eps:
        lo_d, hi_d, best = 1, d, coef
    else:
        lo_d, hi_d, best = d, -1, coef
        while d < d_max:
            d = min(int(np.ceil(1.7 * d)), d_max)
            r, coef, _ = r_star_band(lam_law, w, lo, hi, d)
            if r <= eps:
                hi_d, best = d, coef
                break
            lo_d = d
        if hi_d < 0:
            return -1, coef
    while hi_d - lo_d > 1:
        mid = (hi_d + lo_d) // 2
        r, c, _ = r_star_band(lam_law, w, lo, hi, mid)
        if r <= eps:
            hi_d, best = mid, c
        else:
            lo_d = mid
    return hi_d, best


if __name__ == "__main__":
    rows = []
    for s, n in ((1.4, 6), (2.0, 8), (3.0, 12), (2.0, 12)):
        lt, ll = spectrum(s, n), law_spectrum(s, n)
        lo, hi = float(lt[0]), float(lt[-1])
        kap = hi / lo
        miss = float(np.max([np.min(np.abs(np.log(tv / ll))) for tv in lt]))
        print(f"--- s={s} n={n} kappa={kap:.1f}  law worst per-eigenvalue miss = "
              f"{100*(np.exp(miss)-1):.1f}%", flush=True)
        for w in (0.05, 0.15, 0.30):
            d, coef = min_degree_band(ll, w, lo, hi, 1e-3, max(4 * n, int(0.25 * kap)))
            if d < 0:
                print(f"    w={w:.2f}: > d_max", flush=True)
                continue
            xt = lt / hi
            yt = HEADROOM * np.sqrt(lo / lt)
            pt = core.cheb_design(np.clip(xt, -1, 1), d) @ coef
            err = float(np.max(np.abs(pt - yt) / yt))
            rows.append(dict(s=s, n=n, kappa=kap, w=w, d=int(d), err_on_true=err,
                             law_max_miss=miss, covered=bool(np.exp(miss) - 1.0 <= w)))
            print(f"    w={w:.2f}  d={d:4d}  err_on_TRUE={err:.2e}  "
                  f"bands_cover_truth={np.exp(miss)-1.0 <= w}", flush=True)
            json.dump(rows, open(OUT, "w"))
    print("ALLDONE")
