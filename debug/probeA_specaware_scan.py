"""Probe A -- main scan: spectrum-AWARE vs GENERIC minimal QSVT degree for S^{-1/2}.

Four arms (2 accuracy sets x 2 rescaling conventions), one LP machine
(debug/probeA_specaware_lp.py), on the H2+-style two-center Shibuya-Wulfman metric
spectrum {1 +/- sigma_k} built by geovac.sturmian_sigma_law.

Writes debug/data/probeA_scan.json incrementally.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from geovac.sturmian_sigma_law import sw_cross_block, sigma_spectrum  # noqa: E402
import probeA_specaware_lp as core                                     # noqa: E402

core.BND_PER_DEG = 6          # Ehlich-Zeller guarantee sup|p| <= 1/cos(pi/12) = 1.035
N_GRID = 130
HEADROOM = 0.5                # c = 0.5 * sqrt(lam_min); held fixed across all arms
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "probeA_scan.json")


def spectrum(s: float, n: int) -> np.ndarray:
    sig = sigma_spectrum(sw_cross_block(s, n))
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def bisect_degree(x, y, eps, d_seed, d_max):
    """Smallest d with r*(d) <= eps.  Returns (d, r_at_d, coef, n_lp)."""
    n_lp = 0
    d = max(2, d_seed)
    r, coef = core.min_rel_error(x, y, d); n_lp += 1
    if r <= eps:
        hi, best = d, (r, coef)
        lo = 1
        while hi - lo > 1:                      # search downward
            mid = (hi + lo) // 2
            rm, cm = core.min_rel_error(x, y, mid); n_lp += 1
            if rm <= eps:
                hi, best = mid, (rm, cm)
            else:
                lo = mid
        return hi, best[0], best[1], n_lp
    lo, hi, best = d, -1, (r, coef)
    while d < d_max:
        d = min(int(np.ceil(1.6 * d)), d_max)
        r, coef = core.min_rel_error(x, y, d); n_lp += 1
        if r <= eps:
            hi, best = d, (r, coef)
            break
        lo = d
    if hi < 0:
        return -1, best[0], best[1], n_lp
    while hi - lo > 1:
        mid = (hi + lo) // 2
        rm, cm = core.min_rel_error(x, y, mid); n_lp += 1
        if rm <= eps:
            hi, best = mid, (rm, cm)
        else:
            lo = mid
    return hi, best[0], best[1], n_lp


def d_model(kappa: float, eps: float) -> float:
    return kappa * np.log(kappa / eps)


# (mode, convention) -> (d_max, seed_fn).  Seeds are calibrated slightly ABOVE the
# expected answer from the pilot (d ~ 2.3 kappa generic/p60, ~2.8 sqrt(kappa) shift),
# so the bisection walks DOWNWARD through progressively cheaper LPs.
ARMS = {
    ("aware", "p60"):     (900,  lambda k, nn, e: max(2 * nn + 6, int(0.30 * k))),
    ("generic", "p60"):   (600,  lambda k, nn, e: int((2.7 if e >= 1e-3 else 5.0) * k)),
    ("aware", "shift"):   (400,  lambda k, nn, e: max(6, int(1.6 * np.sqrt(k)))),
    ("generic", "shift"): (400,  lambda k, nn, e: max(8, int((3.2 if e >= 1e-3 else 5.5)
                                                             * np.sqrt(k)))),
}

# per-arm nmax windows (the generic/p60 arm has degree ~2*kappa and is the cost driver)
WINDOW = {
    ("aware", "p60"):     {1.4: (4, 18), 2.0: (4, 18), 3.0: (4, 18)},
    ("generic", "p60"):   {1.4: (4, 6), 2.0: (4, 8), 3.0: (4, 12)},
    ("aware", "shift"):   {1.4: (4, 18), 2.0: (4, 18), 3.0: (4, 18)},
    ("generic", "shift"): {1.4: (4, 18), 2.0: (4, 18), 3.0: (4, 18)},
}

#: the tight-eps leg is run only where the LP stays affordable (kappa cap per arm)
KAPPA_CAP_EPS5 = {("generic", "p60"): 45.0}

NS = [4, 6, 8, 10, 12, 14, 16, 18]
SS = [2.0, 3.0, 1.4]
EPS = [1e-3, 1e-5]


def main():
    rows = []
    if os.path.exists(OUT):
        rows = json.load(open(OUT))
    done = {(r["s"], r["n"], r["mode"], r["conv"], r["eps"]) for r in rows}
    for s in SS:
        for n in NS:
            lam = spectrum(s, n)
            kap = float(lam[-1] / lam[0])
            for (mode, conv), (dmax, seed) in ARMS.items():
                lo_n, hi_n = WINDOW[(mode, conv)][s]
                if not (lo_n <= n <= hi_n):
                    continue
                x, y, c = core.make_problem(lam, mode, conv, n_grid=N_GRID,
                                            headroom=HEADROOM)
                for eps in EPS:
                    key = (s, n, mode, conv, eps)
                    if key in done:
                        continue
                    if eps < 1e-4 and kap > KAPPA_CAP_EPS5.get((mode, conv), 1e18):
                        continue
                    t0 = time.time()
                    d, r, coef, nlp = bisect_degree(x, y, eps, seed(kap, n, eps), dmax)
                    rec = dict(s=s, n=n, dim=2 * n, kappa=kap, mode=mode, conv=conv,
                               eps=eps, d=int(d), r_at_d=float(r),
                               sup=float(core.true_sup(coef)) if d > 0 else None,
                               d_model=float(d_model(kap, eps)), n_lp=nlp,
                               secs=round(time.time() - t0, 1),
                               coef=[float(v) for v in coef] if d > 0 and d <= 400 else None)
                    rows.append(rec)
                    json.dump(rows, open(OUT, "w"))
                    print(f"s={s} n={n:2d} kap={kap:8.2f} {mode:7s} {conv:5s} "
                          f"eps={eps:g}  d={d:4d}  r={r:.2e} sup={rec['sup']:.4f} "
                          f"({rec['secs']}s, {nlp} LPs)", flush=True)


if __name__ == "__main__":
    main()
