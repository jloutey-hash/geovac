"""Grid ladder for Paper 60's 1-norm families at a chosen radial box.

Same rungs and same norm families as debug/sturmian_exact_ladder.py, but built
with the production quadrature engine geovac/sturmian_secular.py at a settable
(R_MAX, N_GRID).  Used to expose the box-truncation bias as a function of K:
the production box R_MAX = 60 is smaller than the 10s-14s orbital extent at
Q ~ 1, so the bias is systematic IN K -- exactly the fitted axis.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

import geovac.sturmian_secular as ss

Z = 2.0


def regrid(rmax: float, ngrid: int) -> None:
    ss.R_MAX = rmax
    ss.N_GRID = ngrid
    ss.r = np.linspace(1e-7, rmax, ngrid)
    ss.dr = ss.r[1] - ss.r[0]
    ss.r2 = ss.r * ss.r
    ss.reset_caches()


def rung(n_max: int, rmax: float, ngrid: int, lmax: int = 3):
    regrid(rmax, ngrid)
    cfgs = ss.build_configs(ss.gen_configs(lmax, {l: n_max for l in range(lmax + 1)}))
    M = ss.build_M(cfgs, Z)
    Rnu = np.array([c.Rnu for c in cfgs])
    t = M - np.diag(Z * Rnu)
    d = np.abs(np.diag(M)).sum()
    tot = np.abs(M).sum()
    return dict(K=len(cfgs), M1=float(tot), M1_diag=float(d),
                M1_off=float(tot - d), T0=float(Z * Rnu.sum()),
                Tp1=float(np.abs(t).sum()),
                Tp1_diag=float(np.abs(np.diag(t)).sum()),
                E0=float(-np.sort(np.linalg.eigvalsh(M))[-1] ** 2 / 2))


if __name__ == "__main__":
    rmax = float(sys.argv[1])
    ngrid = int(sys.argv[2])
    top = int(sys.argv[3])
    out = sys.argv[4] if len(sys.argv) > 4 else None
    rows = []
    for n in range(4, top + 1):
        t0 = time.time()
        r = rung(n, rmax, ngrid)
        rows.append(r)
        print("n=%2d K=%4d  ||M||1=%.5f diag=%.5f off=%.5f T0=%.5f "
              "||T'||1=%.5f  E0=%.6f  (%.0fs)"
              % (n, r["K"], r["M1"], r["M1_diag"], r["M1_off"], r["T0"],
                 r["Tp1"], r["E0"], time.time() - t0), flush=True)
    if out:
        json.dump(rows, open(out, "w"), indent=1)
        print("wrote", out)
