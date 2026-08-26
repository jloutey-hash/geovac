"""Probe A -- optional leg: water (C2v) A1 canonical-correlation block.

Same LP machinery on a spectrum that is NOT of the two-center 1 +/- sigma form:
the C2v-adapted A1 block of water's three-center SW metric (the block that holds
the ground state and the symmetry-irremovable O<->H coupling; Paper 60 sec:molecular).
Tests whether the aware-vs-generic verdict is a property of the 1 +/- sigma
clustering or of finite-dimensionality as such.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from numpy.linalg import eigvalsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import probeA_specaware_lp as core                                # noqa: E402
from probeA_specaware_scan import bisect_degree, HEADROOM, N_GRID, d_model  # noqa: E402
import sturmian_sw_water_conditioning as W                        # noqa: E402

core.BND_PER_DEG = 6
OUT = os.path.join(HERE, "data", "probeA_water.json")

rows = []
for nmax in (3, 4, 6, 8, 10, 12):
    _, a1, _, _ = W.water_blocks(nmax, "S")
    lam = np.sort(eigvalsh(a1))
    kap = float(lam[-1] / lam[0])
    for mode, seed in (("aware", int(max(2 * (2 * nmax) + 6, 0.30 * kap))),
                       ("generic", int(2.7 * kap))):
        x, y, c = core.make_problem(lam, mode, "p60", n_grid=N_GRID, headroom=HEADROOM)
        d, r, coef, _ = bisect_degree(x, y, 1e-3, seed, 600)
        rows.append(dict(nmax=nmax, dim=len(lam), kappa=kap, mode=mode, d=int(d),
                         r=float(r), d_model=float(d_model(kap, 1e-3))))
        print(f"water nmax={nmax:2d} dim={len(lam):2d} kappa={kap:8.2f} {mode:7s} "
              f"d={d:4d} r={r:.2e}  d_model={d_model(kap,1e-3):7.0f}", flush=True)
        json.dump(rows, open(OUT, "w"))
