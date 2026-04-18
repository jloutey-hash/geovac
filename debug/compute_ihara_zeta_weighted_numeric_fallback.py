"""
Numeric-only fallback: compute Ramanujan verdicts and Hashimoto
diagnostics at chosen alpha values for all 4 Dirac-S^3 weighted
configurations, WITHOUT the expensive full symbolic determinant.

Writes to a separate JSON; the full symbolic driver
(compute_ihara_zeta_weighted.py) may still be running in parallel
on the biggest case (Rule B n_max=3).

Run:
  python debug/compute_ihara_zeta_weighted_numeric_fallback.py
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import sympy as sp

from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.ihara_zeta import _count_components, ihara_zeta_bass
from geovac.ihara_zeta_weighted import (
    alpha_sym,
    build_weighted_dirac_adjacency,
    weighted_hashimoto_matrix,
    is_weighted_ramanujan,
    hashimoto_pair_correlation_cv,
)

OUT_JSON = os.path.join(
    os.path.dirname(__file__),
    "data",
    "ihara_zeta_weighted_numeric_fallback.json",
)

ALPHA_PHYS = 1.0 / 137.036
ALPHA_UNIT = 1.0


def main():
    results = {}
    for (n_max, rule) in [(2, "A"), (2, "B"), (3, "A"), (3, "B")]:
        key = f"rule{rule}_nmax{n_max}"
        t0 = time.perf_counter()
        print(f"[numeric] rule={rule} n_max={n_max}", flush=True)
        A01, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
        A_w, _ = build_weighted_dirac_adjacency(n_max, rule)
        V = int(A01.shape[0])
        E = int(A01.sum()) // 2
        c = _count_components(A01)
        r = E - V + c

        # alpha sweep for Ramanujan deviation
        alpha_sweep = []
        for alpha_val in [0.0, 1e-6, 1e-4, 1e-2, 1 / 137.036, 0.05, 0.1, 0.3, 0.5, 1.0]:
            A_num = A_w.subs(alpha_sym, alpha_val)
            is_ram, dev, text = is_weighted_ramanujan(A_num)
            cv = hashimoto_pair_correlation_cv(A_num)
            alpha_sweep.append({
                "alpha": alpha_val,
                "deviation": float(dev),
                "is_ramanujan": bool(is_ram),
                "cv_spacing": cv["cv_spacing"],
                "n_spacings": cv["n_spacings"],
            })

        elapsed = time.perf_counter() - t0
        results[key] = {
            "n_max": n_max,
            "adjacency_rule": rule,
            "V": V,
            "E": E,
            "c": c,
            "r_betti1": r,
            "max_degree": int(deg.max()),
            "alpha_sweep": alpha_sweep,
            "elapsed_seconds": elapsed,
        }
        print(f"  done ({elapsed:.2f}s)", flush=True)

    out = {
        "meta": {
            "module": "geovac.ihara_zeta_weighted",
            "track": "RH-K (numeric fallback)",
            "alpha_physical": ALPHA_PHYS,
            "alpha_unitary": ALPHA_UNIT,
            "weight_convention": (
                "w(a,b) = 1 + alpha^2 * |f_SO(a) + f_SO(b)| / 2, "
                "f_SO(n,kappa) = -(kappa+1) / [4 n^3 l (l+1/2)(l+1)]; "
                "f_SO = 0 for l=0 (Kramers)."
            ),
        },
        "results": results,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
