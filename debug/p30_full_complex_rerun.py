"""Paper 30 Sec. 5.2 reruns on the full plaquette complexes (group5 debt, 2026-11-08).

The printed n_max = 3 monopole fit (c = 8.95) used 200 of the 994 plaquettes
and 80 of the 3342 three-plaquette 2-cycles; the n_max = 4 Monte Carlo Wilson
loops (sigma = -0.080 at beta = 0.3) used 1500 of the 5047 plaquettes.  This
driver reruns both with the caps lifted and every other setting unchanged.

Usage:
    python debug/p30_full_complex_rerun.py mono3    # monopole density, n_max = 3
    python debug/p30_full_complex_rerun.py wilson4  # MC Wilson loops, n_max = 4
Output: debug/data/p30_full_complex_<mode>.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUPPORT = ROOT / "tests" / "wilson_rule_b_support"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SUPPORT))

OUT_DIR = ROOT / "debug" / "data"
FULL = 10 ** 6  # no cap


def mono3() -> dict:
    from xcwg_monopole_density import run_at_nmax
    # Settings as the printed n_max = 3 run, caps lifted.
    return run_at_nmax(
        n_max=3,
        beta_grid=[0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 30.0],
        n_therm=800,
        n_sample=400,
        sample_interval=30,
        cycle_max_count=FULL,
        cycle_time_budget_sec=3600.0,
        rng_seed=43,
        plaq_cap=FULL,
        larger_cycle_cap=3,  # elementary sites only (k_delta = 3)
    )


def wilson4() -> dict:
    from xcwg_full_mc_wilson_loops import run_at_nmax
    # Settings as the printed n_max = 4 run, plaquette cap lifted.
    return run_at_nmax(
        n_max=4,
        beta_grid=[0.3, 1.0, 3.0],
        n_therm=600,
        n_sample=3000,
        sample_interval=10,
        n_loops_per_area=8,
        rng_seed=44,
        A_targets=[1, 2, 3, 4, 5, 6],
        max_clusters_per_area=30,
        plaq_cap=FULL,
    )


def main() -> None:
    mode = sys.argv[1]
    t0 = time.time()
    res = {"mono3": mono3, "wilson4": wilson4}[mode]()
    res["wall_seconds"] = time.time() - t0
    out = OUT_DIR / f"p30_full_complex_{mode}.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
