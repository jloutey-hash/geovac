"""γ_n vs sqrt(n) fit for each of the three Dirichlet series."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from geovac.hp_operator import load_rhm_zeros


def main():
    data_path = ROOT / "debug" / "data" / "spectral_zero_stats.json"
    out = {}
    for which in ("D_full", "D_even", "D_odd"):
        ims = np.sort(load_rhm_zeros(str(data_path), which=which))
        n_arr = np.arange(1, len(ims) + 1, dtype=float)
        # fit ims = a * sqrt(n) + b
        sqrt_n = np.sqrt(n_arr)
        slope, intercept = np.polyfit(sqrt_n, ims, 1)
        fit = slope * sqrt_n + intercept
        rms = float(np.sqrt(np.mean((ims - fit) ** 2)))
        max_dev = float(np.max(np.abs(ims - fit)))
        # try also ims = a * sqrt(n) (no intercept)
        # least squares for a: min sum (ims - a*sqrt_n)^2 → a = (sqrt_n @ ims) / (sqrt_n @ sqrt_n)
        a_only = float(sqrt_n @ ims / (sqrt_n @ sqrt_n))
        fit_noint = a_only * sqrt_n
        rms_noint = float(np.sqrt(np.mean((ims - fit_noint) ** 2)))
        # try ims = a*n + b for comparison
        slope_n, intercept_n = np.polyfit(n_arr, ims, 1)
        rms_n = float(np.sqrt(np.mean((ims - slope_n * n_arr - intercept_n) ** 2)))
        # try ims = a*(n*log(n)) + b
        nlogn = np.where(n_arr >= 2, n_arr * np.log(np.maximum(n_arr, 2)), n_arr * 0.01)
        slope_nl, intercept_nl = np.polyfit(nlogn, ims, 1)
        rms_nl = float(np.sqrt(np.mean((ims - slope_nl * nlogn - intercept_nl) ** 2)))

        print(f"\n=== {which}: {len(ims)} zeros ===")
        print(f"  sqrt(n) fit: g ~ {slope:.4f} * sqrt(n) + {intercept:.4f}   RMS = {rms:.3f}   max_dev = {max_dev:.3f}")
        print(f"  sqrt(n) no-intercept: g ~ {a_only:.4f} * sqrt(n)   RMS = {rms_noint:.3f}")
        print(f"  linear-n fit:   g ~ {slope_n:.4f} * n + {intercept_n:.4f}   RMS = {rms_n:.3f}")
        print(f"  n*log(n) fit:   g ~ {slope_nl:.4f} * (n log n) + {intercept_nl:.4f}   RMS = {rms_nl:.3f}")
        print(f"  best model: {min([('sqrt_n', rms), ('sqrt_n_noint', rms_noint), ('linear', rms_n), ('nlogn', rms_nl)], key=lambda x: x[1])[0]}")

        out[which] = {
            "n_zeros": len(ims),
            "sqrt_n_fit": {
                "slope_a": float(slope),
                "intercept_b": float(intercept),
                "rms": rms,
                "max_dev": max_dev,
            },
            "sqrt_n_no_intercept": {"slope_a": a_only, "rms": rms_noint},
            "linear_fit": {"slope": float(slope_n), "intercept": float(intercept_n), "rms": rms_n},
            "nlogn_fit": {"slope": float(slope_nl), "intercept": float(intercept_nl), "rms": rms_nl},
        }

    # save
    out_path = ROOT / "debug" / "data" / "hp_operator_sqrt_n_fit.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
