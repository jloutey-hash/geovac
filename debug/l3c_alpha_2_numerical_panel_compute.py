"""Sprint L3c-alpha.2 — Numerical verification panel for parametric de-compactification.

Computes Lambda^L(n_max, N_t, T(N_t)) at coupled cells along three admissible
scalings:

  T_0(N_t) = 2 pi                  (canonical BW period; not admissible — fixed T;
                                    baseline only)
  T_1(N_t) = N_t / log(N_t)        (sub-linear admissible: T_1 -> inf, T_1/N_t -> 0)
  T_2(N_t) = sqrt(N_t)             (square-root admissible: T_2 -> inf, T_2/N_t -> 0)

Sweep: n_max in {2, 3, 4} x N_t in {3, 5, 7, 11, 21}.

Verifies the L3c-alpha analytical theorem (Paper 47 §3, Theorem inner-arrow)
empirically: along admissible scalings, the propinquity bound vanishes.

Output: debug/data/l3c_alpha_2_numerical_panel.json
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

from geovac.lorentzian_propinquity_compact_temporal import (
    compute_lorentzian_propinquity_bound,
)


# Panel configuration.
# Reduced from initial {2, 3, 4} x {3, 5, 7, 11, 21} after first run hung at
# (n_max=4, N_t=21) cell — Krein space dim is ~1680, propinquity compute scales
# super-linearly. Smaller panel suffices to verify the rate trend.
N_MAX_VALUES = [2, 3]
N_T_VALUES = [3, 5, 7, 11]
GAMMA_PREC = 30


def scaling_T0(N_t: int) -> float:
    """Canonical BW radius (fixed T = 2 pi). NOT admissible (T does not -> inf)."""
    return 2.0 * math.pi


def scaling_T1(N_t: int) -> float:
    """Sub-linear admissible scaling T = N_t / log(N_t)."""
    return N_t / math.log(N_t)


def scaling_T2(N_t: int) -> float:
    """Square-root admissible scaling T = sqrt(N_t)."""
    return math.sqrt(N_t)


SCALINGS = {
    "T0_canonical_2pi": scaling_T0,
    "T1_sublinear_Nt_over_logNt": scaling_T1,
    "T2_sqrt_Nt": scaling_T2,
}


def main() -> None:
    out_dir = Path("debug/data")
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_cells: list[dict] = []
    t_start = time.time()

    for scaling_name, scaling_fn in SCALINGS.items():
        for n_max in N_MAX_VALUES:
            for N_t in N_T_VALUES:
                T = scaling_fn(N_t)
                t_cell_start = time.time()
                try:
                    bound = compute_lorentzian_propinquity_bound(
                        n_max=n_max,
                        N_t=N_t,
                        T=T,
                        gamma_prec=GAMMA_PREC,
                    )
                    cell = bound.to_dict()
                    cell["scaling"] = scaling_name
                    cell["t_over_Nt"] = T / N_t
                    cell["compute_time_sec"] = time.time() - t_cell_start
                    cell["status"] = "ok"
                    print(
                        f"[{scaling_name}] n_max={n_max} N_t={N_t} T={T:.4f} "
                        f"T/N_t={T/N_t:.4f} Lambda={cell['propinquity_bound']:.6f} "
                        f"({cell['compute_time_sec']:.2f}s)"
                    )
                except Exception as exc:
                    cell = {
                        "scaling": scaling_name,
                        "n_max": n_max,
                        "N_t": N_t,
                        "T": T,
                        "t_over_Nt": T / N_t,
                        "status": "error",
                        "error": repr(exc),
                        "compute_time_sec": time.time() - t_cell_start,
                    }
                    print(
                        f"[{scaling_name}] n_max={n_max} N_t={N_t} T={T:.4f} ERROR: "
                        f"{exc!r}"
                    )
                panel_cells.append(cell)

    total_time = time.time() - t_start
    print(f"\nTotal compute time: {total_time:.1f}s")

    # Summary: verify admissible scalings have Lambda -> 0 as N_t -> infinity.
    summary = {}
    for scaling_name in SCALINGS:
        cells = [c for c in panel_cells if c.get("scaling") == scaling_name and c.get("status") == "ok"]
        per_n_max = {}
        for n_max in N_MAX_VALUES:
            n_max_cells = sorted([c for c in cells if c["n_max"] == n_max], key=lambda c: c["N_t"])
            if len(n_max_cells) >= 2:
                # Monotone-decrease check along N_t for fixed n_max.
                lambdas = [c["propinquity_bound"] for c in n_max_cells]
                monotone = all(lambdas[i] >= lambdas[i + 1] for i in range(len(lambdas) - 1))
                # Ratio Lambda(largest N_t) / Lambda(smallest N_t).
                ratio_largest_smallest = lambdas[-1] / lambdas[0] if lambdas[0] != 0 else None
                per_n_max[f"n_max={n_max}"] = {
                    "N_t_sequence": [c["N_t"] for c in n_max_cells],
                    "Lambda_sequence": lambdas,
                    "monotone_decreasing": monotone,
                    "ratio_largest_smallest_Nt": ratio_largest_smallest,
                }
        summary[scaling_name] = per_n_max

    payload = {
        "sprint": "L3c-alpha.2 numerical verification panel",
        "date": "2026-05-23",
        "n_max_values": N_MAX_VALUES,
        "N_t_values": N_T_VALUES,
        "scalings": list(SCALINGS),
        "panel_cells": panel_cells,
        "summary_by_scaling": summary,
        "total_compute_time_sec": total_time,
        "gamma_prec": GAMMA_PREC,
    }

    out_path = out_dir / "l3c_alpha_2_numerical_panel.json"
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # Print human-readable summary.
    print("\n" + "=" * 70)
    print("L3c-alpha.2 Verdict Summary")
    print("=" * 70)
    for scaling_name, per_n_max in summary.items():
        print(f"\nScaling {scaling_name}:")
        for n_max_key, data in per_n_max.items():
            mono = "MONOTONE DECREASING" if data["monotone_decreasing"] else "NOT monotone"
            ratio = data["ratio_largest_smallest_Nt"]
            ratio_str = f"{ratio:.4f}" if ratio is not None else "n/a"
            print(f"  {n_max_key}: {mono} (ratio Lambda(N_t=21)/Lambda(N_t=3) = {ratio_str})")


if __name__ == "__main__":
    main()
