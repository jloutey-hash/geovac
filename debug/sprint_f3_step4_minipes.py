"""Sprint F3 Step 4 — Mini-PES at NaH max_n=2 with cross-block h1.

5-7 R-points to confirm whether there's an internal minimum and what its
shape looks like. Tests the three architectures (bare, W1c+mz,
W1c+mz+xblockh1) on the same R grid for direct comparison.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, '.')

from debug.sprint_f3_step3_fci import run_arch
from geovac.molecular_spec import nah_spec


def main():
    out_path = Path('debug/data/sprint_f3_step4_minipes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    spec = nah_spec(max_n=2)
    R_grid = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]

    architectures = [
        ('bare',
            dict(screened_cross_center=False, multi_zeta_basis=False,
                 cross_block_h1=False)),
        ('W1c+mz',
            dict(screened_cross_center=True, multi_zeta_basis=True,
                 cross_block_h1=False)),
        ('W1c+mz+xblockh1',
            dict(screened_cross_center=True, multi_zeta_basis=True,
                 cross_block_h1=True)),
    ]

    print("=" * 78)
    print("Sprint F3 Step 4 — mini-PES at NaH max_n=2")
    print("=" * 78)
    print(f"{'arch':>20} {'R':>5}: {'E_gs (Ha)':>14} {'dom_NO':>8} {'Na/H':>11}")

    results = {}
    for label, kwargs in architectures:
        arch_results = []
        for R in R_grid:
            r = run_arch(spec, R, label, **kwargs)
            arch_results.append(r)
            print(f"{label:>20} {R:>5.1f}: {r['E_gs']:>+14.6f} "
                  f"{r['dominant_NO']['occupation']:>8.4f} "
                  f"{r['dominant_NO']['Na_amp_squared']:>5.2f}/"
                  f"{r['dominant_NO']['H_amp_squared']:.2f}")
        results[label] = arch_results

    # Identify R_min for each architecture
    summary = {}
    print(f"\n{'arch':>20} {'R_min':>7} {'E_min':>14} {'D_e (Ha)':>10} {'internal_min':>13}")
    for label in [a[0] for a in architectures]:
        E_arr = [r['E_gs'] for r in results[label]]
        R_arr = [r['R'] for r in results[label]]
        i_min = int(np.argmin(E_arr))
        E_min = E_arr[i_min]
        R_min = R_arr[i_min]
        # internal_min = True if argmin is not at the boundary
        internal_min = (i_min != 0 and i_min != len(R_arr) - 1)
        # D_e: well depth at R_min from dissociation
        E_diss = E_arr[-1]  # at R=10
        D_e = E_diss - E_min
        summary[label] = {
            'R_min': R_min,
            'E_min': E_min,
            'D_e_well_depth_Ha': D_e,
            'internal_min': internal_min,
            'E_at_R3.5': float(E_arr[R_arr.index(3.5)]),
            'E_at_R10.0': float(E_arr[R_arr.index(10.0)]),
        }
        print(f"{label:>20} {R_min:>7.1f} {E_min:>+14.6f} {D_e:>+10.6f} "
              f"{str(internal_min):>13}")

    # Verdict
    f3 = summary['W1c+mz+xblockh1']
    f3_internal_min_OK = f3['internal_min']
    f3_at_R3p5 = abs(f3['R_min'] - 3.566) < 1.0  # within 1 bohr of experimental
    f3_in_window = 0.5 * 0.075 <= f3['D_e_well_depth_Ha'] <= 2.0 * 0.075

    if f3_internal_min_OK and f3_at_R3p5 and f3_in_window:
        verdict = "WALL_CLOSES_GENUINE_BINDING"
    elif f3_internal_min_OK and f3_in_window:
        verdict = "WALL_CLOSES_R_OFF"
    elif f3_internal_min_OK:
        verdict = "PARTIAL_BINDING_MAG_OFF"
    elif f3['D_e_well_depth_Ha'] > 0:
        verdict = "PARTIAL_NO_INTERNAL_MIN"
    else:
        verdict = "WALL_PERSISTS_NO_BINDING"

    print(f"\n{'=' * 78}")
    print(f"Final F3 verdict (mini-PES): {verdict}")
    print(f"  R_min = {f3['R_min']:.2f} bohr (experimental NaH: 3.566)")
    print(f"  D_e = {f3['D_e_well_depth_Ha']:+.4f} Ha (experimental: 0.075)")
    print(f"  Internal minimum: {f3_internal_min_OK}")
    print(f"  R close to R_eq: {f3_at_R3p5}")
    print(f"  D_e in 2x window [0.0375, 0.15]: {f3_in_window}")
    print(f"{'=' * 78}")

    out = {
        'sprint': 'F3 Step 4 — mini-PES',
        'date': '2026-05-23',
        'R_grid': R_grid,
        'architectures': {k: results[k] for k in results},
        'summary': summary,
        'final_verdict': verdict,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
