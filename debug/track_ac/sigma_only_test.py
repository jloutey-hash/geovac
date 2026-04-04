"""Track AC Part 1e: Sigma-only vs sigma+pi comparison.

Runs l_max=4 and l_max=6 with m_max=0 (sigma only) and compares
to Track AA sigma+frozen-pi results.

If the staircase is from selection rules (frozen pi), sigma-only
should show a similar even-odd pattern.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import json
import time
from geovac.level4_multichannel import solve_level4_h2_multichannel

D_E_EXACT = 0.17447  # Ha (Kolos & Wolniewicz)

R_GRID = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0, 10.0])


def run_sigma_only(l_max, method='2d'):
    """Run PES with sigma-only channels."""
    n_coupled = 1 if method == 'adiabatic' else -1
    n_alpha = 60 if method == '2d' else 200
    n_Re = 120 if method == '2d' else 400

    points = []
    t0_total = time.time()
    for R in R_GRID:
        t0 = time.time()
        res = solve_level4_h2_multichannel(
            R=R, l_max=l_max, m_max=0,  # sigma only
            n_coupled=n_coupled, n_alpha=n_alpha, n_Re=n_Re,
            angular_method='spectral', n_basis_angular=10,
            n_quad_angular=100, verbose=False
        )
        dt = time.time() - t0
        points.append({'R': float(R), 'E_total': float(res['E_total']), 'time': float(dt)})
        print(f"  R={R:.1f}: E={res['E_total']:.6f} ({dt:.1f}s)", flush=True)

    R_arr = np.array([p['R'] for p in points])
    E_arr = np.array([p['E_total'] for p in points])
    E_inf = E_arr[-1]
    D_e = E_inf - np.min(E_arr)
    D_e_pct = D_e / D_E_EXACT * 100

    return {
        'l_max': l_max, 'method': method, 'm_max': 0,
        'D_e': float(D_e), 'D_e_pct': float(D_e_pct),
        'E_min': float(np.min(E_arr)),
        'E_inf': float(E_inf),
        'wall_time': float(time.time() - t0_total),
        'points': points,
    }


if __name__ == '__main__':
    results = {}

    for l_max in [2, 3, 4, 5, 6]:
        print(f"\n--- Sigma-only 2D l_max={l_max} ---", flush=True)
        res = run_sigma_only(l_max, method='2d')
        results[f'sigma_2d_lmax{l_max}'] = res
        print(f"  D_e = {res['D_e_pct']:.2f}%", flush=True)

    # Summary
    print("\n" + "=" * 70)
    print("SIGMA-ONLY RESULTS (m_max=0)")
    print("=" * 70)
    print(f"{'l_max':>5} | {'N_ch':>5} | {'D_e%':>8} | {'E_min':>10} | {'E_inf':>10} | {'Time':>8}")
    print("-" * 70)

    from geovac.level4_multichannel import _channel_list
    for l_max in [2, 3, 4, 5, 6]:
        key = f'sigma_2d_lmax{l_max}'
        r = results[key]
        n_ch = len(_channel_list(l_max, homonuclear=True))
        print(f"{l_max:>5} | {n_ch:>5} | {r['D_e_pct']:>8.2f} | {r['E_min']:>10.6f} | "
              f"{r['E_inf']:>10.6f} | {r['wall_time']:>7.1f}s")

    # Compare to Track AA frozen-pi values
    print("\nComparison with Track AA frozen-pi results:")
    print(f"{'l_max':>5} | {'sigma-only':>12} | {'sigma+frozen-pi':>15} | {'diff':>8}")
    print("-" * 50)
    # Track AA 2D values from results.md
    track_aa = {2: 86.44, 3: 87.09, 4: 93.87, 5: 93.95, 6: 95.64}
    for l_max in [2, 3, 4, 5, 6]:
        key = f'sigma_2d_lmax{l_max}'
        sigma = results[key]['D_e_pct']
        frozen = track_aa[l_max]
        print(f"{l_max:>5} | {sigma:>12.2f} | {frozen:>15.2f} | {frozen - sigma:>8.2f}")

    with open('debug/track_ac/sigma_only_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved to debug/track_ac/sigma_only_results.json")
