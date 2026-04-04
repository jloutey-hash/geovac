"""Track AA: Fast convergence scan - adiabatic l_max=1-6, 2D l_max=1-6.

Uses reduced R-grid for expensive 2D runs at high l_max.
"""
import json
import time
import sys
import numpy as np
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.cusp_correction import cusp_correction_h2_pes

D_E_EXACT = 0.17447  # Ha (Kolos & Wolniewicz)

R_GRID_FULL = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0, 10.0])
R_GRID_REDUCED = np.array([1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 10.0])

def run_pes(l_max, method, R_grid):
    n_coupled = 1 if method == 'adiabatic' else -1
    n_alpha = 60 if method == '2d' else 200
    n_Re = 120 if method == '2d' else 400

    points = []
    t_total = 0
    for R in R_grid:
        t0 = time.time()
        res = solve_level4_h2_multichannel(
            R=R, l_max=l_max, m_max=1,
            l_max_per_m={0: l_max, 1: min(l_max, 2)},
            n_coupled=n_coupled, n_alpha=n_alpha, n_Re=n_Re,
            angular_method='spectral', n_basis_angular=10,
            n_quad_angular=100, verbose=False
        )
        dt = time.time() - t0
        t_total += dt
        points.append({'R': float(R), 'E_total': float(res['E_total']), 'time': float(dt)})
        print(f"  R={R:.1f}: E={res['E_total']:.6f} ({dt:.1f}s)", flush=True)

    R_arr = np.array([p['R'] for p in points])
    E_arr = np.array([p['E_total'] for p in points])
    E_inf = E_arr[-1] if R_arr[-1] > 8.0 else -1.0
    D_e = E_inf - np.min(E_arr)
    D_e_pct = D_e / D_E_EXACT * 100

    cusp = cusp_correction_h2_pes(R_arr, E_arr, l_max=l_max)

    return {
        'l_max': l_max, 'method': method,
        'D_e': float(D_e), 'D_e_pct': float(D_e_pct),
        'D_e_cusp_corrected': float(cusp['D_e_corrected']),
        'D_e_pct_cusp_corrected': float(cusp['D_e_pct_corrected']),
        'E_min': float(np.min(E_arr)),
        'R_eq_approx': float(R_arr[np.argmin(E_arr)]),
        'wall_time_total': float(t_total),
        'R_grid': R_arr.tolist(),
        'E_pes': E_arr.tolist(),
        'delta_E_cusp': cusp['delta_E_cusp'].tolist(),
        'points': points,
    }

if __name__ == '__main__':
    all_results = {}

    # Phase 1: All adiabatic (fast: ~4s/point * 60 = ~4 min)
    print("=== PHASE 1: ADIABATIC SCANS ===", flush=True)
    for lm in [1, 2, 3, 4, 5, 6]:
        print(f"\nAdiabatic l_max={lm}:", flush=True)
        res = run_pes(lm, 'adiabatic', R_GRID_FULL)
        all_results[f'adiabatic_lmax{lm}'] = res
        print(f"  -> D_e = {res['D_e_pct']:.2f}%, cusp = {res['D_e_pct_cusp_corrected']:.2f}%", flush=True)

    # Phase 2: 2D scans (slower)
    print("\n=== PHASE 2: 2D SCANS ===", flush=True)
    for lm in [1, 2, 3, 4, 5, 6]:
        # Use full grid for l_max<=4, reduced for l_max>=5
        grid = R_GRID_FULL if lm <= 4 else R_GRID_REDUCED
        print(f"\n2D l_max={lm} ({len(grid)} points):", flush=True)
        res = run_pes(lm, '2d', grid)
        all_results[f'2d_lmax{lm}'] = res
        print(f"  -> D_e = {res['D_e_pct']:.2f}%, cusp = {res['D_e_pct_cusp_corrected']:.2f}%", flush=True)

    # Summary
    print("\n" + "="*95)
    print(f"{'Method':<12} {'l_max':>5} {'N_pts':>5} {'D_e%':>8} {'D_e%+cusp':>10} {'E_min':>10} {'R_eq':>6} {'Time':>8}")
    print("="*95)
    for key in sorted(all_results.keys()):
        r = all_results[key]
        print(f"{r['method']:<12} {r['l_max']:>5} {len(r['R_grid']):>5} {r['D_e_pct']:>8.2f} {r['D_e_pct_cusp_corrected']:>10.2f} {r['E_min']:>10.6f} {r['R_eq_approx']:>6.1f} {r['wall_time_total']:>7.1f}s")

    with open('debug/track_aa/convergence_data.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\nSaved to debug/track_aa/convergence_data.json")
