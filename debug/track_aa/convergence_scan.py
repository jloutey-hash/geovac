"""Track AA: Level 4 H2 convergence scan at l_max=1-6.

Runs adiabatic and 2D solvers, applies cusp correction, saves results.
"""
import json
import time
import numpy as np
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.cusp_correction import cusp_correction_h2_pes

D_E_EXACT = 0.17447  # Ha (Kolos & Wolniewicz)
E_ATOMS = -1.0  # Two isolated H atoms

R_GRID = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0, 10.0])

def run_pes(l_max, method='adiabatic', label=''):
    """Run PES scan at given l_max."""
    n_coupled = 1 if method == 'adiabatic' else -1
    n_alpha = 60 if method == '2d' else 200
    n_Re = 120 if method == '2d' else 400

    results = []
    t_total = 0
    for i, R in enumerate(R_GRID):
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
        results.append({
            'R': float(R),
            'E_total': float(res['E_total']),
            'time': float(dt),
        })
        print(f"  {label} R={R:.1f}: E={res['E_total']:.6f}, {dt:.1f}s")

    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r['E_total'] for r in results])

    # D_e from PES
    E_inf = E_arr[-1] if R_arr[-1] > 8.0 else E_ATOMS
    D_e = E_inf - np.min(E_arr)
    D_e_pct = D_e / D_E_EXACT * 100

    # Cusp correction
    cusp = cusp_correction_h2_pes(R_arr, E_arr, l_max=l_max)

    return {
        'l_max': l_max,
        'method': method,
        'n_channels': results[0].get('n_ch', None),
        'R_grid': R_arr.tolist(),
        'E_pes': E_arr.tolist(),
        'D_e': float(D_e),
        'D_e_pct': float(D_e_pct),
        'D_e_cusp_corrected': float(cusp['D_e_corrected']),
        'D_e_pct_cusp_corrected': float(cusp['D_e_pct_corrected']),
        'delta_E_cusp': cusp['delta_E_cusp'].tolist(),
        'E_min': float(np.min(E_arr)),
        'R_eq_approx': float(R_arr[np.argmin(E_arr)]),
        'wall_time_total': float(t_total),
        'points': results,
    }


if __name__ == '__main__':
    all_results = {}

    # Adiabatic l_max = 1-6
    for lm in [1, 2, 3, 4, 5, 6]:
        label = f"Adiab l_max={lm}"
        print(f"\n=== {label} ===")
        res = run_pes(lm, method='adiabatic', label=label)
        all_results[f'adiabatic_lmax{lm}'] = res
        print(f"  D_e = {res['D_e_pct']:.2f}%, cusp-corrected = {res['D_e_pct_cusp_corrected']:.2f}%")

    # 2D l_max = 1-6
    for lm in [1, 2, 3, 4, 5, 6]:
        label = f"2D l_max={lm}"
        print(f"\n=== {label} ===")
        res = run_pes(lm, method='2d', label=label)
        all_results[f'2d_lmax{lm}'] = res
        print(f"  D_e = {res['D_e_pct']:.2f}%, cusp-corrected = {res['D_e_pct_cusp_corrected']:.2f}%")

    # Summary table
    print("\n" + "="*90)
    print(f"{'Method':<15} {'l_max':>5} {'D_e%':>8} {'D_e%+cusp':>10} {'E_min':>10} {'R_eq':>6} {'Time':>8}")
    print("="*90)
    for key in sorted(all_results.keys()):
        r = all_results[key]
        print(f"{r['method']:<15} {r['l_max']:>5} {r['D_e_pct']:>8.2f} {r['D_e_pct_cusp_corrected']:>10.2f} {r['E_min']:>10.6f} {r['R_eq_approx']:>6.1f} {r['wall_time_total']:>7.1f}s")

    # Save
    with open('debug/track_aa/convergence_data.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\nSaved to debug/track_aa/convergence_data.json")
