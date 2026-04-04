"""
Track BP2: Test whether 2D variational solver eliminates l_max divergence
in LiH composed PES.

Compares adiabatic vs variational_2d at l_max=2, then optionally l_max=3.
"""

import json
import time
import sys
import os
import numpy as np

# Ensure project root is on path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from geovac.composed_diatomic import ComposedDiatomicSolver

# Reference data
REF_R_EQ = 3.015  # bohr
REF_D_E = 0.0920  # Ha

# Coarse R-grid
R_GRID = np.array([2.0, 2.5, 2.8, 3.0, 3.2, 3.5, 4.0, 5.0, 6.0])

N_RE = 200  # reduced from 300 for speed


def run_pes(l_max: int, method: str) -> dict:
    """Run LiH composed PES with given l_max and Level 4 method."""
    print(f"\n{'='*70}")
    print(f"  LiH composed PES: l_max={l_max}, level4_method='{method}'")
    print(f"  pk_channel_mode='l_dependent', n_Re={N_RE}")
    print(f"{'='*70}")

    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=l_max,
        level4_method=method,
        pk_channel_mode='l_dependent',
        verbose=True,
    )

    t0 = time.time()
    solver.solve_core()
    t_core = time.time() - t0
    print(f"  Core solve time: {t_core:.1f}s")

    t0 = time.time()
    pes = solver.scan_pes(R_grid=R_GRID, n_Re=N_RE)
    t_pes = time.time() - t0

    R_eq = pes['R_eq']
    E_min = pes['E_min']
    D_e = pes['D_e']
    wall_times = pes['wall_times']

    R_eq_err = abs(R_eq - REF_R_EQ) / REF_R_EQ * 100
    D_e_err = abs(D_e - REF_D_E) / REF_D_E * 100

    print(f"\n  Summary:")
    print(f"    R_eq = {R_eq:.3f} bohr (ref {REF_R_EQ}, err {R_eq_err:.1f}%)")
    print(f"    E_min = {E_min:.6f} Ha")
    print(f"    D_e = {D_e:.6f} Ha (ref {REF_D_E}, err {D_e_err:.1f}%)")
    print(f"    Total PES time: {t_pes:.1f}s, avg {wall_times.mean():.1f}s/pt")

    return {
        'l_max': l_max,
        'method': method,
        'R_eq': float(R_eq),
        'E_min': float(E_min),
        'D_e': float(D_e),
        'R_eq_err_pct': float(R_eq_err),
        'D_e_err_pct': float(D_e_err),
        'R_grid': R_GRID.tolist(),
        'E_composed': pes['E_composed'].tolist(),
        'E_elec': pes['E_elec'].tolist(),
        'wall_times': wall_times.tolist(),
        'avg_time_per_R': float(wall_times.mean()),
        'total_time': float(t_pes),
        'core_time': float(t_core),
    }


def main():
    results = {}

    # --- l_max=2: adiabatic ---
    results['lmax2_adiabatic'] = run_pes(l_max=2, method='adiabatic')

    # --- l_max=2: variational_2d ---
    results['lmax2_variational_2d'] = run_pes(l_max=2, method='variational_2d')

    # --- Comparison table ---
    print(f"\n{'='*70}")
    print(f"  COMPARISON: l_max=2 adiabatic vs variational_2d")
    print(f"{'='*70}")
    a = results['lmax2_adiabatic']
    v = results['lmax2_variational_2d']

    print(f"\n  {'R':>6s}  {'E_adiab':>12s}  {'E_var2d':>12s}  {'diff(mHa)':>10s}"
          f"  {'t_ad(s)':>8s}  {'t_2d(s)':>8s}")
    print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*8}  {'-'*8}")
    for i, R in enumerate(R_GRID):
        E_a = a['E_composed'][i]
        E_v = v['E_composed'][i]
        diff = (E_v - E_a) * 1000  # mHa
        t_a = a['wall_times'][i]
        t_v = v['wall_times'][i]
        print(f"  {R:6.2f}  {E_a:12.6f}  {E_v:12.6f}  {diff:+10.3f}"
              f"  {t_a:8.1f}  {t_v:8.1f}")

    print(f"\n  Adiabatic:       R_eq={a['R_eq']:.3f}  D_e={a['D_e']:.6f}  "
          f"R_eq err={a['R_eq_err_pct']:.1f}%  D_e err={a['D_e_err_pct']:.1f}%")
    print(f"  Variational 2D:  R_eq={v['R_eq']:.3f}  D_e={v['D_e']:.6f}  "
          f"R_eq err={v['R_eq_err_pct']:.1f}%  D_e err={v['D_e_err_pct']:.1f}%")

    # --- l_max=3 with variational_2d (timing test) ---
    print(f"\n{'='*70}")
    print(f"  l_max=3 variational_2d: timing test (single R-point first)")
    print(f"{'='*70}")

    # Time a single R-point first
    solver3 = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=3,
        level4_method='variational_2d',
        pk_channel_mode='l_dependent',
        verbose=True,
    )
    solver3.solve_core()

    t0 = time.time()
    try:
        E_test = solver3._solve_valence_at_R(3.0, n_Re=N_RE)
        t_single = time.time() - t0
        print(f"  Single R=3.0 point: {t_single:.1f}s, E_elec={E_test:.6f}")

        if t_single > 1800:  # 30 min
            print(f"  SKIP: {t_single:.0f}s > 1800s per point, full scan infeasible")
            results['lmax3_variational_2d'] = {
                'status': 'skipped',
                'reason': f'single point took {t_single:.0f}s (>30 min)',
                'single_point_time': float(t_single),
            }
        else:
            # Run full PES
            results['lmax3_variational_2d'] = run_pes(l_max=3, method='variational_2d')

            # Compute drift
            v2 = results['lmax2_variational_2d']
            v3 = results['lmax3_variational_2d']
            drift = v3['R_eq'] - v2['R_eq']
            print(f"\n  l_max drift (var2d): R_eq(3) - R_eq(2) = {drift:+.3f} bohr")
            results['drift_var2d_lmax2_to_3'] = float(drift)

    except Exception as e:
        t_single = time.time() - t0
        print(f"  FAILED after {t_single:.1f}s: {e}")
        results['lmax3_variational_2d'] = {
            'status': 'failed',
            'error': str(e),
            'time_before_failure': float(t_single),
        }

    # --- Save results ---
    out_path = os.path.join(os.path.dirname(__file__), 'bp2_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to {out_path}")


if __name__ == '__main__':
    main()
