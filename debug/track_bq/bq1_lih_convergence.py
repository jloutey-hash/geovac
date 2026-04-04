"""
Track BQ-1: LiH composed PES l_max convergence study with variational 2D solver.
"""
import sys
import os
import json
import time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# Force unbuffered
sys.stdout.reconfigure(line_buffering=True)

from geovac.composed_diatomic import ComposedDiatomicSolver

R_EQ_REF = 3.015
D_E_REF = 0.0920

R_GRID_FINE = np.array([
    2.4, 2.6, 2.7, 2.8, 2.9, 3.0, 3.05, 3.1, 3.15, 3.2,
    3.3, 3.5, 4.0, 5.0, 6.0, 8.0
])

R_GRID_COARSE = np.array([2.8, 3.0, 3.2, 3.5, 5.0, 8.0])

OUTDIR = os.path.dirname(os.path.abspath(__file__))
ALL_RESULTS = []


def save_intermediate():
    outpath = os.path.join(OUTDIR, 'bq1_results.json')
    output = {
        'description': 'LiH composed PES l_max convergence with variational 2D solver',
        'reference': {'R_eq': R_EQ_REF, 'D_e': D_E_REF},
        'results': ALL_RESULTS,
    }
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2)


def analyze_pes(R_valid, E_valid):
    """Spline fit to find R_eq and D_e."""
    if len(R_valid) >= 4:
        cs = CubicSpline(R_valid, E_valid)
        i_min = np.argmin(E_valid)
        R_lo = R_valid[max(0, i_min - 2)]
        R_hi = R_valid[min(len(R_valid) - 1, i_min + 2)]
        opt = minimize_scalar(cs, bounds=(R_lo, R_hi), method='bounded')
        R_eq = opt.x
        E_min = opt.fun
    else:
        R_eq = R_valid[np.argmin(E_valid)]
        E_min = np.min(E_valid)
    E_dissoc = E_valid[-1]
    D_e = E_dissoc - E_min
    return R_eq, E_min, D_e, E_dissoc


def run_lmax(l_max, r_grid, n_Re=200):
    print(f"\n{'='*60}", flush=True)
    print(f"  l_max = {l_max}", flush=True)
    print(f"{'='*60}", flush=True)

    t_start = time.time()

    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=l_max,
        level4_method='variational_2d',
        pk_channel_mode='l_dependent',
    )
    solver.solve_core()
    t_core = time.time() - t_start
    print(f"  Core solve: {t_core:.1f}s", flush=True)

    result = solver.scan_pes(r_grid, n_Re=n_Re)

    R_valid = result['R_valid']
    E_valid = result['E_valid']
    R_eq, E_min, D_e, E_dissoc = analyze_pes(R_valid, E_valid)

    t_total = time.time() - t_start
    wt = result['wall_times']
    avg_t = float(np.mean(wt[wt > 0]))

    R_eq_err = abs(R_eq - R_EQ_REF) / R_EQ_REF * 100
    D_e_err = abs(D_e - D_E_REF) / D_E_REF * 100

    out = {
        'l_max': l_max,
        'R_eq_spline': float(R_eq),
        'E_min_spline': float(E_min),
        'R_eq_error_pct': float(R_eq_err),
        'D_e': float(D_e),
        'D_e_error_pct': float(D_e_err),
        'E_dissoc': float(E_dissoc),
        'wall_time_total': float(t_total),
        'wall_time_per_pt_avg': float(avg_t),
        'wall_times': wt.tolist(),
        'R_valid': R_valid.tolist(),
        'E_valid': E_valid.tolist(),
    }

    print(f"\n  --- l_max={l_max} Summary ---", flush=True)
    print(f"  R_eq = {R_eq:.4f} bohr  (err {R_eq_err:.2f}%)", flush=True)
    print(f"  D_e  = {D_e:.4f} Ha     (err {D_e_err:.1f}%)", flush=True)
    print(f"  E_min= {E_min:.6f} Ha", flush=True)
    print(f"  Avg t/pt = {avg_t:.1f}s, total = {t_total:.1f}s", flush=True)

    return out


if __name__ == '__main__':
    for l_max in [2, 3, 4]:
        try:
            res = run_lmax(l_max, R_GRID_FINE, n_Re=200)
            ALL_RESULTS.append(res)
            save_intermediate()

            # If l_max=4 was fast enough, try 5
            if l_max == 4 and res['wall_time_per_pt_avg'] < 60:
                print(f"\n  l_max=4 avg {res['wall_time_per_pt_avg']:.0f}s/pt, trying l_max=5", flush=True)
                res5 = run_lmax(5, R_GRID_FINE, n_Re=200)
                ALL_RESULTS.append(res5)
                save_intermediate()
            elif l_max == 4:
                print(f"\n  l_max=4 avg {res['wall_time_per_pt_avg']:.0f}s/pt, skipping l_max=5", flush=True)
        except Exception as e:
            print(f"\n  l_max={l_max} FAILED: {e}", flush=True)
            import traceback; traceback.print_exc()

    # CBS extrapolation
    if len(ALL_RESULTS) >= 3:
        l_v = np.array([r['l_max'] for r in ALL_RESULTS])
        r_v = np.array([r['R_eq_spline'] for r in ALL_RESULTS])
        d_v = np.array([r['D_e'] for r in ALL_RESULTS])
        x = 1.0 / (l_v + 0.5) ** 2
        cr = np.polyfit(x, r_v, 1)
        cd = np.polyfit(x, d_v, 1)
        R_CBS, D_CBS = cr[1], cd[1]
        cbs = {
            'R_eq_CBS': float(R_CBS),
            'R_eq_CBS_err_pct': float(abs(R_CBS - R_EQ_REF) / R_EQ_REF * 100),
            'D_e_CBS': float(D_CBS),
            'D_e_CBS_err_pct': float(abs(D_CBS - D_E_REF) / D_E_REF * 100),
        }
    else:
        cbs = {'feasible': False}

    # Convergence classification
    if len(ALL_RESULTS) >= 2:
        r_diffs = np.diff([r['R_eq_spline'] for r in ALL_RESULTS])
        d_diffs = np.diff([r['D_e'] for r in ALL_RESULTS])
        r_mono = all(d >= 0 for d in r_diffs) or all(d <= 0 for d in r_diffs)
        d_mono = all(d >= 0 for d in d_diffs) or all(d <= 0 for d in d_diffs)
        if r_mono and d_mono:
            if len(r_diffs) >= 2 and abs(r_diffs[-1]) < 0.3 * abs(r_diffs[0]):
                conv = 'PLATEAU'
            else:
                conv = 'MONOTONIC'
        else:
            conv = 'OSCILLATORY'
    else:
        conv = 'INSUFFICIENT_DATA'

    # Final summary
    print(f"\n{'='*80}", flush=True)
    print(f"  LiH Composed PES: l_max Convergence (variational 2D, l-dependent PK)", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"  {'l_max':>5}  {'R_eq':>8}  {'R_err%':>7}  {'D_e':>8}  {'D_err%':>7}  {'t/pt(s)':>8}", flush=True)
    print(f"  {'-'*5}  {'-'*8}  {'-'*7}  {'-'*8}  {'-'*7}  {'-'*8}", flush=True)
    for r in ALL_RESULTS:
        print(f"  {r['l_max']:5d}  {r['R_eq_spline']:8.4f}"
              f"  {r['R_eq_error_pct']:7.2f}  {r['D_e']:8.4f}"
              f"  {r['D_e_error_pct']:7.1f}  {r['wall_time_per_pt_avg']:8.1f}", flush=True)
    print(f"  {'ref':>5}  {R_EQ_REF:8.4f}  {'':>7}  {D_E_REF:8.4f}", flush=True)
    print(f"\n  Convergence: {conv}", flush=True)

    if 'R_eq_CBS' in cbs:
        print(f"\n  CBS (1/(l+1/2)^2):", flush=True)
        print(f"    R_eq_CBS = {cbs['R_eq_CBS']:.4f} (err {cbs['R_eq_CBS_err_pct']:.2f}%)", flush=True)
        print(f"    D_e_CBS  = {cbs['D_e_CBS']:.4f} (err {cbs['D_e_CBS_err_pct']:.1f}%)", flush=True)

    # Final save
    final = {
        'description': 'LiH composed PES l_max convergence with variational 2D solver',
        'reference': {'R_eq': R_EQ_REF, 'D_e': D_E_REF},
        'results': ALL_RESULTS,
        'cbs_extrapolation': cbs,
        'convergence_class': conv,
    }
    outpath = os.path.join(OUTDIR, 'bq1_results.json')
    with open(outpath, 'w') as f:
        json.dump(final, f, indent=2)
    print(f"\n  Saved to {outpath}", flush=True)
