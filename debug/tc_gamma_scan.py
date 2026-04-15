"""
TC gamma parameter scan for He 2D variational solver.

Maps gamma_optimal(l_max) for l_max = 0..5 to understand the sweet spot.

Track DI investigation: why is gamma_optimal=0.11 << Kato value 1.0?
"""

import sys
import os
import time
import json
import numpy as np

# Ensure project root is on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.level3_variational import solve_he_tc_2d, solve_he_variational_2d

E_EXACT = -2.903724377034119598


def run_single_point(l_max, gamma, n_basis_R=20, n_basis_alpha=20):
    """Run TC solver at a single (l_max, gamma) point."""
    result = solve_he_tc_2d(
        Z=2.0,
        n_basis_R=n_basis_R,
        n_basis_alpha=n_basis_alpha,
        l_max=l_max,
        alpha_R=1.5,
        R_min=0.05,
        gamma_tc=gamma,
        symmetry='singlet',
        n_quad_alpha=100,
    )
    return result


def timing_calibration():
    """Time a single point at each l_max to decide scan strategy."""
    print("=" * 70)
    print("TIMING CALIBRATION")
    print("=" * 70)
    timings = {}
    for l_max in range(6):
        t0 = time.time()
        result = run_single_point(l_max, gamma=0.11)
        dt = time.time() - t0
        timings[l_max] = dt
        print(f"  l_max={l_max}: {dt:.2f}s  (dim={result['dim']}, "
              f"E_tc={result['E_tc']:.6f}, err={result['error_tc_pct']:.4f}%)")
    print()
    return timings


def gamma_scan(l_max, gammas, n_basis_R=20, n_basis_alpha=20):
    """Scan gamma at fixed l_max. Returns list of result dicts."""
    results = []
    for g in gammas:
        t0 = time.time()
        try:
            r = run_single_point(l_max, g, n_basis_R, n_basis_alpha)
            dt = time.time() - t0
            r['wall_time'] = dt
            results.append(r)
        except Exception as e:
            dt = time.time() - t0
            results.append({
                'gamma_tc': g, 'error_tc_pct': float('nan'),
                'E_tc': float('nan'), 'wall_time': dt, 'error': str(e),
                'l_max': l_max,
            })
    return results


def find_optimal(results):
    """Find the gamma with minimum TC error from scan results."""
    valid = [(r['gamma_tc'], r['error_tc_pct']) for r in results
             if np.isfinite(r['error_tc_pct'])]
    if not valid:
        return None, float('nan'), float('nan')
    best = min(valid, key=lambda x: x[1])
    # Also get the standard error (same for all gamma at fixed l_max)
    std_err = results[0].get('error_standard_pct', float('nan'))
    return best[0], best[1], std_err


def main():
    print("TC Gamma Parameter Scan for He 2D Variational Solver")
    print("=" * 70)
    print()

    # Phase 1: Timing calibration
    timings = timing_calibration()

    # Decide scan strategy based on timing
    max_time_per_point = max(timings.values())
    print(f"Max time per point: {max_time_per_point:.2f}s")

    if max_time_per_point > 120:
        print("WARNING: Very slow. Using ultra-coarse scan (5 points).")
        gammas = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
    elif max_time_per_point > 30:
        print("Using coarse scan (6 points).")
        gammas = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
    elif max_time_per_point > 10:
        print("Using medium scan (10 points).")
        gammas = np.arange(0.02, 0.32, 0.03)
    else:
        print("Fast enough for fine scan (30 points).")
        gammas = np.arange(0.01, 0.31, 0.01)

    print(f"Gamma values: {gammas}")
    print()

    # Phase 2: Full gamma scan at each l_max
    all_results = {}
    summary = []

    for l_max in range(6):
        print(f"--- l_max = {l_max} ---")
        t0 = time.time()
        results = gamma_scan(l_max, gammas)
        total_t = time.time() - t0
        print(f"  Completed in {total_t:.1f}s")

        g_opt, err_opt, err_std = find_optimal(results)
        print(f"  Standard error: {err_std:.5f}%")
        print(f"  Optimal gamma:  {g_opt:.3f}")
        print(f"  TC error:       {err_opt:.5f}%")
        print(f"  Improvement:    {err_std / max(err_opt, 1e-15):.1f}x")
        print()

        all_results[l_max] = results
        summary.append({
            'l_max': l_max,
            'gamma_optimal': g_opt,
            'error_standard_pct': err_std,
            'error_tc_optimal_pct': err_opt,
            'improvement': err_std / max(err_opt, 1e-15),
            'scan_time_s': total_t,
        })

    # Phase 3: Analysis
    print()
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'l_max':>5} {'gamma_opt':>10} {'std_err%':>10} {'tc_err%':>10} {'improve':>10}")
    print("-" * 50)
    for s in summary:
        print(f"{s['l_max']:>5} {s['gamma_optimal']:>10.4f} "
              f"{s['error_standard_pct']:>10.5f} "
              f"{s['error_tc_optimal_pct']:>10.5f} "
              f"{s['improvement']:>10.1f}x")

    print()
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    # Check for 1/(l_max + c) fit
    g_opts = [s['gamma_optimal'] for s in summary if s['gamma_optimal'] is not None]
    l_vals = [s['l_max'] for s in summary if s['gamma_optimal'] is not None]

    if len(g_opts) >= 3:
        # Fit gamma_opt = A / (l_max + c)
        # Try c = 1, 2, 3, ... and find best linear fit
        best_c = None
        best_r2 = -1
        for c_try in np.arange(0.5, 10.0, 0.1):
            x = 1.0 / (np.array(l_vals) + c_try)
            A_fit = np.sum(np.array(g_opts) * x) / np.sum(x * x)
            pred = A_fit * x
            ss_res = np.sum((np.array(g_opts) - pred) ** 2)
            ss_tot = np.sum((np.array(g_opts) - np.mean(g_opts)) ** 2)
            if ss_tot > 0:
                r2 = 1 - ss_res / ss_tot
            else:
                r2 = 0
            if r2 > best_r2:
                best_r2 = r2
                best_c = c_try
                best_A = A_fit

        print(f"\nBest fit: gamma_opt = {best_A:.4f} / (l_max + {best_c:.1f})")
        print(f"R^2 = {best_r2:.4f}")
        print()

    # Check TC convergence exponent at fixed gamma=0.11
    print("TC convergence at fixed gamma=0.11:")
    print(f"{'l_max':>5} {'std_err%':>10} {'tc_err%':>10}")
    for s in summary:
        print(f"{s['l_max']:>5} {s['error_standard_pct']:>10.5f} "
              f"{s['error_tc_optimal_pct']:>10.5f}")

    # Log-log fit for convergence exponent (skip l=0)
    l_fit = []
    std_fit = []
    tc_fit = []
    for l_max in range(6):
        res_011 = None
        for r in all_results[l_max]:
            if abs(r['gamma_tc'] - 0.11) < 0.005 or abs(r['gamma_tc'] - 0.10) < 0.005:
                res_011 = r
                break
        if res_011 is not None and l_max >= 1:
            l_fit.append(l_max)
            std_fit.append(res_011['error_standard_pct'])
            tc_fit.append(res_011['error_tc_pct'])

    if len(l_fit) >= 3:
        # Log-log fit: log(err) = a + b*log(l_max)
        log_l = np.log(np.array(l_fit))
        log_std = np.log(np.array(std_fit))
        log_tc = np.log(np.array(tc_fit))

        A_std = np.vstack([log_l, np.ones(len(log_l))]).T
        b_std, _ = np.linalg.lstsq(A_std, log_std, rcond=None)[:2]
        b_tc, _ = np.linalg.lstsq(A_std, log_tc, rcond=None)[:2]

        print(f"\nConvergence exponents (log-log fit, l_max >= 1):")
        print(f"  Standard: error ~ l_max^({b_std[0]:.2f})")
        print(f"  TC (gamma~0.11): error ~ l_max^({b_tc[0]:.2f})")

    # Check: at what l_max does TC become counterproductive?
    print("\nTC benefit analysis (optimal gamma at each l_max):")
    for s in summary:
        benefit = "BENEFICIAL" if s['error_tc_optimal_pct'] < s['error_standard_pct'] else "COUNTERPRODUCTIVE"
        print(f"  l_max={s['l_max']}: {benefit} "
              f"(std {s['error_standard_pct']:.5f}% vs TC {s['error_tc_optimal_pct']:.5f}%)")

    # Phase 4: Detailed gamma profile at each l_max
    print()
    print("=" * 70)
    print("DETAILED GAMMA PROFILES")
    print("=" * 70)
    for l_max in range(6):
        print(f"\nl_max = {l_max}:")
        print(f"  {'gamma':>8} {'TC_err%':>12} {'std_err%':>12}")
        for r in all_results[l_max]:
            print(f"  {r['gamma_tc']:>8.3f} {r['error_tc_pct']:>12.5f} "
                  f"{r.get('error_standard_pct', float('nan')):>12.5f}")

    # Save results
    save_data = {
        'summary': summary,
        'gammas': gammas.tolist(),
        'profiles': {},
    }
    for l_max in range(6):
        save_data['profiles'][str(l_max)] = [
            {'gamma': r['gamma_tc'], 'error_tc_pct': r['error_tc_pct'],
             'error_standard_pct': r.get('error_standard_pct', None),
             'E_tc': r['E_tc'], 'E_standard': r.get('E_standard', None)}
            for r in all_results[l_max]
        ]

    out_path = os.path.join(os.path.dirname(__file__), 'data', 'tc_gamma_scan.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nData saved to {out_path}")


if __name__ == '__main__':
    main()
