"""
l_max convergence study: algebraic PK vs Gaussian PK modes for LiH.

Runs LiH PES scans at l_max = 2, 3, 4 (optionally 5) for three PK modes:
  - channel_blind: Gaussian PK applied to all channels
  - l_dependent:   Gaussian PK with delta_{l,0}
  - algebraic:     Rank-1 Phillips-Kleinman projector

Reports R_eq, drift rate, and channel weights for each configuration.
"""

import csv
import json
import os
import sys
import time
import numpy as np
from scipy.stats import linregress

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.composed_diatomic import ComposedDiatomicSolver
from geovac.level4_multichannel import _channel_list

# Constants
R_EQ_EXPT = 3.015  # bohr

# PES scan grid: ~18 points over R in [2.0, 7.0]
R_GRID = np.concatenate([
    np.linspace(2.0, 2.5, 3),
    np.linspace(2.7, 4.0, 10),
    np.linspace(4.5, 7.0, 5),
])

# Solver parameters
N_ALPHA = 100
N_RE = 300

# Output paths
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
CSV_PATH = os.path.join(DATA_DIR, 'lmax_convergence_algebraic.csv')
JSON_PATH = os.path.join(DATA_DIR, 'lmax_convergence_algebraic.json')
REPORT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'algebraic_pk_lmax_report.md')


def make_solver(pk_mode: str, l_max: int) -> ComposedDiatomicSolver:
    """Create a LiH solver for the given PK mode and l_max."""
    if pk_mode == 'channel_blind':
        return ComposedDiatomicSolver.LiH_ab_initio(
            l_max=l_max,
            pk_channel_mode='channel_blind',
            n_alpha=N_ALPHA,
            verbose=False,
        )
    elif pk_mode == 'l_dependent':
        return ComposedDiatomicSolver.LiH_ab_initio(
            l_max=l_max,
            pk_channel_mode='l_dependent',
            n_alpha=N_ALPHA,
            verbose=False,
        )
    elif pk_mode == 'algebraic':
        return ComposedDiatomicSolver.LiH_algebraic_pk(
            l_max=l_max,
            n_alpha=N_ALPHA,
            verbose=False,
        )
    else:
        raise ValueError(f"Unknown pk_mode: {pk_mode}")


def run_single(pk_mode: str, l_max: int) -> dict:
    """Run a single (mode, l_max) configuration and return results."""
    channels = _channel_list(l_max, homonuclear=False)
    n_ch = len(channels)

    print(f"\n{'='*60}")
    print(f"  pk_mode={pk_mode}, l_max={l_max}, n_channels={n_ch}")
    print(f"{'='*60}")

    t0 = time.time()

    solver = make_solver(pk_mode, l_max)
    solver.solve_core()

    t_core = time.time() - t0
    print(f"  Core solved in {t_core:.1f}s")
    print(f"  E_core = {solver.E_core:.6f} Ha")
    if pk_mode == 'algebraic':
        print(f"  PK projector energy_shift = "
              f"{solver.pk_projector['energy_shift']:.4f} Ha")
    else:
        print(f"  Gaussian PK: A={solver.pk_A:.4f}, B={solver.pk_B:.4f}")

    # PES scan
    try:
        pes = solver.scan_pes(R_grid=R_GRID, n_Re=N_RE)
    except Exception as e:
        print(f"  PES SCAN FAILED: {e}")
        wall_time = time.time() - t0
        return {
            'pk_mode': pk_mode,
            'l_max': l_max,
            'n_channels': n_ch,
            'R_eq_bohr': float('nan'),
            'R_eq_error_pct': float('nan'),
            'E_min_Ha': float('nan'),
            'wall_time_s': wall_time,
            'status': f'PES_FAILED: {e}',
            'pes_R': [],
            'pes_E': [],
        }

    # Check if PES has a minimum
    R_valid = pes['R_valid']
    E_valid = pes['E_valid']
    i_min = np.argmin(E_valid)

    if i_min == 0 or i_min == len(E_valid) - 1:
        print(f"  WARNING: Minimum at boundary (i={i_min}), "
              f"R={R_valid[i_min]:.3f})")
        # Still report the grid minimum

    # Morse fit for R_eq
    try:
        spectro = solver.fit_spectroscopic_constants()
        R_eq = spectro['R_eq']
        E_min = spectro['E_min']
    except Exception as e:
        print(f"  Morse fit failed: {e}")
        R_eq = pes['R_eq']  # grid minimum
        E_min = pes['E_min']

    R_eq_err = 100.0 * (R_eq - R_EQ_EXPT) / R_EQ_EXPT
    wall_time = time.time() - t0

    print(f"\n  R_eq = {R_eq:.4f} bohr  (error: {R_eq_err:+.1f}%)")
    print(f"  E_min = {E_min:.6f} Ha")
    print(f"  Wall time: {wall_time:.1f}s")

    result = {
        'pk_mode': pk_mode,
        'l_max': l_max,
        'n_channels': n_ch,
        'R_eq_bohr': R_eq,
        'R_eq_error_pct': R_eq_err,
        'E_min_Ha': E_min,
        'wall_time_s': wall_time,
        'status': 'OK',
        'pes_R': R_valid.tolist(),
        'pes_E': E_valid.tolist(),
    }

    return result


def compute_drift_rates(results: list) -> dict:
    """Compute linear drift rate (bohr/l_max) for each PK mode."""
    drift = {}
    for mode in ['channel_blind', 'l_dependent', 'algebraic']:
        pts = [(r['l_max'], r['R_eq_bohr'])
               for r in results
               if r['pk_mode'] == mode and r['status'] == 'OK'
               and np.isfinite(r['R_eq_bohr'])]
        if len(pts) < 2:
            drift[mode] = {'slope': float('nan'), 'R2': float('nan'),
                           'n_points': len(pts)}
            continue
        lm = np.array([p[0] for p in pts])
        req = np.array([p[1] for p in pts])
        slope, intercept, r_value, p_value, std_err = linregress(lm, req)
        drift[mode] = {
            'slope': slope,
            'R2': r_value**2,
            'intercept': intercept,
            'n_points': len(pts),
            'lm_range': f"{int(min(lm))}-{int(max(lm))}",
        }
    return drift


def write_csv(results: list) -> None:
    """Write results to CSV."""
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(CSV_PATH, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'pk_mode', 'l_max', 'n_channels', 'R_eq_bohr',
            'R_eq_error_pct', 'E_min_Ha', 'wall_time_s',
        ])
        w.writeheader()
        for r in results:
            w.writerow({
                'pk_mode': r['pk_mode'],
                'l_max': r['l_max'],
                'n_channels': r['n_channels'],
                'R_eq_bohr': f"{r['R_eq_bohr']:.4f}" if np.isfinite(r['R_eq_bohr']) else 'NaN',
                'R_eq_error_pct': f"{r['R_eq_error_pct']:.1f}" if np.isfinite(r['R_eq_error_pct']) else 'NaN',
                'E_min_Ha': f"{r['E_min_Ha']:.6f}" if np.isfinite(r['E_min_Ha']) else 'NaN',
                'wall_time_s': f"{r['wall_time_s']:.1f}",
            })
    print(f"\nCSV written to {CSV_PATH}")


def write_report(results: list, drift: dict) -> None:
    """Write diagnostic report."""
    lines = []
    lines.append("# Algebraic PK l_max Convergence Report\n")
    lines.append(f"**Date:** {time.strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"**n_alpha:** {N_ALPHA}")
    lines.append(f"**n_Re:** {N_RE}")
    lines.append(f"**R_grid:** {len(R_GRID)} points over "
                 f"[{R_GRID[0]:.1f}, {R_GRID[-1]:.1f}] bohr")
    lines.append(f"**Reference R_eq:** {R_EQ_EXPT} bohr (experiment)\n")

    # Full data table
    lines.append("## 1. Full Data Table\n")
    lines.append("| PK Mode | l_max | n_ch | R_eq (bohr) | R_eq error (%) "
                 "| E_min (Ha) | Wall time (s) |")
    lines.append("|:--------|:-----:|:----:|:-----------:|:--------------:"
                 "|:----------:|:-------------:|")
    for r in results:
        if r['status'] != 'OK':
            lines.append(f"| {r['pk_mode']} | {r['l_max']} | "
                         f"{r['n_channels']} | FAILED | — | — | "
                         f"{r['wall_time_s']:.0f} |")
        else:
            lines.append(f"| {r['pk_mode']} | {r['l_max']} | "
                         f"{r['n_channels']} | {r['R_eq_bohr']:.4f} | "
                         f"{r['R_eq_error_pct']:+.1f}% | "
                         f"{r['E_min_Ha']:.6f} | "
                         f"{r['wall_time_s']:.0f} |")

    # Drift rate comparison
    lines.append("\n## 2. Drift Rate Comparison\n")
    lines.append("Linear fit: R_eq = slope × l_max + intercept\n")
    lines.append("| PK Mode | Slope (bohr/l_max) | R² | l_max range | "
                 "n_points |")
    lines.append("|:--------|:------------------:|:--:|:-----------:|"
                 ":--------:|")

    ref_drift = {
        'channel_blind': '+0.23',
        'l_dependent': '+0.168',
    }

    for mode in ['channel_blind', 'l_dependent', 'algebraic']:
        d = drift[mode]
        if np.isfinite(d['slope']):
            ref = ref_drift.get(mode, '—')
            lines.append(f"| {mode} | {d['slope']:+.3f} | "
                         f"{d['R2']:.3f} | {d['lm_range']} | "
                         f"{d['n_points']} |")
        else:
            lines.append(f"| {mode} | — | — | — | {d['n_points']} |")

    lines.append("\n**Reference drift rates (Paper 17 / diagnostic arc):**")
    lines.append("- channel_blind: +0.23 bohr/l_max (l_max 2–5, "
                 "diagnostic arc)")
    lines.append("- l_dependent: +0.168 bohr/l_max (l_max 2–7, Paper 17)")
    lines.append("- r_dependent: +0.160 bohr/l_max (l_max 2–7, Paper 17)\n")

    # Assessment
    lines.append("## 3. Assessment\n")

    alg_drift = drift.get('algebraic', {})
    ldep_drift = drift.get('l_dependent', {})
    cb_drift = drift.get('channel_blind', {})

    if np.isfinite(alg_drift.get('slope', float('nan'))):
        alg_slope = alg_drift['slope']
        ldep_slope = ldep_drift.get('slope', float('nan'))

        if np.isfinite(ldep_slope) and ldep_slope > 0:
            ratio = alg_slope / ldep_slope
            if abs(alg_slope) < 0.02:
                verdict = ("**ELIMINATED.** The algebraic PK projector "
                           "eliminates the l_max divergence "
                           f"(drift = {alg_slope:+.3f} bohr/l_max).")
            elif ratio < 0.5:
                verdict = (f"**REDUCED.** The algebraic PK projector "
                           f"reduces the drift by "
                           f"{(1-ratio)*100:.0f}% "
                           f"(from {ldep_slope:+.3f} to "
                           f"{alg_slope:+.3f} bohr/l_max).")
            else:
                verdict = (f"**NOT ELIMINATED.** The algebraic PK "
                           f"projector drift ({alg_slope:+.3f}) is "
                           f"similar to l_dependent ({ldep_slope:+.3f}).")
        else:
            verdict = (f"Algebraic drift = {alg_slope:+.3f} bohr/l_max. "
                       f"l_dependent reference unavailable for comparison.")
    else:
        verdict = ("Insufficient data to compute algebraic drift rate.")

    lines.append(verdict)

    # PES data for each configuration
    lines.append("\n## 4. PES Data\n")
    for r in results:
        if r['status'] == 'OK' and r.get('pes_R'):
            lines.append(f"### {r['pk_mode']}, l_max={r['l_max']}\n")
            lines.append("| R (bohr) | E_composed (Ha) |")
            lines.append("|:--------:|:---------------:|")
            for R_val, E_val in zip(r['pes_R'], r['pes_E']):
                lines.append(f"| {R_val:.3f} | {E_val:.6f} |")
            lines.append("")

    with open(REPORT_PATH, 'w') as f:
        f.write('\n'.join(lines))
    print(f"\nReport written to {REPORT_PATH}")


def main():
    print("=" * 64)
    print("  LiH l_max Convergence: Algebraic PK vs Gaussian PK")
    print("=" * 64)

    l_max_values = [2, 3, 4]
    pk_modes = ['channel_blind', 'l_dependent', 'algebraic']

    results = []

    # Step 1: Validation at l_max=2
    print("\n" + "=" * 64)
    print("  STEP 1: Validation at l_max=2")
    print("=" * 64)

    for mode in ['channel_blind', 'l_dependent']:
        r = run_single(mode, l_max=2)
        results.append(r)

        if r['status'] != 'OK':
            print(f"\n  VALIDATION FAILED: {mode} at l_max=2 did not produce "
                  f"a valid PES. Aborting.")
            write_csv(results)
            return

    # Validate against known results
    cb_req = results[0]['R_eq_bohr']
    ld_req = results[1]['R_eq_bohr']

    print(f"\n  Validation results (l_max=2):")
    print(f"    channel_blind: R_eq = {cb_req:.4f} bohr "
          f"(expect ~3.21, Paper 17)")
    print(f"    l_dependent:   R_eq = {ld_req:.4f} bohr "
          f"(expect ~3.17, Paper 17)")

    # Sanity check: R_eq should be in [2.0, 5.0] range
    for r in results:
        req = r['R_eq_bohr']
        if not (2.0 < req < 5.0):
            print(f"\n  WARNING: {r['pk_mode']} R_eq = {req:.3f} is "
                  f"outside expected range [2.0, 5.0]. Check parameters.")

    # Step 2: Run algebraic at l_max=2
    print("\n" + "=" * 64)
    print("  STEP 2: Algebraic PK at l_max=2")
    print("=" * 64)

    r = run_single('algebraic', l_max=2)
    results.append(r)

    if r['status'] != 'OK':
        print(f"\n  Algebraic PK FAILED at l_max=2.")
        if r.get('pes_R'):
            print(f"  Raw PES data:")
            for R_val, E_val in zip(r['pes_R'], r['pes_E']):
                print(f"    R={R_val:.3f}  E={E_val:.6f}")

    # Step 3: Run all modes at l_max=3, 4
    for l_max in [3, 4]:
        print(f"\n{'='*64}")
        print(f"  STEP 3: l_max={l_max} (all modes)")
        print(f"{'='*64}")

        # Time estimate from previous l_max
        prev_times = [r['wall_time_s'] for r in results
                      if r['l_max'] == l_max - 1 and r['status'] == 'OK']
        if prev_times:
            est = max(prev_times) * 2.5  # rough scaling
            print(f"  Estimated time per mode: {est:.0f}s")

        for mode in pk_modes:
            r = run_single(mode, l_max)
            results.append(r)

            # Save intermediate results
            write_csv(results)

        # Check if l_max=4 took too long for l_max=5
        if l_max == 4:
            t4 = max(r['wall_time_s'] for r in results
                     if r['l_max'] == 4 and r['status'] == 'OK')
            if t4 > 1800:  # 30 minutes
                print(f"\n  l_max=4 took {t4:.0f}s (> 30 min). "
                      f"Skipping l_max=5.")
            else:
                print(f"\n  l_max=4 took {t4:.0f}s. "
                      f"l_max=5 would take ~{t4*2.5:.0f}s.")
                # Optionally run l_max=5
                if t4 < 600:  # Only if l_max=4 was under 10 min
                    print(f"  Running l_max=5...")
                    for mode in pk_modes:
                        r = run_single(mode, l_max=5)
                        results.append(r)
                        write_csv(results)
                else:
                    print(f"  Skipping l_max=5 (too slow).")

    # Step 4: Compute drift rates
    drift = compute_drift_rates(results)

    print("\n" + "=" * 64)
    print("  DRIFT RATE SUMMARY")
    print("=" * 64)
    for mode in pk_modes:
        d = drift[mode]
        if np.isfinite(d['slope']):
            print(f"  {mode:20s}  slope = {d['slope']:+.3f} bohr/l_max"
                  f"  (R² = {d['R2']:.3f})")
        else:
            print(f"  {mode:20s}  insufficient data")

    # Step 5: Save results
    write_csv(results)

    # Save full JSON with PES data
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(JSON_PATH, 'w') as f:
        json.dump({
            'results': results,
            'drift_rates': drift,
            'parameters': {
                'n_alpha': N_ALPHA,
                'n_Re': N_RE,
                'R_grid': R_GRID.tolist(),
                'R_eq_expt': R_EQ_EXPT,
            },
        }, f, indent=2, default=str)
    print(f"\nJSON written to {JSON_PATH}")

    # Write report
    write_report(results, drift)

    print(f"\nDone. Total configurations: {len(results)}")


if __name__ == '__main__':
    main()
