#!/usr/bin/env python3
"""
l_max divergence isolation experiment.

Tests whether the l_max divergence in LiH composed geometry originates from:
  (A) Intrinsic Level 4 solver behavior at asymmetric charges
  (B) Z_eff(r) injection breaking algebraic nuclear coupling
  (C) Static core / dynamic valence mismatch

Runs bare Level 4 solver (no Z_eff, no PK, no composed geometry) on four
systems with increasing charge asymmetry at l_max = 2, 3, 4.

System 1: H2 (Z_A=1, Z_B=1) — symmetric baseline
System 2: HeH+ (Z_A=2, Z_B=1) — 2:1 asymmetric
System 3: Effective LiH valence (Z_A=1, Z_B=1 at LiH distances, midpoint)
System 4: Effective LiH valence with charge-center origin (as if Z_A=3)
"""

import numpy as np
import json
import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.level4_multichannel import (
    solve_level4_h2_multichannel,
    solve_angular_multichannel,
)


def parabolic_fit_Req(R_grid: np.ndarray, E_grid: np.ndarray) -> tuple:
    """Find R_eq via parabolic fit near the minimum.

    Returns (R_eq, E_min, success).
    """
    idx_min = np.argmin(E_grid)

    # Need at least one point on each side
    if idx_min == 0 or idx_min == len(R_grid) - 1:
        return R_grid[idx_min], E_grid[idx_min], False

    # Use 3-5 points around minimum for parabolic fit
    lo = max(0, idx_min - 2)
    hi = min(len(R_grid), idx_min + 3)
    R_fit = R_grid[lo:hi]
    E_fit = E_grid[lo:hi]

    # Quadratic fit: E = a*R^2 + b*R + c
    coeffs = np.polyfit(R_fit, E_fit, 2)
    a, b, c = coeffs

    if a <= 0:  # No minimum (concave)
        return R_grid[idx_min], E_grid[idx_min], False

    R_eq = -b / (2 * a)
    E_min = c - b**2 / (4 * a)

    # Sanity: R_eq should be near the grid minimum
    if R_eq < R_grid[0] * 0.8 or R_eq > R_grid[-1] * 1.2:
        return R_grid[idx_min], E_grid[idx_min], False

    return R_eq, E_min, True


def extract_channel_weights(
    R: float, R_e: float, l_max: int,
    Z_A: float, Z_B: float, z0: float = 0.0,
    n_alpha: int = 200,
) -> dict:
    """Extract channel weight distribution at a specific (R, R_e) point."""
    rho = R / (2.0 * R_e)
    evals, vecs, alpha, channels = solve_angular_multichannel(
        rho, R_e, l_max=l_max, n_alpha=n_alpha,
        Z_A=Z_A, Z_B=Z_B, z0=z0,
    )

    # Ground state eigenvector: shape (n_ch * n_alpha,)
    vec0 = vecs[0]
    n_ch = len(channels)

    weights = {}
    for c_idx, ch in enumerate(channels):
        block = vec0[c_idx * n_alpha : (c_idx + 1) * n_alpha]
        w = np.sum(block**2) * (np.pi / 2) / (n_alpha + 1)
        weights[str(ch)] = float(w)

    total = sum(weights.values())
    # Normalize
    for k in weights:
        weights[k] /= total

    return weights


def run_pes_scan(
    system_name: str,
    Z_A: float,
    Z_B: float,
    R_grid: np.ndarray,
    l_max_values: list,
    origin: str = 'midpoint',
    n_alpha: int = 200,
    n_Re: int = 400,
) -> dict:
    """Run PES scan for a system at multiple l_max values."""
    results = {}

    for l_max in l_max_values:
        print(f"\n{'='*60}")
        print(f"System: {system_name}, l_max={l_max}, origin={origin}")
        print(f"Z_A={Z_A}, Z_B={Z_B}, R range=[{R_grid[0]:.2f}, {R_grid[-1]:.2f}]")
        print(f"{'='*60}")

        E_total = np.zeros(len(R_grid))
        t_start = time.time()

        for i, R in enumerate(R_grid):
            try:
                result = solve_level4_h2_multichannel(
                    R=R,
                    l_max=l_max,
                    Z_A=Z_A,
                    Z_B=Z_B,
                    n_alpha=n_alpha,
                    n_Re=n_Re,
                    origin=origin,
                    verbose=False,
                )
                E_total[i] = result['E_total']
                n_ch = result['n_ch']
            except Exception as e:
                print(f"  FAILED at R={R:.3f}: {e}")
                E_total[i] = np.nan

            print(f"  R={R:.3f}  E_total={E_total[i]:.6f} Ha")

        t_elapsed = time.time() - t_start

        # Extract R_eq
        valid = ~np.isnan(E_total)
        R_eq, E_min, success = parabolic_fit_Req(R_grid[valid], E_total[valid])

        print(f"\n  R_eq = {R_eq:.4f} bohr, E_min = {E_min:.6f} Ha")
        print(f"  Fit success: {success}")
        print(f"  Wall time: {t_elapsed:.1f}s")

        results[l_max] = {
            'R_grid': R_grid.tolist(),
            'E_total': E_total.tolist(),
            'R_eq': float(R_eq),
            'E_min': float(E_min),
            'fit_success': success,
            'n_ch': n_ch,
            'wall_time_s': round(t_elapsed, 1),
        }

    return results


def main():
    l_max_values = [2, 3, 4]

    # Match composed geometry parameters
    n_alpha = 200
    n_Re = 400

    all_results = {}

    # ---- System 1: H2 (symmetric baseline) ----
    R_grid_h2 = np.linspace(0.8, 3.0, 12)
    all_results['H2'] = run_pes_scan(
        'H2 (symmetric)', Z_A=1.0, Z_B=1.0,
        R_grid=R_grid_h2, l_max_values=l_max_values,
        n_alpha=n_alpha, n_Re=n_Re,
    )

    # ---- System 2: HeH+ (2:1 asymmetric) ----
    R_grid_heh = np.linspace(1.0, 4.0, 12)
    all_results['HeH+'] = run_pes_scan(
        'HeH+ (2:1 asymmetric)', Z_A=2.0, Z_B=1.0,
        R_grid=R_grid_heh, l_max_values=l_max_values,
        origin='charge_center',  # Standard for heteronuclear
        n_alpha=n_alpha, n_Re=n_Re,
    )

    # ---- System 3: Effective LiH valence (Z_A=Z_B=1 at LiH distances) ----
    R_grid_lih = np.linspace(2.0, 7.0, 12)
    all_results['eff_LiH_midpoint'] = run_pes_scan(
        'Eff LiH valence (1:1, midpoint)', Z_A=1.0, Z_B=1.0,
        R_grid=R_grid_lih, l_max_values=l_max_values,
        origin='midpoint',
        n_alpha=n_alpha, n_Re=n_Re,
    )

    # ---- System 4: Effective LiH valence with charge-center origin ----
    # Simulate charge-center as if Z_A=3: z0 = R*(3-1)/(2*(3+1)) = R/4
    # But we use Z_A=Z_B=1 for the actual nuclear coupling.
    # The origin parameter with a float value gives explicit z0.
    # For charge-center with "effective Z_A=3, Z_B=1": z0 = R*(3-1)/(2*(3+1)) = R/4.
    # Since z0 depends on R, we run manually.
    all_results['eff_LiH_chargecenter'] = {}
    for l_max in l_max_values:
        print(f"\n{'='*60}")
        print(f"System: Eff LiH valence (1:1, charge-center origin), l_max={l_max}")
        print(f"{'='*60}")

        E_total = np.zeros(len(R_grid_lih))
        t_start = time.time()

        for i, R in enumerate(R_grid_lih):
            z0 = R * (3.0 - 1.0) / (2.0 * (3.0 + 1.0))  # As if Z_A=3, Z_B=1
            try:
                result = solve_level4_h2_multichannel(
                    R=R,
                    l_max=l_max,
                    Z_A=1.0, Z_B=1.0,  # Bare charges still 1:1
                    n_alpha=n_alpha,
                    n_Re=n_Re,
                    origin=z0,  # Explicit origin shift
                    verbose=False,
                )
                E_total[i] = result['E_total']
                n_ch = result['n_ch']
            except Exception as e:
                print(f"  FAILED at R={R:.3f}: {e}")
                E_total[i] = np.nan

            print(f"  R={R:.3f}  z0={z0:.3f}  E_total={E_total[i]:.6f} Ha")

        t_elapsed = time.time() - t_start
        valid = ~np.isnan(E_total)
        R_eq, E_min, success = parabolic_fit_Req(R_grid_lih[valid], E_total[valid])

        print(f"\n  R_eq = {R_eq:.4f} bohr, E_min = {E_min:.6f} Ha")
        print(f"  Fit success: {success}")
        print(f"  Wall time: {t_elapsed:.1f}s")

        all_results['eff_LiH_chargecenter'][l_max] = {
            'R_grid': R_grid_lih.tolist(),
            'E_total': E_total.tolist(),
            'R_eq': float(R_eq),
            'E_min': float(E_min),
            'fit_success': success,
            'n_ch': n_ch,
            'wall_time_s': round(t_elapsed, 1),
        }

    # ---- Channel weights at l_max=4 for HeH+ and reference ----
    print(f"\n{'='*60}")
    print("Extracting channel weights at R_eq for l_max=4...")
    print(f"{'='*60}")

    channel_weights = {}

    # HeH+ at its R_eq
    heh_req = all_results['HeH+'][4]['R_eq']
    z0_heh = heh_req * (2.0 - 1.0) / (2.0 * (2.0 + 1.0))
    # Use R_e ~ 1.5 (typical equilibrium hyperradius for HeH+)
    channel_weights['HeH+'] = extract_channel_weights(
        R=heh_req, R_e=1.5, l_max=4, Z_A=2.0, Z_B=1.0, z0=z0_heh,
    )

    # H2 at its R_eq for comparison
    h2_req = all_results['H2'][4]['R_eq']
    channel_weights['H2'] = extract_channel_weights(
        R=h2_req, R_e=1.5, l_max=4, Z_A=1.0, Z_B=1.0,
    )

    all_results['channel_weights_lmax4'] = channel_weights

    # ---- Save data ----
    outpath = os.path.join(os.path.dirname(__file__), 'data',
                           'lmax_divergence_isolation.json')
    with open(outpath, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nData saved to {outpath}")

    # ---- Print summary table ----
    print(f"\n{'='*70}")
    print("SUMMARY: R_eq vs l_max for all systems")
    print(f"{'='*70}")
    print(f"{'System':<30} {'l_max':>5} {'R_eq':>8} {'n_ch':>5} {'time':>7}")
    print(f"{'-'*30} {'-'*5} {'-'*8} {'-'*5} {'-'*7}")

    for sys_name in ['H2', 'HeH+', 'eff_LiH_midpoint', 'eff_LiH_chargecenter']:
        for l_max in l_max_values:
            d = all_results[sys_name][l_max]
            print(f"{sys_name:<30} {l_max:>5} {d['R_eq']:>8.4f} "
                  f"{d['n_ch']:>5} {d['wall_time_s']:>6.1f}s")
        print()

    # Drift rates
    print(f"\n{'='*70}")
    print("DRIFT RATES: (R_eq[l_max=4] - R_eq[l_max=2]) / 2")
    print(f"{'='*70}")
    for sys_name in ['H2', 'HeH+', 'eff_LiH_midpoint', 'eff_LiH_chargecenter']:
        r2 = all_results[sys_name][2]['R_eq']
        r4 = all_results[sys_name][4]['R_eq']
        drift = (r4 - r2) / 2.0
        print(f"  {sys_name:<30}: {drift:+.4f} bohr/l_max  "
              f"(R_eq: {r2:.4f} -> {r4:.4f})")

    # Composed LiH reference from CSV
    print(f"\n  {'Composed LiH (l_dependent)':<30}: +0.303 bohr/l_max  "
          f"(R_eq: 3.176 -> 3.783)")
    print(f"  {'Composed LiH (channel_blind)':<30}: +0.367 bohr/l_max  "
          f"(R_eq: 3.222 -> 3.956)")

    # Channel weights
    print(f"\n{'='*70}")
    print("CHANNEL WEIGHTS at l_max=4 (R_e=1.5)")
    print(f"{'='*70}")
    for sys_name in ['H2', 'HeH+']:
        weights = channel_weights[sys_name]
        print(f"\n  {sys_name} (R={all_results[sys_name][4]['R_eq']:.3f}):")
        sorted_w = sorted(weights.items(), key=lambda x: -x[1])
        for ch, w in sorted_w:
            if w > 0.001:
                print(f"    {ch}: {w:.4f}")


if __name__ == '__main__':
    main()
