#!/usr/bin/env python3
"""
Fast l_max divergence isolation — abbreviated runs for remaining systems.

Uses n_alpha=100, 8 R points, and skips l_max=3 for symmetric systems
where only drift rate matters. Incorporates full-resolution data from
the initial run for H2 and HeH+ l_max=2,3.
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


def parabolic_fit_Req(R_grid, E_grid):
    idx_min = np.argmin(E_grid)
    if idx_min == 0 or idx_min == len(R_grid) - 1:
        return R_grid[idx_min], E_grid[idx_min], False
    lo = max(0, idx_min - 2)
    hi = min(len(R_grid), idx_min + 3)
    R_fit = R_grid[lo:hi]
    E_fit = E_grid[lo:hi]
    coeffs = np.polyfit(R_fit, E_fit, 2)
    a, b, c = coeffs
    if a <= 0:
        return R_grid[idx_min], E_grid[idx_min], False
    R_eq = -b / (2 * a)
    E_min = c - b**2 / (4 * a)
    if R_eq < R_grid[0] * 0.8 or R_eq > R_grid[-1] * 1.2:
        return R_grid[idx_min], E_grid[idx_min], False
    return R_eq, E_min, True


def run_scan(label, Z_A, Z_B, R_grid, l_max, origin='midpoint', n_alpha=100, n_Re=300):
    print(f"\n{'='*60}")
    print(f"{label}, l_max={l_max}, n_ch=?, origin={origin}")
    print(f"{'='*60}")

    E_total = np.zeros(len(R_grid))
    t0 = time.time()
    n_ch = 0

    for i, R in enumerate(R_grid):
        if isinstance(origin, str) and origin == 'charge_center':
            orig = origin
        elif origin == 'lih_shift':
            # Charge-center origin as if Z_A=3, Z_B=1 but actual charges 1,1
            orig = float(R * (3.0 - 1.0) / (2.0 * (3.0 + 1.0)))
        else:
            orig = origin

        try:
            result = solve_level4_h2_multichannel(
                R=R, l_max=l_max, Z_A=Z_A, Z_B=Z_B,
                n_alpha=n_alpha, n_Re=n_Re, origin=orig, verbose=False,
            )
            E_total[i] = result['E_total']
            n_ch = result['n_ch']
        except Exception as e:
            print(f"  FAILED at R={R:.3f}: {e}")
            E_total[i] = np.nan

        print(f"  R={R:.3f}  E={E_total[i]:.6f} Ha")

    dt = time.time() - t0
    valid = ~np.isnan(E_total)
    R_eq, E_min, ok = parabolic_fit_Req(R_grid[valid], E_total[valid])
    print(f"  => R_eq={R_eq:.4f}, E_min={E_min:.6f}, n_ch={n_ch}, {dt:.0f}s")

    return {
        'R_eq': float(R_eq), 'E_min': float(E_min),
        'n_ch': n_ch, 'wall_time_s': round(dt, 1),
        'fit_success': ok,
        'R_grid': R_grid.tolist(), 'E_total': E_total.tolist(),
    }


def extract_channel_weights(R, R_e, l_max, Z_A, Z_B, z0=0.0, n_alpha=100):
    rho = R / (2.0 * R_e)
    evals, vecs, alpha, channels = solve_angular_multichannel(
        rho, R_e, l_max=l_max, n_alpha=n_alpha, Z_A=Z_A, Z_B=Z_B, z0=z0,
    )
    vec0 = vecs[0]
    n_ch = len(channels)
    weights = {}
    for c_idx, ch in enumerate(channels):
        block = vec0[c_idx * n_alpha : (c_idx + 1) * n_alpha]
        w = np.sum(block**2) * (np.pi / 2) / (n_alpha + 1)
        weights[str(ch)] = float(w)
    total = sum(weights.values())
    for k in weights:
        weights[k] /= total
    return weights


def main():
    results = {}

    # ===== Data from full-resolution run (n_alpha=200, 12 R points) =====
    results['H2'] = {
        2: {'R_eq': 1.4514, 'E_min': -1.158317, 'n_ch': 5, 'wall_time_s': 470.1, 'fit_success': True},
        3: {'R_eq': 1.4582, 'E_min': -1.159687, 'n_ch': 9, 'wall_time_s': 1214.6, 'fit_success': True},
        4: {'R_eq': 1.4912, 'E_min': -1.174004, 'n_ch': 14, 'wall_time_s': 3796.1, 'fit_success': True},
    }
    results['HeH+'] = {
        2: {'R_eq': 1.4640, 'E_min': -2.865305, 'n_ch': 9, 'wall_time_s': 1375.4, 'fit_success': True},
        3: {'R_eq': 1.6513, 'E_min': -2.926623, 'n_ch': 16, 'wall_time_s': 6600.4, 'fit_success': True},
    }

    # ===== HeH+ l_max=4: run with 8 points, n_alpha=100 =====
    R_heh = np.linspace(1.0, 4.0, 8)
    results['HeH+'][4] = run_scan(
        'HeH+ (2:1)', Z_A=2.0, Z_B=1.0, R_grid=R_heh, l_max=4,
        origin='charge_center', n_alpha=100, n_Re=300,
    )

    # ===== System 3: Effective LiH valence (1:1, midpoint, LiH distances) =====
    R_lih = np.linspace(2.0, 7.0, 8)
    results['eff_LiH_midpoint'] = {}
    for lm in [2, 4]:
        results['eff_LiH_midpoint'][lm] = run_scan(
            'Eff LiH (1:1, midpoint)', Z_A=1.0, Z_B=1.0,
            R_grid=R_lih, l_max=lm, origin='midpoint',
        )

    # ===== System 4: Effective LiH valence with charge-center origin =====
    results['eff_LiH_chargecenter'] = {}
    for lm in [2, 4]:
        results['eff_LiH_chargecenter'][lm] = run_scan(
            'Eff LiH (1:1, CC origin)', Z_A=1.0, Z_B=1.0,
            R_grid=R_lih, l_max=lm, origin='lih_shift',
        )

    # ===== Channel weights =====
    print(f"\n{'='*60}")
    print("Channel weights at l_max=4")
    print(f"{'='*60}")

    channel_weights = {}

    # HeH+ at R_eq
    heh_req = results['HeH+'][4]['R_eq']
    z0_heh = heh_req * (2.0 - 1.0) / (2.0 * (2.0 + 1.0))
    channel_weights['HeH+'] = extract_channel_weights(
        R=heh_req, R_e=1.5, l_max=4, Z_A=2.0, Z_B=1.0, z0=z0_heh,
    )

    # H2 at R_eq
    channel_weights['H2'] = extract_channel_weights(
        R=1.4912, R_e=1.5, l_max=4, Z_A=1.0, Z_B=1.0,
    )

    results['channel_weights_lmax4'] = channel_weights

    # ===== Save data =====
    outpath = os.path.join(os.path.dirname(__file__), 'data',
                           'lmax_divergence_isolation.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nData saved to {outpath}")

    # ===== Summary table =====
    print(f"\n{'='*70}")
    print("SUMMARY: R_eq vs l_max")
    print(f"{'='*70}")
    print(f"{'System':<30} {'l_max':>5} {'R_eq':>8} {'n_ch':>5}")
    print(f"{'-'*30} {'-'*5} {'-'*8} {'-'*5}")

    for sys_name in ['H2', 'HeH+', 'eff_LiH_midpoint', 'eff_LiH_chargecenter']:
        for lm in sorted(results[sys_name].keys()):
            if isinstance(lm, int):
                d = results[sys_name][lm]
                print(f"{sys_name:<30} {lm:>5} {d['R_eq']:>8.4f} {d['n_ch']:>5}")
        print()

    # Drift rates
    print(f"\n{'='*70}")
    print("DRIFT RATES")
    print(f"{'='*70}")

    # H2: use l2->l4
    r2 = results['H2'][2]['R_eq']
    r4 = results['H2'][4]['R_eq']
    print(f"  H2 (1:1 sym):              {(r4-r2)/2:+.4f} bohr/l_max  ({r2:.4f} -> {r4:.4f})")

    # HeH+: use l2->l4
    r2 = results['HeH+'][2]['R_eq']
    r4 = results['HeH+'][4]['R_eq']
    print(f"  HeH+ (2:1 asym):           {(r4-r2)/2:+.4f} bohr/l_max  ({r2:.4f} -> {r4:.4f})")

    # Eff LiH midpoint
    r2 = results['eff_LiH_midpoint'][2]['R_eq']
    r4 = results['eff_LiH_midpoint'][4]['R_eq']
    print(f"  Eff LiH (1:1, midpoint):   {(r4-r2)/2:+.4f} bohr/l_max  ({r2:.4f} -> {r4:.4f})")

    # Eff LiH charge-center
    r2 = results['eff_LiH_chargecenter'][2]['R_eq']
    r4 = results['eff_LiH_chargecenter'][4]['R_eq']
    print(f"  Eff LiH (1:1, CC origin):  {(r4-r2)/2:+.4f} bohr/l_max  ({r2:.4f} -> {r4:.4f})")

    print(f"\n  --- Composed LiH reference (from CSV) ---")
    print(f"  Composed (l_dependent PK): +0.303 bohr/l_max  (3.176 -> 3.783)")
    print(f"  Composed (channel_blind):  +0.367 bohr/l_max  (3.222 -> 3.956)")

    # Channel weights
    print(f"\n{'='*70}")
    print("CHANNEL WEIGHTS at l_max=4 (R_e=1.5)")
    print(f"{'='*70}")
    for sys_name in ['H2', 'HeH+']:
        weights = channel_weights[sys_name]
        print(f"\n  {sys_name}:")
        sorted_w = sorted(weights.items(), key=lambda x: -x[1])
        for ch, w in sorted_w:
            if w > 0.001:
                print(f"    {ch}: {w:.4f}")


if __name__ == '__main__':
    main()
