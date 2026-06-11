"""
Phase 6: CI with Relaxed Orbitals
==================================
Combine Eckart SCF orbital improvement with CI correlation:
  1. Optimize Z_eff via Eckart variational
  2. Generate sigma_g and sigma_u orbitals at Z_eff*
  3. Run 2x2 CI with these improved orbitals

Expected: D_e should improve over both frozen CI (0.072 Ha) and
Eckart SCF alone (0.101 Ha). Target: >65% of exact D_e (0.175 Ha).

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import os
if sys.stdout.encoding and sys.stdout.encoding.lower().startswith('cp'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_scf import (
    relaxed_orbital_ci,
    optimize_zeff_eckart,
    fit_spectroscopic_constants,
)


def main():
    print("=" * 70)
    print("  H2 RELAXED-ORBITAL CI")
    print("  Eckart SCF + 2x2 CI on prolate spheroidal lattice")
    print("  Date: 2026-03-13")
    print("=" * 70)

    E_H = -0.5
    E_H2_exact = -1.1745
    R_eq_exact = 1.401
    D_e_exact = 0.1745

    # --- Test 1: Single point at R=1.4 ---
    print("\n" + "=" * 60)
    print("  TEST 1: Relaxed CI at R=1.4")
    print("=" * 60)
    res = relaxed_orbital_ci(R=1.4, N_xi_solve=5000, N_grid=50, verbose=True)

    if not np.isnan(res['E_total']):
        D_e = 2 * E_H - res['E_total']
        print(f"\n  D_e(relaxed CI) = {D_e:.6f} Ha ({D_e / D_e_exact * 100:.1f}% of exact)")
        print(f"  D_e(Eckart HF)  = {2 * E_H - res['E_HF_scf']:.6f} Ha")
        print(f"  D_e(frozen HF)  = {2 * E_H - res['E_HF_frozen']:.6f} Ha")
        print(f"  CI correlation:  {res['E_total'] - res['E_HF_scf']:.6f} Ha")

    # --- Test 2: Comparison table ---
    print("\n" + "=" * 60)
    print("  TEST 2: Method comparison at R=1.4")
    print("=" * 60)
    from debug.stress_test_prolate_h2 import h2_minimal_ci
    frozen = h2_minimal_ci(R=1.4, N_xi_solve=5000, N_grid=50, verbose=False)

    print(f"\n  {'Method':<25s} {'E_total':>10s} {'D_e':>10s} {'%exact':>8s}")
    print(f"  {'-' * 55}")

    methods = [
        ("Frozen CI (Z=1)", frozen['E_total'], 2 * E_H - frozen['E_total']),
        ("Eckart HF", res['E_HF_scf'], 2 * E_H - res['E_HF_scf']),
        ("Relaxed CI (Z_eff*)", res['E_total'], 2 * E_H - res['E_total']),
        ("Exact", E_H2_exact, D_e_exact),
    ]
    for name, E, De in methods:
        pct = De / D_e_exact * 100
        print(f"  {name:<25s} {E:10.6f} {De:10.6f} {pct:7.1f}%")

    # --- Test 3: PES scan ---
    print("\n" + "=" * 60)
    print("  TEST 3: Relaxed CI PES scan")
    print("=" * 60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])
    results = []
    for R in R_vals:
        print(f"\n  --- R = {R:.2f} ---")
        try:
            r = relaxed_orbital_ci(R=R, N_xi_solve=5000, N_grid=40, verbose=False)
            results.append(r)
            De = 2 * E_H - r['E_total']
            print(f"  Z_eff*={r.get('Z_eff_opt', np.nan):.3f}  "
                  f"E_CI={r['E_total']:.6f}  D_e={De:.6f}  [{r['time']:.0f}s]")
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append({'R': R, 'E_total': np.nan, 'E_HF_scf': np.nan,
                            'E_HF_frozen': np.nan})

    # Summary table
    print(f"\n  {'R':>6s} {'E_CI':>12s} {'D_e':>10s} {'Z_eff*':>8s} {'c_u^2':>8s}")
    print(f"  {'-' * 50}")
    for r in results:
        R = r['R']
        E = r.get('E_total', np.nan)
        De = 2 * E_H - E
        z = r.get('Z_eff_opt', np.nan)
        cu = r.get('c_u', np.nan)
        cu2 = cu ** 2 if not np.isnan(cu) else np.nan
        print(f"  {R:6.3f} {E:12.6f} {De:10.6f} {z:8.4f} {cu2:8.4f}")

    # Fit spectroscopic constants
    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r.get('E_total', np.nan) for r in results])
    valid = ~np.isnan(E_arr)

    if np.sum(valid) >= 5:
        fit = fit_spectroscopic_constants(R_arr[valid], E_arr[valid])
        D_e_fit = 2 * E_H - fit['E_min']

        R_err = abs(fit['R_eq'] - R_eq_exact) / R_eq_exact * 100
        D_err = abs(D_e_fit - D_e_exact) / D_e_exact * 100

        print(f"\n  Spectroscopic constants:")
        print(f"    R_eq = {fit['R_eq']:.4f} bohr (exact {R_eq_exact}, err {R_err:.1f}%)")
        print(f"    D_e  = {D_e_fit:.6f} Ha (exact {D_e_exact}, err {D_err:.1f}%)")

        print(f"\n  {'=' * 55}")
        print(f"  RELAXED CI RESULTS:")
        print(f"    D_e: {D_e_fit:.4f} Ha ({D_e_fit / D_e_exact * 100:.0f}% of exact)")
        print(f"    vs frozen CI:  0.072 Ha (41%)")
        print(f"    vs Eckart HF:  0.101 Ha (58%)")
        print(f"  {'=' * 55}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_h2_relaxed_ci.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# H2 Relaxed-Orbital CI on Prolate Spheroidal Lattice\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Method: Eckart Z_eff optimization + 2x2 CI\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_CI':>14s} {'E_HF_SCF':>14s} "
                f"{'D_e':>14s} {'Z_eff':>10s}\n")
        for r in results:
            R = r['R']
            E = r.get('E_total', np.nan)
            Ehf = r.get('E_HF_scf', np.nan)
            z = r.get('Z_eff_opt', np.nan)
            f.write(f"  {R:8.4f} {E:14.8f} {Ehf:14.8f} "
                    f"{2 * E_H - E:14.8f} {z:10.4f}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
