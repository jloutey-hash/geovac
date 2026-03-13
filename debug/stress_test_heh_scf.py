"""
Phase 7: Heteronuclear SCF — HeH+ with Per-Atom Z_eff
======================================================
HeH+ doesn't bind with frozen HeH2+ orbitals because:
  - Orbital concentrated on He (higher Z)
  - J_11 ~ 1.3 Ha — huge electron repulsion
  - No flexibility to delocalize toward H

Fix: optimize independent Z_eff_A (He), Z_eff_B (H) to allow
orbital shape change. Expected:
  - Z_He^eff < 2.0 (screening)
  - Z_H^eff > 1.0 (partial delocalization)
  - D_e > 0 (binding)
  - R_eq ~ 1.5 bohr (exact 1.46)

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

from geovac.prolate_heteronuclear_scf import (
    heteronuclear_scf_energy,
    optimize_heteronuclear_zeff,
    heteronuclear_scf_ci,
)
from geovac.prolate_spheroidal_lattice import fit_spectroscopic_constants


def main():
    print("=" * 70)
    print("  HeH+ HETERONUCLEAR SCF WITH PER-ATOM Z_eff")
    print("  Z_A=2 (He), Z_B=1 (H)  |  Independent orbital exponents")
    print("  Date: 2026-03-13")
    print("=" * 70)

    # Reference values
    E_He = -2.9037   # He atom exact
    E_Hplus = 0.0     # H+ energy
    E_ref = E_He + E_Hplus  # Dissociation limit
    R_eq_exact = 1.46  # bohr
    D_e_exact = 0.075  # Ha

    # --- Test 1: Z_eff scan at R=1.5 ---
    print("\n" + "=" * 60)
    print("  TEST 1: Z_eff landscape at R=1.5")
    print("=" * 60)

    R = 1.5
    print(f"\n  {'Z_eff_A':>8s} {'Z_eff_B':>8s} {'E_HF':>12s} {'J_gg':>10s}")
    print(f"  {'-' * 42}")

    # Scan a few points to show the landscape
    for za in [1.5, 1.7, 2.0]:
        for zb in [0.8, 1.0, 1.2, 1.5]:
            res = heteronuclear_scf_energy(R, 2.0, 1.0, za, zb,
                                           N_xi_solve=3000, N_grid=30)
            print(f"  {za:8.3f} {zb:8.3f} {res['E_HF']:12.6f} "
                  f"{res.get('J_gg', np.nan):10.6f}")

    # --- Test 2: Full optimization at R=1.5 ---
    print("\n" + "=" * 60)
    print("  TEST 2: Per-atom Z_eff optimization at R=1.5")
    print("=" * 60)
    opt = optimize_heteronuclear_zeff(R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
                                     N_xi_solve=5000, N_grid=40, verbose=True)

    if not np.isnan(opt['E_HF_opt']):
        De_opt = E_ref - opt['E_HF_opt']
        De_froz = E_ref - opt['E_HF_frozen']
        print(f"\n  D_e(opt)    = {De_opt:.6f} Ha")
        print(f"  D_e(frozen) = {De_froz:.6f} Ha")
        print(f"  Bound (opt): {De_opt > 0}")
        print(f"  Bound (frozen): {De_froz > 0}")

    # --- Test 3: SCF + CI at R=1.5 ---
    print("\n" + "=" * 60)
    print("  TEST 3: Heteronuclear SCF + CI at R=1.5")
    print("=" * 60)
    res = heteronuclear_scf_ci(R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
                               N_xi_solve=5000, N_grid=40, verbose=True)

    if not np.isnan(res['E_total']):
        De = E_ref - res['E_total']
        print(f"\n  D_e(SCF+CI) = {De:.6f} Ha ({De / D_e_exact * 100:.1f}% of exact)")

    # --- Test 4: PES scan ---
    print("\n" + "=" * 60)
    print("  TEST 4: HeH+ SCF PES scan")
    print("=" * 60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])
    results = []

    for R in R_vals:
        print(f"\n  --- R = {R:.2f} ---")
        try:
            r = heteronuclear_scf_ci(R=R, Z_A_phys=2.0, Z_B_phys=1.0,
                                     N_xi_solve=5000, N_grid=35,
                                     verbose=False)
            results.append(r)
            De = E_ref - r['E_total']
            print(f"  Z_eff=({r.get('Z_eff_A', np.nan):.3f}, {r.get('Z_eff_B', np.nan):.3f})  "
                  f"E_CI={r['E_total']:.6f}  D_e={De:.6f}  [{r['time']:.0f}s]")
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append({'R': R, 'E_total': np.nan, 'E_HF_scf': np.nan})

    # Summary table
    print(f"\n  {'R':>6s} {'E_CI':>12s} {'D_e':>10s} {'Z_eff_A':>8s} {'Z_eff_B':>8s}")
    print(f"  {'-' * 50}")
    for r in results:
        R = r['R']
        E = r.get('E_total', np.nan)
        De = E_ref - E
        za = r.get('Z_eff_A', np.nan)
        zb = r.get('Z_eff_B', np.nan)
        print(f"  {R:6.3f} {E:12.6f} {De:10.6f} {za:8.4f} {zb:8.4f}")

    # Fit spectroscopic constants
    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r.get('E_total', np.nan) for r in results])
    valid = ~np.isnan(E_arr)

    if np.sum(valid) >= 5:
        fit = fit_spectroscopic_constants(R_arr[valid], E_arr[valid])
        D_e_fit = E_ref - fit['E_min']

        R_err = abs(fit['R_eq'] - R_eq_exact) / R_eq_exact * 100

        print(f"\n  Spectroscopic constants:")
        print(f"    R_eq = {fit['R_eq']:.4f} bohr (exact {R_eq_exact}, err {R_err:.1f}%)")
        print(f"    D_e  = {D_e_fit:.6f} Ha (exact {D_e_exact})")

        print(f"\n  {'=' * 55}")
        if D_e_fit > 0:
            print(f"  HeH+ IS BOUND (D_e = {D_e_fit:.4f} Ha, "
                  f"{D_e_fit / D_e_exact * 100:.0f}% of exact)")
        else:
            print(f"  HeH+ NOT BOUND (D_e = {D_e_fit:.4f} Ha)")
        print(f"  {'=' * 55}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_heh_scf.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# HeH+ SCF with Per-Atom Z_eff Optimization\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Z_A_phys=2 (He), Z_B_phys=1 (H)\n")
        f.write(f"# Dissociation limit: He + H+ = {E_ref:.6f} Ha\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_CI':>14s} {'E_HF':>14s} "
                f"{'D_e':>14s} {'Z_eff_A':>10s} {'Z_eff_B':>10s}\n")
        for r in results:
            R_val = r['R']
            E = r.get('E_total', np.nan)
            Ehf = r.get('E_HF_scf', np.nan)
            za = r.get('Z_eff_A', np.nan)
            zb = r.get('Z_eff_B', np.nan)
            f.write(f"  {R_val:8.4f} {E:14.8f} {Ehf:14.8f} "
                    f"{E_ref - E:14.8f} {za:10.4f} {zb:10.4f}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
