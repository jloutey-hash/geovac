"""
Phase 8: Full Grid-Based SCF for H2
=====================================
Instead of optimizing a single Z_eff, solve the Fock equation directly
on the 2D (xi, eta) grid:

  F phi = eps phi
  F = T + V_nuc + V_J/2   (closed-shell, one spatial orbital)

where V_J is the Coulomb potential from rho = 2|phi|^2, computed via
azimuthal averaging on the prolate spheroidal grid.

This captures orbital polarization and shape flexibility beyond
what uniform Z_eff scaling provides.

Expected improvement over Eckart: ~5-15% additional D_e recovery.

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
    grid_scf,
    eckart_scf_energy,
    optimize_zeff_eckart,
    fit_spectroscopic_constants,
)


def main():
    print("=" * 70)
    print("  H2 FULL GRID-BASED SCF")
    print("  Fock operator on 2D prolate spheroidal grid")
    print("  Date: 2026-03-13")
    print("=" * 70)

    E_H = -0.5
    E_H2_exact = -1.1745
    R_eq_exact = 1.401
    D_e_exact = 0.1745

    # --- Test 1: Grid convergence at R=1.4 ---
    print("\n" + "=" * 60)
    print("  TEST 1: Grid convergence at R=1.4")
    print("=" * 60)

    for N_fd in [20, 30, 40]:
        print(f"\n  --- N_fd = {N_fd} ---")
        res = grid_scf(R=1.4, N_xi_fd=N_fd, N_eta_fd=N_fd,
                       xi_max_fd=15.0, verbose=True)
        if not np.isnan(res['E_total']):
            De = 2 * E_H - res['E_total']
            print(f"  D_e = {De:.6f} Ha ({De / D_e_exact * 100:.1f}% of exact)")

    # --- Test 2: Comparison with Eckart at R=1.4 ---
    print("\n" + "=" * 60)
    print("  TEST 2: Grid SCF vs Eckart at R=1.4")
    print("=" * 60)

    res_grid = grid_scf(R=1.4, N_xi_fd=40, N_eta_fd=40,
                        xi_max_fd=15.0, verbose=True)
    opt_eckart = optimize_zeff_eckart(R=1.4, N_xi_solve=5000, N_grid=50)

    if not np.isnan(res_grid['E_total']):
        De_grid = 2 * E_H - res_grid['E_total']
        De_eckart = 2 * E_H - opt_eckart['E_HF_opt']
        De_frozen = 2 * E_H - opt_eckart['E_HF_frozen']

        print(f"\n  {'Method':<20s} {'E_HF':>12s} {'D_e':>10s} {'%exact':>8s}")
        print(f"  {'-' * 52}")
        print(f"  {'Frozen (Z=1)':<20s} {opt_eckart['E_HF_frozen']:12.6f} "
              f"{De_frozen:10.6f} {De_frozen / D_e_exact * 100:7.1f}%")
        print(f"  {'Eckart (Z_eff*)':<20s} {opt_eckart['E_HF_opt']:12.6f} "
              f"{De_eckart:10.6f} {De_eckart / D_e_exact * 100:7.1f}%")
        print(f"  {'Grid SCF':<20s} {res_grid['E_total']:12.6f} "
              f"{De_grid:10.6f} {De_grid / D_e_exact * 100:7.1f}%")
        print(f"  {'Exact':<20s} {E_H2_exact:12.6f} "
              f"{D_e_exact:10.6f} {'100.0':>7s}%")

        improvement = De_grid - De_eckart
        print(f"\n  Grid SCF improvement over Eckart: {improvement:.6f} Ha "
              f"({improvement / D_e_exact * 100:.1f}% of exact)")

    # --- Test 3: PES scan ---
    print("\n" + "=" * 60)
    print("  TEST 3: Grid SCF PES scan")
    print("=" * 60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0])
    results = []

    for R in R_vals:
        print(f"\n  --- R = {R:.2f} ---")
        try:
            r = grid_scf(R=R, N_xi_fd=35, N_eta_fd=35,
                         xi_max_fd=max(15.0, R * 4),
                         max_iter=40, verbose=False)
            results.append(r)
            De = 2 * E_H - r['E_total']
            print(f"  E_HF={r['E_total']:.6f}  D_e={De:.6f}  "
                  f"J_gg={r['J_gg']:.4f}  iters={r['n_iter']}  "
                  f"conv={'Y' if r['converged'] else 'N'}  [{r['time']:.0f}s]")
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append({'R': R, 'E_total': np.nan, 'time': 0})

    # Summary table
    print(f"\n  {'R':>6s} {'E_HF':>12s} {'D_e':>10s} {'J_gg':>10s} "
          f"{'iters':>6s} {'conv':>5s}")
    print(f"  {'-' * 55}")
    for r in results:
        R_val = r['R']
        E = r.get('E_total', np.nan)
        De = 2 * E_H - E
        J = r.get('J_gg', np.nan)
        it = r.get('n_iter', 0)
        cv = 'Y' if r.get('converged', False) else 'N'
        print(f"  {R_val:6.3f} {E:12.6f} {De:10.6f} {J:10.6f} {it:6d} {cv:>5s}")

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
        print(f"  GRID SCF SUMMARY:")
        print(f"    D_e: {D_e_fit:.4f} Ha ({D_e_fit / D_e_exact * 100:.0f}% of exact)")
        print(f"    vs Eckart HF: {0.101:.3f} Ha (58%)")
        print(f"    vs frozen CI: {0.072:.3f} Ha (41%)")
        print(f"  {'=' * 55}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_h2_grid_scf.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# H2 Full Grid-Based SCF on Prolate Spheroidal Lattice\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Method: Fock operator on 2D (xi,eta) FD grid\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_HF':>14s} {'D_e':>14s} "
                f"{'J_gg':>14s} {'iters':>6s}\n")
        for r in results:
            R_val = r['R']
            E = r.get('E_total', np.nan)
            J = r.get('J_gg', np.nan)
            it = r.get('n_iter', 0)
            f.write(f"  {R_val:8.4f} {E:14.8f} {2 * E_H - E:14.8f} "
                    f"{J:14.8f} {it:6d}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
