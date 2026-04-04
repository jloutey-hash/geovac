"""Track AO Part 2: Coupled-channel convergence scan for 4-electron LiH.

Tests 1, 3, 6, 12 coupled channels at l_max=2 with different q_modes.
"""

import numpy as np
import sys
import os
import time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.n_electron_solver import solve_4e_lih_coupled


def run_convergence_scan():
    """Run coupled-channel convergence at l_max=2."""

    R_values = np.array([0.7, 1.0, 1.5, 2.0, 3.0])
    n_channels_list = [1, 3, 6]
    q_mode = 'diagonal'

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact

    print(f"{'='*70}")
    print(f"COUPLED-CHANNEL CONVERGENCE SCAN")
    print(f"l_max=2, q_mode={q_mode}")
    print(f"R values: {R_values}")
    print(f"{'='*70}\n")

    all_results = {}

    for n_ch in n_channels_list:
        print(f"\n{'#'*60}")
        print(f"  n_channels = {n_ch}")
        print(f"{'#'*60}")
        t0 = time.time()

        results_ch = {}
        for R in R_values:
            result = solve_4e_lih_coupled(
                R, n_grid=6, l_max=2, n_channels=n_ch,
                q_mode=q_mode, verbose=True,
            )
            results_ch[R] = result

        t_total = time.time() - t0
        all_results[n_ch] = results_ch

        print(f"\n  {n_ch}-channel scan complete in {t_total:.1f}s")

    # Summary table
    print(f"\n{'='*70}")
    print(f"CONVERGENCE SUMMARY TABLE")
    print(f"{'='*70}")
    print(f"{'n_ch':>5s}  {'R':>6s}  {'E_coupled':>12s}  {'E_adiab':>12s}  "
          f"{'Delta':>10s}  {'D_e_coup':>10s}  {'D_e_adiab':>10s}")
    print(f"{'':>5s}  {'bohr':>6s}  {'Ha':>12s}  {'Ha':>12s}  "
          f"{'Ha':>10s}  {'Ha':>10s}  {'Ha':>10s}")
    print(f"{'-'*5}  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*10}")

    for n_ch in n_channels_list:
        for R in R_values:
            r = all_results[n_ch][R]
            delta = r['E_total'] - r['E_total_single']
            print(f"{n_ch:5d}  {R:6.3f}  {r['E_total']:12.6f}  "
                  f"{r['E_total_single']:12.6f}  {delta:+10.6f}  "
                  f"{r['D_e']:10.6f}  {r['D_e_single']:10.6f}")
        print()

    # Min energy / equilibrium analysis
    print(f"\n{'='*70}")
    print(f"EQUILIBRIUM ANALYSIS")
    print(f"{'='*70}")
    for n_ch in n_channels_list:
        E_vals = np.array([all_results[n_ch][R]['E_total'] for R in R_values])
        E_single_vals = np.array([all_results[n_ch][R]['E_total_single'] for R in R_values])
        D_vals = np.array([all_results[n_ch][R]['D_e'] for R in R_values])
        D_single_vals = np.array([all_results[n_ch][R]['D_e_single'] for R in R_values])

        i_min = np.argmin(E_vals)
        i_min_s = np.argmin(E_single_vals)

        print(f"\n  {n_ch} channels:")
        print(f"    Coupled:   R_eq ~ {R_values[i_min]:.3f} bohr, "
              f"E_min = {E_vals[i_min]:.6f} Ha, "
              f"D_e = {D_vals[i_min]:.6f} Ha ({D_vals[i_min]/D_e_exact*100:.1f}%)")
        print(f"    Adiabatic: R_eq ~ {R_values[i_min_s]:.3f} bohr, "
              f"E_min = {E_single_vals[i_min_s]:.6f} Ha, "
              f"D_e = {D_single_vals[i_min_s]:.6f} Ha ({D_single_vals[i_min_s]/D_e_exact*100:.1f}%)")

    # q_mode comparison at 3 channels
    print(f"\n\n{'='*70}")
    print(f"Q-MODE COMPARISON (3 channels, R=1.0)")
    print(f"{'='*70}")
    R_test = 1.0
    for qm in ['none', 'diagonal', 'full']:
        result = solve_4e_lih_coupled(
            R_test, n_grid=6, l_max=2, n_channels=3,
            q_mode=qm, verbose=False,
        )
        delta = result['E_total'] - result['E_total_single']
        print(f"  q_mode={qm:10s}  E_total = {result['E_total']:.6f}  "
              f"Delta = {delta:+.6f}  D_e = {result['D_e']:.6f} "
              f"({result['D_e_pct']:.1f}%)")

    # Save results
    output_dir = os.path.dirname(__file__)
    print(f"\nResults saved to {output_dir}/")


if __name__ == '__main__':
    run_convergence_scan()
