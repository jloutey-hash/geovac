"""Quick P-matrix + coupled-channel diagnostic for Track AO."""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.n_electron_solver import (
    compute_multichannel_angular_sweep_4e,
    compute_p_matrix_fd_4e,
    solve_4e_lih_coupled,
    CENTRIFUGAL_4E,
)
import time


def main():
    print("=" * 70)
    print("TRACK AO: COUPLED-CHANNEL DIAGNOSTIC")
    print("=" * 70)

    # PART 1: P-matrix diagnostic
    print("\n--- PART 1: P-MATRIX DIAGNOSTIC ---")
    R = 1.0
    n_channels = 3
    R_e_grid = np.concatenate([
        np.linspace(0.5, 1.5, 8),
        np.linspace(1.5, 4.0, 10),
        np.linspace(4.0, 8.0, 5),
    ])
    R_e_grid = np.unique(R_e_grid)

    t0 = time.time()
    angular_data = compute_multichannel_angular_sweep_4e(
        R=R, R_e_grid=R_e_grid,
        Z_A=3.0, Z_B=1.0, z0=0.0,
        n_grid=6, l_max=2, n_channels=n_channels,
        symmetry='s4', verbose=False,
    )
    p_data = compute_p_matrix_fd_4e(angular_data, n_channels=n_channels)
    t1 = time.time()
    print(f"  Angular sweep + P-matrix: {t1-t0:.1f}s")

    P = p_data['P']
    DBOC = p_data['DBOC']
    mu = angular_data['mu']

    print(f"\n  P-MATRIX DIAGNOSTIC TABLE (3 channels, l_max=2, R=1.0)")
    print(f"  {'R_e':>6s}  {'||P||_F':>8s}  {'max|P|':>8s}  "
          f"{'DBOC_0':>8s}  {'gap01':>8s}")
    for i in range(0, len(R_e_grid), 2):
        gap = mu[1, i] - mu[0, i] if n_channels > 1 else 0
        print(f"  {R_e_grid[i]:6.3f}  {p_data['P_frob_norm'][i]:8.4f}  "
              f"{p_data['P_max_offdiag'][i]:8.4f}  "
              f"{DBOC[0, i]:8.5f}  {gap:8.4f}")

    print(f"\n  Peak ||P||_F:  {np.max(p_data['P_frob_norm']):.6f}")
    print(f"  Peak |P_off|:  {np.max(p_data['P_max_offdiag']):.6f}")
    print(f"  Peak DBOC_0:   {np.max(DBOC[0]):.6f}")
    if n_channels > 1:
        print(f"  Peak DBOC_1:   {np.max(DBOC[1]):.6f}")

    # PART 2: Coupled-channel convergence
    print("\n\n--- PART 2: COUPLED-CHANNEL CONVERGENCE ---")
    print("Computing at R = 0.7, 1.0, 1.5, 2.0, 3.0 bohr")

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact

    R_values = [0.7, 1.0, 1.5, 2.0, 3.0]

    for n_ch in [1, 3, 6]:
        print(f"\n  --- {n_ch} channel(s), q_mode='diagonal' ---")
        t_start = time.time()

        for R_val in R_values:
            result = solve_4e_lih_coupled(
                R_val, n_grid=6, l_max=2, n_channels=n_ch,
                q_mode='diagonal', verbose=False,
            )
            delta = result['E_total'] - result['E_total_single']
            print(f"    R={R_val:.1f}: E_coup={result['E_total']:.6f} "
                  f"E_adiab={result['E_total_single']:.6f} "
                  f"dE={delta:+.6f} D_e={result['D_e']:.4f} "
                  f"({result['D_e_pct']:.0f}%)")

        print(f"    Time: {time.time()-t_start:.1f}s")

    # PART 3: q_mode comparison at 3 channels
    print(f"\n\n--- PART 3: Q-MODE COMPARISON (3 ch, R=1.0) ---")
    for qm in ['none', 'diagonal']:
        result = solve_4e_lih_coupled(
            1.0, n_grid=6, l_max=2, n_channels=3,
            q_mode=qm, verbose=False,
        )
        delta = result['E_total'] - result['E_total_single']
        print(f"  q_mode={qm:10s}: E={result['E_total']:.6f} "
              f"dE={delta:+.6f} D_e={result['D_e']:.4f} "
              f"({result['D_e_pct']:.0f}%)")

    print(f"\n  Exact D_e = {D_e_exact:.6f} Ha")
    print(f"\nDone.")


if __name__ == '__main__':
    main()
