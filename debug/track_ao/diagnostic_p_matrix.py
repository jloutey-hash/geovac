"""Track AO Part 1: P-matrix diagnostic for 4-electron coupled-channel.

Computes P-matrix at l_max=2, measures magnitudes, compares to Level 3.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.n_electron_solver import (
    compute_multichannel_angular_sweep_4e,
    compute_p_matrix_fd_4e,
    CENTRIFUGAL_4E,
)


def run_p_diagnostic():
    """Run P-matrix diagnostic at l_max=2."""
    R = 1.0  # Near equilibrium
    n_channels = 6  # Up to 6 channels for diagnostic

    # Dense R_e grid for good FD accuracy
    R_e_grid = np.concatenate([
        np.linspace(0.5, 1.5, 15),
        np.linspace(1.5, 4.0, 20),
        np.linspace(4.0, 10.0, 10),
    ])
    R_e_grid = np.unique(R_e_grid)

    print(f"P-matrix diagnostic at R={R}, l_max=2, {n_channels} channels")
    print(f"R_e grid: {len(R_e_grid)} points from {R_e_grid[0]:.2f} to {R_e_grid[-1]:.2f}")
    print()

    angular_data = compute_multichannel_angular_sweep_4e(
        R=R, R_e_grid=R_e_grid,
        Z_A=3.0, Z_B=1.0, z0=0.0,
        n_grid=6, l_max=2, n_channels=n_channels,
        symmetry='parity', verbose=True,
    )

    p_data = compute_p_matrix_fd_4e(angular_data, n_channels=n_channels)

    P = p_data['P']
    DBOC = p_data['DBOC']

    print(f"\n{'='*70}")
    print(f"P-MATRIX DIAGNOSTIC TABLE")
    print(f"{'='*70}")
    print(f"{'R_e':>8s}  {'||P||_F':>10s}  {'max|P_off|':>10s}  "
          f"{'DBOC_0':>10s}  {'DBOC_1':>10s}  {'DBOC_2':>10s}")
    print(f"{'bohr':>8s}  {'':>10s}  {'':>10s}  "
          f"{'Ha':>10s}  {'Ha':>10s}  {'Ha':>10s}")
    print(f"{'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
    for i in range(len(R_e_grid)):
        print(f"{R_e_grid[i]:8.3f}  {p_data['P_frob_norm'][i]:10.6f}  "
              f"{p_data['P_max_offdiag'][i]:10.6f}  "
              f"{DBOC[0, i]:10.6f}  {DBOC[1, i]:10.6f}  "
              f"{DBOC[2, i]:10.6f}")

    # Peak values
    print(f"\nPeak values:")
    print(f"  ||P||_F max:     {np.max(p_data['P_frob_norm']):.6f}")
    print(f"  |P_offdiag| max: {np.max(p_data['P_max_offdiag']):.6f}")
    for ch in range(min(n_channels, 6)):
        print(f"  DBOC_{ch} max:     {np.max(DBOC[ch]):.6f}")

    # Level 3 comparison
    print(f"\nLevel 3 comparison (He, Z=2, from Track B/D):")
    print(f"  Level 3 ||P||_F peak ~ 0.3-0.5 (at R ~ 1-2 bohr)")
    print(f"  Level 3 DBOC_0 peak  ~ 0.05-0.15")
    print(f"  Level 4N ||P||_F peak = {np.max(p_data['P_frob_norm']):.4f}")
    print(f"  Level 4N DBOC_0 peak  = {np.max(DBOC[0]):.6f}")

    # Eigenvalue gap analysis
    mu = angular_data['mu']
    print(f"\n{'='*70}")
    print(f"EIGENVALUE GAP ANALYSIS")
    print(f"{'='*70}")
    print(f"{'R_e':>8s}  {'mu_0':>10s}  {'mu_1':>10s}  {'mu_2':>10s}  "
          f"{'gap_01':>10s}  {'gap_12':>10s}")
    print(f"{'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
    for i in range(0, len(R_e_grid), 3):
        gap01 = mu[1, i] - mu[0, i] if n_channels > 1 else 0
        gap12 = mu[2, i] - mu[1, i] if n_channels > 2 else 0
        print(f"{R_e_grid[i]:8.3f}  {mu[0, i]:10.4f}  {mu[1, i]:10.4f}  "
              f"{mu[2, i]:10.4f}  {gap01:10.4f}  {gap12:10.4f}")

    # Save results
    output_dir = os.path.dirname(__file__)
    np.savez(os.path.join(output_dir, 'p_matrix_diagnostic.npz'),
             R_e_grid=R_e_grid,
             P=P, DBOC=DBOC,
             P_frob_norm=p_data['P_frob_norm'],
             P_max_offdiag=p_data['P_max_offdiag'],
             mu=mu, R=R)
    print(f"\nSaved to {output_dir}/p_matrix_diagnostic.npz")


if __name__ == '__main__':
    run_p_diagnostic()
