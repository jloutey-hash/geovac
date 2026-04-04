"""l_max=3 PES scan for 4-electron LiH.

Uses reduced R_e grid (20 points instead of 48) to make each R-point ~20 min.
Targets key R values: 0.7, 1.0, 1.5, 2.0, 3.0.
"""
import sys
sys.path.insert(0, 'c:/Users/jlout/OneDrive/Desktop/Project_Geometric')

import numpy as np
import time
import os
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal

from geovac.n_electron_solver import (
    solve_angular_4e_multichannel,
    _enumerate_channels_4e,
    CENTRIFUGAL_4E,
)


def solve_4e_lih_fast(R, l_max=3, n_grid=4, n_Re_angular=20, n_Re_radial=200,
                       R_e_min=0.3, R_e_max=12.0, verbose=True):
    """Fast version of solve_4e_lih_multichannel with reduced R_e grid."""
    Z_A, Z_B = 3.0, 1.0
    z0 = 0.0  # midpoint

    R_A = R / 2.0
    R_B = R / 2.0

    # Reduced R_e grid: concentrate points near the expected minimum
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.5, 5),
        np.linspace(1.5, 4.0, 8),
        np.linspace(4.0, 8.0, 5),
        np.linspace(8.0, R_e_max, 2),
    ])
    R_e_angular = np.unique(R_e_angular)

    n_pts = len(R_e_angular)
    mu_vals = np.zeros(n_pts)

    t0 = time.time()
    for i, R_e in enumerate(R_e_angular):
        rho_A = R_A / R_e
        rho_B = R_B / R_e

        evals, _ = solve_angular_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
            s4_projection=True, symmetry='s4', n_states=1,
        )
        mu_vals[i] = evals[0]

        if verbose:
            elapsed = time.time() - t0
            print(f"  [{i+1}/{n_pts}] R_e={R_e:.3f}: mu={mu_vals[i]:.4f}  ({elapsed:.0f}s elapsed)")

    U_angular = (mu_vals + CENTRIFUGAL_4E) / R_e_angular**2

    # Interpolate and solve radial
    U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

    h_Re = (R_e_max - R_e_min) / (n_Re_radial + 1)
    R_e_radial = R_e_min + (np.arange(n_Re_radial) + 1) * h_Re
    V_radial = U_spline(R_e_radial)

    diag = np.ones(n_Re_radial) / h_Re**2 + V_radial
    off_diag = -0.5 * np.ones(n_Re_radial - 1) / h_Re**2

    evals_r, _ = eigh_tridiagonal(diag, off_diag, select='i', select_range=(0, 0))
    E_elec = evals_r[0]

    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    D_e = E_atoms - E_total

    t_total = time.time() - t0

    if verbose:
        E_exact = -8.0705
        D_e_exact = E_atoms - E_exact
        print(f"  R={R:.3f}: E_total={E_total:.6f}, D_e={D_e:.6f}, time={t_total:.0f}s")
        if D_e_exact > 0:
            print(f"  D_e% = {D_e/D_e_exact*100:.1f}%")

    return {
        'R': R, 'E_total': E_total, 'E_elec': E_elec, 'D_e': D_e,
        'V_NN': V_NN, 'time': t_total,
    }


def main():
    l_max = 3
    n_grid = 4  # dim = 2560

    R_values = [0.7, 1.0, 1.5, 2.0, 3.0]

    print("=" * 70)
    print(f"4-electron LiH PES scan at l_max={l_max}, n_grid={n_grid}")
    print(f"S4 [2,2] projected, dim=2560, 20 R_e points per R")
    print(f"R values: {R_values}")
    print("=" * 70)

    results = []
    t_start = time.time()

    for R in R_values:
        print(f"\n--- R = {R:.3f} bohr ---")
        result = solve_4e_lih_fast(R, l_max=l_max, n_grid=n_grid)
        results.append(result)

        # Save intermediate results
        with open('debug/track_am/lmax3_intermediate.txt', 'w') as f:
            f.write(f"l_max={l_max}, n_grid={n_grid}\n")
            f.write(f"{'R':>8s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>8s}\n")
            for r in results:
                f.write(f"{r['R']:8.3f}  {r['E_total']:12.6f}  {r['D_e']:10.6f}  {r['time']:8.0f}\n")

    t_total = time.time() - t_start

    # Summary
    print("\n" + "=" * 70)
    print(f"SUMMARY: l_max={l_max}, n_grid={n_grid}")
    print("=" * 70)

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact
    R_eq_expt = 3.015

    print(f"{'R':>8s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>8s}")
    print(f"{'bohr':>8s}  {'Ha':>12s}  {'Ha':>10s}  {'s':>8s}")
    for r in results:
        print(f"{r['R']:8.3f}  {r['E_total']:12.6f}  {r['D_e']:10.6f}  {r['time']:8.0f}")

    # Find minimum
    E_tots = [r['E_total'] for r in results]
    i_min = np.argmin(E_tots)
    R_eq = results[i_min]['R']
    has_min = i_min > 0 and i_min < len(results) - 1

    print(f"\nMinimum at R = {R_eq:.3f} bohr (index {i_min})")
    print(f"Has interior minimum: {has_min}")
    if has_min:
        R_eq_err = abs(R_eq - R_eq_expt) / R_eq_expt * 100
        print(f"R_eq error: {R_eq_err:.1f}%")
    print(f"Total wall time: {t_total:.0f}s ({t_total/60:.1f} min)")

    # Save final results
    with open('debug/track_am/lmax3_results_partial.txt', 'w') as f:
        f.write(f"4-electron LiH PES: l_max={l_max}, n_grid={n_grid}\n")
        f.write(f"S4 [2,2] projected, dim=2560\n")
        f.write(f"Total wall time: {t_total:.0f}s\n\n")
        f.write(f"{'R':>8s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>8s}\n")
        for r in results:
            f.write(f"{r['R']:8.3f}  {r['E_total']:12.6f}  {r['D_e']:10.6f}  {r['time']:8.0f}\n")
        f.write(f"\nMinimum at R = {R_eq:.3f}, has_interior_min = {has_min}\n")


if __name__ == '__main__':
    main()
