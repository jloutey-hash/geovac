"""
2D variational solver for 4-electron LiH in mol-frame hyperspherical coordinates.

Track AR: Bypasses the adiabatic approximation by treating the full
(R_e, angular) problem variationally. Assembles the tensor product
Hamiltonian:

    H_2D = T_Re (x) I_ang + [H_ang(R_e) + 99/8 * I] / R_e^2

where H_ang(R_e) is the full angular Hamiltonian (kinetic + centrifugal +
V_ee + V_nuc) evaluated at each radial grid point, and T_Re is the
hyperradial kinetic energy on an FD grid.

This follows the Level 4 solve_direct_2d() pattern in level4_multichannel.py.

Key differences from Level 4:
  - 3 hyperangles (alpha1, alpha2, alpha3) instead of 1 (alpha)
  - S4 [2,2] singlet projection instead of gerade constraint
  - Angular dimension: n_ch_S4 * n_grid^3 instead of n_ch * n_alpha
  - Centrifugal constant: 99/8 instead of 15/8

References:
  - Level 4 2D solver: geovac/level4_multichannel.py (solve_direct_2d)
  - N-electron adiabatic solver: geovac/n_electron_solver.py
  - Track AP scoping: debug/track_ap/analysis.md
  - Track AO failure analysis: debug/track_ao/results.txt
"""

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, eye as sp_eye
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from typing import Tuple, Dict, Optional
import time

from geovac.n_electron_solver import (
    build_angular_hamiltonian_4e_multichannel,
    build_angular_hamiltonian_4e_parity_multichannel,
    CENTRIFUGAL_4E,
)


def solve_4e_lih_2d(
    R: float,
    Z_A: float = 3.0,
    Z_B: float = 1.0,
    n_grid: int = 6,
    l_max: int = 2,
    n_Re: int = 200,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    origin: str = 'midpoint',
    symmetry: str = 'parity',
    verbose: bool = True,
) -> Dict:
    """2D variational solver for 4-electron LiH.

    Builds the full (R_e, angular) tensor product Hamiltonian and finds
    the lowest eigenvalue, bypassing the adiabatic approximation entirely.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle. Angular dim = n_ch_S4 * n_grid^3.
    l_max : int
        Maximum angular momentum per electron.
    n_Re : int
        FD grid points for hyperradius R_e.
    R_e_min, R_e_max : float
        Hyperradial boundaries.
    origin : str
        'midpoint' or 'charge_center'.
    symmetry : str
        'parity' (half-grid pair antisymmetry) or 's4' (channel-space [2,2]).
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
        Contains E_elec, E_total, D_e, timing, etc.
    """
    t0 = time.time()

    # Origin shift
    if origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:
        z0 = 0.0

    # Build angular Hamiltonian at a reference R_e to determine dimensions
    R_A = R / 2.0 - z0
    R_B = R / 2.0 + z0
    R_e_ref = 2.0
    rho_A_ref = R_A / R_e_ref
    rho_B_ref = R_B / R_e_ref

    if symmetry == 'parity':
        H_ang_ref, N_ang, n_ch = build_angular_hamiltonian_4e_parity_multichannel(
            rho_A_ref, rho_B_ref, R_e_ref, Z_A, Z_B, n_grid, l_max,
        )
    else:
        H_ang_ref, N_ang, n_ch = build_angular_hamiltonian_4e_multichannel(
            rho_A_ref, rho_B_ref, R_e_ref, Z_A, Z_B, n_grid, l_max,
            s4_projection=True,
        )

    if N_ang == 0:
        if verbose:
            print("  No S4 [2,2] channels at this l_max. Returning zero energy.")
        return {
            'E_elec': 0.0, 'E_total': 0.0, 'D_e': 0.0,
            'R': R, 'l_max': l_max, 'method': '2d',
        }

    N_total = n_Re * N_ang

    if verbose:
        print(f"4-electron 2D solver at R = {R:.4f} bohr")
        print(f"  Z_A={Z_A}, Z_B={Z_B}, l_max={l_max}")
        print(f"  n_grid={n_grid}, angular_dim={N_ang}, n_ch={n_ch}")
        print(f"  n_Re={n_Re}, R_e range=[{R_e_min:.2f}, {R_e_max:.2f}]")
        print(f"  Total 2D dim: {N_total}")
        print(f"  symmetry={symmetry}, origin={origin}")

    # R_e FD grid
    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_grid = R_e_min + (np.arange(n_Re) + 1) * h_Re

    # Build sparse 2D Hamiltonian in COO format
    # Index: r * N_ang + a  (r = radial index, a = angular index)
    kinetic_diag = 1.0 / h_Re**2
    kinetic_off = -0.5 / h_Re**2

    rows_list = []
    cols_list = []
    vals_list = []

    # Radial kinetic energy: T_Re (x) I_ang
    # Diagonal: kinetic_diag at each (r, a) point
    ang_idx = np.arange(N_ang)
    for r in range(n_Re):
        diag_idx = r * N_ang + ang_idx
        rows_list.append(diag_idx)
        cols_list.append(diag_idx)
        vals_list.append(np.full(N_ang, kinetic_diag))

        if r < n_Re - 1:
            idx1 = r * N_ang + ang_idx
            idx2 = (r + 1) * N_ang + ang_idx
            rows_list.extend([idx1, idx2])
            cols_list.extend([idx2, idx1])
            vals_list.extend([np.full(N_ang, kinetic_off)] * 2)

    t_build_start = time.time()

    # Angular blocks: [H_ang(R_e) + 99/8 * I] / R_e^2 at each R_e
    for r, R_e in enumerate(R_e_grid):
        rho_A = R_A / R_e
        rho_B = R_B / R_e

        if symmetry == 'parity':
            H_ang, dim_check, _ = build_angular_hamiltonian_4e_parity_multichannel(
                rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
            )
        else:
            H_ang, dim_check, _ = build_angular_hamiltonian_4e_multichannel(
                rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
                s4_projection=True,
            )

        # Add centrifugal constant and scale by 1/R_e^2
        np.fill_diagonal(H_ang, H_ang.diagonal() + CENTRIFUGAL_4E)
        H_ang *= (1.0 / R_e**2)

        # Extract nonzero entries and add to COO
        nz_i, nz_j = np.nonzero(np.abs(H_ang) > 1e-15)
        if len(nz_i) > 0:
            offset = r * N_ang
            rows_list.append(offset + nz_i)
            cols_list.append(offset + nz_j)
            vals_list.append(H_ang[nz_i, nz_j])

        if verbose and (r % 50 == 0 or r == n_Re - 1):
            elapsed = time.time() - t_build_start
            print(f"  Building angular block {r+1}/{n_Re} ({elapsed:.1f}s)")

    t_assemble = time.time()
    if verbose:
        print(f"  Angular block assembly: {t_assemble - t_build_start:.1f}s")

    # Assemble sparse matrix
    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)

    H_2d = coo_matrix(
        (all_vals, (all_rows, all_cols)),
        shape=(N_total, N_total),
    ).tocsr()

    if verbose:
        nnz = H_2d.nnz
        fill = nnz / N_total**2 * 100
        print(f"  Sparse matrix: {N_total}x{N_total}, nnz={nnz}, fill={fill:.2f}%")

    # Find lowest eigenvalue
    # Estimate sigma for shift-invert
    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    V_NN = Z_A * Z_B / R
    sigma = E_atoms - V_NN - 0.5  # target below atomic limit

    t_solve_start = time.time()

    try:
        if verbose:
            print(f"  Eigensolve: shift-invert with sigma={sigma:.4f}")
        evals, evecs = eigsh(H_2d, k=3, sigma=sigma, which='LM',
                             maxiter=5000)
        E_elec = np.min(evals)
        if verbose:
            print(f"  Eigenvalues: {np.sort(evals)}")
    except Exception as e:
        if verbose:
            print(f"  Shift-invert failed ({e}), trying SA mode")
        try:
            evals, evecs = eigsh(H_2d, k=6, which='SA', maxiter=10000)
            E_elec = evals[0]
            if verbose:
                print(f"  Eigenvalues (SA): {np.sort(evals[:6])}")
        except Exception as e2:
            if verbose:
                print(f"  SA also failed ({e2}), using dense solver on subblock")
            # Last resort: dense solve if dimension is manageable
            if N_total <= 10000:
                evals_all, _ = eigh(H_2d.toarray())
                E_elec = evals_all[0]
            else:
                raise RuntimeError(f"Eigensolve failed for dim={N_total}: {e2}")

    t_solve = time.time()
    if verbose:
        print(f"  Eigensolve: {t_solve - t_solve_start:.1f}s")

    # Physical results
    E_total = E_elec + V_NN

    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact  # 0.0924 Ha
    D_e = E_atoms - E_total

    t_total = time.time() - t0

    if verbose:
        print(f"\n  === 2D Variational Results (l_max={l_max}) ===")
        print(f"  E_elec  = {E_elec:.6f} Ha")
        print(f"  V_NN    = {V_NN:.6f} Ha")
        print(f"  E_total = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
        print(f"  E_atoms = {E_atoms:.6f} Ha")
        print(f"  D_e     = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        if D_e_exact > 0:
            print(f"  D_e %   = {D_e / D_e_exact * 100:.1f}%")
        print(f"  Total time: {t_total:.1f}s")

    return {
        'E_elec': E_elec,
        'E_total': E_total,
        'V_NN': V_NN,
        'D_e': D_e,
        'D_e_pct': (D_e / D_e_exact * 100) if D_e_exact > 0 else None,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'l_max': l_max,
        'n_grid': n_grid,
        'n_Re': n_Re,
        'n_channels': n_ch,
        'angular_dim': N_ang,
        'total_dim': N_total,
        'E_atoms': E_atoms,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'z0': z0,
        'origin': origin,
        'symmetry': symmetry,
        'method': '2d',
        'time_build': t_assemble - t0,
        'time_solve': t_solve - t_solve_start,
        'time_total': t_total,
    }


def scan_pes_4e_lih_2d(
    R_values: np.ndarray = None,
    Z_A: float = 3.0,
    Z_B: float = 1.0,
    n_grid: int = 6,
    l_max: int = 2,
    n_Re: int = 200,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    origin: str = 'midpoint',
    symmetry: str = 'parity',
    verbose: bool = True,
    output_dir: str = None,
) -> Dict:
    """Scan LiH PES with the 2D variational solver.

    Parameters
    ----------
    R_values : ndarray or None
        Internuclear distances (bohr). Default: standard scan grid.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum.
    n_Re : int
        Radial FD grid points.
    R_e_min, R_e_max : float
        Hyperradial range.
    origin : str
        'midpoint' or 'charge_center'.
    symmetry : str
        'parity' or 's4'.
    verbose : bool
        Print progress.
    output_dir : str or None
        Directory to save results.

    Returns
    -------
    result : dict
    """
    if R_values is None:
        R_values = np.array([0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0,
                             2.5, 3.0, 4.0, 6.0])

    n_R = len(R_values)
    E_total = np.zeros(n_R)
    E_elec = np.zeros(n_R)
    D_e = np.zeros(n_R)
    times = np.zeros(n_R)

    t_start = time.time()

    for j, R in enumerate(R_values):
        if verbose:
            print(f"\n{'='*60}")
            print(f"R = {R:.4f} bohr ({j+1}/{n_R})")
            print(f"{'='*60}")

        result = solve_4e_lih_2d(
            R, Z_A, Z_B, n_grid=n_grid, l_max=l_max,
            n_Re=n_Re, R_e_min=R_e_min, R_e_max=R_e_max,
            origin=origin, symmetry=symmetry, verbose=verbose,
        )
        E_total[j] = result['E_total']
        E_elec[j] = result['E_elec']
        D_e[j] = result['D_e']
        times[j] = result['time_total']

    t_total = time.time() - t_start

    # Find minimum
    i_min = np.argmin(E_total)
    R_eq = R_values[i_min]
    E_min = E_total[i_min]
    D_e_min = D_e[i_min]

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact
    R_eq_expt = 3.015

    has_minimum = (0 < i_min < n_R - 1)

    if verbose:
        print(f"\n{'='*60}")
        print(f"2D VARIATIONAL PES SCAN (l_max={l_max})")
        print(f"{'='*60}")
        print(f"  {'R':>6s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>8s}")
        print(f"  {'bohr':>6s}  {'Ha':>12s}  {'Ha':>10s}  {'s':>8s}")
        print(f"  {'-'*6}  {'-'*12}  {'-'*10}  {'-'*8}")
        for j in range(n_R):
            marker = " <-- min" if j == i_min else ""
            print(f"  {R_values[j]:6.3f}  {E_total[j]:12.6f}  "
                  f"{D_e[j]:10.6f}  {times[j]:8.1f}{marker}")
        print()
        if has_minimum:
            R_eq_err = abs(R_eq - R_eq_expt) / R_eq_expt * 100
            print(f"  EQUILIBRIUM FOUND: R_eq = {R_eq:.3f} bohr")
            print(f"    (expt: {R_eq_expt:.3f}, error: {R_eq_err:.1f}%)")
            print(f"  E_min  = {E_min:.6f} Ha  (exact: {E_exact:.6f})")
            print(f"  D_e    = {D_e_min:.6f} Ha  (exact: {D_e_exact:.6f})")
            if D_e_exact > 0:
                print(f"  D_e %  = {D_e_min / D_e_exact * 100:.1f}%")
        else:
            print(f"  NO EQUILIBRIUM FOUND (min at boundary R={R_eq:.3f})")
        print(f"  Total time: {t_total:.1f}s")

    scan_result = {
        'R': R_values,
        'E_total': E_total,
        'E_elec': E_elec,
        'D_e': D_e,
        'R_eq': R_eq,
        'E_min': E_min,
        'D_e_min': D_e_min,
        'has_minimum': has_minimum,
        'R_eq_expt': R_eq_expt,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'time_per_point': times,
        'time_total': t_total,
        'n_grid': n_grid,
        'l_max': l_max,
        'n_Re': n_Re,
        'n_channels': result['n_channels'],
        'angular_dim': result['angular_dim'],
        'total_dim': result['total_dim'],
        'origin': origin,
        'symmetry': symmetry,
        'method': '2d',
    }

    if output_dir is not None:
        import os
        os.makedirs(output_dir, exist_ok=True)
        np.savez(
            os.path.join(output_dir, f'pes_4e_lih_2d_lmax{l_max}.npz'),
            **scan_result,
        )
        # Save human-readable summary
        with open(os.path.join(output_dir, f'pes_4e_lih_2d_lmax{l_max}.txt'), 'w') as f:
            f.write(f"Track AR: 2D Variational PES for 4-electron LiH\n")
            f.write(f"l_max={l_max}, n_grid={n_grid}, n_Re={n_Re}\n")
            f.write(f"symmetry={symmetry}, origin={origin}\n")
            f.write(f"angular_dim={result['angular_dim']}, total_dim={result['total_dim']}\n\n")
            f.write(f"{'R':>8s}  {'E_total':>14s}  {'D_e':>12s}  {'time':>8s}\n")
            f.write(f"{'bohr':>8s}  {'Ha':>14s}  {'Ha':>12s}  {'s':>8s}\n")
            for j in range(n_R):
                f.write(f"{R_values[j]:8.4f}  {E_total[j]:14.8f}  "
                        f"{D_e[j]:12.8f}  {times[j]:8.1f}\n")
            f.write(f"\nR_eq = {R_eq:.4f} bohr (expt {R_eq_expt:.3f})\n")
            if has_minimum:
                R_eq_err = abs(R_eq - R_eq_expt) / R_eq_expt * 100
                f.write(f"R_eq error = {R_eq_err:.1f}%\n")
            f.write(f"E_min = {E_min:.8f} Ha (exact {E_exact:.4f})\n")
            f.write(f"D_e = {D_e_min:.8f} Ha (exact {D_e_exact:.4f})\n")

    return scan_result
