"""
Compute the H2 potential energy surface using prolate spheroidal methods.

Method 1 (fast): Eckart HF — optimize Z_eff, compute E_HF. ~20s per R point.
Method 2 (accurate): Eckart + 2x2 CI — adds sigma_u for correlation. ~140s per R.

Both approaches use exact H2+ orbitals from the prolate spheroidal lattice
with azimuthally-averaged V_ee (elliptic K integral).
"""

import numpy as np
import time
import sys
import os
from typing import Dict, List

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.prolate_scf import (
    optimize_zeff_eckart, relaxed_orbital_ci,
)

# R values to scan (bohr)
R_VALUES = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]


def compute_pes_hf(
    R_values: List[float] = None,
    N_grid: int = 40,
    xi_max_grid: float = 12.0,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    """Compute H2 PES using Eckart HF (no CI).

    For each R, optimizes Z_eff to minimize E_HF = 2*h + J_gg + 1/R.
    Fast (~20s per R) but misses electron correlation (~20% of D_e).

    Parameters
    ----------
    R_values : list of float
        Bond lengths in bohr.
    N_grid : int
        Quadrature grid size.
    xi_max_grid : float
        Maximum xi for integration.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with 'R', 'E_total', 'Z_eff', 'time' arrays.
    """
    if R_values is None:
        R_values = R_VALUES

    if verbose:
        print(f"Eckart HF PES: N_grid={N_grid}, xi_max={xi_max_grid}")
        print(f"R values: {R_values}")
        print("=" * 60)

    results_R: List[float] = []
    results_E: List[float] = []
    results_Zeff: List[float] = []
    results_time: List[float] = []

    for R in R_values:
        t0 = time.time()
        if verbose:
            print(f"  R = {R:.2f} bohr ...", end='', flush=True)

        try:
            opt = optimize_zeff_eckart(
                R=R, N_xi_solve=3000, N_grid=N_grid,
                xi_max_grid=xi_max_grid, verbose=False,
            )
            E_total = opt['E_HF_opt']
            Z_eff = opt['Z_eff_opt']
        except Exception as e:
            if verbose:
                print(f" FAILED: {e}")
            E_total = np.nan
            Z_eff = np.nan

        dt = time.time() - t0
        results_R.append(R)
        results_E.append(E_total)
        results_Zeff.append(Z_eff)
        results_time.append(dt)

        if verbose:
            print(f" E = {E_total:.6f} Ha, Z_eff = {Z_eff:.4f}, {dt:.1f}s")
            sys.stdout.flush()

    return {
        'R': np.array(results_R),
        'E_total': np.array(results_E),
        'Z_eff': np.array(results_Zeff),
        'time': np.array(results_time),
    }


def compute_pes_ci(
    R_values: List[float] = None,
    N_grid: int = 40,
    xi_max_grid: float = 12.0,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    """Compute H2 PES using Eckart + 2x2 CI (sigma_g/sigma_u).

    More accurate (~5% of D_e) but slower (~120s per R).

    Parameters
    ----------
    R_values : list of float
        Bond lengths in bohr.
    N_grid : int
        Quadrature grid size.
    xi_max_grid : float
        Maximum xi for integration.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with 'R', 'E_total', 'Z_eff', 'time' arrays.
    """
    if R_values is None:
        R_values = R_VALUES

    if verbose:
        print(f"Prolate CI PES: N_grid={N_grid}, xi_max={xi_max_grid}")
        print(f"R values: {R_values}")
        print("=" * 60)

    results_R: List[float] = []
    results_E: List[float] = []
    results_Zeff: List[float] = []
    results_time: List[float] = []

    for R in R_values:
        t0 = time.time()
        if verbose:
            print(f"  R = {R:.2f} bohr ...", end='', flush=True)

        try:
            res = relaxed_orbital_ci(
                R=R, N_xi_solve=3000, N_grid=N_grid,
                xi_max_grid=xi_max_grid, verbose=False,
            )
            E_total = res['E_total']
            Z_eff = res['Z_eff_opt']
        except Exception as e:
            if verbose:
                print(f" FAILED: {e}")
            E_total = np.nan
            Z_eff = np.nan

        dt = time.time() - t0
        results_R.append(R)
        results_E.append(E_total)
        results_Zeff.append(Z_eff)
        results_time.append(dt)

        if verbose:
            print(f" E = {E_total:.6f} Ha, Z_eff = {Z_eff:.4f}, {dt:.1f}s")
            sys.stdout.flush()

    return {
        'R': np.array(results_R),
        'E_total': np.array(results_E),
        'Z_eff': np.array(results_Zeff),
        'time': np.array(results_time),
    }


def compute_pes_hylleraas(
    R_values: List[float] = None,
    j_max: int = 2,
    l_max: int = 2,
    alpha_fixed: float = 1.19,
    l_max_neumann: int = 20,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    """Compute H2 PES using Hylleraas basis with Neumann V_ee.

    Most accurate but slowest. Uses fixed alpha to skip optimization.

    Parameters
    ----------
    R_values : list of float
        Bond lengths in bohr.
    j_max, l_max : int
        Hylleraas basis truncation.
    alpha_fixed : float
        Fixed screening parameter.
    l_max_neumann : int
        Neumann series truncation.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with 'R', 'E_total', 'time' arrays.
    """
    from geovac.hylleraas import (
        generate_basis, build_quadrature_grids, solve_hylleraas,
    )

    if R_values is None:
        R_values = R_VALUES

    grids = build_quadrature_grids(N_xi=30, N_eta=20, N_phi=1)
    basis = generate_basis(j_max=j_max, l_max=l_max, p_max=0, alpha=alpha_fixed)
    n_bf = len(basis)

    if verbose:
        print(f"Hylleraas PES: j={j_max}, l={l_max}, alpha={alpha_fixed} ({n_bf} bf)")
        print("=" * 60)

    results_R: List[float] = []
    results_E: List[float] = []
    results_time: List[float] = []

    for R in R_values:
        t0 = time.time()
        if verbose:
            print(f"  R = {R:.2f} ...", end='', flush=True)

        try:
            res = solve_hylleraas(
                basis, R, grids, verbose=False,
                vee_method='neumann', l_max_neumann=l_max_neumann,
            )
            E_total = res['E_total']
        except Exception as e:
            if verbose:
                print(f" FAILED: {e}")
            E_total = np.nan

        dt = time.time() - t0
        results_R.append(R)
        results_E.append(E_total)
        results_time.append(dt)

        if verbose:
            print(f" E = {E_total:.6f} Ha, {dt:.1f}s")
            sys.stdout.flush()

    return {
        'R': np.array(results_R),
        'E_total': np.array(results_E),
        'time': np.array(results_time),
    }


def save_pes(results: Dict[str, np.ndarray], filename: str) -> None:
    """Save PES results to a text file."""
    with open(filename, 'w') as f:
        f.write("# H2 PES\n")
        f.write(f"# {'R_bohr':>10s}  {'E_total_Ha':>14s}\n")
        for i in range(len(results['R'])):
            f.write(f"  {results['R'][i]:10.4f}  {results['E_total'][i]:14.8f}\n")


if __name__ == '__main__':
    results = compute_pes_hf(verbose=True)
    outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'h2_pes_data.txt')
    save_pes(results, outfile)
    print(f"\nPES data saved to {outfile}")
