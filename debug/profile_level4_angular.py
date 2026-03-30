"""Profile the Level 4 angular solver bottleneck.

Measures:
1. Total angular sweep time for H2 at R=1.4, l_max=2
2. Per-iteration breakdown (assembly vs eigensolve)
3. Matrix dimensions and channel count
"""
import sys
sys.path.insert(0, '.')

import time
import numpy as np
from geovac.level4_multichannel import (
    build_angular_hamiltonian, solve_angular_multichannel,
    compute_adiabatic_curve_mc, _channel_list,
)


def profile_single_solve():
    """Profile a single angular solve at one R_e point."""
    l_max = 2
    n_alpha = 200
    R = 1.4
    R_e = 1.5
    rho = R / (2.0 * R_e)

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    channels = _channel_list(l_max, homonuclear=True)
    n_ch = len(channels)
    N = n_ch * n_alpha

    print(f"=== Single angular solve profile ===")
    print(f"l_max={l_max}, n_alpha={n_alpha}, n_ch={n_ch}, N={N}")
    print(f"R={R}, R_e={R_e}, rho={rho:.4f}")
    print()

    # Time assembly
    t0 = time.perf_counter()
    for _ in range(10):
        H = build_angular_hamiltonian(alpha, rho, R_e, l_max, Z=1.0)
    t_assembly = (time.perf_counter() - t0) / 10
    print(f"Assembly:    {t_assembly*1000:.1f} ms  (matrix {N}x{N})")

    # Time eigensolve
    from scipy.linalg import eigh
    t0 = time.perf_counter()
    for _ in range(10):
        evals, evecs = eigh(H)
    t_eigen = (time.perf_counter() - t0) / 10
    print(f"Eigensolve:  {t_eigen*1000:.1f} ms")
    print(f"Total:       {(t_assembly + t_eigen)*1000:.1f} ms per R_e point")
    print(f"mu_0 = {evals[0]:.6f}")
    print()


def profile_full_sweep():
    """Profile the full 130-point adiabatic sweep."""
    l_max = 2
    n_alpha = 200
    R = 1.4

    # Build R_e grid (same as production)
    grids = [
        np.linspace(0.3, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, 15.0, 20),
    ]
    R_e_grid = np.unique(np.concatenate(grids))
    n_Re = len(R_e_grid)

    channels = _channel_list(l_max, homonuclear=True)
    n_ch = len(channels)

    print(f"=== Full adiabatic sweep profile ===")
    print(f"l_max={l_max}, n_alpha={n_alpha}, n_ch={n_ch}")
    print(f"R={R}, n_Re={n_Re}")
    print()

    t0 = time.perf_counter()
    U = compute_adiabatic_curve_mc(R, R_e_grid, l_max, Z=1.0, n_alpha=n_alpha)
    t_total = time.perf_counter() - t0

    print(f"Total sweep: {t_total*1000:.0f} ms ({n_Re} points)")
    print(f"Per point:   {t_total/n_Re*1000:.1f} ms")
    print(f"U_min = {np.min(U):.6f} Ha at R_e = {R_e_grid[np.argmin(U)]:.2f}")
    print()

    return t_total, n_Re


if __name__ == '__main__':
    profile_single_solve()
    t_total, n_Re = profile_full_sweep()

    # Estimate spectral solver target
    n_basis = 10
    channels = _channel_list(2, homonuclear=True)
    n_ch = len(channels)
    N_spectral = n_ch * n_basis
    print(f"=== Speedup estimate ===")
    print(f"FD matrix:       {n_ch * 200} x {n_ch * 200}")
    print(f"Spectral matrix: {N_spectral} x {N_spectral}")
    print(f"Eigensolve ratio: ({n_ch*200}/{N_spectral})^3 = {(n_ch*200/N_spectral)**3:.0f}x")
