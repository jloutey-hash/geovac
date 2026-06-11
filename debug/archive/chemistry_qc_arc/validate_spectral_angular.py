"""Validate spectral angular solver against FD for Level 4.

Compares eigenvalues at multiple (rho, R_e) points and benchmarks speedup.
"""
import sys
sys.path.insert(0, '.')

import time
import numpy as np
from geovac.level4_multichannel import (
    solve_angular_multichannel, compute_adiabatic_curve_mc, _channel_list,
)
from geovac.level4_spectral_angular import (
    Level4SpectralAngular, solve_angular_spectral,
    compute_adiabatic_curve_spectral,
)


def test_single_point():
    """Compare single-point eigenvalues at several (rho, R_e) values."""
    R = 1.4
    test_points = [
        (0.5, 1.4),   # near minimum
        (1.0, 0.7),   # small R_e
        (2.0, 0.35),  # large rho
        (0.3, 2.33),  # large R_e
    ]

    l_max = 2
    n_alpha = 200
    n_basis = 10

    print("=== Single-point eigenvalue comparison ===")
    print(f"l_max={l_max}, n_alpha={n_alpha} (FD), n_basis={n_basis} (spectral)")
    print()

    # Build spectral solver once
    solver = Level4SpectralAngular(l_max=l_max, n_basis=n_basis, n_quad=100)

    for R_e_val in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
        rho = R / (2.0 * R_e_val)

        # FD reference
        evals_fd, _, _, _ = solve_angular_multichannel(
            rho, R_e_val, l_max, Z=1.0, n_alpha=n_alpha, n_eig=3,
        )

        # Spectral
        evals_sp, _ = solver.solve(rho, R_e_val, n_eig=3)

        delta0 = abs(evals_sp[0] - evals_fd[0])
        print(f"R_e={R_e_val:5.1f}  rho={rho:.4f}  "
              f"mu0_fd={evals_fd[0]:10.5f}  mu0_sp={evals_sp[0]:10.5f}  "
              f"delta={delta0:.2e}")


def test_adiabatic_curve():
    """Compare full adiabatic curves."""
    R = 1.4
    l_max = 2

    # Use a subset of R_e grid for speed
    R_e_grid = np.linspace(0.5, 8.0, 50)

    print("\n=== Adiabatic curve comparison ===")

    # FD
    t0 = time.perf_counter()
    U_fd = compute_adiabatic_curve_mc(R, R_e_grid, l_max, Z=1.0, n_alpha=200)
    t_fd = time.perf_counter() - t0

    # Spectral
    t0 = time.perf_counter()
    U_sp = compute_adiabatic_curve_spectral(
        R, R_e_grid, l_max, Z=1.0, n_basis=10, n_quad=100,
    )
    t_sp = time.perf_counter() - t0

    delta_U = np.abs(U_fd - U_sp)
    rel_err = delta_U / (np.abs(U_fd) + 1e-10)

    print(f"FD time:       {t_fd*1000:.0f} ms")
    print(f"Spectral time: {t_sp*1000:.0f} ms")
    print(f"Speedup:       {t_fd/t_sp:.1f}x")
    print(f"Max |U_fd - U_sp|: {np.max(delta_U):.2e} Ha")
    print(f"Max rel error:     {np.max(rel_err):.2e}")
    print(f"U_min (FD):        {np.min(U_fd):.6f} Ha at R_e={R_e_grid[np.argmin(U_fd)]:.2f}")
    print(f"U_min (spectral):  {np.min(U_sp):.6f} Ha at R_e={R_e_grid[np.argmin(U_sp)]:.2f}")


def test_basis_convergence():
    """Test convergence with n_basis."""
    R = 1.4
    R_e = 1.5
    rho = R / (2.0 * R_e)
    l_max = 2

    # FD reference
    evals_fd, _, _, _ = solve_angular_multichannel(
        rho, R_e, l_max, Z=1.0, n_alpha=200, n_eig=1,
    )
    mu0_fd = evals_fd[0]

    print("\n=== Basis convergence ===")
    print(f"FD reference: mu0 = {mu0_fd:.8f}")
    print()

    for nb in [5, 8, 10, 15, 20, 25]:
        solver = Level4SpectralAngular(l_max=l_max, n_basis=nb, n_quad=100)
        evals, _ = solver.solve(rho, R_e, n_eig=1)
        delta = abs(evals[0] - mu0_fd)
        print(f"n_basis={nb:3d}  dim={solver._total_dim:4d}  "
              f"mu0={evals[0]:12.8f}  delta={delta:.2e}")


def test_timing():
    """Detailed timing comparison."""
    R = 1.4
    l_max = 2

    R_e_grid = np.linspace(0.5, 8.0, 50)

    print("\n=== Detailed timing ===")

    # Time spectral solver components
    solver = Level4SpectralAngular(l_max=l_max, n_basis=10, n_quad=100)
    print(f"Solver init (once): n_ch={solver.n_ch}, dim={solver._total_dim}")

    # Time nuclear coupling computation
    rho = R / (2.0 * 1.5)
    t0 = time.perf_counter()
    for _ in range(100):
        V_nuc = solver._compute_nuclear_spectral(rho, 1.5)
    t_nuc = (time.perf_counter() - t0) / 100
    print(f"Nuclear coupling:  {t_nuc*1000:.2f} ms per point")

    # Time eigensolve
    dim = solver._total_dim
    H = np.random.randn(dim, dim)
    H = (H + H.T) / 2
    t0 = time.perf_counter()
    for _ in range(1000):
        evals, evecs = np.linalg.eigh(H)
    t_eig = (time.perf_counter() - t0) / 1000
    print(f"Eigensolve ({dim}x{dim}): {t_eig*1000:.3f} ms per point")

    # Total per point estimate
    t_total = t_nuc + t_eig
    print(f"Total per point:   {t_total*1000:.2f} ms")
    print(f"FD per point:      ~327 ms")
    print(f"Estimated speedup: {327.0 / (t_total * 1000):.0f}x")


if __name__ == '__main__':
    test_single_point()
    test_basis_convergence()
    test_adiabatic_curve()
    test_timing()
