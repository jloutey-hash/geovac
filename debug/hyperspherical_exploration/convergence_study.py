"""
Convergence study for the hyperspherical adiabatic He solver.

Tests convergence in:
1. N_R (radial grid points)
2. l_max (partial wave truncation)
3. n_alpha (angular FD grid)

Also generates the adiabatic potential curve plot.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import time

# Add project root to path
sys.path.insert(0, r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric')

from geovac.hyperspherical_angular import solve_angular
from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
    plot_adiabatic_curves,
)
from geovac.hyperspherical_radial import solve_helium

E_EXACT = -2.903724


def convergence_n_alpha() -> None:
    """Convergence with angular FD grid points."""
    print("=" * 60)
    print("CONVERGENCE: n_alpha (angular FD grid)")
    print("=" * 60)
    print(f"{'n_alpha':>8} {'E (Ha)':>12} {'Error %':>10} {'Time (s)':>10}")
    print("-" * 45)

    for n_alpha in [50, 80, 100, 150, 200]:
        t0 = time.time()
        result = solve_helium(
            Z=2.0, l_max=3, n_alpha=n_alpha,
            N_R_angular=120, R_min=0.05, R_max=20.0,
            N_R_radial=2000, verbose=False,
        )
        dt = time.time() - t0
        E = result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT) * 100
        print(f"{n_alpha:>8} {E:>12.6f} {err:>10.4f} {dt:>10.2f}")


def convergence_l_max() -> None:
    """Convergence with partial wave truncation."""
    print("\n" + "=" * 60)
    print("CONVERGENCE: l_max (partial waves)")
    print("=" * 60)
    print(f"{'l_max':>8} {'E (Ha)':>12} {'Error %':>10} {'Time (s)':>10}")
    print("-" * 45)

    for l_max in [0, 1, 2, 3, 4]:
        t0 = time.time()
        result = solve_helium(
            Z=2.0, l_max=l_max, n_alpha=150,
            N_R_angular=100, R_min=0.05, R_max=20.0,
            N_R_radial=2000, verbose=False,
        )
        dt = time.time() - t0
        E = result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT) * 100
        print(f"{l_max:>8} {E:>12.6f} {err:>10.4f} {dt:>10.2f}")


def convergence_N_R() -> None:
    """Convergence with radial grid points."""
    print("\n" + "=" * 60)
    print("CONVERGENCE: N_R (radial grid)")
    print("=" * 60)
    print(f"{'N_R':>8} {'E (Ha)':>12} {'Error %':>10}")
    print("-" * 35)

    # Precompute the potential curve once
    from geovac.hyperspherical_radial import solve_radial
    from scipy.interpolate import CubicSpline

    R_grid_ang = np.concatenate([
        np.linspace(0.1, 1.0, 50),
        np.linspace(1.0, 5.0, 50),
        np.linspace(5.0, 20.0, 51),
    ])
    R_grid_ang = np.unique(R_grid_ang)

    mu = compute_adiabatic_curve(R_grid_ang, Z=2.0, l_max=3, n_alpha=200)
    V_eff = effective_potential(R_grid_ang, mu[0])
    V_spline = CubicSpline(R_grid_ang, V_eff, extrapolate=True)

    for N_R in [500, 1000, 2000, 4000, 8000]:
        E, _, _ = solve_radial(V_spline, R_min=0.05, R_max=20.0,
                               N_R=N_R, n_states=1)
        err = abs(E[0] - E_EXACT) / abs(E_EXACT) * 100
        print(f"{N_R:>8} {E[0]:>12.6f} {err:>10.4f}")


def generate_potential_plot() -> None:
    """Generate the adiabatic potential curve plot."""
    print("\n" + "=" * 60)
    print("GENERATING ADIABATIC POTENTIAL PLOT")
    print("=" * 60)

    R_grid = np.concatenate([
        np.linspace(0.2, 1.0, 30),
        np.linspace(1.0, 5.0, 30),
        np.linspace(5.0, 15.0, 20),
    ])
    R_grid = np.unique(R_grid)

    mu = compute_adiabatic_curve(R_grid, Z=2.0, l_max=3, n_alpha=200,
                                 n_channels=1)

    save_path = r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric\debug\plots\hyperspherical_Veff.png'
    plot_adiabatic_curves(R_grid, mu, Z=2.0, save_path=save_path)
    print(f"Plot saved to {save_path}")


if __name__ == "__main__":
    convergence_l_max()
    convergence_N_R()
    convergence_n_alpha()
    generate_potential_plot()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Exact He:    {E_EXACT:.6f} Ha")
    print(f"S^3 lattice: -2.8508 Ha (1.82% error)")
    print(f"Literature single-channel: ~-2.879 Ha (0.85% error)")
    print(f"Our implementation: see tables above")
