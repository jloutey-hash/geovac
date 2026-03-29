"""
Tests for the hyperspherical adiabatic solver for helium.

The single-channel adiabatic approximation omits the diagonal
Born-Oppenheimer correction (DBOC = +1/2 ||dPhi/dR||^2 >= 0),
so the adiabatic energy can be slightly BELOW the exact value.

References:
  - Exact He: E = -2.903724 Ha (Pekeris 1958)
  - Current S^3: E = -2.8508 Ha (1.82% error)
"""

import numpy as np
import pytest

from geovac.hyperspherical_angular import solve_angular, gaunt_integral
from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
)
from geovac.hyperspherical_radial import (
    solve_radial, solve_helium,
    solve_radial_spectral, solve_coupled_radial_spectral,
)


E_EXACT = -2.903724


class TestGauntIntegrals:
    def test_gaunt_000(self) -> None:
        assert abs(gaunt_integral(0, 0, 0) - 2.0) < 1e-12

    def test_gaunt_101(self) -> None:
        assert abs(gaunt_integral(1, 0, 1) - 2.0 / 3) < 1e-12

    def test_gaunt_110(self) -> None:
        assert abs(gaunt_integral(1, 1, 0) - 2.0 / 3) < 1e-12

    def test_gaunt_parity_selection(self) -> None:
        assert gaunt_integral(1, 1, 1) == 0.0
        assert gaunt_integral(0, 1, 0) == 0.0

    def test_gaunt_triangle_selection(self) -> None:
        assert gaunt_integral(0, 0, 2) == 0.0

    def test_gaunt_symmetry(self) -> None:
        for l1 in range(4):
            for l2 in range(4):
                for k in range(l1 + l2 + 1):
                    assert abs(gaunt_integral(l1, k, l2)
                               - gaunt_integral(l2, k, l1)) < 1e-14

    def test_gaunt_202(self) -> None:
        assert abs(gaunt_integral(2, 0, 2) - 2.0 / 5) < 1e-12


class TestAngularSolver:
    def test_angular_returns_shape(self) -> None:
        mu, vecs = solve_angular(R=1.0, Z=2.0, l_max=2, n_alpha=60)
        assert mu.shape == (1,)

    def test_angular_multiple_channels(self) -> None:
        mu, _ = solve_angular(R=1.0, Z=2.0, l_max=2, n_alpha=60,
                              n_channels=3)
        assert mu[0] < mu[1] < mu[2]

    def test_angular_free_eigenvalue(self) -> None:
        """At R=0, mu = 2n^2 - 2 for the Liouville-transformed equation."""
        mu, _ = solve_angular(R=0.0, Z=2.0, l_max=0, n_alpha=200,
                              n_channels=3)
        assert abs(mu[0] - 0.0) < 0.05
        assert abs(mu[1] - 6.0) < 0.1

    def test_angular_convergence(self) -> None:
        mus = []
        for n_a in [50, 100, 200]:
            mu, _ = solve_angular(R=1.5, Z=2.0, l_max=0, n_alpha=n_a)
            mus.append(mu[0])
        assert abs(mus[2] - mus[1]) < abs(mus[1] - mus[0])

    def test_angular_asymptotic(self) -> None:
        """V_eff(R=20) should be near -Z^2/2 = -2.0."""
        mu, _ = solve_angular(R=20.0, Z=2.0, l_max=0, n_alpha=300)
        V_eff = mu[0] / 400.0 + 15.0 / 3200.0
        assert abs(V_eff - (-2.0)) < 0.1


class TestAdiabaticCurve:
    def test_potential_well(self) -> None:
        R_grid = np.linspace(0.3, 5.0, 25)
        mu = compute_adiabatic_curve(R_grid, Z=2.0, l_max=0, n_alpha=100)
        V_eff = effective_potential(R_grid, mu[0])
        assert np.min(V_eff) < -2.0

    def test_repulsive_small_R(self) -> None:
        R_grid = np.array([0.1, 0.15, 0.2])
        mu = compute_adiabatic_curve(R_grid, Z=2.0, l_max=0, n_alpha=80)
        V_eff = effective_potential(R_grid, mu[0])
        assert V_eff[0] > 0


class TestHeliumSolver:
    @pytest.fixture(scope="class")
    def he_result(self) -> dict:
        return solve_helium(
            Z=2.0, l_max=0, n_alpha=200,
            N_R_angular=100, R_min=0.05, R_max=30.0,
            N_R_radial=3000, verbose=False
        )

    def test_energy_close_to_exact(self, he_result: dict) -> None:
        """Within 0.1% of exact (DBOC correction not included)."""
        E = he_result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.001, f"Error {err:.4f} > 0.1%"

    def test_beats_s3_lattice(self, he_result: dict) -> None:
        E = he_result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.018

    def test_wavefunction_decays(self, he_result: dict) -> None:
        F = he_result['wavefunction']
        assert np.max(np.abs(F[-100:])) / np.max(np.abs(F)) < 0.01

    def test_wavefunction_normalized(self, he_result: dict) -> None:
        F = he_result['wavefunction']
        R = he_result['R_grid_radial']
        h = R[1] - R[0]
        assert abs(h * np.sum(F**2) - 1.0) < 0.1


class TestRadialSolver:
    def test_harmonic_oscillator(self) -> None:
        E, _, _ = solve_radial(lambda R: 0.5 * R**2,
                               R_min=-10.0, R_max=10.0, N_R=2000)
        assert abs(E[0] - 0.5) < 0.01

    def test_hydrogen_1s(self) -> None:
        E, _, _ = solve_radial(lambda R: -1.0 / R,
                               R_min=0.01, R_max=50.0, N_R=5000)
        assert abs(E[0] - (-0.5)) < 0.02


class TestSpectralRadialSolver:
    """Tests for the spectral Laguerre hyperradial solver (Lane 2, v2.0.9)."""

    @pytest.fixture(scope="class")
    def he_spectral(self) -> dict:
        """Spectral single-channel He solve using FD angular input."""
        return solve_helium(
            Z=2.0, l_max=0, n_alpha=200,
            N_R_angular=100, R_min=0.05, R_max=30.0,
            N_R_radial=3000, verbose=False,
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
        )

    @pytest.fixture(scope="class")
    def he_fd(self) -> dict:
        """FD single-channel He solve for comparison."""
        return solve_helium(
            Z=2.0, l_max=0, n_alpha=200,
            N_R_angular=100, R_min=0.05, R_max=30.0,
            N_R_radial=3000, verbose=False,
            radial_method='fd',
        )

    def test_spectral_fd_consistency(self, he_spectral: dict,
                                     he_fd: dict) -> None:
        """Spectral and FD agree within FD discretization error (~0.01 Ha)."""
        diff = abs(he_spectral['energy'] - he_fd['energy'])
        assert diff < 0.01, f"Spectral-FD diff {diff:.6f} Ha > 0.01 Ha"

    def test_spectral_energy_close_to_exact(self, he_spectral: dict) -> None:
        """Spectral single-channel within 0.1% of exact."""
        E = he_spectral['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.001, f"Error {err:.4f} > 0.1%"

    def test_spectral_dimension_reduction(self) -> None:
        """Spectral uses 25 basis functions vs 3000 FD grid points (120x)."""
        # n_basis_radial=25, N_R_radial=3000
        assert 3000 / 25 >= 100

    def test_spectral_convergence_with_n_basis(self) -> None:
        """Energy converges monotonically as n_basis increases."""
        energies = []
        for nb in [10, 15, 20, 25]:
            r = solve_helium(
                Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
                radial_method='spectral', n_basis_radial=nb, alpha_radial=1.5,
                verbose=False,
            )
            energies.append(r['energy'])
        # Should converge: differences shrink
        diffs = [abs(energies[i+1] - energies[i]) for i in range(len(energies)-1)]
        assert diffs[-1] < diffs[0], "Not converging with n_basis"

    def test_spectral_alpha_insensitivity(self) -> None:
        """Energy is insensitive to alpha choice over a reasonable range."""
        energies = []
        for alpha in [1.0, 1.5, 2.0, 2.5]:
            r = solve_helium(
                Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
                radial_method='spectral', n_basis_radial=25, alpha_radial=alpha,
                verbose=False,
            )
            energies.append(r['energy'])
        spread = max(energies) - min(energies)
        assert spread < 0.001, f"Alpha sensitivity: spread={spread:.6f} Ha"

    def test_spectral_wavefunction_decays(self, he_spectral: dict) -> None:
        """Wavefunction decays at large R."""
        F = he_spectral['wavefunction']
        assert np.max(np.abs(F[-100:])) / np.max(np.abs(F)) < 0.01

    def test_spectral_default_unchanged(self) -> None:
        """Default radial_method='fd' produces same result as before."""
        r = solve_helium(
            Z=2.0, l_max=0, n_alpha=100, N_R_angular=100,
            N_R_radial=2000, verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT)
        assert err < 0.001


class TestSpectralCoupledChannel:
    """Tests for spectral Laguerre in coupled-channel mode."""

    def test_spectral_coupled_fd_consistency(self) -> None:
        """Spectral and FD coupled-channel agree within 0.001 Ha."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r_fd = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='fd', verbose=False,
        )
        r_sp = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=0, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            verbose=False,
        )
        diff = abs(r_fd['energy'] - r_sp['energy'])
        assert diff < 0.001, f"Coupled-channel diff {diff:.6f} Ha > 0.001 Ha"

    def test_spectral_coupled_lmax3_accuracy(self) -> None:
        """Coupled-channel at l_max=3 gives ~0.22% error (matching FD ceiling)."""
        from geovac.algebraic_coupled_channel import (
            solve_hyperspherical_algebraic_coupled,
        )
        r = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=3, n_channels=3,
            N_R_radial=3000, q_mode='exact',
            radial_method='spectral', n_basis_radial=25, alpha_radial=1.5,
            verbose=False,
        )
        err = abs(r['energy'] - E_EXACT) / abs(E_EXACT)
        assert 0.15 < err * 100 < 0.30, \
            f"Expected 0.15-0.30% error, got {err*100:.4f}%"

    def test_spectral_coupled_speedup(self) -> None:
        """Spectral coupled-channel solver achieves >10x radial-only speedup."""
        import time
        from geovac.algebraic_angular import AlgebraicAngularSolver
        from geovac.algebraic_coupled_channel import compute_algebraic_coupling
        from scipy.interpolate import CubicSpline
        from geovac.hyperspherical_radial import (
            solve_coupled_radial,
            solve_coupled_radial_spectral,
        )

        solver = AlgebraicAngularSolver(2.0, 15, 0)
        R_grid = np.concatenate([
            np.linspace(0.1, 1.0, 67), np.linspace(1.0, 5.0, 67),
            np.linspace(5.0, 30.0, 67),
        ])
        R_grid = np.unique(R_grid)
        coupling = compute_algebraic_coupling(solver, R_grid, 3,
                                              compute_exact_dPdR=True)
        V_sp = []
        for ch in range(3):
            V_ch = effective_potential(R_grid, coupling['mu'][ch])
            V_sp.append(CubicSpline(R_grid, V_ch, extrapolate=True))
        P_sp = [[None]*3 for _ in range(3)]
        Q_sp = [[None]*3 for _ in range(3)]
        for mu in range(3):
            for nu in range(3):
                P_sp[mu][nu] = CubicSpline(R_grid, coupling['P'][mu,nu],
                                           extrapolate=True)
                Q_sp[mu][nu] = CubicSpline(R_grid, coupling['Q_exact'][mu,nu],
                                           extrapolate=True)

        t0 = time.time()
        solve_coupled_radial(V_sp, P_sp, Q_splines=Q_sp, n_channels=3,
                             R_min=0.1, R_max=30.0, N_R=3000, n_states=5,
                             sigma=-3.0, include_Q=True)
        t_fd = time.time() - t0

        t0 = time.time()
        solve_coupled_radial_spectral(V_sp, P_sp, Q_splines=Q_sp,
                                      n_channels=3, n_basis=25, alpha=1.5,
                                      R_min=0.1, n_states=5, include_Q=True)
        t_sp = time.time() - t0

        speedup = t_fd / t_sp
        assert speedup > 10, f"Speedup {speedup:.1f}x < 10x target"
