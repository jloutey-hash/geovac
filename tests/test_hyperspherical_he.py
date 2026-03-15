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
from geovac.hyperspherical_radial import solve_radial, solve_helium


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
