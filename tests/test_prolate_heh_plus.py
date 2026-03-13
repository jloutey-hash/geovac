"""
Tests for HeH+ two-electron CI on prolate spheroidal lattice.

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.stress_test_prolate_heh_plus import (
    get_orbital_on_grid_general,
    heh_plus_ci,
)

R_TEST = 1.5
N_SOLVE = 3000
N_GRID = 30
XI_MAX = 12.0


@pytest.fixture(scope="module")
def ci_result():
    """Run HeH+ CI once for the module."""
    return heh_plus_ci(
        R=R_TEST, N_xi_solve=N_SOLVE, N_grid=N_GRID,
        xi_max_grid=XI_MAX, verbose=False,
    )


class TestHeteronuclearOrbitals:
    """Test orbital generation for Z_A=2, Z_B=1."""

    def test_orbital_exists(self):
        """HeH2+ ground state orbital solves."""
        orb = get_orbital_on_grid_general(
            R=R_TEST, Z_A=2, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        assert not np.isnan(orb['E_elec'])

    def test_ground_below_excited(self):
        """1sigma energy < 2sigma energy."""
        orb1 = get_orbital_on_grid_general(
            R=R_TEST, Z_A=2, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        orb2 = get_orbital_on_grid_general(
            R=R_TEST, Z_A=2, Z_B=1, n_angular=1,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        assert orb1['E_elec'] < orb2['E_elec']

    def test_orbital_normalized(self):
        """Orbital should be normalized."""
        orb = get_orbital_on_grid_general(
            R=R_TEST, Z_A=2, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J = (R_TEST / 2)**3 * (XI**2 - ETA**2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm = np.sum(orb['psi']**2 * J * W_XI * W_ETA) * 2 * np.pi
        assert abs(norm - 1.0) < 0.02


class TestHeHPlusCI:
    """Test HeH+ CI physics."""

    def test_ci_succeeds(self, ci_result):
        """CI calculation should not fail."""
        assert not np.isnan(ci_result['E_total'])

    def test_heh_energy_reasonable(self, ci_result):
        """E_total should be in a physically reasonable range.

        Note: frozen HeH2+ orbitals lack screening, so HeH+ may not be bound
        relative to He + H+ = -2.904 Ha. The electron repulsion (J~1.3 Ha)
        is too large with unrelaxed orbitals. SCF Z_eff optimization needed
        for binding.
        """
        # Should at least be below -2.5 (not wildly wrong)
        assert ci_result['E_total'] < -2.5, \
            f"HeH+ energy too high: E={ci_result['E_total']:.4f}"

    def test_ci_below_hf(self, ci_result):
        """CI should be below HF."""
        assert ci_result['E_total'] < ci_result['E_HF'] + 1e-6

    def test_reasonable_energy(self, ci_result):
        """Energy should be in a reasonable range."""
        assert -4.0 < ci_result['E_total'] < -2.5

    def test_j11_positive(self, ci_result):
        """Coulomb integral should be positive."""
        assert ci_result['J_11'] > 0

    def test_ground_state_dominated(self, ci_result):
        """Ground state should be >50% bonding orbital."""
        assert ci_result['c_1']**2 > 0.5

    def test_equilibrium_exists(self):
        """PES should have a minimum between R=1 and R=3."""
        energies = []
        for R in [1.0, 1.5, 3.0]:
            res = heh_plus_ci(
                R=R, N_xi_solve=N_SOLVE, N_grid=N_GRID,
                xi_max_grid=XI_MAX, verbose=False,
            )
            energies.append(res['E_total'])
        # E(1.5) should be lower than both E(1.0) and E(3.0)
        assert energies[1] < energies[0] and energies[1] < energies[2], \
            f"No minimum: E(1.0)={energies[0]:.4f}, E(1.5)={energies[1]:.4f}, E(3.0)={energies[2]:.4f}"
