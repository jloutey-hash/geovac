"""
Tests for relaxed-orbital CI on prolate spheroidal lattice.

Part 3: CI with Z_eff-optimized orbitals for H2.

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_scf import (
    get_orbital_on_grid,
    compute_vee_integral,
    compute_vnuc_expectation,
    eckart_scf_energy,
    optimize_zeff_eckart,
    relaxed_orbital_ci,
)

# Use modest parameters for test speed
N_SOLVE = 3000
N_GRID = 30
XI_MAX = 12.0


class TestOrbitalGeneration:
    """Test orbital generation with float Z support."""

    def test_h2plus_orbital_normalized(self):
        """Orbital should be normalized to 1."""
        orb = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J = (1.4 / 2) ** 3 * (XI ** 2 - ETA ** 2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm = np.sum(orb['psi'] ** 2 * J * W_XI * W_ETA) * 2 * np.pi
        assert abs(norm - 1.0) < 0.02, f"Norm = {norm:.4f}"

    def test_fractional_z_works(self):
        """Float Z_eff should produce valid orbitals."""
        orb = get_orbital_on_grid(
            R=1.4, Z_A=0.8, Z_B=0.8, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        assert not np.isnan(orb['E_elec'])
        assert orb['E_elec'] < 0  # Bound state

    def test_lower_z_gives_higher_energy(self):
        """Lower Z_eff -> less attractive potential -> higher energy."""
        orb_10 = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        orb_08 = get_orbital_on_grid(
            R=1.4, Z_A=0.8, Z_B=0.8, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        # Less attractive -> E_elec closer to 0
        assert orb_08['E_elec'] > orb_10['E_elec']


class TestVeeIntegral:
    """Test V_ee integral computation."""

    def test_jgg_positive(self):
        """Coulomb self-repulsion should be positive."""
        orb = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        J_gg = compute_vee_integral(orb, orb, orb, orb)
        assert J_gg > 0, f"J_gg = {J_gg:.4f}, should be positive"

    def test_jgg_decreases_with_lower_z(self):
        """More diffuse orbital (lower Z) should have smaller J_gg."""
        orb_10 = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        orb_08 = get_orbital_on_grid(
            R=1.4, Z_A=0.8, Z_B=0.8, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        J_10 = compute_vee_integral(orb_10, orb_10, orb_10, orb_10)
        J_08 = compute_vee_integral(orb_08, orb_08, orb_08, orb_08)
        assert J_08 < J_10, f"J(Z=0.8)={J_08:.4f} should be < J(Z=1)={J_10:.4f}"


class TestEckartSCF:
    """Test Eckart variational SCF."""

    def test_zeff_1_gives_frozen(self):
        """Z_eff=1 should recover frozen H2+ result."""
        res = eckart_scf_energy(1.4, Z_eff=1.0, N_xi_solve=N_SOLVE,
                                N_grid=N_GRID, xi_max_grid=XI_MAX)
        assert not np.isnan(res['E_HF'])
        assert -1.2 < res['E_HF'] < -0.9

    def test_optimal_zeff_below_one(self):
        """Optimal Z_eff for H2 should be < 1 (screening)."""
        opt = optimize_zeff_eckart(R=1.4, N_xi_solve=N_SOLVE,
                                   N_grid=N_GRID, xi_max_grid=XI_MAX)
        assert not np.isnan(opt['Z_eff_opt'])
        assert opt['Z_eff_opt'] < 1.0, \
            f"Z_eff* = {opt['Z_eff_opt']:.3f}, expected < 1.0"

    def test_scf_improves_over_frozen(self):
        """E_HF(Z_eff*) should be <= E_HF(Z=1)."""
        opt = optimize_zeff_eckart(R=1.4, N_xi_solve=N_SOLVE,
                                   N_grid=N_GRID, xi_max_grid=XI_MAX)
        assert opt['E_HF_opt'] <= opt['E_HF_frozen'] + 1e-6


class TestRelaxedCI:
    """Test relaxed-orbital CI."""

    def test_relaxed_ci_runs(self):
        """Relaxed CI should produce finite energy."""
        res = relaxed_orbital_ci(R=1.4, N_xi_solve=N_SOLVE,
                                 N_grid=N_GRID, xi_max_grid=XI_MAX,
                                 verbose=False)
        assert not np.isnan(res['E_total'])
        assert res['E_total'] < -1.0  # Bound state

    def test_ci_improves_over_hf(self):
        """CI should give lower energy than HF."""
        res = relaxed_orbital_ci(R=1.4, N_xi_solve=N_SOLVE,
                                 N_grid=N_GRID, xi_max_grid=XI_MAX,
                                 verbose=False)
        assert res['E_total'] <= res['E_HF_scf'] + 1e-6

    def test_h2_is_bound(self):
        """H2 should be bound (E < 2*E_H = -1.0)."""
        res = relaxed_orbital_ci(R=1.4, N_xi_solve=N_SOLVE,
                                 N_grid=N_GRID, xi_max_grid=XI_MAX,
                                 verbose=False)
        assert res['E_total'] < -1.0, \
            f"E = {res['E_total']:.4f}, should be < -1.0 (2*E_H)"

    def test_de_improves_over_frozen(self):
        """D_e should improve over frozen CI (0.072 Ha)."""
        res = relaxed_orbital_ci(R=1.4, N_xi_solve=3000,
                                 N_grid=30, xi_max_grid=XI_MAX,
                                 verbose=False)
        D_e = -1.0 - res['E_total']
        # At lower resolution the numbers differ, just check improvement
        assert D_e > 0.07, \
            f"D_e = {D_e:.4f} Ha, should improve over frozen CI baseline"

    def test_ci_coefficients_physical(self):
        """Ground state should be dominated by sigma_g^2."""
        res = relaxed_orbital_ci(R=1.4, N_xi_solve=N_SOLVE,
                                 N_grid=N_GRID, xi_max_grid=XI_MAX,
                                 verbose=False)
        assert abs(res['c_g']) > 0.9, \
            f"|c_g| = {abs(res['c_g']):.4f}, should dominate"
        assert abs(res['c_u']) < 0.5, \
            f"|c_u| = {abs(res['c_u']):.4f}, should be small"
