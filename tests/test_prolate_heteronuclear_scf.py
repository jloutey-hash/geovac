"""
Tests for heteronuclear SCF with per-atom Z_eff optimization.

Part 2: HeH+ with independent Z_eff_A, Z_eff_B.

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_heteronuclear_scf import (
    heteronuclear_scf_energy,
    optimize_heteronuclear_zeff,
    heteronuclear_scf_ci,
)

N_SOLVE = 3000
N_GRID = 25
XI_MAX = 12.0


class TestHeteronuclearEnergy:
    """Test heteronuclear one-electron energy computation."""

    def test_physical_charges_give_finite_energy(self):
        """HeH+ with physical charges should give finite energy."""
        res = heteronuclear_scf_energy(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            Z_eff_A=2.0, Z_eff_B=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
        )
        assert not np.isnan(res['E_HF']), f"Energy is NaN"
        assert res['E_HF'] < 0, f"E_HF = {res['E_HF']:.4f}, should be negative"

    def test_inv_ra_inv_rb_positive(self):
        """<1/r_A> and <1/r_B> should both be positive."""
        res = heteronuclear_scf_energy(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            Z_eff_A=2.0, Z_eff_B=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
        )
        assert res['inv_rA'] > 0, f"<1/r_A> = {res['inv_rA']:.4f}"
        assert res['inv_rB'] > 0, f"<1/r_B> = {res['inv_rB']:.4f}"

    def test_he_dominates(self):
        """<1/r_A> (He center) should be larger than <1/r_B> (H center)."""
        res = heteronuclear_scf_energy(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            Z_eff_A=2.0, Z_eff_B=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
        )
        assert res['inv_rA'] > res['inv_rB'], \
            f"<1/r_A>={res['inv_rA']:.4f} should exceed <1/r_B>={res['inv_rB']:.4f}"

    def test_screening_lowers_j(self):
        """Screening Z_eff_A < Z_phys should reduce J_gg."""
        res_phys = heteronuclear_scf_energy(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            Z_eff_A=2.0, Z_eff_B=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
        )
        res_screen = heteronuclear_scf_energy(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            Z_eff_A=1.5, Z_eff_B=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
        )
        if np.isnan(res_phys['J_gg']) or np.isnan(res_screen['J_gg']):
            pytest.skip("Calculation failed")
        assert res_screen['J_gg'] < res_phys['J_gg'], \
            f"J(screened)={res_screen['J_gg']:.4f} should be < J(phys)={res_phys['J_gg']:.4f}"


class TestOptimization:
    """Test per-atom Z_eff optimization."""

    def test_optimization_finds_minimum(self):
        """Optimization should find finite Z_eff values."""
        opt = optimize_heteronuclear_zeff(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        assert not np.isnan(opt['E_HF_opt']), "Optimization failed"
        assert 0.5 < opt['Z_eff_A'] < 3.0, \
            f"Z_eff_A = {opt['Z_eff_A']:.3f} out of range"
        assert 0.3 < opt['Z_eff_B'] < 2.5, \
            f"Z_eff_B = {opt['Z_eff_B']:.3f} out of range"

    def test_optimization_improves_energy(self):
        """Optimized Z_eff should give lower energy than physical charges."""
        opt = optimize_heteronuclear_zeff(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        if np.isnan(opt['E_HF_opt']):
            pytest.skip("Optimization failed")
        assert opt['E_HF_opt'] <= opt['E_HF_frozen'] + 1e-4, \
            f"Opt ({opt['E_HF_opt']:.6f}) not below frozen ({opt['E_HF_frozen']:.6f})"

    def test_he_screening(self):
        """Z_eff_A (He) should be screened below physical Z=2."""
        opt = optimize_heteronuclear_zeff(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        if np.isnan(opt['Z_eff_A']):
            pytest.skip("Optimization failed")
        assert opt['Z_eff_A'] < 2.0, \
            f"Z_eff_A = {opt['Z_eff_A']:.3f}, expected screening (< 2.0)"


class TestSCFCI:
    """Test heteronuclear SCF + CI."""

    def test_scf_ci_runs(self):
        """SCF + CI should produce finite energy."""
        res = heteronuclear_scf_ci(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        assert not np.isnan(res['E_total']), "SCF+CI failed"

    def test_ci_improves_over_hf(self):
        """CI should give lower energy than HF."""
        res = heteronuclear_scf_ci(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        if np.isnan(res['E_total']):
            pytest.skip("SCF+CI failed")
        assert res['E_total'] <= res['E_HF_scf'] + 1e-6

    def test_ground_state_dominates(self):
        """Ground state should be dominated by 1sigma^2 configuration."""
        res = heteronuclear_scf_ci(
            R=1.5, Z_A_phys=2.0, Z_B_phys=1.0,
            N_xi_solve=N_SOLVE, N_grid=N_GRID, xi_max_grid=XI_MAX,
            verbose=False,
        )
        if np.isnan(res.get('c_1', np.nan)):
            pytest.skip("SCF+CI failed")
        assert abs(res['c_1']) > 0.9, \
            f"|c_1| = {abs(res['c_1']):.4f}, should dominate"
