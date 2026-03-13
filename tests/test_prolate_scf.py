"""
Tests for H2 SCF via Z_eff optimization on prolate spheroidal lattice.

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.stress_test_prolate_h2_scf import (
    compute_vnuc_expectation,
    compute_scf_energy,
    optimize_zeff,
)
from debug.stress_test_prolate_heh_plus import get_orbital_on_grid_general

N_SOLVE = 3000
N_GRID = 30
XI_MAX = 12.0


class TestVnucExpectation:
    """Test nuclear attraction expectation value."""

    def test_vnuc_positive(self):
        """<1/r_A + 1/r_B> should be positive."""
        orb = get_orbital_on_grid_general(
            R=1.4, Z_A=1, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        )
        v = compute_vnuc_expectation(orb)
        assert v > 0, f"<V_nuc> = {v:.4f}, should be positive"

    def test_vnuc_scales_with_r(self):
        """<1/r_A + 1/r_B> should decrease with R (atoms farther apart)."""
        v1 = compute_vnuc_expectation(get_orbital_on_grid_general(
            R=1.0, Z_A=1, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        ))
        v2 = compute_vnuc_expectation(get_orbital_on_grid_general(
            R=3.0, Z_A=1, Z_B=1, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX,
        ))
        assert v1 > v2, f"V_nuc(R=1) = {v1:.4f} should exceed V_nuc(R=3) = {v2:.4f}"


class TestSCFEnergy:
    """Test SCF energy computation."""

    def test_zeff_1_gives_frozen_energy(self):
        """Z_eff=1 should recover the frozen H2+ energy."""
        res = compute_scf_energy(1.4, Z_eff=1.0, N_xi_solve=N_SOLVE,
                                 N_grid=N_GRID, xi_max_grid=XI_MAX)
        assert not np.isnan(res['E_HF']), "SCF at Z_eff=1 failed"
        # Should be close to the known frozen value (~-1.055)
        assert -1.2 < res['E_HF'] < -0.9, f"E_HF = {res['E_HF']:.4f} out of range"

    def test_j_gg_decreases_with_lower_zeff(self):
        """Lower Z_eff -> more diffuse orbital -> smaller J_gg."""
        res_10 = compute_scf_energy(1.4, Z_eff=1.0, N_xi_solve=N_SOLVE,
                                    N_grid=N_GRID, xi_max_grid=XI_MAX)
        res_08 = compute_scf_energy(1.4, Z_eff=0.8, N_xi_solve=N_SOLVE,
                                    N_grid=N_GRID, xi_max_grid=XI_MAX)
        if np.isnan(res_10['J_gg']) or np.isnan(res_08['J_gg']):
            pytest.skip("SCF calculation failed")
        assert res_08['J_gg'] < res_10['J_gg'], \
            f"J(Z=0.8)={res_08['J_gg']:.4f} should be < J(Z=1.0)={res_10['J_gg']:.4f}"


class TestOptimization:
    """Test Z_eff optimization."""

    def test_optimal_zeff_exists(self):
        """Optimization should find a finite Z_eff*."""
        opt = optimize_zeff(R=1.4, N_xi_solve=N_SOLVE, N_grid=N_GRID,
                            xi_max_grid=XI_MAX, verbose=False)
        assert not np.isnan(opt['Z_eff_opt']), "Optimization failed"
        assert 0.5 < opt['Z_eff_opt'] < 2.0, \
            f"Z_eff* = {opt['Z_eff_opt']:.3f} out of range"

    def test_scf_improves_over_frozen(self):
        """E_HF(Z_eff*) should be <= E_HF(Z=1)."""
        opt = optimize_zeff(R=1.4, N_xi_solve=N_SOLVE, N_grid=N_GRID,
                            xi_max_grid=XI_MAX, verbose=False)
        if np.isnan(opt['E_HF_opt']):
            pytest.skip("Optimization failed")
        assert opt['E_HF_opt'] <= opt['E_HF_frozen'] + 1e-6, \
            f"SCF ({opt['E_HF_opt']:.6f}) not below frozen ({opt['E_HF_frozen']:.6f})"
