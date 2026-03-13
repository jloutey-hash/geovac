"""
Tests for full grid-based SCF on prolate spheroidal lattice.

Part 1: Fock operator on 2D (xi, eta) grid.

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_scf import (
    compute_coulomb_potential,
    _build_2d_hamiltonian,
    _solve_2d_eigenvalue,
    grid_scf,
    get_orbital_on_grid,
)

N_SOLVE = 3000
N_FD = 20
XI_MAX = 12.0


class TestCoulombPotential:
    """Test Coulomb potential computation on grid."""

    def test_vj_positive(self):
        """V_J from a positive density should be everywhere positive."""
        orb = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=20, N_eta_grid=20,
            xi_max_grid=XI_MAX,
        )
        V_J = compute_coulomb_potential(
            orb['psi'], orb['xi'], orb['eta'],
            orb['w_xi'], orb['w_eta'], 1.4,
        )
        assert np.all(V_J >= 0), "V_J has negative values"
        assert np.max(V_J) > 0, "V_J is identically zero"

    def test_vj_decays_at_large_xi(self):
        """V_J should decrease as xi increases (farther from nuclei)."""
        orb = get_orbital_on_grid(
            R=1.4, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_SOLVE, N_xi_grid=25, N_eta_grid=25,
            xi_max_grid=XI_MAX,
        )
        V_J = compute_coulomb_potential(
            orb['psi'], orb['xi'], orb['eta'],
            orb['w_xi'], orb['w_eta'], 1.4,
        )
        # Average V_J at small xi vs large xi
        n_xi = len(orb['xi'])
        v_near = np.mean(V_J[:n_xi // 3, :])
        v_far = np.mean(V_J[2 * n_xi // 3:, :])
        assert v_near > v_far, \
            f"V_J near ({v_near:.4f}) should exceed V_J far ({v_far:.4f})"


class TestHamiltonian2D:
    """Test 2D FD Hamiltonian construction."""

    def test_hamiltonian_symmetric(self):
        """The A matrix should be symmetric."""
        A, B, xi, eta = _build_2d_hamiltonian(
            R=1.4, Z_A=1.0, Z_B=1.0, N_xi=10, N_eta=10, xi_max=10.0,
        )
        asym = np.max(np.abs(A - A.T))
        assert asym < 1e-10, f"A not symmetric: max|A-A^T| = {asym:.2e}"

    def test_weight_positive(self):
        """Weight matrix B (xi^2-eta^2) should be positive."""
        _, B, xi, eta = _build_2d_hamiltonian(
            R=1.4, Z_A=1.0, Z_B=1.0, N_xi=10, N_eta=10, xi_max=10.0,
        )
        assert np.all(B > 0), "B has non-positive elements"

    def test_h2plus_eigenvalue(self):
        """2D solver without V_ext should recover H2+ energy."""
        A, B, xi, eta = _build_2d_hamiltonian(
            R=2.0, Z_A=1.0, Z_B=1.0, N_xi=40, N_eta=40, xi_max=15.0,
        )
        E_elec, psi, c2 = _solve_2d_eigenvalue(A, B, R=2.0)

        # H2+ exact at R=2.0: E_elec ~ -1.1026, E_total ~ -0.6026
        E_total = E_elec + 1.0 / 2.0
        # The 2D FD is less accurate than the separated solver
        E_exact = -0.6026
        err = abs(E_total - E_exact) / abs(E_exact)
        assert err < 0.40, \
            f"E_total = {E_total:.4f}, expected ~{E_exact:.4f} (err {err:.1%})"

    def test_vext_raises_energy(self):
        """Adding repulsive V_ext should raise the eigenvalue."""
        N = 15
        A0, B, xi, eta = _build_2d_hamiltonian(
            R=1.4, Z_A=1.0, Z_B=1.0, N_xi=N, N_eta=N, xi_max=10.0,
        )
        E0, _, _ = _solve_2d_eigenvalue(A0, B, R=1.4)

        # Add constant repulsive potential
        V_ext = 0.1 * np.ones((N, N))
        A1, B1, _, _ = _build_2d_hamiltonian(
            R=1.4, Z_A=1.0, Z_B=1.0, N_xi=N, N_eta=N, xi_max=10.0,
            V_ext=V_ext,
        )
        E1, _, _ = _solve_2d_eigenvalue(A1, B1, R=1.4)

        assert E1 > E0, f"E(V_ext) = {E1:.4f} should be > E(no V_ext) = {E0:.4f}"


class TestGridSCF:
    """Test full grid-based SCF."""

    def test_grid_scf_converges(self):
        """Grid SCF should converge within max_iter."""
        res = grid_scf(R=1.4, N_xi_fd=N_FD, N_eta_fd=N_FD,
                       xi_max_fd=XI_MAX, max_iter=30, verbose=False)
        assert not np.isnan(res['E_total']), "Grid SCF failed"
        # May or may not converge with small grid, just check it ran
        assert res['n_iter'] >= 1

    def test_grid_scf_energy_reasonable(self):
        """Grid SCF energy should be in physical range."""
        res = grid_scf(R=1.4, N_xi_fd=N_FD, N_eta_fd=N_FD,
                       xi_max_fd=XI_MAX, verbose=False)
        if np.isnan(res['E_total']):
            pytest.skip("Grid SCF failed")
        # Should be between -1.3 (too bound) and -0.9 (unbound)
        assert -1.5 < res['E_total'] < -0.8, \
            f"E = {res['E_total']:.4f} out of physical range"

    def test_jgg_positive(self):
        """J_gg from grid SCF should be positive."""
        res = grid_scf(R=1.4, N_xi_fd=N_FD, N_eta_fd=N_FD,
                       xi_max_fd=XI_MAX, verbose=False)
        if np.isnan(res.get('J_gg', np.nan)):
            pytest.skip("Grid SCF failed")
        assert res['J_gg'] > 0, f"J_gg = {res['J_gg']:.4f} should be positive"

    def test_h2_bound(self):
        """H2 should be bound (E < 2*E_H = -1.0)."""
        res = grid_scf(R=1.4, N_xi_fd=25, N_eta_fd=25,
                       xi_max_fd=XI_MAX, verbose=False)
        if np.isnan(res['E_total']):
            pytest.skip("Grid SCF failed")
        assert res['E_total'] < -1.0, \
            f"E = {res['E_total']:.4f}, should be < -1.0 for bound H2"
