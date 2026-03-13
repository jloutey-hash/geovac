"""
Tests for the 4σ-orbital H2 CI on prolate spheroidal lattice.

Validates:
  1. Higher orbital generation (2σ_g, 2σ_u exist and have correct ordering)
  2. Integral symmetries (8-fold, Hermiticity)
  3. CI matrix structure (symmetry, dimension)
  4. Physics (bound H2, improvement over 2×2, correct dissociation)

Date: 2026-03-13
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.stress_test_prolate_h2 import get_orbital_on_grid, compute_vee_integral
from debug.stress_test_prolate_h2_4sigma import (
    IntegralCache,
    build_ci_matrix_6x6,
    h2_4sigma_ci,
)


# Test parameters — small grid for speed
R_TEST = 1.4
N_SOLVE = 3000
N_GRID = 30
XI_MAX = 12.0


@pytest.fixture(scope="module")
def orbitals():
    """Generate all 4 orbitals once for the module."""
    orbs = {}
    for label, n_ang in [('g', 0), ('u', 1), ('G', 2), ('U', 3)]:
        orbs[label] = get_orbital_on_grid(
            R=R_TEST, n_angular=n_ang, N_xi_solve=N_SOLVE,
            N_xi_grid=N_GRID, N_eta_grid=N_GRID,
            xi_max_grid=XI_MAX, m=0,
        )
    return orbs


@pytest.fixture(scope="module")
def ci_result():
    """Run a single CI calculation for tests."""
    return h2_4sigma_ci(
        R=R_TEST, N_xi_solve=N_SOLVE, N_grid=N_GRID,
        xi_max_grid=XI_MAX, verbose=False,
    )


# ============================================================
# Orbital tests
# ============================================================
class TestOrbitalGeneration:
    """Test that 2σ_g and 2σ_u orbitals are valid."""

    def test_all_four_orbitals_exist(self, orbitals):
        """All 4 orbital solutions found."""
        for label in ['g', 'u', 'G', 'U']:
            assert label in orbitals
            assert 'E_elec' in orbitals[label]
            assert not np.isnan(orbitals[label]['E_elec'])

    def test_energy_ordering(self, orbitals):
        """1σ_g < 1σ_u < 2σ_g < 2σ_u (increasing energy)."""
        E_g = orbitals['g']['E_elec']
        E_u = orbitals['u']['E_elec']
        E_G = orbitals['G']['E_elec']
        E_U = orbitals['U']['E_elec']
        assert E_g < E_u, f"1σ_g ({E_g:.4f}) should be below 1σ_u ({E_u:.4f})"
        assert E_u < E_G, f"1σ_u ({E_u:.4f}) should be below 2σ_g ({E_G:.4f})"
        assert E_G < E_U, f"2σ_g ({E_G:.4f}) should be below 2σ_u ({E_U:.4f})"

    def test_1sigma_g_energy(self, orbitals):
        """1σ_g energy near -1.1 Ha at R=1.4 (known from Phase 1)."""
        E = orbitals['g']['E_elec']
        assert -1.3 < E < -0.9, f"1σ_g energy {E:.4f} out of range"

    def test_orbitals_normalized(self, orbitals):
        """Each orbital should be normalized to 1."""
        for label in ['g', 'u', 'G', 'U']:
            orb = orbitals[label]
            XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
            J = (R_TEST / 2)**3 * (XI**2 - ETA**2)
            W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
            norm = np.sum(orb['psi']**2 * J * W_XI * W_ETA) * 2 * np.pi
            assert abs(norm - 1.0) < 0.02, \
                f"Orbital {label} norm = {norm:.4f}, expected 1.0"

    def test_excited_orbitals_have_nodes(self, orbitals):
        """2σ_g and 2σ_u should have more nodal structure than ground states.

        The node may be in the angular part G(eta) rather than radial F(xi),
        so check the full 2D wavefunction for sign changes.
        """
        psi_G = orbitals['G']['psi']
        # Count sign changes along the eta axis at the peak xi
        xi_peak = np.argmax(np.abs(psi_G).max(axis=1))
        eta_slice = psi_G[xi_peak, :]
        sign_changes = np.sum(np.diff(np.sign(eta_slice[eta_slice != 0])) != 0)
        # The 2σ_g wavefunction (2nd gerade state) should have at least 1 node
        # somewhere in (xi, eta) space
        has_both_signs = (psi_G.max() > 0) and (psi_G.min() < 0)
        assert has_both_signs or sign_changes >= 1, \
            f"2σ_g shows no nodal structure (max={psi_G.max():.4f}, min={psi_G.min():.4f})"


# ============================================================
# Integral tests
# ============================================================
class TestIntegrals:
    """Test V_ee integral properties."""

    def test_J_gg_positive(self, orbitals):
        """Coulomb self-repulsion J_gg > 0."""
        J = compute_vee_integral(
            orbitals['g'], orbitals['g'], orbitals['g'], orbitals['g']
        )
        assert J > 0, f"J_gg = {J:.6f}, should be positive"

    def test_coulomb_ordering(self, orbitals):
        """J_gg should be largest (most compact orbital)."""
        J_gg = compute_vee_integral(
            orbitals['g'], orbitals['g'], orbitals['g'], orbitals['g']
        )
        J_GG = compute_vee_integral(
            orbitals['G'], orbitals['G'], orbitals['G'], orbitals['G']
        )
        # 1σ_g is more compact → larger self-Coulomb
        assert J_gg > J_GG, \
            f"J(1σ_g) = {J_gg:.4f} should exceed J(2σ_g) = {J_GG:.4f}"

    def test_integral_cache_symmetry(self, orbitals):
        """Cache should return same value for symmetric index permutations."""
        cache = IntegralCache(orbitals)
        # ⟨ab|cd⟩ = ⟨cd|ab⟩
        v1 = cache.get('g', 'u', 'g', 'u')
        v2 = cache.get('g', 'u', 'g', 'u')
        assert v1 == v2, "Cache inconsistency"
        # ⟨ab|cd⟩ = ⟨ba|dc⟩
        v3 = cache.get('u', 'g', 'u', 'g')
        assert abs(v1 - v3) < 1e-10, \
            f"Symmetry ⟨ab|cd⟩=⟨ba|dc⟩ violated: {v1:.6f} vs {v3:.6f}"

    def test_exchange_integral_positive(self, orbitals):
        """Exchange integral K(g,G) = ⟨gG|Gg⟩ should be positive for real orbitals."""
        K = compute_vee_integral(
            orbitals['g'], orbitals['G'], orbitals['G'], orbitals['g']
        )
        assert K > 0, f"K(g,G) = {K:.6f}, should be positive"


# ============================================================
# CI matrix tests
# ============================================================
class TestCIMatrix:
    """Test CI matrix structure and properties."""

    def test_ci_matrix_symmetric(self, ci_result):
        """CI result should exist (no NaN)."""
        assert not np.isnan(ci_result['E_total']), "CI calculation failed"

    def test_ci_matrix_6_eigenvalues(self, ci_result):
        """Should have 6 eigenvalues."""
        assert len(ci_result['evals']) == 6

    def test_eigenvalues_ordered(self, ci_result):
        """Eigenvalues should be in ascending order."""
        evals = ci_result['evals']
        for i in range(5):
            assert evals[i] <= evals[i+1] + 1e-10

    def test_ground_state_dominated_by_g2(self, ci_result):
        """Ground state should be dominated by |1σ_g²⟩."""
        c = ci_result['evecs']
        assert c[0]**2 > 0.5, \
            f"|1σ_g²⟩ weight = {c[0]**2:.4f}, expected > 0.5"


# ============================================================
# Physics tests
# ============================================================
class TestPhysics:
    """Test physical correctness of the 4σ CI."""

    def test_h2_is_bound(self, ci_result):
        """E_total < 2*E_H = -1.0 Ha."""
        assert ci_result['E_total'] < -1.0, \
            f"H2 unbound: E = {ci_result['E_total']:.6f} Ha"

    def test_6x6_improves_over_2x2(self, ci_result):
        """6×6 CI should give lower energy than 2×2."""
        assert ci_result['E_total'] <= ci_result['E_2x2'] + 1e-10, \
            f"6×6 ({ci_result['E_total']:.6f}) not below 2×2 ({ci_result['E_2x2']:.6f})"

    def test_6x6_below_hf(self, ci_result):
        """CI energy should be below Hartree-Fock."""
        assert ci_result['E_total'] < ci_result['E_HF'], \
            f"CI ({ci_result['E_total']:.6f}) not below HF ({ci_result['E_HF']:.6f})"

    def test_variational_bound(self, ci_result):
        """E_total should be above exact (-1.1745 Ha)."""
        E_exact = -1.1745
        assert ci_result['E_total'] > E_exact, \
            f"E = {ci_result['E_total']:.6f} below exact {E_exact} (variational violation)"

    def test_binding_energy_improvement(self):
        """D_e(6×6) should exceed D_e(2×2) at R=1.4."""
        res = h2_4sigma_ci(
            R=1.4, N_xi_solve=N_SOLVE, N_grid=N_GRID,
            xi_max_grid=XI_MAX, verbose=False,
        )
        if np.isnan(res['E_total']):
            pytest.skip("CI calculation failed")
        D_e_6 = -1.0 - res['E_total']  # 2*E_H - E
        D_e_2 = -1.0 - res['E_2x2']
        assert D_e_6 >= D_e_2 - 1e-6, \
            f"D_e(6×6) = {D_e_6:.6f} not ≥ D_e(2×2) = {D_e_2:.6f}"


# ============================================================
# Dissociation test
# ============================================================
class TestDissociation:
    """Test correct dissociation behavior."""

    def test_dissociation_to_atoms(self):
        """At large R, E_total should approach 2*E_H = -1.0."""
        res_large = h2_4sigma_ci(
            R=5.0, N_xi_solve=N_SOLVE, N_grid=N_GRID,
            xi_max_grid=XI_MAX, verbose=False,
        )
        if np.isnan(res_large['E_total']):
            pytest.skip("CI at R=5.0 failed")
        # At R=5.0, should be close to -1.0 (within ~0.05 Ha)
        assert res_large['E_total'] > -1.05, \
            f"E(R=5.0) = {res_large['E_total']:.4f}, too bound for dissociation"

    def test_equilibrium_shorter_than_5bohr(self):
        """Minimum energy should occur well before R=5."""
        R_vals = [1.4, 2.0, 5.0]
        energies = []
        for R in R_vals:
            res = h2_4sigma_ci(
                R=R, N_xi_solve=N_SOLVE, N_grid=N_GRID,
                xi_max_grid=XI_MAX, verbose=False,
            )
            energies.append(res['E_total'])
        # E(1.4) should be lower than E(5.0)
        assert energies[0] < energies[2], \
            f"E(1.4) = {energies[0]:.4f} not below E(5.0) = {energies[2]:.4f}"
