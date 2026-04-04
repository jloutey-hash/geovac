"""
Tests for the Sturmian CI solver (Track BU-1).

Validates:
  1. Hydrogenic radial function normalization
  2. Sturmian overlap matrix properties
  3. One-body Hamiltonian consistency
  4. He ground state energy vs standard FCI
  5. Variational bound
  6. Gaunt selection rule sparsity preservation
"""

import pytest
import numpy as np
from scipy.integrate import quad

from geovac.sturmian_solver import (
    hydrogenic_radial, _radial_overlap, _slater_rk, _slater_rk_fast,
    SturmianCI, StandardFCI, GeneralizedSturmianCI,
    compare_he, compare_he_generalized,
)


# ---------------------------------------------------------------------------
# Radial function tests
# ---------------------------------------------------------------------------

class TestHydrogenicRadial:
    """Test hydrogenic radial wavefunctions with variable Z_eff."""

    def test_1s_normalization(self):
        """1s orbital with Z_eff=2 is normalized."""
        val, _ = quad(lambda r: hydrogenic_radial(r, 1, 0, 2.0) ** 2 * r ** 2,
                      0, 30, limit=200)
        assert abs(val - 1.0) < 1e-10

    def test_2s_normalization(self):
        """2s orbital with Z_eff=4 is normalized."""
        val, _ = quad(lambda r: hydrogenic_radial(r, 2, 0, 4.0) ** 2 * r ** 2,
                      0, 30, limit=200)
        assert abs(val - 1.0) < 1e-10

    def test_2p_normalization(self):
        """2p orbital with Z_eff=4 is normalized."""
        val, _ = quad(lambda r: hydrogenic_radial(r, 2, 1, 4.0) ** 2 * r ** 2,
                      0, 30, limit=200)
        assert abs(val - 1.0) < 1e-10

    def test_same_zeff_orthogonality(self):
        """1s and 2s with same Z_eff are orthogonal."""
        Z = 2.0
        val, _ = quad(
            lambda r: hydrogenic_radial(r, 1, 0, Z) * hydrogenic_radial(r, 2, 0, Z) * r ** 2,
            0, 50, limit=300
        )
        assert abs(val) < 1e-10

    def test_sturmian_zeff_nonorthogonal(self):
        """1s(Z1=k) and 2s(Z2=2k) Sturmians are NOT orthogonal."""
        k = 1.5
        val = _radial_overlap(1, 0, 1 * k, 2, 0, 2 * k)
        # Should be nonzero for Sturmians with different Z_eff
        assert abs(val) > 0.01, f"Sturmian overlap too small: {val}"


# ---------------------------------------------------------------------------
# Sturmian overlap and h1 tests
# ---------------------------------------------------------------------------

class TestSturmianMatrices:
    """Test one-particle matrices in the Sturmian basis."""

    def test_overlap_diagonal(self):
        """Overlap diagonal should be 1 (normalized orbitals)."""
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        S = solver._build_overlap(k=1.5)
        np.testing.assert_allclose(np.diag(S), 1.0, atol=1e-8)

    def test_overlap_symmetry(self):
        """Overlap matrix should be symmetric."""
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        S = solver._build_overlap(k=1.5)
        np.testing.assert_allclose(S, S.T, atol=1e-12)

    def test_overlap_positive_definite(self):
        """Overlap matrix should be positive definite."""
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        S = solver._build_overlap(k=1.5)
        eigvals = np.linalg.eigvalsh(S)
        assert np.all(eigvals > 0), f"Non-positive eigenvalue: {eigvals.min()}"

    def test_h1_hydrogen_1s(self):
        """For H (Z=1), 1s energy should be -k²/2 + k*(1/1 - 1) = -k²/2."""
        # Actually: h1_{1s,1s} = k²/2 - Z*k/n = k²/2 - k
        # For Z=1, n=1: h1 = k²/2 - k
        # This is -1/2 at k=1 (correct for H ground state)
        solver = SturmianCI(Z=1, n_electrons=1, max_n=1)
        S = solver._build_overlap(k=1.0)
        h1 = solver._build_h1_sturmian(k=1.0, S=S)
        expected = 0.5 - 1.0  # k²/2 - Z*k/n = 0.5 - 1.0 = -0.5
        assert abs(h1[0, 0] - expected) < 1e-12

    def test_h1_symmetry(self):
        """h1 should be symmetric."""
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        S = solver._build_overlap(k=1.5)
        h1 = solver._build_h1_sturmian(k=1.5, S=S)
        np.testing.assert_allclose(h1, h1.T, atol=1e-12)


# ---------------------------------------------------------------------------
# FCI comparison tests
# ---------------------------------------------------------------------------

class TestHeliumComparison:
    """Compare Sturmian CI and standard FCI for He."""

    @pytest.mark.slow
    def test_he_max_n2_variational_bound(self):
        """Sturmian CI for He at max_n=2 respects variational bound."""
        E_exact = -2.903724
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = solver.optimize_k(k_range=(0.5, 4.0), n_scan=15)
        assert result['energy'] > E_exact, (
            f"Variational bound violated: {result['energy']:.6f} < {E_exact}"
        )

    @pytest.mark.slow
    def test_he_max_n2_comparison(self):
        """Compare Sturmian and standard FCI for He at max_n=2."""
        result = compare_he(max_n=2, verbose=True)

        # Both should give reasonable energies
        assert result['err_standard'] < 20.0, f"Standard FCI error too large: {result['err_standard']:.2f}%"
        assert result['err_sturmian'] < 20.0, f"Sturmian CI error too large: {result['err_sturmian']:.2f}%"

        # Report which is better
        print(f"\n  Standard FCI: {result['err_standard']:.4f}%")
        print(f"  Sturmian CI:  {result['err_sturmian']:.4f}%")

    @pytest.mark.slow
    def test_he_max_n3_comparison(self):
        """Compare Sturmian and standard FCI for He at max_n=3."""
        result = compare_he(max_n=3, verbose=True)

        assert result['err_standard'] < 10.0
        assert result['err_sturmian'] < 20.0

    @pytest.mark.slow
    def test_standard_fci_hydrogen(self):
        """Standard FCI for H (1 electron) should give exact energy."""
        std = StandardFCI(Z=1, n_electrons=1, max_n=2)
        result = std.solve(Z_eff=1.0)
        # H ground state: -0.5 Ha
        assert abs(result['energy'] - (-0.5)) < 1e-4

    @pytest.mark.slow
    def test_sturmian_hydrogen(self):
        """Sturmian CI for H (1 electron) should give exact energy."""
        solver = SturmianCI(Z=1, n_electrons=1, max_n=2)
        result = solver.optimize_k(k_range=(0.5, 2.0), n_scan=10)
        assert abs(result['energy'] - (-0.5)) < 1e-3

    @pytest.mark.slow
    def test_gaunt_sparsity_preserved(self):
        """V_ee sparsity pattern should match between Sturmian and standard."""
        # The Gaunt selection rules are angular — they should be identical
        # in both bases. Only the radial integrals differ.
        solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = solver.solve(k=1.5)

        # At max_n=2 with s+p orbitals, most ERIs are zero due to
        # m-selection rule and parity. Check that we get a reasonable
        # nonzero fraction.
        assert result['eri_count'] > 0, "No nonzero ERIs"
        assert result['vee_sparsity'] < 0.5, "V_ee not sparse enough"


# ---------------------------------------------------------------------------
# Quick smoke test (not marked slow)
# ---------------------------------------------------------------------------

def test_sturmian_solver_import():
    """Module imports without error."""
    from geovac.sturmian_solver import SturmianCI, StandardFCI, compare_he
    solver = SturmianCI(Z=2, n_electrons=2, max_n=2)
    assert solver.n_spatial == 5
    assert solver.n_sd == 45


def test_generalized_sturmian_import():
    """GeneralizedSturmianCI imports and initializes."""
    solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
    assert solver.n_spatial == 5
    assert solver.n_sd == 45
    assert solver._eri_angular_count > 0


# ---------------------------------------------------------------------------
# Fast R^k consistency test
# ---------------------------------------------------------------------------

class TestFastRk:
    """Verify _slater_rk_fast matches _slater_rk."""

    def test_rk_k0_same_z(self):
        """R^0 with same Z_eff: fast vs slow should agree."""
        Z = 2.0
        slow = _slater_rk(1, 0, Z, 1, 0, Z, 1, 0, Z, 1, 0, Z, 0)
        fast = _slater_rk_fast(1, 0, Z, 1, 0, Z, 1, 0, Z, 1, 0, Z, 0)
        assert abs(slow - fast) < 1e-4 * abs(slow), (
            f"R^0 mismatch: slow={slow:.8f}, fast={fast:.8f}"
        )

    def test_rk_k0_mixed_z(self):
        """R^0 with mixed Z_eff for generalized Sturmian case."""
        slow = _slater_rk(1, 0, 1.7, 1, 0, 1.7, 1, 0, 2.5, 1, 0, 2.5, 0)
        fast = _slater_rk_fast(1, 0, 1.7, 1, 0, 1.7, 1, 0, 2.5, 1, 0, 2.5, 0)
        assert abs(slow - fast) < 1e-4 * abs(slow), (
            f"Mixed-Z R^0 mismatch: slow={slow:.8f}, fast={fast:.8f}"
        )

    def test_rk_k1_p_orbitals(self):
        """R^1 with p orbitals."""
        slow = _slater_rk(2, 1, 3.0, 2, 1, 3.0, 2, 1, 3.0, 2, 1, 3.0, 1)
        fast = _slater_rk_fast(2, 1, 3.0, 2, 1, 3.0, 2, 1, 3.0, 2, 1, 3.0, 1)
        if abs(slow) > 1e-12:
            assert abs(slow - fast) / abs(slow) < 1e-3, (
                f"R^1 p-orbital mismatch: slow={slow:.8f}, fast={fast:.8f}"
            )


# ---------------------------------------------------------------------------
# Generalized Sturmian CI tests (Track BU-1b)
# ---------------------------------------------------------------------------

class TestGeneralizedSturmian:
    """Tests for the generalized Sturmian CI solver."""

    def test_beta_computation(self):
        """β values should be distinct for different (n1,n2) configs."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
        betas = solver._compute_betas(E_trial=-2.8)
        unique_betas = sorted(set(np.round(betas, 8)))
        # max_n=2: (1,1), (1,2), (2,2) → 3 distinct β values
        assert len(unique_betas) == 3, (
            f"Expected 3 unique betas, got {len(unique_betas)}: {unique_betas}"
        )

    def test_beta_ordering(self):
        """Configs with higher-n orbitals should get larger β."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
        betas = solver._compute_betas(E_trial=-2.8)
        # β = sqrt(-2E / (Z² × Σ 1/n²))
        # (1,1): Σ=2 → small β;  (2,2): Σ=0.5 → large β
        beta_by_cfg = {}
        for I, sd in enumerate(solver.sd_basis):
            ns = tuple(sorted([solver.states[sp >> 1][0] for sp in sd]))
            beta_by_cfg.setdefault(ns, betas[I])
        assert beta_by_cfg[(1, 1)] < beta_by_cfg[(1, 2)] < beta_by_cfg[(2, 2)]

    @pytest.mark.slow
    def test_he_max_n2_energy(self):
        """Generalized Sturmian He at max_n=2: reasonable energy."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
        result = solver.solve(E_trial=-2.81, max_iter=20, tol=1e-6,
                              damping=0.5, verbose=True)
        E = result['energy']
        E_exact = -2.903724
        # Must respect variational bound
        assert E > E_exact, f"Variational bound violated: {E:.6f}"
        # Should give a reasonable energy (within 10%)
        err_pct = 100 * abs(E - E_exact) / abs(E_exact)
        assert err_pct < 10.0, f"Error too large: {err_pct:.2f}%"

    @pytest.mark.slow
    def test_he_max_n2_convergence(self):
        """Self-consistency should converge within 20 iterations."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
        result = solver.solve(E_trial=-2.81, max_iter=20, tol=1e-6,
                              damping=0.5)
        assert result['converged'], (
            f"Not converged after {result['n_iterations']} iterations"
        )

    @pytest.mark.slow
    def test_he_max_n2_gaunt_preserved(self):
        """Angular ERI sparsity should match Coulomb Sturmian."""
        gen = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=2)
        sturm = SturmianCI(Z=2, n_electrons=2, max_n=2)
        sturm_r = sturm.solve(k=1.5)
        # Angular ERI count should be identical
        assert gen._eri_angular_count == sturm_r['eri_count'], (
            f"ERI count mismatch: gen={gen._eri_angular_count}, "
            f"sturm={sturm_r['eri_count']}"
        )

    @pytest.mark.slow
    def test_he_max_n3_energy(self):
        """Generalized Sturmian He at max_n=3: energy should be reasonable."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=3)
        result = solver.solve(E_trial=-2.844, max_iter=20, tol=1e-5,
                              damping=0.5, verbose=True)
        E = result['energy']
        E_exact = -2.903724
        assert E > E_exact, f"Variational bound violated: {E:.6f}"
        err_pct = 100 * abs(E - E_exact) / abs(E_exact)
        assert err_pct < 10.0, f"Error too large: {err_pct:.2f}%"

    @pytest.mark.slow
    def test_he_max_n3_beta_differentiation(self):
        """At max_n=3, configs should have measurably different β."""
        solver = GeneralizedSturmianCI(Z=2, n_electrons=2, max_n=3)
        result = solver.solve(E_trial=-2.844, max_iter=20, tol=1e-5,
                              damping=0.5)
        betas = result['betas']
        unique_betas = sorted(set(np.round(betas, 6)))
        # max_n=3: 6 distinct (n1,n2) types → 6 unique β
        assert len(unique_betas) == 6, (
            f"Expected 6 unique betas, got {len(unique_betas)}"
        )
        # β range should be substantial (not collapsed to single value)
        beta_range = max(unique_betas) / min(unique_betas)
        assert beta_range > 1.5, (
            f"β range too narrow: max/min = {beta_range:.2f}"
        )
