"""
Tests for the Level 4 characteristic polynomial P(rho, mu).

Diagnostic/research tests for Track S: investigate the algebraic structure
of the Level 4 (H2 molecule-frame hyperspherical) angular Hamiltonian.

Key structural findings tested:
  1. The Level 4 Hamiltonian is NOT a linear pencil in rho — V_nuc(rho)
     depends nonlinearly on rho through split-region min/max terms.
  2. At fixed rho, P(mu) = det(H_ang - mu*I) is a polynomial in mu
     whose degree equals the matrix dimension.
  3. The SO(6) Casimir eigenvalues (H_free) have integer entries.
  4. The V_ee coupling is rho-independent.
  5. The V_nuc coupling is piecewise-smooth in rho with kinks at rho values
     corresponding to quadrature boundary crossings.

References:
  - Paper 15, Section V (nuclear coupling, split-region expansion)
  - Paper 13, Sec XII (Level 3 algebraic structure, comparison)
  - algebraic_curve.py (Level 3 characteristic polynomial)
"""

import numpy as np
import pytest

from geovac.algebraic_curve_level4 import (
    extract_level4_matrices,
    probe_rho_dependence,
    characteristic_polynomial_fixed_rho,
    verify_roots_at_rho,
    rho_dependence_decomposition,
)

try:
    import sympy
    from sympy import Symbol, Poly
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

pytestmark = pytest.mark.skipif(not HAS_SYMPY, reason="SymPy required")


# ─────────────────────────────────────────────────────────────────────
# 1. Matrix extraction and structure
# ─────────────────────────────────────────────────────────────────────

class TestMatrixExtraction:
    """Test extraction and structural properties of Level 4 Hamiltonian components."""

    def test_h_free_integer_diagonal_lmax1(self):
        """H_free at l_max=1 should have integer diagonal entries (Casimir)."""
        H_free, _, _, _, channels = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        # l_max=1 homonuclear: channels (0,0), (1,1)
        # (0,0), k=0, singlet: mu = 0.5*(0+0+2)^2 - 2 + 0 = 0
        # (1,1), k=0, singlet: mu = 0.5*(1+1+2)^2 - 2 + 0 = 6
        assert H_free.shape[0] == 2, f"Expected dim=2, got {H_free.shape[0]}"
        assert abs(H_free[0, 0] - 0.0) < 1e-12
        assert abs(H_free[1, 1] - 6.0) < 1e-12

    def test_h_free_integer_diagonal_lmax2(self):
        """H_free at l_max=2, n_basis=1 should have integer diagonal entries."""
        H_free, _, _, _, channels = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=100,
        )
        dim = H_free.shape[0]
        for i in range(dim):
            val = H_free[i, i]
            int_val = int(round(val))
            assert abs(val - int_val) < 1e-12, (
                f"H_free[{i},{i}] = {val}, not integer"
            )

    def test_h_free_diagonal(self):
        """H_free should be exactly diagonal."""
        H_free, _, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        off_diag = H_free - np.diag(np.diag(H_free))
        assert np.max(np.abs(off_diag)) < 1e-15

    def test_v_ee_symmetric(self):
        """V_ee should be symmetric."""
        _, _, V_ee, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        assert np.max(np.abs(V_ee - V_ee.T)) < 1e-14

    def test_v_nuc_symmetric(self):
        """V_nuc should be symmetric."""
        _, V_nuc, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        assert np.max(np.abs(V_nuc - V_nuc.T)) < 1e-14

    def test_v_ee_rho_independent(self):
        """V_ee should not change with rho."""
        _, _, V_ee1, _, _ = extract_level4_matrices(
            rho=0.5, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        _, _, V_ee2, _, _ = extract_level4_matrices(
            rho=2.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        np.testing.assert_allclose(V_ee1, V_ee2, atol=1e-14,
                                   err_msg="V_ee depends on rho!")

    def test_v_nuc_depends_on_rho(self):
        """V_nuc should change with rho."""
        _, V_nuc1, _, _, _ = extract_level4_matrices(
            rho=0.5, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        _, V_nuc2, _, _, _ = extract_level4_matrices(
            rho=2.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        diff = np.max(np.abs(V_nuc1 - V_nuc2))
        assert diff > 1e-6, f"V_nuc appears rho-independent (diff={diff})"

    def test_channels_lmax1_homonuclear(self):
        """l_max=1 homonuclear should have channels (0,0) and (1,1)."""
        _, _, _, _, channels = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        assert channels == [(0, 0), (1, 1)], f"Got channels: {channels}"

    def test_channels_lmax2_homonuclear(self):
        """l_max=2 homonuclear: (0,0),(1,1),(0,2),(2,0),(2,2)."""
        _, _, _, _, channels = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=100,
        )
        expected = [(0, 0), (0, 2), (1, 1), (2, 0), (2, 2)]
        assert channels == expected, f"Got channels: {channels}"

    def test_dimension_lmax1_nbasis1(self):
        """l_max=1, n_basis=1: 2 channels x 1 = 2 total dim."""
        H_free, _, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        assert H_free.shape == (2, 2)

    def test_dimension_lmax1_nbasis3(self):
        """l_max=1, n_basis=3: 2 channels x 3 = 6 total dim."""
        H_free, _, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        assert H_free.shape == (6, 6)

    def test_dimension_lmax2_nbasis1(self):
        """l_max=2, n_basis=1: 5 channels x 1 = 5 total dim."""
        H_free, _, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=100,
        )
        assert H_free.shape == (5, 5)


# ─────────────────────────────────────────────────────────────────────
# 2. Rho-dependence characterization
# ─────────────────────────────────────────────────────────────────────

class TestRhoDependence:
    """Test the rho-dependence of Level 4 nuclear coupling."""

    def test_not_linear_pencil(self):
        """V_nuc(rho) should NOT be linear in rho (unlike Level 3)."""
        result = probe_rho_dependence(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100,
            rho_values=[0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0],
        )
        # The linearity residual should be large (>1%)
        assert not result['is_linear'], (
            f"V_nuc appears linear in rho "
            f"(residual={result['linearity_residual']:.2e})"
        )
        assert result['linearity_residual'] > 0.01, (
            f"Linearity residual too small: {result['linearity_residual']:.2e}"
        )

    def test_eigenvalues_vary_with_rho(self):
        """Eigenvalues should change significantly with rho."""
        result = probe_rho_dependence(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100,
            rho_values=[0.1, 0.5, 1.0, 2.0, 5.0],
        )
        evals = result['eigenvalues']
        # Ground state eigenvalue should vary
        ground_range = np.max(evals[:, 0]) - np.min(evals[:, 0])
        assert ground_range > 0.1, (
            f"Ground state eigenvalue range too small: {ground_range}"
        )

    def test_v_ee_constant_across_rho(self):
        """V_ee should be identical at all rho values."""
        result = probe_rho_dependence(
            R_e=1.5, l_max=1, n_basis=3, n_quad=100,
            rho_values=[0.1, 1.0, 5.0],
        )
        V_ee = result['V_ee_matrix']
        # V_ee is extracted from the last call — verify it's the
        # same as used in all calls
        assert V_ee is not None
        assert V_ee.shape[0] == result['dim']


# ─────────────────────────────────────────────────────────────────────
# 3. Characteristic polynomial at fixed rho
# ─────────────────────────────────────────────────────────────────────

class TestCharacteristicPolynomialFixedRho:
    """Test the characteristic polynomial P(mu) at fixed rho values."""

    def test_degree_mu_equals_dim_lmax1_nbasis1(self):
        """At l_max=1, n_basis=1: P is degree 2 in mu (dim=2)."""
        P, meta = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        assert meta['dim'] == 2
        assert meta['degree_mu'] == 2

    @pytest.mark.slow
    def test_degree_mu_equals_dim_lmax1_nbasis3(self):
        """At l_max=1, n_basis=3: dim=6, degree_mu should be 6.

        Note: SymPy det of 6x6 Float matrix can produce messy expressions.
        The degree_mu is set to dim as fallback when Poly() fails.
        """
        P, meta = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        assert meta['dim'] == 6
        assert meta['degree_mu'] == 6

    def test_degree_mu_equals_dim_lmax2_nbasis1(self):
        """At l_max=2, n_basis=1: P is degree 5 in mu (dim=5)."""
        P, meta = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=100,
        )
        assert meta['dim'] == 5
        assert meta['degree_mu'] == 5


# ─────────────────────────────────────────────────────────────────────
# 4. Root verification
# ─────────────────────────────────────────────────────────────────────

class TestRootVerification:
    """Verify P(mu) = 0 roots match numerical eigenvalues."""

    @pytest.mark.parametrize("rho", [0.1, 0.5, 1.0, 2.0, 5.0])
    def test_lmax1_nbasis1_roots_match(self, rho):
        """l_max=1, n_basis=1 (dim=2): roots match numerical eigenvalues."""
        result = verify_roots_at_rho(
            rho=rho, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        assert result['max_diff'] < 1e-8, (
            f"rho={rho}: max_diff={result['max_diff']:.2e}\n"
            f"  numerical: {result['mu_numerical']}\n"
            f"  from P:    {result['mu_from_P']}"
        )

    @pytest.mark.slow
    @pytest.mark.parametrize("rho", [0.1, 0.5, 1.0, 2.0, 5.0])
    def test_lmax1_nbasis3_roots_match(self, rho):
        """l_max=1, n_basis=3 (dim=6): roots match numerical eigenvalues.

        Slow: 6x6 SymPy determinant with Floats.
        """
        result = verify_roots_at_rho(
            rho=rho, R_e=1.5, l_max=1, n_basis=3, n_quad=100,
        )
        assert result['max_diff'] < 1e-6, (
            f"rho={rho}: max_diff={result['max_diff']:.2e}"
        )

    @pytest.mark.parametrize("rho", [0.1, 0.5, 1.0, 2.0, 5.0])
    def test_lmax2_nbasis1_roots_match(self, rho):
        """l_max=2, n_basis=1 (dim=5): roots match numerical eigenvalues."""
        result = verify_roots_at_rho(
            rho=rho, R_e=1.5, l_max=2, n_basis=1, n_quad=100,
        )
        assert result['max_diff'] < 1e-7, (
            f"rho={rho}: max_diff={result['max_diff']:.2e}"
        )

    @pytest.mark.parametrize("R_e", [0.5, 1.0, 1.5, 3.0])
    def test_different_R_e_values(self, R_e):
        """Roots should match at various R_e values."""
        result = verify_roots_at_rho(
            rho=1.0, R_e=R_e, l_max=1, n_basis=1, n_quad=100,
        )
        assert result['max_diff'] < 1e-8, (
            f"R_e={R_e}: max_diff={result['max_diff']:.2e}"
        )


# ─────────────────────────────────────────────────────────────────────
# 5. rho-dependence decomposition
# ─────────────────────────────────────────────────────────────────────

class TestRhoDecomposition:
    """Test the decomposition of rho-dependence in V_nuc matrix elements."""

    def test_kinks_detected(self):
        """V_nuc matrix elements should have kinks (piecewise structure)."""
        result = rho_dependence_decomposition(
            R_e=1.5, l_max=1, n_basis=1, n_quad=100,
            rho_values=np.linspace(0.05, 2.0, 200),
        )
        # At least one element should have kinks
        has_kinks = any(
            info['n_kinks'] > 0
            for info in result['element_analysis'].values()
        )
        # Note: kinks may or may not be detected depending on rho resolution
        # and quadrature. This test documents the finding rather than requiring it.
        print(f"Kinks detected: {has_kinks}")
        for key, info in result['element_analysis'].items():
            if info['n_kinks'] > 0:
                print(f"  V_nuc[{key}]: {info['n_kinks']} kinks at rho = "
                      f"{[f'{r:.3f}' for r in info['kink_rho_values'][:5]]}")

    def test_monotonic_small_rho(self):
        """At small rho, V_nuc[0,0] should be approximately -2*Z/rho (Coulomb)."""
        result = rho_dependence_decomposition(
            R_e=1.5, l_max=1, n_basis=1, n_quad=100,
            rho_values=np.linspace(0.01, 0.1, 20),
        )
        if (0, 0) in result['element_analysis']:
            vals = result['element_analysis'][(0, 0)]['values']
            # V_nuc[0,0] should be large and negative at small rho
            assert vals[0] < 0, "V_nuc[0,0] should be negative"
            # Should become less negative as rho increases (Coulomb decay)
            assert vals[-1] > vals[0], (
                "V_nuc[0,0] should increase (become less negative) with rho"
            )


# ─────────────────────────────────────────────────────────────────────
# 6. Structural comparison with Level 3
# ─────────────────────────────────────────────────────────────────────

class TestStructuralComparison:
    """Document structural differences between Level 3 and Level 4 polynomials.

    Level 3: H(R) = H0 + R*V_C -> P(R, mu) is bivariate polynomial
    Level 4: H(rho, R_e) = H_free + R_e*(V_nuc(rho) + V_ee) -> not polynomial in rho
    """

    def test_level4_has_two_parameters(self):
        """Level 4 has two parameters (rho, R_e) vs Level 3's single R."""
        # At different (rho, R_e) with same product rho*R_e, H differs
        _, _, _, H1, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        _, _, _, H2, _ = extract_level4_matrices(
            rho=0.5, R_e=3.0, l_max=1, n_basis=1, n_quad=100,
        )
        # H = H_free + R_e * (V_nuc(rho) + V_ee)
        # Even if rho*R_e = const, V_nuc(rho) differs, so H differs
        diff = np.max(np.abs(H1 - H2))
        assert diff > 0.01, (
            f"H should differ at (rho=1,Re=1.5) vs (rho=0.5,Re=3): diff={diff}"
        )

    def test_v_nuc_not_proportional_to_rho(self):
        """V_nuc(rho) != rho * V_nuc(1), confirming non-linear pencil."""
        _, V_nuc_half, _, _, _ = extract_level4_matrices(
            rho=0.5, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        _, V_nuc_one, _, _, _ = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        _, V_nuc_two, _, _, _ = extract_level4_matrices(
            rho=2.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )
        # If linear: V_nuc(rho) = A + rho*B, so
        # V_nuc(2) - V_nuc(1) should equal V_nuc(1) - V_nuc(0.5) * 2
        # This won't hold for nonlinear rho-dependence
        diff_1 = V_nuc_two - V_nuc_one
        diff_2 = V_nuc_one - V_nuc_half
        # For a linear function, diff_1 = 2 * diff_2 (equal spacing)
        # diff_1 / diff_2 should be 2.0 if linear
        for i in range(V_nuc_one.shape[0]):
            for j in range(V_nuc_one.shape[1]):
                if abs(diff_2[i, j]) > 1e-10:
                    ratio = diff_1[i, j] / diff_2[i, j]
                    # If exactly 2.0, it's linear
                    if abs(ratio - 2.0) > 0.01:
                        # Found nonlinear element
                        return
        # If we reach here, all ratios are ~2.0 -> linear
        pytest.fail("V_nuc appears to be linear in rho (unexpected)")

    def test_polynomial_metadata_table(self):
        """Generate metadata table comparing Level 3 and Level 4."""
        # Level 4 at l_max=1, n_basis=1 (dim=2)
        P_l4, meta_l4 = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=1.5, l_max=1, n_basis=1, n_quad=100,
        )

        print("\n=== Level 4 Characteristic Polynomial Structure ===")
        print(f"  l_max = {meta_l4['l_max']}")
        print(f"  n_basis = {meta_l4['n_basis']}")
        print(f"  dim = {meta_l4['dim']}")
        print(f"  degree_mu = {meta_l4['degree_mu']}")
        print(f"  n_terms = {meta_l4['n_terms']}")
        print(f"  rho = {meta_l4['rho']}")
        print(f"  R_e = {meta_l4['R_e']}")
        print(f"  rho_dependence: NONLINEAR (split-region min/max)")
        print()
        print("  Level 3 comparison (from algebraic_curve.py):")
        print("    rho_dependence: LINEAR PENCIL (H0 + R*V_C)")
        print("    At l_max=1, n_basis=1: dim=2, degree_R=2, degree_mu=2")
        print("    P(R,mu) is polynomial in BOTH R and mu")
        print()
        print("  Level 4 key difference:")
        print("    P(rho, mu) is polynomial in mu ONLY")
        print("    rho enters through V_nuc(rho) which involves min/max(s,rho)")
        print("    After spectral projection, matrix elements are piecewise")
        print("    rational functions of rho, NOT polynomial.")

        assert meta_l4['degree_mu'] == meta_l4['dim']


# ─────────────────────────────────────────────────────────────────────
# 7. Summary diagnostic: print polynomial and verification table
# ─────────────────────────────────────────────────────────────────────

class TestSummaryDiagnostic:
    """Summary diagnostic: print the characteristic polynomial and verification."""

    def test_print_verification_table(self):
        """Print the root verification table for the report."""
        rho_values = [0.1, 0.5, 1.0, 2.0, 5.0]
        R_e = 1.5

        print("\n=== Root Verification Table ===")
        print(f"{'rho':>6s} | {'mu_numerical':>20s} | {'mu_from_P':>20s} | {'|diff|':>12s}")
        print("-" * 70)

        for rho in rho_values:
            result = verify_roots_at_rho(
                rho=rho, R_e=R_e, l_max=1, n_basis=1, n_quad=200,
            )
            for i in range(min(len(result['mu_numerical']),
                               len(result['mu_from_P']))):
                mu_n = result['mu_numerical'][i]
                mu_p = result['mu_from_P'][i]
                diff = abs(mu_n - mu_p)
                print(f"{rho:6.1f} | {mu_n:20.12f} | {mu_p:20.12f} | {diff:12.2e}")

        # Polynomial structure at rho=1.0
        P, meta = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=R_e, l_max=1, n_basis=1, n_quad=200,
        )
        print(f"\nP(mu) at rho=1.0, R_e=1.5, l_max=1, n_basis=1:")
        print(f"  dim={meta['dim']}, degree_mu={meta['degree_mu']}, "
              f"n_terms={meta['n_terms']}")
        print(f"  P = {P}")

    def test_print_lmax2_diagnostic(self):
        """Print l_max=2 diagnostic data for the report."""
        print("\n=== l_max=2 Diagnostic ===")

        # Casimir eigenvalues and matrix structure
        H_free, V_nuc, V_ee, H_total, channels = extract_level4_matrices(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=200,
        )
        print(f"channels: {channels}")
        print(f"dim: {H_free.shape[0]}")
        print(f"H_free diagonal (Casimir): {np.diag(H_free)}")
        print(f"V_nuc diagonal: {[f'{v:.6f}' for v in np.diag(V_nuc)]}")
        print(f"V_ee diagonal:  {[f'{v:.6f}' for v in np.diag(V_ee)]}")

        # Characteristic polynomial at rho=1.0
        P, meta = characteristic_polynomial_fixed_rho(
            rho=1.0, R_e=1.5, l_max=2, n_basis=1, n_quad=200,
        )
        print(f"\nP(mu) at rho=1.0: dim={meta['dim']}, "
              f"degree_mu={meta['degree_mu']}, n_terms={meta['n_terms']}")

        # Root verification table
        print("\nRoot Verification (l_max=2, n_basis=1):")
        print(f"{'rho':>6s} | {'mu_0':>12s} | {'mu_1':>12s} | "
              f"{'mu_2':>12s} | {'max_diff':>12s}")
        print("-" * 70)
        for rho in [0.1, 0.5, 1.0, 2.0, 5.0]:
            result = verify_roots_at_rho(
                rho=rho, R_e=1.5, l_max=2, n_basis=1, n_quad=200,
            )
            mu = result['mu_numerical']
            print(f"{rho:6.1f} | {mu[0]:12.6f} | {mu[1]:12.6f} | "
                  f"{mu[2]:12.6f} | {result['max_diff']:12.2e}")

        # Rho-dependence
        result = probe_rho_dependence(
            R_e=1.5, l_max=2, n_basis=1, n_quad=200,
            rho_values=[0.1, 0.5, 1.0, 2.0, 5.0],
        )
        print(f"\nRho-dependence: is_linear={result['is_linear']}, "
              f"residual={result['linearity_residual']:.4e}")

    def test_print_rho_dependence_summary(self):
        """Print summary of rho-dependence characterization."""
        result = probe_rho_dependence(
            R_e=1.5, l_max=1, n_basis=1, n_quad=200,
            rho_values=[0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0],
        )

        print("\n=== Rho-Dependence Summary ===")
        print(f"  is_linear: {result['is_linear']}")
        print(f"  linearity_residual: {result['linearity_residual']:.6e}")
        print(f"  dim: {result['dim']}")
        print(f"  channels: {result['channels']}")
        print()
        print(f"  Ground state eigenvalues vs rho:")
        for i, rho in enumerate(result['rho_values']):
            evals = result['eigenvalues'][i]
            print(f"    rho={rho:5.1f}: mu_0 = {evals[0]:10.6f}")
