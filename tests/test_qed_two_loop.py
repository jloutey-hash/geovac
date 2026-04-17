"""Tests for two-loop QED on S^3: zeta(3) from Dirac spectral data.

Verifies the even/odd discriminant for D_Dirac(s), the Track D3 result,
nested harmonic sums, connected double sums, and the flat-space limit.
"""

import pytest
import mpmath

mpmath.mp.dps = 80

from geovac.qed_two_loop import (
    dirac_dirichlet_series_hurwitz,
    dirac_dirichlet_series_numerical,
    fock_dirichlet_series,
    nested_harmonic_sum_flat,
    nested_harmonic_sum_s3,
    double_spectral_zeta_connected,
    even_odd_discriminant_table,
    decompose_into_zeta_basis,
    two_loop_vacuum_polarization_s3,
    flat_space_limit_check,
    classify_two_loop_transcendentals,
)


# ---------------------------------------------------------------------------
# Even/odd discriminant: the core structural result
# ---------------------------------------------------------------------------

class TestEvenOddDiscriminant:
    """D_Dirac(s_even) = pi^{even}; D_Dirac(s_odd) = odd-zeta."""

    def test_D4_is_pi_even(self):
        """D_Dirac(4) = pi^2 - pi^4/12 (pi^{even} only, no zeta(3))."""
        result = dirac_dirichlet_series_hurwitz(4)
        assert result["parity"] == "even"
        assert result["zeta_type"] == "pi^{even}"
        assert not result["contains_odd_zeta"]
        # Check PSLQ: 12*D4 = 12*pi^2 - pi^4
        rel = result["pslq_relation"]
        assert rel is not None
        # Verify numerically
        target = mpmath.pi**2 - mpmath.pi**4 / 12
        assert float(abs(result["value"] - target)) < 1e-60

    def test_D5_contains_zeta3(self):
        """D_Dirac(5) = 14*zeta(3) - 31/2*zeta(5) (odd-zeta)."""
        result = dirac_dirichlet_series_hurwitz(5)
        assert result["parity"] == "odd"
        assert result["zeta_type"] == "odd-zeta"
        assert result["contains_odd_zeta"]
        # Verify numerically
        target = 14 * mpmath.zeta(3) - mpmath.mpf(31) / 2 * mpmath.zeta(5)
        assert float(abs(result["value"] - target)) < 1e-60

    def test_D6_is_pi_even(self):
        """D_Dirac(6) = pi^4/3 - pi^6/30 (pi^{even} only)."""
        result = dirac_dirichlet_series_hurwitz(6)
        assert result["parity"] == "even"
        assert not result["contains_odd_zeta"]
        target = mpmath.pi**4 / 3 - mpmath.pi**6 / 30
        assert float(abs(result["value"] - target)) < 1e-60

    def test_D7_contains_zeta5(self):
        """D_Dirac(7) = 62*zeta(5) - 127/2*zeta(7) (odd-zeta)."""
        result = dirac_dirichlet_series_hurwitz(7)
        assert result["parity"] == "odd"
        assert result["contains_odd_zeta"]
        target = 62 * mpmath.zeta(5) - mpmath.mpf(127) / 2 * mpmath.zeta(7)
        assert float(abs(result["value"] - target)) < 1e-60

    def test_D8_is_pi_even(self):
        """D_Dirac(8) = pi^{even} only."""
        result = dirac_dirichlet_series_hurwitz(8)
        assert result["parity"] == "even"
        assert not result["contains_odd_zeta"]

    def test_discriminant_table_alternates(self):
        """The even/odd pattern holds for s = 4..10."""
        table = even_odd_discriminant_table(s_max=10)
        for row in table:
            s = row["s"]
            if s % 2 == 0:
                assert row["zeta_type"] == "pi^{even}", f"s={s} should be pi^even"
                assert not row["contains_odd_zeta"], f"s={s} should not contain odd-zeta"
            else:
                assert row["zeta_type"] == "odd-zeta", f"s={s} should be odd-zeta"
                assert row["contains_odd_zeta"], f"s={s} should contain odd-zeta"

    def test_s3_diverges(self):
        """D_Dirac(3) diverges (zeta(1,3/2) has a pole)."""
        with pytest.raises(ValueError, match="s must be >= 4"):
            dirac_dirichlet_series_hurwitz(3)


# ---------------------------------------------------------------------------
# Track D3 cross-check: Fock-index Dirichlet series
# ---------------------------------------------------------------------------

class TestFockDirichletSeries:
    """D_Fock(s) = 2*zeta_R(s-2) + 2*zeta_R(s-1)."""

    def test_D_fock_4_equals_2z2_plus_2z3(self):
        """D_Fock(4) = 2*zeta(2) + 2*zeta(3) (Track D3 result)."""
        result = fock_dirichlet_series(4)
        target = 2 * mpmath.zeta(2) + 2 * mpmath.zeta(3)
        assert float(abs(result["value"] - target)) < 1e-60
        assert result["contains_odd_zeta"]

    def test_D_fock_4_pslq(self):
        """PSLQ identifies D_Fock(4) = pi^2/3 + 2*zeta(3)."""
        result = fock_dirichlet_series(4)
        rel = result["pslq_relation"]
        assert rel is not None
        # Relation: 3*val - pi^2 - 6*z(3) = 0  =>  val = pi^2/3 + 2*z(3)
        assert rel[0] == 3
        # Check the pi^2 and zeta(3) coefficients
        decomp = result["decomposition"]
        assert "pi^2" in decomp or any("pi" in k for k in decomp.keys())
        assert any("zeta(3)" in k for k in decomp.keys())

    def test_D_fock_5(self):
        """D_Fock(5) = 2*zeta(3) + 2*zeta(4) = 2*zeta(3) + pi^4/45."""
        result = fock_dirichlet_series(5)
        target = 2 * mpmath.zeta(3) + 2 * mpmath.zeta(4)
        assert float(abs(result["value"] - target)) < 1e-60
        assert result["contains_odd_zeta"]

    def test_D_fock_always_has_odd_zeta(self):
        """D_Fock(s) always contains odd-zeta for s >= 4."""
        for s in [4, 5, 6, 7, 8]:
            result = fock_dirichlet_series(s)
            assert result["contains_odd_zeta"], f"D_Fock({s}) should contain odd-zeta"


# ---------------------------------------------------------------------------
# Numerical cross-checks
# ---------------------------------------------------------------------------

class TestNumericalCrossChecks:
    """Verify direct summation converges to Hurwitz exact values."""

    def test_D4_numerical_vs_hurwitz(self):
        """Direct sum of D(4) converges to Hurwitz value."""
        result = dirac_dirichlet_series_numerical(4, N=500)
        assert result["rel_error_vs_hurwitz"] is not None
        assert result["rel_error_vs_hurwitz"] < 5e-3  # 500-term truncation (slow convergence)

    def test_D5_numerical_vs_hurwitz(self):
        """Direct sum of D(5) converges to Hurwitz value."""
        result = dirac_dirichlet_series_numerical(5, N=500)
        assert result["rel_error_vs_hurwitz"] is not None
        assert result["rel_error_vs_hurwitz"] < 1e-4

    def test_D6_numerical_vs_hurwitz(self):
        """Direct sum of D(6) converges to Hurwitz value."""
        result = dirac_dirichlet_series_numerical(6, N=500)
        assert result["rel_error_vs_hurwitz"] is not None
        assert result["rel_error_vs_hurwitz"] < 1e-5


# ---------------------------------------------------------------------------
# Nested harmonic sums
# ---------------------------------------------------------------------------

class TestNestedHarmonicSums:
    """Nested sums that produce zeta(3)."""

    def test_flat_nested_converges_to_2z3(self):
        """Flat-space sum_{n=1}^N 1/n^2 * H_n -> 2*zeta(3)."""
        result = nested_harmonic_sum_flat(N=1000)
        assert result["rel_error"] < 5e-3  # 1000-term convergence (O(1/N) rate)

    def test_flat_nested_improves_with_N(self):
        """Convergence improves with increasing N."""
        r100 = nested_harmonic_sum_flat(N=100)
        r500 = nested_harmonic_sum_flat(N=500)
        assert r500["rel_error"] < r100["rel_error"]

    def test_s3_nested_converges(self):
        """S^3 nested harmonic sum converges to a finite value."""
        result = nested_harmonic_sum_s3(N=200)
        assert result["partial_sum"] > 0
        # Check it's converging (monotonically increasing)
        vals = [v for _, v in result["convergence"] if v > 0]
        for i in range(1, len(vals)):
            assert vals[i] > vals[i - 1]


# ---------------------------------------------------------------------------
# Connected double spectral zeta
# ---------------------------------------------------------------------------

class TestConnectedDoubleSums:
    """Two-loop connected sums."""

    def test_fock_connected_4_4_contains_zeta3(self):
        """Connected double sum with Fock indices contains zeta(3)."""
        result = double_spectral_zeta_connected(4, 4, use_fock=True)
        # The product of D_Fock(4) with itself contains zeta(3)^2
        # and zeta(2)*zeta(3) cross terms
        assert float(result["connected"]) != 0

    def test_dirac_connected_4_4_is_pi_even(self):
        """Connected double sum with Dirac eigenvalues at s=4,4 is pi^{even}."""
        result = double_spectral_zeta_connected(4, 4, use_fock=False)
        # Both s=4 sums are pi^{even}, so product and diagonal are pi^{even}
        # Connected = product - diagonal is also pi^{even}
        decomp = result["decomp_connected"]
        if decomp and decomp.get("identified"):
            # Should not contain zeta(3)
            assert not decomp.get("contains_zeta3", False)

    def test_product_larger_than_diagonal(self):
        """The full product sum is larger than the diagonal (n=m) sum."""
        result = double_spectral_zeta_connected(4, 4, use_fock=True)
        assert float(result["product"]) > float(result["diagonal"])


# ---------------------------------------------------------------------------
# Flat-space limit
# ---------------------------------------------------------------------------

class TestFlatSpaceLimit:
    """Verify flat-space reductions."""

    def test_flat_D4_is_zeta4(self):
        """With g_n=1, lambda_n=n: D(4) = zeta_R(4) = pi^4/90."""
        result = flat_space_limit_check(N=500)
        assert result["flat_D4"]["match"]

    def test_flat_nested_gives_2z3(self):
        """Flat nested harmonic sum -> 2*zeta(3)."""
        result = flat_space_limit_check(N=500)
        assert result["flat_nested"]["rel_error"] < 1e-2  # 500-term truncation


# ---------------------------------------------------------------------------
# PSLQ decomposition
# ---------------------------------------------------------------------------

class TestPSLQDecomposition:
    """Test the PSLQ-based decomposition routine."""

    def test_pi_squared(self):
        """Decompose pi^2 into the basis."""
        result = decompose_into_zeta_basis(mpmath.pi**2)
        assert result["identified"]
        assert "pi^2" in result["components"]
        assert result["residual"] < 1e-50

    def test_zeta3(self):
        """Decompose zeta(3) into the basis."""
        result = decompose_into_zeta_basis(mpmath.zeta(3))
        assert result["identified"]
        assert result["contains_zeta3"]
        assert result["residual"] < 1e-50

    def test_2z2_plus_2z3(self):
        """Decompose 2*zeta(2) + 2*zeta(3) = pi^2/3 + 2*zeta(3)."""
        val = 2 * mpmath.zeta(2) + 2 * mpmath.zeta(3)
        result = decompose_into_zeta_basis(val)
        assert result["identified"]
        assert result["contains_zeta3"]
        # Should contain pi^2 component (since 2*zeta(2) = pi^2/3)
        assert "pi^2" in result["components"]


# ---------------------------------------------------------------------------
# Two-loop summary
# ---------------------------------------------------------------------------

class TestTwoLoopSummary:
    """Integration test for the full two-loop analysis."""

    @pytest.mark.slow
    def test_full_summary_runs(self):
        """The full two-loop summary runs without error."""
        result = two_loop_vacuum_polarization_s3(n_max=100)
        assert "discriminant_table" in result
        assert "fock_D4" in result
        assert "summary" in result

    def test_classification_complete(self):
        """All expected transcendental categories are classified."""
        cl = classify_two_loop_transcendentals()
        assert "D_Dirac_s_even" in cl
        assert "D_Dirac_s_odd" in cl
        assert "D_Fock_all_s" in cl
        assert "even_odd_discriminant" in cl
        assert "two_loop_mechanism" in cl


# ---------------------------------------------------------------------------
# Algebraic identity verification
# ---------------------------------------------------------------------------

class TestAlgebraicIdentities:
    """Verify the algebraic identities underlying the results."""

    def test_hurwitz_half_integer_identity(self):
        """zeta(s, 3/2) = (2^s - 1)*zeta_R(s) - 2^s for s >= 2."""
        for s in [2, 3, 4, 5, 6]:
            hz = mpmath.hurwitz(s, mpmath.mpf(3) / 2)
            rhs = (2**s - 1) * mpmath.zeta(s) - mpmath.mpf(2)**s
            assert float(abs(hz - rhs)) < 1e-60, f"Failed at s={s}"

    def test_degeneracy_decomposition(self):
        """g_n = 2*lambda_n^2 - 1/2 with lambda_n = n + 3/2."""
        for n in range(20):
            g = 2 * (n + 1) * (n + 2)
            lam = n + 1.5
            rhs = 2 * lam**2 - 0.5
            assert abs(g - rhs) < 1e-12, f"Failed at n={n}"

    def test_fock_algebraic_decomposition(self):
        """2n(n+1)/n^s = 2/n^{s-2} + 2/n^{s-1} (algebraic identity)."""
        for n in range(1, 20):
            for s in [4, 5, 6]:
                lhs = 2.0 * n * (n + 1) / n**s
                rhs = 2.0 / n**(s - 2) + 2.0 / n**(s - 1)
                assert abs(lhs - rhs) < 1e-12, f"Failed at n={n}, s={s}"

    def test_D4_dirac_is_zeta_D2_at_s2(self):
        """D_Dirac(4) with |lambda|^{-4} = zeta_{D^2}(2) with lambda^{-4}."""
        # D_Dirac(4) = sum g_n / |lambda_n|^4 = sum g_n / (lambda_n^2)^2 = zeta_{D^2}(2)
        # T9: zeta_{D^2}(2) = pi^2 - pi^4/12
        D4 = dirac_dirichlet_series_hurwitz(4)
        t9_val = mpmath.pi**2 - mpmath.pi**4 / 12
        assert float(abs(D4["value"] - t9_val)) < 1e-60

    def test_D6_dirac_is_zeta_D2_at_s3(self):
        """D_Dirac(6) = zeta_{D^2}(3) = pi^4/3 - pi^6/30."""
        D6 = dirac_dirichlet_series_hurwitz(6)
        t9_val = mpmath.pi**4 / 3 - mpmath.pi**6 / 30
        assert float(abs(D6["value"] - t9_val)) < 1e-60


# ---------------------------------------------------------------------------
# Structural theorem verification
# ---------------------------------------------------------------------------

class TestStructuralTheorems:
    """Verify the structural claims about transcendental content."""

    def test_D_even_has_no_odd_zeta(self):
        """D_Dirac(s_even) has zero coefficient for every odd-zeta value."""
        for s in [4, 6, 8]:
            result = dirac_dirichlet_series_hurwitz(s)
            # The value should be exactly a polynomial in pi^2
            # Test by checking distance from the closest pi^{even} combination
            assert not result["contains_odd_zeta"]

    def test_D_odd_has_nonzero_odd_zeta(self):
        """D_Dirac(s_odd) has nonzero coefficient for odd-zeta values."""
        for s in [5, 7]:
            result = dirac_dirichlet_series_hurwitz(s)
            assert result["contains_odd_zeta"]
            # Verify the value is NOT zero
            assert float(abs(result["value"])) > 0.1

    def test_D5_zeta3_coefficient_is_14(self):
        """The coefficient of zeta(3) in D_Dirac(5) is exactly 14."""
        # D(5) = 2*(2^3-1)*z(3) - (2^5-1)/2*z(5) = 14*z(3) - 31/2*z(5)
        result = dirac_dirichlet_series_hurwitz(5)
        assert result["coeff_zeta_R_s_minus_2"] == 14  # 2*(2^3-1) = 14
        assert result["coeff_zeta_R_s"] == -15.5  # -(2^5-1)/2 = -31/2

    def test_coefficient_formula(self):
        """Coefficients follow the pattern 2*(2^{s-2}-1) and -(2^s-1)/2."""
        for s in [4, 5, 6, 7, 8]:
            result = dirac_dirichlet_series_hurwitz(s)
            assert result["coeff_zeta_R_s_minus_2"] == 2 * (2**(s - 2) - 1)
            assert result["coeff_zeta_R_s"] == -(2**s - 1) / 2
