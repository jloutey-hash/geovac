"""Tests for three-loop QED spectral sums on S^3.

Tests the three-loop iterated sunset diagram: three electron lines connected
by two photon lines in a chain topology, with SO(4) vertex selection rules
at each vertex.

Includes tests for the O(N^3) factorized algorithm (Track Q-2).
"""

import time

import pytest
import mpmath

from geovac.qed_three_loop import (
    three_loop_sunset_s3,
    three_loop_unrestricted,
    three_loop_vertex_restricted,
    three_loop_cg_weighted,
    three_loop_convergence_table,
    decompose_three_loop,
    classify_three_loop_transcendentals,
    three_loop_factorized,
    three_loop_factorized_convergence,
    decompose_three_loop_mzv,
    three_loop_euler_maclaurin_tail,
)


class TestThreeLoopBasic:
    """Sanity checks for three-loop sums."""

    def test_sunset_positive(self):
        r = three_loop_sunset_s3(n_max=3)
        assert r["value_float"] > 0

    def test_sunset_finite(self):
        r = three_loop_sunset_s3(n_max=5)
        assert 0 < r["value_float"] < 1e10

    def test_sunset_has_quintuples(self):
        r = three_loop_sunset_s3(n_max=3)
        assert r["n_quintuples"] > 0

    def test_unrestricted_is_D4_cubed(self):
        u = three_loop_unrestricted(electron_exponent=4)
        D4 = mpmath.mpf(u["D_a_float"])
        expected = D4**3
        assert abs(u["D_a_cubed_float"] - float(expected)) < 1e-10

    def test_unrestricted_is_pi_even(self):
        u = three_loop_unrestricted(electron_exponent=4)
        assert u["transcendental_class"] == "pi^{even}"

    def test_unrestricted_odd_exponent(self):
        u = three_loop_unrestricted(electron_exponent=5)
        assert u["transcendental_class"] == "odd-zeta"


class TestThreeLoopConvergence:
    """Convergence tests."""

    def test_monotone_increasing(self):
        vals = []
        for n in [3, 5, 7]:
            r = three_loop_sunset_s3(n_max=n)
            vals.append(r["value_float"])
        for i in range(1, len(vals)):
            assert vals[i] > vals[i - 1]

    def test_convergence_table(self):
        table = three_loop_convergence_table(n_max_values=[3, 5, 7])
        assert len(table) == 3
        assert table[0]["delta_from_previous"] is None
        assert table[1]["delta_from_previous"] > 0
        assert table[2]["delta_from_previous"] > 0

    def test_cg_weights_differ_from_uniform(self):
        cg = three_loop_sunset_s3(n_max=5, use_cg_weights=True)
        uniform = three_loop_sunset_s3(n_max=5, use_cg_weights=False)
        assert abs(cg["value_float"] - uniform["value_float"]) > 0.01


class TestThreeLoopStructure:
    """Structural tests for the three-loop topology."""

    def test_vertex_restricted_below_unrestricted_at_low_n(self):
        """At low n_max, vertex restriction removes terms."""
        vr = three_loop_vertex_restricted(n_max=5)
        u = three_loop_unrestricted()
        # Vertex-restricted should be less than unrestricted at low n
        assert vr["value_float"] < u["D_a_cubed_float"]

    def test_cg_weighted_vs_unrestricted(self):
        """CG-weighted sum differs from unrestricted D(4)^3."""
        cg = three_loop_cg_weighted(n_max=5)
        assert "unrestricted_float" in cg
        assert "difference_float" in cg
        assert abs(cg["difference_float"]) > 0.01

    def test_classification_has_expected_keys(self):
        c = classify_three_loop_transcendentals()
        assert "unrestricted_three_loop" in c
        assert "vertex_restricted_three_loop" in c
        assert "cg_weighted_three_loop" in c
        assert "depth_k_prediction" in c

    def test_decompose_returns_dict(self):
        r = three_loop_sunset_s3(n_max=3)
        d = decompose_three_loop(r["value"])
        assert isinstance(d, dict)
        assert "identified" in d


class TestThreeLoopCustomExponents:
    """Test non-default exponent choices."""

    def test_higher_electron_exponent(self):
        r = three_loop_sunset_s3(n_max=5, electron_exponent=6)
        assert r["value_float"] > 0
        assert r["electron_exponent"] == 6

    def test_higher_photon_exponent(self):
        r = three_loop_sunset_s3(n_max=5, photon_exponent=2)
        assert r["value_float"] > 0
        assert r["photon_exponent"] == 2


@pytest.mark.slow
class TestThreeLoopSlow:
    """Slow convergence and PSLQ tests."""

    def test_convergence_n15(self):
        r = three_loop_sunset_s3(n_max=15)
        assert r["value_float"] > 0
        assert r["n_quintuples"] > 100000

    def test_vertex_restricted_pslq_at_n15(self):
        """At n_max=15, try PSLQ on vertex-restricted sum.

        This is expected to FAIL (insufficient convergence) — documenting
        that the sum is not yet converged enough for identification.
        """
        vr = three_loop_vertex_restricted(n_max=15)
        # We don't assert identification succeeds — just that it runs
        assert "decomposition" in vr

    def test_cg_weighted_convergence_n20(self):
        """CG-weighted at n_max=20 for convergence assessment."""
        r = three_loop_sunset_s3(n_max=20)
        # Should be approaching D(4)^3 ~ 5.38 from below or above
        assert r["value_float"] > 3.0


# ============================================================================
# Track Q-2: Factorized O(N^3) three-loop sum
# ============================================================================

class TestFactorizedMatchesDirect:
    """Verify the factorized O(N^3) algorithm matches the direct O(N^5) sum."""

    def test_factorized_matches_direct_cg_n5(self):
        """Factorized CG-weighted at n_max=5 matches direct to high precision."""
        direct = three_loop_sunset_s3(
            n_max=5, electron_exponent=4, photon_exponent=1,
            use_cg_weights=True,
        )
        factorized = three_loop_factorized(
            n_max=5, a=4, p=1, use_cg_weights=True,
        )
        # Must match to at least 25 significant digits
        diff = abs(direct["value"] - factorized["value"])
        assert diff < mpmath.mpf(10) ** (-25), (
            f"Mismatch: direct={float(direct['value'])}, "
            f"factorized={float(factorized['value'])}, diff={float(diff)}"
        )

    def test_factorized_matches_direct_cg_n10(self):
        """Factorized CG-weighted at n_max=10 matches direct to 1e-25."""
        direct = three_loop_sunset_s3(
            n_max=10, electron_exponent=4, photon_exponent=1,
            use_cg_weights=True,
        )
        factorized = three_loop_factorized(
            n_max=10, a=4, p=1, use_cg_weights=True,
        )
        diff = abs(direct["value"] - factorized["value"])
        assert diff < mpmath.mpf(10) ** (-25), (
            f"Mismatch at n_max=10: diff={float(diff)}"
        )

    def test_factorized_matches_direct_vertex_restricted(self):
        """Factorized with W=1 (vertex-restricted) matches direct."""
        direct = three_loop_sunset_s3(
            n_max=7, electron_exponent=4, photon_exponent=1,
            use_cg_weights=False,
        )
        factorized = three_loop_factorized(
            n_max=7, a=4, p=1, use_cg_weights=False,
        )
        diff = abs(direct["value"] - factorized["value"])
        assert diff < mpmath.mpf(10) ** (-25), (
            f"Vertex-restricted mismatch at n_max=7: diff={float(diff)}"
        )

    def test_factorized_matches_direct_higher_exponents(self):
        """Factorized matches direct with non-default exponents."""
        direct = three_loop_sunset_s3(
            n_max=5, electron_exponent=6, photon_exponent=2,
            use_cg_weights=True,
        )
        factorized = three_loop_factorized(
            n_max=5, a=6, p=2, use_cg_weights=True,
        )
        diff = abs(direct["value"] - factorized["value"])
        assert diff < mpmath.mpf(10) ** (-25), (
            f"Higher exponents mismatch: diff={float(diff)}"
        )


class TestFactorizedConvergence:
    """Convergence and performance tests for the factorized sum."""

    def test_factorized_convergence_monotonic(self):
        """Value changes monotonically with n_max (for the default sum)."""
        vals = []
        for n in [5, 10, 15]:
            r = three_loop_factorized(n_max=n, a=4, p=1, use_cg_weights=True)
            vals.append(r["value"])
        # The sum should increase monotonically (more terms, all positive contributions)
        for i in range(1, len(vals)):
            assert vals[i] > vals[i - 1], (
                f"Non-monotone: S({[5,10,15][i-1]})={float(vals[i-1])}, "
                f"S({[5,10,15][i]})={float(vals[i])}"
            )

    def test_factorized_speed_n50(self):
        """n_max=50 completes in < 60 seconds."""
        t0 = time.perf_counter()
        r = three_loop_factorized(n_max=50, a=4, p=1, use_cg_weights=True)
        elapsed = time.perf_counter() - t0
        assert r["value_float"] > 0
        assert elapsed < 60, f"n_max=50 took {elapsed:.1f}s (limit: 60s)"

    def test_factorized_convergence_table_structure(self):
        """Convergence table has expected structure."""
        table = three_loop_factorized_convergence(
            n_max_values=[5, 10, 15], a=4, p=1,
        )
        assert len(table) == 3
        assert table[0]["delta_from_previous"] is None
        assert table[1]["delta_from_previous"] > 0
        assert table[2]["delta_from_previous"] > 0
        # Delta should decrease (convergence)
        assert table[2]["delta_from_previous"] < table[1]["delta_from_previous"]
        # All entries have timing
        for row in table:
            assert "time_seconds" in row
            assert row["time_seconds"] > 0


class TestEulerMaclaurinTail:
    """Tests for the Euler-Maclaurin tail correction."""

    def test_tail_correction_positive_direction(self):
        """For this monotonically increasing sum, tail should be positive."""
        tail = three_loop_euler_maclaurin_tail(n_max=15, n_em_terms=5, a=4, p=1)
        # The sum is increasing, so the tail estimate should be positive
        assert float(tail) > 0, f"Tail correction is negative: {float(tail)}"

    def test_tail_correction_finite(self):
        """Tail correction returns a finite value."""
        tail = three_loop_euler_maclaurin_tail(n_max=15, n_em_terms=5, a=4, p=1)
        assert mpmath.isfinite(tail)

    def test_tail_correction_converging_sum(self):
        """For a well-converging sum (higher exponent), tail improves estimate."""
        # With a=6, p=2 the sum converges faster
        r20 = three_loop_factorized(n_max=20, a=6, p=2, use_cg_weights=True)
        r15 = three_loop_factorized(n_max=15, a=6, p=2, use_cg_weights=True)
        gap_raw = abs(float(r20["value"] - r15["value"]))

        tail = three_loop_euler_maclaurin_tail(n_max=15, n_em_terms=5, a=6, p=2)
        corrected = r15["value"] + tail
        gap_corrected = abs(float(r20["value"] - corrected))
        # For this faster-converging sum, the tail should help
        assert gap_corrected < gap_raw, (
            f"Tail correction did not help for fast-converging sum: "
            f"raw gap={gap_raw}, corrected gap={gap_corrected}"
        )


class TestDecomposeMZV:
    """Tests for the extended MZV PSLQ decomposition."""

    def test_decompose_known_pi_even(self):
        """A known pi^{even} value should be identified."""
        # D(4) = pi^2 - pi^4/12, a simple weight-4 pi^{even} value
        from geovac.qed_three_loop import _dirac_D
        D4 = _dirac_D(4)
        result = decompose_three_loop_mzv(D4, n_digits=50)
        assert result["identified"], "Failed to identify D(4)"
        assert result["basis_tier"] == 0

    def test_decompose_returns_dict(self):
        """decompose_three_loop_mzv returns a well-formed dict."""
        val = mpmath.mpf("3.14159265358979")
        result = decompose_three_loop_mzv(val, n_digits=14)
        assert isinstance(result, dict)
        assert "identified" in result

    def test_decompose_rational(self):
        """A rational number should be identified as such."""
        val = mpmath.mpf("42") / mpmath.mpf("17")
        result = decompose_three_loop_mzv(val, n_digits=50)
        assert result["identified"]
        assert "1" in result.get("components", {})


@pytest.mark.slow
class TestFactorizedSlow:
    """Slow factorized tests (n_max=100+)."""

    def test_factorized_n100(self):
        """n_max=100 completes and gives a positive value."""
        r = three_loop_factorized(n_max=100, a=4, p=1, use_cg_weights=True)
        assert r["value_float"] > 0
        assert r["time_seconds"] < 600  # 10 minute limit
