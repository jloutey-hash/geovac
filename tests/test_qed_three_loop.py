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
    """Tests for the Euler-Maclaurin tail correction.

    The p=1 physical regime is a DOCUMENTED dead-end (log-modulated
    (ln N)^2/N^2 increments; CLAUDE.md sec. 3 / Paper 28's retraction
    of the N^-1.31 'power law').  Rewritten 2026-07-04 (cert-run gap
    closure): the old test asserted only tail > 0 at the broken default
    regime -- vacuously true for a wrong value.  The tests now assert
    the estimator WARNS there and PIN the documented failure (the
    vacuous-furry remediation pattern, v4.64.0)."""

    def test_tail_warns_in_log_modulated_regime(self):
        """At the physical a=4, p=1 the estimator emits RuntimeWarning."""
        with pytest.warns(RuntimeWarning, match="log-modulated"):
            tail = three_loop_euler_maclaurin_tail(
                n_max=15, n_em_terms=5, a=4, p=1)
        assert mpmath.isfinite(tail)

    def test_tail_estimate_wrong_at_physical_exponents(self):
        """Documented failure, pinned: at n_max=20, a=4, p=1 the estimate
        (~0.47) undershoots the true remaining tail S_inf - S(20) =
        40.6554... - 12.8366... ~ 27.8 by ~60x (and overshoots ~2.3x at
        n_max=50).  Corrected W, 2026-11-07 (old-count tail 31.5726 -
        7.8719).  The estimator must NOT be trusted at the physical
        exponents; the S^(3) headline uses the decomposition route."""
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore", RuntimeWarning)
            tail = three_loop_euler_maclaurin_tail(
                n_max=20, n_em_terms=5, a=4, p=1)
        # corrected W (2026-11-07): Paper 28 eq:s3_value minus the N=20 row
        true_tail = 40.6553875636 - 12.8365763232
        ratio = float(tail) / true_tail
        assert ratio < 1.0 / 3.0 or ratio > 3.0, (
            f"estimator unexpectedly accurate at the physical exponents "
            f"(ratio {ratio:.3f}) -- if genuinely fixed, update Paper 28 "
            f"and CLAUDE.md sec. 3 before trusting it")

    def test_true_tail_literals_pinned(self):
        """The two literals in true_tail above, backed (Paper 28 eq:s3_value
        and tab:three_loop_convergence, N = 20 row).

        Rejects: the old-count pair 31.5726 / 7.8719 (and any S(20) that the
        module computes with the old channel count, which gives 7.8719).
        The ratio test above passed with the old literals too, so on its
        own it does not see the channel count.

        Route: S(20) = sum_{m<=20} phi(m) L(m)^2, L(m) = 2 sum_{n<=20}
        min(n,m) phi(n), in exact Fraction arithmetic (no geovac), against
        the literal and against the module's factorized sum (which calls
        so4_channel_count).  S_inf against the 25-term closed form
        eq:s3_closed_form is backed in tests/test_s3_decomposition.py; here
        it is checked against the rigorous direct-sum bracket at N = 2000
        (lower = partial sum, upper = partial + 32[(ln N+1)^2+2(ln N+1)+2]/N).
        """
        from fractions import Fraction
        import math

        def phi(n):
            o = 2 * n + 3
            return Fraction(8 * (o * o - 1), o ** 4)
        N = 20
        ph = [phi(n) for n in range(N + 1)]
        exact = sum(ph[m] * (2 * sum(min(n, m) * ph[n] for n in range(N + 1))) ** 2
                    for m in range(N + 1))
        assert abs(float(exact) - 12.8365763232) < 1e-10
        assert abs(float(three_loop_factorized(N)["value_float"])
                   - float(exact)) < 1e-10
        # S_inf literal inside a rigorous bracket (float64 is ample here)
        Nb = 2000
        phs = [8 * ((2 * m + 3) ** 2 - 1) / (2 * m + 3) ** 4 for m in range(Nb + 1)]
        P_inf = math.pi ** 2 - math.pi ** 4 / 12
        head, cum_phi, cum_nphi = 0.0, 0.0, 0.0
        for m in range(Nb + 1):
            A = cum_nphi + m * (P_inf - cum_phi)
            head += phs[m] * (2 * A) ** 2
            cum_phi += phs[m]
            cum_nphi += m * phs[m]
        u = math.log(Nb) + 1
        hi = head + 32 * (u * u + 2 * u + 2) / Nb
        assert head < 40.6553875636 < hi
        assert not (head < 31.5726 < hi)

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
        # Local precision: under xdist an earlier test may have lowered the
        # global mp.dps, and D(4) at 15 digits cannot be identified at 50.
        with mpmath.workdps(60):
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
