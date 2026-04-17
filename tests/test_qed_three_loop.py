"""Tests for three-loop QED spectral sums on S^3.

Tests the three-loop iterated sunset diagram: three electron lines connected
by two photon lines in a chain topology, with SO(4) vertex selection rules
at each vertex.
"""

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
