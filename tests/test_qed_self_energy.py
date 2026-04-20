"""Tests for one-loop electron self-energy on S^3.

Validates the spectral mode sum implementation in geovac.qed_self_energy.

Key physical result discovered during test design:
  Sigma(n_ext=0) = 0 IDENTICALLY due to the vertex selection rule.
  At n_ext=0, the triangle inequality forces q = n_int, and the parity
  rule requires 0 + n_int + q = 2*n_int to be odd. Since 2*n_int is
  always even, there are NO allowed vertex triples. This is a genuine
  physical result: the Dirac ground state on S^3 has zero one-loop
  self-energy from the mode sum.

Test structure:
1. Ground-state zero: Sigma(0) = 0 from vertex selection rule
2. Excited states: Sigma(n_ext >= 1) is nonzero and positive
3. Convergence: monotonic increase with n_max
4. State dependence: different n_ext give different values
5. Unrestricted sum factorizes as Dirichlet product
6. PSLQ decomposition structure
7. Vertex rule consistency with qed_vertex module
8. Vertex correction (stretch) computable at small n_max
9. Table function correctness
10. Edge cases and error handling
"""

import mpmath
import pytest

from geovac.qed_self_energy import (
    self_energy_spectral,
    self_energy_table,
    self_energy_convergence,
    self_energy_transcendental_class,
    self_energy_unrestricted,
    vertex_correction_spectral,
    schwinger_convergence,
)


# ---------------------------------------------------------------------------
# Ground state zero: Sigma(n_ext=0) = 0 from vertex selection rule
# ---------------------------------------------------------------------------

class TestGroundStateZero:
    """The n_ext=0 state has zero self-energy due to vertex selection rule.

    Physics: at n_ext=0, the triangle inequality forces q = n_int,
    and the parity rule requires 0 + n_int + q = 2*n_int to be odd.
    Since 2*n_int is always even, there are NO allowed vertex triples,
    and Sigma(0) = 0 identically.
    """

    def test_sigma_0_is_zero(self) -> None:
        """Sigma(n_ext=0) = 0 exactly."""
        sigma = self_energy_spectral(n_ext=0, n_max=20)
        assert sigma == 0, f"Expected Sigma(0)=0, got {float(sigma)}"

    def test_sigma_0_zero_at_small_n_max(self) -> None:
        """Sigma(0) = 0 even at small n_max (structural, not convergence)."""
        for n_max in [1, 3, 5, 10]:
            sigma = self_energy_spectral(n_ext=0, n_max=n_max)
            assert sigma == 0, (
                f"Expected Sigma(0)=0 at n_max={n_max}, got {float(sigma)}")

    def test_no_allowed_triples_at_n_ext_0(self) -> None:
        """No vertex triples are allowed for n_ext=0."""
        from geovac.qed_self_energy import _vertex_allowed

        for n_int in range(20):
            q = n_int  # triangle forces q = n_int
            if q >= 1:
                assert not _vertex_allowed(0, n_int, q), (
                    f"Should NOT be allowed: (0, {n_int}, {q})")

    def test_first_nonzero_state(self) -> None:
        """n_ext=1 is the first external state with nonzero self-energy."""
        sigma_0 = self_energy_spectral(n_ext=0, n_max=10)
        sigma_1 = self_energy_spectral(n_ext=1, n_max=10)
        assert sigma_0 == 0
        assert sigma_1 > 0, f"Expected Sigma(1)>0, got {float(sigma_1)}"


# ---------------------------------------------------------------------------
# Basic properties (excited states)
# ---------------------------------------------------------------------------

class TestSelfEnergyBasic:
    """Core properties of the self-energy spectral sum for excited states."""

    def test_n_ext_1_nonzero(self) -> None:
        """Self-energy at n_ext=1 should be nonzero and positive."""
        sigma = self_energy_spectral(n_ext=1, n_max=5)
        assert sigma > 0, f"Expected positive self-energy at n_ext=1, got {float(sigma)}"

    def test_n_ext_2_nonzero(self) -> None:
        """Self-energy at n_ext=2 should be nonzero and positive."""
        sigma = self_energy_spectral(n_ext=2, n_max=5)
        assert sigma > 0, f"Expected positive self-energy at n_ext=2, got {float(sigma)}"

    def test_n_ext_1_finite(self) -> None:
        """Self-energy at n_ext=1 should be finite (not inf or nan)."""
        sigma = self_energy_spectral(n_ext=1, n_max=5)
        sigma_f = float(sigma)
        assert not (sigma_f == float('inf') or sigma_f == float('-inf')), (
            f"Self-energy is infinite: {sigma_f}")
        assert sigma_f == sigma_f, "Self-energy is NaN"

    def test_real_valued(self) -> None:
        """Self-energy should be a real mpmath.mpf value."""
        sigma = self_energy_spectral(n_ext=1, n_max=5)
        assert isinstance(sigma, mpmath.mpf), (
            f"Expected mpmath.mpf, got {type(sigma)}")

    def test_positive_for_excited_states(self) -> None:
        """Self-energy should be positive for n_ext = 1, 2, 3."""
        for n_ext in [1, 2, 3]:
            sigma = self_energy_spectral(n_ext=n_ext, n_max=5)
            assert sigma > 0, (
                f"Expected positive self-energy at n_ext={n_ext}, got {float(sigma)}")


# ---------------------------------------------------------------------------
# Convergence
# ---------------------------------------------------------------------------

class TestSelfEnergyConvergence:
    """Convergence properties of the self-energy sum."""

    def test_monotonic_convergence_n_ext_1(self) -> None:
        """Self-energy at n_ext=1 should increase monotonically with n_max.

        Each new internal mode adds a non-negative contribution (all terms
        are non-negative: W >= 0, g > 0, d_T > 0, propagators > 0).
        """
        sigmas = []
        for n_max in [2, 4, 6, 8]:
            sigma = self_energy_spectral(n_ext=1, n_max=n_max)
            sigmas.append(float(sigma))

        for i in range(1, len(sigmas)):
            assert sigmas[i] >= sigmas[i - 1], (
                f"Non-monotonic convergence at n_ext=1: "
                f"n_max={2*(i+1)} ({sigmas[i]}) < n_max={2*i} ({sigmas[i-1]})")

    def test_monotonic_convergence_n_ext_2(self) -> None:
        """Same monotonicity test for n_ext=2."""
        sigmas = []
        for n_max in [2, 4, 6, 8]:
            sigma = self_energy_spectral(n_ext=2, n_max=n_max)
            sigmas.append(float(sigma))

        for i in range(1, len(sigmas)):
            assert sigmas[i] >= sigmas[i - 1], (
                f"Non-monotonic at n_ext=2: n_max step {2*i} -> {2*(i+1)}")

    def test_convergence_function(self) -> None:
        """self_energy_convergence returns valid entries."""
        results = self_energy_convergence(
            n_ext=1, n_max_values=[3, 6, 9, 12])
        assert len(results) == 4

        # First entry has no delta
        assert results[0]["delta"] is None

        # All sigmas should be positive
        for r in results:
            assert r["Sigma_float"] > 0

    def test_convergence_deltas_decrease(self) -> None:
        """Successive deltas should decrease (convergence)."""
        results = self_energy_convergence(
            n_ext=1, n_max_values=[5, 10, 15, 20])
        deltas = [r["delta"] for r in results[1:]]

        for i in range(1, len(deltas)):
            assert deltas[i] <= deltas[i - 1] * 1.05, (
                f"Convergence not improving: delta[{i}]={deltas[i]} "
                f"> delta[{i-1}]={deltas[i-1]}")

    def test_convergence_entry_computed(self) -> None:
        """At least one convergence table entry is computable."""
        results = self_energy_convergence(n_ext=1, n_max_values=[3])
        assert len(results) == 1
        assert results[0]["Sigma_float"] > 0


# ---------------------------------------------------------------------------
# State dependence
# ---------------------------------------------------------------------------

class TestSelfEnergyStateDependence:
    """Self-energy should depend on the external state n_ext."""

    def test_different_excited_states_different_sigma(self) -> None:
        """Different excited states should give different self-energies."""
        sigma_1 = self_energy_spectral(n_ext=1, n_max=10)
        sigma_2 = self_energy_spectral(n_ext=2, n_max=10)
        sigma_3 = self_energy_spectral(n_ext=3, n_max=10)

        assert abs(sigma_1 - sigma_2) > 1e-10, (
            f"n_ext=1 and n_ext=2 give same Sigma: {float(sigma_1)}")
        assert abs(sigma_2 - sigma_3) > 1e-10, (
            f"n_ext=2 and n_ext=3 give same Sigma: {float(sigma_2)}")

    def test_ground_state_vs_excited(self) -> None:
        """Ground state (zero) vs excited (nonzero)."""
        sigma_0 = self_energy_spectral(n_ext=0, n_max=10)
        sigma_1 = self_energy_spectral(n_ext=1, n_max=10)
        assert sigma_0 == 0
        assert sigma_1 > 0
        assert sigma_1 > sigma_0


# ---------------------------------------------------------------------------
# Unrestricted sum matches Dirichlet product
# ---------------------------------------------------------------------------

class TestUnrestricted:
    """The unrestricted (no vertex rule) sum should factorize."""

    def test_unrestricted_positive(self) -> None:
        """Unrestricted sum is positive."""
        sigma_unr = self_energy_unrestricted(n_max=10)
        assert sigma_unr > 0

    def test_unrestricted_factorizes(self) -> None:
        """Unrestricted sum = D_electron * D_photon (by construction)."""
        n_max = 10
        sigma_unr = self_energy_unrestricted(n_max=n_max)

        # Verify by direct double sum
        from geovac.qed_self_energy import (
            _g_n_dirac, _lambda_n, _d_q_transverse, _mu_q)

        total_direct = mpmath.mpf(0)
        for n_int in range(n_max + 1):
            for q in range(1, 2 * n_max + 1):
                total_direct += (_g_n_dirac(n_int) * _d_q_transverse(q)
                                 / (_lambda_n(n_int) ** 4 * _mu_q(q)))

        rel_err = float(abs(sigma_unr - total_direct) / abs(sigma_unr))
        assert rel_err < 1e-25, f"Factorization mismatch: rel_err = {rel_err}"

    def test_restricted_leq_unrestricted(self) -> None:
        """Vertex-restricted sum should be <= unrestricted.

        The vertex rule excludes some (n_int, q) pairs, so the
        restricted sum is a subset of the unrestricted sum.
        Use n_ext=1 (nonzero) to verify strict inequality.
        """
        n_max = 10
        sigma_restricted = self_energy_spectral(n_ext=1, n_max=n_max)
        sigma_unr = self_energy_unrestricted(n_max=n_max)

        assert sigma_restricted <= sigma_unr, (
            f"Restricted ({float(sigma_restricted)}) > unrestricted ({float(sigma_unr)})")
        # Strict inequality because vertex rule is non-trivial
        assert sigma_restricted < sigma_unr, (
            f"Restricted equals unrestricted — vertex rule had no effect")


# ---------------------------------------------------------------------------
# Transcendental classification
# ---------------------------------------------------------------------------

class TestTranscendentalClass:
    """PSLQ decomposition of the self-energy."""

    def test_pslq_n_ext_1(self) -> None:
        """PSLQ decomposition at n_ext=1 produces a result dict."""
        result = self_energy_transcendental_class(n_ext=1, n_max=20)
        assert "Sigma_float" in result
        assert result["Sigma_float"] > 0

        decomp = result["decomposition"]
        assert "identified" in decomp

    def test_pslq_result_structure(self) -> None:
        """PSLQ decomposition returns expected dict keys."""
        result = self_energy_transcendental_class(n_ext=1, n_max=15)
        assert "n_ext" in result
        assert "n_max" in result
        assert "decomposition" in result


# ---------------------------------------------------------------------------
# Table function
# ---------------------------------------------------------------------------

class TestSelfEnergyTable:
    """Tests for the multi-state table function."""

    def test_table_correct_length(self) -> None:
        """Table returns correct number of entries."""
        n_ext_values = [0, 1, 2, 3]
        results = self_energy_table(n_ext_values, n_max=5)
        assert len(results) == 4

    def test_table_ground_zero_excited_positive(self) -> None:
        """Ground state entry is zero; excited state entries are positive."""
        results = self_energy_table([0, 1, 2], n_max=5)
        assert results[0]["Sigma_float"] == 0.0
        assert results[1]["Sigma_float"] > 0
        assert results[2]["Sigma_float"] > 0

    def test_table_excited_entries_distinct(self) -> None:
        """Table entries for different excited states should be distinct."""
        results = self_energy_table([1, 2, 3], n_max=8)
        vals = [entry["Sigma_float"] for entry in results]
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                assert abs(vals[i] - vals[j]) > 1e-10, (
                    f"Entries {i} and {j} give same Sigma")


# ---------------------------------------------------------------------------
# Vertex selection rule consistency
# ---------------------------------------------------------------------------

class TestVertexConsistency:
    """Cross-check vertex rules against qed_vertex module."""

    def test_vertex_channel_count_consistency(self) -> None:
        """Our local _so4_channel_count should match qed_vertex."""
        from geovac.qed_vertex import so4_channel_count
        from geovac.qed_self_energy import _so4_channel_count

        for n1 in range(6):
            for n2 in range(6):
                for q in range(1, 12):
                    local = _so4_channel_count(n1, n2, q)
                    ref = so4_channel_count(n1, n2, q)
                    assert local == ref, (
                        f"Channel count mismatch at ({n1},{n2},{q}): "
                        f"local={local}, ref={ref}")

    def test_n_ext_1_has_allowed_triples(self) -> None:
        """At n_ext=1, there should be allowed triples.

        n_ext=1, n_int=1: q in [0,2], q>=1 so q in [1,2].
          q=1: 1+1+1=3 (odd) => YES.
        n_ext=1, n_int=2: q in [1,3].
          q=2: 1+2+2=5 (odd) => YES.
        """
        from geovac.qed_self_energy import _vertex_allowed

        assert _vertex_allowed(1, 1, 1), "Should be allowed: (1, 1, 1)"
        assert _vertex_allowed(1, 2, 2), "Should be allowed: (1, 2, 2)"
        assert not _vertex_allowed(1, 0, 1), (
            "Should NOT be allowed: (1, 0, 1), parity 1+0+1=2 even")


# ---------------------------------------------------------------------------
# Vertex correction (stretch)
# ---------------------------------------------------------------------------

class TestVertexCorrection:
    """Tests for the stretch-goal vertex correction."""

    def test_vertex_correction_computable(self) -> None:
        """Vertex correction at small n_max is computable."""
        lam = vertex_correction_spectral(n_ext=1, n_max=3)
        assert isinstance(lam, mpmath.mpf)
        assert lam >= 0

    def test_vertex_correction_positive_n_ext_1(self) -> None:
        """Vertex correction at n_ext=1 should be positive."""
        lam = vertex_correction_spectral(n_ext=1, n_max=5)
        assert lam > 0, f"Expected positive vertex correction, got {float(lam)}"

    def test_vertex_correction_zero_at_n_ext_0(self) -> None:
        """Vertex correction at n_ext=0 should be zero (same vertex rule)."""
        lam = vertex_correction_spectral(n_ext=0, n_max=5)
        assert lam == 0, f"Expected zero vertex correction at n_ext=0, got {float(lam)}"

    def test_schwinger_convergence_runs(self) -> None:
        """schwinger_convergence produces results at small n_max."""
        results = schwinger_convergence([2, 3, 4], n_ext=1)
        assert len(results) == 3
        assert results[0]["Lambda_float"] >= 0


# ---------------------------------------------------------------------------
# Edge cases and error handling
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge cases and error handling."""

    def test_n_ext_negative_raises(self) -> None:
        """Negative n_ext should raise ValueError."""
        with pytest.raises(ValueError):
            self_energy_spectral(n_ext=-1, n_max=5)

    def test_n_max_zero(self) -> None:
        """n_max=0 should be computable (may be zero)."""
        sigma = self_energy_spectral(n_ext=1, n_max=0)
        assert isinstance(sigma, mpmath.mpf)

    def test_large_n_ext_small_n_max(self) -> None:
        """Large n_ext with small n_max: few or no allowed triples."""
        sigma = self_energy_spectral(n_ext=10, n_max=2)
        assert isinstance(sigma, mpmath.mpf)
        assert sigma >= 0

    def test_n_max_negative_raises(self) -> None:
        """Negative n_max should raise ValueError."""
        with pytest.raises(ValueError):
            self_energy_spectral(n_ext=1, n_max=-1)

    def test_equal_n_ext_n_max(self) -> None:
        """n_ext = n_max should be computable."""
        sigma = self_energy_spectral(n_ext=3, n_max=3)
        assert isinstance(sigma, mpmath.mpf)
        assert sigma >= 0
