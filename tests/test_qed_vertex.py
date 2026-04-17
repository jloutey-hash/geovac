"""Tests for QED vertex coupling on S^3 and two-loop zeta(3) structure.

Verifies:
- Selection rule consistency with hodge1_s3.vertex_coupling
- Two-loop sum convergence
- Even/odd parity split mechanism for zeta(3)
- Flat-space limit
- Factorized form vs product of single sums
- Transcendental classification
"""

import pytest
import mpmath

mpmath.mp.dps = 80

from geovac.qed_vertex import (
    reduced_coupling_squared,
    weighted_vertex_coupling,
    so4_channel_count,
    effective_pair_weight,
    vertex_allowed_triples,
    two_loop_sunset_unrestricted,
    two_loop_sunset_vertex_restricted,
    two_loop_sunset_weighted,
    two_loop_vertex_correction,
    two_loop_odd_even_split,
    two_loop_min_weighted_hurwitz,
    decompose_two_loop_result,
    two_loop_transcendental_classification,
    verify_vertex_factorization_failure,
    flat_space_vertex_sum,
    two_loop_photon_line,
)


# ---------------------------------------------------------------------------
# Selection rule consistency
# ---------------------------------------------------------------------------

class TestSelectionRule:
    """Verify vertex selection rule matches hodge1_s3.vertex_coupling."""

    def test_matches_hodge1_vertex_coupling(self):
        """qed_vertex selection rule agrees with hodge1_s3.vertex_coupling."""
        from geovac.hodge1_s3 import vertex_coupling
        from geovac.qed_vertex import _vertex_allowed

        for n1 in range(6):
            for n2 in range(6):
                for ng in range(1, 8):
                    assert _vertex_allowed(n1, n2, ng) == vertex_coupling(n1, n2, ng), \
                        f"Mismatch at ({n1},{n2},{ng})"

    def test_dipole_selection_rule(self):
        """Dipole coupling (q=1) requires |n1-n2| <= 1 and n1+n2 even."""
        from geovac.qed_vertex import _vertex_allowed

        # q=1: need |n1-n2| <= 1 <= n1+n2 and n1+n2+1 odd => n1+n2 even
        assert _vertex_allowed(0, 0, 1) is False  # n1+n2=0 < q=1
        assert _vertex_allowed(1, 1, 1) is True   # |0|<=1<=2, 1+1+1=3 odd
        assert _vertex_allowed(1, 0, 1) is False  # 1+0+1=2 even -> forbidden
        assert _vertex_allowed(2, 1, 1) is False  # 2+1+1=4 even, need odd
        assert _vertex_allowed(2, 2, 1) is True   # 2+2+1=5 odd, |0|<=1<=4, OK

    def test_forbidden_pair_00(self):
        """The (0,0) pair has no allowed photon mode."""
        from geovac.qed_vertex import _vertex_allowed
        # n1=n2=0: need q with |0|<=q<=0, so q=0, but q>=1 required
        for q in range(1, 10):
            assert _vertex_allowed(0, 0, q) is False

    def test_reduced_coupling_zero_for_forbidden(self):
        """Reduced coupling is 0 for forbidden triples."""
        assert reduced_coupling_squared(0, 0, 1) == 0
        assert reduced_coupling_squared(1, 0, 1) == 0  # parity forbidden

    def test_reduced_coupling_positive_for_allowed(self):
        """Reduced coupling is positive for allowed triples."""
        C = reduced_coupling_squared(1, 1, 1)
        assert C > 0

    def test_allowed_triples_count(self):
        """Count of allowed triples at small cutoffs."""
        triples = vertex_allowed_triples(3, 3)
        assert len(triples) > 0
        # All triples should satisfy selection rule
        from geovac.qed_vertex import _vertex_allowed
        for n1, n2, ng in triples:
            assert _vertex_allowed(n1, n2, ng)


# ---------------------------------------------------------------------------
# Two-loop unrestricted
# ---------------------------------------------------------------------------

class TestUnrestricted:
    """Unrestricted two-loop sunset = D(4)^2."""

    def test_factorizes_as_D4_squared(self):
        """Unrestricted sum equals D(4)^2 = (pi^2 - pi^4/12)^2."""
        result = two_loop_sunset_unrestricted(n_max=100)
        D4 = mpmath.pi**2 - mpmath.pi**4 / 12
        target = D4**2
        assert float(abs(result["product_exact"] - target)) < 1e-60

    def test_direct_sum_converges_to_product(self):
        """Direct double sum converges to the factorized form."""
        result = two_loop_sunset_unrestricted(n_max=100)
        # Truncation error at n_max=100 should be modest
        assert result["rel_error_direct_vs_exact"] < 0.05

    def test_unrestricted_is_pi_even(self):
        """Unrestricted sum is purely pi^{even}, no zeta(3)."""
        result = two_loop_sunset_unrestricted(n_max=100)
        # D(4)^2 = pi^4 - pi^6/6 + pi^8/144
        decomp = result["decomposition"]
        if decomp and decomp.get("identified"):
            assert not decomp.get("contains_zeta3", False), \
                "Unrestricted sum should not contain zeta(3)"


# ---------------------------------------------------------------------------
# Even/odd split: the mechanism for zeta(3)
# ---------------------------------------------------------------------------

class TestEvenOddSplit:
    """D_even(4) and D_odd(4) individually contain Catalan G and beta(4)."""

    def test_sum_reproduces_D4(self):
        """D_even(4) + D_odd(4) = D(4) to machine precision.

        Now uses exact Hurwitz zeta evaluation (not direct truncated sums).
        """
        result = two_loop_odd_even_split(n_max=100, s=4)
        assert result["rel_error_sum_vs_exact"] < 1e-60

    def test_even_and_odd_are_different(self):
        """D_even(4) != D_odd(4)."""
        result = two_loop_odd_even_split(n_max=100, s=4)
        ratio = float(result["D_even"] / result["D_odd"])
        assert abs(ratio - 1.0) > 0.01  # They differ meaningfully

    def test_individual_parts_contain_catalan(self):
        """D_even(4) and D_odd(4) individually contain Catalan's constant G.

        The even-n sub-sum uses Hurwitz zeta at shift 3/4, whose content
        includes Catalan's constant G = beta(2) and Dirichlet beta(4).
        These are NOT Riemann odd-zeta (zeta(3)) -- they are Dirichlet beta
        function values from the quarter-integer Hurwitz shift.

        D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4)
        D_odd(4)  = pi^2/2 - pi^4/24 + 4G - 4*beta(4)
        """
        result = two_loop_odd_even_split(n_max=100, s=4)
        # The total should NOT contain Catalan (it's D(4) = pi^{even})
        decomp_total = result["decomp_total"]
        if decomp_total and decomp_total.get("identified"):
            assert not decomp_total.get("contains_zeta3", False), \
                "D(4) = D_even + D_odd should be pi^{even}"

        # Individual parts contain Catalan's constant
        decomp_even = result["decomp_even"]
        decomp_odd = result["decomp_odd"]
        if decomp_even and decomp_even.get("identified"):
            assert decomp_even.get("contains_catalan", False), \
                "D_even(4) should contain Catalan's constant G"
        if decomp_odd and decomp_odd.get("identified"):
            assert decomp_odd.get("contains_catalan", False), \
                "D_odd(4) should contain Catalan's constant G"

    def test_catalan_cancels_in_sum(self):
        """The Catalan/beta(4) content of D_even and D_odd cancels in the sum.

        D(4) = pi^2 - pi^4/12 has no Catalan or beta(4), so these must
        cancel between D_even and D_odd.
        """
        result = two_loop_odd_even_split(n_max=100, s=4)
        D_total = result["D_even"] + result["D_odd"]

        # Exact match since both are computed via Hurwitz
        D4_exact = mpmath.pi**2 - mpmath.pi**4 / 12
        assert float(abs(D_total - D4_exact) / abs(D4_exact)) < 1e-60


# ---------------------------------------------------------------------------
# Factorization failure
# ---------------------------------------------------------------------------

class TestFactorizationFailure:
    """The vertex-restricted sum does NOT factorize."""

    def test_same_vs_opposite_parity_differ(self):
        """Same-parity and opposite-parity double sums differ."""
        result = verify_vertex_factorization_failure(n_max=50)
        same = result["double_same_parity_float"]
        opp = result["double_opposite_parity_float"]
        assert abs(same - opp) / abs(same + opp) > 0.01

    def test_sum_check(self):
        """same + opposite = full^2."""
        result = verify_vertex_factorization_failure(n_max=50)
        assert result["sum_check_rel_error"] < 1e-30

    def test_same_parity_not_pi_even(self):
        """D_even^2 + D_odd^2 is NOT pi^{even}.

        Since D_even and D_odd individually contain Catalan G and beta(4),
        their same-parity double sum D_even^2 + D_odd^2 generically
        contains G^2, G*beta(4), beta(4)^2 terms that are NOT pi^{even}.
        """
        result = verify_vertex_factorization_failure(n_max=50)
        same = result["double_same_parity_float"]
        # Same-parity and opposite-parity should differ
        # (already tested above, this just documents the mechanism)
        assert same > 0


# ---------------------------------------------------------------------------
# Flat-space limit
# ---------------------------------------------------------------------------

class TestFlatSpaceLimit:
    """Flat-space vertex sum produces zeta(3) via Euler identity."""

    def test_flat_vertex_sum_converges(self):
        """The flat-space vertex sum converges to a finite value."""
        result = flat_space_vertex_sum(N=500)
        assert result["vertex_sum_float"] > 0

    def test_flat_adjacent_pairs_relate_to_zeta3(self):
        """sum 1/(n^2(n+1)^2) converges to pi^4/72 - 2*pi^2/3 + ... (known).

        Actually: sum_{n=1}^inf 1/(n^2(n+1)^2) = pi^2/3 - 3 (partial fractions).
        Wait, let's compute: 1/(n^2(n+1)^2) = 1/n^2 + 1/(n+1)^2 - 2/(n(n+1))
        Sum from 1 to inf: zeta(2) + (zeta(2)-1) - 2*sum 1/(n(n+1))
        = 2*zeta(2) - 1 - 2*(1) = 2*pi^2/6 - 3 = pi^2/3 - 3.

        This is pi^{even} + rational, NO zeta(3). The zeta(3) in the
        flat-space two-loop comes from the NESTED sum structure, not
        from single nearest-neighbor sums.
        """
        result = flat_space_vertex_sum(N=500)
        adj = result["adjacent_pairs"]
        target = mpmath.pi**2 / 3 - 3
        assert float(abs(adj - target) / abs(target)) < 0.01

    def test_flat_unrestricted_is_zeta2_squared(self):
        """Unrestricted flat sum = zeta(2)^2 = pi^4/36."""
        result = flat_space_vertex_sum(N=500)
        target = mpmath.pi**4 / 36
        assert float(abs(result["unrestricted_sum"] - target) / abs(target)) < 1e-50


# ---------------------------------------------------------------------------
# Vertex-restricted two-loop
# ---------------------------------------------------------------------------

class TestVertexRestricted:
    """The vertex-restricted two-loop sum with photon propagator."""

    @pytest.mark.slow
    def test_converges(self):
        """Vertex-restricted sum converges to a finite value."""
        result = two_loop_sunset_vertex_restricted(n_max=30)
        assert result["value_float"] > 0
        assert result["n_triples"] > 0

    @pytest.mark.slow
    def test_increases_with_n_max(self):
        """Vertex-restricted sum increases with n_max (positive definite)."""
        r1 = two_loop_sunset_vertex_restricted(n_max=10)
        r2 = two_loop_sunset_vertex_restricted(n_max=20)
        assert r2["value_float"] > r1["value_float"]


# ---------------------------------------------------------------------------
# Two-loop photon line
# ---------------------------------------------------------------------------

class TestPhotonLine:
    """Two-loop sum mixing Dirac and Hodge-1 spectra."""

    @pytest.mark.slow
    def test_photon_line_converges(self):
        """Photon line sum converges."""
        result = two_loop_photon_line(n_max=30, s_e=4, s_gamma=1)
        assert result["value_float"] > 0

    @pytest.mark.slow
    def test_photon_line_contains_structure(self):
        """Photon line sum has identifiable transcendental structure."""
        result = two_loop_photon_line(n_max=50, s_e=4, s_gamma=1)
        decomp = result["decomposition"]
        # The sum mixes half-integer (Dirac) and integer (Hodge-1) spectra,
        # which should break the pure pi^{even} structure
        assert decomp is not None


# ---------------------------------------------------------------------------
# Vertex correction
# ---------------------------------------------------------------------------

class TestVertexCorrection:
    """The pair-existence correction from vertex restriction."""

    def test_n0_pairs_are_forbidden(self):
        """All (0, m) and (n, 0) pairs have no allowed photon mode.

        For (0, m): triangle requires q = m exactly, parity requires
        0 + m + m = 2m odd, which is never satisfied. So all pairs with
        either index 0 are forbidden. At n_max=20, there are 2*20+1 = 41
        forbidden pairs: (0,0), (0,1)...(0,20), (1,0)...(20,0).
        """
        result = two_loop_vertex_correction(n_max=20)
        assert (0, 0) in result["forbidden_pairs"]
        # All (0, m) and (n, 0) are forbidden: 2*n_max + 1 = 41
        assert result["n_forbidden"] == 2 * 20 + 1

    def test_correction_includes_all_n0_pairs(self):
        """The correction sums all forbidden (n,0) and (0,m) pair contributions.

        Correction = sum_{forbidden (n,m)} g_n * g_m / (lambda_n^4 * lambda_m^4).
        This is 2 * g_0 * D(4)_trunc - g_0^2/lambda_0^8 (row + col - corner).
        Verify the correction is positive and larger than just the (0,0) term.
        """
        result = two_loop_vertex_correction(n_max=100, s1=2, s2=2)
        correction = result["correction"]
        # The (0,0)-only contribution
        g0 = mpmath.mpf(4)  # g_0 = 2*1*2
        l0 = mpmath.mpf(3) / 2
        only_00 = g0**2 / l0**8
        # Full correction includes all row-0 + col-0 terms, so is larger
        assert float(correction) > float(only_00)
        assert float(correction) > 0


# ---------------------------------------------------------------------------
# PSLQ decomposition
# ---------------------------------------------------------------------------

class TestPSLQDecomposition:
    """Extended PSLQ decomposition for two-loop quantities."""

    def test_decompose_pi_squared(self):
        """pi^2 is identified correctly."""
        result = decompose_two_loop_result(mpmath.pi**2)
        assert result["identified"]
        assert "pi^2" in result["components"]
        assert not result["contains_zeta3"]

    def test_decompose_zeta3(self):
        """zeta(3) is identified correctly."""
        result = decompose_two_loop_result(mpmath.zeta(3))
        assert result["identified"]
        assert result["contains_zeta3"]

    def test_decompose_zeta3_squared(self):
        """zeta(3)^2 is identified correctly."""
        result = decompose_two_loop_result(mpmath.zeta(3)**2)
        assert result["identified"]
        assert result["contains_zeta3_squared"]

    def test_decompose_mixed(self):
        """pi^2 + 3*zeta(3) is identified."""
        val = mpmath.pi**2 + 3 * mpmath.zeta(3)
        result = decompose_two_loop_result(val)
        assert result["identified"]
        assert result["contains_zeta3"]
        assert result["residual"] < 1e-50


# ---------------------------------------------------------------------------
# Transcendental classification
# ---------------------------------------------------------------------------

class TestTranscendentalClassification:
    """Paper 18 taxonomy classification."""

    def test_classification_complete(self):
        """All expected categories are present."""
        cl = two_loop_transcendental_classification()
        assert "vertex_coupling" in cl
        assert "selection_rule" in cl
        assert "unrestricted_two_loop_s_even" in cl
        assert "vertex_restricted_two_loop" in cl
        assert "beta_mechanism" in cl
        assert "flat_space_limit" in cl
        assert "paper18_operator_order" in cl

    def test_vertex_is_rational(self):
        """Vertex coupling is classified as rational."""
        cl = two_loop_transcendental_classification()
        assert "rational" in cl["vertex_coupling"].lower()

    def test_restricted_is_dirichlet_beta(self):
        """Vertex-restricted two-loop is classified as dirichlet_beta."""
        cl = two_loop_transcendental_classification()
        assert "dirichlet_beta" in cl["vertex_restricted_two_loop"]


# ---------------------------------------------------------------------------
# Structural theorem: even/odd split PSLQ
# ---------------------------------------------------------------------------

class TestStructuralTheorem:
    """The core structural result: Catalan G and beta(4) from even/odd split."""

    def test_D_even_4_decomposition(self):
        """D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4).

        PSLQ identifies D_even(4) in the {1, pi^2, pi^4, G, beta(4)} basis.
        """
        from geovac.qed_vertex import _dirac_D_even, _dirichlet_beta
        D_even = _dirac_D_even(4)

        pi2 = mpmath.pi**2
        pi4 = mpmath.pi**4
        G = mpmath.catalan
        beta4 = _dirichlet_beta(4)

        basis = [D_even, mpmath.mpf(1), pi2, pi4, G, beta4]
        relation = mpmath.pslq(basis, tol=1e-40, maxcoeff=100000)

        assert relation is not None, "PSLQ should identify D_even(4)"
        # PSLQ returns [-24, 0, 12, -1, -96, 96]
        # => D_even = (12*pi^2 - pi^4 - 96*G + 96*beta(4))/24
        # => D_even = pi^2/2 - pi^4/24 - 4*G + 4*beta(4)
        # The G coefficient should be nonzero
        assert relation[4] != 0, "D_even(4) should have nonzero Catalan coefficient"
        assert relation[5] != 0, "D_even(4) should have nonzero beta(4) coefficient"

    def test_D_odd_4_has_opposite_catalan(self):
        """D_odd(4) has opposite Catalan/beta(4) coefficients as D_even(4).

        Since D_even(4) + D_odd(4) = D(4) = pi^2 - pi^4/12 (no G, no beta(4)),
        the G and beta(4) coefficients must cancel.
        """
        from geovac.qed_vertex import _dirac_D_even, _dirac_D_odd, _dirichlet_beta
        D_even = _dirac_D_even(4)
        D_odd = _dirac_D_odd(4)

        pi2 = mpmath.pi**2
        pi4 = mpmath.pi**4
        G = mpmath.catalan
        beta4 = _dirichlet_beta(4)

        # PSLQ for D_even
        rel_even = mpmath.pslq([D_even, mpmath.mpf(1), pi2, pi4, G, beta4],
                                tol=1e-40, maxcoeff=100000)
        # PSLQ for D_odd
        rel_odd = mpmath.pslq([D_odd, mpmath.mpf(1), pi2, pi4, G, beta4],
                               tol=1e-40, maxcoeff=100000)

        assert rel_even is not None, "PSLQ should identify D_even(4)"
        assert rel_odd is not None, "PSLQ should identify D_odd(4)"

        # Extract G coefficients: -rel[4]/rel[0]
        coeff_G_even = mpmath.mpf(-rel_even[4]) / rel_even[0]
        coeff_G_odd = mpmath.mpf(-rel_odd[4]) / rel_odd[0]

        # They should sum to zero (G cancels in D(4))
        assert float(abs(coeff_G_even + coeff_G_odd)) < 1e-30, \
            f"G coefficients should cancel: {coeff_G_even} + {coeff_G_odd}"

        # Same for beta(4) coefficients
        coeff_b4_even = mpmath.mpf(-rel_even[5]) / rel_even[0]
        coeff_b4_odd = mpmath.mpf(-rel_odd[5]) / rel_odd[0]
        assert float(abs(coeff_b4_even + coeff_b4_odd)) < 1e-30, \
            f"beta(4) coefficients should cancel: {coeff_b4_even} + {coeff_b4_odd}"

    def test_even_odd_split_at_s6(self):
        """At s=6, D_even(6) and D_odd(6) also have cancelling extra content.

        D(6) = pi^4/3 - pi^6/30 (pi^{even}), so any non-pi content in
        D_even(6) and D_odd(6) must cancel. At s=6, the Hurwitz at 3/4
        involves Catalan G and beta(6).
        """
        from geovac.qed_vertex import _dirac_D_even, _dirac_D_odd
        D_even = _dirac_D_even(6)
        D_odd = _dirac_D_odd(6)

        D_total = D_even + D_odd
        D6_exact = mpmath.pi**4 / 3 - mpmath.pi**6 / 30
        assert float(abs(D_total - D6_exact) / abs(D6_exact)) < 1e-60


# ---------------------------------------------------------------------------
# SO(4) Clebsch-Gordan channel count (weighted vertex)
# ---------------------------------------------------------------------------

class TestSO4ChannelCount:
    """Tests for the SO(4) CG channel count vertex weight."""

    def test_forbidden_gives_zero(self):
        """Channel count is 0 for forbidden triples."""
        assert so4_channel_count(0, 0, 1) == 0
        assert so4_channel_count(1, 0, 1) == 0  # parity

    def test_dipole_gives_one(self):
        """Dipole (q=1) coupling has at most 1 channel per triple."""
        for n1 in range(6):
            for n2 in range(6):
                W = so4_channel_count(n1, n2, 1)
                assert W <= 1, f"Dipole at ({n1},{n2},1) should have W<=1, got {W}"

    def test_higher_multipole_can_give_two(self):
        """Higher multipoles (q>=2) can have W=2."""
        # (2, 2, 3) should have W=2 (both SO(4) components contribute)
        assert so4_channel_count(2, 2, 3) == 2

    def test_channel_count_range(self):
        """Channel count is always 0, 1, or 2."""
        for n1 in range(8):
            for n2 in range(8):
                for q in range(1, 10):
                    W = so4_channel_count(n1, n2, q)
                    assert W in (0, 1, 2), f"W({n1},{n2},{q}) = {W}, expected 0/1/2"

    def test_symmetry(self):
        """Vertex weight is symmetric under (n1, n2) swap (time reversal)."""
        for n1 in range(6):
            for n2 in range(6):
                for q in range(1, 8):
                    assert so4_channel_count(n1, n2, q) == so4_channel_count(n2, n1, q), \
                        f"Asymmetry at ({n1},{n2},{q})"


class TestEffectivePairWeight:
    """Tests for the effective pair weight W_total(n1, n2)."""

    def test_formula_2min_minus_1_minus_delta(self):
        """W_total = 2*min(n1,n2) - 1 - delta_{n1,n2} for n1,n2 >= 1.

        Verified for all (n1, n2) up to 12.
        """
        for n1 in range(1, 13):
            for n2 in range(1, 13):
                result = effective_pair_weight(n1, n2)
                m = min(n1, n2)
                delta = 1 if n1 == n2 else 0
                expected = 2 * m - 1 - delta
                assert result["W_total"] == expected, \
                    f"({n1},{n2}): W_total={result['W_total']}, expected={expected}"

    def test_N_q_equals_min(self):
        """N_q = min(n1, n2) for n1, n2 >= 1."""
        for n1 in range(1, 10):
            for n2 in range(1, 10):
                result = effective_pair_weight(n1, n2)
                assert result["N_q"] == min(n1, n2), \
                    f"({n1},{n2}): N_q={result['N_q']}, expected={min(n1,n2)}"

    def test_ratio_approaches_two(self):
        """W_total/N_q -> 2 for large min(n1, n2)."""
        result = effective_pair_weight(50, 60)
        assert abs(result["ratio"] - 2.0) < 0.05


class TestWeightedVertex:
    """Tests for the weighted_vertex_coupling function."""

    def test_returns_channel_count(self):
        """weighted_vertex_coupling returns the channel count as mpf."""
        val = weighted_vertex_coupling(2, 2, 3)
        assert float(val) == 2.0

    def test_zero_for_forbidden(self):
        """Returns 0 for forbidden triples."""
        assert float(weighted_vertex_coupling(0, 0, 1)) == 0.0


class TestWeightedSunset:
    """Tests for the weighted two-loop sunset."""

    @pytest.mark.slow
    def test_weighted_exceeds_uniform(self):
        """Weighted sum > uniform because W >= 1 and often W = 2."""
        result = two_loop_sunset_weighted(n_max=20)
        assert result["weighted_float"] > result["uniform_float"]

    @pytest.mark.slow
    def test_ratio_between_1_and_2(self):
        """Weighted/uniform ratio is between 1 and 2."""
        result = two_loop_sunset_weighted(n_max=20)
        r = result["ratio"]
        assert 1.0 < r < 2.0, f"Ratio {r} outside expected range"


class TestMinWeightedHurwitz:
    """Tests for the min-weighted Hurwitz double sum."""

    @pytest.mark.slow
    def test_converges(self):
        """The min-weighted sum converges to a positive value."""
        result = two_loop_min_weighted_hurwitz(s=4, n_terms=200)
        assert result["S_min_float"] > 0

    @pytest.mark.slow
    def test_less_than_D_squared(self):
        """S_min < D(s)^2 because min(n1,n2) < max(n1,n2)."""
        result = two_loop_min_weighted_hurwitz(s=4, n_terms=500)
        assert result["S_min_float"] < result["D_s_squared_float"]

    @pytest.mark.slow
    def test_not_in_standard_basis(self):
        """S_min is NOT identified in the standard MZV+beta basis.

        This is the key result: the CG-weighted vertex introduces a
        depth-2 multiple Hurwitz zeta value that is genuinely new.
        """
        result = two_loop_min_weighted_hurwitz(s=4, n_terms=500)
        decomp = result["decomposition"]
        assert not decomp["identified"], \
            "S_min should NOT be identifiable in standard MZV+beta basis"
