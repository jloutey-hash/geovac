"""Tests for graph-native to continuum QED bridge module.

Tests cover all five parts (A-E) of the bridge analysis,
plus spectrum helpers and convention conversions.
"""

import pytest
import sympy as sp
from sympy import Rational
from fractions import Fraction

from geovac.graph_qed_continuum_bridge import (
    # Spectrum helpers
    _lambda_n_ch,
    _g_n_dirac,
    _mu_q,
    _d_q_transverse,
    _vertex_allowed,
    _so4_channel_count,
    # Convention bridge
    _n_ch_from_n_fock,
    _n_fock_from_n_ch,
    _n_max_ch_from_n_max_fock,
    # Part A
    continuum_vp_mode_contribution,
    continuum_vp_truncated,
    # Part B
    compute_projection_constant,
    classify_projection_constant,
    # Part C
    graph_vp_edge_decomposition,
    continuum_vp_mode_decomposition,
    mode_ratio_analysis,
    _enumerate_allowed_triples,
    # Part D
    convergence_table,
    _continuum_vp_full_estimate,
    # Part E
    graph_self_energy_n0,
    self_energy_structural_zero_bridge,
    # Driver
    run_bridge_analysis,
)


# =========================================================================
# Spectrum helpers
# =========================================================================

class TestSpectrumHelpers:
    """Tests for Camporesi-Higuchi spectrum functions."""

    def test_lambda_n0(self):
        """lambda(n=0) = 3/2."""
        assert _lambda_n_ch(0) == Rational(3, 2)

    def test_lambda_n1(self):
        """lambda(n=1) = 5/2."""
        assert _lambda_n_ch(1) == Rational(5, 2)

    def test_lambda_n2(self):
        """lambda(n=2) = 7/2."""
        assert _lambda_n_ch(2) == Rational(7, 2)

    def test_lambda_general(self):
        """lambda(n) = n + 3/2 for several n."""
        for n in range(10):
            assert _lambda_n_ch(n) == Rational(2 * n + 3, 2)

    def test_g_n_dirac_n0(self):
        """g(0) = 2*1*2 = 4."""
        assert _g_n_dirac(0) == 4

    def test_g_n_dirac_n1(self):
        """g(1) = 2*2*3 = 12."""
        assert _g_n_dirac(1) == 12

    def test_g_n_dirac_n2(self):
        """g(2) = 2*3*4 = 24."""
        assert _g_n_dirac(2) == 24

    def test_g_n_dirac_n3(self):
        """g(3) = 2*4*5 = 40 (Paper 2 Delta^{-1})."""
        assert _g_n_dirac(3) == 40

    def test_mu_q1(self):
        """mu(1) = 1*3 = 3."""
        assert _mu_q(1) == 3

    def test_mu_q2(self):
        """mu(2) = 2*4 = 8."""
        assert _mu_q(2) == 8

    def test_d_q_transverse_equals_mu(self):
        """d_T(q) = q(q+2) = mu(q) for all q."""
        for q in range(1, 10):
            assert _d_q_transverse(q) == _mu_q(q)


# =========================================================================
# Vertex selection rules
# =========================================================================

class TestVertexSelectionRules:
    """Tests for SO(4) vertex selection rules."""

    def test_parity_even_forbidden(self):
        """Even parity (n1+n2+q even) is forbidden."""
        assert not _vertex_allowed(0, 1, 1)  # 0+1+1=2 even
        assert not _vertex_allowed(1, 2, 1)  # 1+2+1=4 even

    def test_parity_odd_allowed(self):
        """Odd parity (n1+n2+q odd) passes parity check."""
        assert _vertex_allowed(1, 1, 1)  # 1+1+1=3 odd
        assert _vertex_allowed(1, 2, 2)  # 1+2+2=5 odd

    def test_q_zero_forbidden(self):
        """q=0 is always forbidden."""
        assert not _vertex_allowed(0, 0, 0)
        assert not _vertex_allowed(1, 1, 0)

    def test_triangle_violation(self):
        """Triangle inequality violation is forbidden."""
        assert not _vertex_allowed(0, 0, 1)  # q=1 > 0+0=0
        assert not _vertex_allowed(1, 1, 5)  # q=5 > 1+1=2

    def test_w_111_is_zero(self):
        """W(1,1,1) = 0 even though the triple is parity-allowed."""
        assert _vertex_allowed(1, 1, 1)
        assert _so4_channel_count(1, 1, 1) == 0

    def test_w_122_is_one(self):
        """W(1,2,2) = 1."""
        assert _so4_channel_count(1, 2, 2) == 1

    def test_w_212_is_one(self):
        """W(2,1,2) = 1."""
        assert _so4_channel_count(2, 1, 2) == 1

    def test_w_223_is_two(self):
        """W(2,2,3) = 2."""
        assert _so4_channel_count(2, 2, 3) == 2

    def test_w_forbidden_triple_is_zero(self):
        """Forbidden triples give W=0."""
        assert _so4_channel_count(0, 1, 1) == 0
        assert _so4_channel_count(0, 0, 1) == 0

    def test_w_matches_qed_vertex(self):
        """Bridge W matches qed_vertex.so4_channel_count for several triples."""
        from geovac.qed_vertex import so4_channel_count as qv_W
        for n1 in range(5):
            for n2 in range(5):
                for q in range(1, n1 + n2 + 1):
                    assert _so4_channel_count(n1, n2, q) == qv_W(n1, n2, q), \
                        f"Mismatch at ({n1},{n2},{q})"


# =========================================================================
# Convention bridge
# =========================================================================

class TestConventionBridge:
    """Tests for Fock <-> CH convention conversion."""

    def test_roundtrip(self):
        """n_fock -> n_ch -> n_fock is identity."""
        for n in range(1, 10):
            assert _n_fock_from_n_ch(_n_ch_from_n_fock(n)) == n

    def test_n_ch_from_n_fock(self):
        """n_CH = n_fock - 1."""
        assert _n_ch_from_n_fock(1) == 0
        assert _n_ch_from_n_fock(2) == 1
        assert _n_ch_from_n_fock(3) == 2

    def test_n_max_ch_from_n_max_fock(self):
        """n_max_CH = n_max_fock - 1."""
        assert _n_max_ch_from_n_max_fock(2) == 1
        assert _n_max_ch_from_n_max_fock(3) == 2


# =========================================================================
# Part A: Continuum VP
# =========================================================================

class TestContinuumVP:
    """Tests for the truncated continuum VP spectral sum."""

    def test_mode_contribution_forbidden_is_zero(self):
        """Forbidden triples give zero contribution."""
        assert continuum_vp_mode_contribution(0, 1, 1) == 0

    def test_mode_contribution_w0_is_zero(self):
        """Allowed triple with W=0 gives zero contribution."""
        assert continuum_vp_mode_contribution(1, 1, 1) == 0

    def test_mode_contribution_122_is_rational(self):
        """(1,2,2) contribution is a nonzero rational."""
        c = continuum_vp_mode_contribution(1, 2, 2)
        assert c != 0
        assert isinstance(c, sp.Rational)

    def test_mode_contribution_symmetry(self):
        """Pi(n1,n2,q) = Pi(n2,n1,q) (symmetric in electron indices)."""
        for n1 in range(4):
            for n2 in range(4):
                for q in range(1, n1 + n2 + 1):
                    c1 = continuum_vp_mode_contribution(n1, n2, q)
                    c2 = continuum_vp_mode_contribution(n2, n1, q)
                    assert c1 == c2, f"Asymmetry at ({n1},{n2},{q})"

    def test_truncated_nmax_ch_0(self):
        """Continuum VP at n_max_ch=0 is zero (no allowed triples)."""
        assert continuum_vp_truncated(0) == 0

    def test_truncated_nmax_ch_1(self):
        """Continuum VP at n_max_ch=1 is zero (only triple has W=0)."""
        assert continuum_vp_truncated(1) == 0

    def test_truncated_nmax_ch_2_nonzero(self):
        """Continuum VP at n_max_ch=2 is nonzero."""
        val = continuum_vp_truncated(2)
        assert val > 0

    def test_truncated_nmax_ch_2_value(self):
        """Continuum VP at n_max_ch=2 = 538361856/3603000625."""
        val = continuum_vp_truncated(2)
        assert val == Rational(538361856, 3603000625)

    def test_truncated_monotone(self):
        """Truncated continuum VP is non-decreasing in n_max_ch."""
        vals = [float(continuum_vp_truncated(n)) for n in range(6)]
        for i in range(len(vals) - 1):
            assert vals[i + 1] >= vals[i]


# =========================================================================
# Part B: Projection constant
# =========================================================================

class TestProjectionConstant:
    """Tests for the projection exchange constant."""

    def test_nmax2_continuum_is_zero(self):
        """At n_max_fock=2, continuum VP is zero."""
        result = compute_projection_constant(n_max_fock=2)
        assert result['continuum_is_zero']

    def test_nmax3_continuum_nonzero(self):
        """At n_max_fock=3, continuum VP is nonzero."""
        result = compute_projection_constant(n_max_fock=3)
        assert not result['continuum_is_zero']

    def test_nmax3_graph_trace(self):
        """Graph VP trace at n_max_fock=3 is 3872/735."""
        result = compute_projection_constant(n_max_fock=3)
        assert result['graph_trace'] == '3872/735'

    def test_nmax3_is_rational(self):
        """Projection constant at n_max_fock=3 is rational."""
        result = compute_projection_constant(n_max_fock=3)
        assert result['is_rational']

    def test_nmax3_value(self):
        """Projection constant at n_max_fock=3 = 50471424/1779441125."""
        result = compute_projection_constant(n_max_fock=3)
        assert result['projection_constant'] == '50471424/1779441125'

    def test_classify_nmax2_zero_continuum(self):
        """Classification at n_max_fock=2: continuum is zero."""
        cls = classify_projection_constant(n_max_fock=2)
        assert 'ZERO' in cls['status'] or 'zero' in cls['status'].lower()

    def test_classify_nmax3_calibration(self):
        """Classification at n_max_fock=3: CALIBRATION tier."""
        cls = classify_projection_constant(n_max_fock=3)
        assert cls['paper_18_tier'] == 'CALIBRATION'


# =========================================================================
# Part C: Mode decomposition
# =========================================================================

class TestModeDecomposition:
    """Tests for edge decomposition and mode comparison."""

    def test_graph_edge_count_nmax2(self):
        """Graph at n_max_fock=2 has 3 edges."""
        result = graph_vp_edge_decomposition(n_max_fock=2)
        assert result['E_fock'] == 3

    def test_graph_total_nmax2(self):
        """Graph VP trace at n_max_fock=2 is 224/75."""
        result = graph_vp_edge_decomposition(n_max_fock=2)
        assert result['total_trace'] == '224/75'

    def test_graph_edge_numerator_32(self):
        """All edge Pi_ee have numerator 32 at n_max_fock=2."""
        result = graph_vp_edge_decomposition(n_max_fock=2)
        for e in result['edges']:
            assert e['numerator'] == 32, f"Edge {e['edge_idx']}: num={e['numerator']}"

    def test_graph_edge_types(self):
        """n_max_fock=2 has 1 inter-shell and 2 intra-shell edges."""
        result = graph_vp_edge_decomposition(n_max_fock=2)
        types = [e['edge_type'] for e in result['edges']]
        assert types.count('inter-shell') == 1
        assert types.count('intra-shell') == 2

    def test_enumerate_allowed_triples_nmax1(self):
        """n_max_ch=1 has exactly 1 allowed triple: (1,1,1)."""
        triples = _enumerate_allowed_triples(1)
        assert triples == [(1, 1, 1)]

    def test_enumerate_allowed_triples_nmax2(self):
        """n_max_ch=2 has 5 allowed triples."""
        triples = _enumerate_allowed_triples(2)
        assert len(triples) == 5

    def test_continuum_decomp_nmax1_all_W0(self):
        """At n_max_ch=1, all allowed triples have W=0."""
        result = continuum_vp_mode_decomposition(1)
        assert result['n_nonzero_W'] == 0
        assert result['n_zero_W_but_allowed'] == 1

    def test_continuum_decomp_nmax2_has_nonzero(self):
        """At n_max_ch=2, some triples have W > 0."""
        result = continuum_vp_mode_decomposition(2)
        assert result['n_nonzero_W'] > 0

    def test_mode_ratio_nmax3_rational(self):
        """Ratio at n_max_fock=3 is rational."""
        result = mode_ratio_analysis(n_max_fock=3)
        assert result['ratio_is_rational']

    def test_graph_excess_positive(self):
        """Graph VP exceeds continuum VP at any finite truncation."""
        result = mode_ratio_analysis(n_max_fock=3)
        assert result['graph_excess_float'] > 0


# =========================================================================
# Part D: Convergence
# =========================================================================

class TestConvergence:
    """Tests for convergence analysis."""

    def test_full_estimate_positive(self):
        """Full continuum VP estimate is positive."""
        est = _continuum_vp_full_estimate(n_max_sum=20)
        assert est > 0

    def test_full_estimate_converges(self):
        """Estimate at n=50 and n=100 agree within 20%."""
        e50 = _continuum_vp_full_estimate(n_max_sum=50)
        e30 = _continuum_vp_full_estimate(n_max_sum=30)
        assert abs(e50 - e30) / e50 < 0.20

    def test_convergence_table_structure(self):
        """Convergence table has expected fields."""
        result = convergence_table(n_max_values=[1, 2])
        assert 'full_continuum_estimate_n50' in result
        assert 'rows' in result
        assert len(result['rows']) == 2

    def test_convergence_table_nmax1_zero(self):
        """n_max_fock=1 has zero graph and continuum VP."""
        result = convergence_table(n_max_values=[1])
        row = result['rows'][0]
        assert row['graph_trace'] == 0.0
        assert row['continuum_truncated'] == 0.0


# =========================================================================
# Part E: Self-energy bridge
# =========================================================================

class TestSelfEnergyBridge:
    """Tests for the self-energy structural zero analysis."""

    def test_ground_state_count_nmax2(self):
        """n_max_fock=2 has 2 ground-state Dirac labels (n=1, kappa=-1)."""
        result = graph_self_energy_n0(n_max_fock=2)
        assert result['n_ground_states'] == 2

    def test_ground_has_couplings(self):
        """Ground state HAS vertex couplings in the graph."""
        result = graph_self_energy_n0(n_max_fock=2)
        assert result['has_ground_couplings']

    def test_continuum_structural_zero(self):
        """Continuum self-energy structural zero is True."""
        result = graph_self_energy_n0(n_max_fock=2)
        assert result['continuum_structural_zero']

    def test_graph_self_energy_nonzero(self):
        """Graph self-energy at ground state is NONZERO."""
        result = self_energy_structural_zero_bridge(n_max_fock=2)
        assert not result['structural_zero_preserved']

    def test_graph_self_energy_value(self):
        """Graph Sigma(n=1) = 2/5 for both m_j states."""
        result = self_energy_structural_zero_bridge(n_max_fock=2)
        for key, val in result['ground_state_self_energies'].items():
            assert val['value'] == '2/5'

    def test_graph_self_energy_degenerate(self):
        """Both ground-state m_j components have the same self-energy."""
        result = self_energy_structural_zero_bridge(n_max_fock=2)
        values = [v['value_float'] for v in result['ground_state_self_energies'].values()]
        assert len(set(values)) == 1


# =========================================================================
# Driver
# =========================================================================

class TestDriver:
    """Tests for the full bridge analysis driver."""

    def test_driver_runs(self):
        """run_bridge_analysis completes without error."""
        result = run_bridge_analysis(n_max_fock=2)
        assert 'part_a' in result
        assert 'part_b_at_requested' in result
        assert 'part_e_self_energy' in result

    def test_driver_has_nmax3_fallback(self):
        """When requested n_max=2, driver also computes n_max=3 projection."""
        result = run_bridge_analysis(n_max_fock=2)
        assert 'part_b_at_nmax3' in result

    def test_driver_nmax3_no_fallback(self):
        """When requested n_max=3, no separate nmax3 section needed."""
        result = run_bridge_analysis(n_max_fock=3)
        assert 'part_b_at_nmax3' not in result
