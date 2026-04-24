"""Consistency tests between QED spectral modules and the authoritative
Dirac-on-S³ infrastructure (dirac_s3.py).

The QED modules (qed_vertex, qed_self_energy, qed_three_loop) define local
spectrum helpers for mpmath precision.  These must agree with the authoritative
sympy-Rational definitions in dirac_s3.py.  This test file is the bridge
that catches convention drift between the two tracks.
"""

import pytest
import mpmath

from geovac.dirac_s3 import (
    dirac_eigenvalue_abs,
    dirac_degeneracy,
    fock_to_ch,
    ch_to_fock,
)
from geovac.ihara_zeta_dirac import build_dirac_s3_graph, DiracLabel
from geovac.dirac_matrix_elements import kappa_to_l


# ---------------------------------------------------------------------------
# Collect all QED local spectrum helpers (private, but we test them here)
# ---------------------------------------------------------------------------

from geovac.qed_vertex import (
    _lambda_n as vertex_lambda,
    _g_n_dirac as vertex_g,
    _mu_q as vertex_mu,
    _d_q_transverse as vertex_dq,
    _vertex_allowed as vertex_allowed,
)

from geovac.qed_self_energy import (
    _lambda_n as se_lambda,
    _g_n_dirac as se_g,
    _mu_q as se_mu,
    _d_q_transverse as se_dq,
    _vertex_allowed as se_vertex_allowed,
)

from geovac.qed_three_loop import (
    _lambda_n as tl_lambda,
    _g_n_dirac as tl_g,
    _mu_q as tl_mu,
    _d_q_transverse as tl_dq,
    _vertex_allowed as tl_vertex_allowed,
)


N_MAX_TEST = 10


# ---------------------------------------------------------------------------
# 1. Dirac eigenvalue consistency: QED locals vs dirac_s3 authority
# ---------------------------------------------------------------------------

class TestDiracEigenvalueConsistency:
    """All QED _lambda_n must agree with dirac_s3.dirac_eigenvalue_abs."""

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_vertex_lambda(self, n):
        expected = float(dirac_eigenvalue_abs(n, convention="ch"))
        assert float(vertex_lambda(n)) == pytest.approx(expected, abs=1e-30)

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_self_energy_lambda(self, n):
        expected = float(dirac_eigenvalue_abs(n, convention="ch"))
        assert float(se_lambda(n)) == pytest.approx(expected, abs=1e-30)

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_three_loop_lambda(self, n):
        expected = float(dirac_eigenvalue_abs(n, convention="ch"))
        assert float(tl_lambda(n)) == pytest.approx(expected, abs=1e-30)


# ---------------------------------------------------------------------------
# 2. Dirac degeneracy consistency
# ---------------------------------------------------------------------------

class TestDiracDegeneracyConsistency:
    """All QED _g_n_dirac must agree with dirac_s3.dirac_degeneracy."""

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_vertex_degeneracy(self, n):
        expected = dirac_degeneracy(n, sector="dirac", convention="ch")
        assert int(vertex_g(n)) == expected

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_self_energy_degeneracy(self, n):
        expected = dirac_degeneracy(n, sector="dirac", convention="ch")
        assert int(se_g(n)) == expected

    @pytest.mark.parametrize("n", range(N_MAX_TEST + 1))
    def test_three_loop_degeneracy(self, n):
        expected = dirac_degeneracy(n, sector="dirac", convention="ch")
        assert int(tl_g(n)) == expected


# ---------------------------------------------------------------------------
# 3. Hodge-1 / photon spectrum consistency across QED modules
# ---------------------------------------------------------------------------

class TestHodgeSpectrumConsistency:
    """mu_q and d_q^T must agree across all QED modules."""

    @pytest.mark.parametrize("q", range(1, N_MAX_TEST + 1))
    def test_mu_q_cross_module(self, q):
        vals = [float(vertex_mu(q)), float(se_mu(q)), float(tl_mu(q))]
        assert vals[0] == pytest.approx(vals[1], abs=1e-30)
        assert vals[0] == pytest.approx(vals[2], abs=1e-30)
        assert vals[0] == pytest.approx(q * (q + 2), abs=1e-30)

    @pytest.mark.parametrize("q", range(1, N_MAX_TEST + 1))
    def test_dq_transverse_cross_module(self, q):
        vals = [float(vertex_dq(q)), float(se_dq(q)), float(tl_dq(q))]
        assert vals[0] == pytest.approx(vals[1], abs=1e-30)
        assert vals[0] == pytest.approx(vals[2], abs=1e-30)
        assert vals[0] == pytest.approx(q * (q + 2), abs=1e-30)


# ---------------------------------------------------------------------------
# 4. Vertex selection rule consistency across QED modules
# ---------------------------------------------------------------------------

class TestVertexSelectionConsistency:
    """_vertex_allowed must agree across all QED modules."""

    @pytest.mark.parametrize("n1", range(6))
    @pytest.mark.parametrize("n2", range(6))
    @pytest.mark.parametrize("ng", range(8))
    def test_selection_rule_cross_module(self, n1, n2, ng):
        v = vertex_allowed(n1, n2, ng)
        s = se_vertex_allowed(n1, n2, ng)
        t = tl_vertex_allowed(n1, n2, ng)
        assert v == s, f"vertex vs self_energy disagree at ({n1},{n2},{ng})"
        assert v == t, f"vertex vs three_loop disagree at ({n1},{n2},{ng})"


# ---------------------------------------------------------------------------
# 5. Vertex selection rule vs Rule B adjacency (structural bridge)
# ---------------------------------------------------------------------------

class TestVertexVsRuleBStructure:
    """The QED vertex selection rule and the Dirac graph Rule B encode
    related but DISTINCT selection rules:

    - Vertex: SO(4) gamma-matrix parity on level indices (n1+n2+n_gamma odd)
    - Rule B: E1 dipole selection on individual node labels (Δl=±1)

    Rule B is MORE PERMISSIVE than the vertex rule.  Example: the
    n_CH=0→1 transition with n_gamma=1 gives parity sum 0+1+1=2 (even),
    forbidden by the vertex rule, but Rule B allows individual 1s→2p
    node edges within those levels.

    The correct structural relationship is: vertex-allowed level
    transitions imply Rule-B-connected node pairs (not the reverse).
    """

    def test_vertex_implies_rule_b_connected(self):
        """If vertex_allowed(n1, n2, ng), then there exist nodes at
        Dirac levels n1 and n2 that are connected by Rule B."""
        n_max = 3
        adj, labels, _, _ = build_dirac_s3_graph(n_max, adjacency_rule='B')

        # Build level-to-level connectivity from Rule B adjacency
        level_connected = set()
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if adj[i, j] > 0:
                    n1_ch = fock_to_ch(labels[i].n_fock)
                    n2_ch = fock_to_ch(labels[j].n_fock)
                    level_connected.add((min(n1_ch, n2_ch), max(n1_ch, n2_ch)))

        # Every vertex-allowed (n1, n2) pair should appear in level_connected
        n_max_ch = fock_to_ch(n_max)
        for n1 in range(n_max_ch + 1):
            for n2 in range(n1, n_max_ch + 1):
                for ng in range(1, n1 + n2 + 2):
                    if vertex_allowed(n1, n2, ng):
                        pair = (min(n1, n2), max(n1, n2))
                        assert pair in level_connected, (
                            f"vertex_allowed({n1},{n2},{ng}) but no Rule B "
                            f"edges between levels {n1} and {n2}"
                        )

    def test_rule_b_strictly_more_permissive(self):
        """Rule B allows level pairs that the vertex rule forbids.
        Specifically, n_CH=0→1 (1s→2p) is Rule-B-allowed but
        vertex-forbidden (0+1+1=2 is even)."""
        n_max = 3
        adj, labels, _, _ = build_dirac_s3_graph(n_max, adjacency_rule='B')

        # Find Rule-B-connected level pairs
        rule_b_levels = set()
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if adj[i, j] > 0:
                    n1_ch = fock_to_ch(labels[i].n_fock)
                    n2_ch = fock_to_ch(labels[j].n_fock)
                    rule_b_levels.add((min(n1_ch, n2_ch), max(n1_ch, n2_ch)))

        # Find vertex-allowed level pairs
        n_max_ch = fock_to_ch(n_max)
        vertex_levels = set()
        for n1 in range(n_max_ch + 1):
            for n2 in range(n1, n_max_ch + 1):
                for ng in range(1, n1 + n2 + 2):
                    if vertex_allowed(n1, n2, ng):
                        vertex_levels.add((n1, n2))

        # Rule B should be strictly more permissive
        assert vertex_levels <= rule_b_levels
        extra = rule_b_levels - vertex_levels
        assert len(extra) > 0, "Expected Rule B to be strictly more permissive"
        assert (0, 1) in extra, "n_CH=0→1 should be in Rule B but not vertex"

    def test_vertex_parity_when_allowed(self):
        """When vertex_allowed returns True, the parity sum is odd."""
        for n1 in range(6):
            for n2 in range(6):
                for ng in range(8):
                    if vertex_allowed(n1, n2, ng):
                        assert (n1 + n2 + ng) % 2 == 1


# ---------------------------------------------------------------------------
# 6. Convention sanity: Fock <-> CH round-trip
# ---------------------------------------------------------------------------

class TestConventionSanity:
    """Fock/CH conversion round-trips and QED modules use CH convention."""

    @pytest.mark.parametrize("n_fock", range(1, N_MAX_TEST + 1))
    def test_fock_ch_roundtrip(self, n_fock):
        assert ch_to_fock(fock_to_ch(n_fock)) == n_fock

    def test_qed_convention_is_ch(self):
        """QED n=0 must give λ=3/2 (CH convention, not Fock n=1 -> λ=5/2)."""
        assert float(vertex_lambda(0)) == pytest.approx(1.5, abs=1e-30)
        assert float(se_lambda(0)) == pytest.approx(1.5, abs=1e-30)
        assert float(tl_lambda(0)) == pytest.approx(1.5, abs=1e-30)

    def test_qed_degeneracy_at_zero(self):
        """n_CH=0 has g=4 (2*(0+1)*(0+2)=4), not g=12 (Fock n=1)."""
        assert int(vertex_g(0)) == 4
        assert int(se_g(0)) == 4
        assert int(tl_g(0)) == 4
