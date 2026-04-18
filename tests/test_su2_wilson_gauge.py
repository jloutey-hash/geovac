"""
Tests for geovac/su2_wilson_gauge.py.

Verifies:
  1. Link variables are in SU(2) (det = 1, unitary).
  2. Plaquette holonomy is in SU(2).
  3. Wilson action is gauge-invariant under U_e -> g_{src} U_e g_{tgt}^dagger.
  4. U(1) reduction: diagonal SU(2) link variables give the U(1) Wilson
     action of Paper 25.
  5. Character expansion coefficients converge and reduce to 1 at
     beta = 0 (d_0 = 1), all d_j vanish for j >= 1.
  6. Plaquette enumeration agrees with Ihara zeta primitive-walk count.
"""

import numpy as np
import pytest

from geovac.lattice import GeometricLattice
from geovac.su2_wilson_gauge import (
    OrientedEdge,
    diagonal_su2_from_phase,
    edge_link_variable,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    expectation_wilson_loop,
    gauge_transform,
    is_su2,
    partition_function_character_expansion,
    plaquette_holonomy,
    su2_character,
    su2_character_coefficient,
    su2_from_angle_axis,
    su2_from_quaternion,
    su2_random,
    u1_action_from_su2,
    wilson_action,
)


# ---------------------------------------------------------------------------
# 1. Link variables are in SU(2)
# ---------------------------------------------------------------------------

class TestSU2Membership:
    def test_identity_is_su2(self):
        U = su2_from_quaternion(1.0, 0.0, 0.0, 0.0)
        assert is_su2(U)

    def test_pauli_x_rotation_is_su2(self):
        # pi rotation about x-axis: U = i sigma_x
        U = su2_from_angle_axis(np.pi, [1, 0, 0])
        assert is_su2(U)

    def test_haar_random_is_su2(self):
        rng = np.random.default_rng(42)
        for _ in range(20):
            U = su2_random(rng)
            assert is_su2(U, atol=1e-10)

    def test_quaternion_nonunit_raises(self):
        with pytest.raises(ValueError):
            su2_from_quaternion(1.0, 1.0, 0.0, 0.0)  # norm = sqrt(2)


# ---------------------------------------------------------------------------
# 2. Plaquette holonomy is in SU(2)
# ---------------------------------------------------------------------------

class TestPlaquetteHolonomy:
    def test_trivial_plaquette_is_identity(self):
        # Three identity links around a triangle -> identity holonomy
        links = {
            (0, 1): np.eye(2, dtype=complex),
            (1, 2): np.eye(2, dtype=complex),
            (2, 0): np.eye(2, dtype=complex),
        }
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        U_P = plaquette_holonomy(plaq, links)
        assert np.allclose(U_P, np.eye(2), atol=1e-12)
        assert is_su2(U_P)

    def test_random_plaquette_is_su2(self):
        rng = np.random.default_rng(7)
        links = {
            (0, 1): su2_random(rng),
            (1, 2): su2_random(rng),
            (2, 0): su2_random(rng),
        }
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        U_P = plaquette_holonomy(plaq, links)
        assert is_su2(U_P, atol=1e-10)

    def test_reverse_edge_gives_conjugate(self):
        rng = np.random.default_rng(9)
        U = su2_random(rng)
        links = {(0, 1): U}
        # Forward
        e_fwd = OrientedEdge(0, 1)
        assert np.allclose(edge_link_variable(e_fwd, links), U)
        # Reverse should give U^dagger
        e_rev = OrientedEdge(1, 0)
        assert np.allclose(edge_link_variable(e_rev, links), U.conj().T)


# ---------------------------------------------------------------------------
# 3. Wilson action is gauge-invariant
# ---------------------------------------------------------------------------

class TestGaugeInvariance:
    def test_action_invariant_under_gauge_transform(self):
        """
        Apply a random node-local gauge transformation and verify that
        the Wilson action is invariant.
        """
        rng = np.random.default_rng(123)
        # Triangle graph 0-1-2-0
        edges = [(0, 1), (1, 2), (2, 0)]
        links = {edge: su2_random(rng) for edge in edges}
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        plaquettes = [plaq]
        beta = 2.0
        S_before = wilson_action(plaquettes, links, beta)

        # Random gauge transformation
        gauge = {i: su2_random(rng) for i in range(3)}
        links_gauged = gauge_transform(links, gauge)
        S_after = wilson_action(plaquettes, links_gauged, beta)

        assert np.isclose(S_before, S_after, atol=1e-10)

    def test_identity_gauge_is_noop(self):
        rng = np.random.default_rng(1)
        links = {(0, 1): su2_random(rng), (1, 2): su2_random(rng)}
        gauge = {0: np.eye(2, dtype=complex), 1: np.eye(2, dtype=complex),
                 2: np.eye(2, dtype=complex)}
        gauged = gauge_transform(links, gauge)
        for k in links:
            assert np.allclose(links[k], gauged[k], atol=1e-14)


# ---------------------------------------------------------------------------
# 4. U(1) reduction of SU(2) action
# ---------------------------------------------------------------------------

class TestU1Reduction:
    def test_diagonal_su2_is_su2(self):
        phi = 1.234
        U = diagonal_su2_from_phase(phi)
        assert is_su2(U)

    def test_diagonal_character_is_cos(self):
        # (1/2) Re tr diag(e^{i phi}, e^{-i phi}) = cos(phi)
        for phi in [0.0, 0.3, 1.5, np.pi / 2, np.pi]:
            U = diagonal_su2_from_phase(phi)
            assert np.isclose(su2_character(U), np.cos(phi), atol=1e-12)

    def test_su2_action_reduces_to_u1_at_diagonal(self):
        """
        Set link variables to diagonal SU(2) matrices and compare the
        SU(2) Wilson action to the U(1) action in terms of phases.
        This is the maximal-torus reduction.
        """
        # Triangle plaquette
        phi = {(0, 1): 0.3, (1, 2): 0.7, (2, 0): -0.5}
        links = {k: diagonal_su2_from_phase(v) for k, v in phi.items()}
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        beta = 1.7
        S_su2 = wilson_action([plaq], links, beta)
        S_u1 = u1_action_from_su2([plaq], phi, beta)
        assert np.isclose(S_su2, S_u1, atol=1e-12)


# ---------------------------------------------------------------------------
# 5. Character expansion coefficients
# ---------------------------------------------------------------------------

class TestCharacterExpansion:
    def test_beta_zero_limit(self):
        # d_0(0) = 1 (trivial rep survives), d_j>=1 = 0
        assert np.isclose(su2_character_coefficient(0, 0.0), 1.0)
        assert np.isclose(su2_character_coefficient(1, 0.0), 0.0)
        assert np.isclose(su2_character_coefficient(2, 0.0), 0.0)

    def test_coefficients_decrease_with_j(self):
        # At fixed beta, d_j decreases with j (Bessel tail)
        beta = 1.0
        c_vals = [su2_character_coefficient(j, beta) for j in range(6)]
        for j in range(len(c_vals) - 1):
            assert c_vals[j] > c_vals[j + 1]

    def test_coefficients_positive(self):
        # Character coefficients I_j(beta) > 0 for beta > 0
        beta = 2.5
        for j in range(4):
            assert su2_character_coefficient(j, beta) > 0


# ---------------------------------------------------------------------------
# 6. Plaquette enumeration on S^3 graph
# ---------------------------------------------------------------------------

class TestPlaquetteEnumeration:
    def test_s3_nmax2_has_no_plaquettes(self):
        """At n_max=2 the S^3 Coulomb graph is a forest (beta_1 = 0).
        No closed non-backtracking walks."""
        lat = GeometricLattice(max_n=2)
        A = lat.adjacency.toarray()
        plaqs = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        assert len(plaqs) == 0

    def test_s3_nmax3_has_two_primitive_4cycles(self):
        """At n_max=3, beta_1 = 2, both cycles are 4-cycles in the p-block.
        With both_orientations=False we get 2 unoriented 4-cycles."""
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        # Only 4-cycles
        plaqs_4 = enumerate_plaquettes(A, max_length=4, both_orientations=False)
        length_4_count = sum(1 for p in plaqs_4 if len(p) == 4)
        assert length_4_count == 2

    def test_plaquette_is_closed(self):
        """Every enumerated plaquette must be a closed walk."""
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        for P in plaqs:
            assert P[0].source == P[-1].target, (
                f"Plaquette {P} is not closed"
            )

    def test_plaquette_has_no_backtracking(self):
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        for P in plaqs:
            for i in range(len(P) - 1):
                # Consecutive edges must not be reverses
                assert not (P[i].target == P[i + 1].source
                            and P[i].source == P[i + 1].target)


# ---------------------------------------------------------------------------
# 7. Partition function sanity
# ---------------------------------------------------------------------------

class TestPartitionFunction:
    def test_no_plaquette_Z_is_one(self):
        Z = partition_function_character_expansion([], 1.0, R_max=3)
        assert np.isclose(Z, 1.0)

    def test_Z_is_positive(self):
        # At any finite beta, with positive c_j, Z > 0
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        plaqs = enumerate_plaquettes(A, max_length=4)
        for beta in [0.1, 1.0, 5.0, 20.0]:
            Z = partition_function_character_expansion(plaqs, beta)
            assert Z > 0


# ---------------------------------------------------------------------------
# 8. Wilson loop expectation value
# ---------------------------------------------------------------------------

class TestWilsonLoop:
    def test_unknown_loop_gives_zero(self):
        # A loop not matching any plaquette returns 0 at leading order
        empty_loop = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        plaqs = []
        val = expectation_wilson_loop(empty_loop, plaqs, beta=1.0)
        assert val == 0.0

    def test_plaquette_loop_gives_bessel_ratio(self):
        # A loop that IS one of the plaquettes returns I_2/I_1
        from scipy.special import iv
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        plaquettes = [plaq]
        beta = 2.0
        val = expectation_wilson_loop(plaq, plaquettes, beta)
        expected = iv(2, beta) / iv(1, beta)
        assert np.isclose(val, expected, atol=1e-10)

    def test_small_beta_wilson_loop_small(self):
        # At beta -> 0, I_2/I_1 -> beta/4 -> 0; strong coupling limit.
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        val = expectation_wilson_loop(plaq, [plaq], 0.01)
        assert abs(val) < 0.01

    def test_large_beta_wilson_loop_to_one(self):
        # At beta -> infty, I_2/I_1 -> 1; weak coupling limit, deconfined.
        # The asymptotic expansion I_2/I_1 ~ 1 - 3/(2 beta) + O(beta^-2).
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        val = expectation_wilson_loop(plaq, [plaq], 1000.0)
        assert np.isclose(val, 1.0, atol=0.005)
        # And the value monotonically increases toward 1
        val10 = expectation_wilson_loop(plaq, [plaq], 10.0)
        val100 = expectation_wilson_loop(plaq, [plaq], 100.0)
        assert val10 < val100 < val < 1.0


# ---------------------------------------------------------------------------
# 9. Full smoke test on S^3 graph at max_n = 3
# ---------------------------------------------------------------------------

def test_s3_coulomb_max_n3_endtoend():
    """
    End-to-end test: build the S^3 Coulomb graph at max_n=3, enumerate
    plaquettes, assign random SU(2) link variables, compute the Wilson
    action, check gauge invariance.
    """
    lat = GeometricLattice(max_n=3)
    A = lat.adjacency.toarray()
    V = lat.num_states

    oriented, _ = enumerate_oriented_edges(A)
    plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
    assert len(plaqs) > 0

    rng = np.random.default_rng(11)
    # Assign random link variables (forward direction only)
    links = {}
    for e in oriented:
        # Only set forward direction for each undirected edge
        if e.source < e.target:
            links[(e.source, e.target)] = su2_random(rng)

    beta = 1.5
    S_before = wilson_action(plaqs, links, beta)

    # Random gauge transformation
    gauge = {i: su2_random(rng) for i in range(V)}
    links_gauged = gauge_transform(links, gauge)
    S_after = wilson_action(plaqs, links_gauged, beta)
    assert np.isclose(S_before, S_after, atol=1e-8)
