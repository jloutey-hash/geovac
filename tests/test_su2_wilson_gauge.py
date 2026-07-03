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

        Tolerance 1e-13: measured worst |dS| over 200 seeded Haar trials
        is 1.4e-15 (machine precision), so 1e-13 keeps ~50x margin while
        matching Paper 30's machine-precision claim (previously the assert
        was a loose 1e-10).
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

        assert np.isclose(S_before, S_after, atol=1e-13)

    def test_action_invariant_many_trials_machine_precision(self):
        """Gauge invariance at machine precision across 50 fresh Haar trials.

        Backs Paper 30 Sec. 'Gauge transformation': the demonstrated
        invariance residual is O(1e-15); assert < 1e-13 per trial.
        """
        edges = [(0, 1), (1, 2), (2, 0)]
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        for t in range(50):
            rng = np.random.default_rng(1000 + t)
            links = {edge: su2_random(rng) for edge in edges}
            S_b = wilson_action([plaq], links, 2.0)
            gauge = {i: su2_random(rng) for i in range(3)}
            S_a = wilson_action([plaq], gauge_transform(links, gauge), 2.0)
            assert abs(S_b - S_a) < 1e-13

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
    # Tolerance 1e-13: measured worst |dS| over 100 seeded Haar trials on
    # this graph is 2.7e-15 (machine precision); previously a loose 1e-8.
    gauge = {i: su2_random(rng) for i in range(V)}
    links_gauged = gauge_transform(links, gauge)
    S_after = wilson_action(plaqs, links_gauged, beta)
    assert np.isclose(S_before, S_after, atol=1e-13)


# ---------------------------------------------------------------------------
# 10. Paper 30 Result 2: weak-coupling kinetic term (1/8 coefficient and
#     the plaquette / curl quadratic form) -- mirrors the SU(3) sibling
#     tests/test_su3_wilson_s5.py::TestWeakCouplingKinetic (1/12).
# ---------------------------------------------------------------------------

def _pauli_matrices():
    return [
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]


def _hopf_graph_edge_data(max_n=3):
    """Forward edges, plaquettes, and incidence matrices at n_max=3.

    Returns (forward, eidx, plaqs, B, D) where B is the V x E signed
    node-edge incidence and D is the p x E signed plaquette-edge
    incidence (rows = signed indicator vectors of the Ihara 4-cycles).
    """
    lat = GeometricLattice(max_n=max_n)
    A = lat.adjacency.toarray()
    V = lat.num_states
    oriented, _ = enumerate_oriented_edges(A)
    forward = [(e.source, e.target) for e in oriented if e.source < e.target]
    eidx = {k: i for i, k in enumerate(forward)}
    plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
    B = np.zeros((V, len(forward)))
    for (u, v), i in eidx.items():
        B[u, i] = 1.0
        B[v, i] = -1.0
    D = np.zeros((len(plaqs), len(forward)))
    for pi, P in enumerate(plaqs):
        for e in P:
            if (e.source, e.target) in eidx:
                D[pi, eidx[(e.source, e.target)]] += 1.0
            else:
                D[pi, eidx[(e.target, e.source)]] -= 1.0
    return forward, eidx, plaqs, B, D


class TestWeakCouplingKinetic:
    """Paper 30 Prop. 3 (Result 2): S_W = (beta/8) sum_a <A^a, M A^a> + O(A^4)
    with M = B2 B2^T the plaquette (curl) quadratic form."""

    def test_kinetic_general_formula(self):
        """SU(N_c) coefficient 1/(4 N_c): SU(2) -> 1/8, SU(3) -> 1/12."""
        import sympy as sp
        assert sp.Rational(1, 4 * 2) == sp.Rational(1, 8)
        assert sp.Rational(1, 4 * 3) == sp.Rational(1, 12)
        # ratio SU(2)/SU(3) = 3/2
        assert sp.Rational(1, 8) / sp.Rational(1, 12) == sp.Rational(3, 2)

    def test_quadratic_expansion_numerically(self):
        """
        Verify numerically that for small A, with U_e = exp(i A_e^a sigma_a/2),

            1 - (1/2) Re tr U_P = (1/8) sum_a (sum_e A_e^a)^2 + O(A^4).

        Mirrors the SU(3) sibling (1/12). Measured |diff| <= 2e-7 over
        5 seeded trials at eps = 0.01; assert < 1e-5.
        """
        from scipy.linalg import expm
        sigma = _pauli_matrices()
        rng = np.random.default_rng(42)
        eps = 0.01
        for _ in range(5):
            A_list = [eps * rng.standard_normal(3) for _ in range(3)]
            U_P = np.eye(2, dtype=complex)
            for A in A_list:
                X = sum(A[a] * sigma[a] / 2.0 for a in range(3))
                U_P = U_P @ expm(1j * X)
            A_P = sum(A_list)
            S_predicted = (1.0 / 8.0) * float(np.sum(A_P ** 2))
            S_actual = 1.0 - 0.5 * float(np.real(np.trace(U_P)))
            assert abs(S_actual - S_predicted) < 1e-5

    def test_quadratic_expansion_on_hopf_graph(self):
        """Full Wilson action on the n_max=3 Coulomb graph vs the quadratic
        form (beta/8) sum_a ||D A^a||^2 built from the plaquette incidence.

        Measured relative deviation (seed 7): 2.0e-3 at eps=0.01 and
        7.3e-5 at eps=0.005 -- shrinking with eps as the O(A^4) remainder
        predicts. Assert both, with margin.
        """
        from scipy.linalg import expm
        sigma = _pauli_matrices()
        forward, eidx, plaqs, B, D = _hopf_graph_edge_data(3)
        E = len(forward)
        beta = 1.7
        rng = np.random.default_rng(7)
        rel_diffs = []
        for eps in [0.01, 0.005]:
            Afield = eps * rng.standard_normal((E, 3))
            links = {}
            for (u, v), i in eidx.items():
                X = sum(Afield[i, a] * sigma[a] / 2.0 for a in range(3))
                links[(u, v)] = expm(1j * X)
            S_act = wilson_action(plaqs, links, beta)
            S_quad = (beta / 8.0) * sum(
                float((D @ Afield[:, a]) @ (D @ Afield[:, a]))
                for a in range(3))
            rel_diffs.append(abs(S_act - S_quad) / S_quad)
        assert rel_diffs[0] < 5e-3
        assert rel_diffs[1] < 5e-4
        # O(A^4) remainder: halving eps must shrink the relative deviation
        assert rel_diffs[1] < rel_diffs[0]

    def test_kinetic_form_is_curl_complement_of_L1(self):
        """The kinetic quadratic form M = B2 B2^T is NOT a scalar multiple
        of Paper 25's edge Laplacian L1 = B^T B; it is its curl (up-part)
        complement: rank M = beta_1 = 2 while rank L1 = E - beta_1 = 11,
        and L1 @ B2 = 0 exactly (the plaquette rows span ker L1).
        Together they assemble the full Hodge-1 Laplacian of the 2-complex
        obtained by attaching the Ihara plaquettes as faces.

        This pins the corrected statement of Paper 30 Prop. 3 (the earlier
        'coincides with L1 up to a positive scalar' wording was wrong:
        the two forms have complementary supports).
        """
        forward, eidx, plaqs, B, D = _hopf_graph_edge_data(3)
        E = len(forward)
        assert E == 13
        assert len(plaqs) == 2
        L1 = B.T @ B
        M = D.T @ D
        # Integer matrices: exact comparisons are safe.
        assert np.linalg.matrix_rank(M) == 2          # = beta_1
        assert np.linalg.matrix_rank(L1) == 11        # = E - beta_1
        # Curl rows live exactly in ker L1:
        assert np.allclose(L1 @ D.T, 0.0, atol=0.0)
        # Not proportional (different ranks already prove it; make explicit):
        nz = np.abs(L1) > 0.5
        ratios = M[nz] / L1[nz]
        assert ratios.max() - ratios.min() > 0.5
        # Full Hodge-1 Laplacian of the 2-complex is invertible on the
        # orthocomplement of harmonic space; here beta_1 cycles are all
        # filled by plaquettes, so L1 + M has full rank E.
        assert np.linalg.matrix_rank(L1 + M) == E


# ---------------------------------------------------------------------------
# 11. Paper 30 Result 3: Monte Carlo Wilson-loop table (seeded, slow)
# ---------------------------------------------------------------------------

class TestMonteCarloWilsonLoop:
    """Pins Table 'Wilson-loop expectation' of Paper 30 (Result 3):
    <W>_MC = 0.126 / 0.158 / 0.332 / 0.681 at beta = 0.5 / 1 / 2 / 5,
    n_max = 3, loop = first p-block 4-cycle, 2000 samples after 500
    thermalization steps, seed = 12.

    The Metropolis chain is fully deterministic for a fixed seed
    (verified: two consecutive runs agree bit-exactly), so the pins use
    atol = 1e-6 (guard against libm variation across platforms; the
    paper's numbers are quoted to 3 decimals).
    """

    def _setup(self):
        from geovac.su2_wilson_gauge import monte_carlo_wilson_expectation
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        oriented, _ = enumerate_oriented_edges(A)
        forward = [(e.source, e.target)
                   for e in oriented if e.source < e.target]
        edge_list = [OrientedEdge(u, v) for (u, v) in forward]
        plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
        assert len(forward) == 13   # "13 forward links" (paper)
        assert len(plaqs) == 2
        # The loop is the first primitive 4-cycle in the p-block:
        # (2,1,-1) -> (2,1,0) -> (3,1,0) -> (3,1,-1) -> close.
        loop_states = [lat.states[e.source] for e in plaqs[0]]
        assert loop_states == [(2, 1, -1), (2, 1, 0), (3, 1, 0), (3, 1, -1)]
        return monte_carlo_wilson_expectation, plaqs, edge_list

    @pytest.mark.slow
    def test_mc_wilson_table_seed12(self):
        mc, plaqs, edge_list = self._setup()
        # (beta, mean, stderr) recomputed live 2026-07-02; means round to
        # the paper table 0.126 / 0.158 / 0.332 / 0.681.
        expected = {
            0.5: (0.12576712406708293, 0.011194572527873982),
            1.0: (0.1575325588927505, 0.010766717015609136),
            2.0: (0.3324166614690038, 0.009133638708602723),
            5.0: (0.6812474363347273, 0.00603985232651571),
        }
        paper_table = {0.5: 0.126, 1.0: 0.158, 2.0: 0.332, 5.0: 0.681}
        means = {}
        for beta, (m_exp, s_exp) in expected.items():
            m, s = mc(plaqs[0], plaqs, edge_list, beta,
                      n_samples=2000, n_thermalize=500, seed=12)
            means[beta] = m
            assert abs(m - m_exp) < 1e-6, (beta, m, m_exp)
            assert abs(s - s_exp) < 1e-6, (beta, s, s_exp)
            assert abs(round(m, 3) - paper_table[beta]) < 1e-12
        # Monotonic in beta (strong -> weak coupling)
        bs = sorted(means)
        assert all(means[bs[i]] < means[bs[i + 1]] for i in range(len(bs) - 1))

    @pytest.mark.slow
    def test_mc_seed_determinism(self):
        """Two runs with the same seed agree bit-exactly."""
        mc, plaqs, edge_list = self._setup()
        m1, s1 = mc(plaqs[0], plaqs, edge_list, 0.5,
                    n_samples=2000, n_thermalize=500, seed=12)
        m2, s2 = mc(plaqs[0], plaqs, edge_list, 0.5,
                    n_samples=2000, n_thermalize=500, seed=12)
        assert m1 == m2
        assert s1 == s2
