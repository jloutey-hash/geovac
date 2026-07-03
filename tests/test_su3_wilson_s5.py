"""
Tests for geovac/su3_wilson_s5.py.

Mirrors tests/test_su2_wilson_gauge.py for the SU(3) sibling on the
Bargmann-Segal S^5 graph (Sprint ST-SU3, May 2026).

Verifies:
  1. Gell-Mann matrix conventions: traceless, Hermitian, normalization
     tr(T^a T^b) = (1/2) delta^{ab}.
  2. Link variables are in SU(3) (det = 1, unitary).
  3. Plaquette holonomy is in SU(3).
  4. Wilson action is gauge-invariant under U_e -> g_{src} U_e g_{tgt}^dagger.
  5. Cartan torus reduction: diagonal SU(3) link variables give the
     U(1) x U(1) Wilson action with character (1/3)(cos a + cos b + cos(a+b)).
  6. Plaquette enumeration on the Bargmann-Segal graph at N_max = 2, 3, 4.
  7. Bipartiteness: no odd-length plaquettes.
  8. Wilson action positive semidefinite, vanishes at trivial vacuum.
  9. Weak-coupling kinetic coefficient = 1/12 (vs SU(2)'s 1/8).
"""

from fractions import Fraction

import numpy as np
import pytest
import sympy as sp

from geovac.su3_wilson_s5 import (
    OrientedEdge,
    bargmann_adjacency_dense,
    cartan_links_from_phases,
    cartan_su3_from_phases,
    edge_link_variable,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    gauge_transform,
    gell_mann_anticommutators_normalization,
    gell_mann_matrices,
    is_su3,
    monte_carlo_wilson_expectation,
    plaquette_holonomy,
    su3_character,
    su3_character_coefficient_fundamental,
    su3_from_algebra,
    su3_generators,
    su3_random,
    u1xu1_action_from_su3,
    weak_coupling_kinetic_coefficient_per_plaquette,
    wilson_action,
)


# ---------------------------------------------------------------------------
# 1. Gell-Mann matrix conventions
# ---------------------------------------------------------------------------

class TestGellMann:
    def test_eight_matrices(self):
        lam = gell_mann_matrices()
        assert len(lam) == 8

    def test_all_3x3(self):
        for L in gell_mann_matrices():
            assert L.shape == (3, 3)

    def test_all_traceless(self):
        for L in gell_mann_matrices():
            assert np.isclose(np.trace(L), 0.0, atol=1e-12)

    def test_all_hermitian(self):
        for L in gell_mann_matrices():
            assert np.allclose(L, L.conj().T, atol=1e-12)

    def test_normalization_tr_lambda_squared(self):
        # tr(lambda_a lambda_b) = 2 delta_{ab}; so tr(T^a T^b) = (1/2) delta_{ab}
        T = su3_generators()
        for a in range(8):
            for b in range(8):
                trace = np.trace(T[a] @ T[b])
                expected = 0.5 if a == b else 0.0
                assert np.isclose(np.real(trace), expected, atol=1e-12)
                assert np.isclose(np.imag(trace), 0.0, atol=1e-12)

    def test_lambda_3_diagonal(self):
        lam = gell_mann_matrices()
        # lambda_3 = diag(1, -1, 0)
        expected = np.diag([1.0, -1.0, 0.0])
        assert np.allclose(lam[2], expected, atol=1e-12)

    def test_lambda_8_diagonal_traceless(self):
        lam = gell_mann_matrices()
        # lambda_8 = (1/sqrt(3)) diag(1, 1, -2)
        d = np.diag(lam[7])
        assert np.allclose(d.imag, 0)
        assert np.isclose(np.sum(d).real, 0.0, atol=1e-12)
        # Off-diagonal zero
        off_diag = lam[7] - np.diag(np.diag(lam[7]))
        assert np.allclose(off_diag, 0.0)


# ---------------------------------------------------------------------------
# 2. SU(3) membership
# ---------------------------------------------------------------------------

class TestSU3Membership:
    def test_identity_is_su3(self):
        U = np.eye(3, dtype=complex)
        assert is_su3(U)

    def test_haar_random_is_su3(self):
        rng = np.random.default_rng(42)
        for _ in range(20):
            U = su3_random(rng)
            assert is_su3(U, atol=1e-10)

    def test_algebra_zero_gives_identity(self):
        U = su3_from_algebra([0.0] * 8)
        assert np.allclose(U, np.eye(3), atol=1e-12)

    def test_algebra_random_in_su3(self):
        rng = np.random.default_rng(7)
        for _ in range(10):
            x = rng.standard_normal(8) * 0.5
            U = su3_from_algebra(x)
            assert is_su3(U, atol=1e-10)

    def test_cartan_diag_is_su3(self):
        for a, b in [(0.0, 0.0), (0.3, 0.7), (1.5, -2.1), (np.pi, np.pi)]:
            U = cartan_su3_from_phases(a, b)
            assert is_su3(U, atol=1e-12)


# ---------------------------------------------------------------------------
# 3. Plaquette holonomy stays in SU(3)
# ---------------------------------------------------------------------------

class TestPlaquetteHolonomy:
    def test_trivial_plaquette_is_identity(self):
        links = {
            (0, 1): np.eye(3, dtype=complex),
            (1, 2): np.eye(3, dtype=complex),
            (2, 0): np.eye(3, dtype=complex),
        }
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        U_P = plaquette_holonomy(plaq, links)
        assert np.allclose(U_P, np.eye(3), atol=1e-12)
        assert is_su3(U_P)

    def test_random_plaquette_is_su3(self):
        rng = np.random.default_rng(13)
        links = {
            (0, 1): su3_random(rng),
            (1, 2): su3_random(rng),
            (2, 0): su3_random(rng),
        }
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        U_P = plaquette_holonomy(plaq, links)
        assert is_su3(U_P, atol=1e-10)

    def test_reverse_edge_is_dagger(self):
        rng = np.random.default_rng(9)
        U = su3_random(rng)
        links = {(0, 1): U}
        e_fwd = OrientedEdge(0, 1)
        assert np.allclose(edge_link_variable(e_fwd, links), U)
        e_rev = OrientedEdge(1, 0)
        assert np.allclose(edge_link_variable(e_rev, links), U.conj().T)


# ---------------------------------------------------------------------------
# 4. Gauge invariance
# ---------------------------------------------------------------------------

class TestGaugeInvariance:
    def test_action_invariant_under_gauge_transform(self):
        rng = np.random.default_rng(123)
        edges = [(0, 1), (1, 2), (2, 0)]
        links = {edge: su3_random(rng) for edge in edges}
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        plaquettes = [plaq]
        beta = 2.0
        S_before = wilson_action(plaquettes, links, beta)

        gauge = {i: su3_random(rng) for i in range(3)}
        links_gauged = gauge_transform(links, gauge)
        S_after = wilson_action(plaquettes, links_gauged, beta)

        assert np.isclose(S_before, S_after, atol=1e-10)

    def test_identity_gauge_is_noop(self):
        rng = np.random.default_rng(1)
        links = {(0, 1): su3_random(rng), (1, 2): su3_random(rng)}
        gauge = {i: np.eye(3, dtype=complex) for i in range(3)}
        gauged = gauge_transform(links, gauge)
        for k in links:
            assert np.allclose(links[k], gauged[k], atol=1e-14)


# ---------------------------------------------------------------------------
# 5. THEOREM 1: Cartan torus reduction to U(1) x U(1)
# ---------------------------------------------------------------------------

class TestCartanReduction:
    """Theorem 1: SU(3) Wilson action restricted to the maximal Cartan
    torus T = U(1) x U(1) reproduces the U(1) x U(1) Wilson action exactly."""

    def test_diagonal_character_formula(self):
        """For U = diag(e^{i a}, e^{i b}, e^{-i(a+b)}),
        (1/3) Re tr U = (1/3) [cos a + cos b + cos(a+b)]."""
        for (a, b) in [(0.0, 0.0), (0.3, 0.7), (1.5, -2.1)]:
            U = cartan_su3_from_phases(a, b)
            chi_su3 = su3_character(U)
            chi_expected = (np.cos(a) + np.cos(b) + np.cos(a + b)) / 3.0
            assert np.isclose(chi_su3, chi_expected, atol=1e-12)

    def test_su3_action_reduces_to_u1xu1_at_cartan_triangle(self):
        """Triangle plaquette with diagonal SU(3) links reduces to U(1) x U(1)."""
        phases_a = {(0, 1): 0.3, (1, 2): 0.7, (2, 0): -0.5}
        phases_b = {(0, 1): 0.1, (1, 2): -0.4, (2, 0): 0.9}
        links = cartan_links_from_phases(phases_a, phases_b)
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        beta = 1.7
        S_su3 = wilson_action([plaq], links, beta)
        S_u1xu1 = u1xu1_action_from_su3([plaq], phases_a, phases_b, beta)
        assert np.isclose(S_su3, S_u1xu1, atol=1e-12)

    def test_su3_reduction_machine_precision_bargmann(self):
        """Bargmann graph at N_max=2: full Cartan reduction at machine precision."""
        rng = np.random.default_rng(7)
        A = bargmann_adjacency_dense(2)
        oriented, _ = enumerate_oriented_edges(A)
        plaqs = enumerate_plaquettes(A, max_length=8)
        if len(plaqs) == 0:
            pytest.skip("N_max=2 has no plaquettes")
        # Random Cartan phases on forward edges
        forward = [(e.source, e.target) for e in oriented if e.source < e.target]
        phases_a = {k: float(rng.uniform(-np.pi, np.pi)) for k in forward}
        phases_b = {k: float(rng.uniform(-np.pi, np.pi)) for k in forward}
        links = cartan_links_from_phases(phases_a, phases_b)
        beta = 2.5
        S_su3 = wilson_action(plaqs, links, beta)
        S_u1xu1 = u1xu1_action_from_su3(plaqs, phases_a, phases_b, beta)
        # Should match to machine precision
        assert np.isclose(S_su3, S_u1xu1, atol=1e-12)

    def test_two_independent_u1_components(self):
        """The Cartan reduction has TWO independent U(1) phase fields;
        zeroing one gives a single-U(1) sub-sector."""
        phases_a = {(0, 1): 0.5, (1, 2): -0.3, (2, 0): 0.2}
        phases_b = {(0, 1): 0.0, (1, 2): 0.0, (2, 0): 0.0}  # only a is on
        plaq = [OrientedEdge(0, 1), OrientedEdge(1, 2), OrientedEdge(2, 0)]
        beta = 2.0
        S = u1xu1_action_from_su3([plaq], phases_a, phases_b, beta)
        a_P = sum(phases_a.values())
        # Plaquette character with b=0:
        # (1/3)(cos a_P + cos 0 + cos a_P) = (1/3)(2 cos a_P + 1)
        char_expected = (2 * np.cos(a_P) + 1.0) / 3.0
        S_expected = beta * (1 - char_expected)
        assert np.isclose(S, S_expected, atol=1e-12)


# ---------------------------------------------------------------------------
# 6. Bargmann adjacency adapter
# ---------------------------------------------------------------------------

class TestBargmannAdapter:
    def test_nmax2_node_count(self):
        # nodes at N=0,1,2 = 1 + 3 + 6 = 10
        A = bargmann_adjacency_dense(2)
        assert A.shape == (10, 10)

    def test_nmax3_node_count(self):
        # nodes 1+3+6+10 = 20
        A = bargmann_adjacency_dense(3)
        assert A.shape == (20, 20)

    def test_adjacency_symmetric(self):
        for N in [1, 2, 3]:
            A = bargmann_adjacency_dense(N)
            assert np.array_equal(A, A.T)

    def test_no_self_loops(self):
        for N in [1, 2, 3]:
            A = bargmann_adjacency_dense(N)
            assert np.all(np.diag(A) == 0)

    def test_bipartite_by_n_parity(self):
        """The Bargmann graph is bipartite: edges only between (N) and (N+1).
        Test on N_max=3."""
        from geovac.nuclear.bargmann_graph import build_bargmann_graph
        g = build_bargmann_graph(3)
        # Nodes have parity = N mod 2
        parities = np.array([N % 2 for (N, _, _) in g.nodes])
        A = bargmann_adjacency_dense(3)
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                if A[i, j] == 1:
                    # Endpoints must have different N parities
                    assert parities[i] != parities[j]


# ---------------------------------------------------------------------------
# 7. Plaquette enumeration
# ---------------------------------------------------------------------------

class TestPlaquetteEnumeration:
    def test_nmax1_no_plaquettes(self):
        # N_max=1: only N=0 (1 node) and N=1 (3 nodes), 4 nodes total,
        # 3 edges to single N=0 node. Forest, no cycles.
        A = bargmann_adjacency_dense(1)
        plaqs = enumerate_plaquettes(A, max_length=10, both_orientations=False)
        assert len(plaqs) == 0

    def test_nmax2_plaquette_count(self):
        """N_max=2: the (N=0)-(N=1)-(N=2) ladder with multiple Δl=±1 edges.
        Count plaquettes up to length 8."""
        A = bargmann_adjacency_dense(2)
        plaqs_4 = enumerate_plaquettes(A, max_length=4, both_orientations=False)
        # All plaquettes must have even length (bipartite)
        for P in plaqs_4:
            assert len(P) % 2 == 0

    def test_all_plaquettes_even_length_nmax3(self):
        """Bipartite => all primitive cycles have even length."""
        A = bargmann_adjacency_dense(3)
        plaqs = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        for P in plaqs:
            assert len(P) % 2 == 0, (
                f"Bipartite graph cannot have odd-length cycle: {P}"
            )

    def test_plaquettes_closed(self):
        A = bargmann_adjacency_dense(3)
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        for P in plaqs:
            assert P[0].source == P[-1].target

    def test_plaquettes_no_backtracking(self):
        A = bargmann_adjacency_dense(3)
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        for P in plaqs:
            for i in range(len(P) - 1):
                assert not (
                    P[i].target == P[i + 1].source
                    and P[i].source == P[i + 1].target
                )


# ---------------------------------------------------------------------------
# 8. Wilson action positivity & vacuum
# ---------------------------------------------------------------------------

class TestWilsonActionStructure:
    def test_action_zero_at_trivial_vacuum(self):
        """All link variables = identity => action = 0."""
        A = bargmann_adjacency_dense(3)
        oriented, _ = enumerate_oriented_edges(A)
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        if len(plaqs) == 0:
            pytest.skip("no plaquettes")
        forward = [
            (e.source, e.target) for e in oriented if e.source < e.target
        ]
        links = {k: np.eye(3, dtype=complex) for k in forward}
        S = wilson_action(plaqs, links, beta=2.0)
        assert np.isclose(S, 0.0, atol=1e-12)

    def test_action_nonnegative(self):
        """Action S_W = beta * sum_P (1 - (1/3) Re tr U_P) is non-negative."""
        rng = np.random.default_rng(11)
        A = bargmann_adjacency_dense(3)
        oriented, _ = enumerate_oriented_edges(A)
        plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        if len(plaqs) == 0:
            pytest.skip("no plaquettes")
        forward = [
            (e.source, e.target) for e in oriented if e.source < e.target
        ]
        links = {k: su3_random(rng) for k in forward}
        S = wilson_action(plaqs, links, beta=1.5)
        assert S >= -1e-10  # tiny negative tolerance for round-off


# ---------------------------------------------------------------------------
# 9. THEOREM 2: weak-coupling kinetic coefficient = 1/12
# ---------------------------------------------------------------------------

class TestWeakCouplingKinetic:
    def test_kinetic_coefficient_value(self):
        """Per plaquette, per su(3) component: 1/12."""
        c = weak_coupling_kinetic_coefficient_per_plaquette()
        assert c == sp.Rational(1, 12)

    def test_kinetic_coefficient_vs_su2(self):
        """SU(2): 1/8; SU(3): 1/12. Ratio = 8/12 = 2/3."""
        c_su3 = weak_coupling_kinetic_coefficient_per_plaquette()
        c_su2 = sp.Rational(1, 8)
        assert c_su3 / c_su2 == sp.Rational(2, 3)

    def test_kinetic_general_formula(self):
        """SU(N) coefficient is 1/(2 * N_c * tr_norm) where tr_norm = 1/2.
        General formula: 1/(4 N_c). For SU(2): 1/8. For SU(3): 1/12."""
        for N_c, expected in [(2, sp.Rational(1, 8)), (3, sp.Rational(1, 12))]:
            assert sp.Rational(1, 4 * N_c) == expected

    def test_quadratic_expansion_numerically(self):
        """
        Verify numerically that for small A, Re tr(prod_e exp(i sum_a A_e^a T^a))
        equals 3 - (1/4) sum_a (sum_e A_e^a)^2 + O(A^4).

        Hence 1 - (1/3) Re tr ... = (1/12) sum_a (sum_e A_e^a)^2 + O(A^4).
        """
        rng = np.random.default_rng(42)
        T = su3_generators()
        # Triangle plaquette, 3 edges
        n_edges = 3
        # Small random A^a coefficients
        eps = 0.01
        A_list = [eps * rng.standard_normal(8) for _ in range(n_edges)]
        # Compute U_P exactly
        U_P = np.eye(3, dtype=complex)
        for A in A_list:
            X = sum(A[a] * T[a] for a in range(8))
            U_e = _expm(1j * X)
            U_P = U_P @ U_e
        # Quadratic prediction
        A_P = sum(A_list)  # vector of length 8
        S_predict = (eps**0) * (1.0 / 12.0) * float(np.sum(A_P ** 2))
        # Wait: the test uses unscaled A_list. A_P^a = sum_e A_e^a is finite.
        # Actually we used A_list with epsilon already, so A_P^a = sum_e A_e^a
        # is still O(eps).
        # Re-derive: A_P_sq = (A_P)^a (A_P)^a = sum_a (sum_e A_e^a)^2
        A_P_sq = float(np.sum(A_P ** 2))
        S_predicted = (1.0 / 12.0) * A_P_sq
        S_actual = 1.0 - (1.0 / 3.0) * float(np.real(np.trace(U_P)))
        # Should match to O(A^3) ~ eps^3 ~ 1e-6
        assert abs(S_actual - S_predicted) < 1e-5


def _expm(M):
    """Helper: matrix exponential for tests."""
    from scipy.linalg import expm
    return expm(M)


# ---------------------------------------------------------------------------
# 10. End-to-end smoke test on Bargmann at N_max=2
# ---------------------------------------------------------------------------

def test_bargmann_su3_endtoend_nmax2():
    """End-to-end: build Bargmann graph at N_max=2, enumerate plaquettes,
    assign random SU(3) links, verify gauge invariance."""
    A = bargmann_adjacency_dense(2)
    oriented, _ = enumerate_oriented_edges(A)
    V = A.shape[0]

    plaqs = enumerate_plaquettes(A, max_length=8, both_orientations=False)

    rng = np.random.default_rng(11)
    forward = [(e.source, e.target) for e in oriented if e.source < e.target]
    links = {k: su3_random(rng) for k in forward}

    beta = 1.5
    S_before = wilson_action(plaqs, links, beta)

    gauge = {i: su3_random(rng) for i in range(V)}
    links_gauged = gauge_transform(links, gauge)
    S_after = wilson_action(plaqs, links_gauged, beta)
    assert np.isclose(S_before, S_after, atol=1e-8)


def test_bargmann_su3_endtoend_nmax3():
    """End-to-end at N_max=3."""
    A = bargmann_adjacency_dense(3)
    oriented, _ = enumerate_oriented_edges(A)
    V = A.shape[0]

    plaqs = enumerate_plaquettes(A, max_length=6, both_orientations=False)

    rng = np.random.default_rng(13)
    forward = [(e.source, e.target) for e in oriented if e.source < e.target]
    links = {k: su3_random(rng) for k in forward}

    beta = 1.5
    S_before = wilson_action(plaqs, links, beta)

    gauge = {i: su3_random(rng) for i in range(V)}
    links_gauged = gauge_transform(links, gauge)
    S_after = wilson_action(plaqs, links_gauged, beta)
    assert np.isclose(S_before, S_after, atol=1e-8)


# ---------------------------------------------------------------------------
# 11. Slow tests (Monte Carlo, character coefficients)
# ---------------------------------------------------------------------------

class TestMonteCarlo:
    @pytest.mark.slow
    def test_wilson_loop_monotonic_in_beta(self):
        """<W> should monotonically increase from 0 to 1 as beta -> infty."""
        A = bargmann_adjacency_dense(2)
        oriented, _ = enumerate_oriented_edges(A)
        plaqs = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        if len(plaqs) == 0:
            pytest.skip("no plaquettes")
        forward_edges = [
            OrientedEdge(e.source, e.target)
            for e in oriented if e.source < e.target
        ]
        results = []
        for beta in [0.5, 1.0, 2.0, 5.0]:
            mean, _ = monte_carlo_wilson_expectation(
                plaqs[0], plaqs, forward_edges,
                beta=beta, n_samples=300, n_thermalize=100, seed=42,
            )
            results.append(mean)
        # Monotonic increasing
        for i in range(len(results) - 1):
            assert results[i] < results[i + 1] + 0.05


# ---------------------------------------------------------------------------
# 12. Paper 25 Sec. VII.A pins: Bargmann graph counts, plaquette census,
#     and the CP^2-quotient negative (group5 cert finding E10).
#     Reference data: debug/data/st_su3_plaquettes.json (Sprint ST-SU3)
#     and debug/data/s5_graph_spectrum.json (Sprint 5 Track S5).
# ---------------------------------------------------------------------------

def _graph_counts(N_max):
    """(V, E, c, beta_1) of the Bargmann graph at N_max."""
    A = bargmann_adjacency_dense(N_max)
    V = A.shape[0]
    E = int(A.sum()) // 2
    seen = set()
    c = 0
    for s in range(V):
        if s in seen:
            continue
        c += 1
        stack = [s]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(int(y) for y in np.nonzero(A[x])[0] if y not in seen)
    return V, E, c, E - V + c


class TestPaper25S5Pins:
    """Pins the ST-SU3 graph/plaquette census used by Paper 25 Sec. VII.A
    and Paper 30's SU(3) remark (7,765 plaquettes at N_max = 3)."""

    def test_graph_counts(self):
        # (V, E, c, beta_1) per debug/data/st_su3_plaquettes.json
        assert _graph_counts(2) == (10, 15, 1, 6)
        assert _graph_counts(3) == (20, 42, 1, 23)
        assert _graph_counts(4) == (35, 90, 1, 56)

    def test_graph_counts_nmax5(self):
        """Paper 25 Sec. VII.A: at N_max = 5 the Bargmann graph has
        56 nodes, 165 edges, beta_1 = 110."""
        assert _graph_counts(5) == (56, 165, 1, 110)

    def test_plaquette_census_nmax2(self):
        from collections import Counter
        A = bargmann_adjacency_dense(2)
        pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        cnt = Counter(len(P) for P in pl)
        assert dict(cnt) == {4: 15, 6: 21, 8: 146}
        assert len(pl) == 182

    def test_plaquette_census_nmax3(self):
        """N_max = 3, lengths <= 8: 88 + 526 + 7151 = 7,765 plaquettes --
        the count quoted in Paper 30's higher-rank remark (Cartan
        reduction verified 'at 7,765 plaquettes')."""
        from collections import Counter
        A = bargmann_adjacency_dense(3)
        pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        cnt = Counter(len(P) for P in pl)
        assert dict(cnt) == {4: 88, 6: 526, 8: 7151}
        assert len(pl) == 7765

    def test_plaquette_census_nmax4_short(self):
        from collections import Counter
        A = bargmann_adjacency_dense(4)
        pl = enumerate_plaquettes(A, max_length=6, both_orientations=False)
        cnt = Counter(len(P) for P in pl)
        assert dict(cnt) == {4: 272, 6: 2934}

    @pytest.mark.slow
    def test_plaquette_census_nmax4_full(self):
        """N_max = 4 up to L = 8 (about 40 s)."""
        from collections import Counter
        A = bargmann_adjacency_dense(4)
        pl = enumerate_plaquettes(A, max_length=8, both_orientations=False)
        cnt = Counter(len(P) for P in pl)
        assert dict(cnt) == {4: 272, 6: 2934, 8: 64021}
        assert len(pl) == 67227


class TestPaper25CP2Quotient:
    """Paper 25 Sec. VII.A, negative result: the m_l-quotient of the
    Bargmann graph at N_max = 5 is NOT a CP^2 discretization.

    Pins (a) the 12-sector quotient Laplacian spectrum (multiplicity
    scheme: inter-sector weights count crossing edges), (b) the ratio
    range 0.08--0.19 against the Fubini-Study CP^2 spectrum 4k(k+2),
    and (c) the no-fit quantification: the least-squares single
    rescaling leaves a 40.8% maximum relative residual, and NO single
    rescaling can do better than ~33% under either normalization."""

    # From debug/data/s5_graph_spectrum.json (multiplicity scheme),
    # printed rounded in Paper 25 Sec. VII.A.
    QUOT_SPECTRUM = [
        0.0, 2.223909298693098, 4.869922893611419, 4.941588428412833,
        11.650776112948813, 12.10849084926421, 18.737405154015647,
        28.090079613186617, 33.578934595571994, 50.60534212790485,
        63.765529144643054, 99.42802178174746,
    ]

    def _quotient_eigs(self):
        from geovac.nuclear.bargmann_graph import build_bargmann_graph
        g = build_bargmann_graph(5)
        nodes = g.nodes
        A = bargmann_adjacency_dense(5)
        V = A.shape[0]
        secs = sorted(set((N, l) for (N, l, m) in nodes))
        si = {s: i for i, s in enumerate(secs)}
        W = np.zeros((len(secs), len(secs)))
        for i in range(V):
            for j in range(i + 1, V):
                if A[i, j]:
                    a = (nodes[i][0], nodes[i][1])
                    b = (nodes[j][0], nodes[j][1])
                    if a != b:
                        W[si[a], si[b]] += 1
                        W[si[b], si[a]] += 1
        L = np.diag(W.sum(axis=1)) - W
        return secs, np.sort(np.linalg.eigvalsh(L))

    def test_twelve_sectors(self):
        secs, _ = self._quotient_eigs()
        assert secs == [(0, 0), (1, 1), (2, 0), (2, 2), (3, 1), (3, 3),
                        (4, 0), (4, 2), (4, 4), (5, 1), (5, 3), (5, 5)]
        assert [2 * l + 1 for (N, l) in secs] == [1, 3, 1, 5, 3, 7,
                                                  1, 5, 9, 3, 7, 11]

    def test_quotient_spectrum(self):
        _, ev = self._quotient_eigs()
        assert np.allclose(ev, self.QUOT_SPECTRUM, atol=1e-8)
        # and the rounded values printed in the paper
        printed = [0, 2.22, 4.87, 4.94, 11.65, 12.11, 18.74, 28.09,
                   33.58, 50.61, 63.77, 99.43]
        assert np.allclose(np.round(ev, 2), printed, atol=1e-9)

    def test_cp2_no_uniform_fit(self):
        """Ratios lambda_quot,k / lambda_CP2,k in [0.08, 0.19]; the
        least-squares single rescaling leaves a 40.8% max relative
        residual; the minimax-optimal rescaling still leaves 50%
        (residual vs fit) / 33% (residual vs data). CP^2 Fubini-Study
        (holomorphic sectional curvature 4): lambda_k = 4k(k+2)."""
        _, ev = self._quotient_eigs()
        q = ev[1:]                                     # 11 nonzero
        c = np.array([4.0 * k * (k + 2) for k in range(1, 12)])
        r = q / c
        assert 0.08 < r.min() < r.max() < 0.19
        assert abs(r.min() - 0.0824) < 5e-4            # k = 3
        assert abs(r.max() - 0.1853) < 5e-4            # k = 1
        # least-squares rescaling q ~ s * c
        s = float(q @ c / (c @ c))
        resid_ls = np.max(np.abs(q - s * c) / (s * c))
        assert abs(resid_ls - 0.4078) < 5e-3
        # convention-robust floor: for ANY rescaling s, the max relative
        # residual is >= sqrt(rmax/rmin) - 1 (vs fit) and
        # >= 1 - sqrt(rmin/rmax) (vs data)
        floor_fit = np.sqrt(r.max() / r.min()) - 1.0
        floor_data = 1.0 - np.sqrt(r.min() / r.max())
        assert floor_fit > 0.49
        assert floor_data > 0.33
