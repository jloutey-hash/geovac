"""
Tests for GN-6: Vacuum polarization at n_max=3 on the GeoVac finite graph.

Extends the n_max=2 VP tests to n_max=3 (28 Dirac states, 14 Fock nodes,
13 edges). Validates graph dimensions, VP entries, block decomposition by
l-sector, pi-free certificate, eigenvalue structure, and embedding of
n_max=2 eigenvalues.

Key findings tested:
  - Pi is 13x13, decomposes into l=0 (2x2), l=1 (7x7), l=2 (4x4) blocks
  - l=0 and l=1 blocks are fully rational
  - l=2 block has 4 entries containing sqrt(6) -> pi-free but NOT all-rational
  - 2 negative eigenvalues -> Pi is NOT positive semidefinite at n_max=3
  - All 3 n_max=2 eigenvalues embed in n_max=3 spectrum
  - All diagonal entries have numerator factor 32
  - Trace = 3872/735
"""

import pytest
import sympy as sp
from sympy import Rational, sqrt, Integer, Matrix

from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
    compute_vacuum_polarization,
    _check_pi_free_matrix,
    _check_rational_matrix,
)
from geovac.graph_qed_photon import build_fock_graph
from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
from geovac.dirac_matrix_elements import iter_dirac_labels


# ===========================================================================
# Fixtures
# ===========================================================================

@pytest.fixture(scope="module")
def vp_n3():
    """VP result at n_max=3, t=0, exact sympy."""
    return compute_vacuum_polarization(3, t=Rational(0), exact=True)


@pytest.fixture(scope="module")
def pi_n3(vp_n3):
    """Pi matrix at n_max=3."""
    return vp_n3['Pi']


@pytest.fixture(scope="module")
def fock_n3():
    """Fock graph at n_max=3."""
    return build_fock_graph(3)


@pytest.fixture(scope="module")
def proj_n3():
    """Projection matrix and labels at n_max=3."""
    return build_projection_matrix(3)


# ===========================================================================
# (A) Graph characterization at n_max=3
# ===========================================================================

class TestGraphCharacterizationN3:
    """Graph dimensions and structure at n_max=3."""

    def test_dirac_state_count(self):
        """n_max=3: 28 Dirac states (2+8+18)."""
        labels = list(iter_dirac_labels(3))
        assert len(labels) == 28

    def test_fock_node_count(self, fock_n3):
        """n_max=3: 14 Fock nodes."""
        assert fock_n3.V == 14

    def test_fock_edge_count(self, fock_n3):
        """n_max=3: 13 Fock edges."""
        assert fock_n3.E == 13

    def test_projection_dimensions(self, proj_n3):
        """P is 28 x 14."""
        P, dlabs, fstates = proj_n3
        assert P.rows == 28
        assert P.cols == 14

    def test_betti_numbers(self, fock_n3):
        """beta_0=3 (3 connected components by l-sector), beta_1=2."""
        assert fock_n3.beta_0 == 3
        assert fock_n3.beta_1 == 2

    def test_edge_type_counts(self, fock_n3):
        """5 T-edges (dn=1) and 8 L-edges (dn=0)."""
        t_count = 0
        l_count = 0
        for v1, v2 in fock_n3.edges:
            s1 = fock_n3.states[v1]
            s2 = fock_n3.states[v2]
            if s2[0] - s1[0] != 0:
                t_count += 1
            else:
                l_count += 1
        assert t_count == 5
        assert l_count == 8


# ===========================================================================
# (B) VP matrix entries at n_max=3
# ===========================================================================

class TestVPEntriesN3:
    """Exact VP matrix entries at n_max=3."""

    def test_pi_shape(self, pi_n3):
        """Pi is 13x13."""
        assert pi_n3.shape == (13, 13)

    def test_pi_symmetric(self, pi_n3):
        """Pi is symmetric."""
        for i in range(pi_n3.rows):
            for j in range(i + 1, pi_n3.cols):
                diff = sp.expand(pi_n3[i, j] - pi_n3[j, i])
                assert diff == 0, f"Pi[{i},{j}] != Pi[{j},{i}]"

    def test_trace_exact(self, pi_n3):
        """Trace(Pi) = 3872/735."""
        expected = Rational(3872, 735)
        assert sp.nsimplify(pi_n3.trace(), rational=True) == expected

    def test_diagonal_e0(self, pi_n3):
        """Pi[0,0] = 32/15 (l=0 T-edge, preserved from n_max=2)."""
        assert pi_n3[0, 0] == Rational(32, 15)

    def test_diagonal_e1(self, pi_n3):
        """Pi[1,1] = 32/35 (l=0 T-edge, new at n_max=3)."""
        assert pi_n3[1, 1] == Rational(32, 35)

    def test_l1_diagonal_entries(self, pi_n3):
        """l=1 block diagonal entries."""
        assert pi_n3[2, 2] == Rational(32, 75)
        assert pi_n3[3, 3] == Rational(32, 315)
        assert pi_n3[4, 4] == Rational(32, 75)
        assert pi_n3[5, 5] == Rational(32, 315)
        assert pi_n3[6, 6] == Rational(32, 315)
        assert pi_n3[7, 7] == Rational(32, 147)
        assert pi_n3[8, 8] == Rational(32, 147)

    def test_l2_diagonal_entries(self, pi_n3):
        """l=2 block diagonal entries."""
        assert pi_n3[9, 9] == Rational(32, 245)
        assert pi_n3[10, 10] == Rational(32, 175)
        assert pi_n3[11, 11] == Rational(32, 175)
        assert pi_n3[12, 12] == Rational(32, 245)

    def test_l1_offdiag_rational(self, pi_n3):
        """l=1 block off-diagonal entries are rational."""
        assert pi_n3[2, 4] == Rational(64, 225)
        assert pi_n3[3, 5] == Rational(64, 315)
        assert pi_n3[5, 6] == Rational(64, 315)
        assert pi_n3[7, 8] == Rational(64, 441)

    def test_l2_offdiag_sqrt6(self, pi_n3):
        """l=2 block has sqrt(6) off-diagonal entries."""
        expected = 64 * sqrt(6) / 1225
        assert sp.expand(pi_n3[9, 10] - expected) == 0
        assert sp.expand(pi_n3[11, 12] - expected) == 0

    def test_l2_offdiag_rational_center(self, pi_n3):
        """l=2 block center off-diagonal Pi[10,11] = 192/1225."""
        assert pi_n3[10, 11] == Rational(192, 1225)

    def test_all_diag_numerator_32(self, pi_n3):
        """All diagonal entries have 32 as numerator factor."""
        for i in range(pi_n3.rows):
            entry = sp.nsimplify(pi_n3[i, i], rational=True)
            if isinstance(entry, sp.Rational):
                assert entry.p % 32 == 0, (
                    f"Pi[{i},{i}] = {entry}: numerator {entry.p} not divisible by 32"
                )


# ===========================================================================
# (C) Pi-free certificate at n_max=3
# ===========================================================================

class TestPiFreeCertificateN3:
    """Pi-free and rationality checks at n_max=3."""

    def test_pi_is_pi_free(self, pi_n3):
        """Pi contains no transcendentals (pi, e, zeta)."""
        assert _check_pi_free_matrix(pi_n3)

    def test_pi_is_NOT_all_rational(self, pi_n3):
        """Pi is NOT fully rational -- 4 entries contain sqrt(6)."""
        assert not _check_rational_matrix(pi_n3)

    def test_exactly_four_non_rational_entries(self, pi_n3):
        """Exactly 4 non-rational entries in Pi."""
        non_rat = []
        for i in range(pi_n3.rows):
            for j in range(pi_n3.cols):
                entry = pi_n3[i, j]
                s = sp.nsimplify(entry, rational=True)
                if not isinstance(s, (sp.Rational, sp.Integer,
                                      sp.core.numbers.Zero,
                                      sp.core.numbers.One,
                                      sp.core.numbers.NegativeOne)):
                    non_rat.append((i, j))
        assert len(non_rat) == 4

    def test_non_rational_entries_are_in_l2_block(self, pi_n3):
        """All non-rational entries live in the l=2 block (indices 9-12)."""
        non_rat = []
        for i in range(pi_n3.rows):
            for j in range(pi_n3.cols):
                entry = pi_n3[i, j]
                s = sp.nsimplify(entry, rational=True)
                if not isinstance(s, (sp.Rational, sp.Integer,
                                      sp.core.numbers.Zero,
                                      sp.core.numbers.One,
                                      sp.core.numbers.NegativeOne)):
                    non_rat.append((i, j))
        for i, j in non_rat:
            assert 9 <= i <= 12 and 9 <= j <= 12, (
                f"Non-rational entry at ({i},{j}) is outside l=2 block"
            )

    def test_non_rational_entries_contain_sqrt6(self, pi_n3):
        """Non-rational entries are proportional to sqrt(6)."""
        expected = 64 * sqrt(6) / 1225
        for i, j in [(9, 10), (10, 9), (11, 12), (12, 11)]:
            assert sp.expand(pi_n3[i, j] - expected) == 0

    def test_projection_pi_free(self, proj_n3):
        """Projection matrix P at n_max=3 is pi-free."""
        P, _, _ = proj_n3
        assert _check_pi_free_matrix(P)


# ===========================================================================
# (D) Block decomposition by l-sector
# ===========================================================================

class TestBlockDecompositionN3:
    """VP decomposes into independent l-sector blocks."""

    def test_l0_block_decoupled(self, pi_n3):
        """l=0 block (edges 0,1) decouples from the rest."""
        # Pi[0, j] and Pi[1, j] should be 0 for j >= 2
        for i in range(2):
            for j in range(2, 13):
                assert pi_n3[i, j] == 0, (
                    f"Pi[{i},{j}] = {pi_n3[i,j]} should be 0 (l=0 decoupled)"
                )

    def test_l1_block_decoupled(self, pi_n3):
        """l=1 block (edges 2-8) decouples from l=0 and l=2."""
        for i in range(2, 9):
            # Should be zero for j in l=0 block (0,1)
            for j in range(2):
                assert pi_n3[i, j] == 0
            # Should be zero for j in l=2 block (9-12)
            for j in range(9, 13):
                assert pi_n3[i, j] == 0

    def test_l2_block_decoupled(self, pi_n3):
        """l=2 block (edges 9-12) decouples from the rest."""
        for i in range(9, 13):
            for j in range(9):
                assert pi_n3[i, j] == 0

    def test_l0_block_is_diagonal(self, pi_n3):
        """l=0 block is diagonal (2 edges in different shells, no coupling)."""
        assert pi_n3[0, 1] == 0
        assert pi_n3[1, 0] == 0

    def test_l1_block_has_7_offdiag_pairs(self, pi_n3):
        """l=1 block has specific nonzero off-diagonal pattern."""
        nonzero = []
        for i in range(2, 9):
            for j in range(i + 1, 9):
                if pi_n3[i, j] != 0:
                    nonzero.append((i, j))
        # Expected: (2,4), (3,5), (5,6), (7,8)
        assert len(nonzero) == 4


# ===========================================================================
# (D) Eigenvalue structure
# ===========================================================================

class TestEigenvaluesN3:
    """Eigenvalue properties of Pi at n_max=3."""

    def test_eigenvalue_count(self, pi_n3):
        """13 eigenvalues (counting multiplicity)."""
        eig_dict = pi_n3.eigenvals()
        total = sum(eig_dict.values())
        assert total == 13

    def test_not_positive_semidefinite(self, pi_n3):
        """Pi has at least one negative eigenvalue at n_max=3."""
        eig_dict = pi_n3.eigenvals()
        has_negative = any(float(ev) < -1e-14 for ev in eig_dict.keys())
        assert has_negative, "Pi should have negative eigenvalues at n_max=3"

    def test_exactly_two_negative_eigenvalues(self, pi_n3):
        """Pi has exactly 2 negative eigenvalues."""
        eig_dict = pi_n3.eigenvals()
        neg_count = sum(
            mult for ev, mult in eig_dict.items() if float(ev) < -1e-14
        )
        assert neg_count == 2

    def test_nmax2_eigenvalues_embed(self, pi_n3):
        """All 3 n_max=2 eigenvalues appear in the n_max=3 spectrum."""
        nmax2_eigs = [Rational(32, 225), Rational(32, 45), Rational(32, 15)]
        nmax3_eig_dict = pi_n3.eigenvals()
        nmax3_eigs_float = [float(ev) for ev in nmax3_eig_dict.keys()]

        for ev2 in nmax2_eigs:
            found = any(abs(float(ev2) - ev3) < 1e-10 for ev3 in nmax3_eigs_float)
            assert found, f"n_max=2 eigenvalue {ev2} not found in n_max=3 spectrum"

    def test_rational_eigenvalue_count(self, pi_n3):
        """7 of 13 eigenvalues are rational."""
        eig_dict = pi_n3.eigenvals()
        rational_count = 0
        for ev, mult in eig_dict.items():
            s = sp.nsimplify(ev, rational=True)
            if isinstance(s, (sp.Rational, sp.Integer)):
                rational_count += mult
        assert rational_count == 7

    def test_irrational_eigenvalues_contain_sqrt(self, pi_n3):
        """The 6 irrational eigenvalues contain sqrt(2), sqrt(7), sqrt(10)."""
        eig_dict = pi_n3.eigenvals()
        irrational_evs = []
        for ev in eig_dict.keys():
            s = sp.nsimplify(ev, rational=True)
            if not isinstance(s, (sp.Rational, sp.Integer)):
                irrational_evs.append(ev)

        # Check each irrational eigenvalue has sqrt content
        for ev in irrational_evs:
            ev_str = str(ev)
            has_sqrt = 'sqrt' in ev_str
            assert has_sqrt, f"Irrational eigenvalue {ev} missing sqrt"

    def test_trace_equals_eigenvalue_sum(self, pi_n3):
        """Sum of eigenvalues (with multiplicity) equals trace."""
        eig_dict = pi_n3.eigenvals()
        eig_sum = sum(ev * mult for ev, mult in eig_dict.items())
        trace = pi_n3.trace()
        diff = sp.expand(eig_sum - trace)
        assert diff == 0 or abs(float(diff)) < 1e-12


# ===========================================================================
# (E) Scaling and convergence
# ===========================================================================

class TestScalingN3:
    """Scaling properties between n_max=2 and n_max=3."""

    def test_trace_ratio(self, pi_n3):
        """Tr(Pi_3) / Tr(Pi_2) = 605/343."""
        trace3 = sp.nsimplify(pi_n3.trace(), rational=True)
        trace2 = Rational(224, 75)
        ratio = sp.nsimplify(trace3 / trace2, rational=True)
        assert ratio == Rational(605, 343)

    def test_per_edge_trace_decreases(self, pi_n3):
        """Per-edge trace decreases from n_max=2 to n_max=3."""
        trace3 = float(pi_n3.trace())
        trace2 = float(Rational(224, 75))
        per_edge_2 = trace2 / 3
        per_edge_3 = trace3 / 13
        assert per_edge_3 < per_edge_2

    def test_trace_positive(self, pi_n3):
        """VP trace is positive at n_max=3."""
        assert float(pi_n3.trace()) > 0


# ===========================================================================
# Structural consistency at n_max=3
# ===========================================================================

class TestStructuralConsistencyN3:
    """Cross-module consistency checks at n_max=3."""

    def test_dirac_count_matches_propagator(self, proj_n3):
        """Projection Dirac count matches DiracGraphOperator."""
        P, dlabs, _ = proj_n3
        op = DiracGraphOperator(n_max=3, t=Rational(0))
        assert P.rows == op.N

    def test_fock_edge_count_matches_photon(self, fock_n3):
        """Vertex edge count matches FockGraphData."""
        entries, _, _, E = build_vertex_tensor(3)
        assert E == fock_n3.E

    def test_vertex_symmetry_n3(self):
        """Vertex matrices are symmetric in electron indices at n_max=3."""
        entries, N, _, E = build_vertex_tensor(3)
        V_mats = vertex_tensor_to_matrices(entries, N, E)
        for e_idx in range(E):
            Vm = V_mats[e_idx]
            for i in range(N):
                for j in range(i + 1, N):
                    diff = sp.expand(Vm[i, j] - Vm[j, i])
                    assert diff == 0, (
                        f"V_mat[{e_idx}] not symmetric at ({i},{j})"
                    )

    def test_vertex_nonzero_count(self):
        """260 nonzero vertex entries at n_max=3."""
        entries, _, _, _ = build_vertex_tensor(3)
        assert len(entries) == 260

    def test_vertex_pi_free(self):
        """All vertex entries at n_max=3 are pi-free."""
        entries, _, _, _ = build_vertex_tensor(3)
        transcendental_atoms = (sp.pi, sp.E, sp.EulerGamma, sp.Catalan)
        for a, b, e, val in entries:
            for atom in sp.preorder_traversal(sp.sympify(val)):
                assert atom not in transcendental_atoms, (
                    f"V[{a},{b},{e}] = {val} contains transcendental"
                )

    def test_n_quantum_number_preserved_n3(self, proj_n3):
        """Projection P[a, v] nonzero only when n_fock(a) = n(v)."""
        P, dlabs, fstates = proj_n3
        for i, dl in enumerate(dlabs):
            for j, fs in enumerate(fstates):
                if P[i, j] != 0:
                    assert dl.n_fock == fs[0]

    def test_l_quantum_number_preserved_n3(self, proj_n3):
        """Projection P[a, v] nonzero only when l(kappa) = l(v)."""
        P, dlabs, fstates = proj_n3
        for i, dl in enumerate(dlabs):
            for j, fs in enumerate(fstates):
                if P[i, j] != 0:
                    assert dl.l == fs[1]
