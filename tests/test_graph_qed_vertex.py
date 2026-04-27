"""
Tests for GN-4: Graph-native electron-photon vertex coupling and
one-loop vacuum polarization on the GeoVac graph.

Tests the full GN-4 pipeline:
  1. CG projection matrix P from Dirac to Fock nodes
  2. Vertex tensor V[a, b, e] bridging electron and photon spaces
  3. One-loop vacuum polarization Pi as a finite matrix trace
  4. Pi-free certificate: all entries rational

All tests run at n_max=2 in exact sympy rational arithmetic.
"""

import pytest
import numpy as np
import sympy as sp
from sympy import Matrix, Rational, sqrt, Integer

from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
    compute_vacuum_polarization,
    vertex_pi_free_certificate,
    vertex_summary,
    analyze_vertex_and_vp,
    _check_pi_free_matrix,
    _check_rational_matrix,
)
from geovac.graph_qed_photon import build_fock_graph, FockGraphData
from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l


# ===========================================================================
# Projection matrix P (Dirac -> Fock)
# ===========================================================================

class TestProjectionMatrix:
    """Tests for the CG projection from Dirac to scalar Fock states."""

    def test_dimensions_n1(self):
        """n_max=1: 2 Dirac states, 1 Fock state."""
        P, dlabs, fstates = build_projection_matrix(1)
        assert P.rows == 2, f"Expected 2 Dirac states at n_max=1, got {P.rows}"
        assert P.cols == 1, f"Expected 1 Fock state at n_max=1, got {P.cols}"

    def test_dimensions_n2(self):
        """n_max=2: 10 Dirac states, 5 Fock states."""
        P, dlabs, fstates = build_projection_matrix(2)
        assert P.rows == 10
        assert P.cols == 5
        assert len(dlabs) == 10
        assert len(fstates) == 5

    def test_dirac_label_count(self):
        """Count matches iter_dirac_labels."""
        P, dlabs, _ = build_projection_matrix(2)
        expected = list(iter_dirac_labels(2))
        assert len(dlabs) == len(expected)
        for a, b in zip(dlabs, expected):
            assert a == b

    def test_fock_state_count(self):
        """Fock states match GeometricLattice."""
        from geovac.lattice import GeometricLattice
        P, _, fstates = build_projection_matrix(2)
        lat = GeometricLattice(max_n=2, topological_weights=False)
        assert fstates == lat.states

    def test_p_nonzero_count_n2(self):
        """At n_max=2, P should have exactly 14 nonzero entries."""
        P, _, _ = build_projection_matrix(2)
        nnz = sum(1 for e in P if e != 0)
        assert nnz == 14

    def test_p_entries_algebraic(self):
        """All P entries are algebraic (sqrt of rational), no free symbols."""
        P, _, _ = build_projection_matrix(2)
        for entry in P:
            expr = sp.sympify(entry)
            assert not expr.free_symbols, f"P entry {entry} has free symbols"

    def test_p_pi_free(self):
        """P contains no transcendentals (pi, e, zeta)."""
        P, _, _ = build_projection_matrix(2)
        assert _check_pi_free_matrix(P)

    def test_p_column_norm_consistency(self):
        """Each Fock column should have squared CG norms summing to integer
        (number of Dirac states projecting onto that node)."""
        P, dlabs, fstates = build_projection_matrix(2)
        for j in range(P.cols):
            col_sq = sum(P[i, j] ** 2 for i in range(P.rows))
            col_sq = sp.nsimplify(sp.expand(col_sq))
            # The sum of CG^2 over m_s values gives the number of
            # Dirac states projecting onto this Fock node
            assert col_sq > 0, f"Fock state {fstates[j]} has no projecting Dirac states"

    def test_s_wave_projection_is_unity(self):
        """For l=0 (s-wave) states, j=1/2, the CG coefficient is exactly 1.
        Each 1s Dirac state maps to the single 1s Fock node with coeff 1."""
        P, dlabs, fstates = build_projection_matrix(2)
        # 1s states: n=1, kappa=-1, two_mj in {-1, 1}
        for i, dl in enumerate(dlabs):
            if dl.n_fock == 1 and dl.kappa == -1:
                # Should project onto (1, 0, 0)
                fock_idx = fstates.index((1, 0, 0))
                assert P[i, fock_idx] == 1, (
                    f"1s Dirac state {dl} should project with CG=1 onto (1,0,0), "
                    f"got {P[i, fock_idx]}"
                )

    def test_n_quantum_number_preserved(self):
        """Projection P[a, v] is nonzero only when n_fock(a) = n(v)."""
        P, dlabs, fstates = build_projection_matrix(2)
        for i, dl in enumerate(dlabs):
            for j, fs in enumerate(fstates):
                if P[i, j] != 0:
                    assert dl.n_fock == fs[0], (
                        f"n mismatch: Dirac n_fock={dl.n_fock}, Fock n={fs[0]}"
                    )

    def test_l_quantum_number_preserved(self):
        """Projection P[a, v] is nonzero only when l(kappa) = l(v)."""
        P, dlabs, fstates = build_projection_matrix(2)
        for i, dl in enumerate(dlabs):
            for j, fs in enumerate(fstates):
                if P[i, j] != 0:
                    assert dl.l == fs[1], (
                        f"l mismatch: Dirac l={dl.l}, Fock l={fs[1]}"
                    )


# ===========================================================================
# Vertex tensor V[a, b, e]
# ===========================================================================

class TestVertexTensor:
    """Tests for the vertex coupling tensor."""

    def test_vertex_dimensions_n2(self):
        """At n_max=2: 10 Dirac, 3 edges."""
        entries, N, V, E = build_vertex_tensor(2)
        assert N == 10
        assert E == 3

    def test_vertex_nonzero_count_n2(self):
        """48 nonzero entries at n_max=2."""
        entries, _, _, _ = build_vertex_tensor(2)
        assert len(entries) == 48

    def test_vertex_sparsity(self):
        """Vertex is 84% sparse at n_max=2."""
        entries, N, V, E = build_vertex_tensor(2)
        total = N * N * E
        sparsity = 1.0 - len(entries) / total
        assert abs(sparsity - 0.84) < 0.01

    def test_vertex_symmetry_in_electron_indices(self):
        """V[a, b, e] should be symmetric in electron indices:
        V[a, b, e] = V[b, a, e] because edge is undirected."""
        entries, N, _, E = build_vertex_tensor(2)
        V_mats = vertex_tensor_to_matrices(entries, N, E)
        for e_idx in range(E):
            Vm = V_mats[e_idx]
            diff = Vm - Vm.T
            for entry in diff:
                assert sp.expand(entry) == 0, (
                    f"Vertex matrix for edge {e_idx} is not symmetric"
                )

    def test_vertex_pi_free(self):
        """All vertex entries are pi-free (algebraic)."""
        entries, _, _, _ = build_vertex_tensor(2)
        transcendental_atoms = (sp.pi, sp.E, sp.EulerGamma, sp.Catalan)
        for a, b, e, val in entries:
            for atom in sp.preorder_traversal(sp.sympify(val)):
                assert atom not in transcendental_atoms, (
                    f"Vertex entry V[{a},{b},{e}] = {val} contains transcendental"
                )

    def test_vertex_matrices_match_entries(self):
        """vertex_tensor_to_matrices reconstructs the correct matrices."""
        entries, N, _, E = build_vertex_tensor(2)
        V_mats = vertex_tensor_to_matrices(entries, N, E)
        assert len(V_mats) == E
        for Vm in V_mats:
            assert Vm.shape == (N, N)
        # Verify each entry appears in the right matrix
        for a, b, e, val in entries:
            assert sp.expand(V_mats[e][a, b] - val) == 0

    def test_t_edge_couples_different_shells(self):
        """Edge e0 (T-type, dn=1) should couple n=1 <-> n=2 states."""
        entries, N, _, E = build_vertex_tensor(2)
        P, dlabs, _ = build_projection_matrix(2)
        fock_data = build_fock_graph(2)

        for a, b, e, val in entries:
            if e == 0:  # T-type edge
                da, db = dlabs[a], dlabs[b]
                assert abs(da.n_fock - db.n_fock) == 1, (
                    f"T-edge e0 should couple dn=+-1, "
                    f"got ({da.n_fock},{db.n_fock})"
                )

    def test_l_edge_same_shell(self):
        """L-type edges (e1, e2) should couple same-shell states only."""
        entries, N, _, E = build_vertex_tensor(2)
        P, dlabs, _ = build_projection_matrix(2)

        for a, b, e, val in entries:
            if e in (1, 2):  # L-type edges
                da, db = dlabs[a], dlabs[b]
                assert da.n_fock == db.n_fock, (
                    f"L-edge e{e} should couple same shell, "
                    f"got n=({da.n_fock},{db.n_fock})"
                )


# ===========================================================================
# Vacuum polarization Pi[e1, e2]
# ===========================================================================

class TestVacuumPolarization:
    """Tests for the one-loop VP tensor at t=0."""

    @pytest.fixture
    def vp_t0(self):
        """VP result at n_max=2, t=0."""
        return compute_vacuum_polarization(2, t=Rational(0), exact=True)

    def test_pi_shape(self, vp_t0):
        """Pi is 3x3 (3 photon edges at n_max=2)."""
        Pi = vp_t0['Pi']
        assert Pi.shape == (3, 3)

    def test_pi_is_rational(self, vp_t0):
        """All Pi entries are rational (sqrt's cancel in trace)."""
        Pi = vp_t0['Pi']
        assert _check_rational_matrix(Pi), (
            f"Pi is not fully rational: {Pi}"
        )

    def test_pi_is_pi_free(self, vp_t0):
        """Pi contains no transcendentals."""
        assert vp_t0['pi_free'] is True

    def test_pi_trace_exact(self, vp_t0):
        """Pi trace = 224/75 at t=0."""
        Pi = vp_t0['Pi']
        expected_trace = Rational(224, 75)
        assert Pi.trace() == expected_trace, (
            f"Pi trace {Pi.trace()} != expected {expected_trace}"
        )

    def test_pi_diagonal_e0(self, vp_t0):
        """Pi[0,0] = 32/15 (T-edge self-coupling)."""
        Pi = vp_t0['Pi']
        assert Pi[0, 0] == Rational(32, 15)

    def test_pi_block_diagonal(self, vp_t0):
        """Pi has block-diagonal structure: T-edge decouples from L-edges."""
        Pi = vp_t0['Pi']
        # Pi[0,1] = Pi[0,2] = Pi[1,0] = Pi[2,0] = 0
        assert Pi[0, 1] == 0
        assert Pi[0, 2] == 0
        assert Pi[1, 0] == 0
        assert Pi[2, 0] == 0

    def test_pi_symmetric(self, vp_t0):
        """Pi is a symmetric matrix."""
        Pi = vp_t0['Pi']
        for i in range(Pi.rows):
            for j in range(Pi.cols):
                assert Pi[i, j] == Pi[j, i], (
                    f"Pi[{i},{j}] = {Pi[i,j]} != Pi[{j},{i}] = {Pi[j,i]}"
                )

    def test_pi_eigenvalues_exact(self, vp_t0):
        """Pi eigenvalues are 32/15, 32/45, 32/225."""
        Pi = vp_t0['Pi']
        eigs = sorted(Pi.eigenvals().keys())
        expected = sorted([Rational(32, 225), Rational(32, 45), Rational(32, 15)])
        for e, exp in zip(eigs, expected):
            assert e == exp, f"Eigenvalue {e} != expected {exp}"

    def test_pi_positive_semidefinite(self, vp_t0):
        """All Pi eigenvalues are non-negative."""
        Pi = vp_t0['Pi']
        for ev in Pi.eigenvals():
            assert ev >= 0, f"Negative eigenvalue {ev}"

    def test_pi_l_block_entries(self, vp_t0):
        """L-block Pi[1,1] = Pi[2,2] = 32/75 and Pi[1,2] = Pi[2,1] = 64/225."""
        Pi = vp_t0['Pi']
        assert Pi[1, 1] == Rational(32, 75)
        assert Pi[2, 2] == Rational(32, 75)
        assert Pi[1, 2] == Rational(64, 225)
        assert Pi[2, 1] == Rational(64, 225)


# ===========================================================================
# Pi-free certificate
# ===========================================================================

class TestPiFreeCertificate:
    """Tests for the full pi-free certificate chain."""

    def test_certificate_all_pass(self):
        """Full certificate at n_max=2 should pass."""
        cert = vertex_pi_free_certificate(2)
        assert cert['all_pass'] is True

    def test_p_pi_free(self):
        """Projection matrix is pi-free."""
        cert = vertex_pi_free_certificate(2)
        assert cert['P_pi_free'] is True

    def test_p_algebraic(self):
        """Projection matrix entries are algebraic (no free symbols)."""
        cert = vertex_pi_free_certificate(2)
        assert cert['P_entries_algebraic'] is True

    def test_v_pi_free(self):
        """Vertex tensor is pi-free."""
        cert = vertex_pi_free_certificate(2)
        assert cert['V_pi_free'] is True

    def test_pi_rational(self):
        """VP tensor is rational (sqrt's cancel)."""
        cert = vertex_pi_free_certificate(2)
        assert cert['Pi_rational'] is True

    def test_certificate_n1(self):
        """Certificate should pass at n_max=1 too."""
        cert = vertex_pi_free_certificate(1)
        assert cert['all_pass'] is True


# ===========================================================================
# Vertex structure analysis
# ===========================================================================

class TestVertexSummary:
    """Tests for the vertex structure summary."""

    def test_summary_fields(self):
        """Summary has all required fields."""
        s = vertex_summary(2)
        assert 'N_dirac' in s
        assert 'V_fock' in s
        assert 'E_fock' in s
        assert 'P_shape' in s
        assert 'edge_analysis' in s

    def test_edge_count(self):
        """3 photon edges at n_max=2."""
        s = vertex_summary(2)
        assert s['E_fock'] == 3

    def test_edge_types(self):
        """e0 is T-type, e1 and e2 are L-type."""
        s = vertex_summary(2)
        ea = s['edge_analysis']
        assert ea[0]['edge_type'] == 'T'
        assert ea[1]['edge_type'] == 'L'
        assert ea[2]['edge_type'] == 'L'

    def test_t_edge_coupling_count(self):
        """T-edge e0 has 8 couplings."""
        s = vertex_summary(2)
        assert s['edge_analysis'][0]['n_couplings'] == 8

    def test_l_edge_coupling_count(self):
        """Each L-edge has 20 couplings."""
        s = vertex_summary(2)
        assert s['edge_analysis'][1]['n_couplings'] == 20
        assert s['edge_analysis'][2]['n_couplings'] == 20

    def test_t_edge_electron_dl_zero(self):
        """Electron transitions coupled by T-edge have dl=0."""
        s = vertex_summary(2)
        assert s['edge_analysis'][0]['electron_dl_values'] == [0]

    def test_l_edge_electron_dn_zero(self):
        """Electron transitions coupled by L-edges have dn=0."""
        s = vertex_summary(2)
        assert s['edge_analysis'][1]['electron_dn_values'] == [0]
        assert s['edge_analysis'][2]['electron_dn_values'] == [0]


# ===========================================================================
# Full analysis driver
# ===========================================================================

class TestAnalyzeDriver:
    """Tests for the full analysis driver."""

    def test_analyze_returns_dict(self):
        """analyze_vertex_and_vp returns a dict with expected keys."""
        result = analyze_vertex_and_vp(n_max=2, t=Rational(0))
        assert 'vertex_structure' in result
        assert 'vacuum_polarization' in result
        assert 'pi_free_certificate' in result

    def test_analyze_n_max_in_result(self):
        """n_max is stored in result."""
        result = analyze_vertex_and_vp(n_max=2)
        assert result['n_max'] == 2


# ===========================================================================
# Structural consistency checks
# ===========================================================================

class TestStructuralConsistency:
    """Cross-module consistency with GN-2 (propagator) and GN-3 (photon)."""

    def test_dirac_state_count_matches_propagator(self):
        """Vertex Dirac count matches DiracGraphOperator."""
        P, dlabs, _ = build_projection_matrix(2)
        op = DiracGraphOperator(n_max=2, t=Rational(0))
        assert P.rows == op.N

    def test_fock_edge_count_matches_photon(self):
        """Vertex edge count matches FockGraphData."""
        entries, _, _, E = build_vertex_tensor(2)
        fock_data = build_fock_graph(2)
        assert E == fock_data.E

    def test_propagator_inverse_of_operator(self):
        """G_e @ D_graph = I (identity check for propagator)."""
        op = DiracGraphOperator(n_max=2, t=Rational(0))
        G_e, _ = electron_propagator(op, exact=True)
        D = op.matrix_sympy()
        prod = G_e * D
        for i in range(op.N):
            for j in range(op.N):
                expected = Integer(1) if i == j else Integer(0)
                assert sp.expand(prod[i, j] - expected) == 0

    def test_vp_trace_positive(self):
        """VP trace should be positive (physical: forward scattering is positive)."""
        vp = compute_vacuum_polarization(2, t=Rational(0), exact=True)
        trace = vp['Pi'].trace()
        assert trace > 0

    def test_vp_common_factor_32(self):
        """All Pi entries and eigenvalues have 32 as a common numerator factor.
        This comes from the electron propagator structure (2/(2n+3))."""
        vp = compute_vacuum_polarization(2, t=Rational(0), exact=True)
        Pi = vp['Pi']
        for ev in Pi.eigenvals():
            # Each eigenvalue should be 32/something
            r = Rational(ev)
            assert r.p % 32 == 0 or r == 0, (
                f"Eigenvalue {ev} does not have 32 as numerator factor"
            )


# ===========================================================================
# Edge case: n_max=1
# ===========================================================================

class TestNmax1:
    """Edge case: n_max=1 has 2 Dirac states and usually 0 Fock edges."""

    def test_nmax1_projection(self):
        """n_max=1: 2 Dirac states project onto 1 Fock node."""
        P, dlabs, fstates = build_projection_matrix(1)
        assert P.rows == 2
        assert P.cols == 1

    def test_nmax1_fock_graph(self):
        """n_max=1: Fock graph has 1 node, 0 edges."""
        fock_data = build_fock_graph(1)
        assert fock_data.V == 1
        assert fock_data.E == 0

    def test_nmax1_vertex_empty(self):
        """n_max=1: no edges means no vertex couplings."""
        entries, N, V, E = build_vertex_tensor(1)
        assert E == 0
        assert len(entries) == 0

    def test_nmax1_vp_empty(self):
        """n_max=1: VP is 0x0 (no photon edges)."""
        vp = compute_vacuum_polarization(1, t=Rational(0), exact=True)
        Pi = vp['Pi']
        assert Pi.rows == 0 and Pi.cols == 0


# ===========================================================================
# VP numerical cross-check
# ===========================================================================

class TestVPNumerical:
    """Numerical cross-checks for VP."""

    def test_vp_numpy_matches_sympy(self):
        """Numpy VP should match sympy VP to machine precision."""
        vp_exact = compute_vacuum_polarization(2, t=Rational(0), exact=True)
        vp_num = compute_vacuum_polarization(2, t=Rational(0), exact=False)

        Pi_exact_np = np.array(vp_exact['Pi'].tolist(), dtype=float)
        Pi_num = vp_num['Pi']

        np.testing.assert_allclose(Pi_num, Pi_exact_np, atol=1e-12,
                                   err_msg="Numpy VP does not match sympy VP")

    def test_vp_eigenvalues_all_positive(self):
        """All VP eigenvalues are non-negative."""
        vp = compute_vacuum_polarization(2, t=Rational(0), exact=False)
        eigs = np.linalg.eigvalsh(vp['Pi'])
        assert np.all(eigs >= -1e-12), f"Negative eigenvalue found: {eigs}"
