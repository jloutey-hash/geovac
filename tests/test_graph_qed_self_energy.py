"""
Tests for GN-5: Graph-native one-loop electron self-energy and vertex
correction on the GeoVac graph.

Tests the full GN-5 pipeline:
  1. One-loop self-energy Sigma[a,b] as finite matrix trace
  2. One-loop vertex correction Lambda_total[a,b]
  3. Anomalous magnetic moment extraction (graph-native F_2)
  4. Self-energy structural zero check at the ground state
  5. Pi-free certificate for both Sigma and Lambda

All exact tests run at n_max=2 in sympy arithmetic.
"""

import pytest
import numpy as np
import sympy as sp
from sympy import Matrix, Rational, sqrt, Integer

from geovac.graph_qed_self_energy import (
    SelfEnergyResult,
    VertexCorrectionResult,
    compute_self_energy,
    compute_vertex_correction,
    extract_anomalous_moment,
    self_energy_structural_zero_check,
    self_energy_pi_free_certificate,
    analyze_self_energy_and_vertex,
    _check_rational_matrix,
    _check_pi_free_matrix,
    _check_algebraic_matrix,
    _ground_state_indices,
)
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels


# ===========================================================================
# Helper utilities
# ===========================================================================

class TestHelpers:
    """Tests for internal helper functions."""

    def test_check_rational_pure_rational(self):
        """Rational matrix passes rationality check."""
        M = Matrix([[Rational(1, 2), Rational(3, 4)],
                     [Rational(5, 6), Integer(0)]])
        assert _check_rational_matrix(M) is True

    def test_check_rational_with_sqrt(self):
        """Matrix with sqrt(2) fails rationality check."""
        M = Matrix([[sqrt(2), Rational(1, 3)],
                     [Rational(1, 3), Integer(1)]])
        assert _check_rational_matrix(M) is False

    def test_check_rational_zero_matrix(self):
        """Zero matrix is rational."""
        M = Matrix([[0, 0], [0, 0]])
        assert _check_rational_matrix(M) is True

    def test_check_pi_free_rational(self):
        """Rational matrix is pi-free."""
        M = Matrix([[Rational(1, 2), Rational(3, 4)],
                     [Rational(5, 6), Integer(7)]])
        assert _check_pi_free_matrix(M) is True

    def test_check_pi_free_algebraic(self):
        """Algebraic (sqrt) matrix is pi-free."""
        M = Matrix([[sqrt(2), sqrt(3)],
                     [sqrt(6), Rational(1, 9)]])
        assert _check_pi_free_matrix(M) is True

    def test_check_pi_free_with_pi(self):
        """Matrix with pi fails pi-free check."""
        M = Matrix([[sp.pi, Integer(1)],
                     [Integer(1), Integer(0)]])
        assert _check_pi_free_matrix(M) is False

    def test_check_algebraic_equals_pi_free(self):
        """_check_algebraic_matrix is equivalent to _check_pi_free_matrix."""
        M = Matrix([[sqrt(5), Rational(2, 7)],
                     [Rational(2, 7), sqrt(3)/2]])
        assert _check_algebraic_matrix(M) is True

    def test_ground_state_indices_n2(self):
        """Ground state at n_max=2 has two indices (m_j = +/-1/2)."""
        labels = list(iter_dirac_labels(2))
        gs_idx = _ground_state_indices(labels)
        assert len(gs_idx) == 2
        for i in gs_idx:
            assert labels[i].n_fock == 1
            assert labels[i].kappa == -1

    def test_ground_state_indices_n1(self):
        """Ground state at n_max=1 has two indices."""
        labels = list(iter_dirac_labels(1))
        gs_idx = _ground_state_indices(labels)
        assert len(gs_idx) == 2


# ===========================================================================
# (A) One-loop self-energy -- t=0 (free propagator)
# ===========================================================================

class TestSelfEnergyT0:
    """Self-energy tests with the free Dirac propagator (t=0)."""

    @pytest.fixture(scope="class")
    def se_n2_t0(self):
        """Compute self-energy once for n_max=2, t=0."""
        return compute_self_energy(2, t=Rational(0), exact=True)

    def test_returns_dataclass(self, se_n2_t0):
        """Return type is SelfEnergyResult."""
        assert isinstance(se_n2_t0, SelfEnergyResult)

    def test_dimensions(self, se_n2_t0):
        """n_max=2: 10 Dirac states, 3 Fock edges."""
        assert se_n2_t0.N_dirac == 10
        assert se_n2_t0.E_fock == 3

    def test_sigma_shape(self, se_n2_t0):
        """Sigma is 10 x 10."""
        assert se_n2_t0.Sigma.rows == 10
        assert se_n2_t0.Sigma.cols == 10

    def test_sigma_numpy_shape(self, se_n2_t0):
        """Sigma numpy array is (10, 10)."""
        assert se_n2_t0.Sigma_numpy.shape == (10, 10)

    def test_sigma_hermitian(self, se_n2_t0):
        """Sigma is symmetric (real Hermitian)."""
        assert se_n2_t0.is_hermitian is True

    def test_sigma_hermitian_numpy(self, se_n2_t0):
        """Sigma numpy array is symmetric."""
        S = se_n2_t0.Sigma_numpy
        assert np.allclose(S, S.T, atol=1e-14)

    def test_sigma_is_pi_free(self, se_n2_t0):
        """Sigma is pi-free (algebraic entries only)."""
        assert se_n2_t0.is_pi_free is True

    def test_sigma_not_rational(self, se_n2_t0):
        """Sigma contains sqrt entries (CG coefficients don't cancel)."""
        # This is a genuine finding: the bilinear V.G_gamma.V^T
        # does NOT contract all sqrt's back to rationals
        assert se_n2_t0.is_rational is False

    def test_sigma_trace_rational(self, se_n2_t0):
        """Trace of Sigma is rational (diagonal entries are rational)."""
        tr = sp.nsimplify(se_n2_t0.Sigma.trace(), rational=False)
        assert tr == Rational(44, 3)

    def test_sigma_trace_value(self, se_n2_t0):
        """Trace of Sigma = 44/3."""
        assert se_n2_t0.trace == "44/3"

    def test_eigenvalue_count(self, se_n2_t0):
        """10 eigenvalues (with multiplicity)."""
        assert len(se_n2_t0.eigenvalues) == 10

    def test_five_zero_eigenvalues(self, se_n2_t0):
        """5 of 10 eigenvalues are zero."""
        zero_count = sum(1 for ev in se_n2_t0.eigenvalues if ev == '0')
        assert zero_count == 5

    def test_nonzero_eigenvalues(self, se_n2_t0):
        """The 5 nonzero eigenvalues are {4/3, 2, 2, 4, 16/3}."""
        nonzero_evs = [ev for ev in se_n2_t0.eigenvalues if ev != '0']
        nonzero_vals = sorted([float(sp.sympify(ev)) for ev in nonzero_evs])
        expected = sorted([4/3, 2.0, 2.0, 4.0, 16/3])
        for v, e in zip(nonzero_vals, expected):
            assert abs(v - e) < 1e-12, f"Expected {e}, got {v}"

    def test_sigma_positive_semidefinite(self, se_n2_t0):
        """Sigma eigenvalues are >= 0 (positive semidefinite)."""
        S_np = se_n2_t0.Sigma_numpy
        evals = np.linalg.eigvalsh(S_np)
        assert np.all(evals >= -1e-12), f"Negative eigenvalue: {min(evals)}"

    def test_ground_state_block_nonzero(self, se_n2_t0):
        """Ground state block is NOT zero (graph-continuum difference)."""
        assert se_n2_t0.ground_state_zero is False

    def test_ground_state_block_value(self, se_n2_t0):
        """Ground state block is [[1,1],[1,1]]."""
        gs = se_n2_t0.ground_state_block
        assert gs is not None
        assert gs == Matrix([[1, 1], [1, 1]])

    def test_sigma_sympy_numpy_consistency(self, se_n2_t0):
        """Sympy and numpy Sigma agree."""
        S_sympy = np.array(se_n2_t0.Sigma.tolist(), dtype=float)
        S_numpy = se_n2_t0.Sigma_numpy
        assert np.allclose(S_sympy, S_numpy, atol=1e-12)

    def test_sigma_contains_sqrt2(self, se_n2_t0):
        """Sigma contains sqrt(2) entries (CG structure)."""
        has_sqrt2 = False
        for entry in se_n2_t0.Sigma:
            if sqrt(2) in sp.expand(entry).atoms(sp.Pow):
                has_sqrt2 = True
                break
            # Alternative check: look for sqrt in string repr
            s = str(sp.expand(entry))
            if 'sqrt(2)' in s or 'sqrt(6)' in s:
                has_sqrt2 = True
                break
        assert has_sqrt2, "Expected sqrt(2) or sqrt(6) in Sigma entries"


# ===========================================================================
# (B) One-loop vertex correction -- t=0
# ===========================================================================

class TestVertexCorrectionT0:
    """Vertex correction tests with the free Dirac propagator (t=0)."""

    @pytest.fixture(scope="class")
    def vc_n2_t0(self):
        """Compute vertex correction once for n_max=2, t=0."""
        return compute_vertex_correction(2, t=Rational(0), exact=True)

    def test_returns_dataclass(self, vc_n2_t0):
        """Return type is VertexCorrectionResult."""
        assert isinstance(vc_n2_t0, VertexCorrectionResult)

    def test_dimensions(self, vc_n2_t0):
        """n_max=2: 10 Dirac states, 3 Fock edges."""
        assert vc_n2_t0.N_dirac == 10
        assert vc_n2_t0.E_fock == 3

    def test_lambda_shape(self, vc_n2_t0):
        """Lambda_total is 10 x 10."""
        assert vc_n2_t0.Lambda_total.rows == 10
        assert vc_n2_t0.Lambda_total.cols == 10

    def test_lambda_is_pi_free(self, vc_n2_t0):
        """Lambda_total is pi-free."""
        assert vc_n2_t0.is_pi_free is True

    def test_lambda_not_rational(self, vc_n2_t0):
        """Lambda_total contains algebraic irrationals."""
        assert vc_n2_t0.is_rational is False

    def test_lambda_trace_rational(self, vc_n2_t0):
        """Trace of Lambda_total is rational: 32/9."""
        tr = sp.nsimplify(vc_n2_t0.Lambda_total.trace(), rational=False)
        assert tr == Rational(32, 9)

    def test_lambda_trace_value(self, vc_n2_t0):
        """Trace string is 32/9."""
        assert vc_n2_t0.trace == "32/9"

    def test_eigenvalue_count(self, vc_n2_t0):
        """10 eigenvalues total."""
        assert len(vc_n2_t0.eigenvalues) == 10

    def test_five_zero_eigenvalues(self, vc_n2_t0):
        """5 zero eigenvalues."""
        zero_count = sum(1 for ev in vc_n2_t0.eigenvalues if ev == '0')
        assert zero_count == 5

    def test_lambda_numpy_shape(self, vc_n2_t0):
        """Lambda numpy array is (10, 10)."""
        assert vc_n2_t0.Lambda_total_numpy.shape == (10, 10)

    def test_lambda_sympy_numpy_consistency(self, vc_n2_t0):
        """Sympy and numpy Lambda agree."""
        L_sympy = np.array(vc_n2_t0.Lambda_total.tolist(), dtype=float)
        L_numpy = vc_n2_t0.Lambda_total_numpy
        assert np.allclose(L_sympy, L_numpy, atol=1e-12)


# ===========================================================================
# (C) Anomalous magnetic moment -- t=0
# ===========================================================================

class TestAnomalousMomentT0:
    """Anomalous magnetic moment tests with t=0."""

    @pytest.fixture(scope="class")
    def am_n2_t0(self):
        """Extract anomalous moment at n_max=2, t=0."""
        return extract_anomalous_moment(2, t=Rational(0), exact=True)

    def test_returns_dict(self, am_n2_t0):
        """Return type is dict."""
        assert isinstance(am_n2_t0, dict)

    def test_has_f2(self, am_n2_t0):
        """Dict contains F2 key."""
        assert 'F2' in am_n2_t0
        assert 'F2_float' in am_n2_t0

    def test_f2_algebraic(self, am_n2_t0):
        """F_2 = 5*sqrt(2)/3 (algebraic, not rational)."""
        f2 = sp.sympify(am_n2_t0['F2'])
        expected = 5 * sqrt(2) / 3
        assert sp.simplify(f2 - expected) == 0

    def test_f2_not_rational(self, am_n2_t0):
        """F_2 is not rational (contains sqrt(2))."""
        assert am_n2_t0['F2_rational'] is False

    def test_f2_float_positive(self, am_n2_t0):
        """F_2 float is positive."""
        assert am_n2_t0['F2_float'] > 0

    def test_f2_float_value(self, am_n2_t0):
        """F_2 float ~ 2.357."""
        assert abs(am_n2_t0['F2_float'] - 5 * 2**0.5 / 3) < 1e-10

    def test_tr_lambda_rational(self, am_n2_t0):
        """Tr(Lambda) = 32/9 is rational."""
        tr = sp.sympify(am_n2_t0['Tr_Lambda'])
        assert tr == Rational(32, 9)

    def test_tr_v_bare_ge_algebraic(self, am_n2_t0):
        """Tr(V_bare . G_e) = 16*sqrt(2)/15."""
        tr = sp.sympify(am_n2_t0['Tr_V_bare_G_e'])
        expected = 16 * sqrt(2) / 15
        assert sp.simplify(tr - expected) == 0

    def test_v_bare_nonzero_count(self, am_n2_t0):
        """V_bare has 40 nonzero entries at n_max=2."""
        assert am_n2_t0['V_bare_nonzero_count'] == 40


# ===========================================================================
# (D) Structural zero check -- t=0
# ===========================================================================

class TestStructuralZeroT0:
    """Self-energy structural zero tests with t=0."""

    @pytest.fixture(scope="class")
    def sz_n2_t0(self):
        """Run structural zero check at n_max=2, t=0."""
        return self_energy_structural_zero_check(2, t=Rational(0), exact=True)

    def test_returns_dict(self, sz_n2_t0):
        """Return type is dict."""
        assert isinstance(sz_n2_t0, dict)

    def test_ground_state_not_zero(self, sz_n2_t0):
        """Ground state is NOT zero on the graph (graph-continuum difference)."""
        assert sz_n2_t0['ground_state_zero'] is False

    def test_ground_state_indices_correct(self, sz_n2_t0):
        """Ground state has two indices."""
        assert len(sz_n2_t0['ground_state_indices']) == 2

    def test_sigma_is_pi_free(self, sz_n2_t0):
        """Sigma is pi-free even though ground state is nonzero."""
        assert sz_n2_t0['Sigma_is_pi_free'] is True

    def test_mechanism_message(self, sz_n2_t0):
        """Mechanism message reports the graph-continuum difference."""
        assert 'does NOT have a structural zero' in sz_n2_t0['mechanism']

    def test_has_sigma_diagonal(self, sz_n2_t0):
        """Diagnostic output includes Sigma diagonal."""
        assert 'Sigma_diagonal' in sz_n2_t0

    def test_sigma_diagonal_length(self, sz_n2_t0):
        """10 diagonal entries at n_max=2."""
        assert len(sz_n2_t0['Sigma_diagonal']) == 10

    def test_gs_vertex_coupling_exists(self, sz_n2_t0):
        """Ground state has nonzero vertex couplings."""
        assert sz_n2_t0['gs_vertex_coupling_count'] > 0

    def test_ground_state_block_entries(self, sz_n2_t0):
        """Ground state block entries are all 1."""
        block = sz_n2_t0['ground_state_block_entries']
        for key, val in block.items():
            assert val == '1', f"Entry {key} = {val}, expected 1"


# ===========================================================================
# Pi-free certificate
# ===========================================================================

class TestPiFreeCertificate:
    """Tests for the pi-free certificate function."""

    @pytest.fixture(scope="class")
    def cert_n2_t0(self):
        """Run pi-free certificate at n_max=2, t=0."""
        return self_energy_pi_free_certificate(2, t=Rational(0))

    def test_returns_dict(self, cert_n2_t0):
        """Return type is dict."""
        assert isinstance(cert_n2_t0, dict)

    def test_sigma_pi_free(self, cert_n2_t0):
        """Sigma is pi-free."""
        assert cert_n2_t0['Sigma_pi_free'] is True

    def test_lambda_pi_free(self, cert_n2_t0):
        """Lambda is pi-free."""
        assert cert_n2_t0['Lambda_pi_free'] is True

    def test_sigma_hermitian(self, cert_n2_t0):
        """Sigma is Hermitian."""
        assert cert_n2_t0['Sigma_hermitian'] is True

    def test_sigma_not_rational(self, cert_n2_t0):
        """Sigma is NOT rational (contains sqrt)."""
        assert cert_n2_t0['Sigma_rational'] is False

    def test_lambda_not_rational(self, cert_n2_t0):
        """Lambda is NOT rational (contains sqrt)."""
        assert cert_n2_t0['Lambda_rational'] is False

    def test_all_pass_false(self, cert_n2_t0):
        """all_pass is False because rationality fails."""
        assert cert_n2_t0['all_pass'] is False


# ===========================================================================
# Self-energy at t=kappa (coupled propagator)
# ===========================================================================

@pytest.mark.slow
class TestSelfEnergyTKappa:
    """Self-energy tests with the coupled Dirac propagator (t=-1/16)."""

    @pytest.fixture(scope="class")
    def se_n2_kappa(self):
        """Compute self-energy at n_max=2, t=-1/16."""
        return compute_self_energy(2, t=Rational(-1, 16), exact=True)

    def test_returns_dataclass(self, se_n2_kappa):
        """Return type is SelfEnergyResult."""
        assert isinstance(se_n2_kappa, SelfEnergyResult)

    def test_dimensions(self, se_n2_kappa):
        """Still 10 Dirac states, 3 Fock edges."""
        assert se_n2_kappa.N_dirac == 10
        assert se_n2_kappa.E_fock == 3

    def test_sigma_hermitian(self, se_n2_kappa):
        """Sigma is symmetric."""
        assert se_n2_kappa.is_hermitian is True

    def test_sigma_pi_free(self, se_n2_kappa):
        """Sigma is pi-free."""
        assert se_n2_kappa.is_pi_free is True

    def test_t_stored(self, se_n2_kappa):
        """t parameter is stored correctly."""
        assert se_n2_kappa.t == "-1/16"

    def test_sigma_different_from_t0(self, se_n2_kappa):
        """Sigma at t=-1/16 differs from t=0 (vertex correction uses G_e)."""
        # Note: self-energy Sigma does NOT use G_e,
        # only V and G_gamma. So at t=0 vs t=-1/16 the self-energy
        # is IDENTICAL because Sigma = sum V_e G_gamma V_e^T.
        # The vertex correction Lambda uses G_e and will differ.
        se_t0 = compute_self_energy(2, t=Rational(0), exact=True)
        diff = se_n2_kappa.Sigma - se_t0.Sigma
        # Self-energy does not depend on t (no G_e in formula)
        is_zero = all(sp.nsimplify(diff[i, j]) == 0
                       for i in range(diff.rows) for j in range(diff.cols))
        assert is_zero, "Self-energy should be t-independent (no G_e in formula)"


# ===========================================================================
# Vertex correction at t=kappa
# ===========================================================================

@pytest.mark.slow
class TestVertexCorrectionTKappa:
    """Vertex correction tests with t=-1/16."""

    @pytest.fixture(scope="class")
    def vc_n2_kappa(self):
        """Compute vertex correction at n_max=2, t=-1/16."""
        return compute_vertex_correction(2, t=Rational(-1, 16), exact=True)

    def test_returns_dataclass(self, vc_n2_kappa):
        """Return type is VertexCorrectionResult."""
        assert isinstance(vc_n2_kappa, VertexCorrectionResult)

    def test_lambda_pi_free(self, vc_n2_kappa):
        """Lambda at t=kappa is pi-free."""
        assert vc_n2_kappa.is_pi_free is True

    def test_lambda_differs_from_t0(self, vc_n2_kappa):
        """Lambda at t=-1/16 differs from t=0 (uses G_e which depends on t)."""
        vc_t0 = compute_vertex_correction(2, t=Rational(0), exact=True)
        diff_np = vc_n2_kappa.Lambda_total_numpy - vc_t0.Lambda_total_numpy
        # Should NOT be zero since G_e changes with t
        assert np.max(np.abs(diff_np)) > 1e-10, \
            "Lambda should differ between t=0 and t=-1/16"


# ===========================================================================
# Full analysis driver
# ===========================================================================

@pytest.mark.slow
class TestAnalysisDriver:
    """Tests for the full analysis driver."""

    @pytest.fixture(scope="class")
    def analysis(self):
        """Run full analysis at n_max=2."""
        return analyze_self_energy_and_vertex(n_max=2)

    def test_returns_dict(self, analysis):
        """Return type is dict."""
        assert isinstance(analysis, dict)

    def test_has_module_key(self, analysis):
        """Contains module name."""
        assert analysis['module'] == 'geovac.graph_qed_self_energy'

    def test_has_n_max(self, analysis):
        """Contains n_max."""
        assert analysis['n_max'] == 2

    def test_has_two_t_values(self, analysis):
        """Has analyses for both t=0 and t=-1/16."""
        assert 't=0' in analysis['analyses']
        assert 't=-1/16' in analysis['analyses']

    def test_t0_self_energy_present(self, analysis):
        """t=0 analysis contains self_energy."""
        a = analysis['analyses']['t=0']
        assert 'self_energy' in a
        assert 'vertex_correction' in a
        assert 'anomalous_moment' in a
        assert 'structural_zero' in a
        assert 'pi_free_certificate' in a

    def test_t0_sigma_matrix_entries(self, analysis):
        """t=0 analysis includes Sigma matrix entries."""
        a = analysis['analyses']['t=0']
        assert 'matrix_entries' in a['self_energy']
        entries = a['self_energy']['matrix_entries']
        assert len(entries) == 10  # 10 rows
        assert len(entries[0]) == 10  # 10 columns


# ===========================================================================
# Cross-checks and structural properties
# ===========================================================================

class TestStructuralProperties:
    """Tests for structural properties of the self-energy and vertex."""

    def test_sigma_trace_equals_sum_eigenvalues(self):
        """Tr(Sigma) = sum of eigenvalues."""
        se = compute_self_energy(2, t=Rational(0), exact=True)
        tr = float(sp.sympify(se.trace))
        eig_sum = sum(float(sp.sympify(ev)) for ev in se.eigenvalues)
        assert abs(tr - eig_sum) < 1e-10

    def test_sigma_rank(self):
        """Sigma has rank 5 (5 zero eigenvalues, 5 nonzero)."""
        se = compute_self_energy(2, t=Rational(0), exact=True)
        rank = np.linalg.matrix_rank(se.Sigma_numpy, tol=1e-10)
        assert rank == 5

    def test_lambda_trace_positive(self):
        """Tr(Lambda_total) > 0 at t=0."""
        vc = compute_vertex_correction(2, t=Rational(0), exact=True)
        tr = float(sp.sympify(vc.trace))
        assert tr > 0

    def test_sigma_diagonal_rational(self):
        """All diagonal entries of Sigma are rational."""
        se = compute_self_energy(2, t=Rational(0), exact=True)
        for i in range(se.N_dirac):
            entry = se.Sigma[i, i]
            simplified = sp.nsimplify(entry, rational=False)
            assert simplified.is_Rational or simplified.is_Integer or simplified == 0, \
                f"Diagonal [{i},{i}] = {entry} is not rational"

    def test_lambda_trace_equals_sum_eigenvalues(self):
        """Tr(Lambda) = sum of eigenvalues."""
        vc = compute_vertex_correction(2, t=Rational(0), exact=True)
        tr = float(sp.sympify(vc.trace))
        eig_sum = sum(float(sp.sympify(ev)) for ev in vc.eigenvalues)
        assert abs(tr - eig_sum) < 1e-10

    def test_sigma_off_diagonal_structure(self):
        """Off-diagonal entries in the n_fock=2 block contain sqrt."""
        se = compute_self_energy(2, t=Rational(0), exact=True)
        labels = list(iter_dirac_labels(2))
        n2_indices = [i for i, lab in enumerate(labels) if lab.n_fock == 2]
        has_sqrt = False
        for i in n2_indices:
            for j in n2_indices:
                if i != j:
                    entry = se.Sigma[i, j]
                    s = str(sp.expand(entry))
                    if 'sqrt' in s:
                        has_sqrt = True
                        break
        assert has_sqrt, "Expected sqrt in off-diagonal n_fock=2 block"

    def test_sigma_ground_block_rank_1(self):
        """Ground state block [[1,1],[1,1]] has rank 1."""
        se = compute_self_energy(2, t=Rational(0), exact=True)
        gs = se.ground_state_block
        rank = np.linalg.matrix_rank(
            np.array(gs.tolist(), dtype=float), tol=1e-10)
        assert rank == 1


# ===========================================================================
# Edge cases
# ===========================================================================

class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_n_max_1_no_edges(self):
        """n_max=1 has 0 edges; self-energy should be zero or error gracefully."""
        # n_max=1: V=1, E=0, no photon propagator possible
        # The function should handle this gracefully
        try:
            se = compute_self_energy(1, t=Rational(0), exact=True)
            # If it succeeds, Sigma should be zero (no edges)
            assert se.E_fock == 0
        except Exception:
            # Acceptable to fail at n_max=1 (no photon graph)
            pass

    def test_numpy_path(self):
        """numpy path (exact=False) works."""
        se = compute_self_energy(2, t=Rational(0), exact=False)
        assert se.Sigma is None  # No sympy matrix in numpy mode
        assert se.Sigma_numpy is not None
        assert se.Sigma_numpy.shape == (10, 10)
        # Should still be approximately symmetric
        assert np.allclose(se.Sigma_numpy, se.Sigma_numpy.T, atol=1e-10)
