"""
Tests for Sturmian qubit encoding (Track BU-2).

Validates that:
1. SturmianCI.solve() exposes h1_ortho and eri_4d_ortho
2. JW qubit Hamiltonian eigenvalues match FCI eigenvalues
3. Standard GeoVac encoding reproduces Paper 14 Pauli counts
"""

import numpy as np
import pytest

from geovac.sturmian_solver import SturmianCI

try:
    from openfermion import FermionOperator, jordan_wigner, get_sparse_operator
    HAS_OPENFERMION = True
except ImportError:
    HAS_OPENFERMION = False


def _build_jw_from_integrals(h1, eri_4d, n_spatial, threshold=1e-12):
    """Build JW QubitOperator from orthonormalized integrals."""
    fermion_op = FermionOperator()
    for p in range(n_spatial):
        for q in range(n_spatial):
            val = h1[p, q]
            if abs(val) < threshold:
                continue
            for sigma in range(2):
                fermion_op += FermionOperator(
                    ((2*p+sigma, 1), (2*q+sigma, 0)), val
                )
    for a in range(n_spatial):
        for b in range(n_spatial):
            for c in range(n_spatial):
                for d in range(n_spatial):
                    val = eri_4d[a, b, c, d]
                    if abs(val) < threshold:
                        continue
                    coeff = 0.5 * val
                    for sigma in range(2):
                        for tau in range(2):
                            sp_a = 2 * a + sigma
                            sp_b = 2 * b + tau
                            if sp_a == sp_b:
                                continue
                            fermion_op += FermionOperator(
                                ((sp_a, 1), (sp_b, 1), (2*d+tau, 0), (2*c+sigma, 0)),
                                coeff,
                            )
    return jordan_wigner(fermion_op)


class TestSturmianSolveExposesIntegrals:
    """SturmianCI.solve() must return h1_ortho and eri_4d_ortho."""

    def test_solve_returns_h1_ortho(self):
        sci = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = sci.solve(k=1.812)
        assert 'h1_ortho' in result
        assert isinstance(result['h1_ortho'], np.ndarray)
        assert result['h1_ortho'].shape == (5, 5)

    def test_solve_returns_eri_4d_ortho(self):
        sci = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = sci.solve(k=1.812)
        assert 'eri_4d_ortho' in result
        assert isinstance(result['eri_4d_ortho'], np.ndarray)
        assert result['eri_4d_ortho'].shape == (5, 5, 5, 5)

    def test_h1_ortho_hermitian(self):
        sci = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = sci.solve(k=1.812)
        h1 = result['h1_ortho']
        assert np.allclose(h1, h1.T, atol=1e-12)


@pytest.mark.skipif(not HAS_OPENFERMION, reason="openfermion not installed")
class TestSturmianQubitEigenvalueMatch:
    """JW qubit Hamiltonian must reproduce FCI eigenvalues."""

    def test_qubit_eigenvalue_matches_fci_maxn2(self):
        sci = SturmianCI(Z=2, n_electrons=2, max_n=2)
        result = sci.solve(k=1.812)
        h1_ortho = result['h1_ortho']
        eri_4d_ortho = result['eri_4d_ortho']
        e_fci = result['energy']

        qubit_op = _build_jw_from_integrals(h1_ortho, eri_4d_ortho, sci.n_spatial)
        n_qubits = 2 * sci.n_spatial

        sparse_mat = get_sparse_operator(qubit_op, n_qubits=n_qubits)
        mat = sparse_mat.toarray()
        eigvals = np.linalg.eigvalsh(mat)
        e_qubit = float(eigvals[0])

        assert abs(e_qubit - e_fci) < 1e-8, \
            f"Qubit GS {e_qubit:.8f} != FCI GS {e_fci:.8f}"


@pytest.mark.skipif(not HAS_OPENFERMION, reason="openfermion not installed")
class TestStandardEncodingPaper14:
    """Standard GeoVac encoding must reproduce Paper 14 Pauli counts."""

    def test_pauli_count_maxn2(self):
        from geovac.vqe_benchmark import build_geovac_he
        _, of_op, n_q, _ = build_geovac_he(max_n=2)
        assert len(of_op.terms) == 120

    @pytest.mark.slow
    def test_pauli_count_maxn3(self):
        from geovac.vqe_benchmark import build_geovac_he
        _, of_op, n_q, _ = build_geovac_he(max_n=3)
        assert len(of_op.terms) == 2659
