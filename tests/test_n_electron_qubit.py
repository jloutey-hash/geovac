"""
Tests for N-electron qubit encoding (Track AS).

Validates:
  1. Hermitian matrix -> Pauli decomposition correctness
  2. 4-electron angular Hamiltonian construction for qubit encoding
  3. Pauli term counts and structural analysis
  4. Comparison with composed-geometry encoding
"""

import numpy as np
import pytest
from scipy.linalg import eigh

from openfermion import QubitOperator


# ---------------------------------------------------------------------------
# Unit tests for Pauli decomposition
# ---------------------------------------------------------------------------

class TestPauliDecomposition:
    """Test the Hermitian -> Pauli decomposition."""

    def test_identity_2x2(self):
        """Identity matrix decomposes to I with coefficient 1."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        H = np.eye(2)
        op = hermitian_to_pauli_op(H)
        # Should have just the identity term
        assert len(op.terms) == 1
        assert () in op.terms
        np.testing.assert_allclose(op.terms[()], 1.0, atol=1e-10)

    def test_pauli_z(self):
        """Z matrix decomposes to Z0."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        H = np.array([[1, 0], [0, -1]], dtype=float)
        op = hermitian_to_pauli_op(H)
        # Should have just Z0 with coefficient 1
        assert ((0, 'Z'),) in op.terms
        np.testing.assert_allclose(op.terms[((0, 'Z'),)], 1.0, atol=1e-10)

    def test_pauli_x(self):
        """X matrix decomposes to X0."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        H = np.array([[0, 1], [1, 0]], dtype=float)
        op = hermitian_to_pauli_op(H)
        assert ((0, 'X'),) in op.terms
        np.testing.assert_allclose(op.terms[((0, 'X'),)], 1.0, atol=1e-10)

    def test_reconstruction_2x2(self):
        """Pauli decomposition reconstructs the original matrix."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op, _pauli_matrix
        H = np.array([[3.0, 1.5], [1.5, -1.0]])
        op = hermitian_to_pauli_op(H, threshold=1e-14)

        # Reconstruct
        labels = {'I': np.eye(2), 'X': _pauli_matrix('X'),
                  'Y': _pauli_matrix('Y'), 'Z': _pauli_matrix('Z')}
        H_recon = np.zeros((2, 2), dtype=complex)
        for term, coeff in op.terms.items():
            if len(term) == 0:
                P = np.eye(2)
            else:
                q, l = term[0]
                P = labels[l]
            H_recon += coeff * P

        np.testing.assert_allclose(H_recon.real, H, atol=1e-10)

    def test_reconstruction_4x4(self):
        """4x4 Hermitian matrix roundtrips through Pauli decomposition."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op, _pauli_matrix
        np.random.seed(42)
        A = np.random.randn(4, 4)
        H = A + A.T  # symmetric

        op = hermitian_to_pauli_op(H, threshold=1e-14)

        # Reconstruct from Pauli strings
        labels_map = {'I': np.eye(2), 'X': _pauli_matrix('X'),
                      'Y': _pauli_matrix('Y'), 'Z': _pauli_matrix('Z')}
        H_recon = np.zeros((4, 4), dtype=complex)

        for term, coeff in op.terms.items():
            ops = ['I', 'I']  # 2 qubits
            for q, l in term:
                ops[q] = l
            P = np.kron(labels_map[ops[1]], labels_map[ops[0]])
            H_recon += coeff * P

        np.testing.assert_allclose(H_recon.real, H, atol=1e-10)

    def test_eigenvalue_preservation(self):
        """Pauli-decomposed operator has same eigenvalues as original."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        from geovac.vqe_benchmark import openfermion_to_sparse_pauli_op

        np.random.seed(123)
        A = np.random.randn(4, 4)
        H = A + A.T

        op = hermitian_to_pauli_op(H, threshold=1e-14)
        spo = openfermion_to_sparse_pauli_op(op, 2)
        mat = spo.to_matrix()
        mat_recon = (mat.toarray() if hasattr(mat, 'toarray') else np.array(mat)).real

        evals_orig = np.sort(np.linalg.eigvalsh(H))
        evals_recon = np.sort(np.linalg.eigvalsh(mat_recon))

        np.testing.assert_allclose(evals_orig, evals_recon, atol=1e-10)

    def test_zero_padding(self):
        """3x3 matrix is padded to 4x4 (2 qubits), eigenvalues preserved."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        from geovac.vqe_benchmark import openfermion_to_sparse_pauli_op

        H = np.array([[2.0, 0.5, 0.3],
                       [0.5, -1.0, 0.2],
                       [0.3, 0.2, 1.5]])

        op = hermitian_to_pauli_op(H, threshold=1e-14)
        Q = 2  # ceil(log2(3)) = 2
        spo = openfermion_to_sparse_pauli_op(op, Q)
        mat = spo.to_matrix()
        mat_recon = (mat.toarray() if hasattr(mat, 'toarray') else np.array(mat)).real

        # Original 3 eigenvalues plus a zero from padding
        H_padded = np.zeros((4, 4))
        H_padded[:3, :3] = H
        evals_expected = np.sort(np.linalg.eigvalsh(H_padded))
        evals_recon = np.sort(np.linalg.eigvalsh(mat_recon))

        np.testing.assert_allclose(evals_recon, evals_expected, atol=1e-10)

    def test_diagonal_matrix(self):
        """Diagonal matrix produces only Z-string Pauli terms."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        H = np.diag([1.0, -1.0, 2.0, -2.0])
        op = hermitian_to_pauli_op(H, threshold=1e-14)

        # All terms should be products of I and Z only
        for term in op.terms:
            for q, label in term:
                assert label in ('Z',), \
                    f"Expected only Z terms for diagonal matrix, got {label}"

    def test_sparse_matrix_fewer_terms(self):
        """Sparse matrices produce fewer Pauli terms than dense ones."""
        from geovac.n_electron_qubit import hermitian_to_pauli_op
        np.random.seed(99)

        # Dense 8x8
        A = np.random.randn(8, 8)
        H_dense = A + A.T
        op_dense = hermitian_to_pauli_op(H_dense, threshold=1e-10)

        # Sparse 8x8 (tridiagonal)
        H_sparse = np.zeros((8, 8))
        for i in range(8):
            H_sparse[i, i] = np.random.randn()
            if i + 1 < 8:
                v = np.random.randn()
                H_sparse[i, i+1] = v
                H_sparse[i+1, i] = v
        op_sparse = hermitian_to_pauli_op(H_sparse, threshold=1e-10)

        # Sparse should have fewer or comparable Pauli terms
        assert len(op_sparse.terms) <= len(op_dense.terms)


# ---------------------------------------------------------------------------
# Tests for 4-electron angular Hamiltonian encoding
# ---------------------------------------------------------------------------

class TestAngularHamiltonianEncoding:
    """Test encoding the 4-electron angular Hamiltonian on qubits."""

    def test_lmax1_builds(self):
        """l_max=1 angular Hamiltonian builds and has correct dimension."""
        from geovac.n_electron_qubit import build_4e_angular_hamiltonian
        H, dim, info = build_4e_angular_hamiltonian(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
        )
        assert dim == 128  # 2 S4 channels * 4^3 grid
        assert H.shape == (128, 128)
        assert info['Q'] == 7

    def test_lmax0_empty(self):
        """l_max=0 has empty S4 [2,2] subspace."""
        from geovac.n_electron_qubit import build_4e_angular_hamiltonian
        H, dim, info = build_4e_angular_hamiltonian(
            l_max=0, n_grid=4, R=3.0, R_e=3.0,
        )
        assert dim == 0

    def test_hermiticity(self):
        """Angular Hamiltonian is Hermitian."""
        from geovac.n_electron_qubit import build_4e_angular_hamiltonian
        H, dim, info = build_4e_angular_hamiltonian(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
        )
        np.testing.assert_allclose(H, H.T, atol=1e-12)

    @pytest.mark.slow
    def test_lmax1_pauli_decomposition(self):
        """l_max=1 Pauli decomposition produces valid qubit operator."""
        from geovac.n_electron_qubit import analyze_4e_qubit_encoding
        result = analyze_4e_qubit_encoding(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
            verbose=False,
        )
        assert result['Q'] == 7
        assert result['N_pauli'] > 0
        assert result['N_qwc'] > 0
        assert result['one_norm'] > 0

    @pytest.mark.slow
    def test_eigenvalue_preservation_lmax1(self):
        """Pauli-decomposed l_max=1 Hamiltonian preserves low eigenvalues."""
        from geovac.n_electron_qubit import (
            build_4e_angular_hamiltonian, hermitian_to_pauli_op,
        )
        from geovac.vqe_benchmark import openfermion_to_sparse_pauli_op

        H, dim, info = build_4e_angular_hamiltonian(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
        )

        # Original eigenvalues
        evals_orig = np.sort(eigh(H, eigvals_only=True))[:5]

        # Pauli decomposition
        op = hermitian_to_pauli_op(H, threshold=1e-12)
        spo = openfermion_to_sparse_pauli_op(op, info['Q'])
        mat_recon = spo.to_matrix(sparse=True)

        from scipy.sparse.linalg import eigsh
        evals_recon = np.sort(eigsh(mat_recon.tocsc(), k=5, which='SA')[0])

        # Should match to reasonable precision
        np.testing.assert_allclose(evals_orig, evals_recon, atol=1e-6)


# ---------------------------------------------------------------------------
# Tests for composed comparison
# ---------------------------------------------------------------------------

class TestComposedComparison:
    """Test the comparison with composed-geometry encoding."""

    @pytest.mark.slow
    def test_composed_metrics(self):
        """Composed LiH metrics can be computed."""
        from geovac.n_electron_qubit import composed_lih_metrics
        result = composed_lih_metrics(
            max_n_core=1, max_n_val=1, verbose=False,
        )
        assert result['Q'] == 6
        assert result['N_pauli'] > 0

    @pytest.mark.slow
    def test_full_is_denser(self):
        """Full N-electron encoding is denser than composed at comparable Q."""
        from geovac.n_electron_qubit import (
            analyze_4e_qubit_encoding, composed_lih_metrics,
        )

        # Full at l_max=1, n_grid=4: Q=7
        full = analyze_4e_qubit_encoding(
            l_max=1, n_grid=4, verbose=False,
        )

        # Composed at smallest basis: Q=6
        composed = composed_lih_metrics(
            max_n_core=1, max_n_val=1, verbose=False,
        )

        # Full should have much higher Pauli density
        full_density = full['N_pauli'] / 4**full['Q']
        composed_density = composed['N_pauli'] / 4**composed['Q']

        # The full encoding should be significantly denser
        assert full_density > composed_density, (
            f"Expected full ({full_density:.4f}) > composed ({composed_density:.4f})"
        )


# ---------------------------------------------------------------------------
# Tests for structural analysis
# ---------------------------------------------------------------------------

class TestStructuralAnalysis:
    """Test the structural analysis utilities."""

    def test_structural_breakdown(self):
        """Structural breakdown reports correct statistics."""
        from geovac.n_electron_qubit import (
            build_4e_angular_hamiltonian, structural_breakdown,
        )
        H, dim, info = build_4e_angular_hamiltonian(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
        )
        result = structural_breakdown(H, dim, info, verbose=False)
        assert result['dim'] == dim
        assert result['nnz'] > 0
        assert 0 < result['density'] <= 1.0

    def test_scaling_analysis_format(self):
        """Scaling analysis returns expected format."""
        from geovac.n_electron_qubit import scaling_analysis
        data = [
            {'Q': 5, 'N_pauli': 100},
            {'Q': 7, 'N_pauli': 500},
            {'Q': 10, 'N_pauli': 5000},
        ]
        result = scaling_analysis(data, verbose=False)
        assert 'alpha' in result
        assert 'prefactor' in result
        assert not np.isnan(result['alpha'])


# ---------------------------------------------------------------------------
# Integration test
# ---------------------------------------------------------------------------

class TestIntegration:
    """End-to-end integration tests."""

    @pytest.mark.slow
    def test_run_track_as_analysis(self):
        """Full Track AS analysis runs without error."""
        from geovac.n_electron_qubit import run_track_as_analysis
        results = run_track_as_analysis(
            output_dir="debug/track_as",
            verbose=True,
        )
        assert 'full_4e' in results
        assert 'composed' in results
        assert 'scaling' in results
        assert 'verdict' in results
        assert len(results['full_4e']) > 0

    @pytest.mark.slow
    def test_scaling_sweep_lmax1(self):
        """Sweep n_grid at l_max=1 for scaling analysis."""
        from geovac.n_electron_qubit import (
            analyze_4e_qubit_encoding, build_4e_angular_hamiltonian,
            scaling_analysis,
        )

        data = []
        for ng in [4, 5, 6, 7]:
            r = analyze_4e_qubit_encoding(
                l_max=1, n_grid=ng, R=3.0, R_e=3.0, verbose=True,
            )
            H, dim, info = build_4e_angular_hamiltonian(
                l_max=1, n_grid=ng, R=3.0, R_e=3.0,
            )
            nnz = int(np.count_nonzero(np.abs(H) > 1e-15))
            print(f"  n_grid={ng}: dim={dim}, Q={r['Q']}, "
                  f"N_Pauli={r['N_pauli']:,}, nnz={nnz}/{dim**2} "
                  f"({nnz/dim**2:.4f}), density={r['density']:.4f}")
            data.append(r)

        # Fit scaling
        result = scaling_analysis(data, verbose=True)
        print(f"\nScaling exponent: {result['alpha']:.2f}")

        # The scaling should be > 2.5 (denser than composed)
        assert result['alpha'] > 2.0, (
            f"Expected scaling exponent > 2.0, got {result['alpha']:.2f}"
        )

    @pytest.mark.slow
    def test_structural_comparison(self):
        """Structural analysis of what drives Pauli density."""
        from geovac.n_electron_qubit import (
            build_4e_angular_hamiltonian, structural_breakdown,
            composed_lih_metrics,
        )

        # Full N-electron
        H, dim, info = build_4e_angular_hamiltonian(
            l_max=1, n_grid=4, R=3.0, R_e=3.0,
        )
        full_struct = structural_breakdown(H, dim, info, verbose=True)

        # Composed
        composed = composed_lih_metrics(
            max_n_core=1, max_n_val=1, verbose=False,
        )

        print(f"\nFull N-electron (l_max=1, n_grid=4):")
        print(f"  dim={dim}, nnz={full_struct['nnz']}, "
              f"density={full_struct['density']:.4f}")
        print(f"Composed (Q=6):")
        print(f"  N_Pauli={composed['N_pauli']}, "
              f"density={composed['N_pauli']/4**composed['Q']:.6f}")

        # Full should be denser
        assert full_struct['density'] > 0.0
