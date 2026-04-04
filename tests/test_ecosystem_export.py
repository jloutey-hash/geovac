"""
Tests for ecosystem export — GeoVac -> OpenFermion / Qiskit / PennyLane
======================================================================

Validates:
  - GeoVacHamiltonian properties (n_qubits, n_terms, one_norm)
  - Export round-trips (eigenvalue agreement)
  - Hermiticity of exported operators
  - Published term counts (LiH l_max=2: 334 Pauli terms)
  - 1-norm consistency
  - Multiple systems (He, H2, LiH)

Author: GeoVac Development Team
Date: April 2026
"""

import warnings

import numpy as np
import pytest

from openfermion import QubitOperator

from geovac.ecosystem_export import (
    GeoVacHamiltonian,
    hamiltonian,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def he_hamiltonian() -> GeoVacHamiltonian:
    """He system (small, fast)."""
    return hamiltonian('He', max_n=2, verbose=False)


@pytest.fixture(scope="module")
def h2_hamiltonian() -> GeoVacHamiltonian:
    """H2 system (STO-3G, fast)."""
    return hamiltonian('H2', verbose=False)


@pytest.fixture(scope="module")
def lih_hamiltonian() -> GeoVacHamiltonian:
    """LiH composed system at R=3.015, l_max=2."""
    return hamiltonian('LiH', R=3.015, l_max=2, verbose=False)


# ---------------------------------------------------------------------------
# Basic property tests
# ---------------------------------------------------------------------------

class TestGeoVacHamiltonianProperties:
    """Test n_qubits, n_terms, one_norm properties."""

    def test_he_n_qubits(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        assert he_hamiltonian.n_qubits > 0

    def test_he_n_terms(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        assert he_hamiltonian.n_terms > 0

    def test_he_one_norm_positive(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        assert he_hamiltonian.one_norm > 0.0

    def test_h2_n_qubits(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        # H2 bond-pair (max_n=2): 5 spatial orbitals -> 10 spin-orbitals = 10 qubits
        assert h2_hamiltonian.n_qubits == 10

    def test_h2_n_terms(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        # H2 bond-pair (max_n=2): exactly 112 Pauli terms
        assert h2_hamiltonian.n_terms == 112

    def test_lih_n_terms_published(self, lih_hamiltonian: GeoVacHamiltonian) -> None:
        """LiH at l_max=2 should have exactly 334 Pauli terms (Paper 14)."""
        assert lih_hamiltonian.n_terms == 334

    def test_lih_n_qubits(self, lih_hamiltonian: GeoVacHamiltonian) -> None:
        # LiH composed l_max=2: 15 spatial orbitals -> 30 qubits
        assert lih_hamiltonian.n_qubits == 30

    def test_one_norm_consistency(self, lih_hamiltonian: GeoVacHamiltonian) -> None:
        """one_norm should match manual computation from QubitOperator."""
        of_op = lih_hamiltonian.to_openfermion()
        manual_norm = sum(abs(c) for c in of_op.terms.values())
        assert abs(lih_hamiltonian.one_norm - manual_norm) < 1e-12

    def test_repr(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        r = repr(he_hamiltonian)
        assert 'GeoVacHamiltonian' in r
        assert 'He' in r

    def test_metadata(self, lih_hamiltonian: GeoVacHamiltonian) -> None:
        meta = lih_hamiltonian.metadata
        assert meta['system'] == 'LiH'
        assert meta['R_bohr'] == 3.015


# ---------------------------------------------------------------------------
# OpenFermion export
# ---------------------------------------------------------------------------

class TestOpenFermionExport:
    """Test to_openfermion()."""

    def test_returns_qubit_operator(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        of_op = he_hamiltonian.to_openfermion()
        assert isinstance(of_op, QubitOperator)

    def test_term_count_matches(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        of_op = he_hamiltonian.to_openfermion()
        assert len(of_op.terms) == he_hamiltonian.n_terms

    def test_hermiticity(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """QubitOperator should be Hermitian (all real coefficients for JW)."""
        of_op = h2_hamiltonian.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12, f"Non-real coefficient: {coeff}"


# ---------------------------------------------------------------------------
# Qiskit export
# ---------------------------------------------------------------------------

class TestQiskitExport:
    """Test to_qiskit()."""

    def test_returns_sparse_pauli_op(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        from qiskit.quantum_info import SparsePauliOp
        spo = he_hamiltonian.to_qiskit()
        assert isinstance(spo, SparsePauliOp)

    def test_num_qubits_matches(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        spo = h2_hamiltonian.to_qiskit()
        assert spo.num_qubits == h2_hamiltonian.n_qubits

    def test_hermiticity(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """SparsePauliOp matrix should be Hermitian."""
        spo = h2_hamiltonian.to_qiskit()
        mat = spo.to_matrix()
        if hasattr(mat, 'toarray'):
            mat = mat.toarray()
        assert np.allclose(mat, mat.conj().T, atol=1e-12)

    def test_eigenvalue_roundtrip_h2(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """Eigenvalues from Qiskit export should match OpenFermion."""
        from openfermion import get_sparse_operator

        of_op = h2_hamiltonian.to_openfermion()
        spo = h2_hamiltonian.to_qiskit()

        # OpenFermion eigenvalues
        of_mat = get_sparse_operator(of_op, n_qubits=h2_hamiltonian.n_qubits)
        of_eigvals = sorted(np.real(np.linalg.eigvalsh(of_mat.toarray())))

        # Qiskit eigenvalues
        qk_mat = spo.to_matrix()
        if hasattr(qk_mat, 'toarray'):
            qk_mat = qk_mat.toarray()
        qk_eigvals = sorted(np.real(np.linalg.eigvalsh(qk_mat)))

        np.testing.assert_allclose(of_eigvals, qk_eigvals, atol=1e-10)

    def test_eigenvalue_roundtrip_he(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        """He eigenvalue round-trip (slightly larger system)."""
        from openfermion import get_sparse_operator

        of_op = he_hamiltonian.to_openfermion()
        spo = he_hamiltonian.to_qiskit()
        nq = he_hamiltonian.n_qubits

        of_mat = get_sparse_operator(of_op, n_qubits=nq).toarray()
        of_gs = np.real(np.linalg.eigvalsh(of_mat))[0]

        qk_mat = spo.to_matrix()
        if hasattr(qk_mat, 'toarray'):
            qk_mat = qk_mat.toarray()
        qk_gs = np.real(np.linalg.eigvalsh(qk_mat))[0]

        assert abs(of_gs - qk_gs) < 1e-10


# ---------------------------------------------------------------------------
# PennyLane export
# ---------------------------------------------------------------------------

class TestPennyLaneExport:
    """Test to_pennylane()."""

    def test_returns_hamiltonian(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        import pennylane as qml
        pl_h = he_hamiltonian.to_pennylane()
        assert isinstance(pl_h, qml.Hamiltonian)

    def test_coefficient_count_matches(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        pl_h = h2_hamiltonian.to_pennylane()
        assert len(pl_h.coeffs) == h2_hamiltonian.n_terms

    def test_hermiticity_h2(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """PennyLane Hamiltonian matrix should be Hermitian."""
        pl_h = h2_hamiltonian.to_pennylane()
        mat = _pennylane_to_matrix(pl_h, h2_hamiltonian.n_qubits)
        assert np.allclose(mat, mat.conj().T, atol=1e-12)

    def test_eigenvalue_roundtrip_h2(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """PennyLane eigenvalues should match OpenFermion for H2."""
        from openfermion import get_sparse_operator

        of_op = h2_hamiltonian.to_openfermion()
        pl_h = h2_hamiltonian.to_pennylane()
        nq = h2_hamiltonian.n_qubits

        of_mat = get_sparse_operator(of_op, n_qubits=nq).toarray()
        of_gs = sorted(np.real(np.linalg.eigvalsh(of_mat)))[0]

        pl_mat = _pennylane_to_matrix(pl_h, nq)
        pl_gs = sorted(np.real(np.linalg.eigvalsh(pl_mat)))[0]

        assert abs(of_gs - pl_gs) < 1e-10

    def test_eigenvalue_roundtrip_he(self, he_hamiltonian: GeoVacHamiltonian) -> None:
        """PennyLane eigenvalues should match OpenFermion for He."""
        from openfermion import get_sparse_operator

        of_op = he_hamiltonian.to_openfermion()
        pl_h = he_hamiltonian.to_pennylane()
        nq = he_hamiltonian.n_qubits

        of_mat = get_sparse_operator(of_op, n_qubits=nq).toarray()
        of_gs = sorted(np.real(np.linalg.eigvalsh(of_mat)))[0]

        pl_mat = _pennylane_to_matrix(pl_h, nq)
        pl_gs = sorted(np.real(np.linalg.eigvalsh(pl_mat)))[0]

        assert abs(of_gs - pl_gs) < 1e-10


def _pennylane_to_matrix(pl_hamiltonian: "Any", n_qubits: int) -> np.ndarray:
    """Convert a PennyLane Hamiltonian to a dense matrix for testing."""
    import pennylane as qml
    mat = qml.matrix(pl_hamiltonian, wire_order=list(range(n_qubits)))
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    return np.array(mat)


# ---------------------------------------------------------------------------
# Convenience entry point tests
# ---------------------------------------------------------------------------

class TestHamiltonianEntryPoint:
    """Test the hamiltonian() convenience function."""

    def test_case_insensitive(self) -> None:
        h1 = hamiltonian('he', max_n=2, verbose=False)
        h2 = hamiltonian('He', max_n=2, verbose=False)
        h3 = hamiltonian('HE', max_n=2, verbose=False)
        assert h1.n_terms == h2.n_terms == h3.n_terms

    def test_unknown_system_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown system"):
            hamiltonian('Unobtanium')

    def test_h2_builds(self) -> None:
        h = hamiltonian('H2', verbose=False)
        # Default is now bond-pair encoding (max_n=2): 112 Pauli terms, 10 qubits
        assert h.n_terms == 112
        assert h.n_qubits == 10

    def test_lih_builds(self) -> None:
        h = hamiltonian('LiH', R=3.015, l_max=2, verbose=False)
        assert h.n_terms == 334

    def test_all_three_exports_lih(self, lih_hamiltonian: GeoVacHamiltonian) -> None:
        """LiH should export successfully to all three formats."""
        of_op = lih_hamiltonian.to_openfermion()
        assert len(of_op.terms) == 334

        spo = lih_hamiltonian.to_qiskit()
        assert spo.num_qubits == 30

        pl_h = lih_hamiltonian.to_pennylane()
        assert len(pl_h.coeffs) == 334


# ---------------------------------------------------------------------------
# Cross-format eigenvalue consistency (small system only due to 2^Q cost)
# ---------------------------------------------------------------------------

class TestCrossFormatConsistency:
    """Verify all three exports produce the same ground-state energy."""

    def test_h2_ground_state_all_formats(self, h2_hamiltonian: GeoVacHamiltonian) -> None:
        """H2 ground state energy should match across all three export formats."""
        from openfermion import get_sparse_operator

        nq = h2_hamiltonian.n_qubits

        # OpenFermion
        of_op = h2_hamiltonian.to_openfermion()
        of_mat = get_sparse_operator(of_op, n_qubits=nq).toarray()
        of_gs = np.real(np.linalg.eigvalsh(of_mat))[0]

        # Qiskit
        spo = h2_hamiltonian.to_qiskit()
        qk_mat = spo.to_matrix()
        if hasattr(qk_mat, 'toarray'):
            qk_mat = qk_mat.toarray()
        qk_gs = np.real(np.linalg.eigvalsh(qk_mat))[0]

        # PennyLane
        pl_h = h2_hamiltonian.to_pennylane()
        pl_mat = _pennylane_to_matrix(pl_h, nq)
        pl_gs = np.real(np.linalg.eigvalsh(pl_mat))[0]

        # All should agree to ~1e-10
        assert abs(of_gs - qk_gs) < 1e-10
        assert abs(of_gs - pl_gs) < 1e-10
        assert abs(qk_gs - pl_gs) < 1e-10
