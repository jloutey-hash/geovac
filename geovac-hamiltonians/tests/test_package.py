"""
Integration tests for geovac-hamiltonians package.

Validates:
  - Import works
  - Each system builds without error
  - n_terms matches published counts
  - Hermiticity check
  - Export methods work
"""

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Import test
# ---------------------------------------------------------------------------

def test_import():
    """Package imports successfully."""
    from geovac_hamiltonians import hamiltonian, GeoVacHamiltonian
    assert callable(hamiltonian)
    assert GeoVacHamiltonian is not None


def test_version():
    """Version string is present."""
    import geovac_hamiltonians
    assert hasattr(geovac_hamiltonians, '__version__')
    assert geovac_hamiltonians.__version__ == '0.1.0'


# ---------------------------------------------------------------------------
# H2 (fast, small)
# ---------------------------------------------------------------------------

class TestH2:
    """H2 bond-pair encoding."""

    @pytest.fixture(scope="class")
    def h2(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('H2')

    def test_n_qubits(self, h2):
        assert h2.n_qubits == 10

    def test_n_terms(self, h2):
        assert h2.n_terms == 112

    def test_one_norm_positive(self, h2):
        assert h2.one_norm > 0.0

    def test_hermiticity(self, h2):
        """All JW coefficients should be real."""
        of_op = h2.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12

    def test_metadata(self, h2):
        meta = h2.metadata
        assert meta['system'] == 'H2'

    def test_openfermion_export(self, h2):
        from openfermion import QubitOperator
        of_op = h2.to_openfermion()
        assert isinstance(of_op, QubitOperator)
        assert len(of_op.terms) == 112


# ---------------------------------------------------------------------------
# He (fast, small)
# ---------------------------------------------------------------------------

class TestHe:
    """He atomic encoding."""

    @pytest.fixture(scope="class")
    def he(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('He', max_n=2)

    def test_n_qubits(self, he):
        assert he.n_qubits > 0

    def test_n_terms_positive(self, he):
        assert he.n_terms > 0

    def test_one_norm_positive(self, he):
        assert he.one_norm > 0.0

    def test_hermiticity(self, he):
        of_op = he.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12


# ---------------------------------------------------------------------------
# LiH (composed, the headline result)
# ---------------------------------------------------------------------------

class TestLiH:
    """LiH composed encoding -- the headline benchmark."""

    @pytest.fixture(scope="class")
    def lih(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('LiH', R=3.015, l_max=2)

    def test_n_qubits(self, lih):
        assert lih.n_qubits == 30

    def test_n_terms_exact(self, lih):
        """LiH at l_max=2 must have exactly 334 Pauli terms (Paper 14)."""
        assert lih.n_terms == 334

    def test_one_norm(self, lih):
        assert abs(lih.one_norm - 37.33) < 0.5  # within 0.5 Ha

    def test_hermiticity(self, lih):
        of_op = lih.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12

    def test_metadata(self, lih):
        meta = lih.metadata
        assert meta['system'] == 'LiH'
        assert meta['R_bohr'] == 3.015
        assert meta['l_max'] == 2

    def test_openfermion_export(self, lih):
        from openfermion import QubitOperator
        of_op = lih.to_openfermion()
        assert isinstance(of_op, QubitOperator)
        assert len(of_op.terms) == 334


# ---------------------------------------------------------------------------
# BeH2 (composed)
# ---------------------------------------------------------------------------

class TestBeH2:
    """BeH2 composed encoding."""

    @pytest.fixture(scope="class")
    def beh2(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('BeH2', R=2.502, l_max=2)

    def test_n_qubits(self, beh2):
        assert beh2.n_qubits == 50

    def test_n_terms_exact(self, beh2):
        """BeH2 at l_max=2 must have exactly 556 Pauli terms (Paper 14)."""
        assert beh2.n_terms == 556

    def test_hermiticity(self, beh2):
        of_op = beh2.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12


# ---------------------------------------------------------------------------
# H2O (composed)
# ---------------------------------------------------------------------------

class TestH2O:
    """H2O composed encoding."""

    @pytest.fixture(scope="class")
    def h2o(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('H2O', R=1.809, l_max=2)

    def test_n_qubits(self, h2o):
        assert h2o.n_qubits == 70

    def test_n_terms_exact(self, h2o):
        """H2O at l_max=2 must have exactly 778 Pauli terms (Paper 14)."""
        assert h2o.n_terms == 778

    def test_hermiticity(self, h2o):
        of_op = h2o.to_openfermion()
        for coeff in of_op.terms.values():
            assert abs(coeff.imag) < 1e-12


# ---------------------------------------------------------------------------
# Qiskit export (optional)
# ---------------------------------------------------------------------------

class TestQiskitExport:
    """Test Qiskit export (skipped if qiskit not installed)."""

    @pytest.fixture(scope="class")
    def h2(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('H2')

    def test_qiskit_export(self, h2):
        try:
            spo = h2.to_qiskit()
        except ImportError:
            pytest.skip("qiskit not installed")
        assert spo.num_qubits == 10


# ---------------------------------------------------------------------------
# PennyLane export (optional)
# ---------------------------------------------------------------------------

class TestPennyLaneExport:
    """Test PennyLane export (skipped if pennylane not installed)."""

    @pytest.fixture(scope="class")
    def h2(self):
        from geovac_hamiltonians import hamiltonian
        return hamiltonian('H2')

    def test_pennylane_export(self, h2):
        try:
            pl_h = h2.to_pennylane()
        except ImportError:
            pytest.skip("pennylane not installed")
        assert len(pl_h.coeffs) == 112


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    """Test error cases."""

    def test_unknown_system(self):
        from geovac_hamiltonians import hamiltonian
        with pytest.raises(ValueError, match="Unknown system"):
            hamiltonian('Unobtanium')

    def test_case_insensitive(self):
        from geovac_hamiltonians import hamiltonian
        h1 = hamiltonian('h2')
        h2 = hamiltonian('H2')
        assert h1.n_terms == h2.n_terms

    def test_repr(self):
        from geovac_hamiltonians import hamiltonian
        h = hamiltonian('H2')
        r = repr(h)
        assert 'GeoVacHamiltonian' in r
        assert 'H2' in r
