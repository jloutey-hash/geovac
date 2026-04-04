"""
Tests for H2 bond-pair qubit encoding (Track AZ)
=================================================

Validates:
  - build_h2_bond_pair() produces correct qubit Hamiltonians
  - Hermiticity of the qubit Hamiltonian
  - Pauli term counts at max_n=2,3
  - R-independence of Pauli term count (selection rule structure)
  - 1-norm is positive and finite
  - Qubit count matches 2*M
  - ecosystem_export integration

Author: GeoVac Development Team
Date: April 2026
"""

import numpy as np
import pytest

from openfermion import QubitOperator


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def h2_maxn2():
    """H2 bond-pair at max_n=2, R=1.4."""
    from geovac.composed_qubit import build_h2_bond_pair
    return build_h2_bond_pair(max_n=2, R=1.4, verbose=False)


@pytest.fixture(scope="module")
def h2_maxn3():
    """H2 bond-pair at max_n=3, R=1.4."""
    from geovac.composed_qubit import build_h2_bond_pair
    return build_h2_bond_pair(max_n=3, R=1.4, verbose=False)


# ---------------------------------------------------------------------------
# Basic structure tests
# ---------------------------------------------------------------------------

class TestH2BondPairStructure:
    """Test basic structural properties of the H2 bond-pair encoding."""

    def test_qubit_count_maxn2(self, h2_maxn2: dict) -> None:
        """Q = 2*M = 2*5 = 10 at max_n=2."""
        assert h2_maxn2['Q'] == 2 * h2_maxn2['M']
        assert h2_maxn2['Q'] == 10

    def test_qubit_count_maxn3(self, h2_maxn3: dict) -> None:
        """Q = 2*M = 2*14 = 28 at max_n=3."""
        assert h2_maxn3['Q'] == 2 * h2_maxn3['M']
        assert h2_maxn3['Q'] == 28

    def test_spatial_orbitals_maxn2(self, h2_maxn2: dict) -> None:
        """M = sum_{n=1}^{2} n^2 = 1+4 = 5."""
        assert h2_maxn2['M'] == 5

    def test_spatial_orbitals_maxn3(self, h2_maxn3: dict) -> None:
        """M = sum_{n=1}^{3} n^2 = 1+4+9 = 14."""
        assert h2_maxn3['M'] == 14

    def test_pauli_count_maxn2(self, h2_maxn2: dict) -> None:
        """Pauli term count at max_n=2 should be exactly 112."""
        assert h2_maxn2['N_pauli'] == 112

    def test_pauli_count_maxn3(self, h2_maxn3: dict) -> None:
        """Pauli term count at max_n=3 should be exactly 2627."""
        assert h2_maxn3['N_pauli'] == 2627

    def test_states_list(self, h2_maxn2: dict) -> None:
        """States should be (n,l,m) tuples in canonical order."""
        states = h2_maxn2['states']
        assert len(states) == 5
        assert states[0] == (1, 0, 0)
        assert states[1] == (2, 0, 0)

    def test_h1_diagonal(self, h2_maxn2: dict) -> None:
        """h1 should be diagonal with -1/(2n^2) entries."""
        h1 = h2_maxn2['h1']
        # Off-diagonal should be zero
        assert np.allclose(h1 - np.diag(np.diag(h1)), 0.0)
        # 1s orbital: -1/2
        assert abs(h1[0, 0] - (-0.5)) < 1e-12
        # 2s orbital: -1/8
        assert abs(h1[1, 1] - (-0.125)) < 1e-12

    def test_nuclear_repulsion(self, h2_maxn2: dict) -> None:
        """V_NN = 1/R = 1/1.4."""
        assert abs(h2_maxn2['nuclear_repulsion'] - 1.0 / 1.4) < 1e-12


# ---------------------------------------------------------------------------
# Hermiticity
# ---------------------------------------------------------------------------

class TestH2BondPairHermiticity:
    """Test that the qubit Hamiltonian is Hermitian."""

    def test_real_coefficients(self, h2_maxn2: dict) -> None:
        """All Pauli coefficients should be real (JW of real integrals)."""
        qop = h2_maxn2['qubit_op']
        for coeff in qop.terms.values():
            assert abs(coeff.imag) < 1e-12, f"Non-real coefficient: {coeff}"

    def test_matrix_hermiticity(self, h2_maxn2: dict) -> None:
        """The Hamiltonian matrix should be Hermitian."""
        from openfermion import get_sparse_operator
        qop = h2_maxn2['qubit_op']
        mat = get_sparse_operator(qop, n_qubits=h2_maxn2['Q']).toarray()
        assert np.allclose(mat, mat.conj().T, atol=1e-12)

    def test_eri_symmetry(self, h2_maxn2: dict) -> None:
        """ERI tensor should have chemist symmetry: (pq|rs) = (rs|pq)."""
        eri = h2_maxn2['eri']
        assert np.allclose(eri, eri.transpose(2, 3, 0, 1), atol=1e-12)


# ---------------------------------------------------------------------------
# R-independence of Pauli count
# ---------------------------------------------------------------------------

class TestR_Independence:
    """Pauli term count should be R-independent (selection rule structure)."""

    def test_pauli_count_r_independent(self) -> None:
        """Pauli count at max_n=2 should be the same at multiple R values."""
        from geovac.composed_qubit import build_h2_bond_pair

        R_values = [0.5, 1.0, 1.4, 2.0, 3.0]
        counts = []
        for R in R_values:
            result = build_h2_bond_pair(max_n=2, R=R, verbose=False)
            counts.append(result['N_pauli'])

        # All should be equal
        assert all(c == counts[0] for c in counts), (
            f"Pauli counts vary with R: {dict(zip(R_values, counts))}"
        )

    def test_one_norm_varies_with_r(self) -> None:
        """1-norm should vary with R (V_NN = 1/R changes coefficients)."""
        from geovac.composed_qubit import build_h2_bond_pair
        from geovac.trotter_bounds import pauli_1norm

        r1 = build_h2_bond_pair(max_n=2, R=0.5, verbose=False)
        r2 = build_h2_bond_pair(max_n=2, R=3.0, verbose=False)

        norm1 = pauli_1norm(r1['qubit_op'])
        norm2 = pauli_1norm(r2['qubit_op'])

        assert norm1 != norm2, "1-norm should differ at different R"
        assert norm1 > norm2, "1-norm should be larger at smaller R (larger V_NN)"


# ---------------------------------------------------------------------------
# 1-norm tests
# ---------------------------------------------------------------------------

class TestOneNorm:
    """Test 1-norm properties."""

    def test_positive_and_finite(self, h2_maxn2: dict) -> None:
        from geovac.trotter_bounds import pauli_1norm
        norm = pauli_1norm(h2_maxn2['qubit_op'])
        assert norm > 0.0
        assert np.isfinite(norm)

    def test_one_norm_maxn2(self, h2_maxn2: dict) -> None:
        """1-norm at max_n=2, R=1.4 should be approximately 8.17."""
        from geovac.trotter_bounds import pauli_1norm
        norm = pauli_1norm(h2_maxn2['qubit_op'])
        assert abs(norm - 8.1739) < 0.01


# ---------------------------------------------------------------------------
# Ecosystem export integration
# ---------------------------------------------------------------------------

class TestEcosystemExportH2:
    """Test that ecosystem_export.hamiltonian('H2') uses bond-pair encoding."""

    def test_bond_pair_default(self) -> None:
        """Default H2 should use bond-pair encoding, not STO-3G."""
        from geovac.ecosystem_export import hamiltonian
        h = hamiltonian('H2', max_n=2, verbose=False)
        meta = h.metadata
        assert meta.get('encoding') == 'bond-pair'
        assert meta['Q'] == 10
        assert h.n_terms == 112

    def test_sto3g_fallback(self) -> None:
        """basis='sto-3g' should fall back to Gaussian reference."""
        from geovac.ecosystem_export import hamiltonian
        # The _build_h2 function accepts basis via the internal call
        from geovac.ecosystem_export import _build_h2
        h = _build_h2(basis='sto-3g', verbose=False)
        meta = h.metadata
        assert meta.get('basis') == 'STO-3G'
        assert meta['Q'] == 4
        assert h.n_terms == 15

    def test_r_parameter_passed(self) -> None:
        """R parameter should be passed through to bond-pair builder."""
        from geovac.ecosystem_export import hamiltonian
        h = hamiltonian('H2', max_n=2, R=2.0, verbose=False)
        meta = h.metadata
        assert meta['R_bohr'] == 2.0
