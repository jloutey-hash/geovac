"""
Tests for coupled composition scoping (Track CB).

Validates:
1. Gaunt selection rules preserved in cross-block ERIs
2. PK removal is complete (no PK terms in coupled h1)
3. Coupled Hamiltonian is Hermitian
4. Within-block ERIs unchanged by cross-block addition
5. Pauli count and 1-norm regression values
6. FCI solver consistency
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    lih_spec,
    _enumerate_states,
    _ck_coefficient,
)
from geovac.coupled_composition import (
    build_coupled_hamiltonian,
    coupled_fci_energy,
)


@pytest.fixture
def lih_data():
    """Build LiH composed and coupled Hamiltonians for comparison."""
    spec = lih_spec(max_n_core=2, max_n_val=2, include_pk=True)
    composed = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
    coupled = build_coupled_hamiltonian(spec, verbose=False)
    no_pk = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)
    return {'spec': spec, 'composed': composed, 'coupled': coupled, 'no_pk': no_pk}


class TestGauntSelectionRules:
    """Verify Gaunt selection rules in cross-block ERIs."""

    def test_cross_block_eris_obey_gaunt(self, lih_data):
        """Cross-block ERIs should satisfy m-conservation from Gaunt coupling.

        In chemist notation (pq|rs), the physicist ERI is <pr|qs>.
        m-conservation for physicist <ab|cd>: ma + mb = mc + md.
        Mapping to chemist: a=p, b=r, c=q, d=s → mp + mr = mq + ms.
        """
        coupled = lih_data['coupled']
        no_pk = lih_data['no_pk']
        eri_cross = coupled['eri'] - no_pk['eri']

        # Enumerate states for each sub-block
        states = _enumerate_states(2)

        M = coupled['M']
        n_violations = 0
        for p in range(M):
            for q in range(M):
                for r in range(M):
                    for s in range(M):
                        if abs(eri_cross[p, q, r, s]) < 1e-14:
                            continue
                        mp = states[p % 5][2]
                        mq = states[q % 5][2]
                        mr = states[r % 5][2]
                        ms = states[s % 5][2]

                        # m-conservation: mp + mr = mq + ms
                        # (from physicist <pr|qs> → ma + mb = mc + md)
                        if mp + mr != mq + ms:
                            n_violations += 1

        assert n_violations == 0, f"Found {n_violations} m-conservation violations"

    def test_cross_block_eri_count(self, lih_data):
        """Cross-block ERI count should match BX-3b result (2/3 ratio)."""
        coupled = lih_data['coupled']
        assert coupled['cross_block_eri_count'] == 130
        assert coupled['within_block_eri_count'] == 195
        ratio = coupled['cross_block_eri_count'] / coupled['within_block_eri_count']
        assert abs(ratio - 2.0 / 3.0) < 1e-10


class TestPKRemoval:
    """Verify PK is fully removed from coupled Hamiltonian."""

    def test_no_pk_in_coupled_h1(self, lih_data):
        """Coupled h1 should not contain PK matrix elements."""
        coupled = lih_data['coupled']
        no_pk = lih_data['no_pk']

        # h1 should be identical between coupled and no-PK
        np.testing.assert_allclose(
            coupled['h1'], no_pk['h1'], atol=1e-14,
            err_msg="Coupled h1 differs from no-PK h1 (PK not fully removed)"
        )

    def test_pk_matrix_returned(self, lih_data):
        """PK matrix should still be returned for reference."""
        coupled = lih_data['coupled']
        assert coupled['h1_pk'] is not None
        # PK should have nonzero entries
        assert np.count_nonzero(np.abs(coupled['h1_pk']) > 1e-15) > 0


class TestHermiticity:
    """Verify coupled Hamiltonian is Hermitian."""

    def test_eri_symmetry(self, lih_data):
        """ERI tensor should satisfy (pq|rs) = (rs|pq)."""
        eri = lih_data['coupled']['eri']
        np.testing.assert_allclose(
            eri, eri.transpose(2, 3, 0, 1), atol=1e-14,
            err_msg="ERI tensor not symmetric under (pq|rs) = (rs|pq)"
        )

    def test_h1_symmetric(self, lih_data):
        """h1 should be symmetric (it's diagonal, so trivially true)."""
        h1 = lih_data['coupled']['h1']
        np.testing.assert_allclose(
            h1, h1.T, atol=1e-14,
            err_msg="h1 matrix not symmetric"
        )


class TestWithinBlockUnchanged:
    """Verify within-block ERIs are not modified by cross-block addition."""

    def test_within_block_eris_preserved(self, lih_data):
        """Within-block ERIs should be identical in composed and coupled."""
        composed = lih_data['composed']
        coupled = lih_data['coupled']

        # Core block: indices 0-4
        np.testing.assert_allclose(
            composed['eri'][:5, :5, :5, :5],
            coupled['eri'][:5, :5, :5, :5],
            atol=1e-14,
            err_msg="Core block ERIs modified"
        )

        # Bond center block: indices 5-9
        np.testing.assert_allclose(
            composed['eri'][5:10, 5:10, 5:10, 5:10],
            coupled['eri'][5:10, 5:10, 5:10, 5:10],
            atol=1e-14,
            err_msg="Bond center block ERIs modified"
        )

        # Bond partner block: indices 10-14
        np.testing.assert_allclose(
            composed['eri'][10:15, 10:15, 10:15, 10:15],
            coupled['eri'][10:15, 10:15, 10:15, 10:15],
            atol=1e-14,
            err_msg="Bond partner block ERIs modified"
        )


class TestPauliRegression:
    """Regression tests for Pauli count and 1-norm."""

    def test_composed_pauli_count(self, lih_data):
        """Composed LiH at Q=30 should have exactly 334 Pauli terms."""
        assert lih_data['composed']['N_pauli'] == 334

    def test_coupled_pauli_count(self, lih_data):
        """Coupled LiH at Q=30 should have exactly 854 Pauli terms."""
        assert lih_data['coupled']['N_pauli'] == 854

    def test_pauli_ratio(self, lih_data):
        """Pauli ratio should be ~2.56x."""
        ratio = lih_data['coupled']['N_pauli'] / lih_data['composed']['N_pauli']
        assert abs(ratio - 2.56) < 0.1

    def test_coupled_one_norm(self, lih_data):
        """Coupled 1-norm should be ~85.69 Ha."""
        one_norm = sum(abs(c) for c in lih_data['coupled']['qubit_op'].terms.values())
        assert abs(one_norm - 85.69) < 1.0


class TestFCISolver:
    """Test the sector-restricted FCI solver."""

    def test_fci_hermiticity(self, lih_data):
        """FCI for a simple case should give real eigenvalues."""
        # Restrict to valence subspace for fast test
        h1 = lih_data['coupled']['h1'][5:15, 5:15]
        eri = lih_data['coupled']['eri'][5:15, 5:15, 5:15, 5:15]
        result = {
            'M': 10, 'Q': 20,
            'h1': h1, 'eri': eri,
            'nuclear_repulsion': 0.0,
        }
        fci = coupled_fci_energy(result, 2, verbose=False)
        assert fci['exact_diag'] is True
        # Eigenvalues should be real (check that they are reasonable)
        assert all(isinstance(e, float) for e in fci['eigenvalues'])

    def test_fci_valence_ground_state(self, lih_data):
        """2-electron valence FCI should give a reasonable energy."""
        spec = lih_data['spec']
        nr = spec.nuclear_repulsion_constant

        h1 = lih_data['no_pk']['h1'][5:15, 5:15]
        eri = lih_data['no_pk']['eri'][5:15, 5:15, 5:15, 5:15]
        result = {
            'M': 10, 'Q': 20,
            'h1': h1, 'eri': eri,
            'nuclear_repulsion': nr,
        }
        fci = coupled_fci_energy(result, 2, verbose=False)
        # Should be in the right ballpark for LiH
        assert -10.0 < fci['E_coupled'] < -5.0
