"""
Tests for the general composed-geometry builder (Track BH).

Verifies that ``build_composed_hamiltonian`` with each spec factory
produces results IDENTICAL (to machine precision) to the original
hardcoded builders.
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    build_composed_lih,
    build_composed_beh2,
    build_composed_h2o,
    build_h2_bond_pair,
    lih_spec,
    beh2_spec,
    h2o_spec,
    h2_bond_pair_spec,
    he_spec,
)


# ---------------------------------------------------------------------------
# LiH: general builder vs hardcoded builder
# ---------------------------------------------------------------------------

class TestLiHGeneralBuilder:
    """Verify LiH general builder reproduces the hardcoded builder exactly."""

    @pytest.fixture(scope='class')
    def results(self):
        """Run both builders once and cache results."""
        old = build_composed_lih(
            max_n_core=2, max_n_val=2, verbose=False,
            pk_in_hamiltonian=True,
        )
        spec = lih_spec(max_n_core=2, max_n_val=2)
        new = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        return old, new

    def test_qubit_count(self, results):
        old, new = results
        assert old['Q'] == new['Q'] == 30

    def test_pauli_count(self, results):
        old, new = results
        assert old['N_pauli'] == new['N_pauli'] == 334

    def test_h1_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1'], new['h1'], atol=1e-14)

    def test_eri_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['eri'], new['eri'], atol=1e-14)

    def test_nuclear_repulsion_match(self, results):
        old, new = results
        assert abs(old['nuclear_repulsion'] - new['nuclear_repulsion']) < 1e-12

    def test_h1_pk_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1_pk'], new['h1_pk'], atol=1e-14)


# ---------------------------------------------------------------------------
# BeH2: general builder vs hardcoded builder
# ---------------------------------------------------------------------------

class TestBeH2GeneralBuilder:
    """Verify BeH2 general builder reproduces the hardcoded builder exactly."""

    @pytest.fixture(scope='class')
    def results(self):
        old = build_composed_beh2(
            max_n_core=2, max_n_val=2, verbose=False,
            pk_in_hamiltonian=True,
        )
        spec = beh2_spec(max_n_core=2, max_n_val=2)
        new = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        return old, new

    def test_qubit_count(self, results):
        old, new = results
        assert old['Q'] == new['Q'] == 50

    def test_pauli_count(self, results):
        old, new = results
        assert old['N_pauli'] == new['N_pauli'] == 556

    def test_h1_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1'], new['h1'], atol=1e-14)

    def test_eri_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['eri'], new['eri'], atol=1e-14)

    def test_nuclear_repulsion_match(self, results):
        old, new = results
        assert abs(old['nuclear_repulsion'] - new['nuclear_repulsion']) < 1e-12

    def test_h1_pk_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1_pk'], new['h1_pk'], atol=1e-14)


# ---------------------------------------------------------------------------
# H2O: general builder vs hardcoded builder
# ---------------------------------------------------------------------------

class TestH2OGeneralBuilder:
    """Verify H2O general builder reproduces the hardcoded builder exactly."""

    @pytest.fixture(scope='class')
    def results(self):
        old = build_composed_h2o(
            max_n_core=2, max_n_val=2, verbose=False,
            pk_in_hamiltonian=True,
        )
        spec = h2o_spec(max_n_core=2, max_n_val=2)
        new = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        return old, new

    def test_qubit_count(self, results):
        old, new = results
        assert old['Q'] == new['Q'] == 70

    def test_pauli_count(self, results):
        old, new = results
        assert old['N_pauli'] == new['N_pauli'] == 778

    def test_h1_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1'], new['h1'], atol=1e-14)

    def test_eri_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['eri'], new['eri'], atol=1e-14)

    def test_nuclear_repulsion_match(self, results):
        old, new = results
        assert abs(old['nuclear_repulsion'] - new['nuclear_repulsion']) < 1e-12

    def test_h1_pk_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1_pk'], new['h1_pk'], atol=1e-14)


# ---------------------------------------------------------------------------
# H2 bond pair: general builder vs hardcoded builder
# ---------------------------------------------------------------------------

class TestH2BondPairGeneralBuilder:
    """Verify H2 bond pair general builder reproduces the hardcoded builder."""

    @pytest.fixture(scope='class')
    def results(self):
        old = build_h2_bond_pair(max_n=2, R=1.4, verbose=False)
        spec = h2_bond_pair_spec(max_n=2, R=1.4)
        new = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        return old, new

    def test_qubit_count(self, results):
        old, new = results
        assert old['Q'] == new['Q'] == 10

    def test_pauli_count(self, results):
        old, new = results
        assert old['N_pauli'] == new['N_pauli'] == 112

    def test_h1_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['h1'], new['h1'], atol=1e-14)

    def test_eri_match(self, results):
        old, new = results
        np.testing.assert_allclose(old['eri'], new['eri'], atol=1e-14)

    def test_nuclear_repulsion_match(self, results):
        old, new = results
        assert abs(old['nuclear_repulsion'] - new['nuclear_repulsion']) < 1e-12


# ---------------------------------------------------------------------------
# He: standalone test via general builder
# ---------------------------------------------------------------------------

class TestHeGeneralBuilder:
    """Verify He general builder produces expected Pauli count."""

    @pytest.fixture(scope='class')
    def result(self):
        spec = he_spec(max_n=2)
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)

    def test_qubit_count(self, result):
        assert result['Q'] == 10

    def test_pauli_count(self, result):
        # He via the composed pipeline (single block Z=2, no nuclear repulsion)
        # gives 112 Pauli terms -- different from VQE benchmark's 120 which
        # uses a different code path (build_geovac_he in vqe_benchmark.py).
        assert result['N_pauli'] == 112


# ---------------------------------------------------------------------------
# MolecularSpec dataclass tests
# ---------------------------------------------------------------------------

class TestMolecularSpec:
    """Basic tests for the spec factory functions."""

    def test_lih_spec_blocks(self):
        spec = lih_spec()
        assert spec.name == 'LiH'
        assert len(spec.blocks) == 2
        assert spec.blocks[0].block_type == 'core'
        assert spec.blocks[1].block_type == 'bond'
        assert spec.blocks[1].has_h_partner is True

    def test_beh2_spec_blocks(self):
        spec = beh2_spec()
        assert spec.name == 'BeH2'
        assert len(spec.blocks) == 3
        assert spec.blocks[0].block_type == 'core'
        # Two independent bond blocks
        assert spec.blocks[1].block_type == 'bond'
        assert spec.blocks[2].block_type == 'bond'

    def test_h2o_spec_blocks(self):
        spec = h2o_spec()
        assert spec.name == 'H2O'
        assert len(spec.blocks) == 5
        # core + 2 bonds + 2 lone pairs
        types = [b.block_type for b in spec.blocks]
        assert types.count('core') == 1
        assert types.count('bond') == 2
        assert types.count('lone_pair') == 2

    def test_h2_bond_pair_spec_blocks(self):
        spec = h2_bond_pair_spec()
        assert spec.name == 'H2'
        assert len(spec.blocks) == 1
        assert spec.blocks[0].block_type == 'bond_pair'
        assert spec.blocks[0].pk_A == 0.0

    def test_he_spec_blocks(self):
        spec = he_spec()
        assert spec.name == 'He'
        assert len(spec.blocks) == 1
        assert spec.blocks[0].Z_center == 2.0

    def test_lih_pk_override(self):
        spec = lih_spec(A_pk=10.0, B_pk=11.0)
        assert spec.blocks[1].pk_A == 10.0
        assert spec.blocks[1].pk_B == 11.0

    def test_lih_no_pk(self):
        spec = lih_spec(include_pk=False)
        assert spec.blocks[1].pk_A == 0.0
        assert spec.blocks[1].pk_B == 0.0


# ---------------------------------------------------------------------------
# pk_in_hamiltonian=False test
# ---------------------------------------------------------------------------

class TestPKSeparation:
    """Verify pk_in_hamiltonian=False keeps PK out of h1."""

    def test_lih_pk_not_in_h1(self):
        spec = lih_spec(include_pk=True)
        result_with = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        result_without = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)

        # h1 should differ (PK not added)
        assert not np.allclose(result_with['h1'], result_without['h1'])

        # h1_pk should be the same
        np.testing.assert_allclose(result_with['h1_pk'], result_without['h1_pk'], atol=1e-14)

        # h1_with = h1_without + h1_pk
        np.testing.assert_allclose(
            result_with['h1'],
            result_without['h1'] + result_without['h1_pk'],
            atol=1e-14,
        )
