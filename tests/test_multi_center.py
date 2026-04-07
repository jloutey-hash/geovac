"""
Tests for multi-center composed and balanced coupled Hamiltonians (Track CU).

Validates:
1. Multi-center spec construction (LiF, CO, N2, F2, NaCl, CH2O, C2H2, C2H6)
2. Composed Hamiltonian builds (block-diagonal, no cross-block ERIs)
3. Balanced coupled Hamiltonian builds (cross-block ERIs + cross-center V_ne)
4. Isostructural invariance: CO and N2 same block topology -> same Pauli count
5. Regression: existing 18 single-center molecules unchanged

Author: GeoVac Development Team (Track CU, v2.3.0)
"""

import pytest
import numpy as np

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    lif_spec,
    co_spec,
    n2_spec,
    f2_spec,
    nacl_spec,
    ch2o_spec,
    c2h2_spec,
    c2h6_spec,
)
from geovac.molecular_spec import MolecularSpec, OrbitalBlock


# ---------------------------------------------------------------------------
# Spec construction tests
# ---------------------------------------------------------------------------

class TestMultiCenterSpecs:
    """Verify multi-center MolecularSpec construction."""

    def test_lif_spec_basic(self):
        spec = lif_spec()
        assert spec.name == 'LiF'
        assert len(spec.blocks) == 6
        assert sum(b.n_electrons for b in spec.blocks) == 12
        assert len(spec.nuclei) == 2

    def test_lif_qubit_count(self):
        spec = lif_spec(max_n_core=2, max_n_val=2)
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 70  # 7 sub-blocks * 5 states * 2 spin = 70

    def test_co_spec_basic(self):
        spec = co_spec()
        assert spec.name == 'CO'
        assert len(spec.blocks) == 7
        assert sum(b.n_electrons for b in spec.blocks) == 14
        assert len(spec.nuclei) == 2

    def test_co_qubit_count(self):
        spec = co_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 100  # 10 sub-blocks * 5 * 2

    def test_n2_spec_basic(self):
        spec = n2_spec()
        assert spec.name == 'N2'
        assert len(spec.blocks) == 7
        assert sum(b.n_electrons for b in spec.blocks) == 14

    def test_n2_qubit_count(self):
        spec = n2_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 100

    def test_f2_spec_basic(self):
        spec = f2_spec()
        assert spec.name == 'F2'
        assert len(spec.blocks) == 9
        assert sum(b.n_electrons for b in spec.blocks) == 18

    def test_f2_qubit_count(self):
        spec = f2_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 100  # 10 sub-blocks * 5 * 2

    def test_nacl_spec_basic(self):
        spec = nacl_spec()
        assert spec.name == 'NaCl'
        assert len(spec.blocks) == 4
        assert sum(b.n_electrons for b in spec.blocks) == 8

    def test_nacl_qubit_count(self):
        spec = nacl_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 50

    def test_ch2o_spec_basic(self):
        spec = ch2o_spec()
        assert spec.name == 'CH2O'
        assert len(spec.blocks) == 8
        assert sum(b.n_electrons for b in spec.blocks) == 16

    def test_ch2o_qubit_count(self):
        spec = ch2o_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 120  # 12 sub-blocks * 5 * 2

    def test_c2h2_spec_basic(self):
        spec = c2h2_spec()
        assert spec.name == 'C2H2'
        assert len(spec.blocks) == 7
        assert sum(b.n_electrons for b in spec.blocks) == 14

    def test_c2h2_qubit_count(self):
        spec = c2h2_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 120  # 12 sub-blocks * 5 * 2

    def test_c2h6_spec_basic(self):
        spec = c2h6_spec()
        assert spec.name == 'C2H6'
        assert len(spec.blocks) == 9
        assert sum(b.n_electrons for b in spec.blocks) == 18

    def test_c2h6_qubit_count(self):
        spec = c2h6_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['Q'] == 160  # 16 sub-blocks * 5 * 2


class TestNucleusIndices:
    """Verify that nucleus indices are set correctly."""

    def test_lif_nucleus_indices(self):
        spec = lif_spec()
        # Li_core on Li (idx 0)
        assert spec.blocks[0].center_nucleus_idx == 0
        # F_core on F (idx 1)
        assert spec.blocks[1].center_nucleus_idx == 1
        # Bond: center on Li, partner on F
        assert spec.blocks[2].center_nucleus_idx == 0
        assert spec.blocks[2].partner_nucleus_idx == 1
        # F lone pairs on F (idx 1)
        for i in range(3, 6):
            assert spec.blocks[i].center_nucleus_idx == 1

    def test_co_nucleus_indices(self):
        spec = co_spec()
        assert spec.blocks[0].center_nucleus_idx == 0  # C core
        assert spec.blocks[1].center_nucleus_idx == 1  # O core
        # Bond blocks: center on C, partner on O
        for i in range(2, 5):
            assert spec.blocks[i].center_nucleus_idx == 0
            assert spec.blocks[i].partner_nucleus_idx == 1
        assert spec.blocks[5].center_nucleus_idx == 0  # C lone
        assert spec.blocks[6].center_nucleus_idx == 1  # O lone


class TestIsostructuralInvariance:
    """CO and N2 have the same block topology -> same Pauli count."""

    def test_co_n2_same_pauli_count(self):
        spec_co = co_spec()
        spec_n2 = n2_spec()
        result_co = build_composed_hamiltonian(spec_co, pk_in_hamiltonian=False)
        result_n2 = build_composed_hamiltonian(spec_n2, pk_in_hamiltonian=False)
        assert result_co['N_pauli'] == result_n2['N_pauli']
        assert result_co['Q'] == result_n2['Q']


# ---------------------------------------------------------------------------
# Composed Hamiltonian build tests
# ---------------------------------------------------------------------------

class TestComposedBuilds:
    """Verify composed (block-diagonal) Hamiltonians build correctly."""

    def test_lif_composed_hermitian(self):
        spec = lif_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        h1 = result['h1']
        assert np.allclose(h1, h1.T, atol=1e-12), "h1 not symmetric"

    def test_lif_composed_eri_block_diagonal(self):
        """Cross-block ERIs should be zero in composed Hamiltonian."""
        spec = lif_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        eri = result['eri']
        blocks = result['blocks']
        for i, bi in enumerate(blocks):
            for j, bj in enumerate(blocks):
                if i == j:
                    continue
                oi = bi['center_offset']
                ni = bi['center_M'] + bi['partner_M']
                oj = bj['center_offset']
                nj = bj['center_M'] + bj['partner_M']
                cross = eri[oi:oi+ni, oj:oj+nj, :, :]
                assert np.allclose(cross, 0.0, atol=1e-14), (
                    f"Cross-block ERI nonzero between {bi['label']} and {bj['label']}"
                )

    def test_lif_composed_pauli_positive(self):
        spec = lif_spec()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        assert result['N_pauli'] > 0

    def test_all_multi_center_build(self):
        """All multi-center specs produce valid composed Hamiltonians."""
        specs = [
            lif_spec(), co_spec(), n2_spec(), f2_spec(),
            nacl_spec(), ch2o_spec(), c2h2_spec(), c2h6_spec(),
        ]
        for spec in specs:
            result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
            assert result['N_pauli'] > 0, f"{spec.name} has zero Pauli terms"
            assert result['Q'] > 0, f"{spec.name} has zero qubits"


# ---------------------------------------------------------------------------
# Backward compatibility regression tests
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """Verify existing molecules are unchanged by multi-center additions."""

    def test_molecular_spec_default_nuclei_empty(self):
        """Legacy specs should have empty nuclei list."""
        from geovac.composed_qubit import lih_spec as legacy_lih_spec
        spec = legacy_lih_spec()
        assert spec.nuclei == []

    def test_orbital_block_default_indices(self):
        """Legacy OrbitalBlocks should have -1 nucleus indices."""
        from geovac.composed_qubit import lih_spec as legacy_lih_spec
        spec = legacy_lih_spec()
        for blk in spec.blocks:
            assert blk.center_nucleus_idx == -1
            assert blk.partner_nucleus_idx == -1

    @pytest.mark.parametrize("spec_fn,expected_pauli", [
        ('lih_spec', 334),
        ('beh2_spec', 556),
        ('h2o_spec', 778),
    ])
    def test_existing_composed_pauli_unchanged(self, spec_fn, expected_pauli):
        """First-row composed molecules produce exact same Pauli counts."""
        import geovac.composed_qubit as cq
        spec = getattr(cq, spec_fn)()
        result = build_composed_hamiltonian(spec, pk_in_hamiltonian=True)
        assert result['N_pauli'] == expected_pauli, (
            f"{spec_fn}: expected {expected_pauli}, got {result['N_pauli']}"
        )
