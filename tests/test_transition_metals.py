"""
Tests for transition metal classifier (Z=21-30) and hydride Hamiltonians.

Track CZ/DA extension — v2.8.0
"""

import pytest
import numpy as np


# ---------------------------------------------------------------------------
# _enumerate_states with l_min
# ---------------------------------------------------------------------------

class TestEnumerateStatesLmin:
    def test_l_min_0_unchanged(self):
        """l_min=0 gives the same result as the original function."""
        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(3, l_min=0)
        # n=1: (1,0,0); n=2: (2,0,0),(2,1,-1),(2,1,0),(2,1,1);
        # n=3: 1+3+5=9 => total 1+4+9=14
        assert len(states) == 14

    def test_l_min_2_d_only(self):
        """l_min=2 with max_n=3 gives exactly the 5 d-orbitals."""
        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(3, l_min=2)
        assert len(states) == 5
        assert all(l == 2 for n, l, m in states)
        assert all(n == 3 for n, l, m in states)
        ms = [m for _, _, m in states]
        assert ms == [-2, -1, 0, 1, 2]

    def test_l_min_1_p_and_d(self):
        """l_min=1 excludes s-orbitals."""
        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(3, l_min=1)
        # n=1: nothing; n=2: 3 (p); n=3: 3+5=8 => total 11
        assert len(states) == 11
        assert all(l >= 1 for _, l, _ in states)

    def test_l_min_too_high_empty(self):
        """l_min >= max_n gives empty list."""
        from geovac.composed_qubit import _enumerate_states
        states = _enumerate_states(2, l_min=2)
        assert len(states) == 0


# ---------------------------------------------------------------------------
# Atomic classifier Z=11-30
# ---------------------------------------------------------------------------

class TestAtomicClassifierExtended:
    @pytest.mark.parametrize("Z", range(11, 19))
    def test_second_row_supported(self, Z):
        """Second-row atoms (Z=11-18) are classified and supported."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(Z)
        assert c.supported
        assert c.n_core_electrons == 10
        assert c.period == 3
        assert c.pk_source == 'frozen_core'
        assert c.pk_params is None

    @pytest.mark.parametrize("Z", [19, 20])
    def test_fourth_row_sblock(self, Z):
        """K and Ca are classified with [Ar] core."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(Z)
        assert c.supported
        assert c.n_core_electrons == 18
        assert c.period == 4

    @pytest.mark.parametrize("Z", range(21, 31))
    def test_transition_metals_supported(self, Z):
        """Transition metals (Z=21-30) are classified as type F."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(Z)
        assert c.supported
        assert c.structure_type == 'F'
        assert c.n_core_electrons == 18
        assert c.period == 4
        assert c.group_type == 'd_block'
        assert c.pk_source == 'frozen_core'
        assert c.pk_params is None
        assert c.n_valence_electrons == Z - 18

    def test_cr_anomalous(self):
        """Cr (Z=24) has anomalous 3d5 4s1 configuration."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(24)
        assert '3d5' in c.valence_config
        assert '4s1' in c.valence_config

    def test_cu_anomalous(self):
        """Cu (Z=29) has anomalous 3d10 4s1 configuration."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(29)
        assert '3d10' in c.valence_config
        assert '4s1' in c.valence_config

    def test_z31_unsupported(self):
        """Z=31 (Ga) is beyond the classifier range."""
        from geovac.atomic_classifier import classify_atom
        c = classify_atom(31)
        assert not c.supported


# ---------------------------------------------------------------------------
# Transition metal hydride specs
# ---------------------------------------------------------------------------

class TestTMHydrideSpecs:
    @pytest.mark.parametrize("Z", range(21, 31))
    def test_spec_creation(self, Z):
        """All 10 TM hydride specs can be created."""
        from geovac.molecular_spec import transition_metal_hydride_spec
        spec = transition_metal_hydride_spec(Z)
        assert len(spec.blocks) == 2
        assert spec.blocks[0].block_type == 'bond'
        assert spec.blocks[1].block_type == 'lone_pair'
        assert spec.blocks[1].l_min == 2

    def test_cr_single_bond_electron(self):
        """CrH bond block has 1 electron (from 4s1)."""
        from geovac.molecular_spec import transition_metal_hydride_spec
        spec = transition_metal_hydride_spec(24)
        assert spec.blocks[0].n_electrons == 1  # 4s1
        assert spec.blocks[1].n_electrons == 5  # 3d5

    def test_cu_single_bond_electron(self):
        """CuH bond block has 1 electron (from 4s1)."""
        from geovac.molecular_spec import transition_metal_hydride_spec
        spec = transition_metal_hydride_spec(29)
        assert spec.blocks[0].n_electrons == 1   # 4s1
        assert spec.blocks[1].n_electrons == 10  # 3d10

    def test_convenience_aliases(self):
        """All convenience aliases work."""
        from geovac.molecular_spec import (
            sch_spec, tih_spec, vh_spec, crh_spec, mnh_spec,
            feh_spec, coh_spec, nih_spec, cuh_spec, znh_spec,
        )
        specs = [sch_spec(), tih_spec(), vh_spec(), crh_spec(), mnh_spec(),
                 feh_spec(), coh_spec(), nih_spec(), cuh_spec(), znh_spec()]
        names = [s.name for s in specs]
        assert names == ['ScH', 'TiH', 'VH', 'CrH', 'MnH',
                         'FeH', 'CoH', 'NiH', 'CuH', 'ZnH']


# ---------------------------------------------------------------------------
# Full Hamiltonian builds (slower — require ERI computation)
# ---------------------------------------------------------------------------

class TestTMHydrideHamiltonians:
    @pytest.mark.parametrize("Z", [21, 24, 29, 30])
    def test_build_representative_hydrides(self, Z):
        """Representative TM hydrides build successfully via general builder."""
        from geovac.molecular_spec import transition_metal_hydride_spec
        from geovac.composed_qubit import build_composed_hamiltonian

        spec = transition_metal_hydride_spec(Z, max_n_bond=2, max_n_d=3)
        result = build_composed_hamiltonian(spec, verbose=False)

        assert result['M'] > 0
        assert result['Q'] == 2 * result['M']
        assert result['N_pauli'] > 0

        # d-block has 5 orbitals (l_min=2, max_n=3)
        # Bond block has 5 center + 5 partner = 10 orbitals
        # Total: 15 spatial orbitals, 30 qubits
        assert result['M'] == 15
        assert result['Q'] == 30

    def test_d_block_eri_sparsity(self):
        """d-block ERI is sparser than s+p block at same orbital count."""
        from geovac.composed_qubit import (
            _enumerate_states, _compute_rk_integrals_block, _build_eri_block,
        )
        # d-only block (5 orbitals)
        d_states = _enumerate_states(3, l_min=2)
        assert len(d_states) == 5
        rk_d = _compute_rk_integrals_block(3.0, d_states)
        eri_d = _build_eri_block(3.0, d_states, rk_d)
        n_d_eri = len(eri_d)

        # s+p block (same Z, 5 orbitals: max_n=2 gives 1s,2s,2p*3)
        sp_states = _enumerate_states(2)
        assert len(sp_states) == 5
        rk_sp = _compute_rk_integrals_block(3.0, sp_states)
        eri_sp = _build_eri_block(3.0, sp_states, rk_sp)
        n_sp_eri = len(eri_sp)

        # d-block should be sparser (fewer nonzero ERI from Gaunt rules at l=2)
        assert n_d_eri > 0
        assert n_sp_eri > 0
        # Track CZ result: d-block 4.0% vs s+p 8.9% ERI density
        d_density = n_d_eri / (5 ** 4)
        sp_density = n_sp_eri / (5 ** 4)
        assert d_density < sp_density, (
            f"d-block ERI density {d_density:.3f} should be < "
            f"s+p density {sp_density:.3f}"
        )


# ---------------------------------------------------------------------------
# Ecosystem export
# ---------------------------------------------------------------------------

class TestEcosystemExportTM:
    def test_tm_hydrides_in_registry(self):
        """All 10 TM hydrides are in the system registry."""
        from geovac.ecosystem_export import _SYSTEM_REGISTRY
        for name in ['sch', 'tih', 'vh', 'crh', 'mnh',
                     'feh', 'coh', 'nih', 'cuh', 'znh']:
            assert name in _SYSTEM_REGISTRY

    @pytest.mark.parametrize("name", ['ScH', 'TiH'])
    def test_hamiltonian_api(self, name):
        """TM hydrides are accessible via hamiltonian() API."""
        from geovac.ecosystem_export import hamiltonian
        h = hamiltonian(name)
        assert h.n_qubits == 30
        assert h.n_terms > 0
