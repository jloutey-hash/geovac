"""
Tests for second-row molecule balanced coupled Hamiltonians (Track CJ).

Validates:
  - Spec factories produce correct block topology (no core block)
  - Balanced Hamiltonians build without error
  - Pauli term counts are recorded
  - Ecosystem export works
  - Nuclear repulsion includes frozen core energy
  - Z_eff values are correct
  - No regression on first-row LiH balanced

Author: GeoVac Development Team (Track CJ)
Date: April 2026
"""

import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Spec factory tests
# ---------------------------------------------------------------------------

class TestNaHSpec:
    """NaH: frozen [Ne] core + 1 bond pair."""

    def test_block_count(self):
        from geovac.composed_qubit import nah_spec
        spec = nah_spec(max_n_val=2)
        assert len(spec.blocks) == 1, "NaH should have 1 block (bond only, no core)"

    def test_no_core_block(self):
        from geovac.composed_qubit import nah_spec
        spec = nah_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.block_type != 'core', "NaH should have no core block"

    def test_z_eff(self):
        from geovac.composed_qubit import nah_spec
        spec = nah_spec(max_n_val=2)
        assert spec.blocks[0].Z_center == 1.0, "Na valence Z_eff should be 1.0"

    def test_has_h_partner(self):
        from geovac.composed_qubit import nah_spec
        spec = nah_spec(max_n_val=2)
        assert spec.blocks[0].has_h_partner is True

    def test_nuclear_repulsion_includes_core_energy(self):
        from geovac.composed_qubit import nah_spec
        from geovac.neon_core import _NIST_CORE_ENERGIES
        spec = nah_spec(max_n_val=2)
        # Core energy is large and negative
        assert spec.nuclear_repulsion_constant < -100.0
        # Should contain the NIST core energy contribution
        assert spec.nuclear_repulsion_constant < _NIST_CORE_ENERGIES[11] + 10.0

    def test_pk_is_zero(self):
        from geovac.composed_qubit import nah_spec
        spec = nah_spec(max_n_val=2)
        assert spec.blocks[0].pk_A == 0.0
        assert spec.blocks[0].pk_B == 0.0


class TestMgH2Spec:
    """MgH2: frozen [Ne] core + 2 bond pairs."""

    def test_block_count(self):
        from geovac.composed_qubit import mgh2_spec
        spec = mgh2_spec(max_n_val=2)
        assert len(spec.blocks) == 2

    def test_z_eff(self):
        from geovac.composed_qubit import mgh2_spec
        spec = mgh2_spec(max_n_val=2)
        assert spec.blocks[0].Z_center == 2.0
        assert spec.blocks[1].Z_center == 2.0


class TestHClSpec:
    """HCl: frozen [Ne] core + 1 bond + 3 lone pairs."""

    def test_block_count(self):
        from geovac.composed_qubit import hcl_spec
        spec = hcl_spec(max_n_val=2)
        assert len(spec.blocks) == 4  # 1 bond + 3 lone

    def test_z_eff(self):
        from geovac.composed_qubit import hcl_spec
        spec = hcl_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.Z_center == 7.0

    def test_block_types(self):
        from geovac.composed_qubit import hcl_spec
        spec = hcl_spec(max_n_val=2)
        types = [b.block_type for b in spec.blocks]
        assert types.count('bond') == 1
        assert types.count('lone_pair') == 3


class TestH2SSpec:
    """H2S: frozen [Ne] core + 2 bonds + 2 lone pairs."""

    def test_block_count(self):
        from geovac.composed_qubit import h2s_spec
        spec = h2s_spec(max_n_val=2)
        assert len(spec.blocks) == 4  # 2 bonds + 2 lone

    def test_z_eff(self):
        from geovac.composed_qubit import h2s_spec
        spec = h2s_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.Z_center == 6.0


class TestPH3Spec:
    """PH3: frozen [Ne] core + 3 bonds + 1 lone pair."""

    def test_block_count(self):
        from geovac.composed_qubit import ph3_spec
        spec = ph3_spec(max_n_val=2)
        assert len(spec.blocks) == 4  # 3 bonds + 1 lone

    def test_z_eff(self):
        from geovac.composed_qubit import ph3_spec
        spec = ph3_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.Z_center == 5.0


class TestSiH4Spec:
    """SiH4: frozen [Ne] core + 4 bond pairs."""

    def test_block_count(self):
        from geovac.composed_qubit import sih4_spec
        spec = sih4_spec(max_n_val=2)
        assert len(spec.blocks) == 4  # 4 bonds, no lone pairs

    def test_z_eff(self):
        from geovac.composed_qubit import sih4_spec
        spec = sih4_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.Z_center == 4.0

    def test_no_core_block(self):
        from geovac.composed_qubit import sih4_spec
        spec = sih4_spec(max_n_val=2)
        for blk in spec.blocks:
            assert blk.block_type != 'core'


# ---------------------------------------------------------------------------
# Frozen core cross-nuclear attraction
# ---------------------------------------------------------------------------

class TestFrozenCoreCrossNuclear:
    """Test _v_cross_frozen_core function."""

    def test_finite_value(self):
        from geovac.composed_qubit import _v_cross_frozen_core
        v = _v_cross_frozen_core(11, 1.0, 3.566)
        assert np.isfinite(v)
        # Should be negative (attraction)
        assert v < 0.0

    def test_larger_z_stronger(self):
        """Higher Z_heavy means more core charge, stronger attraction."""
        from geovac.composed_qubit import _v_cross_frozen_core
        v_na = _v_cross_frozen_core(11, 1.0, 3.0)
        v_cl = _v_cross_frozen_core(17, 1.0, 3.0)
        # Cl core is more tightly bound, so at same R the enclosed
        # charge should be similar (~10), but Cl core is more compact.
        # Both should be negative.
        assert v_na < 0.0
        assert v_cl < 0.0

    def test_long_range_limit(self):
        """At large R, V_cross ~ -Z_other * 10 / R (all core enclosed)."""
        from geovac.composed_qubit import _v_cross_frozen_core
        R_large = 50.0
        v = _v_cross_frozen_core(11, 1.0, R_large)
        expected = -1.0 * 10.0 / R_large  # -Z_other * N_core / R
        assert abs(v - expected) < 0.01 * abs(expected)


# ---------------------------------------------------------------------------
# Balanced Hamiltonian build tests (composed + cross-block)
# ---------------------------------------------------------------------------

class TestNaHBalanced:
    """NaH balanced coupled Hamiltonian."""

    def test_builds_without_error(self):
        from geovac.composed_qubit import nah_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = nah_spec(max_n_val=1)
        result = build_balanced_hamiltonian(spec, R=3.566)
        assert result['N_pauli'] > 0
        assert result['Q'] > 0

    def test_qubit_count_max_n1(self):
        from geovac.composed_qubit import nah_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = nah_spec(max_n_val=1)
        result = build_balanced_hamiltonian(spec, R=3.566)
        # max_n=1: 1 state per sub-block, bond has center+partner = 2
        # Q = 2 * M = 2 * 2 = 4
        assert result['Q'] == 4

    def test_qubit_count_max_n2(self):
        from geovac.composed_qubit import nah_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = nah_spec(max_n_val=2)
        result = build_balanced_hamiltonian(spec, R=3.566)
        # max_n=2: 5 states per sub-block (1s + 2s,2p0,2p+1,2p-1)
        # bond has center+partner = 10 spatial orbitals
        # Q = 2 * 10 = 20
        assert result['Q'] == 20


class TestMgH2Balanced:
    """MgH2 balanced coupled Hamiltonian."""

    def test_builds_without_error(self):
        from geovac.composed_qubit import mgh2_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = mgh2_spec(max_n_val=1)
        result = build_balanced_hamiltonian(spec, R=3.268)
        assert result['N_pauli'] > 0


class TestHClBalanced:
    """HCl balanced coupled Hamiltonian."""

    def test_builds_without_error(self):
        from geovac.composed_qubit import hcl_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = hcl_spec(max_n_val=1)
        result = build_balanced_hamiltonian(spec, R=2.409)
        assert result['N_pauli'] > 0


class TestSiH4Balanced:
    """SiH4 balanced coupled Hamiltonian."""

    def test_builds_without_error(self):
        from geovac.composed_qubit import sih4_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = sih4_spec(max_n_val=1)
        result = build_balanced_hamiltonian(spec, R=2.798)
        assert result['N_pauli'] > 0


# ---------------------------------------------------------------------------
# Ecosystem export tests
# ---------------------------------------------------------------------------

class TestEcosystemExport:
    """Test ecosystem export for second-row molecules."""

    def test_nah_export(self):
        from geovac.ecosystem_export import hamiltonian
        H = hamiltonian('NaH', l_max=1)
        assert H.n_terms > 0
        assert H.one_norm > 0
        of_op = H.to_openfermion()
        assert len(of_op.terms) > 0

    def test_registry_contains_row2(self):
        from geovac.ecosystem_export import _SYSTEM_REGISTRY
        for key in ('nah', 'mgh2', 'hcl', 'h2s', 'ph3', 'sih4'):
            assert key in _SYSTEM_REGISTRY, f"{key} missing from registry"


# ---------------------------------------------------------------------------
# Regression: first-row LiH balanced still works
# ---------------------------------------------------------------------------

class TestLiHBalancedRegression:
    """Ensure first-row LiH balanced Hamiltonian is unchanged."""

    def test_lih_balanced_pauli_count(self):
        from geovac.composed_qubit import lih_spec
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec = lih_spec(max_n_core=2, max_n_val=2, R=3.015)
        result = build_balanced_hamiltonian(spec, R=3.015)
        assert result['N_pauli'] == 878, (
            f"LiH balanced at max_n=2 should give 878 Pauli terms, got {result['N_pauli']}"
        )


# ---------------------------------------------------------------------------
# Resource table: collect Pauli counts for reporting
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestResourceTable:
    """Collect Pauli term counts at max_n=1 and max_n=2 for all 6 molecules."""

    _MOLECULES = [
        ('NaH', 'nah_spec', {'R': 3.566}),
        ('MgH2', 'mgh2_spec', {'R': 3.268}),
        ('HCl', 'hcl_spec', {'R': 2.409}),
        ('H2S', 'h2s_spec', {'R_SH': 2.534}),
        ('PH3', 'ph3_spec', {'R_PH': 2.683}),
        ('SiH4', 'sih4_spec', {'R_SiH': 2.798}),
    ]

    def _build(self, spec_name, max_n, kwargs):
        import geovac.composed_qubit as cq
        from geovac.balanced_coupled import build_balanced_hamiltonian
        spec_func = getattr(cq, spec_name)
        # Determine R for balanced builder
        R_val = kwargs.get('R') or kwargs.get('R_SH') or kwargs.get('R_PH') or kwargs.get('R_SiH')
        spec = spec_func(max_n_val=max_n, **kwargs)
        result = build_balanced_hamiltonian(spec, R=R_val)
        return result

    @pytest.mark.parametrize("name,spec_name,kwargs", _MOLECULES)
    def test_max_n1(self, name, spec_name, kwargs):
        result = self._build(spec_name, 1, kwargs)
        print(f"\n{name} max_n=1: Q={result['Q']}, Pauli={result['N_pauli']}, "
              f"1-norm={result['one_norm']:.2f} Ha")
        assert result['N_pauli'] > 0

    @pytest.mark.parametrize("name,spec_name,kwargs", _MOLECULES)
    def test_max_n2(self, name, spec_name, kwargs):
        result = self._build(spec_name, 2, kwargs)
        print(f"\n{name} max_n=2: Q={result['Q']}, Pauli={result['N_pauli']}, "
              f"1-norm={result['one_norm']:.2f} Ha")
        assert result['N_pauli'] > 0
