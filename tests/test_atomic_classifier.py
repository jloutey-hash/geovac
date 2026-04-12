"""
Tests for Atomic Classifier (Track BG, v2.0.30)
=================================================

Validates classify_atom() for Z=1-10, Z^2 PK scaling, edge cases,
and consistency with composed_qubit.py PK defaults.

Author: GeoVac Development Team
Date: April 2026
"""

import pytest

from geovac.atomic_classifier import (
    AtomClassification,
    classify_atom,
    pk_params_z2_scaled,
    pk_params_from_formulas,
)


# ---------------------------------------------------------------------------
# Structure type tests (Paper 16)
# ---------------------------------------------------------------------------

class TestStructureTypes:
    """Each Z=1-10 maps to the correct Paper 16 structure type."""

    @pytest.mark.parametrize("Z, expected", [
        (1, 'A'),
        (2, 'B'),
        (3, 'C'),
        (4, 'D'),
        (5, 'E'),
        (6, 'E'),
        (7, 'E'),
        (8, 'E'),
        (9, 'E'),
        (10, 'B'),
    ])
    def test_structure_type(self, Z: int, expected: str) -> None:
        c = classify_atom(Z)
        assert c.structure_type == expected, f"Z={Z}: expected {expected}, got {c.structure_type}"

    @pytest.mark.parametrize("Z, expected", [
        (1, 'hydrogen'),
        (2, 'noble_gas'),
        (3, 'alkali_metal'),
        (4, 'alkaline_earth'),
        (5, 'p_block'),
        (6, 'p_block'),
        (7, 'p_block'),
        (8, 'p_block'),
        (9, 'p_block'),
        (10, 'noble_gas'),
    ])
    def test_group_type(self, Z: int, expected: str) -> None:
        c = classify_atom(Z)
        assert c.group_type == expected


# ---------------------------------------------------------------------------
# Z_eff and electron counts
# ---------------------------------------------------------------------------

class TestElectronConfig:
    """Electron configuration and Z_eff for each Z."""

    @pytest.mark.parametrize("Z, n_core, n_val, Z_eff", [
        (1, 0, 1, 1.0),
        (2, 0, 2, 2.0),
        (3, 2, 1, 1.0),
        (4, 2, 2, 2.0),
        (5, 2, 3, 3.0),
        (6, 2, 4, 4.0),
        (7, 2, 5, 5.0),
        (8, 2, 6, 6.0),
        (9, 2, 7, 7.0),
        (10, 2, 8, 8.0),
    ])
    def test_electron_counts(self, Z: int, n_core: int, n_val: int, Z_eff: float) -> None:
        c = classify_atom(Z)
        assert c.N_electrons == Z
        assert c.n_core_electrons == n_core
        assert c.n_valence_electrons == n_val
        assert c.Z_eff_valence == pytest.approx(Z_eff)

    @pytest.mark.parametrize("Z, core_cfg, val_cfg", [
        (1, 'none', '1s1'),
        (2, 'none', '1s2'),
        (3, '1s2', '2s1'),
        (4, '1s2', '2s2'),
        (5, '1s2', '2s2 2p1'),
        (6, '1s2', '2s2 2p2'),
        (7, '1s2', '2s2 2p3'),
        (8, '1s2', '2s2 2p4'),
        (9, '1s2', '2s2 2p5'),
        (10, '1s2', '2s2 2p6'),
    ])
    def test_configurations(self, Z: int, core_cfg: str, val_cfg: str) -> None:
        c = classify_atom(Z)
        assert c.core_config == core_cfg
        assert c.valence_config == val_cfg


# ---------------------------------------------------------------------------
# Period, nu, mu_free
# ---------------------------------------------------------------------------

class TestQuantumNumbers:
    """Period, angular quantum number nu, and Pauli centrifugal cost."""

    @pytest.mark.parametrize("Z, period", [
        (1, 1), (2, 1),
        (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (8, 2), (9, 2), (10, 2),
    ])
    def test_period(self, Z: int, period: int) -> None:
        assert classify_atom(Z).period == period

    @pytest.mark.parametrize("Z, nu", [
        (1, 0),   # N=1, nu = max(1-2, 0) = 0
        (2, 0),   # N=2, nu = max(2-2, 0) = 0
        (3, 1),   # N=3, nu = 1
        (4, 2),   # N=4, nu = 2
        (5, 3),
        (6, 4),
        (7, 5),
        (8, 6),
        (9, 7),
        (10, 8),
    ])
    def test_nu(self, Z: int, nu: int) -> None:
        assert classify_atom(Z).nu == nu

    @pytest.mark.parametrize("Z, mu_free", [
        (1, 0.0),
        (2, 0.0),
        (3, 2.0),     # 2 * 1^2
        (4, 8.0),     # 2 * 2^2
        (5, 18.0),    # 2 * 3^2
        (6, 32.0),    # 2 * 4^2
        (7, 50.0),    # 2 * 5^2
        (8, 72.0),    # 2 * 6^2
        (9, 98.0),    # 2 * 7^2
        (10, 128.0),  # 2 * 8^2
    ])
    def test_mu_free(self, Z: int, mu_free: float) -> None:
        assert classify_atom(Z).mu_free == pytest.approx(mu_free)


# ---------------------------------------------------------------------------
# PK parameters
# ---------------------------------------------------------------------------

class TestPKParams:
    """PK parameters match Paper 17 Table 1 / Z^2 scaling."""

    def test_h_no_pk(self) -> None:
        c = classify_atom(1)
        assert c.pk_params is None
        assert c.pk_source == 'none'

    def test_he_no_pk(self) -> None:
        c = classify_atom(2)
        assert c.pk_params is None
        assert c.pk_source == 'none'

    def test_li_computed(self) -> None:
        """Li PK from Paper 17 Table 1."""
        c = classify_atom(3)
        assert c.pk_source == 'computed'
        assert c.pk_params['A'] == pytest.approx(6.93)
        assert c.pk_params['B'] == pytest.approx(7.00)

    def test_be_computed(self) -> None:
        """Be PK from Paper 17 Table 1 (computed, not Z^2 scaled)."""
        c = classify_atom(4)
        assert c.pk_source == 'computed'
        assert c.pk_params['A'] == pytest.approx(13.01)
        assert c.pk_params['B'] == pytest.approx(12.53)

    def test_oxygen_z2_scaled(self) -> None:
        """O (Z=8) PK matches _PK_HELIKE_DEFAULTS in composed_qubit.py."""
        c = classify_atom(8)
        assert c.pk_source == 'z2_scaled'
        assert c.pk_params['A'] == pytest.approx(49.28, abs=0.01)
        assert c.pk_params['B'] == pytest.approx(49.78, abs=0.01)

    @pytest.mark.parametrize("Z, A_expected, B_expected", [
        (5, 21.40, 18.46),   # ab initio (Track BI)
        (6, 31.37, 25.54),   # ab initio (Track BI)
        (7, 43.09, 33.05),   # ab initio (Track BI)
        (9, 71.80, 48.61),   # ab initio (Track BI)
    ])
    def test_computed_pk_values(self, Z: int, A_expected: float, B_expected: float) -> None:
        c = classify_atom(Z)
        assert c.pk_source == 'computed'
        assert c.pk_params['A'] == pytest.approx(A_expected, rel=1e-3)
        assert c.pk_params['B'] == pytest.approx(B_expected, rel=1e-3)

    def test_neon_z2_scaled(self) -> None:
        """Ne (Z=10) still uses Z^2 scaling (not run in Track BI)."""
        c = classify_atom(10)
        assert c.pk_source == 'z2_scaled'
        assert c.pk_params['A'] == pytest.approx(6.93 * (10/3)**2, rel=1e-6)
        assert c.pk_params['B'] == pytest.approx(7.00 * (10/3)**2, rel=1e-6)


# ---------------------------------------------------------------------------
# Z^2 scaling function
# ---------------------------------------------------------------------------

class TestPKScaling:
    """pk_params_z2_scaled() produces correct values."""

    def test_li_anchor(self) -> None:
        """At Z=3 (the anchor), Z^2 scaling reproduces exact values."""
        p = pk_params_z2_scaled(3)
        assert p['A'] == pytest.approx(6.93)
        assert p['B'] == pytest.approx(7.00)

    def test_be_z2(self) -> None:
        p = pk_params_z2_scaled(4)
        assert p['A'] == pytest.approx(6.93 * (4/3)**2, rel=1e-10)
        assert p['B'] == pytest.approx(7.00 * (4/3)**2, rel=1e-10)

    def test_oxygen_z2(self) -> None:
        p = pk_params_z2_scaled(8)
        assert p['A'] == pytest.approx(6.93 * (8/3)**2, rel=1e-10)
        assert p['B'] == pytest.approx(7.00 * (8/3)**2, rel=1e-10)

    def test_scaling_is_quadratic(self) -> None:
        """Verify A(2Z) / A(Z) = 4."""
        p3 = pk_params_z2_scaled(3)
        p6 = pk_params_z2_scaled(6)
        assert p6['A'] / p3['A'] == pytest.approx(4.0, rel=1e-10)
        assert p6['B'] / p3['B'] == pytest.approx(4.0, rel=1e-10)


# ---------------------------------------------------------------------------
# Ab initio PK formula
# ---------------------------------------------------------------------------

class TestPKFromFormulas:
    """pk_params_from_formulas() implements Paper 17 Sec IV."""

    def test_basic(self) -> None:
        p = pk_params_from_formulas(Z=3, E_core=-4.5, r_eff=0.5)
        assert p['A'] == pytest.approx(9.0)
        assert p['B'] == pytest.approx(4.0)

    def test_returns_dict(self) -> None:
        p = pk_params_from_formulas(Z=3, E_core=-1.0, r_eff=1.0)
        assert isinstance(p, dict)
        assert 'A' in p and 'B' in p


# ---------------------------------------------------------------------------
# Consistency with composed_qubit.py defaults
# ---------------------------------------------------------------------------

class TestComposedQubitConsistency:
    """Classifier PK values must match the defaults in composed_qubit.py."""

    def test_li_matches_pk_defaults(self) -> None:
        """Li classifier matches _PK_DEFAULTS[3] in composed_qubit.py."""
        c = classify_atom(3)
        # Values from _PK_DEFAULTS in composed_qubit.py
        assert c.pk_params['A'] == pytest.approx(6.93)
        assert c.pk_params['B'] == pytest.approx(7.00)

    def test_be_matches_pk_defaults(self) -> None:
        """Be classifier matches _PK_DEFAULTS[4] in composed_qubit.py."""
        c = classify_atom(4)
        # Values from _PK_DEFAULTS in composed_qubit.py
        assert c.pk_params['A'] == pytest.approx(13.01)
        assert c.pk_params['B'] == pytest.approx(12.53)

    def test_oxygen_matches_helike_defaults(self) -> None:
        """O classifier matches _PK_HELIKE_DEFAULTS[8] in composed_qubit.py."""
        c = classify_atom(8)
        # Values from _PK_HELIKE_DEFAULTS in composed_qubit.py
        assert c.pk_params['A'] == pytest.approx(49.28, abs=0.01)
        assert c.pk_params['B'] == pytest.approx(49.78, abs=0.01)


# ---------------------------------------------------------------------------
# Edge cases and error handling
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Invalid inputs and unsupported atoms."""

    def test_z_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="Z must be >= 1"):
            classify_atom(0)

    def test_z_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="Z must be >= 1"):
            classify_atom(-1)

    @pytest.mark.parametrize("Z", [11, 12, 15, 20, 50])
    def test_unsupported_z(self, Z: int) -> None:
        c = classify_atom(Z)
        assert c.supported is False
        assert c.support_note != ''
        assert c.pk_params is None

    def test_unsupported_has_correct_N(self) -> None:
        c = classify_atom(11)
        assert c.N_electrons == 11
        assert c.nu == 9  # 11 - 2
        assert c.mu_free == pytest.approx(2.0 * 81)  # 2 * 9^2

    def test_all_supported_atoms_return_true(self) -> None:
        for Z in range(1, 11):
            c = classify_atom(Z)
            assert c.supported is True, f"Z={Z} should be supported"


# ---------------------------------------------------------------------------
# Dataclass structure
# ---------------------------------------------------------------------------

class TestDataclass:
    """AtomClassification is a proper dataclass."""

    def test_is_dataclass(self) -> None:
        c = classify_atom(1)
        assert isinstance(c, AtomClassification)

    def test_fields_accessible(self) -> None:
        c = classify_atom(6)
        # Just verify all fields exist and are the right types
        assert isinstance(c.Z, int)
        assert isinstance(c.N_electrons, int)
        assert isinstance(c.structure_type, str)
        assert isinstance(c.nu, int)
        assert isinstance(c.mu_free, float)
        assert isinstance(c.pk_params, dict)
        assert isinstance(c.supported, bool)

    def test_default_supported_true(self) -> None:
        """Default for supported is True."""
        c = classify_atom(1)
        assert c.supported is True
        assert c.support_note == ''
