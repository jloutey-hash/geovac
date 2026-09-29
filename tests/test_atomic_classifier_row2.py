"""
Tests for Atomic Classifier — Second Row Z=11-18 (Track CH, v2.0.45)
=====================================================================

Validates classify_atom() for Z=11-18 (Na through Ar).  All second-row
atoms use frozen-core PK (pk_source='frozen_core', pk_params=None).

Author: GeoVac Development Team
Date: April 2026
"""

import pytest

from geovac.atomic_classifier import AtomClassification, classify_atom


# ---------------------------------------------------------------------------
# Structure type tests (Paper 16)
# ---------------------------------------------------------------------------

class TestStructureTypesRow2:
    """Each Z=11-18 maps to the correct Paper 16 structure type."""

    @pytest.mark.parametrize("Z, expected", [
        (11, 'C'),   # Na: single valence over closed core
        (12, 'D'),   # Mg: s^2 valence over closed core
        (13, 'E'),   # Al: open p-shell
        (14, 'E'),   # Si
        (15, 'E'),   # P
        (16, 'E'),   # S
        (17, 'E'),   # Cl
        (18, 'B'),   # Ar: closed shell
    ])
    def test_structure_type(self, Z: int, expected: str) -> None:
        c = classify_atom(Z)
        assert c.structure_type == expected, f"Z={Z}: expected {expected}, got {c.structure_type}"


# ---------------------------------------------------------------------------
# Group type tests
# ---------------------------------------------------------------------------

class TestGroupTypesRow2:
    """Each Z=11-18 maps to the correct group type."""

    @pytest.mark.parametrize("Z, expected", [
        (11, 'alkali_metal'),
        (12, 'alkaline_earth'),
        (13, 'p_block'),
        (14, 'p_block'),
        (15, 'p_block'),
        (16, 'p_block'),
        (17, 'p_block'),
        (18, 'noble_gas'),
    ])
    def test_group_type(self, Z: int, expected: str) -> None:
        c = classify_atom(Z)
        assert c.group_type == expected, f"Z={Z}: expected {expected}, got {c.group_type}"


# ---------------------------------------------------------------------------
# Electron configuration tests
# ---------------------------------------------------------------------------

class TestElectronConfigRow2:
    """Electron counts, Z_eff, and configurations for Z=11-18."""

    @pytest.mark.parametrize("Z, n_core, n_val, Z_eff", [
        (11, 10, 1, 1.0),
        (12, 10, 2, 2.0),
        (13, 10, 3, 3.0),
        (14, 10, 4, 4.0),
        (15, 10, 5, 5.0),
        (16, 10, 6, 6.0),
        (17, 10, 7, 7.0),
        (18, 10, 8, 8.0),
    ])
    def test_electron_counts(self, Z: int, n_core: int, n_val: int, Z_eff: float) -> None:
        c = classify_atom(Z)
        assert c.N_electrons == Z
        assert c.n_core_electrons == n_core, f"Z={Z}: n_core expected {n_core}, got {c.n_core_electrons}"
        assert c.n_valence_electrons == n_val, f"Z={Z}: n_val expected {n_val}, got {c.n_valence_electrons}"
        assert c.Z_eff_valence == pytest.approx(Z_eff)

    @pytest.mark.parametrize("Z", range(11, 19))
    def test_core_config_neon(self, Z: int) -> None:
        """All second-row atoms have Ne-like core: 1s2 2s2 2p6."""
        c = classify_atom(Z)
        assert c.core_config == '1s2 2s2 2p6', f"Z={Z}: core_config = {c.core_config}"

    @pytest.mark.parametrize("Z, val_cfg", [
        (11, '3s1'),
        (12, '3s2'),
        (13, '3s2 3p1'),
        (14, '3s2 3p2'),
        (15, '3s2 3p3'),
        (16, '3s2 3p4'),
        (17, '3s2 3p5'),
        (18, '3s2 3p6'),
    ])
    def test_valence_config(self, Z: int, val_cfg: str) -> None:
        c = classify_atom(Z)
        assert c.valence_config == val_cfg, f"Z={Z}: valence_config = {c.valence_config}"


# ---------------------------------------------------------------------------
# Period and quantum numbers
# ---------------------------------------------------------------------------

class TestQuantumNumbersRow2:
    """Period, nu, and mu_free for Z=11-18."""

    @pytest.mark.parametrize("Z", range(11, 19))
    def test_period_3(self, Z: int) -> None:
        """All second-row atoms are period 3."""
        assert classify_atom(Z).period == 3

    # nu = lowest-spin Pauli-filling count (Paper 16 eq:nu_min, corrected
    # 2026-09-28): E(ceil(N/2)) + E(floor(N/2)), degree d holding (d+1)(d+2)/2
    # functions. Values computed BY HAND, not from the code. The retracted
    # 'universal nu = N-2' agrees only through N = 8; N = 9, 10 reject it.
    @pytest.mark.parametrize("Z, nu", [
        (11, 12), (12, 14), (13, 16), (14, 18),
        (15, 20), (16, 22), (17, 24), (18, 26),
    ])
    def test_nu(self, Z: int, nu: int) -> None:
        assert classify_atom(Z).nu == nu

    # mu_free = nu(nu+3N-2)/2 (Paper 16 eq:mu_formula); N=Z.
    @pytest.mark.parametrize("Z, mu_free", [
        (11, 258.0), (12, 336.0), (13, 424.0), (14, 522.0),
        (15, 630.0), (16, 748.0), (17, 876.0), (18, 1014.0),
    ])
    def test_mu_free(self, Z: int, mu_free: float) -> None:
        assert classify_atom(Z).mu_free == pytest.approx(mu_free)


# ---------------------------------------------------------------------------
# PK parameters — frozen core
# ---------------------------------------------------------------------------

class TestPKParamsRow2:
    """All Z=11-18 use frozen_core PK with pk_params=None."""

    @pytest.mark.parametrize("Z", range(11, 19))
    def test_pk_source_frozen_core(self, Z: int) -> None:
        c = classify_atom(Z)
        assert c.pk_source == 'frozen_core', f"Z={Z}: pk_source = {c.pk_source}"

    @pytest.mark.parametrize("Z", range(11, 19))
    def test_pk_params_none(self, Z: int) -> None:
        c = classify_atom(Z)
        assert c.pk_params is None, f"Z={Z}: pk_params should be None, got {c.pk_params}"


# ---------------------------------------------------------------------------
# Support status
# ---------------------------------------------------------------------------

class TestSupportRow2:
    """Z=11-18 are supported; Z>18 is not."""

    @pytest.mark.parametrize("Z", range(11, 19))
    def test_supported(self, Z: int) -> None:
        c = classify_atom(Z)
        assert c.supported is True, f"Z={Z} should be supported"
        assert c.support_note == ''

    @pytest.mark.parametrize("Z", [50])
    def test_unsupported_beyond_row3(self, Z: int) -> None:
        c = classify_atom(Z)
        assert c.supported is False
        assert c.support_note != ''
        assert c.pk_params is None

    @pytest.mark.parametrize("Z", [19, 20])
    def test_row3_s_block_supported(self, Z: int) -> None:
        """Z=19 (K) and Z=20 (Ca) are supported via [Ar] frozen core."""
        c = classify_atom(Z)
        assert c.supported is True

    @pytest.mark.parametrize("Z", [21, 30])
    def test_transition_metals_supported(self, Z: int) -> None:
        """Z=21-30 are now supported (CLAUDE.md v2.8.0 TM hydrides extension)."""
        c = classify_atom(Z)
        assert c.supported is True

    def test_unsupported_has_correct_N(self) -> None:
        """Z=19 returns correct N_electrons and nu."""
        c = classify_atom(19)
        assert c.N_electrons == 19
        assert c.nu == 28  # E(10)+E(9) = 15+13 (retracted N-2 gave 17)
        assert c.mu_free == pytest.approx(1162.0)  # 28*83/2


# ---------------------------------------------------------------------------
# Dataclass integrity
# ---------------------------------------------------------------------------

class TestDataclassRow2:
    """AtomClassification fields for second-row atoms."""

    def test_na_full_classification(self) -> None:
        """Na (Z=11) smoke test for all fields."""
        c = classify_atom(11)
        assert isinstance(c, AtomClassification)
        assert c.Z == 11
        assert c.N_electrons == 11
        assert c.structure_type == 'C'
        assert c.n_core_electrons == 10
        assert c.n_valence_electrons == 1
        assert c.Z_eff_valence == pytest.approx(1.0)
        assert c.nu == 12   # E(6)+E(5) = 7+5
        assert c.mu_free == pytest.approx(258.0)  # 12*43/2
        assert c.pk_params is None
        assert c.pk_source == 'frozen_core'
        assert c.core_config == '1s2 2s2 2p6'
        assert c.valence_config == '3s1'
        assert c.period == 3
        assert c.group_type == 'alkali_metal'
        assert c.supported is True
        assert c.support_note == ''

    def test_ar_full_classification(self) -> None:
        """Ar (Z=18) smoke test for all fields."""
        c = classify_atom(18)
        assert isinstance(c, AtomClassification)
        assert c.Z == 18
        assert c.N_electrons == 18
        assert c.structure_type == 'B'
        assert c.n_core_electrons == 10
        assert c.n_valence_electrons == 8
        assert c.Z_eff_valence == pytest.approx(8.0)
        assert c.nu == 26   # E(9)+E(9) = 13+13
        assert c.mu_free == pytest.approx(1014.0)  # 26*78/2
        assert c.pk_params is None
        assert c.pk_source == 'frozen_core'
        assert c.core_config == '1s2 2s2 2p6'
        assert c.valence_config == '3s2 3p6'
        assert c.period == 3
        assert c.group_type == 'noble_gas'
        assert c.supported is True


# ---------------------------------------------------------------------------
# Heavy rows of Paper 16 tab:mu_values (8th cert backfill, 2026-07-02):
# Kr/Xe/Rn/Og ride the classifier's generic fallback formula path (Xe/Rn are
# not in _ATOM_DATA, supported=False) and had zero test coverage anywhere.
# ---------------------------------------------------------------------------

class TestHeavyRowFormulaPath:
    """Lowest-spin Pauli-filling nu and mu_free on the fallback path, pinned to
    the corrected Paper 16 tab:mu_values heavy rows (retracted N-2 values:
    34/52/84/116)."""

    @pytest.mark.parametrize("Z, nu, mu_free", [
        (36, 78, 7176.0),      # Kr: 2*E(18) = 2*39
        (54, 146, 22338.0),    # Xe (fallback): 2*E(27) = 2*73
        (86, 290, 79170.0),    # Rn (fallback): 2*E(43) = 2*145
        (118, 456, 184224.0),  # Og: 2*E(59) = 2*228
    ])
    def test_heavy_row(self, Z: int, nu: int, mu_free: float) -> None:
        c = classify_atom(Z)
        assert c.nu == nu
        assert c.mu_free == pytest.approx(mu_free)
