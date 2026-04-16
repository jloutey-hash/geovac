"""
Tests for the [Kr] (Z=37-38) and [Xe] (Z=55-56) frozen-core extensions.

Sprint 3 Track HA-A+B (v2.12.0): heavy-atom cores for Sunaga 2025 comparison
(SrH, BaH). Extends the FrozenCore infrastructure previously supporting
[Ne], [Ar], [Ar]3d¹⁰ cores. All cores use Clementi-Raimondi-Reinhardt (1967)
single-zeta Hartree-Fock orbital exponents.

Validates:
  - Kr core auto-detection for Z=37, 38; density integrates to 36 within 1%
  - Xe core auto-detection for Z=55, 56; density integrates to 54 within 1%
  - Z_eff asymptotic behavior at r=0 (→ Z) and r=∞ (→ Z - n_core)
  - Z_eff monotonically decreasing
  - Non-negative density
  - NIST core energies accessible and negative
  - atomic_classifier produces structure type C (alkali) and D (alkaline earth)
    with correct core counts, period assignment, and group labels
  - Unsupported atoms (Z=39-54, 57-86) return supported=False with informative
    messages
  - No regression in Z=11-36 classifications
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.neon_core import FrozenCore, _NIST_CORE_ENERGIES
from geovac.atomic_classifier import classify_atom


# ---------------------------------------------------------------------------
# [Kr] core (36 electrons): Z=37, 38
# ---------------------------------------------------------------------------

@pytest.fixture(params=[37, 38], ids=["Z=37(Rb)", "Z=38(Sr)"])
def kr_core(request):
    fc = FrozenCore(Z=request.param)
    fc.solve()
    return fc


class TestKrCoreBasic:
    """[Kr] core construction and autodetection."""

    def test_autodetect_kr_core(self, kr_core):
        assert kr_core.core_type == 'Kr'
        assert kr_core.n_core_electrons == 36

    def test_raw_density_normalizes_to_36(self, kr_core):
        """Analytically normalized density must integrate very near 36."""
        r_fine = np.linspace(0.001, 15.0, 5000)
        d = kr_core.density(r_fine)
        total = np.trapezoid(d, r_fine)
        # FrozenCore re-normalizes post-hoc, and the hydrogenic-basis density
        # converges to 36 on the internal grid (r_max=20). The integral on a
        # r=15 truncated grid should still be within 1% of 36.
        rel_err = abs(total - 36.0) / 36.0
        assert rel_err < 0.01, (
            f"Z={kr_core.Z}: density integrates to {total:.4f}, "
            f"rel err {rel_err*100:.3f}%"
        )

    def test_density_nonnegative(self, kr_core):
        r_test = np.linspace(0.01, 15.0, 500)
        d = kr_core.density(r_test)
        assert np.all(d >= -1e-10)

    def test_zeff_near_nucleus(self, kr_core):
        """Z_eff(r→0) → Z."""
        z_near = kr_core.z_eff(0.001)
        assert abs(z_near - kr_core.Z) < 0.5

    def test_zeff_at_infinity(self, kr_core):
        """Z_eff(r→∞) → Z - 36."""
        z_far = kr_core.z_eff(100.0)
        expected = kr_core.Z - 36
        assert abs(z_far - expected) < 0.01, (
            f"Z={kr_core.Z}: z_eff(100)={z_far:.6f}, expected {expected}"
        )

    def test_zeff_monotonic(self, kr_core):
        """Z_eff(r) decreases monotonically."""
        r_test = np.linspace(0.01, 15.0, 500)
        z_vals = kr_core.z_eff(r_test)
        diffs = np.diff(z_vals)
        assert np.all(diffs <= 1e-6), (
            f"Z={kr_core.Z}: z_eff not monotonic. "
            f"Max positive diff={np.max(diffs):.2e}"
        )

    def test_energy_is_negative(self, kr_core):
        assert kr_core.energy < 0

    def test_energy_matches_table(self, kr_core):
        """Energy property matches NIST table value."""
        expected = _NIST_CORE_ENERGIES[kr_core.Z]
        assert kr_core.energy == expected


# ---------------------------------------------------------------------------
# [Xe] core (54 electrons): Z=55, 56
# ---------------------------------------------------------------------------

@pytest.fixture(params=[55, 56], ids=["Z=55(Cs)", "Z=56(Ba)"])
def xe_core(request):
    fc = FrozenCore(Z=request.param)
    fc.solve()
    return fc


class TestXeCoreBasic:
    """[Xe] core construction and autodetection."""

    def test_autodetect_xe_core(self, xe_core):
        assert xe_core.core_type == 'Xe'
        assert xe_core.n_core_electrons == 54

    def test_raw_density_normalizes_to_54(self, xe_core):
        r_fine = np.linspace(0.001, 15.0, 5000)
        d = xe_core.density(r_fine)
        total = np.trapezoid(d, r_fine)
        rel_err = abs(total - 54.0) / 54.0
        assert rel_err < 0.01, (
            f"Z={xe_core.Z}: density integrates to {total:.4f}, "
            f"rel err {rel_err*100:.3f}%"
        )

    def test_density_nonnegative(self, xe_core):
        r_test = np.linspace(0.01, 15.0, 500)
        d = xe_core.density(r_test)
        assert np.all(d >= -1e-10)

    def test_zeff_near_nucleus(self, xe_core):
        z_near = xe_core.z_eff(0.001)
        assert abs(z_near - xe_core.Z) < 0.5

    def test_zeff_at_infinity(self, xe_core):
        z_far = xe_core.z_eff(100.0)
        expected = xe_core.Z - 54
        assert abs(z_far - expected) < 0.01

    def test_zeff_monotonic(self, xe_core):
        r_test = np.linspace(0.01, 15.0, 500)
        z_vals = xe_core.z_eff(r_test)
        diffs = np.diff(z_vals)
        assert np.all(diffs <= 1e-6)

    def test_energy_matches_table(self, xe_core):
        expected = _NIST_CORE_ENERGIES[xe_core.Z]
        assert xe_core.energy == expected


# ---------------------------------------------------------------------------
# Atomic classifier: Z=37, 38, 55, 56
# ---------------------------------------------------------------------------

class TestHeavyAtomClassifier:
    """classify_atom() for Z=37, 38, 55, 56 (supported) and neighbors
    (unsupported)."""

    def test_rb_z37(self):
        c = classify_atom(37)
        assert c.supported is True
        assert c.structure_type == 'C'
        assert c.n_core_electrons == 36
        assert c.n_valence_electrons == 1
        assert c.Z_eff_valence == 1.0
        assert c.period == 5
        assert c.group_type == 'alkali_metal'
        assert c.pk_source == 'frozen_core'
        assert c.valence_config == '5s1'

    def test_sr_z38(self):
        c = classify_atom(38)
        assert c.supported is True
        assert c.structure_type == 'D'
        assert c.n_core_electrons == 36
        assert c.n_valence_electrons == 2
        assert c.Z_eff_valence == 2.0
        assert c.period == 5
        assert c.group_type == 'alkaline_earth'
        assert c.pk_source == 'frozen_core'
        assert c.valence_config == '5s2'

    def test_cs_z55(self):
        c = classify_atom(55)
        assert c.supported is True
        assert c.structure_type == 'C'
        assert c.n_core_electrons == 54
        assert c.n_valence_electrons == 1
        assert c.Z_eff_valence == 1.0
        assert c.period == 6
        assert c.group_type == 'alkali_metal'
        assert c.pk_source == 'frozen_core'
        assert c.valence_config == '6s1'

    def test_ba_z56(self):
        c = classify_atom(56)
        assert c.supported is True
        assert c.structure_type == 'D'
        assert c.n_core_electrons == 54
        assert c.n_valence_electrons == 2
        assert c.Z_eff_valence == 2.0
        assert c.period == 6
        assert c.group_type == 'alkaline_earth'
        assert c.pk_source == 'frozen_core'
        assert c.valence_config == '6s2'

    def test_rb_cs_isostructural(self):
        """Rb and Cs should share structure type, group, and valence config
        up to principal quantum number."""
        rb = classify_atom(37)
        cs = classify_atom(55)
        assert rb.structure_type == cs.structure_type == 'C'
        assert rb.group_type == cs.group_type == 'alkali_metal'
        # Both alkalis have 1 valence electron in an s orbital
        assert rb.n_valence_electrons == cs.n_valence_electrons == 1

    def test_sr_ba_isostructural(self):
        sr = classify_atom(38)
        ba = classify_atom(56)
        assert sr.structure_type == ba.structure_type == 'D'
        assert sr.group_type == ba.group_type == 'alkaline_earth'
        assert sr.n_valence_electrons == ba.n_valence_electrons == 2


class TestUnsupportedHeavyAtoms:
    """4d block, 5p block, lanthanides, 5d block, 6p block remain
    out of scope and return supported=False with informative messages."""

    def test_z39_y_unsupported(self):
        """Yttrium (4d block) is out of scope."""
        c = classify_atom(39)
        assert c.supported is False
        assert '4d block' in c.support_note

    def test_z48_cd_unsupported(self):
        """Cadmium (4d block) is out of scope."""
        c = classify_atom(48)
        assert c.supported is False

    def test_z49_in_unsupported(self):
        """Indium (5p block) is out of scope."""
        c = classify_atom(49)
        assert c.supported is False
        assert '5p block' in c.support_note

    def test_z54_xe_unsupported(self):
        """Xenon itself is the Xe core target, not a valence system."""
        c = classify_atom(54)
        assert c.supported is False

    def test_z57_la_unsupported(self):
        """Lanthanum is out of scope (lanthanide)."""
        c = classify_atom(57)
        assert c.supported is False
        assert 'Lanthanide' in c.support_note or 'lanthanide' in c.support_note.lower()

    def test_z79_au_unsupported(self):
        """Gold (5d) is out of scope."""
        c = classify_atom(79)
        assert c.supported is False

    def test_z86_rn_unsupported(self):
        """Radon (6p) is out of scope."""
        c = classify_atom(86)
        assert c.supported is False

    def test_z92_u_beyond_range(self):
        """Uranium beyond supported range."""
        c = classify_atom(92)
        assert c.supported is False
        assert 'beyond the supported range' in c.support_note


class TestNoRegressionInExistingZ:
    """First row (Z=1-10), second row (Z=11-18), and fourth row
    (Z=19-36) classifications must not change."""

    @pytest.mark.parametrize("Z,expected_type,expected_core", [
        (1, 'A', 0),
        (2, 'B', 0),
        (3, 'C', 2),
        (4, 'D', 2),
        (10, 'B', 2),
        (11, 'C', 10),
        (12, 'D', 10),
        (17, 'E', 10),
        (18, 'B', 10),
        (19, 'C', 18),
        (20, 'D', 18),
        (21, 'F', 18),
        (30, 'F', 18),
        (31, 'E', 28),
        (35, 'E', 28),
        (36, 'B', 28),
    ])
    def test_classification_unchanged(self, Z, expected_type, expected_core):
        c = classify_atom(Z)
        assert c.structure_type == expected_type, (
            f"Z={Z}: got {c.structure_type}, expected {expected_type}"
        )
        assert c.n_core_electrons == expected_core, (
            f"Z={Z}: core electrons {c.n_core_electrons}, "
            f"expected {expected_core}"
        )
        assert c.supported is True


class TestFrozenCoreForLegacyNobleGasZ:
    """Existing behavior for Ne/Ar/Ar3d10 cores should work identically."""

    def test_ne_core_still_works(self):
        fc = FrozenCore(Z=11)
        fc.solve()
        assert fc.core_type == 'Ne'
        assert fc.n_core_electrons == 10

    def test_ar_core_still_works(self):
        fc = FrozenCore(Z=19)
        fc.solve()
        assert fc.core_type == 'Ar'
        assert fc.n_core_electrons == 18

    def test_ar3d10_core_still_works(self):
        fc = FrozenCore(Z=31)
        fc.solve()
        assert fc.core_type == 'Ar3d10'
        assert fc.n_core_electrons == 28
