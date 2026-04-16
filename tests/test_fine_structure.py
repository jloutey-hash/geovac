"""Tests for geovac.fine_structure (Track T8 deliverable).

Verifies the Darwin and mass-velocity corrections, the combined
α⁴ fine-structure formula, and the Dirac accidental degeneracy
(2s_{1/2} = 2p_{1/2}).

All verifications are symbolic (sympy exact equality), not numerical.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Integer, Rational, simplify

from geovac.dirac_matrix_elements import (
    DiracLabel,
    alpha_sym,
    Z_sym,
    kappa_to_j,
    kappa_to_l,
)
from geovac.spin_orbit import so_diagonal_matrix_element
from geovac.fine_structure import (
    darwin_diagonal,
    mass_velocity_diagonal,
    fine_structure_total,
    verify_dirac_formula,
)


# ---------------------------------------------------------------------------
# Darwin term
# ---------------------------------------------------------------------------


class TestDarwin:
    """Darwin term is nonzero only for l=0 (s-states)."""

    def test_darwin_l0_n1(self):
        """1s (n=1, κ=-1, l=0): E_D = Z⁴α²/2."""
        val = darwin_diagonal(1, -1, Z=1)
        expected = alpha_sym**2 / Integer(2)
        assert simplify(val - expected) == 0

    def test_darwin_l0_n2(self):
        """2s (n=2, κ=-1, l=0): E_D = Z⁴α²/16."""
        val = darwin_diagonal(2, -1, Z=1)
        expected = alpha_sym**2 / Integer(16)
        assert simplify(val - expected) == 0

    def test_darwin_l0_n3(self):
        """3s: E_D = Z⁴α²/(2·27) = Z⁴α²/54."""
        val = darwin_diagonal(3, -1, Z=1)
        expected = alpha_sym**2 / Integer(54)
        assert simplify(val - expected) == 0

    def test_darwin_zero_for_l1(self):
        """2p (n=2, l=1, κ=-2): Darwin must be exactly 0."""
        val = darwin_diagonal(2, -2, Z=1)
        assert val == 0

    def test_darwin_zero_for_l1_other_j(self):
        """2p_{1/2} (n=2, l=1, κ=+1): Darwin must be exactly 0."""
        val = darwin_diagonal(2, +1, Z=1)
        assert val == 0

    def test_darwin_zero_for_l2(self):
        """3d (n=3, l=2, κ=-3): Darwin must be exactly 0."""
        val = darwin_diagonal(3, -3, Z=1)
        assert val == 0

    def test_darwin_z4_scaling(self):
        """Darwin scales as Z⁴."""
        d1 = darwin_diagonal(2, -1, Z=Integer(1))
        d3 = darwin_diagonal(2, -1, Z=Integer(3))
        ratio = simplify(d3 / d1)
        assert ratio == Integer(81)  # 3⁴

    def test_darwin_nonrelativistic_limit(self):
        """Darwin → 0 as α → 0."""
        val = darwin_diagonal(2, -1, Z=1)
        assert val.subs(alpha_sym, 0) == 0


# ---------------------------------------------------------------------------
# Mass-velocity term
# ---------------------------------------------------------------------------


class TestMassVelocity:
    """Mass-velocity is always negative and applies to all (n, l)."""

    def test_mv_sign_negative_1s(self):
        """MV is negative for 1s."""
        val = mass_velocity_diagonal(1, -1, Z=1)
        val_num = float(val.subs(alpha_sym, Rational(1, 100)))
        assert val_num < 0

    def test_mv_sign_negative_2p(self):
        """MV is negative for 2p_{3/2}."""
        val = mass_velocity_diagonal(2, -2, Z=1)
        val_num = float(val.subs(alpha_sym, Rational(1, 100)))
        assert val_num < 0

    def test_mv_sign_negative_all_n4(self):
        """MV is negative for all states up to n=4."""
        for n in range(1, 5):
            for l in range(n):
                kappa = -(l + 1)
                val = mass_velocity_diagonal(n, kappa, Z=1)
                val_num = float(val.subs(alpha_sym, Rational(1, 100)))
                assert val_num < 0, f"MV positive at n={n}, l={l}: {val}"

    def test_mv_1s_value(self):
        """1s: E_MV = −(Z⁴α²/2) · [1/(1/2) − 3/4] = −(Z⁴α²/2) · 5/4 = −5Z⁴α²/8."""
        val = mass_velocity_diagonal(1, -1, Z=1)
        expected = -Integer(5) * alpha_sym**2 / Integer(8)
        assert simplify(val - expected) == 0

    def test_mv_2s_value(self):
        """2s: E_MV = −(Z⁴α²/32) · [2/(1/2) − 3/4] = −(Z⁴α²/32) · 13/4 = −13α²/128."""
        val = mass_velocity_diagonal(2, -1, Z=1)
        expected = -Integer(13) * alpha_sym**2 / Integer(128)
        assert simplify(val - expected) == 0

    def test_mv_2p_value(self):
        """2p (κ=-2): E_MV = −(α²/32) · [2/(3/2) − 3/4] = −(α²/32)·7/12 = −7α²/384."""
        val = mass_velocity_diagonal(2, -2, Z=1)
        expected = -Integer(7) * alpha_sym**2 / Integer(384)
        assert simplify(val - expected) == 0

    def test_mv_z4_scaling(self):
        """MV scales as Z⁴."""
        mv1 = mass_velocity_diagonal(2, -1, Z=Integer(1))
        mv4 = mass_velocity_diagonal(2, -1, Z=Integer(4))
        ratio = simplify(mv4 / mv1)
        assert ratio == Integer(256)  # 4⁴

    def test_mv_nonrelativistic_limit(self):
        """MV → 0 as α → 0."""
        val = mass_velocity_diagonal(2, -1, Z=1)
        assert val.subs(alpha_sym, 0) == 0


# ---------------------------------------------------------------------------
# Dirac formula verification: the central test
# ---------------------------------------------------------------------------


class TestDiracFormula:
    """E_SO + E_D + E_MV = −(Z⁴α²)/(2n⁴) · [n/(j+½) − ¾]
    must hold exactly for ALL (n, l, j) states.
    """

    def test_hydrogen_n_max_4(self):
        """Verify for hydrogen (Z=1) up to n=4: 16 states."""
        result = verify_dirac_formula(4, Z=Integer(1))
        assert result['all_match'], f"Failures: {result['failures']}"
        assert result['n_states'] == 16

    def test_lithium_like_n_max_3(self):
        """Verify for Z=3 up to n=3: 9 states."""
        result = verify_dirac_formula(3, Z=Integer(3))
        assert result['all_match'], f"Failures: {result['failures']}"
        assert result['n_states'] == 9

    def test_symbolic_Z_n_max_3(self):
        """Verify with symbolic Z up to n=3."""
        result = verify_dirac_formula(3, Z=Z_sym)
        assert result['all_match'], f"Failures: {result['failures']}"

    def test_beryllium_like_n_max_2(self):
        """Verify for Z=4 up to n=2."""
        result = verify_dirac_formula(2, Z=Integer(4))
        assert result['all_match']

    def test_heavy_Z38_n_max_2(self):
        """Verify for Z=38 (Sr) up to n=2."""
        result = verify_dirac_formula(2, Z=Integer(38))
        assert result['all_match']


# ---------------------------------------------------------------------------
# Dirac accidental degeneracy
# ---------------------------------------------------------------------------


class TestDiracDegeneracy:
    """The 2s_{1/2} and 2p_{1/2} have the same fine-structure energy
    despite different l values. This is the Dirac "accidental" degeneracy.
    """

    def test_2s_equals_2p_hydrogen(self):
        """2s_{1/2} (κ=-1) and 2p_{1/2} (κ=+1) have equal E_FS."""
        e_2s = fine_structure_total(2, -1, Z=1)
        e_2p12 = fine_structure_total(2, +1, Z=1)
        assert simplify(e_2s - e_2p12) == 0

    def test_2s_equals_2p_symbolic_Z(self):
        """Same degeneracy with symbolic Z."""
        e_2s = fine_structure_total(2, -1, Z=Z_sym)
        e_2p12 = fine_structure_total(2, +1, Z=Z_sym)
        assert simplify(e_2s - e_2p12) == 0

    def test_3s_equals_3p12_hydrogen(self):
        """3s_{1/2} (κ=-1) and 3p_{1/2} (κ=+1) have equal E_FS."""
        e_3s = fine_structure_total(3, -1, Z=1)
        e_3p12 = fine_structure_total(3, +1, Z=1)
        assert simplify(e_3s - e_3p12) == 0

    def test_3p32_equals_3d32_hydrogen(self):
        """3p_{3/2} (κ=-2) and 3d_{3/2} (κ=+2) have equal E_FS (same j=3/2)."""
        e_3p32 = fine_structure_total(3, -2, Z=1)
        e_3d32 = fine_structure_total(3, +2, Z=1)
        assert simplify(e_3p32 - e_3d32) == 0

    def test_4s_equals_4p12_hydrogen(self):
        """4s_{1/2} = 4p_{1/2}."""
        e_4s = fine_structure_total(4, -1, Z=1)
        e_4p12 = fine_structure_total(4, +1, Z=1)
        assert simplify(e_4s - e_4p12) == 0

    def test_4p32_equals_4d32_hydrogen(self):
        """4p_{3/2} = 4d_{3/2}."""
        e_4p32 = fine_structure_total(4, -2, Z=1)
        e_4d32 = fine_structure_total(4, +2, Z=1)
        assert simplify(e_4p32 - e_4d32) == 0

    def test_4d52_equals_4f52_hydrogen(self):
        """4d_{5/2} = 4f_{5/2}."""
        e_4d52 = fine_structure_total(4, -3, Z=1)
        e_4f52 = fine_structure_total(4, +3, Z=1)
        assert simplify(e_4d52 - e_4f52) == 0

    def test_degeneracy_is_j_dependent(self):
        """States with different j must have DIFFERENT E_FS."""
        e_2s = fine_structure_total(2, -1, Z=1)   # j=1/2
        e_2p32 = fine_structure_total(2, -2, Z=1)  # j=3/2
        diff = simplify(e_2s - e_2p32)
        assert diff != 0


# ---------------------------------------------------------------------------
# Fine-structure splitting values
# ---------------------------------------------------------------------------


class TestFineStructureSplitting:
    """Verify known hydrogen fine-structure splitting values."""

    def test_2p_splitting_alpha2_over_32(self):
        """The 2p fine-structure splitting (j=3/2 minus j=1/2) is α²/32.

        This is the SAME as the SO-only splitting. The Darwin and MV
        terms don't change this splitting because:
        - Darwin is zero for l=1.
        - MV depends on l, not j, so it's the same for both κ=-2 and κ=+1
          at l=1... wait, κ=+1 has l=1 and κ=-2 has l=1, so MV is the same
          for both and cancels in the splitting.
        """
        e_3h = fine_structure_total(2, -2, Z=1)  # j=3/2
        e_1h = fine_structure_total(2, +1, Z=1)  # j=1/2
        splitting = simplify(e_3h - e_1h)
        assert simplify(splitting - alpha_sym**2 / Integer(32)) == 0

    def test_2s_2p32_splitting(self):
        """The 2s_{1/2} − 2p_{3/2} splitting.

        Since 2s_{1/2} = 2p_{1/2} (Dirac degeneracy), this equals the
        2p splitting α²/32.
        """
        e_2s = fine_structure_total(2, -1, Z=1)
        e_2p32 = fine_structure_total(2, -2, Z=1)
        splitting = simplify(e_2p32 - e_2s)
        assert simplify(splitting - alpha_sym**2 / Integer(32)) == 0

    def test_hydrogen_1s_fine_structure(self):
        """1s has E_FS = −(α²/2)·[1/1 − 3/4] = −(α²/2)·(1/4) = −α²/8."""
        val = fine_structure_total(1, -1, Z=1)
        expected = -alpha_sym**2 / Integer(8)
        assert simplify(val - expected) == 0

    def test_3d_splitting(self):
        """3d_{5/2} − 3d_{3/2}: j=5/2 minus j=3/2 at n=3.

        E(j=5/2) = −α²/(2·81)·[3/3 − 3/4] = −α²/162 · 1/4 = −α²/648
        E(j=3/2) = −α²/(2·81)·[3/2 − 3/4] = −α²/162 · 3/4 = −3α²/648 = −α²/216
        Splitting = −α²/648 − (−α²/216) = −α²/648 + 3α²/648 = 2α²/648 = α²/324.
        """
        e_d52 = fine_structure_total(3, -3, Z=1)  # j=5/2
        e_d32 = fine_structure_total(3, +2, Z=1)   # j=3/2
        splitting = simplify(e_d52 - e_d32)
        expected = alpha_sym**2 / Integer(324)
        assert simplify(splitting - expected) == 0


# ---------------------------------------------------------------------------
# Non-relativistic limit
# ---------------------------------------------------------------------------


class TestNonRelLimit:
    """All three corrections vanish as α → 0."""

    def test_all_corrections_zero_at_alpha_zero(self):
        for n, kappa in [(1, -1), (2, -1), (2, -2), (2, 1), (3, -1), (3, -2), (3, 2)]:
            e_d = darwin_diagonal(n, kappa, Z=1)
            e_mv = mass_velocity_diagonal(n, kappa, Z=1)
            e_fs = fine_structure_total(n, kappa, Z=1)
            assert e_d.subs(alpha_sym, 0) == 0
            assert e_mv.subs(alpha_sym, 0) == 0
            assert e_fs.subs(alpha_sym, 0) == 0


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


class TestInputValidation:
    def test_darwin_kappa_zero(self):
        with pytest.raises(ValueError, match="κ = 0"):
            darwin_diagonal(1, 0, Z=1)

    def test_mv_kappa_zero(self):
        with pytest.raises(ValueError, match="κ = 0"):
            mass_velocity_diagonal(1, 0, Z=1)

    def test_darwin_n_too_small(self):
        with pytest.raises(ValueError):
            darwin_diagonal(0, -1, Z=1)

    def test_darwin_l_too_large(self):
        with pytest.raises(ValueError, match="l must be < n"):
            darwin_diagonal(1, -2, Z=1)  # l=1 >= n=1

    def test_mv_l_too_large(self):
        with pytest.raises(ValueError, match="l must be < n"):
            mass_velocity_diagonal(1, -2, Z=1)


# ---------------------------------------------------------------------------
# Regression: T2 + T1 + D1 tests must still pass
# ---------------------------------------------------------------------------


class TestRegression:
    """Verify that T2's SO formula is unchanged by T8 additions."""

    def test_so_2p32_unchanged(self):
        val = so_diagonal_matrix_element(2, -2, Z=1)
        expected = alpha_sym**2 / Integer(96)
        assert simplify(val - expected) == 0

    def test_so_2p12_unchanged(self):
        val = so_diagonal_matrix_element(2, +1, Z=1)
        expected = -alpha_sym**2 / Integer(48)
        assert simplify(val - expected) == 0

    def test_so_splitting_unchanged(self):
        top = so_diagonal_matrix_element(2, -2, Z=1)
        bot = so_diagonal_matrix_element(2, +1, Z=1)
        assert simplify(top - bot - alpha_sym**2 / Integer(32)) == 0

    def test_kramers_unchanged(self):
        assert so_diagonal_matrix_element(1, -1, Z=1) == 0
        assert so_diagonal_matrix_element(2, -1, Z=1) == 0
        assert so_diagonal_matrix_element(3, -1, Z=38) == 0
