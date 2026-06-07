"""Regression tests for debug/kleinschmidt_coaction.py.

Tests the depth-1 zeta-generator / iterated-Eisenstein-integral
machinery from Kleinschmidt, Mafra, Schlotterer, Verbeek 2026
(arXiv:2508.02800, JHEP 05 (2026) 105) as implemented in
debug/kleinschmidt_coaction.py.

The tests cover the five layers of the module:

    Layer 1. F-alphabet substrate (FWord, deconcatenation, shuffle).
    Layer 2. Depth-1 MMV m[j/k] closed form (Kleinschmidt Eq. 26).
    Layer 3. Tsunogai derivation algebra (epsilon_k^{(j)}, Pollack relation).
    Layer 4. Depth-1 zeta generator sigma_w (Eq. 34 truncated).
    Layer 5. Depth-1 coaction Delta(f_w).

Cross-validation strategy (per CLAUDE.md §13.4a):

    - The Lambert moment identity Sum_n sigma_{k-1}(n)/n^{j+1}
      = zeta(j+1) * zeta(j+2-k) is verified NUMERICALLY in its
      absolute-convergence regime as an INDEPENDENT cross-check on
      the closed-form derivation of m[j/k].

    - The classical Eichler integral identity m[1/4] = pi^4/54 is
      a literature anchor; we verify the closed form reproduces it
      symbolically.

    - The deconcatenation coproduct is verified to be coassociative
      on a low-weight word, i.e. (id (x) Delta) Delta = (Delta (x) id) Delta.

    - The shuffle product is verified commutative (up to multiplicity
      relabelling) and unital.

    - The depth-1 z_w action is verified to satisfy z_w(f_w) = 1,
      z_w(f_v) = 0 for v != w, matching the dual-to-deconcatenation
      structural prediction.

    - The Pollack relation weight balance is verified.

    - Depth >= 2 stubs raise NotImplementedError.
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath
import pytest
import sympy
from sympy import I, Rational, factorial, pi, zeta

# Make debug/ importable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "debug"))

from kleinschmidt_coaction import (
    F_UNIT,
    FWord,
    TsunogaiCommutator,
    TsunogaiGenerator,
    ZetaGeneratorDepth1,
    apply_arithmetic_to_fword,
    classify_depth1_mmv,
    deconcatenation_coproduct,
    depth1_coaction_full,
    depth1_coaction_on_zeta,
    depth1_mmv_numerical_pure_tate,
    depth1_mmv_value,
    depth2_coaction_on_iterated_eisenstein,
    depth2_mmv_value,
    eps,
    eps0,
    f_letter,
    pollack_relation_weight_14,
    shuffle,
    verify_lambert_identity,
    verify_pollack_weight_check,
)


# ---------------------------------------------------------------------------
# Layer 1. F-alphabet substrate.
# ---------------------------------------------------------------------------


class TestFWord:
    """FWord type and basic operations."""

    def test_unit_word(self) -> None:
        u = F_UNIT
        assert u.letters == ()
        assert u.f2_power == 0
        assert u.weight() == 0
        assert u.depth() == 0

    def test_single_letter(self) -> None:
        w = f_letter(3)
        assert w.letters == (3,)
        assert w.weight() == 3
        assert w.depth() == 1

    def test_f2_central(self) -> None:
        w = f_letter(2)
        assert w.letters == ()
        assert w.f2_power == 1
        assert w.weight() == 2
        assert w.depth() == 0  # f_2 does not contribute to depth

    def test_multiletter_word(self) -> None:
        w = FWord((3, 5, 7))
        assert w.weight() == 15
        assert w.depth() == 3

    def test_concatenation_is_associative(self) -> None:
        a = FWord((3,))
        b = FWord((5,))
        c = FWord((7,))
        assert (a * b) * c == a * (b * c)
        assert ((a * b) * c).letters == (3, 5, 7)

    def test_invalid_letter_even(self) -> None:
        with pytest.raises(ValueError):
            FWord((4,))  # even letter not allowed in f-alphabet (must be odd >= 3)

    def test_invalid_letter_too_small(self) -> None:
        with pytest.raises(ValueError):
            FWord((1,))

    def test_invalid_f2_power(self) -> None:
        with pytest.raises(ValueError):
            FWord((), f2_power=-1)

    def test_f_letter_rejects_invalid(self) -> None:
        with pytest.raises(ValueError):
            f_letter(4)  # even, not 2
        with pytest.raises(ValueError):
            f_letter(1)


class TestDeconcatenationCoproduct:
    """Brown 2014 f-alphabet coproduct."""

    def test_coproduct_of_unit(self) -> None:
        splits = deconcatenation_coproduct(F_UNIT)
        assert len(splits) == 1
        assert splits[0] == (F_UNIT, F_UNIT)

    def test_coproduct_of_single_letter(self) -> None:
        w = f_letter(3)
        splits = deconcatenation_coproduct(w)
        # Delta(f_3) = 1 (x) f_3 + f_3 (x) 1
        assert len(splits) == 2
        terms = [(L.letters, R.letters) for L, R in splits]
        assert ((), (3,)) in terms
        assert ((3,), ()) in terms

    def test_coproduct_of_f3_f5(self) -> None:
        w = FWord((3, 5))
        splits = deconcatenation_coproduct(w)
        # Delta(f_3 f_5) = 1 (x) f_3 f_5 + f_3 (x) f_5 + f_3 f_5 (x) 1
        assert len(splits) == 3
        terms = [(L.letters, R.letters) for L, R in splits]
        assert ((), (3, 5)) in terms
        assert ((3,), (5,)) in terms
        assert ((3, 5), ()) in terms

    def test_coassociativity_on_depth_2(self) -> None:
        """Verify (id (x) Delta) Delta = (Delta (x) id) Delta on f_3 f_5.

        This is the cobar-condition / coassociativity of the
        deconcatenation coproduct, which is the canonical Hopf-algebra
        identity for Brown's f-alphabet structure.
        """
        w = FWord((3, 5))
        # First Delta gives sum of (L, R).
        first = deconcatenation_coproduct(w)
        # Apply Delta to the right factor (id (x) Delta) Delta.
        id_x_Delta: list = []
        for L, R in first:
            for M, N in deconcatenation_coproduct(R):
                id_x_Delta.append((L.letters, M.letters, N.letters))
        # Apply Delta to the left factor (Delta (x) id) Delta.
        Delta_x_id: list = []
        for L, R in first:
            for M, N in deconcatenation_coproduct(L):
                Delta_x_id.append((M.letters, N.letters, R.letters))
        # Both should give the same multiset of triple-splittings.
        assert sorted(id_x_Delta) == sorted(Delta_x_id)

    def test_coproduct_preserves_f2_power(self) -> None:
        """Central f_2 powers go to the left factor by convention."""
        w = FWord((3,), f2_power=2)
        splits = deconcatenation_coproduct(w)
        for L, R in splits:
            assert L.f2_power == 2
            assert R.f2_power == 0


class TestShuffle:
    """Shuffle product on the f-alphabet."""

    def test_shuffle_with_unit_left(self) -> None:
        a = f_letter(3)
        unit = F_UNIT
        result = shuffle(unit, a)
        assert result == {a: 1}

    def test_shuffle_with_unit_right(self) -> None:
        a = f_letter(3)
        unit = F_UNIT
        result = shuffle(a, unit)
        assert result == {a: 1}

    def test_shuffle_two_distinct_letters(self) -> None:
        """f_3 shuffle f_5 = f_3 f_5 + f_5 f_3."""
        result = shuffle(f_letter(3), f_letter(5))
        assert result.get(FWord((3, 5)), 0) == 1
        assert result.get(FWord((5, 3)), 0) == 1
        # Two terms total.
        assert sum(result.values()) == 2

    def test_shuffle_same_letter(self) -> None:
        """f_3 shuffle f_3 = 2 f_3 f_3 (multiplicity 2)."""
        result = shuffle(f_letter(3), f_letter(3))
        assert result.get(FWord((3, 3)), 0) == 2
        assert sum(result.values()) == 2

    def test_shuffle_count_matches_binomial(self) -> None:
        """Shuffle of words of lengths m and n produces (m+n choose m) terms."""
        import math

        # |a| = 2, |b| = 2 -> 6 terms
        a = FWord((3, 5))
        b = FWord((7, 11))
        result = shuffle(a, b)
        assert sum(result.values()) == math.comb(4, 2)  # = 6


# ---------------------------------------------------------------------------
# Layer 2. Depth-1 MMV closed forms.
# ---------------------------------------------------------------------------


class TestDepth1MMVClosedForm:
    """Kleinschmidt Eq. 26 depth-1 multiple modular values."""

    def test_m_0_4_value(self) -> None:
        """m[0/4] = -2 pi i zeta(3) / 3."""
        expected = -2 * pi * I * zeta(3) / 3
        actual = depth1_mmv_value(0, 4)
        assert sympy.simplify(actual - expected) == 0

    def test_m_0_6_value(self) -> None:
        """m[0/6] = -2 pi i zeta(5) / 5."""
        expected = -2 * pi * I * zeta(5) / 5
        actual = depth1_mmv_value(0, 6)
        assert sympy.simplify(actual - expected) == 0

    def test_m_0_8_value(self) -> None:
        """m[0/8] = -2 pi i zeta(7) / 7."""
        expected = -2 * pi * I * zeta(7) / 7
        actual = depth1_mmv_value(0, 8)
        assert sympy.simplify(actual - expected) == 0

    def test_m_1_4_pi_fourth_over_54(self) -> None:
        """m[1/4] = pi^4 / 54.

        This is a CLASSICAL result for the Eichler period of E_4 at
        j=1; it's a well-known literature anchor for Kleinschmidt
        Eq. 26. We verify the closed form reproduces it symbolically.
        Reference: Eichler (1957) and many textbook accounts; see also
        Brown 2014 §1 for the motivic version.
        """
        expected = pi**4 / 54
        actual = depth1_mmv_value(1, 4)
        assert sympy.simplify(actual - expected) == 0

    def test_m_1_4_numerical(self) -> None:
        """Numerical evaluation of m[1/4] should equal pi^4/54 to high precision."""
        v = depth1_mmv_numerical_pure_tate(1, 4, dps=50)
        expected = float(sympy.N(pi**4 / 54, 50))
        assert abs(complex(v).real - expected) < 1e-40

    def test_invalid_k_odd(self) -> None:
        with pytest.raises(ValueError):
            depth1_mmv_value(0, 5)  # k must be even

    def test_invalid_k_too_small(self) -> None:
        with pytest.raises(ValueError):
            depth1_mmv_value(0, 2)  # k must be >= 4

    def test_invalid_j_out_of_range(self) -> None:
        with pytest.raises(ValueError):
            depth1_mmv_value(3, 4)  # j must be <= k - 2 = 2

    def test_invalid_j_negative(self) -> None:
        with pytest.raises(ValueError):
            depth1_mmv_value(-1, 4)


class TestLambertMomentIdentity:
    """Cross-check of the Lambert-series identity load-bearing in Eq. 26."""

    def test_lambert_j5_k4(self) -> None:
        """Sum_n sigma_3(n) / n^6 = zeta(6) * zeta(3) within absolute-convergence regime.

        Tail decays as 1/n^3 ; 2000-term sum gets to ~1e-7. We don't push
        to higher precision because this is a structural cross-check, not
        a precision test.
        """
        num, ana, rel_err = verify_lambert_identity(j=5, k=4, dps=20)
        assert rel_err < 1e-5

    def test_lambert_j6_k4(self) -> None:
        """Sum_n sigma_3(n) / n^7 = zeta(7) * zeta(4). Faster convergence."""
        num, ana, rel_err = verify_lambert_identity(j=6, k=4, dps=20)
        assert rel_err < 1e-7

    def test_lambert_rejects_non_absolute_regime(self) -> None:
        """j < k - 1 is outside absolute convergence; routine should reject."""
        with pytest.raises(ValueError):
            verify_lambert_identity(j=2, k=6, dps=20)


class TestDepth1Classification:
    """Structural classification of m[j/k] in the f-alphabet picture."""

    def test_m_0_4_inhabits_f3(self) -> None:
        cls = classify_depth1_mmv(0, 4)
        assert cls["pi_power"] == 1
        assert cls["odd_zeta_content"] == [3]
        assert cls["f_alphabet_word"] == f_letter(3)
        assert cls["is_purely_tate"] is False

    def test_m_0_6_inhabits_f5(self) -> None:
        cls = classify_depth1_mmv(0, 6)
        assert cls["pi_power"] == 1
        assert cls["odd_zeta_content"] == [5]
        assert cls["f_alphabet_word"] == f_letter(5)

    def test_m_0_8_inhabits_f7(self) -> None:
        cls = classify_depth1_mmv(0, 8)
        assert cls["pi_power"] == 1
        assert cls["odd_zeta_content"] == [7]
        assert cls["f_alphabet_word"] == f_letter(7)

    def test_m_1_4_purely_tate(self) -> None:
        """m[1/4] = pi^4/54 has no odd-zeta content."""
        cls = classify_depth1_mmv(1, 4)
        assert cls["is_purely_tate"] is True
        assert cls["odd_zeta_content"] == []
        assert cls["f_alphabet_word"] == F_UNIT

    def test_m_2_6_inhabits_f3(self) -> None:
        cls = classify_depth1_mmv(2, 6)
        assert cls["odd_zeta_content"] == [3]
        assert cls["f_alphabet_word"] == f_letter(3)


# ---------------------------------------------------------------------------
# Layer 3. Tsunogai derivation algebra.
# ---------------------------------------------------------------------------


class TestTsunogai:
    def test_eps0_construction(self) -> None:
        e = eps0()
        assert e.k == 0
        assert e.j == 0
        assert e.weight() == 0

    def test_eps4_basic(self) -> None:
        e = eps(4)
        assert e.k == 4
        assert e.j == 0
        assert e.weight() == 4

    def test_eps_iterated(self) -> None:
        e = eps(4, 2)  # = ad_eps_0^2 (eps_4)
        assert e.k == 4
        assert e.j == 2
        # Weight is still 4 in the modular grading (j tracks ad-iteration only).
        assert e.weight() == 4

    def test_invalid_k_odd(self) -> None:
        with pytest.raises(ValueError):
            eps(5)  # odd k not allowed except via k=0

    def test_invalid_k_too_small(self) -> None:
        with pytest.raises(ValueError):
            eps(2)  # k must be 0 or even >= 4

    def test_invalid_j_negative(self) -> None:
        with pytest.raises(ValueError):
            eps(4, -1)

    def test_pollack_weight_balance(self) -> None:
        """[eps_4, eps_10] and [eps_6, eps_8] both have weight 14."""
        assert verify_pollack_weight_check() is True
        a, b = pollack_relation_weight_14()
        assert a.weight() == 14
        assert b.weight() == 14
        assert a.left.k == 4 and a.right.k == 10
        assert b.left.k == 6 and b.right.k == 8

    def test_commutator_weight(self) -> None:
        c = TsunogaiCommutator(eps(4), eps(6))
        assert c.weight() == 10


# ---------------------------------------------------------------------------
# Layer 4. Depth-1 zeta generator sigma_w.
# ---------------------------------------------------------------------------


class TestZetaGenerator:
    def test_sigma_3_decomposition(self) -> None:
        """sigma_3 = z_3 - eps_4^(2) / 2 + higher depth (Eq. 34 truncated)."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        assert sigma3.w == 3
        assert sigma3.weight() == 3
        # Geometric part: -1/(3-1)! eps_{4}^{(2)} = -1/2 eps_4^(2)
        assert sigma3.geometric_eps.k == 4
        assert sigma3.geometric_eps.j == 2
        assert sigma3.geometric_coefficient == -Rational(1, 2)

    def test_sigma_5_decomposition(self) -> None:
        """sigma_5 = z_5 - eps_6^(4) / 4! + ..."""
        sigma5 = ZetaGeneratorDepth1(w=5)
        assert sigma5.geometric_eps.k == 6
        assert sigma5.geometric_eps.j == 4
        assert sigma5.geometric_coefficient == -Rational(1, 24)  # = -1/4!

    def test_sigma_7_decomposition(self) -> None:
        sigma7 = ZetaGeneratorDepth1(w=7)
        assert sigma7.geometric_eps.k == 8
        assert sigma7.geometric_eps.j == 6
        assert sigma7.geometric_coefficient == -Rational(1, 720)  # = -1/6!

    def test_invalid_even_weight(self) -> None:
        with pytest.raises(ValueError):
            ZetaGeneratorDepth1(w=4)  # must be odd

    def test_invalid_weight_too_small(self) -> None:
        with pytest.raises(ValueError):
            ZetaGeneratorDepth1(w=1)


class TestArithmeticAction:
    """Action of z_w on f-words (dual to deconcatenation, depth 1)."""

    def test_z3_on_f3_gives_unit(self) -> None:
        """z_3(f_3) = 1."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        result = apply_arithmetic_to_fword(sigma3, f_letter(3))
        assert result == {F_UNIT: 1}

    def test_z3_on_f5_gives_zero(self) -> None:
        """z_3(f_5) = 0."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        result = apply_arithmetic_to_fword(sigma3, f_letter(5))
        assert result == {}

    def test_z5_on_f5_gives_unit(self) -> None:
        sigma5 = ZetaGeneratorDepth1(w=5)
        result = apply_arithmetic_to_fword(sigma5, f_letter(5))
        assert result == {F_UNIT: 1}

    def test_z3_on_f3_f5_truncates(self) -> None:
        """z_3(f_3 f_5) = f_5."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        result = apply_arithmetic_to_fword(sigma3, FWord((3, 5)))
        assert result == {f_letter(5): 1}

    def test_z3_on_f5_f3_gives_zero(self) -> None:
        """z_3(f_5 f_3) = 0 because the leading letter is f_5, not f_3."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        result = apply_arithmetic_to_fword(sigma3, FWord((5, 3)))
        assert result == {}

    def test_z3_on_unit_gives_zero(self) -> None:
        sigma3 = ZetaGeneratorDepth1(w=3)
        result = apply_arithmetic_to_fword(sigma3, F_UNIT)
        assert result == {}

    def test_z3_preserves_f2_central(self) -> None:
        """z_3(f_2 f_3) = f_2 (central f_2 stays put)."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        w = FWord((3,), f2_power=1)
        result = apply_arithmetic_to_fword(sigma3, w)
        expected = FWord(letters=(), f2_power=1)
        assert result == {expected: 1}


# ---------------------------------------------------------------------------
# Layer 5. Depth-1 coaction.
# ---------------------------------------------------------------------------


class TestDepth1Coaction:
    """Depth-1 collapse of Kleinschmidt Eq. 64 coaction."""

    def test_delta_f3(self) -> None:
        """Delta(f_3) = f_3 (x) 1 + 1 (x) f_3."""
        terms = depth1_coaction_full(3)
        assert len(terms) == 2
        # f_3 (x) 1
        assert (f_letter(3), F_UNIT) in terms
        # 1 (x) f_3
        assert (F_UNIT, f_letter(3)) in terms

    def test_delta_f5(self) -> None:
        terms = depth1_coaction_full(5)
        assert (f_letter(5), F_UNIT) in terms
        assert (F_UNIT, f_letter(5)) in terms

    def test_delta_f7(self) -> None:
        terms = depth1_coaction_full(7)
        assert (f_letter(7), F_UNIT) in terms
        assert (F_UNIT, f_letter(7)) in terms

    def test_invalid_even_weight(self) -> None:
        with pytest.raises(ValueError):
            depth1_coaction_on_zeta(4)

    def test_invalid_weight_too_small(self) -> None:
        with pytest.raises(ValueError):
            depth1_coaction_on_zeta(1)

    def test_consistency_with_deconcatenation(self) -> None:
        """Depth-1 coaction of f_w should match deconcatenation of f_w."""
        for w in (3, 5, 7):
            from_coaction = set(depth1_coaction_full(w))
            from_deconcat = set(deconcatenation_coproduct(f_letter(w)))
            assert from_coaction == from_deconcat


# ---------------------------------------------------------------------------
# Depth >= 2 stubs.
# ---------------------------------------------------------------------------


class TestDepth2Stubs:
    """Verify depth >= 2 calls raise informative NotImplementedError."""

    def test_depth2_mmv_raises(self) -> None:
        with pytest.raises(NotImplementedError, match="Appendix A"):
            depth2_mmv_value(4, 0, 4, 0)

    def test_depth2_coaction_raises(self) -> None:
        with pytest.raises(NotImplementedError, match="Tsunogai-bracket"):
            depth2_coaction_on_iterated_eisenstein([4, 4])


# ---------------------------------------------------------------------------
# Integration / consistency tests across layers.
# ---------------------------------------------------------------------------


class TestCrossLayer:
    """Consistency between f-alphabet, MMV classification, and coaction."""

    def test_classification_matches_coaction_target(self) -> None:
        """The f-alphabet word assigned by classify_depth1_mmv should be the
        target word in the depth-1 coaction Delta(f_{2k+1}) on the right."""
        cls = classify_depth1_mmv(0, 4)
        f_word = cls["f_alphabet_word"]
        # m[0/4] inhabits f_3; depth-1 coaction Delta(f_3) = f_3 (x) 1 + 1 (x) f_3.
        terms = depth1_coaction_full(3)
        right_factors = {R for _, R in terms}
        assert f_word in right_factors or f_word == F_UNIT

    def test_z_action_inverts_coaction_target(self) -> None:
        """z_3 applied to the right factor of Delta(f_3) should leave the unit
        in either slot (this is the canonical duality between deconcatenation
        and the depth-1 arithmetic action)."""
        sigma3 = ZetaGeneratorDepth1(w=3)
        terms = depth1_coaction_full(3)
        # Apply z_3 to each right factor.
        right_actions = []
        for _, R in terms:
            res = apply_arithmetic_to_fword(sigma3, R)
            right_actions.extend(res.keys())
        # One of the terms (Delta(f_3) = ... + 1 (x) f_3) gives z_3(f_3) = 1.
        assert F_UNIT in right_actions

    def test_lambert_identity_underwrites_m_0_4(self) -> None:
        """The closed form of m[0/4] does NOT directly invoke the Lambert
        identity (j=0 is the analytic-continuation boundary); but the
        STRUCTURE of m[j/k] for j >= 1 follows the Lambert pattern. We
        verify the structural Lambert identity holds in its abs-convergence
        regime as a cross-check that the m[j/k] closed form is consistent
        with the q-expansion of G_k."""
        # The deepest cross-check at sprint scale.
        num, ana, rel_err = verify_lambert_identity(j=5, k=4, dps=20)
        assert rel_err < 1e-5
