"""Tests for L3b-2 sub-sprint B: joint Schur multiplier cb-norm.

Verifies the closed-form bound

    ||S_K^joint||_cb = ||S_K^SU(2)||_cb * ||S_K^U(1)||_cb = 2/(n_max + 1) * 1
                    = 2/(n_max + 1)

at panel cells (n_max, N_t) in {(2,3), (2,5), (3,3), (3,5), (4,3), (4,5)}.

Eight test functions covering:

  1. SU(2) factor matches Paper 38 cb-norm bit-identically.
  2. U(1) factor cb-norm = 1 (the closed-form constant).
  3. Joint cb-norm = product of factor cb-norms (factorization).
  4. Panel cells against closed-form bound 2/(n_max + 1).
  5. Asymptotic rate 2/(n_max + 1) = O(1/n_max) verified monotonically.
  6. U(1) cb-norm cross-check via Plancherel symbol max-norm.
  7. Joint Plancherel symbol factorization exact in sympy.
  8. The joint cb-norm depends only on n_max, not on N_t.

Companion to debug/l3b_2_sub_sprint_B_cb_norm_memo.md and
debug/l3b_2_sub_sprint_B_compute.py.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Rational

from geovac.central_fejer_su2 import (
    central_multiplier_cb_norm as cb_su2,
    plancherel_symbol as plancherel_su2,
)
from geovac.central_fejer_compact_temporal import (
    cb_norm_circle,
    joint_cb_norm,
    joint_plancherel_symbol,
    plancherel_symbol_circle,
    verify_plancherel_factorization,
)


# Panel cells from the sub-sprint B dispatch specification.
PANEL_CELLS = [(2, 3), (2, 5), (3, 3), (3, 5), (4, 3), (4, 5)]


# ---------------------------------------------------------------------------
# Test 1: SU(2) factor matches Paper 38 cb-norm bit-identically
# ---------------------------------------------------------------------------


class TestSU2FactorPaper38Consistency:
    """SU(2) cb-norm matches Paper 38 Lemma L2(c) closed form 2/(n_max+1)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4, 5, 6, 8, 10, 20, 50])
    def test_su2_cb_norm_closed_form(self, n_max):
        """||S_K^SU(2)||_cb = 2/(n_max + 1) exactly."""
        result = cb_su2(n_max)
        expected = Rational(2, n_max + 1)
        assert result == expected, (
            f"SU(2) cb-norm at n_max={n_max}: got {result}, expected {expected}"
        )

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_su2_cb_norm_attained_at_jmax(self, n_max):
        """The cb-norm equals Plancherel symbol at j = j_max."""
        j_max = Rational(n_max - 1, 2)
        symbol_at_jmax = plancherel_su2(n_max, j_max)
        cb = cb_su2(n_max)
        assert symbol_at_jmax == cb, (
            f"At n_max={n_max}: Plancherel(j_max)={symbol_at_jmax} != cb={cb}"
        )

    def test_su2_cb_norm_product_eq_2(self):
        """(n_max + 1) * cb_norm = 2 exactly at every n_max."""
        for n_max in range(1, 51):
            product = (n_max + 1) * cb_su2(n_max)
            assert product == Rational(2), (
                f"At n_max={n_max}: (n_max+1)*cb = {product} != 2"
            )


# ---------------------------------------------------------------------------
# Test 2: U(1) factor cb-norm = 1
# ---------------------------------------------------------------------------


class TestU1FactorCBNorm:
    """U(1) Fejér kernel cb-norm equals 1 at every N_t."""

    @pytest.mark.parametrize("N_t", [1, 2, 3, 4, 5, 7, 11, 21, 51])
    def test_u1_cb_norm_equals_one(self, N_t):
        """||S_K^U(1)||_cb = 1 exactly."""
        result = cb_norm_circle(N_t)
        assert result == Rational(1), (
            f"U(1) cb-norm at N_t={N_t}: got {result}, expected 1"
        )

    @pytest.mark.parametrize("N_t", [1, 2, 3, 4, 5])
    def test_u1_cb_norm_is_sympy_rational(self, N_t):
        """cb_norm_circle returns sympy.Rational, not float."""
        result = cb_norm_circle(N_t)
        assert isinstance(result, sp.Rational), (
            f"At N_t={N_t}: cb_norm_circle returned {type(result).__name__}, "
            f"expected sympy.Rational"
        )

    def test_u1_cb_norm_rejects_invalid(self):
        """cb_norm_circle raises ValueError for N_t < 1."""
        with pytest.raises(ValueError, match="N_t must be >= 1"):
            cb_norm_circle(0)
        with pytest.raises(ValueError, match="N_t must be >= 1"):
            cb_norm_circle(-1)


# ---------------------------------------------------------------------------
# Test 3: Joint cb-norm = product of factor cb-norms (factorization)
# ---------------------------------------------------------------------------


class TestJointCBNormFactorization:
    """Joint cb-norm factorizes as product of SU(2) and U(1) factors."""

    @pytest.mark.parametrize("n_max,N_t", PANEL_CELLS)
    def test_joint_eq_product_of_factors(self, n_max, N_t):
        """||S_K^joint||_cb = ||S_K^SU(2)||_cb * ||S_K^U(1)||_cb."""
        joint = joint_cb_norm(n_max, N_t)
        su2 = cb_su2(n_max)
        u1 = cb_norm_circle(N_t)
        product = su2 * u1
        assert joint == product, (
            f"At (n_max={n_max}, N_t={N_t}): joint={joint} != "
            f"su2*u1 = {su2}*{u1} = {product}"
        )

    @pytest.mark.parametrize("n_max,N_t", PANEL_CELLS)
    def test_joint_is_sympy_rational(self, n_max, N_t):
        """joint_cb_norm returns sympy.Rational."""
        result = joint_cb_norm(n_max, N_t)
        assert isinstance(result, sp.Rational), (
            f"At (n_max={n_max}, N_t={N_t}): joint_cb_norm returned "
            f"{type(result).__name__}, expected sympy.Rational"
        )


# ---------------------------------------------------------------------------
# Test 4: Panel cells against closed-form bound 2/(n_max + 1)
# ---------------------------------------------------------------------------


class TestPanelCellsAgainstBound:
    """At each panel cell, joint cb-norm equals 2/(n_max + 1) exactly."""

    @pytest.mark.parametrize("n_max,N_t", PANEL_CELLS)
    def test_joint_eq_bound(self, n_max, N_t):
        """||S_K^joint||_cb = 2/(n_max + 1)."""
        result = joint_cb_norm(n_max, N_t)
        expected = Rational(2, n_max + 1)
        assert result == expected, (
            f"At (n_max={n_max}, N_t={N_t}): joint={result} != "
            f"2/({n_max}+1) = {expected}"
        )

    def test_panel_specific_values(self):
        """Hard-coded panel values match the memo's Table in §8."""
        expected_values = {
            (2, 3): Rational(2, 3),
            (2, 5): Rational(2, 3),
            (3, 3): Rational(1, 2),
            (3, 5): Rational(1, 2),
            (4, 3): Rational(2, 5),
            (4, 5): Rational(2, 5),
        }
        for (n_max, N_t), expected in expected_values.items():
            result = joint_cb_norm(n_max, N_t)
            assert result == expected, (
                f"At (n_max={n_max}, N_t={N_t}): got {result}, expected {expected}"
            )


# ---------------------------------------------------------------------------
# Test 5: Asymptotic rate 2/(n_max + 1) = O(1/n_max) verified monotonically
# ---------------------------------------------------------------------------


class TestAsymptoticRate:
    """Joint cb-norm decreases monotonically to 0 as n_max -> infinity."""

    def test_monotone_decreasing(self):
        """At fixed N_t, joint_cb_norm is strictly decreasing in n_max."""
        N_t = 5
        prev = None
        for n_max in [1, 2, 3, 4, 5, 10, 20, 50, 100]:
            cur = joint_cb_norm(n_max, N_t)
            if prev is not None:
                assert cur < prev, (
                    f"Monotonicity violated at n_max={n_max}, N_t={N_t}: "
                    f"prev={prev}, cur={cur}"
                )
            prev = cur

    def test_rate_constant_2(self):
        """(n_max + 1) * joint_cb_norm = 2 exactly across all n_max and N_t."""
        for n_max in [1, 2, 3, 5, 10, 50, 100]:
            for N_t in [1, 3, 5, 11]:
                product = (n_max + 1) * joint_cb_norm(n_max, N_t)
                assert product == Rational(2), (
                    f"At (n_max={n_max}, N_t={N_t}): "
                    f"(n_max+1)*joint = {product} != 2"
                )

    def test_asymptote_zero(self):
        """joint_cb_norm -> 0 as n_max -> infinity."""
        # At n_max = 10000, the cb-norm should be very small
        val = joint_cb_norm(10000, 5)
        assert float(val) < 1e-3, (
            f"At n_max=10000: cb={val}={float(val):.6f}, expected < 1e-3"
        )


# ---------------------------------------------------------------------------
# Test 6: U(1) cb-norm cross-check via Plancherel symbol max-norm
# ---------------------------------------------------------------------------


class TestU1CBNormCrossCheckViaMaxNorm:
    """Verify U(1) cb-norm = 1 by direct computation of max symbol."""

    @pytest.mark.parametrize("N_t", [1, 2, 3, 5, 7, 11, 21])
    def test_max_symbol_attained_at_zero(self, N_t):
        """max_k hat{K}^U(1)(k) is attained at k = 0."""
        K_max = (N_t - 1) // 2 if N_t % 2 == 1 else N_t // 2
        k_vals = list(range(-K_max, K_max + 1))
        symbols = [(k, plancherel_symbol_circle(N_t, k)) for k in k_vals]
        max_k, max_val = max(symbols, key=lambda kv: kv[1])
        assert max_k == 0, (
            f"At N_t={N_t}: max attained at k={max_k}, expected k=0"
        )

    @pytest.mark.parametrize("N_t", [1, 2, 3, 5, 7, 11, 21])
    def test_max_symbol_equals_one(self, N_t):
        """max_k hat{K}^U(1)(k) = hat{K}(0) = 1."""
        symbol_at_zero = plancherel_symbol_circle(N_t, 0)
        assert symbol_at_zero == Rational(1), (
            f"At N_t={N_t}: hat{{K}}(0)={symbol_at_zero}, expected 1"
        )

    @pytest.mark.parametrize("N_t", [3, 5, 7])
    def test_symbol_decreasing_from_zero(self, N_t):
        """Plancherel symbol is decreasing as |k| increases from 0."""
        K_max = (N_t - 1) // 2 if N_t % 2 == 1 else N_t // 2
        # k = 0 should be larger than k = 1
        s0 = plancherel_symbol_circle(N_t, 0)
        s1 = plancherel_symbol_circle(N_t, 1)
        assert s0 > s1, (
            f"At N_t={N_t}: hat{{K}}(0)={s0} should be > hat{{K}}(1)={s1}"
        )


# ---------------------------------------------------------------------------
# Test 7: Joint Plancherel symbol factorization exact in sympy
# ---------------------------------------------------------------------------


class TestPlancherelFactorizationExact:
    """Joint Plancherel symbol factorizes EXACTLY in sympy rationals."""

    @pytest.mark.parametrize("n_max,N_t", PANEL_CELLS)
    def test_factorization_exact_per_cell(self, n_max, N_t):
        """hat{K}^joint(j, k) = hat{K}^SU(2)(j) * hat{K}^U(1)(k) exactly."""
        ok, details = verify_plancherel_factorization(n_max, N_t)
        assert ok, (
            f"At (n_max={n_max}, N_t={N_t}): factorization failed; "
            f"{details['pairs_match']}/{details['pairs_checked']} pairs match; "
            f"failures: {details['failures'][:3]}..."
        )
        assert details["pairs_match"] == details["pairs_checked"]
        assert len(details["failures"]) == 0

    def test_factorization_individual_value(self):
        """At (n_max=3, N_t=5), check a specific symbol value explicitly."""
        # j = 1, k = 0: should give (2*1+1)/Z_su2 * 1 = 3/6 * 1 = 1/2
        n_max, N_t = 3, 5
        j = Rational(1)
        k = 0
        joint = joint_plancherel_symbol(n_max, N_t, j, k)
        su2 = plancherel_su2(n_max, j)
        u1 = plancherel_symbol_circle(N_t, k)
        assert joint == su2 * u1
        assert joint == Rational(1, 2)


# ---------------------------------------------------------------------------
# Test 8: Joint cb-norm depends only on n_max, not on N_t
# ---------------------------------------------------------------------------


class TestJointCBNormNtIndependence:
    """At fixed n_max, joint cb-norm is constant in N_t (it equals 2/(n_max+1))."""

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_joint_cb_norm_independent_of_Nt(self, n_max):
        """joint_cb_norm(n_max, N_t) is the same for all N_t."""
        ref = joint_cb_norm(n_max, 3)
        for N_t in [1, 2, 3, 5, 7, 11, 21, 51]:
            cur = joint_cb_norm(n_max, N_t)
            assert cur == ref, (
                f"At n_max={n_max}: joint_cb_norm depends on N_t? "
                f"N_t=3 gives {ref}, N_t={N_t} gives {cur}"
            )

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_joint_cb_norm_equals_su2_cb_norm(self, n_max):
        """Since u1_cb = 1, joint = su2 * 1 = su2 at all N_t."""
        su2 = cb_su2(n_max)
        for N_t in [1, 3, 5, 11]:
            joint = joint_cb_norm(n_max, N_t)
            assert joint == su2, (
                f"At (n_max={n_max}, N_t={N_t}): joint={joint} != su2_cb={su2}"
            )


# ---------------------------------------------------------------------------
# Bonus: structural sanity (positivity of Plancherel symbols)
# ---------------------------------------------------------------------------


class TestPlancherelSymbolsArePositive:
    """Both factor Plancherel symbols are non-negative (needed for B3 factorization)."""

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_su2_plancherel_nonneg(self, n_max):
        """hat{K}^SU(2)(j) >= 0 for all j."""
        for k in range(0, 2 * n_max):
            j = Rational(k, 2)
            val = plancherel_su2(n_max, j)
            assert val >= 0, (
                f"hat{{K}}^SU(2)({j}) = {val} at n_max={n_max} is negative!"
            )

    @pytest.mark.parametrize("N_t", [1, 3, 5, 7, 11])
    def test_u1_plancherel_nonneg(self, N_t):
        """hat{K}^U(1)(k) >= 0 for all k."""
        K_max = (N_t - 1) // 2 if N_t % 2 == 1 else N_t // 2
        for k in range(-K_max - 2, K_max + 3):
            val = plancherel_symbol_circle(N_t, k)
            assert val >= 0, (
                f"hat{{K}}^U(1)({k}) = {val} at N_t={N_t} is negative!"
            )
