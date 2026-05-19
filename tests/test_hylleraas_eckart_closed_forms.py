"""Tests for geovac/hylleraas_eckart_closed_forms.py.

Track 1 of the Hylleraas-Eckart double-alpha extension sprint (CLAUDE.md
§2 backlog entry, v2.48.0). Verifies the closed-form Hylleraas-Eckart
master integral

    I_HE^cosh(L, M, N; alpha, B)
        = integral_0^infty ds exp(-2 alpha s) s^L
          * integral_0^s du u^{N+1}
          * integral_{-u}^u dt t^{2M} (s^2 - t^2) cosh(B t)

against the existing single-alpha Hylleraas master integral (B->0 limit),
the scoping prototype, and the explicit (0,0,0) closed form
64 / (4 alpha^2 - B^2)^3.
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile
from fractions import Fraction
from typing import List, Tuple

import pytest

from geovac.hylleraas_eckart_closed_forms import (
    ClosedForm,
    DEFAULT_CACHE_PATH,
    I_HE_cosh_b_zero_check,
    I_HE_cosh_numeric,
    I_HE_cosh_polynomial,
    I_HE_cosh_symbolic,
    SCHEMA_VERSION,
    _CACHE_TABLE,
    load_table,
    precompute_table,
    save_table,
)
from geovac.hylleraas_r12 import hylleraas_master_int


# ---------------------------------------------------------------------------
# Closed-form sanity at (0, 0, 0): the textbook case
# ---------------------------------------------------------------------------

class TestClosedFormSanity:
    """The (0, 0, 0) case has known closed form 64 / (4 alpha^2 - B^2)^3."""

    def test_zero_zero_zero_at_b_eq_zero(self):
        """I_HE^cosh(0,0,0; alpha=1, B=0) = 1 (matches single-alpha I=1/alpha^6)."""
        val = I_HE_cosh_numeric(0, 0, 0, 1.0, 0.0)
        assert abs(val - 1.0) < 1e-14

    def test_zero_zero_zero_canonical_alpha(self):
        """I_HE^cosh(0,0,0; alpha=1.6875, B=0) = 1/1.6875^6 (single-alpha value)."""
        val = I_HE_cosh_numeric(0, 0, 0, 1.6875, 0.0)
        expected = 1.0 / 1.6875**6
        assert abs(val - expected) < 1e-14

    def test_zero_zero_zero_at_b_eq_one(self):
        """I_HE^cosh(0,0,0; alpha=1.5, B=1) = 64/(4*1.5^2 - 1)^3 = 64/8^3 = 0.125."""
        val = I_HE_cosh_numeric(0, 0, 0, 1.5, 1.0)
        assert abs(val - 0.125) < 1e-14

    def test_zero_zero_zero_general_closed_form(self):
        """I_HE^cosh(0,0,0; alpha, B) = 64/(4*alpha^2 - B^2)^3 across panel."""
        for alpha, B in [(2.0, 0.5), (1.5, 0.8), (1.0, 0.5)]:
            val = I_HE_cosh_numeric(0, 0, 0, alpha, B)
            expected = 64.0 / (4 * alpha**2 - B**2) ** 3
            assert abs(val - expected) < 1e-14, (
                f"alpha={alpha}, B={B}: got {val}, expected {expected}"
            )


# ---------------------------------------------------------------------------
# B -> 0 reduction to single-alpha master integral
# ---------------------------------------------------------------------------

class TestBZeroReduction:
    """At B=0, the Hylleraas-Eckart master must reduce bit-exactly to the
    single-alpha Hylleraas master from geovac/hylleraas_r12.py.

    This is the load-bearing regression invariant: every Eckart computation
    at beta=0 must recover the existing single-alpha pipeline.
    """

    PROTOTYPE_CASES = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1)]

    @pytest.mark.parametrize("L,M,N", PROTOTYPE_CASES)
    def test_symbolic_b_zero_reduction(self, L: int, M: int, N: int):
        """Sympy: I_HE^cosh(L,M,N; alpha, 0) - I_single(L,M,N; alpha) = 0 symbolically."""
        diff = I_HE_cosh_b_zero_check(L, M, N)
        assert diff == 0, f"({L},{M},{N}): symbolic diff = {diff}"

    def test_panel_b_zero_numerical(self):
        """Numerical: B=0 panel for L+2M+N <= 5 matches single-alpha at machine precision."""
        alpha = 1.6875
        max_total = 5
        worst_rel = 0.0
        n_cells = 0
        for L in range(max_total + 1):
            for M in range((max_total - L) // 2 + 1):
                for N in range(max_total - L - 2 * M + 1):
                    single = hylleraas_master_int(L, M, N, alpha)
                    eckart = I_HE_cosh_numeric(L, M, N, alpha, 0.0)
                    rel = abs(single - eckart) / abs(single)
                    worst_rel = max(worst_rel, rel)
                    n_cells += 1
        assert worst_rel < 1e-12, (
            f"Worst rel_err {worst_rel:.2e} across {n_cells} cells exceeds 1e-12"
        )

    @pytest.mark.parametrize("L,M,N", PROTOTYPE_CASES)
    @pytest.mark.parametrize("alpha", [1.0, 1.5, 1.6875, 2.0, 3.0])
    def test_b_zero_at_various_alpha(self, L: int, M: int, N: int, alpha: float):
        """B=0 reduction holds at multiple alpha values."""
        single = hylleraas_master_int(L, M, N, alpha)
        eckart = I_HE_cosh_numeric(L, M, N, alpha, 0.0)
        rel = abs(single - eckart) / abs(single)
        assert rel < 1e-13, (
            f"({L},{M},{N}) @ alpha={alpha}: rel_err={rel:.2e}"
        )


# ---------------------------------------------------------------------------
# Polynomial encoding invariants
# ---------------------------------------------------------------------------

class TestPolynomialEncoding:
    """The ClosedForm encoding stores polynomials in (a := alpha, b := B^2)
    with integer-rational coefficients. By t-parity of the integrand, only
    even B-powers appear, which is enforced in _sympy_to_closed_form."""

    def test_zero_zero_zero_evaluates_correctly(self):
        """(0,0,0) closed form evaluates to 64/(4 alpha^2 - B^2)^3 at multiple points.

        We test numerical equivalence rather than exact polynomial encoding
        because the recurrence engine and the sympy oracle can produce
        DIFFERENT but EQUIVALENT polynomial reductions of the same rational
        function (e.g., 64/(4a^2-B^2)^3 vs (192 a^2 - ...) / ((4a^2-B^2)^3 * 2a)).
        Both are valid; only the numerical value must agree.
        """
        cf = I_HE_cosh_polynomial(0, 0, 0)
        for alpha, B in [(1.0, 0.0), (1.5, 0.5), (2.0, 1.0), (1.6875, 0.7)]:
            expected = 64.0 / (4 * alpha**2 - B**2)**3
            assert abs(cf.evaluate(alpha, B) - expected) < 1e-13

    def test_dataclass_is_immutable(self):
        """ClosedForm is frozen; cannot be mutated."""
        cf = I_HE_cosh_polynomial(0, 0, 0)
        with pytest.raises(Exception):
            cf.L = 1  # frozen dataclass

    def test_degree_properties_finite(self):
        """Degree accessors return finite non-negative integers."""
        cf = I_HE_cosh_polynomial(0, 0, 0)
        assert cf.num_degree_alpha >= 0
        assert cf.num_degree_b2 >= 0
        assert cf.den_degree_alpha >= 6  # at least (4 a^2 - B^2)^3
        assert cf.den_degree_b2 >= 3

    def test_exact_evaluation(self):
        """evaluate_exact at rational (alpha, B) gives exact Fraction."""
        cf = I_HE_cosh_polynomial(0, 0, 0)
        # I(0,0,0; alpha=2, B=0) = 64 / (4*4)^3 = 64/4096 = 1/64
        val = cf.evaluate_exact(Fraction(2), Fraction(0))
        assert val == Fraction(1, 64)


# ---------------------------------------------------------------------------
# Numerical evaluation invariants
# ---------------------------------------------------------------------------

class TestNumericalEvaluation:
    def test_even_in_B(self):
        """I_HE^cosh is even in B: f(alpha, B) = f(alpha, -B)."""
        for L, M, N in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 1)]:
            v_pos = I_HE_cosh_numeric(L, M, N, 1.5, 0.7)
            v_neg = I_HE_cosh_numeric(L, M, N, 1.5, -0.7)
            assert abs(v_pos - v_neg) < 1e-14 * max(1.0, abs(v_pos)), (
                f"({L},{M},{N}): not even in B"
            )

    def test_positive_for_physical_regime(self):
        """I_HE^cosh > 0 for B in the convergence interior (|B| < 2 alpha)."""
        for L, M, N in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]:
            for alpha in [1.0, 1.5, 2.0]:
                for B in [0.0, 0.5 * alpha, 1.5 * alpha]:  # all < 2*alpha
                    val = I_HE_cosh_numeric(L, M, N, alpha, B)
                    assert val > 0, (
                        f"({L},{M},{N}) @ alpha={alpha}, B={B}: val={val}"
                    )

    def test_alpha_negative_rejected(self):
        """alpha <= 0 is rejected."""
        with pytest.raises(ValueError):
            I_HE_cosh_numeric(0, 0, 0, 0.0, 0.0)
        with pytest.raises(ValueError):
            I_HE_cosh_numeric(0, 0, 0, -1.0, 0.0)

    def test_lmn_negative_rejected(self):
        """Negative L, M, N rejected at symbolic-call boundary."""
        with pytest.raises(ValueError):
            I_HE_cosh_symbolic(-1, 0, 0)
        with pytest.raises(ValueError):
            I_HE_cosh_symbolic(0, -1, 0)
        with pytest.raises(ValueError):
            I_HE_cosh_symbolic(0, 0, -1)


# ---------------------------------------------------------------------------
# Disk cache round-trip
# ---------------------------------------------------------------------------

class TestDiskCache:
    """The closed-form table can be saved to disk and loaded back, with
    schema-version invalidation."""

    def test_save_load_round_trip(self, tmp_path):
        """save_table + load_table reproduces the in-memory table."""
        path = str(tmp_path / "test_cache.pkl")
        # Populate a small table.
        table = precompute_table(max_total_degree=2, verbose=False)
        save_table(table, path=path)

        # Force reload from disk into a fresh process? We can't easily do
        # that in-process. Instead verify load_table returns a table
        # equal to the source.
        loaded = load_table(path=path, populate_cache=False)
        assert loaded is not None
        assert set(loaded.keys()) == set(table.keys())
        for key in table:
            assert loaded[key].num_coeffs == table[key].num_coeffs
            assert loaded[key].den_coeffs == table[key].den_coeffs

    def test_stale_schema_invalidates(self, tmp_path):
        """A cache file with a stale schema_version is treated as missing."""
        import pickle
        path = str(tmp_path / "stale_cache.pkl")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        bogus_payload = {
            "schema_version": "v0.0-NOT-REAL",
            "table": {},
        }
        with open(path, "wb") as f:
            pickle.dump(bogus_payload, f)

        loaded = load_table(path=path, populate_cache=False)
        assert loaded is None  # stale schema -> treated as absent

    def test_missing_file_returns_none(self, tmp_path):
        """load_table on a nonexistent path returns None silently."""
        path = str(tmp_path / "does_not_exist.pkl")
        loaded = load_table(path=path, populate_cache=False)
        assert loaded is None


# ---------------------------------------------------------------------------
# Convergence boundary
# ---------------------------------------------------------------------------

class TestConvergenceBoundary:
    def test_at_convergence_boundary_raises(self):
        """At 2*alpha = |B| the denominator vanishes (ZeroDivisionError)."""
        alpha = 1.5
        B = 2.0 * alpha  # exactly the convergence boundary
        with pytest.raises(ZeroDivisionError):
            I_HE_cosh_numeric(0, 0, 0, alpha, B)
