"""Tests for geovac/hylleraas_eckart_recurrence.py.

Track 1.5 of the Hylleraas-Eckart sprint — the algebraic recurrence
engine that replaces sympy as the closed-form generator.

Verifies bit-exact ClosedForm equivalence vs the sympy-based oracle in
geovac/hylleraas_eckart_closed_forms.py across the d <= 6 panel, plus
machine-precision numerical evaluations at sample (alpha, B).
"""

from __future__ import annotations

import time
from fractions import Fraction

import pytest

from geovac.hylleraas_eckart_recurrence import (
    I_HE_cosh_polynomial_recurrence,
    _G_2M,
    _F_sinh,
    _F_cosh,
    precompute_table_recurrence,
)
from geovac.hylleraas_eckart_closed_forms import I_HE_cosh_polynomial
from geovac.hylleraas_r12 import hylleraas_master_int


# ---------------------------------------------------------------------------
# Direct closed-form sanity at (0, 0, 0): 64 / (4 alpha^2 - B^2)^3
# ---------------------------------------------------------------------------

class TestClosedFormSanityRecurrence:
    """The (0, 0, 0) case has known closed form 64 / (4 alpha^2 - B^2)^3."""

    @pytest.mark.parametrize("alpha,B,expected_factor", [
        (1.0, 0.0, 1.0 / 1.0**6),
        (1.5, 1.0, 64.0 / 8**3),
        (1.5, 0.5, 64.0 / (4 * 1.5**2 - 0.25)**3),
        (1.6875, 0.0, 1.0 / 1.6875**6),
    ])
    def test_zero_zero_zero(self, alpha: float, B: float, expected_factor: float):
        cf = I_HE_cosh_polynomial_recurrence(0, 0, 0)
        val = cf.evaluate(alpha, B)
        assert abs(val - expected_factor) < 1e-14, (
            f"alpha={alpha}, B={B}: got {val}, expected {expected_factor}"
        )


# ---------------------------------------------------------------------------
# Numerical equivalence to sympy oracle across d <= 6 panel
# ---------------------------------------------------------------------------

class TestRecurrenceMatchesSympy:
    """For every (L, M, N) with L+2M+N <= 6 and at multiple (alpha, B),
    the recurrence-engine ClosedForm must evaluate identically to the
    sympy-based ClosedForm."""

    SAMPLE_AB = [(1.5, 0.5), (2.0, 0.8), (1.0, 0.3), (1.6875, 0.0), (1.2, 1.0)]

    def test_panel_d6(self):
        worst_rel = 0.0
        worst_loc = None
        n_cells = 0
        for L in range(7):
            for M in range((6 - L) // 2 + 1):
                for N in range(6 - L - 2 * M + 1):
                    n_cells += 1
                    cf_rec = I_HE_cosh_polynomial_recurrence(L, M, N)
                    cf_sym = I_HE_cosh_polynomial(L, M, N)
                    for alpha, B in self.SAMPLE_AB:
                        v_rec = cf_rec.evaluate(alpha, B)
                        v_sym = cf_sym.evaluate(alpha, B)
                        rel = abs(v_rec - v_sym) / max(abs(v_sym), 1e-30)
                        if rel > worst_rel:
                            worst_rel = rel
                            worst_loc = ((L, M, N), (alpha, B))
        assert worst_rel < 1e-12, (
            f"Worst rel_err {worst_rel:.2e} at {worst_loc} (across {n_cells} cells)"
        )


# ---------------------------------------------------------------------------
# B -> 0 reduction to single-alpha master integral
# ---------------------------------------------------------------------------

class TestBZeroReductionRecurrence:
    def test_panel_b_zero(self):
        """At B=0, recurrence engine reproduces single-alpha master to machine precision."""
        alpha = 1.6875
        max_total = 6
        worst_rel = 0.0
        for L in range(max_total + 1):
            for M in range((max_total - L) // 2 + 1):
                for N in range(max_total - L - 2 * M + 1):
                    single = hylleraas_master_int(L, M, N, alpha)
                    cf = I_HE_cosh_polynomial_recurrence(L, M, N)
                    eckart = cf.evaluate(alpha, 0.0)
                    rel = abs(single - eckart) / abs(single)
                    worst_rel = max(worst_rel, rel)
        assert worst_rel < 1e-12, f"Worst rel_err at B=0: {worst_rel:.2e}"


# ---------------------------------------------------------------------------
# Performance: recurrence is much faster than sympy
# ---------------------------------------------------------------------------

class TestPerformance:
    @pytest.mark.slow
    def test_recurrence_faster_than_sympy(self):
        """At (3, 1, 2), recurrence must be at least 100x faster than sympy.

        Marked slow because we have to clear the in-process cache and force
        the sympy oracle to do its full ~4s symbolic computation.
        """
        from geovac.hylleraas_eckart_closed_forms import _CACHE_TABLE
        L, M, N = 3, 1, 2
        # Warm caches (lru_cache on G_2M and F_sinh/F_cosh).
        I_HE_cosh_polynomial_recurrence(L, M, N)
        # Time recurrence (warm).
        t0 = time.time()
        for _ in range(10):
            I_HE_cosh_polynomial_recurrence(L, M, N)
        t_rec = (time.time() - t0) / 10
        # Time sympy: must bypass the in-memory cache (which is populated
        # by either recurrence or disk-loaded ClosedForms).
        _CACHE_TABLE.pop((L, M, N), None)
        t0 = time.time()
        I_HE_cosh_polynomial(L, M, N)
        t_sym = time.time() - t0
        # Speedup should be > 100x for (3,1,2).
        speedup = t_sym / max(t_rec, 1e-9)
        assert speedup > 100, (
            f"Speedup {speedup:.1f}x at (3,1,2) below threshold "
            f"(t_rec={t_rec*1000:.2f}ms, t_sym={t_sym*1000:.1f}ms)"
        )


# ---------------------------------------------------------------------------
# Sub-component sanity: G_{2M}(u, B) values
# ---------------------------------------------------------------------------

class TestGRecurrence:
    """Direct sanity on the G_{2M}(u, B) recurrence outputs."""

    def test_G_0_structure(self):
        """G_0 = 2 sinh(Bu) / B  =>  single term coeff=2, b_exp=-1, kind=sinh."""
        G0 = _G_2M(0)
        assert len(G0.terms) == 1
        t = G0.terms[0]
        assert t.coeff == Fraction(2)
        assert t.b_exp == -1
        assert t.s_pow == 0
        assert t.u_pow == 0
        assert t.kind == 'sinh_Bu'

    def test_G_2_structure(self):
        """G_2 = 2 u^2 sinh(Bu)/B - 4 u cosh(Bu)/B^2 + 4 sinh(Bu)/B^3."""
        G2 = _G_2M(1)
        terms_dict = {(t.b_exp, t.u_pow, t.kind): t.coeff for t in G2.terms}
        # 2 u^2 sinh(Bu)/B  =>  (b_exp=-1, u_pow=2, sinh_Bu) -> 2
        assert terms_dict.get((-1, 2, 'sinh_Bu')) == Fraction(2)
        # -4 u cosh(Bu)/B^2 =>  (b_exp=-2, u_pow=1, cosh_Bu) -> -4
        assert terms_dict.get((-2, 1, 'cosh_Bu')) == Fraction(-4)
        # +4 sinh(Bu)/B^3   =>  (b_exp=-3, u_pow=0, sinh_Bu) -> 4
        assert terms_dict.get((-3, 0, 'sinh_Bu')) == Fraction(4)


# ---------------------------------------------------------------------------
# Sub-component sanity: F_sinh, F_cosh
# ---------------------------------------------------------------------------

class TestFRecurrence:
    def test_F_sinh_0(self):
        """F^sinh_0(s, B) = (cosh(Bs) - 1)/B  =>  cosh_Bs/B - 1/B."""
        F = _F_sinh(0)
        terms_dict = {(t.b_exp, t.s_pow, t.kind): t.coeff for t in F.terms}
        assert terms_dict.get((-1, 0, 'cosh_Bs')) == Fraction(1)
        assert terms_dict.get((-1, 0, 'one')) == Fraction(-1)

    def test_F_cosh_0(self):
        """F^cosh_0(s, B) = sinh(Bs) / B."""
        F = _F_cosh(0)
        terms_dict = {(t.b_exp, t.s_pow, t.kind): t.coeff for t in F.terms}
        assert terms_dict.get((-1, 0, 'sinh_Bs')) == Fraction(1)
        assert len([t for t in F.terms if t.coeff != 0]) == 1

    def test_F_sinh_1(self):
        """F^sinh_1(s, B) = s cosh(Bs)/B - sinh(Bs)/B^2."""
        F = _F_sinh(1)
        terms_dict = {(t.b_exp, t.s_pow, t.kind): t.coeff for t in F.terms}
        assert terms_dict.get((-1, 1, 'cosh_Bs')) == Fraction(1)
        assert terms_dict.get((-2, 0, 'sinh_Bs')) == Fraction(-1)
