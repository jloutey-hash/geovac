"""Verification of Paper 46 closed-form constant C_3^op(n_max) = sqrt(1 - 1/n_max).

Per CLAUDE.md §13.4a equation-verification protocol.

Paper 46 (`papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex`)
Lemma L3 (lem:L3op, eq:C3op_closed_form, lines 766-784, 811-844) defines the
envelope-aware Lipschitz constant for the operator-norm seminorm
L_op(a) := ||[D_L, a]||_op on the natural chirality-doubled scalar-multiplier
substrate:

    C_3^op(n_max)
      := sup_{2 <= N <= 2*n_max - 1}  (N - 1) / sqrt(N^2 - 1)
       = sup_{2 <= N <= 2*n_max - 1}  sqrt((N - 1) / (N + 1))
       = sqrt((2*n_max - 2) / (2*n_max))
       = sqrt(1 - 1/n_max).

The per-harmonic constant from Paper 38 (eq:C3_per_harmonic) is

    C_3^(N) = sqrt((N - 1) / (N + 1)) = (N - 1) / sqrt(N^2 - 1),

monotone increasing in N, so the supremum is attained at N = 2*n_max - 1
(the Avery-Wen-Avery envelope max from Paper 44 §3).

Existing `tests/test_lorentzian_propinquity.py` covers the Riemannian limit
(N_t = 1) and panel cells of the joint propinquity bound but does not test
the closed-form C_3^op(n_max) directly. This file fills that gap.

Verification strategy (bit-exact via sympy):

  1. Compute C_3^op(n_max) at n_max in {2, 3, 4, 5, 6, 10, 20, 100, 1000}
     by taking the supremum over N in {2, ..., 2*n_max - 1} of the
     per-harmonic constant (N - 1) / sqrt(N^2 - 1).
  2. Verify each equals sqrt(1 - 1/n_max) exactly under symbolic
     simplification (sympy.simplify result is zero).
  3. Verify monotone strict increase: C_3^op(n_max) < C_3^op(n_max + 1)
     for the panel.
  4. Verify limit -> 1^- : C_3^op(1000) - 1 < 0 with magnitude ~1/(2*1000).
  5. Verify per-harmonic monotonicity in N (load-bearing for the supremum
     argument).
  6. Verify boundary values at small n_max numerically:
        n_max = 2: sqrt(1/2)
        n_max = 3: sqrt(2/3)
        n_max = 4: sqrt(3/4)
        n_max = 5: sqrt(4/5)
"""

from __future__ import annotations

import pytest
import sympy as sp

# Total runtime ~18s on the standard panel (dominated by the n_max in
# {100, 1000} enumeration in test_C3op_supremum_equals_closed_form_large
# and test_C3op_strict_upper_bound_by_one). Mark the module as slow per
# CLAUDE.md §14 so the standard fast regression skips it.
pytestmark = pytest.mark.slow


# ---------------------------------------------------------------------------
# Definitions from Paper 46 eq:C3_per_harmonic and eq:C3op_closed_form.
# ---------------------------------------------------------------------------


def C3_per_harmonic(N: int) -> sp.Expr:
    """Paper 38 / Paper 46 eq:C3_per_harmonic.

    C_3^(N) := (N - 1) / sqrt(N^2 - 1) = sqrt((N - 1) / (N + 1)).

    Defined for integer N >= 2.
    """
    if N < 2:
        raise ValueError(f"C3_per_harmonic requires N >= 2, got N = {N}")
    Ns = sp.Integer(N)
    return (Ns - 1) / sp.sqrt(Ns**2 - 1)


def C3op_from_sup(n_max: int) -> sp.Expr:
    """Compute C_3^op(n_max) directly from the supremum definition.

    Paper 46 eq:C3op_closed_form (first equality):
        C_3^op(n_max) := sup_{2 <= N <= 2 n_max - 1} (N - 1) / sqrt(N^2 - 1).

    By symbolic monotonicity (verified separately), the supremum is attained
    at N = 2 n_max - 1, but we DO NOT assume that here -- this routine
    enumerates all N in the envelope range and picks the symbolic max.
    """
    if n_max < 2:
        raise ValueError(f"C3op requires n_max >= 2 (envelope range non-empty), got {n_max}")
    N_lo, N_hi = 2, 2 * n_max - 1
    # Enumerate exact symbolic values.
    vals = [C3_per_harmonic(N) for N in range(N_lo, N_hi + 1)]
    # Take symbolic max via pairwise sympy.Max, then simplify.
    sup = vals[0]
    for v in vals[1:]:
        sup = sp.Max(sup, v)
    return sp.simplify(sup)


def C3op_closed_form(n_max: int) -> sp.Expr:
    """Paper 46 eq:C3op_closed_form (final closed form): sqrt(1 - 1/n_max)."""
    if n_max < 2:
        raise ValueError(f"closed form requires n_max >= 2, got {n_max}")
    return sp.sqrt(1 - sp.Rational(1, n_max))


# ---------------------------------------------------------------------------
# Tests.
# ---------------------------------------------------------------------------


# Panel: small n_max where supremum enumeration is fast, plus large values
# probing the n_max -> infinity limit.
N_MAX_PANEL_SMALL = [2, 3, 4, 5, 6, 7, 8, 10]
N_MAX_PANEL_LARGE = [20, 100, 1000]
N_MAX_PANEL_ALL = N_MAX_PANEL_SMALL + N_MAX_PANEL_LARGE


def test_C3op_supremum_equals_closed_form_small():
    """For n_max in {2, ..., 10}, supremum-defined C_3^op matches sqrt(1 - 1/n_max) bit-exactly."""
    for n_max in N_MAX_PANEL_SMALL:
        sup_val = C3op_from_sup(n_max)
        closed = C3op_closed_form(n_max)
        residual = sp.simplify(sup_val - closed)
        assert residual == 0, (
            f"n_max = {n_max}: supremum {sup_val} != closed form {closed} "
            f"(residual {residual})"
        )


def test_C3op_supremum_equals_closed_form_large():
    """For large n_max (up to 1000), the supremum-defined C_3^op still matches.

    The enumeration is O(n_max); n_max = 1000 needs 1999 sympy values.
    """
    for n_max in N_MAX_PANEL_LARGE:
        sup_val = C3op_from_sup(n_max)
        closed = C3op_closed_form(n_max)
        residual = sp.simplify(sup_val - closed)
        assert residual == 0, (
            f"n_max = {n_max}: supremum {sup_val} != closed form {closed} "
            f"(residual {residual})"
        )


def test_C3op_boundary_values():
    """Bit-exact match to hand-computed boundary values.

    n_max = 2: sqrt(1/2)
    n_max = 3: sqrt(2/3)
    n_max = 4: sqrt(3/4) = sqrt(3) / 2
    n_max = 5: sqrt(4/5) = 2 / sqrt(5)
    """
    expected = {
        2: sp.sqrt(sp.Rational(1, 2)),
        3: sp.sqrt(sp.Rational(2, 3)),
        4: sp.sqrt(sp.Rational(3, 4)),
        5: sp.sqrt(sp.Rational(4, 5)),
    }
    for n_max, exp in expected.items():
        assert sp.simplify(C3op_closed_form(n_max) - exp) == 0, (
            f"n_max = {n_max}: closed form != {exp}"
        )
        # And the supremum agrees.
        assert sp.simplify(C3op_from_sup(n_max) - exp) == 0, (
            f"n_max = {n_max}: supremum != {exp}"
        )


def test_C3op_supremum_attained_at_envelope_max():
    """The supremum sup_{2<=N<=2 n_max - 1} (N-1)/sqrt(N^2 - 1) is attained at N = 2 n_max - 1.

    This is the load-bearing claim used in the proof of eq:C3op_closed_form
    (Paper 46 §sec:L3 Step C, line 841).
    """
    for n_max in N_MAX_PANEL_SMALL:
        N_envmax = 2 * n_max - 1
        per_envmax = C3_per_harmonic(N_envmax)
        # Compare to all other values in the envelope.
        for N in range(2, N_envmax + 1):
            diff = sp.simplify(per_envmax - C3_per_harmonic(N))
            # diff should be >= 0; verify with sp.simplify().
            assert diff >= 0, (
                f"n_max = {n_max}: per-harmonic({N_envmax}) = {per_envmax} "
                f"is not >= per-harmonic({N}) = {C3_per_harmonic(N)} (diff {diff})"
            )


def test_C3op_per_harmonic_monotone_in_N():
    """Per-harmonic C_3^(N) is strictly monotone increasing in N for N >= 2.

    This underwrites the supremum-at-N=2*n_max-1 fact above.  Verified
    symbolically: C_3^(N+1) > C_3^(N) <=> N(N+2) > N^2 - 1 <=> 2N > -1 (true).
    """
    # Symbolic verification at sample points.
    for N in range(2, 30):
        diff = sp.simplify(C3_per_harmonic(N + 1) - C3_per_harmonic(N))
        assert diff > 0, (
            f"per-harmonic({N + 1}) - per-harmonic({N}) = {diff} is not strictly positive"
        )


def test_C3op_monotone_increase_in_n_max():
    """C_3^op(n_max) is strictly monotone increasing in n_max for n_max >= 2.

    Symbolic verification: sqrt(1 - 1/(n+1)) - sqrt(1 - 1/n) > 0 since
    1 - 1/(n+1) > 1 - 1/n.
    """
    for n_max in N_MAX_PANEL_SMALL:
        diff = sp.simplify(C3op_closed_form(n_max + 1) - C3op_closed_form(n_max))
        assert diff > 0, (
            f"C3op({n_max + 1}) - C3op({n_max}) = {diff} is not strictly positive"
        )


def test_C3op_strict_upper_bound_by_one():
    """C_3^op(n_max) < 1 for every finite n_max (Paper 46 line 1121).

    Used as the multiplier dominance in C_5^joint = max(1, 1, C_3^op, 0) = 1.
    """
    for n_max in N_MAX_PANEL_ALL:
        c = C3op_closed_form(n_max)
        # sqrt(1 - 1/n_max) < 1 <=> 1 - 1/n_max < 1 <=> 1/n_max > 0.  Symbolic.
        assert sp.simplify(1 - c) > 0, (
            f"C3op({n_max}) = {c} is not strictly less than 1"
        )


def test_C3op_asymptotic_to_one():
    """C_3^op(n_max) -> 1^- as n_max -> infinity (Paper 46 eq:C3op_closed_form)."""
    n = sp.Symbol("n", positive=True, integer=True)
    closed_sym = sp.sqrt(1 - sp.Rational(1, 1) / n)
    lim = sp.limit(closed_sym, n, sp.oo)
    assert lim == 1, f"limit is {lim}, expected 1"

    # And the approach is from below at a 1/(2 n_max) rate (asymptotic
    # expansion sqrt(1 - 1/n) = 1 - 1/(2n) + O(1/n^2)):
    series_expansion = sp.series(closed_sym, n, sp.oo, 3).removeO()
    expected_leading = 1 - sp.Rational(1, 2) / n
    diff = sp.simplify(series_expansion - expected_leading)
    # The leading-order behaviour matches 1 - 1/(2n).
    assert diff.subs(n, 1000) < sp.Rational(1, 10**5), (
        f"asymptotic 1 - 1/(2n) leading order off: residual subbed at n=1000 is {diff.subs(n, 1000)}"
    )


def test_C3op_numerical_panel():
    """Numerical sanity-check at the panel cells used in Paper 46 §sec:numerical.

    Paper 46 Table 5.1 (line 1205-1218): panel cells (n_max, N_t) in
    {(2,3), (3,5), (4,7)}.  C_3^op depends only on n_max, so the test
    reduces to three numerical values that should agree with sqrt(1 - 1/n_max).
    """
    expected_numerical = {
        2: float(sp.sqrt(sp.Rational(1, 2))),
        3: float(sp.sqrt(sp.Rational(2, 3))),
        4: float(sp.sqrt(sp.Rational(3, 4))),
    }
    for n_max, ref in expected_numerical.items():
        val = float(C3op_closed_form(n_max))
        assert abs(val - ref) < 1e-15, (
            f"n_max = {n_max}: numerical mismatch {val} vs {ref}"
        )


def test_C3op_envelope_tightness_remark_4_2():
    """Paper 46 Remark 4.2 (rem:env_tightness): at (n_max, N_t) = (3, 5),
    the (N=4, L=1, M=+/-1) generator has LHS/RHS ratio
    C_3^(4) / C_3^op(3) = sqrt(3/5) / sqrt(2/3) approx 0.949.

    Verified bit-exactly via sympy.
    """
    n_max = 3
    N_interior = 4
    ratio_sym = C3_per_harmonic(N_interior) / C3op_closed_form(n_max)
    ratio_expected = sp.sqrt(sp.Rational(3, 5)) / sp.sqrt(sp.Rational(2, 3))
    residual = sp.simplify(ratio_sym - ratio_expected)
    assert residual == 0, f"ratio mismatch: {ratio_sym} vs {ratio_expected} (residual {residual})"
    # Numerical value of the ratio (Paper 46 line 884: ~0.949).
    assert abs(float(ratio_sym) - 0.9486832980505137) < 1e-12, (
        f"numerical ratio {float(ratio_sym)} != 0.9486..."
    )
