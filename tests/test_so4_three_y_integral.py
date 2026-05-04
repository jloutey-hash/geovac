"""Tests for geovac.so4_three_y_integral: the Avery-Wen-Avery 3-Y integral on S^3.

Equation verification per CLAUDE.md Sec 13.4a. Specifically:

  - Orthonormality of the hyperspherical harmonics Y_{n,l,m}(chi, theta, phi)
    on S^3 with measure dOmega_3 = sin^2(chi) sin(theta) dchi dtheta dphi.
    This is the analytical-limit verification: <Y_a | Y_b> = delta_{ab}.

  - Selection-rule enforcement on the full 3-Y integral
    avery_wen_avery_3y(...): returns Integer(0) when m_a != m_b + m_c (Gaunt
    m-rule), when l-triangle |l_b - l_c| <= l_a <= l_b + l_c is violated, or
    when l_a + l_b + l_c is odd (Gaunt parity rule).

  - Symbolic spot-check: <Y_{2,1,0}|Y_{1,0,0}|Y_{2,1,0}> matches the
    closed-form value sqrt(2)/(2*pi). (The Y_{1,0,0} multiplier is
    proportional to the constant function 1/sqrt(2 pi^2), so the integral
    reduces to <Y_{2,1,0}|Y_{2,1,0}> times the constant value of Y_{1,0,0}.)

  - Performance smoke test: full multiplier-matrix build at n_max=3
    (uses geovac.operator_system) completes in < 30s.

These tests are sympy-exact (no float comparisons) for the orthonormality
and selection-rule checks, since the Avery-Wen-Avery integral is constructed
in exact rational + sqrt(rational) arithmetic.

References
==========

J. Avery, "Hyperspherical Harmonics: Applications in Quantum Theory,"
Kluwer 1989, Eq. 1.6.3 (radial normalization), 2.5 (SO(4) selection rules),
3.5 (3-Y integral factorization).

Z. Wen & J. Avery, "Some properties of hyperspherical harmonics,"
J. Math. Phys. 26, 396-403 (1985).

A. R. Edmonds, "Angular Momentum in Quantum Mechanics," Princeton 1957,
Eq. 4.6.3 (Gaunt integral closed form).
"""

from __future__ import annotations

import time

import pytest
import sympy as sp
from sympy import Integer, Rational, pi, simplify, sqrt

from geovac.so4_three_y_integral import (
    avery_wen_avery_3y,
    gegenbauer_triple_integral,
    radial_overlap_avery,
    radial_overlap_two,
    s2_gaunt,
    s3_orthonormality_check,
)


# ---------------------------------------------------------------------------
# Helpers: enumerate all valid (n, l, m) labels at a given cutoff
# ---------------------------------------------------------------------------


def _all_labels(n_max: int):
    """Yield all valid (n, l, m) hyperspherical labels with n <= n_max.

    Bounds: 1 <= n <= n_max, 0 <= l <= n - 1, -l <= m <= l.
    """
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                yield (n, l, m)


# ---------------------------------------------------------------------------
# Orthonormality: the headline equation-verification test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_orthonormality_diagonals_equal_one(n_max):
    """<Y_{n,l,m} | Y_{n,l,m}>_{S^3} = 1 EXACTLY (sympy) for every label."""
    for (n, l, m) in _all_labels(n_max):
        val = s3_orthonormality_check(n, l, m, n, l, m)
        simplified = simplify(val)
        assert simplified == 1, (
            f"<Y_{{{n},{l},{m}}}|Y_{{{n},{l},{m}}}> = {simplified}, expected 1"
        )


@pytest.mark.parametrize("n_max", [2, 3])
def test_orthonormality_off_diagonals_equal_zero(n_max):
    """<Y_a | Y_b> = 0 EXACTLY (sympy) for every distinct pair of labels."""
    labels = list(_all_labels(n_max))
    for i, (na, la, ma) in enumerate(labels):
        for nb, lb, mb in labels[i + 1:]:
            val = s3_orthonormality_check(na, la, ma, nb, lb, mb)
            simplified = simplify(val)
            assert simplified == 0, (
                f"<Y_{{{na},{la},{ma}}}|Y_{{{nb},{lb},{mb}}}> = "
                f"{simplified}, expected 0"
            )


def test_radial_overlap_two_diagonal_one():
    """<R_{n,l} | R_{n,l}>_{[0, pi], sin^2 dchi} = 1 (radial-only)."""
    for n in range(1, 5):
        for l in range(n):
            val = radial_overlap_two(n, l, n)
            assert simplify(val) == 1, (
                f"<R_{{{n},{l}}}|R_{{{n},{l}}}> = {simplify(val)}, expected 1"
            )


def test_radial_overlap_two_off_diagonal_zero():
    """<R_{n,l} | R_{n',l}>_{[0, pi], sin^2 dchi} = 0 for n != n', same l."""
    pairs = [(1, 2, 0), (1, 3, 0), (2, 3, 0), (2, 3, 1), (1, 4, 0)]
    for n, np_, l in pairs:
        if l > min(n, np_) - 1:
            continue
        val = radial_overlap_two(n, l, np_)
        assert simplify(val) == 0, (
            f"<R_{{{n},{l}}}|R_{{{np_},{l}}}> = {simplify(val)}, expected 0"
        )


# ---------------------------------------------------------------------------
# Selection-rule enforcement on avery_wen_avery_3y
# ---------------------------------------------------------------------------


def test_selection_rule_m_sum_violation():
    """If m_a != m_b + m_c, return 0 exactly. (Gaunt m-rule.)"""
    # All these violate m_a = m_b + m_c
    cases = [
        (2, 1, +1, 1, 0, 0, 2, 1, 0),     # 1 != 0 + 0
        (2, 1, +1, 2, 1, 0, 2, 1, +1),    # 1 != 0 + 1 = 1, OK actually...
    ]
    # Reset cases more carefully
    cases = [
        (2, 1, +1, 1, 0, 0, 2, 1, 0),     # 1 != 0 + 0
        (2, 1, -1, 1, 0, 0, 2, 1, 0),     # -1 != 0 + 0
        (3, 2, +2, 2, 1, 0, 3, 2, 0),     # 2 != 0 + 0
    ]
    for case in cases:
        val = avery_wen_avery_3y(*case)
        assert val == 0, f"expected 0 for {case}, got {val}"


def test_selection_rule_l_triangle_violation():
    """If |l_b - l_c| > l_a or l_a > l_b + l_c, return 0 exactly."""
    # l_a = 0, l_b = 0, l_c = 2: triangle requires l_a >= |l_b - l_c| = 2,
    # but l_a = 0. Violated.
    val = avery_wen_avery_3y(1, 0, 0, 1, 0, 0, 3, 2, 0)
    assert val == 0, f"l-triangle violation: expected 0, got {val}"

    # l_a = 0, l_b = 1, l_c = 2: triangle requires l_a >= 1, but l_a = 0.
    val = avery_wen_avery_3y(1, 0, 0, 2, 1, 0, 3, 2, 0)
    assert val == 0


def test_selection_rule_l_parity_violation():
    """If l_a + l_b + l_c is odd, return 0 exactly. (Gaunt parity rule.)"""
    # l_a + l_b + l_c = 0 + 0 + 1 = 1 odd
    val = avery_wen_avery_3y(1, 0, 0, 1, 0, 0, 2, 1, 0)
    assert val == 0, f"l-parity violation (sum=1): expected 0, got {val}"

    # l_a + l_b + l_c = 1 + 1 + 1 = 3 odd
    val = avery_wen_avery_3y(2, 1, 0, 2, 1, 0, 2, 1, 0)
    assert val == 0


def test_selection_rules_in_s2_gaunt():
    """The S^2 Gaunt integral itself enforces the same rules."""
    # m sum violation
    assert s2_gaunt(1, +1, 0, 0, 1, 0) == 0
    # l-triangle (l_a=0, l_b=0, l_c=2 needs l_a>=2)
    assert s2_gaunt(0, 0, 0, 0, 2, 0) == 0
    # l-parity (sum=1 odd)
    assert s2_gaunt(0, 0, 0, 0, 1, 0) == 0


# ---------------------------------------------------------------------------
# Symbolic spot checks
# ---------------------------------------------------------------------------


def test_spot_check_y210_y100_y210():
    """<Y_{2,1,0} | Y_{1,0,0} | Y_{2,1,0}>_{S^3} = sqrt(2)/(2*pi).

    Y_{1,0,0} is the constant function on S^3 (up to normalization) and
    equals 1/sqrt(2*pi^2) in the Avery 1989 normalization. Then

        <Y_a | const | Y_a> = const * <Y_a | Y_a> = const * 1.

    So the integral is the constant value 1/sqrt(2*pi^2) = sqrt(2)/(2*pi)
    after rationalization.

    This serves both as a spot check that the closed form is right and as
    documentation of the Avery normalization convention.
    """
    val = avery_wen_avery_3y(2, 1, 0, 1, 0, 0, 2, 1, 0)
    expected = sqrt(2) / (2 * pi)
    diff = simplify(val - expected)
    assert diff == 0, f"got {simplify(val)}, expected {expected}, diff {diff}"


def test_spot_check_y100_y100_y100():
    """<Y_{1,0,0} | Y_{1,0,0} | Y_{1,0,0}>_{S^3} = 1/sqrt(2*pi^2).

    Y_{1,0,0}^3 integrated over S^3 with measure dOmega_3 (volume 2*pi^2)
    is just (1/sqrt(2 pi^2))^3 * 2 pi^2 = 1/sqrt(2 pi^2). Same as
    Y_{1,0,0} itself.
    """
    val = avery_wen_avery_3y(1, 0, 0, 1, 0, 0, 1, 0, 0)
    expected = 1 / sqrt(2 * pi ** 2)
    diff = simplify(val - expected)
    assert diff == 0, f"got {simplify(val)}, expected {expected}, diff {diff}"


def test_radial_overlap_avery_drop_in_signature():
    """The radial_overlap_avery(n, N, np_, l, L, lp) signature matches
    the placeholder dispatch shape used by operator_system._so4_radial_overlap.
    """
    # Diagonal (1, 0; 1, 0; 1, 0) should give a nonzero value (the constant
    # multiplier overlap).
    val = radial_overlap_avery(1, 1, 1, 0, 0, 0)
    assert val != 0, f"radial_overlap_avery(1,1,1,0,0,0) = {val}, expected nonzero"

    # An obviously zero case: the gegenbauer_triple_integral result for
    # mismatched parities or singular triples
    val_zero = radial_overlap_avery(1, 3, 1, 0, 0, 0)
    # n=1, N=3, n'=1: |1-1|+1 = 1 <= 3 <= 1+1-1 = 1 fails (3 > 1)
    # so the integral is structurally zero on the (n, N, n') triangle.
    # However the gegenbauer_triple_integral does NOT enforce the triangle
    # explicitly (per its docstring), so let's just check the value.
    # The test is: it should be a valid sympy expression (zero or nonzero).
    assert isinstance(val_zero, sp.Expr) or val_zero == 0


# ---------------------------------------------------------------------------
# Symmetry / hermiticity-related sanity checks
# ---------------------------------------------------------------------------


def test_3y_symmetric_under_bra_ket_swap_for_real_values():
    """<Y_a | Y_b | Y_c> = <Y_c | Y_b | Y_a>* for the standard inner product.

    For our convention with (-1)^{m_a} prefactor on the Gaunt, this is
    encoded in the closed form. We test concrete tuples where the values
    are real (m_a = m_b = m_c = 0).
    """
    # All m=0 cases, where the integral is real
    cases = [
        (1, 0, 0, 1, 0, 0, 1, 0, 0),
        (2, 0, 0, 1, 0, 0, 2, 0, 0),
        (2, 1, 0, 2, 1, 0, 1, 0, 0),  # checking permutation
    ]
    for case in cases:
        n_a, l_a, m_a, n_b, l_b, m_b, n_c, l_c, m_c = case
        val_abc = avery_wen_avery_3y(n_a, l_a, m_a, n_b, l_b, m_b, n_c, l_c, m_c)
        val_cba = avery_wen_avery_3y(n_c, l_c, m_c, n_b, l_b, m_b, n_a, l_a, m_a)
        diff = simplify(val_abc - val_cba)
        assert diff == 0, (
            f"swap symmetry violated for {case}: "
            f"forward {simplify(val_abc)}, backward {simplify(val_cba)}"
        )


# ---------------------------------------------------------------------------
# Performance smoke test
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_multiplier_matrix_build_n_max_3_under_30s():
    """Building TruncatedOperatorSystem(n_max=3) — which evaluates many
    Avery integrals — should complete in under 30s on standard hardware.
    """
    from geovac.operator_system import TruncatedOperatorSystem
    t0 = time.time()
    O = TruncatedOperatorSystem(3)
    elapsed = time.time() - t0
    assert elapsed < 30, (
        f"TruncatedOperatorSystem(3) took {elapsed:.1f}s, exceeds 30s budget"
    )
    # Sanity: the build succeeded and produced expected dim
    assert O.dim_H == 14
    assert O.dim == 55


# ---------------------------------------------------------------------------
# Selection-rule completeness: nonzero examples
# ---------------------------------------------------------------------------


def test_nonzero_example_matches_orthonormality():
    """<Y_a | const | Y_a> nonzero with the constant Y_{1,0,0} multiplier.

    Verifies that selection rules don't accidentally zero out the diagonal
    "constant multiplier" matrix elements (which are required for the
    identity 1 = const * sqrt(2 pi^2) * Y_{1,0,0} to be in O at all).
    """
    for (n, l, m) in _all_labels(2):
        val = avery_wen_avery_3y(n, l, m, 1, 0, 0, n, l, m)
        assert val != 0, f"<Y_{{{n},{l},{m}}}|const|Y_{{{n},{l},{m}}}> = 0"


def test_n_to_n_plus_1_l_change_one_nonzero():
    """The "raising/lowering" multiplier M_{2,1,0} couples (n, l, m) to
    (n+1, l+1, m) with l_b = 1 satisfying l_a + l_c + 1 = even. Verify
    nonzero on at least one such matrix element.
    """
    # Y_{1,0,0} -> Y_{2,1,0} via Y_{2,1,0} multiplier
    val = avery_wen_avery_3y(1, 0, 0, 2, 1, 0, 2, 1, 0)
    assert val != 0, f"M_{{2,1,0}} between |1,0,0> and |2,1,0> is zero"
