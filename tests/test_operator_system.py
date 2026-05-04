"""Tests for geovac.operator_system: the Connes-vS truncated operator system
on the Fock-projected S^3.

Verifies, at n_max in {2, 3, 4}:

  - Identity matrix lies in O (the constant function f = 1 is a multiplier).
  - O is *-closed (each generator's conjugate transpose is in O).
  - Witness pair (a, b in O with ab not in O) exists at n_max = 2.
  - dim(O) < N^2 (O is strictly smaller than its C*-envelope).
  - Propagation number prop(O_{n_max}) = 2 at n_max = 2, 3, 4
    (matching Connes-vS Toeplitz S^1 prop = 2 independent of n).
  - Robustness: the propagation number is independent of the radial-overlap
    placeholder choice.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from geovac.operator_system import (
    HyperLabel,
    TruncatedOperatorSystem,
    hilbert_basis,
    hilbert_dim,
    operator_system_dim,
    operator_system_products,
    propagation_number,
    witness_pair,
)


# ---------------------------------------------------------------------------
# Basic sanity
# ---------------------------------------------------------------------------


def test_hilbert_dim_table():
    """N(n_max) = sum n^2 from n=1 to n_max."""
    assert hilbert_dim(1) == 1
    assert hilbert_dim(2) == 5
    assert hilbert_dim(3) == 14
    assert hilbert_dim(4) == 30
    assert hilbert_dim(5) == 55


def test_hilbert_basis_consistency():
    """basis size matches hilbert_dim."""
    for n_max in (1, 2, 3, 4):
        b = hilbert_basis(n_max)
        assert len(b) == hilbert_dim(n_max)


def test_hyperlabel_validation():
    """HyperLabel rejects out-of-range quantum numbers."""
    HyperLabel(n=1, l=0, m=0)        # OK
    HyperLabel(n=2, l=1, m=-1)       # OK
    with pytest.raises(ValueError):
        HyperLabel(n=0, l=0, m=0)
    with pytest.raises(ValueError):
        HyperLabel(n=1, l=1, m=0)    # l = n violates l <= n - 1
    with pytest.raises(ValueError):
        HyperLabel(n=2, l=1, m=2)


# ---------------------------------------------------------------------------
# Operator system construction at small n_max
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_identity_in_O(n_max):
    """The identity matrix lies in O (constant multiplier f = 1).

    The N=1, L=0, M=0 multiplier is the trivial irrep on S^3 (constant
    function), and its matrix in the orthonormal hyperspherical basis is
    proportional to the identity.
    """
    O = TruncatedOperatorSystem(n_max)
    in_O, residual = O.identity_in_O()
    assert in_O, f"identity not in O at n_max={n_max}, residual={residual:.3e}"
    assert residual < 1e-10


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_star_closed(n_max):
    """O is *-closed: every generator's conjugate transpose is in O."""
    O = TruncatedOperatorSystem(n_max)
    star_ok, failures = O.is_star_closed()
    assert star_ok, f"*-closure violated at n_max={n_max}, {len(failures)} failures"


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_O_strictly_smaller_than_envelope(n_max):
    """dim(O) < N^2: O is genuinely an operator system, not its C*-envelope."""
    O = TruncatedOperatorSystem(n_max)
    assert O.dim < O.envelope_dim, (
        f"O at n_max={n_max} has dim {O.dim} = N^2 = {O.envelope_dim}; "
        f"would be the full matrix algebra"
    )
    # Also verify it's strictly bigger than identity-only.
    assert O.dim > 1


def test_dim_O_table():
    """Sanity: dim(O) at small n_max matches the computed values.

    These values were computed with the standard SO(4) selection rules and
    the structural placeholder. They are stable under different placeholder
    choices because the dim is determined by the support pattern, not the
    placeholder values.
    """
    expected = {
        2: 14,
        3: 55,
        4: 140,
    }
    for n_max, dim_expected in expected.items():
        O = TruncatedOperatorSystem(n_max)
        assert O.dim == dim_expected, (
            f"dim(O_{n_max}) = {O.dim}, expected {dim_expected}"
        )


# ---------------------------------------------------------------------------
# Witness pair: a, b in O with ab not in O
# ---------------------------------------------------------------------------


def test_witness_pair_at_nmax_2():
    """At n_max = 2, exhibit a, b in O with ab not in O.

    Take a = M_{N=2, L=1, M=0} (a "raising" multiplier from n -> n+1 with
    Delta l = +/- 1, Delta m = 0). Then b = a^* is the "lowering" partner.
    The product ab represents "raise then lower" — at the top shell n=2 the
    raising step sends the state to n=3 (outside the truncation), where it
    is killed by P_2, producing a deficit on the top-shell diagonal that no
    element of O can supply.
    """
    O = TruncatedOperatorSystem(2)
    a, b, ab, test = witness_pair(O, N_target=2, L_target=1, M_target=0)
    assert a is not None, "M_{2,1,0} generator not present in O"
    in_O, residual = test
    assert not in_O, (
        f"witness pair ab unexpectedly in O at n_max=2 (residual={residual:.3e})"
    )
    # Significant residual means ab is genuinely outside O, not a numerical
    # artifact.
    assert residual > 1e-6, (
        f"ab projection residual {residual:.3e} too small to certify ab not in O"
    )


def test_witness_pair_at_nmax_3():
    """Same witness pair test at n_max = 3 (deeper truncation, smaller
    top-shell deficit ratio but still non-zero)."""
    O = TruncatedOperatorSystem(3)
    a, b, ab, test = witness_pair(O, N_target=2, L_target=1, M_target=0)
    assert a is not None
    in_O, residual = test
    assert not in_O
    assert residual > 1e-6


def test_witness_a_and_b_individually_in_O():
    """Sanity: a = M_{2,1,0} and b = a^* are individually in O even though
    ab is not. This is the defining feature of an operator system: closed
    under * and linear combinations, but not under multiplication."""
    O = TruncatedOperatorSystem(2)
    a, b, ab, _ = witness_pair(O)
    a_in, _ = O.contains(a)
    b_in, _ = O.contains(b)
    assert a_in
    assert b_in


# ---------------------------------------------------------------------------
# Propagation number
# ---------------------------------------------------------------------------


def test_propagation_number_nmax_2():
    """prop(O_2) = 2 (smallest n_max where the test is meaningful)."""
    O = TruncatedOperatorSystem(2)
    prop, dims = propagation_number(O, max_k=4)
    assert prop == 2, f"prop(O_2) = {prop}, expected 2; dims = {dims}"
    assert dims[-1] == O.envelope_dim
    assert dims[0] < O.envelope_dim


def test_propagation_number_nmax_3():
    """prop(O_3) = 2 (Toeplitz S^1 conjecture verification)."""
    O = TruncatedOperatorSystem(3)
    prop, dims = propagation_number(O, max_k=4)
    assert prop == 2, f"prop(O_3) = {prop}, expected 2; dims = {dims}"
    assert dims[0] == 55
    assert dims[1] == 196


def test_propagation_number_nmax_4():
    """prop(O_4) = 2.

    Marked slow because n_max=4 has dim_H = 30 and ~140 generators, and we
    need to compute O^2 = 140 * 140 = 19600 product matrices to check rank.
    Total runtime ~30s.
    """
    O = TruncatedOperatorSystem(4)
    prop, dims = propagation_number(O, max_k=3)
    assert prop == 2, f"prop(O_4) = {prop}, expected 2; dims = {dims}"
    assert dims[0] == 140
    assert dims[1] == 900


# ---------------------------------------------------------------------------
# Robustness: prop is independent of placeholder choice
# ---------------------------------------------------------------------------


def test_propagation_number_robust_to_placeholder():
    """prop is the same under a different (rational) radial-overlap
    placeholder, confirming that the result depends on support pattern only,
    not on the specific values of the radial overlap.

    Note: this test is round-2's invariance check. Sprint WH1-R3.1 swapped
    the default radial-overlap function from the round-2 placeholder to the
    genuine Avery-Wen-Avery integral, so the regression now temporarily
    switches the dispatch point (`_so4_radial_overlap`) to two distinct
    rational placeholders, verifies prop = 2 each time, and restores the
    Avery default.
    """
    from geovac import operator_system as ops
    original = ops._so4_radial_overlap

    def alternate(n, N, np_, l, L, lp):
        if not (abs(n - np_) + 1 <= N <= n + np_ - 1):
            return sp.Integer(0)
        if not L <= N - 1:
            return sp.Integer(0)
        if N == 1:
            return sp.Integer(1)
        # Different rational placeholder
        return sp.Rational(2 * n + 3 * np_ + 5 * N + 7 * l + 11 * L + 13 * lp, 17)

    # Test 1: round-2 over-restrictive placeholder
    ops._so4_radial_overlap = ops._so4_radial_overlap_placeholder
    try:
        for n_max in (2, 3):
            O = ops.TruncatedOperatorSystem(n_max)
            in_I, _ = O.identity_in_O()
            assert in_I
            prop, dims = ops.propagation_number(O, max_k=4)
            assert prop == 2
    finally:
        ops._so4_radial_overlap = original

    # Test 2: alternate over-restrictive rational placeholder
    ops._so4_radial_overlap = alternate
    try:
        for n_max in (2, 3):
            O = ops.TruncatedOperatorSystem(n_max)
            in_I, _ = O.identity_in_O()
            assert in_I
            prop, dims = ops.propagation_number(O, max_k=4)
            assert prop == 2
    finally:
        ops._so4_radial_overlap = original


# ---------------------------------------------------------------------------
# Verify the n_max = 1 trivial case
# ---------------------------------------------------------------------------


def test_nmax_1_is_trivial():
    """At n_max = 1, dim_H = 1, so M_1(C) = C and O = C trivially."""
    O = TruncatedOperatorSystem(1)
    assert O.dim_H == 1
    assert O.envelope_dim == 1
    assert O.dim == 1
    in_I, _ = O.identity_in_O()
    assert in_I
    star_ok, _ = O.is_star_closed()
    assert star_ok
    # prop is trivially 1 here (or 0; the smallest k with O^k = M_1(C) is k=1).
    prop, dims = propagation_number(O, max_k=2)
    assert prop == 1


# ---------------------------------------------------------------------------
# Verify the SO(4) selection rules used by the placeholder are the documented
# Avery (1989) rules.
# ---------------------------------------------------------------------------


def test_selection_rules_match_avery():
    """For each (n, n', N) triple, the placeholder is nonzero iff the SO(4)
    triangle |n - n'| + 1 <= N <= n + n' - 1 holds, and L <= N - 1.
    """
    from geovac.operator_system import _so4_radial_overlap_placeholder
    for n in range(1, 5):
        for np_ in range(1, 5):
            for N in range(1, 8):
                triangle_ok = abs(n - np_) + 1 <= N <= n + np_ - 1
                # Try L=0, l=l'=0: subhood is L <= N - 1, i.e., N >= 1.
                val = _so4_radial_overlap_placeholder(n, N, np_, 0, 0, 0)
                if triangle_ok and N >= 1:  # L=0 always OK if N >= 1
                    assert val != 0, (
                        f"expected nonzero for (n,N,n')=({n},{N},{np_}), got {val}"
                    )
                else:
                    assert val == 0


# ---------------------------------------------------------------------------
# Inspection / regression on dim(O^k) sequences
# ---------------------------------------------------------------------------


def test_dim_O_squared_equals_envelope_at_small_nmax():
    """At n_max in {2, 3, 4}, dim(O^2) = N^2 (full matrix algebra)."""
    for n_max in (2, 3, 4):
        O = TruncatedOperatorSystem(n_max)
        prods = operator_system_products(O.multiplier_matrices, 2)
        dim_O2 = operator_system_dim(prods)
        assert dim_O2 == O.envelope_dim, (
            f"dim(O^2) = {dim_O2}, expected N^2 = {O.envelope_dim} at n_max={n_max}"
        )


# ---------------------------------------------------------------------------
# Check that the contains() method is well-behaved on basis elements
# ---------------------------------------------------------------------------


def test_contains_returns_True_on_basis_elements():
    """For each generator M_i in O, contains(M_i) returns True."""
    O = TruncatedOperatorSystem(2)
    for label, M in O.basis_matrices:
        in_O, residual = O.contains(M)
        assert in_O, f"basis matrix {label} not in own span (residual {residual:.3e})"
        assert residual < 1e-10
