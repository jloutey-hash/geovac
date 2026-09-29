"""Paper 32 (Sec. K decomposition, "B is a finite scalar Casimir trace"): the
normalisation of B = 42.

Paper 32 states B = tr_{A_GV}(C . 1) at n_max = 3, EITHER with C = L^2 (the
SO(3) Casimir l(l+1)) OR with C = (n^2 - 1)/2, HALF the SO(4) Casimir
L^2 + A^2 = n^2 - 1; with the un-halved SO(4) Casimir the same trace is
84 = 2B.  (Paper 2 defines B by the SO(3) route; tests/test_paper2_ingredients.py
pins that route alone.)

    sum_{n<=3} n^2 (n^2 - 1)/2 = 0 + 6 + 36 = 42,
    sum_{n<=3} n^2 (n^2 - 1)   = 84.

WRONG ANSWER REJECTED: B read off the un-halved SO(4) Casimir (B = 84), or a
claim that the SO(3) and SO(4) routes agree without the factor 1/2.

The trace is taken over the node set of geovac.GeometricLattice(max_n=3)
(the production label set), not over a hand-written range.

Fire-test plant: the SO(4) route's halving removed (Fraction(n*n - 1, 2) ->
Fraction(n*n - 1, 1)) -> FIRES.
"""
from __future__ import annotations

from fractions import Fraction

from geovac import GeometricLattice


def _nodes(n_max: int):
    return list(GeometricLattice(max_n=n_max).states)


def _trace_so3(n_max: int) -> Fraction:
    return sum((Fraction(l * (l + 1)) for (n, l, m) in _nodes(n_max)), Fraction(0))


def _trace_half_so4(n_max: int) -> Fraction:
    return sum((Fraction(n * n - 1, 2) for (n, l, m) in _nodes(n_max)), Fraction(0))


def _trace_so4(n_max: int) -> Fraction:
    return sum((Fraction(n * n - 1) for (n, l, m) in _nodes(n_max)), Fraction(0))


def test_node_set_is_the_n_max_3_label_set():
    nodes = _nodes(3)
    assert len(nodes) == 14 and len(set(nodes)) == 14


def test_B_is_42_by_both_routes():
    """SO(3) route and half-SO(4) route both give 42 at n_max = 3."""
    assert _trace_so3(3) == 42
    assert _trace_half_so4(3) == 42


def test_unhalved_so4_trace_is_84_equals_2B():
    """The un-halved SO(4) Casimir trace is 84 = 2B -- NOT B."""
    b = _trace_half_so4(3)
    assert _trace_so4(3) == 84
    assert _trace_so4(3) == 2 * b
    assert _trace_so4(3) != 42


def test_halved_routes_agree_at_every_cutoff():
    """Per shell: sum_l (2l+1) l(l+1) = n^2 (n^2 - 1)/2 exactly, so the two
    routes agree at every cutoff (6, 42, 162, 462 at n_max = 2..5), and the
    un-halved SO(4) trace is exactly twice that everywhere."""
    for n_max, expected in [(2, 6), (3, 42), (4, 162), (5, 462)]:
        assert _trace_so3(n_max) == expected == _trace_half_so4(n_max)
        assert _trace_so4(n_max) == 2 * expected
