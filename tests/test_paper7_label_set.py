"""Paper 7 abstract / item 1 (and Paper 0): the production lattice's node set
IS the Gel'fand--Tsetlin label set {(n, l, m) : 1 <= n <= N, 0 <= l < n,
|m| <= l} of the Peter--Weyl truncation of L^2(S^3), with n^2 labels per
shell.  Until 2026-09-30 (/qa trunk FULL #10) no test asserted this against
`GeometricLattice.states`; the Casimir test builds its own basis.

WRONG ANSWERS REJECTED: a lattice missing a label, duplicating one, or
admitting l = n (n^2 + (2n+1) per shell), or a shell count other than n^2.
"""
from __future__ import annotations

import pytest

from geovac.lattice import GeometricLattice


def _gt_labels(N: int):
    return {(n, l, m) for n in range(1, N + 1) for l in range(n) for m in range(-l, l + 1)}


@pytest.mark.parametrize("N", [1, 2, 3, 5, 8])
def test_node_set_is_gelfand_tsetlin_label_set(N):
    states = list(GeometricLattice(max_n=N).states)
    assert len(states) == len(set(states)), "duplicate labels in the lattice"
    assert set(states) == _gt_labels(N)
    for n in range(1, N + 1):
        assert sum(1 for s in states if s[0] == n) == n * n
