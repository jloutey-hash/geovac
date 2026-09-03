"""Paper 0 sec:VI -- the spin-doubled vertex count |V| = N(N+1)(2N+1)/3 for
all N (earns the [SYMBOLIC PROOF] tag the paper carries; the enumeration in
tests/test_dirac_lattice.py covers N <= 6 only).  Added 2026-09-03 (trunk
FULL run #3, carryforward I.1.17)."""
from __future__ import annotations

import sympy as sp

from geovac.lattice import GeometricLattice


def test_vertex_count_closed_form_symbolic():
    N, n = sp.symbols("N n", positive=True, integer=True)
    spinless = sp.summation(n ** 2, (n, 1, N))          # |V| of the (n, l, m) lattice
    doubled = sp.summation(2 * n ** 2, (n, 1, N))       # spin-doubled
    assert sp.simplify(spinless - N * (N + 1) * (2 * N + 1) / 6) == 0
    assert sp.simplify(doubled - N * (N + 1) * (2 * N + 1) / 3) == 0


def test_vertex_count_matches_production_lattice():
    for N in (1, 2, 3, 5, 8, 12):
        assert GeometricLattice(N).num_states == N * (N + 1) * (2 * N + 1) // 6
