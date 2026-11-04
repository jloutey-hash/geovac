"""Paper 24 thm:pi-free -- the N_max = 5 certificate, pinned by an independent route.

Paper 24 (``papers/group3_foundations/paper_24_bargmann_segal.tex``) prints, in
the abstract and in the "Numerical verification" paragraph after
thm:pi-free, the certificate of ``verify_pi_free(N_max=5)``:

    n_nodes: 56, n_edges: 165, irrationals_encountered: [], pi_free: True

and Sec. "Topology" re-uses 56 x 165 for the signed incidence matrix.

``tests/test_bargmann_graph.py::test_pi_free_larger`` checks only the
``pi_free`` flag and the empty irrational list, i.e. it would accept a graph
with any node or edge count, and its rationality test is an ``isinstance``
check against the same module that built the weights.  This file pins the
certificate by an INDEPENDENT construction:

  * nodes enumerated here from eq:node-count's branching rule
    (l = N, N-2, ..., same parity as N; |m| <= l), not from the module;
  * edges enumerated here from the dipole rules (Delta N = +1,
    |Delta l| = 1, |Delta m| <= 1);
  * every edge weight recomputed from eq:edge-weight with the paper's
    closed-form radial factors and an exact sympy Wigner 3j (a different
    code path from the module's hand-written closed forms), and compared
    to the module's weight as an exact Fraction equality.

Wrong answers this rejects (fire-tested with debug/qa/fire_test.py):
  * the parity filter on l dropped (all 0 <= l <= N): node count 91, not 56;
  * the Delta l = -1 branch dropped: fewer than 165 edges;
  * a float weight entering the build (``float(radial_sq * ang_sq)``):
    not a Fraction, irrational-looking entries in the certificate;
  * a wrong radial factor (n_r + l + 1/2 instead of n_r + l + 3/2): the
    exact weight comparison fails.
"""
from __future__ import annotations

from fractions import Fraction

from sympy import sympify
from sympy.physics.wigner import wigner_3j

from geovac.nuclear.bargmann_graph import build_bargmann_graph, verify_pi_free

N_MAX = 5
PAPER_NODES = 56   # Paper 24 abstract + eq:node-count + certificate
PAPER_EDGES = 165  # Paper 24 abstract + certificate + sec. Topology


def _nodes_independent(n_max: int):
    out = []
    for N in range(n_max + 1):
        for l in range(N % 2, N + 1, 2):
            for m in range(-l, l + 1):
                out.append((N, l, m))
    return out


def _edges_independent(n_max: int):
    """Dipole edges (N,l,m) -> (N+1,l',m'), |l'-l| = 1, |m'-m| <= 1."""
    nodes = set(_nodes_independent(n_max))
    edges = []
    for (N, l, m) in sorted(nodes):
        for lp in (l - 1, l + 1):
            for mp in (m - 1, m, m + 1):
                if (N + 1, lp, mp) in nodes:
                    edges.append(((N, l, m), (N + 1, lp, mp)))
    return edges


def _paper_weight(N: int, l: int, m: int, lp: int, mp: int) -> Fraction:
    """eq:edge-weight: radial^2 (paper's closed forms) x exact 3j^2 (sympy)."""
    n_r = (N - l) // 2
    if lp == l + 1:
        radial = Fraction(2 * n_r + 2 * l + 3, 2)
    else:
        radial = Fraction(n_r + 1, 1)
    q = mp - m
    tj = wigner_3j(l, 1, lp, m, q, -mp)
    tj2 = sympify(tj) ** 2
    assert tj2.is_Rational, (N, l, m, lp, mp, tj2)
    return radial * Fraction(int(tj2.p), int(tj2.q))


def test_certificate_counts_nmax5():
    """56 nodes / 165 edges / zero irrationals, against an independent count."""
    cert = verify_pi_free(N_max=N_MAX)
    nodes = _nodes_independent(N_MAX)
    edges = _edges_independent(N_MAX)

    assert len(nodes) == PAPER_NODES
    assert len(edges) == PAPER_EDGES
    assert cert["n_nodes"] == PAPER_NODES
    assert cert["n_edges"] == PAPER_EDGES
    assert cert["irrationals_encountered"] == []
    assert cert["all_diagonal_rational"] is True
    assert cert["all_adjacency_rational"] is True
    assert cert["pi_free"] is True


def test_every_weight_is_the_exact_paper_rational_nmax5():
    """The module's edge set and every weight equal eq:edge-weight exactly."""
    g = build_bargmann_graph(N_MAX)
    assert g.nodes == _nodes_independent(N_MAX)

    module = {}
    for (i, j), w in g.adjacency.items():
        # type(w) is Fraction, not merely a float equal to a rational
        assert type(w) is Fraction, (i, j, w)
        a, b = g.nodes[i], g.nodes[j]
        lo, hi = (a, b) if a[0] < b[0] else (b, a)
        module[(lo, hi)] = w

    edges = _edges_independent(N_MAX)
    assert set(module) == set(edges)
    for (lo, hi) in edges:
        (N, l, m), (_, lp, mp) = lo, hi
        assert module[(lo, hi)] == _paper_weight(N, l, m, lp, mp), (lo, hi)

    for (N, _, _), d in zip(g.nodes, g.diagonal):
        assert type(d) is Fraction and d == Fraction(2 * N + 3, 2)
