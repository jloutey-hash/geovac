"""
Paper 25 Proposition 1, item 4: the Hopf-quotient S^2 base graph.

Collapsing the magnetic fibers {(n, l, m_l) : m_l = -l..l} of the GeoVac
Hopf graph at n_max = 3 to single (n, l) sectors yields a 6-sector base
graph whose inter-sector edge weights count the S^3 edges crossing
between sectors. Paper 25 (Sec. III.B and Prop. 1 item 4, Track alpha-D)
states its Laplacian spectrum is exactly

    {0, 0, 0, 1, 3, 6}

with three zero modes reflecting the three connected components
(s-chain, p-block, d-chain) inherited from the S^3 graph.

This test builds the quotient live from GeometricLattice(max_n=3) and
pins the spectrum in EXACT sympy arithmetic (group5 cert finding E10:
this spectrum previously had no test).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from geovac.lattice import GeometricLattice


def _quotient_laplacian(max_n: int):
    """Weighted (n,l)-sector Laplacian of the Hopf quotient.

    Sector weight W[s, s'] = number of undirected S^3 edges crossing
    between the m_l-fibers of s and s' (Paper 25 Sec. III.B: 'inter-sector
    edges of the base graph count the number of S^3 edges crossing
    between sectors'; intra-sector angular edges collapse into fiber
    structure and do not appear).
    """
    lat = GeometricLattice(max_n=max_n)
    A = lat.adjacency.toarray()
    states = lat.states
    sectors = sorted(set((n, l) for (n, l, m) in states))
    sidx = {s: i for i, s in enumerate(sectors)}
    k = len(sectors)
    W = sp.zeros(k, k)
    V = len(states)
    for u in range(V):
        for v in range(u + 1, V):
            if A[u, v]:
                su = (states[u][0], states[u][1])
                sv = (states[v][0], states[v][1])
                if su != sv:
                    W[sidx[su], sidx[sv]] += 1
                    W[sidx[sv], sidx[su]] += 1
    deg = sp.diag(*[sum(W[i, j] for j in range(k)) for i in range(k)])
    return sectors, W, deg - W


class TestS2QuotientNmax3:
    def test_sectors_and_fiber_dims(self):
        """Six sectors (Paper 25 Eq. s2_sectors) with fiber dims 2l+1."""
        sectors, _, _ = _quotient_laplacian(3)
        assert sectors == [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
        assert [2 * l + 1 for (n, l) in sectors] == [1, 1, 3, 1, 3, 5]

    def test_quotient_weights(self):
        """Crossing-count weights: s-chain radial edges have weight 1;
        the p-block radial bundle (2,1)-(3,1) has weight 3 (one edge per
        m_l = -1, 0, +1); the d-chain sector (3,2) is isolated."""
        sectors, W, _ = _quotient_laplacian(3)
        s = {sec: i for i, sec in enumerate(sectors)}
        assert W[s[(1, 0)], s[(2, 0)]] == 1
        assert W[s[(2, 0)], s[(3, 0)]] == 1
        assert W[s[(2, 1)], s[(3, 1)]] == 3
        # (3,2) isolated:
        assert all(W[s[(3, 2)], j] == 0 for j in range(6))
        # no other couplings
        assert sum(W[i, j] for i in range(6) for j in range(6)) == 2 * (1 + 1 + 3)

    def test_quotient_spectrum_exact(self):
        """Laplacian spectrum is exactly {0, 0, 0, 1, 3, 6} (Prop 1 item 4).

        Exact sympy eigenvalues; the three zero modes are the three
        connected components (s-chain, p-block, d-chain)."""
        _, _, L = _quotient_laplacian(3)
        eigs = L.eigenvals()  # dict {eigenvalue: multiplicity}
        assert eigs == {sp.Integer(0): 3, sp.Integer(1): 1,
                        sp.Integer(3): 1, sp.Integer(6): 1}

    def test_quotient_charpoly_factors(self):
        """Characteristic polynomial x^3 (x-1)(x-3)(x-6) -- integer roots."""
        _, _, L = _quotient_laplacian(3)
        x = sp.Symbol('x')
        p = L.charpoly(x).as_expr()
        expected = x ** 3 * (x - 1) * (x - 3) * (x - 6)
        assert sp.expand(p - expected) == 0

    def test_numeric_cross_check(self):
        """Float eigenvalues agree with the exact ones."""
        _, _, L = _quotient_laplacian(3)
        Lnp = np.array(L.tolist(), dtype=float)
        evals = np.sort(np.linalg.eigvalsh(Lnp))
        assert np.allclose(evals, [0, 0, 0, 1, 3, 6], atol=1e-12)


def test_full_graph_edge_laplacian_golden_spectrum():
    """Paper 25 Eq. (eq:L1_spectrum): the edge Laplacian L1 = B^T B of
    the FULL n_max=3 Hopf graph (V=14, E=13) has the exact golden-ratio
    spectrum {0 x2, 3/2-sqrt5/2, 1 x2, 5/2-sqrt5/2, 2, 3/2+sqrt5/2,
    3 x3, 5/2+sqrt5/2, 5}.  Added 2026-07-04 (certifying /qa run:
    this equation previously had no test and was cited only to an
    internal memo)."""
    import sympy as sp
    lat = GeometricLattice(max_n=3)
    A = np.asarray(lat.adjacency.toarray())
    V = A.shape[0]
    edges = [(i, j) for i in range(V) for j in range(i + 1, V) if A[i, j]]
    assert (V, len(edges)) == (14, 13)
    B = sp.zeros(V, len(edges))
    for k, (i, j) in enumerate(edges):
        B[i, k] = -1
        B[j, k] = 1
    L1 = B.T * B
    x = sp.symbols("x")
    s5 = sp.sqrt(5)
    expected = {
        sp.Integer(0): 2,
        sp.Rational(3, 2) - s5 / 2: 1,
        sp.Integer(1): 2,
        sp.Rational(5, 2) - s5 / 2: 1,
        sp.Integer(2): 1,
        s5 / 2 + sp.Rational(3, 2): 1,
        sp.Integer(3): 3,
        s5 / 2 + sp.Rational(5, 2): 1,
        sp.Integer(5): 1,
    }
    assert sum(expected.values()) == 13
    poly_expected = sp.expand(sp.prod([(x - k) ** v for k, v in expected.items()]))
    assert sp.expand(L1.charpoly(x).as_expr() - poly_expected) == 0
    # and beta_1 = dim ker L1 = 2 (the two p-block 4-cycles)
    assert expected[sp.Integer(0)] == 2
