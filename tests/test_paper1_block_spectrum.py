"""Papers 0/1/7 -- the GeoVac lattice Laplacian's spectrum is CLOSED FORM, and
the s/p "splitting" is not a robust spectral quantity.

Added 2026-09-03 (trunk DELTA #3, carryforward I.9).  The lattice's edges change
n (T+-) or m (L+-) but never l, so L is block diagonal in l and each block is a
GRID GRAPH: within the l-block the nodes are (n, m) with n = l+1..n_max and
m = -l..l, radial edges along n and azimuthal edges along m, i.e.

    L_l  =  Laplacian of  P_{n_max - l}  x  P_{2l+1}   (Cartesian product),

whose spectrum is the sum of the two path spectra:

    spec(L_l) = { (2 - 2 cos(j pi/(n_max-l))) + (2 - 2 cos(k pi/(2l+1))) }.

Consequences pinned here:

1. The WHOLE spectrum of L is closed form (verified against dense eigh).
2. lambda_max -> 2 d_max = 8 with a PROVEN rate: the deficit is
   8 - lambda_max = (42.6 + o(1)) / n_max^2  (log-log slope -1.984 over
   n_max = 20..320).  The papers' "no rate proven here" understates this.
3. The l = 0 block is the path P_{n_max}: its top eigenvalue is exactly
   2 - 2 cos((n_max-1) pi/n_max) -> 4.
4. **The s/p splitting is ill-conditioned.**  The l = 0 and l = 1 blocks are
   DISCONNECTED components, so there is no 2s/2p degeneracy in the graph
   spectrum to lift; the quantity Paper 1 reports selects, for each of the
   nodes (2,0,0) and (2,1,0), the eigenvector of maximal amplitude there --
   and that selection is a near-tie (top-two |overlap| within 1.2% at
   n_max = 30, the selected eigenvalues differing by 3.0 vs 0.011).  The
   twelve percentages Paper 1 prints are therefore reproducible but not
   robust, and this test records that fact so it cannot be forgotten.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import diags

from geovac.lattice import GeometricLattice


def block_spectrum(n_max: int, l: int) -> np.ndarray:
    """Closed-form spectrum of the l-block = P_{n_max-l} x P_{2l+1}."""
    a = np.array([2 - 2 * np.cos(j * np.pi / (n_max - l)) for j in range(n_max - l)])
    b = np.array([2 - 2 * np.cos(k * np.pi / (2 * l + 1)) for k in range(2 * l + 1)])
    return (a[:, None] + b[None, :]).ravel()


def _dense_laplacian(n_max: int):
    lat = GeometricLattice(n_max)
    A = lat.adjacency.tocsr()
    deg = np.array(A.sum(axis=1)).ravel()
    return lat, (diags(deg) - A).toarray()


@pytest.mark.parametrize("n_max", [4, 7, 10])
def test_closed_form_reproduces_the_whole_spectrum(n_max):
    _, L = _dense_laplacian(n_max)
    direct = np.sort(np.linalg.eigvalsh(L))
    closed = np.sort(np.concatenate([block_spectrum(n_max, l) for l in range(n_max)]))
    assert direct.shape == closed.shape
    assert np.abs(direct - closed).max() < 1e-9


def test_closed_form_wrong_product_is_rejected():
    """Guard: the grid structure is what makes the closed form work -- a wrong
    factorisation (P_{n_max} x P_{2l+1}) must NOT reproduce the spectrum."""
    n_max = 7
    _, L = _dense_laplacian(n_max)
    direct = np.sort(np.linalg.eigvalsh(L))
    wrong = []
    for l in range(n_max):
        a = np.array([2 - 2 * np.cos(j * np.pi / n_max) for j in range(n_max)])
        b = np.array([2 - 2 * np.cos(k * np.pi / (2 * l + 1)) for k in range(2 * l + 1)])
        wrong.append((a[:, None] + b[None, :]).ravel())
    wrong = np.sort(np.concatenate(wrong))
    assert wrong.shape != direct.shape or np.abs(direct - wrong).max() > 1e-3


def test_saturation_rate_is_order_n_squared():
    """8 - lambda_max = (42.6 + o(1)) / n_max^2 -- a proven rate, from the
    closed form (the papers say 'no rate proven here')."""
    ns = [20, 40, 80, 160, 320]
    gaps = [8 - max(block_spectrum(n, l).max() for l in range(n)) for n in ns]
    slope = np.polyfit(np.log(ns), np.log(gaps), 1)[0]
    assert -2.02 < slope < -1.95, slope
    assert 42.0 < gaps[-1] * ns[-1] ** 2 < 43.5


def test_s_wave_block_is_a_path_graph():
    for n_max in (10, 30):
        top = block_spectrum(n_max, 0).max()
        assert abs(top - (2 - 2 * np.cos((n_max - 1) * np.pi / n_max))) < 1e-12
        assert top < 4.0


def test_sp_splitting_identification_is_ill_conditioned():
    """The mode selection behind Paper 1's s/p percentages is a near-tie: at
    n_max = 30 the runner-up eigenvector at the (2,0,0) node is within 1.5% of
    the winner in overlap while its eigenvalue differs by more than 2.  The
    reported splitting is therefore reproducible but not robust."""
    lat, L = _dense_laplacian(30)
    idx = {s: i for i, s in enumerate(lat.states)}
    w, v = np.linalg.eigh(L)
    ov = np.abs(v[idx[(2, 0, 0)], :])
    order = np.argsort(-ov)[:2]
    margin = (ov[order[0]] - ov[order[1]]) / ov[order[0]]
    assert margin < 0.02, margin                       # near-tie
    assert abs(w[order[0]] - w[order[1]]) > 1.0        # wildly different eigenvalues


def test_l_zero_and_l_one_blocks_are_disconnected():
    """There is no 2s/2p degeneracy in the graph spectrum to lift: the two
    nodes live in different connected components."""
    lat = GeometricLattice(12)
    A = lat.adjacency.tocsr()
    st = np.array(lat.states)
    rows, cols = A.nonzero()
    assert np.all(st[rows, 1] == st[cols, 1])
    from scipy.sparse.csgraph import connected_components
    ncomp, labels = connected_components(A, directed=False)
    idx = {s: i for i, s in enumerate(lat.states)}
    assert labels[idx[(2, 0, 0)]] != labels[idx[(2, 1, 0)]]


def test_cg_construction_splitting_does_not_decay():
    """Paper 1 §III's own construction A = |T+|+|T-|+|L+|+|L-|: the s/p lift
    does NOT decay (129/68/84/139% at n_max = 8/10/15/20).  Closes the
    NO-TEST gap on the scoping sentence (DELTA #3 D3).  Same conditioning
    caveat as above applies to the identification."""
    import math
    def build_cg(n_max):
        states = [(n, l, m) for n in range(1, n_max + 1) for l in range(n) for m in range(-l, l + 1)]
        idx = {s: i for i, s in enumerate(states)}
        A = np.zeros((len(states), len(states)))
        for (n, l, m) in states:
            i = idx[(n, l, m)]
            if m < l:
                j = idx[(n, l, m + 1)]; a = math.sqrt((l - m) * (l + m + 1))
                A[i, j] += a; A[j, i] += a
            if (n + 1, l, m) in idx:
                j = idx[(n + 1, l, m)]; a = math.sqrt((n - l) * (n + l + 1) / 4.0)
                A[i, j] += a; A[j, i] += a
        return states, idx, A
    vals = {}
    for n_max in (8, 10, 15, 20):
        states, idx, A = build_cg(n_max)
        L = np.diag(A.sum(axis=1)) - A
        w, v = np.linalg.eigh(L)
        j2s = int(np.argmax(np.abs(v[idx[(2, 0, 0)], :])))
        j2p = int(np.argmax(np.abs(v[idx[(2, 1, 0)], :])))
        vals[n_max] = abs(w[j2p] - w[j2s]) / abs(w[j2s])
    assert vals[8] > 1.0 and vals[20] > 1.0            # >100%, no decay
    assert vals[20] > vals[10]                          # not monotone downward
    for n_max, expected in ((8, 1.2931), (10, 0.6817), (15, 0.8371), (20, 1.3878)):
        assert abs(vals[n_max] - expected) < 5e-3, (n_max, vals[n_max])
