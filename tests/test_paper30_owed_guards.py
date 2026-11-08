"""Owed guards for Paper 30 (papers/group5_qed_gauge/paper_30_su2_wilson.tex),
written as a separate step after the 2026-11-07 merge (CLAUDE.md Sec. 9).

Each guard computes its expected value by a route that does not share the code
path under test:

  * Z(beta) and <W>: seeded Monte Carlo over the RAW link variables of the
    actual label graph (no gauge fixing, no factorization assumed), with an
    SU(2) product written here in (a, b) form, against the closed forms of
    prop:factorization and Sec. 4.1;
  * rank Delta_1: exact integer ranks, with the plaquettes taken from the
    production enumerator and checked against the unit squares of the grids
    (eq:grid_decomposition), built here from the labels;
  * the cubic term: the BCH triple-product coefficient
    c_3 = -(1/8) sum_k sum_{i<j} A_k . (A_i x A_j), derived by hand, against the
    odd part of the production Wilson action;
  * the Rule-B 2-cycles: K_{2,3} subgraphs counted from common neighbours,
    against closed triples found from the plaquettes' signed incidence.

No debug/ import.
"""
from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components
from scipy.special import iv

from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.lattice import GeometricLattice
from geovac.su2_wilson_gauge import (
    diagonal_su2_from_phase,
    enumerate_plaquettes,
    expectation_wilson_loop,
    partition_function_character_expansion,
    u1_action_from_su2,
    wilson_action,
)


# ---------------------------------------------------------------------------
# Independent scalar-graph construction (eq:grid_decomposition)
# ---------------------------------------------------------------------------

def _grid(n_max):
    """States, index map, sorted edge list and unit squares of the scalar
    label graph, built from the labels alone: radial ladders (n, l, m) ~
    (n+1, l, m) and angular ladders (n, l, m) ~ (n, l, m+1).  Each square is
    a cyclic vertex 4-tuple."""
    states = [(n, l, m) for n in range(1, n_max + 1) for l in range(n)
              for m in range(-l, l + 1)]
    idx = {s: i for i, s in enumerate(states)}
    edges = set()
    squares = []
    for (n, l, m) in states:
        if (n + 1, l, m) in idx:
            edges.add(tuple(sorted((idx[(n, l, m)], idx[(n + 1, l, m)]))))
        if (n, l, m + 1) in idx:
            edges.add(tuple(sorted((idx[(n, l, m)], idx[(n, l, m + 1)]))))
        cyc = [(n, l, m), (n, l, m + 1), (n + 1, l, m + 1), (n + 1, l, m)]
        if all(c in idx for c in cyc):
            squares.append(tuple(idx[c] for c in cyc))
    return states, idx, sorted(edges), squares


def _lattice_matches_grid(n_max):
    lat = GeometricLattice(max_n=n_max)
    states, idx, edges, squares = _grid(n_max)
    assert list(lat.states) == states
    A = lat.adjacency.toarray()
    prod_edges = sorted({(min(i, j), max(i, j)) for i, j in zip(*np.nonzero(A))})
    assert prod_edges == edges, "production scalar graph is not the grid union"
    c = int(connected_components(lat.adjacency, directed=False)[0])
    return lat, A, states, idx, edges, squares, c


def _cycle_traversals(cycle, eidx):
    """[(edge index, +1/-1)] for a closed vertex cycle."""
    out = []
    L = len(cycle)
    for k in range(L):
        u, v = cycle[k], cycle[(k + 1) % L]
        if (u, v) in eidx:
            out.append((eidx[(u, v)], 1))
        else:
            out.append((eidx[(v, u)], -1))
    return out


# ---------------------------------------------------------------------------
# SU(2) in (a, b) form: U = a I + i b . sigma, a^2 + |b|^2 = 1
# ---------------------------------------------------------------------------

def _su2_mul(x, y):
    a1, b1 = x[..., 0], x[..., 1:]
    a2, b2 = y[..., 0], y[..., 1:]
    a = a1 * a2 - np.sum(b1 * b2, axis=-1)
    b = a1[..., None] * b2 + a2[..., None] * b1 - np.cross(b1, b2)
    return np.concatenate([a[..., None], b], axis=-1)


def _su2_inv(x):
    y = x.copy()
    y[..., 1:] *= -1.0
    return y


def _haar(rng, shape):
    q = rng.standard_normal(shape + (4,))
    return q / np.linalg.norm(q, axis=-1, keepdims=True)


def _holonomy_halftrace(links, trav):
    """(1/2) tr of the ordered product along `trav`; links shape (N, E, 4)."""
    U = None
    for e, s in trav:
        X = links[:, e, :] if s > 0 else _su2_inv(links[:, e, :])
        U = X if U is None else _su2_mul(U, X)
    return U[:, 0]


def test_su2_ab_product_matches_matrices():
    """Sanity of the (a, b) product used by the MC guards (not a paper claim).
    Rejects the quaternion sign convention (+ b1 x b2), which would not
    represent U = a + i b.sigma."""
    sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
           np.array([[1, 0], [0, -1]], complex)]
    rng = np.random.default_rng(0)
    x, y = _haar(rng, (2,))
    M = lambda q: q[0] * np.eye(2) + 1j * sum(q[k + 1] * sig[k] for k in range(3))
    assert np.allclose(M(_su2_mul(x, y)), M(x) @ M(y), atol=1e-14)
    assert np.allclose(M(_su2_inv(x)), M(x).conj().T, atol=1e-14)


def _f_su2(beta):
    """One-plaquette SU(2) integral e^{-beta} d_0(beta) = 2 e^{-beta} I_1/beta."""
    return 2.0 * np.exp(-beta) * iv(1, beta) / beta


def _f_u1(beta):
    """One-plaquette U(1) integral e^{-beta} I_0(beta)."""
    return np.exp(-beta) * iv(0, beta)


def _mc_mean(sampler, n_batches):
    vals = np.concatenate([sampler() for _ in range(n_batches)])
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(len(vals)))


# ---------------------------------------------------------------------------
# Guard 1: Z(beta) = [2 e^{-beta} I_1(beta)/beta]^{beta_1}  (eq:Z_exact)
# ---------------------------------------------------------------------------

Z_CASES = [(2, 1.0, 1), (3, 1.0, 4), (4, 0.7, 8), (5, 0.4, 8)]


@pytest.mark.parametrize("n_max,beta,n_batches", Z_CASES)
def test_su2_partition_function_closed_form(n_max, beta, n_batches):
    """Paper 30 prop:factorization, eq:Z_exact, at n_max = 2..5.

    Route: Z = E_Haar[exp(-S_W)] over all raw links of the actual graph,
    seeded MC (25k samples per batch), plaquettes = unit squares built from the
    labels.  Also checks the production partition_function_character_expansion
    on the production 4-cycles against the closed form.

    Rejects: the exponent E or V in place of beta_1; the U(1) Bessel form
    [e^{-beta} I_0]^{beta_1}; plaquettes = every closed walk up to length 8."""
    lat, A, states, idx, edges, squares, c = _lattice_matches_grid(n_max)
    V, E = len(states), len(edges)
    beta1 = E - V + c
    assert len(squares) == beta1
    closed = _f_su2(beta) ** beta1
    # production: the 4-cycles, through the production character expansion
    prod_plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
    Z_prod = partition_function_character_expansion(prod_plaqs, beta)
    assert abs(Z_prod / closed - 1.0) < 1e-12, (Z_prod, closed)
    eidx = {e: k for k, e in enumerate(edges)}
    travs = [_cycle_traversals(sq, eidx) for sq in squares]
    rng = np.random.default_rng(3000 + n_max)

    def sampler():
        links = _haar(rng, (25_000, E))
        S = np.zeros(25_000)
        for t in travs:
            S += beta * (1.0 - _holonomy_halftrace(links, t))
        return np.exp(-S)

    Z_mc, se = _mc_mean(sampler, n_batches)
    if beta1 == 0:
        assert Z_mc == 1.0 and closed == 1.0
        assert _f_su2(beta) ** E != 1.0          # exponent E would be wrong here
        return
    assert abs(Z_mc - closed) < 4.0 * se, (Z_mc, closed, se)
    for wrong in (_f_su2(beta) ** E, _f_su2(beta) ** V, _f_u1(beta) ** beta1):
        assert abs(Z_mc - wrong) > 8.0 * se, (wrong, Z_mc, se)


# ---------------------------------------------------------------------------
# Guard 2: the abelian forms (Sec. 4.1)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max,beta", [(3, 1.0), (4, 0.7)])
def test_u1_partition_function_closed_form(n_max, beta):
    """Paper 30 Sec. 4.1: Z = [e^{-beta} I_0(beta)]^{beta_1} on the scalar graph.

    Route: seeded MC over uniform raw link angles theta_e, theta_P the signed
    circulation around the label-built unit squares.

    Rejects: the SU(2) form [2 e^{-beta} I_1/beta]^{beta_1}; the exponent E."""
    lat, A, states, idx, edges, squares, c = _lattice_matches_grid(n_max)
    E = len(edges)
    beta1 = E - len(states) + c
    eidx = {e: k for k, e in enumerate(edges)}
    travs = [_cycle_traversals(sq, eidx) for sq in squares]
    rng = np.random.default_rng(4000 + n_max)

    def sampler():
        th = rng.uniform(-np.pi, np.pi, (100_000, E))
        S = np.zeros(100_000)
        for t in travs:
            S += beta * (1.0 - np.cos(sum(s * th[:, e] for e, s in t)))
        return np.exp(-S)

    Z_mc, se = _mc_mean(sampler, 4)
    closed = _f_u1(beta) ** beta1
    assert abs(Z_mc - closed) < 4.0 * se, (Z_mc, closed, se)
    for wrong in (_f_su2(beta) ** beta1, _f_u1(beta) ** E):
        assert abs(Z_mc - wrong) > 8.0 * se, (wrong, Z_mc, se)


def _ratio_mc(weights_and_obs, n_batches):
    """Reweighted estimator <O> = E[O w]/E[w] with a delta-method stderr."""
    ws, os_ = [], []
    for _ in range(n_batches):
        w, o = weights_and_obs()
        ws.append(w)
        os_.append(o)
    w, o = np.concatenate(ws), np.concatenate(os_)
    R = float(np.sum(o * w) / np.sum(w))
    se = float(np.std(o * w - R * w, ddof=1) / (np.sqrt(len(w)) * w.mean()))
    return R, se


# The loops on the actual lattice.  n_max = 3: P_A (area 1) and the six-edge
# boundary of P_A u P_B (area 2), eq:p_block_plaquettes.  n_max = 4: in the
# l = 1 block (a 3 x 3 grid) the eight-edge boundary of the 2 x 2 square
# (area 4); in the l = 2 block (a 2 x 5 strip) the eight-edge boundary of three
# adjacent squares (area 3).
LOOPS = {
    (3, 1): [(2, 1, -1), (2, 1, 0), (3, 1, 0), (3, 1, -1)],
    (3, 2): [(2, 1, -1), (2, 1, 0), (2, 1, 1), (3, 1, 1), (3, 1, 0), (3, 1, -1)],
    (4, 4): [(2, 1, -1), (2, 1, 0), (2, 1, 1), (3, 1, 1), (4, 1, 1), (4, 1, 0),
             (4, 1, -1), (3, 1, -1)],
    (4, 3): [(3, 2, -2), (3, 2, -1), (3, 2, 0), (3, 2, 1), (4, 2, 1), (4, 2, 0),
             (4, 2, -1), (4, 2, -2)],
}
LOOP_L = {(3, 1): 1, (3, 2): 1, (4, 4): 1, (4, 3): 2}


def _loop_setup(n_max, area):
    """Edges, plaquette traversals and loop traversal restricted to the l-block
    that contains the loop (blocks are disconnected, so this is exact)."""
    lat, A, states, idx, edges, squares, c = _lattice_matches_grid(n_max)
    l = LOOP_L[(n_max, area)]
    keep = [e for e in edges if states[e[0]][1] == l]
    eidx = {e: k for k, e in enumerate(keep)}
    sqs = [sq for sq in squares if states[sq[0]][1] == l]
    travs = [_cycle_traversals(sq, eidx) for sq in sqs]
    loop = [idx[s] for s in LOOPS[(n_max, area)]]
    for k in range(len(loop)):
        assert tuple(sorted((loop[k], loop[(k + 1) % len(loop)]))) in eidx
    # the loop's edge vector is the signed sum of exactly `area` squares
    cvec = np.zeros(len(keep))
    for e, s in _cycle_traversals(loop, eidx):
        cvec[e] += s
    Bp = np.zeros((len(keep), len(sqs)))
    for j, t in enumerate(travs):
        for e, s in t:
            Bp[e, j] += s
    coef = np.linalg.lstsq(Bp, cvec, rcond=None)[0]
    assert np.allclose(Bp @ coef, cvec)
    assert int(np.sum(np.abs(np.round(coef)))) == area
    return keep, eidx, travs, _cycle_traversals(loop, eidx), sqs, idx


SU2_LOOP_CASES = [(3, 1, 2.0), (3, 2, 2.0), (3, 2, 5.0), (4, 4, 2.0), (4, 3, 2.0)]


@pytest.mark.parametrize("n_max,area,beta", SU2_LOOP_CASES)
def test_su2_area_law_on_lattice_loops(n_max, area, beta):
    """Paper 30 prop:factorization, eq:W_exact: <(1/2) tr U_C> = (I_2/I_1)^A for
    loops of the actual lattice (P_A; the six-edge loop around P_A u P_B; at
    n_max = 4 an area-4 square loop and an area-3 strip loop).

    Route: reweighted seeded MC over raw Haar links of the loop's l-block,
    observable = (1/2) tr of the loop's ordered product (no Schur identity, no
    hand-picked matrices).  For P_A also checks the production
    expectation_wilson_loop.

    Rejects: the U(1) ratio (I_1/I_0)^A; a perimeter law (I_2/I_1)^{L/4};
    the wrong area; the withdrawn MC column (0.332 at beta = 2 for P_A)."""
    keep, eidx, travs, ltrav, sqs, idx = _loop_setup(n_max, area)
    E = len(keep)
    rng = np.random.default_rng(5000 + 10 * n_max + area + int(beta))

    def draw():
        links = _haar(rng, (50_000, E))
        S = np.zeros(50_000)
        for t in travs:
            S += beta * (1.0 - _holonomy_halftrace(links, t))
        return np.exp(-S), _holonomy_halftrace(links, ltrav)

    W, se = _ratio_mc(draw, 8)
    expected = (iv(2, beta) / iv(1, beta)) ** area
    assert abs(W - expected) < 4.0 * se, (W, expected, se)
    L = len(ltrav)
    wrongs = [(iv(1, beta) / iv(0, beta)) ** area,
              (iv(2, beta) / iv(1, beta)) ** (area + 1),
              (iv(2, beta) / iv(1, beta)) ** max(area - 1, 0)]
    if L / 4 != area:
        wrongs.append((iv(2, beta) / iv(1, beta)) ** (L / 4))
    if (n_max, area, beta) == (3, 1, 2.0):
        wrongs.append(0.332)
        lat = GeometricLattice(max_n=3)
        P_A = enumerate_plaquettes(lat.adjacency.toarray(), max_length=4,
                                   both_orientations=False)[0]
        assert abs(expectation_wilson_loop(P_A, [P_A], beta) - W) < 4.0 * se
    for wrong in wrongs:
        assert abs(W - wrong) > 6.0 * se, (wrong, W, se)


@pytest.mark.parametrize("n_max,area,beta", [(3, 1, 1.0), (3, 2, 2.0), (4, 4, 2.0)])
def test_u1_area_law_on_lattice_loops(n_max, area, beta):
    """Paper 30 Sec. 4.1: <cos theta_C> = (I_1/I_0)^A on the scalar graph.

    Route: reweighted seeded MC over uniform raw link angles of the loop's
    l-block, observable cos of the loop's signed circulation.

    Rejects: the SU(2) ratio (I_2/I_1)^A; the wrong area."""
    keep, eidx, travs, ltrav, sqs, idx = _loop_setup(n_max, area)
    E = len(keep)
    rng = np.random.default_rng(6000 + 10 * n_max + area)

    def draw():
        th = rng.uniform(-np.pi, np.pi, (100_000, E))
        S = np.zeros(100_000)
        for t in travs:
            S += beta * (1.0 - np.cos(sum(s * th[:, e] for e, s in t)))
        return np.exp(-S), np.cos(sum(s * th[:, e] for e, s in ltrav))

    W, se = _ratio_mc(draw, 4)
    expected = (iv(1, beta) / iv(0, beta)) ** area
    assert abs(W - expected) < 4.0 * se, (W, expected, se)
    for wrong in ((iv(2, beta) / iv(1, beta)) ** area,
                  (iv(1, beta) / iv(0, beta)) ** (area + 1)):
        assert abs(W - wrong) > 6.0 * se, (wrong, W, se)


# ---------------------------------------------------------------------------
# Guard 3: rank Delta_1 = E; Table 1 (tab:scalar_counts)
# ---------------------------------------------------------------------------

def _incidences(A, edges, plaqs):
    eidx = {e: k for k, e in enumerate(edges)}
    V = A.shape[0]
    B = np.zeros((V, len(edges)), dtype=np.int64)
    for k, (a, b) in enumerate(edges):
        B[a, k], B[b, k] = -1, 1
    B2 = np.zeros((len(edges), len(plaqs)), dtype=np.int64)
    for j, P in enumerate(plaqs):
        for e in P:
            u, v = e.source, e.target
            if (u, v) in eidx:
                B2[eidx[(u, v)], j] += 1
            else:
                B2[eidx[(v, u)], j] -= 1
    return B, B2


@pytest.mark.parametrize("n_max", [3, 4, 5, 6])
def test_hodge1_full_rank_on_scalar_graph(n_max):
    """Paper 30 prop:L1_kinetic + eq:h1_defect: Delta_1 = B^T B + B_2 B_2^T has
    full rank E, i.e. dim ker Delta_1 = beta_1 - rank B_2 = 0, at n_max = 3..6
    (the paper prints the n_max = 3 instance; this is the general claim).

    Route: exact integer ranks; B_2 from the production 4-cycles, whose edge
    sets are checked to be exactly the label-built unit squares.

    Rejects: Delta_1 = L_1 alone (rank E - beta_1); a nonzero cohomology
    defect (plaquettes not spanning the cycle space)."""
    lat, A, states, idx, edges, squares, c = _lattice_matches_grid(n_max)
    E, V = len(edges), len(states)
    beta1 = E - V + c
    plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
    prod_sets = sorted(sorted(tuple(sorted((e.source, e.target))) for e in P)
                       for P in plaqs)
    grid_sets = sorted(sorted(tuple(sorted((sq[k], sq[(k + 1) % 4])))
                              for k in range(4)) for sq in squares)
    assert prod_sets == grid_sets
    B, B2 = _incidences(A, edges, plaqs)
    assert not np.any(B @ B2)                            # B B_2 = 0
    L1 = B.T @ B
    Delta1 = L1 + B2 @ B2.T
    rank_B2 = np.linalg.matrix_rank(B2.astype(float))
    assert rank_B2 == beta1 == len(plaqs)
    assert beta1 - rank_B2 == 0                          # eq:h1_defect
    assert np.linalg.matrix_rank(Delta1.astype(float)) == E
    assert np.linalg.matrix_rank(L1.astype(float)) == E - beta1
    if beta1 > 0:
        assert np.linalg.matrix_rank(L1.astype(float)) < E


SCALAR_TABLE = {2: (5, 3, 2, 0, 0), 3: (14, 13, 3, 2, 2),
                4: (30, 34, 4, 8, 8), 5: (55, 70, 5, 20, 20)}


@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_scalar_table_rows(n_max):
    """Paper 30 tab:scalar_counts, every row: (V, E, c, beta_1, N_4), with
    beta_1 = sum_l 2l(n_max - l - 1) (eq:betti_scalar) and N_4 = beta_1.

    Route: V, E, c from the production lattice; N_4 from the production
    enumerator; the closed form summed independently.

    Rejects: N_4 != beta_1 (4-cycles not a cycle basis); a mis-printed row."""
    lat, A, states, idx, edges, squares, c = _lattice_matches_grid(n_max)
    V, E = len(states), len(edges)
    N4 = sum(1 for P in enumerate_plaquettes(A, max_length=4, both_orientations=False)
             if len(P) == 4)
    beta1 = E - V + c
    assert beta1 == sum(2 * l * (n_max - l - 1) for l in range(n_max))
    assert (V, E, c, beta1, N4) == SCALAR_TABLE[n_max]


# ---------------------------------------------------------------------------
# Guard 4: the cubic term (Sec. 3.4, after eq:plaquette_quadratic)
# ---------------------------------------------------------------------------

_SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
        np.array([[1, 0], [0, -1]], complex)]


def _c3_bch(Atrav):
    """BCH: log U_P = eps Z_1 + eps^2 Z_2 + ..., Z_2 = (1/2) sum_{i<j} [X_i, X_j],
    with X = (i/2) A.sigma, [X_i, X_j] = (i/2)(-(A_i x A_j)).sigma.  For
    U = exp((i/2) v.sigma), 1 - (1/2) tr U = |v|^2/8 + O(v^4), so the eps^3
    coefficient is (2 v_1.v_2)/8 with v_1 = sum A_k, v_2 = -(1/2) sum_{i<j}
    A_i x A_j:  c_3 = -(1/8) sum_k sum_{i<j} A_k . (A_i x A_j)."""
    s = 0.0
    n = len(Atrav)
    for k in range(n):
        for i in range(n):
            for j in range(i + 1, n):
                s += float(Atrav[k] @ np.cross(Atrav[i], Atrav[j]))
    c3_expected = -s / 8.0
    return c3_expected


def _P_A():
    lat = GeometricLattice(max_n=3)
    return enumerate_plaquettes(lat.adjacency.toarray(), max_length=4,
                                both_orientations=False)[0]


def _links_for_traversal_fields(P, Atrav, eps):
    """Forward-keyed links such that the k-th traversal of P carries
    exp(i eps A_k.sigma / 2) (a reversed traversal stores the inverse)."""
    links = {}
    for e, a in zip(P, Atrav):
        X = expm(1j * eps * sum(a[k] * _SIG[k] for k in range(3)) / 2.0)
        if e.source < e.target:
            links[(e.source, e.target)] = X
        else:
            links[(e.target, e.source)] = X.conj().T
    return links


def test_su2_cubic_term_nonzero_and_eps_independent():
    """Paper 30 Sec. 3.4: 'The O(A^3) term does not vanish' (cubic vertex
    proportional to epsilon^{abc}).  One 4-plaquette (P_A of the n_max = 3
    graph), random field A ~ N(0,1)^{4x3} with seed 1: the odd part
    [s(eps) - s(-eps)] / (2 eps^3) of s = 1 - (1/2) Re tr U_P, through the
    production wilson_action, equals the BCH coefficient c_3 = -0.14264 at
    eps = 1e-2 and 4e-3 (independent of eps to O(eps^2)).

    Rejects: the retired 'the O(A^3) term vanishes' (c_3 = 0); a reversed path
    ordering (flips the sign of c_3)."""
    P = _P_A()
    assert len(P) == 4
    Atrav = np.random.default_rng(1).standard_normal((4, 3))
    c3 = _c3_bch(Atrav)
    assert abs(c3) > 0.1
    for eps in (1e-2, 4e-3):
        sp_ = wilson_action([P], _links_for_traversal_fields(P, Atrav, eps), 1.0)
        sm_ = wilson_action([P], _links_for_traversal_fields(P, Atrav, -eps), 1.0)
        odd = (sp_ - sm_) / (2.0 * eps ** 3)
        assert abs(odd - c3) < 5e-3 * abs(c3) + 1e-4, (eps, odd, c3)
        # the even (quadratic) part is the curl form: (1/8)|sum A_k|^2
        even = (sp_ + sm_) / (2.0 * eps ** 2)
        assert abs(even - np.sum(Atrav.sum(axis=0) ** 2) / 8.0) < 1e-3


def test_u1_action_has_no_cubic_term():
    """Paper 30 Sec. 4.1: 'the abelian action is even in theta and has no cubic
    term'.  Same plaquette and seed, maximal-torus links: the odd part of the
    SU(2) action on diagonal links and of the U(1) action both vanish to
    rounding, relative to eps^3.

    Rejects: a U(1) (torus) action carrying a cubic term."""
    P = _P_A()
    a = np.random.default_rng(1).standard_normal(4)
    for eps in (1e-2, 4e-3):
        out = {}
        for sgn in (1, -1):
            phases, links = {}, {}
            for e, ak in zip(P, a):
                phi = sgn * eps * ak
                key = (min(e.source, e.target), max(e.source, e.target))
                phases[key] = phi if e.source < e.target else -phi
                links[key] = diagonal_su2_from_phase(phases[key])
            out[sgn] = (wilson_action([P], links, 1.0),
                        u1_action_from_su2([P], phases, 1.0))
        for k in (0, 1):
            odd = (out[1][k] - out[-1][k]) / (2.0 * eps ** 3)
            assert abs(odd) < 1e-5, (eps, k, odd)
        even = (out[1][1] + out[-1][1]) / (2.0 * eps ** 2)
        assert abs(even - a.sum() ** 2 / 2.0) < 1e-3     # (beta/2) theta_P^2


# ---------------------------------------------------------------------------
# Guard 5: Rule-B three-plaquette 2-cycles (Sec. 2.3)
# ---------------------------------------------------------------------------

def _rule_b_labels(n_max):
    """(n, l, 2j, 2m_j) built here, independent of geovac's label iterator."""
    out = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for two_j in ([2 * l + 1] + ([2 * l - 1] if l >= 1 else [])):
                for two_mj in range(-two_j, two_j + 1, 2):
                    out.append((n, l, two_j, two_mj))
    return out


def _rule_b_adjacency_independent(n_max):
    labs = _rule_b_labels(n_max)
    V = len(labs)
    A = np.zeros((V, V), dtype=int)
    for i, j in combinations(range(V), 2):
        a, b = labs[i], labs[j]
        if (abs(a[0] - b[0]) <= 1 and abs(a[1] - b[1]) == 1
                and abs(a[2] - b[2]) <= 2 and abs(a[3] - b[3]) <= 2):
            A[i, j] = A[j, i] = 1
    return labs, A


# Printed counts that this guard confirms: 60 at n_max = 2 and 3342 at
# n_max = 3 (both routes agree; the n_max = 3 count is in the slow test below).
# Rejects at n_max = 3: 80, the retired figure, which was the driver's
# cycle_max_count cap on a 200-plaquette truncated complex
# (tests/wilson_rule_b_support/xcwg_monopole_density.py).
RULE_B_TWO_CYCLES = {2: 60, 3: 3342}
RULE_B_N4 = {2: 44, 3: 994}


def _rule_b_two_cycle_check(n_max):
    """Paper 30 Sec. 2.3: 'The smallest closed 2-cycle consists of three
    plaquettes, each pair sharing two edges, with the six edges each used
    twice; there are 60 such 2-cycles at n_max = 2 [...]'.

    Route A: closed triples found from the plaquettes' signed edge vectors
    (a +-1 combination summing to zero).  Route B: K_{2,3} subgraphs counted
    as sum over vertex pairs of C(#common neighbours, 3).  The two sets of
    edge sets must coincide.  The graph is the production Rule-B graph,
    checked against an E1 adjacency built here from the labels.

    Rejects: a cube-like minimal 2-cycle (six plaquettes, 12 edges, 8
    vertices); 'mutually edge-sharing plaquette triples' in place of closed
    2-cycles; any 2-cycle on fewer than three plaquettes; a mis-printed count.
    Returns the number of closed triples."""
    A, labels, _, _ = build_dirac_s3_graph(n_max, "B")
    labs, Ai = _rule_b_adjacency_independent(n_max)
    key = lambda L: (L.n_fock, -L.kappa - 1 if L.kappa < 0 else L.kappa,
                     2 * abs(L.kappa) - 1, L.two_m_j)
    perm = [labs.index(key(L)) for L in labels]
    assert np.array_equal(A, Ai[np.ix_(perm, perm)])
    V = A.shape[0]
    edges = sorted({(min(i, j), max(i, j)) for i, j in zip(*np.nonzero(A))})
    eidx = {e: k for k, e in enumerate(edges)}
    # 4-cycles from common neighbours (bipartite graph, no triangles)
    nbr = [set(np.nonzero(A[v])[0].tolist()) for v in range(V)]
    plaq = set()
    common = {}
    for a, c in combinations(range(V), 2):
        cn = sorted(nbr[a] & nbr[c])
        common[(a, c)] = cn
        for b, d in combinations(cn, 2):
            plaq.add(frozenset([tuple(sorted(x)) for x in
                                ((a, b), (b, c), (c, d), (d, a))]))
    plaq = sorted(plaq, key=sorted)
    assert len(plaq) == RULE_B_N4[n_max]
    # signed edge vectors: orient each 4-cycle a->b->c->d->a
    vecs = []
    for P in plaq:
        es = list(P)
        verts = sorted({v for e in es for v in e})
        a0 = verts[0]
        nb = sorted(v for e in es for v in e if a0 in e and v != a0)
        cyc = [a0, nb[0]]
        while len(cyc) < 4:
            nxt = [v for e in es for v in e if cyc[-1] in e and v != cyc[-1]
                   and v != cyc[-2]]
            cyc.append(nxt[0])
        vec = np.zeros(len(edges), dtype=np.int64)
        for k in range(4):
            u, w = cyc[k], cyc[(k + 1) % 4]
            vec[eidx[(min(u, w), max(u, w))]] += 1 if u < w else -1
        vecs.append(vec)
    vecs = np.array(vecs)
    # no 2-cycle on one or two plaquettes
    for i, j in combinations(range(len(plaq)), 2):
        assert not (np.array_equal(vecs[i], vecs[j])
                    or np.array_equal(vecs[i], -vecs[j]))
    # Route A: closed triples, searched over pairs sharing >= 1 edge
    share = {}
    for i, j in combinations(range(len(plaq)), 2):
        s = len(plaq[i] & plaq[j])
        if s:
            share[(i, j)] = s
    adj_p = {}
    for (i, j) in share:
        adj_p.setdefault(i, set()).add(j)
        adj_p.setdefault(j, set()).add(i)
    closed = set()
    for (i, j) in share:
        for k in adj_p[i] & adj_p[j]:
            if k <= j:
                continue
            for s2 in (1, -1):
                for s3 in (1, -1):
                    if not np.any(vecs[i] + s2 * vecs[j] + s3 * vecs[k]):
                        closed.add((i, j, k))
    # geometry of every closed triple
    for (i, j, k) in closed:
        for x, y in ((i, j), (i, k), (j, k)):
            assert len(plaq[x] & plaq[y]) == 2
        union = plaq[i] | plaq[j] | plaq[k]
        assert len(union) == 6
        uses = [sum(e in plaq[x] for x in (i, j, k)) for e in union]
        assert uses == [2] * 6
        assert len({v for e in union for v in e}) == 5
    # Route B: K_{2,3} count; same edge sets
    k23 = set()
    for (a, c), cn in common.items():
        for t in combinations(cn, 3):
            k23.add(frozenset(tuple(sorted((x, y))) for x in (a, c) for y in t))
    assert sum(comb(len(cn), 3) for cn in common.values()) == len(k23)
    closed_sets = {frozenset(plaq[i] | plaq[j] | plaq[k]) for (i, j, k) in closed}
    assert closed_sets == k23
    assert len(closed) == len(k23)
    return len(closed)


def test_rule_b_three_plaquette_two_cycles_nmax2():
    """Paper 30 Sec. 2.3 and diagnostic (5): geometry and the count 60 at
    n_max = 2 (see _rule_b_two_cycle_check for routes and 'Rejects:')."""
    assert _rule_b_two_cycle_check(2) == RULE_B_TWO_CYCLES[2]


@pytest.mark.slow
def test_rule_b_three_plaquette_two_cycles_geometry_nmax3():
    """Paper 30 Sec. 2.3 at n_max = 3: every closed three-plaquette 2-cycle is
    a K_{2,3} (pairs share two edges, six edges used twice, five vertices),
    the two routes agree, and the count is the printed 3342.
    Rejects: 80, the retired count (a capped, truncated-complex artifact).
    About 50 s."""
    assert _rule_b_two_cycle_check(3) == RULE_B_TWO_CYCLES[3]
