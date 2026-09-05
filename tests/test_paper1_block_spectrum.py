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
   8 - lambda_max = (C + o(1)) / n_max^2, C = pi^2 (2+2^(1/3))^2
   (1+2^(-2/3))/4 = 42.739654...  (log-log slope -1.984 over
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
    # The naive wrong factorisation has n_max^3 eigenvalues against sum n^2,
    # so a `shape != shape or values differ` assertion short-circuits and never
    # compares a single eigenvalue (found by /qa FULL #4: a decorative negative
    # control).  Compare on a SAME-SHAPE wrong factorisation instead, so the
    # value branch actually runs.
    assert wrong.shape != direct.shape                      # still true, recorded
    same_shape_wrong = []
    for l in range(n_max):
        # radial extent off by one: P_{n_max-l+1} x P_{2l+1} truncated back to
        # the right dimension -- same eigenvalue count, wrong grid
        a = np.array([2 - 2 * np.cos(j * np.pi / (n_max - l + 1))
                      for j in range(n_max - l)])
        b = np.array([2 - 2 * np.cos(k * np.pi / (2 * l + 1))
                      for k in range(2 * l + 1)])
        same_shape_wrong.append((a[:, None] + b[None, :]).ravel())
    same_shape_wrong = np.sort(np.concatenate(same_shape_wrong))
    assert same_shape_wrong.shape == direct.shape, (
        same_shape_wrong.shape, direct.shape)
    assert np.abs(direct - same_shape_wrong).max() > 1e-3


SATURATION_C = np.pi ** 2 * (2 + 2 ** (1 / 3)) ** 2 * (1 + 2 ** (-2 / 3)) / 4


def test_saturation_rate_is_order_n_squared():
    """8 - lambda_max = (C + o(1)) / n_max^2 with C = 42.739654... in closed
    form -- a proven rate, not a fit.

    The previous version asserted `42.0 < gaps[-1]*320**2 < 43.5`, a window so
    wide it could not distinguish the n_max = 320 SAMPLE (42.606) from the
    limit (42.7397).  The papers printed the sample as the limit for that
    reason.  This version pins the closed form and the approach separately.
    """
    ns = [20, 40, 80, 160, 320]
    gaps = [8 - max(block_spectrum(n, l).max() for l in range(n)) for n in ns]
    slope = np.polyfit(np.log(ns), np.log(gaps), 1)[0]
    assert -2.02 < slope < -1.95, slope
    # Derive the constant independently instead of comparing the printed
    # expression to the printed number (both authored together -- /qa DELTA #5
    # planted C + 0.05 with the literal updated and this stayed green).
    # Route: minimise the leading deficit form pi^2[(1-x)^-2 + (2x)^-2] over
    # x = l/n_max numerically, with no reference to the closed form.
    from scipy.optimize import minimize_scalar
    g = lambda x: np.pi ** 2 * (1.0 / (1 - x) ** 2 + 1.0 / (2 * x) ** 2)
    opt = minimize_scalar(g, bounds=(0.05, 0.95), method="bounded",
                          options={"xatol": 1e-12})
    assert abs(opt.x - 0.386488) < 1e-5, opt.x
    assert abs(opt.fun - SATURATION_C) < 1e-6, (opt.fun, SATURATION_C)
    # Every finite sample understates C (one-sided bound).  The approach is NOT
    # monotone -- C_20 = 40.7285 > C_21 = 40.6470 (Paper 0 S VI) -- so there is
    # deliberately no monotone assertion here: the doubling grid below happens
    # to be monotone, and a guard asserting that WAS the retired claim
    # (DELTA #7, CODE-B M6; insert 21 into ns and the old assert failed).
    scaled = [g * n ** 2 for g, n in zip(gaps, ns)]
    assert all(v < SATURATION_C for v in scaled), scaled
    assert abs(scaled[-1] - 42.606) < 0.01, scaled[-1]
    # ... and it really does converge there: one cutoff past the printed data
    n = 4000
    far = min(4 * np.sin(np.pi / (2 * (n - l))) ** 2
              + 4 * np.sin(np.pi / (2 * (2 * l + 1))) ** 2
              for l in range(n))
    assert 42.72 < far * n ** 2 < SATURATION_C, far * n ** 2


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
    # Computed blockwise and BASIS-FREE (2026-09-04).  The previous form
    # took an argmax over individual eigenvector components, which is not
    # well defined inside a degenerate eigenspace -- lambda_2s = 3 has
    # multiplicity 25 here, and /qa DELTA #5 measured a no-op relabelling
    # of the nodes turning a reported 2.70% into 4421%.  Pooling each
    # eigenspace first makes the overlap invariant, and diagonalising the
    # 30 l-blocks instead of the 9455 x 9455 whole takes 1 s rather than 70.
    from geovac import lattice_spectrum as _ls
    from scipy.sparse import diags as _diags

    lat = GeometricLattice(max_n=30)
    A = lat.adjacency.tocsr()
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = _diags(deg, 0, shape=A.shape, format="csr") - A
    idx = {s: i for i, s in enumerate(lat.states)}
    _, blocks = _ls.block_eigh(L, lat.states)
    ov = _ls.eigenspace_overlap(blocks, idx[(2, 0, 0)])
    (w0, o0), (w1, o1) = ov[0], ov[1]
    margin = (o0 - o1) / o0
    assert margin < 0.02, margin                       # near-tie
    # The runner-up overlap is itself an exact tie across several eigenspaces
    # (0.25502003 at lambda = 0.0110, 2.8135, 3.1756 for n_max = 30), so "the"
    # runner-up's eigenvalue is decided by sort stability -- 14 of 25 random
    # relabellings flipped it (DELTA #7, CODE-B M1).  Assert on the SET of
    # runners-up inside the tie, which is relabelling-invariant.
    runners = [w for w, o in ov[1:] if abs(o - o1) <= 1e-9 * max(o1, 1.0)]
    assert any(abs(w0 - w) > 1.0 for w in runners), (w0, runners)


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


# --- the s/p proxy is a residue artifact, and closed form on one branch -----
# Paper 1 Eq. (lambda2s_residue), added by /qa trunk FULL run #4 (2026-09-03).

def _sp_proxy_dense(n_max: int) -> tuple:
    """(relative s/p splitting, lambda_2s) by full dense eigendecomposition of
    the production lattice.  Reference route; O(N^3) in Sum n^2, so it is used
    only for the equivalence cross-check."""
    from geovac.lattice import GeometricLattice
    from scipy.sparse import diags
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency
    deg = np.array(A.sum(axis=1)).flatten()
    L = (diags(deg, 0, shape=A.shape, format="csr") - A).toarray()
    idx = {st: i for i, st in enumerate(lat.states)}
    w, v = np.linalg.eigh(L)
    j2s = int(np.argmax(np.abs(v[idx[(2, 0, 0)], :])))
    j2p = int(np.argmax(np.abs(v[idx[(2, 1, 0)], :])))
    return abs(w[j2p] - w[j2s]) / abs(w[j2s]), w[j2s]


def _sp_proxy(n_max: int) -> tuple:
    """Same quantity, computed per block.

    CORRECTED 2026-09-04 (/qa DELTA #5).  The justification printed here was
    wrong: disconnectedness gives block-DIAGONALITY, not block-supported
    eigenvectors.  lambda_2s = 3 is degenerate ACROSS blocks (multiplicity 25
    at n_max = 30, 15 at 18, 14 at 12), so `eigh` returns an arbitrary basis of
    that eigenspace and a full-spectrum argmax over eigenVECTORS is
    basis-dependent -- permuting the node labelling, a mathematical no-op,
    turns 2.70% into 4421%.

    The well-defined quantity is the basis-free one: argmax over eigenvalue
    PROJECTIONS, lambda -> ||P_lambda e_node||.  That is what the block route
    computes (each block is an invariant subspace, so the projection restricted
    to it is exactly the block's own), and it is what the published numbers
    are.  `test_block_route_matches_projection` checks against that rather than
    against a particular LAPACK ordering.
    """
    from geovac.lattice import GeometricLattice
    from scipy.sparse import diags
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency
    deg = np.array(A.sum(axis=1)).flatten()
    idx = {st: i for i, st in enumerate(lat.states)}

    def block_lambda(node, l):
        rows = [i for st, i in idx.items() if st[1] == l]
        pos = {i: k for k, i in enumerate(rows)}
        sub = A[rows, :][:, rows].toarray()
        Lb = np.diag(deg[rows]) - sub
        w, v = np.linalg.eigh(Lb)
        j = int(np.argmax(np.abs(v[pos[idx[node]], :])))
        return w[j]

    lam_2s = block_lambda((2, 0, 0), 0)
    lam_2p = block_lambda((2, 1, 0), 1)
    return abs(lam_2p - lam_2s) / abs(lam_2s), lam_2s


def _sp_proxy_projection(n_max: int) -> tuple:
    """Basis-free reference: argmax over eigenvalue PROJECTIONS.

    For each distinct eigenvalue lambda, ||P_lambda e_node|| is independent of
    which eigenbasis LAPACK returns, so this is well defined even where the
    eigenvalue is degenerate across blocks -- which lambda_2s = 3 is.
    """
    from geovac.lattice import GeometricLattice
    from scipy.sparse import diags
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency
    deg = np.array(A.sum(axis=1)).flatten()
    L = (diags(deg, 0, shape=A.shape, format="csr") - A).toarray()
    idx = {st: i for i, st in enumerate(lat.states)}
    w, v = np.linalg.eigh(L)

    def best(node):
        e = idx[node]
        # group columns by eigenvalue, then take the projection norm
        order = np.argsort(w)
        lam_groups, cur, out = [], [order[0]], []
        for k in order[1:]:
            if abs(w[k] - w[cur[-1]]) < 1e-9:
                cur.append(k)
            else:
                lam_groups.append(cur)
                cur = [k]
        lam_groups.append(cur)
        for grp in lam_groups:
            amp = float(np.linalg.norm(v[e, grp]))
            out.append((amp, float(w[grp[0]])))
        return max(out)[1]

    lam_2s = best((2, 0, 0))
    lam_2p = best((2, 1, 0))
    return abs(lam_2p - lam_2s) / abs(lam_2s), lam_2s


@pytest.mark.parametrize("n_max", [12, 15])
def test_block_route_matches_projection(n_max):
    """The per-block route equals the BASIS-FREE full-spectrum quantity.

    Replaces a check against `eigh`'s raw eigenvector argmax, which passed only
    because LAPACK happens to return block-supported vectors in the production
    node ordering (/qa DELTA #5: permuting the labelling broke it).
    """
    fast, lam_f = _sp_proxy(n_max)
    ref, lam_r = _sp_proxy_projection(n_max)
    assert abs(fast - ref) < 1e-9, (n_max, fast, ref)
    assert abs(lam_f - lam_r) < 1e-9, (n_max, lam_f, lam_r)


@pytest.mark.parametrize("n_max", [12, 15, 18, 21, 24, 27, 30])
def test_lambda_2s_is_exactly_three_on_multiples_of_three(n_max):
    """lambda_2s = 3 EXACTLY iff n_max = 0 (mod 3).  Measured against a dense
    eigendecomposition, not against the closed form."""
    _, lam = _sp_proxy(n_max)
    assert abs(lam - 3.0) < 1e-9, (n_max, lam)


@pytest.mark.parametrize("n_max", [22, 23, 25, 26, 28, 29])
def test_lambda_2s_is_not_three_off_that_branch(n_max):
    """The contrast that makes the residue structure a finding rather than a
    coincidence: off the branch, lambda_2s is measurably away from 3."""
    _, lam = _sp_proxy(n_max)
    assert abs(lam - 3.0) > 0.05, (n_max, lam)


@pytest.mark.parametrize("n_max", [12, 15, 18, 21, 24, 27, 30])
def test_sp_proxy_closed_form_on_the_branch(n_max):
    """On n_max = 0 (mod 3) the proxy is exactly (2 - 2cos(pi/(n_max-1)))/3."""
    ratio, _ = _sp_proxy(n_max)
    closed = (2 - 2 * np.cos(np.pi / (n_max - 1))) / 3
    assert abs(ratio - closed) < 1e-9, (n_max, ratio, closed)


def test_sp_endpoint_is_selection_biased():
    """The reported endpoint sits on the favourable branch: 0.39% at n_max=30,
    against 1.65% at 28 and 2.58% at 29.  This is what retires the series as
    evidence of convergence -- it tracks cutoff divisibility."""
    at30 = _sp_proxy(30)[0] * 100
    at28 = _sp_proxy(28)[0] * 100
    at29 = _sp_proxy(29)[0] * 100
    assert abs(at30 - 0.391) < 0.01, at30
    assert at28 > 3 * at30, (at28, at30)
    assert at29 > 5 * at30, (at29, at30)
