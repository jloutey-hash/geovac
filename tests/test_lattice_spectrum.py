"""Guards for geovac/lattice_spectrum.py.

Written as its own pass (CLAUDE.md §9: guard-writing is a separate,
separately-reviewed activity), after DELTA #7 recorded the module as having no
dedicated test file and its basis-freeness as unverified.

Every guard below names, in its docstring, THE WRONG ANSWER IT REJECTS, and
each was fire-tested against that answer by planting it in
``geovac/lattice_spectrum.py``.  A guard whose rejected answer cannot be named
is not a guard.

FIRE-TESTED 2026-09-04, 14 plants, 14 FIRED (two only after being fixed):

    closed form: a path factor mis-sized                     FIRED
    closed form: factors multiplied instead of added         FIRED
    spectrum: a block dropped                                FIRED
    lambda_max: max -> min                                   FIRED
    lambda_max: scan restricted to l = 0                     FIRED*
    lambda_max: O(n^3) revert                                FIRED
    path eigenvectors: half-offset dropped                   FIRED
    block_eigh: partition by n instead of l                  FIRED
    eigenspace_overlap: POOLING DISABLED (while False)       FIRED*
    eigenspace_overlap: normalisation slip                   FIRED
    eigenspace_overlap: blockwise != global projector        FIRED
    block_dims: range check removed                          FIRED
    block_dims: factors swapped                              FIRED
    path_spectrum: off-by-one                                FIRED

* did NOT fire on the first attempt, and both reasons were structural rather
  than weak assertions:

  - the lambda_max guard asserted about ``block_spectrum`` and never called
    ``lambda_max``, so a plant on the scan could not reach it;

  - the pooling guard ran on row (2,0,0), which lives in the l = 0 block --
    P_{n_max} x P_1, a bare path, whose spectrum is SIMPLE.
    ``eigenspace_overlap`` reads only the block containing its row, so the
    pooling branch never executed and disabling it changed nothing.  This is
    also why the same plant did not fire against the converted test in
    ``tests/test_paper1_block_spectrum.py``: that test's basis-freeness is
    trivially true for its row.  The guard here now runs on a row in the
    l = 1 block (P_{n_max-1} x P_3, five degenerate pairs), where pooling
    actually happens.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.sparse import diags

from geovac import lattice_spectrum as ls
from geovac.lattice import GeometricLattice

# Cutoffs deliberately different from the ones used while developing the
# module (8, 20, 30) and from CODE-B's re-derivation (6, 11, 17).
CUTOFFS = [3, 5, 9, 14]


def _dense(n_max: int):
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency.tocsr()
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = diags(deg, 0, shape=A.shape, format="csr") - A
    return lat, L


# ---------------------------------------------------------------- spectrum


@pytest.mark.parametrize("n_max", CUTOFFS)
def test_closed_form_reproduces_the_built_laplacian(n_max: int) -> None:
    """REJECTS: a closed form that is not this graph's spectrum.

    Wrong answers this fails on -- each planted and confirmed:
      * either path factor mis-sized (P_{n-l} x P_{2l+1} -> P_{n-l+1} or P_{2l}),
      * the two factors multiplied instead of added,
      * any per-block eigenvalue perturbed.
    Compares the FULL multiset against a dense eigendecomposition of the
    lattice that was actually constructed, not against a second closed form.
    """
    _, L = _dense(n_max)
    ref = np.sort(np.linalg.eigvalsh(L.toarray()))
    got = ls.spectrum(n_max)
    assert got.shape == ref.shape, (got.shape, ref.shape)
    assert np.max(np.abs(ref - got)) < 1e-10


@pytest.mark.parametrize("n_max", CUTOFFS)
def test_spectrum_size_is_the_node_count(n_max: int) -> None:
    """REJECTS: a spectrum with the wrong number of eigenvalues.

    n(n+1)(2n+1)/6 is the node count;  a dropped or double-counted block
    changes it while leaving every individual eigenvalue correct.
    """
    assert ls.spectrum(n_max).size == n_max * (n_max + 1) * (2 * n_max + 1) // 6


@pytest.mark.parametrize("n_max", CUTOFFS)
def test_lambda_max_agrees_with_the_operator_and_the_closed_form(n_max: int) -> None:
    """REJECTS: a lambda_max that is not the true top eigenvalue.

    Three routes must agree: the O(n) shortcut, the full closed-form spectrum,
    and a dense eigendecomposition of the built operator.  Rejects
    ``max`` -> ``min``, a dropped block (the argmax lives in a mid-l block,
    not l = 0), and the factor-top shortcut being taken on the wrong index.
    """
    lat, L = _dense(n_max)
    dense_top = float(np.max(np.linalg.eigvalsh(L.toarray())))
    assert ls.lambda_max(n_max) == pytest.approx(ls.spectrum(n_max)[-1], abs=1e-12)
    assert ls.lambda_max(n_max) == pytest.approx(dense_top, abs=1e-9)
    assert ls.lambda_max_from_operator(L, lat.states) == pytest.approx(dense_top, abs=1e-9)


def test_lambda_max_argmax_is_a_mid_l_block_not_the_s_wave() -> None:
    """REJECTS: the s-wave block being read as the saturating one.

    Paper 0 SS VI: the bound 2 d_max = 8 is attained in a mid-l block;  the
    l = 0 block alone saturates only at 4.  A lambda_max that scanned only
    l = 0 would return ~4 and still 'converge'.
    """
    n = 30
    tops = [ls.block_spectrum(n, l)[-1] for l in range(n)]
    arg = int(np.argmax(tops))
    assert 0 < arg < n - 1, arg
    # the s-wave block's own ceiling is 4; the graph's is 8
    assert ls.block_spectrum(n, 0)[-1] < 4.0 + 1e-9
    # lambda_max must see the mid-l block, so it must exceed the s-wave
    # ceiling and equal the mid-l top.  A scan restricted to l = 0 returns
    # <= 4 and fails here -- the earlier version of this test asserted only
    # about block_spectrum and never called lambda_max, so that plant
    # left it green (DELTA #7 guard pass).
    assert ls.lambda_max(n) > 4.0 + 1e-9
    assert ls.lambda_max(n) == pytest.approx(tops[arg], abs=1e-12)
    assert max(tops) > 7.9


def test_lambda_max_is_linear_time_not_cubic() -> None:
    """REJECTS: an implementation that materialises every block.

    The first version did, and could not finish at n_max = 5000.  The top of
    a Cartesian product is the sum of the factor tops, so this is O(n).  A
    cubic implementation cannot answer at n_max = 20000 in seconds.
    """
    assert ls.lambda_max(20000) == pytest.approx(8.0, abs=1e-6)
    assert ls.lambda_max(20000) < 8.0          # one-sided: always understates


@pytest.mark.parametrize("n_max", [10, 70, 300])
def test_saturation_is_one_sided_and_never_exceeds_the_bipartite_bound(n_max: int) -> None:
    """REJECTS: any lambda_max at or above 2 d_max = 8.

    The bound is one-sided -- this, not monotonicity, is what makes every
    finite sample an understatement (Paper 0 SS VI;  the approach is NOT
    monotone, C_20 = 40.7285 > C_21 = 40.6470).
    """
    assert ls.lambda_max(n_max) < 8.0


# ------------------------------------------------------------ eigenvectors


@pytest.mark.parametrize("m", [1, 2, 4, 7])
def test_path_eigenvectors_are_orthonormal_eigenvectors(m: int) -> None:
    """REJECTS: eigenvectors that do not diagonalise the path Laplacian.

    Rejects the half-offset dropped (cos(j pi i/m) instead of
    cos(j pi (i+1/2)/m)), which is the Dirichlet family, not the Neumann one
    the graph has -- it satisfies neither the eigen-relation nor the boundary
    condition here.
    """
    V = ls.path_eigenvectors(m)
    lam = ls.path_spectrum(m)
    Lp = np.diag([1.0] + [2.0] * (m - 2) + [1.0]) if m >= 2 else np.zeros((1, 1))
    if m >= 2:
        Lp -= np.diag(np.ones(m - 1), 1) + np.diag(np.ones(m - 1), -1)
    assert np.max(np.abs(V.T @ V - np.eye(m))) < 1e-10
    assert np.max(np.abs(Lp @ V - V * lam)) < 1e-10


# -------------------------------------------------------------- block_eigh


@pytest.mark.parametrize("n_max", [5, 9])
def test_block_eigh_reproduces_the_dense_spectrum(n_max: int) -> None:
    """REJECTS: a block decomposition that is not the operator's.

    Rejects block_index keying on n or m instead of l -- the lattice is block
    diagonal in l ONLY (no edge changes l), so any other partition gives a
    different, wrong spectrum.
    """
    lat, L = _dense(n_max)
    ref = np.sort(np.linalg.eigvalsh(L.toarray()))
    got, blocks = ls.block_eigh(L, lat.states)
    assert np.max(np.abs(ref - got)) < 1e-10
    assert sum(sel.size for sel, _, _ in blocks) == len(lat.states)


@pytest.mark.parametrize("n_max", [6, 11])
def test_block_index_partitions_by_l_and_no_edge_crosses_a_block(n_max: int) -> None:
    """REJECTS: a partition that is not by l, or a graph where edges cross it.

    This is the structural fact the whole module rests on.  If any edge
    changed l, block-diagonalising by l would be wrong and every speedup
    with it.
    """
    lat = GeometricLattice(max_n=n_max)
    A = lat.adjacency.tocsr()
    st = np.asarray(lat.states)
    rows, cols = A.nonzero()
    assert np.all(st[rows, 1] == st[cols, 1]), "an edge changes l"
    idx = ls.block_index(lat.states)
    assert sorted(idx) == sorted(set(int(v) for v in st[:, 1]))
    assert sum(v.size for v in idx.values()) == len(lat.states)


# ------------------------------------------------- eigenspace_overlap: the
# property DELTA #7 recorded as unverified


@pytest.mark.parametrize("n_max", [13, 16])
def test_eigenspace_overlap_is_invariant_under_degenerate_basis_rotation(n_max: int) -> None:
    """REJECTS: pooling disabled -- the guard DELTA #7 could not write.

    THE WRONG ANSWER: ``while False:`` in the pooling loop, i.e. treating each
    eigenvalue as its own eigenspace.  Then the returned numbers are
    individual eigenvector amplitudes, which are arbitrary inside a degenerate
    eigenspace: rotating an orthonormal basis of that eigenspace -- a
    mathematical no-op -- changes them.  Pooling first makes ||P_lambda e||
    invariant, and this test measures exactly that difference.

    THE ROW MATTERS.  ``eigenspace_overlap`` reads only the block containing
    the row, and (2,0,0) is in the l = 0 block -- P_{n_max} x P_1, a bare
    path, whose spectrum is SIMPLE.  On that row the pooling branch never
    executes, so disabling it changes nothing and the guard cannot fail;  the
    first version of this test ran there and did not fire (and that is the
    same root cause as the plant CODE-B reported not firing against
    test_paper1_block_spectrum.py's converted test).  This one runs on a row
    in the l = 1 block, P_{n_max-1} x P_3, which has five degenerate pairs.
    """
    lat, L = _dense(n_max)
    idx = {s: i for i, s in enumerate(lat.states)}
    row = next(i for st, i in idx.items() if st[1] == 1)
    # precondition: the chosen row's block really is degenerate
    blk = ls.block_spectrum(n_max, 1)
    assert any(abs(a - b) < 1e-9 for a, b in zip(blk, blk[1:])), \
        "chosen block has a simple spectrum; the pooling branch would not run"
    _, blocks = ls.block_eigh(L, lat.states)
    base = ls.eigenspace_overlap(blocks, row)

    rng = np.random.default_rng(20260904)
    rotated, n_rotated = [], 0
    for sel, w, v in blocks:
        v = v.copy()
        order = np.argsort(w)
        w_s = w[order]
        i = 0
        while i < len(w_s):
            j = i
            while j + 1 < len(w_s) and abs(w_s[j + 1] - w_s[i]) <= 1e-9:
                j += 1
            if j > i:                       # a genuinely degenerate eigenspace
                cols = order[i:j + 1]
                q, _ = np.linalg.qr(rng.normal(size=(j - i + 1, j - i + 1)))
                v[:, cols] = v[:, cols] @ q
                n_rotated += 1
            i = j + 1
        rotated.append((sel, w, v))

    after = ls.eigenspace_overlap(rotated, row)
    assert len(after) == len(base)
    assert max(abs(a[1] - b[1]) for a, b in zip(base, after)) < 1e-9
    assert max(abs(a[0] - b[0]) for a, b in zip(base, after)) < 1e-9
    # Non-vacuity, reported by the loop that did the work: at least one
    # eigenspace really was degenerate and really was rotated.  Without this
    # the test would pass on a spectrum with no degeneracies at all, where
    # pooled and unpooled agree trivially -- and would then not reject the
    # wrong answer it names.
    assert n_rotated > 0, "no degenerate eigenspace was rotated; test is vacuous"
    assert sum(1 for _, o in base if o > 1e-12) < sum(
        len(w) for _, w, _ in blocks), "every eigenvalue distinct: nothing pooled"


def test_eigenspace_overlaps_are_a_partition_of_unity() -> None:
    """REJECTS: overlaps that are not projector norms.

    sum ||P_lambda e||^2 = ||e||^2 = 1 over a complete set of eigenspaces.
    Rejects a normalisation slip or a dropped eigenspace, neither of which
    the invariance test above would catch.
    """
    lat, L = _dense(9)
    idx = {s: i for i, s in enumerate(lat.states)}
    _, blocks = ls.block_eigh(L, lat.states)
    ov = ls.eigenspace_overlap(blocks, idx[(2, 0, 0)])
    assert sum(o ** 2 for _, o in ov) == pytest.approx(1.0, abs=1e-10)


def test_eigenspace_overlap_matches_the_global_dense_projector() -> None:
    """REJECTS: a blockwise overlap that is not the global one.

    The block route must agree with projector norms computed from a dense
    eigendecomposition of the WHOLE operator, including across-block
    degeneracies (lambda = 3 is degenerate across l-blocks).
    """
    n = 9
    lat, L = _dense(n)
    idx = {s: i for i, s in enumerate(lat.states)}
    row = idx[(2, 0, 0)]
    w, v = np.linalg.eigh(L.toarray())
    order = np.argsort(w)
    w, v = w[order], v[:, order]
    glob, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) <= 1e-9:
            j += 1
        glob.append((float(w[i]), float(np.linalg.norm(v[row, i:j + 1]))))
        i = j + 1
    glob.sort(key=lambda t: -t[1])
    _, blocks = ls.block_eigh(L, lat.states)
    got = ls.eigenspace_overlap(blocks, row)
    assert max(abs(a[1] - b[1]) for a, b in zip(glob, got)) < 1e-9


def test_eigenspace_overlap_rejects_a_row_in_no_block() -> None:
    """REJECTS: silently returning something for an out-of-range row."""
    lat, L = _dense(5)
    _, blocks = ls.block_eigh(L, lat.states)
    with pytest.raises(ValueError):
        ls.eigenspace_overlap(blocks, 10_000)


# ------------------------------------------------------------- input guards


def test_block_dims_rejects_out_of_range_l() -> None:
    """REJECTS: silent acceptance of l >= n_max or l < 0, which would return
    a negative or zero path length and a meaningless spectrum."""
    with pytest.raises(ValueError):
        ls.block_dims(5, 5)
    with pytest.raises(ValueError):
        ls.block_dims(5, -1)
    with pytest.raises(ValueError):
        ls.block_dims(0, 0)


def test_block_dims_are_the_grid_factors() -> None:
    """REJECTS: the two factors swapped -- which leaves the SUM of the tops
    unchanged for square blocks and so survives a lambda_max check."""
    assert ls.block_dims(30, 11) == (19, 23)
    assert sum(ls.block_dims(30, l)[0] * ls.block_dims(30, l)[1]
               for l in range(30)) == 30 * 31 * 61 // 6


def test_path_spectrum_endpoints() -> None:
    """REJECTS: an off-by-one in the path spectrum.

    P_m has a zero mode (constants) and top 2 - 2cos((m-1)pi/m) < 4.
    """
    for m in (1, 2, 5, 13):
        s = ls.path_spectrum(m)
        assert s.size == m
        assert s[0] == pytest.approx(0.0, abs=1e-14)
        assert s[-1] == pytest.approx(2 - 2 * math.cos((m - 1) * math.pi / m), abs=1e-14)
        assert s[-1] < 4.0
