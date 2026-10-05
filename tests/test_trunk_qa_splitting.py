"""
TRUNK QA -- Claim 3: the s/p node-amplitude proxy on the BINARY production
lattice (Paper 1 Sec. III + App. B), computed by the BASIS-FREE rule.

Paper 1 (papers/.../paper_1_spectrum.tex, Sec. III "Convergence Analysis" and
App. B, corrected 2026-09-28): for each of the nodes |2,0,0> (2s) and |2,1,0>
(2p), select the eigenvalue lambda of L = D - A whose eigenspace projection
||P_lambda e_node|| is maximal, and report |E_2p - E_2s| / |E_2s|.  This rule
is basis-free; the rule this file used before 2026-09-28 (argmax over the
components of individual eigenvectors) is NOT -- several eigenvalues are
degenerate across the disconnected l-blocks, so the argmax depended on the
arbitrary LAPACK basis inside those eigenspaces (the /qa trunk FULL #9 finding;
measured here: a random rotation inside the degenerate eigenspaces moves the
old rule's n_max = 7 value from 0.65 % to ~1500 %).

Printed values (binary lattice): 37 / 12.7 / 0.65 / 15.6 / 5.1 / 1.7 / 2.7 /
1.7 / 1.1 / 4.1 / 1.8 / 0.39 % at n_max = 5, 6, 7, 8, 9, 10, 12, 15, 18, 20,
25, 30.  The n_max = 30 endpoint is the on-branch closed form
(2 - 2 cos(pi/29))/3 = 0.3908 %.

WRONG ANSWERS REJECTED:
  * a result that changes under a rotation of the eigenbasis inside a
    degenerate eigenspace (the basis-dependent rule);
  * the pre-2026-09-03 waypoints (13 % at 5, 16 % at 10);
  * any off-branch value other than the printed ones (these were "pinned by
    no test" before this rewrite).

Independence of the invariance check (/qa trunk DELTA #14 guard pass): the
rotation test must NOT build its degenerate eigenspaces with the DEG_TOL it is
guarding -- a test that rotates only inside DEG_TOL-clusters cannot see a
DEG_TOL that splits a genuine (numerically ~1e-15-split) eigenspace, and
indeed DEG_TOL = 1e-16 left the previous version all green while the proxy
moved by up to 44 under a genuine rotation.  The rotation clusters are
therefore built with an independent LOOSE_TOL, and a separate test asserts
that DEG_TOL sits well inside the empty gap window of the spectrum (every
consecutive eigenvalue gap is either numerical noise << DEG_TOL or a genuine
gap >> DEG_TOL), so the clustering is not tolerance-sensitive.

Fire-test plants (debug/qa/fire_test.py):
  * _proxy switched to the raw-eigenvector argmax (cl = one cluster per
    eigenvector) -> the rotation-invariance test FIRES;
  * DEG_TOL = 1e-8 -> 1e-16 -> the rotation-invariance test (loose clusters)
    AND the gap-window test FIRE.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import diags

from geovac import GeometricLattice

# Paper 1 printed values, in percent
PRINTED = {5: "37", 6: "12.7", 7: "0.65", 8: "15.6", 9: "5.1", 10: "1.7", 12: "2.7",
           15: "1.7", 18: "1.1", 20: "4.1", 25: "1.8", 30: "0.39"}


def _matches_printed(got_percent: float, printed: str) -> bool:
    """Within half a unit of the printed last digit."""
    decimals = len(printed.split(".")[1]) if "." in printed else 0
    return abs(got_percent - float(printed)) <= 0.5 * 10 ** (-decimals) + 1e-12
DEG_TOL = 1e-8  # eigenvalue clustering; analytically-equal eigenvalues agree to ~1e-13
# Independent, looser tolerance used ONLY by the guards below to define the
# degenerate eigenspaces that are rotated.  Measured (n_max <= 18): numerically
# split degenerate eigenvalues differ by <= ~6e-15; the smallest genuine gap is
# 1.56e-6 (n_max = 10).  1e-7 lies in that empty window, independently of DEG_TOL.
LOOSE_TOL = 1e-7


def _laplacian(max_n: int):
    lat = GeometricLattice(max_n=max_n)
    A = lat.adjacency
    deg = np.array(A.sum(axis=1)).flatten()
    L = (diags(deg, 0, shape=A.shape, format="csr") - A).toarray()
    idx = {s: i for i, s in enumerate(lat.states)}
    return L, idx[(2, 0, 0)], idx[(2, 1, 0)]


def _clusters(w: np.ndarray, tol: float | None = None):
    tol = DEG_TOL if tol is None else tol
    out, start = [], 0
    for i in range(1, len(w) + 1):
        if i == len(w) or w[i] - w[i - 1] > tol:
            out.append((start, i))
            start = i
    return out


def _proxy(w: np.ndarray, v: np.ndarray, i2s: int, i2p: int) -> float:
    """Basis-free: for each node pick the eigenvalue with max ||P_lambda e_node||^2."""
    cl = _clusters(w)

    def pick(node: int) -> float:
        amps = [float(np.sum(np.abs(v[node, a:b]) ** 2)) for a, b in cl]
        a, b = cl[int(np.argmax(amps))]
        return float(np.mean(w[a:b]))

    lam2s, lam2p = pick(i2s), pick(i2p)
    return abs(lam2p - lam2s) / abs(lam2s)


def _splitting(max_n: int) -> float:
    L, i2s, i2p = _laplacian(max_n)
    w, v = np.linalg.eigh(L)
    return _proxy(w, v, i2s, i2p)


def _rotate_within_degenerate(w, v, rng):
    """Random orthogonal rotation inside each degenerate eigenspace, the
    eigenspaces defined by the INDEPENDENT LOOSE_TOL (never DEG_TOL)."""
    vr = v.copy()
    for a, b in _clusters(w, LOOSE_TOL):
        if b - a > 1:
            Q, _ = np.linalg.qr(rng.standard_normal((b - a, b - a)))
            vr[:, a:b] = v[:, a:b] @ Q
    return vr


FAST_GRID = [5, 6, 7, 8, 9, 10, 12, 15, 18]


@pytest.mark.parametrize("n_max", FAST_GRID)
def test_printed_values_reproduce(n_max):
    """Every printed value on the fast grid, to its printed precision."""
    got = 100 * _splitting(n_max)
    assert _matches_printed(got, PRINTED[n_max]), (
        f"n_max={n_max}: {got:.3f}% vs printed {PRINTED[n_max]}%")


@pytest.mark.parametrize("n_max", FAST_GRID)
def test_result_is_invariant_under_rotation_in_degenerate_eigenspaces(n_max):
    """The basis-free rule gives the same number for any orthonormal basis of
    each degenerate eigenspace, at every cutoff of the fast grid (Paper 1
    App. B: n_max = 5..18; the old argmax rule and a too-tight DEG_TOL were
    shown to move at 6, 7, 9, 12, 15, 18).  The eigenspaces rotated are
    defined by LOOSE_TOL, independent of the DEG_TOL under test."""
    L, i2s, i2p = _laplacian(n_max)
    w, v = np.linalg.eigh(L)
    base = _proxy(w, v, i2s, i2p)
    rng = np.random.default_rng(n_max)
    for _ in range(3):
        vr = _rotate_within_degenerate(w, v, rng)
        assert np.allclose(vr.T @ L @ vr, np.diag(w), atol=1e-8)  # still an eigenbasis
        assert abs(_proxy(w, vr, i2s, i2p) - base) < 1e-12, (
            f"n_max={n_max}: result moved under an eigenbasis rotation")


@pytest.mark.parametrize("n_max", FAST_GRID)
def test_deg_tol_sits_inside_the_empty_gap_window(n_max):
    """Every consecutive eigenvalue gap is either numerical noise (a split
    degenerate eigenvalue) or a genuine gap, and DEG_TOL separates the two
    with margin: noise << DEG_TOL << genuine gap.  So the DEG_TOL-clusters are
    exactly the LOOSE_TOL-clusters and the clustering is not tolerance-tuned."""
    L, _, _ = _laplacian(n_max)
    w = np.linalg.eigvalsh(L)
    g = np.diff(w)
    noise = g[g < LOOSE_TOL]
    genuine = g[g >= LOOSE_TOL]
    assert noise.max() < 1e-13, f"n_max={n_max}: split-degeneracy noise {noise.max():.2e}"
    assert genuine.min() > 1e-6, f"n_max={n_max}: genuine gap {genuine.min():.2e}"
    assert 1e3 * noise.max() < DEG_TOL, (
        f"n_max={n_max}: DEG_TOL={DEG_TOL:.0e} does not clear the noise {noise.max():.2e}")
    assert genuine.min() > 1e2 * DEG_TOL, (
        f"n_max={n_max}: inter-cluster gap {genuine.min():.2e} not >> DEG_TOL={DEG_TOL:.0e}")
    assert _clusters(w) == _clusters(w, LOOSE_TOL)


def test_nonmonotone_and_downward():
    """Strongly non-monotone (15.6 % at 8 after 0.65 % at 7) and trending down."""
    vals = [_splitting(n) for n in FAST_GRID]
    assert any(vals[i] > vals[i - 1] for i in range(1, len(vals)))
    assert _splitting(18) < _splitting(8) / 3


def test_pre_2026_09_03_waypoints_do_not_reproduce():
    """13 % at 5 and 16 % at 10 (the list printed before 2026-09-03) are not the measurement."""
    assert abs(_splitting(5) - 0.13) > 0.05
    assert abs(_splitting(10) - 0.16) > 0.05


@pytest.mark.slow
@pytest.mark.xdist_group("heavy_memory")
def test_printed_values_n_max_20_25():
    """The off-branch excursion 4.1 % at 20 and 1.8 % at 25 (dense eigh, ~20 s)."""
    for n in (20, 25):
        got = 100 * _splitting(n)
        assert _matches_printed(got, PRINTED[n]), f"n_max={n}: {got:.3f}%"


@pytest.mark.slow
@pytest.mark.parametrize("n_max", [20, pytest.param(25, marks=pytest.mark.xdist_group("heavy_memory")), pytest.param(30, marks=pytest.mark.xdist_group("heavy_memory"))])
def test_deg_tol_gap_window_slow_values(n_max):
    """The slow printed values (4.1 % at 20, 1.8 % at 25, 0.39 % at 30) rest on
    DEG_TOL too, so the noise << DEG_TOL < genuine-gap check is repeated here.

    Measured 2026-09-28 (/qa DELTA #15 guard pass), and why this is NOT a copy
    of the fast-grid test: the fast test's LOOSE_TOL = 1e-7 split and its
    100 x DEG_TOL margin do not hold at these cutoffs.  At n_max = 25 two
    GENUINELY distinct eigenvalues (0.3971727246 and 0.3971727563, from
    different l-blocks, 13 and 24) sit 3.17e-8 apart -- only ~3 x DEG_TOL, and
    inside LOOSE_TOL, so LOOSE_TOL-clusters would wrongly merge them; at
    n_max = 30 the smallest genuine gap is 4.3e-7 (~43 x DEG_TOL).  So here
    noise and genuine gaps are separated by NOISE_CUT = 1e-10 (numerical noise
    <= ~1e-13, genuine >= 3e-8 at every cutoff measured), and the load-bearing
    check is the one the printed value actually needs: it is IDENTICAL for
    every clustering tolerance from 1e-12 to 1e-6 (measured: 1.7757 % and
    0.39080 % unchanged across that range).
    Rejects: a DEG_TOL below the noise (1e-16 splits numerically-degenerate
    eigenspaces) and a DEG_TOL above a genuine gap."""
    global DEG_TOL
    NOISE_CUT = 1e-10
    L, i2s, i2p = _laplacian(n_max)
    w, v = np.linalg.eigh(L)
    del L
    g = np.diff(w)
    noise = g[g < NOISE_CUT]
    genuine = g[g >= NOISE_CUT]
    assert noise.max() < 1e-12, f"n_max={n_max}: split-degeneracy noise {noise.max():.2e}"
    assert genuine.min() > 1e-8, f"n_max={n_max}: genuine gap {genuine.min():.2e}"
    assert 1e3 * noise.max() < DEG_TOL, (
        f"n_max={n_max}: DEG_TOL={DEG_TOL:.0e} does not clear the noise {noise.max():.2e}")
    assert genuine.min() > DEG_TOL, (
        f"n_max={n_max}: genuine gap {genuine.min():.2e} not above DEG_TOL={DEG_TOL:.0e}")
    assert _clusters(w) == _clusters(w, NOISE_CUT)
    ref = _proxy(w, v, i2s, i2p)
    assert _matches_printed(100 * ref, PRINTED[n_max]), f"n_max={n_max}: {100*ref:.3f}%"
    saved = DEG_TOL
    try:
        for tol in (1e-12, 1e-10, 1e-9, 1e-7, 1e-6):
            DEG_TOL = tol
            got = _proxy(w, v, i2s, i2p)
            assert abs(got - ref) < 1e-12, (
                f"n_max={n_max}: splitting moved with the clustering tolerance "
                f"({tol:.0e}: {100*got:.5f}% vs {100*ref:.5f}%)")
    finally:
        DEG_TOL = saved


@pytest.mark.slow
@pytest.mark.xdist_group("heavy_memory")
def test_endpoint_n_max_30_is_the_on_branch_closed_form():
    """n_max = 30 (9455 nodes, minutes): 0.39 %, within 1e-3 absolute of the
    printed value and equal to (2 - 2 cos(pi/29))/3."""
    s30 = _splitting(30)
    assert abs(s30 - float(PRINTED[30]) / 100) < 1e-3, f"n_max=30: {100*s30:.4f}%"
    closed = (2 - 2 * np.cos(np.pi / 29)) / 3
    assert abs(s30 - closed) < 1e-9
