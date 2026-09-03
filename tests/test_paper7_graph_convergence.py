"""Paper 7, "New contributions" item 1 -- what the discrete graph Laplacian is
MEASURED to do through n_max = 30 (re-scoped 2026-09-03).

HISTORY.  The 2026-09-02 version of this file pinned "the ground-state energy
of the production Hamiltonian against -1/2".  Trunk FULL run #3 showed
(three reviewers independently, PM-recomputed) that this was a false
positive for the paper's convergence claim:  H = kappa (D - A) with
kappa := -0.5/8, so E_0 = -lambda_max(L)/16 IDENTICALLY and "E_0 -> -1/2" is
the same statement as "lambda_max -> 8 = 2 d_max" -- a graph-combinatorial
saturation (bipartite Laplacian bound), rescaled by the constant that was
matched to the target.  Paper 0's abstract had said exactly this ("the
content that could have failed is the spectral bound lambda_max -> 8, not
the target the scale is matched to").  This file now pins the spectral
bound and records the structural facts that limit the reading:

  * the lattice's edges change n (T+-) or m (L+-) but never l, so L splits
    into n_max connected components, one per l, with an n_max-dimensional
    kernel (one constant per block);
  * the l = 0 block is a path graph in n whose top eigenvalue tends to 4
    (3.989 at n_max = 30), i.e. the s-wave block alone gives kappa * 4 = -1/4;
  * lambda_max of the full graph is attained in a mid-l block (l = 11 at
    n_max = 30) whose mode has ~1e-33 weight on the (1,0,0) node;
  * H's spectrum is confined to [-1/2, 0] (L <= 2 d_max), so no level-by-
    level Rydberg ladder is reproduced (six lowest eigenvalues within 0.06%
    of each other at n_max = 30).

Measured sequence (production objects only; kappa = -1/16, Z = 1):

    n_max      :   5        8        10       15       20       30
    lambda_max : 6.618034 7.419972 7.611436 7.821269 7.898179 7.954094
    E_0 = k*lm : -0.41363 -0.46375 -0.47571 -0.48883 -0.49364 -0.49713
    |E_0+1/2|  :  17.27%    7.25%    4.86%    2.23%    1.27%    0.57%

NO rate is asserted (Paper 7 proves none).  The operator convergence
L -> Delta_{S^3} and the S^3 identification of the limit are NOT tested by
anything in tests/ (coverage gap, carryforward I.0.1).  Tier: [MEASURED] for
the bound; the convergence reading is [OBSERVATION].
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

from geovac.atomic_solver import AtomicSolver
from geovac.lattice import GeometricLattice

KAPPA = -1.0 / 16.0
GRID = [5, 8, 10, 15, 20, 30]
MEASURED_LMAX = {5: 6.618034, 8: 7.419972, 10: 7.611436,
                 15: 7.821269, 20: 7.898179, 30: 7.954094}


def _laplacian(max_n: int):
    lat = GeometricLattice(max_n)
    A = lat.adjacency.tocsr()
    deg = np.array(A.sum(axis=1)).ravel()
    return lat, A, diags(deg) - A, int(deg.max())


@pytest.fixture(scope="module")
def lmax() -> dict:
    out = {}
    for n in GRID:
        _, _, L, _ = _laplacian(n)
        out[n] = float(eigsh(L, k=1, which="LA")[0][0])
    return out


def test_lambda_max_reproduces_measured_values(lmax):
    for n in GRID:
        assert abs(lmax[n] - MEASURED_LMAX[n]) < 5e-5, (n, lmax[n])


def test_lambda_max_saturates_bipartite_bound_monotonically(lmax):
    """lambda_max increases toward 2 d_max = 8 at every step and never exceeds it."""
    for n in GRID:
        _, _, _, dmax = _laplacian(n)
        assert dmax == 4
        assert lmax[n] <= 2 * dmax + 1e-9
    for a, b in zip(GRID, GRID[1:]):
        assert lmax[b] > lmax[a]
    assert 0.004 < (8 - lmax[30]) / 8 < 0.0075          # 0.57%, two-sided window


def test_ground_energy_is_kappa_times_lambda_max_by_construction(lmax):
    """E_0 -> -1/2 carries no information beyond lambda_max -> 8."""
    for n in GRID:
        E0 = float(eigsh(AtomicSolver(n, Z=1).H, k=1, which="SA")[0][0])
        assert abs(E0 - KAPPA * lmax[n]) < 1e-9, (n, E0, KAPPA * lmax[n])


def test_laplacian_splits_into_one_block_per_l():
    """Edges change n or m, never l: n_max components, n_max-dim kernel."""
    for n in (6, 12, 30):
        lat, A, L, _ = _laplacian(n)
        st = np.array(lat.states)
        rows, cols = A.nonzero()
        assert np.all(st[rows, 1] == st[cols, 1])           # dl = 0 on every edge
        ncomp, _ = connected_components(A, directed=False)
        assert ncomp == n
        zeros = np.sum(np.abs(np.linalg.eigvalsh(L.toarray())) < 1e-9) if n <= 12 else None
        if zeros is not None:
            assert zeros == n


def test_s_wave_block_saturates_four_not_eight():
    lat, A, L, _ = _laplacian(30)
    st = np.array(lat.states)
    m = st[:, 1] == 0
    L0 = L[m][:, m]
    top0 = float(eigsh(L0, k=1, which="LA")[0][0])
    assert 3.98 < top0 < 4.0
    assert abs(KAPPA * top0 - (-0.25)) < 1e-3


def test_extremal_mode_has_no_1s_weight():
    lat, A, L, _ = _laplacian(30)
    st = np.array(lat.states)
    w, v = eigsh(L, k=1, which="LA")
    vec = v[:, 0] ** 2
    idx_1s = lat.states.index((1, 0, 0))
    assert vec[idx_1s] < 1e-20
    per_l = {l: float(vec[st[:, 1] == l].sum()) for l in range(30)}
    best = max(per_l, key=per_l.get)
    assert best == 11 and per_l[best] > 0.999


def test_spectrum_confined_and_bottom_dense_at_nmax_30():
    """H's spectrum lies in [-1/2, 0]; its six lowest eigenvalues sit within
    0.1% of each other -- not a Rydberg ladder."""
    H = AtomicSolver(30, Z=1).H
    lo = np.sort(eigsh(H, k=6, which="SA")[0])
    hi = float(eigsh(H, k=1, which="LA")[0][0])
    assert lo[0] > -0.5 and hi <= 1e-9
    assert (lo[-1] - lo[0]) / abs(lo[0]) < 1e-3


@pytest.mark.slow
def test_lambda_max_deficit_at_nmax_70_slow():
    """Paper 0 conclusion: 0.11% saturation deficit at n_max = 70 (116,795
    nodes; ~3 min).  Measured 2026-09-03: E0 = -0.499463."""
    _, _, L, dmax = _laplacian(70)
    lm = float(eigsh(L, k=1, which="LA")[0][0])
    assert dmax == 4
    assert 0.0009 < (8 - lm) / 8 < 0.0013


def test_energy_scales_as_z_squared():
    """The He+ row of docs/validation_benchmarks.md.

    DELTA #4 flagged that row as tautological: `AtomicSolver` sets
    `kinetic_scale *= Z**2`, so the RELATIVE saturation deficit is Z-independent
    by construction and no test built the solver at Z != 1 at all.  What is
    genuinely checkable is that the implemented scaling is exactly Z^2 and that
    the graph itself does not depend on Z -- which is what this pins.  It cannot
    detect a spectral error; the H row does that."""
    from geovac.atomic_solver import AtomicSolver
    base = AtomicSolver(max_n=10, Z=1)
    e0_base = float(base.compute_ground_state()[0][0])
    for Z in (2, 3, 7):
        sol = AtomicSolver(max_n=10, Z=Z)
        e0 = float(sol.compute_ground_state()[0][0])
        assert abs(e0 - Z ** 2 * e0_base) < 1e-9
        # and the graph itself carries no Z: same adjacency, same degrees
        assert (sol.lattice.adjacency
                != base.lattice.adjacency).nnz == 0
