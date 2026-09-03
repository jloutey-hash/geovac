"""
TRUNK QA -- Claim 3: 2s/2p splitting on the BINARY production lattice (Papers 1, 7).

Paper 1 sec:convergence (rewritten 2026-09-03) prints the measured binary-lattice
splitting: 37 / 12.7 / 0.65 / 15.6 / 5.1 / 1.7 / 2.7 / 1.7 / 1.1 / 4.1 / 1.8 / 0.39 %
at n_max = 5, 6, 7, 8, 9, 10, 12, 15, 18, 20, 25, 30 (the list printed before --
13 % / 16 % / 0.3 % / 0.005 % -- did not reproduce at 5, 20 or 30).  The decay is
strongly non-monotone and holds on the binary lattice ONLY: on Paper 1's own
CG-magnitude construction A = |T+|+|T-|+|L+|+|L-| the splitting does not decay
(129 / 68 / 84 / 139 % at 8 / 10 / 15 / 20; PM measurement, trunk FULL run #3).

This file computes the splitting from geovac.GeometricLattice's binary graph
Laplacian (the production object AtomicSolver uses), identifying the 2s/2p
eigenstates by FULL-SPECTRUM maximum overlap with |2,0,0> and |2,1,0> (dense
eigh; shift-invert and bottom-k methods mis-identify the states).

Structural claims tested: (i) the splitting is positive and well under 100 %
(perturbation, not force), (ii) it is non-monotone, (iii) it trends downward,
(iv) the n_max = 30 endpoint is sub-percent (slow test).  The xfail records
that the pre-2026-09-03 "16 % at n_max = 10" waypoint does not reproduce.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import diags

from geovac import GeometricLattice


def _splitting_dense(max_n: int) -> float:
    """Relative 2s/2p splitting from L = D - A, identifying 2s/2p eigenstates
    by maximum overlap with bare |2,0,0>, |2,1,0> over the FULL spectrum."""
    lat = GeometricLattice(max_n=max_n)
    A = lat.adjacency
    deg = np.array(A.sum(axis=1)).flatten()
    D = diags(deg, 0, shape=A.shape, format="csr")
    L = (D - A).toarray()

    idx = {s: i for i, s in enumerate(lat.states)}
    i2s, i2p = idx[(2, 0, 0)], idx[(2, 1, 0)]

    w, v = np.linalg.eigh(L)
    j2s = int(np.argmax(np.abs(v[i2s, :])))
    j2p = int(np.argmax(np.abs(v[i2p, :])))
    lam2s, lam2p = w[j2s], w[j2p]
    if lam2s == 0:
        return float("nan")
    return abs(lam2p - lam2s) / abs(lam2s)


# Fast grid (dense eigh up to ~2000 nodes runs in a few seconds total)
N_GRID = [8, 10, 12, 15, 18]


def test_splitting_is_small_and_bounded():
    """At every resolution n_max>=8 the splitting is well under 100%."""
    for n in N_GRID:
        s = _splitting_dense(n)
        assert 0.0 <= s < 1.0, f"n_max={n}: splitting {s*100:.2f}% out of range"


def test_splitting_is_nonmonotone_as_paper_states():
    """Paper 1 itself says the decay 'oscillates'. Verify non-monotonicity:
    at least one n where the splitting INCREASES vs the previous n.
    (A monotone power-law fit would therefore be wrong — the code review
    found exactly this failed fit.)"""
    vals = [_splitting_dense(n) for n in N_GRID]
    increases = [vals[i] > vals[i - 1] for i in range(1, len(vals))]
    assert any(increases), f"expected non-monotone decay; got {vals}"


def test_overall_downward_trend():
    """Late-n splitting much smaller than early-n (artifact decays)."""
    early = _splitting_dense(8)     # ~16%
    late = _splitting_dense(18)     # ~1%
    assert late < early / 3, f"early={early*100:.2f}% late={late*100:.2f}%"


def test_pre_2026_09_03_waypoints_do_not_reproduce():
    """The list Paper 1 printed before 2026-09-03 (13% at 5, 16% at 10, 0.3%
    at 20, 0.005% at 30) does not reproduce: the measured values are 37%,
    1.7%, 4.1% and 0.39%.  Kept as a live guard so the retired waypoints
    cannot return (replaces an xfail that documented a claim the paper no
    longer makes)."""
    assert abs(_splitting_dense(5) - 0.13) > 0.05
    assert abs(_splitting_dense(10) - 0.16) > 0.05


@pytest.mark.slow
def test_splitting_subpercent_by_nmax_30():
    """The convergence-endpoint claim: sub-percent (paper says <0.01%) by
    n_max=30. Dense eigh on 9455 nodes is slow (~3-5 min); marked slow.
    Honest note: we verify <1% (sub-percent), not the paper's <0.01% — the
    precise endpoint magnitude is method-sensitive."""
    s30 = _splitting_dense(30)
    assert s30 < 0.01, f"n_max=30 splitting {s30*100:.4f}% should be <1%"
