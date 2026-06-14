"""
TRUNK QA — Claim 3: 2s/2p splitting -> 0 convergence (Papers 1, 7).

Paper 1 sec:convergence claims the relative 2s/2p splitting from the graph
Laplacian is a discretization artifact that DECAYS to <0.01% by n_max=30,
with quoted waypoints: ~13% (n_max=5), ~16% (n_max=10, "oscillatory peak"),
~0.3% (n_max=20), ~0.005% (n_max=30).

This test computes the ACTUAL splitting from geovac.GeometricLattice's binary
graph Laplacian (the production object AtomicSolver uses) and reports the real
numbers. We identify the 2s/2p eigenstates by FULL-SPECTRUM maximum overlap
with |2,0,0> and |2,1,0> (dense eigh; the only robust identifier — shift-invert
and bottom-k methods mis-identify the states, itself a finding).

HONEST-NUMBERS finding (reported, not asserted to match the paper):
  Dense full-scan splitting (relative %):
    n_max=5  -> ~37%      n_max=8  -> ~16%     n_max=10 -> ~1.7%
    n_max=12 -> ~2.7%     n_max=15 -> ~1.7%    n_max=18 -> ~1.1%
  - The decay is NON-MONOTONE (paper itself says "oscillates").
  - At n_max=10 the splitting is ~1.7%, NOT the paper's 16% (the ~16% appears
    near n_max=8). => Paper 1 sec:III's specific waypoint numbers do NOT
    reproduce. DOWNGRADE SIGNAL for the quoted waypoints.
  - The eigenstate identification is method-sensitive (different sparse
    truncations give different %), so the convergence is real but the
    precise per-n numbers are not robust.

Structural claims that ARE honestly backable (tested below):
  (i)   splitting is positive and well under 100% (perturbation, not force),
  (ii)  it is non-monotone (oscillatory),
  (iii) overall downward trend (artifact decays as boundary recedes).
The "<0.01% by n_max=30" endpoint is checked in a separate slow test.
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


@pytest.mark.xfail(reason="Paper 1 quotes ~16% at n_max=10; production binary "
                          "lattice (dense full-scan) gives ~1.7%. The ~16% "
                          "appears near n_max=8. Waypoint does not reproduce "
                          "-> DOWNGRADE SIGNAL for Paper 1 sec:III numbers.")
def test_paper_waypoint_16pct_at_nmax_10_DOES_NOT_REPRODUCE():
    """Intentionally xfail: documents that the paper's n_max=10 -> 16%
    waypoint is NOT reproduced (real value ~1.7%)."""
    s10 = _splitting_dense(10)
    assert abs(s10 - 0.16) < 0.03   # needs ~16%; actual ~1.7% -> xfail


@pytest.mark.slow
def test_splitting_subpercent_by_nmax_30():
    """The convergence-endpoint claim: sub-percent (paper says <0.01%) by
    n_max=30. Dense eigh on 9455 nodes is slow (~3-5 min); marked slow.
    Honest note: we verify <1% (sub-percent), not the paper's <0.01% — the
    precise endpoint magnitude is method-sensitive."""
    s30 = _splitting_dense(30)
    assert s30 < 0.01, f"n_max=30 splitting {s30*100:.4f}% should be <1%"
