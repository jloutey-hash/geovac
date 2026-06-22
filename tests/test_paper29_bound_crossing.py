"""Genuine backing for Paper 29's Alon-Boppana / Kotani-Sunada bound-crossing claim.

Paper 29 (abstract + RH-D): the GeoVac graph families are Ramanujan
(satisfy the Kotani-Sunada irregular bound |mu_nontrivial| <= sqrt(q_max))
at small sizes, but three of four families CROSS the bound at modest sizes
(V ~ 30-60) -- so the graph-RH for GeoVac is a FINITE-SIZE statement, NOT
an asymptotic one (the dense families fail asymptotically).

This test RECOMPUTES the Hashimoto spectrum and the sub-Ramanujan deviation
    deviation = max|mu_nontrivial| - sqrt(q_max)
directly from the LIVE graph constructors (GeometricLattice / build_bargmann_graph
/ build_dirac_s3_graph) at each size -- it does NOT read the frozen
debug/data/alon_boppana_sweep.json. A regression in any graph constructor or
in hashimoto_matrix changes these signs and fails the test.

The load-bearing content is the finite-size SIGN FLIP (Ramanujan at small
size -> crosses at modest size) for three families, which is exactly the
"finite-size, not asymptotic" claim.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.ihara_zeta import hashimoto_matrix, is_ramanujan
from geovac.lattice import GeometricLattice
from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


def _deviation(A) -> float:
    """sub-Ramanujan deviation = max|mu_nontrivial| - sqrt(q_max), recomputed
    from the Hashimoto spectrum (mirrors debug/compute_alon_boppana_sweep.py)."""
    A = (np.asarray(A) != 0).astype(int)
    deg = A.sum(axis=1).astype(int)
    q_max = int(deg.max()) - 1
    T = hashimoto_matrix(A).astype(float)
    mags = np.abs(np.linalg.eigvals(T))
    rho = float(mags.max())
    nontriv = mags[mags < rho - 1e-9]
    mu_nt = float(nontriv.max()) if nontriv.size else 0.0
    return mu_nt - np.sqrt(max(q_max, 0))


def _bargmann_A(n_max):
    return build_bargmann_graph(n_max).adjacency_dense()


def _dirac_A(n_max, rule):
    A, _labels, _deg, _desc = build_dirac_s3_graph(n_max, rule)
    return A


# ---------------------------------------------------------------------------
# Ramanujan side (small size): deviation < 0, is_ramanujan True  (fast)
# ---------------------------------------------------------------------------

def test_s5_bargmann_small_is_ramanujan():
    """S5 Bargmann N_max=3 (V=20) satisfies the bound (deviation < 0)."""
    A = _bargmann_A(3)
    dev = _deviation(A)
    assert dev < 0.0, f"expected sub-Ramanujan (dev<0), got {dev:+.4f}"
    ram, _, _ = is_ramanujan(A)
    assert ram is True


def test_dirac_B_small_is_ramanujan():
    """Dirac Rule-B n_max=3 (V=28) satisfies the bound (deviation < 0)."""
    A = _dirac_A(3, "B")
    dev = _deviation(A)
    assert dev < 0.0, f"expected sub-Ramanujan (dev<0), got {dev:+.4f}"
    ram, _, _ = is_ramanujan(A)
    assert ram is True


def test_s3_coulomb_small_is_ramanujan():
    """S3 Coulomb max_n=4 (V=30) satisfies the bound (deviation < 0)."""
    A = GeometricLattice(max_n=4).adjacency.toarray()
    dev = _deviation(A)
    assert dev < 0.0, f"expected sub-Ramanujan (dev<0), got {dev:+.4f}"


# ---------------------------------------------------------------------------
# Crossing side (modest size): deviation > 0, is_ramanujan False  -> finite-size
# sign flip.  (slow: V ~ 55-60, Hashimoto 2E-dim eigensolve)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_dirac_B_crosses_at_n4():
    """Dirac Rule-B CROSSES at n_max=4 (V=60): deviation flips positive.
    Pairs with test_dirac_B_small_is_ramanujan -> finite-size sign flip."""
    dev_small = _deviation(_dirac_A(3, "B"))
    dev_cross = _deviation(_dirac_A(4, "B"))
    assert dev_small < 0.0 < dev_cross, \
        f"expected sign flip n3<0<n4, got {dev_small:+.4f}, {dev_cross:+.4f}"
    ram, _, _ = is_ramanujan(_dirac_A(4, "B"))
    assert ram is False


@pytest.mark.slow
def test_s5_bargmann_crosses_at_N5():
    """S5 Bargmann CROSSES at N_max=5 (V=56): deviation flips positive.
    Pairs with test_s5_bargmann_small_is_ramanujan -> finite-size sign flip."""
    dev_small = _deviation(_bargmann_A(3))
    dev_cross = _deviation(_bargmann_A(5))
    assert dev_small < 0.0 < dev_cross, \
        f"expected sign flip N3<0<N5, got {dev_small:+.4f}, {dev_cross:+.4f}"


@pytest.mark.slow
def test_s3_coulomb_crosses_at_n5():
    """S3 Coulomb CROSSES at max_n=5 (V=55): deviation flips positive.
    Pairs with test_s3_coulomb_small_is_ramanujan -> finite-size sign flip."""
    dev_small = _deviation(GeometricLattice(max_n=4).adjacency.toarray())
    dev_cross = _deviation(GeometricLattice(max_n=5).adjacency.toarray())
    assert dev_small < 0.0 < dev_cross, \
        f"expected sign flip n4<0<n5, got {dev_small:+.4f}, {dev_cross:+.4f}"


@pytest.mark.slow
def test_three_of_four_families_cross():
    """The headline: three of four families cross the bound at modest sizes
    (V ~ 30-60), recomputed from live graphs -- the finite-size (not
    asymptotic) Ramanujan statement of Paper 29."""
    crossers = {
        "S3_Coulomb": _deviation(GeometricLattice(max_n=5).adjacency.toarray()),
        "S5_Bargmann": _deviation(_bargmann_A(5)),
        "Dirac_B": _deviation(_dirac_A(4, "B")),
    }
    for fam, dev in crossers.items():
        assert dev > 0.0, f"{fam} expected to cross (dev>0), got {dev:+.4f}"
    assert sum(1 for d in crossers.values() if d > 0) == 3
