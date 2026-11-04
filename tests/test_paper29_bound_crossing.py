"""Backing for Paper 29's finite-size Kotani-Sunada bound statement.

Paper 29 (abstract + RH-D): every family tested at small sizes satisfies the
Kotani-Sunada irregular bound |mu_nontrivial| <= sqrt(q_max); TWO of four
families cross it at modest sizes (scalar S^5 at N_max = 5, Dirac-S^3 Rule B
at n_max = 4), while the scalar S^3 graph -- whose l-blocks are disconnected
components -- stays below the bound on every component through n_max = 5.

The bound is judged per connected component: each component has its own
Perron value, trivial pair and q_max.  A global recipe that drops only the
largest component's Perron value counts the other components' Perron values
as "non-trivial"; on S^3 at n_max = 5 that manufactured the withdrawn
"+0.198 crossing".

Everything here is recomputed from the LIVE graph constructors with a
non-backtracking builder written in this file (not geovac's hashimoto_matrix),
and then compared with geovac.ihara_zeta.is_ramanujan, so a regression in
either the constructors or the library verdict fails the test.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.csgraph import connected_components

from geovac.ihara_zeta import is_ramanujan
from geovac.lattice import GeometricLattice
from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


def _nb_operator(A: np.ndarray) -> np.ndarray:
    """Non-backtracking (Hashimoto) operator on directed edges, built here."""
    n = A.shape[0]
    darts = [(i, j) for i in range(n) for j in range(n) if A[i, j]]
    idx = {dart: k for k, dart in enumerate(darts)}
    T = np.zeros((len(darts), len(darts)))
    for (u, v), k in idx.items():
        for w in np.nonzero(A[v])[0]:
            if w != u:
                T[k, idx[(v, int(w))]] = 1.0
    return T


def _component_deviations(A) -> list:
    """max|mu_nontrivial| - sqrt(q_max) for every component that has a cycle."""
    A = (np.asarray(A) != 0).astype(int)
    n_comp, labels = connected_components(A, directed=False)
    devs = []
    for c in range(n_comp):
        sel = np.where(labels == c)[0]
        Ac = A[np.ix_(sel, sel)]
        if Ac.sum() // 2 < len(sel):
            continue
        mags = np.abs(np.linalg.eigvals(_nb_operator(Ac)))
        nontriv = mags[mags < mags.max() - 1e-9]
        mu = float(nontriv.max()) if nontriv.size else 0.0
        devs.append(mu - np.sqrt(Ac.sum(axis=1).max() - 1))
    return devs


def _s3(n):
    return GeometricLattice(max_n=n).adjacency.toarray()


def _s5(n):
    return build_bargmann_graph(n).adjacency_dense()


def _dirac(n, rule):
    return build_dirac_s3_graph(n, rule)[0]


# Worst-component deviation at each size, pinned to 1e-3.
PINNED = {
    ("S3", 3): -0.1985, ("S3", 4): -0.0658, ("S3", 5): -0.0287,
    ("S5", 3): -0.3568, ("S5", 4): -0.1615, ("S5", 5): +0.6524,
    ("DirA", 3): -0.1092, ("DirA", 4): -0.0424,
    ("DirB", 3): -0.1188, ("DirB", 4): +1.5301,
}
BUILD = {"S3": _s3, "S5": _s5, "DirA": lambda n: _dirac(n, "A"),
         "DirB": lambda n: _dirac(n, "B")}


@pytest.mark.parametrize("fam,n", sorted(PINNED))
def test_worst_component_deviation_pinned(fam, n):
    """Independent builder and is_ramanujan agree with the pinned value."""
    A = BUILD[fam](n)
    dev = max(_component_deviations(A))
    assert abs(dev - PINNED[(fam, n)]) < 1e-3, f"{fam} {n}: {dev:+.4f}"
    ram, lib_dev, _ = is_ramanujan(A)
    assert abs(lib_dev - dev) < 1e-6, f"{fam} {n}: library {lib_dev:+.4f} vs {dev:+.4f}"
    assert ram is bool(dev <= 1e-9)


def test_s3_is_disconnected_into_l_blocks():
    """The premise of the per-component reading: S^3 has n_max components."""
    for n in (3, 4, 5):
        n_comp, _ = connected_components((_s3(n) != 0).astype(int), directed=False)
        assert n_comp == n


def test_s3_never_crosses_through_n5():
    """No S^3 component crosses at n_max <= 5; is_ramanujan agrees."""
    for n in (3, 4, 5):
        assert all(d < 0.0 for d in _component_deviations(_s3(n))), n
        assert is_ramanujan(_s3(n))[0] is True, n


def test_two_of_four_families_cross():
    """Headline: exactly S^5 (N=5) and Dirac Rule B (n=4) cross; S^3 (n<=5)
    and Dirac Rule A (n<=4) do not; each crosser is sub-Ramanujan at n=3."""
    largest = {"S3": 5, "S5": 5, "DirA": 4, "DirB": 4}
    crossers = {f for f, n in largest.items() if max(_component_deviations(BUILD[f](n))) > 0}
    assert crossers == {"S5", "DirB"}
    for f in crossers:
        assert max(_component_deviations(BUILD[f](3))) < 0.0
