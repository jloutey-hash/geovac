"""
Bit-exact verification of the Layer-1 spin-structure moduli + chirality obstruction.

Locks the findings of debug/spin_structure_moduli_memo.md:
  (1) scalar Fock-graph cycle rank beta_1 = sum_l (n_max-l-1)(2l);
  (2) Dirac magnetic reflection m_j -> -m_j is a fixed-point-free (Kramers) graph
      automorphism of the Dirac-S^3 (DiracLattice) graph;
  (3) the genuine-spin (chirality) reflection has NO node bijection: the two
      chirality sectors differ in size by exactly +2 per orbital block, total
      imbalance n_max(n_max+1) -- the discrete Dirac-index obstruction.

All claims are pure-integer / combinatorial (pi-free, Layer-1).
"""
import numpy as np
import pytest

from geovac.fock_graph_hodge import FockGraphHodge
from geovac.dirac_lattice import DiracLattice


# ---------------------------------------------------------------- scalar beta_1
@pytest.mark.parametrize("n_max,expected", [(2, 0), (3, 2), (4, 8)])
def test_scalar_betti1_matches_closed_form(n_max, expected):
    closed = sum((n_max - l - 1) * (2 * l) for l in range(n_max))
    assert closed == expected
    assert FockGraphHodge(n_max).betti_1 == expected


# ---------------------------------------------- Dirac magnetic reflection (Kramers)
@pytest.mark.parametrize("n_max", [2, 3])
def test_dirac_mj_reflection_is_fixed_point_free_automorphism(n_max):
    dl = DiracLattice(n_max, mode="atomic")
    labels = [(lab.n_fock, lab.kappa, lab.two_m_j) for lab in dl.labels]
    idx = {lab: i for i, lab in enumerate(labels)}
    A = dl.adjacency.toarray()
    A = (A != 0).astype(np.int64)

    # m_j -> -m_j is a node permutation (range is symmetric in two_m_j)
    P = np.array([idx[(n, k, -tm)] for (n, k, tm) in labels], dtype=np.int64)

    # fixed-point-free: half-integer m_j is never 0, so two_m_j is odd != -two_m_j
    assert all(P[i] != i for i in range(len(labels)))
    assert all(tm % 2 != 0 for (_, _, tm) in labels)

    # graph automorphism: A[P,P] == A
    assert np.array_equal(A[np.ix_(P, P)], A)


# ---------------------------------------------- chirality obstruction (no bijection)
@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_chirality_sectors_have_no_bijection(n_max):
    dl = DiracLattice(n_max, mode="atomic")
    labels = [(lab.n_fock, lab.kappa, lab.two_m_j) for lab in dl.labels]

    n_pos = sum(1 for (_, k, _) in labels if k < 0)   # chi=+1, j=l+1/2
    n_neg = sum(1 for (_, k, _) in labels if k > 0)   # chi=-1, j=l-1/2

    # no equal-size pairing -> no chirality-reflecting node permutation
    assert n_pos != n_neg
    # global imbalance is exactly n_max(n_max+1) = 2 * (#orbital blocks)
    assert n_pos - n_neg == n_max * (n_max + 1)


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_chirality_block_mismatch_is_plus_two_per_block(n_max):
    dl = DiracLattice(n_max, mode="atomic")
    from collections import defaultdict
    blk = defaultdict(lambda: {"neg": 0, "pos": 0})
    for lab in dl.labels:
        k = lab.kappa
        l = abs(k) - 1 if k < 0 else k
        blk[(lab.n_fock, l)]["neg" if k < 0 else "pos"] += 1

    for (n, l), v in blk.items():
        # chi=+1 (kappa<0, j=l+1/2): 2l+2 states; chi=-1 (kappa>0, j=l-1/2): 2l states
        assert v["neg"] == 2 * l + 2
        assert v["pos"] == 2 * l
        assert v["neg"] - v["pos"] == 2
