"""Paper 60: the M-centre null space is geometry-independent, its RATES are not.

Backs the [MEASURED] scope clause added to `sec:molecular` (CHANGELOG v5.11.4)
after C23 run #2 flagged "spanned by a fixed vector set that does not move with
geometry" plus "a fixed rank-(M-1) rotation removes it" as an over-claim.

The wrong answer this file exists to reject: **"every geometry opens at order
p^2, so one tri(1,2,1) cures any M".**  It does not.  For collinear centres the
order-p^2 form on 1-perp is rank ONE, so only one of the M-1 directions opens at
order 2 and the rest at 4, 6, ..., 2(M-1).

Mechanism (asserted, not just the numbers): j0(pd) = 1 - (pd)^2/6 + O(p^4), so
the order-p^2 behaviour on 1-perp is carried by P D2 P with (D2)_ij = d_ij^2.
Collinear: d_ij^2 = h^2 (i^2 1^T + 1 (j^2)^T - 2 x x^T); P kills the two outer
terms from both sides, leaving -2h^2 P x x^T P, rank one.

Written as a separate pass from the paper edit it protects (CLAUDE.md Sec. 9);
every guard fire-tested via debug/qa/fire_test.py.
"""
from __future__ import annotations

import numpy as np
import pytest

COLLINEAR_3 = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
COLLINEAR_4 = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
EQUILATERAL = [[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3) / 2, 0]]
TETRAHEDRON = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
WATER_BENT = [[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]]


def _j0(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


def _dist(pos) -> np.ndarray:
    p = np.asarray(pos, float)
    return np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)


def _null_orders(pos, p_hi: float = 0.2, p_lo: float = 0.1) -> list[float]:
    """Order in p of each of the M-1 vanishing eigenvalues of [j0(p d_ij)]."""
    D = _dist(pos)
    M = D.shape[0]
    hi = np.sort(np.linalg.eigvalsh(_j0(p_hi * D)))[:M - 1]
    lo = np.sort(np.linalg.eigvalsh(_j0(p_lo * D)))[:M - 1]
    return sorted(float(np.log(hi[k] / lo[k]) / np.log(p_hi / p_lo))
                  for k in range(M - 1))


def _pd2p_rank(pos) -> int:
    D2 = _dist(pos) ** 2
    M = D2.shape[0]
    P = np.eye(M) - np.ones((M, M)) / M
    return int(np.linalg.matrix_rank(P @ D2 @ P, tol=1e-10))


# ---------------------------------------------------- the null SPACE is fixed
@pytest.mark.parametrize("pos", [COLLINEAR_3, COLLINEAR_4, EQUILATERAL,
                                 TETRAHEDRON, WATER_BENT])
def test_null_space_is_the_constants_orthogonal_complement(pos):
    """The half of the paper's claim that IS geometry-independent.

    Rejects: a reading in which the null DIRECTION set moves with geometry.
    At p -> 0 the symbol is the all-ones matrix for every arrangement.
    """
    M = len(pos)
    A = _j0(1e-6 * _dist(pos))
    ev, V = np.linalg.eigh(A)
    small = V[:, :M - 1]
    assert np.max(np.abs(np.ones(M) @ small)) < 1e-6, "null space must be 1-perp"
    assert abs(ev[-1] - M) < 1e-6, "top eigenvalue must be M (all-ones)"


# ---------------------------------------------------- the RATES are not
def test_collinear_geometries_do_not_open_at_order_two():
    """THE LOAD-BEARING GUARD. Rejects: "one tri(1,2,1) cures any M".

    If every direction opened at order 2 this would fail; the orders must be
    2, 4 for three collinear centres and 2, 4, 6 for four.  A preconditioner
    matching a quadratic zero cannot reach a quartic or sextic one.
    """
    assert np.allclose(_null_orders(COLLINEAR_3), [2.0, 4.0], atol=0.05)
    assert np.allclose(_null_orders(COLLINEAR_4), [2.0, 4.0, 6.0], atol=0.08)


@pytest.mark.parametrize("label,pos", [("equilateral", EQUILATERAL),
                                       ("tetrahedron", TETRAHEDRON),
                                       ("water-bent", WATER_BENT)])
def test_non_collinear_geometries_open_entirely_at_order_two(label, pos):
    """Rejects: "the higher orders are generic, so the lever never works".

    They are not generic -- they are the collinear degeneracy.  Water's A_1 is
    in THIS class, which is why the measured table in the paper holds.
    """
    orders = _null_orders(pos)
    assert np.allclose(orders, [2.0] * len(orders), atol=0.05), f"{label}: {orders}"


def test_the_mechanism_is_the_rank_of_the_squared_distance_form():
    """Rejects: "the orders are an empirical fact with no mechanism".

    The paper asserts P D2 P is rank ONE for collinear centres and full rank
    M-1 otherwise; that rank is what predicts how many directions open at
    order 2, so it must be tested and not merely narrated.
    """
    assert _pd2p_rank(COLLINEAR_3) == 1
    assert _pd2p_rank(COLLINEAR_4) == 1
    assert _pd2p_rank(EQUILATERAL) == 2
    assert _pd2p_rank(WATER_BENT) == 2
    assert _pd2p_rank(TETRAHEDRON) == 3
    for pos in (COLLINEAR_3, COLLINEAR_4, EQUILATERAL, TETRAHEDRON, WATER_BENT):
        n_order2 = sum(1 for o in _null_orders(pos) if abs(o - 2.0) < 0.05)
        assert n_order2 == _pd2p_rank(pos), (
            "the count of order-2 directions must equal rank(P D2 P)")
