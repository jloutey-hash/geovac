"""Guard for the atlas's "host symmetry -> survivor" column (docs/atlas.md), at the two doors where the
corpus can build the released operator in the one-centre (n, l, m) basis.

Door 11 (a second nucleus on the z axis): SO(4) -> SO(2) about the axis.  The cross-centre attraction
mixes l at equal m, and the n = 2 Fock multiplet (4 states) splits by |m| into 1 + 2 + 1; the n = 3
multiplet (9 states) into three singles (m = 0) and three pairs (|m| = 1, 2).
Door 24 (a finite nucleus): SO(4) -> SO(3).  The correction is central, so l and m stay; the n = 2
multiplet splits into 1 + 3 (the s level moves, the p level barely), and n mixes (1s with 2s).

Wrong answers this guard rejects: "l survives at door 11" (the l coupling would vanish); "SO(4) survives
at door 11" (the multiplet would stay degenerate); "SO(4) survives at door 24" (2s and 2p would shift
alike); "n survives at door 24" (1s-2s would not couple).  Not a guard of the m rule itself at door 11:
the cross-centre builder imposes m-diagonality, and the rule's own test is
tests/test_paper58_abelian_residue.py.
"""
from __future__ import annotations

from typing import List

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from sympy.physics.hydrogen import R_nl

from geovac.shibuya_wulfman import compute_cross_center_vne

STATES = [(n, l, m) for n in (1, 2, 3) for l in range(n) for m in range(-l, l + 1)]
IDX = {s: i for i, s in enumerate(STATES)}


def multiplicities(eigs: np.ndarray, tol: float = 1e-8) -> List[int]:
    """Sizes of the degenerate clusters of a sorted eigenvalue list."""
    out: List[int] = []
    last = None
    for e in np.sort(eigs):
        if last is not None and abs(e - last) < tol:
            out[-1] += 1
        else:
            out.append(1)
        last = e
    return out


def _shell_block(V: np.ndarray, n: int) -> np.ndarray:
    ids = [IDX[s] for s in STATES if s[0] == n]
    return V[np.ix_(ids, ids)]


@pytest.fixture(scope="module")
def door11():
    return compute_cross_center_vne(1.0, STATES, 1.0, 2.0, 6)


def test_door11_l_released_at_equal_m(door11):
    V = door11
    off_l = max(abs(V[IDX[a], IDX[b]]) for a in STATES for b in STATES if a[2] == b[2] and a[1] != b[1])
    assert off_l > 1e-3, "the second nucleus must couple different l at equal m"


def test_door11_fock_multiplets_split_by_abs_m(door11):
    V = door11
    assert multiplicities(np.linalg.eigvalsh(_shell_block(V, 2))) == [1, 2, 1]
    pattern = sorted(multiplicities(np.linalg.eigvalsh(_shell_block(V, 3))))
    assert pattern == [1, 1, 1, 2, 2, 2]
    # the splitting is not a rounding artefact
    e2 = np.sort(np.linalg.eigvalsh(_shell_block(V, 2)))
    assert e2[-1] - e2[0] > 1e-2


def _finite_nucleus_element(n1: int, l1: int, n2: int, l2: int, Z: int = 1, R: float = 0.05) -> mp.mpf:
    """<n1 l1 | V_sphere + Z/r | n2 l2> with a uniformly charged sphere of radius R; central, so l1 == l2."""
    if l1 != l2:
        return mp.mpf(0)
    r = sp.symbols("r", positive=True)
    f1 = sp.lambdify(r, R_nl(n1, l1, r, Z), "mpmath")
    f2 = sp.lambdify(r, R_nl(n2, l2, r, Z), "mpmath")

    def dV(x):
        return -Z * (3 / (2 * R) - x * x / (2 * R ** 3)) + Z / x

    return mp.quad(lambda x: f1(x) * f2(x) * dV(x) * x * x, [0, R])


def test_door24_so4_breaks_to_so3_and_n_mixes():
    s2 = _finite_nucleus_element(2, 0, 2, 0)
    p2 = _finite_nucleus_element(2, 1, 2, 1)
    assert s2 > 0 and p2 > 0
    assert s2 > 100 * p2, "the s level must move far more than p: the n = 2 multiplet splits 1 + 3"
    assert abs(_finite_nucleus_element(1, 0, 2, 0)) > 1e-4, "n is released: 1s couples to 2s"
    # the (2l+1) rotation multiplets are untouched: the correction is central, so the pattern is 1 + 3
    shell2 = np.array([float(s2)] + [float(p2)] * 3)
    assert sorted(multiplicities(shell2, tol=1e-12)) == [1, 3]
