"""The GeoVac graph Laplacian does NOT converge to the S^3 Laplacian; the label set does carry it.

Backs the 2026-09-25 re-pricing (Papers 0, 1, 7 and their citers):

  WALL   The production edge Laplacian L = D - A of ``GeometricLattice`` (binary default AND the
         ``topological_weights`` = 1/(n_i n_j) option) has no edge that changes l.  Hence it is the
         direct sum of n_max disconnected l-blocks, dim ker L = n_max (the S^3 Laplace-Beltrami
         operator has ONE zero mode), and the heat kernel e^{-tL} is exactly zero between blocks.
         No rescaling of such an operator can converge to Delta_{S^3}, whose heat kernel is
         everywhere positive.

  EXACT  On the node set (n, l, m) -- the Gel'fand-Tsetlin labels of the Peter-Weyl truncation of
         L^2(S^3) -- the S^3 Laplacian acts as the SO(4) Casimir L^2 + A^2 = n^2 - 1 (A = the
         energy-scaled Runge-Lenz vector), diagonal with multiplicity n^2 on EVERY shell including
         the top one, because each shell is a complete irreducible representation.

What wrong answer each check rejects (fire-tested by debug/firetest_graph_s3.py):
  * (a)/(b) reject a lattice with ANY l-changing edge (such an edge merges blocks, so the
    component count, the kernel dimension, the cross-l adjacency count and the cross-l heat
    kernel all change).
  * (c) rejects a wrong Runge-Lenz reduced matrix element: the generators are BUILT here from
    the standard reduced matrix elements (A_z from Pauli's formula, A_+- from the vector-operator
    commutator [L_+-, A_z] = -+ A_+-), squared and summed.  The target diag(n^2 - 1) is never
    used to construct anything.  A wrong element breaks the Casimir (and the so(4) algebra).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from typing import Dict, List, Tuple

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from geovac.lattice import GeometricLattice

State = Tuple[int, int, int]


def _laplacian(n_max: int, topological: bool) -> Tuple[np.ndarray, np.ndarray, List[State]]:
    lat = GeometricLattice(max_n=n_max, topological_weights=topological)
    adj = lat.adjacency.toarray()
    lap = np.diag(adj.sum(axis=1)) - adj
    return adj, lap, list(lat.states)


# --------------------------------------------------------------------------- (a) components / kernel
@pytest.mark.parametrize("topological", [False, True], ids=["binary", "topological_weights"])
@pytest.mark.parametrize("n_max", list(range(2, 13)))
def test_graph_splits_into_nmax_l_blocks(n_max: int, topological: bool) -> None:
    adj, lap, states = _laplacian(n_max, topological)
    ls = np.array([l for (_, l, _) in states])

    # no edge between different l
    cross = (ls[:, None] != ls[None, :]) & (np.abs(adj) > 0)
    assert int(cross.sum()) == 0, f"{int(cross.sum())} adjacency entries join different l"

    # exactly n_max connected components ...
    n_comp, _ = connected_components(adj != 0, directed=False)
    assert n_comp == n_max

    # ... and an n_max-dimensional kernel (S^3 Laplacian: one zero mode)
    ev = np.linalg.eigvalsh(lap)
    tol = 1e-9 * max(1.0, float(np.abs(ev).max()))
    assert int(np.sum(np.abs(ev) < tol)) == n_max


# --------------------------------------------------------------------------- (b) heat kernel
@pytest.mark.parametrize("topological", [False, True], ids=["binary", "topological_weights"])
@pytest.mark.parametrize("n_max", [3, 6, 9])
def test_heat_kernel_vanishes_between_l_blocks(n_max: int, topological: bool) -> None:
    _, lap, states = _laplacian(n_max, topological)
    ls = np.array([l for (_, l, _) in states])
    K = expm(-1.0 * lap)
    cross = ls[:, None] != ls[None, :]
    same = ~cross
    assert float(np.abs(K[cross]).max()) == 0.0, "heat kernel leaks between l-blocks"
    # control: within a block the kernel is strictly positive (connected block), so the zero
    # above is a statement about disconnection, not an underflowed kernel
    assert float(K[same].min()) > 0.0


# --------------------------------------------------------------------------- (c) SO(4) Casimir
def _basis(n_max: int) -> Tuple[List[State], Dict[State, int]]:
    states = [(n, l, m) for n in range(1, n_max + 1) for l in range(n) for m in range(-l, l + 1)]
    return states, {s: i for i, s in enumerate(states)}


def _az_element(n: int, l: int, m: int) -> float:
    """<n, l+1, m | A_z | n, l, m> for the energy-scaled Runge-Lenz vector (Pauli 1926).

    A = n * (Runge-Lenz) in atomic units, so that L^2 + A^2 = n^2 - 1 on shell n.
    """
    return float(np.sqrt((n * n - (l + 1) ** 2) * ((l + 1) ** 2 - m * m) / ((2 * l + 1) * (2 * l + 3))))


def so4_generators(n_max: int) -> Tuple[List[np.ndarray], List[np.ndarray], List[State]]:
    """Build (L_x, L_y, L_z) and (A_x, A_y, A_z) on the (n, l, m) basis from reduced elements."""
    states, idx = _basis(n_max)
    d = len(states)
    Lz = np.zeros((d, d), complex)
    Lp = np.zeros((d, d), complex)
    Az = np.zeros((d, d), complex)
    for (n, l, m), i in idx.items():
        Lz[i, i] = m
        if m < l:  # Condon-Shortley L_+
            Lp[idx[(n, l, m + 1)], i] = np.sqrt(l * (l + 1) - m * (m + 1))
        if l + 1 < n and abs(m) <= l:
            j = idx[(n, l + 1, m)]
            a = _az_element(n, l, m)
            Az[j, i] = a
            Az[i, j] = a
    Lm = Lp.conj().T
    Lx, Ly = (Lp + Lm) / 2, (Lp - Lm) / 2j
    # A is a vector operator:  [L_+-, A_z] = -+ A_+-
    Ap = -(Lp @ Az - Az @ Lp)
    Am = Lm @ Az - Az @ Lm
    Ax, Ay = (Ap + Am) / 2, (Ap - Am) / 2j
    return [Lx, Ly, Lz], [Ax, Ay, Az], states


@pytest.mark.parametrize("n_max", list(range(1, 9)))
def test_so4_casimir_is_n2_minus_1_on_every_shell(n_max: int) -> None:
    Ls, As, states = so4_generators(n_max)
    C = sum(X @ X for X in Ls) + sum(X @ X for X in As)
    target = np.diag([n * n - 1 for (n, _, _) in states]).astype(complex)
    assert np.allclose(C, target, atol=1e-10), f"max dev {np.abs(C - target).max():.3e}"

    # multiplicity n^2 of eigenvalue n^2 - 1 on every shell, the top shell included
    ev = np.round(np.linalg.eigvalsh(C).real, 8)
    for n in range(1, n_max + 1):
        assert int(np.sum(ev == n * n - 1)) == n * n

    # the generators close into so(4):  [A_x, A_y] = i L_z  and  A . L = 0
    Ax, Ay, _ = As
    assert np.allclose(Ax @ Ay - Ay @ Ax, 1j * Ls[2], atol=1e-10)
    assert np.allclose(sum(a @ l for a, l in zip(As, Ls)), 0.0, atol=1e-10)


def test_casimir_is_not_an_edge_laplacian() -> None:
    """The production L does not share the Casimir's spectrum: at n_max = 8 every eigenvalue
    of L lies below 8 (L is bounded in [0, 8]), while the Casimir reaches n_max^2 - 1 = 63."""
    n_max = 8
    _, lap, states = _laplacian(n_max, topological=False)
    ev = np.linalg.eigvalsh(lap)
    assert ev.max() < 8.0
    assert max(n * n - 1 for (n, _, _) in states) == n_max ** 2 - 1 > 8
