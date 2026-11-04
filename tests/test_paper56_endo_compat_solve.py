"""Solve-based guard for Paper 56 `thm:endo_rigidity` (PS-4): the SHAPE claim.

What the paper claims. End_compat(O_{n_max}) -- the Q-linear endomorphisms
phi of O_{n_max} that intertwine the transitions -- is the block-lower-
triangular subalgebra of End_Q(O_{n_max}) for the shell grading, of
dimension sum_{1<=j<=i<=n_max} (i+1)(j+1) = 4, 19, 55, 125, 245
(n_max = 1..5).

Why this file exists. `tests/test_paper56_endo_rigidity.py` pins the closed
form against the SUM only; nothing solved the compatibility conditions. The
cited driver (`debug/compute_q5p_ps4_endo_rigidity.py`) enumerates the
block-lower-triangular elementary matrices and checks each against a
hand-derived "top-right block is zero" criterion -- it assumes the answer.

What this file does instead. For each n_max it treats every entry of phi
(N x N) and of every partner phi_k (N(k) x N(k), 1 <= k < n_max) as an
unknown, writes the intertwining equations

    P_{n_max,k} . phi  -  phi_k . P_{n_max,k}  =  0     (all 1 <= k < n_max)

with P the PRODUCTION `geovac.pro_system.TransitionMap` matrices, and takes
the exact rational nullspace (sympy DomainMatrix over QQ). It then asserts:

  1. the phi-projection of the solution space has dimension 4, 19, 55, 125,
     245 and the partners phi_k are uniquely determined by phi;
  2. its support (positions nonzero in some solution) is EXACTLY the set of
     (row, col) with shell(row) >= shell(col), shells read from
     `sectors_at_cutoff` -- and the dimension equals the size of that set, so
     the solution space IS the block-lower-triangular coordinate subspace,
     not merely inside it;
  3. it is strictly larger than the block-diagonal subalgebra (the off-
     diagonal compatible maps are really there).

Scope note. The equations are imposed with source cutoff m = n_max and every
target k < n_max, using the algebra-level transitions P_{m,k} -- the reading
of "intertwine all transitions" used by the PS-4 driver. The paper writes
Phi_{m,k} (its Hopf-hom lift symbol); on O_{n_max} itself the transition is
P_{m,k}.
"""
from __future__ import annotations

from typing import List, Set, Tuple

import pytest
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

from geovac.pro_system import N_sectors, TransitionMap, sectors_at_cutoff


# The dimensions an exact solve must return (Paper 56 thm:endo_rigidity).
EXPECTED_DIM = {1: 4, 2: 19, 3: 55, 4: 125, 5: 245}


def _solve_compat(n_max: int) -> Tuple[int, int, List[List]]:
    """Exact nullspace of the intertwining system.

    Returns (dim of joint nullspace, dim of its phi-projection, phi-projected
    basis vectors as length-N^2 lists of QQ elements, row-major phi[r][c]).
    """
    N = N_sectors(n_max)
    ks = list(range(1, n_max))
    # Unknown layout: phi (N*N), then phi_k (N(k)^2) for each k in ks.
    offsets = {}
    pos = N * N
    for k in ks:
        offsets[k] = pos
        pos += N_sectors(k) ** 2
    n_unknowns = pos

    rows: List[List] = []
    for k in ks:
        Nk = N_sectors(k)
        P = TransitionMap(n_max, k).matrix  # Nk x N, production object
        assert P.shape == (Nk, N)
        off = offsets[k]
        # Equation (a, b), a < Nk, b < N:
        #   sum_r P[a, r] phi[r, b]  -  sum_s phi_k[a, s] P[s, b]  = 0
        for a in range(Nk):
            for b in range(N):
                row = [QQ(0)] * n_unknowns
                for r in range(N):
                    if P[a, r] != 0:
                        row[r * N + b] += QQ(int(P[a, r]))
                for s in range(Nk):
                    if P[s, b] != 0:
                        row[off + a * Nk + s] -= QQ(int(P[s, b]))
                rows.append(row)

    if rows:
        A = DomainMatrix(rows, (len(rows), n_unknowns), QQ)
        basis = A.nullspace().to_list()  # rows = basis vectors, QQ entries
    else:  # n_max = 1: no transitions, every endomorphism is compatible
        basis = [[QQ(1) if j == i else QQ(0) for j in range(n_unknowns)]
                 for i in range(n_unknowns)]
    joint_dim = len(basis)

    phi_vecs = [v[: N * N] for v in basis]
    if phi_vecs:
        proj = DomainMatrix([list(v) for v in phi_vecs],
                            (len(phi_vecs), N * N), QQ)
        proj_rank = proj.rank()
    else:
        proj_rank = 0
    return joint_dim, proj_rank, phi_vecs


def _lower_block_positions(n_max: int) -> Set[Tuple[int, int]]:
    """(row, col) with shell(row) >= shell(col), shells from the sector labels."""
    shells = [n for (n, _l) in sectors_at_cutoff(n_max)]
    N = len(shells)
    return {(r, c) for r in range(N) for c in range(N)
            if shells[r] >= shells[c]}


def _block_diag_positions(n_max: int) -> Set[Tuple[int, int]]:
    shells = [n for (n, _l) in sectors_at_cutoff(n_max)]
    N = len(shells)
    return {(r, c) for r in range(N) for c in range(N)
            if shells[r] == shells[c]}


@pytest.mark.parametrize("n_max", sorted(EXPECTED_DIM))
def test_end_compat_solved_dimension_and_shape(n_max: int) -> None:
    N = N_sectors(n_max)
    joint_dim, proj_rank, phi_vecs = _solve_compat(n_max)

    # (1) dimension of the phi-solution space, and partners determined by phi
    assert proj_rank == EXPECTED_DIM[n_max]
    assert joint_dim == proj_rank, "phi_k not uniquely fixed by phi"

    # (2) support is exactly the block-lower-triangular positions
    support = {(idx // N, idx % N)
               for v in phi_vecs for idx, x in enumerate(v) if x != 0}
    lower = _lower_block_positions(n_max)
    assert support == lower
    assert proj_rank == len(lower)  # space = the full coordinate subspace

    # (3) strictly larger than block-diagonal for n_max >= 2
    if n_max >= 2:
        assert proj_rank > len(_block_diag_positions(n_max))
        assert support - _block_diag_positions(n_max)
