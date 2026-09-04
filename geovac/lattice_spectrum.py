"""Closed-form and block-decomposed spectra for the GeoVac production lattice.

WHY THIS MODULE EXISTS
----------------------
The production lattice's Laplacian spectrum has been *proved* closed form
(Paper 0 SS VI, Paper 1 SS III):  no edge changes ``l``, so ``L`` is block
diagonal in ``l`` and each block is the Cartesian product of two path graphs,

    L_l  =  P_{n_max - l}  x  P_{2l + 1},

whose Laplacian eigenvalues add:

    spec(L_l) = { 2 - 2 cos(j pi / (n_max - l))  +  2 - 2 cos(k pi / (2l + 1)) }
                for j = 0 .. n_max-l-1,  k = 0 .. 2l.

Until 2026-09-04 that closed form lived as a private helper *inside the test
that proves it*, so no other test could use it -- and fourteen other test
files were still calling ``numpy.linalg.eigh`` on the dense operator.  At
n_max = 30 that is a 9455 x 9455 dense eigendecomposition, measured at 45-69 s,
to obtain numbers this module returns in milliseconds.

MEASURED (n_max = 30, 9455 x 9455)

    route                      time        agreement vs dense eigh
    dense eigh                 69.2 s      --
    block_eigh   (eigenpairs)   1.05 s     max|dlambda| = 1.4e-14   -> 66x
    spectrum     (eigenvalues)  1.7 ms     max|dlambda| = 1.4e-13   -> ~27000x

Both routes are EXACT, not approximations:  ``block_eigh`` diagonalises the
very same operator one block at a time, and ``spectrum`` evaluates the proved
closed form.

DEGENERATE-BASIS CAVEAT -- READ BEFORE USING EIGENVECTORS
---------------------------------------------------------
Eigenvalues are basis-free;  eigenVECTORS at a degenerate eigenvalue are not.
``lambda_2s = 3`` has multiplicity 25 at n_max = 30, so *any* routine --
LAPACK, the block route, or the closed form -- returns an arbitrary basis of
that eigenspace, and a quantity read off individual eigenvectors can change
under a relabelling that is mathematically a no-op.  (Measured during
/qa DELTA #5:  permuting the node labelling turned a reported 2.70% into
4421%.)  So:

  * safe:    eigenvalues;  projector-based quantities such as ||P_lambda e||.
  * unsafe:  argmax over eigenvector components;  "the" eigenvector at an
             eigenvalue;  diagonal fractions in an eigenbasis.

If you need a basis-free overlap, use :func:`eigenspace_overlap`.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np


def block_dims(n_max: int, l: int) -> Tuple[int, int]:
    """Path-graph factor sizes of the ``l``-block:  P_{n_max-l} x P_{2l+1}."""
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    if not 0 <= l < n_max:
        raise ValueError(f"l must satisfy 0 <= l < n_max={n_max}, got {l}")
    return n_max - l, 2 * l + 1


def path_spectrum(m: int) -> np.ndarray:
    """Laplacian spectrum of the path graph ``P_m``:  2 - 2 cos(j pi / m)."""
    if m < 1:
        raise ValueError(f"path length must be >= 1, got {m}")
    return np.array([2.0 - 2.0 * math.cos(j * math.pi / m) for j in range(m)])


def block_spectrum(n_max: int, l: int) -> np.ndarray:
    """Closed-form spectrum of the ``l``-block, sorted ascending."""
    m1, m2 = block_dims(n_max, l)
    a, b = path_spectrum(m1), path_spectrum(m2)
    return np.sort((a[:, None] + b[None, :]).ravel())


def spectrum(n_max: int) -> np.ndarray:
    """Closed-form spectrum of the whole production-lattice Laplacian.

    Sorted ascending, length ``n_max (n_max+1) (2 n_max+1) / 6``.  Agrees with
    a dense ``eigh`` to ~1e-13 and is roughly four orders of magnitude faster
    (see the module docstring).
    """
    return np.sort(np.concatenate(
        [block_spectrum(n_max, l) for l in range(n_max)]))


def lambda_max(n_max: int) -> float:
    """Top Laplacian eigenvalue, in closed form.

    Saturates the bipartite bound ``2 d_max = 8`` from below, and is attained
    in a mid-``l`` block rather than the s-wave block (whose own bound is 4).
    """
    # O(n_max), not O(n_max^3):  the top of a Cartesian product is the sum of
    # the factor tops, so there is no reason to materialise a block spectrum.
    # (The first version of this function built every block in full and could
    # not finish at n_max = 5000;  the profile caught it.)
    best = 0.0
    for l in range(n_max):
        m1, m2 = block_dims(n_max, l)
        top = ((2.0 - 2.0 * math.cos((m1 - 1) * math.pi / m1))
               + (2.0 - 2.0 * math.cos((m2 - 1) * math.pi / m2)))
        if top > best:
            best = top
    return float(best)


def path_eigenvectors(m: int) -> np.ndarray:
    """Orthonormal Laplacian eigenvectors of ``P_m``, columns ordered by j.

    ``v_j(i) = cos(j pi (i + 1/2) / m)``, normalised.
    """
    if m < 1:
        raise ValueError(f"path length must be >= 1, got {m}")
    V = np.array([[math.cos(j * math.pi * (i + 0.5) / m) for j in range(m)]
                  for i in range(m)])
    return V / np.linalg.norm(V, axis=0)


def block_index(states) -> Dict[int, np.ndarray]:
    """Map ``l`` -> array of row indices of ``states`` in that block.

    ``states`` is the lattice's ``(n, l, m)`` list;  membership is read from
    the ``l`` column, which is what makes ``L`` block diagonal.
    """
    st = np.asarray(states)
    return {int(l): np.where(st[:, 1] == l)[0] for l in np.unique(st[:, 1])}


def block_eigh(L, states) -> Tuple[np.ndarray, List[Tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    """Diagonalise ``L`` one ``l``-block at a time.

    Exact:  the same operator, decomposed along the block structure it
    already has.  Returns ``(eigenvalues_sorted, blocks)`` where each block is
    ``(row_indices, block_eigenvalues, block_eigenvectors)``.

    Use this when eigenVECTORS are needed -- it is ~66x faster than a dense
    ``eigh`` at n_max = 30 -- and read the degenerate-basis caveat in the
    module docstring before using the vectors themselves.
    """
    Ld = L.toarray() if hasattr(L, "toarray") else np.asarray(L)
    blocks, evs = [], []
    for l, sel in sorted(block_index(states).items()):
        if sel.size == 0:
            continue
        wl, vl = np.linalg.eigh(Ld[np.ix_(sel, sel)])
        blocks.append((sel, wl, vl))
        evs.append(wl)
    return np.sort(np.concatenate(evs)), blocks


def lambda_max_from_operator(L, states) -> float:
    """Top eigenvalue of the REAL constructed operator, computed blockwise.

    This is the cheap end-to-end route:  it still diagonalises the lattice
    that was actually built, so it does not rest on the closed form -- it is
    an independent check ON the closed form.  Because ``L`` is block diagonal
    in ``l``, the top eigenvalue is the largest of the per-block tops, and
    each block is tiny compared with the whole.

    MEASURED at n_max = 70 (116,795 nodes, largest block 2485):
        global  eigsh(k=1)   70.07 s
        blockwise            0.93 s      -> 75x, agreement 2.1e-14

    Use this where a test should verify the constructed graph;  use
    :func:`lambda_max` where the closed form is the intended subject.
    """
    from scipy.sparse.linalg import eigsh

    best = -np.inf
    for _l, sel in sorted(block_index(states).items()):
        if sel.size == 0:
            continue
        sub = L[sel][:, sel] if hasattr(L, "tocsr") else np.asarray(L)[np.ix_(sel, sel)]
        if sel.size <= 3:
            dense = sub.toarray() if hasattr(sub, "toarray") else sub
            top = float(np.max(np.linalg.eigvalsh(dense)))
        else:
            top = float(eigsh(sub, k=1, which="LA")[0][0])
        best = max(best, top)
    return float(best)


def eigenspace_overlap(blocks, row: int, tol: float = 1e-9) -> List[Tuple[float, float]]:
    """Basis-free overlap of basis vector ``e_row`` with each eigenSPACE.

    Returns ``[(eigenvalue, ||P_lambda e_row||), ...]`` sorted by descending
    overlap.  Eigenvalues within ``tol`` are pooled into one eigenspace first,
    which is what makes the result invariant under the arbitrary choice of
    basis inside a degenerate eigenspace -- unlike an argmax over individual
    eigenvector components.
    """
    for sel, wl, vl in blocks:
        hit = np.where(sel == row)[0]
        if hit.size == 0:
            continue
        amps = vl[int(hit[0]), :]
        order = np.argsort(wl)
        w_s, a_s = wl[order], amps[order]
        out, i = [], 0
        while i < len(w_s):
            j = i
            while j + 1 < len(w_s) and abs(w_s[j + 1] - w_s[i]) <= tol:
                j += 1
            out.append((float(w_s[i]), float(np.linalg.norm(a_s[i:j + 1]))))
            i = j + 1
        return sorted(out, key=lambda t: -t[1])
    raise ValueError(f"row {row} is in no block")
