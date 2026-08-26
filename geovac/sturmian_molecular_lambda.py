"""Paper 60 -- molecular block-encoding 1-norm for two-electron H2 in the two-center
shared-scale Coulomb--Sturmian basis.

This module promotes the Paper-60 sprint driver (``debug/sturmian_h2_ci_1norm.py``) to
tracked code so the paper's honest NEGATIVE resource result is backed by a regression
test rather than a transient ``debug/`` script.  The STANDARD block-encoding 1-norm

    lambda = sum_pq |h_pq| + sum_pqrs |(pq|rs)|

computed in the Loewdin-orthonormalized molecular-Sturmian basis grows POLYNOMIALLY,
``lambda ~ n_orb^2.2``, with NO sublinear behavior.  The atomic
isoenergetic-secular-matrix sublinearity of Paper 60 ``eq:sublinear`` is a single-center
property and does NOT transfer to molecules (``sec:manyelectron``); the standard lambda
below is what actually sets qubitization cost and it is measurable.

Imports the tracked two-center integral engine
:class:`geovac.sturmian_integrals.GoscinskianIntegrals` (does NOT re-implement it).
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np

from geovac.sturmian_integrals import GoscinskianIntegrals, Orb


def build_raw(gi: GoscinskianIntegrals, nmax: int, zeta: float
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One-electron ``h``, overlap ``S`` and ERI tensor for shared-scale s-Sturmians.

    Builds the raw (non-orthogonal) matrices for ``2*nmax`` shared-scale (decay ``zeta``)
    s-Sturmians: principal numbers ``n = 1..nmax`` on each of centers ``A`` and ``B``.
    ``h = T + V`` uses the shared-scale kinetic operator (``kscale = zeta``) and the
    two-center nuclear attraction ``-1/r_A - 1/r_B``.

    Returns ``(h, S, eri)`` with ``h``, ``S`` of shape ``(N, N)`` and ``eri`` of shape
    ``(N, N, N, N)`` in chemist notation ``(ij|kl)``, where ``N = 2*nmax``.
    """
    orbs: List[Orb] = ([(n, "A", zeta) for n in range(1, nmax + 1)]
                       + [(n, "B", zeta) for n in range(1, nmax + 1)])
    N = len(orbs)
    S = np.zeros((N, N))
    h = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            S[i, j] = S[j, i] = gi.overlap(orbs[i], orbs[j])
            hij = gi.kinetic(orbs[i], orbs[j], zeta) + gi.nuclear(orbs[i], orbs[j])
            h[i, j] = h[j, i] = hij
    eri = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(i, N):
            for k in range(N):
                for l in range(k, N):
                    v = gi.eri(orbs[i], orbs[j], orbs[k], orbs[l])
                    for (a, b) in ((i, j), (j, i)):
                        for (c, d) in ((k, l), (l, k)):
                            eri[a, b, c, d] = v
    return h, S, eri


def loewdin(h: np.ndarray, S: np.ndarray, eri: np.ndarray
            ) -> Tuple[np.ndarray, np.ndarray]:
    """Canonical Loewdin (``S^{-1/2}``) orthonormalization of ``h`` and the ERI tensor.

    Drops the null space of ``S`` (eigenvalues ``<= 1e-9``) and rotates the one- and
    two-electron integrals into the orthonormal molecular-orbital basis.  Returns
    ``(hm, em)`` -- the orthonormalized one-electron matrix and ERI tensor.
    """
    w, U = np.linalg.eigh(S)
    keep = w > 1e-9
    X = U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep]))       # canonical S^{-1/2}, drops null
    hm = X.T @ h @ X
    em = np.einsum("ip,jq,kr,ls,ijkl->pqrs", X, X, X, X, eri, optimize=True)
    return hm, em


def h2_lambda(nmax: int, R: float = 1.4, zeta: float = 1.2,
              Lmax: int = 20, nr: int = 2200, nth: int = 150,
              rmax: float = 55.0) -> Tuple[int, float]:
    """Standard block-encoding 1-norm ``lambda`` for two-electron H2 at basis size ``nmax``.

    Builds the shared-scale two-center Sturmian basis (``n_orb = 2*nmax`` orbitals),
    Loewdin-orthonormalizes, and returns ``(n_orb, lambda)`` with

        lambda = sum_pq |h_pq| + sum_pqrs |(pq|rs)|.

    The grid parameters ``(Lmax, nr, nth, rmax)`` are forwarded to
    :class:`~geovac.sturmian_integrals.GoscinskianIntegrals`.  Defaults reproduce the
    Paper-60 driver; a regression test uses a coarser grid (the exponent is stable).
    """
    gi = GoscinskianIntegrals(R=R, Lmax=Lmax, nr=nr, nth=nth, rmax=rmax)
    h, S, eri = build_raw(gi, nmax, zeta)
    hm, em = loewdin(h, S, eri)
    lam = float(np.abs(hm).sum() + np.abs(em).sum())
    return 2 * nmax, lam


def lambda_scaling(nmax_values: Tuple[int, ...] = (1, 2, 3), R: float = 1.4,
                   zeta: float = 1.2, Lmax: int = 20, nr: int = 2200,
                   nth: int = 150, rmax: float = 55.0
                   ) -> Tuple[List[int], List[float], float]:
    """Sweep ``nmax`` and fit ``lambda ~ n_orb^p`` (log-log slope).

    Returns ``(n_orb_list, lambda_list, p)``.  A polynomial (non-sublinear) result has
    ``p`` well above 1; the Paper-60 measurement is ``p ~ 2.2``.
    """
    n_orbs: List[int] = []
    lams: List[float] = []
    for nmax in nmax_values:
        n_orb, lam = h2_lambda(nmax, R=R, zeta=zeta, Lmax=Lmax,
                               nr=nr, nth=nth, rmax=rmax)
        n_orbs.append(n_orb)
        lams.append(lam)
    p = float(np.polyfit(np.log(n_orbs), np.log(lams), 1)[0])
    return n_orbs, lams, p


if __name__ == "__main__":
    orbs, lams, p = lambda_scaling()
    for n_orb, lam in zip(orbs, lams):
        print(f"  n_orb={n_orb:>2}  lambda={lam:>12.4f}")
    print(f"\n  lambda ~ n_orb^{p:.2f}")
