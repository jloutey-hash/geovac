"""Double-factorization (DF) 1-norm of a molecular electronic Hamiltonian.

Definition (Lee et al., PRX Quantum 2, 030305 (2021), arXiv:2011.03494, Appendix C,
following von Burg et al., PRR 3, 033055 (2021), arXiv:2007.14460):

    H = sum_{pq,s} h_pq a+_ps a_qs + 1/2 sum_{pqrs,st} (pq|rs) a+_ps a+_rt a_st a_qs + ecore

with real orthonormal spatial orbitals and chemist-ordered (pq|rs).

* Eq. (3):   T_pq  = h_pq - 1/2 sum_r (pr|rq)      (one-body part of Lee's Eq. (1) form)
* Eq. (A9):  T'_pq = T_pq + sum_r (pq|rr)          (two-body correction folded in, Eq. (C5)-(C9))
* lambda_T' = Schatten 1-norm of T' (sum |eig T'|; text after Eq. (C11), Table XIII caption)
* Eq. (C1):  (pq|rs) = sum_l W^(l)_pq W^(l)_rs  (first factorization; here, as in Lee's
             Table XIII data and openfermion.resource_estimates.df, the eigendecomposition of
             the N^2 x N^2 supermatrix, leaves in descending eigenvalue order)
* second factorization: W^(l) = sum_m f^(l)_m u_m u_m^T  (eigendecomposition, item 2 before C41)
* Eq. (C41): drop eigenvector m of leaf l when (sum_p |f^(l)_p|) |f^(l)_m| < threshold;
             if a leaf keeps nothing, stop and drop all later leaves.
* Eq. (C11): lambda_DF = 1/4 sum_l (sum_{m kept} |f^(l)_m|)^2
* lambda = lambda_T' + lambda_DF  (identity / constant shifts are not counted).

lambda_T' uses the untruncated tensor (Lee: "a truncation of the rank for the two-body
operator but not the one-body operator").

Invariance note: the supermatrix transforms by the orthogonal similarity O (x) O under a
real orbital rotation O, so the f^(l) spectra, the truncation criterion and the spectrum of
T' are all invariant.  This DF lambda therefore does NOT depend on the orbital rotation
within the active space (up to ties among degenerate supermatrix eigenvalues), unlike the
raw Jordan-Wigner Pauli 1-norm.

This implementation works in the packed symmetric pair space (p <= q, off-diagonal pairs
weighted by sqrt 2), an isometric image of the symmetric part of the supermatrix; it does
not call openfermion.  Validated against Lee et al. Table XIII (Reiher FeMoco) in
tests/test_df_lambda.py.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np


def one_body_tprime(h1: np.ndarray, eri: np.ndarray) -> np.ndarray:
    """T'_pq = h_pq - 1/2 sum_r (pr|rq) + sum_r (pq|rr)  (Lee et al. Eqs. (3), (A9))."""
    return h1 - 0.5 * np.einsum('prrq->pq', eri) + np.einsum('pqrr->pq', eri)


def _packed_supermatrix(eri: np.ndarray) -> tuple:
    """Isometric packing of the symmetric pair space: returns (M_packed, iu, w)."""
    n = eri.shape[0]
    iu = np.triu_indices(n)
    w = np.where(iu[0] == iu[1], 1.0, np.sqrt(2.0))
    M = eri[iu[0], iu[1]][:, iu[0], iu[1]] * np.outer(w, w)
    return M, iu, w


def df_leaves(eri: np.ndarray, rank_tol: float = 0.0) -> List[np.ndarray]:
    """First factorization (Eq. (C1)): list of symmetric N x N leaves W^(l) with
    (pq|rs) = sum_l W^(l)_pq W^(l)_rs, in descending supermatrix-eigenvalue order.
    Eigenvalues <= rank_tol are dropped (rank_tol = 0 keeps every positive one)."""
    n = eri.shape[0]
    M, iu, w = _packed_supermatrix(eri)
    evals, evecs = np.linalg.eigh(0.5 * (M + M.T))
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    keep = evals > rank_tol
    leaves = []
    for lam, v in zip(evals[keep], evecs[:, keep].T):
        W = np.zeros((n, n))
        W[iu] = v / w
        W = W + np.triu(W, 1).T
        leaves.append(np.sqrt(lam) * W)
    return leaves


def df_lambda(h1: np.ndarray, eri: np.ndarray, ecore: float = 0.0, threshold: float = 1e-3,
              rank_tol: float = 0.0, return_eri: bool = False,
              reconstruct: bool = True) -> Dict[str, object]:
    """Double-factorized 1-norm lambda = lambda_T' + lambda_DF (Lee et al. Eqs. (A9), (C11),
    truncation Eq. (C41)) of H given by h1 (N x N), eri (chemist (pq|rs), N^4), ecore.

    threshold: the Eq. (C41) second-factorization cut (Hartree^2); 0 keeps everything.
    rank_tol: first-factorization eigenvalue cut (0 = every positive eigenvalue, as in
    openfermion.resource_estimates.df.factorize).
    reconstruct: build the truncated tensor (for max_abs_eri_error / return_eri); False
    skips the O(L N^4) rebuild and reports max_abs_eri_error as nan.

    Returns lambda_total, lambda_one_body (lambda_T'), lambda_two_body (lambda_DF),
    n_leaves (L, leaves with at least one kept eigenvector), rank (first-factorization
    rank before truncation), n_eigenvectors (Xi = total kept eigenvectors), max_abs_eri_error
    (of the truncated reconstruction), ecore (passed through; not part of lambda), and
    eri_df (the reconstructed tensor) if return_eri.
    """
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    n = h1.shape[0]
    if eri.shape != (n, n, n, n):
        raise ValueError('eri must be N^4 with N = h1.shape[0]')
    lam_one = float(np.abs(np.linalg.eigvalsh(one_body_tprime(h1, eri))).sum())
    leaves = df_leaves(eri, rank_tol)
    lam_two = 0.0
    n_leaves = 0
    n_vec = 0
    rec = np.zeros_like(eri)
    for W in leaves:
        f, U = np.linalg.eigh(W)
        s = np.abs(f).sum()
        kept = s * np.abs(f) > threshold
        if not kept.any():
            break
        n_leaves += 1
        n_vec += int(kept.sum())
        lam_two += 0.25 * np.abs(f[kept]).sum() ** 2
        Wk = (U[:, kept] * f[kept]) @ U[:, kept].T
        if reconstruct or return_eri:
            rec += np.einsum('pq,rs->pqrs', Wk, Wk)
    out: Dict[str, object] = {
        'lambda_total': float(lam_one + lam_two), 'lambda_one_body': lam_one,
        'lambda_two_body': float(lam_two), 'n_leaves': n_leaves, 'rank': len(leaves),
        'n_eigenvectors': n_vec, 'max_abs_eri_error': float(np.abs(rec - eri).max()) if (reconstruct or return_eri) else float('nan'),
        'threshold': threshold, 'ecore': float(ecore)}
    if return_eri:
        out['eri_df'] = rec
    return out
