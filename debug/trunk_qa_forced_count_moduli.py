"""Trunk QA — Claim B: genuine parameter-count of the Forced-Count moduli chain.

Computes dim_R M(D_F) at each reduction step of Paper 32 thm:forced_count
(claimed chain 1024 -> 512 -> 128 -> 128) by counting the real dimension of the
solution space of the linear constraints, on the SAME H_F = C^32 / gamma_F / J_F
that geovac.standard_model_triple.StandardModelFiniteTriple uses.

The constraints are real-linear operators on the 2048-real-dim space of 32x32
complex matrices.  Moduli dimension = 2048 - rank(stacked constraint operator).

Order-one [[D,a], J b J^{-1}] = 0 is imposed over a real basis of A_F.
"""
from __future__ import annotations

import numpy as np

from geovac.standard_model_triple import StandardModelFiniteTriple

f = StandardModelFiniteTriple()
GAMMA = np.real(f.chirality_F())          # 32x32 real diag +-1
J_U = f.real_structure_F()                # unitary part U of J = U K
N = 32
RDIM = 2 * N * N                          # 2048


def from_realvec(v):
    re = v[: N * N].reshape(N, N)
    im = v[N * N :].reshape(N, N)
    return re + 1j * im


def realvec(M):
    return np.concatenate([np.real(M).ravel(), np.imag(M).ravel()])


# --- real basis for the input space (columns = standard real coords) ---
_BASIS_INPUTS = [from_realvec(np.eye(RDIM)[:, k]) for k in range(RDIM)]


def constraint_matrix(func):
    """func: complex 32x32 -> complex violation; returns its real matrix (rows=viol)."""
    cols = [realvec(func(M)) for M in _BASIS_INPUTS]
    return np.array(cols).T


def herm(M):
    return M - M.conj().T


def chir(M):
    return GAMMA @ M + M @ GAMMA


def jreal(M):
    return J_U @ np.conj(M) @ J_U.conj().T - M


# --- A_F real basis for order-one ---
def piF(lam, q, m):
    return f.algebra_action(lam, q, m)


def a_f_basis():
    out = []
    for lam in (1.0, 1j):
        out.append(piF(lam, (0, 0, 0, 0), np.zeros((3, 3), complex)))
    for k in range(4):
        q = [0, 0, 0, 0]
        q[k] = 1
        out.append(piF(0, tuple(q), np.zeros((3, 3), complex)))
    for i in range(3):
        for j in range(3):
            for val in (1.0, 1j):
                m = np.zeros((3, 3), complex)
                m[i, j] = val
                out.append(piF(0, (0, 0, 0, 0), m))
    return out


A_BASIS = a_f_basis()


def order_one_matrix():
    """Stack [[D,a],JbJ^-1] over all (a,b) in A_BASIS as a real-linear op in D."""
    JbJ = [J_U @ np.conj(b) @ J_U.conj().T for b in A_BASIS]
    cols = []
    for M in _BASIS_INPUTS:
        parts = []
        for a in A_BASIS:
            Da = M @ a - a @ M
            for Jb in JbJ:
                parts.append(realvec(Da @ Jb - Jb @ Da))
        cols.append(np.concatenate(parts))
    return np.array(cols).T


def nulldim(*mats):
    # rank(A) == rank(A^T A); use the small RDIM x RDIM Gram matrix for speed
    # (A here is rows=violation coords, cols=input coords; we want column rank).
    G = np.zeros((RDIM, RDIM))
    for A in mats:
        G += A.T @ A
    r = np.linalg.matrix_rank(G, tol=1e-7)
    return RDIM - r


if __name__ == "__main__":
    H = constraint_matrix(herm)
    C = constraint_matrix(chir)
    Jc = constraint_matrix(jreal)
    print(f"general complex 32x32 real dim         = {RDIM}")
    print(f"Hermitian only                          = {nulldim(H)}")
    print(f"Herm + chir-anticomm                    = {nulldim(H, C)}")
    print(f"Herm + chir + J-real                    = {nulldim(H, C, Jc)}")
    print(f"A_F real basis size                     = {len(A_BASIS)}")
    O = order_one_matrix()
    print(f"order-one constraint matrix shape       = {O.shape}")
    print(f"Herm + chir + J-real + order-one        = {nulldim(H, C, Jc, O)}")
    print()
    print("Paper 32 thm:forced_count claims chain  : 1024 -> 512 -> 128 -> 128")
