r"""Non-symmetric (bra != ket) ordered-pair coefficient tensors for the LiH pi marriage.

Backs (Paper 19): the CI-relaxed explicit-r12 energy (ci_relaxed.py) -- every new matrix element is
a contraction with a TRANSITION tensor, the generalisation of pi_m16.build_ordered_tensor_compact
from one coefficient matrix (bra = ket = CS) to two (bra X, ket Y):

    <X|O|Y>/4 = sum_{ABCD} T[A,B,C,D] <rho_A(1) rho_B(2) rho_C(3) rho_D(4) O>,
    T = (1/4) sum_{a b a' b'} X[a,b] Y[a',b'] E[(a,a'),(A,B)] E[(b,b'),(C,D)],

ordered pair index A = x*Mk + y with x the BRA orbital and y the KET orbital (rho_A = phi_x^* phi_y).
dA[a]^* dA[a'] (a = (p,q) bra, a' = (p2,q2) ket) expands to the four (left, right) ordered-pair
products  (pp2, qq2) - (pq2, qp2) - (qp2, pq2) + (qq2, pp2)  -- exactly the dp_keys of pi_m16.
Nothing is symmetrised.  Light module (numpy + ref_tensor only): no grid / kernel imports.
"""
from __future__ import annotations

from itertools import combinations
from typing import List, Sequence, Tuple

import numpy as np

import geovac.lih_r12ci.ref_tensor as P0

_SGN = np.array([1.0, -1.0, -1.0, 1.0])


def determinants(ct: np.ndarray, Mk: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    """Kept determinants of a (truncated) CI vector ct in the _dets ordering: alpha-pair index a_I,
    beta-pair index b_I, grouped-order coefficient d_I = CS[a_I, b_I] (psi = sum_I d_I chi_I with
    chi_I = dA[a_I] dB[b_I] / 2 unsigned), and the same-spin pair list."""
    CS, _, _ = P0.fold_cs(ct, Mk)
    ia, ib = np.nonzero(CS)
    apairs = list(combinations(range(Mk), 2))
    return ia, ib, CS[ia, ib].copy(), apairs


def _pair_terms(apairs, Mk: int, a, ap) -> Tuple[np.ndarray, np.ndarray]:
    """(left-electron, right-electron) ordered-pair keys of dA[a]^* dA[a'] (a bra, a' ket); (n, 4)."""
    a = np.asarray(a); ap = np.asarray(ap)
    pr = np.asarray(apairs)
    p = pr[a, 0][:, None]; q = pr[a, 1][:, None]
    p2 = pr[ap, 0][:, None]; q2 = pr[ap, 1][:, None]
    L = np.concatenate([p * Mk + p2, p * Mk + q2, q * Mk + p2, q * Mk + q2], axis=1)
    R = np.concatenate([q * Mk + q2, q * Mk + p2, p * Mk + q2, p * Mk + p2], axis=1)
    return L, R


def transition_stack(bra: Sequence[Tuple[int, int, float]], ket: Sequence[Tuple[int, int, float]],
                     owners: Sequence[int], apairs, Mk: int):
    """Stacked COO of transition tensors: contribution k = (1/4) bra[k].coef ket[k].coef E E, belonging
    to owner owners[k]; bra/ket entries are (alpha-pair index, beta-pair index, coefficient).
    Returns (iA, iB, iC, iD, val, owner) in full ordered-pair indices, duplicates merged per owner."""
    ba = np.array([b[0] for b in bra]); bb = np.array([b[1] for b in bra])
    bc = np.array([b[2] for b in bra], dtype=float)
    ka = np.array([k[0] for k in ket]); kb = np.array([k[1] for k in ket])
    kc = np.array([k[2] for k in ket], dtype=float)
    own = np.asarray(owners, dtype=np.int64)
    LA, RA = _pair_terms(apairs, Mk, ba, ka)                    # alpha: electrons 1, 2
    LB, RB = _pair_terms(apairs, Mk, bb, kb)                    # beta : electrons 3, 4
    coef = 0.25 * bc * kc
    iA = np.repeat(LA, 4, axis=1).reshape(-1); iB = np.repeat(RA, 4, axis=1).reshape(-1)
    iC = np.tile(LB, (1, 4)).reshape(-1); iD = np.tile(RB, (1, 4)).reshape(-1)
    val = (coef[:, None] * np.outer(_SGN, _SGN).reshape(-1)[None, :]).reshape(-1)
    ow = np.repeat(own, 16)
    npair = Mk * Mk
    key = (((ow * npair + iA) * npair + iB) * npair + iC) * npair + iD
    uk, inv = np.unique(key, return_inverse=True)
    v = np.bincount(inv.reshape(-1), weights=val)
    keep = v != 0.0
    uk = uk[keep]; v = v[keep]
    iD = uk % npair; uk = uk // npair
    iC = uk % npair; uk = uk // npair
    iB = uk % npair; uk = uk // npair
    iA = uk % npair; ow = uk // npair
    return iA, iB, iC, iD, v, ow


def split_owners(stack, n_owner: int) -> List[tuple]:
    """The stacked COO as one (iA, iB, iC, iD, val) tuple per owner."""
    iA, iB, iC, iD, val, ow = stack
    order = np.argsort(ow, kind='stable')
    bounds = np.searchsorted(ow[order], np.arange(n_owner + 1))
    return [(iA[s], iB[s], iC[s], iD[s], val[s]) for s in
            (order[bounds[o]:bounds[o + 1]] for o in range(n_owner))]


def bra_det_ket_psi(ia, ib, d) -> Tuple[list, list, list]:
    """(bra, ket, owners) for T_I = transition(bra chi_I, ket psi), owner I."""
    N = len(d)
    bra = [(ia[I], ib[I], 1.0) for I in range(N) for J in range(N)]
    ket = [(ia[J], ib[J], d[J]) for I in range(N) for J in range(N)]
    own = [I for I in range(N) for J in range(N)]
    return bra, ket, own


def bra_det_ket_det(ia, ib) -> Tuple[list, list, list]:
    """(bra, ket, owners) for T_IJ = transition(bra chi_I, ket chi_J), owner I*N+J."""
    N = len(ia)
    bra = [(ia[I], ib[I], 1.0) for I in range(N) for J in range(N)]
    ket = [(ia[J], ib[J], 1.0) for I in range(N) for J in range(N)]
    own = [I * N + J for I in range(N) for J in range(N)]
    return bra, ket, own
