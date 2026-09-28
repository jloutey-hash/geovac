"""Independent build of the TRUE Camporesi-Higuchi truncation on S^3 = SU(2).

Second-route check of the Paper 38 substrate finding A2 (2026-09-28 /qa trunk).
Written from scratch; does not import any geovac module or reviewer script.

Realization (Peter-Weyl + left-invariant framing):
  spinor sections = L^2(SU(2)) (x) C^2,
  L^2(SU(2)) = (+)_j V_j(m) (x) V_j(m'), orthonormal basis e^j_{m m'} = sqrt(2j+1) D^j_{m m'}.
  Left translations act on m, right translations on m'.  Left-invariant vector
  fields generate right translations, so J_R = standard spin-j matrices on m'.
  Dirac operator (unit round S^3, left-invariant trivialization):
        D = 2 sigma . J_R + 3/2      (identity on m).
  On V_j(m') (x) C^2 with K = J_R + sigma/2:  K=j+1/2 -> +(2j+3/2), K=j-1/2 -> -(2j+1/2),
  i.e. spectrum +-(n+3/2), n>=0, multiplicity (n+1)(n+2) (verified numerically below).
  Truncation: H_N = span of eigenspaces with n = 0..n_max-1 (Fock n = 1..n_max),
  i.e. |lambda| <= n_max + 1/2.

Multiplication by a band-N function (spin J=(N-1)/2 matrix coefficient, N^2 of them):
  <e^{j2}_{m2 m2'}| D^J_{M M'} e^{j1}_{m1 m1'}> =
      sqrt((2j1+1)/(2j2+1)) <j1 m1 J M|j2 m2> <j1 m1' J M'|j2 m2'>,
  tensored with identity on C^2.  O = P C(S^3) P (band-limited part).
"""
from __future__ import annotations

import itertools
from functools import lru_cache
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan


@lru_cache(maxsize=None)
def cg(j1: Fraction, m1: Fraction, j2: Fraction, m2: Fraction, J: Fraction, M: Fraction) -> float:
    if m1 + m2 != M:
        return 0.0
    R = lambda x: Rational(x.numerator, x.denominator)
    return float(clebsch_gordan(R(j1), R(j2), R(J), R(m1), R(m2), R(M)))


def spin_mats(j: Fraction) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ms = [j - k for k in range(int(2 * j) + 1)]
    d = len(ms)
    jz = np.diag([float(m) for m in ms]).astype(complex)
    jp = np.zeros((d, d), complex)
    for a, m in enumerate(ms):
        if a > 0:  # |m> -> |m+1> at index a-1
            jp[a - 1, a] = np.sqrt(float(j * (j + 1) - m * (m + 1)))
    jx = (jp + jp.conj().T) / 2
    jy = (jp - jp.conj().T) / (2j)
    return jx, jy, jz


SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def build(n_max: int):
    js = [Fraction(k, 2) for k in range(0, n_max + 1)]  # j = 0 .. n_max/2
    labels: List[Tuple[Fraction, Fraction, Fraction, int]] = []
    for j in js:
        ms = [j - k for k in range(int(2 * j) + 1)]
        for m, mp, s in itertools.product(ms, ms, (0, 1)):
            labels.append((j, m, mp, s))
    idx = {lab: i for i, lab in enumerate(labels)}
    dim = len(labels)
    # Dirac operator, block diagonal in j
    D = np.zeros((dim, dim), complex)
    for j in js:
        d = int(2 * j) + 1
        jx, jy, jz = spin_mats(j)
        blk_r = 2 * (np.kron(jx, SX) + np.kron(jy, SY) + np.kron(jz, SZ)) + 1.5 * np.eye(2 * d)
        blk = np.kron(np.eye(d), blk_r)  # identity on m; ordering (m, m', s)
        ms = [j - k for k in range(d)]
        ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
        D[np.ix_(ids, ids)] = blk
    ev, U = np.linalg.eigh(D)
    return labels, idx, D, ev, U


def mult_op(J: Fraction, M: Fraction, Mp: Fraction, labels, idx) -> np.ndarray:
    dim = len(labels)
    A = np.zeros((dim, dim), complex)
    for (j1, m1, m1p, s) in labels:
        i1 = idx[(j1, m1, m1p, s)]
        for (j2, m2, m2p, s2) in labels:
            if s2 != s or m2 != m1 + M or m2p != m1p + Mp:
                continue
            if not (abs(j1 - J) <= j2 <= j1 + J):
                continue
            c = cg(j1, m1, J, M, j2, m2) * cg(j1, m1p, J, Mp, j2, m2p)
            if c != 0.0:
                A[idx[(j2, m2, m2p, s2)], i1] = np.sqrt(float((2 * j1 + 1) / (2 * j2 + 1))) * c
    return A


def main() -> None:
    for n_max in (1, 2, 3):
        labels, idx, D, ev, U = build(n_max)
        # spectrum check
        spec: Dict[float, int] = {}
        for e in np.round(ev, 9):
            spec[float(e)] = spec.get(float(e), 0) + 1
        print(f"\n=== n_max = {n_max} ===  ambient dim {len(labels)}")
        print("  spectrum (eigenvalue: multiplicity):", {k: v for k, v in sorted(spec.items())})
        keep = np.abs(ev) <= n_max + 0.5 + 1e-9
        V = U[:, keep]
        dimH = V.shape[1]
        expect = 2 * sum((n + 1) * (n + 2) for n in range(n_max))
        print(f"  dim H_N = {dimH} (expected 2*sum (n+1)(n+2) = {expect})")
        Dc = V.conj().T @ D @ V
        # Dirac check: [D, a] must equal geometric derivative -> verify ||[D,a]|| for a=D^{1/2} generic
        O_vecs: List[np.ndarray] = []
        comp_ops: List[np.ndarray] = []
        print("  band N (spin (N-1)/2): rank of compression / N^2")
        for N in range(1, 2 * n_max + 3):
            J = Fraction(N - 1, 2)
            Ms = [J - k for k in range(N)]
            vecs = []
            for M, Mp in itertools.product(Ms, Ms):
                A = mult_op(J, M, Mp, labels, idx)
                T = V.conj().T @ A @ V
                vecs.append(T.reshape(-1))
                comp_ops.append(T)
            Mtx = np.array(vecs)
            r = np.linalg.matrix_rank(Mtx, tol=1e-9)
            print(f"    N={N}: rank {r} / {N*N}")
            O_vecs.extend(vecs)
        Omat = np.array(O_vecs)
        dimO = np.linalg.matrix_rank(Omat, tol=1e-9)
        # orthonormal basis of O
        _, s, Vh = np.linalg.svd(Omat, full_matrices=False)
        basis = Vh[: dimO]
        comm = []
        for b in basis:
            T = b.reshape(dimH, dimH)
            comm.append((Dc @ T - T @ Dc).reshape(-1))
        rc = np.linalg.matrix_rank(np.array(comm), tol=1e-9)
        ker = dimO - rc
        # is the kernel exactly the scalars?  identity in O and commutes with D
        I = np.eye(dimH).reshape(-1)
        in_O = np.linalg.matrix_rank(np.vstack([basis, I]), tol=1e-9) == dimO
        # sanity: multiplication ops self-adjointness structure / Hermitian closure of O
        print(f"  dim O = {dimO};  rank of a->[D,a] on O = {rc};  kernel dim = {ker};  1 in O: {in_O}")
        # smallest nonzero singular value of the commutator map (quantifies L_D is a norm on O/C1)
        sv = np.linalg.svd(np.array(comm), compute_uv=False)
        print(f"  commutator-map singular values: min nonzero = {sv[sv>1e-9].min():.6f}, count zero = {int((sv<=1e-9).sum()) + (dimO-len(sv))}")
        # full-algebra check: D commutes with left translations only; confirm D not proportional to identity
        verdict = "L_D is a NORM on O/C.1" if ker == 1 and in_O else "L_D NOT a norm on O/C.1"
        print(f"  -> {verdict}")


if __name__ == "__main__":
    main()
