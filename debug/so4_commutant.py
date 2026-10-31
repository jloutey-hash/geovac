"""SO(4) generators in the hydrogenic (n, l, m) basis used by geovac.shibuya_wulfman, and the commutant
of an operator inside so(4).  Record SY (docs/predictions/2026-10-30_so4_commutant_door11.md).

    python debug/so4_commutant.py controls     algebra checks + the two controls (invariant, central)
    python debug/so4_commutant.py score        door 11: the second nucleus at R = 2, 4, 1 (run after the record is committed)

Generators: L from the standard ladder algebra (Condon-Shortley, complex Y_lm, as the Gaunt coefficients of
the cross-centre builder); A from Pauli's shell identity P_n r P_n = -(3n/2) P_n A P_n (Z = 1), with the
dipole matrix built from the builder's own radial polynomials and the same Gaunt coefficients.
"""
import sys
from math import factorial, sqrt

import numpy as np
from sympy import N as sN
from sympy.physics.wigner import gaunt

import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.shibuya_wulfman import _hydrogenic_poly_coeffs, _angular_coefficient  # noqa: E402

NMAX = 3
STATES = [(n, l, m) for n in range(1, NMAX + 1) for l in range(n) for m in range(-l, l + 1)]
IDX = {s: i for i, s in enumerate(STATES)}
DIM = len(STATES)


def radial_r3(n1, l1, n2, l2):
    """int R_{n1 l1} R_{n2 l2} r^3 dr with the builder's polynomials."""
    c1, a1 = _hydrogenic_poly_coeffs(1.0, n1, l1)
    c2, a2 = _hydrogenic_poly_coeffs(1.0, n2, l2)
    tot = 0.0
    for i, ci in enumerate(c1):
        for j, cj in enumerate(c2):
            p = i + j + 3
            tot += ci * cj * factorial(p) / (a1 + a2) ** (p + 1)
    return tot


def ang(l1, m1, l2, m2, mu):
    """<l1 m1| sqrt(4 pi / 3) Y_{1 mu} |l2 m2> in the complex basis."""
    return float(sN((-1) ** m1 * gaunt(l1, 1, l2, -m1, mu, m2))) * sqrt(4 * np.pi / 3)


def dipole_block(n):
    """x, y, z restricted to shell n, in the complex basis."""
    sh = [s for s in STATES if s[0] == n]
    k = len(sh)
    X = np.zeros((k, k), complex); Y = np.zeros((k, k), complex); Z = np.zeros((k, k), complex)
    for i, (_, l1, m1) in enumerate(sh):
        for j, (_, l2, m2) in enumerate(sh):
            rad = radial_r3(n, l1, n, l2)
            z = ang(l1, m1, l2, m2, 0)
            yp = ang(l1, m1, l2, m2, 1)
            ym = ang(l1, m1, l2, m2, -1)
            # x = r sqrt(2pi/3) (Y_{1,-1} - Y_{1,1}),  y = i r sqrt(2pi/3) (Y_{1,-1} + Y_{1,1})
            Z[i, j] = rad * z
            X[i, j] = rad * (ym - yp) / sqrt(2)
            Y[i, j] = rad * 1j * (ym + yp) / sqrt(2)
    return sh, X, Y, Z


def generators():
    L = [np.zeros((DIM, DIM), complex) for _ in range(3)]
    A = [np.zeros((DIM, DIM), complex) for _ in range(3)]
    for (n, l, m) in STATES:
        i = IDX[(n, l, m)]
        L[2][i, i] = m
        if m + 1 <= l:
            L[0][IDX[(n, l, m + 1)], i] += sqrt(l * (l + 1) - m * (m + 1)) / 2
            L[1][IDX[(n, l, m + 1)], i] += sqrt(l * (l + 1) - m * (m + 1)) / (2j)
        if m - 1 >= -l:
            L[0][IDX[(n, l, m - 1)], i] += sqrt(l * (l + 1) - m * (m - 1)) / 2
            L[1][IDX[(n, l, m - 1)], i] += -sqrt(l * (l + 1) - m * (m - 1)) / (2j)
    for n in range(1, NMAX + 1):
        sh, X, Y, Z = dipole_block(n)
        ids = [IDX[s] for s in sh]
        for k, M in enumerate((X, Y, Z)):
            A[k][np.ix_(ids, ids)] = -(2.0 / (3.0 * n)) * M
    return L, A


def comm(P, Q):
    return P @ Q - Q @ P


def commutant(V, gens, names, zero_tol):
    """Null space of X -> [V, X] over span(gens); returns (dimension, singular values, basis vectors)."""
    M = np.array([comm(V, G).ravel() for G in gens]).T   # DIM^2 x 6
    u, s, vh = np.linalg.svd(M, full_matrices=True)
    scale = np.linalg.norm(V) * max(np.linalg.norm(G) for G in gens)
    rel = s / scale
    null = [vh[k].conj() for k in range(len(gens)) if (rel[k] if k < len(rel) else 0.0) < zero_tol]
    return len(null), rel, null


if __name__ == "__main__" and (len(sys.argv) < 2 or sys.argv[1] == "controls"):
    L, A = generators()
    names = ["Lx", "Ly", "Lz", "Ax", "Ay", "Az"]
    gens = L + A
    # --- algebra checks (conventions) ---
    eps = lambda i, j, k: np.sign((j - i) * (k - i) * (k - j))  # noqa: E731
    err_LL = max(np.abs(comm(L[i], L[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max()
                 for i in range(3) for j in range(3))
    err_LA = max(np.abs(comm(L[i], A[j]) - 1j * sum(eps(i, j, k) * A[k] for k in range(3))).max()
                 for i in range(3) for j in range(3))
    err_AA = max(np.abs(comm(A[i], A[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max()
                 for i in range(3) for j in range(3))
    cas = sum(G @ G for G in gens)
    cas_expect = np.diag([n * n - 1.0 for (n, l, m) in STATES])
    err_cas = np.abs(cas - cas_expect).max()
    # the builder's own z-angular factor against ours
    err_ang = max(abs(_angular_coefficient(l1, m, l2, m, 1) - ang(l1, m, l2, m, 0))
                  for (_, l1, m) in STATES for (_, l2, m2) in STATES if m2 == m)
    print("[L,L]=iL err", err_LL, " [L,A]=iA err", err_LA, " [A,A]=iL err", err_AA)
    print("L^2+A^2 = n^2-1 err", err_cas, "  builder z-angular vs ours", err_ang)

    # --- control 1: an SO(4)-invariant operator (a function of H): commutant must be all six ---
    H = np.diag([-0.5 / (n * n) for (n, l, m) in STATES]).astype(complex)
    d, rel, _ = commutant(H, gens, names, 1e-10)
    print("control invariant: dim", d, " rel singular values", np.array2string(rel, precision=2))

    # --- control 2: a central operator (finite nucleus, uniform sphere R = 0.05): commutant must be L only ---
    R = 0.05
    r = np.linspace(1e-6, R, 20001)
    dV = -(3 / (2 * R) - r * r / (2 * R ** 3)) + 1 / r
    Vc = np.zeros((DIM, DIM), complex)
    polys = {(n, l): _hydrogenic_poly_coeffs(1.0, n, l) for n in range(1, NMAX + 1) for l in range(n)}
    for (n1, l1, m1) in STATES:
        for (n2, l2, m2) in STATES:
            if (l1, m1) != (l2, m2):
                continue
            c1, a1 = polys[(n1, l1)]; c2, a2 = polys[(n2, l2)]
            R1 = np.exp(-a1 * r) * np.polyval(c1[::-1], r)
            R2 = np.exp(-a2 * r) * np.polyval(c2[::-1], r)
            Vc[IDX[(n1, l1, m1)], IDX[(n2, l2, m2)]] = np.trapezoid(R1 * R2 * dV * r * r, r)
    d, rel, null = commutant(Vc, gens, names, 1e-10)
    print("control central: dim", d, " rel singular values", np.array2string(rel, precision=3))
    for v in null:
        print("   null vector", np.array2string(np.abs(v), precision=3))
    print("control central: max rel commutator among L", max(np.linalg.norm(comm(Vc, G)) for G in L) / np.linalg.norm(Vc),
          " min among A", min(np.linalg.norm(comm(Vc, G)) for G in A) / np.linalg.norm(Vc))


def score():
    from geovac.shibuya_wulfman import compute_cross_center_vne
    L, A = generators()
    gens = L + A
    names = ["Lx", "Ly", "Lz", "Ax", "Ay", "Az"]
    hopf = L[2] + A[2]
    for R_AB in (2.0, 4.0, 1.0):
        V = compute_cross_center_vne(1.0, STATES, 1.0, R_AB, 6).astype(complex)
        d, rel, null = commutant(V, gens, names, 1e-10)
        print("door 11, R = %.1f: commutant dim %d; rel singular values %s" % (R_AB, d, np.array2string(rel, precision=3)))
        for v in null:
            print("   null vector |c| over", names, np.array2string(np.abs(v), precision=4))
        nV = np.linalg.norm(V)
        print("   rel commutator: Lz %.2e  Lz+Az (Hopf) %.3f  Az %.3f  Lx %.3f" % (
            np.linalg.norm(comm(V, L[2])) / nV, np.linalg.norm(comm(V, hopf)) / nV,
            np.linalg.norm(comm(V, A[2])) / nV, np.linalg.norm(comm(V, L[0])) / nV))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "score":
    score()
