"""Guard for the atlas door-11 cell decided by record SY (docs/predictions/2026-10-30_so4_commutant_door11.md):
a second nucleus on the z axis leaves exactly one direction of so(4), the rotation L_z about the axis, and
breaks the Hopf fibre's generator L_z + A_z.

Self-contained (no debug/ import): the six SO(4) generators are rebuilt here in the hydrogenic (n, l, m)
basis of geovac.shibuya_wulfman, L from the ladder algebra and A from Pauli's shell identity
P_n r P_n = -(3n/2) P_n A P_n (Z = 1), and the algebra is checked before use.

Wrong answers this guard rejects: "the Hopf U(1) survives a second nucleus" (its commutator would vanish);
"some Runge-Lenz direction survives" (the commutant would be larger than one); "the generators are wrong"
(the algebra check would fail); and a commutant routine that cannot see a symmetry (the central control
would not return the three rotations).
"""
from __future__ import annotations

from math import factorial, sqrt
from typing import List, Tuple

import numpy as np
import pytest
from sympy import N as sN
from sympy.physics.wigner import gaunt

from geovac.shibuya_wulfman import _hydrogenic_poly_coeffs, compute_cross_center_vne

NMAX = 3
STATES = [(n, l, m) for n in range(1, NMAX + 1) for l in range(n) for m in range(-l, l + 1)]
IDX = {s: i for i, s in enumerate(STATES)}
DIM = len(STATES)


def _radial_r3(n1: int, l1: int, n2: int, l2: int) -> float:
    c1, a1 = _hydrogenic_poly_coeffs(1.0, n1, l1)
    c2, a2 = _hydrogenic_poly_coeffs(1.0, n2, l2)
    tot = 0.0
    for i, ci in enumerate(c1):
        for j, cj in enumerate(c2):
            p = i + j + 3
            tot += ci * cj * factorial(p) / (a1 + a2) ** (p + 1)
    return tot


def _ang(l1: int, m1: int, l2: int, m2: int, mu: int) -> float:
    return float(sN((-1) ** m1 * gaunt(l1, 1, l2, -m1, mu, m2))) * sqrt(4 * np.pi / 3)


def generators() -> Tuple[List[np.ndarray], List[np.ndarray]]:
    L = [np.zeros((DIM, DIM), complex) for _ in range(3)]
    A = [np.zeros((DIM, DIM), complex) for _ in range(3)]
    for (n, l, m) in STATES:
        i = IDX[(n, l, m)]
        L[2][i, i] = m
        if m + 1 <= l:
            c = sqrt(l * (l + 1) - m * (m + 1))
            L[0][IDX[(n, l, m + 1)], i] += c / 2
            L[1][IDX[(n, l, m + 1)], i] += c / (2j)
        if m - 1 >= -l:
            c = sqrt(l * (l + 1) - m * (m - 1))
            L[0][IDX[(n, l, m - 1)], i] += c / 2
            L[1][IDX[(n, l, m - 1)], i] += -c / (2j)
    for n in range(1, NMAX + 1):
        sh = [s for s in STATES if s[0] == n]
        ids = [IDX[s] for s in sh]
        k = len(sh)
        X = np.zeros((k, k), complex)
        Y = np.zeros((k, k), complex)
        Z = np.zeros((k, k), complex)
        for i, (_, l1, m1) in enumerate(sh):
            for j, (_, l2, m2) in enumerate(sh):
                rad = _radial_r3(n, l1, n, l2)
                yp, ym = _ang(l1, m1, l2, m2, 1), _ang(l1, m1, l2, m2, -1)
                Z[i, j] = rad * _ang(l1, m1, l2, m2, 0)
                X[i, j] = rad * (ym - yp) / sqrt(2)
                Y[i, j] = rad * 1j * (ym + yp) / sqrt(2)
        for a, M in enumerate((X, Y, Z)):
            A[a][np.ix_(ids, ids)] = -(2.0 / (3.0 * n)) * M
    return L, A


def _comm(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    return P @ Q - Q @ P


def commutant(V: np.ndarray, gens: List[np.ndarray], zero_tol: float = 1e-10) -> Tuple[int, np.ndarray, List[np.ndarray]]:
    M = np.array([_comm(V, G).ravel() for G in gens]).T
    _, s, vh = np.linalg.svd(M, full_matrices=True)
    rel = s / (np.linalg.norm(V) * max(np.linalg.norm(G) for G in gens))
    null = [vh[k].conj() for k in range(len(gens)) if rel[k] < zero_tol]
    return len(null), rel, null


@pytest.fixture(scope="module")
def gens():
    L, A = generators()
    return L, A


def test_sy_generators_close_the_so4_algebra(gens):
    L, A = gens
    eps = lambda i, j, k: np.sign((j - i) * (k - i) * (k - j))  # noqa: E731
    for i in range(3):
        for j in range(3):
            assert np.abs(_comm(L[i], L[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max() < 1e-12
            assert np.abs(_comm(L[i], A[j]) - 1j * sum(eps(i, j, k) * A[k] for k in range(3))).max() < 1e-12
            assert np.abs(_comm(A[i], A[j]) - 1j * sum(eps(i, j, k) * L[k] for k in range(3))).max() < 1e-12
    cas = sum(G @ G for G in L + A)
    assert np.abs(cas - np.diag([n * n - 1.0 for (n, l, m) in STATES])).max() < 1e-12


def _central_correction(R: float = 0.05) -> np.ndarray:
    """Door 24's operator: the uniform-sphere finite-nucleus correction (radius R), Z = 1."""
    r = np.linspace(1e-6, R, 4001)
    dV = -(3 / (2 * R) - r * r / (2 * R ** 3)) + 1 / r
    Vc = np.zeros((DIM, DIM), complex)
    for (n1, l1, m1) in STATES:
        for (n2, l2, m2) in STATES:
            if (l1, m1) != (l2, m2):
                continue
            c1, a1 = _hydrogenic_poly_coeffs(1.0, n1, l1)
            c2, a2 = _hydrogenic_poly_coeffs(1.0, n2, l2)
            R1 = np.exp(-a1 * r) * np.polyval(c1[::-1], r)
            R2 = np.exp(-a2 * r) * np.polyval(c2[::-1], r)
            Vc[IDX[(n1, l1, m1)], IDX[(n2, l2, m2)]] = np.trapezoid(R1 * R2 * dV * r * r, r)
    return Vc


def test_door24_both_chiral_factors_are_broken(gens):
    """Atlas door 24 (decided 2026-10-31): a finite nucleus breaks both chiral SU(2) factors
    J^+- = (L +- A)/2, the left and right multiplications that make S^3 a group; only their sum L survives.

    Wrong answers rejected: "the chiral SU(2) of Paper 30's S^3 = SU(2) survives a finite nucleus" (a
    commutant inside span J^+ or span J^-); a commutant routine blind inside a 3-dim span (the SO(4)-invariant
    control must keep all three of each factor).
    """
    L, A = gens
    Jp = [(L[a] + A[a]) / 2 for a in range(3)]
    Jm = [(L[a] - A[a]) / 2 for a in range(3)]
    H = np.diag([-0.5 / n ** 2 for (n, l, m) in STATES]).astype(complex)
    for J in (Jp, Jm):
        assert commutant(H, J)[0] == 3
    Vc = _central_correction()
    for J in (Jp, Jm):
        d, rel, _ = commutant(Vc, J)
        assert d == 0, d
        assert min(rel) > 1e-3, rel
    assert commutant(Vc, L)[0] == 3


def test_door24_rotations_fix_a_vector_in_every_shell_and_the_chiral_factor_none(gens):
    """The label-level trace of fixed points against free action: on each shell n >= 2 the rotation group L
    leaves exactly one vector invariant (the s state) and the chiral factor J^+ leaves none.

    Wrong answers rejected: a survivor that acts freely (no invariant vector per shell), and swapping the roles
    of L and J^+.
    """
    L, A = gens
    Jp = [(L[a] + A[a]) / 2 for a in range(3)]
    for n in range(2, NMAX + 1):
        ids = [IDX[s] for s in STATES if s[0] == n]
        for G, expect in ((L, 1), (Jp, 0)):
            M = np.vstack([g[np.ix_(ids, ids)] for g in G])
            s = np.linalg.svd(M, compute_uv=False)
            inv = len(ids) - int(np.sum(s > 1e-10))
            assert inv == expect, (n, expect, inv)


def test_sy_central_control_keeps_exactly_the_rotations(gens):
    L, A = gens
    R = 0.05
    r = np.linspace(1e-6, R, 4001)
    dV = -(3 / (2 * R) - r * r / (2 * R ** 3)) + 1 / r
    Vc = np.zeros((DIM, DIM), complex)
    for (n1, l1, m1) in STATES:
        for (n2, l2, m2) in STATES:
            if (l1, m1) != (l2, m2):
                continue
            c1, a1 = _hydrogenic_poly_coeffs(1.0, n1, l1)
            c2, a2 = _hydrogenic_poly_coeffs(1.0, n2, l2)
            R1 = np.exp(-a1 * r) * np.polyval(c1[::-1], r)
            R2 = np.exp(-a2 * r) * np.polyval(c2[::-1], r)
            Vc[IDX[(n1, l1, m1)], IDX[(n2, l2, m2)]] = np.trapezoid(R1 * R2 * dV * r * r, r)
    d, rel, null = commutant(Vc, L + A)
    assert d == 3
    for v in null:
        assert np.abs(v[3:]).max() < 1e-6, "a central operator breaks every Runge-Lenz direction"
    assert sorted(rel)[3] > 1e-3


def test_sy1_second_nucleus_leaves_only_Lz(gens):
    L, A = gens
    V = compute_cross_center_vne(1.0, STATES, 1.0, 2.0, 6).astype(complex)
    d, rel, null = commutant(V, L + A)
    assert d == 1
    v = np.abs(null[0])
    assert v[2] > 0.999
    assert max(v[0], v[1], v[3], v[4], v[5]) < 0.03
    assert sorted(rel)[1] > 1e-3, "the five broken directions are broken by a wide margin"


def test_sy2_hopf_fibre_generator_is_broken(gens):
    L, A = gens
    V = compute_cross_center_vne(1.0, STATES, 1.0, 2.0, 6).astype(complex)
    hopf = L[2] + A[2]
    assert np.linalg.norm(_comm(V, L[2])) / np.linalg.norm(V) < 1e-12
    assert np.linalg.norm(_comm(V, hopf)) / np.linalg.norm(V) > 1e-3
