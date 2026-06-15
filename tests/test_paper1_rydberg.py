"""
Paper 1 §II — eq:rydberg dynamical-algebra backing test.

Paper 1 §II claims the SU(2) (x) SU(1,1) dynamical-algebra operators and
their Casimir invariants reproduce the exact hydrogen spectrum
E_n = -1/(2 n^2), "from operator eigenvalues, not from solving differential
equations." The 2026-06-14 trunk QA re-check flagged this as a load-bearing
claim with no GeoVac artifact (the only backing was the continuum Fock chain;
the production graph Laplacian L=D-A is positive-semidefinite).

This test builds T+- (radial SU(1,1)) and L+- (angular SU(2)) from the
Biedenharn-Louck matrix elements Paper 1 documents in its appendix:

    <l,m+1|L+|l,m>  = sqrt((l-m)(l+m+1))
    <n+1,l|T+|n,l>  = sqrt((n-l)(n+l+1)/4)

and forms, PURELY from the ladder operators (nothing about n or l is placed
on a diagonal by hand):

    N   := -2 [T+, T-]                       (SU(1,1) number operator)
    L^2 := Lz^2 + 1/2 (L+ L- + L- L+)        (SU(2) Casimir)

The non-tautological content is that

    * N has integer spectrum exactly {n}     (verified by hand: n=1 -> 1,
      n=2 -> 2; the -1/2 commutator eigenvalue times -2 gives n), and
    * L^2 has spectrum exactly {l(l+1)},

so feeding N through the energy-shell relation E = -1/(2 N^2) yields the
exact Rydberg levels. n and l(l+1) emerge as eigenvalues of operators
assembled from the generators, not inserted by hand.

Truncation note: on the top shell n = max_n the raise partner is outside the
finite window, so the commutator is truncated there; the N assertions are
made only on INTERIOR states (those whose |n+1,l,m> partner is present).
L^2 is exact on every state (the m-tower at fixed (n,l) is complete).
"""
from __future__ import annotations

import math

import numpy as np


def build_operators(max_n: int):
    states = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    idx = {s: i for i, s in enumerate(states)}
    dim = len(states)
    Lp = np.zeros((dim, dim))
    Tp = np.zeros((dim, dim))
    Lz = np.zeros((dim, dim))
    for (n, l, m) in states:
        i = idx[(n, l, m)]
        Lz[i, i] = m
        if m < l:  # L+ : <l,m+1|L+|l,m> = sqrt((l-m)(l+m+1))
            j = idx[(n, l, m + 1)]
            Lp[j, i] = math.sqrt((l - m) * (l + m + 1))
        if (n + 1, l, m) in idx:  # T+ : <n+1,l|T+|n,l> = sqrt((n-l)(n+l+1)/4)
            j = idx[(n + 1, l, m)]
            Tp[j, i] = math.sqrt((n - l) * (n + l + 1) / 4.0)
    return states, idx, Lp, Tp, Lz


def test_su11_number_operator_integer_spectrum_n():
    """N := -2[T+,T-] is diagonal with eigenvalue exactly n on interior
    states -- n emerges from the radial ladder algebra, not by hand."""
    max_n = 6
    states, idx, Lp, Tp, Lz = build_operators(max_n)
    Tm = Tp.T
    N = -2.0 * (Tp @ Tm - Tm @ Tp)
    off = N - np.diag(np.diag(N))
    assert np.linalg.norm(off, 2) < 1e-10, "N is not diagonal"
    checked = 0
    for (n, l, m) in states:
        if (n + 1, l, m) in idx:  # interior: commutator exact
            i = idx[(n, l, m)]
            assert abs(N[i, i] - n) < 1e-9, (n, l, m, N[i, i])
            checked += 1
    assert checked > 0


def test_su2_casimir_spectrum_l_l_plus_1():
    """L^2 := Lz^2 + 1/2(L+L- + L-L+) has eigenvalue exactly l(l+1)."""
    max_n = 6
    states, idx, Lp, Tp, Lz = build_operators(max_n)
    Lm = Lp.T
    L2 = Lz @ Lz + 0.5 * (Lp @ Lm + Lm @ Lp)
    off = L2 - np.diag(np.diag(L2))
    assert np.linalg.norm(off, 2) < 1e-10, "L^2 is not diagonal"
    for (n, l, m) in states:
        i = idx[(n, l, m)]
        assert abs(L2[i, i] - l * (l + 1)) < 1e-9, (n, l, m, L2[i, i])


def test_rydberg_levels_from_algebra():
    """E = -1/(2 N^2) with N the algebra's number operator reproduces
    E_n = -1/(2 n^2) exactly on interior states (the §II claim chain)."""
    max_n = 6
    states, idx, Lp, Tp, Lz = build_operators(max_n)
    Tm = Tp.T
    N = -2.0 * (Tp @ Tm - Tm @ Tp)
    for (n, l, m) in states:
        if (n + 1, l, m) in idx:
            i = idx[(n, l, m)]
            E = -1.0 / (2.0 * N[i, i] ** 2)
            assert abs(E - (-1.0 / (2.0 * n**2))) < 1e-12, (n, l, m, E)
