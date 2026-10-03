"""
Paper 1 §II — eq:rydberg dynamical-algebra backing test.

Paper 1 §II claims the SU(2) (x) SU(1,1) dynamical-algebra operators and
their Casimir invariants reproduce the hydrogen QUANTUM NUMBERS (n on
interior shells, l(l+1) on every shell) from operator eigenvalues.  The
energies E_n = -1/(2 n^2) are then an IMPOSED map from n through Fock's
energy shell p0^2 = -2E -- not an eigenvalue output of the algebra, and
not a reproduction of the spectrum by the operators. The 2026-06-14 trunk QA re-check flagged this as a load-bearing
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

so feeding N through the (imposed) energy-shell map E = -1/(2 N^2) gives
the Rydberg levels.  What is tested is that n and l(l+1) emerge as
eigenvalues of operators assembled from the generators, not inserted by
hand; the energy map itself is imposed, not derived.

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


def test_number_operator_exact_arithmetic():
    """[SYMBOLIC PROOF] backing (Paper 1 sec II, upgraded 2026-10-02): in EXACT
    arithmetic (sympy, the T+ weights as surds) N = -2[T+,T-] equals n on every
    interior state and -(n-1-l)(n+l)/2 on the truncated top shell.  The floating
    tests above only bound |N - n| < 1e-9; this one rejects any weight convention
    whose commutator is not exactly n (e.g. w^2 = (n-l)(n+l)/4 gives n - 1/2)."""
    import sympy as sp

    max_n = 6
    states = [(n, l, m) for n in range(1, max_n + 1) for l in range(n)
              for m in range(-l, l + 1)]
    idx = {s: i for i, s in enumerate(states)}
    dim = len(states)
    Tp = sp.zeros(dim, dim)
    for (n, l, m) in states:
        if (n + 1, l, m) in idx:
            Tp[idx[(n + 1, l, m)], idx[(n, l, m)]] = sp.sqrt(sp.Rational((n - l) * (n + l + 1), 4))
    Tm = Tp.T
    N = (-2 * (Tp * Tm - Tm * Tp)).applyfunc(sp.nsimplify)
    for i in range(dim):
        for j in range(dim):
            if i != j:
                assert N[i, j] == 0, (states[i], states[j], N[i, j])
    interior = 0
    for (n, l, m) in states:
        v = sp.nsimplify(N[idx[(n, l, m)], idx[(n, l, m)]])
        if (n + 1, l, m) in idx:
            assert v == n, (n, l, m, v)
            interior += 1
        else:  # top shell: truncated commutator
            assert v == sp.Rational(-(n - 1 - l) * (n + l), 2), (n, l, m, v)
    assert interior > 0


def test_rydberg_levels_from_algebra():
    """The imposed energy-shell map E = -1/(2 N^2), applied to the algebra's
    number operator N, gives E_n = -1/(2 n^2) exactly on interior states;
    the load-bearing content is N's integer spectrum {n}, not the map."""
    max_n = 6
    states, idx, Lp, Tp, Lz = build_operators(max_n)
    Tm = Tp.T
    N = -2.0 * (Tp @ Tm - Tm @ Tp)
    for (n, l, m) in states:
        if (n + 1, l, m) in idx:
            i = idx[(n, l, m)]
            E = -1.0 / (2.0 * N[i, i] ** 2)
            assert abs(E - (-1.0 / (2.0 * n**2))) < 1e-12, (n, l, m, E)
