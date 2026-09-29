"""Paper 16 -- minimum grand angular momentum nu_min of an N-electron atom
(sec:nu_gs, eq:nu_min, tab:mu_values, tab:nu_verified).

Retraction backing (2026-09-28, p16-nu-universal-n-minus-2). Paper 16 previously
claimed nu_min = N - lambda_1 = N - 2 "universally" for every ground state with
S < N/2, and nu = 0 for the fully polarized S = N/2 state. Both are false.

nu_min is the lowest degree K of a harmonic homogeneous polynomial on R^{3N}
(hyperspherical harmonic on S^{3N-1}) whose S_N content contains the spatial irrep
conjugate to the spin diagram [n_up, n_down]. The verified statement is the Pauli
filling count
    nu_min(N, S) = E(n_up) + E(n_down),
    E(m) = sum of the m lowest one-particle polynomial degrees,
where degree d carries (d+1)(d+2)/2 one-particle polynomials (the isotropic 3D
harmonic-oscillator shells 1, 3, 6, 10, ...). It equals N - 2 only for the
lowest-spin state with N <= 8.

Route A (here): exact character generating function
    mult_mu(K) = [t^K] (1 - t^2) (1/N!) sum_sigma chi_mu(sigma) prod_c (1 - t^|c|)^-3,
characters by Murnaghan-Nakayama (checked against the S_3/S_4 tables). For the
L/parity-constrained count, (1 - t^c)^-3 -> 1/det(1 - t^c g^c), g in O(3), with
the O(3) average done by parity projection + the z^L - z^{L+1} coefficient rule.
Route B (here, small N): Young-symmetrizer rank on the explicit monomial space --
no characters used -- as an independent cross-check.
"""
import itertools
from collections import Counter
from fractions import Fraction
from math import factorial

import numpy as np
import pytest


# ------------------------------------------------------------------ helpers
def _partitions(n, maxp=None):
    maxp = n if maxp is None else maxp
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in _partitions(n - k, k):
            yield (k,) + rest


def _mn_char(lam, mu):
    """Murnaghan-Nakayama rule via beta-numbers."""
    lam, mu = list(lam), list(mu)
    if not mu:
        return 1 if sum(lam) == 0 else 0
    r, rest, L = mu[0], mu[1:], len(lam)
    beta = [lam[i] + (L - 1 - i) for i in range(L)]
    bset, total = set(beta), 0
    for b in beta:
        nb = b - r
        if nb < 0 or nb in bset:
            continue
        cnt = sum(1 for x in beta if nb < x < b)
        newbeta = sorted([x for x in beta if x != b] + [nb], reverse=True)
        newlam = [x for x in (newbeta[j] - (L - 1 - j) for j in range(L)) if x > 0]
        total += (-1) ** cnt * _mn_char(newlam, rest)
    return total


def _class_size(mu, n):
    z = 1
    for k, m in Counter(mu).items():
        z *= k ** m * factorial(m)
    return factorial(n) // z


def _spatial_irrep(N, twoS):
    """Conjugate of the spin diagram [n_up, n_down]."""
    nu, nd = (N + twoS) // 2, (N - twoS) // 2
    spin = (nu, nd) if nd > 0 else (nu,)
    return tuple(sum(1 for x in spin if x > j) for j in range(spin[0]))


def _harm_mult_SN(lam, K):
    n = sum(lam)
    acc = [0] * (K + 1)
    for mu in _partitions(n):
        f = [1] + [0] * K
        for c in mu:
            g = [0] * (K + 1)
            for j in range(K // c + 1):
                g[j * c] = (j + 1) * (j + 2) // 2
            f = [sum(f[i] * g[k - i] for i in range(k + 1)) for k in range(K + 1)]
        w = _class_size(mu, n) * _mn_char(lam, mu)
        acc = [a + w * x for a, x in zip(acc, f)]
    acc = [Fraction(a, factorial(n)) for a in acc]
    h = [acc[k] - (acc[k - 2] if k >= 2 else 0) for k in range(K + 1)]
    assert all(x.denominator == 1 for x in h)
    return [int(x) for x in h]


def _kmin_SN(N, twoS, K=30):
    h = _harm_mult_SN(_spatial_irrep(N, twoS), K)
    return next(k for k, x in enumerate(h) if x > 0)


def _harm_mult_SNO3(lam, L, parity, K):
    n = sum(lam)
    tot = [Fraction(0)] * (K + 1)
    for eps in (1, -1):
        acc = [dict() for _ in range(K + 1)]
        for mu in _partitions(n):
            w = _class_size(mu, n) * _mn_char(lam, mu)
            if w == 0:
                continue
            F = [dict() for _ in range(K + 1)]
            F[0] = {0: 1}
            for c in mu:
                G = [dict() for _ in range(K + 1)]
                for m in range(K // c + 1):
                    for b in range(m + 1):
                        for d in range(m - b + 1):
                            e = c * (b - d)
                            G[m * c][e] = G[m * c].get(e, 0) + (eps ** c) ** m
                H = [dict() for _ in range(K + 1)]
                for i in range(K + 1):
                    for j in range(K + 1 - i):
                        for e1, v1 in F[i].items():
                            for e2, v2 in G[j].items():
                                H[i + j][e1 + e2] = H[i + j].get(e1 + e2, 0) + v1 * v2
                F = H
            for k in range(K + 1):
                for e, v in F[k].items():
                    acc[k][e] = acc[k].get(e, 0) + w * v
        wpar = 1 if eps == 1 else parity
        for k in range(K + 1):
            def co(e, kk):
                return acc[kk].get(e, 0) if kk >= 0 else 0
            mL = (co(L, k) - co(L + 1, k)) - (co(L, k - 2) - co(L + 1, k - 2))
            tot[k] += Fraction(wpar * mL, 2 * factorial(n))
    assert all(x.denominator == 1 for x in tot)
    return [int(x) for x in tot]


def _kmin_LP(N, twoS, L, parity, K):
    h = _harm_mult_SNO3(_spatial_irrep(N, twoS), L, parity, K)
    return next((k for k, x in enumerate(h) if x > 0), None)


def _E(m):
    degs, d = [], 0
    while len(degs) < m:
        degs += [d] * ((d + 1) * (d + 2) // 2)
        d += 1
    return sum(degs[:m])


def _fill(N, twoS):
    return _E((N + twoS) // 2) + _E((N - twoS) // 2)


# Hund's-rule ground terms N = 2..12: (2S, L, parity)
GROUND = {2: (0, 0, 1), 3: (1, 0, 1), 4: (0, 0, 1), 5: (1, 1, -1), 6: (2, 1, 1),
          7: (3, 0, -1), 8: (2, 1, 1), 9: (1, 1, -1), 10: (0, 0, 1),
          11: (1, 0, 1), 12: (0, 0, 1)}


# ------------------------------------------------------------------ tests
def test_characters_match_known_tables():
    """Rejects a broken Murnaghan-Nakayama implementation (everything below rests on it)."""
    # S_3, irreps [3],[2,1],[1,1,1] on classes (1^3),(2,1),(3)
    assert [_mn_char(l, (1, 1, 1)) for l in [(3,), (2, 1), (1, 1, 1)]] == [1, 2, 1]
    assert [_mn_char(l, (2, 1)) for l in [(3,), (2, 1), (1, 1, 1)]] == [1, 0, -1]
    assert [_mn_char(l, (3,)) for l in [(3,), (2, 1), (1, 1, 1)]] == [1, -1, 1]
    # S_4 row of [2,2] on classes (1^4),(2,1,1),(2,2),(3,1),(4)
    assert [_mn_char((2, 2), c) for c in [(1,) * 4, (2, 1, 1), (2, 2), (3, 1), (4,)]] \
        == [2, 0, 2, -1, 0]
    for n in range(1, 9):
        assert sum(_mn_char(l, (1,) * n) ** 2 for l in _partitions(n)) == factorial(n)


def test_ground_spin_nu_min_table():
    """tab:nu_verified S_N-only column at the Hund's-rule spin, N = 2..12.
    Rejects nu = N - 2 (wrong for N = 7..12: N 6 not 5, O 7 not 6, Ne 10 not 8, ...)."""
    expected = {2: 0, 3: 1, 4: 2, 5: 3, 6: 4, 7: 6, 8: 7, 9: 8, 10: 10, 11: 12, 12: 14}
    got = {N: _kmin_SN(N, GROUND[N][0]) for N in expected}
    assert got == expected
    assert [N for N in expected if got[N] != N - 2] == [7, 8, 9, 10, 11, 12]


def test_lowest_spin_nu_min_departs_from_N_minus_2_at_N9():
    """Lowest-spin (S = 0 or 1/2) nu_min equals N - 2 through N = 8 only.
    Rejects the 'universal N - 2' reading (N = 9..12 give 8, 10, 12, 14)."""
    got = [_kmin_SN(N, N % 2) for N in range(2, 13)]
    assert got == [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14]
    assert all(got[N - 2] == N - 2 for N in range(2, 9))


def test_fully_polarized_is_not_nu_zero():
    """S = N/2: the spatial irrep is [1^N] (fully ANTISYMMETRIC), so nu_min = E(N) > 0.
    Rejects the retracted 'lambda_1 = N, spatially symmetric, nu = 0' statement."""
    got = [_kmin_SN(N, N) for N in range(2, 13)]
    assert got == [1, 2, 3, 5, 7, 9, 11, 13, 15, 18, 21]
    assert _spatial_irrep(5, 5) == (1, 1, 1, 1, 1)


def test_filling_rule_equals_character_count_all_spins():
    """eq:nu_min: nu_min(N,S) = E(n_up) + E(n_down) for EVERY N <= 12 and every S.
    Rejects any spin-independent rule (e.g. N-2 for all S < N/2): at N = 7 the
    rule gives 5 (S=1/2), 6 (S=3/2), 9 (S=7/2)."""
    for N in range(1, 13):
        for twoS in range(N % 2, N + 1, 2):
            assert _kmin_SN(N, twoS) == _fill(N, twoS), (N, twoS)
    assert [_fill(7, s) for s in (1, 3, 5, 7)] == [5, 6, 7, 9]


def test_production_filling_rule_equals_character_count_all_spins():
    """The PRODUCTION geovac.atomic_classifier.pauli_filling_nu_min (which the
    classifier's nu / mu_free read) equals the Murnaghan-Nakayama character
    count for every N <= 12 and every spin -- not only the hand-pinned Z values.
    Rejects a production shell-capacity error (e.g. the hydrogenic (d+1)^2
    count in place of the oscillator (d+1)(d+2)/2)."""
    from geovac.atomic_classifier import pauli_filling_nu_min
    for N in range(1, 13):
        for twoS in range(N % 2, N + 1, 2):
            n_up, n_down = (N + twoS) // 2, (N - twoS) // 2
            assert pauli_filling_nu_min(n_up, n_down) == _kmin_SN(N, twoS), (N, twoS)


def test_closed_form_and_large_N_table():
    """eq:nu_min closed form E(m) = 3 C(d+2,4) + d (m - C(d+2,3)) and the tab:mu_values
    rows N >= 18 plus the N = 137 (S = 1/2) and Cr values of sec:dirac_limit / Cr paragraph.
    Rejects the retracted 2(N-2)(N-1) numbers (Ar 544, Og 27144, N=137 36720, Cr 1012)."""
    from math import comb
    for m in range(0, 200):
        d = 0
        while comb(d + 3, 3) <= m:
            d += 1
        assert _E(m) == 3 * comb(d + 2, 4) + d * (m - comb(d + 2, 3))
    rows = {18: (0, 26, 1014), 36: (0, 78, 7176), 54: (0, 146, 22338),
            86: (0, 290, 79170), 118: (0, 456, 184224), 137: (1, 570, 279015)}
    for N, (twoS, nu, mu2) in rows.items():
        k = _fill(N, twoS)
        assert k == nu
        assert k * (k + 3 * N - 2) // 2 == mu2 and (k * (k + 3 * N - 2)) % 2 == 0
        assert k != N - 2
    assert _fill(24, 4) == 42 and _fill(24, 6) == 43   # Cr 3d4 4s2 vs 3d5 4s1
    assert _fill(10, 0) * (_fill(10, 0) + 28) // 2 == 190   # Ne isoelectronic series


def test_Li_Be_L_parity_constrained():
    """Physical ground terms: Li 2S (L=0, even) needs nu = 2, not the S_N-only 1
    (degree-1 harmonics are the coordinates themselves: L=1, odd). Be 1S: 2.
    Rejects nu(Li ground state) = 1."""
    assert _kmin_LP(3, 1, 0, 1, 4) == 2
    assert _kmin_LP(3, 1, 1, -1, 4) == 1    # the nu=1 [2,1] harmonics are 2P
    assert _kmin_LP(4, 0, 0, 1, 4) == 2
    assert _harm_mult_SNO3((2, 1), 0, 1, 2)[2] == 2   # r_i.r_j contain [2,1] twice


@pytest.mark.parametrize("N,expected", [(5, 3), (6, 4), (7, 7), (8, 8), (9, 9),
                                        (10, 10), (11, 12), (12, 14)])
def test_L_parity_constrained_B_to_F(N, expected):
    """tab:nu_verified L/parity column, B..Mg ground terms (2P o, 3P e, 4S o, 3P e,
    2P o, 1S e, 2S e, 1S e).  Rejects both N - 2 and the S_N-only count at
    N = 7, 8, 9 (6/7/8 -> 7/8/9), and N - 2 at N = 10, 11, 12 (8/9/10 vs 10/12/14)."""
    twoS, L, par = GROUND[N]
    assert _kmin_LP(N, twoS, L, par, expected + 2) == expected


def test_O3_resolution_is_consistent_with_SN_count():
    """sum_{L,p} (2L+1) m(L,p,K) must equal the S_N-only multiplicity: guards the O(3) code."""
    for lam in [(2, 1), (2, 2), (1, 1, 1)]:
        K = 5
        hs = _harm_mult_SN(lam, K)
        tot = [0] * (K + 1)
        for L in range(K + 1):
            for p in (1, -1):
                tot = [a + (2 * L + 1) * b for a, b in zip(tot, _harm_mult_SNO3(lam, L, p, K))]
        assert tot == hs


# ---------------- Route B: Young symmetrizer on explicit monomials (no characters)
def _young_rank(lam, K):
    N = sum(lam)
    if K < 0:
        return 0
    mons = []
    for c in itertools.combinations_with_replacement(range(3 * N), K):
        e = [0] * (3 * N)
        for i in c:
            e[i] += 1
        mons.append(tuple(e))
    idx = {m: i for i, m in enumerate(mons)}
    tab, k = [], 0
    for r in lam:
        tab.append(list(range(k, k + r)))
        k += r
    cols = [[tab[i][j] for i in range(len(tab)) if j < len(tab[i])] for j in range(lam[0])]

    def group(blocks):
        perms = [list(range(N))]
        for b in blocks:
            perms = [[dict(zip(b, q)).get(i, p[i]) if i in b else p[i] for i in range(N)]
                     for p in perms for q in itertools.permutations(b)]
        return perms

    def sgn(p):
        p, s = list(p), 1
        for i in range(len(p)):
            while p[i] != i:
                j = p[i]
                p[i], p[j] = p[j], p[i]
                s = -s
        return s

    els = Counter()
    for c in group(cols):
        for r in group(tab):
            els[tuple(c[r[i]] for i in range(N))] += sgn(c)
    M = np.zeros((len(mons), len(mons)))
    for j, m in enumerate(mons):
        for p, v in els.items():
            new = [0] * (3 * N)
            for i in range(N):
                for a in range(3):
                    new[3 * p[i] + a] = m[3 * i + a]
            M[idx[tuple(new)], j] += v
    return int(np.linalg.matrix_rank(M))


@pytest.mark.parametrize("lam,kmin", [((2, 1), 1), ((1, 1, 1), 2), ((2, 2), 2),
                                      ((2, 1, 1), 2), ((1, 1, 1, 1), 3)])
def test_route_B_brute_force_agrees(lam, kmin):
    """Independent route (no characters): harmonic multiplicity = rank(P_K) - rank(P_{K-2})
    first becomes nonzero at kmin, and matches Route A there.
    Rejects nu([1,1,1]) = 0 and nu([1^4]) = 0 (the retracted S = N/2 statement)."""
    ranks = [_young_rank(lam, K) for K in range(kmin + 1)]
    harm = [ranks[K] - (ranks[K - 2] if K >= 2 else 0) for K in range(kmin + 1)]
    assert all(h == 0 for h in harm[:kmin]) and harm[kmin] > 0
    assert harm[kmin] == _harm_mult_SN(lam, kmin)[kmin]
