"""Paper 24 sec:large_spin -- releasing spin by the large-spin contraction.

Record: docs/predictions/2026-10-18_large_spin.md (LS1-LS5, committed before the computation in
fb4fdd85). Written as its own step after the result; it does not import the computation script.
Each test names the wrong answer it rejects.

  LS1  A_S^dagger = S_- / sqrt(2S): element sqrt(n+1) sqrt(1 - n/(2S)), [A_S, A_S^dagger] = 1 - n/S
       on |n>, symbolic in S and n. Rejects 1 - n/(2S) and an exact oscillator at finite S.
  LS2  the S_x eigenvector components by an integer (Krawtchouk) three-term recurrence, a route
       independent of Wigner's sum: they solve S_x v = m v, and their squares are rational.
  LS3  p_S(x) = sqrt(S) C(2S, S+m) / 4^S is rational at S = N^2 and tends to pi^(-1/2) e^(-x^2)
       with c_1 = -1/8 + x^2/2 - x^4/6, no S^(-1/2) term, and the record's measured c_2.
  LS4  S^(1/4) <S, S-n | X_S = x> tends to psi_n(x) at rate 1/S (top component positive).
  LS5  the sphere's coherent-state measure contracts to the Bargmann measure, symbolically.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb
from typing import List

import mpmath as mp
import sympy as sp

# ---------------------------------------------------------------- LS1


def test_ls1_scaled_element_and_commutator_symbolic():
    """Rejects a commutator of 1 - n/(2S) (the error of dropping one of the two terms) and the
    claim that the contraction is exact at finite S."""
    S, n = sp.symbols("S n", positive=True)
    m = S - n
    el2 = (S + m) * (S - m + 1) / (2 * S)  # |<n+1| A^dagger |n>|^2, with S_- : m -> m - 1
    assert sp.simplify(el2 - (n + 1) * (1 - n / (2 * S))) == 0
    el2_down = el2.subs(n, n - 1)  # |<n| A^dagger |n-1>|^2 = <n| A^dagger A |n>
    comm = sp.simplify(el2 - el2_down)
    assert sp.simplify(comm - (1 - n / S)) == 0
    assert sp.simplify(comm - (1 - n / (2 * S))) != 0
    assert sp.limit(comm, S, sp.oo) == 1


# ---------------------------------------------------------------- LS2


def _u(S: int, m: int, nmax: int) -> List[int]:
    """(n+1) u_{n+1} = 2m u_n - (2S - n + 1) u_{n-1}, u_0 = 1: integer Krawtchouk values."""
    u = [Fraction(1)]
    prev = Fraction(0)
    for nn in range(nmax):
        nxt = (2 * m * u[nn] - (2 * S - nn + 1) * prev) / (nn + 1)
        prev = u[nn]
        u.append(nxt)
    assert all(v.denominator == 1 for v in u)
    return [int(v) for v in u]


def _comp_sq(S: int, m: int, n: int) -> Fraction:
    """<S, S-n | S_x = m>^2 = C(2S, S+m) u_n^2 / (4^S C(2S, n)), exactly rational."""
    u = _u(S, m, n)[n]
    return Fraction(comb(2 * S, S + m) * u * u, 4**S * comb(2 * S, n))


def _comp(S: int, m: int, n: int) -> mp.mpf:
    u = _u(S, m, n)[n]
    sq = _comp_sq(S, m, n)
    return mp.sign(u) * mp.sqrt(mp.mpf(sq.numerator) / sq.denominator)


def test_ls2_recurrence_vectors_are_sx_eigenvectors_with_rational_squares():
    """Rejects a recurrence or normalization that is not the S_x eigenvector (residual check at
    50 digits). The squares are rational because the recurrence stays in the integers, asserted
    in _u."""
    mp.mp.dps = 50
    for S in (3, 6):
        dim = 2 * S + 1
        for m in range(-S, S + 1):
            v = [_comp(S, m, n) for n in range(dim)]
            norm = mp.fsum(c * c for c in v)
            assert abs(norm - 1) < mp.mpf(10) ** -45
            for n in range(dim):  # (S_x v)_n, state S - n; S_- (S-n -> S-n-1) = sqrt((2S-n)(n+1))
                acc = mp.mpf(0)
                if n > 0:
                    acc += mp.sqrt(n * (2 * S - n + 1)) * v[n - 1] / 2
                if n < dim - 1:
                    acc += mp.sqrt((n + 1) * (2 * S - n)) * v[n + 1] / 2
                assert abs(acc - m * v[n]) < mp.mpf(10) ** -45


# ---------------------------------------------------------------- LS3


def _p(N: int, x: Fraction) -> Fraction:
    S = N * N
    m = x * N
    assert m.denominator == 1
    return Fraction(N * comb(2 * S, S + int(m)), 4**S)


def _fit(Ns: List[int], vals: List[mp.mpf], exps: List[mp.mpf]) -> List[mp.mpf]:
    A = mp.matrix([[mp.mpf(N) ** (-2 * e) for e in exps] for N in Ns])
    sol = mp.lu_solve(A, mp.matrix(vals))
    return [sol[j] for j in range(len(exps))]


NS = [40, 60, 80, 120, 160, 240, 320]


def _c1(x: mp.mpf) -> mp.mpf:
    return mp.mpf(-1) / 8 + x**2 / 2 - x**4 / 6


def _c2(x: mp.mpf) -> mp.mpf:
    return (mp.mpf(1) / 128 - 11 * x**2 / 48 + 19 * x**4 / 48 - 3 * x**6 / 20 + x**8 / 72)


def test_ls3_density_limit_c1_c2_and_no_half_powers():
    """Rejects a different first correction (for example -x^4/12 in place of -x^4/6), a wrong
    printed c_2, a limit other than pi^(-1/2) e^(-x^2), and an S^(-1/2) term."""
    mp.mp.dps = 60
    exps_int = [mp.mpf(j) for j in range(len(NS))]
    for x in (Fraction(0), Fraction(1), Fraction(2)):
        xm = mp.mpf(x.numerator) / x.denominator
        g = mp.e ** (-xm**2) / mp.sqrt(mp.pi)
        vals = [mp.mpf(_p(N, x).numerator) / _p(N, x).denominator / g for N in NS]
        co = _fit(NS, vals, exps_int)
        assert abs(co[0] - 1) < mp.mpf(10) ** -15
        assert abs(co[1] - _c1(xm)) < mp.mpf(10) ** -10
        assert abs(co[2] - _c2(xm)) < mp.mpf(10) ** -6
        # an S^(-1/2) term, if present, would show in a fit that allows it
        exps_half = [mp.mpf(j) / 2 for j in range(len(NS))]
        coh = _fit(NS, vals, exps_half)
        assert abs(coh[1]) < mp.mpf(10) ** -6 and abs(coh[2] - _c1(xm)) < mp.mpf(10) ** -3


# ---------------------------------------------------------------- LS4


def _psi(n: int, x: mp.mpf) -> mp.mpf:
    return mp.hermite(n, x) * mp.e ** (-x**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))


def test_ls4_components_tend_to_hermite_functions_at_rate_one_over_S():
    """Rejects a limit other than psi_n(x) (for example a missing 2^n in the normalization, or
    the argument x / sqrt(2)), and a rate O(1/sqrt(S)): S (ratio - 1) must settle."""
    mp.mp.dps = 50
    x = Fraction(1)
    xm = mp.mpf(1)
    for n in range(4):
        scaled = []
        for N in (40, 80, 160):
            S = N * N
            r = _comp(S, N, n) * mp.sqrt(N) / _psi(n, xm)  # S^(1/4) = sqrt(N); m = x N
            scaled.append((r - 1) * S)
        # O(1/S): (r - 1) S converges; a 1/sqrt(S) term would double it with each doubling of N
        assert abs(scaled[2] - scaled[1]) < abs(scaled[1] - scaled[0])
        assert abs(scaled[2] - scaled[1]) < mp.mpf("0.05") * (1 + abs(scaled[2]))
        assert abs(scaled[2]) < 10


# ---------------------------------------------------------------- LS5


def test_ls5_coherent_measure_contracts_to_bargmann():
    """Rejects a contracted measure other than pi^(-1) e^(-|w|^2), for example e^(-|w|^2 / 2),
    and a limit of C(2S, n)/(2S)^n other than 1/n!."""
    S = sp.symbols("S", positive=True)
    r2 = sp.symbols("r2", positive=True)  # |w|^2
    # (2S+1)/(4 pi) dOmega, with z = w / sqrt(2S) and dOmega = 4 d^2z/(1+|z|^2)^2:
    dens = (2 * S + 1) / (4 * sp.pi) * 4 / (2 * S) / (1 + r2 / (2 * S)) ** 2
    # times the coherent-state weight |<z|n>|^2 / |z|^(2n) C(2S, n)^(-1) -> (1 + |z|^2)^(-2S)
    lim = sp.limit(dens * (1 + r2 / (2 * S)) ** (-2 * S), S, sp.oo)
    assert sp.simplify(lim - sp.exp(-r2) / sp.pi) == 0
    assert sp.simplify(lim - sp.exp(-r2 / 2) / sp.pi) != 0
    for n in range(6):
        c = sp.binomial(2 * S, n) / (2 * S) ** n
        assert sp.simplify(sp.limit(sp.expand_func(c), S, sp.oo) - sp.Rational(1, sp.factorial(n))) == 0
