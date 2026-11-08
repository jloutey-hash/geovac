"""Paper 28 eq:s3_closed_form / eq:s3_value at cutoffs past the existing pins
(2026-11-07 guard pass; CLAUDE.md guard-asymptotics rule).

The 25-term corrected-W closed form for S^(3) (a = 4, p = 1) is checked in
tests/test_s3_decomposition.py against the 60-digit canonical value and an
independent t-value-free route, and its rigorous direct-sum bracket is
checked at N = 20000.  This file re-checks the closed form against the
DEFINITION at two larger cutoffs, N = 2e5 and N = 2e6 (100x past the old
bracket), where the bracket width falls to ~5e-3:

    S^(3) = sum_{m>=0} phi(m) L(m)^2,  L(m) = 2 sum_{n>=0} min(n, m) phi(n),
    phi(n) = 8 (o^2 - 1) / o^4,  o = 2n + 3.

Lower bound: the positive partial sum to N.  Upper bound: the partial sum plus
the integral tail bound 32 [(ln N + 1)^2 + 2 (ln N + 1) + 2] / N (from
L(m) <= 4 (ln m + 1) and phi(m) <= 2/(m+1)^2; Paper 28 sec:three_loop).
L(m) is evaluated with the full-range Hurwitz tail
T(m) = sum_{n>=m} phi(n) = 2 zeta(2, m + 3/2) - zeta(4, m + 3/2)/2, so the
partial sum is the exact definition truncated in m only.  No geovac import.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from mpmath import mp, mpf, nsum, inf, zeta, log, pi
from scipy.special import zeta as hz


def _s3_closed_form_float() -> float:
    """eq:s3_closed_form (25 terms) at 30 dps; zeta(5,3) by a log-free
    Hurwitz sum."""
    with mp.workdps(30):
        z3, z5 = zeta(3), zeta(5)
        z53 = nsum(lambda n: n ** mpf(-5) * (z3 - zeta(3, n)), [2, inf])
        ln2 = log(2)
        F = lambda a, b: mpf(a) / b
        v = (128 * pi**2 * ln2**2 - F(32, 3) * pi**4 * ln2**2
             + 64 * pi**2 * ln2 - F(304, 3) * pi**4 * ln2 + 16 * pi**6 * ln2
             - F(2, 3) * pi**8 * ln2
             - 96 * pi**2 * ln2 * z3 + 16 * pi**4 * ln2 * z3
             - 80 * pi**2 * ln2 * z5
             + 20 * pi**2 * z3**2 - 7 * pi**4 * z3**2 + 70 * pi**2 * z3 * z5
             - 24 * pi**2 * z3 + 40 * pi**4 * z3 - 9 * pi**6 * z3
             + pi**8 * z3 / 2
             - 20 * pi**2 * z5 + 30 * pi**4 * z5 - F(5, 2) * pi**6 * z5
             + 6 * pi**2 * z53
             + F(8, 3) * pi**4 + F(374, 45) * pi**6 - F(8231, 3780) * pi**8
             + F(124571, 680400) * pi**10 - pi**12 / 192)
        return float(v)


def _bracket(N: int):
    m = np.arange(0, N + 1, dtype=np.float64)
    o = 2 * m + 3
    ph = 8 * (o * o - 1) / o ** 4
    T = 2 * hz(2, m + 1.5) - 0.5 * hz(4, m + 1.5)        # sum_{n>=m} phi(n)
    nph = m * ph
    cum_nph = np.concatenate(([0.0], np.cumsum(nph)[:-1]))  # sum_{n<m} n phi(n)
    L = 2 * (cum_nph + m * T)
    lo = math.fsum((ph * L * L).tolist())
    u = math.log(N) + 1
    hi = lo + 32 * (u * u + 2 * u + 2) / N
    return lo, hi


S3_CANONICAL = 40.6553875635854225          # eq:s3_value
S3_OLD_COUNT = 31.5725612075120228          # pre-correction channel count


def test_closed_form_matches_printed_value() -> None:
    """eq:s3_closed_form evaluates to eq:s3_value (40.65538756358542...).
    Rejects: a closed form whose printed coefficients do not sum to the
    printed value (e.g. 374/45 -> 374/44)."""
    assert abs(_s3_closed_form_float() - S3_CANONICAL) < 1e-13


@pytest.mark.parametrize("N", [200_000, 2_000_000])
def test_closed_form_inside_direct_bracket_past_data(N: int) -> None:
    """The closed form lies inside the rigorous direct-sum bracket at
    N = 2e5 and 2e6, and the bracket excludes the old-count value 31.5726
    and every value more than the bracket width away.

    Rejects: the old-count S^(3) = 31.5725612..., and any closed form off
    by more than ~5e-3 (at N = 2e6), i.e. 100x tighter than the
    N = 20000 bracket in tests/test_s3_decomposition.py."""
    lo, hi = _bracket(N)
    cf = _s3_closed_form_float()
    assert lo <= cf <= hi, (N, lo, hi, cf)
    assert not (lo <= S3_OLD_COUNT <= hi)
    # the bracket really is tight at the larger cutoff
    if N == 2_000_000:
        assert hi - lo < 6e-3
        # the true tail S - S(N) is positive and below the bound
        assert 0 < cf - lo < hi - lo
