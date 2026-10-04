"""Paper 38 Lemma L2(d.ii) + Open question (i): the Fejer-moment expansion

    n * gamma_n = c * log n + b + O(log n / n),   c = 4/pi,  b = 4.10932146...

(dual-Coxeter normalisation).  The leading constant c is the Fejer-moment
constant (a kernel constant, not a convergence rate); b is the next
coefficient, closed form [OPEN].

Route: the exact sum rule gamma_n = pi - 4 T_n / (pi Z_n) (Paper 38 eq.
gamma_sum_rule), summed here in O(n) memory by odd offset d, cross-checked
against the production mpmath sum rule, then a least-squares fit of
n * gamma_n over n = 200..6400 with the subleading terms log n/n, 1/n,
log n/n^2 free.  Neither c nor b is put in: both come out of the fit.
"""
from __future__ import annotations

import math

import numpy as np

from geovac.central_fejer_su2 import gamma_n_via_sum_rule

NS = (200, 400, 800, 1600, 3200, 6400)


def _n_gamma(n: int) -> float:
    """n * gamma_n from the exact sum rule, summed over odd offsets d = k2 - k1."""
    k = np.arange(1, n + 1, dtype=float)
    s = np.sqrt(k)
    t = 0.0
    for d in range(1, n, 2):
        a = k[: n - d]
        t += 2.0 * float(np.sum(s[: n - d] * s[d:] * (1.0 / d ** 2 - 1.0 / (2 * a + d) ** 2)))
    z = n * (n + 1) / 2
    return n * (math.pi - 4.0 * t / (math.pi * z))


def _fit() -> np.ndarray:
    y = np.array([_n_gamma(n) for n in NS])
    a = np.array([[math.log(n), 1.0, math.log(n) / n, 1.0 / n, math.log(n) / n ** 2] for n in NS])
    return np.linalg.lstsq(a, y, rcond=None)[0]


def test_odd_offset_sum_matches_production_sum_rule() -> None:
    for n in (7, 50, 200):
        assert abs(_n_gamma(n) - n * float(gamma_n_via_sum_rule(n))) < 1e-9


def test_leading_constant_is_four_over_pi() -> None:
    """c = 4/pi to 1e-5 (measured residual ~5e-7); the circle constant 2/pi and
    Kac's 2 sqrt 2/pi are excluded by orders of magnitude."""
    c = _fit()[0]
    assert abs(c - 4 / math.pi) < 1e-5, c
    assert abs(c - 2 * math.sqrt(2) / math.pi) > 0.3


def test_next_coefficient_b() -> None:
    """b = 4.10932146... (Paper 38 open question (i)), pinned to 2e-5."""
    b = _fit()[1]
    assert abs(b - 4.1093215) < 2e-5, b
