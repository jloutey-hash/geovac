"""Paper 38 lem:lifted_state (a)-(b) on the TRUE spinor window.

(a) xi = (h (x) chi)/sqrt(Z), h = sum_{j <= (n_max-1)/2} sqrt(2j+1) chi_j, lies
    exactly in the span of the complete Camporesi-Higuchi eigenspaces with
    |lambda| <= n_max + 1/2 (no spillover), for any unit spinor chi.
(b) <xi, P M_f P xi> = int f K_{n_max} dg: on the matrix coefficient
    f = D^J_{M M'} the right side is delta_{M M'} c_J / ((2J+1) Z), with
    c_J = sum_{j1, j2 <= j_max, triangle(j1, j2, J)} sqrt((2j1+1)(2j2+1))
    (Paper 38 Lemma L2(b)).

Window, Dirac and multipliers are the Clebsch-Gordan build of
test_p38_true_ch_window.py (Route A).  In its basis e^j_{m m' s} =
sqrt(2j+1) D^j_{m m'} (x) e_s the amplitude h is sum_{j, m} e^j_{m m}.
The right side of (b) is computed from the character triangle rule, not
from the matrices.
"""
from __future__ import annotations

import itertools
from fractions import Fraction

import numpy as np
import pytest

from test_p38_true_ch_window import _build, _mult_op


def _xi(n_max: int, spinor: np.ndarray) -> np.ndarray:
    labels, idx, _, _, _ = _build(n_max)
    jmax = Fraction(n_max - 1, 2)
    v = np.zeros(len(labels), complex)
    for (j, m, mp, s) in labels:
        if j <= jmax and m == mp:
            v[idx[(j, m, mp, s)]] = spinor[s]
    z = n_max * (n_max + 1) / 2
    return v / np.sqrt(z)


def _projector(n_max: int) -> np.ndarray:
    _, _, _, ev, U = _build(n_max)
    V = U[:, np.abs(ev) <= n_max + 0.5 + 1e-9]
    return V @ V.conj().T


def _c(J: Fraction, n_max: int) -> float:
    js = [Fraction(k, 2) for k in range(n_max)]
    return sum(np.sqrt(float((2 * a + 1) * (2 * b + 1))) for a in js for b in js
               if abs(a - b) <= J <= a + b and (a + b + J).denominator == 1)


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_window_inclusion_is_exact(n_max):
    rng = np.random.default_rng(n_max)
    spinor = rng.normal(size=2) + 1j * rng.normal(size=2)
    spinor /= np.linalg.norm(spinor)
    xi = _xi(n_max, spinor)
    assert abs(np.linalg.norm(xi) - 1) < 1e-12
    assert np.linalg.norm(xi - _projector(n_max) @ xi) < 1e-12


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_lifted_state_is_the_fejer_measure(n_max):
    labels, idx, _, _, _ = _build(n_max)
    P = _projector(n_max)
    xi = _xi(n_max, np.array([1.0, 0.0]))
    z = n_max * (n_max + 1) / 2
    checked = 0
    for N in range(1, 2 * n_max + 1):
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        for M, Mp in itertools.product(Ms, Ms):
            A = _mult_op(J, M, Mp, labels, idx)
            lhs = xi.conj() @ (P @ A @ P) @ xi
            rhs = _c(J, n_max) / (N * z) if M == Mp else 0.0
            assert abs(lhs - rhs) < 1e-12, (n_max, J, M, Mp, lhs, rhs)
            checked += 1
    assert checked == sum(N * N for N in range(1, 2 * n_max + 1))
