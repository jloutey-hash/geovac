"""Backing test for Paper 11 Sec. "Heteronuclear extension" -- the eta-label
continuity paragraph added 2026-09-06.

The eta-equation (Paper 11 eq:eta_equation, m = 0) in the normalized Legendre
basis is H = -l(l+1) delta + c^2 <eta^2> + b <eta>.  Its ground eigenvector's
participation deficit D = 1 - max_l |v_l|^2 measures how far the eta-label has
drifted from a pure Legendre (united-atom l) label.

Claim: the heteronuclear term b*eta couples l <-> l+-1 and turns the mixing on
at FIRST order in the coupling (D ~ b^2 ~ R^2), whereas the homonuclear
c^2 eta^2 term couples l <-> l+-2 and turns it on at second order
(D ~ (c^2)^2 ~ R^4).  Rejected wrong answer: both onsets share one exponent.

The R-parametrization uses b = R(Z_B - Z_A) and c^2 = kappa R^2 with a fixed
united-atom kappa; the ONSET EXPONENT is independent of kappa, which is why
this test needs no solver.  The non-saturation of the heteronuclear deficit at
large R (Paper 11 text) is a solver-level statement backed by the CHANGELOG
v5.10.2 chronicle, not by this fast test.
"""

from __future__ import annotations

import numpy as np
from scipy.special import eval_legendre, roots_legendre

_X, _W = roots_legendre(80)


def _legendre_matrices(N: int):
    P = np.array([eval_legendre(l, _X) * np.sqrt((2 * l + 1) / 2.0) for l in range(N)])
    eta1 = (P * _W * _X) @ P.T          # <P_l | eta | P_l'>
    eta2 = (P * _W * _X**2) @ P.T       # <P_l | eta^2 | P_l'>
    return eta1, eta2


def _deficit(c2: float, b: float, N: int = 30) -> float:
    eta1, eta2 = _legendre_matrices(N)
    ls = np.arange(N)
    H = -np.diag(ls * (ls + 1)).astype(float) + c2 * eta2 + b * eta1
    w, V = np.linalg.eigh(H)
    v = V[:, np.argmax(w)]              # ground state = smallest A = largest -A
    return 1.0 - float(np.max(v**2))


def _onset_exponent(kappa: float, dZ: float) -> float:
    Rs = np.array([0.04, 0.08, 0.16])
    D = np.array([_deficit(kappa * R**2, dZ * R) for R in Rs])
    return float(np.polyfit(np.log(Rs), np.log(D), 1)[0])


def test_legendre_matrix_sanity():
    eta1, eta2 = _legendre_matrices(6)
    assert abs(eta1[0, 1] - 1.0 / np.sqrt(3.0)) < 1e-12        # <P0|eta|P1>
    assert abs(eta2[0, 0] - 1.0 / 3.0) < 1e-12                  # <P0|eta^2|P0>
    assert np.allclose(eta1, eta1.T) and np.allclose(eta2, eta2.T)


def test_heteronuclear_onset_first_order_homonuclear_second():
    # HeH2+: Z_A = 2, Z_B = 1 -> b = -R ; united-atom Li2+ E = -4.5 -> c^2 = 2.25 R^2
    p_hetero = _onset_exponent(kappa=2.25, dZ=-1.0)
    # H2+: b = 0 ; united-atom He+ E = -2 -> c^2 = R^2
    p_homo = _onset_exponent(kappa=1.0, dZ=0.0)
    assert abs(p_hetero - 2.0) < 0.15, p_hetero
    assert abs(p_homo - 4.0) < 0.15, p_homo
    assert p_homo - p_hetero > 1.5              # the two onsets are NOT one exponent
    # kappa-independence of the heteronuclear onset (it is the b-term's order)
    assert abs(_onset_exponent(kappa=1.0, dZ=-1.0) - p_hetero) < 0.1
