"""Verification of Paper 51 Theorem thm:dS_vacuum (Section G5).

Claim: the decompactified spectral-action density on S^3_R x R_tau with
Gaussian cutoff is

    s(R, Lambda) = (R^3/2) Lambda^4 - (R/4) Lambda^2,

with extremum R_crit Lambda = 1/sqrt(6) and minimum s_min = -Lambda/(6 sqrt 6).

Rebuilt 2026-11-07 (group5 guard pass).  The previous version hard-coded
the density with the rank-2 coefficients (R^3/4, R/8) and so certified the
withdrawn s_min = -Lambda/(12 sqrt 6) by restating it.  This version does
not take the density from the paper: it sums Tr f(D^2/Lambda^2) over the
Camporesi--Higuchi spectrum (|lambda_n| = (n+3/2)/R, multiplicity
2(n+1)(n+2)) and the continuous tau-momentum, with the trace over the
rank-4 bundle C^2 (x) C^2 of the operator D_{S^3} (x) sigma_1 +
I (x) sigma_2 (-i d_tau).  The coefficients, the extremum and s_min are then
read off the computed numbers.

What wrong answer this rejects: the rank-2 trace (half the density, the
withdrawn -Lambda/(12 sqrt 6)); any other prefactor or sign of either term.
"""

import numpy as np
import sympy as sp
from scipy import integrate


def _sigma_space_trace_factor() -> int:
    """Trace of the identity on the sigma factor, read off D^2.

    D = D3 (x) s1 + I (x) s2 * p.  Its square is D3^2 (x) I + I (x) p^2 (x) I
    because s1 s2 + s2 s1 = 0; so each (S^3 mode, momentum) pair appears
    with multiplicity dim of the sigma space.  The rank-4 bundle itself is
    the paper's input (sec:g2); this documents the trace, it does not test
    the bundle choice.
    """
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    lam, p = 1.7, 0.6  # generic eigenvalue and momentum
    D = lam * s1 + p * s2
    D2 = D @ D
    assert np.allclose(D2, (lam ** 2 + p ** 2) * np.eye(2))
    return int(round(np.trace(D2).real / (lam ** 2 + p ** 2)))


def spectral_action_density(f, R: float, Lam: float, nmax: int = 4000) -> float:
    """S/beta on S^3_R x R_tau, summed over the CH spectrum (beta -> infinity)."""
    rank = _sigma_space_trace_factor()
    total = 0.0
    for n in range(nmax):
        l2 = ((n + 1.5) / R) ** 2
        if l2 / Lam ** 2 > 60.0:
            break
        mult = 2 * (n + 1) * (n + 2)
        half, _ = integrate.quad(lambda w: f((l2 + w * w) / Lam ** 2), 0.0, np.inf, limit=200)
        total += mult * 2.0 * half / (2.0 * np.pi)
    return rank * total


def _coefficients(Lam: float):
    """Fit s(R) = A R^3 + B R from the computed density at two radii."""
    f = lambda x: np.exp(-x)
    R1, R2 = 1.0, 2.0
    s1, s2 = spectral_action_density(f, R1, Lam), spectral_action_density(f, R2, Lam)
    M = np.array([[R1 ** 3, R1], [R2 ** 3, R2]])
    A, B = np.linalg.solve(M, np.array([s1, s2]))
    return A, B


def test_density_coefficients_from_spectrum():
    """Computed A = Lambda^4/2 and B = -Lambda^2/4 (rank-4 trace)."""
    Lam = 8.0
    A, B = _coefficients(Lam)
    assert abs(A / Lam ** 4 - 0.5) < 1e-9, A / Lam ** 4
    assert abs(B / Lam ** 2 + 0.25) < 1e-6, B / Lam ** 2


def test_two_term_exactness_third_radius():
    """The two-term form predicts a third radius it was not fitted to."""
    Lam = 8.0
    A, B = _coefficients(Lam)
    R3 = 0.7
    s3 = spectral_action_density(lambda x: np.exp(-x), R3, Lam)
    assert abs(s3 - (A * R3 ** 3 + B * R3)) < 1e-9 * abs(s3)


def test_dS_extremum_and_smin_from_computed_density():
    """thm:dS_vacuum: R_crit Lambda = 1/sqrt 6, s_min = -Lambda/(6 sqrt 6)."""
    Lam = 8.0
    A, B = _coefficients(Lam)
    R_crit = np.sqrt(-B / (3.0 * A))
    s_min = A * R_crit ** 3 + B * R_crit
    assert abs(R_crit * Lam - 1.0 / np.sqrt(6.0)) < 1e-7
    assert abs(s_min / Lam - (-1.0 / (6.0 * np.sqrt(6.0)))) < 1e-7, s_min / Lam
    # it is a minimum in R
    assert 6.0 * A * R_crit > 0.0


def test_dS_closed_form_algebra():
    """Symbolic algebra on the stated density: extremum and minimum value."""
    R, Lam = sp.symbols("R Lambda", positive=True)
    s = R ** 3 * Lam ** 4 / 2 - R * Lam ** 2 / 4
    crit = [c for c in sp.solve(sp.diff(s, R), R) if c.is_positive]
    assert len(crit) == 1
    assert sp.simplify(crit[0] - 1 / (sp.sqrt(6) * Lam)) == 0
    assert sp.simplify(s.subs(R, crit[0]) + Lam / (6 * sp.sqrt(6))) == 0
