"""
Paper 8 (Bond Sphere / Sturmian) — Theorem 1 (sigma-bond selection rules) backing.

Paper 8 Theorem 1 (`thm:selection`) states two SO(4) Wigner-D identities for the
n=2 representation, under the bond rotation e^{i gamma A_y}:

    D^{(2)}_{(1,0),(1,0)}(gamma) = 1   for all gamma     (Eq. D2_1010)
    D^{(2)}_{(0,0),(1,0)}(gamma) = 0   for all gamma     (Eq. D2_0010)

Physical content: the 2p_0 (sigma) orbital is a TRANSPARENT mode of the SO(4)
bond channel (perfectly transmitted at every bond angle), and the bond rotation
alone cannot mix 2s into 2p_0 -- so sigma-bonds come from p_z alignment, not
from geometric s-p coupling (the s-p hybridization is energy-optimization, not
the Runge-Lenz rotation).

The paper's "Numerical verification" paragraph cites a 10,000-point gamma grid
with both deviations < 1e-14.  That verification originally lived only in a
SPECULATIVE diagnostic, debug/archive/test_harmonic_phase_lock.py -- a probe of
the (failed, no-hits) hypothesis that D-matrix critical points predict LiH R_eq.
The two identities were an incidental byproduct of that probe; this test is the
purpose-built, default-collected backing for Theorem 1, decoupled from the
archived dead-end.

This test reproduces the paper's exact numerical statement on the same grid AND
checks the closed forms from the in-paper proof (cos^2 + sin^2 = 1 for the
diagonal element; the antisymmetric s-p combination cancels identically).

Provenance: SYMBOLIC PROOF (closed-form identity) + MEASURED (10,000-point grid,
machine precision).
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.wigner_so4 import wigner_D_so4


# Paper 8 grid: 10,000 values of gamma in [0.01, pi - 0.01].
GAMMA_MIN = 0.01
GAMMA_MAX = np.pi - 0.01
N_POINTS = 10_000
TOL = 1e-14  # the paper's stated bound

# (n, l', m', l, m) for the two Theorem-1 elements.
ELEM_DIAG = (2, 1, 0, 1, 0)   # D^{(2)}_{(1,0),(1,0)}
ELEM_OFF = (2, 0, 0, 1, 0)    # D^{(2)}_{(0,0),(1,0)}


@pytest.fixture(scope="module")
def gamma_grid() -> np.ndarray:
    return np.linspace(GAMMA_MIN, GAMMA_MAX, N_POINTS)


def test_transparent_mode_identity_D2_1010(gamma_grid):
    """Eq. (D2_1010): D^{(2)}_{(1,0),(1,0)}(gamma) = 1 for all gamma.

    The 2p_0 (sigma) orbital is transmitted without attenuation at every bond
    angle.  Reproduces the paper's |D - 1| < 1e-14 on the 10,000-point grid.
    """
    vals = np.array([wigner_D_so4(*ELEM_DIAG, g) for g in gamma_grid])
    assert np.max(np.abs(vals - 1.0)) < TOL


def test_no_sp_mixing_identity_D2_0010(gamma_grid):
    """Eq. (D2_0010): D^{(2)}_{(0,0),(1,0)}(gamma) = 0 for all gamma.

    The bond rotation alone cannot scatter 2s into 2p_0.  Reproduces the
    paper's |D| < 1e-14 on the 10,000-point grid.
    """
    vals = np.array([wigner_D_so4(*ELEM_OFF, g) for g in gamma_grid])
    assert np.max(np.abs(vals)) < TOL


def test_diagonal_matches_closed_form(gamma_grid):
    """In-paper proof: D^{(2)}_{(1,0),(1,0)} = cos^2(g/2) + sin^2(g/2) = 1.

    Cross-checks the computed element against the closed form term-by-term
    (not just against the constant 1), tying the numerical value to the
    Clebsch-Gordan derivation in the proof of Theorem 1.
    """
    vals = np.array([wigner_D_so4(*ELEM_DIAG, g) for g in gamma_grid])
    closed = np.cos(gamma_grid / 2.0) ** 2 + np.sin(gamma_grid / 2.0) ** 2
    assert np.allclose(vals, closed, atol=TOL, rtol=0.0)


def test_offdiagonal_matches_closed_form(gamma_grid):
    """In-paper proof: D^{(2)}_{(0,0),(1,0)} = cos(g/2)sin(g/2) - sin(g/2)cos(g/2) = 0.

    The antisymmetric s-p combination cancels identically; cross-check the
    computed element against this difference form.
    """
    vals = np.array([wigner_D_so4(*ELEM_OFF, g) for g in gamma_grid])
    half = gamma_grid / 2.0
    closed = np.cos(half) * np.sin(half) - np.sin(half) * np.cos(half)
    assert np.allclose(vals, closed, atol=TOL, rtol=0.0)


def test_identities_are_genuinely_distinct():
    """Guard: the two elements are not the same quantity collapsed -- the
    diagonal is the constant 1, the off-diagonal the constant 0 (so a routine
    returning a single constant for both would fail one of these)."""
    g = 1.234  # arbitrary interior angle
    diag = wigner_D_so4(*ELEM_DIAG, g)
    off = wigner_D_so4(*ELEM_OFF, g)
    assert abs(diag - 1.0) < TOL
    assert abs(off) < TOL
    assert abs(diag - off) > 0.5
