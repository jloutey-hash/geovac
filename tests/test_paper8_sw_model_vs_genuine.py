"""Guard: Paper 8, sec. "The Single-n D-Matrix Model of the SW Integrals".

Paper: papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex,
Eq. (eq:sw_exact) and the paragraph after it.  The paper states that the
single-n D-matrix expression is a MODEL, not an identity:

  * the model gives <1s_A|1/r_B|1s_A> = 1/p0 = 1 at every R (D^(1) = 1);
  * the genuine integral is 1/R - (1 + 1/R) e^{-2R}, printed as
    0.896, 0.610, 0.328, 0.125 at R = 0.5, 1.4, 3.015, 8 bohr (p0 = Z = 1);
  * <1s_A|1/r_B|2s_A>, which the block-diagonal D^(n) sets to zero, is nonzero.

The genuine integrals are computed here by independent radial quadrature
(the exact monopole reduction of 1/|r - R| for a spherically symmetric
product density: int rho(r) / max(r, R) d^3r), not from the closed form,
and then compared with the closed form and the printed values.

Fire-test: return the model value (1/p0, resp. 0) from ``genuine_element``;
every test must fail -- this is the "exact" reading the section retracted.
"""
from __future__ import annotations

import math

import pytest
from scipy.integrate import quad

P0 = 1.0          # shared Sturmian momentum, p0 = Z_A = Z_B = 1
PRINTED = {0.5: 0.896, 1.4: 0.610, 3.015: 0.328, 8.0: 0.125}


def _sturmian_s(n: int, r: float) -> float:
    """L2-normalized s-wave Coulomb-Sturmian radial at scale k = p0.

    Equals the hydrogenic R_n0 at Z = n p0 (rho = 2 p0 r).
    """
    x = 2.0 * P0 * r
    if n == 1:
        return 2.0 * P0 ** 1.5 * math.exp(-x / 2)
    if n == 2:
        z = 2.0 * P0
        return (z / 2) ** 1.5 * (2 - z * r) * math.exp(-z * r / 2)
    raise ValueError(n)


def genuine_element(n_bra: int, n_ket: int, R: float) -> float:
    """<n_bra s_A | 1/r_B | n_ket s_A> by radial quadrature (Y00^2 = 1/4pi)."""
    f = lambda r: _sturmian_s(n_bra, r) * _sturmian_s(n_ket, r) * r * r / max(r, R)
    a, _ = quad(f, 0.0, R, epsabs=1e-13, limit=200)
    b, _ = quad(f, R, math.inf, epsabs=1e-13, limit=200)
    val = a + b
    return val


def _closed_form(R: float) -> float:
    return 1.0 / R - (1.0 + 1.0 / R) * math.exp(-2.0 * R)


def test_normalization_sanity():
    for n in (1, 2):
        v, _ = quad(lambda r: _sturmian_s(n, r) ** 2 * r * r, 0, math.inf)
        assert abs(v - 1.0) < 1e-10


@pytest.mark.parametrize("R", sorted(PRINTED))
def test_1s1s_genuine_matches_closed_form_and_printed(R):
    v = genuine_element(1, 1, R)
    assert abs(v - _closed_form(R)) < 1e-9, (R, v)
    assert abs(v - PRINTED[R]) < 5e-4, (R, v)


def test_1s1s_is_R_dependent_unlike_model():
    vals = [genuine_element(1, 1, R) for R in sorted(PRINTED)]
    model = 1.0 / P0                       # D^(1) = 1 -> constant in R
    assert all(abs(v - model) > 0.1 for v in vals), vals
    assert all(a > b + 0.1 for a, b in zip(vals, vals[1:])), vals


@pytest.mark.parametrize("R", [0.5, 1.4, 3.015])
def test_1s2s_cross_n_nonzero(R):
    v = genuine_element(1, 2, R)
    assert abs(v) > 1e-2, (R, v)            # block-diagonal D^(n) gives 0
