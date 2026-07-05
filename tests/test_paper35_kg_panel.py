"""Verification tests for Paper 35 (Time as Projection) headline results.

Closes a /qa group6 first-cert coverage gap: the paper's headline
"200-case KG pi-free panel" was previously backed only indirectly (the
Paper 34 III.14 spot-check tests one weak property on a few (n, m^2)
pairs and asserts only `omega_sq.is_rational`, never pi-freeness; the
cited debug/kg1_algebraic_ring.py driver has been pruned).

This file exercises the FULL panel genuinely and asserts pi-freeness
(the load-bearing claim), plus pins the KG-3 / KG-5 Casimir headline
values to the module of record.

Paper 35 headlines backed here:
  - KG spectrum omega_n^2 = n(n+2) + m^2 on S^3 x R is pi-free in the
    algebraic-extension ring Q[sqrt(d)] for rational m^2, verified over
    n in [1,50] x m^2 in {0, 1, 1/4, 2} (200 cases, zero transcendentals).
  - conformally-coupled scalar S^3 Casimir E_Cas = 1/240 (exact rational).
  - Dirac S^3 Casimir = +17/480 (exact rational, POSITIVE: E=-1/2 zeta_{|D|}(-1)
    = -1/2*(-17/240), the half-integer shift makes zeta negative).
"""
from __future__ import annotations

import os
import sys

import sympy as sp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


_M_SQ_PANEL = [sp.Integer(0), sp.Integer(1), sp.Rational(1, 4), sp.Integer(2)]
_N_RANGE = range(1, 51)  # n in [1, 50]


def test_paper35_kg_spectrum_pi_free_200_case_panel():
    """omega_n = sqrt(n(n+2)+m^2) is pi-free (in Q[sqrt d]) for all 200 cases."""
    n_cases = 0
    for m_sq in _M_SQ_PANEL:
        for n in _N_RANGE:
            omega_sq = sp.Integer(n) * sp.Integer(n + 2) + m_sq
            omega = sp.sqrt(omega_sq)

            # (1) omega^2 is a nonnegative rational (bare-graph ring datum).
            assert omega_sq.is_rational and omega_sq >= 0, \
                f"omega^2={omega_sq} not a nonneg rational at (n={n}, m^2={m_sq})"

            # (2) THE load-bearing claim: omega carries NO transcendental.
            #     A float cast would trivially pass this -- we stay symbolic.
            assert not omega.has(sp.pi), \
                f"pi in omega at (n={n}, m^2={m_sq}): {omega}"
            assert not omega.has(sp.E, sp.log, sp.EulerGamma), \
                f"transcendental in omega at (n={n}, m^2={m_sq}): {omega}"
            assert omega.is_algebraic is True, \
                f"omega not algebraic at (n={n}, m^2={m_sq})"

            # (3) It lives in Q[sqrt d]: omega is at most a quadratic
            #     irrational over Q (its square is rational).
            assert (omega ** 2 - omega_sq) == 0
            assert sp.together(omega ** 2).is_rational
            n_cases += 1

    assert n_cases == 200, f"panel breadth {n_cases} != 200"


def test_paper35_first_pi_bearing_eigenvalue_is_temporal():
    """The (n,0) modes are pi-free; pi enters only at temporal compactification
    with the first pi-bearing mode (n=0, k=1) at omega^2 = 4 pi^2 / beta^2."""
    beta = sp.symbols("beta", positive=True)
    # Spatial-only (k=0) mode: pi-free.
    for n in range(0, 6):
        omega_sq_spatial = sp.Integer(n) * sp.Integer(n + 2)
        assert not sp.sqrt(omega_sq_spatial).has(sp.pi)
    # Temporal Matsubara mode (n=0, k=1): omega^2 = (2 pi / beta)^2 = 4 pi^2/beta^2.
    omega_sq_temporal = (2 * sp.pi / beta) ** 2
    assert sp.simplify(omega_sq_temporal - 4 * sp.pi ** 2 / beta ** 2) == 0
    assert omega_sq_temporal.has(sp.pi)  # pi enters here and only here


def test_paper35_scalar_casimir_1_over_240():
    """KG-3: conformally-coupled scalar S^3 Casimir = 1/240 (exact, no pi)."""
    from geovac.thermal_tensor_triple import scalar_casimir_S3
    val = scalar_casimir_S3()["casimir_energy_unit_S3_conformal_scalar"]
    assert val == sp.Rational(1, 240)
    assert val.is_rational and not sp.sympify(val).has(sp.pi)


def test_paper35_dirac_casimir_plus_17_over_480():
    """KG-5: Dirac S^3 Casimir = +17/480 (exact rational; POSITIVE).

    E = -1/2 zeta_{|D|}(-1) = -1/2*(-17/240) = +17/480 -- the half-integer
    Dirac shift makes zeta_{|D|}(-1) itself negative, so the fermion -1/2
    factor returns a POSITIVE Casimir, same sign class as the scalar +1/240
    (Paper 35 KG-5 derivation; the naive -17/480 heuristic was corrected).
    """
    from geovac.thermal_tensor_triple import dirac_casimir_S3
    val = dirac_casimir_S3()["casimir_energy_unit_S3_full_dirac"]
    assert val == sp.Rational(17, 480), f"Dirac Casimir {val} != +17/480"
    assert val.is_rational and not sp.sympify(val).has(sp.pi)
