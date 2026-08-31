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


# ---------------------------------------------------------------------------
# Genuine derivations of the two S^3 Casimir constants (added 2026-08-28 by the
# /qa group6 FULL run).  The tests above pin module constants against literal
# targets; these DERIVE the values from Hurwitz/Bernoulli, so a consistent
# mis-transcription into both module and target can no longer pass unnoticed.
# ---------------------------------------------------------------------------
def test_paper35_scalar_casimir_derived_from_zeta():
    """E_Cas(conformal scalar, unit S^3) = (1/2) zeta_R(-3) = 1/240, derived."""
    import sympy as sp

    zeta_m3 = sp.zeta(-3)
    assert zeta_m3 == sp.Rational(1, 120), f"zeta_R(-3) = {zeta_m3} != 1/120"
    E = sp.Rational(1, 2) * zeta_m3
    assert E == sp.Rational(1, 240), f"derived scalar Casimir {E} != 1/240"
    # rational: carries no transcendental at all (the paper's pi-free claim)
    assert not E.has(sp.pi) and E.is_rational


def test_paper35_dirac_casimir_derived_from_bernoulli():
    """E_Cas(Dirac, unit S^3) = -(1/2) zeta_{|D|}(-1) = +17/480, derived.

    zeta_{|D|}(s) = 4 [ zeta_H(s-2, 3/2) - (1/4) zeta_H(s, 3/2) ], and
    zeta_H(-n, a) = -B_{n+1}(a)/(n+1).
    """
    import sympy as sp

    x = sp.Symbol("x")
    B2 = sp.bernoulli(2, sp.Rational(3, 2))
    B4 = sp.bernoulli(4, sp.Rational(3, 2))
    assert B2 == sp.Rational(11, 12), f"B_2(3/2) = {B2} != 11/12"
    assert B4 == sp.Rational(127, 240), f"B_4(3/2) = {B4} != 127/240"

    # zeta_H(-3, 3/2) = -B_4(3/2)/4 ;  zeta_H(-1, 3/2) = -B_2(3/2)/2
    zh_m3 = -B4 / 4
    zh_m1 = -B2 / 2
    zeta_D_m1 = 4 * (zh_m3 - sp.Rational(1, 4) * zh_m1)
    assert zeta_D_m1 == sp.Rational(-17, 240), f"zeta_|D|(-1) = {zeta_D_m1} != -17/240"

    E = -sp.Rational(1, 2) * zeta_D_m1
    assert E == sp.Rational(17, 480), f"derived Dirac Casimir {E} != +17/480"
    # SIGN is the historically fragile part: it must be POSITIVE
    assert E > 0, "Dirac Casimir sign flipped (the 2026-07-04 regression)"
    assert not E.has(sp.pi) and E.is_rational

    # the module constant must agree with the derivation, not merely with itself
    from geovac.thermal_tensor_triple import dirac_casimir_S3
    got = dirac_casimir_S3()["casimir_energy_unit_S3_full_dirac"]
    assert sp.nsimplify(got) == E, f"module {got} disagrees with derived {E}"

    # discrimination guard: the derivation rejects a mis-transcribed Bernoulli
    bad = 4 * (-sp.Rational(127, 241) / 4 - sp.Rational(1, 4) * zh_m1)
    assert -sp.Rational(1, 2) * bad != E, "derivation insensitive to B_4 corruption"


def test_paper35_kg2_seventh_eigenvalue_is_2pi():
    """Paper 35 obs KG-2: at beta=1, m=0 the first seven eigenvalues of the
    compactified KG spectrum omega^2 = n(n+2) + (2 pi k / beta)^2, sorted
    ascending, are 0, sqrt3, 2sqrt2, sqrt15, 2sqrt6, sqrt35, 2pi -- the
    seventh is the FIRST pi-bearing one.  (Cert-2 coverage gap: the ordinal
    was corrected 2026-08-28 but never sorted/counted by a test.)"""
    import math
    vals = set()
    for n in range(0, 12):
        for k in range(0, 4):
            vals.add(n * (n + 2) + (2 * math.pi * k) ** 2)
    first7 = sorted(vals)[:7]
    # first six are the integers n(n+2), n=0..5 (pi-free by construction)
    assert first7[:6] == [0, 3, 8, 15, 24, 35]
    # seventh is exactly (2 pi)^2 -- the (n=0, k=1) Matsubara mode
    assert abs(first7[6] - 4 * math.pi ** 2) < 1e-12
    # and it beats the next spatial mode n=6 -> 48
    assert first7[6] < 48


def test_paper35_kg1_square_free_generator_census():
    """Paper 35 KG-1 census (recounted 2026-08-28): for m^2 in {0, 1/4, 2}
    and n in [1, 50], the union of square-free radicand generators has 136
    distinct integers, from 2 to 10401, and 1 is NOT a generator.  (Cert-2
    coverage gap: the recount was never pinned.)"""
    import sympy as sp

    def squarefree_part(q):
        # q rational > 0: write sqrt(q) = c*sqrt(d), c in Q, d square-free int
        num, den = sp.fraction(sp.Rational(q))
        # sqrt(num/den) = sqrt(num*den)/den -> radicand = squarefree(num*den)
        m = int(num * den)
        d = 1
        for p, e in sp.factorint(m).items():
            if e % 2 == 1:
                d *= p
        return d

    gens = set()
    for n in range(1, 51):
        for m_sq in (sp.Integer(0), sp.Rational(1, 4), sp.Integer(2)):
            om2 = sp.Integer(n) * (n + 2) + m_sq
            d = squarefree_part(om2)
            if d != 0:
                gens.add(d)
    # perfect squares contribute generator 1 only if omega is rational;
    # the paper's census counts nontrivial generators (1 is not in the set)
    gens.discard(1)
    assert len(gens) == 136, f"census {len(gens)} != 136"
    assert min(gens) == 2 and max(gens) == 10401
