"""Category 2 -- two-centre hydrogenic electron-repulsion integrals (Paper 58).

Four closed-form classes from ``geovac/two_center_eri.py``, each evaluated on a
canonical grid and each cross-validated against an INDEPENDENT numerical
quadrature that shares no code with the closed form.

  class            quartet shape                       transcendence seeds
  ---------------  ----------------------------------  --------------------
  aabb             (a b | c d), ab on A, cd on B        {exp}
  hybrid (l=0)     (a b | c d), abc on A, d on B,       {exp}
                   a,b both s-type
  hybrid (l>0)     same, a,b with l > 0                 {exp, E_1, ln}
  exchange         ordered two-centre xi integral       {exp, E_1, ln, gamma}

The closed forms are exact symbolic expressions, so their digit count is limited
only by how far one chooses to evaluate them; the honest limit on what may be
CLAIMED as a verified integral is set by the independent quadrature, a float64
route good to roughly 1e-10..1e-15.  This module therefore reports two separate
things and never conflates them:

  * ``digits_claimed`` = 50, established by evaluating the SAME exact expression
    at working precisions 60 and 90 and requiring the first 50 significant
    digits to be identical.  This certifies that the printed digits are the
    digits of the closed-form expression -- it says nothing on its own about
    whether the closed form is the right integral.
  * ``independent_route_rel_agreement`` = the measured agreement with the
    quadrature route, which is what certifies that the closed form IS the
    integral, at that route's own achievable precision.

Both are reported for every entry.  Neither is silently upgraded into the other.

DOMAIN NOTE (measured while building this table, 2026-08-21): for l_a, l_b > 0
the hybrid closed form goes through the shell route, which requires
Z_B < Z_A strictly.  At Z_B = Z_A it returns NaN (a removable rate coincidence:
approaching Z_B -> Z_A from below converges to the quadrature value), and at
Z_B > Z_A it raises ``AssertionError: Ei reached step 3`` (the exponential-
integral argument turns negative and that branch is not implemented).  The
l > 0 hybrid grid therefore uses (Z_A, Z_B) in {(3,1), (4,2)} rather than
{(1,1), (3,1)}.
"""
from __future__ import annotations

import warnings
from fractions import Fraction
from typing import Any, Dict, List, Tuple

import mpmath as mp
import numpy as np
import sympy as sp
from scipy import integrate

from geovac.two_center_eri import (aabb_closed_form, aabb_quadrature,
                                   hybrid_closed_form, hybrid_quadrature,
                                   ordered_xi_closed)

from ._common import agree_digits, entry

# Canonical grid -------------------------------------------------------------
# R = 3.015 bohr is the balanced-solver LiH equilibrium separation; 1.4 bohr is
# near the H2 equilibrium; 2.0 bohr is a round intermediate point.
R_GRID = [("1.4", sp.Rational(7, 5)),
          ("2.0", sp.Integer(2)),
          ("3.015", sp.Rational(603, 200))]
CHARGE_GRID = [(1, 1), (3, 1)]
CHARGE_GRID_LGT0 = [(3, 1), (4, 2)]           # see DOMAIN NOTE above

CLAIM = 50            # significant digits claimed for a closed-form value
DPS_LO, DPS_HI = 60, 90


def _two_precision_value(expr: sp.Expr) -> Tuple[str, int]:
    """Evaluate an exact expression at two working precisions.

    Returns the value printed to CLAIM significant digits together with the
    number of leading significant digits on which the two precisions agree.
    """
    lo = sp.N(sp.re(expr), DPS_LO)
    hi = sp.N(sp.re(expr), DPS_HI)
    with mp.workdps(DPS_HI + 10):
        agree = agree_digits(mp.mpf(str(lo)), mp.mpf(str(hi)))
    with mp.workdps(CLAIM + 5):
        printed = mp.nstr(mp.mpf(str(sp.N(sp.re(expr), CLAIM + 5))),
                          CLAIM, strip_zeros=False)
    return printed, agree


def _seed_tag(expr: sp.Expr) -> str:
    names = {type(f).__name__ for f in expr.atoms(sp.Function)}
    seeds = []
    if "exp" in names:
        seeds.append("exp")
    if "expint" in names:
        seeds.append("E_1")
    if "log" in names:
        seeds.append("ln")
    if expr.has(sp.EulerGamma):
        seeds.append("gamma")
    return "{" + ", ".join(seeds) + "}"


def _exchange_quadrature(p1: float, p2: float) -> float:
    """Nested quadrature of the ordered-xi exchange object.

        int_1^oo int_1^oo e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2

    with P_0 = 1 and Q_0(x) = (1/2) ln((x+1)/(x-1)).  Shares no code with
    ``ordered_xi_closed``: it never forms an E_1 moment, a log moment or a
    Legendre-Q antiderivative.
    """
    def outer(x1: float) -> float:
        lo, _ = integrate.quad(lambda x2: np.exp(-p2 * x2), 1.0, x1,
                               epsabs=1e-14, epsrel=1e-12, limit=300)
        hi, _ = integrate.quad(
            lambda x2: np.exp(-p2 * x2) * 0.5 * np.log((x2 + 1) / (x2 - 1)),
            x1, np.inf, epsabs=1e-14, epsrel=1e-12, limit=300)
        q0 = 0.5 * np.log((x1 + 1) / (x1 - 1))
        return np.exp(-p1 * x1) * (q0 * lo + hi)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = integrate.quad(outer, 1.0, np.inf, epsabs=1e-13, epsrel=1e-11,
                                limit=300)
    return val


def _independent_note(closed_50: str, ref: float, route: str) -> Tuple[str, float]:
    got = float(closed_50)
    rel = abs(got - ref) / max(abs(ref), 1e-300)
    return (f"independent route ({route}) gives {ref:.16g}; "
            f"relative agreement {rel:.1e}"), rel


# --------------------------------------------------------------------------
def build(mode: str = "full") -> List[Dict[str, Any]]:
    fast = (mode == "fast")
    rows: List[Dict[str, Any]] = []

    # ---- class 1: (AA|BB), s-type ----------------------------------------
    for ZA, ZB in CHARGE_GRID:
        for rlabel, R in R_GRID:
            e = aabb_closed_form(Fraction(ZA), (1, 0, 0), (1, 0, 0),
                                 Fraction(ZB), (1, 0, 0), (1, 0, 0), R)
            val, agree = _two_precision_value(e)
            if fast:
                indep, rel = "(skipped in fast mode)", None
            else:
                ref = aabb_quadrature(Fraction(ZA), (1, 0, 0), (1, 0, 0),
                                      Fraction(ZB), (1, 0, 0), (1, 0, 0),
                                      float(R))
                indep, rel = _independent_note(
                    val, ref, "direct numerical quadrature of the same quartet")
            rows.append(entry(
                f"eri.aabb.Z{ZA}{ZB}.R{rlabel}",
                "two_center_eri",
                f"(1s_A 1s_A | 1s_B 1s_B), Z_A={ZA}, Z_B={ZB}, R={rlabel} bohr",
                val, CLAIM,
                "Exact closed form, geovac.two_center_eri.aabb_closed_form: both "
                "charge distributions sit on one centre each, so the quartet "
                "reduces to a one-electron two-centre problem via the exact (L,M) "
                "multipole decomposition of each density together with the "
                "closed-form radial potential V_L.  No expansion of 1/r12 and no "
                "quadrature anywhere in the route.",
                f"Exact symbolic expression evaluated at working precisions "
                f"{DPS_LO} and {DPS_HI}: first {min(agree, CLAIM)} significant "
                f"digits identical (claiming {CLAIM}).  Independent-route check: "
                f"{indep}.  Backing tests: "
                "tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, "
                "::test_closed_form_centre_swap_consistency.",
                transcendence_class=_seed_tag(e),
                independent_route_rel_agreement=rel,
                two_precision_agree_digits=min(agree, 10 ** 6),
                backing_test="tests/test_two_center_eri_aabb.py",
                geovac_entry_point="geovac.two_center_eri.aabb_closed_form",
            ))

    # ---- class 1b: (AA|BB) with l > 0 on centre A ------------------------
    for ZA, ZB in CHARGE_GRID:
        rlabel, R = R_GRID[1]                       # R = 2.0 only
        e = aabb_closed_form(Fraction(ZA), (2, 1, 0), (2, 1, 0),
                             Fraction(ZB), (1, 0, 0), (1, 0, 0), R)
        val, agree = _two_precision_value(e)
        if fast:
            indep, rel = "(skipped in fast mode)", None
        else:
            ref = aabb_quadrature(Fraction(ZA), (2, 1, 0), (2, 1, 0),
                                  Fraction(ZB), (1, 0, 0), (1, 0, 0), float(R))
            indep, rel = _independent_note(
                val, ref, "direct numerical quadrature of the same quartet")
        rows.append(entry(
            f"eri.aabb_p.Z{ZA}{ZB}.R{rlabel}",
            "two_center_eri",
            f"(2p0_A 2p0_A | 1s_B 1s_B), Z_A={ZA}, Z_B={ZB}, R={rlabel} bohr",
            val, CLAIM,
            "Exact closed form, geovac.two_center_eri.aabb_closed_form, with "
            "l = 1 on the one-centre pair.  The (AA|BB) class stays elementary at "
            "any angular momentum: l enters only through Gaunt coefficients, "
            "which are rational multiples of square roots.",
            f"Exact symbolic expression at working precisions {DPS_LO} and "
            f"{DPS_HI}: first {min(agree, CLAIM)} significant digits identical "
            f"(claiming {CLAIM}).  Independent-route check: {indep}.  Backing "
            "test: tests/test_two_center_eri_aabb.py::"
            "test_closed_form_matches_quadrature_for_l_gt_0_and_M_ne_0.",
            transcendence_class=_seed_tag(e),
            independent_route_rel_agreement=rel,
            two_precision_agree_digits=min(agree, 10 ** 6),
            backing_test="tests/test_two_center_eri_aabb.py",
            geovac_entry_point="geovac.two_center_eri.aabb_closed_form",
        ))

    # ---- class 2: hybrid, s-type one-centre pair (elementary) ------------
    for ZA, ZB in CHARGE_GRID:
        for rlabel, R in R_GRID:
            e = hybrid_closed_form(Fraction(ZA), (1, 0, 0), (1, 0, 0), (1, 0, 0),
                                   Fraction(ZB), (1, 0, 0), R)
            val, agree = _two_precision_value(e)
            if fast:
                indep, rel = "(skipped in fast mode)", None
            else:
                ref = hybrid_quadrature(Fraction(ZA), (1, 0, 0), (1, 0, 0),
                                        (1, 0, 0), Fraction(ZB), (1, 0, 0),
                                        float(R))
                indep, rel = _independent_note(
                    val, ref, "direct numerical quadrature of the same quartet")
            rows.append(entry(
                f"eri.hybrid_s.Z{ZA}{ZB}.R{rlabel}",
                "two_center_eri",
                f"(1s_A 1s_A | 1s_A 1s_B), Z_A={ZA}, Z_B={ZB}, R={rlabel} bohr",
                val, CLAIM,
                "Exact closed form, geovac.two_center_eri.hybrid_closed_form "
                "(direct V_L route).  Three orbitals sit on centre A and one on "
                "centre B.  With an s-type one-centre pair every r_A power is "
                "non-negative, so no exponential-integral seed appears and the "
                "class is elementary.",
                f"Exact symbolic expression at working precisions {DPS_LO} and "
                f"{DPS_HI}: first {min(agree, CLAIM)} significant digits identical "
                f"(claiming {CLAIM}).  Independent-route check: {indep}.  Backing "
                "tests: tests/test_two_center_eri_aabb.py::"
                "test_hybrid_closed_form_matches_quadrature, "
                "::test_hybrid_s_type_is_elementary, "
                "::test_hybrid_two_routes_agree_on_the_s_type_overlap.",
                transcendence_class=_seed_tag(e),
                independent_route_rel_agreement=rel,
                two_precision_agree_digits=min(agree, 10 ** 6),
                backing_test="tests/test_two_center_eri_aabb.py",
                geovac_entry_point="geovac.two_center_eri.hybrid_closed_form",
            ))

    # ---- class 3: hybrid, l > 0 one-centre pair ({E_1, ln}) --------------
    for ZA, ZB in CHARGE_GRID_LGT0:
        for rlabel, R in R_GRID:
            e = hybrid_closed_form(Fraction(ZA), (2, 1, 0), (2, 1, 0), (1, 0, 0),
                                   Fraction(ZB), (1, 0, 0), R)
            val, agree = _two_precision_value(e)
            if fast:
                indep, rel = "(skipped in fast mode)", None
            else:
                ref = hybrid_quadrature(Fraction(ZA), (2, 1, 0), (2, 1, 0),
                                        (1, 0, 0), Fraction(ZB), (1, 0, 0),
                                        float(R))
                indep, rel = _independent_note(
                    val, ref, "direct numerical quadrature of the same quartet")
            rows.append(entry(
                f"eri.hybrid_p.Z{ZA}{ZB}.R{rlabel}",
                "two_center_eri",
                f"(2p0_A 2p0_A | 1s_A 1s_B), Z_A={ZA}, Z_B={ZB}, R={rlabel} bohr",
                val, CLAIM,
                "Exact closed form, geovac.two_center_eri.hybrid_closed_form "
                "(shell route).  With l > 0 on the one-centre pair the minimum "
                "r_A power is -2(l_a + l_b), so the radial integral acquires the "
                "exponential-integral seed E_1; the lower endpoint |r_B - R| "
                "passes through zero and contributes a logarithm whose argument "
                "is a ratio of decay rates, and is therefore R-independent.",
                f"Exact symbolic expression at working precisions {DPS_LO} and "
                f"{DPS_HI}: first {min(agree, CLAIM)} significant digits identical "
                f"(claiming {CLAIM}).  Independent-route check: {indep}.  Backing "
                "tests: tests/test_two_center_eri_aabb.py::"
                "test_hybrid_l_gt_0_via_shells, "
                "::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain "
                "restriction: this route requires Z_B < Z_A strictly (see the "
                "module docstring of benchmarks/certified_reference/"
                "entries_two_center.py).",
                transcendence_class=_seed_tag(e),
                independent_route_rel_agreement=rel,
                two_precision_agree_digits=min(agree, 10 ** 6),
                backing_test="tests/test_two_center_eri_aabb.py",
                geovac_entry_point="geovac.two_center_eri.hybrid_closed_form",
            ))

    # ---- class 4: exchange, ordered-xi closed form -----------------------
    for ZA, ZB in CHARGE_GRID:
        for rlabel, R in R_GRID:
            p1, p2 = sp.Integer(ZA) * R, sp.Integer(ZB) * R
            e = ordered_xi_closed(p1, p2)
            val, agree = _two_precision_value(e)
            if fast:
                indep, rel = "(skipped in fast mode)", None
            else:
                ref = _exchange_quadrature(float(p1), float(p2))
                indep, rel = _independent_note(
                    val, ref, "nested adaptive quadrature of the same double "
                              "integral, forming no E_1 or log moments")
            rows.append(entry(
                f"eri.exchange.Z{ZA}{ZB}.R{rlabel}",
                "two_center_eri",
                f"ordered-xi exchange kernel, rates p1 = {ZA}R, p2 = {ZB}R, "
                f"R = {rlabel} bohr",
                val, CLAIM,
                "Exact closed form, geovac.two_center_eri.ordered_xi_closed: "
                "int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, "
                "the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated "
                "at p_i = Z_i R.  This is an iterated integral over a simplex -- "
                "the shape that defines a period -- and it closes at WEIGHT ONE: "
                "the seeds are exp, E_1, log and Euler's gamma, with no "
                "dilogarithm.",
                f"Exact symbolic expression at working precisions {DPS_LO} and "
                f"{DPS_HI}: first {min(agree, CLAIM)} significant digits identical "
                f"(claiming {CLAIM}).  Independent-route check: {indep}.  Backing "
                "tests: tests/test_two_center_eri_aabb.py::"
                "test_ordered_xi_closed_matches_quadrature, "
                "::test_ordered_xi_is_weight_one; "
                "tests/test_paper59_resurgent_skeleton.py::"
                "test_exchange_class_normal_form.",
                transcendence_class=_seed_tag(e),
                independent_route_rel_agreement=rel,
                two_precision_agree_digits=min(agree, 10 ** 6),
                backing_test="tests/test_two_center_eri_aabb.py",
                geovac_entry_point="geovac.two_center_eri.ordered_xi_closed",
            ))

    return rows
