"""Category 4 -- resurgent and connection data in closed form.

Every corpus result reproduced here is already certified by a backing test; this
module re-evaluates the closed forms so the table carries the actual numbers,
and records which test does the certifying.  Nothing is re-derived.

Background in one paragraph.  Several of the corpus's integrals are given not by
a convergent series but by a divergent (asymptotic) one.  Such a series is still
a complete description of the function provided one also knows its "resurgent"
data: where, in the plane of the Borel transform, the series has singularities,
and with what amplitude.  Those amplitudes are the Stokes constants.  The corpus
result is that across four independent objects the Stokes data comes out
ALGEBRAIC (a root of a polynomial with rational coefficients) up to a fixed power
of pi that is determined by the local type of the singularity -- pole gives
2*pi*i, square-root branch gives pi^0, three-halves branch gives 1/sqrt(pi) --
while any genuinely transcendental content is pushed to the boundary of the
domain, where it appears as a single period.  The entries below are the numbers
in that statement.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List

import mpmath as mp

from ._common import entry, frac_str

CLAIM = 50
DPS = 70

#: reference moduli for the one-mass fibre N(D).  rho = 1/2 is the CM fibre
#: reached by the physical three-centre configuration; rho = 1/5 is the
#: corpus's standard off-CM test point.
RHO_REFS = [(Fraction(1, 5), "1/5"), (Fraction(1, 2), "1/2 (the CM fibre)")]

#: reference decay-rate pairs (a, b) for the exchange class.
RATE_PAIRS = [(Fraction(1), Fraction(2)), (Fraction(3, 2), Fraction(5, 2))]


def _psi(rho):
    """Borel transform of the N(D) series: [(z+2)(rho(1+z)^2 + 1 - rho)]^{-1/2}."""
    return lambda z: 1 / mp.sqrt((z + 2) * (rho * (1 + z) ** 2 + (1 - rho)))


def build(mode: str = "full") -> List[Dict[str, Any]]:
    del mode
    rows: List[Dict[str, Any]] = []

    with mp.workdps(DPS):
        # ---- the E_1 seed --------------------------------------------------
        jumps = []
        for x in (mp.mpf("1.5"), mp.mpf("2.5")):
            eps = mp.mpf("1e-40")
            jumps.append(mp.e1(mp.mpc(-x, -eps)) - mp.e1(mp.mpc(-x, eps)))
        worst = max(abs(j - 2 * mp.pi * mp.mpc(0, 1)) for j in jumps)
        rows.append(entry(
            "resurgent.e1_seed.stokes",
            "resurgent",
            "Stokes constant of the exchange seed e^a E_1(a)",
            "2*pi*i  (rational multiple 1)",
            "exact",
            "The Paper 18 Level-2 exchange seed e^a E_1(a) has Borel transform "
            "1/(1 + zeta), a single simple pole at zeta = -1.  Its Stokes "
            "constant is the discontinuity across the cut, "
            "E_1(-x - i0) - E_1(-x + i0) = 2*pi*i, with rational multiple exactly 1. "
            "This is the rank-1 calibration anchor of the four-object pattern: "
            "pole-type singularity gives the pi-power +1.",
            f"EXACT symbolic identity.  Numerically confirmed here at working "
            f"precision {DPS}: the measured jump differs from 2*pi*i by at most "
            f"{mp.nstr(worst, 3)} at x = 1.5 and x = 2.5, the contour-offset "
            f"floor.  Backing test: tests/test_paper59_resurgent_skeleton.py::"
            "test_e1_seed_stokes_constant.",
            value_kind="algebraic",
            transcendence_class="{pi} (pi-power +1, fixed by the pole type)",
            backing_test="tests/test_paper59_resurgent_skeleton.py",
        ))

        # ---- N(D): real-sector Stokes amplitude ---------------------------
        deltas = []
        for rho_f, _ in RHO_REFS:
            rho = mp.mpf(rho_f.numerator) / rho_f.denominator
            a2 = mp.limit(lambda h, p=_psi(rho): p(-2 + h) * mp.sqrt(h), 0)
            deltas.append(abs(a2 - 1))
        rows.append(entry(
            "resurgent.N.stokes_real",
            "resurgent",
            "N(D) real-sector Stokes amplitude a_*(-2)",
            "1",
            "exact",
            "N(D) is the one-mass fibre of the Paper 59 three-centre integral, a "
            "divergent (Gevrey-1) Bessel series whose Borel transform is the "
            "algebraic function psi(z) = [(z+2)(rho(1+z)^2 + 1 - rho)]^{-1/2}. "
            "The amplitude of its real branch point at z = -2 is the limit of "
            "psi(-2 + h) sqrt(h) as h -> 0, which is exactly 1 -- independent of "
            "the modulus rho.  Square-root branch point, so the pi-power is 0.",
            f"EXACT algebraic value.  Recomputed here at working precision {DPS} "
            f"for both reference moduli rho = 1/5 and rho = 1/2: worst deviation "
            f"from 1 is {mp.nstr(max(deltas), 3)}.  Backing test: "
            "tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.",
            value_kind="algebraic",
            transcendence_class="{} (algebraic over Q(rho); pi-power 0)",
            backing_test="tests/test_routeC_momentum.py",
        ))

        # ---- N(D): complex-sector Stokes amplitude squared ----------------
        for rho_f, rho_label in RHO_REFS:
            rho = mp.mpf(rho_f.numerator) / rho_f.denominator
            omega = mp.sqrt((1 - rho) / rho)
            ac = mp.limit(lambda h, p=_psi(rho), w=omega:
                          p((-1 + 1j * w) + h) * mp.sqrt(h), 0)
            target = mp.mpf(-1) / 2 - (1j / 2) * mp.sqrt(rho / (1 - rho))
            delta = abs(ac ** 2 - target)
            # exact algebraic form of sqrt(rho/(1-rho)) at these moduli
            ratio = Fraction(rho_f, 1 - rho_f)
            sq = mp.sqrt(mp.mpf(ratio.numerator) / ratio.denominator)
            exact_txt = ("-1/2 - (1/4) i" if rho_f == Fraction(1, 5)
                         else "-1/2 - (1/2) i")
            rows.append(entry(
                f"resurgent.N.stokes_complex.rho{rho_f.numerator}_{rho_f.denominator}",
                "resurgent",
                f"N(D) complex-sector Stokes amplitude squared, "
                f"a_*(-1 + i omega)^2, at rho = {rho_label}",
                exact_txt,
                "exact",
                "The same Borel transform psi has a conjugate pair of branch "
                "points at z = -1 +- i omega with omega = sqrt((1-rho)/rho). "
                "The square of the amplitude there is "
                "a_*(-1 + i omega)^2 = -1/2 - (i/2) sqrt(rho/(1-rho)), algebraic "
                f"over Q(rho); at rho = {rho_label} the square root is rational, "
                f"omega = {mp.nstr(omega, 3)}, so the value is a Gaussian rational.",
                f"EXACT algebraic value.  Recomputed here at working precision "
                f"{DPS}: the limit of psi(z) sqrt(z - z_*) squared differs from "
                f"the closed form by {mp.nstr(delta, 3)}.  Backing test: "
                "tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.",
                value_kind="algebraic",
                omega=mp.nstr(omega, 20),
                transcendence_class="{} (algebraic over Q(rho); pi-power 0)",
                backing_test="tests/test_routeC_momentum.py",
            ))

        # ---- N(D): the boundary period ------------------------------------
        rho = mp.mpf(1) / 5
        direct = mp.quad(lambda x: 1 / mp.sqrt((x ** 2 - 1)
                                               * (rho * x ** 2 + (1 - rho))),
                         [1, mp.inf])
        Kval = mp.ellipk(1 - rho)
        rows.append(entry(
            "resurgent.N.boundary_period.rho1_5",
            "resurgent",
            "N(D) boundary period at rho = 1/5: sqrt(c1) N(0) = K(4/5)",
            mp.nstr(Kval, CLAIM, strip_zeros=False),
            CLAIM,
            "All the transcendental content of N(D) sits at the boundary D = 0, "
            "where the Laplace integral degenerates to the complete elliptic "
            "integral of the first kind: sqrt(c1) N(0) = K(1 - rho).  Evaluated "
            "here as mpmath's ellipk with parameter m = 1 - rho = 4/5.",
            f"Two independent routes at working precision {DPS}: mpmath's "
            f"arithmetic-geometric-mean ellipk, and a direct adaptive quadrature "
            f"of the period integral int_1^inf dx / sqrt((x^2-1)(rho x^2 + 1-rho)), "
            f"agreeing to {mp.nstr(abs(direct - Kval), 3)} absolute (the "
            f"quadrature's own endpoint-singularity floor).  Claiming {CLAIM} "
            f"digits from the AGM route, which is the accurate one.  Backing "
            "test: tests/test_routeC_momentum.py::"
            "test_N_stokes_constants_algebraic.",
            value_kind="decimal",
            transcendence_class="elliptic period (weight 1, genus 1)",
            defining_relation="K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta), m = 4/5",
            backing_test="tests/test_routeC_momentum.py",
        ))

        # ---- exchange class: Borel positions, charges, kappa ---------------
        for a_f, b_f in RATE_PAIRS:
            a = mp.mpf(a_f.numerator) / a_f.denominator
            b = mp.mpf(b_f.numerator) / b_f.denominator
            A = a + b
            kappa_f = Fraction(2 * a_f * b_f, a_f + b_f)
            kap = mp.mpf(kappa_f.numerator) / kappa_f.denominator

            sector_res = []
            bnd_res = []
            for R in (mp.mpf("1.3"), mp.mpf("2.9")):
                G = (mp.exp(2 * a * R) * mp.e1(2 * a * R)
                     + mp.exp(2 * b * R) * mp.e1(2 * b * R)
                     - mp.exp(2 * A * R) * mp.e1(2 * A * R))
                lap = mp.quad(lambda s: mp.e ** (-R * s)
                              * (1 / (s + 2 * a) + 1 / (s + 2 * b)
                                 - 1 / (s + 2 * A)), [0, 1, mp.inf])
                sector_res.append(abs(G - lap))
                bnd = mp.euler + mp.log(kap * R)
                lbnd = -R * mp.quad(lambda s: mp.e ** (-R * s) * mp.log(s / kap),
                                    [0, 1, mp.inf])
                bnd_res.append(abs(bnd - lbnd))

            a_txt, b_txt = frac_str(a_f), frac_str(b_f)
            rows.append(entry(
                f"resurgent.exchange.borel.a{a_f.numerator}_{a_f.denominator}"
                f".b{b_f.numerator}_{b_f.denominator}",
                "resurgent",
                f"exchange-class Borel data at decay rates a = {a_txt}, b = {b_txt}",
                f"positions {{0, -{frac_str(2 * a_f)}, -{frac_str(2 * b_f)}, "
                f"-{frac_str(2 * (a_f + b_f))}}}; charges (-1, +1, +1, -1)",
                "exact",
                "The exchange-class closed form has the exact reduced normal form "
                "F = [e^{-AR}(gamma + ln(kappa R)) + e^{(a-b)R} E_1(2aR) "
                "+ e^{(b-a)R} E_1(2bR) - e^{AR} E_1(2AR)] / (2 a b R^2), with "
                "A = a + b and kappa = 2ab/A.  Reading it as a Laplace transform "
                "in R exposes four Borel singularities.  The three E_1 sectors "
                "give simple poles at s = -2a, -2b, -2A with residues +1, +1, -1 "
                "(each e^{2xR} E_1(2xR) = int_0^inf e^{-Rs}/(s + 2x) ds).  The "
                "boundary bundle gamma + ln(kappa R) gives a LOGARITHMIC branch "
                "point at s = 0 with density -ln(s/kappa), i.e. charge -1.  So "
                "the positions are {0, -2a, -2b, -2A} and the integer charges "
                "are (-1, +1, +1, -1).",
                f"EXACT integer data.  The two Laplace representations behind it "
                f"are verified here at working precision {DPS} and two values of "
                f"R: the three-pole density reproduces the E_1 sector sum to "
                f"{mp.nstr(max(sector_res), 3)}, and the logarithmic density "
                f"reproduces gamma + ln(kappa R) to {mp.nstr(max(bnd_res), 3)} "
                f"(both at the quadrature floor).  The normal form itself is an "
                f"exact symbolic identity against "
                "geovac.two_center_eri.ordered_xi_closed, pinned by "
                "tests/test_paper59_resurgent_skeleton.py::"
                "test_exchange_class_normal_form.  Its genuine multivaluedness in "
                "R (unlike the cut-free hybrid class) is pinned by "
                "::test_exchange_class_is_multivalued_in_R.",
                value_kind="integer_tuple",
                transcendence_class="{} (integer charges; algebraic positions)",
                backing_test="tests/test_paper59_resurgent_skeleton.py",
            ))

            rows.append(entry(
                f"resurgent.exchange.kappa.a{a_f.numerator}_{a_f.denominator}"
                f".b{b_f.numerator}_{b_f.denominator}",
                "resurgent",
                f"exchange-class log scale kappa = 2ab/(a+b) at a = {a_txt}, "
                f"b = {b_txt}",
                frac_str(kappa_f),
                "exact",
                "kappa = 2ab/A is the charge-weighted product of the Borel "
                "positions -- the exponents ARE the Stokes charges -- and it is "
                "the argument of the only logarithm in the exchange closed form. "
                "It is a rational function of the decay rates, so at rational "
                "rates it is a rational number.  Euler's gamma appears only "
                "bundled with it, as gamma + ln(kappa R), which is why gamma is "
                "coordinate bookkeeping here and not a Borel-plane "
                "transcendental.",
                "EXACT rational.  Fixed by the exact symbolic normal-form "
                "identity, backing test tests/test_paper59_resurgent_skeleton.py::"
                "test_exchange_class_normal_form (which fails for any other "
                f"kappa).  Decimal form {mp.nstr(kap, 20)}.",
                value_kind="rational",
                decimal_value=mp.nstr(kap, 20),
                transcendence_class="{} (rational)",
                backing_test="tests/test_paper59_resurgent_skeleton.py",
            ))

    return rows
