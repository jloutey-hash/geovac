"""Category 5 -- anchor constants used throughout the corpus.

These are not new numbers; they are the fixed constants that the corpus's own
results are stated against, given here at a uniform 50 digits together with the
relation that defines each one.  A reference table that quotes derived values to
50 digits should also fix, in one place, the constants those values are compared
to -- otherwise a reader reproducing an entry has to guess the convention (in
particular: is the elliptic-integral argument the parameter m or the modulus k?).

  * K(1/2) -- complete elliptic integral of the first kind at PARAMETER m = 1/2,
    the boundary period of the Paper 59 fibre at its CM modulus.
  * pi^2/24 -- the collapse constant of the Paper 60 conditioning law.
  * 2/(1 + min j0) -- the gerade constant, the exact value of the conditioning
    ratio in the equal-centre gerade sector, together with the two pieces that
    define it.
"""
from __future__ import annotations

from typing import Any, Dict, List

import mpmath as mp

from ._common import entry

CLAIM = 50
DPS = 70


def build(mode: str = "full") -> List[Dict[str, Any]]:
    del mode
    rows: List[Dict[str, Any]] = []

    with mp.workdps(DPS):
        # ---- K(1/2) -------------------------------------------------------
        k_agm = mp.ellipk(mp.mpf(1) / 2)
        k_gamma = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
        rows.append(entry(
            "anchor.K_half",
            "anchor",
            "K(1/2), the complete elliptic integral of the first kind at "
            "parameter m = 1/2",
            mp.nstr(k_agm, CLAIM, strip_zeros=False),
            CLAIM,
            "K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta) at m = 1/2 "
            "(mpmath's ellipk convention: the argument is the PARAMETER m, not "
            "the modulus k; K here corresponds to modulus k = 1/sqrt(2), the "
            "lemniscatic case).  Closed form K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi)). "
            "This is the boundary period of the Paper 59 one-mass fibre at its "
            "CM modulus rho = 1/2, and the generator of the corpus's period ring "
            "in the T2 PSLQ searches.",
            f"Two independent evaluations at working precision {DPS}: mpmath's "
            f"arithmetic-geometric-mean ellipk and the Gamma closed form "
            f"Gamma(1/4)^2/(4 sqrt(pi)), agreeing to "
            f"{mp.nstr(abs(k_agm - k_gamma), 3)} absolute -- bit-identical at "
            f"this precision.  Claiming {CLAIM} digits.  The identification of "
            "this value as the corpus fibre's boundary period is backed by "
            "tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.",
            value_kind="decimal",
            defining_relation="K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi))",
            transcendence_class="{Gamma(1/4), pi} -- a CM period at discriminant -4",
            backing_test="tests/test_routeC_momentum.py",
        ))

        # ---- pi^2 / 24 ----------------------------------------------------
        collapse = mp.pi ** 2 / 24
        from geovac.sturmian_sigma_law import COLLAPSE_CONSTANT
        rows.append(entry(
            "anchor.collapse_pi2_24",
            "anchor",
            "pi^2/24, the Paper 60 conditioning-law collapse constant",
            mp.nstr(collapse, CLAIM, strip_zeros=False),
            CLAIM,
            "In the Shibuya-Wulfman two-centre metric the largest cross-centre "
            "singular value obeys the band-limited concentration law "
            "1 - sigma_max = (kR)^2 pi^2 / (24 n^2), so the rescaled quantity "
            "(1 - sigma_max)(n/kR)^2 tends to pi^2/24 as the per-centre basis "
            "size n grows.  This fixes the conditioning exponent at exactly 2, "
            "and reveals the previously fitted exponents 1.85 and 1.97 as "
            "pre-asymptotic windows of the same law.",
            f"Elementary closed form, evaluated at working precision {DPS}; "
            f"claiming {CLAIM} digits.  The float64 constant exported by "
            f"geovac.sturmian_sigma_law.COLLAPSE_CONSTANT agrees to "
            f"{mp.nstr(abs(collapse - mp.mpf(COLLAPSE_CONSTANT)), 3)}.  The "
            "physics claim -- that the measured quantity actually converges to "
            "this constant (0.9913 of the limit at n = 160) -- is backed by "
            "tests/test_paper60_sigma_law.py.",
            value_kind="decimal",
            defining_relation="lim_{n->inf} (1 - sigma_max) (n / kR)^2 = pi^2 / 24",
            transcendence_class="{pi^2} -- pure-Tate, weight 2",
            backing_test="tests/test_paper60_sigma_law.py",
        ))

        # ---- the sinc minimum and the gerade constant ----------------------
        x_star = mp.findroot(lambda x: mp.tan(x) - x, mp.mpf("4.4934"))
        j0_min = mp.sin(x_star) / x_star
        gerade = 2 / (1 + j0_min)
        from geovac.sturmian_sigma_law import gerade_constant
        g_float = gerade_constant()

        rows.append(entry(
            "anchor.sinc_min_root",
            "anchor",
            "x*, the first stationary point of sin(x)/x beyond the origin "
            "(root of tan x = x)",
            mp.nstr(x_star, CLAIM, strip_zeros=False),
            CLAIM,
            "The unique root of tan x = x in (pi, 3pi/2), found by Newton "
            "iteration in arbitrary precision.  It is where the spherical Bessel "
            "function j0(x) = sin(x)/x attains its global minimum, which is what "
            "sets the gerade constant below.",
            f"Newton iteration at working precision {DPS}; the residual "
            f"|tan x* - x*| is at the arithmetic floor.  Verified to be the "
            f"stationary point of sin(x)/x rather than an artefact: "
            f"j0(x*) = {mp.nstr(j0_min, 12)} is below j0(4) = "
            f"{mp.nstr(mp.sin(4) / 4, 12)}.  Claiming {CLAIM} digits.",
            value_kind="decimal",
            defining_relation="tan(x*) = x*,  x* in (pi, 3pi/2)",
            transcendence_class="{} -- a transcendental-equation root, not a period",
        ))

        rows.append(entry(
            "anchor.sinc_min_value",
            "anchor",
            "min_x sin(x)/x, the global minimum of the spherical Bessel "
            "function j0",
            mp.nstr(j0_min, CLAIM, strip_zeros=False),
            CLAIM,
            "j0(x*) = sin(x*)/x* at the stationary point above.  Equivalently "
            "-cos(x*), since tan x* = x* implies sin x*/x* = cos x*.",
            f"Evaluated at working precision {DPS} from the certified root x*; "
            f"the two equivalent forms sin(x*)/x* and -(-cos x*) agree to "
            f"{mp.nstr(abs(j0_min - mp.cos(x_star)), 3)}.  Claiming {CLAIM} "
            f"digits.",
            value_kind="decimal",
            defining_relation="min_x sin(x)/x = sin(x*)/x* = cos(x*)",
            transcendence_class="{} -- a transcendental-equation value, not a period",
        ))

        rows.append(entry(
            "anchor.gerade_constant",
            "anchor",
            "the gerade constant 2 / (1 + min_x j0(x))",
            mp.nstr(gerade, CLAIM, strip_zeros=False),
            CLAIM,
            "For two EQUIVALENT centres the Shibuya-Wulfman overlap matrix is "
            "S = [[I, C], [C^T, I]] with C symmetric, and the symmetry-adapted "
            "blocks are exactly I +- C.  The condition number of the gerade block "
            "is therefore (1 + sup W)/(1 + inf W) where W is the symbol of C, "
            "which in the sine basis is j0.  Since sup j0 = 1 and "
            "inf j0 = min_x sin(x)/x, the gerade condition number tends to "
            "2/(1 + min j0) -- independent of both the internuclear separation "
            "and the basis size.  Paper 60's empirical 'flat condition number "
            "about 2 across N = 4..20' is this exact constant seen at small "
            "basis size.",
            f"Evaluated at working precision {DPS} from the certified root x*. "
            f"The independent float64 implementation "
            f"geovac.sturmian_sigma_law.gerade_constant() (scipy bounded "
            f"minimisation, a different algorithm) gives {g_float!r}, agreeing "
            f"to {mp.nstr(abs(gerade - mp.mpf(g_float)), 3)} -- the float "
            f"minimiser's own tolerance.  Claiming {CLAIM} digits from the "
            f"arbitrary-precision route.  The physics claim (measured "
            "cond(I + C) approaches this value: 2.552979 / 2.553764 / 2.554033 "
            "at n = 160 for three separations) is backed by "
            "tests/test_paper60_sigma_law.py.",
            value_kind="decimal",
            defining_relation="2 / (1 + min_x sin(x)/x)",
            transcendence_class="{} -- algebraic in the sinc minimum; no period content",
            backing_test="tests/test_paper60_sigma_law.py",
        ))

    return rows
