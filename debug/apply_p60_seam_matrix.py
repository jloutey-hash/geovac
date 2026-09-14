"""Add the two contraction-seam rows to docs/claim_test_matrix.md.

Appended after the sec:resource legitimacy row (the last Paper 60 row).
Idempotent.
"""
from __future__ import annotations

import sys

DOC = "docs/claim_test_matrix.md"
MARKER = "contraction window"
ANCHOR = ("Fire-tested with the tempting shortcut `X = P^-1/2` alone "
          "(the factor carrying the DST): FIRES on both halves |\n")

ROWS = (
    "| 60 | §molecular [SYMBOLIC + MEASURED] — the transcendental tagging owed "
    "under CLAUDE.md §4: the two π's of this section come from opposite sides of "
    "Paper 18's compactness boundary. The `π^2` of eq:sigma_law is TRUNCATION-side "
    "— the KMS constant c_1, the first Dirichlet eigenvalue of -d^2/dx^2 on the "
    "unit interval, which appears with NO Bessel function present anywhere "
    "(`n^2 min<θ^2>` -> π^2, Richardson 9.86949 vs 9.86960) and whose minimum the "
    "near-null direction ATTAINS (`n·rms(θ)` -> π, Richardson 3.14158). The "
    "`(2π)^-1/2` and `π/4` of eq:chirp_decay are CONTINUUM-side (Bessel "
    "normalisation + branch phase at the p->∞ pole). Both calibration-tier M2 — "
    "the `π^2·Q` half and the `√π·Q` half respectively | "
    "`tests/test_paper60_contraction_window.py```::test_band_minimum_is_pi_squared_"
    "with_no_bessel_present`` + ``::test_the_constant_tracks_the_zero_order_not_the_"
    "value_pi_squared`` + ``::test_near_null_direction_attains_the_band_minimum`` + "
    "``::test_factorisation_holds_on_the_near_null_direction`` | tracked "
    "`geovac/sturmian_sigma_law.py` (`theta2_band_matrix`, closed form — no "
    "quadrature) | **NEW 2026-09-12** | BACKED-SOUND. The Bessel-free route is the "
    "whole evidential point and it is guarded against being vacuous: a companion "
    "guard raises the symbol's zero from quadratic to quartic and requires the "
    "constant to COLLAPSE, so the π^2 cannot be arriving by construction. The "
    "near-null guard asserts the top direction reaches π/n AND that a lower "
    "singular direction is >1.5x worse — without the second half it would accept "
    "\"every direction looks the same\". Fire-tested 8 ways (corrupt the closed "
    "form; second singular vector; wrong curvature 1/24->1/12; quartic->quadratic; "
    "and the four weight plants below). Second route for c_1: the exact tridiagonal "
    "spectrum, `tests/test_paper60_kms_attribution.py` (v5.10.18) — two "
    "representations, one constant. rests on: eq:sigma_law; Paper 18 "
    "§\"Compactness as the source of discreteness\" (the tier vocabulary and the "
    "toll-per-decompactified-axis reading); the C23 re-attribution of the π/4 as a "
    "branch phase (row 548) |\n"
    "| 60 | §molecular [MEASURED] — the conditioning law is carried by the "
    "TRANSLATION, not by the metric that carries it: deform the metric by any "
    "smooth positive radial weight `W` on the Fock sphere (intra block = finite "
    "section of `W`, cross block = that of `W j_0`) and the generalised symbol is "
    "the quotient `W j_0 / W = j_0`, which does not see `W`. Measured at kR=2, "
    "n=160: collapse 0.4072 / 0.4123 / 0.4107 / 0.4164 for W = 1, 1+0.8cosχ, "
    "2+sinχ, e^-χ against π^2/24 = 0.4112, exponents -1.98..-2.01 | "
    "`tests/test_paper60_contraction_window.py```::test_smooth_weights_preserve_the_"
    "collapse_constant`` (4 weights) + ``::test_vanishing_weight_control_moves_the_"
    "constant`` + ``::test_weight_independence_holds_at_the_exponent_too`` + "
    "``::test_unit_weight_reproduces_the_tracked_sw_block`` | tracked "
    "`geovac/sturmian_sigma_law.py` (`weighted_blocks`, `generalized_sigma_max`) | "
    "**NEW 2026-09-12** | BACKED-SOUND, and the CONTROL is the load-bearing half: "
    "a weight that VANISHES at the degeneracy (`W = 1+cosχ`) moves the constant to "
    "0.828 while leaving the exponent at -1.97, so the agreement is a measurement "
    "rather than an insensitivity — the smooth-weight guard alone would pass on a "
    "pipeline that ignored `W` entirely. Anchored to tracked code by requiring "
    "`W = 1` to reproduce `(I, sw_cross_block)` to 1e-9, so the deformed machinery "
    "cannot be a self-consistent separate universe. Fire-tested: `weighted_blocks` "
    "ignoring its weight FIRES the control; making the control secretly smooth "
    "FIRES; `generalized_sigma_max` ignoring the intra block FIRES the smooth "
    "guard; a detuned kR FIRES the anchor. **SCOPE, stated in-paper because the "
    "overstatement is close by:** this is NOT independence of `V_0` — a "
    "position-space-local `V_0` acts on momentum space by CONVOLUTION, not "
    "multiplication, and leaves this class entirely. Consistent with the one "
    "measured out-of-class case (v5.10.15: the L^2 overlap degrades at the same "
    "exponent with a constant ~1.4x larger). Whether any position-space `V_0` "
    "preserves the constant is OPEN. rests on: eq:sigma_law |\n"
)


def main() -> int:
    with open(DOC, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED -- marker present; nothing done.")
        return 1
    if text.count(ANCHOR) != 1:
        print(f"ANCHOR not unique (count={text.count(ANCHOR)}); aborting.")
        return 2
    text = text.replace(ANCHOR, ANCHOR + ROWS)
    with open(DOC, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("applied: 2 claim-matrix rows added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
