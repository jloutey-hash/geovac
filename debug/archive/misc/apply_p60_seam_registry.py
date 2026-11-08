"""Register the contraction-seam numerals in debug/qa/numeric_registry.py.

Three keys, inserted before "p60_stateprep_overlap_exc":

  p60_window_richardson_pi2     n^2 min<theta^2> -> pi^2, measured with no
                                Bessel function present -- the evidence that
                                the pi^2 of eq:sigma_law is truncation-side.
  p60_window_rms_richardson_pi  n * rms(theta) -> pi on the near-null direction.
  p60_weighted_collapse_control the W-independence control (W vanishing at the
                                degeneracy), carried with the four smooth-W
                                values as a MATCHED SET in aliases.

Idempotent: refuses to run twice.
"""
from __future__ import annotations

import sys

REG = "debug/qa/numeric_registry.py"
ANCHOR = '    "p60_stateprep_overlap_exc": dict(\n'
MARKER = "p60_window_richardson_pi2"

NEW = '''    "p60_window_richardson_pi2": dict(
        value=9.86949, convention="constant: first-order Richardson limit in 1/n "
                                  "of n^2 * min<theta^2>, the minimal mean-square "
                                  "spread in the Fock polar angle theta = pi - chi "
                                  "over span{sin(a chi)}_{a<=n}; target pi^2 = "
                                  "9.8696044",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_window_constant_probe.py test F2, "
                   "from n = 20..640. This is the Kac-Murdock-Szego constant c_1 "
                   "recomputed in a SECOND representation -- a band-limited "
                   "concentration problem with NO Bessel function anywhere in it -- "
                   "which is the whole evidential point: it shows the pi^2 of "
                   "eq:sigma_law is a truncation constant, not Bessel content. "
                   "First route was the exact tridiagonal spectrum (v5.10.18, "
                   "tests/test_paper60_kms_attribution.py). Two routes, one "
                   "constant, per the independent-route rule.",
        aliases={9.84287: "raw n=640, unextrapolated",
                 9.81624: "raw n=320", 9.76331: "raw n=160"}),
    "p60_window_rms_richardson_pi": dict(
        value=3.14158, convention="constant: first-order Richardson limit in 1/n of "
                                  "n * rms(theta) for the near-null (top singular) "
                                  "direction of the SW cross block at kR = 2; "
                                  "target pi = 3.1415927",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_window_constant_probe.py test F1, "
                   "from n = 20..640. Says the degeneracy direction ACHIEVES the "
                   "band-limited minimum of p60_window_richardson_pi2 (same "
                   "constant, squared), i.e. the near-null direction is the optimal "
                   "concentrator at the p = 0 pole. Paired with F3, which checks "
                   "1 - sigma_max = (kR)^2 <theta^2>/24 on that direction to 0.3% "
                   "at n = 320 across kR = 0.5..4.",
        aliases={3.13733: "raw n=640, unextrapolated", 3.13309: "raw n=320"}),
    "p60_weighted_collapse_control": dict(
        value=0.828, convention="constant: the collapse (1-sigma_max)(n/kR)^2 at "
                                "kR=2, n=160 under the CONTROL weight W = 1+cos(chi), "
                                "which VANISHES at the degeneracy chi = pi -- the "
                                "value that must differ from pi^2/24 = 0.4112 for "
                                "the weight-independence measurement to mean "
                                "anything",
        q=None,
        provenance="MEASURED 2026-09-12, debug/p60_contraction_seam_probe.py test E. "
                   "MATCHED SET -- this control and the four smooth-weight values in "
                   "aliases must move together; quoting the agreement without the "
                   "control would report an insensitivity as a measurement. Fitted "
                   "exponent stays -1.97 here, so the control moves the CONSTANT "
                   "only, not the exponent.",
        aliases={0.40718: "smooth W = 1 (SW reference), n=160",
                 0.41229: "smooth W = 1 + 0.8 cos(chi), n=160",
                 0.41065: "smooth W = 2 + sin(chi), n=160",
                 0.41635: "smooth W = exp(-chi), n=160"}),
'''


def main() -> int:
    with open(REG, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED -- marker present; nothing done.")
        return 1
    if text.count(ANCHOR) != 1:
        print(f"ANCHOR not unique (count={text.count(ANCHOR)}); aborting.")
        return 2
    text = text.replace(ANCHOR, NEW + ANCHOR)
    with open(REG, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("applied: 3 registry keys added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
