"""Version note for the literal-registration pass (v5.11.11). Idempotent."""
from __future__ import annotations
import sys

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.10] - 2026-09-13\n"
MARKER = "## [v5.11.11]"
ENTRY = """## [v5.11.11] - 2026-09-13

**Paper 60 literal-registration pass -- FULL-run cert-blocker 1 partly discharged, partly reclassified, and it surfaced a finding.** Five declared-debt literals were MEASURED from tracked code (`debug/fullrun_measure_literals.py`) and registered under C21 with `\\gvq` annotations: `p60_l2_inflation_hydrogenic` (1.19), `p60_l2_inflation_sturmian` (3.33), `p60_sw_cond_exponent` (1.85, window N=12..16; full-range 1.80, asymptote N^2), `p60_l2_overlap_exponent` (**1.71** -- the paper printed 1.70, corrected to the measured value), `p60_water_a1_exponent` (1.96, the raw `_water_A1` "N^1.96 here" value). C21 green, compile clean, registry self-test 18/18.

**The finding:** three of the debt literals do NOT reproduce from the tracked test suite, which the earlier "correct; unregistered" label had hidden -- the Gaussian ratio exponents (tracked `_gaussian_metric` gives ratio-1.6 -> N^6.30 with κ~1.5e6 at N=16, not the paper's ~1e5, and ratio-3 -> N^1.54 not 1.4; the paper's own caveat already calls these "illustrative rather than a fixed factor"), the water N^1.97 sec:molecular probe (a different three-center-SW construction from `_water_A1`), and the `tab:resource` row (`[RESOURCE MODEL]`, O(1) cancels in ratios). These need their drivers promoted to `tests/` or reclassification, deferred to the group2 review. **Sec.15 rule 3 honoured: every registered value was measured; the three that could not be measured from tracked code were NOT registered rather than guessed.** Paper 60 full cert deferred to the group2 review (PI direction).

"""


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("skip (already applied)"); return 0
    if t.count(ANCHOR) != 1:
        print(f"MISS anchor count={t.count(ANCHOR)}"); return 3
    with open(CL, "w", encoding="utf-8") as fh:
        fh.write(t.replace(ANCHOR, ENTRY + ANCHOR))
    print("ok CHANGELOG v5.11.11")
    return 0


if __name__ == "__main__":
    sys.exit(main())
