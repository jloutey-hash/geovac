"""Update the DoD declared-debt table after the 2026-09-13 registration pass.

Registered (measured from tracked code, C21 annotated): C8.1 (Q^3.33/Q^1.19),
C8.6 (SW N^1.85 window + L2 N^1.71), and the water 'N^1.96 here' value.

NOT registrable from the tracked suite -- a finding, promote-or-reclassify at
the group2 review: the Gaussian ratio (illustrative by the paper's own caveat;
tracked helper gives 6.30/1.54 not 6/1.4), the water N^1.97 sec:molecular probe
(a different three-center-SW construction), the tab:resource row (RESOURCE MODEL,
modelled), and the H2+/H2 toy figures (not yet measured this pass).

Idempotent.
"""
from __future__ import annotations

import sys

DOD = "docs/qa/paper_60.done.md"

OLD = """| locus | literal | status |
|---|---|---|
| branch criterion 4 | Gaussian ratio-dependence (`N^6` / `N^1.4`) | correct; unregistered |
| W5 | water A1 `cond∼N^1.97` | correct; unregistered |
| W8 | H2+ / H2 toy-validation figures | correct; unregistered |
| C8.1 | naive L2 inflation `Q^3.33` vs `Q^1.19` | correct; unregistered |
| C8.6 | SW `cond(S)∼N^1.85` (vs L2 `N^1.70`) | correct; unregistered |
| C8.7 | the `tab:resource` row | correct; unregistered |
| C8.9 | water gerade-lever-fails numbers | correct; unregistered |"""

NEW = """| locus | literal | status (updated 2026-09-13) |
|---|---|---|
| C8.1 | naive L2 inflation `Q^3.33` vs `Q^1.19` | **REGISTERED** `p60_l2_inflation_sturmian` (3.33, meas 3.3328) / `p60_l2_inflation_hydrogenic` (1.19, meas 1.1909); annotated eq:blowup |
| C8.6 | SW `cond(S)∼N^1.85` (vs L2 `N^1.70`) | **REGISTERED** `p60_sw_cond_exponent` (1.85, meas 1.8518 on window N=12..16; full-range 1.80, asymptote N^2) / `p60_l2_overlap_exponent` (meas **1.71**, paper's 1.70 corrected); annotated sec:molecular |
| W5 / C8.9 | water A1 `cond∼N^1.96` (raw, `_water_A1`) | **REGISTERED** `p60_water_a1_exponent` (1.96, meas 1.9580); annotated sec:resource. The **N^1.97 sec:molecular probe (19.9→698 over N=6..36) is a DIFFERENT construction, still unregistered** — see note |
| branch criterion 4 | Gaussian ratio-dependence (`N^6` / `N^1.4`) | **NOT REGISTRABLE — finding.** The tracked `_gaussian_metric` gives ratio-1.6 → N^6.30 (κ~1.5e6 at N=16, not the paper's ~1e5) and ratio-3 → N^1.54 (not 1.4). The paper's own caveat calls these "illustrative rather than a fixed factor"; treat as illustrative, or promote the original driver |
| W8 | H2+ / H2 toy-validation figures | unregistered (not measured this pass; toy validations) |
| C8.7 | the `tab:resource` row | `[RESOURCE MODEL]` — O(1) cancels in ratios per the caption; modelled, not a measured canonical value |

**Registration pass 2026-09-13 (FULL-run cert-blocker 1).** Five literals were
measured from tracked code (`debug/fullrun_measure_literals.py`) and registered
with `\\gvq` annotations; C21 green. **The pass surfaced a finding the earlier
"correct; unregistered" label hid:** three of the debt literals — the Gaussian
ratio exponents, the water `N^1.97` sec:molecular probe, and the `tab:resource`
row — do NOT reproduce from the tracked test suite (the first two trace to
constructions that live only in `debug/` drivers or are modelled). They cannot
be registered as *measured* values without promoting their drivers to `tests/`
first, which is a group2-review task. So the literal-registration debt is now
**partly discharged (5 keys) and partly reclassified** (2 need driver promotion,
1 is modelled) — no longer a flat "register 7 before cert."
"""

OWED_OLD = ("**Owed before the next FULL\ncertifying run:** register each under C21 or a C17 family (Sec.15 rule 3 applies —\n"
            "measure or cite, never guess), then delegate. Until then a reviewer treats them\n"
            "as literals and checks them against the body.")
OWED_NEW = ("**Status 2026-09-13:** the tracked-reproducible five are registered (above);\n"
            "the three that are not tracked-reproducible (Gaussian ratio, water N^1.97 probe,\n"
            "tab:resource) need their drivers promoted to `tests/` or reclassification, deferred\n"
            "to the group2 review. Sec.15 rule 3 was honoured — every registered value was\n"
            "MEASURED, and the three that could not be measured from tracked code were NOT\n"
            "registered rather than guessed.")


def main() -> int:
    with open(DOD, encoding="utf-8") as fh:
        t = fh.read()
    n = 0
    if "REGISTERED** `p60_l2_inflation_sturmian`" in t:
        print("  skip table (already applied)")
    elif t.count(OLD) == 1:
        t = t.replace(OLD, NEW)
        n += 1
        print("  ok    debt table updated")
    else:
        print(f"  MISS table anchor count={t.count(OLD)}")
        return 3
    if OWED_NEW[:40] in t:
        print("  skip owed-note")
    elif t.count(OWED_OLD) == 1:
        t = t.replace(OWED_OLD, OWED_NEW)
        n += 1
        print("  ok    owed-note updated")
    else:
        print(f"  MISS owed-note count={t.count(OWED_OLD)}")
    with open(DOD, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
