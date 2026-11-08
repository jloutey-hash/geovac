"""Append the follow-on section to the contraction-seam sprint memo."""
from __future__ import annotations

import sys

M = "debug/sprint_contraction_seam_memo.md"
MARKER = "## 8. Follow-on items"

SECTION = """
---

## 8. Follow-on items (2026-09-12, PI-directed)

### 8.1 Paper 34 candidate (e) logged in §VIII -- DONE, still PI-gated

Verified first that the gap is real: the string "translat" appears **nowhere**
in Paper 34, and the only Shibuya-Wulfman mentions (§III.22, multipole
expansion) name the SW *basis* expansion as a **contrast** -- something that
does NOT terminate -- never as a projection. §III.11 covers Wigner-`D`
*rotation* between centres and explicitly preserves rationality up to
`Q[sqrt2,sqrt3,sqrt6]`.

Candidate (e) is now logged in §VIII's open-question list alongside (a)-(d),
**not** written into §III, per the tag-transcendentals STOP rule. It carries:
the two-sided M2 signature (truncation prices `pi^2`, symbol asymptotics price
`sqrt(pi)`); the caution that the contraction reading is prior art and would be
cited rather than claimed, with the antipodal parity caveat; and a clause
reconciling it with the Sprint-3 "structurally complete" verdict -- that verdict
is scoped to the master Mellin engine's *accounting*, and (e) is M2 on both
sides, so what it adds is a **slot, not a mechanism**.

Promotion to a §III entry remains a PI call. group6 gates all PASS.

### 8.2 A gate for the splice class -- BUILT, PROPOSED, not adopted

`debug/qa/check_prose_continuity.py` (proposed **C24**). The 2026-09-12
self-catch was that an applier's anchor ended mid-sentence, and **every gate
passed** -- C10 compiles references, not prose. The class is not new: the
criteria document records that C20's own registration block was spliced into
the middle of the C19 sentence. Script-driven paragraph insertion is the
corpus's standard editing mechanism, so the class is systematic.

Detects two independent signals, each on its own: a prose paragraph that ends
without terminal punctuation, and one that opens with a lowercase ordinary
word. Neighbour-aware, so prose running INTO a display equation and prose
resuming AFTER one are both exempt, as is prose introducing a displayed
theorem.

**Measured: 0 findings across all 70 papers in 9 scopes**, while firing on both
halves of the real defect when it is re-planted into Paper 60.

*Two corrections it needed, both found by `debug/qa/_prose_gate_probe.py`, which
re-plants the real defect and requires the gate to fire.* The first draft
reported PASS on three scopes with zero findings and **would have shipped as a
gate that examines nothing**:

1. Environment tracking counted `\\begin{document}`, which never closes until
   the last line, so every paragraph in the body was skipped as "inside an
   environment".
2. The conjunction was the wrong shape. The draft required one seam to BOTH end
   unterminated and be followed by a lowercase opener, on the reasoning that
   this keeps false positives near zero. But when a paragraph is spliced into
   the middle of a sentence, the two halves land at **opposite ends of the
   insertion**. The real defect triggered neither half.

Refinement then took the corpus from 107 findings to 0 without losing the
detection: inline `\\begin{smallmatrix}` was disqualifying whole prose
paragraphs; preamble macro blocks were being judged as prose; trailing
`\\checkmark` left whitespace the closer-strip did not remove.

**Adding a QA criterion is a gate change**, hence a minor version and a PI
decision (the precedent is C23 at v5.11.0). Proposed, not adopted; the script
and its probe stand on their own until then.
"""


def main() -> int:
    with open(M, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    with open(M, "a", encoding="utf-8") as fh:
        fh.write(SECTION)
    print("appended: follow-on section")
    return 0


if __name__ == "__main__":
    sys.exit(main())
