"""Record DELTA-verification run #3 (v5.11.9). Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
CL_ANCHOR = "## [v5.11.8] - 2026-09-13\n"
CL_MARKER = "## [v5.11.9]"

ENTRY = """## [v5.11.9] - 2026-09-13

**`/qa paper_60` DELTA-verification #3 = DEFECTS, remediated. The cleanest run of this lineage: two of four dimensions CLEAN, zero LARGE, zero mathematical, and all three defects were one claim not reaching two summary surfaces.** Four dimensions over the claim-impact set, unseeded. Deterministic layer 14/14.

### Two dimensions came back CLEAN — a first for this paper

- **External citations: CLEAN-DELTA.** All six v5.11.8 citation edits verified at primary source, none over-corrected: the withdrawn Shibuya-Wulfman quotation (1965 abstract still unreachable, verbatim string gone, paraphrase claims only what Wulfman-Takahata 1967 and Red-Weatherford 2004 independently establish -- both re-confirmed at source, the 1967 reference's abstract naming E4/R5/O(4,1) exactly); the added `grochenig_leinert2006` inverse-closedness cite; the completed `lowdin1950` title and `rokob2008` venue; the `gslw2019` arXiv-numbering note; the KMS smoothness-hypothesis pointer.
- **Code/test-backing: CLEAN-DELTA** (narrow scope, no full-suite re-run). All three post-DELTA-#2 guards pass and fire on the specific retired value each names: the graded-mesh value guard rejects both 1.29e-7 and 1.5676e-7; its two-sided convergence leg rejects a 0.5x collapse and accepts a 1.010x plateau; the uniform-banding guard rejects the 0.98 exponent. The two c=3 docstrings tell one consistent story (converged value on the graded mesh; order-only on the still-climbing uniform ladder), and claim-matrix coverage is intact.

### The three defects, all remediated, all the same claim

The v5.11.8 water-control fix (the null-direction rotation is load-bearing; the preconditioner alone does not suffice) reached the paper body, abstract, conclusion and synthesis -- but not two summary surfaces, and left a stale number inside the corrected body sentence.

1. **`papers/group2_quantum_chemistry/paper_60...tex` L1338 (claims dimension) -- a stale number inside the corrected sentence.** The rewrite fixed the exponent to N^1.94 but kept the displayed low value `2766`, which is the N=48 INTERIOR point, not the low endpoint; `2766 -> 42008` over the table's N=12..192 reproduces the discredited N^0.98. **The claims reviewer proposed 729.2 as the fix -- that was itself wrong**, the N=24 value from a different measurement grid. Measured at the paper's own N grid (N=12/48/192): uniform 194.0 / 2766.0 / 42007.7, so the low endpoint is 194 and the exponent N^1.94 (raw N^1.96 on this grid). Corrected to `194 -> 42008 ... N^1.94`, with the interior 2766 kept as a named aside. **Verifying the reviewer's number before printing it prevented the sixth wrong value on this one claim.**
2. **`papers/INDEX.md` L88 (claim-impact) -- an omitted requirement.** The status map said the preconditioner reaches water's A_1 block with no mention of the rotation. C16 had been widened to this file and caught its token-matchable defects, but an OMITTED requirement is an absence, not a token. Caveat added.
3. **`docs/walls/register.md` L111 (claim-impact) -- a stale control argument.** The operative dispatch register cited the BLIND uniform control (`blockdiag(P,P)`, "leaves the growth intact, so the rotation is doing the work") as the discriminating evidence -- the exact defect the paper carried before v5.11.8. That control commutes with the rotation (3e-13) and cannot discriminate it. Replaced with the valid selective-unrotated control (N^3.79, 106x worse than untreated).

### The category nobody had nominated

**Paraphrase-level requirement-omission surviving token gates in already-nominated files.** Both claim-impact defects sit in files C16 was widened to in v5.11.8. The token gate fixed their string-matchable staleness (the "Proposition D" label, the "6.44" literal) but by construction cannot see an omitted requirement or a stale argument stated in the file's own words. The lesson, now demonstrated twice: **when C16 is widened to a file for a token, the claim-impact reviewer must RE-READ that file's Paper-60 sentence for direction and omission, not trust the token gate.**

### The process win, stated because last round's process was the finding

DELTA #2's worst defect was the target moving under the reviewers five times. This round the tree was frozen at dispatch and **no edit was made until all four reviewers returned** -- so every verdict covers the state it was formed on. That defect did not recur.

### The trajectory

DELTA #1: 2 LARGE + many SMALL, two leaked categories. DELTA #2: 0 LARGE but a remediation that never landed, a false plateau, five corrections to one number. DELTA #3: 0 LARGE, 0 mathematical, 3 SMALL all on one claim, two dimensions fully clean, and the one reviewer-proposed number that was wrong was caught before it shipped. The defects are converging toward mop-up of a single fix. **This is still DEFECTS, not a clean delta** -- the three fixes are summary-surface edits whose own cleanliness a fourth delta would confirm -- so the certifying FULL run stays locked. But the distance to a clean delta is now three summary edits, not a mathematical question.

### Verdict

**DEFECTS.** Remediated. A delta cannot return PASS; a clean delta remains the precondition for the certifying FULL run.

"""

DOD = "docs/qa/paper_60.done.md"
DOD_ANCHOR = "## Change log\n"
DOD_MARKER = "2026-09-13 — **DELTA #3"
DOD_NEW = """## Change log
- 2026-09-13 — **DELTA-verification #3 = DEFECTS, remediated. Cleanest of the lineage.**
  Four dimensions; **citations and code both CLEAN-DELTA**; 0 LARGE, 0 mathematical.
  Three SMALL, all the v5.11.8 water-control fix (rotation is load-bearing) not reaching
  two summary surfaces plus a stale number in the corrected body sentence: paper L1338
  (`2766` was the N=48 interior point, not the low endpoint — measured 194 at the paper's
  own grid; **the claims reviewer's proposed 729.2 was itself wrong and was caught before
  printing**), `papers/INDEX.md` (omitted rotation caveat), `docs/walls/register.md` (cited
  the blind uniform control as discriminating evidence). New category: paraphrase-level
  requirement-omission survives token gates even in files C16 was widened to. **Process:
  the tree was frozen and no edit made until all four reviewers returned — DELTA #2's
  moving-target defect did not recur.** Public-web layer stale but a known PI-gated gap,
  unchanged. **Next: a clean delta before any FULL certifying run.**
"""


def main() -> int:
    n = 0
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) == 1:
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
        n += 1
        print("  ok   CHANGELOG v5.11.9")
    else:
        print("  skip CHANGELOG (already applied or anchor moved)")
    with open(DOD, encoding="utf-8") as fh:
        d = fh.read()
    if DOD_MARKER not in d and d.count(DOD_ANCHOR) == 1:
        with open(DOD, "w", encoding="utf-8") as fh:
            fh.write(d.replace(DOD_ANCHOR, DOD_NEW))
        n += 1
        print("  ok   paper_60.done.md")
    else:
        print("  skip done.md (already applied or anchor moved)")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
