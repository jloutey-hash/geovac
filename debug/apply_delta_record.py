"""Record the DELTA-verification run (v5.11.7). Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
DOD = "docs/qa/paper_60.done.md"
CL_ANCHOR = "## [v5.11.6] - 2026-09-12\n"
CL_MARKER = "## [v5.11.7]"

ENTRY = """## [v5.11.7] - 2026-09-13

**`/qa paper_60` DELTA-verification = DEFECTS, remediated. NOT a clean delta, so the certifying FULL run stays locked.** Four dimensions over the claim-impact set (not the byte diff), unseeded. Deterministic layer 14/14. **172 tests pass; zero mathematical defects in any dimension.**

### Scope, and why two LARGE findings sat outside it

Nominated from the three deterministic sources: the `\\cite` graph (2 documents), a topical-vocabulary sweep (7 further papers), and the declared `rests on:` edges (8 matrix rows). **Both LARGE findings landed outside that set, and by CATEGORY rather than accident:**

1. **The generated-artifact layer.** The 2026-09-11 Kac-Murdock-Szego credit was applied to `docs/certified_reference_values.md` and never to the generator. Measured: 1 mention in the generated doc, **0 in `entries_anchors.py`, 0 in the distributed JSON**. So the reference data this corpus offers outward already presented a 1953 theorem as an internal finding, and the next `generate_table` would have silently deleted the human-readable credit. Fixed in the generator and regenerated; the credit now survives in all three layers by construction. **Rule for future sweeps: nominate the GENERATOR whenever you nominate a generated doc.**
2. **The recall layer.** `memory/avery_method_and_prior_art_gaps.md` auto-loads into every session and carried four pre-2026-09-12 readings -- one as an OPERATIVE INSTRUCTION ("do claim the identification") folding in the translation reading that is prior art on three counts; the WITHDRAWN frames reading asserted as live paper content; the demoted "Proposition D" label; and `west_ruedenberg2013` described as "the one named source that uses an SVD/principal-angle construction", which its own abstract contradicts. This layer sits outside the repo and outside C16/C21/C22, and the same file had been self-caught for this exact class one day earlier.

### The over-correction, and then the over-correction of the correction

**The 2026-09-12 chain fix added "values are truncated, not rounded". That is false**, provably without measuring anything: the 1s^2 rung is exactly `-(27/16)^2 = -2.84765625`, which truncates to `-2.8476` while the paper prints `-2.8477`. All four rungs round. The convention was read off the old THREE-decimal display and carried over to four. Corrected, registry convention text with it.

**Then the reverse.** The code reviewer filed the c=3 box-rule column as a quadrature floor that never converges, making the appendix's `~1e-7` only a bound. Its three tabulated points showed 12-13x falls per resolution step; the PM verified those points, agreed, and re-tiered the claim matrix and the test docstring to "consistent but not confirmed". **The reviewer then extended its own ladder two points further, found a plateau, and withdrew the finding.** Independently re-measured here, full ladder at `n_max=8`:

| npts | c=1 | c=2 | c=3 | step |
|--:|--:|--:|--:|--:|
| 12 000 | 3.6130e-02 | 6.5629e-04 | 2.5709e-05 | -- |
| 40 000 | 3.6166e-02 | 6.8727e-04 | 2.1107e-06 | 12.18x |
| 100 000 | 3.6169e-02 | 6.8984e-04 | 1.5592e-07 | 13.54x |
| 250 000 | 3.6170e-02 | 6.9025e-04 | 1.5676e-07 | **0.99x** |
| 600 000 | 3.6170e-02 | 6.9031e-04 | 2.0598e-07 | 0.76x |

It PLATEAUS. The paper's figure is right, the original "VINDICATED" wording was right, and the correction was the error. Restored with the full ladder, which is stronger evidence than the original had. `test_c3_sits_below_the_quadrature_floor` -- a guard that PASSED while naming a false mechanism -- replaced by one asserting convergence.

**The lesson is the PM's:** a measurement-based finding was accepted on the reviewer's sample instead of extending the sample. Third grid-convergence trap of this arc, and the only one that reached the record.

### The rest

**Claims (7 SMALL):** the s-only sector label reached the conclusion and not the abstract or body -- the same silent-basis-mix defect the chain fix was written to close, one locus away; the third-lever scope omitted `s`-sector at two loci ("s-sector" occurs exactly ONCE in the paper); "perfectly-conditioned gerade sector" survived at its own locus; the water-control sentence contradicted itself (clause A says preconditioning is not doing the work, clause B measures it halving the exponent 1.96 -> 0.98 -- the rotation supplies BOUNDEDNESS); "complete set at one scale" unqualified inside a molecular comparison; the criteria file froze two RETIRED values as goalposts.

**Citations (2 SMALL):** "Bernstein" survived uncited at one locus, naming a mechanism the cited source does not use; and Theorem 73 was described as constraining polynomial approximation when it lower-bounds the block-encoding QUERY COUNT. Both corrected conservatively. **All three new bibitems CONFIRMED at source, and the FULL run's recommendation to drop the Wulfman-Takahata prior-art count was wrong** -- the work exists and its abstract names E4, R5, O(4,1) verbatim. The three-count concession is fully supported.

**Claim-impact (5 SMALL beyond the two LARGE):** the walls register stated the breach scope as "M = 2, 3", covering the collinear case the paper declines -- and BeH2 and CO2 are collinear and in this corpus's own library; the same register still named a 1954 result as our proposition; `code_architecture.md` used "derived", the word the owning module's docstring forbids; and the synthesis and claims register quoted the accuracy floor as a single value where the paper says the BRACKET is what it quotes -- with 6.4 sitting BELOW the bracket's own lower endpoint. Found independently by two reviewers.

**Code (1 MATERIAL, 6 NITs):** `tests/test_paper60_full_shell_family.py` was cited by ZERO claim-matrix rows -- the 2026-09-12 remediation closed the test leg of "no test, no registry key, no C17 family" and left the registration leg open, so deleting the file would have left every gate green (C22 checks rows->tests, never tests->rows). Row added. NITs cleared: a residual `x == x` that had been supplemented rather than replaced; a commutation pin whose docstring claimed protection it cannot provide; a superlinear guard that `fit > 1.02` would let 1.5 through, now bracketed, plus the paper's 0.893 global fit asserted for the first time.

### What the delta found nothing wrong with

158 test executions across five orderings, 0 failures. 8 of 8 guards fired. Both replaced guards genuinely discriminate -- the null direction is specific, not "any rotation": sweeping the angle, only 1 of 200 random angles reaches the positive test's bound. Every closed form, table and exponent recomputed by an independent route matched, `theta2_band_matrix` by two routes neither of which is the in-repo anchor's trapezoid. The `1e-300` clamp in `generalized_sigma_max` never engages. The per-state threshold sits 4.6x below the true separation and infinitely above an exactly-zero noise floor.

### Verdict

**DEFECTS.** A delta cannot return PASS by construction, and this one was not clean, so **the certifying FULL run remains locked** until a delta comes back clean.

"""

DOD_ANCHOR = "## Change log\n"
DOD_MARKER = "2026-09-13 — **DELTA-verification"
DOD_NEW = """## Change log
- 2026-09-13 — **DELTA-verification = DEFECTS, remediated. NOT a clean delta.** Four
  dimensions over the claim-impact set, unseeded; deterministic 14/14; 172 tests pass;
  **zero mathematical defects**. 2 LARGE, both OUTSIDE the nominated set by category —
  the generated-artifact layer (prior-art credit in the output but not the generator, so
  the distributed JSON lacked it and regeneration would erase the markdown copy) and the
  auto-loading recall layer (four stale readings, one an operative instruction to claim
  what was given away). 14 SMALL across claims/citations/claim-impact, plus one MATERIAL
  code finding (the new full-shell test cited by no row). **Two self-inflicted errors
  recorded:** the 2026-09-12 chain fix added "values are truncated, not rounded", false
  against the exact 1s² value; and the PM accepted a reviewer's under-sampled convergence
  finding and re-tiered a CORRECT claim, which the reviewer then withdrew and an extended
  ladder refuted (c=3 plateaus at ~1.6e-7 by 250k points). **Next: a clean delta before
  any FULL certifying run.**
"""


def main() -> int:
    n = 0
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) == 1:
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
        n += 1
        print("  ok   CHANGELOG.md")
    with open(DOD, encoding="utf-8") as fh:
        d = fh.read()
    if DOD_MARKER not in d and d.count(DOD_ANCHOR) == 1:
        with open(DOD, "w", encoding="utf-8") as fh:
            fh.write(d.replace(DOD_ANCHOR, DOD_NEW))
        n += 1
        print("  ok   paper_60.done.md")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
