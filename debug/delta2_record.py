"""Record DELTA-verification run #2 (v5.11.8). Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
CL_ANCHOR = "## [v5.11.7] - 2026-09-13\n"
CL_MARKER = "## [v5.11.8]"

ENTRY = """## [v5.11.8] - 2026-09-13

**`/qa paper_60` DELTA-verification #2 = DEFECTS, remediated. Still not a clean delta, so the certifying FULL run stays locked.** Four dimensions over the claim-impact set, unseeded. Deterministic layer **14/14**. **Zero mathematical defects in any dimension**, and every external theorem re-derived to its source text. Three findings are the PM's own from the previous round.

### The remediation that was recorded as done and was not in the corpus

`debug/delta_fix_04_last_three.py` accumulates edits in memory and `return 2`s from inside its edit loop on a single stale anchor -- **before** its write loop. So one bad anchor discards every edit in the script, silently, while the run reports a failure nobody reads as total. Three fixes written 2026-09-13 (the walls register's "Proposition D" at three loci; the accuracy-floor bracket in the synthesis and the claims register) **were never in the corpus**, and `CHANGELOG.md` v5.11.7 recorded them as remediated. Its sibling `delta_fix_03c_remaining.py` documents the identical failure one step earlier, in its own docstring.

Found by the claim-impact reviewer, which checked whether the fixes were *present* rather than whether they had been *written*. Rerun: all five landed. **Every applier written since writes first and reports misses**, so a miss can no longer discard a match.

### The control that could not test what it was named after

Paper 60 claimed *"a control confirms that it is the rotation that bounds the growth ... the naive `blockdiag(P,P)` without it leaves the growth unbounded (2766 -> 42008 ... an exponent of 0.98, halved from the raw column's 1.96)."* Both halves are wrong, and the repo already knew the first:

* `blockdiag(P,P)` is `I2 (x) tri(1,2,1)` and the rotation is `V (x) I`. **They commute** -- measured to `3e-13` at `N=192` -- so that control returns the identical spectrum in either frame and is blind to the rotation by construction. The 2026-09-12 FULL run had recorded exactly this, and the fix reached the TEST DOCSTRING and not the paper.
* The exponent is **not** halved. Measured:

| column | N=24 | N=192 | exponent |
|---|--:|--:|--:|
| raw | 698.8 | 41 700 | N^1.967 |
| uniform `blockdiag(P,P)` | 729.2 | 42 008 | N^1.950 |
| selective `blockdiag(P,I)`, unrotated | 1 690 | 4 437 000 | N^3.791 |
| selective `blockdiag(P,I)`, rotated | 42.4 | 44.1 | N^0.018 |

Banding alone does essentially nothing. **The discriminating control is the selective preconditioner in the UNROTATED frame, and it ends 106x WORSE than leaving the metric untreated.** The conclusion the paragraph wanted is true and the real evidence is far stronger than what was printed. Rewritten from measurement, and the caveat now travels with the claim to the abstract, the conclusion and the synthesis -- where the gain had been stated three times with no control at all.

### The plateau that was not one -- third grid trap of this arc, and mine

v5.11.7 reverted a re-tiering of the box rule and restored "VINDICATED", on the strength of a `0.99x` step at `npts=250000`. **The code reviewer sampled `npts=400000` -- a point nobody had taken -- and the column climbs straight past it.** Independently re-measured here, digit for digit:

| npts | c=3 relative | step |
|--:|--:|--:|
| 100 000 | 1.55917e-07 | -- |
| 250 000 | 1.56764e-07 | 1.005x |
| 400 000 | 1.93056e-07 | 1.232x |
| 600 000 | 2.05980e-07 | 1.067x |

**The v5.11.7 table already said so.** Its last row rises 31% over the row marked `0.99x`, printed as a "0.76x step" directly beneath the word PLATEAU, unreconciled. The flat step is a **crossing**: the quantity is relative to a `c=5` box whose own truncation (`7.69e-08` absolute, against `c=3`'s `1.29e-07`) is comparable and partially cancels there.

**The conclusion survives; the mechanism and the third digit do not.** App. A's `~1e-7` is vindicated as an ORDER by every route -- uniform mesh `1.29e-7`, graded mesh `2.16e-7`, in-repo relative `1.6-2.1e-7` -- and the revert was right. But the two discretizations differ by ~1.7x, so no converged value should be quoted and v5.11.7 quoted one. The guard built on the false mechanism pinned `(100000, 250000)`, **the only adjacent pair in the ladder under a 10% window** (measured: 0.54%, 23.15%, 6.69%), and the reviewer fired it by planting MORE resolution -- a change with no physics in it. That is verbatim the defect `delta_fix_05` was written to remove from the guard before it. Replaced by one asserting the ORDER and the NON-collapse, which is what is robust, naming both wrong answers.

**The lesson is the PM's, and it is the mirror of the previous one:** v5.11.7 recorded accepting a reviewer's under-sampled finding; this round recorded publishing an under-sampled finding of its own, in the opposite direction, in the correction of that correction.

### Two instrument defects, both of which were hiding live findings

1. **The C16 entries added 2026-09-12 shipped with a NARROWER `files` list than the older Paper-60 entries beside them** -- missing `docs/claims_register.md`, `docs/code_architecture.md`, `docs/topic_to_paper_lookup.md`, `papers/INDEX.md` and `docs/qa/paper_60.done.md`. That is why the walls register's "Proposition D" and the single-value floor survived two runs. Widening the six took C16 from **PASS to FAIL with 7 live loci**, all now closed (two were legitimate withdrawal records missing the standardized token; one was genuinely stale). The widened entry is proven to discriminate: planting the retired wording FAILS, restoring PASSES. Correctly-flagged occurrences went 36 -> 40.
2. **`generate_table --check` compared only the `value` field**, so the prior-art credit DELTA #1 found missing from the generator could be deleted again from `method`/`provenance` with the gate green. **A guard added after a defect that cannot see the field the defect lived in is not guarding it.** Widened to attribution fields and fire-tested by deleting that exact credit: 1 drift detected, 0 when restored. `render_markdown` also gained the `provenance` branch it never had.

### Two affected CATEGORIES nobody has nominated

DELTA #1's two LARGE findings sat outside the nominated set by category. So do these.

1. **The public web-rendering layer** (`viz/public/papers/*.html`, `index.html`, `sitemap.xml`) -- the outward-facing sibling of DELTA #1's generated-artifact LARGE. The group2 synthesis renders an abstract saying *"the nine Group 2 papers"* where the source says twelve, and **Papers 58, 59, 60 and 61 have no public page at all**, so not one of the nine changed claims has any public rendering. Built from a stale `debug/data/zenodo_manifest.json`. **DECLARED GAP, not fixed:** rebuilding that manifest touches the Zenodo/DOI distribution surface, which is PI-gated.
2. **`docs/walls/register.md` is unreachable by every gate** -- the string `docs/walls` appears nowhere in `debug/qa/` -- yet it carries an operative dispatch rule. Same "operative instructions outside the gate" shape as the recall layer that was DELTA #1's second LARGE, but inside the repo. Closed by the widening above.

### Citations: 0 LARGE, and one quotation withdrawn

Both nominated fixes verified correct at the primary source: **Theorem 73 of `gslw2019` IS "Lower bound for eigenvalue transformation"** and bounds applications of the block-encoding, and GSLW themselves pair it with Corollary 67 in their own proof paragraph -- so the chain is the source's, not ours. The Bernstein removal is complete. **`wulfman_takahata1967` is real** (JCP 47(2), 488-498) and its abstract names E4, R5, O(4,1) verbatim -- the FULL run's drop recommendation was wrong and overturning it was right, now confirmed a second time independently.

Fixed: the Bernstein removal left the inverse-closedness mechanism **asserted with no cite**, while `grochenig_leinert2006` -- which is exactly that mechanism -- sat in the bibliography never cited; `lowdin1950` carried no title; `rokob2008`'s venue was unresolvable; the `gslw2019` bibitem leads with the STOC abridgement while the theorem numbers are the arXiv version's. **Withdrawn:** a verbatim quotation attributed to Shibuya-Wulfman's 1965 abstract. The bibliographic record is exact, but the abstract is paywalled and three independent retrievals failed to reach the quoted string -- so it is now paraphrased and marked unverified rather than asserted. The prior-art surrender is unaffected: it rests on three counts and the other two were verified at source.

### The rest

**Claims (7 further):** the He chain's *"every rung at the same `n_max=10` family"* is false for the leading `1s^2` single configuration; *"the largest computed basis `K=452`"* is contradicted by the paper's own `K=514` loci **and by the registry's own alias** (it is the largest basis where BOTH roots were computed); the previous round's completeness fix left a clause with no antecedent pointing the wrong way; the synthesis dropped the `s`-sector scope the paper carries at all three of its own loci; and "derived" -- the prohibited word for the conditioning law -- survived at two loci above the `[PRIOR ART]` paragraph.

**Claim-impact (4 further):** `papers/INDEX.md` quoted the floor as `6.4` mHa, **below the bracket's own lower endpoint**, and still described molecular conditioning as *"grows with basis"*, the pre-breach reading; the `eq:W_diagonal` row stated a `5e-11` entrywise check as the claim's evidence where the DoD says flatly it is an implementation check and not the evidence (an UNDERCLAIM); and the **variational bound proved by Sylvester inertia has no claim-matrix row at all** though every excited-state number in Sec.4 rests on it -- flagged in DELTA #1 and now a declared coverage gap raised to the PI.

**Code (4 further):** a debt this delta had PAID was still recorded as owed; the global-fit guard pinned a value and not the window it names (K=35-220 -> 0.8865, K=84-220 -> 0.8995, K=56-165 -> 0.8885 all passed a +-0.01 band, while grid and box variants move it by <1e-4 -- tightened to +-0.004); the box-rule row omitted the very test backing its own note; and App. A's exponent-convergence sentence has no test, now declared.

### The PM's own process defect

**The target moved under the reviewers.** The code reviewer finished its mandated run at 08:25; the PM then applied six appliers between 08:30 and 08:43, rewriting the paper, five documents, the synthesis and a test file that gained a 13th test the run had never collected. **A DELTA verdict formed on a moving tree certifies a state that no longer exists.** The `/qa` protocol already says to snapshot into a worktree when the PM will keep editing, and the PM did not. The reviewer caught it itself, re-ran, and independently reproduced every new number -- raw `N^1.9666`, uniform `N^1.9497`, selective `N^3.7913`, `106.4x` at `N=192` -- so the content is sound and the process is the finding. **Standing rule for the next run: snapshot before dispatch, or do not remediate until every reviewer has returned.**

### What the delta found nothing wrong with

Every closed form, table and exponent recomputed by an independent route matched. The generator is idempotent -- 66/66 entries, 0 field diffs, `render_markdown` byte-identical. The two thresholds tightened in the previous round were **correct and not overshot**: across four grid/box variants the `T'` fit reads 1.0727 every time (spread `1e-4` against an 0.08-wide bracket) and the global fit 0.8932 every time. The commutation-pin docstring's honesty rewrite is accurate, and the discrimination it points to genuinely lives in the selective control. 95 tests pass with `--slow`; the deterministic layer is 14/14.

### Verdict

**DEFECTS.** A delta cannot return PASS by construction, and this one was not clean. **The certifying FULL run remains locked** until a delta comes back clean.

"""


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER in t:
        print("already applied")
        return 0
    if t.count(CL_ANCHOR) != 1:
        print(f"  MISS anchor count={t.count(CL_ANCHOR)}")
        return 3
    with open(CL, "w", encoding="utf-8") as fh:
        fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
    print("  ok   CHANGELOG v5.11.8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
