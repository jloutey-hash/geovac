"""Record the Monkhorst-Jeziorski verdict in the CHANGELOG and the sprint memo."""
from __future__ import annotations

import sys

CL_ANCHOR = "### Owed (PI items)\n"
CL_MARKER = "### Monkhorst-Jeziorski 1979"
CL_NEW = """### Monkhorst-Jeziorski 1979: a standing corpus claim retracted, and the better answer underneath it

PI-directed read of the source flagged earlier the same day. Memo `debug/lit_scan/monkhorst_jeziorski_1979_memo.md`. **Verdict SECONDARY-QUOTED** -- abstract and bibliographic record verified at source on two independent indexes and re-verified here via Crossref; the **two-page body is closed with zero repository copies anywhere and is UNREAD**, so the mechanism below is reconstructed from the lineage and is labelled as such in the paper. (The DOI the dispatch supplied was wrong; the correct one is `10.1063/1.438337`, vol. 71(12), 5268-5269.)

**The abstract is more direct than the title.** For a many-center **one-electron** system, eigenvalues follow from "diagonalizations of simple overlap matrices", and "the problems of many-center integrals and instabilities due to overcompleteness of basis sets do not appear at all."

**Verdict on the two suspended memory claims, which had been travelling together:**

- **"No prior art for secular-matrix norm growth" SURVIVES**, narrowed. No condition number, spectrum, smallest eigenvalue or growth law is in evidence; the reference list carries no numerical-linear-algebra source; no citing paper cites it for a conditioning result. Rests on abstract + reference list + page count, not a read.
- **"Conditioning is a blind spot in the whole Avery canon" DOES NOT SURVIVE and is RETRACTED as phrased.** An abstract naming "instabilities due to overcompleteness" is explicit published engagement. The defensible remainder is the first claim: nobody *priced* it.

**The resolution is worth more than the attribution, and it inverts the expected answer.** The PM's predicted resolution -- that their claim lives in the `V_0` metric and ours in `L^2` -- was **wrong**. It is the *same pencil*: same translation-phase matrix, and `sigma_max -> 1` at `p = 0` is a property of the symbol, indifferent to who is looking. The difference is *extraction*. Letting the overlap enter only multiplicatively -- a determinantal root search in the scale, or direct diagonalization -- **inverts nothing**, so a near-null direction yields a harmless spurious branch instead of amplified error; and the basis is exactly orthonormal in the metric actually used, which is this corpus's own measured "the intra-center block is exactly the identity", so the `L^2` Gram never enters. **They do not remove the degeneracy; they keep it out of the denominator.** GeoVac inverts because a block-encoding wants a standard Hermitian eigenproblem -- so **the conditioning exposure is created by the ENCODING REQUIREMENT, not by the basis and not by the metric.** Captured in Paper 60 `sec:obstruction` with a verified bibitem.

**The lever, priced honestly and left open.** A determinantal re-posing would remove the conditioning multiplier but reinstates the outer nonlinear search over the scale that `eq:secular` exists to eliminate. Which is cheaper at scale is an unasked, answerable resource question.

**A self-catch inside the fix.** The retraction was first written as a new item at the foot of the memory file while the original claim stayed live at its own locus 170 lines above -- the corrected-in-the-owner-left-standing-in-the-citer class, committed inside the edit that was fixing it. Both loci now carry the marker, and the file's `description:` line, which drives recall, was itself asserting the retracted half and is corrected.

**Two owed threads.** (i) The method is credited in the later literature to **Novosadov**, not to Monkhorst-Jeziorski -- unverified, and an untouched Russian-language thread on overcompleteness continuing to 1986. Deliberately NOT cited in the paper, since an inline attribution with no bibitem is a C20 defect. (ii) The lineage obtains three-centre **one-electron** integrals with no three-centre evaluation, via a completeness sum over Fock-sphere states. **Not a breach of the three-centre wall** -- this corpus's wall is the two-electron ERI, and the source concedes the advertisement weakens at three or more centres -- but a diagnostic probe candidate.

"""

MEMO_ANCHOR = "\n---\n\n## 8. Follow-on items (2026-09-12, PI-directed)\n"
MEMO_MARKER = "## 7f. Monkhorst-Jeziorski"
MEMO_NEW = """
## 7f. Monkhorst-Jeziorski 1979 -- verdict, and the resolution

**SECONDARY-QUOTED.** Abstract + record verified (Crossref, re-verified by the
PM); two-page body closed, zero repository copies, UNREAD. Correct DOI
`10.1063/1.438337`, JCP **71**(12), 5268-5269. Memo
`debug/lit_scan/monkhorst_jeziorski_1979_memo.md`.

| suspended claim | verdict |
|:--|:--|
| "no prior art for secular-matrix norm growth" | **SURVIVES**, narrowed; rests on abstract + ref list + page count |
| "conditioning is a blind spot in the whole Avery canon" | **RETRACTED as phrased** |

**The predicted resolution (Sec. 7e) was WRONG.** The PM expected an
`L^2`-vs-`V_0` metric distinction. It is the *same pencil*; the degeneracy is in
their matrix too. The difference is that they never **invert**: the overlap
enters only multiplicatively, so a near-null direction gives a spurious branch
rather than amplified error, and the basis is exactly orthonormal in the metric
used. **They keep the degeneracy out of the denominator.** GeoVac inverts
because a block-encoding wants a standard Hermitian eigenproblem -- so the
conditioning exposure belongs to the ENCODING REQUIREMENT, not the basis or the
metric. This is a better answer than the one that was predicted, and it is the
single most useful thing the whole prior-art arc produced.

**Open lever:** determinantal re-posing removes the conditioning multiplier and
reinstates the outer nonlinear scale search `eq:secular` exists to eliminate.
Unpriced.

**Impact set (Sec. 7e) verified against the verdict:** the memory items were the
real exposure and both moved; Paper 60 `sec:quantum`'s novelty claim was
insulated exactly as predicted (it is scoped to *quantum* algorithms, which a
1979 classical note cannot reach); `sec:obstruction` gained the resolution
rather than needing a scoping clause. Pre-enumeration held.
"""


def main() -> int:
    n = 0
    with open("CHANGELOG.md", encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) >= 1:
        with open("CHANGELOG.md", "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, CL_NEW + CL_ANCHOR, 1))
        n += 1
        print("  ok   CHANGELOG.md")
    with open("debug/sprint_contraction_seam_memo.md", encoding="utf-8") as fh:
        m = fh.read()
    if MEMO_MARKER not in m and m.count(MEMO_ANCHOR) == 1:
        with open("debug/sprint_contraction_seam_memo.md", "w",
                  encoding="utf-8") as fh:
            fh.write(m.replace(MEMO_ANCHOR, MEMO_NEW + MEMO_ANCHOR))
        n += 1
        print("  ok   sprint memo")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
