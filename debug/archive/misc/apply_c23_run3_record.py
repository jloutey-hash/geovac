"""Record C23 run #3 in the CHANGELOG entry and the sprint memo. Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
MEMO = "debug/sprint_contraction_seam_memo.md"

CL_ANCHOR = "### Owed (PI items)\n"
CL_MARKER = "### C23 run #3"
CL_NEW = """### C23 run #3 (at-authorship): the tagging claim was audited the day it was written, and it needed correcting

The at-authorship trigger adopted in v5.11.2 fired on the `[SYMBOLIC + MEASURED]` paragraph written earlier the same day. Memo `debug/lit_scan/c23_run_003_m2_tagging_memo.md`. Three verdicts, and the scan's value was not the attributions:

- **T1 (`c_1 = pi^2` read as a Dirichlet eigenvalue): PRIOR ART**, in the source the paper ALREADY CITES. Boettcher-Widom present `c_alpha` as the least eigenvalue of `(-1)^a u^(2a) = lam u` on `[0,1]` with clamped ends; at `a = 1` that is the Dirichlet problem. Verified here: bibitem present, title reads "...to higher-order **Wirtinger-Sobolev** inequalities".
- **T2 (the band-limited second-moment minimum): PRIOR ART**, and the priority risk resolved opposite to the one flagged — it is NOT Slepian/Landau/Pollak (their functional is energy concentration, eigenvalues transcendental in the time-bandwidth product, never `pi^2`) but the second-moment/**Wirtinger** branch. Which means the "no Bessel present" measurement is the `b == 1` case of the same theorem: a change of REPRESENTATION, not an independent route. The independent-route claim is withdrawn.
- **T3 (the two-sided taxonomy itself): ABSENT** as a taxonomy — but two of its three legs are someone else's, and its corollary was wrong.

**The finding that mattered was a defect in prose written hours earlier**, not an attribution. The paragraph's operational corollary — truncation-side prices are matrix-level and reachable, continuum-side prices are symbol-level and untouchable — is **false in both halves**, and the withdrawal is recorded above. Verified locally before editing (line 1084: the matching polynomial is chosen to *share the symbol's zero*, per Serra). Corrected at the owner and swept to all three dependents that restated it (claim matrix, this entry, the sprint memo).

**One reported defect did NOT hold and was not "fixed".** The scan flagged the inline `c_1 = pi^2` attribution as a C20 bibitem-less case; checked, and the bibitem exists and the sentence cites it — the flagged phrase is a second mention inside the same sentence. C20 passes in scope. **One did hold:** `bottcher_widom2005` was cited arXiv-only; published coordinates added after Crossref verification (Birkhaeuser, 2007, pp. 73-87, doi:10.1007/978-3-7643-7980-3_4). The series volume number was NOT added — it is search-level only.

**Flagged, not resolved:** `debug/lit_scan/toeplitz_finite_section_memo.md` records BOTH "the conditioning and the non-locality are one fact" and "they are two independent facts" at different points. The corrected paragraph adopts the second and now says so.

"""

MEMO_ANCHOR = "\n---\n\n## 8. Follow-on items (2026-09-12, PI-directed)\n"
MEMO_MARKER = "## 7b. C23 run #3"
MEMO_NEW = """
## 7b. C23 run #3 -- the at-authorship trigger paid for itself immediately

Memo `debug/lit_scan/c23_run_003_m2_tagging_memo.md`. Run on the tagging
paragraph the day it was written, per the v5.11.2 scope change.

| claim | verdict |
|:--|:--|
| T1 `c_1 = pi^2` as a Dirichlet eigenvalue | **PRIOR ART** -- Boettcher-Widom, the source already cited |
| T2 band-limited second-moment minimum | **PRIOR ART** -- Wirtinger branch, NOT Slepian/Landau/Pollak |
| T3 the two-sided taxonomy | **ABSENT** as a taxonomy, but its corollary was wrong |

**The scan's value was a defect in prose written hours earlier.** The
removability corollary was false in both halves (see Sec. 3, withdrawn).
Verified locally before editing. Corrected at the owner and swept to all three
dependents.

One reported defect did NOT hold: the inline `c_1` attribution was flagged as
C20 bibitem-less, but the bibitem exists and the sentence cites it. Not
"fixed". One did: `bottcher_widom2005` was arXiv-only; published coordinates
added after Crossref verification, **without** the series volume number, which
was search-level only.

*Open direction the scan could not close:* Serra-Capizzano's *Practical Band
Toeplitz Preconditioning and Boundary Layer Effects* was unreachable (Springer
IDP redirect) and is the source most likely to state a truncation-vs-symbol
separation in the literature's own words. T3's ABSENT verdict should be re-run
against it before being leaned on.
"""


def main() -> int:
    n = 0
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t:
        if t.count(CL_ANCHOR) < 1:
            print("CL anchor missing")
            return 2
        t = t.replace(CL_ANCHOR, CL_NEW + CL_ANCHOR, 1)
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(t)
        n += 1
        print("  ok   CHANGELOG.md")
    with open(MEMO, encoding="utf-8") as fh:
        m = fh.read()
    if MEMO_MARKER not in m:
        if m.count(MEMO_ANCHOR) != 1:
            print(f"memo anchor count={m.count(MEMO_ANCHOR)}")
            return 3
        m = m.replace(MEMO_ANCHOR, MEMO_NEW + MEMO_ANCHOR)
        with open(MEMO, "w", encoding="utf-8") as fh:
            fh.write(m)
        n += 1
        print("  ok   sprint memo")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
