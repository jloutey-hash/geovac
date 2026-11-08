"""Propagate the primaries-scan corrections to every dependent record.

Three claims moved in Paper 60 (their owner); the documents whose argument
rests on them move here, per CLAUDE.md Sec. 9's retraction->dependents rule.

Also records the two owed items the scan surfaced, one of which threatens a
STANDING corpus claim.
"""
from __future__ import annotations

import sys

EDITS = [
    # ---- claim matrix: the weight row is no longer "ours" ----
    ("docs/claim_test_matrix.md",
     "| 60 | §molecular [MEASURED] — the conditioning law is carried by the "
     "TRANSLATION, not by the metric that carries it:",
     "| 60 | §molecular [MEASURED + PRIOR ART] — the conditioning law is "
     "carried by the TRANSLATION, not by the metric that carries it "
     "(**re-tiered 2026-09-12: this is a KNOWN theorem, in a stronger form — "
     "Ahmad, Al-Aidarous, Alrehaili, Ekström, Furci & Serra-Capizzano, *Numer. "
     "Algorithms* **78**(3), 867–893 (2018), whose Eq. (22) gives the "
     "preconditioned-pencil eigenvalues EXACTLY for the τ/DST-I algebra — i.e. "
     "Toeplitz minus Hankel, this paper's own structure — with the weight "
     "cancelling identically. The C23 run #2 verdict of ABSENT was reached "
     "without this source and is superseded; the measurement confirms a "
     "theorem rather than establishing one**):"),

    ("docs/claim_test_matrix.md",
     "Whether any position-space `V_0` preserves the constant is OPEN. "
     "rests on: eq:sigma_law |",
     "Whether any position-space `V_0` preserves the constant is OPEN. "
     "**Hypothesis gap, stated flatly:** our chirp symbol is not a "
     "trigonometric polynomial and its ratio has infinitely many sign changes, "
     "so Ahmad et al.'s Theorem 1 hypotheses FAIL for it; only their τ-algebra "
     "identity and the L¹ localization results apply unconditionally. "
     "rests on: eq:sigma_law |"),

    # ---- CHANGELOG ----
    ("CHANGELOG.md",
     "### The one new result: the law is carried by the TRANSLATION, not the metric",
     "### The law is carried by the TRANSLATION, not the metric -- and it is a known theorem"),

    ("CHANGELOG.md",
     "**C3 (the weight-independence above): ABSENT** -- the one thing here that is ours.",
     "**C3 (the weight-independence above): ABSENT** -- which was the verdict "
     "*on the sources run #2 reached*, and **it did not survive the primaries "
     "pass the same day** (below)."),
]

CL_ANCHOR = "### Owed (PI items)\n"
CL_MARKER = "### Primaries pass"
CL_NEW = """### Primaries pass: two sources two earlier scans could not reach, and both moved a claim

Memo `debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md`. What unlocked it was route, not persistence alone: direct `curl` plus local `pdftotext` reached publisher abstracts and PDFs that the summarising fetch path 403s on.

**1. The translation identification is NOT ours, on three counts.** Shibuya and Wulfman's own 1965 abstract builds the molecular `p0` operator from the united-atom one "by a sum of unitary transformations, one for each nucleus in the molecule" -- one unitary per centre, the structure in substance. Wulfman & Takahata gave the explicit continuous-group formulation two years later. Weatherford & Red titled 2002-2004 papers on representing that operator in a Coulomb-Sturmian basis, and even plot against "the translation distance multiplied by the screening parameter" -- our own `kR`. **What survives as ours is the SYMBOL**: that in the sine basis on the Fock polar angle the operator is the finite section of multiplication by `j0(kR cot(chi/2))`. Corrected in place. *Provenance, stated because it matters:* the Royal Society page is 403 from here as it was for two prior scans, so the abstract quotation is RELAYED from the scan's direct read; the paper BODY is still unread. The correction is conservative either way -- it gives a claim away.

**2. The weight-independence result is prior art, in a stronger form.** Ahmad, Al-Aidarous, Alrehaili, Ekstroem, Furci & Serra-Capizzano, *Numer. Algorithms* **78**(3), 867-893 (2018) -- Crossref-verified here for title, six authors, volume, issue, pages and year -- give the preconditioned-pencil eigenvalues in almost closed form. Their Eq. (22) is an **identity** for the tau (DST-I) algebra, i.e. Toeplitz minus Hankel, which is exactly this paper's structure; the weight cancels identically. So GeoVac sits closer to the exactly-solvable core than to the asymptotic Toeplitz statement. Re-tiered from `[MEASURED]` to `[MEASURED + PRIOR ART]`. **Hypothesis gap stated flatly:** our chirp symbol is not a trigonometric polynomial and its ratio has infinitely many sign changes, so their Theorem 1 hypotheses fail for it; only the tau identity and the L^1 localization results apply unconditionally.

**3. `eq:sigma_law`'s unexplained ~1% residue is mostly the grid convention.** The tau statement places eigenvalues at `j pi/(n+1)`; collapsing with `(n+1)^2` moves the residue at `n=160, kR=2` from **-0.99% to +0.26%**, a fourfold reduction with a sign flip. Re-measured here across `n = 40..320` before editing. A genuine `O(1/n)` term survives in both conventions, so the sentence says "dominated by", not "is".

**Two reported defects did NOT hold** and were not "fixed": the `avery2004` bibitem is already correct at vol. **100**, 121 (2004), and the inline `c_1` attribution does carry its `\\cite`. **One bibliographic correction recorded:** the Serra-Capizzano paper the corpus declined to cite as "LAA 270 (1998)" is the wrong reference entirely -- that is Tyrtyshnikov-Zamarashkin; the intended ones are LAA **267** (1997) 139-161 and LAA **282** (1998) 161-183.

"""

MEMO_ANCHOR = "\n---\n\n## 8. Follow-on items (2026-09-12, PI-directed)\n"
MEMO_MARKER = "## 7c. Primaries pass"
MEMO_NEW = """
## 7c. Primaries pass -- both targets moved a claim

Memo `debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md`.

| target | verdict | consequence |
|:--|:--|:--|
| Shibuya-Wulfman 1965 | **SECONDARY-QUOTED** (abstract + full reference list; body still unread) | the **translation** identification is prior art on three counts; only the **symbol** survives as ours |
| Toeplitz pencil / ratio symbol | **REACHED** | Ahmad et al., *Numer. Algorithms* **78**(3) 867-893 (2018), Eq. (22): an identity for the tau/DST-I algebra. Our weight-independence result is prior art, re-tiered |

Plus: `eq:sigma_law`'s ~1% residue is **dominated by the `n -> n+1` grid
convention** (-0.99% -> +0.26% at `n=160`, `kR=2`), re-measured locally before
editing; an `O(1/n)` term survives in both conventions.

**Two reported defects did not hold** (the `avery2004` volume is already right;
the inline `c_1` attribution does cite). Recording this because the pattern
across three scans is consistent: **the attribution findings are reliable, the
"live citation defect" findings are roughly half right**, and every one must be
checked against the current file before acting.

### OWED, and it threatens a standing claim

**Monkhorst & Jeziorski, "No linear dependence or many-center integral problems
in momentum space quantum chemistry", J. Chem. Phys. 71, 5268 (1979).**
Unopened. A 1979 title of that exact shape bears directly on
`memory/avery_method_and_prior_art_gaps.md` items 2 and 3 -- "there is NO prior
art for secular-matrix norm growth" and "conditioning is a blind spot in the
whole Avery canon" -- which are load-bearing for Paper 60's novelty framing and
have already survived three scans. **Read it before repeating either claim.**

Also owed: Serra-Capizzano's *Practical Band Toeplitz Preconditioning and
Boundary Layer Effects* (Springer IDP redirect, unreachable), the source most
likely to state a truncation-vs-symbol separation in the literature's own
words. C23 run #3's T3 ABSENT verdict should be re-run against it.
"""


def main() -> int:
    applied = 0
    for path, old, new in EDITS:
        with open(path, encoding="utf-8") as fh:
            t = fh.read()
        if new[:60] in t:
            continue
        if t.count(old) != 1:
            print(f"  MISS {path}: count={t.count(old)} :: {old[:70]!r}")
            continue
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new))
        applied += 1
        print(f"  ok   {path}")

    with open("CHANGELOG.md", encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) >= 1:
        with open("CHANGELOG.md", "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, CL_NEW + CL_ANCHOR, 1))
        applied += 1
        print("  ok   CHANGELOG.md (primaries section)")

    with open("debug/sprint_contraction_seam_memo.md", encoding="utf-8") as fh:
        m = fh.read()
    if MEMO_MARKER not in m and m.count(MEMO_ANCHOR) == 1:
        with open("debug/sprint_contraction_seam_memo.md", "w",
                  encoding="utf-8") as fh:
            fh.write(m.replace(MEMO_ANCHOR, MEMO_NEW + MEMO_ANCHOR))
        applied += 1
        print("  ok   sprint memo (primaries section)")
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
