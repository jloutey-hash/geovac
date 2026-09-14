"""Record the Paper 61 DELTA-verification (v5.11.12). Idempotent."""
from __future__ import annotations
import sys

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.11] - 2026-09-13\n"
MARKER = "## [v5.11.12]"
ENTRY = """## [v5.11.12] - 2026-09-13

**`/qa paper_61` DELTA-verification = DEFECTS, remediated. The paper's own mathematics and prose are CLEAN; both defects sat outside it.** Four dimensions over the seam scope (Paper 59 + Paper 61 + both syntheses), unseeded, tree frozen at dispatch. Deterministic layer 14/14. This is the confirming clean-delta the PI-ratified DoD requires before Paper 61's FULL run; the paper is unchanged since its 2026-09-08 remediation, so the delta verified that remediation held corpus-wide.

### Dimension results

| dimension | verdict |
|---|---|
| deterministic (14 gates) | PASS |
| claims / prose (C3/C5/C6/C8) | **CLEAN-DELTA** — 0 defects; all 4 fixes-of-fixes from the 09-08 same-day delta hold; every Sp4 statement is monodromy-in-Sp4(Z)/Galois-in-Sp4(C); W0=pi^2/rho^2 [SYMBOLIC] derivation stands alone; 20+ honest-scope sentences all point the limiting direction |
| code / test-backing (C1-C2, --slow) | **CLEAN-DELTA** — 75 passed / 0 skipped / 0 failed, reproduced twice; all four DoD branch criteria hold; the W0 [OPEN]->[SYMBOLIC] upgrade fire-tested (excludes W0=1); no restricted-evaluation false positive |
| claim-impact | **DEFECTS** (1 NIT) |
| external citations (C4) | **DEFECTS** (1 MATERIAL) |

### The 09-08 remediation held

All three REGISTERED retractions verified corpus-wide with no survivor: `Gal in Sp4(Z)` (now monodromy/Galois-correct in every paper, doc, test, driver, synthesis, and the paper_59.done.md goalpost that had certified it), "transcendence cancels in the determinant" (W0 now [SYMBOLIC]), and "every CM fibre" (all four previously-live loci including the self-reseeding memory file). Both C7 cross-paper edges (Paper 35/WH7, Paper 56 seam T-2) are non-stale.

### The two defects, both outside the paper's prose, both remediated

1. **Citations MATERIAL — Chowla-Selberg had no bibitem.** The CM-period Gamma-value attribution is made by name at three load-bearing loci (the abstract [MEASURED] claim, L117, L182 [MEASURED]) with no `\\bibitem` and no `\\cite` anywhere — the bibitem-free layer the C20 gate cannot see. The attribution is CORRECT. Added `chowla_selberg1967` (A. Selberg and S. Chowla, J. Reine Angew. Math. 227, 86 (1967), author order and details verified at de Gruyter/EUDML) and cited it at the two body loci. All 17 other bibitems were verified GROUNDED at source; graph balanced; every DoD-flagged distinction (K(1/2)/lemniscate, Sp4(Z)/(C), Broadhurst-Mellit determinant vs Broadhurst-Roberts quadratic) correct in the paper.

2. **Claim-impact NIT — the retired K(1/2)=lemniscate naming survived in a tracked driver.** `debug/routeC_pslq_fit.py:6` labelled K(1/2)=varpi "the lemniscate constant" (the classical lemniscate constant is sqrt(2)*varpi, a different number). The paper fixed this in its own prose, but retired claim #6 had been given NO `check_retracted_terms.py` entry (a paper-prose-only fix), so it was never swept corpus-wide -- the locus-by-locus-remediation class. Fixed the comment AND added the missing C16 entry `p61-k12-is-lemniscate`, fire-tested both directions: it FIRES on the planted conflation and stays SILENT on every correct usage (the paper's sqrt2 distinction, the drivers' true-lemniscate definition, the done.md record, the lit-memo negatives). The exemption was tightened after a first version wrongly exempted the conflation via the `sqrt(2)` in a neighbouring disc-8-period definition.

### Verdict and what remains

**DEFECTS**, remediated. A delta that found defects is not a clean delta, so per the ratified DoD the FULL certifying run stays locked. The remaining surface is now three mechanical fixes (a bibitem+2 cites, a driver comment, a fire-tested gate entry); a confirming clean delta over them -- foldable into the group3 review -- then the FULL. **The paper itself passed claims and code cleanly**, so the distance to certification is citation-completeness and gate-hygiene, not a mathematical question.

Pre-existing declared items surfaced for the FULL run (not delta regressions): the 64-digit PSLQ negative has no pytest (test-backed part is the narrower disc-4 wt<=2 negative); `test_paper59_diagonal_A` PSLQ leg has a decoy but no positive control; T2 digits 22-66 are quoted from a permanent artifact, not recomputed; `docs/claims_register.md` has no row for Papers 54-61.

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
    print("ok CHANGELOG v5.11.12")
    return 0


if __name__ == "__main__":
    sys.exit(main())
