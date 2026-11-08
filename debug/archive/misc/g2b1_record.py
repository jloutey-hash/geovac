"""Record group2 baseline FULL run — Batch 1 (Papers 8, 11, 12). v5.11.14. Idempotent."""
from __future__ import annotations
import sys, pathlib

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.13] - 2026-09-13\n"
MARKER = "## [v5.11.14]"
ENTRY = """## [v5.11.14] - 2026-09-13

**`/qa group2` baseline FULL run — Batch 1 (guardrail Paper 8 + Papers 11, 12) = FAIL, remediated.** Five reviewers (3 code, 1 chunked claims, 1 chunked citations), tree frozen. The baseline's payoff is a real pre-existing LARGE that no delta would have reached.

### The LARGE — a guardrail negative overstated at its summary surfaces (PI-approved fix)

The **Sturmian Structural Theorem** (Papers 8-9 guardrail) was scoped down in v4.72.x (PI-reviewed): the no-go holds only for the single-n / shared-p0 / block-diagonal APPROXIMATION; the genuine cross-n Coulomb-Sturmian overlap couples n, escapes the theorem, and binds (supplying the R-dependence algebraically via Shibuya-Wulfman). That correction reached the paper bodies, CLAUDE.md §3.5, and `paper_fci_molecules` — but NOT the summary surfaces. **Two reviewers independently (claims chunk + code Paper 8) caught the pre-correction, over-strong reading still live at Paper 8's abstract, section-intro, and conclusion — the conclusion asserting "binding reintroduces continuous spatial geometry (Corollary cor:binding)" while that corollary says the opposite — and at four loci in Paper 11.** This is exactly the Summary-Surface Reading Rule class. **PI-confirmed 2026-09-13: propagate the approved scope.** Seven loci rewritten (3 in Paper 8, 4 in Paper 11) mirroring the wording already approved in `paper_fci_molecules`; the self-contradiction is resolved.

### MATERIAL — Kolos-Wolniewicz reference value (Paper 12 citations)

The exact-H2 benchmark was printed `-1.17475`/`0.17475` (dropped digit) at three loci; true KW 1968 value is `-1.174475`/`0.174475`. The D_e% column already used the correct denominator (the numerical column 79.6/80.1 matches 0.174475), so only the printed literal was corrected; the column and 92.4% headline are the paper's own computed values, untouched.

### Verdicts and the run-#1 regression class

Per-paper: **Paper 8 code FAIL** (the LARGE + one SMALL/BUG — a backing test NaN'd, see below); **Paper 11 code PASS** (216 tests, 0 failed); **Paper 12 code PASS** (32 tests); **claims chunk FAIL** (the LARGE + NITs); **citations chunk FAIL** (the KW MATERIAL + SMALLs). **Crucially, the 2026-06-26 run-#1 "headlines backed by DELETED code" regression is CLOSED for both at-risk papers:** Paper 11's spectral solver was deleted in v2.7.0 and restored 2026-06-27 (now live, 168 core tests); Paper 12's "no integrals" contradiction is remediated (the paper now names its one residual B_l quadrature in every load-bearing surface). No current deleted-code regression.

### Also remediated

- **A failing backing test:** `test_paper8_overlap_cross_n::test_overlap_conserves_m` NaN'd under `--slow` (m!=0 overlap quadrature broken near the coordinate corner; the test was also vacuous — it never built the cross-m element its docstring claimed). Non-load-bearing (the paper uses only m=0 overlaps). Rewritten onto the m=0 path with the m!=0 limitation documented; now 3 passed. **Caught only because the reviewer corrected its own earlier misread of `tail`'s exit code for pytest's** (the never-pipe-verification rule).
- **Two stale claim-matrix rows:** Paper 8's Structural Theorem (still NO-TEST) and Paper 12's full-V_ee/92.4% (still NO-TEST/tautological-B_l) both moved to BACKED-SOUND — the backing was backfilled and is now genuine (independent anchors, two-route V_ee exactness).
- **NITs:** Paper 12 abstract "92.4% with 27" -> 92.2% (92.4% is the N=72 plateau); Paper 11 abstract "computed algebraically, eliminating quadrature entirely" -> "can be computed" (the default path is Gauss-Laguerre quadrature); Paper 8 Fiedler "gives the same object" -> "a discrete analog"; Herbst-Avery-Dreuw scope softened (atomic CS-HF, not molecular binding).

### Carried to the group review (declared NITs, not blocking)

Three orphaned bibitems (`Boys1970` P8, `James1933`/`loutey_fci`/`loutey_paper10` P11); a bare uncited He `-2.9037` in Paper 8; Paper 11's `<0.1% one-electron atoms` is trunk-dependent (C7) and awaits the trunk re-cert; Paper 11 R_eq 2.005-vs-2.001 signposting and a headline inline-provenance tag; the "single 1D quadrature" wording (defensible as-is).

### Deterministic + verdict

Whole-group deterministic layer 14/14 (established v5.11.13). Batch-1 papers: post-remediation deterministic PASS, all 13 group2 papers compile. **Batch 1 verdict: FAIL, remediated** — 1 LARGE (guardrail, PI-approved fix applied), 1 MATERIAL (KW literal), the rest SMALL/NIT; zero fabricated results, zero live deleted-code regressions. **Batches 2 (Papers 13/15/17), 3 (Paper 19 + FCI), 4 (Paper 58 + synthesis C9 + completeness) remain for the baseline.**

"""


def main() -> int:
    p = pathlib.Path(CL); t = p.read_text(encoding="utf-8")
    if MARKER in t:
        print("skip"); return 0
    if t.count(ANCHOR) != 1:
        print(f"MISS anchor count={t.count(ANCHOR)}"); return 3
    p.write_text(t.replace(ANCHOR, ENTRY + ANCHOR), encoding="utf-8")
    print("ok CHANGELOG v5.11.14")
    c = pathlib.Path("CLAUDE.md"); ct = c.read_text(encoding="utf-8")
    ct = ct.replace("**Version:** v5.11.13 (September 13, 2026)", "**Version:** v5.11.14 (September 13, 2026)", 1)
    c.write_text(ct, encoding="utf-8"); print("CLAUDE.md v5.11.14")
    return 0


if __name__ == "__main__":
    sys.exit(main())
