"""Record group2 baseline FULL run — Batch 2 (Papers 13, 15, 17). v5.11.15. Idempotent."""
from __future__ import annotations
import sys, pathlib

CL = "CHANGELOG.md"
ANCHOR = "## [v5.11.14] - 2026-09-13\n"
MARKER = "## [v5.11.15]"
ENTRY = """## [v5.11.15] - 2026-09-13

**`/qa group2` baseline FULL run — Batch 2 (Papers 13, 15, 17, the natural-geometry hierarchy) = FAIL, remediated. Much cleaner than Batch 1: no LARGE.** Five reviewers, tree frozen. Every code dimension PASSED and all headline numbers reproduce from current code.

### The three run-#1 regressions are all CLOSED

- **Paper 15 spectral speedups** (16× radial / 20× / 269× angular) were "backed by v2.7.0-deleted code" in run #1; the radial solver was restored from the pre-deletion commit and the angular test un-archived — both now exist, run, and reproduce (128 tests, 0 failed).
- **Paper 15 "exceeds Paper 12"** was a run-#1 FALSE-POSITIVE (passing via the disavowed adiabatic solver); now genuinely backed by the 2D variational solver, whose <100% band actively rejects the adiabatic 105% artifact.
- **Paper 17 BeH₂ 11.7%** was flagged stale in run #1 (→19.7% under a PK-blind bug); the bug was fixed 2026-06-27 and 2.80 bohr / 11.7% now reproduces exactly. Every Paper-17 composed R_eq was re-run: no drift.

### Remediated

- **H₂O reconciled (paper right, docs stale).** The code reviewer confirmed H₂O 19.4% (R_eq 1.459) reproduces from a passing test. The claim-matrix, claims-register, and group2 DoD still said **26%**; all three corrected to 19.4%. The paper was already current.
- **A self-contradiction in Paper 17 §VIII.B:** "molecular equilibrium exists without PK ... 5× overbinding" is the ADIABATIC reading, contradicting the paper's own authoritative 2D-variational result (UNBOUND, D_e<0). Reconciled: the minimum is a feature of the adiabatic PES, not a bound molecule.
- **Citations (the FAIL driver, all SMALL):** a genuine **WRONG-ID** — `Mitnik2021` cited to Comput. Phys. Commun. 269/108145, which resolves to an unrelated QCD code; corrected to the verified Mol. Phys. 119(8), e1881179 (2021), Mitnik/López/Ancarani (arXiv:2006.06616), third author Gasaneo→López. A **MISATTRIBUTION** — a 3.29 doubly-excited ¹S autoionization-width ratio credited to Madden–Codling, who measured the dipole-allowed ¹Pᵒ series; reframed. And the same dropped-digit class as Batch 1 — Paper 15's exact-H₂ D_e printed 0.17447 (truncation of 0.174475); corrected.
- **Prose:** Paper 15's abstract now states the 96.0% includes a ~1pp Schwartz cusp correction (the pure-variational value ~95% already exceeds 92.4%); Paper 17's "first realization of a category ... morphisms" (an unverified category-theory novelty claim) reduced to the fiber-bundle framing the paper already proves.
- **Two claim-matrix tiers upgraded:** the l_max=2/structural-divergence keystone (was print-only NO-TEST → BACKED-SOUND, `TestLmaxDivergenceMonotone` reproduces the monotone drift) and H₂O (NO-TEST/26% → BACKED-SOUND/19.4%).

### Coverage gaps — two CI-infeasible headlines, disclosed not fixed

Papers 13 (graph-native 0.19% at n_max=7) and 15 (96.0% D_e at l_max=6) each carry a headline that is too expensive to test at its stated truncation (l_max=6 ≈ 754 s/point; n_max=7 CI-uncomputable). Both are on monotone converging sequences whose lower anchors ARE tested (l_max=4 = 94.3%; n_max=5 = 0.25%) and both are already disclosed NO-TEST in the matrix. The code Paper 15 reviewer flagged the 96.0% "LARGE by principle" but well-mitigated; recorded as a standing disclosed coverage gap, not a defect. **Owed for the eventual cert:** an in-paper NO-TEST/extrapolation footnote on each (Paper 13's 0.19% should read as an extrapolation from the 0.25%@n=5 anchor, which the 2026-08-29 ERI fix moved it to).

### Carried to the group review (declared NITs)

Seven orphaned bibitems (KlarKlar1980 P13; macek1968/james1933/morse1953 P15; paper12/paper16/fci_a P17-internal); citation metadata slips (Kereselidze2016 title/authors, Abdouraman2016 initials, bishop1977 citekey year, huber1979 "ab initio" vs experimental caption); Paper 13's S³-vs-hydrogenic integral provenance and its 0.05%-over-0.022% table billing; Paper 15's l_max=2-jump endpoint numbers; the speedup multipliers stated flat (want "measured, hardware-dependent"); Paper 17 tab:pk B-value drift (7.00→6.80) and loose-range regression guards; stale test docstrings.

### Verdict

**Batch 2: FAIL, remediated** — driven by SMALL citation defects (one wrong-ID, one misattribution) plus doc staleness and one prose contradiction; **zero LARGE, zero live deleted-code regressions, every headline number reproduces.** Deterministic layer green, all 13 group2 papers compile. **Batches 3 (Paper 19 + FCI-atoms + FCI-molecules) and 4 (Paper 58 + synthesis C9 + completeness-critic) remain.**

"""


def main() -> int:
    p = pathlib.Path(CL); t = p.read_text(encoding="utf-8")
    if MARKER in t:
        print("skip"); return 0
    if t.count(ANCHOR) != 1:
        print(f"MISS anchor count={t.count(ANCHOR)}"); return 3
    p.write_text(t.replace(ANCHOR, ENTRY + ANCHOR), encoding="utf-8")
    print("ok CHANGELOG v5.11.15")
    c = pathlib.Path("CLAUDE.md"); ct = c.read_text(encoding="utf-8")
    ct = ct.replace("**Version:** v5.11.14 (September 13, 2026)", "**Version:** v5.11.15 (September 13, 2026)", 1)
    c.write_text(ct, encoding="utf-8"); print("CLAUDE.md v5.11.15")
    return 0


if __name__ == "__main__":
    sys.exit(main())
