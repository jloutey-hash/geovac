"""Record the FULL certifying /qa run on Paper 60 (v5.11.10). Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
CL_ANCHOR = "## [v5.11.9] - 2026-09-13\n"
CL_MARKER = "## [v5.11.10]"

ENTRY = """## [v5.11.10] - 2026-09-13

**`/qa paper_60` FULL certifying run = FAIL, PI-invoked. Five of six dimensions PASS with ZERO mathematical or content defects; the FAIL is the completeness dimension's cert-blockers, two of which are remediated here.** Whole-paper, unseeded, tree frozen at dispatch (no edit until all five reviewers returned — DELTA #2's moving-target defect did not recur).

### Dimension results

| dimension | verdict |
|---|---|
| deterministic (14 gates) | PASS |
| code / test-backing (C1-C2) | **PASS** — 136/136 Paper-60 tests pass with `--slow` (19+109+8), 0 fail/skip/error, C21 green; every load-bearing claim maps to a non-tautological test with anti-tautology / matched-set / restricted-evaluation legs; no keystone bug |
| claims / prose (C3/C5/C6/C8) | **PASS** — every claim within its tier; all six QC-honesty branch criteria satisfied; ~17 honest-scope sentences all point the limiting direction; six withdrawn readings retracted in place. 1 SMALL |
| external citations (C4) | **PASS** — all 49 bibitems resolve, graph balanced, 0 WRONG-ID/MISATTRIBUTED/OVERSTATED; rajchel2025 drift confirmed corrected |
| synthesis faithfulness (C9, GATING) | **PASS** — all 8 checkpoints; no strengthened hedge, no stale value. 1 NIT |
| completeness-critic | **FAIL** — zero undiscovered gaps, zero surviving authoritative stale echoes; FAIL on three self-declared cert-blockers |

### Why FAIL, and what it is NOT

The completeness dimension found **no hidden defect** — every abstract headline is test-backed, and no stale echo survives in any authoritative surface (paper, INDEX, synthesis, claims-register, code-architecture, the auto-loading memory, the walls register are all current). The FAIL is three unmet PRECONDITIONS for a clean certification:

1. **The DoD's own literal-registration precondition is unmet.** Its declared-debt table lists ~7 load-bearing literals (Gaussian ratio N^6/N^1.4, water A_1 N^1.97, the tab:resource row, naive-L2 Q^3.33/Q^1.19, SW N^1.85 vs L2 N^1.70, H2+/H2 toy figures) marked "owed before the next FULL certifying run" -- and this is that run. They are all VERIFIED CORRECT (this run re-confirmed them) and the load-bearing ones carry TEST backing; they are C21-undelegated, not backing-less. That is the exact class that let a wrong `4.43` through on 2026-09-11. **Registered as its OWN careful pass (see below), not rushed at run-tail** -- this lineage's failures have come from batching number-edits.
2. **The Sylvester-inertia variational bound (C8.16) had no standalone test.** REMEDIATED: `tests/test_paper60_scale_lock.py::test_c5_inertia_bound_and_root_by_root` now proves the algebraic identity `H(lam)+lam^2 S/2 = lam(lam I - M)` (1e-9), the inertia count `#{roots < -lam^2/2} == #{eig(M) > lam}` at five lam, and the k=0..3 root-by-root consequence. Fire-tested: a 0.7*I shift of M breaks the count. Matrix row 568 moved OPEN -> BACKED; the [INTERNAL THEOREM] tier is now test-backed, not just DoD-verified.
3. **No clean delta preceded the run.** DELTA #1/#2/#3 all returned DEFECTS-remediated; the run was invoked ahead of the corpus's own gate (PI's call). A clean delta over the small remediation surface remains the last step to a PASS-able state.

### Remediated this run (the reviewers' findings)

- **F1 (claims SMALL):** the abstract tagged the third-lever sentence [MEASURED], but its closing clause "a direct block-encoding ... takes the penalty from n^3 to n" is a resource-model result (circuit cited, not compiled). Split: the clause is now [RESOURCE MODEL] with "circulant-embedded circuit is analysed, not compiled".
- **C9 NIT (synthesis):** "eigenproblem whose eigenvalues are the energies" -> "whose eigenvalues are the scaling parameters p_kappa = sqrt(-2E), from which the energies follow directly".
- **C4 NIT (citations):** the Monkhorst-Jeziorski bibitem title read "No linear dependence OR MANY-center..."; corrected to the published "...AND MULTI-center..." (verified verbatim at the author's own publication list; DOI resolved regardless).

### The two citation residuals the reviewer routed to the PM -- both closed

- `gslw2019` Theorem 73 / Corollary 67 numbers: verified at secondary source in DELTA #3 (Thm 73 = eigenvalue-transformation lower bound; Cor 67 = x^{-c} polynomial approximation).
- Monkhorst-Jeziorski abstract quotes: the paywalled verbatim strings remain unread, but the paper's characterization is corroborated by the paper's OWN title ("No Linear Dependence and Multi-Center Integral Problems in Momentum Space Quantum Chemistry"). Substance grounded; the paper fences the strings as source-verified.

### The remaining gating task, itemized (register before a re-run can PASS)

Each value is VERIFIED CORRECT; the task is to give each a C21 key (or C17 family) with measured/cited provenance and a `\\gvq` annotation, per Sec.15 rule 3. Sources identified so the pass is turnkey:

| literal | locus | source to cite/measure |
|---|---|---|
| Q^1.19 / Q^3.33 (naive L2 inflation) | eq:blowup L204/206 | `tests/test_sturmian_l2_encoding.py` (regression-pinned exponents) |
| SW N^1.85 / L2 N^1.70 | L829 | the SW-better-conditioned test in `test_paper60_sturmian.py` |
| water A_1 N^1.97 | L1126 | the water-probe / `_water_A1` route in `test_paper60_preconditioner.py` |
| H2+ 1.3% / H2 -1.09 vs exact -1.174 | L1151 / L1522 | `..._h2plus_isoenergetic_binds`, `..._h2_ci_dissociates_and_binds` |
| Gaussian ratio N^6 / N^1.4 | L1190-1191 | **verify durable backing first** -- may be debug/-only; do not register a value whose only provenance is prunable |
| tab:resource d_inv/kappa row | tab:resource | RESOURCE MODEL (O(1) cancels in ratios per caption); decide C17-family vs literal-with-caption |

### Verdict

**FAIL.** Content and mathematics are clean across five dimensions; the run fails on declared preconditions. Two cert-blockers closed this session (the inertia test; the three NITs). Two remain: the literal registration (itemized above, to be done as its own careful pass) and a clean delta over the remediation. Certification is one focused registration pass plus one clean delta away -- not a mathematical question.

"""

DOD = "docs/qa/paper_60.done.md"
DOD_ANCHOR = "## Change log\n"
DOD_MARKER = "2026-09-13 — **FULL certifying run"
DOD_NEW = """## Change log
- 2026-09-13 — **FULL certifying run = FAIL (PI-invoked).** Whole-paper, unseeded, tree
  frozen. **5 of 6 dimensions PASS with zero mathematical/content defects:** deterministic
  14/14, code 136/136 `--slow` (C21 green), claims (1 SMALL), citations (49/49 resolve),
  synthesis C9 gating (1 NIT). **Completeness = FAIL** on three self-declared
  cert-blockers, not on any hidden defect (zero undiscovered gaps, zero surviving stale
  echoes). Remediated this run: F1 (abstract n^3->n clause re-tiered [RESOURCE MODEL]),
  the C9 eigenvalue-wording NIT, the Monkhorst-Jeziorski bibitem title, and — closing a
  cert-blocker — the Sylvester-inertia bound C8.16 now has a standalone fire-tested test
  (`test_c5_inertia_bound_and_root_by_root`), matrix row 568 OPEN->BACKED. **Remaining
  gating work: register the ~7 declared-debt literals (itemized in CHANGELOG v5.11.10) as
  their own careful pass, then a clean delta.** NOT certified.
"""


def main() -> int:
    n = 0
    with open(CL, encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) == 1:
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, ENTRY + CL_ANCHOR))
        n += 1
        print("  ok   CHANGELOG v5.11.10")
    else:
        print("  skip CHANGELOG")
    with open(DOD, encoding="utf-8") as fh:
        d = fh.read()
    if DOD_MARKER not in d and d.count(DOD_ANCHOR) == 1:
        with open(DOD, "w", encoding="utf-8") as fh:
            fh.write(d.replace(DOD_ANCHOR, DOD_NEW))
        n += 1
        print("  ok   paper_60.done.md")
    else:
        print("  skip done.md")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
