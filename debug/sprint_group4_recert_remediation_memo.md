# Sprint: /qa group4 re-cert (FAIL) + remediation of all 4 findings — canonical memo

**Date:** 2026-06-28 · **Version:** v4.55.0 · **PI direction:** `/qa group4` (re-cert,
2nd run); then **"fix all 4 now, but hold on recert."**
Predecessors: v4.52.0 (pre-work), v4.53.0 (A/B framing rule + disclosure), v4.54.0 (first
cert FAIL→remediated, 9 findings). Re-cert run notes: this memo. Answer key (gitignored):
`debug/qa/group4_seed_key.json`.

## 1. The re-cert run (held for record; the fixes below are the action)

Second cert of group4 (Papers 14/16/20/23 + synthesis), whole-group, against the frozen DoD
(`docs/qa/group4.done.md`) at commit **7834faf** (v4.54.0). Detached seeded worktree, **8
fresh/varied calibration seeds** (4 code S2/S3, 1 framing-zombie, 1 status-overstatement, 1
wrong-arXiv-ID, 1 synthesis-zombie), fresh path-pinned panel (4 code 1/paper, 2 claims chunks,
1 citation, 1 synthesis; opus).

**Verdict: FAIL — fully calibrated (trustworthy).** Sensitivity **8/8**, specificity **5/5**
(the v4.54.0-disclosed `sec:eri_rule` + the now-correct Swain cite both correctly NOT flagged).
Worktree removed, leak-scan clean (all 8 seeds = 0 in the real corpus).

**Why the re-cert mattered:** it caught that the v4.54.0 citation remediation (finding A below)
was itself partly wrong — the value of the fresh-adversary gate (qa.md principle 2).

## 2. The 4 genuine MATERIAL findings (seeds excluded) — all remediated

| # | Finding (dimension) | Fix |
|---|---|---|
| A | **`Swain2022` misattribution** — the v4.54.0 fix (`Sunaga2025`→`Swain2022` arXiv:2211.06907) was wrong. The RaH-18q benchmark is **Chawla et al. arXiv:2406.04992 = PRA 111, 022817 (2025)**; **47,099 = two-electron integrals, not Pauli; the relativistic Pauli count is 12,556** (non-rel 2,740). — citation/C4 | Re-keyed `Swain2022`→**`Chawla2024`** (P14 bibitem + P20 `@article` + all `\cite`); RaH-18q 47,099→**12,556** (rel) / 2,740 (non-rel); recomputed the P20 "Chawla ratio" column on denominator **12,556** (×3.751) + P14 native-ratio column; reworded "Swain"→"Chawla et al." throughout |
| B | **trenev2025 range still attributed** — the v4.54.0 fix kept "Trenev report the Q^3.9–4.3 range"; deeper primary check: Trenev is a **vibrational-spectra** study with **no electronic Gaussian scaling exponents at all** — citation/C4 | Trenev = **methodology cite only**; the Q^3.92/4.25 exponents AND the Q^3.9–4.3 range they bracket are GeoVac's own OpenFermion recompute. Fixed P14 (7 loci + bibitem note), P20 caption + bib note, synthesis |
| C | **P14 §origin framing zombie** (§origin, ~L2183) — "structural sparsity advantage that grows with angular complexity" presented the pair-diagonal d-block density (4.0% vs 8.9%) as a genuine saving — claims (C16-class) | Disclosed: it is a **pair-diagonal-approximation** property; under the exact global-$M_L$ rule the $d$-block is **denser**, so the "advantage that grows with angular complexity" is an artifact. (Also dropped the wrong "same three multipole orders as $l=1$" — $l=1$ has two.) |
| D | **P20 STO-3G market-test convention mix** — the LiH STO-3G "907 @ Q=12" is **raw** JW, but every other Gaussian count in the table (LiH 6-31G/cc-pVDZ, all H₂O rows; Q=20,36,12,46 are the reduced sizes) is **2-qubit-reduced** (matching `GAUSSIAN_LIH_PUBLISHED`, which lists STO-3G LiH at 276 @ Q=10) — claims/code | **Disclosure** (not a number reversal): caption + market-test prose now state the LiH multiplier is raw-vs-raw (334@Q30 vs 907@Q12) and **narrows toward parity under a uniform 2-qubit-reduced convention** (276 vs 334); raw-vs-raw and reduced-vs-reduced bracket the genuine advantage — parallel to the CF-1 ERI-rule disclosure |

**Worth-doing NITs swept:** P14 "30 systems (Z=1–36)" v2.4.0 note → **37 systems** (matches the 7 other library loci); P20 §heading "Propinquity-derived basis-truncation error bound" → **"GH-convergence-derived…"** + legacy-field-name note (Paper 38 is state-space GH, not propinquity); P20 L382 "propinquity error bounds" → "GH-convergence error bounds".

**Deferred (with rationale, non-blocking):** the `propinquity_bound` **code field** rename (breaking public-API change; legacy name now disclosed in prose); the P14 `debug/sprint_df_multipole_lift_memo.md` + other `debug/` cites (the policy-deferred corpus-wide C14-advisory sweep, CLAUDE.md §9 standing debt); the 11.10 vs 11.11 ±identity convention (adequately disclosed at P14:63 "exact, non-identity terms"); two `docs/*_proposal.tex` pre-production drafts still carrying the old Sunaga attribution (out of the authoritative corpus — flag for a future docs sweep).

## 3. Verification

- Deterministic gates C11/C13/C14/C15/C16 `--gate group4`: **PASS** (C16: the 3 pair-diagonal
  occurrences correctly carry withdrawal flags; no live retracted claim).
- No active `Swain2022`/`Sunaga2025` cites remain (only the intentional re-attribution
  history-comment block); `Chawla2024` cites resolve (P14 bibitem + P20 `@article`).
- All 5 papers compile three-pass clean: **P14 26pp, P16 7pp, P20 12pp** (+1pp from the D
  disclosure), **P23 11pp, synthesis 4pp** — zero undefined refs/cites; P20 bibtex 0 warn/0 err.
- **No production code edited this turn** (paper/cite/framing edits only) → no `/regression`
  needed; the v4.54.0 test suite (eri_rule 3, scaling, classifier 197, balanced_row2 39) stands.

## 4. Honest scope

- **Theorem grade:** none.
- **Process verdict (trustworthy):** the re-cert was fully calibrated (8/8, 5/5), so the FAIL
  is trustworthy and every fix traces to primary text, the paper's own data, or the
  web-verified Chawla publication.
- **Biggest catch:** finding A — the re-cert caught my own v4.54.0 citation fix as a
  misattribution. This is exactly why the fresh-adversary re-cert exists.
- **D is a disclosure, not a headline reversal:** the LiH 2.7× market test is internally
  consistent raw-vs-raw; the fix discloses its convention-sensitivity (it narrows to parity
  under uniform qubit-reduction) rather than asserting a reversed number on uncertain reasoning.
  **Flagged for PI awareness** as a market-test-robustness disclosure (parallel to CF-1).
- **Re-cert HELD per PI direction** ("hold on recert"): the certified-PASS confirmation run is
  the next `/qa group4` invocation, on the PI's timing — NOT auto-fired.
