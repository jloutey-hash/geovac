# Sprint: /qa group4 confirmation cert (3rd run, FAIL) + remediation — canonical memo

**Date:** 2026-06-29 · **Version:** v4.56.0 · **PI direction:** released the v4.55.0 recert hold
(`/qa group4`); then **"go ahead"** (remediate); then **"diagnose first, then release"** (the rel-λ conflict).
Predecessors: v4.52 (pre-work), v4.53 (A/B framing rule), v4.54 (1st cert FAIL→rem, 9 findings),
v4.55 (re-cert FAIL→rem, 4 findings). Answer key (gitignored): `debug/qa/group4_seed_key.json`.

## 1. The confirmation cert (3rd cert run @ 4092dcd/v4.55.0)

Whole-group, frozen DoD. 8 fresh/varied seeds, fresh path-pinned panel (4 code 1/paper, 2 claims
chunks, 1 citation, 1 synthesis; opus). **Verdict: FAIL — fully calibrated.** Sensitivity **7/7**
valid seeds (c1 excluded as an inert/dud seed — removing a redundant reversal-assert didn't degrade
the test); specificity **6/6** (incl. the citation-reviewer *confirming* the now-correct Chawla2024
cite, and the eri_rule + STO-3G disclosures accepted). Worktree removed, leak-scan clean.

The run peeled a **deeper layer** than the first two certs (the claims reviewer ran a full
relativistic-table internal-consistency cross-check that prior runs' chunking never did).

## 2. The findings (all remediated)

| # | Finding | Fix |
|---|---|---|
| **M1** | **P14 relativistic-section staleness.** `tab:sunaga` GeoVac native = stale **805/805/534** vs the paper's own `tab:spinor_resource` **1413/1413/942**; "16–24× / 40–67×" headline rode the stale counts; obs-3 QWC 6571 (8.3×) stale. (v4.55.0's A-fix had recomputed P14's ratio col on the stale 805 while P20 used 1413 → P14↔P20 inconsistent.) | counts→1413/942, ratios→0.113/0.075, advantage→**9–13×**, projected→**17–32×** (honest O(Q^2.5)), obs-3→**11865 (15.0×)**; all code-validated |
| **M2** | **P20 abstract "binds LiH at R_eq=3.015"** vs the body's own computed **R_eq=3.227 (7% high)**; 0.20% is an n_max=3/84q single-point at the experimental geometry, mis-paired with the n_max=2/878-Pauli/30q resource | abstract + tab:scope reworded to separate the n=2 resource from the n=3 accuracy; R_eq→3.227 |
| **M3** | **0.20% n_max=3 accuracy headline has no test** (deferred in test_paper19) | logged as a coverage gap (84q FCI single-point test is heavy); M2 fixed the abstract conflation |
| **M4** | **synthesis "~1/M²"** vs P14's "~1/M" | → ~1/M (M^{-0.85}) |
| **M5** | **P23 §4 title "First Nuclear Qubit Hamiltonian"** (unqualified priority) | → "First **Two-Species** Nuclear Qubit Hamiltonian" (matches the abstract) |
| **+ proj.** | P14↔P20 projected mismatch (P20 said 21×) | both → honest O(Q^2.5) **17–32×** |

## 3. The rel λ_ni conflict — "diagnose first" (the 6th, code-vs-paper finding)

The M1 agent surfaced that the code's **relativistic 1-norm** for *first-row* LiH/BeH matched
**neither** paper's table (BeH n=2: code **143.96** vs table **40.26** = 3.6×), while frozen-core
CaH/SrH/BaH matched **exactly**. PI: "diagnose first." Diagnosis (conclusive — **stale table, NOT a
code regression**):
- **Physical impossibility in the table:** BeH scalar λ = **139.12**, so the table's BeH *rel* λ
  (40.26) claimed the relativistic 1-norm — with 4.24× more Pauli terms — was **3.5× BELOW its own
  scalar**. Impossible. The code's 143.96 (≈ scalar +3.5%) is sensible.
- **Pattern:** code rel λ is a uniform modest increase over scalar (+3.5–12.5% at n=2, +20–30% at
  n=3) across all three molecules — CaH's code value (which equals the table) follows the same
  pattern; the table's BeH (rel < scalar) was the lone outlier.
- **Git:** the table λ dates to v2.15.0; the first-row composed-λ path changed afterward (chemistry
  re-entry v3.54/v3.56). The Pauli *counts* (1413/942) stayed (test-pinned, match); only the λ
  *coefficients* drifted, and the table's first-row λ column was never refreshed. (`include_breit`
  defaults False, so Breit is not the cause — the reported λ is the non-Breit build the papers describe.)
- **Resolution:** code is authoritative. Updated `tab:spinor_resource` (P14) + `tab:paper20_tier2`
  (P20) + obs-2 + the P20 prose bullet to the code values; reworded obs-2 from the (doubly-wrong)
  "decreases 30–40%, favorable" to the honest "**modest ~20–30% increase, no order-of-magnitude
  inflation**"; added the missing pinning test **`tests/test_paper14_rel_lambda.py`** (the coverage
  gap that let it drift — λ was never test-guarded).

## 4. Verification

- New test `test_paper14_rel_lambda.py` (pins λ_ni^rel n=2 + a rel≥0.9·scalar sanity that catches
  the stale-table mode): **4 passed**.
- Deterministic gates C11/C13/C14/C15/C16 `--gate group4`: **PASS**.
- All 5 papers compile clean (P14 26 / P16 7 / P20 12 / P23 11 / synth 4 pp; zero undefined; P20 bibtex clean).
- P14↔P20 now consistent (rel λ 40.59/143.96/218.78/413.90; native 0.113/0.075; projected 17–32×; advantage 9–13×). No residual 805/534/6571/8.3×/21×.

## 5. Deferred (flagged NITs, non-blocking)

Citation NITs I will NOT auto-fix (my citation fixes were wrong twice this session — Sunaga→Swain):
**Pachucki "2023" vs bibitem-2018** (needs a source-check on whether the 2018 paper supports the FW
two-particle claim or a genuine 2023 paper was meant); **rocca2024** 4th author ("P.D. Johnson"→
"P.J. Ollitrault" per the reviewer, but P.D. Johnson is a real QC author — verify); **caesura**
"A. Pol"→"W. Pol". And **ScH 277/278** (convention-ambiguous — 278 with-identity is consistent with
the 334/778 table convention). All flagged for a careful source-check pass.

## 6. Honest scope

- **Theorem grade:** none. **Process verdict (trustworthy):** the cert was fully calibrated (7/7, 6/6).
- **Biggest catch:** the rel λ_ni stale-table drift — a real, physically-impossible table value
  (BeH rel < scalar) that survived because λ had no pinning test. Now fixed + guarded. This is the
  6th genuine finding the cert arc surfaced beyond what pre-work + 2 prior certs caught (peel-one-layer).
- **Self-implication:** M1 also caught that v4.55.0's own A-remediation left P14↔P20 inconsistent
  (ratio col on stale 805). The fresh-adversary gate working as designed, a third time.
- **Re-cert HELD** per the gate (PI timing) — the certified-PASS confirmation is the next `/qa group4`.
