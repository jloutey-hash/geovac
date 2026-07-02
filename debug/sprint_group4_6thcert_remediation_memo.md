# /qa group4 — 6th cert (whole-group) — FAIL → remediated (v4.60.0, 2026-07-01)

PI-fired confirmation run after five FAIL→remediated cycles (v4.54.0–v4.58.0/v4.59.0).
Criteria frozen at v4.58.0 (verified: no DoD change since; the v4.59.0 M-A call touched
papers only). Worktree seeded at HEAD 1c049ce; answer key `debug/qa/group4_seed_key.json`.

## Calibration scorecard — panel FULLY CALIBRATED

**Sensitivity 8/8** (every planted defect caught), **specificity 6/6** (no known-good
control flagged MATERIAL; several explicitly two-way verified sound: the eri_rule
disclosure, Trenev Table 5, Chawla RaH numbers, the 342.2 MeV pin, five-type A/B/C/D/E,
R_eq 3.227).

| Seed | Class | Catcher | Note |
|---|---|---|---|
| s1 rel-λ pin → regenerate-from-builder tautology | S2 | code-P14 ✓ (+ code-P20 independently) | "the pin can never fail; hardcoded 40.59/143.96/18.68 dead" |
| s2 δ_1s divergence asserts → ordering-only | S3/weaker | code-P16 **reviewer #2** ✓ (reviewer #1 MISSED — de-calibrated, discarded; fresh strength-matching re-dispatch caught it with a capped-δ counterfactual) | the run's one calibration wrinkle |
| s3 library-37 pin → registry self-tautology | S2 | code-P20 ✓ (live demo: 3 fake systems still pass) | |
| s4 deuteron 1-norm tol 1.0→120 ("BLAS headroom") | S3 | code-P23 ✓ (\|227.3−342.2\|=114.9<120 — re-admits the historical drift) | |
| s5 O(Q^2.5) "derived exactly … independent of any fit" | S4 | claims-{14,20} ✓ | |
| s6 magic numbers "predicted … no empirical input" | nuclear-honesty | claims-{16,23} ✓ | |
| s7 Thompson1977 → NPA 302,159 (1978) wrong locator | S1 | citation ✓ (correct A286,53 (1977) identified) | |
| s8 abstract "2.7× at matched accuracy" zombie | S9 | synthesis ✓ | the abstract-clause class |

**Verdict: FAIL** — calibrated panel + 3 verified genuine MATERIALs (below). All
remediated same-session. All dimensions exercised (deterministic ×7 PASS pre- and
post-remediation; 4 code + 2 claims + 1 citation + 1 synthesis + completeness-critic +
2 focused gap-closure re-dispatches on the critic's regions).

## Verified genuine MATERIAL findings (all remediated)

- **M1 (C8, second-locus propagation):** the v4.58.0 M-C "Z=1–36 → Z=1–56" fix never
  reached P20's conclusion (:1048) or four P14 loci (60, 1381, 1523, 1543). Registry
  probe re-confirmed SrH(38)/BaH(56) among the 37. All five loci → Z=1–56 (H–Ba).
  Independently converged on by claims-{14,20} + synthesis + code-P20.
- **M2 (C2/C3):** P14 §ℓ-parity claimed "Verified bit-exact spectrum preservation …
  (test tests/test_extended_tapering.py)" but that file's `_spectrum_lowest` was defined
  and NEVER CALLED — no eigenvalue comparison existed. Fixed by making the claim true:
  new `test_extended_hopf_ell_spectrum_preserved_h2` (naive JW vs min-over-sectors
  extended-tapered H₂; PASSES, agreement <1e-10) + paper sentence now cites the actual
  assertion.
- **M3 (C8, v4.52.0 propagation, PM-caught):** P14 tab:multi_center still carried the
  three organics (CH₂O/C₂H₂/C₂H₆) dropped from the library in v4.52.0 — table said 8
  rows while the prose below said "these five"; registry has no organics. Rows + intro
  ("eight … up to eight nuclei") + two prose count loci (1517, 2815 ethane example)
  removed/replaced with shipping systems. New C16 registry entry `organics-in-library`
  backstops recurrence.

## NITs fixed on sight (carve-out classes)

Papers: P14 residual "~1/M_b²"→"~1/M_b" (:1374); LiH λ 33.3→32.6 + 0.97×→0.95× (P14
:2783/:2789; P20 worked example :161, tab:resources :445, :469 — live
`hamiltonian('LiH',R=3.015).one_norm`=32.59, diagnosed stale drafting-era value, not a
convention difference); P23 "O(Q^4) single-center"→sparsity-reduction-below-naive-O(Q^4)
+ measured 3.15 (:1137); P23 Pauli-round-trip sentence reworded to matrix-path
source-of-truth + sign-limitation \ref (:409, new `subsec:known_limitation` label); P23
ne-scales ~3e5→~5e5; P16 journal-target header comment deleted (off-policy); P20
\date{\today}→fixed date; P14 eri_rule now states the m-swap per-vertex factor
|c²(1,±1;1,∓1)|=√6/5≈0.49 (sympy-verified), sourcing the synthesis's illustration.

Citations (PM-verified per the literature discipline): FriarPayne2005 title↔venue
mismatch → **PRA 56, 5173 (1997)** (web-confirmed; key→FriarPayne1997); Goings-et-al
P450 THC numbers were cited to `lee2021` → new `goings2022` bibitem (PNAS 119,
e2203533119 (2022), arXiv:2202.01244) + repoint; caesura "A. Pol"→"W. Pol"; duplicate
`caesura2025_tensor_factorization` bibitem deleted + cite repointed; burkat title →
"Efficient Circuits with Applications"; Chawla title += "trapped ion";
Pachucki2023→Pachucki2018 key (2 cites); §V.G Migdalek–Bylicki "PRA 57, 3456 (1998)"
WRONG-ID (DOI 404; real joint papers PRA 22/24, 1980–81) → reworded to
"Migdalek–Bylicki-type model potential, e.g. PRA 24, 649 (1981), evaluated here with…";
"Drake 1971 §IV" unverifiable locator → generic Drake-reviews pointer; Brink–Satchler
App.~5 suspect sub-locator dropped (Varshalovich §5.17 primary).

Code/tests: test_paper16_dirac_metric now pins ALL six tab:metric δ rows (was 2);
test_nuclear_electronic 10¹³-ratio guard >1e9→>1e12; test_ibm_quantum_demo shot-noise
bar 10→30 mHa (smoke test, backs no claim, 24.5 mHa observed once — the OPPOSITE
disposition from seed s4, which guarded a §C8 headline); atomic_classifier module
docstring += Type F; composed_qubit sweep docstring/print mark Q^4.60 + 631-anchor as
historical/superseded. Matrix: 2 stale class-name pointers fixed; P16 Z=137 row
NO-TEST→BACKED-SOUND (gap #10 CLOSED); He-4 1.20×/12.25× BACKED-WEAK→**BACKED-SOUND
upgrade** (endpoints exactly pinned — the two-way verdict). DoD watch-notes synced
(~1/M² → ~1/M; "38 molecules Z=1–36" → 37 systems Z=1–56).

## Coverage closure

Completeness-critic found P14 §V.G (1957–2151, the Breit–Pauli/He 2³P benchmark) as a
panel blind spot + smaller regions (fig captions, MPO-bond-rank prose, P20 tier2
table-notes, P16 §I-B). Two focused re-dispatches: **claims — CLEAN across all 8
regions, zero MATERIAL** (the §V.G symbolic claims are genuinely backed by
test_fine_structure.py + test_breit_integrals.py); **citations — the Migdalek–Bylicki
WRONG-ID above** + unverifiable sub-locators (fixed).

## Post-remediation state

All 6 deterministic gates PASS; 4 papers compile 0-errors; 122+21+1(slow spectrum)
affected tests green; worktree removed, branch deleted, zero seed leakage (grep +
git-status verified).

## Held for PI (not auto-decided)

1. **CLAUDE.md §1.5** (PM-locked): "LiH electronic-only 1-norm matches STO-3G at 0.97x
   with 2.7x fewer Pauli terms" — needs 0.97→0.95 AND the un-demoted "2.7x fewer"
   phrasing predates the v4.59.0 M-A call. §2 already synced (0.95×).
2. **MPO bond-rank (Thm 3.2.A) empirical companions** backed only by transient debug/
   drivers — candidate for the tests/-migration pattern (rank2_rate_support /
   wh7_support precedent). Named follow-on.
3. **"Intrinsic tier" label** on the Drake combining coefficients (P14 §V.G ~:2022) —
   focused reviewer judged within-tier but flagged for verification against Paper 18's
   intrinsic-tier definition + Drake 1971.
4. **Recert timing** — the 7th (certifying) run; also note the account spend-limit
   interruption mid-run (agents recovered via transcript-resume; two focused
   re-dispatches burned twice).

Per-run chronicle: CHANGELOG v4.60.0. Seed key: debug/qa/group4_seed_key.json.
