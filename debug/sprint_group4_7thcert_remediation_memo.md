# /qa group4 — 7th cert (whole-group, PI-fired) — FAIL → remediated (v4.62.0, 2026-07-02)

The certifying recert after six cycles. Criteria frozen (verified: since the v4.58.0
freeze the DoD changed only in the dated v4.60.0 watch-note syncs + chronicle entries;
criteria.md untouched). Worktree at HEAD 51a8b1a (v4.61.0); answer key
`debug/qa/group4_seed_key.json` (run-7 set — all classes varied from run 6).

## Calibration scorecard — panel FULLY CALIBRATED

**Sensitivity 8/8** — every seed caught by its own agent; three also caught
cross-dimension (r7s5 by claims+synthesis+code-P20; r7s6 by claims+code-P23+synthesis;
r7s1 by code-P14+code-P20). **Specificity 6/6** — no control flagged MATERIAL
(Trenev App. B independently re-confirmed; the near-parity market lead verified honest;
342.2/592 pins, five-type, the new MPO universal-core form, R_eq 3.227 all protected).

| Seed | Class | Catcher | Note |
|---|---|---|---|
| r7s1 scaling band [2.2,2.8]→[1.8,3.5] | S3 over-wide tolerance on a C8 exponent | code-P14 ✓ | best analysis of the run: the band admits 3.15 — the single-geometry value the composed headline differentiates from; "BLAS variation" bogus for integer counts; matrix contradicts |
| r7s2 μ_free expected-from-classifier | S2 tautology | code-P16 ✓ | counterfactual: the historical 2ν² bug passes the seeded test, fails the row2 cross-check |
| r7s3 library ==37 → >=35 floor | S3 floor guard | code-P20 ✓ (graded NIT — severity variance noted; detection precise) | |
| r7s4 He-4 Coulomb sign check → magnitude-only | S3 direction kill | code-P23 ✓ | reviewer out-analyzed the planter: the direction is transitively pinned by the two non-overlapping absolute 1-norm tests + dE>0 elsewhere — seed partially inert in situ |
| r7s5 GH rate → "bounds the energy directly" | S9 regression-zombie of the v4.58.0 #2 fix | claims-{14,20} ✓ | flagged as self-contradicting the surviving hedge 9 lines below |
| r7s6 P23 abstract "reasonable agreement with experiment… nuclear-structure calculation" | nuclear-honesty flip | claims-{16,23} ✓ | all three contradicting body loci quoted |
| r7s7 Fock1935 Z.Phys. 98→93 | S1 wrong volume on THE foundational cite | citation ✓ | correct locator independently confirmed (ADS/Springer); P23's copy noted correct |
| r7s8 synthesis: κ-parity negatives flipped to working | S8 §3-dead-end zombie | synthesis ✓ | exact §3 quotes (ΔQ=0; three-sprint closure) |

## Verdict: FAIL → all remediated same-session; **3 verified genuine MATERIALs**

- **M1 (P16:292-293, wrong printed identity):** μ_free/N² expansion read "2 − 8/N + 8/N²";
  2(N−2)(N−1)/N² = **2 − 6/N + 4/N²** (the paper's own table column, 1.44 at Ne, matches the
  correct form). PM-verified algebra + real-corpus grep. Fixed 8→6, 8→4.
- **M2 (P14:2781, headline-floor drift):** "(190×–1,712×, Sec. sec:composed_gaussian)" — the
  cited table's own floor is **51×** (51/746/1712) and P20:113 says 51–1,712×. Fixed →51×.
- **M3 (C2 test-validity on the P23 magic-number headline):** `verify_magic_ordering`'s
  min_gap=1e-10 is a necessary-condition check nine OoM below the real gap scale — it could
  never catch a magic-gap collapse. PM probe (reproduced + extended the reviewer's):
  at eq:optimal (0.1709/0.0211) the six lower magics are the SIX LARGEST gaps ≤126, but
  **126 is non-dominant** (0.107ℏω; ten non-magic boundaries larger: 6,14,16,32,40,56,78,
  90,104,118), and boundaries above 126 (136, 154) are n_max=7 truncation-edge artifacts.
  The paper's "recovered as gaps" prose was literally honest; the GUARD was the defect.
  Fixed: new `tests/test_paper23_magic_gaps.py` (3 passed — seven gap/ℏω pins + dominance-
  of-six + 126-non-dominant band, region-scoped ≤126 with the truncation-edge rationale) +
  an honest gap-magnitude disclosure sentence in P23 (126 non-dominant; edge caveat).
  Production scan code untouched (changing its default would move eq:optimal — result-changing).

## Completeness-critic → the balanced-table λ sweep (the run's biggest cleanup)

The critic asked for the never-verified `tab:molecules` balanced cells. Live sweep of all
12: **every Pauli count exact** (+identity convention), but the **λ_ni column was stale at
drafting vintage in 9 of 12 cells** (0.3–7%; AsH₃ worst 886.2→824.2; NaH/HCl/… exact —
per-cell vintage mix, the v4.56.0 stale-table class, NOT a code regression). All cells
re-synced to the live builder; **new `tests/test_paper20_balanced_lambda.py`** pins every
cell (the balanced-table analog of test_paper14_rel_lambda); caption carries the sync note.
Also caught: the run-6 KH λ fix (27.5→28.15) had been memo-listed but never applied — the
PM's own propagation miss. Landing it exposed a second layer: the run-6 probe value 28.15
is UNREPRODUCIBLE at any candidate geometry (live: 31.57 at the factory experimental
R=4.243; 29.56 at R=3.566) — an unreproduced reviewer number I had imported without
recomputing. Final authoritative cell: **31.6** (factory-default convention, consistent
with all other rows), pinned. Lesson re-learned: reproduce before syncing, even from a
calibrated reviewer.

## NITs fixed on sight

GH "truncation-error bounds" phrasing → "truncation-convergence rates" (P20:56/385/1071 +
synthesis ×2; the claims-vs-code severity disagreement on this cluster recorded — adjudicated
NIT per the claims reviewer + six prior calibrated panels, but ALL loci tightened);
Navratil2000 title↔venue mixed pair → PRL title restored; synthesis 13× QWC given the
raw/pair-diagonal qualifier; P14 symmetry-adapted sentence reworded to what the test proves
(ext vs ext+hidden; naive comparison carried by the H₂ ℓ-parity test); P23 deuteron
export-round-trip docstring made honest (FCI match informational, sign-limitation \ref);
dirac-metric divergence test given a local magnitude assert; matrix — 878-pin upgrade
(BACKED-WEAK→SOUND; stale "not pinned" note), Z=1–36→56, O(Q^2.5) row "never fit" note
stale→LiH-fit BACKED-WEAK, magic row →BACKED-SOUND via the new gap test, two §VII rows
added (cross-register 69 tests + magnetization 57 tests were backed-but-unlogged).

## Honest ceiling

P20's heavy full-library slow suites subset-run (11.10/9.23 confirmed via exact
representative pins); z2_tapering + breit suites ran green but not fully enumerated
(back the §V Breit aside, cleared by run-6's focused pass); R²=0.997 remains a logged
secondary gap (fit_pauli_scaling returns prefactor, not R²); inline debug/ cites remain
C14-advisory standing debt. Run interrupted once by ECONNRESET + once by the account
spend limit; both recovered via transcript-resume, no loss.

## State for the 8th (certifying) run

All 3 MATERIALs + the λ column + all NITs remediated; 6 deterministic gates PASS
post-remediation; 5 papers compile 0-errors; worktree removed, zero seed leakage
(8/8 signatures grepped absent). Per-run detail: CHANGELOG v4.62.0.
