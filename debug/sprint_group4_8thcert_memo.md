# /qa group4 — 8th cert (the FULL certifying run) — 2026-07-02

**Shape:** FULL (the run-shapes protocol's final certifying pass; precondition met — delta-1
CLEAN, v4.62.2). Worktree seeded at 906d2dc. First full run under the v4.62.1 cost package:
Sonnet-tiered code+citation reviewers (two seeds each), terser reports, C17 live.

## Calibration — the panel

**Sensitivity 13/13** (every seed caught by its own agent). **Specificity 6/6** (no control
flagged; the citation reviewer handled c1/Trenev exactly right — UNVERIFIABLE-this-session
deferring to the documented triple-verification, not a false flag).

| Agent (tier) | Seeds | Caught | Notes |
|---|---|---|---|
| code-P14 (sonnet) | r8s1a MPO ==→<=, r8s1b JW 1e-4→1e-1 | 2/2 | graded s1a SMALL-mitigated via the JSON-pin chain (fair); one confused git-archaeology note (attributed the seed to a past commit — the catch stands) |
| code-P16 (sonnet) | r8s2a row2 tautology, r8s2b ==36720→>0 | 2/2 | s2a with the 2ν²-regression counterfactual + Mg-row exposure analysis |
| code-P20 (sonnet) | r8s3a count-assert deleted, r8s3b R_eq band | 2/2 | verified the seed's own justifying comment was FALSE; rebuilt all 12 systems |
| code-P23 (sonnet) | r8s4a 0.005→0.05, r8s4b floor>400 | 2/2 | read the production source to prove "no grid/interpolation exists"; also caught r8s6 + r8s7b cross-dimension |
| claims-{14,20} (opus) | r8s5 "modestly ahead" inversion | 1/1 | internal-consistency forcing worked as designed |
| claims-{16,23} (opus) | r8s6 Fock-rigidity flip | 1/1 | six contradicting loci + theorem + §3 marshaled |
| citation (sonnet) | r8s7a Tung 134→138, r8s7b Bargmann 14→17 | 2/2 | DOI-level confirmation both |
| synthesis (opus) | r8s8 1.69→1.49 | 1/1 | cross-referenced body + P14; noted 1.47/1.57 are different quantities |

**The Sonnet tier held**: all 10 tiered-agent seeds caught, several with deeper forensics than
prior Opus panels (the P23 reviewer's no-interpolation-exists proof; the P20 reviewer's
false-comment check). The tiering-with-calibration design is validated at full-run scale.

**Protocol hardening (2 items, adopted):** (1) two reviewers legitimately localized seeds via
`git diff` in the worktree (uncommitted seeds are enumerable) → future runs commit seeds onto
the worktree branch so its git state is clean; (2) the delta-1 internal-consistency lesson was
already baked into qa.md and demonstrably carried (r8s5/r8s6 catches).

## Genuine findings (verified against primary text/code) → ALL remediated in-run

**MATERIAL:**
1. **M2 (C7/C8, P14)** — §d_orbitals + Note-v2.4.0 stale at the two-hydride/no-automation
   v2.4.0 state vs P14's own abstract + all of P20 (all-ten). Registry-verified truth: 37
   systems incl. sch…znh; classify_atom(23–30) type-F supported; 35 = 37 − He − H₂ exact.
   → both loci updated to the current state (benchmark honesty retained: resource-level only).
2. **D2 (C1, P14)** — §hopf_tapering "6-system machine-precision spectrum preservation
   verified directly" + "254 total savings" had NO test (only H₂; numbers lived in a transient
   debug memo) — the identical hollow-citation class fixed one subsection later in v4.60.0.
   → new `TestLibraryTaperingHeadlines` (test_z2_tapering.py): He sector-sweep spectrum test
   (2nd continuously-verified system) + per_block delta_Q library sum == 254 pin (75s slow);
   paper sentence now states exactly what is continuously verified (H₂+He) vs v2.6.0-measured.
3. **D3 (C1/C8, P20)** — FCIDUMP round-trip "seven-system sample across LiH, BeH₂, H₂O, NaH,
   KH, MgH₂, CaH₂" vs the actual test list (CO, ScH, H₂ in; KH/MgH₂/CaH₂ never tested).
   → sentence corrected to the tested set + inline test cite.
4. **M3/N1 (C8/§1.5, P14+P20)** — "at matched qubit count" mislabeling the Q30-vs-Q12 (LiH)
   and Q70-vs-46 (H₂O) comparisons → "matched raw-JW convention" with the qubit ratios stated.

**§BeH₂–H₂O re-verification (the P20-code coverage-gap finding, PM-probed live):**
λ convention = identity-INCLUDED. Live: balanced BeH₂ 306.4 (printed 304.7), H₂O 1,511.1
(printed 1,509), composed-H₂O 28,055 (printed 28,053) — vintage rounding drift; but composed
BeH₂ **354.9 was the deprecated legacy-builder value → 373.4 live** (the ~5% gap the legacy
path's own docstring documents). The QPE-direction claim SURVIVES (18% lower, was 14%);
0.86×→0.82× (P14). tab:resources BeH₂ electronic 67.4→66.0 (P14 already carried 66.0; live
confirms). All six cells now pinned (`test_paper20_library.py::test_beh2_h2o_qpe_regime_one_norms_pinned`
+ the 32.6 LiH pin). New C17 family `beh2-h2o-qpe-onenorm-vintage` registers the retired values.

**Other genuine (SMALL/NIT), all fixed:**
- P23:567 body ~3×10⁵ → ~5×10⁵ (second locus of my own run-7 table fix).
- P23 multipole termination "verified at ℓ∈{0,1,2}" but tests instantiated ℓ≤1 → live-verified
  tight at L=4 (V₄≡V₁₀ machine-exact; L=4 contributes 8.5e-5 on d×d) + new test point.
- P16 heavy rows (Kr/Xe/Rn/Og, fallback path, zero coverage) → `TestHeavyRowFormulaPath` pins.
- P16 five orphaned bibitems removed (all real works, never cited).
- Citations: Tung author Wei-Chia→Wei-Cheng (+ the seed-confirmed vol 134), Burkat–Fitzpatrick
  initials → J./N., Chawla title synced to arXiv v2.
- P14 §111 channel ERI counts (12/24/16 summed 52 vs the sentence's own 40) → rewritten to the
  counts forced by the paper's own 7²+4² combinatorics (24 k=0 non-direct + 16 k=1 = 40);
  Pauli split 16/24/16 (tested, =56) unchanged.
- P20 "D_e×2.4"→"~75% vs experiment" (baseline named), 74%→75%, "guarantees"→"implies";
  balanced-MPO panel got a live-recompute drift guard (LiH, bit-exact vs shipped JSON).
- DoD:108 provenance note 3.015→3.227-computed (secondary-provenance carve-out).

**Logged, not remediated (unchanged tier):** BeH₂/H₂O FCI accuracy cells (vintage,
accuracy-tier), VQE openfermion-vs-qiskit paragraph (demo-class), D_e 0.161 pin, rel-λ
n_max=1/3 rows, 0.20%@n=3 (PI-deferred), P16 inline-tier-labels editorial note (7 calibrated
runs accepted the hedged prose), "propinquity_bound" metadata key naming (code-key NIT class).

## Deterministic + compiles
7/7 gates PASS post-remediation (incl. C17 with the new family); 4 papers recompile 0 errors
(synthesis untouched this round). All new/affected tests green: 30-passed slow batch + He
spectrum + 254-sum + heavy rows + ℓ=2 + magic gaps + balanced-λ.

## Completeness-critic + focused gap-closure (loop cycle 2)
Critic named ~12 regions un-enumerated BY THIS RUN's panel (most verified in prior calibrated
runs — e.g. the NaH snippet was a run-5 fix). Focused gap-closure reviewer on all twelve:
**11 SOUND/NIT, 1 MATERIAL** — **eq:dirac_fs printed (Zα)⁴ where its own backing code
(`fine_structure.py:259`), its own test docstring, and its own spin-orbit term (L1796, α²)
all carry Z⁴α²** (the Hartree-unit form; (Zα)⁴ needs the mc² factor the equation lacks).
The same wrong-printed-equation class as run-7's P16 M1 — surfaced only because the
final-run loop-until-dry sent a reviewer into a region seven calibrated panels never
enumerated. Fixed → Z⁴α² + "(Hartree units)"; fine-structure suite 43/43. Sound findings
en route: the Hund/Cr arithmetic verified (141=105+36 vs 136=91+45), the P23 sign-bug
disclosure's count/1-norm-invariance reasoning confirmed, the NaH snippet reconciled
(composed-vs-balanced, full-vs-nonid) + annotated, "T3's spin-orbit"→Track T2, TC vintage
2–3-term mismatch disclosed. New logged gaps: H-set 6j percentages (no test), DF-rank-7
pin migration (debug-memo-backed), TC table cells (vintage).

## In-run re-scan of the remediation diff (loop cycle 3 — dry)
One reviewer over all 18 remediation hunks, live-corpus verification, one seed (the
eq:dirac_fs hunk pasted as α⁴): **caught** (verified the live file directly, α², and
flagged the paste as wrong-if-live). 17/18 hunks SOUND against code+tests; **1 genuine
tail defect introduced by my own h1 sync** — the sibling sentence P14:1036–37 ("BeH₂ 355,
H₂O 361") left contradicting the corrected line above and carrying the legacy with-PK
value in an electronic-only list. Fixed → LiH 33 / BeH₂ 66 / H₂O 359 (all live-probed +
pinned); P14 recompiles 0 errors; C17/C13 PASS. The remaining fix was a two-number sync
PM-verified against the already-pinned live values — loop declared dry.

## Verdict — **PASS (group4 CERTIFIED)**
Every gating dimension exercised + calibrated (13/13 sensitivity incl. all 10 Sonnet-tier
seeds, 6/6 specificity) + clean at run close (zero remaining verified MATERIAL after the
in-run loop-until-dry: panel → remediation → critic → gap-closure → remediation → re-scan
→ tail fix). Deterministic ×7 PASS (C17 + the new vintage family); 4 papers 0 errors;
all new pins green (He spectrum, 254-sum, QPE 1-norms ×4+2, heavy rows ×4, ℓ=2, 32.6,
MPO balanced drift guard). Worktree removed; zero seed leakage (12 signatures).

**Honest ceiling.** PASS = survived calibrated detectors for the seeded defect classes +
the pre-registered criteria across three in-run convergence cycles — not provably perfect.
Named unexercised/deferred items: 0.20%@n=3 (PI-deferred, heavy), BeH₂/H₂O FCI accuracy
cells + VQE-comparison paragraph + D_e 0.161 (vintage/demo-class, matrix-logged), rel-λ
n_max=1/3 rows, H-set 6j %s, DF-rank pin migration, the P16 inline-tier-label editorial
note, the corpus-wide debug-ref sweep (standing §9 debt). Protocol hardening adopted for
future runs: commit seeds onto the worktree branch (git-diff enumerability); the run-4/
delta-1 internal-consistency forcing demonstrably carried.

**Cost note (the v4.62.1 package's first full outing):** panel 8 agents ≈ 1.90M subagent
tokens (5 Sonnet-tiered ≈ 1.20M of it at Sonnet rates) + critic 259k + gap-closure 193k +
re-scan 146k ≈ **2.49M subagent total** — comparable token count to prior full runs but at
a substantially lower dollar cost (≈half the panel volume billed at the Sonnet tier), with
MORE seeds (13 vs 8) and three convergence cycles instead of one.
