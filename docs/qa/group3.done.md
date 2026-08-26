# Group 3 (Foundations) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group3-specific scope + deltas + bite tracking.

**Scope (non-trunk group3):** Papers **18, 22, 24, 31, 54, 55, 56, 57** + the
**group3 foundations synthesis**. Trunk papers **0, 1, 7** are taken as
already-certified (`/qa trunk` PASS) and not re-litigated except where a group3
claim restates them (C7).

**Deterministic `--gate`:** `group3`.

## Bite tracking (this branch certifies in bites)

- **Bite 1 — PASS** (`/qa group3` run #7, v4.18.0): Papers **22, 24, 31** + the
  synthesis's coverage *of those three*. Calibrated panel, 6/6 sensitivity,
  zero false positives.
- **Bite 2 — run 1 (2026-06-16) = FAIL.** Papers **18, 54, 55, 56, 57** + the
  **full** group3 synthesis (C9 across the whole paper set). Verified material
  defects across all five papers + synthesis (tautological keystone backing,
  conjecture-env K-labels, Paper 38 propinquity drift, count contradictions, no
  matrix rows). Synthesis dimension uncalibrated this run (reviewer cross-read
  the real corpus, missed its seed) → re-run with worktree-pinned reviewer.
  Findings + remediation checklist: `debug/sprint_qa_group3bite2_memo.md`.
  Branch criteria unchanged — not relaxed.
- **Re-touch (bite 2 must re-confirm):** **Paper 24** was edited in v4.19.0
  (exchange-constant citation synced "four types" → six-tier list) *after* its
  bite-1 cert; bite 2 re-confirms Paper 24's C4/C8 on that one edit. **Paper 18**
  was edited in v4.19.0 (taxonomy reconciled to six tiers + `\Z` compile-bug
  fix) — it is in bite 2's primary scope, so this is covered.

## Branch deltas (the only non-inherited content)

- **C4 high-fabrication surface (watch).** The periods / Tannakian apparatus of
  Papers 55/56 (Deligne–Milne, Brown, Fathizadeh–Marcolli, mixed-Tate /
  cosmic-Galois literature) is the branch's highest fabrication-risk surface —
  every cited theorem/def number verified.
- **C6 (watch).** Bite-2 spectrum statements in Papers 18, 24, 54 (and 22/31
  already cleared in bite 1).
- **C7 (trunk-dependent status).** Where a group3 paper cites a trunk result:
  Paper 38/WH1 as PROVEN scoped to the van Suijlekom state-space GH distance;
  κ = −1/16 as an Observation; no overstatement of a trunk keystone.
- **C8 (headline honesty), per-paper (bite-2 papers).**
  - **Paper 18** = the exchange-constant taxonomy is **six tiers** (intrinsic,
    conformal/calibration, embedding, algebraic-implicit, composition,
    inner-factor input data — the v4.19.0 reconciliation); the master Mellin
    engine M1/M2/M3 + the three §VIII theorems stated at their tiers.
  - **Paper 54** = two-body selection rules from the tensor-product triple stated
    as structural; radial coupling **NOT forced** (the honest negative).
  - **Paper 55** = every GeoVac period = **cyclotomic mixed-Tate at level ≤ 4**
    as a *classification* (tier-appropriate), not a derivation of physics.
  - **Paper 56** = Tannakian/cosmic-Galois reconstruction as a **theorem-grade
    abelianized homomorphism (Reading A) at finite cutoff** — explicitly **NOT a
    closed immersion** (M3 column rank 2 via χ₋₄ vertex parity, capped ≪ N;
    parity-blind restriction rank-1) — with honest scope (infinite-cutoff
    equality NOT claimed). *(A reviewer presenting the injection as a genuine
    closed immersion / faithful / rank-N is a MATERIAL overstatement.)*
  - **Paper 57** = forced/free P5 packing-reachability is an **internal-consistency
    check**, not an independent discriminator — the packing-reachable tag is the F/C
    label by construction, so its **98.3%** is an agreement rate (not a discovered
    accuracy), with the I3 ambiguous case noted, not "100%" (Paper 57 §6.1). The
    genuine discriminator signal is the two-family decomposition. *(A reviewer
    presenting P5 as a genuine discriminator / "98.3% accuracy" / a principle that
    "decides forcedness" is a MATERIAL overstatement.)*
  - *(Bite-1: Paper 22 angular sparsity theorem; Paper 24 π-free Bargmann–Segal
    lattice + Coulomb/HO six-layer asymmetry; Paper 31 A/D partition — certified.)*
- **No branch-specific C14+** (no DESCOPED/PARTIAL papers in group3 scope).

## Change log
- 2026-06-15 — DRAFTED, then FROZEN (PI-confirmed). Bite 1 = Papers 22/24/31.
- 2026-06-16 — bite 1 = **PASS** (`/qa group3` run #7, v4.18.0).
- 2026-06-16 — **slimmed to a profile** (criteria → `docs/qa/criteria.md`, no
  criterion changed); **bite 2 defined + FROZEN** (Papers 18, 54–57 + full
  synthesis), PI-authorized. Paper 24 re-touch flagged.
- 2026-06-16 — **bite 2 run 1 = FAIL** (calibrated code/claims/citation/det
  dimensions; synthesis dimension uncalibrated, re-run needed). Findings +
  remediation checklist: `debug/sprint_qa_group3bite2_memo.md`.
- 2026-06-16 — **bite 2 run-1 findings REMEDIATED** (v4.20.0): all material
  defects fixed — 3 keystone corrections (P56 C4 refuted → abelianized; P40 →
  state-space GH; P57 P5 → consistency-check) + the PM-fixables; deterministic
  gate (C10–C13) GREEN. **Remediated, not re-certified** — a fresh `/qa group3`
  (synthesis dimension worktree-pinned) certifies. Memo passes 1–2 +
  post-remediation honest scope: `debug/sprint_qa_group3bite2_memo.md`.
- 2026-06-17 — **bite 2 re-cert run 2 = FAIL** (`/qa group3`, post-v4.20.5).
  Calibrated panel (sensitivity 5/5, specificity 5/5; synthesis dimension
  worktree-pinned this run). 8 verified MATERIAL defects — all *missed
  instances* of already-fixed classes: C6 graph-produces-spectrum (P18 l.149/648,
  P55 l.145); P56 "closed sub-pro-algebraic group" zombie (l.1539/1676/1805);
  P57 bare-128 (l.77) + "38 forced" (l.175); P18 α²-Ihara false-badge + κ–B
  no-test; P56 5,864-residual headline absorbing the refuted C4; P54
  angular-theorem tier. Answer key `debug/qa/group3_recert_seed_key.json`.
- 2026-06-17 — **run-2 findings REMEDIATED** (v4.20.6); re-run pending.
- 2026-06-17 — **bite 2 re-cert run 3 = FAIL** (`/qa group3`, on v4.20.6).
  Calibrated panel (sensitivity 5/5 — S5/S2/S1/S9/S7 all caught; specificity
  6/6 — M1–M6 incl. the v4.20.6 C6/κ/abelianized fixes all held). 4 reviewers
  died mid-run on a spend limit, re-dispatched after tokens restored → all 5
  dimensions ultimately exercised + calibrated (the first code-56 missed S2;
  the fresh code-56 caught it — panel redundancy recovered). The 8 remediated
  classes are confirmed FIXED (controls passed); the FAIL is a **deeper,
  pre-existing defect** the gate peeled to: **Paper 56 `thm:injection_g4`
  C1/C2 legs** — C1 multiplicativity test is tautological (`simplify(lhs −
  product)` with identical operands), C2 coproduct test is vacuous (`assert
  … is not None`), and C2's sole literature support, **"Brown 2017
  Proposition 5.2," appears nonexistent** (Brown ICM-2014 §5 = depth-filtration
  / modular-forms conjecture, eqs 5.1/5.2 + Conjecture 5.1, no Prop 5.2;
  convergent run-2 "unverifiable" + run-3 "nonexistent" + ToC topic-mismatch).
  Plus SMALL: P56 §open_g4 boxed-theorem leads "closed at theorem grade …
  closed subgroup" before the abelianized qualifier. Seed key
  `debug/qa/group3_recert2_seed_key.json`. Remediation pending PI direction.
- 2026-06-17 — **run-3 findings REMEDIATED** (v4.20.7): C2 regrounded on
  Cartier–Milnor–Moore, C1/C2 tests genuine-ified, box reworded.
- 2026-06-17 — **bite 2 re-cert run 4 = FAIL** (`/qa group3`, on v4.20.7).
  Calibrated panel (sensitivity 5/5 — S4/S2/S1/S8/S7 all caught; specificity
  6/6 — M1–M6, incl. M5 confirming the v4.20.7 C2/CMM fix is ACCEPTED, all
  held). No spend-limit deaths. The v4.20.7 fixes confirmed; FAIL is **two
  more genuine, deeper, distinct defects** (4th consecutive calibrated FAIL,
  each on different papers — severity now down to a worked-example typo + a
  broken code-ref): (1) **Paper 18 Thm 1(2) worked example D(4) = 2ζ(2)+2ζ(3)
  is WRONG** (l.1661) — that is the *Fock-index* value; Dirac D(4) = π²−π⁴/12
  (π-even, no ζ(3)) per the paper's own formula + `test_D4_is_pi_even`; the
  odd-zeta example should be D(5). The validated test contradicts the prose.
  (2) **Paper 55 `thm:jlo_depth2_reading_A` cites nonexistent
  `geovac/jlo_chi.py`** (l.1305) — Reading-A disambiguation backed only by a
  debug/ script + a broken code-module ref (a class C13 does not cover — it
  checks test refs, not code-module refs). Seed key
  `debug/qa/group3_recert3_seed_key.json`. Remediation pending PI direction.
- 2026-06-17 — **C14 added** (`debug/qa/check_file_refs.py`, paper↔file
  reference integrity; `tests/test_file_ref_check.py`): closes the C13 gap
  (gates `geovac/`/`benchmarks/`/`demo/` code-artifact refs that C13's
  tests/-only scope missed — the jlo_chi.py class). Corpus-wide serious-class
  sweep = **CLEAN** after fixing the one remaining hit
  (`geovac/hyperspherical.py`→`hyperspherical\_adiabatic.py` in paper_34);
  **group3 PASSES C14** (0 serious-class misses). **443 `debug/` dangling
  pointers logged ADVISORY** (transient clean-room dir, pruned by design —
  a separate hygiene-sweep decision for the PI, not a cert blocker).
- 2026-06-17 — **bite 2 re-cert run 5 = PASS** ✅ (`/qa group3`, on v4.21.0).
  **First calibrated PASS** after 4 remediation rounds. Calibrated panel
  (sensitivity 5/5 — S5/S2/S1/S9/S7 all caught; specificity 6/6 — M1–M6,
  incl. M4 the v4.20.8 D(4) fix and M5 the v4.20.7 C2/CMM fix, both ACCEPTED
  by the fresh blind panel). Zero verified MATERIAL defects (criteria test:
  no target result changes; `thm:injection_g4` soundly backed by the
  periodmap rank-1 driver). Residual NITs (recommended polish, non-blocking):
  slow-variant C1 test still tautological + unused imports in
  `test_paper56_injection_g4.py`; "5,864 residuals" headline conflates count
  bookkeeping with computed zeros for the C1/C2 portion; worked-Example-1 C6
  wording in P18; peripheral Paper-51 one-loop-gravity cite; NA-1/S⁵
  open-section theorems backed by debug/ not tests/. Convergence arc: keystone
  corrections → missed-instances → C2-leg → D(4)+broken-ref → test-hygiene
  NITs. **BITE 2 CERTIFIED.** Seed key `debug/qa/group3_recert4_seed_key.json`.
- 2026-07-05 — **C8 Paper 57 headline tightened** (PI direction, post-Topos-4 v4.70.0
  delta-QA): "P5 packing-reachability **discriminator** at 98.3%" → "P5 is an
  **internal-consistency check, not a discriminator** (98.3% = agreement rate, not
  accuracy; Paper 57 §6.1)". Brings the DoD into line with Paper 57's own body (the
  §6.1 tautology resolution) and the corrected field-guide seam paragraph; the
  looser "discriminator/accuracy" wording was flagged by the Topos-4 synthesis
  reviewer as the source of the field-guide overstatement. No criterion changed;
  the headline statement is made tier-honest (a deliberate, dated edit — not a
  mid-review goalpost move; bite-2 cert stands).
- 2026-07-09 — **C8 Paper 56 headline reconciled** (PI direction, post-v4.74.1
  delta-QA): "theorem-grade **closed immersion** at finite cutoff" → "theorem-grade
  **abelianized homomorphism** (Reading A), NOT a closed immersion (M3 rank 2 via
  χ₋₄ vertex parity, capped ≪ N; parity-blind restriction rank-1)". The old wording
  predated the v4.20.0 bite-2 remediation that refuted the closed-immersion reading
  (C4 refuted → abelianized); the DoD line was never synced to the certified paper
  body. A deliberate, dated edit bringing the DoD into line with the certified paper
  (same pattern as the 2026-07-05 P57 reconciliation) — not a mid-review goalpost
  move; bite-2 cert stands. Surfaced by the v4.74.1 (M3 rank-1→2 flip) delta-QA,
  which ran CLEAN-DELTA (calibrated) against the corrected framing.

---

## Re-review OWED (2026-08-22, v5.0.0) — group3

**Papers in this target changed after certification.** Logged here, at the
owning source, so the certified status is not read as covering text that
post-dates it.

- **Paper 24** (`paper_24_bargmann_segal.tex`): the two-fermion
  entanglement-rigidity corollary is **RETRACTED** — claim, mechanism and its
  INTERNAL THEOREM tier all withdrawn. Section rewritten with the corrected
  measurements (S = 0.0671/0.0716/0.0833; commutator 0.74/0.63/0.67;
  E₀ = 21.6538/21.6279/21.5442). The π-freeness results are untouched.
- C16 registry entry `p24-entanglement-rigidity` added (scope group3 + group4).
- Backing tests rewritten (`test_paper24_ho_entropy.py`), plus a new tripwire
  `test_paper24_n_tot_guard_stays_removed`.

**Discharge condition:** a CLEAN DELTA over the changed loci (diff-scoped, per
the qa.md run-shape rule), which is also the standing precondition for the next
FULL certifying pass. The delta must re-test *these specific defects* rather
than trust this entry (qa.md hard rule, added the same day).

**Deterministic layer already re-run and GREEN on this target:** C10 (compiles,
with the aux-clean fix), C13, C14, C16, C17, C18, C19, C5/C12 (now corpus-wide
after the scope fix). What is owed is the LLM-judgment layer: claims, and
synthesis where the target has one.

**DISCHARGED — CLEAN-DELTA, 2026-08-24 (v5.0.9).** Delta scope = the changed loci
(P24 `sec:entanglement-rigidity` retraction + the group3 synthesis's rigidity
passages). Pasted-hunks delta (seeds only in prompts). Two dimensions dispatched, both
**exercised + calibrated**: claims on P24 (opus, 1/1 seed — an inverted "negligibly
small" clause, caught) and synthesis on the group3 HO-rigidity passage (opus, 1/1 seed
— an H²(S⁵)→H²(S³) manifold slip, caught); specificity clean. **Zero verified
MATERIAL** — the P24 retraction is honest and internally consistent, and the group3
synthesis carries only the *valid parent* HO rigidity theorem (unique central potential
from the Euler operator on H²(S⁵)), NOT the retracted entropy corollary. Seed key
`debug/qa/group36_delta_seed_key.json`.

**Honest correction to the note above:** the v5.0.0 "deterministic layer GREEN" claim was
NOT accurate for this target — C10 (run correctly) had **pre-existing broken/cross-document
`\ref`s** in paper_18 and paper_22 (`sec:theorem`, `sec:universal_partition`, +cross-doc to
p14/p18/p28) and a bad `\cite{paper7}`, all fixed on sight this run. Deterministic layer is
green **now**. Branch stays CERTIFIED; a clean delta is the precondition for the next FULL
certifying pass.

---

## FULL certifying pass — CERTIFIED ✅ (2026-08-24, v5.0.9)

Fired after the clean delta above (qa.md run-shape rule: only a FULL run emits PASS).
Completes the directed group4 → group6 → **group3** retraction re-cert sequence.

**Panel:** 9 dispatched reviewers (4 claims chunks, 1 citation, 1 synthesis, 3 code)
+ 1 completeness-critic. **All 9 calibrated** — every planted seed caught (P18
Observation→Theorem, P56 injection-flip, P7 discrete/continuum, P24 single-Slater,
2× P55 arXiv transpositions, synthesis Fock-rigidity uniqueness-flip, 3× code
b1/b2 tautologies); zero false-positives on the known-good controls. The critic
independently re-caught the P24 seed (blind) and confirmed the P31 §9 remediation scope.
Seed key `debug/qa/group3_fullcert_seed_key.json` (scratch deleted; seed-leak check clean).

**One cert-blocking defect — found & remediated (analog of the group6 synthesis zombie):**
- **P31 §9 (`sec:sig_l2_verification`) descope-zombie** — the withdrawn Lorentzian
  "literal identification at finite cutoff / Lorentzian closure complete" reading
  survived because the C16 entry `lorentzian-literal-identification-krein` was
  scoped **group6-only** (files = P34 + group6 synth); the group6 §III.29 correction
  never propagated to P31 (group3). Remediated to the corrected structural-
  correspondence / signature-blind / compact-boost verdict (Riemannian closure kept);
  **C16 widened to `group3 group6` + P31 file + two single-line tells, discrimination
  proven** (fires on old wording, exempt on corrected).

**Fix-on-sight (all remediated):** κ tier-visibility (P31 "set by Jacobian"→coincides/
Observation; P31 π-table yes→no; P18 §kappa +Observation caveat); P22 O(Q^2.5)/51×–1712×
universality→hedged to Corollary 1; P24 Minnesota −0.81→−0.55 MeV (paper+test); P57
"magic numbers"→HO-closures; P22 ~4×→~3×; P7 abstract + P0 cross-ref + P7 N>1
discrete/continuum precision; 2 citation titles; synth date; a pre-existing C18
duration ("year of internal work").

**Backing gaps — PI directed "write the tests first" (all CLOSED with genuine tests):**
- P55 M1 (`test_paper55_m1_pure_tate.py`): rewrote the tautological F-theorem sub-test
  to call production `qed_two_loop.{scalar,dirac}_F_theorem()` (retires the resurrected
  Paper-50 typed-in-F false-positive) + added a discrimination guard (rejects ζ(3),
  ζ(3)/π, Catalan, ζ(5)).
- P55 M3-on-S⁵ (`test_paper55_m3_s5.py`, NEW): direct spectral sum vs Hurwitz closed
  form + χ₋₄ parity discriminant, ~40-digit, + wrong-degeneracy guard.
- P55 M2 Grothendieck (`test_paper55_grothendieck.py`, NEW): derives the explicit
  L-polynomial class from the defining equation, general-n symbolic identity + guard.
- P22 Breit rank-2 (`test_paper22_breit_angular.py`, NEW): reproduces the Z=4
  Coulomb/SS/SOO density table (l_max 0–3) from Wigner-3j + potential-independence + guard.
- P54 (zero tests/ cites → cited): added inline `test_paper31_two_body.py` citations to
  `thm:angular_structure` + connected-fraction; removed the `debug/` DF↔multipole citation.
- Papers now cite the tracked tests (S⁵/Grothendieck/Breit debug/ citations removed);
  `docs/claim_test_matrix.md` rows synced (54/55/56/22-Breit/24-count).

**Deterministic layer GREEN (whole-target):** C5, C10 (12 docs compile), C11, C13, C14,
C16, C17, C18, C19. **Regression:** 139 passed / 9 slow-skipped over the touched +
baseline suite (incl. the 18 symbolic S³ proofs).

**Verdict: group3 CERTIFIED (FULL pass PASS).** The v5.0.0 retraction re-review OWED is
fully discharged across group3/group4/group6.
