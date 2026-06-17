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
    closed immersion at finite cutoff** with honest scope (infinite-cutoff
    equality NOT claimed; Reading A).
  - **Paper 57** = forced/free P5 packing-reachability discriminator at its
    **true 98.3%** with the I3 ambiguous misclassification noted, not "100%".
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
