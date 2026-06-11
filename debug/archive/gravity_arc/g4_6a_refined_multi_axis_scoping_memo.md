# Sprint G4-6a refined — multi-axis exploration scoping

**Date:** 2026-05-29
**Path:** Multi-task thread 9, Track β''. Honest scoping memo for the genuinely multi-month G4-6a refined sub-sprint.
**Verdict:** **GO-WITH-HONESTLY-REVISED-SCOPE.** The N_φ-axis clean negative (thread 8 Track a') and the spectral-slope re-examination finding (this thread Track α'') together indicate that G4-6a refined is **genuinely multi-month, requires multi-axis substrate exploration, AND requires a prior audit of v3.19.0 driver conventions before committing**. Three named scoping options. Honest verdict: continuing G4-6a refined toward theorem-grade UV recovery is **not the highest-value next move at this point** given the day's findings; **alternative G4-6 closure paths via scope-limitation are equally viable**.

## 1. The starting state

After threads 1-8 + thread 9 Track α'':

**Closed sprint-scale:**
- G4-6d (spectral azimuthal): bit-exact F6, production tests pass, 160× UV improvement over FD.
- G4-6b (IR boundary): B_substrate = 0.163 at R≥10, within 2.3% of continuum +1/6.

**Substantive new structural findings:**
- N_φ-axis refinement (thread 8 Track a'): A recovery DEGRADES with N_0 — single-axis refinement does NOT close A extraction.
- Spectral-FD substrate-discretization invariance (thread 9 Track α''): both substrates give identical slope at α=2 (within 4 decimal places).
- N_ρ-axis stability (thread 9 Track α''): slope CV = 3.4% across N_ρ ∈ {100, 200, 400, 600}.
- **Sign discrepancy with v3.19.0 Möbius (thread 9 Track α''): direct slope measurement at α=2 gives +0.052, opposite sign from v3.19.0 -0.056 measurement. Substrate-level Möbius identification needs re-examination.**

## 2. The G4-6a refined open question

The original G4-6a sub-sprint asked: can the discrete substrate's measured A coefficient converge to A_cont = 1/(24π) via systematic substrate refinement?

After the day's findings, the question structurally splits:

**Q1.** Does the substrate's A coefficient converge to A_cont under MULTI-AXIS refinement (vary $a$, $N_\rho$, $N_\phi$ jointly)?

**Q2.** Or is the substrate's UV content structurally limited at substrate values reachable in compute budgets?

**Q3.** Or is the "A coefficient" itself a substrate-specific observable that doesn't correspond to the continuum target (A_cont) in a clean way?

The N_φ single-axis result (Track a' thread 8) suggests Q2 or Q3 is the likely answer. The honest multi-axis exploration would test all three:

- If Q1 holds (multi-axis convergence): G4-6a refined closes at multi-month effort.
- If Q2 holds (compute-budget limit): G4-6a refined needs alternative observable or substrate.
- If Q3 holds (observable mismatch): the original A_cont target may not be the right closure goal.

## 3. Three scoping options

### Option 1 — Multi-axis substrate sweep (genuine multi-month commitment, 3-6 months)

**Approach:**
- 4-axis sweep: $(a, N_\rho, N_\phi, t)$ jointly varied at substrate-feasible refinements.
- Substrate sizes: $a \in \{0.025, 0.05, 0.1\}$, $N_\rho \in \{100, 200, 400, 800\}$, $N_\phi \in \{60, 120, 240, 480\}$, $t \in \{a^2, 2a^2, 5a^2, 10a^2, 50a^2\}$.
- 192 substrate cells × 5 t-values = 960 heat-trace computations.
- Per-cell compute scales as $O(N_\phi · N_\rho^2)$ — largest cell (a=0.025, N_ρ=800, N_φ=480) is ~10 minutes.
- Total estimated compute: ~10-20 hours wall on a single workstation.

**Expected outcomes:**
- Identify whether A recovery improves under multi-axis refinement (Q1 closes positively).
- Or characterize the structural limit (Q2 or Q3).

**Risk:** even with 4-axis sweep, the substrate may not converge to A_cont. The compute is substantial but the closure is uncertain.

**Effort:** 3-6 months including the data analysis, multi-axis Richardson extrapolation methodology, paper documentation.

### Option 2 — Audit v3.19.0 + scope-limit Paper 51 (sprint-scale, ~1-2 weeks)

**Approach:**
- Audit `debug/alpha_gt_1_ansatz_test.py` and related v3.19.0 drivers to identify the convention difference behind the Möbius sign discrepancy.
- Resolve whether the substrate-level Möbius identification in Paper 51 is correct.
- Update Paper 51 with honest scope: either re-affirm the Möbius (after convention clarification) or remove the substrate-level identification and document the N_ρ-stable +0.052 slope behavior instead.

**Expected outcomes:**
- Clean Paper 51 documentation reflecting actual substrate behavior.
- Möbius mechanism question reframed honestly (either supported or removed).

**Risk:** the audit may reveal that the v3.19.0 substrate-level identification was an artifact, requiring substantial Paper 51 revision.

**Effort:** 1-2 weeks for audit + Paper 51 update.

### Option 3 — Accept observation-rigor closure + scope-limited Paper 51 (sprint-scale, ~1 week)

**Approach:**
- Accept the day's findings as the realistic G4-6 closure state.
- Update Paper 51 to reflect: (i) two sub-sprints closed (G4-6d + G4-6b); (ii) substrate-level Möbius identification flagged for re-examination; (iii) honest acknowledgment that A coefficient extraction at substrate values is multi-axis-dependent and not closed.
- Do not pursue Option 1's multi-month commitment.
- Move on to other research directions.

**Expected outcomes:**
- Paper 51 documents the gravity arc at observation rigor + the structural-skeleton-scope reading sharpened by today's findings.
- G4-6 multi-month commitment is honestly scope-limited.

**Risk:** the theorem-grade closure that the original G4-6 plan envisioned is not delivered. The compensating finding is the structural-skeleton-scope sharpening (compression-pattern-has-limits) which is publication-grade in its own right.

**Effort:** ~1 week for Paper 51 update + closure narrative.

## 4. Honest recommendation

**Option 2 is the highest-value immediate next move.**

Three reasons:

1. **The substrate-level Möbius identification flag is load-bearing.** If Paper 51 currently documents an incorrect substrate-level identification (per the Track α'' sign discrepancy), this should be resolved before any further G4-6 work. The audit is sprint-scale and the resolution is structurally definitive.

2. **Option 2's outcome informs Option 1's viability.** If the audit reveals that v3.19.0 was using a convention that genuinely produces the Möbius slope, the substrate-level mechanism stands and Option 1 multi-month commitment is justified. If the audit reveals the Möbius was a convention artifact, Option 1's foundation is undermined and Option 3 becomes the cleaner closure.

3. **Day's pattern continues to favor closing structural questions before committing to multi-month work.** The compression pattern that worked for G4-6d and G4-6b was: do the focused diagnostic first, see what the substrate actually does, then commit. Option 2 follows this discipline.

**Option 1 (multi-month G4-6a refined) is recommended only AFTER Option 2 resolves the substrate-level Möbius status.**

**Option 3 (scope-limited closure) becomes the cleanest path if Option 2 reveals the Möbius substrate-level identification was an artifact.**

## 5. Decision gate for Option 1 (multi-axis G4-6a refined)

If the PI commits to Option 1, the decision gate is:

- After 2-3 months of multi-axis sweep, can the substrate's A coefficient be Richardson-extrapolated to A_cont with sub-10% error?
- If YES → G4-6a refined closes; continue to G4-6e (Mellin moment theorem grade) and G4-6f (final).
- If NO → revisit Q3 (observable mismatch); the substrate may not have the UV content; scope-limit Paper 51 closure.

This is a real fork. The compute budget is substantial; commit only with PI direction.

## 6. Connection to the day's compression-pattern findings

The compression pattern from threads 1-7 worked at 10-50× ratios for tractable structural questions. Thread 8 (N_φ-sweep) and Thread 9 Track α'' (slope re-examination) reveal where the pattern hits limits:

- **Compresses well:** problems that decompose into focused diagnostics tractable in main session (G4-6d, G4-6b).
- **Doesn't compress:** problems that require multi-axis substrate exploration to reveal structural behavior (G4-6a refined Q1) OR audit of previous convention choices to disambiguate measurements (the Möbius slope question).

The honest framework-state reading is:
- **Structural-skeleton work compresses.** The framework's selection rules, transcendental classification, propinquity convergence, etc. — all are tractable in main session.
- **Calibration / direct numerical convergence at substrate values does NOT compress.** Multi-axis substrate exploration genuinely takes multi-month effort.

## 7. Honest scope

This Track β'':
- **Documents** the multi-axis G4-6a refined option honestly.
- **Names** three scoping options (1/2/3) with effort estimates and risk assessments.
- **Recommends** Option 2 (audit v3.19.0 + scope-limit Paper 51) as the highest-value immediate next move.
- **Identifies** Option 3 as the cleanest closure path if Option 2 reveals an artifact.

Does NOT:
- Execute Option 1, 2, or 3.
- Make the PI strategic choice.

## 8. Cross-references

- `debug/g4_6a_refined_v3_nphi_sweep_memo.md` — thread 8 Track a' (N_φ clean negative)
- `debug/sprint_moebius_reading_b_nrho_sweep_memo.md` — thread 9 Track α'' (slope re-examination)
- `debug/g4_6_scoping_memo.md` — original G4-6 plan (needs update with this thread's findings)
- Paper 51 §subsubsec:g4_5_v3_20_followon — Möbius documentation flagged for re-examination

## 9. Files

- `debug/g4_6a_refined_multi_axis_scoping_memo.md` (this)
- No driver (scoping work only)
