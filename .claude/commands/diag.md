---
description: Diagnostic-before-engineering — force a ~1-day diagnostic step before launching an implementation sprint
---

Stop. Before launching an implementation sprint, run a diagnostic pass. Reference:\ `memory/feedback_diagnostic_before_engineering.md`. The rule:\ when ≥ 2 honest negatives accumulate on the same wall, do a diagnostic-only sprint before another engineering attempt. Many recent W1c-arc sprints (F2 kernel substitution, F4 bonding PK, F5 explicit-core Hartree, F6 basis enlargement, Schmidt, core correlation) closed in hours-to-days because diagnostic pre-passes ruled them out cheaply.

**Step 1. Name the closure target.** What specific quantity, threshold, or behaviour would the implementation be trying to achieve? State it as a decision gate (GO / BORDERLINE / STOP) with explicit thresholds.

**Step 2. Name the structural falsifier.** What would have to be true (about the underlying physics / math / architecture) for the implementation to work? State it as a single-sentence claim that can be tested cheaply.

**Step 3. Cheap test.** What is the smallest computation that exercises the falsifier? Examples from recent practice:
- Algebraic check (sympy identity, closed-form differential)
- 1–2 point computation at the smallest tractable system
- Step-1-of-3-step gate (algebraic → eigenvalue → mini-PES, with each step gating the next)
- Diagnostic decomposition (e.g.\ per-config CI contributions, per-mechanism residual breakdown)

**Step 4. Predict before computing.** What do you expect the cheap test to return if the falsifier holds? What do you expect if it doesn't? Pre-commit this prediction — write it to memo / data file BEFORE running the test, per the W3 protocol.

**Step 5. Run the cheap test.**

**Step 6. Verdict.**
- GO:\ falsifier survives, implementation sprint is worth the cost.
- STOP:\ falsifier killed, implementation would not have closed the gate.
- REFINE:\ falsifier partially holds, but the closure mechanism is structurally different from what was originally proposed. Refine the closure target before any implementation.

**Step 7. If STOP or REFINE — name the new closure target.** What is the next thing worth attempting? Or is this wall genuinely structural and worth recording as a closed dead end (§3) / multi-focal-composition-wall instance?

**Discipline.** The user's "diagnostic-before-engineering" rule has converted multi-week sprints into one-day diagnostics multiple times (NaH max_n=3 prediction test, Phase 1B-E interior closure, L3b-2a candidate validation). Pay the diagnostic cost; it dwarfs the cost of the implementation that does not close.
