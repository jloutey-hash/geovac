# Sprint: /qa all-dimensions enforcement (runs #2–#4)

**Date:** 2026-06-15 · **Version:** v4.16.2 · **Verdict:** trunk PASS (after one in-protocol remediation); gate hardened.

## Why this sprint

`/qa trunk` run #2 (2026-06-14) reported **PASS** while exercising only 2 of the gate's review dimensions — **claims** + **citations**. The **code** (C1–C2) and **synthesis** (C9) dimensions were unexercised and footnoted in the "honest ceiling." PI flagged it: *"how can we ensure /qa attempts all three dimensions if preceding gates pass."* That question exposed a structural hole — a clean early dimension was short-circuiting the rest, and an unexercised dimension was being downgraded to a footnote instead of forcing INCONCLUSIVE.

## What was done

### Run #3 — ran the two skipped dimensions
- **Code dimension** (code-reviewer ×4: P7, P1, P32, P38): calibrated (S2 tautological test, S3 false-positive test both caught; controls clean). Genuine corpus tests **sound** → PASS.
- **Synthesis dimension** (claims-reviewer on group3 synthesis): calibrated (S8 zombie, S9 status-overstatement caught) AND independently surfaced **2 real defects** the paper-only runs structurally could not see:
  - **κ "derivable from the Fock projection rather than fitted"** (C8) — contradicts Paper 7's own "observation, not a derivation."
  - **K = π(B+F−Δ) "conjectural" ×4 sites** (C5 hard prohibition) — stale vs the 2026-06-14 conjecture→Observation downgrade.
- → synthesis dimension **FAIL**; the trunk's run-#2 PASS was wrong.

### Fixes (cea711f, group3 synthesis)
- κ recast → **Observation/coincidence** (geometric 1/16 is a real computation; no bridge derives κ's calibration role from it; aligns with Paper 7 §III).
- K → **"an Observation"** at all 4 sites; "never a conjecture or derivation."
- Nits: six-vs-five Coulomb/HO layer count → six; `λ_n=−(n²−1)` graph → S³ eigenvalues (C6). Recompiles ERRORS=0.

### Gate hardening (qa.md, trunk.done.md, seed_defects.md)
- Step 4 names the **review dimensions** (code / paper-claims / citations / synthesis / deterministic), reviewer types, and gated criteria; **all mandatory in one invocation**, no early termination.
- Step 7: **per-dimension scorecard + AND roll-up**; new Hard rule "all dimensions every run"; **unexercised gating dimension ⇒ INCONCLUSIVE, never PASS**.
- Step 3 + seed_defects.md: **per-dimension seed coverage**.
- trunk.done.md: Review-dimensions table.

### Run #4 — first full all-dimensions certification (on cea711f)
- 8 seeds, ≥1 per gating dimension; 7-agent panel + deterministic layer.
- **Sensitivity 8/8, specificity 0 false positives.** Synthesis re-cert **PASSED** (κ/K fixes confirmed clean; only S8/S9 flagged).
- **But the 1st claims-reviewer(P32) MISSED planted S5** ("K is now derived as a theorem", line 5532) — it *sampled* K-appearances and was fooled by the self-contradicting hedge two lines down ("does not promote K to a theorem"). Under the new rule → **INCONCLUSIVE**.
- Remediated per step 7 ("fix the reviewer prompt/agent, and re-run"): re-dispatched a sharpened claims-reviewer that **enumerates-and-quotes every K-sentence** → caught S5 exactly (left the compliant non-selection *theorem* and hypotheticals unflagged). Claims-prose dimension calibrated → **trunk PASS**.

### Durable fix from the run-#4 miss
`.claude/agents/claims-reviewer.md` step 5 now **requires exhaustive enumerate-and-quote of every K-sentence** + "a self-contradicting hedge doesn't cure a tripwire," with the run-#4 lesson cited inline. (Parallel to the run-#1 → C11 title lesson. Stronger option offered to PI: a **deterministic K-candidate flagger** to give C5 the C11 treatment — not yet built.)

## Lessons (durable)
1. **All dimensions every run.** A multi-axis QA verdict is an AND over independent dimensions; a clean axis never short-circuits another; an unexercised gating axis is INCONCLUSIVE, not a footnoted PASS.
2. **The K hard-prohibition is missable by an LLM reviewer under hedging.** Run-#2 caught S5; run-#4's first reviewer didn't (variance). Tripwires need exhaustive enumeration, not sampling — and ideally a deterministic backstop.
3. **The gate keeps catching its own blind spots** (run-#1 missed S7 → deterministic title check; run-#4 missed S5 → exhaustive-K-enumeration). The three-way verdict (INCONCLUSIVE) is what makes that visible instead of a silent false PASS.

## Artifacts
- Fixes: `cea711f` (synthesis + gate). Hardening: `.claude/agents/claims-reviewer.md`.
- Run key: `debug/qa/trunk_seed_key.json` (run #4). Ledger: `debug/honest_review_2026_06_14_ledger.md`.
