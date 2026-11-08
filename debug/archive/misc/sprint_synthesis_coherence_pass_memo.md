# Sprint memo — synthesis-layer coherence pass (group1 + group3 + field guide + Paper 18/24)

**Date:** 2026-06-16. **Version:** v4.19.0.
**Verdict:** The three existing synthesis documents and the foundation papers they lean on were brought into internal consistency, primarily around three drift classes the recent QA never reached: the **Lorentzian descope** (group1), the **κ "derivable" overclaim + stale bookkeeping** (field guide), and the **exchange-constant taxonomy count** (Paper 18, with a pre-existing compile bug fixed on the way). group3 synthesis audited **clean** (recent first-bite cert held). All edited `.tex` compile clean (`-halt-on-error`, exit 0).

> Companion: this was a **coherence pass** (step 1 of the §9 branch-QA protocol — "synthesis update"), NOT a `/qa` certification. No calibration seeds, no fresh-adversary panel verdict. Two independent `claims-reviewer` agents did the group3-synthesis + field-guide diagnostics; PM verified every finding against primary text before acting.

## 1. Scope & method

Triggered by the group1 synthesis coherence work (the Lorentzian arc was the document's weakest part and we had just reconstructed its current state in conversation). Extended to "the same review" for group3 and the trunk. Finding: the trunk has **no separate synthesis paper**; its story lives in `geovac_field_guide.tex` (the KEYSTONE front-door doc), so the trunk review = the field-guide review.

Method per document: read current state → (for group3 synth + field guide) dispatch a `claims-reviewer` agent with the defect taxonomy used on group1 → PM verifies each hit against primary text → fix in-place → compile-check.

## 2. group1 synthesis — Lorentzian descope coherence (8 edits)

The doc had been patched at the **top** (abstract + a June-2026 Status note) but the descope was only ~80% threaded; the holes were in the worst places.

- **Self-contradiction:** the "Honest scope" subsection called strong-form Lorentzian propinquity **open**, while the §46 closure and the outlook declared it **closed by Paper 46**. Same document, opposite verdicts.
- **Stale tail (worst offender):** the outlook's bold header read *"Strong-form Lorentzian propinquity (Q1) — CLOSED by Paper 46"* with "the free-upgrade reading holds" — the reader's last impression was a retracted theorem presented as a triumph.
- **Stale Headline 4 title:** "first Lorentzian propinquity convergence theorem" (its own body retracts it).
- **Abstract overclaim:** Papers 48/49 described as triumphant closures (body descopes them); twin-paradox deficit miscredited to "Uhlmann's relative-entropy monotonicity" — the body correctly uses **Datta's max-divergence chain inequality**, and the project specifically learned ordinary relative entropy does NOT satisfy a chain inequality.

Fixes (8): Headline-4 title → "closed by a degeneracy theorem"; abstract 48/49 reframed to descoped + Uhlmann→Datta; top Status-note hedge tightened (body is now coherent throughout); §46 closure descope note + "closes"→"claimed to close"; §47 closure note (outer norm-resolvent arrow survives, inner propinquity arrow descoped); outlook lead-in (convergence credited to Paper 38, Lorentzian convergence open); outlook Q1 "CLOSED"→"OPEN (Paper 46 descoped)"; outlook state-side softened.

## 3. group3 synthesis — audited clean (0 edits)

`claims-reviewer` found **zero** high-severity defects: no descoped-as-live, no internal contradiction, no κ/K/discrete-continuum violations, tail clean. The K=π(B+F−Δ) hard prohibition handled correctly throughout (every appearance "candidate"/"Observation at the rule level"); κ correctly "coincides…an Observation, not a derivation"; six-layer asymmetry + six-tier taxonomy correct. This is the doc that just went through 7 `/qa` rounds — the cert held. Only two optional wording NITs, not actioned.

## 4. field guide — front-door fixes (8 edits)

Handled the *scary* stuff perfectly (Lorentzian descope told correctly; K hard-prohibition respected everywhere). The defects were the *quiet* drift:

- **κ overclaim, in the keystone identity section:** "κ = −1/16, **derivable** from the packing construction" — contradicts Paper 0 itself ("coincides…an observation, not a derivation") and the v4.13.0 κ→Observation downgrade. → "coincides with a conformal-geometry quantity — a numerical observation, not a derivation."
- **Discrete-vs-continuum:** "the graph Laplacian…has eigenvalues −(n²−1)" → reworded to point at the conformal equivalence (the next section already explains it).
- **Stale tail (again):** gravity entry "S_BH = A/4 derivation **pending** the L6 sub-sprint (**closed** in May 2026)" → "closed by the L6 sub-sprint."
- Stale paper-counts (Groups 1/3/5 advertised 14/7/7 → 16/11/11); "Papers 38–49 = eleven" miscount → "Papers 42–49" matching the cite list; garbled U(1)/Paper-25 cross-ref; stray `''`; dropped "independent" from "18 symbolic proofs" (Paper 7's own caveat).

## 5. Paper 18 / Paper 24 — taxonomy reconciliation + compile-bug fix

The cross-file taxonomy count was inconsistent: Paper 24 said **four** types, Paper 18 §IV intro/closing prose said **five** (old names "calibration"/"flow"), Paper 18 abstract + line 2141 ("a sixth named tier") + the developed body + the certified group3 synthesis said **six**. **Six is canonical** (the paper's own committed abstract). Resolved at the owning source (Paper 18), not deferred — per the PI principle "deferral is churn":

- Paper 18 §IV: intro prose → six tiers with canonical names; **added the missing 6th enumeration item** ("Inner-factor input data"), mirroring the paper's own inner-factor definition (lines 2144–2150) and its honest non-selection framing (admits Higgs structurally, does NOT autonomously select Yukawas; external calibration data disjoint from the M1/M2/M3 outer-factor periods); closing prose → canonical names + inner-factor; the definition's disjointness list → canonical names.
- **Pre-existing compile bug fixed on sight:** Paper 18 line 1052 used `\Z` (undefined control sequence; every other instance uses `\mathbb{Z}`). It was halting the entire paper's compile. Fixed → Paper 18 now compiles `EXITCODE=0`, which also confirms the taxonomy edits are error-free.
- Paper 24: stale "four types: intrinsic, calibration, embedding, flow" citation → Paper 18's current six-tier list.

**Deliberately left visible (not silently changed, not vaguely deferred):** Paper 18 line 2065's heading calls the inner factor "a fourth tier" — that is the *Mellin-engine* counting scheme (M1/M2/M3 + inner = 4), a different classification from the six-tier exchange taxonomy. Internally correct in its own scheme; "tier" doing double duty is a wording-judgment for the author, flagged with the exact line for the full group3 cert.

## 6. Honest scope

- **No theorem-grade or numerical-observation claims were created or changed.** This sprint moved *prose tiers and bookkeeping into line with already-established results* — no physics, no new equations, no new benchmarks.
- **Closed (coherence):** group1 synthesis now tells the Lorentzian descope consistently end-to-end (abstract/body/outlook agree); field guide front-door drift removed; Paper 18's six-tier taxonomy is internally consistent (intro/enumeration/closing) with its own abstract; Paper 24's citation synced; INDEX + CLAUDE.md §2 layer count synced to Paper 24's authoritative six.
- **Fixed bug:** Paper 18 `\Z` (was blocking compile).
- **Hard-prohibition status:** NOT violated — two near-violations *fixed* (INDEX Paper-2 K "conjectural"→Observation; field-guide κ "derivable"→observation).
- **Named open follow-ons (NOT closed here):**
  1. Paper 18 line 2065 "fourth tier" vs "sixth tier" — the two "tier" schemes (Mellin vs exchange-taxonomy) want disambiguation; for the full group3 cert.
  2. Full group3 cert (Papers 18/54/55/56/57 + full synthesis) still pending; this pass touched 18/24 only on the taxonomy axis.
  3. The structural anti-churn fix for cross-file duplicated facts (single-source-of-truth and/or a deterministic count-consistency check — the ratchet, cf. `check_internal_titles.py`) — proposed, not built.
  4. group1 synthesis is coherence-passed, not `/qa`-certified; the field guide likewise.
- **Process note:** "try again" mid-sprint was a failed-API retry, misread as a critique → a memory was deleted and re-created; net no harm, the memory `feedback_deferral_is_churn.md` stands.
