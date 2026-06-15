# GeoVac Authoring Conventions

The standards every paper is checked against in the **§9 Branch QA Review Protocol** (the conformance lens on the paper-adversarial reviewer). Project-wide rules apply to all papers; per-group sections add the audience-specific norms ("a lot of different areas with different audiences"). Reviewers use this as a checklist; authors use it as a style guide.

## Project-wide (all papers)

1. **Rhetoric (§1.5).** Dual-description framing; no ontological-priority claims (neither the graph nor the continuum is asserted "more fundamental"). Lead with concrete computational advantages, not interpretation.
2. **Discrete vs continuum precision.** State clearly whether a claim is about the **discrete graph** or its **continuum limit**. The −(n²−1) spectrum is a *continuum* (S³ Laplace–Beltrami) property; the discrete graph Laplacian is positive-semidefinite and *converges* to it (degeneracy recovery), it does not *have* that spectrum. Do **not** write "the graph produces the spectrum" or "the graph IS the Schrödinger equation" when the code shows convergence/matching. *(Added 2026-06-14 from the trunk code review — the framework's most code-exposed framing.)*
3. **Benchmarking.** Compare to the strongest available baseline (cc-pVTZ+ for atoms, explicitly correlated for molecules), not just STO-3G; if the comparison is unfavorable, say so and state what the framework offers instead.
4. **Conjectural labels (§13.5).** K = π(B+F−Δ) stays CONJECTURAL. Any numerically-matched-but-underived identity is an **Observation**, not a derivation.
5. **Transcendental tagging.** Every π / ζ / Catalan G / log is tagged to its Paper 18 tier + Paper 34 projection + master-Mellin M1/M2/M3 slot.
6. **Claim → artifact.** Every load-bearing claim maps to a backing test in `docs/claim_test_matrix.md`. **"derived" / "SYMBOLIC PROOF" requires a test that proves the derivation** — a matching-condition or convergence test backs "matched/converges," not "derived."
7. **Citations.** Every `\cite` resolves to a real publication that says what we claim; verified by the citation-grounding pass (`.claude/agents/citation-reviewer.md`).
8. **Status vocabulary.** Use the claims-register tiers accurately (SYMBOLIC PROOF / MEASURED / PANEL-VERIFIED / INTERNAL THEOREM / CONDITIONAL / OBSERVATION / CONJECTURE / RETRACTED). "SYMBOLIC PROOF" is reserved for exact symbolic/machine-verified derivations.
9. **Provenance visibility (§9 QA principle 1).** Every load-bearing claim wears its rule-8 tier *inline* — next to the claim in the paper, not only in `docs/claims_register.md` — so a reader (or future-self) cannot mistake an observation for a theorem. The prose may assert no more than the tier: "derived" / "SYMBOLIC PROOF" requires a *derivation* test (rule 6), never a matching / convergence test. The κ episode is the cautionary case — a prose promotion ("no longer fitted") with no tier and no artifact survived ~50 versions. The code-reviewer flags any load-bearing claim whose prose exceeds its backing, in either direction (an *under*-stated claim is upgraded — §9 QA principle 3).

## Per-group (audience-specific)

### group1 — Operator algebras / NCG (read as a math.OA referee)
- Theorems with explicit hypotheses; distinguish **proven / finite-cutoff-panel / imported**.
- Correct NCG lineage (Connes, van Suijlekom, Latrémolière, Marcolli) — exact theorem names/numbers, verified by the citation pass.
- **vS state-space GH distance** vs **Latrémolière propinquity** stated precisely (don't claim propinquity when the proved object is the state-space distance).

### group2 — Quantum chemistry (read as a skeptical quantum chemist)
- Strongest-baseline benchmarking; report error vs reference **and the reference source**.
- Honest negatives (the chemistry walls, W1e) stated as such.

### group3 — Foundations (read as a skeptical mathematical physicist)
- Rule 2 (discrete-vs-continuum precision) is the load-bearing one here.
- **Forced vs consistent:** distinguish "forced by the packing axiom" from "matches known physics" (WH3 falsifier).

### group4 — Quantum computing (read as a NISQ/VQE practitioner)
- Resource counts in standard units (Pauli terms, qubits, 1-norm, QWC groups); **matched-qubit vs matched-accuracy** caveat co-stated.
- NISQ/VQE positioning, not fault-tolerant QPE.

### group5 — QED / gauge / SM (read as a HEP referee)
- α combination rule conjectural; gauge structure dual-description.

### group6 — Precision / observations (read as an AMO precision physicist)
- Reference precision and its source stated; focal-length decomposition for multi-component observables (Lamb shift, hyperfine, fine structure).
