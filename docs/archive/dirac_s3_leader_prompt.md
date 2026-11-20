# Leader Agent Prompt — Dirac-on-S³ Tier 1 Sprint Design

**To be used by:** pasting `agents/LEADER.md` + `CLAUDE.md` + this file into a fresh Claude session.

---

## Context

The α sprint series (Phases 4B–4H, documented in CLAUDE.md §2 backlog) has produced three structural identifications for the components of the α combination rule K = π(B + F − Δ):

- **B = 42** — finite Casimir trace on S³ at m=3; κ↔B link established via Fock-weight identity (Phase 4B, α-C).
- **F = π²/6 = D_{n²}(d_max)** — Dirichlet series of Fock degeneracies at the packing exponent (Phase 4F, α-J).
- **Δ = 1/40 = g₃^Dirac(S³)** — third single-chirality Dirac mode degeneracy on unit S³ (Phase 4H, SM-D).

The third identification is the key: it puts **Dirac structure inside the α combination rule**, without that having been the hypothesis. Seven non-Dirac mechanisms have been eliminated (Phases 4B–4G plus SM-A/B/C/F). The natural next research direction is to formalize the Dirac-on-S³ sector and ask whether the combination rule becomes a theorem rather than a three-tier coincidence.

This is also the natural next step for the quantum-computing track (Papers 14, 20): a spin-ful / relativistic upgrade to the composed qubit encoding is a bounded refinement (qubit count ~2× or 4×, Pauli scaling ~Q^2.5 preserved, JW mapping unchanged) that opens the heavy-atom relativistic-chemistry regime where Gaussian methods are weakest.

The PI has committed to a three-tier program:

- **Tier 1 (months 1–2):** Dirac-on-S³ formalized; Fock rigidity theorem extended to Dirac; α combination rule revisited as potential theorem.
- **Tier 2 (months 3–4):** Spin-ful composed qubit encoding; fine-structure splittings for He/Li/Be.
- **Tier 3 (months 5–6):** Heavy-atom (Au or Hg) qubit Hamiltonian with relativistic corrections baked into the geometry.

This prompt asks the Leader to scope **Tier 1 only**.

## Your task

Produce a Strategic Brief for **Tier 1: Dirac-on-S³ formalization**. Specifically:

1. **Minimal viable sprint.** What is the smallest well-defined sprint (3–6 tracks, 2–4 weeks of work) that could either (a) promote the α combination rule K = π(B + F − Δ) from a three-tier coincidence to a theorem, or (b) cleanly falsify that possibility and identify what would be needed instead? "Cleanly falsify" is an acceptable outcome — we want a decisive answer, not a partial one.

2. **Required infrastructure.** What mathematical objects must be constructed first before the combination-rule question can be posed precisely? At minimum I expect: the Dirac operator on unit S³, its spinor spherical harmonics indexed by (n,l,m,spin) compatible with our existing (n,l,m) graph labels, its spectrum with exact degeneracies, and a version of the Fock projection that maps the Dirac-Coulomb problem on R³ to free Dirac on S³. Identify any additional objects, and for each say whether it's off-the-shelf (cite the reference) or needs to be built.

3. **Ranked tracks.** Propose 3–6 tracks ranked by expected information gain per unit effort. For each track give: (a) one-sentence goal, (b) what success looks like, (c) what failure looks like, (d) what papers/memos it would produce or affect, (e) dependency on other tracks.

4. **Guardrail check.** Which existing Section 3 (Approaches That Failed) entries are adjacent to this sprint? The SM-origin sprint (Phase 4H) is the most obvious. Are there any non-Dirac α mechanisms that a naive "just extend to Dirac" approach would accidentally re-derive? Flag them so the track prompts can avoid them.

5. **What the Leader does NOT do.** Do not propose Tier 2 or Tier 3 work. Do not propose writing production code. Do not propose paper drafts. You propose the sprint; the PI approves it; the PM decomposes it into sub-agent dispatches; the workers execute.

6. **Open questions for the PI.** At the end, list 2–4 questions whose answers would materially change the sprint design. Example shape: "Should the Dirac operator be two-component (Weyl) or four-component (full Dirac) at Tier 1? The α identification uses single-chirality modes, but the Fock rigidity extension may need the full Dirac to match the Coulomb spectrum."

## Output format

```
# Strategic Brief — Dirac-on-S³ Tier 1

## One-paragraph summary
[what the sprint accomplishes and why now]

## Infrastructure required
[numbered list; for each, off-the-shelf vs. to-build]

## Proposed tracks (ranked)
Track D1: [name] — [one-line goal]
  Success: ...
  Failure: ...
  Deliverable: ...
  Depends on: ...
Track D2: ...
[etc.]

## Guardrail adjacency
[which Section 3 entries border this sprint, and how to not re-derive them]

## Open questions for the PI
1. ...
2. ...
```

Keep the brief under 1500 words. Prioritize actionable scoping over literature review — the Explorer agent will be invoked separately for the literature pass.
