---
description: Classify a new transcendental appearance under Paper 18 + Paper 34 taxonomies
---

A transcendental (π, π², ζ(3), Catalan G, log 2, etc.) has appeared in a computation. Tag it under both the Paper 18 transcendence hierarchy AND the Paper 34 projection inventory. This enforces the "transcendentals must be pinned to a projection" rule from CLAUDE.md §4.

**Step 1. State the appearance cleanly.** Where did the transcendental show up? In what observable? At what coefficient? Reference the computation file / paper section.

**Step 2. Paper 18 tier classification.** Six tiers as of v2.29.0:

- **Intrinsic** (rational / ℤ / ℚ — π-free; from graph topology, selection rules, Gaunt 3j)
- **Calibration** (π, π², π³, … — from Hopf-base measure, Seeley-DeWitt, conformal projection)
- **Embedding** (1/r₁₂, cusp content — from continuum geometry the discrete graph cannot absorb)
- **Algebraic-implicit** (μ(R) polynomial roots, Cardano radicals — algebraic over ℚ(π) but defined implicitly)
- **Composition** (PK pseudopotential — from cross-block ERIs that the composed-geometry factorization cannot reach directly)
- **Inner-factor input data** (Yukawa Dirichlet ring — parameter-tied, categorically disjoint from the other five)

Which tier? Justify in one sentence.

**Step 3. Master Mellin engine sub-mechanism.** If the transcendental is in the calibration tier, identify which of M1/M2/M3 it comes from (per Paper 32 §VIII case-exhaustion theorem):

- **M1 (Hopf-base measure)** — k = 0, π or 2π or Vol(S²)/4 = π type
- **M2 (Seeley-DeWitt)** — k = 2, √π·ℚ ⊕ π²·ℚ ring
- **M3 (vertex-parity Hurwitz)** — k = 1, quarter-integer Hurwitz / Dirichlet L / Catalan G

Per Paper 18 §III.7, master Mellin engine is `𝓜[Tr(D^k · e^{-tD²})]` with k ∈ {0, 1, 2}.

**Step 4. Paper 34 projection chain.** Which of the 28 projections (§III.1 through §III.28) does this transcendental traverse? State the chain explicitly:\ source layer → projection 1 → projection 2 → … → observable. The transcendental should be attributable to a specific projection in the chain.

**Step 5. Three-axis tag.** Per Paper 34 §IV:
- **Variable** axis (which physical variable does the chain promote / project?)
- **Dimension** axis (what is the dimension of the observable?)
- **Transcendental** axis (which ring does the value live in?)

**Step 6. Catalogue update.** If the appearance is new (not yet in Paper 34 §V or §V.B), draft a row for the empirical-matches catalogue. If it's an off-precision match, tag with the error source (T / B / A / C / S per §V.B). If it's a 16th–28th projection candidate that doesn't fit, flag for §VIII open-question review.

**Step 7. If the transcendental cannot be pinned to a Paper 34 projection — STOP.** Anonymous transcendentals are not allowed in production code or papers (CLAUDE.md §4). Either (a) extend Paper 34 with a new named projection (requires PI direction), (b) identify a missing algebraic structure that would absorb the transcendental, or (c) record it as an open-question entry. Do not let it sit unattributed.
