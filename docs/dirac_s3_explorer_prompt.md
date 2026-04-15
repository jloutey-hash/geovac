# Explorer Agent Prompt — Dirac-on-S³ Infrastructure Scan

**To be used by:** pasting `agents/EXPLORER.md` + `CLAUDE.md` (sections 4–6 and 12 are the most relevant) + this file into a fresh Claude session. Run in parallel with the Leader prompt (`dirac_s3_leader_prompt.md`).

---

## Context

GeoVac is a spectral-graph-theoretic framework for computational quantum chemistry built on the Fock projection from the hydrogen atom in R³ to the free Laplacian on the 3-sphere S³ (Paper 7, 18 symbolic proofs). The framework is currently non-relativistic (Schrödinger). We are extending to Dirac-on-S³ for two reasons:

1. **α combination rule (Paper 2, conjectural):** Phase 4H identified Δ = 1/40 = g₃^Dirac(S³), the third single-chirality Dirac mode degeneracy on unit S³, as a component of the α formula K = π(B + F − Δ). This puts Dirac structure inside the combination rule.

2. **Quantum-computing heavy-atom regime:** Relativistic basis methods (Dirac-Hartree-Fock, X2C) are the dominant approach for Z > 50. The composed qubit encoding needs a spin-ful / relativistic upgrade to enter this regime.

Before scoping the Tier 1 sprint, we need to know what's off-the-shelf in the mathematical-physics literature. The goal of this Explorer pass is to **minimize the amount of foundational infrastructure GeoVac has to build from scratch**.

## Your task

Produce a structured literature and infrastructure scan covering the following domains. For each find, give: full citation, one-paragraph summary of what's usable, and a compatibility assessment with GeoVac's (n,l,m) graph labeling and Fock projection machinery.

### Domain 1: Dirac operator on S³ (free)

Camporesi–Higuchi is already cited in CLAUDE.md (Phase 4H, SM-D). Confirm it and find:

- **Spectrum and degeneracies** — the |λ_n| = (2n+3)/2, g_n = 2(n+1)(n+2) result used in Phase 4H. Verify with a second independent source.
- **Spinor spherical harmonics on S³** — explicit basis functions, preferably in a form indexed by (n, l, j, m_j) or (n, l, m, spin) compatible with our (n, l, m) graph nodes.
- **Chirality decomposition** — the single-chirality Weyl vs. full 4-component Dirac distinction. Phase 4H uses single-chirality; is that the right object for an α theorem, or does Fock rigidity require the full Dirac spectrum?
- **Matrix elements of standard operators** (r, 1/r, ∂_r, angular momentum) in the spinor basis — if computed in the literature, saves a lot of Decomposer work later.

### Domain 2: Dirac-Coulomb on curved space (specifically, the Fock analog)

The original Fock result (1935) mapped the non-relativistic hydrogen atom in R³ to free Laplacian on S³ via stereographic projection with p₀ = √(−2E). Is there a known analog for Dirac-Coulomb? I.e., does the relativistic hydrogen atom with its Sommerfeld fine-structure spectrum map to free Dirac on S³ (or some deformation of S³) under a similar conformal projection?

Literature to probe:
- **Dirac-hydrogen accidental symmetries** — is there an SO(4,1) or SO(2,2) that plays the role SO(4) plays for Schrödinger-hydrogen? The "Johnson-Lippmann operator" and the "Biedenharn-Temple" operator are two classical references.
- **Dirac-Fock-Sommerfeld projection** — papers by Biedenharn, Martin, Pratt, Rose from the 1960s–80s, and more recent work by Khachidze-Khelashvili or Zhang-Silbar. Does such a projection exist cleanly, or is it only partial?
- **What deforms** — if the Fock analog exists, what replaces p₀ = √(−2E)? Likely something with both E and mc² in it. And what's the resulting "sphere" — still S³, or a deformation?

### Domain 3: Dirac spectra on other spheres (S⁵, S⁷)

Phase 4 already found that the Coulomb/HO asymmetry is structural (Paper 24). The HO lives on S⁵ via Bargmann-Segal. For the Dirac step:
- Dirac on S⁵ and S⁷ spectra — is there a pattern in single-chirality degeneracies g_n^Dirac(S^{2k+1}) that extends the SM-D identification?
- Any connection between Dirac-on-S³ and Hopf bundle structure used in Paper 2 — the U(1) fiber of the Hopf map and the U(1) gauge group of electromagnetism are the same group; is this exploited anywhere in the mathematical-physics literature as a "Kaluza-Klein"-style mechanism for α?

### Domain 4: Graph-theoretic / discrete Dirac operators

Our framework is fundamentally a graph object. The non-relativistic graph Laplacian replaces ∇². What's the analog for Dirac? Literature to probe:
- Combinatorial Dirac operators (Friedman, Tillich, and others)
- Dirac-Kähler discretization on simplicial complexes
- Lattice fermion formulations (staggered, Wilson, overlap) — these are from lattice QCD but the discrete-geometry lessons may transfer
- Specifically: is there a discrete Dirac operator whose spectrum on a finite approximation to S³ converges to Camporesi-Higuchi?

### Domain 5: Relativistic quantum chemistry qubit Hamiltonians

For Tier 2/3 downstream, but useful to know now:
- Has anyone mapped Dirac-Hartree-Fock or X2C Hamiltonians to qubit form for VQE/QPE?
- Known resource estimates (Pauli count, 1-norm) for relativistic Gaussian bases at heavy atoms?
- This defines the market we'd be entering — we don't need to reproduce it, just know what exists so GeoVac's positioning is accurate.

## Output format

```
# Explorer Report — Dirac-on-S³ Infrastructure Scan

## Executive summary
[3–5 bullets: what's off-the-shelf, what needs to be built, biggest surprise]

## Domain 1: Dirac on S³ (free)
### Finding 1.1 — [short title]
Citation: [full ref]
What's usable: [paragraph]
GeoVac compatibility: [paragraph — (n,l,m) labeling? rational arithmetic? fiber-bundle compatible?]

[etc., organized by domain]

## Gaps
[list of objects Tier 1 will have to build from scratch because no usable literature exists]

## Recommendations for Tier 1 scoping
[2–4 concrete suggestions based on what you found, addressed to the Leader agent]
```

Target 2000–3000 words. Prioritize finds that directly enable the α theorem attempt over broad surveys. If a domain is too thin to be useful, say so — a short "negative" section is more valuable than padded coverage.

## Constraints

- **Do not propose code or sprint tracks.** Those are the Leader's and PM's job.
- **Do not propose paper drafts.** That's the Reviewer's job downstream.
- **Do cite specifically** — "Camporesi-Higuchi 1996" is usable; "the literature suggests" is not.
- **Flag anything adjacent to GeoVac's Section 3 failed approaches** — if a literature technique smells like "single-center molecular" or "shared-exponent unified basis," note it so we don't re-derive dead ends in Dirac clothing.
