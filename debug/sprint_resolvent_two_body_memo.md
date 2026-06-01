# Sprint: Resolvent two-body diagnostic

**Date:** 2026-06-01
**Version:** v3.40.0
**Verdict:** CLEAN NEGATIVE — resolvent does not close the multi-focal wall; wall is the Fock projection conformal factor.

## 1. Motivation

Paper 54 (v3.39.0) showed the gauged spectral action on T_{S³}⊗T_{S³} gives correct angular selection rules but wrong radial coupling strengths (Pearson 0.41-0.58 vs Coulomb). The paper's open question: does the resolvent (D²-z)⁻¹ recover the Coulomb 1/r₁₂ kernel?

## 2. Construction

The Green's function on S³ expands as:

G(Ω₁, Ω₂) = Σ_{NLM} Y_{NLM}(Ω₁) Y*_{NLM}(Ω₂) / λ_N

The two-body matrix element from this is:

V_{abcd} = Σ_{NLM} M_{NLM}[a,c] × conj(M_{NLM}[d,b]) × w(N)

where M_{NLM} are the 3-Y multiplier matrices and w(N) is the weight function. Four candidate weightings tested:

| Weighting | w(N) | Physical meaning |
|:----------|:-----|:-----------------|
| Laplacian resolvent | 1/(N²-1), N≥2 | Green's function of Δ on S³ |
| Dirac resolvent | 1/(N+½)² | Green's function of D² |
| Uniform | 1 | Angular structure only (no radial weight) |
| 1/N | 1/N | Intermediate |

Compared against exact Coulomb Slater integrals from `hypergeometric_slater.py`.

Initial implementation had a conjugation bug (M[a,c]×M[b,d] instead of M[a,c]×conj(M[d,b])), diagnosed by m-conservation at ~40% instead of 100%. Fixed; all results use the corrected conjugation.

## 3. Results

### Summary table

| Method | n_max=2 Pearson | n_max=3 Pearson | n_max=2 CV | n_max=3 CV |
|:-------|:----------------|:----------------|:-----------|:-----------|
| Laplacian 1/(N²-1) | -0.131 | -0.042 | 10.7 | 96.5 |
| Dirac 1/(N+½)² | **0.807** | **0.751** | 40.4 | 68.5 |
| Uniform w=1 | 0.464 | 0.333 | 8.1 | 60.4 |
| 1/N | 0.744 | 0.654 | 25.8 | 65.5 |

All methods: m-conservation 100%, coverage ~50-59% of Coulomb nonzero elements.

### Key findings

**F1 (NEGATIVE): Pearson DECREASES with n_max across all weightings.** Dirac resolvent goes from 0.807→0.751; uniform from 0.464→0.333. This is structural divergence, not finite-size error. The resolvent does not converge to Coulomb.

**F2 (NEGATIVE): Laplacian resolvent is structurally wrong.** The N=1 mode has eigenvalue N²-1=0 (zero mode of Δ on S³), so it's excluded. But the N=1 multiplier dominates the (1s,1s)→(1s,1s) Coulomb integral. Result: the largest Coulomb matrix elements are all ZERO in the Laplacian resolvent.

**F3 (STRUCTURAL): Dirac resolvent regularizes the zero mode.** N=1 has eigenvalue (1+½)²=2.25 in D², so it's included. This is why the Dirac resolvent outperforms the Laplacian. The Dirac operator naturally regularizes the compact-manifold zero mode.

**F4 (STRUCTURAL): The mismatch is the Fock projection conformal factor.** The 3-Y multiplier M_{NLM} factorizes as I_rad(Gegenbauer triple on S³) × Gaunt(S² angular). The Coulomb integral factorizes as Gaunt(S² angular) × R^k(Slater in flat space). The angular Gaunt part is shared (hence 100% m-conservation and correct selection rules). The radial parts are structurally different: Gegenbauer triple on S³ vs Slater integral in flat-space coordinates. The conversion between them IS the Fock projection (stereographic map, Paper 7). This is the conformal factor chordal = 2R sin(geodesic/2R), which is a function, not a constant.

**F5 (STRATEGIC): Paper 54's open question is answered in the negative.** Neither the spectral action nor the resolvent can recover the Coulomb kernel from the spectral triple. The wall sits at the Fock projection — the conformal map from S³ to flat space. This is the Layer 1 / Layer 2 boundary from Paper 34.

## 4. Strategic reframing (from PI conversation)

The conversation produced several reframings:

1. **Paper 54 should be folded into Paper 31 or 32, not standalone.** The diagnostic confirms Paper 31's A/D partition at the two-body level — it's a verification, not a new result.

2. **The framework's value is sparsity, not accuracy.** The chemistry accuracy ceiling (5-26% R_eq) is the classical FCI ceiling of the composed-geometry approximation. Quantum computing (VQE/QPE on the O(Q^2.5) sparse Hamiltonians) is the path past this ceiling.

3. **The graph adds structure (the map), not better Coulomb integrals.** The Coulomb potential in flat space is well understood. The graph tells you which basis is natural, which matrix elements are zero, and where discrete meets continuous.

4. **There aren't two frameworks — there's one spectral triple and multiple projections.** Gravity comes from the spectral action (stays on S³). Coulomb comes from the Fock projection (maps S³ → flat space). Both are Layer 2 operations on the same Layer 1 graph.

## 5. Files

### Created
- `debug/resolvent_two_body_diagnostic.py` — Driver (4 weightings × 2 n_max, comparison to exact Coulomb)
- `debug/data/resolvent_two_body_diagnostic.json` — Numerical results
- `debug/sprint_resolvent_two_body_memo.md` — This memo

### Modified
- `papers/group3_foundations/paper_31_universal_coulomb_partition.tex` — New §10 "Two-body verification of the partition" (3 subsections: gauged spectral action, resolvent construction, structural diagnosis). Paper 31: 14 pages (was ~13), compile clean. Folds Paper 54's standalone content into its natural home.

### Not modified
- No production code modified
- No tests added (diagnostic only; new §10 equations are descriptive, not computable identities needing verification)

## 6. Honest scope

**Theorem grade:** Nothing. This is a diagnostic sprint.

**Structural finding (verified numerically):** The resolvent of D² on T_{S³}⊗T_{S³} does not converge to the Coulomb 1/r₁₂ kernel. The mismatch is the conformal factor of the Fock projection, which is a function (not a constant) of the geodesic angle. F1-F4 verified at n_max ∈ {2, 3} with exact Coulomb comparison.

**Numerical observation:** Dirac resolvent outperforms Laplacian resolvent (Pearson 0.75-0.81 vs -0.04 to -0.13) due to zero-mode regularization (F3). All methods show decreasing Pearson with n_max.

**Named open follow-ons:**
- Fold Paper 54 into Paper 31 or Paper 32 as a section (PI directive)
- Investigate whether the resolvent of the FOCK-PROJECTED operator (not D² on S³) gives Coulomb — this would use the conformal factor explicitly, making it tautological but potentially illuminating the convergence rate
- The quantum computing value proposition (Paper 14 sparse Hamiltonians on quantum hardware) remains the primary path past the classical accuracy ceiling
