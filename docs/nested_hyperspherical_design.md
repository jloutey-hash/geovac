# Nested Hyperspherical Architecture: Design Document

**Track:** DF Sprints 1-5 (COMPLETE)
**Date:** 2026-04-07
**Status:** INVESTIGATION COMPLETE — MIXED. Novel sparsity results, compact encoding validated, molecular geometry unsolved. Heterogeneous nested NEGATIVE (Sprint 5).

---

## 1. Executive Summary

This document maps the Aquilanti group's nested (recursive) hyperspherical
decomposition onto GeoVac's existing infrastructure and evaluates coupling
scheme choices for preserving Gaunt selection rule sparsity.

**Key finding:** The H-set (pairwise) coupling scheme is recommended for LiH
because it preserves pair-angular-momentum as a good quantum number at each
nesting level, which is exactly the quantum number that Gaunt integrals
constrain. The K-set (sequential) coupling mixes pair angular momenta across
nesting levels, weakening the selection rules.

**Relationship to GeoVac hierarchy:** The nested approach is a strict
generalization of the existing Level 3 → Level 4 progression. Level 3 (He on
S⁵) is the N=2 case of the hierarchy. Level 4 (H₂ on S⁵ in molecule-frame) is
Level 3 with nuclear coordinates added. The nested approach extends this to
N=3,4,... by recursively introducing pair hyperangles that parametrically
couple subsystems, maintaining the full S^(3N-1) topology throughout.

---

## 2. The Aquilanti Nested Hierarchy

### 2.1 Core Idea

For N electrons in 3D, the configuration space is R^(3N). Separating the
hyperradius R = sqrt(Σ r_i²) leaves (3N-1) angular coordinates on S^(3N-1).
The Aquilanti approach organizes these angular coordinates into a binary tree
of nested pair hyperangles.

**Standard (flat) approach:** Treat all (3N-1) angles simultaneously. This is
what the existing Level 4N solver (`n_electron_solver.py`) does for N=4 on
S¹¹ — a 3D finite difference grid in (α₁, α₂, α₃) plus channel labels
(l₁, l₂, l₃, l₄). Full SO(12) problem, ~13,000 channels at l_max=2.

**Nested (tree) approach:** Organize the hyperangles hierarchically. First
couple electrons into pairs (defining pair hyperangles), then couple pairs
together (defining inter-pair hyperangles). At each level of the tree, the
problem is a 1D hyperangular eigenvalue problem coupled to the levels above
and below.

### 2.2 Coordinate Definitions

For 4 electrons (LiH, Be), the Jacobi tree has two levels:

```
Level 2 (inter-pair):     α₂ = arctan(ρ₃₄/ρ₁₂)
                         /                    \
Level 1 (intra-pair): α₁ = arctan(r₂/r₁)   α₃ = arctan(r₄/r₃)
                     /    \                  /    \
Electrons:          e₁    e₂              e₃    e₄
```

where ρ₁₂ = sqrt(r₁² + r₂²) and ρ₃₄ = sqrt(r₃² + r₄²).

The fractional radii are:
- s₁ = cos(α₂)cos(α₁), s₂ = cos(α₂)sin(α₁)
- s₃ = sin(α₂)cos(α₃), s₄ = sin(α₂)sin(α₃)

This is *exactly* the coordinate system already used in `n_electron_solver.py`
(lines 7-18). The difference is not in the coordinates but in the
*angular basis* — how the angular eigenfunctions are expanded.

### 2.3 Two Coupling Schemes

**H-set (pairwise, "democratic"):**
Each pair of electrons (1,2) and (3,4) first has its own pair angular momentum
l₁₂ and l₃₄ coupled from individual l₁,l₂ and l₃,l₄ via standard CG
coefficients. Then the pair angular momenta l₁₂, l₃₄ are coupled to give total
L. The quantum numbers at each level are:

```
Individual:    l₁, l₂, l₃, l₄         (per-electron angular momentum)
Pair-coupled:  l₁₂ ∈ |l₁-l₂|..l₁+l₂   l₃₄ ∈ |l₃-l₄|..l₃+l₄
Total:         L ∈ |l₁₂-l₃₄|..l₁₂+l₃₄
```

The hyperspherical harmonic basis element is:
Y_{ν; l₁l₂l₁₂, l₃l₄l₃₄, L}(Ω) = 
    [Y_{l₁}(Ω₁) ⊗ Y_{l₂}(Ω₂)]_{l₁₂} ⊗ [Y_{l₃}(Ω₃) ⊗ Y_{l₄}(Ω₄)]_{l₃₄}]_L
    × Φ_{ν}^{l₁₂,l₃₄}(α₁, α₂, α₃)

where Φ is the hyperangular part (Jacobi polynomials in the nested scheme).

**K-set (sequential, "hierarchical"):**
Electrons are coupled sequentially: (1,2) → 12, then (12,3) → 123, then
(123,4) → 1234. Quantum numbers:

```
Step 1:  l₁₂ ∈ |l₁-l₂|..l₁+l₂
Step 2:  l₁₂₃ ∈ |l₁₂-l₃|..l₁₂+l₃
Step 3:  L ∈ |l₁₂₃-l₄|..l₁₂₃+l₄
```

### 2.4 Relationship to Existing GeoVac Levels

| GeoVac Level | N | Angular Space | Aquilanti Equivalent |
|:-------------|:-:|:--------------|:---------------------|
| Level 1 (H)  | 1 | S² (SO(3))    | Trivial: single Y_l |
| Level 3 (He) | 2 | S⁵ (SO(6))    | H-set = K-set (only one pair, identical) |
| Level 4 (H₂) | 2+nuc | S⁵ + R      | Same as Level 3 with nuclear parametric coupling |
| Level 4N (LiH) | 4 | S¹¹ (SO(12)) | H-set or K-set binary tree |
| Level 5 (composed) | 4 | S⁵ × S⁵ + PK | Factorized approximation to nested |
| **Nested (proposed)** | 4 | S¹¹ (SO(12)) | Full tree with Gaunt-truncated basis |

**Critical insight:** The composed architecture (Level 5) is an *approximation*
to the nested hierarchy where:
1. The inter-pair hyperangle α₂ is frozen (pairs are independent)
2. Inter-pair coupling is replaced by PK pseudopotential
3. Each pair gets its own Z_eff (different basis functions)

The nested approach unfreezes α₂ and restores inter-pair coupling directly
through V_ee and V_ne matrix elements.

---

## 3. Where Gaunt Integrals Enter

### 3.1 Electron-Electron Repulsion

For a pair (i,j), the multipole expansion gives:
1/r_{ij} = Σ_k (r_<^k / r_>^{k+1}) P_k(cos θ_{ij})

In hyperspherical coordinates, r_i = R·s_i(Ω), so:
1/r_{ij} = (1/R) Σ_k [min(s_i,s_j)/max(s_i,s_j)]^k / max(s_i,s_j) · P_k(cos θ_{ij})

The angular matrix element of P_k(cos θ_{ij}) between states with angular
momenta (l_i, l_j) and (l_i', l_j') factorizes into Gaunt integrals:

⟨l_i m_i|P_k|l_i' m_i'⟩ = G(l_i, k, l_i') × δ_{m_i, m_i'}

The Gaunt selection rule requires:
- |l_i - l_i'| ≤ k ≤ l_i + l_i'  (triangle inequality)
- l_i + k + l_i' even  (parity)

### 3.2 The Key Question: Does Nesting Preserve This?

In the **flat** (Level 4N) approach, V_ee between pair (1,2) produces a
coupling matrix in the 4-channel space (l₁, l₂, l₃, l₄) that is:
- Gaunt-restricted in (l₁, l₂) → (l₁', l₂')
- Diagonal in (l₃, l₄) [electrons 3,4 not involved]

This gives block-diagonal structure. The question is whether this structure
is preserved when we change basis from the uncoupled (l₁,l₂,l₃,l₄) to the
coupled (l₁₂, l₃₄, L) representation.

**Answer: YES, for the H-set coupling.** The Gaunt integral selection rules
on (l₁, l₂) → (l₁', l₂') translate to selection rules on the coupled pair
angular momentum l₁₂:

For V_ee between electrons 1 and 2:
- l₃, l₄ are unchanged (diagonal)
- Therefore l₃₄ is unchanged (diagonal in pair angular momentum)
- l₁₂ changes by |Δl₁₂| ≤ 2k (from recoupling CG coefficients)
- But k is bounded by the Gaunt triangle inequality

The critical point: in the H-set scheme, l₃₄ is a *good quantum number*
for V_ee(1,2). This means the coupling matrix is block-diagonal in l₃₄.
Similarly, V_ee(3,4) is block-diagonal in l₁₂. The inter-pair V_ee terms
(1,3), (1,4), (2,3), (2,4) couple both l₁₂ and l₃₄, but the Gaunt
restrictions on each individual electron's angular momentum still constrain
the transitions.

**For the K-set coupling:** l₁₂₃ is not a pair angular momentum — it's a
3-electron cumulative coupling. V_ee(3,4) does not leave l₁₂₃ diagonal
(electron 3 is mixed into l₁₂₃). This destroys the block-diagonal structure
for cross-pair interactions.

### 3.3 Selection Rules in the Nested H-Set Basis

For 4 electrons with H-set coupling at L=0, the complete selection rules are:

**Intra-pair V_ee (pairs (1,2) or (3,4)):**
- Same pair angular momentum: Gaunt restricted (|Δl_i| ≤ k, parity)
- Other pair: diagonal in l (and hence l_paired)
- Hyperangular part: new integrals in α₁ or α₃ (1D, computable)

**Inter-pair V_ee (pairs (1,3), (1,4), (2,3), (2,4)):**
- Both pair angular momenta change: l₁₂ → l₁₂', l₃₄ → l₃₄'
- Individual electron Gaunt restrictions: |l_i - l_i'| ≤ k, |l_j - l_j'| ≤ k
- Hyperangular part: 3D integral in (α₁, α₂, α₃) — the expensive part

**Nuclear V_ne (each electron):**
- Only the electron's own l changes: |Δl_i| ≤ L_multipole
- Paired partner's l unchanged: l_partner unchanged
- For single-center (Be): V_ne is diagonal in all l (monopole only)
- For two-center (LiH): V_ne from far nucleus requires multipole expansion
  (already implemented in `shibuya_wulfman.py`)

### 3.4 ERI Density Estimate (Analytical)

At l_max=1 (l ∈ {0,1}) for 4 electrons, the H-set coupled basis at L=0:

Possible pair angular momenta:
- l₁₂ = 0 (from l₁=l₂=0, or l₁=l₂=1 coupled to 0)
- l₁₂ = 1 (from l₁=0,l₂=1 or l₁=1,l₂=0, or l₁=l₂=1 coupled to 1)
- l₁₂ = 2 (from l₁=l₂=1 coupled to 2)

Similarly for l₃₄. L=0 requires l₁₂ = l₃₄.

Counting the coupled basis states (l₁, l₂, l₁₂ = l₃₄, l₃, l₄) at L=0:

| l₁₂ = l₃₄ | (l₁,l₂) pairs | (l₃,l₄) pairs | States |
|:-----------|:--------------|:--------------|:-------|
| 0 | (0,0), (1,1) = 2 | same = 2 | 4 |
| 1 | (0,1), (1,0), (1,1) = 3 | same = 3 | 9 |
| 2 | (1,1) = 1 | same = 1 | 1 |
| **Total** | | | **14** |

Compare to uncoupled basis at M=0: (l_max+1)⁴ = 16 channels (but with
S₄ [2,2] projection, ~3 independent). The coupled basis has 14 states before
S₄ projection, reducing to a similar count.

**Intra-pair V_ee sparsity:** For V_ee(1,2), l₃₄ is conserved. Of the 14
states, coupling is block-diagonal with blocks of size 2 (l₃₄=0), 3 (l₃₄=1),
1 (l₃₄=2). Within each block, the Gaunt selection rules restrict l₁₂
transitions. ERI density from intra-pair: comparable to Level 3 He (~8.9%).

**Inter-pair V_ee sparsity:** For V_ee(1,3), both l₁₂ and l₃₄ can change,
but each individual electron's Gaunt constraint limits allowed transitions.
The 6j/9j recoupling coefficients introduce additional selection rules beyond
the bare Gaunt restriction. This is where the nested approach may achieve
*better* sparsity than the flat approach — the recoupling coefficients can
vanish even when individual Gaunt integrals don't.

**Estimated ERI density:** 5-10% at l_max=1, potentially *lower* than the
flat uncoupled basis (8.9% for s/p blocks) due to additional recoupling
selection rules. This needs to be verified numerically in Sprint 2.

---

## 4. Antisymmetry in the Nested Framework

### 4.1 The Composed Architecture's Master Obstruction

CLAUDE.md Section 3 documents three failed attempts at inter-group
antisymmetry in the composed architecture. The root cause: composed geometry
factorizes R^(3N) into independent coordinate systems (S⁵_core × S⁵_valence),
and the permutation operator P_{ij} that swaps a core electron with a valence
electron maps *between* these coordinate systems. Since the coordinates are
incompatible (different Z_eff, different hyperradii), the matrix element
⟨core state|P_{ij}|valence state⟩ cannot be evaluated within either
coordinate system alone.

### 4.2 Why Nested Avoids This Obstruction

The nested approach works in a *single* S^(3N-1) coordinate system. All N
electrons share the same hyperradius R and the same angular manifold. The
permutation operator P_{ij} acts on this single space by permuting electron
labels:

P_{12}: (α₁, α₂, α₃, Ω₁, Ω₂, Ω₃, Ω₄) → (π/2-α₁, α₂, α₃, Ω₂, Ω₁, Ω₃, Ω₄)

For within-pair exchange (P_{12} or P_{34}), this is a simple reflection in
the pair hyperangle. For cross-pair exchange (P_{13}, etc.), it involves a
transformation mixing all three hyperangles — more complex but still
well-defined within the single coordinate system.

The S_N irrep projection machinery from `n_electron_solver.py` (lines 119-239)
applies directly: S₄ acts on the angular channel labels (l₁,l₂,l₃,l₄) by
permutation, and the [2,2] projector selects the correct singlet subspace.
In the coupled basis, S₄ acts via 6j recoupling:

P_{13}: |l₁l₂l₁₂; l₃l₄l₃₄; L⟩ → Σ_{l₁₃,l₂₄} W(l₁l₂l₃l₄; l₁₂l₁₃...) |...⟩

where W involves Racah (6j) coefficients. This is computable, not
obstructed.

### 4.3 Cost of Antisymmetry

The S₄ [2,2] projection reduces the basis by a factor of ~1/6 at l_max=1
(from 16 uncoupled channels to ~3 independent [2,2] channels). In the
coupled basis, the reduction factor is similar but the surviving states
are different linear combinations. The key advantage: antisymmetry is
*exact* in the nested scheme, whereas it is *approximated* (via PK) in
the composed scheme.

---

## 5. Relationship to Balanced Coupled Approach (Paper 19)

The balanced coupled approach (Track CD) is an intermediate between composed
and nested:

| Feature | Composed (Level 5) | Balanced Coupled | Nested (Proposed) |
|:--------|:-------------------|:-----------------|:------------------|
| Coordinate system | S⁵ × S⁵ (separate) | S⁵ × S⁵ (separate) | S¹¹ (single) |
| Inter-pair α₂ | Frozen | Frozen | Dynamic variable |
| Cross-center V_ne | Via PK or SW multipole | SW multipole (exact) | SW multipole (exact) |
| Inter-pair V_ee | None (PK proxy) | Cross-block ERIs | Direct (all 6 pairs) |
| Antisymmetry | PK approximation | Occupation number | Exact S₄ projection |
| Pauli terms (LiH) | 334 | 878 | TBD (Sprint 3) |
| FCI energy error | ~5% (PK) | 1.8% (n_max=2) | TBD |

The balanced coupled approach already showed that PK-free accuracy is
achievable (0.20% at n_max=3). The nested approach would maintain this
accuracy while potentially achieving better sparsity through the hierarchical
angular truncation.

**Critical question for Sprint 2:** Is the nested Pauli count closer to 334
(composed, block-diagonal but approximate) or 3,288 (full 4N, exact but dense)?
The answer depends entirely on whether the Gaunt selection rules in the
coupled basis produce enough sparsity to offset the larger Hilbert space.

---

## 6. Infrastructure Reuse

### 6.1 Existing Code That Transfers Directly

| Module | What Transfers | What Changes |
|:-------|:---------------|:-------------|
| `n_electron_solver.py` | Coordinates (α₁,α₂,α₃), S₄ permutation matrices, [2,2] character table | Angular basis: product Y_lm → coupled Y_{l₁₂,l₃₄,L} |
| `n_electron_scope.py` | Casimir eigenvalues, dimension counting | Channel enumeration in coupled scheme |
| `hyperspherical.py` | Level 3 angular solver (template for pair sub-problems) | Used as pair-level building block |
| `shibuya_wulfman.py` | Cross-center V_ne multipole + Wigner D rotation | Operates on angular matrix elements (basis-agnostic) |
| `composed_qubit.py` | JW encoding, Pauli counting | Feed different 1-body/2-body integrals |
| `ecosystem_export.py` | Export API | Add `method='nested'` option |

### 6.2 New Infrastructure Required

| Component | Description | Estimated Complexity |
|:----------|:------------|:---------------------|
| Coupled angular basis | H-set CG coupling of Y_lm products | Medium (Wigner 3j already available) |
| Recoupling coefficients | 6j/9j symbols for cross-pair V_ee in coupled basis | Medium (can use sympy.physics.wigner) |
| Nested hyperangular integrals | 1D integrals in α₁, α₂, α₃ for Jacobi polynomial basis | Low (analogous to Level 3 Gegenbauer integrals) |
| Cross-pair V_ee matrix elements | Full 3D hyperangular integrals for inter-pair repulsion | High (most complex new component) |

---

## 7. Coupling Scheme Recommendation

### 7.1 H-Set Advantages

1. **Pair angular momentum is conserved by intra-pair V_ee.** This creates
   block-diagonal structure in the coupling matrix — the same structural
   sparsity that makes the composed architecture sparse.

2. **Natural match to GeoVac's existing pair decomposition.** The composed
   architecture already groups electrons into pairs (core pair, bond pairs,
   lone pairs). The H-set coupling formalizes this pairing within a single
   coordinate system.

3. **6j coefficients for inter-pair recoupling are well-known.** The
   transformation between H-set and K-set is a Racah W-coefficient (6j
   symbol), which is algebraic and has its own triangle inequality selection
   rules that can only *reduce* the number of nonzero matrix elements.

4. **Compatible with S₄ [2,2] projection.** The [2,2] irrep of S₄ for 4
   electrons in the H-set basis has a clean decomposition: P₁₂ (within-pair
   exchange) acts on α₁ by reflection, and P₁₃ (cross-pair exchange) acts
   via 6j recoupling. Both are computable.

### 7.2 K-Set Disadvantages

1. **No conserved pair angular momentum for cross-pair V_ee.** When electron 3
   is coupled to the (1,2) subsystem first, V_ee(3,4) does not leave l₁₂₃
   diagonal — it changes l₃ and hence l₁₂₃.

2. **Asymmetric treatment of electrons.** The K-set privileges the order of
   coupling, making the permutation symmetry less transparent.

3. **Worse for core-valence systems.** In LiH, the natural pairing is
   core(1,2) + valence(3,4). The H-set respects this; the K-set mixes core
   and valence electrons in intermediate coupling steps.

### 7.3 Decision

**Use H-set (pairwise) coupling for LiH.** The pair angular momentum
conservation under intra-pair V_ee directly produces the block-diagonal
ERI structure that Gaunt selection rules enforce. This is the coupling
scheme most likely to preserve O(Q^2.5) or better Pauli scaling.

---

## 8. Sprint 2 Plan: Analytical Sparsity Verification

With the H-set coupling scheme selected, Sprint 2 must answer:

1. **Count nonzero ERI matrix elements** in the coupled (l₁₂, l₃₄, L=0)
   basis at l_max=1 and l_max=2 for 4 electrons. Separate intra-pair
   (block-diagonal, expected sparse) from inter-pair (coupled, need to
   evaluate 6j restrictions).

2. **Compute ERI density** = (nonzero ERIs) / (total possible ERIs) and
   compare to:
   - Composed s/p block density: 8.9%
   - d-block density: 4.0%
   - Full Level 4N uncoupled density: TBD

3. **Prove or disprove block-diagonality** of the full ERI tensor in the
   coupled basis. Specifically: does V_ee(1,3) in the H-set basis have
   additional vanishing matrix elements (from 6j selection rules) beyond
   the individual Gaunt restrictions?

4. **Go/no-go gate:** ERI density > 15% → STOP and archive.

---

## 9. Risk Assessment Update

After literature synthesis, the risk levels from the investigation prompt
are updated:

| Risk | Original | Updated | Rationale |
|:-----|:---------|:--------|:----------|
| Gaunt sparsity destroyed by nesting | Medium | **Low** | H-set coupling preserves pair angular momentum conservation; intra-pair V_ee is block-diagonal by construction |
| Angular basis too large at ν_max=2 | High | **Medium** | S₄ [2,2] projection reduces by ~1/6; coupled basis further compresses via L=0 constraint |
| SW integrals incompatible with nested basis | Medium | **Low** | SW operates on angular matrix elements ⟨l,m|P_L|l',m'⟩ — these are basis-agnostic |
| Antisymmetry as hard as composed | Medium | **Very Low** | Nested works in single S¹¹; S₄ acts by label permutation + 6j recoupling; no coordinate boundary crossing |
| 1-norm inflation | Medium | **Medium** | Unchanged; must measure empirically |
| **NEW: Inter-pair V_ee integrals expensive** | — | **High** | 3D hyperangular integrals for cross-pair V_ee in coupled basis may require significant new infrastructure |

---

## Sprint 4B Result: Two-Center Charge-Center Origin — NEGATIVE

**Date:** 2026-04-07
**Hypothesis:** Centering hydrogenic orbitals at the charge center (between Li
and H) with V_ne from both nuclei via Shibuya-Wulfman multipole expansion would
fix the R_eq problem (Sprint 4: R_eq = 2.0 bohr, 33.7% error) by distributing
orbital density more symmetrically between the two nuclei.

**Result: NEGATIVE.** The charge-center approach catastrophically degrades the
core electron description. Energy error increases from 11.4% (single-center on
Li) to 48-52% (charge-center) regardless of Z_orb choice.

### Data

**Z_orb scan** (charge-center origin, max_n=2, R=3.015 bohr):

| Z_orb | E (Ha) | Error (%) | Pauli | 1-norm_ni (Ha) |
|:------|:-------|:----------|:------|:---------------|
| 1.0 | -3.891 | 51.8 | 120 | 4.43 |
| 2.0 | -3.933 | 51.3 | 120 | 9.00 |
| 3.0 | -4.069 | 49.6 | 120 | 21.43 |
| 4.0 | -4.122 | 48.9 | 120 | 37.31 |

For comparison: Sprint 4 single-center (Z=3 on Li): E=-7.153, error=11.4%.

**Center position sensitivity** (Z_orb=3, max_n=2, R=3.015):

| Displacement from Li | d_A (bohr) | Error (%) |
|:---------------------|:-----------|:----------|
| 0% (on Li) | 0.000 | 11.4 |
| 2% | 0.060 | 14.4 |
| 5% | 0.151 | 25.3 |
| 10% | 0.302 | 40.2 |
| 25% (charge center) | 0.754 | 49.6 |

### Root Cause

The 1s² core electrons on Li have characteristic orbital radius ~a₀/Z = 0.33
bohr. Moving the orbital center even 0.06 bohr (2% of R) from Li makes the 1s
hydrogenic orbital unable to capture the nuclear cusp at Li. At max_n=2 (5
spatial orbitals), there are not enough higher-n orbitals to compensate.

**Why Paper 15's midpoint centering works but this doesn't:** Paper 15 uses a
continuous hyperspherical angular basis (complete in the angular sector at any
ρ-value), which can represent electron density at any position relative to the
origin. The second-quantized approach uses a truncated hydrogenic orbital basis
that is inherently localized around its center. Midpoint centering in
hyperspherical coordinates ≠ midpoint centering in orbital space.

### All Sprint 4B gates: FAIL

| Gate | Target | Actual | Status |
|:-----|:-------|:-------|:-------|
| R_eq error | < 15% | N/A (energy too wrong for PES) | FAIL |
| FCI error | — | 48.2% (best Z_orb=3.5) | FAIL |
| Pauli terms | ≤ 150 | 120 (PASS trivially) | PASS |
| max_n=3 FCI | feasible | OOM at Q=28 | FAIL |

### Implications

1. The nested architecture remains viable for single-center atoms (Be) and
   single-center molecular approximations (Sprint 4 on Li)
2. For molecular PES with correct R_eq, the nested approach is limited to
   fixed-geometry calculations at externally-specified R
3. The composed architecture (separate blocks with PK) remains the production
   approach for molecular PES
4. A direct hyperspherical approach (Paper 15-style, not second-quantized)
   would avoid this obstruction, but requires the full Level 4N angular
   machinery rather than the orbital-based qubit encoding

### Key Files

- `debug/track_df_lih_twocenter.py` — Two-center builder and comparison
- `debug/data/track_df_lih_twocenter.json` — Saved data
- `tests/test_nested_hyperspherical.py` — TestTwoCenterLiHConstruction, TestTwoCenterLiHFCI

---

## Appendix A: Quantum Number Dictionary

| Symbol | GeoVac Name | Aquilanti Name | Definition |
|:-------|:------------|:---------------|:-----------|
| l_i | Per-electron angular momentum | Orbital angular momentum | Eigenvalue of l̂_i² |
| l₁₂ | Pair angular momentum | Intermediate angular momentum (H-set) | |l₁-l₂| ≤ l₁₂ ≤ l₁+l₂ |
| L | Total angular momentum | Total orbital angular momentum | |l₁₂-l₃₄| ≤ L ≤ l₁₂+l₃₄ |
| ν | Grand angular momentum | Grand angular momentum | SO(3N) Casimir: μ = ν(ν+3N-2)/2 |
| n_ν | Hyperangular excitation | Hyperangular quantum number | Excitation within ν sector |
| R | Hyperradius | Hyperradius | sqrt(Σ r_i²) |
| α₁, α₃ | Pair hyperangles | Kinematic rotation angles (pair) | arctan(r₂/r₁), arctan(r₄/r₃) |
| α₂ | Inter-pair hyperangle | Kinematic rotation angle (inter-pair) | arctan(ρ₃₄/ρ₁₂) |

## Appendix B: Key Literature Mapping

| Aquilanti Concept | GeoVac Equivalent | Status |
|:------------------|:------------------|:-------|
| Hyperspherical harmonics on S^(3N-1) | SO(3N) Casimir eigenfunctions (Paper 7 Sec V) | Proven, eigenvalues known |
| Alternative coupling schemes (H/K) | Channel labeling in n_electron_solver.py | Currently uncoupled; H-set proposed |
| Hyperquantization algorithm | Not applicable (GeoVac is graph → continuous) | Different starting point, same destination |
| Generalized Sturmian method | Tested and characterized (Paper 8, Track BU) | Negative for qubit encoding (1-norm inflation) |
| Shibuya-Wulfman integrals | `shibuya_wulfman.py` (Paper 19) | Implemented for cross-center V_ne |
| Fock projection to S³ | Paper 7 (18 symbolic proofs) | Foundation of GeoVac |

---

## Sprint 5: Heterogeneous Nested — COMPLETE (NEGATIVE)

**Date:** 2026-04-07
**Result:** Löwdin orthogonalization between different-Z orbital sets destroys sparsity and causes BSSE-like overbinding.

### Architecture

Two orbital sets (core at Z_core=3.0, valence at Z_val), both centered on Li, Löwdin-orthogonalized into a single basis for second quantization. Full cross-Z ERIs computed via general Neumann expansion with Gaunt selection rules.

### Z_val Optimization

| Z_val | E_FCI (Ha) | N_Pauli | S_min | 1-norm_ni (Ha) |
|:-----:|:----------:|:-------:|:-----:|:--------------:|
| 0.50 | -7.531 | 1,711 | 0.093 | 28.3 |
| 0.80 | -7.532 | 1,711 | 0.006 | 42.0 |
| 1.00 | -7.665 | 1,711 | 0.006 | 60.2 |
| 1.20 | -7.720 | 1,711 | 0.011 | 57.5 |
| 1.50 | -7.751 | 1,711 | 0.015 | 60.6 |
| **1.70** | **-7.752** | **1,711** | **0.013** | **64.0** |
| 2.00 | -7.738 | 1,711 | 0.009 | 70.0 |
| 2.50 | -7.690 | 1,711 | 0.002 | 81.0 |

Optimal Z_val = 1.70 (matches composed Z_eff for Li valence).

### FCI PES at Z_val = 1.70

| R (bohr) | E_FCI (Ha) | 1-norm_ni (Ha) |
|:--------:|:----------:|:--------------:|
| 1.50 | -7.585 | 63.0 |
| 2.00 | -7.714 | 63.5 |
| **2.50** | **-7.755** | **63.8** |
| 3.015 | -7.752 | 64.0 |
| 4.00 | -7.691 | 64.3 |
| 5.00 | -7.619 | 64.5 |
| 6.00 | -7.569 | 64.7 |
| 8.00 | -7.513 | 65.0 |

R_eq = 2.50 bohr (17.1% error, exact 3.015). D_e = 0.241 Ha (162% error, exact 0.092).

### Comparison Table

| Metric | Uniform | Hetero | Composed+PK | Balanced |
|:-------|:------:|:------:|:-----------:|:--------:|
| Q | 10 | 20 | 30 | 30 |
| N_Pauli | 120 | 1,711 | 334 | 878 |
| 1-norm_ni (Ha) | 13.49 | 63.95 | 37.23 | 74.1 |
| R_eq error (%) | 33.7 | 17.1 | ~5 | 7.0 |
| D_e error (%) | 4.7 | 162.4 | -- | -- |

### Success Criteria Assessment

| Criterion | Target | Actual | Verdict |
|:----------|:-------|:-------|:--------|
| R_eq error | < 15% | 17.1% | **FAIL** |
| D_e error | < 10% | 162% | **FAIL** |
| Pauli terms | ≤ 150 | 1,711 | **FAIL** |
| 1-norm_ni | < 30 Ha | 64.0 Ha | **FAIL** |
| ERI density | ≤ 12% | < 12% | PASS |

### Root Cause Analysis

1. **Löwdin destroys Gaunt sparsity.** The overlap matrix S between Z_core=3 and Z_val=1.7 orbitals is large (1s-1s overlap ~0.89). The Löwdin transform S^{-1/2} mixes radial indices within each (l,m) sector, creating nonzero ERI elements in the transformed basis that were zero in the original angular-momentum-labeled basis. N_Pauli inflates from ~120 (expected for Q=10 nested) to 1,711 at Q=20 — a ~2.5× inflation factor on top of the Q-scaling.

2. **BSSE-like overbinding.** The Löwdin-mixed basis provides more variational flexibility at short R (where both orbital sets overlap significantly with the nuclear region) than at long R (where the valence set decouples). This R-dependent basis quality causes artificial energy lowering at short range, manifesting as D_e = 0.241 Ha (2.6× exact). This is the same mechanism as the v0.9.x LCAO BSSE (cf. Section 1.6 of investigation record).

3. **Near-singular overlap.** S eigenvalues as low as 0.006 (at Z_val=0.8) indicate near-linear dependence. The 1/sqrt(0.006) = 12.9× amplification factor in the Löwdin transform inflates all transformed integrals, particularly the 1-norm.

### Approaches That Failed (Section 3 entry)

| Approach | Sprint | Failure Mode | Key Lesson |
|:---------|:------:|:-------------|:-----------|
| Heterogeneous nested (per-pair Z_eff + Löwdin) | 5 | 1,711 Pauli (5× composed), D_e 162% error (BSSE-like), R_eq 17.1% (still above target) | Löwdin orthogonalization between different-Z hydrogenic orbital sets destroys Gaunt sparsity and introduces BSSE-like overcounting. The composed architecture's block-diagonal factorization is structurally superior to Löwdin mixing. |

### Files

| File | Description |
|:-----|:------------|
| `tests/test_nested_hyperspherical.py` | 14 Sprint 5 tests (helpers + construction + FCI + comparison) |
| `debug/track_df_heterogeneous.py` | Z_val scan + PES investigation script |
| `debug/data/track_df_heterogeneous.json` | Complete numerical results |

---

## Track DF: Final Assessment

**Investigation complete (6 sprints).** The nested hyperspherical architecture has been thoroughly characterized.

**What worked:**
- 6j recoupling sparsity theorem (Sprint 2) — novel, publishable
- Compact encoding: 3× fewer qubits for atoms (Sprint 3)
- Correct D_e for uniform nested (Sprint 4)

**What failed:**
- Single-center molecular: R_eq 33.7% (Sprint 4)
- Two-center: catastrophic energy degradation (Sprint 4B)
- Heterogeneous: Löwdin destroys sparsity, BSSE-like overbinding (Sprint 5)

**The composed architecture remains the production path.** The nested approach cannot compete on any resource metric (Pauli terms, 1-norm, qubit count) once molecular geometry is required. The root obstruction is that core and valence electrons at different effective charges require either (a) separate Hilbert spaces (composed), which preserves sparsity, or (b) a mixed basis with orthogonalization, which destroys it.

**Total tests:** 50 (15 Sprint 3 + 11 Sprint 4 + 10 Sprint 4B + 14 Sprint 5).
