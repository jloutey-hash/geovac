# Sprint 2: Dirac Cusp + SS/SOO Breit Interaction

**Version target:** v2.13.0 (with Sprint 1 follow-ups) → v2.14.0 on Sprint 2 completion
**Date:** April 2026
**Tracks:** Two parallel tracks (independent, no blocking dependencies)

Sprint 2 has two goals, both downstream of the Dirac-on-S³ Tier 1-3 infrastructure:

1. **Track DC (Dirac Cusp):** Does the Dirac-Coulomb coalescence condition converge faster in the angular basis than Kato's non-relativistic cusp? Tests whether the relativistic wavefunction has a "naturally softer" cusp on S³.

2. **Track BR (Breit Rank-2):** Does the Breit interaction (magnetic spin-spin + spin-other-orbit, rank-2 tensor) preserve the Gaunt selection rules that underlie Paper 22's angular sparsity theorem? Addresses the 66-211% fine-structure gap from T8 honest negative.

---

## Track DC: Dirac-Coulomb Cusp Characterization

**Goal:** Determine whether the Dirac-Coulomb two-electron coalescence is structurally different from the Schrödinger-Coulomb (Kato) cusp in a way that changes angular basis convergence on S³. This is the physics question behind the memory file's claim that "full cusp regularization requires one-loop QED."

**Principle:** Transcendental Cataloging + Algebraic Deconstruction

### Background

The non-relativistic Kato cusp condition at r₁₂ = 0 is:
```
    (∂ψ/∂r₁₂)|_{r₁₂=0} = (1/2) ψ(r₁₂=0)
```

This is a first-derivative condition on the many-body wavefunction. The Schwartz partial-wave expansion of the cusp converges as l⁻⁴, which is the bottleneck for graph-native CI at Z≥4.

The Dirac-Coulomb equation has a different coalescence structure:
- Dirac wavefunctions are 4-component spinors, not scalars
- The singular part of the electron-electron interaction includes not just 1/r₁₂ but also Breit corrections (α² magnetic)
- The coalescence condition mixes large and small components via the γ matrices

Paper 18 §II.B classifies the 1/r₁₂ cusp as an "embedding exchange constant" with value "removable by similarity transformation" (TC Jastrow). But CUSP-3 proved TC is dead in the composed basis. The open question: does the Dirac structure give a natural regularization that TC doesn't?

### Sub-tracks

**DC-A: Derive the Dirac-Coulomb coalescence condition (algebraic)**

1. Start from the two-electron Dirac-Coulomb Hamiltonian:
   ```
   H = h_D(1) + h_D(2) + 1/r₁₂
   ```
   where h_D is the single-particle Dirac operator.

2. Compute the coalescence condition by matching the wavefunction near r₁₂ = 0. For pure 1/r₁₂ (no Breit), the relativistic analog of Kato is:
   - Known in the literature: Kutzelnigg 1988, Salomonson & Öster 1989. 
   - The relativistic cusp condition depends on the spin state (singlet vs triplet) and on γ = √(1-(Zα)²).
   - Singlet: ψ'/ψ|_{r=0} = 1/(2γ) for 1s²
   - Triplet: different coefficient

3. Compute the partial-wave expansion of the Dirac cusp. Does the leading convergence rate differ from Schwartz l⁻⁴?

4. All of this should be pure algebra — Dirac matrix elements, spinor algebra, and partial-wave expansion. Use the Tier 1-3 infrastructure (`dirac_s3.py`, `dirac_matrix_elements.py`, spinor harmonics).

**DC-B: Numerical test — He 1s² cusp on the composed Dirac graph**

1. Build a two-electron Dirac-Coulomb Hamiltonian for He in the composed qubit encoding (leverage Tier 2 T3 machinery in `composed_qubit_relativistic.py`).

2. At n_max=2, 3, 4: compare the angular-basis convergence of:
   - Schrödinger-Coulomb FCI (existing, via `casimir_ci.py` graph-native CI)
   - Dirac-Coulomb FCI (new, using spinor basis)

3. Fit the error vs l_max. If Dirac gives a faster convergence rate (l⁻⁶ or l⁻⁸?), that's the positive result. If same rate (l⁻⁴), the Dirac structure doesn't help — the cusp is an angular-basis-completeness issue, not a spin-structure issue.

4. For Z=4 (Be²⁺) where the graph validity boundary is suppressed: the comparison is cleanest.

**DC-C: Structural interpretation (synthesis)**

1. If DC-B positive: document the mechanism. Is it the γ = √(1-(Zα)²) correction changing the short-distance behavior? The spinor coupling mixing singlet/triplet? The Breit magnetic term?

2. If DC-B negative: document why the Dirac cusp is "same as Schrödinger cusp" for this purpose. This is the clean result — the cusp is genuinely an angular basis issue, not a relativistic one. Paper 18 updated accordingly.

3. Update Paper 18 §II.B: whether 1/r₁₂ as "embedding exchange constant" is Schrödinger-specific or universal across 2nd-order operator sectors.

### Success criteria

- Dirac-Coulomb coalescence condition derived symbolically (sympy exact)
- Partial-wave convergence rate for Dirac cusp (analytic + numerical)
- He n_max=2-4 FCI with Dirac-Coulomb (convergence table)
- Clean verdict: Dirac cusp faster, same, or slower than Kato

### Files to read first
- `geovac/dirac_s3.py` — Dirac spinor infrastructure
- `geovac/dirac_matrix_elements.py` — Szmytkowski matrix elements, κ↔(l,σ) bridge
- `geovac/composed_qubit_relativistic.py` — Tier 2 T3 relativistic composed builder
- `geovac/casimir_ci.py` — graph-native CI with Slater integrals
- `geovac/hypergeometric_slater.py` — exact rational R^k integrals
- CLAUDE.md §2 cusp characterization bullets (v2.9.0-v2.9.2)
- Paper 18 §II.B (1/r₁₂ as embedding exchange constant)

### Failed approaches to avoid
- TC (transcorrelated) — decisively dead at every n_max (CUSP-3). Don't attempt.
- TC angular gradient — 100×+ cost/benefit worse (Track BX-4).
- Alpha-only cusp factor on S⁵ — negative (Track U).
- Graph absorption of 1/r₁₂ into S⁵ — structural obstruction (Track W).

---

## Track BR: Breit Interaction and Rank-2 Gaunt Selection Rules

**Goal:** Determine whether the Breit interaction (spin-spin dipole-dipole + spin-other-orbit) preserves the Gaunt selection rules that give Paper 22's angular sparsity theorem. Fine-structure accuracy (currently 66-211% error on He/Li/Be 2p splittings) is the downstream test — but the structural question about selection rules is the primary science.

**Principle:** Natural Geometry Search + Transcendental Cataloging

### Background

T8 established that Darwin + mass-velocity don't improve multi-electron 2p splittings because both 2p states share l=1 (Darwin=0 for l≥1; MV identical for same l). The residual 66-211% error comes from the two-body Breit interaction, which has two parts:

1. **Spin-spin (SS):** magnetic dipole-dipole, ∝ α² [σ₁·σ₂/r₁₂³ − 3(σ₁·r̂)(σ₂·r̂)/r₁₂³]
2. **Spin-other-orbit (SOO):** current-current, ∝ α² L₁·σ₂/r₁₂³ type terms

Both have **rank-2 tensor structure** in spin space (vs Coulomb's rank-0 and dipole's rank-1). The question: does this preserve the Gaunt angular selection rules that make Paper 22's theorem work?

Paper 22's angular sparsity theorem currently covers:
- Rank-0 (Coulomb): 1.44% density at l_max=3
- Rank-1 (dipole): verified in earlier work

Rank-2 (Breit) is the next test. If selection rules hold, Paper 22 extends cleanly. If they don't, we learn the boundary of angular sparsity.

### Sub-tracks

**BR-A: Angular decomposition of Breit operator (pure algebra)**

1. Write the SS and SOO operators in tensor form:
   ```
   H_SS = (α²/r₁₂³) Σ_q T²_q(σ₁) · T²_q(σ₂) · Y²_{-q}(r̂₁₂)
   H_SOO = similar structure with L₁·σ₂
   ```

2. Compute the angular matrix elements in the (κ, m_j) basis using the existing Szmytkowski machinery (`dirac_matrix_elements.py`).

3. Verify: are the only nonzero matrix elements those satisfying |Δl| ≤ 2, |Δj| ≤ 2, |Δm_j| ≤ 2? This is the Gaunt rank-2 selection rule. Use sympy Wigner 3j/6j symbols.

4. Compute the angular sparsity density d_Breit(l_max) for l_max = 0..5. Compare with d_scalar and d_spinor from Paper 22.

**BR-B: Radial part of Breit matrix elements**

1. The radial part ⟨1/r₁₂³⟩ is more singular than Coulomb 1/r₁₂. In partial-wave expansion, it has a different structure.

2. Use the Dirac radial functions (from T7, when available) or hydrogenic approximations (from T1).

3. Check: are the Breit radial integrals algebraic (Bethe-Salpeter rational)? Or do they introduce new transcendental content?

4. For the T7-deferred n_r ≥ 1 states with s = -2, -3: the Kramers-Pasternak recursion is needed. Either implement it (extension to T7) or restrict to n_r = 0 states for BR.

**BR-C: He 2³P fine-structure benchmark**

1. Implement two-body Breit (SS + SOO) in the composed relativistic He Hamiltonian.

2. Diagonalize for the 2³P multiplet (3 levels: J=0, 1, 2).

3. Current error: 66% on the 2³P span. Target: <20%.

4. If target hit: Paper 14 §V and Paper 20 Tier 2 table updated with SS/SOO row. Fine-structure benchmark table adds Breit column.

5. If not hit: honest negative — identify what's missing (higher-order QED? multi-electron correlations beyond Breit?).

**BR-D: Paper 22 extension**

1. If BR-A positive (Gaunt rank-2 selection rules hold): add "spinor + Breit" row to Paper 22's sparsity density table.

2. The extended theorem: ERI density depends on l_max and tensor rank, not on V(r). Applies to Coulomb (rank-0), dipole (rank-1), and Breit (rank-2).

3. If BR-A negative: document which selection rule breaks and why.

### Success criteria

- Angular sparsity density d_Breit(l_max) for l_max = 0..5 computed exactly
- Gaunt rank-2 selection rule verified (or counterexample found)
- He 2³P fine-structure benchmark: target <20% error on the multiplet span
- Paper 22 extended (positive) or obstruction documented (negative)

### Files to read first
- `geovac/dirac_matrix_elements.py` — angular matrix elements, Wigner symbols
- `geovac/spin_orbit.py` — T2 Breit-Pauli H_SO implementation pattern
- `geovac/composed_qubit_relativistic.py` — where to add Breit terms
- `papers/core/paper_22_angular_sparsity.tex` — existing theorem statement
- `docs/tier3_verdict.md` T8 section — honest negative and deferred SS/SOO

### Failed approaches to avoid
- Attempting Breit in the composed basis without first verifying Gaunt selection rules — would produce dense matrices with no sparsity advantage
- Using scalar (n,l,m) basis for Breit — must use (κ,m_j) spinor basis

---

## Sprint 2 PM Prompt

```
Read CLAUDE.md, docs/sprint2_tier4_plan.md, debug/q1_reframing_memo.md
(for Sprint 1 context), and the files listed in each track.

Dispatch Track DC and Track BR as parallel sub-agents. The tracks are
independent: DC tests the cusp convergence question, BR tests the Breit
selection rule question. Both leverage Tier 1-3 infrastructure.

Track DC (Dirac Cusp):
  Sub-agent 1 (DC-A): Derive Dirac-Coulomb coalescence condition
    symbolically. Partial-wave convergence rate analytically.
    Pure algebra via sympy. Reference Kutzelnigg 1988 / Salomonson
    & Öster 1989.

  Sub-agent 2 (DC-B): Numerical test. Build Dirac-Coulomb FCI for
    He at Z=4 (Be²⁺ — avoid Z=2 graph validity boundary) at
    n_max=2,3,4. Compare convergence rate to Schrödinger-Coulomb
    FCI. Fit error vs l_max for both.

  Sub-agent 3 (DC-C): Synthesis. Interpret positive or negative
    result. Update Paper 18 §II.B.

Track BR (Breit):
  Sub-agent 4 (BR-A): Angular decomposition of SS and SOO.
    Compute d_Breit(l_max) for l_max=0..5 using sympy Wigner symbols.
    Verify Gaunt rank-2 selection rules.

  Sub-agent 5 (BR-B): Radial Breit matrix elements. Use hydrogenic
    approximation for diagonal states; flag Kramers-Pasternak as
    optional extension.

  Sub-agent 6 (BR-C): He 2³P fine-structure benchmark. Implement
    Breit in composed_qubit_relativistic.py. Compare to NIST.

  Sub-agent 7 (BR-D): Paper 22 extension (conditional on BR-A positive).

Algebraic-first: before any numerical run over a few minutes, confirm
there's no closed form. Angular matrix elements should be exact sympy
rationals; Breit radial integrals should be Bethe-Salpeter-like.

Exit criteria:
- DC: Dirac cusp convergence rate verdict (faster / same / slower than Kato)
- BR: Gaunt rank-2 verdict + He 2³P <20% error (or honest negative)
- Papers 18, 22, 14, 20 updated as appropriate
```

---

## Sprint 3 Dependency

Sprint 3 (heavy atoms, [Kr]/[Xe] frozen cores) is independent of both DC and BR. Can start anytime. The Sunaga 2025 matched-Q comparison is the market test for the relativistic pipeline, with BR's fine-structure accuracy (if achieved) as a bonus.

---

## What This Sprint Resolves

1. **Cusp connection to Dirac/QED:** Either (a) the Dirac cusp converges faster than Kato, validating the intuition that relativistic structure helps, or (b) same rate, meaning the cusp is purely an angular basis issue and QED won't help directly. Either answer is valuable.

2. **Paper 22 scope:** Does the angular sparsity theorem extend to rank-2 tensor operators? This determines whether Breit-corrected Hamiltonians preserve GeoVac's structural sparsity advantage.

3. **Fine-structure honest negative:** Close the 66-211% gap or document the next layer below Breit.
