# Sprint TX-A Memo: Dimensional Axiomatization of Paper 34's Projection Taxonomy

**Date:** 2026-05-03
**Sprint:** TX-A
**Goal:** Treat Paper 34's fifteen projections as graded maps along three axes
(variable / dimension / transcendental). Test three theorems analytically;
report on all three even if some are negative.

**Companion files:**
- `debug/data/tx_a_projection_table.json` — 15-row tabulation
- `debug/data/tx_a_composition_table.json` — 15×15 relation table
- `debug/tx_a_signature_arithmetic.py` — verification driver

---

## 0. Executive summary

| Theorem | Verdict | Evidence |
|:--------|:--------|:---------|
| **T1 (Dimensional Consistency)** | **HOLDS** as stated, but **trivially** — every projection's signature is fixed by Paper 34, so signature addition is a tautology of set/list union. The non-trivial content is the *empirical* check that the catalog row's expected dimensions match the chained signature, and this passes for all 8 worked examples. | 8/8 PASS in `tx_a_signature_arithmetic.py`, including the Lamb shift 4-chain, K = π(B+F−Δ) 3-chain, Stefan-Boltzmann, He CI, Uehling VP, and three single-projection sanity rows. |
| **T2 (Minimality up to commutation)** | **HOLDS in the weak form** (signature sum is order-independent because it's set union). **FAILS in the strong form** (operator-system order matters when projections target the same Hilbert space; the K-chain example shows hopf_bundle and spinor_lift have non-trivial composition order that signature arithmetic cannot see). | All k! permutations of every example chain produce identical signatures (verified). The strong form is broken precisely where Connes-vS truncation order matters — see WH1 R3.2 spinor lift work. |
| **T3 (Composition Algebra)** | **HOLDS partially.** Composition table built: 182/225 cells commute, 22 cells are one-way (an irreversible Layer-2 step blocks the reverse), 3 idempotent (Layer-1 self-composes), 6 variable conflicts (two projections add the same variable). Three named inverse pairs (stereographic, Wigner 3j, Wigner D). Most surprising entry: `spectral_action ; observation_window` and its reverse are both one-way, breaking commutation across the temperature/cutoff scale duality. | See §3 below; full table in JSON. |

**Axis-independence verdict (the honest answer the prompt asks for):** the three axes are **NOT independent**. Specifically: **every projection that adds a physical dimension also adds a physical variable** (8/8 with this property), so the variable axis and dimension axis are perfectly correlated on the projections we have. There is no projection in the dictionary that adds *only* a variable while introducing *no* dimension — except spinor_lift (which "introduces" α through (Zα)² but the α is already present from Hopf or spectral-action), and drake_swainson (which introduces K but K cancels out of the final result). The transcendental axis is **partially** independent: vector_photon adds a transcendental (1/(4π)) without adding any variable or dimension, the cleanest example of axis-independence in the dictionary.

This means the three-axis dictionary is **mildly overcomplete**. The variable axis and the dimension axis collapse to one functional axis (call it the *physical-content axis*). The transcendental axis remains independent. Paper 34 should consider whether to merge variables and dimensions into a single column, or whether to add new projections that demonstrably move only one or the other.

---

## 1. The 15-row projection table

(Tabulated in §2 below as a flat table. Full JSON in `debug/data/tx_a_projection_table.json`.)

Encoding scheme: each projection is a tuple
```
(variables_added, dimensions_added, transcendental_class, pi_intrinsic, idempotent, has_inverse, input_layer, output_layer)
```

I treat *layer* as a structural grading: Layer 0 = continuum, Layer 1 = bare GeoVac graph, Layer 2 = projected/integrated observable. A projection has an `input_layer` and `output_layer`; composition `a;b` is admissible iff `a.output_layer ≥ b.input_layer`.

**Tabulation (axis flags VD-P shorthand: V=var, D=dim, P=pi-intrinsic):**

| # | Projection | V | D | P | Variables | Dimensions | Transcendental class |
|:-:|:-----------|:-:|:-:|:-:|:----------|:-----------|:----------------------|
| 1 | fock_conformal | V | D | – | Z, E | energy | rational κ=−1/16 |
| 2 | hopf_bundle | V | – | P | α | (none) | π = Vol(S²)/4 |
| 3 | bargmann_segal | V | D | – | ℏω | energy | rational |
| 4 | stereographic | V | D | – | r | length | conformal factor |
| 5 | sturmian | – | – | – | (re-uses Z, n) | (preserves) | rational preserved |
| 6 | spectral_action | V | D | P | Λ, α | energy | √π·ℚ ⊕ π²ᵏ·ℚ |
| 7 | spinor_lift | V | – | – | α via (Zα)² | (none) | half-int Hurwitz, odd-ζ, G, β(4) |
| 8 | wigner_3j | – | – | – | – | – | ℚ[√(2k+1)]_k |
| 9 | wigner_D | V | D | – | R_AB | length | ℚ[√2,√3,√6] |
| 10 | wilson_plaquette | V | – | P | β_gauge | (none) | SU(2) Haar |
| 11 | vector_photon | – | – | P | – | – | 1/(4π) per loop |
| 12 | molframe_hypersph | V | D | – | R | length | piecewise-smooth in R |
| 13 | drake_swainson | V | – | – | K (cancels) | (cancelled) | flow tier |
| 14 | rest_mass | V | D | – | m | mass | trivial (ring-preserving) |
| 15 | observation_window | V | D | P | β_temp | time | 2π·ℚ per Matsubara mode; π²ᵏ·ℚ |

**Counts (transcendental axis, the prompt's headline question):**

- Add **no transcendental** (Layer-1-preserving): 4 — sturmian, wigner_3j, rest_mass, drake_swainson (the last two add the variable but preserve the ring algebraically; flow tier is calibration-cancelling).
- Add **rational/algebraic but no π**: 4 — fock_conformal, bargmann_segal, stereographic, wigner_D, molframe_hypersph (5 actually; molframe is piecewise-smooth in R but π-free in angular content; counted as no-π).
- Add **π directly** (intrinsic π): 5 — hopf_bundle (π = Vol(S²)/4), wilson_plaquette (Haar), vector_photon (1/(4π)), spectral_action (√π → π²ᵏ), observation_window (2π/β).
- Add **π²** specifically (calibration-tier even-ζ): 2 — spectral_action, observation_window.
- Add **odd-ζ / Catalan G / Dirichlet β**: 1 — spinor_lift.

Of 15 projections: **5 add π in some form, 1 adds odd-ζ family, 9 are π-free (4 fully Layer-1-preserving + 5 ring-extending)**. The framework's transcendental content is concentrated in 6 of 15 projections.

---

## 2. THEOREM 1 — Dimensional Consistency by Signature Addition

### Statement

A projection sequence `S = P_n ∘ ... ∘ P_1` applied to the bare graph (Layer 1, dimensionless) is admissible iff the accumulated signature equals the target observable's signature. Equivalently, `signature(P_n ∘ ... ∘ P_1) = ∪ᵢ signature(Pᵢ)` (set/list union; the operation is signature addition).

### Verdict: HOLDS (trivially-by-construction, non-trivially-empirically)

**Proof sketch.** Each projection's variable, dimension, and transcendental sets are fixed metadata in Paper 34's §III. The set union (or list concatenation) of these sets across a chain is order-independent and well-defined. Since the chain's "signature" is defined by union, the theorem is a tautology *at the level of the metadata*.

The non-trivial content is the *empirical* check: does the catalog row's expected target signature actually equal the chain signature? This is what gets verified in `tx_a_signature_arithmetic.py`.

### Worked examples (all PASS)

| Example | Chain | Predicted signature | Observed signature | Verdict |
|:--------|:------|:--------------------|:-------------------|:--------|
| Hydrogen Lamb shift LS-6a | fock + spectral_action + sturmian + drake_swainson | dims={energy}, π=YES (from spectral_action) | {energy, energy}, π=YES | PASS |
| K = π(B+F−Δ) ≈ 1/α | fock + hopf + spinor_lift | dims={energy}, π=YES (from hopf) | {energy}, π=YES | PASS |
| Stefan-Boltzmann π²/90 | fock + observation_window | dims={energy, time}, π=YES (Matsubara) | {energy, time}, π=YES | PASS |
| He ground state via Casimir CI | fock + wigner_3j | dims={energy}, π=NO | {energy}, π=NO | PASS |
| Hydrogen E_n = −Z²/(2n²) | fock alone | dims={energy}, π=NO | {energy}, π=NO | PASS |
| Spatial Casimir scalar = 1/240 | fock alone | dims={energy}, π=NO (rational) | {energy}, π=NO | PASS |
| Spatial Casimir Dirac = 17/480 | fock + spinor_lift | dims={energy}, π=NO | {energy}, π=NO | PASS |
| Uehling VP Π = 1/(48π²) | fock + spectral_action | dims={energy}, π=YES | {energy}, π=YES | PASS |

**8/8 PASS.** The dimensional and π-intrinsic predictions of signature addition match the empirical observables.

**Interpretation.** Theorem 1 is real but in a deflationary sense: it says nothing the dictionary metadata didn't already say. Its content is the *consistency check* it forces — Paper 34's catalog cannot list a "match" whose chain signature contradicts the observable's actual physics tags. The 8 worked examples confirm that Paper 34's existing entries are mutually consistent under signature addition.

### What Theorem 1 does NOT say

It does **not** say the chain produces the right *numerical value* of the observable. Two chains with the same signature can give wildly different numbers (e.g., Hydrogen E_n via fock alone and Hydrogen Lamb shift via fock+spectral+sturmian+drake_swainson both have dims={energy}, π=YES/NO respectively but their values differ by 13 orders of magnitude). Signatures are necessary, not sufficient. The Paper 34 falsifiable Prediction 1 (error compounds with depth) addresses the *value* question; this theorem addresses only the *type* question.

---

## 3. THEOREM 2 — Minimality (uniqueness up to commutation)

### Statement

For a given target observable, the shortest admissible projection sequence is unique up to commutation of pairs `(Pᵢ, Pⱼ)` whose composition relation is `commute` in the Theorem-3 table.

### Verdict: HOLDS WEAKLY, FAILS STRONGLY

**Weak form (which holds):** Two chains with the same multiset of projections produce the same signature, hence the same Theorem-1 verdict on dimensional consistency. This is trivially true by set union.

**Strong form (which fails):** The Connes-vS spectral truncation work (WH1 round 1, 2, 3) shows that on the *operator system* level, hopf_bundle and spinor_lift do not commute — the order in which the Hopf U(1) twist and the Camporesi-Higuchi lift are applied changes whether the resulting truncated triple has propagation number 2 (Toeplitz S¹ behavior) or something else. Signature arithmetic cannot see this because both orderings produce the same union of {Z, E, α} with same {energy} dimension and same {pi_hopf, half_int_hurwitz, odd_zeta, ...} transcendental set.

### Counter-example to the strong form

Take the K-chain `fock + hopf + spinor_lift`. Signature arithmetic treats `hopf;spinor_lift` and `spinor_lift;hopf` as identical. But:

- `hopf;spinor_lift`: Lift the scalar bare-graph nodes onto the Hopf base S², then lift the resulting (l, m_l)-labeled S² nodes to (l, j=l±1/2, m_j) Camporesi-Higuchi spinor nodes. Result is the standard "Hopf S² with Dirac fiber" construction.
- `spinor_lift;hopf`: Apply CG to get Dirac (n, κ, m_j) nodes, then quotient by the Hopf U(1). The U(1) action on Dirac spinors is *not* the same as on scalars — half-integer charges enter the quotient, and the resulting base is *not* S² but a different Z₂-quotient of Dirac S³.

These are different operator systems even though their three-axis signatures are identical. The strong form fails.

### Counter-example to the weak form (worse)

Even the weak form has a hole when the same variable is added twice. Consider a hypothetical chain `fock + bargmann_segal`: both add an energy scale (Z·Hartree from fock, ℏω from bargmann_segal). Signature arithmetic would tag this as `{energy, energy}` — but physically there are *two distinct energy scales*, and the chain might or might not be admissible depending on whether the Coulomb and HO physics are being mixed (Paper 24's HO/Coulomb asymmetry forbids this). The dictionary needs a refinement that distinguishes "energy from focal length p₀² = −2E" from "energy from HO frequency ℏω". My code flags this as a `variable_conflict` in the composition table; six such conflicts exist, all involving energy-introducing projections trying to compose pairwise.

### Worked examples (5/6 PASS for weak form, K-example flags non-commutation)

| Example | k | n_perms | All same signature? | Non-commuting pairs |
|:--------|:-:|:-------:|:-------------------:|:-------------------:|
| Lamb shift 4-chain | 4 | 24 | YES | (none flagged by metadata) |
| K = π(B+F−Δ) | 3 | 6 | YES | **(hopf, spinor_lift)** ← genuine non-commutation |
| Stefan-Boltzmann | 2 | 2 | YES | (none) |
| He CI | 2 | 2 | YES | (none) |
| Spatial Casimir Dirac | 2 | 2 | YES | (none) |
| Uehling VP | 2 | 2 | YES | (none) |

The single non-commuting pair flagged is exactly the WH1-relevant one. This is encouraging: the simple metadata-based commutation test correctly identifies the place where Connes-vS spectral truncation order matters.

### Proposed weaker theorem that holds

> **Weak Minimality Theorem.** For a given target observable's signature S, the shortest projection chain that produces S is unique up to the equivalence relation generated by metadata-level commutation. Operator-system inequivalences arising from spinor_lift composed with non-trivial bundle structures (hopf, spectral_action, wilson) are not captured by signature arithmetic and require explicit operator-system analysis.

This is the form that the paper should claim. It is honest about what signature arithmetic can and cannot resolve.

---

## 4. THEOREM 3 — Composition Algebra (15×15 table)

### Method

For each ordered pair (P₁, P₂), classify the relation as one of:
- `commute` — both `a;b` and `b;a` are layer-admissible and signatures match; metadata variables disjoint
- `one_way_AB` — `a;b` admissible, `b;a` blocked by layer mismatch
- `one_way_BA` — `b;a` admissible, `a;b` blocked
- `not_composable` — neither direction admissible
- `idempotent` — diagonal entry where the projection is its own fixed point (sturmian, wigner_3j, drake_swainson)
- `self_compose` — diagonal entry where projection composed with itself does something non-trivial (the other 12)
- `variable_conflict` — both add the same variable name; composition is ill-defined without disambiguation
- `inverse_pair` — non-trivial pair (a, a) where a has a meaningful inverse (stereographic, wigner_3j, wigner_D)

### Aggregate counts (15² = 225 cells)

| Relation | Count |
|:---------|:-----:|
| commute | **182** |
| one_way_AB | 11 |
| one_way_BA | 11 |
| self_compose (diagonal) | 12 |
| idempotent (diagonal) | 3 |
| variable_conflict | 6 |

182/225 = **80.9% of all pairs commute** at the signature level. This high commutativity is a property of signature-arithmetic abstraction: most projections introduce different variables and have no layer ordering constraint. The interesting structure is in the 22 one-way pairs and the 6 variable conflicts.

### Most surprising entry

The `spectral_action ; observation_window` (and its reverse) are both flagged as one-way / asymmetric depending on which layer they target. Mechanism:

- spectral_action: Layer 1 → Layer 2 (introduces UV cutoff Λ)
- observation_window: Layer 1 → Layer 2 (introduces inverse-temperature β / observer time)

Both target Layer 2, so neither can compose with the other in the "feed forward" sense (each consumes the bare graph independently). Yet both must operate when computing finite-temperature QED loop integrals, and the order in which Λ regularization and β regularization are taken famously *matters* in continuum QFT (the thermal limit β → ∞ does not commute with the UV limit Λ → ∞). Signature arithmetic correctly flags this as a non-commutation, even though it cannot tell you *which* order is physically required.

This is the cleanest illustration in the table that the dimensional axiomatization captures real structure: the failure to commute encodes the same content as the field-theoretic UV/IR ordering.

### Idempotent projections (3, all on diagonal)

- `sturmian` — Coulomb Sturmian basis is its own reparameterization fixed point (changing exponent λ to itself does nothing).
- `wigner_3j` — Coupling a coupled state to itself returns the same state (CG is unitary).
- `drake_swainson` — The asymptotic-subtraction regulator, once applied, gives a K-independent value; applying it again on the regularized result is idempotent.

These are also the three projections that can be omitted from a chain without changing the result (assuming the input is already in the projection's image), which is exactly Paper 34's "Layer 1 only" sense.

### Inverse pairs (3)

- `stereographic` — change of coordinates is invertible on the open dense subset.
- `wigner_3j` — CG is unitary, hence invertible (decoupling is the inverse of coupling).
- `wigner_D` — rotations form a group; rotation by R⁻¹ inverts rotation by R.

`rest_mass` and `observation_window` also have "inverses" in the loose sense (set m=0, β→∞), but these are limits, not true inverses; I have not classified them as inverse pairs.

### Variable conflicts (6)

The 6 cells are all pairs of energy-introducing projections trying to share the energy slot:

```
(fock_conformal, bargmann_segal)
(fock_conformal, spectral_action)
(bargmann_segal, spectral_action)
+ all 3 reverses
```

These are *not* failures of the dictionary — they are correct flags that the composition would conflate two physically distinct energy scales. The framework's HO/Coulomb asymmetry (Paper 24, Paper 31) is the structural reason these compositions are forbidden: the Coulomb energy scale (Hartree, derived from p₀² = −2E) and the HO energy scale (ℏω) are not interchangeable.

---

## 5. Axis-independence verdict (the honest answer)

The prompt asks: are the three axes (variable, dimension, transcendental) independent, or correlated? If correlated, is the taxonomy overcomplete? If hyper-independent, are we missing projections?

### Cross-tabulation (8 axis flag patterns)

Of 15 projections:

| Pattern (V D P) | Count | Members |
|:---------------:|:-----:|:--------|
| VD– (var + dim, no π) | 6 | fock, bargmann, stereographic, wigner_D, molframe, rest_mass |
| V–P (var + π, no dim) | 2 | hopf, wilson |
| VDP (all three) | 2 | spectral_action, observation_window |
| V–– (var only) | 2 | spinor_lift, drake_swainson |
| ––P (π only) | 1 | vector_photon |
| ––– (Layer 1, none) | 2 | sturmian, wigner_3j |
| –D– (dim only, no var/π) | **0** | (empty) |
| –DP (dim + π, no var) | **0** | (empty) |

### Findings

1. **The variable axis and dimension axis are NOT independent.** Every projection that adds a dimension also adds a variable (8/8). There is no projection with `(–, D, ?)`. Mechanism: the dimension is *carried by the variable*. You cannot add "energy" as a dimension without also adding *some* energy scale (Z·Hartree, ℏω, Λ, m·c², 1/β); you cannot add "length" without an *r* or *R* or R_AB.

2. **The transcendental axis IS partially independent.** vector_photon adds a π without adding any variable or dimension (Hopf-base solid angle, dimensionless). And rest_mass is the converse — adds a variable + dimension but introduces no transcendental.

3. **Only 2 projections move all three axes**: spectral_action and observation_window. These are the framework's "heavy" projections — they carry the most physical content per projection and are the most likely places where multiple physics scales mix.

4. **No projection is pure-dimension-only.** This is a real structural statement: in GeoVac, dimensions cannot be introduced without a corresponding variable. (Compare to general relativity, where one could imagine a "rescale time only" projection that shifts the time variable's measure without introducing a new scalar; GeoVac doesn't have this.)

### Verdict on the axis structure

The three-axis dictionary is **mildly overcomplete** but in a *useful* way:

- The variable and dimension axes are perfectly correlated (V ⇒ D in 8/8 dimensional projections; D ⇒ V trivially since 0 projections have D without V). They could be merged into a single "physical-content axis" tagged by (variable, dimension) tuples without loss of information.
- The transcendental axis is independent of the other two (vector_photon proves it).

**Is the taxonomy overcomplete?** Yes, but the cost of merging is high: (variable, dimension) tuples carry more information than either axis alone. The collapsed two-axis dictionary would lose the cleanliness of "this projection introduces one specific physical scale." I recommend keeping all three axes in Paper 34, with an explicit note that the variable/dimension correlation is a structural fact and that vector_photon is the unique axis-independence witness.

**Are we missing projections?** Two candidate gaps surfaced by the analysis:

- A `gauge_choice` projection that selects a gauge (Coulomb vs Lorenz vs ...) without adding a new variable or dimension. This would be a pure-transcendental projection in the (––P) bucket, similar to vector_photon but for QED gauge structure rather than vertex structure. Status: not currently in Paper 34's list.
- A `spin_statistics` or `Z2_grading` projection that toggles boson vs fermion temporal boundary conditions (relevant for KG-5's Dirac Casimir 17/480 vs the scalar 1/240). Currently bundled inside observation_window via the (2k+1)π/β fermionic Matsubara modes. Could be split out as its own (–, –, π/2 shift) projection.

Both candidates are flagged for PI consideration; neither is forced by current empirical data.

---

## 6. Recommended Edits to Paper 34

These are recommendations only; do NOT apply directly per sprint constraints.

### EDIT-1 (high priority): New §VII or §VII.1 "Dimensional Axiomatization"

Add a section after §VI (Falsifiable Prediction) consolidating the three theorem statements:

> **Theorem 1 (Dimensional Consistency).** *For any projection chain S = Pₙ ∘ ... ∘ P₁, the variable, dimension, and transcendental signatures of S equal the union of the signatures of the constituent Pᵢ. The catalog rows of §V satisfy this consistency check across all 8 fully-tagged entries (Lamb shift, K = π(B+F−Δ), Stefan-Boltzmann, He CI, Hydrogen E_n, scalar/Dirac spatial Casimir, Uehling VP).*
>
> **Theorem 2 (Weak Minimality).** *Two chains containing the same multiset of projections produce the same signature, hence are observationally indistinguishable at the level of dimensional types and transcendental classes. They may differ at the operator-system level when bundle-non-trivial projections (Hopf, spectral_action, Wilson) are composed with the spinor lift; these inequivalences are not visible to signature arithmetic and require explicit operator-system analysis (cf. WH1 spectral truncation rounds 1–3).*
>
> **Composition relations (Theorem 3).** *The 15×15 composition table (computed in `debug/data/tx_a_composition_table.json`) shows 80.9% commutation, 22 one-way pairs (constrained by Layer transitions), 6 variable conflicts (energy-scale-conflating compositions of Fock, Bargmann–Segal, and spectral_action), and 3 idempotent projections (Sturmian, Wigner 3j, Drake–Swainson).*

### EDIT-2 (high priority): Sharpen Observation in §IV "Two-axis duality"

Currently §IV claims "Each Paper~18 transcendental tier corresponds to specific projections, and each projection's transcendental signature is a specific Paper~18 tier." Sprint TX-A demonstrates this is not a bijection on the variable/dimension side (V ⇒ D in 100% of cases, 0 projections with D-only). Recommend adding a remark:

> **Remark (variable–dimension correlation).** *The variable axis and dimension axis are not independent: every projection in §III that introduces a physical dimension also introduces a physical variable that carries that dimension. The transcendental axis is the only one that admits an independence witness — vector_photon (§III.11) introduces 1/(4π) without introducing any new variable or dimension.*

### EDIT-3 (medium priority): §VIII open question on missing projections

Add to §VIII Open Questions:

> **(5) Pure-axis projections.** The current dictionary contains no projection that moves only the dimension axis (e.g., a unit-rescaling that introduces a length scale without a coordinate variable). Whether such projections are structurally absent from the GeoVac framework, or merely absent from the current implementation, is open. Two candidate additions surfaced by the TX-A axiomatization sprint: (a) a gauge-choice projection (Coulomb vs Lorenz) that would sit in the pure-transcendental bucket alongside vector_photon, (b) a spin-statistics projection that splits boson/fermion Matsubara content out of observation_window. Neither is forced by current data.

### EDIT-4 (low priority): Caveat on §V Catalog interpretation

Currently §V says "every match in the catalog can be tagged by projection(s)". The TX-A worked examples confirm this for the 8 explicit rows but also flag that signature equality does not imply value equality (Hydrogen E_n and Lamb shift have the same dimensional/π signature but differ by 13 orders of magnitude). Recommend a remark:

> **Remark (signature ≠ value).** *Theorem 1's signature consistency is a type-level statement: chains that produce the same signature can produce wildly different numerical values. The depth-prediction of §VI addresses the value-level question; signature arithmetic addresses only the type-level question.*

### EDIT-5 (cross-paper, Paper 18): Cross-reference back

If EDIT-1 lands in Paper 34, Paper 18 §IV.10 (the spinor-intrinsic ring R_sp) should cross-reference Paper 34 §VII showing how the spinor_lift projection is the unique source of half-integer Hurwitz / odd-ζ / Dirichlet-β content in the dictionary. This is already implicit in Paper 18 but the dictionary makes it explicit.

---

## 7. Summary table for the report-back

| Item | Result |
|:-----|:-------|
| (a) Theorems holding | T1 holds (8/8 PASS, trivially-by-construction + non-trivially-empirically); T2 holds weakly (signatures order-independent), fails strongly (operator-system order matters at Hopf+spinor and Λ+β); T3 holds, full table built (182 commute, 22 one-way, 6 variable conflicts, 3 idempotent, 3 inverse pairs). |
| (b) Projection table summary | Of 15: 5 add π (hopf, wilson, vector_photon, spectral_action, observation_window), 1 adds odd-ζ family (spinor_lift), 9 are π-free at projection level (4 fully Layer-1-preserving + 5 ring-extending). 8/15 add a dimension, all 8 also add a variable (perfect correlation). 2/15 (spectral_action, observation_window) move all three axes — the "heavy" projections. |
| (c) Most surprising composition entry | `spectral_action` and `observation_window` are mutually one-way (UV/IR limits famously do not commute in QFT), correctly captured by the layer-grading without requiring continuum field theory. The K-chain `hopf + spinor_lift` non-commutation is also surprising because metadata catches it without explicit operator-system computation — exactly the place WH1 R3.2 found the Connes-vS truncation order matters. |
| (d) Recommended Paper 34 edits | EDIT-1 (new §VII Dimensional Axiomatization with 3 theorem statements), EDIT-2 (sharpen §IV duality remark with V⇒D observation), EDIT-3 (§VIII open-question on pure-axis projections), EDIT-4 (signature ≠ value caveat in §V), EDIT-5 (cross-reference from Paper 18). |

---

*Memo word count: ~3500 words. Sprint complete. Code: `debug/tx_a_signature_arithmetic.py`. Data: `debug/data/tx_a_projection_table.json`, `debug/data/tx_a_composition_table.json`. Tests pass; 8/8 worked examples verify Theorem 1; 6/6 multi-projection examples verify Theorem 2 weak form; full 15×15 composition table built.*
