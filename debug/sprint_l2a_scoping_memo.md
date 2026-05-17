# Sprint L2-A — BBB-Krein-Lift Scoping Memo

**Date:** 2026-05-16
**Sprint:** L2-A (literature audit + sprint planning, NO code, NO paper edits)
**Status deliverable:** memo + structured bibliography JSON + sequential sprint plan
**Builds on:** Sprint L0 audit (`debug/lorentzian_l0_audit_memo.md`), Paper 42 §3 Nieuviarts NO-GO, lorentz-boost NO-GO (`debug/lorentz_boost_scoping_memo.md`), lit update of 2026-05-16 (`debug/lorentzian_literature_update_2026_05_16.md`)
**Companion files:** `debug/data/sprint_l2a_bibliography.json`, `debug/sprint_l2a_sprint_plan.md`

---

## §1. Executive summary

**Verdict: GO-WITH-CAVEATS on the BBB Krein-lift path. Reachable as a 4-6 month sprint sequence if Sprint L3 (Lorentzian propinquity) is deferred. L3 itself is 6-12 months of original NCG-mathematics and is not blocking for the structural sub-goals (L2-B, L2-C, L2-D, L2-E).**

The audit changes the named path. The literature-update trigger of 2026-05-16 reported that the Nieuviarts twist-morphism "may shorten the 6-12 month Lorentzian-NCG-math sub-budget significantly." This turned out to be wrong for GeoVac: Nieuviarts Definition 2.2 (twist by grading) and its companion in arXiv:2402.05839 are **explicitly restricted to even-dimensional manifolds** — the author of those papers states this himself in §2.2 of arXiv:2502.18105 around eq (2.12): *"the twist by the grading is the only way to minimally twist the spectral triple of even-dimensional manifolds based on C∞(M) ⊗ ℂ². The fact that no such procedure for twisting odd-dimensional manifold's spectral triple have been found will be one of the justifications to focus on the study of even-dimensional manifolds."* Paper 42 §3.2 already records this NO-GO; this audit confirms it from the original.

What the audit found instead: the **load-bearing concrete prescription is NOT in BBB 2018 itself.** BBB 2018 is the *classification* — it tells you which (m, n) ∈ (ℤ/8)² signatures admit Connes-axiom-compatible structures, gives the sign tables, and proves tensor-product additivity. It does NOT contain a Riemannian → Lorentzian lift procedure or a worked example on a compact spatial slice × ℝ. The QED example in BBB §8 is on Lorentzian 4-D Minkowski (m,n) = (4, 6), not S³ × ℝ.

The actual lift prescription that GeoVac would use is **van den Dungen 2016 Proposition 4.1** (arXiv:1505.01939, *Math. Phys. Anal. Geom.* 19, 4): given an even-dimensional pseudo-Riemannian spin manifold (M, g) of signature (t, s) with a spacelike reflection r, the Wick-rotated metric g_r = g(rv, w) makes (M, g_r) Riemannian, and the data (C∞_c(M), L²(S), i^t D̸, J_M) is an even Krein spectral triple of Lorentz-type iff t is odd. For GeoVac, take M = S³ × ℝ (4-dim, even spacetime; spatial S³ is 3-dim odd but doesn't need its own grading); t = 1; this gives a Lorentz-type triple at signature (3, 1).

Three caveats keep the verdict at GO-WITH-CAVEATS rather than clean GO:

1. **No published worked example of S³ × ℝ as a Lorentzian / Krein spectral triple.** The Prop 4.1 construction is generic over (M, g); no paper instantiates it for the GeoVac spatial slice or any compact-spatial × ℝ case. Sprint L2-C will need to actually write this down.
2. **No Lorentzian propinquity exists.** Latrémolière propinquity used in Paper 38 is fundamentally Riemannian (Lipschitz seminorm requires positive-definite inner product; Krein indefinite-norm fails the propinquity axioms). Sprint L3 is genuine original NCG-mathematics on the timescale of 6-12 months. The good news: Sprint L3 is NOT required for closing the (3,1) structural sub-goals at finite n_max (Sprints L2-B, L2-C, L2-D, L2-E land theorems-at-finite-cutoff without needing propinquity convergence to be established).
3. **S³ × S¹_β is structurally inadmissible as a (3,1) globally hyperbolic spacetime.** Compact Lorentzian manifolds without boundary contain closed timelike curves; the natural Matsubara compactification Sprint TD Track 1 / Track 4 uses on the Riemannian side does NOT Wick-rotate to a globally hyperbolic Lorentzian spectral triple. The Lorentzian version uses S³ × ℝ (non-compact in time); thermal physics enters via KMS states on the wedge sub-algebra rather than at the spectral-triple level. This is consistent with Paper 42's BW-α / BW-γ unified Wick-rotation theorem at the metric-functional level; just blocks the obvious "tensor with the Matsubara circle" shortcut.

Sprint L2 is **reachable, structurally honest, contains genuinely new mathematics (no derivative application of any single existing paper), and lands a clean Lorentzian reading of the WH1 PROVEN result at finite cutoff.** Estimated cost is 4-6 months if L3 is deferred (with named L3-prerequisite caveat on the GH-convergence statement), or 9-12 months if L3 is run in parallel.

---

## §2. BBB 2018 prescription, decoded

**Verified citation:** Bizi, Brouder, Besnard, *J. Math. Phys.* 59, 062303 (2018), arXiv:1611.07062. arXiv version is v2, dated 16 Oct 2017.

### §2.1 The classification

BBB classify any algebra carrying a quadruple (𝓗, χ, η, J) of [complex Hilbert space, self-adjoint chirality χ²=1, self-adjoint fundamental symmetry η²=1, antiunitary charge conjugation J with J†J = 1 and J² = ε = ±1] satisfying

- J² = ε
- Jχ = ε'' χ J
- Jη = εκ η J
- ηχ = ε''κ'' χ η

with four sign choices (ε, ε'', κ, κ''). The KO-dimension n ∈ ℤ/8 of Connes corresponds to (ε, ε''); BBB introduce a *second* dimension m ∈ ℤ/8 corresponding to (κ, κ''). The pair (m, n) classifies the data.

**The sign table (BBB Table 1):**

| m, n | 0 | 2 | 4 | 6 |
|:-----|:-:|:-:|:-:|:-:|
| κ, ε   | +1 | −1 | −1 | +1 |
| κ'', ε'' | +1 | −1 | +1 | −1 |

m and n are independently in ℤ/8 but BBB Table 1 only lists even values (this is the foundation of the "even-dim-only" restriction discussed in §2.3 below).

**Tensor product additivity (BBB Section 3):** for two mod-8-spacetime representations 𝒮₁, 𝒮₂ with signatures (m₁, n₁), (m₂, n₂),
- 𝒮 = 𝒮₁⊗̂𝒮₂ has signature (m, n) = (m₁ + m₂, n₁ + n₂) mod 8.
- The graded tensor product is constructed via χ = χ₁⊗̂χ₂, J = J₁χ₁^{|J₂|} ⊗̂ J₂χ₂^{|J₁|}, η = i^{|η₁||η₂|} η₁χ₁^{|η₂|} ⊗̂ η₂χ₂^{|η₁|}.

**(p, q) ↔ (s, t) ↔ (m, n) translation table (BBB Table 3, condensed):**
- m = t + s mod 8 ("total dimension")
- n = t − s mod 8 ("KO-dim")
- (j, k) and (j+4, k+4) give the same (m, n) (this is the Clifford algebra isomorphism Cℓ(s, t+8) ≃ Cℓ(s+8, t) ≃ Cℓ(s+4, t+4)).

### §2.2 The Krein spectral triple definition

BBB Section 5: an **even-dimensional real indefinite spectral triple** consists of (𝒜, 𝒦, D, J, χ, η) where:

1. 𝒜 is a *-algebra represented on a Krein space 𝒦 with Hermitian form (·, ·) and fundamental symmetry η, π(a*) = π(a)^× (Krein-adjoint).
2. χ is a chirality operator (χ² = 1, χ commutes with 𝒜).
3. J is antilinear with J†J = 1.
4. (ε, ε'', κ, κ'') describes the (m, n) signature.
5. D is Krein-self-adjoint (D^× = D), JD = DJ, χD = −Dχ.

The Krein-adjoint and Hilbert-adjoint are related by T† = η T^× η.

### §2.3 The dimensional restriction (HEADLINE)

BBB classification is **even-dimensional only.** The (m, n) sign table (Table 1) gives only even entries; the construction via gamma matrices uses Cℓ(p, q) with p+q = 2ℓ even; the explicit indefinite spectral triple in §5 is "even-dimensional real indefinite spectral triple." Section 9 of BBB acknowledges:

> "For applications to topological insulators, it would be desirable to extend these results to the case of odd-dimensional algebras."

This is the same odd-dim obstacle that affects Nieuviarts. **For GeoVac at signature (3, 0) (spatial S³, KO-dim 3, odd), this means BBB does NOT classify the spatial-only triple.** But this is NOT fatal: GeoVac does not need the spatial-only triple at (3, 0) to fit BBB; what GeoVac needs is the *spacetime* triple at (3, 1) — which is m + n combined = 4-dimensional (even). The spatial S³ inherits its description from the spacetime restriction.

### §2.4 What BBB does NOT contain

- No Riemannian → Lorentzian lift procedure. BBB gives the (m, n) sign table and the abstract Connes-axiom signature, but no procedure that takes a Riemannian spectral triple at signature (m, 0) as input and produces a Lorentzian spectral triple at signature (m', 1) as output.
- No (3, 1) example. The QED example in BBB §8 is on 4-D Minkowski with manifold signature (s, t) = (3, 1) but the *spectral-triple-level* signature (m₁, n₁) = (4, 6) — i.e., it sits at the "even-only" entries of the table. No (3, 1) signature triple is constructed.
- No compact spatial × time example. Nothing of the form N × ℝ with N a compact spatial manifold.
- No propinquity, GH-convergence, or finite-cutoff truncation discussion. BBB is an axiomatic classification paper.

**Net assessment of BBB 2018 for GeoVac:** the (m, n) classification is needed for naming the signature of the GeoVac (3, 0) → (3, 1) lift, but BBB does not give the construction. The construction recipe must come from elsewhere.

---

## §3. Literature landscape: where the construction recipe actually is

### §3.1 The load-bearing reference: van den Dungen 2016

Van den Dungen, *Math. Phys. Anal. Geom.* 19, 4 (2016), arXiv:1505.01939. **Proposition 4.1** is the clean Riemannian → Lorentzian (Krein) lift:

> Let (M, g) be an n-dimensional time- and space-oriented pseudo-Riemannian spin manifold of signature (t, s). Let r be a spacelike reflection on TM (TM = E_t ⊕ E_s decomposition; r acts as −1 on E_t, +1 on E_s). Let g_r = g(rv, w) be the Wick-rotated metric; (M, g_r) is Riemannian. Then (C∞_c(M), L²(𝕊), i^t D̸, J_M) is an even Krein spectral triple, with grading Γ_M; Lorentz-type (i.e., J_M odd) iff t is odd.

This is the construction GeoVac needs. For (M, g) = (S³ × ℝ, ds²_S³ − dt²) at signature (s, t) = (3, 1):
- Wick-rotated metric g_r is ds²_S³ + dt² → (M, g_r) = (S³ × ℝ, standard Riemannian) which is well-defined.
- t = 1 is odd → Lorentz-type triple.
- The Dirac operator is i D̸ where D̸ is the standard Dirac on the Wick-rotated Riemannian S³ × ℝ.
- The fundamental symmetry 𝒥_M = γ(e₀) where e₀ is the unit timelike vector.
- The grading Γ_M is the spacetime chirality γ^{5} = i γ⁰γ¹γ²γ³.

**Dimensional verification for GeoVac:** spacetime is 4-dim (even); t = 1 (odd) → Lorentz-type. The spatial S³ is 3-dim (odd) but does NOT need its own grading; the spacetime grading Γ_M lives on the full S³ × ℝ.

**Caveat:** van den Dungen 2016 gives the abstract construction. It does NOT instantiate it for any specific N × ℝ with N a compact spatial manifold. Sprint L2-C will need to actually instantiate this for N = S³.

### §3.2 The Kasparov-product alternative: van den Dungen–Rennie 2016

Van den Dungen, Rennie, *Ann. Henri Poincaré* 17, 3255 (2016), arXiv:1503.06916. Constructs the Lorentzian spectral triple as a *family* of Riemannian spectral triples indexed by a foliation of globally hyperbolic spacetime by spacelike hypersurfaces. The Kasparov product of the family gives an indefinite Kasparov module which is the Lorentzian analog of a spectral triple.

This approach is closer in spirit to what GeoVac does at finite n_max: each spatial slice S³_t carries a Camporesi–Higuchi spectral triple T_{n_max}, and the Lorentzian object is the "stack" of these over the time direction. The abstract construction is in 1503.06916; for GeoVac it would have to be instantiated.

Three concrete examples covered in the paper (per abstract): (1) Dirac on a pseudo-Riemannian spin manifold; (2) the harmonic oscillator; (3) Kasparov-product construction from a family. The S³ × ℝ case is structurally compatible but not worked out explicitly.

### §3.3 The temporal-Lorentzian alternative: Franco–Eckstein 2014

Franco, Eckstein, *Rev. Math. Phys.* 26, 1430007 (2014), arXiv:1210.6575. Adds a 3+1 decomposition and a global-time element T ∈ 𝒜 to a pseudo-Riemannian spectral triple. Provides a Lorentzian distance formula

d(p, q) = sup_{a ∈ 𝒜 causal, [D, a] Krein-bounded} (a(q) − a(p)).

The 3+1 decomposition is *required* — it's part of the definition. For GeoVac, the natural 3+1 decomposition is the obvious one (spatial S³ + temporal ℝ).

This framework is the right venue for the Lorentzian-distance / Connes-distance-at-(3,1) question, which corresponds to the R2.3-leg work of WH1.

### §3.4 Foundations: Strohmaier 2006

Strohmaier, *J. Geom. Phys.* 56, 175 (2006), arXiv:math-ph/0110001. The foundational paper defining pseudo-Riemannian spectral triples. Shows that for noncommutative tori and pseudo-Riemannian spin manifolds, dimension, signature, and integral are recoverable from spectral data. Construction is abstract; no compact-spatial × ℝ example.

Strohmaier is the upstream reference cited by van den Dungen, Franco–Eckstein, BBB, and all subsequent work. Its role for L2 is naming the canonical definitions.

### §3.5 Twisted-emergence alternatives: Devastato–Lizzi–Martinetti 2018, Nieuviarts 2024–2026

Devastato, Farnsworth, Lizzi, Martinetti, JHEP 03, 089 (2018), arXiv:1710.04965. Shows that twisting the SM spectral triple by an inner-factor twist naturally yields the Krein space associated with Lorentzian signature. The twist acts on the inner factor (the SM finite spectral triple A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ)) — NOT on the outer manifold. This means for GeoVac, the twist mechanism could in principle apply to the inner factor of the AC extension (Sprint H1), but not to the outer S³ side.

Nieuviarts 2024–2026 (arXiv:2402.05839, 2502.18105, 2512.15450) extends this: the twist-morphism produces a pseudo-Riemannian spectral triple from a Riemannian twisted spectral triple "as an algebraic alternative to Wick rotation." The signature change is determined by a unitary K with K² = 1.

**Direct check of the odd-dim caveat in arXiv:2502.18105 §2.2:** The author explicitly states *"the twist by the grading is the only way to minimally twist the spectral triple of even-dimensional manifolds based on C∞(M) ⊗ ℂ². The fact that no such procedure for twisting odd-dimensional manifold's spectral triple have been found will be one of the justifications to focus on the study of even-dimensional manifolds."* For GeoVac at S³ = SU(2) (odd-dim), the twist-morphism is NOT applicable directly. This confirms the Paper 42 §3.2 footnote NO-GO.

**Possible extension:** the twist-morphism is applicable on the S³ × ℝ spacetime (4-dim, even). If Sprint L2 wants to use the twist-morphism, it would be on the spacetime, not on the spatial S³. The twist mechanism for spacetime SM almost-commutative geometry could in principle be applied — this is a parallel path that should be flagged for evaluation.

### §3.6 The propinquity gap: Latrémolière 2022 and follow-ons

Latrémolière, *Adv. Math.* 404, 108393 (2022), arXiv:1811.10843. **STRICTLY RIEMANNIAN.** The Lipschitz seminorm L(a) = ‖[D, a]‖ used to define the propinquity requires the operator norm derived from a positive-definite Hilbert inner product. The Krein indefinite norm does NOT satisfy this — there is no positive seminorm that captures the same Lipschitz-class information.

Recent (2024–2025) extensions — Toyota 2025 (analytic-path-of-Riemannian-metrics, arXiv:2504.11715), Hekkelman–McDonald 2024 (compact-quantum-metric-spaces truncations, arXiv:2412.00628) — all extend the propinquity along the Riemannian side. **No Lorentzian propinquity exists.** Constructing one is an original NCG-mathematics contribution.

This is the hard wall for Sprint L3. The good news: Sprint L2-B / L2-C / L2-D / L2-E do not require Lorentzian propinquity to land their structural results at finite n_max.

### §3.7 Compact-spatial-times-time: not in the published literature

The audit found **no published example** of a compact spatial manifold N (with N = S^n for n ≥ 2, or any compact n-manifold) being instantiated as a (n, 1) Lorentzian / Krein spectral triple via N × ℝ. The closest published worked examples are:
- Rennie–Strung 2014 (arXiv:1411.3578): globally hyperbolic 2-D surfaces (i.e., (1, 1) signature, very minimal).
- van den Dungen–Paschke–Rennie 2013 (arXiv:1207.2112): harmonic oscillator on (n, 1) pseudo-Riemannian manifolds, but the abstract treatment.
- BBB 2018: QED on Minkowski (3, 1) (non-compact spatial).

None of these is S³ × ℝ or any compact-spatial × ℝ at (n, 1). **Sprint L2-C is therefore writing down a new explicit construction**, not applying an off-the-shelf one. This is structurally honest scope, and a natural mathematical contribution.

### §3.8 Globally hyperbolic obstruction (THE FUNDAMENTAL OBSERVATION)

A Lorentzian manifold without boundary that is compact contains closed timelike curves (Geroch's theorem). This means:

- **S³ × ℝ_t** is globally hyperbolic (S³ is a compact Cauchy surface; ℝ_t is non-compact). ✓ admissible.
- **S³ × S¹_β** is NOT globally hyperbolic; closed timelike curves exist. ✗ inadmissible at the spectral-triple level for a *Lorentzian* signature.

**Implication for GeoVac.** Sprint TD Track 1 / Track 4 use M = S³ × S¹_β with β interpreted as inverse temperature on the Riemannian side. Wick-rotating S¹_β to S¹_t gives the same compact-time topology, which would have CTCs in Lorentzian signature. Hence the natural Matsubara compactification GeoVac uses on the Riemannian side does NOT have a direct Lorentzian spectral-triple lift.

**The correct reading,** consistent with Paper 42's four-witness Wick-rotation theorem (Hawking + Sewell + Bisognano-Wichmann + Unruh), is: the Lorentzian object is M = S³ × ℝ_t (non-compact in time), and *thermal* physics enters via a KMS state on the wedge sub-algebra (the modular flow / Bisognano-Wichmann reading). This is exactly the operator-system-level extension that Paper 42 sets up.

This is a structural caveat for the sprint plan, not a blocker. Sprint L2-B (Krein space construction) and L2-C (Lorentzian Dirac) build on S³ × ℝ_t. Sprint L2-E (literal Wick-rotation identification at Krein level) re-derives the Paper 42 unified four-witness theorem on this Lorentzian object.

---

## §4. Technical obstacles for GeoVac specifically

Walking through each Connes axiom and structural ingredient of GeoVac at the (3, 0) → (3, 1) lift:

### §4.1 Algebra 𝒜_GV

Function algebra on the Fock-projected S³ graph (Paper 32 §III). Tensor with 𝒜_t = C^∞(ℝ_t) (or compactly supported sections) to get 𝒜_GV ⊗ 𝒜_t. **Wick-rotation-trivial.** No GeoVac-specific obstacle.

### §4.2 Hilbert space → Krein space

Current 𝓗_GV is the truncated Camporesi–Higuchi spinor space at n_max (Paper 32 §III). Krein lift: tensor with L²(ℝ_t), use fundamental symmetry 𝒥 = γ⁰ on the spinor bundle of S³ × ℝ. Becomes Krein space 𝒦 = (𝓗_GV ⊗ L²(ℝ_t), η = γ⁰).

**GeoVac-specific question:** is the Camporesi–Higuchi (3, 0) spinor bundle on S³ structurally compatible with the (3, 1) gamma matrix algebra Cℓ(3, 1)? Cℓ(3, 1) ≃ M_4(ℝ) but Cℓ(3, 0) ≃ M_2(ℂ) ⊕ M_2(ℂ) (or ℍ ⊕ ℍ in real form). The spinor dimensions match (4 = 2+2), but the embedding of the spatial gamma matrices is not unique — needs to be chosen consistently with the Camporesi–Higuchi spectral decomposition.

**Sprint L2-B deliverable.** This is concrete construction work, 2-3 weeks.

### §4.3 Dirac operator

Current D_GV is the truncated Camporesi–Higuchi Dirac with spectrum |λ_n| = n + 3/2 on round S³. Lorentzian D = γ⁰ ∂_t + D_{S³} on S³ × ℝ.

**GeoVac-specific subtlety #1:** the Lorentzian Dirac has continuous spectrum because of the non-compact ℝ_t direction. Compact resolvent / summability of the Riemannian spatial part is preserved, but the *overall* compact resolvent fails.

Resolution (per van den Dungen 2016 Prop 4.1 and van den Dungen–Rennie 2016 framework): work with the *family* of (3, 0) spatial spectral triples indexed by t ∈ ℝ, and use Kasparov product / temporal compactification techniques. The full (3, 1) "spectral triple" is not literally a Connes spectral triple; it's an indefinite Kasparov module.

**Sprint L2-C deliverable.** Construct the explicit Lorentzian Dirac operator D_L on (S³ × ℝ, ds² − dt²) using the Wick-rotation lift D_L = i^t D̸_g_r where D̸_g_r is the Riemannian Dirac on (S³ × ℝ, ds² + dt²). 3-4 weeks.

**GeoVac-specific subtlety #2:** the four-witness Wick-rotation theorem in Paper 42 (Hawking + Sewell + BW + Unruh) restricts to a hemispheric wedge of S³ at finite cutoff. The Lorentzian analog is the wedge sub-algebra 𝒜_W ⊂ 𝒜_L of S³ × ℝ, with modular flow given by the boost generator K_α. Sprint L2-E will re-derive the unified-strong theorem at the Krein level.

### §4.4 Real structure J

GeoVac Riemannian J has J² = −I, JD = +DJ at KO-dim 3 (Paper 32 §IV verified bit-exact). Under (3, 0) → (3, 1):
- BBB Table 1 entry at (m, n) = (4, 6) (corresponding to (s, t) = (3, 1) per Table 3 East-coast convention) gives ε = +1 (J² = +I), ε'' = −1 (Jχ = −χJ), κ = +1 (Jη = +ηJ), κ'' = −1.
- This is the standard (3, 1) charge conjugation: J² = +I in Lorentzian East-coast, J² = −I in Lorentzian West-coast.

**Sign convention is non-trivial.** Need to fix West-coast vs East-coast (BBB §2) at the outset. GeoVac's current Riemannian J² = −I is more consistent with West-coast convention if the Lorentzian lift uses J² = +I at (3, 1), but consistency check needed at sign-table level.

**Sprint L2-D deliverable:** redo Paper 32 §IV Connes axiom audit at signature (3, 1) per BBB sign table, verify all four signs (ε, ε'', κ, κ'') and the modified Connes axioms (J², JD, χD) at finite n_max ∈ {1, 2, 3}. 2-3 weeks.

### §4.5 Order-one / order-zero

Currently trivial because 𝒜_GV is commutative. Krein-side, order-one becomes [[D, π(a)], π(b)°] = 0 with π(b)° = η π(b)* η. For commutative 𝒜, the property is structurally preserved. **No GeoVac-specific obstacle.**

### §4.6 Latrémolière propinquity (the WH1 PROVEN axiom)

**This is the hard wall.** Paper 38's five-lemma propinquity-convergence proof is Riemannian. No Lorentzian propinquity exists; constructing one is Sprint L3 (6-12 months).

**Sprint L2-B/C/D/E do NOT require a Lorentzian propinquity** at finite n_max. They produce structural theorems at finite cutoff. The Lorentzian analog of WH1 PROVEN (convergence T_{n_max,L} → T_{S³ × ℝ}^L in Lorentzian propinquity) would be the Sprint L3 deliverable. Defer.

### §4.7 Master Mellin engine M1 / M2 / M3

- **M1 (Hopf-base measure 2π = Vol(S¹)):** Wick-rotation-trivial via the four-witness theorem. M1 sub-mechanism transfers cleanly.
- **M2 (Seeley–DeWitt heat kernel):** Tr e^{-tD²} requires elliptic D. Lorentzian D² is hyperbolic; heat kernel becomes Schwinger proper-time integral ∫₀^∞ dt e^{-itD² + iεt²/2}. Coefficient structure changes; the Seeley–DeWitt √π·ℚ ring is not the right ring for Lorentzian-side coefficients without rederivation.
- **M3 (vertex-parity Hurwitz / Catalan G):** Riemannian eigenvalue parity condition on the Dirac spectrum. Krein-side analog: track 𝒥-positive vs 𝒥-negative sectors. **Sprint L2-D prediction:** on chirality-symmetric spectrum at (3, 1), M3 trivializes (the (3, 1) BBB Table 2 entry flips the {J, γ_5} anticommutation table, which has the effect of making the M3 vertex-parity sum identically zero on the chirality-symmetric truncation).

This M3 trivialization prediction is **a clean Sprint L2-E falsifier**: the case-exhaustion theorem (Paper 32 §VIII) predicts M1/M2/M3 are the only π-sources; if M3 trivially closes at (3, 1) on the chirality-symmetric truncation, the Lorentzian-side case-exhaustion theorem has only M1/M2 active, which is a structurally simpler statement than the Riemannian one. Verifiable as a finite-n_max computation.

### §4.8 Roothaan multipole termination L_max = 2 ℓ_max

Currently exact via Wigner-3j triangle inequality. Under a Lorentz boost, angular momentum mixes with linear momentum (Wigner rotations); spatial L is no longer conserved.

**Status:** the question "what is the Lorentzian analog of Roothaan termination?" is open in the published literature. For *stationary observers* (no boost), the termination is preserved because the spatial S³ structure is preserved. For boosted observers, the termination is broken.

GeoVac's current precision-catalogue work (Sprint MH, Sprint HF, etc.) uses static atomic configurations — i.e., stationary observers — and the Roothaan termination remains valid. Lorentzian extension at the boost level (as opposed to just the signature level) is beyond Sprint L2 scope.

### §4.9 Hopf base measure 2π

Riemannian, unchanged. The Hopf base S² is spatial; Wick rotation doesn't touch the spatial fibration. **No GeoVac-specific obstacle.**

### §4.10 Summary table of obstacles

| Structural piece | Wick verdict | GeoVac specific issue | Sprint deliverable |
|:----------------|:-----------:|:----------------------|:-------------------|
| Algebra 𝒜_GV | trivial | none | — |
| Hilbert → Krein space | extends | spatial Camporesi–Higuchi spinor lift via Cℓ(3, 1) | L2-B |
| Dirac operator | extends | i^t D̸ via Prop 4.1 lift on S³ × ℝ | L2-C |
| Real structure J | extends | BBB Table 1 sign verification at (m,n) per (3,1) | L2-D |
| Bounded commutators | trivial | none | — |
| Compact resolvent | extends | family-of-spatial-triples / Kasparov | L2-C |
| Order-one / zero | trivial | none | — |
| **Latrémolière propinquity** | **BLOCKS** | original NCG-math (B1) | **L3, deferred** |
| M1 master Mellin | trivial | Paper 42 four-witness already covers | — |
| M2 master Mellin | blocks | Lorentzian heat-kernel rederivation | L2-D (predicted minimal) |
| M3 master Mellin | extends | trivializes on chirality-symmetric (3,1) | L2-D / L2-E |
| Roothaan termination | extends (stationary) | preserved for stationary observers | — |
| Hopf base measure | trivial | none | — |
| Modular Hamiltonian / Wick rotation literal identification | extends | re-derive Paper 42 unified-strong at Krein level on wedge of S³ × ℝ | L2-E |

---

## §5. Recommended sprint plan (L2-B, L2-C, L2-D, L2-E, L2-F)

Sequential dispatch with named risks. Full PM briefs in `debug/sprint_l2a_sprint_plan.md`.

### §5.1 L2-B — Krein space construction on S³ × ℝ at signature (3, 1)

**Scope:** Construct the Krein space 𝒦_{n_max} = 𝓗_GV^{n_max} ⊗ L²(ℝ_t) with fundamental symmetry 𝒥 = γ⁰ (the temporal Dirac matrix), and the corresponding indefinite inner product. Verify Krein-axiom basics (𝒥² = 𝒥*𝒥 = I; indefinite inner product is hermitian; positive-definite splits 𝒦 = 𝒦⁺ ⊕ 𝒦⁻). At finite n_max compactify the temporal direction to ℝ_t → small bounded interval for numerical concreteness (NOT S¹_β — that's CTC-violating).

**Deliverable:** memo + new module `geovac/krein_space_construction.py` (~250 lines) + ~30 tests verifying Krein axioms at finite n_max ∈ {1, 2, 3}.

**Estimated wall-clock:** 3 weeks (1 week design, 1.5 weeks implementation, 0.5 weeks verification).

**Prerequisites:** none (independent of L2-C/D/E).

**Named risks:** (R1) the spatial gamma matrix embedding into Cℓ(3, 1) might not be unique up to similarity, and the choice could affect downstream sign computations. Resolution: fix West-coast convention upfront (BBB §2).

### §5.2 L2-C — Lorentzian Dirac operator on S³ × ℝ at finite cutoff

**Scope:** Implement the van den Dungen 2016 Proposition 4.1 lift D_L = i D̸_g_r on S³ × ℝ with Wick-rotated Riemannian metric g_r = ds²_S³ + dt². At finite n_max truncate the spatial direction (Camporesi–Higuchi at n_max) and the temporal direction (small grid on ℝ_t). Verify Krein-self-adjointness D_L^× = D_L and the anticommutation χ D_L = −D_L χ.

**Deliverable:** memo + new module `geovac/lorentzian_dirac.py` (~300 lines) + ~30 tests at finite n_max ∈ {1, 2, 3}, also verifying that the Riemannian limit recovers Paper 32 §III's D_GV bit-identically.

**Estimated wall-clock:** 4 weeks (2 weeks construction, 1 week verification, 1 week documentation).

**Prerequisites:** L2-B (Krein space).

**Named risks:** (R2) the temporal-direction truncation is non-canonical; choices of boundary conditions affect spectrum. Need to document the convention and verify the static / Wick-rotation invariance separately. (R3) Continuous spectrum is unavoidable in the ℝ_t direction; finite-cutoff computations need to handle this carefully (avoid spectral-density artifacts).

### §5.3 L2-D — Connes axiom audit at signature (3, 1)

**Scope:** Verify the Paper 32 §IV Connes axiom audit (J² = ?I, JD = ?DJ, χD = −Dχ, order-one, order-zero) at signature (3, 1) using BBB Table 1 sign conventions at the relevant (m, n) entry. Determine which entry — there are two candidates per convention:
- East-coast: (s, t) = (3, 1) → (m, n) = (4, 6) → ε = +1, ε'' = −1, κ = +1, κ'' = −1.
- West-coast: (p, q) = (1, 3) → (m, n) = (4, 6) → same signs.

The sign-table entry is the same in both conventions, so this is unambiguous. Verify J²(3,1) = +I bit-exact at finite n_max. **Predict M3 sub-mechanism trivialization** on chirality-symmetric Dirac spectrum at (3, 1) (per Sprint L0 prediction) and verify computationally.

**Deliverable:** memo + new module `geovac/connes_axiom_audit_31.py` (~200 lines) + ~25 tests + paper-edit recommendation for Paper 32 §IV (extend the audit table from KO-dim 3 only to KO-dim 3 + (m,n)=(4,6) Lorentzian).

**Estimated wall-clock:** 2 weeks.

**Prerequisites:** L2-B + L2-C.

**Named risks:** (R4) the new sign at J² = +I changes the structural reading of Paper 32 §IV; need to update CLAUDE.md §1.7 WH1 entry to reflect that Riemannian J² = −I and Lorentzian J² = +I are the two natural cases. Not a blocker, just bookkeeping.

### §5.4 L2-E — Operator-system-level Wick-rotation literal identification (Krein-level Paper 42 theorem)

**Scope:** Re-derive the Paper 42 unified-strong four-witness theorem (Hartle-Hawking + Sewell + Bisognano-Wichmann + Unruh) at the Krein level on a hemispheric wedge of S³ × ℝ. The Riemannian-side construction gives BW-α (geometric K_α) and BW-γ (Tomita-Takesaki K_TT); the Lorentzian-side construction should give the corresponding Lorentzian modular Hamiltonian K_L_α and K_L_TT, with σ_{2π}^L(O) = O bit-exact at finite n_max under the BW choice of local Hamiltonian H_local = K_L_α^W / β.

**Deliverable:** memo + new module `geovac/modular_hamiltonian_lorentzian.py` (~400 lines) + ~50 tests verifying σ_{2π}^L = id at finite n_max ∈ {1, 2, 3} (mirroring the Paper 42 Riemannian-side verification at machine precision); paper-edit recommendation for Paper 42 §10 open question O1 (Lorentzian extension to (3,1)) and Paper 32 §VIII.D Lorentzian-side closure.

**Estimated wall-clock:** 4 weeks (mirrors the Paper 42 Riemannian L1 sprint which was ~5 days end-to-end, but the Krein-side construction is genuinely new and will take longer).

**Prerequisites:** L2-B, L2-C, L2-D.

**Named risks:** (R5) The H_local = K_α^W / β choice in Paper 42 was a derived structural finding, not a uniqueness theorem; in the Lorentzian case, the natural local Hamiltonian could differ (spectral-action paradigm vs thermal-time paradigm; this is Paper 42 §7.2 / O3). Sprint L2-E might land a verdict that the Lorentzian H_local differs from the Riemannian one, which would be a refinement of Paper 42's load-bearing scope finding.

### §5.5 L2-F — Falsifiers (parallel, lightweight)

**Scope:** A 1-week documentation pass: codify the named falsifiers from L2-B/C/D/E so they can be tested as the sprints land. Includes:
- L2-B falsifier: do Krein axioms hold at finite n_max bit-exactly? (If not: the Camporesi–Higuchi spinor lift to (3, 1) is structurally incompatible.)
- L2-C falsifier: does the Riemannian limit of D_L recover D_GV bit-identically? (If not: the Wick-rotation lift fails.)
- L2-D falsifier: does J²(3,1) = +I as predicted by BBB? Does M3 trivialize on chirality-symmetric spectrum? (If not: a structural error in BBB application, or M3 has more residual content than expected.)
- L2-E falsifier: does σ_{2π}^L = id hold at finite n_max bit-exactly? (If not: the Lorentzian extension of the unified-strong theorem fails, which would be a clean structural negative.)

**Deliverable:** falsifier memo `debug/sprint_l2_falsifiers.md` + cross-reference table.

**Estimated wall-clock:** 1 week (concurrent with L2-B).

**Prerequisites:** none.

### §5.6 Sequencing summary

```
Week 1-3:    L2-B (Krein space)  +  L2-F (falsifiers, parallel)
Week 4-7:    L2-C (Lorentzian Dirac)
Week 8-9:    L2-D (Connes axiom audit at (3,1))
Week 10-13:  L2-E (Krein-level Paper 42 redo)
Week 14:     Synthesis memo + paper edits + standalone Lorentzian paper draft outline.
```

**Total wall-clock: 14 weeks (~3.5 months) for L2-B through L2-E + synthesis.** Add 1-2 weeks of buffer for unforeseen complications → **4 months realistic estimate.**

### §5.7 Cross-track risk

The biggest cross-track risk is that L2-B reveals a structural incompatibility between the Camporesi–Higuchi spatial spinor bundle and the Cℓ(3, 1) gamma matrix embedding — i.e., that the spatial gamma matrices used for D_GV cannot be promoted to spacetime gamma matrices in a way consistent with both the (3, 0) Riemannian KO-dim 3 signs and the BBB (4, 6) Lorentzian signs simultaneously. This would force a different spinor-bundle choice on the Riemannian side, which would re-open WH1 PROVEN.

**Mitigation:** L2-B includes a Riemannian-limit check (Section 5.1 above): the Riemannian limit of the (3, 1) Krein space must recover Paper 32 §III's 𝓗_GV bit-identically. If this fails at L2-B, escalate immediately — the verdict is no longer GO-WITH-CAVEATS but something needs to be rethought.

---

## §6. Decision on Sprint L3 (Lorentzian propinquity)

**Recommendation: DEFER Sprint L3.**

L3 is the construction of a Lorentzian / Krein analog of the Latrémolière propinquity. This is original NCG-mathematics on the 6-12 month timescale. The audit confirms that:

1. No published Lorentzian propinquity exists as of May 2026.
2. The Lipschitz seminorm definition that grounds Latrémolière propinquity is fundamentally Riemannian.
3. Recent (2024–2026) follow-ons (Toyota, Hekkelman–McDonald) all extend the Riemannian side.

**Critical observation:** L3 is NOT required for Sprints L2-B/C/D/E. Those four sprints all produce *finite-cutoff theorems at signature (3, 1)*:
- Krein space at finite n_max.
- Lorentzian Dirac at finite cutoff.
- Connes axioms verified at finite n_max.
- Operator-system-level Wick-rotation literal identification at finite cutoff.

The Lorentzian analog of WH1 PROVEN — convergence T_{n_max}^L → T_{S³ × ℝ}^L in some Lorentzian propinquity — would be the L3 deliverable. Without it, the Lorentzian-side structural theorems are at finite cutoff, with the GH-convergence claim sitting as named open question (Paper 42 §10 O1 already named).

**This is consistent with how Sprint L0 / Sprint L1 already operate**: they land structural correspondence at the metric-functional level and operator-system literal identification at finite cutoff, without claiming a continuum limit. Sprint L2 does the same on the Lorentzian side.

**Standalone Lorentzian-extension paper** drafted after L2-B/C/D/E land would have the same scope: structural theorems at finite cutoff at signature (3, 1), with the Lorentzian-propinquity convergence question named as open (and load-bearing for Sprint L3, but separate from the paper's main results).

**Revisit trigger for Sprint L3:** if and only if (a) Lorentzian propinquity is published independently by Latrémolière / Toyota / Hekkelman-McDonald or similar; or (b) the Lorentzian extension of WH1 PROVEN becomes the load-bearing deliverable for some specific application.

---

## §7. Honest unknowns

The following items could not be answered from the literature alone; they will need to be worked out as Sprint L2-B/C/D/E proceed.

1. **Whether the spatial Camporesi–Higuchi spinor bundle on S³ is structurally compatible with the (3, 1) gamma matrix embedding** is not a literature question — no paper does S³ × ℝ as a Lorentzian spectral triple. This will be discovered during Sprint L2-B.

2. **Whether the Lorentzian Dirac D_L at signature (3, 1) on S³ × ℝ has the spectral structure needed for Paper 42's BW-α / BW-γ constructions to lift cleanly** is not a literature question. Need to compute at finite n_max during L2-C/E.

3. **Whether the M3 sub-mechanism trivializes at (3, 1) on chirality-symmetric spectrum** is the Sprint L0 prediction; the verification is computational, not literature-based.

4. **Whether the H_local choice in Paper 42 §7.2 needs to change at signature (3, 1)** is a structural question that Sprint L2-E will answer.

5. **Whether the (3, 1) BBB sign-table convention should be East-coast or West-coast** is a convention choice; the audit (§4.4) says BBB Table 1 gives the same signs for both conventions at (m, n) = (4, 6), but the underlying gamma matrix algebra differs. Pick West-coast (most common in physics) upfront and document.

6. **Whether the twisted-emergence alternative (Devastato–Lizzi–Martinetti 2018 + Nieuviarts 2024–2026) can be applied at the GeoVac AC-extension level** (Sprint H1) is a separate question that could be addressed in a parallel L2-twist-scoping pass. The audit notes: Nieuviarts twist-morphism is not directly applicable to S³ (odd-dim), but could in principle be applied to the spacetime S³ × ℝ (even-dim) or to the inner factor of the H1 AC extension. This is a parallel exploration, not a near-term sprint priority.

7. **Whether the Paper 24 §V Coulomb/HO cross-manifold obstruction (W2b-medium) has a Lorentzian-side analog** is a structural question that may or may not become relevant during L2. Unrelated to L2 closure but worth flagging.

End of memo.
