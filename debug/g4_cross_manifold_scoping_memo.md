# G4 Cross-Manifold Scoping: Tensor-Product Spectral Triple T_{S³} ⊗ T_{S⁵}

**Sprint:** G4 scoping (parallel to G3 sprint), Paper 32 §VIII.B dominant gap.
**Date:** 2026-05-06.
**Author:** PM scoping fork.
**Status:** Scoping only. No closure claimed.

---

## §0. Executive verdict

**FAR**, with dominant obstruction **GEOMETRIC** (Dirac asymmetry).

The cross-manifold tensor product T_{S³} ⊗ T_{S⁵} cannot be made a spectral triple in the canonical Connes / Marcolli–van Suijlekom sense without one of three structural concessions, each of which gives away some content the construction was meant to preserve. The dominant obstruction is the four-layer Coulomb/HO asymmetry of Paper 24 §V, surfacing here as: T_{S³} carries a Camporesi–Higuchi *Riemannian spinor* Dirac with KO-dim 3, while T_{S⁵} as Paper 24 builds it carries a *first-order complex Euler operator* on the Hardy sector — these are different kinds of operators, and the natural NCG tensor product requires both factors to be Riemannian spectral triples in the same sense.

**Splitting the gap.** G4 cleanly factorizes into two sub-gaps with very different reachability:

- **G4a (formal Connes unification on T_{S³} alone)**: take A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) and tensor *only* with the S³ side. Reachable in 1–2 month sprint extending `geovac/almost_commutative.py`. Closes G4 formally if it works, at the cost of making the Sprint ST-SU3 Bargmann Wilson construction structurally redundant for SM purposes.
- **G4b (genuine cross-manifold S³ ⊗ S⁵ unification)**: requires a non-Riemannian spectral-triple framework, or a unifying ambient manifold containing both S³ Coulomb and S⁵ Bargmann as natural sub-structures. Multi-month, and probably requires structural mathematics that does not currently exist in the Connes-NCG literature.

The recommendation is to pursue **G4a** as the cross-manifold sprint after G3 lands, and to record G4b as a structural open question rather than a sprint target.

---

## §1. Construction attempt

### §1.1 The two factors

**T_{S³}** is the GeoVac truncated metric spectral triple from Paper 32 §III, the operator-system flavor R2 of WH1, with:

- **Algebra**: A_{S³} = O_{n_max} = P_{n_max} C^∞(S³) P_{n_max}, the Connes–vS truncated operator system on the Fock-projected S³ graph. *-closed but not multiplicatively closed; propagation number prop = 2 at n_max ∈ {2, 3, 4} (matching Connes–vS Toeplitz S¹ Prop 4.2 verbatim).
- **Hilbert space**: H_{S³} = ⊕_{n ≤ n_max} (n_fock, l, m_l, m_j) basis, with dim_H scaling as the cumulative shell sum 4 × Σ (n+1)(n+2) for the full Dirac sector.
- **Dirac**: D_{S³} = Camporesi–Higuchi Dirac on round S³, |λ_n| = n + 3/2, with degeneracies g_n^Dirac = 2(n+1)(n+2). KO-dimension 3 (J² = −I, JD_{S³} = +D_{S³}J), verified to machine precision in `geovac/real_structure.py`.
- **Real structure**: J_{S³} the standard Riemannian spinor charge conjugation on S³.
- **Status**: WH1 PROVEN as of 2026-05-06. The truncated triple converges to the round-S³ Camporesi–Higuchi spectral triple in the Latrémolière propinquity, with rate constant 4/π = Vol(S²)/π² (the M1 Hopf-base measure signature of the master Mellin engine).

**T_{S⁵}** does not exist as a metric spectral triple in the same sense as T_{S³}. The Bargmann–Segal lattice of Paper 24 is built from the *holomorphic Hardy sector* H²(S⁵) restricted to symmetric SU(3) irreps (N, 0). Its natural diagonal operator is the Euler / number operator N̂ = Σ z_i ∂/∂z_i, which:

- Is **first-order complex-analytic** (Folland–Stein, Kohn Laplacian on H²(S⁵)), not second-order Riemannian.
- Has **linear** spectrum N (in units of ℏω), not the nonlinear stereographic spectrum n² − 1 of S³.
- The Bargmann graph adjacency encodes **dipole transitions** (ΔN = ±1, Δl = ±1), not the spectrum: A is spectroscopically informative but spectrally inert (Paper 24 §III).
- **Has no canonical Riemannian Dirac**. The round-S⁵ Riemannian spinor Dirac (5-dimensional manifold, KO-dim 5) is a different object from the Bargmann Euler operator and uses a different graph (one that would discretize the round-S⁵ Laplace–Beltrami, not the Hardy-sector Bargmann construction).

This is the four-layer Coulomb/HO asymmetry of Paper 24 §V, reformulated at the spectral-triple level. Layer (i) (spectrum-computing role of L_0) and layer (ii) (calibration π) are the surface manifestations. The deeper structural fact is that S⁵ in the Bargmann construction is *complex-analytic* in a way S³ in the Coulomb construction is not — the two sides live in different categories of geometric structure.

### §1.2 Two candidate T_{S⁵}

To form T_{S³} ⊗ T_{S⁵} we must first commit to a T_{S⁵}. There are two natural choices, each with severe drawbacks:

**Option A (Riemannian round-S⁵, abandons Bargmann content).**
Take T_{S⁵} = (C^∞(S⁵), L²(S⁵, S), D_{S⁵}^{spinor}) with D_{S⁵}^{spinor} the standard Riemannian spinor Dirac on round S⁵. This is a textbook Connes spectral triple with KO-dim 5 (J² = −I, JD = −DJ in the conventional Riemannian sign table for d ≡ 5 mod 8).

| Pros | Cons |
|------|------|
| Standard, well-defined, KO-arithmetic clean | Discards the Bargmann (N, 0)-tower structure — the very thing that made Paper 24 a "natural" S⁵ construction |
| Tensor product T_{S³} ⊗ T_{S⁵}^{Option A} is a textbook NCG object | The SU(3) Wilson construction of Sprint ST-SU3 lives on the Bargmann graph, NOT on round-S⁵. Option A throws ST-SU3 away. |
| KO arithmetic: 3 + 5 = 8 ≡ 0 mod 8 (J² = +I, JD = +DJ on the combined object) | The Bargmann graph's bit-exact π-free certificate (Paper 24 §III) is irrelevant to round-S⁵; round-S⁵ has Vol(S⁵) = π³ etc. via standard Weyl asymptotics. Calibration π reappears on the S⁵ side under Option A. |

**Option B (Bargmann-Euler-based, breaks the spectral-triple formalism).**
Take T_{S⁵}^{Option B} with algebra A_{S⁵} = polynomial truncation on Bargmann nodes, Hilbert space H_{S⁵} = ⊕_{N ≤ N_max} (N, 0) symmetric SU(3) irrep states, and "Dirac-like" operator built from a chiral doubling of the Euler operator N̂.

| Pros | Cons |
|------|------|
| Preserves Sprint ST-SU3 SU(3) Wilson content on the Bargmann graph | The Euler operator is first-order complex-analytic, not Riemannian. Promoting it to a "Dirac" requires a chirality grading on the (N, 0) tower that does not naturally exist (the (N, 0) symmetric irreps are real / non-chiral) |
| Stays in the bit-exact π-free regime of Paper 24 | KO-dim is not classically defined (the Hardy-sector setup is *Type II / B-field* in NCG language; see Hawkins 2000 for the closest published analog, on the unit disk rather than S⁵, and even there the formalism is non-standard) |
|  | No published precedent for tensoring a Riemannian spectral triple with a Hardy-sector first-order complex object. The literature stays inside Riemannian (Connes) or stays inside Berezin-Toeplitz (Hawkins, Hekkelman). |

**Net.** Neither candidate preserves both (a) standard NCG axioms and (b) the structural content (Camporesi–Higuchi Dirac on S³, SU(3) Wilson on Bargmann) that Papers 25 / 30 / ST-SU3 established as the GeoVac record on the two manifolds. We can have a clean spectral triple (Option A) at the cost of the Bargmann content, or we can keep the Bargmann content (Option B) at the cost of leaving the spectral-triple formalism as published.

### §1.3 Tensor product KO arithmetic (Option A)

Under Option A we can carry through the standard Connes–Marcolli tensor product (Connes–Marcolli 2008 §13.4):

- T_{S³} has KO-dim 3, ε = −, ε' = +. (J² = −I, JD = +DJ.)
- T_{S⁵}^{Option A} has KO-dim 5, ε = +, ε' = −. (J² = +I, JD = −DJ.)
- **Combined KO-dim = 3 + 5 = 8 ≡ 0 mod 8.** Combined ε = (−)(+) = −1, but the JD relation is (+)(−) = −, so JD = −DJ on the combined.

Wait — KO-dim 0 mod 8 requires J² = +I and JD = +DJ in the Connes–Marcolli sign table. The naive ε product gives J² = −I, which would put us at KO-dim 1 mod 8 instead. The correct way to compute this is to use the proper sign-table multiplication:

| KO-dim | J² | JD = ±DJ | ε | ε' |
|--------|----|----|---|-----|
| 0 | +I | +DJ | + | + |
| 1 | −I | +DJ | − | + |
| 2 | −I | −DJ | − | − |
| 3 | +I | −DJ | + | − |
| 4 | +I | +DJ | + | + |
| 5 | −I | +DJ | − | + |
| 6 | −I | −DJ | − | − |
| 7 | +I | −DJ | + | − |

(Connes–Marcolli 2008 Table 13.1.) But Paper 32 §IV uses Euclidean convention KO-dim 3 with JD = +DJ, ε = −, ε' = +, which differs from the table above (table has KO-3 with JD = −DJ). The Euclidean convention is from Connes–vS 2021 / Hekkelman 2024, and consistent with that convention:

- T_{S³} (Euclidean KO-dim 3): ε = −, ε' = +, JD = +DJ
- T_{S⁵}^{Option A} (Euclidean KO-dim 5): ε = −, ε' = +, JD = +DJ
- **Combined**: ε = (−)(−) = +, ε' = (+)(+) = +, KO-dim 3 + 5 = 8 ≡ 0 mod 8 with J² = +I, JD = +DJ.

This is structurally clean. The combined tensor product T_{S³} ⊗ T_{S⁵}^{Option A} is a well-defined Riemannian spectral triple with KO-dim 0.

**But**: this is the tensor product of round-S³ and round-S⁵, neither of which is the GeoVac-specific construction. It is a textbook NCG object that any two compact symmetric-space spectral triples would produce. The GeoVac specificity is in the *truncations* — the Fock-projected S³ Coulomb graph and the Bargmann-Segal Hardy-sector S⁵ graph — and Option A discards the truncations on the S⁵ side.

### §1.4 What can survive the tensoring

The honest accounting at the construction level is:

- **Option A** preserves the spectral-triple axioms and KO arithmetic but loses the SU(3) Wilson content of Sprint ST-SU3 (the Wilson links live on the Bargmann graph, not on round-S⁵). What we get is a vanilla S³ × S⁵ Riemannian product, which is published mathematics and does not constitute a GeoVac result.
- **Option B** preserves the SU(3) Wilson content but exits the Connes spectral-triple formalism, requiring extension of the NCG framework to non-Riemannian / Hardy-sector cases that does not currently exist as standard mathematics.

Neither candidate gives a "T_{S³} ⊗ T_{S⁵}" that is *both* a clean spectral triple *and* a faithful representation of the GeoVac content on the two sub-manifolds. This is the geometric obstruction at its sharpest.

---

## §2. SM gauge content from inner fluctuations

Suppose we accept Option A and compute the inner-fluctuation gauge content of the AC extension

T_{AC} = T_{S³} ⊗ T_{S⁵}^{Option A} ⊗ T_F

with A_F a finite almost-commutative algebra. There are two natural choices for A_F:

### §2.1 A_F = ℂ ⊕ ℍ (electroweak only, à la Sprint H1)

Sprint H1 already showed this works on T_{S³} alone, recovering U(1) × SU(2) gauge fluctuations. Tensoring with T_{S⁵}^{Option A} adds a second Riemannian factor but does not change the inner-fluctuation algebra: A_{S³} ⊗ A_{S⁵} ⊗ A_F still has the same C ⊕ H finite factor, and inner fluctuations still produce U(1) × SU(2) at the gauge level. The Higgs sector is also unchanged (Sprint H1 verdict: admits Higgs structurally, doesn't autonomously emit Yukawa).

The S⁵ factor is a passive spectator at the gauge level. It contributes additional Hilbert-space dimensions (which in a continuum spectral-action computation would affect Seeley–DeWitt coefficients, hence vacuum energies and cosmological constant), but it does not introduce new gauge groups.

### §2.2 A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) (Connes' SM convention)

This is the Pati–Salam style finite algebra of Connes' Standard Model construction. Adding M_3(ℂ) to A_F produces an SU(3) gauge sector via the inner-automorphism subgroup of M_3(ℂ) — *independently of whether S⁵ is in the construction at all*.

This is structurally significant: **the Connes-SM SU(3) does not need S⁵.** Connes' SM works on a single 4-dimensional spacetime manifold tensored with a finite AC algebra; the SU(3) comes from the M_3(ℂ) factor, not from any internal sub-manifold structure.

If GeoVac adopts this convention, then:

- U(1) on T_{S³} (Paper 25, recovered as Cartan torus reduction of SU(2) per Paper 30 §5.1)
- SU(2) on T_{S³} (Paper 30, gauge piece of inner fluctuation with A_F = ℂ ⊕ ℍ)
- SU(3) from the M_3(ℂ) factor of A_F (NOT from S⁵)

This recovers the three SM gauge groups on a single manifold (S³), with all three coming from inner fluctuations of D_{S³} on different finite AC factors. **Sprint ST-SU3's SU(3) Wilson construction on the Bargmann S⁵ graph becomes structurally redundant for SM-unification purposes** in this reading: the M_3(ℂ) finite factor handles SU(3) without needing the Bargmann graph.

### §2.3 What Sprint ST-SU3 was and was not

The ST-SU3 result is the SU(3) Wilson lattice gauge theory on the Bargmann S⁵ graph, with verified Cartan reduction and 1/12 = 1/(4·N_c) kinetic coefficient. This is a real and structural result about **what kind of non-abelian gauge theories the GeoVac graph machinery supports**. It is *not* a unification claim — Sprint ST-SU3 §5.4 already noted that matter coupling rediscovers the Sprint 5 CG obstruction, so ST-SU3's SU(3) is "pure gauge" without natural matter.

For SM unification, Connes' M_3(ℂ) convention is structurally cleaner: the SU(3) comes with a built-in prescription for color-triplet quark matter (the M_3(ℂ) algebra acts on a 3-dim color factor in H_F), and there is no CG-intertwiner obstruction because the matter content is in a fixed irrep across all shells.

### §2.4 The honest gauge-content conclusion

Under Connes' M_3(ℂ) convention plus Option A on the S⁵ side:

- The cross-manifold tensor product **does** produce U(1) × SU(2) × SU(3) gauge content at the inner-fluctuation level.
- The construction is structurally **standard** Connes–Chamseddine SM extended to a higher-dimensional base manifold (S³ × S⁵ instead of M^4).
- ST-SU3's Wilson SU(3) on the Bargmann graph and the Connes-SM SU(3) on M_3(ℂ) are structurally **distinct** SU(3) homes for the same Lie group. The construction picks one or the other.
- The Higgs sector remains in the Sprint H1 status (admits structurally, no autonomous Y selection). Adding the M_3(ℂ) factor to A_F gives a richer Higgs structure (off-diagonal blocks between M_3(ℂ), ℂ, ℍ), but Yukawa entries are still free inputs.

This is **G4a** — the formal Connes unification — and it is reachable as a 1–2 month sprint extending `geovac/almost_commutative.py`. It closes G4 formally **if** we accept that the SU(3) lives on the finite M_3(ℂ) factor rather than on the Bargmann graph.

---

## §3. Where the cross-manifold obstruction sits

### §3.1 Three layers of obstruction

The cross-manifold gap has at least three layers, ordered by increasing depth:

**Layer 3 (least fundamental): matter-coupling rediscovery.** Sprint ST-SU3 §5.2 showed that coupling Wilson SU(3) on Bargmann to natural matter (the (N, 0) shells) requires Clebsch–Gordan intertwiners between irreps of growing dimension — exactly Sprint 5 Track S5's obstruction. This persists at the tensor-product level: if we use Option B for T_{S⁵}, the matter coupling on the S⁵ side stays in (N, 0)-rep CG-intertwiner mode regardless of what we do on the S³ side.

**Layer 2 (medium-fundamental): algebraic mismatch between operator systems.** A_{S³} (Connes–vS truncated operator system, prop = 2) is *-closed but not multiplicatively closed. The natural Bargmann-side algebra A_{S⁵} (polynomial truncation in the (N, 0) basis) has its own multiplicative structure inherited from the Hardy-space polynomial ring, but the two algebras' multiplicative-closure properties are not parallel: Connes–vS truncations break multiplicative closure by index cutoff, while Bargmann polynomial truncation respects multiplicative closure (polynomials of degree ≤ N times polynomials of degree ≤ M give polynomials of degree ≤ N + M, which a hard cutoff may exclude). The two truncations are structurally different categories of object.

**Layer 1 (most fundamental, the dominant obstruction): the Dirac asymmetry.** T_{S³} carries a Riemannian spinor Dirac (D_{S³} = Camporesi–Higuchi, KO-dim 3, square is the spinor Laplacian on S³). T_{S⁵} as Paper 24 builds it carries a *first-order complex Euler operator* (Hardy-sector setup), not a Riemannian Dirac. The standard NCG tensor product requires both factors to be in the same category — both Riemannian (Connes), or both Berezin–Toeplitz (Hawkins / Hekkelman), or both Type-II / B-field (some operator-algebraic generalizations). Mixing categories has no published prescription.

### §3.2 Q3 answer: where does the obstruction sit?

**The dominant obstruction is GEOMETRIC (Layer 1).** The matter-coupling obstruction (Layer 3) and the algebraic mismatch (Layer 2) are surface manifestations of the deeper geometric asymmetry. If we resolve Layer 1 by going to Option A (round-S⁵ Riemannian Dirac), Layers 2 and 3 partially heal: A_{S⁵}^{Option A} = round-S⁵ truncated operator system can be made parallel to A_{S³}, and the matter-coupling problem is replaced by a different choice (matter content on round-S⁵ is determined by the spinor bundle, not by Bargmann (N, 0) shells).

If we resolve Layer 1 by going to Option B, Layer 2 becomes worse (Hardy-sector polynomial algebra vs Connes–vS truncated operator system are categorically different, and there is no published bridging) and Layer 3 stays unchanged (CG intertwiners on (N, 0) tower).

### §3.3 The deeper structural reading

The four-layer Coulomb/HO asymmetry of Paper 24 §V identifies what kind of geometric object S³ Coulomb is and what kind S⁵ Bargmann is. They are in different categories:

- **S³ Coulomb** is a *Riemannian discretization*: the round-S³ Laplace–Beltrami discretizes onto the Fock-projected graph; the spinor Dirac on round-S³ discretizes onto the Camporesi–Higuchi Dirac on the same graph. Both are second-order / first-order Riemannian operators. Paper 32 establishes these as a metric spectral triple.
- **S⁵ Bargmann** is a *Hardy-space discretization*: the holomorphic sector H²(S⁵) restricts to (N, 0) symmetric irreps, the Euler operator gives a linear spectrum, the dipole adjacency encodes complex-analytic transition amplitudes. The natural operator-algebraic framework for this is Berezin–Toeplitz / Hawkins, not Connes' Riemannian spectral triples.

These are two *different* discretization paradigms in NCG. The cross-manifold gap is not just "two manifolds, one construction" — it is **two structurally different categories of geometric object**, joined only by the universal Wilson–Hodge combinatorial vocabulary (Paper 24 §V, layer iv).

The headline answer to Q3: **the dominant obstruction is that we are trying to tensor two spectral triples that live in different NCG categories.** Connes' axioms apply to Riemannian; Hawkins / Berezin–Toeplitz apply to complex-analytic; the cross-category tensor product is not a published construction. This is genuinely outside the existing NCG literature, not just outside the GeoVac codebase.

### §3.4 Cross-validation: Paper 24 §V (the four-layer asymmetry) re-read

Re-reading Paper 24 §V with this scoping in mind, the four layers map onto the cross-manifold gap as follows:

| Layer | Coulomb side (S³) | HO side (S⁵) | Cross-manifold consequence |
|-------|-------------------|--------------|----------------------------|
| (i) Spectrum-computing role of L_0 | Yes (D − A on graph = Rydberg) | No (HO spectrum is diagonal) | Tensor-product Dirac D_{S³} ⊗ I + γ ⊗ D_{S⁵} requires both factors to *be* Dirac operators; the S⁵ side's natural operator (Euler) does not square to a Riemannian Laplacian |
| (ii) Calibration π | Yes (Vol(S²)/4 = π Hopf measure) | No (bit-exact rational) | Tensor product would need a unified calibration; either Option A reintroduces π on round-S⁵ side (loses Paper 24 π-free), or Option B keeps π-free on S⁵ but breaks the spectral triple |
| (iii) Wilson gauge with natural matter coupling | Yes (Paper 30 SU(2) on S³ = SU(2)) | Pure-gauge only (ST-SU3, CG obstruction on matter) | Matter coupling to SU(3) Wilson on Bargmann fails uniformly — adding a tensor factor doesn't fix the (N, 0)-rep CG structure |
| (iv) Universal Wilson–Hodge vocabulary | Yes | Yes | The combinatorial vocabulary IS unified, but it doesn't carry NCG content by itself |

Layer (i) is the critical obstruction at the spectral-triple level; layer (ii) is the transcendental signature of (i); layers (iii)–(iv) are downstream consequences. The four-layer asymmetry of Paper 24 is therefore **identical to the cross-manifold gap of Paper 32 §VIII.B G4**, just expressed at different levels of structure.

---

## §4. Concrete computational test for closing G4 in 1–2 months

### §4.1 The G4a sprint: formal Connes unification on T_{S³} alone

The cleanest sprint-scale test that could close G4a (the formal-Connes unification flavor) in 1–2 months is:

> **Extend `geovac/almost_commutative.py` to A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) on T_{S³} alone, compute the gauge piece of the inner fluctuation, and verify that it produces U(1) × SU(2) × SU(3) at the structural level.**

**Implementation outline:**

1. **Algebra extension.** Replace the current `ElectroweakFiniteTriple` (A_F = ℂ ⊕ ℍ, dim_C(H_F) = 8 doubled) with a Pati–Salam-style `SMFiniteTriple` that has A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ). The Hilbert space H_F doubles to dim_C(H_F^{matter}) = 8 (lepton + quark sectors with appropriate color triplet), then matter–antimatter doubles to 16. KO-dim of the finite factor is 6 (Connes' SM convention).

2. **Combined KO-dim arithmetic.** With T_{S³} at KO-dim 3 (Euclidean) and T_F at KO-dim 6, the combined is 3 + 6 = 9 ≡ 1 mod 8, with J² = −I, JD = +DJ. This matches the Sprint H1 §1.2 computation for the C ⊕ H case; the M_3(ℂ) extension doesn't change the KO arithmetic, only the dimension of H_F.

3. **Gauge inner fluctuation.** For ω = Σ a_i [D, b_i] with a_i, b_i ∈ A_GV ⊗ A_F, the gauge piece is

   ω_gauge = Σ a_i^{GV} [D_{GV}, b_i^{GV}] ⊗ a_i^F b_i^F

   The product structure of A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) splits this into three independent gauge blocks (one per summand of A_F), giving U(1) × SU(2) × SU(3) at the inner-automorphism level. Each block can be computed independently using the same machinery as Sprint H1.

4. **Verification targets.**
   - The U(1) block matches the Cartan torus reduction of the SU(2) block (Paper 30 §5.1).
   - The SU(2) block reproduces Paper 30 / Sprint H1 results.
   - The SU(3) block produces a new SU(3) sector that is structurally **distinct** from the Sprint ST-SU3 Bargmann Wilson construction (this distinction is itself the diagnostic).
   - All three Higgs blocks (off-diagonal between summands of A_F) are non-trivial under any non-zero Yukawa, but no GeoVac structure selects the Yukawa (Sprint H1 verdict generalizes).

5. **Falsifier.** The strong falsifier (analog of Sprint H1 §3): show that for every A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) Yukawa structure derivable from GeoVac data alone, the inner fluctuation produces only gauge-1-forms and never a Higgs. This is expected to fail (Higgs sector exists for any non-zero Y), confirming the positive-thin verdict at the SM scale.

**Estimated effort.** 3–6 weeks for the algebra extension and inner-fluctuation calculation; 2–4 weeks for the Higgs falsifier sweep at n_max ∈ {2, 3} on the larger H_F. Total 1.5 months.

**Outcome (predicted).** G4a closes positive-thin in the same sense as Sprint H1: the Connes-SM construction works on top of the GeoVac S³ triple, producing U(1) × SU(2) × SU(3) gauge content; no GeoVac mechanism selects the SM-distinguishing data (Yukawas, hypercharge assignments, generation count). This would be a valuable **negative on autonomous SM emergence** more than a positive on SM unification.

### §4.2 The G4b open question (out of sprint scope)

For the genuine cross-manifold question (T_{S³} ⊗ T_{S⁵} faithful to both factors), the missing structural ingredient is:

> An NCG framework that allows tensoring a Riemannian spectral triple with a Hardy-sector / first-order complex spectral triple, preserving the structural content of both (KO arithmetic, π-free certification on the complex side, calibration π on the Riemannian side).

This is **not** GeoVac-specific. It is a question about NCG itself. As far as I can tell from the published literature (Connes & van Suijlekom 2021 / 2024; Hekkelman 2022 / 2024; Hekkelman–McDonald 2024; Latrémolière propinquity series; Hawkins 2000; Berezin–Toeplitz literature):

- Connes-style spectral triples on round odd-dim spheres (S^{2n+1}) are well-understood for the *Riemannian* side, but the literature has not addressed the Bargmann–Segal Hardy-sector reformulation as a spectral triple.
- The Hawkins / Berezin–Toeplitz framework on S² and on the unit disk has been developed but not extended to S⁵.
- No paper I have located addresses cross-category tensor products.

Hence: G4b would require either contributing a structural mathematical result to NCG (new framework, multi-month or multi-year), or accepting Option A and doing a vanilla tensor product that loses the Bargmann content.

Recommendation: **Do not open G4b as a sprint target.** Record it as a structural open question in the paper record (Paper 32 §VIII.B G4 entry already does this). If GeoVac develops a research direction that needs cross-category tensor products (e.g., a unified S³ Coulomb + S⁵ Bargmann Hamiltonian for a specific physical system), that will be the motivation to revisit; without such a motivator, G4b is open-ended in the way the α-derivation program is.

---

## §5. Recommendation

### §5.1 Verdict

| Metric | Value |
|--------|-------|
| Reachability | **G4a: 1–2 month sprint, FAR for G4b (multi-month / multi-year, missing NCG framework)** |
| Dominant obstruction | **GEOMETRIC** — Dirac asymmetry between the Riemannian S³ Camporesi–Higuchi triple and the Hardy-sector S⁵ Bargmann construction. Surfaces as the four-layer Coulomb/HO asymmetry of Paper 24 §V at the spectral-triple level. |
| Sealed? | **No**, but probably *categorically resistant*. The obstruction is at the level of NCG framework choice, not at the level of computational hardness. |

### §5.2 Sprint-timing recommendation

**One-sentence recommendation:** G4 should be split into G4a (formal Connes unification with M_3(ℂ) on T_{S³} alone, reachable in 1–2 month sprint after G3 closes) and G4b (genuine cross-manifold S³ ⊗ S⁵ unification, not currently sprint-scale, requires structural NCG mathematics that does not exist in published form).

The G4a sprint, if pursued after G3, would close the formal SM-unification question at the same positive-thin level as Sprint H1: the construction works, but no GeoVac mechanism autonomously selects the SM-distinguishing data. This is a genuine and useful structural result — it would establish that GeoVac sits inside the Marcolli–vS / Connes-SM lineage at the level of "all three SM gauge groups arise as inner fluctuations on a single S³ spectral triple", with the honest caveat that the GeoVac autonomy claim does not extend beyond the gauge-group structure.

The G4b question stays as a structural open question. The obstruction is well-named (Layer 1 in §3.1, the Riemannian-vs-Hardy-sector asymmetry), and the missing ingredient (cross-category NCG tensor product) is also well-named. Whether to pursue this depends on a research direction that motivates it, which is currently absent.

### §5.3 What would change this recommendation

The G4b verdict would shift if any of the following emerged:

1. **A new NCG construction in the literature** that bridges Connes' Riemannian framework with Berezin–Toeplitz / Hardy-sector triples in a tensor-product compatible way. Worth periodically scanning the NCG literature for this; arXiv listings under math.OA and math-ph.
2. **A unifying ambient manifold** containing both S³ Coulomb and S⁵ Bargmann as natural sub-structures with a common operator-algebraic description. Speculative candidates: U(2) (= S³ × S¹ as a real manifold; not obviously containing S⁵), or some homogeneous space of a larger Lie group with both S³ and S⁵ as orbits. No obvious candidate from GeoVac's working palette.
3. **A "Hardy-sector spectral triple" framework** for the Bargmann construction itself — Paper 24 already names the holomorphic Hardy space H²(S⁵) as the natural setting; if someone (us or others) developed this as a Type-II / B-field spectral triple in the Hawkins lineage, the cross-category tensor product question could be reformulated.

### §5.4 Net structural reading

The cross-manifold gap is not a computational gap. **It is a categorical gap in the NCG framework itself**, surfaced by the GeoVac construction's choice to discretize Coulomb (S³ Riemannian) and HO (S⁵ Hardy) in their respective natural categories. Resolving it requires either (a) flattening the asymmetry by going Riemannian on both sides (Option A; loses the bit-exact π-free certificate of Paper 24), (b) developing new structural NCG (multi-year), or (c) accepting that GeoVac as currently constituted produces *two* spectral-triple objects that are not naturally tensorable in the canonical sense.

The GeoVac record can name this gap honestly without claiming closure. Paper 32 §VIII.B G4 already does this; no edits to the paper body are recommended from this scoping pass.

---

## §6. Files referenced

- `papers/synthesis/paper_32_spectral_triple.tex` (esp. §VIII.B G4 entry)
- `papers/observations/paper_30_su2_wilson.tex`
- `papers/synthesis/paper_25_hopf_gauge_structure.tex`
- `papers/core/paper_24_bargmann_segal.tex` (esp. §V Coulomb/HO asymmetry)
- `papers/core/paper_31_universal_coulomb_partition.tex`
- `geovac/almost_commutative.py` (Sprint H1 single-manifold AC extension; the natural extension target for G4a)
- `geovac/operator_system.py` (WH1 R2 truncated operator system, A_{S³} = O_{n_max})
- `geovac/full_dirac_operator_system.py` (R3.5 full Dirac sector on T_{S³})
- `geovac/real_structure.py` (J on truthful CH, used by H1)
- `geovac/su3_wilson_s5.py` (Sprint ST-SU3 Wilson SU(3) on Bargmann)
- `geovac/nuclear/bargmann_graph.py` (Paper 24 implementation)
- `debug/h1_ac_extension_memo.md`
- `debug/st_su3_wilson_memo.md`
- `debug/sm_unified_gauge_synthesis_memo.md`
- `debug/wh1_round1_connes_vs_pdf_verification.md`

External references: Connes & van Suijlekom (2021, CMP 383); Marcolli–van Suijlekom (2014, J. Geom. Phys. 75, arXiv:1301.3480); Perez-Sanchez (2024, arXiv:2401.03705 + 2025 arXiv:2508.17338); Connes–Marcolli (2008, Noncommutative Geometry, Quantum Fields and Motives); Hawkins (2000, J. Reine Angew. Math.); Folland & Stein (1974, CPAM 27, 429).

---

**End of memo.** Word count: ~5,100 words.
