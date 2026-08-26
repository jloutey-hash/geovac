# Literature verification scout: KMS/thermal ↔ CM (Seam 1) and the conductor-4 arithmetic web (Seam 2)

**Date:** 2026-08-21
**Type:** read-only literature verification. **No papers modified.**
**Guard applied:** every arXiv ID below was fetched and its title/authors/abstract read; every
"absent" claim is a `pdftotext` word-count on the downloaded PDF, not a recollection. Numerical
identities were re-verified locally at 30–40 dps. IDs that did *not* resolve, or resolved to a
different paper than claimed, are called out explicitly.

---

## 0. Bottom line

Two things the project is about to write about are **already theorems, and must be cited**:

* CMR 2005 proves a KMS ↔ CM statement — but about *ground states labelled by CM points of a
  Shimura variety*, not about the modular **flow** being a Hodge circle.
* "CM Hodge structure ⇒ periods are Γ-values" is the Chowla–Selberg / **Gross–Deligne** theorem-
  and-conjecture complex. P56's `prop:hodge_cm_point` is an *instance* of the classical
  "Hodge group of a CM abelian variety is a torus", not new mathematics.

One thing is **not found anywhere** after four differently-phrased search passes:

* the identification of a **Bisognano–Wichmann / thermal modular flow** with the **Hodge circle of
  a Mumford–Tate (norm-1 CM) torus**, with μ₄ torsion at the quarter-periods and the spinor/scalar
  discriminant sitting at the half-period. That is the candidate-novel core of `_kms_mt_torus.py`.

And one **corpus citation defect**: `connes_marcolli2004` (used in P18 §, P56 §) is
arXiv:math/0409306 *"Renormalization and motivic Galois theory"*. It contains **no** KMS states,
**no** Q-lattices, **no** complex multiplication. It cannot carry any Seam-1 weight. The QSM/CM
results need `math/0404128` (GL(2) system) and `math/0501424` (CMR) added as separate keys.

---

## 1. Per-item verdict table

Legend: **THM** = theorem in the literature, cite it. **FOLK** = folklore-adjacent / textbook
analogy, safe to assert with a standard reference but not a quotable theorem. **NF** = not found;
candidate novel (search-negative, not a proof of novelty).

### Seam 1 — KMS/thermal ↔ complex multiplication

| # | Item | Verdict | Reference / evidence |
|:--|:-----|:-------:|:---------------------|
| 1.1 | A QSM system whose **KMS_β (β>1) extremal states are parameterized by invertible K-lattices ≅ 𝔸\*_{K,f}/K\***, partition function ζ_K, idele-class symmetry, values in K^ab, class-field-theory intertwining with Gal(K^ab/K) | **THM** | Connes–Marcolli–Ramachandran, arXiv:math/0501424, **Theorem 5.1** (quoted §3 below). Verified: exists, 29 pp, Selecta Math. |
| 1.2 | The CM system is a **specialization of the GL(2) system at CM points τ ∈ ℍ** ("non-generic ground states") | **THM** | CMR math/0501424, Introduction ¶3 (verbatim in §3 below) |
| 1.3 | GL(2) QSM system: KMS_∞ states settle onto the **Shimura variety Sh(GL₂, ℍ^±)**; arithmetic subalgebra ≈ modular Hecke algebra; Galois group of the **modular field** acts on ground-state values; analyzed for **generic (transcendental j)** | **THM** | Connes–Marcolli, arXiv:math/0404128 (6 Apr 2004), Part I. Verified title/abstract. CMR §3 restates: "In the range β > 2 the set of extremal KMS states is given by the invertible ℚ-lattices, namely by the Shimura variety Sh(GL₂, ℍ^±)" |
| 1.4 | Anything in CMR / the GL(2) papers connecting a **BW or Tomita modular flow** to a **Mumford–Tate torus / Hodge structure** | **NF** | Word counts on the extracted CMR text (1634 lines, 71× "KMS"): **Hodge 0, Mumford 0, Tomita 0, Bisognano 0, "modular flow" 0, "modular automorphism" 0, Gamma 0, Chowla 0** |
| 1.5 | Γ-values / Chowla–Selberg content inside CMR | **NF** | same word counts (Gamma 0, Chowla 0). The Γ-value side of CM lives in Gross–Deligne (item 2.1), not in the QSM literature |
| 1.6 | **Hodge group (special MT group) of a CM elliptic curve = norm-1 torus of the CM field** | **THM (classical)** | Mumford/Shimura; standard in Moonen's *Notes on Mumford–Tate groups*. ⚠ Terminology nit for P56 `prop:hodge_cm_point`: the ℚ-Zariski closure of the **Hodge circle** h(U(1)) is the **Hodge group / special MT group**; the full MT group additionally carries the weight cocharacter, so for a CM elliptic curve MT = Res_{K/ℚ}𝔾_m (dim 2) while Hg = K¹ (dim 1). P56's proof computes the closure of the Hodge circle, i.e. Hg. The conclusion (a 1-dim torus, proper in SL₂, CM field ℚ(i)) is unaffected; the *name* is off by one notion |
| 1.7 | **Connes–Rovelli thermal time** carrying any arithmetic/CM structure in later literature | **NF** | Searched thermal-time ∩ {CM, motivic, arithmetic}. Returns are philosophy-of-physics (Swanson; Chua; Paetz) and the Connes–Marcolli "cooling procedure" for emergent geometry. No arithmetic realization of thermal time found |
| 1.8 | **Wick rotation = passing between real forms of the complexified structure group** (signature = choice of real form) | **THM** | Helleland–Hervik, *Wick rotations and real GIT*, arXiv:1703.04576 (12 Mar 2017), J. Geom. Phys. Verified. Uses real GIT on real forms of the complexified structure group; discusses O(p,q) real forms of O(n,ℂ) and G₂ ↔ split-G₂ |
| 1.9 | The specific **1-dim torus** case: SO(2) (compact, norm-1 torus of an imaginary quadratic field) vs SO(1,1) (split 𝔾_m) as the two real forms | **FOLK** | Textbook algebraic-group fact (the two real forms of 𝔾_{m,ℂ}); widely used in the higher-spin/AdS literature as "the Wick rotation maps 𝔰𝔬(1,1) → 𝔰𝔬(2)". Not singled out as a *torus* statement in 1703.04576 |
| 1.10 | **"Wick rotation = base change to ℚ(i)"** — an *arithmetic* (motivic / field-of-definition) reading of signature change | **NF** | No hit in the motivic or NCG literature. Nearest neighbours are geometric (1703.04576, real forms) or analytic (Kontsevich–Segal allowable metrics; Zhou's contour rotation, item 2.7) — none arithmetic |
| 1.11 | A **physical flow realizing the Hodge circle** (closest prior art) | **FOLK / adjacent** | Angius–Volpato, arXiv:2605.30418 (28 May 2026; v2 10 Jul 2026), *Hodge Loci and Complex Multiplication via Generalized Symmetries in Calabi–Yau sigma models*: "the Hodge decomposition is determined by the **U(1)×U(1) R-charges**", rational structure from BPS boundary states, polarization from the open-string Witten index, CM number fields embedding at special loci. **But**: R-symmetry, not a thermal/modular flow; no KMS, no Tomita, no BW. This is the closest thing found to "a physical circle action *is* the Hodge circle" |
| 1.12 | **BW modular flow = Hodge circle of a ℚ(i) MT torus**; μ₄ torsion at quarter-periods; e^{iπK} = −I for half-integer m_j vs +I for integer m_l ("spin double cover inside the thermal circle") | **NF — candidate novel** | Four independent search framings (BW∩MT; Tomita∩Deligne-torus; modular-flow∩Weil-operator; Deligne-torus∩QFT) all empty. See §4 caveats |
| 1.13 | BW normalization: Δ^{is} = boost, KMS at β = 2π | **THM (textbook)** | Bisognano–Wichmann 1975/76; Sewell 1982. Universally stated as U(s) = e^{2πisK₁}. The project's β = 2π is the standard normalization, not a free choice |

### Seam 2 — the conductor-4 arithmetic web

| # | Item | Verdict | Reference / evidence |
|:--|:-----|:-------:|:---------------------|
| 2.1 | **CM Hodge structure ⇒ periods are Γ-values** (with exponents set by the Hodge decomposition) | **THM (Chowla–Selberg) + conjecture (Gross–Deligne)** | Chowla–Selberg 1949 PNAS / 1967 J. reine angew. Math. **227**, 86–110. Gross, *On the periods of abelian integrals and a formula of Chowla and Selberg*, Invent. Math. **45** (1978) 193–212 (Deligne supplied the general formulation). Modern: Fresán, *Periods of Hodge structures and special values of the gamma function*, arXiv:1403.4105, Invent. Math. (2017) — proves an **alternating variant** for smooth projective varieties with finite-order automorphisms (via Saito–Terasoma; improves Maillot–Rössler) |
| 2.2 | **K(1/√2) = Γ(1/4)²/(4√π)** | **THM (classical)** | First singular value k₁; Legendre/Gauss, and the CM case of Chowla–Selberg. **Verified locally to 40 dps: exact agreement** |
| 2.3 | **ϖ = Γ(1/4)²/(2√(2π)) = π/M(1,√2)** (Gauss, 30 May 1799) | **THM (classical)** | Canonical historical treatment: D. A. Cox, *The arithmetic-geometric mean of Gauss*, L'Enseign. Math. **30** (1984) 275–330. **Verified locally: π/agm(1,√2) = Γ(1/4)²/(2√(2π)) to 40 dps** |
| 2.4 | **G = β(2) = L(2, χ₋₄)** | definitional | χ₋₄ = the odd primitive character mod 4; disc ℚ(i) = −4 |
| 2.5 | **y² = x³ − x has conductor 32**, newform f = η(4τ)²η(8τ)² = q − 2q⁵ − 3q⁹ + 6q¹³ + 2q¹⁷ + …, weight 2 level 32; CM by ℤ[i] (j = 1728) | **THM** | Verbatim from Moerman, arXiv:2008.06749 §1 (verified). CM by ℤ[i] classical (Deuring) |
| 2.6 | **Is L(E,2) of the conductor-32 CM curve known to involve G or ϖ?** | **NO — and this is a clean negative** | (a) Moerman arXiv:2008.06749 (J. Number Theory, DOI 10.1016/j.jnt.2021.09.013): full text scanned — **Catalan 0, lemniscate 0, Γ(1/4) 0**. L(E,2) and L(E,3) are known **as periods** (Zudilin, *Period(d)ness of L-values*, Springer Proc. Math. Stat. **43** (2013) 381–395), L(E,4) as an explicit 4-cube period (Moerman Thm 1); the reported L(E,2) form is a **θ₂/θ₃ theta-integral** — a Γ(2)-level object, i.e. the project's own modular home. (b) Structurally, Bloch proved the Beilinson conjecture for **L(E,2) of CM elliptic curves** — so L(E,2) is a **regulator** (elliptic dilogarithm) times an algebraic factor, not a Γ- or G-expression. New proof via hypergeometrics: Ito, arXiv:1605.01145. ⚠ Beware: some literature writes the conductor-32 curve as y² = x³ + 4x (isogenous/twist); L(E,s) is the same for any conductor-32 curve over ℚ (Moerman states this explicitly) |
| 2.7 | **Catalan G / β(2) in the Bessel-moment literature** | **NOT PRESENT** (four canonical sources, all zero) | Word counts on downloaded full texts: Broadhurst arXiv:1604.03057 *Feynman integrals, L-series and Kloosterman moments* — **Catalan 0**; Fresán–Sabbah–Yu arXiv:2006.02702 — **Catalan 0**; Zhou arXiv:1706.08308 — **Catalan 0**; Bailey–Borwein–Broadhurst–Glasser arXiv:0801.0891 *Elliptic integral evaluations of Bessel moments* — **Catalan 0** (but 18× "singular value", 9× AGM, 14× "modular"). Extraction sanity-checked (43× "Bessel", 5× "Clausen"). Their conductors/levels are 3, 6, 8, 15 …, not 4 |
| 2.8 | **χ₋₄ (conductor 4) in the Broadhurst–Dorigoni resurgent-Lambert framework** | **THM — present, but in the topological-string sector, not the Feynman sector** | Broadhurst–Dorigoni, arXiv:2607.14020, PoS **LL2026** 020 (submitted 15 Jul 2026) — verified. Feynman sunrise/banana use χ₃,₂ and χ₆,₅ on Γ₀(6). **§5.1** (spectral trace of local ℙ^{a,b}, conductor N = a+b+1): at (a,b) = (2,1) "we encounter the odd character χ_{4,3}(n) = ±1 for n = ±1 mod 4", i.e. **χ₋₄**, in the combination L₁(χ₄,₃; τ) − √−4 · L₁(χ₄,₃; −1/4τ), requiring **directional Borel resummation** at both couplings. Catalan G itself is **not named** (word counts: Catalan 0, "complex multiplication" 0, Γ(2) 0). Companion: arXiv:2507.21352 *Resurgent Lambert series with characters* (28 Jul 2025; rev. 30 Jul 2026) — verified |
| 2.9 | **∫₀¹ K(k) dk = 2G** provenance | **FOLK / classical, with a clean citable secondary source** | Elementary: swap the order of integration → ∫₀^{π/2} θ/sin θ dθ. That integral = 2G is classical and is stated as "one of the most important equivalent forms", eq. (4), in **G. Jameson & N. Lord, "Integrals evaluated in terms of Catalan's constant", Math. Gazette 101 (March 2017)**. The moment family lives in Borwein–Borwein–Glasser–Wan, *Moments of Ramanujan's generalized elliptic integrals and extensions of Catalan's constant*, arXiv:1101.1132, JMAA (2011) (e.g. ∫₀¹ arctan(x)K(x) dx/x = G). **Both verified locally: ∫₀¹K dk − 2G = 0 at 40 dps; ∫₀^{π/2}θ csc θ dθ − 2G ≈ −2e−31 at 30 dps.** No single "first" paper — do **not** attribute it to Ramanujan or to any one author |
| 2.10 | **Theta constants are weight-1/2 forms on the metaplectic double cover** | **THM** | Shimura, *On modular forms of half integral weight*, Ann. Math. **97** (1973) 440–481: the automorphy factor for half-integral weight is *defined* via the Jacobi theta function |
| 2.11 | **Metaplectic double cover ↔ spin double cover** ("the Weil rep is the symplectic analogue of the spin rep") | **FOLK (textbook analogy)** | Weil 1964; standard in every treatment. Sharp caveats found: π₁(SO(n)) = ℤ₂ so Spin(n) is the *universal* cover, whereas π₁(Sp(2n,ℝ)) = ℤ so Mp is only *one* of many covers; and Mp(2n,ℝ) has no faithful finite-dimensional representation while Spin(n) does. Safe as an analogy, **not** as an isomorphism claim |
| 2.12 | **Half-integer weight read as a genuine spin structure** (not merely an analogy) | **THM — and this is the sharp citation the project actually wants** | **M. F. Atiyah, "Riemann surfaces and spin structures", Ann. Sci. ÉNS (4) 4 (1971) 47–62**: theta characteristics (square roots of the canonical bundle; 2^{2g} of them in genus g) **are** spin structures on the Riemann surface. This is the theorem-level content behind "half-integer modular weight is a spin double-cover shadow" — the theta *multiplier* is carried by the choice of theta characteristic = spin structure |
| 2.13 | Weil representation ↔ spin **as an arithmetic statement** about GeoVac-style half-integer m_j | **NF** | Nothing found linking a physical half-integer angular-momentum sector to the metaplectic weight-1/2 story |

---

## 2. Corpus-citation audit (the four keys named in the brief)

| Key | Resolves to | Verified content | Does it cover the new seams? |
|:----|:------------|:-----------------|:-----------------------------|
| `connes_marcolli2004` (P18, P56) | arXiv:**math/0409306**, *Renormalization and motivic Galois theory*, IMRN 2004(76) 4073–4091 | Renormalization ↔ motivic Galois group; Tannakian; mixed Tate | **NO.** Zero KMS, zero CM, zero Q-lattices. Cannot support any Seam-1 claim. Add `math/0404128` (GL(2) QSM) and `math/0501424` (CMR) as *new, separate* keys |
| `zhou_wick2018` (P59) | arXiv:**1706.08308**, CNTP **12** (2018) 127–192 | Bessel moments via contour deformation + integration over modular forms; evaluates IKM as explicit constants or **critical values of modular L-series**; verifies Broadhurst conjectures | Partially. ⚠ **"Wick rotation" there means a 90° contour rotation converting I₀K₀-moments into J₀Y₀-moments** (§2.2, §5.1), *not* a modular transformation and *not* arithmetic. Do **not** cite it for "Wick rotation = arithmetic base change". Catalan 0 |
| `broadhurst_dorigoni2026` (P59) | arXiv:**2607.14020**, PoS **LL2026** 020 | Resurgent transseries for Lambert series Σ a(n)q^n/(1−q^n) with a(n) = χ(n)/n^s; Feynman sector on Γ₀(6); topological-string sector at conductor N = a+b+1 | Yes for the resurgence route. **Adds a conductor-4 datapoint the project did not know it had**: χ₋₄ at local ℙ^{2,1} (§5.1) |
| `fresansabbahyu2023` (P59) | arXiv:**2006.02702**, Algebra Number Theory **17** (2023) 541–602 | Quadratic relations between Bessel moments via the period pairing on Sym^k of the Kloosterman connection; Deligne's conjecture made explicit | Yes for the Bessel-moment framework. Catalan 0, Γ(1/4) 0 |

---

## 3. The five most load-bearing exact references

**(R1) Connes–Marcolli–Ramachandran, "KMS states and complex multiplication", arXiv:math/0501424 (24 Jan 2005), Selecta Math.**
*Theorem 5.1* (structure, from the extracted text): for the C\*-dynamical system (𝒜_K, σ_t) built from
commensurability of 1-dimensional K-lattices, K imaginary quadratic, **arbitrary class number**:
(i) unique KMS state for 0 < β ≤ 1; (ii) for β > 1, extremal KMS states are parameterized by
invertible K-lattices, 𝔸\*_{K,f}/K\* ≅ ℰ_β, with a **free and transitive** idele-class-group action;
(iii) explicit Gibbs form φ_{β,L}(f) = ζ_K(β)⁻¹ Σ_{J ideal in 𝒪} f(J⁻¹L, J⁻¹L) n(J)^{−β};
(iv) the extremal KMS states evaluated on the arithmetic subalgebra 𝒜_{K,ℚ} **take values in K^ab**,
and the class-field-theory isomorphism intertwines the 𝔸\*_{K,f}/K\* symmetry with Gal(K^ab/K).
*Introduction, ¶3 (verbatim):* "it is also a specialization of the GL₂-system of [9] to elliptic
curves with complex multiplication by K. In this case the ground states can be related to the
non-generic ground states of the GL₂-system, associated to points τ ∈ ℍ with complex multiplication,
and the group of symmetries is the Galois group of the maximal abelian extension of K."
→ **This is the "KMS ↔ CM" theorem. It is about ground states *labelled by* CM points of a Shimura
variety. It says nothing about the modular flow σ_t itself being a Hodge circle, and contains no
Hodge / Mumford–Tate / Γ-value content whatsoever.**

**(R2) Connes–Marcolli, "From Physics to Number Theory via Noncommutative Geometry, Part I: Quantum Statistical Mechanics of Q-lattices", arXiv:math/0404128 (6 Apr 2004).**
The GL(2) system. Abstract (verbatim): "The system at zero temperature settles onto a classical
Shimura variety, which parameterizes the pure phases of the system… It acts on values of the ground
states at the rational elements via the Galois group of the modular field." Ground states are
analyzed **for the generic case of transcendental j-invariant**; CM points are the *non-generic*
case handed to (R1). For β > 2 the extremal KMS states = Sh(GL₂, ℍ^±).

**(R3) Gross, Invent. Math. 45 (1978) 193–212 (formulation with Deligne) + Fresán, arXiv:1403.4105, Invent. Math. (2017).**
The Gross–Deligne conjecture: *periods of geometric Hodge structures with multiplication by an
abelian number field are products of Γ-values at rational arguments, with exponents determined by
the Hodge decomposition.* Fresán proves an alternating variant for smooth projective varieties with
finite-order automorphisms.
→ **This is the theorem-level home of "GeoVac sits at a CM point, therefore its periods are
Γ-values / Chowla–Selberg". Present the Γ(1/4)² appearance as an *instance* of
Chowla–Selberg / Gross–Deligne, never as an independent discovery.**

**(R4) Moerman, "L-values for conductor 32", arXiv:2008.06749, J. Number Theory (2021), DOI 10.1016/j.jnt.2021.09.013.**
Fixes the curve **E : y² = x³ − x**, conductor 32, newform f(τ) = q Π_m (1−q^{4m})²(1−q^{8m})²
= q − 2q⁵ − 3q⁹ + 6q¹³ + 2q¹⁷ + …, weight 2 level 32, with L(E,s) = L(f,s) valid for **any**
conductor-32 curve over ℚ. Theorem 1 gives L(E,4) as an explicit integral over [0,1]⁴; L(E,2) and
L(E,3) as periods are due to Zudilin (Springer Proc. Math. Stat. 43 (2013) 381–395), the L(E,2)
representation being a **θ₂θ₃ theta-integral**.
→ **Answers the brief's question directly: L(E,2) of the lemniscatic CM curve is known — as a period
/ Beilinson regulator (Bloch's theorem for CM curves) — and is *not* known in terms of G or ϖ. No
Catalan, no lemniscate constant, no Γ(1/4) anywhere in the paper.**

**(R5) Atiyah, "Riemann surfaces and spin structures", Ann. Sci. ÉNS (4) 4 (1971) 47–62.**
Theta characteristics (square roots of the canonical class; 2^{2g} of them in genus g) are precisely
the spin structures. Together with Shimura, Ann. Math. **97** (1973) 440–481 (half-integral weight
defined via the Jacobi theta automorphy factor), this is the citable backbone for reading
half-integer modular weight as spin-structure data — substantially stronger than the
metaplectic-vs-Spin *analogy*, which is textbook but only an analogy.

---

## 4. What is genuinely open / candidate-novel, with honest caveats

**Candidate novel (Seam 1 core).** The identification
`BW modular generator K = diag(2m_j) ⟹ e^{itK} conjugates (in the Kramers frame) to
cos t·I + sin t·J = the Hodge circle of the ℚ(i) norm-1 torus`, together with (i) β = 2π closure ⇔
compactness, (ii) **e^{iπK} = −I only in the half-integer (spinor) sector**, +I in the integer
(scalar) sector — the spin double cover living *inside* the thermal circle — and (iii) the four
quarter-points {I, J, −I, −J} = μ₄ torsion. No prior art found. Caveats that must be stated in any
paper:

1. **Search-negative ≠ novel.** Four framings, all empty, is decent evidence but not exhaustive; the
   Hodge-theory and AQFT literatures use disjoint vocabularies, and a hit could be phrased in either.
2. **The two halves are individually classical.** "Hodge circle of a CM elliptic curve generates the
   norm-1 torus" (item 1.6) and "BW modular flow is the boost, KMS at β = 2π" (item 1.13) are both
   textbook. The claim is the *identification*, and it is a 2×2 computation — a referee will read it
   as an Observation, not a theorem, unless it is shown to do work.
3. **Terminology.** Fix Hodge group / special Mumford–Tate group vs Mumford–Tate group (item 1.6)
   before the claim goes into prose.
4. **Closest prior art exists and should be cited defensively:** Angius–Volpato arXiv:2605.30418
   realizes the Hodge decomposition by a *physical* U(1)×U(1) R-symmetry flow, with CM at special
   loci. Different flow, same shape of idea. Not citing it would look like a literature gap.
5. **The compactness argument is one WH7 already owns.** e^{2πiK} = I holds because K has integer
   spectrum at finite cutoff — the same fact that made the truncated Bisognano–Wichmann boost compact
   and the Lorentzian signature metrically invisible (CLAUDE.md §1.7, WH7). This is a consistency win,
   but it also means the μ₄/Hodge-circle statement is a *re-reading* of an already-registered
   structural fact, not independent evidence for it.

**Candidate novel (Seam 2).** The observation that GeoVac's period content sits at **conductor 4**
(G = β(2) = L(2,χ₋₄), Γ(1/4)², ϖ, K(1/√2), ∫₀¹K = 2G) while the **Bessel-moment / Feynman-period
literature it is adjacent to sits at conductors 3, 6, 8, 15** — with Catalan literally absent from
all four canonical Bessel-moment sources — is, as far as this scout can establish, unremarked. Two
qualifications: (a) the individual conductor-4 identities are all classical, so the novelty is the
*disjointness observation*, not the constants; (b) Broadhurst–Dorigoni arXiv:2607.14020 §5.1 **does**
reach χ₋₄, via the local ℙ^{2,1} topological-string spectral trace — so the honest statement is
"absent from the Feynman/Bessel-moment sector, present in the topological-string sector of the same
resurgence framework", not "absent from the literature".

**Not novel — do not claim.** (i) CM Hodge structure ⇒ Γ-values (Chowla–Selberg / Gross–Deligne).
(ii) Hodge group of a CM elliptic curve is the norm-1 torus. (iii) KMS ground states ↔ CM points
(CMR). (iv) Wick rotation ↔ real forms (Helleland–Hervik). (v) K(1/√2) = Γ(1/4)²/(4√π),
ϖ = π/M(1,√2), ∫₀¹K dk = 2G — all classical. (vi) Theta constants as weight-1/2 metaplectic forms
(Shimura 1973); theta characteristics = spin structures (Atiyah 1971).

---

## 5. Verification log (local recomputation, mpmath)

| Identity | Result |
|:---------|:-------|
| ∫₀¹ K(k) dk vs 2G | 1.831931188354438030109207029864768221548 both; **diff = 0** at 40 dps |
| ∫₀^{π/2} θ/sin θ dθ vs 2G | diff ≈ −2.0e−31 at 30 dps |
| K(1/√2) vs Γ(1/4)²/(4√π) | 1.854074677301371918433850347195260046218 both; exact at 40 dps |
| π/agm(1,√2) vs Γ(1/4)²/(2√(2π)) | 2.622057554292119810464839589891119413683 both; exact at 40 dps (the direct ∫₀¹dt/√(1−t⁴) agrees to 22 dps — endpoint singularity in the quadrature, not a discrepancy) |

---

## 6. arXiv IDs verified in this scout (all resolve to the claimed paper and content)

math/0404128 · math/0409306 · math/0501424 · 0801.0891 · 1005.2941 · 1101.1132 · 1403.4105 ·
1604.03057 · 1605.01145 · 1703.04576 · 1706.08308 · 2006.02702 · 2008.06749 · 2507.21352 ·
2605.30418 · 2607.14020.

Non-arXiv references verified via publisher/repository page: Atiyah, Ann. Sci. ÉNS (4) **4** (1971)
47–62 (Numdam/EuDML) · Shimura, Ann. Math. **97** (1973) 440–481 (Annals) · Gross, Invent. Math.
**45** (1978) 193–212 (Springer/EuDML) · Cox, L'Enseign. Math. **30** (1984) 275–330 ·
Jameson & Lord, Math. Gazette **101** (March 2017) · Zudilin, Springer Proc. Math. Stat. **43**
(2013) 381–395 (via Moerman ref [22]).

**No fabricated IDs encountered.** One misuse found in the existing corpus: `connes_marcolli2004`
is a real paper cited under a real key, but it is the *renormalization* paper and carries none of
the KMS/CM content the new seam needs.
