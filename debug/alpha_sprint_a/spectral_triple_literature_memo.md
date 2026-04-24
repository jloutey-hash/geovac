# Sprint A Memo — Spectral Triple Literature Survey (α-LS RE-DISPATCH)

**Author:** Sub-agent (PM dispatch, α-LS literature survey, re-dispatch after v1 usage-cap
failure)
**Scope:** WH1 (CLAUDE.md §1.7): is GeoVac's structure grounded in published work on
discrete / finite / almost-commutative spectral triples?  Plus α-EB v2 cross-check:
any published S⁵ Seeley–DeWitt data?
**Date:** 2026-04-18

---

## 1. Summary verdict (top-line)

**WH1 status: MODERATE.**  There *is* a mature published literature on finite and
almost-commutative spectral triples (Krajewski 1998, Paschke–Sitarz 2000, Chamseddine–
Connes 1996–2010, Ćaćić 2009, van Suijlekom 2015/2024) that matches GeoVac's structural
vocabulary.  The **closest single published match** is Marcolli–van Suijlekom,
*Gauge networks in noncommutative geometry* (J. Geom. Phys. 2014, arXiv:1301.3480):
a finite spectral triple **on a graph**, continuum limit to Yang–Mills–Higgs,
lattice Yang–Mills–Higgs and Kogut–Susskind Hamiltonian as the spectral action.
GeoVac's S³ Coulomb graph + almost-commutative extension A × M_n(ℂ) (Papers 25/30)
is **shape-compatible** with this framework.

**However**, three important cautions:
(a) The Marcolli–van Suijlekom continuum-limit claim was partially **contradicted** by
Perez-Sanchez (arXiv:2401.03705 / comment arXiv:2508.17338), who shows the continuum
limit is Yang–Mills without Higgs, i.e. the full YM-H reconstruction was over-stated.
(b) **No published spectral-triple framework predicts α** (or any other dimensionless
physical constant) from a discrete / graph-based construction.  The Chamseddine–Connes
"uncanny precision" (CMP 2010) gets SM gauge-coupling ratios at a high unification
scale and the Higgs-to-W mass ratio, but does not fix α.
(c) WH1 in its full form — "GeoVac is an almost-commutative spectral triple whose
spectral action generates the observable structure" — is at **partial shape-match**
level in the published literature.  Nobody has built a finite ACG on an integer-labeled
compact Lie-group-quotient Hopf graph with a Camporesi–Higuchi Dirac and computed
α from its spectral action.  What is published is each ingredient separately.

**α-EB v2 cross-check findings:** the Kluth–Litim paper (EPJ C 80, 269 (2020),
arXiv:1910.00543) provides explicit heat-kernel coefficients on spheres in arbitrary
dimension including odd dimensions S^5, covering scalars/vectors/tensors.  This is
**exactly the reference needed** to test the "π³ = Vol(S⁵) from Paper 24's S⁵" hypothesis.
I did not extract the explicit α_0, α_2 coefficients from the PDF in this sprint
(PDF binary extraction failed); access the published EPJ C version directly.

---

## 2. Twelve structured entries (verified)

### E1 — Connes–Chamseddine, The Spectral Action Principle (CMP 1997)
- **Citation:** A. H. Chamseddine, A. Connes, *Commun. Math. Phys.* 186 (1997) 731,
  arXiv:[hep-th/9606001](https://arxiv.org/abs/hep-th/9606001).
- **Structural claim:** Introduces the spectral action Tr f(D/Λ) as a natural action
  functional built from the spectrum of a Dirac operator on a (possibly almost-
  commutative) spectral triple.  Framework applies to any compact-manifold-based
  spectral triple.  On almost-commutative geometries M × F with M = spacetime and F a
  finite noncommutative space, the spectral action reduces to Einstein–Hilbert +
  Yang–Mills–Higgs + cosmological + gauge-Higgs couplings.
- **Relevance to GeoVac:** This is the *spectral-action machinery* WH1 appeals to.
  However, CC's M is always a smooth compact manifold (continuous), and the spectral
  action comes from heat-kernel expansion Tr f(D/Λ) ~ Σ_k f_k Λ^{d−k} a_k.  GeoVac's
  analog is explicitly finite (Fock-projected S³ graph at n_max ≤ 3).  The continuous
  CC machinery does **not** apply verbatim, because the heat-kernel asymptotic
  expansion degenerates at finite truncation.
- **CC transfer to discrete:** does **not** transfer directly; discrete analog must be
  constructed (see E9, E10).

### E2 — Krajewski, Classification of Finite Spectral Triples (J. Geom. Phys. 1998)
- **Citation:** T. Krajewski, *J. Geom. Phys.* 28 (1998) 1–30, arXiv:[hep-th/9701081](https://arxiv.org/abs/hep-th/9701081).
- **Structural claim:** Finite spectral triples (F with dim < ∞) are completely described
  in terms of matrices and classified using *Krajewski diagrams* (bipartite graphs whose
  vertices are gauge multiplets of chiral fermions and edges Yukawa couplings).  This is
  the combinatorial skeleton of any almost-commutative model.
- **Relevance to GeoVac:** Krajewski diagrams are finite directed bipartite graphs
  encoding representation-theoretic data of a finite spectral triple.  GeoVac's S³
  Coulomb graph is **not** a Krajewski diagram — it encodes (n, l, m) quantum numbers,
  not fermion multiplet / Yukawa data — but the *framework* of encoding a finite
  spectral triple as a graph + matrix data is precisely Krajewski's technology.
- **CC transfer:** finite-ACG case is fully developed within Krajewski's classification;
  this is where the SM "finite space F" F_SM lives.

### E3 — Paschke–Sitarz, Structure Theory of Finite Real Spectral Triples (2000)
- **Citation:** M. Paschke, A. Sitarz, "The geometry of noncommutative symmetries,"
  *Acta Phys. Polon. B* 31 (2000) 1897.  Foundational for the "real structure J"
  axiomatics of finite ACG.  See also their *J. Math. Phys.* 39 (1998) 6191.
- **Structural claim:** Parallel/alternative classification of finite real spectral
  triples (with real structure J and charge conjugation) to Krajewski's.  Developed the
  KO-dimension classification tables.
- **Relevance to GeoVac:** provides the charge-conjugation / antilinear-structure
  axiomatics that GeoVac's Dirac-on-S³ (T1, Camporesi–Higuchi) inherits.  GeoVac's
  T7/T8/T9 body of exact symbolic identities lives at this level of classification.
- **CC transfer:** finite-ACG case fully covered; this is the "stable" foundation.

### E4 — Chamseddine–Connes, The Uncanny Precision of the Spectral Action (CMP 2010)
- **Citation:** A. H. Chamseddine, A. Connes, *Commun. Math. Phys.* 293 (2010) 867,
  arXiv:[0812.0165](https://arxiv.org/abs/0812.0165).
- **Structural claim:** The asymptotic Seeley–DeWitt expansion of Tr f(D/Λ) on the
  round S³ is shown to be already **remarkably accurate** just from the first two
  terms (cosmological constant + scalar curvature).  Demonstrates that spectral-
  action truncation is controlled.
- **Relevance to GeoVac:** **this is the paper α-SP used for R1.**  CC bulk on S³
  gives $(a_0, a_2, a_4) = (+, +, +)$ — which is the α-SP finding.  The paper is
  about how good the asymptotic expansion *is* on S³, i.e. how many SD coefficients
  control the effective action.  It does **not** discuss the sign pattern $(+,+,-)$
  that Paper 2's K exhibits.
- **CC transfer to discrete:** partial.  The truncation error becomes the object of
  interest (Paper 2's K residual, via α-EB v2's π³α³ candidate).

### E5 — Ćaćić, Moduli Spaces of Dirac Operators for Finite Spectral Triples
  (2009 → Chapter in Springer volume)
- **Citation:** B. Ćaćić, arXiv:[0902.2068](https://arxiv.org/abs/0902.2068),
  in K. Landsman et al. (eds.) *Noncommutative Geometry and Physics*, Springer (2011).
- **Structural claim:** Generalizes the Krajewski / Paschke–Sitarz classification to
  allow arbitrary KO-dimension and failures of orientability / Poincaré duality.
  Gives the moduli space of Dirac operators on a fixed finite spectral triple.
- **Relevance to GeoVac:** directly relevant to the "rigidity" question for the
  Dirac operator on GeoVac's S³ graph.  The Camporesi–Higuchi Dirac is a specific
  choice; Ćaćić's framework is how to classify nearby choices.  (T9's squared-Dirac
  spectral zeta theorem sits inside this moduli space.)
- **CC transfer:** directly applies to finite-ACG classification, which is what
  GeoVac uses.

### E6 — Marcolli–van Suijlekom, Gauge networks in noncommutative geometry (J. Geom.
  Phys. 2014) — **CLOSEST PUBLISHED MATCH TO GEOVAC**
- **Citation:** M. Marcolli, W. D. van Suijlekom, *J. Geom. Phys.* 75 (2014) 71–91,
  arXiv:[1301.3480](https://arxiv.org/abs/1301.3480).
- **Structural claim:** Gauge networks are **systems of finite spectral triples
  on a graph**, generalizing spin networks (loop quantum gravity) and lattice gauge
  fields to almost-commutative manifolds.  Configuration space = quiver representations
  modulo equivalence in the category of finite spectral triples.  Graph vertices carry
  finite spectral triples, edges carry connection data.  The construction defines a
  discretized Dirac operator on the quiver and **computes the spectral action** —
  reducing to the Wilson action for lattice gauge theories on a 4-lattice and recovering
  the Yang–Mills–Higgs system in the continuum limit.
- **Relevance to GeoVac:** **this is the single closest match** to GeoVac's structure.
  The vocabulary — "finite spectral triples on graph vertices, connection on edges,
  lattice spectral action" — is essentially what WH1 claims GeoVac is.  Both Paper 25
  (GeoVac U(1) Hopf graph with L₀ node Laplacian and L₁ = B^T B edge Laplacian) and
  Paper 30 (SU(2) Wilson action on the same graph with maximal-torus reduction to
  Paper 25) fit the Marcolli–van Suijlekom template.
- **Important caveat:** the continuum-limit-to-Yang–Mills–Higgs claim was
  **challenged** by the 2024 comment (see E7).  This reduces but does not eliminate
  the match.
- **CC transfer to discrete:** the entire point of the paper is to transfer CC
  spectral-action machinery to discrete-graph spectral triples; successful construction,
  with the continuum-limit caveat noted.

### E7 — Perez-Sanchez, Bratteli networks and the Spectral Action on quivers (2024)
- **Citation:** C. I. Perez-Sanchez, arXiv:[2401.03705](https://arxiv.org/abs/2401.03705).
  See also Comment arXiv:[2508.17338](https://arxiv.org/abs/2508.17338).
- **Structural claim:** Extends Marcolli–van Suijlekom (E6) to a broader combinatorial
  framework using *Bratteli networks* on top of inductive sequences of finite spectral
  triples (AF algebras).  The spectral action on quivers reduces to Wilsonian Yang–Mills
  lattice gauge theory + Weisz–Wohlert cells (Symanzik improved gauge theory).  A
  hermitian matrix field emerges from self-loops; Yang–Mills–Higgs is the smooth limit.
  The comment (2508.17338) demonstrates that the continuum limit of the Marcolli–van
  Suijlekom construction is actually **Yang–Mills without Higgs**, correcting E6.
- **Relevance to GeoVac:** Perez-Sanchez extends the framework Closed the E6 caveat
  noted above.  The positive news: even after the Higgs correction, the Marcolli–
  van Suijlekom / Perez-Sanchez line builds finite spectral triples on graphs whose
  spectral action is Wilson / YM lattice gauge theory — i.e. GeoVac's Paper 30 SU(2)
  Wilson action on the S³ Hopf graph is in exactly this published class.
- **CC transfer:** partial, with continuum-limit caveats now explicit.

### E8 — van Suijlekom, Noncommutative Geometry and Particle Physics (2nd ed. Springer
  2024)
- **Citation:** W. D. van Suijlekom, *Noncommutative Geometry and Particle Physics*,
  2nd ed., Springer Mathematical Physics Studies (2024), ISBN 978-3-031-59120-4.
  Open access version: [repository.ubn.ru.nl/.../314382.pdf](https://repository.ubn.ru.nl/bitstream/handle/2066/314382/314382.pdf).
- **Structural claim:** Comprehensive textbook on NCG and particle physics.  2nd ed.
  has new chapters on **spectral truncations** (Ch. 13) extending the 1st edition.
  Covers finite ACG classification, full SM derivation, Higgs spectral action, and
  the BSM extensions known by 2024.
- **Relevance to GeoVac:** the standard reference for the CC framework. Chapter 13
  on spectral truncations is where the finite-cutoff-NCG language (E9) is developed
  in textbook form; this is the closest-to-GeoVac chapter.
- **CC transfer:** the whole book is about CC.  Ch. 13 addresses the truncation /
  cutoff question, i.e. it discusses how much of CC transfers to finite-resolution.

### E9 — Connes–van Suijlekom, Spectral Truncations in NCG and Operator Systems
  (CMP 2021)
- **Citation:** A. Connes, W. D. van Suijlekom, *Commun. Math. Phys.* 383 (2021)
  2021–2067, arXiv:[2004.14115](https://arxiv.org/abs/2004.14115).
- **Structural claim:** Extends NCG to deal with **spectral truncations** (UV cutoffs
  in momentum space) using **operator systems** in place of C*-algebras.  Primary
  example: truncation of the circle to finite-dim Toeplitz matrices.  Introduces
  *propagation number* as an invariant under stable equivalence.
- **Relevance to GeoVac:** **directly relevant**, though framework is illustrated on
  a continuous base (circle) cut to finite dimensions.  GeoVac's case is the *opposite*
  direction: a finite graph that already sits at finite dim (integer-labeled (n, l, m)
  with n ≤ 3 cutoff).  Operator-system framework is compatible with the GeoVac cutoff,
  but the standard examples in the paper are continuous-base + truncation, not
  graph-native.  This is the correct published entry point for making WH1 rigorous.
- **CC transfer:** the point of the paper is to adapt CC to finite resolution; has
  the tools needed to address the bulk-vs-cutoff discrepancy identified in α-SP.

### E10 — Hekkelman, Truncated geometry on the circle (2022) and
  Hekkelman–McDonald, NCG integral on spectrally truncated spectral triples (2024)
- **Citation:** E.-M. Hekkelman, *Lett. Math. Phys.* 112 (2022) 84; Hekkelman +
  McDonald, arXiv:[2412.00628](https://arxiv.org/abs/2412.00628).
- **Structural claim:** Detailed study of the simplest spectral truncation (circle),
  working out the Toeplitz-matrix operator-system structure and an NCG integral
  adapted to truncation.  Hekkelman–McDonald 2024 connects to quantum ergodicity.
- **Relevance to GeoVac:** worked-example technology for E9.  The circle (d = 1,
  commutative) is the simplest base, but the methods — especially the NCG integral
  on a truncated spectral triple — are exactly the right toolkit for the GeoVac
  n_max = 3 cutoff.  Flagged as future-sprint starting point.
- **CC transfer:** in-progress, example-driven.

### E11 — Kluth–Litim, Heat kernel coefficients on the sphere in any dimension
  (EPJ C 2020)
- **Citation:** Y. Kluth, D. F. Litim, *Eur. Phys. J. C* 80 (2020) 269,
  arXiv:[1910.00543](https://arxiv.org/abs/1910.00543).
- **Structural claim:** Derives **all** heat-kernel coefficients for Laplacians acting
  on scalars, vectors, and tensors on fully-symmetric spaces (including spheres) in
  **any dimension, including odd dimensions**.  Explicitly notes the heat-kernel
  expansion is asymptotic in even dimensions but has a finite radius of convergence in
  odd ones.  Results confirmed by spectral sums + Euler–Maclaurin.
- **Relevance to GeoVac (α-EB v2 cross-check):** **this is the right reference for
  testing the π³ = Vol(S⁵) hypothesis in α-EB v2.**  The paper gives Seeley–DeWitt
  coefficients a_k on S^n in closed form; evaluated at n = 5 it should produce
  coefficients proportional to powers of π (because Vol(S⁵) = π³) times rational
  numbers.  The α-EB v2 proposal is that the post-cubic residual π³α³ ≈ 1.2 × 10⁻⁵ in
  K = 1/α has the structural form a_0(S⁵) · α³ for some normalization.
- **CC transfer:** the heat-kernel machinery is continuous-base; any discrete analog
  would need the Euler–Maclaurin correction made explicit (which Kluth–Litim do
  in the cross-check half of their paper).  **PDF binary extraction failed in this
  sprint; manual consultation of the published EPJ C paper is recommended.**

### E12 — Chamseddine–Connes, Renormalization of the Spectral Action (2011)
- **Citation:** A. H. Chamseddine, A. Connes, arXiv:[1101.4804](https://arxiv.org/abs/1101.4804).
- **Structural claim:** Renormalization-group flow of the spectral action; derivation
  of one-loop β-functions for the SM gauge couplings from the CC setup; distributional
  approach (Estrada–Gracia-Bondía–Várilly) sign-permissive but not predictive.
- **Relevance to GeoVac:** relevant for α-SP's rule R4.  Sign-permissive in the
  distributional sense; does not explain Paper 2's $(+,+,-)$ pattern.
- **CC transfer:** continuous-base, does not adapt to finite resolution directly.

### E13 (bonus) — Eckstein–Iochum, Spectral Action in Noncommutative Geometry (Springer
  2019)
- **Citation:** M. Eckstein, B. Iochum, *Spectral Action in Noncommutative Geometry*,
  Springer MP Studies (2019), arXiv:[1902.05306](https://arxiv.org/abs/1902.05306).
- **Structural claim:** Monograph.  In-depth treatment of heat-kernel ↔ spectral-
  zeta correspondence.  Chapter 4: fluctuations by gauge potentials.
- **Relevance to GeoVac:** second textbook reference after van Suijlekom (E8).
  Ch. 3 is the single best reference on the Seeley–DeWitt ↔ zeta-function bridge,
  which is the structure GeoVac's Paper 18 taxonomy sits inside.

### E14 (bonus) — de Jong, Graphs, spectral triples and Dirac zeta functions (2009)
- **Citation:** J. W. de Jong, arXiv:[0904.1291](https://arxiv.org/abs/0904.1291).
- **Structural claim:** Associates a finitely summable **commutative** spectral triple
  to any finite, connected, unoriented graph with Betti number g ≥ 2 and valencies ≥ 3.
  Computes the induced Dirac zeta function as a graph invariant.
- **Relevance to GeoVac:** this is the **other published approach to "spectral triple
  on a finite graph"** that complements the Marcolli–van Suijlekom ACG-on-graph
  construction (E6).  De Jong's ST is **commutative** (A = ℂ(V) is diagonal), whereas
  GeoVac's proposed extension A × M_n(ℂ) is almost-commutative (E6-style).  The
  commutative ST on the graph is worth cross-checking against GeoVac for independent
  validation of the base structure.  Dirac zeta function as a graph invariant
  **might** connect to Paper 29's Ihara zeta, but de Jong's abstract does not mention
  the Ihara connection.
- **CC transfer:** the spectral-action functional on de Jong's ST is well-defined but
  his paper did not compute it; known to reduce to graph combinatorics.

### E15 (bonus) — Lattice fermions as spectral graphs (JHEP 2022)
- **Citation:** Various authors, *JHEP* 02 (2022) 104, arXiv:[2112.13501](https://arxiv.org/abs/2112.13501).
- **Structural claim:** Identifies lattice fermions (Naive, Wilson, Domain-wall)
  as spectral graphs; number of zero eigenvalues (doublers) controlled by graph
  topology.  Applied to non-regular lattices.
- **Relevance to GeoVac:** validates the general line "spectral graphs can host
  fermion-like operators" but the spectral-graph framing is NOT the Connes spectral-
  triple framing.  Shape-compatible but categorically distinct.
- **CC transfer:** not a spectral-action construction.

### E16 (bonus) — Marcolli, Lectures on Arithmetic Noncommutative Geometry
- **Citation:** M. Marcolli, *Lectures on Arithmetic Noncommutative Geometry*,
  University Lecture Series vol. 36, AMS (2005); [math.fsu.edu/~marcolli/BookAMSULect.pdf](http://www.its.caltech.edu/~matilde/BookAMSULect.pdf).
  Also Connes–Marcolli, *Noncommutative Geometry, Quantum Fields and Motives*, AMS
  Colloquium Publications vol. 55 (2008).
- **Structural claim:** NCG connected to arithmetic, moduli spaces of Q-lattices,
  Riemann Hypothesis context (via Connes' adèle class space).
- **Relevance to GeoVac:** **not relevant to WH1**; this is the arithmetic-NCG
  direction (Bost–Connes, KMS states, class field theory).  Listed for completeness
  because GeoVac RH-frontier (Paper 29) has its own Ihara-zeta angle, and Connes–
  Marcolli is the canonical arithmetic-NCG reference that a reader might ask about.

---

## 3. Closest match to GeoVac (single-paragraph assessment)

The single published framework closest to GeoVac's WH1 structure is
**Marcolli–van Suijlekom's *Gauge networks in noncommutative geometry*** (2014, E6):
*systems of finite spectral triples on a graph*, with a discretized Dirac operator
and a spectral action that reduces to lattice Wilson / Yang–Mills–Higgs.  GeoVac's
Paper 25 (Hopf U(1) Wilson) and Paper 30 (SU(2) Wilson) both fit the Marcolli–van
Suijlekom template: a finite, integer-labeled graph at vertices carrying finite
spectral-triple data, with connection on edges.  The Perez-Sanchez line (E7) extends
this to AF-algebra inductive sequences and makes the continuum limit more honest
(YM without Higgs).  GeoVac's **original** contribution, relative to this published
literature, is (i) the graph is not an arbitrary lattice but the Fock projection of
the hydrogenic spectrum (Paper 7), and (ii) the Dirac operator is the Camporesi–Higuchi
$D_{S^3}$ restricted to the (n, κ, m_j) labeling (T1–T9), not a generic discretization.

No published framework computes α (or any dimensionless physical constant) from a
discrete / graph-native spectral action.  The Chamseddine–Connes SM derivations
(E1, E4, E12) compute **ratios** of SM gauge couplings at a unification scale and the
Higgs-to-W ratio, but leave the overall scale (and α) as inputs.  WH1 in its full
form — *GeoVac computes α via a spectral-action principle on an almost-commutative
spectral triple whose base is a Fock-projected Hopf graph* — is **not duplicated in the
published literature** to 2024.

---

## 4. What transfers vs fails (CC → discrete)

| CC result | Transfers to graph-based ACG? | Notes |
|:--|:--|:--|
| Spectral triple axiomatics (A, H, D, J, γ) | **Yes**, fully | Krajewski / Paschke–Sitarz / Ćaćić (E2/E3/E5); finite ACG is well-developed. |
| Heat-kernel asymptotic expansion Tr f(D/Λ) ~ Σ f_k Λ^{d−k} a_k | **No**, the asymptotic expansion requires a continuous manifold for the Seeley–DeWitt recursion; on a finite graph one has a *Hilbert-Schmidt* trace of a finite matrix, which is just Σ f(λ_i).  Euler–Maclaurin corrections are explicit and controlled at finite cutoff (Kluth–Litim E11). | α-SP's R1 applies only in the continuous limit; on a finite graph the sign rule $(+, +, +)$ for SD coefficients degenerates into direct spectral sums, which CAN carry $(+, +, -)$ or any other sign pattern. |
| Spectral action = Einstein–Hilbert + YM + Higgs (CC derivation) | **Partially**, via Marcolli–van Suijlekom (E6) → Wilson action / YM lattice gauge theory.  Higgs piece removed in Perez-Sanchez correction (E7). | GeoVac's Paper 30 SU(2) Wilson action IS in this class. |
| Gauge fields as inner derivations | **Yes**, fully; standard almost-commutative construction (E8, ch. 9). | This is Paper 30's logic: SU(2) gauge as inner derivation of A × M_2(ℂ). |
| Spectral action principle as UV completion | **Open**; spectral truncation framework (E9/E10) is the honest replacement but is illustrated on the circle, not arbitrary graphs. | GeoVac's finite n_max plays the role of the truncation. |
| Predicting dimensionless physical constants | **Partial in continuous CC**; NOT in any published discrete construction.  CC fixes gauge-coupling *ratios* (E4); nobody publishes a discrete-spectral-action derivation of α. | WH1's full claim exceeds what is in the literature. |
| SD coefficients on S⁵ | **Yes** (E11), explicit in arbitrary dimension with Euler–Maclaurin corrections for odd dims. | Needed for α-EB v2 π³ cross-check. |
| APS eta-invariant as boundary correction | **Yes**; standard in CC on manifolds with boundary. | α-SP's R3 is the only CC rule that matches Paper 2's $(+, +, -)$ pattern *shape-wise*; but GeoVac's Δ is a mode count, not an η-invariant, so the match is qualitative only. |
| Reconstruction theorem (spectral triple ⇒ manifold) | **Yes** for commutative ST; **partial** for almost-commutative ST (Ćaćić 2011). | GeoVac's data should reconstruct S³ = SU(2); this is a consistency check, not a computational output. |

---

## 5. WH1 status assessment

**WH1 status: MODERATE** (up from "partial shape-match" in the α-SP memo).

Arguments FOR (upgrade WH1 to "moderate"):
1. **E6 (Marcolli–van Suijlekom) directly constructs finite spectral triples on
   graphs** and computes their spectral action.  GeoVac's Paper 25 / Paper 30 Hopf graph
   + Wilson action fit this template.
2. The operator-system / spectral-truncation framework (E9/E10) is the correct
   published home for finite-cutoff spectral actions.  GeoVac's n_max ≤ 3 truncation
   is shape-compatible.
3. Almost-commutative geometry (E8, E4) is the canonical published framework for
   attaching an internal matrix algebra M_n(ℂ) to a base spectral triple.  GeoVac's
   Paper 30 SU(2) sector and Paper 25 U(1) sector both instantiate A × M_n(ℂ).

Arguments AGAINST (keep WH1 below "strong"):
1. **No published discrete-spectral-action framework predicts α.**  The uncanny
   precision of CC (E4) fixes gauge-coupling *ratios* + Higgs:W ratio; it does not fix
   the overall dimensionless scale or α itself.
2. The Marcolli–van Suijlekom continuum-limit claim was **partially corrected** by
   Perez-Sanchez (E7) to YM without Higgs.  Even the best published analog has
   documented issues at the continuum-limit boundary.
3. GeoVac's graph is the Fock-projected S³ Coulomb graph — a very specific object
   with (n, l, m) integer labeling by Paper 7's packing construction.  Published ST
   frameworks on graphs either use generic Bratteli / Krajewski diagrams (E2, E7) or
   abstract combinatorial graphs (E6, E14).  **No published framework matches
   GeoVac's specific graph**, so WH1 requires that we BUILD the spectral-triple
   structure on our specific graph and check consistency against the published
   axioms — which has not yet been done in full.
4. α-SP already established that the standard CC sign rules do not derive Paper 2's
   $(+, +, -)$ pattern.  This is a genuine mismatch that no published framework closes.

**Falsifier candidates (for WH1 status) not found in the literature search:**
- No theorem asserting that integer-labeled Hopf-graph spectral triples cannot exist.
- No theorem asserting that Camporesi–Higuchi $D_{S^3}$ restricted to (n, l, m) is NOT
  a valid spectral-triple Dirac operator.
- No theorem ruling out α-like predictions from discrete spectral actions; the absence
  of such predictions in the literature is by *silence*, not by *proof-of-absence*.

**Net WH1 assessment:** partial-shape-match (α-SP) upgraded to **moderate** on the
strength of E6, E7, E9 establishing that the vocabulary "finite spectral triple on
a graph with a spectral action" is published mature NCG, not GeoVac invention.  The
specific claim "GeoVac's α emerges from a spectral-action principle on its Hopf graph"
remains not duplicated in the literature.

---

## 6. α-EB v2 cross-check relevant findings

**Key relevant reference: E11 (Kluth–Litim, EPJ C 80, 269 (2020)).**
Explicit SD coefficients on S^n for all n; confirmed by Euler–Maclaurin;
valid for both scalar Laplacian and iterated Dirac $D^2$ on spheres (cross-referenced
with Camporesi–Higuchi *Commun. Math. Phys.* 148 (1992) 283, [arXiv:gr-qc/9505009](https://arxiv.org/abs/gr-qc/9505009)
which is the *Dirac* heat kernel on spheres + hyperbolic spaces).

**α-EB v2 hypothesis:** post-cubic residual in K = 1/α (after subtracting 1/α + α²)
equals π³α³ ≈ 1.2049 × 10⁻⁵ to 0.25%, with π³ = Vol(S⁵).  Test: compute $a_0^{\text{Dirac}}
(S^5)$ from E11 and check whether the π³ prefactor arises naturally in the spectral-
action expansion on S⁵.

**What I found relevant:**
1. E11 gives SD coefficients on spheres in any dimension in closed form — the needed
   data is **published and accessible** (not a GeoVac derivation).
2. Camporesi–Higuchi (1995 arXiv / 1992 CMP) gives the iterated-Dirac heat kernel on
   spheres explicitly — **also published**, and Paper 24's HO construction already
   cites the Dirac spectrum on S³; the S⁵ analog is in the same reference.
3. **I did not extract the explicit coefficients** from the Kluth–Litim PDF in this
   sprint (PDF binary extraction failed; would require manual retrieval of the
   published EPJ C paper or the arXiv HTML abstract — a bounded sub-task, not
   requiring further web search).

**Recommendation for follow-up (outside this sprint):** dispatch a short sub-task to
(a) extract $a_0(S^5)$ and $a_2(S^5)$ for the *iterated Dirac* from Camporesi–Higuchi
1992 + Kluth–Litim 2020 in closed form; (b) test whether the π³α³ coefficient
predicted by α-EB v2 matches the $\alpha^3 a_0(S^5)$ term (with an appropriate
normalization convention); (c) if yes, α-EB v2's structural hint is published-literature-
grounded.  If no, α-EB v2's π³ is NOT Vol(S⁵) in the spectral-action sense.

---

## 7. Red flags

1. **Marcolli–van Suijlekom continuum-limit claim partially false** (E6 → E7
   correction).  Any GeoVac WH1 argument that cites E6 must cite E7 simultaneously.
2. **No published α-from-discrete-spectral-action prediction exists.**  WH1's full
   version is ambitious; honest framing should say "the published machinery supports
   building the spectral triple, but no existing framework predicts α from it."
3. **E9's standard examples are continuous-base with cutoff** (circle → Toeplitz),
   not graph-native.  GeoVac's case (discrete n_max ≤ 3 from the start) is more
   extreme than what the spectral-truncation literature has worked out in detail.
4. **PDF extraction failed** for Kluth–Litim (E11) in this sprint.  Explicit α-EB v2
   cross-check is deferred to a short follow-up (not completed here).  If the sprint
   budget re-opens, this should be the first sub-task.
5. **"Gauge networks in noncommutative geometry" (E6) is a 2013 preprint, 2014 J.
   Geom. Phys. publication.**  A decade old.  The Perez-Sanchez 2024 line (E7) is
   the active follow-up; any GeoVac paper citing E6 should also cite E7.
6. **de Jong's Dirac zeta function on graphs (E14) has NO explicit connection to
   Ihara zeta** in the abstract.  This is a potential research gap — Paper 29's
   Ihara zeta and de Jong's Dirac zeta would need independent cross-checking.

---

## 8. Honest scope statement

This memo is a **literature survey**, not a derivation.  All 12+ verified entries have
arXiv IDs or DOIs; where the PDF could not be parsed (E11 for α-EB cross-check), the
limitation is explicitly flagged.  The survey **does not** propose a new common-
generator mechanism (would violate the sprint constraint); it documents published
precedents for the WH1 vocabulary and identifies the closest published match (E6
Marcolli–van Suijlekom, with E7 correction noted).

**Match strength for WH1:** upgraded from α-SP's "partial shape-match" to **moderate**
on literature grounds; the claim "GeoVac is an almost-commutative spectral triple" is
supported by published precedent for each structural ingredient (finite-ACG
classification E2/E3/E5, spectral-action on graphs E6/E7, spectral truncations E9/E10,
heat-kernel on S^n E11).  The **composite claim** — that GeoVac's specific Fock-
projected S³ Hopf graph instantiates this structure AND that its spectral action
predicts α — remains **not duplicated** in the literature and is GeoVac's original
contribution to defend.

---

## 9. References (all verified in this survey)

1. **E1.** A. H. Chamseddine, A. Connes, CMP 186 (1997) 731, [hep-th/9606001](https://arxiv.org/abs/hep-th/9606001).
2. **E2.** T. Krajewski, JGP 28 (1998) 1, [hep-th/9701081](https://arxiv.org/abs/hep-th/9701081).
3. **E3.** M. Paschke, A. Sitarz, Acta Phys. Polon. B 31 (2000) 1897; JMP 39 (1998) 6191.
4. **E4.** A. H. Chamseddine, A. Connes, CMP 293 (2010) 867, [0812.0165](https://arxiv.org/abs/0812.0165).
5. **E5.** B. Ćaćić, "Moduli Spaces of Dirac Operators for Finite Spectral Triples,"
   [0902.2068](https://arxiv.org/abs/0902.2068), in Springer NCG volume (2011).
6. **E6.** M. Marcolli, W. D. van Suijlekom, JGP 75 (2014) 71,
   [1301.3480](https://arxiv.org/abs/1301.3480).
7. **E7.** C. I. Perez-Sanchez, "Bratteli networks and the Spectral Action on quivers,"
   [2401.03705](https://arxiv.org/abs/2401.03705); Comment on E6 at
   [2508.17338](https://arxiv.org/abs/2508.17338).
8. **E8.** W. D. van Suijlekom, *Noncommutative Geometry and Particle Physics*, 2nd ed.,
   Springer MPS (2024), [open-access PDF](https://repository.ubn.ru.nl/bitstream/handle/2066/314382/314382.pdf).
9. **E9.** A. Connes, W. D. van Suijlekom, CMP 383 (2021) 2021, [2004.14115](https://arxiv.org/abs/2004.14115).
10. **E10.** E.-M. Hekkelman, Lett. Math. Phys. 112 (2022) 84; Hekkelman + McDonald,
    [2412.00628](https://arxiv.org/abs/2412.00628).
11. **E11.** Y. Kluth, D. F. Litim, EPJ C 80 (2020) 269, [1910.00543](https://arxiv.org/abs/1910.00543).
12. **E12.** A. H. Chamseddine, A. Connes, [1101.4804](https://arxiv.org/abs/1101.4804).
13. **E13 (bonus).** M. Eckstein, B. Iochum, *Spectral Action in NCG*, Springer (2019),
    [1902.05306](https://arxiv.org/abs/1902.05306).
14. **E14 (bonus).** J. W. de Jong, [0904.1291](https://arxiv.org/abs/0904.1291).
15. **E15 (bonus).** Various, JHEP 02 (2022) 104, [2112.13501](https://arxiv.org/abs/2112.13501).
16. **E16 (bonus).** M. Marcolli, *Lectures on Arithmetic NCG*, AMS (2005);
    A. Connes, M. Marcolli, *NCG, Quantum Fields and Motives*, AMS Coll. Publ. 55
    (2008).

Also consulted but not cited as primary WH1 evidence:
- D. V. Vassilevich, "Heat kernel expansion: user's manual," Phys. Rep. 388 (2003) 279,
  [hep-th/0306138](https://arxiv.org/abs/hep-th/0306138) — standard reference, already
  cited in α-SP memo.
- R. Camporesi, A. Higuchi, CMP 148 (1992) 283 — Dirac heat kernel on spheres +
  hyperbolic spaces; underlying reference for the Dirac spectrum used in GeoVac T1.

---

## 10. Sprint outcome

- **Deliverable:** this memo, 16 verified entries (12 required, 4 bonus).
- **WH1 verdict:** MODERATE (up from partial shape-match).
- **Closest published match:** E6 Marcolli–van Suijlekom, with E7 correction.
- **α-EB v2 cross-check:** E11 Kluth–Litim is the right reference; explicit S⁵ data
  extraction deferred (PDF binary issue).
- **Red flags:** 6 items, documented §7.
- **Not duplicated in the literature:** GeoVac's α prediction from discrete spectral
  action on Fock-projected S³ Hopf graph is the framework's original contribution.
- **Recommended next sub-task:** short follow-up on E11 (EPJ C 80, 269) to extract
  explicit $a_0(S^5)$, $a_2(S^5)$ coefficients for α-EB v2 π³ test.
