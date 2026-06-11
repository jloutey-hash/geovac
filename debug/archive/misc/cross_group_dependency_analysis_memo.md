# Cross-Group Dependency Analysis — GeoVac Paper Series

**Date:** 2026-05-22
**Scope:** Internal citation dependencies across the six audience-targeted paper groups under `papers/` (synthesis folder empty; archive treated separately). Read-only analysis; no .tex files modified.
**Coverage:** ~85% of internal GeoVac series citations. Methodology combines bibitem citations (`\cite{...}` with the conventions `paper_N`, `Paper_N`, `loutey_paper_N`, `GeoVac_Paper_N`) with prose references (`Paper~N` / `Paper N`). Narrative dependencies were merged into the per-paper citation lists because Paper 34 in particular relies heavily on prose references without expanding the bibitem table to match. Counts are structural indicators, not precise bibliometrics.

---

## Section 1: Executive Summary

The most central paper of each group, ranked by inbound GeoVac-internal citations:

- **Group 1 (math.OA):** Paper 32 (spectral-triple synthesis, 14 inbound) — the operator-algebras program's hub.
- **Group 2 (quantum chemistry):** Paper 11 (prolate spheroidal H₂⁺, 13 inbound) — the foundational diatomic result.
- **Group 3 (foundations):** Paper 7 (S³ dimensionless vacuum, 31 inbound) — by far the most-cited paper in the entire series.
- **Group 4 (quantum computing):** Paper 14 (qubit encoding, 22 inbound) — the headline computational result.
- **Group 5 (QED / gauge):** Paper 2 (α observation, 15 inbound) — the conjectural anchor for QED on S³.
- **Group 6 (precision / observations):** Paper 34 (projection taxonomy, 8 inbound) — the living catalogue.

The highest-weight cross-group edges are (i) **Group 1 → Group 3 (17 edges)**: math.OA papers depend on Papers 7, 24, 18 for the underlying spectral-triple substrate; (ii) **Group 3 → Group 2 (22) and Group 2 → Group 3 (27)**: the foundations–chemistry coupling is mutual and dense, anchored by Papers 7, 0, 1 cited from chemistry, and Papers 11, 13, 17 cited from Paper 18's exchange-constant taxonomy; (iii) **Group 5 → Group 3 (27 edges)**: QED on S³ rests on Papers 7, 18, 22, 24; (iv) **Group 6 → Group 3 (17 edges)**: precision catalogues (Paper 34, 35) cite Papers 18, 22, 24, 7, 31 as the structural backbone; (v) **Group 1 → Group 1 (34 edges, intra-group)**: the math.OA papers form a tight forward chain 38→39→40→42→43→44→45 with Paper 32 as connector. Paper 7 is the single most universal dependency in the framework — cited by 31 of 39 papers.

---

## Section 2: 6×6 Group-to-Group Citation Matrix

Rows = source group, columns = target group. Cell counts include both bibitem and narrative dependencies; intra-group counts on the diagonal.

|        | G1 (math.OA) | G2 (chem) | G3 (found) | G4 (QC) | G5 (QED) | G6 (obs) |
|--------|:---:|:---:|:---:|:---:|:---:|:---:|
| **G1 (math.OA)** | 34 | 0 | 17 | 3 | 11 | 3 |
| **G2 (chem)** | 0 | 23 | 27 | 8 | 0 | 0 |
| **G3 (found)** | 7 | 22 | 26 | 11 | 13 | 6 |
| **G4 (QC)** | 3 | 8 | 10 | 3 | 1 | 2 |
| **G5 (QED)** | 7 | 1 | 27 | 10 | 17 | 3 |
| **G6 (obs)** | 7 | 5 | 17 | 7 | 10 | 6 |

**Row totals (citations made by each group):** G1=68, G2=58, G3=85, G4=27, G5=65, G6=52.
**Column totals (citations received by each group):** G1=58, G2=59, G3=124, G4=42, G5=52, G6=20.

Key observations:
- **Group 3 (foundations) is the gravitational center of the series**: 124 inbound citations across 9 papers (≈14 per paper average), far above any other group.
- **Group 2 (chemistry) is structurally isolated from G1, G5, G6**: zero edges out of G2 into G1/G5/G6. The chemistry papers were written first and don't refer to the math.OA / QED / observations program; the dependency is one-way (later papers reach back).
- **Group 6 (observations) makes citations everywhere but receives only 20 inbound** — it is a sink for catalogue-style writeups that synthesize earlier work; its main upstream readers are the math.OA papers (G1 → G6 = 3) and the QED papers (G5 → G6 = 3, both citing Paper 34/35 as taxonomy).
- **G4 (quantum computing) is the smallest source and target group**, but Paper 14 alone is one of the most-cited papers (22 inbound) — Group 4's role in the series is concentrated in a single hub paper.

---

## Section 3: Per-Paper Centrality Within Each Group

Inbound counts from all other GeoVac papers in the series. Each entry lists the paper, the inbound count, and a one-line reason it gets cited.

### Group 1 — Operator algebras

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **32** | 14 | Spectral-triple synthesis; canonical reference for (A_GV, H_GV, D_GV) construction, Connes axiom audit, case-exhaustion theorem §VIII |
| **38** | 11 | WH1 PROVEN; SU(2) propinquity convergence theorem with five-lemma proof, 4/π universal rate |
| **40** | 8 | Unified GH-convergence on compact Lie groups; structural generalization of Paper 38 to all rank-r |
| **39** | 7 | Tensor-product propinquity convergence; cross-focal-length theorem |
| **42** | 6 | Tomita-Takesaki modular Hamiltonian, four-witness Wick-rotation literal identification at finite cutoff |
| **43** | 6 | Lorentzian extension to Krein space at signature (3,1); H_local ≠ D_W structural finding |
| **29** | 5 | Ramanujan Hopf graphs, Ihara zeta; graph-RH program |
| **44** | 1 | Operator-system substrate at signature (3,1); cited only by Paper 45 |
| **45** | 0 | Most recent (K⁺-weak-form Lorentzian propinquity); no inbound yet |

### Group 2 — Quantum chemistry

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **11** | 13 | Prolate spheroidal H₂⁺; foundational diatomic with spectral Laguerre; cited by every chemistry paper and several others |
| **13** | 11 | Hyperspherical He; level 3 architecture; cited by Papers 7, 14, 17, 31 |
| **17** | 10 | Composed geometries (LiH/BeH₂/H₂O); cited by Papers 14, 19, 20, 26, 27, 31, 34 |
| **15** | 8 | Level 4 molecular geometry (H₂); cited by Papers 0, 17, 31, FCI-M |
| **19** | 7 | Coupled composition; cited by Papers 14, 17, 20, 23, 27 |
| **12** | 6 | Neumann V_ee expansion; cited by chemistry papers 13, 15, 17 and Paper 18 |
| **8** | 4 | Bond sphere Sturmian (GUARDRAIL); cited by Papers 11, FCI-A, FCI-M, 36 |
| **FCI-A** | 0 | Atomic FCI benchmark suite; standalone, not cited by name |
| **FCI-M** | 0 | LCAO molecular FCI dead-end (GUARDRAIL); standalone reference |

### Group 3 — Foundations / spectral graph

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **7** | 31 | S³ dimensionless vacuum proof; cited by 31 of 39 papers in the series — most universal dependency |
| **18** | 22 | Exchange constants taxonomy; cited from every group, anchors transcendental classification |
| **0** | 20 | Geometric packing axioms, K = −1/16; cited from every group |
| **24** | 18 | Bargmann-Segal lattice, HO rigidity, Coulomb/HO asymmetry; cited heavily by G1 and G5 |
| **22** | 13 | Angular sparsity theorem; cited by G4 (Paper 14), G5, G6 |
| **1** | 12 | Spectral graph theory; cited by chemistry and Paper 31 |
| **31** | 6 | Universal/Coulomb partition; cited by Papers 32, 33, 34, 35 |
| **21** | 1 | Synthesis paper; cited by Paper 25 only |
| **6** | 1 | Quantum dynamics; cited by Paper 0 only |

### Group 4 — Quantum computing

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **14** | 22 | Qubit encoding, O(Q²·⁵) Pauli scaling, composed sparsity; flagship computational result |
| **23** | 11 | Nuclear shell model qubit Hamiltonians; cited by Papers 18, 22, 24, 27, 30, 31, 32, 34, 35, 36 |
| **20** | 5 | Resource benchmarks; cited by Papers 14, 17, 19, 21, 23 |
| **16** | 4 | Chemical periodicity / S_N rep theory; cited by Papers 17, 19, 21, 25 |

### Group 5 — QED / gauge

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **2** | 15 | α conjecture; cited by every group except quantum computing, anchors transcendental program |
| **28** | 12 | QED on S³; cited by Papers 18, 25, 29, 30, 32, 33, 34, 35, 36, 41, 43 |
| **25** | 11 | Hopf gauge structure synthesis; cited by Papers 18, 24, 28, 29, 30, 33, 34, 41 |
| **30** | 9 | SU(2) Wilson lattice gauge; cited by Papers 18, 24, 25, 28, 31, 33, 34, 41 |
| **36** | 3 | Bound-state QED Lamb shift one-loop closure; cited by Papers 18, 34, 35 |
| **33** | 2 | 1+6+1 QED selection rule partition; cited by Papers 34, 35 |
| **41** | 0 | Rule B Wilson U(1) — most recent QED paper, no inbound yet |

### Group 6 — Precision / observations

| Paper | Inbound | Reason cited |
|---|:---:|---|
| **34** | 8 | Projection taxonomy (28 projections, three-axis tagging); cited by Papers 18, 23, 27, 31, 32, 35, 36, 42 |
| **35** | 6 | Klein-Gordon spectrum, π-free graph, observation-vs-rest-mass split; cited by Papers 18, 31, 32, 33, 34, 36, 41 |
| **27** | 4 | Entropy as projection artifact; cited by Papers 18, 24, 26, 34 |
| **26** | 2 | Entanglement structure; cited by Papers 27, 34 |

---

## Section 4: Cross-Group Edge Catalogue

Enumeration of non-diagonal cells with ≥3 edges; each line is "Source paper (Group) → Target paper (Group): brief reason / load-bearing class". Tags: **S** = structural (theorem reused), **M** = methodological (technique borrowed), **C** = contextual (acknowledgement only).

### G1 → G3 (17 edges)

Math.OA papers consistently cite the foundational substrate.

- Paper 29 → Paper 7 (G3): S³ Fock projection underlies the Hopf-graph construction (**S**).
- Paper 29 → Paper 18 (G3): exchange-constant taxonomy, π-free certificate (**S**).
- Paper 29 → Paper 24 (G3): Bargmann-Segal S⁵ comparison (**S**).
- Paper 29 → Paper 1 (G3): spectral graph methodology (**M**).
- Paper 32 → Paper 7 (G3): D_GV constructed on the S³ Fock graph (**S**).
- Paper 32 → Paper 22 (G3): angular sparsity theorem feeds into operator-system propagation (**S**).
- Paper 32 → Paper 24 (G3): Coulomb/HO asymmetry, four-layer §VIII.B/§VIII.C (**S**).
- Paper 32 → Paper 31 (G3): universal/Coulomb partition framing (**S**).
- Paper 38 → Paper 7 (G3): Camporesi-Higuchi Dirac on S³ comes from Paper 7's S³ proof (**S**).
- Paper 38 → Paper 24 (G3): Coulomb/HO asymmetry sets scope for SU(2) GH convergence (**M**).
- Paper 39 → Paper 24 (G3): cross-manifold W2b blocker explicitly out of scope (**C**).
- Paper 42 → Paper 24 (G3): scope cross-reference for Lorentzian extension (**C**).
- Paper 43 → Paper 18 (G3): exchange-constant taxonomy hosts the Pythagorean M1 prefactor (**S**).
- Paper 43 → Paper 24 (G3): Pythagorean orthogonality four-layer Coulomb/HO connection (**S**).
- Paper 44 → Paper 24 (G3): same context (**C**).
- Paper 45 → Paper 24 (G3): same context (**C**).
- Paper 42 → Paper 31 (G3): KO-dimension argument cites universal-Coulomb partition (**M**).

### G1 → G5 (11 edges)

- Paper 29 → Paper 25 (G5): Hopf gauge structure shares the Ramanujan-graph Hopf bundle (**S**).
- Paper 29 → Paper 28 (G5): vertex parity / QED on S³ (**S**).
- Paper 32 → Paper 2 (G5): K = π(B+F−Δ) three-sector reading (**S**).
- Paper 32 → Paper 25 (G5): U(1) Wilson construction as Cartan-subalgebra projection (**S**).
- Paper 32 → Paper 28 (G5): QED on S³ as projection of the GeoVac triple (**S**).
- Paper 32 → Paper 30 (G5): SU(2) Wilson as non-abelian sibling on S³=SU(2) (**S**).
- Paper 32 → Paper 36 (G5): bound-state QED closure target (**S**).
- Paper 38 → Paper 2 (G5): the α conjecture motivates the propinquity proof (**C**).
- Paper 42 → Paper 2 (G5): same (**C**).
- Paper 42 → Paper 25 (G5): cross-reference (**C**).
- Paper 43 → Paper 28 (G5): QED on S³ as Lorentzian-extension target (**M**).

### G3 → G2 (22 edges)

Foundations cite chemistry results extensively because the natural geometry hierarchy is anchored in concrete molecular benchmarks.

- Paper 0 → Papers 11, 13, 15, 17 (G2): packing axioms grounded in concrete level results (**S** all).
- Paper 7 → Paper 11 (G2): H₂⁺ benchmark (**S**).
- Paper 18 → Papers 11, 12, 13, 15, 17, 19 (G2): every chemistry paper appears in the exchange-constant inventory (**S** all).
- Paper 21 → Papers 11, 13, 15, 17, 19 (G2): synthesis paper cites the full natural geometry hierarchy (**S** all).
- Paper 22 → none (no G2 cites).
- Paper 24 → none (G2 cites only via Paper 27).
- Paper 31 → Papers 11, 13, 15, 17 (G2): universal-Coulomb partition uses Level 2-4 benchmarks (**S** all).

### G3 → G3 (26 edges, intra-group)

Foundations papers form a tight cluster. The largest sub-hub is Paper 18's pull from all other foundations papers.

### G3 → G4 (11 edges)

- Paper 18 → Paper 14 (G4): qubit encoding hosts the spinor-block exchange constants (**S**).
- Paper 18 → Paper 23 (G4): nuclear shell qubit Hamiltonian as universal-Coulomb partition case study (**S**).
- Paper 22 → Paper 14 (G4): O(Q²·⁵) Pauli scaling is the angular-sparsity theorem applied (**S**).
- Paper 22 → Paper 23 (G4): nuclear validation of angular sparsity (**S**).
- Paper 24 → Paper 14 (G4): Bargmann-Segal qubit applications (**M**).
- Paper 24 → Paper 23 (G4): nuclear shell extension (**S**).
- Paper 31 → Paper 14 (G4): qubit work as universal-sector application (**S**).
- Paper 31 → Paper 23 (G4): nuclear as second universal-sector test (**S**).
- Paper 7 → Paper 14 (G4): qubit encoding derived from S³ proof (**S**).
- Paper 0 → Paper 14, 16 (G4): packing/periodicity (**S** both).

### G3 → G5 (13 edges)

- Paper 7 → Paper 2 (G5): conformal equivalence underpins α conjecture (**S**).
- Paper 18 → Papers 25, 28, 29, 30 (G5): exchange-constant taxonomy hosts QED on S³ (**S** all).
- Paper 21 → Paper 2 (G5): synthesis (**C**).
- Paper 22 → none (G5 cites G3 heavily, not the reverse here).
- Paper 24 → Papers 2, 25, 30 (G5): Coulomb/HO asymmetry constrains gauge content (**S** all).
- Paper 31 → Papers 2, 25, 28, 30 (G5): universal-Coulomb partition organizes Wilson constructions (**S** all).

### G3 → G6 (6 edges)

- Paper 18 → Papers 27, 34, 35 (G6): exchange-constant taxonomy is the substrate of the projection catalogue (**S** all).
- Paper 24 → Paper 27 (G6): HO entanglement rigidity (**S**).
- Paper 31 → Papers 34, 35 (G6): universal/Coulomb partition feeds projection taxonomy (**S** both).

### G5 → G3 (27 edges)

QED papers consistently cite the foundational substrate. The pattern is uniform: every QED paper cites Papers 7, 18, plus a subset of {22, 24}.

- Paper 2 → Papers 7, 18, 22, 23, 24 (5 edges total).
- Paper 25 → Papers 0, 7, 18, 21, 22, 24 (6 edges).
- Paper 28 → Papers 7, 18 (2 edges).
- Paper 30 → Papers 0, 7, 18, 22, 24 (5 edges).
- Paper 33 → Papers 7, 18, 22, 24, 31 (5 edges).
- Paper 36 → Papers 7, 18 (2 edges).
- Paper 41 → Papers 7, 18, 24 (3 edges).

All structural/methodological; foundations provide the operator-algebraic substrate.

### G5 → G4 (10 edges)

- Paper 2 → Paper 14 (G4): qubit encoding hosts the α-conjecture spectral data (**S**).
- Paper 25 → Paper 14 (G4): Hopf gauge in qubit encoding (**M**).
- Paper 28 → Paper 14 (G4): QED scaling validated against Pauli scaling (**S**).
- Paper 30 → Paper 14 (G4): same (**M**).
- Paper 33 → Paper 14 (G4): selection rules in qubit space (**S**).
- Paper 36 → Paper 14 (G4): bound-state QED with composed encoding (**S**).
- Paper 41 → Paper 14 (G4): Wilson U(1) qubit application (**M**).
- Paper 28 → Paper 23 (G4): nuclear cross-reference (**C**).
- Paper 36 → Paper 23 (G4): nuclear context (**C**).
- Paper 33 → Paper 14 / Paper 23 (additional cross-reference, **C**).

### G5 → G5 (17 edges, intra-group)

Tight QED cluster: every QED paper cites Paper 2 and Paper 28; Papers 25, 28, 30, 29 (G1) form the gauge sub-cluster.

### G6 → G3 (17 edges)

Observations consistently anchor on foundations.

- Paper 26 → Paper 7, 22 (G3) (2 edges).
- Paper 27 → Papers 0, 18, 24 (G3) (3 edges; also G6's Papers 5(archived), 23, 26).
- Paper 34 → Papers 0, 7, 18, 22, 24, 31 (G3) (6 edges).
- Paper 35 → Papers 0, 7, 18, 22, 24, 31 (G3) (6 edges).

All structural; the foundations provide the projection-taxonomy substrate.

### G6 → G5 (10 edges)

- Paper 27 → none G5.
- Paper 34 → Papers 2, 25, 28, 30, 33, 36 (G5) (6 edges).
- Paper 35 → Papers 2, 28, 33, 36 (G5) (4 edges).

The catalogue papers (34, 35) extensively reference QED/gauge results as named projection mechanisms.

### G6 → G4 (7 edges)

- Paper 26 → Paper 14, 17 (G4: 14 only? note: 17 is G2). Paper 26 cites Paper 14, 22, 27 internally; G4 hit is Paper 14.
- Paper 27 → Paper 23 (G4) (1).
- Paper 34 → Papers 14, 16, 20, 23 (G4) (4 edges).
- Paper 35 → Paper 14, 23 (G4) (2 edges).

### G6 → G6 (6 edges, intra-group)

Paper 26 → 27; Paper 27 → 26; Paper 34 → 26, 27, 35; Paper 35 → 34.

### G2 → G3 (27 edges)

Chemistry papers consistently cite the foundational triad (Papers 0, 1, 7) plus relevant substrate.

- All 9 chemistry papers cite Paper 7 (**S** all).
- 7 of 9 cite Paper 0 (**S** all).
- 6 of 9 cite Paper 1 (**S/M** all).
- Paper 13 cites Paper 18 (one of the few G2 → G3-18 edges) (**M**).
- Paper 11 → Paper 18 not present in this dataset; Paper 13 → Paper 18 is the main connection.

### G2 → G4 (8 edges)

- Paper 8 → Paper 14 (G4) (**M**).
- Paper 15 → Paper 14 (G4) (**M**).
- Paper 17 → Papers 14, 16, 20 (G4) (**M** all).
- Paper 19 → Papers 14, 16, 20 (G4) (**M** all).

### G4 → G3 (10 edges) and G4 → G2 (8 edges)

- Paper 14 → Papers 7, 13, 17, 18, 22 — anchors qubit encoding in foundations + chemistry.
- Paper 16 → Papers 0, 1, 7, 13 (mostly G3).
- Paper 20 → Papers 14, 17, 19 — applications paper depends on QC + chemistry.
- Paper 23 → Papers 0, 2, 7, 14, 17, 18, 22, 32, 34, 38, 39 — Paper 23 is the broadest G4 paper, drawing from every group.

---

## Section 5: Reading Prerequisites Per Group

Recommended 3–5 paper prerequisite set for someone diving into each group cold. Order = suggested reading order.

### Prereq for Group 1 (math.OA)

1. **Paper 7** (G3): the S³ Fock projection establishes the underlying manifold and Camporesi-Higuchi Dirac spectrum.
2. **Paper 18** (G3): exchange-constant taxonomy, master Mellin engine framing, transcendental classification — without this, propinquity-rate constants like 4/π read as mysterious.
3. **Paper 24** (G3): Bargmann-Segal lattice and the four-layer Coulomb/HO asymmetry — the scope boundary for all math.OA-side cross-manifold work.
4. **Paper 32** (G1): GeoVac spectral triple construction and Connes axiom audit; the synthesis paper for the math.OA program (intra-group prereq).
5. **Paper 2** (G5): the α conjecture — motivates K = π(B+F−Δ) and explains why the math.OA papers exist.

### Prereq for Group 2 (chemistry)

1. **Paper 0** (G3): packing axioms, K = −1/16.
2. **Paper 1** (G3): spectral graph theory and O(N) eigenvalue method.
3. **Paper 7** (G3): S³ proof; SO(4) symmetry.
4. **Paper 11** (G2): foundational diatomic (prolate spheroidal H₂⁺) — first natural geometry beyond S³.
5. **Paper 17** (G2): composed geometries — the production architecture (LiH/BeH₂/H₂O).

### Prereq for Group 3 (foundations)

The foundations are largely self-contained. For a cold reader:

1. **Paper 0** (G3): packing axioms.
2. **Paper 7** (G3): S³ dimensionless vacuum proof.
3. **Paper 11** (G2): H₂⁺ benchmark (referenced by every G3 paper).
4. **Paper 13** (G2): hyperspherical He (referenced by Papers 1, 7, 18, 21, 22, 31).
5. **Paper 18** (G3): exchange-constant taxonomy, which is the conceptual hub of the foundations group.

### Prereq for Group 4 (quantum computing)

1. **Paper 7** (G3): S³ Fock projection.
2. **Paper 17** (G2): composed geometry that all qubit encodings sit on.
3. **Paper 22** (G3): angular sparsity theorem — explains O(Q²·⁵) Pauli scaling.
4. **Paper 14** (G4): the headline result (intra-group prereq).
5. **Paper 18** (G3): exchange-constant taxonomy classifying transcendentals in Pauli counts.

### Prereq for Group 5 (QED / gauge)

1. **Paper 7** (G3): S³ proof.
2. **Paper 18** (G3): exchange-constant taxonomy + master Mellin engine (essential).
3. **Paper 2** (G5): α conjecture — the structural motivation (intra-group prereq).
4. **Paper 22** (G3): angular sparsity theorem — needed for QED on S³ Pauli counts.
5. **Paper 24** (G3): Bargmann-Segal / HO comparison; Coulomb/HO asymmetry frames Wilson lattice gauge scope.

### Prereq for Group 6 (precision / observations)

1. **Paper 7** (G3): S³ Fock projection.
2. **Paper 18** (G3): exchange-constant taxonomy — Paper 34 explicitly classifies projections using Paper 18's tier system.
3. **Paper 31** (G3): universal/Coulomb partition.
4. **Paper 24** (G3): Coulomb/HO asymmetry.
5. **Paper 32** (G1): GeoVac spectral triple — for Paper 34's §VIII open questions and for the master Mellin engine / case-exhaustion connection.

---

## Section 6: Load-Bearing Dependency Chains

Five chains that span multiple groups and are critical for the framework's coherence.

### Chain 1: WH1 PROVEN closure

Paper 0 (packing) → Paper 7 (S³ proof) → Paper 18 (exchange constants) → Paper 32 (GeoVac spectral triple) → Paper 38 (SU(2) propinquity convergence) → Paper 40 (unified GH on compact Lie groups) → Paper 42 (Tomita-Takesaki / four-witness Wick rotation) → Paper 43 (Lorentzian extension at (3,1)) → Paper 45 (K⁺-weak-form Lorentzian propinquity).

**Span:** G3 → G3 → G3 → G1 → G1 → G1 → G1 → G1 → G1.
**Load-bearing because:** each math.OA paper relies on the previous one's theorem statement. If Paper 7's S³ proof or Paper 18's master Mellin engine identification is wrong, every downstream propinquity result loses its substrate. The 4/π rate constant in Papers 38 / 40 IS Paper 18 §III.7's M1 Hopf-base-measure signature.

### Chain 2: K = π(B+F−Δ) three-sector decomposition

Paper 7 (S³ Casimir) → Paper 28 (QED on S³, Dirac mode count g_3^Dirac = 40 = Δ⁻¹) → Paper 18 (exchange-constant taxonomy as classification grid) → Paper 32 §VIII (case-exhaustion theorem: every π in any finite chain is M1/M2/M3) → Paper 2 (α observation, conjectural).

**Span:** G3 → G5 → G3 → G1 → G5.
**Load-bearing because:** Paper 2's combination rule K is the framework's most-discussed open question; the structural decomposition is split across G3 (taxonomy), G5 (QED/Dirac structure), and G1 (case-exhaustion theorem). Paper 2's status as conjecture rests on this triangulation. If Paper 28's structural-zero theorem or Paper 32's case-exhaustion theorem changes, the K combination rule's interpretation changes immediately.

### Chain 3: Composed-architecture quantum encoding

Paper 7 (S³) → Paper 13 (hyperspherical) → Paper 15 (Level 4) → Paper 17 (composed geometries) → Paper 22 (angular sparsity theorem) → Paper 14 (qubit encoding, O(Q²·⁵) Pauli scaling) → Paper 23 (nuclear extension) → Paper 20 (resource benchmarks).

**Span:** G3 → G2 → G2 → G2 → G3 → G4 → G4 → G4.
**Load-bearing because:** the quantum-computing pitch (Paper 14 headline) is anchored in five upstream chemistry/foundation results. The angular-sparsity claim in Paper 14 cannot be evaluated without Paper 22; Paper 22 cannot be evaluated without Paper 17's composed architecture; Paper 17 requires Paper 15's molecule-frame solver; Paper 15 sits on Paper 13's hyperspherical formulation; Paper 13 sits on Paper 7's S³ proof.

### Chain 4: Precision-catalogue dictionary

Paper 7 (S³) → Paper 18 (transcendental taxonomy) → Paper 31 (universal/Coulomb partition) → Paper 32 §VIII (case-exhaustion theorem) → Paper 34 (28-projection taxonomy) → Paper 35 (observation-vs-rest-mass split, π enters at exactly one step) → Paper 36 (bound-state QED Lamb shift one-loop closure, Layer-2 inputs).

**Span:** G3 → G3 → G3 → G1 → G6 → G6 → G5.
**Load-bearing because:** the Paper 34 catalogue's §V/§V.B/§V.C/§V.D format depends on the master Mellin engine + universal/Coulomb partition + case-exhaustion theorem being formally stated upstream. Paper 35's Falsifiable Prediction 1 ("π enters iff continuous integration over a temporal/spectral parameter") is a sharpened restatement of Paper 32 §VIII's case-exhaustion theorem in observable language.

### Chain 5: Wilson lattice gauge across SM groups

Paper 7 (S³ = SU(2)) → Paper 24 (Bargmann-Segal S⁵, HO rigidity, four-layer asymmetry) → Paper 25 (Hopf gauge structure, U(1) Wilson on S³) → Paper 30 (SU(2) Wilson on S³ = SU(2), Cartan-subalgebra reduction to Paper 25) → Paper 32 §VIII.B (unified U(1)×SU(2)×SU(3) on three sub-manifolds with 1/(4·N_c) coefficient) → Paper 41 (Rule B Wilson U(1), seven-witness 3D compact U(1) closure on Dirac-S³ graph).

**Span:** G3 → G3 → G5 → G5 → G1 → G5.
**Load-bearing because:** the SM-gauge-content-saturated claim (Bertrand × Hopf-tower truncation at n ≤ 3) requires all three Wilson constructions to exist on GeoVac sub-manifolds. Paper 32 §VIII.B is the synthesis claim; Papers 25, 30, 41 are the gauge-construction primitives. If any one of those constructions falls (e.g., if SU(3) on S⁵ Bargmann is structurally blocked beyond the gauge level — Sprint ST-SU3 verdict), the synthesis sharpens to "gauge YES, matter NO" but doesn't break.

---

## Section 7: Recommendations

### R1. Paper 38 should explicitly cite Paper 18 (master Mellin engine)

Paper 38 cites Papers 2, 7, 24, 32, 40 but not Paper 18. The 4/π rate constant IS Paper 18 §III.7's M1 Hopf-base-measure signature (recorded in CLAUDE.md §1.7 WH1 and memory `l2_quantitative_rate_4_over_pi.md`). The structural identification is currently implicit; making the citation explicit would tighten the math.OA → foundations connection for math.OA readers who don't already know Paper 18.

### R2. Paper 14 (qubit encoding) is undercited from G6

Paper 14 has 22 inbound citations across the series but the G6 catalogue papers (34, 35) cite it only via prose, not bibitem in some places. Paper 34's projection taxonomy explicitly references qubit encoding at multiple projections (§III.6 spectral action, §III.21 multipole/Gaunt, §III.23 symmetry/Young) — bibitem-level citations would help. Currently the G6 → G4 edge count is 7 across 4 papers, modest given Paper 14's centrality.

### R3. Paper 24's four-layer Coulomb/HO asymmetry deserves a load-bearing citation from Paper 39 / Paper 45

Papers 39 and 45 both flag cross-manifold $\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$ as out of scope (G3 = Paper 24 §V four-layer asymmetry), but Paper 39's cite is only contextual. A structural citation pointing to the specific layer (i)-(iv) being blocked would tighten the scope statement and avoid the "out of scope" reading scanning as hand-wave.

### R4. Paper 23 (nuclear shell) is a hidden hub

Paper 23 has 11 inbound citations — second-most in Group 4 — and is cited from every other group (G1, G3, G5, G6). Its role as a universal/Coulomb-partition test case (alongside Paper 24's HO rigidity) is broader than the "nuclear application" framing suggests. Consider promoting Paper 23 in the §6 inventory's loading guide from "On-topic" to "Always" given its cross-group reach.

### R5. Paper 21 (synthesis) is underused

Paper 21 has only 1 inbound citation (from Paper 25) despite being the framework's synthesis paper. The synthesis content (S³ proof chain, N-electron generalization, exchange constant taxonomy) is duplicated and extended in Papers 31 and 32 — Paper 21 may be functionally archived without being formally moved. Consider whether Paper 21 should be (a) archived to `papers/archive/`, (b) folded into Paper 31's universal/Coulomb partition, or (c) cited more explicitly from G5/G6 papers as the entry-point synthesis. Current state has it lingering.

### R6. The synthesis/ folder is empty

`papers/synthesis/` exists but contains no files. Either the folder is for future use (in which case a placeholder README would help), or Paper 21 belongs there (currently sits in `group3_foundations/`), or the folder should be deleted.

### R7. Two-way edges from chemistry → math.OA / observations are zero

G2 → G1, G2 → G5, G2 → G6 are all zero. Chemistry papers were written first and don't refer to the math.OA / QED / observations program; this is historically defensible but creates an artificially one-way dependency. If a future chemistry paper extends Paper 17 to second-row chemistry (currently flagged with W1c-residual wall), it should explicitly cite Paper 34 §V.C autopsies and Paper 32 §VIII spectral-triple framing to close the loop. The G2 ← G1/G5/G6 directions have substantial traffic; making G2 → G1/G5/G6 nonzero would tighten the series.

### Bibliographic surprises (no recommendation, structural notes only)

- **Paper 7 is cited by 31 of 39 papers** — the most universal dependency in the framework, matching its "Always" tier in CLAUDE.md §6. No surprise.
- **Paper 18 (exchange constants) is the second-most-cited paper at 22 inbound**, despite sitting in foundations and not being the headline result. Its taxonomy function is genuinely load-bearing.
- **Paper 0 (packing) is third at 20 inbound** — the axiomatic role is well-recognized.
- **Paper 45 has zero inbound** because it is the most recent (2026-05-18); same for Paper 41 (zero). These are leaf-of-DAG papers, as expected.
- **Paper 6 (Quantum Dynamics) has 1 inbound and 1 outbound** — it is essentially an orphan in the citation network. Its content (O(V) dynamics, Rabi/AIMD/spectroscopy) is computational machinery that other papers don't reach for. Either it serves a standalone purpose for the dynamics community or it should be considered for archive review.

---

**End of memo.**
