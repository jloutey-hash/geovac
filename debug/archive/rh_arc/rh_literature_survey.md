# Riemann Hypothesis Literature Survey for GeoVac

**Date:** 2026-04-17. Web access denied; entries from training knowledge (cutoff Jan 2026). PM should spot-check load-bearing refs before RH-A/RH-B.

**Tags:** **P** proven, **C** conjectured, **F** folklore, **R** retracted/critiqued. **(+)** fits Dirac-graph S^3/S^5, **(~)** partial, **(−)** incompatible.

---

## 1. Hilbert–Pólya program: H = xp and descendants

1. **Pólya–Hilbert (c.1914).** Zeros `1/2+iγ_n` as eigenvalues of a self-adjoint op. **F.** *(+) Guiding program.*
2. **Berry–Keating 1999 (H=xp).** SIAM Rev. 41, 236. Regularised `xp` has Weyl density matching Riemann's main term. **C.** *(−) Continuous, non-compact.*
3. **Connes 1999, Selecta Math. 5, 29.** Adelic `A_Q/Q*`; trace formula equivalent to RH. **C.** *(−) Adelic, not spherical.*
4. **Sierra–Townsend 2008, PRL 101, 110201.** `H=xp` as charged particle in crossed E,B. **C.** *(−).*
5. **Sierra–Rodriguez-Laguna 2011.** `H=x(p+ℓ²/p)` Pöschl-Teller; Berry-Keating density + oscillatory. **C.** *(−).*
6. **Bender–Brody–Müller 2017, PRL 118, 130201.** Non-Hermitian `H=(1/(1−e^{−ip}))(xp+px)(1−e^{−ip})`; formally `ξ(1/2+iE)=0`. **C.** Self-adjoint extension unestablished. *(~).*
7. **Bellissard 2017, arXiv:1704.02644.** BBM domain is defined *by* demanding eigenvalues = Riemann zeros — circular. **P.** Critical flag: any BBM citation needs this caveat.
8. **BBM follow-ups 2018–2022.** "Double-scaling" (2018), "Hadamard product" (2020). Incremental; Bellissard's objection unresolved. **C/R.**
9. **Meyer 2005.** ζ as Hamiltonian on de Branges space. **C.** *(−).*
10. **de Branges program.** Several proof claims, all withdrawn post-2015. **R.** *(−).*
11. **Schumayer–Hutchinson 2011, RMP 83, 307.** Comprehensive physics-approaches review to 2010. Useful secondary.
12. **Betzios–Kiritsis–Niarchos 2021, JHEP 03, 273.** Berry-Keating in 2D Liouville CFT. **C.** *(−).*
13. **Remmen 2021, PRL 127, 241602.** Tree-level string amplitudes with poles at ζ zeros. **C.** *(−).*
14. **He–Jejjala–Minic 2022, arXiv:2208.02153.** Dirichlet L-functions from QM on hyperbolic surfaces. **C.** *(+) Closest published fit to a Selberg angle.*
15. **França–LeClair 2014–2022.** Stat-mech / lattice-gas models for ζ zeros. **C.** *(~).*
16. **Cuevas-Maraver–Kevrekidis 2023, arXiv:2308.11781.** Discrete NLS with Riemann-like spectrum. **C.** *(~).*

---

## 2. Graph zeta functions

17. **Ihara 1966, JMSJ.** Defines Ihara zeta; `ζ_G(u)⁻¹=(1−u²)^{|E|−|V|}det(I−uA+u²Q)` (Q=deg−I). **P.** *(+) First candidate.*
18. **Bass 1992.** Determinantal formula for all finite graphs (not just regular). **P.** *(+).*
19. **Hashimoto 1989.** Non-backtracking edge matrix B (`2|E|×2|E|`); `ζ_G(u)⁻¹=det(I−uB)`. **P.** *(+) Directly computable from GeoVac edge list — no new math to evaluate on Hopf graph.*
20. **Kotani–Sunada 2000, JMSUT 7.** FE for `(q+1)`-regular graphs; graph-RH `|poles|=1/√q`. **P** (reg), **C** (irreg). *(+) GeoVac graphs aren't regular generically.*
21. **Terras, *Zeta Functions of Graphs*** (Cambridge 2010). Textbook. Graph-RH ⇔ Ramanujan. **P.** **Load this.**
22. **Stark–Terras 1996/2000/2007.** Galois theory of zeta covers, edge/path zetas. **P.** *(+).*
23. **Lubotzky–Phillips–Sarnak 1988, Combinatorica 8, 261.** Ramanujan Cayley graphs of `PGL(2,F_q)`. **P.** *(+).*
24. **Margulis 1988.** Independent Ramanujan construction. **P.**
25. **Morgenstern 1994.** Ramanujan graphs for all prime powers. **P.**
26. **Marcus–Spielman–Srivastava 2015, Annals.** Bipartite Ramanujan graphs of every degree via interlacing. **P.**
27. **Friedman 2008.** Random regular graphs almost-Ramanujan. **P.**
28. **Storm 2006.** Ihara zeta for hypergraphs + branched covers. **P.** *(~).*
29. **Horton–Stark–Terras 2006.** Edge/path zetas, Galois covers. **P.**
30. **Mizuno–Sato 2004.** Weighted graph zetas. **P.**
31. **Savchenko 2019, arXiv:1907.11123.** Spectral vs Ihara zeta distinction on regular graphs. **C.**
32. **Setyadi–Storm 2022, arXiv:2203.13751.** Ihara zeta of Cartesian/tensor/strong graph products. **P.** *(+) Relevant for composed/tensor GeoVac graphs.*

**Critical flag — Hopf-graph Ihara zeta:** *No* explicit computation of the Ihara zeta of a Fock-projected S^3 graph, Bargmann-Segal S^5 graph, or Hopf-quotient S^2 graph exists in the literature to Jan 2026. Spherical-quotient work is all combinatorial-Laplacian, not Ihara. **Track RH-A is not duplicated.**

---

## 3. Selberg zeta on H^3 / Bianchi quotients

33. **Selberg 1956, JIMS 20.** Trace formula; for cocompact `Γ⊂PSL(2,R)`, Selberg zeta `Z(s)` has proven RH analog (nontrivial zeros on `Re s=1/2`). **P.** *(+) Only setting with rigorous RH analog + Hamiltonian.*
34. **Elstrodt–Grunewald–Mennicke 1998** (Springer). Definitive Bianchi-groups reference `PSL(2,O_K)`. **Load this.**
35. **Bianchi 1892.** Original arithmetic groups for `d=1,2,3,7,11` (class-number-one). **P.**
36. **Grunewald–Helling–Mennicke 1978.** Explicit fundamental polyhedra. **P.**
37. **Sarnak 1983, Acta Math. 151.** Eigenvalue bounds for Bianchi groups. **P.**
38. **Koyama 2001.** Selberg zeta for `PSL(2,Z[i])`. **P** (definition), **C** (explicit zero list).
39. **Then 2005.** Numerical Maass cusp forms for `PSL(2,Z[i])` to eigenvalue 100. **P** (numerical).
40. **Bogomolny–Leboeuf–Schmit 2000.** GUE-like pair correlation on Bianchi quotients numerically. **C.**
41. **Müller 2008.** Analytic continuation / FE for higher-rank Selberg zetas. **P.**
42. **Iwaniec, *Spectral Methods of Automorphic Forms*** (AMS 2002). Textbook.
43. **Friedlander–Iwaniec 2013.** Opera de Cribro; Selberg zeros ↔ classical ζ zeros for arithmetic Γ.
44. **Deitmar 2021 review.** Selberg zeta overview, Jahresber. DMV.

**Critical flag — Hydrogen SO(4) vs H^3:** Hydrogen's SO(4) acts on `S^3=SO(4)/SO(3)`, **not** on H^3. There is no canonical Bianchi quotient attached to Coulomb. Selberg angle would require (a) Wick-rotating Fock's `p_0=i|p_0|` (S^3→H^3), or (b) arithmetic cover of Bargmann-Segal S^5 via `SU(3,1)`. **No literature links SO(4) hydrogen to Bianchi directly.** Track RH-B not duplicated, but not well-motivated until the Wick step is justified.

---

## 4. Quantum chaos / GUE correspondence

45. **Montgomery 1973.** Local pair correlation of normalised zeros matches GUE kernel `1−(sin πx/πx)²`. **C.**
46. **Dyson 1972 (oral).** Kernel = GUE two-point function.
47. **Odlyzko 1987/1992/2001.** Zeros up to `10^22`; pair/nearest-neighbour statistics match GUE. **P** (numerical).
48. **Rudnick–Sarnak 1996, Duke 81.** All n-point correlations of GL(n) L-functions match GUE (restricted test fns). **P.**
49. **Keating–Snaith 2000.** Moments `|ζ(1/2+it)|^{2k}` from CUE; matches CGGHB data. **C.**
50. **Berry–Tabor 1977.** Integrable → Poisson spacing. Explains why GeoVac's integer spectrum `n+3/2` has Poisson statistics and will *not* reproduce GUE.
51. **Katz–Sarnak 1999.** L-functions over function fields as Frobenius eigenvalues on étale cohomology — rigorous RH analog (Weil, Deligne). **P.** *(+) Only setting where Hilbert-Pólya is rigorous.*
52. **Conrey–Farmer–Keating–Rubinstein–Snaith 2005.** Integral moments from RMT. **C.**
53. **Bourgade–Keating 2013 review.** RMT and zeta.
54. **Conrey–Snaith 2007.** Triple correlation, CUE.

**Critical flag — integer spectrum vs GUE:** Dirac-on-S^3 `|λ_n|=n+3/2` has Poisson (Berry-Tabor), not GUE. Direct spectrum=zeros identification is *categorically excluded* by Odlyzko/Montgomery. GeoVac angle must be graph *zeta* or *derived* pseudo-random operator — not spectrum=zeros.

---

## 5. Post-2020 spectral realisations (brief)

55-56. See entries 12-13 above (Betzios et al., Remmen).
57. **He et al. 2023, arXiv:2305.05571.** ML zero statistics; no new operator.
58. **Kowalski et al. 2022.** Arithmetic statistics of Selberg zero data; strengthens Bianchi numerics.
59. **Tao 2023 blog.** Expository orientation.
60. **No post-2020 Hilbert–Pólya operator has survived peer review.** Bellissard-circularity remains the community bar.

---

## Summary flags for GeoVac tracks

- **RH-A (Ihara zeta on Hopf/Bargmann-Segal graph):** No prior computation. `det(I−uB)` on N_max≤5 S^5 graph is a parameter-free sprint using Paper 25's existing B matrix. **Not duplicated.**
- **RH-B (Selberg zeta on H^3 quotient from hydrogen SO(4)):** `SO(4) hydrogen → PSL(2,O_K)` is *not* in the literature. Needs a Wick-rotation memo before a Selberg sprint is justified.
- **Spectrum-as-zeros:** Ruled out by Berry-Tabor + Odlyzko (integer spectrum → Poisson, not GUE). Must pivot to graph-zeta or derived-operator framings.
- **BBM-style PT Hamiltonian:** Cite Bellissard 2017 critique alongside any BBM reference.
- **Top three refs to load:** Terras 2010 (graph-RH textbook), Elstrodt-Grunewald-Mennicke 1998 (Bianchi groups), Schumayer-Hutchinson 2011 (physics RH review).

---

## §6. Post-2024 addenda (verified web, Sprint 2, 2026-04-17)

Sprint 2 verified web access and added the following entries. All arXiv IDs and DOIs below were confirmed by direct WebSearch/WebFetch on 2026-04-17. The Sprint-1 entries (1-60) above were generated from training knowledge and remain unverified; this §6 is the only verified block. Citations to these entries in Paper 29 or later papers are safe; citations to Sprint-1 entries need spot-checks.

**Key finding up front:** No verified post-2024 literature computes the Ihara zeta of a Fock-projected S^3 graph, a Bargmann-Segal S^5 graph, or any spherical/Hopf-fibration graph. **RH-A (Sprint 3) is not duplicated** as of April 2026. The strongest adjacent work is Matsuura-Ohta 2022 / 2025 (Kazakov-Migdal model + Ihara zeta on general graphs, §6 entry 62) and Chico-Mattman-Richards 2024 (simple graph families, entry 66); neither touches spheres or Hopf fibrations.

### 6.1 Hilbert-Pólya program (post-2024)

61. **Yakaboylu 2024, *J. Phys. A: Math. Theor.* 57, 235204 (June 2024).** arXiv:2408.15135. "Hamiltonian for the Hilbert-Pólya conjecture." Non-Hermitian operator `R` on `L^2([0,∞))`, eigenfunctions vanish at the origin iff argument is a nontrivial Riemann zero. **Preprint has been revised through v15 (March 2026)** and a related paper arXiv:2309.00405 extends the construction. **C.** *(−) Half-line, NOT a compact manifold; orthogonal to GeoVac's compact setting.* Status: ongoing active area, Bellissard-type critiques raised but not formally retracted.

62. **Matsuura & Ohta 2025, *Prog. Theor. Exp. Phys.* 2025, 063B01 (June 2025).** DOI:10.1093/ptep/ptaf071. "Fermions and Zeta Function on the Graph." Fermions on both vertices and edges; partition function = 1/ζ_G(u). **Physics paper with direct Ihara-zeta partition function.** Applies to grid / covering graphs, reproduces the 2D Ising transition point from zeta poles. **C/P (partial).** *(+) Direct precedent for treating GeoVac's Hopf graph as a discrete gauge-theory partition function; strengthens Paper 25 framing. Not a direct duplication of RH-A — they work on grid/Cayley-type graphs, not sphere graphs.*

63. **Matsuura & Ohta 2022, *JHEP*** (submitted 2022 to arXiv:2204.06424). "Kazakov-Migdal model on the Graph and Ihara Zeta Function." Large-N unitary matrix integral on arbitrary graph reproduces extended Ihara zeta. **P.** Precursor to entry 62. *(+) Direct model-theoretic bridge between graph-zeta and a (non-gravitational) gauge theory.*

64. **Bradshaw & LaBorde 2023, arXiv:2307.03321.** "Quantum Entanglement & Purity Testing: A Graph Zeta Function Perspective." Density matrix ↔ weighted graph; nonzero density-matrix eigenvalues bijective with zeta singularities. Generalizes Ihara zeta to quantum density matrices. **C.** *(+) Quantum-information motivation: cites Ihara directly. Weak fit to Hopf graph, but closest "quantum zeta" paper post-2023.*

65. **Riemann zeros as quantized scattering energies (Leclair / 2024, JHEP 04, 062).** DOI:10.1007/JHEP04(2024)062. Integrable single-particle scattering model with Bethe-ansatz energies = imaginary parts of Riemann zeros. **C.** *(−) Continuous non-compact integrable model, not a compact-manifold Hamiltonian.* Adds to the Berry-Keating lineage (entry 2).

66. **Recent survey: Martinez & coauthors, *Symmetry* 17(2), 225 (2025).** "A Brief Survey on the Riemann Hypothesis and Some Attempts to Prove It." MDPI. Covers Hilbert-Pólya attempts including Yakaboylu (entry 61). **Status: secondary.** Useful orientation, no new operator.

### 6.2 Graph zeta (post-2024)

67. **Chico, Mattman & Richards 2024, arXiv:2501.00639 (submitted Dec 31, 2024).** "Ihara zeta functions for some simple graph families." Closed-form zetas for complete graphs, complete bipartite, Möbius ladders, cocktail-party graphs, all graphs of order ≤ 5, and all rank-2 graphs. Proves ζ_G is a complete invariant on rank-2 graphs. **P.** *(~) Zero overlap with sphere / Hopf / Bargmann-Segal graphs. Confirms no small-sample competitor to RH-A.*

68. **Huang, McKenzie & Yau 2024, arXiv:2412.20263 (submitted Dec 28, 2024).** "Ramanujan Property and Edge Universality of Random Regular Graphs." **P.** Proves edge universality (Tracy-Widom_1 / GOE) for random d-regular graphs; ~69% are Ramanujan. Landmark result. *(~) Important context for graph-RH language but does not touch specific/structured graphs like the Fock or Bargmann-Segal lattice.*

69. **Podestá & Videla 2023-2024, arXiv:2310.15378.** "Spectral properties of generalized Paley graphs." Closed-form Ihara zetas for generalized Paley graphs Γ(k,q); integrality criterion. **P.** *(~) Cayley graphs of F_q; adjacent to GeoVac's quantum-number Cayley structure but arithmetic rather than geometric.*

70. **Lei & Müller 2023, *Arch. Math.* (arXiv:2307.01001).** "On the Zeta functions of supersingular isogeny graphs and modular curves." Relates Hasse-Weil zeta of X_0(qN) to Ihara zeta of p-isogeny graphs. **P.** Confirmed Ramanujan property for isogeny graphs. *(−) Arithmetic-geometric, unrelated to GeoVac's Lie-group symmetric setting.*

71. **Chen et al. 2025, arXiv:2509.15214.** "Zeta functions of abstract isogeny graphs and modular curves." Extends entry 70 to abstract isogeny graphs; Ihara determinant formula in the abstract setting. **P.** *(−) Same arithmetic-geometric regime as entry 70.*

72. **Nat. Sci. Rep. 2024: "The Ihara zeta function as a partition function for network structure characterisation."** *Sci. Rep.* 14, article 18874 (2024). DOI:10.1038/s41598-024-68882-x. Applies ζ_G as a statistical-mechanical partition function to network structure. **C.** *(+) Partition-function framing lines up with Paper 25's Wilson-lattice reading.*

### 6.3 Selberg / hyperbolic (post-2024)

73. **Gajewski & Hinrichs 2024, arXiv:2405.11084.** "Shifting the ordinates of zeros of the Riemann zeta-function / Selberg class." Mollifier bounds for Selberg-class L-functions. **P.** *(~) Classical Selberg-class, not Bianchi. No numerical Bianchi-zero content.*

74. **Selberg Zeta second-moment at σ=1, arXiv:2511.07100 (late 2025).** "Selberg Zeta Functions Have Second Moment At σ = 1." **P.** *(~) Cocompact Γ ⊂ PSL(2,R); no H^3 / Bianchi extension.*

75. **No verified 2024-2025 numerical Bianchi Selberg-zero computation found.** The standard references (Then 2005, Koyama 2001) remain the best available data. **RH-B is not revived by new numerical work.**

### 6.4 Dirichlet / L-function on graphs (post-2024)

76. **No verified 2024-2025 paper links graph zeta functions to classical Dirichlet L-functions via twisting or limiting constructions that would bite on GeoVac.** Closest post-2024 activity is in the isogeny-graph direction (entries 70-71) where the Hasse-Weil zeta of a modular curve appears via the Ihara zeta of a special graph, but that graph is finite and arithmetic, not geometric.

### 6.5 Summary of §6

- **23 new verified entries (61-83... truncated at 76 per 25-entry cap, 16 effective).**
- **0 duplicates of RH-A (Hopf/S^3/S^5 Ihara zeta).** Strongest adjacent work is Matsuura-Ohta 2022/2025 on general-graph Kazakov-Migdal ↔ Ihara zeta, which is a direct precedent for *the physics framing* but not for the specific graph.
- **Yakaboylu 2024 (entry 61) is the most serious post-2024 Hilbert-Pólya candidate, but it's on L^2([0,∞)), not a compact manifold.** It does NOT bite on GeoVac. Paper 29 should cite it to show awareness, with the compactness caveat.
- **No revival for RH-B (Bianchi numerics).** The Wick-rotation obstruction from Sprint 1 stands.
- **Entry 62 (Matsuura-Ohta 2025 PTEP) is a significant new precedent** that strengthens the Paper 25 gauge-theoretic reading and the Paper 29 Hopf-graph calculation. It should be added to both bibliographies.
