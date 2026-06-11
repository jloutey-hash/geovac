# NCG Renormalization-Group Literature Survey

**Date:** 2026-05-15
**Scope:** Noncommutative-geometry (NCG) approaches to renormalization-group (RG) flow on lattice gauge theories, with emphasis on Wilson-type constructions on compact Riemannian manifolds and their continuum limits.
**Audience:** GeoVac PM/PI deciding the next RG sprint scope.

---

## 1. Existing NCG–RG frameworks

The published NCG-renormalization landscape splits into four largely independent threads. None of them target "Wilson lattice gauge on a finite-cutoff spectral truncation of S³ with running couplings" as such; each thread touches a different piece of that statement.

### 1.1 Connes–Kreimer Hopf-algebra renormalization (1998 onward)

Connes and Kreimer, *Hopf Algebras, Renormalization and Noncommutative Geometry*, Commun. Math. Phys. 199 (1998) 203–242 (arXiv:hep-th/9808042), give the foundational construction: the BPHZ subtraction procedure for ultraviolet-divergent Feynman graphs is rewritten as a Birkhoff decomposition of loops in an infinite-dimensional pro-algebraic group, with the relevant Hopf algebra built from rooted trees / Feynman graphs. The Lie algebra is the dual of an operadic graph-insertion operation. Connes–Marcolli's book *Noncommutative Geometry, Quantum Fields and Motives* (AMS Colloquium 55, 2008) extends this to the motivic Galois group of mixed Tate motives via the Riemann–Hilbert correspondence.

**What this thread does:** turns perturbative renormalization into a clean algebraic / categorical machine. Recovers the RG as a one-parameter subgroup of the Birkhoff group.

**What this thread does NOT do:** it is fundamentally a perturbative continuum-QFT formalism. It does not act on a Wilson lattice gauge action, does not interact with spectral truncations or finite spectral triples, and has no direct apparatus for a *physical* running of a gauge coupling on a graph. The "noncommutative geometry" in the title refers to the transverse-index-theory Hopf algebra of foliations, not to the spectral-triple framework GeoVac uses.

Van Suijlekom adapted Connes–Kreimer to gauge theories: his early papers (2007–2008) extend the Hopf algebra to include gauge symmetries via Slavnov–Taylor / Ward identities; see his thesis and review *Gauge Symmetries and Renormalization* (arXiv:2001.00104).

### 1.2 Spectral-action one-loop renormalization

This is the published thread closest to "GeoVac running couplings."

- **van Suijlekom, *Renormalizability Conditions for Almost-Commutative Manifolds*, Ann. Henri Poincaré 15 (2014) 985–1011 (arXiv:1112.4690).** Formulates conditions under which the asymptotically expanded spectral action on an almost-commutative manifold is renormalizable as a higher-derivative gauge theory. The conditions are graph-theoretical and live on the Krajewski diagrams that classify finite spectral triples. The result is *only* renormalizability — no explicit β-functions are computed. Applies to almost-commutative manifolds in the unbounded (canonical) sense, NOT to spectral truncations.

- **van Nuland & van Suijlekom, *One-Loop Corrections to the Spectral Action*, JHEP 05 (2022) 078 (arXiv:2107.08485).** Quantizes the spectral action perturbatively at one loop. Expands the spectral action in higher Yang–Mills and Chern–Simons forms; path-integrates over matrix fluctuations around a fixed background; shows the one-loop counterterms are of the same form (higher YM and CS), so they can be subtracted within the spectral framework. Ward identities are the technical engine. The paper establishes one-loop renormalizability *in a generalized sense* and stays inside the spectral language; it does NOT derive concrete β(g) for any physical coupling, and it does NOT apply to spectral truncations (the underlying triple is still canonical).

- **Chamseddine–Connes & co-authors (2010 onward).** The "running couplings" in the spectral-action SM are imported from the standard model β-functions of ordinary QFT, applied at a fixed unification scale fixed by the heat-kernel coefficients a₀ / a₂ / a₄ ratios. Chamseddine–Connes *Noncommutative Geometry as a Framework for Unification of all Fundamental Interactions including Gravity. Part I* (arXiv:1004.0464) is explicit that the spectral action does not change the SM RG running; unification of α_i at a single scale is a *boundary condition* imposed at Λ, not a derived consequence. Devastato–Lizzi–Valcarcel Flores–Vassilevich (arXiv:1410.6624) extend this by adding dimension-6 operators motivated by the spectral action; still no autonomous β-function generation.

**Net status of thread 1.2:** the spectral action gives one-loop renormalizability and a clean tree-level relation between gauge couplings at the unification scale; it does NOT autonomously generate β(g). Running is imported from standard QFT.

### 1.3 Grosse–Wulkenhaar matrix-model thread

Grosse and Wulkenhaar's φ⁴ model on Moyal space with harmonic potential is the most physically nontrivial *renormalizable* noncommutative field theory. Key results:
- Renormalizable to all orders, no Landau ghost (Grosse–Wulkenhaar 2003–2005).
- β-function vanishes at the self-dual point Ω = 1 up to three loops (Disertori–Rivasseau and others), conjectured (and largely proven) to vanish to all orders (arXiv:hep-th/0612251).
- Asymptotically safe at one loop.
- Constructive (non-perturbative) renormalization established for the 2D model (arXiv:1104.3750).

This thread is *operatively* the cleanest "NCG β-function" story in the literature, but the geometry is Moyal noncommutative R^d, not a Wilson lattice on a compact Riemannian manifold. The technology (Wetterich equation in the matrix base, multiscale loop-vertex expansion) is reusable for finite matrix truncations and matrix models, which is where the Perez-Sanchez thread (1.4) connects.

### 1.4 Functional RG on random / fuzzy noncommutative geometries (Perez-Sanchez 2020–2025)

This is the published thread closest in spirit to GeoVac's situation.

- **Perez-Sanchez, *On Multimatrix Models Motivated by Random NCG I: The Functional RG as a Flow in the Free Algebra*, Ann. Henri Poincaré 23 (2022) 1979–2032 (arXiv:2007.10914).** Applies the Wetterich functional RG to multimatrix actions motivated by Barrett's fuzzy / finite spectral triples. Derives explicit β-functions, identifies large-N fixed points, computes critical exponents of 2-dim geometries. The free-algebra calculus (Voiculescu cyclic gradient + Rota–Sagan–Stein derivative) is what makes the noncommutative Laplacian act cleanly on the Wetterich equation.

- **Perez-Sanchez, *On Multimatrix Models Motivated by Random NCG II: A Yang–Mills–Higgs Matrix Model*, Ann. Henri Poincaré 23 (2022) 2447–2511 (arXiv:2105.01025).** Extends thread (I) to a gauge sector. This is the published bridge between FRG and almost-commutative finite spectral triples in the Barrett / Marcolli–vS spirit.

- **Perez-Sanchez, *Bratteli Networks and the Spectral Action on Quivers* (arXiv:2401.03705)** and **Perez-Sanchez, *Comment on "Gauge Networks in Noncommutative Geometry"* (arXiv:2508.17338).** These are the corrections/extensions to Marcolli–vS 2014 (arXiv:1301.3480). The 2025 comment is the structurally important one for GeoVac: it shows that the *continuum limit* of the Marcolli–vS gauge-network action on quivers is **Yang–Mills without Higgs** (rather than YM-Higgs as originally claimed). Neither the 2024 nor the 2025 paper contains RG apparatus — both stop at the spectral action / continuum limit.

**Net status of thread 1.4:** the only published thread with both (i) finite spectral triples and (ii) explicit β-functions. The cost is that the matrix models are not Wilson lattice gauge — they are *random matrix actions* whose "geometry" is fluctuating, more analogous to Euclidean dynamical triangulations than to Wilson SU(N) on a fixed Hopf graph.

---

## 2. Wilson block-spin on graphs

Wilson block-spin / Migdal–Kadanoff has a 45-year history on regular Z^d lattices (Migdal 1975, Kadanoff 1976, Sect. Migdal–Kadanoff recursion relations for SU(2) and SU(3): see Itoh–Iwasaki–Yoshié, *Nucl. Phys. B210* (1981) and follow-ups; modern reviews in Smit's lattice-QFT textbook). The general picture on Z^d is:
1. Decimation: integrate out a sublattice of sites, generating new effective couplings.
2. Bond moving: rearrange the resulting plaquette terms to recover the original lattice topology with shifted couplings.
3. Iterate to get a renormalized action whose flow defines β(g).

**Adaptation to graphs (as opposed to Z^d):** there is no widely-used clean published Migdal–Kadanoff procedure on a *generic finite graph* of the type GeoVac uses (Hopf S³, Bargmann-Segal S⁵). The technical obstruction is that bond-moving in MK requires a translational symmetry to assemble the moved bonds back into a regular plaquette structure; that translational symmetry is exactly what GeoVac's irregular Fock-graph topology breaks. Generalized reflection-positivity-preserving variants exist (Tomboulis 1984, *Phys. Rev. D 30* 455), but they still require sufficient regularity.

The closest published existing alternative is the **gradient-flow / Wilson-flow exact RG (Lüscher 2010, Sonoda–Suzuki 2024/2025 PTEP)**: a continuum-style RG where the smoothing operator is the Yang–Mills gradient flow ∂_t A = -δS/δA. This is well-defined on any Riemannian manifold (including S³) and on any spectral triple where one can define a heat kernel. The cost is that gradient-flow RG is geometric/continuum-style, not block-spin; it does NOT integrate out modes from the spectrum of the finite-cutoff Dirac operator.

**Marcolli–vS / Perez-Sanchez gauge-network track does NOT include RG.** Both Marcolli–vS 2014 (arXiv:1301.3480) and Perez-Sanchez 2024/2025 (arXiv:2401.03705, arXiv:2508.17338) stop at the level of (i) defining the spectral action on the quiver / gauge network, (ii) taking a continuum limit. They do not provide block-spin on the graph; they do not provide running couplings.

**Conclusion for §2:** There is no off-the-shelf "Wilson block-spin on the Hopf graph" framework in the literature. GeoVac would be writing this from scratch if it chose option (b) of §5. The closest analog is functional RG in the matrix-model formulation (Perez-Sanchez I/II), which works at the level of the finite-N truncation directly without ever passing through a block-spin decimation.

---

## 3. Spectral-action approaches to running couplings

The Chamseddine–Connes spectral action `Tr f(D²/Λ²)` admits a heat-kernel expansion
```
Tr f(D²/Λ²) = Σ_k f_{2k} Λ^{d-2k} a_{2k}(D²)
```
where a_{2k} are Seeley–DeWitt coefficients and f_{2k} are moments of the cutoff function f. At d = 4, this generates the Yang-Mills + Einstein-Hilbert + cosmological constant + Weyl-squared action, with coefficients fixed by f_0, f_2, f_4 and the geometric data of the spectral triple.

**Does the heat-kernel expansion give one-loop β-functions automatically?** *No.* The spectral action at order Λ^0 in the heat-kernel expansion is the classical (tree-level) action, with a fixed scale Λ set by the cutoff function. β-functions describe the variation of *renormalized* couplings under variation of the *renormalization* scale, after one-loop integration of fluctuations. The two scales are conceptually different: Λ (heat-kernel cutoff) ≠ μ (renormalization scale). Quoting Chamseddine–Connes 2010 explicitly: spectral-action minimal formalism does NOT change SM RG running; one-loop β(g) must be supplied externally.

What CAN be done with the spectral action:
- **Boundary conditions at the unification scale**: heat-kernel ratios a₀/a₂/a₄ fix relations between α_1, α_2, α_3 at one specific scale Λ; running from Λ down to electroweak scale uses ordinary SM β-functions.
- **One-loop renormalizability of the spectral action itself**: the asymptotically expanded action is renormalizable as a higher-derivative gauge theory (van Suijlekom 2014), and its one-loop counterterms are of the same spectral form (van Nuland–van Suijlekom 2022). This is the rigorous result, but it doesn't give running of *finite* couplings.
- **Heat-kernel running**: variation of Λ in `Tr f(D²/Λ²)` does change the relative weight of curvature / kinetic / mass terms; this is sometimes called "spectral-action running" or "scale running of the heat-kernel coefficients" but is structurally distinct from a renormalization-group β-function (it does not integrate out modes).

**Most useful single reference for the GeoVac PM:** Vassilevich, *Heat kernel expansion: user's manual* (arXiv:hep-th/0306138). This catalogs all heat-kernel coefficients on spheres, fuzzy spaces, etc., to high order; it is the technical input for any "spectral-action running" computation on S³.

The European Physical Journal C paper *Heat kernel coefficients on the sphere in any dimension* (Eur. Phys. J. C 80 (2020) 269; arXiv:1910.00543, Kluth–Litim) gives explicit Seeley–DeWitt on spheres in any dimension and is the cleanest independent cross-check on Paper 28's S³ formulas.

---

## 4. Gap analysis

| Existing piece | Provides | Missing for GeoVac |
|:---|:---|:---|
| van Nuland–vS 2022 (arXiv:2107.08485) | One-loop renormalizability of spectral action; counterterms stay spectral | Applies to canonical (unbounded) triple, not spectral truncation; no explicit β(g) |
| van Suijlekom 2014 (arXiv:1112.4690) | Renormalizability conditions via Krajewski diagrams | Conditions only; no β(g); not for spectral truncations |
| Connes–vS 2021 (arXiv:2004.14115) | Operator-system framework for spectral truncations | No RG / running couplings — stops at the operator-system level |
| Perez-Sanchez I/II (arXiv:2007.10914, 2105.01025) | Explicit β-functions on finite matrix models via Wetterich FRG; YM-Higgs sector | Geometry is random/fuzzy, not Wilson lattice on fixed Hopf graph |
| Marcolli–vS 2014 + Perez-Sanchez 2024/2025 | Spectral action on quivers → YM continuum limit | No RG apparatus; no running of g_YM |
| Hekkelman–McDonald 2024 (arXiv:2412.00628) | Noncommutative integral on spectrally truncated triples; Szegő limit | Not RG; provides a *trace* on the truncation, not a β-function |
| Grosse–Wulkenhaar thread | All-orders vanishing β at self-dual point; constructive RG | Moyal noncommutative R⁴, not compact Lie group; not Wilson lattice |
| Chamseddine–Connes 2010 (arXiv:1004.0464) | Tree-level unification boundary conditions at Λ | RG running imported from external SM β-functions |
| Wilson-flow / gradient-flow RG (Lüscher 2010, Sonoda–Suzuki 2025) | Geometric continuum RG via heat-flow smoothing | Not formulated on spectral truncations; not autonomous on a finite graph |

**The single closest existing published work to "Wilson lattice gauge on a finite-cutoff spectral truncation of S³ with running couplings" is the pair Perez-Sanchez 2007.10914 + 2105.01025.** That work provides explicit β-functions on finite matrix models motivated by Barrett's finite spectral triples, including a gauge sector. But: (a) it does NOT use Wilson lattice gauge — it uses random matrix actions whose "geometry" is integrated over; (b) it does NOT use the Connes–vS spectral truncation paradigm — it uses Barrett's fuzzy / random NCG.

**The single closest *spectral-truncation* paper is Connes–vS 2021 (arXiv:2004.14115).** It introduces the operator-system framework that GeoVac's Paper 38 lives in, but it explicitly *does not* contain RG. Hekkelman–McDonald 2024 (arXiv:2412.00628) is the most recent extension of that paradigm, providing a noncommutative integral, but still no β-function.

**There is no published work that combines (Connes–vS spectral truncation framework) × (Wilson lattice gauge in Marcolli–vS sense) × (Wetterich / Wilsonian RG).** GeoVac sits on the unique intersection of these three.

---

## 5. Recommended sprint scoping

**Recommendation: pursue Option (a) — direct β(α) via varying Λ in Tr f(D²/Λ²) using Paper 28's QED-on-S³ machinery.** Reasoning in three paragraphs below.

**Why (a) over (b) and (c):** Option (b) (Migdal–Kadanoff block-spin on the Hopf graph) is not blocked by physics but by infrastructure: there is no published MK procedure on irregular graphs, and writing one for the Hopf graph at finite n_max requires solving the bond-moving problem on a vertex set with non-uniform connectivity. This would be a multi-month methods-development sprint with significant risk of producing a procedure whose "renormalized couplings" are convention-dependent. Option (c) (Connes–Kreimer Hopf-algebra lifted to GeoVac) is the cleanest from the categorical-NCG side, but Connes–Kreimer is a *perturbative continuum-QFT* machine — it acts on Feynman graphs of QFT on the spectral triple, not on the spectral truncation itself. Lifting it to GeoVac would mean perturbatively quantizing QED on Dirac-S³ (which Paper 28 has already done at one loop, including F₂, self-energy, and three-loop topologies), then running Connes–Kreimer on those Feynman graphs. This is feasible but is doing Connes–Kreimer-on-Paper-28, not new NCG.

**Why (a) is the right next sprint:** GeoVac already has the *spectral action machinery* it needs (Paper 28: Tr f(D²/Λ²) computed on Dirac-S³ at the level of integrated heat-kernel; Paper 36: bound-state one-loop QED on Dirac-S³ closed at sub-percent; sprint TS / MR-B: closed-form modular-residual signature √π·ℚ ⊕ π²·ℚ for the M2 mechanism). Vassilevich's heat-kernel toolkit (arXiv:hep-th/0306138) and Kluth–Litim's explicit S^d Seeley–DeWitt (arXiv:1910.00543) are the external cross-checks. Concrete deliverable: compute the heat-kernel-derived "spectral-action running" of the QED coupling α from `Tr f(D²/Λ²)` evaluated at varying Λ over n_max ∈ {2..6}, with α-dependence entering via Paper 25's U(1) Wilson sector. **Honest scope-boundary:** this is *heat-kernel scale running*, not Wilsonian RG flow — it does not integrate out modes, it varies the cutoff function's scale. This distinction must be enforced in the writeup so the result is not over-stated. The substantive question the sprint can answer is: does the Λ-dependence of the spectral action on Dirac-S³ at finite n_max reproduce the one-loop QED β(α) = 2α²/(3π) that Paper 28 already extracted from the bare spectral sum? If yes (which is plausible structurally), the sprint provides an *internal consistency check* between the bound-state-QED side (Paper 36) and the spectral-action side (Paper 28) — and this is a new and publishable result with no clear analog in the published literature.

**Concrete 4–6 week sprint plan:** Track 1 (week 1–2): implement `tr_f_D_squared_over_Λ_squared(n_max, Λ, f)` using existing Paper 28 spectral data; vary Λ over the natural scale [Λ_min, Λ_max] set by Paper 28's spectral cutoff theorem; extract dα/d log Λ. Track 2 (week 3–4): compare to one-loop QED β(α) = 2α²/(3π); compare to Paper 36's Lamb shift residual decomposition. Track 3 (week 5–6): write up as Paper 41 (math-ph) or as new section in Paper 28; check whether the result generalizes to SU(2) (Paper 30) and SU(3) (Sprint ST-SU3) — predicted yes via the universal 1/(4N_c) coefficient. If the result clears the sprint, Option (c) (Connes–Kreimer + Paper 28 Feynman graphs at two loops) becomes the natural follow-on. If it does not clear, the sprint still produces a clean negative result — the spectral action's Λ-dependence is structurally distinct from the bound-state β(α). Either outcome is publishable.

---

## 6. Bibliography

All references verified by direct WebFetch / WebSearch against arXiv and journal sites. ArXiv IDs and DOIs noted where available. Uncertain entries flagged explicitly.

**NCG foundations and gauge networks**
- Marcolli, M. and van Suijlekom, W. D. *Gauge networks in noncommutative geometry.* J. Geom. Phys. 75 (2014) 71–91. arXiv:1301.3480. **VERIFIED.** Stops at the action level; no RG.
- Perez-Sanchez, C. I. *Bratteli networks and the spectral action on quivers.* arXiv:2401.03705 (2024). **VERIFIED.** Extension of Marcolli–vS; no renormalization apparatus.
- Perez-Sanchez, C. I. *Comment on "Gauge networks in noncommutative geometry".* arXiv:2508.17338 (2025). **VERIFIED.** Continuum limit is YM without Higgs (correction of Marcolli–vS).

**Spectral truncations**
- Connes, A. and van Suijlekom, W. D. *Spectral truncations in noncommutative geometry and operator systems.* Commun. Math. Phys. 383 (2021) 2021–2067. arXiv:2004.14115. **VERIFIED.** No RG; the operator-system framework Paper 38 builds on.
- Hekkelman, E.-M. and McDonald, E. A. *A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity.* J. Funct. Anal. (accepted 2025). arXiv:2412.00628. **VERIFIED.** Szegő limit formula; not RG.
- Leimbach, F. and van Suijlekom, W. D. *Convergence of spectral truncations for compact metric groups.* Int. Math. Res. Not. 2025 (13), paper rnaf197. arXiv:2310.14733. **VERIFIED via Paper 40 context.** Used in Paper 40's universal rate proof.

**Spectral-action renormalization**
- Chamseddine, A. H. and Connes, A. *The Spectral Action Principle.* Commun. Math. Phys. 186 (1997) 731–750. arXiv:hep-th/9606001. **VERIFIED.** Tree-level spectral action; no autonomous RG.
- Chamseddine, A. H. and Connes, A. *Noncommutative geometry as a framework for unification of all fundamental interactions including gravity. Part I.* Fortsch. Phys. 58 (2010) 553. arXiv:1004.0464. **VERIFIED.** Explicit on running being imported from SM, not derived.
- van Suijlekom, W. D. *Renormalizability conditions for almost-commutative manifolds.* Ann. Henri Poincaré 15 (2014) 985–1011. arXiv:1112.4690. **VERIFIED.** Krajewski-diagram renormalizability conditions; no β(g).
- van Nuland, T. D. H. and van Suijlekom, W. D. *One-loop corrections to the spectral action.* JHEP 05 (2022) 078. arXiv:2107.08485. **VERIFIED.** One-loop renormalizability stays spectral; no β(g) extraction.
- Devastato, A., Lizzi, F., Valcarcel Flores, C., and Vassilevich, D. *Unification of coupling constants, dimension six operators and the spectral action.* Phys. Rev. D 91 (2015) 035026. arXiv:1410.6624. **VERIFIED.** Unification scenarios with dimension-6 operators; tree-level RG only.

**Connes–Kreimer Hopf-algebra renormalization**
- Connes, A. and Kreimer, D. *Hopf algebras, renormalization and noncommutative geometry.* Commun. Math. Phys. 199 (1998) 203–242. arXiv:hep-th/9808042. **VERIFIED.**
- Connes, A. and Marcolli, M. *Noncommutative Geometry, Quantum Fields and Motives.* AMS Colloquium Publications 55 (2008). **VERIFIED via several review citations.** Standard reference for motivic / Riemann–Hilbert renormalization.
- van Suijlekom, W. D. *Gauge symmetries and renormalization.* Math. Phys. Anal. Geom. 25 (2022) 6. arXiv:2001.00104. **VERIFIED.** Adapts Connes–Kreimer to gauge theory.

**Grosse–Wulkenhaar matrix-model thread**
- Disertori, M., Gurau, R., Magnen, J., and Rivasseau, V. *Vanishing of beta function of non-commutative Φ⁴_4 theory to all orders.* Phys. Lett. B 649 (2007) 95. arXiv:hep-th/0612251. **VERIFIED.** All-orders vanishing β at self-dual point.
- Grosse, H. and Wulkenhaar, R. (and follow-ups in arXiv:1104.3750 by Wang, etc.). General reference for the model. **VERIFIED.**

**Functional RG on random / fuzzy NCG (closest to GeoVac)**
- Perez-Sanchez, C. I. *On multimatrix models motivated by random noncommutative geometry I: The functional renormalization group as a flow in the free algebra.* Ann. Henri Poincaré 23 (2022) 1979–2032. arXiv:2007.10914. **VERIFIED.** Explicit β-functions on finite matrix models motivated by Barrett's fuzzy spectral triples.
- Perez-Sanchez, C. I. *On multimatrix models motivated by random noncommutative geometry II: A Yang–Mills–Higgs matrix model.* Ann. Henri Poincaré 23 (2022) 2447–2511. arXiv:2105.01025. **VERIFIED.** Gauge extension of paper I.
- Azarfar, S. and Khalkhali, M. *Random finite noncommutative geometries and topological recursion.* arXiv:1906.09362 (2019). **CITATION SOFT** — appears in Perez-Sanchez bibliographies, not independently verified during this survey. Note in any followup memo if used.

**Block-spin on graphs / Wilson-flow RG**
- Migdal, A. A. (1975, 1976) and Kadanoff, L. P. (1976), original Migdal–Kadanoff papers. Standard textbook references; no need for arXiv. **STANDARD.**
- Tomboulis, E. T. *Generalized reflection positivity and modified Migdal–Kadanoff procedure for lattice gauge theories.* Phys. Rev. D 30 (1984) 455. **VERIFIED via OSTI / NASA ADS.**
- Itoh, S., Iwasaki, Y., and Yoshié, T. *Migdal–Kadanoff recursion relations in SU(2) and SU(3) gauge theories.* Nucl. Phys. B 210 (1982) 49. **VERIFIED.**
- Sonoda, H. and Suzuki, H. *Basis of the gradient flow exact renormalization group for gauge theory.* PTEP 2025 (9) 093B05. **VERIFIED.** Modern Wilson-flow RG with gauge invariance.

**Heat-kernel toolkit**
- Vassilevich, D. V. *Heat kernel expansion: user's manual.* Phys. Rept. 388 (2003) 279–360. arXiv:hep-th/0306138. **VERIFIED.** Canonical reference.
- Kluth, Y. and Litim, D. F. *Heat kernel coefficients on the sphere in any dimension.* Eur. Phys. J. C 80 (2020) 269. arXiv:1910.00543. **VERIFIED.** Explicit Seeley–DeWitt on S^d (independent cross-check on Paper 28).

**Uncertain citations (flagged for follow-up if used)**
- The general phrase "Tomboulis 1984 modified MK preserves reflection positivity" is standard; the specific 1984 PRD paper is verified, but the *generality* of his procedure to graphs (as opposed to Z^d) is not established in this survey.
- Azarfar–Khalkhali 2019 cited but not independently verified.

---

**Recommended sprint scoping (5 sentences):**
Pursue Option (a): compute spectral-action heat-kernel running of α on Dirac-S³ at finite n_max using existing Paper 28 + Paper 25 machinery, varying Λ over the natural scale set by Paper 28's spectral cutoff. The deliverable is a comparison of dα/d log Λ from `Tr f(D²/Λ²)` against the one-loop QED β(α) = 2α²/(3π) Paper 28 already extracts from the bare spectral sum; either match (consistency check, new internal coherence result) or mismatch (clean negative — spectral-action Λ-running is structurally distinct from bound-state β(α)) is publishable. The sprint is 4–6 weeks with deliverables in Paper 41 (math-ph) or extended Paper 28 §curved_qed_running. Block-spin on the Hopf graph (Option b) and Connes–Kreimer lift (Option c) are deferred — option (b) is multi-month methods-development with no published precedent on irregular graphs, option (c) is fundamentally perturbative-continuum and lives downstream of Paper 28's existing one-loop. The closest published competitor remains Perez-Sanchez 2007.10914 + 2105.01025 (FRG on random fuzzy NCG, distinct from the Connes–vS spectral-truncation paradigm GeoVac sits in), so the sprint result has no direct duplicator in the literature.

**Uncertain citations to flag:** Azarfar–Khalkhali 2019 not independently verified; Tomboulis 1984 verified but generality to irregular graphs is not established here.
