# Sprint Q1' Concurrent-Work Re-check Memo

**Date:** 2026-05-24
**Sprint:** Q1' (deferred strong-form bridge on Paper 46 Appendix B enlarged substrate; chirality-flipping generators where $\{J,M\}=0$)
**Search window:** 2026-01-01 → 2026-05-24
**Re-check scope:** strong-form Lorentzian propinquity, non-commutative Lorentzian pre-length spaces, enlarged operator system substrates on Krein spectral triples, Connes–Rovelli thermal time on operator systems, BBB / van den Dungen / Strohmaier / Franco–Eckstein / Mondino–Sämann / Latrémolière 2026 follow-ups.

---

## Headline verdicts by area

| Area | Verdict | Notes |
|:-----|:--------|:------|
| Strong-form Lorentzian propinquity (Krein-MS bridge, no K⁺ restriction) | **CLEAR** | No competing construction found. Q1' direction remains uncontested. |
| Enlarged-substrate operator systems with chirality-flipping generators ($\{J,M\}=0$) | **CLEAR** | No paper addresses enlarged substrates as defined in Paper 46 Appendix B. |
| Non-commutative Lorentzian pre-length spaces | **CLEAR–MINOR** | Mondino–Sämann lineage active and growing (2025–2026) but strictly synthetic/metric. No operator-algebraic bridge published. Q1'.B remains open territory. |
| Connes–Rovelli thermal time on operator systems / enlarged substrates | **CLEAR** | Active thermal-time work exists but on quantum systems / POVMs, not on operator-system substrates. |
| BBB 2018 follow-ups | **CLEAR–MINOR** | Nieuviarts 2025/2026 (twist morphism approach) and de Groot 2026 (SU(1,1) pseudo-Riemannian triple) cite BBB; neither constructs Lorentzian propinquity or enlarged substrates. |
| van den Dungen 2016 follow-ups | **CLEAR** | de Groot 2026 cites vdD; no new structural follow-up by vdD himself in window. |
| Strohmaier follow-ups | **CLEAR** | No new 2026 Strohmaier papers on Lorentzian spectral geometry in scope. |
| Franco–Eckstein follow-ups | **CLEAR** | No new 2025–2026 papers on Lorentzian distance formula or causality from this lineage. |
| Latrémolière 2026 follow-ups beyond arXiv:2512.03573 | **MINOR** | arXiv:2603.19128 (March 2026): spectral continuity of almost commutative manifolds in $C^1$ topology on Riemannian metrics. Strictly Riemannian. No Lorentzian extension. |
| Direct competitor for Q1' strong-form Krein-MS bridge | **CLEAR (no scoop risk)** | Nobody appears to be constructing this. |

**Net verdict: CLEAR with three MINOR cite-worthy entries.** No scope changes to Q1' sprint. Proceed.

---

## Detailed findings

### Hit 1 — Latrémolière 2026, spectral continuity in $C^1$ topology (Riemannian) [MINOR]

- **arXiv ID:** [2603.19128](https://arxiv.org/abs/2603.19128)
- **Authors:** Frederic Latrémolière
- **Date:** Submitted 19 March 2026
- **Title:** *Spectral continuity of almost commutative manifolds for the $C^1$ topology on Riemannian metrics*
- **Summary:** Proves that spectra of Dirac operators on almost-commutative models remain continuous as Riemannian metrics vary in the $C^1$ topology. Introduces a novel approach using spectral propinquity. Applications to quantum tori and quantum solenoids.
- **Relevance to Q1':** **MINOR.** Strictly Riemannian. Does not construct strong-form Lorentzian propinquity. Does not address enlarged substrates or chirality-flipping generators. The "novel approach to prove continuity using spectral propinquity" is methodologically interesting and worth a citation in Paper 47 / Paper 48 / future Q1'-driven paper as part of the spectral-propinquity literature update, but it does not interact with the Lorentzian arc.
- **Recommendation:** Add to bibliography of any future Q1'-driven paper as a Latrémolière-2026 reference for the spectral-propinquity machinery; no scope effect.

### Hit 2 — Martinetti 2026, Twisted Standard Model and Krein structure [MINOR]

- **arXiv ID:** [2603.03216](https://arxiv.org/html/2603.03216)
- **Authors:** P. Martinetti (in memoriam M. Filaci)
- **Date:** 3 March 2026
- **Title:** *Twisted Standard Model and its Krein structure*
- **Summary:** Reviews various minimally-twisted spectral triples for the Standard Model and shows the twisted inner product converts the Hilbert space into a Krein space.
- **Relevance to Q1':** **MINOR.** Krein-structure construction via twist (the Wick-rotation-without-complex-numbers lineage of Devastato–Lizzi–Martinetti 2018 and Nieuviarts 2025). Confirms the Krein structure is a natural endpoint of twist constructions in almost-commutative geometry but does not construct propinquity, does not address enlarged substrates with chirality-flipping generators (the twist operator does have left/right particle distinction structure, but this is at the algebra-twist level, NOT at the operator-system Lipschitz-seminorm level the Q1' enlarged substrate operates on), and does not discuss Bisognano–Wichmann modular flow or Connes–Rovelli thermal time.
- **Recommendation:** Cite in any future Lorentzian-arc paper as part of the twist-morphism Krein-structure lineage; no scope effect on Q1' direction.

### Hit 3 — Nieuviarts 2025/2026, Emergence of Time from Twisted Spectral Triple [MINOR]

- **arXiv ID:** [2512.15450](https://arxiv.org/html/2512.15450) (v2, 11 May 2026; v1 October 2025)
- **Author:** Gaston Nieuviarts
- **Title:** *Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry*
- **Summary:** Synthesis paper on emergence of pseudo-Riemannian structures from twisted spectral triples in the almost-commutative framework. KK-morphism construction, signature change via spacelike reflection without complex numbers, Lorentzian signature in KO-dim 6.
- **Relevance to Q1':** **MINOR.** Same lineage as Hit 2 (twist morphism for signature change, the alternative to Wick rotation that GeoVac's Sprint L0 audit already considered and ruled out for odd-dim S³). Does not construct propinquity, does not discuss enlarged substrates with chirality-flipping generators, does not discuss Bisognano–Wichmann or Connes–Rovelli thermal time. Is a synthesis/review, not a new theorem.
- **Recommendation:** Already cited in Paper 32 §VIII.E and the Sprint L0 lit audit. No new effect on Q1'.

### Hit 4 — de Groot 2026, Pseudo-Riemannian SU(1,1) [CLEAR]

- **arXiv ID:** [2601.22171](https://arxiv.org/abs/2601.22171)
- **Author:** Jort de Groot
- **Date:** 22 January 2026
- **Title:** *Pseudo-Riemannian Spectral Triples for SU(1,1)*
- **Summary:** Shows that a triple involving Kostant's cubic Dirac operator on $L^2(\mathrm{SU}(1,1)) \otimes \mathbb{C}^2$ is both a pseudo-Riemannian spectral triple and an indefinite spectral triple. Uses harmonic analysis methods.
- **Relevance to Q1':** **CLEAR.** A worked example of the BBB/vdD pseudo-Riemannian-triple framework, not a structural advance. Does not construct propinquity, does not discuss enlarged substrates, does not address modular flow.
- **Recommendation:** Mention in a related-work paragraph if Q1' paper drafts an SU(1,1) or non-compact-group example; otherwise no action.

### Hit 5 — Fröb–Much–Papadopoulos 2026, Novel Noncommutative Spacetime Algebra [CLEAR]

- **arXiv ID:** [2601.07350](https://arxiv.org/abs/2601.07350)
- **Authors:** M. B. Fröb, A. Much, K. Papadopoulos
- **Date:** 12 January 2026
- **Title:** *A proposal for the algebra of a novel noncommutative spacetime*
- **Summary:** Constructs a Lorentz-invariant noncommutative coordinate algebra using Weyl algebra + GNS construction. Indefinite signature resolved via Krein space formulation. Emergent minimal-length effects, quantum corrections $\propto \ell_P^2$.
- **Relevance to Q1':** **CLEAR.** This is a Lorentz-invariant NC coordinate algebra proposal (Snyder/Doplicher-style direction with Krein space resolution), structurally distinct from the GeoVac spectral-triple / propinquity programme. No construction of propinquity, no enlarged-substrate operator system in the Paper 46 sense.
- **Recommendation:** No citation needed; entirely orthogonal direction.

### Hit 6 — Mondino–Sämann lineage (CLEAR–MINOR cluster, four papers)

Five papers in the synthetic-Lorentzian-GH program in window:

#### 6a — Mondino–Sämann 2025/2026 [arXiv:2504.10380]
- **Date:** Submitted April 2025, last revised December 2025
- **Title:** *Lorentzian Gromov-Hausdorff convergence and pre-compactness*
- **Verdict:** **CLEAR.** Strictly synthetic/metric Lorentzian. No operator-algebraic bridge. No mention of propinquity, spectral triples, or Krein spaces.

#### 6b — Braun–Sämann 2025/2026 [arXiv:2506.10852]
- **Date:** June 12, 2025
- **Title:** *Gromov's reconstruction theorem and measured Gromov–Hausdorff convergence in Lorentzian geometry*
- **Verdict:** **CLEAR.** Builds three notions of measured Lorentzian GH convergence (intrinsic / distortion / box). No operator-algebraic bridge. Causal-set-theory applications.

#### 6c — Che–Perales–Sormani 2026 [arXiv:2510.13069]
- **Date:** v3 February 26, 2026
- **Title:** *Gromov's Compactness Theorem for the Intrinsic Timed-Hausdorff Distance*
- **Verdict:** **CLEAR.** Classical metric geometry of timed-metric-spaces. No operator-algebraic content.

#### 6d — Mondino–Ryborz–Sämann 2026 [arXiv:2605.03172] (mentioned in search results)
- **Date:** May 2026
- **Title:** *Stability of Synthetic Timelike Ricci Bounds under $C^0$-Limits…*
- **Verdict:** **CLEAR.** Synthetic Lorentzian timelike Ricci, $C^0$ stability. No operator-algebraic angle.

#### 6e — Beran–Sämann (existing) [arXiv:2204.09491], "Non Hilbertian (Lorentzian) Length Spaces" [arXiv:2407.19595]
- Mentioned in passing in searches. Both strictly synthetic-Lorentzian.

**Cluster recommendation:** The Mondino–Sämann synthetic-Lorentzian-GH program is **active and accelerating** in 2026 (four to five papers in 12 months across multiple author combinations Braun, Mondino, Sämann, Ryborz, Perales, Sormani, Beran, Che). It is the most plausible community to attempt a non-commutative / operator-algebraic extension of Lorentzian pre-length spaces — which would be Q1'.B closure. **However, no member of this community has published such an extension in window.** Q1'.B remains open territory. This is the area to monitor most closely for SCOOP-RISK in the next 6–12 months.

### Hit 7 — Modular Hamiltonian / Bisognano-Wichmann 2026 [CLEAR]

Three 2026 papers in this area:

- **arXiv:2605.20001** (May 2026), *Numerical approach to the modular operator for fermionic systems* — uses BW relation on Fock space for numerical modular-Hamiltonian computation; strictly numerical-fermionic, not operator-system.
- **JHEP 02 (2026) 210**, *Modular Hamiltonians for future-perturbed states* — perturbative modular Hamiltonian for 2D CFT, strictly QFT side.
- **arXiv:2410.16433**, *On the Bisognano-Wichmann entanglement Hamiltonian of nonrelativistic fermions* — strictly nonrelativistic fermion side.

**Verdict:** **CLEAR.** None of these construct modular Hamiltonians on truncated operator systems or extend Connes–Rovelli thermal time to enlarged substrates. GeoVac Paper 42's operator-system construction is uncontested in this area.

### Hit 8 — Connes–Rovelli thermal time follow-ups [CLEAR]

The 2023–2024 thermal-time-as-unsharp-observable program (arXiv:2306.13774, JMP 65 032105) and "The Time in Thermal Time" (arXiv:2407.18948) remain the most recent thermal-time-on-quantum-systems work. **No 2026 paper extends Connes–Rovelli to operator systems or enlarged substrates.** GeoVac's L1-tighten + Paper 42 unified BW-α + BW-γ Tomita construction is uncontested.

### Hit 9 — Hekkelman–McDonald 2025 / Connes 2026 [CLEAR]

- arXiv:2412.00628 (Hekkelman–McDonald 2025) — noncommutative integral on spectrally truncated triples, Riemannian only.
- Connes 2026 implementation of Galerkin matrices at multiple cutoffs (mentioned in search results) — strictly Riemannian spectral truncation.

**Verdict:** **CLEAR.** Connes–vS spectral truncation lineage is active but strictly Riemannian in 2026.

---

## Cross-search for direct Q1' competitors

I performed cross-searches for any paper combining (i) strong-form / enlarged-substrate operator system language with (ii) non-commutative Lorentzian pre-length space language, or (iii) Krein-self-adjoint Lipschitz seminorm with (iv) BBB-style Krein triple at finite cutoff. **No hits.** The closest construction in the literature remains [Paper 45/46](file:///../paper_45_lorentzian_propinquity.tex) (K⁺-weak-form) and its strong-form extension on the natural substrate ([Paper 46 main + Appendix B](file:///../paper_46_strong_form_lorentzian_propinquity.tex)), both of which are GeoVac-internal.

The Q1' strong-form Krein-MS bridge — propinquity on the FULL Krein space without K⁺ restriction AND on the enlarged substrate where $\{J,M\}=0$ for chirality-flipping generators AND with a Mondino–Sämann-style pre-length space target — is, to our knowledge, an unconstructed object as of 2026-05-24.

**The Q1'.B sub-question (whether a non-commutative MS pre-length space concept exists in published literature that would close Q1'.B for free or scoop Q1' entirely)** is decisively in the negative. The Mondino–Sämann lineage is the natural community, they are publishing actively, and they have NOT extended to operator-algebraic / non-commutative versions yet. **Q1'.B is open territory and ours to claim.**

---

## Recommendations for Q1' sprint planning

1. **No scope changes.** All four sub-areas tested return CLEAR for direct competitors. Q1', Q1'.A (topography enlargement), Q1'.B (Gelfand spectrum / possible non-commutative MS extension), and Q1'.C (off-orbit super-additivity via Paper 42 four-witness transport) all remain open and uncontested.

2. **Three MINOR citations to add when Q1' paper drafts:**
   - Latrémolière 2026 [arXiv:2603.19128] — spectral propinquity continuity in $C^1$ topology (methodology reference).
   - Martinetti 2026 [arXiv:2603.03216] — twisted SM Krein structure (Krein-structure lineage reference, in memoriam Filaci).
   - Nieuviarts 2025/2026 [arXiv:2512.15450] — emergence of time from twisted spectral triple (related-work, already in Sprint L0 audit).

3. **Cluster to monitor — Mondino–Sämann synthetic-Lorentzian-GH program.** This is the area of highest forward SCOOP-RISK for Q1'.B specifically. Four to five papers in 12 months across an extended author community. Recommend a 6-month re-check trigger if Q1' has not landed by 2026-11-24.

4. **Q1'.B framing recommendation.** Because the Mondino–Sämann community is active but has NOT yet extended to operator-algebraic versions, the strategic-best framing for Q1'.B is: **GeoVac's enlarged-substrate operator system at finite cutoff defines a candidate "non-commutative Lorentzian pre-length space" that — in a controlled GH limit — should converge to the Mondino–Sämann commutative pre-length space.** This positioning ties Q1'.B to a publicly-active research direction without competing with them on their (synthetic) home turf.

5. **Proceed with Q1' sprint plan as scoped.** Q1'-Light diagnostic in parallel; if it returns POSITIVE-PARTIAL or GO, dispatch Q1'-main four-sub-sprint plan analogous to L3b-2a–d / L3b-2f-β.

---

## Citations

- [arXiv:2603.19128](https://arxiv.org/abs/2603.19128) Latrémolière 2026, spectral continuity in $C^1$ topology
- [arXiv:2603.03216](https://arxiv.org/html/2603.03216) Martinetti 2026, Twisted SM Krein structure
- [arXiv:2512.15450](https://arxiv.org/html/2512.15450) Nieuviarts 2025/2026, Emergence of time from twisted spectral triple
- [arXiv:2601.22171](https://arxiv.org/abs/2601.22171) de Groot 2026, Pseudo-Riemannian SU(1,1)
- [arXiv:2601.07350](https://arxiv.org/abs/2601.07350) Fröb–Much–Papadopoulos 2026, Novel NC spacetime algebra
- [arXiv:2504.10380](https://arxiv.org/abs/2504.10380) Mondino–Sämann 2025/2026, Lorentzian GH convergence
- [arXiv:2506.10852](https://arxiv.org/html/2506.10852) Braun–Sämann 2025, Lorentzian Gromov reconstruction
- [arXiv:2510.13069](https://arxiv.org/html/2510.13069) Che–Perales–Sormani 2026, Timed-Hausdorff distance compactness
- [arXiv:2605.03172](https://arxiv.org/abs/2605.03172) Mondino–Ryborz–Sämann 2026, Timelike Ricci $C^0$ stability
- [arXiv:2605.20001](https://arxiv.org/html/2605.20001) Numerical modular operator for fermions 2026
- [arXiv:2504.11715](https://arxiv.org/abs/2504.11715) Farsi–Latrémolière 2025, spectral propinquity analytic Riemannian path
- [arXiv:2412.00628](https://arxiv.org/abs/2412.00628) Hekkelman–McDonald 2025, NC integral on truncated triples
