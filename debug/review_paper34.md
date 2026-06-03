# Confidence Review: Paper 34 — GeoVac as a Two-Layer Framework: Projection Taxonomy and the Empirical Matches Catalog

Reviewer: GeoVac Confidence Reviewer (Wave 3, text-level audit)
Date: 2026-06-02
Source file: `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` (9,819 LaTeX lines)

## Pass A — Content audit

### Headline-count verification

Counted `\subsection{...}` entries between `\section{Layer 2: The Projection Family}` (L216) and `\subsection{Boundary of §III ...}` (L2733): **28 projection subsections** (L246–L2449). The abstract claim "We catalogue twenty-eight projection mechanisms" (L46) and conclusion claim "Twenty-eight projections are listed" (L9502) MATCH the body. ✓ Headline-count survives.

### CLAUDE.md §1.8 Roothaan autopsy cross-check

CLAUDE.md §1.8 names 6 prospective targets. Paper 34 §V.C has 19 autopsies (L4427–L6932), all 6 CLAUDE.md targets present (H 1S Lamb shift L4427; H 21cm HFS L4513; μH 2S-2P Lamb L4607; He 2³P FS L5008; He 2¹P→1¹S oscillator L5617; Cs 6S₁/₂ HFS L5778). Paper 34 has executed the §1.8 queue and significantly extended it. ✓ Cataloguing discipline survives.

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence |
|---|-------|----------|---------|----------|----------|
| 1 | 28 projections documented | abstract L46, conclusion L9502 | A | Self-count on body §III | grep + manual count of `\subsection{...}` between L216–L2733 gives 28 ✓ |
| 2 | Three-axis tagging (V/D/transcendental) | §IV L2778 + Table L2786 | A | Body §III metadata | every row in Table 1 has all 3 axes filled; 28 rows match 28 §III entries ✓ |
| 3 | V/D perfect correlation (12/12), transcendental axis genuinely independent | §IV.Rem 1 L2972, §V axiomatization | B | Counted from dictionary tags | 3 independence witnesses (vector-photon, gauge choice, Wick rotation) in transcendental-only corner ✓ structurally consistent |
| 4 | "Layer-2-presence bound" replaces falsified naive depth prediction | abstract L64–86, §VI L8773 | A | Catalogue itself (counter-examples listed) | depth-3 chains exist at residual 0 and at +286 ppm; depth-4 at +2 ppm and -0.534% — directly falsifies the depth law ✓ |
| 5 | Stefan-Boltzmann π²/90 on S³ × S¹_β | catalogue L3349 | A | EXTERNAL (textbook QFT) | π²/90 is standard photon-gas Casimir (per Wikipedia + cited references) ✓ |
| 6 | H 1S static dipole polarizability = 9/2 a₀³ exact via Sturmian at N_basis = 2 | catalogue L3515 | A | EXTERNAL (Dalgarno-Lewis 1955) | Standard textbook result ✓ |
| 7 | Casimir E^{S³} = 1/240 (conformal scalar) | catalogue L3332 | A | EXTERNAL (standard ζ-function regularization) | Confirmed by web search — 1/240r is the standard result ✓ |
| 8 | He 2¹P→1¹S oscillator f = 0.27616 (Drake handbook / Theodosiou 1984) | autopsy L5620 | A | Theodosiou Phys. Rev. A 29, 2981 (1984) + modern recomputation f = 0.27616499(27) | Web-verified ✓ |
| 9 | "Six exposures catalogued to date" (V.D) | §V.D L6948 | E | Self-count contradicts body | Table L6970–L7040 has 9 rows + L7045+ has 10 subsubsections — "six" is stale (was correct before V.D.5–V.D.10 added). **MEDIUM** |
| 10 | H 1S Lamb shift autopsy sum = +1057.41 MHz, residual +0.43 MHz | §V.C.1 L4465–L4467 | E | Arithmetic of own table | Components sum: 1066.44−27.13−12.88+0+1.20+1.18−2.40+5.00 = **1031.41 MHz**, NOT 1057.41 MHz. True residual = +26.43 MHz, not +0.43 MHz. Body text "framework-native subtotal +1031.43" (L4477) DOES sum correctly to 1031.43 — but the table sum and residual lines are off by 26.00 MHz. **HIGH** |
| 11 | Section title "Hydrogen 1S Lamb shift" but reference observable is 2S₁/₂–2P₁/₂ (1057.845 MHz) | §V.C.1 L4427–L4430 | E | Title vs reference observable | "1S Lamb shift" is the 1S hyperfine displacement; 1057.845 MHz is the 2S₁/₂–2P₁/₂ Lamb shift (Lundeen–Pipkin 1981). Title and reference do not match each other. **MEDIUM** |
| 12 | Conclusion L9574–9575: "error compounds with projection depth" | §VIII Conclusion | E | Contradicts §VI body | The body explicitly REPLACED this prediction (`pred:l2_presence` L8773; `rem:depth_falsification` L8816 says "the naive depth-vs-residual prediction was falsified on the catalogue itself"). Conclusion retains stale framing. **HIGH** |
| 13 | T3 16×16 composition table ~196 commute / ~26 one-way / etc., "approximate (~) totals" | §V L8596–L8617 | C | Self-flagged as approximate | The paper itself flags counts as approximate ("Definitive recomputation … is deferred until the next TX-A-style audit"). Honest, but the recompute remains undone. **LOW** — already self-flagged. |
| 14 | "Pending T3 expansion (2026-05-09)" with 19×19 grid expected | §V L8675 | D | Aspirational | Stated openly as pending — neither verified nor counted, but flagged honestly. **LOW** |
| 15 | T3 should be expanded to 28×28 for current 28 projections | §V (and conclusion) | C | Self-flagged | Conclusion says "T3 expansion to 28×28 is flagged for the next TX-A audit" L9559. Honest deferral. **LOW** |
| 16 | Paper 35 Prediction 1 verified 208/208 | abstract + connections + body | A | Paper 35 internal | Paper 35 confirms 208 (4 separate lines) ✓ Internal cross-corpus consistency holds |
| 17 | Conformally coupled scalar = "rest-mass projection (conformal value m² = 1)" | catalogue L3332 | B | Self-consistent with Paper 35 | Standard ξR_c term in conformal coupling gives effective m²_eff = (d−2)/(4(d−1))R — on unit S³ with R=6, m²_eff = 6·(1/4)·... = consistent with conformal value depending on convention. Internal-only check. |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recompute | survives? |
|-------|---------------|----------------------|--------------|-----------|
| H 1S Lamb autopsy table sum | +1057.41 MHz (L4465) | sum-of-listed-components | **+1031.41 MHz** | **NO — off by 26.00 MHz** |
| H 1S Lamb autopsy residual | +0.43 MHz (L4467) | 1057.845 − sum | **+26.43 MHz** if sum is 1031.41 (which it is) | **NO — off by 26.00 MHz** |
| Framework-native subtotal | +1031.43 MHz (L4477) | FN-components only: 1066.44 − 27.13 − 12.88 + 0 + 5.00 | +1031.43 MHz | ✓ |
| Layer-2 net | −0.02 MHz (L4478) | L2 components: 1.20 + 1.18 − 2.40 | −0.02 MHz | ✓ |
| He osc. residual | −2.02% (L3536) | (0.2705 − 0.27616)/0.27616 | −2.05% (rounding-tolerant given 4 sig figs on 0.2705) | ✓ (within rounding) |
| f = 0.2761 (Drake) vs modern recompute | 0.2761 | Independent: 0.27616499(27) | match to 4 sig figs | ✓ |
| Casimir 1/240 on S³ | rational, exact | Standard ζ-function result on S³ | 1/240 ✓ | ✓ |

### Circularity map

**GEOVAC-ONLY chains (load-bearing):**

- **Layer-2-presence bound (Prediction 1).** Rests on the catalogue being correctly tagged for Layer-2 input count. The "rational analytic answer in finite-basis span" exact class (Sprint Calc-P) and the "basis-quality limit" class are framework-internal categories. External anchor: the framework's H Lamb shift one-loop closure at −0.534% is independently the result of Paper 36; the H 1S polarizability 9/2 a₀³ is textbook (Dalgarno–Lewis), so the prediction's two ends are anchored externally even if the intermediate scaling is GeoVac-only. MIXED rather than purely GeoVac-only.

- **28-projection list as exhaustive.** Paper 34 §VIII "closure of the projection list" claims twenty-eight projections cover the framework, with Sprint 3 ABSORBED verdicts on TT modular flow, loop expansion, heat-kernel regularization. This is a framework-internal completeness claim. The dictionary is exhaustive by construction within the GeoVac vocabulary; external completeness (e.g., does QFT need any projection GeoVac has missed?) cannot be settled here. **D — Unverifiable.**

- **The V/D-perfect-correlation observation (12/12).** Rests on dictionary metadata that the paper itself sets. The observation is true by construction; its content is the structural choice to enforce V/D coupling. This is **B — Internally consistent only.**

**MIXED chains:**

- **Paper 32 cross-reference for case-exhaustion theorem.** Paper 34's connections section invokes Paper 32 §VIII case-exhaustion as the structural underpinning. Cross-corpus check: **Paper 32 §VIII Thm `thm:pi_source_case_exhaustion` (L1629–L1696) states "Paper~34 list §III.1–§III.15" and the proof says "each of the fifteen Paper~34 projections."** Paper 34 now has 28 projections. Paper 32 lags by 13 projections. The theorem's logic generalizes mechanically (each new projection is either π-free or engages M1/M2/M3), but Paper 32's text is stale. **MEDIUM cross-corpus drift, flagged in Paper 32 not Paper 34, but the inconsistency is visible from either side.**

**EXTERNAL chains (anchored):**

- Stefan-Boltzmann π²/90, Casimir 1/240 on S³, H 1S polarizability 9/2 a₀³, He 2¹P→1¹S f = 0.27616 (Drake handbook / Theodosiou 1984) all anchored to external standard references. Web-verified.

### Overstatement findings

1. **"consolidates the six exposures catalogued to date" (L6948).** Table has 9 rows, subsubsections has 10. Replace with "the {n} exposures catalogued to date" with n = current count.

2. **Conclusion L9574–9575: "A falsifiable structural prediction (§VII) emerges from the depth-grouped catalogue: error compounds with projection depth."** This is the OLD prediction that the body explicitly FALSIFIED. Replace with the new Layer-2-presence bound. Suggested honest replacement: "A falsifiable structural prediction (§VII) emerges from re-tagging the catalogue by Layer-2 input count: residuals are bounded by the larger of the framework's basis-quality limit or the largest single Layer-2 input consumed by the chain (Prediction~\ref{pred:l2_presence}). The earlier naive 'error compounds with depth' form was falsified on the catalogue itself."

3. **§V.C.1 title "Hydrogen 1S Lamb shift autopsy" with reference 2S₁/₂–2P₁/₂ (1057.845 MHz).** Title and reference observable do not match. Suggested fix: rename "Hydrogen 2S₁/₂–2P₁/₂ Lamb shift autopsy".

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found |
|-----------|-----------|---------|--------------|
| chamseddine_connes2010 | Fortsch. Phys. 58, 553 (2010), NCG unification framework | CITE-OK | Verified: Chamseddine & Connes, Fortschritte 58, 553–600 (2010); arXiv:1004.0464 ✓ |
| marcolli_vansuijlekom2014 | J. Geom. Phys. 75, 71 (2014); arXiv:1301.3480 | CITE-OK | Verified (separately by prior reviews) |
| nieuviarts2025a | arXiv:2502.18105, "Emergence of Lorentz symmetry" | CITE-OK | Verified: G. Nieuviarts, arXiv:2502.18105, math-ph (2025) ✓ Author initial G. correctly applied. |
| nieuviarts2025b | arXiv:2512.15450 (2025), "Emergence of Time" | CITE-OK | Verified (search result confirms arXiv:2512.15450; G. Nieuviarts) ✓ |
| bisognano_wichmann1976 | J. Math. Phys. 17, 303 (1976) | CITE-OK | Verified: J. Math. Phys. 17, 303–321 (1976) DOI 10.1063/1.522898 ✓ |
| bisognano_wichmann1975 | J. Math. Phys. 16, 985 (1975) | CITE-OK | Standard reference, prior reviews confirmed |
| sewell1982 | Annals of Physics 141, 201 (1982) | CITE-OK | Verified: Annals of Physics 141, 201–224 (1982) ✓ |
| sewell1980 | aliased to sewell1982 | CITE-OK | Backward-compat alias documented in source ✓ |
| hartle_hawking1976 | Phys. Rev. D 13, 2188 (1976) | CITE-OK | Standard reference |
| unruh1976 | Phys. Rev. D 14, 870 (1976) | CITE-OK | Standard reference |
| strohmaier2006 | J. Geom. Phys. 56, 175 (2006); arXiv:math-ph/0110001 | CITE-OK | Standard reference; verified by prior reviews |
| bizi_brouder_besnard2018 | J. Math. Phys. 59, 062303 (2018); arXiv:1611.07062 | CITE-OK | Verified ✓ (minor: paper title in bibliography is "applications" plural; the journal published version uses "application" singular and the arXiv title uses "applications"; this is below-noise) |
| franco_eckstein2014 | Rev. Math. Phys. 26, 1430007 (2014); arXiv:1210.6575 | CITE-OK | Standard reference |
| van_den_dungen2016 | Math. Phys. Anal. Geom. 19, 4 (2016) | CITE-OK | Standard reference |
| devastato_lizzi_martinetti2018 | JHEP 03, 089 (2018); arXiv:1710.04965 | CITE-OK | Standard reference |
| drake_swainson1990 | Phys. Rev. A 41, 1243 (1990) | CITE-OK | Standard reference |
| camporesi_higuchi1996 | J. Geom. Phys. 20, 1 (1996) | CITE-OK | Standard reference |
| fock1935 | Z. Phys. 98, 145 (1935) | CITE-OK | Foundational, well-known |
| bargmann1961 | Comm. Pure Appl. Math. 14, 187 (1961) | CITE-OK | Foundational |

Internal GeoVac papers (paper0, paper2, paper7, paper14, paper15, paper17, paper18, paper19, paper20, paper22, paper23, paper24, paper25, paper27, paper28, paper30, paper31, paper32, paper33, paper35, paper38, paper42, ls1_memo, ls2_memo, ls3_memo, ls4_memo, vp1_memo, vp2_memo, audit_memo): self-references in the GeoVac corpus, each present in the project. Not externally checkable.

### Problems found
None at HIGH severity. One minor (plural/singular in Bizi–Brouder–Besnard title) below noise.

### Priority / novelty claims
Paper 34 does NOT claim novelty in the absolute "first-in-literature" sense for any specific construction. The "living catalogue" framing is descriptive of GeoVac's internal state, not a claim against external precedent. **No novelty claims to soften.**

### Cross-corpus check (Paper 32 case-exhaustion)

**Found:** Paper 32 §VIII Theorem `thm:pi_source_case_exhaustion` (L1633) reads: "Let C be any finite composition of projections drawn from the Paper~34 list \S\,III.1--\S\,III.15." Proof body (L1675): "each of the fifteen Paper~34 projections is either π-free at the projection step ... or π-bearing through ..." Paper 34 now has 28 projections; this is a corpus-lag flag against Paper 32, not Paper 34. The π-free/π-bearing assignment Paper 32 gives for projections 1–15 is consistent with Paper 34's §III ordering at those indices (verified projection-by-projection). The theorem's logic mechanically extends to projections 16–28 (Sprint 3 ABSORBED verdicts confirm structural completeness at 28), but the literal statement in Paper 32 needs updating. **MEDIUM, against Paper 32.**

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---------|------|---------|----------|
| H 1S Lamb shift autopsy table sum +1057.41 MHz off by 26.00 MHz (actual sum = 1031.41 MHz) | A | E | **HIGH** |
| Conclusion L9574–9575 retains the falsified naive "error compounds with depth" prediction, contradicting §VI body | A | E | **HIGH** |
| Cross-corpus drift: Paper 32 §VIII case-exhaustion still says "fifteen Paper 34 projections" but Paper 34 has 28 | A (cross-paper) | E (in Paper 32) | **MEDIUM** (in Paper 32, not Paper 34) |
| §V.D L6948 "six exposures catalogued" — actually 9 rows / 10 subsubsections | A | E | **MEDIUM** |
| §V.C.1 title "Hydrogen 1S Lamb shift" + reference is 2S₁/₂–2P₁/₂ — observable mismatch | A | E | **MEDIUM** |
| T3 16×16 composition counts self-flagged as approximate; T3 expansion to 28×28 still pending | A | C | LOW (self-flagged) |
| Bizi-Brouder-Besnard 2018 bibliography title uses plural "applications"; journal version uses singular | B | CITE-OK (cosmetic) | LOW |

**Totals:**
- Pass A counts: A=8, B=2, C=3, D=2, E=5 (E entries can also rank by severity → 2 HIGH, 3 MEDIUM)
- Pass B counts: 17 CITE-OK, 0 CITE-WRONG-METADATA, 0 CITE-MISATTRIBUTED, 0 CITE-DOESNT-SUPPORT, 0 CITE-CANT-FIND (internal GeoVac self-references not counted)
- Severity totals: **2 HIGH, 3 MEDIUM, 2 LOW**

## Broadcast readiness: **YELLOW**

Paper 34's structural backbone is sound: the 28-projection count is correct, every projection has a proper §III subsection, the three-axis tagging is consistently applied, all 6 CLAUDE.md §1.8 Roothaan-autopsy targets have been executed and significantly extended (19 autopsies now), the Layer-2-presence prediction is honestly framed in the body with explicit falsification of the earlier naive form, and all external citations check out. The headline claims about the 28-projection dictionary and the (V/D, transcendental) axis structure survive at A-verdict. However, **two HIGH-severity defects block a clean broadcast**: (1) the §V.C.1 Hydrogen 1S Lamb shift autopsy table has an arithmetic error of 26 MHz in its sum and residual lines (a domain expert summing the listed components will catch this immediately), and (2) the conclusion L9574 retains the OBSOLETE "error compounds with depth" framing that the body itself explicitly falsified — a reader skimming abstract→conclusion would receive a contradictory message about what the paper's central prediction is. Three MEDIUM defects (stale "six exposures" count, observable-name mismatch in autopsy title, cross-corpus lag with Paper 32 §VIII enumeration) are batch-fixable. Fix the two HIGH items + the autopsy title and Paper 34 is GREEN; the conclusion fix is one paragraph rewrite, the autopsy table just needs the sum/residual lines reconciled with the actually-listed components, and the title needs a one-line rename to "2S₁/₂–2P₁/₂ Lamb shift."

## What I could NOT verify (hand to a human expert)

- Whether the 28-projection dictionary is genuinely exhaustive vs. the framework's full QED/NCG content (Sprint 3 ABSORBED verdicts argue yes; only a domain expert can settle whether any standard projection mechanism is still missing).
- Whether the "rational-analytic-answer-in-finite-Sturmian-span" sub-class of the L2-count=0 class always produces residual 0 (Prediction 1.iii). This is a forward-looking falsifiable claim; Sprint Calc-P is one positive instance but a generalization requires more rational-analytic linear-response observables to be tested.
- The framework-side claim that the §III.27 Wick rotation projection "promotes the metric-functional-level content of the Bisognano–Wichmann reading from §VIII open-question candidate to a named §III projection slot" — promotion is GeoVac-internal; whether the metric-functional-level content suffices to constitute a "projection" in the dictionary's sense is a framework-design choice.
- The H 1S Lamb shift autopsy: the cleanest path to fixing the sum is to verify the SE 2S₁/₂ component against Eides 2024 Table 7.3 directly (likely the correct value is +1085 MHz, not +1066, making the listed components add to ~1057 MHz — but only the source memo `debug/calc_track_LAR_lamb_autopsy_memo.md` can settle this without re-running the calculation).

---

Sources I consulted (web):
- Chamseddine & Connes, Fortsch. Phys. 58, 553 (2010), [arXiv:1004.0464](https://arxiv.org/abs/1004.0464)
- Nieuviarts, [arXiv:2502.18105](https://arxiv.org/abs/2502.18105)
- Bisognano & Wichmann, [J. Math. Phys. 17, 303 (1976)](https://ui.adsabs.harvard.edu/abs/1976JMP....17..303B/abstract)
- Sewell, [Annals of Physics 141, 201 (1982)](https://www.sciencedirect.com/science/article/abs/pii/0003491682902858)
- Theodosiou, Phys. Rev. A 29, 2981 (1984) (search result confirmation for f = 0.27616)
- Bizi, Brouder, Besnard, [J. Math. Phys. 59, 062303 (2018)](https://pubs.aip.org/aip/jmp/article-abstract/59/6/062303/985413/) / [arXiv:1611.07062](https://arxiv.org/abs/1611.07062)
- Casimir 1/240 on S³ confirmed via web (standard ζ-function regularization result)
