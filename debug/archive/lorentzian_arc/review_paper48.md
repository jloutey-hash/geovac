# Confidence Review: Paper 48 — Krein-pointed proper QMS and the bridge to Mondino-Sämann Lorentzian pre-length spaces, with application to the GeoVac hemispheric wedge

**Reviewer**: GeoVac Confidence Reviewer (combined auditor + citation checker), Wave 2 Lorentzian-arc audit.
**Source audited**: `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` (text-only, .tex source).
**Date**: 2026-06-02.

## Calibration check

Not a calibration run. Audit performed against the stated standards in `agents/CONFIDENCE_REVIEW.md`.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence |
|:--|:------|:---------|:--------|:---------|:---------|
| 1 | Wick-rotation functor $\Wfun: \KreinMM \to \LorPLG$ is well-defined and contravariant by Gelfand duality | abstract, Def 4.5 | B | GEOVAC-ONLY (Papers 42, 44, 45, 47) plus external (Latrémolière 2025, Mondino-Sämann v4, Connes-Rovelli 1994) | Both load-bearing externals exist (verified arXiv 2512.03573 and 2504.10380v4); functor construction is novel and internally consistent. No external check of categorical correctness. |
| 2 | Bridge Theorem 6.4 properties B1/B2/B3/B4 hold at theorem-grade rigor | abstract, Thm 4.7 | B | GEOVAC-ONLY at proof level (Papers 42/44/45/47 are GeoVac-internal); chains into external Latrémolière/Mondino-Sämann results | All proofs in Paper 48 are "Proof sketch" — none reach formal-proof level in this paper. Domain-expert verification needed for full theorem-grade claim. |
| 3 | Numerical panel $\Lprop \in \{2.0746, 1.6101, 1.3223\}$ inherited bit-identically from Paper 45 | Thm 5.6 (T3), §6.2 | B | GEOVAC-ONLY (Paper 45 Sub-Sprint D) | Cannot recompute from this paper alone; inheritance assertion not verifiable here. |
| 4 | Off-orbit super-additivity Case (iii) is structurally EMPTY at K⁺-weak-form via M-diagonal topography (Lemma 4.10) | §4.4 | A | EXTERNAL (basic set-theoretic exhaustion: $(N_x,L_x)=(N_y,L_y)\neq(N_z,L_z)=(N_x,L_x)$ is a contradiction) | The case-exhaustion proof in §4.4 (Lemma 4.10) is correct symbolic logic. |
| 5 | "First quantitative pLGH-convergence panel in the Mondino-Sämann literature from operator-algebraic input" | abstract, Thm 5.6, §6.2, §6.5 | C | The novelty claim itself: cannot confirm absence | Per protocol §"honest ceiling," a clean search downgrades "first in literature" → "to our knowledge, the first." Search returned no prior art (mondino_samann2025 §10 lists only continuous spacetime examples — Chruściel-Grant, warped products). Hence: claim survives in softened form. |
| 6 | T6 G2 metric-level closure for GeoVac wedge — workaround for published-open Latrémolière non-compact propinquity | abstract, Thm 5.9 | B | GEOVAC-ONLY (rests on Paper 47 + bridge functor + Theorems B2/B4) | T6 is conditional on Bridge Theorem (claim 2). The "wedge-specific workaround" framing is honest (paper explicitly distinguishes wedge from general carrier in §5.8). |
| 7 | F2 mismatch resolved via R3 Connes-Rovelli thermal time | abstract, §4.1.3 | B | MIXED — Connes-Rovelli 1994 EXTERNAL verified; specific identification $t_{thermal} = -i t_{geometric}/\kappa_g$ on the BW wedge is GeoVac-internal interpretation | Connes-Rovelli citation is real and on-topic. The functorial-not-isometric framing is honestly stated. |
| 8 | Six BW witnesses reduce to two synthetic-side equivalence classes by $\kappa_g$ rescaling (T4) | §5.5 (Thm 5.7) | B | GEOVAC-ONLY (Paper 42 six-witness collapse) | Honest derivation if Paper 42 holds. |
| 9 | Pythagorean orthogonality with $1/\pi^2$ M1 signature transports to synthetic side (T5) | §5.6 (Thm 5.8) | B | GEOVAC-ONLY (Paper 43 §10.2) | Internal transport — assumes Paper 43 closure. |
| 10 | Riemannian-limit recovery at $N_t = 1$ bit-exact | §4.6 (Cor 4.13) | B | GEOVAC-ONLY (Paper 45 Sub-Sprint D Lem 5.1) | Inherited claim. |
| 11 | Honest scope: K⁺-weak-form only; Q1' strong-form deferred | abstract, §1.4, §7.1 | A | EXTERNAL methodological standard | Scope properly fenced — verbatim "the bridge is K⁺-weak-form only" and explicit Q1' three-step decomposition. |

### Numbers I recomputed

No numerical claims in Paper 48 are independently verifiable from this paper alone — every number ($2.0746$, $1.6101$, $1.3223$, $\kappa_g = 1/4, 1/8, 1, 2$, dim $W_L = 16, 60, 160$, propagation $= 2$, slab $2\pi$) is asserted as inherited from Paper 42/43/44/45/47. The recompute task would require those papers' code/proofs in front of me.

**Independent reference for the Mondino-Sämann pLGH framework**: arXiv:2504.10380v4 verified via WebFetch — Def 3.6 (LGH convergence), Def 3.8 (covered LPLS), Def 3.12 (pLGH) exist as cited. **The structural ground that Paper 48 builds on is real.**

**Independent reference for Latrémolière hypertopology**: arXiv:2512.03573 verified — pointed proper QMS framework with hypertopology and inframetric is real. (4-point relaxed forward inequality not explicitly confirmed via abstract; published in December 2025.)

### Circularity map

**GEOVAC-ONLY chains** (most at risk):

- **Claim 2 (Bridge Theorem)** rests on Papers 42, 44, 45, 47 for every load-bearing input; if any of those has an upstream error, B1–B4 cascade.
- **Claim 3 (numerical panel)** rests entirely on Paper 45 Sub-Sprint D; this paper makes no independent computation.
- **Claim 6 (T6 G2 closure)** rests on the conjunction of Bridge Theorem (claim 2) + Paper 47 three-carrier identification — two GeoVac chains.
- **Claim 8 (T4 synthetic four-witness)** rests on Paper 42 six-witness collapse.
- **Claim 9 (T5 Pythagorean)** rests on Paper 43 §10.2.

**MIXED chains** (good external core + GeoVac interpretive layer):

- **Claim 1 (functor existence)** — Latrémolière + Mondino-Sämann frameworks are EXTERNAL and verified; the specific Wick-rotation functor construction is GeoVac-internal.
- **Claim 7 (F2 via Connes-Rovelli)** — Connes-Rovelli 1994 thermal-time hypothesis is established and on-topic; the wedge-specific identification is GeoVac-internal.

**EXTERNAL** (solid):

- **Claim 4 (Case (iii) empty)** — direct set-theoretic exhaustion verifiable by inspection.
- **Claim 11 (scope fencing)** — methodological standard, applied properly.

### Overstatement findings

- **Abstract**: "To the author's knowledge this is the first quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input." → **Already softened** ("to the author's knowledge"). This is the correct epistemic level per protocol. No change needed. (For comparison, §5.6 and §6.5 use stronger "the first quantitative pLGH-convergence panel" — these should be softened to match the abstract's "to the author's knowledge" framing.)

  **Suggested replacement**:
  - §5.6 last sentence "This is the **first quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input.**" → "This is, to our knowledge, the **first** quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input."
  - §6.2 first sentence "the GeoVac wedge family furnishes the first quantitative pLGH-convergence panel" → "...furnishes, to our knowledge, the first quantitative pLGH-convergence panel."

- **Abstract**: "Bridge Theorem 6.4 holds at theorem-grade rigor" → Paper 48 contains only "Proof sketch" environments for all four properties (B1, B2, B3, B4). The claim "theorem-grade rigor" is therefore overstated when read against the proofs actually written into this paper, even if the underlying proofs exist in the cited Phase-A memos.

  **Suggested replacement**: "Bridge Theorem 6.4 holds at theorem-grade rigor (full proofs in Phase-A sub-sprint memos; this paper presents proof sketches with detailed references)." Or alternatively, expand the proof sketches into full proofs in the paper itself before broadcast.

- **§1.4 / §5.8**: "the substantive deepest content" / "the deepest substantive new content" — the framing is fine but slightly self-congratulatory; consider softening or removing the superlative, since it is hard for the reader to evaluate.

- **§1.5 last paragraph**: "the bridge constructed here is novel" — same priority issue as the abstract. Already implicitly softened by "to the author's knowledge" earlier in the paragraph, so this is acceptable.

## Pass B — Citation and novelty

### Citation table

| key | claimed as | verdict | what I found |
|:---|:----------|:-------|:------------|
| `allen_burtscher2022` | "Properties of the null distance and spacetime convergence," IMRN 2022 | CITE-OK (assumed; not independently verified — low-risk) | n/a |
| `bertozzini_conti_lewkeeratiyutkul2009` | arXiv:0901.4031 | **CITE-MISATTRIBUTED** | arXiv:0901.4031 is "Convergence Radii for Eigenvalues of Tri-diagonal Matrices" by Adduci-Djakov-Mityagin (Jan 2009); has nothing to do with Bertozzini et al. The actual Bertozzini-Conti-Lewkeeratiyutkul "Modular theory, non-commutative geometry and quantum gravity" is in SIGMA 2010 (arXiv:1007.4094). **The arXiv ID is wrong.** |
| `bisognano_wichmann1976` | J. Math. Phys. 17 (1976), 303-321 | CITE-OK (assumed; standard reference) | n/a |
| `bizi_brouder_besnard2018` | "Spectral action in Lorentzian signature," Class. Quantum Grav. 35 (2018) 175004. arXiv:1611.07062 | **CITE-WRONG-METADATA / DOESNT-SUPPORT** | arXiv:1611.07062 has title "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics" by Bizi-Brouder-Besnard, published in *J. Math. Phys.* 59 (2018) 062303, NOT *Class. Quantum Grav.* 35 175004. **Both title and journal/article number are wrong.** (A separate Bizi-Brouder-Besnard paper "Lorentz signature and twisted spectral triples" is in JHEP 2018 — also not the cited venue.) |
| `camporesi_higuchi1996` | J. Geom. Phys. 20 (1996), 1-18 | CITE-OK (assumed standard, used across the math.OA series) | n/a |
| `chamseddine_connes1997` | Comm. Math. Phys. 186 (1997), 731-750 | CITE-OK (assumed standard) | n/a |
| `che_perales_sormani2025` | DGA 103 (2026), arXiv:2510.13069 | CITE-OK | Verified: "Gromov's Compactness Theorem for the Intrinsic Timed-Hausdorff Distance" by Che-Perales-Sormani, 2025 submission. |
| `connes1994` | Noncommutative Geometry, Academic Press, 1994 | CITE-OK (book reference) | n/a |
| `connes_rovelli1994` | Class. Quantum Grav. 11 (1994), 2899-2917 | CITE-OK (standard, load-bearing) | Verified |
| `connes_vs2021` | CMP 383 (2021), arXiv:2004.14115 | CITE-OK | Verified: title and authors match. |
| `farsi_latremoliere2024` | "Quantum GH propinquity for crossed products...," J. Funct. Anal. 286 (2024), 110293 | CITE-OK (assumed; consistent with related citation 2025) | n/a |
| `farsi_latremoliere2025` | arXiv:2504.11715 | CITE-OK | Verified |
| `franco_eckstein2014` | "Causality in noncommutative geometry," Class. Quantum Grav. 30 (2013), 135007 | CITE-WRONG-METADATA (minor: year says "(2014)" in key but body says (2013); this is internal inconsistency only) | Plausibly OK on author/topic; not independently verified. |
| `hartle_hawking1976` | Phys. Rev. D 13 (1976), 2188-2203 | CITE-OK (standard reference) | n/a |
| `hekkelman_mcdonald2024` | arXiv:2412.00628, "The semi-classical limit of Connes spectral triples on a torus" | **CITE-WRONG-METADATA** | arXiv:2412.00628 actual title is "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity" by Hekkelman-McDonald. **Title is wrong.** Authors and arXiv ID match the real paper, but the title attribution is fabricated/wrong. |
| `hekkelman_mcdonald2024b` | arXiv:2411.04566, "Noncommutative integration on non-compact spectral triples" | **CITE-MISATTRIBUTED (Fursaev-Solodukhin failure mode)** | arXiv:2411.04566 is "Exponential improvements to the average-case hardness of BosonSampling" by Bouland-Datta-Fefferman-Hernandez (quantum complexity theory). It has NOTHING to do with Hekkelman-McDonald or spectral triples. **The arXiv ID is fabricated/wrong.** I could not find a published Hekkelman-McDonald paper with title "Noncommutative integration on non-compact spectral triples"; the only Hekkelman-McDonald 2024 paper appears to be arXiv:2412.00628 (the same one mis-cited at `hekkelman_mcdonald2024`). |
| `ketterer2026` | arXiv:2605.11271, "$\ell$-convergence on metric measure spaces and applications" | **CITE-WRONG-METADATA** | arXiv:2605.11271 actual title is "Convergence of Lorentzian spaces and curvature bounds for generalized cones" by Christian Ketterer (May 11, 2026). Author and submission date match; **title is wrong**. |
| `kubota2026` | arXiv:2605.09101, "A Lorentzian coarea inequality" | CITE-OK (essentially correct; actual is "Lorentzian coarea inequality" without "A") | Verified |
| `kunzinger_samann2018` | Ann. Global Anal. Geom. 54 (2018), 399-447; arXiv:1711.08990 | CITE-OK | Verified |
| `kostant1965` | Bull. AMS 75 (1965), 627-642 | Date inconsistency: Kostant 1965 paper in vol 75 was published in 1969, not 1965. The volume number is wrong if year is 1965; vol 71 is correct for 1965. **Minor metadata issue** — not load-bearing in this paper (Kostant is not cited in the body). | Possibly stale or vestigial bibitem. |
| `latremoliere_metric_propinquity_2015` | J. Math. Pures Appl. 103 (2015), 303-351 | CITE-OK (standard) | n/a |
| `latremoliere2018` | Trans. AMS 368 (2016), 365-411; arXiv:1302.4058 | CITE-OK | Verified |
| `latremoliere_pp_2018af` | J. Math. Anal. Appl. 469 (2019), 378-404 | CITE-OK (assumed) | n/a |
| `latremoliere_metric_st_2017` | "The dual GH propinquity for metric spectral triples," Adv. Math. 411 (2023), 108790; arXiv:1811.10843 | **CITE-WRONG-METADATA (minor)** | The arXiv 1811.10843 title is "The Gromov-Hausdorff propinquity for metric Spectral Triples" — without "dual". Paper 48 adds "dual" to the title. Author, journal, year align with related work. (Latrémolière's "dual" variant is a separate 2015 paper.) Minor title attribution error. |
| `latremoliere2025_hypertopology` | arXiv:2512.03573, "Quantum GH hypertopology on pointed proper QMS" | CITE-OK | Verified, load-bearing, central to Paper 48. |
| `leimbach_vs2024` | Leimbach-van Suijlekom, "Convergence of Peter-Weyl truncations of compact quantum groups," Adv. Math. 439 (2024), 109496. arXiv:2409.16698 | **CITE-WRONG-METADATA** | arXiv:2409.16698 has single author Malte Leimbach (no van Suijlekom co-author), title is correct, submitted Sept 2024 / revised April 2025, to appear in J. Noncommut. Geom. (NOT Adv. Math.). **Authorship and journal are wrong.** (Note: CLAUDE.md repeatedly cites this work as "Leimbach-van Suijlekom 2024 (Adv. Math. 439, 109496)" — the same erroneous citation appears throughout the GeoVac corpus.) |
| `marcolli_vs2014` | "Gauge networks in NCG," J. Geom. Phys. 75 (2014), 71-91; arXiv:1301.3480 | CITE-OK | Verified |
| `minguzzi_suhr2024` | "Bounded Lorentzian metric spaces," Lett. Math. Phys. 114 (2024) 73, arXiv:2209.14384 | **CITE-WRONG-METADATA** | arXiv:2209.14384 actual title is "Lorentzian metric spaces and their Gromov-Hausdorff convergence" by Minguzzi-Suhr, Lett. Math. Phys. 114 (2024). Author, journal, year align; **title is wrong.** |
| `mondino_ryborz_samann2025` | "Stability of synthetic timelike curvature bounds," arXiv:2605.03172 (May 2026) | CITE-WRONG-METADATA (minor) | Actual title is "Stability of Synthetic Timelike Ricci Bounds under $C^0$-Limits and Applications to Impulsive Gravitational Waves" — Paper 48's title is a shorter paraphrase. Authors, arXiv ID, date all match. |
| `mondino_samann2025_pointed` | "Pointed Lorentzian GH convergence of covered pre-length spaces," arXiv:2504.10380v4 (Dec 2025) | CITE-WRONG-METADATA (title) but CITE-OK on content | The arXiv abstract title is "Lorentzian Gromov-Hausdorff convergence and pre-compactness" — the paper's actual title in v4 does NOT say "pointed" or "covered" in the title. **However, content verified**: Defs 3.6, 3.8, 3.12 (LGH convergence, covered LPLS, pLGH) exist in v4 as cited. The substantive ground IS real; only the title attribution is off. **This is the most load-bearing citation — flagging because precision matters here.** |
| `muller2022` | Müller, "Lorentzian GH convergence and Cauchy slabs," CMP 391 (2022), 855-882 | CITE-OK (assumed standard) | n/a |
| `nieuviarts2025_v2` | "Twisted spectral triples and pseudo-Riemannian spectral geometry, v2," arXiv:2512.15450 | **CITE-WRONG-METADATA** | arXiv:2512.15450 actual title is "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry." Authors and v2 date (May 11, 2026) match; **title is wrong.** |
| `reed_simon_iv` | Reed-Simon vol IV | CITE-OK (standard book) | n/a |
| `sakovich_sormani2024` | "Timed GH convergence of spacetimes with timelike curvature bounds," arXiv:2410.16800 | **CITE-WRONG-METADATA** | arXiv:2410.16800 actual title is "Introducing Various Notions of Distances between Space-Times" by Sakovich-Sormani, 2024 (latest v3 Oct 2025). Authors, arXiv ID match; **title is wrong.** |
| `sewell1982` | Ann. Phys. (NY) 141 (1982), 201-224 | CITE-OK (standard) | n/a |
| `sormani_vega2016` | Class. Quantum Grav. 33 (2016), 085001 | CITE-OK (assumed standard) | n/a |
| `strohmaier2006` | J. Geom. Phys. 56 (2006), 175-195 | CITE-OK (assumed) | n/a |
| `unruh1976` | Phys. Rev. D 14 (1976), 870-892 | CITE-OK (standard) | n/a |
| `vandungen2016` | Math. Phys. Anal. Geom. 19 (2016), 4; arXiv:1505.01939 | CITE-OK | Verified |
| `paper24`-`paper47` | GeoVac internal preprints | n/a (internal — circularity flagged in Pass A) | n/a |

### Problems found

**CITE-MISATTRIBUTED (HIGH severity — load-bearing or Fursaev-Solodukhin failure mode):**

1. **`hekkelman_mcdonald2024b` → arXiv:2411.04566**: this arXiv ID is "Exponential improvements to the average-case hardness of BosonSampling" by Bouland et al., NOT a Hekkelman-McDonald spectral-triple paper. The Hekkelman-McDonald paper with the claimed title "Noncommutative integration on non-compact spectral triples" appears not to exist; the only 2024 Hekkelman-McDonald paper is arXiv:2412.00628 (already cited separately as `hekkelman_mcdonald2024`). **Fix: either correct the arXiv ID or remove the bibitem** (Paper 48 cites this only in §1.5 as part of the adjacent-work paragraph; it is not load-bearing for any theorem, but the misattribution would embarrass on broadcast).

2. **`bertozzini_conti_lewkeeratiyutkul2009` → arXiv:0901.4031**: this arXiv ID is "Convergence Radii for Eigenvalues of Tri-diagonal Matrices" by Adduci-Djakov-Mityagin, NOT a Bertozzini et al. paper. The actual Bertozzini-Conti-Lewkeeratiyutkul work appears in arXiv:1007.4094 (SIGMA 2010). **Fix: correct the arXiv ID to 1007.4094 and update venue.**

**CITE-WRONG-METADATA (MEDIUM severity — fixable but cosmetic if low-load-bearing):**

3. `bizi_brouder_besnard2018` — title "Spectral action in Lorentzian signature" and journal "Class. Quantum Grav. 35 (2018) 175004" both wrong. Actual: "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics," *J. Math. Phys.* 59 (2018) 062303. Used in §2.4 (Lorentzian Dirac via vdD 2016) — cited in support of Krein construction. Load-bearing: not directly; the cite serves as adjacency, not as proof input. Still HIGH-MEDIUM because wrong-venue/wrong-title looks careless on broadcast.

4. `hekkelman_mcdonald2024` (arXiv:2412.00628) — title wrong (real: "A noncommutative integral on spectrally truncated spectral triples..."). MEDIUM.

5. `ketterer2026` — title wrong (real: "Convergence of Lorentzian spaces and curvature bounds for generalized cones"). MEDIUM.

6. `leimbach_vs2024` — authorship wrong (single author Leimbach, not Leimbach-van Suijlekom) and journal wrong (J. Noncommut. Geom., not Adv. Math.). **NOTE: this same erroneous attribution appears throughout CLAUDE.md and the GeoVac corpus** — this is corpus-wide, not just Paper 48. MEDIUM-HIGH on its own; HIGH if PI cares about corpus-wide propagation.

7. `minguzzi_suhr2024` — title wrong (real: "Lorentzian metric spaces and their Gromov-Hausdorff convergence"). Author/journal/year correct. MEDIUM.

8. `mondino_samann2025_pointed` — title wrong (real: "Lorentzian Gromov-Hausdorff convergence and pre-compactness"). The substantive content (covered LPLS, pLGH) IS in the paper (Defs 3.6, 3.8, 3.12 confirmed in v4 HTML). MEDIUM — but this is the most load-bearing cite of Paper 48, so the title should be correct.

9. `nieuviarts2025_v2` — title wrong (real: "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry"). MEDIUM.

10. `sakovich_sormani2024` — title wrong (real: "Introducing Various Notions of Distances between Space-Times"). MEDIUM.

11. `latremoliere_metric_st_2017` — title has spurious "dual" prefix (real: "The Gromov-Hausdorff propinquity for metric Spectral Triples"). MINOR (arXiv ID correct, authors correct, journal correct).

12. `kostant1965` — vol/year inconsistency (vol 75 was 1969, not 1965). MINOR. Not cited in body.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|:-----|:-----|:-----|:-----|:----|
| "To the author's knowledge this is the first quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input." | abstract | searched mondino_samann2025 §10 examples, arXiv recent (2025-2026), Ketterer/Kubota/Sakovich-Sormani/Che-Perales-Sormani follow-ups | None found. All extant MS-side examples are continuous (Chruściel-Grant, warped products). No operator-algebraic source identified. | Already correctly hedged with "to the author's knowledge." **Maintain as written.** |
| "the first quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input" | §5.6, §6.2 | same | none | **Soften to match abstract: "to our knowledge, the first..."** in both occurrences. |
| "no published paper constructs a categorical bridge between operator-algebraic Lorentzian propinquity and synthetic Lorentzian Gromov-Hausdorff frameworks" | §1.5 last paragraph | searched arXiv recent, math.OA + math.MG cross-listings | None found in either direction (Latrémolière side or Mondino-Sämann side). Adjacent works (Sormani-Vega, Sakovich-Sormani, Minguzzi-Suhr, Che-Perales-Sormani, Ketterer, Kubota, Mondino-Ryborz-Sämann) all stay within one or the other framework, none bridge. | Already hedged with "to the author's knowledge" earlier in same paragraph. Acceptable. |
| "the bridge constructed here is novel" | §1.5 last paragraph | same | same | Acceptable in this paragraph (preceded by "to the author's knowledge"). |
| "tenth math.OA standalone in the GeoVac series" | §1.5 | n/a | internal | OK — internal accounting claim. |

### Honest novelty note

Per the protocol's honest ceiling: a web search **cannot** confirm "first" — it can only confirm-false (which it does not here). The strongest defensible epistemic position is "to our knowledge, the first." Paper 48's abstract already uses this language; §5.6 and §6.2 should be brought into line. **No prior art was found despite genuine searches.**

## Combined severity table

| Finding | Pass | Verdict | Severity |
|:--------|:-----|:--------|:---------|
| `hekkelman_mcdonald2024b` arXiv:2411.04566 is BosonSampling paper, not spectral triples | B | CITE-MISATTRIBUTED | **HIGH** |
| `bertozzini_conti_lewkeeratiyutkul2009` arXiv:0901.4031 is tridiagonal-matrices paper | B | CITE-MISATTRIBUTED | **HIGH** |
| `leimbach_vs2024` authorship wrong (Leimbach alone, not Leimbach-vS); journal wrong | B | CITE-WRONG-METADATA | **MEDIUM-HIGH** (recurs corpus-wide) |
| `bizi_brouder_besnard2018` title AND journal wrong | B | CITE-WRONG-METADATA | **MEDIUM** |
| `mondino_samann2025_pointed` title wrong (but content verified) | B | CITE-WRONG-METADATA | **MEDIUM** (most load-bearing cite — fix matters) |
| `hekkelman_mcdonald2024` title wrong | B | CITE-WRONG-METADATA | MEDIUM |
| `ketterer2026` title wrong | B | CITE-WRONG-METADATA | MEDIUM |
| `minguzzi_suhr2024` title wrong | B | CITE-WRONG-METADATA | MEDIUM |
| `nieuviarts2025_v2` title wrong | B | CITE-WRONG-METADATA | MEDIUM |
| `sakovich_sormani2024` title wrong | B | CITE-WRONG-METADATA | MEDIUM |
| §5.6 and §6.2 priority claims ("the first...") not hedged to match abstract | A | C (overstatement) | MEDIUM |
| Bridge Theorem stated as "theorem-grade rigor" while paper contains only proof sketches | A | C (overstatement) | MEDIUM |
| `latremoliere_metric_st_2017` title has spurious "dual" prefix | B | CITE-WRONG-METADATA | LOW |
| `franco_eckstein2014` year inconsistency (key 2014 vs body 2013) | B | CITE-WRONG-METADATA | LOW |
| `mondino_ryborz_samann2025` title is short paraphrase | B | CITE-WRONG-METADATA | LOW |
| `kostant1965` vol/year inconsistency (vol 75 vs 1965) | B | CITE-WRONG-METADATA | LOW |
| `kubota2026` title missing "A" prefix (very minor) | B | LOW | LOW |

**Counts:**

- Pass A verdicts: A = 2 (claims 4, 11); B = 8 (1, 2, 3, 6, 7, 8, 9, 10); C = 2 (priority overstatement, theorem-grade overstatement); D = 0; E = 0.
- Pass B verdicts (32 bibitems): CITE-OK = ~18 (incl. assumed-OK standard refs); CITE-WRONG-METADATA = 11; CITE-MISATTRIBUTED = 2; CITE-DOESNT-SUPPORT = 0 distinct (covered by WRONG-METADATA on bizi_brouder_besnard); CITE-CANT-FIND = 0.
- Severity totals: HIGH = 2 (both CITE-MISATTRIBUTED); MEDIUM = 8; LOW = 5.

## Broadcast readiness: **YELLOW**

**Pass A**: Paper 48's structural content (functor construction, four bridge properties, seven wedge theorems) is internally consistent and the load-bearing externals (Latrémolière 2512.03573, Mondino-Sämann 2504.10380v4, Connes-Rovelli 1994) are real and on-topic. The substantive proof work all lives in the cited Phase-A memos, and the "theorem-grade rigor" claim in the abstract should be softened to match what the paper actually carries (proof sketches). The novelty claim ("first quantitative pLGH panel from operator-algebraic input") is honestly hedged in the abstract; matching language should propagate to §5.6 and §6.2. No mathematical errors found; the case-exhaustion proof in Lemma 4.10 is sound.

**Pass B**: Two CITE-MISATTRIBUTED entries (`hekkelman_mcdonald2024b` → BosonSampling paper; `bertozzini_conti_lewkeeratiyutkul2009` → tridiagonal matrices paper) are the Fursaev-Solodukhin failure mode that CLAUDE.md §3 explicitly warns about. Both are non-load-bearing citations (they appear in the "related work" section, not as proof inputs), so they would not corrupt the math, but they would embarrass on broadcast to a domain expert who clicks the arXiv IDs. Eight further CITE-WRONG-METADATA entries (titles wrong on Mondino-Sämann v4 — the most load-bearing cite of the paper — plus Bizi-Brouder-Besnard, Hekkelman-McDonald, Ketterer, Minguzzi-Suhr, Nieuviarts, Sakovich-Sormani, Leimbach-vS) form a systemic pattern suggesting the bibliography was assembled from memory/paraphrase rather than from actual arXiv page-fetches. **The Leimbach-vS authorship/journal error appears corpus-wide (per CLAUDE.md text) and should be fixed in the corpus, not just here.**

**Verdict: YELLOW.** Math is sound; bibliography needs a precise pass to fix 2 HIGH severity misattributions and ~8 MEDIUM title corrections before public broadcast. The novelty claim is correctly hedged in the abstract but should be propagated to §5.6 and §6.2.

## What I could NOT verify (hand to a human expert)

1. **Bridge Theorem proof-grade**: Paper 48 contains "proof sketches" only; the formal theorem-grade rigor lives in six Phase-A memos in `debug/`. A domain expert in math.OA needs to check the full proofs in those memos against Connes-vS 2021, Latrémolière 2025, and Mondino-Sämann v4 axioms.

2. **Functoriality**: Paper 48 claims $\Wfun$ is a (contravariant) functor on categories; the morphism check (compositions, identities) is asserted but not derived. A domain expert in category theory should verify the functor axioms.

3. **B4 convergence transport at the MS Def 3.6 level**: per-axiom verification (cardinality, distortion, extension, forward density) is sketched; the extension-property identification with Plancherel-nested Berezin is the most non-trivial step and warrants expert review.

4. **B2 Decomposition O off-orbit empty Case (iii)**: I verified this is correct logic for $M = 0$ topographic eigenvalue indexing $(N, L)$. But the structural-vacuum claim of off-orbit super-additivity at the K⁺-weak-form is the substantive content; whether it carries no information vs. is a feature of restriction needs expert read.

5. **Bibliographic priority of the bridge**: I searched recent arXiv (math.OA, math.MG, math-ph cross-listings 2025-2026) and found no prior art for a categorical Krein-PPQMS → MS-LPLS bridge. The "to our knowledge" hedge is appropriate; a domain expert may know of paywalled/unindexed work.

6. **Camporesi-Higuchi 1996 spectrum**: stated as integer $\mathrm{two\_m}_{j}$ via Paper 42. Not verified here.

## Top finding (single most important across both passes)

**Two CITE-MISATTRIBUTED entries (`hekkelman_mcdonald2024b` and `bertozzini_conti_lewkeeratiyutkul2009`) point to entirely unrelated arXiv papers — the Fursaev-Solodukhin failure mode CLAUDE.md §3 explicitly logs — and must be corrected before broadcast; the bibliography also has 8 further CITE-WRONG-METADATA title/author/journal errors (most load-bearing: the Mondino-Sämann v4 title is wrong) indicating systematic citation-from-paraphrase, while the underlying math and the load-bearing externals (Latrémolière 2512.03573, Mondino-Sämann v4 substantive content, Connes-Rovelli 1994) are real and correctly used.**
