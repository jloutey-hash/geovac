# Sprint Phase A.3' — Concurrent-Work Re-check Memo

**Sprint date:** 2026-05-24
**Author:** Concurrent-work re-check sub-agent (PM-dispatched, parallel to Phase A.3' bridge construction)
**Search window:** 2026-03-01 → 2026-05-24 (~3 months)
**Predecessor audits:** Phase A.1 (L3e-scoping, ~2026-05-22), Sprint 1E concurrent-work re-check (2026-05-24 morning)
**Method:** WebSearch + WebFetch against arXiv math.MG / math.DG / math.OA / math-ph / gr-qc, plus direct author-page checks (Mondino, Sämann, Latrémolière).

---

## Headline verdict by priority area

| Priority area | Verdict | Notes |
|:--------------|:-------:|:------|
| 1. math.MG Lorentzian GH (post-Mondino-Sämann 2504.10380 follow-ups) | **MINOR** | Two new follow-on papers (Kubota 2605.09101 coarea; Ketterer 2605.11271 ℓ-convergence) — synthetic-side internal extensions; no bridge to operator-algebra side. |
| 2. math.DG synthetic Lorentzian curvature / pre-length spaces | **MINOR** | Mondino-Ryborz-Sämann 2605.03172 (May 2026, 91pp) is the only Mondino paper in the window — synthetic-side internal, no bridge. Beran-Sämann 2204.09491 v4 (Feb 2026 revision) — no scope change. |
| 3. math.OA Lorentzian propinquity / Krein spectral triple convergence | **CLEAR** | Latrémolière has no new papers post-2512.03573. Full math.OA April 2026 + May 2026 listings audited — no Lorentzian-side propinquity work. |
| 4. math-ph + gr-qc modular flow / thermal time Lorentzian NCG | **CLEAR** | No new work in window. Earlier Bagarello arXiv:2502.12750 (Feb 2025, "Thermal time of NC Minkowski") not new and does not bridge to synthetic frameworks. |
| 5. Mondino direct check | **MINOR** | One paper in window: 2605.03172 (synthetic-side stability, no operator-algebraic content). |
| 6. Sämann direct check | **MINOR** | Same paper (2605.03172) plus Beran-Sämann revision; no bridge work. |
| 7. Sormani / Cavalletti direct check | **MINOR** | Che-Perales-Sormani 2510.13069 (v3 Feb 2026) timed-Hausdorff compactness, synthetic-side, no bridge. Cavalletti has no new papers in window matching the bridge target. |
| 8. McCann | **CLEAR** | No new papers in window matching scope. |
| 9. Vega | **CLEAR** | No new papers in window. |
| 10. Cross-citation (papers citing both Mondino-Sämann 2504.10380 AND Latrémolière 2512.03573) | **CLEAR** | No such paper found in any subfield. |
| 11. Twisted spectral triple / pseudo-Riemannian Dirac (Nieuviarts lineage) | **MINOR** | Nieuviarts arXiv:2512.15450 has a v2 revision (May 11, 2026) — same scope; still KO-dim 3 restriction documented in L3e Phase A.1 audit. |

**Overall verdict: CLEAR-WITH-MINOR-FLAGS.** No paper in the window constructs the bridge between operator-algebraic Lorentzian propinquity and synthetic Lorentzian Gromov-Hausdorff frameworks (the F2 forward-vs-reverse triangle mismatch resolution). The Phase A.3' sprint can proceed without scope or scoop concerns. Three minor items warrant citation acknowledgement in the merged Paper 48 bibliography.

---

## New papers found (in priority order)

### Tier 1 — Direct hits in the search window, worth assessing in detail

#### Ketterer arXiv:2605.11271 — "Convergence of Lorentzian spaces and curvature bounds for generalized cones"

- **Posted:** May 11, 2026
- **Authors:** Christian Ketterer (solo)
- **Subjects:** math.DG, math-ph, math.MG
- **Contribution:** Introduces "ℓ-convergence," a new notion extending Mondino-Sämann 2504.10380 Lorentzian GH convergence. Proves stability of timelike curvature bounds under measured ℓ-convergence; establishes convergence of generalized Lorentzian cones; precompactness theorems for smooth generalized cones with uniform Ricci bounds.
- **Relevance to A.3' bridge target:** **MINOR.** Ketterer's ℓ-convergence is an internal refinement of the Mondino-Sämann synthetic framework. **The paper does NOT cite Latrémolière 2512.03573** (verified via PDF fetch on full reference list). No operator-algebraic content. No metametric vs reverse-triangle discussion. The ℓ-convergence definition is a measure-theoretic strengthening of pre-length space convergence, NOT a metametric.
- **Recommendation:** Cite in merged Paper 48 bibliography as concurrent synthetic-side convergence-notion refinement. No scope change. **One sentence in related-work** noting that Ketterer extends Mondino-Sämann 2504.10380 internally on the synthetic side, parallel to the GeoVac Phase A.3' bridge construction operating across the synthetic/operator-algebraic divide.

#### Kubota arXiv:2605.09101 — "Lorentzian coarea inequality"

- **Posted:** May 9, 2026 (v1); revised May 19, 2026 (v3)
- **Authors:** Hikaru Kubota (solo)
- **Subjects:** math.DG / math.MG
- **Contribution:** Introduces "locally uniformly d-controlling map between Lorentzian pre-length spaces preserving diameters of causal diamonds." Establishes coarea inequality for Lorentzian Hausdorff measure. Derives covering lemma under "local causal enlargement property."
- **Relevance to A.3' bridge target:** **MINOR.** Synthetic-side measure theory on Lorentzian pre-length spaces. No operator-algebraic content. The diameter-of-causal-diamonds technology is internal to the Mondino-Sämann framework. Could in principle be cited if the merged Paper 48 needs Hausdorff-measure-on-causal-diamonds technology, but not load-bearing.
- **Recommendation:** Note in A.3' memo if Hausdorff measure on causal diamonds becomes relevant for the bridge construction; otherwise skip.

#### Mondino-Ryborz-Sämann arXiv:2605.03172 — "Stability of Synthetic Timelike Ricci Bounds under C⁰-Limits and Applications to Impulsive Gravitational Waves"

- **Posted:** May 4, 2026
- **Authors:** Andrea Mondino, Vanessa Ryborz, Clemens Sämann
- **Subjects:** math.DG, gr-qc, math.MG (91 pages)
- **Contribution:** Proves stability of synthetic TCD^e_p(K,N) curvature-dimension condition under C⁰-convergence of smooth Lorentzian metrics. Applies to impulsive gravitational waves arising from nonlinear interactions.
- **Relevance to A.3' bridge target:** **MINOR.** This is the Mondino-Sämann program's flagship May-2026 paper. **PDF audit confirms it does NOT cite Latrémolière 2512.03573 and contains no operator-algebraic content.** Convergence framework is C⁰-Lorentzian-metric convergence + optimal transport (Sturm-Lott-Cavalletti lineage). No engagement with the metametric vs reverse-triangle pre-length space distinction (assumes the reverse-triangle pre-length framework throughout).
- **Recommendation:** Cite in merged Paper 48 as the May-2026 state-of-the-art Mondino-Sämann-line paper, documenting that the synthetic-side program continues to expand internally without bridging to the operator-algebraic side. Confirms the bridge gap remains open at the literature level.

#### Che-Perales-Sormani arXiv:2510.13069 — "Convergence of Timed-Metric Spaces and Causality" (revised v3 Feb 26, 2026)

- **Posted:** Oct 15, 2025 (v1); revised v3 Feb 26, 2026
- **Authors:** Mauricio Che, Raquel Perales, Christina Sormani
- **Subjects:** math.DG, math.MG
- **Contribution:** Compactness theorem for intrinsic timed-Hausdorff convergence of timed-metric spaces. Triangle inequality for timed-Hausdorff distance. Introduces "addresses" as conceptual tool. Arzelà-Ascoli theorem for Lipschitz functions on converging metric spaces.
- **Relevance to A.3' bridge target:** **MINOR.** Synthetic-side; timed-metric framework distinct from both Mondino-Sämann pre-length spaces AND Latrémolière metametric. The v3 revision (Feb 2026) sits just inside the search window edge. This was already caught in Sprint L3b-2 pre-submission hardening (Paper 45 §1.4 cross-reference). No new content beyond what was already noted there.
- **Recommendation:** Already handled in Paper 45 hardening pass. No new action needed; if the A.3' bridge memo cites Mondino-Sämann synthetic GH lineage, mention Che-Perales-Sormani as the timed-metric alternative.

#### Nieuviarts arXiv:2512.15450 — "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" (v2 May 11, 2026)

- **Posted:** Dec 17, 2025 (v1); revised May 11, 2026 (v2)
- **Author:** Gaston Nieuviarts (solo)
- **Subjects:** math-ph, hep-th
- **Contribution:** Almost-commutative-framework synthesis of twisted spectral triple → pseudo-Riemannian structure morphism, as alternative to Wick rotation. v2 (May 11, 2026) is a revision of the December v1 (same scope; size dropped slightly 50KB → 48KB, suggesting tightening rather than expansion).
- **Relevance to A.3' bridge target:** **MINOR.** Already documented in CLAUDE.md §2 and L3e Phase A.1 audit: KO-dim 3 restriction blocks direct applicability to GeoVac's S³ truncation. v2 May 2026 revision does not change the KO-dim 3 restriction (confirmed via fetch — no new content listed). Sprint 1E re-check already caught this.
- **Recommendation:** Already covered. No action.

### Tier 2 — Adjacent but not in direct scope

#### Manzano-Mosani-Sämann-Zoghlami arXiv:2512.05842 — "Conformal transformations of metric spaces and Lorentzian pre-length spaces" (Dec 5, 2025)

- **Posted:** Dec 5, 2025 (just outside primary window, but in the broader Phase A audit envelope)
- **Authors:** Miguel Manzano, Karim Mosani, Clemens Sämann, Omar Zoghlami
- **Contribution:** First consistent notion of conformal length for low-regularity Lorentzian spaces; conformal invariance of angles and causality; extends Nomizu-Ozeki theorem.
- **Relevance to A.3' bridge target:** **CLEAR.** Pure synthetic-side conformal-transformation work. No operator-algebraic content. PDF audit explicitly confirms no Krein / Latrémolière / metametric citations.
- **Recommendation:** No action. The Phase A.1 audit picked this up structurally.

#### Erös-Gieger arXiv:2506.22197 — "A synthetic Lorentzian Cartan-Hadamard theorem" (June 27, 2025; revised Jan 21, 2026)

- **Posted:** June 2025; v2 revision Jan 21, 2026
- **Authors:** Darius Erös, Sebastian Gieger
- **Contribution:** Synthetic Cartan-Hadamard for Lorentzian; existence + uniqueness of timelike geodesics under global hyperbolicity + future one-connectedness.
- **Relevance to A.3' bridge target:** **CLEAR.** Pre-existing synthetic-side internal extension. No operator-algebraic content.
- **Recommendation:** No action.

#### Beran-Erös-Ohta-Rott arXiv:2509.26196 — "Concavity of spacetimes"

- **Posted:** Sep 2025 (revised in window)
- **Authors:** Tobias Beran, Darius Erös, Shin-ichi Ohta, Felix Rott
- **Contribution:** Local concavity of time separation functions on Finsler spacetimes (Lorentzian Busemann-convexity analog).
- **Relevance to A.3' bridge target:** **CLEAR.** Internal synthetic Finsler-Lorentzian; no operator-algebraic content.
- **Recommendation:** No action.

#### Braun-Sämann arXiv:2506.10852 — "Gromov's reconstruction theorem and measured Gromov-Hausdorff convergence in Lorentzian geometry"

- **Posted:** June 12, 2025 (revised Aug 16, 2025)
- **Authors:** Mathias Braun, Clemens Sämann
- **Contribution:** Lorentzian-signature Gromov reconstruction theorem; three convergence notions for normalized bounded Lorentzian metric measure spaces. Applications to causal-set theory spacetime reconstruction.
- **Relevance to A.3' bridge target:** **MINOR (already caught in Phase A.1).** PDF audit confirms no operator-algebraic content and no F2-mismatch engagement. This was caught in Phase A.1.
- **Recommendation:** Mention in A.3' memo only if Gromov reconstruction technology is load-bearing for the bridge construction; otherwise the Mondino-Sämann 2504.10380 + Latrémolière 2512.03573 pair remains the central reference axis.

### Tier 3 — Searched but no hits

- **arXiv math.OA April 2026 listing (50+ papers):** zero hits on Lorentzian propinquity, Krein spectral triple convergence, modular flow noncommutative Lorentzian, or any bridge to synthetic frameworks. Closest adjacent paper (Ponge 115pp, "Noncommutative Geometry, Spectral Asymptotics, and Semiclassical Analysis") is Riemannian-signature only based on its title and category. The math.OA April list is dominated by C*-algebra structure / groupoid / von Neumann / K-theory work, none Lorentzian.
- **arXiv math.OA May 2026 listing (65+ papers):** zero hits with same scope. Same Riemannian-side focus.
- **Latrémolière author page:** no new papers post-2512.03573 (Dec 2025). Confirmed by direct author-page fetch.
- **Sämann author page (per arXiv listing query):** only 2605.03172 (already covered) and Beran-Sämann 2204.09491 v4 Feb 2026 (no scope change) in the search window.
- **Mondino author page:** only 2605.03172 (already covered) in window.
- **Cavalletti, McCann, Vega, Sormani direct searches:** no new bridge-construction work in window. Cavalletti has nothing new in Lorentzian/synthetic in the search window beyond his ongoing Mondino collaboration. McCann: no new Lorentzian work in window. Vega: no new work in window.
- **Connes-Rovelli thermal time + synthetic Lorentzian cross-search:** Bagarello arXiv:2502.12750 (Feb 18, 2025; outside window) on κ-Minkowski thermal time. No new work in the March-May 2026 window connecting thermal-time framework to synthetic Lorentzian geometry. **The R3 candidate (thermal-time framework applied to synthetic Lorentzian geometry, flagged in Phase A.3' scope) is NOT scooped.**

---

## F2 forward-vs-reverse triangle bridge — explicit scoop assessment

**The central question:** Has any paper in the March 2026 → May 2026 window constructed the bridge between Latrémolière metametric (forward triangle, operator-algebraic) and Mondino-Sämann reverse-triangle pre-length space (synthetic, signature-flipped triangle inequality)?

**Verdict: NO. The bridge is not scooped.**

Evidence:

1. **Ketterer 2605.11271** is the closest May-2026 synthetic-side convergence-notion refinement; PDF audit confirms zero Latrémolière citations and zero metametric engagement. Ketterer's ℓ-convergence is a measure-theoretic strengthening internal to the Mondino-Sämann reverse-triangle framework — it does NOT cross the F2 mismatch.

2. **Mondino-Ryborz-Sämann 2605.03172** is the flagship May-2026 synthetic-side paper; PDF audit confirms zero Latrémolière citations and zero operator-algebraic content. Uses C⁰-Lorentzian + Sturm-Lott-Cavalletti optimal transport machinery, all synthetic-side.

3. **No new math.OA papers in April 2026 or May 2026** address Lorentzian propinquity, Krein spectral triple convergence, or Lorentzian noncommutative metric geometry. The Latrémolière metametric program (Dec 2025) has no follow-up by Latrémolière or co-authors in window.

4. **No paper in the window cites BOTH Mondino-Sämann 2504.10380 AND Latrémolière 2512.03573 simultaneously.** Direct cross-citation search returned zero hits.

5. **The Nieuviarts twisted-spectral-triple line (2512.15450 v2, May 11, 2026)** does not cross the bridge either. Its Wick-rotation-alternative scope is algebraic (almost-commutative framework), not metric-convergence-theoretic. KO-dim 3 restriction documented in L3e Phase A.1 audit still applies.

**Implication for Phase A.3' sprint:** the bridge target is open. The Phase A.3' agent can construct the bridge (or document why it cannot be constructed cleanly) without scoop concerns. The merged Paper 48 program retains its full priority on the F2 mismatch resolution.

---

## Differential update relative to Sprint 1E (2026-05-24 morning)

Sprint 1E found Latrémolière 2504.11715 (Apr 2025, Riemannian analytic continuity), Mondino-Sämann 2605.03172 (May 4, 2026, synthetic-side stability), and Nieuviarts 2512.15450 (Dec 2025 v1). This re-check **adds three substantive items**:

1. **Ketterer 2605.11271** (ℓ-convergence, May 11, 2026) — new synthetic-side convergence-notion refinement.
2. **Kubota 2605.09101** (Lorentzian coarea inequality, May 9-19, 2026) — synthetic-side measure-theoretic addition.
3. **Nieuviarts 2512.15450 v2** (May 11, 2026 revision) — same scope, KO-dim 3 restriction preserved.

None of the three address the F2 bridge target. The combined Sprint 1E + this re-check gives a complete picture of the search window with high confidence the bridge has not been constructed.

---

## Recommendations for Phase A.3' memo + merged Paper 48 bibliography

### To Phase A.3' memo (running in parallel):

- **No scope change required.** Bridge target remains open.
- **Acknowledge in concurrent-work section** that Ketterer 2605.11271 is a synthetic-side parallel convergence-notion refinement (one sentence, "Ketterer's ℓ-convergence sharpens the Mondino-Sämann 2504.10380 framework internally on the synthetic side, parallel to but distinct from the cross-framework bridge construction here").
- **No engagement with Kubota 2605.09101** required unless Hausdorff measure on causal diamonds becomes load-bearing for the bridge construction.

### To merged Paper 48 bibliography (when drafted):

Add the following five bibitems to the related-work / concurrent-work paragraph:

1. **Ketterer 2605.11271** — synthetic-side ℓ-convergence refinement, May 2026.
2. **Kubota 2605.09101** — synthetic-side Lorentzian coarea inequality, May 2026.
3. **Mondino-Ryborz-Sämann 2605.03172** — May 2026 flagship synthetic-side TCD stability paper.
4. **Braun-Sämann 2506.10852** — synthetic-side Gromov reconstruction in Lorentzian metric measure spaces, August 2025 (Phase A.1 inheritance).
5. **Nieuviarts 2512.15450 v2** — twisted-spectral-triple alternative to Wick rotation (May 2026 revision; KO-dim 3 restriction documented).

Recommended related-work paragraph framing (one paragraph, ~150 words): document that the synthetic side (Mondino-Sämann program + Ketterer + Kubota + Braun) and the operator-algebraic side (Latrémolière 2512.03573 metametric program) have continued to expand internally in parallel during the year leading up to the merged Paper 48, but no published work bridges the two. The F2 forward-vs-reverse triangle mismatch resolution constructed here is the first such bridge to our knowledge.

### To CLAUDE.md §2:

This re-check produced no findings that change WH1 PROVEN, the Paper 38/45/47 status, or the L3e roadmap. Standard one-liner sprint entry sufficient when the Phase A.3' sprint lands. No standalone update needed for this re-check.

---

## Honest scope and caveats

- arXiv indexing has a 2-3 day lag; papers posted May 22-24, 2026 may not be fully indexed at the time of this re-check.
- Some recent papers may be in author archives or institutional preprint servers not surfaced by search. The Latrémolière, Mondino, and Sämann author-page checks reduce this risk to near-zero for the three central lineages.
- "No bridge found" is a statement about the published literature in the search window. It is logically consistent with there being unpublished work in progress in another group; the Phase A.3' agent's parallel work would scoop or be scooped by any such hypothetical work depending on timing. This is a fundamental limit of concurrent-work re-checks.
- The math.OA April-May 2026 listing audit covered ~115 papers but relies on title-and-abstract scanning; a paper with a Riemannian-sounding title but Lorentzian content in body could in principle be missed. Risk judged low given the lineage-author direct checks.

---

## Sources

- Mondino-Sämann arXiv:2504.10380 — Lorentzian Gromov-Hausdorff convergence and pre-compactness
- Latrémolière arXiv:2512.03573 — Quantum Gromov-Hausdorff Hypertopology on pointed Proper Quantum Metric Spaces
- **Ketterer arXiv:2605.11271** — Convergence of Lorentzian spaces and curvature bounds for generalized cones (NEW, May 11, 2026)
- **Kubota arXiv:2605.09101** — Lorentzian coarea inequality (NEW, May 9-19, 2026)
- **Mondino-Ryborz-Sämann arXiv:2605.03172** — Stability of Synthetic Timelike Ricci Bounds under C⁰-Limits (May 4, 2026)
- **Manzano-Mosani-Sämann-Zoghlami arXiv:2512.05842** — Conformal transformations of metric spaces and Lorentzian pre-length spaces (Dec 5, 2025)
- **Braun-Sämann arXiv:2506.10852** — Gromov's reconstruction theorem in Lorentzian geometry (Phase A.1 inheritance)
- **Nieuviarts arXiv:2512.15450 v2** — Emergence of Time from a Twisted Spectral Triple (May 11, 2026 revision)
- **Beran-Sämann arXiv:2204.09491 v4** — Hyperbolic angles in Lorentzian length spaces (Feb 4, 2026 revision)
- **Erös-Gieger arXiv:2506.22197** — Synthetic Lorentzian Cartan-Hadamard theorem
- **Beran-Erös-Ohta-Rott arXiv:2509.26196** — Concavity of spacetimes (Finsler Lorentzian)
- **Che-Perales-Sormani arXiv:2510.13069 v3** — Convergence of Timed-Metric Spaces and Causality (Feb 26, 2026 revision)
- arXiv math.OA April 2026 listing — audited (50+ papers, zero Lorentzian/bridge hits)
- arXiv math.OA May 2026 listing — audited (65+ papers, zero Lorentzian/bridge hits)
- arXiv math.MG April 2026 listing — audited (50 papers visible, one timed-metric paper, zero bridge hits)
- arXiv math.MG May 2026 listing — audited (Kubota + Ketterer + Mondino-Ryborz-Sämann the three Lorentzian-relevant hits, all synthetic-side)
- Mondino author page (arxiv.org/a/mondino_a_1.html) — direct verification, only 2605.03172 in window
- Sämann author page (arxiv.org/a/saemann_c_1.html) — direct verification, only 2605.03172 + Beran-Sämann revision in window
- Latrémolière author page (arxiv.org/a/latremoliere_f_1.html) — direct verification, no new papers post-2512.03573
