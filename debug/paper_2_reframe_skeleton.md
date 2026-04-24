# Paper 2 Reframe Skeleton (Sprint A Target)

**Status (2026-04-18, Phase 2 complete):** Sprint A reframe edits applied to Paper 2. File grew from 1428 → 1596 lines (+168, +12%). All major Sprint A findings now integrated. Title and §VII Links table left unchanged pending PI review. Sprint B (g−2) plan pending Kluth-Litim cross-check completion and PI sign-off.

**Edits applied (2026-04-18):**
1. ✅ §II.A "Spectral-triple setting and published lineage" (new subsection, ~180 words) — Marcolli-vS lineage framing per α-LS finding
2. ✅ Abstract extension (+~100 words) — three Sprint-A structural features + projection-constant working hypothesis
3. ✅ Bibliography +4 entries — Marcolli-vS 2014, Perez-Sanchez 2024, Perez-Sanchez 2025 Comment, alpha_sprint_a_memo
4. ✅ §III.B "Triple selection at n_max=3" (new subsection, ~220 words) — α-MI triple coincidence finding
5. ✅ §III.D APS-shape paragraph (~180 words) — α-SP structural hint on (-) sign
6. ✅ §IV "Structural home of external π" (~140 words) — α-PI Hopf-measure identification with Eq. labeled `eq:pi_hopf_measure`
7. ✅ Open Question #3 π³α³ footnote (~240 words + numerics) — α-EB v2 + α-X narrow-scope insertion per agent recommendation

**Edits NOT applied (deferred to PI judgment):**
- Title change: kept as-is. Current "from Spectral Geometry of the Hopf Fibration" is descriptive; the implicit "from" weakness is now contextualized by the new §II.A. PI may revise.
- §VII Links table: kept as-is. The Sprint A findings give partial structural support to gaps but don't decisively flip any × to ✓; updating the table symbols is a judgment call.
- §VII underived choices list (5 items at lines ~927-936): not edited — the new sections in §III and §IV implicitly address gaps #1, #2 (sign pattern via APS-shape; π via Hopf measure), but the explicit list still reads as fully open. PI may want to softly update.
- Conclusion: kept as-is. Already ends with "conjectural" acknowledgment; the abstract carries the Sprint A framing.

**Verification needed:**
- LaTeX compile (no automated check yet — cross-refs and \cite keys verified by inspection only)
- ~~α-X v2 (Kluth-Litim S⁵ SD verification) still in flight~~ — **DONE 2026-04-18: MATCH.** Kluth & Litim 2020 (EPJ C 80:269, arXiv:1910.00543v2) Table 1 (scalar) + Eq. 13 transformation (Dirac via Lichnerowicz E_KL = −R/4 × dim_spin = 4) reproduce α-X's six sympy values exactly: scalar (1, 10/3, 16/3) on R=20 unit S⁵; Dirac (4, −20/3, 14/3). Two trivial normalization conventions reconciled (KL's b is volume-normalized; sign of Lichnerowicz E). α-X PARTIAL-NEGATIVE verdict stands. Paper 2 §IV.G footnote needs NO softening.

---
**Owner:** PM
**Date started:** 2026-04-18
**Governed by:** CLAUDE.md §1.7 (WH1, WH5), §13.5 (conjectural label protection for combination rule)

---

## 1. Current state assessment (2026-04-18)

Paper 2 v1.x (`papers/conjectures/paper_2_alpha.tex`, 1428 lines) has absorbed most of the Phase 4B-4I machinery:

- Three-homes theorem (§IV.A) — stated cleanly
- Three obstructions to unification (§IV.B) — documented
- Connection to rigidity theorems (§IV.C) — Papers 7/23/24 cross-referenced
- ζ(3) aside (§IV.D) — transcendental content flagged for Paper 18
- Closure of structural-decomposition program (§IV.E) — explicit
- Statistical validation (§V) — p = 5.2×10⁻⁹
- Circulant / Z_3 Hermiticity (§VI) — complete
- Links 1-5 assessment (§VII) — honest, Link 3 currently marked ×

The honesty is already there. The abstract already says "the combination rule remains conjectural." The paper already states what is and isn't established.

**What's missing, in order of Sprint A's work:**

1. Spectral-triple framing (WH1) — the paper uses "spectral data" and "spectral invariants" language but doesn't explicitly place itself in the Connes-Chamseddine almost-commutative spectral triple context. Adding this is a §II.A or §IV.B.5 move, ~1-2 paragraphs.

2. Projection-constant interpretation (WH5) — the working hypothesis that α is the conversion constant between three structurally distinct spectral regimes, not a quantity derivable from a single principle. Lives in the abstract and in a new paragraph near §IV.E closure.

3. Sign-pattern structural justification (Sprint A α-SP) — if positive, a new §IV.B.4 subsection "Structural origin of the +/- pattern." If negative, an honest paragraph in §IV.E saying "we checked this, the sign pattern is not forced by known heat-kernel rules."

4. m-independence characterization (Sprint A α-MI) — K(m=2,3,4,5) values + structural interpretation. New §IV.F or table extension to §III/§IV.

5. π-prefactor justification (Sprint A α-PI) — if positive, resolves §VII Link 2 gap #2. If negative, documented in §VII.

6. Error budget (Sprint A α-EB) — new §IV.G or §VII.A "Structure of the residual at 8.8×10⁻⁸."

7. (If Sprint B lands positive) electron g-2 from QED-on-S³ — new §VIII "Auxiliary prediction: anomalous magnetic moment." This is the promotion-to-core lever.

## 2. Target changes

### Title (tentative)

**From:** "The Fine Structure Constant from Spectral Geometry of the Hopf Fibration"
**To:** "α as a Three-Sector Spectral Coincidence on the Fock-Projected S³"
*or:* "Spectral Homes of α on the Fock-Projected S³: A Structural Decomposition"

Rationale: remove "from," which implies derivation. Current title claims more than current content.

### Abstract (target edit)

Keep the first three sentences (K formula, numerical agreement, p-value). Replace the fourth sentence block:

*Current:* "However, the combination rule K = π(B + F − Δ) remains conjectural: this is an empirical formula with structural support, not a first-principles derivation."

*Replace with:* "Each of B, F, and Δ has a canonical spectral home proved as a theorem (Theorem 1): B is the finite Casimir truncation of the scalar Laplace-Beltrami spectrum at m=3; F is the infinite Dirichlet series of the scalar Fock-shell degeneracy at the packing exponent s = d_max = 4, yielding F = ζ_R(2); Δ⁻¹ is the single-level Dirac degeneracy g_3^Dirac = 40. Nine mechanisms eliminated across Phases 4B-4I confirm these three objects do not share a common spectral generator. The combination rule K = π(B + F − Δ) remains conjectural at the level of the sum; we offer the projection-constant interpretation (α as the conversion factor between three structurally distinct spectral regimes) as the working hypothesis."

### Structure migration plan

No section deletion. Only additions and abstract/introduction edits:

- §I (Introduction): add one paragraph naming the projection-constant interpretation as the framing
- §II (Hopf Bundle): add §II.A "Spectral-triple setting" placing the bundle inside the GeoVac spectral triple (WH1 framing)
- §IV (The Formula and its Three Homes): extend with §IV.B.4, §IV.B.5, §IV.F, §IV.G from Sprint A outputs
- §VII (What Is and Is Not Established): update Link 3 entry based on Sprint A results (× → ◐ partial if sign-pattern matches)
- §VIII (new, Sprint B contingent): "Auxiliary prediction: electron g-2 on S³"
- §IX (new epilogue): "Conditions for promotion to core" — explicit decision rule

### Status migration

Conjectural → core, with the conjectural label migrating to a narrower scope:

- Paper 2 moves from `papers/conjectures/` to `papers/core/`
- The conjectural label stays on the combination rule itself (Link 3 in the §VII assessment)
- CLAUDE.md §13.5 hard prohibition updated: "never remove the 'conjectural' label from the combination-rule observation in Paper 2, even when the surrounding theorems are promoted"
- CLAUDE.md §6 paper inventory tier updated: Paper 2 moves to Core (always load)
- MEMORY.md unchanged (memory doesn't track paper tier)

Promotion trigger: Sprint A passes (sign-pattern or m-independence gives structural support for at least one of the five §VII underived choices) AND Sprint B either reproduces Schwinger at quoted precision OR gives honest diagnosed miss that doesn't falsify WH1.

If Sprint A passes and Sprint B misses cleanly, Paper 2 still strengthens enough to consider promotion — but the PI decides, not the PM.

If Sprint A fails (sign pattern is ad hoc AND m-independence shows m=3 is a single-point fit with no structural preference), hold at conjectural. The reframe still applies; the promotion doesn't.

## 3. Citation cleanup list (to run after reframe lands)

Papers that cite Paper 2 in its "derivation attempt" framing — reference wording needs updating to "structural decomposition" or "three-homes theorem" as appropriate:

- Paper 18 (exchange constants) — multiple refs; the spine of this paper is about α transcendental content
- Paper 21 (synthesis) — Sec VI on Hopf α
- Paper 25 (Hopf gauge structure) — multiple refs to Paper 2's (B, F, Δ) as "gauge-theoretic invariants"
- Paper 29 (Ramanujan Hopf) — references in §6 and conclusion
- Paper 30 (SU(2) Wilson) — §8 analog-of-Paper-2 framing
- README / Zenodo description — if they mention Paper 2
- CLAUDE.md §2 α-decomposition narrative — already mostly aligned; check for "derivation" language

Estimated: 8-12 touch-ups. Can be batched into a single cleanup pass as a sub-agent dispatch after Paper 2 reframe is finalized.

## 4. Integration points with Paper 18 v2.0 restructure (WH2)

Paper 18 v2.0 will promote the Seeley-DeWitt + ζ-invariant decomposition from observation to theorem-with-evidence. That content overlaps Paper 2's §IV.D (ζ(3) aside) and should:

- Paper 18 takes the three-axis grid as its thesis (operator order × bundle × vertex topology)
- Paper 2 references Paper 18 v2.0's grid when discussing B/F/Δ placement
- Avoid duplication: Paper 2 mentions placement briefly, Paper 18 develops it
- Planned as a coordinated paper-pair sprint after Sprint A lands

## 5. What this skeleton is NOT

- Not a patch to apply now
- Not a commitment to specific wording
- Not a decision to promote (PI decides after Sprint A/B)
- Not a replacement for the Sprint A memos — those generate the actual content
- Not a plan for new computations — Sprint B (g-2) lives in a separate plan doc

---

## 6. Sprint A integration checklist (in-progress)

Two of five Sprint A tracks complete as of 2026-04-18.

### α-PI (π-prefactor) — POSITIVE PARTIAL

Source memo: `debug/alpha_sprint_a/pi_prefactor_memo.md`

**Finding:** External π in K = π(B + F − Δ) identified as the Hopf-bundle topological measure Vol(S²)/4 = Vol(S³)/Vol(S¹) = π. Three equivalent readings (base solid-angle per chirality, spinor half-period of fiber, bundle-volume per fiber-length) collapse by the Hopf submersion identity Vol(S³) = ½·Vol(S²)·Vol(S¹). Independent heat-kernel candidate a₀² = π. WH1-natural split: external π = Hopf base measure, internal π² inside F = ζ(2) = fiber zeta.

**Integration into Paper 2 reframe:**
- New §II.A subsection "Hopf-bundle measure factor" (~1 paragraph + the submersion identity).
- §VII Link 2 updated: gap #2 ("Why π as the overall prefactor?") moves from unresolved to *partially resolved* (structural candidate identified; derivation status still pending CC spectral-action computation).
- §IV.5 or new §IV.C.2 "Sector-by-sector π content" paragraph: external π (Hopf measure, first-order sector) vs internal π² (F = ζ(2), second-order Jacobi-θ sector) vs rational sectors (B, Δ). Matches Paper 28 T9 taxonomy.
- Promotion-readiness update: this resolves one of the five §VII underived choices.

### α-MI (m-independence) — NEGATIVE-CONFIRMED-BUT-SHARPENING

Source memo: `debug/alpha_sprint_a/m_independence_memo.md`
Data: `debug/data/alpha_m_independence.json`
Code: `debug/compute_alpha_m_independence.py`

**Finding:** K(m) grows as Θ(m⁵); no series structure; m = 3 is a single-point coincidence, NOT a partial sum of a convergent object. **However**, m = 3 is a *triple* coincidence: (a) B(m)/N(m) = dim(S³) = 3 selection principle (already known, Paper 2 §III), (b) the (+,+,−) sign pattern uniquely hits 1/α by 3.5 orders of magnitude over the next-best pattern, (c) the two independent canonical forms of Δ (Paper 2 boundary product |λ_m|·N(m−1) vs Phase 4H Dirac degeneracy g_m^Dirac/2) agree ONLY at m = 3. Three independent mechanisms converging at one integer is new structural support.

**Also:** F is m-independent (depends on packing axiom d_max = 4, not truncation cutoff). B and Δ are m-dependent. The three K-ingredients live on two distinct cutoff variables.

**Integration into Paper 2 reframe:**
- New §IV.F "The triple selection at m = 3" subsection. Tabulate all three mechanisms and note that they converge at the same integer.
- Sign-pattern table (8 sign patterns at m=3 with rel errors) added as a new table.
- Two-Δ-forms-agree-only-at-m=3 noted as a new structural finding.
- F's m-independence vs B, Δ m-dependence documented in a short paragraph (links to d_max packing axiom from Paper 0).
- §VII updated: Link 2 gap #1 ("Why additive B + F − Δ") strengthened by the sign-pattern uniqueness at m=3.
- Promotion readiness: strong positive signal — the triple selection softens the "single-point coincidence" worry.

### α-SP (sign-pattern) — HONEST NEGATIVE (with weak APS-shape partial)

Source memo: `debug/alpha_sprint_a/sign_pattern_memo.md` (~2200 words)

**Finding:** Standard Connes-Chamseddine bulk expansion gives all-positive signs (+,+,+) on S³, NOT the (+,+,−) of K. The only standard heat-kernel rule producing a (−) next to bulk (+) is the Atiyah-Patodi-Singer (APS) eta-invariant boundary subtraction. Paper 2's Δ is *shape-compatible* with APS-style boundary subtraction (Δ⁻¹ = g_3^Dirac is a single-level Dirac mode count, in the right structural family) but is NOT literally an APS eta-invariant. No spectral-action functional generates K as its leading three terms.

**Decisive structural mismatches:** (1) Paper 2 has no cutoff function f or moments f_k; (2) the three pieces live in different operator sectors at different degrees of regularization (finite Casimir trace, ζ-regularized infinite sum, single-level degeneracy); (3) Phase 4E α-I already eliminated the closest spectral-action candidate (S⁵ SD coefficients). 

**Sources consulted (now in Paper 2 reference pool):** Connes-Chamseddine 1997, Chamseddine-Connes 2010 ("Uncanny Precision"), Vassilevich 2003 (heat kernel manual), van Suijlekom 2024 (textbook 2nd ed), Connes-van Suijlekom CMP 2020 (spectral truncations), Estrada-Gracia-Bondía 2011, APS 1975.

**New candidate framework flagged:** Connes-van Suijlekom "spectral truncations" (CMP 2020, arXiv:2004.14115) — finite-resolution spectral truncation with a "tolerance" remainder term. Not yet tested for Paper 2's Δ. Recommended as candidate for any future Sprint addressing Link 3.

**Integration into Paper 2 reframe:**
- New §III final paragraph: honest acknowledgment that (+,+,−) is shape-compatible with APS-flavor boundary subtraction but NOT derivable from any specific spectral-action functional. Quote the agent's recommended paragraph wording verbatim (good calibration).
- New §IV.G "Open question: spectral truncation": optional one-paragraph note pointing to Connes-van Suijlekom as a candidate for any future α-program reopening.
- §VII Link 2 gap #1 ("Why additive combination?"): partially addressed — the additive structure has a partial structural home (APS shape), not derivable.
- WH1 STATUS UPDATE in CLAUDE.md §1.7: from "shape-match only, no explicit derivation" to "shape-match for π-prefactor (α-PI POSITIVE PARTIAL) + APS-shape for sign pattern (α-SP WEAK PARTIAL); no derivation; literature survey for discrete spectral triples pending (α-LS hit usage cap)."

### α-EB v2 (error budget) — POSITIVE SUBSTANTIVE (with cross-paper hint)

Source memo: `debug/alpha_sprint_a/error_budget_memo_v2.md` (~840 words)
Data: `debug/data/alpha_error_budget_v2.json`
Code: `debug/alpha_sprint_a/compute_error_budget_v2.py`

**v1 stalled at watchdog after partial finding; v2 verified the partial at 80 dps and corrected one definitional conflation.**

**Key verified numbers (80 dps):**
- K − 1/α_codata = 6.5330×10⁻⁵ ≈ α² (trivial cubic offset, as expected from the cubic α³ − Kα + 1 = 0)
- R_predict = K − 1/α − α² = 1.2079×10⁻⁵ (substantive post-cubic residual)
- π³ · α³ = 1.2049×10⁻⁵
- **R_predict / (π³α³) = 1.00251** — within 0.25% of unity at 80 dps

**π³ uniquely picked out:** π² gives ratio 3.15, π⁴ gives 0.32. Only π³ is order unity.

**Structural reading:** π³ = Vol(S⁵). Connects Paper 2's residual to Paper 24's Bargmann-Segal S⁵ lattice. Three-loop QED reading is ruled out by sign (3-loop QED gives α³/π³, opposite sign).

**Next-order coefficient:** R_predict = π³α³(1 + Cα + ...) gives C ≈ 0.344. Closest closed forms: π/9 (1.4% gap, preferred), 1/3 = 1/dim(S³) (3.2% gap, suggestive). Single CODATA data point cannot uniquely fix C.

**Integration into Paper 2 reframe (post α-X cross-check, 2026-04-18):**
- α-X verdict is PARTIAL (negative-leaning). π³ shape-match survives at the Seeley-DeWitt coefficient level (every S⁵ SD coefficient factorizes as π³ × rational; scalar a₀ = π³, a₂ = 10π³/3, a₄ = 16π³/3; Dirac variants have analogous factors). BUT no S⁵ spectral-action coefficient produces a contribution at order α³ structurally. The CC expansion on d=5 is in odd powers of Λ, and standard CC uses a₀/a₄ ratios to fix α once, not α³.
- **DO NOT integrate α-EB v2 as a substantive §IV.G "Structure of the residual" section** in Paper 2. The π³α³ match is a structural hint, not a derivation, and the S⁵ machinery does not produce it at the right order.
- **DO integrate as a narrowly scoped footnote** in Paper 2's Open Questions or §VII. Suggested framing: "the post-cubic residual matches π³α³ to 0.25%; π³ appears naturally in S⁵ Seeley-DeWitt coefficients, but no specific S⁵ spectral-action coefficient produces the α³ structure; this is a suggestive numerical pattern, not a derivation."
- Cross-reference to Paper 24 in this footnote is appropriate but NARROW: "Vol(S⁵) = π³ is the integration measure of S⁵ heat-kernel traces; Paper 24's π-free certificate on the discrete Bargmann-Segal graph is consistent with this — the π-free claim applies to the discrete spectrum and edge weights, NOT to continuum integration measures."
- Paper 24 itself stays unchanged. No new cross-paper section needed.
- §VII Link 3 entry: still ×, unchanged. The α-EB v2 finding is too thin to move it to ◐.

**Cross-paper hint status:** the Paper 24 ↔ Paper 2 connection via π³ = Vol(S⁵) is real at the shape level (SD coefficients factorize) but not derivational (no coefficient produces α³). Does NOT reshape Sprint B priorities — the agent's original plan (g−2 from QED-on-S³) stands.

**WH implication:** WH5 (α as projection constant) is *slightly strengthened* by this cross-check: even after extending the machinery to S⁵, no common generator appears. WH1 is neutral (the SD calculation is consistent with the spectral-triple framing, but the lack of α³ coefficient means no new positive evidence either). WH4 (four-way S³ coincidence) is neutral — the S⁵ connection is a shape-only extension that doesn't promote the coincidence to a five-way pattern.

### α-LS (literature survey) — MODERATE (re-dispatch succeeded)

Source memo: `debug/alpha_sprint_a/spectral_triple_literature_memo.md` (16 verified entries)

**Finding:** WH1 literature grounding is MODERATE — each structural ingredient of the spectral-triple framing has published precedent.

**Closest published match to GeoVac:** Marcolli-van Suijlekom 2014, "Gauge networks in noncommutative geometry" (J. Geom. Phys. 75, arXiv:1301.3480). Finite spectral triples on graph vertices with connection on edges; spectral action = Wilson lattice gauge theory. Papers 25 (U(1)) and 30 (SU(2)) directly fit this template. **This is close enough that Papers 25 and 30 should cite it explicitly.**

**Important correction:** Perez-Sanchez 2024 (arXiv:2401.03705, Comment arXiv:2508.17338) partially corrects Marcolli-vS — the continuum limit is Yang-Mills *without* Higgs, not YM-H as originally claimed. Any GeoVac citation of Marcolli-vS should pair with Perez-Sanchez.

**Published precedent for each WH1 ingredient:**
- Finite ACG classification: Krajewski 1998 (hep-th/9701081), Paschke-Sitarz 2000, Ćaćić 2009 (arXiv:0902.2068)
- Graph spectral triples: Marcolli-vS 2014, de Jong 2009 (arXiv:0904.1291)
- Spectral truncations / operator systems: Connes-van Suijlekom CMP 2021, Hekkelman 2022+2024
- SM spectral action: Chamseddine-Connes 1997 (hep-th/9606001), 2010 (arXiv:0812.0165)
- van Suijlekom textbook 2024 2nd edition (open access)

**Critical finding:** NO published framework predicts a physical constant like α from a discrete spectral action. This means GeoVac's Paper 2 claim (α from K = π(B + F − Δ) on the Fock-projected S³) is *genuinely novel* and not duplicated in the literature. Paper 2's ambitious claim remains GeoVac's own.

**Red flags raised:**
1. Marcolli-vS continuum-limit claim is overstated; Perez-Sanchez correction must accompany any citation.
2. Published spectral-truncation examples (Hekkelman) are continuous-base-with-cutoff (circle → Toeplitz), NOT graph-native like GeoVac. GeoVac is more discrete than the canonical spectral-truncation setting.
3. No published α-from-spectral-action precedent: WH1's strongest form (α derivable from CC spectral action on the right triple) *exceeds* the literature. Framing should acknowledge this.

**α-EB v2 / α-X follow-up lead:** Kluth-Litim EPJ C 80:269 (2020) (arXiv:1910.00543) gives explicit Seeley-DeWitt coefficients on spheres in all dimensions including odd. α-LS agent couldn't extract the PDF; deferred to a short follow-up. This would independently verify α-X's sympy-computed S⁵ SD values.

**Integration into Paper 2 reframe:**
- Add Marcolli-van Suijlekom citation to Paper 2's references as the published framework closest to GeoVac's graph-spectral-triple reading.
- Add Perez-Sanchez 2024 alongside as the correction.
- Paper 2's §II.A "Spectral-triple setting" should reference these as "GeoVac is a specific graph-spectral-triple in the Marcolli-vS lineage."
- Do NOT overstate: the α prediction remains GeoVac's claim, not derivable from Marcolli-vS as written.

**Cross-paper updates flagged (separate touch-up sprint):**
- Paper 25 (Hopf gauge): add Marcolli-vS 2014 + Perez-Sanchez 2024 citations. Paper 25 already looks like Marcolli-vS but doesn't cite it.
- Paper 30 (SU(2) Wilson): same — add Marcolli-vS 2014 + Perez-Sanchez 2024.
- Paper 24 (Bargmann-Segal): optional add of Kluth-Litim 2020 if we verify S⁵ SD values in a follow-up.

**WH1 STATUS UPDATE:** from "partial shape-match" to **MODERATE**. Concrete evidence: each WH1 structural ingredient has published precedent. The ambitious part of Paper 2 (predicting α) remains novel. CLAUDE.md §1.7 WH1 Status updated accordingly.

Three agents still running in background. Will update this section as they complete.

### Running assessment of promotion readiness

After α-PI + α-MI + α-SP (3 of 5 tracks complete; α-EB v2 running, α-LS deferred):

- Link 1 (Hopf as arena): ✓ unchanged
- Link 2 gap #1 (additive combination): partially addressed — APS-shape gives qualitative structural origin; sign-pattern uniqueness at m=3 (α-MI) gives finite-N selection support. Still no derivation.
- Link 2 gap #2 (π prefactor): partially resolved — Hopf-measure structural candidate identified (α-PI Vol(S²)/4 = π).
- Link 2 gap #3-5 (ζ(2) for fiber, 1/(|λ_3|N(2)) for boundary, det M = −1): not yet touched by Sprint A.
- Link 3 (combination rule): still × (not established). Constraints have tightened — three independent m=3 selections converge (α-MI), the (+,+,−) sign pattern is shape-compatible with APS (α-SP), and the π prefactor has a Hopf-measure home (α-PI). But no spectral-action functional generates K (α-SP decisive structural mismatch).
- Link 4 (cubic from Z_3 circulant): ✓ unchanged
- Link 5 (α from cubic): ✓ unchanged

**Net after three tracks:** Paper 2's framing strengthens from "empirical formula with structural support" to "three-mechanism selection at m=3 with Hopf-measure π-prefactor and APS-shape sign pattern." That's a nontrivial upgrade. **However:** α-SP delivered an honest negative on direct spectral-action derivation, which means WH1 in CLAUDE.md §1.7 should be downgraded from "strong shape-match" to "partial shape-match — π-prefactor and three-objects pattern fit; sign pattern requires APS-flavor subtraction not derivable from CC alone."

**Implication for promotion decision:** Sprint A's net effect is to *narrow* the open question from "where do B, F, Δ come from + why this combination" (where most of the original conjectural mass sat) to "why does the APS-shape combination at the unique B/N=3 selection point hit α⁻¹ to 8.8×10⁻⁸." The narrower question is the genuine residual mystery. Promoting Paper 2 to core under this framing would mean: the paper is core-grade as a *structural decomposition with named open question*, not as a derivation. The promotion test (Layer 2 in the original plan: out-of-sample prediction via electron g−2) becomes more important now, not less, because Layer 1 internal consistency landed mixed (3 partials, no decisive positive).

**WH1 status update (CLAUDE.md §1.7):** I recommend updating WH1's "Status:" field after α-EB v2 lands and before any further reframe action. Current language "shape-match only, no explicit derivation of K from a spectral-action functional" was already accurate; α-SP confirms it; α-PI partially upgrades the π-prefactor cell; α-MI doesn't change WH1 directly but strengthens WH5.

**WH5 status update (CLAUDE.md §1.7):** WH5 ("α as projection constant") is *strengthened* by Sprint A. The triple m=3 selection (α-MI) and the α-SP decisive negative on common-generator spectral action both support the projection-constant reading. WH5's "Status:" field should be updated to note the new support.
