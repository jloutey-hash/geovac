# Sprint Modular Propinquity — five-track synthesis

**Date:** 2026-05-23.
**Sprint position:** Synthesis of the 5-track parallel modular propinquity sprint (Track 0 lit audit + M-X Sturmian-FCI + M-Y NaH bimodule + M-Z Bethe log + M-H1 Higgs dual-modular). All five tracks landed in one session.

**Outcome (one line):** **The PI's modular bet pays substantively on the chemistry / bound-state QED side (3 of 4 application tracks return new content) but does not pay on the SM-unification side (M-H1 vocabulary-only). Three substantive new structural findings; one concrete implementation path with magnitude prediction.**

---

## §1. Track 0 literature audit — canonical references locked

Located six Latrémolière modular papers via WebSearch + primary-PDF extraction:

| arXiv | Year / Journal | What it provides |
|:------|:--------------|:-----------------|
| 1608.04881 | 2019 / *Dissertationes Math.* 544 | Modular propinquity, pseudo-metric on MQVBs (4-axiom definition: domination, closed unit-ball, inner + modular quasi-Leibniz) |
| 1811.04534 | 2021 / *J. Noncommut. Geom.* 15 | Dual modular propinquity, COMPLETE on MQVBs; relaxed MQVB Def 2.8 + metrical-QVB Def 2.12 |
| 1811.10843 | 2022 / *Adv. Math.* 404 | Spectral propinquity, natural fit for Camporesi-Higuchi substrate |
| 1703.07073 / 1803.06601 | 2017–2020 | Heisenberg modules over quantum 2-tori as MQVBs (structural parallel) |
| 2112.11000 | 2023 / *Math. Ann.* | Dirac spectrum continuity |
| 2512.03573 | Dec 2025 | Pointed-QLCMS (non-unital, non-compact, pin state — already deep-read in Phase A.2') |

**Critical correction Track 0 surfaced:** the "dual" in **dual modular propinquity** is **bridge ↔ tunnel construction duality**, NOT Hilbert-module duality, NOT Morita equivalence. My earlier discussion response to the PI used "Morita-equivalence-respecting" loosely — this was incorrect framing. The correct framing: "dual" refers to two complementary topology constructions (embeddings INTO a containing algebra via *bridges* vs surjections FROM a common algebra via *tunnels*); 1811.04534 proves completeness, which is the central reason for the dual construction.

**Roothaan-autopsy overlap verdict (PI's flagged question):** **SURFACE-LEVEL VOCABULARY SIMILARITY ONLY.** Track 0's six-row comparison: different "dualization" mechanisms, different decomposition structure (additive in autopsy; infimum-over-tunnels in math.OA), different output type (decomposed observable + residual in autopsy; single distance value in propinquity), different convention handling. The structures are complementary at *different layers*: Roothaan = Layer-2 itemization across Paper 34 projections; modular propinquity = Layer-1 spectral-triple substrate. Do not claim the Roothaan autopsy has a math.OA equivalent in propinquity vocabulary.

## §2. Three convergent threads across M-X / M-Y / M-Z

### Thread A: Sturmian basis is a MODULE basis (M-X + M-Z agree)

Both Sub-sprint X's verdict and Sub-sprint Z's bound-state QED reading independently identified that the Coulomb-Sturmian basis at $\lambda = Z/n_*$ is naturally a **module basis** (under the $1/r$-weighted Sturmian inner product) on the Hilbert C*-module $\Mcal = L^2(\R^3) \otimes \C^{2j+1}$ over the algebra $\Acal = C_0(\R^3) \otimes M_d$. The Sturmian projector $P_{n_{\max}}^{\mathrm{Sturmian}}$ truncates the *module* without changing the *algebra*; the algebra acts on the truncated module via natural compression.

**Why this matters:** the X/Z framing divergence from Sub-sprint X (Coulomb-mult algebra vs spectral algebra readings of where Sturmians live) **dissolves under the modular framework**. Both prior readings were forcing module-side structure into algebra-only language. The natural modular reading places Sturmians at the module level cleanly.

### Thread B: R1 gradient-Dirac is canonical, not a kludge (M-X)

Sub-sprint X's R1 workaround for Schrödinger Leibniz failure was: replace $L(f) = \|[D_S, M_f]\|$ (which fails Leibniz for second-order $D_S$) with $L(f) = \|\nabla f\|_\infty$. In the modular framework, this is **canonical**:

$$D_S = -\tfrac{1}{2}\nabla^2 - Z/r = \tfrac{1}{2}\nabla^* \nabla + V_{\mathrm{Coulomb}}$$

is the **Bochner Laplacian** of the gradient connection $\nabla$ on the trivial line bundle, plus multiplication by the potential. The modular D-norm lives at the first-order connection $\nabla$, not at the second-order Bochner Laplacian. Modular Leibniz $D(a\xi) \le \|a\|_\Acal D(\xi) + L_\Acal(a) \|\xi\|_\Mcal$ is immediate from the gradient product rule.

**Why this matters:** R1 was previously framed as a workaround for a structural defect (Schrödinger Leibniz failure). Under modular propinquity, R1 *is* the structural choice — the Schrödinger operator is correctly placed at the *second-order* Bochner Laplacian level (where it doesn't satisfy Leibniz, nor should it), while the *Lipschitz seminorm* is correctly placed at the *first-order* connection level (where Leibniz holds automatically). The framework's atomic FCI calculations have a uniform Latrémolière-compatible Lipschitz seminorm regardless of whether the operator of interest is Schrödinger or Dirac.

### Thread C: Multi-focal-composition wall has internal structure (M-Z headline)

Sub-sprint M-Z's substantive new finding: the multi-focal-composition wall (CLAUDE.md §1.7) splits cleanly into two pieces under the modular framework:

- **Bimodule cross-shifts (HANDLED by the framework, production code already composes them):**
  - HF-3 recoil (cross-register $V_{eN}$, `geovac/cross_register_vne.py`)
  - HF-4 Zemach radius (§III.18, `geovac/magnetization_density.py`)
  - HF-5 multi-loop $a_e$ (Parker-Toms first-order curvature, Paper 28 §curved_qed)
  - Finite nuclear size (Foldy-Friar, §III.17)
  - Hyperfine averaging
  
- **Module endomorphisms (NOT HANDLED — the LS-8a wall is precisely this piece):**
  - LS-8a renormalization counterterms ($Z_2 - 1$, $\delta m$)
  - Sprint H1 Yukawa selection
  - Inner-factor calibration data

**Why this matters:** previously the multi-focal wall was a single observation ("framework couples discrete labels cleanly; doesn't compose multiple Fock projections at once"). Under modular propinquity, the wall is a precise structural taxonomy: **the framework DOES compose bimodule cross-shifts (production code shows this empirically: HF-3 at sub-ppm, HF-4 at 0.012% of Zemach LO shift); it does NOT compose module endomorphisms (LS-8a renormalization). The LS-8a wall is the module-endomorphism piece, structurally distinct from cross-shift composition.** This sharpens the structural-skeleton-scope statement and makes the LS-8a follow-up sprint scope cleaner.

## §3. M-Y substantive new content — two-axis L/R decomposition

Sub-sprint M-Y returned the strongest single new-content finding of the sprint. Three substantive new diagnostics beyond Y's original structural ordering:

### (a) Two-axis bimodule distance reveals right-action dominance

The bimodule structure for NaH at the W1c level is:
- $\Acal_L$ = H-centered multiplication algebra (the H-valence carrier)
- $\Acal_R$ = Na-centered with FrozenCore structure (the Na-valence carrier)
- $\Mcal$ = valence bimodule (dim 5 at $\max_n = 2$)

The natural two-axis decomposition $d_\Mcal^\text{LR}(\xi_a, \xi_b) = (d_L^2 + d_R^2)^{1/2}$ applied to the four Y-candidates reveals:
- $d_\Mcal^\text{LR}((i), (iii)) \approx 0.67$ Ha
- $d_R / d_L \approx 6.7$ — **the W1c residual is right-action (Na-side wavefunction shape) dominant**
- H-side hydrogenic Z=1 is already correct; Na-side is the problem

### (b) PK was hitting the wrong axis

Track 3's Phillips-Kleinman cross-center barrier ($\Delta H_{pq}^{PK} = \sum_c (E_v - E_c) S_{pc} S_{cq}$) acts ONLY on the H-side valence via cross-overlaps against Na core. **In bimodule language: PK is a left-action shift parameterized by Na-core data.** But the W1c residual is right-action dominant. The 14.6% PK improvement is the subleading left-side correction; the dominant right-side correction needs a different mechanism. M-Y surfaces a new candidate **Path B (bilateral PK)** not in Track 3.

### (c) Magnitude prediction + alkali-hydride scaling

Predicted binding-recovery magnitude from the bimodule analysis: $\sim 0.1$ Ha (vs experimental NaH $D_e \approx 0.075$ Ha — right order of magnitude). Testable scaling prediction across the alkali-hydride series: $d_R \sim r_\text{valence}^\text{phys}(M) - 1.5$ bohr → LiH small (consistent: only binding alkali hydride in the catalogue), NaH/KH/RbH/CsH increasingly large.

**Why this matters:** M-Y converts the structural-skeleton-scope statement from "framework gives skeleton; calibration external" to "framework gives skeleton + two-axis decomposition + magnitude prediction + uniform scaling law across a chemistry series." The implementation work to close the NaH wall is now **bimodule-grounded** (Path A: one-sided Na-only physical wavefunction substitution; Path B: bilateral PK). Has a clean testable next-sprint deliverable.

## §4. M-H1 negative — PI's deeper bet on SM unification doesn't pay

Sub-sprint M-H1 tested the PI's hypothesis: does dual modular propinquity's duality requirement impose a non-trivial constraint on the Yukawa $Y$ that Sprint H1's algebra-level analysis missed?

**Verdict: NO — SELF-ADJOINTNESS-REPHRASING with one structural sharpening of the falsifier.**

Three findings:
1. The dual modular duality requirement reduces to Hermiticity of $D_F$ + matter-antimatter doubling, both already imposed by Connes' axiom system. No new constraint on $Y$.
2. The Higgs cross-block $\Phi$ lives entirely in the fiber sub-bimodule $\Acal_F = \C \oplus \mathbb{H}$, structurally decoupled from GeoVac index. **No GeoVac-side bimodule data exists to constrain it.** Consistent with the inner-factor Mellin engine reading (CLAUDE.md §1.7 WH1).
3. The Morita-triviality constraint ("could $\Phi = u^* v$ for GeoVac-derived $u, v$?") **does not survive matter-antimatter doubling** — and Track 0 audit independently confirmed that the "dual" in dual modular propinquity is NOT Morita anyway, so the Morita route was structurally barred.

**Net:** H1 POSITIVE-THIN verdict survives unchanged. G3 NEGATIVE and G4a POSITIVE-THIN survive. G4b cross-manifold obstruction (Paper 24 §V four-layer asymmetry) gets cleaner naming but no resolution.

**Why this matters:** clean information about where the modular framework's leverage lies. **The framework gains substantive new content on the chemistry / bound-state QED side (M-X, M-Y, M-Z), not on the SM-unification side (M-H1).** The Yukawa-undetermined verdict is now precisely characterized: it lives in the fiber sub-bimodule, has no GeoVac-side index data, and no propinquity-side construction can constrain it. This is consistent with the structural-skeleton-scope statement at the SM side AND with Sprint TD Track 3's empirical refinement ($Y_F > 0$ forced by electroweak phase transition at $T_C \approx 160$ GeV — empirical, not derived).

## §5. Net structural picture

Summarizing the modular framework's leverage on GeoVac's open questions:

| Open question | Modular framework verdict | New content |
|:--------------|:--------------------------|:------------|
| X/Z framing divergence (Sturmian) | DISSOLVED | Module-truncation reading, R1 canonical |
| Schrödinger Leibniz failure | RESOLVED | Bochner Laplacian of gradient connection |
| Non-compact Paper 38 $4/\pi$ rate | STRUCTURAL REASON IDENTIFIED | M1 Hopf-base mechanism has no Sturmian analog on $\R^3$; closed-form rate remains open multi-month frontier |
| W1c NaH binding wall | TWO-AXIS DIAGNOSTIC | Right-action dominant; PK wrong axis; ~0.1 Ha prediction; Path A/B candidates; alkali-hydride scaling |
| Multi-focal-composition wall (Sprint HF) | SPLIT INTO TWO PIECES | Bimodule cross-shifts handled; module endomorphisms not (LS-8a wall is the endomorphism piece) |
| Sprint H1 Yukawa $Y$ unspecified | NO NEW CONSTRAINT | Cross-block lives in fiber sub-bimodule, decoupled from GeoVac index |
| G3 chirality co-location | NEGATIVE SURVIVES | Two $\Z_2$'s on independent tensor factors |
| G4a Connes SM ($M_3$) | POSITIVE-THIN SURVIVES | Higher-value 1-2 month sprint (independent of dual modular) |
| G4b cross-manifold | OBSTRUCTION CLEANER-NAMED | Paper 24 §V fourth layer = bimodule-level modular-Hamiltonian obstruction; not solved |

Convergent reading: **modular propinquity provides the right structural vocabulary for GeoVac's bound-state-QED / atomic-FCI / chemistry work** (M-X canonical R1; M-Y bimodule diagnostic; M-Z multi-focal wall split). It does NOT provide new constraints on the framework's open SM-side questions (M-H1). This is informative: the framework's near-term physics-side leverage is on chemistry + precision atomic physics + bound-state QED, NOT on Connes SM electroweak unification.

## §6. Next sprint options ranked

Three options surfaced from the synthesis. Ranked by expected new content per week of work:

### Option α (HIGHEST EXPECTED VALUE): M-Y.1 implementation sprint — ~2 weeks

**Goal:** implement the M-Y two-axis L/R decomposition numerically for NaH and verify the predicted ratio $d_R / d_L \approx 6.7$. Open Track 3 Option 1 (Path A: one-sided Na-only physical wavefunction substitution) with bimodule-grounded justification. Test alkali-hydride scaling on LiH and KH.

**Predicted deliverable:** either binding-recovery closure of the W1c NaH wall at $\sim 0.1$ Ha (validating the bimodule prediction) OR a clean negative that sharpens the diagnostic (telling us which assumption to revise). Either outcome is high-value. Also: testable LiH/KH scaling that, if confirmed, lifts the bimodule diagnostic to a uniform tool across the chemistry catalogue.

**Why first:** M-Y is the only sub-sprint with a magnitude prediction. Implementation work is well-scoped (architecture exists: `screened_valence_basis.py`, `phillips_kleinman_cross_center.py`); the bimodule reading specifies which axis to address. Highest probability of closing an open wall.

### Option β (PAPER DRAFT): Light synthesis paper + Paper 18 / Paper 38 cross-references — ~1-2 weeks

**Goal:** Draft a math.OA-style synthesis paper (or extension of existing arc) capturing the three convergent threads (Thread A module basis, Thread B R1 canonical, Thread C cross-shift vs endomorphism partition). Could be a new standalone paper (Paper 48 candidate) OR a substantial extension of an existing paper (Paper 18 §III.7 master Mellin engine; Paper 38 cross-reference; Paper 32 §VIII bimodule remark).

**Predicted deliverable:** synthesis writeup; conservative paper-update set across §IV.6 (Paper 18 master Mellin engine extension); §V.D (Paper 34 multi-focal-wall sharpening — bimodule cross-shifts ARE composed, module endomorphisms are NOT); §1.7 WH1 multi-focal-wall taxonomy refinement.

**Why second:** locks the structural insights in writing before they decay. Conservative; no autonomy claim beyond what's already established. Paper drafts can run in parallel with Option α.

### Option γ (MULTI-MONTH FRONTIER): Non-compact modular propinquity rate theorem — ~3-6 months

**Goal:** Prove an explicit non-compact analog of Paper 38's $4/\pi$ universal rate for Sturmian truncations on $L^2(\R^3)$. Originally an open question after M-X verdict (Latrémolière 1608.04881 / 1811.04534 give existence of complete metric but no rate).

**Predicted deliverable:** a new math.OA standalone paper (Paper 49 candidate) on Sturmian convergence rate. High mathematical originality if it lands.

**Why third:** multi-month timeline; structurally distinct from the framework's near-term physics deliverables; comparable scope to Papers 47/40 — substantial NCG math investment. The question is open ("is there a Plancherel-weighted asymptote for Sturmian on $\R^3$?") and Latrémolière hasn't answered it.

### Cross-option observation

Options α and β can run in parallel (different work streams). Option γ is the only one requiring multi-month commitment; it should be deferred unless Option α produces a clean negative that requires fundamentally new math to close.

**Default recommendation:** Open Option α + Option β in parallel. Option α gives a concrete physics deliverable (NaH binding closure or sharpened diagnostic); Option β locks the structural findings in paper form. Total time ~2-3 weeks.

## §7. Paper edits suggested (conservative)

Three small-paste edits flagged for review:

1. **Paper 18 §III.7 master Mellin engine** — add one-paragraph cross-reference: "Structural reason why M1's $4/\pi$ rate does not transport to non-compact Coulomb (Sturmian basis): the Hopf-base measure mechanism is compact-Lie-quotient specific; Sturmian's continuum content dilutes Plancherel weight to 1. See `debug/sprint_modular_propinquity_mZ_bethe_log_memo.md`."

2. **Paper 34 §V.D running catalogue** — add new entry: §V.D.7 "Multi-focal-composition wall splits into bimodule cross-shifts (handled) vs module endomorphisms (LS-8a wall). Per Sub-sprint M-Z." Conservative, one-paragraph entry; refers to existing memos.

3. **Paper 32 §VIII.C addendum to Sprint H1 verdict** — add one-paragraph remark: "Bimodule reformulation of H1 falsifier (Sub-sprint M-H1, 2026-05-23) confirms POSITIVE-THIN verdict at the operator-system level. Higgs cross-block lives in fiber sub-bimodule $\Acal_F$, structurally decoupled from GeoVac index. No GeoVac-side data constrains $Y$ at the bimodule level." Locks the negative finding.

No paper edits applied in this sprint per the constraint that the sprint was a research / verification task, not an implementation task. PI to decide whether to apply.

## §8. CLAUDE.md updates suggested

Two CLAUDE.md additions flagged:

1. **§1.7 multi-focal-wall taxonomy refinement** — add one-paragraph clarification: "The multi-focal-composition wall splits into (a) bimodule cross-shifts (handled by production code: HF-3, HF-4, HF-5, etc.) and (b) module endomorphisms (LS-8a renormalization, Sprint H1 Yukawa, inner-factor calibration — NOT handled by the framework). Per Sub-sprint M-Z, 2026-05-23."

2. **§2 active development frontier** — add one-paragraph sprint summary capturing the three convergent threads + the M-Y two-axis prediction + the M-H1 negative.

## §9. Honest scope

This synthesis:
- IS a synthesis of five independent sub-sprint verdicts (Track 0 + M-X + M-Y + M-Z + M-H1)
- IS NOT a new structural result (the structural findings are from the sub-sprints individually)
- DOES surface the convergent threads across M-X / M-Y / M-Z (substantive)
- DOES propose three sprint options with magnitude estimates and risk assessment

**Confidence:**
- HIGH on Track 0 audit findings (six papers cited with arXiv numbers, primary-PDF extraction)
- HIGH on convergent threads (M-X / M-Z agree independently on Sturmian-as-module basis; M-X / M-Z agree independently on Schrödinger as Bochner Laplacian)
- HIGH on M-Y's two-axis L/R decomposition (concrete numerical estimate $d_R / d_L \approx 6.7$)
- MEDIUM on M-Y's magnitude prediction ($\sim 0.1$ Ha) — order-of-magnitude estimate from the Frame−Physical-3p differential, not a rigorous calculation
- HIGH on M-H1 negative (cross-block in fiber sub-bimodule, no GeoVac-side data)
- MEDIUM on the M-Z multi-focal wall split — the partition is structurally motivated but the framework verification is at the structural-sketch level, not a rigorous theorem

**Files:**
- `debug/sprint_modular_propinquity_synthesis_memo.md` (this memo, ~4500 words)
- Cross-references:
  - `debug/sprint_modular_propinquity_literature_audit.md` (Track 0)
  - `debug/sprint_modular_propinquity_mX_sturmian_memo.md` (M-X)
  - `debug/sprint_modular_propinquity_mY_pinstate_memo.md` (M-Y)
  - `debug/sprint_modular_propinquity_mZ_bethe_log_memo.md` (M-Z)
  - `debug/sprint_modular_propinquity_mH1_higgs_memo.md` (M-H1)
  - `debug/sprint_l3e_p3_synthesis_memo.md` (prior X/Y/Z synthesis at 2026-05-23 before the modular extension)

## §10. Session summary

**Sprint outcome:** the modular propinquity bet pays substantively on the chemistry / bound-state QED side (3 of 4 application tracks return new content) and gives clean information about its non-applicability on the SM-unification side (M-H1 vocabulary-only). Three substantive new findings: (i) Sturmian-as-module-basis dissolves the X/Z framing divergence; (ii) M-Y two-axis L/R decomposition gives bimodule-grounded NaH binding prediction at $\sim 0.1$ Ha; (iii) M-Z multi-focal-composition wall splits into bimodule cross-shifts (handled) vs module endomorphisms (LS-8a wall). PI's flagged Roothaan-autopsy vocabulary similarity is confirmed surface-level only.

Recommended: parallel Option α (M-Y.1 implementation, ~2 weeks) + Option β (paper draft + conservative paper edits, ~1-2 weeks).

Pending PI direction.
