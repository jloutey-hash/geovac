# Sprint C-arc closure — E6 + D5/D6 + chemistry η-trivialization analog (2026-06-08)

## TL;DR

**Three theorems / structural targets added to Paper 32 §VIII in one bundled closure sprint, completing the path-#1 vein C arc.** Paper 32 §VIII now contains the C-arc terminal-state remark (`rem:c_arc_terminal_state`) tallying eight theorem-grade non-selection results across four sectors:

1. **E6 (combination rule K):** `thm:no_single_mechanism_K` + `rem:e6_scope`. No single morphism in $\mathcal{A}$ generates $K = \pi(B + F - \Delta)$. Three spectral homes (B Casimir, F Fock Dirichlet, Δ Dirac mode count) live in two Mellin sub-rings; twelve mechanisms empirically eliminated. Combination rule for $\alpha^{-1}$ remains conjectural per §13.5; theorem strengthens empirical non-derivability into structural impossibility.

2. **D5/D6 (cutoff function):** `thm:cutoff_function_external` + `cor:cc_fine_tuning_external` + `rem:d5_d6_scope`. Cutoff function $f$ and its Mellin moments $\varphi(s)$ are external test-function data; no axiom in $\mathcal{A}$ selects moment values. Wald forces RELATIONS among moments; values are external. CC fine-tuning $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ is external moment selection.

3. **Chemistry-side η-trivialization analog:** `rem:chemistry_eta_analog`. Named multi-month NCG-research target: a structural property of the FrozenCore $Z_{\rm eff}(r)$ pipeline mirroring $\{\gamma_F, D_F\} = 0$ on the chemistry inner-factor side. Sprint-scale path unclear; closure is by documentation, not by theorem. H6, H7 remain empirical-only pending.

Plus a final closure remark **`rem:c_arc_terminal_state`** tallying the full corpus state: eight theorem-grade non-selection results across inner-factor / multi-focal-composition / α combination / gravity sectors, one multi-month named target (chemistry analog), all subsidiary residuals (F4-F7, I-entries) accounted for.

Paper 32 compiles three-pass clean at 83 pages (up from 80 in v3.102.0). The C-arc is in stable terminal state.

## Sprint context

The user wants to "wrap up these 3 before I reset context." The framing is closure of the C-arc, not new research. E6 and D5/D6 are sprint-scale theorem upgrades; the chemistry analog is a multi-month target documented at the structural-pointer level.

This sprint is therefore three deliverables of distinct character:
- **E6 theorem:** sprint-scale upgrade, follows C3/F3 pattern (crystallise existing structural content into Theorem block)
- **D5/D6 theorem:** sprint-scale upgrade, parallels the inner-factor input-data theorems on the gravity sector
- **Chemistry analog:** scoping/documentation closure, names the multi-month target precisely

The bundling is honest: all three are closures of named follow-on targets from the v3.99.0-v3.102.0 sprint memos.

## Theorem 1 — E6 combination rule K

**Statement (paraphrase, see `thm:no_single_mechanism_K`).** Under $\mathcal{A}$ and the master Mellin engine case-exhaustion theorem, the three spectral homes that appear in the empirical combination rule $K = \pi(B + F - \Delta)$ for $\alpha^{-1}$ (Paper 2) live in two distinct Mellin sub-rings:
- $B = 42$ — Casimir trace at $m = 3$, $M_2$ sub-ring
- $F = \pi^2/6$ — Fock-degeneracy Dirichlet at $d = 3$, $M_2$ sub-ring
- $\Delta = 1/40$ — Dirac mode degeneracy $g_3^{\rm Dirac}$, $M_1$ sub-ring

No morphism in $\mathcal{A}$ produces the combination $K = \pi(B + F - \Delta)$ as an output of one Mellin mechanism or a composition of two/three with the specific weights $(+1, +1, -1)$ and external $\pi$ prefactor. Twelve mechanisms empirically eliminated (Phases 4B-4I + Sprint A α-SP + Sprint K-CC).

**Proof structure.** Composition of three facts:
1. Master Mellin engine case-exhaustion theorem (Paper 32 §VIII `thm:pi_source_case_exhaustion`): every transcendental classified into $M_1, M_2, M_3$; sub-rings don't inter-translate within $\mathcal{A}$.
2. Three independent spectral homes for $B, F, \Delta$ (Paper 2 §IV + Phases 4B-4I).
3. Twelve mechanism eliminations (Phases 4B-4I: nine candidates; Sprint A α-SP: CC direct spectral-action derivation; Sprint K-CC: PSLQ at 100 dps + natural Λ search + unified CC heat-kernel blocked by T9).

**Honest scope per §13.5.** The theorem identifies the structural reason for the 12-mechanism elimination but does NOT derive $K$. The conjectural label on the combination rule is preserved (per CLAUDE.md §13.5). Paper 2's status as Observations-folder paper (post-2026-05-02 curve-fit-audit migration) is unchanged. Sibling axioms producing $K$ from a single mechanism remain logically possible but unexhibited.

## Theorem 2 — D5/D6 cutoff function

**Statement (paraphrase, see `thm:cutoff_function_external`).** The CC spectral action $S(D, \Lambda, f) = \mathrm{Tr}\, f(D/\Lambda)$ takes $f$ as a Schwartz-class test function not constrained by $\mathcal{A}$. The Mellin moments $\varphi(s) = \int_0^\infty x^{s-1} f(x) dx$ are external calibration data. Wald's two-term-exactness $\Rightarrow$ pure Einstein $\Rightarrow$ action-$G$ $\equiv$ entropy-$G$ argument (Paper 51 §G7) forces RELATIONS among the moments via the sector-wise Mellin moment map (tip $\leftrightarrow \varphi(0)$, Einstein-Hilbert $\leftrightarrow \varphi(1)$, $\Lambda_{cc} \leftrightarrow \varphi(2)$); the relations are theorem-grade within $\mathcal{A}$. The moment VALUES are not constrained.

**Corollary (`cor:cc_fine_tuning_external`).** The CC problem $\Lambda_{cc} \cdot G_{\rm eff} \approx 36\pi \cdot \varphi(2)/\varphi(1)^2 \approx 10^{-122}$ requires $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ with $\varphi(0), \varphi(1) = O(1)$. No natural cutoff function produces this fine-tuned ratio; no axiom in $\mathcal{A}$ constrains the moment values to it.

**Honest scope.** The theorem formalises the precise location of the CC problem in the framework: gravitational *relations* are forced (Wald, theorem-grade); cutoff-function *moments* are external. The CC problem is a moment-selection problem on $f$, structurally similar to the Higgs-Yukawa selection problem on a different calibration object. Sibling-axiom direction (RG-flow, regulator independence, asymptotic safety) named but unexhibited.

## Closure 3 — Chemistry-side η-trivialization analog

**Scoping closure (`rem:chemistry_eta_analog`).** The SM-side η-trivialization (Paper 18 `thm:eta_trivialization`) forces the inner-factor period ring to $M_1 \cup M_2$ via the chirality grading axiom $\{\gamma_F, D_F\} = 0$. The chemistry-side analog would be a structural property of the FrozenCore $Z_{\rm eff}(r)$ pipeline mirroring this chirality grading.

The most concrete candidate direction: identify a discrete grading on the FrozenCore radial profile pipeline (e.g., parity of multi-zeta exponents, partition of core orbital basis into chirality-graded blocks) that produces an η-trivialization-style theorem for chemistry inner-factor data.

Sprint-scale path is unclear. Unlike the SM-side case where the axiom $\{\gamma_F, D_F\} = 0$ is established in Krajewski 1998 and canonically adopted in CCM 2010, the chemistry-side analog would need to identify the analog grading axiom from scratch. Multi-month NCG-research target.

**Catalogue status of H6, H7 unchanged:** empirical-only pending the structural-research target's completion. The W1e period-class diagnostic (`debug/sprint_w1e_period_class_memo.md`, 2026-06-04) establishes that H6, H7 are not in $M_1 \cup M_2 \cup M_3$ at audit ceiling 100; the open question is whether they sit in a structurally distinct chemistry-side ring forced by a hypothetical analog axiom.

**Closure type.** Documentation, not theorem. The C-arc's sprint-scale run closes at v3.103.0 with this target named precisely but not attempted. Multi-month engagement is the right next move only if PI commits to a chemistry-side NCG-research arc.

## The C-arc terminal state

**`rem:c_arc_terminal_state` consolidates the full final tally.** Eight theorem-grade non-selection results in the corpus across four sectors:

| # | Theorem | Sector | Sprint | Version |
|:-:|:--------|:-------|:-------|:-------:|
| 1 | Forced-Count | Inner-factor | pre-2026-06-08 | — |
| 2 | H1 Yukawa non-selection | Inner-factor | pre-2026-06-08 | — |
| 3 | $N_{\mathrm{gen}}$ non-selection | Inner-factor | C3 | v3.99.0 |
| 4 | Inner KO-dim non-selection | Inner-factor | F3 | v3.100.0 |
| 5 | Single-cutoff spectral action | Multi-focal renormalization | E7/E8 | v3.101.0 |
| 6 | Spatial-composition radial wall | Multi-focal spatial | G1/G2/G5 | v3.102.0 |
| 7 | No-single-mechanism K | α combination | E6 | v3.103.0 |
| 8 | Cutoff function external | Gravity | D5/D6 | v3.103.0 |

**Sectors closed at theorem-grade:**
- Inner-factor structural-skeleton boundary (canonical-rep level): theorems 1-4
- Multi-focal-composition wall: theorems 5-6 (full closure both sub-sectors)
- α combination rule single-mechanism generation: theorem 7
- Gravity cutoff-function external selection: theorem 8

**Named open structural-research target:**
- Chemistry-side η-trivialization analog (`rem:chemistry_eta_analog`): multi-month NCG-research; H6, H7 remain empirical-only pending.

**Subsidiary residuals (no separate theorem needed):**
- F4 (Higgs VEV): subsidiary calibration scale inheriting from F1 (Yukawa non-selection)
- F5/F6 (CKM, PMNS): subsidiary cross-generation Yukawa-Higgs structure
- F7 (Neutrino masses): subsidiary (inherits F1 + F4)
- I1 (α value): inherits E6 (combination rule conjectural)
- I2 (Born rule): foundational calibration via Gleason via Hilbert inheritance; framework does not improve on Gleason
- I3 (Higgs direction $\hat n$): conditional, depends on Hopf-base identification

## What this closure means

The forced/free seam has gone from "scattered observations across CLAUDE.md §1.7 + §3 + sprint memos" (pre-v3.97.0) to "characterised by eight theorem-grade non-selection results across four sectors with one named multi-month follow-on" (post-v3.103.0) in one day of work.

**Three things this corpus state now supports:**

1. **A clean structural story.** The framework's structural-skeleton scope can now be stated with precision: it forces the inner-factor moduli space + dimension + cross-domain unity in selection rules + universal propinquity rate + master Mellin engine classification + … and it is structurally silent on Yukawa values, $N_{\mathrm{gen}}$, inner KO-dim, multi-loop counterterms, spatial-composition radial profiles, the α combination rule, and the cosmological-constant cutoff moments. The silence is theorem-grade in eight cases.

2. **A clean entry point for outreach (path #2).** Each theorem is a concrete handle for an external mathematician (Brown, Marcolli, Glanois, Kleinschmidt, vS, etc.) to engage with. The paper-grade entry tickets are now complete for the Brown-Kleinschmidt ERC engagement.

3. **A clean named multi-year research target.** The chemistry-side η-trivialization analog is the singular remaining structural-research target. If the project ever commits to a chemistry-side NCG-research arc, this is where it goes.

## Decision gate

The C-arc has reached its terminal state per `rem:c_arc_terminal_state`. PI choice:

1. **Release v3.103.0** and pause. The C-arc is decisively closed; further work on the forced/free seam either tackles the chemistry-analog multi-month target or is residual subsidiaries.
2. **Pivot to path #2 (Brown-Kleinschmidt outreach).** Eight theorem-grade results across four sectors provides a strong entry-ticket profile.
3. **Pivot to a different vein within path #1.** Vein A (Tannakian extension), Vein B (period classification sharpening), or Vein D (Hain-Brown / JLO continuation) all remain.

Honest read: this is a natural stopping point for the project, not just the C-arc. The corpus state is in its strongest single-day-net-gain configuration in the project's history.

## Files

- **`papers/group1_operator_algebras/paper_32_spectral_triple.tex`** — three new theorems (`thm:no_single_mechanism_K`, `thm:cutoff_function_external`) + one corollary (`cor:cc_fine_tuning_external`) + four remarks (`rem:e6_scope`, `rem:d5_d6_scope`, `rem:chemistry_eta_analog`, `rem:c_arc_terminal_state`) added to §VIII. Three-pass clean compile at 83 pages.
- **`debug/sprint_c_arc_closure_e6_d5d6_chemistry_memo.md`** — this canonical bundled memo.
- (Pending) Paper 57 §6.3 update with final tally.
- (Pending) CLAUDE.md + CHANGELOG updates.

## Honest scope

- **Two theorem-grade upgrades** (E6, D5/D6) follow the existing C3/F3/E7-E8/G1-G5 pattern: composition of existing structural content into Theorem blocks. No new mathematics in the proofs.
- **One scoping closure** (chemistry analog): documentation-grade closure of a multi-month NCG-research target. Names the target precisely; does not attempt the theorem.
- **One terminal-state remark** (`rem:c_arc_terminal_state`): tally + sector map of the final C-arc state. Documentation-grade closure of the arc.
- **All conditional on canonical CCM rep** (where applicable). Unconditional versions are multi-year NCG-research targets.
- **Sibling-axiom directions named** for E6 (sibling axiom that generates K from single mechanism), D5/D6 (RG-flow / regulator-independence / asymptotic-safety axiom), chemistry analog (chirality grading on FrozenCore pipeline). All unexhibited.
- **§13.5 protection preserved** — E6 theorem strengthens structural impossibility, does not promote $K$ from observation to theorem-grade derivation. Paper 2 conjectural label intact.

## Cross-references

**Strengthens / consolidates:**
- Paper 2 §IV.G (12-mechanism elimination)
- Paper 18 `thm:eta_trivialization` (SM-side analog)
- Paper 32 `thm:pi_source_case_exhaustion` (case-exhaustion)
- Paper 51 §G7 (Wald's theorem) + §G4-5 (sector-wise Mellin moment map)
- `memory/cc_phi_moment_fine_tuning_statement.md`
- `memory/inner_factor_mellin_engine.md`
- `debug/sprint_w1e_period_class_memo.md` (chemistry-side empirical wall)

**Updates:**
- Paper 32 §VIII (three theorems + four remarks + terminal-state)
- (Pending) Paper 57 §6.3 final tally
- (Pending) CLAUDE.md + CHANGELOG for v3.103.0
