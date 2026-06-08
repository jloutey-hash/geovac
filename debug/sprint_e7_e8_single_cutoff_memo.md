# Sprint E7/E8 — Single-cutoff spectral action theorem (2026-06-08)

## TL;DR

**Structural theorem closing the LS-8a renormalization wall.** Paper 32 §VIII gains `def:multi_cutoff`, `thm:single_cutoff_spectral_action`, `cor:multi_loop_renormalization_wall`, and `rem:single_cutoff_scope`. The Chamseddine--Connes spectral action axiom in $\mathcal{A}$ specifies a single-cutoff functional $S(D, \Lambda) = \mathrm{Tr}\, f(D/\Lambda)$; no morphism derivable from $\mathcal{A}$ produces a multi-cutoff structure. Consequently the framework cannot autonomously generate multi-loop QED renormalization counterterms ($Z_2, \delta m$, etc.), formalising the LS-8a wall observed empirically in Sprint LS-8a (May 2026).

Closes **four catalogue entries** (E7, E8, G3, G4) under one structural theorem. Paper 32 compiles three-pass clean at 79 pages.

## Sprint context

Sprint F3 (this same close-of-day, v3.100.0) named E7/E8 as the next tractable target after the four canonical-rep-level inner-factor non-selection theorems (Forced-Count, H1 Yukawa, $N_{\mathrm{gen}}$, KO-dim). E7/E8 is structurally distinct from those four: it's in the *spectral-action* sector rather than the *inner-factor* sector. The wall is the LS-8a multi-loop QED renormalization gap observed empirically in May 2026 (Paper 36 LS-8a).

The structural reason for the LS-8a wall was already documented in the multi-focal-composition wall pattern memo: "framework reproduces the bare integrand but cannot generate counterterms because the spectral action axiom is single-cutoff." This sprint formalises the structural reason as a theorem in Paper 32 §VIII.

## Theorem statement (paraphrase; see Paper 32 `thm:single_cutoff_spectral_action`)

**Definition (Multi-cutoff structure).** A multi-cutoff spectral action on a real spectral triple is a functional depending on $k \ge 2$ independent energy scales in a manner not reducible to a function of a single scale.

**Theorem.** The framework's current axiom set $\mathcal{A}$ specifies the CC spectral action $S(D, \Lambda) = \mathrm{Tr}\, f(D/\Lambda)$ as a single-cutoff functional. No morphism derivable from $\mathcal{A}$ produces a multi-cutoff spectral action structure.

**Corollary (LS-8a wall).** Multi-loop QED renormalization counterterms ($Z_2, \delta m$, etc.) require a multi-cutoff $(\Lambda_{UV}, \mu_R)$ structure. By the theorem, $\mathcal{A}$ contains no morphism producing such a structure. Therefore multi-loop QED counterterms are external calibration data, classifying catalogue entries E7, E8, G3, G4 as multi-cutoff renormalization calibration data.

**Proof structure.** Four-bullet structural inspection of each axiom in $\mathcal{A}$:
1. CCM spectral action is single-cutoff by stipulation.
2. Paper 0 packing + standard real-spectral-triple axioms + Hopf-tower + Bertrand + Upgrade B are algebraic, kinematic, topological — introduce no energy scale.
3. Inner fluctuations preserve single-Dirac structure.
4. Tensor products preserve single-cutoff structure with combined $D$.

Composition: no axiom in $\mathcal{A}$ introduces an independent second scale.

## What this is and is not

**What it is.** A theorem-grade structural impossibility result for multi-loop QED renormalization within the framework's current axiom set. Closes four catalogue entries (E7, E8, G3, G4) under one structural argument.

**What it isn't.** A claim that multi-cutoff structure is impossible in any extension of the framework. A sibling axiom adding renormalization-group flow machinery would extend the axiom set and could produce multi-cutoff observables. Such a sibling axiom is logically possible and named as the analog of the $N_{\mathrm{gen}}$ sibling-axiom direction.

**Scope clarification: this theorem covers only the RENORMALIZATION sub-sector of the multi-focal-composition wall.** Catalogue entries G1 (HF-3 recoil cross-register $V_{eN}(\hat r_e, \hat R_N)$), G2 (HF-4 Zemach magnetization density), and G5 (W1e chemistry inner-region correlation) are spatial $\times$ spatial composition walls, structurally distinct from the multi-cutoff renormalization sub-sector. They lie outside this theorem's scope. A general multi-focal-composition wall theorem unifying both sub-sectors remains the named open structural-research target (Paper 57 §6.4).

## What's left in the C-arc multi-focal sector

After Sprint E7/E8:

| Entry | Coverage | Status |
|:------|:---------|:-------|
| E7 (Z_2, δm counterterms) | this theorem | THEOREM-grade |
| E8 (multi-loop QED beyond LS-7) | this theorem | THEOREM-grade |
| G3 (HF-5 multi-loop hyperfine) | this theorem | THEOREM-grade |
| G4 (LS-8a, = E7) | this theorem | THEOREM-grade |
| G1 (HF-3 recoil V_eN) | spatial-composition wall, separate | EMPIRICAL |
| G2 (HF-4 Zemach magnetization) | spatial-composition wall, separate | EMPIRICAL |
| G5 (W1e chemistry) | spatial-composition / inner-region wall | EMPIRICAL |
| G6 (= F1 Yukawa) | covered by H1 | THEOREM-grade |
| E6 (combination rule K) | special; 12-mechanism-elimination | EMPIRICAL |

Remaining empirical-only entries in the multi-focal-composition wall family that could be upgraded:
- **G1 / G2 spatial-composition walls** — the framework's Hilbert-space tensor product gives single-particle one-body operators on each register but no native two-body operator. Paper 54 (drafted 2026-06-01) addressed this: tensor-product spectral action gives correct angular selection rules but wrong radial weights for two-body Coulomb (Pearson 0.58 / 0.41 decreasing with $n_{\max}$; the Fock-projection conformal-factor wall). A structural theorem might formalise "the framework's natural composition product does not reproduce two-body Coulomb because the Fock-projection conformal factor breaks the natural composition." Sprint-scale, harder than E7/E8.
- **G5 chemistry-side W1e** — chemistry-side analog of Yukawa non-selection / multi-focal walls. The chemistry-side $\eta$-trivialization analog (named in W1e period-class memo) would be the natural theorem-grade upgrade. Multi-month NCG-research.

## Decision gate

After Sprints C3 + F3 + E7/E8 (this sprint), the corpus carries **five theorem-grade non-selection results** total:
1. Forced-Count Theorem (inner factor moduli dimension)
2. H1 Yukawa non-selection (inner factor diagonal slice)
3. $N_{\mathrm{gen}}$ non-selection (inner factor multiplicity)
4. Inner KO-dim non-selection (inner factor real-structure signature)
5. Single-cutoff spectral action (multi-loop QED renormalization sub-sector)

The first four close the inner-factor structural-skeleton boundary at the canonical-rep level. The fifth closes the renormalization sub-sector of the multi-focal-composition wall.

PI choice:
1. **Continue** to E6 (combination rule $K$ structural impossibility, medium-hard).
2. **Continue** to G1/G2 spatial-composition theorem (medium, uses Paper 54 substrate).
3. **Continue** to D5/D6 CC fine-tuning structural impossibility (medium).
4. **Pause** the C-arc here. Five theorem-grade results in one day; the C-arc is decisively closed for the major sectors. Remaining targets are smaller residuals.
5. **Pivot** to a different vein or path #2 outreach.

## Files

- **`papers/group1_operator_algebras/paper_32_spectral_triple.tex`** — `def:multi_cutoff` + `thm:single_cutoff_spectral_action` + `cor:multi_loop_renormalization_wall` + `rem:single_cutoff_scope` added to §VIII; `\bibitem{paper57}` added. Three-pass clean compile at 79 pages.
- **`debug/sprint_e7_e8_single_cutoff_memo.md`** — this memo.
- (Pending) Paper 57 §6.3 update with E7/E8 closure note.
- (Pending) CLAUDE.md + CHANGELOG updates.

## Honest scope

- **Theorem-grade within the framework's current axiom set $\mathcal{A}$.** The proof inspects each axiom and shows none introduces an independent second scale. Standard structural-inspection argument.
- **No new mathematics.** The structural reason for the LS-8a wall was already named in the multi-focal-wall-pattern memo (single-cutoff spectral action). This sprint formalises the structural reason as a theorem block in Paper 32 §VIII.
- **Theorem covers only the renormalization sub-sector.** Spatial-composition walls (G1, G2, G5) lie outside its scope. A general multi-focal wall theorem unifying both sub-sectors remains open.
- **Sibling-axiom direction.** The theorem does not rule out an extension of $\mathcal{A}$ with multi-cutoff machinery. This is logically possible but requires new primitive content (RG-flow axiom).
- **Cross-paper compatibility.** Paper 36 (Bound-State QED on Dirac-S^3) §LS-8a empirically established the renormalization wall in May 2026. This theorem formalises the structural reason behind that empirical wall. No conflict; the theorem strengthens the existing scope statement.

## Cross-references

**Strengthens / consolidates:**
- Paper 36 §LS-8a (empirical wall, May 2026)
- Paper 54 (drafted 2026-06-01; two-body Coulomb via tensor-product spectral action — addresses spatial-composition sub-sector, outside this theorem's scope)
- `memory/multi_focal_wall_pattern.md`
- `memory/ls8a_two_loop_renormalization_gap.md`

**Updates:**
- Paper 32 §VIII (new theorem + corollary + remark)
- (Pending) Paper 57 §6.3 catalogue closure update
- (Pending) CLAUDE.md + CHANGELOG
