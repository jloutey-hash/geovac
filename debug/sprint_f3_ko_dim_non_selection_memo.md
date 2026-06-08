# Sprint F3 — Inner KO-dim non-selection theorem (2026-06-08)

## TL;DR

**Companion to the C3 N_gen non-selection theorem.** Paper 32 §VIII gains `thm:ko_dim_non_selection` + `rem:full_inner_factor_boundary`. The inner KO-dimensional signature is external calibration data under the canonical CCM SM representation choice. The argument is more direct than for N_gen: packing is kinematic (Paper 0 §VII.B, explicit), KO-dim depends on $(D_F, J_F, \gamma_F)$ data, packing produces none of these. Therefore the framework's axiom set $\mathcal{A}$ cannot autonomously select an inner KO-dim signature; the CCM canonical rep (KO-dim 6) is one choice among consistent alternatives, supplied as external input. Paper 32 compiles three-pass clean at 77 pages.

With this theorem the corpus now carries **four theorem-grade non-selection results** characterising the full free-side content of the inner-factor structural-skeleton boundary at the canonical-representation level: Forced-Count, H1 Yukawa, $N_{\mathrm{gen}}$, inner KO-dim.

## Sprint context

Sprint C3 named F3 (inner KO-dim) as the most tractable next upgrade target after F2 (N_gen): "sprint-scale, structurally analogous to F2." This sprint executes that handoff.

The structural argument is cleaner than for N_gen because the Direction 2 seam-packing scoping memo (`debug/seam_packing_scoping_memo.md`, 2026-06-03) provides the precise composition:

> "KO-dimension is a real-structure/Dirac datum, doubly out of packing's reach. Door 4f T1 / Door 4b Q3: inner KO-dim is a property of $(H_F, D_F, J_F, \gamma_F)$. Packing produces none of $D_F, J_F, \gamma_F$."

Paper 0 §VII.B "What the construction does not provide" supplies the kinematic-scope statement explicitly. Combining the two gives a one-line composition theorem.

## Theorem statement (paraphrase; see Paper 32 `thm:ko_dim_non_selection`)

**Theorem.** Let $\mathcal{A}$ denote the framework's current axiom set. The Paper 0 packing axiom is kinematic, producing labels and graph topology but no real-structure data. The inner KO-dimensional signature is a property of $(D_F, J_F, \gamma_F)$, none of which is in packing's output. The remaining axioms constrain $(D_F, J_F, \gamma_F)$ to be consistent with the outer GeoVac triple's KO-dim 3 and the standard real-spectral-triple axioms, but do not uniquely select an inner KO-dim signature. Different choices yield different signatures, each consistent with $\mathcal{A}$. The canonical CCM SM rep (KO-dim 6) is one such choice, supplied as external input.

**Proof structure.** Two-fact composition:
1. Paper 0 §VII.B explicit scope statement (packing is kinematic).
2. Door 4f T1 / Door 4b Q3 (KO-dim is real-structure data; packing produces no such data).
Composition: $\mathcal{A}$ has no morphism producing $(D_F, J_F, \gamma_F)$ autonomously, and KO-dim depends on this data. Therefore $\mathcal{A}$ does not autonomously select KO-dim.

The conditional is the same as for the N_gen theorem: canonical CCM SM representation. The unconditional question (whether any non-canonical rep that respects SM phenomenology yields a different KO-dim) remains open.

## What this is and is not

**What it is.** A conditional structural non-selection theorem for inner KO-dim under canonical CCM rep, matching the N_gen theorem's structure but with a more direct argument (no representation-distribution analysis needed; the kinematic-scope statement directly closes the question).

**What it isn't.** An unconditional impossibility theorem. The CCM canonical rep IS specified for KO-dim 6 specifically; alternative reps would have alternative KO-dims by definition. Whether any non-canonical rep with KO-dim ≠ 6 admits SM-phenomenologically-consistent fermion content is the open question.

**Companion theorem (rem:full_inner_factor_boundary).** Paper 32 §VIII now has a remark consolidating the four theorem-grade non-selection results: Forced-Count Theorem forces the moduli space and its dimension (128 per generation); H1 forces 8 free Yukawa parameters per generation (the point in the moduli space); $N_{\mathrm{gen}}$ non-selection forces the multiplicity to be external; KO-dim non-selection forces the real-structure signature to be external. The four together characterise the full inner-factor structural-skeleton boundary at the canonical-rep level.

## What's left on the C-arc

Three of the four theorem-grade non-selection results in the corpus were added today (H1 pre-existed; Forced-Count was June 3; F2 N_gen and F3 KO-dim are today). The remaining empirical / diagnostic-grade non-selection entries that could be upgraded:

- **E6** combination rule $K = \pi(B + F - \Delta)$ for $\alpha^{-1}$: 12 mechanisms eliminated empirically; structural impossibility theorem would assert "no single morphism in $\mathcal{A}$ generates $K$ as a combination of the three independent spectral homes." Medium difficulty; would require formalising "single morphism" precisely.
- **E7 / E8** LS-8a multi-loop QED counterterms: structural theorem would assert "the spectral action axiom is single-cutoff; no morphism in $\mathcal{A}$ produces multi-cutoff renormalization counterterms." Medium difficulty.
- **G1–G6** multi-focal composition wall: would unify six independent empirical instances under one theorem ("no multi-focal composition theorem in $\mathcal{A}$"). Hardest of the candidates; would require precise definition of "multi-focal composition theorem" as a structural object.
- **D5 / D6** CC fine-tuning: the $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ requirement is currently an empirical statement; structural-grade upgrade would characterize the cutoff-function moments as external input categorically.
- **Chemistry-side $\eta$-trivialization analog** (H6 / H7): multi-month NCG-research per existing flagging.

## Decision gate

C3 arc has produced two theorem-grade upgrades today (F2 N_gen, F3 KO-dim) in addition to the pre-existing H1 and Forced-Count. Together they close the inner-factor structural-skeleton boundary at the canonical-rep level. The user can:

1. **Pause** the C-arc here. Four theorem-grade non-selection results in the corpus; the inner-factor structural-skeleton boundary is fully theorem-bound at the canonical-rep level. Remaining candidates are outside the inner-factor sector (gravity calibration, QED counterterms, multi-focal composition, chemistry analogs).
2. **Continue** to E7/E8 (LS-8a single-cutoff argument, medium).
3. **Continue** to E6 (combination rule, medium-hard).
4. **Continue** to G1–G6 multi-focal composition unification theorem (hardest).
5. **Pivot** to a different vein within path #1 or to path #2.

PI choice.

## Files

- **`papers/group1_operator_algebras/paper_32_spectral_triple.tex`** — `thm:ko_dim_non_selection` + `rem:full_inner_factor_boundary` added to §VIII, three-pass clean compile at 77 pages.
- **`debug/sprint_f3_ko_dim_non_selection_memo.md`** — this memo.
- (Pending) Paper 57 §3.1 + §6.3 update with F3 reference.
- (Pending) CLAUDE.md + CHANGELOG updates.

## Honest scope

- **Theorem-grade conditional on canonical CCM rep.** Same conditional structure as F2.
- **No new mathematics.** Two-fact composition: Paper 0 §VII.B explicit scope + Door 4f T1 / Door 4b Q3 multiplicity-invisibility result. Both established in 2026-06-02 and 2026-06-03 work. This sprint crystallises the existing composition into a Theorem block.
- **Cleaner structural argument than N_gen.** The N_gen theorem required Read 2's three-obstruction analysis of why the standard CCM rep cannot be reorganised to make algebra factors = generations. The KO-dim argument is more direct: packing is kinematic (Paper 0 §VII.B explicit), KO-dim is real-structure data (definition), packing produces none → packing cannot reach KO-dim. One-line composition.
- **The unconditional version** (whether any non-canonical SM-phenomenologically-consistent rep yields a different KO-dim) remains the multi-year question. The Read 2 sharpest falsifier for N_gen has a natural analog for KO-dim: exhibit a non-canonical rep that admits a different KO-dim while preserving SM fermion-doubling and Majorana neutrino-mass structure.
- **Verified at theorem level by reference, not bit-exactly.** No new computational verification was needed; the proof composition cites existing memos and Paper 0 §VII.B directly.
