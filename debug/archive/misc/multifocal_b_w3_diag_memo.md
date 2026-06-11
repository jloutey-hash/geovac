# Phase B Sub-Sprint W3-diag — Inner-factor parameter selection wall

**Date:** 2026-05-07
**Track:** Phase B-W3-diag of the multi-focal-composition sprint
**Author:** PM (diagnostic; no production files modified, no paper edits)
**Sources:** Phase A synthesis §6 Q2 + §1 W3 row; Phase A internal audit §D + entry B12; CLAUDE.md §1.7 WH register (WH3, WH4, WH5); CLAUDE.md §2 Bertrand × Hopf-tower entry, multi-focal wall pattern, structural-skeleton scope pattern; Sprint H1 + G3 + G4a memos; `inner_factor_mellin_engine_memo.md`; `sm_gauge_content_forcing_memo.md`; Paper 18 §IV (six tiers); Paper 32 §VIII.C; memory files `bertrand_sm_gauge_truncation`, `geovac_structural_skeleton_scope_pattern`, `inner_factor_mellin_engine`.

---

## §1. The diagnostic question and W3's place in the taxonomy

Phase A's refined wall taxonomy distinguishes three structural classes among the six sub-walls W1a/b/c, W2a, W2b, W3:

- **W1a/b/c** — *cross-register two-body operator* failures. Architecture permits the operator; it has not been built. Tooling-addressable internally.
- **W2a/W2b** — *frontier-of-field* failures. UV/IR renormalization (W2a) and tensor product of two infinite metric spectral triples (W2b) are open in the broader spectral-action / NCG programs, not GeoVac-specific defects.
- **W3** — *inner-factor parameter selection*. Operator class is well-defined and the construction works, but no GeoVac-internal data selects a value within the class.

W3's load-bearing instances: (i) Sprint H1's Yukawa-undetermined verdict (the off-diagonal $\mathbb{C} \leftrightarrow \mathbb{H}$ block of $D_F$ produces a Higgs sector for any non-zero $Y$, but no GeoVac data couples to that block); (ii) G3's commuting-$\mathbb{Z}_2$ theorem ($\|\gamma_{\mathrm{GV}} \otimes I_F - I_{\mathrm{GV}} \otimes \gamma_F\|_{\mathrm{op}} = 2$ exactly); (iii) the Sprint H1 / G3 collapse — G2 ("no autonomous Yukawa") and G3 ("no chirality co-location") are one structural fact; (iv) the Bertrand × Hopf-tower upper bound's *missing lower bound*; (v) Paper 18 §IV.6's inner-factor input data tier (parameter-tied Dirichlet ring $\mathbb{Q}[y_1^{-2s}, \ldots, y_n^{-2s}]$, structurally orthogonal to outer-factor exchange constants).

The PI has named this the *second packing axiom* question: the Paper 0 packing axiom forces the discrete labels (n, l, m, s) and S³ topology, producing the *outer* spectral triple. If a sibling axiom existed that forced calibration content (Yukawa values, magnetization profiles, generation count, hypercharge assignments, multi-focal compositions), the framework would close. Phase B-W3-diag tests whether that sibling axiom has any concrete content.

---

## §2. Q1: Catalogue of every "second packing axiom" speculation

The phrase "second packing axiom" appears explicitly in three locations:

1. **`debug/inner_factor_mellin_engine_memo.md` §7–§8** (Sprint Angle 2, 2026-05-07): "Does not provide a 'second packing axiom.' GeoVac's packing axiom is a 2D→3D lift producing (n,l,m,s). No analogous axiom is here for color or generations." §8 verdict: "The 'second packing axiom' question is genuinely unsolved by this analysis."

2. **Paper 18 §IV at line 1854**: "This is a clean negative on the 'second packing axiom' question: GeoVac's machinery says where $A_F$'s Dirichlet ring sits in the taxonomy but does not pick out which finite spectral triple realizes it."

3. **`debug/multifocal_audit_internal_memo.md` §D + B12 + §E**: B12's "what's missing structurally" is "A separate axiom or structure that supplies inner-factor data ($A_F$ choice, Yukawa values, hypercharges, generation count, renormalization counterterms). The 'second packing axiom' question." §E follow-up: "Whether a 'second packing axiom' candidate has emerged anywhere in the project record... the audit didn't find any concrete proposal."

Beyond the explicit phrase, related framings:

4. **CLAUDE.md §1.7 WH3** ("the lattice exists a priori; match to physics is evidence, not derivation"). Origin-story framing that the existing axiom does more work than retrospective fitting could explain. *Does not propose a second axiom.*

5. **CLAUDE.md §1.7 WH4** ("the four-way S³ coincidence is one structure expressing itself four times"). Strongest coincidence in the project; suggests structural unity but proposes no axiom. (See §5.)

6. **CLAUDE.md §1.7 WH5** ("α is a projection constant, not a derivable number"). K = π(B + F − Δ) accepted as structural coincidence with three identified spectral homes; twelve mechanisms eliminated for K = single-principle derivation. *Implicitly accepts non-derivation; does not propose what would need to be added.*

7. **`debug/sm_gauge_content_forcing_memo.md` §8** (2026-05-06). Three things would upgrade Bertrand × Hopf-tower from "compatible with SM" to "forces SM": (i) generation tripling structure, (ii) closed forcing argument that S⁷ cannot exist in GeoVac, (iii) co-located SM gauge content on a single triple (G4a). *Structural targets, not axioms.*

8. **`debug/geovac_structural_skeleton_scope_pattern.md`** (memory file). Pattern observation across four convergent findings (Bertrand, inner Mellin, H1/G3/G4a, LS-8a). Names the second packing axiom question as "the speculative frontier."

9. **CLAUDE.md §2 Sprint TS-C #14 rest-mass framing**: "rest-mass: CONFIRMED as central deformation $D^2 \to D^2 + m^2 \cdot \mathbf{1}$. New tier in Paper 31 §VII.2.5: parametric / central sector." Closest structural-place-for-parameter, but rest-mass is a Layer-2 input, not a derivation.

10. **The Paper 18 sixth tier itself.** Adding "inner-factor input data" as a categorically disjoint tier is not a *derivation*; it is a *structural place*. The taxonomy distinguishes intrinsic / calibration / embedding / flow / composition / inner-factor-input — six tiers; what populates the sixth is empirical.

11. **`debug/g4a_connes_sm_scoping_memo.md` calibration list**. G4a names eight specific calibration choices that remain free even after Connes' SM construction is set up: number of generations; four Yukawa entries per generation; CKM mixing; PMNS mixing; hypercharge normalization; $A_F$ choice at the level above. *Catalogue of "what GeoVac does not pick" at the Connes-SM level.*

12. **Phase 4H Track SM-B** (April 2026). Tested $\Sigma_f N_c Q_f^2 = 8 = |\lambda_3|$ as a structural map between SM charge content and S³ shell content. NEGATIVE — per-generation 8/3 is charge-universal while every per-shell S³ invariant varies with n. (CLAUDE.md §3.)

13. **Phase 4H Track SM-C** (April 2026). Tested whether $\Delta = 1/40$ embeds naturally in SU(5)/E_8/Spin(10) GUT structures. NEGATIVE — closest hit is tautological renaming using integers 5 and 8 already present in Paper 2. (CLAUDE.md §3.)

14. **Sprint G4a deferred candidate** (Paper 32 §VIII.D). Not yet attempted; predicted positive-thin per Sprint H1 / G3. Closest concrete sprint scope for closing within-spectral-triple electroweak unification, but does not derive Yukawas.

**What the catalogue contains and does not contain.** Many instances of *naming* the question (1, 2, 3, 6, 8, 10), several structural targets that would close gaps the axiom currently fills (4, 7, 11), two specific second-axiom candidates that were tested and falsified (12, 13), one tractable sprint scope that does not derive parameters (14), one mathematical target whose closure would constitute "lower bound forcing" in the gauge sector (item 7's "no S⁷ theorem"). **Zero instances of a concrete proposal for a second packing axiom.** The phrase is consistently used as a placeholder for structure that does not yet exist, not as an active research target.

---

## §3. Q2: W3-distinct vs W3 ⊂ W2a — argued explicitly

### §3.1 Argument for W3-distinct (Phase A's recommendation)

The failure modes are structurally different:

- **W3 (Sprint H1).** Operator class well-defined; the off-diagonal $D_F$ block is a structurally permissible Hermitian operator. Failure is *value selection within the class*. G3's commuting-$\mathbb{Z}_2$ theorem makes this concrete: $\gamma_{\mathrm{GV}}$ and $\gamma_F$ are independent commuting $\mathbb{Z}_2$'s with operator-norm residual exactly 2. Construction works for any $Y \neq 0$; negative is silence on which $Y$.

- **W2a (LS-8a).** Operator class would be "iterated CC spectral action with renormalization counterterms," but the bare iterated CC action *cannot generate the counterterms*. Framework reproduces the divergent integrand correctly (right prefactor, sign, divergence ~$N^{3.43}$) but lacks the subtraction machinery. Failure is *machinery-not-present*, not value-selection.

These are different in mathematical content. W3 is "operator exists, value not selected." W2a is "operator class incomplete; subtraction machinery not in the framework." G3's commuting-$\mathbb{Z}_2$ is *proven*; no analogous proof exists for W2a (only a power-law divergence diagnostic).

### §3.2 Argument for W3 ⊂ W2a

Every multi-loop renormalization counterterm IS an inner-factor input under Paper 18's six-tier sense. The Yukawa $Y$ in the AC inner factor and the QED counterterm $\delta m$ both sit in the parameter-tied Dirichlet ring (Paper 18 §IV.6) — both are numerical inputs that determine specific physical observables but are not derivable from other framework axioms. Under this reading, W2a is W3 restricted to RG counterterms.

The Phase A internal audit §D acknowledges this explicitly: "Whether to call this a third wall (W3) or a sub-wall of W2 is a framing choice."

### §3.3 Resolution

The cleanest reading distinguishes them by *what the framework does at the operator level*:

- **W2a** is *operator-class incompleteness* — the natural CC iteration produces a divergent series; closing W2a requires adding the subtraction operator at the *operator* level.
- **W3** is *parameter-value-selection-within-a-complete-operator-class*. Closing W3 requires data outside the framework's input.

Under this resolution, **W2a counterterms become W3 inputs once they are admitted as parameters in a complete operator class** — i.e., once the renormalized operator class is set up with the counterterms as free parameters, the value-selection question is genuinely a W3 question. The "running" reading is that W2a is the *upstream* wall that creates a W3 problem downstream once it is closed by importing renormalization machinery.

**W3-distinct is the right framing for the current state**; W3 ⊂ W2a is the right framing if/when W2a is closed by adding subtraction machinery. The two readings are not mutually exclusive at different stages of framework development. **Phase A's recommendation (W3-distinct) is confirmed.**

---

## §4. Q3: Concrete proposal hunt

A "concrete proposal" must be expressible as a structural axiom or derivation target whose success/failure is testable. Evaluating each catalogued speculation:

| Item | Concreteness |
|:-----|:-------------|
| 1 (`inner_factor_mellin_engine_memo.md`) | Named open question; explicitly "genuinely unsolved" — **not concrete** |
| 2 (Paper 18 §IV) | Named open question in published-text vehicle — **not concrete** |
| 3 (`multifocal_audit_internal_memo.md`) | Named follow-up — **not concrete** |
| 4 (WH3) | Origin-story framing — **not concrete (different framing)** |
| 5 (WH4) | Coincidence observation — **not concrete (suggests target, see §5)** |
| 6 (WH5) | Compatible with W3 stance — **not concrete** |
| 7 (`sm_gauge_content_forcing_memo.md` §8) | Three structural targets — **concrete-but-targets-not-axioms** (§6) |
| 8 (skeleton pattern memo) | Pattern observation; explicitly "speculative frontier" — **not concrete** |
| 9 (TS-C rest-mass tier) | Structural-place-for-parameter — **not concrete** |
| 10 (Paper 18 sixth tier) | Taxonomy without content — **not concrete** |
| 11 (G4a calibration list) | Catalogue of free parameters — **concrete-as-catalogue, no derivation proposal** |
| 12 (Phase 4H SM-B) | Specific structural map; tested — **concrete but FALSIFIED** |
| 13 (Phase 4H SM-C) | Specific GUT embedding; tested — **concrete but FALSIFIED** |
| 14 (G4a sprint) | Constructs SM architecture given $A_F$; does not derive — **concrete-but-architecture-not-derivation** |

**Net finding: no concrete proposal for a second packing axiom exists in the project record.** Two specific second-axiom candidates (Tracks SM-B and SM-C, Phase 4H) were tested and falsified. No alternative concrete proposal has been formulated since.

The closest "concrete content" is:

(a) **The Bertrand × Hopf-tower lower-bound gap** (item 7): a structural argument fixing the lower bound — i.e., why the SM gauge content saturates the upper limit rather than being a subset — would constitute partial second-axiom content for the gauge sector. Testable in principle (§6).

(b) **The four-way S³ coincidence as derivation target** (item 5, WH4): if S³'s four roles were forced by a single structural statement, that statement would have second-axiom-ish character. (§5.)

(c) **Krajewski–Mellin compatibility audit** (Option B from `inner_factor_mellin_engine_memo.md`, ~2-3 weeks): systematically compute Tr($D_F^k \cdot e^{-tD_F^2}$) for the small Krajewski classification (dim $H_F \le 16$, KO-6 only) and rank candidates by GeoVac-compatibility. If GeoVac-compatibility metrics are picky enough to select a small subclass, that would be partial second-axiom content for inner-factor selection. Flagged as follow-up; not yet executed.

None of (a)/(b)/(c) is a concrete second-axiom *proposal* — they are *targets* whose closure would constitute partial second-axiom content. The status: the second packing axiom question is *not a concrete research target*; it is a placeholder for a class of structural results that, if they existed, would close the catalogued gaps.

---

## §5. Q4: WH4 four-way S³ coincidence — derivation target or speculation?

WH4: the four roles of S³ (Fock projection image, Hopf base, Dirac spin carrier, SU(2) gauge manifold) are *one structure expressing itself four times*. S³ = SU(2) is the unique rank-1 non-abelian compact Lie group where all four roles coincide.

Can WH4 be promoted from coincidence-observation to derivation-target?

**Argument that WH4 is a tractable derivation target:**

1. The Fock projection image (Paper 7) is determined by SO(4) = SU(2) × SU(2) symmetry of the Coulomb Hamiltonian → base is forced to be $S^3 = SO(4)/SO(3)$.
2. The Camporesi–Higuchi spinor bundle (Paper 32 / Sprint TS) is canonically defined on $S^3$ by parallelizability → Dirac carrier role.
3. The Hopf base $S^2$ and circle fiber $S^1$ are determined by the maximal-torus reduction of SU(2) on itself, automatic from $S^3 = SU(2)$.
4. SU(2) gauge (Paper 30) is the strict-natural Wilson construction on $S^3 = SU(2)$.

**Three of the four roles are automatically forced once the first one is established.** The forcing chain: (a) Coulomb $\to S^3$; (b) $S^3 = SU(2)$ parallelizable $\to$ Dirac spinor bundle; (c) $S^3 \to S^2$ Hopf is the maximal-torus reduction; (d) Wilson SU(2) on $S^3 = SU(2)$ is strict-natural.

**The fourfold coincidence is structurally a one-way coincidence (Coulomb gives $S^3$) plus three forced consequences.** Accepting Bertrand's theorem (Coulomb is one of two central potentials with all bound orbits closed), "outer GeoVac is hosted on $S^3$" is forced by the Hamiltonian, and WH4 is just a *reading* of that fact.

**However**, this reading does not extend to inner-factor selection. The forcing chain (a)–(d) determines the *outer* spectral triple structure; it says nothing about $A_F$, Yukawa values, or generation count.

**Verdict on Q4: WH4 is a derivation target for the *outer* spectral triple's structural unity (essentially completed by Sprint TS-D plus the master Mellin engine reading), but is NOT a candidate axiom for inner-factor selection.** The S³ structure is universally Coulomb-side; Bargmann S⁵ (Paper 24) is HO-side; nothing in the existing chain reaches $A_F$. WH4 may deserve deflation to "outer triple's structural unity is forced by Fock projection from $-Z/r$" — see §7 closing recommendation.

---

## §6. Q5: Bertrand × Hopf-tower lower bound — is there a structural argument?

The Bertrand × Hopf-tower argument establishes:

- **Upper bound (forced).** GeoVac produces only $S^3$ and $S^5$ as natural sub-manifolds, by Bertrand's theorem. The complex-Hopf tower $S^{2n-1} \to SU(n)$ truncated at $n \le 3$ is the natural Wilson gauge content. Higher SU(n>3), Sp(n), G_2, exceptionals all require sub-manifolds GeoVac doesn't produce.

- **Lower bound (NOT forced).** Why the SM saturates the upper bound — i.e., why all three of U(1), SU(2), SU(3) are realized rather than a proper subset — is not forced by GeoVac.

A structural argument fixing the lower bound would have to: (i) force HO physics to be present, justifying $S^5$; (ii) force Coulomb physics to be present, justifying $S^3$; (iii) force the U(1) sub-bundle to be relevant alongside SU(2) on $S^3$. Each of these is essentially a second packing axiom in disguise — empirical observations that the universe has both atoms and nuclei, that Coulomb-class bound systems exist, that the abelian sub-construction is independently relevant rather than a derived feature.

**There is no structural argument in the project record that fixes the lower bound on gauge content.** The Bertrand × Hopf-tower memo §8 explicitly names the closest target: a closed forcing argument that S⁷ cannot exist in GeoVac. "Currently the absence of S⁷ is empirical (no mechanism has produced it). A structural theorem ('any GeoVac packing produces only spheres of dim 3 and 5') would upgrade the truncation argument from 'currently no candidate' to 'no candidate exists.' Bertrand's theorem is the natural starting point but doesn't directly give this."

**Verdict on Q5: NO structural argument fixes the lower bound on gauge content.** Closing the lower bound would require a second axiom or equivalent. The closest-tractable target ("any GeoVac packing produces only spheres of dim 3 and 5") is named as a desired theorem but has no proof attempt. This is W3-class content.

---

## §7. Q6: Verdict

The verdict options:

- **(a)** W3 is its own wall, structurally distinct from W2a, with at least one concrete proposal worth a Phase C sprint.
- **(b)** W3 is its own wall but no concrete proposal exists; W3 is a "framework's open question" to be honestly named in Paper 18 §IV.6 or Paper 32 §VIII.D rather than attacked.
- **(c)** W3 ⊂ W2a; the wall taxonomy should fold W3 into W2a.
- **(d)** W3 is the wrong abstraction; the parameter-selection question lives in CalibrationPhysics-not-StructurePhysics, end-of-story.

**My verdict: (b).** W3 is its own wall, structurally distinct from W2a (per §3 — value-selection-within-permissible-class is categorically different from operator-class-incompleteness; G3's commuting-$\mathbb{Z}_2$ theorem is the structural fingerprint), but no concrete proposal exists in the project record. The fourteen catalogued items (§2) are exhaustive of what the project has named; none is a proposal in the strict sense. Two specific second-axiom candidates (Tracks SM-B, SM-C) were tested in Phase 4H and falsified.

W3 is best handled as a *framework's open question* with structural content (the inner-factor input data tier in Paper 18 §IV.6 already does this; Paper 32 §VIII.C's H1 verdict already does this; the Bertrand × Hopf-tower lower-bound gap already does this). Phase C-W3 should be scoped as *cataloguing and framing*, not closure — verdict-conditional:

- If a concrete proposal emerges later (e.g., from a Krajewski–Mellin compatibility audit picking a small subclass), upgrade to verdict (a).
- If W2a is closed by adding renormalization machinery, the W3/W2a boundary will need re-examining (§3.3).
- If verdict (d) becomes the operating stance, W3's status should be reframed as a *scope statement* rather than an open frontier.

The strongest argument against (a) is that no concrete proposal exists despite a year of investigation; against (c) is the structural difference in failure modes (§3.1); against (d) is that calling parameter selection "end-of-story" is premature without the Krajewski–Mellin audit being run. **(b) is the only honest stance.**

---

## §8. Verdict (b) implementation: draft language for Paper 18 §IV.6 / Paper 32 §VIII.D

If the PI accepts verdict (b), the natural venue is a brief addition to **Paper 18 §IV.6** (the inner-factor input data tier) naming the open question explicitly. Draft text (~150 words; not applied, awaiting PI direction):

> **Open question: a second packing axiom for inner-factor data.** The Paper~0 packing axiom forces the outer spectral triple's discrete labels $(n, l, m, s)$ and the $S^3$ topology, but is silent on inner-factor content $A_F$. Sprint H1's $A_{\mathrm{GV}} \otimes (\mathbb{C} \oplus \mathbb{H})$ extension produces a Higgs sector for any non-zero Yukawa~$Y$, but no GeoVac data couples to that selection (the G3 commuting-$\mathbb{Z}_2$ theorem $\| \gamma_{\mathrm{GV}} \otimes I_F - I_{\mathrm{GV}} \otimes \gamma_F \|_{\mathrm{op}} = 2$ is the structural fingerprint). The Bertrand × Hopf-tower truncation provides only an upper bound on gauge content; the lower-bound forcing (why SU(2) and U(1) appear alongside SU(3), why generation count = 3, why specific Yukawa ratios) is calibration-side. Whether a *sibling* axiom — distinct from the Paper~0 packing axiom — generates this data is an open question. Two candidate routes (charge-universal generation map; SU(5)/$E_8$/Spin(10) topological embedding) were tested in Phase~4H and falsified (Tracks SM-B, SM-C). The honest scope statement is that GeoVac's machinery places inner-factor data in this taxonomic tier but does not derive it.

Alternative venue: a one-paragraph addendum to **Paper 32 §VIII.D** (when next touched, after Paper 38 lands). Paper 32 already has the H1 verdict in §VIII.C; a §VIII.D paragraph naming the lower-bound gauge gap and the inner-factor selection gap as *one structural question* (the second packing axiom question) would be appropriate. Both venues are acceptable; Paper 18 §IV.6 is preferred because the six-tier taxonomy is the natural place for taxonomy-level scope statements, and Paper 32 §VIII is currently in equilibrium per the Phase A synthesis Q4 recommendation.

---

## §9. Honest scope and uncertainty

**Catalogue scope.** I read the listed memos and CLAUDE.md sections directly. The catalogue is exhaustive within the search terms (`second packing axiom`, `sibling axiom`, `inner.factor.input`, `calibration.side`) plus the structural context from WH3/4/5 and Bertrand × Hopf-tower. I did *not* search the entire debug/ corpus exhaustively; if a concrete proposal exists under different terminology (e.g., "calibration packing," "color packing"), it would be missed.

**W3 ⊂ W2a framing uncertainty.** §3 argued that W3-distinct is cleaner at the *current state*, with W3 ⊂ W2a being right after W2a closure. The PI may reasonably prefer a more compressed taxonomy (W1/W2 only). The structural difference in failure modes is sharp enough to warrant a separate label, but this is defensible rather than forced.

**Krajewski–Mellin audit payoff uncertainty.** The audit (~2-3 weeks) has not been executed. If GeoVac-compatibility metrics are picky enough, the verdict could shift (b) → (a). If not, (b) is reinforced. I have not predicted the outcome.

**WH4 deflation.** §5 argued the four-way coincidence is structurally a one-way Fock-projection statement plus three forced consequences. This deflates WH4 — it is less a "deep coincidence" than a structural unfolding from Bertrand-then-Fock-then-parallelizability-then-strict-naturalness. The PI may want to re-examine whether WH4 deserves separate WH-register status, or fold it into WH3 or WH1.

**(b) vs (d) strategic choice.** (b) treats W3 as an open frontier (does not foreclose future second-axiom result); (d) treats it as end-of-story (framework explicitly accepts calibration is empirical). (d) is more decisive but also more committal. The choice is strategic-framing rather than technical.

**Follow-ups.** (i) Krajewski–Mellin compatibility audit as the one tractable lever that could shift (b) → (a). (ii) Whether WH4 should be deflated to "outer triple's structural unity is forced by Fock projection from $-Z/r$." (iii) Whether the Bertrand × Hopf-tower lower-bound gap deserves its own diagnostic in Paper 18, since it is the most structurally specific instance of a missing forcing argument in the gauge sector.

---

## §10. Files

This memo: `debug/multifocal_b_w3_diag_memo.md` (~3500 words, 9 sections + verdict implementation).

No production files modified. No paper edits made. No code changes. Existing memos and papers cited verbatim with quotation marks where used.

**Verdict: (b).** W3 is its own wall, structurally distinct from W2a, but no concrete proposal exists in the project record; W3 is best handled as a "framework's open question" via a Paper 18 §IV.6 addition (draft language in §8). Phase C-W3 should be scoped as cataloguing-and-framing rather than closure.

---

**End of Phase B-W3-diag memo.**
