# GeoVac Honest Review — Findings Ledger (2026-06-14)

Adversarial review in dependency order (roots first). Severity:
**SOUND** / **FIX** (constructive, non-blocking) / **OVERCLAIM** / **INCONSISTENCY** / **VERIFY** (needs PM confirmation before action).
PM verifies every load-bearing finding + citation before it is treated as actionable.

---

## Paper 7 — Dimensionless Vacuum (TAPROOT hinge) — VERDICT: SOUND + honest

- **SOUND** — §2 "What is New vs Known" attribution exemplary (Fock/Bargmann/Barut known; the 18 proofs *verify* the continuous identities, explicitly "do not claim to derive them for the first time"). §5.2 disavows ontological priority ("a mathematical equivalence, not a proof of physical priority") → §1.5 rhetoric-compliant. §5.7 honestly self-flags convergence as empirical (no rate, to n_max=30).
- **SOUND** — κ=−1/16 is derived not fitted: the universal s-wave Fock coupling c²(n,0)=1/16 from the Chebyshev/Gegenbauer three-term recurrence (amplitude ¼; p₀-independent). Three equivalent readings (§5.4).
- **FIX-1 (forward-ref)** — §5.7 self-flags "no convergence rate proven (empirical to n_max=30)"; Paper 38 (June 2026) later proved rate (4/π)·log n/n. Add forward-ref in §3.3/§5.7 → Paper 38. **VERIFY first** (Paper 38 agent): Paper 7's object = scalar graph-Laplacian spectral convergence (s/p degeneracy recovery); Paper 38's object = operator-system GH/propinquity. Confirm these are the *same* convergence before wording the cross-ref.
- **FIX-2 (framing nit)** — line 111 quotes the weakest κ reading (1/Ω⁴(0), needs Ω(0)=2 ⇒ p₀=1); lead instead with the universal c²(n,0)=1/16 (reading #3).
- **EDGE** — Δ=1/40 in Paper 2's α formula K=π(B+F−Δ) = Paper 7's Fock coupling c²(4,3)=1/40 (§5.4, lines 719–723). Real dependency 2←7; note for Paper 2 adversarial review (the most exotic α ingredient has a derived home in the root).
- **CORRECTION to FIX-1** (from Paper 38 verification): Paper 38 is the **Dirac/spinor** triple; it does NOT rigorize Paper 7's **scalar** Laplacian→Δ_S³ convergence. The scalar state-space GH result is **prior art** (Gaudillot-Estrada–van Suijlekom, arXiv:2310.14733). Corrected cross-ref: "scalar = GE-vS; Paper 38 = the Dirac/spinor analog." (register row 8 confirms.)

---

## Paper 0 — Geometric Packing (AXIOM root) — VERDICT: SOUND (most honest paper on forcing-vs-consistency)

- **SOUND** — packing genuinely FORCES the angular structure (annular area A_k ∝ 2k−1 ⇒ 2ℓ+1 capacity; the identity k²−(k−1)²=2k−1 is construction-internal, not reverse-engineered; WH3 falsifier clears for the angular sector). The forcing-vs-consistency seam is *explicitly drawn by the author* (§IV.C, §VII.C): n-grouping, spin factor 2, radial spacing flagged "natural but not uniquely forced." Rhetoric §1.5-compliant ("a choice of interpretation, not a consequence of the mathematics"). No zombie citations.
- **Load-bearing insight** — Paper 0's NEGATIVE scope (§VII.B "what it does NOT provide": energy, κ, α — "purely kinematic") is MORE load-bearing than its positive claim: Paper 32 §VIII non-selection theorems build on "packing is kinematic ⇒ no real-structure/inner-factor data." [PM-verified: abstract L36–39 defers spectrum to Papers 1/7/18; §VII.B L682–695 disclaims energy/κ/α.]
- **FIX-F (register precision)** — claims-register row 2 attributes the −(n²−1)/κ=−1/16 SYMBOLIC PROOF to "Papers 0, 7", but Paper 0 disclaims deriving them (the test-file proofs are Paper 7's). Retag row 2 "Where" → "Paper 7 (κ + spectrum); Paper 0 (kinematic labels only)". [PM-verified.]
- **FIX (housekeeping)** — Paper 0 dated Feb 2026 but cited as "2025" in paper_34/35 + group3_synthesis. Normalize.

## Paper 32 — Spectral Triple (NCG root, WH1 keystone) — VERDICT: SOUND content, STATUS-DRIFT defects

- **SOUND** — axiom audit cleanly separates finite-cutoff / continuum-limit / imported (Table). All 8 §VIII non-selection theorems conditionalized (𝒜-scope + canonical-CCM-rep), no drift to unconditional. K=π(B+F−Δ) conjectural label preserved at every touch-point (§13.5 honored). Rhetoric clean. NO zombie citations (Krein four-witness explicitly finite-cutoff-only; does not use withdrawn K⁺).
- **FIX-A (INCONSISTENCY, status drift)** — WH1 "MODERATE--STRONG" at L201, L6769 (STALE) vs "WH1 PROVEN" at L1173, 5915, 6026, 6290. Self-contradictory. → "PROVEN (2026-06-10)". [PM-verified via grep.] The §13.11-rule-9 failure mode exactly.
- **FIX-B (OVERCLAIM in thm statement, VERIFY-first)** — `thm:gh_convergence` reportedly stated in "Latrémolière propinquity"; status remark concedes proven object = vS state-space GH (register row 8 correct). Align statement to remark. [PM: grep did not line-confirm exact statement wording — confirm before editing.]
- **FIX-C (propagated)** — Paper 39 citation (~L1375) inherits the same propinquity-vs-state-space gap, undisclosed at point of use. Add one-line caveat.
- **FIX-D (register narrowing)** — case-exhaustion theorem restricted to Paper 34 §III.1–15 (15 of 28 projections); register row 17 implies full catalogue. Narrow row 17.
- **FIX (hygiene)** — H1-theorem cross-refs point to §sec:open (actually §sec:higgs_h1); "§VIII" is shorthand with no actual \section VIII. Low stakes.

## Paper 38 — SU(2) Propinquity Convergence (WH1 PROOF object) — VERDICT: QUALIFIED-SOUND

- **THE SUBSTANTIVE FINDING** — "unconditional" is REAL (not relabeling: G1 genuinely dissolved, G2 genuinely proven as a structural reduction), BUT for a weaker object than the one-word headline implies. Convergence is vS state-space GH under the **translation seminorm**; the paper itself PROVES (rem:dirac_degeneracy) the **truthful Dirac-commutator seminorm does NOT metrize at any finite cutoff** (kernel 10/14 at n=2, 26/55 at n=3). Honest, disclosed, test-frozen — two residual qualifiers behind "unconditional":
  - (i) METRIC: translation seminorm (= Dirac seminorm only in the continuum, never finite-n);
  - (ii) GAP-STILL-OPEN: the kernel condition's base-case non-vanishing (AWA stretched 3-Y element ≠ 0 ∀N) is closed-form + numerics-to-n=5, NOT a general proof (Schur reduces it to one nonzero element/band — small but open).
- **SOUND** — G1 dual reach genuinely DISSOLVED (lifted-state defects collapse to the same Fejér moment γ_n; exact-fit spinor window band-limited by construction). 4/π rate derived (Σ_{d odd}1/d²=π²/8; Hopf-base reading post-hoc, flagged). Prior-art scoping CLEAN (scalar=GE-vS; contribution = spinor transport + per-band injectivity + 4/π). Rhetoric SOUND. No zombie citations. `tests/test_p38_action_seminorm.py` PASSES (3).
- **FIX-H** — abstract should carry the one-clause Dirac-seminorm-degeneracy disclosure (in rem:dirac_degeneracy, not abstract). Minor: "4 digits"→"≈3 digits" at a_1600; cosmetic (n+1)(n+2) vs n(n+1) dimension label in Lemma lifted_state(a) (does not break the conclusion).

---

## CONSOLIDATED — 4 roots (7, 0, 32, 38)

**Content verdict: SOUND.** No physics overclaim, no zombie citations, rhetoric-compliant, K=π(B+F−Δ) conjectural honored, §VIII theorems conditionalized. The foundations are honest.

**Defect pattern: STATUS-DRIFT + cross-ref hygiene + 2 honest-scope sharpenings — NOT soundness.** The content is honest; the *metadata* (status labels, register attributions, cross-refs) drifted under correct-in-place editing. This validates the new Current-State Check + §13.11-rule-9 process additions — they target exactly this failure mode.

**One substantive load-bearing finding:** WH1's "PROVEN unconditional" rests on two small residual qualifiers (translation-seminorm metric; base-case non-vanishing to n=5). Sound, but the one-word "unconditional" is slightly stronger than the proof's own two footnotes. The claims register row 8 already carries the translation-seminorm qualifier; the WH1 §1.7 headline and the "unconditional" rhetoric should match it.

### Actionable update batch (proposed; NOT yet applied)

| ID | Target | Change | Severity | Verify-first? |
|----|--------|--------|----------|---------------|
| A | Paper 32 (L201, L6769) | WH1 "MODERATE-STRONG" → "PROVEN (2026-06-10)" | INCONSISTENCY | no (PM-verified) |
| B | Paper 32 `thm:gh_convergence` | "Latrémolière propinquity" → "vS state-space GH" | OVERCLAIM | **YES** (line-confirm) |
| C | Paper 32 ~L1375 | vS-distance caveat on Paper 39 closure | INCONSISTENCY | YES |
| D | register row 17 | narrow to "15 of 28 projections; 16–28 pending" | precision | no |
| E | register row 2 | "Where" → Paper 7 (κ+spectrum); Paper 0 (labels) | precision | no (PM-verified) |
| F | Paper 7 §3.3/§5.7 | convergence cross-ref: scalar=GE-vS, Dirac analog=Paper 38 | FIX | no (register-confirmed) |
| G | Paper 7 L111 | lead κ with universal c²(n,0)=1/16 | nit | no |
| H | Paper 38 abstract | Dirac-seminorm-degeneracy disclosure; "4 digits"→"≈3" | FIX | no |
| I | WH1 §1.7 + register row 8 | note 2 residual qualifiers behind "unconditional" | **substantive** | PI direction |
| J | Paper 0 + 4 citing files | normalize date → 2026 | housekeeping | no |

### Applied 2026-06-14 (safe batch, this session)

- **A ✓** Paper 32 L201/L6769: WH1 `MODERATE--STRONG` → `PROVEN` (internal self-contradiction resolved).
- **D ✓** claims_register row 17: narrowed to "first 15 of 28 projections (§III.1–15; 16–28 pending)".
- **E ✓** claims_register row 2: attribution → "Paper 7 (κ + spectrum); Paper 0 (kinematic node labels only)".
- **F ✓** Paper 7 §5.7: convergence cross-ref added (scalar = GE-vS arXiv:2310.14733; Dirac analog = Paper 38). Text references, no `\cite` → no undefined-ref risk.
- **G ✓** Paper 7 §VI.C κ line: now leads with the universal s-wave coupling c²(n,0)=1/16 (p₀-independent).
- **J ✓** paper_34 (L9610) + paper_35 (L1398) bibitems: Paper 0 "(2025)" → "(2026)" (authoritative \date = Feb 22 2026; group3_synthesis already 2026; archive files left untouched).
- **H (abstract clause): NO EDIT NEEDED** — already present in Paper 38 abstract L87–88 ("the truthful Dirac-commutator seminorm itself degenerates at finite cutoff; a remark quantifies this") + prior-art scoping L117–122. Agent finding #4 overstated; verified.

### Verify-batch — RESOLVED 2026-06-14

- **I — RESOLVED (no change; agent overstated).** `lem:band_injectivity` proof (L351–371): SO(4)-equivariance + irreducibility ⇒ Schur dichotomy (zero-or-injective, *general*) + non-vanishing of the **stretched** AWA 3-Y element (standardly nonzero = positive product of factorials, *general*), with n≤5 as supplementary confirmation (also checks the stronger full-N²-rank). NOT numerics-only-to-n=5. **WH1 "unconditional" stands; no footnote.** Optional polish (non-blocking): add a one-line in-text "stretched ⇒ positive ⇒ nonzero for all N" so the n≤5 verification isn't misread as the sole support.
- **B ✓ APPLIED** — Paper 32 `thm:gh_convergence` (L3098): "Latrémolière propinquity" → "van Suijlekom state-space GH distance" + Remark~\ref{rem:gh_status_caveat} pointer; displayed bound symbol Λ → d_GH (L3100). Aligns the statement to the paper's own `rem:gh_status_caveat`(iii) + Paper 38 source + claims-register row 8.
- **C ✓ APPLIED** — Paper 32 L1393 (Paper 39 tensor-product citation): same "Latrémolière propinquity" → "vS state-space GH distance" + Remark pointer.
- **H-nit — STILL OPEN (lone remaining item).** Paper 38 L1599 "converging to 4 digits at n=1600": needs the closed-form γ_n extracted + a_1600 recomputed before substituting a digit count. Low value (deep-appendix numerical-verification sentence; the abstract is already honest). Optional.

**Roots review COMPLETE (7, 0, 32, 38): content sound; all status-drift/precision fixes applied except the H-nit. Compile-before-release: B touched a keystone theorem statement — run three-pass compile on Paper 32/7/34/35 before any release.**

### Compile verification (2026-06-14, v4.12.1 checkpoint)

- **Paper 32 — three-pass GATE: PASS** (0 errors / 0 undef / 0 multidef). The keystone `thm:gh_convergence` edit + new `\ref{rem:gh_status_caveat}` are clean.
- **Papers 7 / 34 / 35 — pre-existing undefined refs (NOT from this session's edits):** Paper 7 `\cite{loutey_paper2}` (L721, untouched by me — Paper 7's bib lacks the bibitem); Paper 34 ~5 undefined `\ref`s (`sec:matches`, `sec:curvature_coefficients`, `sec:layer2_d`, `tab:catalog_off`, `sec:proj_pk`, pp.53–66); Paper 35 one. A year change in a bibitem + text additions cannot create undefined refs.
- **LATENT FINDING (new):** the corpus carries pre-existing undefined `\cite`/`\ref` in at least Papers 7, 34, 35. Candidate for a dedicated bib/cross-ref hygiene sweep (out of scope for this review's substantive root pass). Not introduced here; pre-existing.

---

## Foundations BRANCH review (2026-06-14) — synthesis update + adversarial pass

Two parallel opus agents (update-proposal + adversarial). Both verdicts: synthesis + branch papers content **SOUND**; defect = same metadata-drift as the roots. The synthesis is **"one branch short"** — it predates the 54–57 reconvergence.

### VERIFIED CONFLICT — Coulomb/HO layer count = SIX (not four)

Agents disagreed (update agent: four, per Paper 57; adversarial: six, per Paper 24 body). **PM-verified Paper 24 L612–649:** body enumerates **Layers 1–6** (L6 = master Mellin chirality-parity, marked "(new)", June 2026); but Paper 24's own intro prose says "fourth layer" (L612) and a footnote says "five" — Paper 24 is internally drifted (4→5→6 as the count grew). **Canonical current count = 6.** The update agent's "four" was stale; do NOT propagate it.

### Synthesis update plan (Agent 1 draft, PM-reconciled)

- **NEW `\section{...}\label{sec:reconvergence}`** (between §8 `sec:partition` and §9 `sec:open`): integrates Papers 54 (two-body), 55 (periods), 56 (Tannakian/cosmic-Galois injection), 57 (forced/free) as the arc's *reconvergence*. Full draft prose returned by Agent 1 (ready to place).
- **Abstract + §1 intro:** "nine papers" → eleven (0,1,7,18,22,24,31,54,55,56,57); Papers 6/21 = archive lineage; add the taproot→fan-out→reconvergence framing.
- **§3 convergence:** empirical → theorem (Paper 38 Dirac vS state-space GH, rate 4/π; scalar = GE-vS prior art — keep DISTINCT per ledger FIX-F) + Paper 40 universality.
- **§6 (`sec:bs`):** layer count 4↔5 internal contradiction → **SIX** (verified Paper 24 body).
- **§7:** add the Paper 55 period-ring home (M1⊂ℚ[π,π⁻¹]; M2⊂⊕π^{2k}ℚ; M3⊂MT(ℤ[i,1/2]) level≤4); case-exhaustion = "15 of 28" (register row 17).
- **§9 open:** reframe W3 under Paper 57's forced/free seam; add cosmic-Galois equality (multi-year, **NOT claimed**) + M2/M3 unification.
- **Bib + date:** add bibitems for 54–57 + Paper 38/40; bump `\date` to 2026-06-14.

### Adversarial branch findings (Agent 2) — source-paper fixes (follow-on hygiene pass)

- **F1 (INCONSISTENCY, VERIFY canonical convention):** D_pd vs physical-D — Paper 31/synthesis/register row 3 quote the pair-diagonal D_pd (1.44% @ l_max=3) as "the" angular density; Paper 22's physical Coulomb D = 6.06%. Decide canonical convention across Paper 22/31/synthesis/register.
- **F2 (INCONSISTENCY):** layer-count drift across Paper 24 (intro "fourth"/footnote "five"/body "six"), Paper 31 (5 claimed, 3 described), synthesis (4-vs-5). Propagate **SIX**.
- **F3/F4 (FIX, cross-ref hygiene):** undefined `\ref` in Paper 31 (`sec:universal_sector`, `sec:potential_specific_sector` → `sec:universal`/`sec:specific`) + Paper 22 (`sec:theorem`, `sec:universal_partition` + cross-paper Paper-14 refs). Adds Papers 22/31 to the bib-hygiene sweep (with 7/34/35).
- **F5 (FIX, register):** Paper 18 states case-exhaustion unrestricted; add "15 of 28" caveat (register row 17).
- **F6 (OVERCLAIM, borderline):** Paper 55 L1372 "GeoVac IS the abelianisation" too bald → soften to "realizes the abelianised image (depth-2 cocycle cocommutative; equality multi-year)".
- **F7 (INCONSISTENCY, register):** row 18 residual count stale 2,960 → Paper 56 abstract 5,864. Update register.
- **F8 (FIX):** "Latrémolière propinquity" in Paper 54 L67 + Paper 31 B1 → vS state-space GH (same as B/C; the §signature_partition Lorentzian-propinquity uses are legitimately different — do NOT touch those).
- **F9/F10 (minor):** Paper 56 N double-def (`Nsec`=n(n+1)/2 vs `N`=n(n+3)/2); abstract 𝒰₄ vs body 𝒰₄^ab → add `^ab`.
- **Clean:** zero zombie citations; α combination rule conjectural honored everywhere; cosmic-Galois scope correct (injection into 𝒰₄^ab, equality NOT claimed).

### Synthesis update APPLIED 2026-06-14 (v4.12.2) — three-pass GATE: PASS

- NEW `\section{...}\label{sec:reconvergence}` (Papers 54–57); abstract reframed (taproot→fan-out→reconvergence); §1 intro → eleven papers + dependency framing; §3 convergence → theorem (Paper 38 Dirac / GE-vS scalar, kept distinct); §6 layer count → **SIX** (added gravity L5 + master-Mellin L6; fixed "Paper 25 synthesis group" → QED/gauge group); `\date` → 2026-06-14.
- **OPTIONAL-REMAINING (synthesis):** §9 open-questions reframe (W3 under Paper 57; the cosmic-Galois equality frontier is already stated in the reconvergence section — redundant) and §7 period-ring forward pointer (now covered by the reconvergence Periods subsection). Low value; deferred.
- **NEXT BATCH — branch-paper hygiene (Agent 2 F1–F10, source papers):** D_pd convention (22/31/register row 3); layer-drift inside Paper 24 (intro "fourth"/footnote "five" vs body "six") + Paper 31; undefined `\ref` in 22/31; Paper 18 "15 of 28" caveat; Paper 55 L1372 "IS the abelianisation" softening; register row 18 residual 2,960→5,864; Paper 56 N double-def + abstract 𝒰₄→𝒰₄^ab.

**Compile note:** all edits are text/table only (F uses text refs, no new macros/`\cite`) — low risk; run three-pass compile before any release of these.

---

## Trunk code review + QA execution (2026-06-14, v4.13.0)

5-agent per-paper code review (Papers 0,1,7,32,38) + 2-agent test-writing pass. Full claim→test→code matrix: `docs/claim_test_matrix.md`. PM-verified κ and 4/π against primary text.

**QA infrastructure created:** CLAUDE.md §9 Branch QA Review Protocol; `.claude/agents/code-reviewer.md` + `citation-reviewer.md`; `docs/authoring_conventions.md`; `docs/claim_test_matrix.md`.

**Dispositions APPLIED:**
- **κ "derived" → OBSERVATION** [PI-approved]: matching κ=−1/16 *coincides* with geometric 1/16 (c²(n,0), 1/Ω⁴(0)); counterfactual (`test_trunk_qa_kappa.py`) shows no bridge. Reframed Paper 7 (L111+§VII.D), Paper 0 §VI, register row 2 (SYMBOLIC PROOF→OBSERVATION).
- **4/π UPGRADED** (Paper 38): circular `test_asymptotic_constant_value` superseded by `test_trunk_qa_fejer_4_over_pi.py` (derives from γ_n sum rule, rejects 2/π decoy).
- **Paper 1 §III downgraded:** degree table / 16% waypoint / 1.77 gearing not reproduced (CG construction absent); QA-correction note; qualitative stands.
- **Paper 32 Forced-Count downgraded:** full-axiom moduli = 260 not 128 (128=matter-sector); "45 tests verify" re-scoped; correction note.
- **5 BACKED artifacts accepted:** annular→2ℓ+1, c²(4,3)=1/40=Δ, [T₊,L₊]=0, splitting-convergence, 4/π.

**Key structural finding:** L=D−A is positive-semidefinite — −(n²−1) is a *continuum* (S³) property, not the graph's. Machinery (axiom audit + negative controls, O(V), continuum proofs, sum rules, K-conjectural) genuinely backed; soft spot = headline derived/exact constants.

**Compile:** Papers 0, 32 GATE PASS; Papers 1, 7 fail only on PRE-EXISTING undefs (`Note1` L76; `loutey_paper2`). All 8 `test_trunk_qa_*` green.

---

## Trunk QA RE-CHECK + K→observation (2026-06-14, v4.14.0)

Full chronicle: CHANGELOG [4.14.0]. Summary only here.

**Fresh re-check (§9 principle 2):** 5 fresh opus agents (one per trunk paper) re-derived pass-1 verdicts independently. *Every* downgrade held under independent reproduction (Paper 1 agent recomputed 1.9/3.8, 1.7%@n10, ratio 2.0 from scratch). Both Paper 32 hard-prohibition flags CLEAR (K-conjectural intact pre-downgrade; WH1 vS-GH not over-stated).

**Upgrades (two-way verdict, principle 3):** per-band injectivity (P38) BACKED-WEAK→BACKED-SOUND (Schur+AWA all-N); c²(4,3)=1/40 (P7)→INTERNAL THEOREM (two routes). 4/π refined to derived-**numerics-pinned** (doubling estimator, not symbolic E–M) — corrects my v4.13.0 "PROVEN" over-statement.

**New finds → fixed:** Paper 0 |V|=2n_max² BUG → Σ2n² (code was right; pass-1's "spinless" dispute was wrong — reconciled vs primary text + Paper 41). 2 NO-TEST→BACKED tests written (`test_paper1_rydberg`: N=−2[T₊,T₋]→{n}, L²→l(l+1); `test_paper7_vee_s3`: §V exact integrals). Paper 1 §III QA-note under-scoping + L143 contradiction fixed. Paper 32 stale Q2 (J²=+1 → resolved J²=−1), 2 zombie propinquity titles, vacuous `test_J_preserves_O` (print→assert). Paper 38 circular test annotated.

**K = π(B+F−Δ) conjecture→observation** [PI direction]: 14 corpus sites (Paper 2 ×8, Paper 32 ×6) + register + matrix + CLAUDE.md (WH5/§2/§6/§13.5/§13.8). §13.5/§13.8 prohibition reframed: "never stronger than an Observation" (floor=observation). Three derived spectral homes unaffected.

**3 QA principles** memorialized: §9 + code-reviewer/citation-reviewer agents + authoring_conventions rule 9.

**Verification:** Papers 0,2,32 GATE PASS; 1,7 ERRORS=0 (same pre-existing undefs). 214 passed / 4 skipped (new + edited tests + 18 topological proofs).

---

## Trunk citation-grounding pass (2026-06-14, v4.14.1)

Full chronicle: CHANGELOG [4.14.1]. 5 opus web agents (one/paper); PM re-verified every external WRONG-ID vs primary source (AIP/ScienceDirect/EMS/arXiv) before editing.

**Headline: NO fabricated arXiv IDs / nonexistent theorem numbers on the trunk.** The worst historical failure mode (Fursaev–Solodukhin class) is absent. Every hit = metadata slip on a *real* work + one missing internal bibitem.

**Fixes (8):**
- P38 `avery_wen_avery1986` (LOAD-BEARING): 3 authors→2, title, vol 27→26, year 1986→1985, 402→403 → Z.-Y. Wen & J. Avery, JMP 26, 396–403 (1985). In-text "Avery–Wen–Avery"→"Wen–Avery". Proof verifies 3-Y non-vanishing computationally (n≤5) so cite-correction is proof-neutral.
- P38 `latremoliere2016`: Banach JMA 10 (2016) → JMPA 103, 303–351 (2015).
- P32 `glanois2015`: JNT 182, 36–90 (2018) → 160, 334–384 (2016).
- P32 `latremoliere2021dual`: title +"Gromov–Hausdorff".
- P7 `loutey_paper2`: missing bibitem added (clears undef).
- P1: SU(1,1) attribution narrowed off Biedenharn–Louck (SU(2)) → Barut; 2 orphans (`berry1984`,`chung1997`) grounded by citing.
- P0: CLEAN (5/5 grounded).

**Bonus:** the two long-standing trunk undefs cleared — P7 `loutey_paper2` (bibitem) + P1 `Note1` (bibtex cycle). All 4 edited papers now compile GATE PASS / 0 undefs.

**Residual (low-risk):** P32 `connes_vs2021` theorem-numbers (Def 2.39/Prop 4.2/4.3) consistent w/ abstract, not byte-verified vs PDF. **Flagged for PI (separate):** Paper 2 title "…from Spectral Geometry…" reads derivation-y post K→observation — a retitle question, not a cite fix.
