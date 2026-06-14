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

**Compile note:** all edits are text/table only (F uses text refs, no new macros/`\cite`) — low risk; run three-pass compile before any release of these.
