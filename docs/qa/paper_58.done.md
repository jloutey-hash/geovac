# Paper 58 (Abelian Residue) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-58 scope + deltas + watch-notes.

> **STATUS: CERTIFIED ✅ 2026-08-15 — FULL certifying run PASS.** First single-paper
> `/qa` target. Certified by a **full calibrated run over the final combined text**
> (Paper 58 + the group2 synthesis's Paper-58 content): the panel was **fully calibrated**
> (sensitivity **6/6** planted seeds caught — 2 code, 2 citation, 1 claims, 1 synthesis;
> specificity clean, no known-good control false-flagged), with **zero verified material
> defects** beyond the seeds across all four dimensions (code / claims / citations /
> synthesis) plus the completeness-critic, on an all-green deterministic layer (C10–C18).
> This **discharges the 2026-08-14 honest cap**: the earlier cert was on a clean delta with
> the full run PI-waived; the full run has now been run and PASSED, so no single-FULL-run
> gap remains. Only NIT-tier residue survives (see change log). Paper 58 (Aug 10) postdates
> the 2026-06-28 group2 cert; this is its own cert record.

**Scope (single-paper):**
- **Paper 58** — `papers/group2_quantum_chemistry/paper_58_abelian_residue.tex`
  ("Angular Sparsity Is an Atomic-Sector Property: The Abelian Residue of
  Multi-Center Bonding").
- **C9 target** = the **Paper-58 promotions in the group2 synthesis**
  (`papers/synthesis/group2_quantum_chemistry_synthesis.tex`, the 2026-08-14 edits:
  the new "abelian residue" subsection + the W1e/NaH corrections + the
  open-questions/abstract/bib updates). Only the changed regions are the C9 scope;
  the rest of the synthesis rests on its 2026-06-28 cert.
- **Out of scope:** the other 9 group2 papers (unchanged since the 2026-06-28 cert).
  Trunk papers (0/1/7/14/18) canonical; in scope only where Paper 58 restates them (C7).

**Deterministic `--gate`:** `group2` (Paper 58 lives under `group2_quantum_chemistry/`).

## Dimensions exercised (ALL, one invocation — FULL run, first cert of a fresh target)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1 on Paper 58. Map + **RUN**:
  `tests/test_two_center_eri_aabb.py` (closed forms + the `(AA|BB)` decided census);
  the NaH binding panel; the bra/ket certificate `eq:braket_certificate` + its
  two-corruption diagnostic; the Gaussian-corroboration cross-checks. For each: does
  the test genuinely *prove* the claim (not tautological / weaker than prose)?
- **Paper claims / prose (C3, C5, C6, C8)** — `claims-reviewer` ×1 on Paper 58,
  enumeration-forced (every table row, every tier-label, every headline).
- **External citations (C4)** — `citation-reviewer` ×1 on Paper 58 (Mulliken/Roothaan
  `A_n`/`B_n`, Ruedenberg exchange, McMurchie–Davidson, Shibuya–Wulfman, and the
  internal Papers 8–9/11/14/17/19/20 cross-refs).
- **Synthesis faithfulness (C9)** — `claims-reviewer` ×1 on the group2-synthesis
  Paper-58 promotions: do they faithfully reflect Paper 58's tiers, and introduce no
  claim Paper 58 does not support?
- **Deterministic (C10–C18)** — the step-1 scripts, `--gate group2`.
- **Completeness-critic** ×1 (FULL run).

## Branch-defining criterion (inherited from `group2.done.md`): benchmarking + guardrail-negative honesty

The reviewers must verify ALL of: (1) **Benchmarking rule** — every accuracy headline
names its baseline and uses the strongest relevant one; a favorable number vs a weak
baseline without strong-baseline context is MATERIAL. (2) **Guardrail negatives stay
negative** — the Papers 8–9 Sturmian Structural Theorem (no shared-p₀ Sturmian binds a
heteronuclear pair) and the FCI-M graph-concatenation negative are presented as
negatives, never re-asserted as method. (3) **PK / l_max ceilings honest.**
(4) **Positioning** — research instrument, not production-chemistry replacement;
zero-parameter asserted as exactly that.

## Paper-58-specific watch-notes (the risk surface — ranked)

- **W1 — exact ≠ accurate [HIGHEST RISK].** The closed-form ERI engine and the decided
  census bought **decidability, not accuracy**. No sentence may frame
  closed-form / π-free / weight-one, or the `(AA|BB)` decision, as an accuracy
  improvement, or imply any molecule is made more accurate by it. Accuracy is a basis
  (max_n) axis. Any "exact ⇒ better energy" reading = MATERIAL.
- **W2 — no QC advantage from the engine [framing-zombie].** QC-1 tested NEGATIVE
  (contraction is free in qubit terms ⇒ no qubit/Pauli advantage; N3b sparsity also
  negative). The paper must not claim the engine reopens sparsity or yields a
  qubit-resource edge. The honest ceiling is the "advantages contract with system size"
  small-system line. Any qubit/Pauli/sparsity advantage attributed to the closed form
  = MATERIAL.
- **W3 — counted vs decided [C8/C3].** The `g`-row OVERALL tier stays **COUNTED**; only
  `(AA|BB)` + one-center (20.5%) are **DECIDED** (195/195, 0 accidental, 0 missed). The
  80% cross bulk stays COUNTED pending an `{E₁, ln, γ}` independence argument. No prose
  may upgrade the whole `g`-row to decided/verified.
- **W4 — hydrogenic, not Coulomb-Sturmian [C3].** The builder is hydrogenic (Z_eff per
  block), NOT shared-p₀ Coulomb-Sturmian. The abstract fix (≈line 31) must hold; no
  stale "Coulomb-Sturmian describes the builder" anywhere in the body.
- **W5 — NaH demo forfeits sparsity + is minimal-basis quality [benchmarking].** NaH
  R_eq +4.8% / 91.1% in-basis FCI is a PANEL-VERIFIED binding *demonstration* that
  explicitly forfeits the sparsity and makes no sparsity claim; it is minimal-basis
  quality, not production accuracy. Must read as such.
- **W6 — [INTERNAL THEOREM] tier + novelty honesty [C3].** The abelian residue is an
  INTERNAL THEOREM (only axial m survives two centers, by the axial-U(1) intersection).
  The group theory is elementary (diatomic Λ); the contribution is wiring it to the
  qubit-sparsity story. Assert exactly that — no "proven / first / breakthrough" beyond
  internal tier. The abstract's "We settle … necessary" is the assertive edge: verify it
  claims only the group-theoretic necessity, nothing empirical.
- **W7 — polyatomic numbers are scope, not a result [C8/benchmarking].** BeH₂
  +10.0%→+1.4%, H₂O +14.9%→+2.5%, bond angle 106.3° vs 104.5° (≈lines 188–198) are
  Gaussian-evaluator RHF figures, explicitly "reported as scope, not as a result of this
  paper." Verify the caveat is intact and the numbers are backed. The closed-form engine
  is TWO-CENTER ONLY — no claim may imply it reaches three centers.
- **W8 — transcendental tagging [Paper 18].** The `{E₁, ln, γ}` seeds are embedding-tier,
  weight-one, π-free (ln IS the Neumann Q₀ kernel; γ not independent of E₁). Verify the
  tagging is present and correct; no anonymous transcendental.

## C8 headlines (enumerated, with tiers — the frozen goalposts)

1. **Abelian residue [INTERNAL THEOREM].** Cross-center: Gaunt l-rule fails, axial m-rule
   survives; mechanism = axial U(1) is the common subgroup (Pontryagin-dual ℤ,
   additive charge conservation).
2. **Census [MEASURED].** Cross-center one-body 0.19–0.80 Ha (bond scale); permitted ERI
   density inflates 13.8× over the builder at n_max=2 (corroborated 29.8% vs 29.4%
   counted) / 15.0× at n_max=3 (counted only).
3. **`(AA|BB)` [DECIDED].** 195/195 permitted entries genuinely nonzero, 0 accidental,
   0 missed — counted density = true density on that block. `g`-row OVERALL = COUNTED.
4. **Closed forms.** All four classes closed; weight-one, π-free over `{E₁, ln, γ}`;
   ordered-ξ integral agrees with quadrature ~1e-13 / 6.3e-16, σ≠0 to 2.5e-14.
5. **NaH [PANEL-VERIFIED].** R_eq = 3.736 a₀ (+4.8%), 91.1% of in-basis FCI binding,
   3 fragment-native determinants; demo forfeits sparsity.
6. **Tapering [OBSERVATION].** Equivalent-atom-swap reduces the qubit saving
   (BeH₂ 12→9, N₂ 22→19, F₂ 22→15), inflates Pauli 2.3–3.6×.
7. **Löwdin re-reading [OBSERVATION].** The 17.9× inflation was largely already owed —
   genuine tensor inflates ~15× before any transformation.
8. **H₂ native FCI** = −1.106556606091 Ha (1s-per-centre basis); Gaussian converges onto
   it as the fit deepens (a basis-limited value, NOT the exact H₂ energy — W1).
9. **Polyatomic (SCOPE, not a result):** BeH₂ +10.0%→+1.4%, H₂O +14.9%→+2.5%, bond
   angle 106.3° vs 104.5°.

## Seeding plan (worktree only; never touches the real corpus)

K ≈ 5 planted defects, ≥1 catchable by each dimension (code / prose / citation /
synthesis), spanning the watch-notes — e.g. an **exact⇒accurate** seed (W1), a
**whole-g-row-decided** overstatement (W3), a **QC-advantage** re-assertion (W2), a
citation splice (C4), a synthesis-promotion drift (C9). M ≈ 5 known-good controls (the
verified C8 headlines) that must NOT be flagged. Tiered agents (code + citation =
Sonnet) get 2 seeds each. Answer key → `debug/qa/paper_58_seed_key.json`.

## Change log
- 2026-08-14 — **DRAFTED** by PM for PI freeze. First single-paper `/qa` target
  (Paper 58 postdates the 2026-06-28 group2 cert). Inherits criteria.md C1–C18 + the
  group2 benchmarking/guardrail deltas.
- 2026-08-14 — **DELTA-verification run (Avery-scoping edits) = DEFECTS → remediated.**
  Paste-don't-point, seeded. Deterministic whole-target green. Claims delta CALIBRATED
  (caught W3 seed) + clean on the real lineage paragraph + abstract hook. Citation delta
  CALIBRATED (caught the wrong-year + wrong-title seeds) AND surfaced a **genuine**
  defect: the `avery2004` bibitem is a fabricated composite (2006 book title + 2004
  article venue/vol/year + wrong page) — verified vs primary; **inherited verbatim from
  the CERTIFIED group2 synthesis**. Fixed in BOTH Paper 58 and the synthesis →
  `avery2006` (Avery & Avery, *Generalized Sturmians and Atomic Spectra*, World
  Scientific 2006); both recompile clean. Missing-cites sweep (advisory) returned a
  ranked list; load-bearing gaps: Löwdin 1955 (the actual non-orthogonal cofactor
  method; `lowdin1950` is a *different* Löwdin paper), the Head-Gordon NOCI pair,
  Barnett–Coulson 1951.
- 2026-08-14 — **Missing cites ADDED (PI: 3 load-bearing + top optional), each
  PM-verified vs primary sources before insertion:** Löwdin 1955 (Phys. Rev. 97, 1474 —
  the non-orthogonal cofactor method; note the already-cited `lowdin1950` is a
  *different* Löwdin paper), Thom–Head-Gordon 2009 (JCP 131, 124113) + Sundstrom–
  Head-Gordon 2014 (JCP 140, 114103) [NOCI], Barnett–Coulson 1951 (Phil. Trans. R. Soc.
  A 243, 221), Filter–Steinborn 1978 (Phys. Rev. A 18, 1). Each placed at its supporting
  claim (cofactor eval; NOCI demo; two-center Slater-integral machinery); compiles clean,
  all 5 cited (no orphans). Citation delta now clean.
- 2026-08-14 — **CERTIFIED (PI direction: "good enough for certification").** Final FULL
  certifying run WAIVED by PI; cert rests on the clean delta + the base FULL run + green
  deterministic layer (see STATUS honest cap). Committed to `work/sparsity-boundary`.
  Standing follow-up (PI-noted, not blocking): sweep the corpus's Avery/hyperspherical
  citations (Papers 16/34, other syntheses) for the same `avery2004` composite that this
  run found had propagated from the certified group2 synthesis.
- 2026-08-15 — **FULL certifying run = PASS** (the PI-waived full run, now run). Fresh
  seeded worktree over the final combined text; 6 seeds (2 code / 2 citation / 1 claims /
  1 synthesis) all caught, specificity clean; **zero verified material defects** beyond
  seeds across code / claims / citations / synthesis + completeness-critic; deterministic
  C10–C18 all green. Honest cap DISCHARGED. NIT-tier residue (fix-on-sight, non-blocking):
  (a) the H2 native-FCI value and the 195 census total are backed only by OWED `debug/`
  drivers (matrix-logged), not pytest asserts — recommend pinning both in `tests/`;
  (b) uncited provisional benchmarks (D_e = 1.961 eV flagged provisional; Clementi
  −161.12); (c) the `(ab|cd)` closed-form display equation lacks a `\label` (§13.4a);
  (d) "conservation law" / "readable in advance" rhetoric slightly ahead of its hedge.
  Seed key `debug/qa/paper_58_full_seed_key.json`; worktree removed, no seed leaked.
- 2026-08-14 — **FROZEN (PI: Path A + "freeze as written").** Review begins:
  deterministic `--gate group2` on the real corpus; seeded worktree + panel for the
  four LLM dimensions + completeness-critic.
- 2026-08-14 — **FIRST-CERT FULL run = FAIL (trustworthy).** Panel FULLY CALIBRATED:
  sensitivity **6/6** seeds caught (code S1/S2 by recompute; claims S5 W1; citation
  S3/S4; synthesis S6 W3), specificity clean (**0/5** controls false-flagged). Three
  LLM dimensions (claims / synthesis / citations) clean beyond seeds. Genuine
  verified defects: (a) **C11 deterministic** — Paper 58 bibitem titles for P19+P20
  wrong → **FIXED at source, C11 re-run PASS**; (b) **code MATERIAL** — §I.A
  polyatomic-extension numbers (BeH₂/H₂O, M=19, 106.3°, +1.4%/+2.5%) have **no
  backing artifact** anywhere in the repo (verified incl. uncommitted `debug/poly3*`)
  → **RAISED TO PI**. SMALL (logged): bra/ket certificate "holds on every entry"
  under-tested (panel 6-case test is sound); Thm 1(iii) no dedicated test
  (argument-tier, hedged, already in claim_test_matrix). Citation BeH⁺-vs-BeH₂ flag
  DISMISSED (Paper 17 title genuinely "BeH⁺"). Seed key
  `debug/qa/paper_58_seed_key.json`; worktree removed, no seed leaked. **Path to
  cert:** PI resolves the polyatomic finding → delta-verification re-run → FULL
  certifying run.
