# Corpus Accessibility & External Validation Plan

**Created:** 2026-06-09. **Approved:** 2026-06-09 (PI, all four gates; Option C on the Phase-1 finding).

**Status (2026-06-10):**
- **Phase 1 — COMPLETE.** The refutation pass falsified Paper 45's main theorem (annihilation theorem; fabricated citation; L2 conflation). Option C executed: all 8 affected papers corrected **in place** (PI directive 2026-06-10: single source of truth, no v2/erratum scaffolding; git/Zenodo is the version record). Falsifier frozen in `tests/test_p45_kplus_degeneracy.py`.
- **Phase 2 — COMPLETE.** 2a: N1 + N2 front-door notes (`docs/outreach/`, both compile clean). 2b: `docs/vocabulary_translation.md` (33 rows). 2c: `docs/claims_register.md` (20 rows). 2d: standalone audit (`debug/p38_p45_standalone_audit.md`) + all rewords applied. 2e: README (claims-status section, reading paths, stale claims fixed), field guide (conditional/degeneracy framing, rate display corrected), Group-1 synthesis (status note + 13 corrections), `.zenodo.json` description corrected. Bonus: Paper 28 Theorem 3 proof factor-2 error caught while drafting N2, fixed, verified to 40 digits.
- **Phase 3 — COMPLETE (2026-06-10).** CLAUDE.md §1 viability sentence reworded (PI-directed); README badge "Production"→"Research" + matched-qubit caveats at both headline sites; Paper 14 abstract caveat inserted (Paper 20 already compliant, 6/6 sites); corpus-wide α-mention audit CLEAN (zero fixes needed — §13.5 discipline held).
- **Phase 4 — PREPARED, sends pending (2026-06-10).** Full send kit at `docs/outreach/phase4_send_kit.md`: recipient ladders (N1: vS → Latrémolière → early-career; N2: Brown/Kleinschmidt), four email drafts, pre-registered outcomes, send checklist. The S1-vs-S2 decision flagged in the kit is now **moot in the gap direction**: the P38 gap-closure landed (see next line), so the strongest send posture — repaired unconditional theorem in hand — is available. N1 was updated by the cascade to state the unconditional result. Release precedes any send.
- **Phase 5 — ACTIVE; P38 gap-closure COMPLETE (2026-06-10, v3.109.0).** The standing-exception work closed both named gaps via the translation-seminorm reframing (G2: truthful-substrate kernel condition via Schur + per-band injectivity; G1: exact-fit spinor lifted state, defects = Fejér smoothings). Paper 38's theorem is **unconditional**; cascade applied to P45/P32/claims register/README/field guide/N1; falsifier `tests/test_p38_action_seminorm.py`. The arc freeze otherwise holds; the remaining standing exception is the Lorentzian repair target (Toeplitz temporal algebra + non-compressing Krein seminorm).
**Goal:** Convert the corpus from internally-coherent to externally-checkable, and run the first external-validation experiments. Addresses the four standing concerns from the 2026-06-09 honest assessment: (1) unverified mathematical superstructure, (2) corpus unreadability for outsiders, (3) tool-positioning honesty, (4) zero external expert contact.

---

## 0. Design principles

- **Shrink the claim surface; never grow it during this arc.** No new research arcs while Phases 1–4 run (Phase 5 freeze). Exceptions only for gaps named by the adversarial pass or by an external replier.
- **One front door per audience, ≤ 4 pages, zero project-internal vocabulary.** No reader should need CLAUDE.md, the field guide, or 57 papers to evaluate one claim.
- **The flagship gets hardened adversarially before anyone external sees it.** Finding a gap in Phase 1 is a success outcome, not a failure.
- **Honesty tiering made externally legible.** The internal verification culture (bit-exact panels, falsifiers, dead-end registry) becomes a visible artifact (claims register), not a private practice.
- **Success metric for the whole arc:** one content-engaged reply from one domain expert. "This is wrong because X" counts as success. "This is known, see Y" counts as success (calibration). Only silence is failure.

---

## 1. Phase 1 — Adversarial hardening of the flagship pair (Papers 38 + 45) [~3–5 sprint-days]

The flagship is the propinquity chain: Paper 38 (Riemannian GH-convergence, the foundation) + Paper 45 (Lorentzian propinquity, the first-in-literature novelty claim). Paper 46 joins if 45 survives.

**1a. Referee-grade proof write-out.** Expand every step previously labeled "straightforward" or "bookkeeping":
- P38 L5 (Latrémolière tunneling-pair assembly) — the reach/height bookkeeping in full.
- P45 L2 (Bożejko–Fendler cb-norm equality on the amenable product) — full proof, not citation-by-analogy.
- P45 L5 K⁺-preservation under Latrémolière 2017/2023 §5 — the named gap, closed line-by-line (Phase 1B-A claimed this 2026-05-24; re-derive fresh).
- Rule: no "clearly," no "straightforward," no "verbatim" survives in a load-bearing lemma. Each step gets a full argument or an explicitly named gap.

**1b. Adversarial refutation pass.** Fresh-context reviewers prompted to *break* each lemma, not confirm it:
- Hunt the known failure mode: citing an imported theorem whose hypotheses were never checked against our objects line-by-line (quasi-Leibniz constants, *-closure of the operator system, completeness of the Lip-pair, admissibility conditions on tunnels).
- Attempt small-n_max counterexamples against each lemma's statement as written.
- Separately: a "does the theorem statement claim more than the proof delivers" pass (qualitative-rate vs quantitative-rate, K⁺ scope, substrate scope).

**1c. Hypothesis checklist appendix.** A table per flagship paper: every imported theorem (Latrémolière 2013/2017/2023, Connes–vS 2021, Bożejko–Fendler, Leimbach–vS 2024) × every hypothesis × where our construction satisfies it. Doubles as verification and as a referee-trust device.

**1d. Freeze numeric falsifiers as tests.** Mostly exists (Riemannian-limit bit-exact recovery, panel monotonicity); confirm coverage and pin in `tests/`.

**Gate:** if 1b finds a real gap → fix or descope honestly before any outreach. Outreach does not proceed on a flagship that hasn't survived its own refutation pass.

---

## 2. Phase 2 — The accessibility layer [~3–4 sprint-days, partly parallel with Phase 1]

**2a. Front-door notes** (2–4 pages each, self-contained, standard notation only):
- **N1 (math.OA):** "Gromov–Hausdorff convergence of truncated spectral triples on SU(2), and a Lorentzian propinquity theorem." Theorem statements, the K⁺ trick in one paragraph, precise novelty claim relative to Connes–vS 2021's explicit deferral and the Riemannian-only Latrémolière lineage, pointers to Papers 38/45/46. No "GeoVac," no "Fock graph," no "M1/M2/M3," no "WH1."
- **N2 (periods — Brown/Kleinschmidt):** two or three self-contained, ten-minute-checkable identities — candidates: the χ₋₄ closed form D_even(s) − D_odd(s) = 2^(s−1)(β(s) − β(s−2)); ζ_{D²}(2) = π² − π⁴/12; the pure-Tate Seeley–DeWitt sub-ring statement on S³ — followed by ONE precise question. Not the motivic-Galois framing; an expert can supply their own framing if the identities interest them.
- **N3 (QC/chemistry, optional, later):** honest encoding-theory positioning — structural sparsity results + the accuracy walls, presented together.

**2b. Vocabulary translation table** (`docs/vocabulary_translation.md`). Internal term → standard term, e.g.: "Fock graph" → truncated Sturmian/Peter–Weyl basis; "master Mellin engine M1/M2/M3" → heat-kernel Mellin mechanisms at operator order k = 0, 1, 2; "focal length" → projection/scale parameter; "Layer 1/Layer 2" → combinatorial data vs projected observables; "PROVEN" → internally verified lemma chain, unreviewed. Used to scrub N1/N2 and the flagship intros.

**2c. Claims register** (`docs/claims_register.md`, surfaced from README). Every headline claim an outsider would meet first (~15 rows to start) → status ∈ {symbolic proof, panel-verified at finite n_max, internal theorem (unreviewed), observation, conjecture} × falsifier × where verified in `tests/`. This converts the internal honesty culture into the externally legible artifact that earns trust.

**2d. Standalone audit of the flagship papers.** Confirm P38/P45 are readable with zero corpus context: all notation defined in-paper, GeoVac motivation demoted to a remark, no load-bearing references to internal memos or other corpus papers.

**2e. Field guide + README restructure.** Lead with the defensible core in this order: conformal equivalence (Paper 7) → sparsity theorems (Papers 22/14) → propinquity arc (38/45–49) → transcendental-entry observation (18/34/35). The packing origin story and α move to an honest back-section ("origins and open coincidences"). Add per-audience reading paths and surface the claims register at the top of README.

---

## 3. Phase 3 — Repositioning edits [~1 sprint-day]

- **Viability-case language** (README + CLAUDE.md §1): "usable, benchmarked computational tool" → discretization/encoding-structure research program with a tool artifact. **Touches §1/§1.5 — PI approval required; PM drafts the diff only.**
- **Papers 14/20 pass:** every sparsity headline co-located with its accuracy caveat (benchmarking rule mostly does this; verify and close gaps).
- **α references audit:** every corpus mention of K = π(B + F − Δ) leads with observation status + twelve eliminated mechanisms. No mention carries motivational weight for other arcs.

---

## 4. Phase 4 — The external experiment [calendar-gated; sends only after Phase 1 passes]

**4a. math.OA recipient ladder (N1):**
1. W. van Suijlekom (Radboud) — spectral truncations are his program; P45 closes a deferral his own paper states.
2. F. Latrémolière — propinquity is his framework.
3. Early-career rung: Leimbach, Hekkelman, McDonald — higher reply base-rates.
Email shape: ≤ 5 sentences + N1 attached + one paper PDF. The ask: "Is Theorem X known, or wrong? We believe it is the first Lorentzian instance; happy to be corrected."

**4b. Brown/Kleinschmidt (N2):** send the identities note, not the corpus. One question, checkable in minutes. Recommendation: send ~2 weeks after the math.OA rung — any engagement (even critical) from 4a strengthens this contact, and Brown is a one-shot high-value contact not to be spent on an unhardened corpus.

**4c. arXiv (math.OA) for P38 + P45.** Timestamped, discoverable, citable; the de facto external venue for this material (Zenodo DOIs have near-zero math.OA discoverability). Requires endorsement as an unaffiliated author — the endorsement request is itself a micro-instance of the experiment. Not a journal submission; consistent with the Zenodo-first policy. **PI decision.**

**4d. Pre-registered outcomes.** Content-engaged reply (incl. refutation) = success. "Known, see X" = success (calibration). 3 weeks silence → next rung. Outcomes logged like sprint negatives in CHANGELOG.

---

## 5. Phase 5 — Freeze discipline

No new research arcs while Phases 1–4 run. Exceptions: gaps named by Phase 1b, requests from external repliers, and mechanical maintenance. Rationale: each additional internally-verified paper currently subtracts external credibility by growing the surface a reviewer must distrust.

---

## 6. Sequencing

```
Week 1:    Phase 1a/1b (flagship hardening)  ||  Phase 2a-N1 draft + 2b translation table
Week 1-2:  Phase 1c/1d  ||  Phase 2c claims register + 2d standalone audit
Week 2:    Phase 2e field guide/README  +  Phase 3 repositioning (PI gate on §1 language)
Week 2+:   Phase 4a sends (if Phase 1 passed) → 4b ~2 weeks later → 4c in parallel
Ongoing:   Phase 5 freeze
```

First concrete sprint: **Phase 1a/1b on Paper 45's L2 + L5 lemmas** (highest-risk load-bearing steps), with N1 drafted in parallel.

---

## 7. PI decision gates

| # | Decision | Recommendation |
|---|----------|----------------|
| 1 | Flagship + outreach order | P38+P45 flagship; math.OA rung first, Brown/Kleinschmidt ~2 weeks later |
| 2 | arXiv posting (P38, P45) | Yes — discoverability is the point; endorsement hurdle is informative |
| 3 | §1/§1.5 viability-case rewording | Approve PM-drafted diff (PM cannot edit §1.5 unilaterally) |
| 4 | External-facing status vocabulary | Rename "PROVEN" → "internally verified (unreviewed)" in all outward docs; internal register unchanged |
