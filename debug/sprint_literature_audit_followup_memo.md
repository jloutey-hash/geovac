# Sprint literature-audit-followup — canonical memo

**Date:** 2026-06-03
**Version:** v3.45.2
**Scope:** PI-requested comprehensive literature audit (4 parallel phases) followed by full follow-through — verification gate against hallucinated citations, 21 citation additions across 14 papers, 5 substantive framing edits, 2 incidental pre-existing-bug fixes (Paper 51 `\nmax`, Paper 36 `\to` in `\text{}`).

---

## Why this sprint

PI asked: comprehensive sweep of the literature to (a) check whether GeoVac findings can be attributed to prior work, (b) identify connections worth building, (c) lay precedence stakes for the most novel claims. No prior comprehensive sweep had been done — only per-paper or per-sprint surveys (e.g., RH-I 2026-04, WH1 Marcolli-vS lineage tracking).

## Methodology — 4 parallel phases + verification gate

**Phase 1 (precedence):** Papers 45-50 — math.OA standalones with "first published" claims. Aggressive 2024-2026 sweep of arXiv math.OA / math-ph for Lorentzian propinquity, Krein spectral triple convergence, OSLPLS / synthetic-Lorentzian-NCG follow-ups. Targeted lookback to Latrémolière / Hekkelman / Mondino-Sämann / van den Dungen lineages.

**Phase 2 (connections):** Master Mellin engine, projection taxonomy, spectral triple construction, case-exhaustion theorem. 2020-2026 aggressive on Connes-Chamseddine / Marcolli-vS / Perez-Sanchez lineage; targeted on motivic/period theory.

**Phase 3 (physics-side):** Paper 14 qubit scaling, Paper 22 angular sparsity, Paper 36 Lamb shift, gravity arc (Papers 51, 53). 2018-2026 on qubit/gravity; 2010-2026 on Lamb compilations.

**Phase 4 (long-tail):** Bargmann-Segal lattice (Paper 24), nuclear shell (Paper 23), Hopf-graph Ramanujan (Paper 29), SU(2)/SU(3) Wilson (Papers 25/30/41). 2015-2026 broad sweep.

**Phase 2a verification gate (added at PI request after Phase 1):** independent WebFetch verification of every proposed arXiv ID against arxiv.org abstract page — checks paper exists, title matches, authors match. Defensive gate before any citation lands in published `.tex` files.

## Headline finding: hallucination signature

**Six of 33 verified citations (18%) had MISMATCHED authors but topically-adjacent content.** Pattern:
- Pre-2024 IDs (10 in sample): 100% VERIFIED clean
- 2024-2026 IDs: 6 of 23 had wrong author attribution
- **Every mismatch resolves to a real paper on the same topic** — the hallucination signature is "grab a nearby plausible arXiv ID, invent an author attribution that fits the topic"
- Three of six are famous-name-swaps (Chirco-Josset-Rovelli, Iazzi-Glaser, "Stetcu group") to better-known researchers in the relevant field
- Pre-flagged 2604.08349 (suspicious month-code) confirmed as one of the six

**Practical takeaway:** content-only citation checks miss this entirely. Author-level verification is required when sub-agents propose citations near or past their training cutoff. The verification gate caught all six before they landed in publications.

## Per-phase results

### Phase 1 — precedence audit

Memo: `debug/lit_survey_phase1_precedence_memo.md`

- All 6 Papers 45-50 novelty claims survive the June 2026 audit
- Most important new citation: **Mondino-Sämann arXiv:2504.10380 v4 (Dec 2025)** — synthetic Lorentzian GH convergence; categorically separate but Paper 48 explicitly bridges to MS framework
- Single most-relevant concurrent work: **Latrémolière arXiv:2512.03573 (Dec 2025)** — already correctly cited in Papers 45, 47, 48 (verified)
- Paper 49 Datta max-divergence fix verified IN PLACE (30+ references throughout the .tex)
- One framing recommendation: Paper 48 §1 explicit "what we take from MS / what we contribute to MS" paragraph

### Phase 2 — connection-finding

Memo: `debug/lit_survey_phase2_connections_memo.md`

- **Highest cite-to-cost ratio addition:** Perez-Sanchez arXiv:2508.17338 (Aug 2025) — establishes Marcolli-vS continuum limit is Yang-Mills without Higgs, matching Sprint H1's discrete-substrate verdict. Already in CLAUDE.md WH1 but missing from Paper 32 §VIII.B
- **Sprint candidate surfaced (1 week):** are GeoVac discrete SD coefficients on S³ mixed-Tate periods (Fathizadeh-Marcolli arXiv:1611.01815 lineage)? Either outcome publishable
- Natural collaboration cluster: Hekkelman / van Suijlekom / Latrémolière / Farsi (4 directly-adjacent papers in 2024 alone)
- Theme 2 (projection taxonomy) returned weaker matches than expected — Paper 35 §VI is uncommon territory; Tener 2025 (BW for non-unitary CFTs) is the closest theorem-grade analog

### Phase 3 — physics-side audit

Memo: `debug/lit_survey_phase3_physics_memo.md`

| Paper | Verdict | Headline |
|:------|:--------|:---------|
| 14 (qubit scaling) | NOVEL-BUT-NARROW | Defensible vs raw JW; not benchmarked vs THC/DF at fault-tolerant scale. Paper's own honest disclosure is the model |
| 22 (angular sparsity) | NOVEL | Machinery is textbook; universality across 5 potentials is genuinely new |
| 36 (Lamb shift) | NEEDS-STRONGER-BASELINE | 0.534% is 3 OoM looser than precision-QED's kHz; structural-not-precision framing is right; sharpen it |
| 51, 53 (gravity arc) | NOVEL mechanism / DUPLICATED observation | Chamseddine-Connes 2008 already had two-term exactness "up to astronomically small corrections" — GeoVac sharpens, doesn't discover. Paper 51 framing risk |

Most-critical missing citation: **Yerokhin et al. PRL 2024 (arXiv:2411.12459)** — state-of-the-art two-loop self-energy at low Z; corpus-wide gap for Paper 36.

### Phase 4 — long-tail connection scout

Memo: `debug/lit_survey_phase4_longtail_memo.md`

- **Strongest single bridge candidate:** Garofalo et al. arXiv:2311.15926 + 2503.03397 — SU(2) Hamiltonian lattice gauge theory built by partitioning S³ ≅ SU(2). Same Wilson Hamiltonian on the same group manifold as Paper 30, different community. Urbach (Bonn) / Jansen (DESY) are the natural conversation partners
- **Open-question hit:** Higgs-Pickrell arXiv:2503.23549 (March 2025) — Pickrell's open question for $d \ge 3$ on SO(d)-invariant HO is exactly Paper 24's $d=5$ case
- Direct deuteron VQE competitors for Paper 23: Sarma-Stevenson (Surrey) arXiv:2510.02124 and Pillai-Ramanan-Balakrishnan-Lakshmibala arXiv:2509.08948

### Phase 2a — verification gate

Table: `debug/citation_verification_table.md`

Ship-list after verification:
- 23 cleanly VERIFIED as claimed
- 4 borderline-VERIFIED with cleanups (#16 Cavalletti-Mondino: add arXiv:2004.08934; #17 Farsi-Latrémolière-Packer: add arXiv:2301.00274 distinct from 2302.09117; #21 Caesura et al.: cite 11-author list not "PsiQuantum-BI" shorthand; #24 Garofalo: title is "S₃ partitionings" not "S³ ≅ SU(2)")
- 6 MISMATCH with substituted attributions (#3 → 2507.01482 Behrndt-Holzmann-Stelzer-Landauer; #4 Rotondo 2026 not Chirco-Josset-Rovelli; #13 Tener 2025 not "Camassa et al."; #14 Morinelli 2024 not "Crisand"; #20 Homšak-Veroni 2024 not Iazzi-Glaser; #27 Sarma-Stevenson 2025 not "Stetcu group")

## Citation additions (21 across 14 papers)

Dispatched as 3 parallel sub-agents (A: math.OA + gravity; B: foundations + nuclear + Ramanujan; C: qubit + gauge + projection + Lamb). All three reported clean compiles. Substituted attributions applied per verification table.

| Paper | Group | New cites |
|:-----:|:------|:----------|
| 32 | g1 | van Suijlekom 2409.02773 (other 3 already present) |
| 47 | g1 | Behrndt-Holzmann-Stelzer-Landauer 2507.01482 |
| 48 | g1 | Braun-Sämann 2506.10852 (MS 2504.10380 already present) |
| 49 | g1 | Rotondo 2604.08349 |
| 51 | g5 | Homšak-Veroni 2404.11670 |
| 18 | g3 | Fathizadeh-Marcolli 1611.01815, Rejzner 1603.02748, Connes-Marcolli math/0409306 |
| 23 | g4 | Sarma-Stevenson 2510.02124, Pillai et al. 2509.08948 |
| 24 | g3 | Higgs-Pickrell 2503.23549 |
| 29 | g1 | Matsuura-Ohta 2403.07385 PTEP 2024 |
| 14 | g4 | Caesura et al. 2501.06165 |
| 30 | g5 | Garofalo et al. 2311.15926, Jakobs et al. 2503.03397 |
| 34 | g6 | Tener 2025 2506.10625, Cavalletti-Mondino 2004.08934 |
| 35 | g6 | Morinelli 2024 2412.20410, Tener 2506.10625 (forward pointer) |
| 36 | g5 | Yerokhin-Harman-Keitel 2411.12459 |

Bibitem-key collision avoided in Paper 14 (existing `caesura2025` was a different paper; used `caesura2025_tensor_factorization`).

## Framing edits (5 papers)

All edits made in main session with prose deltas shown to PI before commit:

1. **Paper 51 abstract** — promoted Chamseddine-Connes 2008 "uncanny precision" into abstract; reframed contribution as "Bernoulli-mechanism sharpening" not de novo discovery. Mitigates reviewer-pushback risk identified by Phase 3.
2. **Paper 48 §1** — added explicit "what we take from MS / what we contribute to MS" paragraph at start of `\subsection{Main result}`. Most-precedence-sensitive paper of the six (Lorentzian propinquity space is GeoVac-only as of June 2026, but MS team is publishing at sprint cadence).
3. **Paper 50 abstract** — added explicit sentence reframing as verification-with-structural-insight (M2/M3 orthogonal decomposition is novel; KPS 2011 F-theorem values are not). Defends against reviewer reading (A) as new F-theorem.
4. **Paper 35 §7 Falsifiable prediction** — added forward-pointer to BW reading immediately after the prediction statement; named Bisognano-Wichmann theorem (with Morinelli 2024 + Tener 2025) as closest theorem-grade analog.
5. **Paper 36 §1 Introduction** — added explicit `\paragraph{Scope of the claim}` acknowledging 3-OoM gap to precision-QED kHz frontier; sharpens existing structural-not-precision framing.

## Incidental bug fixes (pre-existing, surfaced during compile verification)

- **Paper 51:** added `\newcommand{\nmax}{n_{\max}}` to preamble (matching project-standard from Papers 38, 39, 42, 43). Paper had a pre-existing fatal compile error from undefined `\nmax` macro — likely existed for some time, surfaced when Agent A tried to verify a citation insertion.
- **Paper 36:** fixed pre-existing math-mode error `\Delta_\text{LS-1 \to LS-6a}` → `\Delta_{\text{LS-1} \to \text{LS-6a}}`. Same family of bug, surfaced when verifying Phase 3.5 edit.

Both are defensible 1-line cleanups; both were necessary to verify edits in their respective papers; both are flagged here for visibility.

## Compile verification

All 5 papers touched in Phase 3 compile clean three-pass after Phase 4 fixes:
- Paper 51: 37 pages
- Paper 48: 29 pages
- Paper 50: 17 pages
- Paper 35: 19 pages
- Paper 36: 15 pages

Plus the 9 additional papers touched by citation sub-agents (Phase 2b) — all reported clean compiles at agent completion.

## Carved-out — deliberately not in this sprint

- **Mixed-Tate period test** (1-week sprint candidate). Are GeoVac discrete SD coefficients on S³ mixed-Tate periods? Either outcome (yes → Paper 32 §VIII M2 inherits as corollary; no → discrete sector more restricted than continuum) is publishable. PI decision pending.
- **Collaboration outreach** (Garofalo/Urbach group at Bonn/DESY; Hall at Notre Dame for Bargmann-Segal; Stetcu/Johnson nuclear VQE comparison; Hekkelman cluster). PI direction needed.
- **Paper 14 §sec:ft_gaussian** explicit BLISS-THC comparison row (Caesura et al. cited but full comparison-table population would need additional resource estimate; defer).
- **Pre-existing math errors elsewhere** in the corpus (Paper 14 missing figures, Paper 30 remark env errors) — flagged by Phase 2b agents, not in scope here.

## Hallucination-pattern learnings (preserve as discipline)

1. **Sub-agents past their training cutoff hallucinate arXiv IDs with high prevalence (18% in this sweep).** Always run a verification gate when sub-agents propose post-training-cutoff citations.
2. **Hallucinations are topical, not random.** Content-checks alone will miss them; author-level verification (WebFetch the abstract page, cross-check author names) is required.
3. **The pattern is consistent enough to predict:** for any future citation sweep with sub-agents, expect ~15-20% of recent-year IDs to need correction. Budget verification time accordingly.
4. **Substitute papers (right topic, corrected attribution) usually serve the citation purpose** — the agent's content-matching was right even where author attribution was wrong. Verification + substitution is more efficient than re-search.

## Status

| Phase | Status |
|:------|:-------|
| Phase 1 (verification of existing fixes) | DONE — both Datta and Latrémolière already in place |
| Phase 2a (citation verification gate) | DONE — 6 mismatches caught + corrected |
| Phase 2b (citation additions) | DONE — 21 cites across 14 papers, all compile clean |
| Phase 3 (framing edits) | DONE — 5 papers, prose deltas shown to PI |
| Phase 4 (release wrap) | IN PROGRESS — this memo + CHANGELOG + CLAUDE.md §2 + version bump + commit + tag + (push pending PI confirmation) |

## Files modified

26 papers + 5 lit-survey memos + 1 verification table + this sprint memo. Full list in git status; aggregate diffstat is ~37 file changes, ~600 line insertions (excluding regenerated PDFs).
