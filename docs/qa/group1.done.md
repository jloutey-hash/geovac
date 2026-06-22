# Group 1 (Operator algebras / NCG) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group1-specific scope + deltas + the branch-specific
> criterion C14.

> **STATUS: DRAFTED 2026-06-16 — NOT YET FROZEN.** Sequenced *after* group3
> bite 2. Awaiting PI confirmation of scope + first bite before FREEZE.

**Scope (non-trunk group1):** Papers **29, 39, 40, 42, 43, 44, 45, 46, 47, 48,
49, 50, 52, 53** + the **group1 operator-algebras synthesis**. Trunk papers
**32, 38** taken as already-certified (`/qa trunk` PASS); in scope only where a
group1 paper restates them (C7).

**Deterministic `--gate`:** `group1`.

## Branch deltas (the only non-inherited content)

- **C14 — Descope / partial status accuracy (branch-specific; the defining
  criterion).** Every DESCOPED or PARTIAL paper presents withdrawn claims at
  current status:
  - Papers **45, 46 = DESCOPED** — the "first Lorentzian propinquity
    (convergence) theorem" is **retracted**; Paper 45's headline is the K⁺
    annihilation theorem; no body text presents weak/strong-form Lorentzian
    propinquity convergence as established.
  - Papers **47, 48, 49 = PARTIAL** — surviving arrows stated as surviving (47
    norm-resolvent + three-carrier; 48 bridge as conditional design; 49
    cocycle-deficit / OSLPLS algebra); descoped metric-level / Λ-inheritance
    claims flagged descoped.
  - The **product-carrier S³×S¹ convergence** is presented as
    **signature-agnostic (Euclidean), NOT a Lorentzian claim**.
  - Each carries an in-place Status note; no pre-descope claim cited as live.
  *(If a second branch later needs C14, promote it into `criteria.md`.)*
- **C4 high-priority for this branch.** The math.OA citation apparatus
  (Connes–vS, Latrémolière, van Suijlekom state-space GH, Marcolli,
  Camporesi–Higuchi, Bizi–Brouder–Besnard, Mondino–Sämann, Datta) is the
  corpus's highest fabrication surface and *already produced a real defect here*
  (the nonexistent "Latrémolière Thm 5.5 / Def 3.4" refs that drove the Paper 45
  descope; the Paper 38 `avery_wen_avery` v4.14.1 fix). Every cited theorem/def
  *number* is verified.
- **C7 (trunk-dependent status).** Paper 38/WH1 as PROVEN scoped to the van
  Suijlekom state-space GH distance; Paper 32 axioms / κ at current tiers.
- **C8 (headline honesty), per-paper.** Paper 29 = Ramanujan + integer-algebraic
  Ihara zeros; Paper 40 = 4/π rate **universal** (derived-numerics-pinned + PRV
  tightness, not full symbolic proof beyond the verified panel); Paper 50 =
  **bit-exact** CFT-on-sphere F-theorem match (arithmetic fact, not a CFT
  derivation); Paper 45 headline = the **annihilation theorem** (a NEGATIVE
  result), never a convergence theorem.

## Proposed first bite (PI to confirm at FREEZE)

- **(recommended) Papers 29, 40, 50 + synthesis** — the ACTIVE math.OA backbone
  (Ramanujan / universal 4/π rate / CFT F-theorem); calibrate the panel on the
  non-descoped core, C9 scoped to the synthesis's coverage of those three plus
  its Lorentzian-descope narrative (coherence-passed v4.19.0).
- **Bite 2 (proposed):** the **Lorentzian cluster 42–49** — the C14-heavy,
  highest-risk subset, after the panel is calibrated.
- *Alternative: attack the Lorentzian cluster first if you'd rather hit the risk
  head-on while the descope is fresh.*

## Change log
- 2026-06-16 — **DRAFTED** by PM for PI review (third pre-registered `/qa`
  target). Adds C14 (descope/partial status); flags C4 high-priority. Scope
  excludes trunk-certified 32/38. **Awaiting PI freeze** (sequenced after group3
  bite 2).
- 2026-06-17 — **Paper 29 (single-paper bite) = FAIL → FIXED** (v4.21.3).
  PI-scoped to one paper (token slow-roll); dimensions deterministic + code +
  claims + citation (synthesis excluded, single paper). Calibrated panel
  (sensitivity 4/4, specificity 4/4). Science sound; fixed a degree-arithmetic
  error (84→80), two C4 citation conflations (matsuura → JHEP 09(2022)178,
  yakaboylu → H-P-Hamiltonian arXiv:2309.00405), the McKenzie initial, + closed
  the S⁵-N3 closed-form coverage gap (new `test_s5_N3_zeta_matches_paper_closed_form`).
  Seed key `debug/qa/group1_seed_key.json`.
- 2026-06-17 — **Bite A (Papers 40, 50 + synthesis) = FAIL → REMEDIATED**
  (v4.22.0). Calibrated panel (sensitivity 4/4, specificity 5/5). Genuine
  defects fixed: 2 wrong-ID + 1 fabricated citation (beccaria→1406.3542,
  hartman→1807.11401, hekkelman-T^d removed→Leimbach-vS); paper_50 F-theorem
  false-positive backing → genuine `test_paper50_f_theorem.py` (ζ′(0) from
  framework spectrum → KPS 1e-82); 3 synthesis/paper-50 Lorentzian/propinquity
  zombies; predictive-CFT + over-rigor calibration. **Paper 40 rank≥2
  universality**: pruned `l2_universal_rate_memo` resurrected from git history
  → genuine computation but **fit-sensitive** extraction; calibrated to rank-1
  rigorous + rank≥2 A-over-B robust (named gap), backing
  `test_paper40_universal_rate.py`. Seed key `group1_A_seed_key.json`, memo
  `debug/sprint_qa_group1_biteA_memo.md`. **Remaining: Bite B (39, 42–49, 52,
  53) — Lorentzian cluster, smaller sub-bites.**
- 2026-06-18 — **Bite B sub-bite 1 (Papers 42, 43, 44 + synthesis) = FAIL →
  REMEDIATED** (v4.23.0). 10-agent panel calibrated (sensitivity 5/5,
  specificity 5/5). **Two reviewer LARGE findings OVERTURNED by PM verification**
  (§9 reconcile): paper_42's "derived finding is false" (reviewer omitted the
  β factor — D_W generator is 2π·D_W, does NOT close; finding CORRECT, no
  reframe) + §10 "descoped-zombie" (over-flag — genuine finite-cutoff Krein
  closure; fixed the stale intro disclaimer instead). Genuine fixes: 6 wrong-ID/
  fabricated citations (hekkelman_mcdonald2024 fabricated ×3, hekkelman2022,
  zhu_casini2020 authors, latremoliere2018 ×3, avery, devastato) + theorem-
  number verification (van den Dungen Prop 4.1 + Nieuviarts Def 2.2 GROUNDED) +
  Paper-38 propinquity→state-space-GH labels + P39 zombie title + synthesis
  Paper-45 degeneracy framing + paper44 bibitem. All compile clean; C11/C13/C14
  PASS. Seed key `group1_B1_seed_key.json`, memo
  `debug/sprint_qa_group1_biteB1_memo.md`. **Remaining: sub-bite 2 (45–49),
  sub-bite 3 (39, 52, 53).**
- 2026-06-18 — **Bite B sub-bite 2 (Papers 45–49) = FAIL → REMEDIATED** (full
  sweep, PI direction; v4.24.0). 13-agent panel calibrated (sensitivity 5/5,
  specificity 5/5). Heaviest FAIL: a systematic **"Status-note-only descope"** —
  all 5 papers correctly recorded the descope in a Status note but re-asserted
  the descoped Lorentzian-metric convergence as established throughout theorem
  bodies/abstracts/appendices (p46 Appendix B "closes G1′", p47 header/§7, p48
  T3/T6/§2.2, p49 B4′/aggregate). Propagated each paper's own Status-note framing
  into ~20 zombie spots; fixed citations (devastato ×2, hekkelman ×3, bizi
  confabulated title/venue p47, latremoliere ×; the nonexistent Latrémolière Thm
  5.5 re-pointed in p46); theorem-numbers (van den Dungen/Nieuviarts/Mondino–
  Sämann) verified GROUNDED. All 5 compile clean; C11/C13/C14 PASS. Deferred:
  NIT-grade citation-title slips (real IDs). Seed key `group1_B2_seed_key.json`,
  memo `debug/sprint_qa_group1_biteB2_memo.md`. **Remaining: sub-bite 3 (39, 52,
  53).**
- 2026-06-18 — **Bite B sub-bite 3 (Papers 39, 52, 53) = FAIL → REMEDIATED**
  (v4.25.0; **completes group1**). 9-agent panel, calibrated (sensitivity 5/5,
  specificity held). Two SERIOUS findings remediated **rescue-first** (PI
  direction), both PM-verified vs code/drivers: **(p39)** the boxed C₃<1
  Pythagorean keystone is FALSE + footnote fabricated — real-operator probe
  gives Pyth-ratio 2.0 (maximally violated) / triangle-ratio 1.0 (tight), so
  C₃→√2 (triangle); convergence theorem SURVIVES, "<1" downgraded; §10 k-fold/
  master theorems "k-INDEPENDENT"→√k-sketch; new `test_paper39_triangle_tight.py`
  + `c3_full_triangle_bound`; 4 citation fixes. **(p53)** disk positivity(s≥2)/
  rate Λ^{−1.30} "verified numerically" are refuted/phantom (drivers: min eig
  −0.13, Λ^{+0.07}) — rescue SUCCEEDS via the boundaryless plane (Bochner–Riesz
  s≥½, Λ^{−0.6..−0.9}); new `test_paper53_plane_bochner_riesz.py`; thm:interior
  (b),(c) withdrawn on the finite disk. p52: 4 C7 mislabels + §5 title.
  Synthesis: Pythagorean→triangle + 52/53 coverage gap closed. **Bonus (PI
  mid-sprint):** Lorentzian-propinquity rescue probe — P45 K⁺ annihilation is
  spatial-multiplier-only (Toeplitz temporal survives K⁺=ω_q) but signature-blind
  robustly (v2 Wick involution = ω_q too) → WH7 weakens-to-convention confirmed;
  `test_lorentzian_toeplitz_kplus.py` (6/6). C5/C11/C13/C14 PASS. Seed key
  `group1_B3_seed_key.json`, memo `debug/sprint_qa_group1_biteB3_memo.md`.
  **All bites remediated** (per-paper FAIL→fixed) — but see the re-cert below.
- 2026-06-19 — **RE-CERTIFICATION begun (PI), batched ~5 papers, Lorentzian-first.**
  Re-running the calibrated `/qa` gate on the remediated corpus to convert
  FAIL→remediated cycles into a trustworthy PASS. **Batch 1 (Papers 45–49 +
  synthesis) = FAIL → REMEDIATED** (v4.27.0; re-run deferred to after Batches
  2–3 per PI). 12-agent panel, calibrated (sensitivity 4/5 — citation missed an
  inline-arXiv transposition seed; specificity 5/5). **Key finding: the v4.24.0
  (B2) full-sweep was incomplete** — residual C14 "Status-note-only descope"
  zombies (p46 thm:enlarged_main, p47 thm:g2_metric, p48 T6 ×4, p49 B4′,
  synthesis l.1593+4/π), C7 (p48 ×4 "Riemannian propinquity hypertopology"), and
  C4 (bizi/leimbach in p48+p49 — B2 fixed p47 only; p49 hekkelman title; p49
  Datta "Theorem 11" → elementary + framework-verified, citation-fix). All
  remediated, all compile clean, C5/C11/C13/C14/**C15** PASS. New deterministic
  gate **C15** (`check_inline_arxiv.py`) closes the citation calibration gap
  (inline-arXiv transposition class). Seed key `group1_rc1_seed_key.json`, memo
  `debug/sprint_qa_group1_rc1_memo.md`. **Remaining: re-cert Batch 2 (42/43/44/53),
  Batch 3 (29/39/40/50/52), then re-run all to certify.**
- 2026-06-19 — **Re-cert Batch 2 (Papers 42, 43, 44, 53 + synthesis) = FAIL →
  REMEDIATED** (v4.28.0; re-run deferred to after Batch 3 per PI). 14-agent panel,
  calibrated (sensitivity 4/4, specificity 5/5); deterministic C5/C11/C13/C14/C15
  PASS. Genuine fixes: p42 ×3 result-level C7 (Paper-38 "Latrémolière
  propinquity"→state-space GH; the L5-lemma-name mention kept — matches P38's own
  naming, reconcile), p44 C8 two-way UPGRADE (P38 "qualitative-rate"→unconditional
  rate 4/π), **p53 status-note-only descope fixed in the theorem** (thm:interior
  (b),(c) → plane-only; §boundary phantom Λ^{−1.30} corrected; stale "interior
  reconstruction" ×2 → plane), p43 hekkelman2022 wrong-ID 2206.13744→2111.13865,
  p53 stempak impossible-year→ETNA 14 (2002) 223–235, connes_vs Prop 4.2
  sub-number softened, corpus-wide ack Edward→Eva-Maria Hekkelman (p42/43 + 38/39/40).
  3 reviewer over-flags reconciled (synth phantom-52/53-coverage = summary exists
  l.136–144; synth L5 = faithful to P38; vdD Prop 4.1 = GROUNDED per citation-42).
  All compile errors=0/undef=0. Seed key `group1_rc2_seed_key.json`, memo
  `debug/sprint_qa_group1_rc2_memo.md`. **Remaining: re-cert Batch 3 (29, 39, 40,
  50, 52), then re-run all to certify.**
- 2026-06-19 — **Re-cert Batch 3 (Papers 29, 39, 40, 50, 52 + synthesis) = FAIL →
  REMEDIATED** (v4.29.0; completes the three re-cert batches; certifying re-run
  deferred per PI). 15-agent panel, calibrated (sensitivity 4/4, specificity 5/5);
  C5/C11/C13/C14/C15 PASS. Genuine fixes: **p40** C8 twin (rank-uniform 4/π proof
  "now complete"→named gap) + **C4 LARGE wrong-ID** (hekkelman_mcdonald_vs2024_ucp
  = arXiv:2410.15454 is actually Bhattacharyya et al. → redirected to Rieffel/Hawkins,
  fabricated bibitem removed) + 3 bibitem metadata (bizi title, latremoliere2016
  venue, latremoliere2018 year/vol); **p39** "master theorem" zombie→sketch +
  asymmetric-sup bug (corner→interior (8,5)) + Λ^full code-comment; **p50** C7 ×7 +
  Henningson–Skenderis "inheritance"→taxonomic classification + "Weyl"→single
  (3d) Dirac + lei initials; **p52** C7 ×12 (incomplete B-series body relabel);
  **p29** 3-prism graph-ID bug (V=6/E=7/β₁=2, not C₃×P₂) + cross-paper dangling
  \ref→Paper 28; **synth** Sachdev→Safdi. 2 over-flags reconciled (synth 52/53
  coverage; vdD N/A). Coverage gaps logged, claims verified-true (p40 PRV C₃=1 not
  wired into tests/; p29 bound-crossing JSON smoke test). All 6 papers compile
  errors=0/undef=0. Seed key `group1_rc3_seed_key.json`, memo
  `debug/sprint_qa_group1_rc3_memo.md`. **All three batches remediated — the
  certifying re-run (45–49 / 42–44+53 / 29+39+40+50+52) is the next step.**
- 2026-06-21 — **CERTIFYING RE-RUN begun. Batch 1 (Papers 45–49 + synthesis) =
  FAIL → REMEDIATED** (v4.30.0). The re-run (to convert FAIL→remediated into a
  PASS) itself FAILed: calibrated (sensitivity 4/4, specificity 5/5),
  C5/C11/C13/C14/C15 PASS, and it caught rc1-residuals **plus two
  keystone-correctness issues all prior passes missed** — **p49 strict
  super-additivity was FALSE** ("strict whenever pairwise distinct"; distinct
  commuting states saturate the D_max chain at deficit 0 → restated "generically")
  and **p49's "Datta max-divergence" panel was computed with Umegaki** (relabeled
  illustrative + pointed to the real 96/96 D_max test); **p47's three surviving
  keystones had no test** (marked analytic proof-by-argument). Plus C14 zombies
  (p46 "closed convergence-theorem leg"→degeneracy; p48 §6.3 "closure
  exists"→conditional/descoped), C7 (p46 l.1289), C4 (latremoliere "dual"→404
  (2022); farsi_latremoliere2024 wrong-ID removed from p47). All 4 papers compile
  errors=0/undef=0. Seed key `group1_rr1_seed_key.json`, memo
  `debug/sprint_qa_group1_rerun1_memo.md`. **NOTE: Batch 1 re-run is FAIL→
  remediated, NOT yet a clean PASS — certification requires a clean calibrated
  pass (convergence). Re-run-2 of Batch 1 + Batches 2/3 re-runs still pending.**
- 2026-06-21 — **Re-run-2 of Batch 1 (Papers 45–49 + synthesis) = FAIL →
  REMEDIATED + dedicated C4 sweep** (v4.31.0). 2nd re-run calibrated (sensitivity
  4/4, specificity 5/5); **the v4.30.0 keystone fixes HELD** (strictness
  "generically", analytic-proof framing). Claims/code converged to SMALL residuals
  (p47 Remark 7.6 + cb-norm; p49 ×4 Umegaki-summary propagation); C4 surfaced 2 new
  LARGE (p47 hekkelman misattribution, p48 muller fabrication). Per PI, ran a
  **dedicated exhaustive citation sweep** of 45–49 (every bibitem web-verified) →
  drained the tail in one pass: ~25 bibitem corrections (wrong titles, years,
  authors; muller→GRG 54 (2022) 117, farsi→2404.00240, martinetti title; 2×
  hekkelman_mcdonald2024b dup removed). p48 MS-defs PI-caveat resolved (p49
  independently verified MS Def 2.3/3.8/4.4 exist). All 5 compile errors=0/undef=0;
  C5/C11/C13/C14/C15 PASS. Seed key `group1_rr1b_seed_key.json`, memo
  `debug/sprint_qa_group1_rerun1b_sweep_memo.md`. **C4 tail now drained + claims/
  code converged — the next clean re-run should certify Batch 1; Batches 2/3 still
  pending (shared bibliography, tail likely reduced).**
- 2026-06-21 — **Re-run of Batch 2 (Papers 42, 43, 44, 53 + synthesis) = FAIL →
  REMEDIATED + folded-in C4 sweep** (v4.32.0). Calibrated (sensitivity 4/4,
  specificity 5/5); exhaustive citation verify folded in (drained Batch-2 C4 in one
  pass). Genuine fixes: **p42 C5/K hard-prohibition** ("conjectural observation"→
  "Observation") + C8 abstract/body Lorentzian alignment; **p43 C7 ×2** (Paper-38
  "qualitative-rate propinquity"→state-space GH); **p44 N_t>1 propagation coverage
  gap** (stale docstring prop=∞→prop=2; new test passes); **p53 LARGE**
  (latremoliere2025 = compact-only 1811.10843 → pointed-proper hypertopology
  2512.03573) + stale header-comment corrected + stein_weiss kernel→means; +11
  SMALL citation-tail. 3 over-flags reconciled. All 4 compile errors=0/undef=0;
  C5/C11/C13/C14/C15 PASS. Seed key `group1_rr2_seed_key.json`, memo
  `debug/sprint_qa_group1_rerun2_memo.md`. **Batch 2 C4 drained + claims/code
  converged — next clean re-run should certify Batch 2; Batch 3 re-run pending.**
- 2026-06-21 — **Re-run of Batch 3 (Papers 29, 39, 40, 50, 52 + synthesis) = FAIL →
  REMEDIATED + folded-in C4 sweep** (v4.33.0). 15-agent panel, calibrated
  (sensitivity 4/4, specificity 5/5); exhaustive citation verify folded in. Genuine
  fixes: **p40 C4 LARGE** (vinberg1990 misattribution — the Weyl-chamber dominance
  lemma backing C₃=1 tightness cited Vinberg's shift-of-argument paper → standard
  dominance-order fact, Humphreys §13) + abstract/theorem tier caveat (rank-uniform
  4/π = named gap); **p52 C14 LARGE** (Q4 Lorentzian zombie — descoped 45/46 +
  partial 47 shown as established → descope-accurate) + C7 ×2; **p39 code C2**
  (`c3_full_pythagorean_bound_symbolic` corner-only → full-grid; interior max
  √(121/87) at (8,5); new regression test); **p50** C7 propinquity-bound→state-space
  GH + squashed-S³ misattribution (JKPS→Hama–Hosomichi–Lee) + Weyl-3d-Dirac wording;
  **p29** kotani "Thm 1.3 Ξ-completion"→proposal + bander title. Coverage gaps logged
  (p29 Obs-2 block-decomp + alon-boppana JSON-smoke; p50 wedge-KMS + CHM —
  verified-true, weak/no test; raise to PI). All 5 compile errors=0/undef=0;
  C5/C11/C13/C14/C15 PASS. Seed key `group1_rr3_seed_key.json`, memo
  `debug/sprint_qa_group1_rerun3_memo.md`. **All three batch re-runs now done; each
  needs one CLEAN re-run to convert FAIL→remediated into a certified PASS.**
- 2026-06-21 — **CERTIFYING re-run #1 of Batch 1 (Papers 45–49 + synthesis) = FAIL →
  REMEDIATED** (v4.35.0). First clean cert re-run; 16-agent panel, calibrated
  (sensitivity 5/5, specificity 5/5) — but surfaced genuine defects beyond the seeds
  (the convergence reality), so still FAIL→remediated, not yet PASS. Genuine fixes:
  **p46 C14 LARGE** (informal Main theorem re-asserted withdrawn strong-form Lorentzian
  quantum-GH convergence as a live metric distance → inline descope) + l.1424
  "convergence theorem"→"(descoped) rate-formula result"; **p48 C14** (§3 "nine
  Latrémolière axioms at theorem-grade rigor" → categorical/structural lift,
  degenerate-seminorm descope at §1.2 + §3 opener); **p47 C1** (thm:outer/thm:three_carriers
  registered proof-by-argument + three_carriers inline-tagged analytic; matrix Batch-1
  rows added — 45–49 were unregistered); **C4 ×3** (farsi_latremoliere2024 "crossed
  products" fabrication→Collapse/2404.00240 in p48+p49; kubota T.→H.; bertozzini
  title→"Modular Theory, NCG and Quantum Gravity", all web-verified). All 4 edited
  papers compile errors=0/undef=0; C5/C11/C12/C13/C14/C15 PASS. Seed key
  `group1_cert1_seed_key.json`, memo `debug/sprint_qa_group1_cert_rerun1_memo.md`.
  **A clean calibrated re-run (zero genuine defects) is still required to certify Batch 1.**
- 2026-06-21 — **CERTIFYING re-run #2 of Batch 1 (Papers 45–49 + synthesis) = FAIL →
  REMEDIATED, ONE SMALL defect** (v4.35.1). 16-agent panel, calibrated (5/5 sens with
  5 FRESH seeds, 5/5 spec); **the v4.35.0 cert1 fixes HELD (none re-flagged)**. Converged
  sharply: exactly one genuine new defect — p49 §11 C4 author misattribution
  (bousso_casini_fisher_maldacena2020 prose "Bousso, Casini, Fisher, Maldacena" →
  web-verified "Bousso, Chandrasekaran, Rath, Shahbazi-Moghaddam"; non-load-bearing,
  bibitem already correct). Fixed. p49 compiles errors=0/undef=0; C11/C13/C15 PASS.
  Seed key `group1_cert2_seed_key.json`. **Trajectory: cert1 = 4 defects (2 LARGE) →
  cert2 = 1 SMALL. One more clean re-run (zero defects) certifies Batch 1.**
- 2026-06-22 — **CERTIFYING re-run #3 of Batch 1 (Papers 45–49 + synthesis) = FAIL →
  REMEDIATED, 3 defects** (v4.36.0). 16-agent panel, calibrated (5/5 sens fresh seeds,
  5/5 spec); cert1+cert2 fixes HELD. Convergence NOT monotone (a deeper Bridge-Theorem
  read found a new zombie): **p48 C14 LARGE** (B4 convergence-transport + T2
  synthetic-compactness asserted "theorem-grade rigor" on the degenerate Krein metric →
  inline descope to abstract B4 item + thm:convergence_transport + thm:synthetic_compactness);
  **synthesis C4 MATERIAL** (wrong-ID — bousso_etal2020 = arXiv:2008.03319 resolves to
  Akers–Penington, a different paper → repointed to BCRS "Gravity dual of Connes cocycle
  flow" 2007.00230 matching p49); **p46 C8 SMALL** (Appendix L1' prop=1 → proof-sketch tag).
  All compile errors=0/undef=0; C11/C12/C15 PASS. Seed key `group1_cert3_seed_key.json`.
  **Systemic finding: the Lorentzian cluster has Status-note-only descopes on many theorems;
  a dedicated inline descope-tagging sweep of P46/P47/P48 is recommended before cert4
  (efficient vs blind re-runs).**
- 2026-06-22 — **Descope-tagging sweep DONE** (v4.37.0, PI-directed). Enumerated every
  theorem/prop in P45–49; tagged the 6 metric-level-convergence statements with
  Status-note-only descopes (p47 thm:inner; p48 thm:bridge-B4 / relaxed_triangle /
  coincidence / inframetric; p46 prop:reach_height_op). Surviving results (norm-resolvent,
  N_t=1, prop=2, B1/B3, D_max chain, B3' cardinality bound) left untagged. All compile clean;
  C5/C11/C13/C15 PASS. Memo `debug/sprint_qa_group1_descope_sweep_memo.md`. **The
  metric-convergence statements now each carry an inline descope tag — cert4 (the next
  calibrated re-run) should converge to a clean certified PASS.**

- 2026-06-22 — **cert4 FAIL→remediated** (v4.38.0). Calibrated 5/5 sens, 5/5 spec; the
  v4.37.0 descope tags HELD (no theorem-statement zombie re-surfaced). But cert4 found the
  NEXT class — descope zombies in §1-Introduction PROSE (p47) + abstract enumerations
  (p48 B4), a thm:compact_agreement substrate-vs-rate over-reach (joint rate carries the
  U(1) T/4 offset; only substrate + SU(2) factor are bit-exact), and a p49 dangling driver
  ref (→ tests/test_lorentzian_propinquity_foundation.py). All remediated; +NIT docstring
  retraction banner. che_perales venue OVER-FLAG dropped (web-confirmed DGA 103 (2026)).
  Compiles clean; C5/C11/C13/C14/C15 PASS. Memo cert_rerun1_memo §cert4. **Still not a
  clean PASS — the prose/abstract descope layer is the remaining tail; cert5 next.**

- 2026-06-22 — **prose/abstract descope sweep DONE + CLAUDE.md compaction round 4** (v4.39.0,
  PI-directed). Follow-on to cert4's finding that descope zombies survive in §1-prose/abstracts.
  Read all five P45–49 abstracts + §1 Introductions in full; 2 genuine fixes of 10 regions:
  **p48 abstract B2** (listed flatly under "theorem-grade rigor" though its own §1 calls B2
  "structurally vacuous off-orbit" → restructured B1/B3 structural vs B2 structural/vacuous vs
  B4 descoped) + **p49 §1** honest-scope bullet ("Bridge Theorem holds at theorem-grade rigor"
  → added per-leg "structural+thermal-stack legs survive; B4′ descoped"). P45/46/47 abstracts +
  P45/46/48 §1 + P49 abstract all already clean. Compiles clean; C5/C11/C13/C14/C15 PASS.
  Both the theorem-statement layer (round 1, v4.37.0) and the prose/abstract layer (round 2) are
  now swept; any remaining zombie would be in the paper *body*. CLAUDE.md 165→135 KB (round-4
  §2 compaction: granular cert/historical one-liners → docs/development_frontier_archive.md).
  **cert5 is the fresh-adversary validation.**

- 2026-06-22 — **cert5 FAIL→remediated** (v4.40.0). Calibrated 5/5 sens, 5/5 spec.
  cert5 hypothesis confirmed — remaining zombies were in the paper BODY: p47 §5.4
  (descoped inner-arrow panel "confirms the theorem at full machine precision") +
  §1-prose residuals (p48 §1.2 B1–B4, §1-roadmap substrate-vs-rate) + synthesis
  Umegaki-misattribution. LARGE: **p46 C3^op fabrication** — diagnostic confirmed the
  per-harmonic √((N−1)/(N+1)) / envelope √(1−1/n_max) is op-norm-false + mis-attributed
  to Paper 38 L3 (which is C₃=1); reworked to C₃=1 corpus-wide in p46, envelope RANGE
  kept, pure-sympy test archived, operator-system test added (13/13). All compile clean;
  gates PASS. Memo cert_rerun1 §cert5. **Theorem + abstract/§1 + body layers all now swept;
  cert6 = first real shot at a clean PASS.**

- 2026-06-22 — **rem1 targeted remediation pass** (v4.41.0, PI-directed). Re-ran ONLY the
  3 live dimensions (claims+code+synthesis+deterministic), dropping the converged-clean
  per-paper citation dimension → REMEDIATION PASS, not a cert (citation unexercised ⇒
  INCONCLUSIVE-for-cert by rule). Calibrated 4/4 sens, specificity held; all cert5 fixes
  confirmed holding. Caught TWO genuine cross-paper misses from the cert5 p46 rework:
  (1) **p46 Table 1 column** still printed √(1−1/n_max) (cert5 swept prose/eqs, missed the
  table) → fixed to C₃=1; (2) **Paper 45 carries the IDENTICAL C3^op fabrication** (eq:C3_joint_bound
  √→√(1−1/n_max) "Paper 38 L3 verbatim" + false-positive test) — never swept (cert5 was p46-only).
  code-45 caught it by running the machinery; **claims-45 had accepted the prose as sound** —
  the two-dimension save the code dimension bought. Transported the p46 fix to Paper 45 (C₃=1
  corpus-wide: eq + rem:envelope_v2 withdrawal + Q4 + production docstring + test reframe).
  All compile clean; gates PASS; tests green. **Lesson: the targeted set found a real miss the
  full cert hadn't, at 11 agents; keeping CODE in mattered. cert6 (full, w/ citation) still certifies.**

- 2026-06-22 — **rem2 second targeted remediation pass** (v4.42.0, PI-directed). Claims+code+
  synthesis+deterministic (citation dropped); remediation pass, not a cert. Calibrated 4/4 sens,
  specificity held. Validated the rem1 Paper-45 C3 transport and **found it INCOMPLETE — 3 residual
  sites**: p45 L4(d) proof (l.1131 C_3^(N)→C₃=1), p45 Table 1 (l.1570-1575 + caption → C₃=1, the
  same Table-1 class rem1 fixed in p46 but didn't check in p45), and production c_lipschitz_joint
  wiring (build() stored the √-value → now C₃=1; test → ==1.0). Plus 3 pre-existing descope-zombies
  prior sweeps missed: synthesis Headline-4 ("propinquity satisfies [→0]" lead-in, retraction not
  inline), p46 rem:limit_identification_op, p47 rem:bit_exact_g2 (metric-value remarks w/o inline
  descope). +NIT: p46 p38-test-citation overstatement softened. All compile clean; gates PASS;
  test 33/2. **The CODE dimension again caught what claims accepted (production wiring). cert6
  (full, w/ citation) still certifies; next targeted pass re-validates rem2.**

- 2026-06-22 — **rem3 CODE-DIMENSION-ONLY validation pass** (v4.42.1, PI-directed). 5 code-reviewers,
  4 vacuous-tolerance seeds (one per paper's backing test). Calibrated 4/4 sens (graded by materiality);
  specificity held (rem2 code fixes confirmed holding: production c_lipschitz_joint=1.0, √-formula
  unwired, envelope-max tight). **ZERO new MATERIAL code defects — first MATERIAL-clean dimension.**
  Drained 3 recurring test-hygiene NITs: gamma_su2_residual tautology comment (3rd flag → reframed as
  determinism check), paper_38_bound √-mis-attribution (→ per_harmonic_bound + clarified), asymptotic_rate
  docstring (C₃=1). Logged for next claims pass: p49 "joint propinquity rate" label = SU(2) factor.
  Tests pass; gates PASS. **Three laps (rem1→rem2→rem3) converged the code dimension to MATERIAL-clean,
  each lap narrower than the last. cert6 (full, w/ citation) certifies.**
