# Sprint: `/qa group1` re-certification Batch 3 (Papers 29, 39, 40, 50, 52 + synthesis) — 2026-06-19 (v4.29.0)

PI-invoked `/qa group1` re-certification, **Batch 3 = the results backbone**
(29 Ramanujan/Ihara, 39 tensor convergence, 40 universal 4/π, 50 CFT₃
F-theorem, 52 Category III) + synthesis. PI directive (standing): remediate now,
defer re-run. Completes the three re-cert batches.

## 1. Verdict: FAIL (calibrated) → remediated

15-agent panel (claims ×5 incl. synthesis, citation ×5, code ×4 [29/39/40/50]),
path-pinned to worktree `../geovac-qa-seed-group1-rc3` (removed; no seed leaked;
real corpus confirmed clean on all 4 anchors). Seed key
`debug/qa/group1_rc3_seed_key.json`.

- **Sensitivity 4/4 (all seeds caught):** S-claims-C8 (p40 named-gap rank-uniform
  proof flipped to "now complete") → claims-40; S-citation-C4 (p29 Matsuura–Ohta
  bibitem 2204.06424→2210.06424) → citation-29; S-code-C2 (p50 KPS scalar
  bit-exact assertion `10**(-40)`→`10**(40)`) → code-50; S-synthesis-C9 (Paper 50
  "reproduces bit-exactly"→"derives from a CFT calculation") → synthesis.
- **Specificity 5/5:** all controls (M1 p29 Ramanujan/π-free, M2 p39 √2-triangle,
  M3 p40 rank-1 + PRV C₃=1, M4 p50 KPS bit-exact, M5 p52 Category-III) affirmed.
- **All gating dimensions exercised + calibrated** → verdict trustworthy.
  Deterministic C5/C11/C13/C14/C15 group1 PASS. (code-42/44 spend-limit issue was
  Batch 2; this batch's code reviewers all ran.)

## 2. Genuine MATERIAL defects remediated (all PM-verified vs primary text/code)

**Paper 40:**
- **C8** (the seed's genuine twin): l.2019 "proves $c(G)=4/\pi$ rigorously at all
  ranks" → "supplies the *mechanism* … rigorous at rank 1, a proof sketch at
  rank ≥ 2 (the rank-uniform write-up is a named gap)", matching the theorem's own
  proof-sketch header.
- **C4 LARGE wrong-ID:** `hekkelman_mcdonald_vs2024_ucp` cited arXiv:2410.15454,
  which resolves to **Bhattacharyya–Duhan–Pradhan** "GH convergence of metric
  spaces of UCP maps" (no S²/Kähler content), NOT Hekkelman–McDonald–vS. The
  Kähler-S²-Berezin–Toeplitz content it was cited for is genuinely covered by
  **Rieffel (CMP 401, 2023)** + **Hawkins (CMP 215, 2000)**, both already in p40's
  bib → redirected `\cite` and removed the fabricated bibitem.
- **C4 ×3 metadata** (real-paper, B-series canonical forms): bizi_brouder_besnard
  title → "Space and time dimensions of algebras…"; latremoliere2016 venue
  Banach JMA 10 (2016)→JMPA 103 (2015) 303–351; latremoliere2018 →"*quantum* GH
  propinquity," TAMS 368 (2016) (not 370/2018); + latremoliere_metric_st_2017
  preprint-year 2017→2018 (NIT).

**Paper 39:**
- **C8 zombie:** Paper-45 bibitem (l.1476) called eq:master_theorem_v "the master
  theorem" — but §open(v) descoped it to a **sketch** (k≤2 only). Reworded to
  "the $k=2$ instance of the multi-factor form (general-$G$ is a sketch)".
- **prose bug:** "Asymmetric supremum" remark said the (10,4) grid max is ≈1.17 at
  corner (11,5); the actual sup is **≈1.179 at interior (8,5)** = √(121/87)
  (the code grid-scan is right, the remark was wrong) → corrected.
- **C7:** ack "single-factor S³ propinquity convergence" → "state-space GH".
- **Λ^full table:** the paper already documents Λ^full = C₃·(1+2√2)·max(γ)
  (closed form) and attributes only Λ^fact to the code; added a note that the
  code's `propinquity_bound_full` field is the looser conservative 2γ variant +
  fixed that field's stale "tighten to 1" docstring (withdrawn Pythagorean
  direction) to the √2 correction.

**Paper 50:**
- **C7 cluster ×7** (Paper 38 result = state-space GH, not Latrémolière
  propinquity): l.254, 1069, 1155 (comparison table), 1161, 1171, 1181, 1255.
- **C8 over-framing:** the Henningson–Skenderis "structural inheritance … both
  compute the same operator-algebraic object … reproduces structurally" softened
  to a **taxonomic M2-ring classification** ("the framework does not compute the
  HS holographic anomaly; the content is the taxonomic placement").
- **terminology:** F_D = log2/4 + 3ζ(3)/(8π²) ≈ 0.219 is the *single (3d) Dirac*
  value; "Weyl (single chirality) sector" mislabel (odd dim has no Weyl) → "single
  (3d) Dirac fermion … the doubled four-component Dirac".
- **C4 NIT:** lei_van_leuven initials J./W.D. → Y. Lei / S. van Leuven.

**Paper 52:**
- **C7 cluster ×12** (incomplete B-series relabel — §5 title was fixed, body was
  not): the at-a-glance taxonomy table (l.327) + ~10 body sentences + ack all
  named Paper 38's result "propinquity convergence/limit" → "state-space GH".
  Kept the L5-proof / propinquity-style-metric / Connes-vS / Paper-45-bibitem uses.

**Paper 29:**
- **C1 graph-ID bug:** the S³ ℓ=1 component (degree-14 zeta ⟹ E=7; n=2,3 × m∈{−1,0,1}
  = 6 vertices ⟹ β₁=2) was mis-identified as "the bipartite 3-prism C₃×P₂" — but
  C₃×P₂ has E=9/β₁=4/deg-18 and *has triangles* (not bipartite); the description
  was self-contradictory. Rewritten to the correct V=6/E=7/β₁=2 invariants.
- **dangling \ref:** eq:D_diff_closed (referenced an equation in **Paper 28**, a
  cross-paper \ref that can never resolve) → "the depth-1 identity of
  Paper 28~\cite{paper28}". p29 now compiles undef=0.

**Synthesis:** Klebanov–Pufu–**Sachdev** → **Safdi** (l.137).

## 3. Reconcile catches (reviewer over-flags — NOT acted on)
- **synthesis "Paper 52 coverage broken" (rated LARGE) = over-flag.** The 52/53
  summary exists inline (synthesis l.138–144); the "fourteen papers" core list vs
  "five further papers (48,49,50,52,53)" is a deliberate core/extra tiering, not
  an inconsistency. (Same conclusion as the Batch-2 reconcile of this finding.)
- **van den Dungen "Prop 4.1"** — N/A to this batch.

## 4. Coverage gaps logged (claim true, weak/no test — raised to PI)
- **p40 C₃(G)=1 PRV asymptotic-tightness:** no test in `tests/` (the
  `verify_dirac_triangle` driver lives in `debug/qa/_resurrected/`, passes on
  spot-check); the claim is analytically backed by the cited PRV/Kumar/Vinberg
  theorems. Wiring a `tests/` row is the named follow-up.
- **p29 bound-crossing (abstract result #1):** backed by a JSON smoke test that
  re-checks arithmetic, not a live spectral recompute. code-29 recomputed live
  (S³ n=5 dev +0.198, S⁵ N=5 +0.652, matching the JSON) — **claim verified-true**;
  the test architecture should be upgraded to a live recompute.
- **p39 Λ^full assembled column:** the code field of the same name computes the
  looser 2γ; the table is reproducible from the documented closed form (now
  clarified in prose).

## 5. Honest scope / deferred
- **Material findings all fixed.** All 6 edited papers compile errors=0/undef=0
  (p29 14pp, p39 16pp, p40 26pp, p50 17pp, p52 11pp, synth 24pp). C5/C11/C13/C14/C15
  group1 PASS. Touched code: `gh_convergence_tensor.py` comment-only.
- **Deferred (NIT, B2/rc1/rc2 precedent cosmetic sweep):** p29 Setyadi–Storm-2022
  inline attribution (unverifiable) + Kotani–Sunada "Thm 1.3" number; p39 Connes
  1994/leimbach-title; p40 polo/branson key-year + leimbach title; p50 cite-key
  years (beccaria/cardy/maldacena). Logged for the cosmetic sweep.
- **Re-run deferred** (PI): the three batches (1, 2, 3) now all remediated; the
  certifying re-run is the next step.

## 6. Next
Re-run all three batches (45–49, 42/43/44/53, 29/39/40/50/52) on the remediated
corpus to convert the FAIL→remediated cycles into a trustworthy PASS → certify
group1.
