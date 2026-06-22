# Sprint: `/qa group1` certifying re-run of Batch 2 (rr2) + folded-in C4 sweep — 2026-06-21 (v4.32.0)

Certifying re-run of Batch 2 (Papers 42, 43, 44, 53 + synthesis), with an
**exhaustive citation verify folded in** (citation reviewers checked every
bibitem) so one run both calibrates the verdict and drains the Batch-2 C4 tail.

## 1. rr2 verdict: FAIL (calibrated) → remediated

12-agent panel (claims ×4, citation ×4 exhaustive, code ×3, synthesis ×1),
worktree removed (no seed leaked). Seed key `debug/qa/group1_rr2_seed_key.json`.
- **Sensitivity 4/4:** S-claims-C14 (p43 continuum propinquity "remains open"→"is
  now established")→claims-43; S-citation-C4 (p42 Camporesi–Higuchi JGP 20 (1996)→
  24 (1998))→citation-42; S-code-C2 (p42 four-witness period-closure `<1e-12`→
  `<1e12`)→code-42; S-synthesis-C9 (p53 "recovered only on the plane"→"recovered
  on the finite disk")→synthesis.
- **Specificity 5/5:** controls affirmed; the v4.28.0 rc2 fixes held.

## 2. Genuine MATERIAL defects remediated

- **p42 — C5/K hard-prohibition:** "Paper 2's status as a **conjectural**
  observation" → "as an **Observation**" (K=π(B+F−Δ)/Paper 2 is an Observation,
  never conjectural — §13.5).
- **p42 — C8 abstract/body alignment:** abstract said the Lorentzian (3,1)
  extension "remains the named multi-month follow-up", but §9 delivers a
  finite-cutoff Krein four-witness theorem → reworded: the finite-cutoff Krein
  closure is given in §9, the *continuum* (Lorentzian-propinquity-rate) closure is
  the follow-up.
- **p43 — C7 ×2:** "Paper 38's qualitative-rate **propinquity** convergence" →
  "**state-space GH** convergence" (l.1822, l.1935); the Lorentzian-propinquity
  *target* mentions (the open thing) correctly kept as "propinquity".
- **p44 — N_t>1 propagation coverage gap:** the headline prop_ach=2 "for all N_t"
  was tested only at N_t=1, and the code docstring asserted the OPPOSITE (prop=∞
  at N_t>1). Fixed the stale docstring (achievable envelope = dim_Weyl²·N_t,
  prop=2 at all N_t — matching the code's own `achievable_envelope_dim` and the
  paper) + added `test_propagation_number_nmax_2_Nt_3_achievable_envelope`
  (prop=2, dims 42→192; passes).
- **p53 — LARGE citation misattribution:** `latremoliere2025` (co-cited 5× for the
  "pointed/proper locally-compact propinquity machinery" of the non-compact plane
  assembly) pointed to arXiv:1811.10843 — a **compact/unital-only** paper.
  Repointed to the actual pointed-proper-locally-compact work, the hypertopology
  paper arXiv:2512.03573 (key year "2025" now correct; co-citation {survey 2016 +
  hypertopology 2025} both genuinely support the assembly).
- **p53 — stale header-comment** (non-rendered, but harvestable): the FIRST-DRAFT
  block still asserted disk positivity / Λ^{−1.30} / disk-as-carrier → replaced
  with the CORRECTED status (disk obstructed; plane is the carrier).
- **p53 — stein_weiss overstatement:** "the Bochner–Riesz **kernel** is positive
  above the critical index" → "the Bochner–Riesz **means** are positivity-preserving
  … confirmed numerically, min g ≈ −2×10⁻⁴" (the kernel is sign-changing).

**C4 tail (11 SMALL, drained via the exhaustive sweep):** p42 strohmaier label
2006→2000 + latremoliere2018 label →2016; p43 leimbach "for tori" + connes 1994 +
franco label 2013; p44 bizi "applications" + bykov "…unbounded case" title +
latremoliere_metric_st_2017 preprint 2018 + latremoliere2018 label →2016 + leimbach
"for tori" + toyota R.

## 3. Reconcile catches (over-flags / not-acted)
- synthesis "Papers 52–53 have no body section / promise unfulfilled" — **over-flag**
  (same as Batch-2/3 reconcile): the 52/53 summary is the inline abstract/intro
  block (synth l.136–144); "summarized below" = those sentences, no `\cite` → no
  bibitem needed. Not acted.
- p43 l.401–402 "an analog of Paper~38 at signature (3,1)" — NOT a C7 violation
  (names the *Lorentzian* propinquity target, correct), left as-is.
- p42 six-witness "collapse" tested via residual-equality (code-42 NIT) — the Δ
  bit-identity is true (reviewer verified by hand); coverage NIT, not acted.

## 4. Honest scope / deferred
- All 4 papers compile errors=0/undef=0 (p42 28pp, p43 25pp, p44 18pp, p53 9pp);
  new p44 test passes; touched code (`operator_system_lorentzian.py` docstring)
  comment-only. C5/C11/C13/C14/C15 group1 PASS. Worktree removed; no seed leaked.
- Deferred NIT (cosmetic): stale cite-KEY names whose displayed bibitems are
  correct (p42 `ucp_maps_2024`=Rieffel 2004, `zhu_casini2020`=Zhang et al.; p53
  `stempak` key, `cheeger1983` orphan-but-correct) — reader-invisible.
- **Certification status:** Batch 2's C4 tail is now drained (exhaustive verify)
  and claims/code converged → the next clean re-run should certify Batch 2.

## 5. Next
Batch 3 re-run (29, 39, 40, 50, 52), then clean re-runs to certify.
