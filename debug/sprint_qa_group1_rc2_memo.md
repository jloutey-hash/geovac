# Sprint: `/qa group1` re-certification Batch 2 (Papers 42, 43, 44, 53 + synthesis) — 2026-06-19 (v4.28.0)

PI-invoked `/qa group1` re-certification, Lorentzian-first batching, **Batch 2 =
the Lorentzian foundations + disk** (42 modular four-witness, 43 finite-cutoff
Wick rotation, 44 BBB Krein operator system, 53 disk-with-cone) + synthesis
coverage. PI directive (standing from rc1): **remediate now, defer re-run.**

## 1. Verdict: FAIL (calibrated) → remediated

14-agent panel (claims ×5 incl. synthesis, citation ×4, code ×3 [p42/p44/p53]),
path-pinned to worktree `../geovac-qa-seed-group1-rc2` (removed; no seed leaked;
real corpus confirmed clean on all 4 anchors). Seed key
`debug/qa/group1_rc2_seed_key.json`. (code-42/code-44 hit a mid-run spend limit
and returned empty; re-dispatched after the limit cleared — both SOUND.)

- **Sensitivity 4/4 (all seeds caught):** S-claims-C14 (p43 descoped continuum
  Lorentzian propinquity flipped to "is now established") → claims-43;
  S-citation-bibitem-C4 (p44 van den Dungen bibitem 1505.01939→1505.04821) →
  citation-44; S-code-C2 (p53 plane-convergence `assert p > 0.3` gutted to
  `> -1.0e9`) → code-53; S-synthesis-C9 (synthesis Λ-inheritance "descoped"→
  "established") → synthesis.
- **Specificity 5/5:** all controls (M1 p42 four-witness, M2 p43 finite-cutoff
  Wick, M3 p53 plane Bochner–Riesz, M4 connes_rovelli/camporesi_higuchi cites,
  M5 v4.26.1 Lorentzian-closure cross-refs) affirmed sound — no false positive.
- **All gating dimensions exercised + calibrated** → verdict trustworthy.
  Deterministic C5/C11/C13/C14/**C15** group1 PASS.

## 2. Genuine MATERIAL defects remediated (all PM-verified vs primary text/code)

**C7 (Paper 38 = van Suijlekom state-space GH, not Latrémolière propinquity):**
- p42 ×3 **result-level** mislabels — "converge … in the Latrémolière quantum
  GH propinquity … established in Paper 38" (l.201, l.468, l.2213) →
  "in van~Suijlekom's state-space Gromov–Hausdorff distance". Reconcile narrowed
  the reviewer's count: **l.486 is CORRECT** (names Paper 38's *L5 lemma*
  "Latrémolière propinquity assembly" — Paper 38 itself names it so; the C7 rule
  bites the *result/distance*, not the proof *machinery*).
- synthesis Headline-2 **header** "Riemannian propinquity convergence" →
  "Riemannian state-space Gromov–Hausdorff convergence" (its own body l.221
  already said state-space GH). Synthesis **Lemma L5** left intact — it mirrors
  Paper 38's own L5/Λ structure (faithful, not a zombie).

**C8 two-way UPGRADE (p44):** conclusion l.1481 *understated* keystone Paper 38
as "qualitative-rate state-space GH" → "state-space GH (unconditional, rate
constant 4/π)", matching p44's own §9 l.1282 and Paper 38's UNCONDITIONAL status.

**C8/C14 status-note-only descope (p53 — the heaviest finding):** the disk
`thm:interior` (b),(c) were still asserted "for all Λ" in the **theorem body**,
repaired only by the following Remark (the run-#4 / B2 failure pattern). Fixed
*in the theorem statement*: lead-in now scopes (a),(d) to the finite disk and
(b),(c) to the boundaryless plane only ("withdrawn on the finite disk"); item
headers marked "(plane only)". Proof of (c) (l.421) re-scoped to plane
(γ_Λ→0 on plane; Λ^{+0.07} no-decay on disk). §boundary phantom **Λ^{−1.30}**
(l.507, self-contradicting its own rem:numerics) → "even the interior
coefficient does not decay (Λ^{+0.07}); convergence recovered only on the
plane". Two stale "interior reconstruction" (l.482 cor proof, l.492 remark) →
"plane reconstruction".

**C4 citations:**
- p43 `hekkelman2022` bibitem WRONG-ID `arXiv:2206.13744` (resolves to a
  Kerr–Melvin black-hole paper) → `arXiv:2111.13865`, "Truncated geometry on the
  circle", LMP 112 (2022) 20 (matches p42/p44's correct form; orphan bibitem but
  canonical failure mode → fixed).
- p53 `stempak1989` bibitem impossible year (ETNA began 1993) → added
  ETNA **14 (2002) 223–235**.
- connes_vs `[Prop.~4.2]` sub-number (UNVERIFIABLE from source, flagged by 3
  reviewers) softened to plain `\cite{connes_vs2021}` in p42 (l.197) + p53
  (l.323); the prop=2 result itself is framework-verified (Paper 32 §III).

**Corpus-wide factual name fix (deferral-is-churn):** acknowledgments "Edward
Hekkelman" → "Eva-Maria Hekkelman" (the real arXiv:2412.00628 author) in p42,
p43, **and** p38/p39/p40 (out of batch but fixed at discovery).

## 3. Reconcile catches (reviewer over-flags — NOT acted on)

- **synthesis Finding 1 "phantom Paper 52/53 coverage" (rated LARGE) = OVER-FLAG.**
  The abstract's "Papers 52–53 are summarized below" IS satisfied: the summary
  sits immediately below at l.136–144 (P52 Category III; P53 disk-obstructed /
  plane-genuine). 52/53 are prose-named (no `\cite`) → no bibitem required → no
  compile issue. B3's coverage-gap fix is working; "fixing" it would be churn.
- **van den Dungen "Prop 4.1 content" UNVERIFIABLE** (citation-43/44): citation-42
  verified the i^t Wick recipe GROUNDED *verbatim* against arXiv:1505.01939, and
  B1 confirmed it — no action (the two reviewers' PDFs wouldn't decode).
- **code dimension** (p42/p44): all load-bearing claims map to genuine,
  non-tautological, passing tests (four-witness β=2π closure reproduced
  bit-exactly from the real CH integer-doubled spectrum; BBB envelope-dependent
  prop_ach=2/prop_full=∞ verified non-vacuous). SOUND.

## 4. Honest scope / deferred
- **Material findings all fixed.** All 5 structural-edit papers (+38/39/40 ack)
  compile errors=0 / undef=0 (p42 28pp, p43 25pp, p44 18pp, p53 9pp, synth 24pp).
  C5/C11/C13/C14/C15 group1 PASS.
- **Deferred (NIT, B2/rc1 precedent):** cosmetic citation slips with correct
  arXiv IDs — Latrémolière 2017/2018 display-year labels, strohmaier year,
  toyota "M."→"R.", bykov/leimbach/colzani title paraphrases. Plus a flagged
  re-verify of the BBB (m,n)=(4,6) Table-1 sign attribution against the BBB PDF
  (the period closures are DL-independent, so not load-bearing). Logged for the
  cosmetic title/year sweep.
- **Re-run deferred** (PI direction): Batch 2 re-run after Batch 3.

## 5. Next
Batch 3 (Papers 29, 39, 40, 50, 52 — results backbone), then re-run all three
batches to certify.
