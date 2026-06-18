# Sprint: `/qa group1` Bite A (Papers 40, 50 + synthesis) + remediation — 2026-06-17 (v4.22.0)

PI-invoked `/qa group1`, second bite (after Paper 29 single-paper, v4.21.3).
Scope: the active math.OA backbone **Papers 40 (universal 4/π rate), 50 (CFT₃
F-theorem)** + the **group1 synthesis**. C4 (citation apparatus) high-priority
for this branch.

## 1. Verdict: FAIL (calibrated) → remediated

Panel **calibrated**: sensitivity **4/4** (S-claims C8, S-citation C4, S-code
C2, S-synthesis C9 all caught), specificity **5/5** (controls held). The
calibrated panel surfaced substantial **genuine** defects, vindicating the C4
high-priority flag (two fabricated/wrong-ID citations beyond the seed).

Seed key: `debug/qa/group1_A_seed_key.json`. Worktree
`../geovac-qa-seed-group1` removed; no seed leaked (verified).

## 2. Genuine defects found + fixed

**LARGE — citations (all verified vs arXiv by PM):**
- **paper_50 `beccaria_tseytlin2017`** — arXiv:1702.02325 = Smith, Goldfeld
  conjecture (math.NT); cited title ("N=4 conformal supergravity dilatation
  operator") matches no B–T paper. → **Beccaria, Bekaert, Tseytlin, "Partition
  function of free conformal higher spin theory," JHEP 08 (2014) 113,
  arXiv:1406.3542** (genuinely covers free-field/Maxwell partition functions).
- **paper_50 `hartman…2019`** — arXiv:1902.10893 = Caputa–Datta–Shyam. Title +
  venue (JHEP 03 (2019) 004) correct → arXiv **1807.11401** (verified).
- **paper_40 `hekkelman_mcdonald2024`** ("Spectral truncations of T^d") —
  **fabricated**: no such paper; arXiv:2403.18619 = an OpenMP shortest-path
  paper (cs.DC). The flat-tori work is Leimbach–vS (already cited). →
  re-attributed to Leimbach–vS, bibitem removed.

**LARGE — code backing:**
- **paper_50 F-theorem was FALSE-POSITIVE-backed**: the only "coverage" was the
  literal cross-product `F_D−2F_s` in `test_paper55` (typed-in F-values, never
  computed). → new `tests/test_paper50_f_theorem.py` (8 tests) +
  `qed_two_loop.py::{scalar_F_theorem, dirac_F_theorem,
  conformal_scalar_zeta_s3}` — computes ζ′(0) from the framework S³ spectrum →
  KPS bit-exact (1e-82); Dirac via the shipped T9 Hurwitz form, scalar via the
  binomial continuation; two independent methods + spectrum-tie cross-checks.
- **paper_40 rank≥2 universality had NO genuine test** (only SU(2) = Paper 38).
  → see §3.

**LARGE — synthesis zombies:**
- l.1797 "extended in turn to the Lorentzian K⁺-weak-form by Paper 45" (Paper 45
  is the degeneracy theorem, not an extension) → corrected to descope framing.
- §9.1 (Paper 48) re-asserted "G2 metric-level **closure**" + theorem-grade
  bridge inside a descoped subsection → conditioned to "proposed route".
- paper_50 theorem + abstract asserted "Latrémolière **propinquity**
  convergence … via Paper 45" → state-space GH (Paper 38 descope); 5 bare
  "paper38, paper45" cites pinned to the surviving product-carrier rate.

**SMALL:** paper_50 "predictive classification of CFT physics" → "reproduction,
not prediction"; paper_40 Paper-39 cross-ref title synced (propinquity →
state-space GH); paper_40 over-rigor calibration (see §3).

## 3. Paper 40 rank≥2 universality — the load-bearing finding

The "4/π universal across compact Lie groups" headline at rank ≥ 2 was claimed
as a "rigorous theorem at the rigor level of Paper 38 Appendix A," with detailed
bookkeeping deferred to `l2_universal_rate_memo` — a **pruned** debug memo.

**De-risk** (`debug/qa/su3_rate_prototype.py`): a from-scratch Weyl-integration
reproduces SU(2) bit-exactly but lands SU(3) off ~2.9× under naive conventions;
tuning a factor to force 4/π would be circular.

**Resurrection** (PI's question — "can't we just resurrect it?"): YES. Recovered
from git history (`56d29dd^`) into `debug/qa/_resurrected/`:
`l2_universal_rate_verify.py`, `sp2_g2_rate_constant.py` (+ `dirac_triangle_
extended_verify.py` dep), proof memo, data JSONs. The genuine computation
**runs + is reproducible** (Haar = 1.0). BUT it reveals the rate extraction is
**fit-sensitive**: the leading Stein–Weiss constant scatters with fit order +
min-L cut (Sp(2) 3-param across cuts: 3.12 / 1.67 / 0.73 / 0.13 / −0.64; G2:
1.97 / 4.27 / 2.06). The clean table values (SU(3) 1.243, Sp(2) 1.087, G2 1.177)
come from a specific 2-param fit + the generic→canonical Λ rescaling
(2.175/2 = 1.087; 2.88/√6 = 1.176).

**Resolution (PI chose calibrate; resurrection confirmed it):**
- Genuine test `tests/test_paper40_universal_rate.py`: SU(2) → 4/π (rigorous
  rank-1 anchor); rank-2 machinery Haar = 1.0; **G2 A-over-B discrimination**
  (extracted c ≪ Reading-B 24/π² ≈ 2.43, closer to A = 4/π) — the robust
  content, NOT the exact value.
- Paper 40 prose calibrated: rank-1 rigorous; rank≥2 mechanism = proof sketch,
  numerically confirmed, **extraction fit-sensitive**, robust claim = A-over-B,
  full rank-uniform proof = **named gap**. Removed all "rigorous theorem at
  Paper 38 rigor at all ranks" + the pruned-memo "detailed bookkeeping" cites
  (orphan bibitem deleted).

## 4. Files

- **New:** `tests/test_paper50_f_theorem.py`, `tests/test_paper40_universal_rate.py`;
  `debug/qa/_resurrected/*` (recovered drivers + README), `debug/qa/su3_rate_prototype.py`,
  `debug/qa/validate_recovered_rate.py`, `debug/qa/group1_A_seed_key.json`.
- **Code:** `geovac/qed_two_loop.py` (+3 F-theorem functions).
- **Papers:** paper_40, paper_50, group1 synthesis (.tex + .pdf); all compile 0
  errors / 0 undefined cites.
- **Docs:** `docs/claim_test_matrix.md` (group1 Bite-A section); `docs/qa/group1.done.md` (log).

## 5. Honest scope
- **Theorem-grade:** Paper 50 F-theorem (F_s, F_D bit-exact from framework
  spectrum). Paper 40 SU(2) rank-1 = 4/π.
- **Numerical / calibrated:** Paper 40 rank≥2 universality = A-over-B robust,
  exact-4/π fit-sensitive (named gap for the rank-uniform proof).
- **Pre-existing (out of scope, noted):** `test_qed_two_loop.py` has a global
  `mpmath.mp.dps` test-isolation fragility (7 failures only when co-run with
  `test_two_loop_self_energy.py`; passes alone — reproduces with my changes
  stashed out). Hygiene follow-on, not a Bite-A defect.
- **Deferred (PI-acknowledged):** the 443 debug/-citation corpus sweep (paper_40
  still cites `sp2_g2_rate_memo`/`su4_rate_memo` — advisory debug-refs).

## 7. Follow-on (v4.22.1): resurrect-pruned rule + walked-past S⁵ gap

PI direction: make "resurrect pruned artifacts for context before concluding"
a standing rule (`memory/feedback_resurrect_pruned_artifacts.md`), and look back
for a walked-past case.

**Walked-past case found + closed:** Paper 50's **S⁵ F-theorem** (Thm scalar_S5
/ dirac_S5). The Bite-A S³ test noted the `ads_track_a_s5` backing was pruned and
left S⁵ untested. Resurrecting it surfaced a **factor-4 multiplicity bug in the
original scalar memo** (prefactor 1/3 → degeneracy 4 at n=0; standard S⁵
harmonic count is 1). **The paper value is correct** (1/12, deg 1,6,20) —
recomputed bit-exact: F_s^{S⁵} = −log2/128 − ζ(3)/128π² + 15ζ(5)/256π⁴;
D′(0)^{S⁵} verified |diff|~1e-52 (Dirac memo + paper agree). New
`qed_two_loop.{scalar,dirac}_F_theorem_s5` + S⁵ block in
`test_paper50_f_theorem.py` (11 tests total). No paper edit (values already
correct). S⁷ negative-finding memo also recovered → S⁷ catalogue-row NIT
available for the next pass.

**Lesson:** the rule's first application immediately caught a real (latent) bug
+ confirmed the paper over a buggy pruned source — exactly the "back /
calibrate-with-evidence / build-fresh-from-original" payoff.

## 6. Next
Bite B (PI: smaller sub-bites) — Lorentzian cluster 39, 42–49, 52, 53
(C14-heavy). Suggested: 42–44, then 45–49, then 39/52/53.
