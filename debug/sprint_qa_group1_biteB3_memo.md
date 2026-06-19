# Sprint: `/qa group1` Bite B sub-bite 3 (Papers 39, 52, 53) — 2026-06-18 (v4.25.0)

The last group1 sub-bite. PI-invoked `/qa group1 Bite B sub-bite 3`. Verdict
**FAIL (calibrated)**; remediated via two PI-directed **rescue-first** efforts
(both headlines kept where the object was real, downgraded where it wasn't) plus
a spontaneous Lorentzian-propinquity rescue probe the PI requested mid-sprint.

## 1. Verdict + calibration

9-agent panel (claims/citation/code per paper + synthesis), path-pinned to
worktree `../geovac-qa-seed-group1-B3` (removed; no seed leaked). **Calibrated:
sensitivity 5/5** (S-claims-C7→claims-39, S-claims-C8→claims-52,
S-citation-C4→citation-52, S-synthesis-C9→synthesis, S-code-C2→code-39),
**specificity held** (honest controls affirmed SOUND). Seed key
`debug/qa/group1_B3_seed_key.json`.

Two SERIOUS genuine findings (both PM-verified against code/drivers per §9
reconcile rule before acting):

## 2. p39 — false C₃<1 keystone (rescue attempted, downgraded)

**code-39 F1.** The boxed operator-norm **Pythagorean identity**
`‖[D,M_f⊗M_g]‖² = ‖[Dₐ,M_f]‖²‖M_g‖² + ‖M_f‖²‖[D_b,M_g]‖²` (the basis of
`C₃<1`) is **false**, and its footnote ("verified via direct sympy
diagonalisation … machine-precision identity") is **fabricated** — the cited
`c3_full_pythagorean_bound` does pure per-irrep arithmetic, no diagonalization.

**PM verification (rescue attempt, PI: "try to rescue first").** Built the real
graded joint commutator from the Camporesi–Higuchi Dirac + graded multipliers
(`full_dirac_operator_system`): over every single-harmonic pair at
(2,2),(2,3), **max ‖A+B‖²/(‖A‖²+‖B‖²) = 2.0000 exactly** (Pythagorean maximally
violated — the two Leibniz terms add constructively) and **max ‖A+B‖/(‖A‖+‖B‖) =
1.0000 exactly** (triangle TIGHT). So C₃<1 is unrescuable; the rigorous constant
is the tight triangle **C₃ = √(((N_a−1)+(N_b−1))²/(N_a²+N_b²−2)) ↗ √2** (= 1 at
(2,2), > 1 beyond (3,3)). **The convergence theorem survives** (√2 is finite →
propinquity → 0); only the cosmetic "<1" sharpness died.

**Remediation:** corrected Lemma L3-T (Pythagorean→triangle + Remark
"Withdrawn Pythagorean refinement"), eq:per_irrep_ratio_T, eq:C3_2_full, the
two remarks, abstract, intro, main theorem, §1 ingredient list, conclusion,
numerical table (C₃ + Λ_full columns recomputed), the §10 k-fold + master
theorems (F2: "k-INDEPENDENT via Pythagorean →1⁻" corrected to "√k via
triangle", downgraded from CLOSED to sketch). Code:
`c3_full_triangle_bound` (new, correct formula) + deprecated alias
`c3_full_pythagorean_bound`; symbolic twin fixed. Tests:
`test_paper39_triangle_tight.py` (NEW — verifies Pyth ratio→2, triangle→1 on
real harmonics; the diagonalization the footnote falsely claimed) +
`TestTriangleC3Bound` rewrite of `test_gh_convergence_tensor.py`.
Citations: hekkelman2022 (2206.13744→2111.13865 "Truncated geometry on the
circle"), hekkelman_mcdonald2024 (fabricated 2403.18619 → repointed to real
2024b + bibitem removed), avery_wen_avery1986 (→ Wen & Avery JMP 26 (1985) 396),
latremoliere2016/2018 orphan bibitems (→ dual JMPA 103 (2015) / quantum TAMS 368
(2016), web-verified).

## 3. p53 — phantom disk numerics (rescue SUCCEEDED via the plane)

**code-53.** The cited backing tests are S³ (wrong carrier); the disk
Markov–Cesàro positivity (s≥2) and interior rate Λ^{−1.30} are claimed "verified
numerically" but are **refuted/phantom**.

**PM verification.** Ran all three drivers:\ `b3_markov_berezin.py` (the paper's
disk construction) gives **min eig = −0.13** (positivity FAIL) and rate
**Λ^{+0.07}** (no decay, not Λ^{−1.30}); `b3_disk_heat_berezin.py` also
Λ^{+0.10} FAIL; but `b3_plane_bochner_riesz.py` (the boundaryless plane =
de-compactified cone, already the paper's stated cure) **converges**:
positivity at Bochner–Riesz s≥½, rate **Λ^{−0.6 to −0.9}**. So the genuine
convergent carrier is the plane (rescue succeeds); the finite disk is obstructed
by its Dirichlet boundary.

**Remediation:** abstract + §results + honest-scope corrected (disk:\ unitality
only; positivity + rate fail = boundary obstruction; plane:\ genuine);
thm:interior got a correction Remark (parts (b),(c) withdrawn on the finite
disk); rem:numerics replaced with actual driver output; cor:apex_backbone routed
through the plane; Q2 (phantom rate) → plane rate; new genuine carrier test
`test_paper53_plane_bochner_riesz.py`; C7 propinquity→state-space-GH relabels;
Paper 45/46 lineage zombie fixed; latremoliere2025 grab-bag → verified
Adv. Math. 404 (2022).

## 4. p52 + synthesis

p52:\ 4 C7 mislabels (Paper 38 "Latrémolière propinquity" → state-space GH:\
abstract, §1.4, §examples, §6 table) + §5 title. Synthesis:\ Paper-39
Pythagorean→triangle + the k-fold/master-theorem "k-INDEPENDENT" → √k-sketch
corrections; **52/53 coverage gap closed** (added a summary paragraph + the 50/52/53
one-liners).

## 5. Lorentzian-propinquity rescue probe (PI mid-sprint question)

PI: "can we rescue the Lorentzian propinquity?" → "run the probe now" +
"scope v2 now". Separate memo `debug/sprint_lorentzian_toeplitz_kplus_probe_memo.md`;
result: **rescuable from annihilation** (P45 K⁺ annihilation is
spatial-multiplier-only — the Toeplitz temporal algebra survives K⁺ = ω_q), but
**signature-blind robustly** (v2 genuine Wick involution J_t=sign(D_t) leaves
the surviving seminorm = ω_q bit-exact) → confirms WH7 "weakens-to-convention".
Recorded:\ P45 abstract + §open Q1 update, CLAUDE.md §1.7 WH7, frozen test
`test_lorentzian_toeplitz_kplus.py` (6/6).

## 6. Honest scope
- **Theorem-grade kept:** p39 tensor state-space GH convergence (now C₃≤√2,
  triangle-tight, backed); p53 plane Bochner–Riesz convergence (driver +
  test-backed); p50 F-theorem (Bite-A); the K⁺/Toeplitz v1+v2 facts (bit-exact).
- **Downgraded (honest):** p39 "C₃<1" → "C₃≤√2"; p39 k-fold/master theorems
  CLOSED → sketch (√k, not k-independent); p53 finite-disk convergence →
  withdrawn (plane only).
- **Named open follow-ons:** p39 k>2 / general-G master theorem (sketch);
  p53 sharp closed-form plane rate; genuine indefinite-signature Lorentzian
  metric (needs a state/modular-level device, not a Lipschitz seminorm).
- All 5 papers compile errors=0/undefined=0; C5/C11/C13/C14 group1 PASS.

## 7. Closes group1
With bites A (29/40/50), B1 (42/43/44), B2 (45–49), and B3 (39/52/53) all
FAIL→remediated, the **group1 branch QA is complete**.
