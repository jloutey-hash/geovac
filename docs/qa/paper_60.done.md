# Paper 60 (Generalized-Sturmian Secular Equation) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-60 scope + deltas + watch-notes.

> **STATUS: CERTIFIED ✅ 2026-08-18 — final FULL certifying run = PASS.** Over the post-engine-sprint
> text: panel **FULLY CALIBRATED** (sensitivity **6/6** fresh seeds — code FC1 gutted-band + FC2
> gutted-tolerance; citation FT1 gslw2019-wrong-ID + FT2 lowchuang2019-vol/year; claims FP1
> beats-DF/THC + FP2 novelty-tier-flip), **specificity clean** (every reviewer independently recomputed
> the paper's numbers and confirmed they reproduce — no wrong-value false-flags), **zero verified
> cert-blocking MATERIAL** across code / claims / citations + completeness-critic. The code reviewer
> re-ran the full tracked suite (**28 passed**) and confirmed the engine sprint's backing is genuine
> (1.19/3.33, K^0.842, cond 4→3672.6, n_orb^2.22, N^1.97 all reproduce exactly; engine `validate()`
> non-circular vs the symbolic route). Only fix-on-sight NITs remained — one loose L²-divergence
> tolerance (**tightened to a band around 3673 during the run**), plus acknowledged texture NITs
> (collective-scale RSS identity; un-anchored inline resource numbers; 3 Avery-source UNVERIFIABLE
> book/thesis citations that are inaccessible primary sources, not misattributions). Deterministic
> C10–C18 green; paper compiles clean (7 pp). Final seed key `debug/qa/paper_60_final_seed_key.json`;
> worktree removed, no seed leaked. **HONEST CEILING:** "CERTIFIED" = survived the calibrated detectors
> for the seeded defect classes + the pre-registered criteria; the Avery unpublished-thesis/2006-book
> content attributions remain UNVERIFIABLE (primary source not web-accessible), and the exact isoenergetic
> secular-matrix (−p_κ𝟙) form is primary-source-gated to Avery's formalism — both honest standing limits,
> not coverage gaps in GeoVac's own machinery (now fully tracked + tested).
>
> *(superseded)* **PASS-ELIGIBLE 2026-08-18 — weak-backing closed by the engine sprint; awaiting final PI stamp.**
> Arc: first-cert FULL = FAIL (trustworthy) → remediation → **CLEAN-DELTA** → full certifying run
> (panel calibrated **6/6** seeds, specificity clean, **zero wrong §C8 headline value**) → the full run's
> only genuine findings were *weak-backing/coverage* (exact exponents + large-R/L²-divergence/eq:secular
> tracked only by untracked drivers) → **ENGINE SPRINT (v4.96.0)** promoted the machinery to 4 tracked
> `geovac/sturmian_*.py` modules with faithful ≥4-pt tests, closing every one of those gaps. The
> substantive weak-backing that held the PASS is resolved; the residual is fix-on-sight NITs (all applied).
> A final formal certifying pass over the strengthened text is the only remaining step for the CERTIFIED
> stamp. Full-cert seed key `debug/qa/paper_60_full2_seed_key.json`.
>
> *(superseded)* **FIRST-CERT FULL run = FAIL (trustworthy) 2026-08-18 — remediated.**
> Panel **FULLY CALIBRATED**: sensitivity **6/6** planted seeds caught (2 code SC1/SC2, 2
> citation ST1/ST2, 2 claims SP1/SP2), specificity **5/5** (0 known-good controls
> false-flagged), across all three exercised dimensions (code / claims / citation; C9
> out-of-scope). Deterministic C5/C10–C18 all green (C17 vacuous for Paper 60 — no families).
> **Verified non-seed MATERIAL:** W1 — abstract `K^0.78` vs canonical `eq:sublinear K^0.84`
> (the claims enumeration cross-check caught it independently; touches C8#4 headline).
> **SMALL genuine:** `baek2023` author "J."→"U." (Unpil Baek, verified vs arXiv:2205.09039).
> **Load-bearing coverage gaps:** `λ~n_orb^2.2` (C8#10) + `Q^3.33/Q^1.19` (C8#1) have NO
> backing test; `K^0.84` (C8#4) backed only qualitatively. W2 (abstract gerade ratio) = NIT
> (conservative under-statement). W7 (`rajchel2025`) confirmed already-clean (v4.95.0 fix).
> Seed key `debug/qa/paper_60_seed_key.json`; worktree removed, no seed leaked. **Path to
> cert:** PM fixes mechanical (W1/W2/baek2023/tolerances) + adds C17 families; PI directs the
> coverage-gap decisions (add tests vs accept) → delta-verification re-run → FULL certifying
> run. Third single-paper `/qa` target (after Paper 58 CERTIFIED 2026-08-15/17, Paper 59
> CERTIFIED 2026-08-17), capstone of the Avery arc (58→59→60). Paper 60 (v4.89–v4.95, 2026-08-18) postdates the
> 2026-06-28 group2 cert and has **no synthesis footprint** (grep of the group2 synthesis
> for paper-60 / isoenergetic / secular / metric-free = empty), so C9 is OUT OF SCOPE this
> run (as for Paper 59). Backing suite `tests/test_paper60_sturmian.py` = **17/17 green with
> `--slow`** (9 pass / 8 slow-skipped by default — the slow ones carry the load-bearing
> headlines, so the code dimension MUST run `--slow`). **Two genuine abstract↔body headline-
> number tensions were found while drafting these criteria (W1, W2 below) — flagged for PI
> decision before freeze.**

**Scope (single-paper):**
- **Paper 60** — `papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`
  ("The Generalized-Sturmian Secular Equation Block-Encodes Without a Metric:
  A Quantum-Algorithm Analysis, Metric-Free for Atoms and a Conditioning Frontier
  for Molecules").
- **C9 (synthesis): OUT OF SCOPE** — Paper 60 has no group2-synthesis footprint yet
  (verified by grep). If the PI wants a synthesis promotion written first, that is a
  separate task; it is NOT part of this cert.
- **Out of scope:** the other group2 papers (unchanged since 2026-06-28 cert). Papers
  58/59 (`loutey_paper58`/`loutey_paper59`) are CERTIFIED; in scope only as cross-refs
  (C7). Trunk papers (0/1/7/14/18) canonical; in scope only where Paper 60 restates them.

**Deterministic `--gate`:** `group2` (Paper 60 lives under `group2_quantum_chemistry/`).
**C17 note:** the headline-number registry currently has **NO Paper-60 families** — they
must be ADDED on freeze (K^0.84 sublinear exponent, Q^3.33/Q^1.19 inflation, N^1.85/N^1.70
SW/L² conditioning, N^1.97 water, n_orb^2.2 molecular 1-norm, the tab:resource d_inv/kappa
row family, He −2.847/−2.897/−2.90372).

## Dimensions exercised (ALL, one invocation — FULL run, first cert of a fresh target)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1 on Paper 60. **RUN
  `pytest tests/test_paper60_sturmian.py --slow` (17 tests; `--slow` MANDATORY — the 8
  default-skipped tests carry the sublinear-1-norm, gerade-conditioning, water-probe, and
  H2-CI headlines).** Map each headline to its test and audit whether the test *proves* the
  claim (not tautological / weaker than prose): potential-weighted orthonormality (diag);
  L2 ill-conditioning growth; F0=5/8; single-config He variational −2.847/−2.84766;
  multiconfig lowering + **sublinear 1-norm** (`multiconfig_lowers_and_sublinear_onenorm`,
  `lgt0_angular_correlation_and_sublinear`); SW intra-center = identity; SW better-
  conditioned than L2; **gerade flat conditioning** kappa≈2; SW beats Gaussian metric; H2+
  isoenergetic binds; many-electron metric-does-not-compound; two-center ERI vs exact
  closed form; interacting-H2 binds; **collective scale = root-sum-of-squares** (the −4 vs
  −8 gate that caught the T0 bug); H2 CI dissociates+binds; **water gerade-lever-fails-for-
  inequivalent-center**.
- **Paper claims / prose (C3, C5, C6, C8)** — `claims-reviewer` ×1 on Paper 60,
  **enumeration-forced**: every row of `tab:resource` AND the SW-conditioning table AND the
  He convergence chain, every `\textbf{[TIER]}` label, every headline number. Special
  charge: adjudicate the two abstract↔body number tensions (W1, W2) — a self-inconsistent
  headline number is a defect; quote the abstract clause and the body/labelled-equation
  clause side by side.
- **External citations (C4)** — `citation-reviewer` ×1 on Paper 60. Cites: `avery1989`,
  `avery2006`, `averyphd`, `averymsc` (Avery generalized-Sturmian / isoenergetic /
  Goscinskian); `babbush2018`, `su2021` (plane-wave first-quantization "compute-don't-load"
  precedent); `liang2022`, `baek2023`, `rajchel2025` (quantum generalized-eigenvalue);
  `cks2017`, `gslw2019`, `lowchuang2019` (QSVT / qubitization primitives); `calderini2012`;
  `herbst2019` (Herbst–Avery–Dreuw CS-HF). **Verify `rajchel2025` (W7 — known title/author
  drift).**
- **Synthesis faithfulness (C9)** — **NOT EXERCISED (out of scope, no footprint).** Recorded
  as a scope exclusion, NOT a skipped gating dimension: C9 does not gate this verdict
  because Paper 60 makes no synthesis claim to audit.
- **Deterministic (C10–C18)** — the step-1 scripts, `--gate group2`. **Add Paper-60 C17
  families first** (see C17 note).
- **Completeness-critic** ×1 (FULL run).

## Branch-defining criterion (Paper-60-specific, PI to confirm): QC-resource honesty

Paper 60's single highest risk is overselling a quantum-resource advantage it does not have.
The reviewers must verify ALL of:
1. **The atomic sublinearity is a CONFIGURATION-count (K) statement, not a qubit-count (Q)
   one.** ‖M‖1∼K^0.84 is sublinear in CI configuration count for a single atom; the paper
   itself flags it as not-yet-mapped to Q. No prose may imply a qubit-count sublinearity or
   a many-electron win. **exact ≠ accurate (inherited Paper 58 W1):** metric-free / pi-free /
   pure-number structure buys *encoding cost*, not accuracy.
2. **Molecules are polynomial, and the paper must say so.** The N-electron interacting
   molecular block-encoding 1-norm is ∼n_orb^2.2 (standard second-quantization ballpark, no
   advantage over DF/THC); the paper concedes "not a sublinear matrix." Any prose implying a
   molecular many-electron 1-norm advantage = MATERIAL.
3. **Benchmarking rule (strongest baseline).** DF / THC / sparse-qubitization / plane-wave
   first-quantization are the mature FT baselines; Paper 60 does NOT beat them at molecular
   scale and must not imply it does. Honest ceiling = "novel metric-free ATOMIC secular
   equation + conceptual kinship with compute-don't-load," not a competitive molecular lambda.
4. **Metric-vs-metric ratios carry the ratio-dependent caveat.** The 10^2–10^4× vs Gaussian
   is ratio-dependent (Gaussian ratio-1.6 → cond∼N^6; ratio-3 → N^1.4). Caveat must be intact.
5. **No hardware.** Resource/structural analysis; validation is textbook/variational, not
   production chemistry. Any "demonstrated on a device" reading = MATERIAL.
6. **Novelty scoped to the exact triple.** The [OBSERVATION] novelty is "no existing quantum
   algorithm combining (i) Sturmian basis, (ii) quantum eigenvalue routine, (iii) isoenergetic
   inversion" — not "first quantum-chemistry algorithm" or any broader claim.

## Paper-60-specific watch-notes (the risk surface — ranked)

- **W1 — abstract K-exponent drift [HIGHEST, headline-number].** Abstract (line 46)
  ‖M‖1∼K^0.78; body `eq:sublinear` (labelled, line 247) ‖M‖1∼K^0.84; CLAUDE §2 + the
  lit-comparison memo both use 0.84. **Canonical = the labelled body equation, 0.84.** The
  abstract figure must match. Mismatch = MATERIAL (C8/C17).
- **W2 — gerade-vs-Gaussian ratio drift [HIGH, headline-number].** Abstract (line 62)
  attributes "10^2–10^3× smaller than Gaussian" to the *gerade*; body (line 409) gives
  full-metric = 10^2–3×10^2 and **gerade = 10^3–10^4×** — consistent with `tab:resource`
  (3.3×10^5/18 ≈ 1.8×10^4 gerade; 3.3×10^5/1050 ≈ 3×10^2 full). The abstract appears to have
  mis-attributed the full-metric range to the gerade. **Canonical = the table.** Mismatch =
  MATERIAL.
- **W3 — sublinear-axis conflation [framing-zombie].** The atomic K^0.84 (config-count
  LCU-λ) and the plane-wave "sublinear in basis size N" (Babbush 2019, Toffoli-in-N) are
  DIFFERENT sublinearities on different axes/mechanisms; the lit-memo flags reader-conflation
  risk. The paper must not blur them.
- **W4 — "metric-free"/"no metric" is ATOMS-ONLY.** Molecules re-introduce the SW metric
  (generalized eigenproblem). No prose may present the metric-free result as general.
- **W5 — the gerade lever is EQUIVALENT-CENTER-only.** [MEASURED] water probe: cond(A1)∼N^1.97,
  the lever fails for a symmetry-unique heavy center (O↔H coupling is the sole driver). Must
  not be presented as a general polyatomic property.
- **W6 — He is DELIBERATELY low-accuracy.** −2.897 (spdf, K=164) vs exact −2.90372 (~7 mHa)
  is worse than STO-3G-class; the Goscinskian basis is deliberately poor for the He GS. The
  validation proves the *machinery is correct*, NOT that it is accurate. "Accurate helium"
  reading = MATERIAL.
- **W7 — `rajchel2025` citation drift [C4].** Lit-memo: bibitem title ("…via a singularity
  test") ≠ real arXiv title ("Quantum algorithm for solving generalized eigenvalue problems
  with application to the Schrödinger equation"); initials "K. Plis"/"A. Zak" ≠ real (Szymon
  Plis; Emil Zak). Fix-on-sight; verify vs arXiv:2506.13534.
- **W8 — H2+/H2 are TOY validations.** H2+ 1.3% s-only underbound; H2 CI −1.09 vs exact
  −1.174. They validate structure (dissociation, binding, metric behavior), not accuracy.
- **W9 — `tab:resource` is a [RESOURCE MODEL] on measured kappa.** The d_inv values are modelled
  (O(1) constant cancels in cross-metric ratios, not absolute counts — per the caption). Not
  measured circuit counts.
- **W10 — transcendental / pi-free tagging.** The load-bearing "pure numbers" claim (T′
  independent of p_kappa/E/Z; F0=5/8; sqrt2 in the He scale) is the paper's pi-free content; verify
  it is exact-rational/algebraic as stated, no anonymous transcendental.

## C8 headlines (enumerated, with tiers — the frozen goalposts; canonical = BODY values)

1. **Naive L2 inflation [MEASURED].** Shared-scale Coulomb-Sturmian is L2-non-orthogonal;
   JW/Löwdin LCU 1-norm inflates λ∼Q^3.33 vs Q^1.19 hydrogenic; driven by overlap
   ill-conditioning.
2. **Isoenergetic metric-free ATOMIC secular equation [ESTABLISHED, from Avery].**
   [diag(Z R_nu)+T′−p_kappa·1]B=0; eigenvalues p_kappa=sqrt(−2E) = the energies directly (no outer
   loop); T′ = a matrix of pure numbers, independent of p_kappa, E, and Z.
3. **Helium validation [MEASURED].** Single Goscinskian config −2.847 Ha = textbook
   variational (−2.84766; the pure number 5/8·sqrt2^−1); bare (no V′) matrix −4.0 (two He+ 1s,
   exact non-interacting); multiconfig −2.847(1s^2)→−2.873(s)→−2.894(+p)→−2.897(spdf,K=164) →
   exact −2.90372; ~7 mHa residual = basis incompleteness (deliberately poor basis; Avery &
   Avery reach −2.90250 with 102 configs), every point above exact (no overshoot).
4. **Atomic sublinear 1-norm [MEASURED].** ‖M‖1∼**K^0.84** (`eq:sublinear`, full s+p+d+f) in
   configuration count K — sublinear, opposite of the naive inflation. (CANONICAL 0.84;
   abstract 0.78 = the W1 drift.)
5. **Novelty [OBSERVATION].** No prior quantum algorithm combines (i) Sturmian/hyperspherical
   basis, (ii) a quantum eigenvalue routine, (iii) the isoenergetic inversion; the two
   documented cost risks (metric conditioning, outer energy search) are both ABSENT in the
   atomic case.
6. **SW molecular metric [MEASURED].** Better-conditioned than L2 (intra-center block =
   EXACT identity; better at large R) but still grows: cond(S)∼N^1.85 (≈ L2 N^1.70;
   power-law not exponential; lambda_max<2, all growth in lambda_min→0), inter-center-coupling
   driven (off-block ~0.4–0.6). Cross-checked 3 ways (2D position grid / closed-form momentum
   quadrature / 3D momentum integral) to ~10^−10.
7. **Gerade lever + resource model [RESOURCE MODEL / MEASURED].** H2+ sigma_g is gerade; gerade
   sector flat kappa≈2 across N=4–20; ground-state d_inv≈16–18 basis-independent → ~60× vs full
   SW metric (d_inv 1050→18). `tab:resource` (N=16): atomic 23q/kappa1/d_inv0; SW gerade
   26q/kappa2.3/d_inv18; SW full 31q/kappa91.8/d_inv1050; Gaussian LCAO(r2) 31q/kappa1.99×10^4/d_inv3.3×10^5.
   SW vs Gaussian: 10^2–3×10^2 (full) / 10^3–10^4× (gerade) [body/table]; ratio-dependent caveat.
   H2+ binds 1.3% s-only.
8. **Metric cost is ONE-electron [MEASURED].** k-electron config overlap = k-th compound
   matrix of the 1e overlap; SW-orthonormal molecular Sturmians → identity metric at every k;
   metric discharged once at the 1e MO construction — does not compound with electron number.
9. **Water gerade-lever-FAILS [MEASURED].** C2v water A1 ground-state block cond∼N^1.97
   (R^2=1.000; 19.9→698 over N=6→36), NOT flat; zeroing the O↔H block → cond≈2. The lever is
   equivalent-center-only, not symmetry-as-such.
10. **N-electron 1-norm NOT sublinear [MEASURED; the OPEN question answered negative].**
    Two-center CI (H2 dissociates correctly, E→−1.0 at R=6; underbinds −1.09 vs exact −1.174
    at R_eq) has standard block-encoding λ∼n_orb^2.2 — polynomial, no sublinear behaviour;
    the atomic sublinearity rides on the single-center diagonal T0=Z R_nu. Real molecular
    savings = the metric levers (gerade / large-R / one-electron), not a sublinear matrix.

## Seeding plan (worktree only; never touches the real corpus)

K ≈ 5–6 planted defects, ≥1 catchable by each EXERCISED dimension (code / prose / citation;
C9 not exercised), spanning the watch-notes — e.g. an **accurate-helium** overclaim (W6), a
**"molecular many-electron 1-norm is sublinear"** reversal (branch-crit #2 / C8#10), a
**"metric-free in general"** over-generalization (W4), a **"gerade lever holds for water/
polyatomics"** reversal (W5/C8#9), a **citation splice** on `rajchel2025` or `babbush2018`
(C4/W7), a **K^0.84 → "derived/proven-optimal"** tier overclaim (C3/S4). M ≈ 5 known-good
controls (verified C8 headlines: F0=5/8, He −2.847, cond(S)∼N^1.85, a `tab:resource` row,
n_orb^2.2) that must NOT be flagged. Tiered agents (code + citation = Sonnet) get 2 seeds
each. Answer key → `debug/qa/paper_60_seed_key.json`.

## Change log
- 2026-08-18 — **FINAL FULL certifying run = PASS → CERTIFIED ✅** (PI: "run one final formal certifying
  pass"). Fresh git-worktree over the post-sprint text (full `geovac` package so the tracked suite runs),
  6 fresh seeds across all 3 dimensions. **Calibration 6/6, specificity clean, zero cert-blocking MATERIAL.**
  Code re-ran the tracked suite (28 passed) and independently recomputed every headline number (all
  reproduce). One fix-on-sight applied mid-run (L²-divergence tolerance `>100` → band `2000<c<6000` around
  the real 3672.6). Completeness-critic re-caught both claims seeds, no new material. Worktree removed, no
  leak. Honest ceiling: the Avery book/thesis content attributions stay UNVERIFIABLE (inaccessible primary
  source); the exact `−p_κ𝟙` secular form is Avery-formalism primary-source-gated — standing limits, not
  GeoVac coverage gaps.
- 2026-08-18 — **ENGINE SPRINT (v4.96.0) — the full-cert weak-backing CLOSED.** PI-directed after
  a "talk about the engine" discussion. Promoted the Goscinskian machinery from untracked `debug/`
  drivers to four tracked, type-hinted, closed-form-validated `geovac/sturmian_*.py` modules
  (`sturmian_integrals`, `sturmian_l2_encoding`, `sturmian_secular`, `sturmian_molecular_lambda`) +
  faithful ≥4-pt tests, via a main-session engine foundation + **three parallel target agents**
  (Q/K/N). Exact exponents now regression-protected: Q^3.33/Q^1.19 (measured **3.333/1.191**, JW LCU
  1-norm), K^0.84 (**0.842**), n_orb^2.2 (**2.22**); + L²-divergence (cond 4.1→**3672.6**), eq:secular
  T'-Z-indep (**2.2e-16**), large-R lever (cond→2.2 @R=10), and the manyelectron SW-side de-decorated.
  Fix-on-sight applied: babbush2018 (not first-quantized), eq:blowup caveat, single-ζ wording. Full
  regression **28 passed** + topological **18 passed**; paper recompiles clean (7 pp); deterministic
  C10–C18 green (C17 guards K^0.84/n_orb^2.2). The residual to CERTIFIED is a final formal certifying pass.
- 2026-08-18 — **FIRST-CERT FULL run = FAIL (trustworthy).** Panel FULLY CALIBRATED: sensitivity **6/6**
  (code FC1 false-positive-tol + FC2 F0-tautology; citation FT1 liang2022-wrong-ID + FT2 cks2017-vol;
  claims FP1 chemically-accurate-overclaim + FP2 tier-flip), specificity clean (all reviewers verified
  the paper's numbers reproduce closely — no wrong-value false-flags). **Zero wrong §C8 headline value.**
  Genuine non-seed findings were all fix-on-sight NITs (babbush2018 mischaracterization; textbook-single-ζ
  uncited; Acknowledgments blanket) + *weak-backing/coverage* (exact exponents + large-R/L²-divergence/
  eq:secular driver-only) → raised to PI (§9) → resolved by the engine sprint above. Seed key
  `debug/qa/paper_60_full2_seed_key.json`; worktree removed, no seed leaked.
- 2026-08-18 — **DELTA-verification run = CLEAN-DELTA ✅.** Diff-scoped over the remediation
  edits, fresh seeded delta worktree, all three affected dimensions. **Calibration: 5/5
  seeds caught** (claims DC1 tier-overclaim; citation DT1 wrong-initial + DT2 wrong-year;
  code DK1 gutted-ordering + DK2 useless-band), **specificity clean**. **Real edits verified
  correct across every dimension:** claims — Fixes 1–4 FIX-OK (abstract `K^0.84`↔`eq:sublinear`
  consistent, gerade `10³–10⁴`↔tab:resource, `N^1.85` consistent, Acknowledgments honestly
  scoped, no self-contradiction/hard-prohibition drift); citation — `baek2023` = `U.~Baek`
  (2023) confirmed correct vs arXiv:2205.09039 + PRX Quantum 4, 030307; code — the tight
  assertions (`p_sh>p_hy+0.5`, band `0.6<p<0.95`, cond(S) 3.0/5.8/13.9/32.2 sequence, He
  variational) confirmed genuine strengthenings. **Zero genuine (non-seed) defect.** Delta key
  `debug/qa/paper_60_delta_seed_key.json`; worktree removed, no seed leaked. **A clean delta is
  the precondition for the FULL certifying run** (the next and final step for a PASS).
- 2026-08-18 — **REMEDIATION applied (PI-directed), toward delta re-run.**
  **Mechanical fixes:** W1 (`K^0.78`→`K^0.84` headline + s-only clarified), W2 (gerade
  `10²–10³`→`10³–10⁴` + ratio-dependent note), `baek2023` `J.`→`U.` (verified vs
  arXiv:2205.09039), `N^1.8`→`N^1.85`, Acknowledgments umbrella softened, `averyphd` cap;
  + a sec:atomic note that the quoted exponents are numerical fits. Paper recompiles clean
  (7 pp, 0 undefined). **Test strengthenings (17/17 `--slow` green):** L²-inflation → pins
  the mechanism (cond(S) 3.0/5.8/13.9/32.2 sequence + shared>hydrogenic Löwdin-λ ordering,
  Q^1.86 vs Q^0.76), self-contained; atomic-sublinear → ≥5-pt full s+p+d+f band fit
  `0.6<p<0.95` (measured 0.77). **C17:** added `paper60-atomic-sublinear-exponent` (K^0.84
  headline-form; ignores bare s-only 0.78) + `paper60-molecular-lambda-exponent` (n_orb^2.2)
  families — both discriminating (fire on wrong value, PASS on corrected paper). Deterministic
  C10–C18 `--gate group2` all green. **PI decision (exponent backing):** the EXACT
  K^0.84/n_orb^2.2/Q^3.33/Q^1.19 are accepted as driver-measured + logged in
  `docs/claim_test_matrix.md` with an engine-promotion follow-up (promote the untracked
  Goscinskian λ engine to tracked code); the tracked tests pin the regime/band/mechanism/
  ordering. Next: delta-verification re-run (diff-scoped, seeded) → FULL certifying run.
- 2026-08-18 — **DRAFTED** by PM for PI freeze. Third single-paper `/qa` target; capstone of
  the Avery arc (58→59→60). Inherits criteria.md C1–C18 + the group2 benchmarking/guardrail
  deltas; adds a Paper-60-specific QC-resource-honesty branch criterion. C9 out of scope (no
  synthesis footprint). Backing suite 17/17 green (`--slow`). **Two genuine abstract↔body
  headline-number tensions surfaced during drafting — RAISED TO PI (fix-at-source before
  freeze vs let-the-run-catch):** W1 (K^0.78 abstract vs K^0.84 labelled body/`eq:sublinear`)
  and W2 (gerade-vs-Gaussian 10^2–10^3× abstract vs 10^3–10^4× body/`tab:resource`). C17
  registry has NO Paper-60 families yet — must be added on freeze.
