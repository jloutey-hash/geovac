<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED

> **SUPERSESSION NOTE (2026-09-11, v5.10.18).** This record verified `cond(S)~N^1.85` as a C8 headline and treated `eq:sigma_law` as a result derived in this corpus. The *values* stand; the *attribution* does not. The law is the Kac-Murdock-Szego extreme-eigenvalue asymptotic (c_1 = pi^2, 1953), identified by a literature scan on 2026-09-11; Paper 60 now cites it and claims only the identification of the Shibuya-Wulfman metric as such a finite section. Do not read this record as ratifying the originality of `eq:sigma_law`. Backing: `tests/test_paper60_kms_attribution.py`.

> This record certifies the state as of **2026-08-18**. Since then **1 `.tex` changed**: paper_60_sturmian_secular_quantum.tex.
>
> **The CERTIFIED verdict below is therefore historical, not current.** Do not cite it as present-tense status.
>
> Re-measure rather than trusting this banner — it is itself a snapshot and will go stale the same way:
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Paper 60 (Generalized-Sturmian Secular Equation) — `/qa` profile

> **RUN 2026-09-11 — DELTA verification, unseeded: DEFECTS, remediated. NOT certified.**
> Five dimensions. Deterministic layer 12/12. Citations CLEAN. **2 LARGE + 25
> MATERIAL-SMALL + 5 upgrade candidates, and zero mathematical defects** — every
> finding was prose, attribution, staleness or coverage. Five adversarial passes,
> one of which re-derived every keystone independently and fire-tested the guards
> it was handed, found nothing wrong with the arithmetic.
>
> **LARGE 1 (corpus surface).** `CLAUDE.md:119` asserted four retired values
> (`K^0.84`, "local slope 0.906", "T⁰ is the sublinear block at K^0.70", "T′ is
> SUPERlinear at K^1.05") with no supersession marker, in the file loaded by every
> session and every subagent dispatch — two lines below the bullet that supersedes
> it. §13.11 rule 9 was not applied: four newer bullets were appended and this one
> left standing. Replaced; superseded text relocated to the frontier archive.
>
> **LARGE 2 (coverage) — CLOSED 2026-09-11**, as its own pass per §9. Six
> abstract-level `[MEASURED]` families had driver-only backing in prunable
> `debug/`; all six are now recomputed from tracked `geovac/` code in
> `tests/test_paper60_resource_ladder.py` (7 tests, all `@pytest.mark.slow`,
> 14 min). Four were fire-tested against the specific wrong answer each excludes
> — including the K=105/K=136 substitution that was this run's one wrong number,
> which FIRES. The seventh family, the floor **bracket**, is **PARTIAL and
> declared**: the claim form (fit-from-below vs Shanks-from-above) is backed on
> the spdf ladder K=74–244; the endpoint values [6.47, 6.62] / [1.647, 1.676]
> need the full ladder to K=452 (~20 min) and remain driver-backed, recorded as
> PARTIAL in `docs/claim_test_matrix.md`.
>
> Two things the writing found that the review had not. **(i)** The
> bracket-direction claim is SECTOR-SPECIFIC: written first on the cheap s-only
> ladder, it failed, because there the windowed fit *falls* (4.3098 → 4.3059 →
> 4.3035) where on the spdf ladder it *rises*. A cheap proxy in the wrong sector
> would have reported the claim backed while measuring something that behaves
> oppositely. **(ii)** A seventh locus of the "ill-conditioned" cluster, in
> `docs/claim_test_matrix.md` and in a backing test **named**
> `test_paper60_l2_overlap_illconditioned_grows` — the test-asserts-the-zombie
> sub-flavour catalogued in v4.43.5. Its assertions were sound and are unchanged;
> the name and framing were not. C16's pattern was widened twice to reach the
> bare adjective (post- and pre-nominal), and re-proved silent on the paper's own
> denial.
>
> **The one wrong number in the paper:** the s-only span deficit read `4.43` at
> K=136; the driver's own output gives `4.40376` there and `4.43413` at K=105 — a
> mismatched pair, in the sentence carrying the floor's re-attribution. It was an
> *unregistered* literal, so C21 was blind to it, and it was absent from the
> declared-debt table that asserted its list was all correct. Now registered as a
> matched pair (`p60_span_deficit_sonly_locked` / `_free`).
>
> **Instrument findings (three), all fixed before any content edit:** no C16 entry
> could reach `CLAUDE.md`, `claims_register.md`, `code_architecture.md` or
> `geovac/sturmian_variational.py`; `p60-sublinear-as-regime` stated its own
> CORRECTION in values retired the next day; and its `exempt_if_nearby` contained
> `nuclear diagonal|T^0` — **the vocabulary of the retired mechanism**, so a locus
> asserting "the nuclear diagonal T⁰ is the sublinear block" exempted itself. Two
> new entries (`p60-tprime-superlinear`, `p60-l2-metric-diverges`) use the
> standardized per-entry marker and `exempt_if_nearby = (?!)`, because CLAUDE.md:119
> sits two lines from a bullet containing both "WITHDRAWN" and "2026-09-08" — any
> plausible nearby vocabulary would have exempted the LARGE. Both proven to
> discriminate two ways: 6/6 fire on retired wording, 8/8 silent on corrected
> wording, including two loci where an early draft of my own pattern fired on the
> *denial* ("Neither T′ grouping is superlinear") and would have made the right
> answer unwritable.
>
> **Two upgrades** (C8.15/C8.16 below) — the backing proved more than the prose on
> both results this week rests on.
>
> Full sweep: 4 live zombies fixed, 16 chronicle loci marked, 6-locus
> "ill-conditioned" cluster swept claim-wide, `reverses to superlinear` corrected at
> owner + citer + register, the floor/bracket contradiction resolved, the
> Acknowledgments' four surrendered theorems reclaimed, and the `cor:dual_p0`
> overstatement fixed in the auto-loaded memory and in `paper_fci_molecules`.



> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-60 scope + deltas + watch-notes.

> **STATUS: 2026-09-07 — FAIL, TWICE; NOT re-certified.** Remediated on the
> second pass. The paper's central claim CHANGED, it was not corrected.
>
> The morning run found `eq:sublinear` overclaimed (a window fit stated as a
> regime, with the mechanism attributed to the wrong block) and remediated it.
> **That remediation was then itself refuted the same day** — its replacement
> numbers were measured on the same truncated radial domain as the claim they
> replaced. `geovac/sturmian_secular.py` fixes `R_MAX = 60.0` while the
> Goscinskian mixed pairs carry orbitals reaching 100–200 bohr, and
> `hyd_radial` renormalises the truncated stub to unit norm — manufacturing a
> compact pseudo-orbital whose coupling is O(1) where the true one decays as
> `n^-2`. The error therefore GROWS with K, the fit variable (×1.011 → ×1.205
> across the fitted range), which is what produced the apparent rising slope.
>
> **Settled by two independent routes agreeing to 4 decimal places** — exact
> grid-free Slater algebra (dps 60–70) and converged quadrature under the derived
> rule `R_MAX ≥ 3n_max²`. Both also reproduce the paper's published numbers when
> run on the old domain, so this is not a pipeline difference.
>
> | leg | window | published | **converged** |
> |---|---|---|---|
> | `‖M‖₁` | 74–164 | 0.84 | **0.8193** |
> | `‖M‖₁` | 74–340 | 0.854 | **0.8058** |
> | `‖T′‖₁^off` | 74–340 | 1.05 | **0.9368** |
> | `‖T′‖₁` full | 74–340 | (prose: superlinear) | **0.8750** |
> | local slopes | 100–514 | rising to 0.906 | **falling to 0.766** |
>
> **What the paper now says.** The sublinearity belongs to the *basis-growth
> rule* (fixed `l_max`), not to the isoenergetic construction — grown as full
> hydrogenic shells the same construction gives `K^1.07` and a *rising* exponent.
> And the cheap rule carries an accuracy floor of **6.44 mHa**.
>
> **Corrected 2026-09-13:** the floor is quoted as a BRACKET, **[6.47, 6.62] mHa**,
> not as a single fitted value — the paper says so in its own honest-scope sentence
> ("it is the bracket we quote"), because the fit family approaches the floor from
> below and a single number would be a lower estimate rather than a central one.
> The 6.44 above sits BELOW the bracket's own lower endpoint. This record ratified it.
>
> **Superseded 2026-09-08 — the floor's mechanism and universality.** This
> record previously ratified "4.0× chemical accuracy, at any basis size. Cost
> growth and attainable accuracy are one fact." The floor VALUE stands; its
> attribution and its universality do not. (i) The mechanism is the **scale
> lock**, not `l_max`: metric-free holds iff `E = -λ²/2`, hence
> `λ = p_κ`, and that is not the variational optimum — freeing `λ`
> over the *identical* span reaches 1.28 mHa at K=130 against 7.46 locked, and
> 0.15 mHa of the independently known s-limit at K=136. The He `l≥4`
> partial-wave tail is 0.37–0.53 mHa, an order below the floor, so angular
> truncation cannot be it. (ii) The floor is **ground-state specific**: at
> K=202, with `‖M‖₁` identical because it does not depend on which root
> is extracted, the ground state sits 4.49× above chemical accuracy and
> 2¹S sits 1.12×; the posing cost falls 3–4× per rung up the ¹S
> ladder. (iii) The price of freeing the scale is the whole encoding advantage,
> `‖·‖₁` from `K^0.72` to `K^2.75` (s-only ladder, K=21..136 — not the
> headline window). Backing:
> `tests/test_paper60_scale_lock.py`; drivers `debug/p60_{variational_probe,
> scale_scan,freescale_resource,posing_cost_by_state,excited_ladder}.py`. `‖T⁰‖₁` is now `[SYMBOLIC]` and is not a
> power law at all — `Z√(2K)·ln(K/2)`, asymptotic exponent 1/2.
>
> **Second LARGE, a different claim.** `cond(S)` “4 → 3673” is entirely a domain
> artifact (converged: 4.07 → 16.0 → 23.5, growing ~0.12·K — an ordinary Gram
> matrix). It switches on exactly where `n_max²` first exceeds `R_MAX`, and every
> domain agrees to 4 digits below that point. **`test_sturmian_secular.py`
> asserted `2000 < cond < 6000`** — it pinned the artifact and would have failed
> on repair. Withdrawn from the paper; test replaced and fire-tested. The
> SEPARATE `eq:blowup` cost result (L² Löwdin LCU 1-norm ~Q^3.33 vs Q^1.19,
> through a faithful Jordan–Wigner LCU) is a different quantity on a different
> route and is **untouched** — the paper still has its obstruction.
>
> **Two guards retired for defending withdrawn claims:** the `cond(S)` test
> above, and `test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal`,
> written that same morning, which asserted `p_off > 1.0`.
>
> ---
>
> ### Branch-defining criterion (REPLACED 2026-09-07)
>
> It is no longer “does the exponent name its window.” It is:
>
> **Every quantity in this paper must name its evaluation domain, its nuclear
> charge, and its basis-growth family — and a run must verify the domain
> satisfies `R_MAX ≥ 3n_max²` (App. A) before accepting any exponent.**
>
> The three axes are not decoration. `p_total` moves 0.947 (Z=½) → 0.717 (Z=30),
> so **“0.84” was a helium number**; it moves ±0.05 with `l_max`; and it reverses
> its drift under a different growth rule. Any quoted exponent missing one of the
> three is a defect.
>
> **Restricted-evaluation is the top watch-note for this paper**, not a general
> mandate: two successive remediations were both defeated by the same undeclared
> restriction, and `E(1s²) = −729/256` is bit-identical at every domain, so the
> module's own calibration cannot see it.
>
> *(superseded — historical)* **RE-RUN 2026-09-07 — FAIL (one LARGE), remediated.**
> `/qa` (PI-invoked).
>
> **LARGE — `eq:sublinear` was stated as a regime, and its mechanism was
> backwards.** Re-measured with the paper's own `gen_configs`/`solve` past its
> largest fitted point (K = 164) to K = 340:
>
> | | window K ≤ 164 | extended K ≤ 340 |
> |---|---|---|
> | total ‖M‖₁ | **K^0.840** (reproduces the published 0.84) | K^0.854 |
> | local slope | 0.838, 0.840, 0.841 | 0.850, 0.868, 0.882, **0.906** |
> | nuclear diagonal T⁰ = Z·ΣR_ν | — | **K^0.704**, local slope FALLING |
> | off-diagonal T′ (the "pure numbers") | — | **K^1.049**, SUPERlinear |
>
> So (i) 0.84 is a **window fit**, not an asymptotic regime — the exponent
> trends toward 1; and (ii) the sublinearity is carried **entirely by the
> nuclear diagonal**, while the pure-number block the paper credited is the
> superlinear one. §7 (molecular) already said this correctly ("rides on the
> clean diagonal T⁰"), so §4 contradicted §7 inside one document.
> *What survives:* ‖M‖₁ does grow more slowly than the matrix dimension over
> every computable basis, and the contrast with the L² superlinear inflation is
> real. The encoding claim stands; its asymptotic reading and mechanism did not.
>
> Note for future runs: the corpus carried **three different windows** for this
> one exponent — K = 9..100 (0.842, `test_sturmian_secular.py`), K ≤ 24 (0.77,
> the self-contained sweep), K = 74..164 (0.84, the paper) — each reported as
> "the" value. That is why the window is now part of the registered object.
>
> **Remediated:** abstract / intro / `eq:sublinear` / mechanism / §5 / conclusion
> restated; new `eq:sublinear_split` records the T⁰/T′ split; group2 synthesis
> block rewritten (it carried the bare "sublinear 1-norm" — the LARGE's second
> locus) together with two SMALLs there ("the one genuine lever" → the paper
> names **three**; the isoenergetic posing is **Avery's** and was uncredited);
> 4 numeric-registry entries + 5 `\gvq` annotations (C21 was examining ZERO
> here, now 8); C16 entry `p60-sublinear-as-regime` with `cited_by`, proven to
> discriminate both ways; claim-matrix row re-tiered, split row added, three
> uses of "asymptote" corrected; backing test
> `test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal` written as a
> separate activity and fire-tested in both directions.
>
> *(superseded — historical)* **CERTIFIED ✅ 2026-08-18 — final FULL certifying run = PASS.** Over the post-engine-sprint
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
- **C9 (synthesis): GATING (corrected 2026-09-07).** This record said until then:
  *"OUT OF SCOPE — Paper 60 has no group2-synthesis footprint yet (verified by grep)."*
  **True when frozen (2026-08-18), false now.** The group2 synthesis gained
  `\subsection{The isoenergetic secular equation as a quantum algorithm (Paper~60)}`
  (L692–716) plus L43/L76/L259/L778/L942–945/L981–982 and bibitem L1168 on 2026-09-06
  (v5.10.6) — and that un-certified block already carries a defect the FULL run caught
  from the other side: L703–706 attaches "that a symmetry-unique heavy atom reinstates"
  to the *lever* where Paper 60 says such a center reinstates the **metric cost**, i.e.
  a documented negative inverted into a positive by a relative clause.
  C9 is therefore a **gating dimension** for Paper 60 and was **not exercised** in the
  FULL run of 2026-09-07. See the identical correction in `paper_59.done.md`.
- **Out of scope:** the other group2 papers (unchanged since 2026-06-28 cert). Papers
  58/59 (`loutey_paper58`/`loutey_paper59`) are CERTIFIED; in scope only as cross-refs
  (C7). Trunk papers (0/1/7/14/18) canonical; in scope only where Paper 60 restates them.

**Deterministic `--gate`:** `paper_60`.
  Note the scope is the SINGLE-PAPER one, not `group2`: `qa_scopes.py` deliberately excludes 58/59/60 from `group2` (they are their own cert targets), so the `--gate group2` this record carried until 2026-09-07 examined NOTHING for the paper it certifies.
**C17 note (CORRECTED 2026-09-12 — this note froze two RETIRED values as
goalposts, the same class the file records as having stopped the 2026-09-11
run at protocol step 1):** two Paper-60 families EXIST
(`paper60-atomic-sublinear-exponent`, `paper60-molecular-lambda-exponent`,
added 2026-08-18). Do NOT register `K^0.84` (retired 2026-09-07 →
`p60_onenorm_exponent` = 0.82) or He `−2.897` (retired 2026-09-12 →
`p60_he_chain_spdf_k164` = −2.8964). Original note, kept for the record: the
registry has NO Paper-60 families — they
must be ADDED on freeze (K^0.84 sublinear exponent, Q^3.33/Q^1.19 inflation, N^1.85/N^1.70
SW/L² conditioning, N^1.97 water, n_orb^2.2 molecular 1-norm, the tab:resource d_inv/kappa
row family, He −2.847/−2.897/−2.90372).

## Dimensions exercised (ALL, one invocation — FULL run, first cert of a fresh target)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1 on Paper 60. **RUN
  `pytest tests/test_paper60_sturmian.py --slow` (17 tests; `--slow` MANDATORY — the 8
  default-skipped tests carry the sublinear-1-norm, gerade-conditioning, water-probe, and
  H2-CI headlines).** Map each headline to its test and audit whether the test *proves* the
  claim (not tautological / weaker than prose): potential-weighted orthonormality (diag);
  L2 overlap NON-ORTHOGONALITY growth (cond is modest — do NOT audit this as
  ill-conditioning, that reading is withdrawn); F0=5/8; single-config He
  variational −2.847/−2.84766;
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
- **Synthesis faithfulness (C9)** — **GATING (corrected 2026-09-07; this line said
  “NOT EXERCISED, out of scope, no footprint” until then).** The premise was false by
  the time it was relied on: the group2 synthesis carries a Paper-60 block, and on
  2026-09-07 that block held the LARGE's second locus — the bare “sublinear 1-norm”
  with none of the paper's qualification — plus two SMALLs (“the one genuine lever”
  where the paper names three; Avery uncredited). A dimension recorded as having no
  footprint held three defects. See the C9 entry in the criteria block above; the
  scope now carries `synthesis/group2_quantum_chemistry_synthesis.tex`.
- **Deterministic (C10–C18)** — the step-1 scripts, `--gate paper_60`. **Add Paper-60 C17
  families first** (see C17 note).
- **Completeness-critic** ×1 (FULL run).

## Branch-defining criterion (Paper-60-specific, PI to confirm): QC-resource honesty

Paper 60's single highest risk is overselling a quantum-resource advantage it does not have.
The reviewers must verify ALL of:
1. **The atomic sublinearity is a CONFIGURATION-count (K) statement, not a qubit-count (Q)
   one.** The exponent of `eq:sublinear` is sublinear in CI *configuration* count for a
   single atom — value owned by C21 key `p60_onenorm_exponent`, and it is a WINDOW fit, not
   a regime. No prose may imply a qubit-count sublinearity or a many-electron win.
   **exact ≠ accurate (inherited Paper 58 W1):** metric-free / pi-free / pure-number
   structure buys *encoding cost*, not accuracy.
2. **Molecules are polynomial, and the paper must say so.** The N-electron interacting
   molecular block-encoding 1-norm is polynomial (standard second-quantization ballpark, no
   advantage over DF/THC) — value owned by C17 family `paper60-molecular-lambda-exponent`;
   the paper concedes "not a sublinear matrix." Any prose implying a molecular many-electron
   1-norm advantage = MATERIAL.
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

## C8.15 / C8.16 — the two results upgraded from MEASURED to proved (2026-09-11)

**C8.15 — `eq:W_diagonal` is [SYMBOLIC], not measured.** The one-body Coulomb
metric is exactly diagonal, by three cases: angular orthogonality when
`l_μ ≠ l_ν`; hermiticity of `(T − E)` on the two Sturmian equations giving
`(Q_μ − Q_ν)·W_μν = 0` when the `l` agree and the roots differ; and disjoint
n-multisets in the degenerate branch, which first occurs at `n_max = 35` (the
largest basis computed is 17). The `5×10⁻¹¹` entrywise agreement is a check on
the implementation, **not the evidence for the claim** — stating it as the
evidence is now an UNDERCLAIM and MATERIAL. The equation is stated at unit scale;
at general λ the diagonal is `λ R_ν`, and dropping that qualifier is MATERIAL
(it makes the labelled equation false by a factor λ).

**C8.16 — the variational bound is proved by inertia, and proves more.** Because
W is diagonal, `H(λ) + ½λ²S = λ(λ𝟙 − M)` **exactly** (verified to 1.4e-17 —
machine precision, i.e. an algebraic identity, not a fit). S ≻ 0, so by
Sylvester's law of inertia the number of pencil roots below `−½λ²` equals
`#{k : λ_k(M) > λ}`: zero at `λ_max` (the bound), exactly k at `λ_k` (the
**root-by-root correspondence**, on which every excited-state number in §4
depends). Verified against a direct pencil solve at k=0,1,2,3. The paper must
NOT revert to asserting "E_iso is the lowest root" as an unproved intermediate —
it is a consequence, not a premise.

## Un-delegated literals — DECLARED DEBT (2026-09-11)

The criteria sections below still write these numbers as literals, because **no
C21 registry key and no C17 family owns them yet** (this paper has exactly two
C17 families: `paper60-atomic-sublinear-exponent`, `paper60-molecular-lambda-exponent`).
Delegating to a key that does not exist would leave the DoD pointing at nothing —
worse than the literal. They are recorded here so the exposure is **visible and
dated** rather than silent:

| locus | literal | status (updated 2026-09-13) |
|---|---|---|
| C8.1 | naive L2 inflation `Q^3.33` vs `Q^1.19` | **REGISTERED** `p60_l2_inflation_sturmian` (3.33, meas 3.3328) / `p60_l2_inflation_hydrogenic` (1.19, meas 1.1909); annotated eq:blowup |
| C8.6 | SW `cond(S)∼N^1.85` (vs L2 `N^1.70`) | **REGISTERED** `p60_sw_cond_exponent` (1.85, meas 1.8518 on window N=12..16; full-range 1.80, asymptote N^2) / `p60_l2_overlap_exponent` (meas **1.71**, paper's 1.70 corrected); annotated sec:molecular |
| W5 / C8.9 | water A1 `cond∼N^1.96` (raw, `_water_A1`) | **REGISTERED** `p60_water_a1_exponent` (1.96, meas 1.9580); annotated sec:resource. The **N^1.97 sec:molecular probe (19.9→698 over N=6..36) is a DIFFERENT construction, still unregistered** — see note |
| branch criterion 4 | Gaussian ratio-dependence (`N^6` / `N^1.4`) | **NOT REGISTRABLE — finding.** The tracked `_gaussian_metric` gives ratio-1.6 → N^6.30 (κ~1.5e6 at N=16, not the paper's ~1e5) and ratio-3 → N^1.54 (not 1.4). The paper's own caveat calls these "illustrative rather than a fixed factor"; treat as illustrative, or promote the original driver |
| W8 | H2+ / H2 toy-validation figures | unregistered (not measured this pass; toy validations) |
| C8.7 | the `tab:resource` row | `[RESOURCE MODEL]` — O(1) cancels in ratios per the caption; modelled, not a measured canonical value |

**Registration pass 2026-09-13 (FULL-run cert-blocker 1).** Five literals were
measured from tracked code (`debug/fullrun_measure_literals.py`) and registered
with `\gvq` annotations; C21 green. **The pass surfaced a finding the earlier
"correct; unregistered" label hid:** three of the debt literals — the Gaussian
ratio exponents, the water `N^1.97` sec:molecular probe, and the `tab:resource`
row — do NOT reproduce from the tracked test suite (the first two trace to
constructions that live only in `debug/` drivers or are modelled). They cannot
be registered as *measured* values without promoting their drivers to `tests/`
first, which is a group2-review task. So the literal-registration debt is now
**partly discharged (5 keys) and partly reclassified** (2 need driver promotion,
1 is modelled) — no longer a flat "register 7 before cert."


**Every one listed above was verified correct on 2026-09-11** — none is a retired value.  *But the table was not exhaustive, and that is the lesson:* the DELTA run found an eighth un-delegated literal, the s-only locked span deficit, written as `4.43` when the measured K=136 value is `4.40` (the 4.43 is the K=105 row).  It was not in this table, so the table's own reassurance did not cover it.  An enumeration offered as complete is a stronger claim than the literals it lists;  it is now registered (`p60_span_deficit_sonly_locked` / `_free`) and this table asserts only what it enumerates. The risk is
structural, not present: each is a literal that will rot the next time its
measurement moves, exactly as `K^0.84` did in W1. **Status 2026-09-13:** the tracked-reproducible five are registered (above);
the three that are not tracked-reproducible (Gaussian ratio, water N^1.97 probe,
tab:resource) need their drivers promoted to `tests/` or reclassification, deferred
to the group2 review. Sec.15 rule 3 was honoured — every registered value was
MEASURED, and the three that could not be measured from tracked code were NOT
registered rather than guessed.

## Paper-60-specific watch-notes (the risk surface — ranked)

- **W1 — K-exponent agreement [HIGHEST, headline-number].** The abstract's K-exponent,
  the labelled body equation `eq:sublinear`, and C21 key `p60_onenorm_exponent` must **all
  three agree**. Any disagreement = MATERIAL (C8/C17/C21). *No value is written here by
  design* — see "DoD criteria name the OWNING GATE" in `criteria.md`. This note previously
  froze `0.84` as canonical and was still asserting it after that value was retired
  (2026-09-07), which is what stopped the 2026-09-11 run at protocol step 1.
- **W2 — gerade-vs-Gaussian ratio drift [HIGH, headline-number].** Abstract (line 62)
  attributes "10^2–10^3× smaller than Gaussian" to the *gerade*; body (line 409) gives
  full-metric = 10^2–3×10^2 and **gerade = 10^3–10^4×** — consistent with `tab:resource`
  (3.3×10^5/18 ≈ 1.8×10^4 gerade; 3.3×10^5/1050 ≈ 3×10^2 full). The abstract appears to have
  mis-attributed the full-metric range to the gerade. **Canonical = the table.** Mismatch =
  MATERIAL.
- **W3 — sublinear-axis conflation [framing-zombie].** The atomic exponent of
  `eq:sublinear` (config-count LCU-λ, C21 `p60_onenorm_exponent`) and the plane-wave
  "sublinear in basis size N" (Babbush 2019, Toffoli-in-N) are
  DIFFERENT sublinearities on different axes/mechanisms; the lit-memo flags reader-conflation
  risk. The paper must not blur them.
- **W4 — "metric-free"/"no metric" is ATOMS-ONLY.** Molecules re-introduce the SW metric
  (generalized eigenproblem). No prose may present the metric-free result as general.
- **W5 — the gerade lever is EQUIVALENT-CENTER-only.** [MEASURED] water probe: cond(A1)∼N^1.97,
  the lever fails for a symmetry-unique heavy center (O↔H coupling is the sole driver). Must
  not be presented as a general polyatomic property.
- **W6 — He is DELIBERATELY low-accuracy, and the RESIDUAL IS THE SCALE LOCK.** The
  helium ladder is worse than STO-3G-class; the validation proves the *machinery is correct*,
  NOT that it is accurate. "Accurate helium" reading = MATERIAL. **Superseded 2026-09-08:**
  this note used to attribute the residual to *basis incompleteness*. It is not — a
  variational CI over the identical span reaches far closer (C21 key
  `p60_span_deficit_spdf`), so the residual is the scale lock, `eq:scale_lock`. Any prose
  still calling it basis incompleteness = MATERIAL.
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

## C8 headlines (enumerated, with tiers — the frozen goalposts)

> **Values are NOT written here.** Each headline names the claim, its tier, and the
> **gate that owns its number** (C21 registry key / C17 family / equation label). See
> "DoD criteria name the OWNING GATE" in `criteria.md`. The pre-2026-09-11 version of this
> block froze three claims that were later retired or withdrawn — `K^0.84` as canonical,
> "residual = basis incompleteness", and the Avery 102-configuration comparison — and would
> have failed a corrected paper.

1. **Naive L2 inflation [MEASURED].** Shared-scale Coulomb-Sturmian is L2-non-orthogonal;
   JW/Löwdin LCU 1-norm inflates λ∼Q^3.33 vs Q^1.19 hydrogenic. **Mechanism corrected
   2026-09-11:** this line read "driven by overlap ill-conditioning" — the
   characterization §2 of the owner WITHDREW on 2026-09-07 (converged cond(S) is
   ordinary, ~0.12·K). The inflation is driven by the *density* of Löwdin's
   S^{-1/2} and the resulting spread of the coefficient distribution, not by
   numerical instability. Asserting ill-conditioning here = MATERIAL; C16
   `p60-l2-metric-diverges` guards it.
2. **Isoenergetic metric-free ATOMIC secular equation [ESTABLISHED, from Avery].**
   [diag(Z R_nu)+T′−p_kappa·1]B=0; eigenvalues p_kappa=sqrt(−2E) = the energies directly (no outer
   loop); T′ = a matrix of pure numbers, independent of p_kappa, E, and Z.
3. **Helium validation [MEASURED].** Single Goscinskian config reproduces the textbook
   single-exponent variational value from the pure number `5/8·sqrt2^-1`; the bare (no V′)
   matrix returns two non-interacting He+ 1s exactly; the multiconfiguration ladder descends
   monotonically toward the exact non-relativistic energy with **every point above it**.
   **The residual is the SCALE LOCK, not basis incompleteness** (`eq:scale_lock`; span
   deficit owned by C21 `p60_span_deficit_spdf`). **The Avery 102-configuration comparison is
   WITHDRAWN** — `eq:no_selection` proves that figure unreachable in the locked posing at any
   K or selection, and a relayed consultation attributes it to a scale-optimized
   (ordinary variational CI) calculation. Any prose reinstating either = MATERIAL.
4. **Atomic sublinear 1-norm [MEASURED].** `eq:sublinear` — ‖M‖₁ grows more slowly than
   the configuration count K on a **converged radial domain**, opposite of the naive
   inflation. Exponent owned by C21 `p60_onenorm_exponent` (s-only companion
   `p60_onenorm_exponent_sonly`); split owned by `p60_T0_asymptotic_exponent` /
   `p60_Tprime_full_exponent`. **It is a WINDOW fit, not an asymptotic regime** — the local
   slope falls monotonically and no window value is stable; prose asserting a regime, or any
   value measured on a truncated domain, = MATERIAL.
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
    at R_eq) has a standard POLYNOMIAL block-encoding λ — exponent owned by C17 family
    `paper60-molecular-lambda-exponent` — with no sublinear behaviour;
    the atomic sublinearity rides on the single-center diagonal T0=Z R_nu. Real molecular
    savings = the metric levers (gerade / large-R / one-electron), not a sublinear matrix.


11. **Metric-free ⟺ the scale is locked to the eigenvalue [INTERNAL THEOREM].**
    `eq:W_diagonal` (the one-body Coulomb metric is exactly diagonal) and `eq:scale_lock`
    (metric-free ⟺ E = −λ²/2 ⟺ λ = p_κ). The posing IS the variational problem of its own
    span at the one scale where the L² metric cancels. Backing
    `tests/test_paper60_scale_lock.py`. **Consequences that must travel with it:** the
    variational bound is automatic, not fortunate; freeing the scale reaches chemical
    accuracy but hands back BOTH advertised cost risks and the encoding advantage
    (`p60_freescale_set_sonly`, a MATCHED SET on one ladder — quoting one member against a
    value from another ladder = MATERIAL); and the floor is a **ground-state pathology**
    (`p60_posing_cost_ground` vs `p60_posing_cost_exc`, `p60_gnd_ratio_k452` vs
    `p60_exc_ratio_k452`), not a property of the method.
12. **No selection rescues the locked posing [INTERNAL THEOREM].** `eq:no_selection` — T′
    being pure numbers makes any sub-family's secular matrix exactly a principal submatrix,
    so Cauchy interlacing gives E(A) ≥ E(M). Measured corollary owned by C21
    `p60_best102_locked`. Backing `tests/test_paper60_no_selection.py`. The property that
    makes the encoding attractive is what supplies the bound.
13. **The general-V₀ form, and what the molecular metric IS [INTERNAL THEOREM].**
    `eq:general_v0` — V·C = V₀·B·C, every L² overlap cancelling for ANY local V₀, orthonormal
    configurations or not; atomic specialisation reproduces `eq:secular`. Backing
    `tests/test_paper60_general_v0.py`. **Provenance split that must stay intact:** the
    Shibuya–Wulfman integrals are Avery's; the identification of that matrix *as* the
    V₀-weighted overlap is ours, and the variational bound for the fixed-scale metric-free
    problem is not in his canon. Crediting either to Avery = MATERIAL.
14. **The accuracy mechanism is AVERY'S, the price in qubits is OURS [ESTABLISHED, from
    Avery].** In-out radial correlation; split-shell 1s1s′ carries two independent exponents;
    a Goscinskian 1s² pins both electrons to one exponent for *any* weighting potential.
    Presenting this mechanism as a GeoVac discovery = MATERIAL.

15. **The conditioning lever: a symbol zero of known order and location [SYMBOLIC +
    MEASURED].** The ill-conditioning is a zero of known order and location, so a
    band-Toeplitz preconditioner in Serra's sense removes it; the matching polynomial is
    EXACTLY tridiagonal (its Hankel part vanishes identically) and DST-I diagonalizable in
    closed form; `cond(G)` is FLAT in `n` against the raw `n^2` growth. Backing
    `tests/test_paper60_preconditioner.py`. **Legitimacy leg, which must stay attached:**
    any `X` with `X^T S X = I` preserves the generalized spectrum, so this is not a change
    of problem. **Scope that must stay attached:** `s`-sector shared-scale, `M = 2, 3`; and
    it does NOT recover `l`-selection — the block-diagonal congruence result
    (Löwdin / Slater-Koster 1954) is untouched. Prose implying the lever
    restores sparsity = MATERIAL.
16. **The lever reaches water's `A_1`, and the rotation is what does the work [MEASURED].**
    The symmetry-inequivalent-centre case that defeats the gerade lever (C8#9); raw growth
    reproduces `N^1.97` independently while the preconditioned column is bounded. **The
    CONTROL is load-bearing, but the UNIFORM one does NOT isolate the rotation**
    (corrected 2026-09-13): `blockdiag(P,P)` = `I2 (x) tri(1,2,1)` commutes with the
    rotation `V (x) I`, so it returns the same spectrum in either frame, and its
    exponent is `N^1.950` against the raw `N^1.967`. The discriminating control is the
    SELECTIVE `blockdiag(P,I)` in the UNROTATED frame: `N^3.79`, ending 106x WORSE
    than untreated. Reporting the gain without THAT control = MATERIAL. Backing `tests/test_paper60_preconditioner.py`.
17. **The M-centre null space is geometry-independent; its RATES are not [MEASURED,
    SCOPE].** At the symbol point every block symbol tends to `j0(0)`, so the null space is
    the constants' orthogonal complement for every arrangement — but the order at which each
    direction opens is governed by `rank(P D2 P)`, which is ONE for collinear centres, so a
    linear polyatomic opens at orders 2, 4, ..., 2(M-1). The lever is established for `M = 2`
    and NON-COLLINEAR `M = 3`; **claiming it for a linear polyatomic = MATERIAL.** Backing
    `tests/test_paper60_mcentre_orders.py`. External: Batenkov-Demanet-Goldman-Yomdin.
18. **The amplitude floor [SYMBOLIC].** `eq:amplitude_floor` — any `X` with `X^T S X = I`
    satisfies `||X|| = ||S^{-1/2}||` EXACTLY, independent of the factorization, so no
    whitening can lower the block-encoding subnormalization. Backing
    `tests/test_paper60_preconditioner.py`. Prose implying a factorization buys amplitude =
    MATERIAL.
19. **The direct block-encoding of `G` [SYMBOLIC + MEASURED].** `eq:ratio_symbol` — `G`'s
    symbol is a bounded ratio of two symbols with zeros of the same order, and its sup
    COINCIDES with `||G||`, so a circulant-embedded Toeplitz-minus-Hankel encoding carries
    `O(1)` subnormalization and the metric penalty scales as `n` rather than `n^3`. Backing
    `tests/test_paper60_direct_encoding.py`. **Two honest limits that must stay attached:**
    the circuit is CITED, not compiled (a resource model, not a gate count); and `B != G` by
    a stated operator-norm fraction that does not grow with `n`, harmless only because `B`
    enters as a whitening and the amplitude floor (C8#18) makes any such `X`
    spectrum-preserving. Presenting this as a compiled circuit = MATERIAL.
20. **Overcompleteness is ONE DIRECTION, not a property of the basis [MEASURED].** The
    Bessel deficit of a displaced Sturmian against the one-centre span PLATEAUS and does not
    tend to zero, so the one-centre set is measurably far from complete IN THE MOLECULAR
    METRIC. **The 2026-09-11 frames reading — "overcompleteness is the price of one-centre
    completeness" — is WITHDRAWN [retracted 2026-09-12: p60-frames-completeness];**
    re-asserting it = MATERIAL. Ron-Shen is the surviving
    mechanism. The paired claim form is load-bearing: the gap must collapse WHILE the deficit
    does not, which no single-sided bug satisfies. Backing
    `tests/test_paper60_one_direction.py`. **Consequent caveat that must stay attached:** the
    basis's motivating completeness is in the ATOMIC metric, not the molecular one.
21. **`eq:sigma_law` is Kac-Murdock-Szego [PRIOR ART].** Not derived here; the corpus claims
    only the IDENTIFICATION of the metric as such a finite section, Toeplitz minus Hankel
    with the stated symbol. Re-claiming the asymptotic = MATERIAL. Backing
    `tests/test_paper60_kms_attribution.py`. **Residue leg (2026-09-12):** the ~1% figure is
    DOMINATED by the `n -> n+1` grid convention, with a genuine `O(1/n)` term surviving both
    conventions; calling it purely an asymptotic tail = MATERIAL-SMALL.
22. **The `l`-selection loss is NOT a conditioning effect [SYMBOLIC].** The
    block-diagonal congruence result — a
    block-diagonal congruence cannot orthogonalize a metric that is not block diagonal, at
    EVERY `cond(S) > 1`, and it does not relax as `cond(S) -> 1+`. **Re-attributed 2026-09-12
    (C23 run #1): this is Löwdin symmetry-preservation specialized to the `l` grading, known
    since Slater-Koster (1954); what the paper claims is the `l`-vs-`m` application.**
    Presenting it as a new proposition of ours = MATERIAL
    [retracted 2026-09-12: p60-prop-d-as-new — the LABEL is retired; the
    `l`-vs-`m` application is what the paper claims]. Backing
    `tests/test_paper60_kms_attribution.py`.
23. **`eq:chirp_decay` is a Bessel asymptotic [SYMBOLIC + PRIOR ART].** Constant AND phase
    from DLMF, no stationary-phase argument needed; the `pi/4` is the BRANCH phase of the
    square-root prefactor, not a stationary-phase signature; `sum|c_j|` CONVERGES while
    `sum j|c_j|` diverges, which is the Böttcher-Widom hypothesis and a different condition.
    Backing `tests/test_paper60_kms_attribution.py`. Asserting a stationary-phase origin, or
    conflating the two sums, = MATERIAL-SMALL.
24. **Transcendental tagging of this section [SYMBOLIC + MEASURED].** Both constants are
    calibration-tier M2 and the tagging is PROVENANCE ONLY. **Prior art that must be
    credited (C23 run #3, 2026-09-12):** the Dirichlet-eigenvalue reading of the constant,
    the extremal Wirtinger-Sobolev problem behind it, and the independence of the constant
    from the rest of the symbol are ALL Böttcher-Widom's, in a source the paper already
    cites. **The Bessel-free measurement is a change of REPRESENTATION, not an independent
    route** — presenting it as independent corroboration = MATERIAL. **The removability
    corollary is WITHDRAWN [retracted 2026-09-12: p60-removability-corollary]** ("truncation-side prices are matrix-level and reachable;
    continuum-side prices are symbol-level and untouchable"): both halves are false, because
    the preconditioner is built FROM the symbol and preconditioning IS a congruence.
    Re-asserting it = MATERIAL. Backing `tests/test_paper60_contraction_window.py`.
25. **The minimiser carries the antipodal parity [MEASURED].** The band minimiser is the
    Dirichlet ground state in the band index ONLY with the alternating factor; without it
    the two are EXACTLY ORTHOGONAL, so the bare statement is not an approximation of the
    right one. Stating it bare = MATERIAL-SMALL. Backing
    `tests/test_paper60_contraction_window.py`.
26. **The law is carried by the TRANSLATION, not the metric [MEASURED + PRIOR ART].** The
    generalized symbol is the quotient, so a smooth positive radial weight cancels; the
    vanishing-weight CONTROL is the load-bearing half and reporting the agreement without it
    = MATERIAL. **Prior art (2026-09-12): Ahmad et al. give this EXACTLY for the tau/DST-I
    algebra — the structure this paper works in — so the measurement confirms a theorem
    rather than establishing one;** claiming it as novel = MATERIAL. **Scope that must stay
    attached:** this is NOT `V_0`-independence, since a position-space-local `V_0` acts by
    convolution and leaves the class. Backing `tests/test_paper60_contraction_window.py`.
27. **The translation identification is NOT ours [PRIOR ART].** Shibuya-Wulfman's own
    abstract builds the molecular operator from one unitary per nucleus; the explicit
    group formulation and the Coulomb-Sturmian translation-operator reading are both later
    published work. **What survives as ours is the SYMBOL.** Re-claiming the translation
    reading = MATERIAL. **Companion resolution that must stay intact:** Monkhorst-Jeziorski's
    "no linear dependence" and this paper's measured conditioning are the SAME pencil and
    both true — they never INVERT the overlap, and GeoVac inverts because a block-encoding
    wants a standard Hermitian eigenproblem, so the exposure belongs to the ENCODING
    REQUIREMENT. **Provenance cap:** that paper's two-page body is UNREAD (closed, no
    repository copy); the mechanism is reconstructed from the lineage and the paper says so.
    Dropping that cap = MATERIAL-SMALL.

> **Pre-registration completed 2026-09-12 (PI-confirmed) for v5.11.0–v5.11.4.** Headlines
> 15–27 cover the preconditioner lever, the water `A_1` transfer, the M-centre rate scoping,
> the amplitude floor, the direct block-encoding, the one-direction correction, and the
> five prior-art re-tierings of 2026-09-12. This ADDS claims to be checked; no existing
> goalpost was relaxed. **Standing caution carried forward:** three of these (20, 24, 26)
> record a WITHDRAWN reading, and a withdrawn reading re-surfacing is the corpus's most
> frequent defect class — C16 entries are owed for 24's removability corollary and 20's
> frames reading.

## Seeding plan (worktree only; never touches the real corpus)

K ≈ 5–6 planted defects, ≥1 catchable by each EXERCISED dimension (code / prose / citation;
C9 **is** exercised and GATING --- corrected 2026-09-07, and again 2026-09-12: this line still said "not exercised" after the file had corrected that premise twice, and the dimension then held a LARGE defect), spanning the watch-notes — e.g. an **accurate-helium** overclaim (W6), a
**"molecular many-electron 1-norm is sublinear"** reversal (branch-crit #2 / C8#10), a
**"metric-free in general"** over-generalization (W4), a **"gerade lever holds for water/
polyatomics"** reversal (W5/C8#9), a **citation splice** on `rajchel2025` or `babbush2018`
(C4/W7), a **K^0.84 → "derived/proven-optimal"** tier overclaim (C3/S4). M ≈ 5 known-good
controls (verified C8 headlines: F0=5/8, He −2.847, cond(S)∼N^1.85, a `tab:resource` row,
n_orb^2.2) that must NOT be flagged. Tiered agents (code + citation = Sonnet) get 2 seeds
each. Answer key → `debug/qa/paper_60_seed_key.json`.

## Change log
- 2026-09-13 — **FULL certifying run = FAIL (PI-invoked).** Whole-paper, unseeded, tree
  frozen. **5 of 6 dimensions PASS with zero mathematical/content defects:** deterministic
  14/14, code 136/136 `--slow` (C21 green), claims (1 SMALL), citations (49/49 resolve),
  synthesis C9 gating (1 NIT). **Completeness = FAIL** on three self-declared
  cert-blockers, not on any hidden defect (zero undiscovered gaps, zero surviving stale
  echoes). Remediated this run: F1 (abstract n^3->n clause re-tiered [RESOURCE MODEL]),
  the C9 eigenvalue-wording NIT, the Monkhorst-Jeziorski bibitem title, and — closing a
  cert-blocker — the Sylvester-inertia bound C8.16 now has a standalone fire-tested test
  (`test_c5_inertia_bound_and_root_by_root`), matrix row 568 OPEN->BACKED. **Remaining
  gating work: register the ~7 declared-debt literals (itemized in CHANGELOG v5.11.10) as
  their own careful pass, then a clean delta.** NOT certified.
- 2026-09-13 — **DELTA-verification #3 = DEFECTS, remediated. Cleanest of the lineage.**
  Four dimensions; **citations and code both CLEAN-DELTA**; 0 LARGE, 0 mathematical.
  Three SMALL, all the v5.11.8 water-control fix (rotation is load-bearing) not reaching
  two summary surfaces plus a stale number in the corrected body sentence: paper L1338
  (`2766` was the N=48 interior point, not the low endpoint — measured 194 at the paper's
  own grid; **the claims reviewer's proposed 729.2 was itself wrong and was caught before
  printing**), `papers/INDEX.md` (omitted rotation caveat), `docs/walls/register.md` (cited
  the blind uniform control as discriminating evidence). New category: paraphrase-level
  requirement-omission survives token gates even in files C16 was widened to. **Process:
  the tree was frozen and no edit made until all four reviewers returned — DELTA #2's
  moving-target defect did not recur.** Public-web layer stale but a known PI-gated gap,
  unchanged. **Next: a clean delta before any FULL certifying run.**
- 2026-09-13 — **PROCESS DEFECT, the PM's: THE TARGET MOVED UNDER THE REVIEWERS.**
  The code reviewer finished its mandated run at 08:25, then the PM applied
  `debug/delta2_fix_01..06_*.py` between 08:30 and 08:43 — rewriting the paper's
  resource-control paragraph, the claim matrix, the claims register, the walls
  register, `papers/INDEX.md`, the group2 synthesis, and
  `tests/test_paper60_preconditioner.py`, which gained a 13th test that did not
  exist when that run collected its 52 items. **A DELTA verdict formed on a
  moving tree certifies a state that no longer exists.** The `/qa` protocol
  already says how to avoid this — *"snapshot it into a worktree first if the PM
  will keep editing during the run"* — and the PM did not. The reviewer caught it
  itself, re-ran against the new state (12 passed / 1 skipped) and independently
  reproduced the new numbers (raw N^1.9666, uniform N^1.9497, selective N^3.7913;
  106.4x at N=192), so the CONTENT is sound and the PROCESS is the finding.
  **Standing rule for the next run: snapshot before dispatch, or do not remediate
  until every reviewer has returned.**
- 2026-09-13 — **DELTA #2 — two NEW affected categories, and a remediation that
  never landed.** (1) **The applier class.** `delta_fix_04_last_three.py` returns
  from inside its edit loop on a stale anchor, BEFORE its write loop, so one bad
  anchor discards every edit in the script silently — while CHANGELOG v5.11.7
  recorded the remediation as applied. Three fixes (the walls register's retired "Proposition D" label
  [retracted 2026-09-12: p60-prop-d-as-new] at three loci, the
  floor bracket in the synthesis and claims register) were written 2026-09-13 and
  were not in the corpus. `delta_fix_03` had failed the same way one step earlier
  and said so in its own docstring. Appliers now WRITE FIRST and report misses.
  (2) **The public web-rendering layer** (`viz/public/papers/*.html`,
  `index.html`, `sitemap.xml`) — the outward-facing sibling of DELTA #1's
  generated-artifact LARGE. The group2 synthesis renders an abstract saying "the
  nine Group 2 papers" where the source says twelve, and **Papers 58, 59, 60 and
  61 have no public page at all**, so none of the nine changed claims has any
  public rendering. Built by `debug/build_paper_pages.py` from a stale
  `debug/data/zenodo_manifest.json`. **DECLARED GAP, not fixed here:** rebuilding
  the manifest touches the Zenodo/DOI distribution surface, which is PI-gated.
  (3) **`docs/walls/register.md` is unreachable by every gate** — the string
  `docs/walls` appears nowhere in `debug/qa/`, yet the file carries an operative
  dispatch rule. Same "operative instructions outside the gate" shape as the
  recall layer that was DELTA #1's second LARGE, but inside the repo. Closed by
  widening the six 2026-09-12 C16 entries, which had shipped with a NARROWER
  `files` list than the older p60 entries beside them; widening took C16 from
  PASS to FAIL with 7 live loci, all now closed, and the widened entry is proven
  to discriminate both ways.
- 2026-09-13 — **DELTA-verification = DEFECTS, remediated. NOT a clean delta.** Four
  dimensions over the claim-impact set, unseeded; deterministic 14/14; 172 tests pass;
  **zero mathematical defects**. 2 LARGE, both OUTSIDE the nominated set by category —
  the generated-artifact layer (prior-art credit in the output but not the generator, so
  the distributed JSON lacked it and regeneration would erase the markdown copy) and the
  auto-loading recall layer (four stale readings, one an operative instruction to claim
  what was given away). 14 SMALL across claims/citations/claim-impact, plus one MATERIAL
  code finding (the new full-shell test cited by no row). **Two self-inflicted errors
  recorded:** the 2026-09-12 chain fix added "values are truncated, not rounded", false
  against the exact 1s² value; and the PM accepted a reviewer's under-sampled convergence
  finding and re-tiered a CORRECT claim, which the reviewer then withdrew and an extended
  ladder refuted (c=3 plateaus at ~1.6e-7 by 250k points). **Next: a clean delta before
  any FULL certifying run.**
- 2026-09-12 — **FULL run (PI-invoked) = FAIL, remediated. NOT certified.** Six dimensions,
  unseeded. C8 extended to headlines 15–27 first (PI-confirmed) so v5.11.0–v5.11.4 were
  pre-registered rather than unmeasured. **Zero mathematical defects**; every published
  converged constant independently reproduced by both code reviewers; 46/46 bibitems resolve.
  Five load-bearing findings: a retired `-2.897` He chain endpoint in the abstract (violating
  the paper's own interlacing bound, now `-2.896` and registered); abstract + conclusion
  pre-breach against `sec:resource`; the water control provably blind (it commutes with the
  rotation it names); the τ prior-art surrender over-scoped (our matrices are not in that
  algebra — measured 1.1% off-diagonal at n=160); and an `x == x` guard behind an abstract
  headline. Two reviewer findings OVERTURNED by PM verification (a "missing" 1967 reference
  exists; a wrong-basis-point finding did not apply to the abstract). Instrument findings:
  three gates scoped to `trunk` on first invocation; C16 clean on files it scopes while
  zombies live in them; a claim-matrix row citing a deleted test; the staleness banner
  measuring the paper and not the synthesis on a GATING dimension; C21 blind to the molecular
  half. Remediated in six passes, content before guards; four new guard assertions all
  fire-tested. **Next: a delta-verification run over the remediation.**
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
