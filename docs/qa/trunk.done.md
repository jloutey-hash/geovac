<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED — AND NOT CLOSE
> The CERTIFIED verdict below is **historical**. Do not cite it as
> present-tense status.
>
> **Run record (all adverse, all remediated):**
>
> | run | date | verdict |
> |---|---|---|
> | FULL #1 | 2026-09-01 | FAIL (calibrated 13/13, 0 FP) |
> | FULL #2 | 2026-09-01/02 | FAIL (19/21 seeds, 0/8 FP) |
> | DELTA #1 | 2026-09-02 | DEFECTS (9/9 seeds, 0/8 FP) |
> | DELTA #2 | 2026-09-02 | DEFECTS (8/9, one void seed) |
> | FULL #3 | 2026-09-02/03 | FAIL (unseeded; 55 MATERIAL) |
> | DELTA #3 | 2026-09-03 | DEFECTS |
> | DELTA #4 | 2026-09-03 | DEFECTS |
> | FULL #4 | 2026-09-03 | FAIL (4 dimensions, 13 gates) |
> | DELTA #5 | 2026-09-03/04 | DEFECTS |
> | DELTA #6 | 2026-09-04 | DEFECTS |
> | DELTA #7 | 2026-09-04 | DEFECTS (40 MATERIAL, 10 LARGE; P39 re-priced CONDITIONAL) |
> | DELTA #8 | 2026-09-05 | DEFECTS (5 MATERIAL, 1 LARGE; mirror-staleness) |
> | FULL #5 | 2026-09-05 | FAIL (seeded re-measure; 10/12 seeds, 0 FP; 9 MATERIAL-SMALL) |
> | DELTA #9 | 2026-09-05 | CLEAN-DELTA (v5.9.5 remediation material-clean; 3 NITs swept) |
> | FULL #6 | 2026-09-05 | FAIL (certifying; 2 LARGE + 6 SMALL) -> remediated v5.10.0 (claim-scope, gate-first) |
> | DELTA #10 | 2026-09-05 | DEFECTS (first claim-impact run) -> trunk items fixed v5.10.1; group1 cluster owed |
> | FULL #7 | 2026-09-06 | FAIL (certifying; 4 SMALL, 0 LARGE; all keystones re-derived SOUND) -> content remediated same day; guard pass + C21-on-trunk owed |
> | DELTA #11 | 2026-09-06 | CLEAN-DELTA (FULL #7 remediation + guard pass + C21 closure all verified; 3/3 dimensions clean) -> unlocks FULL #8 |
> | FULL #8 | 2026-09-06 | FAIL (certifying; 2 SMALL, 0 LARGE; fresh adversary caught 2 pre-existing P1 s/p paraphrase zombies FULL #7 missed) -> remediated + C16-backstopped same day; C21 symbolic-fraction constants owed |
> | DELTA #12 | 2026-09-06 | DEFECTS (claims CLEAN; C21 symbolic instrument built but its tight-match was inert -> guard fixed + fire-tested) -> FULL #8 remediation + C21-on-trunk verified |
> | DELTA #13 (claim-impact) | 2026-09-27/28 | DEFECTS (~19 LARGE, one family: graph/κ credited with the hydrogen spectrum) -> v5.17.1 |
> | re-checks #1–#3 | 2026-09-28 | DEFECTS each (residue paraphrases, then ~13 LARGE pre-existing debt) -> v5.17.2, list-closing pass v5.17.3 |
> | FULL #9 | 2026-09-28 | FAIL (unseeded; ~15 LARGE + ~35 SMALL; RESEARCH-GRADE NOT earned). Keystone mathematics SOUND. LARGE cluster: P32's diagonal finite triple ([D,π(f)]≡0) and its consequences; finite four-witness "literal identification" trivial (K=2J_z); order-zero/one residuals grow; E6 proof; P38 substrate A1/A2; P7 Stage-3 tier; P1 title/α. PI approved all recommendations; remediation in progress |
>
> **No run has produced PASS.** Only a FULL run can, and the last four FULL
> runs each found material defects — several *introduced by the previous
> run's own remediation*, which is why a clean DELTA is the precondition for
> firing another FULL.
>
> Remediation scope: [`trunk.carryforward.md`](trunk.carryforward.md),
> Parts F–O. Open at Part N.5: **the guard pass**, deliberately deferred as
> its own separately-reviewed unit (CLAUDE.md §9) with six carried items.
>
> Spill outside trunk: group1 (Papers 38/39/40) and group3 were both edited
> by trunk remediation, so their own records are staler than their banners say.
>
> *Re-measure rather than trusting this banner — it is a snapshot and has
> already gone stale once, stopping at FULL #2 while eight further runs
> landed:*
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Trunk — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only trunk-specific scope + deltas. C1–C13, the verdict rule,
> the review-dimensions map, and the hard rules live in `criteria.md`.

> **STATUS: NOT CERTIFIED; RESEARCH-GRADE NOT earned** (FULL #9, 2026-09-28 = FAIL).
> The historical "certified PASS" (run #4, 8/8 seeds, v4.16.2, 2026-06) is superseded
> by every run since; see the banner.

**Scope:** Papers **0, 1, 7** (group3 foundations) + **32, 38** (group1
operator-algebras) + the **group3 foundations synthesis**. The trunk is the
foundation every branch depends on; it is QA'd first so a finding at a root
re-prices everything above it.

**Deterministic `--gate`:** `trunk` (the default scope for the check scripts).

## Branch deltas (the only non-inherited content)

- **C7 (trunk-dependent status).** The trunk *is* the foundation, so C7 reduces
  to WH1 self-consistency: Paper 38 / WH1 is PROVEN **scoped to the van
  Suijlekom state-space GH distance** (translation-seminorm metrization), with
  no residual "Latrémolière propinquity" overclaim where the proved object is
  the state-space GH distance.
- **C8 (headline honesty), per-paper.** κ = −1/16 is an **Observation**
  (coincidence, no bridge), not a derivation; **4/π** is **derived
  (numerics-pinned)**, not asserted as full symbolic proof; the **Forced-Count**
  moduli chain is stated at its **full-axiom** count, not the matter-sector
  subcount.
- **C6/C8 — the graph→S³ re-pricing (added 2026-09-27, pre-registered for the
  post-v5.16.0 DELTA).** The edge Laplacian D−A of the (n,l,m) lattice does
  **NOT** converge to Δ_{S³}. It preserves l, so it splits into n_max disconnected
  blocks with n_max zero modes, and its spectrum is algebraic and bounded in [0,8]
  (`tests/test_graph_s3_nonconvergence.py`). The positive anchor is the exact
  label set, with the S³ Laplace–Beltrami operator acting as the **diagonal SO(4)
  Casimir** n²−1, and Fock's continuum projection carrying the equivalence.
  **π-free** means the graph's entries are rational, so its eigenvalues are
  *algebraic*, **not** "integers or rationals". The pure-integer ladder −(n²−1) is
  the continuum spectrum. κ = −1/16 maps the extremal eigenvalue (a band edge),
  not the spectrum.
  - **LARGE:** any live locus asserting any of these —
    - that the graph Laplacian converges to, reproduces, or equals Δ_{S³} or the
      Rydberg ladder;
    - an integer or rational graph spectrum;
    - that `tests/test_paper7_graph_convergence.py` backs such a convergence (it
      measures λ_max → 2 d_max only).
    The C16 entries `graph-laplacian-converges-to-s3` and
    `graph-eigenvalues-integer-or-rational` are the deterministic backstop; the
    claim-impact reviewer covers paraphrases.
  - **Also in claim-impact scope:** WH3's graph-limit leg, and Coulomb/HO Layer 1
    (withdrawn: 7 catalogued, 6 standing), wherever trunk papers restate them.
- **C7 extension — which object converges (pre-registered 2026-09-28 after FULL #9,
  PI-approved; applies from the next run).** Paper 38's convergence theorem is about
  the **true Camporesi–Higuchi truncation** (Dirac eigenspaces up to the cutoff), not
  the doubled Fock spin-orbital space. **LARGE:** any live locus that describes the
  convergent object inconsistently with the proof, or that motivates the translation
  seminorm by a Dirac-seminorm degeneracy that does not hold on the true window.
- **C8 extension — Paper 32's finite axiom table (pre-registered 2026-09-28 after
  FULL #9, PI-approved; applies from the next run).** The finite triple's audited
  status is a headline. **LARGE:** any live locus asserting that the diagonal
  spectral form D_GV carries a metric or a nonzero commutator with π(A); that the
  finite truncation is odd (it is even via chirality doubling; no KO label at finite
  cutoff); that the order-zero/one failures decay; that the finite four-witness 2π
  period is a literal identification of a Bisognano–Wichmann modular flow (it is a
  trivial symbolic fact, K = 2J_z); or the E6 "no single mechanism" statement at
  any tier above Observation. C16 backstops: `diagonal-commutator-lipschitz`,
  `truncation-inherits-odd`, `order-zero-one-decaying-artifact`,
  `finite-four-witness-literal-identification`, `e6-no-single-mechanism-theorem`.
- **No branch-specific C14+.**

## Change log
- 2026-09-28 (v5.17.6 @ 83b731c, DELTA #15, claim-impact, unseeded) — **DEFECTS;
  RESEARCH-GRADE NOT earned.** ~13 LARGE + ~40 SMALL (down from ~16 at DELTA #14),
  every LARGE PM-verified. Mathematics all CONFIRMED: P16 Theorem 2 proof SOUND (a third
  independent route, dual Jacobi-Trudi); the true-window oddness, J relations and
  order-zero values reproduced by a third construction (polynomials on R^4, no Wigner/CG
  code); P42 §11 merge proofs SOUND; citations 0 LARGE, no fabricated IDs. Dimensions:
  deterministic PASS; citations 0 LARGE; code 1 LARGE (provenance: the second "route" for
  the true-window oddness was cited to a test that does not contain it -- PM's own
  error); claims P32/38 3 LARGE + 1 in doubt (unscoped "finite truncation is even" x6
  after the true window was found odd -- PM's own propagation miss; thm:GV_triple scope
  clause garbled; retracted cb-norm 2/(n+1) live; P38 §marcolli_vs gauge-network
  anchor); claims P0/1/7+synthesis 2 LARGE (synthesis §s3 pre-remediation Stage-3
  account; P7 calls the Coulomb kernel the Green's function of Delta_{S^3}, contradicted
  by its own 2pi^2/n eigenvalue); claim-impact 6 LARGE (P24 'D_W off-diagonal' after
  P42 corrected it; P45 Q1 conditionality not mirrored in P42/P34; P50 CHM geometry
  CFT3 vs CFT4 (S^2 vs S^1 entangling surface) undercutting headline (B); P42 test
  attribution 13 cells vs 7; group1 synthesis LvS O(1/n) rate; P16 novelty over-credited
  vs Rittenhouse Thm 1 (any non-interacting N) and the K-harmonics Young-scheme
  literature). Dominant class again: a correction not yet propagated to every locus
  that restates it -- now largely inside the remediation's own new content.
- 2026-09-28 (v5.17.5 working tree, DELTA #14, claim-impact, unseeded) — **DEFECTS;
  RESEARCH-GRADE NOT earned.** ~16 LARGE + ~43 SMALL, every LARGE PM-verified in the
  frozen worktree. Dimensions: deterministic PASS (all 14 gates); citations CLEAN of
  LARGE (0 LARGE / 5 SMALL); code — P7/P16/classifier CLEAN (0 LARGE, all numbers
  re-derived by a second route, 10/10 plants fired), guard pass 1 LARGE (the rewritten
  s/p splitting guard cannot see a wrong cluster tolerance: DEG_TOL 1e-8 -> 1e-16 left
  19/19 green, PM-reproduced); claims P0/1/7 1 LARGE family; claims P32/38 5 LARGE + 1
  pre-existing candidate; synthesis 2 LARGE; claim-impact 7 LARGE. Dominant class: a
  correction applied in the owner (P7 Stage 3 kernel, P38 true window / A2, P32
  sections, P16, P56) not yet propagated to every locus that summarizes it — including
  the mirror direction (owner STRENGTHENED, summaries kept "imposed"). One LARGE is a
  remediation error of this cycle (P50's CHM "both spectra integer-valued" is false for
  the continuum modular Hamiltonian). List brought to the PI per the standing agreement.
- 2026-09-28 (v5.17.3, FULL #9) — **FAIL; RESEARCH-GRADE NOT earned.** Unseeded;
  standing calibration record cited. Dimensions: deterministic PASS (C5, C10–C24);
  code FAIL (P7 1 LARGE, P32 3 LARGE, P38 2 LARGE; P0/P1 sound); claims FAIL
  (P1 title/α/"unexplored"; P32 gauge/Lorentzian); citations FAIL (Leimbach–van
  Suijlekom credited with a rate they do not give); completeness FAIL (P32 prop:odd
  contradicts its own measured grading; 14 SMALL; this DoD's status line was stale).
  Every LARGE verified by the PM against text or code. Keystone mathematics SOUND
  (P38 theorem and constants, L5 refutation; P32 Forced-Count 32 by two routes, sign
  triple, prop = 2; P7 wall, Casimir, λ_max, F⁰). PI approved every recommendation
  (P32 rescope, P38 A2 independent build first, P7 kernel test, E6 → Observation,
  P1 retitle, the two DoD extensions above). Instrument fix: `fire_test.py` now
  takes a per-target lock (two concurrent reviewers had corrupted a target).
- 2026-09-27/28 (v5.17.0 → v5.17.1) — **claim-impact DELTA (unseeded) = DEFECTS;
  RESEARCH-GRADE NOT earned** (about 19 LARGE, about 16 SMALL). Dimensions:
  - deterministic: green, after a C20 fix;
  - code (Opus): CLEAN;
  - citations (Opus): CLEAN;
  - claims: 7 LARGE;
  - synthesis: 2 LARGE;
  - claim-impact: 12 LARGE families.

  Every LARGE was one family: the graph or κ credited with the hydrogen spectrum.
  A v4.76.0 audit confirmed the refutation (the model has been unchanged since v4).
  Remediated gate-first:
  - two C16 entries (`kappa-maps-graph-spectrum`, `graph-hamiltonian-carries-spectrum`);
  - about 45 papers edited;
  - Paper 7 retitled;
  - the AtomicSolver default changed to the physical H.
  **A claim-impact re-check is owed** before RESEARCH-GRADE can be stamped. Its
  recently-changed loci are those listed in CHANGELOG v5.17.1. Standing calibration
  record cited.
- 2026-09-27 (v5.17.0) — pre-registration for the post-v5.16.0 DELTA: the graph→S³
  re-pricing added to C6/C8 as a LARGE class (above). RESEARCH-GRADE status
  available (criteria.md); this target is OWED.
- 2026-06-14 — created (co-authored PM + PI) as the first pre-registered `/qa`
  target.
- 2026-06-14 — added the all-dimensions-mandatory rule after run #3 found
  runs #1–2 had exercised only the claims + citation dimensions; the code
  (C1–C2) and synthesis (C9) dimensions surfaced 2 real synthesis defects
  (κ "derivable", K "conjectural"). Run #4 = **PASS** (8/8 seeds).
- 2026-06-16 — **slimmed to a profile**; C1–C13 + rules moved verbatim to
  `docs/qa/criteria.md` (no criterion changed). PASS status unaffected.
- 2026-09-01 — FULL run #1 = **FAIL** (see banner). C3 inline-tier pass,
  reverse-citation pass, and the C11 could-not-fail fix followed the same day
  (run record §§1c–1g).
- 2026-09-02 — FULL run #2 (dispatched 2026-09-01, PI-directed FULL rather
  than the DELTA the carryforward asked for) = **FAIL**. Criteria frozen as at
  `dbe7ae2`; nothing added or relaxed. Two C16 registry entries added with
  two-way discrimination proofs (`latremoliere-propinquity-named-for-gh-rate`,
  `p45-kplus-compression-theorem-live`). Code dimension certified only on the
  seed classes each Sonnet agent demonstrably caught; delta-run seeding
  requirements in the carryforward.
- 2026-09-02 (v5.3.0) — the five PI items from run #2 executed (Part F
  F1.1/F1.2/F2.1/F2.2/F6). **Gate change (PI-authorized, criteria.md +
  qa.md):** CODE dimension → Opus on the trunk roots and on any file whose
  tests assert a convergence endpoint; three seeds per Sonnet code agent
  elsewhere. Three more C16 entries, two-way proven (registry 30 → 33). The
  FAIL stands; the rest of Part F and the DELTA are still owed.
- 2026-09-02 (v5.3.1) — Part F remediation complete (F1.3–F4.6, 30 rows,
  each with its DONE note) plus seven findings the remediation surfaced
  (F7.10–F7.16: Paper 40 semisimple scope, the dim_H = g₃ coincidence at
  five loci, c²(4,3), wall-entry arithmetic, two order-of-magnitude
  restatements, a bibitem-key collision, and the widened C16 entry finding
  eight more loci in Paper 18 + the group1 synthesis). All deterministic
  gates PASS on `trunk` and the `group1` subset. The FAIL stands until the
  DELTA run under the F-seeding rule returns CLEAN-DELTA.
- 2026-09-02 (v5.3.1, later) — DELTA run #1 under the F-seeding rule:
  calibrated 9/9 (0/8 FP), verdict **DEFECTS** (26 genuine items in the
  Part F remediation; carryforward Part G), remediated the same day; all
  gates green again. DELTA #2 (fresh seeds) is the precondition for the
  FULL certifying run. The FAIL stands.
- 2026-09-02 (v5.3.1, DELTA #2, PI-invoked) — scope = DELTA #1's
  remediation; calibrated 8/9 (one void seed, CLAIMS-1 INCONCLUSIVE), 0/8
  FP; verdict **DEFECTS** (15 items, carryforward Part H), remediated the
  same day; two PI adjudications (Paper 40 title; finite-triple KO label,
  measured (−,+,+) with the production grading). DELTA #3 owed. FAIL stands.
- 2026-09-02 (v5.4.0, PI minor bump) — gate change: blind seeding is
  OPT-IN (`/qa <target> seeded`); default runs are unseeded and cite the
  standing calibration record. PI adjudications applied: Paper 40 retitled
  (semisimple), finite-triple KO label dropped for the measured sign triple
  (−,+,+), Paper 7 NO-TEST row closed by `tests/test_paper7_graph_convergence.py`,
  Paper 40 §L5 reframed as state-space GH. All deterministic gates green on
  trunk + group1 (C10 trunk + P40, C11, C13, C14, C15, C16, C17, C18, C19,
  C20, C21, C22, C5). DELTA #3 skipped by PI direction; FULL run #3
  (unseeded) is next. FAIL stands until it reports.
- 2026-09-03 (v5.4.0 -> remediation) — **FULL run #3 (unseeded) = FAIL**: 55
  verified MATERIAL findings across all LLM dimensions (code x6, claims x3,
  citations P32, synthesis); deterministic 13/13 PASS but two C16 instrument
  defects found. Five re-pricings raised to the PI (carryforward Part I.0):
  the graph->S^3 convergence claim measures lambda_max -> 2 d_max; Paper 38's
  rate constant is 2/pi in the paper's own unit-S^3 metric (4/pi = rotation-
  angle metric); the Forced-count endpoint 260 is a degenerate-sample artefact
  (correct: 32); Paper 32 Theorem 1 names the wrong algebra for the limit;
  prop = 2 is generic. Remediation in progress; FAIL stands.
- 2026-09-03 (v5.4.1) — FULL run #3 remediated (carryforward Part I.8): five
  re-pricings applied with dated correction notes, one production bug fixed
  (SM representation), three new tests + six rewritten, C16 widened. FAIL
  stands until an unseeded DELTA and FULL run #4.
- 2026-09-03 (v5.4.2) — **DELTA #3 (unseeded) = DEFECTS**, remediated same day
  (carryforward Part J): the FULL-run-#3 remediation had been applied
  locus-by-locus, so the descoped readings survived at twelve loci; swept
  claim-wide. Two new findings: the s/p splitting is a node-amplitude proxy,
  not a spectral gap (disconnected l-blocks, near-tie mode selection); and
  the block spectrum is closed form with a proven O(n_max^-2) saturation rate
  (upgrade). Citations CLEAN-DELTA. FAIL stands; DELTA #4 owed before FULL #4.
- 2026-09-03 (v5.4.3) — self-caught: the L5 panel check added in v5.4.1 was
  itself wrong (the panel height rises toward 1 while gamma falls, so the
  inequality holds only for n_max <= 5). Reframed as a measurement in the
  test, the module, Paper 32 and a new panel note in Paper 38; the panel-side
  quantity is a named open check. C10 trunk PASS, C19 clean.
- 2026-09-03 (v5.4.4) — the L5 panel crossing is now MEASURED (margin +0.053
  at n_max = 6, −0.069 at 7, so the crossing is between 6 and 7); Papers 32/38
  and the test carry the measured table in place of the extrapolated
  "near n_max ≈ 6", and a slow test pins the crossing.
- 2026-09-03 (v5.5.0, PI direction) — the Paper 40 normalisation cross-check
  is SETTLED: the dual-Coxeter rule Cas(ad) = h∨ on su(2) gives the rotation
  angle as its geodesic distance, so Paper 38's moment is a dual-Coxeter
  moment and 4/π is the canonical value. Paper 38 restated on that sphere
  (Option 1); Paper 40's named check closed at rank 1. Consequence, after
  investigation: the constant is the unit-sphere quotient Vol(S²)/Vol(S³) =
  2/π, so M1's volume content and period ring stand; only the "Hopf base"
  label is a misnomer (that base carries half the radius). No sweeping change.
  **The word "canonical" in this entry was itself wrong — see the v5.5.2 line.**
- 2026-09-03 (v5.5.2) — **DELTA #4 (unseeded) = DEFECTS**, remediated same day
  (carryforward Part L). Two reviewers, dispatched blind on the same diff,
  independently produced the same two LARGE findings. (i) "4/π is canonical"
  is an over-claim: Cas(ad) = h∨ is the corpus's declared rule, not the
  field-standard one — Kac's basic form (θ|θ) = 2 gives Cas(ad) = 2h∨, radius
  √2 and 2√2/π. What survives is that P38 and P40 use the *same* rule, so
  the rank-1 cross-check is genuinely settled; what fails is the word
  "canonical" and Paper 40's attribution of its rule to the literature.
  (ii) The Hopf-base retraction had again been applied locus-by-locus and
  missed ~10 loci, including a theorem body; swept claim-wide and registered
  in C16 so the class is now mechanical. Three new guards could not fail
  (h∨ a literal; the radius test an identity on its own literals; the γ ratio
  pin kernel-blind) — all three rebuilt with routes the claim does not supply,
  and each re-planted to prove it fires. Also closed: the §13.4 equation gate
  on eq:seminorm_normalisation, and the He⁺ benchmark row that asserted a
  check its own constructor imposes. Twelve deterministic gates PASS.
  FAIL stands; the DELTA is now clean, so FULL run #4 is unlocked.
- 2026-09-03 (v5.6.0) — **FULL run #4 = FAIL**, remediated same day in four
  tranches (carryforward Part M). Eleven reviewers, all four dimensions, 13
  deterministic gates. Headline: **two gates were examining nothing on trunk**
  — C17 had zero trunk families and C21 zero `\gvq` annotations — so seven
  headline re-pricings this cycle all went to C16, a blocklist that cannot
  catch a correct value drifting to a second locus. Two published readings
  withdrawn, both PI-approved: the s/p evidence leg is an artifact of cutoff
  divisibility (λ₂ₛ = 3 exactly iff n_max ≡ 0 mod 3; the 0.39% endpoint sits on
  the favourable branch), and Paper 38's L5 height bound is FALSE by a
  finite-band witness (height_B ≡ 1, so it fails for every n_max ≥ 6).
  *Corrected 2026-09-04:* this line said “exactly the v5.4.4 crossing”, conflating two different events — γ falls below 1 between n_max 5 and 6, while the PANEL rises past γ between 6 and 7. **WH1 and
  `thm:main` are untouched**; the unconditional theorem uses no height
  constituent. Both withdrawals are also upgrades — each replaced an asserted
  trend with a proved closed form. Three papers named operators they do not
  use, and in all three the code was right and said so in a docstring. Five
  guards that could not fail, two of them written the previous day. FAIL
  stands; a DELTA on the remediation is owed before FULL #5.
