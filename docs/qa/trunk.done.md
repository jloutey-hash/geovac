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
> | FULL #6 | 2026-09-05 | FAIL (certifying; 2 LARGE + 6 SMALL, all prose-zombie class; math sound) |
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

> **STATUS: FROZEN — certified PASS** (`/qa trunk` run #4, 8/8 seeds, v4.16.2).

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
- **No branch-specific C14+.**

## Change log
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
