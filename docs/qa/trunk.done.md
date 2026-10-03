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
  chirality-doubled MODEL truncation is odd (it is even via chirality doubling), or that the TRUE window (the converging object) is even (it is odd: commutant of O = scalars; rescoped 2026-09-29 after DELTA #16, PI-approved); no KO label at finite
  cutoff); that the order-zero/one failures decay; that the finite four-witness 2π
  period is a literal identification of a Bisognano–Wichmann modular flow (it is a
  trivial symbolic fact, K = 2J_z); or the E6 "no single mechanism" statement at
  any tier above Observation. C16 backstops: `diagonal-commutator-lipschitz`,
  `truncation-inherits-odd`, `order-zero-one-decaying-artifact`,
  `finite-four-witness-literal-identification`, `e6-no-single-mechanism-theorem`.
- **No branch-specific C14+.**

## Change log
- 2026-10-02 (v5.18.1 @ 6c4d162, DELTA #25, claim-search DELTA, unseeded; 4 owner-in-full + 4 claim-search + code + citations, all Opus) — **DEFECTS; RESEARCH-GRADE NOT earned.** ~10 OWNER-LARGE (P38 Berezin pair live inside refuted L5 + M1 placement as fact; P32 Yukawa 'forced/structural reason' + 'M1 owns GH rates'; P18/P55 'M1 produces the 4/pi rate'; P55 'all GeoVac transcendentals' contradicted by P61's elliptic period; P40 C3(G)*gamma display + 'closes it positively'; P57 false equation Delta = c^2(3,2)), ~15 SURFACE-LARGE (README K 'conjectural' (section 13.5) + 'eight theorem-grade'; P18 Three-homes clause; group1 synthesis unsynced since the P38 re-pricing; group3 synthesis ladder tier; field guide; cited docs/ catalogues; dated docs), 3 instrument (C25 comma-list escape on the pi-source keystone; dead p40 overclaim; README/docs off C25's front door) + C20 natbib-blind and silent re-baseline. Citations clean bar 2 SMALL. Remediated (PI-approved incl. A: archive pitch + splinter, banners, re-sync forced_free_seam; B: P1 ladder -> [SYMBOLIC PROOF]; C: P55 'census'): 18 idea-level C16 entries (198), C25 hardened + 2 claims (19) + every claim probed, C20 natbib + refusing re-baseline, C17 README scope. Two pre-checkpoint self-checks: 10 then 6 LARGE siblings (incl. README's whole retired quantitative surface, P2 circulant 'explains', depth = loop order, P39 Gaudillot-Estrada credit, S_min 'open' contradiction introduced by wave 2), all fixed. Guards fire-tested. All gates PASS; C10 whole corpus. **Self-check loop to zero LARGE (new rule), 15 rounds: 10, 6, 9, 4, 7, 6, 5, 3, 3, 2, 2, 3, 1, 3, 0** -- headline results: Paper 39 Lemma L3-T upgraded to C_3^(2) = 1 by compression on the Clifford-doubled module (k-fold and master forms with constant 1, rate the joint moment / sum); Paper 40 per-summand bridge withdrawn and (INT) scoped by summand class; Paper 45 Fejer rate O(T log N_t / N_t) and panel / Riemannian-limit on the spatial leg; propinquity credit removed from Papers 32/39/42 and three syntheses; code layer re-labelled (max(gamma) a reference value). C16 at 226 entries. Debt: docs/qa/trunk.debt.md.
- 2026-10-02 (v5.18.1, trunk debt pass, PI-directed) — all five debt items closed (P34 longtable; P32 19 attributions verified + cited, CCM year fixed; Kac normalisation computed; P7 novelty withdrawn as Bander-Itzykson prior art; P32 slot-vs-ring clarified). Instrument defect: C20 skipped 52 of 57 unbaselined papers; now fails on them and the baseline covers all 57 (147 corpus tags recorded). a fresh Opus reviewer found 0 LARGE, 1 SMALL, 5 NIT, all fixed: the new slot/ring sentence glossed the slot as "which power of D enters the trace", which would file its own M2 example (zeta_{D^2}, a D^0 trace) under k = 0 (now "the Mellin slot of the master engine"); a GeoVac inductive system "hosts" the pro-unipotent factor (now "would host"); "Tannakian theorem" for the Hopf-algebra/group-scheme anti-equivalence (renamed, Waterhouse 1.4 restored); the Adams cite placed on the H-space clause only; Paper 7 item retitled "organizing use only". It also flagged Paper 18's master-Mellin list still assigning M1 the domain "state-space GH rates" (now volume-slot normalizations, with the GH-rate constant sharing the value as a reading). Adams, Bander-Itzykson, the longtable, and the Kac test were confirmed correct. All gates PASS; C10 whole corpus.
- 2026-10-02 (v5.18.0 @ a5bab2e, DELTA #24, first claim-search DELTA, unseeded; 3 owner-in-full + 4 claim-search + code + citations, all Opus) — **DEFECTS; RESEARCH-GRADE NOT earned.** 5 OWNER-LARGE (P1 L53 '[PANEL-VERIFIED] bit-exactly' after the PM re-tiered every other locus to [MEASURED] 1e-9; P1 'higher-l states have higher connectivity' stated as verified, false on both constructions; P32 Door-4d 'closes the n=2 fork by selecting H'; P32 Forced-Count 'the framework forces the dimension' (borderline, PM call OWNER); P2 'Three homes' Theorem asserting 'no spectral construction generates all three', the negative P32 withdrew), ~20 SURFACE-LARGE (group1 synthesis and field guide crediting P38 with the almost-commutative/gauge-network convergence + retired C3*gamma display + 'holds for all n'; group3 synthesis 'M1 governs the GH rate'; P55/P18/field guide 'complete stratification'; group5 synthesis 'largely forced'; P57 abstract factor count 'autonomously generated'; F^k 'on S^3' at 7 loci; claims-register row 6 He 0.004% as variational (caught by C25 pre-run); P25 'This is a prediction' about K; P38->P40 'resolved'), 2 instrument LARGE in C25 itself (owner checks that could not fail, incl. the K hard prohibition; qualifiers satisfied by their negations: unconditional/almost-commutative/half-integer), ~45 SMALL. Citations 10/12 CONFIRMED (CvS 'outlook' misattribution SMALL; Kac 2h-dual UNVERIFIABLE -> debt). Remediated (PI-approved, incl. the P2 Three-homes clause moved out of the Theorem): C25 hardened (adjacency-only negation, arrow exemption only toward the current tier, word-bounded qualifiers, brace-balanced \gvc parser, macro normalisation, specific owner anchors + owner-overclaim check), 5 claims registered (14 total); 17 idea-level C16 entries (180 total, discrimination-tested) enumerated 85 live loci, all fixed; ~25 further SMALLs fixed; 3 SMALL to debt. Guard pass: E* guard rebuilt to model the real -11/8 error (half the total) + paper-text check; P1 plateau (n=5) and CG non-monotone tests; C25 fire tests for the six escaped strings, the fake K owner and nested braces. Pre-checkpoint self-check (fresh Opus claims reviewer on the full diff + summary surfaces): 5 LARGE + 9 SMALL + 6 NIT, all siblings the sweep had not reached, all fixed -- Paper 55's own abstract/conclusion/Sec. 5 still said "classification" / "every physical observable decomposes" after the citers were scoped (now census + conditional); the field guide contradicted itself ten lines apart; the group1 synthesis paragraph after the edited sentence still said the agreement "identifies" the M1 signature; Paper 20's body still called the S^3-triple bound a basis-truncation rate for the Hamiltonians; Paper 2 still claimed a "complete enumeration" of mechanisms. SMALL: P32 deep-wall "composition of these two facts", Door-4d wording vs Upgrade B, thermal paragraph overreach (now "disfavours", with the gauge/scalar thermal masses named); Gaudillot-Estrada credit added without \cite in P44 (bibitem added) and P40; P42 recap still printing C3*gamma; P57 abstract verb and a fact-stated M1 reading; gh_convergence "rate proven unconditionally" (test now pins the honest tier); Hekkelman-McDonald metadata in P40/P42/P44/P45/group1 synthesis; README (P22 1.44% -> 6.06%, nuclear 592/712 -> 688/828, P12 99.09% -> 99.97%) and the onepager's whole retired pair-diagonal paragraph (334 / Q^2.5 / 11.11 x Q, rewritten from registered exact-rule figures only). Four C25 false-fail risks it found are closed and pinned by selftest probes: owner windows were not LaTeX-normalised, the wh1 overclaim fired on the corpus's correct form, the K pattern matched the Mellin index k and "theorem-grade", and a negation did not cover a list. All gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.17 @ b429344, FULL #11, certifying, unseeded; Opus panel + completeness critic) — **FAIL; RESEARCH-GRADE NOT earned.** 5 LARGE, no keystone math defect for the third FULL run running: P7 He 'do not bracket' (they do: 0.004% is an extrapolation below exact); P7 Parseval single-integral credited to S^3 at 3 summary loci; P32 Door-4 'mostly forced' / 'relative theorem of necessity' / 'disjointness forbids selection' / 'the one inner datum the construction does reach'; P1 [PANEL-VERIFIED] 'higher l has more connections' (false: degree peaks at l=1); claims_register row 17 case-exhaustion as INTERNAL THEOREM (owner CONDITIONAL). ~65 SMALL. Remediated (PI-approved) gate-first: 4 idea-level C16 entries (162 total after this step; 163 with the self-check entry; discrimination-tested, 0 live); row 17 + row 8 + INDEX L11 fixed under the new C25 claim registry, which found them on its first run. ~60 SMALL fixed on sight (P32 x17 incl. CODATA edition behind 0.25% (2022 gives 0.52% the other side), prop=infinity vs CvS Def 2.39, U(1) not the N_c=1 limit, G3 closed negative, decouple-exactly scoped, Forced-Count 'three routes'; P38 x10 incl. 4/pi twice the quotient, refuted Berezin-pair refinement, (d.iii) MEASURED, L4 tier; P7 x12; P0 x7; P1 x4; synthesis x5; code x6 incl. an off-block refusal in lattice_spectrum.lambda_max_from_operator); Malkin-Man'ko + Fronsdal SO(4,2) priority, de Broglie, Brouwer-Haemers, Kac, Hekkelman-McDonald (J. Funct. Anal. 289(12) 111154) cited. One SMALL to debt (P32's 19 bibitem-less inline attributions; blocked on primary-source verification). Guard pass (separate, fire-tested): tests/test_trunk_qa_full11_guards.py (E* = -729/256; degree not monotone in l; off-block refusal). Instrument change (v5.18.0, MINOR): OWNER-LARGE / SURFACE-LARGE split in the RESEARCH-GRADE rule; DELTA reshaped to a claim-search DELTA; new deterministic gate C25 (claim registry). Pre-checkpoint self-check (fresh Opus claims reviewer on the diff + summary surfaces): 4 LARGE + 5 SMALL + 7 NIT, all in the PM's own remediation, all fixed -- (1) P0 conclusion made the linear spacing a consequence of the 2k-1 capacity it is an input to (circular; now input + converse); (2) the (d.iii) re-tier left P38's main-theorem statement and 4 restatements (P32, synthesis x3) at [PANEL-VERIFIED] for all n_max -- now [MEASURED] finite-range, new C16 entry p38-uniform-bound-panel-verified; (3) the P32 Door-4 fix contradicted the unchanged seam paragraph ('free because it generates its own ring', 'the seam already forbids selection') -- reworded to an Observation-level reading, C16 pattern widened; (4) claims_register row 16 still gave K's match as 8.8e-8 (that is the cubic's root; K is 4.8e-7). SMALL: the P38 M1 'tagging convention' insertion was transplanted from the B = 42 fix and misattributed to Paper 55 (reverted to its M1 witness table); P7 'own results' framing; Hekkelman-McDonald metadata unified across P32/P38/P39; criteria/qa wording (claim-search DELTA; 'every gating dimension' defined for a DELTA's dispatched scope; stale 'most runs found zero LARGE' dated); the P1 guard now asserts the printed minima and the top-shell truncation. The reviewer confirmed the He bracket, the Parseval re-crediting, the degree values, L4's [INTERNAL THEOREM] tier, prop = infinity vs CvS Def. 2.39 and the K tripwire clean. All 15 gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.16 @ 834c78f, DELTA #23, claim-impact, unseeded; claims + claim-impact(+citations) + code, Opus) — **DEFECTS; RESEARCH-GRADE NOT earned.** 1 LARGE (10 -> 1 -> 1) + 1 borderline, BOTH the PM's own prior-round residue rather than corpus debt: (1) P32 L2431 (keystone proof sketch) and P35 L846 'engaging a mechanism does not force pi' stated the converse FALSE against the 'not established' ruling (siblings: P55 theorem box, test docstring 'B = 42 (proven)'); (2) the DELTA #22 relabel called the SU(2)xU(1) value 2/(n+1)*1 a 'joint mass maximum' at ~6 code/test loci, when Paper 45 states it is a HYBRID (SU(2) mass max x U(1) symbol max), neither a mass maximum nor a cb-norm. Math clean (465 tests pass); the joint-UCP guard is sound but its positivity leg could be evaded by an aliased negative lobe. ~18 SMALL (residual cb-norm/Plancherel-symbol conflation in berezin_reconstruction, central_fejer_su2, Cesaro-2 test, gh_convergence; Bozejko-Fendler credited in code + P40 + orphan bibitems; P57 L347 'witness derivation'; synthesis 4/pi 'quotient' -> 'twice the quotient'; claim-matrix 'structural corollary'). Remediated gate-first: 3 new idea-level C16 entries over papers + tests + geovac (pi-source-converse-stated-false, lorentzian-hybrid-labelled-mass-max, bozejko-fendler-as-support; 158 total, all discrimination-tested, 0 live). Guard repaired: positivity now ||K||_{L1} = 1 by quadrature (equivalent to K >= 0 for unit mass), fire-tested on the exact aliasing kernel that escaped. Two 119-382 s tests marked @slow. Pre-checkpoint self-check: 0 LARGE, 4 SMALL + 3 NIT, all fixed (a regex in the PM's own relabel had turned 'joint mass maximum' into 'joint hybrid value' inside the new hybrid docstring; lorentzian field docstrings; joint_plancherel_symbol hybrid note incl. a U(1) support-convention discrepancy vs P45 Lemma L2(b), flagged; residual Bozejko-Fendler in central_fejer_su2; synthesis orphan bibitem). All 14 gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.15 @ 635c9f8, DELTA #22, claim-impact, unseeded; claims + claim-impact + code + citations, Opus) — **DEFECTS; RESEARCH-GRADE NOT earned.** 1 LARGE (10 -> 1): the cb-norm correction had not reached the CODE layer -- tests and geovac modules behind Papers 39/45 still labelled the Plancherel mass maximum 2/(n+1) / 4/((n_a+1)(n_b+1)) a cb-norm, and P39's joint cb-norm = 1 was untested. Claims 0 LARGE (every corrected result re-derived; K tripwire clean); claim-impact 0 LARGE (all nine changed claims clean corpus-wide by idea); citations 0 LARGE. ~25 SMALL. PI rulings: the gauge group stays in Paper 57's forced bucket 'relative to the adopted Upgrade B axiom'; the [SYMBOLIC PROOF] tag on the 18 checks stays. Remediated gate-first: a new code-layer C16 entry over tests/ + geovac/ (155 total) caught 4 further mislabelled loci the reviewer had not listed; docstrings/comments/test prose relabelled across ~13 files (identifiers kept, NAMING NOTEs added). Paper SMALLs: Bozejko-Fendler withdrawn (discrete groups) at the group1 synthesis, P32, P39 and the code; Paulsen cited (P39); P18 F = half zeta_{S^1}(1); P55 abstract 'uniquely identifies' -> 'indexes'; P35 section title; case-exhaustion scope residue; P45 header+caption; P50/P56/P57/field guide/INDEX wording. Guard pass: joint UCP guard (Paper 39 L2-T(c)). Pre-checkpoint self-check: 1 LARGE -- my own new wording called M1's B = 42 a 'proven witness' so that the converse 'does not hold in general'; B's M1 membership is Paper 55's tagging convention, not a proof, so the converse is restored to 'withdrawn: not established' -- plus 7 SMALL (formula lines still written as ||.||_cb, residual Bozejko-Fendler in code, P39 L1405, P50, P57, P45 caption) and the first joint guard was near-tautological (a normalisation it could not see); it was rebuilt to integrate the actual kernel against Haar measure and fire-tested on three wrong answers. All 14 gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.14 @ 1cae4e8, DELTA #21, claim-impact, unseeded; claims + claim-impact + code + citations, Opus) — **DEFECTS; RESEARCH-GRADE NOT earned.** ~10 LARGE loci in 6 families, all PM-verified; no keystone math defect; citations 0 LARGE; every new guard fired; the re-pointed WH1 falsifier genuinely covers the top band. The claim-impact reviewer searched by IDEA across the whole corpus and reached pre-existing restatements no prior round had touched, so the count rose rather than fell. Families: (1) pi-source case-exhaustion restated as a proven/complete EQUIVALENCE (group1 synthesis thm, P35, P55 thm, field guide 'classification is complete'); (2) the FULL #10 re-scope was itself wrong -- M2 at half-integer s (zeta_{D^2}(5/2) = 14 zeta(3) - 31/2 zeta(5)) is a second pi-free witness -- so, PI-approved, the converse (iii)=>(i) is WITHDRAWN outright (theorem now (i)<=>(ii) empirical + (ii)=>(iii) over the first fifteen projections, [CONDITIONAL]; witnesses M1 rational core B = 42, M2 half-integer s, M3 at convergent s = 4: 8(beta(4) - G)); (3) retired cb-norm 2/(n+1) and 'Plancherel symbol' in the group1 synthesis, and P39 Lemma L2-T (b)/(c) + remark + two tables stating 4/[(n_a+1)(n_b+1)] as a cb-norm (UCP => cb-norm 1; these are mass-distribution maxima; no P39 bound used them); (4) F as the Fock-Dirichlet evaluation 'on the scalar Laplacian' (P55); (5) 'discreteness is compactness' biconditional (group6 synthesis abstract); (6) '18 independent proofs' (P11). ~35 SMALL fixed on sight (M2 integer-s scope at 3 table cells + P24/P18; P57/field guide 'theorem' labels; P57 gauge group out of 'autonomously generates'; 'disjoint rings' P34/P50; Paper 2 F = half zeta_{S^1}(1), value unchanged; Cesaro-1 label circle-only; EMN/Deligne citations; 2017 -> 2018; a self-inflicted matrix error, Forced-Count test collects 7 not 47). Remediated gate-first: 5 new idea-level C16 entries (154 total), all discrimination-tested. Guard pass (fire-tested): two converse-withdrawal witness tests; Paper 0 'linear spacing required' upgraded [OBSERVATION] -> [SYMBOLIC PROOF] (telescoping test). Pre-checkpoint self-check: 0 LARGE, 10 SMALL sibling loci, all fixed (notably the M2 witness vs 'M2 is pure-Tate' at P32/P56, now scoped to SD coefficients and integer s; 'no known pi' in place of 'pi-free'). Process note: two heredoc backslash-halving slips during remediation (both caught; one left a missing comma fixed directly). All 14 gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.13 working tree, FULL #10, certifying, unseeded; 10 Opus reviewers + completeness critic) — **FAIL; RESEARCH-GRADE NOT earned.** 10 LARGE, all PM-verified against the files, none a mathematical defect in a keystone (every keystone re-derived independently: P38 bound + 4/pi to 6 digits by a fit on computed sums; 2pi^2/n to n = 12 by Funk-Hecke; C_sat; block spectrum; forced count; B = 42; the Casimir ladder). Code + citations: 0 LARGE for P0/P1/P7. LARGE: (1) P38 intro torus-rate / 'Cesaro-2' residue of a withdrawn comparison; (2) the WH1 falsifier test_p38_action_seminorm runs on the DOUBLED model (cannot see the top band) yet P38 L645, matrix rows 70/71/150 and CLAUDE.md named it the injectivity backing (true backing: test_p38_true_ch_window); (3) P32 Berezin map in the refuted N/Z, N <= n_max form; (4) P32 retracted 1:2:3 / 'trend toward correlation' still listed as robust (+ 86.6020 != 50 sqrt3); (5) P32 'categorically disjoint' (withdrawn strong form, false: both rings contain Q); (6) P32 Catalan irrationality asserted (open); (7) P32 F = zeta_{Delta}(2) FALSE (= 0.885; F = zeta_{1-Delta}(2)); (8) P32 case-exhaustion (iii)=>(i) fails for M3 (G - 1/2, no pi); (9) P32 'seven theorem-grade non-selection results' (all proof-by-argument, one with an empirical leg; Forced-Count is a dimension count); (10) group3 synthesis T9 ring unscoped (half-integer s carries odd zeta). About 60 SMALL fixed on sight (incl. P7 master formula credited as Parseval [PRIOR ART], P0 block spectrum credited as the textbook path-graph spectrum, Pearson 0.81 -> 0.846/0.796 at P31/P32/P57, P1 'curvature' -> discrete gradient, WH3 register language out of the synthesis, stale test docstrings). Remediated gate-first: 11 new idea-level C16 entries (149 total, all discrimination-tested), mirror sweeps into P18/P31/P40/P55/P57/field guide. Guard pass (separate, all fire-tested): test_paper32_f_zeta_and_m3_scope (3), order-one counts 43/196 -> 1337/3025 pinned, test_paper7_label_set, linear-spacing counterexample in test_trunk_qa_annular. Pre-checkpoint self-check: 5 LARGE sibling loci (P32 theorem statement still 'equivalent'; P57 unscoped classification + gauge group 'autonomously generated' + retired discreteness<->compactness biconditional; P40 2Vol(base)/Vol(group) contradicting P38; P18 disjoint-ring sibling) + 7 SMALL, all fixed; it also exposed that the claims reviewer's 'zeta(odd) is not M3' SMALL contradicted the owner Paper 55 (zeta(3) in M3 cap MT(Z)), so that edit was reverted. Run note: one reviewer fire-tested a shared test file while another ran it (one transient failure; tree clean) -- fire tests go on copies from now on. P32's slow Lorentzian panel timed out with two uncaptured failures: INCOMPLETE, not a defect. Owed: P40's per-summand L3 route is un-re-audited against P38's 2026-09-28 withdrawal (C3 <= 1 itself now rests on compression). All 14 gates PASS; C10 whole corpus.
- 2026-09-30 (v5.17.12 @ 9c7561d, DELTA #20, claim-impact, unseeded) — **DEFECTS; RESEARCH-GRADE NOT earned.** 4 LARGE (12 -> 8 -> 5 -> 4 -> 4) in three families: (1) Hopf-base M1 label in wordings the prior gate missed (P32 obs:pi_coulomb 'Hopf-base Weyl constant'; group1 synthesis 'the Hopf bundle produces Vol(S^2)=4pi', false: the base S^2(1/2) has area pi); (2) Paper 28 T9 -- the OWNER theorem -- false at half-integer s ('no odd zeta at any s'; 'sqrt(pi) populated at half-integer s'): zeta_{D^2}(5/2) = 14 zeta(3) - (31/2) zeta(5), PM-verified to 30 digits; copied into P34/P55 by the prior round; (3) 'from M1 alone' -- GeoVac presented as deriving the Hawking/Unruh beta (P32, P34 x6, P35 x3) against P34's own honest scope. Code CLEAN (grading guard sound; sign blind spot fixed with a shifted-D check); citations CLEAN (BW/Sewell/Unruh/HH and the Hopf S^3(r)->S^2(r/2) geometry confirmed). Two of the three families PRE-DATE the remediation (the owner-theorem error; the M1-alone phrasing): the DELTA sweeps now excavate adjacent pre-existing debt, so PI directed the next round be a trunk FULL run. Remediated gate-first (2 new idea-level C16 entries + hopf-label widened; 28 live loci -> 0); T9 rewritten as an integer-s theorem (residues -1/4, 1 and heat-trace coefficients PM-verified by hand); new guard test_paper28_t9_half_integer.py. Pre-checkpoint self-check: 2 LARGE (M1-alone stragglers P32 L3722, P34 L10281; gate widened) + 6 SMALL fixed (incl. half-integer odd-zeta scoped to s >= 5/2 -- negative half-integers are rational; P33 Hopf-base area; P34 'every pi' softened). All 14 gates PASS; C10 whole corpus.
- 2026-09-29 (v5.17.11 @ 96778b5, DELTA #19, claim-impact, unseeded) — **DEFECTS; RESEARCH-GRADE NOT earned.** 4 LARGE (12 -> 8 -> 5 -> 4) in two families, both twins of the prior round's fixes at unmoved loci: 'the framework's M1 2pi IS the modular-flow / boost-orbit period' (P32, P34 x4, P35) and the retired sqrt(pi) Q (+) pi^2 Q ring (P32, P55). Code CLEAN (control genuine at n_max=2; grading leg never exercised -> new synthetic test); citations CLEAN (GK per-component F_f, Rittenhouse->Cavagnero, BBT rename confirmed). PI-approved APPROACH CHANGE: idea-level C16 entries (`framework-2pi-is-boost-period`, `sqrtpi-ring-anywhere`) keyed on the claim's subject, not its phrasing, with EVERY corpus hit swept (continuum Sewell/Unruh statements deliberately silent); and the deferred 'Hopf-base measure' label swept corpus-wide (~77 loci, 16 papers) with a gate (`hopf-base-measure-label`) that then caught two line-wrapped survivors. Pre-checkpoint self-check: 0 LARGE, no overcorrection; 6 SMALL fixed incl. a PM-seeded error ('a quarter of the area of the Hopf base' is false under the corpus's own S^3(r)->S^2(r/2) convention; pi is the WHOLE area of S^2(1/2)) and P55's a_k(D^2) in pi*Q (normalised coefficients are 4pi^2, -2pi^2). All 14 gates PASS; C10 whole corpus.
- 2026-09-29 (v5.17.10 @ c81ae22, DELTA #18, claim-impact, unseeded) — **DEFECTS; RESEARCH-GRADE NOT earned.** 5 LARGE (12 -> 8 -> 5), all one-sentence residue at unmoved loci of changed claims: P32 L7486 ('H_local IS the BW vacuum Hamiltonian'); group3 synthesis L320 + P1 L414 ('not an output of these operators' next to the Casimir); P18 thm:k_decomposition defining K inside a Theorem; P34 'IS already Lorentz-boost physics'. Code CLEAN (0 LARGE; Richardson combination verified exact to 1/n^2; pinned P M_g P values match the paper verbatim). Citations CLEAN (0 LARGE; Rieffel Thm 2.3, BBB odd-m extension, GK, Rittenhouse all confirmed from primary texts). Remediated gate-first (5 new C16 entries, 9 live loci -> 0; SMALL sweep: C_sat rename, Beccaria key, proof parenthetical sharpened, P53 Stein-Weiss credit, two new fire-tested guards). Pre-checkpoint self-check: 1 LARGE (P34 'already captures the Lorentz-boost content', gate widened) + 5 SMALL, all fixed; incl. a guard that pinned SVD noise at n_max=1 (grading > commutant, impossible) -- now asserts only the robust collapse there. Classical Bochner-Riesz facts re-tiered [PRIOR ART]. All 14 gates PASS; C10 whole corpus.
- 2026-09-29 (v5.17.9 @ 47c6bcb, DELTA #17, claim-impact, unseeded) — **DEFECTS; RESEARCH-GRADE NOT earned.** 8 LARGE (16 -> 13 -> 12 -> 8), all retired forms at loci whose bytes did not move: P0 abstract (E_n 'not a function of the Casimir'); group3 synthesis abstract ('framework's reproduction' of CFT values); P32 'alpha observation at 8.8e-8' (x3); P32 BW section 'boost-class generator' / 'analog of the BW vacuum' (x3); INDEX (P50 M2/M3, P2 8.8e-8); P2 'alpha is a root of a polynomial over Q(pi)'; P16 novelty (body + group4 synthesis); P2 K match inside thm:three_homes (touches sec. 13.5). Code CLEAN (0 LARGE; 138 tests; every aimed fire-test FIRED; P M_g P negativity verified at the operator level). Citations CLEAN (0 LARGE; GK eq. 2.1 covers S^7; GE-vS Prop. 14 is the lifted-state bound). P53 first review: proof and every constant re-derived by two reviewers, 0 LARGE. Remediated gate-first (6 new C16 entries + 1 widened, 26 live loci -> 0, 12 dependents stamped) by three disjoint sweeps; PM caught a false counterexample in sweep B's new P2 text. Pre-checkpoint self-check: 3 LARGE (P18 twin of the counterexample + pi^2 for pi^3; the cubic's-root transcendence is PROVABLE, one-line Lindemann argument, was 'likely (not proven)'; group1 synthesis 'bit-exact' body) + 6 SMALL, all fixed; new C16 entry `cubic-root-transcendence-understated` and guard `tests/test_paper2_cubic_root_transcendence.py` (fire-tested). All 14 gates PASS; C10 whole corpus.
- 2026-09-29 (v5.17.7 @ ba34fc8, DELTA #16, claim-impact, unseeded) — **DEFECTS;
  RESEARCH-GRADE NOT earned.** ~12 LARGE (16 -> 13 -> 12). Code CLEAN (0 LARGE; every
  new guard fires; the n_max=4 J O J^-1 guard rejects a block-phased J that passes the
  sign tests); citations CLEAN (0 LARGE; Rittenhouse/Cavagnero/FdlR/BBB/BDF/CHM all
  confirmed from primary sources). LARGE by family: (1) Fock <-> Delta_{S^3} framing
  (P7 roadmap, P0 x2, synthesis x5) -- PM finding: the v5.17.6/7 wording "not a
  free-particle problem / Delta is not the link" OVERCORRECTED DELTA #15's genuine
  Green's-function defect (H_bound = -(Z^2/2)(1-Delta_{S^3})^{-1}; eigenstates are
  the S^3 harmonics; the kernel ~ (1-Delta)^{-1/2} fixes p0) -- PI decision; (2) P32
  rate tier (thm:gh_convergence "proven unconditionally") and rem:wh1_proven
  (structural AC claim "proved"; only the commutative factor converges) -- WH1
  wording, PI decision; (3) P18 K "equals alpha^-1 to 8.8e-8" (true gap 4.8e-7; 8.8e-8
  is the cubic root) and "alpha is a root of the cubic" (pre-existing); (4) P50
  one-sentence summary "realizes ... modular structure"; (5) group1 synthesis HS
  closed form stated at N_t = 11 (owner: N_t = 1 only); P34 live falsified M3
  prediction; P34/P31 BBB (3,0)/(3,1) residue; P42 abstract untagged signature-blind
  clause. INSTRUMENT: this file's C8 extension ("LARGE: asserting the finite
  truncation is odd") now contradicts the owner (true window odd) -- rescope owed
  after the run (criteria frozen during it).
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
