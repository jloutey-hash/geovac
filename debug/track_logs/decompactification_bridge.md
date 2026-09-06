# Track log — decompactification bridge (PM holds this thread)

**Opened:** 2026-09-05. **PI direction:** dispatch agents, PM maintains context.

## The conceptual frame (do not lose this)
Two open items from the propinquity retraction:
- **NC Mondino-Samann extension** — Lorentzian (time-separation function). NOT the bridge target; scoping correction made to PI.
- **Metric-level de-compactification** (Paper 47 T->infinity limit) — Riemannian. THIS is the bridge target.

PI's intuition: decompactification feels Boolean, but sliding two atoms together at constant R-rate should decompactify *predictably*. Resolution offered: the LABEL is Boolean (m exact / l dies), the COUPLING is continuous. Precedent already in corpus: composition arc found irreducibility is a knife-edge (Boolean, about center count) while the nested-commutator magnitude decays smoothly (continuous, the bonding signal). "Two invariants, not one."

Wavefunction answer given to PI: prolate-spheroidal eta-equation deforms l continuously (Paper 0 "Level 2"); phi/m never decompactifies; Paper 8 bond-sphere angle gamma(R): pi (united atom) -> 0 (separated). Chemistry name = Mulliken correlation diagram, governed by the non-crossing rule (same m cannot cross) = Paper 58 abelian residue drawn as a picture.

Per-shell front (READING, untested, from an Aug memory note only): two-center character peaks at r_n ~ R i.e. n ~ sqrt(ZR); as R shrinks the window sweeps toward the core; core decompactifies last.

Bridge shape: Paper 47 has inner compact rate + outer de-compactification rate, coupled by an admissible scaling ruled "convention". Two-center family has the same two dials (atomic cutoff, R). Prediction (untested): a finite-truncation Lipschitz seminorm sees only the compact part (phi circle + atomic spheres) and treats R like it treated signature -> "decompactification is convention at every computable level" showing up a third time. Every released axis costs one transcendental seed; R enters through the seed argument = Paper 18 compactness thesis as an integral.

## Dispatched 2026-09-05
1. **Explore (read-only):** does any synthesis / Paper 58 already state m <-> correlation-diagram / non-crossing rule? PI: "I hope that's in a synthesis somewhere."
2. **general-purpose (measurement):** the R-sweep diagnostic. Deliverable 1 = per-shell front (window vs R, fit vs sqrt(ZR)). Deliverable 2 = four aligned curves vs R: gamma(R) [P8], ||[P_A,P_B]|| + principal angles [comp arc, anchor to LiH 0.50 / 7.6-44.7-67.3 deg], aggregate l-mixing, transcendental seed a(R) [neumann_vee.py]. Gate: GO continuous / BORDERLINE scaling-off / STOP discontinuity.

## Status: agents running. PM awaiting reports, then synthesize + decide paper capture.

## Explore result (2026-09-05): m <-> correlation-diagram = GENUINE GAP
Connection stated NOWHERE. "correlation diagram"/"non-crossing"/"Walsh"/"Mulliken correlation" = 0 hits in papers/.
- Paper 58 states the abelian-residue theorem (Thm 1, L594-613; proof via axial U(1)/Pontryagin L627-652) but never connects to the correlation diagram / non-crossing rule.
- group2 synthesis restates it (L588-619), same stop-short.
- Paper 0 Level 2 (L489-501) closest: "united-atom limit", "m is exact" — but no diagram, framed as H2+ separation.
- Field guide: doesn't cover P58 abelian residue at all.
BEST HOME: Paper 58 § "Mechanism: the abelian residue", remark after proof (L652) before the shibuya_wulfman corroboration (L654). Secondary: group2 synthesis after L617-619 italic summary.
CAVEAT to honor when writing: non-crossing rule is same-m AND same g/u (homonuclear) AND same spin; state precisely, don't overclaim.
STATUS: capture DEFERRED — fold into one decision once R-sweep lands.

## R-sweep result (2026-09-05): BORDERLINE — continuous, but scaling corrected
Continuous crossover confirmed (no discontinuity in any curve) => "Boolean label / continuous coupling" reading HOLDS. Not Boolean.
CORRECTION to PM's hypothesis: front scales n*(R) ∝ R^1.0 (fit 0.98/1.02 two routes), NOT √(ZR)=R^0.5.
Sharper truth: front set by orbital EXPONENTIAL DECAY LENGTH ℓ_n=n/Z (tail reach), not mean radius r_n=n²/Z.
  - R*/ℓ_n = [1.56,2.19,2.19,2.13,2.12,2.14,2.15,2.14] FLAT (~2.14, <3% drift n≥2)
  - R*/r_n  = [1.56,1.10,0.73,0.53,0.42,0.36,0.31,0.27] DRIFTS
  Physics: two-center overlap is a TAIL phenomenon -> governed by decay length, not where bulk density sits. Correct + cleaner than √(ZR).
Four curves all continuous: gamma(R) monotone 157→9.5°; principal angles sweep smoothly 0.8/2.1/7.4° → 44.6/88.2/88.6° (||[P,P]|| itself capped at ½, carries no scaling); Σchar_n single smooth peak; seed a(R)=R linear, e^a E_1(a) smooth; bonus η-eq l-mixing rises 0→~0.5 (metric bound).
VALIDATION: LiH anchor reproduced EXACTLY (7.6/44.7/67.3°); Mulliken closed-form vs grid quadrature <1.3e-8; η-eigenproblem matches prolate solver bit-for-bit.
Files: debug/decompactification_R_sweep.py, debug/data/decompactification_R_sweep.json, debug/sprint_decompactification_R_sweep_memo.md (1293w).
Caveats: single-electron H2+ by design; e-e correlation untested; m-exactness asserted (structural φ δ_mm'), not re-measured.

## CAPTURE DECISION (pending PI): two candidates, recommend but do NOT auto-edit (user in explore/think-aloud mode, no approval given).
(A) Paper 58 remark: m = correlation-diagram organizing label / non-crossing rule (gap confirmed). Safe conceptual bridge. Honor g/u + spin precision.
(B) NEW finding: decompactification is a continuous R-crossover set by decay length n/Z (linear front), not √(ZR). Single-electron diagnostic. Recommend PI review before baking into a paper claim.

## 2026-09-06: capture (A) DONE — Paper 58 correlation-diagram remark
Placed `\paragraph{The correlation-diagram reading.}` [OBSERVATION] between Thm 1 proof end and the shibuya_wulfman corroboration. +2 bibitems (neumann_wigner1929 Phys.Z. 30, 467; herzberg1950 vol I 2nd ed Ch. VI). Matrix row added (rests on: Paper 11 eta-eq united-atom limit). Compile 0 err / 0 undefined / 11pp. Gates: C19 PASS, C16 group2 PASS, C17 group2 PASS, C13 PASS. NOT committed (PI did not ask).
Capture (B) — decay-length front — still HELD for PI review.

## Machinery survey for PI Q: multi-electron + unequal Z (NOT mass — mass doesn't enter electronic problem; BO)
Unequal Z (heteronuclear):
- 1e exact: geovac/prolate_spheroidal_lattice.py takes Z_A, Z_B (heteronuclear term b = R(Z_B - Z_A)) -> eta-mixing front directly available for HeH2+ etc.
- overlaps: Mulliken A_n/B_n + two_center_grid_lm.py take per-orbital exponent -> decay-length front with l_A = n/Z_A != l_B = n'/Z_B. Sweep already anchored on LiH (unequal Z) for principal angles.
- Papers 8-9 Cor dual_p0: NO shared p0 for heteronuclear -> Paper 8's gamma has no single heteronuclear form -> heteronuclear = TWO focal lengths (cf. Paper 39 tensor-product "two focal lengths"; multi_focal_wall fires at 2). Lead, not result.
Multi-electron:
- THEOREM already in corpus (Paper 60 sec "the metric does not compound"): k-electron config overlap = k-th compound matrix of the 1e overlap -> overlap/metric-level front is ONE-ELECTRON by identity; N-electron adds nothing at that level.
- What e-e CAN change: effective decay length via screening (Z_eff per block, composed_qubit.py) -> l_n = n/Z_eff. Testable: does the front shift by Z->Z_eff and nothing else?
- Solvers with R-sweep: Level 4 (H2; HeH+ tests in test_level4_multichannel.py; united-atom limit verified P15), NOCI NaH ladder (noci_engine.py rows by R), balanced_coupled, composed. Paper 60 2e H2 Sturmian CI (dissociates correctly).
Mass: only hyperfine_a_constant / cross_register_vne (nuclear-electronic PoC, Track NI) / magnetization. Isotope/mass = nuclear-motion axis, orthogonal to electronic decompactification.

## 2026-09-06: PI "let's start it" -> DISPATCHED correlation-vs-front diagnostic
Question: does the decompactification front move by Z -> Z_eff and NOTHING else?
Ladder: (0) 1e heteronuclear law (HeH2+ eta-mixing + unequal-decay-length overlap front: is R* = c(l_A + l_B)?) -> (1) HF on H2/HeH+ (screening only) -> (2) FCI (screening + correlation). Residual beyond Z_eff = correlation touching compactness.
Hard control: basis must have radial flexibility so l_eff is an OUTPUT (basis-doubling control), else tautological.
Gate: GO = 1e law at empirical Z_eff matches within ~10% both molecules; BORDERLINE = one molecule only; FINDING = systematic residual.

## 2026-09-06 ~01:05: correlation ladder REPORTED — GO (principal angle) + FINDING (coherence)
GO: occupied-one-body-space front (M1, occupation-weighted principal angle, = Paper 60 compound-matrix object) matches 1e 1s-1s law at EMPIRICAL Z_eff(R): H2 HF -0.7% / FCI -1.2%; HeH+ HF 0.0% / FCI -0.2%. Both rungs, both molecules, 3 bases, 2 ERI grids. => screening moves the front, nothing else, to ~1%.
FINDING: signed cross-center coherence (M2, bond-order-weighted overlap) front moves INWARD at FCI: H2 -23% (0.987 vs 1.285), HeH+ -9% (0.690 vs 0.758). Mechanism: antibonding NO enters per-center density with weight (1+S)/(1-S) ~5 at S=0.67, so even n_u=0.025 pulls coherence from S (MO limit) toward S^2 (Heitler-London). OCCUPATION-NUMBER effect (bond-order collapse), NOT a length scale. Two-orbital formula 0.591 vs measured 0.596 — PM re-did arithmetic: checks.
Rung 0 surprises: (1) absolute |S|=1/sqrt2 front EXISTS only for exponent ratio t < 2.746 (united-atom overlap must exceed it); (2) tail-reach front of two unequal tails = GEOMETRIC MEAN 0.79*2sqrt(l_A l_B) (max dev 2%), NOT additive (23%) nor max-law (62%). HeH2+ eta-mixing turns on FIRST order (b*eta, ~R^2.1) vs homonuclear second order (~R^4.0); heteronuclear deficit does NOT saturate near 0.5 — label lost to LOCALIZATION on He focus.
Hard control PASSED: zeta_eff is an output (1.46->0.83 across R; HF vs FCI differ 12% at R=4); basis doubling dR* <=0.8%. Route-A/B zeta agree 0.3% at front. Known route limitation: single-center multipole (aa|aa) self-repulsion of steepest B function at R>=4 is a 1-3% defect not cured by grid/L_max; harmless (2e-5 on all fronts) but RECORDED not certified away.
Caveats: s-only basis (no p-sigma control); 2e only (N>2 untested); M3 subspace measure is occupation-blind, its FCI crossing = threshold artifact.
STATUS: agent's final driver re-run (PID 22660, started 01:05:12) still running when reported; JSON on disk = 01:03 (pre-final-edit). PM waiting for exit, then re-verify numbers vs memo before any capture decision.
VERIFIED 01:0x by PM from on-disk JSON (01:03:19): H2 HF M1 1.2828/pred 1.2918 (-0.70%); H2 FCI M1 1.2696/pred 1.2845 (-1.16%), M2 0.9873, M3 3.1508; HeH+ HF M1 0.7571/pred 0.7569 (+0.03%); HeH+ FCI M1 0.7567/pred 0.7583 (-0.21%), M2 0.690. Memo table == data. Two-orbital M2 formula re-derived by hand: 0.591 OK.
Driver PID 22660 (started 01:05:12) still alive at check; JSON not yet overwritten. Background waiter bmyiqe680 will re-invoke PM on exit -> re-diff JSON vs memo then. Verdict reported to PI on the verified 01:03 data.
CLOSED 01:08: driver re-run exited 0; JSON regenerated 01:08:05 (172 s). All four M1/M2/pred fronts identical to the verified 01:03 values (worst |delta| = 0.0000); only addition = 'lmax_control' block (24->40, <=2e-5 on every reported quantity, as memo states). Verdict stands on final data. No agents running. Three measured structural results now in memos awaiting PI placement: (1) decay-length front R* ~ 2.14 n/Z (v. 09-05 memo); (2) correlation moves only the coherence front, via antibonding occupation (1+S)/(1-S), not any length; (3) unequal tails combine as geometric mean 0.79*2sqrt(l_A l_B). Paper 58 correlation-diagram remark WRITTEN, gates pass, NOT committed.

## 2026-09-06 CAPTURE IN PROGRESS (PI: "let's put em in papers")
- Paper 58: new Sec. continuous side (front, exponent, geometric-mean tail law, t_c) + backing-table rows; gvq macro. Paper 60: sec:manyelectron two new paragraphs (front moves only via Z->zeta_eff; coherence front via occupation, eq:coherence_two_orbital). Paper 11: eta-label drift paragraph.
- Tests: test_paper58_decompactification_front.py (5 fast), test_paper60_coherence_front.py (2 fast + 1 slow ladder), test_paper11_eta_mixing_onset.py (2 fast). Registry +11 keys. Matrix +5 rows. CHANGELOG v5.10.2. CLAUDE.md v5.10.2 + one-liner. Memory axis-map corrected.
- INDEPENDENT-ROUTE CATCH: agent t_c=2.7456 was a scan-grid artifact; closed form + quadrature give 2.664 (S(R) monotone in R for all t). Fixed in paper/registry/CHANGELOG/matrix; memo carries a correction note.
- Pre-existing, NOT mine: Paper 11 `fig:pes` undefined at L759 (noted).

## 2026-09-06 CAPTURE COMPLETE (pending PI checkpoint)
- Papers 58/60/11 compile 3-pass: 0 errors, 0 undefined (P11 `fig:pes` L759 is a PRE-EXISTING missing figure, not mine -- reported).
- Fast tests 9/9 green; slow ladder test green (9.5 s); registry self-test 17/17; fire_test: 5/5 guards fire (58a needed a physics plant -- the first quadrature-parameter plant was a no-op, not a sleeping guard).
- Gates: C19 PASS, C16 group2 PASS, C17 group2 PASS, C13 PASS, C22 PASS. C21 was FAIL on ONE pre-existing zombie: Paper 28 L158 "O(Q^2.5) Pauli scaling" (retired pair-diagonal figure, group5, untouched by this arc) -> fixed on discovery to the exact-rule linear statement (27.90 Q). Re-verify C21 next.
- Version cursor: CLAUDE.md v5.10.2 + section-2 one-liner; CHANGELOG [v5.10.2]. NOT committed (PI runs /checkpoint).
- Untracked new files to stage: 3 tests, 2 drivers, 3 JSONs (incl. _quick), 2 memos, this track log.
