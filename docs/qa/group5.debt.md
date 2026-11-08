# group5 QA debt

Created 2026-11-07 by `/qa group5` FULL (incremental, whole target; snapshot 99cd34bf).
Verdict FAIL, remediated (CHANGELOG v6.5.0). This file holds what the remediation did not
close: owed computations (PI-approved "re-tier now, owe sprints"), owed guards for the
guard pass, PI decisions outstanding, and SMALL findings not fixed on sight.

## PI decisions (resolved 2026-11-07, all on PM recommendation)

1. Retitles applied: Paper 33 "... Graph Topology, Angular Momentum, and the Diagonal alpha_z Vertex";
   Paper 36 "Bound-State QED on S^3: The Leading One-Loop Hydrogen Lamb Shift"; Paper 41 "... Compatibility
   with Compact U(1) on a Non-Cubic Substrate, Including Monopole-Density and Polyakov-Loop Diagnostics".
   Bibitems updated (Papers 7, 32, group5 synthesis); old titles blocked by C16 `g5-retired-titles-33-36-41`.
2. The 1/(24 pi t) criteria amendment is PI-confirmed.
3. MINOR bump v6.5.0 confirmed.

## Owed computations (pre-register before running)

- **P28 vector-photon vertex.** Selection rule DONE 2026-11-07 (VV stage A, GO; now Paper 33
  Sec. 2): the ground state couples, Theorem 4 is a model result. Still owed (stage B): the
  normalized Sigma and vertex sums with the physical rule, gauge sector, C for Furry.
- **P28 anomalous-moment normalization.** Derive the absolute normalization (e^2, loop
  measure, harmonic normalization on S3) so F2 is a moment in units of alpha, then the
  flat-space limit R -> infinity.
- **P33 physical-vertex census.** Re-score rules 2 and 7 with the full alpha coupling and
  transverse vector harmonics.
- **P36 order-by-order Lamb budget** from Eides-Grotch-Shelyuto: A50, recoil, finite size,
  two-loop B40/B50, before any residual attribution; converge the native 2S acceleration-form
  Bethe log (N=40 is a crossing).
- **P41 discriminating tests.** A 4D (Z^4) control on the same diagnostics; monopole data above
  beta_c ~ 1.01; <|P_lat|> on larger spatial graphs; fixed-L Wilson-loop ensemble.
- **P41 Polyakov-rate mismatch** (c_sigma^MC vs c_rho/2): open.
- **P51 sector-wise Mellin (iii)** re-scoped to the bare 2D tip (phi(0)) vs the 4D cigar's
  leading phi(1); re-check the empirical F14 "sub-5%" claim against this split.

## Guards (guard pass DONE 2026-11-07, v6.5.1)

Rebuilt (each fire-tested against the retired value; 22/22 plants fire; independent review 0 LARGE / 9 SMALL,
7 fixed): `test_paper51_dS_vacuum.py` (from the CH spectrum, rank-4 trace), `test_paper51_cutoff_mellin.py`
(direct spectral sum -> G_eff, Lcc; Laplace-convention cigar tip incl. f(0) != 1), `test_paper51_replica_sbh.py`
(action-route G == entropy-route G), `test_su2_wilson_gauge.py` (cycle basis at n_max = 3, 4, 5; Haar I2/I1 tied
to `wilson_action`; Schur factorization), `test_paper41_wilson_witnesses.py` (paper estimator passes 3; fixed-band
gap; centre symmetry; <|P_lat|>). New: `test_paper2_dirac_cumulative.py`, `test_paper36_a50_term.py`,
`test_paper51_tt_heat_trace.py`.

Still owed (backing for printed numbers, not pins of retired values):
- P30: long-chain MC values (0.232-0.249, 0.429-0.452, three 20 000-sample chains) have no test; the area-2 loop
  check uses the Schur identity on hand-picked matrices, not the lattice's actual 6-cycle.
- P51: continuum tip = 1/6 at t/R^2 = 0.01, 0.02, 0.05 (exact Bessel zeros); a0 continuum 1.888 at R=10, t=0.1;
  B4 error ratio -> sqrt2. The rank-4 bundle and the total-spinor cone coefficient 1/12 are inputs no test checks.
  Possible naming SMALL: the TT spectrum m^2 - 3 is labelled 'Lichnerowicz' (verify the -nabla^2 vs Delta_L
  convention on unit S^3).
- P41: correlator normalization (1/800); centre-symmetry and P_lat now backed.
- P33: diagonal-current M1 projection (-4/3 for the ground state, fg = 1).
- P28: physical-vertex coupling at q = n1+n2+1 (vertex_amp_polarized nonzero for 0 -> 0 at q = 1).
- P36: native-log Lamb ~1055.6; 2S acceleration-form N=50 residual +2.98% (A50 now backed).

Code defects found (not fixed in this run; production code):
- `tests/wilson_rule_b_support/xcwg_polyakov_loop.py:617-618` divides the correlator by n_sample twice.
- `geovac/qed_self_energy.py:72` sets mpmath.mp.dps = 50 at import (breaks
  test_paper55_m3_cyclotomic_mixed_tate under xdist).
- `geovac/vector_qed.py` imposes V(a,a) = 0 by hard-code (now disclosed in Paper 33).
- `geovac/gravity/warped_dirac.py:322` comment cites the withdrawn 1/(24 pi t) target.
- `tests/test_graph_qed_self_energy.py --slow` did not complete in 35 min (pendant sweep at n_max >= 3).
- Paper 28 vertex gate (`vertex_allowed`, `_so4_channel_count`) is copied into qed_vertex,
  qed_three_loop, qed_anomalous_moment, dirac_vertex, graph_qed_continuum_bridge.

## SMALL findings not fixed on sight

- P2: triple-selection circularity, narrative non-accident, per-manifold attribution, APS premise,
  "correct object", Hopf ontology wording, twelve-mechanism list; K = 137.036064 and 4.8e-7 not
  pinned; loose 8.8e-8 gate; "top ~10" is 11th; p-value not regression-gated.
- P25: n=4 claimed but tested to n=3; L_z test at R=2 only; CP2 pairing multiplicity-blind;
  covariance untested; arithmetic-only hodge1 tests.
- P28: C_VP/C_SE exponent contradictions (n^0.93 vs n^1.48); stale caption "C not computed at
  n_max = 5,6"; Neumann series "terminates by Cayley-Hamilton"; F2 evenness as matrix identity;
  "<0.1%" tracking (0.19%, 0.13%); Connes-vS boundary-mechanism reading; census tier drift
  (charge conjugation intrinsic vs vector); "Sigma Hermiticity" rule 8; Bh conical prop vs Wald
  caveat; graviton DOF assertion; "every stated value backed by a test"; depth bookkeeping
  (Li4(1/2)); "first comprehensive" (deleted); CLAUDE.md/memory citations; NITs (counts, garbles,
  dangling refs, Rosner attribution, Paper 29 for GUE).
- P30: d1 = B^T mislabel; missing tier tags; 1/(4 N_c) and Cartan reduction provable (upgrade).
- P33: slow pendant test; native-Dirac rule tautology (now disclosed in the verdict).
- P36: "structural" denominator wording; Friar -3.842 vs -3.696; 1057.845 MHz unsourced; Eides
  title and section 3.2; 3D Bethe log converges to 0.033% (upgrade).
- P41: S2 table <S(4)> 0.2477 vs 0.2500; R^2 range; W3 band fitted to own data; dim ker(d1)
  mislabel at remaining loci; MK grid 2001 vs 301; W4 "permanent confinement" leading-order
  conclusion; README test count; Guth misattribution; monopole sentence; Perez-Sanchez 2024.
- P51: zeta guards weak (FT1a); L6a envelope band; replica-literal test; cutoff_mellin tautology;
  G4-2 sign pair (scalar-sign SC term with S = -dI/dalpha); upgrades U1 (all odd spheres),
  U2 (scalar a_k SYMBOLIC), U4 (tip convergence MEASURED), U6 (a1 = -4 pi/3 SYMBOLIC).
- Citations: Eides title; Robertson; Luscher T4; Berg-Luscher; Aharony; Ikeda-Taniguchi;
  Rosner/Laporta not two-loop; Fegan-Menotti-Onofri naming; Barrow edition; missing
  Deligne/Glanois/Camporesi bibitems; orphans.

## Guard escape recorded

- C16 `coulomb-ho-calibration-pi-layer` (2026-10-19) stamped Paper 2 as remediated, but Paper 2's
  "pi in K is a specific consequence of the Coulomb S3 Riemannian character" survived: its wording
  matched no pattern. New entry `g5-p2-calibration-pi-coulomb` covers it.

## Added after self-check round 2

- Paper 36 muonium HFS: an itemization attributed to Karshenboim 2005 (all-positive terms "closing" a positive
  +199 ppm residual) withdrawn; the residual is unattributed. Verify the source and build the Mu HFS budget.
- Paper 34 HF-2 row now reads "Schwinger a_e" (was "Parker-Toms a_e"); PM to confirm the HF-2 chain used the plain
  Schwinger value.
- SMALL left: Paper 30 says non-abelian couplings enter at A^3 in one place and A^4 in two; Paper 51 gives the spinor
  coefficient as +(1/12)(1/alpha - alpha) at one locus vs -(1/12)(...) in its abstract (sign convention); Paper 28 c2
  "8-digit" match is about 6.3 digits given the n_int = 50 shift.

## Restructure carry-forward (2026-11-07, PI-approved group5 advisory pass)

Group5 is now Papers 2, 28, 30, 33, 51; Paper 36 moved to group6; Papers 25 and 41 archived
into Paper 30 (rows in docs/retired_papers.md). The certification stamps of this run predate the
rewrite, so every rewritten unit is CHANGED for the next /qa plan.

**PI calls.**
- Paper 2 mechanism count: the historical record gives 9 (Phases 4B-4I) + 1 + 3 in one place
  and 8 + 1 + 3 in Paper 32; Paper 2 now follows 8 + 1 + 3 with the two fibre rows merged and
  says so. Confirm.
- Paper 34's 13th projection is named "Drake-Swainson asymptotic subtraction", but the K-split it
  named is cut from Paper 36. Rename ("Drake-Swainson normalization of Bethe-log sums") or retire
  (changes the count of 28). PM recommends the rename.
- `docs/qa/group5.done.md` criterion 3 and l.118 cite C_VP/C_SE ~ 3/29 as an example; with the
  corrected channel count the ratio is 0.137 and no simple rational. The criteria file is
  PI-confirmed, so the PM has not edited it.
- INDEX QA table is folder-based: Paper 36 now counts under group6 (OWED) though FULL #1 ran under
  group5.
- `docs/qa/group5.done.md` and `group6.done.md` still list the old paper sets; `debug/qa/qa_scopes.py`
  is updated (group5 = 2, 28, 30, 33, 51; group6 adds 36). The done files are PI-confirmed: amend at
  the next /qa freeze.
- debug/data/zenodo_manifest.json carries the old paths and titles; regenerate at the next deposit.
- Seam rows E2/E3 (`docs/forced_free_seam.md`): statuses kept so the gated counts (43 F) don't move.
  E2's P = M1 rested on the cut per-loop 1/(4 pi) reading (PM recommends "none"); E3 now carries the
  Paper 33 pendant-edge value in place of the cut graph-native F2 = 5 sqrt2/3 (keep, or retire the row).
- `papers/INDEX.md`: Paper 33 status set OBSERVATION -> ACTIVE (it now carries thm:physical_rule). Confirm.
- Paper 18 Theorem 2 (axis A3) was rescoped by the sweep to a parity restriction (it had attributed the
  odd-parity clause to gamma^mu, which Paper 33 refutes). Owner theorem statement: review the wording.
- Paper 34's axis-grid argument (~l.2335-2360; tables ~l.3173, 3425; l.3259/3302/3353/3365/9447/9758/9931)
  still makes vector-photon promotion's "1/(4 pi) per loop" the gauge-sector witness of the
  transcendental-only corner; the physical vertex does not carry that factor (it is the model's S^2
  harmonic normalization). Not rewritten: owner call.
- Version: committed as a patch; the PM recommends MINOR (reorganization of the paper series and
  a correction that moves published numbers).

**Owner items found by the citer sweeps (not fixed).**
- Paper 34: hydrogen Lamb rows still tagged Fock o Sturmian o spectral action and "framework-native"
  (~l.4068-4100, 4856-4866, He+ ~5321; group6 synthesis 427-429), feeding the composition theorem's
  "eight tested rows", while Paper 36 says S^3 never enters the number; the muonium -0.6% row and
  "framework total 1047.63 MHz" unverified (~l.4317, 7896); l.4183 "Eides Tab. 7.4 / Karshenboim
  2005" citation unverified (the same one Paper 36 dropped).
- Paper 42 l.730: "axis fixed canonically by alignment with the Hopf-base direction" lost its
  support (archived Paper 25 Sec. V.3); drop "canonically".
- Paper 29: "the only sub-action of the azimuthal U(1)": an m-reflection is in O(2), not U(1).
- Paper 32: Cor. cc_fine_tuning_external says 10^-122 where Paper 51 says ~10^-120; the U(1)
  coefficient 1/4 in tab:unified_gauge has no test.
- Paper 36: the per-level muonic Uehling values it printed (-210.59 / -5.59 meV) do not reproduce
  (driver -219.58 / -14.58; same difference 205.007); not carried to Paper 34.
- Paper 18: the Sec. V Hopf exchange-constant reading of K was cut (it rested on the false S^2
  attribution of B), removing two of Paper 18's Observations. Reviewed by the sweep, flagged here.

**Owed tests (claims kept, no test).**
- Paper 33: `tests/test_paper33_physical_vertex.py` porting stage A (C1-C4, SA1-SA5, the 9j closed
  form, m conservation, Ward-Takahashi on route 2); the L1 cross-sector fractions 20% / 57%.
- Paper 28: c2(n), c3(n), D2, D3, the zeta(3) complementarity (Theorem 5), the S_{r,s} Euler-sum
  tables; D4-D6 pinned only to frozen decimals.
- Paper 30: Z = [2 e^-beta I1/beta]^beta1 at general n_max; the abelian Z and <cos theta_C> =
  (I1/I0)^A; rank Delta_1 = E at n_max >= 4; the three-plaquette 2-cycle geometry.
- Paper 51: B1, B2, B4; continuum tip 1/6 at Bessel zeros; a0 = 1.888; alpha-uniformity CV 1.3%;
  disk propagation 2; the S^5 / S^7 odd-sphere identity; the rank-4 bundle and the total-trace cone
  coefficient (inputs). test_paper51_replica_sbh.py docstring still uses the retired sign pair.
- Paper 2: Delta forms equal only at m = 3; B = B_formal only for d = 3; the S^2 trace 36; the
  8 sign-pattern errors; K's 4.77e-7; CODATA 2022 differences; a reduced deterministic slice of the
  look-elsewhere search (N p, rank 11, 27 at K precision, 3 distinct values); the obstruction
  identities; the twelve mechanism rows.
- Paper 36: 2S acceleration at N = 50 (+2.98%); 2P at N = 55 (+1.24%); residual vs Bezginov
  (-5.64 MHz); LS-8a totals at n_max = 2, 6 and the local slope 2.8; muonic-H driver pins; the
  (2l+1) normalization and I_v(2S) -> 2 Z^4/n^3.
- Registry keys to create (C21): alpha_k_value, alpha_k_rel_err, alpha_cubic_root_inv,
  alpha_cubic_rel_err, p2_search_n, p2_search_p, p2_chance_count; lamb_1loop_mhz, lamb_exp_lp81,
  lamb_exp_bezginov19, a50_2s_mhz.

**Tests whose claims were cut (kept, not deleted; retire or re-home in the guard pass).**
test_paper2_corrections::circulant_hermiticity; test_paper2_ingredients::TestCirculant::
su2_structure_constant_analog; test_su2_wilson_gauge::TestMonteCarloWilsonLoop; several
test_paper41_wilson_witnesses witnesses (W1 Weyl d_W, W2 MK accuracies, W4 large-beta, V5 dual
graph); test_paper25_hopf_grouping (now backs Paper 2's scope section); test_su3_wilson_s5::
TestPaper25S5Pins / TestPaper25CP2Quotient; the LS-7 prefactor tests in test_two_loop_self_energy
(back Paper 32/35 citers only); test_qed_flat_limit, test_qed_vacuum_polarization and the graph-QED
trace tests (test_graph_qed_*), which no longer back Paper 28 text; test_vector_qed's 7/8 tallies
and the rule-6 proxy (now back only Paper 33's footnote e).

## Guard pass (DONE 2026-11-07, v6.6.0) and PI rulings on the carry-forward

PI: "I go with all recommendations." Applied: MINOR (v6.6.0); Paper 2 mechanism count 8 + 1 + 3 kept;
Paper 33 ACTIVE kept; seam E2 P -> none; Paper 18 Theorem 2 wording reviewed (sound; the axis still
carries the name "diagram topology", SMALL); the done-file amendments wait for the next /qa freeze.
Paper 34 owner pass done: projection 13 renamed "Drake-Swainson normalization of Bethe-log sums" (count 28
kept); the gauge sector no longer witnesses the transcendental-only corner (witnesses now the Hopf and Wick
volume slots); the hydrogen Lamb rows are textbook assembly, so the composition theorem has 7 tested rows,
not 8; the muonium row's experiment was Mariam 1982, an HFS measurement (now Mu-MASS, Ohayon et al. PRL 128,
011802 (2022), residual +0.04%, not attributed); "Eides Tab. 7.4" dropped. Chronicle:
docs/chronicles/paper34_owner_pass_cut.md.

Guards: new test_paper33_physical_vertex.py (+ support), test_paper2_owed_guards.py, test_paper30_owed_guards.py,
test_paper36_owed_guards.py, test_paper51_owed_guards.py, test_paper28_euler_sum_derivation.py,
test_paper28_s3_closed_form_bracket.py; four asleep W-pins rebuilt with siblings; every guard fire-tested
(plants recorded in the agent reports, CHANGELOG v6.6.0); independent review 1 LARGE (a Paper 51 citation of
an unguarded panel, closed by a new guard) and 14 SMALL (closed). Owner defects the guards found and the PM
fixed: Paper 2 row 9 ("all-positive" Seeley-DeWitt terms; the Dirac trace has two terms, signs + and -);
Paper 30 Rule-B 2-cycles at n_max = 3 (3342, not the capped 80; Monte Carlo truncation disclosed); Paper 24
CP^2 residual (40.8% least squares, 38.4% minimax, not 24.98%); Paper 51 B2 t = 10 sign change (a finite
alpha-step artifact; exact derivative converges from below to 0.166520) and "1/6 at every t"; small numerics
in Papers 2, 36, 51. Tests whose claims were cut: TestMonteCarloWilsonLoop -> tests/_archive/superseded/; the
LS-7 literature-target and verdict tests -> tests/_archive/diagnostics/; the slow C_2S ratio test re-anchored
to C(4) = 40.78.

**Still owed.**
- PI call: Paper 34 Thm. 1's Uehling row uses the chain fock o spectral_action for Pi = 1/(48 pi^2), while
  Paper 34's own catalogue enters Pi as the flat-space value; if the row has no chain the count falls to 6.
- Paper 34, same pattern, not yet done: the muH, mud and muHe+ Lamb rows still tagged Fock o SA o Sturmian o
  DS "framework-native"; "Eides Tab. 7.3" (~l.3721, 4794, 4917) and "Tab. 7.6" (He+) use non-existent
  numbering; "Eides Sec. 3.2" vs Paper 36's Sec. VIII.A; the 3P Bethe log (+6.3%) has no owner.
- Paper 30: rerun the n_max = 3 monopole fit and the n_max = 3, 4 Monte Carlo on the full plaquette sets
  (994 and 5047); the printed fits are for truncated actions.
- Paper 51: the alpha coefficient of variation 1.3% and propagation number 2 have no test.
- Test docstring re-anchors: test_paper34_projection_spot_checks.py l.46-48, _batch1 l.11, 295-360, _batch2
  l.765-780 (old projection-13 name, per-loop 1/(4 pi) reading); tests/paper36_precision_support/
  precision_catalogue_muonium_lamb.py (Mariam as the Lamb reference); test_paper41_wilson_witnesses cut asserts
  (W1 d_W, W2 MK accuracies, W4 large-beta, V5 dual graph) to trim and re-anchor to Paper 30 Sec. 5.2;
  test_su3_wilson_s5::TestPaper25S5Pins census tests -> diagnostics (keep the n_max = 5 counts, anchor Paper 24).
- TX-A metadata (debug/tx_a_signature_arithmetic.py, tx_a_projection_table.json) still gives Hopf the variable
  alpha, Drake-Swainson K_subtract and a Lamb chain.

## Carry-forward worked (2026-11-08, v6.6.1; PI: "take your recommendations and work the carry-forward")

Every "Still owed" item above is closed. Two of them changed printed results.
- **Paper 30 full-complex reruns** (`debug/p30_full_complex_rerun.py`, settings as printed, caps lifted; JSONs
  archived in tests/wilson_rule_b_support/data/). n_max = 4 Wilson loops on all 5047 plaquettes: joint-fit
  sigma = -0.0043(15), -0.0006(7), -0.0005(1) vs mu = 0.022, 0.0058, 0.0020 at beta = 0.3, 1, 3. The printed
  -0.080(2) at beta = 0.3 was the 1500-plaquette truncation, not the cluster construction; no area law at
  n_max = 4 either. n_max = 3 monopoles on 994 plaquettes and all 3342 2-cycles: c = 22.9, R^2 = 0.935 on
  0.05-0.5, a two-decade drop between beta = 0.1 and 0.2, zero from 0.7; the printed c = 8.95 (R^2 0.981) was
  the truncated action. Paper 30 (5) and (6) rewritten; the 1.196 / 1.242 rates now say n_max = 2. Pins in
  test_paper41_wilson_witnesses moved to the full-complex JSONs.
- **Paper 34 TX-A** (PI call taken on recommendation: Uehling row out, six tested rows). The driver moved to
  tests/paper34_support/ and its metadata matched to the paper (Hopf adds no variable; projection 13 renamed,
  no variable; Lamb and Uehling examples removed). Theorem 3's 16x16 counts were never computed and the driver
  could not return an inverse pair; replaced by the computed 15x15 table (186 commute; 22 one-way, all involving
  Drake-Swainson; 2 conflicts, spectral action <-> spinor lift; 3 idempotent; 12 self-compose). The "metadata
  correctly flags (hopf, spinor_lift)" claim came only from the withdrawn Hopf variable alpha; both orderings
  are now stated as operator-level / semantic, as Weak Minimality already said.
- **Eides pointers** (citation reviewer, primary sources): the 2007 book's Ch. 7 is the muonic-atom Lamb
  shift, so all 19 "Tab. 7.3"/"Tab. 7.6" loci for H Lamb, H HFS, D HFS, He+ and multi-loop budgets were wrong;
  re-pointed to sections (6.1, 3.2-3.5, Ch. 11, 11.1.1, 11.1.3, 7.5, 12.1.5, 12.1.8, 12.2.2) in Paper 34 and
  Paper 23. The "Eides -39.500 ppm" Zemach reference had no source: the closed form -2 Z alpha m_e r_Z gives
  -39.4953 ppm, which the operator reproduces; the "0.012% residual" claim is gone (Papers 34, 23). "Eides
  Sec. 3.2" pointers are correct (Sec. 3.2.1; Phys. Rep. Sec. VIII.A Eq. (40)). Code comments re-pointed.
- **Paper 34 muonic Lamb rows** (agent owner pass, chronicle appended): muH, mud, mu4He+ Uehling / Eides SE /
  closed-form Friar are textbook assembly (TA), finite-size rows L2; Uehling chains "entered Pi o Fock o
  rest-mass". Also the Ps structural reading, the He+ Foldy status (FN -> L2), "proton" (not deuteron)
  polarizability for muH, and the muH row's stray "flow" class.
- **3P Bethe log** reproduces by the 2P/3D route: -0.0405950 at N = 40, +6.30% vs Drake-Swainson; guarded.
- **Paper 51**: the 1.3% CV is at t = 10 (both routes, 1.35%); 0.23% at t = 1; paper now states the t. Both owed
  measurements guarded (propagation number 2 ported from debug/l6_disk_prop2.py).
- **Test re-anchors**: Paper 34 spot-check docstrings; muonium driver to Mu-MASS 1047.2(2.3)(1.1) MHz; Paper 41
  witness file trimmed of Paper 41-only asserts (d_W, MK accuracies, large-beta sigma, dual-graph counts) and
  anchored to Paper 30 sec:ruleb; the ST-SU3 census -> tests/_archive/diagnostics/test_su3_s5_plaquette_census.py
  (the n_max = 5 counts stay live as TestPaper24S5GraphCounts).
- C16 `g5-carry-forward-2026-11-08` (fire-tested: 11 retired forms fire, 13 corrected/near forms silent).

**Still owed (new, found on the way).**
- Paper 34: figures attributed to the old "Tab. 7.3" (D HFS Zemach -150 ppm, the +12 to +18 ppm and -18.5 ppm
  H HFS bands) are not verified against any table; the book's list of tables was behind authentication.
- Paper 34: "Eides Eq. 2.35" (7 loci) for the Friar/Foldy closed form; the book's Ch. 2 is six pages of general
  features, so this is probably wrong (UNVERIFIABLE without the book). geovac/two_loop_self_energy.py comments
  cite "Eides 2001 Tab. 7.3" for C_2S.
- Paper 34: the "Three-projection chain" heading still holds the muH 202.17 row; the "III.14 transports the
  Eides-bracket machinery verbatim" sentence in the mu4He+ Z-scaling paragraph is the same overclaim, mild.
- Paper 34 observation obs:axes_correlated cross-tabulates 28 projections by hand; the computed table covers 15,
  and in it Drake-Swainson is also a transcendental-only projection (flow tier), which the observation omits.
- Paper 30 abstract and the group5 synthesis say the area coefficient is consistent with zero "at n_max = 2";
  true, and n_max = 4 now agrees in substance (the beta = 3 value is -0.0005(1), small but not zero).

## Follow-ons worked (2026-11-08, v6.7.0, MINOR)

The "Still owed (new)" list above is closed except as noted. Primary sources read (Phys. Rep. 342 PDF text; Friar-Sick
arXiv record); full account in CHANGELOG v6.7.0. The Paper 34 Calc-rZG fit is withdrawn (circular H point, unsourced D
budgets); D HFS re-sourced to Friar-Sick 2004 and Eq. (414); the 21 cm residual is unattributed with sourced
magnitudes for scale; "Eq. 2.35" -> Eq. (158); the LS-8a code target -> B_40 = 0.53895; the "Three-projection" row
and the mu4He+ transport sentence fixed; the axes observation left as is (Drake-Swainson's tag in tab:projection_axes
is "none", so it is not a transcendental-only witness; the driver's stale "flow" tag was the defect).

**Still owed.**
- "Friar moment" used for the leading (Z alpha)^4 <r^2> term (~15 loci in Paper 34); SMALL terminology sweep.
- A Zemach-radius extraction with sourced Layer-2 budgets; a complete hydrogen 21 cm budget.
- V.D.N exposure IDs now skip V.D.1 (kept stable rather than renumbering ~15 loci); Pattern C narrative at n = 5 needs
  an owner re-read at the next /qa.
- TX-A driver tags vector_photon as transcendental-only (one_over_4pi), which the paper's observation rejects; the
  driver's axis printout is not cited, so left.
- debug/ top-level file count (2680 vs 600): archive sweep is a PI decision.

