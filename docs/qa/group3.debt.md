# group3 — debt (SMALL findings, owed guards, completeness gaps)

Owner of this list: the PM. SMALLs here are fixed in one scheduled sweep outside a cert cycle
(edit bar, `docs/qa/criteria.md`). Guards are written as their own activity and fire-tested.
Source run: `/qa group3` FULL (incremental, first cert), 2026-11-03, snapshot 7727c540
(`refs/qa/group3-full-2026-11-03`).

## Owed guards (separate activity; fire-test each against the wrong answer it excludes)

These are OWNER-LARGE "test proves less than the prose / no test" findings. The prose was not
changed where it is true (re-derived by the reviewers); the owner units stay unstamped until the
guard exists.

- **Paper 54** `thm:angular_structure` and `tests/test_paper31_two_body.py::_angular_audit`: the
  Gaunt / k = 0 checks pass for any direct-diagonal tensor (plant `if True:` did not fire). Rewrite
  to assert the support is exactly {[a,b,a,b]} (25, 196) and that the 82 / 3604 Coulomb-only
  elements are absent; drop the "Gaunt" label.
- **Paper 54** `thm:free_factorization` + `tab:heat`: no test of {D, gamma} = 0, the heat-trace ratio,
  the table values, or "4 of 14 / 29 of 55".
- **Paper 55** `thm:jlo_depth2_reading_A` (40 zero residuals; Paper 56's Reading A rests on it),
  `thm:na1_diagonal_collapse`, `thm:m3_gammaP_closed_form`, `thm:na1_trace_functional_collapse`:
  debug-only backing. Reviewer re-derived M3^{gammaP}(s) = 2^{2s-3}(beta(2s-1) - beta(2s-3)) at
  s = 2..6 to 1e-70.
- **Paper 55** M1 / M2 witness panels are tautological (plants 17 pi^8 -> 16, 4/pi -> 3/pi did not
  fire); S^(3) definitional tie only ~4 digits (independent extrapolation 31.5725612075).
- **Paper 56** C1-C3 checks cannot fail (C1 plant F4 did not fire; C3 is det g = 1 on the panel);
  claim-test matrix row 119 "BACKED-SOUND (real det computation)" is wrong; no test of the dual Hopf
  map. PS-4 End_compat = block-lower-triangular has no test (exact solve confirms 4, 19, 55, 125).
- **Paper 57** 98.3% / Table I: no test; `debug/principle_hunt_audit.py` hard-codes a stale catalogue
  (B6 128, H2 11.10 x Q, C9 outside). Gate the table against `docs/forced_free_seam.md`.
- **Paper 61** (open since 2026-09-07): the weight <= 3, 64-digit negative has no test (reviewer:
  negative holds to height 100, control recovered); `test_paper61_stokes_amplitude.py` covers one point
  at 3.5% vs the prose's four at ~4e-4; `test_paper59_coarea_reduction.py` asserts 5e-3 / 1e-8 vs the
  prose's 2.2e-9 / 5e-29; `test_paper59_diagonal_A` PSLQ basis omits 1, K0, K1, E1 (full ring negative
  to height 100).
- **Paper 18** `thm:ac_factorization`: only a debug toy; `eq:dirac_modular_residual` guard asserts
  dict literals (plant 2 m^2 -> 3 m^2 did not fire; formula true to 1e-40). Upgrade-ready: kappa-B
  weights {7/16, 5/8, 3/8, 5/8, 17/32, 11/32} sum 63/4; alpha^2-Ihara degree 44 / -35/648.
- **Paper 22**: assert radial_zero == 0 and density == D(l_max) in the potential-independence test;
  n_max >= 2 independence; sign-holonomy numbers looser than prose (13.3x asserted > 10x; Delta 43-48x
  and 2.5-3.2x unpinned; n_max = 4 log-only).
- **Paper 24**: assert 56 nodes / 165 edges at N_max = 5.

## Routed to another target

- **Paper 28 Theorem 3 "Deligne-Glanois descent from level 4 to level 2"** (owner, group5), restated
  in Paper 55 (L237, L644, L1034, L1164, prop:vertex_parity_descent, abstract, conclusion), Paper 18
  L1160 and Paper 32 L2481. Reviewer CL-B: no level-2 object appears in the vertex-parity
  decomposition (parts are level 4, the sum level 1). Not a self-contradiction (4->2->1 = 4->1), so
  not graded LARGE here; the group5 run must check the reading at its owner.

## SMALL sweep (edit bar) -- OPEN

- **Paper 18.** Extra "tiers" vs the six-tier headline (L714, L1905, L2281; chi_-4 "calibration-tier"
  L761 vs "beyond calibration" L2281); algebraic numbers called "transcendental" (L717, L775); "only
  when non-compact" (L116 vs L1331); kappa "calibration constant" (L2855, L3139 vs L3375); Class S vs
  P (L3059 vs L3252); Paper 7 equivalence widened (L3273); untiered L1756, L1761/L3352, L1785 (v2.9.1
  graph eigenbasis, check hop era), L2361; TS-E3 (L1182); NITs L294, L701, L1779, L2446 vs L2487;
  L2662 "structural prediction"; abstract strong-form headline (self-hedged; standing since the
  2026-09-14 scope note); L3693 "validated through Level 5" vs L570. Tests: d_even4/d_odd4 cite
  should point at `tests/test_spectral_chi_neg4.py`; "35/35 tests" (L2273) now 50; aabb no-pi guard;
  elliptic-kernel backing -> Paper 59 tests. Citations: goncharov1999 (3-manifold volumes are K_3 /
  zeta_F(2)); perez_sanchez2024 says Yang-Mills-Higgs (drop at L1373); Rovelli quotes are from the SEP
  entry (cite it), not quant-ph/9609002 (CLAUDE.md WH7 status corrected 2026-11-03, v6.1.4); rejzner2016 a survey;
  liu_noui "H. Liu"; adamchik2003 -> Z. Anal. Anwend. 21 (2002) 817-826, uncited; karamata1962
  reprint unverifiable.
- **Paper 22.** Theorem 2 converse heuristic; no tier tags on Theorems 1-3 (upgrade: Theorem 3
  densities are exact, independently re-derived); L423 jj "identically" vs L528; NITs 43-48x vs
  43-50x, 2.5-3.2x vs 2.44, M_L = 0 vs singlet space. Citations: KlassenTerhal2019 is not the hardness
  result (Marvian-Lidar-Hen); Hangleiter2020 does not define the stoquastic gap.
- **Paper 24.** L1105 "both biconditional ... standard" vs L1047 caveat; Theorem 2 untiered; L754
  calibration pi in K without Observation label; Wilson layer unscoped (L773, L785); Layers 5-6
  untagged. Citations: hall1994 credit (L1066 -> Bargmann 1961); iachello2006 chapter.
- **Paper 31.** False displayed equation D^2 = Delta_LB + 3/2 on the scalar sector (L248-255; the
  Lichnerowicz term is spinorial); L786 "same reason" / plural heading L769; abstract L59-63 Wilson
  "decided" drops "argued"; pi in K as fact (L112, L857); two-body numbers untagged, Paper 54 uncited,
  vs L179 "no new computations"; nuclear coefficients untagged (L668); CLAUDE.md / debug cites;
  "decreases across all four weightings" (only Dirac tested; Laplacian at n_max = 2 bound 0.15);
  "18-orbital set" vs test's 17.
- **Paper 54.** "Six observables" vs bib "Five" vs corpus 7; L390 "Fock projection maps the graph onto
  S^3" (section 1.5); section 7.2 "bit-exactly" for 1e-14, NO-TEST; 27.90 not `\gvq`-annotated; no tier
  tags; CV 16.46 / 20.32 unasserted.
- **Paper 55.** L1963 "three theorem-grade" incl. a null PSLQ; abstract "S_min fully resolved" (PSLQ
  220 dps, conjectural bases); Thm NA-1 needs per-chirality (L1272); F(s) shift-parity vs tie-break
  (L2756 vs L1514); depth of pi^2k 0 vs 1; JLO "Theorem" on spot checks; "closes" / "theorem-grade
  closure" / "L1 is forced" on finite panels; "candidate" dropped for U* (L2842, 2919, 2957, 3133);
  A1-Matsubara open vs closed; abstract lists archived Paper 43; S^5 / Yukawa listed open while
  closed; "irreducible" L1163; Delta -> M3, F -> M2 unhedged; memory/ cite L1468; D_even labelling
  n_Fock vs n_CH (L326, L813, L1039 vs L766/L907). Citations: Hoffman "Table 2"; FM "Theorem 7.5" is
  Prop. 7.5; "Lemma 7.1.5" = Lemma 7.1(5); karamata1962 reprint; CCM "explicit formulas" overstated;
  FM "admits MZV" is our reading.
- **Paper 56.** Reading A vs B "open" (L1474) vs "A wins" (L2139); equality conjecture's depth >= 2 vs
  gr^1 (L2076 vs L1428); "rank 2 per weight"; M3 rank presented as computed (L1271, L1566) though set
  by the Definition; no theorem tier tags; L1084 conjectured vs L1042 uniform; KO-3 owner (L1659);
  stale test docstrings (closed immersion; rank-1 skip reasons; "GENUINE" C1). Citations:
  Cartier-Milnor-Moore not in Brown section 5; glanois2015 -> "groupoid", 2016.
- **Paper 57.** "seven non-selection theorems" vs the tally (L598); Family 1's 17 double-counts;
  "Accuracy" column header with P5 bolded; L113 "motivates" vs Paper 0 "counts"; 9 debug-memo
  bibitems (C14 advisory).
- **Paper 61.** "transcendentals placed" universal (7/e, K0(1)/K1, 1/sqrt pi unplaced); L128-131
  "so the arithmetic lives" upgrades a bounded negative; stale W0 clause (L617 vs L619); resummation
  unconditional (L503 vs L481); abstract Eisenstein/CM untagged; "tier tags carried over unchanged"
  false; L387 2.2e-9 and L369 1e-32 vs tests; NITs L425 vs L453, L652. Tests: b5 literal;
  bd_jacobian 1e-20 vs 1e-32; Euclidean corner; theta decay identity; "lemniscatic" variable;
  maximal chain unchecked. Profile `paper_61.done.md` stale: C16 now has p61-* entries; 12 matrix rows,
  not eleven. Citation: broadhurst2016 uncited (PI note).
- **Synthesis.** "ecosystem" (L876); nonlinear stereographic (L1037); EMN converse (L1440);
  cusp-form (L1481); "generator" of pi-content (L1683); integer spectrum (L440); 98.3% labelling date
  (L1530); Wilson "argued" (L1307); Paper 7 PRIOR ART (L380); Paper 57 paragraph omits W and
  Observation entries; date L29; "Four more recent papers" (L1346).

## Self-check (2 rounds, closed 2026-11-03)

Round 1: 6 LARGE (P18 L505 mu transcendental; P22 L647 wrap-up; P31 heading L1226 + Cor 1
restatements; P54 heat-kernel explanation x2; P56 surjectivity onto G_4 / G_m and Ga^inf x| SL_2
identified with a sub-quotient of G_4; synthesis "derives selection rules" + heat-kernel reason) ->
all fixed. Round 2 (last): 1 SURFACE-LARGE (P31 L1460 "correct selection rules from both routes")
-> fixed after the last round, so that unit is unreviewed and stays unstamped. Round-2 SMALLs, open:
- **PI item / possible OWNER-LARGE (Paper 56, pre-existing):** Definition `def:period_map` gives
  pi: H_GV -> periods, but a homomorphism U* -> G_4 dualises from a Hopf map O(G_4) -> H_GV; and (C2)
  checks only depth-1 coproduct compatibility, which supports at most a map into U_4^ab. Check against
  the primary math (Deligne-Milne Thm 2.11 direction) before the next group3 run.
- P56 codomain stated as G_4 in the theorem (L1257) but U_4^ab elsewhere; L1419 "sits inside gr^1_D
  U_4" (a quotient); equality conjecture "iff every depth >= 2 MZV is realised" against the depth-1
  target U_4^ab (L2080, L1314; synthesis L1449 "each weight"); L240 / L1157 "motivic Galois group on
  the mixed-Tate period ring" should be its abelianized unipotent part; L1478 Reading A vs B "open"
  vs L2142 / L2198.
- Seam doc L183 "all at outside" vs E6, I1 (M1+M2); L187 lists 12 vs 17; P57 L236 same; P57 table
  "Family 1 \checkmark" for P2 though E6, I1 are missed.
- P18 L148, L2192 "operator order determines the motivic weight" unscoped.
- P31 table L654 Casimir under the Coulomb block; P31 L433, L516, L537 "unique central potential"
  unconditional beside the [CONDITIONAL] theorem.
- P54 L376 "determines the multipole hierarchy"; L92-100, L391 "now verified" / "explains"; L495
  "derives which nuclear channels couple".

## Completeness-critic gaps (2026-11-03)

- K = pi(B+F-Delta): the claims reviewers enumerated every K locus per paper (CL-A P18, CL-B P55, CL-C
  P57, CL-D P24/P31, SYN synthesis), no prohibited tier word; the critic's list was against the PM's
  summary ledger, not the reports. C5/C12 PASS.
- Calibration-pi status clauses: P24 L754 ("The Coulomb S^3 side carries calibration pi", untied from
  Coulomb two sentences later), P18 L906, P31 L712, synthesis L526 / L953-955 / L1168: not
  clause-read for the withdrawn Layer 2; P24 L984, L1183, section 1120 bodies unread.
- Claim-matrix row 88 (P22 footnote) was stale: re-checked and rewritten 2026-11-03 (routine computes
  the global D; `test_angular_zero_count_computes_global_headline_density`).
- Cross-paper pairs not compared: P24 Cor L1228 vs P18 L1946 (bundle type vs spinor diagnostic; P18 now
  scoped); P24 L740-752 six-tier list vs P18; P56 rem:paper59_cm (L1834-1880) vs P61 and the section-3
  "Q(i) seam as a shared mechanism" negatives; P55 S^5 sections (L1077, L1872, L1940) vs P24; P55 vs
  P56 level statements (rem:why_level_4 L1397).
- Unread regions: P57 L408-606 (Bohr-site), L156-181 (J doors), L661-689; P31 L877-943 (TS-C),
  L1015-1225 (four-bucket, blockers, falsified (3,1)); P24 L430, L677-736, L1257-1429, L1460; P22
  L496, L521, L739-906, L908; P54 L436-509; P55 L331 (case-exhaustion), L1489-1561, L1607, L1793,
  App. L3254; P56 L796-1160 (TC-2a-2d) and remarks L470, L577, L622, L776, L872, L1073, L1414, L1433
  (honest-scope direction), L1482, L1731, L1768; P18 L227-500, L789-983, L3424-3540; P61 L134-360;
  synthesis L785 ("Quantum Dynamics on the Graph", at risk after the graph->S^3 refutation).
- Tables unread: P24 table* L574; P31 L513, L635 body; P56 L483 ("1,771 bit-exact zero residuals",
  same class as the 5,864 headline); P22 L350 body; synthesis L705.
- P18 L3162-3163 "eigenvalues, degeneracies, kappa -- integer or rational": defensible in context (band
  edge 8, label Casimir); SMALL wording risk under C16 `graph-eigenvalues-integer-or-rational`.
- C3 inline tiers nearly absent: P56 0, P55 4 (~3,300 lines), P22 2, P54 2, P57 6, P31 5.
- Matrix rows absent: P24 (Bargmann rationality L366, Prop L656, Cor L1228, Cor L1362, KS frame); P55
  (case-exhaustion, thm:m3_cyclotomic_mixed_tate, descent, JLO); P54 (free_factorization, Prop L338,
  Cor L401); P31 (signature partition); P56 (TC-2a-2d, rem:paper59_cm); P18 (ac_factorization,
  eta-trivialization).
- Docstrings unread: geovac/nuclear/potential_sparsity.py, nuclear/bargmann_graph.py,
  extended_tapering.py, tannakian.py, pro_system.py, operator_system.py.

## Upgrade candidates (PI, tier changes)

- Paper 22 Theorem 3 densities: exact enumeration, independently reproduced -> INTERNAL THEOREM.
- Paper 61: weight <= 3 negative and the cusp coefficient A negative both hold to height 100;
  lambda = 1 - rho also for rho > 1.
- Paper 55: S_min closed form via Hoffman's t-value parity theorem.
- Paper 18: kappa-B identity guard; alpha^2-Ihara degree 44 / -35/648.
