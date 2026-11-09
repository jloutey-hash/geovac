# /qa group6 seeded FULL 2026-11-08 -- running findings log (PM)

## Calibration
- CODE-35: CO-35a caught (graded NIT, abs() redeemed by L149, fire-tested), CO-35b caught (D7). 2/2.
- SYNTH: SY-1 caught (STALE #2), SY-2 caught (STRONGER #3). 2/2.

## Genuine findings (UNVERIFIED until PM check)
### CODE-35
- D1 OWNER-LARGE: P35 L790-807 "208 checks, 100% match, load-bearing principle" backed by tests asserting what they construct; tally hard-coded (predictions:198-211); no iff-direction test.
- D2 OWNER-LARGE: Heisenberg-Euler pi content convention-dependent (HL units: 2 alpha^2/(45 m^4), no pi); "Dunne 2012 Eq 2.3" UNVERIFIABLE; 1/pi form matches neither convention.
- D3 SMALL: SB guard pins wrong sign (predictions:50, :58).
- D4 SMALL: Dirac degeneracy label 2(n+1)(n+2) vs prefactor 4 (P35 L1224; module s3_dirac_spectrum).
- D5 SMALL: KG panel self-constructed; "bare-graph ring" names wrong object (n(n+2) is continuum).
- D6 SMALL: "iff m^2 in same ring" false (m^2=sqrt2 counterexample) P35 L255-256.
- D7 SMALL: module docstring "1/2 zeta(-3) = -1/240 ... sign depends" false (thermal_tensor_triple.py:394-404).
- NIT: falsification section concedes K open counterexample vs "no exceptions" L793.
### SYNTH
- S-1 LARGE: Paper 36 absent (synthesis "four Group 6 papers" :47, :134).
- S-4 STALE: ERI fill 42%->100% (:384-386) vs owner 17.1%, "graded".
- S-5 STALE: core-valence MI 2e-3 / Z>=5 (:398-399) vs owner 4e-3 / Z>=7.
- S-6 STALE: gamma_inf 1.93-1.96 (:346-347, :508) vs owner 1.96-1.97.
- S-7 probable STRONGER: "any reduced density matrix" (:294-297) vs owner "single-particle".
- SMALL: N_Pauli=27.90xQ unannotated (:389-390).
- Owner notes: P34 says focal-length program "advances" (L4720, L8586) vs CLAUDE.md 1.8 RETIRED -> PI; P34 L7543 lists Karshenboim among exposures while row withdrawn; multi-focal six vs seven instances; P35 L799 "load-bearing principle" vs L863 "structural Observation".

### CODE-34 (seeds: CO-34a caught as D10 NIT, CO-34b caught as D9 NIT; also caught CL-B as D2) 2/2
- D1 OWNER-LARGE: Prediction 1' L2-count=1 class (tex 9339-9340) rests on H 21cm (+18.4 ppm now itemized as omitted textbook terms, r_Z moves 0.34 ppm) and D +23.0 (chain-truncation); test_paper34_l2_presence_bound.py excludes D BF row for chain truncation (:123-134) but keeps H 21cm row (:111) hidden by fitted L2_CLASS_SCALE (:168).
- D3 SMALL: H21cm reduced-mass printed 0.998366 (tex 4817) vs 0.9983679 -> chain 1420.42907 vs printed 1420.4318 (1.9 ppm); chain guard tolerance 5e-3 MHz too loose.
- D4 SMALL: "operator-level" Zemach = definitional identity (magnetization_density.py:567, profile width set so M_1=r_Z); "machine precision" by construction; matrix rows 505, 932.
- D5 SMALL: r_Z(p)=1.00(11) uncertainty omits kernel/Table XIX-r_Z systematic (~0.07 fm).
- D6 SMALL: "remaining +7.2" is a plug; "term by term" wording.
- D7 SMALL NO-TEST: D +23.0 vs +25.3 not tested; omits D recoil.
- D8 NIT: magnetization_density docstrings 0.10124 stale; -6657 vs -6665.2; eq 95/96.

### CODE-36 (seeds: CO-36b caught (D2); CO-36a sign-blind abs MISSED -> PARTIAL 1/2; also caught CL-E as D1) -> re-dispatch focused on S3 class
- Verified: all 2S literals vs arXiv REVIEW.TEX; 2P formulas; Table XX at 0.805/0.862/0.891.
- D3 SMALL: budget's source-anchored guard (Table XX closure, 2 kHz) weak: mutants pass it, caught only by self-pins; add 0.805/0.891 rows + per-formula checks.
- D4 SMALL SURFACE: matrix missing P36 rows (erratum, A50, 2P/3D logs, 1S crossing).
- N1: FNS scales rounded 145 (closed form 145.25 -> 1057.8247); N2: "every term the review tabulates" omits ~5 Hz terms; N3: Barker-Glover not in Table VIII; N4: two_loop_self_energy docstrings "derive C_2S / projecting"; N5 lamb_chain docstring.
- U1: two-loop divergence provable (nonneg summands, g/lambda^2 -> 2) -> internal theorem.

### CODE-26 (seeds: CO-26a caught F5, CO-26b caught F4) 2/2
- F1 OWNER-LARGE: "essentially unique / isolated sparse basis" false: global SO(3) Wigner-d rotation (mixes m) keeps ERI count 107/625 exactly; real-harmonic mixing 109. Abstract L51-56, III L298-316, L402-407, concl L753-758.
- F2 OWNER-LARGE: O ground multiplet 9-fold (3P), not 7; C 9-fold; Lanczos k=8 artifact; paper L471, test L321 pins artifact; "combinatorial not LS" false.
- F3 OWNER-LARGE: hub migration NO-TEST; C irreproducible (2s 4.03% vs 3.7%, top pair 0.941 vs 1.393); B/C member-dependent undisclosed.
- F6 SMALL: H2 bond R-independence by construction (MEASURED tag overstates); 40x needs pk_in_hamiltonian=False (default gives ~57); "different electron counts" explanation false.
- F7 SMALL: eq:icv pairwise vs table bipartite; He 0.751 reproduces in neither (0.7638/0.7633).
- F8 SMALL: stale test docstrings (many); F9 SMALL: matrix P26 rows stale (233, 235, 236, 245) and FREEZE C8 P26 block stale (0.3033139/49.6/4/7/6).
- F10 NIT list. U1 core closure N/O/F at n_max=2 provable by L,S symmetry; U2 107/625 exact Gaunt count.

### CLAIMS-34A (seed CL-A L2922 "proven for every spectral triple" MISSED -- listed 'proven (L2922)' as checked; bonus: flagged out-of-range CL-C L9632) -> UNCALIBRATED, re-dispatch
- L1 OWNER-LARGE (C5-adjacent): L3072-3075 "graph-derived ... the foundational alpha" vs L157, L3031-3037.
- L2 SURFACE-LARGE: L1188-1192 "confirming PK is not the bottleneck" inverts Paper 17 L1172-1175 (and P34 L1878-1879).
- L3 SURFACE-LARGE: -(11/48) attributed to Karshenboim 2005 S9 at L831, 873, 3086, 3309 vs own citation_status (framework's reduction); "sub-MHz closure" circular (calibration input); 0.54% not 0.57%; 14 MHz not sub-MHz.
- L4 SURFACE-LARGE: Bethe log "flow tier" L680, 3152, 3264-3267 vs Paper 18 (calibration; flow not a tier).
- SMALL: S1 L558-562 alpha^n finite parts framework-native vs two-loop divergence; S2 L2902-2903 PSLQ "not algebraic combinations" overclaim; S3 L1510 uniform precision false; S4 L1148-1150 balanced+PK mislabel; S5 L1404/1421 vs L1316; S6 L2057 vs L2093; S7 V/D correlation inconsistency; S8 "only" superlatives L3272, L1022; L389-391 "match experiment exactly"; L3030 WH8 cited as authority, L3084 WH5 for Yukawas; abstract L41-43 single-origin unscoped.
- NITs incl L410/L3548 "sixth" multi-focal, L2372 Lamb precision wording, L515 LS-2 2-3%.

### CLAIMS-35/36 (seed CL-E caught as L1) 1/1
- L2 SURFACE-LARGE: P35 L1340-1349 "structurally ready for two-loop QED ... additive rather than blocked, ~8 sprints" vs P36 L386-391 divergence; also self-contradicted L1387. SMALL: L1360-1364 "LS-6a predicts ... 5.8x reduction" drops "not a GeoVac computation".
- L3 OWNER-LARGE: P35 rem:graduation L796-802 + concl L1442-1444 "graduates ... load-bearing principle", 208 checks inflated (200 trivial) vs L862-867 Observation. (same as CODE-35 D1 -> converged)
- S1 iff false (L254-256) [=CODE-35 D6]; S2 "Hurwitz 3/2 drives all three Casimirs" false for scalar (L755-760, 923, 966, 1450); S3 tally ten vs 12/10 (L956-966); S4 "Paper 34 is not amended" stale (L610-617, L1427); S5 projection-depth + broken \ref{sec:open} (L822-823); S6 refined Prediction 1 overreach (L1389-1396); S7 "newly computed -31/60480" and -17/3840 priority unverified (L429, L752, L1446).
- NITs: P36 6.9 vs 7.1 "alone"; title "Convention Erratum"; P35 phantom XIV.A/B; 7pi^2/720 per-Weyl; etc.

### CLAIMS-26/27 (seed CL-D caught as L1) 1/1
- L2 SURFACE-LARGE: P27 abstract L44 "equally large edge ... carries none" vs L540-541 (|V| half), L1247 "comparably large".
- L3 OWNER-LARGE (verify): P27 sec:scope L202-234 Li/Be examples are hop-model (EP-2j) values unlabelled; physical H1 has 2s/2p degenerate -> "open-shell non-degenerate" violation mode collapses.
- L4 OWNER-LARGE: P26 "essentially unique / isolated" (converges with CODE-26 F1).
- SMALL P26: S1 L643-647 40x "consistent with Z^-2.52" (predicts 16) + false electron-count excuse [=CODE-26 F6]; S2 L493 "theorem" vs provenance; S3 L47-50/668-670 "<0.2%"/"no two-body" vs L255 7%; S4 L663-665 unqualified separable->no entanglement; S5 L57/389-391 Z-independence automatic; S6 L723-724 "diverges at dissociation" vs ln2; S7 L396-407 cutoff dependence undisclosed; S8 L772-777 attribution.
- SMALL P27: S9 L1046/1105-1107 "bare-graph ground state" entropy; S10 L1101-1103 wrong concavity bound; S11 concl L1172-1173 limit misdescribed; S12 L803-805 vs L990-996 composed vs balanced; S13 L583-585 "irrational" tagging; S14 L36/432-433 "rigorously recovered" circular; S15 title lacks qualifier.
- NITs list (graph eigenbasis, P27 -2.56 stale L840-842, etc).

### CITE-B (seeds CI-B1 Kempf caught F3, CI-B2 Bezginov interval caught F1) 2/2
- F2 near-OWNER-LARGE: P35 L1238-1239 "17/480 matches textbook value of Camporesi-Higuchi 1996 (Eq. 5.27)" -- Eq 5.27 is a root-system measure; no Casimir in paper (arXiv gr-qc/9505009). Other two supports UNVERIFIABLE.
- F4 SMALL: P36 L182 ln k0(2P) = -0.0300167089 vs Drake-Swainson (via Eides Table I) -0.030016708630213 -> -0.0300167086.
- F5 SMALL: P35 L1354 "Eides 3.2" pointer nonexistent.
- SMALL: P27 L1357 Moshinsky1996 bibitem uncited.
- UNVERIFIABLE load-bearing: Dowker-Critchley PRD 13 3224 for 1/240 (title is de Sitter!) -> PM check; Bytsenko 4.5; Dowker-Altaie value; Lieber states; Huff SO(4,2); P36 CODATA constants unattributed.

### CLAIMS-34A re-dispatch (L2380-3616; caught CL-A as L1) -> P34A recalibrated (class: status/scope overstatement)
- L2 = alpha graph-derived (converges with CLAIMS-34A L1).
- S1 PSLQ disjointness L2901-2903, 2830, 3055-3056; S2 BO composition "exercised operationally" L2983 never run; S3 "all four faces" sympy L2502 (only 2); S4 "every spectral-side transcendental" vs 15 L2831; S5 rest-mass "only" L3272; S6 Breit dimension energy vs dimensionless (counts 12->11, [E] 3->2); S7 corner classification L3333-3339; S8 "Lorentzian extension has since been carried out" L2545-2547 vs open/withdrawn; S9 "bit-identically" with 3e-10 L2881.
- PM-check: L3092-3094 "calibration pi in M1" (calibration-pi layer withdrawn 2026-10-19) STALE?; L3475-3477 "empirical separation across Lamb/21cm/D" may rest on withdrawn rZG fit.

### CLAIMS-34C (L6501-10261; seed CL-C caught as L1) 1/1
- L2 OWNER-LARGE: He hop zombie L7308, 7335, 7343 ("~30" vs ~3), 7346.
- L3 OWNER-LARGE: convention-exposure "unique diagnostic instrument" framing (L7536-7547, 8516-8519, 8555) contradicted by every row.
- L4 OWNER-LARGE: pred:l2_presence failure forms omit direct violation; D row violation explained away (L9312); Ps epsilon inconsistent (L9275); "sharpest test" cannot fail (L9319-9323); muH -7227 ppm row omitted. (converges with CODE-34 D1)
- L5 SURFACE-LARGE: -(11/48) "canonical Karshenboim closed-form" L6719 (converges with CLAIMS-34A L3).
- L6/L8 OWNER-LARGE: muH HFS "64 ppm convention split resolvable / above any plausible residual budget" L7413-7418, 7428 vs -0.72% (7250 ppm) residual.
- L7 SURFACE-LARGE: L9457 "calibration pi" live vs Paper 31 withdrawn (and P34 L9678).
- SMALL: mud arithmetic L7889-7892 (+0.12% -> ~0); Lyman-alpha sqrt2 called transcendental L7139/7151; Breit tier L9215 vs L6683; W1e list L8443-8448 counts superseded HF-3/4; counts L7548 eight exposures, L10228 thirteen vs 14; V.D.3 live at L7973, 8639; alkali cliff reason false L8546-8549; Ps closure by construction L6610; muHe "close the gap" L6923, 7064 vs -1.44%; L6534 "muH 1S HFS +2 ppm" mislabel; H 2P FS overclaim L7170, 7232, 7218-7221; He L7349 unbacked; Lorentzian "no structural cost" L8830, 8898, 9022 vs L8958-8970; conclusion L10115 universality vs L10181; L9318 +18.4 comparator (should be +25.3).
- NITs list; abstract "depth-3 -211" vs rem "-456" PM check.

### CLAIMS-34B (L3617-6500; seed CL-B caught as S10) 1/1
- L1 OWNER-LARGE: He+ Lamb autopsy residual misattributed (L5194-5246, 5280-5290, 5325-5327, row 4273-4278): omitted one-loop alpha(Za)^5 scales Z^5 -> ~+228 MHz at Z=2 accounts for -200 MHz deficit; caption "closed by Drake-Yan" vs -181.55 residual; "-0.55%" stale.
- L2 SURFACE-LARGE: catalogue rows tag textbook Uehling/SE with Fock chains (L3699, 3732, 3736-3738, 3752-3753, rem 3903-3904) vs autopsy "no GeoVac projection enters".
- L3 OWNER-LARGE: FN label definition (L4714-4716) violated by BF / reduced mass / Schwinger rows in 21cm, D, T, He-3, Li-7 tables; HF-2 "graph-native F2" vs "taken as input" L4086-4088; "framework-native subtotal" headlines.
- SMALL: T "inside +12+-5" false (L5740, 5791, 4220); Mu HFS itemization doesn't add + misplaced caveat (L4291-4304) + C-REF 4.703; Lyman-alpha residual = reduced mass (L3849-3852); Hylleraas-Eckart honest-scope reversed (L4638-4640 vs 6280-6290; 6266-6271); Schwinger row stale "33%" (L4003); HF-4 "no magnetization operator" stale (L4097); He 2^3P "same architecture closes H Lamb" stale (L5447, 5344-5345, 5422); He-3 reading stale (L6094-6100); D/T overclaim (L5692, 5884); D table ppm conventions mixed (L5565-5586); 2S Bethe log as match (L3784; 4011, 4014); row 3651 "Fock + kappa"; C-REF: mud -0.50% vs Krauth route mismatch (L5079), muHe 1666.305, muH FS -4.3% irreproducible (L4955, 4242), Ps HFS Ishida 203.3942 not 203.389 (L4130), Mu 1S-2S ref; Cs F_R conflation + sum 1233.6 (L6334-6407).
- Also flagged H21cm 0.998366 (=CODE-34 D3).

### CODE-36 re-dispatch (caught CO-36a as N1, CO-36b as M1) -> P36 code recalibrated 2/2
- M2 SMALL: _pct helper sign-blind (owed_guards:50-53) though "Signed"; backstopped.
- M3 SMALL: G_SE,7 2P term (0.12 kHz) unguarded; budget docstring overpromises.
- NITs N2-N6 (redundant/identity asserts; N4 sympy identity not code).

### CODE-27 (CO-27a caught as D4; CO-27b rows[-1] runner-up MISSED -> PARTIAL; also caught CL-D as D5) -> re-dispatch assertion enumeration
- D1 OWNER-LARGE: cusp control edge (1s,2p0)<->(1s,3p0) is parity-forbidden (odd sector, amplitudes 1e-17) -> guard locus:111-119 cannot fail; in-sector edge (1s2s)<->(1s3s) |V| 0.49x moves S by -14% -> "mass at cusp vertex, not anywhere" false. Abstract L44-45, L539-545.
- D2 OWNER-LARGE: gamma_inf: sub-2 gap is the von Neumann log (S ~ eps ln 1/eps, eps exponent -> 2), present at n_max=2, not multi-shell aggregation; Z->inf limit is 2 (slope minimum ~Z 200, climbs back); conclusion "as Z, n_max -> inf" wrong. Abstract L60-68, L788-792, L1172-1175.
- D3 SMALL NO-TEST physical gamma_inf; matrix row L242 stale.
- D6 SMALL: sec:scope Li/Be numbers hop-model and gauge-dependent (converges with CLAIMS-27 L3); Be 1.23 guard pins LAPACK gauge.
- D7 SMALL NO-TEST list. NITs: commutator normalization header; ep2i docstring; papers/core path.
- U: pred2 "asymptotic to 2nd-order RS" earned.

### CITE-34 (CI-A1 Wineland vol caught; CI-A2 Hellwig digits caught) 2/2
OWNER-LARGE (C-REF):
- C1 Lyman-alpha NIST f-values printed 0.13881/0.27764/0.41645 vs NIST 0.13880/0.27760/0.41641 -> residual -41 ppm not +64; "+64 ppm Dirac+radiative" component built on wrong digits (L7085-7086, 7144-7146).
- C2 mud polarizability: CGV 2014 = 2.01(74) meV (not 1.690); 1.690(20) is Hernandez 2014; 1.941(19) Friar 2013 -> 1.697 after Coulomb distortion; "13% spread" exposure verdict falls (L5068-5069, 5147-5155, 7867-7914).
- C3 He 2^3P exposure: PY 2009 quote "eliminated this disagreement by reevaluating relativistic Bethe log" NOT in source; discrepancy is theory-vs-experiment (L7797-7849, 5426-5429).
- C4 He+ Lamb sourcing: 14041.18(13) credited to Pachucki PRA 63 042503 (a hydrogen paper); Layer-2 rows credited to Drake-Yan PRA 46 2378 (helium Rydberg states) (L5177-5180, 5207-5209, 5284-5290, 7798).
MATERIAL: muH 202.3706 credited to CREMA 2010 (is Antognini 2013) L4917, 7647; He-3 2^3S1 credited to Schuessler/Prior-Wang (is Rosner-Pipkin 1970) L5892, 5935; tritium "+12+-5 ppm Karshenboim budget" not in source + "Greene 2017" phantom (L5740-5834, 5705, 5739, 3802); Antognini -3.84 meV not in source (L7702-7720) -> "-4.3%" spurious; Li-7 digits inconsistent (L6422 vs 7930; source 803504086.6(1.0)); HD Komasa-Pachucki 2020 JCP wrong (PRA 102 012814; 224.51 not there) L6766-6862; Mu 1S-2S Mu-MASS Crivelli 2018 is proposal (Meyer 2000) L3772.
SMALL: drake_swainson title; strohmaier title; Fee 1993 authors; Bonilla Tab IV -> Table 3; Nieuviarts 2025; citation_status thirteen/14 & certifies Wineland.
UNVERIFIABLE: Eides Table XIX items (ar5iv incomplete) -- PM verified earlier from PDF text; Diepold; muHe 1378.521; etc.

### CODE-27 re-dispatch (caught CO-27b as F2) -> P27 code recalibrated
- F1 = control edge false positive (CONVERGED independently with CODE-27 D1): control nodes outside GS connected component; in-sector edge -14.0% (n3), -12.5% (n4).
- F3 OWNER-LARGE: P27 Sec V physical commutator 10.6%/9.3% and block fractions 0.945/0.926 (tex 494-501) untested; tests pin hop model k=1; verification list L1228-1235 claims "Verified".
- F4 SMALL: "saturates not vanishes" test can't reject decay (two points); paper L501 stronger than backing.

### COMPLETENESS-CRITIC
- (a) Unmeasured: FREEZE C8 P27/P26 headline lists stale (gamma 1.93-1.96 vs 1.96-1.97; edge share 42.5/40.0 vs 58.5/54.1; rho 0.94 vs 0.87; S_HO values; P26 0.3033/49.6 vs 0.330/~40); criterion 4.4 audits a retracted claim; C3 inline tiers near-absent in P34/P35/P36 (unmeasured); P36 C-REF measured by CODE-36 (REVIEW.TEX fetched).
- b5 OWNER-LARGE candidate: P34 III.28 (L2836-2847) apparatus-identity / PSLQ entropies are hop-model values (0.040811 vs physical 0.0377).
- b6 SURFACE-LARGE: P34 L9526-9527 "PSLQ negative proved ... categorically distinct" vs L2851.
- b7 STALE-LARGE: P35 L822-823 depth analysis live.
- b9 SURFACE: claim_test_matrix rows 233, 236, 241, 242, 245 stale (BACKED-SOUND on retracted values).
- SMALL: P35 L1425/1432 bare-graph ring; P26 title/INDEX row 133 decoupling headline; batch2 III26 test docstring false identity; P27 L842.
- Not closed: P34 He-3 Zemach input under review; batch1/3 bodies; INDEX outside group6.

## PM verification (2026-11-08)
- Verified against primary text in the real corpus: P35 graduation (L796-802, 1444), two-loop ready (L1343); P34 alpha graph-derived (L3075), PK inverted (L1186-1192 vs P17 L1172), flow tier (L680, 3152), He hop zombie (L7308, 7335, 7343), Lyman-alpha (L7085: NIST A-values imply f 0.13881/0.27760/0.41641; framework+mu 0.416423 is +35 ppm ABOVE NIST, residual sign flips), mud pol (CGV 2014 abstract 2.01(74) meV), PY quote absent (ar5iv 0905.3046: nu01 3 sigma discrepancy with experiment), PRA 63 042503 = hydrogen two-loop paper (Crossref), P27 gamma limit text (L1170-1175, L786-793), P26 O/C 9-fold (dense eigh: C 9, N 4, O 9, F 6).
- Converged independently (two reviewers): P26 uniqueness false; P27 control edge parity false positive; P27 sec:scope hop values; P34 pred:l2_presence; -(11/48) attribution; P35 graduation; P35 iff false.

## Calibration
21 seeds: first pass 18/21 (missed CL-A, CO-27b, CO-36a); all three recovered by focused re-dispatch (claims P34A status/scope; code P27 and P36 assertion enumeration). Controls 6/6 clean.
