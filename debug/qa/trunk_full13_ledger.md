# FULL #13 ledger (HEAD ac6a79c)
Gates: 16/16 PASS (trunk).

## cite: synthesis + P0 + P1 — 0 LARGE
- AUDIT-SMALL P1 L70 chung1997 attributes combinatorial D-A / kinetic-energy reading to Chung (normalized Laplacian book). -> SMALL sweep
- NIT orphans: synthesis loutey_paper23 (L1904), P0 paper18 (L1168)
- UNVERIFIABLE (not defects): synth kac_infinite1990 locus; bander_itzykson1966 "explicit" scale invariance
- NIT: synth L1461 mixed-elliptic target could cite Hain–Matsumoto 1512.03975

## claims P38 — 0 LARGE (reviewer); PM verified loci
- SMALL edit-bar history/dates: L126-127, L311 (rendered retraction token!), L521, L1828-31, L1854-55, L1270-82, L417/784/794/812
- SMALL insider vocab: L1540-41/L1559 "WH1 keystone"; L1051 "not a decoy"; L1744 "Sprint H1"
- SMALL C3 untagged lemmas continuum_lip L366, band_injectivity L389, kernel_condition L428, lifted_state L449
- BORDERLINE (SMALL vs LARGE) C3 L303-304 "a tapered amplitude makes O(1/n)" untagged; siblings L273/531/1601 tagged MEASURED. PM verified text. Fix: add [MEASURED].
- SMALL (AUDIT unit) L1' L948-964 PANEL-VERIFIED "every cross-shell pair" vs n_max 2,3 evidence; footnote L983-4 robustness clause unbacked; not load-bearing.
- SMALL stale "Why SU(2) harder" L1663-94 complications not used by proof.
- NITs: "rate" wording for Fejér moment; "five lemmas" vs "four"; L796 0-indexing; thm:main proof L3 mention; chi notation; L493 q.
- Upgrades: any complex central amplitude; L1' n_max=2 backed by test_p45_kplus_degeneracy.

## claims P32-B — reviewer 0 LARGE; PM re-grade:
- **LARGE (PM, doubt rule) OWNER** M4: cor:multi_loop_renormalization_wall L4905-4906 "LS-8a sprint empirically confirmed this structural verdict" (inside CONDITIONAL corollary) + rem L4942-4944 "Paper 57 ... record this theorem as the renormalization sub-sector closure" — contradicts L5090/L5108 "classified conditionally" and P57 L613 ("conditional structural reading"). Same family as FULL #12 OWNER-LARGE (single-cutoff re-tier). Also thm env renders "Theorem ... [OBSERVATION]". Fix: "is consistent with this reading"; "closure" -> "conditional classification". Gate: C16 entry for "empirically confirmed ... structural verdict"/"sub-sector closure".
- SMALL M1 G4 header "G4a reachable"/"Predicted outcome" L3921-3943 (gate g4a-live-path misses "reachable").
- SMALL M2 L3947 "already verified by H1" (H1 gauge content exhibited; P25/30 identification Observation)
- SMALL M3 L4197-4203 Door 4b "exactly two" needs per-factor-unimodular convention stated
- SMALL M5 L3274-3285 "adds to it ... sharpens, rather than adds to" self-contradiction
- SMALL M6 L4015-4018 "modular-Hamiltonian Pythagorean orthogonality" -> posited wedge-generator
- NITs: dated test counts; untiered novelty L3305; L3458 G2 electroweak; "bit-exact" L4252; CLAUDE.md §1.5 cited in paper L4185/L4605; L4582 "confirmation"; hekkelman year L4002/5442; v3.103.0 L5299; debug/ refs; history-note loci list (L3278..L5148)

## synthesis C9 — 0 LARGE
- SMALL S1 L1169-1176 M1 "collapses to ... 4/π" vs P18 (M1 = Vol(S^d); 4/π shares value, reading). PM verified: SMALL (same bullet says reading; not rate).
- SMALL S2 L1711-1720 conclusion "each ... frontier-of-field" vs P32 L5385 (3 of 5 walls internal)
- SMALL S3 L1652/1675 A/D partition untagged (OBSERVATION)
- SMALL S4 L1263 27.90×Q unscoped (main-group), no \gvq
- SMALL S5 history prose: abstract L67-70, L108-110, L149-150, L234-236, L481-483, L515-516, L799-801, L944-947, L1246
- NITs: L141 "directed"; Layer 4 L978-983; Layer 5 L995 SYMBOLIC tag scope; abstract L76 PANEL "every cutoff"; L1670 "is"; L1579; L1146; L649-656 withdrawn s/p listed as contribution

## claims P0/1/7 — 0 LARGE, 0 AUDIT-LARGE
- SMALL S1 P1 L378 "nonzero Berry phase requires complex amplitudes" (false: Z2 flux) AUDIT unit
- SMALL S2 P7 L825 Limitations: relativistic corrections understated; Paper 38 called spectral-triple paper (it's P32); P36 has Dirac-on-S3 bound-state corrections
- SMALL S3 history prose: P0 L732, L46; P1 L14, L262-267, L283-4, L390-1; P7 L17, L113, L693-700, L829, L833
- SMALL S4 P7 L735 eigenvalues "count modes" (degeneracy does)
- NITs: P7 L825 "extends"/"two rates proved"; P7 L686 analogy; P1 L376 list; P1 L403 k≈2; P1 L444; P7 L1050 "Sprint 3D"
- Upgrades: P1 ||L+||/||T+||=2 provable; P1 O(1/n) rate sketch closable; P7 L85 nnz≤4V linear

## claims P32-A
- **SURFACE-LARGE (PM verified)** M1 P32 conclusion L5705-5707 "convergence of that commutative factor is proved (explicit asymptotic rate [INTERNAL THEOREM + MEASURED])" — retired 4/pi-as-rate family (FULL #12); owner P38 + P32 L2830-2837 say bound d<=gamma_n, rate O(1/n) MEASURED. Fix: "is proved, with explicit bound d_GH <= gamma_n (Paper 38)". Gate: C16 "explicit asymptotic rate" near convergence/proved.
- SMALL M2 thm:axiom_complete L1504-1510 W_d paired with unit sphere (dual-Coxeter vs unit homothety)
- SMALL M3 rem:gh_status_caveat L2858 "Four gaps ... are named" vs relocated assembly (only iii stands)
- SMALL M4 rem:wh1_proven L2945 order-zero/one called open/downstream (they fail, measured)
- SMALL M5 L2495 TS-E3 "showing" M1-M3 sub-cases (reading)
- SMALL M6 L2505 "K's chain" -> ingredients' chains
- SMALL M7 rendered retraction tokens L1717, L2083 + history notes list (L922-943, L2029, L2256/2267/2338 inside theorem, ... L3089)
- NITs: "three structurally distinct sectors" L2047/L5722; table "≥2" vs measured range; L2776; L1772; P39 untiered L1707; L1848 uncited attribution; pi-source counterexamples vs (ii); L4835 theorem env; L5724 obs:pi_coulomb ref
- Upgrades: F,G orthogonally similar; order-zero diagonal row Holds; (ii)=>(iii) case-check

## cite P32 — 0 LARGE
- SMALL D1 L1023-1024 GH continuum clause cited to connes_vs2021/hekkelman_mcdonald (CvS defer GH; P32 L786 says so) -> cite vansuijlekom2021_fejer (orphaned)
- SMALL D2 L2874 bare arXiv:2005.08544; L2875 Latrémolière uncited (latremoliere2021dual orphan); L2870 Leimbach–vS no bibitem
- SMALL D3 L5285 "established in Krajewski 1998" overcredit
- SMALL orphans: 24 external + 2 internal bibitems uncited after relocation
- NITs L5769 title; L590 TeX line numbers

## cite P38+P7 — CLEAN
- NITs toyota "compact Lie groups" -> non-abelian; key/year mismatches; P38 L1388/L1688 "no integrable almost-complex structure" (S3 odd-dim: none at all) -> claims NIT; orphan paper7 in P38

## code P0 — 0 LARGE; 154 passed (--slow), exit 0
- SMALL (PM verified, borderline) P0 abstract L23 + L100-106 SYMBOLIC PROOF "(l,m) emerge ... central, mathematically rigorous" vs proof = count identity; L345 caveat "bijection is a labelling choice" not carried to abstract. Graded SMALL: body caveat present, label set = counted positions modulo naming.
- SMALL test_trunk_qa_kappa.py T4 :158-194 decorative asserts (hard-coded 8, unfalsifiable counterfactual)
- SMALL docs/claim_test_matrix.md L152 "force ... from packing alone" vs P0 "motivated, not forced"; L155 kappa counterfactual overclaim
- SMALL geovac/lattice.py L26-28 "Edges represent allowed transitions"; L217 "Requires separate V=-Z/r" stale
- SMALL coverage gap: algebraic-not-rational spectrum has no test (verified by reviewer via charpoly)
- Fire tests: 4 plants rejected. Re-derivation all match.
- Upgrade: one-sided bound provable; 8-λmax = C_sat/(n+1/2)^2 + O(n^-4)

## code synthesis — 114+32 passed, exit 0
- **SURFACE-LARGE (PM verified)** F2 tests/test_rabi_oscillation.py L7-8 "Proves that the GeoVac lattice produces quantitatively correct coherent quantum dynamics", L74 "Build the hydrogen system", L588 "It is a Time Machine." — runs hamiltonian="graph" (non-physical excited states). v5.17.1 family (graph credited w/ hydrogen). Fix docstring; C16 entry.
- **SURFACE-LARGE (PM, doubt rule)** F5 geovac/lattice.py L454-475 stitch_lattices "n_bridges ... H2 optimal 8-24", "Default: 16 (experimentally validated for H₂)" — superseded LCAO graph-concatenation (FCI-M guardrail). Fix: mark superseded, drop "validated". C16 entry.
- SMALL F1 synth L829-831 "Bloch–Siegert-corrected two-level period to 0.41%": coded BS formula (Ω_R/2ω) vs standard (Ω_R/4ω) (convention-sensitive; BS term 0.05%, ~8x below residual); guard misses a 0.6% Ω error (fire-tested). Claim-matrix row 156 BACKED-SOUND overstated. Integrator itself matches exact propagation ~5e-5.
- SMALL F3 geovac/dynamics.py L197-198 "exact ... hydrogen radial integrals"; coded constants wrong (1s-2p 0.7449 vs 1.2903 etc.)
- SMALL F4 unitarity test starts from eigenvector (unrestricted random state also 1.4e-13: not FP)
- SMALL F6 synth L513 cite test_paper1_block_spectrum for C_sat; L487 "same test" Casimir vs Δ_S3; lattice.py L69 node_weights "(potential energy)"
- SMALL 18-check eigenvalue coverage n=2,3 l=0 (already caveated); fock harmonics ↔ momentum wavefns checked for 4 states only

## code P7 — 0 LARGE; fast 513 pass/11 skip, slow 5+4 pass
- diff 85bb318..HEAD on test_fock_laplacian / fock_kernel: docstrings only, nothing weakened; 18 checks still 18
- SMALL M1 P7 "Hierarchy of Spheres": SO(6) free eigenvalues "reproduced ... to 1e-10 (test_algebraic_angular)" — tautological (casimir typed into diagonal, algebraic_angular.py L192/L312). Claim true (prior art).
- SMALL M2 test_paper7_gegenbauer_coupling docstring + matrix row 47 "SYMBOLIC + MEASURED" (no symbolic step)
- SMALL M3 P7 4.2x/6.6x vs test >3x/>5x
- gaps G1 chordal-ansatz table; G2 Li 2S; G3 form factors; G4 FCI basis invariance
- Upgrades: kernel eigenvalue all n via Funk-Hecke (verified k<=6); exact c^2 values
- Fire tests: 6 plants, all fired
- NOTE: concurrent agents clobbered shared scratchpad out.txt; serial reruns are the record

## code P38 — 338 pass + 9 OOM under xdist (17 pass serially, exit 0)
- **OWNER-LARGE (PM verified, backing-test validity)** tests/test_p38_tapered_amplitude.py:68-76 test_tapered_moment_is_order_one_over_n admits a c·log n/n rate with c<~0.4 (only rejects (4/pi)log n; vals[100]-vals[25]<0.6 ⇒ c·ln4<0.6 ⇒ c<0.43). Fire-tested: 0.75·taper+0.25 (genuine log rate c≈0.28) passes all claim-level asserts. Claim (P38 abstract, MEASURED O(1/n)) true per reviewer's exact sum rule (doubling estimator 0.113,0.059,0.030,0.0155,0.0078 → 0) but unbacked. Fix: assert doubling estimator n·γ(2n)-n·γ(n) halves/→0 out to n=1600 via exact sum rule; pin 50/200/400 values; fire-test the eps=0.25 plant. (Guard-writing separate pass.)
- SMALL lifted_state (a)/(b) spinor-window freeze missing (reviewer verified leak ≤2e-16; scratch incl.py)
- SMALL prop:kernel_condition joint injectivity needs one line (SO(4)-inequivalent bands)
- SMALL L3 cites test_r25_l3_lipschitz_bound (doubled model) not test_paper40_l3_compression
- SMALL b≈4.109 unpinned (reviewer fit 4.1093)
- SMALL infra: 6 memory-heavy tests not xdist-safe
- Re-derived: γ1=π, γ2=2.074551, Fejér doubling → 4/π from above; taper 6.3526..6.8644 match; conventions match; App A Step 3 OK. Fejér guard fires on 2/π and 2√2/π plants.

## code P1 — 0 LARGE; 216 pass (MemoryErrors under load pass serially)
- SMALL S1 O(V) appendix describes shift-invert eigsh(σ=-1, k=20); test times eigsh k=1 SA
- SMALL S2 log-holonomy test rebuilds weights locally; production plant asleep (P10)
- SMALL S3 Casimir n²−1 inside abstract [SYMBOLIC PROOF] block; backing float n_max≤8 + Pauli prior art
- SMALL S4 generic O(1/n) rate untested (self-disclosed)
- NITs N1-N3; Upgrades U1 Frobenius closed form, U2 nnz<4V theorem, U3 gearing=2 symbolic
- Fire: 11/12 fire; all values reproduce

## code P32 — 706 pass + 8 env-OOM (pass serially), 7/7 slow pass; FULL#12 OOM tests now 0.3-0.5s
- **SURFACE-LARGE (PM verified, doubt rule; trunk.done P32 LARGE list "finite four-witness literal BW identification")** tests/test_modular_hamiltonian_lorentzian.py L563 docstring "BW canonical strong identification" + asserts verdict=="STRONG_IDENTIFICATION_LORENTZIAN" (L198, L568) as the claim; module L9 says label is legacy/closure-only. Also tests/test_modular_hamiltonian.py L229 "Expected: STRONG_IDENTIFICATION". Out of trunk code units (Paper 42 tests), now fast-run. Fix: assert residuals, reword docstrings; C16 widen finite-four-witness-literal-identification to tests/ docstrings.
- SMALL S1 test_real_structure.py:246 print-only (item-4 offdiag residuals unpinned)
- SMALL S3 thm:forced_count proof "at n_max=2" spurious
- SMALL S4 abstract L120-122 "doubled" for Weyl model (undoubled Weyl sector is odd)
- SMALL S5 thm:single_cutoff in theorem env with [OBSERVATION] (also found by P32-A/B)
- Re-derived 8 constants all match; fire 6 plants rejected
- Upgrades: U1 order-zero fails for every J (commutant scalar); U2 closed forms 1-n/(n+1)^2

## completeness critic
- STAMP LEAK (PM verified): P32 ::Conclusion unit stamped current (last changed v5.17.8) carries the SURFACE-LARGE at L5705. Record; double audit fraction next run.
- Unexercised: C23 inverse-citation (candidates: C_sat closed form; gearing 2√6/3; taper O(1/n) (Jackson/Fejér–Riesz); c²(n,l)); no dedicated arithmetic-audit pass (in-panel arithmetic only); C9-currency mechanical check not recorded.
- Gate scope holes: C16 doesn't scan demo/, benchmarks/scripts/, non-trunk tests/, docs/*.md; C25 front door excludes README.
- Siblings PM-VERIFIED SURFACE-LARGE:
  * tests/test_universal_constant_origin.py L360-363 graph ratios "converge to 1:4:9:16 (Rydberg) ... topological property of graph's continuum limit (unit S3)" (+L2 "Universal Constant", L371 "derivation of -1/16")
  * papers/group3_foundations/paper_55_periods_of_geovac.tex L411 "4/π (asymptotic state-space GH rate)", L1629 "universal state-space GH rate", L1631
  * papers/group4_quantum_computing/paper_20_resource_benchmarks.tex L1054 "with asymptotic rate 4/π"
  * tests/test_paper34_projection_spot_checks_batch1.py L83-87 "universal topological constant kappa ... map graph eigenvalues to Rydberg energies"; misquotes CLAUDE.md §8
- Siblings critic-reported, NOT PM-verified: demo/demo_h2.py (LCAO 0.43% "accuracy achieved", optimized R, universal constant); tests/advanced_benchmarks.py + bulk_physics_puzzles.py (α "first-principles", fitted n=5; §14 archive; PI item); claims_register row 8 backing attribution test_p38_action_seminorm vs true window; README L49 "worked out as a discrete almost-commutative spectral triple" (borderline); test_p38_metric_convention "rate constant" naming + one-sided assert (borderline)
- SMALL siblings: many 4/π "rate" docstrings (test_trunk_qa_fejer_4_over_pi L1, test_paper40_universal_rate L1, modular_hamiltonian L1524/1534, central_fejer_compact_temporal L57, ecosystem_export L238); demo_h2_dirac; benchmark.py Li; docs/topic_to_paper_lookup L5 "Universal constant"; INDEX L59; track_nk_sprint2 L242; P32 L3853/3879/3902/L2219; P7 L755 rendered token + INTERNAL THEOREM on s/p proxy; P32 L2083 token, L2112 withdrawn L0 row; P32 L5535/L254 "forced-graph claim"; P32 L1441 |D| uncited; README badge 5.19.0

---
# REMEDIATION STATUS (resume point, 2026-10-03)
Done (uncommitted at time of writing):
- OWNER: P32 single-cutoff -> observation env; corollary LS-8a "confirmed" + remark "closure" fixed. Taper guard rebuilt (exact sum rule + doubling halving; fire-tested vs 0.75 taper+0.25 plant: FIRED).
- SURFACE: P32 conclusion; P55 (table, canonical ex. 3, L1667, L2442); P20 L1054-1067; P41 L2430; P57 L619 + bibitem L793; test_rabi_oscillation; test_universal_constant_origin; test_paper34_projection_spot_checks{,_batch1}; test_modular_hamiltonian{,_lorentzian} docstrings; geovac/lattice.py (header, class doc, stitch_lattices SUPERSEDED, 3 SMALLs); benchmarks/scripts/h2_dissociation.py; docs/topic_to_paper_lookup.md; docs/claims_register.md row 8.
- Archived: demo/demo_h2*.py -> demo/_archive/superseded/ (+README); tests/advanced_benchmarks.py, bulk_physics_puzzles.py -> tests/_archive/dead_ends/.
- C16: new fejer-constant-as-gh-rate, multi-loop-corollary-confirmed-closure, lcao-bridge-graph-validated, graph-model-credited-with-hydrogen-dynamics; widened kappa-universal-topological-constant (+tests/demo/benchmarks/docs) and finite-four-witness-literal-identification (+tests). All discrimination-tested (HEAD fires, now silent).
- SMALL sweep list appended to docs/qa/trunk.debt.md.
- Self-check round 1: 4 SURFACE-LARGE (P20, P41, P57, test/bench docstrings) -> fixed. Round 2 dispatched.
Remaining: self-check to zero LARGE; all gates green; CHANGELOG v5.20.2 + trunk.done.md change-log line + CLAUDE.md §2 one-liner + version; cert_stamps --stamp --run "FULL #13" --unchanged-since ac6a79c; /checkpoint (PI).

Self-check rounds LARGE: 4, 4, 2, 1, 0 (round 5 ZERO LARGE). Remaining: final gates, stamps, /checkpoint (PI).
