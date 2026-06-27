# `/qa group2` — whole-group run #1 notes (2026-06-26)

**Verdict: FAIL** (calibrated panel + many verified MATERIAL defects, several structural). First QA of the quantum-chemistry branch; group2 has never been QA'd before.

## Setup
- DoD `docs/qa/group2.done.md` FROZEN (PI 2026-06-26). Scope: 9 papers (8,11,12,13,15,17,19,FCI-A,FCI-M) + new synthesis. `--gate group2`.
- Worktree `../geovac-qa-seed-group2` (detached HEAD d0c9fcf) + uncommitted cert-target files copied in (fixed paper_15, new synthesis). 10 seeds planted (key `debug/qa/group2_seed_key.json`). Worktree removed + leak-scan clean.

## Calibration scorecard — panel FULLY CALIBRATED
- **Sensitivity 10/10.** Every seed caught by its expected dimension: seed1 (claims-A C3), seed2 (claims-B C6), seed3 (claims-B C8b), seed4 (claims-C C8b + code-FCI-M), seed5 (cite-1 C4), seed6 (cite-2 C4), seed7 (code-P8 C2), seed8 (code-P17 C2), seed9 (code-P13+code-FCI-A C2), seed10 (synthesis C9).
- **Specificity clean.** Controls M1–M6 all treated SOUND; zero false positives.
- All gating dimensions exercised + calibrated → the FAIL is trustworthy.

## Deterministic layer
- C10 (compile): CLEAN content across all 10 docs **after** fixing a genuine pre-existing bug in paper_15 (frontmatter `\begin{document}` was AFTER `\title/\author` → revtex author-macro cascade; + a `dcolumn` math-mode break at the `$^\dagger$` cell). The only residual `!` lines are **6 missing generated-figure assets** (P13×1, P15×3, fci_atoms×2) → draft-mode fallback; figures carry no physics claim (all numbers in tables/text). Asset-hygiene item (the figure analog of the debug/-ref debt; fci_atoms figures exist at `debug/plots/` but are referenced as `figures/`).
- C11/C12/C13/C14/C15/C16: all PASS (`--gate group2`).

## "Code decides" (PI directive) — He numbers
code-P13 recomputed: 2D-var raw l_max=7 = **0.0217%**, cusp l_max=4 = **0.0035%**, graph-native n_max=5 = **0.2496%** (trend→0.19% at n_max=7). The code **confirms the paper** (0.022%/0.004%/0.19%) and the v4.49.x **§2 best-results fix**. CLAUDE.md **§5** + `docs/claims_register.md` still carry the stale 0.019%/0.20%(n_max=9) → §5 needs PI sign-off (not PM-editable), claims_register PM-fixable.

## Converged genuine MATERIAL findings (seeds set aside)

### TIER 1 — STRUCTURAL: headlines backed by deleted/retired code (PI strategic decision)
- **A. Paper 11 — entire spectral headline backed by code DELETED in v2.7.0** (commit 8d692a0). Spectral Laguerre solver (0.0002%, 250×/270×/5000×, σ-algebraic, π/δ single seed e^aE₁(a)) — live `ProlateSpheroidalLattice` has no `radial_method`/`n_basis`; the methods are gone; tests archived + ImportError. Paper cites the removed API. Live backing = FD solver at ~1% + qualitative CI only.
- **B. Paper 15 — spectral solver (16×/20×/269×) REMOVED v2.7.0**, presented as live (`angular_method='spectral'` → TypeError); the 96.0% D_e (l_max=6) and the 2D-variational D_e have NO working test (only-2D test ERRORs on stale kwargs); "exceeds Paper 12" tests pass only via the adiabatic solver the paper disavows (σ-only test asserts the OPPOSITE of the paper); "σ–π completely decoupled" contradicted by its own backing test + code.
- **C. FCI-atoms & FCI-molecules — backed by RETIRED `MolecularLatticeIndex`.** LiH D_e^CP=0.110 "first heteronuclear" + the entire FCI-M pipeline: sole test archived + import-broken. The keystone **graph-concatenation NEGATIVE** (monotonic PES, no equilibrium) has NO collected test (only docs + an ignored debug script). The Sturmian Structural Theorem (P8) likewise has no matrix-level test.

### TIER 2 — STALE / NON-REPRODUCIBLE headline numbers (verify + correct)
- **D. Paper 17 BeH₂ 11.7% → code gives 19.7%** (status-drift; test gate `<30%` masks it). l_max=2-optimal / structural-divergence = NO-TEST (`assert True`). H₂O 26% = NO accuracy test.
- **E. Paper 19 balanced LiH 0.20% energy NON-REPRODUCIBLE** — recompute E=−15.21 (paper −7.924), D_e flips sign (over-binds 0.158 vs paper under-binds 0.037). Stale post the 2026-06-07 V_NN double-count fix. NO test pins the accuracy.
- **F. Paper 13 H₂ rovibrational (ν₀₁=4157, 0.1%) NON-REPRODUCIBLE** — cited PES (`benchmarks/ab_initio_nuclear/`) recorded FAILED; working pipeline gives +10.5%. graph-native 0.19% @ n_max=7 = NO-TEST. The 0.004% cusp result is BELOW exact (−2.903826 < −2.903724) → NOT "properly variational" as the paper emphasizes.
- **G. Paper 12 "no numerical integration of any kind" is FALSE** — B_l computed by `scipy.integrate.quad` (paper candid in §III.D but headline/abstract/conclusion overclaim); the B_l backing test is tautological (quad vs identical quad); the 92.4% D_e headline + Table I have NO backing test; "12–20pp improvement" test mis-asserted (passes if Neumann is worse).

### TIER 3 — CLAIMS / §3-suppression (CONFIRMED vs primary text)
- **H. l_max-divergence §3-suppression + internal contradiction.** P17:1068–1072 ("divergence entirely from adiabatic, not PK; 2D is the correct fix path") **contradicts** P17:1173–1175 ("intrinsic to the PK architecture, not adiabatic; l_max=2 optimal") — the §3-correct version. P15:1031–1035 propagates "2D is the path to eliminating the l_max divergence in composed geometries" — contradicts the §3 record (2D gives IDENTICAL drift; l_max-via-2D is a documented dead-end). Reconcile corpus-wide to the structural-ceiling framing.
- **I.** fci_atoms LiH "minimum near R=2.5, contracted 17%, traces to basis incompleteness" (475–477) re-asserts a spurious equilibrium + a FALSIFIED mechanism, contradicting the corrected companion FCI-M (LARGE). P11 "zero free parameters" on the Eckart-optimized-exponent H₂ (SMALL). P13 Li quasi-Coulomb "<0.01%" is a tautological back-fit presented as predictive accuracy + a trend (SMALL).

### TIER 4 — Citations (genuine, non-seed)
- **J.** `shull_lowdin` spliced wrong vol/page/year (P8 + FCI-M); `Avery2004` spliced title (P8); `bishop1979` (HeH⁺ paper) cited as source of BeH⁺ reference numbers (P17, wrong-attribution); `koga1987` common-vs-different-exponent misattribution (P19); `hehre1969` dangling `\cite` with no `\bibitem` (FCI-A). + NITs (Aquilanti page/initials, herbst orphan, hoggan chapter title).

### Deterministic-fixed this run
- paper_15 frontmatter `\begin{document}` reorder + `\dagger` dcolumn fix (committed-pending). 6 missing figures flagged (asset hygiene).

## Disposition
Far beyond "PM fixes small issues directly." TIER 1 (A/B/C) requires a PI strategic decision: **restore the deleted spectral / MolecularLatticeIndex code, or honestly descope/reframe the affected paper headlines** (the v2.7.0 removal + the MolecularLatticeIndex retirement happened after these papers were written; the papers were never reconciled). TIER 2–4 are substantial but mechanical-ish once TIER 1 strategy is set. Recommend PI direction before remediation.

---

## Remediation (2026-06-27, PI-directed: investigate TIER-1 first + start independent fixes)

### Independent mechanical fixes — APPLIED (all papers recompile content-clean)
- **TIER-4 citations (6):** P8+FCI-M shull_lowdin 25,1035(1956)→30,617(1959); P8 Avery2004 retitled to match its 108,8848 coords; P17 bishop1979→huber1979 for the BeH⁺ values (HeH⁺ paper was wrong) + unused bibitem removed; P19 koga1987 common-vs-different-exponent prose; FCI-A missing hehre1969 bibitem added. (P13 Pekeris + P19 Shibuya were seeds — real corpus correct.)
- **TIER-3 §3-suppression:** P17:1068 + P15:1031 reconciled to the structural-ceiling framing (fixed P17's 1068-vs-1185 self-contradiction; matches §3 record).
- **TIER-3 P12 "no integration of any kind":** abstract/intro/conclusion reframed to "single 1D quadrature for the log-singular B_l" (matches §III.D).
- **claims_register** He row + **P11/P13 wording** (zero-free-params on Eckart H2; Li quasi-Coulomb tautology) corrected.

### TIER-1 restores (PI: full restore plan) — code restored + VERIFIED GREEN
- **P11 spectral-radial → production:** prolate_spheroidal_lattice.py 364→982 lines (pre-removal 8d692a0^ restored; FD path preserved). Tests green: test_prolate_h2plus (FD gate), test_paper11_spectral_radial, test_paper11_associated_laguerre, test_paper11_associated_kinetic. H2+ 0.0002% re-grounded.
- **P15 spectral-angular:** un-archived test_paper15_spectral_angular.py (module was always live) — GREEN.
- **P15 spectral-radial (16×) → production:** level4_multichannel.py 1406→2033 lines; hyperspherical_radial.py rewired; test_paper15_spectral_radial.py GREEN.
- **MolecularLatticeIndex → geovac/_archive/superseded/molecular_lattice_index.py** (4916 lines, NOT production — guardrail preserved) + un-archived test_fcim_lih_molecular.py. [verification pending]
- **18 symbolic S³ proofs (test_fock_*) PASS** after all restores. level4/hyperspherical consumer regression + FCI-M in progress.

### Still pending
- TIER-2 stale-number recompute (P17 BeH2 11.7%-vs-19.7%-vs-28%, P19 balanced-LiH 0.20%, P13 H2 rovib) — each needs a verification recompute before editing.
- §5 He-number sync (line 351) — exact diff prepared, awaiting PI sign-off (§5 not PM-editable).
- Re-run /qa group2 to convert FAIL→remediated into a certified PASS (after all tiers land).

### TIER-1 restore — FINALIZED + VERIFIED (2026-06-27)
- **P11 spectral-radial** (prolate_spheroidal_lattice 364→982, FD byte-identical): test_paper11_spectral_radial 36 + associated 119 + test_prolate_h2plus 11 PASS; H2+ spectral 0.00000% (n_basis=20). PM-verified independently.
- **P15 angular**: test_paper15_spectral_angular 18 PASS (module was always live; archive note was false). **P15 16× radial**: hyperspherical_radial +~328 (genuine spectral helper resurrected — NOT rewired to FD, which would have been tautological; agent's correct deviation from my prompt), level4_multichannel +~450; test_paper15_spectral_radial 13 PASS; 16× reproduced, spectral↔FD <0.001 Ha.
- **MolecularLatticeIndex → geovac/_archive/superseded/** (4915 lines, NOT production): test_fcim_lih_molecular 47 PASS --slow / skip-by-default; cross-atom J within 1% of reference (bit-faithful); LiH bound + monotonic-PES + BSSE<0 asserted.
- **Regression**: 18 S³ proofs PASS; level4_multichannel 45 PASS (FD clean); hyperspherical consumers 105 PASS; PM's own consumer regression 133 PASS. No speed regressions (additive spectral behind non-default flags).

### Restore flags (carry to re-cert / claim_test_matrix)
1. **FCI-M magnitude coverage gap**: paper D_e^CP=0.110/BSSE=0.115 is **nmax=3** (367,290 dets = C(56,4)); the un-archived test asserts only *signs* at nmax=2 (gives 0.2875/−0.2184). Exact-determinant-count match confirms same code; the nmax=3 magnitude has no value-asserting test (367k-SD run too slow for CI). FCI-M is a negative-result paper — the no-equilibrium/monotonic headline IS tested; only the 0.110/0.115 magnitude is gap.
2. Restoring solve_radial_spectral re-enables non-default spectral paths in n_electron_solver.py + algebraic_coupled_channel.py (were ImportError/NotImplementedError stubs); FD defaults unaffected (additive). Spectral paths there untested.
3. sw_form_factor not restored (used only in MolecularLatticeIndex use_dmatrix/use_sturmian paths, not reached by default FCI).
4. Generated artifacts (not committed): geovac/_archive/superseded/cache/, tests/_durations.json.

### TIER-2 recompute results (2026-06-27)
- **P19 balanced-LiH (nmax=2, correct V_NN convention)** [PM driver, fine R-grid]: R_eq=**3.20** bohr (6.1%; paper 3.226/7.0% — consistent, fine), D_e=**0.161 Ha** (well depth to R=5.0; matches CHANGELOG/test ~0.158). Absolute E_coupled=−15.21 is an offset internal quantity (NOT the −8.07-comparable total — the ~7.3 Ha core/convention offset the reviewer saw), so D_e (difference) + R_eq (argmin) are the decisive convention-robust numbers. **MATERIAL: paper D_e=0.037 (under-binding) is stale → correct ≈0.16 (OVER-binding vs exact 0.092); binding direction flips.** R_eq fine. Single-point energy claims (1.8%/0.20%) are offset-sensitive accounting not reproduced by this driver — likely unaffected by the shape-corrector fix; flagged for the paper edit. → HOLD paper_19 edit for PI sign-off (material narrative change).
- **P17 BeH2 (full-1RDM-exchange, l_max=2)** [PM driver]: R_eq=3.00 bohr (grid 0.2), R_ref=2.507 → error **19.66%**. **MATERIAL: paper 11.7% → current 19.66%** (~2× worse; test-comment "~28%" also stale; the number "drifted as PK and cross-block h1 evolved" per the test comment). → HOLD paper_17 edit for PI sign-off.
- **Pattern**: group2 composed/balanced chemistry headline numbers are stale-optimistic vs the current (post-v2.7.0 / post-2026-06-07-V_NN-fix / post-PK-evolution) code. "Code decides" → update to current values (worse), honest framing, PI sign-off on the material ones (BeH2 11.7→19.66%; LiH D_e under→over-binding).
- **P13 H2-rovibrational** [agent a6fda diagnosis]: **NOT a defect — stale benchmark FILE only.** Hylleraas+Neumann solver runs clean today (E=−1.160955 Ha @R=1.4, reproduces archived data); paper §IX ω_e=4435/ν₀₁=4157/B_e=59.5 REPRODUCE to ~1% (residual = grid: paper used linspace(1.0,2.0,8); benchmark was 0.2-spaced). The "ALL FAILED" rows were junk written at v1.2.0 (d9568c2) overwriting good α-opt data from 27693f5; generating script pruned. Reviewer's +10.5% was the WRONG (prolate-CI) pipeline. **NO paper change.** FIX (PM, benchmarks/): restored h2_neumann_pes.txt from 27693f5 (un-FAILed; 0 FAILED rows). Optional polish: re-run on linspace(1.0,2.0,8) to land exactly on 4435 (deferred). No geovac/ edit (latent dead-helper arity bug in neumann_vee._compute_vee_numerical_p0 noted, not on Neumann path, untouched).

### TIER-2 CORRECTION after "investigate first" (2026-06-27)
- **P19 balanced-LiH — REPRODUCES; NO paper edit needed.** Under the paper's TC convention (E_min(TC)=E_coupled−E_core_off, E_core_off=_FIRST_ROW_CORE_ENERGY[3]=−7.2799): E_min(TC)=**−7.933 Ha (1.71%)** ≈ paper −7.924/1.8% ✓; R_eq=**3.230 (7.1%)** ≈ paper 3.226/7.0% ✓; D_e=0.057 (to R=4.2; paper 0.037 to farther ref — same small-under-binding regime, reference-dependent) ✓. **My + the code-19 reviewer's "non-reproducible E=−15.21 / D_e flip to over-binding 0.16" was an ACCESSOR ERROR** (raw E_coupled without the E_core convention; D_e to R=5.0 vs paper's nearer ref). The −15.21 is expected (E_core+V_NN baked into the Hamiltonian constant), NOT a regression. **Disposition: P19 numbers stand; residual = §9 backing gap (archived driver chemistry_solver_retest_lih.py, not a permanent test) → log coverage gap, backfill a permanent test replicating the TC-convention computation.** code-19 P19 "non-reproducible" finding = FALSE POSITIVE (wrong accessor).
- **P17 BeH2 — CODE REGRESSION (not stale paper).** Current code: monopole 4.00/59.6% (paper 4.06/62% — unchanged), exchange 3.00/19.7% (paper 3.01/20% — unchanged), **full_exchange 3.00/19.7% (paper 2.80/11.7% — collapsed EXACTLY onto the exchange value)**. The off-diagonal 1-RDM exchange contribution evaporated; specific to the full_exchange path (other modes unchanged) → regression introduced during PK/cross-block-h1 evolution. **Paper's 11.7% was correct; fix the CODE (inter_fiber_coupling/composed_triatomic full_exchange path), not the paper.** code-17 BeH2 "11.7%→stale" finding reframed: it's a regression to repair.
- **NET TIER-2: no group2 paper headline needs downgrading.** P13 (benchmark restored), P19 (reproduces, accessor-error false-positive), BeH2 (code regression to fix). The "papers stale-optimistic" hypothesis is overturned; the PI's "investigate first" prevented two wrong downgrades + surfaced one real code regression.

### BeH2 regression FIXED + VERIFIED (2026-06-27)
- **Root cause:** v2.7.0 composed_qubit refactor dropped `pk_potentials` from the `solve_angular_multichannel` call in `inter_fiber_coupling.extract_channel_data` → inter-fiber angular eigenvectors (vec_2d) PK-blind while radial F was PK-aware (v3.56.0) → off-diagonal 1-RDM exchange (full_exchange) collapsed onto the diagonal exchange minimum.
- **Fix:** restore `pk_potentials=pk_potentials` propagation (geovac/inter_fiber_coupling.py, +15/−7). Physically: vec_2d now PK-consistent with F. composed_qubit passes None → bit-unchanged.
- **Verified:** full_exchange R_eq 2.800/11.69% (paper 2.80/11.7% ✓), exchange 3.000/19.66% (unchanged ✓), monopole 4.00/59.6% (unchanged ✓); fast guard `test_extract_channel_data_propagates_pk` PASS; tightened slow `test_full_exchange_reduces_req_error` (<15% + full<exchange); 18 S³ proofs PASS; test_composed_h2o 17 PASS.
- **Consumer-safe:** LiH composed (ComposedDiatomicSolver) does NOT use extract_channel_data — unaffected. H2O (ComposedWaterSolver) DOES + passes non-None pk_potentials → bond_channel_data now PK-consistent (more correct); H2O 26% R_eq (untested coverage gap) may shift slightly toward correct — recompute in the H2O backfill. **Paper 17 BeH2 11.7% reproduces; no paper edit.**

### TIER-2 CLOSED — no group2 paper headline downgraded
P13 (benchmark restored, reproduces), P19 (reproduces under TC convention; accessor-error false-positive), BeH2 (code regression FIXED, 11.7% recovered). The "papers stale-optimistic" hypothesis fully overturned.

### Backfill phase (4 agents, 2026-06-27) — partial
- **P19 permanent test** `test_paper19_balanced_lih.py`: 5 PASS — E_min(TC)=−7.933 (1.71%), R_eq=3.223 (6.9%); §9 archived-driver backing gap CLOSED (recompute-from-framework TC-convention test).
- **P17 H2O — MATERIAL IMPROVEMENT (flag for PI):** post-BeH2-fix (PK-consistent extract_channel_data now feeds H2O bond_channel_data), H2O composed R_eq = **1.459 bohr (19.4%)** vs paper 1.34 bohr (26%). 19.4% is CLOSER to experimental 1.809 → the fix improved H2O. Paper 17 H2O 26% → 19.4% is a material paper update (an improvement, not a downgrade) requiring PI sign-off. Test updated to current value; CONFIRMED 2026-06-27: H2O R_eq=1.4589 bohr (19.4%), test passes.
- **P15** (Agent 3, DONE): false-positive tests fixed (2D solver: σ-only 87.3% below P12, σ+π 94.1% above, <100 variational discriminator); l_max=4 2D+cusp 94.3% test + honest NO-TEST note for l_max=6/96.0%. Flag: tests/test_cusp_correction_2d.py still broken (stale E_atoms kwarg) — separate cleanup.

### Backfill COMPLETE (4 clusters + P12 + cusp, 2026-06-27) — all verified
- **P8** `test_paper8_sturmian_structural.py` (7✓): Sturmian H∝S + R-independent eigenvalues (the headline NEGATIVE) now matrix-level tested.
- **P12** `test_neumann_vee.py` +7 (32✓ --slow): B_l tautology replaced by recurrence + B0 closed-form + Q1 anchors; V_ee machine-precision vs independent Legendre-reduced assembly (~2.8e-11) + literal-4D-quad normalization; 92.4%-D_e recompute (E=−1.160960, 92.11%, matches N=27).
- **P13** `test_direct_ci.py` (8✓): vacuous E<0 FALSE-POSITIVE tightened to pin −2.845; FCI-A He hybrid-config + graph-native (`test_casimir_ci.py` 95✓).
- **P15** `test_level4_multichannel.py`: "exceeds-P12" FALSE-POSITIVE fixed (2D solver: σ-only 87.3% below, σ+π 94.1% above + <100 variational discriminator); l_max=4 2D+cusp 94.3% + honest NO-TEST note for l_max=6/96.0%.
- **P17** `test_composed_diatomic.py` l_max-divergence (TestLmaxDivergenceMonotone, 2 PASS) + `test_composed_h2o.py` H2O R_eq=1.459/19.4% (CONFIRMED; paper updated 26%→19.4% + provenance, PI-signed-off).
- **P19** `test_paper19_balanced_lih.py` (5✓): permanent TC-convention test (E_min(TC)=−7.933/1.71%, R_eq=3.223/6.9%) — §9 archived-driver backing gap CLOSED.
- **cusp**: test_cusp_correction_l4 E_atoms-fixed (9✓); test_cusp_correction_2d ARCHIVED (stale E_atoms-reference API; 2D+cusp covered by P15 test_lmax4_2d_cusp_de).
- P17+synthesis recompile clean (0 content errors) after H2O edits. No geovac/ bugs found in backfill (all production verified correct). NEXT: re-run /qa group2 (re-cert).

---

## §6 Honest scope (sprint-close)

- **Nothing new at theorem grade.** This was a *remediation*, not new physics — the value is in restoring/verifying existing claims, not proving new ones.
- **Restorations (code re-grounded), VERIFIED:** P11 spectral-Laguerre (H2+ 0.00000% recompute, FD byte-identical); P15 spectral angular (18✓) + 16× radial (13✓, genuine spectral helper, not FD-masquerade); FCI MolecularLatticeIndex→_archive (cross-atom J within 1% bit-faithful; sign-level @n_max=2).
- **Regression fix, VERIFIED:** BeH₂ full_exchange off-diagonal restored (2.80/11.69% recovered; fast + slow guards green). Physically-principled (PK-consistency), not a fudge.
- **Improvement, VERIFIED:** H2O 26%→19.4% (config-faithful uncoupled recompute; closer to experiment; from the same PK-consistency fix).
- **Reproductions (numerical, code-confirmed = paper):** P19 balanced-LiH TC-convention E_min=−7.933/1.71%, R_eq 3.223/6.9%; P13 H2-rovib ω_e/ν₀₁ to ~1%; P12 92.4%-D_e = 92.11%; He 0.022/0.004% & 0.19%.
- **Honest caps / named open follow-ons:**
  - FCI-M LiH D_e^CP=0.110 / BSSE=0.115 is **n_max=3 (367k-SD), NO-TEST in CI** — too slow; sign-level @n_max=2 is the guard (determinant-count C(56,4)=367,290 confirms same code). Documented, not faked.
  - P15 96.0% D_e @ l_max=6 — **NO-TEST** (uncomputable in CI; l_max=4 2D+cusp 94.3% tested + in-paper note).
  - P13 graph-native 0.19% @ n_max=7 — **NO-TEST** (n_max=5=0.25% confirmed; n_max=7 impractical).
  - The PK-consistency fix changed PK-aware paths through `extract_channel_data`; verified BeH₂ + H2O; its other consumers (`algebraic_slater`, `lone_pair`) not re-checked for value shifts (FD defaults unaffected).
  - **Re-cert PENDING** — the certified PASS requires a fresh calibrated `/qa group2` run on the committed remediated corpus (this release lands the remediation; the re-cert is the next step).
