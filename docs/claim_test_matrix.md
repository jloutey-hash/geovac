# GeoVac Claim → Test Matrix

**Purpose:** the granular source-of-truth mapping each load-bearing *paper claim* to its *backing test* and *code module*, with an adversarial verdict on whether the test actually **proves** the claim. Companion to the reader-facing `docs/claims_register.md` (21 headline claims). This matrix is the per-paper, per-claim QA artifact of the **§9 Branch QA Review Protocol**.

**Populated by:** the per-paper `code-reviewer` agent (`.claude/agents/code-reviewer.md`) during each branch's QA cycle. Built in dependency order — **trunk first** (Papers 0, 1, 7, 32, 38), then branches.

**Classification:**

| Tag | Meaning |
|---|---|
| BACKED-SOUND | test exists, passes, and genuinely proves the claim |
| BACKED-WEAK | test passes but proves *less* than the prose (special case / necessary-not-sufficient) |
| NO-TEST | claim has no backing test (coverage gap) |
| FALSE-POSITIVE | test passes for the wrong reason |
| BUG | test or code is wrong |

**Disposition rule:** a load-bearing claim that is NO-TEST / BACKED-WEAK / FALSE-POSITIVE / BUG is **raised to the PI**. Small backings (minor claims, easy missing tests) the PM fixes directly. New equations follow the §13.4a naming convention (`test_paper{N}_*`).

---

## Trunk (Papers 0, 1, 7, 32, 38) — COMPLETE (code review + test-writing, 2026-06-14)

5-agent per-paper code review + 2-agent test-writing pass; PM-verified the two headline constants (κ, 4/π) against primary text. 5 new artifacts accepted, 5 claims downgraded; all `test_trunk_qa_*` green (48 passed, 2 slow-skip, 2 xfail).

**Inline-reference pass (2026-06-15):** the trunk papers now cite their backing tests *inline* (previously the mapping lived only in this matrix). Added: Paper 0 → `test_trunk_qa_annular.py` (angular capacity), `test_dirac_lattice.py` (|V|); Paper 1 → `test_ov_scaling_rigorous.py` (O(V)); Paper 7 → `test_fock_projection.py` / `test_fock_laplacian.py` (the 18-proof appendix modules), `test_trunk_qa_c2_delta.py` (c²(4,3)=1/40); Paper 38 → `test_central_fejer_su2.py` (γ_n rate) / `test_trunk_qa_fejer_4_over_pi.py` (4/π). Paper 32 already cited its tests inline (20 refs). All four edited papers compile ERRORS=0.

| Paper | Claim | Backing test | Code module | Status | Verdict |
|---|---|---|---|---|---|
| 0 | Annular area → 2ℓ+1 angular capacity | `test_trunk_qa_annular.py` | `lattice.py` | BACKED-SOUND | ✅ accepted (new) |
| 0 | Cumulative shell count Σ(2k−1)=n² | (exact, agent-verified) | `lattice.py` | BACKED-SOUND | sound |
| 0,7 | **κ = −1/16 "derived from Fock projection"** | `test_trunk_qa_kappa.py` | `lattice.py`, `atomic_solver.py` | **OBSERVATION** (was "SYMBOLIC PROOF") | ⬇ **downgraded** — matching scale *coincides* with geometric 1/16; counterfactual shows no bridge |
| 0 | \|V\| spin-doubled vertex set = Σ2n² = n(n+1)(2n+1)/3 | `test_dirac_lattice.py` (state-count, spinor-doubling, per-shell) | `dirac_lattice.py` (spin-doubled — NOT spinless; 1st-pass map was wrong) | BACKED-SOUND (code) | ⬇ paper formula BUG fixed 2026-06-14: §VI had 2n_max² (top shell only) → Σ2n²; code was always correct |
| 1 | O(V)/O(N) sparse eigen-scaling | `test_ov_scaling_rigorous.py` | `lattice.py` | BACKED-SOUND | sound (multi-point fit ~1.0) |
| 1,7 | Spectrum = −(n²−1) | `test_fock_laplacian.py` (continuum, n≤3) | symbolic S³ L–B | BACKED-SOUND (continuum only) | ⚠ the **discrete graph L=D−A is positive-semidefinite**; −(n²−1) is a *continuum* property, not the graph's |
| 1,7 | 2s/2p splitting → 0 (convergence) | `test_trunk_qa_splitting.py` | `lattice.py` | BACKED (oscillatory decay) | ✅ accepted (new) |
| 1 | §III degree table D_2s=0.854, D_2p=3.416 | `test_trunk_qa_degree_table.py` | — (CG construction absent) | DOWNGRADED | ⬇ faithful build gives 1.9/3.8; numbers not reproduced |
| 1 | §III 16% splitting @ n_max=10 | `test_trunk_qa_splitting.py` | `lattice.py` | DOWNGRADED | ⬇ real ~1.7% @ n=10 (~16% lands near n=8) |
| 1 | Gearing ratio ‖L₊‖/‖T₊‖ → 1.77 | `test_trunk_qa_gearing.py` | — | DOWNGRADED | ⬇ not convergent (spectral norm = 2.0) |
| 1 | [T₊, L₊] = 0 | `test_trunk_qa_gearing.py` | `lattice.py` | BACKED-SOUND | ✅ accepted (bit-exact) |
| 7 | 18 symbolic Fock-S³ proofs (continuum geometry + S³ eigenvalues n≤3) | `test_fock_projection.py`, `test_fock_laplacian.py` | symbolic | BACKED-SOUND | sound; ~4 of 18 weak/tautological (volume Jacobian, etc.) — soften "18 independent" |
| 7 | c²(4,3) = 1/40 = Δ (Paper 2 boundary term) | `test_trunk_qa_c2_delta.py` | `dirac_s3.py` | BACKED-SOUND | ✅ accepted — two independent routes to 1/40 (genuine bridge) |
| 7 | V_ee Slater integrals F⁰ on S³ (§V closed form eq:f0_s3) | `test_paper7_vee_s3.py` (NEW; exact ∫ → 5π/32, 17π/324, 77π/2048) | `algebraic_slater.py` | BACKED-SOUND | ✅ NEW — §V closed form now has a live symbolic test (direct test was archived w/ removed code in v2.7.0); `test_algebraic_exchange.py` tests a different (Laguerre-grid) method |
| 1 | §II Rydberg E_n=−1/2n² "from operator eigenvalues" (SU(2)⊗SU(1,1) algebra) | `test_paper1_rydberg.py` (NEW; N=−2[T₊,T₋]→{n}, L²→l(l+1), bit-exact) | algebra from documented matrix elements | BACKED-SOUND | ✅ NEW — n & l(l+1) emerge from the ladder algebra (non-tautological); was NO-TEST (continuum-only backing, PSD graph) |
| 32 | Connes axiom audit (J²=−1, JD=+DJ, order, Lorentzian signs) | `test_real_structure.py`, `test_connes_axiom_audit_31.py` | `real_structure.py` | BACKED-SOUND | sound — bit-exact + negative controls |
| 32 | K = π(B+F−Δ) | `test_paper2_corrections.py` (arithmetic only) | — | COMPLIANT (observation) | ✅ §13.5 honored — no test asserts derivation; conjecture→observation downgrade 2026-06-14 |
| 32 | Forced-Count moduli chain 1024→512→128→128 | `test_trunk_qa_forced_count_moduli.py` | `standard_model_triple.py` | DOWNGRADED | ⬇ full-axiom = 260, not 128; 128 = matter-sector; "45 tests verify" mis-pointed |
| 32 | §VIII non-selection theorems (N_gen, KO-dim, H1 Yukawa, single-cutoff, …) | none | structural prose | PROOF-BY-ARGUMENT | not test candidates; verify the argument (paper review) |
| 38 | "Unconditional" state-space GH convergence | `test_p38_action_seminorm.py` | `full_dirac_operator_system.py` | QUALIFIED-SOUND | honest to the translation-seminorm scope |
| 38 | Per-band injectivity (full N² rank) | `test_p38_action_seminorm.py`, `p38_g1g2_band_diagnostics.py` | `spinor_operator_system.py` | BACKED-SOUND (general leg = Schur irreducibility + AWA non-vanishing; finite n≤5 corroborates) | ⬆ upgraded from BACKED-WEAK (2026-06-14 re-check) — the proof is all-N, not finite-cutoff |
| 38 | **4/π asymptotic constant** | `test_trunk_qa_fejer_4_over_pi.py` | `central_fejer_su2.py` | BACKED (derived, **numerics-pinned**) | ✅ rejects 2/π decoy; subcoefficients pinned by doubling estimator (not executed symbolic E–M, so not full SYMBOLIC PROOF); **supersedes** circular `test_asymptotic_constant_value` (still live, now annotated as a constant lock) |
| 38 | γ_n closed-form sum rule | `test_central_fejer_su2.py` | `central_fejer_su2.py` | BACKED-SOUND | sound |

**Pattern:** the trunk's machinery (axiom audit, O(V), continuum symbolic proofs, sum rules, K-conjectural) is genuinely backed; the soft spot was specifically the *headline derived/exact constants* — κ (coincidence, downgraded) and 4/π (genuinely derivable, upgraded from a hardcoded circular test). Plus Paper 1 §III quantitative content (from a CG construction never in production).

## Foundations branch (Papers 18, 22, 24, 31, 54–57) — IN PROGRESS (group3 first bite: 22/24/31, 2026-06-16)

group3 `/qa` first-pass test-backing (3 coder agents, 2026-06-16) closed the three NO-TEST C1 gaps run #3 surfaced. New tests all green (28 passed, 7 slow-skip).

| Paper | Claim | Backing test | Code module | Status | Verdict |
|---|---|---|---|---|---|
| 22 | Thm 3 headline global-M_L Coulomb density D = 14.84/8.52/6.06/4.83/3.99% (l_max 1–5) | `test_paper22_density.py` (global_D exact-sympy + production _gaunt_ck) | `casimir_ci.py::_gaunt_ck` | BACKED-SOUND | exact integer counts (sympy 3j) + production _gaunt_ck; all 5 l_max bit-match |
| 22 | Thm 3 pair-diagonal density D_pd = 7.81/2.76/1.44/0.90/0.62% | `test_paper22_density.py` (pair_diagonal exact + production) | `potential_sparsity.py::angular_zero_count` | BACKED-SOUND | exact integer counts; all 5 l_max bit-match |
| 22 | Thm 2 potential-independence (zero pattern bit-identical across 5 potentials) | `test_paper22_density.py::test_potential_independence_lmax2` (slow) | `potential_sparsity.py::compute_eri_tensor` | BACKED-SOUND | ERI mask bit-identical + uniform density on common-orbital set |
| 22 | Table II spinor (jj-coupled) density: d_sp^FG = 25.00/8.59/6.46/5.17/4.30/3.68%, pair-diag d_sp = 25.00/4.30/2.11/1.23/0.81/0.57% (l_max 0–5) | `test_paper22_spinor_density.py` (exact jj 3j counts; l_max 3–5 slow) | self-contained jj selection rule (Dyall §9.3, Grant §7.5) | BACKED-SOUND | exact integer counts; both conventions bit-match Table II; +Q=2(l+1)², subset & d_sp≤d_sc checks. Closes run-#4 NO-TEST gap |
| 22 | §III footnote "angular_zero_count computes the global-M_L D" | `test_paper22_density.py::test_angular_zero_count_computes_pair_diagonal_not_global` | `potential_sparsity.py::angular_zero_count` | MISLABEL (CF-1) | routine computes D_pd, not global D; footnote inaccurate; regression-pinned, NOT fixed — see `docs/qa/group4.carryforward.md` CF-1 |
| 24 | Two-fermion HO entanglement rigidity: S_full=0 (single Slater det) and ‖[H_HO,V]‖_F/‖H_HO‖_F < 1e-15 (Moshinsky–Talmi N_tot conservation) | `test_paper24_ho_entropy.py` (N_max∈{2,3}) | `nuclear/ho_two_fermion.py` | BACKED-SOUND | S=0.0, rel-commutator ~2.4e-16/2.7e-16; self-contained (no archived-debug import) |
| 31 | §two_body gauged spectral action: connected fraction 77%/32% (n_max 2/3) | `test_paper31_two_body.py` (spectral_action) | `operator_system.py`, `casimir_ci.py` | BACKED-SOUND | bit-reproduces 76.7%/32.4% (double-sum inner fluctuation) |
| 31 | §two_body angular part: 100% m-conserving, pure k=0 monopole | `test_paper31_two_body.py` | `operator_system.py` | BACKED-SOUND | exact Coulomb selection-rule match (A/D positive leg) |
| 31 | §two_body radial Pearson 0.41–0.58, decreasing with n_max (honest-negative) | `test_paper31_two_body.py::...nmax3_decreasing` (slow) | `operator_system.py`, `casimir_ci.py` | BACKED-SOUND | pins 0.58→0.41 + decreasing trend (CLAUDE.md §3 dead-end) |
| 31 | §two_body Dirac resolvent Pearson 0.81/0.75 (n_max 2/3) | `test_paper31_two_body.py::test_resolvent_two_body_dirac_pearson` | `operator_system.py`, `casimir_ci.py` | BACKED-SOUND | 0.807/0.751; vs production Coulomb 0.825/0.75 (§3 dead-end) |
| 31 | §two_body Laplacian resolvent uncorrelated (N=1 zero-mode kills (1s,1s)); Dirac regularizes | `test_paper31_two_body.py::test_resolvent_laplacian_zero_mode_decorrelation` | `operator_system.py`, `casimir_ci.py` | BACKED-SOUND | \|r_lap\|≪r_dir; paper's "r<0.05" is cutoff-specific (exact at n_max=3) |

## Other branches (chemistry / QC / QED-gauge / NCG-OA / precision) — PENDING

_Pre-populated during the group3 sweep where a finding spilled into a future branch._

| Paper | Claim | Backing test | Code module | Status | Verdict |
|---|---|---|---|---|---|
| 8 | Thm 1 σ-bond selection rules: D²_(1,0),(1,0)(γ)≡1 (transparent mode), D²_(0,0),(1,0)(γ)≡0 (no s–p mixing) | `test_paper8_sigma_bond_selection.py` (10k-pt grid <1e-14 + closed forms) | `wigner_so4.py::wigner_D_so4` | BACKED-SOUND | machine-precision on the paper's grid + matches in-paper cos²+sin² / antisym closed forms. Replaces stale `debug/test_harmonic_phase_lock.py` cite (archived dead-end probe) |
