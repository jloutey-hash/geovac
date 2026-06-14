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

| Paper | Claim | Backing test | Code module | Status | Verdict |
|---|---|---|---|---|---|
| 0 | Annular area → 2ℓ+1 angular capacity | `test_trunk_qa_annular.py` | `lattice.py` | BACKED-SOUND | ✅ accepted (new) |
| 0 | Cumulative shell count Σ(2k−1)=n² | (exact, agent-verified) | `lattice.py` | BACKED-SOUND | sound |
| 0,7 | **κ = −1/16 "derived from Fock projection"** | `test_trunk_qa_kappa.py` | `lattice.py`, `atomic_solver.py` | **OBSERVATION** (was "SYMBOLIC PROOF") | ⬇ **downgraded** — matching scale *coincides* with geometric 1/16; counterfactual shows no bridge |
| 0 | \|V\|=2n² spin-doubled vertex set | none | `dirac_lattice.py` (`n²`, spinless) | NO-TEST (code/paper cardinality mismatch) | flag — reconcile |
| 1 | O(V)/O(N) sparse eigen-scaling | `test_ov_scaling_rigorous.py` | `lattice.py` | BACKED-SOUND | sound (multi-point fit ~1.0) |
| 1,7 | Spectrum = −(n²−1) | `test_fock_laplacian.py` (continuum, n≤3) | symbolic S³ L–B | BACKED-SOUND (continuum only) | ⚠ the **discrete graph L=D−A is positive-semidefinite**; −(n²−1) is a *continuum* property, not the graph's |
| 1,7 | 2s/2p splitting → 0 (convergence) | `test_trunk_qa_splitting.py` | `lattice.py` | BACKED (oscillatory decay) | ✅ accepted (new) |
| 1 | §III degree table D_2s=0.854, D_2p=3.416 | `test_trunk_qa_degree_table.py` | — (CG construction absent) | DOWNGRADED | ⬇ faithful build gives 1.9/3.8; numbers not reproduced |
| 1 | §III 16% splitting @ n_max=10 | `test_trunk_qa_splitting.py` | `lattice.py` | DOWNGRADED | ⬇ real ~1.7% @ n=10 (~16% lands near n=8) |
| 1 | Gearing ratio ‖L₊‖/‖T₊‖ → 1.77 | `test_trunk_qa_gearing.py` | — | DOWNGRADED | ⬇ not convergent (spectral norm = 2.0) |
| 1 | [T₊, L₊] = 0 | `test_trunk_qa_gearing.py` | `lattice.py` | BACKED-SOUND | ✅ accepted (bit-exact) |
| 7 | 18 symbolic Fock-S³ proofs (continuum geometry + S³ eigenvalues n≤3) | `test_fock_projection.py`, `test_fock_laplacian.py` | symbolic | BACKED-SOUND | sound; ~4 of 18 weak/tautological (volume Jacobian, etc.) — soften "18 independent" |
| 7 | c²(4,3) = 1/40 = Δ (Paper 2 boundary term) | `test_trunk_qa_c2_delta.py` | `dirac_s3.py` | BACKED-SOUND | ✅ accepted — two independent routes to 1/40 (genuine bridge) |
| 7 | V_ee Slater integrals F⁰ on S³ | `test_algebraic_exchange.py` | `algebraic_slater.py` | BACKED-SOUND | sound |
| 32 | Connes axiom audit (J²=−1, JD=+DJ, order, Lorentzian signs) | `test_real_structure.py`, `test_connes_axiom_audit_31.py` | `real_structure.py` | BACKED-SOUND | sound — bit-exact + negative controls |
| 32 | K = π(B+F−Δ) | `test_paper2_corrections.py` (arithmetic only) | — | COMPLIANT (conjectural) | ✅ §13.5 honored in code — no test asserts derivation |
| 32 | Forced-Count moduli chain 1024→512→128→128 | `test_trunk_qa_forced_count_moduli.py` | `standard_model_triple.py` | DOWNGRADED | ⬇ full-axiom = 260, not 128; 128 = matter-sector; "45 tests verify" mis-pointed |
| 32 | §VIII non-selection theorems (N_gen, KO-dim, H1 Yukawa, single-cutoff, …) | none | structural prose | PROOF-BY-ARGUMENT | not test candidates; verify the argument (paper review) |
| 38 | "Unconditional" state-space GH convergence | `test_p38_action_seminorm.py` | `full_dirac_operator_system.py` | QUALIFIED-SOUND | honest to the translation-seminorm scope |
| 38 | Per-band injectivity (full N² rank) | `test_p38_action_seminorm.py`, `p38_g1g2_band_diagnostics.py` | `spinor_operator_system.py` | BACKED-WEAK (finite cutoff n≤5) | general leg = Schur + AWA closed-form |
| 38 | **4/π asymptotic constant** | `test_trunk_qa_fejer_4_over_pi.py` | `central_fejer_su2.py` | BACKED (derived, non-circular) | ✅ accepted — rejects 2/π decoy; **replaces circular `test_asymptotic_constant_value`** |
| 38 | γ_n closed-form sum rule | `test_central_fejer_su2.py` | `central_fejer_su2.py` | BACKED-SOUND | sound |

**Pattern:** the trunk's machinery (axiom audit, O(V), continuum symbolic proofs, sum rules, K-conjectural) is genuinely backed; the soft spot was specifically the *headline derived/exact constants* — κ (coincidence, downgraded) and 4/π (genuinely derivable, upgraded from a hardcoded circular test). Plus Paper 1 §III quantitative content (from a CG construction never in production).

## Foundations branch (Papers 18, 22, 24, 31, 54–57) — PENDING

## Other branches (chemistry / QC / QED-gauge / NCG-OA / precision) — PENDING
