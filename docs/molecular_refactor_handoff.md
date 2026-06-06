# Molecular subsystem refactor — handoff packet

> **RETIRED 2026-06-05 (PI directive, same-day).** The planned molecular refactor is **not needed.** Sprint chemistry-test-rot reconciliation (v3.56.0, `debug/sprint_lih_binding_fix_memo.md`) closed the load-bearing items: §2.1 LiH binding collapse resolved by F.1 V_NN(R) fix + F.2 PK pseudopotential restore (LiH ab initio R_eq 2.82%, balanced 6.93%, 943/943 tests pass); §2.2 nuclear-electronic commutativity and §2.3 BeH₂ cross-block h1 dissolved on diagnosis (G + H sub-tracks). The §4 shim audit and §3.1 Paper 17 reconciliation remain as low-priority follow-ons but do not require a sprint-scale refactor. This document is preserved for institutional memory; treat it as a historical snapshot of the v3.52.0–v3.55.0 problem statement, not as an active to-do.

**Prepared:** 2026-06-05 (close of v3.52.0 sprint cycle)
**Audience:** the next session that picks up the planned molecular refactor.
**Purpose:** hand off the chemistry-test-rot situation cleanly so the refactor inherits a problem statement rather than having to re-diagnose.

## Context

The v3.52.0 sprint cycle did a substantial test-cleanup pass (`debug/sprint_test_cleanup_memo.md`) on the chemistry subsystem. Roughly 40 test files were fixed in-place, 11 were archived to `tests/_archive/superseded/` (or `dead_ends/`), and 10 production-code shims were added as backward-compat bandaids. **The non-molecular project is green:** 313/313 paper-equation tests pass, every spectral-triple / QED / GH-convergence / topo / Wigner test passes, and the partial full-suite scan (first 46%) showed exactly 1 failure. The rot is entirely concentrated in the molecular chemistry subsystem (composed_qubit / balanced_coupled / LCAO legacy / PK pipeline / frozen-core builds).

Because a molecular refactor is planned, fixing tests against soon-to-be-replaced APIs is wasted motion. The molecular cleanup was therefore stopped after the sub-agent's pass and the remaining work is being handed off here. **What this doc preserves is the things that would otherwise go silent.**

Three categorically different things are kept distinct below; conflating them is the failure mode this doc exists to prevent.

---

## §1. Test-rot candidates (~30–60 files)

Tests that fail because the production API moved under them. Pure mechanical cleanup. The refactor decides what to keep, rewrite against the new API, or delete.

The sub-agent's per-file table is the authoritative source:\ `debug/sprint_test_cleanup_memo.md` §"Per-file triage table" (lines 114–168).

The dominant patterns were classified A–E in that memo:

| Pattern | What | Treatment | Count |
|:-------:|:-----|:----------|:-----:|
| A | API drift / removed kwarg | Fix-in-place (drop kwarg) or production shim (kwarg as no-op) | ~12 files |
| B | Removed function / class import | Redirect to current API or archive | ~10 files |
| C | Block-layout / result-dict shape change | Update test to use new shape | ~5 files |
| D | Numerical drift from v3.x sprints | Update expected value or relax tolerance, NAMED if physics-real | ~6 files |
| E | NAMED PRODUCTION REGRESSION (LiH R_eq) | Tests skipped pending PI decision | ~15 test methods across 4 files |

**Archived in v3.52.0 cleanup** (don't re-archive; just verify the new API doesn't restore them):

```
tests/_archive/superseded/
  test_algebraic_pk_superseded.py            # AbInitioPK.algebraic_projector removed
  test_algebraic_zeff_superseded.py          # zeff_method= kwarg removed
  test_h2_energy_decomposition_superseded.py # LCAO MoleculeHamiltonian path
  test_lih_validation_superseded.py          # LCAO MoleculeHamiltonian path
  test_locked_shell_superseded.py            # MolecularLatticeIndex retired
  test_locked_shell_hyperspherical_superseded.py  # MolecularLatticeIndex retired
  test_molecular_bugs_superseded.py          # LCAO bridge connectivity tests
  test_prolate_active_superseded.py          # LockedShellMolecule + removed function
  test_spectral_pes_wiring_superseded.py     # radial_method='spectral' retired
  test_track_r_benchmark_superseded.py       # multiple removed v2.0.x APIs

tests/_archive/dead_ends/
  test_spectral_pk_superseded.py             # spectral_rank_k_projector retired
                                              # (one of the 6 documented PK failures, CLAUDE.md §3)
```

**Two TBD files** the cleanup sub-agent didn't get to (token-limit cutoff): `test_level4_multichannel.py`, `test_prolate_heteronuclear_scf.py`. The refactor should run these first as a "is this still relevant" gate.

---

## §2. Production bugs guarded by rotting tests

These are physics-correctness regressions, NOT test rot. If the refactor archives the tests without separately addressing the underlying bugs, the bugs go silent until they cause a paper claim to fail downstream.

### §2.1 LiH binding collapse in both `balanced_coupled` and `LiH_ab_initio` (Pattern E)

**Symptom.** Both production paths give monotone-decreasing PES with no equilibrium minimum near experimental R_eq = 3.015 bohr:

| Path | E(R=1.5) | E(R=3.015) | E(R=7.0) | PES minimum |
|:-----|:--------:|:----------:|:--------:|:------------|
| `ComposedDiatomicSolver.LiH_ab_initio(l_max=2)` | −1.821 | −1.323 | (n/r) | At R=1.5 (panel boundary) |
| `build_balanced_hamiltonian(lih_spec, R)` | **−16.040** | −15.205 | −14.172 | At R=1.5 (panel boundary) |
| `build_composed_hamiltonian(lih_spec, R)` | −13.138 | −14.143 | −14.838 | Trivial: E_elec ≈ −15.138 constant; total tracks only −3/R |

Composed_qubit's R-independence is by-design (no cross-center V_ne). The balanced_coupled regression is the real issue.

**Diagnostic drivers and data:**
- `debug/diag_lih_balanced_pes.py` + `debug/data/diag_lih_balanced_pes.json`
- `debug/diag_lih_composed_qubit_pes.py` + `debug/data/diag_lih_composed_qubit_pes.json`
- `/tmp/lih_probe.log` (legacy `LiH_ab_initio` probe)

**Tests currently skipped pending PI decision** (per sub-agent memo §"Pattern E"):
- `tests/test_ab_initio_pk_v2.py::test_lih_r_eq_near_experiment`
- `tests/test_ab_initio_pk.py::test_r_eq_physical`
- `tests/test_ab_initio_pk.py::test_beh_r_eq_physical`
- `tests/test_l_dependent_pk.py::TestLiH{L2,L3,L4}*::test_r_eq_*` (3 skipped)
- `tests/test_composed_diatomic.py::TestLiHL2Backward, TestLiHL3Accuracy, TestLiHL3SigmaPi, TestLiHL4SigmaPi, TestBeHPlusPipeline, TestArchitectureGenerality, TestComparisonTable, TestLmaxConvergenceReport` (8 class-level skips)

**Diagnosis question for the refactor:** git-bisect through v3.x to identify which sprint broke balanced_coupled LiH binding. The Pauli count (878) is preserved across the regression, so the structural architecture is intact; the broken piece is in the cross-center V_ne assembly or in the (h1, eri) integration. Track CD (v2.0.39+) documented "LiH 878 Pauli, 1.8% energy, 7.0% R_eq; only bound 4e config" — the regression is *from* that state.

**Anti-pattern to avoid in the refactor:** rebuilding on top of the current balanced_coupled code without verifying the binding curve is restored. If the new architecture inherits the same (h1, eri) machinery without revisiting it, the regression carries forward silently.

### §2.2 `nuclear_electronic` composed-matrix commutativity (0.375 Ha residual)

**Symptom.** `geovac/nuclear/nuclear_electronic.py::build_deuterium_composed_matrix(include_hyperfine=False)` produces a Hamiltonian whose `[H_nuc, H_el]` commutator has max coefficient ~0.375 Ha. Previously it was bit-zero (< 1e-6 Ha residual) when hyperfine was off.

**Test tolerance:** `tests/test_nuclear_electronic.py` had its tolerance relaxed from 1e-6 → 1.0 Ha by the cleanup pass to mask the failure. The test now passes silently; the production bug is undocumented in `geovac/`.

**Hypothesis:** `include_hyperfine=False` is not actually zeroing all cross-register coupling. Either: (a) a coupling block is added unconditionally and the `include_hyperfine` flag is checked too late, or (b) the cross-register block has a non-zero contribution that survives the flag.

**Diagnosis question for the refactor:** audit `build_deuterium_composed_matrix` and confirm `include_hyperfine=False` truly returns a block-diagonal Hamiltonian. If yes, the regression is in some auxiliary function being called unconditionally. If no, the production code has been silently coupling hyperfine even when callers asked it not to.

### §2.3 BeH2 cross-block h1 disagreement (~5%)

**Symptom.** The old `build_composed_beh2(pk_in_hamiltonian=...)` builder and the new `build_composed_hamiltonian(beh2_spec)` produce h1 matrices that disagree by ~5% on diagonal entries.

**Test:** `tests/test_composed_beh2.py::test_h1_match` currently skipped.

**Question:** is the new builder's cross-block h1 addition (W1d arc per CLAUDE.md §1.7) the correct one, or is the old builder correct? If the new is correct → remove the deprecated `build_composed_beh2` entirely. If the old is correct → the cross-block addition is a bug.

This is a clean either-or; pick one and remove the other.

---

## §3. Paper claims at risk

### §3.1 Paper 17 §V Table II line 1154: "5.3% R_eq for LiH via Adiab. + l-dep PK"

The legacy `ComposedDiatomicSolver.LiH_ab_initio(l_max=2)` path that produced this number is broken (see §2.1). The number is currently NOT reproducible from any path in the codebase. Three options for the refactor:

| Option | Action | Paper-honesty cost |
|:-------|:-------|:-------------------|
| (a) Restore | Debug the regression, restore v2.x behavior, re-verify 5.3% | Zero — paper stays as published |
| (b) Reframe | Add Paper 17 footnote noting the v3.x refactor changed the result | Moderate |
| (c) Erratum / archive the claim | Remove or pin the 5.3% claim to a historical commit | High — paper claim becomes unmoored |

The cleanup sub-agent's memo §"Named follow-ons" item 1 calls out this decision explicitly. Until the refactor decides, the tests stay skipped with explanatory markers and Paper 17 stays as published. **The refactor should make this decision in the same sprint** where it touches LiH binding — not as a separate cleanup pass — so the paper stays consistent with the code at the point the new architecture lands.

### §3.2 Paper 27 HO ground-state energy: 14.898 → 22.185 MeV (Minnesota V_S/V_T fix)

The v3.38.0 Minnesota fix changed the expected GS energy. The cleanup sub-agent took 22.185 MeV as the new authoritative value and updated the test; the value was not independently verified against published Minnesota interaction parameters.

**Diagnosis question for the refactor:** if the molecular refactor touches the nuclear subsystem at all, the 22.185 MeV value should be verified against Minnesota literature parameters. If the published-parameter calculation gives something else, the v3.38.0 fix has its own bug.

### §3.3 BeH2 monopole R-flatness tolerance relaxed 15% → 30%

`tests/test_composed_beh2.py::test_monopole_r_flatness` had its tolerance widened from 15% to 30% to pass. This is a qualitative regression (still flat-ish) but the original tight bound was a meaningful structural property. The refactor should consider whether to restore the original ~15% behavior.

---

## §4. Backward-compat shims that should be audited / removed

The cleanup sub-agent added 10 production-code shims as bandaids to keep test fixtures running. These should be audited and either:\ (a) preserved as legitimate backward compatibility, (b) removed because the underlying caller pattern is gone in the new architecture, or (c) recognized as silently-ignored kwargs that should error instead.

| File:line | Shim added | What it actually does | Risk of leaving in place |
|:----------|:-----------|:----------------------|:-------------------------|
| `geovac/inter_fiber_coupling.py:321, 178` | `extract_channel_data` / `extract_origin_density` no longer propagate `pk_potentials` to `solve_angular_multichannel` | PK info silently dropped | Production callers (`composed_water`, `composed_triatomic`) think they're using PK; they're not |
| `geovac/level4_multichannel.py:1125` | `solve_level4_h2_multichannel(pk_potentials=...)` accepted as no-op | Kwarg silently ignored | Same — callers think PK is applied |
| `geovac/algebraic_coupled_channel.py:43-45` | `solve_radial_spectral` / `solve_coupled_radial_spectral` import shim | Function-not-found at call time | Caller path is dead; safe to remove if no caller remains |
| `geovac/level4_spectral_angular.py:39-47` | `compute_core_screening_analytical` / `compute_pk_pseudopotential` import shim | Same — NotImplementedError at call time | Same |
| `geovac/cusp_factor.py:53-61` | Same import-shim issue | Same | Same |
| `geovac/rho_collapse_cache.py:348-352` | `build_angular_hamiltonian` was being called with removed `Z_A_func=` and `n_theta=` kwargs | Real bug fix — caller now matches current signature | Keep |
| `geovac/frozen_core.py:586+` | `MolecularLatticeIndex` (LCAO) `NotImplementedError` stub | Errors with clear message if invoked | Keep — but remove entire path if LCAO is gone |
| `geovac/locked_shell.py:147` | Same `MolecularLatticeIndex` shim | Same | Same |
| `geovac/composed_qubit.py:1872, 2443` | `build_composed_beh2/_h2o` accept `pk_in_hamiltonian=` kwarg | Silently ignored if the new builder doesn't honor it | Verify behavior matches old builder |
| `geovac/hyperspherical_radial.py` | `_laguerre_moment_matrices` and `_build_laguerre_SK_algebraic` restored | Real restore — CLAUDE.md §12 algebraic-registry entries depend on these | Keep |

**Highest-risk shims:** the first three. They silently absorb a `pk_potentials` argument and drop it on the floor. Production callers in `composed_water.py` and `composed_triatomic.py` still pass `pk_potentials=self.pk_potentials` thinking it has effect; it doesn't. The refactor should decide whether to:\ (a) make these kwargs raise an explicit error (so callers update consciously), or (b) actually wire PK through to `solve_angular_multichannel` again (restoring the original architecture).

---

## §5. What the refactor should NOT need to re-do

These are settled in v3.52.0 and the refactor inherits clean:

- **Z_2 Hopf-U(1) tapering** (`geovac/z2_tapering.py`, Paper 14 §sec:hopf_tapering, `tests/test_z2_tapering.py` 17/17 green). Both `global` and `per_block` modes work for all 37 molecules. Backward-compatible `tapered={None, 'global', 'per_block'}` flag on `ecosystem_export.hamiltonian()`.
- **`extract_channel_data` and `extract_origin_density`** no longer crash when `pk_potentials` is passed. (Bug-fix retained as a feature in §4 above — but the silent-drop semantics need reconsideration.)
- **All test files in `tests/` that did NOT move to `tests/_archive/`** import without error and have at least their structural / fixture-setup tests passing.
- **Non-molecular tests** (QED, spectral triple, GH-convergence, paper-equation, topo, Wigner, Ihara, Berezin, real-structure, hypergeometric_slater, Breit) confirmed green during the v3.52.0 cycle. The refactor shouldn't need to revisit any of these.

---

## §6. Suggested refactor sprint shape

1. **Diagnose §2.1 (LiH binding)** via git-bisect on `balanced_coupled` LiH PES. This is the load-bearing decision; everything else cascades. Likely 1–3 days.
2. **Make the Paper 17 call** in the same sprint where (1) lands. Don't let the codebase and paper drift apart even one more release cycle.
3. **Audit §4 high-risk shims** as part of the architecture review. Decide silent-drop vs explicit-error vs restore-original.
4. **Run §2.2 and §2.3 audits** as quick checks (each ~30 min). Likely independent of (1).
5. **Re-test the §3.2 Paper 27 GS energy** against literature if the refactor touches nuclear physics at all.
6. **Refresh the test suite** against the new API: unarchive any tests that the new API restores, write fresh tests for any new architectural surface.

After steps 1–5 land, the molecular subsystem rejoins the non-molecular project at "clean green." Step 6 is the test-rewrite work this handoff was created to defer; it's the right work to do *once*, against the new API, not against the current rotting one.
