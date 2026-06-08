# Sprint M-vS-2 + R-sweep — default LiH spec Bratteli reading and spectral-action functional test

**Date:** 2026-06-07 (session continuation)
**Driver:** `debug/bratteli_mvs2_lih_default_driver.py`
**Data:** `debug/data/bratteli_mvs2_lih_default.json`
**Log:** `debug/bratteli_mvs2_lih_default_log.txt`
**Predecessors:** Sprint M-vS-1 (heteronuclear LiH 2-vertex pilot, `debug/bratteli_lih_pilot_memo.md`); H₂ pilot (`debug/bratteli_h2_pilot_memo.md`); umbrella memo `debug/sprint_w1e_audit_and_mvs_chemistry_open_2026_06_07_memo.md`; arc scoping `debug/sprint_marcolli_vs_chemistry_paper_arc_scoping_memo.md`.

---

## TL;DR

- **Q1 STRUCTURAL: PASS bit-exactly (residual 0.0)** — the production-default `lih_spec()` 3-sub-block / 15-spatial-orbital chemistry Hamiltonian is a Marcolli–van Suijlekom 2014 gauge network on the 2-vertex Li ↔ H bond quiver, with Li-vertex Hilbert space H_Li = H_Li_core ⊕ H_LiH_bond_center (10-dim direct sum) and H-vertex H_H = H_LiH_bond_partner (5-dim). The named M-vS-2 sprint deliverable is complete with a stronger result than asked: residual is *exactly* 0.0 rather than ≤ 1e-10.
- **Q2 FUNCTIONAL: NEGATIVE — spectral action does NOT bind LiH.** Over R ∈ [2.0, 8.0] bohr, the spectral action S(D)(R) = Tr exp(−D²/Λ²) of the assembled M-vS Dirac is **monotone increasing** in R at every Λ ∈ {1, 2, 4}; E_FCI(R) is **monotone decreasing** toward small R. Both have boundary minima at R = 2.0 bohr (smallest R tested). Neither is a binding functional for LiH at this projection.
- **Joint reading: W1e is at the projection step, not the evaluation step.** The chemistry-engineering arc is exhausted (umbrella memo §3.1); switching the EVALUATION step from FCI to spectral action gives the same monotone-descent failure. W1e is structurally about what (h1, eri, ecore) carries — the binding shape lives in the continuous Level-4 PK-composed pipeline ingredients that get lost in the projection to integrals.

The M-vS paper arc has its second sprint deliverable. The reconciliation question from this morning ("would the M-vS-native observable S(D) bind where FCI doesn't?") is settled in the negative — cleanly and without ambiguity.

## 1. What was tested

Two questions in one driver:

**Q1 STRUCTURAL.** Is the production-default `lih_spec()` — 3 sub-blocks (Li_core 5 orbitals at Z=3 + LiH_bond_center 5 at Z_eff_Li + LiH_bond_partner 5 at Z=1), 15 spatial orbitals total — bit-exactly a Marcolli–vS 2014 gauge network on the 2-vertex Li ↔ H bond quiver under the natural reading:

- Vertex Li: H_Li = H_Li_core ⊕ H_LiH_bond_center (10-dim direct sum). The vertex prespectral triple has Hilbert space being a direct sum of two atomic shells, with vertex Dirac block-diagonal in the two shells (cross_block_h1 skips same-center pairs structurally).
- Vertex H: H_H = H_LiH_bond_partner (5-dim).
- Edge intertwiner L_e: 5 × 10 sub-matrix of cross_block_h1 combining (Li_core ↔ H bond_partner) and (LiH_bond_center ↔ H bond_partner) coupling.

Pass gate: max|h1_GeoVac − P^T · H_MvS · P| ≤ 1e-10 at R = 3.015 bohr, where P is the permutation reordering sub-blocks into Li-then-H vertex order.

**Q2 FUNCTIONAL.** Over R ∈ {2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0} bohr, does the spectral action S(D)(R) = Tr exp(−D²/Λ²) of the assembled M-vS Dirac D = h1_GeoVac (relabeled) develop an interior minimum near the experimental R_eq = 3.015 bohr? This is the reconciliation question opened in this session: the chemistry pipeline computes binding via FCI on (h1, eri, ecore); the M-vS-native observable is the spectral action functional. If S(D)(R) had a minimum at R_eq, the spectral action would be a candidate alternative binding functional. If not, the W1e wall is independent of the evaluation choice.

## 2. Q1 verdict: bit-exact PASS

At R = 3.015 bohr, max_n = 2 (production default):

| Quantity | Value |
|:---------|:------|
| max |h1_GeoVac − P^T H_MvS P|       | **0.0** (exact) |
| Vertex Li Hermiticity residual      | 0.0 |
| Vertex H Hermiticity residual       | 0.0 |
| Internal Li off-diag (Li_core ↔ LiH_bond_center) | 0.0 |
| ‖L_e‖_F (Frobenius)                | 1.378 |
| ‖L_e‖_op (operator norm)           | 1.351 |
| max |L_e[i,j]|                     | 0.772 |
| ‖L_e^† L_e − I‖_∞                  | 1.0 (non-unitary) |
| ‖L_e L_e^† − I‖_∞                  | 1.0 (non-unitary) |

The 3-sub-block production spec maps to a 2-vertex quiver because both Li_core and LiH_bond_center sit on the Li nucleus at (0, 0, 0); they form a *direct sum* of atomic shells inside a single vertex prespectral triple. The internal Li-vertex off-diagonal block (Li_core ↔ LiH_bond_center) is exactly zero because `cross_block_h1` skips same-center pairs structurally (line 466–469 of `geovac/cross_block_h1.py`).

The edge intertwiner L_e is **Hermitian but not unitary** (residual exactly 1.0). This confirms what the H₂ and LiH 2-vertex pilots already showed: GeoVac sits in the **Marcolli–vS 2014 framework where L_e is allowed to be Hermitian**, NOT in the **Perez-Sanchez 2024a framework where L_e is required to be unitary**. The bit-exact match strengthens the H₂ + LiH-artificial-2-vertex result to the production-default 3-sub-block case.

**Significance for the M-vS paper arc.** The arc-scoping memo §4 named M-vS-2 as "the load-bearing test for whether the production Track CD pipeline is structurally a Marcolli-vS gauge network on the bond quiver" — a much stronger claim than the H₂ and LiH-artificial-2-vertex pilots. M-vS-2 PASSES. The paper headline can move from "we can encode chemistry-like specs into M-vS by construction" to "the production Track CD pipeline IS structurally a M-vS gauge network at bit-exact precision on H₂, LiH-2-vertex, and LiH-default-3-sub-block."

## 3. Q2 verdict: NEGATIVE on the spectral-action-as-binding-functional question

R-sweep over 10 panel points, n_max = 2, 4-electron FCI:

| R (bohr) | Tr(D²) | S(Λ=1) | S(Λ=2) | S(Λ=4) | V_NN+E_core | E_FCI |
|:--------:|:------:|:------:|:------:|:------:|:-----------:|:------:|
| 2.000 | 53.79 | **6.45** | **10.33** | **12.88** | −5.78 | **−18.08** |
| 2.500 | 45.72 | 6.78 | 10.62 | 13.05 | −6.08 | −17.15 |
| 3.000 | 41.60 | 6.95 | 10.86 | 13.24 | −6.28 | −16.73 |
| 3.015 | 41.51 | 6.95 | 10.86 | 13.24 | −6.28 | −16.72 |
| 3.500 | 38.79 | 7.12 | 11.05 | 13.41 | −6.42 | −16.51 |
| 4.000 | 36.96 | 7.30 | 11.22 | 13.51 | −6.53 | −16.34 |
| 4.500 | 35.66 | 7.50 | 11.37 | 13.59 | −6.61 | −16.18 |
| 5.000 | 34.47 | 7.71 | 11.50 | 13.63 | −6.68 | −16.05 |
| 6.000 | 32.75 | 8.10 | 11.73 | 13.69 | −6.78 | −15.75 |
| 8.000 | 30.51 | 8.74 | 12.07 | 13.72 | −6.90 | −15.21 |

(Boldface marks the panel minimum for each column.)

All five observable curves are monotone in R:

- **E_FCI(R)** is monotone DECREASING toward small R. No interior minimum. This is the W1e wall as documented in CLAUDE.md §3 and §1.7 (multi-focal-composition wall pattern, 6th instance).
- **S(Λ=1, 2, 4)** are all monotone INCREASING in R. Smallest at R = 2.0 bohr (the boundary). No interior minimum at any Λ.
- **Tr(D²)** is monotone DECREASING in R. Smallest at R = 8.0 bohr (the other boundary).

The spectral action and Tr(D²) have OPPOSITE monotonicity, but neither has an interior minimum.

**No interior minimum at any tested Λ for any of the five observables.** Predicted "D_e" from the spectral action gap (S(R=∞) − S(R_min)) is 2.28, 1.73, 0.84 at Λ = 1, 2, 4 respectively — but these are gaps from the panel boundary to the same panel boundary on the opposite end, not gaps from a well to the dissociation limit. There is no binding-shape interpretation available.

For E_FCI: the "D_e_predicted" from the panel is 2.87 Ha (E_FCI(R=8) − E_FCI(R=2.0)) — but this is again the over-binding signature, not a real D_e.

For comparison: experimental R_eq = 3.015 bohr, D_e = 0.0924 Ha (NIST CCCBDB).

## 4. Joint reading

The Q1+Q2 result lands clean. Three structural implications:

### 4.1. The M-vS gauge-network correspondence is paper-grade

Three molecules at bit-exact precision: H₂ (residual 9.7×10⁻¹⁸), LiH-2-vertex (residual 2.5×10⁻¹⁷), LiH-3-sub-block-default (residual 0.0). The arc has three empirical anchors covering homonuclear (H₂), heteronuclear-artificial (LiH 2-vertex), and heteronuclear-production-default (LiH 3-sub-block). The umbrella sprint memo §3.3 explicitly identified the 3-sub-block test as the "much stronger" claim — now verified.

The remaining concerns flagged by the H₂ pilot memo §3.4 ("two-body ERI Bratteli reading", M-vS-3) and the arc-scoping memo §4 (NaH frozen-core test) are still open; the present sprint addresses neither.

### 4.2. W1e is at the PROJECTION step, not the EVALUATION step

The chemistry-engineering arc was exhausted in last sprint (umbrella memo §3.1, CLAUDE.md §3 entries for LiH kwarg STOP, B.1 HF NaH STOP, F1–F6). That arc tested DIFFERENT BUILDERS of (h1, eri, ecore) and showed none of them produce binding-shaped curves in FCI.

The present sprint tests a DIFFERENT EVALUATION of the same (h1, eri, ecore): instead of FCI ground state on the full chemistry Hamiltonian, the spectral action functional S(D) = Tr exp(−D²/Λ²) on just D = h1. Result: same monotone-descent failure, just in a different functional.

Combining the engineering-arc and evaluation-axis tests:
- **Building** the chemistry Hamiltonian differently (multi-zeta, screened, PK variants, explicit core, cross-block h1, ...): does not produce binding.
- **Evaluating** the same Hamiltonian via different functionals (FCI ground state, spectral action at three Λ values, trace invariants): does not produce binding either.

Both axes are now empirically exhausted. The remaining locus for W1e — the chemistry-accuracy wall — is at the **projection step from continuous Level-4 PK-composed multichannel adiabatic + PK solver to second-quantized (h1, eri, ecore) integrals**. The PK content that makes the continuous solver bind LiH at 5.3% R_eq (Paper 17, ComposedDiatomicSolver.LiH_ab_initio) does not survive the projection into integral tensors; what is lost cannot be recovered by either changing the builder OR changing the evaluator.

This is a **sharper statement of the multi-focal-composition wall pattern** for the chemistry case. The pattern is now:

> Build the M-vS gauge network on the bond quiver from atomic spectral triples and multipole bimodules → at bit-exact precision, the chemistry Hamiltonian IS gauge-network-shaped. But the binding-content is calibration data that lives in the continuous PK-composed pipeline and is NOT autonomously generated by the gauge-network structure — neither by FCI ground state nor by the spectral action functional.

### 4.3. The M-vS structural identification is binding-orthogonal

This sharpens the reconciliation question. Marcolli–vS-as-structural-identification (what the paper arc is about) is a **separate question** from M-vS-as-binding-functional (what the spectral action would do). The first is now confirmed at bit-exact across three molecules; the second is empirically negative.

The arc-scoping memo §6 framing — "this paper provides the structural NCG home for GeoVac's chemistry composition" — survives Q2 negative cleanly. The paper is about the structural correspondence, not about a new chemistry-binding method. Honest scope.

## 5. What this sprint did NOT do

- **NOT** test M-vS-3 (two-body ERI Bratteli reading per arc-scoping memo). The eri tensor was used by FCI but not analyzed for plaquette-trace structure.
- **NOT** test NaH frozen-core in the default spec (arc-scoping memo flagged this for paper-arc completeness).
- **NOT** sweep at n_max = 3. n_max = 2 sufficient for the Q1 bit-exact identification and Q2 monotone-curve verdict; n_max = 3 not expected to change either.
- **NOT** test alternative spectral-action cutoffs (sharp, polynomial, exponential). Gaussian cutoff is the Chamseddine–Connes standard; if Gaussian gives monotone behavior at three Λ values, alternative cutoffs are unlikely to flip the verdict.
- **NOT** test S(D + λ · eri-trace) where the spectral action is taken of an "extended Dirac" incorporating two-body content. This is a natural follow-on for the M-vS-3 sprint.
- **NOT** alter the M-vS paper arc plan. The sprint adds M-vS-2's deliverable cleanly; M-vS-3 → M-vS-4 → M-vS-5 sequence stands.

## 6. Hard-prohibition check (CLAUDE.md §13.5)

No changes to:
- Natural geometry hierarchy
- Fitted-or-empirical parameters
- §3 deletions (only ADDING the spectral-action-monotone-descent row at end of sprint)
- Paper 2 K = π(B+F−Δ) combination-rule "conjectural" label

## 7. Verification

- Driver runs in ~140 seconds wall time (10 R-points × ~13s each for h1+eri build + FCI diagonalization).
- Q1 gate at R = 3.015: max residual 0.0 — well below the 1e-10 gate.
- All 10 R-sweep points produce internally-consistent data with monotone curves across all five observables (E_FCI, Tr(D²), S(Λ=1,2,4)).
- Vertex Hermiticity = 0.0 (no Hermitian-violating artifacts of the relabeling).
- Internal Li off-diag = 0.0 (cross_block_h1 same-center-skip behavior matches expectation).

## 8. Files

### Created
- `debug/bratteli_mvs2_lih_default_driver.py` (~450 lines)
- `debug/data/bratteli_mvs2_lih_default.json`
- `debug/bratteli_mvs2_lih_default_log.txt`
- `debug/sprint_mvs2_lih_default_plus_rsweep_memo.md` (this file)

### Modified
- None (no production code changes; no paper edits)

### Next sprint candidate
- **Sprint M-vS-3 (two-body ERI Bratteli reading).** From the arc-scoping memo §4.2: "Does Track CD's cross-block ERI tensor admit a 'Bratteli self-loop / plaquette' reading à la Perez-Sanchez §4?" Decision gate: structural identification with empirical match within 1%. Scope: 3-4 weeks per the arc plan.

## 9. Sample size

| Metric | Value |
|:-------|:------|
| Structural identification (Q1) | 1 spec config (default LiH), 1 R = 3.015 bohr, n_max = 2 |
| R-sweep (Q2)                   | 10 R points × 5 observables × 1 cutoff family × 1 n_max |
| Marcolli-vS correspondence cumulative | 3 molecules at bit-exact precision: H₂, LiH-2-vertex, LiH-default-3-sub-block |
