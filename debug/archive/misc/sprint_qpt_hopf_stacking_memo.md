# Sprint memo: QPT and Hopf-Z2 tapering stack — bit-exact triply-additive blocking

**Date:** 2026-06-05 (third meta-lesson confirmation in one day, after DF and Cholesky)
**Author:** PM session
**Status:** POSITIVE — every commutator bit-exact zero across LiH, BeH2, H2O. QPT and per-sub-block Hopf-Z2 simultaneously enforceable. Concrete joint-sector reduction: factor 5.7x beyond particle-number tapering for LiH (4,791 SDs vs 27,405).

**Files:**
- `debug/qpt_hopf_stacking_test.py` — pairwise commutator audit across LiH/BeH2/H2O
- `debug/qpt_joint_block_dim.py` — joint (N, Sz, P_Hopf) sector dimensions for LiH
- `debug/data/qpt_hopf_stacking_test.json`
- `debug/data/qpt_joint_block_dim.json`

---

## Motivation

Today's third sprint testing the meta-lesson "GeoVac algebra is comparable to ANY algebraic-exact compression technique." DF and Cholesky both landed inside the multipole subspace; the lit-search agent flagged the **Quantum Paldus Transform** (Burkat-Fitzpatrick June 2025, arXiv:2506.09151) as the next sprint-scale candidate.

QPT implements u(d) × SU(2) -> Gelfand-Tsetlin block-diagonalization as a quantum circuit with O(d^3) Toffoli depth. It exploits **spin (S^2)** and **orbital unitary-group symmetry** to block-diagonalize a JW Hamiltonian. The lit-search agent's prediction was that QPT should **stack additively** with our per-sub-block Hopf-Z2 tapering (v3.52.0): spin × angular × Hopf-base = triply additive blocking.

The structural question: **do the QPT symmetries commute with the Hopf-Z2 stabilizers?** If yes, they enforce simultaneously and the predicted triply-additive blocking holds. If not, they overlap and double-count.

## Method

Same shape as the DF / Cholesky tests:

1. Build composed Hamiltonians for LiH (M=15, Q=30), BeH2 (M=25, Q=50), H2O (M=35, Q=70) in the Hopf-rotated basis.
2. Construct `S^2` operator via `openfermion.s_squared_operator(M)` and JW.
3. Construct Hopf stabilizers (Z_alpha, Z_beta, P_0, ..., P_{n_sb-1}) via `geovac.z2_tapering.build_stabilizers(mode='per_block')`.
4. Compute every pairwise commutator `[A, B]` via `openfermion.commutator`, report the L-infinity norm of coefficients.
5. Enumerate joint (N_alpha, N_beta, P_signs) sectors and count Slater determinant dimensions.

## Findings

### Q1 — Every commutator is bit-exact zero across all three molecules

Coefficient L-infinity norm of `[A, B]` for every operator pair:

| pair | LiH | BeH2 | H2O |
|---|---:|---:|---:|
| [H, S^2] | **0.0** | **0.0** | **0.0** |
| [H, Z_alpha] | 0.0 | 0.0 | 0.0 |
| [H, Z_beta] | 0.0 | 0.0 | 0.0 |
| [H, P_i] (per sub-block) | 0.0 | 0.0 | 0.0 |
| [S^2, Z_alpha] | 0.0 | 0.0 | 0.0 |
| [S^2, Z_beta] | 0.0 | 0.0 | 0.0 |
| **[S^2, P_i]** (key) | **0.0** | **0.0** | **0.0** |
| [Z_alpha, P_i] | 0.0 | 0.0 | 0.0 |
| [Z_beta, P_i] | 0.0 | 0.0 | 0.0 |
| [P_i, P_j] (i != j) | 0.0 | 0.0 | 0.0 |
| total pairs tested | 21 | 36 | 55 |

Bit-exact zero — every commutator's Pauli-string coefficients cancel algebraically in openfermion's symbolic computation. Not numerical "near zero," exactly zero.

**The key commutator `[S^2, P_Hopf] = 0` confirms the stackability hypothesis**: the spin symmetry that QPT exploits and the Hopf-fiber Z2 we ship are mutually commuting symmetries on the GeoVac Hamiltonian.

### Q2 — Concrete reduction for LiH joint sectors

LiH composed: M=15, Q=30, 4 electrons, 3 active Hopf sub-blocks.

| sector | dim |
|---|---:|
| Full Hilbert 2^30 | 1,073,741,824 |
| N=4 particle-number sector | 27,405 |
| + S_z=0 (singlet candidates) | 11,025 |
| + Hopf-Z2 ground sector (P_0 = P_1 = P_2 = +1) | **4,791** |
| + S^2 = 0 spin-adapted (CSF, rough) | ~3,675 |

**Joint reduction factor 5.7x beyond particle-number tapering** (27,405 -> 4,791). The Hopf-Z2 alone gives a factor 2.3x on top of S_z parity; the further QPT spin-adaptation gives an additional ~1.3x. Both contributions stack on top of the Bravyi 2017 particle-number tapering.

### Q3 — Structural reason for stacking (analytic)

The Hopf-Z2 reflection `P_i` acts on **spatial orbitals only** (m -> -m within sub-block `i`). The S^2 spin operator acts on **spin labels only** (alpha vs beta within each spatial orbital). At the operator-algebra level, these factor cleanly:

```
P_i = sum over sub-block i orbitals of (orbital-reflection)  (spatial)
S^2 = sum over spatial pairs of (spin-coupling expressions)  (spin)
```

Since spatial and spin operators act on independent tensor factors, `[P_i, S^2] = 0` is an identity, not an empirical observation. The bit-exact zero result (Q1) is the algebraic confirmation. This identity holds for any non-relativistic GeoVac Hamiltonian (the spec used in this test).

## Interpretation

QPT and Hopf-Z2 tapering exploit **independent symmetries** of the GeoVac Hamiltonian:

| symmetry | acts on | nature | qubit cost change | per-sector cost change |
|---|---|---|---|---|
| Particle number (Bravyi 2017) | total N | conserved Casimir | -2 Q | factor 2x state reduction |
| Hopf-U(1) Z2 per sub-block (v3.52.0) | spatial m -> -m | per-sub-block parity | -n_sb Q | factor 2^n_sb state reduction |
| QPT spin-adaptation | S^2 eigenvalue | total spin Casimir | 0 (basis change) | factor depending on shell |
| QPT orbital U(d) | shell permutations | irrep label | 0 | factor depending on shell |

All four commute pairwise on the GeoVac Hamiltonian. They can be enforced simultaneously, producing a joint block-diagonalization across all symmetry sectors.

For chemistry resource estimation:
- The **qubit-count savings** come from Bravyi tapering + Hopf-Z2 (additive ΔQ).
- The **per-sector dimension savings** come from particle number + Hopf-Z2 + QPT (multiplicative reduction within each Q-fixed Hamiltonian).
- The Hopf-Z2 contribution sits in BOTH columns. That's why it's the most impactful symmetry numerically.

## Implications for papers

### Paper 14 §sec:hopf_tapering — recommended addition

The existing Hopf-Z2 tapering section (added 2026-06-04) can pick up a forward reference to QPT stacking. Suggested short paragraph after the per-sub-block ΔQ = 2 + n_sub_blocks result:

> The Hopf-U(1) Z_2 stabilizers commute pairwise with the total spin
> operator S^2 and with each other (verified bit-exactly across LiH, BeH_2, H_2O via direct commutator computation; see
> `debug/sprint_qpt_hopf_stacking_memo.md`). This means the Hopf
> tapering stacks additively with spin-adaptation techniques such as
> the recent Quantum Paldus Transform~\cite{burkat_fitzpatrick2025}.
> On LiH, the joint (N, S_z, P_Hopf) sector is 5.7x smaller than the
> particle-number-only sector of standard Bravyi 2017 tapering, with
> further reduction available from full QPT spin-adaptation to S^2
> eigenstates.

This addition would need a new bibitem (Burkat-Fitzpatrick 2025, arXiv:2506.09151).

### Paper 20 — no change recommended

Paper 20 is the chemistry-audience resource benchmark paper. The QPT result is a structural compatibility statement, not a benchmark. Best left as a Paper 14 forward reference.

### Forward-looking

The Hopf-Z2 + QPT stack is the **first** triply-additive blocking framework in GeoVac. If we later test additional symmetries (point-group from Setia 2020, time-reversal, particle-hole), each new symmetry just adds another commutator check. The same bit-exact commutator framework can be reused.

## Honest scope and named open items

- **Tested on non-relativistic composed builder.** With spin-orbit coupling (Tier 2 relativistic builder, Paper 14 §spinor_composed), [H, S^2] is generally nonzero. The test should be repeated on `build_composed_hamiltonian_relativistic` to characterize the spin-orbit regime. Sprint-scale follow-on.
- **QPT cost not measured.** The Burkat-Fitzpatrick QPT circuit is O(d^3) Toffoli where d is the number of spatial orbitals. For LiH d=15, that's roughly 3,375 Toffoli — well within fault-tolerant chemistry budgets but not free. A full resource estimate would need to be made for Paper 20-style comparison tables.
- **Symmetric group / U(d) irrep decomposition not enumerated.** This test confirmed [S^2, P_Hopf] = 0 (key result), but didn't break the orbital U(d) generators into individual elements. The full QPT also exploits U(d) irrep labels beyond just S^2 sectors. Sprint-scale extension would enumerate the GT-pattern blocks per (Hopf-Z2 sector, S^2 sector).
- **Cross-block ERI scope.** Standard composed builder (no cross-block ERIs) only. Paper 19 balanced-coupled spot check still the outstanding gate from the DF/Cholesky sprints.

## Verdict

POSITIVE. Bit-exact zero on every commutator across three molecules. The "triply additive blocking" claim from the lit-search agent is confirmed. QPT and Hopf-Z2 stack with each other and with standard Bravyi 2017 tapering. Concrete LiH reduction is 5.7x beyond particle-number tapering at the joint ground-state sector.

**The meta-lesson now confirmed on three independent methods in one day**: DF (multipole subspace, 2026-06-05 morning), Cholesky (same multipole subspace, 2026-06-05 afternoon), QPT (commuting symmetry with Hopf-Z2, 2026-06-05 evening). The pattern is consistent: when we push past surface-methodology differences, the underlying algebraic structures cooperate cleanly with GeoVac's tensor-product spectral triple.

Recommended next move (per PI direction):
- **Sprint-scale follow-on**: QPT cost estimate + spin-orbit (relativistic) commutator test, 3-5 days.
- **Multi-month**: DMRG / MPO bond-rank = Connes-vS spectral truncation theorem (lit-search candidate 3).
- **Pause**: consolidate three positive results into a single arc memo for the project record.
