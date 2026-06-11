# Papers index — status map

One row per paper. This is the newcomer's map of the corpus; the agent-facing loading
guide lives in CLAUDE.md §6, and per-claim verification tiers live in
[`docs/claims_register.md`](../docs/claims_register.md).

**If you read only five things:**
1. `synthesis/geovac_field_guide.tex` — what GeoVac is, as a story (10 pp)
2. [`docs/claims_register.md`](../docs/claims_register.md) — every headline claim + its verification tier + falsifier
3. `group3_foundations/Paper_7_Dimensionless_Vacuum.tex` — the S³ equivalence (the theoretical foundation)
4. `group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` — the discrete→continuum convergence theorem (unconditional)
5. `group4_quantum_computing/paper_14_qubit_encoding.tex` — the headline computational result

## Status legend

| Status | Meaning |
|:-------|:--------|
| **KEYSTONE** | Load-bearing identity paper — the framework stands on it |
| **ACTIVE** | Current and accurate; on-topic reference |
| **OBSERVATION** | Honest empirical/structural observation; interpretive derivation explicitly NOT claimed |
| **GUARDRAIL** | Proven negative result that constrains future approaches (read before working in its domain) |
| **PARTIAL** | Carries an in-place Status note descoping part of its claims; surviving content stated in the paper |
| **DESCOPED** | Main theorem withdrawn/corrected in place (2026-06-09 hardening pass); paper documents the corrected state |
| **DRAFT** | First draft, not yet hardened |
| **HISTORICAL** | Archived — superseded or orphaned; content remains a valid record |

## Synthesis (`synthesis/`)

| File | Status | One-liner |
|:-----|:------:|:----------|
| `geovac_field_guide.tex` | **KEYSTONE** | The "what is GeoVac" document — packing puzzle → spectral triple → periods → forced/free seam |
| `group1_operator_algebras_synthesis.tex` | ACTIVE | Narrative of the math.OA arc (Papers 29–49) |
| `group3_foundations_synthesis.tex` | ACTIVE | Narrative of the foundations arc (Papers 0–31) |

## Group 3 — Foundations (`group3_foundations/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 0 `Paper_0_Geometric_Packing.tex` | **KEYSTONE** | The packing axiom; K = −1/16; nodes are quantum numbers |
| 1 `paper_1_spectrum.tex` | **KEYSTONE** | Spectral graph theory; O(N) eigenvalue methods |
| 7 `Paper_7_Dimensionless_Vacuum.tex` | **KEYSTONE** | S³ conformal equivalence, 18 symbolic proofs — most-cited paper in the corpus |
| 18 `paper_18_exchange_constants.tex` | ACTIVE | Transcendental taxonomy; master Mellin engine M1/M2/M3 |
| 22 `paper_22_angular_sparsity.tex` | **KEYSTONE** | Potential-independent angular sparsity theorem (underwrites the Pauli scaling) |
| 24 `paper_24_bargmann_segal.tex` | **KEYSTONE** | π-free HO lattice on S⁵; HO rigidity; Coulomb/HO asymmetry (4 layers) |
| 31 `paper_31_universal_coulomb_partition.tex` | **KEYSTONE** | Universal vs Coulomb-specific (A/D) split of the spectral triple |
| 54 `paper_54_tensor_product_two_body.tex` | ACTIVE | Two-body selection rules from the tensor-product triple; radial coupling NOT forced |
| 55 `paper_55_periods_of_geovac.tex` | ACTIVE | Periods of the Mellin sectors: pure-Tate / level-4 cyclotomic mixed-Tate |
| 56 `paper_56_tannakian_substrate.tex` | ACTIVE | Tannakian/cosmic-Galois reconstruction at finite cutoff (5,864 bit-exact residuals) |
| 57 `paper_57_forced_free_seam.tex` | ACTIVE | 60-entry forced/free catalogue; P5 packing-reachability discriminator (98.3%) |

## Group 1 — Operator algebras / NCG (`group1_operator_algebras/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 29 `paper_29_ramanujan_hopf.tex` | ACTIVE | Hopf graphs are Ramanujan; graph-RH is a finite-size statement |
| 32 `paper_32_spectral_triple.tex` | **KEYSTONE** | The GeoVac spectral triple; Connes axiom audit; §VIII case-exhaustion + non-selection theorems |
| 38 `paper_38_su2_propinquity_convergence.tex` | **KEYSTONE** | **Unconditional** state-space GH convergence at rate (4/π)·log n/n (2026-06-10) |
| 39 `paper_39_tensor_propinquity_convergence.tex` | ACTIVE | Tensor-product convergence (two focal lengths) |
| 40 `paper_40_unified_propinquity_convergence.tex` | ACTIVE | 4/π rate universal across compact Lie groups |
| 42 `paper_42_modular_hamiltonian_four_witness.tex` | ACTIVE | Four-witness Wick-rotation theorem at finite cutoff (Riemannian) |
| 43 `paper_43_lorentzian_extension.tex` | ACTIVE | Krein (3,1) extension at finite cutoff; Pythagorean orthogonality with 1/π² M1 signature |
| 44 `paper_44_lorentzian_operator_system.tex` | ACTIVE | Lorentzian operator-system substrate; prop = 2 |
| 45 `paper_45_lorentzian_propinquity.tex` | **DESCOPED** | Main theorem withdrawn (K⁺ seminorm ≡ 0, annihilation theorem); spatial statement unconditional; product-carrier S³×S¹ convergence rebuilt in the action-seminorm framework (2026-06-10, signature-agnostic) |
| 46 `paper_46_strong_form_lorentzian_propinquity.tex` | **DESCOPED** | Strong-form claims pending product repair; Lemma 3.2 = the degeneracy diagnosis |
| 47 `paper_47_two_rate_hybrid_convergence.tex` | PARTIAL | Norm-resolvent arrow + three-carrier identification stand; propinquity arrow descoped |
| 48 `paper_48_krein_ms_bridge.tex` | PARTIAL | Bridge design conditional on repair; T3/T6 descoped |
| 49 `paper_49_oslpls_strong_form_bridge.tex` | PARTIAL | Λ inheritance descoped; cocycle-deficit / TICI algebra survives |
| 50 `paper_50_cft3_partition_function.tex` | ACTIVE | Bit-exact F-theorem match on S³ and S⁵ (Klebanov–Pufu–Safdi) |
| 52 `paper_52_category_iii_correspondence.tex` | DRAFT | Positioning: spectral-triple discretization as Category III (non-holographic) |
| 53 `paper_53_disk_propinquity.tex` | DRAFT | First manifold-with-boundary carrier; plane Bochner–Riesz reconstruction |

## Group 2 — Quantum chemistry (`group2_quantum_chemistry/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 8–9 `Paper_8_Bond_Sphere_Sturmian.tex` | **GUARDRAIL** | Sturmian structural theorem: single-center molecular encodings are R-independent (proven dead end) |
| 11 `paper_11_prolate_spheroidal.tex` | ACTIVE | H₂⁺ at 0.0002% via spectral Laguerre |
| 12 `paper_12_algebraic_vee.tex` | ACTIVE | Algebraic V_ee (Neumann expansion) |
| 13 `paper_13_hyperspherical.tex` | ACTIVE | He at 0.019%; graph-native CI at 0.20% with zero parameters |
| 15 `paper_15_level4_geometry.tex` | ACTIVE | H₂ at 96.0% D_e (molecule-frame hyperspherical) |
| 17 `paper_17_composed_geometries.tex` | ACTIVE | Composed geometry: LiH R_eq 5.3%; the production molecular architecture |
| 19 `paper_19_coupled_composition.tex` | ACTIVE | Balanced coupled builder; PK-free cross-center V_ne |
| `paper_fci_atoms.tex` | ACTIVE | He/Li/Be FCI benchmarks |
| `paper_fci_molecules.tex` | **GUARDRAIL** | LCAO graph-concatenation fails (R-independent kinetic energy) — why natural geometry was necessary |

## Group 4 — Quantum computing (`group4_quantum_computing/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 14 `paper_14_qubit_encoding.tex` | **KEYSTONE** | The headline: O(Q^2.5) Pauli scaling, 51×–1,712× vs Gaussian baselines |
| 16 `paper_16_periodicity.tex` | **KEYSTONE** | Chemical periodicity from S_N representation theory |
| 20 `paper_20_resource_benchmarks.tex` | ACTIVE | Resource benchmarks, FCI PES, pip install |
| 23 `paper_23_nuclear_shell.tex` | **KEYSTONE** | Nuclear qubit Hamiltonians; Fock rigidity theorem (cross-group hub) |

## Group 5 — QED / gauge / SM (`group5_qed_gauge/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 2 `paper_2_alpha.tex` | **OBSERVATION** | K = π(B+F−Δ) matches 1/α at 8.8×10⁻⁸; combination rule **conjectural** (12 mechanisms eliminated) |
| 25 `paper_25_hopf_gauge_structure.tex` | OBSERVATION | Hopf graph as discrete lattice-gauge structure |
| 28 `paper_28_qed_s3.tex` | ACTIVE | QED on S³: 5 theorems; S_min irreducible at 200 digits |
| 30 `paper_30_su2_wilson.tex` | ACTIVE | SU(2) Wilson gauge on S³ = SU(2); Paper 25 is its maximal torus |
| 33 `paper_33_qed_selection_rules.tex` | OBSERVATION | 1+6+1 partition of QED selection rules |
| 36 `paper_36_bound_state_qed.tex` | ACTIVE | Hydrogen Lamb shift at −0.534% at one loop, no fits; two-loop wall named |
| 41 `paper_41_rule_b_wilson_u1.tex` | ACTIVE | Wilson U(1) on the Dirac Rule-B graph; seven-witness compatibility |
| 51 `paper_51_gravity_arc.tex` | ACTIVE | Gravity from the spectral action; two-term exactness; J-blindness theorem |

## Group 6 — Precision / observations (`group6_precision_observations/`)

| Paper | Status | One-liner |
|:------|:------:|:----------|
| 26 `paper_26_entanglement.tex` | ACTIVE | Energy–entanglement decoupling; basis-intrinsic sparsity |
| 27 `paper_27_entropy_projection.tex` | **KEYSTONE** | Entropy as projection artifact; HO zero-entropy rigidity |
| 34 `paper_34_projection_taxonomy.tex` | ACTIVE | **Living catalogue**: 28 named projections — where physics enters the graph |
| 35 `paper_35_time_as_projection.tex` | ACTIVE | π enters exactly at temporal compactification (Matsubara); falsifiable prediction |

## Archive (`archive/`) — HISTORICAL, load on request

| Paper | Note |
|:------|:-----|
| 3 `paper_3_holography.tex` | Holographic machinery partially retracted |
| 4 `Paper_4_Universality.tex` | Proton-radius claims overtaken by measurement; kernel absorbed into Papers 34/35 |
| 5 `Paper_5_Geometric_Vacuum.tex` | Early synthesis |
| 6 `Paper_6_Quantum_Dynamics.tex` | Valid tool paper; citation-orphan in the corpus |
| 10 `paper_10_nuclear_lattice.tex` | Early nuclear draft |
| `paper_18_exchange_constants_v1.tex` | Superseded by current Paper 18 |
| 21 `paper_21_geometric_vacuum_synthesis.tex` | Superseded by the two group syntheses |
