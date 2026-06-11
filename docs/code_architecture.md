# Code architecture — full entry-point catalogue (extracted from CLAUDE.md S7, live document)

> Extracted verbatim from CLAUDE.md on 2026-06-10 (v3.110.0 state) during compaction round 2. This file is the canonical archive; the CLAUDE.md section now holds only the compact working form.

## S7 (update here, not in CLAUDE.md)

## 7. Code Architecture

### Key Entry Points

| Task | Module | Entry Point |
|:-----|:-------|:------------|
| Atomic lattice | `geovac/lattice.py` | `GeometricLattice(Z, max_n)` |
| Atomic Hamiltonian | `geovac/hamiltonian.py` | `GraphHamiltonian(lattice)` |
| Multi-electron FCI | `geovac/lattice_index.py` | `LatticeIndex(Z, n_electrons, max_n)` |
| Direct CI (large systems) | `geovac/direct_ci.py` | `DirectCISolver(lattice_index)` |
| Molecular FCI (LCAO) | `geovac/lattice_index.py` | `MolecularLatticeIndex(atoms, R)` |
| Prolate spheroidal (H2+) | `geovac/prolate_spheroidal.py` | Prolate spheroidal lattice for diatomics |
| Hyperspherical (He) | `geovac/hyperspherical.py` | Two-electron hyperspherical solver |
| Level 4 (H2) | `geovac/level4_multichannel.py` | Molecule-frame hyperspherical |
| Core screening | `geovac/core_screening.py` | `CoreScreening(Z).solve()` |
| Ab initio PK | `geovac/ab_initio_pk.py` | `AbInitioPK(core, n_core)` |
| Composed diatomic | `geovac/composed_diatomic.py` | `ComposedDiatomicSolver.LiH_ab_initio(l_max)` |
| Rho-collapse cache | `geovac/rho_collapse_cache.py` | `AngularCache`, `FastAdiabaticPES` |
| Quantum dynamics | `geovac/dynamics.py` | O(V) time evolution |
| Hopf bundle | `geovac/hopf_bundle.py` | S3 embedding, Hopf projection, fiber analysis |
| VQE benchmark | `geovac/vqe_benchmark.py` | `build_geovac_he()`, `run_vqe()`, `collect_static_metrics()` |
| Gaussian reference | `geovac/gaussian_reference.py` | `h2_sto3g()`, `he_sto3g()`, `he_cc_pvdz()`, `he_cc_pvtz()` |
| Measurement grouping | `geovac/measurement_grouping.py` | `qwc_groups()`, `analyze_measurement_cost()` |
| Trotter bounds | `geovac/trotter_bounds.py` | `pauli_1norm()`, `trotter_steps_required()` |
| Composed qubit | `geovac/composed_qubit.py` | `build_composed_lih()`, `build_composed_beh2()`, `build_composed_h2o()`, `build_h2_bond_pair()` |
| Algebraic angular solver | `geovac/algebraic_angular.py` | `AlgebraicAngularSolver(Z, l_max)` |
| Algebraic coupled-channel | `geovac/algebraic_coupled_channel.py` | `solve_hyperspherical_algebraic_coupled()` |
| Hyperradial solver | `geovac/hyperspherical_radial.py` | `solve_radial()`, `solve_radial_spectral()` |
| Spectral angular (Level 4) | `geovac/level4_spectral_angular.py` | `SpectralAngularSolver`, `angular_method='spectral'` |
| Perturbation series (Level 3) | `geovac/algebraic_angular.py` | `perturbation_series_mu()`, `pade_approximant()` |
| Algebraic curve (Level 3) | `geovac/algebraic_curve.py` | `characteristic_polynomial()`, `algebraic_veff_matrix_lmax0()` |
| Algebraic curve (Level 4) | `geovac/algebraic_curve_level4.py` | `extract_level4_matrices()`, `probe_rho_dependence()` |
| Algebraic Z_eff | `geovac/algebraic_zeff.py` | `LaguerreZeffExpansion`, `fit_density_laguerre()` |
| Algebraic Slater integrals | `geovac/algebraic_slater.py` | `slater_fk_integral_algebraic()`, `slater_f0_algebraic()` |
| Cusp factor (Level 4) | `geovac/cusp_factor.py` | `CuspFactorSolver` (negative result — alpha-only factor insufficient) |
| Cusp graph analysis | `geovac/cusp_graph.py` | `analyze_cusp_graph()` (theory: 1/r₁₂ cannot be absorbed into S⁵ graph) |
| Cusp correction | `geovac/cusp_correction.py` | `he_cusp_correction()`, `h2_cusp_correction()`, `cusp_pes_correction()` |
| Cusp angular basis | `geovac/cusp_angular_basis.py` | `verify_so6_casimir()`, `extract_channel_coefficients()` (Track Y negative result) |
| Geometric elevation | `geovac/geometric_elevation.py` | Blow-up, Lie algebra, S³×S³ analysis (Track Z negative result) |
| PK R_ref derivation | `geovac/pk_rref.py` | `compute_rref_candidates()`, `select_rref()` (Track AD negative result) |
| N-electron solver | `geovac/n_electron_solver.py` | `scan_pes_4e_lih_multichannel()` |
| N-electron scoping | `geovac/n_electron_scope.py` | `four_electron_channel_count_molecular()` |
| N-electron spectral | `geovac/n_electron_spectral.py` | Spectral compression analysis |
| N-electron 2D solver | `geovac/n_electron_2d.py` | 2D variational solver for 4-electron system (Track AR) |
| N-electron qubit encoding | `geovac/n_electron_qubit.py` | Full N-electron quantum encoding comparison (Track AS) |
| Ecosystem export | `geovac/ecosystem_export.py` | `GeoVacHamiltonian`, `hamiltonian()`, `.to_openfermion()`, `.to_qiskit()`, `.to_pennylane()` |
| PK partitioning | `geovac/pk_partitioning.py` | `pk_classical_energy()`, `reconstruct_1rdm_from_statevector()`, `validate_pk_partitioning()` |
| Atomic classifier | `geovac/atomic_classifier.py` | `classify_atom(Z)`, `AtomClassification`, `pk_params_z2_scaled()` |
| Molecular spec | `geovac/molecular_spec.py` | `MolecularSpec`, `OrbitalBlock` |
| General composed builder | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec)`, `lih_spec()`, `beh2_spec()`, `h2o_spec()`, `hf_spec()`, `nh3_spec()`, `ch4_spec()` |
| TC integrals | `geovac/tc_integrals.py` | `compute_tc_integrals_block()`, `build_tc_composed_hamiltonian()` |
| Coupled composition | `geovac/coupled_composition.py` | `build_coupled_hamiltonian()`, `coupled_fci_energy()`, `run_coupled_scoping()` (Track CB, negative result) |
| Shibuya-Wulfman V_ne | `geovac/shibuya_wulfman.py` | `compute_cross_center_vne()`, `compute_cross_center_vne_element()` |
| Balanced coupled builder | `geovac/balanced_coupled.py` | `build_balanced_hamiltonian(spec, nuclei)` |
| Frozen core (Ne-like) | `geovac/neon_core.py` | `FrozenCore(Z)` |
| NaH spec | `geovac/composed_qubit.py` | `nah_spec()` |
| MgH₂ spec | `geovac/composed_qubit.py` | `mgh2_spec()` |
| HCl spec | `geovac/composed_qubit.py` | `hcl_spec()` |
| H₂S spec | `geovac/composed_qubit.py` | `h2s_spec()` |
| PH₃ spec | `geovac/composed_qubit.py` | `ph3_spec()` |
| SiH₄ spec | `geovac/composed_qubit.py` | `sih4_spec()` |
| KH spec | `geovac/composed_qubit.py` | `kh_spec()` |
| CaH₂ spec | `geovac/composed_qubit.py` | `cah2_spec()` |
| GeH₄ spec | `geovac/composed_qubit.py` | `geh4_spec()` |
| AsH₃ spec | `geovac/composed_qubit.py` | `ash3_spec()` |
| H₂Se spec | `geovac/composed_qubit.py` | `h2se_spec()` |
| HBr spec | `geovac/composed_qubit.py` | `hbr_spec()` |
| SrH spec ([Kr] frozen core, Sprint 3) | `geovac/molecular_spec.py` | `srh_spec()`, `srh_spec_relativistic()` |
| BaH spec ([Xe] frozen core, Sprint 3) | `geovac/molecular_spec.py` | `bah_spec()`, `bah_spec_relativistic()` |
| Alkaline-earth monohydride builder | `geovac/molecular_spec.py` | `_alkaline_earth_monohydride_spec(Z, name, ...)` |
| LiF spec (multi-center) | `geovac/composed_qubit.py` | `lif_spec()` |
| CO spec (multi-center) | `geovac/composed_qubit.py` | `co_spec()` |
| N₂ spec (multi-center) | `geovac/composed_qubit.py` | `n2_spec()` |
| F₂ spec (multi-center) | `geovac/composed_qubit.py` | `f2_spec()` |
| NaCl spec (multi-center) | `geovac/composed_qubit.py` | `nacl_spec()` |
| CH₂O spec (multi-center) | `geovac/composed_qubit.py` | `ch2o_spec()` |
| C₂H₂ spec (multi-center) | `geovac/composed_qubit.py` | `c2h2_spec()` |
| C₂H₆ spec (multi-center) | `geovac/composed_qubit.py` | `c2h6_spec()` |
| Multi-center nuclear repulsion | `geovac/composed_qubit.py` | `_compute_multi_center_nuclear_repulsion()` |
| TM hydride spec (Z=21-30) | `geovac/molecular_spec.py` | `transition_metal_hydride_spec(Z)`, `sch_spec()`, ..., `znh_spec()` |
| General composed builder | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec)` |
| General-l Wigner D rotation | `geovac/shibuya_wulfman.py` | `_build_rotation_matrix_real_sh(l, theta, phi)`, `_build_block_rotation_matrix()` |
| Nested hyperspherical | `geovac/nested_hyperspherical.py` | `build_nested_hamiltonian()`, H-set coupled basis |
| Level 3 variational 2D | `geovac/level3_variational.py` | `solve_he_variational_2d()`, `solve_he_precision()` |
| Algebraic Casimir CI | `geovac/casimir_ci.py` | `build_fci_matrix()`, `solve_self_consistent()`, `solve_variational()`, `convergence_table()` |
| Graph-native CI | `geovac/casimir_ci.py` | `build_graph_native_fci()`, `build_graph_consistent_fci()`, `_build_graph_h1()` |
| Hypergeometric Slater integrals | `geovac/hypergeometric_slater.py` | `compute_rk_float()`, `compute_rk_exact()`, exact Fraction-arithmetic R^k evaluator for arbitrary n_max |
| DUCC downfolding | `geovac/downfolding.py` | Exact (2J-K) core potential, PK divergence root-cause analysis |
| Physical constants | scattered next to their modules | `C_LIGHT` in `geovac/dirac_hamiltonian.py`; `ALPHA` in `geovac/two_loop_self_energy.py`; `KAPPA_SCALAR = Rational(-1,16)` in `geovac/graph_qed_propagator.py`. No central `geovac/constants.py` module. `-1/16` per CLAUDE.md §8 may be used directly. |
| Breit retarded integrals | `geovac/breit_integrals.py` | `compute_radial()`, exact Fraction/sympy r_<^k / r_>^{k+3} kernel |
| Breit composed pipeline | `geovac/composed_qubit.py` | `build_composed_hamiltonian(spec, include_breit=True)` passthrough to relativistic builder |
| QED vacuum polarization | `geovac/qed_vacuum_polarization.py` | `seeley_dewitt_coefficients_s3()`, `vacuum_polarization_coefficient()`, `beta_function_qed()`, `spectral_zeta_massive()`, `spectral_zeta_derivative_at_zero()`, `classify_transcendental_content()` |
| QED vertex coupling | `geovac/qed_vertex.py` | `two_loop_odd_even_split()`, `two_loop_sunset_vertex_restricted()`, `two_loop_transcendental_classification()`, `vertex_allowed_triples()`, `reduced_coupling_squared()` |
| QED flat-space limit | `geovac/qed_flat_limit.py` | `weyl_density_check()`, `two_loop_sunset_R_dependent()`, `flat_space_limit_scaling()`, `two_loop_vacuum_energy_R_dependent()`, `extract_zeta3_coefficient()`, `verify_weyl_zeta3()` |
| QED self-energy | `geovac/qed_self_energy.py` | `self_energy_spectral()`, `self_energy_table()`, `self_energy_convergence()`, `self_energy_transcendental_class()`, `vertex_correction_spectral()`, `schwinger_convergence()` |
| QED three-loop | `geovac/qed_three_loop.py` | `three_loop_sunset_s3()`, `three_loop_unrestricted()`, `three_loop_vertex_restricted()`, `three_loop_cg_weighted()`, `three_loop_factorized()`, `three_loop_factorized_convergence()`, `decompose_three_loop_mzv()`, `three_loop_euler_maclaurin_tail()` |
| Graph-native Fock graph/Hodge | `geovac/fock_graph_hodge.py` | `build_fock_graph()`, `compute_hodge_decomposition()`, `verify_hodge_identity()` |
| Graph-native electron propagator | `geovac/graph_qed_propagator.py` | `DiracGraphOperator`, `electron_propagator()`, `neumann_electron_propagator()` |
| Graph-native photon propagator | `geovac/graph_qed_photon.py` | `build_fock_graph()`, `compute_photon_propagator()`, `compute_vacuum_polarization()` |
| Graph-native vertex tensor | `geovac/graph_qed_vertex.py` | `build_projection_matrix()`, `build_vertex_tensor()`, `vertex_tensor_to_matrices()` |
| Graph-native self-energy/vertex | `geovac/graph_qed_self_energy.py` | `compute_self_energy()`, `compute_vertex_correction()`, `extract_anomalous_moment()`, `self_energy_structural_zero_check()`, `self_energy_pi_free_certificate()` |
| Graph-native continuum bridge | `geovac/graph_qed_continuum_bridge.py` | `compute_projection_exchange_constant()`, `compute_continuum_bridge()`, `bridge_self_energy()` |
| Vector-photon QED | `geovac/vector_qed.py` | `build_electron_states()`, `build_photon_modes()`, `vertex_coupling()`, `compute_self_energy()`, `check_selection_rules()`, `check_selection_rules_dirac()` |
| Connes-vS truncated operator system | `geovac/operator_system.py` | `TruncatedOperatorSystem(n_max)`, `.contains(M)`, `.identity_in_O()`, `.is_star_closed()`, `.compute_propagation_number()` |
| Circulant-S³ falsification comparator | `geovac/circulant_s3.py` | `CirculantS3Truncation(n_points)`, `.compute_propagation_number()`, `circulant_for_geovac(n_max)`, `compare_to_geovac(n_max)` |
| Connes distance SDP | `geovac/connes_distance.py` | `compute_connes_distance(op_sys, v, w, dirac_mode='offdiag')`, `compute_distance_matrix(op_sys)`, gauge-fixed cvxpy SDP with operator-norm LMI |

### Solver Methods

| Method | Access | Use Case |
|:-------|:-------|:---------|
| Mean-field | `LatticeIndex(method='mean_field')` | Quick atomic energies |
| Full CI (matrix) | `LatticeIndex(method='full_ci')` | Small systems (N_SD < 5000) |
| Full CI (direct) | `DirectCISolver` or `fci_method='direct'` | Large systems (N_SD >= 5000) |
| Auto | `fci_method='auto'` | Switches at N_SD = 5000 |
| Frozen core | `FrozenCoreLatticeIndex` | Active-space CI for core-valence |
| Locked shell | `LockedShellMolecule` | Extreme SD reduction for heavy atoms |

---

