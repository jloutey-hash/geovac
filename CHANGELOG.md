# Changelog

All notable changes to GeoVac will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.11] - 2026-03-29

### Paper 18: Spectral-Geometric Exchange Constants

New observations paper cataloging the transcendental constants that appear when discrete algebraic structures are projected onto continuous coordinate systems across the GeoVac hierarchy. Identifies them as instances of Weyl's law and the Selberg trace formula.

#### Paper content
- **Section II.B expanded with Track J evidence:** Full partial-fraction decomposition of centrifugal singularity, associated Laguerre basis L_n^{|m|}(x), algebraic M_{-1} via DLMF 18.9.13, Stieltjes recurrence seeded by single transcendental e^a·E₁(a). Key observation: continuous math compressed from function (quadrature integrand) to constant (Stieltjes seed).
- **Basis adaptation observation:** Identifying the exchange constant and absorbing it into the basis improves convergence (associated Laguerre stable by N=10 vs N≈20 for ordinary). The transcendental tells you where the curvature mismatch lives.
- **π as founding example (Section III):** π is the primordial exchange constant (S¹ eigenvalue count ↔ circumference); all GeoVac constants are descendants.
- **Computational cost observation (Section VI):** Exchange constants are conceptual tolls, not computational bottlenecks. κ: one multiply; seed: one function call; μ(R): millisecond diagonalization.
- **α transcendence strengthened (Section V):** Two formal observations — structural match of K = π(B + F − Δ) to exchange constant pattern; α likely transcendental as root of cubic with transcendental coefficients.

#### MD integration
- Paper 18 added to CLAUDE.md: observations inventory, reading guide (item 11), topic-to-paper lookup (4 entries)
- Paper 18 added to README.md Paper Series table

**Files changed:**
- `papers/observations/paper_18_exchange_constants.tex` — new paper
- `CLAUDE.md` — Paper 18 registration, version bump
- `README.md` — Paper 18 registration, version bump
- `CHANGELOG.md` — this entry

## [2.0.10] - 2026-03-29

### Algebraic Laguerre Matrix Elements, Spectral Level 4 Hyperradial, & Associated Laguerre m≠0

#### Track H — Level 3 Algebraic Laguerre Matrix Elements (COMPLETE)

Replaced Gauss-Laguerre quadrature with algebraic three-term Laguerre recurrence for the hyperradial overlap (S) and kinetic (K) matrices in the Level 3 spectral solver.

- **Overlap S:** Pentadiagonal M2 moment matrix from three-term recurrence x·Lₙ = -(n+1)Lₙ₊₁ + (2n+1)Lₙ - nLₙ₋₁. S = M2/(8α³). Machine-precision agreement with quadrature (< 1e-14 relative error).
- **Kinetic K:** Derivative kernel Bₙ = -n/2 Lₙ₋₁ + 1/2 Lₙ + (n+1)/2 Lₙ₊₁ (tridiagonal expansion). K = BᵀB/(4α). Machine-precision agreement (< 1e-14 relative error).
- **11× matrix build speedup:** 52 μs (algebraic) vs 585 μs (quadrature)
- **V_eff stays quadrature:** Effective potential V_eff(R) = μ(R)/R² + 15/(8R²) is transcendental (Track G proven). Quadrature required — this is the algebraic boundary for Level 3.
- **Energy agreement:** < 1e-14 Ha vs quadrature (end-to-end validation)
- **Coupled-channel ceiling unchanged:** 0.15%–0.25% at l_max=3 (same physics, faster build)
- `matrix_method='algebraic'` parameter in spectral radial solver
- 11 new tests, 41/41 passing

**Files changed:**
- `geovac/hyperspherical_radial.py` — `_build_laguerre_matrices_algebraic()`, `matrix_method` parameter
- `geovac/algebraic_coupled_channel.py` — `matrix_method` parameter wiring
- `tests/test_hyperspherical_he.py` — 11 new tests

#### Track I — Spectral Laguerre for Level 4 Hyperradial (COMPLETE)

Applied spectral Laguerre pattern to the Level 4 molecule-frame hyperspherical solver, replacing FD R_e grid for all three solution pathways (adiabatic, coupled-channel, 2D variational).

- **16× dimension reduction (adiabatic):** N_FD=400 → n_basis=25
- **5× dimension reduction (2D):** N_FD×N_ang=16,000 → n_basis×N_ang=3,200
- **Spectral-FD consistency:** < 0.0005 Ha (adiabatic), < 0.002 Ha (2D)
- **D_e% preserved:** 90.0% for both spectral and FD at l_max=2 (H₂)
- **No wall time speedup — structural finding:** Angular sweep (building H_ang at ~130 ρ-values) dominates 99% of Level 4 cost. FD radial uses `eigh_tridiagonal` (O(N), sub-ms) — already optimally fast, unlike Level 3's FD which used N_R=3000 with `eigsh`. Spectral radial substitution does not improve total runtime.
- **Spectral value:** accuracy (< 0.0005 Ha with 16× fewer points), memory (400 → 25 grid), parameterization (just n_basis and α), extensibility (same basis for adiabatic, DBOC, and 2D)
- **n_basis=20 optimal:** Convergence plateau at n_basis=15-20. Mild overlap matrix conditioning degradation at n_basis≥25 (~5e-5 Ha shift).
- `radial_method='spectral'` parameter in Level 4 solver
- 13 new tests, 58/58 passing (45 existing + 13 new)

**Files changed:**
- `geovac/level4_multichannel.py` — `radial_method`/`n_basis_radial`/`alpha_radial` parameter wiring, spectral radial solve
- `tests/test_spectral_level4.py` — new (13 tests)

#### Track J — Level 2 Algebraic m≠0 Associated Laguerre (COMPLETE)

Eliminated the last remaining quadrature in the Level 2 prolate spheroidal radial solver. The centrifugal singularity m²/(ξ²−1) → m²/(x(x+2β)) for m≠0 states (π, δ, ...) was previously handled by Gauss-Laguerre quadrature. The associated Laguerre basis L_n^{|m|}(x) with weight x^|m|·e^{-x} absorbs the singularity, and the partial-fraction decomposition yields two algebraic components:

- **Lowered moment M_{-1}^(α):** Computed from DLMF 18.9.13 summation identity L_n^(α)(x) = Σ L_k^(α−1)(x). Cumulative sum of h_p^(α−1) = Γ(p+α)/p!. Dense, positive definite.
- **Stieltjes integral J^(α)(a):** Three-term recurrence seeded by single transcendental value S_0^(0)(a) = e^a·E₁(a). Miller's backward recurrence for numerical stability. Hybrid forward/backward assembly for matrix rows.
- **Associated moment matrices M₀, M₁, M₂:** Generalization of ordinary Laguerre moments to weight x^α·e^{-x}. Band structure preserved (diagonal, tridiagonal, pentadiagonal).
- **Associated kinetic matrix:** Derivative kernel G₀ (bidiagonal) + G₁ (diagonal) with weighted moment matrices P_k = M_k + 2β·M_{k-1}.

**Results:**
- Associated basis converges faster than ordinary for |m|=1: stable by N=10 (vs non-monotonic at N=40 for ordinary)
- |m|=2 (δ states): algebraic-quadrature agreement ~5e-9 at N=20
- PES scan (12 R-points, m=1): max algebraic-quadrature discrepancy 6e-4 (basis difference, not error); R_eq matches
- Algebraic method is MORE ACCURATE than FD: converges below all FD grid values
- `matrix_method='algebraic'` now works for ALL m values (auto-detected, no user change needed)
- 119 tests passing: 77 associated Laguerre moment/stability tests + 42 kinetic/eigenvalue/PES tests

**Qualitative distinction:** The m=0 algebraic construction is fully rational. For m≠0, the single transcendental seed e^a·E₁(a) is irreducible — the Stieltjes integral inherently involves an exponential integral. This is the algebraic boundary for Level 2.

**Files changed:**
- `geovac/prolate_spheroidal_lattice.py` — `_associated_laguerre_moment_matrices()`, `_lowered_moment_matrix()`, `_stieltjes_matrix()`, `_build_laguerre_matrices_algebraic_associated()`
- `tests/test_associated_laguerre.py` — new (77 tests)
- `tests/test_associated_kinetic.py` — new (42 tests)
- `papers/core/paper_11_prolate_spheroidal.tex` — Sec V.D expanded: associated Laguerre subsection, convergence table

## [2.0.9] - 2026-03-28

### Five-lane sprint (2 rounds, 3 parallel PMs)

#### Round 1

##### Lane 1 / Track C: Spectral PES Production Wiring (COMPLETE)

Wired spectral Laguerre radial solver into production PES scan and wavefunction reconstruction for H₂⁺.

- **`scan_h2plus_pes(radial_method='spectral')`:** 287× speedup (0.53s vs 153s, 41 R-points), 5000× E_min accuracy (0.0005% vs 1.01%), 160× D_e accuracy (0.04% vs 5.9%)
- **`solve_with_wavefunction(radial_method='spectral')`:** Spectral wavefunction reconstruction via Galerkin eigenvector on 500-point grid
- R_eq: 0.38% (slightly worse than FD's 0.15% due to coarse 41-point grid + polynomial fit — underlying energies are far more precise)
- Default behavior unchanged (`radial_method='fd'`)
- 20 new tests (test_spectral_pes_wiring.py), 37 total, all passing

**Files changed:**
- `geovac/prolate_spheroidal_lattice.py` — `radial_method`/`n_basis` params in `scan_h2plus_pes()`, spectral wavefunction reconstruction
- `tests/test_spectral_pes_wiring.py` — new (20 tests)

##### Lane 2 / Track F: Spectral Laguerre for Level 3 Hyperradial (COMPLETE)

Applied spectral Laguerre pattern to the Level 3 hyperradial solver, replacing 3000-point FD grid.

- **Basis:** φₙ(R) = (R-R_min)·exp(-α(R-R_min))·Lₙ(2α(R-R_min)), Gauss-Laguerre quadrature
- **120× dimension reduction:** 3000 FD grid → 25 spectral basis functions (converged at n_basis=15)
- **95× coupled-channel speedup:** 561 ms → 5.9 ms (3-channel, l_max=0)
- **Spectral-FD consistency:** < 0.00003 Ha across all l_max values tested (0-3)
- **Coupled-channel physics preserved:** l_max=3 error 0.221% (spectral) vs 0.220% (FD) — same ceiling
- Alpha-insensitive: energy spread < 0.001 Ha over α=[0.5, 3.0]
- Default behavior unchanged (`radial_method='fd'`)
- 10 new tests (test_hyperspherical_he.py), 131 total, all passing

**Files changed:**
- `geovac/hyperspherical_radial.py` — `solve_radial_spectral()`, `solve_coupled_radial_spectral()`, `_build_laguerre_matrices_dirichlet()`
- `geovac/algebraic_coupled_channel.py` — `radial_method`/`n_basis_radial`/`alpha_radial` param wiring
- `tests/test_hyperspherical_he.py` — 10 new tests

##### Lane 3 / Track A: 2D Solver Integration into Composition Pipeline (COMPLETE)

Integrated variational 2D solver (Paper 15) into composed diatomic pipeline for LiH.

- **Clean integration:** 2D solver already accepts PK potentials through `build_angular_hamiltonian()` pathway. No interface gap.
- **4× drift reduction:** l_max divergence rate +0.400 → +0.100 bohr/l_max (adiabatic → 2D)
- **75% of divergence is adiabatic** (as diagnosed in v2.0.6); **25% residual is non-adiabatic** — likely PK-related or intrinsic to composed-geometry separation. This is new information.
- **Adiabatic wins at l_max=2** (2.8% vs 6.1%) due to error cancellation; **2D wins at l_max=3** (9.5% vs 16.1%)
- Wall time comparable (2D sometimes faster — avoids 130-point R_e angular sweep)
- Default behavior unchanged (`level4_method='adiabatic'`)
- 77/77 tests passing, zero regressions

**Files changed:**
- `geovac/composed_diatomic.py` — `level4_method` parameter (`'adiabatic'` default, `'variational_2d'`)

#### Round 2

##### Lane A / Track C+: Algebraic Laguerre Matrix Elements (COMPLETE)

Replaced Gauss-Laguerre quadrature with algebraic three-term recurrence for all Level 2 radial matrix elements (m=0).

- **All matrix elements algebraic** for m=0 (σ states): overlap S, kinetic K, potential V computed via Laguerre moments M_k[i,j] = ∫x^k·Lᵢ·Lⱼ·e^{-x}dx using three-term recurrence x·Lₙ = -(n+1)Lₙ₊₁ + (2n+1)Lₙ - nLₙ₋₁
- **Machine-precision agreement:** overlap max diff 4.4e-15, Hamiltonian max diff 4.1e-12, energy diff ~1e-14 Ha
- **1.6× speedup** over quadrature (eliminates GL root computation)
- **m≠0 not algebraicizable:** m²/(ξ²-1) produces 1/x singularity in Laguerre space with no finite moment expansion. Falls back to quadrature automatically. Documented and tested.
- Default behavior unchanged (`matrix_method='quadrature'`)
- 19 new tests, 56 total, all passing

**Files changed:**
- `geovac/prolate_spheroidal_lattice.py` — `_laguerre_moment_matrices()`, `_build_laguerre_matrices_algebraic()`, `matrix_method` parameter
- `tests/test_spectral_radial.py` — 19 new tests

##### Lane B / Track G: Algebraic Adiabatic Curves — Perturbation Series (COMPLETE)

Rayleigh-Schrödinger perturbation series for Level 3 angular eigenvalues μ(R), exploiting linear matrix pencil H(R) = H₀ + R·V_C.

- **a₁ validated exactly:** -5.590189, matches Paper 13 Table II to ~10⁻¹⁵. Z-scaling formula confirmed for Z=1,2,3,5.
- **Perturbation series to order 20:** coefficients decay rapidly (|a_k| < 10⁻⁴ for k≥7)
- **Raw series convergence radius: R ≈ 2 bohr.** Padé [7/7] extends to R ≈ 3 bohr (< 1% error). All Padé orders fail beyond R ≈ 5 bohr.
- **μ(R) proven transcendental:** transitions from O(R) to O(R²) asymptotic behavior. No finite rational function spans both regimes. Confirms Paper 13 Sec XII.B.
- **P-matrix series** converges even more slowly (useful only R < 1 bohr)
- **Practical implication:** Point-by-point R-grid diagonalization cannot be eliminated globally. Series provides exact analytical derivatives at R=0 (spline seeding), validated algebraic coefficients, and definitive convergence boundary.
- 8 new tests, 83 total, all passing

**Files changed:**
- `geovac/algebraic_angular.py` — `perturbation_series_mu()`, `evaluate_perturbation_series()`, `pade_approximant()`, `evaluate_pade()`, `perturbation_series_P_matrix()`
- `tests/test_algebraic_angular.py` — 8 new tests

## [2.0.8] - 2026-03-28

### Three-lane concurrent sprint (first 3-PM parallel session)

#### Lane 1 / Track C: Spectral Radial Solver for Level 2 (COMPLETE)

Spectral Laguerre basis replaces 5000-point finite-difference radial ξ-solver in prolate spheroidal coordinates.

- **Method:** φₙ(ξ) = exp(-α(ξ-1)) · Lₙ(2α(ξ-1)) with α = √c² adapted to wavefunction decay. Gauss-Laguerre quadrature for matrix elements. Generalized eigenvalue problem Hc = λSc inside existing Brent c² root-finding loop.
- **250× dimension reduction:** 5000 FD grid points → 20 spectral basis functions
- **270× speedup:** 4.3s → 0.016s per energy evaluation at R=2.0 bohr
- **5000× accuracy improvement:** 1.01% → 0.0002% error at R=2.0 (FD error was masquerading as physics)
- **Spectral convergence:** n_basis=5 already at machine-level accuracy
- **H₂⁺ PES:** R_eq = 2.005 bohr (0.38% error vs exact 1.997), E_min matches exact to 0.01%
- FD solver preserved as `radial_method='fd'` (default unchanged); spectral via `radial_method='spectral'`
- 17 new tests (test_spectral_radial.py), 11/11 existing FD tests unchanged

**Files changed:**
- `geovac/prolate_spheroidal_lattice.py` — added `_radial_top_eigenvalue_spectral()`, `radial_method`/`n_basis` params
- `tests/test_spectral_radial.py` — new (17 tests)

#### Lane 2: Paper 14 Qubit Encoding — Expanded Benchmarks

Strengthened Paper 14 with genuine Gaussian baselines, fault-tolerant resource narrative, and composability discussion.

- **Genuine Gaussian baselines:** He cc-pVDZ (Q=10, 156 terms, λ=42.95) and cc-pVTZ (Q=28, 21,607 terms, λ=530.47) added to Table 1 and 1-norm table. Equal-qubit head-to-head: 1.3× fewer terms / 3.8× lower λ at Q=10; 8.1× fewer terms / 6.8× lower λ at Q=28.
- **Composed 1-norm table:** Trotter step counts for LiH (r≈4,500 at Q=84) and BeH₂ (r≈16,500 at Q=140) at ε=10⁻³
- **New Section V.D (Composability):** GeoVac structural sparsity composes additively with qubit tapering, Pauli grouping, and low-rank factorization (orthogonal axes)
- **New Section V.E (Fault-Tolerant Resources):** Explicit Trotter step estimates at ε=10⁻³ and 10⁻⁴; qubitization query complexity argument (Q^1.69 λ-scaling → ~3× query reduction at Q=110, ~10× at Q=1000)
- Commutator Trotter bounds (v2.0.3) verified present and correct (Q^1.47, 7× at Q=60)
- Synthetic Gaussian entries marked with † to distinguish from validated computed integrals
- 25/25 test_composed_qubit.py, 5/5 test_paper14_revision.py pass

**Files changed:**
- `papers/core/paper_14_qubit_encoding.tex` — expanded Sections III.G, IV.C, new V.D and V.E
- `tests/test_paper14_revision.py` — path fix (was pointing to nonexistent _revised.tex)

#### Lane 3: Level 3 Coupled-Channel Convergence (COMPLETE — structural ceiling identified)

Extended He coupled-channel convergence study to l_max=5. Definitive characterization of the adiabatic approximation floor.

- **Convergence through l_max=5 (q_mode='exact', n_channels=3):**
  | l_max | Error % | Wall Time |
  |------:|--------:|----------:|
  | 0 | 1.105 | 0.6s |
  | 1 | 0.311 | 0.8s |
  | 2 | 0.239 | 1.0s |
  | 3 | 0.220 | 1.4s |
  | 4 | 0.214 | 1.6s |
  | 5 | 0.211 | 2.0s |
- **Convergence rate:** Algebraic (~l_max⁻²), NOT exponential. Structural floor at ~0.19–0.20%.
- **Sub-0.1% NOT reachable** with coupled-channel P+Q machinery. Bottleneck is adiabatic approximation, not any numerical parameter.
- **3 channels sufficient:** n_channels=5 gives only 0.0006% improvement over n_channels=3
- **n_basis sensitivity (l_max=3):** n_basis=10: 0.283%, n_basis=30: 0.195%. Converges to ~0.19%.
- **Radial grid converged:** N_R_radial=3000 vs 5000 gives identical results.
- 2 new tests (17 total in test_algebraic_coupled_channel.py)

**Files changed:**
- `tests/test_algebraic_coupled_channel.py` — 2 new tests (tests 16-17)
- `debug/lmax_convergence_exact_v2.py` — convergence study script
- `debug/nbasis_sensitivity.py` — basis sensitivity script
- `debug/profile_coupled_channel.py` — profiling script
- `debug/data/lmax_convergence_exact_v2.json` — convergence data
- `debug/data/nbasis_sensitivity.json` — sensitivity data
- `debug/track_logs/lane3_coupled_channel_convergence.md` — track log

## [2.0.7] - 2026-03-28

### Track C: Level 2 Algebraic Audit
- Angular solver confirmed 100% algebraic (Legendre spectral basis)
- Radial ξ-solver identified as algebraicization target (5000-point FD → ~15 spectral functions)
- Literature review: Mitnik et al. (2021) and Kereselidze et al. (2016) provide proven spectral Sturmian bases
- Audit report: debug/level2_algebraic_audit.md

### Track D: Exact Q-Matrix
- Algebraic dP/dR from Hellmann-Feynman (verified vs FD to <1e-9)
- q_mode='exact' reduces coupled-channel error by 14-19% at l_max≥1
- Best result: 0.22% at l_max=3 (down from 0.27% with closure)
- l_max=0 overcorrection (1.1%) identified as structural (diagonal DBOC dominance)
- 5 new tests (15 total in test_algebraic_coupled_channel.py)

### First multi-agent PM session
- Two parallel tracks dispatched and completed via sub-agent architecture
- 54/54 tests pass, zero regressions

## [2.0.6] - 2026-03-28

### Track B: Algebraic Angular Solver and Coupled-Channel Integration

Algebraic spectral basis replaces finite-difference angular solver at Level 3. Coupled-channel integration resolves l_max convergence problem.

#### Algebraic angular solver
- Gegenbauer spectral basis: matrix dimension 200-800 (FD) → 10-30 (algebraic)
- Exact Hellmann-Feynman P-matrix from R-independent dH/dR
- He l_max=0: 0.16% single-channel (algebraic), 0.095% P-only coupled
- All angular matrix elements algebraic (Casimir, Gaunt, partial harmonic sums)

#### Sturmian variant (research investigation)
- Nuclear monopole >100% of V_coupling (validates Sturmian premise)
- 1.6-3.6x faster convergence at small n_basis for l_max=0
- Counterproductive at l_max≥1 (converges faster to wrong adiabatic limit)
- Preserved as research artifact, not recommended for production

#### Coupled-channel integration
- Algebraic coupling matrices feed into coupled-channel radial solver
- l_max convergence FIXED: single-channel diverges (0.16%→0.65%), coupled converges (0.37%→0.27%)
- P-only and P+Q bracket exact energy from below and above
- Q closure overcorrects at l_max=0 (1.1%) — improvement pathway identified (exact dP/dR)
- DBOC analysis: +0.035 Ha correction, 97% cancelled by off-diagonal coupling

#### Files added
- `geovac/algebraic_angular.py` — Algebraic angular solver (spectral basis)
- `geovac/algebraic_angular_sturmian.py` — Sturmian variant (research artifact)
- `geovac/algebraic_coupled_channel.py` — Coupled-channel integration
- `tests/test_algebraic_angular.py` — 28 tests
- `tests/test_algebraic_angular_sturmian.py` — 12 tests
- `tests/test_algebraic_coupled_channel.py` — 10 tests

#### Paper updates
- Paper 13: new section on algebraic angular solver and coupled-channel convergence

### Track A: l_max Divergence Root Cause Isolation

Definitive diagnostic identifying the source of the l_max divergence in composed LiH geometry.

#### Key findings
- HeH⁺ bare Level 4 drifts at +0.262 bohr/l_max with no Z_eff, no PK, no composition — divergence is intrinsic to the adiabatic multichannel expansion with asymmetric nuclear charges
- Mechanism: asymmetric nuclear coupling populates odd-l channels forbidden in homonuclear systems; these channels preferentially stabilize large-R configurations
- PK-induced symmetry breaking: composed LiH (Z_A=Z_B=1, symmetric) drifts at +0.303 because PK on atom A breaks A↔B symmetry, making the system effectively heteronuclear
- Channel weight evidence: H₂ has 94% in (0,0); HeH⁺ has only 68%, with 20% in odd-l channels

#### Negative results (new dead ends)
- Algebraic PK projector: rank-1 |core⟩⟨core| too weak; drift 5.6× worse than Gaussian PK (+1.697 vs +0.303 bohr/l_max)
- Enhanced Z_eff (PK-free): shape mismatch (1/r² vs 1/r) requires extreme A_rep; odd-l channel population 30.7% — worse than bare HeH⁺ (20%)
- Single-channel DBOC: 97% cancellation with off-diagonal P-coupling; overcorrects by 23×

#### Path forward identified
- The l_max divergence is an adiabatic approximation artifact
- The variational 2D solver (Paper 15) avoids this approximation and should eliminate the divergence
- Priority: reduce Level 4 iteration times to enable rapid testing of 2D solver in composition pipeline

#### Diagnostic scripts
- `debug/lmax_divergence_fast.py`, `debug/enhanced_zeff_quick.py`, `debug/coupled_channel_test.py`

## [2.0.5] - 2026-03-27

### Asymmetric Bond Diagnostic Arc (Tasks 1–11)

**Problem investigated:** l_max divergence for asymmetric bonds in composed geometry framework. R_eq drifts linearly at +0.23 bohr per l_max unit with l-dependent PK (no saturation).

#### Key Findings
- Angular wavefunction spreading is 100% nuclear (zero e-e contribution)
- Eigenvalue convergence is exponential at fixed geometry (β = 0.42/l_max)
- Eigenchannel rotation and prolate spheroidal basis provide no compression
- The Level 4 solver itself converges; divergence is a composed-geometry (PK) property

#### New Feature: R-Dependent PK Scaling
- `pk_mode='r_dependent'` option in `ComposedDiatomicSolver`
- Formula: w_PK(R) = δ_{l,0} × min(cap, R/R_ref)
- Reduces LiH R_eq error from 13.5% to 2.0% at l_max=4
- Advantage erodes at l_max ≥ 5; asymptotic drift rate unchanged
- R_ref = R_eq(expt) — ab initio derivation remains open

#### Also Tested (Negative Results Documented)
- Self-consistent PK: 50% drift reduction, intermediate performance
- Projected PK (one-shot): PES collapses, no equilibrium
- Channel-blind PK: diverges fastest (baseline)

#### Paper Updates
- Paper 17: expanded l_max divergence section, added R-dependent PK equation, updated conclusion
- Paper 15: added channel weight distribution data, clarified composed-geometry boundary

#### Files Added
- 13 diagnostic scripts in `debug/`
- 11 data files in `debug/data/`
- 36 plots in `debug/plots/`
- `diagnostic_arc_report.md` (comprehensive technical report)

## [2.0.4] - 2026-03-26

### H₂O Composed Solver

First triatomic PES from the composed natural geometry framework. Five-block decomposition: O 1s² core (Level 3, Z=8) + 2 O–H bond pairs (Level 4 with charge-center origin, Z_eff=6, Z_B=1) + 2 lone pairs (Level 3, Z_eff=6). Zero free parameters.

#### Architecture
- New module `composed_water.py`: `ComposedWaterSolver` with full pipeline (core → lone pair → PES → Morse fit)
- New module `lone_pair.py`: `solve_lone_pair()` — Level 3 hyperspherical solver for non-bonding electron pairs, `extract_channel_data_level3()` for coupling interface, `compute_inter_fiber_overlap()` for cross-fiber overlap
- Inter-fiber coupling at arbitrary bond angles via P_{l₁}(cos θ) · P_{l₂}(cos θ) Wigner rotation phases — backward compatible (θ=π recovers linear molecule formula)
- 6 coupling pairs: bond-bond (1, at θ_HOH), bond-lone (4, at θ_bond_lone), lone-lone (1, at θ_lone_lone)
- `coupling_pairs` parameter: 'all' (default), 'bond_bond' (validated only), 'none' (uncoupled)
- Wall time ~29s/point uncoupled, ~34s/point coupled

#### PES Results (uncoupled, charge-center origin)
- R_eq = 1.34 bohr (26% error vs 1.809 expt)
- Bound PES with correct qualitative shape (repulsive wall, minimum, dissociation)
- Bond-bond coupling is small (~0.5 Ha) and barely affects R_eq — consistent with BeH₂

#### Tests
- 21 tests in test_composed_water.py (core, bond pair, lone pair, PES scan, R_eq, no-coupling, coupling, nuclear repulsion, constructor, charge-center)
- 29 tests in test_bond_angle.py (bond pair charge asymmetry diagnostic)
- 13 tests in test_lone_pair.py (lone pair solver validation)

### Bond Pair Charge Asymmetry Diagnostic

Systematic study of Level 4 bond pair accuracy vs charge ratio Z_A/Z_B.

| Z_A/Z_B | R_eq error (midpoint) | R_eq error (charge-center) |
|:--------:|:---------------------:|:--------------------------:|
| 1:1 (H₂) | 0.1% | — |
| 2:1 (BeH) | 18% | — |
| 6:1 (OH) | 56% | 17% |

The charge-center origin z0 = R(Z_A−Z_B)/(2(Z_A+Z_B)) dramatically improves R_eq for asymmetric bonds by centering coordinates near electron density rather than at the geometric midpoint.

### Inter-Fiber Coupling at Arbitrary Bond Angles

Generalization of the P_{l}(-1) = (-1)^l linear-molecule rotation to P_{l₁}(cos θ) · P_{l₂}(cos θ) for general bond angle θ. Applied in `_channel_rotation_phases()` and `full_exchange_inter_fiber_energy()`. For θ = π, recovers existing formula exactly.

### Lone Pair Coupling — Negative Result

Coupling decomposition at R = 1.81 bohr (experimental R_eq):

| Pair type | Energy (Ha) | % of total |
|:----------|:-----------:|:----------:|
| Bond-Bond (1 pair) | −0.59 | 1.4% |
| Bond-Lone (4 pairs) | −27.59 | 63.5% |
| Lone-Lone (1 pair) | −15.27 | 35.1% |
| **Total** | **−43.46** | 100% |

The S·F⁰ approximation breaks down at Z_eff=6: Slater integrals scale as ~Z, producing coupling energies that exceed the total electronic energy. Bond-bond coupling (~0.5 Ha) remains physically reasonable and consistent with BeH₂ validation. Lone pair coupling disabled for production use.

### Origin Audit

- Charge-center origin applied to H₂O bond pairs (6:1 asymmetry): R_eq error 56% → 17%
- BeH₂ (2:1) and LiH (1:1 effective) left at midpoint origin where the effect is small
- `origin='charge_center'` parameter in `solve_level4_h2_multichannel()` — available since v2.0.3

#### Paper Updates
- **Paper 17:** Section VI.A gains H₂O outlook paragraphs (architecture, results, coupling decomposition, bottleneck identification)

## [2.0.3] - 2026-03-25

### Commutator-Based Trotter Bounds

First-order commutator Trotter error bound via symplectic Pauli algebra. Exploits the fact that selection-rule-sparse GeoVac Hamiltonians have most Pauli term pairs commuting, yielding a dramatically tighter simulation cost bound.

#### Implementation
- `pauli_commutator_bound()` in trotter_bounds.py: first-order commutator bound via symplectic inner product on Pauli binary representations
- `_pauli_string_to_binary()`, `_symplectic_inner_product_matrix()`: algebraic Pauli commutativity check (O(N²Q), no Hilbert-space matrices)
- `analyze_commutator_cost()`: multi-epsilon comparison of commutator vs 1-norm bounds
- Chunked computation for large N (memory-safe up to N ~ 30,000)

#### Results

| System | Q | N_terms | Anticomm % | r_comm | r_1norm | Tightening |
|:-------|:--|:--------|:-----------|:-------|:--------|:-----------|
| He nmax=2 | 10 | 120 | 28.0% | 60 | 253 | 0.056 |
| He nmax=3 | 28 | 2,659 | 21.1% | 303 | 1,753 | 0.030 |
| He nmax=4 | 60 | 31,039 | 12.2% | 821 | 5,849 | 0.020 |

- Commutator Trotter scaling: r ~ Q^{1.47} (vs Q^{1.69} 1-norm bound)
- Anticommuting fraction: 28% at Q=10, 12% at Q=60, scales as Q^{-0.45}
- 7× fewer Trotter steps at Q=60 (821 vs 5,849 at ε=10⁻³)
- nmax=5 skipped (N×N matrix ~ 52 GB); 3-point fit sufficient for scaling

#### Tests
- 8 new tests in test_trotter_bounds.py (24 total): binary conversion, symplectic product, commutator bound correctness, tightening vs 1-norm, anticommuting fraction decrease with system size

#### Paper Updates
- **Paper 14:** Sec III.F updated with commutator bound derivation, results table, and Q^{1.47} scaling fit; Sec III.H notes structural origin of advantage; Limitations section updated (commutator bounds computed, no longer pending); scaling summary table gains commutator row

## [2.0.2] - 2026-03-25

### Full 1-RDM Inter-Fiber Exchange

Off-diagonal one-particle density matrix exchange closes the gap between diagonal S·F⁰ and the fitted exchange model, with zero free parameters. Kinetic orthogonalization tested and found negligible.

#### Full Exchange Implementation
- `full_exchange_inter_fiber_energy()` in inter_fiber_coupling.py: contracts off-diagonal channel 1-RDM with parity-transformed Slater F⁰ matrix via Eq. K_AB = -Σ [γ_{l₁,l₁'}]² · (-1)^{l₁+l₁'} · F⁰_{l₁,l₁'}
- `extract_channel_1rdm()`: traces electron 2 from Level 4 eigenvector to obtain γ_{l₁,l₁'}, |F(R_e)|²-weighted over hyperradius
- `_compute_l1_indexed_f0_matrix()`: aggregates per-channel Be-centered densities into l₁-marginal F⁰
- `kinetic_orthogonalization_energy()`: Löwdin kinetic correction from channel overlap × centrifugal T
- `interbond_mode='full_exchange'` and `'full_exchange_kinetic'` in composed_triatomic.py

#### BeH₂ Results

| Method | R_eq (bohr) | Error vs 2.507 |
|:-------|:-----------:|:--------------:|
| Block-diagonal (no coupling) | 3.28 | 31% |
| Monopole F⁰ (Phase 1) | 4.06 | 62% |
| Exchange S·F⁰ (Phase 3) | 3.01 | 20% |
| Full 1-RDM exchange | 2.80 | 11.7% |
| Fitted model −K·S (diagnostic) | 2.81 | 12% |

- Off-diagonal γ elements are 68–121% of diagonal (grows with R as channel mixing increases)
- Off-diagonal exchange terms are subtractive (~10% of diagonal at R_eq), providing correct R-dependence
- Kinetic orthogonalization: +0.03 Ha uniform, no R_eq effect (negative result)
- Residual 11.7% attributed to basis truncation (l_max=2) and adiabatic approximation

#### LiH Transferability (Not Applicable)
- Full 1-RDM exchange is specific to bond–bond (Level 4 + Level 4) topology
- LiH core–valence fibers have different channel structures (Level 3 vs Level 4); cannot share γ matrices
- PK pseudopotential already encodes core–valence orthogonality
- Cross-level overlap formalism deferred to future work

#### Tests
- 13 tests in test_1rdm_exchange.py: 1-RDM properties (diagonal consistency, Hermiticity, trace, off-diagonal), full exchange (sign, magnitude, R-dependence, R_eq improvement, parity symmetry), kinetic orthogonalization (positive, R-dependent, reasonable magnitude), combined mode

### Algebraic Methods & l-Dependent Pseudopotential

Algebraic reformulation of core screening, inter-fiber coupling, and nuclear coupling audit. l-dependent PK pseudopotential improves LiH R_eq from 6.9% to 5.3%. Higher Slater multipoles confirmed as dead end.

#### Algebraic Core Density (Phase 1)
- **CoreScreening** gains `method='algebraic'` (now default): replaces histogram binning with coordinate-transform identity using Level 3 channel coefficients
- **Z_eff closed-form discovery:** Z_eff(r) ≈ 1 + 2(1 + ar + br²)exp(-cr) fits Li⁺ at 0.7% accuracy; deviations from hydrogenic ratios quantify electron correlation
- New function: `compute_core_density_algebraic()` — accepts pre-solved Level 3 result
- All 21 core_screening tests pass, all 33 composed_diatomic tests pass

#### l-Dependent Phillips-Kleinman Pseudopotential (Phase 2)
- **ComposedDiatomicSolver** gains `pk_channel_mode='l_dependent'`: scales PK barrier per-electron by δ_{l_i, 0}
- LiH R_eq improves from 6.9% to 5.3% at l_max=2 (channel-blind PK is a bug, not a feature — 1s² core has zero overlap with l>0)
- l_max divergence reduced but not eliminated: 14.7% at l3 (vs 16.9% channel-blind), 26.2% at l4 (vs 32.8%)
- 8 new unit tests, 6 new integration tests in test_l_dependent_pk.py

#### Channel-Decomposed Inter-Fiber Coupling (Phase 3)
- **inter_fiber_coupling.py**: new algebraic pathway (`method='algebraic'`) eliminates density-extraction + shell-theorem + cumulative-charge pipeline
- `extract_channel_data()` — single angular solver pass shared by density, overlap, and F⁰ computations
- `compute_channel_fk_matrix(k)` — channel-channel Slater Fᵏ decomposition at arbitrary multipole order
- Direct exchange formula: E_direct = −Σ_ch (−1)^{l₁+l₂} F⁰_ch (avoids S·F⁰ factorization)
- `multipole_exchange_inter_fiber_energy(k_max)` — generalized multipole exchange
- 8 new algebraic tests, all 32 backward-compatibility tests pass

#### Higher Slater Multipoles — Negative Result
- k=1 (dipole) is repulsive (+0.009 Ha), k=2 (quadrupole) is negligible
- Total k=0+1+2 gives 20.5% R_eq error (worse than monopole alone at 20.0%)
- Added to CLAUDE.md Section 3 (Approaches That Failed)

#### Nuclear Coupling Audit
- Confirmed `compute_nuclear_coupling` in level4_multichannel.py already uses algebraic split-region Legendre expansion terminated by Wigner 3j triangle inequality — no quadrature
- Gauss-Legendre quadrature scoped to Z_eff screening correction only

#### Paper Updates
- **Paper 15:** Nuclear coupling section corrected (was "Gauss-Legendre quadrature", now describes actual algebraic method with split-region Legendre expansion and Gaunt integrals)
- **Paper 17:** Algebraic core density, l-dependent PK, channel-decomposed exchange, multipole negative result, path forward narrowed to 1-RDM exchange + kinetic orthogonalization
- **Paper 13:** Cross-reference to nuclear coupling algebraization added after Gaunt integral definition
- **CLAUDE.md:** Updated best results (LiH 5.3%), new failed approach (higher Slater multipoles), algebraic structure note in Section 5

## [2.0.1] - 2026-03-24

### BeH₂ Polyatomic Diagnostic & Failed Approach Documentation

Systematic validation of the composed natural geometry framework on BeH₂ (simplest polyatomic), identifying the block-diagonal approximation boundary.

#### New Module
- **composed_triatomic.py**: BeH₂ composed PES solver with 5-block decomposition (core + 2 bond pairs × [Be-side + H-side]), symmetric stretch scan, inter-bond repulsion model

#### BeH₂ PES Results
- Qualitatively correct PES: bound state, repulsive wall, correct dissociation limit — a positive finding for a zero-parameter method
- R_eq = 3.31 bohr (32% error vs 2.507 expt), D_e overestimated by 19×
- Root cause identified: block-diagonal approximation is insufficient for valence-valence coupling between bond pairs sharing the same nucleus

#### Failed Approaches Tested and Documented
- **Z_eff partitioning** (4 schemes: equal, valence_scaled, sqrt, per-bond): all worsen R_eq, best result 103% error. Charge redistribution cannot substitute for missing inter-bond interaction
- **Classical inter-bond repulsion**: point-charge V_inter(R) between bond pairs has wrong R-dependence, pushes R_eq outward (wrong direction) at all scaling factors (0.25–1.5)

#### Key Finding
- Block-diagonal approximation validated for core-valence separation (LiH 6.4%) where energy scales differ, but breaks down for valence-valence separation where bond pairs compete for nuclear attraction
- Orbital-level inter-bond coupling (exchange, orthogonalization) is the active research problem
- Pauli scaling (Q^2.5) confirmed unaffected — scalar corrections don't change qubit Hamiltonian structure

#### Phase 1 Monopole Inter-Fiber Coupling (Negative Result)

Implemented monopole (k=0) Slater F⁰ inter-fiber Coulomb coupling between BeH₂'s two bond-pair fibers.

- **New module: inter_fiber_coupling.py** — standalone pipeline: extract one-electron density from Level 4 wavefunction (re-runs angular solver at sample R_e points), transform to Be-centered coordinates via shell theorem convolution, compute Slater F⁰ integral using O(N) cumulative charge method
- **Functions:** `extract_origin_density()`, `transform_to_center_density()`, `slater_f0_integral()`, `monopole_inter_fiber_energy()`
- **Integration:** `composed_triatomic.py` gains `interbond_mode='monopole'` option for PES scan
- **Result:** Monopole energy is ~1.6 Ha and nearly R-independent near equilibrium. R_eq worsens from 32% to 62% error — same failure mode as classical inter-bond repulsion (CLAUDE.md Section 3). Scalar Coulomb corrections cannot fix R_eq because they lack differential R-dependence
- **6 new tests** in test_composed_beh2.py: monopole positive, decreasing with R, flat at short R, electron count ≈ 2, F⁰ hydrogenic validation (5Z/8 = 0.625), norm-preserving coordinate transform

#### Inter-Fiber Overlap Diagnostic (Positive Result)

Computed the inter-fiber channel overlap S_avg(R) = Σ (−1)^{l₁+l₂} |c⁰_{l₁l₂}|² as a diagnostic for whether exchange coupling can provide the missing R-dependence.

- **New function: `compute_overlap_diagnostic()`** in inter_fiber_coupling.py — extracts channel expansion coefficients from Level 4 angular eigenvectors, computes parity-weighted overlap, returns |F|²-weighted average over hyperradius
- **Diagnostic script:** debug/overlap_diagnostic.py — scans 14 R-points (2.0–6.0 bohr), produces table + 4-panel plot
- **Key result:** S_avg has strong R-dependence (0.47 at R=2.0, 0.24 at R=6.0, range 0.225). Slope dS/dR = −0.16 at R_eq. The (0,0) channel weight drops from 61% to 20% as odd-parity channels grow with R
- **Exchange feasibility fit:** E_exch = −K·S(R) with K=4.08 Ha reduces R_eq error from 31% to 12% (spline-interpolated). Both S^1 and S^2 models pass the < 15% target. Phase 3 exchange implementation is warranted
- **Data:** debug/data/overlap_diagnostic.json, debug/plots/overlap_vs_R.png

#### Phase 3 Exchange Inter-Fiber Coupling (Partial Success)

Implemented exchange inter-fiber coupling using the channel overlap diagnostic to modulate the monopole Coulomb integral.

- **New function: `exchange_inter_fiber_energy()`** in inter_fiber_coupling.py — computes E_exch(R) = −S_avg(R) · F⁰(R) from pre-computed overlap and monopole data
- **Integration:** `composed_triatomic.py` gains `interbond_mode='exchange'` option for PES scan
- **Result:** R_eq improves from 31% (block-diagonal) to 20% (exchange). The ab initio coupling captures ~40% of the fitted correction (K_eff ≈ 1.6 vs K_fit = 4.08 Ha)
- **Remaining gap (~60%):** off-diagonal 1-RDM exchange integrals (not the direct F⁰), higher multipoles k≥1, kinetic energy corrections from Löwdin-type orthogonalization
- **4 new tests** in test_composed_beh2.py: exchange energy negative, decreasing with R, R_eq improvement over block-diagonal, K_eff magnitude

**Progression summary:**

| Method | R_eq (bohr) | Error vs 2.507 |
|:-------|:-----------:|:--------------:|
| Block-diagonal (no coupling) | 3.28 | 31% |
| Monopole F⁰ (Phase 1) | 4.06 | 62% |
| Exchange S·F⁰ (Phase 3) | 3.01 | 20% |
| Fitted model −K·S (diagnostic) | 2.81 | 12% |

**Key finding:** The inter-fiber channel overlap S(R) — a purely algebraic quantity from angular eigenvector coefficients and parity phases — captures the R-dependent coupling between bond-pair fibers. This validates the fiber bundle framework for polyatomics.

#### LiH l_max Sigma+Pi Negative Result

Attempted extending LiH composed solver to l_max=3,4 with sigma+pi (m=1) channels, following Paper 17's prediction that this would reduce R_eq error from 6.4% to ~2–3%. The result is the opposite — R_eq diverges monotonically from experiment:

| Configuration | R_eq (bohr) | Error vs 3.015 |
|---------------|:-----------:|:--------------:|
| l_max=2 sigma (manual PK) | 2.97 | 1.5% |
| l_max=2 sigma (ab initio PK) | 3.21 | 6.4% |
| l_max=3 sigma+pi (ab initio PK) | 3.52 | 16.9% |
| l_max=4 sigma+pi (ab initio PK) | 4.00 | 32.8% |

- **Root cause:** PK pseudopotential is channel-blind (derived from 1s² core). Higher-l channels add correlation preferentially at large R (diffuse valence cloud), shifting PES minimum outward
- **Code changes:** `m_max` and `l_max_per_m` parameters added to `ComposedDiatomicSolver` (backward-compatible)
- **Next direction:** l-dependent or m-dependent pseudopotential, or multi-channel adiabatic corrections (DBOC)

#### Documentation Updates
- CLAUDE.md: three new entries in Section 3 (Approaches That Failed), updated active frontier
- README.md: polyatomic accuracy limitation, BeH₂ benchmark row, corrected l_max prediction
- Paper 14: new limitation item in Section IV.E (valence-valence coupling)
- Paper 17: l_max convergence table and root cause analysis in Outlook section, BeH₂ diagnostic note

#### Tests
- 36 tests passing (20 unit + 6 monopole/overlap + 10 integration) in test_composed_beh2.py and test_composed_triatomic.py
- 11 new tests in test_composed_diatomic.py for l_max=3,4 sigma+pi (TestLiHL3SigmaPi, TestLiHL4SigmaPi, TestLmaxConvergenceReport)

## [2.0.0] - 2026-03-23

### Composed Multi-Center Qubit Hamiltonians

Block-diagonal qubit Hamiltonians for composed multi-center molecules (LiH, BeH₂, H₂O) with universal Q^2.5 Pauli scaling, representing 51×–1,712× fewer Pauli terms than published Gaussian baselines.

#### New Module
- **composed_qubit.py**: Composed multi-center qubit Hamiltonians via pseudopotential block decomposition
  - `build_composed_lih()` — 3-block LiH (core + Li-valence + H-valence)
  - `build_composed_beh2()` — 5-block BeH₂ (core + 2 bond blocks × [Be-side + H-side])
  - `build_composed_h2o()` — 7-block H₂O (core + 2 bond blocks + 2 lone pair blocks)
  - Block-diagonal ERIs from composed natural geometry pseudopotential (cross-block ERIs exactly zero)
  - Hydrogenic orbitals with per-block Z_eff screening (Z_eff = Z - n_core)

#### Key Results: Universal Composed Scaling
| Molecule | Blocks | Q (nmax=2) | Pauli Terms | Scaling α | vs Gaussian |
|----------|:------:|:----------:|:-----------:|:---------:|:-----------:|
| LiH | 3 | 30 | 334 | 2.50 | 51× fewer |
| BeH₂ | 5 | 50 | 556 | 2.51 | 460× fewer |
| H₂O | 7 | 70 | 778 | 2.52 | 1,712× fewer |

Exponent spread of 0.02 across three molecules demonstrates topology-independent universality.

#### Published Gaussian Baselines
- Gaussian Pauli term counts from Trenev et al. (Quantum, 2025): LiH Q^4.25, H₂O Q^3.92
- At equal qubit counts, composed GeoVac advantage grows with system size due to ~1.5 exponent gap

#### Paper 14 Revised
- New Section IV: Composed Systems (block decomposition, scaling analysis, Gaussian comparison)
- Added Trenev et al. (2025) reference for published molecular Gaussian Pauli counts
- Fixed H₂O Z_eff typo: 2 → 6 (= Z_O - n_core = 8 - 2)
- Old version archived to `papers/archive/paper_14_qubit_encoding_v1.9.0.tex`

#### Tests
- **test_composed_qubit.py**: 25 tests (LiH block structure, Pauli terms, cross-block ERIs, scaling)
- **test_composed_beh2.py**: 14 tests (BeH₂ 5-block structure, Pauli terms, bond symmetry)
- **test_composed_h2o.py**: 17 tests (H₂O 7-block structure, Pauli terms, lone pair blocks)
- **test_paper14_revision.py**: 5 tests (Paper 14 revision validation)
- Total: 61 new tests

---

## [1.9.0] - 2026-03-23

### VQE Benchmark Infrastructure & Validated Gaussian Baselines

End-to-end VQE benchmark pipeline comparing GeoVac lattice encodings against genuine Gaussian-basis Hamiltonians at equal qubit counts. Replaced estimated Pauli term counts with actual computed values from a custom single-center Gaussian integral engine.

#### New Modules
- **vqe_benchmark.py**: VQE benchmark pipeline with OpenFermion→Qiskit conversion
  - `build_geovac_he(max_n)`, `build_gaussian_he()`, `build_gaussian_h2()` — system builders
  - `run_vqe()` — COBYLA-based VQE with EfficientSU2 ansatz
  - `collect_static_metrics()` — CX count, circuit depth, QWC groups without VQE execution
  - `format_comparison_table()`, `save_results()` — formatted output and JSON persistence
- **debug/compute_he_gaussian_integrals.py**: Single-center Gaussian integral engine for He
  - No PySCF dependency — uses BSE basis set parameters + scipy/numpy
  - Slater R^k radial integrals via incomplete gamma function (analytical inner integral)
  - Real spherical harmonic Gaunt coefficients via numerical quadrature
  - Validated against published FCI energies (Woon & Dunning, JCP 100, 2975, 1994)

#### New Benchmark
- **benchmarks/vqe_head_to_head.py**: Head-to-head VQE comparison
  - He GeoVac nmax=2,3 vs Gaussian STO-3G (original VQE comparison)
  - H2 Gaussian STO-3G VQE benchmark
  - Equal-qubit accuracy comparison at Q=2, 4, 10, 28
  - Pauli term scaling crossover analysis with plot (`debug/plots/pauli_crossover.png`)
  - Headline accuracy-per-qubit summary

#### Updated gaussian_reference.py
- **`he_cc_pvdz()`**: Hardcoded computed MO integrals for He cc-pVDZ [2s1p]
  - 5 spatial orbitals, 10 qubits, 156 Pauli terms
  - FCI energy: -2.8876 Ha (published: -2.8877 Ha, within 0.0001 Ha)
- **`he_cc_pvtz()`**: Cached MO integrals for He cc-pVTZ [3s2p1d]
  - 14 spatial orbitals, 28 qubits, 21,607 Pauli terms
  - FCI energy: -2.9002 Ha (published: -2.9003 Ha, within 0.0001 Ha)
  - Cache file: `geovac/cache/he_cc_pvtz_mo_integrals.npz`

#### Key Results: Equal-Qubit Comparison
| Q | GeoVac | Err% | Pauli | Gaussian | Err% | Pauli |
|:-:|--------|:----:|:-----:|----------|:----:|:-----:|
| 10 | nmax=2 | 0.56% | 120 | cc-pVDZ | 0.56% | 156 |
| 28 | nmax=3 | 0.45% | 2,659 | cc-pVTZ | 0.12% | 21,607 |

- At 10 qubits: comparable accuracy, GeoVac uses 1.3× fewer Pauli terms
- At 28 qubits: GeoVac uses 8.1× fewer Pauli terms; Gaussian is 3.8× more accurate
- GeoVac advantage grows with system size due to structural sparsity

#### Correction: Gaussian Atomic Scaling
- Old estimate (H₂ molecules only): Gaussian Pauli terms scale as Q^4.60
- Actual (atoms + molecules combined): Q^3.56 for atoms (spherical symmetry zeros many ERIs)
- Estimated counts were ~5× too high: 814 → 156 (cc-pVDZ), 91,923 → 21,607 (cc-pVTZ)
- All comparisons now use actual computed Pauli counts, not scaling extrapolations

#### Tests
- **test_gaussian_reference.py**: 23 tests (was 12)
  - 6 new cc-pVDZ tests: metadata, h1/ERI symmetry, FCI energy, Pauli count (156), p-orbital degeneracy
  - 5 new cc-pVTZ tests: metadata, h1/ERI symmetry, FCI energy, Pauli count (21,607)
- **test_vqe_benchmark.py**: 12 tests
  - System builders, VQE execution, metric collection, OpenFermion→Qiskit conversion
- All 35 new tests passing

## [1.8.0] - 2026-03-22

### Quantum Simulation Cost Analysis (Paper 14 Expansion)

Comprehensive expansion of Paper 14 (qubit encoding) with three new code modules, benchmark infrastructure, and three new results sections quantifying measurement and simulation costs.

#### New Modules
- **gaussian_reference.py**: Genuine Gaussian baseline Hamiltonians from published integrals
  - `h2_sto3g(R)` — H₂ STO-3G from Szabo & Ostlund Table 3.15 (Hehre/Stewart/Pople JCP 1969)
  - `he_sto3g()` — He STO-3G from analytical single-center integrals
  - `build_qubit_hamiltonian()` — JW transform for Gaussian systems
- **measurement_grouping.py**: Qubit-wise commuting (QWC) measurement group analysis
  - `qwc_groups(qubit_op)` — Greedy sorted-insertion QWC partitioning
  - `analyze_measurement_cost(qubit_op)` — Full measurement cost metrics
  - `compare_measurement_cost()` — Side-by-side GeoVac vs Gaussian comparison
- **trotter_bounds.py**: First-order Trotter error bounds and simulation cost
  - `pauli_1norm(qubit_op)` — Pauli coefficient 1-norm (lambda)
  - `trotter_steps_required()` — Steps for target accuracy from Campbell PRL 2019
  - `analyze_trotter_cost(qubit_op)` — Full Trotter cost analysis

#### New Benchmark
- **benchmarks/qubit_scaling_sweep.py**: Systematic sweep across nmax=2–5
  - Collects Pauli terms, QWC groups, 1-norms, Trotter steps
  - Fits power-law scaling exponents with R² values
  - Outputs `qubit_scaling_data.json` and `qubit_scaling_summary.txt`

#### Scaling Exponents (GeoVac He, 2e)

| Metric | Exponent | R² | Gaussian baseline |
|--------|:--------:|:--:|:-----------------:|
| Pauli terms | Q^3.15 | 0.9995 | Q^4.60 (conventional) |
| QWC measurement groups | Q^3.36 | 1.0000 | — |
| Pauli 1-norm (λ) | Q^1.69 | 0.9972 | — |
| Trotter steps (ε=10⁻³) | Q^1.69 | 0.9972 | — |

#### Paper 14 Expanded (~590 → ~850 lines)
- **Sec III.E:** Measurement cost — QWC grouping analysis, O(Q^3.36) scaling
- **Sec III.F:** Simulation cost — Pauli 1-norm, first-order Trotter bounds, O(Q^1.69) scaling
- **Sec III.G:** Genuine Gaussian comparison — published integrals replace placeholder estimates
- **Updated:** Abstract, discussion, limitations, and conclusion with new metrics
- **New bibliography:** Campbell PRL 2019 (Trotter bounds), Hehre/Stewart/Pople JCP 1969 (STO-3G basis)

#### ERI Density Decay
- 10.4% → 3.9% → 2.1% → 1.3% over nmax=2–5
- Confirms ~1/M² selection rule from angular momentum conservation

#### Tests
- 48 new tests across 3 files:
  - `tests/test_gaussian_reference.py` — 12 tests (integral validation, JW consistency)
  - `tests/test_measurement_grouping.py` — 20 tests (QWC correctness, grouping, scaling)
  - `tests/test_trotter_bounds.py` — 16 tests (1-norm, Trotter steps, cost analysis)

## [1.7.7] - 2026-03-22

### Epistemic Calibration & Rhetorical Tightening

No code changes. Rhetoric calibrated across the core paper series to properly distinguish known results from new contributions, and to use consistent terminology for κ = −1/16.

#### Papers Revised
- **Paper 0:** Rewrote Section I introduction — replaced top-down axiomatic presentation with bottom-up discovery narrative (naive spatial discretization produces angular momentum quantum numbers unsolicited). Axioms now presented as distillation of the construction, not the starting point.
- **Paper 7:** Calibrated abstract to attribute Fock's 1935 S³ equivalence properly, foregrounding discrete graph convergence as this paper's new contribution. "Dimensionless Vacuum Hypothesis" → "Scale-Invariance Property." Added Fock attribution in Coulomb potential discussion (Sec VI.C). Concluding remarks now explicitly distinguish known (Fock's scale invariance) from new (graph convergence + organizing principle).
- **Paper 16:** Title changed from "Chemical Periodicity as S_N Representation Theory" to "The Periodic Table Recast as S_N Representation Theory" — signals reframing rather than derivation.
- **FCI-A (paper_fci_atoms):** Added benchmark scope paragraph acknowledging STO-3G is a minimal basis (1969), not a production standard; positions comparison honestly as minimal-basis test.
- **Paper 6:** Removed dead TODO comment block from abstract (content was already written).

#### Terminology Audit (All Papers)
- **"universal kinetic scale" / "universal topological constant" → "kinetic scale factor"** throughout all .tex files (Papers 0, 2, 5, 7, 8-9, 14, FCI-A)
- κ = −1/16 now consistently described as determined by Rydberg matching, not derived from topology alone

## [1.7.6] - 2026-03-22

### Paper 2: Universal Algebraic Identity & Hopf Fibration Generalization

No code changes. Paper 2 updated with Phase 5 (residual analysis) and Phase 6 (Hopf generalization) findings.

#### Paper 2 Updated
- **Proven:** Universal algebraic identity B_formal(m)/N(m) = d at m = 2 for ALL round spheres S^d, via quadratic m² + dm − 2(d+2) = 0
- **Proven:** B_formal = B_branch ONLY for d = 3 (unit embedding index of SO(3) ⊂ SO(4)); this is the true source of S³ specificity
- **Established (negative):** Quaternionic Hopf (S³→S⁷→S⁴) gives 1/α ≈ 970 — no match to weak coupling α_W ≈ 1/29.5
- **Established (negative):** Octonionic Hopf (S⁷→S¹⁵→S⁸) gives 1/α ≈ 7164 — no match to strong coupling α_s ≈ 1/8.5
- **Established (negative):** Formula residual (8.8×10⁻⁸) is intrinsic — matches no QED correction, spectral invariant, or mathematical constant
- **Revised:** Sec VIII.E "S³ specificity" with sharper B_formal vs B_branch distinction and three-point representation-theoretic obstruction
- **Added:** Quaternionic/octonionic Hopf negative results to Sec VIII.E and conclusion
- **Added:** Bibliography entry for Phase 6 investigation

#### Investigation Scripts (`debug/kappa_investigation/`)
- `phase5/phase5_residual.py` — Residual structure analysis (50-digit precision, QED corrections, CODATA vintages)
- `phase5/phase5_summary.md` — Phase 5 summary: residual is intrinsic, no known correction matches
- `phase6/phase6_quaternionic_hopf.py` — All four Hopf fibrations, B_formal vs B_branch, representation theory
- `phase6/phase6_summary.md` — Phase 6 summary: universal identity, quaternionic negative result, obstruction analysis

## [1.7.5] - 2026-03-22

### Paper 2: Kappa Investigation — Transition Operators, Second Selection, Circulant Hermiticity

No code changes. Paper 2 updated with findings from three-phase kappa investigation.

#### Paper 2 Updated
- **Added:** Sec VIII — "Transition Operators and the Hopf Geometry": L± = Hopf fiber motion, T± = radial S3 motion, d_max = 4 from 2 motion types × 2 directions
- **Added:** Sec VIII — "Second Selection Principle": bipartite spectral radius 2·d_max = 8 matches |λ₃| = n²−1 = 8, independently selecting n_max = 3
- **Proven:** Circulant coupling (b, c) are necessarily complex conjugates; M is Hermitian, not real symmetric (Sec VI D)
- **Proven:** Physical root at θ ≈ 3π/2 (nearly pure imaginary coupling); determinant condition fixes cos(3θ) = −1/(2(K/3)^{3/2})
- **Proven:** Self-referential property (b+c satisfies same cubic) is the k=0 DFT mode of any circulant — structural, not coincidental
- **Established (negative):** Best additive bridges to B = 42 are π/72 (0.002%) and 1/24 (0.003%), neither derived from first principles
- **Established (negative):** Circulant does not factor K into independent spectral invariants of individual bundle components
- **Updated:** Link 4 in derivation chain with footnote on forced Hermiticity result
- **Updated:** Link 2 footnote noting second selection principle as independent support
- **Added:** Two new open questions (common origin of selection principles; circulant phase ↔ Hopf connection)

#### Investigation Scripts (`debug/kappa_investigation/`)
- `kappa_investigation.py` — Phase 1: lambda_max analysis, d_max counting, kappa self-consistency
- `phase2/phase2_spectral_determinants.py` — Phase 2: spectral determinant bridge search, d_max correction
- `phase3/phase3_circulant.py` — Phase 3: circulant polar parameterization, Hermiticity proof, DFT analysis

## [1.7.4] - 2026-03-22

### Paper 2: Algebraic Selection Principle & Spectral Geometry Survey

No code changes. Paper 2 updated with proven results and systematic negative results.

#### Paper 2 Revised
- **Proven:** Selection principle B/N = 3 selects n_max = 3 uniquely (closed form B/N = 3(m+2)(m-1)/10, quadratic (m+4)(m-3) = 0)
- **Proven:** Per-shell Casimir identity ⟨C₃⟩_n = C₄(n) = (n²-1)/2 (SO(3) average equals SO(4) Casimir)
- **Proven:** S³ specificity — selection principle has no analogue on S^d for d ≠ 3
- **Established (negative):** Spectral determinants, analytic torsion, heat kernels, and APS invariants cannot derive K or K/π
- **Established (negative):** Other cubic roots match no known Standard Model parameters at 1% level
- **Established (negative):** QED running confirms formula gives α at Thomson limit; residual is intrinsic
- **Fixed:** det'(S²) formula (added missing 1/2 term in exponent)
- **Fixed:** B(4), B(5) values in selection principle table (132→162, 310→462)
- **Upgraded:** Link 2 from ◐ to ● (largely established) in derivation chain
- **Updated:** Conclusion reflects algebraic proof and spectral geometry negative results

#### New Debug Scripts (`debug/`)
- `alpha_cubic_roots.py` — Cubic root survey and QED running analysis
- `alpha_spectral_bridge.py` — Spectral determinant bridge search (near-miss investigation)
- `alpha_torsion_heat.py` — Analytic torsion, heat kernel, and index theory survey
- `alpha_rep_theory.py` — Representation-theoretic attack on B = 42 and selection principle

#### Principles Applied
- All mathematical content preserved — proven results added, negative results documented
- Paper remains in `papers/conjectures/` (combination rule still not derived)
- Epistemic standards: "we prove" used only for algebraically proven results

## [1.7.3] - 2026-03-22

### Paper 2 Rewrite: Statistical Validation & Circulant Structure

Complete rewrite of Paper 2 (fine structure constant from Hopf bundle) based on four-phase computational audit. No code changes.

#### Paper 2 Revised
- **Fixed:** Boundary term formula (Δ = 1/(|λ₃|×N(2)) = 1/40, where N(2) = 5 is cumulative state count, not g₂ = 4 as previously stated)
- **Fixed:** Half-trace notation removed; B = 42 is the full degeneracy-weighted Casimir trace
- **Added:** Section V — Statistical Validation (p-value = 5.2×10⁻⁹ from combinatorial search over 1.92×10⁹ formulas)
- **Added:** Section VI — Circulant Structure (cubic is characteristic polynomial of traceless Z₃-symmetric circulant matrix)
- **Added:** Section VII — What Is and Is Not Established (five-link derivation chain with explicit gap identification)
- **Added:** Spectral determinant near-miss as open question (det'₁·det'₃/π ≈ 41.957, 0.1% from B = 42)
- **Removed:** QED charge renormalization analogy (α² does not match any perturbative order, confirmed Phase 4)
- **Removed:** "Single chirality" subsection and factors of 1/2
- **Removed:** Overclaiming language ("greatest damn mystery has a geometric answer")
- **Tone:** "We observe" throughout; abstract states "empirical formula with structural support, not a first-principles derivation"

#### Computational Audit (`debug/alpha_audit/`)
- Phase 1: Arithmetic verification, boundary term fix, sensitivity analysis
- Phase 2: Combinatorial search (1.92B formulas, p = 5.2×10⁻⁹)
- Phase 3: Cubic origin (circulant matrix, spectral determinants, heat kernel, topology, QED comparison)
- Phase 4: Derivation chain assessment (Links 1,5 ✓; 2,4 ◐; 3 ✗)

#### Principles Applied
- All mathematical content preserved — only framing, evidence, and structural analysis added
- Epistemic standards match v1.7.2 tightening (Papers 0, 1, 16)
- Paper remains in `papers/conjectures/` (combination rule not derived)

## [1.7.2] - 2026-03-22

### Documentation Review & Epistemic Tightening

Comprehensive review of paper framing and project documentation. No code changes.

#### Papers Revised
- **Paper 0:** Added "What is New vs. What is Known" section; softened ontological claims (removed "Quantum mechanics may not be fundamental. Geometry is."); qualified axiom presentation to acknowledge structural choices beyond stated axioms; added forward reference to Paper 7 confirming energy eigenvalue conjecture
- **Paper 1:** Removed alpha-at-n=5 Platonic solids speculation; fixed Berry phase section (correctly describes log-holonomy, removes retracted relativistic scaling claim from abstract); softened wavefunction-as-statistical-distribution to "one possible interpretation"; cleaned up conclusion to match Paper 7's epistemic standard
- **Paper 16:** Qualified central claim from "derives" to "recasts in representation-theoretic language"; added "What is New vs. What is Known" section; strengthened Hund's rule acknowledgment (Level 1 topology cannot determine ground-state configuration); softened "periodic law IS" to "periodic law corresponds to"

#### Project Documentation
- **CLAUDE.md:** Restructured from ~544 to ~315 lines. Absorbed PAPERS.md. Added "Approaches That Failed" table (7 dead ends), "Natural Geometry Hierarchy" table, "Current Development Frontier" section, and "Project Context" paragraph. Removed emoji headers, directory structure rules, file naming conventions, common commands, documentation workflow, AdS/CFT section, and CLAUDE.md version history.
- **README.md:** Added "Development Methodology" section (AI-augmented workflow disclosure). Added paper acknowledgment template in Citation section.
- **PAPERS.md:** Retired. Content absorbed into CLAUDE.md Section 6.

#### Principles Applied
- All mathematical content preserved exactly — only framing, rhetoric, and interpretation changed
- "What is New vs. What is Known" sections added to Papers 0 and 16 (Paper 7 already had one)
- Downstream papers verified: no references to Papers 0/1 echo removed claims

## [1.7.1] - 2026-03-21

### Paper 2 Rewrite: Fine Structure Constant from Hopf Bundle Geometry
- **Cubic equation:** α³ - Kα + 1 = 0 where K = π(42 + ζ(2) - 1/40)
- **Precision:** 8.8×10⁻⁸ relative error (500× improvement over v1 symplectic approach)
- **Zero free parameters** — all terms derived from spectral geometry
- **Selection principle:** n=3 cutoff from ⟨Casimir⟩/state = dim(S³)

### New Modules
- **hopf_bundle.py**: S³ embedding, Hopf projection, fiber analysis, discrete volume sums
- **Debug scripts**: alpha_eigenvalue_search, alpha_zeta_analysis, alpha_derivation_attempt,
  alpha_fiber_analysis, alpha_precision_analysis, hopf_convergence_analysis

### Key Findings
- 42 = degeneracy-weighted Casimir trace Tr(L²)|_{n≤3}
- ζ(2) = π²/6 = S¹ fiber spectral zeta (single chirality)
- 1/40 = 1/(λ₃ × g₂) = boundary term at n=3, l=2 cutoff
- Cubic structure suggests bare coupling + α² self-interaction (renormalization)

### Paper 2 Superseded
- Old: "Geometric Impedance: Symplectic Framework" (0.0045% error, fitted pitch)
- New: "Spectral Geometry of the Hopf Fibration" (8.8×10⁻⁸ error, no parameters)

## [1.7.0] - 2026-03-21

### Paper 17: Composed Natural Geometries for Core-Valence Diatomics
- **Fiber bundle architecture:** G_total = G_nuc ⋉ G_core(R) ⋉ G_val(R, core_state)
- **LiH:** R_eq = 3.21 bohr (6.4% error, improved from 17% LCAO), ω_e = 1471 cm⁻¹ (4.6%)
- **BeH⁺:** Bound molecule, same code path, zero per-molecule tuning
- **Ab initio pseudopotential:** A from core-valence energy gap, B from ⟨1/r²⟩-weighted core radius
- **Zero experimental molecular data** in pseudopotential derivation

### New Modules
- **CoreScreening** (`core_screening.py`): Z_eff(r) from hyperspherical 2e core solve
- **AbInitioPK** (`ab_initio_pk.py`): Parameter-free Phillips-Kleinman pseudopotential
- **AngularCache / FastAdiabaticPES** (`rho_collapse_cache.py`): ρ-collapsed cache, ~48s PES
- **CoreValenceProjector** (`pauli_projector.py`): Pauli exclusion projection
- **ComposedDiatomicSolver** (`composed_diatomic.py`): General 2+2 diatomic solver

### Performance
- **ρ-collapse:** Angular eigenvalues cached as function of ρ=R/(2R_e), eliminating redundant solves
- **Per-R-point cost:** ~2 seconds (1D tridiagonal radial solve)
- **Full pipeline:** ~4 minutes per molecule on commodity hardware

### Modified
- `level4_multichannel.py`: Z_A_func callable support for screened nuclear charges

### Tests
- 102 new tests across 7 test files (core_screening, rho_collapse, z_eff_injection, lih_composed, composed_diatomic, ab_initio_pk, ab_initio_pk_v2)

## [1.6.1] - 2026-03-18

### Added
- **FrozenCoreLatticeIndex**: Freeze core orbitals, solve active space FCI
  - LiH: 257x SD reduction, 14x speedup, 1.33% error

- **LockedShellMolecule**: Lock complete shells as single states
  - LiH: 2,400x SD reduction, 0.35s runtime, 1.36% error
  - BeH: 4,681x reduction
  - BH: 10,610x reduction
  - Realizes Paper 16 insight: Type C hierarchical structure is computationally separable

### Known Limitations
- Cross-atom integral engine slow for Z > 10 (NaCl blocked by integrals, not SD count)

## [1.6.0] - 2026-03-18

### Chemical Periodicity from Representation Theory

#### Added
- **Paper 16:** "Chemical Periodicity as S_N Representation Theory on S^(3N-1)"
  - Derives periodic table structure from representation theory
  - Universal formula: μ_free = 2(N-2)² for all ground states with S < N/2
  - Theorem 1: λ₁ = 2 universally (spatial Young diagram first-column length)
  - Structure types A/B/C/D/E from Young diagram shape
  - Periodic law as irrep sequence: C → D → E → B within each period
  - Two-level separation: topology (N-dependent) vs metric (Z-dependent)
  - Dirac limit (Z ≈ 137) identified as metric singularity, not topological
  - Extended analysis: H through Og (Z=118), transition metals (Sc-Zn)

- **Debug scripts** for periodic table analysis:
  - `debug_beryllium_periodic.py`: Be [2,2] irrep analysis, H/He/Li/Be comparison
  - `debug_periodic_extended.py`: Na, Mg, Ne, Ar structure types and patterns
  - `debug_transition_metals.py`: Cr/Cu anomalies, Hund's rule at Level 2
  - `debug_superheavy_dirac.py`: Superheavy extrapolation, Dirac limit as metric singularity

#### Theoretical Advances
- S_N representation theory determines chemical periodicity
- Alkali metals (Type C): [N-1,1] irrep — hierarchical, low IE
- Noble gases (Type B): [2,2,...,2] irrep — democratic, high IE
- Per-pair angular momentum contribution → 4 as N → ∞
- Excitation fraction ν/(3N-2) → 1/3 (ground state never squeezed)
- No topological limit to periodic table — ends from nuclear physics (Z ~ 126)
- Relativistic metric correction: δ = (Zα)²/(1-(Zα)²) diverges at Z ~ 137

#### Corrected
- Falsified claim that Hund's rule arises from ν minimization
- Falsified claim that half-filled shell stability is topological
- Cr/Cu anomalies: both configs have same ν = N-2, difference is Level 2 (exchange pairs)
- All corrections applied to `debug_periodic_extended.py` for consistency

## [1.5.0] - 2026-03-18

### Algebraic Structure & SO(3N) Generalization

#### Added
- **Paper 7, Section VI:** "Generalization to N Electrons"
  - SO(3N) isometry group for N-electron atoms on S^(3N-1)
  - General Casimir formula: μ_free = ν(ν + 3N - 2)/2
  - Fermionic statistics determine ground-state irrep
  - Table: H (SO(3)), He (SO(6)), Li (SO(9)), Be (SO(12))

- **Paper 13, Section XII:** "Algebraic Structure of the Angular Problem"
  - A. SO(6) Casimir eigenvalues verified to <0.001%
  - B. Exact first-order formula: a₁ = (8/3π)(√2 - 4Z), match to 9×10⁻⁶
  - C. Selection rules: Δν = 0 mod 4 from S₂ exchange symmetry
  - D. Algebraic/transcendental boundary: a_n algebraic, μ(R) transcendental
  - E. Topology flow S⁵ → S³ × ℝ: IPR, effective dimension, R_c ≈ 4.7 bohr
  - F. Lithium generalization: SO(9) Casimir, [2,1] of S₃, ν=1 ground state
  - G. Emerging pattern: quasi-Coulomb accuracy improves with electron count

- **Debug scripts** for algebraic structure analysis:
  - `debug_so6_casimir.py`: Casimir verification, free eigenvalue tables
  - `debug_algebraic_coefficients.py`: Matrix element computation
  - `debug_a2_final.py`: Second-order perturbation convergence
  - `debug_closed_form_final.py`: Partial harmonic sum derivation (definitive)
  - `debug_crossover_analysis.py`: Padé, interpolation, continued fractions
  - `debug_curvature_correlation.py`: Conformal factor and Ricci flow tests
  - `debug_topological_transition.py`: S⁵ → S³ × ℝ surgery characterization
  - `debug_flow_comparison.py`: Symplectic impedance vs topology flow
  - `debug_lithium_algebraic.py`: 3-electron SO(9) analysis

#### Theoretical Advances
- Perturbation coefficients are individually algebraic (rational × π^{-n})
- Full eigenvalues are transcendental (topology flow integral)
- Quasi-Coulomb accuracy improves with electron count:
  He (2e): 5.5% error → Li (3e): ~0% error
- Fermionic statistics as mechanism: higher irreps have larger μ_free
- Surgery correction Δ = −0.16 Ha is dynamical, not algebraic

#### Changed
- Paper 7: Added Section VI before Discussion (~1.5 pages)
- Paper 13: Added Section XII with 7 subsections (~3 pages)
- Paper 13: Updated "Open questions" to reference lithium results

## [1.4.0] - 2026-03-17

### Papers 14-15, Level 4 Solver, Heteronuclear Extension

#### Added
- Paper 14: Structurally sparse qubit Hamiltonians (O(Q^3.15) Pauli scaling)
- Paper 15: Level 4 natural geometry for two-center two-electron molecules
  - H₂: 94.1% D_e with variational 2D solver (29 channels, σ+π)
  - HeH⁺: 93.1% D_e with charge-center origin (57 channels, σ+π)
- `geovac/level4_multichannel.py` — Molecule-frame hyperspherical solver
- Algebraic nuclear coupling via multipole expansion (replaces quadrature)
- π-channel support (m_max parameter) for angular correlation
- Heteronuclear molecules (Z_A, Z_B parameters, charge-center origin)
- Variational 2D solver with shift-invert eigensolver
- Wigner 3j symbols for generalized angular coupling
- `tests/test_level4_multichannel.py` — Convergence and heteronuclear tests

#### Changed
- Nuclear coupling: 1019× speedup from algebraic vs quadrature
- Channel enumeration: extended to (l₁, m₁, l₂, m₂) tuples

#### Fixed
- 2D solver ghost modes at large matrix sizes (shift-invert with sigma)
- Adiabatic approximation diagnosed (~11% systematic overestimate)

## [1.2.0] - 2026-03-15

### The QA & Corrections Release

#### Paper 1: Berry Phase Retraction
- Section IV rewritten: "Geometric Curvature" -> "Geometric Phase Structure of the Lattice"
- Retracted k = 2.113 claim (arg() = 0 for real SU(2)/SU(1,1) operators)
- Documented log-holonomy as valid quantity (k = 1.0 exact)
- Updated Discussion, Conclusion, and Appendix A/B references
- Added erratum noting the retraction

#### Paper 13: Autoionization & Adiabatic Limits
- New Section X: Autoionization channel classification from angular topology
- New Section XI: Limits of the adiabatic approximation for quantitative widths
- Updated Discussion and Conclusion with new open questions

#### New Modules
- `geovac/berry_phase.py` — Log-holonomy computation on geometric lattice plaquettes
- `geovac/hyperspherical_complex_scaling.py` — Exterior complex scaling for resonances
- `geovac/hyperspherical_coupling.py` — Coupled-channel adiabatic solver
- `geovac/hyperspherical_resonances.py` — Resonance detection and analysis

#### QA Sprint
- `debug/qa_sprint/berry_phase_reconciliation.md` — Full Berry phase investigation
- `debug/qa_sprint/benchmark_audit.md` — Verified He energy, Neumann D_e, kappa, symbolic proofs
- `debug/qa_sprint/test_health.md` — Test suite health check (528 pass, 0 fail, 1 xfail)
- `debug/qa_sprint/cross_document_consistency.md` — 7/8 claims consistent; Berry phase sole issue (now corrected)

#### Test Fixes
- `tests/test_lih_validation.py`: Added missing assertions to 2 tests
- `tests/test_hylleraas.py`: Made `test_energy_below_atoms` unconditional (enlarged basis to j_max=2)
- `benchmarks/ab_initio_nuclear/results.md`: Annotated stale prolate CI numbers vs Paper 13 Hylleraas PES

#### Test Suite
- 528 passed, 0 failed, 1 xfailed (Sturmian structural theorem — expected)

## [1.1.0] - 2026-03-15

### The Multi-Particle Release

#### Paper 12: Algebraic Two-Electron Integrals
- `geovac/neumann_vee.py` — Neumann expansion of 1/r_12 in prolate spheroidal
  coordinates: replaces 5D numerical quadrature with sums of 1D algebraic integrals
- Auxiliary integral tables: C_l (eta moments), A_l (xi exponential moments),
  B_l (second-kind Legendre), X_l (2D xi integrals via integration by parts)
- H2 at R=1.4011 bohr: 92.4% D_e with 27 basis functions (vs 79.6% numerical)
- Uniform 12-20 percentage point improvement over numerical V_ee at all basis sizes
- Remaining 7.6% gap diagnosed as electron-electron cusp (basis completeness limit)
- `tests/test_neumann_vee.py` — convergence and accuracy tests

#### Paper 13: Hyperspherical Lattice for Two-Electron Atoms
- `geovac/hyperspherical_angular.py` — Angular eigenvalue solver:
  - Liouville substitution u = sin(a)cos(a)f for standard Schrodinger form
  - FD discretization with Dirichlet BCs
  - Gaunt integral coupling via Wigner 3j symbols
  - L=0 partial-wave reduction for He ground state
- `geovac/hyperspherical_adiabatic.py` — Adiabatic potential curves:
  - mu(R) computed on R grid, V_eff(R) = mu/R^2 + 15/(8R^2)
  - Asymptotic: V_eff -> -Z^2/2 (He+ threshold)
- `geovac/hyperspherical_radial.py` — Hyperradial solver:
  - Self-adjoint FD + cubic spline interpolation
  - solve_helium() orchestrator: angular -> adiabatic -> radial
- **He ground state:** E = -2.9052 Ha (0.05% error vs exact -2.9037 Ha)
- Single-channel adiabatic, 600x600 sparse matrix, ~3 seconds runtime
- First non-trivial fiber bundle in GeoVac (angular structure varies with R)
- `tests/test_hyperspherical_he.py` — 20 tests (Gaunt, angular, adiabatic, solver)

#### Ab Initio Molecular Spectroscopy
- Full pipeline: electron lattice -> PES -> Morse fit -> nuclear lattice -> spectrum
- H2 Morse parameters: R_e = 1.42 bohr (+1.2%), omega_e = 4435 cm-1 (+0.8%)
- Rovibrational transitions: nu_01 = 4157 cm-1 (-0.1%), B_e = 59.5 cm-1 (-2.2%)
- Zero experimental spectroscopic input at any stage
- `benchmarks/ab_initio_nuclear/` — full pipeline scripts and results

#### Qubit Encoding Benchmarks
- GeoVac Pauli terms: O(Q^3.15) vs conventional O(Q^4.60)
- ERI density: 10% (nmax=2) -> 1.3% (nmax=5), drops as ~1/M^2
- Structural advantage from angular momentum selection rules (Wigner 3j)
- `benchmarks/qubit_encoding/` — comparison scripts and results

#### Zenodo Publication Cluster
- `papers/zenodo/` — Papers 7, 11, 12, 13 prepared for upload
- `.zenodo.json`, `CITATION.cff` — metadata files
- All four papers compile cleanly with pdflatex
- Cross-reference audit: all inter-paper references resolved

#### Documentation
- Paper 13 written: `papers/core/paper_13_hyperspherical.tex` (8 pages, revtex4-2)
- README updated: v1.1.0 benchmarks, paper series, quick start examples
- CLAUDE.md updated: Paper 12/13 in paper list, new modules, fiber bundle concept

## [1.0.1] - 2026-03-13

### Prolate Spheroidal Stress Tests Complete

#### Two-Electron H2 CI
- `geovac/prolate_scf.py` — Eckart SCF, grid SCF, relaxed-orbital CI
- V_ee integrals via azimuthal averaging with complete elliptic integral K(k)
- **Frozen CI:** D_e = 0.077 Ha (44% of exact), R_eq = 1.77 bohr
- **Relaxed CI (Eckart Z_eff):** D_e = 0.101 Ha (58% of exact), R_eq = 1.40 bohr
- Pi orbital selection rule proven: <sigma^2|pi^2> = 0 (fundamental symmetry)
- 4-sigma CI (6x6) adds only 0.1 mHa — confirms minimal basis is the bottleneck

#### Heteronuclear Extension
- `geovac/prolate_heteronuclear_scf.py` — Per-atom Z_eff optimization
- HeH+ binding recovered with independent Z_eff_A, Z_eff_B
- Frozen HeH2+ orbitals unbound (J_11 ~ 1.3 Ha); optimized Z_eff restores binding

#### Grid-Based SCF
- `geovac/prolate_scf.py` — 2D finite-difference Fock solver (proof of concept)
- Non-uniform xi grid via quadratic transformation (critical for accuracy)
- H2 binding confirmed on 2D FD grid

#### Tests
- 95/95 stress tests pass across 8 phases
- New test files: test_prolate_stress.py, test_prolate_h2_4sigma.py,
  test_prolate_heh_plus.py, test_prolate_scf.py, test_prolate_relaxed_ci.py,
  test_prolate_heteronuclear_scf.py, test_prolate_grid_scf.py

#### Documentation
- Paper 11 Sec. VIII.D updated: H2 CI results, pi selection rule, HeH+ extension
- README updated: H2 prolate CI benchmark, scope corrections
- `debug/PROLATE_STRESS_TESTS.md` — Phases 6-8 documented

## [1.0.0] - 2026-03-12

### The Natural Geometry Release

v1.0.0 marks the completion of the GeoVac theoretical framework: from atoms
on S3 (Fock projection) to molecules on prolate spheroids (molecular Fock
projection), unified by the natural geometry principle.

#### Prolate Spheroidal Lattice (Paper 11)
- `geovac/prolate_spheroidal_lattice.py` — `ProlateSpheroidalLattice` class:
  separated prolate spheroidal solver for one-electron diatomics.
- Spectral angular solver (Legendre basis, machine precision) + self-adjoint
  FD radial solver + Brent root-finding in c^2.
- **H2+ results:** R_eq = 2.001 bohr (0.21% error), E_min = -0.5984 Ha
  (0.70% error), correct dissociation to H + p. Zero free parameters.
- HeH2+ heteronuclear extension validated.
- 11 tests in `tests/test_prolate_h2plus.py` (all passing).

#### Nuclear Lattice (Paper 10)
- `geovac/nuclear_lattice.py` — graph-based vibrational/rotational solver.
- `geovac/coupled_en_lattice.py` — coupled electronic-nuclear framework.
- Tests: `tests/test_nuclear_lattice.py`, `tests/test_coupled_en.py`.

#### Paper Series (11 papers)
- Paper 11 written: "The Molecular Fock Projection: Diatomic Molecules as
  Prolate Spheroidal Lattices" (7 pages, revtex4-2).
- Paper 10 written: "The Nuclear Lattice: Discrete Graph Structures for
  Molecular Vibration and Rotation" (6 pages).
- Forward references to Paper 11 added to Papers 7, 8-9, FCI, LCAO-FCI, 10.
- All 6 updated/new papers compile cleanly (pdflatex, no errors).
- `papers/PAPER_STATUS.md` updated with full inventory.

#### Cleanup
- `run_all.py` moved from root to `benchmarks/scripts/`.
- `figures/` contents moved to `debug/plots/` and `debug/data/`.
- `nul` (Windows artifact) removed.
- `setup.py` description updated, version synced.
- README.md rewritten for v1.0.0: updated benchmarks, project structure,
  paper series table, prolate spheroidal quick-start example.

#### Unchanged
- All production code in geovac/ (lattice.py, hamiltonian.py, lattice_index.py,
  direct_ci.py, etc.) unchanged from v0.9.37.
- All existing tests pass.

## [0.9.0 – 0.9.37] - 2026-02-22 to 2026-03-11

### LCAO/Graph-Concatenation Diagnostic Arc (29 versions, superseded)

This version range represents the LCAO molecular diagnostic arc — a systematic
29-version investigation of graph-concatenation approaches to molecular electronic
structure. The arc is **superseded** by the natural geometry framework (Papers 11,
13, 15, 17) starting at v1.0.0.

#### Key Results
- **Root cause identified:** Graph Laplacian kinetic energy is R-independent in the
  LCAO basis, producing monotonically attractive PES with no equilibrium geometry.
  The Mulliken cross-nuclear diagonal is too strong at short R.
- **LiH benchmark (best):** D_e_CP = 0.083 Ha (10% error) at R=3.015 via
  Boys-Bernardi counterpoise correction (v0.9.9–0.9.11). R_eq ~2.5 bohr
  (expt 3.015) — a discretization artifact.
- **Sturmian arc (v0.9.20–0.9.34):** 15 versions testing Sturmian alternatives.
  Structural theorem proven: Sturmian H proportional to S, eigenvalues
  geometry-independent. All variants failed. See `docs/sturmian_attempts.md`.
- **Orbital exponent relaxation (v0.9.36):** Shifts PES uniformly, not
  differentially. R_eq unchanged. Not a mechanism for equilibrium geometry.

#### What Was Built (still in codebase)
- `MolecularLatticeIndex` — heteronuclear FCI with separate hydrogenic lattices
- `DirectCISolver` — excitation-driven sparse FCI (O(N_SD^<2) scaling)
- `compute_bsse_correction()` — Boys-Bernardi counterpoise
- `wigner_so4.py` — SO(4) Wigner D-matrix elements
- Slater F0 V_ee with disk cache, hybrid h1 method
- 3-electron and 4-electron FCI (Li 1.07%, Be 0.90%)

#### Dead Ends (full list)
See **CLAUDE.md Section 3 ("Approaches That Failed")** for the complete table
of negative results from this period, including LCAO graph concatenation,
Sturmian basis, orbital relaxation, Mulliken cross-nuclear, and Berry phase.

#### Final State (v0.9.37)
Codebase locked. Documentation-only release. All production code frozen pending
transition to natural geometry at v1.0.0.

---

#### Paper 8 Narrative Complete
- Added Sturmian Structural Theorem (Sec. XII): H proportional to S,
  eigenvalues geometry-independent. Binding requires beta_k(R) from PDE.
- Added SO(4) Selection Rules (Sec. XIII): D2_10_10 = 1 and D2_00_10 = 0
  for all gamma. The 2p0 orbital is a transparent mode of the bond channel.
- Added Harmonic Phase Locking negative result (Sec. XIV): no D-matrix
  critical points map to R_eq ~ 3.015 bohr. R_eq is emergent variational
  balance, not geometric fixed point.
- Revised abstract and conclusion to reflect three classes of results:
  positive structural theorems, negative Sturmian theorem, negative
  phase locking observation.

#### Codebase Cleanup
- Created `debug/archive/` and moved 31 one-off diagnostic scripts from
  the v0.9.20-v0.9.36 experimental arc.
- Created `docs/sturmian_attempts.md` — complete narrative of the Sturmian
  arc (v0.9.20-v0.9.34) with structural theorem summary.

#### New Paper Scaffold
- Created `papers/core/paper_geovac_lcao_fci.tex` — scaffold for the
  LCAO FCI results paper with section structure, key claims, and
  bibliography. Working title: "Topological Full Configuration Interaction
  for Heteronuclear Diatomics: LiH as a Benchmark."

#### Scientific Results (Final, nmax=3)
- D_e_CP = 0.093 Ha (1.0% error vs expt 0.092 Ha)
- R_eq ~ 2.5 bohr, identified as nmax discretization artifact
- Sturmian structural theorem proven (29 diagnostic versions)
- All specific R_eq hypotheses tested and eliminated

#### Unchanged
- All production code in geovac/ untouched.
- All tests pass.


## [0.1.0 – 0.8.0] - 2026-02-01 to 2026-02-21

### Foundation Era (pre-natural-geometry)

Early development of the GeoVac framework before the natural geometry principle
was established. These versions built the core infrastructure that later versions
depend on, but the molecular approaches from this era are superseded.

#### Key Milestones
- **v0.1.0** (Feb 1): Initial release — `GeometricLattice`, `HeliumHamiltonian`,
  graph Laplacian kinetic energy, sparse eigenvalue solver.
- **v0.2.0–0.2.1** (Feb 12–13): Discovery of universal constant kappa = -1/16.
  H2+ control experiment (0% error). Classification as topological Hartree-Fock.
  `MoleculeHamiltonian` class for LCAO bonding.
- **v0.3.1–0.3.2** (Feb 13–14): Multi-solver architecture (mean-field, DFT, FCI,
  Dirac). `AtomicSolver` class. Z^2-scaling for hydrogenic ions.
- **v0.4.0–0.4.2** (Feb 15): Isoelectronic scaling breakthrough (Three Laws:
  conformal, Coulomb, torsion). Li+ 0.03%, Be2+ 0.15% error.
- **v0.5.0** (Feb 15): Distance-dependent bridges. Alpha-metric bond constants.
  H2 R_eq = 1.40 bohr (experimental match with fitted parameters).
- **v0.6.0** (Feb 15): Schwarzschild torsion metric. Full periodic table coverage
  (Z=1 to Z=92+).
- **v0.7.0** (Feb 15): Quantum dynamics engine (`TimePropagator`). Rabi oscillation
  validation.
- **v0.8.0** (Feb 21): Dynamics & thermodynamics release. Delta-kick spectroscopy,
  PES scanning, AIMD, Langevin thermostat. Paper 6 written.

#### Infrastructure Still in Use
- `GeometricLattice`, `GraphHamiltonian`, `LatticeIndex`, `DirectCISolver`
- Sparse eigenvalue solver, Slater F0 integrals, disk cache
- `TimePropagator` for quantum dynamics (Paper 6)
- 18 symbolic S3 proofs (`test_fock_projection.py`, `test_fock_laplacian.py`)

---

## Version Naming Convention

- **Major (X.0.0):** Architectural changes (e.g., LCAO -> natural geometry)
- **Minor (0.X.0):** New features, completed diagnostic arcs, paper additions
- **Patch (0.0.X):** Bug fixes, documentation updates
