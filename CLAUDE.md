# GeoVac Project Guidelines

## 1. Project Identity

**Name:** GeoVac (The Geometric Vacuum)
**Version:** v2.0.11 (March 29, 2026)
**Mission:** Spectral graph theory approach to computational quantum chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (unit S3) that is mathematically equivalent to the Schrodinger equation via Fock's 1935 conformal projection. This equivalence is exploited computationally to replace expensive continuous integration with O(N) sparse matrix eigenvalue problems.

**Authoritative source rule:** The papers in `papers/core/` and `papers/observations/` are the authoritative source for all physics. If any documentation (README, CHANGELOG, code comments) conflicts with the papers, the papers win. Flag the conflict to the user rather than silently resolving it.

**Project context:** GeoVac is an independent research project with no institutional affiliation, developed using an AI-augmented agentic workflow. The principal investigator provides scientific direction and quality control; implementation and documentation drafting are performed collaboratively with LLMs (Anthropic Claude). The primary dissemination channel is GitHub + Zenodo (DOI-stamped releases). The papers in `papers/core/` and `papers/observations/` are written to academic standards but are not submitted to traditional journals. The project's viability case rests on producing a usable, benchmarked computational tool. Do not suggest formatting papers for specific journals or pursuing traditional peer review unless the user asks.

---

## 1.5. Positioning & Framing

GeoVac is a discretization framework that exploits the natural geometry of separable quantum systems to produce sparse graph Hamiltonians. It is **not** a new foundation for quantum mechanics — it works because the graph Laplacian converges to the known continuous Laplace-Beltrami operators in the continuum limit (Paper 7).

**Rhetoric rule:** The GeoVac framework is demonstrably conformally equivalent to standard continuous quantum mechanics at every level where it has been tested (Paper 7 provides 18 symbolic proofs for the S³ case; subsequent papers verify operator convergence for each natural geometry). The papers should present this conformal equivalence as the primary result. The discrete graph topology and the continuous Laplace–Beltrami operators are dual descriptions connected by proven conformal maps — the mathematics supports both readings, and the interpretation of which is more fundamental is left as a choice for the reader. Avoid language that asserts ontological priority of either description. Lead with concrete computational results (sparsity, scaling, accuracy, structural insight) rather than interpretive claims about the nature of quantum mechanics.

**Lead with concrete advantages:** O(V) sparsity, angular momentum selection rules baked into the basis, zero-parameter construction from nuclear charges and geometry alone, efficient qubit encodings (Paper 14). These structural properties are the framework's actual selling points — not philosophical claims about the nature of quantum mechanics.

**Benchmarking rule:** When comparing to other methods, always use the strongest available baseline (cc-pVTZ or better for atoms, explicitly correlated methods for molecules), not just STO-3G. If the comparison is unfavorable, say so honestly and identify what the framework offers instead (sparsity, scaling, structural insight). Position the framework as a computationally principled alternative, not as a replacement for production quantum chemistry.

---

## 2. Current Development Frontier

**Best results by system type:**
- Atoms: He at 0.05% error (hyperspherical coordinates, Paper 13)
- One-electron molecules: H2+ at 0.0002% (spectral Laguerre, Paper 11)
- Two-electron molecules: H2 at 94.1% D_e (Level 4 mol-frame hyperspherical, Paper 15)
- Core-valence molecules: LiH R_eq 5.3% error with l-dependent PK at l_max=2 (composed geometry, Paper 17)
- Polyatomic molecules: BeH₂ R_eq 11.7% error (full 1-RDM exchange, composed geometry, Paper 17)
- Triatomic molecules: H₂O R_eq 26% error, uncoupled with charge-center origin (composed geometry, 5-block Level 3+4, zero parameters, Paper 17)
- Qubit encoding: O(Q^3.15) Pauli terms, O(Q^3.36) QWC groups, O(Q^1.69) 1-norm (Paper 14)
- VQE benchmark: 1.3x fewer Pauli terms at Q=10, 8.1x at Q=28 vs validated Gaussian cc-pVDZ/cc-pVTZ (Paper 14, v1.9.0)
- Composed qubit encoding: universal Q^2.5 Pauli scaling across LiH/BeH2/H2O (exponent spread 0.02), 51x-1,712x advantage over published Gaussian baselines (Paper 14, v2.0.0)
- Fine structure constant: alpha from Hopf bundle at 8.8x10^-8, zero free parameters, p-value 5.2x10^-9, universal algebraic identity B_formal/N = d, Hopf generalization negative result, circulant Hermiticity, second selection principle (Paper 2, conjectural)

**Active frontier:**
- l-dependent PK pseudopotential for LiH composed geometry: reduces R_eq error from 6.9% to 5.3% at l_max=2 (higher l_max with channel-blind PK diverges — see Section 3). Remaining l_max divergence is from differential angular correlation. Diagnostic arc (v2.0.5) confirmed: divergence is linear (+0.23 bohr/l_max), angular spreading is 100% nuclear, Level 4 solver converges (β=0.42/l_max). R-dependent PK w_PK(R) = δ_{l,0} × min(cap, R/R_ref) achieves 2.0% at l_max=4 but requires empirical R_ref; ab initio R_ref derivation is the next milestone.
- Polyatomic PES: BeH2 full 1-RDM exchange achieves 11.7% R_eq error (down from 20% diagonal S·F⁰, matches 12% fitted model with zero free parameters). Kinetic orthogonalization tested: +0.03 Ha uniform, no R_eq effect (negative result). Remaining 11.7% attributed to basis truncation (l_max=2) and adiabatic approximation. Next steps: l_max convergence with full exchange, non-adiabatic corrections, cusp corrections.
- H₂O composed solver: 5-block architecture (O 1s² core + 2 O–H bond pairs + 2 lone pairs), R_eq 26% error with charge-center origin. Coupling architecture validated (bond-bond ~0.5 Ha, consistent with BeH₂), but lone pair coupling unphysical at Z_eff=6 (S·F⁰ produces −28 Ha bond-lone, −15 Ha lone-lone — exceeds total electronic energy). Bottleneck is Level 4 angular basis at 6:1 charge asymmetry (R_eq error scales: 0.1% at 1:1, 18% at 2:1, 17% at 6:1 with charge-center). Next step: improved angular basis or non-uniform grid for highly asymmetric bonds.
- Algebraic audit confirmed: Level 4 nuclear coupling (Paper 15) uses split-region Legendre expansion with Gaunt integrals — not quadrature. Core screening (Paper 17) has algebraic density from channel coefficients. Paper 15 Section V.A corrected to match production code.
- Improving Level 4 D_e recovery beyond 94% (channel convergence, cusp corrections)
- Chemical periodicity as representation theory (Paper 16) -- computational exploitation of hierarchical structure
- Quantum simulation cost: equal-qubit Gaussian comparison validated with computed integrals (v1.9.0); commutator-based Trotter bounds: r ~ Q^{1.47} (vs Q^{1.69} 1-norm), 7x fewer steps at Q=60 (Paper 14, v2.0.3)
- Track B — Hyperradial algebraicization: COMPLETE. Algebraic spectral basis replaces FD angular solver (matrix dimension 200-800 → 10-30). Exact Hellmann-Feynman coupling matrices (R-independent dH/dR). Single-channel adiabatic approximation identified as bottleneck at l_max≥1 (energy dives below exact). Coupled-channel integration fixes convergence direction: error decreases monotonically with l_max (0.37% → 0.27% at l_max 1-3). Q-matrix closure overcorrects at l_max=0 (1.1%); improvement pathway identified (exact dP/dR). Sturmian variant investigated — modest improvement at l_max=0, counterproductive at l_max≥1 (converges faster to wrong adiabatic limit). Preserved as research artifact. (Paper 13, v2.0.6)
- Track A — l_max divergence: 2D SOLVER INTEGRATED, PARTIAL FIX. Variational 2D solver (Paper 15) wired into composed LiH pipeline (v2.0.9). Reduces l_max drift by 4× (from +0.400 to +0.100 bohr/l_max). 75% of divergence is from adiabatic approximation (as diagnosed v2.0.6); 25% residual survives the 2D solver — likely PK-related or intrinsic to composed-geometry separation. Adiabatic still wins at l_max=2 (2.8% vs 6.1% R_eq error) due to error cancellation; 2D wins at l_max=3 (9.5% vs 16.1%). Default remains `level4_method='adiabatic'`. Next investigation: test with `pk_mode='none'` to isolate residual drift source. (v2.0.9)
- Track C — Level 2 algebraicization (prolate spheroidal): COMPLETE. Spectral Laguerre basis replaces 5000-point FD radial ξ-solver. 250× dimension reduction, 270× speedup, 5000× accuracy improvement (1.01% → 0.0002%). Production wiring complete: `scan_h2plus_pes(radial_method='spectral')` gives 287× speedup with 5000× E_min accuracy (v2.0.9). Algebraic Laguerre matrix elements via three-term recurrence eliminate ALL quadrature for m=0 σ states — the Level 2 radial solver is now fully algebraic (`matrix_method='algebraic'`, machine-precision agreement with quadrature). (v2.0.9)
- Track D — Level 3 exact Q-matrix: COMPLETE. Algebraic dP/dR derived from Hellmann-Feynman quantities (verified vs FD to <1e-9, zero finite differences needed). Exact Q = PP + dP/dR replaces closure approximation. Results: 14-19% improvement over closure at l_max≥1, best coupled-channel error 0.22% at l_max=3. l_max=0 overcorrection (1.1%) is structural — dominated by diagonal DBOC, not fixable by Q-matrix improvement alone. q_mode='exact' recommended for l_max≥1 production use. (v2.0.6)
- Track E — Level 3 coupled-channel convergence ceiling: COMPLETE. Extended to l_max=5 with q_mode='exact'. Error converges algebraically (~l_max⁻²) to structural floor of 0.19-0.20%, set by adiabatic channel truncation. Sub-0.1% requires variational 2D solver (Paper 15). 3 channels sufficient (4-5 contribute < 0.001%). n_basis and radial grid fully converged. (v2.0.8)
- Track F — Level 3 spectral hyperradial: COMPLETE. Spectral Laguerre basis replaces 3000-point FD hyperradial grid. 120× dimension reduction (3000 → 25 basis functions, converged at n_basis=15), 95× coupled-channel speedup (561 ms → 5.9 ms). Spectral-FD consistency < 0.00003 Ha. Coupled-channel ceiling unchanged (0.221% at l_max=3 — same as FD's 0.220%). Alpha-insensitive. Default preserved (`radial_method='fd'`). (v2.0.9)
- Track G — Level 3 perturbation series: COMPLETE (definitive negative result). Rayleigh-Schrödinger perturbation series for μ(R) exploiting linear matrix pencil H(R) = H₀ + R·V_C. a₁ validated to 10⁻¹⁵ against Paper 13 Table II. Raw series converges for R < 2 bohr; Padé extends to R ≈ 3 bohr; all orders fail beyond R ≈ 5 bohr. Root cause: μ(R) is transcendental (O(R) → O(R²) regime transition). Point-by-point R-grid diagonalization cannot be eliminated globally. Series provides exact R=0 derivatives for spline seeding. Confirms Paper 13 Sec XII.B. (v2.0.9)
- Track H — Level 3 algebraic Laguerre matrix elements: COMPLETE. Overlap S and kinetic K fully algebraic via three-term Laguerre recurrence (pentadiagonal M2 for S, tridiagonal derivative expansion for K). Machine-precision agreement (< 1e-14 relative error). 11× matrix build speedup (52 μs vs 585 μs). V_eff stays quadrature (transcendental, Track G). `matrix_method='algebraic'` in spectral solver. (v2.0.10)
- Track I — Level 4 spectral Laguerre hyperradial: COMPLETE. Spectral Laguerre basis replaces FD R_e grid for all three Level 4 pathways (adiabatic, coupled-channel, 2D variational). 16× dimension reduction (400 → 25), < 0.0005 Ha FD agreement (adiabatic). No wall time speedup — angular sweep dominates 99% of Level 4 cost (structural finding: FD radial is already O(N) tridiagonal, sub-ms). Spectral value is accuracy, memory, and parameterization reduction. n_basis=20 optimal; mild conditioning at n_basis≥25. `radial_method='spectral'` in Level 4 solver. (v2.0.10)
- Track J — Level 2 algebraic m≠0 (associated Laguerre): COMPLETE. Associated Laguerre basis L_n^{|m|}(x) with weight x^|m|·e^{-x} absorbs the 1/x centrifugal singularity. Partial-fraction decomposition splits centrifugal term into lowered moment M_{-1} (algebraic, DLMF 18.9.13 summation identity) plus Stieltjes integral J (three-term recurrence seeded by e^a·E₁(a)). All matrix elements algebraic except single transcendental seed. Associated basis converges faster than ordinary Laguerre for m=1 (stable by N=10 vs non-monotonic at N=40). m=2 delta states: algebraic-quadrature agreement ~5e-9. PES shape preserved (R_eq matches). `matrix_method='algebraic'` now works for ALL m values. 119 tests passing (77 associated Laguerre + 42 kinetic). Paper 11 Sec V.D updated. (v2.0.10)

**Backlog:**
- Q-matrix improvement for Level 3 coupled-channel — DONE (v2.0.6)
- Track C implementation + production wiring + algebraic matrix elements — DONE (v2.0.9)
- Track F spectral hyperradial — DONE (v2.0.9)
- Track G perturbation series — DONE, definitive negative (v2.0.9)
- Rebuild composed-geometry Hamiltonians with spectral Level 2 solver (expect accuracy improvement from elimination of FD error)
- Track A next step: test composed LiH with `pk_mode='none'` to isolate residual 25% l_max drift source
- Apply algebraic Laguerre matrix elements to Level 3 hyperradial — DONE (Track H, v2.0.10). Overlap S and kinetic K fully algebraic via three-term Laguerre recurrence (< 1e-14 relative error, 11× build speedup). Potential V_eff stays quadrature (transcendental, Track G). Energy agreement < 1e-14 Ha. `matrix_method='algebraic'` in spectral solver.
- Apply spectral Laguerre to Level 4 hyperradial — DONE (Track I, v2.0.10). Three pathways (adiabatic, coupled-channel, 2D) all wired. 16× dimension reduction, < 0.0005 Ha FD agreement. No wall time speedup — angular sweep dominates 99% of Level 4 cost (structural finding). `radial_method='spectral'` in Level 4 solver.
- Angular sweep caching/acceleration for Level 4 (identified as 99% bottleneck, Track I)
- Dense spectral PES scan (200+ R-points) for precision H₂⁺ spectroscopic constants (R_eq, ω_e, B_e, ν₀₁)

**Architecture locked:** The LCAO/graph-concatenation approach (v0.9.x series) is superseded. All molecular work uses natural geometry (Papers 11, 13, 15, 17).

---

## 3. Approaches That Failed

Critical institutional memory. Do not re-derive these dead ends.

| Approach | Why It Fails | Resolution | Reference |
|:---------|:-------------|:-----------|:----------|
| LCAO graph concatenation for molecules | Graph Laplacian kinetic energy is R-independent -> monotonically attractive PES, no equilibrium | Prolate spheroidal lattice (Paper 11) or composed geometry (Paper 17) | FCI-M, 29-version diagnostic arc |
| Sturmian basis with shared p0 | H proportional to S -> eigenvalues R-independent | Prolate spheroidal separation introduces R-dependence via beta_k(R) | Papers 8-9, Structural Theorem |
| Berry phase from lattice plaquettes | arg() = 0 identically for real SU(2)/SU(1,1) operators | Log-holonomy Theta(n) = -2ln((n+1)/n) ~ 1/n is the valid geometric quantity | Paper 1 v1.2.0 erratum |
| Numerical V_ee quadrature on prolate spheroid | Coulomb singularity at r12=0 causes slow convergence, saturates at ~80% D_e | Algebraic Neumann expansion eliminates quadrature entirely | Paper 12 |
| Single-S3 molecular encoding | One p0 cannot encode R-dependent bonding physics | Natural geometry principle: use coordinate system where separation occurs | Papers 8-9, 11 |
| Orbital exponent relaxation (zeta) | Shifts PES uniformly, not differentially; R_eq unchanged | Not a mechanism for equilibrium geometry | v0.9.36 |
| Mulliken cross-nuclear diagonal | Too strong at short R, drives R_eq inward | Bond sphere (Paper 8) or natural geometry approach | v0.9.35 diagnosis |
| Z_eff partitioning for polyatomics | Dividing screened charge among bond pairs weakens bonds without adding inter-bond repulsion; R_eq goes from 32% to 103% error | Inter-bond orbital coupling needed (exchange/orthogonalization), not charge redistribution | v2.0.0 BeH2 diagnostic |
| Classical inter-bond repulsion | Point-charge V_inter(R) has wrong R-dependence; pushes R_eq outward (wrong direction) at all scaling factors | Orbital-level inter-bond coupling required; scalar energy corrections cannot capture wavefunction modification | v2.0.0 BeH2 diagnostic |
| Higher l_max / sigma+pi channels for LiH composed | PK pseudopotential is channel-blind (derived from l=0 core); higher-l channels add correlation preferentially at large R, pushing R_eq outward monotonically (l2: 6.4%, l3: 16.9%, l4: 32.8%) | l-dependent or m-dependent PK pseudopotential needed | Paper 17 Sec VI.A, v2.0.1 |
| Monopole (F⁰) inter-fiber Coulomb coupling for BeH₂ | Monopole energy ~1.6 Ha, nearly R-independent near equilibrium; uniform repulsion pushes R_eq outward (32% → 62% error) | Exchange coupling via inter-fiber channel overlap S(R) — has correct R-dependence (20% R_eq error) | v2.0.1 |
| Higher Slater multipoles (k≥1) for BeH₂ inter-fiber exchange | k=1 dipole is repulsive (worsens R_eq); k=2 quadrupole is negligibly small; total k=0+1+2 gives 20.5% error (worse than monopole alone at 20.0%) | Missing exchange is from off-diagonal 1-RDM terms, not higher angular multipoles | v2.0.1 |
| Kinetic orthogonalization for BeH₂ | Löwdin correction is +0.03 Ha, nearly R-independent (peaks at R~4.0 where l≥1 channel growth counteracts overlap decay); no effect on R_eq | Residual 11.7% error is from basis truncation (l_max=2) and adiabatic approximation, not missing orthogonalization | v2.0.2 |
| Lone pair inter-fiber coupling at Z_eff=6 (H₂O) | S·F⁰ produces −28 Ha bond-lone and −15 Ha lone-lone coupling — exceeds total electronic energy. Slater integrals scale as ~Z, unphysical at Z_eff=6 (validated only at Z_eff=2 for BeH₂) | Disable lone pair coupling; bond-bond only is validated (~0.5 Ha, consistent with BeH₂) | v2.0.4 H₂O diagnostic |
| Midpoint origin for asymmetric bond pairs (Z_A ≫ Z_B) | Default midpoint origin places coordinate center far from electron density; R_eq error 56% at Z_A/Z_B=6:1 vs 17% with charge_center | Use origin='charge_center' with z0 = R(Z_A−Z_B)/(2(Z_A+Z_B)) | v2.0.4 H₂O diagnostic |
| Eigenchannel rotation for l_max compression | Rotation to eigenchannels of coupling matrix does not reduce number of channels needed; all channels contribute at large R | Divergence is not from channel mixing but from PK coupling to higher-l channels | v2.0.5 diagnostic arc |
| Prolate spheroidal basis compression | Spheroidal harmonics at 6:1 aspect ratio require 85% of spherical channels to capture 99% weight | Spheroidal basis provides no compression for asymmetric systems | v2.0.5 diagnostic arc |
| Self-consistent PK iteration | Reduces drift rate by ~50% but does not eliminate l_max divergence; intermediate performance | R-dependent PK scaling is more effective | v2.0.5 diagnostic arc |
| Projected PK (one-shot from core density) | PES collapses — no equilibrium geometry | PK must retain variational character; projection destroys it | v2.0.5 diagnostic arc |
| Enhanced Z_eff (PK-free core exclusion) | 1/r² vs 1/r shape mismatch requires extreme A_rep (Z_eff(0)=−38); penetration effect (Z_eff>1) overwhelms repulsive dip; any Z_eff modification destroys homonuclear symmetry (31% odd-l channels vs 0% for PK) | Decoupled architecture (constant Z_A + PK) is essential; core exclusion cannot be encoded in the nuclear potential | v2.0.5 diagnostic arc, v2.0.6 |
| Algebraic PK projector (rank-1 |core⟩⟨core|) | Rank-1 too weak; valence avoids core in function space while penetrating in coordinate space. Drift 5.6× worse than Gaussian PK (+1.697 vs +0.303 bohr/l_max) | Gaussian PK provides coordinate-space exclusion that rank-1 projector cannot. The l_max divergence root cause is the adiabatic approximation, not PK form | v2.0.6 algebraic PK diagnostic |
| Single-channel DBOC correction | 97% cancelled by off-diagonal P-coupling; overcorrects by 23×. Requires full coupled-channel solver | Use coupled-channel hyperradial solver instead of single-channel + DBOC | v2.0.6 DBOC diagnostic |

---

## 4. The Dimensionless Vacuum Principle

The graph Laplacian is dimensionless and scale-invariant, topologically equivalent to the unit three-sphere S3. Physical energies emerge only through the energy-shell constraint p0^2 = -2E, which acts as a stereographic focal length. The 1/r Coulomb potential is not an input force law -- it is the coordinate distortion from stereographic projection (chordal distance identity). The universal kinetic scale kappa = -1/16 maps graph eigenvalues to the Rydberg spectrum. Eigenvalues of the Laplace-Beltrami operator on unit S3 are pure integers: lambda_n = -(n^2 - 1).

For molecules, the natural geometry shifts from S3 (atoms) to prolate spheroidal coordinates (Paper 11), hyperspherical coordinates (Paper 13), molecule-frame hyperspherical (Paper 15), or composed fiber bundles (Paper 17). The choice of geometry is determined by where separation of variables occurs. At Level 1, the separated coordinates are all angular (on S3), making the problem fully discrete. At Levels 2+, separation produces continuous radial-like coordinates alongside discrete angular channels — a consequence of SO(4) symmetry breaking by multi-center or multi-electron potentials (Papers 8-9, negative theorem).

**Prime directive (level-aware):** The framework preserves two distinct kinds of structure, and the rules differ by level:

*Angular/symmetry structure (discrete at all levels):* The quantum number labeling (n, l, m), the per-shell degeneracy 2l+1, the selection rules from Gaunt integrals, and the channel structure are combinatorial invariants derived from the packing construction (Paper 0). These must never be modified, approximated, or bypassed at any level. They are the framework's core identity.

*Radial/parametric structure (level-dependent):* At Level 1, Fock's SO(4) symmetry converts the radial coordinate into an angular coordinate on S3, making the entire problem discrete. The graph Laplacian is the exact object and must not be modified to artificially recover continuous differential terms (like 1/r or nabla^2). At Levels 2-4, the SO(4) symmetry is broken by multi-center or multi-electron physics, and the radial-like coordinates (internuclear distance R, hyperradius R_e, prolate spheroidal xi) currently require continuous numerical methods. However, this is not necessarily a fundamental feature of the physics — it may reflect algebraic structures that haven't been found yet. The project has already demonstrated this replacement in several cases: Paper 12's Neumann expansion replaced numerical quadrature of V_ee with algebraic recurrence relations; Gaunt integrals replace angular integration at all levels; the split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. In each case, what appeared to require continuous integration was replaced by algebraic evaluation once the right structure was identified.

The places where continuous numerical methods currently survive are: Z_eff screening quadrature in composed geometries, spline caching for the rho-collapse, finite-difference radial grids in hyperspherical solvers, and the adiabatic potential curves U(R). These are computational tools for evaluating coupling matrix elements between discrete channel structures. Standard numerical methods (grid refinement, basis convergence, quadrature improvement) are appropriate for improving accuracy in these areas. Finding algebraic replacements for these numerical steps — as Paper 12 did for V_ee — is an ongoing goal of the project.

The framework's long-term aspiration is: for every natural geometry, there exists an algebraic structure that computes all coupling matrix elements without spatial integration. This has been achieved for Level 1 (S3 graph eigenvalues), for angular couplings at all levels (Gaunt integrals), and for electron-electron repulsion in prolate spheroidal coordinates (Neumann expansion). It has not yet been achieved for the hyperradial coupling in Levels 3-4 or for Z_eff screening in Level 5. What must be preserved at every level is the discrete channel structure and selection rules, not the specific numerical method used to evaluate radial amplitudes within those channels.

*Practical test:* If a proposed modification changes which quantum numbers label the states or which transitions are allowed, it violates the prime directive. If it changes how accurately the radial amplitude is computed within a given channel, it is a legitimate numerical improvement.

**Topological integrity tests:** The 18 symbolic proofs in `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py` validate S3 conformal geometry. These must never be broken or bypassed. Run before any release.

---

## 5. Natural Geometry Hierarchy

The core organizational principle of the project. Each electron configuration has a natural coordinate system where the physics separates.

| Level | System | Natural Geometry | Best Result | Paper |
|:-----:|:-------|:-----------------|:------------|:-----:|
| 1 | H (1-center, 1e) | S3 (Fock) | < 0.1% | 7 |
| 2 | H2+ (2-center, 1e) | Prolate spheroid | 0.0002% (spectral) | 11 |
| 3 | He (1-center, 2e) | Hyperspherical | 0.05% | 13 |

*Level 3 note:* The 0.05% result uses the original FD adiabatic solver at l_max=0. The algebraic spectral basis (v2.0.6) gives 0.16% single-channel at l_max=0 (removing lucky FD error cancellation). However, the algebraic solver enables coupled-channel integration, which shows convergent behavior: 0.37% → 0.27% at l_max 1-3, compared to divergent single-channel (0.56% → 0.65%). The coupled-channel solver with exact Q-matrix converges algebraically (~l_max⁻²) to a structural floor of 0.19-0.20% (v2.0.8, l_max=0-5 study). Sub-0.1% requires the variational 2D solver (Paper 15) to bypass the adiabatic approximation.
| 4 | H2 (2-center, 2e) | Mol-frame hyperspherical | 94.1% D_e | 15 |
| 5 | LiH (core+valence) | Composed (Level 3 + 4) | R_eq 5.3% | 17 |
| 5 | BeH₂ (polyatomic) | Composed (Level 3 + 4) + exchange | R_eq 11.7% | 17 |
| 5 | H₂O (triatomic) | Composed (Level 3 + 4) + lone pairs | R_eq 26% | 17 |

The composed geometry (Level 5) is a fiber bundle: G_total = G_nuc semi-direct G_core(R) semi-direct G_val(R, core_state). Each electron group gets its own natural coordinate system, coupled via Z_eff screening and Phillips-Kleinman pseudopotential.

**Algebraic structure:** At every level, angular matrix elements are computed from quantum number labels and Wigner 3j symbols (via Gaunt integrals), with no spatial quadrature. The split-region Legendre expansion (Paper 15) terminates exactly via the 3j triangle inequality. At Level 2, the radial solver is fully algebraic for all m: σ states (m=0) use ordinary Laguerre three-term recurrence with zero numerical integration (v2.0.9); π/δ states (m≠0) use associated Laguerre basis L_n^{|m|}(x) with partial-fraction decomposition and Stieltjes integral recurrence, reducing non-algebraic content to a single transcendental seed e^a·E₁(a) (Track J, v2.0.10). At Level 3, the angular problem is fully algebraic, and the hyperradial overlap S and kinetic K matrices are algebraic via three-term Laguerre recurrence (Track H, v2.0.10: pentadiagonal M2 for S, tridiagonal derivative expansion for K, < 1e-14 relative error, 11× build speedup). The adiabatic eigenvalues μ(R) are proven transcendental (O(R) → O(R²) regime transition, v2.0.9 Track G) — the potential V_eff(R) must stay quadrature, point-by-point diagonalization is irreducible, though spectral radial solvers achieve 95-120× speedups. At Level 4, the spectral Laguerre basis achieves 16× dimension reduction and < 0.0005 Ha FD agreement for the hyperradial coordinate; wall time is dominated by the angular sweep (99% of total cost), not the radial solve. n_basis=20 optimal (mild conditioning at n_basis≥25). Spatial quadrature enters for: Level 3 hyperradial potential matrix elements (V_eff transcendental), Level 4 angular eigenvalue sweeps (μ(ρ) transcendental), Z_eff screening in composed geometries, and rho-collapse spline caching.

---

## 6. Paper Series

### Reading Guide

1. **Start here:** Paper 7 (the theoretical foundation -- graph Laplacian = S3 = Schrodinger)
2. **Atoms:** Papers 0, 1 (graph construction, eigenvalue methods), then FCI-A (multi-electron)
3. **Multi-electron atoms:** Paper 13 (hyperspherical lattice, He at 0.05%, fiber bundle)
4. **Dynamics:** Paper 6 (time evolution, spectroscopy, AIMD on graph Hamiltonians)
5. **Molecules -- the problem:** Papers 8-9 (bond sphere geometry, why single-S3 fails)
6. **Molecules -- the solution:** Paper 11 (prolate spheroidal lattice, H2+ zero free params)
7. **Molecules -- two-electron:** Paper 12 (algebraic V_ee, Neumann expansion, H2 92.4% D_e)
8. **Molecules -- natural geometry:** Paper 15 (Level 4 hyperspherical, H2 94.1% D_e)
9. **Molecules -- core-valence:** Paper 17 (composed geometry, LiH 6.4%, BeH+ bound)
10. **Periodicity:** Paper 16 (periodic table from S_N representation theory on S^(3N-1))
11. **Exchange constants:** Paper 18 (why transcendentals appear at each level, Weyl's law connection)
12. **Ab initio spectroscopy:** Paper 13 Sec IX (PES -> Morse -> nuclear lattice)
13. **Quantum computing:** Paper 14 (qubit encoding, O(Q^3.15) Pauli scaling)

### Paper Inventory

#### Core (`papers/core/`)

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Paper 0 | `Paper_0_Geometric_Packing.tex` | Active | Universal constant K = -1/16 |
| Paper 1 | `paper_1_spectrum.tex` | Active | Spectral graph theory, O(N) eigenvalue methods. Berry phase retracted v1.2.0 |
| Paper 6 | `Paper_6_Quantum_Dynamics.tex` | Active | O(V) dynamics: Rabi, spectroscopy, AIMD |
| Paper 7 | `Paper_7_Dimensionless_Vacuum.tex` | Active | S3 proof (18/18 symbolic), Schrodinger recovery, SO(3N) generalization |
| Paper 10 | `paper_10_nuclear_lattice.tex` | Draft | Rovibrational spectra from SU(2) algebraic chains |
| Paper 11 | `paper_11_prolate_spheroidal.tex` | Draft | Prolate spheroidal lattice: H2+ 0.70% E_min |
| Paper 12 | `paper_12_algebraic_vee.tex` | Active | Neumann V_ee: H2 92.4% D_e, cusp diagnosis (7.6% gap) |
| Paper 13 | `paper_13_hyperspherical.tex` | Active | Hyperspherical lattice: He 0.05%, fiber bundle, ab initio spectroscopy |
| Paper 14 | `paper_14_qubit_encoding.tex` | Active | Structurally sparse qubits: O(Q^3.15) atoms, O(Q^2.5) composed; Trenev et al. Gaussian baselines |
| Paper 15 | `paper_15_level4_geometry.tex` | Active | Level 4: H2 94.1% D_e, HeH+ 93.1% D_e |
| Paper 17 | `paper_17_composed_geometries.tex` | Active | Composed geometry: LiH R_eq 6.4%, BeH+ bound, ab initio PK |
| Papers 8-9 | `Paper_8_Bond_Sphere_Sturmian.tex` | Draft | Bond sphere (positive), Sturmian structural theorem (negative), SO(4) selection rules |
| FCI-A | `paper_fci_atoms.tex` | Draft | He 0.35%, Li 1.10%, Be 0.90% |
| FCI-M | `paper_fci_molecules.tex` | Scaffold | LCAO LiH results and diagnostic arc |

#### Observations (`papers/observations/`)

| Paper | File | Status | Key Result |
|:------|:-----|:------:|:-----------|
| Paper 16 | `paper_16_periodicity.tex` | Active | Periodic table from S_N representation theory on S^(3N-1) |
| Paper 18 | `paper_18_exchange_constants.tex` | Draft | Spectral-geometric exchange constants: Weyl-Selberg identification of κ, e^a E₁(a), μ(R); α connection |

#### Conjectures (`papers/conjectures/`)

| Paper | File | Key Topic |
|:------|:-----|:----------|
| Paper 2 | `paper_2_alpha.tex` | Fine structure constant from Hopf bundle (8.8x10^-8, p = 5.2x10^-9, circulant Hermiticity, second selection, universal B_formal/N = d identity, Hopf generalization negative result, zero params) |
| Paper 3 | `paper_3_holography.tex` | Holographic entropy, spectral dimension, central charge |
| Paper 4 | `Paper_4_Universality.tex` | Mass-independence, universality, muonic hydrogen |
| Paper 5 | `Paper_5_Geometric_Vacuum.tex` | Comprehensive geometric vacuum framework (synthesis) |

---

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
| Composed qubit | `geovac/composed_qubit.py` | `build_composed_lih()`, `build_composed_beh2()`, `build_composed_h2o()` |
| Algebraic angular solver | `geovac/algebraic_angular.py` | `AlgebraicAngularSolver(Z, l_max)` |
| Algebraic coupled-channel | `geovac/algebraic_coupled_channel.py` | `solve_hyperspherical_algebraic_coupled()` |
| Hyperradial solver | `geovac/hyperspherical_radial.py` | `solve_radial()`, `solve_radial_spectral()` |
| Perturbation series (Level 3) | `geovac/algebraic_angular.py` | `perturbation_series_mu()`, `pade_approximant()` |
| Physical constants | `geovac/constants.py` | `HBAR`, `C`, `ALPHA`, etc. |

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

## 8. Coding Standards

### Sparse vs Dense: Context-Dependent

- **Hamiltonian and CI matrices (N > 100):** Always `scipy.sparse` (csr_matrix, coo_matrix). Never densify.
- **Hot-loop lookup tables (ERI, h1 in direct CI):** Use dense NumPy when array fits in memory (n_spinorb <= ~300). `scipy.sparse._validate_indices` overhead (~24us/call) is prohibitive at 100K+ lookups.

Rule of thumb: sparse for the physics matrix (N_SD x N_SD), dense for orbital-index lookup tables (n_spinorb x n_spinorb or n_spatial^4).

### Type Hints Required

All function signatures must have type hints.

```python
def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    ...
```

### Physical Constants

Import from `geovac.constants` or define at module top. No hardcoded magic numbers.
**Exception:** `-1/16` is the universal topological constant (can be used directly).

### Vectorization Over Loops

Avoid Python loops for graph operations; use NumPy masking/vectorization.

```python
mask = (n_values >= 1) & (l_values < n_values)
states_filtered = states[mask]
```

---

## 9. Workflow Protocols

### Theory Check Rule

Before implementing new physics:
1. Check `papers/` for the derivation
2. If code contradicts paper -> flag it and ask user
3. If changing physics in code -> prompt user to update papers

### Benchmarking Rule

After any modification to `hamiltonian.py`, `lattice.py`, or `solver.py`:
1. Run topological integrity: `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py -v`
2. Run validation: `pytest tests/advanced_benchmarks.py`
3. Verify 18/18 symbolic proofs pass (topological foundation)
4. Verify H2+ < 0.1% error (topological control)
5. Verify H2 Full CI < 1.0% error (accuracy control)
6. Report any speed regression > 10%

### Clean Room Rule

- Generated plots -> `debug/plots/` or `papers/figures/`
- Generated data -> `debug/data/`
- Scripts -> `debug/`, `demo/`, `tests/`, or `benchmarks/` (never root)
- Documentation -> `docs/`

### Changelog Protocol

**Changelog granularity:** Each CHANGELOG entry should correspond to a releasable unit of work — a new feature, a completed diagnostic arc, a paper update, or a benchmark result. Do not create a new version entry for intermediate debugging steps, parameter sweeps, or exploratory runs. If a task requires multiple iterations to complete, document it as a single entry when the task concludes, noting the key findings and negative results. Intermediate data belongs in `debug/` with timestamped filenames, not in the CHANGELOG.

**Version numbering:** Patch versions (x.y.Z) for bug fixes and documentation; minor versions (x.Y.0) for new features, completed diagnostic arcs, or paper additions; major versions (X.0.0) for architectural changes. A diagnostic arc that tests 10 hypotheses and finds 9 negative results is one minor version, not 10 patches.

---

## 10. Validation Benchmarks

| Test | Max Error | Purpose |
|:-----|:---------:|:--------|
| Symbolic proofs (18 tests) | 0 failures | Topological foundation |
| H (hydrogen) | < 0.1% | Basic validation |
| He+ (helium ion) | < 0.1% | Z-scaling check |
| H2+ (ionized H2, FD) | < 0.1% | Topological control |
| H2+ (prolate spheroidal, spectral) | < 0.001% | Spectral Laguerre accuracy control |
| He (hyperspherical) | < 0.1% | Multi-electron control |
| H2 Full CI | < 1.0% | Accuracy control |
| H2 Neumann V_ee | 92.4% D_e | Algebraic integral accuracy |
| H2 Level 4 (2D solver) | 94.1% D_e | Molecule-frame hyperspherical |
| HeH+ Level 4 | 93.1% D_e | Heteronuclear extension |
| LiH Composed (ab initio PK) | R_eq 6.4% | Composed geometry |
| BeH+ Composed | Bound, physical | Transferability |
| Hyperspherical (20 tests) | 0 failures | Angular + adiabatic + radial |
| Muonic H energy ratio | < 0.01% | Mass-independence |
| V_ee S3 overlap (1s-1s, 1s-2s, 2s-2s) | < 0.01% | Topological integrity |
| Direct CI vs matrix CI | < 1e-8 Ha | Algorithmic consistency |
| JW eigenvalue consistency | < 1e-4 Ha | Qubit encoding correctness |
| H2 STO-3G Pauli terms | exactly 15 | Gaussian reference validation |
| He cc-pVDZ Pauli terms | exactly 156 | Computed Gaussian baseline |
| He cc-pVTZ Pauli terms | exactly 21,607 | Computed Gaussian baseline |
| He cc-pVDZ FCI energy | < 0.001 Ha vs published | Integral engine validation |
| He cc-pVTZ FCI energy | < 0.0002 Ha vs published | Integral engine validation |
| QWC grouping correctness | 0 violations | Measurement group integrity |
| LiH composed Pauli terms (Q=30) | exactly 334 | Composed qubit validation |
| BeH2 composed Pauli terms (Q=50) | exactly 556 | Composed qubit validation |
| H2O composed Pauli terms (Q=70) | exactly 778 | Composed qubit validation |
| Composed cross-block ERIs | exactly 0 | Block-diagonal integrity |
| Algebraic angular Casimir (R=0) | < 10⁻⁶ | SO(6) eigenvalue correctness |
| Algebraic a₁ perturbation | < 10⁻⁸ | First-order perturbation theory |
| Algebraic GL vs quad consistency | < 10⁻¹⁰ | Quadrature correctness |
| Algebraic l_max monotonic convergence | monotonic | Centrifugal singularity elimination |
| Spectral vs FD consistency | < 0.01 Ha | Spectral and FD agree within FD error |
| Spectral dimension reduction | > 100× | 20 basis functions vs 5000 FD grid |
| Spectral speedup | > 10× | Wall time reduction control |
| l_max=5 coupled-channel floor | 0.15%–0.25% | Adiabatic ceiling characterization |
| n_channels=5 vs 3 consistency | < 1 mHa | Channel truncation validation |
| Spectral PES (H₂⁺) | R_eq < 0.5%, E_min < 0.001% | Spectral PES scan accuracy |
| Algebraic vs quadrature matrices | < 1e-12 elementwise | Laguerre recurrence correctness |
| Algebraic PES energy | < 1e-14 Ha vs quadrature | End-to-end algebraic validation |
| Level 3 spectral-FD consistency | < 0.0001 Ha | Hyperradial spectral solver |
| Level 3 spectral coupled-channel | 0.15%–0.25% at l_max=3 | Coupled-channel ceiling preserved |
| Level 3 spectral dimension reduction | > 100× | 25 basis functions vs 3000 FD grid |
| Level 3 spectral coupled speedup | > 50× | Coupled-channel wall time |
| Perturbation a₁ vs Paper 13 | < 10⁻⁹ | RS series first-order validation |
| Perturbation series convergence | converges R < 2 bohr | Convergence radius characterization |
| 2D solver composed LiH | drift < +0.15 bohr/l_max | 2D vs adiabatic drift comparison |
| Level 3 algebraic vs quadrature S, K | < 1e-10 elementwise | Laguerre recurrence correctness (Track H) |
| Level 3 algebraic energy consistency | < 1e-10 Ha | Algebraic S+K gives identical physics |
| Level 3 algebraic ceiling unchanged | 0.15%–0.25% at l_max=3 | Coupled-channel ceiling preserved with algebraic |
| Level 4 spectral vs FD consistency | < 0.001 Ha | Level 4 spectral hyperradial solver |
| Level 4 spectral D_e% match | within 0.5% | D_e% preserved under spectral substitution |
| Level 4 spectral convergence plateau | n_basis 20-30 within 0.0001 Ha | Convergence plateau verified |
| Speed regression | < 10% | Performance control |

---

## 11. Topic-to-Paper Lookup

| Topic | Paper | Section | Tier |
|:------|:-----:|:--------|:----:|
| Universal constant -1/16 | 0 | Sec 2 | Core |
| Graph Laplacian method | 1 | Sec 3 | Core |
| O(V) quantum dynamics | 6 | All | Core |
| Rabi oscillations | 6 | -- | Core |
| Delta-kick spectroscopy | 6 | -- | Core |
| AIMD / Langevin thermostat | 6 | -- | Core |
| Dimensionless vacuum proof | 7 | All | Core |
| Schrodinger recovery | 7 | Sec 4 | Core |
| S3 conformal geometry | 7 | Sec 3 | Core |
| V_ee on S3 (node overlap) | 7 | Sec VI | Core |
| Slater F0 master formula | 7 | Sec VI.B | Core |
| SO(3N) generalization | 7 | Sec VI | Core |
| Nuclear rovibrational spectra | 10 | All | Core |
| Prolate spheroidal lattice | 11 | All | Core |
| Neumann V_ee expansion | 12 | Sec III-V | Core |
| Prolate spheroidal CI (H2) | 12 | Sec VI | Core |
| Cusp diagnosis (7.6% gap) | 12 | Sec VII | Core |
| Hyperspherical coordinates | 13 | Sec II | Core |
| Angular eigenvalue (Gaunt) | 13 | Sec III | Core |
| Adiabatic potential curves | 13 | Sec IV | Core |
| Fiber bundle structure | 13 | Sec VII | Core |
| Natural geometry hierarchy | 13 | Sec VIII | Core |
| Ab initio spectroscopy | 13 | Sec IX | Core |
| Algebraic structure (SO(6)) | 13 | Sec XII | Core |
| Algebraic angular solver | 13 | Sec XIII | Core |
| Coupled-channel convergence | 13 | Sec XIII.C | Core |
| Hellmann-Feynman P-matrix | 13 | Sec XIII.A | Core |
| Adiabatic bottleneck (l_max) | 13 | Sec XIII.B | Core |
| Perturbation series μ(R) | 13 | Sec XII-XIII | Core |
| Transcendental boundary (μ(R)) | 13 | Sec XII.B | Core |
| Spectral Laguerre (Level 2) | 11 | Sec V | Core |
| Algebraic Laguerre moments | 11 | Sec V | Core |
| Spectral Laguerre (Level 3) | 13 | — | Core |
| Algebraic Laguerre S, K (Level 3) | 13 | — | Core |
| Spectral Laguerre (Level 4) | 15 | Sec VI.H | Core |
| 2D solver in composition | 15, 17 | Sec VI.D, — | Core |
| Qubit Pauli scaling | 14 | All | Core |
| Structural sparsity | 14 | Sec III | Core |
| QWC measurement groups | 14 | Sec III.E | Core |
| Pauli 1-norm scaling | 14 | Sec III.F | Core |
| Trotter error bounds | 14 | Sec III.F | Core |
| Genuine Gaussian baselines | 14 | Sec III.G | Core |
| VQE head-to-head benchmark | 14 | Sec III.G | Core |
| Equal-qubit Gaussian comparison | 14 | Sec III.G | Core |
| Composed qubit Hamiltonians | 14 | Sec IV | Core |
| Composed Pauli scaling (Q^2.5) | 14 | Sec IV.B | Core |
| Trenev et al. Gaussian baselines | 14 | Sec IV.C | Core |
| Level 4 mol-frame hyperspherical | 15 | All | Core |
| Mol-frame charge function | 15 | Sec III | Core |
| Multichannel expansion | 15 | Sec V | Core |
| Heteronuclear extension | 15 | Sec V.D | Core |
| Variational 2D solver | 15 | Sec VI.D | Core |
| HeH+ convergence | 15 | Sec VI.E | Core |
| Double-adiabatic fiber bundle | 15 | Sec VII.C | Core |
| Chemical periodicity (S_N reps) | 16 | All | Observation |
| Structure types A/B/C/D/E | 16 | Sec III | Observation |
| Hierarchical decomposition | 16 | Sec IV | Observation |
| Dirac limit (Z~137) | 16 | Sec VI | Observation |
| Exchange constants (Weyl-Selberg) | 18 | All | Observation |
| Stieltjes seed e^a E₁(a) | 18 | Sec II.B | Observation |
| Transcendence hierarchy | 18 | Sec IV, VI | Observation |
| α as exchange constant | 18 | Sec V | Observation |
| Composed natural geometries | 17 | All | Core |
| Core-valence fiber bundle | 17 | Sec II | Core |
| Ab initio Phillips-Kleinman | 17 | Sec IV | Core |
| Rho-collapse cache | 17 | Sec V | Core |
| l_max divergence root cause | 17 | Sec VI.A | Core |
| LiH/BeH+ benchmarks | 17 | Sec VI | Core |
| Bond sphere theory | 8-9 | All | Core |
| Sturmian structural theorem | 8-9 | Sec IV | Core |
| SO(4) selection rules | 8-9 | Sec III | Core |
| Fine structure alpha (Hopf bundle) | 2 | Sec 3-5 | Conjecture |
| Hopf fibration spectral invariants | 2 | Sec 3 | Conjecture |
| Cubic self-consistency (alpha) | 2 | Sec 4 | Conjecture |
| Spectral dimension d_s | 3 | Sec 4 | Conjecture |
| Holographic entropy S | 3 | Sec 5 | Conjecture |
| Central charge c | 3 | Sec 6 | Conjecture |
| Mass-independence | 4 | Sec 3-4 | Conjecture |
| Muonic hydrogen | 4 | Sec 5 | Conjecture |
| Contact geometry | 4 | Sec 5 | Conjecture |
| Comprehensive framework | 5 | All | Conjecture |

---

## 12. Algebraic Registry

Tracks which matrix elements at each level are computed algebraically vs numerically. Status: **algebraic** (closed-form from quantum numbers), **algebraic-pending** (algebraic route identified but production code still uses quadrature), **numerical-required** (no known algebraic replacement).

### Level 3 (Hyperspherical — He)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| SO(6) Casimir eigenvalues | algebraic | Exact: ν(ν+4)/2 from representation theory |
| Centrifugal barrier | algebraic | Diagonal in Gegenbauer spectral basis (confirmed Track B, v2.0.6) |
| Nuclear coupling | algebraic | Partial harmonic sums (Eqs. 31-32, Paper 13) |
| V_ee coupling | algebraic-pending | Split-region Legendre structure confirmed; GL quadrature still used in production |
| Centrifugal matrix elements | algebraic | Diagonal in Gegenbauer basis (v2.0.6) |
| Hyperradial overlap (S) | algebraic | Pentadiagonal M2 moment matrix from three-term Laguerre recurrence. S = M2/(8α³). Machine-precision agreement with quadrature (< 1e-14 relative). (Track H, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Derivative kernel B_n = -n/2 L_{n-1} + 1/2 L_n + (n+1)/2 L_{n+1} (tridiagonal). K = bbᵀ/(4α). Machine-precision agreement (< 1e-14 relative). 11× build speedup. (Track H, v2.0.10) |
| Hyperradial potential (V_eff) | numerical-required | V_eff(R) = μ(R)/R² + 15/(8R²) is transcendental (Track G proven). Quadrature required. Analogous to Level 2 m≠0 centrifugal singularity. |
| Hellmann-Feynman P-matrix | algebraic | Exact from R-independent dH/dR (v2.0.6) |
| Q-matrix (second derivative coupling) | algebraic | Exact Q = PP + dP/dR computed from Hellmann-Feynman quantities (v2.0.6) |
| Coupled-channel radial solve | numerical-required | Coupled ODE integration on R grid |
| Adiabatic potential curves U(R) | numerical-required | Small-matrix diagonalization at each R (10-30 dim, algebraic solver) or FD grid (200-800 dim, original) |

### Level 2 (Prolate Spheroidal — H₂⁺)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Angular η-eigenfunctions | algebraic | Legendre spectral basis (n_basis=50), same as Mitnik et al. angular component |
| Azimuthal m | algebraic | Exact separation of variables, integer quantum number |
| Radial ξ-solver (m=0) | algebraic | Ordinary Laguerre basis (n_basis=20, α-adapted). 250× dimension reduction vs FD. Three-term recurrence for ALL matrix elements — zero quadrature. Machine-precision agreement with quadrature (< 1e-14). (v2.0.9) |
| Radial ξ-solver (m≠0) | algebraic | Associated Laguerre basis L_n^{|m|}(x) with weight x^|m|·e^{-x}. Partial-fraction decomposition of centrifugal 1/x singularity into lowered moment M_{-1} (algebraic, DLMF 18.9.13) + Stieltjes integral J (three-term recurrence). Single transcendental seed e^a·E₁(a); all other elements algebraic. Associated basis converges faster than ordinary for m=1 (stable by N=10). (Track J, v2.0.10) |
| Separation parameter c² | numerical-required | Iterative root-finding (Brent method) for self-consistency between angular and radial equations |
| Coupled-channel ceiling | characterized | Error floor 0.19-0.20% from adiabatic approximation; algebraic convergence ~l_max⁻²; 3 channels sufficient; n_basis and R-grid converged. Sub-0.1% requires non-adiabatic (2D variational) solver. (v2.0.8) |

### Level 4 (Mol-Frame Hyperspherical — H₂)

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Hyperradial overlap (S) | algebraic | Same Laguerre three-term recurrence as Level 3. Pentadiagonal M2 moment matrix. (Track I, v2.0.10) |
| Hyperradial kinetic (K) | algebraic | Same derivative kernel as Level 3. Tridiagonal expansion. (Track I, v2.0.10) |
| Hyperradial potential (U_eff) | numerical-required | Adiabatic curve from angular sweep. μ(ρ) transcendental (same root cause as Level 3, Track G). |
| Angular eigenvalue sweep | numerical-required | Point-by-point diagonalization of H_ang at ~130 ρ-values. Dominates 99% of Level 4 wall time. |
| 2D tensor product assembly | numerical-required | H_ang evaluated at each quadrature point for (R_e, α) tensor product. Dense kronecker assembly. |
| Nuclear coupling (split-region Legendre) | algebraic | Gaunt integrals, exact via 3j triangle inequality (Paper 15). |

---

## 13. Multi-Agent Protocol

The GeoVac project uses an AI-augmented agentic workflow with a formalized planner–worker architecture for Claude Code sessions.

### 13.1 Architecture

Three layers:

1. **Plan mode (human + chat):** Strategic direction, framing decisions, result synthesis. Conceptual changes to papers and CLAUDE.md originate here. This layer is not automated.
2. **PM agent (main Claude Code session, full context):** Reads CLAUDE.md, README, and all papers relevant to the current track. Decomposes tasks into sub-agent work units. Evaluates sub-agent results against verification checklists. Does NOT write production code itself — it plans, delegates, and synthesizes.
3. **Worker sub-agents (scoped context):** Execute specific coding, computation, or drafting tasks. Each receives only the files it needs. Reports results back to the PM agent.

### 13.2 PM Session Kickoff

Every PM session begins by reading CLAUDE.md and then executing the following:

1. Identify the current track(s) and relevant papers from the plan mode directive
2. Read those papers and any results from the previous session
3. Check the failed approaches table (Section 3) for relevant dead ends
4. Decompose the session goal into sub-agent tasks
5. Draft sub-agent prompts using the standard template (13.3)
6. Identify which papers need updating based on the session's results (see 13.8)

### 13.3 Sub-Agent Prompt Template

Every sub-agent dispatch uses this format:

```
CONTEXT FILES: [minimal set of files to read — list explicitly]
TASK: [one clear deliverable, stated in one sentence]
CONSTRAINTS:
  - Failed approaches to avoid: [list relevant entries from Section 3]
  - Structures that must be preserved: [quantum numbers, selection rules, etc.]
  - Do NOT modify: [list any files that are off-limits]
SUCCESS CRITERIA:
  - Tests to pass: [specific test files or assertions]
  - Numerical targets: [specific error thresholds if applicable]
  - Consistency check: [which papers or results to verify against]
OUTPUT FORMAT: [what to report back — specific data, not just pass/fail]
PAPER UPDATES:
  - Papers affected: [list paper numbers]
  - Autonomous changes: [list what can be committed directly per 13.8]
  - Escalation changes: [list what needs plan-mode review, with proposed diffs]
```

### 13.4 Verification Gates

Before the PM agent accepts a sub-agent result, it checks:

1. **Test gate:** Do all relevant tests pass? (Non-negotiable.)
2. **Dead-end gate:** Does the approach match any entry in the failed approaches table? If so, reject unless the sub-agent explicitly explains what is different this time.
3. **Prime directive gate:** Does the result modify any discrete structure — quantum number labeling, selection rules, channel structure, Gaunt integral coupling? If so, do NOT accept. Escalate to plan mode for human review.
4. **Consistency gate:** Does the result contradict any established result in the papers? If uncertain, flag for PM review rather than accepting.

### 13.5 Escalation Rules

The following changes require escalation to plan mode (human review) and must NOT be made by sub-agents or the PM agent autonomously:

- Any modification to CLAUDE.md
- Any change to the natural geometry hierarchy (new levels, changed coordinate systems)
- Any new entry in the failed approaches table
- Any result that would change the "Best Result" column in the hierarchy table
- Introduction of any fitted or empirical parameter
- Paper changes in the **escalation tier** (see 13.8)

### 13.8 Paper Update Policy

Papers are the authoritative source for all physics (Section 1). Code that outpaces the papers creates documentation drift. PMs are expected to keep papers in sync with code results, subject to a tiered autonomy rule.

#### Autonomous tier (PMs may commit directly)

PMs may modify papers in `papers/core/`, `papers/observations/`, or `papers/conjectures/` WITHOUT plan-mode approval for changes that are strictly additive and factual:

- Adding or updating benchmark tables (new numerical results, updated error percentages)
- Adding rows to existing comparison tables
- Updating numbers that have improved (e.g., error percentages, speedup factors, dimension counts)
- Adding a subsection documenting a new solver variant, algorithm, or computational method with its results — provided it implements physics already described in the paper or in CLAUDE.md
- Adding references to new code modules or test files
- Fixing typos, grammatical errors, or LaTeX formatting issues

**The test:** Does the change add evidence or detail for claims the paper already makes? If yes, autonomous. If it changes what the paper claims, escalate.

#### Escalation tier (must go through plan-mode review)

The following paper changes must NOT be made autonomously. The PM should draft the proposed changes and include them in the lane report under a **"Paper updates for review"** section with full diffs. These are reviewed in plan-mode chat at sprint-end.

- Changing, softening, or strengthening any claim or conclusion
- Reframing a result (e.g., from "demonstrated" to "proven", or vice versa)
- Adding or modifying theoretical arguments or derivations
- Changing the abstract or introduction framing
- Adding a new section that introduces a new concept, principle, or interpretation
- Documenting a negative result that changes the paper's narrative arc
- Modifying or retracting any previously published result
- Any change the PM is uncertain about

#### PM report format for escalation-tier changes

When a PM has escalation-tier paper changes, the lane report should include:

```
## Paper updates for review

### Paper N — [title]

**Type:** [new subsection / claim modification / narrative change / etc.]
**Rationale:** [why this change is needed, in 1-2 sentences]

**Proposed diff:**
[exact LaTeX to add/change, with enough surrounding context to locate it]
```

These are reviewed in the plan-mode sprint-end chat. The human approves, modifies, or rejects each proposed change.

### 13.6 Track Management

Active work is organized into tracks. Each track has:

- A name and one-sentence goal
- A list of relevant papers and code modules
- A current status (active / blocked / complete)
- A log of sub-agent dispatches and results (maintained by the PM in `debug/track_logs/`)

The PM agent maintains a brief track status file at `debug/track_logs/STATUS.md` that is updated at the end of each session. This file is read at the start of the next PM session to restore context.

### 13.7 General Guidance

This protocol is a starting point and will be refined based on experience. The overriding principles are: **code changes can be autonomous** (the test suite catches errors); **factual paper updates are autonomous** (adding results, tables, and new-method subsections keeps papers in sync with code); **but conceptual changes cannot be autonomous** (claims, framing, and theoretical arguments require human judgment in plan mode). When in doubt, escalate rather than proceeding.
