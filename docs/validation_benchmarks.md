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
| H2 Level 4 (2D solver) | 96.0% D_e | Molecule-frame hyperspherical |
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
| H2 bond-pair Pauli terms (Q=10) | exactly 112 | Bond-pair qubit validation |
| H2 bond-pair R-independence | 112 at all R | Selection rule sparsity |
| Composed cross-block ERIs | exactly 0 | Block-diagonal integrity |
| HF composed Pauli terms (Q=60) | exactly 667 | New molecule validation |
| NH₃ composed Pauli terms (Q=80) | exactly 889 | New molecule validation |
| CH₄ composed Pauli terms (Q=90) | exactly 1000 | New molecule validation |
| General builder exact reproduction | < 1e-14 | Refactored builder matches hardcoded |
| Atomic classifier Z=1-10 | 101 tests pass | First-row classification |
| TC radial-only LiH Pauli (Q=30) | exactly 562 | TC composed radial validation (BX-3) |
| TC angular Gaunt selection rules | |ΔL|=1, |Δm|≤1 | Angular gradient preserves Gaunt rules |
| TC angular He max_n=1 identity | rad = full | No l>0 orbitals → identical |
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
| Level 4 spectral angular vs FD U_min | < 2e-5 Ha | Angular spectral solver accuracy (Track K) |
| Level 4 spectral angular speedup | > 200× | Angular sweep speedup control (Track K) |
| Level 4 spectral angular dimension | 20× reduction | 1000 → 50 matrix dimension (Track K) |
| Cusp factor baseline reproduction | < 1e-14 eigenvalue diff | f=1 matches standard solver (Track U) |
| Cusp factor D_e degradation | D_e decreases with γ | Negative result verified (Track U) |
| S⁵ Green singularity order | 1/d³ (not 1/d¹) | Dimensionality mismatch proof (Track W) |
| He cusp correction sign | ΔE < 0 | Cusp lowers energy (Track X) |
| He cusp correction l_max=2 | error < 0.15% | Breaks through 0.19-0.20% floor (Track X) |
| H₂ cusp correction R-dependent | varies with R | Required for D_e correction (Track X) |
| Cusp correction convergence | → 0 as l_max → ∞ | Basis-dependent, not systematic error (Track X) |
| 4-electron LiH equilibrium at l_max≥2 | structural | PK-free equilibrium validation (Track AJ) |
| S₄ [2,2] channel reduction | verified | N-electron antisymmetry machinery (Track AJ) |
| N-electron spectral compression | ≥ 100× | FD-to-spectral ratio (Track AK) |
| N-electron 2D variational bound | E_min > exact | Variational principle respected (Track AR) |
| N-electron 2D D_e sign | D_e < 0 (unbound) | Adiabatic overcounting confirmed as artifact (Track AR) |
| Full vs composed Pauli terms | full > composed | Full encoding categorically denser (Track AS) |
| Balanced coupled LiH Pauli terms (Q=30) | exactly 878 | Balanced coupled qubit validation |
| Balanced coupled LiH 4e FCI bound | D_e > 0 | Only bound 4e config |
| Balanced coupled LiH n_max=3 energy | -8.055 Ha (0.20% err) | Convergence validation |
| Balanced coupled variational bound | E > exact at all R | Variational principle |
| Balanced coupled BeH₂ Pauli terms (Q=50) | exactly 2,652 | Polyatomic balanced validation |
| Balanced coupled BeH₂ 1-norm < composed | 304.7 < 354.9 Ha | 1-norm advantage verification |
| Balanced coupled BeH₂ 4e FCI | 10.7% error | Polyatomic energy validation |
| Balanced coupled H₂O Pauli terms (Q=70) | exactly 5,798 | Non-collinear balanced validation |
| Balanced coupled H₂O 1-norm | 1,509.3 Ha | Non-collinear 1-norm |
| Balanced coupled non-collinear V_ne | direction=(0,0,-1) matches nuc_parity to 1e-14 | Wigner D rotation validation |
| FrozenCore Z_eff asymptotic | Z_eff(0)≈Z, Z_eff(∞)≈Z-10 | Ne-like screening validation |
| FrozenCore density normalization | integral = 10 ± 1% | Core electron count |
| Second-row Pauli scaling | Q^2.50 | O(Q^2.5) universality |
| NaH balanced Pauli (Q=20) | exactly 239 | Second-row qubit validation |
| Atomic classifier Z=11-18 | 97 tests pass | Second-row classification |
| Wigner D l=2 orthogonality | R^T R = I to 1e-12 | l=2 rotation validation |
| Wigner D l=2 determinant | det(R) = +1 | Proper rotation check |
| Wigner D l=1 legacy match | < 1e-13 | General vs legacy consistency |
| Block rotation l=0,1,2 | R^T R = I to 1e-12 | Mixed-l block validation |
| NaH n_max=3 build | succeeds | l=2 rotation unblock |
| NaH n_max=3 FCI bound | E(3) < E(2) at all R | Variational convergence |
| LiF composed Pauli terms (Q=70) | exactly 778 | Multi-center qubit validation |
| CO composed Pauli terms (Q=100) | exactly 1111 | Multi-center qubit validation |
| N₂ composed Pauli terms (Q=100) | exactly 1111 | Multi-center isostructural invariance |
| CO = N₂ Pauli count | identical | Isostructural invariance (multi-center) |
| F₂ composed Pauli terms (Q=100) | exactly 1111 | Multi-center qubit validation |
| NaCl composed Pauli terms (Q=50) | exactly 556 | Mixed frozen-core multi-center |
| CH₂O composed Pauli terms (Q=120) | exactly 1333 | Multi-center polyatomic validation |
| C₂H₂ composed Pauli terms (Q=120) | exactly 1333 | Multi-center polyatomic validation |
| C₂H₆ composed Pauli terms (Q=160) | exactly 1777 | Multi-center polyatomic validation |
| Composed Pauli/Q ratio | 11.11 ± 0.1 | Universal linear scaling law |
| All TM hydrides composed Pauli terms (Q=30) | exactly 277 (non-identity) | Transition metal qubit validation (isostructural: all 10 identical). Note: Track CZ/DA reported 278 including the identity term; v2.8.0 standardized on excluding identity from N_pauli across all builders. |
| SrH / BaH composed Pauli terms (Q=20) | exactly 222 (non-identity) | Heavy-atom alkaline-earth monohydride validation ([Kr], [Xe] frozen cores); isostructural with KH, NaH, CaH (Sprint 3 HA-C, v2.12.0). |
| SrH_rel / BaH_rel composed Pauli terms (Q=20) | exactly 942 (non-identity) | Relativistic heavy-atom monohydride validation; isostructural with CaH_rel (post-TR, Sprint 4 v2.15.0; pre-TR was 534). |
| SrH_rel / BaH_rel / CaH_rel λ_ni | bit-identical 13.87 Ha | Cross-species relativistic 1-norm invariance (frozen core screens Z, spin-orbit sees Z_eff=2 uniformly). |
| SrH_rel / BaH_rel / CaH_rel QWC | bit-identical 52 groups | Cross-species QWC structural invariance. |
| d-only block ERI density | 4.0% | d-orbital Gaunt sparsity |
| d-block Pauli/Q (composed) | 9.23 (< main-group 11.11) | d-orbital sparsity advantage |
| Nested Be Pauli terms (Q=10) | exactly 112 | Nested encoding qubit validation |
| Nested Be 1-norm < composed | 18.95 < 121.35 Ha | PK elimination 1-norm advantage |
| H-set ERI density < uncoupled | 9.2% < 12.5% (l_max=1) | 6j recoupling sparsity |
| LiH/BeH₂/H₂O Pauli unchanged | 334/556/778 | Backward compatibility regression |
| He 2D variational bound | E > exact at all l_max | Variational principle |
| He 2D l_max monotonic | E decreasing with l_max | Convergence validation |
| He 2D breaks adiabatic floor | < 0.10% at l_max=5 | Non-adiabatic improvement |
| He 2D cusp-corrected | < 0.01% at l_max=4 | Track DI Phase 1 target |
| He 2D radial converged | stable at n_R ≥ 20 | Radial basis validation |
| Casimir CI n_max=1 k* | 9/4 (exact) | SC Hartree screening |
| Casimir CI n_max=1 E_var | -729/256 (exact) | Variational Hartree screening |
| Casimir CI polynomial structure | residual < 1e-10 | H(k) = Bk + Ck² |
| Casimir CI variational bound | E_var > exact at all n_max | Variational principle |
| Casimir CI n_max=3 error | < 2% (variational) | Algebraic CI accuracy |
| Graph-native CI n_max=7 error | < 0.20% | Graph-native CI accuracy |
| Graph-native CI n_max=8 error | 0.207% (2,262 configs) | Exact algebraic float integrals |
| Graph-native CI n_max=9 error | 0.201% (3,927 configs) | Exact algebraic float integrals |
| He 2D variational self-consistent cusp | < 0.020% | Self-consistent cusp correction (l_max=7, n_R=35) |
| Hylleraas-Eckart He 1¹S (ω=4, 22 basis) | < 0.001% | Hylleraas-Eckart double-α, Track 1 closure |
| Hylleraas-Eckart He 2¹S-2³S splitting | -1.4% | Eckart 1933, Track 3 closure |
| Hylleraas-Eckart He 2¹P→1¹S f-value | -2.02% vs Drake 0.2761 | Full Schwartz 1961 two-channel, Track 5 closure |
| Hylleraas-Eckart hydrogen 1S→2P f | < 1% vs 0.4162 | Wigner-Eckart factor-of-2 sanity check |
| P-state quadrature kinetic vs algebraic | < 1.5e-5 (sym×sym) | Universal kinetic validation |
| P-state kinetic Hermiticity all channels | < 1e-8 | T(p,q) = T(q,p) at β=0.3 |
| P-state antisym kinetic at β=0 | 0.0 identically | Basis vanishes (sinh→0) |
| P-state cross-sector kinetic at β=0 | 0.0 identically | Cross-sector vanishes (sinh→0) |
| Algebraic Slater R^k machine precision | < 1.5e-12 | `compute_rk_float()` vs exact Fraction |
| casimir_ci F²(2p,2p) corrected | 45/512 (was 43/512) | Typo fix validated by hypergeometric evaluator |
| Graph-native CI beats diagonal | error_graph < error_diag | Graph topology dominance |
| Graph-consistent = graph-native | < 1e-10 Ha | FCI basis invariance |
| Graph-consistent full spectrum | all evals match to 1e-8 | FCI orbital rotation invariance |
| Graph-native CI Z_c crossover | Z_c ≈ 1.84 at n_max=7 | Graph validity boundary |
| Graph-native H⁻ over-binding | E_CI < E_exact at n_max≥2 | Non-variational characterization |
| Standard FCI H⁻ variational | E_CI > E_exact at all n_max | Variational bound verification |
| He energy decomposition verified | <h1>+<V_ee>=E to 1e-12 | Decomposition correctness |
| V_ee full rank (graph eigenbasis) | rank = dim at all n_max | Cusp rank characterization |
| Diagonal V_ee converged by n_max=3 | < 0.01 mHa change | Mean-field convergence |
| PsH bound (adiabatic) | V_min < -0.5 Ha | Exotic atom binding |
| PsH energy l_max=3 | 4.1% error | Exotic atom accuracy |
| 111 Pauli per s/p block | exactly 111 | Composed Pauli derivation |
| 111 = 55 + 56 decomposition | 55 direct + 56 exchange | Pauli channel decomposition |
| TC gamma_opt(l=3) ≈ 0.10 | sweet spot | TC optimization |
| Dirac-on-S³ π-free certificate (Weyl sector) | 0 non-rationals, n_max=6 | Analog of Paper 24 §III Bargmann-Segal certification; every \|λ_n\| is exact sympy Rational (n + 3/2), every g_n^Weyl is positive int ((n+1)(n+2)) |
| Dirac-on-S³ π-free certificate (Dirac sector) | 0 non-rationals, n_max=6 | Every \|λ_n\| is Rational, every g_n^Dirac is positive int (2(n+1)(n+2)) |
| Dirac-on-S³ label generator | exactly g_n labels per level | `spinor_labels_at_n` generates exactly (n+1)(n+2) Weyl or 2(n+1)(n+2) Dirac labels; stronger invariant than bare π-free |
| Dirac Δ⁻¹ identity | g_3^Dirac = 40 exactly | Phase 4H SM-D identity (Δ = 1/(g_3^Dirac)) reproduced in D1 API as `delta_inverse_identity() == (40, Rational(1,40))` |
| Fock ↔ CH convention conversion | invertible at all n | `fock_to_ch(ch_to_fock(n)) == n` for n = 1..10; label compatibility with scalar graph (n_Fock = n_CH + 1) |
| D2 cumulative Dirac trace closed form | exact Rational | Σ g_n^Dirac = (N+1)(N+2)(2N+3)/3 symbolic identity, verified for N = 0..5 |
| D2 \|λ_{m-1}\|·g_{m-1}^Weyl = 42 at m=3 | (m-1)(m+2) = 10 uniquely | Single-point Dirac/Weyl coincidence with B = 42 occurs only at m=3, sympy exact |
| D3 Dirac Dirichlet at s=4 | 2ζ(2) + 2ζ(3) | `summation(2*m*(m+1)*m**(-4), (m,1,oo)) == pi**2/3 + 2*zeta(3)` symbolic |
| D3 Weyl Dirichlet at s=4 | ζ(2) + ζ(3) | Factor of 2 difference from Dirac, sympy exact |
| D4 Hopf charge partition sum | Σ mult = g_n^Dirac | Sympy-exact at n_CH = 0..5; 40 = 20 half-integer-charge + 20 integer-charge at n_CH = 3 |
| D4 Dirac Hurwitz spectral zeta at s=4 | π² − π⁴/12 | `summation(2*(n+1)*(n+2)/(n+Rational(3,2))**4, (n,0,oo))` symbolic closed form |
| T0 d_spinor pair-diag at l_max=0..5 | exact rational | 1/4, 11/20, 553/724, 101/118, 2533/2820, 9611/10396 from sympy Wigner 3j |
| T0 d_spinor full-Gaunt at l_max=0..5 | exact rational | 1/4, ..., 0.923 (sec'dry table, physically correct Coulomb rule) |
| T0 d_spinor ≤ d_scalar | monotonic ∀ l_max | Spinor basis sparsity bounded by scalar sparsity |
| T0 scalar density reproduces Paper 22 | bit-exact | Pair-diag convention matches Paper 22 Table §III |
| `dirac_matrix_elements` module tests | 108 tests pass | Angular (Szmytkowski) + diagonal radial (Bethe-Salpeter) + κ↔(l,σ) bridge + Kramers-Pasternak direct integration |
| σ·r̂ reduction identity | `(−κ, m_j, −1)` exact | Szmytkowski Eq. 2.7 verification |
| Diagonal ⟨1/r³⟩_{n,l} hydrogenic | Z³/[n³·l(l+½)(l+1)] | T1 closed form, verified to sympy machine precision |
| T1 ⟨1s\|r\|2s⟩ at Z=1 | −32√2/81 | Off-diagonal radial sympy integration sanity |
| `spin_orbit` module tests | 22 tests pass | H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)] closed form |
| SO Kramers cancellation at l=0 | H_SO = 0 exact | κ=−1 forces numerator zero before ⟨1/r³⟩ evaluated |
| SO Z⁴ scaling | (Z/Z_ref)⁴ symbolic | Verified at Z ∈ {1, 3, 4, 38} |
| 2p doublet splitting (Z=1) | α²/32 exact | Breit-Pauli fine-structure benchmark |
| `spin_ful_composed` module tests | 13 tests pass | Full composed-rel pipeline regression |
| LiH/BeH₂/H₂O scalar regression | 334/556/778 Pauli preserved | Bit-exact scalar path unchanged when relativistic=False |
| LiH relativistic Pauli at n_max=1 | exactly 9 | Matches scalar (no spin-orbit at l=0) |
| LiH rel/scalar Pauli ratio at n_max=2 | ∈ [3.7, 4.9] | Pinned at 4.24× in regression suite (post-TR, Sprint 4 v2.15.0) |
| Spinor FCI at α=0 matches scalar FCI | \|ΔE\| < 1e-10 Ha | TR regression test (Sprint 4): jj reduced-matrix-element phase fix |
| Spinor composed block-diagonal ERI | zero cross-block entries | Factorization preserved in relativistic path |
| α → 0 zeroes H_SO diagonal | exact zero | Non-relativistic limit verification |
| `spinor_certificate` module tests | 25 tests pass | Ring R_sp = ℚ(α²)[γ]/(γ²+(Zα)²−1) enforcement |
| Contamination rejection (π, π², ζ(3), log, E₁, unregistered) | raises SpinorTaxonomyError | Six negative controls |
| T3 H_SO block R_sp membership | passes n_max ≤ 4 | Every κ-branch coefficient in ring |
| Sunaga RaH-18q baseline | 47,099 Pauli (published) | Single calibrated cell from PRA 111, 022817 |
| GeoVac rel/Sunaga RaH-18q ratio | 0.011×–0.017× native-Q | Resource advantage verified |
| Fine-structure Li 2²P splitting | sign + OoM correct | Breit-Pauli + Z_eff sanity |
| Fine-structure He 2³P span | sign + OoM correct | 66% relative error (accepted) |
| Fine-structure Be 2s2p ³P span | sign + OoM correct | 78% relative error (accepted) |
| Sprint 5 CP: Li 2²P doublet splitting | < 20% | +8.89% err with std conv (Z_val=1, Z_eff=1) |
| Sprint 5 CP: Be 2s2p ³P span (P₀-P₂) | < 20% | +2.76% err with std conv + Slater Z_eff=1.95 |
| Sprint 5 CP: Be 2s2p ³P individual (3 splittings) | < 20% | all three pass (+18.9%, +6.3%, +2.8%) |
| Sprint 5 CP: 2s2p ↔ 2p² parity-forbidden coupling | exactly 0 | symbolic parity verification |
| Sprint 5 CP: Li core polarization worsens accuracy | err_cp > err_bare | MB α_d/r_c increases Δζ in wrong direction |
| Breit zero at α=0 | breit_eri_count = 0 | α² prefactor kills all Breit terms when alpha_num=0 |
| Breit 1-norm order of magnitude | rel_shift < 1% | Breit 1-norm shift is O(α²) ~ 10⁻⁴ relative to Coulomb |
| Breit Pauli count unchanged | N_pauli(breit) ≤ 2× N_pauli(no breit) | Same Gaunt selection rules; bounded increase |
| Breit block diagonality | cross_block_eri_count = 0 | Breit ERI remains block-diagonal |
| Breit Hermiticity | max(|imag(coeff)|) < 1e-10 | All qubit operator coefficients are real |
| Breit He 2³P splittings | O(10⁻⁵) Ha total | Drake J-pattern physical, nonzero splittings |
| Breit radial Z³ scaling | ratio = 8.0 (Z=2/Z=1) | Retarded integrals scale as Z³ |
| Kramers-Pasternak matches Pochhammer n_r=0 | < 1e-10 | Direct integration reproduces T7 Pochhammer results |
| Kramers-Pasternak matches Hellmann-Feynman ⟨1/r⟩ | < 1e-10 | Direct integration matches all-state HF formula |
| Kramers-Pasternak ⟨r⁰⟩ = 1 normalization | exact | Normalization integral equals 1 for all states |
| Seeley-DeWitt a₀ on unit S³ | √π (exact sympy) | Heat kernel coefficient from Dirac D² spectrum |
| Seeley-DeWitt a₁ on unit S³ | √π (exact sympy) | Curvature correction (R_scalar/6 = 1) |
| Seeley-DeWitt a₂ on unit S³ | √π/8 (exact sympy) | Vacuum polarization coefficient source |
| Vacuum polarization coefficient | 1/(48π²) (exact) | Standard Dirac fermion VP from S³ spectral data |
| β(α) QED one-loop | 2α²/(3π) (exact) | Reproduces standard QED beta function |
| No odd-zeta at one loop | structural theorem | T9 guarantees ζ_{D²}(s) is polynomial in π² |
| D_even(4) decomposition | π²/2 − π⁴/24 − 4G + 4β(4) (PSLQ 80 digits) | Vertex parity exposes Catalan G and Dirichlet β(4) |
| D_odd(4) opposite sign | π²/2 − π⁴/24 + 4G − 4β(4) (PSLQ 80 digits) | Opposite Dirichlet content in odd-n sub-sum |
| D_even(4) + D_odd(4) cancellation | = D(4) to 1e-60 | Dirichlet L-values cancel in full (unrestricted) sum |
| Vertex selection rule consistency | matches hodge1_s3 exactly | SO(4) triangle + parity: n₁+n₂+n_γ odd |
| Fine-structure Dirac formula (n<=4) | exact symbolic | Dirac formula verified for all (n,l,j) through n=4 |
| Dirac accidental degeneracy | 6/6 pairs confirmed | E_FS depends on (n,j) only |
| gamma radial <1/r> NR limit | Z/n^2 | Hellmann-Feynman all-state formula |
| gamma radial <1/r^2> n_r=0 NR limit | Z^2/(n^3(l+1/2)) | Pochhammer ratio |
| gamma radial <1/r^3> n_r=0 NR limit | Z^3/(n^3 l(l+1/2)(l+1)) | Pochhammer ratio |
| zeta_{D^2}(2) | pi^2 - pi^4/12 | Squared Dirac spectral zeta (T9) |
| zeta_{D^2}(s) pi^{even} only | theorem | No odd-zeta content at any integer s |
| Weyl density Dirac S³ | O(1/n_max) convergence | Spectral-to-momentum correspondence |
| Sunset R-scaling power | exact R^10 at (s1,s2,p)=(2,2,1) | Dimensional analysis verification |
| D(s,R)/R^s R-independence | ratio constant to 1e-30 | Transcendental class is R-independent |
| D(5) = 14ζ(3) − 31/2·ζ(5) | exact to 80 dps | ζ(3) structural identification |
| Even-s stays π^{even} at all R | structural theorem | One-loop parity persists |
| Screened radial solver hydrogenic limit | < 1% at Z_eff=const | Reproduces analytical ⟨1/r³⟩ |
| Screened SO enhancement Ca/Sr/Ba | 12×/61×/144× over Z_eff=2 | Core penetration quantified |
| Screened SO CaH/SrH/BaH splitting | within 70% of physical | Leading-order Breit-Pauli |
| NaH balanced FCI overattraction | no equilibrium | Frozen-core cross-V_ne limitation |
| MgH₂ balanced FCI overattraction | no equilibrium | Same mechanism as NaH |
| Self-energy Σ(n_ext=0) structural zero | exactly 0 | Vertex parity selection rule proof |
| Self-energy Σ(n_ext=1) sign | < 0 | Physical (ground state protected, excited states shifted) |
| Self-energy convergence monotonic | monotonic with n_max | Spectral sum convergence |
| Factorized matches direct three-loop | < 1e-25 at n_max=15 | O(N³) vs O(N⁵) algorithmic equivalence |
| Factorized n_max=50 speed | < 60 seconds | O(N³) performance benchmark |
| Factorized convergence monotonic | vals increasing with n_max | Three-loop sum convergence |
| S_min 200 digits verified | 3 independent methods agree | mpmath.nsum, Euler-Maclaurin, direct sum |
| S_min PSLQ irreducibility | 15 failures across 100+ basis | Extended basis including Tornheim-Witten, colored MZV |
| c₂ cross-invariant 8-digit match | |c₂_cross - c₂_apparent| < 2e-8 | Paper 2↔28 bridge verification |
| c₂ symbolic identity | 19/100 - 41π²/25200 (exact sympy) | Rational+π² decomposition |
| c₃ from n_int=0..50 | -5.946(3)×10⁻⁷ at 200.9 sigma | Nonzero; expansion does not terminate at c₂ |
| c₂ T9 consistency | rational + rational·π² only | One-loop π^{even} constraint |
| Paper 28 theorem count | ≥ 4 theorems with proofs | Q-4 exit criterion |
| D₅ c₅(n) exact match n=1..15 | rational equality | Sommerfeld closed-form verification |
| D₅ PSLQ decomposition | residual < 1e-190 | Weight-9 MZV decomposition |
| D₅ ζ(2)ζ(7) cancellation | coefficient = 0 | Product survival rule verification |
| D₆ analytical assembly | 62 matching digits | Weight-11 MZV decomposition from 7 Euler sums |
| D₆ ζ(2)ζ(9) cancellation | coefficient = 0 | Product survival rule verification at p=6 |
| D₆ surviving products | exactly 3 | z(3)z(8), z(4)z(7), z(5)z(6) as predicted |
| Product survival rule | max(0, floor((2p-5)/2)) | Verified D₂..D₆ |
| K/π not in Q-span of D₂..D₆ | PSLQ null | K-Sommerfeld structural separation |
| Graph-native QED tests (non-slow) | 72 pass | GN-5 self-energy + vertex pipeline |
| Graph-native Σ trace (t=0, n_max=2) | exactly 44/3 | Self-energy trace validation |
| Graph-native Λ trace (t=0, n_max=2) | exactly 32/9 | Vertex correction trace validation |
| Graph-native F₂ (t=0, n_max=2) | exactly 5√2/3 | Anomalous moment (irrational, π-free) |
| Graph-native Σ ground-state block ≠ 0 | [[1,1],[1,1]] | Broken structural zero (CG opens couplings) |
| Graph-native Σ eigenvalues (t=0) | {0(×5), 4/3, 2(×2), 4, 16/3} | Self-energy spectrum |
| Graph-native VP Π trace (n_max=2) | exactly 32 | Vacuum polarization trace (GN-4) |
| Graph-native VP Π trace (n_max=3) | exactly 3192/5 | VP extension (GN-6) |
| Graph-native projection C (n_max=3) | 50471424/1779441125 | Rational projection exchange constant (GN-7) |
| Graph-native all π-free | structural | All graph QED quantities algebraic, no transcendentals |
| GN-6 n_max=3 VP tests | 45 pass | Extended VP at n_max=3 |
| GN-7 continuum bridge tests | 63 pass | Projection exchange constant verification |
| F₂(t) even-function | c₁ = c₃ = c₅ = 0 exact | Graph-native anomalous moment parity |
| F₂(κ) convergence | 2.353, 1.873, 1.589 | Monotonic decrease at n_max=2,3,4 |
| F₂(κ) at n_max=5 | 1.39581063 | Extended convergence (numpy) |
| F₂(κ) at n_max=6 | 1.25321124 | Extended convergence (numpy) |
| F₂ power-law exponent | -0.573 (R²=0.99990) | F₂ ~ 3.507 × n^(-0.573) |
| Pendant-edge n_max=5 | Σ(GS) = 1.60000000 = 8/5 | Exact match 2(n-1)/n |
| Pendant-edge n_max=6 | Σ(GS) = 1.66666667 = 5/3 | Exact match 2(n-1)/n |
| F₂(t) n_max=3 rational degree | 16/16 | Even function over ℚ(√2,√3,√5,...) |
| F₂ successive ratios → 1 | 0.796, 0.849, 0.878, 0.898 | Power-law consistency |
| Selection rule census | 1/8 survives | Only Gaunt/CG sparsity survives on graph |
| C × F₂ divergence | grows with n_max | VP projection ≠ vertex projection |
| Dirac graph Rule B selection rules | 4/8 survive | Spinor-recoverable vs vector-required partition |
| Dirac graph GS NOT pendant | degree 2 (A), 5 (B) | Constant across n_max=2,3,4 |
| Dirac graph Σ(GS) ≠ 0 | nonzero both rules | Structural zero not recovered |
| Dirac graph Σ strictly PSD | 0 zero eigenvalues | Unlike scalar graph (5 zeros) |
| Dirac graph π-free | ℚ[√2,√17,√41,√881] | Algebraic, no transcendentals |
| Dirac graph Tr(Σ_A) | exactly 17/2 | Exact algebraic trace |
| Dirac graph Tr(Σ_B) | exactly 103/20 | Exact algebraic trace |
| Dirac graph Furry recovered | tadpole = 0 | Off-diagonal identity vertex parity |
| VP/SE projection ratio constant | CV < 1% across n_max=3,4,5 | Two-tier calibration structure |
| VP/SE ratio ≈ 3/29 | 0.1035 ± 0.0009 | Closest simple rational |
| C_F2 exponent matches F₂ | +0.57 vs −0.57 | Complementary scaling |
| Catalan/β from graph parity | 24 PSLQ null | Clean negative |
| Per-mode projection topological | ρ varies 0 to 1.66 | Not multiplicative |
| Transverse Σ_T(GS) = 0 (n_max=3) | < 1e-16 | Plaquette-recovered GS protection |
| Transverse Σ_T(GS) = 0 (n_max=4) | < 1e-16 | Plaquette-recovered GS protection |
| Co-exact eigenvalue q=3 match (n_max=3) | exact 15.0 = 3×5 | Transverse mode spectrum |
| Triangle inequality modes q=3,4 (n_max=4) | 0% violations | Low-q mode-resolved SR test |
| Triangle inequality modes q=5,6 (n_max=4) | 63-86% violations | High-q structural negative |
| Ward identity ||[Σ_T,H₀]||/||Σ_T|| | ~0.58 | Non-diagonal (calibration needed) |
| Vector QED scalar+vector 7/8 rules | 7 of 8 PASS | Vector photon selection rule recovery |
| Vector QED Dirac+vector 8/8 rules | 8 of 8 PASS | Dirac spinor phase constraint recovers Furry |
| Vector QED GS structural zero | Σ(GS) = 0 exactly | Angular momentum GS protection |
| Vector QED l=1 self-energy at n_max=2 | exactly 1/(4π) | S² Weyl calibration exchange constant |
| Vector QED Furry scalar FAIL | nonzero odd-q diagonal | Scalar electrons: no spinor phase constraint |
| Vector QED Furry Dirac PASS | tadpole = 0 exactly | Dirac spinor phase constraint: V(a,a)=0 |
| Vector QED scalar ≠ Dirac rules | 7/8 vs 8/8 | Dirac spinor phase recovers Furry |
| L₁ non-block-diag by (q,m_q) | cross-sector 20% | Photon QN are transitions, not conserved |
| Vector QED tests | 99 pass | Full vector QED pipeline |
| Operator system tests | 24 pass | Truncated O_{n_max} construction, *-closure, witness pair, propagation number |
| prop(O_{n_max}) at n_max=2,3,4 | exactly 2 | Matches Connes-vS Toeplitz S¹ Prop 4.2 verbatim |
| dim(O_{n_max}^2) = N² at n_max=2,3,4 | exact | Operator system squared spans full envelope (envelope-collapse verified) |
| O_{n_max} witness pair (n_max=2) | residual ≥ 1e-2 | M_{2,1,0}² has 14.9% residual against O — exhibits ab ∉ O concretely |
| 𝟙 ∈ O_{n_max} | residual ≤ 1e-13 | Constant function f=1 acts as unit of O |
| Circulant-S³ tests | 35 pass | Falsification comparator: A_circ_N is commutative C*-subalgebra, prop=1 intrinsic / ∞ ambient |
| prop(A_circ_N) intrinsic | exactly 1 | Verified at N ∈ {2, 5, 14, 30, 55, 140} — categorically different from GeoVac prop=2 |
| prop(A_circ_N) ambient | infinity | Commutative subalgebra cannot generate non-commutative M_N(ℂ) |
| Connes distance tests | 17 pass | SDP framework, gauge fixing, distance matrix, m-reflection forced zeros |
| Connes distance d(v,v) = 0 | exact | Diagonal positivity (sanity check) |
| Connes distance triangle inequality | holds | d(u,w) ≤ d(u,v) + d(v,w) verified at n_max=2,3 |
| Connes distance m-reflection | d(|n,l,-m⟩,|n,l,+m⟩)=0 | SO(3) m-reflection symmetry — physical, not pathology |
| Connes distance n_max=2 ratios | 1:2:3 with 50√3 | Algebraic fingerprint: 28.87 : 57.73 : 86.60, 86.6020 = 50√3 to 5dp |
| Speed regression | < 10% | Performance control |
| Paper 2 cubic root 1/α | 137.036011 (above CODATA, 8.8e-8) | Table II/abstract value correction (mpmath) |
| Paper 2 §VIII.D identity (k=0) | B_formal/N=d, root m=2 | Summation-index correction (sympy) |
| Propinquity bound γ_{n_max=2,3,4} | 2.0746 / 1.6101 / 1.3223 (system-independent) | Paper 38 main theorem evaluation on `GeoVacHamiltonian.propinquity_bound` metadata; first quantitative basis-truncation error estimate at chemistry-consumer API (Target A, v3.85.0) |
| FCIDUMP LiH composed n_max=2 round-trip | max\|h_1 diff\| = 0.0, max\|eri diff\| = 0.0 | `to_fcidump()` / `read_fcidump()` bit-exactness (P1, v3.85.0) |
| FCIDUMP 7-system sample (LiH/BeH₂/H₂O/NaH/KH/MgH₂/CaH₂) | < 1e-10 relative | Multi-system round-trip integrity (P1, v3.85.0) |
| LiH composed n_max=2 R=3.015 bohr qubit FCI | −14.143 Ha at Q=30, 333 Pauli | Publication-grade headline from FCIDUMP→DMRG path (R3-A, v3.85.0) |
| LiH composed/balanced qubit FCI PES | monotone-descending across R ∈ [2.5, 5.0] bohr | Operational confirmation of W1e wall at projection step (R3-A, v3.85.0) |
| NaH balanced n_max=3 W1e over-binding at R_e^exp | +1.50 Ha vs experimental D_e = 0.0713 Ha (21× over) | W1e baseline (P4, v3.85.0) |
| NaH balanced DMRG vs P4 direct-FCI baseline | bit-identical at every R (max diff ~6e-13 Ha) | DMRG-falsifier STOP verdict (R3-B, v3.85.0) |
| H₂ Q=10 composed n_max=2 R=1.4 openfermion UCCSD | error 8.6e-13 mHa (165 evals, 0.92 s) | Openfermion-native UCCSD vs L-BFGS-B from HF init (R3-C, v3.85.0); 5 OoM tighter than published STO-3G UCCSD literature |
| Openfermion vs qiskit-nature UCCSD at H₂ Q=10 | ~3000× speedup per energy evaluation | Empirical validation that openfermion is the correct VQE stack for GeoVac (R3-C, v3.85.0) |

---
