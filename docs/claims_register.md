# GeoVac Claims Register

**Purpose:** one line of verification truth per headline claim, for external readers and internal discipline. This register makes the project's internal status vocabulary externally legible. It is maintained alongside the papers (single source of truth: this file + the cited tests; history in git/CHANGELOG).

**Status vocabulary** (strongest → weakest):

| Status | Meaning |
|---|---|
| SYMBOLIC PROOF | Exact symbolic derivation, machine-verified (sympy/exact rationals), frozen in `tests/` |
| MEASURED | Computed result vs. an external reference value; number and reference stated |
| PANEL-VERIFIED | Exact/bit-exact at finite cutoffs on a stated panel; general-cutoff statement is induction from the panel unless noted |
| INTERNAL THEOREM | Proof written and internally checked; **no external expert review** |
| CONDITIONAL | Internal theorem with explicitly named gaps; gaps listed |
| OBSERVATION | Numerical/structural coincidence; no derivation claimed |
| CONJECTURE | Stated as such |
| RETRACTED | Withdrawn claim; retained here so old DOI'd copies resolve to the truth |

**No claim in this register is externally peer-reviewed as of June 2026.** The project's dissemination is GitHub + Zenodo (DOI-stamped releases).

---

| # | Claim | Where | Status | Verification / falsifier |
|---|---|---|---|---|
| 1 | Fock 1935 conformal equivalence: the discrete S³ graph Laplacian construction reproduces hydrogenic structure (18 symbolic proofs) | Paper 7 | SYMBOLIC PROOF | `tests/test_fock_projection.py`, `tests/test_fock_laplacian.py` (run before every release) |
| 2 | Unit-S³ Laplace–Beltrami eigenvalues = −(n²−1) (continuum operator, symbolic n≤3). Separately: κ = −1/16 is the energy-matching scale, which *coincides* with the geometric 1/16 (Fock coupling c²(n,0), inverse Jacobian 1/Ω⁴(0)) — a coincidence, NOT a derivation | Paper 7; Paper 0 | eigenvalues **SYMBOLIC PROOF**; **κ-derivation OBSERVATION** (corrected 2026-06-14 — was mislabeled "derived / SYMBOLIC PROOF") | `test_fock_laplacian.py` (eigenvalues); `test_trunk_qa_kappa.py` (κ — counterfactual shows no geometric→matching bridge) |
| 3 | Angular sparsity theorem: two-electron-integral density depends only on l_max, not on the potential V(r) | Paper 22 | INTERNAL THEOREM + PANEL-VERIFIED | verified across Coulomb/HO/Woods-Saxon; 1.44% at l_max = 3 |
| 4 | Composed qubit Hamiltonians: N_Pauli ≈ 11.11 × Q across 40 molecules; O(Q^2.5) scaling; 51×–1712× fewer Pauli terms than Gaussian-basis baselines at matched qubit counts | Papers 14, 20 | MEASURED | library benchmarks; **caveat co-stated:** comparison is at matched qubits, not matched accuracy (see row 5) |
| 5 | Chemistry accuracy ceilings: LiH R_eq 5.3% (composed), H₂O R_eq 26%, NaH does not bind; the wall (W1e) is at the Hamiltonian-specification level and is structurally tied to the sparsity advantages | Papers 17, 19; CLAUDE.md §3 | MEASURED (honest negative) | DMRG=FCI falsifier sprints (R3-A/B); F1–F6 + kwarg-sweep dead-end records |
| 6 | He ground state: 0.019% (2D variational + cusp), 0.20% (zero-parameter graph-native CI) | Paper 13 | MEASURED | vs. NIST reference |
| 7 | Hydrogen Lamb shift at one loop: −0.534% vs. experiment, no fitted parameters | Paper 36 | MEASURED | 1052.19 MHz vs 1057.845 MHz; two-loop NOT closed (renormalization gap, LS-8a) |
| 8 | State-space GH convergence of truncated Camporesi–Higuchi (Dirac) triples on S³, rate (4/π)log n/n, translation-seminorm metrization | Paper 38 | INTERNAL THEOREM (unconditional) | Former gaps closed 2026-06-10: kernel condition proved on the truthful substrate (Schur multiplicity-one + per-band injectivity; `tests/test_p38_action_seminorm.py`); dual direction via exact-fit spinor lifted state — both almost-inverse defects are Fejér smoothings at γ_n, no inverse estimate remains. Scalar-sector convergence for compact groups is **prior art** (Gaudillot-Estrada–van Suijlekom, arXiv:2310.14733); contribution = spinor transport + per-band injectivity + explicit 4/π rate |
| 9 | 4/π rate constant universal across compact Lie groups (dual-Coxeter normalization) | Paper 40 | INTERNAL THEOREM | inherits row 8's L2 correction; rate moment is a kernel property, unaffected |
| 10 | K⁺ annihilation theorem: compressing a Krein-self-adjoint product Dirac to the Krein-positive subspace annihilates the spatial Dirac; the compressed Lipschitz seminorm vanishes identically on the momentum-diagonal operator system | Paper 45 | INTERNAL THEOREM + PANEL-VERIFIED (bit-exact) | `tests/test_p45_kplus_degeneracy.py`; `debug/p45_kplus_seminorm_check.py` (42/42, 275/275 kernel) |
| 11 | "First Lorentzian propinquity convergence theorem on truncated Krein spectral triples" | formerly Paper 45 | **RETRACTED** (2026-06-09) | the bounded quantity was degenerate (row 10); imported theorem numbers were nonexistent; see Paper 45 History remark + CHANGELOG v3.106.0 |
| 12 | Strong-form / bridge / hybrid Lorentzian metric claims built on the row-11 quantity (Λ^strong equality, pLGH panel transport, G2-metric closure) | Papers 46–49 | RETRACTED/DESCOPED (2026-06-09) | Status notes in each paper; norm-resolvent outer arrow of Paper 47 and the cocycle/deficit algebra of Paper 49 survive as stated there |
| 13 | χ₋₄ identity: D_even(s) − D_odd(s) = 2^{s−1}(β(s) − β(s−2)) for the S³ Dirac parity sectors; D(4) = π² − π⁴/12 | Paper 28 Thm 3 | SYMBOLIC PROOF | 3-line Hurwitz proof; verified to 40 digits vs direct spectral sums (proof display corrected 2026-06-10: factor-2 convention) |
| 14 | ζ_{D²}(−k) = 0 for all k ≥ 0 on unit S³ (Bernoulli mechanism) ⟹ spectral action exactly two-term; scalar heat trace a_k = 2π²/k! all orders | Paper 51 | SYMBOLIC PROOF | `tests/test_paper51_*.py` |
| 15 | GeoVac Hopf graphs are Ramanujan — as a finite-size statement; three of four families cross the bound at V ≈ 30–60 | Paper 29 | MEASURED | Alon-Boppana sweep data; NOT an asymptotic claim |
| 16 | α observation: K = π(B + F − Δ) matches 1/α to 8.8×10⁻⁸; B, F, Δ each have derived spectral homes | Paper 2 | OBSERVATION | combination rule **conjectural by standing policy**; 12 single-mechanism derivations eliminated; no derivation claimed |
| 17 | Master Mellin engine / π-source case-exhaustion: every π entering the framework's projections arises from Tr(D^k e^{−tD²}) at k ∈ {0,1,2} | Papers 18, 32 §VIII | INTERNAL THEOREM | case-checked over the first 15 of 28 Paper 34 projections (§III.1–15; 16–28 pending the theorem's extension); status caveat in Paper 32 §VIII applies to the GH-theorem entry |
| 18 | Tannakian substrate: Aut^⊗(ω) reconstruction closed at finite cutoff (2,960 zero residuals) | Paper 56 | PANEL-VERIFIED | exact-rational panels at n_max ≤ 6; infinite-cutoff identification with motivic Galois is open (multi-year), **not claimed** |
| 19 | Quantum-simulation positioning: structural sparsity is basis-intrinsic and survives tapering/grouping; market test = 1-norm parity with STO-3G (0.97×), 2.7× fewer Pauli | Papers 14, 20 | MEASURED | honest framing: parity with the weakest standard basis; accuracy caveats of row 5 apply |
| 20 | Negative-results corpus: 100+ documented dead ends (PK modifications, cusp treatments, single-center molecular encodings, relativistic-Z₂ tapering, K⁺ compression, …) | CLAUDE.md §3; papers | MEASURED (negatives) | each row carries a memo path; guardrail papers (8–9, FCI-M) protect against re-derivation |
| 21 | Product-carrier convergence on S³×S¹_T in the action-seminorm framework: truncated joint system → C(S³×S¹_T) in vS state-space GH at additive rate γ_n^{SU(2)} + γ_K^{U(1)}; the temporal direction carries genuine metric weight (Toeplitz algebra: L(S_q) = 2πq/T exact at every finite window; the 2026-06-09 invisibility was the momentum-diagonal algebra, not time) | Paper 45 `prop:product_action_seminorm` (2026-06-10) | INTERNAL THEOREM + PANEL-VERIFIED | `tests/test_wh7_b1_joint.py`, `tests/test_wh7_toeplitz_temporal.py`; **signature-agnostic** (Wick-rotated/thermal carrier) — explicitly NOT a Lorentzian claim; row 11 stays RETRACTED; the Lorentzian target is reopened Q1 (boost/modular-flow seminorm candidate) |

---

*Maintained since 2026-06-10 (Phase 2c of `docs/corpus_accessibility_plan.md`). When a claim's status changes, update the row and record the change in CHANGELOG.md; do not fork the register.*
