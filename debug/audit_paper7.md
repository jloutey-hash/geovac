# Confidence Audit: Paper 7 — The Dimensionless Vacuum: Recovering the Schrödinger Equation from Scale-Invariant Graph Topology

## Calibration check
Not a calibration run. Cross-corpus check (Paper 18 §sec:kappa, §sec:kappa_derivation; Paper 0 §kinetic-scale; Paper 13 §sec:casimir; Paper 1 §discretization-artifacts) was applied before each finding was promoted to a defect. Two candidate findings were demoted to "lags-corpus" framing after cross-checking against Paper 18, and one numerical inconsistency in §VI.B survived cross-check as a genuine math error.

## Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | Fock stereographic projection sends ℝ³ to unit S³; ∑n_i²=1 for all p₀>0 | Abstract, §IV.B, Eq. (eq:unit_sphere) | A | EXTERNAL (Fock 1935) | sympy verified ∑n_i² = 1 identically |
| 2 | Chordal distance identity ‖n−n'‖² = Ω·Ω'·‖p−p'‖² | §IV.C, Eq. (eq:chordal) | A | EXTERNAL (standard Fock projection algebra) | sympy LHS−RHS = 0 identically |
| 3 | Conformal Laplacian formula Δ_{S³}f = Ω⁻²[Δ_flat f + ∇(ln Ω)·∇f] in n=3 | §V.A, Eq. (eq:conformal_laplacian), Appendix B | A | EXTERNAL (standard Riemannian-geometry result) | Appendix B derivation reproduces standard Wald/Lee Christoffel transform; numerically verified Δ_{S³}(cos χ) = −3 cos χ and Δ_{S³}(4cos²χ − 1) = −8(4cos²χ−1) |
| 4 | S³ Laplace-Beltrami spectrum λ_n = −(n²−1) (pure integers) | §V.B–C, Eq. (eq:s3_spectrum) | A | EXTERNAL (Bargmann 1936; standard sphere harmonics) | sympy verified n=2,3 directly via Gegenbauer C¹_{n−1}; index matches U_{n−1} Chebyshev second-kind identity (sympy C¹_{n−1} − U_{n−1} = 0) |
| 5 | "We demonstrate that the GeoVac discrete graph Laplacian converges to Δ_{S³} in the continuum limit" | Abstract, §I.3rd-paragraph, §III.C, §VII.A Stage 1 | **C** | GEOVAC-ONLY for the operator convergence claim | The paper provides ONE convergence indicator (the 2s/2p splitting decay to <0.01% at n_max=30) and an SO(4) isometry identification. This is empirical spectral evidence, not a proven operator-convergence theorem. The corpus contains a real theorem of this kind (Paper 38, WH1 PROVEN), but Paper 7 does not cite it and does not match that level of rigor. Abstract overstates: "We demonstrate ... converges" → suggest "We provide spectral evidence (degeneracy recovery) consistent with convergence". |
| 6 | "κ = −1/16 is the unique constant that maps the dimensionless graph eigenvalues to the physical hydrogen spectrum" | §III.A Eq. (eq:hamiltonian), §VII.D | **C** | GEOVAC, but cross-corpus-corrected | Cross-corpus check (Paper 18 §sec:kappa line 167): "Per-shell energies E_n = −Z²/(2n²) emerge from the full graph Hamiltonian H = κℒ + W (with node weights W_nn = −Z/n²), not from κ alone." Paper 7 states H = κL with no node weights, then asserts κ alone produces the Rydberg spectrum. **Paper 7 LAGS the corrected version in Paper 18.** A single constant cannot map −(n²−1) onto −1/(2n²) for all n simultaneously (different functions of n). Suggest acknowledging W in §III.A or adding a footnote pointing to Paper 18 §sec:kappa. |
| 7 | The free SO(6) Casimir eigenvalues for He (S⁵) are "ν(ν+4)/2 = 0, 3, 8, 15, …" | §VI.B (immediately after Table I) | **E** | Math error | ν(ν+4)/2 for ν = 0,1,2,3 gives 0, 5/2, 6, 21/2 — NOT 0, 3, 8, 15. The numbers 0, 3, 8, 15 equal n²−1 for n=1,2,3,4 (the S³ spectrum from §V.B). **Doubly wrong:** even if quoted correctly as ν(ν+4)/2, the ¹S singlet (Paper 7's stated case) restricts to even ν only, giving 0, 6, 16, 30, 48 — see Paper 13 Eq. (eq:mu_free) and surrounding text. Suggested replacement: "the free eigenvalues ν(ν+4)/2 = 0, 5/2, 6, 21/2, …, with ¹S singlet selecting even ν only (0, 6, 16, 30, …)". |
| 8 | F⁰(1s,1s) = 5Z/8; F⁰(1s,2s) = 17Z/81; F⁰(2s,2s) = 77Z/512 | §VI Eqs. (eq:slater_exact), (eq:phi1s)–(eq:phi2s), three verification subsubsections | A | EXTERNAL (Slater 1929; standard atomic-structure result) | sympy reproduced all three rational values from the t-substitution exactly; the algebraic chain (q = 2Zt → Φ_a(t) → ∫₀^∞ Φ_a Φ_b dt × 4Z/π) is internally tight and gives bit-exact rationals |
| 9 | κ = 1/Ω⁴(0) with Ω(0) = 2 | §VII.D | **C** (notational sloppiness, not error) | GEOVAC + standard projection | Test 5 in Paper 7's own appendix correctly says "Ω(0) = 2/p₀". §VII.D drops the p₀-dependence and writes "Ω(0) = 2", which is only true at p₀=1 (the n=1 ground state). The identification is correct when read with p₀=1 implicit, but as written the equation Ω(0)=2 is inconsistent with the rest of the paper. Suggest "Ω(0) = 2/p₀ = 2 evaluated at p₀=1 (ground state)". |
| 10 | c²(n,l) = (1/16)[1 − l(l+1)/(n(n+1))], universal at l=0, c²(4,3) = 1/40 | §VII.D Eq. (eq:fock_coupling) | A | Internal arithmetic, externally checkable | Direct evaluation confirms 1/40 for (n,l)=(4,3). The "1/40 = Δ from Paper 2" connection in the same paragraph rests on Paper 2 (which CLAUDE.md §13.5 keeps explicitly conjectural at the combination-rule level) — this side claim is GEOVAC-ONLY interpretation, separate from the c² formula itself. |
| 11 | He graph-native FCI at n_max=7: −2.8983 Ha, 0.19% error vs exact −2.903724 Ha | Appendix "Note added v2.6.0" | A | EXTERNAL reference (Drake / Pekeris-style high-precision He value −2.903724 Ha) | Recomputed: 100·\|−2.8983 − (−2.903724)\|/\|−2.903724\| = 0.187%. Paper's 0.19% matches. |
| 12 | 18 symbolic proofs verify "every algebraic step of the conformal transformation chain" | Abstract; §II.B item 2; §VII.A Stage chain; Appendix A | B | GEOVAC tests of GEOVAC code, but verifying STANDARD identities | The 18 items in Appendix A all verify externally-established results (Fock projection, conformal Laplacian, Gegenbauer eigenvalues). Re-derivable trivially from textbook math; tests are correct but they check known identities, not new theorems. Honest framing already present in §II.B: "Our symbolic proofs in the Appendix verify these identities; they do not claim to derive them for the first time." Paper 7 itself is appropriately careful here; flag only because abstract reads more strongly. |
| 13 | Bargmann lithium ²S claim: spatial irrep [2,1] of S₃ first appears at ν=1, lowest free eigenvalue μ_free = 1·8/2 = 4 | §VI.C | D | GEOVAC + group theory | The arithmetic ν(ν+7)/2 at ν=1 = 4 is consistent with Paper 7's own Table I formula. The S₃ representation-theory claim (mixed-symmetry [2,1] first appears at ν=1) is a standard hyperspherical harmonic decomposition result; not verified here. Hand to a domain expert if needed. |

## Numbers I recomputed

| Claim | Paper's figure | Independent reference | Recomputed value/error | Survives? |
|---|---|---|---|---|
| F⁰(1s,1s) | 5Z/8 | Slater 1929; Condon-Shortley | sympy exact 5Z/8 | YES |
| F⁰(1s,2s) | 17Z/81 | Standard atomic-structure tables | sympy exact 17Z/81 | YES |
| F⁰(2s,2s) | 77Z/512 | Standard atomic-structure tables | sympy exact 77Z/512 | YES |
| Δ_{S³}(cos χ) | −3 cos χ | Standard S³ Laplace-Beltrami spectrum (textbook) | sympy bit-exact | YES |
| Δ_{S³}(4cos²χ−1) | −8(4cos²χ−1) | Standard | sympy bit-exact | YES |
| ν(ν+4)/2 at ν=0,1,2,3 | "0, 3, 8, 15" (Paper 7 §VI.B) | Definition + Paper 13 Eq. (eq:mu_free) | 0, 5/2, 6, 21/2 (or 0, 6, 16, 30 for even-ν singlet) | **NO** |
| He n_max=7 error | 0.19% | Drake/Pekeris −2.903724 Ha | 0.187% | YES |
| c²(4,3) | 1/40 | Paper 7's own formula | 1/40 exact | YES |

## Circularity map

Claims grounded EXTERNAL:
- Fock projection algebra (Claim 1, 2, 4): external math
- Conformal Laplacian formula (Claim 3): external Riemannian geometry
- Slater F⁰ values (Claim 8): external atomic-structure tables
- He reference energy −2.903724 Ha (Claim 11 denominator): external

Claims grounded GEOVAC-ONLY (require care if quoted as proof outside the corpus):
- **Graph Laplacian operator convergence to Δ_{S³}** (Claim 5): supported here only by the 2s/2p splitting decay; the rigorous operator-convergence theorem lives in Paper 38 and is not cited from Paper 7. If a reader pulled Paper 7 alone, the "we demonstrate convergence" claim is empirically supported, not proven.
- **κ as "the unique constant"** (Claim 6): the corrected version in Paper 18 says κ + node weights W are both needed. Paper 7 LAGS Paper 18 here.
- **κ = 1/Ω⁴(0) with Ω(0)=2** (Claim 9): valid only at p₀=1; presented without the p₀=1 specification.
- **c²(4,3) = 1/40 = Δ in Paper 2** (Claim 10 footnote): the connection to Paper 2 is GEOVAC-internal and inherits Paper 2's "conjectural" status (CLAUDE.md §13.5).

Claim 12 (the 18 symbolic proofs): MIXED. Internal tests of internal code, but the identities being tested are externally established. Paper 7 itself acknowledges this in §II.B. Read in context the framing is honest; only the abstract's "we verify these claims through 18 independent symbolic proofs" reads slightly louder than the body warrants.

## Overstatement findings

1. **Abstract**: "We demonstrate that the GeoVac discrete graph Laplacian converges to this same operator in the continuum limit, establishing a precise mathematical bridge"
   → Suggested replacement: "We provide spectral evidence (degeneracy recovery, eigenvalue identification) that the GeoVac discrete graph Laplacian is consistent with convergence to this operator; a rigorous Gromov-Hausdorff propinquity proof for the related Camporesi-Higuchi spectral triple appears in Paper 38 [cite]."

2. **§III.A**: "κ = −1/16 is the kinetic scale factor---the unique constant that maps the dimensionless graph eigenvalues to the physical hydrogen spectrum"
   → Suggested replacement: "κ = −1/16 is the kinetic scale factor that calibrates the graph-Laplacian energy unit. The full per-shell hydrogen spectrum E_n = −Z²/(2n²) emerges from H = κℒ + W with node weights W_nn = −Z/n² (Paper 18 §III.A); κ alone calibrates the overall scale."

3. **§VI.B** (math error, not overstatement): "the free eigenvalues ν(ν+4)/2 = 0, 3, 8, 15, ..."
   → Required replacement: "the free eigenvalues ν(ν+4)/2 = 0, 5/2, 6, 21/2, 16, ..., with the ¹S singlet selecting even ν only and giving 0, 6, 16, 30, 48, ... (Paper 13 §sec:casimir)". (The "0.001%" numerical agreement claim immediately following also needs to be checked against the corrected list — Paper 13 confirms the even-ν singlet values 0, 6, 16, 30 are what is actually verified numerically.)

4. **§VII.D**: "$1/16 = 1/\Omega^4(0)$ where $\Omega(0) = 2$ is the stereographic conformal factor at $p = 0$"
   → Suggested replacement: "$1/16 = 1/\Omega^4(0)$ evaluated at the ground state $p_0 = 1$, where $\Omega(0) = 2/p_0 = 2$." (The current wording elides the p₀=1 specialization; cf. Test 5 in Appendix A, which correctly writes Ω(0) = 2/p₀.)

## Broadcast readiness: YELLOW

Paper 7 is a foundational document whose mathematical core (the Fock projection algebra, the conformal Laplacian eigenvalue chain, the S³ density-overlap Slater integrals) is solid — every symbolic claim that the paper makes about continuous identities reproduces bit-exact in sympy, and the 18 symbolic proofs do what they advertise. The He graph-native FCI number checks against a standard external reference. The Slater F⁰ values are all standard and reproduce exactly via the paper's own t-substitution.

The defect that prevents a clean GREEN is the §VI.B arithmetic error: the listed numbers "0, 3, 8, 15" are not ν(ν+4)/2 and are not the ¹S singlet sequence either. Because Paper 7 is the most-cited paper in the corpus (31 inbound citations) and §VI is the entry point for the N-electron generalization, the error propagates to any reader trying to learn the framework from Paper 7. This is a HIGH-severity, isolated correction — one number-list edit fixes it.

Two MEDIUM items both have the same shape: the abstract and §III.A read more strongly than the body actually supports, and the corpus already contains the corrected version. Paper 7 LAGS the corrected version (Paper 18 §sec:kappa for κ + W; Paper 38 for the actual GH-convergence theorem) without flagging that lag. These are paper-text-only fixes (no new physics), but they matter because Paper 7 is the load-bearing foundation.

LOW items: §VII.D's "Ω(0) = 2" notation, Appendix A Test 12 phrasing, and one or two stylistic items.

Net: with the §VI.B fix and two suggested abstract/§III.A softenings (plus the §VII.D notation tightening), Paper 7 is broadcast-safe. Until those are applied, an expert reader sees a math error on the He generalization, an over-strong "we demonstrate convergence" in the abstract, and an under-stated κ derivation in §III.A.

## What I could NOT verify (hand to a human expert)

- **Claim 13** (lithium [2,1] irrep first appearing at ν=1): The arithmetic ν(ν+7)/2 at ν=1 = 4 is consistent with Paper 7's own Table I. The S₃ rep-theory claim that [2,1] first appears at ν=1 in the hyperspherical-harmonic decomposition is a standard but non-trivial group-theory statement; not re-derived here. Hand to a domain expert if it matters for any downstream paper.
- **The actual operator-norm convergence rate** of L_graph → Δ_{S³}: Paper 7 only demonstrates spectral convergence of one eigenvalue gap. The full operator-norm or resolvent-norm convergence on the Hilbert space lives in Paper 38, which I treated as a separate document and did not cross-verify here.
- **The "FCI basis invariance < 3×10⁻¹⁵ Ha" claim** in the Appendix Note (v2.6.0): unverifiable from .tex; rests on the cited Sprint 3D test, which the citation checker can confirm exists.

## Summary line
**Title:** Paper 7 — The Dimensionless Vacuum (S³ proof). **Verdict:** YELLOW. **A/B/C/D/E:** A=6, B=1, C=4, D=1, E=1. **HIGH/MEDIUM/LOW:** HIGH=1, MEDIUM=2, LOW=2. **Top finding:** §VI.B states "the free eigenvalues ν(ν+4)/2 = 0, 3, 8, 15" — the listed numbers are wrong: ν(ν+4)/2 gives 0, 5/2, 6, 21/2 for ν=0..3 (and 0, 6, 16, 30 for the ¹S even-ν restriction Paper 7 explicitly cites). The numbers 0, 3, 8, 15 are the S³ spectrum n²−1, mistakenly transcribed into the S⁵ row. Single-line fix.
