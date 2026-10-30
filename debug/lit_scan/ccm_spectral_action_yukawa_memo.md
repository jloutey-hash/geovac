# CCM spectral-action Yukawa constraint, and whether a truncated S³ outer factor moves it

Date 2026-10-29. Read-only literature + corpus scan. No computation run beyond text extraction of downloaded PDFs (`python -I`). Downloads: `debug/lit_scan/dl_ccm_yukawa/{ccm2007,cc2012,ccvs2013}/` (untrusted data; text extracted into the session scratchpad).

Sources read in full text (PRIMARY):
- [CCM07] Chamseddine, Connes, Marcolli, "Gravity and the standard model with neutrino mixing", Adv. Theor. Math. Phys. 11, 991-1089 (2007), arXiv:hep-th/0610241.
- [CC12] Chamseddine, Connes, "Resilience of the spectral standard model", JHEP 1209 (2012) 104, arXiv:1208.1030.
- [CCvS13] Chamseddine, Connes, van Suijlekom, "Beyond the spectral standard model: emergence of Pati-Salam unification", JHEP 1311 (2013) 132, arXiv:1304.8050.

Abstract only (PRIMARY, abstract page):
- [CCvS15] same authors, "Grand unification in the spectral Pati-Salam model", JHEP 1511 (2015) 011, arXiv:1507.08161.
- [CC08] Chamseddine, Connes, "The uncanny precision of the spectral action", Commun. Math. Phys. 293, 867-897 (2010), arXiv:0812.0165.

---

## 1. The relation (CCM07, verified against the text)

**Inner traces.** With the finite Dirac's Yukawa matrices Y_(↑3), Y_(↓3) (up/down quarks), Y_(↑1), Y_(↓1) (neutrino/electron) and the Majorana matrix Y_R, CCM07 defines (eq. 3.32 and §3.7.1)

    a = Tr( Y*_(↑1)Y_(↑1) + Y*_(↓1)Y_(↓1) + 3( Y*_(↑3)Y_(↑3) + Y*_(↓3)Y_(↓3) ) ),
    b = Tr( (Y*_(↑1)Y_(↑1))² + (Y*_(↓1)Y_(↓1))² + 3(...)² ),
    c = Tr(Y*_R Y_R),   d = Tr((Y*_R Y_R)²),

with e a mixed neutrino/Majorana trace. The factor 3 is colour. (CCM07 gives a explicitly at eq. 3.32 via Tr(M(φ)²) = 4a|φ|². I took b's form from the quartic term, eq. 3.51, f0/(8π²)·Tr((D_(0,1))⁴) = f0/(2π²)·b|φ|⁴, and from CC12 eq. 4. The explicit b formula sits in CCM07 eq. 3.15, which I did not transcribe.)

**Spectral action, Λ⁰ (f0) terms** (Thm 3.13 and §3.7.3-3.7.6), with f0 = f(0):
- Higgs kinetic: f0/(8π²) Tr(M(D_μφ)²) = (f0/2π²) · a · |D_μφ|²  (eq. 3.49)
- Higgs quartic: (f0/2π²) · b · |φ|⁴  (eq. 3.51)
- Yang-Mills: (f0/2π²)(g3² G² + g2² F² + (5/3) g1² B²). The factors come from the fermion-content traces 24, 24 and 80 over H_F (§3.7.6).
- Mass term: (1/π²)(−2a f2 Λ² + e f0)|φ|²  (eq. 3.50). This is the only Higgs term carrying a different moment, f2 = ∫₀^∞ f(u)u du.

**Normalization ⇒ constraints at the cutoff Λ.** Normalizing the gauge kinetic terms gives (eq. 4.10)

    g3² f0 / (2π²) = 1/4,     g3² = g2² = (5/3) g1².

Normalizing the Higgs kinetic term needs H = (√(a f0)/π) φ (eq. 4.9). The physical Yukawa matrices are then k_x = (π/√(a f0)) Y_x (eq. 4.13). That gives the constraint (eq. 4.17)

    Tr( k*_(↑1)k_(↑1) + k*_(↓1)k_(↓1) + 3(k*_(↑3)k_(↑3) + k*_(↓3)k_(↓3)) ) = π²/f0 = 2 g²,

which in mass and SM-Yukawa form reads (eqs. 4.16, 5.25-5.26)

    Σ_σ [ (m_ν^σ)² + (m_e^σ)² + 3(m_u^σ)² + 3(m_d^σ)² ] = 8 M_W²,
    **Y₂(S) ≡ Σ_σ [ (y_ν^σ)² + (y_e^σ)² + 3(y_u^σ)² + 3(y_d^σ)² ] = 4 g²   at Λ.**

The quartic relation (eq. 5.6 and Table 1) is

    λ̃(Λ) = g3² · b / a²     (equivalently λ0 = g² b/a²),

and ξ0 = 1/12 (the R|H|² coupling).

CCM07's own remark after Lemma 4.6 is load-bearing for Part (b): "g² appeared in the same way on both sides and drops out." Both relations are ratios of f0 terms.

## 2. Predictions and their status

**Top mass (CCM07 §5.4).** Keeping only the top in Y₂ gives y_t(Λ) = (2/√3) g ≈ 0.597, with g(10¹⁷ GeV) ≈ 0.517 (eqs. 5.27-5.28). They run this down with the **one-loop** equation dy_t/dt = (1/16π²)[ (9/2)y_t³ − ((17/12)g1² + (9/4)g2² + 8g3²) y_t ] (eq. 5.9). That gives y_t(M_Z) ≈ 1.102 and m_t ≈ 173.683·y0 GeV, about 191 GeV. Their see-saw argument sets y_ντ = y_t at Λ, so Y₂ ≈ 4y_t², which gives y_t ≈ g (eq. 5.29). That lowers y_t(M_Z) to about 1.04, around 181 GeV, which CCM07 calls "acceptable" given the neglected couplings. (The GeV products are my arithmetic on their printed numbers.)

**Higgs mass (CCM07 §5.2).** With top dominance, λ̃(Λ) ≈ (4/3)πα3(Λ) (eq. 5.10), λ0 = 4λ̃ ≈ 0.356 at Λ = 10¹⁷ GeV. One-loop RG running (eqs. 5.13-5.14) gives λ(M_Z) ≈ 0.241 and m_H ≈ 170 GeV (eq. 5.15). Including ν_τ changes 4/3 to 1 and gives about 168 GeV (Remark 5.1). The predicted range was 160-180 GeV.

**Status: falsified by m_H ≈ 125 GeV; repaired in CC12.** CC12 says plainly that the 160-180 GeV prediction is contradicted. They add that a 125 GeV Higgs destabilizes the SM potential below the unification scale, which breaks the "big desert" assumption. Their fix is the real singlet σ coming from the Majorana entry k_νR, "already present ... wrongly neglected". They take top and neutrino dominance with k_ν = √n k_u. With f0/(2π²) = 1/(4g²), the boundary conditions at Λ are (CC12 eqs. 13-15, 22)

    λ_h = (n²+3)/(n+3)² · 4g²,   λ_hσ = 2n/(n+3) · 4g²,   λ_σ = 2·4g² = 8g²,   k_t = √(4/(n+3)) g.

The Yukawa normalization is the same structure as CCM07, with Y₂ = 4g² now split between top and ν by n. One-loop RG over unification scales of 6.5×10¹²-1.4×10¹⁷ GeV gives m_H ≈ 125.5 GeV along a curve n(u), so one free parameter is traded for the Higgs mass. CC12 states the top mass comes out "a few percents lower than the experimental value" and hopes two-loop running fixes it. So the 2012 result is a consistency, not a prediction: n is fitted.

**Pati-Salam.** CCvS13 drops the first-order condition and obtains an SU(2)_R × SU(2)_L × SU(4) model. Its top and Higgs masses are left as "the subject of future work" (§VIII, text read). CCvS15 (abstract) shows that the Pati-Salam gauge couplings unify at one loop for several scalar contents. The abstract gives no Higgs or top numbers. So no updated Yukawa-relation prediction exists in these two papers.

SECONDARY (not verified here, flagged only): later analyses with the full singlet sector and higher loops are widely reported to move the top prediction. I have not read them, so nothing above rests on them.

## 3. Outer-factor dependence: verdict

**In the continuum: the outer factor cancels exactly, on any compact Riemannian spin 4-manifold. The coefficients 4 (in Y₂ = 4g²), 1 (in λ0 = g²b/a²), and 5/3 (in g3² = g2² = (5/3)g1²) come only from the finite factor.**

Reasoning, from the CCM07 computation above:
1. The three terms that enter are the Yang-Mills term, the Higgs kinetic term and the quartic term. All three are Λ⁰ terms of the same Seeley-DeWitt coefficient a₄ and carry the same moment f0. Each is a *local density* multiplied by the same √g d⁴x. The outer geometry enters a₄ only through curvature invariants (R, C², R*R*, and ξ0 R|H|² with ξ0 = 1/12 fixed). None of these multiplies the F², |Dφ|² or |φ|⁴ densities with a manifold-dependent coefficient. Normalizing fields removes f0, and the relation g3² f0/2π² = 1/4 then removes it from the Yukawa constraint (CCM07's "g² drops out").
2. The spinor dimension (4) multiplies all three terms and cancels.
3. What survives is (i) the inner traces a, b over H_F with colour weights, (ii) the gauge-generator traces 24, 24, 80 over the fermion content, and (iii) the fixed ratio between the Ω²/12 and E²/2 terms of a₄. Item (iii) is a universal Gilkey coefficient. **My inference, not stated by CCM07:** (iii) is universal in any dimension and Clifford traces scale with dim S. So the same numerical coefficients would appear on a d-dimensional outer factor at the same heat-kernel order, though for d ≠ 4 the couplings become dimensionful and the shared moment becomes f_{(d−4)/2}Λ^{d−4}.
4. The only Higgs-sector quantity that depends on the outer moments is the mass term μ0² = 2f2Λ²/f0 − e/a. It involves f2/f0 and is free in CCM07 as well.

So the Yukawa relation is a property of the inner factor, *provided* (A) the full operator is a product with a factorized square and (B) the asymptotic a₄ expansion is used. CC08 (abstract) reports that on round S³ the asymptotic expansion reproduces the full spectral action "with remarkable accuracy". On S³×S¹ the two-term sum is exact up to an "astronomically small correction", and a₄ and a₆ cancel. That supports (B) on a *smooth round* S³. It says nothing about a truncation.

**On GeoVac's truncated S³: the cancellation is not automatic, and in the current construction the relation cannot even be posed.** There are three separate reasons:

1. **The kinetic terms vanish on the diagonal algebra.** Paper 32 Prop. (bounded commutators) and eq. (omega_decomp) show [D_GV, b^GV] = 0 for the diagonal A_GV. So the gauge one-form vanishes, and the Higgs field is constant on the outer factor. Neither F² nor |Dφ|² exists, and the normalizations behind Y₂ = 4g² have nothing to act on. A nonzero gauge piece appears only for the operator-system multipliers (FullDiracTruncatedOperatorSystem), which is the representative `geovac/almost_commutative.py` uses. The relation can be formulated only on that representative.
2. **The module's Dirac does not factorize.** `almost_commutative.py` builds D = D_GV⊗1 + sign(D_GV)⊗D_F. Paper 32 (Door-4 seam paragraph) records the cross term in D² as 2|D_GV|⊗D_F ≠ 0, with norms 4.3, 13.2, 28.0 at n_max = 1, 2, 3. The CCM split "outer moment × inner trace" holds only for a graded product, where D² = D_GV²⊗1 + 1⊗D_F² (e.g. D_GV⊗γ_F + 1⊗D_F, which the odd outer factor requires per Dabrowski-Dossena, or the chirality-doubled γ_GV = σ_x⊗1). With the ungraded D, odd powers of D_F enter the spectral action weighted by |D_GV| moments. Then even the constant-Higgs potential's coefficients depend on the outer factor in a way CCM's do not.
3. **No local density at finite cutoff.** Even with a graded product, the truncated spectral action is a finite trace and not an a₄ integral. The gauge and Higgs kinetic coefficients become two *different* finite-trace functionals. One is quadratic in the operator-system "curvature" of a_i[D,b_i]; the other is quadratic in [D_GV, φ] for a multiplier φ. Nothing forces their finite-n_max ratio to equal the continuum ratio. If the truncated triple converges (Paper 38 proves state-space GH convergence for the commutative factor; convergence of the AC extension is OPEN per WH1), both should approach their continuum values. Their ratio would then approach the continuum one and give Y₂ = 4g² back with an n_max-dependent correction. That correction would shift the *coefficient* at finite cutoff while leaving the *inner traces* a, b, 24, 80 untouched.

**Verdict.** The constraint's coefficients belong to the inner factor in the continuum, and the outer factor cancels as a ratio of same-order a₄ densities. At finite n_max, GeoVac's outer factor *could* move the coefficient through a cutoff-dependent ratio of two finite-trace kinetic functionals. This holds only after GeoVac switches to (i) a graded product Dirac and (ii) the operator-system multipliers. In the current module the relation is undefined (diagonal algebra) or unfactorized (sign(D_GV) product).

## 4. What H1 did

H1 did **not** impose or test the CCM Yukawa relation. From Paper 32 §sec:higgs_h1 and Theorem thm:h1_yukawa_non_selection:
- It tested kinematics only: Hermiticity, the J sign pair (−,+), the sign triple (−,+,+) on the chirality-doubled model, order-zero on the finite factor, the order-one/chirality moduli count (rank 16, dimension 32), and whether the Higgs piece of the inner fluctuation is nonzero iff Y ≠ 0.
- Its spectral-action content is two checks, Tr(D_F²) > 0 and Tr(D_F⁴) > 0. Paper 32 itself labels these automatic for any nonzero Hermitian D_F. The module's `heat_kernel_trace`, `spectral_action_eigenvalues` and `yukawa_dependence_scan` compute only Tr e^{−tD²}, eigenvalues and Tr(D²) against y. No kinetic normalization, no f0 moment, and no relation between a, b and g appear anywhere.
- "Non-selection" means the *kinematic* axioms leave the Yukawa free. CCM agree: the axioms leave Y free, and the *dynamical* spectral-action normalization supplies one constraint (Y₂ = 4g²) plus λ0 = g²b/a². That constraint is one equation on a weighted trace, not a selection of individual Yukawa values. So H1's theorem is compatible with CCM: imposing the CCM relation would add one trace constraint and would still not select the values. Paper 32 should not be read as having tested this constraint.

## 5. Proposed minimal pre-registrable GeoVac test

**Question:** at finite n_max, does the GeoVac outer factor change the coefficient c in Y₂ = c·g² (continuum c = 4)? Equivalently, does it change the ratio of Higgs to gauge kinetic normalizations?

**Construction (no new physics, all existing pieces):**
1. Outer: the CH full-Dirac D_GV at n_max = 2, 3, 4 (and 5 if feasible), with the operator-system multipliers of `FullDiracTruncatedOperatorSystem`.
2. Inner: minimal C⊕H finite triple (`FiniteElectroweakTriple`) with one nonzero Yukawa y. Inner traces a = y², b = y⁴ in the 1-generation, lepton-only convention, fixed exactly.
3. Product: the **graded** D = D_GV⊗γ_F + 1⊗D_F, chosen so D² factorizes. Check this numerically first, as a sanity gate and not a scored quantity.
4. Two one-parameter fluctuations built from the *same* lowest nontrivial multiplier pair (a, b) of O_{n_max}:
   - gauge: A_ε = ε·a[D_GV,b]⊗T (T a U(1) or SU(2) generator of A_F), symmetrized with J;
   - Higgs: Φ_ε = ε·(multiplier φ)⊗[D_F, q] (the off-diagonal piece).
5. Spectral action S(ε) = Tr exp(−D_ε²/Λ²) at several Λ. Read off the ε² coefficients K_A(n_max, Λ) and K_Φ(n_max, Λ) by exact second-order perturbation, e.g. Duhamel or an exact finite difference in rational ε, not a fit.
6. Compute the same two coefficients on the continuum round S³, either analytically from a₄ with the same harmonic content or on a large n_max reference. Define δ_A = K_A/K_A^cont and δ_Φ = K_Φ/K_Φ^cont. The finite-cutoff coefficient is then c(n_max, Λ) = 4·δ_A/δ_Φ.

**Pre-registered outcomes:**
- GO (outer cancels at finite cutoff): δ_A/δ_Φ = 1 exactly at every n_max. This would be a finite-level theorem candidate.
- BORDERLINE (cancels in the limit): δ_A/δ_Φ ≠ 1 but → 1 monotonically as n_max → ∞ at fixed Λ/n_max, with the rate stated.
- STOP (outer factor moves the coefficient): δ_A/δ_Φ tends to a limit ≠ 1, or does not converge.
- Controls, measured before scoring per the standing rule: (a) the Y = 0 case gives K_Φ = 0; (b) δ_A and δ_Φ against the n_max-reference continuum at the largest n_max; (c) the ungraded sign(D_GV) product as a negative control. It should show an outer-dependent Tr(D_F) cross-moment in the constant-Higgs potential.

**Not scored here:** the finite-factor traces a, b, 24, 80. They are inner by construction. Any GeoVac "prediction" of top or Higgs mass would inherit the CCM07 falsification and the CC12 fitted n, and should not be attempted. Guardrails: WH1 says AC-extension convergence is OPEN, so a BORDERLINE or GO is about the commutative-factor cutoff, not about the AC extension. Paper 32's H1 non-selection theorem is unaffected by any outcome, since the CCM relation is one trace constraint and selects no values.
