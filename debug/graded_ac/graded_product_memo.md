# Graded almost-commutative product on the true CH window: instrument and controls

Date 2026-10-29. Prototype in `debug/graded_ac/`: `graded_product.py` (the instrument), `controls.py` (controls only), `controls.log` (output, about 110 s wall). The task came from `debug/lit_scan/ccm_spectral_action_yukawa_memo.md` §5. **The scored quantity was not formed.** Nothing evaluated the heat-coefficient functionals with the S³ outer factor. On S³ the controls measure only exact algebra and one-form norms.

## 1. Construction

**Outer factor.** The outer factor is the true Camporesi–Higuchi window: complete Dirac eigenspaces |λ| ≤ n_max + 1/2, realized on L²(SU(2)) ⊗ C² with D = 2σ·J_R + 3/2. This ports Route A of `tests/test_p38_true_ch_window.py`, with one change: the window eigenbasis is built by Clebsch–Gordan coupling of (J_R, spin), not by `eigh`. The basis is exactly orthonormal (to 2e-16) and diagonalizes D exactly. So D_GV = diag(±(n+3/2)) is exact over Q. Window dimensions are 4, 16 and 40, and dim O = 5, 30 and 91, matching the frozen test. The real structure J_GV = V^T U V comes out as an exact signed permutation and is held in Q(i) arithmetic (`SparseQ`).

**Finite factor.** The finite factor is CCM doubled, as in `geovac/almost_commutative.py`: H_F = C⁸, D_F = diag(M, M̄), γ_F = diag(1,1,−1,−1,−1,−1,1,1), J_F = σ_x⊗1·K, A_F = C⊕H acting on matter. It is re-implemented exactly (`FiniteExact`). The controls use rational toy Yukawas y_ν = 1/2 and y_e = 1, plus a complex check y_e = 1/3 + i/4. These are not a physical assignment.

**Products.** Three Dirac operators are built, exactly:

| name | D | D² factorizes? | one-forms split? | J signs |
|---|---|---|---|---|
| `dirac_dd` (recommended) | D_GV⊗γ_F + 1⊗D_F | yes, exact zero | yes | J_GV⊗J_Fγ_F: (+,−) = KO-1 |
| `dirac_swap` | D_GV⊗1 + γ_s⊗D_F | yes, exact zero | **no** | (−,+,−) with grading γ_s⊗γ_F |
| `dirac_sign` (production) | D_GV⊗1 + sign(D_GV)⊗D_F | **no**: cross term 2\|D_GV\|⊗D_F | yes | (−,+) |

**One-forms.** One GV multiplier pair (a, b) supplies both fluctuations:
- gauge: (a⊗T)[D, b⊗1] = a[D_GV,b] ⊗ Tγ_F, with T = π_F(1,1) = P_matter;
- Higgs: (a⊗1)[D, b⊗β] − (a⊗β)[D, b⊗1] = (ab) ⊗ [D_F, β], with β = π_F(0,1) = P_L^matter.

Both identities hold identically for the DD product because β is even. The code takes the commutators of the full D and checks them against the closed forms (residuals 0 to 1e-16). Self-adjoint, J-symmetrized fluctuations D + ω + ε′J′ωJ′⁻¹ with ε′ = −1 preserve J′D = −DJ′ to 2e-16.

**Machinery.**
- `eps2_heat_coefficient(D, F, t)`: the exact second-order Duhamel coefficient of Tr e^{−t(D+εF)²}.
- `kinetic_duhamel(D, C, t)`: (t²/2)∫₀¹ Tr(C e^{−stD²} C† e^{−(1−s)tD²}) ds.
- `eps2_trace_power(D, F, n)`: the ε² coefficient of Tr(D+εF)ⁿ. It is exact when the inputs are `SparseQ`.
- `exact_trace_D2k` and `factorized_trace_D2k`: exact Tr(D^{2k}).

The S³ multipliers carry Clebsch–Gordan square roots, so fluctuated S³ quantities are float64. The torus toy is exact throughout.

## 2. What was blocked, stated plainly

**Blocker 1: the γ_GV⊗D_F form does not exist on the true window.** A grading of the right kind would have to commute with the operator system and anticommute with D_GV. None exists. I re-measured this here: the commutant of O^CH has dimension 1 and the grading space has dimension 0 at n_max = 1 and 2 (Paper 32 prop:odd(c); n_max = 3 rests on the frozen test). This is the oddness of the truncation.

Dropping the requirement to commute with the algebra changes the picture. An involution anticommuting with D_GV then does exist. It is γ_s, the left/right (group-inversion) swap |+, n, m, M⟩ ↔ |−, n, M, m⟩ between the ±(n+3/2) eigenspaces, whose multiplicities are both (n+1)(n+2). γ_s² = 1, {γ_s, D_GV} = 0 and J_GVγ_s = +γ_sJ_GV, all exactly. So `dirac_swap` has an exactly factorized D². **What fails:** γ_s fails to commute with every nonscalar element of O (4/5, 29/30 and 90/91 basis elements at n_max = 1, 2, 3; max ‖[γ_s,b]‖/‖b‖ = 2). As a result the one-form does not split. The gauge one-form becomes a[D_GV,b]⊗T + a[γ_s,b]⊗TD_F. The extra piece carries D_F and has norm 1.86 at n_max = 2 and 2.50 at n_max = 3, against 4.80 and 8.70 for the clean piece. A "gauge" fluctuation that contains the Yukawa operator cannot be normalized separately from the Higgs. The swap product is therefore unusable for the relation.

**The honest alternative is the Dabrowski–Dossena product D_GV⊗γ_F + 1⊗D_F.** This is the product their §4.2 prescribes when the first factor is odd. The grading is taken from the even (finite) factor, so the oddness of the window is no obstruction. Its D² factorizes exactly and its one-forms split exactly. Its cost is in the real structure. With Connes' J_F the sign ε′ is **inconsistent** (exact test: neither +1 nor −1), reproducing DD Table 5, row 3₊, column 6₊. With J_F γ_F the result is KO-1, (ε, ε′) = (+, −), exactly at n_max = 1, 2 and 3. The product is odd: 1⊗γ_F does not anticommute with D. Order zero and order one still fail on the outer factor (Paper 32; not re-measured).

**Blocker 2: the diagonal algebra is resolved by the true-window O.** The fluctuations use compressed functions P f P from O^CH, not the diagonal A_GV. Gauge and Higgs pieces are nonzero (control c).

## 3. Control results (from `controls.log`)

**(a) Cross term, exact.**

| n_max | DD | swap | ungraded: nnz, ‖·‖_F (toy Y) | ungraded = 2\|D_GV\|⊗D_F exactly |
|---|---|---|---|---|
| 1 | 0 | 0 | 32, 13.416 | yes |
| 2 | 0 | 0 | 128, 40.988 | yes |
| 3 | 0 | 0 | 320, 86.948 | yes |

For DD, Tr D^{2k} equals Σ_i C(k,i) Tr D_GV^{2i} Tr D_F^{2(k−i)} exactly for k = 1, 2, 3. At n_max = 3, for example, Tr D⁴ = 36674 and Tr D⁶ = 886177/2. The ungraded product departs from this from k = 2 on: at n_max = 3, Tr D⁴ = 44234 and Tr D⁶ = 1415593/2. Tr D² agrees for all three products, because the cross term is traceless. The complex Yukawa gives the same exact zero for DD.

**(b) Signs, exact.**

| triple | (ε, ε′, ε″) |
|---|---|
| J_F | (+, +, −), KO-6 |
| J_Fγ_F | (−, −, −) |
| J_GV on the window | (−, +) at n_max = 1–3 |
| DD with J_GV⊗J_F | (−, *undefined*) |
| DD with J_GV⊗J_Fγ_F | (+, −) |
| swap with grading γ_s⊗γ_F | (−, +, −) |
| sign product | (−, +), matching Paper 32 |

**(c) One-form pieces, DD, toy Y.** The pair is a = P D^{1/2}_{−½,+½} P and b = P D^{1/2}_{+½,+½} P.

| n_max | ‖gauge‖ | ‖Higgs‖ | Higgs at Y = 0 | gauge at b = 1 |
|---|---|---|---|---|
| 2 | 4.797118 | 1.590869 | 0 | 0 |
| 3 | 8.702251 | 2.964635 | 0 | 0 |

**(d) Toy with a known answer.** The toy is the flat 2-torus [0,2π]² with the antiperiodic spin structure. Its window is |k| ≤ R, a union of complete Dirac eigenspaces, and its multipliers are compressed plane waves e^{iq·x} and e^{ip·x}. I chose it for two reasons. It has the same structure as the S³ case: an outer factor of complete eigenspaces and compressed function multipliers. And its continuum answer can be computed by a route that shares no code with the matrix machinery: a symbolic Gaussian k-integral (sympy) followed by an mpmath s-integral. By Poisson summation that route equals the infinite-lattice sum up to O(e^{−π²/t}). Two answers are known in closed form.
- **Pure gauge (q = 0) has continuum k_A = 0.** D + iε[D,χ] is unitarily equivalent to D, so k_A vanishes. The independent route gives 0, 0 and −2.8e-16 for three p values. This also validates the continuum route.
- **c must return 4 when finite equals continuum.** For p = (1,0), q = (0,1) and t = 0.1, the continuum values are k_A^c = 0.402477639612 and k_Φ^c = 0.607787540617. At fixed t, c(R²) runs −87.9, −56.6, −24.5, 1.66, 3.97974, 3.99999985, then **4.000000000000** at R² = 400. Over the same steps the pure-gauge k_A(R) falls from −8.2 to −8.5e-14. The machinery reproduces the continuum 4.

Machinery cross-checks: Duhamel against a Richardson finite difference agrees to 1.6e-8 (finite-difference-limited). The exact ε² moments of Tr(D+εA)ⁿ are 88, 1392 and 16278 for n = 2, 4, 6, equal to the float values. Through the full DD product with toy D_F, the gauge coefficient factorizes exactly: K_A^full = Tr_F(T²e^{−tD_F²})·k_A^out, to 3.6e-15. The Higgs kin-kin coefficient through the full product does **not** factorize: 1.03369 against 1.12488. It factorizes exactly once D_F is removed from the heat weights.

**Non-vacuity, and a warning for the scoring regime.** With the cutoff tied to the window (fixed τ = tR² = 4), the toy's c is −3.48, −10.66, −15.83 and −40.01 at R² = 16, 36, 64 and 144. The fixed-t run shows the source: the compression breaks gauge invariance, and the pure-gauge k_A is large and negative whenever the heat cutoff reaches the window edge. This toy result is not a statement about S³. It does show that the machinery can return c ≠ 4, and that a fixed-Λ/n_max score is exposed to gauge-breaking contamination.

## 4. The quantity a pre-registration should score

Fix the following before any S³ evaluation:

1. Outer: `true_window(n_max)`, D = D_GV (exact diagonal), n_max ∈ {2, 3, 4, 5}.
2. Pair: a = P D^{1/2}_{−½,+½} P and b = P D^{1/2}_{+½,+½} P (`TrueWindow.multiplier`). These are the controls' pair, and da∧db ≠ 0 in the continuum.
3. Gauge field A = ω + ω†, with ω = a[D_GV, b]. Higgs outer field φ = a·b, the *compressed product* that the DD Higgs one-form produces. C = [D_GV, φ].
4. Outer functionals:
   - k_A(n_max, t) = `eps2_heat_coefficient(D_GV, A, t)`;
   - k_Φ(n_max, t) = `kinetic_duhamel(D_GV, C, t)`.

   These are outer-only by design. The gauge one is exactly the full-product coefficient divided by Tr_F(T²e^{−tD_F²}), shown in (d). The Higgs one is the leading-a₄ kin-kin term and is D_F-free. The full-product kin-kin carries D_F weights at O(tY²). **No D_F enters c.** The 4 is supplied by CCM's inner traces, as in the memo.
5. Continuum reference: k_X^∞(t) = lim over n_ref of k_X(n_ref, t) at fixed t, with the same construction. Stop at the first n_ref with a relative step below 1e-10. On S³ the Toeplitz residual is a boundary term (`tests/test_true_window_real_structure.py`), so edge suppression e^{−t(n_ref+3/2)²} applies as on the toy. **This convergence must be measured and logged first, as a control, before scoring.**
6. Cutoff: t = τ/(n_max + 1/2)², with τ ∈ {1, 4} fixed in advance.
7. **Scored:** c(n_max; τ) = 4·[k_A(n_max,t)/k_A^∞(t)] / [k_Φ(n_max,t)/k_Φ^∞(t)].

Measure these controls before scoring:
- (i) the same pipeline at fixed t = 0.1 with growing n_max must approach 4, as on the toy;
- (ii) the pure-gauge companion k_A^{pg}, taken with a = 1, must have k_A^{pg,∞} = 0 to the reference tolerance. Report ρ_pg(n_max) = k_A^{pg}/k_A at the scored t;
- (iii) Y = 0 gives a zero Higgs one-form;
- (iv) the ungraded product's exact cross term is nonzero.

**Gate proposal.**
- **GO** (the outer factor cancels at finite cutoff): |c/4 − 1| < 1e-10 at every n_max and both τ.
- **BORDERLINE**: c ≠ 4, but |c/4 − 1| decreases monotonically over n_max = 2–5 at both τ, with the rate stated.
- **STOP** (the outer factor moves the coefficient): c is non-monotone, changes sign, grows, or settles away from 4 at either τ.
- **VOID**: control (i) or (ii) fails.

Prior from the toy, to be stated in the registration and not taken as a prediction for S³: at fixed τ the toy is STOP-type. Its c diverges negatively because the compressed gauge functional is not gauge invariant. A gauge-breaking-subtracted variant must be defined *before* scoring, or not at all.

Scope (WH1): a verdict concerns the commutative-factor cutoff on the DD product only, not AC-extension convergence or Paper 32's H1 theorem.

## 5. Continuum reference (added 2026-10-29, PI decision)

The reference k_X^inf(t) is the smooth-sphere quantity, not a large window (`continuum_ref.py`, log `continuum_ref.log`, 36 s). (R1): uncompressed matrix elements of a, b and the true product ab between exact CG-coupled Dirac eigenstates, all levels n <= 51, plus a rigorous Gaussian tail bound (algebraic norms <= 6, level jump <= 3). Tails are <= 6.5e-30 relative at the eight scored t; 9.9e-16 at t = 0.02. (R2): Gilkey a_4 on unit S^3, int|F|^2 = 8 vol, int|dphi|^2 = (4/3) vol, so k_A ~ (2/3) sqrt(pi) t^(1/2), k_Phi ~ (1/3) sqrt(pi) t^(1/2), ratio 2 (constants checked on the flat torus to 4e-4 at t = 0.001). R(t) = 2.352, 2.168, 2.082, 2.040 at t = 0.16..0.02; (R-2)/t -> about 2. Pure gauge k^inf <= 3.6e-15. Window c at fixed t = 0.1: 3.942, 3.9987, 3.99998753 at n_max = 10, 12, 14. The fixed-tau toy prior and the draft pre-registration are in `toy_prior.log` and `prereg_draft.md`. Section 4 item 5 (large-n_ref reference) is superseded by this.
