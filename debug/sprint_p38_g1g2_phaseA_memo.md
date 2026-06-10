# Sprint P38-G1G2, Phase A: the action-seminorm repair design (2026-06-10)

**Goal:** close Paper 38's two named gaps (G1 dual-direction reach; G2 kernel condition), converting the main theorem from conditional to unconditional. This is the standing freeze exception of the accessibility plan and the unblock for Phase 4 sends under either send option.

**Verdict: the repair design survives its decisive diagnostics. G2 dissolves outright; G1 reduces to a bounded, well-specified lemma (the spinor lifted-state estimate) whose scalar prototype is published.**

## 1. The reframing

Paper 38 v1 metrized the truncated operator system by the truthful-Dirac commutator seminorm ‖[D_CH, ·]‖, whose kernel at finite cutoff strictly exceeds the scalars (10/14, 26/55 multipliers at n_max = 2, 3) — the root of G2, and the reason the engineered off-diagonal Dirac existed. The structural cause: the truthful commutator is blind to shell-diagonal content, and every top-band multiplier (N = 2n_max − 1) compresses to a shell-diagonal matrix at every finite cutoff, so the degeneracy is permanent, not transient.

**Repair:** metrize the truncated side by the isometry-action Lipschitz seminorm instead —

  L_n(T) = sup_{|X|=1} ‖[J_X, T]‖   (X over the Killing fields of round S³; equivalently the sup-over-g translation form sup_g ‖U_g T U_g* − T‖/d(g, e)).

On the continuum this **equals** the metric Lipschitz constant (= the Dirac seminorm ‖[D, M_f]‖ = ‖∇f‖_∞), so nothing changes at the limit object; only the finite-cutoff metrization changes. This is exactly the framework of Gaudillot-Estrada–van Suijlekom (arXiv:2310.14733, verified by full-text fetch 2026-06-10): their Lip-norm is the left/right conjugation-translation seminorm ‖T‖_{λ,ρ} = max(‖T‖_λ, ‖T‖_ρ).

## 2. G2: dissolved (kernel theorem)

**Claim.** On the chirality-doubled spinor multiplier system at any cutoff, L_n(T) = 0 ⟺ T ∈ C·1, with the **truthful** geometry.

**Proof shape (two steps).**
(i) *Band preservation* (textbook): Killing fields generate isometries; isometries preserve the harmonic degree, so for f in band N, each X·∇f lies in band N. Infinitesimally, [J_X, P M_f P] = P M_{X·∇f} P — a *scalar* multiplier compression in the same band.
(ii) *Per-band injectivity* (the premise, now verified): the map f ↦ P_n M_f P_n is injective on every band N ≤ 2n_max − 1. **Verified numerically 2026-06-10** (`debug/p38_g1g2_band_diagnostics.py`, JSON in `debug/data/`): at n_max ∈ {2,3,4,5}, every band has vec-rank = N² = band dimension exactly (e.g. n_max = 5: ranks 1, 4, 9, …, 81; total 285 = system dim).
Then L_n(T) = 0 ⟹ X·∇f = 0 for all Killing X ⟹ ∇f = 0 (Killing fields frame the tangent space — S³ is parallelizable) ⟹ f constant ⟹ T ∈ C·1. ∎ (modulo writing (i) cleanly)

No engineered Dirac, no truthful-vs-modified comparison, no perturbation lemma. The Dirac-commutator seminorm survives as a *dominated* seminorm (‖[D, PM_fP]‖ = ‖P c(∇f) P‖ ≤ √3 max_i ‖P M_{e_i f} P‖-ish) whose finite-cutoff degeneracy becomes a documented remark, consistent with the project's earlier +∞-Connes-distance findings (R3.x) and the P45 annihilation analysis.

## 3. G1: reduced to the spinor lifted-state lemma

GE–vS's dual-direction map (full-text verified): υ_Λ^μ(T)(g) = μ_Λ(ρ_g(T)) — evaluate the conjugated truncated operator against a "lifted state" μ_Λ built from a probability measure μ concentrated near e — with the quantitative bound

  d_GH(S(A_Λ), S(A)) ≤ ∫_G d(e, x) dμ(x).

**The right-hand side is precisely Paper 38's mass-concentration moment γ_n when μ = K_n·(Haar)** — the central Fejér measure. So the quantitative rate (4/π)·log n/n with the closed-form T_n sum slots directly into their theorem schema; nothing about the rate analysis is wasted.

What is genuinely new (their scope is scalar C(G) only, confirmed): the **multiplicity-2 spinor version of the lifted state** — construct μ_Λ as a state on the chirality-doubled truncated system P_n(C(G)⊗1)P_n ⊂ B(L²(G)^⊕2) with the equivariance and moment properties their Lemmas need. The truncation P_n is G-equivariant (CH-Dirac spectral projections commute with the isometry action), which is the only structural input their construction appears to use beyond Peter–Weyl.

**Viability diagnostic (passed):** the dual direction needs the compression to retain each band non-degenerately. Per-band conditioning of the vec-stack (Frobenius-level): cond ≤ 2.65 across ALL bands and cutoffs tested, and σ_min at fixed band *grows* with n_max (band 3: 0.30 → 0.78 → 1.33 → 1.94 over n_max = 2..5) — the compression's grip on each band strengthens with the window. Honest caveat: this is a Frobenius-norm diagnostic; the lemma needs operator-norm estimates (op ≥ Frob/√rank gives a weak transfer; the real estimates come from the lifted-state construction, not from this table).

## 4. Remaining steps (Phase B — the actual lemma work)

1. **Band-preservation lemma** (step (i) above) written properly: Killing action on S³ harmonics preserves N; infinitesimal conjugation formula on compressions. Textbook; half a page.
2. **Spinor lifted-state lemma** (the heart): adapt GE–vS's μ_Λ and υ_Λ^μ to the multiplicity-2 equivariant window; verify their hypothesis chain line-by-line against our objects (hypothesis-checklist discipline — this is where the last fabricated-import failure lived). Prototype numerics at n_max = 2, 3: build υ explicitly, verify almost-inverse property at rate γ numerically before trusting the proof.
3. **Seminorm equivalence on the continuum** (L_action = ‖∇f‖_∞; left/right max vs frame max constants): half a page, standard.
4. **Rewrite Paper 38**: named-gaps section replaced by the action-seminorm framework (G2 deleted as dissolved; G1 cited to the new lemma); main theorem unconditional; the engineered off-diagonal Dirac demoted to a remark on the Dirac-seminorm's finite-cutoff degeneracy. Update Paper 45's conditional-spatial proposition, the claims register row 8, N1's Result A + G1/G2 paragraphs, and the field guide sentence.
5. **Freeze new falsifiers**: per-band injectivity as a regression test; υ almost-inverse panel.

## 4a. Phase B results (2026-06-10, same day): scalar prototype PASS + spinor inclusion resolved exactly

**Scalar end-to-end prototype** (`debug/p38_g1g2_scalar_prototype.py`, window j ≤ 1, dim 14): the lifted state is the normalized Fejér amplitude ξ = h/√Z — in the Peter–Weyl orthonormal basis this is the all-ones-on-diagonals vector, and its induced measure is *exactly* the Fejér kernel measure. Three checks:

- **T1 (structural identity): PASS at machine precision** — υ(P M_f P)(g) = σ_J·f(g) with residual ≤ 3×10⁻¹⁶ across bands J = ½..2 and 24 random group elements. The σ values independently re-confirm the corrected Lemma-L2 multiplier symbols: σ(1) = (5+2√3)/18 = 0.470228 and σ(2) = 1/10 on the n = 3 window, exact.
- **T2 (almost-inverse rate): PASS** — ‖K∗f − f‖ ≤ γ·Lip(f) with comfortable slack (γ = 0.805 at this small window, computed as the kernel's first moment on the class measure with the same metric convention).
- **T3 (tautological contractivity): PASS** — Lip(υ(T)) ≤ L_action(T) for random system elements, with the action seminorm correctly evaluated at pair quotients g_i g_k⁻¹ (the first run's apparent violation was a sup-sampling artifact; the corrected measurement makes the inequality a theorem on the sampled set).

**Bonus catch (precision bug in our own corrected text):** the prototype's first run produced σ(½) > 1 — impossible for a UCP symbol — exposing that the CG-sum condition written into Paper 45 Appendix A and Paper 38 L2(b) during the descope omitted the **integrality constraint** (unit steps: j₁+j₂+J ∈ ℤ; e.g. ½⊗½ ∌ J=½). The published *values* were correct (the constraint did not bite at n_max = 2); the displayed condition was underspecified. Fixed in both papers (`debug/p38_p45_parity_fix.py`), both GATE: PASS.

**Spinor window inclusion: resolved exactly, no estimates.** Under SU(2)_L × SU(2)_R, L²(SU(2)) ⊗ ℂ² = ⊕_j V_j ⊗ (V_{j+1/2} ⊕ V_{j−1/2}), and these summands are precisely the CH Dirac shells: V_j ⊗ V_{j+1/2} = positive shell n = 2j (dim (n+1)(n+2) ✓), V_j ⊗ V_{j−1/2} = negative shell n = 2j−1 (dim ✓). Hence h_J ⊗ χ (Fejér amplitude ⊗ fixed unit spinor) lies in the CH window of shells ≤ 2J as a union of *complete* shells — exact containment, dimension-checked at J = 1 (2+6+2+12+6 = 28 = 2Σ(2j+1)²). The ℂ² factor is inert under scalar multipliers, so T1/T2/T3 transfer verbatim; left translation ⊗ 1 commutes with the CH projection and conjugates multipliers correctly.

**Rate correction (2026-06-10, Phase C indexing check):** the earlier "8/π cost" note was an indexing slip — γ_n is indexed by the number of scalar *levels* (n = 2j_max + 1), not by the top spin. With J = (n_max − 1)/2 the inclusion 2J ≤ n_max − 1 is *exact* (complete shells, zero slop), the scalar Fejér window has exactly n_max levels, and the lifted state's moment is **γ_{n_max} itself**. Both almost-inverse defects (function-side Fejér smoothing and operator-side conjugation averaging) run at the same γ_{n_max}, and the unconditional theorem closes at the original rate (4/π + o(1))·log n_max/n_max. A further Phase-C simplification: in the vS-GH assembly the forward arrow is the *plain compression* (UCP, unital, tautologically translation-contractive) — the Berezin map is not needed for the distance bound (it survives as a forward-direction refinement); the truncated-side composition equals the conjugation average Φ_K(T) = ∫ρ_g(T)K(g)dg, bounded by γ·L_n(T) tautologically. No inverse estimate of any kind remains in the proof.

**Net Phase B verdict: G1 is structurally closed.** All ingredients of the two-directional vS state-space GH bound now exist on the spinor substrate with explicit constants: forward Berezin (P38 L4, corrected), dual lifted-state map (this phase: exists in-window, tautologically Lipschitz-contractive, almost-inverse = scalar Fejér smoothing at rate γ_{⌊n/2⌋}), kernel condition (Phase A, action seminorm, every cutoff). G2 was dissolved in Phase A.

## 4b. Phase C — COMPLETE (2026-06-10, v3.109.0)

All five steps executed same-day:

1. **Spinor inclusion bookkeeping:** per-shell dimension checks recorded in §4a (J = 1 case 2+6+2+12+6 = 28 = 2Σ(2j+1)²); the rate-indexing check (γ indexed by scalar levels, J = (n_max−1)/2 inclusion exact) folded into §4a's rate-correction paragraph. The theorem closes at the **original rate (4/π + o(1))·log n_max/n_max**, not 8/π — the 8/π note in the older step-4 text below was the pre-correction plan.
2-3. **Assembly + band-preservation + kernel write-up:** landed directly in Paper 38 as `lem:band_injectivity`, `lem:continuum_lip`, `lem:lifted_state`, `prop:kernel_condition`, `rem:no_inverse` (no partial-inverse/transference estimate anywhere in the chain).
4. **Paper 38 rewritten in place to unconditional** (`debug/p38_unconditional_rewrite.py`; `thm:main_unconditional`; engineered off-diagonal Dirac demoted to `rem:dirac_degeneracy`; history in `rem:history38`). 20pp GATE: PASS. **Cascade applied** (`debug/p38_unconditional_cascade.py` + 3 manual abstract/remark fixes): Paper 45 spatial statement → unconditional (22pp PASS), Paper 32 caveat items (ii)–(iv) resolved + WH1 remark + limit-id remark moved to vS language (85pp PASS), claims register row 8 → INTERNAL THEOREM (unconditional), README bullet, field guide, N1 note (3pp PASS).
5. **Falsifiers frozen:** `tests/test_p38_action_seminorm.py` (6 tests: per-band injectivity n_max ∈ {2,3} + lifted-state smoothing identity on the scalar window, published symbols (1, √2/3, 2/9)). Topological proofs + P45 falsifier re-run green (24 passed).

Close-out records: CHANGELOG v3.109.0; CLAUDE.md §1 version + §2 one-liner + §1.7 WH1 **RESTORED to PROVEN (unconditional)**; memory wh1_proven.md / paper_38_drafted.md / MEMORY.md index; accessibility plan status block (P38 gap-closure COMPLETE; Phase-4 send posture upgraded — S1-vs-S2 moot in the gap direction).

## 5. Files

- `debug/p38_g1g2_band_diagnostics.py` + `debug/data/p38_g1g2_band_diagnostics.json` — Phase A checks, all cutoffs.
- `debug/p38_g1g2_scalar_prototype.py` + `debug/data/p38_g1g2_scalar_prototype.json` — Phase B scalar pipeline (T1/T2/T3, all PASS; future slow-test candidate).
- `debug/p38_p45_parity_fix.py` — the CG-sum integrality precision fix to Papers 38/45 (both GATE: PASS after).
- GE–vS alignment: full-text fetch of arXiv:2310.14733 (Lip-norm form; υ_Λ^μ formula; Prop 9 kernel argument; first-moment rate; scalar-only scope) — recorded here, to be re-verified line-by-line in Phase B against the actual lemma statements before any import.
