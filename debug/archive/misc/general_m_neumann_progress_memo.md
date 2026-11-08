# Quadrature-free / μ=2-stable general-m Neumann V_ee — build progress memo

**Goal.** Make the general-m (μ>0) prolate-spheroidal V_ee stable at μ=2 (currently
the `d^4 Q_l` instability blocks the δ-channels) and quadrature-free, extending
Paper 12's algebraic V_ee from σ-only to general m. This is the top-ranked
accuracy path (`memory/polyatomic_state_of_play.md` §5); the azimuthal work
(v5.11.18/19) reached H₂ 99.09% but only via quadrature/Gaussian, with μ=2
unstable.

**Target to replace.** `debug/prolate_ci_general_m.py` `vee_matrix()`: its radial
`Xtab` uses `legendre_deriv_poly(l,m)` (=d^m P_l) and `q_deriv(l,m,ξ)` (=d^m Q_l)
— **differentiation** — whose coefficients reach ~1e10 and overflow at large l
(its own comment, line ~365). The angular `Ytab` is already algebraic + exactly
selection-rule-terminated (l > Q+2s−m ⇒ 0). The instability is entirely radial.

## Findings (this session)

1. **Pivot = GO.** The associated-Legendre l-recurrence
   `(l−m+1)M_{l+1}=(2l+1)ξM_l−(l+m)M_{l−1}` never differentiates, so it avoids the
   overflow/cancellation. `debug/general_m_moments.py`, `general_m_moments2.py`.

2. **The bare Q-moment ∫ξ^p Q_l^m DIVERGES for m≥2** (Q_l^m ~ (ξ²−1)^{−m/2} at ξ=1).
   The physical moment carries the basis (ξ²−1)^μ (in the driver, the combined
   (ξ²−1)^s, s=(μ_i+μ_j+m)/2). **This weight cannot be expanded into ξ-powers** —
   the pieces diverge and only their sum converges (endpoint cancellation). So the
   moment tables must carry (ξ²−1)^s **intact** (unlike `neumann_vee.compute_Xl`,
   which expands ξ^p P_l into monomials).

3. **Regularized Q-moment B_l^{m=2,μ=1} = ∫ξ^p(ξ²−1)^μ e^{−αξ}Q_l^m: forward
   recurrence is mpmath-stable to 1.6e-15** through l=12 (`general_m_moments2.py`).
   Backward/Miller fails here (the regularized moment is not recessive).

4. **Reduced P via forward recurrence: solved.** Matches direct differentiation to
   1e-15 at all l≤24, m≤4 in **float64**, fixing the documented d^mP_l overflow.
   `general_m_reduced_recur.py`.

5. **Reduced Q *pointwise on the grid* recurrence is float64-UNSTABLE both
   directions** (Q is the recessive solution; error → 1.0 by l≈m+8). So a drop-in
   "evaluate Q_l^m by recurrence on the grid" does NOT work. The robust Q route is
   the **moment level** (finding 3, mpmath), not pointwise float64.

## The route (next increment)

Build a moment-level 2D-ordered radial engine, NOT a pointwise grid patch:

- `A_l^{m,s}(p) = ∫₁^∞ ξ^p (ξ²−1)^s e^{−αξ} P_l^m(ξ) dξ` and the B analog for Q_l^m,
  via **forward** l-recurrence (mpmath), seeded at l=m, m+1 by stable evaluation
  (the (ξ²−1)^s weight kept intact). s runs over (μ_i+μ_j+m)/2.
- The 2D ordered integral X_l^{m,s}(p1,p2) assembled from A/B at α and 2α by the
  same IBP split as `neumann_vee.compute_Xl` — but the inner partial integral
  S_l(ξ₀)=∫₁^{ξ₀}… must also keep the intact weight (no monomial expansion).
- Swap into `vee_matrix`'s `Xtab`; keep the algebraic `Ytab`.

**Validation gates.** (a) reproduce `debug/prolate_ci_general_m.py` V_ee at μ≤1 to
≥8 digits; (b) STABLE at μ=2 (no d^4 blow-up), reaching ~99.1% with δ-channels;
(c) exact reduction to `geovac/neumann_vee.py` at μ=0. Then `test_paper12_*`,
Paper 12 algebraic-V_ee claim σ-only → general-m, CHANGELOG.

**Open numerical question to settle first in the next increment:** whether the
low-l seeds A_m^{m,s}, B_m^{m,s} close in the Paper-18 {e^a E₁(a), ln, γ} ring
(truly quadrature-free) or are seeded by stable 1D mpmath quadrature (still kills
the differentiation instability, matching how μ=0 `neumann_vee` seeds B_l by
`scipy.quad`). The latter is sufficient for the μ=2-stability win; the former is
the quadrature-free ideal.

**Artifacts.** `debug/general_m_moments.py` (P forward-stable, Q bare divergent),
`general_m_moments2.py` (regularized Q-moment forward-stable 1.6e-15),
`general_m_reduced_recur.py` (reduced-P recurrence solves overflow; reduced-Q
pointwise unstable).
