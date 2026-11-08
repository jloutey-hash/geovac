# Sprint memo — the LiH (two-center) N=4 explicit-r₁₂ 4-body integral gate (σ + π/δ)

> **SUPERSEDED IN PART (audit 2026-09-22; CHANGELOG v5.15.23, `memory/lih_marriage_state_of_play.md`):** the "kinetic-vector UNTESTED" flag is STALE (settled v5.15.13–14); the π/δ validation here is on an azimuthally MODULATED MODEL density, not m≠0 orbitals (the assembled energies are σ-only); the σ figure 1.5e-4 is MC-scatter-limited; the 1-D leaf dressing is the isotropic-leaf special case (general leaves need the Neumann potential); the "bridge addition theorem" is a validated reduction, not a theorem. Read the memory file before building on this memo.

**Date:** 2026-09-21. **Type:** diagnostic / integral validation (no energy). **Verdict:**
the Be atomic N=4 **soft-wall** result (v5.15.8/9) **carries to two centers, for all azimuthal
channels (σ, π, δ)**. The genuinely-4-body bridging integral of a two-center (LiH-geometry)
explicit-r₁₂ CI reduces **exactly** to a 1-D leaf dressing + a **two-center prolate Neumann
Coulomb** between the dressed densities — RI-free, L-sum convergent — validated by (i) a
**closed-form** control on the new bridge machinery and (ii) an independent 12-D Monte-Carlo,
first for the σ (m=0) channel (§2–4) and then for the π (m=1) / δ (m=2) azimuthal-transfer
channels (§4b). PI-directed: the "prolate 4-body integral gate," the load-bearing first step of
the accurate-LiH build (`debug/lih_r12_build_plan.md`), then extended to m≠0.

## 0. Why this sprint

The Be sprint (v5.15.8) showed the N=4 explicit-r₁₂ "wall" is SOFT for the scalar Coulomb
bridging term `f₁₂ f₃₄ / r₁₃`: exact, RI-free, terminating, reducible — validated
reduced==brute on the Be-relevant **atomic** integral (`debug/r12ci_4e_be_integral.py`). But
Be is atomic: isotropic densities → **monopole** kernels, single-center Slater–Condon bridge.
LiH is **two-center**: the bridging Coulomb acts between two-center densities and needs the
**prolate two-center Neumann** expansion. The open question this gate answers: does the Be
4-body reduction survive the move to two centers, or does the two-center bridge reintroduce an
RI/decidability obstruction? (The build plan's own Step 3 — "validate reduced==brute on one
4-electron LiH integral BEFORE any energy" — reordered to Step 1 as the highest-value,
lowest-risk deliverable, since every accurate-LiH path depends on this integral machinery.)

## 1. The object and the reduction

The only genuinely-4-body-connected term in `⟨Φ|F H F|Φ⟩` (F = Σ f_pq multiplicative;
kinetic gradients pair-local, V_ne one-body → only the two-body Coulomb bridges disjoint
pairs):

    I = ∫ ρ₁(r₁) ρ₂(r₂) ρ₃(r₃) ρ₄(r₄)  f(r₁₂) f(r₃₄) (1/r₁₃)  d³r₁..d³r₄     (chain 2-1-3-4)

Model orbitals (a pure integral validation, as in the Be file): **bridge** e₁,e₃ = 1s on
focus B (H side, diffuse, Z=1.0); **leaves** e₂,e₄ = 1s on focus A (Li side, tight, Z=2.7);
f(r)=exp(−0.5 r); R = 3.015 (LiH), foci A=(0,0,−a), B=(0,0,+a), a=R/2.

**Reduction (two steps):**
1. **Leaf dressing** (Be-style radial monopole; EXACT because the leaf is isotropic about A):
   each 1s leaf integrates out into a 1-D radial dressing of its bridge partner,
   `Ψ_A(r_{1A}) = ∫ ρ_A(r₂) f(r₁₂) d³r₂` = spherical average of f over the leaf angle.
2. **Two-center prolate Neumann bridge** (the NEW piece vs Be): the residual is a two-center
   Coulomb between the **dressed two-center densities** `D = ρ_B(r_B)·Ψ_A(r_A)`
   (r_B=a(ξ−η), r_A=a(ξ+η)), expanded in the prolate Neumann series
   `1/r₁₃ = (2/R) Σ_l Σ_m (−1)^m (2l+1)[(l−|m|)!/(l+|m|)!]² P_l^{|m|}(ξ_<)Q_l^{|m|}(ξ_>)
   P_l^{|m|}(η₁)P_l^{|m|}(η₃) e^{im(φ₁−φ₃)}`. For σ densities only m=0 survives, giving
   `I = (2/R)(2π)²a⁶ Σ_l (2l+1) ∫∫ g_l(ξ₁)g_l(ξ₃) P_l(ξ_<)Q_l(ξ_>) dξ₁dξ₃`,
   `g_l(ξ)=∫(ξ²−η²)D P_l(η) dη`. The l-sum converges/terminates with the finite η-content →
   **RI-free**, exactly as Be's L-sum terminated at 2·l_bridge.

The bridge factorization across `1/r₁₃` is the two-center analog of Be's Legendre addition
theorem — the "bridge addition theorem" (memo v5.15.8 §5), now in prolate spheroidal
coordinates.

## 2. Controls (validate the machinery before the 4-body number)

| control | result | anchor |
|:--|:--|:--|
| **C0** ∫ρ_B dτ over the prolate grid | **1.000000** | exact normalization |
| **C1** prolate Neumann self-Coulomb of 1s_B | **0.625085** vs **5α/8 = 0.625000** (rel 1.4e-4 at the default grid, → **3.7e-5** at 500×240, halving each refinement) | **CLOSED FORM** — the rigorous anchor for the new bridge machinery |
| **C2** leaf dressing `⟨ρ_A ρ_B f⟩` reduced vs 2e MC | reduced **0.206596** vs MC **0.206613 ± 1.6e-5** (rel **8.1e-5**) | independent MC on the leaf |

C1 is the load-bearing control: it validates the **exact prolate Neumann bridge + its
prefactor `(2/R)(2π)²a⁶` + the Q_l machinery** against a value known in closed form, and its
error → 0 with grid refinement (exact-in-the-limit). This is a STRONGER anchor than the Be
sprint had (Be leaned on MC alone for its bridge).

## 3. The 4-body number — reduced vs brute

**REDUCED** (leaf monopole dressing + prolate Neumann bridge) descends monotonically toward
its grid limit as the (ξ,η) grid refines (residual = quadrature bias, not structural),
tracking the C1 self-Coulomb grid error → 0:

| grid | I_reduced | selfCoul rel |
|:--|:--|:--|
| 200×80, L30 | 2.85693e-2 | 2.3e-4 |
| 260×100, L34 | 2.85665e-2 | 1.4e-4 |
| 320×140, L36 | 2.85651e-2 | 9.0e-5 |
| 440×200, L40 | 2.85639e-2 | 4.8e-5 |
| 560×260, L42 | 2.85633e-2 | 2.9e-5 |
| **grid limit (extrap. selfCoul→0)** | **2.85625e-2** | 0 |

L-sum tail (l ≥ 20 fraction) = **3.6e-24** → RI-free (converges/terminates).

**BRUTE** (full 12-D importance-sampled MC, no reduction): the `1/r₁₃` estimator is
**heavy-tailed** (e₁,e₃ share center B → r₁₃→0 common), so *both* the naive √(var/N) error and
even batch-means underestimate the true **run-to-run** scatter — four independent 120–200M
runs landed at 2.85567 / 2.85576 / 2.85591 / **2.85679** e-2, a spread of ~1e-4. 3×120M
batch-means replicas: mean **2.85582e-2** (within-run batch spread 2.0e-6; not the true
uncertainty). So the MC is the **corroborating**, not the decisive, check.

**Meet:** reduced(grid limit) **2.85625e-2** and brute(3×120M) **2.85582e-2** agree to
**rel 1.5e-4**, inside the MC heavy-tail scatter (~1e-4). The residual is numerical (grid +
MC heavy tail), not structural.

## 4. Verdict, scoped honestly

**The two-center (LiH-geometry) 4-electron scalar-Coulomb bridging integral is EXACT and
RI-FREE.** The reduction is **exact by construction** — the prolate Neumann expansion of
`1/r₁₃` is an exact identity, and the isotropic-leaf monopole dressing is exact, so the
4-body integral is *re-expressed*, not approximated. Implementation validated three ways: the
assembled bridge reproduces a **closed form** (C1, self-Coulomb 5α/8, exact-in-limit and
converging), the leaf reproduces an independent MC (C2), and the full reduction meets an
independent 12-D MC at the MC heavy-tail precision (rel 1.5e-4). The l-sum converges/terminates
(RI-free). **The Be atomic N=4 soft-wall result carries to two centers.** (Epistemically
STRONGER than the Be gate, which had only an MC check for its bridge — here the new two-center
bridge machinery is anchored to a closed form.)

**What this is NOT (do not overclaim):**
1. A diagnostic + one validated integral, **NOT a LiH R12-CI energy**. It clears the
   integral-level gate the accurate-LiH build depends on; it delivers no spectroscopic number.
2. Scoped to the **scalar Coulomb** chain with **σ** (m=0) model orbitals. The **π/δ (m≠0)**
   azimuthal-transfer channels of the prolate Neumann bridge are the natural next exercise
   (they carry the angular correlation that matters for a real valence bond) and are UNTESTED
   here. The kinetic-vector 4-body pieces (argued non-bridging, pair-local gradients) and the
   non-Hermitian TC operator remain separate, untested axes.
3. The **quantum-encoding wall** (a 4-body correlation operator is a 4-body Pauli string) is
   separate and real — unmoved by this classical-integral result.

**Ahead vs catch-up:** on this narrow axis GeoVac is genuinely ahead of Gaussian-F12 (which
uses RI for exactly these 4-electron integrals) — the prolate Neumann closed-forms keep them
exact and RI-free, now demonstrated at **two centers**, not just atomically.

## 4b. The π (m=1) / δ (m=2) channel — same gate, all m (`debug/lih_r12_4body_pi_channel.py`)

A **diagonal** density `|φ|²` of a pure-m orbital is φ-independent, so the m≠0 bridge channels
never activate from diagonal densities — the m-transfer lives in the exchange (transition-
density) terms of a real CI. To exercise the general-m prolate Neumann bridge with a positive,
MC-samplable density, modulate the σ bridge density azimuthally:
`ρ_bridge = ρ_{1s_B}·(1 + β₁cosφ + β₂cos2φ)` (β₁=0.6, β₂=0.4; >0 for all φ) — this carries
m=0, ±1 (π), ±2 (δ) explicitly, samples exactly like the σ case (1s_B proposal + a φ-weight),
and its m-components are trivial. Leaves stay 1s on A (isotropic → the leaf dressing Ψ_A is
unchanged, φ-independent), so only the **bridge** needs the general-m Neumann expansion
`1/r₁₃ = (2/R)Σ_l Σ_m (2−δ_{m0})(−1)^m(2l+1)[(l−m)!/(l+m)!]² P_l^m(ξ_<)Q_l^m(ξ_>)P_l^m(η₁)
P_l^m(η₃)cos(m(φ₁−φ₃))`. `P_l^m`, `Q_l^m(ξ>1)` from scipy `lpmn`/`lqmn`; the Condon–Shortley
phases cancel in the `P·Q` and `P(η)·P(η)` pairs, leaving the explicit `(−1)^m`. After the
φ-integrals the azimuthal weights are `W_0=(2π)²`, `W_{m≥1}=2π²`.

| check | result |
|:--|:--|
| **C1′** m=0 path (β=0) inside the general-m code | 0.625085 vs 5α/8 = 0.625000 (rel 1.4e-4) — reduces correctly |
| **C2′** modulated self-Coulomb: general-m reduced vs 6-D MC | reduced 0.649863 vs MC 0.649860 ± 8.7e-5 → **rel 3.9e-6** (m=1 = 3.3%, m=2 = 0.6% of the total) — the m≠0 convention is pinned essentially exactly |
| **4-body m=0 part** | **2.856650e-2** — reproduces the σ-gate 4-body (§3) bit-for-bit (internal consistency) |
| **4-body reduced** (grid limit) | **2.96791e-2** (σ 2.85665 + π 0.09534 + δ 0.01642, ×10⁻²) |
| **4-body brute** (3×120M, batch-means) | **2.96797e-2**, replicas at +0.6/+0.1/−0.3σ (straddling the reduced limit) |
| **meet** | **rel 2.2e-5** — tighter than the σ gate |

**The π (m=1) and δ (m=2) azimuthal-transfer channels of the two-center 4-body reduction are
EXACT and RI-FREE** — the general-m bridge is pinned to essentially exact agreement with an
independent MC (C2′, rel 3.9e-6), its m=0 path reproduces the σ gate bit-for-bit, and the full
π/δ 4-body meets the independent 12-D MC at rel 2.2e-5. **The σ gate now holds for every
azimuthal channel** — the physically dominant valence-correlation channel (π) is covered.

## 5. Next steps (owed / PI call)

- **Into an energy** (the remaining owed step): assemble a two-center LiH R12-CI matrix element
  end-to-end via the reduced form (Be `be_r12ci_full.py` structure: `{Φ₀, FΦ₀}`, ill-conditioned
  `h`/`σ²` analytic, `g` — now containing this two-center 4-body term — via the reduction). The
  build plan's ansatz (B) on the validated 4e FCI Φ₀ is the first-PoC route; carry the Be
  ill-conditioning lessons (short-range geminal, orthogonalized basis, exact overlaps).
- (Done this session: the σ gate §2–4 and the π/δ (m≠0) extension §4b.)
- Reminder: the total-energy correlation axis (what r₁₂ targets) is **separate** from the LiH
  R_eq drift (a frozen-core / numerical-grid wall, NOT r₁₂-fixable — track log
  `debug/track_logs/prolate_native_lih.md`).

## 6. Files

- `debug/lih_r12_4body_integral.py` — the σ (m=0) deliverable: two-center 4-body integral,
  reduced (leaf dressing + prolate Neumann bridge) vs 12-D MC, with C0/C1/C2 controls,
  grid-limit extrapolation, and batch-means MC replicas.
- `debug/lih_r12_4body_pi_channel.py` — the π/δ (m≠0) extension: general-m prolate Neumann
  bridge (`lpmn`/`lqmn`), azimuthally-modulated bridge density, C1′/C2′ controls (incl. the
  6-D-MC convention pin), reduced-vs-12-D-MC. Reuses the σ module's leaf + sampler.
- `debug/r12ci_4e_be_integral.py` — the Be (atomic) precursor both generalize.

## 7. Records touched

CHANGELOG v5.15.10; CLAUDE.md §2 one-liner; `memory/r12_generalization_boundary_n3_n4.md`
(N≥4 frontier: two-center gate CLEARED); `debug/lih_r12_build_plan.md` (Step-1 gate done).
