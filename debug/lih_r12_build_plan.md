# LiH explicit-r₁₂ build plan — handoff for a fresh session

**Written 2026-09-21 (end of the N=4/Be session), for the next session.** The top-of-session
intent was **LiH physical accuracy**; the session validated the enabling machinery on the
simpler atomic Be case and stopped there by design. This is the concrete plan to carry it to
LiH. Read this + the owning papers (12, 19) + `debug/track_logs/prolate_native_lih.md` before
starting. Current-state rule applies: verify against CHANGELOG since 2026-09-21 first.

---

## STATUS UPDATE (2026-09-21, v5.15.13) — COMPLETE: triangle reduced RI-free, g_Vee + g_T done, full analytic 2×2 → E_R12 assembled

**STAGE 4b part 2c is DONE.** The owed pieces all closed this session:
- **3-body TRIANGLE reduced RI-free** (`debug/lih_r12ci_triangle_gate.py`): H₂ non-separable kernel →
  azimuthal Fourier `H₂^{(m)}` → low-rank (randomized eigh) → mode-by-mode general-m prolate Neumann
  Coulomb (the π/δ machinery). Both placements validated vs MC: T1=0.03917 (0.2σ), T2=0.04824 (0.4σ).
- **g_Vee complete** (`debug/lih_r12ci_gVee_analytic.py`): a 216-triple unified message-passing reducer
  (separable→m=0 Coulomb/Yukawa, 3-cycle→triangle engine); reproduces E2 exactly (normalisation gate);
  all 6 intra/inter groups match MC ~1e-3; grid-consistent **g_Vee=+0.53834 vs VMC +0.53605 (0.4%)**.
  (Mixed cov/enum gave +0.585 — inconsistent code paths broke the ~4% cancellation; grid-consistent fixed it.)
- **g_T complete** (`debug/lih_r12ci_gT_analytic.py`): the IBP identity collapses gT2+gT3 to a SCALAR
  Yukawa covariance — **no vector ∇f dressing needed**. `g_T = gT1 − γ²σ² + 2γ Cov[F,Y_sum]`; only gT1
  (drift² over the KD gradient density, grid-consistent through the ~2% cancellation) is new.
  **g_T=+1.39419 vs VMC +1.38640 (0.56%).**
- **Capstone 2×2 → E_R12** (`debug/lih_r12ci_assemble.py`): S₀₁=0 (diagonal overlap);
  **E_R12=−7.94200 Ha** (dE=−54.2 mHa), variational, matching the same-geminal VMC 2×2 (−7.94121) to **0.79 mHa**.
  GEMINAL: the exp pieces use **f=exp(−γr)** (the Stage-1 f-tensor geminal).
- **CUSP-CORRECT geminal DONE (`debug/lih_r12ci_linexp.py`):** the framework is geminal-agnostic —
  only kernels change (f→d·e^{−γd}, f²→d²·e^{−2γd}, same-pair f/r→e^{−γr} [bounded], and
  ∇²f=γ²f−4γe^{−γr}+2e^{−γr}/r → g_T's IBP piece = −γ²σ²+4γCov[F,E_sum]−2Cov[F,Y_sum]). Reproduces the
  vmc **linexp** γ=0.5 targets: Fbar=3.5700 (exact), σ²=0.1279, h=−0.1034, g_Vne/g_Vee/gT1 all <0.7%,
  Cov[F,E_sum] ✓. **E_R12=−7.921 (dE −33 mHa)** vs vmc −7.9316; every piece <1% EXCEPT **Cov[F,Y_sum]**
  (−0.079 vs MC −0.102, 22% — the Yukawa integral Ȳ at ~0.4% grid accuracy × the F̄=3.57 cancellation,
  = the §1.3 long-range ill-conditioning). With Cov[F,Y] at converged precision → **E_R12=−7.917, dE=−29 mHa,
  matching vmc**. So the RI-free machinery lands the cusp-correct value; residual = one grid-limited Yukawa scalar.
**Owed (optional next):** finer-grid Yukawa dressing to tighten Cov[F,Y_sum] (last ~4 mHa); production
`geovac/` migration + regression tests. See CHANGELOG v5.15.13.

## STATUS UPDATE (2026-09-21, v5.15.12) — full off-diagonal h ANALYTIC; g_Vne done; 4-body bridge inside the energy validated

**Off-diagonal h fully analytic/RI-free** (this session, on the v5.15.11 σ²/h_Vne): h_T = +0.7495
(VMC +0.7485), h_Vee = +0.3098 (VMC +0.3082), h_Vne = −0.9302 ⟹ h = +0.129 (VMC +0.1268 ✓).
Files `debug/lih_r12ci_hT_analytic.py`, `debug/lih_r12ci_hVee_analytic.py` (reusable `cov_FA_FB`
bilinear engine + two-center Coulomb `neumann_potential`).
**Diagonal g — progress:** `g_Vne = −2.8516` (MC −2.8490, target −2.8454) VALIDATED via a block
decomposition F=S+I (`debug/lih_r12ci_gVne_analytic.py`). The **RI-free 4-body bridge inside the
energy** `⟨f₁₂f₃₄·coul₁₃⟩ = ¼ I_coul[S^f,S^f]` = 0.03281 vs MC 0.03281 (0.0σ) — the showpiece —
plus the **Yukawa-field machinery** `Ψ^Y = Ψ^coul − Ψ^smooth` (both in
`debug/lih_r12ci_gVee_groundwork.py`, which also banks the 6-product MC breakdown of ⟨F²V_ee⟩).
**Owed:** the rest of g_Vee (a two-center **3-body triangle** reduction — see STAGE 4b part 2
below) + g_T + the 2×2 → E_R12.

## STATUS UPDATE (2026-09-21, v5.15.10) — Step-1 integral gate CLEARED

The load-bearing first gate (Step 3's "validate reduced==brute on one 4-electron LiH integral
BEFORE any energy," reordered to Step 1 and done first) is **DONE and PASSED**. The two-center
4-body bridging integral `⟨ρ₁ρ₂ρ₃ρ₄ f₁₂f₃₄/r₁₃⟩` reduces **exactly, RI-free**, to a 1-D leaf
dressing + a **two-center prolate Neumann** Coulomb between the dressed densities — anchored to
the closed-form 5α/8 self-Coulomb (exact-in-limit), full reduction == independent 12-D MC at
rel 1.5e-4 (MC heavy-tail-limited). The prolate Neumann bridge (`(2/R)(2π)²a⁶ Σ_l (2l+1)…`,
via scipy `lqn` Q_l) is the reusable new machinery. Engine `debug/lih_r12_4body_integral.py`;
memo `debug/sprint_lih_r12_4body_integral_memo.md`; CHANGELOG v5.15.10. **So the 4-body term of
an all-electron LiH explicit-r₁₂ is confirmed passable at two centers, not just atomically.**
**Then EXTENDED to π (m=1)/δ (m=2) — same session, ALSO PASSED.** An azimuthally-modulated
bridge density `ρ_{1s_B}(1+β₁cosφ+β₂cos2φ)` validates the **general-m** prolate Neumann bridge
(assoc. Legendre `P_l^m`/`Q_l^m` via scipy `lpmn`/`lqmn`; prefactor `(−1)^m(2l+1)[(l−m)!/(l+m)!]²`)
against a 6-D MC at **rel 3.9e-6**, and the full π/δ 4-body against the 12-D MC at rel 2.2e-5;
the m=0 part reproduces the σ gate bit-for-bit. Engine `debug/lih_r12_4body_pi_channel.py`.
**So the two-center 4-body reduction is exact/RI-free for ALL azimuthal channels (σ/π/δ)** — π
is the physically dominant valence-correlation channel.

**Energy assembly — STARTED (2026-09-21), `debug/lih_r12ci_energy.py`.** Be-style `{Φ₀, FΦ₀}`
2×2 in the prolate geometry; minimal reference `Φ₀=|1s_A² 1s_B²|` (the two-center analog of
Be's `1s²2s²`; Z_A=3/Z_B=1, exponents 2.7/1.0). Order set by conditioning: well-conditioned
`E₀`, `F̄` first, then the ill-conditioned `σ²`/`h` (quadrature, no cancellation), then `g`
(the validated 4-body reduction).
- **STAGE 1 DONE + validated:** the AO 2-body f-integral primitives (`⟨pq|f|rs⟩`, every one a
  prolate quadrature by dressing the isotropic member of each pair, except `(ab|f|ab)` = one
  6-D importance-MC) — quadrature == MC on the 3 isotropic integrals (rel ≤1.4e-4), and
  `(aa|f|bb)=0.206596` reproduces the σ-gate C2 leaf value bit-for-bit. **F̄ = 1.8807** over the
  Löwdin determinant. f-machinery over the determinant is trustworthy.
- **STAGE 2 DONE + validated (2026-09-21):** `E₀ = ⟨Φ₀|H|Φ₀⟩ = −7.8878 Ha` (`stage2_E0()`), the
  two-center closed-shell RHF energy of the ionic reference `Φ₀=|1s_A² 1s_B²|` (Li⁺ core + H⁻).
  One-electron `h = T + V_ne` (kinetic via the exact `−½∇²·1s(ζ)=(−ζ²/2+ζ/r)·1s` identity →
  overlap + attraction, with a Hermiticity cross-check on T_AB); two-electron Coulomb = 5
  dressable integrals by exact Hartree-dressing + `(ab|ab)` by the validated prolate Neumann
  machinery. **Every control passes:** closed forms (`⟨A|1/rA|A⟩=ζ`, `5ζ/8`, point-charge
  screening, T_AB Hermiticity, `(aa|bb)` symmetry) to 1e-11..1e-13; `(ab|ab)` Neumann vs
  importance-MC to 4.7e-3 (~1.3σ). **E₀ > exact −8.070 (variational ✓);** the +182 mHa gap is
  basis (crude single-ζ, no Li-2s/p polarization) + correlation (~83 mHa) → recovered by Stages 3-4.
- **STAGE 3 (part 1) DONE — σ² GROUND TRUTH + conditioning verdict (2026-09-21,
  `debug/lih_r12ci_sigma2_mc.py`):** σ² is exactly `Var_{|Φ₀|²}[F]`, a positive quantity, so a
  well-conditioned VMC (Metropolis on the block-factorized `|Φ₀|²=D_p(1,2)D_p(3,4)`) gives it
  with NO cancellation. **σ² = 0.13699 ± 2.3e-4** (geminal `exp(−0.5r)`). **Conditioning ratio
  σ²/F̄² = 0.039 → HEALTHY, NOT the Be trap** (Be's `f→1` gave `<<1e-3`); so `⟨F²⟩−F̄²` is safe
  in float64 and **no geminal switch is forced** (short-range `r·e^{−r}` only marginally better,
  0.055; geminal is now an accuracy question for Step 4, not conditioning). **Bonus — Stages 1&2
  independently cross-validated** by whole-determinant VMC: `⟨F⟩_MC=1.8812` vs analytic
  Fbar=1.8807 (1.1σ), `⟨V_ee⟩_MC=3.6169` vs analytic E2=3.6144 (1.2σ). Caveat: Metropolis
  acceptance 0.21 (tight Li core); cross-checks at ~1σ confirm equilibration, but sampler bias
  is the first suspect if the analytic port disagrees with σ²=0.137.
- **PoC ENERGY GATE PASSED via VMC (2026-09-21, `debug/lih_r12ci_vmc.py`).** The Stage-3
  conditioning verdict (σ² healthy) justified getting the ENERGY by the standard VMC linear
  method (one Jastrow parameter `c` on χ=F−F̄), all 2×2 elements well-conditioned `|Φ₀|²`
  averages, kinetic via the BOUNDED gradient form (`½Σ⟨|∇(χΦ₀)|²⟩`, no Laplacian cusp spikes).
  **E₀_VMC = −7.902 ± 0.014 reproduces the Stage-2 analytic −7.888 (1.0σ)** — validates the
  kinetic machinery + V_NN. Geminal scan (f=r·e^{−γr}, correct cusp sign): **best γ=0.5 →
  E_R12 = −7.932 ± 0.014 Ha, correlation captured dE = −29.5 ± 3.1 mHa** (~36% of LiH's ~83 mHa,
  ahead of Be's 19% PoC). **Variational both ways** (E₀ ≥ E_R12; E_R12 > exact −8.070). So
  explicit r₁₂ lowers the two-center LiH energy correctly-signed — the Be result at two centers.
  Traps hit + fixed: (1) forgot V_NN (a constant → shifts all eigenvalues equally, cancels in
  dE); (2) initial single-step Metropolis under-mixed the A↔B basins (acc 0.21) → mixture
  proposal (small steps resolve the tight core + occasional ~R hops) → acc 0.50, E₀ correct.
  **Scope (honest):** the 4-body content here is handled STOCHASTICALLY (implicit in the
  covariance estimators), NOT via the analytic RI-free reduction; single geminal, crude single-ζ
  ionic Φ₀ (not spectroscopic). The VMC confirms the PHYSICS + validates the stack; the analytic
  RI-free energy assembly below is the remaining distinctive (ahead-of-F12) deliverable.
- **STAGE 3 (part 2) DONE — the analytic RI-free σ² (2026-09-21, `debug/lih_r12ci_sigma2_analytic.py`).**
  **σ²_analytic = 0.13723 vs VMC 0.13699 (rel 1.7e-3, within the VMC's ±2e-4) — no MC, no RI.**
  Key reduction: the block density is SEPARABLE, `D_p² = P00⊗P11 + P11⊗P00 − 2 P01⊗P01`
  (`P_pq=m_p m_q`), so `⟨F²⟩` collapses to grid integrals of one-electron densities and their
  f-DRESSINGS `Ψ^f_h = ∫h(r')f(|r−r'|)dr'`. With orthonormal MOs: marginal `ρ=m0²+m1²`, block
  norm N=2 exact. `σ² = 2α₂+4β₂+8γ_sh+4γ_dj+16δ − 2α₁²−16α₁β₁−16β₁²` (verified: f=const→0).
  α/β/γ_dj are **scalar** f-interactions `I_f[P_pq,P_rs]=c·W·c` (W = 3×3 AO-pair f-tensor =
  Stage 1); only 3-body δ, γ_sh need the dressing FIELDS. **The one new piece — the two-center
  f-dressing `Ψ^f_{ab}`** (2D axially-symmetric convolution with the angle-averaged exp kernel;
  exp has NO singularity, unlike 1/r) — VALIDATED against Stage 1: gate W matches the f-tensor to
  1e-6..1e-8 (incl. `I_f[ab,ab]`), and F̄=2α₁+4β₁ reproduces 1.8807 to 2.5e-5. This is the
  RI-free machinery the "ahead-of-F12" claim rests on, now demonstrated inside a σ² number.
- **STAGE 4a (part 1) DONE — analytic h_Vne (2026-09-21, in `lih_r12ci_sigma2_analytic.py`).**
  `h_Vne = ⟨V_ne(F−F̄)⟩ = −0.9302` (analytic, RI-free) vs **VMC −0.9299 ± 0.012 (0.03σ) ✓**.
  Reduction: `⟨V_ne F⟩=4⟨v₁F⟩`, `⟨v₁F⟩=⟨v₁f₁₂⟩+⟨v₁⟩⟨f₃₄⟩+2⟨v₁f₁₃⟩+2⟨v₁f₂₃⟩`, all grid integrals
  of v-weighted densities × the existing f-dressings — **no new dressing** (V_ne is a known grid
  multiplier). Gates: `⟨V_ne⟩=−20.86474` matches Stage 2 exactly; constant-f gate h_Vne→0 verified
  analytically. VMC validation targets for the analytic port (geminal exp(−0.5r), from
  `lih_r12ci_vmc.py`): **h_T = +0.7485 ± 0.004, h_Vee = +0.3082 ± 0.002, g = −0.9303 ± 0.004**,
  σ² = 0.13696, h_total = +0.1268.
- **STAGE 4a (part 2) DONE — analytic h_T (kinetic×f) (2026-09-21, `debug/lih_r12ci_hT_analytic.py`).**
  `h_T = +0.749453` (analytic, RI-free) vs **VMC target +0.7485 (dev 0.0010) ✓.** Two parts:
  **Part A** = `½⟨χΣ|∇Φ₀|²⟩` reduces exactly like `σ²` but with **gradient one-electron densities**
  `G_pq=∇m_p·∇m_q` (grad-orb dot products via `cosAB=(ξ²+η²−2)/(ξ²−η²)`) alongside `P_pq`, reusing
  the existing `Ψ^f` dressings — `PartA = KI + ∫κ·Ψ^f_ρ + (α₁−F̄)(T₀₀+T₁₁)`, no new field. Gate:
  `T₀₀+T₁₁=8.36746=⟨T⟩_Φ₀` to rel 9e-14. **Part B** = the `∇F·∇Φ₀` coupling — the anticipated
  vector ∇f-dressing is UNNECESSARY: an IBP/Hermiticity identity collapses it to a SCALAR,
  `PartB = −½⟨Σ_{i<j}∇²f(r_ij)⟩ = −½γ²F̄ + γ·Ȳ`, where `Ȳ=⟨Σ e^{−γr}/r⟩` is a **Yukawa**
  (screened-Coulomb) two-electron expectation (`∇²f=γ²f−2γ·e^{−γr}/r`). Ȳ via a 3×3 Yukawa
  AO-pair matrix — isotropic aa/bb use the closed radial Yukawa potential (no grid singularity),
  ab,ab one importance-MC. Trap fixed: the Yukawa W-matrix MUST be in AO-pair order **[aa,ab,bb]**
  (matching `cvec`/`W^f`), not [aa,bb,ab] — the ordering bug flipped α₁^Y negative (MC caught it).
- **STAGE 4b (part 1) DONE — analytic h_Vee (Coulomb×f) (2026-09-21, `debug/lih_r12ci_hVee_analytic.py`).**
  `h_Vee = +0.309790` (analytic, RI-free) vs **VMC target +0.3082 (dev 0.0016) ✓** (MC cross-check
  +0.30863). `h_Vee = Cov[F,V_ee]` = the SAME separable bilinear as `σ²` but with the second kernel
  = Coulomb: a reusable **`cov_FA_FB(A,B,C)` engine** (SAME-pair uses the product kernel C=A·B;
  SHARE/DISJOINT use A- and B-dressings), self-checked to reproduce `σ²=0.13723` for `cov[f,f]`
  (GATE 2). Same-pair f·(1/r)=**Yukawa** (reuses h_T's W^Y); cross-pairs need the **two-center
  Coulomb dressing field** `Ψ^coul` — aa/bb = closed Hartree (`_hartree_1s`), transition-density ab
  = **prolate-Neumann POTENTIAL** `neumann_potential()` (extracted from the v5.15.10 energy
  expansion as `∫∫D·D'/r12=∫D·V_{D'}`). Gate 1: Neumann pot of ρaa vs closed Hartree rel 2.4e-3;
  `W^coul[aa,aa]=5ZA/8`, `[bb,bb]=5ZB/8` exact. **⟹ the full off-diagonal h = h_T+h_Vne+h_Vee =
  +0.749−0.930+0.310 = +0.129 is now analytic/RI-free** (VMC h_total +0.1268 ✓).
- **STAGE 4b (part 2a) DONE — analytic g_Vne (2026-09-21, `debug/lih_r12ci_gVne_analytic.py`).**
  `g_Vne = ⟨χ²V_ne⟩ = −2.8516` vs **MC −2.8490 (0.24σ), target −2.8454 ✓.** Method: block
  decomposition **F = S + I** (S=f₁₂+f₃₄ intra, I=inter f), `⟨F²V_ne⟩ = ⟨S²V⟩+2⟨SIV⟩+⟨I²V⟩ = A+B+C`,
  each a grid integral of **v-weighted dressed densities** (V_ne one-body → a `v`-weight on one
  electron; block-independence factorizes the rest). All three MC-validated (A=−3.942/MC−3.939,
  B=−26.96/MC−26.94, C=−49.24/MC−49.18). The −80.15→−2.85 cancellation is benign in float64.
  Reduction primitives: v-weighted dressings `Ψ^f_{vP}`, slices `S^f`,`S^{f,v}`, `ρ_vpartner`,
  down-field `G2`. This validates the block-decomposition trilinear machinery.
- **STAGE 4b (part 2b) DONE — the 4-body bridge INSIDE the energy + Yukawa-field machinery
  (2026-09-21, `debug/lih_r12ci_gVee_groundwork.py`).** (1) `⟨f₁₂f₃₄·coul₁₃⟩ = ¼ I_coul[S^f,S^f]`
  (S^f the f-dressed intra slice; two-center Coulomb between the dressed densities via
  `neumann_potential`) = **0.032814 vs block-sampler MC 0.032805 (0.0σ)** — the RI-free 4-body
  bridge, now demonstrated inside a real LiH correlation quantity (the showpiece). (2) the
  **Yukawa dressing field** `Ψ^Y_h = Ψ^coul_h − Ψ^{smooth}_h` (kernel `(1−e^{−γd})/d` non-singular,
  no Yukawa-Neumann needed), validated vs closed `yukawa_pot_iso` (rel 4.4e-3) — the enabler for
  the same-pair-Yukawa dressing slices `S^Y`. 6-product MC breakdown of ⟨F²V_ee⟩ banked in-file.
- **STAGE 4b (part 2c) OWED — the rest of g_Vee (the 3-body TRIANGLE), g_T, and the 2×2 → E_R12.**
  **The obstacle (genuine, not bookkeeping):** the full g_Vee's `⟨I²C⟩` products contain two-center
  **3-body triangle** terms, e.g. `⟨f₁₃·f₂₃·coul₁₂⟩` (two inter geminals sharing down-vertex 3,
  closed by intra `coul₁₂`). Down-integration gives a **non-separable 2-point kernel**
  `H₂(1,2)=∫ρ(3)f(r₁₃)f(r₂₃)d³r₃` (= `Kf · diag(geo·ρ) · Kf` as a matrix), and the residual
  `⟨…⟩ = ¼∫coul(r₁₂)·D_u(1,2)·H₂(1,2) d1 d2` is a 2-electron Coulomb of a NON-separable density —
  beyond the disjoint-pair bridge. **Identified path:** low-rank (SVD/eigen) expansion
  `H₂ ≈ Σ_l λ_l φ_l(1)φ_l(2)`, so `D_u·H₂` becomes rank-(3×modes) separable → the Coulomb reduces
  to `Σ I_coul[u_k φ_l, w_k φ_l]` via `neumann_potential` per mode (validate the truncation).
  The 6-product analytic build (P1=⟨S²C_S⟩, P2=⟨S²C_I⟩ both derived: P2 = 2∫Ψ^c_ρ S^{f²} +
  2∫S^f Ψ^c_{S^f} [the bridge]; P3=2⟨SIC_S⟩ = 4[∫Ψ^f_ρ S^Y + ∫Ψ^f_{S^coul} S^f]; P4,P5,P6 have the
  triangles) validates per-product against the banked 6-product MC targets. **g_T** = kinetic of
  `χΦ₀`, three sub-moments: `½⟨χ²Σ|∇Φ₀|²⟩` (=`⟨χ²·KD-density⟩`, reuses `KD`/`G_pq` from h_T Part A;
  target drift²=+1.168), `⟨χΣ∇F·v⟩` (cross, −0.048), `½⟨Σ|∇F|²⟩` (pair part `γ²(2α1^{f²}+4β1^{f²})`
  + a vector `Σ_j f'(r_ij)r̂_ij`-field 3-body part; target +0.266). Then the analytic 2×2
  [[E0,h],[h,g]] / [[1,S01],[S01,σ²]] → E_R12, compare to the VMC −7.932 (closing the loop two ways).
  **Conditioning:** `g_Vne = ⟨F²V⟩−2F̄⟨FV⟩+F̄²⟨V⟩ ≈ −80.15+151.10−73.80` — ~3–4% residual, benign
  in float64 (confirmed this session). `⟨F²V⟩` wants ~1e-3 relative accuracy; the current 72×44
  grid delivered g_Vne to 0.24σ.
The integral machinery (σ + π/δ, + the Yukawa field + the 4-body bridge inside the energy) is built and validated.

---

## 0. The goal, stated precisely

**Primary: LiH TOTAL-ENERGY correlation accuracy** (the original prompt — "more physically
accurate"). Add explicit r₁₂ correlation to a two-center LiH CI and measure how close the
total energy gets to the exact non-relativistic value (LiH exact ≈ −8.070 Ha; the "physical
accuracy" target is chemical accuracy, 1.6 mHa, or better). This is the axis the Be R12-CI
PoC pointed at.

**Secondary: LiH R_eq (bond length).** Best so far 5.3% (composed, Paper 17). The prolate
variational-core work (v5.15.7) confirmed the frozen-core-is-culprit / variational-core-is-cure
mechanism but the clean R_eq is blocked by a **numerical (grid) wall**, NOT the 4-body wall —
so explicit r₁₂ does not directly fix R_eq. Keep these two goals separate; aim at total energy
first.

**Do NOT conflate with the Be atom result.** Be (this session) is a 1-center atom; LiH is a
2-center molecule. The Be number (E_R12 = −14.557, ~19% of correlation) was a PoC that the
4-body RI-free machinery works inside an energy — it is nowhere near physical accuracy and was
never meant to be.

---

## 1. What this session ESTABLISHED that unblocks LiH (the load-bearing facts)

1. **The N=4 explicit-r₁₂ wall is SOFT** (v5.15.8, `docs/walls/register.md`). The only
   genuinely-4-body term in ⟨Φ|F H F|Φ⟩ is the scalar Coulomb chain `f₁₂ f₃₄ / r₁₃` (kinetic
   gradients are pair-local; V_ne is one-body → only the two-body Coulomb bridges disjoint
   pairs). It is **exact, RI-free, terminating** (bridge multipole sum stops at L ≤ 2·l_bridge)
   and **reducible** (Legendre addition theorem factorizes the chain across the bridge into a
   finite (L,M) contraction of vertex kernels). So an **all-electron LiH explicit-r₁₂** with
   the 4-body term is NOT blocked at the integral level — this is exactly the "N=4 explicit-r₁₂
   wall" the prior LiH work (v5.15.6) said the all-electron build would hit. It is passable.
2. **A Be R12-CI energy works** with the ill-conditioned pieces done analytically
   (`debug/be_r12ci_full.py`): E_R12 = −14.5572 Ha, −18 mHa correlation, variational.
3. **Ill-conditioning lesson (CRITICAL — carry this to LiH):** the linear-geminal energy is
   pathological with a long-range geminal (f = 1−e^{−r} → basis 99.8% parallel, F̄≈5, the
   coupling h a 0.03 residual of ~71-magnitude numbers → needs ~1e-4 precision, defeats all MC).
   **Fix: short-range geminal** f = r·e^{−γr} (cusp f′(0)=1, →0 at large r; parallelism
   0.994→0.86, F̄ 4.9→0.3), + **orthogonalized basis** G = (F−F̄)Φ₀ with every small quantity
   computed directly as a variance/covariance (no cancellation), + **exact analytic overlaps**
   (σ² must be exact). Then h and σ² analytic-exact, g well-conditioned by MC.

---

## 2. What TRANSFERS vs what is NEW

**Transfers (reuse the ideas):**
- The 4-body reduction *structure*: term enumeration by up/down block + cross, the
  bridge-factorization of the disjoint-pair chain, the "kinetic gradients don't bridge, only
  Coulomb does" argument. See `debug/be_r12ci_full.py` (TWOPAIR framework) and
  `debug/r12ci_4e_wall_diagnostic.py`.
- The orthogonalized-basis energy formulation + short-range-geminal lesson (Sec. 1.3).
- The analytic-overlaps-are-mandatory lesson (the ill-conditioning is in the overlap; make S
  exact).

**NEW — the actual build (the Be engine does NOT transfer):** Be is *atomic* — all-s orbitals,
isotropic densities, **monopole** kernels, spin-block factorization. LiH is **two-center** —
densities are not isotropic, so none of the monopole reductions apply. LiH needs the
**prolate spheroidal two-center explicit-r₁₂ machinery** with **Neumann kernels** (A K₀ − B K₁),
which already EXISTS for 2 electrons and must be extended to 4.

---

## 3. Reusable artifacts (the substrate to build ON)

| File | What it is | Role in the LiH build |
|:--|:--|:--|
| `debug/prolate_r12_mpf.py` | **HeH⁺ 2-electron prolate explicit-r₁₂ engine** (`assemble_hetero(basis_p,R,alpha,Z_A,Z_B,...)`, `vne_hetero_mpf`, `vee_r12_odd_mpf`, `kinetic_p1p1_mpf`, the odd-r₁₂ `_odd_K0/_odd_K1/_make_godd` machinery) | **THE substrate.** The 2-electron heteronuclear prolate r₁₂ integrals are done here. Extend to 4 electrons. |
| `debug/prolate_allelectron_fci.py` | All-electron (N=4) variational-core prolate FCI, **pairwise V_ee** (no r₁₂), generalized-m π ERIs | The 4-electron LiH CI *structure* without r₁₂ (HF==eckart validated). Add r₁₂ to this, OR use it as the CI reference the r₁₂ correction sits on. **Note its numerical wall (Sec. 5).** |
| `debug/be_r12ci_full.py` | Be analytic block-reduction engine + orthogonalized-basis energy | Template for the reduction structure + the ill-conditioning handling. |
| `geovac/transcorrelated_sturmian.py`, `debug/r12ci_3e_{triangle_kernel,vertex_rules}.py` | N≤3 **atomic** r₁₂ angular rules (RULE A / RULE B / TRIANGLE) | Reference for the term inventory; the angular rules need the prolate analog (the *topology* is the same, kernels differ). |
| `debug/heh_converge.py` | HeH⁺ convergence driver | Pattern for driving the 2-center r₁₂ CI + the lesson "the lever is the ANGULAR basis, not the exponent." |
| `debug/track_logs/prolate_native_lih.md` | Full LiH prolate history | The numerical-wall details, the frozen-core vs variational-core story, HeH⁺ reference (Kolos-Peek). |

---

## 4. The build path (concrete)

**Step 0 — decide the ansatz.** Two options:
- **(A) Explicit-r₁₂ CI (Hylleraas-style, r₁₂ in the basis)** with the 4-body term. This is the
  direct analog of the Be R12-CI and the HeH⁺ engine; the 4-body term is now known soft. Cleanest
  demonstration of "the RI-free 4-body evaluation inside a LiH energy."
- **(B) CI + r₁₂ correction:** take `prolate_allelectron_fci.py`'s pairwise-V_ee CI and add an
  explicit-r₁₂ correction (à la the Be orthogonalized G = (F−F̄)Φ₀). Less integral work; may be
  the faster PoC.
  Recommend **(B) for the first PoC** (reuses the validated 4e CI), then **(A)** for accuracy.

**Step 1 — 2→4 electron extension of the prolate r₁₂ integrals.** `assemble_hetero` does the
2-electron heteronuclear prolate r₁₂ blocks (overlap even, V_ee odd, V_ne hetero, kinetic). For
4 electrons you additionally need: the **3-body** terms (shared-vertex; RULE-A/TRIANGLE analogs
in prolate kernels) and the **one 4-body** disjoint-pair chain `f₁₂ f₃₄ / r₁₃` — evaluated by the
**bridge factorization** (validated on Be) with the prolate **Neumann** kernels replacing the
atomic monopoles. The bridging Coulomb multipole sum terminates (soft wall), so it is finite.

**Step 2 — handle ill-conditioning from the start** (do NOT skip): short-range geminal
f = r·e^{−γ r}; orthogonalized basis; exact analytic overlaps (the two-center overlap/S_11
analog). Compute the small quantities (h, σ²) as variances/covariances. This is the difference
between a trustworthy number and noise (see Be Sec. 7b of `sprint_n4_wall_diagnostic_memo.md`).

**Step 3 — PoC gate:** at a FIXED LiH geometry (R = R_eq), does adding r₁₂ lower the total energy
correctly (variational, correctly-signed), like Be? Validate the r₁₂ integrals reduced==brute on
one 4-electron LiH integral first (as `be_r12ci_4body_exchange.py` did for Be). Then quote the
correlation captured.

**Step 4 — accuracy push (only after the PoC):** bigger basis, higher r₁₂ powers (p≥2 — the Be
`p≤1` plateaued at ~19%; H₂ needed p≥2 for µHa), multiple length scales. Target chemical accuracy
on the LiH total energy at fixed geometry.

---

## 5. Traps (things this and prior sessions hit — do not re-hit)

- **Ill-conditioning (Sec. 1.3).** The #1 killer. Short-range geminal + orthogonalized basis +
  exact analytic overlaps, from the start.
- **The numerical/grid wall (v5.15.7, the R_eq blocker).** The tight Li 1s² core needs
  real-space grid resolution that, in the prolate coordinate r_A=(R/2)(ξ+η), is tied to R — so an
  isolated Li²⁺ core swings 0.36 Ha across R (grid artifact, NOT BSSE), swamping the 0.088 Ha bond
  and collapsing R_eq scans inward. This bites **R_eq scans**, less so a **single-geometry total
  energy** (the target here). Use the mpf engine + a core-adapted/graded grid; if scanning R,
  budget for this. A dedicated atom-centered tight-core treatment is the real fix (scoped, unbuilt).
- **Frozen core reintroduces the R_eq drift.** A frozen Li 1s² gives the +5.5% outward rigidity
  (v5.15.6); the cure is a variational (all-electron) core, which hits the grid wall above. For
  TOTAL ENERGY at fixed geometry this is less acute; for R_eq it is the crux.
- **Don't chase R_eq with r₁₂.** r₁₂ improves correlation/total energy; the R_eq error is a
  core-screening/grid problem. Different axes.

---

## 6. Success criterion & honest scope

- **PoC success:** adding explicit r₁₂ to a two-center LiH CI lowers the total energy, variationally
  and correctly-signed, with the 4-body term evaluated by the soft-wall reduction (reduced==brute
  validated on one integral). This is the direct LiH analog of the Be result.
- **Physical-accuracy success:** LiH total energy within chemical accuracy (1.6 mHa) of exact
  (≈ −8.070 Ha) at a fixed geometry — needs the accuracy push (Step 4).
- **Honest scope:** this is a MAJOR build, **bigger than Be** (two-center, 4-electron prolate
  r₁₂). Expect multiple sessions. First fresh session: aim for the PoC gate (Step 3), NOT immediate
  physical accuracy. Do not start it at the tail of a long session (this session's near-misses were
  all budget-driven).

---

## 7. First moves for the fresh session

1. Read: this file, `debug/track_logs/prolate_native_lih.md`, Papers 12 & 19, `memory/
   r12_generalization_boundary_n3_n4.md`, `memory/polyatomic_state_of_play.md`. Verify current
   state (CHANGELOG since 2026-09-21).
2. Open `debug/prolate_r12_mpf.py` (`assemble_hetero`) and `debug/prolate_allelectron_fci.py` —
   understand the 2-electron r₁₂ blocks and the 4-electron pairwise-V_ee CI.
3. Pick ansatz (B recommended for PoC), and do Step 1's 2→4 extension + Step 3's single-integral
   reduced==brute validation BEFORE any energy. Short-range geminal + orthogonalized basis from
   the start.
