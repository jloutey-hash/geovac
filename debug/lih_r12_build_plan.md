# LiH explicit-r₁₂ build plan — handoff for a fresh session

**Written 2026-09-21 (end of the N=4/Be session), for the next session.** The top-of-session
intent was **LiH physical accuracy**; the session validated the enabling machinery on the
simpler atomic Be case and stopped there by design. This is the concrete plan to carry it to
LiH. Read this + the owning papers (12, 19) + `debug/track_logs/prolate_native_lih.md` before
starting. Current-state rule applies: verify against CHANGELOG since 2026-09-21 first.

---

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
- **STAGE 2–4 owed (in order):** `E₀` (determinant energy: two-center T + V_ne + Coulomb J/K +
  V_NN — a mini two-center HF, well-conditioned); `σ² = ⟨F²⟩−F̄²` and `h = ⟨HF⟩−F̄E₀` (the
  **ill-conditioned** pieces — small residuals of large numbers, so quadrature/near-exact, per
  the Be lesson; **do these with fresh care, not at a session tail**); `g = ⟨G|H|G⟩` with the
  f₁₂(1/r₁₃)f₃₄ chain via the validated σ+π/δ reduction. Then the 2×2 → E_R12.
The integral machinery (σ + π/δ) the 4-body term needs is built and validated.

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
