# Next-session build plan — the rigorous variational LiH energy (VMC-over-FCI)

**Written 2026-09-22 (v5.15.21), for the next session. PI-directed.** Self-contained so a fresh
context can execute without relearning. Current-state rule applies: before starting, read the
owning paper section (Paper 19 "Fixed-geometry energy versus well shape") + CHANGELOG since
2026-09-22, in case anything moved.

---

## THE GOAL
Turn the **additive-F12 estimate** LiH energy (−8.062, near-chemical) into a **fully rigorous,
variational** number. Target: a variational E ≈ −8.06 ± a few mHa that (a) beats the pure-orbital
ceiling −8.032, (b) stays **above** the exact −8.0705 (variational — a value below exact is a
BUG, the frozen falsifier), and (c) is double-counting-free by construction. Success confirms
(or sharpens) the −8.062 estimate as a proven number and closes the LiH energy arc.

## WHY THIS PATH (VMC-over-FCI), not the analytic route
The rigorous-energy feasibility fork (2026-09-22, `debug/sprint_lih_rigorous_f12_memo.md`)
found BOTH rigorous paths are multi-session, and chose one:
- **(A) analytic F12-on-correlated-reference — DO NOT ATTEMPT.** `geovac/lih_r12ci` is welded to
  the 2-MO |1s_A²1s_B²| reference at the *formula* level (P00/P11/P01, c00/c11/c01/crho as
  module constants; DP_COMPS/FkVne/FkSv2 4-body enumerators loop over exactly those; E0
  hardcoded). Retargeting = re-expressing σ²/h/g over the FCI's 1–4-particle RDMs (F12-on-MRCI);
  the ⟨f_ij f_kl⟩ 4-body over an arbitrary MO basis is the piece that doesn't generalize.
- **(B) VMC-over-FCI — THE PLAN.** Variational and double-counting-free by construction
  (E_VMC[Ψ_CI · Jastrow] ≥ exact). Sidesteps the F12-on-MRCI re-derivation. Tractable but a
  from-scratch VMC stack.

## CURRENT STATE (the load-bearing facts the build rests on)
- **Geometry SOLVED:** Route C prolate all-electron LiH R_eq = +0.2% (π-converged, corpus best).
  Do NOT touch it. The energy is the only open axis.
- **Orbital energy ceiling −8.032** (M>16 exhausted; `debug/prolate_float_eri.py push`). A finite
  orbital basis cannot represent the e-e cusp — this is the wall VMC's Jastrow must cross.
- **Additive-F12 estimate −8.062** = −8.032 + core cusp 30.1 mHa (+ ~8 mHa budget-pinned valence).
  The core cusp 30.1 mHa (He-like, validated vs known He to 1.4 mHa) is the number the VMC must
  recover **variationally** (cross-check).
- **Exact LiH BO energy = −8.0705 Ha** (Paper 19's value; 3.015 bohr).
- **Productionized engine (fast, validated):** `debug/fci_fast.py` (sparse connected-pair FCI,
  bit-exact vs dense, ~14× at M=16) + `debug/prolate_float_eri.py` (float64 ERI, energy-exact to
  ~5 µHa). These made M>16 reachable and are the substrate for the FCI eigenvector + energy.
- **The validated geminal (Jastrow factor u):** the cusp-correct `f = r·e^{−γr}` (linexp, γ≈0.5)
  and `f = e^{−γr}` (exp), both in `geovac/lih_r12ci` (the 2-MO PoC gave dE −29.5 mHa VMC).

## THE BUILD (6 steps; validate each before the next)
1. **FCI eigenvector.** Flip `fci_fast.fci_energy_fast` to `return_eigenvectors=True`
   (eigsh already computes it — just return it) + the orthonormal-MO coefficients from
   `assemble_rebased` (the X @ ... transform). GATE: the returned vector's Rayleigh quotient ==
   the eigenvalue.
2. **Arbitrary-point prolate orbital evaluator + real-space gradient.** Only fixed-grid
   evaluators exist (`get_orbital_on_grid`/`_orb_on_grid`). Need φ(r) and ∇φ(r) at ANY r, via the
   chain rule through (ξ, η, φ_azimuth) — r = (R/2)(ξ+η) etc. GATE: finite-difference ∇φ matches
   the analytic ∇φ; the evaluator on the grid points matches the existing grid evaluator.
3. **Multi-determinant Ψ_CI(r₁..r₄) evaluator + ∇Ψ_CI/Ψ_CI.** Slater determinants from step-2
   orbitals, weighted by the step-1 CI coefficients. GATE: ⟨Ψ_CI|Ψ_CI⟩ and ⟨Ψ_CI|H|Ψ_CI⟩ by VMC
   (no Jastrow) reproduce the FCI energy −8.032 within statistics — see the MANDATORY gate below.
4. **Jastrow J = exp(Σ_{i<j} u(r_ij))** with u = the validated linexp/exp geminal + ∇J.
5. **VMC:** Metropolis on |Ψ_CI · J|²; the **gradient-form** kinetic local energy
   (½Σ|∇(Ψ)|²/|Ψ|², no Laplacian cusp spikes — the corpus's bounded form); linear-method (or a
   1-param scan) to optimize γ in u.
6. **THE MANDATORY GATE (do NOT skip):** VMC of Ψ_CI *without* Jastrow must reproduce −8.032
   within statistics BEFORE trusting the Jastrow result. A half-validated VMC is an untrustworthy
   number — the task forbids reporting one. Then correlated sampling for the correlation dE
   (±few mHa, as the 2-MO VMC achieved).

## DECISION GATE
- **GO:** rigorous variational E ≤ −8.05, > −8.0705, gate-6 passed → the rigorous near-chemical
  LiH number. Report it; it confirms/sharpens −8.062.
- **BORDERLINE:** E in [−8.05, −8.04] → partial; report + diagnose (likely γ/basis under-optimized).
- **STOP:** if the VMC stack can't be validated to the −8.032 gate in budget → report the specific
  blocker; do NOT report an ungated VMC number.

## VALIDATION ANCHORS
−8.032 (must beat; the gate-6 no-Jastrow target) · −8.062 (should land near) · −8.0705 (must stay
above — variational, frozen falsifier) · core cusp 30.1 mHa (the correlation the Jastrow recovers).

## TRAPS / DO-NOT-REDERIVE (institutional memory)
- Analytic F12-on-MRCI is BLOCKED (above) — do not attempt.
- The additive −8.062 estimate is already CONFIRMED defensible — do not re-derive it; the VMC's
  job is the *rigorous* number, expected near it.
- The geometry (+0.2%) is solved — do not touch.
- Variational integrity: E > −8.0705 always. Below = a bug (kinetic sign, normalization,
  double-counting). See `memory/feedback_bit_exactness_rule.md`.
- Single-det + Jastrow caps below −8.032 (~HF+cusp) — Ψ_CI MUST be multi-determinant to beat it.

## KEY FILES
- `debug/sprint_lih_rigorous_f12_memo.md` — the fork's full feasibility verdict + this recipe.
- `debug/sprint_lih_additive_f12_memo.md` — the additive estimate (−8.062) + the 30.1 mHa cusp.
- `debug/track_logs/prolate_native_lih.md` — the full LiH prolate chronicle.
- `debug/{prolate_allelectron_c4,prolate_float_eri,fci_fast,lih_additive_f12,lih_r12_ceiling_probe}.py`
  — the engines (orbital C4, float ERI, sparse FCI, additive-F12, the He-like cusp probe).
- `geovac/lih_r12ci/` — the validated geminals (the Jastrow u).
- Paper 19 "Fixed-geometry energy versus well shape" — the owning paper section.

## HONEST SCOPE
A from-scratch VMC stack (orbital evaluator + gradients, multi-det sampler, Jastrow, Metropolis,
optimizer) — multi-session and bug-risky. The gate-6 −8.032 reproduction is the trust anchor:
without it, any VMC number is untrustworthy. Do NOT start it at the tail of a long session.
Expected outcome: a rigorous E ≈ −8.06 ± few mHa, confirming the additive estimate.
