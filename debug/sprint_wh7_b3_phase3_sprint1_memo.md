# Sprint B3 Phase 3, Sprint 1 — state-level cost structure on modular orbits (2026-06-10, v3.115.0)

**Goal:** lay the foundation for the state-level Lorentzian construction fixed by the
Phase-2 negative: intervals from the modular flow, costs from entropy production,
super-additivity from divergence monotonicity. Cost functional = Datta max-divergence
D_max per the standing Paper 49 lesson (memory `m1_datta_max_divergence_replacement`:
Umegaki chain fails on cocycle triples; D_max's chain is a theorem).

**Verdict: FOUNDATION-LAID — one universal fact, one reproducible conditional
structure, one honest caveat.**
Driver `debug/wh7_b3_phase3_state_intervals.py` (deterministic seeds, JSON in
`debug/data/`), falsifier `tests/test_wh7_b3_phase3.py` (7/7), Paper 45 Q1 extended,
24 pp GATE: PASS.

## Setup
j ≤ 1 PW window (dim 14); K = 2J_z (Paper 42 `two_m_j` convention); wedge KMS state
ρ = e^{−βK}/Z, β = 1 (KMS verified directly; β/κ_g conventions documented here, not
load-bearing). Orbits ω_t = U_t ω₀ U_t† of a generically perturbed state; triples
σ₁ = ω₀, σ₃ = ω_1, midpoint kicked by e^{iεG} with G a Hermitized causal-class
representative from Phase 2.

## Results

1. **Anchors (bit-level):** KMS residual 1.7×10⁻¹⁵ (after catching a σ_{iβ} sign
   error in the first run's check — the check was wrong, not the state); flow
   invariance 5×10⁻¹⁷; orbit period 2π at 2.5×10⁻¹⁶; orbit injectivity on (0, 2π)
   (min trace-distance 0.099). Thermal-time intervals are well-defined mod 2π on
   orbits — circular time, as a compact (Euclidean/thermal) carrier must give.

2. **Universal (the twin direction):** the D_max chain inequality
   c(σ₁,σ₃) ≤ c(σ₁,σ₂) + c(σ₂,σ₃) held in **96/96** cells (7 classes × 2 generators
   × 4 ε + 40 random kicks) — "a detour never costs less," verified on the substrate.
   Umegaki also showed 0/96 violations on THIS panel (the Paper 49 ~5% failure was on
   TICI cocycle triples specifically; not reproduced here, panels differ — noted, not
   contradictory).

3. **Conditional structure (reproducible, mechanism open):** the kick-attributable
   excess cost scales bimodally with the kick generator's boost weight:
   p ≈ 2 (quadratic) for |m′| ≤ 1/2 classes; p ≈ 1.1–1.3 for |m′| ≥ 1 classes —
   and the null-class kick carries the largest coefficient (C ≈ 1.09, ~2× others).
   **Not a reparametrization confound:** projecting out the flow-tangent component
   of every kick shifts the exponents by < 0.1 (tangent overlaps all ≤ 0.15).
   Suspected mechanism: D_max's max-eigenvalue non-smoothness + a weight selection
   rule on the top eigenspace; identification = Sprint 2.

4. **Honest caveat (frozen as a test):** the SIGN of the excess over the on-orbit
   baseline is reference-state dependent (T5: fresh seeds produce zero/negative
   excess in some class cells, e.g. (2,0) at seed 202). Whenever the excess is
   positive, the p-split reproduces ((2,0) → 1.97/2.01; (2,2) → 1.23/1.26). The
   baseline deficit itself is large (2.375 nats at this configuration) — the on-orbit
   chain is far from tight, consistent with D_max measuring distinguishability, not
   arc length. **Implication for the interval design:** D_max is the right COST
   functional (its chain is the universal twin direction) but is not itself the
   interval; the interval must be the flow parameter, with costs entering as the
   deficit/penalty layer — matching the Paper 48 bridge design (ℓ = κ_g·τ_mod on
   orbits) rather than a divergence-as-distance reading.

## Sprint 2 (named)
- Ensemble statistics over reference states (the excess-sign question becomes a
  distributional statement; identify the state classes where positivity holds).
- Cost-functional comparison: D_max vs fidelity/Bures vs symmetrized divergences —
  which produces a state-independent positive excess (if any)?
- Mechanism of the bimodal p (top-eigenspace weight selection rule).
- Wedge restriction (HemisphericWedge of `geovac/modular_hamiltonian.py`) as the
  proper substrate; then the interval functional ℓ = flow-parameter on orbits with
  cone-graded admissibility, and the convergence statement (the Phase-3 prize).
