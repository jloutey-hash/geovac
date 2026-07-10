# Sprint GB-4 — the discriminator: two-window witnessing + the ladder's boundary

**Date:** 2026-05-29. **Type:** diagnostic. **Verdict:** PARTIAL POSITIVE with a sharp boundary. The two-window (functional-equation) structure is witnessed at B₄, B₆ from both sides by physically distinct observables (Casimir + thermal); the combinatorial packing core is cleanly OFF the ladder. NOT numerology, NOT new physics — a precise map of where the ladder lives.

## 1. The test

The Casimir landings (GB-3) were all ζ at *negative* integers (skeleton window). The ladder's core claim is that each rung B_{2k} has *two* windows glued by the functional equation: ζ(1−2k) skeleton AND ζ(2k) transcendental. Sharp test: is the **positive-even (M2) window** lit, by *different physics* than the Casimir, on the same rungs? Plus a negative control: the genuinely-GeoVac combinatorial core (from the packing axiom, not generic sphere QFT) must be OFF, or the ladder is just "everything carries a zeta."

## 2. Results (exact sympy)

**Thermal window lit by Stefan–Boltzmann (Bose–Einstein integral I_d = ∫x^d/(e^x−1) = Γ(d+1)ζ(d+1)):**
- S³×S¹ (4D blackbody): I₃ = π⁴/15 → ζ(4) = π⁴/90 = **B₄, M2 window** (cross-checks Sprint TD Track 1's −π²/90 T⁴).
- S⁵×S¹ (6D blackbody): I₅ = 8π⁶/63 → ζ(6) = π⁶/945 = **B₆, M2 window**.

**Two-window table — each rung lit from both sides by distinct physics:**

| rung | skeleton ζ(1−2k) | from | M2 ζ(2k) | from |
|:----:|:----------------:|:-----|:--------:|:-----|
| B₄ | ζ(−3)=1/120 | S³ Casimir (vacuum) | ζ(4)=π⁴/90 | S³×S¹ Stefan–Boltzmann (thermal) |
| B₆ | ζ(−5)=−1/252 | S⁵ Casimir (vacuum) | ζ(6)=π⁶/945 | S⁵×S¹ Stefan–Boltzmann (thermal) |

The two windows ARE the functional equation = **temperature-inversion duality** (T↔1/T on the thermal circle): vacuum Casimir (negative-integer ζ) ↔ blackbody (positive-even ζ). Same Bernoulli rung, both windows, physically distinct observables.

**Negative control — combinatorial core is OFF the ladder:**
- B = 42 (α Casimir trace): not on a rung.
- Δ = 1/40 (α): not on a rung.
- κ = −1/16 (Fock Jacobian): not on a rung.
- (Sprint TD Track 5: GeoVac correlation entropy S_full(GS) PSLQ-null off the engine.)

## 3. Honest verdict — what it does and does NOT prove

- **Beats numerology:** the negative control is clean. The ladder is NOT "everything in GeoVac carries a zeta" — the packing-combinatorial core is off it.
- **Does NOT beat "sphere QFT carries zeta" in the strongest sense:** both Casimir and Stefan–Boltzmann are sphere-QFT objects. A hard skeptic can call the positive landings textbook (Planck blackbody + Casimir + Cardy temperature inversion). Acknowledged.
- **What it genuinely establishes:** the ladder is the **temperature-inversion / functional-equation structure of GeoVac's spectral-geometry sector** (gravity Casimirs + thermal sector), witnessed at three rungs from BOTH windows by distinct physics — and it is **bounded** away from the combinatorial packing core.

## 4. Cross-connection: this re-derives the Phase 4G no-common-generator finding

Phase 4G found K = π(B + F − Δ) has "no common generator": B = combinatorial Casimir trace, F = ζ(2) arithmetic-transcendental, Δ = combinatorial boundary product. GB-4's ladder/off-ladder split IS that split: **F = ζ(2) is on the ladder (B₂, M2 window); B and Δ are off it.** The α-conjecture is anomalous precisely because it sums one on-ladder quantity (F) with two off-ladder combinatorial ones (B, Δ). The ladder boundary and the no-common-generator finding are the same fact from two directions.

## 5. The two-layer picture this sharpens

GeoVac has two structurally distinct layers:
- **Spectral-geometry / zeta-ladder layer:** Casimir + thermal + heat-kernel, organized by the quadratic Weyl multiplicity onto the Bernoulli ladder (both functional-equation windows). This is also where the no-FE/no-RH wall lives (GB-2).
- **Combinatorial packing layer:** B, Δ, κ, graph zeta — off the ladder, the genuine α-conjecture territory.

F = ζ(2) is the single bridge between them (it is both a packing-degeneracy Dirichlet value, D_{n²}(4), and a ladder rung).

## 6. Documentation applied
- Paper 35: `obs:two_window_duality` added (Casimir↔thermal = temperature inversion = functional equation; negative control; boundary).

## Files
- `debug/sprint_gb4_discriminator.py`
