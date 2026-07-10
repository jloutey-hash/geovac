# Sprint GB-3 — the B₆ rung is populated (S⁵ conformal-scalar Casimir)

**Date:** 2026-05-29. **Type:** diagnostic / prediction-test. **Verdict:** B₆ RUNG POPULATED — the Bernoulli ladder now has three forced rungs (B₂, B₄, B₆), with the spreading independently re-confirming the GB-2 mechanism.

## 1. Principled target (not cherry-picked)

The conformal-scalar Casimir on S^{2k−1} leads with ζ(−(2k−1)): S¹→ζ(−1) (B₂), S³→ζ(−3) (B₄), **S⁵→ζ(−5) (B₆)**. S⁵ is a real GeoVac manifold (Bargmann–Segal, Paper 24), so the rung is reachable *within* the framework. Method: E = (1/2)Σ ω_l·deg_l (zeta-regularized), ω_l = l+(d−1)/2, deg_l = degree-l harmonic dimension on S^d.

## 2. Results (exact sympy)

- **Control (validates method):** S³ conformal scalar Casimir = (1/2)Σn³ = (1/2)ζ(−3) = **1/240**. Reproduces Paper 35 exactly.
- **Target:** S⁵ conformal scalar Casimir = **−31/60480** = **(1/24)(ζ(−5) − ζ(−3))**.
  - ζ(−5) = −1/252 (B₆ rung) present with clean coefficient **+1/24**.
  - ζ(−3) = 1/120 (B₄) present with coefficient **−1/24**.
  - ζ(−1) (B₂) **absent** (specific structural content — the spectrum polynomial (n⁵−n³)/12 has no n¹ term).

## 3. Why it spreads (GB-2 mechanism, independent confirmation)

The S⁵ harmonic degeneracy is quartic (n²(n²−1)/12), so ω·deg is degree-5 in n → regularizes to ζ(−5) AND ζ(−3). On S³ the degeneracy is n² → ω·deg is degree-3 → single clean rung (ζ(−3)). This is the *same* multiplicity-spreading that GB-2 found obstructs the functional equation — here seen on the Casimir (skeleton) side, from a completely different observable. Two findings, one mechanism.

## 4. Self-similar M2/M3 split at the rung

- Integer (conformal scalar) → clean M2 ladder: −31/60480 = (1/24)(ζ(−5) − ζ(−3)).
- Half-integer (HO / Bargmann S⁵, TX-B) → M3 sibling: −17/3840 (the "17" echoes Dirac-S³ 17/480).
- M2-window partner of ζ(−5): ζ(6) = π⁶/945 — predicted π⁶ transcendental in the S⁵ sector's observable window (Stefan–Boltzmann-class on S⁵×S¹).

The integer/half-integer = scalar/spinor = M2/M3 split that appeared at the B₄ rung (S³ scalar 1/240 vs Dirac 17/480) recurs at B₆ (S⁵ scalar −31/60480 vs HO −17/3840). The ladder structure is self-similar across rungs.

## 5. Audit (honest bound)

- Free params: zero; forced computation.
- Control: S³ reproduces Paper 35's 1/240.
- **Skeptic's point (acknowledged):** "S⁵ Casimir carries ζ(−5)" is partly textbook spectral geometry. The GeoVac-specific content is NOT the bare ζ(−5) appearance but: (a) S⁵ is in the framework, so the rung is reached internally; (b) the *specific* rung-content (B₆ and B₄, weights ±1/24, B₂ absent) is a forced prediction; (c) the spreading independently re-derives GB-2 from a second observable; (d) the M2/M3 self-similarity across rungs.
- Robustness: exact.

## 6. Net

The Bernoulli ladder promotes from a two-rung observation to a **three-rung structure with a mechanism (GB-2) that predicts the spreading and is independently confirmed**. The reading (two-layer split = ζ functional-equation reflection, per Bernoulli order) is now supported across three rungs and two observable types (functional-equation obstruction + Casimir). It remains a structural correspondence, not a forcing theorem; classical RH stays closed (GB-2). Next falsifiable rung: S⁷ conformal Casimir predicted to spread across ζ(−7), ζ(−5), ζ(−3) (degree-7).

## 7. Documentation applied
- Paper 32 §VIII `rem:bernoulli_ladder`: "B₆ empty" line updated → B₆ populated by S⁵ Casimir.
- Paper 35: S⁵ conformal-scalar Casimir + dimensional ladder added.

## Files
- `debug/sprint_gb3_b6_rung_hunt.py`
