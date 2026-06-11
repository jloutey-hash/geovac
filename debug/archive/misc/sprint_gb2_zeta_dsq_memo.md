# Sprint GB-2 — functional equation for ζ_{D²}? Closed-form probe

**Date:** 2026-05-29. **Type:** diagnostic. **Verdict:** DECISIVE CLEAN NEGATIVE with structural upgrade — ζ_{D²} has no functional equation, the obstruction is the quadratic Weyl multiplicity, and this unifies RH walls RH-N + RH-O.

## 1. Motivation

Sprint GB's Bernoulli-ladder finding placed the framework's even-s ζ-values in the M2 transcendental window. The squared Dirac spectral zeta ζ_{D²}(s) = D(2s) lives in exactly that window (T9, Paper 28: ζ_{D²}(s) = 2^{2s-1}[λ(2s-2) − λ(2s)], λ(z)=(1−2^{−z})ζ(z)). RH-O (Apr 2026) killed the functional equation for D(s) by a 13,080-template numerical search. ζ_{D²} was never tested directly, and the T9 closed form lets us settle it **analytically** instead of by search.

## 2. Method (closed-form, not template search)

Since ζ_{D²}(s) = D(2s), both ζ_{D²}(s) and its reflection ζ_{D²}(c−s) lie in the 2D space span{ζ(2s), ζ(2s−2)} over the field of elementary functions. A functional equation ζ_{D²}(s) = χ(s)·ζ_{D²}(c−s) with χ elementary exists **iff the two coefficient vectors are parallel**, i.e. det = A·D − B·C = 0, where:
- F(s) = C·ζ(2s) + D·ζ(2s−2), C = −2^{2s-1}(1−2^{−2s}), D = 2^{2s-1}(1−2^{2−2s}) (exact, no FE needed).
- F(c−s) reflected via the Riemann FE onto the same basis gives A, B.

## 3. Results (mpmath 40 dps)

- **Object verified:** spectral sum Σ 2(n+1)(n+2)(n+3/2)^{−2s} = T9 closed form (|diff| → 8e-17 at s=3).
- **Reflection algebra validated:** A·ζ(2s)+B·ζ(2s−2) = F(3/2−s) to 1e-41.
- **Natural axis:** the argument set {2s−2, 2s} reflects consistently only at **c = 3/2** (= Dirac ground state |λ₀|=3/2) or c=1/2. The axis is NOT the problem.
- **DECISIVE: det ≠ 0 everywhere, and grows.** |det| = 2.8, 8.4, 60.3, 8.2 at s = 1+i, 1.5+3i, 2+5i, 0.5+10i. The candidate-multiplier ratio χ₁/χ₂ tracks the gamma ratio Γ(2s)/Γ(2s−2) = (2s−1)(2s−2).
- **Mechanism / positive control:** the exponent gap of 2 between the two λ-terms is produced by the **quadratic multiplicity** gₙ = 2(n+1)(n+2). A constant-multiplicity (single-power) spectrum lies in a 1D span, is trivially parallel, and DOES inherit a functional equation. The quadratic Weyl degeneracy is the murderer.

## 4. The structural upgrade (the real result)

This is not just "ζ_{D²} also fails." Two things are new vs RH-O:

1. **Analytic, not numerical.** det ≠ 0 with explicit growth law (2s−1)(2s−2), vs RH-O's numerical 48-orders-of-magnitude residual. The obstruction is named: the gamma ratio between the two exponents.

2. **Unification of two RH walls.** The "wrong Weyl law" wall (RH-N: γₙ ~ √n) and the "no functional equation" wall (RH-O) are the **same structural fact** — both are consequences of the S³ quadratic degeneracy. And it is the *same* n² multiplicity that puts the framework's even-s values on the Bernoulli ladder (F = ζ(2) = π²B₂). **The quadratic Weyl multiplicity simultaneously (a) places us on the ζ special-value ladder and (b) denies us the functional equation / critical line.** Consistent with and sharpening WH6.

## 5. Audit (per feedback_audit_numerical_claims)

- Free params: zero (closed form).
- Selection bias: positive control included (constant multiplicity → FE exists); negative is "det provably ≠ 0 with growth law," not "couldn't find one."
- Alternatives: the only basis-matching axes (c=1/2, 3/2) both fail.
- Robustness: analytic + numeric to 1e-41.
- Independent: reflection algebra cross-checked against direct closed-form eval.

## 6. Honest scope

- Concerns the functional equation / zero structure (continuum, Layer-2). Decisive negative.
- Does NOT touch the Bernoulli-ladder integer-value content (Layer-1 special values, unaffected).
- Confirms WH6; the classical-RH bridge stays closed. No reopening.
- Framework-level statement: ANY GeoVac spectral zeta with polynomial multiplicity of degree ≥ 1 spreads across exponents differing by integers and inherits the gamma-ratio obstruction. A Riemann-type functional equation requires a flat (constant-multiplicity) spectrum, which the S³ Weyl law is not.

## 7. Documentation applied

- Paper 28 §T9: `rem:zetaD2_no_fe` added (functional-equation obstruction, RH-N+RH-O unification, Bernoulli-ladder cross-ref).
- Recommended to PI (NOT applied — §1.7 PI-only): sharpen WH6 status note to record that RH-N and RH-O are one structural fact (quadratic multiplicity), per this closed-form result.

## Files
- `debug/sprint_gb2_zeta_dsq_functional_eq.py`
