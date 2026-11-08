# Increment 3c / the periods question — the ordered ξ integral closes, at **weight 1**

**Date:** 2026-08-12 · branch `work/sparsity-boundary`
**Drivers:** `debug/inc3c_weight_probe.py` (the two new moment families),
`geovac/two_center_eri.py` §"Increment 3c" (the closed form)

## The question, stated properly

The ordered double integral `∫∫ P_τ(ξ_<) Q_τ(ξ_>)` was the last numerical step in
the whole two-center engine. It is also an **iterated integral over a simplex** —
the `ξ_< / ξ_>` ordering *is* the simplex — and that is precisely the shape that
defines a period. So "does it close?" and "where does it sit in the transcendence
hierarchy?" are one question, framed by weight:

| weight | objects |
|---|---|
| 0 | rationals, algebraic |
| 1 | `ln`, `E₁`, γ |
| 2 | Li₂, ζ(2) = π²/6 |

Iterated integrals of weight-1 objects generically land at **weight 2**. That was
the natural expectation, and it would have connected this integral directly to the
M1/M2/M3 period work and Paper 18's ladder.

## Answer: it closes at weight 1

Verified, not argued. The τ = 0 ordered integral assembled in closed form:

| p | closed form | quadrature | rel |
|---|---|---|---|
| 1.5 | 0.01506424613994 | 0.01506424613994 | 0.0 |
| 2.0 | 0.00359621732846 | 0.00359621732846 | 6.0e-16 |
| 2.5 | 0.00093942473885 | 0.00093942473885 | 0.0 |
| 3.0 | 0.00026013637617 | 0.00026013637617 | 6.3e-16 |

Symbolic function content: **`{exp, expint, log}` plus Euler's γ.** No
dilogarithm, no polylog, no ζ(2). Pinned by
`test_ordered_xi_integral_closes_at_weight_one`.

**So the obstruction was assembly, not transcendence.** I had been calling this
"the hard part" in a way that implied a transcendental barrier. There isn't one.

## The move that does the work

Substitute `t = ξ − 1` on the outer integral and split

    ln((t+2)/t) = ln(t+2) − ln(t)

Both halves diverge as ξ → 1 and the divergences cancel. Taken *separately* on
[0,∞) the divergence is **never formed** — which is why no higher-weight object
is generated. This is the same lesson as increments 1c and 2, in a third
costume: a formulation that never creates the singularity beats one that creates
and cancels it.

That leaves exactly two new moment families, both weight 1, validated to ~1e-15:

    log_moment(n, c)       = ∫₀^∞ tⁿ e^{−ct} ln t dt = (n!/c^{n+1})(H_n − γ − ln c)
    log_shift_moment(n,c,s)= ∫₀^∞ tⁿ e^{−ct} ln(t+s) dt

plus the E₁ moments already built for increment 2.

**Where γ comes from, structurally.** `log_moment` is the carrier: ψ(n+1) = −γ +
H_n. It is the same γ Phase 0-e found at the ξ = 1 endpoint — but reached without
ever forming the divergence. And ξ = 1 is the *degenerate* ellipse, i.e. the
internuclear axis itself. So γ enters at the axis, which is an observation-side
(Layer 2) feature, not a skeleton one.

## Scope — read before generalising

Established at **σ = 0**. For σ ≠ 0 the derivatives `d^σ Q_τ/dξ^σ` put poles of
order up to σ at ξ = ±1. Increment 3a's parity fact makes the net exponent there
exactly 0 (regular) — but that is **argued, not verified**, and it is the one
place a higher-weight object could still enter. Given that this session has twice
found a claim verified in an easy corner and wrongly generalised, that gap is
named rather than assumed closed.

Also still assembly-only: general (τ, j, H). The primitives cover them; nobody has
built the loop.

## Tagging obligation — OWED

`feedback_tag_transcendentals` requires every transcendental be classified against
Paper 18 + Paper 34.

- **E₁** — already tagged (Paper 18 §"Level 2", the Stieltjes seed). Discharged by
  citation.
- **ln and γ** — **NOT yet tagged.** They are new to this build. Flagged in the
  v4.77.0 records and still owed. The structural hint above (γ enters at the
  degenerate ellipse = the internuclear axis) is a starting point, not a
  classification.

Do not let an exchange result reach a paper before that is done.

## What this changes

1. **Increment 3 is finishable.** The last numerical step has a closed form at
   σ = 0; the rest is bookkeeping over primitives that exist and are validated.
2. **The periods connection is real but shallower than hoped.** The object *is* an
   iterated integral, but it stays at weight 1, so it does not reach the
   dilogarithm/mixed-Tate territory the M1/M2/M3 work lives in. The seed ladder
   across the four classes (none → {E₁,ln} → {E₁,ln,γ}) is a weight-1 filtration
   throughout.
3. It does **not** revive the quantum-resource case, which QC-1 tested negative
   for independent reasons (contraction is free in qubit terms). Closing this
   integral makes the engine faster and makes vanishing decidable; it does not
   make it cheaper on a device.
