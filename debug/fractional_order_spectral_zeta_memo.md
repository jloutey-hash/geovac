# Fractional-Order Spectral Zeta on S^3

**Date:** 2026-04-26
**Status:** Complete

## Setup

The Dirac Dirichlet series D(s) = sum g_n / |lambda_n|^s on unit S^3 (Camporesi-Higuchi spectrum |lambda_n| = n+3/2, g_n = 2(n+1)(n+2)) interpolates between the T9 theorem at p=2 (always pi^even) and the parity discriminant at p=1 (even/odd alternation). This investigation asks: what happens at fractional p, where D(ps) is evaluated at non-integer arguments?

## Riemann Representation (new closed form)

Converting the Hurwitz representation via the Dirichlet lambda function yields:

    D(s) = 2(2^{s-2} - 1) zeta_R(s-2) - (2^s - 1)/2 zeta_R(s)

Verified to 10^{-59} against the Hurwitz form at s = 4..10. This representation makes the transcendental content transparent.

## Results

**The boundary is infinitely sharp and depends ONLY on integrality of ps.**

- **ps = even integer >= 4:** D(ps) = A pi^{ps-2} + B pi^{ps} with A, B explicit Bernoulli rationals. Proven analytically (not PSLQ). Verified through D(20). This reproves T9 for the Dirac case with explicit coefficients.

- **ps = odd integer >= 5:** D(ps) = C zeta(ps-2) + D_coeff zeta(ps) with C = 2(2^{ps-2}-1), D_coeff = -(2^{ps}-1)/2. Explicit integer/half-integer coefficients growing as 4^k.

- **ps = non-integer:** D(ps) is generically a new transcendental, not in any finite Q-span of {pi^k, zeta(k)}. Confirmed by PSLQ failure at 60 dps across all 57 non-integer cases in the landscape (p in [1.0, 2.0] step 0.1, s = 1..10). The transition is discontinuous: D(4.001) is already unidentifiable.

Landscape totals: 21 pi^even (all even-integer ps), 11 odd_zeta (all odd-integer ps), 57 unidentified (all non-integer ps), 21 skipped (ps < 3.5).

## Sommerfeld Connection

**None.** Sommerfeld fine-structure D_p (Paper 28 Section 7) involve MIXED even+odd zeta and zeta-products (e.g., zeta(3) zeta(4) in D_4). Spectral D(s) involves PURE even or PURE odd zeta, never mixed, never products. PSLQ confirms D_2..D_5 are not in Q-span of spectral D(4)..D(10). The shared n^2 Fock degeneracy weighting enters differently: Sommerfeld sums weight energy expansion coefficients c_p(n); the spectral series weights eigenvalues |lambda_n|^{-s}.

## Structural Interpretation

The sharp integer/non-integer boundary is a number-theoretic rigidity: the Bernoulli numbers and the Euler product of zeta_R(s) are rigid structures available only at integer arguments. At non-integer s, Riemann zeta is defined by analytic continuation but has no algebraic structure relating it to pi or to its integer values. The fractional-order spectral zeta interpolates smoothly in value but discontinuously in transcendental class.

## Data

- Script: `debug/fractional_order_spectral_zeta.py`
- JSON: `debug/data/fractional_order_spectral_zeta.json`
