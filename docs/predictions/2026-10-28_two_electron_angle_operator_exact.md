# Pre-registered prediction: the bond-angle ODE of the two-electron three-centre integral, from exact data (AO3)

**Registered:** 2026-10-28, before the operator search is run on any two-electron exact data. PI-directed ("Agreed.
Checkpoint then step 1", step 1 being the climb to the two-electron integral with the exact method; "Sounds good").
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** AO (`2026-10-28_angle_operator.md`, STOP: float guessing masked) and NA
(`2026-10-28_nuclear_attraction_angle_operator.md`, BORDERLINE: the one-electron integral has an order-8 operator;
x* the only finite branch point, smallest non-integer exponent 3/2 = the Poisson power; the stronger claim that every
other singular point is apparent missed). The predictions below are NA's held claims carried up one rung, plus the
contact form.

## Object

(XY|XZ) over 1s Slater functions, X at the origin (zX), H functions on Y and Z (zH), |XY| = D1, |XZ| = D2,
x = cos(gamma); f(x) = sum_L F_L P_L(x).

## Method

**Exact coefficients** (`debug/ao_exact.py`): F_L / (4 zX^3 zH^3) = sum_c q_{L,c} c, rationals times named constants,
built bottom-up as in NA plus one more integration layer (inner antiderivatives Phi_n, outer J_n and K_n recursions,
finite parts at r = 0). Disclosed, measured on the unscored case S1 = (1.8, 2.6, 1.3, 1.0) before registering: 31
constants, all exponentials times 1, E1 or a logarithm (the K base integrals, E1 products and gamma cancel).
**Search** (`debug/na_modguess.py`, unchanged): generic sequence with random weights, mod 2^61 - 1 and
4611686018427387847; ansatz sum a_ij x^j T^i; for each R the smallest Q with nullity >= 1 by bisection, rows at least
columns + 30; L = 0..2000; R ascending from 1 to 40, stopping four orders after the first find.
**FOUND** = same smallest (R_min, Q0) on both primes, nullity 1 at (R_min, Q0) on both.
**Analysis** (`debug/na_modanalyze.py`, unchanged, both primes): exponents at x* (rational) against candidates
-3..30 and the half-integers -5/2..41/2; every irreducible factor of the leading coefficient; Fuchs at infinity.

## Cases

| Case | role | D1 | D2 | zX | zH | x* |
|:--|:--|:--|:--|:--|:--|:--|
| S1 | unscored (AO's scored case; used for controls only) | 1.8 | 2.6 | 1.3 | 1.0 | 125/117 |
| W1 | scored | 1.7 | 2.5 | 1.4 | 1.1 | 457/425 |
| W2 | scored | 1.6 | 2.2 | 1.0 | 1.2 | 185/176 |

## Controls measured before registering

| Control | Measured | Bound |
|:--|:--|:--|
| (a) exact x 4 zX^3 zH^3 vs the arb engine `debug/ao_radial2.py` (n = 90, 600 bits), L = 0..30, S1 | 7.6e-84 | 1e-60 |
| (a') same vs the AO 110-digit data (n = 140), L = 0..30, evaluated at 250 digits | 2.7e-105 | |
| (b) NA positive control (exact order 4) and decoy, unchanged code | as in NA | must reproduce |
| (c) a 2030 x 2000 nullspace mod p | 1.3 s | |
| (d) modular coefficient cost | 21 s to L = 200 (about L^2) | |

## Predictions

- **AO3-0 (controls; run FIRST).** (a) on W1 and W2 at L = 0..30 <= 1e-60; (b) reproduces.
- **AO3-1 (the gate).** FOUND in both scored cases within L <= 2000, R <= 40. Probability 0.45.
- **AO3-2.** x* is the only finite point of the operator with a non-integer exponent (both cases, both primes).
  Probability 0.7.
- **AO3-3 (the contact form).** 7/2 is an exponent at x*, and no smaller non-integer exponent occurs there (no 1/2,
  3/2, 5/2): the leading branch is the cusp-cusp (pi^2/630) g_Y g_Z |Y - Z|^7 term. Probability 0.7.
- **AO3-4.** Infinity is irregular. Probability 0.85.
- **AO3-5.** W1 and W2 give the same (R_min, Q0). Probability 0.55.

**Reported, not scored:** R_min, Q0, the order-degree trade; all exponents everywhere; the constant count per case and
the components' operators (if affordable).

**Gate.** STOP if a control fails or AO3-1 fails in both cases. GO if AO3-1 holds in a live case and AO3-2 to AO3-4
hold there. BORDERLINE otherwise.

## Novelty

As in AO and NA: no bond-angle ODE for three-centre Slater integrals found in the literature.

## Caps

1s functions, two geometries, L <= 2000, R <= 40. NOT FOUND says nothing beyond the budget. The linear independence of
the constants over Q is assumed.

## Disposition (appended 2026-10-28, after scoring; nothing above edited)

**Deviation disclosed.** The first launch of the four scans was stopped during the coefficient stage, before any search
output existed: the per-coefficient cost grew like L^2.2 (17 hours projected), from exact big-fraction powers reduced
mod p. The fix computes those powers in the field (`Field.pw`, `Field.ifact` in `debug/na_exact.py`; uses in
`debug/ao_exact.py`). Checked before relaunching: the mod-p coefficients equal the exact rationals reduced mod p (0
mismatches, W1, L <= 30); exact vs the arb engine unchanged (W1 1.1e-83); NA N1 exact vs numeric 1.3e-93. What is
computed is unchanged; only the speed (51 minutes per case and prime to L = 2000).

**Controls, all held.** AO3-0(a): W1 1.1e-83, W2 3.9e-84 (bound 1e-60), 31 constants in each; (b) the NA control and
decoy reproduce on both primes.

**AO3-1 MISSED in both cases, both primes. Gate: STOP.** For every order R = 1..40, no operator up to the largest
degree the rows allow (Q = 985 at R = 1 down to 45 at R = 40): the nullity is 0 at the top degree, so at every lower
one. Reading: within the ansatz sum a_ij x^j T^i, the bond-angle operator of the two-electron integral has more than
about 1970 unknowns at every order up to 40; the one-electron operator had 288 (order 8, degree 31). AO3-2 to AO3-5
unscored. The integral is D-finite in x, so the operator exists; it lies beyond this budget.
Data: `debug/data/ao3_scan_{W1,W2}_{p}.{json,log}`.

**Reported only: the components** (`debug/ao3_components.py`, W1, L = 0..1000, prime 2^61 - 1, smallest Q at the
smallest R <= 24; `debug/data/ao3_comp_W1.{json,log}`). 26 of the 31 component sequences have an operator; their
(R, Q0) take seven distinct values: (6, 16) [2 log components], (8, 30) [4 exp*E1], (10, 72) [5 exp*E1 and 3 log],
(12, 46) [1 exp], (12, 60) [6 exp*E1], (13, 59) [4 exp*E1], (13, 61) [2 exp]. Five pure-exponential components have no
operator up to R = 24 at this length. Reading: the generic operator is a common left multiple of at least seven
distinct component operators of orders 6-13 (sum 74) plus the five unfound ones; in NA the six components collapsed to
order 8, here they do not, which accounts for AO3-1's miss.
