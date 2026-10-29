# Pre-registered prediction: the bond-angle ODE of the one-electron three-centre integral, from exact data (NA)

**Registered:** 2026-10-28, before the operator search is run on any three-centre data. PI-directed (the stepping stone,
"Yeah, 2 seems fine"; the exact route, "Sure go ahead with that", "Love it").
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** AO (`2026-10-28_angle_operator.md`, STOP), whose disposition and novelty check explain the change of
method: from floating-point coefficients, fake low-order fits matched a known order-4 control to the data's precision
(1e-115 at 200 digits), so no float search could see the true equation; the guessing literature works from exact data.

## Object

V(x) = int phi_X(r) phi_H(|r - Y|) / |r - Z| d^3r, 1s Slater functions on X (zX) and Y (zH), X at the origin,
|Y| = D1, a point charge at Z, |Z| = D2, x = cos(angle Y X Z). V(x) = sum_L F_L P_L(x).

## Method

**Exact coefficients** (`debug/na_exact.py`, built bottom-up from the labels). With rational parameters,
F_L / (2 pi N) = sum_c q_{L,c} c exactly: c runs over constants e^{e0} {1, E1(z), ln v, gamma} (e0, z, v rational);
q_{L,c} rational, computed from the three-term rules of i_L, k_L, exact integration by parts, and a Hadamard finite part
on the piece at r = 0. Computed mod p for the search. On both test geometries the set is six constants (four
exponentials and two exponential times E1); the gamma and log terms cancel.

**Sequences searched.** Per case: the generic sequence S_L = sum_c w_c q_{L,c} mod p with random weights w_c (the
integral with generic constants; it carries the operator of V if the six constants are linearly independent over Q,
which is assumed and not tested), and each component q_{L,c} alone (reported only).

**Search** (`debug/na_modguess.py`). Ansatz sum_{i<=R, j<=Q} a_ij x^j T^i, T = (1 - x^2) d/dx, local on Legendre
coefficients, all mod p; exact null space (flint nmod_mat). For R = 1..12: the smallest Q with nullity >= 1, rows at
least columns + 30. L = 0..1000. Two primes, 2^61 - 1 and 4611686018427387847. **FOUND** = the same smallest
(R_min, Q0) on both primes, with R_min the smallest order found, and at (R_min, Q0) nullity 1 on both.

**Analysis** (`debug/na_modanalyze.py`), of the null vector at (R_min, Q0), on both primes: standard form; for each
irreducible factor f of the leading coefficient over F_p, the orders of the c_k at its roots, regularity, and the
non-negative integer exponents (indicial polynomial in F_p[t]/f); at the rational point x* the indicial polynomial
tested at candidate exponents; at infinity the Fuchs condition.

## Cases

| Case | role | D1 | D2 | zX | zH | x* |
|:--|:--|:--|:--|:--|:--|:--|
| probe | unscored | 2 | 1.5811388300841898 | 1 | 1 | |
| nb | unscored | 1.85 | 2.65 | 1.3 | 1.0 | 2089/1961 |
| N1 | scored | 1.8 | 2.6 | 1.3 | 1.0 | 125/117 |
| N2 | scored | 1.5 | 2.4 | 1.0 | 1.2 | 89/80 |

## Controls measured before registering

| Control | Measured | Bound |
|:--|:--|:--|
| (a) exact decomposition x 2 pi N vs the numeric engine `debug/ao_na.py` (n = 100, 660 bits), L = 0..40, on nb | 6.8e-89 | 1e-60 |
| (a') same on an N1-like run, L = 0..12 (n = 80, 500 bits) | 1.8e-75 (L = 1, the numeric engine's weak spot) | |
| (b) exact positive control f = sum_{k=1..4} x^k e^{c_k x}, c_k = k/3 + 1/2 (known minimal order 4), L = 0..300, generic weights, both primes | order 1-3: nothing up to the row limit; order 4 at Q = 15 (nullity 1); order 5 at Q = 12, order 6 at Q = 11 (nullity 2); identical on both primes | must reproduce |
| (c) decoy: (b) with one coefficient (L = 150) changed by 1 | nullity 0 at (4, 15); no order-4 operator at any Q up to 60 | must reproduce |
| (d) analysis on (b) at (4, 15) | x = 0 exponents 1, 2, 3, 4; x = +-1: 0, 1, 2, 3; x = -3 (multiplicity 3, apparent): 0, 2, 4, 6; infinity irregular | as found |

The floating-point route of AO, run on the same order-4 control at 200 digits and 501 coefficients, found nothing:
its spectrum fell smoothly to 1e-120 with no gap above 4.9 decades. That is the reason for the change of method.

## Predictions (derived before computing; none uses three-centre search output)

- **NA0 (controls; run FIRST).** (a) on N1 and N2 against `ao_na.py` at L = 0..40: <= 1e-60; (b), (c) reproduce.
- **NA1 (the gate).** FOUND in both scored cases within R <= 12, L <= 1000. Probability 0.8.
- **NA2 (the Landau map).** x* is a root of the leading coefficient, and every other irreducible factor of the leading
  coefficient (x - 1 and x + 1 included) is regular with R_min distinct non-negative integer exponents (apparent
  type). Falsifier: any other factor irregular or without R_min integer exponents. Probability 0.65.
- **NA3 (the size of the cusp term).** 3/2 is an exponent at x*. Reason: Poisson's equation turns the density's
  |r - Y| cusp term into a |Z - Y|^3 term of the potential, (2 D1 D2 (x* - x))^(3/2); the coefficients then fall like
  rho^L L^(-2). Bouferguene and Jones (JCP 109, 5718, Eq. 4.24, as read in the WA novelty check) would give
  (l + 1/2)^(-3) for 1s, i.e. 5/2; if 5/2 is an exponent and 3/2 is not, this derivation fails. Probability 0.8.
- **NA4.** Infinity is irregular. Probability 0.85.
- **NA5.** N1 and N2 give the same (R_min, Q0). Probability 0.6.

**Reported, not scored:** R_min, Q0, the order-degree trade (smallest Q at each R); every exponent at x*; the operators
of the separate components and whether their orders add up to R_min.

**Gate.** STOP if a control fails or NA1 fails in both cases. GO if NA1 holds in a live case and NA2-NA4 hold there.
BORDERLINE otherwise.

## Novelty

Per the AO novelty check: no ODE in the bond angle for three-centre Slater integrals, one- or two-electron, was found;
guessing from exact data is standard (Bluemlein-Kauers-Klein-Schneider 2009; Bostan et al. 2009). The exact split of
each Legendre coefficient into a finite constant basis is the Barnett-Coulson / Lowdin multipole route evaluated in
closed form (prior art for the route; the six-constant count is ours, not searched).

## Caps

1s functions, two geometries, R <= 12, L <= 1000. The linear independence of the six constants over Q is assumed.
NOT FOUND says nothing beyond the budget.
