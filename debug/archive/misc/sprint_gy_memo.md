# GY1-GY3: the Gaussian geminal on the full triangle (memo)

**Record:** `docs/predictions/2026-10-23_gaussian_geminal_triangle.md` (committed b84e1231 before computing; not edited).
Scripts `debug/gy_*.py`, data and logs `debug/data/gy_*`, everything exact as strings in `debug/data/gy_pieces.json`
(assembler `debug/gy_11_assemble.py`). Nothing outside `debug/` was touched. Nothing was committed.

**Outcome: GY1 HELD (with an extra exponent), GY2 HELD (its point clause cannot discriminate at gamma = 1), GY3 FAILED
within stated bounds.** Every registered control passed before scoring.

## Object, loop and kite

Y_G(t) = <rho_i rho_j rho_k g(r_ik) g(r_jk)/r_ij>, g = exp(-gamma r^2), gamma = 1, scales (1, 1, t), t on the apex,
conventions of record QG.

- **Loop (derived, `gy_0_kernel.py`).** int d^3l g(l-x1) g(l-x2)/|l-x3|^2 = e^{-s12/(8 gamma)} 4 gamma pi^(3/2)
  Daw(R/sqrt(2 gamma))/R, R^2 = (s13 + s23)/2 - s12/4 (Daw = Dawson's function; s-integral symbolic). Against direct
  float64 3-D quadrature at five random triples: worst 6.7e-12 (gate 1e-8).
- **Kite (route B).** The pair needs no Feynman parameter. With the Coulomb line in position space:
  Phi~(w) = sqrt(pi/8) gamma^(-3/2) e^{-z/(8 gamma)} int_0^1 du sqrt(u(1-u)) M^-3 K_G.
  - M = sqrt(1 + u(1-u) z) and kappa^2 = (u - 1/2)^2 z.
  - K_G = [H(M - i kappa) - H(M + i kappa)]/(2 i kappa), with H(b) = I(b)(1 - M b/gamma).
  - I(b) = sqrt(pi/(2 gamma)) e^{b^2/(2 gamma)} erfc(b/sqrt(2 gamma)).
  - K_G agrees with direct r-quadrature to 1e-40 at six points (including imaginary kappa and kappa = 0).
  - phi_0 = (pi^2/32) e^{1/2} erfc(1/sqrt 2) = 0.16135464128884115015... (hand value matches to 2e-39).
- **Split.** The phi_n carry two constants: phi = pi^(3/2)/(2 sqrt 2) [E Ah(z) + Bh(z)], Ah, Bh in Q[[z]],
  E = I(1) = sqrt(pi/2) e^{1/2} erfc(1/sqrt 2) (ODEs symbolic, identity to 1e-38).

## Routes and controls (all PASS)

| Control | Result |
|:--|:--|
| Route A (Dawson loop) vs route B, three w | 0, 2.6e-23, 2.6e-23 (gate 1e-14) |
| Cauchy radius 1 vs 2 (phi, phi_b; phi_a) | 5e-244, 6e-244; 3e-231. The only finite singularity is z = -4 |
| Positive: exp pair in the same series driver vs stored QG gamma = 3/2 series | 1.6e-76 over 100 coefficients (120-digit run) |
| Exact mod-p split series vs 600-digit floats | 41 phi_a coefficients recognised as rationals, 0 mismatches |
| F_3 map control on the stored Slater gamma = 2 factor | recovers the stored CF6 map, t (x-4)^-4, degree 20; 1/1000 breaks it |
| Newton code on the Slater L9 (gamma = 2) | Fuchsian (height 0); extra: GA2's L_G gives slope -1, c = 1/4 |
| Jordan-block code at x = 0 | F_3 gives [3]; the Slater G4 (gamma = 1) gives [3, 1] |
| Twisted order-one search on the Slater G4 | finds its algebraic factor at c = 0 |
| Negative: 1/200 on the operator, applied to Y_G (both normalisations) | 2.3e-3 to 1.4 and 5.4e-2 to 1.1, against 1e-27 unperturbed |

## The operator

- **Data.** 800 Taylor coefficients of phi, phi_a, phi_b at 600 digits; Ah, Bh exactly mod p (600 terms, 1 prime; 400, 10).
- **L9G** (both split parts of Y_log stacked, 5 primes, CRT): order 9, x-degree 13, integers up to 10^22; leading
  coefficient 64 (x-4)^2 (27661 x^7 + ...). Residual on all 800 float coefficients 1.9e-362 (phi), 1.3e-362 (phi_b);
  1.2e-82 on phi_a's 240 usable coefficients.
- **On Y_G** (32 digits, t = 0.5 to 3.0, six points): 1.2e-27 to 2.4e-27 of the largest monomial (2.7e-22 to
  1.4e-26 of the largest group); quadrature h vs 2h 3e-27 to 2e-28. No source term.
- **Minimality over Q(x).** No operator of order 1 to 8 up to x-degree 100, exact mod p on both parts (nullity 0 mod p
  implies none over Q). Upper bound: the Mellin route from the kite gives order 17, and L9G right-divides it exactly.
- **Assumption.** Killing both parts separately assumes E is irrational (not proved). Independent of it: the
  phi-alone float search reproduces order 9 at degree 13 (1e-27 on Y_G); its order 5-8 hits drift in degree with the
  tolerance, and the order-5 one leaves 0.02-1.6 on Y_G (an artifact of the factorially vanishing phi_a).
- **Split parts alone.** The E-part Y_a has the order-4 annihilator F1 (exact mod p). The b-part alone needs order 9.
- **Kite.** Order 5, z-degree 12 (none of order <= 4 to degree 100); factors 1 (solution z^(1/2)(z+4)^(-3/2)), 2, 2,
  both order-2 factors irregular at infinity; genuine finite point z = -4 only.

## Newton polygon at x = infinity (L9G)

Points (k, deg c_k): (0,8), (1..5,13), (6,12), (7,11), (8,10), (9,9); height 4 (irregular); slope -1 edge over
k = 5..9, characteristic polynomial 27661 C^5 (2C-1)^2 (4C-1)^2. Formal solutions: two exp(t^2/4) (in F1), two
exp(t^2/8) (in F2), five moderate. The kite has e^{-z/4}, e^{-z/8}; the sign flips through the apex transform.

## Composition series (exact pseudo-division, chain-first peeler)

L9G = F2 o F1. Both inclusions are exact, and the top of the chain equals L9G.

- **F1:** order 4, x-degree 5; t = 0 exponents {1, 3, 3, 5}, as for the Slater gamma = 1 G4 = R3 o R1 (F1 is not G4);
  x = 4 exponents {-2, 0, 1, 2}.
- **F2:** order 5, x-degree 18; x = 4 exponents {-9/2, 0, 1, 2, 3}. All other leading factors apparent.
- **Irreducibility (bounded).** Frobenius probes of F and F* at x = 0, 1, 4: no factor to x-degree 220. No order-one
  factor exp(c x) t^f (1 - x/4)^e P (c from the Newton polygon, deg P <= 120) of F or F*, the latter with poles allowed at
  the apparent points. lib4's exponent-bounded split searches do not apply to irregular pieces and were not run.
- **Jordan blocks at x = 0:** L9G [3, 3, 1, 1, 1]; F1 [3, 1]; F2 [3, 1, 1].
- **Y_G needs F2:** F1 annihilates Y_a but not Y_b.

## F_3 map searches (for GY3)

A = t^s (x-4)^-b app^-c sum_{j<3} N_j theta^j, s in {1, -1, ..., -11, 3, 5}, b <= 6, c <= 2, deg N <= 128, into
F1, F1*, F2, F2*, L9G and L9G*: **none found**. Other denominators cannot occur, so the limit is the degree and pole box.

## Against the registration

| Prediction | Verdict | Evidence |
|:--|:--|:--|
| GY1: irregular point at infinity; c = 1/8 | **HELD**; value clause **PARTLY** | The minimal (Q(x)) annihilator has height 4, slope -1, so the cap on non-minimal operators does not bite. c = 1/8 (double) occurs, but the dominant growth is exp(t^2/4) (double), not predicted. It sits in F1, the annihilator of the E-part. |
| GY2: genuine finite points in x > 0 exactly {4}; x = 1 regular; no point moves with gamma; no order-2 Legendre or Fuchsian factor | **HELD** (point clause not discriminating) | Genuine points: x = 4 only (degree-7 factor apparent: exponents 0..7, 9, log-free); x = 1 is an ordinary point (exponents 0..8). Composition factors have orders 4 and 5, both irregular at infinity. **Cap:** at gamma = 1 the Slater thresholds 2 gamma and 1 + gamma also sit at x = 4, and 1 - gamma at 0, so {4} would hold for Slater too. "No point moves with gamma" rests on the derived kernel, where only M can vanish, not on a second gamma. |
| GY3: a factor equivalent to F_3 twisted by t, needed by Y_G | **FAILED** (bounded) | No order-3 composition factor (orders 4 and 5, irreducible by the bounded tests above). No F_3 map into any sub, quotient or adjoint within the box, while the same code recovers the stored Slater map. The 1+3+1 escape in F2 is closed by the order-one searches. |

Unscored reading: the K3 rung's local trace at t = 0 survives (a size-3 Jordan block in each factor), but globally
it is fused with exponential solutions into irreducible irregular factors; the confluence reaches the K3 piece too.
Novelty was not re-searched; with GY3 failed, the two-part FINDING condition is not met.

## Deviations from the record

1. Mod-p series used though the phi_n carry two constants (split derived and verified first; see the assumption on E).
2. Route B is one-dimensional, since the Gaussian pair needs no Feynman parameter.
3. The Y w-grid step was refined from h = 1/16 (unconverged, 1e-8) to 1/128.
4. The exp-pair control was run at 120 digits on 100 coefficients.
5. An N = 1200 mod-p bound and a stacked float search died in an interruption; not rerun (the N = 600 bounds supersede them).
6. Defect fixed: acb square roots of balls straddling the negative axis gave NaN (exact z = -rho node,
   kappa = sqrt(W) sqrt(z)); no reported number depends on it.

**Closing.** At gamma = 1 the Gaussian triangle's minimal annihilator is order 9 = 5 + 4, irregular at infinity with
exp(t^2/4) and exp(t^2/8) growth, and genuinely singular only at t = 2. It has no elliptic piece, and no
F_3-equivalent piece within the searched box.
