# Pre-registered prediction: the Gaussian geminal on the full triangle (GY1-GY3)

**Registered:** 2026-10-23, before any computation in this file (PI-directed: "Let's keep knocking
down these open items").
**Prompted by:** record GA2 (`2026-10-22_geminal_shapes.md`), which showed at the level of the one
function that a Gaussian geminal's avatar int K0(tx) K0(x) e^{-x^2} dx has an irregular singular point
at infinity, and Paper 59 `sec:density_ring`, which caps that result: "at the level of the one
function". The LiH engine uses Slater shapes; Gaussian geminals are the F12 standard. This file asks
what the Gaussian shape does to the whole correlated triangle at the fibre.
**Extends:** `2026-10-22_geminal_shapes.md` (GA2), `2026-10-17_ring_extensions.md` (QG),
`2026-10-23_gamma1_mechanism.md` (the kite's singular points as pair thresholds).
**Commit status:** committed on its own, before any computation in this file.

## The question

Y_G(t) = <rho_i rho_j rho_k g(r_ik) g(r_jk) / r_ij> with the Gaussian geminal g(r) = exp(-gamma r^2),
gamma = 1, scales (1, 1, t), t on the apex k, at the fingerprint fibre. Which pieces of the Slater
triangle survive, which leave, and where does the Gaussian put its singular point?

## What I know before computing (disclosed)

- The Slater triangle at (1, 1, t): order nine, pieces 1, 1, 1, 1, 2, 3; the algebraic pieces branched
  at t = 2 (density pair) and 2 gamma (geminal pair); the Legendre piece with z = (t^2 + 1 - gamma^2)/(2t),
  thresholds t = 1 +- gamma; F_3's operator twisted by t, singular at x = t^2 = 0, 4 only, independent of
  gamma. The kite's singular points are minus the squared pair masses (record MG).
- g~(k) = (pi/gamma)^(3/2) exp(-k^2/(4 gamma)). The geminal pair with external momentum w needs no
  Feynman parameter: exp(-l^2/4gamma) exp(-(w-l)^2/4gamma) = exp(-(l - w/2)^2/(2 gamma)) exp(-w^2/(8 gamma)).
  So the kite carries the entire factor exp(-z/(8 gamma)), z = w^2, in place of the algebraic N: no
  branch point at any finite z from the geminal pair, and an essential singularity at z = infinity.
- The loop with two Gaussian lines and the Coulomb line, int d^3l g(l - x1) g(l - x2) / |l - x3|^2, is
  exp(-|x1 - x2|^2/(8 gamma)) times a Gaussian integral against 1/|l - q|^2, which is elementary in a
  Dawson or error function of sqrt(alpha) |q| (to be derived and checked against direct 3-D quadrature,
  not taken from memory).
- Landau reading: a Gaussian line has no mass shell, so every cut through a geminal line gives no
  threshold. The cuts that survive are the density pair alone, z = -4, and the cut giving the K3 piece
  (t = 0, 2). The thresholds that bounded the Legendre piece on the Slater triangle, t = 1 +- gamma, and
  the geminal-pair threshold 2 gamma have gone to infinity, where the essential singularity sits.
- Derived by hand, nothing computed: the apex transform of exp(-z/8gamma) z^(1/2) (z + t^2)^(-3/2)
  produces terms like exp(t^2/8gamma) times incomplete gamma functions of t^2/8gamma, so a companion
  solution of Y_G's operator should grow like exp(t^2/(8 gamma)) = exp(t^2/8) at gamma = 1, up to
  powers of t. GA2's avatar companion grew like exp((t-1)^2/4) with its own normalisation.

## Model

The kite of route B with the density pair as before (M^2 = 1 + u(1-u) z, weight sqrt(u(1-u))) and the
Gaussian pair explicit, the Coulomb line in position space by the convolution theorem of
`debug/ga1_0_kernel.py`; or route A (momentum triangle with the closed-form Gaussian loop). Both routes
required, agreeing to 1e-14 at three w; the loop checked against direct 3-D quadrature at five triples
to 1e-8. Taylor coefficients of phi_G by Cauchy's formula at 600 digits (radius 2, the only finite
singularity being z = -4; to be confirmed by a second radius) and, if the coefficients are rational
multiples of a fixed constant, exactly modulo primes; operator guessing in x = t^2 for
Y_log = x^(3/2) sum binom(-3/2, n) phi_n x^n; exact composition series (`debug/ga3_lib.py`); the K3 map
search of `debug/ga3_10_k3map.py`; annihilation on Y_G at 32 digits by the w-grid route; the Newton
polygon of the accepted operator at x = infinity (slopes from the degrees of the Q_i); a numerical
companion by analytic continuation if cheap.

## Predictions

**GY1 (main: the Gaussian shape leaves the Fuchsian class at the level of Y).** The minimal annihilator
of Y_G in x = t^2 has an irregular singular point at x = infinity: its Newton polygon there has a slope
of positive height, and a formal solution grows like exp(c x) with c > 0. Predicted c = 1/(8 gamma) = 1/8
(moderate confidence on the value; the irregularity is the claim).
- Falsifier: a Fuchsian minimal annihilator, checked on Y_G with its controls passing.

**GY2 (no elliptic piece, no geminal-dependent finite point).** The genuine finite singular points of
Y_G's minimal annihilator in x > 0 are exactly {4}; in particular x = 1, the threshold of a zero-momentum
two-density piece L(t, 1; 0), is regular, and no point moves with gamma. No composition factor of order
two is equivalent to a Legendre operator, and no irreducible Fuchsian factor of order two occurs at all.
- Falsifier: a genuine finite point other than 4; or such an order-two factor.
- Confidence: moderate. This is the Landau reading applied to a line without a mass shell.

**GY3 (the K3 piece survives).** Y_G's annihilator has a composition factor equivalent by an exact map
to F_3's operator twisted by t, which Y_G needs (the right factor below it does not annihilate Y_G).
- Falsifier: no such factor, or Y_G annihilated without it.
- Confidence: moderate. On the Slater triangle the K3 piece's singular points never depended on gamma.

Reading registered with the predictions, not scored: the Gaussian is the confluent member of the
family. As the geminal's mass-type thresholds leave for infinity, the elliptic piece's two singular
points coalesce there into the irregular point, and the price of the geminal changes from an elliptic
period to an exponential period. The K3 rung, set by the densities and the Coulomb line, is unchanged.

## Novelty status

- GY1: PROBABLY KNOWN in form. Gaussian-geminal many-electron integrals reduce to Boys and error
  functions (Persson-Taylor, Theor. Chem. Acc. 1997; Barca-Loos, J. Chem. Phys. 2017, arXiv:1704.08051:
  "known in closed form in a few special cases"); such functions are exponential periods. A hit is a
  CHECK.
- GY2 and GY3: NOT FOUND. Searched 2026-10-23 (web): "Gaussian geminal three-electron integrals closed
  form error function Boys Singer Persson Taylor"; "explicitly correlated Slater geminal three-electron
  integrals elliptic Picard-Fuchs". No source gives a differential-operator or period analysis of a
  three-electron geminal integral of either shape. A hit on both is a **FINDING**: at the fibre the
  geminal's shape decides the price (Slater: one elliptic period beside the K3; Gaussian: an exponential
  period beside the same K3, and no elliptic one). The search is repeated at scoring.

## Controls

- Positive: with the Gaussian pair replaced by the exponential pair the same code reproduces the stored
  gamma = 3/2 series of record QG (`debug/data/ring_ext_QG_*`) coefficient by coefficient; the Gaussian
  loop against direct 3-D quadrature; the F_3 map control recovers the stored map on the slice's
  exponential operator at gamma = 2 before it is run on Y_G.
- Negative: a 1/200 change in one coefficient of the accepted annihilator leaves a residual above 1e-6
  of the largest term on Y_G; a 1/1000 change in the K3 map breaks it.
- The Newton polygon is computed on the operator by the same code applied to the Slater operator L9
  (which must come out Fuchsian: every slope zero) before it is applied to Y_G's.

## Caps stated in advance

One gamma, one placement (t on the apex), one Gaussian shape; nothing is run inside the LiH energy. The
confluence reading is a reading. A STOP for feasibility (no annihilator with controls passing) is
reported and not scored. The Newton-polygon clause of GY1 is conclusive only on a minimal annihilator,
since a left multiple of a Fuchsian operator need not be Fuchsian (record GA2); on a non-minimal one the
irregularity rests on the growth of a companion solution, as in GA2, and this is said so.

## Disposition (appended 2026-10-24, after computing; the predictions above are unchanged)

Computed by a sub-agent from this record; scripts `debug/gy_*.py`, data and logs `debug/data/gy_*`, everything
exact as strings in `debug/data/gy_pieces.json`, memo `debug/sprint_gy_memo.md`. PM recheck on an independent code
path, `debug/gy_pm_check.py` (log `debug/data/gy_pm_check.log`).

**Object.** The Gaussian loop is e^{-s12/(8 gamma)} 4 gamma pi^(3/2) Daw(R/sqrt(2 gamma))/R,
R^2 = (s13 + s23)/2 - s12/4, against direct 3-D quadrature at five triples to 6.7e-12. Routes A and B agree to
2.6e-23. The Taylor coefficients of the kite carry two constants, phi = pi^(3/2)/(2 sqrt 2) [E A(z) + B(z)] with A, B
rational series and E = sqrt(pi/2) e^{1/2} erfc(1/sqrt 2); the split was derived and checked to 1e-38 before use.
Controls: the same series driver with the exponential pair reproduces the stored record-QG series at gamma = 3/2 to
1.6e-76 on 100 coefficients; Cauchy radii 1 and 2 agree to 5e-244; the F_3 map search recovers the stored Slater map
at gamma = 2, and a 1/1000 change breaks it; the Newton-polygon code gives Fuchsian for the Slater L9 and slope -1 for
GA2's L_G. **PM recheck of the object:** Y_G by a position-space reduction written independently (the pair closes on
the apex as exp(-2 gamma |x3 - m|^2 - gamma |x12|^2/2), m = (x1 + x2)/2) agrees with the agent's values at t = 0.8,
1.3, 2.5 to 3.6e-14, the float64 floor (grid self-consistency 1.5e-14).

**Operator.** L9G, order nine, x-degree 13, integer coefficients, from exact series modulo five primes; on 800 float
coefficients 1.9e-362; on Y_G at six t at 32 digits 1.2e-27 to 2.4e-27 of the largest monomial; the 1/200 change
leaves 2.3e-3 to 1.4 (GA1 normalisation) and 5.4e-2 to 1.1 (dominant monomial). Minimal over Q(x): no operator of
order one to eight to x-degree 100, exact mod p on both split parts. The Mellin route from the kite (order five) gives
order 17, which L9G right-divides exactly. Requiring both split parts to be annihilated assumes E irrational (not
proved); independently, the float search on phi alone returns the same order nine at degree 13.

**GY1: HELD.** On the minimal annihilator the Newton polygon at x = infinity has height 4 and an edge of slope -1
with characteristic polynomial 27661 C^5 (2C - 1)^2 (4C - 1)^2 (PM recheck exact), so two formal solutions grow like
exp(t^2/8), as predicted, and two like exp(t^2/4), not predicted. Since the operator is minimal, the cap on
non-minimal operators does not apply. The value clause is PARTLY: c = 1/8 occurs, but the dominant rate is 1/4.
Class: **CHECK**.

**GY2: HELD, and not discriminating as registered.** The only genuine finite singular point in x > 0 is x = 4;
x = 1 is an ordinary point (exponents 0 to 8); x = 4 carries the exponent -3/2; the degree-seven factor of the leading
coefficient is apparent (PM recheck exact at x = 1, 4). There is no factor of order two. **Registration error,
disclosed:** gamma = 1 is the equal-mass value, where the Slater triangle also has its thresholds only at t = 2 and
0 and also has no elliptic piece (records QG and MG). So neither clause of GY2 separates the Gaussian shape from the
Slater shape at this gamma. "No point moves with gamma" rests on the derived kernel, where only M can vanish, not on a
second gamma. Class: **CHECK**, of low weight.

**GY3: FAILED, within bounds.** L9G = F2 o F1 exactly, with F1 of order four (x-degree 5) and F2 of order five
(x-degree 18), both irregular at infinity (PM recheck: F1 right-divides L9G with zero remainder; the quotient is F2
times (5x^2 + 48x + 48)^-5). Irreducibility rests on bounded searches: Frobenius probes of F and F* at x = 0, 1, 4 to
x-degree 220, and order-one factors exp(c x) t^f (1 - x/4)^e P with c from the Newton polygon and deg P <= 120, which
also closes a 1 + 3 + 1 split of F2 around an F_3 piece. By Jordan-Holder the composition factors are fixed up to
equivalence, and neither has order three, so no factor is equivalent to F_3's operator. A direct map search from F_3
into F1, F2, L9G and their adjoints (t-powers -11 to 5, (x - 4)^-6, the apparent factor squared, numerator degree 128)
also finds nothing. Each factor keeps a Jordan block of size three at t = 0, the local trace of the K3 piece, inside an
irreducible irregular factor. F1 is exactly the annihilator of the E-part of Y_G. Class: **MISS**.

**Reading (not scored).** The registered confluence reading held for the elliptic rung and failed for the K3 rung:
with Gaussian geminals the K3 piece does not survive as a Fuchsian factor; its local data at t = 0 survive, fused with
exponential solutions. So the K3 rung of the Slater triangle is not set by the densities and the Coulomb line alone;
the geminal's mass shell enters it.

**Novelty.** With GY3 failed the record's FINDING condition (hits on both GY2 and GY3) is not met; the search was not
repeated.

**Deviations.** Exact mod-p series used although phi_n carries two constants (split derived first); route B is
one-dimensional; the Y grid step refined from 1/16 (unconverged) to 1/128; the exponential-pair control at 120 digits
on 100 coefficients; an N = 1200 bound run and a stacked float search died in a rate-limit interruption and were not
rerun (the N = 600 exact bounds cover them); a NaN defect in complex square roots at the branch cut was found and
fixed, and no reported number depends on it.
