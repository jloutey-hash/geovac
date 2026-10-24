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
