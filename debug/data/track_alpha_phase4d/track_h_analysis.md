# Track alpha-H: Base-Fiber Tensor Trace — Analysis

**Phase:** 4D alpha sprint
**Date:** 2026-04-10
**Goal:** Take the Phase 4B Track alpha-C Fock-weight base data on (n,l) and
tensor it with a *continuum* S^1 fiber via the Hopf bundle. Test whether a
regularized trace reproduces any component of
K = pi*(B + F - Delta) = pi*(42 + pi^2/6 - 1/40) = 137.036064...

## Base data (from alpha-C, Fock weight w = (p^2 + p0^2)^-2 at p0 = 1)

| (n,l) | w(n,l) | 2l+1 | l(l+1) | (2l+1) l(l+1) w |
|-------|--------|------|--------|-----------------|
| (1,0) | 7/16   | 1    | 0      | 0               |
| (2,0) | 5/8    | 1    | 0      | 0               |
| (2,1) | 3/8    | 3    | 2      | 9/4             |
| (3,0) | 5/8    | 1    | 0      | 0               |
| (3,1) | 17/32  | 3    | 2      | 51/16           |
| (3,2) | 11/32  | 5    | 6      | 165/16          |

Sum of Hopf-Casimir weights = **63/4** = 6*B*|kappa| with B = 42, |kappa| = 1/16.
(The task prompt stated 21/8 but the actual sum from the alpha-C JSON is 63/4.
This is consistent: the identity 6*B*|kappa| = 6*42/16 = 252/16 = 63/4.)

Also needed:
- raw w sum = sum_{(n,l)} w(n,l) = 7/16 + 5/8 + 3/8 + 5/8 + 17/32 + 11/32 = **47/16**
- Hopf-weighted sum = sum_{(n,l)} (2l+1) w(n,l) = **49/16**

## Targets

| Quantity            | Value              |
|---------------------|--------------------|
| K = pi*(B+F-Delta)  | 137.036064         |
| K/pi                | 43.619900...       |
| B + F               | 43.644934...       |
| F = pi^2/6          | 1.644934...        |
| B                   | 42                 |
| Delta               | 1/40 = 0.025       |

## Variant results

### Variant (a): Uniform fiber
Every (n,l) sector carries the same S^1 zeta, zeta_{S^1}(1) = pi^2/3 (sum 2/k^2 over k >= 1).

    T_a = (63/4) * (pi^2/3) = 21 pi^2 / 4 = 51.815...

Closest target: B+F = 43.645, relative error **18.7%**.
No match. T_a is pure pi^2, no additive rational, wrong magnitude.

### Variant (b): Scaled fiber, circumference proportional to 2l+1
With L_l = 2l+1, zeta_{fiber}(1; l) = L_l^2 / 12, and
    T_b = sum_{n,l} (2l+1) l(l+1) w (2l+1)^2 / 12 = 409/16 = 25.5625

Closest target: B = 42, relative error 39.1%.
Also tried L_l = l(l+1): T_b = 131/4 = 32.75 (21.9% error vs B).
Also tried L_l = sqrt(l(l+1)) -> same as L_l^2 = l(l+1) formula: 97/16 = 6.0625.
All pure rationals, no pi, no match.

### Variant (c): Fiber-averaged base x continuum S^1 zeta
No Hopf-Casimir weighting; just multiply the raw Fock-weight sum by pi^2/3.

    T_c_raw = (47/16) * (pi^2/3) = 47 pi^2 / 48 = 9.664...
    T_c_hopf = (49/16) * (pi^2/3) = 49 pi^2 / 24 = 20.150...
    T_c_avg = T_c_raw / 6 = 47 pi^2 / 288 = 1.611...

Closest targets:
- T_c_avg = 1.611 vs F = pi^2/6 = 1.645 -> 2.1% error
  (But this is trivial: T_c_avg = (47/288) pi^2 vs (1/6) pi^2 = (48/288) pi^2,
  both are rational multiples of pi^2, they're bound to be close in the ~1%
  region for any sensible rational.)
- No match for B, B+F, K.

### Variant (d): Mellin / heat-kernel representation (THE CRITICAL VARIANT)

The Fock weight (p^2 + p0^2)^-2 has the Schwinger representation
    (p^2 + p0^2)^-2 = (1/Gamma(2)) int_0^inf t exp(-(p^2 + p0^2) t) dt.

Build the full trace at p0 = 1:
    I(1) = int_0^inf t exp(-t) theta_S2(t) theta_S1(t) dt
with heat kernels
    theta_S2(t) = sum_{l >= 0} (2l+1) exp(-l(l+1) t)
    theta_S1(t) = sum_{k in Z} exp(-k^2 t).

**Exact closed form.** Expand theta_S2, swap sums, integrate t exp(-a t) dt = 1/a^2:
    I(1) = sum_{l >= 0} (2l+1) * T(1 + l(l+1))
    T(b) = sum_{k in Z} 1/(b + k^2)^2

Mittag-Leffler differentiation (sum_k 1/(b+k^2) = pi coth(pi sqrt(b))/sqrt(b)):
    T(b) = (pi^2 csch^2(pi sqrt(b))) / (2 b) + (pi coth(pi sqrt(b))) / (2 b^{3/2}).

Numerical evaluation at 50 dps:
    I(1) = 3.93746913838782558297385155871...

**Asymptotic expansion (small-t Weyl analysis).** The S^2 heat kernel expands as
    theta_S2(t) ~ 1/t + 1/3 + t/15 + 4 t^2/315 + t^3/315 + ...
and theta_S1(t) ~ sqrt(pi/t) (Jacobi inversion, the *only* place pi enters from the fiber).
Combining with the t*exp(-t) weight and integrating term by term:

| S^2 coef | t^{e}       | Integral = sqrt(pi) * a_k * Gamma(e + 3/2) |
|----------|-------------|--------------------------------------------|
| 1        | 1/t         | pi                                         |
| 1/3      | 1           | pi/6                                       |
| 1/15     | t           | pi/20                                      |
| 4/315    | t^2         | pi/42                                      |
| 1/315    | t^3         | pi/72                                      |

Sum: **I_asym = pi * (1 + 1/6 + 1/20 + 1/42 + 1/72) = 2119 pi / 1680 = 3.96252073...**

Every asymptotic term is **linear in pi**; the Gamma((2m+1)/2) * sqrt(pi)
identity collapses sqrt(pi) * sqrt(pi) * (rational) to a rational multiple of pi.
**No pi^2 appears at any asymptotic order.**

**Residual.** The residual I(1) - I_asym = -0.02505566... does not PSLQ-identify
with any simple basis {1, pi, pi^2, pi^3, zeta(3), log 2, coth(pi)} at tolerance
1e-25. The residual encodes the *non-asymptotic* exp(-l(l+1) t) contributions,
which are csch^2(pi sqrt(1+l(l+1))) and coth(pi sqrt(1+l(l+1))) terms --
genuinely transcendental functions with no algebraic relation to pi^2/6.

**Product with base sum.**
- I(1) * (63/4) = 62.016... (nowhere near K/pi = 43.62)
- I(1) * (47/16) = 11.566... (nowhere near any target)

**Verdict for variant (d):** NEGATIVE.

The continuum Mellin/heat-kernel trace produces:
(i) a leading asymptotic series purely linear in pi (via sqrt(pi) * sqrt(pi) = pi),
(ii) a residual from csch^2 / coth of pi*sqrt(integer), which is genuinely
    transcendental and does not factor as pi^k * Q.

pi^2 = 2 * zeta_R(2) does NOT emerge because the integrand has weight 2*zeta_R(2 * 1)
* Gamma(1) = pi^2/3 times a GAUSSIAN weight -- the sqrt(pi/t) inversion lowers the
power of pi by half, cancelling the zeta_R(2) structure.

### Variant (e): Resolvent at z = -p0^2 = -1

Sum_{k in Z} 1/(k^2 + 1) = pi coth(pi).
    R_e = (63/4) * pi * coth(pi) = 49.665...

Closest target: B+F = 43.645, relative error 13.8%.
Positive-k version: (63/4) * (pi coth(pi) - 1)/2 = 16.958, no match.

## Summary table

| Variant | Expression            | Numeric  | Closest target | Rel err |
|---------|-----------------------|----------|----------------|---------|
| a       | 21 pi^2 / 4           | 51.815   | B+F = 43.645   | 18.7%   |
| b1      | 409/16                | 25.563   | B = 42         | 39.1%   |
| b2      | 131/4                 | 32.750   | B = 42         | 22.0%   |
| b3      | 97/16                 | 6.062    | F = 1.645      | 268%    |
| c1      | 47 pi^2 / 48          | 9.664    | F = 1.645      | 488%    |
| c2      | 49 pi^2 / 24          | 20.150   | B = 42         | 52.0%   |
| c3      | 47 pi^2 / 288         | 1.611    | F = pi^2/6     | 2.1%    |
| d_asym  | 2119 pi / 1680        | 3.963    | F = 1.645      | 141%    |
| d_exact | (closed form)         | 3.937    | F = 1.645      | 139%    |
| d*63/4  | exact I(1) * 63/4     | 62.015   | B+F = 43.645   | 42.1%   |
| e       | (63/4) pi coth(pi)    | 49.665   | B+F = 43.645   | 13.8%   |

**Cleanest "near-miss":** variant (c3), T_c_avg = 47 pi^2 / 288 vs F = pi^2/6, at
2.1% relative error. This is a trivial numerical coincidence between two
rationals (47/288 vs 48/288), not a structural match.

**No variant lands within 1% of any target.** No variant reproduces the
additive B = 42 structure. No variant produces pi^2 with the right coefficient
next to a rational.

## Verdict

**NEGATIVE.**

All five variants of the base-fiber tensor trace fail to reproduce any
component of K = pi*(B + F - Delta). The critical Mellin/heat-kernel variant
(d), which has the cleanest topological motivation (continuum fiber, Jacobi
inversion introducing sqrt(pi) from theta_S1), yields:

- an asymptotic Weyl series that is **purely linear in pi** (because
  sqrt(pi) * sqrt(pi) * Q = pi * Q, collapsing any sqrt(pi) content);
- an exact residual encoding csch^2(pi*sqrt(1+l(l+1))) and
  coth(pi*sqrt(1+l(l+1))), which are transcendental functions NOT of the form
  pi^k * Q.

**pi^2 = 2*zeta_R(2) does not emerge from the base-fiber trace.** The fiber
contribution through Jacobi inversion lowers the pi power by 1/2 at each
asymptotic order, converting the expected zeta_R(2) factor into a linear-pi
factor. The "F = pi^2/6" target cannot be reached through this construction.

## Structural interpretation

Combined with alpha-E's negative result (Fock weight is m-trivial, discrete
S^1 fibers give only rationals), this closes a second entry point to F via the
Hopf bundle: **neither a discrete nor a continuum S^1 fiber tensored with the
(n,l) base data reproduces F = pi^2/6**.

The B = 42 component from alpha-C is a pure Fock-weight property on the
(n,l) base (it survives as Sigma (2l+1) l(l+1) w = 6B/16 = 63/4, independent
of any fiber structure). The F = pi^2/6 component must enter the alpha
formula from a **different mathematical object** than the Hopf-fiber trace --
most likely the *continuous* zeta_R(2) structure of a higher-dimensional
spectral zeta (e.g. zeta_{S^5}(2), a Dedekind eta at a special point, or a
Eisenstein series) rather than the S^1 fiber.

## Recommendation

Move on from Hopf-fiber-trace attempts to identify F. The next entry point
for F should be:
1. The **full S^5 spectral zeta** (since the Bargmann-Segal lattice in Paper
   24 lives on S^5), evaluated at s = 2 where zeta_R(2) = pi^2/6 appears in
   the Epstein-zeta structure;
2. Or a **Dirichlet L-function** associated to the (n,l) lattice's symmetry
   group, where pi^2/6 appears as L(2, chi_0) = zeta_R(2);
3. Or accept that F is a *calibration* exchange constant (Paper 18) that
   enters via the embedding metric and has no microscopic graph origin
   (shelve the alpha combination rule derivation permanently).

The alpha-C result B = 42 from the Fock weight on the (n,l) base is robust
and does not depend on any of the five variants tested here. The alpha-H
sprint closes the "continuum S^1 fiber" avenue.
