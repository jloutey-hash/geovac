# Pre-registered prediction: the dictionary bounds the rung for s-type trees (RB)

**Registered 2026-11-06** (PI: "I agree with your suggestions. Let's start working them.", in reply to the
recommended diagnostic sprint on the owed rung bound), before any scored case is computed.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 63 eq:tree (the holomorphic period of a tree of n linked densities is int prod_v K0(a_v x) dx)
and prop:rung's fibre identities (an s-type density weight is a polynomial in c, c^-1 and c d/dc applied to the 1s
weight); `/qa paper_63` FULL #1 (v6.4.1), which withdrew "rung at most the 1s rung" as unproved.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope (no molecular encoding).

## Why

The dictionary puts every exponential-orbital integral in the D-module of its graph's 1s integral, but the rung is
read on a one-parameter fibre at special values of the other parameters, where that membership does not by itself
bound the composition factors. The /qa run left the bound as an observation. This record tests a proof of it for
the case the ladder is defined on.

## The argument (derived before any computation)

On the one-scale slice (t; 1, ..., 1) of a tree of n linked s-type densities, the holomorphic period is
G(t) = int w_1(tx) prod_{v>=2} w_v(x) dx, each weight an Euler-operator image of 2 c^-1/2 K0(x / sqrt c) at c = 1:
a combination, with coefficients polynomial in x^2, of K0(x) and x K1(x). Call x^(2k) K0^a (x K1)^b, a + b = n - 1,
an even section.

1. **Generation.** With theta = x d/dx: theta K0 = -x K1 and theta (x K1) = -x^2 K0, so
   theta(x^(2k) K0^(a+1) (xK1)^(b-1)) = 2k S - (a+1) x^(2k) K0^a (xK1)^b - (b-1) x^(2k+2) K0^(a+2) (xK1)^(b-2).
   Solving for the middle term (a + 1 <= n - 1, never zero) expresses every even section as Q(x^2, theta) K0^(n-1),
   by recursion on b.
2. **Integration by parts.** Against K0(tx), theta_x transposes to -(1 + theta_x), and theta_x K0(tx) =
   theta_t K0(tx); x^2 K0(tx) = t^-2 (theta^2 K0)(tx) = t^-2 theta_t^2 K0(tx). Boundary terms vanish (logarithmic at
   0 times x, exponential at infinity). So int K0(tx) T(x) dx = R(t^(+-1), theta_t) F_n(t), F_n = int K0(tx) K0^(n-1) dx.
3. **The scaled leaf.** (tx) K1(tx) = -theta_t K0(tx), so its weights are D_t-images too.
4. **Conclusion.** G = E F_n with E a differential operator with Laurent-polynomial coefficients; D_t G is a quotient
   of D_t F_n, so its composition factors are among F_n's: the general s-type tree sits at or below rung n. No
   Mellin transform or contour shift is used.

Scope: s-type densities; trees; the one-scale slice and its holomorphic period. Not covered: l > 0 (the momentum
factors carry directions), rings, the real three-dimensional reduction, other slices (distinct fixed scales give a
tensor product of Bessel modules, not a symmetric power), attainment.

## Method (fixed before scoring)

Driver `debug/rung_bound_trees.py` (exact Fraction recursion for R; mpmath quadrature at 40 digits for both sides,
theta_t^m F_n evaluated as int (theta^m K0)(tx) K0^(n-1) dx). Cases: n = 3, sections (k,a,b) = (0,2,0), (0,1,1),
(0,0,2), (1,2,0), (1,1,1); n = 4: (0,2,1), (0,1,2), (0,0,3), (1,2,1); t in {0.5, 1.5, 2.5} for n = 3 and
{0.5, 2, 2.5} for n = 4 (off the singular points t^2 = 4 and t^2 = 1, 9).

## Controls (measured 2026-11-06, before registering)

- Quadrature: F_n at 40 vs 60 digits agree to 1.3e-40..5.7e-39 at all six (n, t).
- Planted wrong recursion (denominator a + 2 for a + 1): every case with b >= 1 misses by 0.25 to 0.80 relative;
  b = 0 cases do not use the recursion and are unaffected.

## Gates

- **RB1 (identity):** every registered case satisfies g_T = R F_n to relative 1e-30 or better. GO if all 27 do.
- **RB2 (no hidden restriction):** at least one case per n has b >= 2 (two Euler images in one product) and at least
  one has k >= 1 (an x^2 weight). Built into the case list.
- STOP if any case misses 1e-30: the argument or its implementation is wrong, and the observation stands.
- On GO: Paper 63 restates the rung bound for s-type trees on the one-scale slice at the holomorphic-period level as
  [SYMBOLIC PROOF] (this argument) + [MEASURED] (this check); everything outside that scope stays module membership.
  A guard is written as its own activity and fire-tested.

## Disposition

**Scored 2026-11-06 (driver sha256 prefix 8c06bdfc9f54d1aa, committed e6c68d5b before scoring): RB1 GO, RB2 met.**
All 27 cases satisfy g_T = R F_n: worst relative 3.8e-38 (n = 4, (xK1)^3, t = 2.5); n = 3 worst 1.0e-39. Cases with
b = 0 (and some k = 1 cases) agree to 0.0 at 40 digits because x^2 K0(tx) = t^-2 (theta^2 K0)(tx) pointwise, so both
sides integrate the same integrand; the integration-by-parts step is exercised by every b >= 1 case (9 of 9 at
k = 0 agree to 1e-38..1e-40), against 0.25..0.80 for the planted recursion. The operator sizes (1 to 6 terms) are
recorded in the driver output. The argument stands as written; Paper 63 restates the bound in the registered scope.
