# Cold review of debug/exch_reverse_memo.md (exchange order 9 for a > b)

Date 2026-11-19. Adversarial review; nothing under review edited; no paper/test/doc/git change; nothing SEALED opened.
Scripts `debug/exch_reverse_review_{closed,hp,compare,det,wrong,neart1,diagcf,t99,f0check}.py`; data
`debug/data/exch_reverse_review_*.json`; logs `debug/exch_reverse_review_scratch/`.

## Verdict: SOUND WITH CAVEATS

The mathematics holds on every attack. The caveats concern wording and a bug found in the reference evaluator, not
in the claim.

| attack | result |
|---|---|
| 1. continuation and the i pi bookkeeping | tried, held (by hand and numerically) |
| 2. closed form vs exchange_hp at new a > b points, and near t = 1 | tried, held (36.8 to 59 digits). One apparent break at t = 0.99 traced to a bug in the reference evaluator, `geovac/two_center_eri.radial_norm`, not in the closed form |
| 3. decay verdicts at t < 1 | tried, held (two keys checked by hand, audit table re-read) |
| 4. independence lemma for a > b | tried, held (the lemma is symmetric) |
| 5. determinant gcd at new R | tried, held (2ppi, R = 13/4 and 2/3: gcd constant; exact det nonzero at t = 1/10, 1/3, 1/2, 9/10) |
| 6. the coincidence a = 3b | tried, held |

## 1. Analytic continuation

The holomorphy holds. With a = 1 and b = t, the integrand's modulus is e^{-Re t (r1B + r2B)} times terms that do
not depend on t. The integral therefore converges locally uniformly on Re t > 0, so K is holomorphic there (Morera
plus dominated convergence). The normalisation factor t^3 is entire.

The poles at t = 1 do no harm, and here is why. The nine functions, continued from (1, oo), are holomorphic on
D = {Re t > 0} minus (0, 1], because their cuts lie on (0, 1]. The coefficients are rational in t, with poles only
at 0 and +-1. So F is holomorphic on D, which is connected and contains (1, oo). By the identity theorem F = K on
D. Because K is continuous, K(s) = lim F(s + i eps) for s in (0, 1). The limit is taken at s < 1 and never at
t = 1, so the (t - 1)^-k poles never enter.

The bookkeeping, done by hand:
- Ei_cont(w) = gamma + Log w - Ein(-w) for keys positive at t > 1.
- Ei_cont(w) = gamma + Log(-w) - Ein(-w) = -E1(-w) for keys negative at t > 1.
- LC = Log(1 + t) - Log(t - 1).
- psi_b, limit from above: +i pi - 2 i pi + i pi = 0.
- psi_a: -w = 2(t - 1)R lies in the upper half-plane, so it picks up the same +i pi - 2 i pi + i pi = 0.
- From below, every sign flips and the total is again 0.

Each psi bracket is branch-free by itself. f0's integrand FM1 = c(psi_a - psi_b + g7 - g6) is too.

One slip in the memo: it writes Ei(z) = gamma + Log z + Ein z. The correct form is Ein(-z) with a minus sign. The
slip does not matter, because Ein is entire and only the Log carries the branch.

Can the real branch differ from the continuation? No. Every slip in the bookkeeping is a multiple of i pi times a
real coefficient, so it is purely imaginary. The real branch is the continuation exactly, not just up to a term
that happens to vanish.

Numerical check (`exch_reverse_review_closed.py`). The review's own evaluator works from integer numerator and
denominator polynomials of the frozen Q(t) strings. It shares no code with exch_reverse_core. I evaluated F at
t = s + i eps with the continued functions, for all four cases, s in {1/2, 1/3, 9/10}, R = 7/5. The real part
differs from the real-branch value by O(eps^2): about 1e-50 at eps = 1e-24, 1e-26 at eps = 1e-12. The imaginary
part is O(eps), which is just eps times dK/dt. The lower half-plane gives the same numbers. At t = 2 the continued
form equals the real form. The evaluator reproduces the memo's (2,1) values for 2psig and 3dsig, and its
(3/2,4/5) value for 2ppi, to all 50 printed digits.

## 2. Closed form vs exchange_hp at new points

exchange_hp was run at the review's points. It is the Neumann sum and never sees the closed form.
`exch_reverse_review_compare.py`, `exch_reverse_review_neart1.py`:

| case | (a, b) | t | R | tau_max | abs diff | digits | hp last term |
|---|---|---|---|---|---|---|---|
| 2ppi | (5/2, 1) | 0.4 | 7/5 | 20 | 5.5e-48 | 45.9 | -7.1e-45 |
| 3dsig | (5/2, 1) | 0.4 | 7/5 | 18 | 1.6e-38 | 36.8 | 1.5e-35 |
| 2psig | (5/2, 1) | 0.4 | 2 | 18 | 3.5e-36 | 34.4 | 1.8e-33 |
| 2ppi | (1, 9/10) | 0.9 | 2 | 20 | 1.9e-50 | 48.2 | -6.3e-80 |
| 3dsig | (1, 9/10) | 0.9 | 2 | 18 | 8.1e-45 | 43.1 | 7.1e-65 |
| 3ddel | (1, 9/10) | 0.9 | 7/5 | 20 | 1.6e-50 | 47.9 (closed form at dps 90; 42.5 at dps 60) | 3.0e-81 |
| 2ppi | (1, 97/100) | 0.97 | 3/2 | 14 | 4.3e-60 | 58 | -1.6e-64 |
| 2ppi | (1, 49/50), (1, 199/200), (2, 99/50) at R = 3/4, (100/101, 1) | 0.98, 0.995, 0.99, 1.0101 | 3/2 | 14 | rel <= 1.1e-58 | 58-58.6 | |

Three rows need a comment.
- 3dsig (5/2, 1) and 2psig (5/2, 1) are limited by truncation of the hp tau sum: the difference is below the size
  of the dropped tail. The 2psig tail ratio is about 2.3e-3 per step, so the predicted tail is about 4e-36, which
  is the measured difference.
- Near t = 1 the closed form loses roughly 18 digits to cancelling (t - 1)^-9 poles at t = 0.9, and needs
  dps >= 90 there to show its full agreement. This is a numerical caveat that applies on both sides of t = 1.

**The apparent break at t = 0.99, and its cause.** At (a, b) = (1, 99/100), hp and the closed form differ by a
constant relative 1.76e-18 at R = 1 and at R = 3/2, for 2ppi and for 2psig, at hp dps 50 and 100. The closed form
was checked three ways there: the review's evaluator, the memo's exact-coefficient route, and a fresh exact solve
at (1, 99/100), which is equal to the stored form. All three agree to 1e-72.

The same physical point computed as (2, 99/50) at R = 3/4 agrees with the closed form to 58.6 digits. The fault
is therefore in hp's input handling. Cause: `geovac/two_center_eri.radial_norm` returns
`1/sqrt(sp.nsimplify(s))`, and sympy's `nsimplify(Rational(250000, 970299))` returns the wrong number
2^(251/327) 3^(88/109) ..., with relative error 1.76e-18. So the hydrogenic 1s normalisation at Z = 99/100 is
wrong in the 18th digit. A scan of the Z values used in this memo and the earlier memos (1, 2, 3, 4, 9/2, 4/5, 5,
15/2, ...) finds only Z = 99/100 bad. **No published or memo number is affected**, but it is a latent production
bug, outside this review's scope; raise it to the PM. Fix: drop the `nsimplify`, since s is already an exact
Rational.

## 3. Decay at t < 1

I re-read `exch_reverse_decay.json` (13 keys, 4 cases, denominators t, t +- 1 only) and checked two keys whose Ei
argument changes sign:
- e^{-2R} Ei(2(1 - t)R): at t < 1 the argument is positive and Ei(x) ~ e^x/x, so the rate is
  -2 + 2(1 - t) = -2t < 0. It decays, slowly as t -> 0, but on all of (0, 1).
- e^{-2tR} Ei(2(t - 1)R): the argument is negative, so Ei decays like e^{2(t-1)R}. The rate is -2 < 0.
- Also e^{2R} Ei(-2(1 + t)R): the rate is 2 - 2 - 2t = -2t.

The homogeneous solutions (A + B ln r, A + B r) do not decay, so the decaying solution is unique. The fresh solves
at t = 1/3 report 0 free unknowns, 0 residual failures, and no non-decaying key (2psig, 3ddel JSONs). My fresh solve
at t = 99/100 for 2ppi repeats this: 159 solves, 0 free, 0 failures, equal to the stored form.

## 4. Independence

The proof is symmetric in a <-> b. The log coefficients at R = 0 are the same for either sign of the Ei argument,
because the principal value is gamma + ln|x| + an entire part. Step 3 needs {0, +-2a, +-2b} distinct and
c != 0; Step 4 needs {a + b, 2a, 2b} distinct. All three hold for any a != b. Confirmed.

## 5. Determinant (`exch_reverse_review_det.py`, 2ppi)

At R = 13/4: degree 146, multiplicities (59, 34, 22), residual degree 31, one root in (0, 1), in
(0.4845, 0.4846). At R = 2/3: residual degree 31, no root in (0, 1). The gcd of the two is constant.

I also computed exact rational determinants of the 9 x 9 matrix directly, by `fraction_det` at
t in {1/10, 1/3, 1/2, 9/10} and both R. All are nonzero. The degrees and multiplicities match the memo's table.

## 6. a = 3b

Ei(-(a+b)R) and Ei(2(b-a)R) become the same function. But psi_b uses it multiplied by e^{-2bR}, and the keys
(lambda, mu) stay distinct. In the lemma, g2's log coefficient A_2 = 1 is removed at Step 3 by the distinct rates,
before any Ei identity is used. The fresh solve at (3, 1) is clean, and the direct det at t = 1/3 is nonzero.
Held.

## Required caveats

1. Read the real branch explicitly. The paper's psi bracket "+ ln((b-a)/(b+a))" is the log of a negative number
   when a > b. It must become -ln((a+b)/|b-a|), and Ei of a positive argument must be the principal value.
2. Numerical evaluation near t = 1 needs extra working precision, about 18 digits lost at t = 0.9. This holds for
   either ordering.
3. The single-R determinant now has real zeros on both sides of t = 1 (for example near b/a = 0.4846 at
   R = 13/4). The paper's "real zeros b/a > 1" must become "real zeros b/a != 1".
4. Separate from the claim: the radial_norm/nsimplify bug in the reference evaluator (Sec. 2). Report it; it does
   not touch any reported value.

## Proposed Paper 58 wording (minimal)

- Line 442-443: replace "for every $0<a<b$, where in ... the $l>0$ orbital carries the smaller exponent $a$ (the
  reverse ordering is not covered)." with "for every $a\neq b$, $a,b>0$, with $a$ the exponent of the $l>0$
  orbital; for $a>b$ the same closed forms hold with the real branch (principal-value $\mathrm{Ei}$ and
  $\ln\tfrac{|b-a|}{b+a}$), being the analytic continuation in $b/a$ around $b/a=1$."
- Line 441: "real zeros $b/a>1$" becomes "real zeros $b/a\neq1$".
- Line 421: the psi bracket gets the absolute value inside the logarithm, $\ln\tfrac{|b-a|}{b+a}$.
- Lines 454-456 and 459-460 (m != 0): "for every $0<a<b$ ... with $2p$ or $3d$ carrying the smaller exponent $a$"
  becomes "for every $a\neq b$ ... with $a$ the exponent of the $2p$ or $3d$ orbital".
- Line 542-543 (Caps) and table row 1710: "(in the mixed cases with the $l>0$ orbital carrying $a$)" and "every
  $0<a<b$" become "every $a\neq b$". The 1s, $2p\sigma2p\sigma$ and $2p_{+1}2p_{+1}$ cases are symmetric under
  reflection, so they already cover $a>b$.

## Wrong answers a guard should reject (`exch_reverse_review_wrong.py`; 2ppi and 3dsig at (5/2, 1), R = 7/5)

- LC sign flipped inside psi (+LC instead of -LC): off by a factor of 99 (2ppi) and 1.15e4 (3dsig).
- The a <-> b relabelled value K(R; b, a), the "symmetry" shortcut: relative error 0.84 and 0.46.
- Taking the real part of a mis-branched continuation: NOT rejectable. Any i pi slip is purely imaginary, so a
  guard must check that Im = 0 before taking the real part. Otherwise it accepts an expression whose bookkeeping is
  wrong.
- An hp comparison at Z = 99/100-type inputs fails for the wrong reason until radial_norm is fixed.
- A determinant guard that checks roots only in (1, oo), the old criterion: it would pass a gcd with roots in
  (0, 1).
