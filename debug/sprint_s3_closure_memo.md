# Sprint S^(3) closure — follow-ons 1+2 of the decomposition memo (2026-06-11)

Charter: close the two named open items of
`debug/sprint_s3_decomposition_memo.md` §5 — (1) the trailing-argument-1
evaluator fix + the four invalidated PSLQs, (2) the global
decomposition-vs-anchor check — then freeze the falsifier and make the
k=3 content paper-grade.

Drivers: `debug/s3_closure_trailing1.py` (evaluator fix + gates + PSLQs),
`debug/s3_closure_anchor.py` (rigorous bracket), `debug/s3_closure_recompute.py`
(corrected canonical value).
Data: `debug/data/s3_closure_{trailing1,anchor,recompute}.json`.

## Headline

**Both stage-1 evaluations of S^(3) were wrong, for the same root cause,
and the rigorous re-evaluation closes both open items at once.**
The true value is

    S^(3) = 31.5716 +- 0.0012   (rigorous bracket [31.57063, 31.57300])

— NOT the stage-1 decomposition figure 30.6154, and NOT the stage-1
Levin anchor 30.2197.  The "~0.4 anchor mismatch" of the stage-1 memo
was actually two independent errors of the same disease:
**mpmath nsum(method='levin') silently mis-converges on summands with
logarithmic factors.**

Three victims of the one disease:
1. The t3(b1,b2,1) evaluator (digamma branch, ln j factor): the four
   table entries t3(2,1,1), t3(2,3,1), t3(4,1,1), t3(4,3,1).  The two
   entries NOT covered by any stuffle witness — t3(2,1,1), t3(4,1,1),
   precisely the b2=1 slowest-decay cases (~ln j / j^2) — were off far
   beyond the ~1e-12 the R2/R6 residuals suggested, poisoning the
   stage-1 decomposition value through coefficients up to 3072.
2. The route-B factorized anchor (phi(m) L(m)^2 with L ~ 4 ln m): Levin
   gave 30.2197, sumem diverged to -1.8e+403.
3. (Same family as the S_min lesson: acceleration trusted outside its
   convergence class, cf. CLAUDE.md §3 grab-bag-PSLQ row.)

## Follow-on 2 (global anchor): CLOSED

`debug/s3_closure_anchor.py` evaluates the factorized anchor with NO
series acceleration anywhere:

- All terms positive => partial sums are a rigorous LOWER bound.
- Explicit inequalities L(m) <= 4(ln m + 1), phi(m) <= 2/(m+1)^2
  (verified numerically on sampled m up to 10^6, derivations in the
  driver docstring) give the rigorous integral tail bound
  tail(N) <= 32[(ln N + 1)^2 + 2(ln N + 1) + 2]/N.

Gates (all PASS):
- GA: incremental A(m) vs the digamma/Hurwitz closed form at
  m in {1, 2, 5, 100, 1000}: <= 4e-58.
- GB: exact-Fraction cube-truncated triple sum at N=300 vs the mpf
  implementation: 2.0e-59.  (This is the memo's named "exact-rational
  partial triple sum at moderate N"; the analytically-controlled tail
  is the integral bound above — no raw Levin anywhere.)
- GC: bound validity on sampled m: PASS at m = 10, 100, 10^4, 10^6.

Brackets (head = partial sum, rigorous):

| N       | lower (head)     | upper (head + tail bound) |
|---------|------------------|---------------------------|
| 10^4    | 31.27965825754   | 31.68500779755            |
| 10^5    | 31.52754770527   | 31.58629943475            |
| 10^6    | 31.56614854032   | 31.57418471230            |
| 4x10^6  | 31.57062945557   | 31.57300467231            |

Every bracket EXCLUDES both 30.6154 and 30.2197 already at N = 10^4.

## Follow-on 1 (trailing-1 evaluator): fix + the two implementation lessons

Fix per the stage-1 memo's own prescription: Abel summation over the odd
harmonic weight,

    t3(b1,b2,1) = sum_k (2k-1)^{-1} A(k),  A(k) = sum_{j>=k} a_j,
    a_j = (2j+1)^{-b2} 2^{-b1} zeta(b1, j+3/2),

with A(1) from one Levin-safe nsum, the Hurwitz forward recurrence
zeta(b1, (j+1)+3/2) = zeta(b1, j+3/2) - (j+3/2)^{-b1} building a_j at
one mpf power per term, and A(k) = A(1) - cumsum (absolute error stays
at ulp level).  Both the inner sum and the outer summand
A(k)/(2k-1) ~ k^{-(b1+b2-1)} are pure algebraic decay (no logs).

Two implementation iterations, both instructive:

**v1 (killed after 2.5 CPU-hours):** correctness gate recomputed the OLD
digamma-branch evaluator for comparison — but that evaluator thrashes by
construction (Levin retrying endlessly on the log-modulated summand is
the disease itself), and the backward-seed ceiling-extension design
recomputed the whole array from a fresh Levin seed per extension.
Lesson: never use the broken method as a reference; gate against a
rigorous bracket instead.

**v2 (caught by the G1 brackets in 6 s/value):** term arrays frozen at
the calling precision violate mpmath nsum's PRECISION CONTRACT — the
Levin transform re-requests terms at elevated internal precision, and
its cancellations amplify the frozen rounding error (observed final
errors 1e-4 .. 1e-11 at 40 dps, worst on the slowest series; every
value landed BELOW its rigorous lower bound, which is impossible for a
correct evaluation of a positive series).  Fix (v3): term arrays cached
PER WORKING PRECISION, rebuilt at the current precision on demand.
v3 smoke test: t3(4,3,1) matches its bracket floor to 26 digits;
t3(2,1,1) = 0.13981227 in bracket.

**G1 bracket diagnosis of the OLD stage-A/B values** (the numbers the
stage-1 sprint actually used; bracket = direct partial sum to J=1e5 +
integral tail bound, both rigorous):

| value | old error (abs) | witness that bounded it |
|-------|-----------------|--------------------------|
| t3(2,1,1) | -2.896e-4 | NONE (b2=1 blind spot) |
| t3(2,3,1) | -4.237e-12 | R2 |
| t3(3,2,1) | -2.118e-12 | R2 |
| t3(4,1,1) | -1.412e-12 | NONE (b2=1 blind spot) |
| t3(4,3,1) | -3.91e-20 | R6 |
| t3(3,4,1) | -5.87e-20 | R6 |

Exact witness decomposition: R2's stage-1 residual 6.4e-12 =
4.24e-12 + 2.12e-12 (its two trailing-1 members); R6's 9.8e-20 =
3.9e-20 + 5.9e-20.  The stuffle residuals were measuring precisely
these evaluator errors; the two b2=1 entries had no witness, and
t3(2,1,1) — the slowest series, hence the worst Levin mis-convergence —
carries coefficient ~3072 in the linear table: 3072 x 2.896e-4 ~ 0.89
of the 0.956 decomposition shift (exact residual check below once the
corrected assembly lands).

## Follow-on 1 results (v3 run, all gates PASS)

- **G1** (rigorous brackets, 40 dps): all six new values in-bracket;
  old-value errors quantified (table above).
- **G2** (stuffle witnesses, 220 dps): R2 residual 5.4e-216 (was
  6.4e-12); R6 residual 7.5e-217 (was 9.8e-20).
- **G3** (two-working-precision self-consistency): 2.6e-216.
- **G4** (the four pending PSLQs — ALL ACCEPT, complete conjectural
  level-2 weight-homogeneous bases):

      t3(2,1,1) = (1/2) Li4(1/2) + (1/48) ln^4 2 + (1/24) pi^2 ln^2 2
                  - (19/5760) pi^4                      [w4; res 2.6e-216]
      t3(4,1,1) = (1/192) pi^4 ln^2 2 - (1/64) pi^2 ln2 z3
                  - (11/322560) pi^6 - (1/2) t2(5,1)
                  - (7/128) z3^2                        [w6; res 1.1e-218]
      t3(3,2,1) = (3/64) pi^2 ln2 z3 - (1/9216) pi^6
                  - (1/2) t2(5,1) - (49/256) z3^2       [w6; res 3.4e-219]
      t3(3,4,1) = (1/768) pi^4 ln2 z3 + (5/128) pi^2 ln2 z5
                  + (61/2903040) pi^8 - (1/256) pi^2 z3^2
                  + (1/2) t2(5,3) - (1/2) t2(7,1)
                  - (217/512) z3 z5                     [w8; res 7.3e-339]

  Dimension counting was exactly right: all four reduce, the w6 pair
  introducing t2(5,1) and the w8 entry both t2(5,3) and t2(7,1) — every
  component a CATALOGUED level-2 object, coefficients +-1/2 on the
  depth-2 generators.

## Canonical value + exact consistency closure

Corrected assembly (`s3_closure_recompute.py`, 200 dps from the fixed
220-dps cache; S_min via its identified closed form, cross-checked
against the Levin anchor at 2.2e-199):

    S^(3) = 31.57256120751202275476192954506898337195758446957719412870...

- INSIDE the rigorous bracket [31.57063, 31.57300]  (anchor CLOSED).
- Shift from the stale figure: 0.95718067869252275476.
- Exact accounting: the stale 30.6154 came from the 130-dps numerics
  run, whose own t3 errors (stored in s3_decomp_numerics.json) differ
  from the 220-dps cache errors — Levin mis-convergence is
  precision-dependent (t3(2,1,1): -4.674e-4 at 130 dps vs -2.896e-4 at
  220 dps).  Sum of coeff x (stale - corrected) over the stale run's
  own values = -0.95718067869252104699 vs needed
  -0.95718067869252275476: agreement to 1.7e-15 (residual = the
  sub-1e-30 diffs not enumerated).  The shift is fully accounted for,
  dominated by t3(2,1,1) at table coefficient +2048.

## Falsifier

`tests/test_s3_decomposition.py` — 13 tests, all passing (91 s with
--slow): chain-collapse Fraction identity + CG cross-check, trailing-1
evaluator brackets + stale-value exclusion (regression guard for the
Levin/log failure mode), the four identifications, R2/R6 witnesses,
anchor bracket containing the canonical value and excluding both stale
figures.

## Honest scope (sprint close, 2026-06-11)

- **Exact / theorem-grade:** chain-collapse factorization (Fraction
  bit-identity, N=25 + CG cross-check N=6); the bracket inequalities
  (derivations in the driver docstrings); stuffle identities R2/R6 as
  algebraic statements.
- **Rigorous-numerical:** the anchor bracket [31.57063, 31.57300]
  (monotone partial sums + integral tail bound, no acceleration); the
  trailing-1 brackets.
- **PSLQ-pinned (complete conjectural level-2 bases, 220/340 dps):**
  the four trailing-1 identifications — conditional on the standard
  level-2 basis conjectures (BBV data mine covers w<=11).
- **Numerical:** the canonical S^(3) value (200 dps assembly from
  Levin-safe evaluators, anchored by the rigorous bracket at 4 digits).
- **OPEN (named, unchanged):** the w10 graded piece
  (-2048 t3(4,3,3) - 1024 t3(4,4,2)) — numerically included in the
  canonical value, motivically unidentified (level-2 w10 basis dim 89;
  parity heuristic predicts depth <= 2); the full symbolic closed form
  of S^(3) modulo w10 (mechanical assembly, not yet executed).
- **Realized-depth verdict at k=3: <= 2** on every identified
  component (w10 predicted <= 2).  Loop tower: depth grows strictly
  slower than loop order; nothing new minted.

## Method lesson (for CLAUDE.md section 3)

mpmath nsum(method='levin') silently mis-converges on summands with
logarithmic factors, and the error is precision-dependent (so two runs
at different dps disagree — a usable tripwire).  Three victims in the
stage-1 S^(3) work: the trailing-1 evaluator, the route-B anchor, and
the published "270 N^-1.31 power law" (really (ln N)^2/N^2 increments;
local log-slope 2 - 2/ln N over the fit window).  Fix pattern: remove
the log by Abel summation / use rigorous brackets with no acceleration.
Secondary lesson (v2->v3): cached term arrays violate nsum's precision
contract — terms must be served at the CURRENT working precision (the
Levin transform re-requests terms at elevated internal precision).

## Status: SPRINT CLOSED

Both follow-ons closed; falsifier frozen; paper edits applied
(Paper 28 sec three_loop + prop:depth_k + abstract + open questions;
Paper 55 k=3 rows).  Remaining named follow-ons: w10 room, full
symbolic assembly modulo w10, S^(4) (still undefined corpus-wide).
