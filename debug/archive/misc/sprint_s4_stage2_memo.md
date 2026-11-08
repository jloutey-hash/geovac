# Sprint S^(4) stage-2 — symbolic depth verdict (2026-06-13)

Charter: symbolic identification of S^(4) (the k=4 rung), targeting the
realized-depth question (prop:depth_k predicts depth <= k-1 = 3) and the
deferred trailing constants — WITHOUT the high-precision evaluation wall
that blocked stage-1 numerics.

Drivers: `debug/s4_stage2_scope.py`, `s4_stage2_stuffle.py`,
`s4_stage2_assembly.py`. Data: `debug/data/s4_stage2_{scope,stuffle,assembly}.json`.

## HEADLINE: realized depth <= 3 for S^(4) — SETTLED

The depth verdict is closed. Structural fact (o-space relation): depth-4
content enters S^(4) ONLY through the C4 core (coefficient +1; Ce, Cm, C
are all <= depth-3). So S^(4)'s depth-4 content == C4's depth-4 content,
and the question is whether C4's depth-4 part reduces to depth <= 3.

**Weights 9, 11, 13 — PROVEN (stuffle, numerics-free, conclusive).**
At each, the quasi-shuffle relation system spans the ENTIRE depth-4 word
space (stuffle rank on D4 = #depth-4 words: 35/35, 84/84, 165/165), so
C4's depth-4 vector — indeed any depth-4 vector — collapses to depth <= 3.
Verified by exact mod-p (Mersenne 2^61-1) streaming-rank linear algebra.
This is the bulk of the depth-4 content (15+13+4 = 32 nonzero C4 coeffs).

**Weights 5, 7 — SETTLED by the low-weight depth filtration.** The
stuffle system is rank-deficient here (w=5 rank 0; w=7 rank 7 < 10) for a
structural reason: depth-4 stuffle relations need depth-2 x depth-2
products, and a depth-2 multiple-t has weight >= 3, so 3+3 > 5 leaves
weight 5 with NO depth-4 relations at all. But weights 5 and 7 are below
the threshold for ANY depth->=3 irreducible: in the level-2 (odd /
multiple-t) algebra these live in, every weight-5 value reduces to single
polylog constants {Li_5(1/2), ln^5 2, pi^2 ln^3 2, lambda(5),
lambda(2)lambda(3), ...} (depth <= 1) and every weight-7 value to depth
<= 2. The k=3 closure already exhibited this one weight down:
t3(2,1,1) (w4, depth 3) = (1/2)Li_4(1/2) + ln^4 2/48 + pi^2 ln^2 2/24
- 19 pi^4/5760 — a depth-<=1 level-2 combination. The v4.7.0 sprint
verified ALL GeoVac w<=7 multiple-t values are verbatim Hoffman App-A
rows. Hence the w=5,7 depth-4 atoms (t4(2,1,1,1) at w5; seven atoms at
w7) reduce to depth <= 2.

The seven w=7 atoms: t4(2,1,1,3), t4(2,1,2,2), t4(2,1,3,1), t4(2,2,1,2),
t4(2,2,2,1), t4(2,3,1,1), t4(4,1,1,1) (C4 coeffs +-32768/+-16384).
The w=5 atom: t4(2,1,1,1) (C4 coeff 32768).

**Verdict: realized depth <= 3 for S^(4)** — prop:depth_k (depth <= k-1)
holds at k=4, settled WITHOUT the heavy CH machinery and WITHOUT the
high-precision evaluation wall. (k=3 needed CH to find the depth-2
collapse; k=4's depth-4 collapse falls out of stuffle rank at the bulk
weights, and of low-weight reducibility at w=5,7.)

## Track 1 — per-weight stuffle closure (s4_stage2_stuffle.py)

Quasi-shuffle engine ports directly from the W10 driver. Per weight 5-9:
the stuffle-only system does NOT individually reduce depth-4 atoms (they
survive as generators) — exactly the k=3 W10 saturation pattern (stuffle
alone is rank-deficient vs the regularized double shuffle). What settles
the verdict is not individual reduction but (a) full D4-span at w>=9 and
(b) low-weight filtration at w=5,7.

## Track 2 — assembly-cancellation (s4_stage2_assembly.py)

The depth verdict via C4's depth-4 reduction, mod-p exact rank. Table:

  w   | depth-4 words | C4 nonzero | stuffle rank(D4) | verdict
  5   | 1             | 1          | 0                | low-weight filtration
  7   | 10            | 7          | 7                | low-weight filtration
  9   | 35            | 15         | 35               | DEPTH<=3 (stuffle, proven)
  11  | 84            | 13         | 84               | DEPTH<=3 (stuffle, proven)
  13  | 165           | 4          | 165              | DEPTH<=3 (stuffle, proven)

Performance note: initial sympy exact-Q rank was intractable at w=11,13
(hundreds-to-thousands of words). Switched to mod-p (2^61-1) streaming
Gaussian elimination with D4-projection during row construction — runs in
seconds at all weights. (Rank over a large prime = rank over Q with
overwhelming probability; standard for MZV rank computations.)

## Honest scope

- **Theorem-grade / numerics-free:** the w>=9 depth-4 collapse (exact
  mod-p rank = full D4 span; cross-checked at TWO primes 2^61-1 and
  2^31-1, agreement at every weight — removes the unlucky-prime caveat).
- **Literature-grounded structural:** the w=5,7 settlement (low-weight
  depth filtration; cited Deligne-Goncharov depth theory + the v4.7.0
  Hoffman App-A verification of all GeoVac w<=7 t-values). NOT an explicit
  in-framework reduction.
- **NOT done (named follow-on):** explicit closed forms of the w=5,7
  depth-4 atoms via regularized double shuffle / CH-at-depth-4 (the
  s3_w10_ident.py machinery extended to m=4, with leading-1
  regularization). Would make the w=5,7 settlement self-contained rather
  than literature-cited. Also NOT done: the full closed form of S^(4)
  (the depth-<=3 survivors at each weight) — the depth VERDICT is settled
  but the explicit value/closed-form is the larger remaining stage-2 goal.
- **Validation caveat (unchanged from stage-1 close):** no high-precision
  canonical S^(4) value exists (3-digit bracket only); the trailing
  constants resist high-precision numerics. The depth verdict here is
  purely algebraic (rank + filtration), so it does NOT depend on that
  wall — its strength is the algebra, not numerics.

## Status: depth verdict CLOSED; explicit closed form OPEN (follow-on)
