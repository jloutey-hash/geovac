# S1 design: are the moving singular points of (AA|BC)'s bond-angle operator carried by single pieces?

*Design note for a record the PM may register (inventory row 6, singular-points column). It computes no scored
quantity. Controls are in `debug/s1_controls.py`, with data in `debug/data/s1_ctrl_*.json`, collected in
`debug/data/s1_controls.json`. The proposed scored driver is `debug/s1_run.py` (written, not run).*

## Read first: two control findings that change the object

**W0-1. AB's order-17 operator is not minimal.** On both AK primes and under two weight vectors, AB's search
(first null vector at the smallest findable order) gives order 17, degree 48, as in AB. But the gcrd of the
generic sequence's independent annihilators (`intclass_lib.minimal_operator`, the routine IC and AK apply to the
pieces) has **order 14, degree 111**. The order-17 operator is a left multiple of the order-14 one. Order 14 is not
found directly because degree 111 is beyond the row budget at L <= 1000. This is the order-degree trade that
`minimal_operator`'s docstring warns about.

**W0-2. AB's "moving points with non-integer exponents" come from the left multiplier.** In the order-17 operator,
every lc factor with non-integer or irrational exponents other than x* moves when the *weights* change, with zeta_A
held fixed (4 factors on one prime, 8 on the other). AB's mu-shift reused the same weights, so it could not tell
the two kinds of movement apart. The same goes for ±1 carrying irrational exponents. In the minimal operator,
identical on both primes and both weights, the non-apparent points are:

| point | exponents | class |
|:--|:--|:--|
| x* = 125/117 | 0..10, 3/2, 7/2, 9/2 | M (the only non-integer point) |
| ±1 | repeated integers; 10 of 14 log-free | L (logarithms forced) |
| two linear points | -3, 0..12; 13 of 14 log-free | L |
| one linear point | -2, 0..12; 13 of 14 log-free | L |

Every other lc factor is apparent-type and moves with the weights. The exponent sets suggest that the two "-3"
points are AB's mu-moving linear points (exponent -3 there too) and that the "-2" point is AB's mu-fixed
-179401/179400. **No mu shift was computed**, so this is an identification, not a measurement.

**Consequence for the cell.** At A1 there is no mu-moving singular point with a non-integer exponent. The candidates
are two points with integer exponents and a forced logarithm, so their local monodromy is nontrivial but unipotent.
The cell's [ARGUMENT] ("each one with nontrivial local monodromy is singular for some single piece") then applies to
logarithmic points, not to the branch points it implies. The cell also quotes "order 17, degree 48" as the bond-angle
operator. AB's record and the inventory cell may need an erratum on two counts: the minimal order (14), and P3a's
reading of non-integer exponents at moving points. I have not edited either.

## Object

The function is f(x) = (1s_A 1s_A | 1s_B 1s_C) in x = cos(angle C B A), at AB's cases A1 and A2. Each piece is one
of AB's 21 coefficient sequences q_{L,c}, L = 0..1000, built by the frozen engine, unchanged. For every piece, and
for the generic combination sum_c w_c q_{L,c}, the operator is the **minimal operator** (gcrd), computed mod p on
AK's primes 1073741441 and 1073741381.

Each F_p-irreducible factor of a leading coefficient gets one class (`s1_controls.points`):

- **M**: some exponent is not an integer.
- **X**: irregular point.
- **L**: all exponents are integers and a logarithm is forced (`intclass_lib.log_free_count`, linear factors only).
- **P**: integer exponents including a negative one, with no logarithm (a pole).
- **A**: apparent.
- **I**: integer exponents at a factor of degree above 1, where the log test is not run.

A factor is "genuine" if its class is M, X, L or P. Two sets are defined from these classes:

- **WF(zeta_A)**: the lc factors of the full operator that are common to both weight vectors (seeds 2026 and 7).
- **G**: the genuine factors of WF(zeta_A) that are absent from WF(zeta_A + 1/997). G is the object S1 scores.

## Proposed clauses (cases A1 and A2, both primes)

- **S1a (carriage).** Every factor in G has a carrier piece. A monodromy point (class M, X or L) needs a carrier
  with monodromy at that point (M, X or L). A pole (P) needs a carrier of any genuine class. For linear factors, the
  full operator's exponents mod 1 must appear among the carriers' exponents mod 1. Bar: no exceptions in any of the
  four runs.
- **S1b (attribution, reported).** For each factor in G, the list of carrier pieces, and the carried degree per
  piece. The per-piece carried degree must agree between the two primes.
- **S1c (exactness).** U_mov = G, where U_mov is the set of genuine piece factors that are absent from
  WF(zeta_A + 1/997).
- **S1d (non-empty).** G is non-empty in both cases. Expected from W0: two L-points at A1.

## Can it fail?

**S1a is a theorem, given correct operators.** The minimal operator of a sum right-divides the LCLM of the
summands' operators, so every solution of the full operator is a sum of piece solutions. If every piece is analytic
and single-valued at a point, the full operator is too there, and the full operator's local monodromy eigenvalues
are among the pieces'. No Qbar-independence of the constants is needed, because the pieces are rational sequences
(AB P2). **An S1a miss therefore indicts the instrument** (a wrong operator, a mod-p accident or a class misread),
so the proposed gate scores it **VOID**, not STOP. T5 shows the instrument can report a miss.

**The live content is S1b, plus S1c to a lesser degree.** S1c can fail only through isomorphic summands that remove
the non-analytic directions. That is rare, but it does happen with collapse, as T2 shows. A cheaper alternative to
running S1: cite the D-module fact, upgrading the cell's [ARGUMENT] to a cited statement, and replace the cell's
order and points with the order-14 data from W0.

**Disclosed degeneracy.** In AK, all twelve E1 and ln components have the same (R0, Q0) = (4, 28) and the same
trade table, so they probably share one operator up to gauge. A point carried by all twelve shows twelve carriers.

## Gate (proposed)

- **VOID:** any control fails; a piece or full operator is not found at order <= 24; the full order differs between
  weights or primes; x* is not a weight-fixed M point; G's total degree differs between primes; or S1a misses.
- **GO:** S1a, S1c and S1d all hold in all four runs.
- **BORDERLINE:** S1c or S1d fails.
- **STOP:** none. The PM may promote an S1c miss to STOP.

## Controls (measured; no A1/A2 piece operator was built; nothing was compared with the full operator's moving points)

| control | 1073741441 | 1073741381 |
|:--|:--|:--|
| T1: toy pieces t^L, L t^L, (2L+1) s^L. Known points 5/3 and 29/20 read as class M with exponent 1/2 mod 1; the zero of L t^L at 1/3 reads as apparent | pass | pass |
| T2: generic combination of the three under two weights. Order 2 (t^L and L t^L collapse); weight-fixed genuine set {5/3, 29/20}; the weight-moving quadratic is apparent | pass | pass |
| T3: attribution. 5/3 is carried by {t^L, L t^L} and 29/20 by {(2L+1) s^L}; the union of the pieces equals the full operator's set | pass | pass |
| T4: moving point (s -> 3/8). The moving genuine set is {29/20} | pass | pass |
| T5: fire test. u^L is in the full combination but in no piece; S1a reports 29/21 as uncarried | pass | pass |
| T6: log class. For t^{L+1}/(L+1), x = ±1 read as L (x = -1 is on the other branch of the square root) and 5/3 as M | pass | pass |
| N3: NA at a fresh geometry (D1 7/5, D2 11/5, zeta_B 6/5, zeta_C 1). Generic minimal operator of order 8 under both weights, with no weight-moving factor | pass | pass |
| N3: genuine points x* = 85/77 (M), ±1 and -22633/21175 (L). x* is carried by the three exponential pieces only, not by the two E1 pieces or e^{-7/5}; the L points are carried by all six. Union of pieces = full | pass | pass |
| W0: A1 full operator only, as in the table above | order 14/111 | same |

## Cost per prime (measured)

| step | seconds |
|:--|:--|
| A1 coefficients, L = 1000 (two needed per run: zeta_A, and the shift) | 628 (zeta_A cached) |
| NA coefficients, N3 | 312 |
| full minimal operator, one weight vector | ~26 |
| one piece operator (NA pieces 7-9 s; AK's A1 components 10-15 s) | ~10 |
| one run (case, prime) | ~1500 (A2), ~900 (A1) |

All four runs, in parallel pairs: about 45 min wall.

## Caps

Mod p, two primes; F_p factors (cross-prime checks use total degrees); log test on linear factors only; gcrd
minimality not certified; W0 ran A1 only.
