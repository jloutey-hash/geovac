# Pre-registered: the minimal bond-angle operator of (AA|BC) has no moving point with a non-integer exponent (record MO)

**Registered:** 2026-11-17, before any scored computation. Owed item: inventory row 6, singular-points column (the
"S1 per-piece test"), re-aimed by its design pass.
**Owner:** Paper 59 `sec:scope`, record AB paragraph, verbatim: "One registered clause failed: the bond-angle operator,
of order 17, has further singular points that move with zeta_A and no 5/2 exponent at x*"; caps "minimal operators not
certified". The reading that the moving points carry non-integer exponents is record AB's (its "New against
Silverstone" paragraph and P3a), not the paper's.
Cold review before commit: 4 MUST-FIX (MO2 read only the base point set, so it added nothing beyond MO1 and could not
fail at A1; GO could be met with no moving point; shift-side consistency checks missing; instrument outcomes, unaccounted
exponents and irregular points, mapped to STOP), 5 SHOULD-FIX (MO1 read one weight; the div control had no code; owner
misquoted; the order-14 finding goes to the PI whatever the gate; the class-I exclusion undisclosed), 3 NIT; all
applied, the driver re-frozen.

## Why the object changed

The S1 design pass (`debug/s1_design_memo.md`) was asked for a per-piece test of the row-6 [ARGUMENT]: every
zeta_A-moving singular point of AB's operator is singular for some single constant piece. Its control W0, run at case
A1 on record AK's two primes, found something earlier:

- AB's own search (`aabc_run`) reproduces order 17, degree 48 under two weight vectors. The greatest common right
  divisor of the generic sequence's annihilators (`intclass_lib.minimal_operator`, the routine IC and AK use) has
  **order 14, degree 111**. AB's search could not reach it: degree 111 needs more rows than L <= 1000 gives.
- In the order-17 operator every factor with non-integer exponents other than x* moves when only the random weights
  change. AB's zeta_A shift reused the same weights, so it could not tell weight-movement from zeta_A-movement.
- In the order-14 operator, at A1 on both primes and both weights, the only finite point with a non-integer exponent
  is x* = 125/117 (exponents mod 1 {0, 1/2}). The other genuine points are +-1 (forced logarithms) and three linear
  points with integer exponents and one forced logarithm: two with lowest exponent -3 and one with -2 (20665/15246).
  Every other factor is apparent and moves with the weights.

So the A1 reading is already seen (it is a control, not a prediction). This record tests it where it is not seen:
case A2, and the zeta_A shift at both cases, under weights that are varied.

The per-piece carriage clause S1a of the design pass is a theorem once the operators are right (the minimal operator
of a sum right-divides the LCLM of the pieces' operators; the pieces are rational sequences). It is kept here as an
instrument check (MO3, a miss is VOID), not as a prediction.

## Object

Record AB's frozen engine (`debug/aabc_exact.py`, `pieces`/`combine`, unchanged): the 21 constant pieces of the
Legendre coefficients of (1s_A 1s_A | 1s_B 1s_C), L <= 1000, at AB's cases A1 (x* = 125/117) and A2 (x* = 89/80),
reduced modulo AK's primes 1073741441 and 1073741381. The generic sequence is a random combination of the pieces
(two weight vectors, seeds 2026 and 7). Each run builds it at zeta_A and at zeta_A + 1/997 (geometry and the other
exponents fixed). "Genuine" means non-apparent and fixed under the change of weights.

## Predictions

The scored run is

    python debug/s1_run.py controls 1073741441;  python debug/s1_run.py controls 1073741381
    python debug/s1_run.py run A1 <p>;  python debug/s1_run.py run A2 <p>   (both primes)
    python debug/s1_run.py summary          -> debug/data/s1_run.json

- **MO0 (instrument).** Failure on any item is VOID:
  - controls T1-T7 and N3 pass on both primes (T7: a toy point that moves with a parameter and has exponent -1/2 is
    reported as moving with a non-integer exponent, so MO2 can fail);
  - the shifted coefficient sequences differ from the base ones (the shift is applied);
  - a minimal operator is found at order <= 24 in every run;
  - its (order, degree) is the same across the two weights and the two primes, per case, at the base zeta_A and at
    the shift;
  - AB's search operator for the same generic sequence is right-divisible by it (pseudo-remainder zero);
  - x* is a weight-fixed point with a non-integer exponent, at the base and at the shift;
  - no scored point has exponents the instrument cannot account for, and no genuine point is irregular (class X).
  The div control (an AB-style search operator is right-divided by the minimal one; an unrelated operator is not) is
  rerun under `controls`.
- **MO1 (A2, fresh; scored).** In A2's minimal operator, the only finite point with a non-integer exponent is
  x* = 89/80, read identically under both weights. Its exponents mod 1 are {0, 1/2}, and its half-integer exponents
  lie in {3/2, 7/2, 9/2}. Probability 0.75. (Reported for A1, not scored: A1's base side is already seen in W0.)
- **MO2 (both cases; scored).** A genuine point moves if it is in the base operator's genuine set and not the shifted
  one's, or the reverse; moved points are read at both positions. Every genuine moving point other than x* and +-1
  has only integer exponents. The live content is the shifted side at both cases and the base side at A2: at A1 the
  base side is the W0 control. n_moving (per case) is reported. Probability 0.7.
- **MO3 (instrument).** Every genuine point of the minimal operator is carried by some piece's operator, with the right
  class: a point with monodromy by a piece with monodromy there, a pole by any genuine piece point, exponents mod 1
  covered on linear factors. A miss is VOID (it would point at the instrument, not the argument).
- **MO4 (reported).** The attribution table (which pieces carry which point), and that each piece's carried degrees
  agree between primes. Record AK found 12 E1/ln components with identical shape; they may share one operator, so a
  point they carry can show twelve carriers.

**Reported, not scored:** the minimal order and degree at each case, the moving points' positions and exponents, the
A1 MO1 check.

## Decision gate

- **VOID** (overrides everything): any MO0 or MO3 item fails, or a step has not finished 48 h after launch.
- **GO:** MO1 and MO2 hold, and n_moving >= 1 in at least one case.
- **STOP:** MO2 misses (some genuine zeta_A-moving point of the minimal operator has a non-integer exponent). AB's
  reading then stands, restated for the minimal operator.
- **BORDERLINE:** otherwise: MO2 holds and MO1 misses, or no genuine point moves in either case (then the paper's "points
  that move with zeta_A" would itself need restating).

## What a GO would change (for the PI; no paper edit is made by this record)

Whatever the gate, W0's finding goes to the PI now: Paper 59's "of order 17" names a non-minimal operator (the
minimal one at A1 has order 14, degree 111, mod p). Beyond that, Paper 59's AB paragraph and the inventory row 6 cell say the operator has order 17 and points that move with zeta_A
with non-integer exponents. On a GO the recommendation to the PI is to restate it: the minimal operator has order 14
(degree 111, two primes, two weights, not certified over Q); its only non-integer exponents are at x*; the extra
points move with zeta_A but have integer exponents and a forced logarithm, so the cell's [ARGUMENT] about moving
points with nontrivial local monodromy then concerns logarithmic (unipotent) monodromy only. That is a change to an
existing claim and goes to the PI.

## Caps

- 1s functions, two geometries, two primes; mod-p readings, nothing certified over Q.
- The constants are assumed linearly independent over Q-bar (as AB), so the pieces' sequences are rational.
- "Moves" is read from one shift of zeta_A by 1/997; a point could move with zeta_A only through a factor that this
  shift leaves fixed modulo the prime (measure zero, not excluded).
- Positions mod p are compared as roots of factors; rational positions are reconstructed only where reported.
- Factors of degree above 1 with integer exponents are classed I (the logarithm test is not run on them) and are left
  out of the genuine set, MO2 and MO3: a moving point of that kind is not scored.

## Controls (measured by the design pass; `debug/data/s1_controls.json`, `debug/data/s1_ctrl_*.json`)

| id | what | result |
|:--|:--|:--|
| T1 | three toy sequences with known points: recovered as non-integer points, exponent 1/2 mod 1; an apparent zero read as apparent | pass, both primes |
| T2 | two toy pieces collapse: generic order 2, a weight-moving apparent factor appears | pass |
| T3 | attribution and union equality | pass |
| T4 | a point that moves when its piece changes is detected | pass |
| T5 | fire test: a planted point carried by no piece is reported uncarried | pass |
| T6 | logarithmic points at +-1 detected | pass |
| T7 | a toy point moving with a parameter, exponent -1/2, is reported moving with a non-integer exponent | pass, both primes |
| div | an AB-style search operator (order 3, degree 7) is right-divided by the minimal one (order 3, degree 4); an unrelated operator is not (`controls` rerun, both primes) | pass |
| N3 | NA (record NA's object) at a fresh geometry: generic order 8 under both weights; x* = 85/77 carried by the 3 exponential pieces only; log points by all six; union equals full | pass |
| W0 | A1 at zeta_A only (the finding above): AB search order 17, degree 48; minimal order 14, degree 111; points as listed | both primes, both weights |

## Frozen (this commit; unchanged at scoring)

| file | sha256 |
|:--|:--|
| `debug/s1_run.py` | `fbaef39ba20360c3b2b7ea5f77346cdc800389dd1bbf8f6ad29dd12b9d6a8296` |
| `debug/s1_controls.py` | `1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687` |
| `debug/aabc_exact.py` | `218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30` |
| `debug/na_exact.py` | `80acd793e8537b0e8a731ce7d0c47df2ed6ac5ab6564d4f7bbc14c0598be1b61` |
| `debug/intclass_lib.py` | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `debug/aabc_run.py` | `6ec038953ebec4dd4ab8d229fb10ce4041678c5bdc49f9429c878eb361d39529` |
| `debug/data/s1_controls.json` | `adf21b1ee456af41cfd3df5871a9a8e92d314b2365d0e0290a3473d23cbbfd8a` |
| `debug/data/s1_run_controls_1073741441.json` | `85ee42b0530b91bea5dffb1579e941a6d2865b5287ee65bbb707ef0d63d1f6d0` |
| `debug/data/s1_run_controls_1073741381.json` | `6c4f2d459aaa0cd0167f5a68e58305c265d5927875a2ff2cfd11855480c5cc74` |

Environment: Python 3.14.0, python-flint 0.9.0, mpmath 1.3.0.

## Phase 2 (after this record is committed and the commit verified)

The commands under Predictions, in parallel pairs (~40 min).

## Results (scored 2026-11-17; record text above unchanged)

The scored run was launched after commit 56f10c90 was verified; the frozen files' sha256 match. Four runs (A1, A2 x two
primes), about 17 min (A1) and 27 min (A2) each.

**Verdict: GO (MO0, MO1, MO2, MO3 held; no VOID flag; MO4 reported).**

- **MO0 held.** Controls T1-T7, div and N3 pass on both primes. In every run the minimal operator has **order 14,
  degree 111**, at both cases, both weights, both primes, base and shift; AB's search gives order 17, degree 48 at both
  cases and is right-divisible by it; x* is a weight-fixed point with a non-integer exponent at base and shift.
- **MO1 held (A2, fresh).** The only finite point with a non-integer exponent is x* = 89/80, identically under both
  weights: exponents 0..10 and 3/2, 7/2, 9/2 (mod 1 {0, 1/2}). A1 reads the same at 125/117.
- **MO2 held.** n_moving = 4, 2 (A1, the two primes) and 2, 4 (A2). Every moving point is linear with exponents
  {-3, 0, 1, ..., 12} (integers, one forced logarithm). The one fixed extra point has lowest exponent -2 (A2: 2509/1980
  on both primes).
- **MO3 held.** Every genuine point is carried by some piece of the right class. The moving points are carried by the
  zeta_A-carrying pieces (the E1, ln and exponential components sharing one exponential prefactor); x* by mu-free and
  mu-carrying exponential pieces.
- **MO4 (reported).** Per-piece carried degrees agree between primes at A1, not at A2.

Readings, not scored:
- Positions are mod-p roots. Where the driver printed a rational it is a reconstruction on one prime: A1's fixed -2 point
  reads 20665/15246 on one prime and 21025/14886 on the other, so it is not a small-height rational (A2's 2509/1980
  agrees on both). The moving points are likewise read mod p.
- The count of moving points differs between primes (4 vs 2). The plausible reading is that the -3 points are roots of
  factors of degree above 1 over Q that split on one prime and not the other; unsplit, they are class I and excluded
  (the cap stated above). The MO4 mismatch at A2 fits the same reading. Neither is measured further here.

Reading: the minimal bond-angle operator of (AA|BC) at both geometries has order 14; its only non-integer exponents are
at x* (3/2, 7/2, 9/2); the points that move with zeta_A are logarithmic points with integer exponents. The non-integer
"moving" points of AB's order-17 operator came from its non-minimal left factor (weight-dependent). For the paper and
the inventory, see "What a GO would change": these go to the PI.
