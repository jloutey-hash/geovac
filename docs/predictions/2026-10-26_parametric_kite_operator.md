# Pre-registered prediction: the kite's operator with a parameter left free, and a blind second point (PK1-PK8)

**Registered:** 2026-10-26, before any computation that scores a prediction in this file. PI direction ("Let's do both with sub agents"), on the PM's proposal after
record OE: the derivation of the operator is now the larger part of the cost of the operator route, and it is
redone for every parameter set. Can the operator be had with a parameter left free, or predicted at a new point?
**Prompted by:** `2026-10-26_operator_evaluator.md` (OE1-OE4 and its disposition).
**Extends:** `2026-10-24_ibp_kite.md` (IB4), `2026-10-26_operator_evaluator.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

The off-slice kite phi(z), z = w^2, at (a, b, g1, g2) = (1, 2, 3/2, 3), has an exact order-11 operator over Q(z)
(`debug/data/ib_pf_op_Q.json`). It was built from the 18-master differential system dM/dz = A(z) M, which the
reducer `debug/ib_kite_de.py` samples at 80 random z modulo a prime and rebuilds as rational functions of z. Three
primes were run, the system was lifted to Q from two and checked on the third, and the operator was built from it
by the cyclic vector. Every new (a, b, g1, g2) needs all of this again.

Two questions:

1. Is the operator's coarse structure (singular points, local exponents, order, degree) predictable from
   (a, b, g1, g2) by a stated rule? If so, it can be checked blind at a new point.
2. What is the cheapest way to get the system, and so the operator, with a parameter left free, so that a new value
   of that parameter costs a specialization and a cyclic vector instead of a derivation?

## What I know before computing (disclosed)

Facts read from existing data or from cheap probes run for this record. The probes are listed so that none of them
can be mistaken for a scored result. Scripts `debug/pko_*.py`, logs `debug/data/pko_*.log`.

- **Base-point operator, leading coefficient** (`ib_pf_op_Q.json`, factored over Z): z^5, then (4z+1), (z+4),
  (4z+49), (z+16) each to the power 1, then (z+1), (z+9), (4z+9), (4z+81) each to the power 3, then one irreducible
  factor of degree 92. Total 5 + 4 + 12 + 92 = 113. The d/dz-form degrees are 102, 103, ..., 113.
- **Base-point local exponents** (`debug/ib_exponents_Q.py`, exact over Q):
  - z = 0: -4, -1/2, 0, 1/2, 1, 3/2, 2, 5/2, 3, 4, 5.
  - z = -1/4, -4, -49/4, -16: -1/2 and 0, 1, ..., 9.
  - z = -1, -9: -1, 0, 0, 1, 1, 2, 3, ..., 7.
  - z = -9/4, -81/4: -2, -1, 0, 0, 1, 2, ..., 7.
  - z = -25/4, -25: 0, 1, ..., 10, and the leading coefficient is nonzero there (ordinary points).
  - z = infinity (`debug/pko_exponents_inf.py`): solutions behave like z^rho with rho = -13/2, -6, -11/2 (twice),
    -5, -9/2 (twice), -4, -7/2, -3, -1.
  Whether the repeated and integer-spaced exponents carry logarithms was not checked.
- **Base-point system** (`ib_pf_p5_op.json`, entries of c and A mod p5; `debug/pko_probe_system_degrees.py`): 108
  nonzero entries; numerator degree in z at most 9, denominator at most 10, total at most 19. The least common
  multiple of all denominators has degree 15, with roots z = 0 (power 3), -9/4 and -81/4 (power 2), and -1/4, -1,
  -4, -9, -49/4, -16, -25/4, -10 (power 1). So the system has two roots the operator does not: -25/4 and -10.
  At the base point -25/4 is both -(a+g1)^2 and -(b^2+g1^2), and -10 is -(a^2+g2^2).
- **Coincidences at the base point.** (a-g1)^2 = (b-g1)^2 = 1/4 and (b-g2)^2 = (a-b)^2 = 1, so the pseudo-thresholds
  of the adjacent pairings sit on top of genuine points, and the base point cannot tell whether they are genuine.
  Also a + b = g2, b = 2a and g2 = 2 g1. Only the adjacent thresholds -25/4 and -25 were seen to be ordinary.
- **One reduction sample, timed** (`debug/pko_probe_sample_time.py`, prime p5, two probes run at the same time on a
  machine where another job was running): 27.8 s processor time, single process. Almost all of it is seeding (2016
  seeds, the Laporta elimination); the reduction of the target and the closure under d/dz are negligible. 18
  masters.
- **The same sample with g2^2 replaced by a random residue mod p5** (a, b, g1 unchanged): 27.6 s, 18 masters. So the
  reducer already accepts a generic parameter value. Only g2^2 enters the integration-by-parts identities, not g2.
- **The logged derivation times are pool wall times.** `ib_kite_de.py` samples with `Pool(6)` (unchanged since
  the runs, per git). The logs of p5 and p13 end 12 minutes apart, while each run took 921 to 957 s, so prime runs
  overlapped. 80 samples in about 943 s on six workers is about 71 s per sample under that load, against 27.8 s
  measured here. So "2836 s" in record OE is wall time on six workers per run, not time on one process; see Risks.
- **Generic point chosen by arithmetic only** (`debug/pko_choose_point.py`): P2 = (a, b, g1, g2) = (1, 3/2, 7/3, 8/3).
  At P2 all candidate values below are distinct and nonzero, no mass equals another, and none equals a sum or a
  difference of two others. Among points with a = 1 and the others in (1/2)Z or (1/3)Z up to 4, it has the second
  smallest largest candidate (25), and it keeps the base ordering a < b < g1 < g2.
- Nothing about the operator or the system at P2, or at any other parameter value, has been computed. No direct
  series at P2 exists.

## Model

**The cut rule (from Paper 59 sec. density_ring, made explicit).** Route w in at the apex and out at the centre.
The lines are 1 = density i (k, mass a, index 3/2), 2 = density j (w-k, b, 3/2), 3 = geminal i-k (k-q, g1, 2),
4 = geminal j-k (w-k+q, g2, 2), 5 = Coulomb (q, 0, 1). The cuts that separate w are the line sets whose momenta add
to w: {1,2}, {3,4}, {1,4,5}, {2,3,5}. Each gives a threshold z = -(sum of masses)^2 and a pseudo-threshold
z = -(difference)^2; the massless Coulomb line adds nothing. So the genuine finite nonzero singular points are

  -(a+b)^2, -(a-b)^2, -(g1+g2)^2, -(g1-g2)^2, -(a+g2)^2, -(a-g2)^2, -(b+g1)^2, -(b-g1)^2.

The adjacent pairings {1,3} and {2,4} are not cuts of w (their momenta add to q or w+q), so -(a+-g1)^2 and
-(b+-g2)^2 are not singular. This is the standard Landau analysis of normal thresholds and pseudo-thresholds of a
two-point graph. At the base point it reproduces the eight points found.

**The exponent rule.** At a threshold or pseudo-threshold of an n-line cut with line indices nu_i in D = 3, one
solution has the exponent beta = (n-1)(D-1)/2 - 1 - sum_i (nu_i - 1), the phase-space power near threshold with
each index raised above 1 lowering it by one. The other exponents are the integers of an analytic block. This
gives beta = -1 for {1,2}, -2 for {3,4}, -1/2 for the three-line cuts, as at the base point. The rule does not
involve the masses, so the full exponent sets should carry over from the base point cut by cut. The exponents at
z = 0 and at infinity are set by power counting of regions (small and large w), also mass-free.

**Order and apparent points.** The order is the size of the cyclic subspace of the kite's vector in the 18-master
system, generic in the masses. If the exponents above hold at every genuine point and at infinity, and every
apparent point stays simple (exponents 0, 1, ..., 9, 11), Fuchs's relation fixes the number of apparent points, so
the leading coefficient keeps degree 113 with an apparent factor of degree 92.

**Routes to a parametric operator.** Interpolate the system, not the operator: the system's entries have total
degree at most 19 in z, while the operator has degree 113 and 565-bit coefficients. A new parameter value then
costs a specialization of the system and the cyclic vector over Q(z), the step that took 12.8 s in record OE.

- (a) Reduction over F_p(s)(z) with s = g2^2 symbolic. Needs the reducer's arithmetic rewritten from residues to
  rational functions in s, with intermediate expression swell in the elimination of 2016 seeds. Not chosen. The
  standard tools sample over finite fields instead, as in (b).
- (b) Sample s mod p as well as z, and reconstruct each entry as a rational function of (z, s). Needs no change to
  the reducer (probe above), only a new driver. **Chosen.** Plan: at each of N_s random values of s, 36 random z
  (26 to fit, 10 held out, enough for total z-degree up to 22); normalize each entry's denominator in z; interpolate
  every coefficient in s by a univariate rational reconstruction with 3 held-out s values; lift to Q from two
  primes; check on a third prime at 3 s values x 5 z. Samples: 2 x 36 x N_s + 15.
- (c) A differential equation in the parameter. Raising line 4's index gives dM/ds = B(z, s) M at the same
  samples. But B does not determine A at a new s without knowing A or B as functions of s, which is route (b)
  again. It evolves values, not operators. Not chosen.

**Which parameter first.** s = g2^2 with (a, b, g1) = (1, 2, 3/2) fixed. The family passes through the base point,
so the three saved base systems (p5, p10, p13) are an exact specialization control at s = 9. g2 enters four of the
eight cuts ({3,4} and {1,4,5}), so the s-dependence is not trivial. Caveat: along this family (a-g1)^2 = (b-g1)^2
and (a-b)^2 = 1 for every g2, so the family cannot test the adjacent pseudo-thresholds; P2 does that.

**Cost model.** By homogeneity the entries are rational in (z, a^2, b^2, g1^2, g2^2) of fixed weight, so their
degree in s at fixed z is at most their total weighted degree, which for most entries is the z-degree. The
g2-dependent denominator factors seen at the base point are (z+(g1+-g2)^2) squared, (z+(a+-g2)^2) and (z+a^2+g2^2),
of total degree 7 in s. Estimate d_s (numerator plus denominator degree in s) at most 19, so N_s = d_s + 2 + 3 at
most 24, and the family costs at most 2 x 36 x 24 + 15 = 1743 samples, about 1250 if d_s is near 12. At 27.8 s per
sample alone and about 70 s per sample under load on eight workers, that is about 3.5 to 4.5 h wall on eight
workers. A single point with the existing pipeline is 3 x 80 = 240 samples; trimmed to 2 x 36 + 5 = 77 samples.
So the family pays after about 5 to 7 new values of g2 against the existing pipeline, or about 16 to 23 against a
trimmed one. A family in all four masses (three ratios, by homogeneity) multiplies the count by about d^2 more and
is out of reach of this Python reducer.

## Predictions

Stage 1 is the blind test at P2 = (1, 3/2, 7/3, 8/3), with the unchanged pipeline (`ib_kite_de.py`, 80 samples per
prime, Pool(6)), lifted from two primes and checked on a third. It is scored before any family work.
Predicted values at P2:

| cut | -(sum)^2 | -(difference)^2 | predicted exponents |
|:--|:--|:--|:--|
| {1,2} densities | -25/4 | -1/4 | -1, 0, 0, 1, 1, 2, ..., 7 |
| {3,4} geminals | -25 | -1/9 | -2, -1, 0, 0, 1, 2, ..., 7 |
| {1,4,5} | -121/9 | -25/9 | -1/2, 0, 1, ..., 9 |
| {2,3,5} | -529/36 | -25/36 | -1/2, 0, 1, ..., 9 |
| adjacent {1,3}, not a cut | -100/9 | -16/9 | not singular |
| adjacent {2,4}, not a cut | -625/36 | -49/36 | not singular |

**PK1 (main: the cut rule, blind).** At P2 the operator's genuine singular points other than 0 and infinity are
exactly the eight values in the first four rows. The four adjacent values are not genuine singular points. Every
other root of the leading coefficient is apparent.
- Falsifier: a predicted point is not genuine, an adjacent value is genuine, or any other genuine point appears
  (genuineness by `debug/ib_apparent.py` and exact exponents).
- Strong form, reported but not scored: the adjacent values are not roots of the leading coefficient at all.
- Confidence: moderate to high. The rule is standard; the risk is an anomalous threshold the base point hid.

**PK2 (size).** At P2 the reduction closes on 18 masters, the operator has order 11, its d/dz-form degrees are
102, ..., 113, and its leading coefficient is z^5 times the four three-line factors to the power 1, the four
two-line factors to the power 3, and apparent factors of total degree 92.
- Falsifier: any of these differs.
- Confidence: high for 18 masters (the g2 probe above found 18 at a random g2^2; at P2 all four masses differ).
  Moderate to high for order 11. Moderate for the degrees, which rest on PK3 and on simple apparent points.

**PK3 (exponents).** At P2 the exponents at each of the eight points are those in the table, cut by cut the same as
at the base point. At z = 0 they are -4, -1/2, 0, 1/2, 1, 3/2, 2, 5/2, 3, 4, 5. At infinity they are -13/2, -6,
-11/2 (twice), -5, -9/2 (twice), -4, -7/2, -3, -1.
- Falsifier: any exponent set differs.
- Confidence: moderate. The singular exponents follow from the rule; the rest of each set is carried over, not
  derived.

**PK4 (the extra points of the system, exploratory).** At P2 the least common multiple of the system's denominators
has z-roots 0 (power 3), -25 and -1/9 (power 2), the other six predicted points (power 1), and -(a^2+g2^2) = -73/9
and -(b^2+g1^2) = -277/36 (power 1), and nothing else. In particular no adjacent value.
- Falsifier: any other root, a missing one, or a different power.
- Confidence: low. This reads the base point's -25/4 as b^2+g1^2, not as (a+g1)^2, from one coincident example.

Stage 2 is the family in s = g2^2 at (a, b, g1) = (1, 2, 3/2). It starts only after the stage-1 disposition is
appended and the identification control of stage 1 has passed.

**PK5 (degree in the parameter).** At one fixed random z, every entry of c and A is a rational function of s whose
numerator and denominator degrees add to at most 19. Measured on one prime from 30 values of s plus 5 held out
(35 samples), before the full family is run.
- Falsifier: any entry above 19, or no reconstruction with 30 values.
- Confidence: moderate.

**PK6 (the family is exact and specializes).** The family system over Q(z, s), lifted from two primes and checked
on a third, (i) at s = 9 equals the saved base systems at p5, p10 and p13 entry by entry, exactly; (ii) at the
held-out g2 = 7/3 (s = 49/9, not used in the reconstruction) agrees with direct reductions at 5 fresh z; (iii) at
g2 = 7/3 gives an operator whose recurrence reproduces the direct Taylor coefficients 6 to 40 from 0 to 5, to 1e-40
relative. g2 = 7/3 makes no new coincidence among the cut values (11 distinct of 12, the 12th the family's own).
- Falsifier: any mismatch in (i) or (ii), or worse than 1e-40 in (iii).
- Confidence: moderate to high, if PK5 holds.

**PK7 (cut rule along the family).** The operator at g2 = 7/3 from the family has its genuine points at the cut
values of (1, 2, 3/2, 7/3) by the same rule: -9, -1, -529/36, -25/36, -100/9, -16/9, -49/4, -1/4.
  The adjacent values -25/4 and -169/9 and -1/9 are not genuine (the fourth, -1/4, coincides with a genuine one).
- Falsifier: as PK1.
- Confidence: as PK1, and it is not independent of it.

**PK8 (cost).** The family takes at most 1743 reduction samples in total. Once it exists, the operator at a new g2
comes from it in at most 60 s on one process (specialization plus the cyclic vector over Q(z)), against 240
samples for a derivation by the existing pipeline.
- Falsifier: more samples, or more than 60 s per new value.
- Confidence: moderate. The sample count rests on PK5 and on two primes being enough for the heights.

## Novelty status

- Landau singularities, normal thresholds and pseudo-thresholds, and the phase-space exponent near a threshold:
  KNOWN, standard (Landau, Nucl. Phys. 13 (1959) 181; Eden, Landshoff, Olive and Polkinghorne, The Analytic
  S-Matrix, 1966). Recalled, not read here; to be checked against the sources before any disposition cites them.
  PK1 to PK3 are checks of standard analysis on this family.
- Integration by parts with symbolic masses, and multivariate rational reconstruction over finite fields: KNOWN,
  standard in particle physics. Recalled, not read here: von Manteuffel and Schabinger, Phys. Lett. B 744 (2015)
  101, arXiv:1406.4513; Peraro, JHEP 12 (2016) 030, arXiv:1608.01902, and FiniteFlow, JHEP 07 (2019) 031,
  arXiv:1905.08019; Klappert and Lange, FireFly, Comput. Phys. Commun. 247 (2020) 106951, arXiv:1904.00009;
  Kira 2.0, arXiv:2008.06494; FIRE6, Smirnov and Chuharev, Comput. Phys. Commun. 247 (2020) 106877,
  arXiv:1901.07808. To be checked against the sources before any disposition cites them. Route (b) is their method,
  written small. PK5, PK6 and PK8 are measurements on this family, at most a modest finding.
- The D = 3 kite with half-integer density lines and squared geminal lines: as in record IB, no prior treatment
  found in a brief search, not to the standard of a novelty claim.

## Controls

- **Identification at P2 (gate for scoring PK1 to PK4).** The P2 operator's recurrence reproduces the direct
  Taylor coefficients of the P2 kite, 6 to 40 from 0 to 5, to 1e-40 relative. The direct series comes from
  `debug/pf_4b_series_direct.py` run through a new wrapper with the P2 masses. If this fails, PK1 to PK4 are not
  scored until the cause is found.
- **Wrapper control.** The same wrapper at the base masses reproduces the stored 280-digit series, coefficients
  0 to 40, to 1e-50. (`pf_common` binds the masses as default arguments at import, so patching module constants
  alone may not move them; this control catches that.)
- **Negative controls.** The base operator applied to the P2 series fails the recurrence by many orders of
  magnitude. The P2 operator with one coefficient of P_5 times 201/200 fails it too.
- **Family, positive.** PK6 (i), exact specialization to the saved base systems.
- **Family, negative.** The family specialized at s = 9 + 1/1000 fails PK6 (i) at every prime.
- Every rebuilt rational function passes held-out points, in z and in s.

## Caps stated in advance

- One topology, D = 3, one parameter free in stage 2. A family in all four masses is not attempted.
- Not rigorous. Genuineness is by local series counts modulo a prime with a root of the factor; logarithms are not
  checked.
- The cut rule is checked at one generic point and one family member. It is not proved.
- Times on one shared machine, with another job running; processor times per sample are reported alongside wall
  times.

## Compute budget and STOP rules

- At most 8 worker processes at any time. Processes I did not start are never touched.
- Stage 1: at most 4 primes x 80 samples (an extra prime only if the held-out prime disagrees), plus the two series
  runs, in at most 2.5 h wall. If the fourth prime still disagrees: STOP, report, do not score PK1 to PK4.
- Stage 2a (PK5): 35 samples on one prime.
- Stage 2b (the family) runs only if the projected count 2 x 36 x (d_s + 5) + 15 is at most 1743 and the projected
  wall time at the throughput measured in stage 2a is at most 5 h. Otherwise STOP and report the projection; PK5 is
  scored, PK6 to PK8 are not run.
- If two primes do not reach the heights and a third full prime would push the total past 2400 samples: STOP and
  report.
- Whole phase: at most 9 h wall. Exceeding it is a STOP, reported with what was reached.
- Output files: new names only (`debug/pko_*.py`, `debug/data/pko_*`, `debug/data/ib_pko_*`). The existing lift
  script writes `debug/data/ib_pf_op_Q.json` and is not rerun; a copy with a new output name is used.

## Risks

- **Record OE's cost line mixes units.** The three prime runs were wall times of six-worker pools that overlapped,
  so 2836 s is not time on one process. Their processor time is up to about 6 x 2823 s, roughly 17,000 s, which is
  above OE3's threshold of a tenth of the quadrature's bound (11,841 s). Measured per sample on an emptier machine,
  240 x 27.8 s is about 6,700 s, which is below it. The quadrature's bound was taken under load too. OE3 and the
  Paper 59 sentence ("about 2900 s, about a fortieth") need restating on like-for-like terms. This record does not
  edit either; it flags it.
- An anomalous threshold could exist at generic masses that the base point's coincidences hid. That is what PK1
  tests.
- Heights at P2 may exceed what two primes reconstruct (denominators 6 instead of 2). The held-out prime catches it;
  the budget allows a fourth.
- The direct series at P2 may need more nodes than at the base point; the wrapper control sets the node count.
- The family's degree in s could exceed 19, which would push the cost past the cap. PK5 measures it first.

## Disposition, stage 1 (appended 2026-10-26; predictions above unchanged)

**PK1 MISSED. PK2 MISSED. PK3 MISSED. PK4 MISSED.** The cut rule held on every point it names, and so did the exponent
rule. What failed is a ninth genuine point that no cut predicts, at z = -111/1120, together with what it changes.
Identification gate passed. Scripts `debug/pko_run_p2.sh`, `debug/pko_lift.py`, `debug/pko_series.py`,
`debug/pko_recurrence.py`, `debug/pko_landau_leading.py`; data `debug/data/ib_pko_p2_p{0,1,2,3}_op.json`, the
operator `debug/data/ib_pko_p2_op_Q.json`, logs `debug/data/pko_p2_*.log`, `debug/data/pko_landau_p2*.log`.

- **Derivation.** Four prime runs, one after another, Pool(6) each, 80 samples each: 576, 531, 643 and 802 s wall
  (320 samples, 2552 s). Two primes did not reconstruct the system. Its heights reach 72 bits, against 37 at the
  base point. So the system was lifted from three primes and checked on the fourth, which agreed. This is the
  extra prime the budget allowed. The exact operator over Q(z) was built in 46 s: order 11, coefficients up to
  779 bits, agreeing with the operator at 4 of 4 primes.
- **Identification gate.** The operator's recurrence gives Taylor coefficients 6 to 40 from 0 to 5, agreeing with
  the direct P2 series to 2.3e-69 relative, the precision of that series. Deviation, disclosed: the first try fed
  0 to 5 at 70 digits, and the recurrence amplifies rounding about 100 times per coefficient here, so it reached
  1.6 at n = 40, although n = 6 agreed to 9e-63. Inputs at 170 digits (two node counts, agreeing to 9e-170) and a
  250-digit recurrence fixed it. Wrapper control: at the base masses the wrapper reproduces the stored 280-digit
  series, coefficients 0 to 40, to 2.6e-69. Node control at P2: 120x120 against 160x160 nodes, 7.9e-70.
- **Negative controls.** The base operator on the P2 series: 2.9e49. The P2 operator with its lowest coefficient of
  P_5 times 201/200: 7.6e97; P_11 the same way: 1.8e82. Deviation, disclosed: my first choice of coefficient
  (z^41 in P_5) cannot enter the first 40 coefficients, so that control could not fail. It was replaced by the
  lowest nonzero coefficient, before the replacement was run.
- **PK1, MISSED.** Factoring the leading coefficient over Z and testing each factor (`debug/ib_apparent.py`): all
  eight predicted points are genuine. The four adjacent values are not roots of the leading coefficient, so the
  strong form holds too. But z = -111/1120 is a ninth genuine point, a simple root with exponents -7/2, 0, 1, ...,
  9. The other roots, one factor of degree 98, are apparent.
- **PK2, MISSED.** 18 masters and order 11, as predicted. The degrees are 109 to 120, not 102 to 113. The leading
  coefficient is z^5, the four three-line factors to the power 1, the four two-line factors to the power 3, as
  predicted, then (1120z + 111) and an apparent factor of degree 98, not 92.
- **PK3, MISSED.** At each of the eight points the exponents are exactly those in the table. At infinity they are
  exactly the base point's. At z = 0 they are -1/2, 0, 1/2, 1, 3/2, 2, 5/2, 3, 7/2, 4, 5. The base point's -4 is
  replaced by 7/2.
- **PK4, MISSED.** The system's denominators have z-roots -25 and -1/9 to the power 2, the other six predicted
  points to the power 1, and -73/9 = -(a^2+g2^2) and -277/36 = -(b^2+g1^2) to the power 1, as predicted, with no
  adjacent value. So the base point's -25/4 was b^2+g1^2, as read. But z = 0 has power 1, not 3, and
  (1120z + 111) appears to the power 3.
- **The ninth point.** It is not a leading Landau singularity of the first type. With sympy Groebner bases, the
  Landau equations dF/dalpha_i = 0 with every alpha and U nonzero have no solution at P2, for the kite and for
  every graph with one line contracted. With two lines contracted, only the two-line cuts come back. The test is
  degenerate for the massless Coulomb line. A regulator mass^2 of 1e-6 on that line still gives none. Its origin
  is open. A guess, recorded here and not registered as a prediction: at the base point it could sit on z = 0,
  where the base point has the unexplained exponent -4 and the system has z to the power 3. Stage 2's family runs
  through the base point and should show how this point moves with g2.
- **Novelty.** The rule that held is standard Landau and threshold analysis, so those clauses are CHECKs. A
  genuine point outside the normal thresholds was not sought in the literature.
- **Reading.** At a generic point the cut rule gives every genuine point of a normal-threshold type, with the right
  exponents, and no adjacent pairing. It is not the whole singular locus. The base point's coincidences hid one
  more genuine point. Its location is a rational function of the masses that is not yet known. Without that, a
  blind prediction of the leading coefficient, and so of the operator's degree, is not available.

## Disposition, stage 2 (appended 2026-10-26; predictions above unchanged)

**PK5 MISSED. STOP for stage 2b by the record's rule. PK6, PK7 and PK8 NOT RUN.** Script `debug/pko_family.py`;
data `debug/data/pko_fam_deg_p5.json`, log `debug/data/pko_fam_deg_p5.log`. The fitting and specialization script
for stage 2b, `debug/pko_family_fit.py`, was written but not run.

- **PK5, MISSED.** At one fixed random z modulo p5, 35 random values of s = g2^2, with (a, b, g1) = (1, 2, 3/2): 30
  to fit, 5 held out. All 108 nonzero entries were rebuilt as rational functions of s and passed their held-out
  points, with no failed fits. 18 masters at every sample. The largest numerator degree plus denominator degree is
  26, against the predicted at most 19. It comes from the vector c, the kite written in the masters, at 26
  (numerator 12, denominator 14). The matrix A, the differential system itself, stays at 16 (numerator at most 7,
  denominator at most 9), inside the prediction.
- **STOP.** With d_s = 26 the planned family needs 2 x 36 x 31 + 15 = 2247 samples, above the cap of 1743. At the
  measured rate, 35 samples in 297 s wall on six workers, that is about 5.3 h, above the cap of 5 h. Either one
  stops stage 2b under the rule. Stage 1 also showed that the heights can need a third prime, which would make it
  about 3350 samples.
- **What the measurement says about the route.** The system A depends on g2 with low degree. The extra degree is in
  c. So a family of A alone would cost 2 x 36 x 21 + 15 = 1527 samples on two primes, but the operator also needs c.
  c could be had from the same reducer in other ways, for example by choosing a different master basis or by
  reducing the kite only at the new value of s. That would cost a few samples per new value, not a derivation.
  This was not tested and is not claimed.
- **Compute used, whole phase.** Stage 1: 320 reduction samples, 2552 s wall on Pool(6), plus a 46 s lift and about
  2 min of series runs. Stage 2a: 35 samples, 297 s wall. All within the 9 h cap, never more than 7 worker processes
  of mine at once.
