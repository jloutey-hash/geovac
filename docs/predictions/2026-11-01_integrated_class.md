# Pre-registered prediction: does the integrated three-centre integral leave the two-centre class? (IC)

**Registered 2026-11-01** (PI: "Let's start the first two still open."), before the instrument below is run on any
component of (XY|XZ) or of the one-electron three-centre integral.
**Commit status:** to be committed on its own, before any computation it scores.
**Rests on:** Paper 59 sec:scope (each Legendre coefficient of (XY|XZ) is a rational combination of 31 two-centre-class
constants; the order-8 one-electron operator NA, confirmed bottom-up BU; the order-41 recurrence B2 and the
rho^L L^-4 decay B2X; AO3's generic operator beyond about 1970 unknowns and its seven component operators);
sec:f12 (the two-centre exchange class builds the same fibre family and still closes at genus zero); sec:obstruction
(the radial core is a period of an IRREDUCIBLE order-4 IRREGULAR operator, eq:pf); records NA, BU, AO3, B2, B2X;
`docs/qa/v5.20.18_claim_delta.debt.md` (the open question: "whether the bond-angle operator carries an elliptic
factor"); `debug/lit_scan/prior_art_bond_angle_series_memo.md` (Silverstone 1967: the channels are two-centre class;
no source decides the integrated class).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope; no molecular encoding is proposed.

## Why

Every channel of the bond-angle series is in the two-centre class, so any higher-class content of the integrated
value must come from the infinite channel sum. Paper 59 leaves the integrated class OPEN and names the test: the
bond-angle operator. AO3 could not find the generic operator (more than about 1970 unknowns); its 31 constant
components have operators of orders 6 to 13 (26 found at L <= 1000, R <= 24; five pure-exponential ones not found).

## What is decidable (the design argument)

1. **The value at one geometry: not decidable here.** Whether one number lies in a ring of periods is a transcendence
   question no finite computation settles. A PSLQ of (XY|XZ) at one rational geometry against products of two-centre
   seeds (instrument iii of the brief) is weak both ways: the seed basis is not finite or known (arguments are linear
   forms in D1, D2 and |Y - Z| with exponent-dependent coefficients, E1 products, logs), a negative says only that
   the chosen basis failed, and the positive control (XY|XY) validates a basis of two-centre arguments only. Not used.
2. **The function of the angle: decidable one way, and the other way up to weight.** Fix rational D1, D2, zX, zH; let
   w = |Y - Z|, so x = cos(gamma) = (D1^2 + D2^2 - w^2)/(2 D1 D2). A two-centre-class closed form is a finite sum of
   products of rational functions of w, exponentials e^{lambda w} and iterated integrals of these (E1 and log of
   linear forms in w). Such a function is annihilated by an operator that is TRIANGULAR over Qbar(w): all its
   composition factors have order one. So:
   - **NOT TRIANGULAR => leaves the two-centre class** (indeed the class of iterated integrals of hyperexponentials
     over Qbar(w)).
   - **TRIANGULAR => compatible with the two-centre class, not proof of it.** Weight-two polylogarithms are also
     triangular (control N4 below), so weight one is not certified.
3. **Why "elliptic factor" is not tested literally.** In an exponential variable elliptic content need not show as an
   order-2 Fuchsian factor: the radial core's own operator eq:pf is irreducible of order 4 and irregular, with no
   order-2 factor (control P6). And a genus-zero object can leave the Liouvillian class too: Bessel K0 is irreducible
   of order 2 and non-Liouvillian (control B1). So the decisive question is triangularity; "elliptic versus Bessel" is
   a second, classificatory reading of whatever piece blocks it (rule below, validated on controls).
4. **Components carry the answer.** f = sum_j c_j G_j, G_j = sum_L r_{L,j} P_L(x) (formal; some r grow factorially).
   (a) If every G_j's operator is triangular, their LCLM annihilates f and is triangular: f is triangular. No
   assumption on the constants. (b) f is holonomic in (D1, D2, x), so it has an annihilator in Q(x)[d/dx] at rational
   D; applied to f it gives sum_j c_j (rational series) = 0, so under Qbar-linear independence of the 31 constants
   it annihilates every G_j and f's own minimal operator carries every component's composition factors. So a
   non-triangular component means a non-triangular f, under that assumption (Paper 59 assumes Q-independence; this
   needs Qbar).

## Instrument (debug/intclass_lib.py, frozen with this record)

- **Minimal operator of a component** (`minimal_operator`; in the scorer `operator_sets`): the AO3 search
  (`debug/na_modguess.py`, unchanged) gives the smallest order findable within the rows, which can exceed the minimal
  order. Found in the controls: for P8 it returned order 13, degree 42, where the minimal operator has order 12; the
  extra order is a left multiplier M o L_min, whose own factors would sit on top. AO3's nullity-2 entries (for
  example (13, 59, 2)) are of this kind. Remedy: the gcrd of independent annihilators. Set A: every null vector at
  the smallest findable order R0 and the minimal-degree operators at R0+1, R0+2. Set B: R0+1, R0+2, R0+3.
- **Pullback** to w (x = a0 + a2 w^2), exact mod p.
- **Order-one peel** (`peel`): search a hyperexponential solution y = e^{lambda w} prod f_i^{e_i} N(w). lambda runs
  over the F_p roots of the slope-one edge polynomial at infinity and 0. e_i runs over the small-rational indicial
  classes at each F_p-irreducible factor f_i of the leading coefficient (computed in F_p[t]/f_i; pole bounds from the
  most negative integer root). N is a polynomial with degree bound from the indicial polynomial at infinity, found by
  an exact null space. Then right-divide L = Q o (d - y'/y) exactly (the remainder is checked to be 0) and repeat on
  Q. An order-1 remainder is itself an order-one factor. TRIANGULAR if the order reaches 1 or 0. STUCK if an
  operator of order >= 2 has no hyperexponential solution: its minimal submodule then has dimension >= 2, so an
  irreducible composition factor of order >= 2 exists (Jordan-Hölder).
- **Clean:** a STUCK counts only if no flag was raised at the stuck step (edge polynomial split over F_p, every
  exponent a small rational, degree bound <= 6000).
- **Component verdict:** TRIANGULAR if either set's operator peels (a left multiple can only add factors). STUCK if
  both sets are clean STUCK with the same residual order. Otherwise STUCK-UNCONFIRMED, or NOT-FOUND.
- **Classification of a stuck residual** (`classify_residual`, reported):
  - fuchsian-elliptic-type: order 2, one exponential type, a forced logarithm at a finite rational point. Irreducible
    order 2 with unipotent local monodromy means the Galois group contains SL2.
  - confluent-genus1-type: four or more distinct lambda (eq:pf-like).
  - confluent-genus0-type: order 2, two lambda (Bessel-like).
  - fuchsian-higher: one exponential type, order > 2.
  - other.

## Controls (measured before registering; `debug/intclass_controls.py`, data `debug/data/intclass_controls_<p>.json`)

Legendre controls are run at the unscored S1 geometry (D1, D2) = (9/5, 13/5), through the full pipeline (sequence ->
gcrd minimal operator -> pullback -> peel). None touches (XY|XZ) or NA.

| Control | object | expected | measured (p = 1073741789, 1073740783, 1073741833) |
|:--|:--|:--|:--|
| P1 | sum rho^L C(2L,L) 4^-L P_L = period of y^2 = t(rho-t)(1-2xt+t^2) | STUCK | STUCK, residual 2, fuchsian-elliptic-type, all three |
| P3 | P1 + N1 | STUCK | STUCK after 3 peels, residual 2, fuchsian-elliptic-type |
| P7 | P1 + N6 (order 7) | STUCK | STUCK after 5 peels, residual 2, fuchsian-elliptic-type |
| P8 | P1 + N8 (order 12, w-degree 119) | STUCK | STUCK after 10 peels, residual 2, fuchsian-elliptic-type, clean |
| P4 | P1 twisted by e^{3w} | STUCK | STUCK, residual 2, fuchsian-elliptic-type (lambda = 3 double) |
| P5 | Y's minimal operator M13 (K3 R_V + elliptic z1, z2 + order-one F's) | STUCK | 5 order-one peels, STUCK at 8 = 4 + 2 + 2, fuchsian-higher |
| P6 | the radial core eq:pf, rho = 1/5 | STUCK | STUCK, residual 4, confluent-genus1-type (4 lambda); on 1073740783 (= 3 mod 4) flagged "2 lambda outside F_p", so the scoring primes are = 1 mod 4 |
| B1 | Bessel K0 | STUCK | STUCK, confluent-genus0-type |
| N1 | sum rho^L/(L+1) P_L (log of a rational function of w) | TRIANGULAR | TRIANGULAR |
| N4 | sum rho^L/(L+1)^2 P_L (weight-two polylog) | TRIANGULAR | TRIANGULAR |
| N6 | sum_k k rho^L/(L+k) + 3 rho^L/(L+2)^2 (order 5) | TRIANGULAR | TRIANGULAR |
| N2a-c | the two-centre Coulomb integral (1s_Y 1s_Y \| 1s_Z 1s_Z) vs the angle: its 3 constant components (two grow factorially) | TRIANGULAR | TRIANGULAR, each |
| N3a-b | e^{-2w}/w: its 2 components | TRIANGULAR | TRIANGULAR |
| N5 | E1(2w): its 4 components (two E1 constants, two exponentials) | TRIANGULAR | TRIANGULAR (the E1 components reduce to order 1: weak) |
| N8 | N2b + Yukawa(3) components + N5 exp + N6 (order 10, w-degree 79) | TRIANGULAR | TRIANGULAR |

All pass on all three primes, about 45-65 s per prime. The Coulomb and E1 splits were checked in floats against
the closed forms (300-400 digits): Coulomb residual 4e-18 at L <= 60 (truncation), E1 to the truncation.

## What I know before comparing (disclosed)

- AO3 (W1, one prime, L <= 1000, R <= 24, reported only): the 31 constants, and per constant the smallest findable
  (R, Q0, nullity). There are 26 found, at (6,16), (8,30), (10,72), (12,46), (12,60), (13,59), (13,61), with
  nullity 2 at (12,46) and (13,59). Five pure-exponential constants had none: e^{-141/25}, e^{-147/25},
  e^{-197/25}, e^{-469/50}, e^{-581/50}. The generic combination has no operator below about 1970 unknowns (W1, W2).
- NA: order 8, degree 31 at N1 and N2; x* its only finite branch point (exponents 3/2, 5/2); x = +-1 and one rational
  point not apparent; infinity irregular. B2/B2X: integer beta only at x*; the slowest term is (x* - x)^{7/2}.
- The peel has NOT been run on any operator of (XY|XZ) or NA, nor have any W2 components been split. The one run of
  the scorer's own code on real coefficients was a mechanics dry run at the unscored S1 geometry with L = 40,
  Rmax = 2. All 31 components came out NOT-FOUND by construction (too few rows for any operator), and it wrote to a
  scratch directory. Timing: S1 coefficients to L = 200 took 19 s; a bisection at L = 2000 takes about 3 s for an
  unfound order (random sequence).
- Engine blobs (git hash-object, unchanged since AO3): ao_exact.py 60432979, na_exact.py 4841b06f,
  na_modguess.py bf5d7f89, na_modanalyze.py b63b4ed6.

## Cases

| Case | role | D1 | D2 | zX | zH |
|:--|:--|:--|:--|:--|:--|
| W1 | scored, (XY\|XZ) | 1.7 | 2.5 | 1.4 | 1.1 |
| W2 | scored, (XY\|XZ) | 1.6 | 2.2 | 1.0 | 1.2 |
| N1 | scored (IC-4), one-electron NA | 1.8 | 2.6 | 1.3 | 1.0 |
| N2 | scored (IC-4), one-electron NA | 1.5 | 2.4 | 1.0 | 1.2 |
| S1 | controls geometry only | 1.8 | 2.6 | 1.3 | 1.0 |

Scoring primes 1073741789 and 1073741833 (both = 1 mod 4). (XY|XZ): L <= 2000, R <= 40. NA: L <= 1200, R <= 24.

## Predictions

- **IC-0 (controls; run first).** `python debug/intclass_controls.py <p>` on both scoring primes: ALL PASS.
  Falsifier: any control misread. Probability 0.97.
- **IC-1 (the gate).** W1 is NOT TRIANGULAR: at least one component is STUCK (clean, confirmed by both operator sets)
  on both primes. Falsifier: every W1 component TRIANGULAR on both primes (a decided TRIANGULAR verdict). An
  INCOMPLETE outcome leaves it UNDECIDED. Probability 0.6.
- **IC-2.** W2 gets the same case verdict as W1. Falsifier: different decided verdicts, or one decided and one
  INCOMPLETE. Probability 0.8.
- **IC-3 (the elliptic reading; scored only if IC-1 holds).** At least one STUCK W1 residual is
  fuchsian-elliptic-type or confluent-genus1-type. Falsifier: every stuck residual is confluent-genus0-type,
  fuchsian-higher or other. Probability 0.35.
- **IC-4.** The one-electron integral NA is TRIANGULAR in both N1 and N2. Falsifier: a STUCK component in either.
  Probability 0.6.
- **IC-5.** Every component verdict that is decided on both primes is the same on both (W1, W2). Probability 0.95.
- **IC-6.** All 31 W1 components have an operator within L <= 2000, R <= 40. Probability 0.45.

**Reported, not scored:** per component the minimal order and degree, the peel steps (lambda, exponent gauges, degree
of N), every stuck residual's classification, lambda set, forced-log points and adjoint test, and the run times.

## Decision gate

- **VOID** if IC-0 fails.
- **GO** if IC-1 and IC-2 hold: the integrated (XY|XZ), as a function of the bond angle, leaves the two-centre
  (indeed the triangular-Liouvillian) class at two geometries, on two primes, under Qbar-independence.
- **STOP** if IC-1 is missed by a decided TRIANGULAR verdict: the prediction failed. Read: the function of the angle
  is compatible with the two-centre class (weight not certified).
- **BORDERLINE** otherwise (INCOMPLETE, unconfirmed STUCKs only, or IC-2 missed).

## Phase 2 (after this record is committed)

```
python debug/intclass_controls.py 1073741789 ; python debug/intclass_controls.py 1073741833      # IC-0, ~1 min each
python debug/intclass_score.py run W1 1073741789      # ~2-2.5 h: coefficients ~0.5-0.9 h, 31 searches ~1.3 h, peels
python debug/intclass_score.py run W1 1073741833
python debug/intclass_score.py run W2 1073741789
python debug/intclass_score.py run W2 1073741833
python debug/intclass_score.py na N1 1073741789 ; python debug/intclass_score.py na N1 1073741833  # ~10-20 min each
python debug/intclass_score.py na N2 1073741789 ; python debug/intclass_score.py na N2 1073741833
python debug/intclass_score.py score
```
One thread each; memory well under 4 GB (the largest search cache is about 2000 x 1900 entries). Each run writes
`debug/data/intclass_<case>_<p>.{json,log}` incrementally. The peel time on real components is not measured. On
controls of the same order (12-13, w-degree up to 120) it is under 1 s; a peel slower than 2 h stops at its budget
(verdict BUDGET, counted as not found).

## Caps

- Function-level class in the angle at fixed rational lengths and exponents; 1s functions; two geometries; two
  primes. The value at a geometry is not decided.
- TRIANGULAR does not certify weight one (polylogs of any weight are triangular).
- NOT TRIANGULAR needs the 31 constants to be linearly independent over Qbar (assumed, as in Paper 59, there over Q).
- The minimal operator is a gcrd of annihilators, not certified minimal. A spurious left multiplier can only make a
  verdict STUCK, never TRIANGULAR. The A/B confirmation guards against this but does not exclude a multiplier common
  to both sets.
- The hyperexponential search covers small-rational exponents and exponential types that split over F_p. Anything
  else is flagged, and flagged STUCKs do not count.
- The classification is a reading, validated on the controls above. It is not a proof that a residual is a
  pullback of an elliptic family.
- NOT-FOUND means only that the operator lies beyond L <= 2000, R <= 40.

## Disposition

Scored 2026-11-01 by `python debug/intclass_score.py score` (output `debug/data/intclass_score.{json,out}`), after
the Phase 2 runs below. Instrument files unchanged since this record's commit (sha256 as in its message).
**Decision: GO.** IC-0, IC-1, IC-2, IC-3, IC-5, IC-6 held; IC-4 missed.

| Case, prime | components | STUCK (clean, A = B) | STUCK-UNCONFIRMED | stuck residuals (order, reading) | time |
|:--|:--:|:--:|:--:|:--|:--:|
| W1, 1073741789 | 31 | 6 | 25 | 6 x (6, other) | 5451 s |
| W1, 1073741833 | 31 | 31 | 0 | 6 x (6, other); 25 of order 10-16, confluent-genus1-type | 5484 s |
| W2, 1073741789 | 31 | 6 | 25 | 6 x (6, other) | 5238 s |
| W2, 1073741833 | 31 | 6 | 25 | 6 x (6, other) | 5326 s |
| N1, both primes | 6 | 6 | 0 | orders 4-6, other / confluent-genus1-type | ~400 s each |
| N2, both primes | 6 | 0 | 6 | flagged: lambda edge does not split over F_p | ~400 s each |

- **IC-0 HELD.** All controls pass on both scoring primes.
- **IC-1 HELD.** Six W1 components are clean STUCK on both primes, each with a residual of order 6. So W1 is NOT
  TRIANGULAR. Every other component was unconfirmed on the first prime only because some exponential rates at
  infinity do not lie in that prime field (flag "lambda edge has k roots outside F_p"). On the second prime all 31
  are clean STUCK.
- **IC-2 HELD.** W2 is NOT TRIANGULAR on the same evidence: six clean STUCK components on both primes, order-6
  residuals. The other 25 are flagged on both primes.
- **IC-3 HELD, with a scope note.** The genus-one readings all come from the second prime. The six components
  that are STUCK on both primes read "other" on both. The reading is the record's own classifier, not a proof
  (Caps).
- **IC-4 MISSED.** N1 is NOT TRIANGULAR: all six components are clean STUCK on both primes, with no factor peeled.
  N2 is INCOMPLETE: all six are flagged. The one-electron integral, as a function of the angle, also leaves the
  triangular class. Its rates at infinity are lambda = 0 and lambda^2 = -69/100 = -(1.3^2 - 1.0^2), set by the two
  orbital exponents of the N1 case.
- **IC-5 HELD.** No component is decided differently on the two primes.
- **IC-6 HELD.** All 31 W1 components have an operator within the caps (minimal orders 6 to 20).

**Erratum to Caps (found after scoring, before this disposition).** The Caps say the hyperexponential search
flags everything outside its coverage. That is false. The search assigns one exponent to each F_p-irreducible
factor of the leading coefficient. A solution whose exponents differ at the conjugate roots of a factor of degree
two or more is missed without a flag. Planted control: 4(c - w^2) y'' - 4 w y' + y = 0, with solutions
(w - s)^(1/2) and (w + s)^(1/2), s^2 = c. It is TRIANGULAR, and the instrument returns it clean STUCK when c is not
a square mod p (c = 5 at 1073741833). With c = 2, a square there, it peels in one step. Rescaled, this is
Chebyshev's equation of order 1/2. The gap and its cure, working over the extension field, are standard
(Cluzeau and van Hoeij, J. Symb. Comput. 38 (2004) 1043; their recurrence paper, Sec. 8). The miss needs a factor
of degree two or more with a non-integer exponent class, or exponents outside F_p. **Check on the scored
residuals** (rebuilt at L = 600; probe scripts in the session scratchpad, not committed): no such factor and no
exponent outside F_p in any of W1's six both-prime STUCK residuals (both primes), W2's six (second prime) or
N1's six (second prime). In every probed residual all exponents are integers. In W1 and W2, the two ln
components have only degree-one singular factors, where the search is complete. So no scored
STUCK is this miss, and every verdict above stands. Every verdict above depends on that check, not on the Caps
statement as registered.
