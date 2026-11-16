# Record IN: deciding record IC's undecided one-electron geometry N2 (inventory row 4)

*Registered 2026-11-15, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14).
A fresh reviewer read the draft cold before commit; its findings were applied: a pass-side control with a
degree-2 factor, clause (a) of the erratum check and its control E, IN-3 in the gate and UNDECIDED on missing
runs, the reported fields and incremental writes in the driver, IN-0 relabelled a precondition.
Rests on record IC (`docs/predictions/2026-11-01_integrated_class.md`), its frozen instrument and its erratum.*

## Why

Record IC scored the one-electron three-centre integral NA as a function of the bond angle at two geometries. N1
was decided: all six constant components are clean STUCK on both primes, so NA leaves the triangular
(two-centre-type) class there. N2 (D1 = 1.5, D2 = 2.4, zX = 1.0, zH = 1.2) came out INCOMPLETE: all six components
were STUCK but flagged "lambda edge has 2 roots outside F_p" on both primes, so none counted.

**The flag, read before computing (an expectation, reported not scored).** The N2 residuals' rates at infinity
include lambda = 0 and +-6/5 (found) and two roots outside F_p. IC's N1 rates were lambda^2 = -69/100 =
zH^2 - zX^2 (zX = 1.3, zH = 1.0). The same form at N2 gives lambda^2 = 1.44 - 1 = 11/25, lambda = +-sqrt(11)/5,
and 11 is a quadratic non-residue modulo both IC primes 1073741789 and 1073741833 (Legendre symbol -1 at both,
checked). If that is the missing pair, primes in which 11 is a square split the edge, as IC chose primes = 1 mod 4
for P6. Both new primes are 1 mod 4, so 11 and -11 are squares together and the sign of lambda^2 does not matter.
No instrument change. If the edge still fails to split on the new primes, IN-1 is UNDECIDED, not missed.

## Instrument (frozen; unchanged since IC)

`debug/intclass_lib.py`, `debug/intclass_score.py` (`run` in `na` mode, L = 1200, R <= 24), `debug/na_exact.py`,
`debug/na_modguess.py`, `debug/na_modanalyze.py`, `debug/intclass_controls.py`, at this commit. The new driver
`debug/intclass_in_check.py` (frozen with this record) recomputes each component with those functions and adds
only the erratum check below; it does not change any verdict logic.

**Primes:** 1073741717 and 1073741689: both prime, both = 1 mod 4 (P6's +-2i splits), 11 and -69 squares at both
(N2's and N1's edges split).

## Controls (measured before registering)

IC's full control suite (`python debug/intclass_controls.py <p>`; IC's table, all expectations unchanged) on both
new primes: **ALL PASS** on 1073741717 (63 s) and 1073741689 (58 s); `debug/data/intclass_controls_<p>.json`. P6
splits with four lambda and reads confluent-genus1-type, clean, on both.

**The erratum check.** A STUCK residual PASSES iff (a) at every F_p-irreducible factor of its leading coefficient
no indicial exponent lies outside F_p, and (b) at every factor of degree >= 2 all exponents are integers. Clause (b)
is IC's erratum. Clause (a) was added on the cold review: exponents such as +-sqrt(c) at a rational point are
hyperexponential over Qbar(w), and IC's search returns them STUCK without a flag (control E below). IC's own scored
verdicts are not touched by it: its addendum found integer exponents at every point of the probed residuals.

Erratum-check controls (`python debug/intclass_in_check.py control <p>`, `debug/data/intclass_in_controls_<p>.json`),
**ALL PASS** on both primes (c below is the smallest non-square: 2 at 1073741717, 7 at 1073741689):
- IC's planted 4(c - w^2) y'' - 4 w y' + y, c a non-square: STUCK, check FAILS (degree-2 factor, exponent 1/2).
  With c a square (6, 2) the search peels it (TRIANGULAR).
- Pass side with a degree-2 factor: Legendre nu = -1/2 (elliptic K type), (c - w^2) y'' - 2 w y' - y/4, c a
  non-square: STUCK, check PASSES (exponents 0, 0 at w = +-sqrt(c)).
- E: w^2 y'' + w y' - c y (solutions w^{+-sqrt(c)}, TRIANGULAR over Qbar(w)), c a non-square: STUCK, check FAILS
  (degree-1 factor, two exponents outside F_p).
- No-factor case, Bessel K0-type w y'' + y' - w y: STUCK, check PASSES (trivially: no factor of degree >= 2).
The check is in `debug/intclass_in_check.py`, written before this commit and run only on these controls; it calls
IC's `operator_sets`, `component_verdict` and `peel` unchanged.

## Predictions

- **IN-0 (precondition, measured before commit, not a forecast).** IC's controls and the erratum-check controls
  pass on both new primes (above).
- **IN-1 (the gate).** N2 is NOT TRIANGULAR: at least one component is clean STUCK, confirmed by both operator sets
  (A, B) with the same residual order, on both new primes, AND both sets' residuals pass the erratum check (above)
  on both primes. Falsifier: every N2 component TRIANGULAR on both primes. Probability 0.8 (N1's six
  are all STUCK; the N2 residuals were STUCK, only flagged).
- **IN-2 (positive control on the primes).** N1, rerun on both new primes, is again NOT TRIANGULAR with all six
  components clean STUCK and passing the erratum check. Probability 0.9.
- **IN-3 (prime agreement).** No component of N1 or N2 is decided differently on the two new primes, nor
  differently from a decided verdict on IC's primes. Probability 0.95.

**Reported, not scored:** per component the minimal order and degree, the residual order, lambda set and reading,
the erratum-check details (factor degrees, exponents) and the run times.

## Decision gate

- **VOID** if IN-0 fails.
- **GO** if IN-1, IN-2 and IN-3 hold: NA as a function of the bond angle leaves the triangular class at both of
  IC's geometries (two primes each, under the IN erratum check and Qbar-independence of its constants, IC's caps).
  IN-3 is UNDECIDED (not HELD) if any run is missing.
- **STOP** if IN-1 is missed by a decided TRIANGULAR verdict: the two geometries disagree.
- **BORDERLINE** otherwise.

## Caps

IC's caps, unchanged: function-level class in the angle at fixed rational lengths and exponents; 1s functions;
TRIANGULAR would not certify weight one; NOT TRIANGULAR needs Qbar-independence of the constants; the minimal
operator is a gcrd of annihilators, not certified minimal (a spurious left multiplier can only make a verdict
STUCK; the A/B confirmation guards against it, not against a multiplier common to both); the hyperexponential
search's blind spot is the erratum's, which IN-1 checks explicitly.

## Phase 2 (after this record is committed)

```
python debug/intclass_in_check.py run N2 1073741717 ; python debug/intclass_in_check.py run N2 1073741689
python debug/intclass_in_check.py run N1 1073741717 ; python debug/intclass_in_check.py run N1 1073741689
python debug/intclass_in_check.py score     # IN-0..IN-3 and the gate -> debug/data/intclass_in_score.json
```
About 400 s per run (IC's timing). IN-0 is scored from the four control files above. IN-3 compares against IC's
`debug/data/intclass_NA_<case>_<p>.json`.
