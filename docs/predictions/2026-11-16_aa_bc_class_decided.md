# Record AK: record AB's class clause (P4) decided on primes that split its rates

*Registered 2026-11-16, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14) and the
PI's direction to work down the owed items (2026-11-16: inventory row 6, "class as a function of the angle
UNDECIDED"). Rests on record AB (`docs/predictions/2026-11-14_aa_bc_three_centre.md`), record IC's frozen instrument
(`docs/predictions/2026-11-01_integrated_class.md`) and record IN's erratum check
(`docs/predictions/2026-11-15_na_n2_decided.md`). The precedent is IN, which decided NA's N2 the same way. A fresh reviewer read the draft and the driver cold
(0 LARGE; applied: AK-3's code now accepts STUCK with any flag as the text says (it had required clean STUCK),
AK-2 checks every component found, the scorer returns MISSING/UNDECIDED instead of crashing, IC's caps restated,
AK-2's falsifier names the unread components; the rate algebra and the Legendre symbols re-derived by hand, held).*

## Why

AB's P4 ran IC's instrument on the 15 mu-carrying components of (AA|BC) at case A1 on IC's primes 1073741789 and
1073741833. Every component was STUCK from both operator sets at the same residual order (order 4: eight E1 and four
ln components; order 5: e^{-162/25}; order 6: e^{-189/50}, e^{-513/50}), but every residual carried the flag
"lambda edge has 4 roots outside F_p", so none counted: P4 UNDECIDED.

**The rates, read before registering (disclosed; a diagnostic, no verdict computed).** On 1073741789, the slope-one
edge of each component's w-operator (A set) was factored by a scratch script equivalent to the driver's `edge` mode
(all 15 components: each edge has two irreducible quadratics with zero linear term; their constants read for two
components, e^{-162/25} and e^{-189/50}, and equal). Besides the split
rates (0; +-9/5 = +-mu on the two pure exponentials), the edge carries two irreducible quadratics lambda^2 - c1,
lambda^2 - c2 with c1 + c2 = 51/20 and c1 c2 = 81/25:

    lambda^4 - (51/20) lambda^2 + 81/25 = 0,   lambda^2 = (51 +- 3 i sqrt(287))/40,   lambda = (+-sqrt(615) +- i sqrt(105))/20.

All four roots lie in F_p iff 615 and -105 are squares mod p. Neither AB prime qualifies (1073741789: neither;
1073741833: -105 only), which accounts for the flag at both. The new primes **1073741441** and **1073741381** are the
two largest primes below 2^30 (sympy.isprime) that are 1 mod 4 with 615, -105, -69 and 11 all squares (the last two also split NA's rates at
both of IC's geometries; checked by Legendre symbols). The origin of these rates is not derived here.

No instrument change: IC's verdict logic, AB's frozen engine (`debug/aabc_exact.py`, unchanged since AB's commit)
and IN's erratum check (`debug/intclass_in_check.py`, unchanged since IN's registration) are called as they are.

## Instrument

Driver `debug/aabc_ak_check.py` (this commit): `run <p>` builds A1's coefficients to L = 1200 mod p exactly as
`aabc_run.py p4` does, runs IC's two operator sets and the peel on each mu-carrying component, and records IN's
erratum check on each STUCK residual. `score` writes `debug/data/aabc_ak_score.json`.

**Controls (AK-0), measured before registering on both new primes** (logs `debug/data/aabc_ak_ctrl_<p>.log`, json
`intclass_controls_<p>.json`, `intclass_in_controls_<p>.json`):

| control set | 1073741441 | 1073741381 |
|:--|:--|:--|
| IC's 21 (P1, P3-P8, B1 expected STUCK; N1-N8 expected TRIANGULAR) | 21/21 as expected | 21/21 as expected |
| IN's erratum: planted non-square STUCK, check FAILS; square peels; K0 STUCK, PASSES; Legendre-K STUCK, PASSES; Euler w^{+-sqrt c} STUCK, FAILS | 5/5 | 5/5 |

AK-0 is therefore re-checked from these files at scoring, not re-run.

## Predictions (A1 only, both new primes)

- **AK-0 (gate).** IC's and IN's controls pass on both new primes. Failure is VOID.
- **AK-1 (the class).** At least one mu-carrying component is clean STUCK from both operator sets, with the erratum
  check passing, on both new primes: the screened pieces leave the triangular (two-centre) class on their own, not
  only through NA. Falsifier: every component TRIANGULAR on both primes. Mixed or flagged outcomes leave it
  UNDECIDED. Probability 0.8.
- **AK-2 (the edge accounted for).** All 15 components are found and no residual carries a "lambda edge ... outside
  F_p" flag on either prime. Falsifier includes: a component whose edge constants differ from the two read (13 of
  15 were not read) stays flagged on the new primes. Probability 0.8.
- **AK-3 (consistency with AB).** Every component is STUCK on both new primes (any flag) at the same residual order
  as in AB's P4 on both of AB's primes. Probability 0.75.

**Gate:** VOID if AK-0 fails; GO if AK-1, AK-2 and AK-3 hold; STOP if AK-1 is missed (every component
TRIANGULAR); BORDERLINE otherwise.

## Caps

- IC's caps, unchanged: NOT TRIANGULAR rests on the Qbar-independence of the constants (E1, ln and the
  exponentials), and the gcrd is not certified minimal (agreement of sets A and B does not exclude a common
  left multiplier).
- One geometry (A1), as AB's P4. Case A2 is not run.
- A clean STUCK is evidence of leaving the class modulo the instrument's search (small-rational exponents,
  exponential types splitting over F_p, IN's erratum clauses); it is not a proof over Qbar(x).
- The whole-f verdict follows from NA's components (AB, "Already known"); it is not scored here.
- The rates' origin (why |lambda|^2 = 9/5 numerically) is not explained and not claimed.

## Results (2026-11-16; appended after scoring, text above unchanged)

Run `python debug/aabc_ak_check.py run <p>` on both primes after the registration commit (9ccbe13f) landed;
coefficients ~970 s per prime, 10-15 s per component. Data `debug/data/aabc_ak_<p>.json`, score
`debug/data/aabc_ak_score.json`.

**Verdict: GO.**

| clause | 1073741441 | 1073741381 | outcome |
|:--|:--|:--|:--|
| AK-0 controls (IC 21, IN erratum 5) | all pass | all pass | **HELD** |
| AK-1 clean STUCK, both sets, erratum check passing | 15 of 15 | 15 of 15 | **HELD** |
| AK-2 15 found, no lambda-edge flag | yes, none | yes, none | **HELD** |
| AK-3 every component STUCK at AB's residual order | yes | yes | **HELD** |

- Residual orders as in AB's P4: 4 (eight E1 and four ln components), 5 (e^{-162/25}), 6 (e^{-189/50}, e^{-513/50});
  every residual type "confluent-genus1-type", clean, no flags.
- **Reported.** Each residual's rates at infinity are the four roots found in F_p plus 0 (e^{-162/25}) or +-9/5
  (the order-6 pair); the four satisfy lambda^4 - (51/20) lambda^2 + 81/25 = 0 exactly mod p on both primes (sum
  zero), the edge read before registering.
- **Consequence for inventory row 6.** As a function of the bond angle, at A1, the screened (mu-carrying) part of
  (AA|BC) leaves the triangular (two-centre) class on its own, not only through NA's components, under IC's caps
  (Qbar-independence of the constants; minimal operators not certified) and the instrument's search caps.
