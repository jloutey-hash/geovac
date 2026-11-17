# Pre-registered: record EX's missing instrument checks at (1, 3/2), with more primes (EXb)

**Registered:** 2026-11-16, before any EXb computation (PI-approved follow-up to record EX's VOID). Cold review
before commit: 3 MUST-FIX (instrument outcomes mapped to STOP/BORDERLINE; now VOID as in EX), 4 SHOULD-FIX
(coefficient-degree screen, reconstruction ladder, 32 primes, disclosure of the height bound), 3 NIT; all applied.
**Parent:** `docs/predictions/2026-11-16_exchange_unequal_zeta.md` (record EX, scored VOID at 384e2bdb).
**Inventory:** row 3c of `docs/chemistry_period_inventory.md`. Owner: Paper 58.

## Why

Record EX was VOID on its registered instrument clause. At (a, b) = (1, 3/2), CRT plus rational reconstruction of the
mod-p operator failed with the four run primes (M ~ 2^103). So `heldout 1 3/2` produced nothing, and `numeric 1 3/2`
could not run. Every mod-p clause of EX held at both pairs, and (1, 2) passed every instrument check. The failure is
CRT headroom. EXb supplies the headroom and scores only what EX could not.

## What is already known (disclosed)

- **Scored and frozen by EX:** the mod-p readings at (1, 3/2), on four primes.
  - Operator: order 9, degree 16.
  - Finite points: apparent or exponents_only.
  - R = 0: exponents {-1 (x2), 0, 1, 2, 3, 4, 6, 9}, with ln R.
  - R^9 ln R coefficient: 9/500 by ratrec mod p.
  - Rates at infinity: {-5/2, 3, 2, 0 x2, -2 x2, -3 x2}.
- **Not looked at:** the reconstructed operator over Q and any value of it. The four-prime failure does bound the
  coefficient heights from below: the normalised operator's numerator times denominator exceeds about 2^102, so
  heights above about 2^51. `exch_heldout_1_3_2_EXS.json` and `exch_numeric_1_3_2_EXS.json` exist and contain only
  the error strings.
- **Untouched:** the sealed files `debug/data/exch_unequal_SEALED_*.npz` stay unopened. The driver neither reads nor
  regenerates them.

## Instrument

`debug/exchb_run.py` (new, frozen here) changes no frozen file. It works in three steps:

1. **Candidate runs.** It calls the frozen `debug/exch_run.py run 1 3/2 p --suffix EXBC` on candidate primes, in a
   fixed order:
   - the four EXS primes first: 67108859, 33554393, 67108837, 67108819;
   - then the primes below 2^26 in descending order, skipping P3 = 16777213.
2. **Acceptance.** A candidate is accepted when its operator has order 9 and degree 16, with per-coefficient
   degrees equal to the EXS run on 67108859 ([14, 15, 16 x 8]; this screens primes that drop a leading coefficient),
   annihilates all complete rows, and passes the valuation check. The first 32 accepted primes are used (about
   2^830); the cap is 48 candidates. Accepted runs are copied to suffix `EXB`, which is all that `heldout` reads.
   The driver refuses to start if any `_EXB` run file already exists.
3. **Instrument checks and frozen steps.** It records determinism (the EXBC operator on the four EXS primes equals
   the EXS one; a crash is recorded, not raised) and a stability ladder: reconstruction on the first k accepted
   primes for k = 18, 20, ..., 32; k0 is the smallest k that succeeds, and k0 + 4 must also succeed with the same
   operator. It then runs the frozen
   `heldout 1 3/2 --suffix EXB` and `numeric 1 3/2 --suffix EXB`, unchanged from EX:
   - `heldout`: exact residual over Q, plus a third-prime fit with 60 held-out orders;
   - `numeric`: Frobenius from the CRT operator vs `exchange_hp` at 50 and 80 digits, at R = 1.4 and 3.015.

The summary goes to `debug/data/exchb_summary.json`.

## Predictions

- **EXb0 (gate).** Determinism holds on all four EXS primes. At least 32 candidates are accepted within the cap. The
  frozen files match their sha256 below. Failure is VOID.
- **EXb1 (reconstruction).** Some k0 <= 28 on the ladder succeeds, k0 + 4 succeeds, and the two operators are equal.
  Failure is VOID (instrument). Probability 0.85.
- **EXb2 (EX1's missing parts at (1, 3/2)):**
  - **Exact residual:** `heldout` finds an exact residual of zero over Q on every order it checks (about -9..111).
  - **Third prime:** fit nullity 1, proportional to the CRT operator, and the CRT operator kills every row including
    the 60 held out.
  - **Numeric:** `numeric` agrees with `exchange_hp` (80 digits) to 1e-30 at R = 1.4 and 3.015, and the 50- and
    80-digit `exchange_hp` values agree to 1e-35.

  Probability 0.85.
- **EXb3 (EX4b exactly).** The exact-Q series' ln R part starts at R^9 ln R with coefficient exactly **9/500** =
  (8/3375) a^5 b^5. Probability 0.9. A miss is reported; it does not gate.

## Decision gate

- **VOID** (instrument, as in EX) if any of these happens:
  - EXb0 or EXb1 fails;
  - the exact Q residual is nonzero, or the third-prime check fails. The CRT operator already kills the series
    modulo 32 primes to N = 450, so either outcome means a wrong reconstruction, a mismatch between the Q and mod-p
    generators, or a bad third prime;
  - `numeric` errors, or its 50- and 80-digit `exchange_hp` values disagree beyond 1e-35;
  - a step has not finished 48 h after launch.
- **GO** if EXb0, EXb1 and EXb2 hold (EXb3 may miss). EX's own Results line stays VOID. The joint verdict, EX + EXb
  = GO for row 3c, goes in EXb's Results and in the inventory, as AX2b did for AX: every EX clause held, EX2-EX5 at
  (1, 3/2) on mod-p readings, and the instrument now covers both pairs.
- **BORDERLINE** if `numeric` completes, its 50 vs 80 digits agree to 1e-35, everything else holds, and Frobenius
  vs `exchange_hp` misses the 1e-30 bar.
- No STOP outcome: EX's scientific clauses were scored there.

## Caps

These are EX's caps: 1s only, two exponent pairs, mod-p singular-point readings, gcrd not certified minimal, and
constants treated as formal. EXb adds only CRT headroom.

## Frozen (this commit; unchanged at scoring)

The new driver, plus EX's frozen list (record EX, section "Frozen"; its sha256 table applies unchanged).

| file | sha256 |
|:--|:--|
| `debug/exchb_run.py` | `a3c0a3b98e7468fa5615a150db9a43657a3ad6312408dc3a3c9d1e6feb68cc49` |
| `debug/exch_run.py` | `41e582548d38066976444090062ffdc4ff3d3ac1ad5806bbb9b66a44f0f21c00` |
| `debug/exch_gen.py` | `5a17efc6fc9d88883d9e790a04fca3271dcbea3b3690eb44458b138e949996d5` |
| `debug/intclass_lib.py` | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `geovac/qfd_core.py` | `43317e4ecc0d6b1f47eb5221f954833cb114e8ee9ca108559334fdca4c74be24` |

## Phase 2 (after this record is committed and the commit verified)

`python debug/exchb_run.py > debug/data/exchb.out 2>&1`
