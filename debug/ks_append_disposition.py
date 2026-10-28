from pathlib import Path

FN = Path(__file__).resolve().parent.parent / "docs" / "predictions" / "2026-10-28_kite_speed_test.md"
TEXT = r"""
## Disposition (appended 2026-10-28; predictions above unchanged)

**Gate at T12: STOP.** R_marg is about 500: with the operator given, float64 quadrature is about 500 times faster.
KS2 FAILED, KS3 HELD, KS3b HELD, KS4 HELD, KS5 UNCLEAR (not scored as registered, because C7 failed at T9).
Scripts `debug/ks_common.py`, `debug/ks_route_q.py`, `debug/ks_route_o.py`, `debug/ks_time.py`,
`debug/ks_time_q_loop.py`; data `debug/data/ks_outer_{T12,T9}.json`, `debug/data/ks_timing.{json,log}`,
`debug/data/ks_timing_qloop.{json,log}`. One child process at a time; phase 2 took well under the 2 h cap.

**Shared outer grid.** T12: 40 nodes, w from 4.9e-6 to 653. T9: 36 nodes, w from 7.5e-5 to 149. With the exact
60-digit phi, stage 2 on these grids gives Y to 6.5e-14 to 3.8e-13 (T12) and 1.2e-13 to 4.6e-13 (T9).

**Ladders (accuracy only, chosen before any clock was read).**
- Route Q, max relative error on Y over the four t: h = 1/4: 3.7e-9 (fails both targets); h = 1/8: 4.6e-13 (T12),
  5.3e-13 (T9). Scored rung h = 1/8 at both targets.
- Route O at T12: 20, 25 and 30 digits give nan or inf; 40 digits gives 1.8e278; 60 digits gives 3.8e-13. At T9:
  20 to 40 digits blow up (inf, inf, 8.7e263, 4.8e138); 60 digits gives 4.6e-13. Scored rung 60 digits at both.
  The failure below 60 digits is a blow-up, not a gradual loss. Its cause was not diagnosed.
- Route O used 76 main-path steps (80 in all, 2 branch walks) at T12, 69 (73) at T9, against 248 in record OE.

**Controls.**
- C1 PASS (above). C2 PASS, at 60 digits only. C3 PASS: O and Q agree on Y to 1.4e-14 to 7.4e-14 at all four t, and
  both meet the target against the reference.
- C4 PASS: perturbed operator, Y wrong by 1.4e-2 to 6.8e-2 (both targets). C5 PASS: slice assignment, Y wrong by
  6.9e-2 to 1.1e-1.
- C6 PASS: Q at h = 1/16 over h = 1/8, median ratio 5.0 (registered run); 4.94 in the looped run below.
- C7: T12 PASS (O 1.14, Q 1.00). T9 FAIL: after ten repeats Q's max/min was 1.5 (O 1.10). As registered, the T9
  timings are UNCLEAR and the T9 gate is not scored.

**Timer resolution (not foreseen; disclosed).** On this machine `time.process_time()` ticks at 15.625 ms. Route Q's
whole evaluation is 2 to 3 ticks: every registered Q repeat read 0.03125 s, one read 0.046875 s. That is the T9 C7
failure. As a supplement, not part of the registered scoring, Q was repeated 50 times inside each of five children
(`debug/ks_time_q_loop.py`): median 33.1 ms (T12, spread 1.08), 163.8 ms at h = 1/16 (spread 1.13), 30.9 ms (T9,
spread 1.07).

**Timings (processor time, one process, medians).**

| | T12 | T9 |
|:--|:--|:--|
| Route O load (operator, roots) | 0.08 s | 0.08 s |
| Route O input (six coefficients, 70 digits) | 0.22 s | 0.23 s |
| Route O walk (60 digits) | 16.08 s | 16.02 s |
| Route O total, spread | 16.39 s, 1.14 | 16.34 s, 1.10 |
| Route Q total, registered (tick-limited) | 0.031 s, 1.00 | 0.031 s, 1.50 |
| Route Q total, looped supplement | 0.0331 s, 1.08 | 0.0309 s, 1.07 |
| R_marg, registered / looped | 525 / 495 | 523 / 528 |
| R_tot (with 9,200 s derivation), registered / looped | 2.9e5 / 2.8e5 | 2.9e5 / 3.0e5 |

Imports before each run: about 0.30 s (O) and 0.23 s (Q), excluded as registered. Stage 2 was below one tick for
both routes.

**Scores.**
- **KS2 FAILED.** Route O needs 60 digits at both targets; nothing at 40 or fewer reaches either. **MISS.**
- **KS3 HELD.** R_marg = 525 at T12. **UNCLEAR** (novelty not searched at scoring).
- **KS3b HELD.** R_marg > 10 by a factor of about 50. **UNCLEAR**, same reason.
- **KS4 HELD.** R_tot = 2.9e5; route Q costs 0.03 s, not 9 s. **UNCLEAR**, same reason.
- **KS5 UNCLEAR.** Not scorable as registered (C7 failed at T9). In the looped supplement the T9 outcome is the same
  STOP (R_marg 528).

**Departures from the model, all disclosed.**
- Route O's detour around a real apparent root is walked only if 0.65 rho is below the largest node's z. The shorter
  grids end before the later roots, and walking past the last node would only add steps. Otherwise oe_engine and the
  oe_kite walker were imported unchanged, including its 1200-bit floor for the local operator coefficients.
- The input coefficients were computed in-process by `pf_4b_series_direct.u_chunk` (80 x 80 'sub', one process),
  not by a subprocess as in record OE, so that process_time counts them.
- The complex128 port of the walk was not attempted. The flint walk fails at 40 digits, so a 16-digit port cannot
  reach either target.
- The looped Q timing is an addition to the plan, forced by the timer tick. The registered scoring stands without it.

**Reading.** The 1/13 of record OE is a property of 60-digit arithmetic, where the rival was 60-digit ball
quadrature. At chemistry precision the kite's integrand is a smooth 2-D integral, and a graded float64 rule gets Y to
5e-13 in about 33 ms. The operator walk cannot drop below 60 digits on this path, and at 60 digits it costs about
16 s, nearly all of it in the walk. Counting the derivation, which is owed afresh for every parameter set, the
quadrature is about 3e5 times cheaper. This says nothing about Chen et al., whose integrals are a different object.
"""

s = FN.read_text(encoding="utf-8")
assert "## Disposition" not in s
FN.write_text(s.rstrip("\n") + "\n" + TEXT, encoding="utf-8")
print("appended")
