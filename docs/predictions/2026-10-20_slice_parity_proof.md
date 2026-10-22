# Pre-registered prediction: the parity proof that the slice and the chain lie on P (MP3)

**Registered:** 2026-10-20, before any computation that scores it. PI-directed under the same request
("Love proofs. Let's do it"). Found while scoring MP1 and MP2, which failed
(`docs/predictions/2026-10-20_slice_membership_proof.md`).
**Commit status:** committed on its own, before any computation that scores it.

## The argument (derived by hand, nothing computed)

1. R1 S = 1 mod R3 means R1 S = 1 + N R3 for some operator N. Then
   R1 P = R1 (theta+2) - R1 S R1 (theta+2) = E - (1 + N R3) E = -N W, with E = R1 (theta+2).
   So for every y in Sol(W), R1 (P y) = 0, and P y = alpha(y) h, with h the solution of R1:
   h = t/(1-t^2)^4 for the slice, h = t^3/(1-t^2)^2 for the chain.
2. Both h are odd in t. At t -> infinity, h ~ t^-7 (slice) and h ~ -t^-1 (chain): odd powers of t, that
   is, half-integer powers of x = t^2.
3. P has coefficients rational in x. It maps a convergent expansion in integer powers of 1/x and
   powers of log x to another one of the same kind.
4. G's Mellin transform is 2^(s-2) Gamma((s-1)/2) Gamma((s+1)/2) M[J^2](2-s). J is a sum of odd powers
   of x times powers of log x near x = 0, so M[J^2](z) continues with poles only at z = -2, -4, ...
   The Gamma factors have no poles with Re s > 1. So G's large-t expansion has only even powers of t,
   times logs. Since G solves W and infinity is a regular singular point, that expansion is the
   convergent local one. For the chain, chi is a sum of odd powers of x times logs, and
   C(t) = (8 t^3/pi^3) int chi^2 K0(tx) dx gives t^3 times t^-(3+2k): again even powers only.
5. So P G has only even powers at infinity, and alpha h has only odd ones. Hence alpha(G) = 0 and
   P G = 0 exactly. The same for the chain.
6. Corollary: ker P is exactly the set of solutions of W whose expansion at infinity has only even powers
   of t. So P's exponents at infinity are W's integer group there.

## Predictions

**MP3.**
- (a) R1 P = -N W holds exactly, with N rational, for the slice and for the chain.
- (b) W's exponents at x = infinity (in w = 1/x) are {1, 2, 7/2, 4, 6} for the slice and
  {-1, 0, 1/2, 1, 2} for the chain, and P's are the integer ones: {1, 2, 4, 6} and {-1, 0, 1, 2}.
- (c) J near 0 has only odd powers of x times logs, and so does chi. Checked by series to order x^9.
- *Falsifier:* any of (a) to (c) fails. (b) failing would mean an error in steps 1 to 5.

Confidence: high, since (b) follows from (a) and steps 2 to 5.

## Novelty status

- **MP3: NOT FOUND** in the literature for these objects. As for MP1 and MP2, the outcome it proves was
  measured first. A hit is scored **CHECK**, a tier upgrade from measured to proved. A failure is a
  **MISS**.

## Controls

- The exponent code at infinity must reproduce the recorded infinity rows of the bound tables for R3 R1
  (slice {2, 7/2, 4, 6}, chain {0, 1/2, 1, 2}).
- P built from the perturbed section must fail (a).
- Numerically, G t^4 at t = 10, 20, 40 must approach 6 pi with no odd-power drift. This is a sanity check.

## Honest caps

- Step 4 uses the standard Mellin continuation and contour shift, with the decay of the Gamma factors on
  vertical lines. It is used, not re-proved in full.
- One slice and one chain fibre only.

## Result

(To be appended after scoring; the predictions above are never edited.)

## Disposition (scored 2026-10-20)

Scripts `debug/slice_parity_E.py` (stopped after parts (a) and (b), whose output was saved; its sympy
Bessel series was too slow) and `debug/slice_parity_F_series.py`. Data
`debug/data/slice_parity_{E,F}.json`, logs `debug/slice_parity_{E,F}.log`.

**Controls.** The infinity exponents of R3 R1 reproduce the bound tables: slice {2, 7/2, 4, 6}, chain
{0, 1/2, 1, 2}. P from the perturbed section leaves a nonzero remainder in (a), for both objects. The
truncated series reproduce J and chi at z = 0.3 to 2.5e-15 and 2.6e-16. G t^4 = 18.4207, 18.7355,
18.8204 at t = 10, 20, 40, against 6 pi = 18.8496. The gaps 0.429, 0.114, 0.029 fall by a factor near 4
per doubling, as t^-2 times a log does. An odd power would give 2 or 8.

**MP3 (a): HELD.** R1 P = -N W exactly, with
N = (x-4)(42x^4 - 45x^3 + 77x^2 - 75x + 33) / (315 x (x-1)^3 (2x+7)) for the slice and
N = (x-4)(20x^2 - 7x + 5) / (15 x (x-1)(4x+11)) for the chain.

**MP3 (b): HELD.** At infinity W is {1, 2, 7/2, 4, 6} and P is {1, 2, 4, 6} for the slice; W is
{-1, 0, 1/2, 1, 2} and P is {-1, 0, 1, 2} for the chain. P drops exactly the half-integer exponent.

**MP3 (c): HELD.** J and chi carry only the powers x^1, x^3, x^5, x^7, x^9 below x^11, times logs.

**Class: CHECK** for all three, as decided before computing: a tier upgrade. P G = 0 on the Slater slice
and P C = 0 for the chain at the fingerprint fibre go from measured to proved, with the cap on the Mellin
contour shift stated above. One detail was added while scoring and is not a change of argument: on t > 1,
Sol(R1) is spanned by h, so P G = alpha h there; the large-t comparison gives alpha = 0 on t > 1, and
P G = 0 then holds on all of t > 0 by analytic continuation.

**Correction 2026-10-20 (/qa paper_59 claim-search DELTA, three reviewers independently).** The
argument starts from y in Sol(W), and that G solves W (and C its W) is measured, not proved: the Mellin
route gives the order-seven N2 W, and N2 is shown spurious only numerically (2e-31 at five t). So the
proof reduces two measured facts to one: P G = 0 and P C = 0 follow from the measured W G = 0 and
W C = 0. The tier upgrade above is partial. Paper 59 now says so. A closure was suggested by the code
reviewer and is not attempted here: N2 is irregular at infinity (formal solutions like e^(+-i sqrt6 t)),
so a nonzero W G could not carry G's power-log expansion; that would also need L_G G = 0 to be exact,
which is itself derived by the Mellin route and checked numerically. Any such attempt is a new registration.

**Forward pointer (2026-10-21).** The premise W G = 0 (and W C = 0) named in the correction above is now
derived: `docs/predictions/2026-10-21_wg_closure.md` (WG1-WG3 held). The one cap left is the Mellin
contour shift.
