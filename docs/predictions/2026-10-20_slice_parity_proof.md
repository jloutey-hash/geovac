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
