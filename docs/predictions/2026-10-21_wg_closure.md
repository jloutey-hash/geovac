# Pre-registration: closing W G = 0 for the Slater slice and the chain

**Registered:** 2026-10-21, before any computation that scores it.
**Prompted by:** the /qa paper_59 claim-search DELTA of 2026-10-20 (`docs/qa/v5.21.2_claim_delta.debt.md`,
"Open: closing W G = 0"); PI-directed ("I'll go with your recommendation ... Let's work those").
**Extends:** Paper 59 `sec:density_ring`, the parity proof (`docs/predictions/2026-10-20_slice_parity_proof.md`
and its correction note), which derives P G = 0 and P C = 0 from the measured W G = 0 and W C = 0.
**Commit status:** committed on its own, before computing.

## The question

The Mellin route gives an order-seven operator L that kills G (the slice) and factors exactly as
N2 W, with W = R3 R1 (theta + 2) (chain: N2 W, W = R3 R1 (theta - 2)). That W G = 0, so that N2 is not
needed, is measured (2e-31), not proved. Can it be proved? Two gaps: (a) that L G = 0 holds with no source
term, which the Mellin route does not guarantee (it produced the source -3 t^2 for the Coulomb triangle);
(b) that y = W G, which solves N2 y = 0 by (a), must vanish.

## What I know before computing (disclosed)

- Numerically L G and W G are below 2e-31 at five t (Paper 59; `tests/test_paper59_chain_slater_slice.py`).
- The derivation (`debug/density_loop_B_ode.py`, `mellin_t_ode`): J^2 obeys theta^n F = sum_j v_j(z) theta^j F
  (z = x^2, theta acting first); A(sigma) = M[J^2](sigma) gets the recurrence sum_k P_k(-(sigma + 2k)) A(sigma + 2k) = 0
  with boundary terms dropped; Phi(s) = M[G](s) = A(2 - s) E(s), E = M[K1], E(s + 2) = (s^2 - 1) E(s); a common
  polynomial factor in s is cancelled; the t-operator is read off with all terms moved to one line.
- Hand reasoning, not computed. (i) The A recurrence has no source if v_0 = 0 (the relation kills constants):
  then every term carries a theta, theta J^2 decays exponentially at x = infinity, so every Mellin integral
  converges on Re sigma > -2 with no boundary term, and B = M[theta J^2] = -sigma A continues as a meromorphic
  identity. J^2 tends to a constant at infinity, so constants are in the kernel; whether v_0 = 0 in the code's
  representation is not yet checked. (ii) The only other place a source can enter is the contour shift: the
  term t^(2i) Q_i(theta) G is an integral on the line c - 2i, c in (2, 4); moving it to the line c crosses
  the poles of Phi at s = 4, 6, 8 (from J's odd series at x = 0, J^2 = sum x^(2k+2) log^m x). The source is the
  sum of those residues, a finite combination of t^2, t^0, t^-2 times powers of log t, whose coefficients
  use only the singular Laurent parts of Phi there: exact numbers built from J's series coefficients
  (rationals, gamma, log 2) and the Taylor data of E at 4, 6, 8 (Gamma and polygamma values). (iii) N2 for the
  slice (read from the test file) is (t^2 - 4) theta^2 + 6(3t^2 + 28)/(2t^2 + 7) theta + (24t^8 + ...)/(2t^2 + 7)^2.
  In the form y'' + f y' + g y = 0, g tends to 6 and f to 1/t, so the formal solutions at infinity are
  e^(+-i sqrt6 t) t^(-1/2) (the code reviewer's note gives e^(+-i sqrt6 t); chain e^(+-i (2 sqrt5/3) t)).
  With y = v exp(-(1/2) int f), v'' + q v = 0 with q -> 6. If q >= q0 > 0 on [T, infinity) and q' is
  integrable there, Sturm comparison gives every nonzero real solution infinitely many zeros, and
  E = v'^2 + q v^2 obeys E' = q' v^2 >= -(|q'|/q) E, so E stays above a positive bound and limsup |y| t^(1/2) > 0.
  A function with a power-log expansion at infinity either has a nonzero leading term, and then a fixed sign
  for large t, or is O(t^-K) for every K; neither is possible for a nonzero real solution. G is real and W has
  real coefficients, so y is real. (iv) G has a power-log expansion at infinity (poles at s = 4, 6, 8, ...;
  the parity record), and so does y = W G. So (a) and (b) together give W G = 0 on t > T, hence on t > 0.

## Model

- Step 1 (homogeneity). Recompute the slice's L by `mellin_t_ode` (it must equal the test file's QG) and
  check v_0 = 0 exactly. Compute J's series at x = 0 exactly to the order needed, the singular Laurent parts
  of A at sigma = -2, -4, -6 and the Taylor data of E at s = 4, 6, 8, and the residue sum S(t). Exact arithmetic
  over Q with symbols for pi, gamma, log 2 and any zeta value the polygammas bring.
- Step 2 (oscillation). From the exact N2: f, g, q = g - f'/2 - f^2/4 as rational functions. Find T beyond the
  real singular points with q >= 5 on [T, infinity) (sign of the numerator of q - 5 by a Sturm sequence or real
  root isolation) and q' = O(t^-3).
- Step 3. The same for the chain, C = (8 t^3/pi^3) int chi^2 K0(t x) dx, chi = int_0^x s K1, with A = M[chi^2],
  E = M[K0], the t^3 factor as in `_on_chain`, and its N2 from `_chain_pieces`.

## Controls (scored first; a failed control voids the predictions that use it)

- **C1 (positive, KNOWN: Paper 59 eq. for the ring).** The same residue bookkeeping on the Coulomb triangle's
  Mellin route gives its known source, L_G G = -3 t^2, exactly. If the triangle's Mellin data cannot be set up
  in the same code, C1 is replaced by a constructed control: a function with a known nonzero source, stated
  in the disposition before its result is read.
- **C2 (oscillation test, KNOWN).** Step 2's code reports oscillation for Bessel's equation
  theta^2 y + (t^2 - 1) y = 0 (q -> 1) and reports NOT oscillatory for the modified equation
  theta^2 y - (t^2 + 1) y = 0 (q -> -1, which has the solution K1, O(t^-K) for every K).
- **C3 (representation).** The recomputed slice operator equals QG exactly; for the chain, the recomputed
  operator equals QC exactly. A failed C3 voids that object's predictions.

## Predictions

- **WG1 (slice, no source).** v_0 = 0, and the residue sum S(t) is identically zero. *Falsifier:* v_0 != 0
  with a nonzero boundary contribution, or any nonzero coefficient in S.
- **WG2 (slice, N2 oscillates).** There is T <= 10 with q >= 5 on [T, infinity), and q' = O(t^-3).
  *Falsifier:* no such T up to 10, or q' not integrable.
- **WG3 (chain).** WG1 and WG2 hold for the chain, with q -> 20/9 and the bound q >= 2 on [T, infinity), T <= 10.
  *Falsifier:* any clause fails.
- **Consequence, scored with them.** If WG1 and WG2 hold, W G = 0 is derived, and P G = 0 with it, with one
  cap left: the Mellin contour shift (vertical decay of Phi, the power-log expansion at infinity), as now.
  Same for the chain with WG3.

## Novelty status

- WG1-WG3: internal statements about Paper 59's own operators. The tools are textbook (Mellin residues;
  Sturm comparison, e.g. Hartman, Ordinary Differential Equations, Ch. XI). No literature search applies;
  a hit is a tier upgrade inside the corpus, class **CHECK**, as for MP3.
- C1, C2: KNOWN (Paper 59; Bessel's equations).

## Caps, stated in advance

- The Mellin contour shift and the vertical decay of Phi are used, not proved here (the paper's existing cap).
- Exactness of S = 0 is over the symbols named in Step 1; if an independent transcendental appears whose
  vanishing coefficient needs a relation among constants, that coefficient is reported, not assumed zero.

## Disposition, part 0: the replacement control, stated before computing it (2026-10-21)

The triangle cannot be set up in the same code. Its kernel psi solves an inhomogeneous equation,
(theta^2 - 1 - x^2) psi = -2/pi, and its Mellin factor M[psi] needs psi's series at infinity as well as at
zero; the engine built for this record handles K-Bessel kernels, with exact Gamma recurrences and series
data at x = 0 only. So C1 is replaced, by the rule above, by:

- **C1' (constructed positive control).** G1(t) = int_0^inf x K1(t x) J(x) dx, the slice with J in place of
  J^2. Its operator L1 comes from `mellin_t_ode` on J's own relation; its strip is 2 < Re s < 3 and its
  crossed poles are at s = 3, 5, 7, .... The engine's predicted source S1(t) must equal a direct quadrature of
  L1 G1 at t = 0.7 and 1.9 to 1e-20 of the largest term, at 30 digits. C1' counts as a positive control only
  if S1 is not identically zero. If it is, the next variant is G0(t) = int_0^inf x K0(t x) J(x)^2 dx, then
  G2(t) = int_0^inf x K1(t x) J(x)^3 dx, in that order, until one has a nonzero source.
- **C1'' (negative control).** The same engine, run on the C1' function with one crossed pole left out,
  must disagree with the quadrature by more than 1e-6 of the largest term at one of the two t.

## Disposition (appended 2026-10-21, after computing; the predictions above are unchanged)

Driver `debug/wg_closure_A.py`, data `debug/data/wg_closure_A.json`, log `debug/wg_closure_A.log`. Exact
arithmetic in sympy throughout; quadrature at 30 digits for the controls.

**Disclosed code error, fixed before scoring.** The first run reported C3 FAILED for the chain and crashed
in WG2. Both had one cause: the test files build their symbols t and u with different assumptions, so
sympy treated them as different symbols. With the symbols rebound by name the recomputed chain operator
equals `QC` term by term; nothing else changed between the runs.

**Controls.**
- C3 PASSED: the recomputed slice operator equals `QG`, and the chain's equals `QC`, exactly.
- C1' PASSED as a positive control: for G1 = int x K1(t x) J dx the engine predicts L1 G1 = -16/t, and
  quadrature gives -22.857142857142857142857143 at t = 0.7 and -8.4210526315789473684210526 at t = 1.9, equal
  to the engine to all 25 printed digits. The source is nonzero, so C1' counts.
- C1'' PASSED: with its one crossed pole (s = 3) left out the engine gives 0, against -22.86.
- C2 PASSED: Bessel's equation, q = 1 - 3/(4t^2), oscillation reported (T = 1.2247); the modified equation,
  q -> -1, reported not oscillatory.
- Extra check, not registered: C1' has only simple poles, and the slice has a double pole at s = 8. For
  G0 = int x K0(t x) J^2 dx, which has the same poles, the engine gives L0 G0 = -9216 t^2, and quadrature
  gives -4515.84 and -33269.76 at t = 0.7 and 1.9, equal to all printed digits. Without the double pole the
  engine gives 9216(-t^4 - 3)/t^2, so that branch is exercised.

**Predictions.**
- **WG1: HELD.** v_0 = 0 for J^2 (and for J and chi^2). The six crossed residues (s = 4 for i = 1, 2, 3;
  s = 6 for i = 2, 3; s = 8, a double pole, for i = 3) sum to zero exactly. The cancellation is not
  term by term: leaving out the residues at s = 4, 6 or 8 gives 8640 pi + 55440 pi/t^2, -8640 pi - 10080 pi/t^2
  and -45360 pi/t^2. So L_G G = 0 with no source. Class **CHECK**.
- **WG2: HELD.** q -> 6, q = 6 - 59/(4 t^2) + ..., q' = O(t^-3); q >= 5 on [T, infinity) with T = 3.4179
  (the largest real root of q - 5's numerator, denominator and N2's leading coefficient). Class **CHECK**.
- **WG3: HELD.** Chain: v_0 = 0 for chi^2; the six crossed residues (s = 0, 2, 4, with a triple pole at 4)
  sum to zero exactly; q -> 20/9, q = 20/9 - 85/(12 t^2) + ..., q' = O(t^-3), q >= 2 on [5.3156, infinity).
  Class **CHECK**.
- **Consequence.** W G = 0 and W C = 0 now follow, by the argument in "What I know": y = W G solves N2 y = 0
  (WG1), is real, and has a power-log expansion at infinity; every nonzero real solution of N2 has
  infinitely many zeros past T and amplitude bounded below by a multiple of t^(-1/2) (WG2), so y = 0 on
  t > T, hence on t > 0. With the parity record, P G = 0 and P C = 0 follow. The one cap left is the Mellin
  contour shift (vertical decay of Phi, the power-log expansion at infinity), as before. The tier of the
  slice and chain membership moves from "derived from the measured W G = 0" to derived, with that cap.
