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
