# Deriving the two matched steps of the OP52 counting rule

Diagnostic, 2026-11-18. Nothing registered; no paper, record or git touched; no SEALED file opened. One heavy
process (about 2 min of piece-operator guessing on MO's cached (AA|BC) sequences).
Scripts: `debug/op52_deriv_{lemmas,local,endpoint,nextorder,aabc,twopieces,twopieces_b}.py`.
Data: `debug/data/op52_deriv_{local,endpoint,nextorder,aabc_A1_1073741441,aabc_A2_1073741381,twopieces}.json`.

## 0. Bottom line

Both steps follow from one mechanism, and the mechanism predicts exact rational Taylor coefficients that the
operators reproduce. **The half-integer local solution of each endpoint piece at x\* is t^{1/2} times the endpoint
part of the Landau-pinch germ of the potential that piece sees** (t = x\* - x; the endpoint part is taken at the sphere
through C for D1, through A for D2). From this:

- M1: a piece carries the germ of one potential at one endpoint. The germ is built from local data only: the
  potential's odd Taylor coefficients, C's cusp, and e^{-zB r} on the segment. So it does not depend on which
  exponential constant multiplies the piece, and the two halves of a family carry **the same** branch series (exactly,
  to 14 terms).
- M2 within one endpoint: the germ is linear in the potential. Two potentials at one endpoint therefore differ by the
  germ of their difference, so v = (p_diff - p_lead)/2.
- M2 between endpoints: the D1 and D2 germs of one potential agree at leading order and differ at the next by
  rho_D2 - rho_D1. This is an explicit rational function of the geometry, nonzero at all four geometries computed,
  so v = 1 off a proper hypersurface of parameter space.

**Verified exactly.** 3 to 6 Taylor coefficients per object, mod two primes, at A1 and A4 for f' and at A1 and A2 for
(AA|BC). The predicted rationals have up to about 27 digits, so a chance match has probability about 1e-9 per
coefficient. As a control, the D2-endpoint value in place of the D1 one fails.

**Not proved.**
- That each piece's (-1)-sector is one-dimensional. This is argued (germ rigidity) and checked exactly mod p. It is
  not automatic: see the counterexample in Section 1.
- That the endpoint-germ identification holds to all orders. It is checked to relative order |a|^2 (one subleading
  order).
- That the pinch expansion converges (L52 cap (b)).

## 1. The language (Lemmas 1-3, proved)

**Lemma 1 (counting).** Let L be regular singular at x\* with indicial roots in (1/2)Z. Local monodromy T has
eigenvalue e^{2 pi i rho} for each root rho. Hence:

- The number of half-integer roots, counted with multiplicity, is dim W, where W is the generalized (-1)-eigenspace
  of T on the local solution space.
- If W is log-free, then W = t^{1/2} H with H in C{t}, and the half-integer roots are {1/2 + ord h : h in H}. This set
  has exactly dim H elements, namely the distinct pivots of a basis of H in echelon form by valuation.

Exponent values are read from that echelon basis, not from any spanning set. Logs were checked absent at every piece,
block and full operator used here (`piece_logs`, `P4_*_logs` all false). The check is done by Frobenius recursion at
resonances, not read from the indicial polynomial.

**Lemma 2 (sums).**
- W(LCLM(L_1, ..., L_k)) = sum of the W(L_i).
- For the generic combination f = sum c_i f_i, Sol(L_f) is a subspace of the sum of the Sol(L_i), possibly a proper
  one.
- If L_f annihilates every x\*-carrying f_i, then W(L_f) = sum of the W(L_i).

**Lemma 3 (collision).** Take normalized branch series h1 and h2 of equal valuation. Their span has exponents
{e, e + v}, where v = ord(h2 - h1) - ord(h1). So a collision resolves at the valuation of the difference, not at "+1"
by fiat. Sympy check (`op52_deriv_lemmas.py`): v = 1, 2, 3 give {3/2, 5/2}, {3/2, 7/2} and {3/2, 9/2}.

**M1 is not a lemma.** A single function f = (x\*-x)^{1/2} log(x - x1) has minimal operator of order 2 with **two**
half-integer exponents {1/2, 3/2} at x\*. Its continuation around x1 adds 2 pi i (x\*-x)^{1/2} (sympy, same script).
So M1 holds exactly when a piece's (-1)-sector at x\* is one-dimensional, and that is a property of the piece. The
rank-one toys of OP52 K6b have it by integration by parts: the branch is a rational multiple of one square root. The
K6c rank-2 toys read {3/2, k} by Lemma 1 at once: H = span{(x_t-x), (x_t-x)^{k-1/2}(x_u-x)^{1/2}}.

## 2. The endpoint-germ mechanism

**Germ to next order (derived here, extending L52).** Consider a potential term s^p against C's cusp term
|r-C|^q, with density factor g = (-zC)^q/q! e^{-zB R}. Write |a|^2 = 2 D1 D2 t and
R0^2 = (1-u) D2^2 + u D1^2. Then

    I_pq ~ C_pq |a|^{p+q+3} int_0^1 u^{(p+1)/2} (1-u)^{(q+1)/2}
             [ g(R0) + u(1-u)|a|^2 ( zB/(2 R0) g(R0) - lap g / (2(p+q+5)) ) ] du.

The two correction terms come from:

- the x-dependence of the segment point: R^2 = R0^2 - u(1-u)|a|^2;
- the isotropic t-Hessian: int t^2 (t^2+M^2)^{-lambda} d^3t equals (3/2)/(lambda - 5/2) M^2 times the leading
  integral, with lap g = (zB^2 - 2 zB/R0) g.

Substituting r = R0 makes every u-integral elementary, int_{D1}^{D2} e^{-zB r} poly(r, 1/r) dr = F(D2) - F(D1). The
D1 part is -F(D1) and the D2 part is F(D2), each with its exponential removed. This split is canonical, because an
antiderivative of the form e^{-zB r} P(r) has a unique P. No E1 term arises (asserted in the code). This matches the
measurement that no E1 piece carries x\*. The constants C_pq are pi times a rational, so all predicted ratios are
rational.

**Claim (verified, not proved): the piece at endpoint E carrying constant c has as its half-integer local solution
t^{1/2} G^E_V, the E-part of sum_{p,q} v_p(V) I_pq.** V is the potential that piece sees: 1/s for the NA D1 piece;
e^{-mu s}/s (f') or S (f) for the two mu-halves at D1; V' (f') or V_A (f) for the D2 pair. The normalized series of
each piece are computed from its own embedded or guessed operator by Frobenius recursion.

| quantity (predicted exactly, matched mod both primes) | f' A1 | f' A4 | (AA\|BC) A1 | (AA\|BC) A2 |
|:--|:--|:--|:--|:--|
| P2: the two mu-halves have identical normalized branch series (14 terms) | yes | yes | yes (7 terms) | yes |
| P2: the two D2 halves identical | yes | yes | yes | yes |
| D1: [t^{5/2}] of B_mu - B_N | -55222849332/5770496875 | -294034080676/66094864875 | **0** (v_1(V_A) = 0) | **0** |
| D1: [t^{7/2}] of B_mu - B_N | match | match | -1220992024002715071/14440668429687500 | -300547039232/3604453125 |
| D1: [t^{9/2}] of B_S - B_N | -- | -- | match | match |
| D2: rho_D2 = next/leading coefficient | 325288437828681/30021182916725 | match | 342361942152264847/30637470212302525 | match |
| rho_D2 - rho_D1 (the collision coefficient) | -1.33e20/7.95e20 (nonzero) | nonzero | 2.88e26/2.46e28 (nonzero) | nonzero |
| control: the D2-endpoint value in place of the D1 one | fails | fails | -- | -- |

"Match" means the predicted rational reduces to the measured residue, on both primes for f' and on the one prime run
for each (AA|BC) case (A1: p = 1073741441; A2: p = 1073741381). Full values are in the JSON files.

**What this derives.**

1. **M1 for families.** The two mu-halves of Y (or of Y with dY) at D1, and the two e^{+-zC D1} halves at D2, carry
   the same germ. The germ contains no mu- or zC-exponential; those enter only through the constants. So a family
   adds one dimension, not two. Previously this was matched; it is now an exact identity of series, and its value is
   predicted.
2. **M2 within an endpoint.** By linearity, B_V2 - B_V1 is the germ of V2 - V1 divided by the common lead. Its
   valuation is (p_lead(V2 - V1) + 3)/2. That gives v = 1 for f' (V' has an s^1 term) and v = 2 for (AA|BC) (V_A
   starts at s^3; the [t^{5/2}] coefficient is exactly 0). The amplitude of the new exponent is predicted as well.
3. **M2 between endpoints.** The D1 and D2 germs of one potential share their leading exponent. Their difference
   after reduction is (rho_D2 - rho_D1) times the leading coefficient times t^{(next)}. Since rho_D2 - rho_D1 is a
   rational function of (D1, D2, zB, zC, zA), and is nonzero at the four geometries, it vanishes only on a proper
   hypersurface. So the between-endpoint collision resolves at +1 generically. This is no longer an assumption: it
   is a checkable non-coincidence of two explicit rationals, and it was checked.
4. **Blocks and full operator.** By Lemma 2, the half-integer local space of the D1 block, the D2 block and the full
   operator (both weight seeds) equals the span of the pieces' series exactly, to 14 terms (f', both cases, both
   primes; `P4_*_span_equals_sum`). The exponents {3/2, 5/2} / {5/2} / {3/2, 5/2, 7/2} then follow by echelon
   reduction. For (AA|BC), the same germs give {3/2, 7/2} / {7/2} / {3/2, 7/2, 9/2}.

## 3. Why the full operator misses 2 of 21 pieces, yet x\* is weight-independent

By Lemma 2, the minimal operator of a generic sum has Sol(L_f) inside the sum of the Sol(L_i), and the inclusion can
be proper when the pieces' solution spaces are linked. The full operator then need not annihilate every summand. The
two pieces it misses are the same pure exponentials at every case and prime: A1 e^{-162/25} and e^{-9/5}; A4
e^{-239/50} and e^{-32/25}. Each has order 5, and x\* is not a singular point of either, so W = 0 for both. Missing
them can change only the integer sector and other points, never the x\* half-integers.

The guard-pass parenthetical said the two pieces "differ by a rational factor". I could not confirm that. There is no
relation b = r0(x) a + r1(x) T a with deg r <= 4 against any other piece, and no b_L = R(L) a_L with deg R <= 6
(`op52_deriv_twopieces*.py`). The link exists, since the solution space is a proper subspace, but it is unidentified.
The parenthetical should be weakened to that.

## 4. Status

**Proved (lemmas, unconditional).**
- Lemma 1: the half-integer count is dim of the generalized (-1)-eigenspace; values are the echelon valuations.
- Lemma 2: LCLM and generic-sum behaviour.
- Lemma 3: a collision resolves at the difference valuation.
- M1 is not automatic (counterexample).
- The next-order germ formula (derivation above, given L52's conditions (a) and (b)).

**Reduced to an exact, checked condition.**
- The between-endpoint "+1" holds if and only if rho_D2(V) != rho_D1(V). These are explicit rationals, nonzero at
  f' A1 and A4 and at (AA|BC) A1 and A2.
- The within-endpoint step holds given the lowest odd power of the difference potential (sympy table).

**Verified exactly mod p, argued, not proved.**
- The endpoint-germ claim: each piece's x\*-local solution is the E-part of the pinch germ of its potential. This is
  checked to relative order |a|^2 (3 to 6 coefficients per object, 4 geometries).
- The reason it holds is an argument. The engine's per-constant split and the segment integral's endpoint split come
  from the same radial antiderivatives e^{b r} poly(r).
- The pieces are factorially divergent formal series, so "the germ of a piece" is a statement about its operator.

**Still assumed.**
- That each piece's (-1)-sector is one-dimensional. The argument is germ rigidity: any continuation that branches at
  x\* does so by the same local pinch, so it gives the same germ up to scale. This is measured exactly (one simple
  half-integer root, log-free, per piece) but not proved.
- That no other piece carries x\* (measured).
- Minimality of the guessed operators, and agreement over Q (the record's caps).

## 5. Proposed Paper 59 sentence and tier

Replace "the rule's two steps are supported by this prediction, not derived" with:

> Both steps follow from one mechanism, verified exactly: the half-integer local solution of each endpoint piece at
> x\* is the endpoint part, at the sphere through C or through A, of the Landau-pinch germ of the potential that piece
> sees, computed to one order beyond the leading term. A family of pieces therefore carries one branch. Two potentials
> at one endpoint differ by the germ of their difference, so the second exponent sits at the valuation of the
> difference potential. The two endpoints' germs of one potential differ at the next order by an explicit rational,
> nonzero at every geometry computed. These predicted coefficients match the operators exactly, modulo two primes, at
> two geometries for each object. That each piece has a one-dimensional branch is argued, not proved.

**Tier: unchanged, [ARGUMENT + MEASURED].** It is not [SYMBOLIC PROOF], because:

- the operators are guessed and read mod p;
- the germ identification is checked to finite order;
- per-piece one-dimensionality is argued.

The content moves from "matched" to "mechanism with exactly predicted coefficients". The PI decides; this memo makes
no edit.

**Guard.** The existing guard checks exponent sets and the counting rule as a pure function, not the mechanism. A new
fast test (seconds) would recompute, from the embedded f' piece operators and MO's (AA|BC) operators, the
normalized branch series by Frobenius recursion, and assert:

- the predicted rationals for [t^{5/2}], [t^{7/2}] (D1 difference), rho_D2 and rho_D2 - rho_D1;
- that the (AA|BC) [t^{5/2}] is zero;
- that the halves are equal.

It would reject the D2 germ substituted at D1, a "+1" collision within D1 for (AA|BC), a zero between-endpoint
residual, and two branches per family.

The (AA|BC) half needs either MO's per-piece operators embedded or about 1 min of guessing per case. Per Section 9 of
CLAUDE.md, this test should be written and fire-tested as a separate step.
