# Cold review: Pachucki exchange bound, general exponents / a -> b / l > 0

Date 2026-11-18. Reviews `debug/pachucki_extension_memo.md`. Reviewer scripts: `debug/pachucki_ext_review_*.py`.
Nothing else touched, except that re-running `pachucki_ext_symbolic.py` rewrote the untracked
`debug/data/pachucki_ext_symbolic_1s.json` (identical apart from `time_s`).

**Verdict: SOUND WITH CAVEATS.** No error found in any claimed result. The caveats: one check the proof needs
and the memo ran only for K, and wording.

## What I ran

| check | script | result |
|---|---|---|
| Part A re-run | author's `pachucki_ext_symbolic.py` | reproduces every Part A line |
| own integrals | `review_neumann.py` (Neumann expansion, no geovac import, no closed form) | helper's 1s value at (1,2), R = 1.4, to 30 digits |
| l > 0 map at a new point | `review_closedform.py`, (a,b,R) = (4/5, 3/2, 23/10) | 1s, 2p-sigma, 2p0-2p0, 2p+1-2p+1, 3d0: closed form = my integral, relative difference <= 5e-48 |
| scaling, stored Q(t) coefficients | `review_scaling.py`, t = 15/8, aR = 46/25 | a K(aR;1,t) = my K(R;a,b) to 48 digits, five cases; a^2 or a^0 is 20-25% off |
| FM1 sign | `review_closedform.py` | f0 by label-60 quadrature = direct integral to 50 digits; sign flipped: 11% off |
| a = b | `review_scaling.py` | Sugiura form = my integral at a = b = 1, R = 7/5, to 55 digits |
| every t > 1 | `review_qt.py`, 1s and 2p-sigma | all f_alpha, F_alpha, d_alpha p^2 denominators are t, t-1, t+1; every key decays for all t > 1 |
| jets at repeated exponents | `review_jets.py` | jet sources = sympy derivatives of the helper's formulas, exactly, at (1,2) and (4/5,3/2), 16 multi-indices each, up to (3,1,3,1) |

## Answers

**a.** The Q(t) argument is a proof for every real t > 1, with one addition. Pivots and the `_neg` heuristic do
not matter: the residual is checked after the solve, as an identity in Q(t), and it specialises wherever no
coefficient has a pole. Uniqueness by decay is used at every multi-index, so the proof needs the denominators
and the decay of every intermediate f_alpha. The memo reports K's denominators only. I checked the intermediates
for 1s and 2p-sigma; the three 576-index cases were not checked (about an hour each).
"Order <= 9" needs no rank: ten coefficient vectors in a nine-dimensional space over Q(t0)(R) are dependent at
every t0. "Exactly 9" is not proved at any t; it needs the nine functions independent over C(R). The memo says so.
*Wrong answer accepted by the zero residual:* any wrong source, since the solver solves whatever F it is given.
Source correctness rests on the helper (subsets only), the jet check above and numerics.

**b.** K(R;a,b) = a K(aR;1,b/a) holds for normalised orbitals of any l. f0 scales as a^-2. I re-derived the
four normalisations. Every numeric check in the author's memo has a = 1, so the power of a was untested there;
it is tested above.

**c.** Correct. With A at the origin and B at z = R, r_B^2 = r_A^2 - 2zR + R^2, so
z_A = (r_A^2 + R^2 - r_B^2)/(2R), z_B = z_A - R. Each density is a polynomial in r_A^2, r_B^2 with Laurent
coefficients in R, and r^n e^{-ur} = (-1)^{n+1} d_u^{n+1}[e^{-ur}/r] gives the map. End-to-end: table above.

**d.** The integral bounded is int int conj(2p+1_A)(1) 2p+1_B(1) (1/r12) conj(2p+1_B)(2) 2p+1_A(2), the
exchange integral of the complex pair; the orbital product is (1/2) N_A N_B (x^2+y^2) e^{-a r_A - b r_B}, hence
the 1/4. It is EQ's `pppi0`. It is not the real 2p_x-2p_x exchange integral, whose density x^2 also has an
m = +-2 part, (x^2-y^2)/2, which is not covered.

**e.** No circularity. The guessed operators come from exact Neumann-series data mod primes; the derived ones
from the closed form; `exchange_hp` is a Neumann sum and never sees the closed form. My integrals use the same
textbook expansion and share no code with either.

**f.** Legitimate. K is continuous in b at b = a and equals the closed form for t > 1, so the eps^0 coefficient
is K at a = b. The bound 4 follows from the explicit limit lying in a four-function span closed under d/dR, not
from a limit of operators. Nothing beyond a continuity check of Sugiura's formula is claimed.

**g.** The engine's FM1 is the reviewed helper's, unchanged; the fix was in the limit script's own re-entry.
Published Eq. (63) (PDF of arXiv:0905.1635v3) reads
F^(-1) = [2(u^2-w^2)]^-1 {(e^{2rw} - e^{2ru}) Ei(-2r(u+w)) - e^{-2rw}[...] + e^{-2ru}[...]}, which both now carry.
The numeric check above does not depend on the paper. Not tuned.

## Required changes to the claims

1. Sentences 1 and 3: the proof needs the intermediate denominators and decay, not only K's. Either run that
   check for 3d0, 2p0-2p0 and 2p+1-2p+1, or say the identity holds outside a finite set of t and extends by
   continuity (K's own denominators are t, t-1, t+1).
2. Sentence 3: name the 2p+1 case as the exchange integral of the complex pair and say real 2p_x-2p_x exchange
   is not covered.
3. Sentence 3: "(9, 34; 9, 53; 9, 39)" is in a different order from the integrals listed. Attach each degree to
   its case.
4. Sentence 3 tier: [SYMBOLIC PROOF] rests on sources at repeated exponents, which the reviewed helper does not
   cover. Keep the tier only with the jet-versus-sympy check as a guard.
5. Sentence 2: the equal-exponent bound is prior art and is a consistency check; keep "recovered".
6. The memo cites arXiv labels (60, 59, 28); the published numbers are 63, 62, 29. Say which.
7. m != 0 plan: Pachucki's Eqs. (64)-(66) give the r12 recursion at the exchange point with no 0/0. It holds
   only on the exchange locus, so separate exponent derivatives are still missing, but the plan should start
   there. His Eq. (67) is the first exponent derivative at the exchange point in closed form.

## What a guard must pin

* Jet sources = sympy derivatives of the helper's formulas at repeated-exponent multi-indices, exact (rejects a
  jet rule right at order 1 and wrong at order 2; the zero residual cannot see it).
* Denominators of every f_alpha in {t, t-1, t+1} and every key decaying for all t > 1 (rejects a hidden
  exceptional t: a pole at t = 3, a key e^{(b-3a)R}).
* K_1s(23/10; 4/5, 3/2) = 0.0955959302201830564664864485839835 against an independent integral (rejects a wrong
  power of a; no a = 1 pin can).
* f0(23/10; 4/5, 3/2) = 0.0022558102214482983799805712969021 (rejects the FM1 sign flip).
* Same point, l > 0: 2p-sigma 0.171366050015632245442926777969856; 2p0-2p0
  0.0427216693803001655774249080965235; 2p+1-2p+1 0.0797334869400471005133449208692647; 3d0
  0.0850202908833205952575807791934367 (rejects a shifted alpha, a sign (-1)^|n|, 1/4 read as 1/2).
* The author's list otherwise stands.
