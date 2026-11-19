# Cold review: `debug/pachucki_m_nonzero_memo.md` (Pachucki route for m != 0 exchange)

Date 2026-11-18. Reviewer did not write the memo. No paper, test, doc, CHANGELOG, CLAUDE.md or git state touched;
files under review not edited; SEALED files not opened; `debug/exact9_m_*` not touched.

## Verdict: SOUND WITH CAVEATS

The mathematics held under every attack I ran. The caveats are about wording and guard design. None of them
blocks a paper sentence that claims the bound (order <= 9, the nine-function span for every 0 < a < b, and operator
equality at (1,2)). Two items must not be stated as established: "exactly 9" for general a < b, and operator
equality away from (1,2). The memo already lists both as open.

## Attack table

| # | attack | result |
|---|---|---|
| 1 | Re-derive (**) from label 39 | **tried, held.** See below. |
| 1' | Uniqueness of the decaying solution / sign of c/q | **tried, held.** At the exchange point c = 0, so every jet solve has leading operator q0 d^2/dr^2. The homogeneous solutions are A + B r, never e^{+-p r}, and no decaying solution exists. Off the exchange point the c-jets sit on the right-hand side (Leibniz), so the sign of c/q never enters. Instrumented: `_solve_conds` left **0 free unknowns in all 248 solves** at (1,2) and at (2/3,7/4). The uniqueness argument is analytic (true jet minus found solution = A + B r, decaying, so 0) and does not depend on solver rank. |
| 2a | Are the controls circular? Planted wrong coefficients | **tried; held, with one wording defect.** Table below. |
| 2b | Is `exchange_hp` independent? New point | **tried, held.** `exchange_hp` is the prolate-spheroidal Neumann (tau) expansion in `geovac/qfd_core.py`, with no Pachucki equation. At the new point (a,b,R) = (2/3, 7/4, 1.9) the closed form agrees with exchange_hp (tau_max 16, dps 32, tail ~3e-34) to every printed digit (30): 2ppi -0.00515968878885911303770792463990781 vs -0.00515968878885911303770792463991, and 3ddel 0.000355893627169501579264221099213784 vs 0.000355893627169501579264221099214. |
| 3 | Can a small error hide under the 1e-6 quadrature? | **tried, held.** Every planted coefficient error moves the value by O(1) relative or more (table below). The quadrature is a secondary check. The binding numerical checks are 47-48 digits (memo) plus 30 digits at my new point. A 1e-6 check would catch every plant, though not an error of order 1e-8. Nothing suggests such an error exists, and the exact controls rule out the coefficient classes that would produce one. I did not re-audit the quadrature script; it imports the engine only for the comparison value. |
| 4 | Normalisation and m conventions | **tried, held.** The orbital tuples in `pachucki_m_hpref.py` equal `REF_ORBITALS` in `tests/test_paper58_exchange_l1.py` (2ppi) and `tests/test_paper58_exchange_d_pp.py` (3ddel) character for character. The normalisations re-derive by hand: N_2p^2 = z^5/pi, N_1s^2 = z^3/pi, and the integral of \|x+iy\|^4 e^{-2zr} is 12 pi/z^7, so N_3d^2 = z^7/(12 pi). With CS phases, 2ppi = -(2p_x 1s\|2p_x 1s). **Independent of the memo's operator comparison:** I read the committed operators from the test files as text and applied them to the rebuilt closed form as an exact identity in the function algebra. OP_P_1_2 annihilates 2ppi and OP_3DDEL_1_2 annihilates 3ddel. Each cross-family operator leaves 13 nonzero keys. Neither committed record stores a value; the committed tests compare against exchange_hp, which I used too. |
| 5 | Q(t) claims | **tried, held.** I specialised intermediate Q(t) jets from the memo's state pickles at new t and compared them with a fresh Q-field solve at (1,t), a separate arithmetic path. All of these are equal exactly: 3ddel at t = 7/3, f_4(1111), f_2(3131), f_0(3333); 2ppi at t = 11/4, f_2(1111), f_2(1011). For decay near t = 1, the memo's checker requires each linear rate e(t) to satisfy e(1) <= 0 with slope <= 0, not both zero. That gives strict decay for every t > 1, pointwise but not uniformly. The denominators t, t+-1 have no zero in t > 1. The pinned f0 coefficient of 2ppi re-derives exactly at (1,2) and at (2/3,7/4); the second pair tests the a-dependence. |
| 6 | 2p_x 2p_x sigma = 2 part | **tried, held.** Algebra: x1^2 x2^2 averaged over the common azimuth gives rho1^2 rho2^2 [1/4 + (1/8) cos 2(phi1-phi2)], so (xx\|xx) = pppi0 + ppdel/2. Numerically, I expanded each real 2p_x = (2p_{-1} - 2p_{+1})/sqrt2 into the 16 complex exchange_hp terms. The expansion is valid whatever slots exchange_hp conjugates, because x is real. At (1,2), R = 1.4, tau 8: 16-term 0.15072359406719310298797 vs pppi0 + ppdel/2 0.15072359406719310281957, diff 1.7e-19 (tau-truncation level). ppdel = 3(b/a)^2 3ddel follows from the normalisations above. |
| 7 | sympy heugcd patch | **tried, held.** I forced the fallback on 300 random and structured univariate integer pairs, including (t-1)^5, (t+1) and content factors. The fallback gcd equals sympy's Poly.gcd up to sign, and the cofactors are exact. A Q(t) field computation with the fallback forced matches sympy `cancel`. Exactness of the cofactor division is asserted in the patch itself. Maximality affects only reduction, not correctness. The dense-list orientation (reversed both ways) is consistent. |

## Attack 1 in detail: independent derivation of (**)

From label 12, with u1 = t: sigma = sigma0 + sigma2 w1^2 + t^2 w1^4, where sigma0 = c - q t^2 matches label 18, and
sigma2 = t^4 - Sig t^2 + (u2^2 u3^2 + w2^2 w3^2 - u2^2 w2^2 - u3^2 w3^2) = t^4 - Sig t^2 + e. All of these match
the memo. With g = sum (-1)^n G_n w1^n/n! (label 40), the w1^m coefficient of
sigma g' + (1/2) sigma' g + P, multiplied by m!(-1)^{m+1}, is:

sigma0 G_{m+1} + sigma2 G_{m-1}[m(m-1) + m] + t^2 G_{m-3}[m(m-1)(m-2)(m-3) + 2m(m-1)(m-2)] + (-1)^{m+1} P_m = 0

The brackets give m^2 and m(m-1)^2(m-2), and every sign is +. Inverse Laplace for r > 0: t^2 G -> f'' and
t^4 G -> f''''. The boundary terms are polynomials in t, so they map to deltas at r = 0, as do the polynomial-in-t
parts of P_m. The result is exactly (**).

At the exchange point (q = -(u^2-w^2)^2, c = 0, e = 0, Sig = 2(u^2+w^2)), (**) reproduces label 61 differentiated
twice. That includes the n^2 and (n-2)(n-1)^2 n coefficients, so label 61 independently confirms both m-dependences
at the base jet.

I re-checked the four-logarithm grouping of label 13 by hand: the t-carrying arguments, the signs of A3..A6 / D3..D6,
and the t*w1 coefficients of the two eps = 1 logs. L^{-1}[ln(t+c)] = -e^{-cr}/r, and the k-th w1-derivative of the
eps = 1 logs gives (-1)^{k-1} r^{k-1} e^{-cr}. Both match `RJets.source`.

The Leibniz implementation in `solve_fn` (q, c on U != 0; Sig, e on all U; d^4 and f_{m-3}'' on U = 0 only; binomials
on derivative values) is correct.

The termwise inverse Laplace is legitimate because f_n jets are integrals with positive exponents. They are analytic
in the exponents and in w1 near 0, bounded at r = 0, and decay at infinity.

## Planted-coefficient table ((1,2), R = 1.4; `pachucki_m_review_plants_1_2.json`)

"caught" means the exact identity fails, or the value moves.

| plant | label61 n=1 | control (a) [f_2] | label61 n=3 | control (a') [f_4] | 2ppi value | 3ddel value | in 9-span |
|---|---|---|---|---|---|---|---|
| none (baseline) | ok | ok | ok | ok | memo value | memo value | yes |
| m^2 -> m | ok | **ok (misses)** | caught | caught | unchanged | -0.2376 | no |
| drop f_{m-3}'' | ok | **ok (misses)** | caught | caught | unchanged | 0.1057 | no |
| f_{m-3} coef 12 -> 6 | ok | **ok (misses)** | caught | caught | unchanged | 0.0535 | no |
| S_m -> -S_m | caught | caught | caught | caught | -0.1155 | -0.917 | no |
| c -> -c | **ok (misses)** | caught | **ok (misses)** | caught | 0.863 | -82.9 | no |
| e dropped | **ok** | caught | **ok** | caught | 0.00227 | 0.156 | no |
| Sig -> -Sig | caught | caught | caught | caught | -0.79 | 11.3 | no |
| d^4 term dropped | caught | caught | caught | caught | 0.308 | 5.59 | no |
| c-jets with \|U\| >= 4 dropped | ok | caught | ok | caught | **unchanged** | -7.49 | no |

Readings:
- **The memo's guard list is wrong on one item.** It says control (a) over Q(t) rejects "m^2 -> m". Control (a)
  exercises only m = 1, where m^2 = m, so it cannot. The same holds for every f_{m-3} plant. Only the rung-up control
  (a') (nabla^2 r12^3 = 12 r12, f_4 jets) or label 61 at n = 3 catches them. The memo's line "This control cannot
  pass with ... a wrong m-dependence" (control (a) section) is false for control (a) as stated; it is true for (a)
  together with (a'). The mathematics is unaffected, since (a') was run and passed.
- The label-61 checks see only the base jet, where c = e = 0, so they miss the c and e plants. Control (a) catches
  both.
- The 2ppi value is blind to all the m >= 3 plants and to high-|U| c-jet errors. A guard that pins only 2ppi does not
  protect 3ddel; the memo correctly pins both.

## Independent numbers

| quantity | value |
|---|---|
| 2ppi (2/3,7/4), R = 1.9, closed form (dps 60) | -0.0051596887888591130377079246399078133586026001414902 |
| same, exchange_hp tau 16 dps 32 | -0.00515968878885911303770792463991 |
| 3ddel (2/3,7/4), R = 1.9, closed form | 0.00035589362716950157926422109921378413308955590832378 |
| same, exchange_hp | 0.000355893627169501579264221099214 |
| (2p_x 2p_x\|2p_x 2p_x) (1,2), 1.4, 16-term vs pppi0 + ppdel/2 | 0.15072359406719310298797 vs ...281957 (tau 8) |
| (1,2) baseline values reproduced by my own driver | 2ppi, 3ddel at R = 1.4: identical to the memo's 50 digits |
| controls (a), (a'), label 61 n=1,3 at (2/3,7/4) | all hold exactly; 0 free unknowns |
| committed operators on rebuilt K | each annihilates its own family exactly, not the other |

## Caveats required in a paper sentence

1. Name the established level: the bound (order <= 9, span of the nine functions) for every 0 < a < b, plus
   operator equality with the committed (9,25)/(9,35) **at (1,2) only**. "Exactly 9" and minimality for general
   a < b remain unproved, because independence of the nine functions over C(R) is not shown. This blocks any
   sentence that says "order exactly 9 for all a < b".
2. The equation (**) is derived from Pachucki's label 39 by termwise inverse Laplace transform. That step is
   standard, and its consequences are certified by label 61 on the locus, exact Laplacian controls off it, and
   high-precision agreement with the independent Neumann route. It is not a line-by-line theorem in the source, so
   cite it as derived here from label 39, not as Pachucki's equation.
3. The 2p_x 2p_x statement covers the sigma = 2 part by identity with 3ddel; the m = 0 part rests on the earlier
   extension result.

## Wrong answers a regression guard should reject

- **m^2 -> m, a dropped f_{m-3}'' term, or 12 -> 6:** these are caught only by the f_4 rung (control (a')) or by a
  3ddel value. Control (a) at f_2 and any 2ppi pin pass all three. The guard must include (a') or a 3ddel pin.
- **c -> -c, e dropped, high-|U| c-jet terms dropped:** label 61 misses these (base jet only). Control (a) catches them.
- **S_m -> -S_m, Sig -> -Sig, d^4 dropped:** everything catches these.
- **Wrong family operator:** OP_3DDEL_1_2 applied to 2ppi, or the reverse, must leave a nonzero residual (13 keys).
- **A non-unique solve:** assert `last_free == 0` after each `solve_D2`. The current code sets free unknowns to 0
  silently. That is harmless today (0 free in every solve), but a future case could pass a residual check with an
  arbitrary solution.
- **Values at a pair not used in building:** (2/3,7/4), R = 1.9 above, which tests the a-power and the b-dependence
  together.

## Files

`debug/pachucki_m_review_plants.py`, `debug/pachucki_m_review_hp.py`, `debug/pachucki_m_review_qt.py`,
`debug/pachucki_m_review_gcd.py`, `debug/pachucki_m_review_annihilate.py`, `debug/pachucki_m_review_f0coef.py`.
Data: `debug/data/pachucki_m_review_{plants_1_2,plants_2-3_7-4,hp_2ppi_2-3_7-4_1.9_T16,hp_3ddel_2-3_7-4_1.9_T16,hp_xx_1_2_1.4_T8,qt_3ddel_7-3,qt_2ppi_11-4,gcd,annihilate,f0coef}.json`
and the logs `debug/data/pachucki_m_review_hp_*.log`.
