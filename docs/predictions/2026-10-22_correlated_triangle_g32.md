# Pre-registered prediction: the correlated triangle at gamma = 3/2 (GT1-GT3)

**Registered:** 2026-10-22, before any computation in this file (PI-directed: "Let's do those next steps",
item three of the 2026-10-21 overview).
**Prompted by:** Paper 59 `sec:density_ring`: "At gamma = 3/2 a quotient of order three, found modulo a
prime, has the local data of F_3's operator twisted by t ... the explicit map was not found, and a factor
above K3 is not excluded there." This record tests the second clause and tries to close the first.
**Extends:** `2026-10-17_ring_extensions.md` (QG), `2026-10-17_density_loop.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

Y(t) = <rho_i rho_j rho_k f(r_ik) f(r_jk) / r_ij>, f = exp(-gamma r), gamma = 3/2, fingerprint fibre,
scales (1, 1, t), t on the apex k, as in QG. Does Y's annihilator climb above the K3 rung at gamma = 3/2?

## What I know before computing (disclosed)

From QG (`debug/sprint_ring_ext_QG_memo.md`, data `debug/data/ring_ext_QG_6g_L9_g1.5*.json`):
- Y_log has an exact integer operator L9 of order 9 and t^2-degree 25; no order-8 operator up to degree
  43. L9 annihilates Y itself to 7e-28 at t = 0.8, 1.3, 1.7 with no source.
- Singular points t = 1/2, 2, 5/2, 3 (|1 - gamma|, 2, 1 + gamma, 2 gamma), plus 0 and infinity, plus a
  degree-19 apparent polynomial in t^2.
- Exact: L9 = M7 o R2, with R2 of order two and reducible (a rational and an algebraic solution).
- Mod p: an order-two right factor through the t = 3 branch (exponent -3/2), not identified, with no
  order-one factor through that solution; an order-three quotient N3 of L9* with F_3's local data twisted
  by t (maximally unipotent at 0, a pseudo-reflection at t = 2, no other true singular points). The map to
  F_3 was not found with denominators (apparent)^(1,2) (t^2-4)^(0,1) and numerator degree up to 146;
  the apparent factor of N3 has degree 66 (10 at gamma = 1).
- At gamma = 1 the composition factors are algebraic pieces and one K3 piece, nothing else.
- The kite phi has order 6 at gamma = 3/2 (2 at gamma = 1).
- Nothing in this file has been computed.

## Model

Start from L9 as stored. Exact right-factor peeling over Q(t), with mod-p searches for right factors
(local-exponent bounds, multi-prime CRT, rational reconstruction) and an exact division check before any
factor is accepted. Irreducibility of a piece of order two or three by local exponents and Fuchs'
relation where it applies, otherwise by an exact or mod-p search within the leading-term bound, with a
planted split as control. Numerical checks on Y by the QG quadrature at 30 digits or more.

## Predictions

**GT1 (no climb).** L9 has a composition series over Q(t) in which every factor has order at most three,
and exactly one factor has order three; that factor has N3's local data.
- Pass: the series exhibited, each division exact.
- Falsifier: an exactly exhibited piece of order four or more with no right factor of lower order within
  the leading-term bound (mod p, with the planted split found by the same call). That is a climb.
- Anything else (series incomplete, no falsifier): UNCLEAR.

**GT2 (the K3 map).** An exact map, rational in t, sends the solutions of F_3's operator, twisted by t as
at gamma = 1, into the solutions of the order-three factor.
- Pass: the map exhibited and the image annihilated exactly, modulo F_3's operator.
- Not found within the searched box: UNCLEAR, not a falsifier, since rigidity says the map exists if the
  local data are right.
- Falsifier: the order-three factor's local data contradict F_3's at some point (an exponent difference
  at t = 2, 0 or infinity that F_3's twisted operator cannot have).

**GT3 (Y needs it).** The right factor of L9 below the order-three piece does not annihilate Y: applied
to Y it is nonzero at two values of t, by more than 1e-5 of the largest term.
- Falsifier: it annihilates Y to 1e-20 (then Y sits below K3 at gamma = 3/2).

Not scored: the identification of any irreducible order-two factor (elliptic or otherwise) is reported
as exploratory.

## Novelty status

- GT1, GT2: NOT FOUND, as for QG. A hit on both is a **FINDING** of moderate weight: it extends the
  gamma = 1 placement, which is already a FINDING, to a second gamma with moving singular points.
- GT3: a CHECK.

## Controls

- Positive: the same pipeline at gamma = 1 recovers G4 = R3 R1 and the exact F_3 map stored in
  `debug/data/ring_ext_QG_7_k3map_g1.json`.
- Positive: the split search finds a planted split built from L9's known pieces.
- Negative: a 1/200 perturbation of one coefficient of the GT3 operator spoils the GT3 residual pattern.

## Caps stated in advance

One gamma value beyond 1, the fingerprint fibre, the exponential geminal, apex placement. L9's minimality
is bounded (no order 8 up to degree 43), not proved. Mod-p searches are bounded by their boxes.

## Feasibility rule

If the factorization or the map search cannot be completed within the session's memory, the computation
stops and the record says STOP for feasibility, unscored, with the stage reached.

## Disposition

**Computed:** 2026-10-22, after the registration commit 4b1be301. Scripts `debug/correlated_triangle_g32_*.py`, data `debug/data/correlated_triangle_g32_*` (exact pieces in `..._series_final.json`, the map in `..._8_k3map.json`). Nothing committed. The feasibility rule did not fire.

**GT1: PASS (FINDING, with GT2).** An exact composition series of L9 over Q(t). Each inclusion is an exact right division.
- Chain: Sol(ra) < Sol(R2) < Sol(V4) < Sol(V5) < Sol(K6) < Sol(L9), orders 1, 2, 3, 4, 6, 9.
- Factors: four of order one, one of order two (C2b), one of order three (F). The order-one factors are algebraic: R2's rational solution; exponent -5/2 at t = 2; exponent -7/2 at t = 3; t^5 (t^2-9)^(-9/2) times a sextic in t^2.
- K6 was found mod p per prime from N3 (the L9* right factor, x-degree 34), then 16-prime CRT. K6 has x-degree 19. L9 = F o K6 exactly.
- F has N3's local data. MUM at t = 0 (exponents 1, 5, 9; one log-free solution). Pseudo-reflection at t = 2 (exponents -11/2, 0, 1). One Jordan block at infinity. Regular at t = 1/2, 5/2, 3. Apparent part of x-degree 19 + 2 x 14.
- Irreducibility: F has no right factor of order one or two. C2b has none of order one. Method: exponent-sum bounds, valid because each quotient's solutions are analytic away from L9's true singular points; every admissible exponent choice was tested mod p, nullity 0.
- QG's piece (b) is reducible. LCLM(R2, B2) has order 3: B2 contains R2's algebraic solution.
- The order-three piece K6/V4 also splits, as C2b o rho1, found by the same search.

**GT2: PASS (FINDING, with GT1).** A = t (x-4)^(-4) sum_j N_j(x) theta^j, x = t^2, deg N_j <= 20. F o A = 0 modulo F_3's operator, exactly (8 primes). The denominator has no apparent factor. Box: t^s with s odd in [-11, 1], (x-4)^(-b) with b <= 6, numerator x-degree <= 80; the hit came at 24 and was minimised to 20. This box differs from QG's (apparent factors, (t^2-4)^(0,1)). The map was sought into the left quotient F, whose solutions are analytic at its apparent points, so no apparent factor is needed. This widens the (t^2-4) power from 1 to 6.

**GT3: PASS (CHECK).** K6 Y relative to the largest group: 1.46, 2.43, 1.31 at t = 0.8, 1.3, 1.7. Relative to the largest monomial: 0.27, 0.78, 1.21. Reference: L9 Y is 7e-28, 8e-29, 2e-28. Quadrature h against 2h on K6 Y: 2e-16 to 7e-15.

**Controls.**
- Positive, gamma = 1: PASS. R1 recovered. G4 = F o R1 exactly. F equals QG's R3 up to the gauge (t^2-4)^2 (9t^4+312t^2-12560). The stored map was recovered exactly (constant ratio). A map into F needs only t (x-4)^3 in the denominator (numerator degree 7), exact.
- Positive, planted split: PASS. On V4 = B2q o rb o ra the search finds the order-one factors and the order-two factor R2. On B2 it finds R2's algebraic solution. The order-two test first failed this control (a missing 1/c_r in the adjoint-Wronskian relation). It was fixed and rerun on every piece.
- Negative, 1/200: NOT MET as worded. The largest theta^0 scalar coefficient of K6, times 1 + 1/200, moves K6 Y by 0.37 to 0.64 percent. The pattern (nonzero, order one) is unchanged; a nonzero residual cannot be spoiled. Supplementary, not registered: the same perturbation of L9 moves 7e-28 to between 3e-4 and 7e-3.

**Not scored.** C2b is an irreducible order-two factor: unipotent at t = 0, +-1/2, +-5/2 and infinity, regular at +-2 and +-3. It is not identified. At gamma = 1 there was no such factor.

**Caps.** L9's minimality is bounded, not proved. This does not affect GT3: Y's own minimal operator is a right factor of L9 not contained in K6, so it carries F. The irreducibility tests used one prime (2^61-1); they rest on Fuchsian regularity and the standard exponent facts. One gamma, the fingerprint fibre, apex placement.

**Deviations.** The order-three factor tested is F, the exact left quotient, not N3; they are the same composition factor. Irreducibility used exponent sums, not a leading-term bound. The negative control did not discriminate.

**PM check (2026-10-22).** Rerun with the sympy helpers of `tests/test_paper59_generic_slice.py`, a code
path independent of the agent's flint arithmetic (log `debug/data/correlated_triangle_g32_pm_check.log`):
L9 / K6 has zero remainder and the quotient is F up to a rational left factor; each step of the chain
ra, R2, V4, V5, K6 right-divides the next exactly; F o A = 0 modulo F_3's operator exactly; a 1/1000
change in one coefficient of A breaks it. A is injective on Sol(F_3), since F_3's operator is irreducible
and A is a nonzero operator of order two. Scored in the record's words: GT1 HELD, GT2 HELD (FINDING
together), GT3 HELD (CHECK). The negative control was ill-posed at registration, not by the computation:
it perturbs an operator whose residual on Y is meant to be nonzero, so it cannot discriminate. The
unscored C2b is unipotent at exactly four points in x = t^2 (0, 1/4, 25/4, infinity), the pattern of a
Beauville-type elliptic family; that is a lead, not a result.
