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
