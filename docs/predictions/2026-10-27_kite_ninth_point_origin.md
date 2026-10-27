# Pre-registered prediction: where the kite's ninth singular point comes from, and its massive branch (ND1-ND4, NM0-NM4)

**Registered:** 2026-10-27, before any computation that scores a prediction in this file. PI-directed follow-up
(PM brief, "Two follow-ups on the kite's ninth singular point, under pre-registration").
**Prompted by:** `2026-10-26_kite_ninth_point.md` (NP1-NP4 and its Disposition: the closed form is fitted, not
derived; NP2 void because both positive controls failed).
**Extends:** `2026-10-26_kite_ninth_point.md`, `2026-10-26_parametric_kite_operator.md`,
`debug/lit_scan/topology_rung_ninth_point_memo.md`.
**Commit status:** committed on its own, before any computation in this file.

## The questions

Write A = a^2, B = b^2, G1 = g1^2, G2 = g2^2, u = A - G1, v = B - G2. The measured, fitted result is

  z9 = (A G2 - B G1)(A - B - G1 + G2) / ((A - G1)(B - G2)),

where the Lee-Pomeransky critical-point count of the kite drops from 22 to 21.

- **Question D (derive it).** Which algebraic condition does z9 solve? Where do the two numerator factors and the
  denominator come from?
- **Question M (the massive branch; replaces NP2).** Give the Coulomb line mass^2 eps. Map the locus in the (z, eps)
  plane where the critical-point count drops, near (z9, 0), directly from the count, not from eps-poles of rebuilt
  functions. Does a drop curve pass through (z9, 0), and what is its local shape?

Kite conventions as in `debug/lit_scan_euler_kite.py`: lines 1 (k, a), 2 (w-k, b), 3 (k-q, g1), 4 (w-k+q, g2),
5 (q, m5), Euclidean, z = w^2, eps = m5^2, G = U + F.

## What I know before computing (disclosed)

- **A hand model, made while designing this record.** Put all five lines on shell, with q^2 = -eps. Then
  k.q = -(u + eps)/2 and (w-k).q = (v + eps)/2 hold identically. In D = 3 the three momenta k, p = w - k and q
  are generically independent. They are coplanar exactly when their Gram determinant vanishes:

    det [[-A, K, c1], [K, -B, c2], [c1, c2, -eps]] = 0,   K = k.p,  c1 = -(u+eps)/2,  c2 = (v+eps)/2,  z = 2K - A - B.

  - At eps = 0 this equation is linear in K. Its root gives z = B u/v + A v/u - A - B, which is the measured
    z9.
  - With r = u/v it factors as z9 = (r - 1)(B r - A)/r. The numerator factor A - B - G1 + G2 is r = 1, which
    means w.q = 0. The factor A G2 - B G1 is r = A/B. Both put the point at z = 0. The denominator is r = 0 or
    r = infinity. There, k.q = 0 or p.q = 0 would be needed with q null and in the plane. That is impossible
    for nonzero q unless the line is massless, so the configuration runs off to infinite K.
  - This Gram determinant is the maximal-cut Baikov polynomial of the kite. The kite has five propagators, no
    irreducible numerators, L = 2 and E = 1. With unit indices the maximal cut goes like
    P_B^((D-L-E-1)/2) = P_B^(-1/2) in D = 3. Each index above 1 lowers that power by one:
    -1/2 - (1/2 + 1/2 + 1 + 1) = -7/2. That is the exponent measured at P2. This is a postdiction, not a scored
    item.
  - Why the parametric Landau test missed it. When U(alpha) != 0, the loop equations fix k and q as multiples
    of w, so the momenta are collinear. A coplanar, non-collinear solution needs a kernel of the loop matrix,
    so U(alpha) = 0. `pko_landau_leading.py` saturated by U, so it excluded such solutions by construction.
- **Checked by computer in phase 1** (`debug/nq_model_numbers.py`, log `debug/data/nq_model_numbers.log`). This
  script evaluates the model only: no reduction, no Groebner basis of the kite, no count.
  - The model's eps = 0 root equals the measured formula, as an identity in the rates. So the identity
    "Gram-on-shell = z9" is the hypothesis this record tests. It is not one of its predictions.
  - The model's curve at P2: C(z, eps) = -1296 eps^2 z + 1296 eps z^2 + 20484 eps z + 2700 eps - 28000 z - 2775 = 0.
    It is a cubic curve through (z9, 0) = (-111/1120, 0).
  - Slope at that point: dz/deps = kappa = 53517321/2195200000 = 0.0243792460824.
  - Model roots at P2:

    | eps | near root | far root |
    |:--|:--|:--|
    | 1/400 | -0.0990460557660772 | 8626.27129914219 |
    | 1/1600 | -0.0990918971533117 | 34552.1953959095 |
    | 1/6400 | -0.0991033330579241 | 138255.898642299 |

  - Threshold control values: -(c +- m)^2 for c in {11/3, -5/3} and m = 1/20, 1/40, 1/80, listed in the log.
- **Cost probe** (`debug/nq_probe_cost.py`, logs `debug/data/nq_probe_cost_*.log`). Four single-process runs at
  once, at random residues z and eps only, at P2, with p = 2^61 - 1 and p = 32003.
  - Massless: count 22, Groebner basis 58-61 s processor time. The six traces of multiplication by t and by
    x1..x5 take 2 s more.
  - Massive (random eps): count 30, Groebner basis 100-101 s, traces 4.5 s.
  - The prime size does not change the cost.
  - Nothing was evaluated at z9, at a threshold or on the model curve. The trace values at one random point
    carry no information about poles.
- **Already measured and not re-scored:**
  - Counts at P2, eps = 0: 22 at random z; 21 at z9 and at -121/9; 17 at the two-line cuts; 8 at z = 0; 22 at
    -73/9, -277/36 and -100/9.
  - With m5^2 = 1/100, the count is 30 at z9 and at a random z.
  - The face scan of the principal-Landau-determinant type found no nondegenerate face at z9.
  - Parametric leading Landau equations, with U nonzero: no solution at z9.

## Route chosen for Question D, and why

- **Not route (i), the symbolic principal A-determinant or Euler discriminant of G.**
  - At P2 the principal-Landau-determinant face scan already missed z9.
  - The symbolic principal A-determinant over 189 faces, with four rate parameters and z, is out of reach of
    sympy.
  - Helmer-Papathanasiou-Tellander's Whitney stratification needs Macaulay2, which is unavailable.
  - What route (i) can still give is a hypothesis-free list of drop points. I get that from the count itself, as
    poles of traces: ND0 below, and Question M.
- **Route (ii), done in momentum space: a leading Landau singularity whose multipliers have U = 0.** This is the
  hand model above: all five lines on shell, the three momenta coplanar and not collinear. It is first-type in
  momentum space and invisible to the parametric equations with U nonzero. The scored items test it as the
  mechanism, by routes independent of the Gram identity:
  - an elimination in components;
  - the Landau multipliers;
  - the escape direction of the lost critical point;
  - the massive continuation.

## Methods

**Trace-pole scan (hypothesis-free drop locus in z at fixed eps).**
- At fixed eps, take random z mod p. Compute the Groebner basis of the critical-point ideal: the five equations
  nu_i G - s x_i dG/dx_i, plus t x1...x5 G - 1, in grevlex order, with the same nu and s as
  `lit_scan_euler_kite.py`.
- From the basis, form the traces T_v(z) of multiplication by v in {t, x1, ..., x5} on the quotient.
- T_v is a rational function of z. It has a pole exactly where a critical point leaves the very affine variety:
  x_i -> 0 or infinity, or G -> 0. Those are the points where the count drops.
- Each T_v is rebuilt with `fit_rational` (the copy in `ib_kite_de.py`), with 6 held-out samples.
- D_eps(z) = the monic least common multiple of the six denominators. Its roots are the drop points.
- Lift: coefficients are lifted to Q by rational reconstruction from p0 = 2^61 - 1 and accepted only if 6 samples
  at p1 = 2305843009213693921 agree. D_eps is then factored over Q.

**Direct counts.** `lit_scan_euler_kite.py` logic, mod p0, at chosen (z, eps). It is used for controls and for
the drop size.

**Components elimination (ND1).**
- Set w = (omega, 0), k = (k1, k2), q = (q1, q2), z = omega^2. Impose the five on-shell equations with eps = 0.
- Compute a sympy lex Groebner basis that eliminates k1, k2, q1, q2, omega, saturated by k2 (non-collinear).
- Run it at P2 and at two random rational rate points R1, R2, chosen by a seeded generator before running.
- No Gram identity is used.

**Landau multipliers (ND2).**
- At P2, take the exact coplanar on-shell configuration from ND1.
- Solve the two vector loop equations for (alpha_1, ..., alpha_5):
  - alpha1 k - alpha2 (w-k) + alpha3 (k-q) - alpha4 (w-k+q) = 0
  - -alpha3 (k-q) + alpha4 (w-k+q) + alpha5 q = 0
- These are four scalar linear equations, solved exactly over the field of the configuration.

**Escape direction (ND4).**
- Numerical critical points at z = z9 + delta, delta in {1e-2, 1e-3, 1e-4}: Newton in log coordinates from random
  starts, mpmath at 40 digits, deduplicated.
- Scored only if exactly 22 distinct points are found at every delta.
- The escaping point is the one whose norm or projective position moves without bound as delta -> 0.
- Two (nu, s) sets are used.

## Predictions

**ND0 (complete drop locus at P2, massless; also the method control).**
- D_0(z) has exactly the roots 0, z9 = -111/1120, and the eight cut values -25/4, -1/4, -25, -1/9, -121/9, -25/9,
  -529/36, -25/36. Nothing else.
- In particular none of -100/9, -16/9, -625/36, -49/36, -73/9, -277/36.
- *Method gate:* if any of the eight cut values or z9 is missing, the trace-pole method is void, and NM1-NM3 are
  not scored.
- Falsifier: an extra root, which would be a drop point that no known mechanism names. Or a missing root, which
  triggers the gate.
- Probability 0.6 for "exactly these ten". The gate is expected to pass, at 0.9.

**ND1 (the elimination in components).**
- At P2, R1 and R2, the z-eliminant of the non-collinear on-shell system is the model factor
  z(A-G1)(B-G2) - (AG2-BG1)(A-B-G1+G2), evaluated at those rates, times at most a power of the {1,2} cut factor
  (z + A + B)^2 - 4AB. There is no other factor.
- Falsifier: Delta9 absent, or another factor present.
- Probability 0.85. This is mostly a check of the hand algebra by an independent route.

**ND2 (it is a genuine pinch with U = 0).**
- At P2 on the coplanar solution, the loop equations have exactly a one-dimensional solution space alpha*.
- All five components of alpha* are nonzero and finite, and U(alpha*) = 0. So it is a Landau singularity with all
  five lines on shell, not a contraction, and it is excluded by U != 0.
- Reported, not scored: F(alpha*), and whether alpha* is real.
- Falsifier: no nonzero solution, a two-dimensional solution space, some alpha_i = 0, or U(alpha*) != 0.
- Probability 0.7. U = 0 is forced if the loop equations hold with k not parallel to w; the risk is a vanishing
  alpha.

**ND3 (the numerator and denominator are not components; they are where z9 meets z = 0 or z = infinity).**
Direct counts at three rate points, given in squared rates:
- R_a = (1, 9/4, 49/9, 49/4), on AG2 = BG1, where z9 = 0.
- R_b = (1, 9/4, 49/9, 241/36), on A - B - G1 + G2 = 0, where z9 = 0. Only squared masses enter G.
- R_c = (1, 9/4, 1, 64/9), on A = G1, where z9 = infinity.

Predictions:
- (a) At z = -37/101 the count is 22 at all three. None of the factors is a component of the Euler discriminant on
  its own.
- (b) At R_a and at R_b the count at z = 0 is at most 7, below P2's 8 at z = 0, because two drop components meet
  there.
- Falsifiers: (a) any count other than 22; (b) a count of 8 or more at z = 0.
- Probabilities: (a) 0.65. The denominator locus could itself change the number of masters. (b) 0.5.

**ND4 (the count drop is this pinch: the escape direction).**
- As delta -> 0, exactly one of the 22 critical points escapes.
- Its projective direction tends to alpha* from ND2, the same for both (nu, s).
- The projective distance falls at least tenfold from delta = 1e-2 to 1e-4.
- Reported, not scored: whether |x| -> infinity or -> 0, and the rate exponent.
- Falsifier: no escaping point, more than one, or an escape direction away from alpha*.
- Not scored if fewer or more than 22 points are found at any delta, or if ND2 fails (no alpha*).
- Probability 0.5.

**NM0 (positive controls; run FIRST; Question M is void if any fails).** Direct counts at eps = 1/400 (m = 1/20):
- At z = -(11/3 + 1/20)^2 = -49729/3600 and at z = -(5/3 - 1/20)^2 = -9409/3600, the count is below 30.
- At each of those values + 1/1000, and at z = -37/101, it is 30.
- In the trace scans, D_eps(z) must contain the four roots -(11/3 +- m)^2 and -(5/3 +- m)^2 at eps = m^2,
  for every eps scanned.
- If any of these fails, NM1-NM4 are not scored.
- Expected to pass, 0.85. This is a check of the method, not a prediction about z9.

**NM1 (a drop curve passes through (z9, 0), and it is the coplanar pinch).**
- At eps = 1/400, 1/1600 and 1/6400, the model quadratic Q_eps(z) = 1296 eps z^2 + (20484 eps - 1296 eps^2 - 28000) z
  + 2700 eps - 2775 divides D_eps(z) at p0 and at p1.
- Negative control: Q_eps with eps replaced by eps(1 + 1/1000) must not divide.
- Falsifier: Q_eps does not divide at some eps.
- Probability 0.65.

**NM2 (local shape: transverse and analytic in eps = m5^2).**
- In the lifted D_eps, the only root within |z - z9| < 0.01 is the model's near root in the table above.
- The roots move linearly: (z_j - z9)/eps_j = 0.02443, 0.02439, 0.02438, tending to kappa = 0.0243792...
- So the curve is eps = (z - z9)/kappa + O((z - z9)^2), transverse to the axis eps = 0, and analytic in m5^2.
  It is not the tangent, sqrt-type shape of the threshold controls, eps ~ (z - z_c)^2/(4c^2).
- Falsifier:
  - a root in the window other than the model's;
  - no root in the window;
  - (z_1 - z9)/(z_2 - z9) near 2 (sqrt scaling) instead of near 4.
- The window clause is not scored if the lift to Q fails. Then only NM1's divisibility stands.
- Probability 0.6.

**NM3 (the second branch goes to infinity).**
- D_eps also has the model's far root, about 8626.27, 34552.2 and 138255.9.
- It escapes as z ~ uv/eps = (1750/81)/eps when eps -> 0. That is why only one drop point survives at eps = 0.
- Falsifier: far root absent. If NM1 holds this is automatic, so it is scored with NM1 but reported separately.
- Probability 0.65.

**NM4 (drop size).**
- At a point of the model curve mod p0 (eps = 1/400, z a root of Q_eps in F_p0, or in a prime from
  `ib_primes.txt` where the discriminant is a square), the count is 29.
- At a random z with the same eps it is 30.
- Falsifier: a count of 30, or a drop by more than 1.
- Probability 0.55 for exactly 29. A drop of 1 or more is expected at 0.65.

## Novelty status

- **Leading singularities and symbol letters from Gram determinants on cuts in the Baikov representation:
  KNOWN.**
  - Jiang, Liu, Xu and Yang, "Symbol letters of Feynman integrals from Gram determinants", arXiv:2401.07632.
  - "Alphabet of one-loop Feynman integrals", arXiv:2201.12998.
  - Dlapa, Li and Zhang, leading-singularity protocol in Baikov representation, arXiv:2103.04638.
  - These three are [V] at search-listing level only, not abstract pages. To be checked against the abstract
    pages before any disposition cites them.
  - Baikov-representation cuts: Frellesvig and Papadopoulos, arXiv:1701.07356 [UNV].
- **Landau equations in momentum space vs parametric form, and Landau components missed by principal Landau
  determinants:** as in the NP record's novelty section (Fevola-Mizera-Telen; Helmer-Papathanasiou-Tellander;
  Klausen; Mizera-Telen; Vergu).
- **If ND1, ND2 and NM1 hold, the ninth point is a standard object:** the zero of the kite's maximal-cut Baikov
  polynomial in D = 3, a coplanar pinch. Its novelty is then at most that this D = 3 kite's point was first seen as
  "unpredicted" by this program's own cut rule.
  - The derivation would be a CHECK of known machinery, not a finding.
  - The statement "it belongs to the massless line" in Paper 59 would need restating. The pinch exists at every
    eps. At eps = 0 its second branch leaves to infinity.
- No search beyond the above was made, so nothing here is a novelty claim.

## Controls

- **ND0 gate:** the eight cut values and z9 must be seen by the trace-pole scan at eps = 0.
- **NM0 gate:** the threshold branches through -121/9 and -25/9 must be seen, by direct count and in D_eps, before
  anything near z9 is read. The scan script prints the control verdicts before it prints any other root.
- **Fit discipline:** every T_v passes its 6 held-out z. Every lift is confirmed at p1.
- **NM1 negative control:** the perturbed quadratic must not divide.
- **ND3 reference:** the P2 counts already logged (22 at -37/101, 8 at 0).
- **ND4 reproducibility:** two (nu, s) sets must agree on the escape direction.

## Caps stated in advance

- One graph, D = 3, one rate point (P2) for the massive branch.
- Genuineness is by counts modulo a prime and by numerical critical points. It is not rigorous.
- Only the exponent's postdiction is offered for the -7/2. Exponents are not re-measured.
- No closed form is proposed beyond the model already written above. No small-integer searches.

## Compute budget and STOP rules

- **Processes.** At most 8 of mine at any time: Pool(7) plus the parent for scans; at most 6 direct counts at once
  when no pool runs. Processes I did not start are never touched. Every multiprocessing script has an
  `if __name__ == "__main__":` guard.
- **Per-sample cost estimates** (probe, 4 concurrent): 63 s processor time massless, 105 s massive. Allow 1.5x
  under 8 concurrent.
  - STOP rule: if the median of the first 7 samples of a scan exceeds 180 s, STOP that scan.
- **Trace scans:** eps = 0, then 1/400, 1/1600, 1/6400 at p0. For each:
  - 140 samples to start (fit total degree up to 132, with 6 held out). One extension to 220 if no fit; then STOP
    that eps (not scored).
  - 6 samples at p1 for the lift.
  - At most 4 x 226 = 904 samples, about 2.8 h wall at the estimated rate.
- **Direct counts** (NM0, NM4, ND3): at most 16, at most 40 min wall.
- **ND1 / ND2 sympy eliminations:** at most 30 min each. If one has not finished, STOP it; the item is not scored.
- **ND4 numerics:** at most 45 min wall.
- **Whole phase:** at most 8 h wall. Exceeding it is a STOP, reported with what was reached.
- **Output files:** new names only, `debug/nq_*.py` and `debug/data/nq_*`. No saved file is overwritten. No papers,
  CHANGELOG or CLAUDE.md edits in this record's runs.

## Decision gate

- **GO (derived):** ND1, ND2 and NM1 hold. z9 is then derived as the coplanar (Baikov-polynomial) pinch, and the
  massive branch is identified.
- **BORDERLINE:** ND1 holds but ND2 or NM1 fails or is void. The identity holds, but the mechanism is not
  established.
- **STOP:** ND0 gate fails, or the caps are exceeded before ND1 and the scans are read.

The other items (ND0, ND3, ND4, NM2-NM4) are scored separately and do not change the gate.

## Risks

- The traces could have a cancelling pole: two escaping points whose divergent parts cancel in the sum. Six
  different traces make this unlikely; ND0's gate tests it on known points.
- Rebuilt degree in z may be larger than 132. The extension and STOP rule handle it.
- At eps small, the near root sits close to the eps = 0 values. That does not matter mod p. The window clause uses
  the lift over Q.
- ND3's rate points have their own coincidences. For example, R_c makes the pairing {1,3} equal-mass, which
  may change masters. That is exactly what (a) tests.
- ND4 numerics may not find all 22 points. In that case it is not scored, by rule.
