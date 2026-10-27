# Pre-registered prediction: where the kite's ninth singular point comes from, as a function of the rates (NP1-NP4)

**Registered:** 2026-10-26, before any computation that scores a prediction in this file. PI direction ("Let's do
those things"), after the literature scan.
**Prompted by:** `2026-10-26_parametric_kite_operator.md` (PK1-PK8, both Dispositions) and
`debug/lit_scan/topology_rung_ninth_point_memo.md` (Q2).
**Extends:** `2026-10-26_parametric_kite_operator.md`.
**Commit status:** committed on its own, before any computation in this file.

## The question

At the rates P2 = (a, b, g1, g2) = (1, 3/2, 7/3, 8/3) the kite's order-11 operator (`debug/data/ib_pko_p2_op_Q.json`)
has a genuine simple singular point z9 = -111/1120 that none of its cuts predicts. Its exponents are -7/2, 0, 1, ..., 9.
The leading Landau equations of the kite and its contractions have no solution there. The critical-point count of the
Lee-Pomeransky likelihood drops from 22 to 21 at z9, in the top sector only, and does not drop there once the Coulomb
line has mass^2 1/100. Three questions:

1. How does z9 move along the family (a, b, g1) = (1, 2, 3/2), s = g2^2 free? Does it reach z = 0 at g2 = 3, the
   base point?
2. Is z9 the limit of a singular point of the kite with a massive Coulomb line as its mass m5 goes to 0? Or does it
   exist only at m5 = 0?
3. Is there a closed form for z9 in the rates? It counts only if exact reconstruction fixes it from more points than
   it has free coefficients, and it then passes held-out points.

## What I know before computing (disclosed)

Facts read from existing data or from one cost probe run for this record. None of them is a scored result.

- **Where z9 sits in the P2 system** (read from `debug/data/ib_pko_p2_p0_op.json`, prime p0 = 2^61 - 1, during the
  design of this record). The 18 masters are the same labels, in the same order, as at the base point. The common
  denominator of A has degree 12 with simple linear factors, one of them z - z9. In c, z - z9 has power 3. 32
  entries of (c, A) have z - z9 in their denominator. The cheapest ones are in the row of the top master
  (0, 0, 0, 0, 1). The entries A[(0,0,0,0,1)][j] for j = 6, 11, 13, 15, 16 are a constant over (z - z9), with total
  z-degree 1. Entry j = 17 has numerator degree 1 and denominator degree 2. This is what makes a cheap readout
  possible: at a rate point, a few reduction samples at random z fix the pole of those entries.
- **What I did not look at.** I did not look at the same entries at the base point, at any family member, or at any
  massive point. One script listed the factor multiplicities of the base system's A and c without testing any
  particular entry. They are as already disclosed in the PK record: the common denominator of A has simple linear
  factors, and c has one linear factor to the power 3, which is z^3 there.
- **Massive-Coulomb cost probe** (`debug/np_probe_massive_time.py`, log `debug/data/np_probe_massive_time.log`). One
  reduction sample at P2 with m5^2 a random residue mod p0, at one random z. The zero sectors were changed to "lines
  3, 4 and 5 all absent", since with m5 nonzero the sectors with line 5 alone on the q-loop are no longer scaleless.
  Result: 2226 seeds, 30 masters, 31.0 s processor time. 30 equals the generic critical-point count with
  m5^2 = 1/100 reported in the lit-scan memo.
- **Massless cost**, from the PK record: 27.8 s processor time per sample, alone. Stage 2a ran 35 samples in 297 s
  wall on Pool(6), about 8.5 s wall per sample.
- **Primes**: `debug/data/ib_primes.txt`, all about 2^61. Rational reconstruction from one prime recovers a value
  whose numerator and denominator are below about 10^9.
- **The IBP identities involve only squared masses** (`ib_kite.kite_family`; confirmed in the PK record by the g2^2
  residue probe). So z9 is an algebraic function of (a^2, b^2, g1^2, g2^2) of mass dimension 2. If it is rational
  (see NP3), it is homogeneous of weight 1 in the squared rates.
- Nothing about z9 at any rates other than P2 has been computed. No massive-Coulomb system has been reconstructed.

## Method

**Readout of z9 at one rate point (massless).** Take r = 4 reduction samples at random z mod p. In each of the five
carriers A[(0,0,0,0,1)][j], j in {6, 11, 13, 15, 16}, fit y = C/(z - r0) from two samples. Then check the fit on the
other two samples. All five carriers must agree on r0. Then z9 = r0. If any check fails, take 8 more samples. In that
case rebuild the carriers by `fit_rational`, factor their common denominator, and divide out the known polynomials:
z, the cut quadratics (z + X + Y)^2 - 4XY for the pairs (a,b), (g1,g2), (a,g2), (b,g1) in squared rates, and the
system-only factors z + a^2 + g2^2 and z + b^2 + g1^2. What is left is the z9 factor. At most 12 samples per point.
The master basis is matched by label, not by index. At a rate point with rational rates, z9 is lifted to Q by
rational reconstruction and confirmed on a second prime.

**Pilot (gate for NP4 and the cheap readout).** One random s on the family, 40 samples at p0, the full system rebuilt
in z by the `ib_kite_de.py` fitter (35 to fit, 5 held out). The structure of the carriers is read off. Only the
40-sample pilot is used to score NP4.

**Stage A, the family F1: (a, b, g1) = (1, 2, 3/2), s = g2^2 free.**
- 36 random residues s at p0 to fit z9(s) with `fit_rational`. Its own redundancy is 2 points, so total degree
  (numerator plus denominator) is at most 32. Then 6 random s held out at p0 and 4 random s at a second prime p1.
- The coefficients are lifted to Q by rational reconstruction from p0, with monic denominator. The lift is accepted
  only if it reproduces the 4 points at p1.
- Exact rational held-out points, each measured at p0 and p1: g2 = 5/2, 7/2, 4.
- The base point g2 = 3 is checked directly from the three saved base systems (`ib_pf_p5/p10/p13_op.json`), with no
  new samples.
- Points registered: 36 to fit, 6 + 4 + 3 held out. The formula is accepted only if its fitted total degree t
  satisfies t + 2 <= 34 (at least two fit points to spare) and every held-out point agrees exactly.
- **Genuineness along F1** (independent route): the critical-point count (`debug/lit_scan_euler_kite.py`, mod 32003, or
  65521 if a denominator vanishes mod 32003) at g2 = 5/2 and at g2 = 7/2. At each, count at z9(g2), at
  z9(g2) + 1/1000, and at one random z.

**Stage B, the massive limit, at P2.** Use the probe's reducer, with eps = m5^2 a random residue. Sample at five fixed
z: z9 = -111/1120; two positive controls, -121/9 = -(a+g2)^2 and -25/9 = -(a-g2)^2, the {1,4,5} threshold and
pseudo-threshold, which involve line 5; and two random z, the negative controls. At each z there are 64 random eps:
58 to fit and 6 held out. Every entry of (c, A) is rebuilt as a rational function of eps. For each entry, record
v(z0), the order of vanishing at eps = 0 of its denominator, and the other eps-roots of the denominator. A singular
branch of the massive system passes through (z0, eps = 0) when some entry has v(z0) greater than its v at both
random z. A pre-check: if a sample at z = z9 with eps random gives fewer than 30 masters, z9 is singular for every
eps. That is recorded as its own outcome.

**Stage C, a closed form in all four rates.** It runs only if NP3a holds at the Stage-A pilot, that is if the z9
factor is linear over F_p(s).
- Use the readout at random residues (b^2, g1^2, g2^2) with a = 1.
- Step C1: one random line t -> (b^2, g1^2, g2^2) = u + t v at p0. 36 values of t to fit, 6 held out. This gives the
  total degrees (dn, dd).
- Step C2: a dense ansatz z9 = N/D in the three squared rates, with total degrees dn and dd. The number of unknowns
  is U = C(dn+3, 3) + C(dd+3, 3) - 1 after normalization. Solve a linear system from U + 12 random points. Its
  solution space must be one-dimensional, and all U + 12 equations must hold. Then 12 random points held out at p0
  and 6 at p1.
- Coefficients are lifted to Q by rational reconstruction from p0 and accepted only if the p1 points agree. The
  formula is homogenized with a^2 afterwards.
- Exact held-out checks that it must pass: P2 gives -111/1120, a value measured in the PK record and not used in the
  fit; the F1 formula from Stage A is reproduced identically; and the scaling control below.

## Predictions

**NP1 (the family F1).**
- (a) At generic s on F1, the system has exactly one root outside the known list: the cut quadratics, z, and the two
  system-only factors. It is a simple pole of the five top-row carriers and has power 3 in c, as at P2. It is
  genuine: at g2 = 5/2 and at g2 = 7/2 the critical-point count at z9(g2) is one less than at z9(g2) + 1/1000 and at
  the random z, which agree with each other.
- (b) z9(s) is a rational function of s over Q, reconstructed under the rule above, and passes all 13 held-out
  points.
- (c) z9 reaches z = 0 at the base point: the formula gives z9(9) = 0, and at the base point the five carriers
  A[(0,0,0,0,1)][j] have denominator exactly z (all three saved primes).
- Falsifiers: (a) no extra root, more than one, a different power or form, or no drop of the count (or a drop at
  the control); (b) no reconstruction within total degree 32 at 36 points, or any held-out mismatch; (c) z9(9) is
  not 0, or a base carrier's denominator is not z.
- Confidence: (a) moderate, since F1 keeps the coincidences (a-g1)^2 = (b-g1)^2 and (a-b)^2 = 1, which could alter
  or hide z9; (b) moderate to high, if (a) holds; (c) moderate, about 60 percent. The evidence for (c) is the power 3
  of z at the base against the power 3 of (z - z9) at P2, and the base's unexplained exponent -4 at z = 0.

**NP2 (the massive limit).** z9 exists only at m5 = 0. At z = z9, every entry's denominator has the same order of
vanishing at eps = 0 as at both random z, and no other eps-root of any denominator at z9 tends to 0. The reduction at
z9 with eps random has the generic 30 masters. The positive controls must detect a branch: at -121/9 and at -25/9
some entry's v must exceed its generic value, with the extra eps-roots predicted by the threshold quadratic
(z + A^2 + eps)^2 - 4 A^2 eps, namely eps = 4(a+g2)^2 = 484/9 at -121/9 and eps = 4(a-g2)^2 = 100/9 at -25/9.
- Falsifier: some entry has v(z9) above its generic value (a massive branch ends at z9, so z9 is a continuation), or
  the master count drops at z9 for random eps (z9 does not depend on m5).
- Not scored if a positive control fails. In that case the method cannot see a known branch, and Stage B is reported
  as void.
- Confidence: low to moderate, about 55 percent. For it: the lit-scan locates z9 in a stratum that is degenerate
  because line 5 is massless, and the count at z9 no longer drops at m5^2 = 1/100. Against it: a massive branch that
  moves away with the mass would show exactly that count result too.

**NP3 (closed form).**
- (a) z9 is a rational function of the squared rates. At the pilot its factor is linear over F_p(s), and at every
  Stage-C point the readout succeeds with a linear factor.
- (b) Stage C reconstructs z9 = N/D in (a^2, b^2, g1^2, g2^2), homogeneous of weight 1, within the cap U <= 250, under
  the rule above. It passes all held-out points, reproduces -111/1120 at P2, and restricts to the F1 formula.
- Falsifiers: (a) a factor of degree 2 or more, or a readout failing at more than 5 percent of points; (b) U above
  the cap, or any held-out mismatch, or the P2 value missed.
- Confidence: (a) moderate to high, since z9 is rational at P2 and the factor there is linear over Q; (b) low to
  moderate. The degree is unknown, and U grows as the cube of it.

**NP4 (cheap readout).** At the pilot, the five carriers keep the form C/(z - z9) of total z-degree 1. At Stage-A and
Stage-C points the 4-sample readout succeeds without the 8-sample extension at at least 95 percent of points.
- Falsifier: a different form at the pilot, or more than 5 percent of points needing the extension.
- Confidence: moderate.

## Novelty status

- Reading singular points from the denominators of a differential system, and reconstruction over finite fields: KNOWN,
  standard. The finite-field references are as in the PK record (Peraro arXiv:1608.01902 and 1905.08019; Klappert and
  Lange arXiv:1904.00009; von Manteuffel and Schabinger arXiv:1406.4513). They are recalled and not re-verified here.
- Landau components missed by leading Landau equations and by the principal Landau determinant, found by a drop of the
  Euler characteristic: KNOWN. Helmer, Papathanasiou and Tellander, "Landau Singularities from Whitney
  Stratifications", arXiv:2402.14787. The abstract was fetched for this record: it states that the full singular set
  is read off a Whitney stratification, including strata of kinematic limits such as setting masses to zero. The
  parachute example with the drop of the Euler characteristic from 19 to 18 is as read in full by the lit-scan memo.
  Principal Landau determinant: Fevola, Mizera and Telen, arXiv:2311.16219 and 2311.14669; Mizera and Telen,
  arXiv:2109.08036. Principal A-determinants and second-type singularities: Klausen, arXiv:2109.07584. Massless lines
  make the on-shell spaces singular: Vergu, arXiv:2512.12763. Euler characteristic equals the number of masters:
  Lee and Pomeransky, arXiv:1308.6676; Bitoun, Bogner, Klausen and Panzer, arXiv:1712.09215. All of these are [V] in
  the lit-scan memo at abstract level or better.
- The point z9 of this D = 3 kite, its formula and its massless origin: not found in the scan. No search beyond the
  memo's was made, so this is not a novelty claim. NP2's question, continuation from a massive branch versus birth at
  m5 = 0, is the one HPT's stratification of mass limits addresses in general. Any answer here is a measurement on
  one graph.

## Controls

- **Pilot identification.** The pilot's full system has 18 masters and the same labels as P2.
- **Readout identity at P2.** The readout run fresh at P2 (4 samples, p0) must return -111/1120. This is a positive
  control of the readout code, with a known answer.
- **Scaling control.** The readout at (2, 3, 14/3, 16/3), twice P2, must return 4 x (-111/1120) = -111/280.
- **Negative control, permuted labels, which must fail.** Stage A's fit rerun with the 36 values of z9 assigned to
  randomly permuted s must give no reconstruction within total degree 32, or must fail the 6 held-out points.
- **Negative control, critical-point count.** At z9(g2) + 1/1000 the count must not drop (see NP1a).
- **Stage B positive controls** (-121/9, -25/9) and negative controls (two random z), as in NP2. If a positive
  control fails, Stage B is void.
- **Stage B code control.** The massive reducer run with eps = 0 and the massless zero-sector rule must reproduce
  the saved P2 p0 entries of the top-master row at 2 random z.

## Caps stated in advance

- One graph, D = 3. Not rigorous: genuineness is by the critical-point count modulo a small prime, and along F1 the
  exponents at z9 are not measured.
- The F1 formula is a statement about one line in rate space. A four-rate formula is claimed only if Stage C passes.
- No guessing. No closed form is proposed from small-integer searches, from P2's value, or from any fit that has not
  passed its held-out points.

## Compute budget and STOP rules

- At most 8 processes of mine at any time: reductions on Pool(7) plus the parent; critical-point counts only when no
  pool is running, at most 6 at once. Processes I did not start are never touched.
- Per-sample cost: if the pilot's median processor time per sample exceeds 90 s (massless), or Stage B's first 7
  exceed 120 s (massive), STOP that stage and report.
- Stage A: at most 40 + 46 x 12 + 6 x 12 = 664 samples (the planned count is about 250), at most 1.5 h wall. If the
  fit fails at total degree 32: one extension to 60 fit points, then STOP with motion mapped only (BORDERLINE).
- Stage B: at most 5 x 64 + 2 = 322 samples, plus one extension to 96 eps per z if a fit fails, at most 1.5 h wall.
  If it still fails: STOP, NP2 not scored.
- Stage C: C1 at most 42 x 12 samples. C2 runs only if U <= 250 and the projected wall time at the measured rate is
  at most 3 h. Otherwise STOP, report U and the projection, and NP3b is not scored.
- Critical-point counts: at most 6 counts, 1 h wall.
- Whole phase: at most 8 h wall. Exceeding it is a STOP, reported with what was reached.
- Output files: new names only (`debug/np_*.py`, `debug/data/np_*`). No saved file is overwritten.

## Decision gate

- **GO:** z9(g2) on F1 reconstructed exactly under the rule, passing every held-out point.
- **BORDERLINE:** z9's motion on F1 mapped at sampled points, with no accepted formula.
- **STOP:** a per-sample or degree cap above is exceeded before the motion is mapped.

NP2 and NP3b are scored separately and do not change the gate.

## Risks

- F1's own coincidences may move z9 to a degenerate place or remove it. NP1a tests this. If it fails, Stage A
  becomes a map of whatever replaces it, and Stage C still runs at random rates.
- The carriers may have other factors at generic rates. The 12-sample fallback handles this.
- At z = z9, the reducer with eps nonzero might meet a degenerate seed system. That would show as a master count
  other than 30, which is an outcome named in NP2.
- The order of vanishing at eps = 0 compares rebuilt functions of eps. A numerator and a denominator that vanish
  together cancel in the rebuilt form. The positive controls test whether a known branch survives this.
- Rational reconstruction from one 61-bit prime may not reach the coefficient heights. The second prime then
  disagrees, and a CRT lift from two primes doubles the fit samples, within the caps.
