# Pre-registered prediction: a geometric route to the pieces of the off-slice triangle Y (GR1-GR5)

Registered 2026-10-25, before the comparison. PI-directed: the PI asked whether a geometric route, built from the
Fock projection rather than the Feynman-integral toolkit, predicts the structure of Y better than the Landau-cut
reading behind PF1-PF4 (record `docs/predictions/2026-10-23_price_formula_slots.md`).

## Protocol

- One fresh agent (Opus) derived the prediction blind. It was told the object, the order-19 bound and the eight
  published singular points, and pointed at Paper 59's density-topology rule. It was barred from the PF record, every
  `debug/pf_*` and `debug/ib_tetra_*` file and their data, any CHANGELOG entry above v5.21.10, and from computing,
  guessing or factoring any operator for Y or its kite, or any numerical monodromy. It reports that it read Paper 59
  at the committed HEAD and no banned file. Its memo, scripts and data: `debug/geo_y_route_memo.md`,
  `debug/geo_y_k3.py`, `debug/geo_y_pieces.py`, `debug/data/geo_y_k3.*`, `debug/data/geo_y_pieces.*`.
- The agent saw neither the PF predictions nor any result about the pieces. It nevertheless names the same three
  sources PF1 and PF3 name: its E_A is the Legendre family at z = (t^2 - 8)/(2t), PF1's z_2; its E_B is at
  z = (4t^2 + 7)/(16t), PF1's z_1; its K3 piece is the period of V3(1, 2, t), PF3's R_V. This is recorded as a fact
  about the two routes, not as evidence for either.

## What the scorer already knew when this was recorded (disclosed)

The PM had already computed, and the agent had not seen:
- An exact composition chain of L (order 19): seven order-one factors, then one order-12 factor F_8, on which the
  probes of `debug/pf_lib.py` found no proper right or left factor. The probes are bounded (Frobenius exponents,
  degree at most 140), so this is not a proof of irreducibility.
- The local exponents of F_8 and of the order-one factors, and that the order-7 right divisor does not annihilate Y.
- A numerical monodromy run whose irreducibility tests did not separate the planted reducible control, so it decides
  nothing.

The scoring rules below are therefore written by someone who knows part of the answer. They are the plain reading of
each answer and use only checks named here. The PI may override any score.

## Scoring rules

Each of GR1-GR5 is the agent's answer to Q1-Q5 below, scored HIT, MISS, PARTIAL or UNDECIDED.
- GR1 (Q1): HIT if the composition factors of Y's minimal annihilator are one of order 4 equivalent to R_V, two of
  order 2 equivalent to the Legendre pullbacks at z_1 and z_2, and three or four of order one with algebraic
  solutions branched as stated. MISS if there is an irreducible factor of order 5 or more, or the order-2 or order-4
  factors are absent. The test is an exact map search from each source into F_8 and its adjoint, with planted
  controls, after which F_8's remaining part is factored by the same probes.
- GR2 (Q2): the local monodromy of the order-4 factor at each point, read from its exponents and logarithms. PARTIAL
  if some points match.
- GR3 (Q3): the points where that factor's local monodromy is trivial.
- GR4 (Q4): the order of Y's minimal annihilator. HIT if 11; PARTIAL if within 10-13; MISS otherwise.
- GR5 (Q5): HIT if the order-4 factor is equivalent to R_V by an exact map and its monodromy has a symmetric
  invariant form. If GR1 misses because the order-4 factor is absent, GR2, GR3 and GR5 are UNDECIDED, not missed.

### The prediction, verbatim

(Variable x = t^2, function Y/t^3. For Y itself, all eigenvalues at x = 0 and ∞ flip sign.)

**Q1. Composition.**
- **Prediction:** Six to eight irreducible pieces:
  - one of order 4, the three-density K3 piece, with no algebraic solutions;
  - two of order 2, the elliptic pieces E_A (at x = 4, 16) and E_B (at x = 1/4, 49/4), with no algebraic solutions;
  - three or four of order 1, with algebraic solutions, branched at 81/4, at 9/4, and possibly at {1, 9} (each possibly also at 0 and ∞).

  The largest irreducible piece has order 4.
- **Confidence:** moderate (the order-1 count is low).
- **Reason:** Each cut polygon carries the period of its Fock sides: three sides give a K3, two give an elliptic curve at imaginary momentum set by the geminal pole, one gives genus 0. Unequal scales make the two K3 conjugates non-isogenous (F_p counts), so Sym^2 ⊕ ∧^2 fuses to an irreducible rank 4, while E_A and E_B separate.

**Q2. Local monodromy of the rank-4 piece.**
- **Prediction:**
  - x = 0: eigenvalues 1, 1, 1, 1; Jordan blocks 2 + 2; rank(T − I) = 2; two log blocks, with first powers of the log only.
  - x = 1 and x = 9: eigenvalues 1, 1, 1, −1; rank(T − I) = 1; no logs (a reflection).
  - x = 1/4, 9/4, 4, 49/4, 16, 81/4: trivial.
  - x = ∞: eigenvalues −1, −1, −1, −1; Jordan blocks 3 + 1; rank(T − I) = 4; one log block, of size 3 (log^2).
- **Confidence:** high at 1 and 9; moderate at 0 and ∞ (low for the signs there); moderate for the six trivial points.
- **Reason:** It is the tensor induction of Legendre(λ = (u+1)/(3−u)) along the conic w^2 = (x−1)(x−9). The conjugates are exchanged at its branch points, where concurrent branch lines make a D4 point. Over x = 0 one conjugate degenerates (type II); over x = ∞ both do (type III).

**Q3. Points the rank-4 piece does not see.**
- **Prediction:** x = 1/4, 9/4, 4, 49/4, 16, 81/4. It sees only 0, 1, 9 and ∞.
- **Confidence:** moderate.
- **Reason:** The six-line double plane contains a, b and t only. The geminal scales enter other pieces through residues, not this surface.

**Q4. Order of the minimal annihilator.**
- **Prediction:** 11 (plausible range 10–13).
- **Confidence:** low.
- **Reason:** 4 + 2 + 2 + about 3 order-1 pieces. This assumes multiplicity 1 for each piece and no new top-sector period; off the slice there is no leading (dual-tetrahedron) singularity.

**Q5. Geometric type of the largest piece.**
- **Prediction:** The transcendental H^2 of the K3 double plane y^2 = (p^2+1)(q^2+4)((p+q)^2+x), branched along six lines; generic Picard number 18. Its rank is 4 and its Hodge type is (1, 2, 1), weight 2. It is isomorphic to H^1(E_{λ(u)}) ⊗ H^1(E_{λ(−u)}) (a Shioda–Inose partner of E+ × E−). The monodromy preserves a symmetric form: orthogonal, in O(2,2) with reflections, not symplectic.
- **Confidence:** high for this K3 piece; moderate that it is the largest piece of Y.
- **Reason:** Bailey's moduli become Möbius in the conic parameter, and the F_p character sums match ±a_p(E+)a_p(E−) in 117 of 117 rows. The tensor of two symplectic H^1's is orthogonal.

## Disposition (appended 2026-10-25; prediction and scoring rules above unchanged)

**GR1: MISSED. GR2, GR3, GR4, GR5: UNDECIDED.** The data are those of disposition addendum 3 of
`docs/predictions/2026-10-23_price_formula_slots.md`. In short: L = F_8 o F_7 o ... o F_1 with seven order-one
factors and one factor F_8 of order 12; R_7 does not annihilate Y; F_8 has no order-one right or left factor; none of
R_V, the Legendre pullbacks at z_1, z_2 and the direct d_1, d_2 is equivalent to a right or left factor of F_8 or
of L within the stated box, with planted controls found at the full box.

- **GR1: MISSED.** The prediction says every composition factor of Y's minimal annihilator M is the order-4 K3
  piece, one of the two elliptic pieces, or an algebraic order-one piece. M has a composition factor equivalent to an
  irreducible right factor of F_8, and that factor is none of these. Cap: the box.
- **GR2, GR3, GR5: UNDECIDED.** The order-4 piece is not found at either end of F_8. Its local monodromy and form
  cannot be read.
- **GR4: UNDECIDED.** The order of M is not decided.

Reading, not scored: what the geometry got right. These parts can be checked against the data without the pieces:
- the algebraic pieces branched at 81/4, at 9/4 and at {1, 9}: all three groups occur (three, three and one
  order-one factors of L);
- the hard part blind to 9/4 and 81/4: F_8 has trivial local monodromy at both;
- eigenvalue -1 and no logarithm at 1 and 9: F_8 has exactly that, with two such exponents, not one;
- one logarithm at each of 1/4, 4, 49/4, 16: F_8 has exactly that.

What it got wrong: the transcendental part sits in a block of order 12 with no predicted piece at either end and no split found, not
three separate pieces of orders 4, 2 and 2. The agent named the risk itself. The top sector's kernel is not rational
on the Fock curves, so the density-topology rule does not govern it, and a top-sector piece could not be excluded.
Both routes, the Landau reading and the geometry, located the same thresholds and the same candidate pieces, and both
missed that block. Whether it is truly one irreducible piece is not proved.
