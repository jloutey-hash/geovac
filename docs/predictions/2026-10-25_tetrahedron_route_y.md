# Pre-registered prediction: the whole tetrahedron, a second blind geometric route to Y (TG1-TG5)

Registered 2026-10-25, before the comparison. PI-directed ("Then launch the tetrahedron geometry agent"). The first
geometric route (`docs/predictions/2026-10-25_geometric_route_y.md`) named its own blind spot, the top sector; this
agent was asked to price the whole tetrahedron including it.

## Protocol

- One fresh agent (Opus), blind. It read Paper 59 as of commit ddd88961, before the pieces of Y were written in, and
  the first geometric memo `debug/geo_y_route_memo.md` with its scripts. It was barred from both PF and GR records,
  the atlas, the claim matrix, CHANGELOG above v5.21.10, the working copies of Paper 59 and its synthesis, every
  `debug/pf_*` and `debug/ib_tetra_*` file and their data, and git history after ddd88961; and from computing any
  operator or numerical monodromy for Y. Memo `debug/geo_tet_route_memo.md`; scripts `debug/geo_tet_*.py`; data
  `debug/data/geo_tet_*`.

## What the scorer already knew when this was recorded (disclosed)

The factorisation of disposition addendum 3 of the PF record: seven order-one factors and an order-12 factor F_8,
with no order-one factor at either end of F_8, and no R_V or Legendre pullback at either end within a map box whose
denominator allowed poles only at the eight cuts. On reading this prediction the scorer found that the box omits
poles at the sources' own apparent singular points (for the Legendre pullbacks x = B^2 - gamma^2), so that negative
may be an artifact of the search. The re-run with the fixed box is the comparison.

## Scoring rules

TG1-TG5 are the agent's T1-T5, scored HIT, MISS, PARTIAL or UNDECIDED.
- TG1: HIT if the transcendental composition factors of Y's minimal operator are exactly R_V (order 4) and the two
  crossed Legendre pullbacks (order 2 each), each once, and every other factor has order one with an algebraic
  solution. MISS if an irreducible factor of order 5 or more, or of order 2 to 4 of another type, is shown.
- TG2: HIT if the largest irreducible factor has order 4.
- TG3: the local monodromy of the order-4 factor from its exponents and logarithms; PARTIAL if some points match.
- TG4: HIT if the order of Y's minimal operator is 13; PARTIAL if 11 to 16; MISS otherwise.
- TG5: HIT if the order-4 factor is equivalent to R_V by an exact map.
If the order-4 factor is not found, TG3 and TG5 are UNDECIDED.

### The prediction, verbatim

(Variable x = t², function Y/t³. For Y itself, flip the sign of the eigenvalues at 0 and ∞.)

**T1. How many irreducible transcendental pieces, and their orders.**
- **Prediction:** Three:
  - K3, of order 4;
  - E_A, of order 2 (x = 4, 16);
  - E_B, of order 2 (x = 1/4, 49/4).

  Each appears once, as a separate piece; none sits inside a larger irreducible piece. **The top sector contributes no transcendental piece.** It contributes only order-one pieces with algebraic solutions: a character branched at 9/4 and 81/4, plus Tate structure at x = 0, where its leading singularity sits.
- **Confidence:** moderate-high for the top sector; moderate for "each appears once".
- **Reason:** The apex line is opposite the massless Coulomb line, so the top cut is linear in d_t and binary-quadratic in (P, Q). Its fibres are conics and its point counts are Artin–Tate (187 of 187 rows). Every other stratum loses one Fock factor per geminal line, so the census closes on the six-line double plane and two genus-one curves.

**T2. The order of the largest irreducible piece.**
- **Prediction:** 4.
- **Confidence:** moderate-high.
- **Reason:** It is the middle convolution of the rank-2 Legendre system of the density bubble, twisted by √z. Katz's rank formula and the exact algorithm both give 4. Nothing in the census is larger, and the rank-4 "geminal K3" is excluded by the Artin–Tate top cut.

**T3. Local monodromy of the largest piece.**
- **Prediction:**

  | x | eigenvalues | Jordan blocks | rank(T − I) | logarithms |
  |---|---|---|---|---|
  | 0 | 1, 1, 1, 1 | 2 + 2 | 2 | two, first power only |
  | 1 | 1, 1, 1, −1 (a reflection) | — | 1 | none |
  | 9 | 1, 1, 1, −1 (a reflection) | — | 1 | none |
  | ∞ | −1, −1, −1, −1 | 3 + 1 | 4 | one block of size 3 (log²) |

  The piece is trivial at x = 1/4, 9/4, 4, 49/4, 16 and 81/4.
- **Confidence:**
  - at x = 1 and 9: high;
  - at x = 0 and ∞: moderate-high (two independent routes now agree);
  - triviality at the six points: moderate-high.
- **Reason:** The input is Legendre monodromy at v = 1 and 9, −1 at 0 from √z, and unipotent at ∞, pushed through MC₋₁. The slice version of the same computation reproduces Paper 59's single J3 block at t = 0. The geminal scales never enter the density bubble's local system.

**T4. The order of Y's minimal operator.**
- **Prediction:** 13. Plausible range 11–16, and predicted strictly below 19.
- **Confidence:** low for 13; moderate-low that 19 is not minimal.
- **Reason:** The transcendental part is 8. The sector Euler characteristics (sunrise 4, {g1,C} 3, {g2,C} 3, {g1,g2} one character class, top 5) bound the function space at 16. The count of order-one pieces is not decided by the geometry.

**T5. Geometric type of the largest piece.**
- **Prediction:** The transcendental H² of the K3 double plane y² = (p² + 1)(q² + 4)((p+q)² + x), branched along six lines. It has rank 4 and Hodge type (1, 2, 1), weight 2, and is ≅ H¹(E₊) ⊗ H¹(E₋) for Bailey's conjugate curves. The monodromy preserves a **symmetric** form: orthogonal, O(2,2) with reflections, not symplectic.
- **Confidence:** high.
- **Reason:** The algorithm's invariant form is unique and symmetric, which is the middle convolution of a symplectic system by a quadratic character. It agrees with the previous memo's F_p identity S_p − 2 ≡ ±a_p(E₊)a_p(E₋) (mod p) in 117 of 117 rows.

## Disposition (appended 2026-10-25; prediction and scoring rules above unchanged)

**TG1: MISSED. TG2, TG3, TG4, TG5: UNDECIDED.** Data: addenda 3 and 4 of
`docs/predictions/2026-10-23_price_formula_slots.md`. With the denominator gap closed and at twelve primes, none of
R_V and the four Legendre pullbacks is equivalent to a right or left factor of F_8, and F_8 has no order-one end.
Y's minimal operator has a factor equivalent to an irreducible right factor of F_8, which is none of the types TG1
allows. Cap: the box.

- TG2 (largest piece of order 4): not decided, since the irreducibility of F_8 is not proved.
- TG3, TG5: the order-4 factor is not found, so UNDECIDED by the rule.
- TG4 (order 13): the order of the minimal operator is not decided.

Reading, not scored. The prediction agrees with the first geometric route on every piece, by two new computations
(an Artin-Tate count of the top cut in 187 of 187 rows, and a Katz middle convolution of the Legendre bubble). The
data agree with both routes on the local types at most points and disagree on the pieces. One place the data and the
prediction can be compared directly: at x = 1 and x = 9 F_8 has two eigenvalues -1, and the predicted K3 piece has
one. A rank-4 piece induced from E+ over the conic w^2 = (x - 1)(x - 9), rather than tensor-induced, would have two.
That is a hypothesis for the next test, not a result.

Addendum (2026-10-25): the induced-system hypothesis above was tested (PF record, addendum 5). The induced system
and its three twists are at neither end of F_8, at numerator degree up to 800. The scores stand.

## Disposition addendum 2 (2026-10-26; prediction, scoring rules and the disposition above unchanged)

**Rescored. TG1: HIT. TG2: HIT. TG3: PARTIAL. TG4: UNDECIDED. TG5: HIT.** Data: PF record, addendum 6. The disposition
above rested on a left-end map search that was void. Y's operator of order 15 has composition factors R_V, the two
crossed Legendre pullbacks and seven algebraic order-one pieces, each identification an exact map, and Y needs the
three transcendental ones.
- TG1: the transcendental composition factors of Y's minimal operator are exactly R_V and the two crossed pullbacks,
  each once; every other factor has order one with an algebraic solution. The top sector contributes two order-one
  pieces branched at 9/4 and 81/4, as predicted.
- TG2: the largest irreducible factor has order 4. Cap: R_V's irreducibility rests on a bounded probe.
- TG3: every read feature of the table matches (reflections at 1 and 9; trivial at the six points; two blocks at 0
  and at infinity, with the eigenvalues predicted). The block sizes 2 + 2 and 3 + 1 are not read.
- TG4: the order is between 8 and 15, not decided.
- TG5: the order-4 factor is equivalent to R_V by an exact map (56 primes, checked by exact division).
The induced-system hypothesis of the disposition is moot: the two eigenvalues -1 at 1 and 9 in F_8 are its two R_V
copies.

## Disposition addendum 3 (2026-10-26; prediction, scoring rules and earlier dispositions unchanged)

**TG1: HIT. TG2: HIT. TG3: HIT. TG4: HIT. TG5: HIT.** Data: PF record, addendum 7. Convention fixed before computing:
minimal over Q(x).
- TG1 now holds on Y's minimal operator itself: R_V, z_1 and z_2 once each, and five other factors, each of order one
  with an algebraic solution.
- TG3: the Jordan blocks are 2 + 2 at x = 0 and 3 + 1 at infinity, as in the table. Every feature matches.
- TG4: the order of Y's minimal operator is 13.
- Reading, not scored: T5's symmetric form holds (orthogonal, with reflections at 1 and 9; no alternating form).
  The Hodge type and the signature are not tested.
