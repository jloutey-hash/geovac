# Geometric route to the off-slice correlated triangle Y(t): a blind prediction

Date: 2026-10-25. Scripts: `debug/geo_y_k3.py` and `debug/geo_y_pieces.py`. Data: `debug/data/geo_y_k3.{json,log}` and `debug/data/geo_y_pieces.{json,log}`.

**Blindness.** I read Paper 59 from the committed HEAD (`git show HEAD:...`), not from the working-tree copy, which has uncommitted edits. Its last paragraph on Y says "the pieces of Y off the slice ... are not yet in hand", so nothing I read gives the answer away. I read three Paper 59 guards: `test_paper59_three_electron_k3.py` (Bailey's formula), the outline of `test_paper59_generic_slice.py`, and the header of `test_paper59_correlated_triangle_g32.py`. I read no file on the banned list, and I computed no operator for Y, for its kite, or for any subsector. I consulted no literature. The standard facts used are Kodaira fibre types, Picard–Lefschetz theory, Kulikov degeneration types, and the rule that a triple point of a branch curve gives a D4 point; they are stated from general knowledge.

## 1. Setting: lines, Fock factors and the dual polygon

In position space Y is a tetrahedron on the vertices O, 1, 2, 3, with six lines:

- **Density lines** (index 3/2), each a Fock dispersion factor (k^2+m^2)^(-1/2) times rational terms: O–1 with a = 1, O–2 with b = 2, and the apex line O–3 with mass t.
- **Geminal lines** (index 2), rational: 1–3 with g1 = 3/2, and 2–3 with g2 = 3.
- **Coulomb line** 1–2 (index 1, massless), rational.

Only the integer-index lines can be shrunk by integration by parts. A density line has a half-integer index, so it is never removed. The subsectors that depend on t are therefore the top sector, the three sectors with one rational line shrunk (S_12, S_13, S_23), and the sunrise sector with masses (a, b, t).

**The geometric dictionary I use.** Momentum conservation at a vertex closes a triangle whose sides are the masses of the three lines meeting there. A cut closes a polygon in the same way. This is the dual tetrahedron, with masses as edge lengths, the same reading as the triangle of sides t, 1, γ that Paper 59 gives for its elliptic piece. On this reading, t is singular exactly when a flat polygon can be made from the apex side and the other sides of a cut. The massless Coulomb side has length 0, so it lets a four-sided cut flatten as a triangle.

**Computation 1** (`geo_y_pieces.py`, all cuts and all sign choices). The cuts that contain the apex line give exactly the eight given points and no others:

| cut (side containing O) | lines crossing the cut | x = t^2 | Fock factors | threshold exponent β | rung carried |
|---|---|---|---|---|---|
| {O} | a, b, t | 1, 9 | 3 | −1/2 (eigenvalue −1) | K3 surface |
| {O,1} | b, t, g1, Coulomb | 1/4, 49/4 | 2 | 0 (log) | elliptic curve E_B |
| {O,2} | a, t, g2, Coulomb | 4, 16 | 2 | 0 (log) | elliptic curve E_A |
| {O,1,2} | t, g1, g2 | 9/4, 81/4 | 1 | −3/2 (eigenvalue −1) | genus 0 (algebraic) |

The "adjacent" triangles (t, a, g1) and (t, b, g2) are not cuts. They would put points at x = 25/4 and 25, which are exactly the two given ordinary points.

The exponent β is the near-threshold behaviour of a k-line cut in D = 3: β = ((k−1)(D−1) − 2)/2 + Σ(1 − ν_i). The cut with three density lines and the vertex star at 3 come out half-integer, so the local monodromy there has eigenvalue −1. The two four-line cuts come out integer, so the monodromy there is unipotent, with logarithms.

**No leading singularity off the slice.** For the vacuum tetrahedron, the leading Landau-type locus is the closed dual tetrahedron with edge lengths equal to the masses. Its Coulomb edge has length 0, so vertices 1 and 2 of the dual coincide. Closure then needs a = b and g1 = g2. Off the slice that is impossible, so the top sector has no singular point of its own. This matches the given list, which contains only cut values.

## 2. The rule: one Fock factor per side

Paper 59's density-topology rule says that two Fock factors build the curve y^2 = (c1 k^2+1)(c2 k^2+1), and that rational kernels are rational functions on it. Extended to Y, each cut polygon carries the period of the Fock factors on its sides:

- **Three Fock sides (a, b, t).** The double plane y^2 = (p^2+a^2)(q^2+b^2)((p+q)^2+t^2) is branched along six lines. It is a K3 surface, and its holomorphic period is the three-scale moment V3(1,2,t) = ∫K0(s)K0(2s)K0(ts) ds (Paper 59 eq:tree with n = 3).
- **Two Fock sides plus a geminal pole.** Consider the 1-D Fourier form of a sector with one geminal edge. The geminal factor 2g/(k^2+g^2) has a pole at k = ig. Its residue evaluates the two-Fock elliptic family at imaginary momentum, which is Paper 59's own reading. The curve is y^2 = (k^2+m^2)((k−ig)^2+t^2). This gives E_A (m = a, g = g2) and E_B (m = b, g = g1).
  - **Adjacency.** Why does the density a pair with the geminal on the opposite edge? The apex momentum is w = −(p+q). A pole of the geminal at vertex 1 fixes p, the same variable that Fock_a depends on, so the two cannot compose. The geminal at vertex 2 fixes q and shifts w. This is the origin of the ordinary points 25/4 and 25.
- **One Fock side (t).** Both geminal residues fix w, which leaves (t^2 − (g2 ± g1)^2)^(-1/2): genus 0.

## 3. The three-Fock piece: computations

**Computation 2: Bailey's moduli in x** (`geo_y_k3.py`). For masses (1, 2, t), both sorting regimes give the same discriminant:

  λ± = (x − 3 ± √((x−1)(x−9))) / (2x).

For general masses (A, B) the discriminant is (x − (A+B)^2)(x − (A−B)^2)/x^2. The two curves E± are therefore conjugate over the conic w^2 = (x−1)(x−9).

**Computation 3: parametrizing the conic.** Put x = (9 − u^2)/(1 − u^2). Then exactly:

  λ+ = (u+1)/(3−u),  λ− = λ+(−u).

So each conjugate is the Legendre family, pulled back by a Möbius map. The K3 piece is the tensor induction, from the u-line to the x-line, of the rank-2 Legendre local system. That is, its solution space is H^1(E_{λ(u)}) ⊗ H^1(E_{λ(−u)}) of rank 4, twisted by x^(−1/2) through Bailey's prefactor π/(2c).

Where the conjugates are singular:

- E+ is singular at u = −1, 1, 3 (λ = 0, 1, ∞); E− is singular at u = 1, −1, −3.
- The preimages u = ±3 lie over x = 0, and u = ±1 lie over x = ∞.
- The conic's branch points u = ∞ and u = 0 lie over x = 1 and x = 9. Both conjugates are smooth there (λ = −1 and λ = 1/3).

**Computation 4: Bailey against quadrature.** Agreement is ≤ 5×10^-31 at t = 0.5, 1.5, 2.5 and 4, which covers all three sorting regimes.

**Computation 5: the swap at x = 9 and x = 1.**

- Near t = 3: the swap-odd product K(k−)K(k+') − K(k+)K(k−') behaves like √ε. Divided by √ε it is −2.72071794, −2.72069924, −2.72069905 at ε = 10^-4, 10^-6, 10^-8, so it converges. The swap-even sum (which is V3) is analytic there.
- Near t = 1: the same ratio converges to −π, to 8 digits by ε = 10^-8.

So the local monodromy at x = 1 and x = 9 is the exchange E+ ↔ E−: eigenvalues (1, 1, 1, −1), with rank(T − I) = 1.

The same follows from geometry. At t = a ± b the three branch lines p = ±ia, q = ±ib, p+q = ±it become concurrent. The triple point makes a D4 rational double point on the double plane. The K3 does not degenerate, it gains one algebraic class, and the monodromy is the reflection in that (−2)-class.

**Computation 6: point counts over F_p.** Let S_p = Σ_{P,Q∈F_p} χ((P^2+1)(Q^2+4)((P+Q)^2+x)), the character sum of the double plane. I tested u = 1/2, 5/7, 7/3 and 2/5, which give x = 35/3, 52/3, −4/5 and 221/21, at every good prime from 11 to 139.

- In 117 of 117 rows, S_p − 2 ≡ ±a_p(E+)·a_p(E−) (mod p); in 100 of those rows the right side is not 0 mod p. The constant 2 comes from the lines at infinity. A random residue would match with probability about 2/p.
- So the unit-root part of H^2 is H^1(E+) ⊗ H^1(E−) of rank 4. This is an arithmetic check that does not use Bailey's formula.
- Also a_p(E+)^2 ≠ a_p(E−)^2 at many primes, so the conjugates are not isogenous, not even up to a quadratic twist. (The exception is x = −4/5, where λ+ = 5 and 1 − 1/λ− = 5 make them twists of each other. That is an isolated special fibre, not a physical one.)
- The same holds over Q(x) by reduction type: E+ has good reduction at u = −3 while E− has a pole of j there.

Therefore the K3 piece is irreducible of rank 4. Its monodromy lies in O(2,2), generated by the two SL2 factors and the swap reflections. The tensor of two alternating forms gives a symmetric form, so the monodromy is orthogonal, not symplectic. The Hodge type is (1, 2, 1) and the generic Picard number is 18.

**Local monodromy of the rank-4 piece**, from the Legendre data at the preimages:

- **x = 0.** At u = 3, E+ has λ = ∞ (monodromy −J2) and E− is smooth. That gives −(J2 ⊕ J2); the x^(−1/2) twist flips it to +(J2 ⊕ J2).
  - Geometric picture: the pair p+q = ±it collapses to a double line that passes through no other concurrency point. That is a type II degeneration, N ≠ 0 and N^2 = 0.
  - Independent check: Y/t^3 = ∫K0(tr) r^2 Φ(r) dr is A(x) log x + B(x) with A, B analytic at 0. On the slice the K3 factor is the head of the composition chain (top of the series ra < R2 < V4 < V5 < K6 < L9). If it is also the head off the slice, its eigenvalue at x = 0 must be +1.
- **x = ∞.** At u = 1, λ+ = 1 and λ− = 0, so each conjugate is J2 and the product is J2 ⊗ J2 = J3 ⊕ J1. The twist gives −(J3 ⊕ J1).
  - Geometric picture: the collapsing pair now merges with the line at infinity, which contains the two other concurrency points. That is type III (maximally unipotent).
  - V3 ~ t^(−1) log^2 t confirms both the log^2 and the sign.
- **The six other cut points.** The double plane does not depend on g1 or g2 at all, so its monodromy is trivial there.

**Slice control.** At a = b we get λ+ + λ− = 1. Then E− is E+ up to a twist, and the tensor splits as Sym^2 ⊕ ∧^2, of ranks 3 and 1. The rank-3 part is F3's operator, which is the slice's order-three piece. Its singular points move to t = 2 and t = 0. The ∧^2 part is a rank-1 character branched at t = 2, a candidate for one of the slice's "algebraic pieces branched at t = 2". Off the slice that character is absorbed into the irreducible rank-4 piece.

## 4. The elliptic and algebraic pieces

**Computation 7** (`geo_y_pieces.py`). The j-invariant of the shifted two-Fock quartic equals Legendre's, with z = (t^2 + m^2 − g^2)/(2mt), identically.

- E_A: j = 16(x^2−4x+64)^3 / (x(x−4)^2(x−16)^2), and λ = −(t−4)(t+2)/(4t). It has I_2 fibres at x = 4 and 16, and simple poles of j at x = 0 and ∞. The Euler number is 2+2+1+7 = 12, so it is a rational elliptic surface, like the slice's.
- E_B: j = (16x^2+824x+49)^3 / (4x(4x−1)^2(4x−49)^2), with I_2 fibres at x = 1/4 and 49/4.

The two have disjoint multiplicative supports, and neither shares a support with E± (whose j has poles only over x = 0 and ∞). So E_A, E_B and the K3 piece are pairwise non-isomorphic. On the slice E_A and E_B coincide, and the singular set reduces to Paper 59's {|1−γ|, 2, 1+γ, 2γ}.

Local monodromy in the Y/t^3 normalization:

- At each of the piece's two cut points: J2 (one logarithm).
- At x = 0: J2. The model C = ∫K0(tr)K0(mr)cosh(gr)dr is even in t plus a log.
- At x = ∞: −J2, since C ~ x^(−1/2) log x. This is the slice surface's I_1^*/I_1 pattern with 0 and ∞ exchanged by the x^(3/2) normalization.

The genus-0 pieces are characters √(x − 81/4) (from the residue at the geminal-pair pole in S_12) and √(x − 9/4) (from the top-sector vertex star). Each may carry a √x twist.

## 5. Order count, and the top sector

The slice at γ = 3/2 has order 9 = 3 (K3) + 2 (elliptic) + 4×1 (algebraic). Off the slice the geometry gives:

- K3: 3 → 4, because ∧^2 joins Sym^2.
- Elliptic: 2 → 2 + 2, because E_A and E_B separate.
- Algebraic: 4 → about 3. The ∧^2 character is absorbed; the two geminal-pair characters move to 81/4 and 9/4; one density-pair character either stays at {1, 9} or disappears.

That gives about 11, within the integration-by-parts bound of 19.

What this account assumes: the top sector adds no transcendental piece of its own. Two things support it.

- Off the slice there is no leading singularity, so a new top-sector period would have to live on cut points only.
- On the slice the top sector contributed nothing new.

But the top-sector kernel, the arctan loop with √(det G), is not a rational function on the Fock curves. The rational-kernel premise of the density-topology rule therefore fails there, and the geometry does not exclude a top-sector piece that was antisymmetric under the 1↔2 swap on the slice. That swap-odd possibility is the main risk to Q1 and Q4.

## 6. Caps

- The cut decomposition is derived and slice-checked, not proved off the slice; multiplicities above 1 are not excluded.
- The signs at x = 0 and ∞ assume the Y/t^3 normalization and that the K3 factor is the head.
- The algebraic count is a slice-based estimate.
- Point counts: four values of u, primes up to 139.

## PREDICTION

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
