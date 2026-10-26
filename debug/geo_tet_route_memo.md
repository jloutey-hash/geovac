# Geometric route to Y(t), part 2: the whole tetrahedron and its top sector

Date: 2026-10-25. Scripts: `debug/geo_tet_symbolic.py`, `debug/geo_tet_fp_top.py`, `debug/geo_tet_katz.py`, `debug/geo_tet_euler2.py`. Data and logs: `debug/data/geo_tet_*`.

**Blindness.** I read `debug/geo_y_route_memo.md` in full and Paper 59 at commit ddd88961 (`git show`): sec:momentum, sec:genus, sec:reduction and sec:density_ring, plus the docstring of `debug/ring_ext_QG_common.py`. I read nothing on the banned list. I computed no differential operator, for Y or for any piece of it, and no numerical monodromy of any operator. Every monodromy below comes from geometry: Legendre data, vanishing cycles and Katz's algorithm on matrices. I searched no literature. The standard facts I used, all from memory, are: the Baikov representation; the Cayley–Menger determinant; the Van Oosterom–Strackee solid-angle formula; Gauss–Bonnet for spherical triangles; Katz's middle convolution and the Dettweiler–Reiter algorithm; and the Euler-characteristic count of master integrals.

## 1. The variety, including the top sector

**Baikov form.** For three loops in D = 3 with no external momenta, the Baikov representation integrates over the six squared momenta d_e of the lines. Its Gram factor is raised to the power (D − L − 1)/2 = −1/2:

  Y ∝ ∫_{CM ≥ 0} ∏ dd_e · CM(d)^{−1/2} · ∏_e (d_e + m_e²)^{−ν_e}.

The six d_e are the squared edge lengths of the **dual tetrahedron**, whose four vertices are the faces of the position-space tetrahedron O,1,2,3, and CM is that tetrahedron's Cayley–Menger (Gram) determinant. The family over x is the double cover

  **W_x : y² = CM(d) · (d_a + 1)(d_b + 4)(d_t + x)**  over A⁶,

with integer-order poles on d_g1 = −9/4, d_g2 = −9 and d_C = 0. The masses enter only through these hyperplanes. Opposite edges in the dual tetrahedron are (C, t), (a, g2) and (b, g1), the same pairing as in position space.

**The same integral with the loop solved.** Integrating the ring of geminal, geminal and Coulomb lines first gives

  Y = 8π² ∫_{λ(p²,q²,w²) ≤ 0} pqw · ρ̃_a(p) ρ̃_b(q) ρ̃_t(w) · ∂²_{g1²g2²} Ω(p,q,w) dp dq dw.

The kernel Ω is π² arctan(√det M / D)/√det M. This is the Van Oosterom–Strackee formula for the solid angle of a spherical triangle whose vertex vectors e_i have lengths m_i = (g1, g2, 0) and inner products M_ij = (m_i² + m_j² + s_ij)/2. The Gram boundary λ = 0 is a boundary of the chain, not a branch locus, because the outer two-loop Baikov exponent is 0 in D = 3.

**The spherical reading.** A spherical triangle is two-dimensional, so Gauss–Bonnet makes its area an angle excess, a sum of arccosines of Gram ratios. The top-sector kernel therefore has weight one: it is a logarithm on the double cover defined by √det M. It cannot carry an elliptic or K3 period by itself. Such periods arise only where the Fock factors are integrated against it.

The massless Coulomb line makes e₃ a null vector, so the triangle has a vertex on the null cone. This is what degenerates the geometry (S2 and S4 below).

## 2. The top sector is Artin–Tate

**S1 (exact).** The coefficient of d_e² in det G is −d_opp(e)/4. The apex line t is opposite the Coulomb line, so on the Coulomb shell d_C = 0 the determinant is **linear in d_t**.

**S2 (exact).** On the top-sector maximal cut (C, g1, g2 on shell), with P = d_a + g1² and Q = d_b + g2²,

  4 det G = g2²P² + g1²Q² − (g1² + g2² + d_t)PQ.

This is a binary quadratic form in (P, Q). At the apex Fock zero T = d_t + x = 0 it splits into two lines Q = μ±P through the point P = Q = 0. The discriminant of μ is (x − (g1+g2)²)(x − (g1−g2)²), which vanishes at x = 81/4 and 9/4.

With a massive line opposite the apex, the cut is quadratic in d_t, with d_t² coefficient m_C²/4. Its fibre over (A, B) is then an elliptic curve, not a conic.

**S6 (exact).** On the cut the apex-line form T^{−3/2}(c − βT)^{−1/2} dT is exact: it equals d[−2√(c − βT)/(c√T)]. The other basis form has closed period π/√β, which does not depend on c and hence not on x. So x enters the top-sector cut only through the two lines c = 0. On each line the double cover y² ∝ (P − α)(μP − β′) is a conic.

**Computation A: point counts on the cut** (`geo_tet_fp_top.py`). I computed S = Σ_{A,B,T ∈ F_p} χ(A B T · 4 det G_cut), at every good prime from 11 to 157 and x ∈ {7/3, 5, 13/2, 31/5, 2/7, −3/4}. In **187 of 187 rows** it equals the closed form derived from S2:

  S = −χ(−1) + p · Σ_μ [−χ(−1) − χ(−μ(g1² − a²)(g2² − b²))],

where μ runs over the F_p-roots of g1²μ² − (g1² + g2² − x)μ + g2². This involves only finite-order characters times powers of p. **No Frobenius trace of any curve or surface appears.**

Two further checks:
- **Number of top masters.** At primes where the μ-lines split, #{ABT·det G ≠ 0} = p³ − 4p² + 8p − 10 exactly, at x = 7/3 and at x = 5. So χ = −5: the top sector has five master integrals, all with Artin–Tate homogeneous solutions.
- **Massive control.** The massive-Coulomb control did not show a curve trace at the level of p either: its p-level digits stay within ±7. So the point counts do not by themselves isolate the masslessness as the cause; that mechanism rests on S1 and S2.

**S4 (exact): the kernel's singular loci.** For the D = 3 ring (g1, g2, 0):
- D = (g2P′ + g1Q′)/2, with P′ = p² + g1² and Q′ = q² + g2².
- D² + det M = P′Q′(w² + (g1+g2)²)/4. The log divisor is three planes: the two geminal poles and the geminal-pair threshold.
- det M equals minus the top-cut quadric of S2. The kernel's own branch locus is the same pair of planes as the Baikov top cut, found by an independent route.
- On the Gram boundary (collinear momenta, w = p + q), 4 det M = (pq(p+q) + g2²p + g1²q)² is a perfect square. The two factors 2D ± i√(4 det M) split into three linear forms each: ±i(p ∓ ig1)(q ∓ ig2)(p + q ± i(g1+g2)). So on the boundary the solid angle is a sum of three arctangents of linear forms, which is the degenerate spherical triangle.

## 3. The top sector's own singular point hides at x = 0

**S3 (exact).** The leading Landau singularity is the flat dual tetrahedron: all six lines on shell and det G = 0. It is linear in x, with root

  x_LS = −(a g2 − b g1)(a g2 + b g1)(a² − b² − g1² + g2²) / ((a² − g1²)(b² − g2²)).

At (a, b, g1, g2) = (1, 2, 3/2, 3) we have a g2 = b g1 = 3. Then 4 det G = 25x/4, so x_LS = 0. At a control point with g2 = 5/2 the same root is x = 44/45, a ninth genuine point that this parameter choice hides.

On the cut, x_LS is where the branch point P = α of one μ-line meets the other branch point β′/μ: μ = (g2² − b²)/(g1² − a²) = 4. This explains why the given singular list contains only cut values. It also means that x = 0 carries top-sector (Tate) structure on top of the transcendental pieces.

## 4. Census of strata (S5, two-loop form)

Every place where a period can live is an intersection of these loci: the Gram boundary λ = 0, the kernel's log planes, its branch planes, and the Fock zeros.

| stratum | what remains of y² = (p²+1)(q²+4)((p+q)²+x) | piece |
|---|---|---|
| λ = 0 (whole boundary) | all six Fock lines | K3 (six-line double plane) |
| λ = 0, p = ig1 | Fock_b, Fock_t at w = q + ig1: genus 1 | E_B (x = 1/4, 49/4) |
| λ = 0, q = ig2 | genus 1 | E_A (x = 4, 16) |
| λ = 0, p + q = i(g1+g2) | genus 1, but independent of x | constant elliptic period × (x − 81/4)^{−1/2} |
| triple point p = ig1, q = ig2 | concurrent with p + q = i(g1+g2) | character at 81/4 |
| point p = −ig1, q = ig2 | 25(4x − 9)/16 | character at 9/4 |
| interior, P′ = 0 or Q′ = 0 | det M is a perfect square there; branched on two lines | rational (Tate) |
| interior, w² = −(g1+g2)² | det M is a perfect square (threshold) | Tate × (x − 81/4)^{−1/2} |
| interior, top (ABT · det M) | conic fibres over (A, B) | Artin–Tate (Computation A) |

Each geminal line removes exactly one Fock factor, because its direction is a coordinate direction (p, q or p + q). So no stratum produces a curve of genus 2 or more. The list of transcendental pieces is closed: **K3, E_A and E_B**.

## 5. Assembly by middle convolution (`geo_tet_katz.py`)

Y/t³ = const · ∫₀^∞ √z (z + x)^{−3/2} Φ(z) dz, where Φ is the kite. This is an Euler transform whose kernel has monodromy −1. So each irreducible local system W of √z·Φ maps to the Katz middle convolution MC₋₁(W), which is irreducible. Rank-one pieces can also be created at the boundary.

I implemented the Dettweiler–Reiter algorithm in exact arithmetic and validated it two ways:
- MC₋₁ of x^{−1/2}(x − 1)^{−1/2} is Legendre's local system: J2 at 0 and 1, −J2 at ∞.
- MC₋₁ ∘ MC₋₁ returns the original character.

**The K3 piece.** The D = 3 density bubble is P′_{−1/2}((a² + b² − v)/(2ab)) (Paper 59's function 𝓛(a,b;c) at c = L). Its local system is Legendre's: unipotent at v = 1 and 9, −J2 at ∞, trivial at v = 0. The factor √z adds −1 at v = 0. MC₋₁ then gives:

| x | Jordan form | rank(T − I) |
|---|---|---|
| 0 | unipotent, blocks 2 + 2 | 2 |
| 1 | eigenvalues (1, 1, 1, −1) | 1 |
| 9 | eigenvalues (1, 1, 1, −1) | 1 |
| ∞ | −(J3 ⊕ J1) | 4 |

- **Rank 4**, with exactly one invariant bilinear form, and it is **symmetric**.
- This is the previous memo's prediction exactly, reached by a route independent of Bailey's formula and of its point counts.

**Slice control (a = b = 1).** The same algorithm gives rank 3: a single unipotent J3 at x = 0, a reflection at x = 4 (t = 2), and −J3 at ∞, with a symmetric form. This matches Paper 59's order-3 K3 piece, singular only at t = 0, ±2, ∞ and keeping its Jordan block of size 3 at t = 0. The collision of the pseudo-threshold (b − a)² with x = 0 is what turns 2 + 2 into a single block of 3.

**E_A.** MC₋₁ of the rank-one cut character (−1 at v = 0, 4 and 16) has rank 2, with J2 at x = 0, 4 and 16. At ∞ it gives −I, whereas the previous memo's direct curve gives −J2. That detail stays open; it does not affect the largest piece.

**The geminal-pair curve.** The kite's top-sector cut is a curve with I₁ fibres at 9/4 and 81/4 and I₄* at ∞. Its image depends on the twist at v = 0:
- with no twist, the image is **rank 1**, a character with −1 at 9/4 and 81/4;
- with a √z twist, it would be a rank-4 "geminal K3".

The vacuum top cut is Artin–Tate (Computation A), which selects the rank-1 case.

## 6. Order count

Sector master counts, from Euler characteristics:
- **Sunrise: 4.** χ = −4 analytically, and point counts at split primes give p³ − 4p² + 6p − 7, exact at 35 of 36 primes. This equals the rank of the K3 piece.
- **Sector {g1, C}: 3.** In D = 3 it collapses to a Yukawa line, with χ = 3 (p² − 3p + 5, exact at 34 of 36 primes): E_B plus one more.
- **Sector {g2, C}: 3**, by the same argument with g2 in place of g1.
- **Sector {g1, g2}: 3 masters**, all rational multiples of (x − 81/4)^{−1/2}.
- **Top: 5.**

The sectors {g1} and {g2} are the sunrise multiplied by a tadpole, so they are the same functions.

Hence the space of functions containing Y and its x-derivatives has dimension at most 4 + 3 + 3 + 1 + 5 = 16. On this count the integration-by-parts order 19 is not minimal. Integration by parts treats the factorized copies and the D = 3 bubble collapses as independent.

The transcendental part is fixed at 8 (K3 4, E_A 2, E_B 2), and at least one character must carry 9/4 and 81/4. The number of order-one pieces lies between 1 and 8, and the geometry does not fix it. Starting from the slice's four order-one pieces, which persist off the slice, and allowing the pseudo-threshold points to separate gives about 13.

## 7. Caps

- **Twisted cohomology dimensions.** The Euler-characteristic counts assume the exponents are non-resonant. D = 3 is resonant for the sub-bubbles, which is why the two-loop collapsed forms were used for the middle sectors.
- **The bound of 16.** It assumes Y's x-derivatives stay at the integer indices (2, 2, 1), which is true. It also assumes the collapsed sector counts are complete, which is not proved.
- **The massive control.** It was inconclusive.
- **E_A at ∞.** The monodromy at infinity of E_A and E_B is open between −I and −J2.
- **The Katz map.** It assumes the density content of the kite enters as rational multiples of the bubble, which is checked for the sectors with fewer than three integer lines, and covered by the variation-of-parameters argument for the top sector.

## PREDICTION

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
