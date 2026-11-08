# Geometry Exploration Report: tetrahedron and simplex geometry for the tetrahedral integrals Y(t)

Explorer agent, 2026-10-25. Read-only, apart from this memo. Corpus state was checked against Paper 59 (sec. three_electron, sec. four_density, the correlated-triangle and IBP paragraphs). I did not open the restricted predictions records, the atlas, the claim matrix, recent CHANGELOG entries or the pf_/ib_tetra/geo_ debug files.

## Search summary

About 45 searches and fetches. Every citation below was checked on an arXiv, publisher or journal page unless marked UNVERIFIED. The corpus already cites Davydychev–Delbourgo, Fromm–Hill, Remiddi, BBBG and Fresán–Sabbah–Yu, and already identifies the six-line double cover of the three-density chain (Paper 59, sec. three_electron).

**Bridge, noted and not built on.** One-loop N-point integrals in N dimensions are volumes of spherical or hyperbolic (N−1)-simplices: A. I. Davydychev and R. Delbourgo, J. Math. Phys. 39, 4299–4334 (1998), arXiv:hep-th/9709216. A three-point loop in D=3 is the area of a spherical triangle, which is our solid-angle kernel.

**One derived fact used throughout.** Put the nucleus at a vertex and the three electrons at the others. Then
d³r₁ d³r₂ d³r₃ = (8π²/3) · ∏ₑ ℓₑ dℓₑ / V(ℓ),
over the cone where the six edge lengths form a real tetrahedron (Cayley–Menger 288V² > 0). Mirror images are counted. I derived this and checked the single-vertex step by quadrature, 0.7566 against 0.7548 ± 0.002. A raw 6-D Monte Carlo is unreliable here because 1/V gives it infinite variance. K₄ is self-dual, so the momentum-space form is the same measure on the *dual* tetrahedron, with Fourier-transformed edge functions. Y(t) is therefore an integral over the space of tetrahedra, in both representations.

---

## Ranked candidates

### 1. Double covers of projective space branched along hyperplane arrangements (Aomoto–Gelfand E(k,2k) at exponent −1/2)

**What.** The integral of a form against ∏ ℓⱼ^{−1/2}, over the complement of 2k hyperplanes in P^{k−1}, is a period of the double cover w² = ∏ℓⱼ. It is Calabi–Yau: for k=2 an elliptic curve (Legendre), for k=3 a K3 surface, for k=4 a "double octic" threefold. The k=3 system E(3,6) and its monodromy are classical.

**Sources (verified).**
- K. Matsumoto, T. Sasaki, M. Yoshida, "The monodromy of the period map of a 4-parameter family of K3 surfaces and the hypergeometric function of type (3,6)", Internat. J. Math. 3(1), 1–164 (1992).
- S. Cynk, D. van Straten, "Picard–Fuchs operators for octic arrangements I (the case of orphans)", Commun. Number Theory Phys. 13(1), 1–52 (2019), arXiv:1709.09752.
- Cynk–van Straten, "Periods of double octic Calabi–Yau manifolds", arXiv:1709.09751 (journal not checked).
- S. Hosono, B. H. Lian, H. Takagi, S.-T. Yau, "K3 surfaces from configurations of six lines in P² and mirror symmetry I", arXiv:1810.00606.
- I. M. Gel'fand, "General theory of hypergeometric functions", Dokl. Akad. Nauk SSSR 288(1), 14–18 (1986). The English Soviet Math. Dokl. pagination is UNVERIFIED.

**What it computes for us.**
- **The tree ladder.** In Paper 59's tree identity (eq:tree) each density contributes Δᵥ² = cᵥLᵥ² + 1 = (√cᵥ Lᵥ + i)(√cᵥ Lᵥ − i). That is a pair of parallel hyperplanes. So n linked densities give exactly a double cover of P^{n−1} branched along 2n hyperplanes. The CY dimension n−2 then gives the rungs: elliptic at n=2, K3 at n=3, CY3 at n=4.
- **A geometric proof of tree-independence.** For any tree the normals Lᵥ come from the incidence matrix: their only relation is ΣLᵥ = 0, and any n−1 of them are independent. So every tree is GL(n−1)-equivalent to every other, and the scales cᵥ are the only moduli. This is why the period depends only on n, an independent route to the sec:four_density statement.
- **"Elliptic squared" for the chain (prediction).** The six lines are special in two ways. The three intersection points of the pairs are collinear, all on the line at infinity, and the configuration is centrally symmetric (k → −k). The family therefore has 2 moduli inside MSY's 4. That predicts Picard number at least 18 and a transcendental lattice of rank at most 4, which is a Shioda–Inose-type link to a product of two elliptic curves. This would explain BBBG's K(k₋)K(k₊′) + K(k₊)K(k₋′), with the two curves of moduli k±. Not derived; testable.
- **The four-density operator.** It is a double octic branched along eight planes in four parallel pairs, so a candidate member of the Meyer arrangement families. Prediction: AESZ #58 is a double-octic Picard–Fuchs operator. If it is, Hodge numbers, monodromy and rigid (modular) fibres are already tabulated.
- **Geminals.** They act as imaginary translations L → L + iγ of a hyperplane pair. That keeps the pairs parallel but breaks the central symmetry. The prediction is a lower Picard number, a higher transcendental rank, and new irreducible pieces, which fits the order jump from 4 to 9 when γ ≠ 1. The Coulomb line has integer exponent, so it enters as a pole divisor. Y is then a period of the *open* cohomology H^m(X \ D_pole), a mixed Hodge structure. The solid angle and the algebraic pieces are its extension data.
- **Caveat.** For the full tetrahedron I have not established a fibre form with linear L's. If the density lines stay rank-3 quadrics, the branch locus is quadrics and the right tool is candidate 2.

**Fit: high.** The arrangement is purely combinatorial (exact Layer-1 data), π enters only through periods, and the method is algebraic-first. **Cost:** low for the AESZ #58 lookup; low-to-moderate for point counts on the chain; moderate-to-high for Y.

### 2. Aomoto–Machida: hypersphere arrangements, Cayley–Menger determinants, and the Schläfli function

**What.** Aomoto studies integrals ∫ ∏ fⱼ(x)^{λⱼ} ϖ in which every fⱼ = |x|² + 2(αⱼ,x) + αⱼ₀ is a hypersphere, with arbitrary complex exponents. This is the Gauss–Manin theory of generalized simplex volumes. He gives the dimension of the twisted cohomology and an explicit connection whose coefficients are Cayley–Menger minors. Its singularities lie where principal minors vanish. A(jk)=0 means spheres j and k touch; in general, the spheres in J share a tangent affine subspace. Schläfli's volume formula is the origin of the variational formula.

**Sources (verified).**
- K. Aomoto, Y. Machida, "Hypergeometric integrals associated with hypersphere arrangements and Cayley–Menger determinants", Hokkaido Math. J. 49(1), 1–85 (2020), arXiv:1709.09329.
- Aomoto–Machida, "Generalization of Schläfli formula to the volume of a spherically faced simplex", J. Math. Soc. Japan 72(1), 213–249 (2020), arXiv:1710.10735.
- Aomoto–Machida, Proc. Japan Acad. Ser. A 91(6), 77–81 (2015). Title and pages are from search snippets only.
- K. Aomoto, "Analytic structure of Schläfli function", Nagoya Math. J. 68, 1–16 (1977). It writes the Schläfli function through Chen-type hyperlogarithms (Thm 1) and relates it to Mellin–Sato hypergeometric functions (Thm 2).

**What it computes for us.**
- **The one-loop sub-blocks.** A propagator (|q−p|² + m²)^{−ν} is exactly an Aomoto hypersphere, with center p and radius im, and λ = −ν (−3/2 for densities, −2 for geminals, −1 for Coulomb). Their rank formula, dim H = Σ_{ν=1}^{n} C(m,ν) + C(m−1,n), gives 7 for n=3, m=3 and 15 for n=3, m=4. These are predictions to compare with the IBP master counts per sector.
- **The "flat polygon".** The singular locus is the vanishing of Cayley–Menger minors, that is tangency of the spheres, which is the threshold and pseudo-threshold condition. This puts the PI's "flat polygon from the masses" in classical form.
- **The solid-angle kernel.** The JMSJ Schläfli-type formula, the unit-exponent limit, should reproduce its mass derivative in closed form.
- **Multi-loop.** The three-loop tetrahedron is a fibration of such arrangements. The inner q₃ integral is a 3-sphere arrangement whose centers are 0, q₁, q₂.

**Fit: high.** Hyperspheres are stereographic images of hyperplane sections of S³, and the coefficients are exact Cayley–Menger polynomials. **Cost: moderate.** The paper runs to 85 pp. The n=3 case is direct; the multi-loop fibration is new work.

### 3. Shape-space integrals on S³ and squared 6j symbols (Freidel–Louapre; Ponzano–Regge; Roberts; Taylor–Woodward)

**What.** Freidel and Louapre write
{6j}² = ∫_{SU(2)⁴} ∏_{I<J} χ_{l_IJ}(g_J g_I⁻¹) [dg],
and reduce it to an integral over the space of *spherical tetrahedra in S³*, with measure (2/π⁴) ∏ sinθ dθ / √det[cos θ_IJ]. This is the exact spherical analogue of our ∏ℓ dℓ/V; its flat limit is the measure above (my check: √det G → 6V). Roberts proved the Ponzano–Regge asymptotics. Taylor–Woodward showed that U_q(sl₂) 6j symbols at roots of unity are governed by spherical tetrahedra.

**Sources.**
- L. Freidel, D. Louapre, Class. Quantum Grav. 20, 1267–1294 (2003), arXiv:hep-th/0209134. Formula and measure read from the paper.
- J. Roberts, "Classical 6j-symbols and the tetrahedron", Geom. Topol. 3, 21–66 (1999), arXiv:math-ph/9812013.
- Y. U. Taylor, C. T. Woodward, "6j symbols for U_q(sl₂) and non-Euclidean tetrahedra", arXiv:math/0305113. Published in Selecta Math.; volume and pages UNVERIFIED.
- Taylor–Woodward, "Spherical tetrahedra and invariants of 3-manifolds", arXiv:math/0406228.
- S. Garoufalidis, R. van der Veen, "Asymptotics of classical spin networks", Geom. Topol. 17, 1–37 (2013), arXiv:0902.3113: rational generating series, G-function asymptotics.
- G. Ponzano, T. Regge, "Semiclassical limit of Racah coefficients", in Spectroscopic and Group Theoretical Methods in Physics (North-Holland, 1968), pp. 1–58. Pages from a secondary listing.

**What it computes for us.**
- **Coulomb kernels.** Fock's kernel 1/|ω−ω′|² = 1/(2−2cos θ) = Σ_{n≥1} χₙ(θ), the Chebyshev-U generating function at r=1. For points off the sphere, 1/(1−2r cos θ + r²) gives weights r^{n−1}. A tetrahedron of Fock Coulomb kernels on S³ is therefore Σ_{n_e} ∏ r_e^{n_e−1} {6j}²(n_e), a generating series of *rational* numbers ({6j}² is rational) indexed by GeoVac's own label n = 2j+1. The differential equations would come from Wilf–Zeilberger recurrences on the 6j's (the Schulten–Gordon recursion, J. Math. Phys. 16, 1961–1970 (1975)). Singularities would come from the Ponzano–Regge volume.
- **Masses.** In this picture masses become radial depths r inside B⁴ ≅ H⁴, the Davydychev–Delbourgo hyperbolic picture.
- **Wall.** One kernel needs one focal length. Lines with different masses do not share an S³ (the multi-focal wall, §3). This computes an exact single-focal analogue of Y, not Y itself.

**Fit: high** (S³, SO(4)/SU(2) characters, discrete labels, rational coefficients). **Cost: moderate.**

### 4. Polygon spaces and the phase space of a tetrahedron (Kapovich–Millson; Bianchi–Haggard; Aquilanti et al.)

**What.** The moduli space of spatial polygons with fixed side lengths is singular exactly when a collinear polygon exists, on the walls Σ εᵢrᵢ = 0. For 4-gons it is a 2-sphere, the reduced phase space of the 6j symbol. On it the tetrahedron volume generates closed orbits, whose period Bianchi–Haggard give through the complete elliptic integral K.

**Sources (verified).**
- M. Kapovich, J. J. Millson, J. Differential Geom. 44(3), 479–513 (1996).
- E. Bianchi, H. M. Haggard, Phys. Rev. Lett. 107, 011301 (2011), arXiv:1102.5439.
- Bianchi–Haggard, "Bohr–Sommerfeld quantization of space", Phys. Rev. D 86, 124010 (2012), arXiv:1208.2228.
- V. Aquilanti, H. M. Haggard, A. Hedeman, N. Jeevanjee, R. G. Littlejohn, L. Yu, J. Phys. A 45, 065209 (2012), arXiv:1009.2811.

**What it computes for us.**
- **Singular points.** The walls ε·m = √x are exactly thresholds and pseudo-thresholds, the "flat polygon" condition. A wall exists only if the corresponding polygon is a cut of the dual diagram. That is a selection rule to test against Paper 59's two ordinary points, z = −25/4 and −25.
- **An elliptic family from pure tetrahedron geometry.** Level sets of V² on the 4-gon sphere form one. Test: compare its j-invariant and monodromy with the Legendre piece, whose modulus λ = (γ²−(t−1)²)/(4t) = sin²(C/2), C being the angle of the mass triangle opposite the geminal.

**Fit: moderate.** These are SU(2) coadjoint orbits, which quantize to discrete spins. **Cost: low to moderate.**

### 5. Constant-curvature simplex volumes, used as a sorting rule

**What.** Schläfli: dV = ½ Σ ℓₑ dθₑ in 3-D, and d(area) = Σ dαᵢ for spherical triangles. Closed volume formulas for spherical and hyperbolic tetrahedra exist in dihedral angles and in edge lengths, built from Lobachevsky and dilogarithm functions. Goncharov: hyperbolic volumes are Borel-regulator values, periods of mixed Tate motives. Aomoto dilogarithms (BGSV) are the motivic version for pairs of triangles.

**Sources (verified unless marked).**
- H. S. M. Coxeter, Quart. J. Math. os-6, 13–29 (1935).
- R. Kellerhals, Math. Ann. 285 (1989); pages 541–569 UNVERIFIED (DOI 10.1007/BF01452047).
- Y. Cho, H. Kim, Discrete Comput. Geom. 22, 347–366 (1999).
- J. Murakami, M. Yano, Comm. Anal. Geom. 13, 379–400 (2005).
- J. Murakami, A. Ushijima, J. Geom. 83, 153–163 (2005), arXiv:math/0402087 (edge lengths).
- A. Ushijima, arXiv:math/0309216, in Non-Euclidean Geometries (2006), pp. 249–265; pages from a secondary listing.
- D. A. Derevnin, A. D. Mednykh, Russian Math. Surveys 60(2), 346–348 (2005).
- J. Murakami, Proc. AMS 140, 3289–3295 (2012), arXiv:1011.2584 (spherical tetrahedra).
- N. Abrosimov, A. Mednykh, survey, arXiv:1302.4919.
- A. B. Goncharov, J. AMS 12, 569–618 (1999), arXiv:alg-geom/9601021.
- Beilinson–Goncharov–Schechtman–Varchenko, Grothendieck Festschrift I, Progr. Math. 86, 135–172 (1990).
- J. Milnor, "The Schläfli differential equality", Collected Papers vol. 1 (Publish or Perish, 1994); title wording varies between sources.

**What it computes for us.** A firm negative. Anything that is a 2-D or 3-D constant-curvature simplex volume is mixed Tate of weight at most 2 (arctan/log or dilog/Cl₂), with unipotent local monodromy. That covers the solid-angle kernel and the one-centre dilog values (Remiddi, Fromm–Hill). It *excludes* simplex volumes as the source of the Legendre and K3 blocks. Cheap test: check that the kernel's mass derivative obeys the 2-D Schläfli relation, with the angles taken from the Gram data.

**Fit: moderate. Cost: low.**

### 6. The Regge symmetry and its geometry

**What.** Regge found extra linear symmetries of the 6j symbols. Together with the tetrahedral symmetries they generate a group isomorphic to S₃×S₄, of order 144. They preserve volume and dihedral data and are scissors congruences. Doyle–Leibon extend the group to 23040 symmetries of generalized hyperbolic tetrahedra. Boalch ties the Regge symmetry to Okamoto's symmetry of Painlevé VI.

**Sources.**
- Roberts (above): each generic Euclidean tetrahedron has 12 scissors-congruent, non-congruent partners.
- Y. Mohanty, Algebr. Geom. Topol. 3, 1–31 (2003), arXiv:math/0301318.
- A. Akopyan, I. Izmestiev, Bull. LMS 51(5), 765–775 (2019), arXiv:1903.04929.
- P. Doyle, G. Leibon, arXiv:math/0309187.
- P. Boalch, Comm. Math. Phys. 276, 117–130 (2007), arXiv:math/0603398.
- T. Regge, Nuovo Cimento 11, 116 (1959): UNVERIFIED.

**What it computes for us.** A cheap covariance test: apply the 144 maps to the six edge masses and see whether Y's singular set or operator transforms with them. Prior is low, since Regge maps mix line types; look at γ=1, the equal-mass case.

**Fit: low to moderate. Cost: very low.**

### 7. Integral geometry of random tetrahedra

**What.** Moments of random-simplex volumes. Miles, Adv. Appl. Probab. 3(2), 353–382 (1971); Zinani, Monatsh. Math. 139, 341–348 (2003), where the mean volume in a cube is 3977/216000 − π²/2160, the π² coming from the boundary; Ribando, Discrete Comput. Geom. 36, 479–487 (2006), solid angles as hypergeometric series in inner products.

**What it computes for us.** Y is an expectation over the shape measure, but the known results are moments of V, not products of exponential or K₀ edge weights. Use is limited to the Gaussian-geminal arc (Miles) and to expanding the solid-angle kernel in Gram entries (Ribando).

**Fit: low. Cost: low.**

### 8. Positive geometries and the amplituhedron

**What.** A positive geometry carries a canonical form with logarithmic poles on its boundaries. For polytopes, that form is the volume of the dual polytope.

**Sources.**
- N. Arkani-Hamed, Y. Bai, T. Lam, "Positive geometries and canonical forms", JHEP 11 (2017) 039, arXiv:1703.04541.
- N. Arkani-Hamed, E. Y. Yuan, "One-loop integrals from spherical projections of planes and quadrics", arXiv:1712.09991 (no journal reference listed).

**Assessment.** The resemblance is concrete in one place only: the one-loop solid-angle kernel, where "simplex relative to a quadric, with Gram data" is literally a spherical-simplex volume (Davydychev–Delbourgo, Arkani-Hamed–Yuan). Beyond it, it is empty. Polytope canonical forms have dlog poles and give only mixed Tate polylogarithms, never K3 or elliptic blocks; and the amplituhedron is massless, planar and four-dimensional, while our lines are massive, half-integer and in D=3. The real neighbour is candidate 2 (sphere-bounded simplices with arbitrary exponents), whose half-integer exponents give candidate 1's double covers, a step positive geometry does not take.

---

## Negative results (do not re-explore)

- Simplex volumes and Schläfli functions alone cannot give the elliptic or K3 pieces (candidate 5).
- The Cayley–Menger variety has no K3 of its own. CM is quadratic in each squared edge, so the measure ∏ dℓₑ²/√CM is genus 0 edge by edge. The transcendence comes from the edge functions (K₀ lines, which correspond to index 3/2), not from the shape space.
- Fromm–Hill (all edges e^{−wr}/r, and mass derivatives of them) stays at dilogarithms. Y climbs because its density lines are K₀ type, not Yukawa type.

## Connections to the existing framework

Paper 59's six-line K3 and tree identity are candidate 1's k=3 and general cases; new here are the E(k,2k) framing, the tree-independence proof and the Picard-rank prediction. The multi-focal wall (§3) is what blocks candidate 3 from reaching Y.

---

## Recommendation

**Route A: arrangement double covers (candidate 1). Most promising: purely algebraic, cheap first step.**
First computation, two cheap parts:
(a) Look up Paper 59's four-density operator (AESZ #58) in the Cynk–van Straten / Meyer double-octic operator lists. If it is there, the Hodge numbers, monodromy group and rigid (modular) fibres are already tabulated, and the four-density rung is identified as a double octic branched along eight planes in four parallel pairs.
(b) For the three-density chain, count points on w² = (c_b p²+1)(c_c q²+1)(c_a(p+q)²+1) (homogenized) over F_p for p ≤ 200, at two or three rational scale triples. Extract the Frobenius traces on H², and test the predictions: Picard number at least 18, and a transcendental part matching the product of the two BBBG elliptic curves with moduli k±. This is exact integer arithmetic.
Only after both: find out whether Y has a fibre form ∏(linear)^{−1/2} × rational. If it does, repeat the count for Y's order-twelve block, and read Y as a period of the open double cover X \ D_Coulomb.

**Route B: Aomoto–Machida hypersphere arrangements (candidate 2).**
First computation: treat the kite's inner loop (D=3, three or four spheres, exponents −3/2, −2, −1) as an Aomoto–Machida arrangement.
- Compare their rank formula (7 or 15) with the master counts per sector from Paper 59's IBP code.
- Evaluate their Cayley–Menger singular locus and check that it gives Y's eight genuine singular points and leaves z = −25/4 and −25 ordinary.
- Check that the unit-exponent limit (JMSJ 2020) reproduces the solid-angle kernel's Schläfli derivative.
