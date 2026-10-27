# Literature scan: link topology -> period rung, and the kite's ninth singular point

Date: 2026-10-26. Scope: Paper 59 `sec:density_ring` + abstract; `docs/predictions/2026-10-26_parametric_kite_operator.md`
(Disposition, stage 1). Literature scan plus a small exact computation for Q2 (scripts in `debug/lit_scan_*.py`,
logs in `debug/data/lit_scan_*`). No papers or records edited. Nothing committed.

Verification key: **[V]** = abstract page (or full text where noted) fetched this session; **[V-prior]** = verified
in `debug/lit_scan/elliptic_eri_feynman_memo.md`; **[UNV]** = not verified against the primary source.

---

## Q1. Does a graph's topology fix the rung of its period?

### What the literature has

1. **Feynman integrals are periods; the motive belongs to the graph.** Bloch-Esnault-Kreimer, *On motives
   associated to graph polynomials*, CMP 267 (2006), arXiv:math/0510011 [V]. Bogner-Weinzierl, *Periods and
   Feynman integrals*, arXiv:0711.4863 [V]: every Laurent coefficient at rational Euclidean kinematics is a period.
   This is the general frame in which "set by the graph, not the basis or method" holds: the integral is a period
   of the graph hypersurface complement, whatever representation computes it. It says nothing about *which* rung.

2. **The nearest theorem in shape: Brown, *On the periods of some Feynman integrals*, arXiv:0910.0114** [V,
   full text read]. In the paper's words, "the set of periods of any graph ... is minor monotone", so "the set of
   graphs whose periods lie in a fixed ring (e.g. the ring of multiple zeta values) is determined by a finite set
   of forbidden minors". Theorem 2 + Thm 118: vertex width at most 3 implies the period is an MZV (mixed Tate),
   and a non-Tate example of vertex width 4 shows the bound cannot be improved. Scope: massless, primitive
   phi^4 graphs in D=4. This is a theorem that a topological invariant (a minor-closed class, vertex width) bounds
   the rung. It is the only genuinely theorem-shaped "rung by topology" statement found.

3. **First non-Tate rung is K3, set by graph combinatorics.** Brown-Schnetz, *A K3 in phi^4*, Duke 161 (2012),
   arXiv:1006.4064 [V]: graph-hypersurface point counts given by a modular form of a singular K3. Brown-Schnetz,
   *Modular forms in QFT*, arXiv:1304.5342 [V]: point counts mod q^3 (the c2 invariant; the term is not in the
   abstract, so the c2 wording is [UNV] at abstract level) match modular forms of weight <= 8, level <= 17, up to 10
   loops.

4. **Loop number sets the rung for bananas.** Bloch-Kerr-Vanhove, arXiv:1406.2664 [V]: three-loop banana = regulator
   periods of a family of K3s. Bonisch-Fischbach-Klemm-Nega-Safari, arXiv:2008.10574 [V]: the l-loop banana
   corresponds to a Calabi-Yau (l-1)-fold (GKZ, arbitrary masses). Bailey-Borwein-Broadhurst-Glasser,
   arXiv:0801.0891 [V-prior]: Bessel moments of K_0^n, the D=2 banana in position space.

5. **A loop-order bound on the rung.** Bourjaily-McLeod-von Hippel-Wilhelm, arXiv:1810.07689 [V]: rigidity of
   massless D=4 integrals bounded by 2(L-1) at L loops (argued, not a theorem). The arXiv title is "A (Bounded)
   Bestiary of Feynman Integral Calabi-Yau Geometries". Paper 59's bib uses the PRL title, which I did not check
   [UNV]; worth a look.

6. **Two-loop topology -> motive, case by case.** Doran-Harder-Pichon-Pharabod-Vanhove, *Motivic geometry of two-loop
   Feynman integrals*, arXiv:2302.14840 [V]. Planar two-loop graphs split into mixed Tate pieces plus families of
   hyperelliptic, elliptic or rational curves, *depending on D*. The double box is elliptic, the non-planar
   tardigrade is a K3 of Picard rank 11, and the ice-cream-cone family gives pairs of sunset CYs. This is the
   closest in spirit to a "topology table", and it shows the rung depends on D as well as the graph.

7. **Mechanisms that lower the rung.**
   - *Uniqueness / star-triangle*: a loop integrates out exactly when the indices are unique (sum = D at a vertex,
     D/2 around a triangle). Classical (D'Eramo-Peliti-Parisi 1971; Vasil'ev-Pis'mak-Honkonen 1981; Kazakov 1983/84,
     all [UNV] at the primary source). Collected with the general-D constant in Preti, arXiv:1811.04935 (cited in
     P59) and Preti, *The Game of Triangles*, arXiv:1905.07380 (search listing only, [UNV]). The literature uses it
     to *compute*. I found no source that states it as a rung-lowering rule.
   - *Hidden involution*: Marzucca-McLeod-Page-Pogel-Weinzierl, *Genus drop in hyperelliptic Feynman integrals*,
     arXiv:2307.11497 [V]. A hidden involution maps genus g to ceil(g/2) or floor(g/2) (crossed box, massive).
   - *Equal masses*: equal-mass sunrise is modular (Bloch-Vanhove arXiv:1309.5865, cited in P59). This is the same
     degeneration as P59's gamma = 1 drop.
   - *Dimension/index*: Doran et al.'s D-dependence; and see the remark on index 3/2 below.

### Mapping to Paper 59's cases

Read P59's integrands as **position-space graphs**. Densities are edges to the nucleus vertex, kernels are edges
between electron vertices, and at one centre all density edges meet one external vertex 0.
- **Chain of n densities**: a tree of kernels plus n density edges to 0. P59 already ties it to the banana ladder by
  an exact Hadamard product (sec. four_density). The rung-by-n statement is therefore BFKNS's loop-number rule
  transported. That much is prior art; the transport map and the operator identifications are P59's own.
- **Coulomb ring of three** = K4, the tetrahedron or wheel with three spokes: three Coulomb edges plus three density
  edges. In massless D=4 phi^4 this graph's period is 6 zeta(3), the first non-rational rung. P59's one-centre
  weight-two value fits the same "K4 adds weight" pattern, but that is only an analogy (different D, masses and
  indices). The fibre drop to elliptic happens because the Coulomb triangle is **unique** in D=3: three unit-index
  1/r edges, star-triangle. That is P59's reading, and it is the standard uniqueness identity.
- **Correlated triangle Y** = K4 again, with three massive edges (two geminals and a density pair) and one Coulomb
  edge. No index set is unique, so no loop collapses. The inner kite is a two-loop self-energy, and closing it with
  the third density gives the three-loop tetrahedron vacuum graph. For equal masses or zero masses in D=4 the
  massive tetrahedron is polylogarithmic at sixth roots of unity (Broadhurst, arXiv:hep-th/9803091 [V], "10 mass
  cases"). I found no generic-mass or D=3 classification of the tetrahedron that would place it on K3.
- **Remark (index 3/2 in D=3)**: (k^2 + a^2)^{-3/2} in D=3 Fourier-transforms to a K_0(a r) kernel, which is why the
  D=2 Bessel-moment/banana ladder (BBBG) appears in a D=3 problem. This is a standard dimension-shift fact, worth
  one sentence in the paper if it is not already there.

### Prior art vs new

- **Prior art**: rung set by graph + kinematics + D (BEK, Bogner-Weinzierl, Doran et al.); loop number fixes the
  banana rung (BFKNS, BKV, BBBG); a minor-closed topological class bounds the ring (Brown 0910.0114); uniqueness
  collapses loops (classical); the equal-mass drop (Bloch-Vanhove).
- **Not found anywhere**: (i) the reading of molecular integrals over linked densities as these graphs, with the
  chain/ring/kite ladder and its operators; (ii) the specific Coulomb-triangle operator L_G and its factorization,
  the four-density AESZ #58 entry, and the correlated triangle on K3; (iii) uniqueness stated as a rule that lowers
  the rung.
- **Important tension**: Brown's minor monotonicity says the *period set* (all dressings) of a graph contains that
  of its minors. The chain is a minor of the triangle (delete one edge), yet the Coulomb triangle sits *below* the
  chain at the fibre. This is no contradiction. Brown's theorem is about the ring generated over all dressings and
  regularizations, while uniqueness is an index-specific identity at one point of that family. But it shows the
  "rung by link topology" rule cannot be a topological statement alone. The theorem-shaped version is "rung bounded
  above by a minor-closed class of the graph; equalities at special indices (uniqueness) or masses (equal-mass) can
  lower it". P59's own sentence ("a property of uniqueness, not of rings") already points this way. Neither
  Brown's theorem nor anything found covers massive, D=3, non-integer-index graphs, so the rule is **not a special
  case of a published theorem**. It is closest to Brown's minor-monotone bound, at Observation tier, with an
  analogue in BFKNS for the tree (banana) column.

---

## Q2. The ninth singular point z0 = -111/1120

### Literature: singularities that standard Landau analysis misses

- **Second-type singularities** (pinches at infinite loop momentum, U = 0 in parametric form): Fairlie-Landshoff-
  Nuttall-Polkinghorne, J. Math. Phys. 3, 594 (1962), bibliographic data confirmed via Cambridge listing, text
  [UNV]. Klausen, arXiv:2109.07584 [V]: principal A-determinants give the singular locus and cover second-type
  singularities.
- **Principal Landau determinant (PLD)**: Fevola-Mizera-Telen, arXiv:2311.16219 and arXiv:2311.14669 [V]; Mizera-
  Telen, arXiv:2109.08036 [V]. PLD = the product over faces of Newt(U+F) of the face discriminants. In their words it
  "captures the singularities even in the presence of massless particles". Faces whose discriminant vanishes
  identically are discarded.
- **Whitney stratification**: Helmer-Papathanasiou-Tellander, arXiv:2402.14787 [V, html read]. They find
  components that PLD misses. For the slashed box with massless edges, PLD finds 8 of 9 components and misses
  p4^2 - s - t. For the **parachute**, an extra discriminant is "not captured by either direct principal
  A-determinant calculation nor the principal Landau determinant". On it "the Euler characteristic drops from the
  generic value 19 to the value 18". The diagnostic is the **drop of the Euler characteristic**: |chi| = the number
  of masters (Lee-Pomeransky arXiv:1308.6676 [V]; Bitoun-Bogner-Klausen-Panzer arXiv:1712.09215 [V]).
- **Massless lines**: Vergu, arXiv:2512.12763 [V]. Massless propagators make the on-shell spaces singular, so
  Pham's general results do not apply directly. They need resolution of singularities.

### Computation (exact, this session)

Kite with lines 1 (k,a), 2 (w-k,b), 3 (k-q,g1), 4 (w-k+q,g2), 5 (q,0), at P2 = (1, 3/2, 7/3, 8/3). G = U + F as in
`debug/pko_landau_leading.py`.

**(a) PLD-style face scan** (`debug/lit_scan_pld_kite.py`, log `data/lit_scan_pld_kite_p2.log`): Newt(G) has 30
monomials, 10 facets and 189 faces. For each face, does G_face = x_i dG_face/dx_i = 0 have a torus solution at z?
- Controls hit: -25/4, -1/4, -25, -1/9 each on three faces (the bubble faces), and -121/9 on the sunrise face
  {x1^2 x5, x1 x4 x5, x4^2 x5}.
- **z0: no non-degenerate face**.
- 49 faces are degenerate (a torus solution at every z), including the whole polytope. Singular points of
  {G=0} lie on a curve for every z (`data/lit_scan_sing_kite_p2.log`).
- So **naive PLD misses z0**, and z0 can only live in a degenerate stratum, the class Helmer et al. show PLD misses.

**(b) Euler-characteristic drop** (`debug/lit_scan_euler_kite.py`): count the critical points of
sum nu_i log x_i - s log G (generic rational nu, s; Groebner mod 32003, multiplicities counted).

| z | count |
|:--|:--:|
| random -37/101, -111/1121, -1/10, -2/7 | 22 |
| **-111/1120** | **21** |
| -121/9 (three-line cut with the Coulomb line) | 21 |
| -25/4, -1/4, -25, -1/9 (two-line cuts) | 17 |
| 0 | 8 |
| -73/9, -277/36 (system-only roots) | 22 |
| -100/9 (adjacent pairing, not a cut) | 22 |

- Repeat with prime 65521 and different nu, s: 22 at -5/13, 21 at z0 (`data/lit_scan_euler_kite_p2_alt.log`).
- **z0 is a genuine component of the Landau locus in the Whitney/Euler sense.**
- The controls behave as the operator says: genuine points drop, apparent and adjacent ones do not.
- The drop size (1) and the half-integer exponent match the three-line cuts that pass through the massless line.

**(c) Localization** (`debug/lit_scan_euler_sectors.py`, `debug/lit_scan_euler_face5.py`):
- Every sector with 1 to 3 lines contracted has the same count at z0 as at random z, so the drop is in the **top
  sector only**.
- The soft face x5 -> infinity (graph with line 5 deleted, q -> 0) has count 8 at both, so z0 is not soft.
- With a regulator mass m5^2 = 1/100 the generic count rises to 30 and z0 no longer drops: **the point exists
  because line 5 is massless** (it moves or disappears under a mass).
- The leading first-type Landau system at m5^2 = 1/100 is still the unit ideal, so z0 does not continue from a
  massive first-type leading singularity.

**(d) Curve probe** (`debug/lit_scan_curve_kite.py`): a lex Groebner basis with z symbolic, to find where the
singular curve of {G=0} meets x_i = 0 or is itself singular. It did not finish in about 25 min of sympy and was
stopped. Extracting the z-polynomial is the open step (recipe 2).

### Verdict and recipe

**Most likely class**: a top-sector, massless-line Landau component of mixed/second type. It lies in a degenerate
stratum of the U+F hypersurface: the singular locus of {G=0} is a curve for every z because line 5 is massless.
That makes it invisible to the first-type leading Landau equations and to the PLD face discriminants. This is the
exact situation of the Helmer-Papathanasiou-Tellander parachute example, which was found by the same Euler-drop
diagnostic.

Recipe, in increasing cost:
1. **Euler drop at a second rate point** (done at P2). Repeat at the base point and at one more rate set, scanning
   candidate z for a drop. Each count takes about 2-5 min in sympy.
2. **Get the polynomial.** Run HPT's `WhitneyStratifier` (Macaulay2; their paper's code) on the map
   {G=0} -> z-line with masses fixed. Or compute the Euler discriminant by homotopy continuation of the
   critical-point system (`HomotopyContinuation.jl`): track the 22 critical points along z and record where one
   escapes to the boundary (x_i -> 0/infinity or G -> 0). Applied at 3-4 rate sets, the escape location gives the
   rational function of the rates. PLD.jl (FMT) alone will not find it. Run it as a negative control.
3. **Regulator continuation.** Track the drop location z0(m5) as m5 -> 0. If z0(m5) -> -111/1120 smoothly, then
   identify the massive branch (a first-type point off the torus? the 1/100 test says not in the torus). If it
   appears only at m5 = 0, it is intrinsically a massless/IR-collinear stratum (Vergu's setting).
4. **Exponent check**: compare -7/2 with the second-type strength formula for the stratum once its face is known.

No closed formula in the rates was evaluated: the polynomial is not known yet, and guessing would be curve-fitting.

---

## Sources

arXiv:math/0510011, 0711.4863, 0910.0114, 1006.4064, 1304.5342, 1406.2664, 2008.10574, 0801.0891, 1810.07689,
2302.14840, 2307.11497, hep-th/9803091, 1811.04935, 2109.07584, 2311.16219, 2311.14669, 2109.08036, 2402.14787,
1308.6676, 1712.09215, 2512.12763, 1607.01571 (Adams-Bogner-Schweitzer-Weinzierl, kite to all orders in elliptic
polylogs, D=4-2eps [V]). Unverified: Fairlie et al. 1962 text; Kazakov 1983/84; D'Eramo-Peliti-Parisi 1971;
Vasil'ev et al. 1981; arXiv:1905.07380; the PRL title of 1810.07689.
