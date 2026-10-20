# GT label polytope and positive geometry: literature and corpus scan

Date: 2026-10-19. Read-only scan; this file is the only output. Exploratory (PI: "an itch").
Every polytope fact below was derived here and checked with a short sympy script, which
lived in the session scratchpad and is not kept. Literature status is marked
**[V]** (verified in the primary text or on the arXiv abstract page),
**[S]** (seen only in a search-engine summary), or **[C]** (verified only through a corpus bibitem).

## 1. Which polytope?

**The labels.** Shell n is the SO(4) irrep (n-1, 0) of the harmonics on S^3. Under the chain
SO(4) > SO(3) > SO(2), its Gel'fand-Tsetlin (GT) pattern is (l, m), with 0 <= l <= n-1 and
|m| <= l. Corpus anchor: Paper 0 l.824, Paper 7 l.120 and l.747, `tests/test_paper7_label_set.py`.

**Shell polytope (coupled, GT basis).** Put k = n-1. The shell polytope is
T_k = {(l,m) : l <= k, l-m >= 0, l+m >= 0} = k*T, where T = conv{(0,0), (1,1), (1,-1)}.
- **Dimension 2. It is a lattice triangle, so a 2-simplex.** Area k^2; there are 4k boundary
  points and (k-1)^2 interior points.
- **Ehrhart polynomial** L(k) = k^2 + 2k + 1 = (k+1)^2 = n^2, with h* = (1, 1) and normalized
  volume 2. Ehrhart reciprocity holds: L(-k) = (k-1)^2 is the interior count.
- **Integral: yes. Reflexive: no**, since T has no interior point. Among the dilates, only
  2T (shell n = 3) is reflexive: a lattice polygon is reflexive exactly when it has a single
  interior lattice point, here (1, 0). From n = 4 on, a dilate has (n-2)^2 >= 4 interior points.
- **Toric identification.** The map (l, m) -> (l, l+m) has determinant 1 and sends T to
  conv{(0,0), (1,0), (1,2)}. That triangle is unimodularly equivalent to conv{(0,0), (2,0), (0,1)},
  the moment polygon of the weighted projective plane **P(1,1,2)** (the quadric cone) with O(2).
  Check: the dimension of the degree-2k part of C[x,y,z] with deg z = 2 is (k+1)^2.
  The vertex (l,m) = (0,0) has edge vectors (1,1) and (1,-1), with determinant -2. **It is
  the polygon's only non-unimodular vertex (index 2), and it is the cone point of P(1,1,2).**
  The other two vertices are smooth.

**Shell polytope (uncoupled, SU(2)xSU(2) basis, which is also the Kustaanheimo-Stiefel (KS)
lift).** The labels (m1, m2) in [-j, j]^2 with j = (n-1)/2 fill a square: the moment polytope
of P^1 x P^1 = S^2 x S^2, with the same Ehrhart polynomial (k+1)^2.
- S^2 x S^2 is the smooth quadric Q^2 = SO(4)/P. H^0(Q^2, O(k)) is the space of degree-k
  harmonics, which is shell n = k+1. It is also the classical reduced Kepler space, with factors
  J+- = (L +- A)/2 [S, survey arXiv:2301.11759].
- With KS quanta (n1, n2, n3, n4) and the constraint n1+n2 = n3+n4 = n-1, the square is the
  mid-slice of the 4D-oscillator 3-simplex, which is circle reduction seen on polytopes.
- The square (4 vertices) and the triangle (3 vertices) share an Ehrhart polynomial but are not
  unimodularly equivalent.
- **So passing from the uncoupled basis to the GT basis is the toric degeneration
  Q^2 ~> P(1,1,2),** the textbook degeneration xy - zw -> xy - z^2 of the quadric surface.
- Symplectically, (|L|, L_z) on S^2 x S^2 maps onto the triangle. Over the singular vertex
  |L| = 0 the fibre is the antidiagonal 2-sphere J- = -J+, which is Lagrangian. This is derived
  here. It looks like the standard non-toric fibre of Gelfand-Cetlin systems, but no primary
  source was located: PROBABLY KNOWN.

**Truncated lattice up to n_max.** Put K = n_max - 1. The truncated lattice is
S_K = {(k,l,m) : 0 <= l <= k <= K, |m| <= l} = K*S, where
S = conv{(0,0,0), (1,0,0), (1,1,1), (1,1,-1)} is the pyramid over T. **It is a 3-simplex.**
- Ehrhart polynomial L(K) = K^3/3 + 3K^2/2 + 13K/6 + 1 = n_max(n_max+1)(2n_max+1)/6, as required.
- Ehrhart series (1+t)/(1-t)^4, so h* = (1, 1, 0, 0); normalized volume 2.
- Interior points: (K-1)(K-2)(2K-3)/6, which is 0 at K = 2 and 1 at K = 3 (n_max = 4).

**Is it a GT polytope in the toric-degeneration literature's sense?** In type A, GT polytopes
give the Gonciulea-Lakshmibai / Kogan-Miller toric degeneration of flag varieties. Kogan-Miller,
Adv. Math. 193 (2005) 1, arXiv:math/0303208 [S, title and abstract via search]. Alexeev-Brion
generalize this to spherical varieties ("Toric degenerations of spherical varieties",
arXiv:math/0403379 [S]). Other context:
- Littelmann, "Cones, crystals, and patterns", Transform. Groups 3 (1998) 145 [S]: in type A, GT
  patterns are string-cone patterns for a nice reduced word.
- Berenstein-Zelevinsky, Invent. Math. 143 (2001) 77, arXiv:math/9912012 [S]: string/BZ
  polytopes via tropicalization of totally positive varieties.
- Kiritchenko-Smirnov-Timorin, "Schubert calculus and Gelfand-Zetlin polytopes" [S, title only].
- Guillemin-Sternberg built the Gelfand-Cetlin integrable system for U(n) and for SO(n). This is
  known only via citing papers [S]; Harada, arXiv:math/0404485, says the same construction works
  for SO(n).

The GeoVac case is the orthogonal SO(4) > SO(3) chain on a *partial* flag variety (the quadric).
**I found no primary source that writes down the SO(4) shell triangle as a GT / Newton-Okounkov
polytope of Q^2.** It is elementary and almost certainly folklore, so PROBABLY KNOWN. The two
searches run found only type-A and general Newton-Okounkov material.

## 2. Positive geometry and canonical forms

**Polytopes are positive geometries.** Arkani-Hamed-Bai-Lam, "Positive Geometries and Canonical
Forms", JHEP 11 (2017) 039, arXiv:1703.04541 **[V: abstract page and extracted full text]**.
From the full text:
- Sec. 5.6 proves that positive parts of projectively normal toric varieties are positive
  geometries, and that the canonical form has a simple pole on each facet, with residue equal to
  the facet's canonical form.
- Sec. 7.3.2 and App. H obtain a polytope's canonical form as the push-forward of the toric
  variety's form under the algebraic moment map.
- Sec. 5.8 says that (G/P, (G/P)>=0) being a positive geometry "was essentially established in
  [19]", with all scalar normalizations only "expected" to equal 1.
- Lam, "An invitation to positive geometries", arXiv:2208.05407 [V: abstract only] is a short
  mathematical introduction.

**Applied to GeoVac.** The shell polytope and the truncated polytope are both **simplices**. A
simplex's canonical form is unique and has a constant numerator: the adjoint polynomial
(Warren's) is non-trivial only for non-simplicial polytopes. Explicitly:
- Omega(T_k) = 2k dl^dm / [(k-l)(l-m)(l+m)]
- Omega(S_K) = 2K dk^dl^dm / [(K-k)(k-l)(l-m)(l+m)]

In both, the prefactor is the determinant of the homogeneous facet matrix; signs follow
orientation. The uncoupled square's form is the product of two segment forms. **These forms
carry no information beyond the facet list.** Their only "values" are rational dual volumes.
Polytope canonical forms are rational by construction, so "rational skeleton, transcendentals
on integration" holds here trivially. It is not evidence for anything.

**Has anyone computed a GT polytope's canonical form or linked one to an amplitude?**
- An extended search ("Gelfand-Tsetlin polytope canonical form positive geometry") found
  nothing. It surfaced stringy canonical forms (arXiv:1912.08707), open associahedra
  (arXiv:2008.06418) and GT combinatorics, all unconnected [S].
- A second extended search (GT plus amplitudes or Feynman integrals) also found nothing.
- **Verdict: NOT FOUND.** The likely reason is that the nearby positive geometries (positive flag
  varieties) are studied directly, not through their GT degenerations.

## 3. Real bridges between representation-theory polytopes and amplitudes

**Real and cited:**
- **Associahedron and the bi-adjoint scalar.** Arkani-Hamed-Bai-He-Yan, "Scattering Forms and the
  Positive Geometry of Kinematics, Color and the Worldsheet", JHEP 05 (2018) 096,
  arXiv:1711.09102 [V: abstract]. The tree amplitude is the canonical form of an associahedron
  in kinematic space. The associahedron is a generalized permutohedron, not a GT polytope.
- **Cosmological polytopes.** Arkani-Hamed-Benincasa-Postnikov, arXiv:1709.02813 [V: abstract].
  Each Feynman graph's contribution to the wavefunction comes from a polytope.
- **Positive Grassmannian / positroid cells and total positivity.** Covered in ABL Secs. 5.5/5.8
  [V]. **The real common ancestor of amplitude geometry and crystal / canonical-basis polytopes
  is Lusztig's total positivity.** BZ 2001 [S] obtains string polytopes by *tropicalizing*
  totally positive varieties; amplituhedra live *inside* totally positive Grassmannians. This is
  a genuine shared root, not a map from one to the other.
- **MV polytopes.** Kamnitzer, Ann. Math. 171 (2010) 245, arXiv:math/0501365 [S]: MV polytopes
  are moment images of MV cycles and are in bijection with Lusztig's canonical basis. No
  amplitude link was found.
- **Feynman integrals and polytopes.**
  - Davydychev-Delbourgo, J. Math. Phys. 39, 4299 (1998), hep-th/9709216 [C: Paper 59 bibitem,
    sec:density_ring]: one-loop N-point integrals as volumes of non-Euclidean simplices.
  - Arkani-Hamed-Hillman-Mizera, arXiv:2202.12296, "Feynman polytopes" = Newton polytopes of
    Symanzik polynomials [S].
  - de la Cruz-Novichkov-Vanhove, arXiv:2512.10518 [V: abstract; to appear in JHEP]: classify
    Fano and reflexive polytopes from quasi-finite Feynman integrals via (bivariate) Ehrhart
    counts. They find reflexivity rare and tied to Calabi-Yau periods.
  - None of these polytopes is a representation-label polytope.
- **Lattice-point sums versus integrals.** Brion-Vergne, "Lattice points in simple polytopes",
  JAMS 10 (1997) 371 [S]. Khovanskii-Pukhlikov is the unimodular case [S, via Brion-Vergne's
  abstract]. For polynomial weights on a simple lattice polytope, the lattice sum is a Todd
  operator applied to the integral over facet-shifted polytopes, an exact rational identity.
  T_k and S_K are simple, but T has the index-2 vertex, so the Brion-Vergne form is needed
  rather than Khovanskii-Pukhlikov.

**Analogy only, no citation behind it.** "GT patterns <-> positroid cells" and "crystal bases <->
amplituhedron cells". Both connect only through total positivity, as stated above.

## 4. Corpus

The grep over `papers/` (archive excluded) and `docs/` for Gelfand / Tsetlin / Cetlin / Zetlin /
Ehrhart / polytope / canonical form / positive geometry / amplituhedron returned:

- **The GT label-set statement**, made identically in Paper 0 l.824, Paper 7 l.79/120/747,
  Paper 25 l.82 and the group3 synthesis l.154/476: the nodes are the GT labels of the
  Peter-Weyl truncation of L^2(S^3), with n^2 per shell. It is backed by
  `tests/test_paper7_label_set.py` (claim_test_matrix row 787, BACKED-SOUND). **No corpus
  document treats the label set as a polytope, gives an Ehrhart polynomial, or mentions positive
  geometry, canonical forms (in the ABL sense) or the amplituhedron.**
- **Only near-hit:** `docs/track_nk_bargmann_segal_sprint1.md` l.200 says the HO shell labels are
  a simplex (the SU(3) GT pattern) that "replaces the paraboloid".
  - In polytope terms both shells are 2-simplices. The HO shell is the standard simplex dilate:
    unimodular, P^2. The Coulomb shell is the P(1,1,2) triangle, with one index-2 vertex.
  - In the shared (l, m) chart, the HO shell is the points of the same triangle T_N with
    l = N (mod 2).
- **Unrelated "canonical form" hits:** Paper 59 l.349 (a normal form), Paper 18 l.2885, and
  C23 (Halmos).
- **Lattice-sum ingredient, relevant to Sec. 5:** Paper 2's B(n_max) = Sum (2l+1) l(l+1) is a
  polynomial-weight lattice sum over S_K. The selection principle B/N = 3 (Paper 2
  Sec. "Selection principle") is the *lattice-point average* of l(l+1) over S_K, which equals 3
  at K = 2. The continuum (volume) average over S_2 is 11/5 (computed here); Brion-Vergne's Todd
  operator relates the two exactly.

## 5. Verdict

**Plainly: at the level of the label polytope, nothing with an unknown outcome is
well-defined.** The polytope is a simplex in both the shell and the truncated form. A simplex
is the trivial positive geometry: its canonical form is fixed by its facets, its numerator is
constant, and every "volume" it yields is rational. The facts in Sec. 1 (Ehrhart polynomial,
P(1,1,2) identification, index-2 vertex at l = 0, the toric degeneration Q^2 ~> P(1,1,2), the
reflexive dilate at n = 3) are all now computed. A pre-registration of them would be predicting
known outcomes, which docs/predictions/README.md rules out.

**Numerology hazard, flagged so nobody chases it.** 2T, shell n = 3, is the unique reflexive
shell dilate, and Paper 2's selection picks n_max = 3. The two statements concern different
objects: one shell versus the cumulative tetrahedron, whose single-interior-point dilate falls
at n_max = 4, not 3. No mechanism connects them. Under `feedback_audit_numerical_claims.md` this
is a coincidence of the integer 3, not a lead.

Two candidates, offered only because the brief asked, with honest labels.

**C1 -- bookkeeping CHECK (an atlas entry, not a finding).**
- **Statement:** the GT shell polytope is (n-1)*T, unimodularly the moment triangle of
  P(1,1,2), with Ehrhart polynomial n^2, h* = (1,1), and exactly one non-unimodular vertex (the
  s-state vertex l = m = 0, index 2). The uncoupled / KS label set is the S^2 x S^2 square, and
  the basis change is the toric degeneration Q^2 ~> P(1,1,2). The truncated set is a 3-simplex
  with h* = (1,1,0,0).
- **Falsifier:** a lattice-count mismatch at any n <= 50, a second singular vertex, or a GT
  label outside T_k.
- **Novelty:** PROBABLY KNOWN. Two searches turned up no explicit SO(4) source, but the
  Kogan-Miller / Alexeev-Brion framework and the quadric-cone degeneration make it folklore.
- **Cost:** one test file.
- **Value:** it would let the atlas record the label set as FORCED with an explicit
  toric/positive-geometry type. That closes the itch rather than opening a programme.

**C2 -- the only route to non-trivial canonical-form content; conditional, and currently NOT
well-defined.**
- The spinor (Dirac) label set (n, kappa, m_j) is the first GeoVac label set whose natural hull
  is not a simplex. In a chart such as (l, 2(j-l), 2m_j) it is a stacked pair of triangles with
  5-6 vertices. A non-simplicial polytope has a non-constant adjoint polynomial, so its canonical
  form would carry real content: for example, whether the adjoint's zero locus is SO(4)-forced.
- **Blocker:** embedding half-integer labels in a lattice is a convention choice, and the
  polytope and adjoint change with the chart. **Not pre-registerable until a
  representation-forced chart is fixed.** The obvious candidate is the GT pattern of Spin(4) on
  the spinor bundle, for the CH Dirac eigenspaces of Paper 32/38.
- **Falsifier, once a chart is fixed:** the polytope is a simplex or product of simplices
  (trivial adjoint), or its adjoint depends on the choice of chart.
- **Novelty:** NOT FOUND. No search found GT or label polytopes of spinor representations in a
  positive-geometry context, though that search was not run separately.
- **Recommendation:** pursue it only if the PI wants the itch scratched at the Dirac level.

**Bottom line.** The GeoVac label set is the lattice-point set of a simplex: the dilated
P(1,1,2) triangle per shell, which is the GT toric degeneration of the quadric
S^2 x S^2 = reduced Kepler space, and the pyramid over it cumulatively. Polytope positive
geometry is real (ABL 2017). Amplitude bridges are real for associahedra, cosmological polytopes
and the positive Grassmannian, but none passes through GT or label polytopes. The only genuine
meeting point in the literature is Lusztig total positivity as a shared root. At the spinless
level there is nothing to test.
