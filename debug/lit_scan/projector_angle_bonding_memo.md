# Adversarial prior-art scan: projector angles / commutator / compound matrices for bonding

**Date:** 2026-09-06 · **Type:** LITERATURE SCAN (no code/paper edits) · **Scope:** the three
projector-algebra claims in Paper 60 §manyelectron + Paper 58 (decompactification front) + Paper 32
(Rem. `rem:four_center_modulus`), assessed against EDA / localization / nonorthogonal-CI literature.

Confirmed live: every citation below was checked via WebSearch/WebFetch in this session (2026-09-06);
none is asserted from memory alone. No arXiv ID/DOI/title is invented — anything I could not confirm is
marked UNVERIFIED, and there are none left unmarked below.

## Ranked references, by claim

### Claim 1 — principal (canonical) angles between two atoms' occupied subspaces as a bonding coordinate

1. **Jordan, C., "Essai sur la géométrie à n dimensions," Bull. Soc. Math. France 3, 103 (1875).**
   Origin of principal/canonical angles between subspaces (pure math, pre-dates chemistry entirely).
2. **King, H.F.; Stanton, R.E.; Kim, H.; Wyatt, R.E.; Parr, R.G., "Corresponding Orbitals and the
   Nonorthogonality Problem in Molecular Quantum Mechanics," J. Chem. Phys. 47, 1936 (1967),
   doi:10.1063/1.1712221.** Proves any two nonorthogonal orbital sets can be rotated (independently on
   each side) to a *diagonal* overlap — i.e. SVD of the inter-set overlap block, which is exactly
   cos(principal angles). This IS the corresponding-orbital transformation (generalizing Amos–Hall) and
   is the direct ancestor of every "SVD between two orbital spans" method below. Claim 1's core linear
   algebra.
3. **Amos, A.T.; Hall, G.G., "Single Determinant Wave Functions," Proc. R. Soc. Lond. A 263, 483 (1961),
   doi:10.1098/rspa.1961.0175.** First "corresponding orbitals" construction (diagonalizing a rectangular
   overlap by two unitary transformations = principal angles) for single-determinant wavefunctions.
4. **West, A.C.; Schmidt, M.W.; Gordon, M.S.; Ruedenberg, K., "A comprehensive analysis of
   molecule-intrinsic quasi-atomic, bonding, and correlating orbitals. I. Hartree-Fock wave functions,"
   J. Chem. Phys. 139, 234107 (2013), doi:10.1063/1.4840776.** The valence MO space is resolved by SVD
   against the free-atom minimal-basis space to build "oriented quasi-atomic orbitals" — i.e. principal
   angles between the *occupied molecular* subspace and a *free-atom* subspace, used exactly as a bonding
   diagnostic (bonding/lone-pair/radical character read off the resulting angles/occupations). Closest
   single match to "principal angle between an atom-like subspace and the molecular occupied space
   measures how shared/bonded it is" — but static, at one geometry, not tracked continuously across R.
5. **Wu, Q. et al., "Fragment Aligned Molecular Orbital Analysis" (FAMO), J. Chem. Theory Comput.
   (2024), doi:10.1021/acs.jctc.4c00456; PubMed 39046803.** Explicitly uses orthogonal-Procrustes
   alignment (SVD of the cross-overlap between the molecule's occupied MOs and the constituent
   fragments' occupied MOs) — Procrustes alignment IS the principal-angle machinery restated as a
   rotation problem. Confirmed 2024, active research line.
6. **"Natural Fragment Bond Orbital" (NFBO) method, J. Am. Chem. Soc. (2024/2025),
   doi:10.1021/jacs.4c12421.** Derives interfragment bonding/antibonding orbitals from canonical MOs;
   same family (fragment-subspace decomposition of the occupied space).
7. **Khaliullin, R.Z.; Cobar, E.A.; Lochan, R.C.; Bell, A.T.; Head-Gordon, M., "Unravelling the Origin of
   Intermolecular Interactions Using Absolutely Localized Molecular Orbitals," J. Phys. Chem. A 111,
   8753 (2007), doi:10.1021/jp073685z.** ALMO-EDA: partitions the interaction into frozen/polarization/
   charge-transfer by comparing block-localized (fragment-restricted) vs. fully relaxed occupied spaces —
   conceptually the same "how much does letting the two atoms' occupied spaces mix change things"
   question, but measured energetically, not via an explicit principal-angle/R* law.
8. **"The Role of Principal Angles in Subspace Classification," arXiv:1507.04230 (2015).** Confirms the
   general-purpose statistics/ML fact used above: principal angles quantify subspace separation
   (θ_k = π/2 for all k ⇔ maximal separation) — standard tool, unconnected to internuclear-distance
   bonding physics in this source.

**Not found anywhere searched:** a paper that (i) tracks the principal angle between two ATOM-CENTERED
occupied subspaces as an explicit function of internuclear distance R, (ii) identifies the 45° crossing
as coincident with peak `||[P_A,P_B]||`, or (iii) states a decay-length scaling law R* ~ c·n/Z for that
crossing. The mathematical machinery (SVD of the cross-block overlap = principal angles) is 60 years old
and is the explicit working tool in Ruedenberg/QUAO, King et al., and the 2024 FAMO/Procrustes line; the
*R-dependent front / united-atom↔separated-atom decompactification* reading of it is not attested in the
literature this scan reached.

### Claim 2 — the commutator [P_A, P_B] as the exact obstruction to factoring a molecule into a product of atoms

9. **Halmos, P.R., "Two subspaces," Trans. Amer. Math. Soc. 144, 381 (1969).** Foundational operator-
   theory result: any two closed subspaces (equivalently, their orthogonal projections) admit a
   2×2-block "generic position" decomposition parameterized exactly by the principal angles; [P,Q] = 0
   iff the subspaces are in a compatible (simultaneously block-diagonalizable) position.
10. **Böttcher, A.; Spitkovsky, I.M., "A gentle guide to the basics of two projections theory," Linear
    Algebra Appl. 432, 1412 (2010), doi:10.1016/j.laa.2009.11.002.** Modern survey of exactly this
    algebra — the commutator, anticommutator, and norm of functions of two orthogonal projections, all
    expressed through principal angles. This is where the *pure-math* identity behind claim 2
    (`||[P_A,P_B]||` as a function of the σ_k) already lives, decades before any chemistry application.
11. **Manby, T.F.; Miller, T.F. III et al., projector-based embedding — e.g. Lee, S.J.R.; Welborn, M.;
    Manby, F.R.; Miller, T.F. III, "Projection-Based Wavefunction-in-DFT Embedding," Acc. Chem. Res. 52,
    1359 (2019) [confirmed via search of the Manby–Miller level-shift projector P_B = μ S_AB D_B S_BA];
    Wesolowski, T.A.; Warshel, A., frozen-density embedding, J. Phys. Chem. 97, 8050 (1993).** These are
    the chemistry-native statements of the SAME fact: fragment subsystems fail to factor additively
    (non-additive kinetic energy / need for an explicit orthogonality-enforcing projector) *because*
    S_AB ≠ 0, i.e. P_A and P_B do not commute. They diagnose and fix the non-additivity directly through
    the overlap block S_AB and a level-shift projector, not through an explicit `||[P_A,P_B]||` operator
    norm or an "exact obstruction" theorem statement.
12. **Density Matrix Embedding Theory (DMET) non-additivity discussion** (general fragment-embedding
    literature; e.g. reviews of DMET/quantum embedding surfaced this session note "the non-additive
    kinetic term is conspicuously absent" when fragmenting) — same phenomenon, same informal diagnosis
    (orbitals on different fragments overlap ⇒ energy doesn't add), never phrased as a commutator-norm
    theorem.

**Verdict basis:** the underlying linear algebra (non-commuting projectors ⇔ non-orthogonal/overlapping
subspaces ⇔ additive/product decomposition fails) is exactly Halmos's 1969 two-subspaces theorem and is
implicit in every embedding/EDA method that has ever had to deal with S_AB ≠ 0. No chemistry source
found states it explicitly as "`||[P_A,P_B]||` is *the* exact obstruction to factoring the molecular
Hilbert space/wavefunction into a tensor product of atomic factors" — that specific operator-algebraic
packaging (as opposed to working directly with S_AB or a level-shift projector) looks like a genuine,
if small, re-framing rather than new mathematics.

### Claim 3 — k-electron overlap = k-th compound matrix of the one-electron overlap; bond order/collapse under correlation

13. **Löwdin, P.O., "Quantum Theory of Many-Particle Systems. I. Physical Interpretations by Means of
    Density Matrices, Natural Spin-Orbitals, and Convergence Problems in the Method of Configuration
    Interaction," Phys. Rev. 97, 1474 (1955).** Establishes that the overlap of two N-electron Slater
    determinants is det(S) with S the one-electron overlap block — the N×N (top) compound matrix — and
    lays the groundwork for the general cofactor/minor calculus used by every nonorthogonal-CI method
    since.
14. **Prosser, F.; Hagstrom, S., "Simplified evaluation of matrix elements...," Int. J. Quantum Chem. 2,
    89 (1968).** The classical "cofactor"/compound-matrix method for matrix elements between
    nonorthogonal Slater determinants — k-electron quantities built from k×k minors of the one-electron
    overlap, generalized to density matrices of all orders ("super-cofactor" strategy). This is the
    textbook precedent for claim 3's central identity.
15. **King, Stanton, Kim, Wyatt, Parr (1967), op. cit. (ref. 2).** Same paper: once orbitals are rotated
    to corresponding-orbital form (diagonal overlap), the k-electron overlap collapses to a product of
    the k largest singular values — literally the diagonal compound-matrix statement, used to generalize
    the Slater–Condon rules.
16. **Burton, H.G.A., "Generalized nonorthogonal matrix elements: Unifying Wick's theorem and the
    Slater–Condon rules," J. Chem. Phys. 154, 144109 (2021), arXiv:2101.10944, and Part II, J. Chem.
    Phys. 157, 204109 (2022), arXiv:2208.10208.** Modern, general treatment explicitly built on cofactors
    and **compound matrices** (exterior powers) of the one-electron overlap for arbitrary-order
    nonorthogonal matrix elements — confirms the identity is standard, current, and named "compound
    matrix" in the literature (not just an equivalent object under another name).
17. **Mayer, I., "Charge, bond order and valence in the ab initio SCF theory," Chem. Phys. Lett. 97, 270
    (1983); "Bond order and valence indices: a personal account," J. Comput. Chem. 28, 204 (2007),
    doi:10.1002/jcc.20494.** Bond order M_AB = Σ(DS)_μν(DS)_νμ built from the one-electron density and
    overlap — the standard "bond order from overlap" object claim 3 says "drops straight out"; well
    known this is a functional of the one-electron overlap/density, though Mayer's papers do not
    themselves invoke compound-matrix/exterior-algebra language.
18. **Wiberg, K.B., Tetrahedron 24, 1083 (1968).** Original overlap-squared bond index (semiempirical
    precursor to Mayer's ab initio generalization).

**Verdict basis:** "the k-electron overlap is the k-th compound matrix (k×k minors) of the one-electron
overlap" is not a new mathematical fact — it is the Cauchy–Binet identity applied to Slater determinants,
textbook since Löwdin (1955)/Prosser–Hagstrom (1968)/King et al. (1967), and is the explicit modern
vocabulary of Burton (2021). What is genuinely being *added* in our framing (per the decompactification
memo) — that natural-orbital occupation numbers reweight the compound-matrix entries so that a small
antibonding occupation n_u disproportionately pulls the signed coherence toward the Heitler–London limit
— is a specific quantitative reading of that old identity, not the identity itself.

## Verdicts

| Claim | Verdict | Why |
|---|---|---|
| **(1) Principal angles between two atoms' occupied subspaces, tracked over R, 45° = peak commutator = decompactification front, R* ~ 2.14 n/Z** | **ADJACENT** | The linear algebra (SVD of cross-block overlap = principal angles) is 60-year-old machinery, actively used today for essentially this purpose (Ruedenberg/QUAO refs 4; FAMO/Procrustes ref 5, 2024). But no source found tracks it as a continuous function of R with an explicit 45°/peak-commutator/decay-length law tying it to a united-atom↔separated-atom transition. The *tool* is COLLISION-grade familiar; the *R-dependent front reading* of the tool is not attested. |
| **(2) `‖[P_A,P_B]‖` as THE exact obstruction to factoring a molecule into a product of atoms** | **ADJACENT** | The pure math (Halmos 1969; Böttcher–Spitkovsky 2010) is COLLISION — completely standard operator theory, older than the project. The chemistry-native statements of the same fact (Manby–Miller projector embedding, frozen-density non-additive kinetic energy, DMET non-additivity) work directly with S_AB rather than an explicit commutator-norm "obstruction" theorem. Packaging as an operator-algebra statement (rather than an overlap-matrix / non-additive-functional statement) looks like a genuine, modest re-framing — not new mathematics, and not previously phrased this way in the chemistry literature reached by this scan. |
| **(3) k-electron overlap = k-th compound matrix of one-electron overlap ⇒ bond order + its correlation collapse** | **COLLISION** | Textbook since Löwdin (1955) and Prosser–Hagstrom (1968); the corresponding-orbital form is King et al. (1967); the exact term "compound matrix" for this object is current literature (Burton 2021/2022). The identity itself must be cited, not re-derived as new. Only the specific natural-orbital-occupation-weighting argument connecting it to the decompactification front (Paper 60 §manyelectron / the 2026-09-06 memo) is a fresh application. |

## The blunt line

Two of three claims rest on genuinely standard machinery that must be cited, not presented as new:
principal angles/SVD-of-overlap for fragment bonding (Amos–Hall 1961, King et al. 1967, Ruedenberg/QUAO
2013, Procrustes-FAMO 2024) and the compound-matrix/exterior-power identity for k-electron overlaps
(Löwdin 1955, Prosser–Hagstrom 1968, Burton 2021/2022) are both mature, named, and actively published.
The commutator-norm framing of the "composition wall" (claim 2) sits on equally standard pure math
(Halmos 1969; Böttcher–Spitkovsky 2010) that the embedding/EDA literature already uses implicitly via
S_AB, without stating it as an explicit operator-norm obstruction theorem. **What looks genuinely ours is
narrow and specific: (a) reading the SVD/principal-angle machinery as a continuous function of
internuclear distance R and identifying the 45° crossing with the peak of `‖[P_A,P_B]‖` and a measured
decay-length law R* ≈ 2.14·n/Z (Papers 58/60, CHANGELOG v5.10.2); (b) the explicit "obstruction" framing
of the commutator as *the* factorization obstruction, stated as an operator identity rather than worked
through S_AB; and (c) the natural-orbital-occupation argument for why correlation moves the compound-
matrix-based coherence front but not the principal-angle front.** None of that is a new mathematical
object — every piece of hardware (SVD of overlap, compound matrices, projector commutators) must be
attributed to the literature above. The papers should cite refs 2–4, 9–10, and 13–16 explicitly rather
than presenting principal angles / projector commutators / compound matrices as GeoVac inventions.

## Search log (queries run, 2026-09-06)

Principal angles + bonding/dissociation; ALMO-EDA; IBO (Knizia); QUAO (Ruedenberg); Mayer/Wiberg bond
order; compound matrix/exterior power + Slater determinants; non-commuting projectors + fragments;
Grassmann/subspace distance + molecular similarity; Amos–Hall / King–Stanton–Kim–Wyatt–Parr corresponding
orbitals; principal angles + UA/SA dissociation limit; embedding-projector commutator/non-additivity;
Löwdin 1955; Manby–Miller frozen-density projector; Böttcher–Spitkovsky two-projections algebra; Burton
2021/2022 generalized nonorthogonal matrix elements; Prosser–Hagstrom cofactors; FAMO/NFBO/Procrustes
fragment-orbital alignment (2024); avoided-crossing/diabatization + principal angles. No arXiv ID, DOI,
or title in this memo was invented; all were confirmed via live search/fetch this session.
