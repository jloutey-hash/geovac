# Prior art: tori, rings and multiple quantum-number sets (math / HEP side)

Date: 2026-11-01. Read-only scan for the PI's ring/nanotorus prediction. Downloads are in
`debug/lit_scan/dl_prior_ring_hep/` (git-ignored, confirmed). None of the items below is cited
anywhere in `papers/` (grep, 2026-11-01). Manton 1987, Park 1960, Tarter 1970 and the
Mardoyan et al. book are already cited in Paper 25, so they are not re-reported here.
Confidence tags: VERIFIED = read in the downloaded text; ABSTRACT = abstract or bibliographic
record only; SECONDARY = reported by another source.

## Headline

**The prediction already has an exact, published version on S³ itself.** The free Laplacian on S³
separates in six orthogonal coordinate systems. Each one gives a different complete set of quantum
numbers on the same eigenspace, and one of them, the Hopf ("cylindrical") system, labels states by
two winding numbers on the Clifford-torus family. Paper 25's two groupings (by rotation, and by the
Hopf/chiral factor) are two of these six. A warning comes with it: in one of the six (prolate),
quantum numbers **cannot be assigned globally** (quantum monodromy), and this has been shown for
hydrogen itself.

## 1. Separability in toroidal coordinates

- **C. P. Boyer, E. G. Kalnins, W. Miller Jr., "Symmetry and separation of variables for the
  Helmholtz and Laplace equations," Nagoya Math. J. 60, 35–80 (1976), doi:10.1017/S0027763000017165.**
  VERIFIED. Quotes: "(0.1) separates in eleven coordinate systems"; "there are 17 types of cyclidic
  coordinate systems which permit R-separation of variables in the Laplace equation"; "The
  existence of these additional systems is related to the conformal symmetry of (0.2)"; and "all
  known R-separable orthogonal coordinate systems for (0.2) are characterized by a pair of commuting
  second-order symmetry operators in the enveloping algebra of so(4,1)." This is the primary source
  for the conformal (SO(4,1)) origin of R-separation. Toroidal coordinates are one of the 17
  cyclidic systems. Moon–Spencer (*Field Theory Handbook*, Springer 1961) and Morse–Feshbach
  (*Methods of Theoretical Physics*, 1953) are the classic tables. I did not open them. That they
  list toroidal coordinates as R-separable for Laplace but not for Helmholtz is SECONDARY (via
  Pereira–Schmidt and Wikipedia).
- **Boyer, Kalnins, Miller, "R-separable coordinates for three-dimensional complex Riemannian
  spaces," Trans. AMS 242, 355–376 (1978).** ABSTRACT. This is the group-theoretic classification.
- **R. Prus, A. Sym, Phys. Lett. A 336, 459–462 (2005); arXiv:nlin/0402035.** ABSTRACT: proves that
  cyclidic coordinates can be "non-regular (non-Stäckel) R-separable" for Helmholtz.
- **M. E. Pereira, A. G. M. Schmidt, "Schrödinger equation is R-separable in toroidal coordinates,"
  arXiv:2511.08646 (2025), a preprint.** VERIFIED (text). It claims "for the first time, exact
  solutions for the Schrödinger equation in Moon and Spencer's toroidal coordinates" using
  *irregular* R-separation, which gives "fractional angular momentum eigenvalues" ("because of the
  torus's hole"). It concedes that ordinary separability is still impossible ("This is still
  true"). It is not peer-reviewed, and fractional angular momentum is a red flag for a
  single-valuedness problem. Treat it as unvetted.
- **L. G. Mardoyan et al. (already cited), arXiv:2310.17336, Introduction.** VERIFIED. Hydrogen
  separates in "four coordinate systems, namely spherical, sphero-conal, parabolic and prolate
  spheroidal." **Toroidal coordinates are not among them.** No Coulomb problem is known to separate
  in R³ toroidal coordinates.
- **Problems that are separable on a ring or torus:** a particle on a 1D ring (trivially, with
  Aharonov–Bohm flux), on a flat 2D torus, or on an embedded torus surface of revolution, where the
  toroidal angle separates exactly and the poloidal equation carries a curvature potential. Encinosa
  & Mott, PRA 68, 014102 (2003) treat the last case numerically (SECONDARY, via Pereira–Schmidt
  refs). These are constrained-surface problems, not 3D Coulomb problems.

**Elementary check (computed here, not a literature claim):** stereographic projection from S³ maps
Hopf coordinates (x₁+ix₂ = cos ψ e^{iφ₁}, x₃+ix₄ = sin ψ e^{iφ₂}) onto R³ toroidal coordinates,
with tanh η = cos ψ and azimuth = φ₁. Verified numerically to 1e-10 at three values of ψ. The
Clifford torus is ψ = π/4. The R-factor √(cosh η − cos ξ) of toroidal Laplace separation is the
conformal weight, which is why separation needs a *conformally* invariant equation.

**Caveat for the side idea (load-bearing).** Fock's S³ is the stereographic image of **momentum**
space, not position space. A ring of nuclei in position space therefore does *not* become a circle
on Fock's S³. The position-space picture (a ring in R³ maps to a circle in S³, toroidal coordinates
are Hopf coordinates) works for Laplace and conformal problems, not for the Coulomb Schrödinger
equation. That said, a two-centre field on Fock's S³ keeps only L_z (SY result, v5.21.26). A
C_N-symmetric ring keeps a discrete C_N ⊂ that SO(2), a finite subgroup of one of the two torus
U(1)s.

## 2. Hopf tori, Clifford torus, harmonic analysis on the maximal torus

- **S. R. Dawson, H. R. Dullin, "Quantum integrable systems arising from separation of variables on
  S³," arXiv:2405.08778 (2024).** VERIFIED. This is the most relevant find. "All 6 separable
  orthogonal coordinates on the 3-sphere: ellipsoidal, prolate, oblate, Lamé, spherical and
  cylindrical"; "Cylindrical (or Hopf, or doubly cylindrical) coordinates"; the cylindrical system
  is "a toric system on S²×S² with a square as Delzant polytope", and the separated equation is
  Jacobi. Also: "For the prolate system ... there exists a defect in the joint spectrum which
  prohibits a global assignment of quantum numbers." The space of separable systems is a Stasheff
  pentagon, and the six labellings are its vertices, edges and interior.
- **H. R. Dullin, H. Waalkens, "Defect in the joint spectrum of hydrogen due to monodromy," PRL 120,
  020507 (2018); arXiv:1612.00823.** ABSTRACT. Hydrogen in prolate spheroidal separation shows
  quantum monodromy, so quantum numbers cannot be assigned globally near threshold. This is the
  two-centre (Coulson–Joseph) labelling, directly relevant to multi-centre rings.
- **U. Pinkall, "Hopf tori in S³," Invent. Math. 81, 379–386 (1985), doi:10.1007/BF01389060.**
  SECONDARY quote: "the inverse image of any closed curve on S² will be an immersed torus in S³
  which will be called a Hopf torus." The geometry is classical: any loop on the Hopf base lifts to
  a torus. This is a natural "ring ↔ torus" map, but it is pure geometry with no spectral content.
- **A. F. Rañada, "A topological theory of the electromagnetic field," Lett. Math. Phys. 18, 97–106
  (1989).** ABSTRACT: the Hopf index equals the linking number of field lines. **W. T. M. Irvine,
  D. Bouwmeester, "Linked and knotted beams of light," Nature Phys. 4, 716–720 (2008).** SECONDARY:
  the stereographically projected Hopf map gives field lines on nested tori. **L. D. Faddeev,
  A. J. Niemi, Nature 387, 58 (1997)** (knotted solitons, Hopfions) is UNVERIFIED in this scan and
  is cited from memory only.
- Harmonics on the maximal torus U(1)×U(1) ⊂ SO(4) are the (m_L, m_R) weight basis. On Fock's S³
  this is the hydrogen parabolic basis. The identification "parabolic basis = Hopf/cylindrical
  coordinates on S³" is standard group theory (Park's addition-of-angular-momentum argument,
  VERIFIED as described in Mardoyan Introduction). I found no single source that states it in the
  Hopf-coordinate language. That wording is my inference, consistent with Dawson–Dullin's square
  polytope.

## 3. Inequivalent quantum-number sets on a torus

- **K. Sasaki, Y. Kawazoe, Prog. Theor. Phys. 112, 369 (2004); arXiv:cond-mat/0307339.** VERIFIED.
  Nanotorus π-states carry two labels (around and along the tube axis). "Wave vectors k and k + Δk
  are equivalent" under a reciprocal lattice of the chiral and translational vectors, so the integer
  pair is a *choice* of labelling. A twist "can be thought of as an effective gauge field," which can
  halve the flux period. This is the closest condensed-matter analogue of the PI's prediction, and
  it is already worked out.
- **C. T. White, D. H. Robertson, J. W. Mintmire, PRB 47, 5485 (1993).** ABSTRACT: classifies tubules
  by "helical and rotational symmetries". The same states can be labelled by a screw (helical)
  quantum number or by a pure rotation plus translation. These are two valid groupings, which is an
  exact analogue of "two groupings of the same states".
- **A. Ceulemans, L. F. Chibotaru, S. A. Bovin, P. W. Fowler, J. Chem. Phys. 112, 4271 (2000).**
  ABSTRACT: closed-form tight-binding spectra of polyhex tori. Metallic iff chiral and twist indices
  ≡ 0 mod 3. **M. F. Lin, D. S. Chuu, PRB 57, 6731 (1998)** (page from memory; volume and issue
  ABSTRACT): persistent currents in nanotori, hc/e periodic.
- **N. Byers, C. N. Yang, PRL 7, 46 (1961).** SECONDARY: single-valuedness restricts gauge
  transformations to integer flux quanta. This is the large-gauge, winding-number origin of
  ring flux periodicity. **Physics Letters A 96, 365 (1983), "Josephson behavior in small normal
  one-dimensional rings"** (normally attributed to Büttiker, Imry, Landauer). The bibliographic
  record is VERIFIED, but the author list was not confirmed on the DOI page.
- **G. 't Hooft, Nucl. Phys. B 153, 141 (1979)**, twisted boundary conditions on the torus.
  SECONDARY: the flux sectors are integers mod N. **A. Giveon, M. Porrati, E. Rabinovici, Phys. Rep.
  244, 77 (1994), hep-th/9401139.** ABSTRACT: T-duality exchanges momentum and winding on a circle
  or torus. This is the cleanest HEP statement that "two inequivalent quantum-number sets label one
  spectrum", with an O(d,d;ℤ) or SL(2,ℤ) group relabelling the lattice of integer pairs.

**What "multiple valid quantum numbers on a torus" means mathematically.** The torus labels form a
lattice ℤ², and any GL(2,ℤ) change of basis cycles relabels the same states. A *physical*
distinction between labellings appears only when the Hamiltonian (lattice, twist, flux) picks a
sublattice or a direction. For example, a C_N ring makes angular momentum m defined only mod N
(Umklapp).

## 4. Toroidal plasma (low priority)

- Safety factor q (rotational transform ι = 1/q) is a winding number, and rational surfaces ι = m/n
  carry closed field lines and islands (standard, SECONDARY). MHD modes are labelled by the pair
  (m, n).
- **A. Thompson, J. Swearngin, A. Wickes, D. Bouwmeester, PRE 89, 043104 (2014); arXiv:1310.4229.**
  ABSTRACT: plasma configurations "defined by toroidal and poloidal winding numbers"; (1,1) is the
  Kamchatnov–Hopf soliton (**A. M. Kamchatnov, Sov. Phys. JETP 55, 69 (1982)**, ABSTRACT); field lines
  lie on nested tori. This is a genuine Hopf-fibration plasma object.
- Quasi-symmetry of stellarators: |B| depends on Mθ − Nφ, so the choice of which torus U(1) is the
  symmetry is the plasma's version of a "grouping choice" (Landreman & Paul, PRL 128, 035001 (2022),
  arXiv:2108.03711; Helander, Rep. Prog. Phys. 77, 087001 (2014); both ABSTRACT).
- **Verdict: shared mathematics only.** No GeoVac-style discrete label skeleton (Casimir, Gaunt rules)
  appears in plasma physics. The winding numbers there are field-line topology, not eigenvalue labels.
  No predictive bridge.

## Bearing on the ring study

1. Pre-register the prediction against the S³ benchmark: the six Dawson–Dullin labellings already
   give "more than one valid set" for the single-centre problem. A ring result is new only if it
   shows something *beyond* S³'s multiseparability: a labelling forced by C_N, or one broken by it.
2. Expect only one label to survive for a planar C_N ring: the ring index k mod N, a discrete
   subgroup of the surviving L_z. The tube or second winding label has no exact symmetry
   counterpart in a finite molecule. It may survive as a *grouping* (like Paper 25's) or as an
   approximate band index (Sasaki–Kawazoe), not as a conserved quantity.
3. Check for monodromy. Multi-centre labellings can fail globally (Dullin–Waalkens), and a
   pre-registered "labels assignable everywhere" clause would be falsifiable.
4. The "ring → circle on S³" idea is valid for position-space conformal problems but not for
   Fock's momentum S³. State this before building on it.
