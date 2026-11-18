# Cosmological-polytope match: three-vertex chain check

Date 2026-11-18. Follow-up to `2026-11-18_positive_geometry_memo.md`. Script
`debug/posgeom_chain_check.py` (runs in under a minute), data `debug/data/posgeom_chain_check.json`.
No paper, record or git state touched.

## Verdict

The three-site case holds exactly. It is the same statement as Paper 63's shell-theorem radial
recursion, written in the wavefunction-coefficient vocabulary, and the rational closed forms it gives
are textbook atomic physics. The dictionary is worth at most one credited remark in Paper 63; nothing
for Paper 58.

## 1. Two-site identity, normalisation fixed

With the full measures (no extra 1/r factors):

    I2(a,b;mu) = int int d^3r1 d^3r2 e^{-a r1} e^{-b r2} e^{-mu r12}/r12 = 16 pi^2 d_a d_b psi2,
    psi2 = 1/[(a+b)(a+mu)(b+mu)]                      (1709.02813 Eq. (13)/(22))
         = 32 pi^2 (a^2+3ab+b^2+2(a+b)mu+mu^2) / [(a+b)^3 (a+mu)^2 (b+mu)^2].

The earlier memo's "(pi/4)" belongs to its own expression, which is I2/(64 pi) ("up to 4 pi factors");
the constant for the d^3r1 d^3r2 integral is 16 pi^2.

Mechanism: the average of e^{-mu r12}/r12 over the direction of r1 is G_mu(r1,r2)/(r1 r2), with
G_mu = (e^{-mu|r1-r2|} - e^{-mu(r1+r2)})/(2 mu) the half-line Dirichlet Green's function (quadrature
check 3e-22). Each leaf then contributes r^2/r = r, i.e. one derivative.

- Symbolic: half-line integral of e^{-a r1 - b r2} G_mu equals psi2 identically (termwise, sympy).
- Numeric (a,b,mu) = (1.3, 0.7, 0.45): quadrature 67.3787737670912735305 vs formula
  67.3787737670912735306, relative difference 1.7e-21.
- Coulomb limit: 32 pi^2 (a^2+3ab+b^2)/[a^2 b^2 (a+b)^3]; quadrature agrees to 2.5e-21; at a = b = 2Z
  it is the textbook 5 pi^2/(8 Z^5).
- geovac: `compute_rk_exact(1,0,1,0,1,0,1,0,0)` = 5/8 and `compute_rk_exact(1,0,1,0,2,0,2,0,0)` = 17/81,
  both reproduced exactly from derivatives of psi2 at mu = 0 (the 2s polynomial becomes powers of -d_b).

## 2. Three-site chain

    I3 = int d^3r1 d^3r2 d^3r3 e^{-a r1 - b r2 - c r3} e^{-mu r12}/r12 e^{-nu r23}/r23
       = 64 pi^3 d_a d_c psi3(a,b,c;mu,nu),
    psi3 = (a+2b+c+mu+nu) / [(a+b+c)(a+mu)(c+nu)(b+mu+nu)(a+b+nu)(b+c+mu)].

Derivative structure, derived: a vertex of degree d with density r^k e^{-x r} carries radial weight
r^{k+2-d}. Leaves (d = 1) get one derivative; the middle vertex (d = 2) gets none.

- Symbolic: the half-line integral with two Green's functions equals psi3 identically (sympy, termwise).
- Numeric (1.3, 0.7, 2.1; 0.45, 0.9): quadrature 22.1143153979156071534 vs formula
  22.1143153979156071535, relative difference 3.4e-21. Wrong-answer controls: adding d_b misses by 25%;
  dropping the total-energy factor (a+b+c) misses by 45%.
- Coulomb limit agrees with quadrature to 3.4e-21.
- Corpus exact values: Paper 63's one-centre chain limits (`tests/test_paper63_three_electron_k3.py`)
  23/54, 133/144, 289/225 at zeta = (1,1,1), (1,2,3), (2,1,3), central zeta first, are reproduced exactly
  by zeta-cubed-product x 64 d_a d_c psi3 at mu = nu = 0.
- Poles of I3: a+mu, c+nu, b+mu+nu, a+b+nu, b+c+mu, a+b+c. These are the vertex-exponent sums over
  connected subgraphs, plus the masses of the lines leaving each subgraph.

Source caveat: the psi3 formula is quoted from 1709.02813; the ar5iv extraction placed it in Eq. (13)
but returned five denominator factors there and six in the Introduction. The six-factor form is the one
verified here by direct integration; the equation number is not confirmed against the PDF.

## 3. Scope

Covered: one-centre trees with spherical densities and Yukawa or Coulomb lines. The integral is
(4 pi)^n times a product of (-d_x)^{k+2-d} acting on the flat-space tree coefficient.

- Rational only when k+2-d >= 0 at every vertex. A 1s vertex of degree 3 has weight 1/r, which is an
  integral over the exponent, not a derivative, and gives logarithms. Checked: the four-density Coulomb
  star at unit zeta is 0.3098044805149232..., with 192 S = -71 + 3840 ln 2 - 2304 ln 3 (integer-relation
  search, 30 digits; no rational found alone). So "one-centre s-type trees are rational" is false as
  stated in the earlier memo; chains of 1s densities are.
- Rings: no leaf, the angular average does not factor (Paper 63 sec:density_ring; dilogarithms).
- l > 0: the partial-wave kernel r_<^l/r_>^{l+1} replaces G; not a Dirichlet propagator. Not covered
  without extra structure.
- Two centres: no common radial variable. Not covered.

Corpus statements. Paper 58 (the one-centre paragraph after thm:abelian_residue) says each one-centre
integral with rational exponents is one rational times one square root; that is about Gaunt and
normalisation factors at fixed exponents, not about exponent dependence. Paper 63 states the shell-
theorem radial recursion for every s-type tree, and the rational one-centre chain values. Neither states
the pole set as a function of the exponents. The dictionary would add that one sentence; the mechanism
is already Paper 63's.

## 4. Prior art

- The angular reduction is Calais and Lowdin, J. Mol. Spectrosc. 8, 203 (1962) (r12 as integration
  variable; bibliographic record confirmed, full text not read), and for the Coulomb line Newton's
  shell theorem. Ohrn and Nordling (J. Chem. Phys., 1963) extend it to r12, r13 products (search
  record only).
- Closed forms for one-centre three- and four-electron s-type integrals with unlinked r_ij: Fromm and
  Hill 1987, Remiddi 1991, Harris 1997, King, Padhy arXiv:1701.05701 (search records; Fromm-Hill,
  Remiddi and Padhy are already cited in Paper 63).
- No source found that states the identity as "derivative of a wavefunction coefficient", or the
  subgraph pole rule, for Slater integrals. Two searches; an absence from a search.

So the rational-in-sums-of-exponents form is standard; only the dictionary is unrecorded.

## 5. Recommendation

Paper 58: nothing. Paper 63: optional single remark after the shell-theorem proof, tier
[SYMBOLIC + PRIOR ART], with a guard test if added (the three exact chain values and the star logarithm
are ready in the script):

> For 1s densities at one centre the radial recursion is the time-integral representation of the
> flat-space wavefunction coefficients of Arkani-Hamed, Benincasa and Postnikov with radius in place of
> conformal time: a Yukawa or Coulomb chain is a product of exponent derivatives of their rational
> tree coefficient, so its poles sit at sums of exponents over connected sub-chains; the angular
> reduction is that of Calais and Lowdin.

The equation number in 1709.02813 must be checked against the PDF before any citation carries one.
