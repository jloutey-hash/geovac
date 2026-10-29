# Prior-art scan: the Landau collinear-pinch size of the water bond-angle singularity

Date: 2026-10-28. Target: Paper 59, paragraph "At unequal lengths the singular point moves to x*..."
(lines ~2009-2026) and record `docs/predictions/2026-10-28_water_landau_gaussian_fresh.md`.
Method: about 15 web searches (standard + extended), 4 fetches. Two PDFs downloaded into isolated folders
(`debug/lit_scan/dl_pinch_landau_bootstrap/`, `debug/lit_scan/dl_pinch_homeier/`, both added to
`.git/info/exclude`), text extracted with `pdftotext`, nothing executed from inside them.
Reading levels: **FULL** = text read from the extracted PDF; **ABSTRACT** = abstract / front matter only;
**SECONDARY** = search-engine summary only, not a primary read. No content inferred from a title.

## Verdict table

| Item | Verdict | Best source | Level |
|:--|:--|:--|:--|
| 1. Landau / pinch reading of a molecular integral's singularities in nuclear geometry | **NOT FOUND** | none. Nearest: Kawai 1979 (position-space Landau variety), Harris 1997 (Fromm-Hill branch tracking) | SECONDARY |
| 2. Complex-angle singular point x* = (D1^2+D2^2)/(2 D1 D2) | **PARTIAL** (location standard; its rate already in the one-electron analogue) | Sack 1964; Bouferguene & Jones 1998 | ABSTRACT |
| 3. Size = simplex (segment) average of the cusp functions | **PARTIAL** (the general mechanism is standard asymptotics; the molecular use was not found) | Hannesdottir-McLeod-Schwartz-Vergu 2024, App. B; Morse-Bott stationary phase (Brumley-Templier Prop. 8.3) | FULL / SECONDARY |
| 4. The endpoint (cusp-only) value fails | **PARTIAL** (it is the simple-pinch rule used where the pinch is not isolated) | Hannesdottir et al. 2024 (simple vs non-simple pinches; Pham 1967; Berghoff-Panzer 2022) | FULL (the sentences); refs not read |

Overall: **PARTIAL.** The mechanism is standard: a non-isolated (non-simple) pinch, where the leading
coefficient is an integral over the critical set rather than a value at one point. Its use to derive the size of
a bond-angle singularity of a multicentre Coulomb integral was **not found**. Paper 59's sentence "The pinch
mechanism is standard; its use for this term was not searched" can now read "was not found in a search", with
the caveats below.

## Item 1. Landau reading of geometry singularities of molecular integrals: NOT FOUND

Searches tried: "Landau singularities" + "molecular integrals"; complex bond angle / complex internuclear
distance + analytic continuation; Slater / Coulomb integrals as "position-space Feynman integrals" + Landau;
branch points of multicentre integrals; two-electron integral as a function of bond angle + Legendre tail.
None returned a quantum-chemistry source that treats an integral's dependence on nuclear geometry by Landau or
pinch analysis.
- Kawai (1979), as summarised by the search engine: the singularity spectrum of position-space Feynman amplitudes
  "is confined to the extended Landau variety associated with the diagram in which all internal masses are set
  to zero" (SECONDARY; full reference not retrieved). This is a QFT structural statement: no molecular
  integrals, no coefficient.
- Harris, Phys. Rev. A 55, 1820 (1997), on the Fromm-Hill three-electron integral: "eliminates the necessity of
  branch tracking in the complex plane, and avoids numerical instability through a formal cancellation of the
  individual-term singularities" (SECONDARY, abstract via search). These are spurious branch singularities of a
  closed form. It is not a Landau analysis and not about geometry.
- The 2026-09-06 scan `elliptic_eri_feynman_memo.md` (same directory) had already found no chemistry source that
  reads the three-centre ERI as a Feynman-type integral. This scan does not change that.

## Item 2. The complex-angle singular point x*: PARTIAL

The location is |Y-Z|^2 = 2 D1 D2 (x* - x) = 0. That is the singular point of every Legendre expansion of a
function of the distance |r1 - r2| between two points.
- Sack, "Generalization of Laplace's Expansion to Arbitrary Powers and Functions of the Distance between Two
  Points", J. Math. Phys. 5, 245-251 (1964). Abstract (via scilit): the expansion is in Legendre polynomials of
  the angle between the two points, with radial coefficients given by Gauss hypergeometric functions (ABSTRACT).
  The Legendre coefficients of |r1-r2|^nu, and therefore of (x*-x)^{7/2}, are classical. Their asymptotics are
  the standard Darboux / Bernstein-ellipse result, with rate x* - sqrt(x*^2-1) = D</D>.
- Bouferguene & Jones, "Convergence analysis of the addition theorem of Slater orbitals and its application to
  three-center nuclear attraction integrals", J. Chem. Phys. 109, 5718 (1998), doi 10.1063/1.477194. The
  front matter seen through search says "the determination of the convergence type of such series is of great
  importance since it allows one to choose adequately the convergence accelerator" (ABSTRACT). Paper 59 already
  cites it for the geometric D</D> rate in the one-electron analogue. I did not read its Eq. 4.24 here, and I did
  not check whether it gives an amplitude as well as the rate. **PM should check this.**
- Homeier, "On Convergence Acceleration of Multipolar and Orthogonal Expansions", Internet J. Chem. 1, Art. 28
  (1998) (FULL, text extracted). It treats acceleration of Legendre/multipole series for a three-centre
  nuclear-attraction integral. It remarks on behaviour "in the vicinity of the singularities of the Legendre"
  expansion, but gives no singular-term amplitude and does not mention a cusp or Landau.
- Not found: the two-electron cusp-cusp version, i.e. that the 1/r12 meeting of two hydrogen cusps sets an
  l^-4 x (D</D>)^l tail.

## Item 3. The simplex-average size formula: PARTIAL (general mechanism standard)

- Hannesdottir, McLeod, Schwartz, Vergu, "The Landau Bootstrap", arXiv:2410.02424 (FULL, text extracted).
  - Simple pinch, main text: "In the simplest cases, solutions to the Landau equations describe a simple pinch.
    This means that the Landau equations uniquely determine the value of the internal momenta q_e and the
    parameters alpha_e that describe how the integration contour is being pinched ..., and additionally that the
    Hessian matrix associated with this singularity is negative definite."
  - Coefficient rule, App. B, after Eq. (B17): "we have replaced the alpha_i variables in the numerator by the
    value they take at the critical point, which is valid at leading order."
  - Non-simple pinches: "On the other hand, non-simple pinches--which also appear ubiquitously in Feynman
    integrals--require more care to analyze (for more details, see [9, 37])", where [9] = F. Pham, Ann. IHP
    Phys. Théor. 6, 89 (1967) and [37] = M. Berghoff and E. Panzer, arXiv:2212.06661 (references not read).
- The general asymptotic statement that fits Paper 59's Phi: for a phase or Laplace exponent whose critical set
  is a submanifold W (Morse-Bott, or "clean"), the leading coefficient is the integral over W of the amplitude,
  weighted by |det of the transverse Hessian|^(-1/2). Brumley & Templier, "Large values of cusp forms on GL(n)",
  arXiv:1411.4317, Prop. 8.3, as summarised by the search engine (SECONDARY). Standard textbook homes are
  Hörmander and Duistermaat-Guillemin (not checked here).
- My reading, not a source claim: Paper 59's critical set is the H-H segment (a one-parameter family of
  collinear pinch configurations). The weight 24 s u on the simplex plays the part of the transverse factor, and
  Phi is the critical-manifold integral. Its applying to a Slater three-centre ERI, with the cusp functions G_Y,
  G_Z as the amplitude, was not found.

## Item 4. Failure of the endpoint value: PARTIAL

The endpoint value C is the simple-pinch rule ("numerator at the critical point", Hannesdottir et al., quoted
above) applied at one point of a non-isolated critical set. The literature says that rule needs a simple
(isolated, non-degenerate) pinch. Nothing found states or measures its failure for a molecular integral. The
6.5-25% gap and the 3e-5 match are Paper 59's own.

## Caveats
- Search coverage of pre-1980 S-matrix literature (Eden et al., Pham's book, Coleman-Norton 1965) was not
  re-checked. The PM said it need not be. Potential-scattering anomalous-threshold papers found
  (arXiv:1808.08695, 1905.13198, 2212.06157) discuss branch points of partial-wave projections, but not a
  simplex-average amplitude (SECONDARY).
- Steinborn / Weniger / Sharma / Löwdin alpha-function convergence studies: no large-l amplitude result surfaced
  (SECONDARY). Bouferguene, J. Phys. A 38, 2899 (2005), on numerical Barnett-Coulson/Löwdin functions, was not
  read.

## Top 3 for the PM to verify
1. Hannesdottir et al., arXiv:2410.02424, App. B (Eq. B17) and the simple vs non-simple pinch paragraph. Then
   Pham 1967 and Berghoff-Panzer arXiv:2212.06661: check whether a critical-manifold coefficient formula for
   non-simple pinches is stated there. If it is, credit it for the mechanism.
2. Bouferguene & Jones, J. Chem. Phys. 109, 5718 (1998), full text near Eq. 4.24: check whether it gives the
   singular amplitude (not just the rate) for the one-electron analogue, and whether that amplitude is an
   endpoint or a line average.
3. Sack, J. Math. Phys. 5, 245 (1964), plus one standard Morse-Bott / clean stationary-phase reference
   (Hörmander Vol. I §7.7, or Duistermaat-Guillemin): the classical sources for x* and for the
   critical-manifold coefficient formula.
