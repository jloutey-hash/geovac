# Prior-art scan: cusp contact form of the water bond-angle series (2026-10-28)

Target: Paper 59, lines ~1980-2035. Five items: (1) the constant pi^2/630; (2) the slow l^-8 tail coming only from the
two H cusps meeting; (3) g as the Kato coefficient of the product density; (4) zero contact removes the l^-8 term; (5) the
Tan-contact reading in this setting. I did not re-report items already credited (Bouferguene-Jones 1998 Eq. 4.24,
LaBudde-Sahni 1960, Lowdin, Schwartz/Hill/Kutzelnigg-Morgan, Hofmann-Barth-Zwerger arXiv:1304.2891, Kimball).

Reading levels used below: **FULL** = full text read; **ABS** = publisher abstract or landing page only; **SEC** =
secondary source, or a search-engine summary of a source I did not open myself; **TITLE** = bibliographic record only.
No file was downloaded into the repo. One arXiv PDF fetch went to the harness cache, outside the repo, and I did not
parse it. No `debug/lit_scan/dl_contact_*` folders were created, so `.git/info/exclude` was not touched.

## Verdict table

| Item | Verdict | Best source | Level |
|:--|:--|:--|:--|
| 1. pi^2/630 universal constant, size (pi^2/630) g_Y g_Z \|Y-Z\|^7 | **NOT FOUND** as stated; PARTIAL for the two-centre one-centre-density case (Silverstone 1967) | Silverstone, JCP 46, 4377 (1967) | ABS |
| 2. Tail comes only from the cusp-cusp meeting | **PARTIAL** (general mechanism known: a cusp sets the k^-4 Fourier tail; the two-electron, bond-angle form was not found) | Benesch-Smith (via secondary); Lowdin (via Bouferguene-Jones) | SEC |
| 3. g as Kato's coefficient of the density | **PARTIAL** (Kato/Bingel density cusp and the molecular non-isotropic cusp are standard; their use as the size of a multi-centre-integral tail was not found) | Fournais, M. & T. Hoffmann-Ostenhof, Sorensen, AHP 8, 731 (2007) | ABS |
| 4. Zero contact removes the l^-8 term | **NOT FOUND** for nuclear cusps; analogue only (electron-electron r12 correction moves l^-4 to l^-8) | Kutzelnigg, TCA 68, 445 (1985) | SEC |
| 5. Tan-contact reading for multi-centre integrals / molecules | **NOT FOUND** (Coulomb short-distance relations are known for the electron gas and atoms, and already credited) | Hofmann-Barth-Zwerger 1304.2891 (credited) | ABS |

## Item 1: the constant

- **Silverstone, "Series Expansion for Two-Center Noninteger-n Coulomb Integrals", J. Chem. Phys. 46, 4377-4380
  (1967), doi:10.1063/1.1840555.** ABS only. No free copy was found on the AIP page, Scilit (403), author lists (exaly) or
  citing works. The AIP landing page states: "There are in general two types of terms." Its flattened superscripts give
  the families as R^(2N+lambda), with lambda running over |l1-l2| ... l1+l2, and **R^(n1+n2+3+N)**. Most coefficients come
  "from the overlap series" through the Poisson relation, and logarithms appear when n1 and n2 are noninteger but their
  sum is an integer.
  - *Inference, not reading:* if the charge distributions are r^(n-1) e^(-zeta r), noninteger n gives a large-k Fourier
    tail k^-(n+2), and the overlap family R^(n1+n2+1+N) in the companion paper (JCP 46, 4368) matches that. For integer n
    = 1 (1s-type densities) the low coefficients would vanish and the first surviving term would be R^7, consistent with
    Paper 59's (2/315)R^7. If so, Silverstone gives the R^7 coefficient for **two one-centre STO densities**, as a
    function of (n, zeta). Whether he gives it in cusp form (g_1 g_2 times a universal number) cannot be told from the
    abstract. This is the single most important unread text.
- **Silverstone and Kay, "Convergence Properties of Certain Formulas for Multicenter Electron Repulsion Integrals
  Obtained from the Bipolar Expansion of r12^-1", J. Chem. Phys. 50, 5045-5047 (1969), doi:10.1063/1.1671015.** TITLE only:
  the AIP page shows no abstract ("only available via PDF"), and its reference list cites Erdelyi, *Asymptotic
  Expansions*. **New lead.** It is the multicentre two-electron convergence analysis closest in kind to the bond-angle
  rate, and it is unread.
- The Ruedenberg two-centre closed form contains the R^7 term (Paper 59 already checks it), but nowhere in cusp form.
- Verdict: NOT FOUND for the three-centre cusp-cusp size and the universal pi^2/630. The two-centre one-centre-density
  coefficient is probably implicit in Silverstone 1967 (PARTIAL, unverified).

## Item 2: the cusp-cusp origin

- **Benesch and Smith** (the X-ray form factor / momentum density of a Kato-cusp density). SEC: the search result was a
  summary of an explicitly-correlated-Gaussian paper (arXiv:2106.10355, opened via ar5iv). It reported that "Benesch and
  Smith showed that an atomic charge density satisfying Kato's cusp conditions should correspond to a momentum density
  asymptotically decaying as n(p) ∝ |Z rho(0) p^-4 + O(p^-6)|^2". When I opened the ar5iv text I could not find that
  sentence or its reference entry in the part I read (100k of 127k characters). So the exact citation is unconfirmed;
  it is probably R. Benesch and V. H. Smith Jr., Acta Cryst. A26 (1970).
  - The same paper (ar5iv text, read): "If the nth derivative corresponds to a Dirac delta function, then the first n−2
    derivatives are continuous", and the transform "decays polynomially with 1/k^n". This is the standard
    Fourier-asymptotics fact behind item 2.
- **Lowdin**, as credited by Bouferguene-Jones (FULL, earlier PM read), ties slow multipole convergence to the cusp in the
  one-electron case.
- What is missing from all of these: two cusps on different centres, meeting through 1/r12, as the *only* non-analytic
  source in the bond angle, and the l^-8 law that follows. Verdict PARTIAL: the mechanism is textbook Fourier asymptotics
  plus the Kato cusp; its application here was not found.

## Item 3: g as the Kato coefficient

- **Fournais, M. Hoffmann-Ostenhof, T. Hoffmann-Ostenhof, Sorensen, "Non-isotropic cusp conditions and regularity of the
  electron density of molecules at the nuclei", Ann. Henri Poincaré 8, 731-748 (2007), arXiv:math-ph/0603019.** ABS. The
  abstract (paraphrased; it is too long to quote): the molecular density is written as mu times the exponential of an
  explicit function F of the nuclear charges and positions; mu is C^{1,1}; and first- and second-order non-isotropic
  cusp conditions generalising Kato are proved. This supports using a local |r - Y| coefficient of a molecular
  (product) density. It says nothing about integrals.
- Kato (1957) and Bingel (1963) for the density cusp: SEC, via the same search summaries ("His result was extended to a
  molecule by Bingel"). Verdict PARTIAL: the definition is standard; its role as the size of a multi-centre tail was not
  found.

## Item 4: zero contact

- **Kutzelnigg, Theor. Chim. Acta 68, 445 (1985)** (r12-dependent terms). SEC: a search summary of the Kutzelnigg85 PDF
  on git.irsamc.ups-tlse.fr, which I did not open myself, reports: "The energy increments corresponding to the partial
  wave expansion of the correction go asymptotically as l^-8, while conventional partial wave increments go as l^-4."
  This is the electron-electron analogue (removing the cusp steepens the partial-wave law). It is not the nuclear-cusp
  case, and the l^-4 to l^-8 step is not the same exponent shift as Paper 59's l^-8 to l^-9/l^-10. The PM should open
  the PDF and quote it before citing.
- That Gaussian fits copy a cusp only down to about 1/sqrt(alpha_max) is basis-set lore, and Paper 59 already calls it
  that. Verdict NOT FOUND for the nuclear zero-contact removal and its (1-x)^{9/2} / (1-x)^{11/2} next terms.

## Item 5: the Tan analogy

- Searches on Tan contact, Coulomb, molecule, electron-nucleus contact and momentum tail found only nuclear-physics
  generalized contacts (Weiss-Bazak-Barnea, arXiv:1503.07047; TITLE/ABS) and cold-atom work. Nothing applies contact
  relations to multi-centre integrals or to bond-angle series. Verdict NOT FOUND beyond the credited
  Hofmann-Barth-Zwerger and Kimball.

## Other items checked

- **Joy and Handler, JCP 43, S252 (1965).** ABS only: it studies inter- and intrashell convergence for the 1s sigma_g and
  2p sigma_u states of H2+ and says convergence laws were deduced. SEC: "Accurate Single-Center Expansions with
  Slater-Type Orbitals: Hydrogen Atom and Hydrogen Molecule-Ion" (Scilit record; authors and year not confirmed) reports
  the Joy-Handler rule as energy increments falling as l^-4. These are one-centre *wavefunction* energy increments, not
  integral coefficients, so they are not prior art for items 1-4.
- Ozdogan-Ruiz arXiv:1209.3755 (already cited) gives term counts, not a decay law. Kennedy-Zhao, IJQC 71, 1 (1999),
  "complementary convergence characteristics" of zeta-function and ellipsoidal expansions: ABS. It is not about a cusp
  rate.

## Top three for the PM to verify

1. **Silverstone, JCP 46, 4377 (1967)**, and its overlap companion JCP 46, 4368. Is the 1s-type R^(n1+n2+3+N) family
   first non-zero at R^7, and is its coefficient written in a form that factors into one cusp number per centre? If
   yes, item 1 becomes PARTIAL to FOUND for the two-centre case.
2. **Silverstone and Kay, JCP 50, 5045 (1969).** The convergence of bipolar-expansion multicentre ERIs is the nearest
   two-electron analogue of the bond-angle rate. It is unread, and its abstract is not even on the page.
3. **Benesch and Smith (Acta Cryst. A26, c. 1970; exact reference unconfirmed).** The cusp sets the k^-4 Fourier tail
   with coefficient Z rho(0), which is the one-centre ingredient of the pi^2/630 derivation. Cite it for the mechanism,
   alongside Kimball, if it checks out.
