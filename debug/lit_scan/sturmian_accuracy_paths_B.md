# Literature scan, axis B: explicit correlation + hyperspherical methods in exponential-type settings

Date: 2026-09-13. Read-only scan. Scope B only (correlated hyperspherical harmonics;
hyperspherical adiabatic / hyperspherical elliptic; exponential-Hylleraas molecular;
R12/F12 and transcorrelated in exponential-type settings). Other agents cover the
molecular Sturmian record, STO/B-function integral technology, and two-centre
spheroidal Sturmians.

Budget used: 25 web calls (cap 25). Sources reached at source: 4 (3 with abstract
verbatim, 1 metadata-only). Everything else is tagged UNVERIFIED and its metadata
comes from search-result text, not from a publisher page.

**Coverage gaps (stated, not filled in):** APS, AIP, ScienceDirect, IOP and PubMed all
refused WebFetch (403 / cookie wall), and the NIST Hy-CI review PDF returned binary
metadata instead of text. So the Hy-CI/E-Hy-CI molecular record, the Cohen-Alavi 2019
abstract, the Ten-no STG paper, the Kievsky correlated-HH papers, the Korobov
three-body papers and the Tung-Pavanello-Adamowicz LiH curve were NOT read at source.
Their bibliographic data below is concordant across independent search results, but
every one of them is tagged UNVERIFIED and must be re-verified against the publisher
page before it is written into a paper. No content was lost to the interruption.

**Question asked of every path:** does it reach ~1.6 mHa on a MOLECULE (two or more
centres), what is the enabling move, and where is the body-count / centre-count wall?

---

## 1. Summary table

| # | Path | Best demonstrated accuracy + system | Enabling move | What GeoVac would forfeit | Scaling wall (N_e, M centres) | Verdict | Source |
|:--|:-----|:--|:--|:--|:--|:--|:--|
| B1 | Exponential explicitly correlated diatomic (James-Coolidge / Heitler-London with explicit r12 powers) | H2 ground Born-Oppenheimer potential to ~1e-15 Ha over R = 0.1-20 au; E(1.4011) = -1.1744759314002167(3) | Closed one-dimensional integral representation + recurrence relations for the **two-centre two-electron** integral over explicitly correlated exponentials, with *arbitrary* powers of r12; arbitrary-precision arithmetic | The angular closure entirely. No Gaunt/3j algebra survives: r12 is carried as an explicit factor, not expanded, so the matrix element is a 1D integral representation plus recursions, not a selection-rule-structured sum | **N_e = 2 at M = 2.** The technology *is* a two-centre two-electron object. Pachucki's own reach statement is "an arbitrary **diatomic** molecule" - no polyatomic, no third electron | **GO on accuracy, ceiling at 2e/2c** | Pachucki, "Born-Oppenheimer potential for H2", Phys. Rev. A **82**, 032509 (2010), DOI 10.1103/PhysRevA.82.032509 - VERIFIED; and Pachucki, "Correlated exponential functions in high precision calculations for diatomic molecules", Phys. Rev. A **86**, 052514 (2012), DOI 10.1103/PhysRevA.86.052514 - VERIFIED |
| B2 | Hy-CI / E-Hy-CI (one r_ij per configuration state function; exponential correlation factor r_ij exp(-w r_ij)) | Atoms: He to picohartree (E-Hy-CI, >13 digits); Li 2S = -7.4780603234519 Ha; Be isoelectronic sequence 10-20 nHa. Molecules: H2 by Hy-CI with STO and GTO (sub-mHa, not at the B1 level) | **Restrict to ONE r_ij per CSF.** That caps the required correlated integrals at three- and four-electron, which have closed forms / recursions for s-type Slater orbitals - the full many-r_ij problem is never entered | Angular closure inside the correlated block only; the uncorrelated part of the CI keeps its ordinary structure. This is the cheapest forfeit of any path here | **N_e = 4 at M = 1.** The four-electron correlated integral is the published bottleneck (Sims-Hagstrom 2015); 5+ electrons unsolved. Molecular Hy-CI beyond H2 essentially absent | **BORDERLINE** (atoms + H2 only), and the single most transferable move on this list | Sims & Hagstrom, "Mathematical and computational science issues in high precision Hy-CI variational calculations: III. Four-electron integrals", J. Phys. B **48**, 175003 (2015), DOI 10.1088/0953-4075/48/17/175003 - UNVERIFIED; part I (three-electron), J. Phys. B **37**, 1519 (2004), DOI 10.1088/0953-4075/37/7/012 - UNVERIFIED; review: Sims, Ruiz et al., "High-precision Hy-CI and E-Hy-CI studies of atomic and molecular properties", Adv. Quantum Chem. (2021), ScienceDirect S0065327621000125 - UNVERIFIED (403) |
| B3 | Exponentially correlated ATOMIC integrals at 3 and 4 electrons | Closed forms for 2-, 3- and 4-electron atomic integrals with products of s-type Slater orbitals and r_ij exp(-lambda r_ij) | Fourier representation of exp(-lambda r12)/r12, bypassing both the spherical-harmonic addition theorem and the Feynman trick | n/a (this is the integral layer, not a method) | **Single-centre only.** Four electrons is the published maximum; no two-centre generalisation is addressed anywhere I reached | STOP as a molecular path; it is the *evidence* for B1/B2's wall | Padhy, "Analytic Evaluation of some 2-, 3- and 4- Electron Atomic Integrals Containing Exponentially Correlated Functions of r_ij", Asian J. Spectrosc. Special Issue (2012) 157-162, arXiv:1609.00269 - VERIFIED |
| B4 | R12/F12 with Slater-type **geminal** | Chemical accuracy on real polyatomics, routinely, in production codes | RI/CABS: represent all three- and four-electron integrals as sums of two-electron integrals over a large auxiliary basis; correlation factor = Slater-type geminal (Ten-no 2004), itself fitted by a short Gaussian expansion | Everything. The orbital basis is Gaussian; the RI auxiliary basis is large and dense; nothing in the machinery preserves l-selection | No body-count wall - RI removes it. The wall for GeoVac is **basis type**: there is no production F12 over exponential-type *orbitals*; the STG is the correlation factor, not the basis | **STOP** - the enabling move is exactly the densification GeoVac's composition wall names | Ten-no, "Initiation of explicitly correlated Slater-type geminal theory", Chem. Phys. Lett. **398**, 56-61 (2004) - UNVERIFIED; Kutzelnigg-Klopper RI lineage, reviewed in Ten-no, WIREs Comput. Mol. Sci. **2**, 114 (2012), DOI 10.1002/wcms.68 - UNVERIFIED |
| B5 | Transcorrelated / Jastrow similarity transformation | Chemically accurate with small basis sets; ionization potentials, atomization energies, dissociation curves of first-row atoms and molecules | Similarity-transform the Hamiltonian by a Jastrow factor; solve the non-Hermitian eigenproblem projectively (FCIQMC / TC-CCSD) | Hermiticity, and the effective Hamiltonian **contains three-body interactions** by construction | Three-body operator; GeoVac has already measured that it does not collapse under Gaunt/6j (non-abelian), and that xTC keeps atomic sparsity but loses it at two centres | **STOP - already in the failed ledger** (TC Jastrow adiabatic ~46%; TC second-quantized plateau ~3.4%; three TC 2D-variational cusp attempts; cheap-cusp x3; per-pair gamma) | Cohen, Luo, Guther, Dobrautz, Tew, Alavi, "Similarity transformation of the electronic Schrodinger equation via Jastrow factorization", J. Chem. Phys. **151**, 061101 (2019), DOI 10.1063/1.5116024, arXiv:1908.02882 - UNVERIFIED at source (metadata concordant across two independent search hits + arXiv id) |
| B6 | Correlated hyperspherical harmonics (Kievsky-Viviani-Marcucci); potential harmonics (Fabre de la Ripelle) | Nuclear: A = 3 "completely satisfactory"; A = 4 satisfactory with a pair-correlation factor; A = 16 bosons; 4-alpha "rather slowly convergent". **No molecular electronic-structure result at chemical accuracy anywhere I reached** | Multiply the HH basis by a correlation factor built as a product of one-dimensional functions fixed by a *two-body Schrodinger equation*. Note: the factor multiplies the **basis** (Hermitian, variational); it does not transform the Hamiltonian | The angular closure: correlated HH matrix elements need quadrature over the correlation factor. Potential-harmonic truncation partially restores structure but is what fails at A > 3 | **Degeneracy of the HH basis at fixed order grows prohibitively with particle number** - the stated, field-acknowledged wall. Pair correlation + potential basis insufficient for A > 3 | **BORDERLINE-to-STOP**: no molecule, ever; but the *variational* correlated-basis move is NOT the same move as the TC row in GeoVac's ledger (see Section 3) | Kievsky, Rosati, Viviani et al., correlated-HH lineage: "Correlated hyperspherical-harmonic expansion for three-nucleon systems", Few-Body Systems, DOI 10.1007/BF01077669 - UNVERIFIED; "Correlated hyperspherical-harmonic calculations for three- and four-body systems", Nuovo Cimento A, DOI 10.1007/BF02731979 - UNVERIFIED |
| B7 | Hyperspherical elliptic coordinates (HSE) for the three-body Coulomb problem | Structural, not energetic: an approximate separability of the hyperspherical adiabatic eigenvalue problem, each state labelled by a pair of HSE quantum numbers generalizing both the spheroidal quantum numbers of diatomics and the Herrick-Lin numbers of two-electron atoms | An **additional approximate integral of motion**, specific to the Coulomb interaction, that generalizes *both* the Laplace-Runge-Lenz vector of the two-body Coulomb problem *and* the separation constant of the two-centre Coulomb problem | Nothing yet - it is a coordinate/symmetry statement, not a method with a demonstrated molecular energy | Three-body. No four-body or molecular-electron extension found | **BORDERLINE on accuracy; the one structurally novel object on this axis** (see Section 4) | Tolstikhin, Watanabe, Matsuzawa, "Hyperspherical elliptic coordinates and three-body Coulomb problem", Phys. Rev. Lett. **74** (18), 3573-3576 (1995), DOI 10.1103/PhysRevLett.74.3573 - VERIFIED (metadata at source via Semantic Scholar API; abstract field null, so the physics statements in this row come from secondary text and are UNVERIFIED). Follow-up: Phys. Rev. A **63**, 062705 (2001), DOI 10.1103/PhysRevA.63.062705 - UNVERIFIED |
| B8 | Fully exponential few-body (Korobov; Frolov; four-body nonadiabatic exponential) | H2+ / HD+ nonrelativistic energies to 1e-15 - 1e-24 a.u. (three-body, one electron); H2 as a four-body system (2e + 2p) to sub-nHa; positronium hydride PsH; HeH+ with two-centre correlated orbitals (2 electrons) | Basis exp(-alpha_n R - beta_n r1 - gamma_n r2) with complex nonlinear exponents + multiprecision arithmetic; all interparticle distances explicit | Every structural asset: no angular labels, no selection rules, no sparsity. The basis is a nonlinearly optimized soup | **4 bodies = 2 electrons + 2 nuclei.** Frolov's own comment: semi-exponential expansions "cannot compete" and need substantial improvement to reach 7-9 stable digits on PsH | BORDERLINE (H2 / H2+ / PsH class only) | Korobov, Coulomb three-body bound-state variational lineage (Phys. Rev. A, 2000 onward) - UNVERIFIED; Frolov comment arXiv:physics/0503116 - UNVERIFIED; Pachucki & Komasa, "Relativistic Correction from the Four-Body Nonadiabatic Exponential Wave Function", J. Chem. Theory Comput. (2024), DOI 10.1021/acs.jctc.4c00861 - UNVERIFIED |
| B9 | (Control, outside axis B, but it is the answer to "has anything reached chemical accuracy on a 4-electron diatomic from a few-body correlated basis") Explicitly correlated **Gaussians** with shifted centres | **LiH ground-state potential energy curve, absolute accuracy not exceeding 0.3 cm^-1 (~1.4 microhartree), R = 1.8-40 bohr** - a genuine 4-electron 2-centre molecule far better than chemical accuracy | Simultaneous variational optimization of ALL nonlinear parameters (exponents *and* centre shifts) using analytic energy gradients | Everything GeoVac has. Zero symmetry structure, zero sparsity, every matrix element dense; the closed forms are Gaussian ones (pi-laden) | Practical ceiling ~4-6 electrons; cost is nonlinear optimization of thousands of parameters, not integral evaluation | **STOP as a transfer** - but it is the honest answer to the accuracy question and should be the benchmark GeoVac quotes for LiH | Tung, Pavanello, Adamowicz, "Very accurate potential energy curve of the LiH molecule", J. Chem. Phys. **134**, 064117 (2011), DOI 10.1063/1.3554211 - UNVERIFIED at source (metadata + the 0.3 cm^-1 figure from concordant search text) |

---

## 2. The wall, stated once

Across every exponential-type explicitly correlated method reached in this scan, the
ceiling has the same shape, and it is an **integral** ceiling, not a convergence
ceiling:

- **One centre:** correlated integrals are solved to 4 electrons (Hy-CI's three- and
  four-electron integrals; Padhy's closed forms), with the restriction of at most one
  r_ij per configuration.
- **Two centres:** correlated integrals are solved to **2 electrons** (Pachucki's
  two-centre two-electron integral with arbitrary r12 powers). Nothing beyond.
- **Three centres:** nothing at all. The entire exponential explicitly correlated
  literature is *diatomic*; Pachucki's own forward-looking claim is "an arbitrary
  diatomic molecule".

That last line matters for GeoVac: **the three-centre wall is not a GeoVac artifact.**
Paper 59's genus-1 elliptic Bessel moment for (XY|XZ) is GeoVac's local statement of a
ceiling the exponential-basis community has also never crossed. GeoVac is not behind
the field at three centres; the field has no three-centre exponential correlated
integral either.

The field's route past that wall is B4 (RI/CABS over a large Gaussian auxiliary basis)
and B5 (Jastrow similarity transformation) - and both are moves GeoVac has already
priced as densifying, or ledgered as failed.

---

## 3. Things that qualify or contradict a corpus statement

1. **"Exact but not accurate" has a sharper form on this axis.** The only route to
   microhartree accuracy on a real two-centre molecule with exponential functions (B1)
   reaches it by *abandoning* the closed angular form, not by improving it: r12 is
   carried as an explicit factor with arbitrary powers and the matrix element becomes a
   one-dimensional integral representation plus recursions. GeoVac's Neumann/Legendre
   expansion (Paper 12) is the opposite trade - keep the angular closure, pay in
   partial-wave convergence (the L^-3 wall). The literature says the first trade is the
   one that has ever reached microhartrees at two centres.

2. **F12 is not an exponential-basis method and must not be cited as precedent for
   one.** The Slater-type object in F12 is the *geminal correlation factor* (Ten-no
   2004), and it is itself fitted by a short Gaussian expansion; the orbital basis is
   Gaussian and the RI/CABS auxiliary basis is Gaussian. I found no production F12 or
   R12 implementation over Slater-type or other exponential-type *orbitals*. If any
   GeoVac document leans on F12 as evidence that explicit correlation composes with
   exponential bases, that support does not exist.

3. **The ledger's TC-in-hyperspherical row is narrower than the CHH literature.**
   GeoVac's failed row is "TC Jastrow in the adiabatic hyperspherical solver (~46%)" -
   a *similarity transformation of the Hamiltonian*, non-Hermitian, with the known
   consequence that the multiplicative V_ee is replaced by a first-derivative operator
   and adiabatic separation breaks. Kievsky-Viviani correlated HH is a different move:
   the correlation factor multiplies the **basis**, the problem stays Hermitian and
   variational, and the factor is fixed by a two-body Schrodinger equation rather than
   optimized. That version works at A = 3 and is usable at A = 4 in nuclear physics.
   A scope distinction to record, not a contradiction - but the ledger row should not
   be read as excluding the variational correlated-basis form.

4. **The Hy-CI "one r_ij per CSF" restriction is the field's answer to precisely the
   obstruction GeoVac measured.** GeoVac found the TC three-body operator does not
   collapse under Gaunt/6j (non-abelian). Hy-CI never generates a three-body operator:
   by capping at one r_ij per configuration it caps the integral order at three- and
   four-electron and stays Hermitian. GeoVac's own R12-CI He proof of concept
   (0.80 mHa with ~6 functions) is already an instance of this ansatz at one centre.

5. **The LiH benchmark GeoVac should be quoting is 0.3 cm^-1 (ECG, shifted centres),
   not an STO-3G or cc-pVDZ number.** Under the Section 1.5 benchmarking rule ("always
   use the strongest available baseline"), the strongest published LiH PEC is the
   explicitly correlated Gaussian one. Nothing in axis B threatens it; nothing in
   GeoVac approaches it. Worth one honest sentence wherever LiH accuracy is discussed.

---

## 4. The single most promising untraveled path

**Impose the Hy-CI restriction - at most one r12 per configuration - on GeoVac's
two-centre basis, and evaluate that one r12 with Pachucki's two-centre two-electron
exponential recursions rather than with a Neumann/Legendre expansion.**

Why this and not the others:

- It is the only device on this axis that buys correlated accuracy while **capping the
  integral order**. It never produces the three-body operator that GeoVac has already
  measured does not collapse under Gaunt/6j, and it stays Hermitian, so it sidesteps
  the non-Hermitian quantum-axis problem already in the ledger (Q12 strong
  orthogonality, symmetrising the TC operator, the cheap non-Hermitian dressing - all
  ledgered).
- The hard integral is **already solved and published**. Pachucki (2012) gives a
  compact one-dimensional integral representation plus recurrence relations for the
  two-centre two-electron integral over explicitly correlated exponentials with
  *arbitrary* powers of r12, including the awkward case of integer r12 powers with a
  vanishing nonlinear parameter. GeoVac would not be inventing integral technology; it
  would be importing a known one into a basis that already lives in the right
  coordinates (Level 2, prolate spheroidal).
- GeoVac has a working one-centre instance already: R12-CI He at 0.80 mHa with ~6
  functions. The untraveled step is the second centre, not the ansatz.
- It is not excluded by any GeoVac wall I can identify. The scale lock (Paper 60) is a
  Sturmian-posing constraint and does not bind a CI-with-one-r12 ansatz; the
  Lowdin/l-selection theorem is about orthogonalizing mixed-exponent bases, not about
  an added correlated block; the three-centre genus-1 wall does not bind a diatomic;
  PK is not involved, because the correlated block is not a pseudopotential.

**What it costs, stated honestly so the PI can price it:** the r12-carrying block has
no Gaunt/l closure - it is a dense correction block bolted onto a sparse Hamiltonian.
So the deliverable is *accuracy on diatomics as validation*, not a sparsity claim; the
qubit-resource story would have to be told about the uncorrelated part, with the
correlated block as a classical reference. And it inherits the field's ceiling exactly:
**2 electrons per correlated pair at 2 centres, and nothing at 3 centres.** It would put
a chemical-accuracy H2 / HeH+ number on the board with a nameable enabling move. It
would not move LiH, and it would not move water.

**Structural companion, worth one diagnostic day and no more:** the
Tolstikhin-Watanabe-Matsuzawa hyperspherical elliptic coordinate system. Its content is
an approximate additional integral of motion for the three-body Coulomb problem that
generalizes *both* the Laplace-Runge-Lenz vector (GeoVac's Level 1 / Fock SO(4) object)
*and* the separation constant of the two-centre Coulomb problem (GeoVac's Level 2 /
prolate spheroidal object), with a single pair of quantum numbers interpolating between
Herrick-Lin two-electron-atom labels and diatomic spheroidal labels. That is a candidate
bridge between the two lowest rungs of the Section 5 hierarchy, which the corpus
currently treats as separate natural geometries chosen by where separation occurs. It
has no chemical-accuracy result attached and must not be sold as one. Caveat: the PRL
abstract itself was not retrievable, so this characterization rests on secondary text.

---

## 5. Search trail

Queries run (WebSearch unless marked FETCH):

1. Pachucki Komasa H2 Born-Oppenheimer James-Coolidge nanohartree -> arXiv 1007.0322,
   fuw.edu.pl excited-Sigma preprint, JCTC 4c00861, arXiv 1307.6065.
2. Frolov exponential variational four-body LiH HeH+ -> arXiv physics/0503116 (comment),
   arXiv 2405.08818 (PsH), arXiv 1707.04936 (two-centre molecular ions), arXiv 0901.3942.
   Dead end for LiH/HeH+ specifically: the four-body exponential record is positronium
   and hydrogen-class, not LiH.
3. Tolstikhin Watanabe Matsuzawa hyperspherical elliptic -> PRL 74, 3573; PRA 63, 062705.
4. Kutzelnigg Klopper R12 Slater-type three-electron integrals -> Ten-no CPL 2004;
   Ten-no WIREs 2012; Gaussian three-/four-electron recurrence papers. This is where the
   "STG is the factor, not the basis" finding came from.
5. Sims Hagstrom Hy-CI LiH -> only atomic hits (Li, Be sequence, He). **Dead end: no
   Hy-CI LiH exists.** Redirected to the Hy-CI review chapter.
6. Kievsky Viviani correlated HH A=4 A=6 -> Few-Body Systems BF01077669, Nuovo Cimento
   BF02731979; the A = 3,4 nucleon + A = 16 boson statement; the 4-alpha slow-convergence
   statement.
7. Transcorrelated Cohen Luo Alavi 2019 -> TC-CCSD arXiv 2109.11182, TC pseudopotentials
   JCTC 2025, Jastrow optimization.
8. ECG LiH BH Cencek Rychlewski Bubin Adamowicz -> RMP 85, 693 (Mitroy et al. review);
   JCP 134, 064117 (LiH PEC, 0.3 cm^-1); JCP 114, 3393 (non-BO LiH).
9. FETCH arXiv 1007.0322 -> VERIFIED Pachucki PRA 82, 032509, abstract + the "arbitrary
   diatomic molecule" reach statement.
10. FETCH link.aps.org PRL 74.3573 -> **403**, APS blocks WebFetch.
11. FETCH tsapps.nist.gov Hy-CI chapter PDF -> returned binary / Illustrator metadata,
    not text. **Dead end.**
12. F12 with Slater-type ORBITALS (not geminals) -> no such implementation found; every
    hit is STG-as-correlation-factor over Gaussian orbitals. This is a negative result.
13. FETCH pubmed 10058239 -> cookie wall, no content. **Dead end.**
14. E-Hy-CI H2 Sims Hagstrom four-electron bottleneck -> J. Phys. B 48, 175003 (2015)
    four-electron integrals; J. Phys. B 37, 1519 (2004) three-electron; the "one r_ij per
    CSF -> 10-20 nHa on the Be sequence" statement; Hy-CI H2 with STO and GTO.
15. FETCH pubs.aip.org JCP 151/061101 -> **403**, AIP blocks WebFetch.
16. Tung Pavanello Adamowicz LiH PEC -> JCP 134, 064117, 0.3 cm^-1, R = 1.8-40 bohr.
17. FETCH sciencedirect S0065327621000125 (Hy-CI review chapter) -> **403**.
18. Cohen et al. exact DOI -> 10.1063/1.5116024, arXiv:1908.02882, "non-Hermitian and
    contains three-body interactions".
19. Three-electron two-centre exponentially correlated integrals -> arXiv 1209.1258
    (Pachucki diatomic), HeH+ two-centre correlated orbitals, arXiv 1609.00269 (atomic
    2-, 3-, 4-electron), plus the explicit statement that generalizing beyond 2 electrons
    is "fairly difficult" because the integrands contain several interelectronic
    distances. **This is the load-bearing negative.**
20. Fabre de la Ripelle potential harmonics -> the HH-degeneracy-becomes-prohibitive
    statement; integrodifferential reformulation.
21. FETCH arXiv 1209.1258 -> VERIFIED Pachucki PRA 86, 052514, abstract verbatim.
22. FETCH arXiv 1609.00269 -> VERIFIED Padhy, atomic only, 4 electrons max.
23. Hyperspherical adiabatic for molecular electronic structure -> only nuclear-motion
    uses (H2O vibrations, He + H2 surfaces) and the He atom. **Dead end: hyperspherical
    adiabatic has no molecular electronic-structure accuracy result.**
24. Korobov exponential three-body -> exp(-alpha R - beta r1 - gamma r2), 1e-15 to 1e-24
    a.u. on H2+ / HD+.
25. FETCH Semantic Scholar API for PRL 74.3573 -> VERIFIED metadata (title, authors,
    venue, 74(18), 3573-3576, 1995, DOI); abstract field null.

**Not pursued (out of scope B, assigned to other agents):** molecular Coulomb-Sturmian
ERI technology (an Avery chapter surfaced in query 5 and was left alone), STO/B-function
integral machinery, two-centre spheroidal Sturmians.
