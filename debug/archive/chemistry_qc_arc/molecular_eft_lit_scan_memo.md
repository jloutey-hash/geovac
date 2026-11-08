# Literature scan: effective field theory for molecules

**Date:** 2026-10-20. Read-only scan. This file is the only output. Not a /walls or /qa pass.
**Question (PI):** is there an EFT version of the finite-nucleus matching for molecules, especially
two-centre problems in prolate spheroidal coordinates? How does the molecule affect the nuclear radius
one extracts? Where does this sit with the multi-focal composition wall (Paper 57: the label skeleton
has no native matching step)?

**Grading.** CONFIRMED = I fetched the abstract or text and it says what is attributed.
UNVERIFIABLE = not fetched, or seen only in a search-result snippet. WRONG = fetched and it does not
say it. No item below was found WRONG. Identifiers are given as fetched.

**Short verdict.** (A) and (B) are KNOWN and mature. The molecule never changes the radius. It changes
only the electronic factor that multiplies it, and the factor is a matrix element of a local operator
at each nucleus. Born-Oppenheimer has been written as an EFT, with a power counting, for H2+ itself.
(C) is the thin one. Pseudopotentials are matched to all-electron data, but no paper found here does it
with a counting rule and an error estimate per order in chemistry. (D) is empty: nothing found does EFT
or matching in prolate spheroidal coordinates. One concrete open item turned up in (A): the order
m alpha^6 finite-size term for HD+ was estimated with an LCAO wavefunction, because its constant is
state dependent.

---

## A. The nucleus inside a molecule

**What exists.**

- [KK21] Korobov, Karr, "Ro-vibrational spin-averaged transitions in the hydrogen molecular ions",
  arXiv:2107.14497, Phys. Rev. A 104, 032806 (2021), DOI 10.1103/PhysRevA.104.032806. **CONFIRMED**
  (abstract and full text via ar5iv). The nuclear terms are written as multiples of the contact
  density at one nucleus, <pi delta(r_d)>. Their Eq. (2) is a mα^5 term with (R_d/a0)^3 and a
  coefficient C_eta (a Friar-type term). Their Eq. (3) is a mα^6 term,
  proportional to (R_d/a0)^2 (Z alpha)^2 (C_theta - ln(Z R_d/a0)) <pi delta(r_d)>, with C_theta = 0.38(4).
  The paper says this term "is not a 'state-independent' term" and "would require an independent
  calculation". They **estimate it with an LCAO wavefunction**, (psi_1s(r_p) + psi_1s(r_d))/sqrt2,
  borrowing the hydrogen value. They also add the deuteron polarizability, -21.37(8) <pi delta(r_d)> kHz,
  which shifts the fundamental line by 0.33 kHz against a 0.5 kHz theory uncertainty.
- [KHKK16] Karr, Hilico, Koelemeij, Korobov, "Hydrogen molecular ions for improved determination of
  fundamental constants", arXiv:1605.05456, Phys. Rev. A 94, 050501(R) (2016),
  DOI 10.1103/PhysRevA.94.050501. **CONFIRMED** (abstract). H2+ and HD+ are "sensitive to the radii of
  the proton and deuteron" at the level of the proton-size puzzle. They need about 1e-12 accuracy.
- [KK23] Karr, Koelemeij, arXiv:2212.11683 (v3 2023). **CONFIRMED** (abstract). HD+ frequencies are
  fitted as functions of masses, Rydberg constant and proton and deuteron radii; the gain is in the
  electron mass and m_p/m_e, not the radii. Journal and CODATA 2022 use: **UNVERIFIABLE** (snippets).
- [KHK17] Korobov, Hilico, Karr, arXiv:1703.07972, PRL 118, 233001 (2017). **CONFIRMED** (abstract):
  NRQED expansion in alpha for H2+, D2+, HD+.
- [BKTV18] Brambilla, Krein, Tarrus Castella, Vairo, "The Born-Oppenheimer approximation in an
  effective field theory language", arXiv:1707.09647, Phys. Rev. D 97, 016016 (2018),
  DOI 10.1103/PhysRevD.97.016016. **CONFIRMED** (abstract and ar5iv text). The matching for H2+ is
  done to O(mα^5). Nuclei are point charges. The only contact term is electron-nucleus, with
  coefficients c_D and d_2, entering at mα^4. The ultrasoft (Lamb-shift) term at mα^5 depends on the
  electron density at the nuclear positions, as a function of internuclear distance. No nuclear-size
  term is included. In NRQED a nuclear radius would enter as a shift of that same contact coefficient
  (standard; Hill-Paz, already cited in Paper 57).
- [H19] Holsch et al., "Benchmarking theory with an improved measurement of the ionization and
  dissociation energies of H2", arXiv:1902.09471, PRL 122, 103002 (2019). **CONFIRMED** (abstract).
  D0 of ortho-H2 has an uncertainty "three times smaller than" the ~1 MHz finite-size contribution
  of the proton.
- [KPCLP19] Komasa, Puchalski, Czachorowski, Lach, Pachucki, arXiv:1909.13826, Phys. Rev. A 100, 032519
  (2019). **CONFIRMED** (abstract). The finite-size effect for H2 is included "within the
  Born-Oppenheimer approximation only".
- Heavy molecules. [ASS16] Almoukhalalati, Shee, Saue, "Nuclear size effects in vibrational spectra",
  PCCP 18, 15406 (2016), DOI 10.1039/c6cp01913g. Bibliography **CONFIRMED** (Crossref). Content
  **UNVERIFIABLE** (RSC 403; snippet): field shifts via the effective density and its R-derivatives;
  contact density in its place costs about 10%. The rotational field shift goes back to Tiemann et al.
  (Chem. Phys. 67, 133; 68, 13; 68, 21 (1982)), as cited in a Saue abstract I fetched.
- [U21] Udrescu et al., "Isotope Shifts of Radium Monofluoride Molecules", arXiv:2105.10549, PRL 127,
  033001 (2021). **CONFIRMED** (abstract): high nuclear-size sensitivity, framed as a test of the
  electron density inside the nucleus.
- [AN24] Athanasakis-Kaklamanakis, Neyens, arXiv:2402.12157. **CONFIRMED** (abstract): finite size
  shifts the rotational constant at the Hz level; proposes model-independent matter radii, HD+ as
  benchmark.

**What it settles.** The EFT logic carries over to molecules unchanged. The nuclear coefficients
(radius, Friar moment, polarizability) are local and molecule-independent. The molecule enters only
through matrix elements of local operators at each nucleus: <delta(r_a)> at leading order. They are
taken at fixed R, then averaged over vibration. So the molecule does not change the radius one
extracts. It changes the coefficient and adds correlated theory errors (hence KK23's global fit).
Heavy-atom field shifts are R-dependent and break Born-Oppenheimer mass scaling (Tiemann).

**What looks open.** (1) The state-dependent constant of the mα^6 finite-size term in H2+/HD+.
KK21 estimated it with LCAO. I found no later exact two-centre value. Zhong, Yang, Korobov, Li, Shi,
arXiv:2603.07858, Phys. Rev. A 114, 032814 (2026), on mα^6 for HD+ beyond Born-Oppenheimer, does not
mention finite size in its abstract (**CONFIRMED** abstract). (2) Finite-size effects beyond
Born-Oppenheimer in H2 are not included (KPCLP19). Both are tiny.

## B. Nuclear motion: Born-Oppenheimer as an EFT

**What exists.**

- [BKTV18] (above). Power counting from the ar5iv text: m/M ~ alpha^(3/2). Vibrational energy
  ~ mα^2 sqrt(m/M). Radial nuclear gradient ~ mα (M/m)^(1/4). Nonadiabatic coupling is suppressed by
  (m/M)^(1/4). The internuclear distance stays a dynamical variable, expanded harmonically about
  r0 ~ 1/(mα) only to estimate scales. Spherical coordinates are used for the relative nuclear motion.
  **CONFIRMED.**
- [BBMV24] Berwein, Brambilla, Mohapatra, Vairo, "One Born-Oppenheimer Effective Theory to rule them
  all", arXiv:2408.04719, Phys. Rev. D 110, 094040 (2024). **CONFIRMED** (abstract): the QCD version,
  coupled equations with nonadiabatic terms.
- [PST03] Panati, Spohn, Teufel, "Space-Adiabatic Perturbation Theory", arXiv:math-ph/0201055,
  Adv. Theor. Math. Phys. 7, 145 (2003). **UNVERIFIABLE** (abstract seen only in a search summary).
  The rigorous all-orders version, with Born-Oppenheimer as an example.
- [MT19] Matyus, Teufel, arXiv:1904.00042, DOI 10.1063/1.5097899. **CONFIRMED** (abstract): nth-order
  effective nuclear Hamiltonians, with mass corrections from outlying electronic states.
- [PK08] Pachucki, Komasa, "Nonadiabatic corrections to the wave function and energy",
  arXiv:0710.2999, DOI 10.1063/1.2952517. **CONFIRMED**: "the expansion parameter being the
  electron-nucleus mass ratio to the power 3/4". This nonadiabatic perturbation theory (NAPT) is the
  H2 production route (KPCLP19).
- [PW14] Papenbrock, Weidenmueller, arXiv:1307.1181, Phys. Rev. C 89, 014334 (2014). **CONFIRMED**
  (abstract): EFT of rotational bands on vibrations; molecules named in scope.

**What it settles.** Born-Oppenheimer is an EFT, with a power counting, and it has been matched for
H2+. The fixed-R electronic energies are the static potentials (the matching coefficients) of that EFT.
GeoVac's Level 2 prolate solutions for H2+ are exactly such static energies. The (m/M)^(1/4) the PI named
is BKTV18's nonadiabatic suppression. PK08 quotes a 3/4 power for its expansion. Both are powers of
the same hierarchy.

**What looks open.** Nothing structural for diatomics. No paper found combines the BOEFT counting with
finite nuclear size.

## C. Core versus valence

**What exists.**

- [L97] Lepage, "How to Renormalize the Schrodinger Equation", arXiv:nucl-th/9706029 (lectures, 1997).
  **CONFIRMED** (ar5iv text). Section 2 is the template: a Coulomb atom plus a hidden short-range
  potential V_s. The EFT replaces V_s by smeared contact terms. The errors in the s-state binding
  energies fall "order-by-order in a^2" (Sec. 2.4, Fig. 2). The 2S level comes out to better than 1%
  with a about a quarter of the atom's radius. The text I read does not use the word "pseudopotential".
  This is one-centre and one-electron.
- [B18] Bennett et al. (Mitas group), "New generation of effective core potentials from correlated
  calculations: 2nd row elements", arXiv:1805.00607, DOI 10.1063/1.5038135. **CONFIRMED** (abstract).
  ccECPs are matched by near-isospectrality of valence excitations, refined on dimer curves. Average
  errors about 0.03 eV and 3 mA, but "the errors for Al and Si oxide molecules at short bond lengths
  are notably larger". Read as EFT, that is the breakdown when the second centre enters the core.
- [CWC26] Cai, Wang, Chen, "Kohn-Sham Hamiltonian from Effective Field Theory: Quasiparticle Band
  Narrowing from Frozen Core Dynamics", arXiv:2604.25199 (April 2026). **CONFIRMED** (abstract and
  HTML text). The small parameter is eps_F/Delta E_c <~ 0.05. Matched contributions are subtracted
  order by order with counterterms. A formally exact order-by-order expansion is "beyond the scope of
  this work". Diagrams are given to second order only. Phillips-Kleinman is only cited. Solids, not
  molecules.
- [QCM26] Quinzi, Chiarotti, Marzari, "Dynamical pseudopotentials", arXiv:2605.04673 (May 2026).
  **CONFIRMED** (abstract). Tracing out the core gives an energy-dependent embedding potential. They
  derive generalized norm conservation and reproduce all-electron scattering over wide energy ranges.
- [LNC15] Lloyd-Williams, Needs, Conduit, "Pseudopotential for the electron-electron interaction",
  arXiv:1508.02004, Phys. Rev. B 92, 075106. **CONFIRMED** (abstract). It smooths the e-e Coulomb
  short-distance part and cuts the cost about 30 times. It is applied to Li and Be atoms.
- [E14] Evangelista, arXiv:1406.0114 (driven SRG). **CONFIRMED** (abstract). SRG downfolding, not
  core-specific.
- Mueller, Flesch, Meyer, J. Chem. Phys. 80, 3297 (1984), core polarization potentials. A dipole
  polarizability of the core with one cutoff parameter. **UNVERIFIABLE** (citation seen in a reference
  list and snippet only). This is the long-distance side of a core multipole expansion: the
  polarizability is a Wilson coefficient in EFT terms.

**What it settles.** Matching a core to all-electron data is standard (isospectrality, multi-energy
scattering, polarization tails), and so are the EFT-type limits: short bonds, energy dependence.

**What looks open.** I found no chemistry paper that (i) writes the core as a series of local
operators per angular momentum with a counting rule (core radius over valence length, or valence over
core energy), (ii) matches the coefficients in the atom, and (iii) shows the molecular error falling
by the predicted power as orders are added. CWC26 comes closest, and it is for solids and stops at
second order. This is the one gap that touches GeoVac's bottleneck. Paper 17's PK barrier,
V = A e^{-Br^2}/r^2, is one fixed-shape operator. Its A and B come from the core. It has no order
parameter, so there is no way to say what the next term is or how big the error should be.

## D. EFT or matching in two-centre or prolate spheroidal coordinates

**What exists.**

- [MTSP15] Mironova, Tupitsyn, Shabaev, Plunien, arXiv:1411.2519, Chem. Phys. 449, 10 (2015).
  **CONFIRMED** (abstract). For U2^183+ the critical distance is 38.42 fm (point nuclei) and 34.72 fm
  (extended): the two-centre analogue of the Dirac door. One-centre Dirac-Sturm basis, not prolate.
- [KHK13] Korobov, Hilico, Karr, arXiv:1304.4543, Phys. Rev. A 87, 062506 (2013). **CONFIRMED**
  (abstract): relativistic Bethe logarithm for the two-centre problem.
- [MN24] Mochizuki, Nishida, arXiv:2408.06011, Phys. Rev. C 110, 064001 (2024). **CONFIRMED**
  (abstract). Zero-range plus Coulomb via a generalized Bethe-Peierls condition, applied to three
  charges in a Born-Oppenheimer picture. The closest multi-centre "contact plus Coulomb" result found.
- Brambilla et al. (B) use spherical coordinates for the nuclear motion, not prolate ones.

**What it settles.** Nothing prolate-specific. Matching is local, so it is done near each nucleus in
one-centre coordinates. The two-centre solution enters only through its contact density.

**What looks open.** Prolate-specific work, mostly because nobody needs it. My own reasoning, not a
citation: a finite nucleus breaks prolate separability. The interior
potential is a function of r_a^2 = (R/2)^2 (xi + eta)^2, which is not of the separable form
[f(xi) + g(eta)]/(xi^2 - eta^2). So the separated labels (n_xi, n_eta) are released, while m and,
for equal nuclei, g/u parity survive. That is the two-centre version of J14's "labels l, m survive,
n released".

---

## What GeoVac could test that is not already known

Almost all of A and B is KNOWN: NRQED plus Born-Oppenheimer, both mature. The atlas point is
framing. Composition across a
separated scale *is* done natively in EFT, by locality: atomic coefficients times molecular contact
densities. So the multi-focal wall in its nuclear instances (J14, HF-4) should stay read as "GeoVac
has no native matching step" (Paper 57's wording), not as "the two focal lengths cannot be combined".
A molecule adds a third length (the bond). By the same locality it composes the same way, with
corrections suppressed by powers of (nuclear radius / bond length). That is negligible for H2.

Candidate tests, in order of how likely they are to say something new:

1. **Core EFT for the composed builder (C, the real gap).** Replace the single PK barrier with a
   local series per l: a smeared delta, then its Laplacian, and so on, in the style of L97. Fix the
   coefficients on the atom alone (quantum defects and their energy slopes per l). Then predict LiH and
   BeH2 at several R. Pre-register the counting rule and the predicted error ratio per order. Then
   check whether the R_eq error falls by that power. If it does, PK gets an order parameter and an
   error bar it lacks today. If it fails at short R as ccECPs do, that pins the breakdown scale. The
   pieces are known; a clean order-by-order test in molecules was not found. Two cautions. CLAUDE.md
   §3 lists six failed PK modifications, including an l-dependent one, so the new content must be
   the counting rule and the measured error scaling, not a reshaped barrier. And the core's real job
   is orthogonality, which is energy dependent (QCM26); a local, energy-independent series may need
   energy-derivative terms per order to do it.
2. **The state-dependent mα^6 finite-size constant in H2+/HD+ (A).** KK21 used LCAO for it. GeoVac's
   prolate H2+ solver (Paper 11, machine precision) could give the exact two-centre value at fixed R,
   then average over vibration. It needs the near-nucleus regular part of the wavefunction beyond
   <delta>, the molecular analogue of the atomic psi(n), ln n terms. Open as far as this scan found.
   Impact is small (kHz scale), and the Korobov group may have it in a paper I missed.
3. **Label survival under a finite nucleus in prolate coordinates (D).** An atlas row, not new
   physics: its leading term is <delta(r_a)> times the atomic coefficient.
4. **Not worth doing:** Born-Oppenheimer EFT counting (B). BKTV18 did it.

## Search queries run (36, abbreviated)

HD+ NRQED finite size and radii (4 variants, 2016-2026); BO EFT Brambilla/Vairo (2) and BOEFT 2024;
H2 dissociation finite size; NAPT Pachucki-Komasa; (m/M)^(1/4) expansion; Matyus-Teufel;
Panati-Spohn-Teufel; Lepage renormalized Schroedinger; EFT pseudopotential / ECP matching / RG power
counting (4, three extended); ccECP; Evangelista DSRG; molecular contact density and field shift;
Saue nuclear-size papers (3); RaF isotope shifts; CODATA 2022 HD+; two-centre Bethe logarithm;
two-centre Dirac critical distance; prolate two-centre EFT contact; Coulomb plus zero-range two
centres; Rydberg core multipole EFT; Mueller-Flesch-Meyer CPP; NR EFT chemistry review (extended);
Papenbrock-Weidenmueller (2); Visscher-Dyall nuclear models.
