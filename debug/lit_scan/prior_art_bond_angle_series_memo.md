# Prior art: bond-angle (Legendre) series of three-centre Slater integrals (2026-11-01)

Scope: Paper 59 sec:scope, the water passage (~L1960-2140). Already cited there and not re-reported: Matuzas 2026 (both Zenodo preprints), Bouferguene-Jones 1998, Silverstone 1967 (JCP 46, 4377), Hofmann-Barth-Zwerger, Barnett-Coulson, Ozdogan-Ruiz, Bluemlein et al., Bostan et al., Lesiuk-Moszynski. Earlier memos not repeated: `angle_ode_prior_art_memo.md`, `water_contact_prior_art_memo.md`, `matuzas2026_bond_angle_prior_art_memo.md`.
Downloads: `debug/lit_scan/dl_prior_angular_series/*` (git-excluded). Read levels: **VERIFIED** = read from the text; **ABSTRACT** = abstract only (Crossref or publisher); **SECONDARY** = known only through a citing source.

## Headline

1. **Silverstone, PNAS 58, 34-36 (1967): VERIFIED, credit owed.** It is a three-centre 1s two-electron integral written as a Legendre series in the angle between the two bond vectors. Each coefficient is a finite expression in exponentials, the exponential integrals E_n (logarithms enter through their "entire part"), and elementary modified spherical Bessel functions. This is Paper 59's "each channel stays in the two-centre class", published 59 years earlier. The integral type differs, (1s_a1s_a|1s_b1s_c) rather than (XY|XZ). Silverstone-Kay III (JCP 1968) then states the same channel form for (XY|XZ)-type integrals (ABSTRACT).
2. **Weniger, math-ph/0505018 (2005): VERIFIED, credit owed** for generating p functions by differentiating a scalar-orbital integral in a nuclear coordinate (Hobson's theorem; reduced-Bessel / B functions).
3. **No source found** for the l^-10 law (p_sigma at equal lengths) or the rho^l l^-5 law (p_pi), for any bond-angle ODE, or for any proof that a three-centre value is in, or outside, the two-centre class.
4. **Silverstone 1967 (JCP 46, 4377): still unread** (paywalled, not open access). Its abstract (VERIFIED) is consistent with the R^7 term being there, but the coefficient was not seen.

## Q1. Decay laws for single-centre / Legendre expansions; higher l; the derivative trick

**Silverstone, H. J., "Analytical evaluation of three- and four-center integrals of r12^-1 with Slater-type orbitals," PNAS 58, 34-36 (1967), doi:10.1073/pnas.58.1.34. VERIFIED** (all three page scans in PMC335592 read).
- The object is I = ∫∫ 1s_a(r1)² r12⁻¹ 1s_b(r2) 1s_c(r2−R−ℛ), taken as "the simplest example of a three-(noncolinear)-center integral".
- Eq. (6): I = N Σ_l (2l+1)(−1)^l P_l(cos Θ) × R^l(...)ℛ-terms {i_l(ζ_c ℛ) I_1^(l) + k_l(ζ_c ℛ) I_2^(l)}, with cos Θ = R·ℛ/(Rℛ). Eqs. (7)-(8) are finite expressions in exp, E_1, E_2 and E_2l of linear combinations of the exponents times distances, plus finite double sums.
- Definition used: Ê_n(x) = E_n(x) + (−x)^(n−1) ln x/(n−1)!, the "entire part" of E_n.
- Advantage (2) of the formula: "The internuclear angular coordinates enter in a transparent manner." Page 1 says multicentre integrals are usually done numerically "except in the special case of linear triatomic molecules" (Barker-Eyring 1954).
- **Bearing on Paper 59.** It anticipates the structure of "Each Legendre coefficient is a finite sum of rational numbers times fixed constants ... all exponentials times 1, E_1 or a logarithm, so each channel stays in the two-center class." Each Silverstone channel is visibly in that class.
- What Silverstone does **not** have: the l-independent finite basis of 31 constants with rational coefficients, the guessing done modulo primes, and any operator or ODE.
- Credit owed in that sentence. The current citations, Barnett-Coulson and Ozdogan, are the numerical-route references.

**Silverstone 1968, and Silverstone-Kay 1968-1970, JCP "Analytical Evaluation of Multicenter Integrals of r12^-1 with STOs" I-VI.** ABSTRACT (Crossref).
- I: (1-2)-type three-centre, JCP 48, 4098, doi:10.1063/1.1669744.
- II: three-centre nuclear attraction, JCP 48, 4106, doi:10.1063/1.1669745: "an infinite sum in which the internuclear angles appear in spherical harmonics, and the internuclear distances appear in modified spherical Bessel functions and exponential-type integrals."
- III (with Kay): the (2-2)-type three-centre integral, "each electron described by a two-center product", JCP 48, 4108, doi:10.1063/1.1669746. This is the (XY|XZ) family. Its result is "an infinite sum involving spherical harmonics for the internuclear angular coordinates and exponential integral and spherical Bessel-type functions for the internuclear distances", for general n, l, m, ζ.
- **Bearing.** These are the general-l versions of the channel-class statement. III should be cited beside the PNAS note for (XY|XZ). II is the one-electron analogue of Paper 59's six-constant J11 statement.
- Not checked: whether III's per-channel radial pieces are finite or infinite expansions. Paper V, the four-centre paper, says its radial integrals are "convergent infinite expansions".

**Silverstone and Kay, "Convergence properties of certain formulas for multicenter electron repulsion integrals obtained from the bipolar expansion of r12^-1," J. Chem. Phys. 50, 5045-5047 (1969), doi:10.1063/1.1671015.** TITLE only: no abstract, not open access.
- A Bar-Ilan record gives the pages as 4942-4949. Crossref gives 5045-5047, which is the version to cite.
- It is the closest unread two-electron convergence analysis.

**LaBudde and Sahni, "Study of molecular integrals. II. Three-center one-electron and two-electron integrals," J. Chem. Phys. 33, 1022-1027 (1960), doi:10.1063/1.1731326.** ABSTRACT: "Estimates for the convergence rates of the infinite series which arise in all cases are derived." Part I is Sahni-LaBudde, JCP 33, 1015 (1960), on two-centre exchange.
- SECONDARY, via Matuzas Sec. 9.1 Table 12 (read here): printed envelopes for J11 and V, with component powers 2-6, and exchange power 6.
- Matuzas reconciles these to the integrated powers −4 and −8 through cancellation.
- Bouferguene-Jones (copy read here; see provenance note) also say LaBudde-Sahni "established the asymptotic form of many integrals, including the one under study".
- **Bearing.** These are the first published rate estimates for exactly Paper 59's 1s objects. Paper 59 credits Matuzas for the laws, and Matuzas credits LaBudde-Sahni for the envelope structure.
- The 1960 rates are component envelopes, not Paper 59's integrated rates. A one-clause credit is owed in "the decay laws l^-8 and l^-4 times the distance ratio below". PI call.

**Bouferguene and Jones, JCP 109, 5718 (1998).** Already cited. VERIFIED from a third-party-hosted copy of the publisher PDF (ornaverum.org/reference/pdf/611.pdf). Legality is unclear, so it was used only to check two sentences.
- Its large-l analysis is for **s-type** orbitals only.
- It states: "since N2 appears as a power of (l+1/2), the convergence of the series 4.4 is expected to improve with high values of this parameter". So the power of l rises with the principal quantum number of the orbital.
- That is the s-orbital, n-dependent analogue of Paper 59's steeper p_sigma law. It is not the l^-10 or l^-5 laws, and it does not treat p_pi.

**Weniger, "The spherical tensor gradient operator," Collect. Czech. Chem. Commun. 70, 1225 (2005), arXiv:math-ph/0505018. VERIFIED** (abstract and conclusion).
- "it is in principle sufficient to consider only multicenter integrals of scalar functions: Higher angular momentum states can be generated by differentiation with respect to the nuclear coordinates."
- Hobson's theorem makes reduced Bessel functions "particularly compact" under Y_l^m(∇).
- Paper 59's f(s) = (s/ζ + 1/ζ²)e^{−ζs} is a reduced Bessel function, k̂_{3/2}(ζs)/ζ², and its gradient is the 2p Slater orbital. This is the B-function construction of Weniger-Steinborn (1983; SECONDARY).
- **Credit owed** at "so every case is a derivative of the 1s integral in a bond length". What stays new: applying the trick to the singular germ to obtain the decay powers.

**Matuzas 2026 (doi:10.5281/zenodo.23448184), Sec. 9-10 re-read.** It confirms that higher angular momenta are out of its scope. The further convergence literature it cites (Sidi-Hoggan on high-order Barnett-Coulson-Löwdin asymptotics, Petersson-McKoy, Klahn-Morgan, Hill, Kutzelnigg-Morgan, Flygare) is SECONDARY, and none of it gives a p or d law as far as Matuzas's text shows.

**Joy and Handler, "Convergence in one-center expansions. H2+," J. Chem. Phys. 43, S252 (1965), doi:10.1063/1.1701501. VERIFIED** (OSTI report version, osti.gov/servlets/purl/12777511, legally open).
- The finding is about energy convergence of a one-centre **wavefunction**: "a shellwise convergence law of E′ ~ 1/l^4", with exponential convergence within each l.
- **Bearing: none.** It concerns wavefunction energy increments, not integral coefficients. Paper 59 can drop it from its list of unread sources.

Also checked, no decay law: Fernández Rico-López-Ramírez (JCP 1989/1991; ABSTRACT); Bishop 1967 and Ali-Meath 1977 (SECONDARY). No p or d law was found for any multicentre expansion, Slater or Gaussian.

## Q2. Holonomic / D-finite treatment in geometry; Landau singularities

Not found, consistent with `angle_ode_prior_art_memo.md`. New items:
- **Kawai, T., "Extended Landau variety and the singularity spectrum of position-space Feynman amplitudes," Publ. RIMS 15, 767-772 (1979), doi:10.2977/prims/1195187875.** ABSTRACT/SECONDARY. The singularity spectrum of position-space amplitudes lies in an extended Landau variety, which is "always holonomic" (Kashiwara-Kawai). This is the general framework for Paper 59's "collinear pinch of the Landau equations" in position space with fixed external points, and the natural citation for "The pinch mechanism is standard". Molecular (exponential) lines are not treated.
- **Erdoğan and Sterman, "Coordinate-space singularities of massless gauge theories," arXiv:1312.0058.** ABSTRACT: pinch singularities analysed directly in coordinate space. Massless lines only. An analogue.
- **Koutschan, Paule and Suslov** (creative telescoping for relativistic Coulomb integrals) and **Koutschan-Zeilberger** (Pekeris He): ABSTRACT. The HolonomicFunctions tool has been applied to one-centre radial integrals with recurrences in quantum numbers, never with a geometric parameter.
- **Lesiuk-Moszyński Part II:** has auxiliary ODEs in the exponent parameters for two-focus Neumann functions (SECONDARY via Matuzas). Two centres only.

Verdict: no differential equation in a bond angle, and no D-finite or creative-telescoping treatment of a multicentre integral in its geometry, was found.

## Q3. Is a three-centre value in the two-centre class?

- No proof in either direction was found.
- Every exact treatment found is an infinite channel sum whose terms are in the two-centre class. Sources: Silverstone PNAS 1967 (VERIFIED); Silverstone 1968 and Silverstone-Kay III (ABSTRACT); Bosanac-Randić, JCP 56, 337 (1972), three-centre nuclear attraction as multiple sums of A_n, B_m and E_n (ABSTRACT).
- The modern literature states only "no analytically closed form", e.g. Bağcı-Hoggan arXiv:1503.03690 (SECONDARY), without proof.
- **Hirschfelder, Eyring and Rosen, JCP 4, 121 (1936), doi:10.1063/1.1749798. ABSTRACT:** for linear symmetric H3 with equal-exponent 1s functions, "All of the difficult three-center integrals for the repulsions between pairs of electrons were calculated explicitly in ellipsoidal coordinates without the use of infinite series."
  - This is the only claim of a series-free treatment of a genuine three-centre 1s two-electron integral. It is collinear only, and "explicitly" may mean quadrature (full text not open access).
  - Later linear treatments still give series: Barker-Eyring, JCP 22, 114 and 1177 (1954), "convergent series" / "rapidly converging series" (ABSTRACT).
  - Hirschfelder 1938 (JCP 6, 795) and Hirschfelder-Weygandt (JCP 6, 806) went to the differential analyzer for bent geometries and tabulated K(c,ab) = (1/π)∫ r_c⁻¹ e^{−r_a−r_b} dτ numerically (ABSTRACT).
  - Read as evidence on the open class question, this is weak: a possible closed form at the linear endpoint, where Paper 59 says the function is analytic.
- Nothing found bears on whether the infinite channel sum produces elliptic content. The question is open in the literature as well.

## Silverstone 1967 (JCP 46, 4377): the R^7 question

- Not open access (Unpaywall). Not read.
- Abstract (VERIFIED, Crossref): two families, R^{2N+λ} and R^{n1+n2+3+N}. Most coefficients come "from the coefficients of the overlap series by exploiting ∇²(Coulomb integral) = −4π (overlap integral)".
- For 1s-type distributions (n = 1) the second family starts at R^5. That Paper 59's first odd term is R^7 requires the R^5 (and R^6) coefficients to vanish at integer n. This is consistent with a Gamma-function zero, but it is an inference, not a reading.
- Note: Silverstone's Poisson route, Coulomb from overlap, is the two-centre counterpart of Paper 59's "Poisson's equation turns the density's |r−Y| cusp term into a |Z−Y|^3 term". It is worth a pointer once read.
- **Answer: still unchecked.** Paper 59's sentence "whether it gives the two-center R^7 coefficient in this form is unchecked" stays true.

## Top build-on items

1. Silverstone PNAS 1967 eqs. (6)-(8) and Silverstone-Kay III. A ready-made exact channel formula for (aa|bc) and (ab|ac): it is an independent cross-check of Paper 59's 31-constant coefficients, and a second route for the recurrence of order 41.
2. Kawai 1979 and the holonomic extended Landau variety, as the framework for deriving x* and its exponents a priori in position space.
3. Weniger's Y_l^m(∇)/B-function machinery. It would extend the p results to d functions by Hobson's theorem without new germs (prediction: each Y_l(∇) adds a fixed power per cusp-free derivative).
