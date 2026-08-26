# Explorer sweep: elliptic-period-adapted correlation basis?

**Date:** 2026-08-23
**Question:** Does any published work build a *computational basis* (for two-electron
integrals or for accelerating correlation-energy convergence) from the **analytic /
modular / elliptic-period structure of the integral itself**, rather than fitting
Gaussians or generic radial functions? I.e. does anyone bridge the elliptic-Feynman-integral
side (elliptic polylogs, Gamma(2) MMVs, Bessel moments, Broadhurst/Vanhove/Zerbini) to
explicitly-correlated / Sturmian quantum chemistry (F12/R12, Aquilanti-Avery, Hylleraas)
so that the period/modular structure *chooses or accelerates* the basis?

---

## 1. VERDICT

**STOP (unmade), with a BORDERLINE-leaning scaffold.** No published method builds a
correlation/two-electron basis from the period, modular, or elliptic-curve structure of the
integral. The two communities are cleanly disjoint. The single closest scaffold -- and it is
*already GeoVac's own continuous body* -- is the **Aquilanti-Avery momentum-space /
Fock-projection Coulomb-Sturmian ERI machinery**, which lives on exactly the S3/S4
hyperspherical-harmonic geometry where GeoVac's elliptic period sits, but which stops at
special-function evaluation and never touches the period/modular content. Nobody outside
GeoVac has taken the next step.

---

## 2. CANDIDATES

### A. Aquilanti-Avery momentum-space Coulomb-Sturmian / hyperspherical-harmonic ERIs -- CLOSEST
- **Who:** James E. Avery & John S. Avery (Copenhagen); Vincenzo Aquilanti (Perugia) and coworkers.
- **Construction:** ERIs for molecular Coulomb-Sturmians via **Fock projection of momentum
  space onto the unit S3** (Fourier transforms of Coulomb-Sturmians = 4D hyperspherical
  harmonics), densities expanded in 2k-Sturmians, angular part done with hyperspherical-harmonic
  (Gaunt-type) algebra. Benchmarked ~40 ns/ERI, "competitive with Gaussians, 3-4 OoM faster
  than STO-in-Gaussians."
- **Accelerates or evaluates?** *Evaluates* (fast, exact special-function evaluation). The
  basis is chosen for completeness + integral cheapness + the SO(4) symmetry of the hydrogenic
  problem -- **not** from period structure.
- **Primary cites (resolved):**
  - Avery & Avery, "Molecular Integrals for Exponential-Type Orbitals Using Hyperspherical
    Harmonics," *Adv. Quantum Chem.* (2014), ScienceDirect S0065327614000057.
  - Avery & Avery, "Fast Electron Repulsion Integrals for Molecular Coulomb Sturmians,"
    *Adv. Quantum Chem.* (2013), ScienceDirect B9780124115446000066.
  - Aquilanti et al., "Hyperspherical harmonics as Sturmian orbitals in momentum space: a
    systematic approach to the few-body Coulomb problem" (review).
  - J. Avery, *Hyperspherical Harmonics and Generalized Sturmians* (Kluwer, 2000) /
    *Hyperspherical Harmonics and Their Physical Applications* (World Scientific, 2018).
- **Distance to bridge:** ONE structural step. The Fock-projected S3 momentum picture is
  literally the same geometry on which GeoVac identifies the ERI as a Legendre/Gamma(2) elliptic
  period. Avery's lineage supplies the basis + the momentum machinery; it simply never asks
  whether the *transcendence/period class* of the resulting integral should inform the basis.
  This is the GeoVac corpus's own observation that "the chemistry-elliptic-Feynman bridge is
  unmade" -- confirmed from the chemistry side: the scaffold exists, the period step is absent.

### B. Explicitly-correlated F12/R12 (Klopper-Tew-Kutzelnigg / Ten-No) -- accelerates, cusp-driven
- **Who:** Kutzelnigg, Klopper, Tew, Ten-No, Noga, Valeev, and the broad F12 community.
- **Construction:** Augment the wavefunction with a geminal correlation factor f(r12). The
  factor is chosen to satisfy the **Kato electron-electron cusp** (a *local* analytic
  condition at coalescence), giving L^-7 correlation-energy convergence vs L^-3 for orbitals
  (Kutzelnigg-Morgan). Ten-No's Slater geminal (1-e^{-g r12})/g is the de facto standard, in
  practice expanded in a few Gaussians for integral convenience.
- **Accelerates or evaluates?** *Accelerates convergence* -- but from the **cusp** (local
  short-range analytic behavior), NOT from the integral's global period/modular structure.
- **Primary cites:** Kutzelnigg & Morgan, *J. Chem. Phys.* 96, 4484 (1992) [L^-7 law];
  Ten-No, *Chem. Phys. Lett.* 398, 56 (2004) [Slater geminal]; Tew-Klopper-Kutzelnigg review
  "Explicitly Correlated R12/F12 Methods for Electronic Structure," *Chem. Rev.* 112 (2012).
- **Distance to bridge:** FAR in motivation, adjacent in spirit. F12 is the community that
  *does* choose basis functions from the analytic structure of the correlation problem -- but
  the structure it uses is the cusp, a local singularity, not the elliptic period of the whole
  integral. GeoVac's own F12-in-Fock-momentum probe already showed the Slater geminal does not
  lower the genus, so F12's analytic lever and the period lever are orthogonal.

### C. Feynman-integral period/modular machinery (Broadhurst, Vanhove, Weinzierl, Broedel-Duhr) -- evaluates fixed diagrams
- **Who:** Broadhurst; Bailey-Borwein-Broadhurst-Glasser; Vanhove; Zerbini; Weinzierl and
  Bogner-Mueller-Stach; Broedel-Duhr-Dulat-Penante-Tancredi.
- **Construction:** Express families of *master integrals* on elliptic curves / K3 in an
  epsilon-form (canonical differential-equation) basis of **elliptic polylogarithms / iterated
  integrals of modular forms**; evaluate via q-expansion after a modular transformation that
  makes the nome small (fast numerics). Bessel-moment closed forms via contour integration of
  lattice Green functions.
- **Accelerates or evaluates?** *Evaluates* a fixed diagram at fixed external kinematics to
  high precision. Their "basis" is a basis of **master integrals** (a finite DE system), a
  categorically different object from a variational Hilbert-space basis over continuous
  molecular geometry.
- **Primary cites (resolved):**
  - Bailey, Borwein, Broadhurst, Glasser, "Elliptic integral evaluations of Bessel moments,"
    *J. Phys. A* 41 (2008) 205203, arXiv:0801.0891.
  - Broedel, Duhr, Dulat, Penante, Tancredi, "Elliptic polylogarithms and Feynman parameter
    integrals," *JHEP* 05 (2019) 120, arXiv:1902.09971.
  - "Modular transformations of elliptic Feynman integrals," *Nucl. Phys. B* (2021),
    arXiv:2011.07311 [the small-nome fast-evaluation trick].
  - (GeoVac-given, not re-verified here) Bogner-Mueller-Stach-Weinzierl arXiv:1907.01251;
    Broadhurst-Dorigoni arXiv:2607.14020.
- **Distance to bridge:** FAR. This community has no notion of a variational basis; it prices
  and evaluates diagrams. It supplies the *value-side* hand-off (already the GeoVac Route-C
  named endpoint), not a basis-construction method.

### D. Special-function/analytic-structure-adapted correlation expansions -- adjacent, not period
- **Who:** e.g. the Fourier-Legendre expansion of the two-electron-atom density matrix
  (arXiv:0909.3992); Hooke's-atom / harmonium closed-form studies; exotic-harmonium work
  (arXiv:2504.18118).
- **Construction:** Expand correlation quantities (density matrix, pair density) in
  special-function bases (Legendre, etc.) and exploit their analytic structure for
  smoothness/convergence.
- **Accelerates or evaluates?** Mixed -- analytic-structure-aware, but the structure exploited
  is smoothness / known special functions, never the period or modular class of an integral.
- **Distance to bridge:** FAR. Illustrates that "analytic-structure-adapted basis" exists as a
  habit of mind in chemistry, but always at the level of elementary special functions, not
  motives/periods.

---

## 3. CLOSEST EXISTING BRIDGE

**The Aquilanti-Avery momentum-space Fock-projection Coulomb-Sturmian ERI method (Candidate A).**

Precise gap between it and an "elliptic-period-adapted correlation basis":
- Avery's method already *lives on the S3/S4 hyperspherical geometry* obtained by Fock
  projection -- the exact substrate on which the two-center/three-center bond ERI is a
  Legendre/Gamma(2) elliptic period. The basis (Coulomb-Sturmians / hyperspherical harmonics),
  the momentum representation, and the Gaunt/angular algebra are all in place.
- What is missing is the **use of the period/modular class as a design input**: Avery chooses
  Sturmians for completeness + SO(4) symmetry + fast special-function integrals, and treats the
  resulting radial integral as "a number to compute," not as an object whose Gamma(2)/elliptic
  structure could dictate which basis functions to keep, how to grade them, or how to
  resum/accelerate the correlation series. The modular structure is never fed back into the
  basis. That feedback loop is precisely the unmade bridge -- and precisely GeoVac's Route-C /
  Paper 59 territory.

Second-closest: **F12 (Candidate B)** is the only community that *does* choose basis content
from the correlation integral's analytic behavior -- but from the **cusp**, orthogonal to the
period lever (GeoVac's own probe: the geminal does not change the genus).

---

## 4. HONEST NEGATIVE CONTENT -- why each community stops short

- **Quantum chemists (Avery/Aquilanti, F12, Hylleraas):** Basis functions are selected so that
  the *many* two-electron integrals needed across *many* nuclear geometries in an SCF/CI/CC
  loop are individually cheap and analytically tractable -- Gaussians (product theorem), Slater
  geminals, Coulomb-Sturmians (Fock-projection recurrences). The governing analytic constraint
  they build to is the **Kato cusp** (local, at coalescence), which drives the L^-7 acceleration.
  The *global* transcendence class of a single integral (elliptic period, Gamma(2) modular
  content) is irrelevant to their cost model: they never need the closed-form *value* of one
  integral, only fast repeated numerics, so the period structure carries no computational payoff
  for them and is never computed. The genus/period is invisible to the chemistry cost function.

- **Feynman-integral number theorists (Broadhurst, Vanhove, Weinzierl, Broedel-Duhr):** They
  work with a *fixed* diagram at *fixed* external kinematics and seek its closed form or a
  fast high-precision evaluation. Their "basis" is a canonical/epsilon-form basis of **master
  integrals** -- a finite DE system, not a spanning set of a variational Hilbert space
  parameterized by a continuously varying physical geometry. There is no electron, no
  Rayleigh-Ritz, no convergence-with-basis-size notion. Their modular machinery is a *value*
  tool (evaluate the period), so it hands off a number, never a basis.

- **Structural consequence:** The two cost functions never meet. Chemistry wants cheap
  repeated integrals and cusp-driven convergence; period theory wants the exact value of one
  integral. "Use the period class to choose the basis" pays off in neither community's own
  objective -- which is exactly why it is unmade, and exactly the seam GeoVac's Route-C /
  Paper 59 program sits on.

---

## Citations flagged / not independently re-verified
- Bogner-Mueller-Stach-Weinzierl arXiv:1907.01251 and Broadhurst-Dorigoni arXiv:2607.14020
  were supplied as settled by the task and NOT re-resolved in this sweep.
- ScienceDirect landing pages (Avery ERI chapters) returned HTTP 403 to automated fetch;
  their titles/authors/method are confirmed from multiple independent search snippets
  (Copenhagen research portal, ResearchGate, ScienceDirect abstracts) but the full text was
  not opened. Treat the ~40 ns/ERI benchmark as reported-in-abstract.
