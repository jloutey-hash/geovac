# Sprint memo — TC / F12 / explicitly-correlated prior-art scout (2026-08-23)

Read-only literature scout grounding two questions for GeoVac's electron-cusp
quantum-simulation fork (TC on the Coulomb-Sturmian basis; native geminal
integrals). All load-bearing claims below are verified against the primary
source (arXiv/DOI). No code, papers, CHANGELOG, or CLAUDE.md touched.

GeoVac state being grounded (from CLAUDE.md v5.0.2-5):
- (i) native NON-HERMITIAN TC operator is cheap (kappa_V = O(1), QEVE-affordable),
  but its accuracy with a simple one-parameter Slater geminal is
  non-variationally fragile (accurate only at a hand-picked geminal width);
- (ii) variational explicitly-correlated R12-CI is accurate (He 0.45-0.80 mHa)
  and compact (q ~ 6-8), but its geminal is a non-orthogonal correlated basis fn.

===============================================================================
Q1 — Making TC principled: how is the Jastrow/geminal width fixed so the
     NON-variational TC energy is accurate and robust? Cheap internal criterion,
     or full VMC?
===============================================================================

VERDICT: A cheap, principled internal criterion EXISTS and is now STANDARD in the
Alavi-Cohen line: **minimise the VARIANCE of the TC reference energy** (NOT the
variational energy). There is now a fully DETERMINISTIC (VMC-free, analytic-
derivative) implementation of exactly that criterion. => GeoVac's observed
fragility with a hand-picked gamma is EXPECTED for an energy/eyeballed criterion,
and is FIXABLE. An internal-criterion negative would be SURPRISING given this
literature.

Verified sources:

- Haupt, Hosseini, Lopez Rios, Dobrautz, Cohen, Alavi — "Optimizing Jastrow
  factors for the transcorrelated method." arXiv:2302.13683; J. Chem. Phys. 158,
  224105 (2023).
  What it actually says: Jastrow factors obtained by **minimising the variance of
  the TC reference energy** give better, more consistent results than those from
  minimising the variational energy. Energy-minimised Jastrows yield
  non-variational TC energies that "converge very slowly to the basis set limit";
  variance-minimised ones are stable. This is precisely GeoVac's symptom
  (energy-good only at one width) and its named fix. **Cost caveat:** in THIS
  paper the variance minimisation is done by stochastic VMC (CASINO; 44 Jastrow
  params for atoms/homonuclear dimers, 80 for CN), with e-e and e-n cusp
  conditions imposed exactly and the rest optimised. So the *criterion* (variance)
  is cheap-in-principle but the *implementation* here is still full VMC.

- Lopez Rios / Alavi et al. — "Deterministic Optimisation of Jastrow Factors."
  arXiv:2506.04895; J. Chem. Phys. 163, 084107 (2025).
  What it actually says: a DETERMINISTIC alternative to VMC that minimises the
  variance of the TC reference energy **in a standard (orbital) basis set**, with
  **analytic expressions for the derivatives of the TC Hamiltonian matrix
  elements** derived and implemented. Noise-free, reproducible; can optimise from
  scratch or refine a VMC guess; yields Slater-Jastrow variances "almost as low
  as VMC variance-opt" with energies "comparable to energy-minimisation VMC." This
  IS the cheap internal criterion — no stochastic sampling, computed from the same
  integrals a TC calculation already needs.

- Cohen et al. — "Performance of a one-parameter correlation factor for
  transcorrelation: Li-Ne total energies and IPs." arXiv:2202.10443 (builds on
  J. Chem. Phys. 154, 084107 (2021)).
  What it actually says: a SINGLE-parameter correlation factor; benefit is
  simplicity ("efficient numerical-analytical schemes can be set up"). Directly
  relevant to GeoVac's one-parameter Slater geminal — a single width is a known,
  usable regime, but the parameter is set per-system/per-basis, not universal.
  (Full functional-form / transferability detail not extractable from the abstract
  page; the PDF did not render — flagged as partially-verified.)

- (Context, F12 side) Motta et al. (below, Q2) tune the Slater-geminal length
  gamma in F12 = -e^{-gamma r12}/gamma "for a given orbital basis set" by
  minimising a CT-F12/CCSD proxy energy — i.e. gamma is basis-set-adaptive and
  tuned, NOT a universal transferable constant. Ten-no's F12 correlation factor
  (1 - e^{-gamma r12})/gamma is "the de facto standard" geminal form
  (Chem. Rev. 112, 75 / cr200204r).

Q1 bottom line: The fragility is not a GeoVac defect — it is the generic
energy-criterion pathology the Alavi line documented. The principled cheap fix is
variance-of-TC-reference-energy minimisation, available deterministically
(arXiv:2506.04895) using analytic derivatives of the TC matrix elements GeoVac
already forms. Recommended GeoVac test: set the Slater-geminal gamma by minimising
the TC reference-energy variance (deterministic) rather than by scanning energy;
if that ALSO proves fragile on the Coulomb-Sturmian basis, THAT would be the
surprising, report-worthy result.

===============================================================================
Q2 — Encoding explicitly-correlated methods on quantum computers. Is
     "geminal-augmented (non-orthogonal) R12-CI as a NOQE / quantum eigenvalue
     problem" done, or a gap?
===============================================================================

VERDICT: The three ingredients exist SEPARATELY on QC, but the specific GeoVac
combination — a variational, real-space r12/geminal correlated basis encoded as a
NOQE-style generalised eigenvalue problem, on a Coulomb-Sturmian basis with native
(closed-form) geminal integrals — is NOT pre-empted. Closest single prior arts:
Uvarov-Izmaylov (non-Hermitian TC + cheap eigenvalue estimation) and NOQE
(non-orthogonal generalised eigenproblem machinery) — neither combines explicit
r12 correlation with the non-orthogonal-CI encoding.

The three separate ingredients (all verified):

(a) NON-HERMITIAN TC on QC + cheap eigenvalue estimation — the direct GeoVac
    analogue, but Gaussian:
    - McArdle & Tew, arXiv:2006.11181 (2020): first TC-on-QC; non-Hermitian TC,
      right eigenvector obtained non-variationally (imaginary-time); reduces qubit
      count. Foundational.
    - Sokolov, Dobrautz, Luo, Alavi, Tavernelli — arXiv:2201.03049; Phys. Rev.
      Research 5, 023174 (2023): "exact" (non-Hermitian) TC on QC; handles the
      non-Hermiticity (which breaks the variational principle) via **ansatz-based
      quantum imaginary-time evolution (QITE)**. Orders-of-magnitude overhead
      reduction.
    - Dobrautz et al. (Chalmers) — arXiv:2303.02007; JCTC 20, 4146 (2024):
      "Ab initio TC ... near-term hardware" — chemical accuracy on real hardware
      with smaller basis / fewer qubits.
    - **Uvarov & Izmaylov, arXiv:2511.21867 (2025)** — CLOSEST prior art to
      GeoVac's non-Hermitian-TC + QEVE fork: applies a quantum eigenvalue
      estimation algorithm (QEVE) valid for non-Hermitian, real-spectrum
      Hamiltonians; uses the **xTC** approximation (same family as GeoVac's
      v5.0.3-5 xTC). Result: TC in minimal STO-6G more accurate than standard
      cc-pVQZ for Li/Be (worse than cc-pVDZ for O/F/Ne); T-gate count between
      cc-pVTZ and cc-pVQZ qubitisation; **2.5x fewer qubits**. Basis is Gaussian
      (STO-6G / cc-pVnZ). **No Sturmian / Coulomb-Sturmian anywhere** — verified.

(b) Hermitianised F12 on QC — the "F12-on-QC" prior art, but a fixed similarity
    transform, not a variational correlated basis:
    - Motta, Gujarati, Rice, Kumar, Masteran, Latone, Lee, Valeev, Takeshita —
      "Quantum simulation of electronic structure with a transcorrelated
      Hamiltonian," arXiv:2006.02488 (2020), PRX Quantum. Uses **canonical
      transcorrelated F12 (CT-F12)** to build a **Hermitian**, singularity-free,
      2-body TC Hamiltonian solved with **q-UCCSD VQE**. Slater geminal
      F12 = -e^{-gamma r12}/gamma, gamma tuned per basis by a CCSD proxy.
      Gaussian bases only; no Sturmian. This is the canonical "explicitly-
      correlated-on-QC" reference.

(c) Non-orthogonal generalised-eigenvalue machinery on QC (NOQE) — the encoding
    GeoVac wants, but applied to determinants, NOT geminals:
    - Baek, Head-Gordon, et al. — "Say NO to Optimization: A Non-Orthogonal
      Quantum Eigensolver," arXiv:2205.09039; PRX Quantum 4, 030307 (2023).
      VERIFIED from full text: solves the generalised eigenvalue problem
      **H c = E S c**, with off-diagonal H_IJ and overlap S_IJ between
      non-orthogonal states **measured on device via a modified Hadamard test**.
      BUT the non-orthogonal states are UHF **Slater determinants** transformed by
      **UCC-type cluster operators** (multireference); correlation is multi-
      reference. Its "Jastrow" is an ORBITAL-SPACE number-operator factor
      (J_pq n_p n_q) — **NOT a real-space r12 geminal**. So the machinery matches
      GeoVac's R12-CI structure, but has never carried an explicit-correlation
      basis function.
    - Follow-on: "Towards Heisenberg Scaling: Measurement-Efficient NOQE,"
      arXiv:2606.01589 (measurement cost of the same H c = E S c encoding).

Adjacent geminal-on-QC (different kind of geminal — pairing, not e-e cusp):
    - AGP / JAGP on QC: "Correlating AGP on a quantum computer" (arXiv:2008.06138);
      "State Preparation of AGP without Number Projection" (JPCA 2023). These are
      BCS-type PAIRING geminals (number-projected BCS), NOT real-space r12 cusp
      geminals. Not the same object as an R12-CI.

Gaps confirmed (searched, nothing found):
    - F12/R12 + quantum subspace / quantum Krylov / **[2]_R12 with CABS on QC**:
      no direct combination in the literature. UNDONE.
    - Coulomb-Sturmian basis on QC: only at HF basis-construction level
      (arXiv:1811.05777; Phys. Rev. A 99, 012512 (2019)). NOT combined with
      explicit correlation, nor with qubit encoding of a correlated Hamiltonian.

Q2 bottom line — NOVEL vs pre-empted:
    - Non-Hermitian TC + cheap quantum eigenvalue estimation is PRE-EMPTED
      (Uvarov-Izmaylov 2511.21867, Sokolov 023174) — but only on Gaussian bases.
      GeoVac's differentiator here is the BASIS (Coulomb-Sturmian, native
      closed-form geminal integrals), not the algorithm.
    - "Encode a geminal-augmented (real-space r12) non-orthogonal R12-CI as a
      NOQE-style H c = E S c quantum eigenvalue problem" is a GAP. NOQE supplies
      the encoding but only for determinants; F12/TC-on-QC supplies explicit
      correlation but not the non-orthogonal-CI encoding. No paper puts an r12/
      geminal correlated basis function into the NOQE generalised eigenproblem.
    - GeoVac's angle (Coulomb-Sturmian + native geminal integrals + R12-CI-as-
      quantum-eigenvalue-problem) is therefore genuinely unclaimed.

HONEST CAVEAT for GeoVac's Q2 pivot (consistent with CLAUDE.md non-orthogonal-
encoding guardrail row): in a NOQE-style encoding the correlated-basis
NON-ORTHOGONALITY relocates the cost to device-measured many-body overlaps S_IJ
(Hadamard tests), exactly the "cost never disappears, it moves into a measured S"
pattern. GeoVac computing the geminal integrals in closed form helps the CLASSICAL
matrix assembly, but does not by itself remove the on-device S_IJ measurement
cost. Novelty is real; the encoding-cost question is open and is the thing to
measure, not assume.

===============================================================================
Refs (verified against primary source unless flagged)
===============================================================================
Q1:
- arXiv:2302.13683 / JCP 158, 224105 (2023) — variance-of-TC-ref-energy criterion (VMC impl).
- arXiv:2506.04895 / JCP 163, 084107 (2025) — DETERMINISTIC variance criterion, analytic derivs.
- arXiv:2202.10443 (+ JCP 154, 084107 2021) — one-parameter correlation factor [abstract-only, PDF didn't render].
- Chem. Rev. 112, 75 (cr200204r) — Ten-no F12 (1-e^{-gr})/g de-facto-standard geminal.

Q2:
- arXiv:2006.11181 (2020) — McArdle-Tew, first non-Hermitian TC on QC.
- arXiv:2201.03049 / PhysRevResearch 5, 023174 (2023) — Sokolov et al., exact TC via QITE.
- arXiv:2303.02007 / JCTC 20, 4146 (2024) — Dobrautz et al., ab initio TC on hardware.
- arXiv:2511.21867 (2025) — Uvarov-Izmaylov, non-Hermitian TC + QEVE, xTC, Gaussian [closest to GeoVac fork].
- arXiv:2006.02488 / PRX Quantum (2020) — Motta et al., CT-F12 (Hermitian) + q-UCCSD VQE [closest F12-on-QC].
- arXiv:2205.09039 / PRX Quantum 4, 030307 (2023) — Baek et al., NOQE (determinant multiref, H c = E S c) [verified full text].
- arXiv:2606.01589 — measurement-efficient NOQE follow-on.
- arXiv:2008.06138; JPCA 2023 — AGP/JAGP on QC (pairing geminals, not r12).
- arXiv:1811.05777 / PhysRevA 99, 012512 (2019) — Coulomb-Sturmian basis on QC, HF-level only.
