# Geometry Exploration Report: e-e-cusp / explicitly-correlated methods natural for GeoVac (momentum-space, Sturmian, qubit)

**Date:** 2026-08-23 · **Agent:** Explorer (read-only literature sweep; no computation, no production edits)

## Search Summary
Swept the transcorrelated (TC), F12/explicitly-correlated, and momentum-space/Sturmian literature (arXiv, PRResearch, JCTC, JCP, PCCP, PRA, Int. Rev. Phys. Chem.), ~30 sources examined, ~14 verified. Search axes: (1) TC Hamiltonians for quantum computing (VQE/QITE/QPE/DMRG); (2) momentum-space / Fock / plane-wave / Sturmian treatments of the e-e cusp; (3) explicitly-correlated methods on Slater/Sturmian (non-Gaussian) bases with closed-form geminal integrals. Every citation below was resolved to an arXiv ID or journal DOI; the one paper not independently re-fetched is flagged.

## Verdict: **GO**
The transcorrelated method is genuinely natural for GeoVac's setting, and specifically for its two distinctive assets: the **momentum-space/Fock representation** (where the cusp's three-body byproduct is provably tamest) and the **qubit Hamiltonian** (which TC was designed to shrink). The one obstacle — the three-body operator from the similarity transform — is *least* severe exactly in momentum space (Luo-Alavi: plane-wave momentum conservation collapses TC to two-body only), and has two independent tamings (xTC 1-RDM contraction; ISDF low-rank factorization) that are compatible with a structured/sparse ERI tensor. The honest caveat: no one has built TC on a Coulomb-Sturmian / S³-Fock basis — the two-body collapse is proven for *translation-invariant* plane waves, not yet for the SO(4)/Gaunt convolution of Sturmians. That transfer is precisely the GeoVac-native research question, not a settled win.

---

## Candidate methods (ranked by GeoVac-naturalness)

### Candidate 1 — Momentum-space (plane-wave) transcorrelated Hamiltonian: two-body collapse by momentum conservation
**Who:** Luo & Alavi (foundational, HEG); Sokolov, Dobrautz, Luo, Alavi, Tavernelli (exact TC + QITE, quantum computer).
**Construction:** Jastrow similarity transform e^{-τ} H e^{τ} performed in a **plane-wave (momentum) basis**. The generic real-space TC transform makes an irreducible three-body operator; in a plane-wave basis, momentum conservation at each vertex forces the effective Hamiltonian to contain "**up to two-body operators only**," with almost no overhead vs the bare Hamiltonian.
**Result:** CBS convergence of the HEG total energy improves O(M⁻¹) → **O(M⁻⁵ᐟ³)** (M = # orbitals). The exact-TC quantum-computing paper drives this into a compact non-Hermitian Hubbard Hamiltonian solved with Ansatz-based quantum imaginary-time evolution (QITE), reporting "orders of magnitude" accuracy gains at fixed basis.
**Citations (verified):**
- Luo & Alavi, *JCTC* **14**, 1403-1411 (2018), DOI 10.1021/acs.jctc.7b01257 — arXiv:1712.07524.
- Sokolov, Dobrautz, Luo, Alavi, Tavernelli, *Phys. Rev. Research* **5**, 023174 (2023) — arXiv:2201.03049.
**GeoVac composition:** *This is the closest fit.* GeoVac's ERIs are already momentum-space objects (Paper 59: the ERI is the k-space density product ∫ρ̃*·(4π/k²)·ρ̃), and an F12 integral is the same object with the kernel swapped (Paper 59 §f12). A momentum-space Jastrow rides GeoVac's *native* representation, and the two-body collapse is exactly the property GeoVac needs to keep its qubit Hamiltonian two-body. **Open, load-bearing caveat:** plane-wave momentum conservation is continuous translation invariance; on GeoVac's S³/Coulomb-Sturmian basis the analog is SO(4) + Gaunt/6j selection, not a simple momentum delta — so whether the three-body term collapses the same way must be *checked*, not assumed. This is the single highest-value GeoVac-native question in this report.

### Candidate 2 — Non-Hermitian TC Hamiltonian for qubit algorithms: VarQITE (NISQ) and QEVE/QPE (fault-tolerant)
**Who:** McArdle & Tew (VarQITE); Uvarov & Izmaylov (quantum eigenvalue estimation / QPE).
**Construction:** Keep the *exact* non-Hermitian TC Hamiltonian (Jastrow similarity transform, real spectrum, right-eigenvector target). Because it is non-Hermitian, standard variational VQE and standard qubitization do not apply; two routes handle it:
 - **VarQITE / Ansatz-based QITE** — imaginary-time projection onto the right eigenvector (McArdle-Tew, tested on Hubbard; the compactness of the TC ground state is the payoff).
 - **QEVE (quantum eigenvalue estimation)** — a fault-tolerant algorithm for non-Hermitian Hamiltonians with real spectra (Low-Su FOCS 2024), applied by Uvarov-Izmaylov to TC electronic Hamiltonians of the **second-row atoms Li–Ne**.
**Result:** McArdle-Tew: reduced circuit resources at fixed accuracy. Uvarov-Izmaylov: TC in a *minimal* STO-6G basis reaches a **T-gate count between** those of standard qubitization at cc-pVTZ/cc-pVQZ — i.e. the small-basis TC buys the accuracy of a much larger Gaussian basis at lower gate cost; light atoms (Li, Be) beat cc-pVQZ, heavier ones (O, F, Ne) accumulate error beyond cc-pVDZ.
**Citations (verified):**
- McArdle & Tew, "Improving the accuracy of quantum computational chemistry using the transcorrelated method," arXiv:2006.11181 (2020).
- Uvarov & Izmaylov, "Accuracy and resource advantages of quantum eigenvalue estimation with non-Hermitian transcorrelated electronic Hamiltonians," arXiv:2511.21867 (2025) — *JCTC* accepted, DOI 10.1021/acs.jctc.6c00274.
- Foundational transform: Cohen, Luo, Guther, Dobrautz, Tew, Alavi, *J. Chem. Phys.* **151**, 061101 (2019) — arXiv:1908.02882.
**GeoVac composition:** GeoVac already ships qubit Hamiltonians (OpenFermion/Qiskit/PennyLane export) and prices resources in Pauli count / LCU 1-norm. A GeoVac TC Hamiltonian would slot directly into this pipeline; VarQITE (already GeoVac's tighter VQE stack per memory `vqe_via_openfermion`) is the NISQ route, QEVE the FT route. The non-normal block-encoding 1-norm is a *new* resource axis GeoVac would need to measure — GeoVac's existing 1-norm machinery (Paper 14/60) is Hermitian-only.

### Candidate 3 — Canonical Transcorrelated F12 (CT-F12): a *Hermitian, two-body* TC Hamiltonian for VQE
**Who:** Motta, Gujarati, Rice, Kumar, Masteran, Latone, Lee, Valeev, Takeshita (IBM/Valeev), building on Yanai-Shiozaki canonical transcorrelation.
**Construction:** An approximate similarity transformation by an explicitly-correlated *unitary* two-body operator (F12 geminal correlator), truncated so the result is **Hermitian and contains no more than two-particle interactions** and is singularity-free. Sidesteps non-Hermiticity entirely — usable in ordinary q-UCCSD VQE.
**Result:** On small molecules, **qubit count ÷3 and CNOT count ÷~100** to reach cc-pVTZ quality.
**Citation (verified):** Motta et al., *Phys. Chem. Chem. Phys.* **22**, 24270-24281 (2020) — arXiv:2006.02488.
**GeoVac composition:** The most drop-in of the three — a Hermitian two-body Hamiltonian composes with everything GeoVac already does (tapering, Gaunt sparsity, QWC grouping) with no new algorithmic layer. Cost: it is an *approximation* (the three-body content is folded/dropped into effective two-body), so it inherits the accuracy ceiling Uvarov-Izmaylov saw for heavier atoms. Best first target for a GeoVac proof-of-concept precisely because it avoids the non-Hermitian machinery.

### Candidate 4 — Closed-form two-center **Slater-geminal** integrals (the Slater-basis F12 that composes with GeoVac's exact ERIs)
**Who:** Lesiuk & Moszynski (analytic two-center geminal integrals); Ten-no (Slater-type-geminal F12 origin).
**Construction:** Evaluate the nonstandard two-electron integrals over Slater geminals e^{-γr₁₂} (× products of Slater orbitals) in **closed form**: an inhomogeneous 4th-order ODE for a master integral + a new special-function family (series about regular singular points) + open-ended recursions to raise powers of r₁₂.
**Result:** Analytic formulas for the full two-center geminal-integral class, validated numerically — i.e. an *exact* Slater-basis F12 integral engine at two centers.
**Citations (verified):**
- Lesiuk & Moszynski, *Phys. Rev. A* **86**, 052513 (2012) — arXiv:1209.0839.
- Ten-no, "Initiation of explicitly correlated Slater-type geminal theory," *Chem. Phys. Lett.* **398**, 56-61 (2004) — resolves via ScienceDirect (DOI not independently re-fetched; see flag).
**GeoVac composition:** This is the Slater/Sturmian-basis counterpart to GeoVac's own closed-form two-center ERIs. GeoVac already knows the geminal has an exact *regular* Fourier transform 8πγ/(k²+γ²)² (Paper 59 §f12), and that a two-center geminal pair is genus-0 elementary. Lesiuk-Moszynski is the position-space closed form for the same objects — a direct cross-check and a ready-made **two-body TC integral supply** for a Sturmian TC Hamiltonian. It composes with GeoVac's exact ERIs because both are two-center, closed-form, and share the elliptic/genus grading.

### Candidate 5 — Three-body-term tamings compatible with a structured/sparse tensor: xTC and ISDF
**Who:** Christlmaier, Kats, Alavi (xTC); ISDF-for-TC authors (2026).
**Construction:**
 - **xTC** — contract the three-body operator with the reference one-particle density matrix into an *effective two-body* correction to the 0/1/2-body integrals, restoring formal two-body scaling with no explicit three-body storage.
 - **ISDF-for-TC** — compress the grid-evaluated TC integrals with interpolative separable density fitting into a low-rank representation (+ xTC), cutting storage/integration by orders of magnitude, for flexible multi-center correlators.
**Result:** Both remove the three-body bottleneck that would otherwise densify any TC Hamiltonian.
**Citations (verified):**
- Christlmaier, Kats, Alavi, "xTC: An efficient treatment of three-body interactions in transcorrelated methods," *J. Chem. Phys.* **159**, 014113 (2023).
- "Interpolative Separable Density-Fitting for Transcorrelated Hamiltonians," arXiv:2607.17314 (2026).
**GeoVac composition:** These are the fallback if the Sturmian three-body term does *not* collapse by selection rules (Candidate 1's caveat). xTC's 1-RDM contraction is basis-agnostic and would apply to GeoVac's integrals directly; ISDF's low-rank factorization is the natural partner for GeoVac's structured-sparse ERI tensor (analogous to the density-fitting GeoVac already benchmarks against in Paper 60).

### Candidate 6 — TC-DMRG / TC-MPS (classical validation substrate)
**Who:** Baiardi, Lesiuk, Reiher (transcorrelated matrix product operators); Baiardi & Reiher (TC-DMRG).
**Construction:** Represent the (non-Hermitian) TC Hamiltonian as a matrix product operator and solve by DMRG; a VarQITE-for-MPS variant bridges to quantum-circuit ansätze.
**Citation (verified):** Baiardi, Lesiuk, Reiher, *JCTC* (2022), DOI 10.1021/acs.jctc.2c00167 — arXiv:2202.08709.
**GeoVac composition:** Not a geometry source; a classical cross-check for a GeoVac TC Hamiltonian before circuit deployment (the role GeoVac's FCI benchmarks already play for its qubit Hamiltonians).

### Candidate 7 — Momentum-space Sturmians / Fock lineage (the basis GeoVac already lives in)
**Who:** Aquilanti, Cavalli, Coletti, Di Domenico, Grossi.
**Construction:** Hyperspherical harmonics as Sturmian orbitals in momentum space — Fock's S³ projection extended into a systematic (d+1)-dimensional HH basis for the few-body Coulomb problem, with configuration-space Kepler-Coulomb Sturmians dual to momentum-space HH.
**Citation (verified):** Aquilanti et al., *Int. Rev. Phys. Chem.* **20**(4), 673-709 (2001).
**GeoVac connection:** This *is* GeoVac's substrate (Fock S³ + Coulomb-Sturmians), so it is background rather than a new import — but it is the reference that a Sturmian TC construction would build on, and it confirms the momentum-space/Sturmian duality that makes Candidate 1's momentum-space Jastrow conceivable on GeoVac's basis.

---

## Closest existing bridge + the precise gap
**Closest bridge:** Luo-Alavi momentum-space TC (Candidate 1) + Motta CT-F12 qubit result (Candidate 3) + Lesiuk-Moszynski closed-form two-center geminal integrals (Candidate 4). Together they show: (a) the cusp folds into a *two-body* Hamiltonian in momentum space, (b) that Hamiltonian shrinks the qubit/CNOT footprint by integer/two-order factors, and (c) the two-body TC integrals exist in closed form on a Slater/Sturmian basis.

**The precise gap to a GeoVac-native e-e-cusp treatment — two concrete unknowns:**
1. **Does the two-body collapse survive on S³/Sturmians?** Plane-wave TC is two-body because momentum conservation is a continuous delta at each vertex. GeoVac's basis convolves via Gaunt/6j/SO(4) selection rules, not a momentum delta. The open question: does a momentum-diagonal (or SO(4)-adapted) Jastrow make the GeoVac TC three-body operator vanish/sparsify by the *same* selection-rule mechanism, or must it be tamed by xTC/ISDF (Candidate 5)? No one has computed this — it is the natural first GeoVac calculation.
2. **The two-body TC integral supply on the Fock basis.** GeoVac has native momentum-space F12 (kernel-swapped ERI, Paper 59 §f12) and Lesiuk-Moszynski has the position-space closed forms. What is missing is assembling these into the *transformed* one/two-body TC integrals (the ∇²·∇ Jastrow-gradient terms, not just the geminal itself) in the Coulomb-Sturmian basis. GeoVac's existing two-center ERI machinery is the right engine; the TC transform adds gradient/kinetic pieces that need the same closed-form treatment.

If both resolve favorably, the deliverable is a **two-body, small-basis, non-Gaussian TC qubit Hamiltonian** that folds the Coulomb hole into GeoVac's already-sparse structure — the missing piece for GeoVac's accuracy wall, expressed in its own representation.

---

## Honest negatives (where each family stops short for us)
- **The geminal does not lower the integral genus (already in-corpus, confirmed externally).** TC/F12 buys *compactness* (fewer basis functions ⇒ fewer qubits), not cheaper per-integral evaluation. Multicenter TC integrals inherit the same elliptic (genus-1) wall as GeoVac's three-center ERIs (Paper 59). TC helps the qubit-count axis, not the polyatomic three-center-integral axis.
- **Direct momentum-space r₁₂ is hard beyond He.** Fourier-transforming a Hylleraas/explicitly-correlated wavefunction to momentum space is computationally challenging even for helium (King lineage) — so "build r₁₂ correlation directly in momentum space" as a *wavefunction* method does not scale. TC's advantage is that it moves r₁₂ into the *Hamiltonian*, avoiding this.
- **CT-F12's Hermitian two-body form is an approximation.** It folds/drops three-body content; accuracy ceilings appear for heavier second-row atoms (error > cc-pVDZ for O/F/Ne per Uvarov-Izmaylov). Exact accuracy requires the non-Hermitian route + its extra machinery.
- **Non-Hermitian TC breaks the standard qubit toolchain.** Variational VQE and standard qubitization assume Hermiticity; the non-normal Hamiltonian needs QITE or QEVE, and its block-encoding 1-norm / spectral conditioning is a genuinely different (and un-benchmarked-for-GeoVac) resource axis. GeoVac's current 1-norm/Pauli accounting is Hermitian-only.
- **The Sturmian two-body collapse is unproven.** The most GeoVac-native claim (Candidate 1) rests on a plane-wave property that may or may not transfer to S³. Treat it as the research question, not an established result — overclaiming it would be a §1.5 rhetoric violation.
- **Jastrow optimization is an added moving part.** Every TC result depends on a chosen/optimized correlator; recent work (Modular Jastrow construction, arXiv:2512.07530) shows this is still being systematized. GeoVac would inherit that non-uniqueness.

---

## Citation flags
- **Ten-no (2004), Chem. Phys. Lett. 398, 56** — resolves in search (ScienceDirect) and is the canonical STG-F12 origin, but I did not independently re-fetch the DOI page in this sweep. Verify DOI 10.1016/j.cplett.2004.09.041 before it enters a paper.
- **Uvarov & Izmaylov (2025), arXiv:2511.21867** — arXiv abstract verified; the JCTC DOI (10.1021/acs.jctc.6c00274) is listed as accepted/in-press and should be re-checked at paper-writing time.
- All other arXiv IDs and journal DOIs above were resolved during the sweep.

## Recommended next step
Hand Candidate 1's open question to the Decomposer: **does an SO(4)/momentum-diagonal Jastrow make the GeoVac transcorrelated three-body operator collapse to two-body by Gaunt/6j selection rules (the S³ analog of Luo-Alavi plane-wave momentum conservation), or must it be reduced by xTC 1-RDM contraction?** That single determination decides whether GeoVac gets a *native* two-body TC qubit Hamiltonian (a headline result) or an xTC-reduced one (still useful, less clean). Start from CT-F12 (Candidate 3) for a Hermitian proof-of-concept, using Lesiuk-Moszynski (Candidate 4) closed forms as the two-body TC integral supply and GeoVac's FCI as the classical benchmark.
