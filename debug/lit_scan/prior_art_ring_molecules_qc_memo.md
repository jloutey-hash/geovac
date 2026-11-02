# Prior art: ring molecules (H_N rings to nanotori), quantum-chemistry side

2026-11-01. Downloads: `debug/lit_scan/dl_prior_ring_qc/` (git-excluded). Tags: VERIFIED (read in text), ABSTRACT ONLY, SECONDARY, UNVERIFIED.

## 1. Hydrogen chains and rings as correlation benchmarks

**1a. Motta et al. (Simons Collaboration), "Towards the solution of the many-electron problem in real materials: equation of state of the hydrogen chain with state-of-the-art many-body methods", Phys. Rev. X 7, 031059 (2017), DOI 10.1103/PhysRevX.7.031059, arXiv:1705.01608; erratum PRX 11, 029901 (2021).** VERIFIED (arXiv text). Open linear chains, energies per atom, R = 1.0–3.6 bohr. Table II, H10 in STO-6G, FCI per atom: −0.542439 Eh at R = 1.8 bohr (the minimum region) and −0.538436 at R = 1.6. Table III, H10 at the complete-basis-set (CBS) limit: MRCI+Q −0.56655(5) and AFQMC −0.56644(4) at R = 1.8. Table IV gives MRCI+Q per atom: cc-pVDZ −0.561486 at R = 1.8. *Bearing:* the gold-standard basis ladder, but for a chain; finite-N ring and chain energies differ.

**1b. Motta et al., "Ground-state properties of the hydrogen chain: dimerization, insulator-to-metal transition, and magnetic phases", PRX 10, 031058 (2020), DOI 10.1103/PhysRevX.10.031058, arXiv:1911.01618.** ABSTRACT ONLY (PDF downloaded, not mined). Phase map (Mott/AF, dimerization, insulator-to-metal at short spacing); choose spacing against it.

**1c. Stair & Evangelista, "Exploring Hilbert space on a budget: novel benchmark set and performance metric for testing electronic structure methods in the regime of strong correlation", arXiv:2005.11349v2 (2020).** VERIFIED (arXiv text). The journal reference (believed to be J. Chem. Phys. 2020) is UNVERIFIED. Table I gives the **H10 1D ring in STO-6G** with nearest-neighbour distance r in Å: E_FCI = −5.151378 (r = 0.75), **−5.422958 (r = 1.00)**, −5.048052 (1.50) and −4.794398 Eh (2.00), all total energies. Key finding: "DMRG generally gives the most efficient wave function representation … particularly in the 1D chain with a localized basis. For the 2D and 3D systems, all methods (except DMRG) perform best with a delocalized basis." *Bearing:* a validation point (the role STO-3G water plays in Paper 62), and direct prior art for the "two label sets" prediction: site and ring-momentum orbitals give the same FCI energy; which compresses best depends on the method.

**1d. Magoulas, Shen & Piecuch, "Addressing strong correlation by approximate coupled-pair methods with active-space and full treatments of three-body clusters", Mol. Phys. (2022), DOI 10.1080/00268976.2022.2057365, arXiv:2111.13787.** VERIFIED (arXiv text). The journal DOI is SECONDARY. It treats D6h H6 and D10h H10 rings, R_H–H = 0.6 to 2.5 Å in 0.1 Å steps, with exact FCI from GAMESS. **H6 ring, cc-pVTZ FCI** (Table I): −3.410069 Eh at 1.0 Å (the minimum), −3.394673 (1.1), −3.236119 (1.5) and −3.074787 (2.0). **H10 ring, Dunning DZ [4s/2s] FCI** (Table IV): −5.538852 at 1.0 Å (the minimum), −5.287756 (1.5) and −5.079254 (2.0). *Bearing:* strongest exact ring references beyond minimal basis; CCSD fails when stretched (−24.7 mEh at 2.0 Å, H6).

**1e. Wouters, Jiménez-Hoyos, Sun & Chan, "A practical guide to density matrix embedding theory in quantum chemistry", JCTC 12, 2706 (2016), DOI 10.1021/acs.jctc.6b00316, arXiv:1603.08443.** VERIFIED. It runs H10 and H90 rings in STO-6G with FCI as the fragment solver. Quote: "for short bond lengths larger differences are expected with respect to DMRG results due to the difference in finite size effects in a ring vs a chain."

**1f. FCI study of Peierls dimerization in H_N rings and chains, J. Chem. Phys. 138, 074315 (2013), DOI 10.1063/1.4792197.** SECONDARY; authors UNVERIFIED (the catalogue author field was blank). STO-3G FCI, N = 6–14, localization tensor. Also Fertitta, Paulus, Barcza & Legeza, PRB 90, 245129 (2014), arXiv:1406.7038, DMRG on Li/Be rings (ABSTRACT ONLY).

## 2. Slater/Sturmian bases on rings; cyclic-symmetry adaptation

**2a. Reinhardt, Toulouse, Assaraf, Umrigar & Hoggan, "Quantum Monte Carlo facing the Hartree-Fock symmetry dilemma: the case of hydrogen rings", in *Advances in Quantum Monte Carlo*, ACS Symp. Ser. 1094, ch. 5 (2012), DOI 10.1021/bk-2012-1094.ch005, arXiv:1107.3537.** VERIFIED (arXiv text). This is **direct prior art for a Slater vs Gaussian comparison on H rings.** The rings are H_{4n+2} with r = 0.74747 Å, built from one 1s STO per H (ζ = 1.17, integrals from SMILES) or a minimal 5-Gaussian contraction. Quotes:
- "in a minimal basis the orbital coefficients are fixed by the cyclic symmetry"
- "using a Slater basis set, instead of a Gaussian basis set, delays the onset of HF symmetry breaking until larger rings": H50/H62 with Slater against H46/H54 with Gaussian
- for H42 the symmetry-adapted energy per two-atom cell is "−0.950252 hartree with the Gaussian basis and −0.997003 hartree with the Slater basis"
- the symmetry-adapted (ring-momentum) determinant gives the lowest VMC/DMC energy.

*Bearing:* this is an equal-size Slater-beats-Gaussian result on rings, but only at minimal basis and Hartree–Fock/QMC level. Paper 62's DZP-level, correlated, equal-size test on a ring appears **not to have been done** (no hit found). Credit this chapter.

**2b. Elesin, Podlivaev & Openov, "On the momentum-space approach to calculation of one-electron energy spectra and wave functions of atomic clusters", arXiv:quant-ph/0403195 (2004).** VERIFIED; journal version UNVERIFIED. It extends the Fock / Shibuya–Wulfman / Koga momentum-space Sturmian method to three or more centres: "the many-center problem is basically different from the two-center one in that there appear three-center overlap integrals". It is one-electron only, on *linear chains* of He nuclei (He₄⁷⁺, He₆¹¹⁺). *Bearing:* closest Sturmian many-centre precedent; no ring or correlated case found.

**2c.** I found no Coulomb-Sturmian calculation on a ring or periodic system. One many-centre Sturmian paper ("integrals and iteration schemes", reportedly with Avery) is UNVERIFIED.

**2d. Cyclic adaptation in codes.** Standard codes adapt only to D2h subgroups; Stair ran D2h, Reinhardt a custom C_N ring code (SECONDARY). Hückel annulenes (Hückel, Z. Phys. 70, 204 (1931), SECONDARY) are the one-electron limit, E_k = α + 2β cos(2πk/N).

## 3. Carbon rings and nanotori

**3a. Kaiser, Scriven, Schulz, Gawel, Gross & Anderson, "An sp-hybridized molecular carbon allotrope, cyclo[18]carbon", Science 365, 1299 (2019), DOI 10.1126/science.aay1914, arXiv:1908.05904.** VERIFIED. Polyynic (D9h) by AFM. Quote: "sp-hybridized carbon allotropes … possess two perpendicular π-conjugated electron systems … Hoffmann predicted, in 1966, that two orthogonal ring currents would be formed in C18, causing double aromatic stabilization." *Bearing:* the out-of-plane and in-plane π systems give C18 **two ring-momentum ladders on one ring**, which is a known instance of "more than one label set".

**3b. Nanotori.** Lin & Chuu, PRB 57, 6731 (1998), persistent currents, AB period hc/e (ABSTRACT ONLY); Ceulemans, Chibotaru, Bovin & Fowler, J. Chem. Phys. 112, 4271 (2000), DOI 10.1063/1.480972, analytic tight-binding tori of any chirality/twist (SECONDARY).

**3c. Toroidal quantum numbers. Encinosa & Mott, "Curvature-induced toroidal bound states", Phys. Rev. A 68, 014102 (2003), DOI 10.1103/PhysRevA.68.014102, arXiv:quant-ph/0303083.** VERIFIED. For a particle on a torus with α = a/R, the ansatz exp[imφ] separates the around-ring angle exactly. The tube angle θ then obeys an ODE with (1 + α cos θ) coefficients and is solved as a Fourier series Σ c_n e^{inθ}. Only parity θ → −θ survives: "the solutions … can be split into odd and even parity eigenfunctions". *Bearing (load-bearing for the prediction):* **the around-ring winding m is an exact quantum number; the tube winding n is not.** Tube winding is a good label only in the thin-torus limit α → 0, and at finite α it mixes. Atomistic tori: exact labels are the finite symmetry group's irreps.

## 4. Ring of centres mapped to a single centre

I found no treatment that maps a ring of nuclei to a one-centre problem (by conformal inversion, Kelvin transform or stereographic map). SECONDARY classical facts (Moon & Spencer; Morse & Feshbach): Laplace's equation is R-separable in toroidal coordinates; Helmholtz/Schrödinger at E ≠ 0 is not. A 2025 preprint (Pereira & Schmidt, arXiv:2511.08646) claims an irregular R-separable Helmholtz solution (ABSTRACT ONLY, unvetted). The Kelvin transform preserves harmonicity only (E = 0), so it can map a ring of charges to a line of charges only for zero-energy/electrostatic problems. (Consistent with the corpus's E = 0 door.) In Fock's momentum space, a ring of N centres enters only through the phases exp(−ip·R_k). The C_N rotation is a subgroup of SO(3) ⊂ SO(4) acting on S³, so ring momentum is a natural label on the S³ side, but no single-centre reduction is known. OPEN (none found ≠ none exists).

## Design implications

1. **The prediction is partly known and needs sharpening.** Site and ring-momentum labels are related by a discrete Fourier transform, so FCI energies at fixed one-particle space are identical; ring momentum is the exact C_N label, the site label a basis choice. They differ only in (a) truncated methods (DMRG/sCI compressibility, Pauli counts after tapering) and (b) bases that span *different* spaces, for example ring-centred one-centre Sturmians against atom-centred STOs. On a torus, tube winding is approximate (Encinosa–Mott): a falsifiable sub-prediction.
2. **Paper 62's "equal size, correlated, Slater vs Gaussian" on a ring appears new.** The minimal-basis HF/QMC precedent (Reinhardt et al.) must be credited.
3. **Integral feasibility.** An H ring needs four-centre ERIs. Paper 62's Gaussian-fit evaluator handles them; the closed-form route does not.

## Recommended first system

**H6 ring, D6h, R_H–H = 1.0 Å (≈ 1.89 bohr, near the FCI minimum), 6 electrons.** Run three steps:
- **Validation:** reproduce Magoulas et al.'s cc-pVTZ FCI of **−3.410069 Eh** (84 functions; ~10⁹-CSF FCI, so check via selected-CI/CCSD(T) at sub-mEh). Also reproduce Stair & Evangelista's **H10 ring STO-6G FCI of −5.422958 Eh at 1.0 Å** (10 functions, trivial).
- **Equal-size test:** compare cc-pVDZ (5 functions per H, 30 total, 6e FCI trivial) with a Slater DZP of the same shape [2s1p], exponents set on the free H atom/H₂ per Paper 62's rule. Score both against the cc-pVTZ FCI and its own CBS extrapolation. Add a C6-adapted ring-momentum rotation (same space, so the same FCI; compare the Pauli count / 1-norm) and a ring-centred one-centre Sturmian set of equal size (a different space).
- **Then** H10 at 1.0 Å against the Magoulas DZ FCI (**−5.538852 Eh**) and the STO-6G FCI. Add a stretched point at 1.5 Å, where CCSD fails, to probe strong correlation.
