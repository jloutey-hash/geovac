# Prior art: atom-set correlating exponents and atom-to-molecule transfer (Paper 62 / SL-3d)

Date 2026-11-01. Read-only scan, about 20 searches and fetches. Downloads are in
`debug/lit_scan/dl_prior_atom_correlating/` (OSTI pages, Crossref and Semantic Scholar JSON, the Ye-Berkelbach PDF).
Already credited in Paper 62 sec:prior and not re-reported here: Dunning 1989 as the "standard premise" of
atom-optimized families, Ema 2003 as precedent for the design, Forster-Visscher 2021, Lesiuk 2015, Nemec 2010,
Rico 2004 and Koridon 2021.

## Q1. Is the SL-3d rule Dunning's rule? Yes, in substance and nearly in wording.

**F1. T. H. Dunning Jr., J. Chem. Phys. 90, 1007 (1989), doi:10.1063/1.456153.** VERIFIED from the abstract (OSTI
biblio 6593601). The full text was not reached: AIP blocked the fetch and S2 has the abstract elided. The rule, verbatim:
> "compact sets of primitive Gaussian functions effectively and efficiently describe correlation effects if the
> exponents of the functions are optimized in atomic correlated calculations, although the primitive (sp) functions
> for describing correlation effects can be taken from atomic Hartree--Fock calculations if the appropriate primitive
> set is used. Test calculations on oxygen-containing molecules indicate that these primitive basis sets describe
> molecular correlation effects as well as the ANO sets of Almloef and Taylor."

That is exactly the split SL-3d uses: s and p exponents from the free atom's HF energy, the correlating d from the free
atom's correlated energy. Dunning even worked it out on oxygen. The differences are the Slater functions, fc-CCSD(T)
on a UHF triplet in place of CISD (HF+1+2), and one d function.
**Bears on:** the sec:limits sentence "A ladder with the $3d$ set by an atomic correlated criterion would be a different
ladder and has not been run". That different ladder is the cc-pVXZ prescription applied to Slater functions, and the
record and the paper should say so. The SL-3d record's "Why this criterion" item (i) ("a correlated objective is the
only atomic one that sees it") restates Dunning's reason. Credit is owed in sec:prior and in the record.

**F2. D. E. Woon and T. H. Dunning Jr., J. Chem. Phys. 98, 1358 (1993), doi:10.1063/1.464303.** VERIFIED from the
abstract (OSTI 6942130). Verbatim:
> "The exponents for the polarization functions (dfgh) are systematically optimized for a correlated wave function
> (HF+1+2). The (sp) correlation functions are taken from the appropriate HF primitive sets; it is shown that these
> functions differ little from the optimum functions."
This is the same rule restated as a method. Its second sentence supports leaving s/p at their HF values in SL-3d.

**F3. I. Ema, J. M. García de la Vega, G. Ramírez, R. López, J. Fernández Rico, H. Meissner, J. Paldus, J. Comput.
Chem. 24, 859 (2003), doi:10.1002/jcc.10227.** VERIFIED from the abstract (Crossref). Paper 62 already cites it, but
only for "HF- and CI-optimized exponents". The abstract gives the specific Slater-type rule:
> "the exponents for the VB1 bases were obtained by minimizing the sum of the Hartree--Fock (HF) and frozen-core singles
> and doubles configuration interaction (CISD FC) energies of the respective atoms in their ground state. For H and Li,
> we minimized the sum of the HF and CISD FC energies of the corresponding diatoms (i.e., of H2 or Li2) plus the
> ground-state energy of the atom."
VB1 for B to Ne is [5,3,1], so it has a single polarization d. That makes it the closest prior art found: an STO d set
on the free atom's frozen-core correlated energy. Two differences: Ema optimizes all exponents jointly on HF plus CISD,
while SL-3d keeps s/p at HF and moves only the 3d; and their H rule uses H2, which parallels R5H (H exponents from H2,
though at RHF in R5H). **Credit owed** in sec:prior, by sharpening the ema2003 clause, and in sec:ladder at "except on
the last rung, where $1s,2s,2p$ are each set by H$_2$ RHF".

**F4. B. P. Prascher, D. E. Woon, K. A. Peterson, T. H. Dunning Jr., A. K. Wilson, Theor. Chem. Acc. 128, 69 (2011),
doi:10.1007/s00214-010-0764-0.** Citation verified via Crossref; the content is SECONDARY, from Ye-Berkelbach (F5).
Hydrogen has no atomic correlation energy, so its polarization functions are set by H2's correlation energy, and
atoms with one valence electron use the dimer. This is the cc-family source for the R5H-type rule. I did not see what
Dunning 1989 itself did for H, so that is UNVERIFIED.

**F5. H.-Z. Ye and T. C. Berkelbach, "Correlation-consistent Gaussian basis sets for solids made simple",
arXiv:2112.05824.** VERIFIED from the full text. My guess at the journal (JCTC 18, 1595 (2022)) is UNVERIFIED. It is a
clear secondary statement of the rule, verbatim:
> "The valence basis has s and p primitive orbitals ... whose exponents are determined by minimizing the HF ground
> state energy of a free atom ... The polarization functions are primitive orbitals of d angular momentum or higher
> ... whose exponents are determined by minimizing the correlation energy of a free atom."
It also reports that CCSD-optimized polarization exponents "in general agree very well with those in the all-electron
cc-pVXZ basis sets", which were CISD-optimized. The method choice (CCSD(T) vs CISD) therefore moves the exponent
little, so the SL-3d criterion is Dunning's rule up to method.

**ADF / van Lenthe-Baerends, J. Comput. Chem. 24, 1142 (2003), doi:10.1002/jcc.10255.** ABSTRACT ONLY. The exponents
are optimized for ZORA/DFT and tested on atoms and diatomic oxides. Nothing in the abstract says how the polarization
exponents were set, so it is not a precedent for a correlated atomic criterion (consistent with nemec2010's reading).

## Q2. Atom-to-molecule transfer: what is known

- **No quantitative per-function transfer ratio was found.** I found no study of the form "an exponent's atomic
  correlation gain reaches the molecule at ratio rho" that measures rho rung by rung. Dunning (F1) states transfer
  qualitatively, as "as well as the ANO sets". The search was limited, so this is not a proof of absence.
- **The fact that the water-set 3d (2.756) lands within 0.02-0.05 of the atom-set values is Dunning's empirical
  premise,** which is why atomic optimization works for cc sets. It is not new. The record already discloses it.
- **SL3d-H (water RHF gains from the tighter 3d while the atom's UHF barely moves) restates a textbook point.**
  F. Jensen, J. Chem. Phys. 115, 9113 (2001), doi:10.1063/1.1413524 (abstract VERIFIED; the body is SECONDARY via
  search snippets) built the pc-n sets because the atomic HF energy is insensitive to polarization functions. The pc-n
  polarization exponents were therefore chosen on representative molecules, while the s/p exponents were set on
  atoms. Jensen also warns that molecule-set exponents can be biased toward the molecules used. A small credit is
  owed in the R3 diagnosis and SL-3d discussion wherever "polarization the atom does not need" appears.
- **The ANO lineage builds the molecular role into the atomic objective.** J. Almlöf and P. R. Taylor, J. Chem. Phys.
  86, 4070 (1987), doi:10.1063/1.451917 (abstract VERIFIED): "atomic natural orbitals provide an excellent basis for
  molecular electronic structure calculations ... Polarization functions can be included using the same approach."
  Widmark, Malmqvist and Roos, Theor. Chim. Acta 77, 291 (1990), doi:10.1007/BF01120130, is citation-only (no
  abstract). Its density-matrix averaging over ions and, in later papers of the series, the field-polarized atom is
  SECONDARY. Roos, Lindh, Malmqvist, Veryazov and Widmark, J. Phys. Chem. A 108, 2851 (2004), doi:10.1021/jp031064+
  (abstract VERIFIED via the Lund record): the ANOs come "from the average density matrix of the ground and lowest
  excited states of the atom, the positive and negative ions, and the dimer at its equilibrium geometry."
- **Petersson CBS** (Petersson et al., J. Chem. Phys. 89, 2193 (1988), doi:10.1063/1.455064; Nyden and Petersson, J.
  Chem. Phys. 75, 1843 (1981), doi:10.1063/1.442208; both abstracts VERIFIED). These extrapolate molecular pair
  energies through the N^-1 PNO asymptotics. That is not atom-to-molecule transfer, so it does not anticipate SL3d-T.
  It is relevant only as a model of increments.
- **Bardo and Ruedenberg, J. Chem. Phys. 60, 918 (1974), doi:10.1063/1.1681168** (abstract VERIFIED). They optimize
  even-tempered primitives in molecules, water among them, and report "Transferability of the optimal even-tempered
  primitive and contracted bases between molecules is shown". That is molecule-to-molecule at the SCF level, not
  atom-to-molecule correlated transfer.

**Verdict on SL-3d:** the design is Dunning's (and Ema's) rule. The specific predictions are not anticipated by
anything found: SL3d-T (a per-rung transfer ratio in [0.844, 1.156]), SL3d-M (non-monotone R3 persists) and SL3d-C
(the correlation-part ratio below 1, offset by an RHF gain). Paper 62's existing claim of novelty, "the pre-registered
rung-by-rung measurement of the ratio", stands.

## Q3. Build-on items

1. **Clementi-Roetti / Bunge STO tables as an independent check of the O s/p exponents.** E. Clementi and C. Roetti,
   At. Data Nucl. Data Tables 14, 177 (1974), doi:10.1016/S0092-640X(74)80016-1; C. F. Bunge, J. A. Barrientos and A.
   V. Bunge, ADNDT 53, 113 (1993), doi:10.1006/adnd.1993.1003; E. Clementi and D. L. Raimondi, J. Chem. Phys. 38, 2686
   (1963), doi:10.1063/1.1733573 (abstract VERIFIED: minimal-STO exponents "optimized as to ensure the energy
   minimum"). All are citation-verified, but I did not read the tables, so no numbers are quoted. The O ³P DZ and
   minimal sets give an external check of R2-R5's atom-HF exponents. One caveat: ours minimize the
   average-of-configuration energy in a non-standard primitive pattern (1s,1s,2s / 2p,3p), so only the comparable
   rungs should be checked.
2. **Even-tempered Slater rules to set further exponents with no water input.** K. C. Raffenetti, J. Chem. Phys. 59,
   5936 (1973), doi:10.1063/1.1679962 (abstract VERIFIED): even-tempered *exponential* (Slater) bases for atomic SCF,
   with "no deleterious effect on the performance", and "Optimal bases for many atoms are determined and tabulated
   for use in molecular computation". Bardo and Ruedenberg IV, J. Chem. Phys. 59, 5966 (1973), doi:10.1063/1.1679965
   (pseudoscaling for molecules), and Feller and Ruedenberg, Theor. Chim. Acta 52, 231 (1979),
   doi:10.1007/BF00547681 (citation only), extend these. A two-parameter (alpha, beta) family per l would fix whole
   rungs from the atom. Dunning's own sequence of correlating d/f exponents is near-geometric in practice; that point
   is SECONDARY.
3. **A correlated-atomic criterion with a molecular-environment term (ANO-RCC style) for the polarization role.** The
   controls show water's RHF gain from the 3d is something the free-atom objective cannot see. An atomic objective
   that adds the O atom in a weak field, or O⁻/O²⁻, captures that role and stays water-free, as in the
   Widmark-Malmqvist-Roos / Roos 2004 averaging and the Ema VB2 sum over ions. The R3 diagnosis already ran an O²⁻
   probe, so the ingredients exist. A pre-registered variant of this kind would test whether the compensation in
   SL3d-T (RHF gain, correlation loss) can be removed by design.

## Credits owed (summary)
- sec:limits, the "different ladder" sentence, and the SL-3d record's criterion paragraph: cite Dunning 1989 (F1) and
  Woon-Dunning 1993 (F2) as the source of the s/p-from-HF, d-from-correlated-atom split.
- sec:prior, the ema2003 clause: state that their single-d STO sets for B-Ne took exponents from the atoms' HF plus
  frozen-core CISD, and H from H2 (F3). Add Prascher 2011 (F4) for the H2 rule behind R5H.
- Wherever the "polarization the atom does not need" (SL3d-H) point is written up: Jensen 2001.
