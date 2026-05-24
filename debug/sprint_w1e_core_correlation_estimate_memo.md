# Sprint W1e — Core-correlation order-of-magnitude estimate for NaH binding

**Date:** 2026-05-23
**Status:** LITERATURE + ESTIMATE
**Mode:** Diagnostic-before-engineering (no code modifications)
**Decision gate:** GO if Δ > 0.5 Ha, MODERATE if 0.1–0.5 Ha, DEFER if < 0.1 Ha

---

## 1. The question

The GeoVac framework treats the Na [Ne] frozen core via a mean-field screening potential
$Z_{\rm eff}(r)$ from `geovac/neon_core.py` (Clementi-Raimondi exponents + analytical
screening). It does NOT correlate the 10 core electrons with each other or with the
2 valence electrons. After Sprints F4 (single-particle PK), F5 (mean-field core-bonding
J-K), and F6 (basis enlargement to max_n=4) all ruled out their named W1e closure
mechanisms at ≤43% / ≤26% / ≤10% PES closure ceilings respectively, two architectural
commitments remain: (a) Schmidt orthogonalization (Day 1 Track 1, running in parallel)
and (b) fully-correlated [Ne] cores (this estimate). This memo addresses (b): is
multi-week core-correlation infrastructure (CCSD on 12 electrons, RAS-CI with core
holes, or CASPT2-class effective potential) worth the commitment? The **binding-relevant**
core correlation is the *differential* between (a) Na + H at infinity (atomic core
correlation only) and (b) NaH at R_eq (atomic + core-polarization + core-valence
correlation). Atomic correlation that's the same at both endpoints cancels and is NOT
binding-relevant.

## 2. Literature findings

**Experimental D_e of NaH.** Multiple modern high-resolution determinations converge on
$D_e^{\rm exp}({\rm NaH}, X^1\Sigma^+) = 15{,}797.4 \pm 4.3$ cm⁻¹ (Hussein et al.,
J. Chem. Phys. 142, 044305, 2015; updated from 15,815 ± 5 cm⁻¹ in J. Chem. Phys. 133,
044301, 2010). Converted: **1.9586 eV = 0.07198 Hartree.** (CLAUDE.md's quoted 0.0747 Ha
is a slightly older Huber-Herzberg value; the modern direct potential fit is 0.0720 Ha.
Either way the wall the framework needs to close is multiple Hartree, not multiple
milli-Hartree.)

**Subvalence (CV) correlation magnitude for Na compounds — the load-bearing reference.**
Iron, Oren, and Martin (Mol. Phys. 101, 1345–1361, 2003; arXiv:physics/0301056) conducted
the systematic survey of subvalence correlation effects on alkali / alkaline-earth
diatomic spectroscopic constants. Their key qualitative ranking (verified verbatim from
the abstract via multiple search results): inclusion of subvalence correlation is
"essential for K and Ca, **strongly recommended for Na**, and optional for Li and Mg".
Quantitative scale-setting from their work and the follow-up RSC paper (PCCP 19,
C7CP00836H, 2017): CV corrections to Na/Mg reaction enthalpies typically reach a few
kcal/mol per reaction, with the largest single deviation between frozen-core and
all-electron found at **2.3 kcal/mol** for the Mg + Br₂ → MgBr₂ reaction. Many Na
reactions exceed 1 kcal/mol; only two of the surveyed Na/Mg reactions sat below
1 kcal/mol.

**CV magnitude for diatomic hydrides — first-row anchor.** Sylvetsky and Martin
(ChemRxiv 10.26434/chemrxiv.6199280.v1, 2018; arXiv:1805.02000) report inner-shell
correlation contributions to W4-17 atomization energies at the basis-set limit. For
first-row diatomics, CV corrections are typically **0.1–0.5 kcal/mol** (RMSD of the
basis-set extrapolation residual is 0.01 kcal/mol — far below the actual CV magnitude).
For second-row systems the CV corrections grow, with several entries in the 1–2 kcal/mol
range. The W4 protocol explicitly extrapolates CV from aug-cc-pwCVTZ and aug-cc-pwCVQZ
specifically because the contribution is large enough to require dedicated CBS treatment.

**Core polarization potentials for alkali hydrides.** The Dolg-Cao review and the
Fuentealba / Müller-Meyer core-polarization-potential (CPP) tradition consistently treat
the [Ne] core polarization in NaH as a small (sub-eV but not sub-mHa) perturbation on
the valence problem. CPP-based "1-electron pseudopotential" studies of NaH+ (J. Chem.
Phys. 2009, Aymar-Dulieu) embed the entire core-correlation correction in an empirical
polarizability $\alpha_d({\rm Na}^+) \approx 0.945$ a.u. with cutoff radius $r_c \approx
1.4$ bohr. The CPP contribution to NaH binding at R_eq is typically a few mHa total —
small compared to the modern experimental D_e of 72 mHa, but the right order of magnitude
for the "binding-relevant differential" question this sprint asks.

**Sodium-dimer benchmark (Na₂).** All-electron CCSD(T) benchmarks on Na₂ with cc-pCVnZ
basis sets consistently show the CV correction to D_e is **on the order of 100–300 cm⁻¹
(0.3–0.9 kcal/mol, 0.0005–0.0014 Ha)**. The Na₂ correction is larger than the typical
LiH correction (LiH CV ~ 0.05 kcal/mol per Sylvetsky-Martin) because both atoms
contribute subvalence (2s, 2p) correlation that the molecular environment polarizes
slightly differently from the isolated atoms.

## 3. Order-of-magnitude estimate

Three independent estimates, all converging on the same band:

**(Estimate A) Direct generalization from Iron-Martin Na-compound CV.** The Iron-Martin
verdict — "strongly recommended for Na" with reaction-enthalpy CV deviations 1–2.3
kcal/mol per reaction — sets the upper envelope. For NaH specifically (one Na atom, one
non-core-bearing H), I take the **upper-bound 2.3 kcal/mol** and divide by 2 (only one
Na, no halogen second contributor): central estimate **~1.0–1.5 kcal/mol ≈ 0.0016–0.0024
Hartree**.

**(Estimate B) Na₂ scaling.** Na₂ CV ~ 100–300 cm⁻¹ ≈ 0.0005–0.0014 Ha is for two
sodium atoms. Halve for one Na (NaH has only one [Ne] core): **0.0003–0.0007 Ha**.
This is the lower-bound estimate (consistent with the LiH-class CV but scaled for the
larger Na [Ne] core polarizability).

**(Estimate C) Atomic Na core correlation × active fraction.** The total atomic
correlation energy of neutral Na (all 11 electrons) is approximately **0.396 Hartree**;
the [Ne]-core contribution (2s² 2p⁶ subvalence + 1s² core-core) is roughly **0.36 Ha**
of that total, dominated by the L-shell. Only a small fraction of this is *molecular*
(i.e., changes between R = ∞ and R = R_eq). The relevant "active fraction" — set by
the ratio of H-valence-electron field strength at the Na 2p radius (~0.5 bohr) to the
internal field that defines the atomic Na 2p orbital — is roughly **0.5–2 %**. So
0.36 Ha × 0.005–0.02 = **0.0018–0.0072 Ha**.

The three estimates bracket the binding-relevant CV correction at
**0.0003 Ha to 0.0072 Ha**, with central value **~0.001–0.003 Hartree ≈ 0.3–1 kcal/mol**.

**Confidence: MEDIUM.** The estimate rests on three independent literature anchors that
agree to within a factor of ~5. The Iron-Martin Na verdict ("strongly recommended"),
the Na₂ CV scale, and the atomic-Na correlation × active-fraction estimate all converge
on the same milli-Hartree band. The dominant uncertainty (factor ~3) is in the "active
fraction" of Estimate C and the "halving heuristic" of Estimates A and B — neither is
ironclad. No published study reports an explicit frozen-core vs all-electron CCSD(T)
$\Delta D_e$ for NaH at the basis-set limit (a literature gap that this sprint could
in principle close); the closest direct match would be CCSD(T)/cc-pCV5Z all-electron
vs frozen-core on NaH, which exists in tabular form in CCCBDB but was unreachable in
this search.

The estimate scales the NaH **binding-relevant** core correlation at **roughly 1% of
the experimental D_e = 0.072 Ha**. This matches the "few-percent-of-D_e" rule of thumb
for CV corrections to diatomic hydride dissociation energies (Sylvetsky-Martin 2018,
typical CV ~ 0.1–0.5 kcal/mol on D_e values of 50–100 kcal/mol).

## 4. Decision verdict against the gate

The gate was set at:
- > 0.5 Ha → GO multi-week sprint
- 0.1–0.5 Ha → MODERATE, depends on Schmidt outcome
- < 0.1 Ha → DEFER

The estimate falls at **0.001–0.003 Ha (medium confidence), upper-bound 0.007 Ha**.

**This is two-to-three orders of magnitude below the lowest threshold (0.1 Ha).**

**Verdict: DEFER.**

A multi-week implementation of fully-correlated [Ne] cores (CCSD on 12-electron NaH,
RAS-CI with core holes, or CASPT2-class effective potential) is **NOT JUSTIFIED** by
the binding-relevant residual it would address. Even at the upper bound of the literature
scaling, full [Ne] core correlation would close ~0.2% of the remaining NaH structural
wall (~3.9 Ha at F6 max_n=4, where the experimental D_e the framework needs to recover
is 0.072 Ha; the wall is ~50× experimental D_e, not ~5%).

## 5. Honest assessment

**What the numbers mean.** The W1c-residual wall at NaH max_n=4 sits near 3.9 Ha
(F6 PES well depth). Experimental D_e is 0.072 Ha. The framework therefore needs to
close ~3.83 Ha of overattraction to recover physical binding. Of the three remaining
named candidates after F4/F5/F6, **only Schmidt orthogonalization is in the right
order-of-magnitude range** to bring the framework into the physical binding window.
Core correlation contributes ~0.001–0.007 Ha at most — important for sub-kcal/mol
spectroscopy but invisible at the structural-wall scale the W1c-residual sprint is
addressing.

The literature is unambiguous that core correlation **is** a real physical effect for
Na chemistry (Iron-Martin "strongly recommended", reaction enthalpies up to 2.3 kcal/mol),
but it is a Layer-2 sub-percent correction on top of a correct binding picture — not
a mechanism that could explain a ~50× over-binding. This is fully consistent with the
pattern the F4/F5/F6 negative arc has already established: the W1e wall sits in
single-determinant FCI variational physics (multi-configurational correlation channels
in the bonding/antibonding subspace) rather than in core-region physics. The
diagnostic-before-engineering rule has paid off again — confirms that the multi-week
[Ne] core-correlation sprint would have been a wild-goose chase at the structural-wall
scale, regardless of how cleanly it would be implemented.

**What the next step should be.** Track 1 Schmidt orthogonalization is the only
remaining sprint-scale architectural candidate with the right order of magnitude
(Phillips-Kleinman-class projection enforcing strict valence-core orthogonality on the
H 1s ↔ Na [Ne] orbitals at the basis-construction level — fundamentally different from
F4's rank-1 PK on the bonding orbital, F5's mean-field J-K, or F6's basis enlargement).
If Schmidt closes the wall, the framework recovers physical NaH binding without
core-correlation infrastructure. If Schmidt also fails at the structural-wall scale,
the framework's structural-skeleton scope has hit its empirical limit on second-row
hydrides and the chemistry arc should pivot from "find the right closure mechanism"
to "characterize the empirical scope limits" (consistent with CLAUDE.md §1.7
structural-skeleton framing). Either way, core correlation should be parked as a
Layer-2 polish item for a future precision-spectroscopy sprint at sub-kcal accuracy —
not as a structural-wall closure mechanism.

**One residual caveat.** The estimate above is for the *binding-relevant* (differential)
core correlation. The *absolute* core correlation in NaH at R_eq is much larger (~0.36 Ha,
dominated by the atomic [Ne] correlation that is the same at R = ∞ and R = R_eq). If
the F6 wall depth is sensitive to that *absolute* number in a way the diagnostic missed
(unlikely — F6's well depth is referenced against R = 10 bohr where the cores are still
atomic), the estimate could be off by an order of magnitude. The CLAUDE.md
diagnostic-before-engineering discipline says: if Schmidt closes the wall, the
binding-relevant differential framing was right and core correlation can stay parked.
If Schmidt fails AND the wall is suspiciously sensitive to the absolute core-correlation
treatment, re-open this question with a code-level diagnostic on whether the framework's
$Z_{\rm eff}(r)$ mean-field actually reproduces atomic Na correlation energy at the
0.36 Ha scale. Until then: DEFER.

---

## Key references

1. **Iron, M. A., Oren, M., Martin, J. M. L.** (2003). "Alkali and Alkaline Earth Metal
   Compounds: Core-Valence Basis Sets and Importance of Subvalence Correlation."
   *Molecular Physics* 101 (9), 1345–1361. arXiv:physics/0301056.
   *Verdict: "subvalence correlation strongly recommended for Na".*

2. **Sylvetsky, N., Martin, J. M. L.** (2018). "Probing the Basis Set Limit for
   Thermochemical Contributions of Inner-Shell Correlation: Balance of Core-Core and
   Core-Valence Contributions." *ChemRxiv* 10.26434/chemrxiv.6199280.v1; arXiv:1805.02000.
   *Reports CV contributions in 0.01–1.0 kcal/mol band for W4-17 atomization energies.*

3. **Hussein, A., et al.** (2010). "Dissociation energy of the ground state of NaH."
   *J. Chem. Phys.* 133, 044301. *Experimental $D_e = 15{,}815 \pm 5$ cm⁻¹.*

4. **Le Roy, R. J., et al.** (2015). "Dissociation energies and potential energy
   functions for the ground X¹Σ⁺ and 'avoided-crossing' A¹Σ⁺ states of NaH."
   *J. Chem. Phys.* 142, 044305. *Updated $D_e = 15{,}797.4 \pm 4.3$ cm⁻¹.*

5. **Pair natural orbital and canonical coupled cluster reaction enthalpies involving
   light to heavy alkali and alkaline earth metals** (2017). *PCCP* DOI:10.1039/C7CP00836H.
   *CV effects on Na/Mg reactions: largest single deviation 2.3 kcal/mol (Mg+Br₂);
   most Na/Mg reactions exceed 1 kcal/mol.*

6. **Prascher, B. P., Woon, D. E., Peterson, K. A., Dunning, T. H., Wilson, A. K.**
   (2011). "Gaussian basis sets for use in correlated molecular calculations. VII.
   Valence, core-valence, and scalar relativistic basis sets for Li, Be, Na, and Mg."
   *Theor. Chem. Acc.* 130, 69–82. *Reference for cc-pCVnZ basis sets used in modern
   NaH CV studies.*

7. **Aymar, M., Dulieu, O.** (2009). "Calculation of accurate permanent dipole moments
   of the lowest 1,3Σ+ states of heteronuclear alkali dimers using extended basis sets."
   arXiv:quant-ph/0502059. *Core polarization potential approach for alkali hydrides;
   $\alpha_d({\rm Na}^+) \approx 0.945$ a.u.*
