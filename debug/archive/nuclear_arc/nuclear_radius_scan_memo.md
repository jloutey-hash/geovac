# Nuclear radius scan: can the nucleus's own skeleton supply what the finite-nucleus doors consume?

Date: 2026-10-19. Read-only literature and corpus scan; one memo. Nothing in `geovac/`, `tests/` or `papers/` was modified. Scratch computations (exact sympy and mpmath) live in the session scratchpad, not the repo. They are reproducible from the formulas below.

**Verdict in one paragraph.** The PM's hunch is **half right**. In a harmonic-oscillator (HO) Slater determinant, including the Tassie–Barker centre-of-mass correction, every **even** charge moment is a rational multiple of b^{2k}. So ratios such as <r^4>/<r^2>^2 are rational and forced. The **odd** moments that the doors actually consume at odd order (the Friar moment <r^3>_(2) and the Zemach radius <r>_(2)) are **not rational**. They are (rational) x sqrt(rational) x pi^{-1/2} in units of b, because they come from the self-convolution of a Gaussian. This is the same pi^{-1/2} that Paper 23 already measured for its Gaussian profile. The scale b (or hbar*omega) is **free**. Shell filling at fixed hbar*omega gives R ∝ A^{1/6}, not A^{1/3}. The A^{1/3} law needs one consumed saturation number (r0, or equivalently rho0). No first-principles derivation of that number exists: chiral EFT consumes low-energy constants (LECs), and lattice QCD cannot reach saturated matter. Against data, the HO ratios hold at the **1–2% level for odd moments of s-shell nuclei** and are **~6% off for <r^4> of 4He**, because <r^4> weights the tail and the HO tail is Gaussian rather than exponential. Every candidate prediction is KNOWN or PROBABLY KNOWN. The honest door is an *accounting* door, like J14 and J16, not a finding.

---

## 1. Corpus state (Question 4)

- **Paper 23 `sec:finite_nucleus`** (lines 326–369) consumes, as inputs: <r^2> (leading defect -(2/3)Z^2<r^2>), the Friar moment (1/3)Z^3<r^3>_(2) at x^3, then <r^4> with a profile-dependent log constant at x^4, and for the hyperfine factor the Zemach radius r_Z. It already computes these for three model profiles: the uniform sphere (<r^4> = 3R^4/7, Friar 32R^3/21), a thin shell (Friar 16R^3/5), and **a Gaussian, whose Friar term carries pi^{-1/2}** (measured, 1.6e-6). That Gaussian is exactly the HO 0s density. So the corpus has already met the number class that an HO nucleus would bring, without naming it as the nucleons' own shell model. The paper itself notes that "the shared form is noted, not read": the electron's interior oscillator in a uniform sphere is not the nucleons' HO.
- **Seam rows J14 / J16** (`docs/forced_free_seam.md` l.163–164) and **atlas doors 24 and 26** (`docs/atlas.md` l.48, 50) price the radius and the joint moments as **C / FREE**, with "no label of either skeleton fixes them".
- **Paper 23 l.465–468:** the deuteron rms radius "~2.1 fm" is quoted from experiment, **not computed**. It serves only to explain why N_shells=2 at hbar*omega=10 MeV is unconverged.
- **Paper 24 l.1147:** b = 2.0364 fm at hbar*omega = 10 MeV, consistent with b^2 = (hbar c)^2/(m c^2 · hbar*omega) = 41.47/10 fm^2. This is a unit conversion, not a radius.
- **Code.** `geovac/nuclear/form_factor.py` hard-codes `R_PROTON_FM = 0.8414` and treats R as an input. `geovac/magnetization_density.py` takes r_Z as a calibrated input (`r_Z_bohr`, "calibrating the profile width"). `geovac/nuclear/harmonic_shell.py` has `harmonic_radial_wavefunctions` but **no moment or form-factor routine**. A grep for charge-radius or rms identifiers across `geovac/` and `tests/` finds only the electron-side PK and proton-constant uses.

**Conclusion:** GeoVac computes **no** nuclear charge moment from its own shell model. The machinery to do so (HO radial functions, Moshinsky brackets) is present, and the computation is a short script.

---

## 2. HO closed-shell charge moments (Question 1)

### 2.1 Exact results (computed here, sympy, exact arithmetic)

Point-proton density from filled HO orbitals with length b, form factor F(q) = (polynomial in q^2 b^2) x e^{-q^2 b^2/4}. The intrinsic (Tassie–Barker) form factor multiplies by e^{+q^2 b^2/(4A)}. Moments come from the Taylor coefficients of F. Zemach and Friar moments come from the self-convolution, whose Fourier transform is F^2.

| Nucleus (proton occupancy) | F_shell(q)·e^{q²b²/4} | <r²>_shell/b² | <r²>_intr/b² | <r⁴>/<r²>² shell / intr | zero(s) of F, qb |
|---|---|---|---|---|---|
| 4He (0s²) | 1 | 3/2 | 9/8 | 5/3 / 5/3 | none |
| 12C (0s² 0p⁴) | 1 − q²b²/9 | 13/6 | 49/24 | 255/169 / 3575/2401 | 3 |
| 16O (0s² 0p⁶) | 1 − q²b²/8 | 9/4 | 69/32 | 40/27 / 775/529 | 2√2 |
| 40Ca (to 1s0d) | 1 − q²b²/4 + q⁴b⁴/80 | 3 | 237/80 | 17/12 / 26405/18723 | q²b² = 10 ∓ 2√5 |

General closed shells up to N_max = M (one species): **<r²>/b² = 3(M+2)/4**, with A_species = (M+1)(M+2)(M+3)/3. Intrinsic: subtract 3b²/(2A).

Odd moments (intrinsic density, with charge distribution = magnetization distribution for r_Z):

| Nucleus | <r>_(2)/b (Zemach) | <r³>_(2)/b³ (Friar) | r_Z/r_E | Friar/r_E³ |
|---|---|---|---|---|
| 4He | √6/√π | 3√6/√π | **4/√(3π) = 1.30294** | **32√3/(9√π) = 3.47451** |
| 12C | 1337√66/(3267√π) | 643√66/(297√π) | 1.31276 | 3.40153 |
| 16O | 281√30/(450√π) | 139√30/(40√π) | 1.31410 | 3.39150 |
| 40Ca | 1109227√195/(3855735√π) | 239563√195/(109850√π) | 1.31682 | 3.36957 |

**Mechanism of the class split.** Any HO density is (polynomial) x Gaussian, and its self-convolution is (polynomial) x Gaussian of width √2·b. An even moment is a Gamma function at integer-plus-1/2 divided by Γ(3/2), which is rational. An odd moment is Γ(integer)/Γ(3/2), which is rational/√π, times 2^{odd/2}. The CM factor rescales the Gaussian width by (A−1)/A, so odd moments pick up a further √((A−1)/A). In Paper 18 / Paper 34 terms, π^{-1/2} is the Gaussian's normalisation surviving an odd radial projection. For the 4He 0s Gaussian it is the same π^{-1/2} that Paper 23 measured for its Gaussian profile. This is a reading, to be tagged at Paper 18's tier if it is ever written up.

### 2.2 What the literature says, and verification status

- **Gaussian checks against a primary source (VERIFIED, Sick, Phys. Rev. C 90, 064002 (2014), arXiv:1412.2603, read in full text).** Sick quotes, "for Gaussian (exponential) densities ... for the same rms radius", Friar(4He) = 16.50 (17.99) fm³. My closed form 32√3/(9√π)·1.681³ = **16.504 fm³** reproduces his Gaussian figure. For 3He he quotes a Gaussian Zemach radius "with both radii set to the experimental charge radius" of 2.570 fm. My 4/√(3π)·1.973 = **2.5706 fm** reproduces it. So the Gaussian odd-moment closed forms are KNOWN, and the code here is validated against them.
- **HO p-shell form factors.** The 0p radial integral J₀ = (1 − 2z/3)e^{−z}, z = q²b²/4, and the Tassie–Barker Gaussian factorisation were seen **only via search summaries of citing works** (Iraqi J. Phys. papers on 12C/16O HO form factors). The textbook sources named in the task (Bohr–Mottelson Vol. I; de Shalit–Talmi; Brussaard–Glaudemans; Elton; deForest–Walecka; Donnelly–Sick) were **not accessed**, and I give **no equation numbers** for them. The closed forms in the table are standard. Treat them as textbook-KNOWN, with the exact primary location owed.
- **Tassie–Barker.** The CM factor e^{q²b²/(4A)} is exact for a pure HO Slater determinant (0ħω configuration). It is a Gaussian, so it keeps even ratios rational and odd ratios in the √·π^{-1/2} class. The Kirson form of the charge radius (seen via citing work, Xayavong–Lim arXiv:2311.14967 Eq. (3), quoting Kirson, Nucl. Phys. A781, 350 (2007)) is <r²>_ch = <r²>_pt − 3b²/(2A) + r_p² + (N/Z)r_n² + 3ħ²/(4m²c²) + Δ_ls. This is **verified in Xayavong–Lim's text**, not in Kirson's.
- **Nucleon finite size breaks rationality.** Folding with the proton and neutron charge distributions adds new consumed numbers: r_p, r_n², and the higher nucleon moments. The ratios are then functions of r_p/b. They stay rational only if the nucleon density is itself modelled as a Gaussian, and the physical proton is close to dipole/exponential (<r⁴>/<r²>² = 5/2). So *at the level of the measured charge density*, no ratio is a pure rational. The forced statement lives on the **point-nucleon, intrinsic** density.

---

## 3. The scale (Question 2)

- **Practice.** The scale is fitted to radii. Blomqvist–Molinari, Nucl. Phys. A106, 545 (1968): **b² = 0.90 A^{1/3} + 0.70 fm²**. Kirson 2007: b_π² = 0.983(4) A^{1/3} + 0.373(23) fm². Both are seen via Xayavong–Lim arXiv:2311.14967, Eqs. (2) and (4), which I read; the originals were not read. Xayavong–Lim themselves refit b to the Angeli–Marinova (ADNDT 99, 69 (2013)) radii, with Coulomb and symmetry terms added.
- **The ħω ≈ 41 A^{-1/3} MeV rule.** This is a matching, not a derivation. Set the closed-shell HO <r²> = (3/4)(3A/2)^{1/3} b² (large-A limit of 3(M+2)/4, N=Z) equal to the uniform sphere's (3/5) r0² A^{2/3}. That gives b² = (4/5)(2/3)^{1/3} r0² A^{1/3} = 0.699 r0² A^{1/3}, and with (ħc)²/mc² = 41.47 MeV fm², **ħω = 41.2 A^{-1/3} MeV at r0 = 1.2 fm** (derived here). Bohr–Mottelson is cited for this relation (Eq. 2-157, "≈ 40 A^{-1/3}") by arXiv:nucl-th/0402046; I saw that **only via a search summary** and did not verify it in Bohr–Mottelson.
- **What shell filling alone forces.** At fixed ħω, filled HO shells give <r²> = b²·3(M+2)/4 ∝ A^{1/3}, so **R ∝ A^{1/6}**. The empirical A^{1/3} (saturation) is **not** a shell-filling consequence. It enters through ħω ∝ A^{-1/3}, which is one consumed number (r0 ≈ 1.2 fm, or ρ0 ≈ 0.16 fm^{-3}) plus the assumption of saturation. Given that rule, the *shell-structure* fluctuations around the A^{1/3} trend (kinks at closures, the 3(M+2)/4 staircase, the 3b²/(2A) CM term) **are** forced. The trend and its coefficient are not.
- **First principles.**
  - Chiral EFT reproduces radii only with fitted LECs. NNLOsat (Ekström et al., Phys. Rev. C 91, 051301(R) (2015), arXiv:1502.04682) fitted two- and three-nucleon forces simultaneously, including radii of selected C and O isotopes, and then got radii up to 40Ca and a nearly empirical saturation point. Earlier chiral interactions underpredicted the radii of 12C and 16O (the "radius problem"). Seen via search summary of the abstract.
  - Lattice EFT: a four-parameter SU(4)-type interaction reproduces neutron matter to saturation and ground states up to Ca (Lu et al., arXiv:1812.10928, seen via search summary).
  - Lattice QCD: saturated nuclear matter is out of reach because of the sign problem (Cohen et al., arXiv:1312.3339, seen via search summary). Its role is to fix EFT LECs.
  - **No first-principles derivation of r0 or ρ0 without consumed parameters was found.**
- **For the electron's door.** ZR in atomic units is roughly Z·(R m_π)·(m_e/m_π)·α. Every door number is then a product of one QCD-side dimensionless number (R m_π ≈ 0.85 A^{1/3}, fitted or LEC-derived) and two free constants the corpus already consumes (α, m_e/m_π). A nuclear-skeleton door cannot reduce this below one consumed nuclear scale.

---

## 4. Do the ratios match data? (Question 3)

Here "HO" means the 0s Gaussian (intrinsic, CM-corrected), with b fixed by the measured charge radius. It is folded with a dipole proton (r_p = 0.8409 fm) and a Galster neutron, using mpmath and Sick's Friar integral formula. Data are Sick 2014 (VERIFIED, primary).

| Quantity | Pure Gaussian | HO + nucleon size | Measured | HO deviation |
|---|---|---|---|---|
| 4He r_ch | (input) 1.681 | (input) 1.681 | 1.681(4) fm | — |
| 4He <r⁴>/r_ch⁴ | 5/3 = 1.667 | 1.696 | 14.35(11)/1.681⁴ = **1.797(15)** | **−5.6%** (~7σ) |
| 4He Friar/r_ch³ | 3.4745 | 3.486 | 16.73(10)/1.681³ = **3.522(21)** | **−1.0%** (~1.7σ) |
| 3He r_Z (r_M = r_E) | 2.5706 fm | — | **2.528(16) fm** | +1.7% |
| 3He r_M/r_E | 1 (HO, point nucleons) | — | 1.976(47)/1.973(14) = 1.002(25) | consistent |

For contrast, the uniform sphere gives <r⁴>/<r²>² = 25/21 = 1.19 and Friar/R³ = 3.28. The measured 4He sits between the Gaussian (HO) and the exponential, closer to the HO for the Friar moment.

**Reading.** The HO captures the odd moments of s-shell nuclei at the 1–2% level, because they are dominated by the bulk. It misses <r⁴> by about 6%, because higher moments weight the tail. The physical tail is exponential with decay constant √(2 m S) set by the separation energy S, which is another consumed number, and the HO tail is Gaussian. Reinhard, Nazarewicz and Garcia Ruiz (Phys. Rev. C 101, 021301 (2020), arXiv:1911.00699; ar5iv text read) tie <r⁴> to surface thickness through the Helm model. In this language that is a second, independent scale which a one-scale HO cannot carry. So **"forced ratios" is a model statement for <r⁴>/<r²>²**: the forced value is falsified by data at the 6% level. **For Friar/r_E³ and r_Z/r_E it is approximately physical**, within a few percent.

**Owed (not located this session).** Measured <r⁴> and form-factor zero positions for 12C, 16O and 40Ca. Neither the Fricke et al. / de Vries–de Jager–de Vries (ADNDT 36, 495 (1987)) compilations nor a table of <r⁴> values was reached. HO with b fitted to r_ch gives first zeros at **12C: 1.83 fm⁻¹; 16O: 1.61 fm⁻¹; 40Ca: 1.19 fm⁻¹**. The r_ch values used (2.4702, 2.6991, 3.4776 fm) are recalled from Angeli–Marinova and were **not verified** this session. The search returned conflicting unsourced claims for the measured 16O and 40Ca minima, which I am not reporting.

---

## 5. Candidate pre-registerable predictions (Question 5)

The door would read "the nucleus's own HO skeleton, Paper 23 Sec. shell plus Paper 24's π-free lattice, as the deeper theory for J14/J16's matching numbers". Every candidate below is KNOWN or PROBABLY KNOWN. This would be an accounting door, consistent with how J14 and J16 were priced.

**NR1 — Number-class split on the HO skeleton.** For any HO Slater determinant (closed shell or not, CM-corrected by Tassie–Barker), the intrinsic point-nucleon <r^{2k}>/b^{2k} ∈ ℚ, and the Zemach and Friar moments are ∈ √ℚ · π^{-1/2} · ℚ in units of b. So J14's x² and x⁴ coefficients are rational in Zb, and the x³ (Friar) coefficient carries exactly one π^{-1/2}. Quantity: exact symbolic moments for 4He, 12C, 16O, 40Ca, plus one open-shell case (7Li).
*Falsifier:* any even moment irrational in b, or an odd moment outside √ℚ·π^{-1/2}.
*Controls:* positive, reproduce Sick's Gaussian 16.50 fm³ and 2.570 fm (done here); negative, the exponential profile, whose odd moments are rational.
*Novelty:* **KNOWN** for the Gaussian (Sick 2014, Sec. "Moments for 4He", verified). **PROBABLY KNOWN** for general HO shells: it is Γ-function arithmetic, and no explicit statement of the class split was found in the searches here (Sick, the Reinhard et al. <r⁴> paper, generic HO form-factor searches). A hit is a CHECK.

**NR2 — Forced dimensionless ratios vs data, s-shell.** HO (intrinsic Gaussian, b from r_ch, nucleon folding with r_p consumed) predicts 4He Friar/r_ch³ = 3.486 ± (nucleon-model spread) and 3He r_Z/r_E = 4/√(3π) at the Gaussian-nucleon level.
*Falsifier:* deviation from Sick 2014 larger than 3% for either.
*Pre-disclosed:* the computation here already gives −1.0% and +1.7%, so this is no longer blind.
*Novelty:* **KNOWN**, since Sick quotes the Gaussian estimates beside his data. CHECK only.

**NR3 — Scale-free form-factor invariant for p-shell and sd-shell closures.** The first zero q₀ times the intrinsic point rms is forced: q₀²<r²>_intr = **147/8** (12C), **69/4** (16O), **(10 − 2√5)·237/80** (40Ca). The scale b cancels, so this is the one forced, scale-free observable of the HO skeleton visible directly in electron-scattering data. Proton folding does not move zeros, but it does move the extracted point rms.
*Falsifier:* deviation larger than 5% from the de Vries compilation's Fourier–Bessel zeros and r_pt. I expect it to hold for 12C and 16O and fail for 40Ca, where HO form factors are known to be poor.
*Novelty:* the HO form-factor zeros are **KNOWN** (textbook p-shell form factor). The b-free product with the CM-corrected rms is **NOT FOUND** stated as such in this session's searches (HO form-factor and diffraction-minimum queries only). That search was shallow, so "PROBABLY KNOWN" is the safer label. Requires the 1987 compilation, which is owed.

**NR4 — What shell filling forces in the A-dependence.** At fixed ħω, R ∝ A^{1/6}, the staircase 3(M+2)/4, and a −3b²/(2A) CM term are forced. The A^{1/3} trend needs exactly one consumed number (r0 or ρ0) plus saturation. Quantity: the residuals of Angeli–Marinova radii for doubly-closed nuclei against [one-parameter b² = c·A^{1/3}] versus Blomqvist–Molinari's two-parameter fit.
*Falsifier:* the one-parameter (pure saturation) form fits doubly-magic radii as well as the two-parameter form. That would mean the 0.70 fm² offset is not needed.
*Novelty:* **KNOWN** (Blomqvist–Molinari 1968; Kirson 2007; Xayavong–Lim 2023 all do this fit). CHECK only.

**Recommendation.** If the PI wants a door, NR1 is the cleanest atlas entry. It turns J14's "the profile's number class" into a *forced* class for the skeleton's own profile: even orders rational, Friar order π^{-1/2}, with the single consumed scale b. NR3 is the only candidate with a plausible, if thin, novelty claim, and it needs the de Vries compilation before registration. None of the four derives R. The literature's own answer is that R's scale is a consumed saturation number, not forced by any skeleton or by first-principles QCD at present.

## Sources and verification status

- VERIFIED, primary text: Sick, PRC 90, 064002 (2014), arXiv:1412.2603; Xayavong & Lim, arXiv:2311.14967 (for its Eqs. 2–4 and bibliography); Reinhard, Nazarewicz & Garcia Ruiz, arXiv:1911.00699 (Helm/surface-thickness relation only).
- Via citing work only: Blomqvist & Molinari, NPA 106, 545 (1968); Kirson, NPA 781, 350 (2007); Angeli & Marinova, ADNDT 99, 69 (2013).
- Via search summary only: Bohr–Mottelson Vol. I Eq. 2-157 (cited by arXiv:nucl-th/0402046); Ekström et al., PRC 91, 051301 (2015); Lu et al., arXiv:1812.10928; Cohen et al., arXiv:1312.3339; Diepold et al., Ann. Phys. 396, 220 (2018); HO p-shell radial integral and Tassie–Barker factorisation.
- Not accessed: de Shalit–Talmi, Brussaard–Glaudemans, Elton, deForest–Walecka, Donnelly–Sick, the de Vries compilations, and Tassie–Barker's original paper. No equation numbers are claimed for any of them.
