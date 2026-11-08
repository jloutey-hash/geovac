# Confidence Review: Paper 23 — Nuclear Shell Model Hamiltonians on the Hyperspherical Lattice

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | HO shell closures 2, 8, 20, 40, 70, 112 from (N+1)(N+2) degeneracies | Sec II.A, after eq (1) | **A** | EXTERNAL (standard QM) | Python sum of (N+1)(N+2) for N=0..5 = 2,8,20,40,70,112 exactly. |
| 2 | Magic numbers 2,8,20,28,50,82,126 recovered with v_ls/ℏω≈0.171, d_ll/ℏω≈0.021 | Eq (3), Table I | **B** | GEOVAC-ONLY (the scan), but values consistent with Suhonen 2007 / Talmi 1993 | Could not independently re-run the scan; level ordering in Table I matches standard Mayer-Jensen ordering. |
| 3 | Fock projection rigidity theorem (S³ unique to -Z/r) | Sec III, theorem statement | **A** | EXTERNAL (Bander-Itzykson 1966; SO(4) of hydrogen) | Standard textbook result, cited correctly. |
| 4 | Deuteron: 16 qubits, 592 non-identity Pauli, 1-norm ≈227 MeV | Table II | **B** | GEOVAC-ONLY (numerical output of the framework's builder) | Cannot independently verify Pauli count; consistent with Track NI Zenodo memo (also 592 Pauli claim). |
| 5 | He-4: 16 qubits, 712 Pauli, Hilbert dim 784 (= C(8,2)²) | Table III | **B** + **A** for dim | EXTERNAL for dim, GEOVAC-ONLY for Pauli count | C(8,2)² = 28² = 784 verified bit-exact; 784/64 = 12.25 verified; 712/592 = 1.20 verified. |
| 6 | Deuteron experimental B = 2.2246 MeV; rms charge radius ~2.1 fm | Sec IV.A | **A** | EXTERNAL | PDG-consistent; modern values 2.224566 MeV, 2.128 fm. |
| 7 | He-4 experimental B = 28.296 MeV | Sec V opening | **A** | EXTERNAL | Standard nuclear data table value. |
| 8 | Composed nuclear-electronic: 26 qubits, 614 Pauli (592+10+12), 13-OoM scale ratio | Sec VI.B, Table IV | **B** | GEOVAC-ONLY | Internally consistent with Track NI memo. Cannot independently verify Pauli count. |
| 9 | A_hf(H 21cm) ≈ 1420.4 MHz ≈ 2.16×10⁻⁷ Ha; max shift 3A/4 ≈ 1.62×10⁻⁷ Ha | Sec VI.B | **A** | EXTERNAL | Recomputed: 2.1588×10⁻⁷ Ha, 0.75·A = 1.6191×10⁻⁷ Ha — verified. |
| 10 | Finite-size shift ΔE₁ₛ = (2/5) Z⁴ R²/n³ ≈ 1.01×10⁻¹⁰ Ha at R = 0.8414 fm | Sec VI.B | **A** | EXTERNAL (Foldy-Friar formula) | Recomputed: 1.0113×10⁻¹⁰ Ha — verified. |
| 11 | Roothaan J₀(λₑ,λₙ) closed form, gives 5/8 at λ=1, symmetric | Eq (15) | **A** | EXTERNAL (Roothaan 1951) | Sympy verification: J₀(1,1) = 5/8 exact; symmetry verified to zero residual. |
| 12 | Taylor expansion 1 − 2ε² + 5ε³ − 9ε⁴ + 14ε⁵ − 20ε⁶ + 27ε⁷ − 35ε⁸ | Eq (17) | **A** | self-contained algebra | Reproduced bit-exact in sympy. |
| 13 | λₙ = 2√(M_p/m_e) ≈ 85.7 bohr⁻¹; ΔE_recoil = +2.65×10⁻⁴ Ha; BS = +2.72×10⁻⁴ Ha; residual 2.86% | Sec VII.B | **A** | EXTERNAL + numerical | Recomputed bit-exact: λ_n = 85.7007, ΔE_cross = 2.6453×10⁻⁴, BS = 2.7231×10⁻⁴, residual −2.86%. |
| 14 | Δν_Z/ν_F = -2 Z α mₑ r_Z gives -39.500 ppm at r_Z=1.045 fm | Eq (22), Sec VIII | **A** | EXTERNAL (Eides convention) | Recomputed using m_e c/ℏ as inverse-Compton-length: -39.50 ppm exact. **CAVEAT**: the formula's "mₑ" is shorthand for inverse Compton wavelength m_e c/ℏ; reading "mₑ" as mass-in-atomic-units gives -0.29 ppm. Eides convention is standard but notationally non-obvious. |
| 15 | Deuteron BF prediction 327.3975 MHz at +40 ppm vs Wineland-Ramsey 1972 | Sec X | **A** with **C-flag** | EXTERNAL | Numerical value verified bit-exact. **C**: eq (24) writes the result times "Ha" at the end, which is dimensionally awkward — actually means "the Hartree-as-frequency conversion factor 6.5797×10¹⁵ Hz". A reader recomputing eq (24) literally gets 5×10⁻⁸ Hz, not 327 MHz. |
| 16 | Cross-register V_eN at 0.012% residual vs Eides | Sec VIII | **B** | GEOVAC-ONLY (operator code) | Cannot verify the operator-level numerics independently. |
| 17 | "First two-species qubit encoding for the deuteron" | abstract | **D** with **C-flag** | NOVELTY (can't confirm) | I cannot positively verify priority. Roggero et al. 2020 and Romero et al. 2022 both built nuclear qubit Hamiltonians (deuteron-related). Recommend softening to "to our knowledge, the first GeoVac-style two-species encoding for the deuteron." |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed | Survives? |
|---|---|---|---|---|
| HO shell closures | 2,8,20,40,70,112 | (N+1)(N+2) sum, standard QM | 2,8,20,40,70,112 | yes |
| C(8,2)² | 784 | combinatorics | 784 | yes |
| Hilbert ratio 784/64 | 12.25 | arithmetic | 12.25 | yes |
| Pauli ratio 712/592 | 1.20 | arithmetic | 1.2027 | yes |
| g_d^atomic = 2·μ_d/μ_N | 1.7148 | CODATA μ_d/μ_N = 0.857 | 1.7149 | yes |
| Roothaan J₀(1,1) | 5/8 | exact rational | 5/8 | yes |
| Roothaan Taylor coefficients | {-2, 5, -9, 14, -20, 27, -35} | sympy | bit-exact match | yes |
| λ_n = 2√(M_p/m_e) | ~85.7 | CODATA M_p/m_e = 1836.15 | 85.7007 | yes |
| ΔE_recoil cross | +2.65×10⁻⁴ Ha | Z·(J₀-λ_e) | 2.6453×10⁻⁴ | yes |
| BS leading | +2.72×10⁻⁴ Ha | |E_1| m_e/m_p | 2.7231×10⁻⁴ | yes |
| Recoil residual | 2.86% | direct | 2.86% | yes |
| Finite-size ΔE₁ₛ | 1.01×10⁻¹⁰ Ha | (2/5)Z⁴R²/n³ | 1.0113×10⁻¹⁰ | yes |
| A_hf in Ha | 2.16×10⁻⁷ | 1420.4 MHz × 10⁶ / Ha_to_Hz | 2.1588×10⁻⁷ | yes |
| 3A/4 in Ha | 1.62×10⁻⁷ | 0.75 × A | 1.6191×10⁻⁷ | yes |
| Δν_Z/ν_F | -39.500 ppm | -2 α (m_e c/ℏ) r_Z (SI) | -39.50 ppm | yes (with unit caveat) |
| D BF HFS prediction | 327.3975 MHz | (3/2)(4/3) g_d α² (m_e/m_d) × Ha_to_Hz | 327.3975 MHz | yes |
| D BF residual vs Wineland-Ramsey | +40 ppm | direct | +40.1 ppm | yes |
| r_Z in bohr at 1.045 fm | 1.974×10⁻⁵ | division by a_0 = 52917.72 fm | 1.9748×10⁻⁵ | yes |

### Circularity map

**GEOVAC-ONLY chains (house-of-cards risk):**

1. **Pauli counts 592 / 712 / 614** — these are outputs of the framework's qubit builder code. No independent reproduction; I cannot verify whether the builder is correct. They are *internally consistent* across Paper 23 §IV/§V/§VI and the Track NI memo, but that is intra-corpus consistency, not external truth. Severity: standard; this is the framework's central output. A bug in the JW encoding (especially given the **acknowledged sign issue in Sec VI.D**, see below) could affect these counts.

2. **1-norm ≈ 227 MeV / 557 MeV / 552 MeV** — same status; outputs of the builder.

3. **Magic-number scan v_ls ≈ 0.171, d_ll ≈ 0.021** — the existence of *some* (v_ls, d_ll) reproducing the seven magic numbers is well-established (Mayer-Jensen 1949), and the numerical values are within the standard range. The midpoint-of-valid-region procedure is GEOVAC-specific but not load-bearing.

**MIXED chains:**

4. **Cross-register V_eN at 2.86% / 0.012% Eides match** — leading-order physics is EXTERNAL (Bethe-Salpeter recoil, Eides Zemach), the recovery of these values via the operator construction is GEOVAC-ONLY. The 2.86% recoil residual is now decomposed analytically (half-integer-power tower), which is a real algebraic finding — that piece is **A**.

### Overstatement findings

**C-1.** Abstract: "We build the **first** two-species qubit encoding for the deuteron." → **Softening recommended**: "We build a two-species qubit encoding for the deuteron." Roggero et al. 2020 (PRD 101, 074038) already computed deuteron-related quantum-circuit constructions; Romero et al. 2022 implemented ADAPT-VQE for nuclear structure. I cannot positively verify priority, and per the CONFIDENCE_REVIEW honest-ceiling rule, "first" is a claim that cannot be confirmed by web search.

**C-2.** Eq (24): the formula `ν_HFS = (3/2)·(4/3) g_d α² (m_e/m_d) Ha` is dimensionally non-obvious — a literal reading gives a number times a Hartree (energy). The intended meaning is "the Hartree expressed as a frequency", i.e., a conversion via E = hν. **Recommended fix**: replace the trailing "Ha" with an explicit "× Ha_to_Hz" or rewrite as a dimensionless ratio (or add a footnote stating the unit convention).

**C-3.** Sec VIII eq (22): `Δν_Z/ν_F = -2 Z α m_e r_Z`. The "m_e" here is shorthand for the inverse electron Compton wavelength m_e c/ℏ, not m_e as a mass in atomic units. Standard in Eides notation but a reader recomputing literally will get an answer off by α (factor ≈10⁻². The verification at -39.50 ppm only works with m_e read as inverse Compton length. **Recommended fix**: a footnote stating "in natural units m_e c/ℏ" or explicit conversion factor.

**C-4.** Sec VI.D acknowledges a "known sign issue" in the off-diagonal one-body Pauli encoding: `a†_i a_j` encoded as XY − YX instead of XY + YX. The paper claims this "does not affect" the deuteron FCI spectrum because the diagonal kinetic term is zero. **Caveat**: this is a real code-level bug acknowledged in the paper, and the assertion that it "does not affect Pauli counts" is a GEOVAC-internal claim. A referee would (legitimately) ask whether He-4 with antisymmetrization is also unaffected, or whether the cross-register hyperfine bypass (acknowledged in the paper) is the only saving grace. The disclosure is honest but I want to flag this as a known fragility.

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found |
|---|---|---|---|
| Mayer1949 | PR 75, 1969 (1949) — On Closed Shells in Nuclei II | **CITE-OK** | doi:10.1103/PhysRev.75.1969. |
| Jensen1949 | PR 75, 1766 (1949) — Magic Numbers in Nuclear Structure (Haxel/Jensen/Suess) | **CITE-OK** | doi:10.1103/PhysRev.75.1766.2. |
| Suhonen2007 | Springer textbook From Nucleons to Nucleus | **CITE-OK** | ISBN 978-3-540-48859-0. |
| Talmi1993 | Harwood textbook Simple Models of Complex Nuclei | **CITE-OK** | standard reference. |
| Thompson1977 | Nucl. Phys. A 286, 53 — Minnesota potential (Resonating-Group Method) | **CITE-OK** (likely) | Could not access full text behind paywall; Tang/Lemere/Thompson collaboration at Univ. of Minnesota in this era is well-attested; specific V_R=200, V_T=-178, V_S=-91.4 MeV are widely cited Minnesota parameters in subsequent literature. |
| Moshinsky1959 | Nucl. Phys. 13, 104 — Transformation Brackets | **CITE-OK** | Classic well-cited paper. |
| Talmi1952 | Helv. Phys. Acta 25, 185 — HO basis nuclear spectroscopy | **CITE-OK** (likely) | Pre-web journal; standard early reference. |
| Roggero2020 | PRD 101, 074038 — Quantum Computing for Neutrino-Nucleus Scattering | **CITE-OK** | doi:10.1103/PhysRevD.101.074038 confirmed. |
| Hagen2014 | RPP 77, 096302 — Coupled-Cluster Computations of Atomic Nuclei | **CITE-OK** | doi:10.1088/0034-4885/77/9/096302 confirmed. |
| Romero2022 | PRC 105, 064317 | **CITE-OK** | doi:10.1103/PhysRevC.105.064317 confirmed. |
| Navratil2000 | PRL 84, 5728 (Navratil-Vary-Barrett 12C) | **CITE-WRONG-METADATA (minor)** | Real title is "Properties of 12C in the Ab Initio Nuclear Shell Model" (not "Large-Basis Ab Initio NCSM and Its Application to 12C"). Authors and venue correct. Cosmetic only. |
| BarrettNavratilVary2013 | PPNP 69, 131 — Ab Initio NCSM review | **CITE-OK** | confirmed. |
| Fock1935 | Z. Phys. 98, 145 | **CITE-OK** | doi:10.1007/BF01336904 confirmed. |
| Bander1966 | RMP 38, 330 — Group Theory and Hydrogen Atom | **CITE-OK** | doi:10.1103/RevModPhys.38.330 confirmed. |
| Bargmann1961 | Commun. Pure Appl. Math. 14, 187 | **CITE-OK** | doi:10.1002/cpa.3160140303 confirmed. |
| Foldy1958 | RMP 30, 471 — Neutron-Electron Interaction | **CITE-OK** | doi:10.1103/RevModPhys.30.471 confirmed. |
| Friar1979 | Ann. Phys. 122, 151 | **CITE-OK** | doi:10.1016/0003-4916(79)90300-2 confirmed. |
| Roothaan1951 | J. Chem. Phys. 19, 1445 | **CITE-OK** | doi:10.1063/1.1748100 confirmed. |
| EidesGrotchShelyuto2007 | Springer textbook Theory of Light Hydrogenic Bound States | **CITE-OK** | ISBN 978-3-540-45269-8 standard reference. |
| ChamseddineConnes2010 | Fortschr. Phys. 58, 553 | **CITE-OK** | doi:10.1002/prop.201000069 standard reference. |
| WinelandRamsey1972 | PRA 5, 821 — Atomic Deuterium Maser | **CITE-OK** | doi:10.1103/PhysRevA.5.821 confirmed; experimental value 327.384352522 MHz confirmed. |
| **Pachucki2023** | **PRL 130, 023004 — "Three-Photon Exchange Nuclear Structure Correction in Hydrogenic Systems"** | **CITE-MISATTRIBUTED** | **HIGH SEVERITY.** PRL 130, 023004 (2023) is actually "Recoil Corrections to the Energy Levels of Hydrogenic Atoms" by Pachucki & Yerokhin (no Patkóš). The "Three-Photon Exchange" paper (Pachucki/Patkóš/Yerokhin) is arXiv:1803.10313 and was published in **PRL 121, 073001 (2018)**. The paper's use of Pachucki2023 for "Foldy-Wouthuysen-reduced two-particle Hamiltonian" matches the actual content of PRL 130, 023004 (which is about recoil corrections), so the **content** is sort of okay, but the bibitem title and author list are wrong. Fix by either (a) correcting the bibitem title/authors to "Recoil Corrections..." by Pachucki-Yerokhin, or (b) keeping the title and changing the venue to PRL 121, 073001 (2018). |
| **FriarPayne2005** | **PRA 72, 014501 — Higher-order nuclear-size corrections** | **CITE-WRONG-METADATA** | **MEDIUM SEVERITY.** The actual Friar-Payne 2005 paper on nuclear-size / hyperfine corrections is **Phys. Rev. C 72, 014002** (not PRA 72, 014501). Wrong journal AND wrong article number. The actual title is "Nuclear corrections to hyperfine structure in light hydrogenic atoms" (not "Higher-order nuclear-size corrections in atomic hydrogen"). The cited deuteron Zemach radius r_Z(D) = 2.593(16) fm is consistent with the literature value but its specific attribution to this Friar-Payne paper requires verification. |
| **PachuckiYerokhin2010** | **PRL 104, 070403 — "Theoretical hyperfine splitting in hydrogenic atoms with deuteron nucleus"** | **CITE-MISATTRIBUTED — HIGH** | **HIGH SEVERITY.** PRL 104, 070403 (2010) is actually **"Fine Structure of Heliumlike Ions and Determination of the Fine Structure Constant"** by Pachucki & Yerokhin. This is the Fursaev-Solodukhin / `hep-th/9512134` failure mode: the DOI exists but points to a different paper entirely. I could not find a published 2010 Pachucki-Yerokhin paper specifically titled "Theoretical hyperfine splitting in hydrogenic atoms with deuteron nucleus" or similar. The cited deuteron polarizability contribution +44 ppm needs an actual source. |
| **Eides2024** | **Phys. Lett. B (2024), doi:10.1016/j.physletb.2024.139049** | **CITE-MISATTRIBUTED — HIGH** | **HIGH SEVERITY.** The DOI 10.1016/j.physletb.2024.139049 resolves to "Dynamical non-locality in the near-horizon region of a black hole with quantum time" by Hadi & Akbarieh (PLB 2024). This is NOT an Eides paper on hyperfine splitting. I could not find a published 2024 Eides paper with this exact title. The cited paper for Eides Tab. 7.3 hyperfine-budget reference is almost certainly the 2007 Springer book or the Eides-Grotch-Shelyuto Physics Reports 342, 63 (2001). Fix by replacing this bibitem with the correct primary source. |
| TrackNiZenodo | GeoVac Track NI Zenodo memo | **CITE-OK (internal)** | Internal cross-reference to `track_ni_spectral_triple_zenodo.md`; file present and consistent. |
| SprintHF | GeoVac Sprint HF memos | **CITE-OK (internal)** | Internal sprint provenance. |
| Paper0/2/7/14/17/18/22/32/38/39 | GeoVac internal | **CITE-OK (internal)** | Cross-paper consistency; loading priority per CLAUDE.md §6. Paper 39 marked "(In preparation)" — fine. Paper 2 marked "(Conjectural.)" — preserves required label per CLAUDE.md §13.5. |

### Problems found (CITE-MISATTRIBUTED / DOESNT-SUPPORT / CANT-FIND)

**The three HIGH-severity citation problems** (Pachucki2023, PachuckiYerokhin2010, Eides2024) are the headline finding for this review. All three are **CITE-MISATTRIBUTED** — the DOI / volume / page resolves to a real paper, but it is a different paper than the one cited:

1. `Pachucki2023` — PRL 130, 023004 (2023) is "Recoil Corrections..." by Pachucki-Yerokhin, NOT the "Three-Photon Exchange..." paper by Pachucki-Patkóš-Yerokhin (which is PRL 121, 073001 (2018)).

2. `PachuckiYerokhin2010` — PRL 104, 070403 (2010) is "Fine Structure of Heliumlike Ions..." by Pachucki-Yerokhin, NOT a deuterium hyperfine paper. The deuterium polarizability number +44 ppm needs a different source.

3. `Eides2024` — The DOI points to a black-hole physics paper by Hadi-Akbarieh, NOT an Eides hyperfine paper. The +18 ppm reference budget is almost certainly attributable to Eides-Grotch-Shelyuto 2007 (Springer) or Eides 2001 Phys. Rep. 342, 63.

This is the same failure mode flagged in CLAUDE.md §3 (Fursaev-Solodukhin / hep-th/9512134). A referee can verify a DOI in 30 seconds; finding that three load-bearing recent-physics citations have wrong attributions would be a real embarrassment.

**One MEDIUM-severity problem:**

4. `FriarPayne2005` — wrong journal (PRA → PRC), wrong article number (014501 → 014002), and the title in the bibitem is paraphrased rather than the actual published title. The content cited (r_Z(D) = 2.593 fm) is consistent with the literature.

**One LOW-severity problem:**

5. `Navratil2000` — bibitem title is paraphrased ("Large-Basis Ab Initio NCSM and Its Application to 12C") rather than the actual title ("Properties of 12C in the Ab Initio Nuclear Shell Model"). Authors and venue correct. Cosmetic.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|---|---|---|---|---|
| "We build the first two-species qubit encoding for the deuteron" | abstract | yes (Roggero 2020, Romero 2022, quantum computing nuclear) | no exact prior art for two-species deuteron encoding in GeoVac style; Roggero 2020 did quantum computing for neutrino-nucleus scattering using triton; Romero 2022 did ADAPT-VQE for nuclear structure | Per CONFIDENCE_REVIEW honest-ceiling: cannot confirm "first." Downgrade to "to our knowledge, the first" or simply remove "first." |
| Fock rigidity theorem (S³ unique to Coulomb) | Sec III | already standard SO(4)/Bander-Itzykson | well-established | The novelty is in framing this as a constraint on the GeoVac framework's transfer to nuclei; the theorem itself is not new. Paper labels it correctly as a stated theorem with elementary proof, not a new result. **A**. |
| Track NI "Connes-style real-space multi-particle spectral triple" | Sec VI.E | already cross-referenced to Track NI memo | the Zenodo memo itself acknowledges "does not claim the construction is the first of its kind in the broader chemistry-physics literature" | The paper's claim is properly hedged. **A**. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `Pachucki2023` bibitem: PRL 130, 023004 (2023) is "Recoil Corrections", not "Three-Photon Exchange" | B | CITE-MISATTRIBUTED | **HIGH** |
| `PachuckiYerokhin2010` bibitem: PRL 104, 070403 (2010) is "Fine Structure of Heliumlike Ions", not a deuterium HFS paper | B | CITE-MISATTRIBUTED | **HIGH** |
| `Eides2024` bibitem DOI resolves to a black-hole physics paper, not an Eides paper | B | CITE-MISATTRIBUTED | **HIGH** |
| `FriarPayne2005` bibitem: wrong journal (PRA→PRC) and wrong article number (014501→014002) | B | CITE-WRONG-METADATA | **MEDIUM** |
| Abstract "first two-species qubit encoding for the deuteron" — unverifiable novelty claim | A+B | C (overstatement) | **MEDIUM** |
| Eq (24) deuteron HFS dimensional notation: trailing "Ha" is a frequency-conversion shorthand | A | C | **MEDIUM** |
| Eq (22) Eides Zemach formula: "m_e" is shorthand for m_e c/ℏ inverse Compton wavelength | A | C | **MEDIUM** |
| Sec VI.D acknowledged sign issue in JW off-diagonal encoding: claim of "does not affect" is GEOVAC-only | A | B | **MEDIUM** (already honestly disclosed; consider tightening) |
| `Navratil2000` title slightly paraphrased | B | CITE-WRONG-METADATA | **LOW** |
| Pauli counts 592/712/614 are GEOVAC-only outputs (standard but worth flagging) | A | B | **LOW** |

**Severity totals:**
- HIGH: 3
- MEDIUM: 5
- LOW: 2

**Pass A verdict counts:** A: 12, B: 4 (or part-B), C: 4 (overstatements/conventions), D: 1 (priority claim), E: 0.

**Pass B verdict counts:** CITE-OK: 21, CITE-WRONG-METADATA: 2, CITE-MISATTRIBUTED: 3, CITE-DOESNT-SUPPORT: 0, CITE-CANT-FIND: 0.

## Broadcast readiness: **RED**

The numerical content of Paper 23 is solid — every quantitative claim I recomputed survives, the Roothaan closed form and Taylor expansion are correct to bit-exact precision, the recoil and finite-size and HFS numerics all check out, and the Fock projection rigidity theorem is correctly framed and properly delimits where the framework transfers vs does not. The framework's Pauli counts and 1-norms are GEOVAC-only outputs but internally consistent across Paper 23 / Track NI memo. **Pass A is essentially clean** (A counts dominate, the C-class findings are notational rather than substantive, and the one B-class JW-sign issue is honestly disclosed in §VI.D).

**Pass B is what blocks broadcast.** Three load-bearing citations to recent precision-physics literature (Pachucki2023, PachuckiYerokhin2010, Eides2024) have DOIs / volumes / pages that resolve to completely different papers than the bibitem titles claim. This is the **Fursaev-Solodukhin / `hep-th/9512134` failure mode** explicitly flagged in CLAUDE.md §3 — a referee can verify a DOI in 30 seconds, and any one of these three errors would seriously embarrass the broadcast. The cross-references to deuteron polarizability, three-photon exchange, and the Eides hyperfine budget are precisely the parts a precision-AMO referee would zoom in on, because Paper 23 explicitly positions itself against these references (the +286 ppm chain breakdown in §X is *defined* by these citations). A fourth (FriarPayne2005) has wrong journal/article-number and would also be caught.

This is a citation-hygiene problem, not a physics problem. All four can be fixed by a one-sprint errata pass (correct titles, correct DOIs, correct authors), and Paper 23 then becomes broadcast-ready. Without that fix, broadcast is RED.

## What I could NOT verify (hand to a human expert)

1. The actual numerical values of the Pauli counts (592, 712, 614) and the 1-norms (~227, ~557, ~552 MeV) — these are framework outputs requiring full reproduction of the qubit builder. Inter-corpus consistency with Track NI memo holds.
2. Whether the JW off-diagonal "XY − YX vs XY + YX" sign issue in §VI.D really has *zero* impact on the He-4 antisymmetrized Pauli count (claimed but not externally verified).
3. The Minnesota V_R = 200, V_T = -178, V_S = -91.4 MeV parameter values — these are widely cited in the literature but the primary Thompson-Lemere-Tang 1977 paper is behind a paywall.
4. The deuterium polarizability +44 ppm claim attributed to PachuckiYerokhin2010 — needs an actual source given the cited DOI is wrong.
5. The Eides Tab. 7.3 reference for the +18 ppm hydrogen-21cm budget and the -39.5 ppm Zemach value — almost certainly Eides-Grotch-Shelyuto 2007 (Springer) or Eides 2001 Phys. Rep. 342, 63, but Paper 23's Eides2024 citation is wrong.
6. Priority claim ("first two-species qubit encoding for the deuteron") — per CONFIDENCE_REVIEW honest-ceiling, cannot be positively verified. A nuclear-quantum-computing domain expert could settle this.
