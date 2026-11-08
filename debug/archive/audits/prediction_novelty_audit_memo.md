# Prediction novelty audit: ten pre-registration records (2026-10-19)

**Question.** For every scored prediction in ten records under `docs/predictions/` (the four nuclear-door records excluded), was the outcome already known in the literature before we computed it? The real hit rate is the number of confirmed predictions whose outcome was not already known.

**Method.** I read each record in full, with its disposition, and the owning paper sections: Paper 59 `sec:literature`, `sec:four_centre`, `sec:three_electron`, `sec:four_density` and `sec:density_ring`; Paper 7 `sec:continuum`; Paper 24 `sec:large_spin` and `sec:ks_frame`. I checked load-bearing sources in their primary text where I could get it: arXiv PDFs and the Seaton review, text-extracted locally. "Verified" below means I read the stated location myself. "Via citing work" means I saw it only in a work that cites it, or in a search summary.

**Search depth.** The searches were targeted, not exhaustive. Novelty claims for the period-ladder items rest mainly on Paper 59's own targeted search (`sec:literature`) plus one further web search of my own on three-electron/two-centre Slater integrals with K3, Calabi-Yau or elliptic geometry. That search found only false friends: elliptical *coordinates*, and series or quadrature treatments. I did no independent search for the Struve-moment operator of the triangle. Nothing under `debug/ising_class_lit_downloads/` was used.

**Classes.**
- FINDING: held, and not found in the literature.
- CHECK: held, and already known.
- MISS: failed.
- UNCLEAR: held only as evidence, or the outcome is ambiguous.
- UNSCORED: not tested, or stopped.

FINDINGs carry a weight tag. *Routine* means the outcome is guaranteed by elementary calculus on a new object. *Substantive* means it could have come out otherwise.

## Table

### 2026-10-13 four-centre genus (Paper 59 `sec:four_centre`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| P1 | A four-centre (XY\|ZW) integral keeps only two dispersion roots, so its genus is at most 1 | HELD | **PARTLY KNOWN.** The momentum-space form is the Averys' (*Hyperspherical Harmonics*, App. B; the record checked the book, I did not). Paper 59 found no genus or elliptic reading of it anywhere (`sec:literature`). The corpus credits the Averys. | FINDING (low weight: the paper itself calls it "a corollary ... a calibration rather than a hard test") |
| P2 | Geometry enters only through the coefficients W, D1, D2 | HELD | NOT FOUND as a statement. It is immediate from the same reduction. | FINDING (routine) |
| P3 | Genus 2 first appears in three-electron integrals (hunch) | not tested; superseded by P3a–P3e | — | UNSCORED |

### 2026-10-13 three-electron K3 (Paper 59 `sec:three_electron`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| P3a | The angle integral is elementary, and its only roots are the three Δ's | HELD (stronger than predicted: elementary for every D) | NOT FOUND for this object. It is a standard angular reduction of Yukawa-type kernels. | FINDING (routine) |
| P3b | The fibre is a double cover of P² branched on six lines (a K3), with non-isotrivial genus-1 slices | HELD | **PARTLY KNOWN.** That a six-line double cover is K3 is classical (Matsumoto–Sasaki–Yoshida, as the record notes; not re-read). The three-loop banana is K3 (Bloch–Kerr–Vanhove; Paper 59 cites it). Neither Paper 59's search nor mine found this Slater chain placed on a K3. | FINDING (substantive, moderate) |
| P3c | The one-centre limit is elementary, no K3 | HELD (rational) | **PARTLY KNOWN / routine.** With one-centre 1s densities the potentials are elementary, so the integral is elementary by inspection. Closed forms for one-centre three-electron Slater integrals exist (Fromm–Hill 1987; Remiddi 1991; neither read). | CHECK |
| P3d | The holomorphic period is a sum of products of elliptic periods ("elliptic squared") | HELD | **KNOWN.** Bailey–Borwein–Broadhurst–Glasser, J. Phys. A 41, 205203 (2008), arXiv:0801.0891, eq. (33): (2c/π) V̄₃ = K(k₋)K(k₊′) + K(k₊)K(k₋′), from W. N. Bailey's second paper on infinite Bessel integrals. **Verified in the primary text.** Credited in Paper 59 ("the closed form is that of Ref. bbbg2008"). The reduction from the chain to V̄₃ is the corpus's own. | CHECK |
| P3e-1 | At equal scales the true fibre collapses onto K3 periods plus Tate numbers | HELD (216 G = 80 c₃₀ − 6600 c₃₂ + 81π³, stable at 40 and 70 digits) | NOT FOUND. BBBG-type recursions reduce equal-mass moments ∫xᵏK₀ᵃK₁ᵇ to c₃₀, c₃₂, so the mechanism is standard. The object, which carries the iterated piece ∫K₀, and the identity are not found. | FINDING (moderate-low) |
| P3e-2 | At generic scales the fibre is not in the K3 + residue-elliptic + Tate span | Registered test INCONCLUSIVE (cap too low). Holds as evidence under a powered rerun amended after the control and before G was examined. | NOT FOUND. Filed by `/walls` as an open-leaning §3 row. | UNCLEAR (a null PSLQ at one mass triple, with a basis known to be incomplete) |

### 2026-10-16 linking ceiling (Paper 59 `sec:four_density`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| PL1 | Tree identity: P_n = 2ⁿ/(π√∏c) ∫∏K₀(a_v x)dx, for any tree shape | HELD | **PARTLY KNOWN.** This is the standard position-space (Fourier) representation of banana-type integrals: a tree maps unimodularly onto the hyperplane ΣL = 0. The record itself calls it "a one-line Fourier argument". No source states it for trees of linked densities. | FINDING (routine) |
| PL2 | F₄ is not built from elliptic pieces (no quadric, order 4) | HELD | **KNOWN.** The operator is entry #58 of Almkvist–van Enckevort–van Straten–Zudilin (math/0507430), a Hadamard product (Almkvist–Zudilin, math/0402386), with the same holomorphic coefficients. An order-4 Calabi–Yau operator with a non-constant Yukawa coupling has Sp₄ monodromy, so the outcome follows. Broadhurst (2009, "Bessel moments, random walks and Calabi–Yau equations") ties W_N moments to Calabi–Yau equations. Found and credited by the corpus after the result (Paper 59 [PRIOR ART]). The tables were not re-read by me. | CHECK |
| PL3 | Consequence: the ceiling is set by the method, not the molecule | not separately tested | — | UNSCORED |
| PL4 | F₅ is not built from lower rungs (hunch) | PARTLY (one quadric, orthogonal rank 5; not elliptic-built) | **KNOWN.** The F₅ operator is printed with AESZ entry #39 (exterior-square type), per Paper 59's own [PRIOR ART]. The record states the "open question ... is answered by the literature". | CHECK |

### 2026-10-17 density loop (Paper 59 `sec:density_ring`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| PT1 | The loop integrates out by the star-triangle relation, with constant π³ | HELD | **KNOWN.** The uniqueness relation goes back to D'Eramo–Peliti–Parisi and Symanzik; the general-D constant is collected by Preti 2018. Credited in Paper 59. Not re-read by me. | CHECK |
| PT2 (main) | The Coulomb triangle does not climb past the chain | HELD, but in a fourth, unregistered way: it lands *below* the chain, on the elliptic rung (L_G = M(θ − λ), M ≅ Legendre) | NOT FOUND. Paper 59 says no literature search was made for the operator `eq:ring_op`. The record's check of Avery's book found nothing. I did not search the Struve-moment form. | FINDING (substantive, moderate) |
| PT2 (side) | Singular points lie in {0, ±2, ∞} | FAILED (t = ±1 genuine) | — | MISS |
| PT3 | The one-centre triangle is not rational at (1,2,3) and (2,1,3) | HELD (dilogarithms) | **KNOWN.** Remiddi, Phys. Rev. A 44, 5492 (1991): closed form in logs and dilogs. Reproduced to 2e-21 via the PPR 2004 form. Credited. The record noted, before computing, that differentiating in the exponents could lower the weight, so this was a genuine (if likely) prediction. Its outcome is nevertheless contained in Remiddi's formula. | CHECK |
| PT4 | Consequence: a ring of three sits with the chain of three | scored "not supported" (the ring sits lower) | — | MISS |

### 2026-10-17 Ising class (Paper 59 `sec:four_density`, last paragraphs)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| IS1 | The singular t² of F_n's operator are exactly the (n−1)-banana thresholds (n−1−2j)², so t = 1 is singular iff n is even (n = 6, 7, 8 predicted) | HELD (operators checked on F_n to 1e-42) | **PARTLY KNOWN.** That equal-mass banana singular points are the threshold and pseudo-thresholds (m ± m ± … ± m)² is standard. Pögel–Wang–Weinzierl, arXiv:2212.08908, and Bönisch et al., arXiv:2008.10574, state it; I saw this only via search summaries. Combined with the Almkvist–Zudilin Hadamard structure, the n ≤ 8 singular set and the parity follow. The explicit n = 6–8 operators and the exponent pattern are not found. | CHECK (outcome foreseeable from known banana thresholds) |
| IS2 | Odd moments stay on their rung: c₃,₀ is elliptic and not Tate; c₅,₀ and c₇,₀ have no lower relation | HELD within bounds | **KNOWN / consistent with the known.** BBBG: c₃,₀ = 3Γ(1/3)⁶/(32π2^{2/3}) (verified, BBBG l.660ff of the extracted text). Lim–Tu–Yu (arXiv:2609.13576; verified to exist and read) place the n = 5, 7 moments inside 2×2 determinants of periods, so no single lower-basis relation is expected. | CHECK |
| IS3 | Even moments drop: c₄,₀ is a weight-4 L-value and not Tate; c₆,₀ drops too | HELD for the registered value cases | **KNOWN.** Lim–Tu–Yu eq. (43): ∫K₀⁴ = 4π² L(f_8.4.a.a, 2) (credited there to BBBG (47) and Zagier/AMNT). Eq. (44): ∫K₀⁶ ≟ 108π² L(f_24.6.a.a, 4), numerically. **Both verified in the primary text.** Credited in Paper 59. | CHECK |
| IS3' | c₄,₀ is attached to the conifold fibre at x = 1/4 | not tested | — | UNSCORED |

### 2026-10-17 ring extensions (Paper 59 `sec:density_ring`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| QG1 | The 3-D one-loop triangle with masses (γ, γ, 0) is elementary | HELD (solid angle, arctan form) | **KNOWN.** Davydychev & Delbourgo, J. Math. Phys. 39, 4299 (1998), hep-th/9709216, §5.1, eqs. (5.1)–(5.5): the three-dimensional three-point function with arbitrary masses is the solid angle Ω⁽³⁾ (a spherical-triangle area) over √D⁽³⁾, in arccos/arcsin form. They state that it corresponds to Nickel, J. Math. Phys. 19, 542 (1978). **D–D verified in the primary text; Nickel only via D–D.** **Not credited by the corpus** (see credit problems). | CHECK |
| QG2 | The geminal triangle stays on the elliptic rung | FAILED (K3 at γ = 1) | — | MISS |
| QR1 | A ring of four sits on the K3 rung | STOP for feasibility | — | UNSCORED |
| QR2 | The one-centre square is not rational | not reported | — | UNSCORED |

### 2026-10-18 hydrogen continuum (Paper 7 `sec:continuum`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| HC1 | n → iη; Casimir → −(1+η²); zonal-function continuation; l, m survive | HELD | **KNOWN.** Bander & Itzykson 1966. Credited ([PRIOR ART] tag). | CHECK |
| HC2 | The contact density is Z/π on the bound side and (Z/π)/(1 − e^{−2πη}) in the continuum | HELD | **KNOWN.** The Sommerfeld factor (Bethe–Salpeter). Credited. | CHECK |
| HC3 (main) | One universal factor: df/dE = Q_i(iη)/(1 − e^{−2πη}) for every initial state and every l | HELD | **KNOWN.** Seaton, Rep. Prog. Phys. 46, 167 (1983): eq. (2.31) A(ε,l) = ∏_{p=0}^{l}(1 + p²ε); §2.7.3, eq. (2.51), B(ε,l) = A(ε,l)[1 − exp(−2π/k)]⁻¹ for ε > 0, with the energy-normalised s, c built from B^{1/2} times functions f, g that are analytic in ε (§2.5.3); §2.11, the bound s_n normalised per unit energy (∫s_n² = n³/2). Together these say the continuum energy density is the analytic continuation of the bound density per unit energy times 1/(1 − e^{−2π/k}), l-independent. Seaton attributes the mathematical theory to Ham (1955), and bound–free continuation to Burgess–Seaton (1960). **Seaton verified in the primary text** (some equation bodies were lost in text extraction; eq. (2.51) and the definitions are legible). Stobbe's formula is credited; **the one-factor statement is not.** | CHECK |
| HC4 | The price: rational after stripping; no Coulomb phase; e^{πη} is the only new class | HELD | **KNOWN** in substance: Stobbe/Bethe–Salpeter squared quantities contain no arg Γ; it follows from HC3. The tag against Papers 18/34 is an internal classification, not a literature question. | CHECK |

### 2026-10-18 large spin (Paper 24 `sec:large_spin`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| LS1 | Contraction matrix element √(n+1)√(1 − n/2S); the commutator equals 1 − n/S | HELD | **KNOWN.** Holstein–Primakoff 1940; Inönü–Wigner 1953. Credited. | CHECK |
| LS2 | Label data algebraic, with rational squares, at finite S | HELD | **KNOWN** (trivial: Wigner d at π/2 is built from square roots of binomials). | CHECK |
| LS3 | p_S(x) = π^{−1/2}e^{−x²}[1 + c₁/S + …] with c₁ = −1/8 + x²/2 − x⁴/6; no half-integer powers | HELD | **KNOWN** (the Stirling/Edgeworth correction to de Moivre–Laplace for the symmetric binomial; the record derived c₁ by hand from Stirling). Not located in a specific primary text. The unregistered c₂ was measured, not predicted. | CHECK |
| LS4 | Eigenvector components tend to Hermite functions at rate O(1/S) | HELD | **KNOWN.** Krawtchouk → Hermite (Koekoek–Lesky–Swarttouw). Credited. | CHECK |
| LS5 | The sphere's coherent-state measure contracts to the Bargmann measure; π enters only in the limit | HELD | **KNOWN** (Radcliffe 1971; Arecchi et al. 1972; Wallis). Credited. The M1 tag is a reading. | CHECK |

### 2026-10-18 zero energy (Paper 7 `sec:continuum`, "Zero energy")

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| ZE1 | At E = 0 the Runge–Lenz algebra contracts to E(3), with A² = Z² | HELD | **KNOWN.** Immediate from Pauli's 1926 relations at H = 0 (cited). That E = 0 is the third Bander–Itzykson case is noted in Kibler, quant-ph/0503039, fn. 15 ("reduces to a free-particle system for E = 0"; verified). | CHECK |
| ZE2 | n³R_nl²/Z² → 2J_{2l+1}(√(8Zr))²/r from both sides, in whole powers of E | HELD | **KNOWN.** Seaton 1983 §2.5.3 and §2.7.2: f is analytic in ε with convergent ε-expansions, and "Ham (1955) shows that the functions f₍ₙ₎ ... can be expressed in terms of Bessel functions" (verified in Seaton; Ham not read). Whole powers of E with the same coefficient on both sides is exactly that analyticity. | CHECK |
| ZE3 | The threshold function has rational Taylor coefficients and no π | HELD | **KNOWN** (trivial from the 0F1 series). | CHECK |

### 2026-10-19 Kustaanheimo–Stiefel lift (Paper 24 `sec:ks_frame`)

| ID | Statement | Outcome | Literature | Class |
|:--|:--|:--|:--|:--|
| KS1 | The Q = 0 sector has dimension n², with l = 0..n−1 once each | HELD | **KNOWN.** Hydrogen as a constrained 4-D oscillator: Kustaanheimo–Stiefel 1965; Boiteux 1973; Kibler–Négadi, Lett. Nuovo Cim. 37, 225 (1983), J. Phys. A 16, 4265 (1983), PRA 29, 2891 (1984). Seen via Kibler, quant-ph/0503039 §3.3 (verified there; the primaries not read). Paper 24 states "Part (i) is classical". | CHECK |
| KS2 | The Casimir 2(J₊² + J₋²) = n² − 1 on Q = 0 | HELD | **KNOWN**, as KS1 (the constraint selects j₁ = j₂ of SU(2)×SU(2)). | CHECK |
| KS3 | The pushed-down states are Coulomb Sturmians with Z = n p₀ | HELD | **KNOWN.** Coupling-constant metamorphosis (Hietarinta et al. 1984) and Coulomb Sturmians (Rotenberg). Credited. | CHECK |
| KS4 | The R⁴/R³ inner-product constant is π/4, one fibre volume | HELD | **KNOWN** (the KS measure relation d⁴u ∝ r⁻¹ d³x dφ; KS 1965, Stiefel–Scheifele; not re-read). | CHECK |
| KS5 | π² appears at even argument in all four spectral-zeta ladders | HELD | **KNOWN / routine** (partial fractions and Hurwitz-zeta identities). Its value was as an internal correction: it withdrew Paper 24's calibration-π layer. That is an atlas correction, not a literature finding. | CHECK |

## Totals

| Class | Count | Items |
|:--|:-:|:--|
| FINDING | 7 | P1, P2, P3a, P3b, P3e-1, PL1, PT2 (of which routine: P2, P3a, PL1; low weight: P1) |
| CHECK | 27 | P3c, P3d, PL2, PL4, PT1, PT3, IS1, IS2, IS3, QG1, HC1–HC4, LS1–LS5, ZE1–ZE3, KS1–KS5 |
| MISS | 3 | PT2-side, PT4, QG2 |
| UNCLEAR | 1 | P3e-2 |
| UNSCORED | 5 | P3, PL3, IS3', QR1, QR2 |

**Real hit rate.** 34 of the 41 scored items held (FINDING + CHECK); 3 failed; 1 is unclear. Of the 34 that held, 27 were already known. **Only 7 confirmed predictions have outcomes not found in the literature, and only 3 of those (P3b, PT2, P3e-1) are substantive**; P1 is a corollary-level calibration. All four door records (continuum, large spin, zero energy, KS) are 100% CHECK, and they say so themselves in their "honest caps". That is as designed: they are accounting exercises on textbook mathematics. All FINDINGs come from the Paper 59 period-ladder records.

**A common dependency.** Every FINDING inherits Paper 59's base novelty claim, that no one has placed a molecular Slater integral on the elliptic/Bessel-moment Feynman machinery (`sec:literature`, itself "an absence from a targeted search"). If an obscure source has made that bridge, P1, P2, PL1 and probably P3a fall to CHECK together. P3b, PT2 and P3e-1 would likely survive as statements about specific new objects. All FINDINGs are also fibre (D = 0) statements.

## FINDINGs ranked by confidence of novelty

1. **PT2: the Coulomb triangle at the two-centre fibre lands on the elliptic rung, below the chain.** The Struve-moment form `eq:ring_struve`, the operator `eq:ring_op` with L_G = M(θ − λ_G) and M ≅ Legendre, and the drop below the three-density chain. The outcome was not anticipated even by the registration. No search was made for the operator (Paper 59 says so), so the confidence is moderate-high that it is unpublished, given how specific the object is. The ingredients are classical (star-triangle; Legendre).
2. **P3b: the three-electron two-centre chain fibre is a K3** (six lines in P², no four concurrent, non-isotrivial slices). The geometry type is classical; the placement of a chemistry integral on it was not found. Moderate.
3. **P3e-1: the equal-scale collapse 216 G = 80 c₃₀ − 6600 c₃₂ + 81π³.** A specific identity for a new object, by a standard IBP mechanism. Moderate-low as a contribution, but likely new.
4. **P1: no genus 2 at four centres.** New only as an extension of Paper 59's own reading; the momentum form is the Averys'. Low weight.
5. **PL1, P2, P3a: routine.** They are new only because the objects are new.

## Credit problems: the paper home presents a CHECK as its own

1. **Paper 7 `sec:continuum`, eq. `continuum_factor`, and the "Zero energy" paragraph (HC3, ZE2). Clearest case.** The [PRIOR ART] tag covers only the n → iη continuation. The sentence "What the continuum adds beyond the continued bound states is one factor", the l-independence argument, and "approached in whole powers of E from each side, with the same coefficient on both" carry no attribution. All three are quantum-defect theory: Seaton 1983 eqs. (2.31), (2.51) and §2.5.3/§2.11, with Ham (1955) and Burgess & Seaton (1960) for the bound–free continuation. The record's own disclosure ("continuity ... across threshold (Seaton's quantum-defect theory)") shows the source was known before registration, but it did not reach the paper. Suggested fix: add a [PRIOR ART] sentence citing Seaton 1983 §2.7.3/§2.11 and Ham 1955, and keep "the accounting on the seam" as the contribution.
2. **Paper 59 `sec:density_ring` (QG1).** The paper says: "No literature search was made for the solid-angle form of the massive loop; it rests on its numerical check." That form is Davydychev & Delbourgo, J. Math. Phys. 39, 4299 (1998), §5.1 (eqs. 5.1–5.5), who attribute the three-dimensional triangle to Nickel (1978). The ring-extensions record marks Nickel "NOT VERIFIED"; D–D confirm it as the source. Fix: cite both and drop the no-search sentence.
3. **Paper 59 `sec:four_density`, minor (IS1).** "The threshold reading is what the ladder adds." The singular set is the standard threshold/pseudo-threshold locus of the equal-mass banana, (Σ±m)². The paper names the points "thresholds" but cites no source for the locus, and the even-n parity is immediate from it. Suggest citing a banana-threshold source (e.g. Bönisch et al. 2021, already cited elsewhere in the paper as `bonisch2021`; or Pögel–Wang–Weinzierl), and narrowing the novelty to the explicit n = 6–8 operators and the exponent pattern.
4. **Paper 59 `eq:tree`, optional (PL1).** It is tagged [SYMBOLIC + MEASURED] with no note that it is the standard x-space representation of banana integrals. The text does later relate F_n to the banana; a half-sentence would close it.

No credit problem found in Paper 24: both `sec:large_spin` and `sec:ks_frame` say "classical" and cite the sources. Kibler–Négadi or Boiteux could be added for the label-set and Casimir statements (KS1, KS2), which are cited there only as "classical". The ising-class, linking-ceiling and density-loop records match their paper homes on credit.

## Sources read

Primary text: Seaton 1983 (§§2.5–2.11); BBBG arXiv:0801.0891 (eq. 33, c₃,₀); Lim–Tu–Yu arXiv:2609.13576 (eqs. 43–44); Davydychev–Delbourgo hep-th/9709216 (§5.1); Broadhurst arXiv:1604.03057 (§7.3: level 6 belongs to the n = 6 sunrise problem, not c₆,₀); Kibler quant-ph/0503039 (§3.3, fn. 15). Search summaries only: arXiv:2212.08908, arXiv:2008.10574. Not re-read: AESZ, Remiddi, Avery, Ham, Kibler–Négadi, Nickel, Holstein–Primakoff, Koekoek et al.
