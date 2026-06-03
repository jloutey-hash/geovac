# Confidence Review: Paper 27 — Entropy as a Projection Artifact: One-Body Operators are Entanglement-Inert on Sparse Lattices

Reviewer: GeoVac Confidence Reviewer (Wave 1 re-fire, text-level audit of `.tex` source)
Source: `papers/group6_precision_observations/paper_27_entropy_projection.tex`
Date of review: 2026-06-01

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|------|----------|---------|----------|---------------------|
| C1 | "the von Neumann entropy of any single-particle RDM derived from a non-degenerate GS of a purely one-body fermionic Hamiltonian on a finite lattice is identically zero" | abstract, §II prop. | **A** | EXTERNAL (Peschel 2003, Ghirardi–Marinatto 2004 free-fermion / Slater-determinant fact) | Standard result: a single Slater determinant has Gaussian correlation matrix with eigenvalues in {0,1}; the paper correctly cites both. The proposition is well-formulated and the scope (§II.C) explicitly carves out the degenerate / open-shell exceptions with concrete numerical examples (Li 0.637 nats antisymmetrization floor, Be 1.23 nats degenerate-ensemble floor). Honest scope language. |
| C2 | "S_kin/S_full ~ 10^{-14} for He at n_max=2,3" | abstract, Table I | **A** | MIXED (external Peschel fact + GeoVac code) | Recomputed from `debug/data/ep1_entropy_partition.json`: n_max=2 ratio = -5.6e-15, n_max=3 ratio = +1.9e-14. **Bit-exact match to Table I.** |
| C3 | "the area-law A_n ∝ n^4 is recovered as pair counting: A_n = g_n^2 = 4 n^4, factor 4 = two-body signature" | abstract, §III, Eq. (4) | **A** | EXTERNAL (symbolic combinatorics on Paper 0 degeneracies) | Symbolic sympy verification: g_n = 2n^2 ⇒ g_n^2 = 4 n^4 (exact), g_n(g_n-1)/2 = 2n^4 - n^2, g_n g_{n-1} = 4n^4 - 8n^3 + 4n^2 (all three agree at O(n^4)). Cumulative N(n) = n(n+1)(2n+1)/3 ~ (2/3) n^3 with spin — paper formula correct. The "Paper 5 statement was mislabeled in one-particle language" framing is honest and well-presented. |
| C4 | "V_(1s,1s),(1s,1s) = 5/8 exactly; 42% of next-largest; relative commutator 6.1% at n_max=3, 5.3% at n_max=4; F-diagonality 0.920/0.892" | §IV (cusp section) | **A** | GEOVAC-internal numerics, externally reproducible | From `debug/data/energy_graph_nmax3.json`: V_(1s,1s) = 0.625 = 5/8 (bit-exact); rel commutator = 0.06087 (paper 6.1%); V_diag_fraction = 0.9203 (paper 0.920). From `energy_graph_nmax4.json`: rel commutator = 0.0532 (paper 5.3%). All four numerical claims survive. |
| C5 | "HO zero-entropy theorem: S_HO = 0 identically for ANY central V(r12), two-fermion closed-shell GS on Bargmann–Segal lattice" | abstract (4th result), §VI.A | **A** | EXTERNAL (Moshinsky–Talmi total-quanta conservation, standard nuclear-structure fact) + internal verification | Mechanism is sound: M–T transformation makes N_tot = N_rel + N_CM conserved under any central V(r12), so the (0s)^2 GS sits in 1D subspace ⇒ single Slater det ⇒ S=0. Numerical confirmation `debug/data/ep2b_ho_two_fermion.json` shows S_full = 0.0 exactly with occupations (2.0, 0, 0, 0) at both N_max=2 and N_max=3. The structural claim ("for ANY central V") is the right scope: M–T mixes only within the conserved N_tot block, but the GS sits in the lowest N_tot block which is 1D for the closed-shell two-fermion case. |
| C6 | "S_B = A (w̃_B/δ_B)^γ universal scaling; γ_∞ ≈ 1.96 in the joint Z, n_max → ∞ limit" | abstract (5th result), §VI.B, Table III | **C (overstated)** | GEOVAC-internal numerics + external 2nd-order RS expectation (γ=2) | Local-slope numerical claims for individual n_max all survive bit-exact: my recomputed local slopes at Z=15→30 are γ=1.992, 1.976, 1.966, 1.959 for n_max=2,3,4,5 (paper claims identical to 4 digits). EP-2c power-law fit at n_max=3 reproduces γ=2.383, A=8.16, R²=0.998 bit-exactly. **However:** Aitken Δ² extrapolation on (γ_3, γ_4, γ_5) = (1.9755, 1.9662, 1.9589) yields γ_∞ ≈ 1.931; linear-in-1/n extrapolation yields 1.935. **Neither reaches the claimed 1.96.** The paper's "γ_∞ ≈ 1.96" effectively just quotes the n_max=5 value as if it were converged. The text says "Richardson extrapolation on n_max=3,4,5 local slopes" — but Richardson is not actually performed here. Soft overstatement; see overstatement-findings below. |
| C7 | "fits at R²=0.990 with γ=2.228, A=0.157 across 22 calibration points" | §VI.B end of first paragraph | **B** | GEOVAC-only (not re-fit independently here) | The number is plausible given the He-like sub-fit (γ=2.383 at fixed-δ_B=0.20), but I did not re-aggregate the 22-point mixed dataset; rely on `debug/data/ep2g_*.json`. Internally consistent only. |
| C8 | "Be analytical 3×3 degenerate-PT: GS = 0.98|2s²⟩ − 0.14|2p_0²⟩ − 0.14|(2p_-1 2p_+1)_{S=0}⟩, occupations (2.000, 1.924, 0.020, 0.037, 0.020), S_full = 0.794 nats" | §II.C EP-2N paragraph | **A** | derivable from claimed coefficients | Re-derived: 2·0.98² = 1.9208 (paper 1.924, agrees within roundoff if true coefficient is slightly larger than 0.98); 2·0.14² = 0.0392 (paper 0.037 — fine, true coefficient slightly smaller than 0.14); single-orbital S from listed occupations = 0.795 nats (paper 0.794). **S_full claim survives.** S_kin = 1.358 claim from "equal mixture" of 3 degenerate states does not match my naïve recompute (1.229 from listed (2, 0.76, 0.76, 0.47, 0)), but the paper says these are the H_1-only ensemble occupations, not equal-weight; without the explicit ensemble construction I cannot fully reproduce. **Minor: not load-bearing.** |
| C9 | "Sprint TD Track 2: S_thermo(β) = k_B · S_microstate-info(β) at machine precision (residual 4.8×10^{-14}, 90 test points)" | §VI.C | **B** | GEOVAC-only verification of a textbook identity | This identity is essentially the standard statement that the Gibbs ensemble's thermodynamic entropy equals the Shannon entropy of the Boltzmann population over energy eigenstates (the canonical-ensemble identity, e.g., Landau–Lifshitz Vol. 5). The Track 2 sprint verifies it holds for GeoVac eigenstates. The paper's framing ("Sprint TD Track~2 verifies this operationally") is appropriate — does not claim novelty. |
| C10 | "PSLQ tested S_full(GS) against 12,312-form mechanical basis at 150 dps; clean negative; S_full joins K=π(B+F−Δ), L₂ constant c, Wolfenstein parameters on the list of natural framework constants outside the master Mellin ring" | §VI.C | **D** (unverifiable here) | GEOVAC-only PSLQ negative | The claim is internally honest (negative result clearly stated as negative, not novelty). I did not re-run a 12,312-form PSLQ at 150 dps; rely on the linked Sprint TD Track 5 memo. The framing — placing S_full alongside K, c, and Wolfenstein parameters as transcendentals outside the Mellin ring — is consistent with Paper 18 §IV and the memory rule on transcendental tagging. **Not flagged.** |
| C11 | "rigidity statement structurally parallel to the Fock and Bargmann–Segal theorems of Paper 24" (HO zero-entropy theorem) | abstract, §VI.A end | **C (mild)** | MIXED | The parallel is structurally apt at a heuristic level (HO N_tot conservation ↔ Fock S^3 closure ↔ Bargmann–Segal Hardy-sector π-freeness), but "structurally parallel to a theorem" needs more careful language. Paper 24 has named rigidity theorems for two different things (Fock projection rigidity = S^3 unique to −Z/r, HO rigidity = bit-exact π-free in rational arithmetic); Paper 27's HO zero-entropy theorem is genuinely a third rigidity statement of its own, and "structurally parallel to" is a reasonable framing, but the abstract could more precisely state what it shares. **Recommend softening — see below.** |

### Numbers I recomputed

| claim | paper's figure | independent reference (+ source) | my recomputed value/error | survives? |
|-------|----------------|----------------------------------|----------------------------|-----------|
| S_full/S_kin ratio at He n_max=2 | -5.6e-15 | `debug/data/ep1_entropy_partition.json` (regenerated would require pipeline; raw data only) | -5.61e-15 from JSON | **YES** |
| Same at n_max=3 | 1.9e-14 | same | 1.90e-14 | **YES** |
| EP-2c γ He-like fit | γ = 2.383, A = 8.16, R² = 0.998 | my numpy lstsq on the JSON | γ = 2.3830, A = 8.160, R² = 0.9982 | **YES (bit-exact)** |
| EP-2c Z-scaling α_Z | -2.613 | same | -2.6134 | **YES** |
| Local slope γ_loc at n_max=2,3,4,5 (Z=15→30) | 1.992, 1.976, 1.966, 1.959 | `debug/data/ep2L_nmax5_gamma.json` | 1.9921, 1.9755, 1.9662, 1.9589 | **YES** |
| γ_∞ extrapolated | ≈ 1.96 | Richardson / Aitken on n_max=3,4,5 | Aitken: 1.931; linear-in-1/n: 1.935 | **NO — paper's γ_∞ is essentially γ_5, not extrapolated** |
| V_(1s,1s),(1s,1s) at n_max=3 | 5/8 = 0.625 | `debug/data/energy_graph_nmax3.json` | 0.625 (top_hubs[0]/diagonal) | **YES (exact)** |
| Frobenius diag fraction at n_max=3 | 0.920 | same | 0.9203 | **YES** |
| Rel commutator at n_max=3 | 6.1% | same | 0.06087 | **YES** |
| Rel commutator at n_max=4 | 5.3% | `debug/data/energy_graph_nmax4.json` | 0.05319 | **YES** |
| Li open-shell S | 0.637 nats | analytical: -(2/3)ln(2/3) - (1/3)ln(1/3) | 0.6365 | **YES** |
| Be S_full (analytical occupations) | 0.794 nats | analytical from stated occs | 0.795 | **YES** |
| Be coefficient norm² | implied = 1 | 0.98² + 0.14² + 0.14² | 0.9996 | **YES (rounding artifact)** |
| Area-law A_n = 4 n^4 | g_n² = 4n^4 | sympy: (2n²)² | 4 n^4 (exact) | **YES** |
| Unordered pair count | 2 n^4 + O(n²) | sympy: 2n²(2n²-1)/2 | 2 n^4 - n² | **YES** |
| Adjacent-shell pair | 4 n^4 + O(n³) | sympy: 2n²·2(n-1)² | 4n^4 - 8n³ + 4n² | **YES** |

### Circularity map

Most claims rest on a healthy mix of EXTERNAL + GEOVAC-internal:

- **C1, C3, C5** (one-body inert + area-law combinatorics + HO zero entropy) are **EXTERNAL at the math level**: free-fermion Slater-det Gaussian correlation structure (Peschel/Ghirardi–Marinatto), elementary combinatorics on Paper 0 degeneracies, and Moshinsky–Talmi N_tot conservation are all standard. The GeoVac-specific layer is the numerical realization on the S^3 / Bargmann–Segal lattices, which is internal but explicit.
- **C2, C4, C6, C7** (numerical entropy values, EP-2c fits, energy-graph cusp characterization) are **GEOVAC-ONLY** in the sense that they rest on GeoVac's CI pipeline and exact-Fraction Slater integral evaluator. They are NOT externally cross-validated against another quantum-chemistry code on the same basis. However: (i) the EP-1 ratio is the ratio of *two* GeoVac calculations done with V_ee on vs off, so the kinetic-only result S_kin~0 is essentially a self-consistency check of the FCI code against a textbook identity, not an external accuracy claim. (ii) Paper 26's Z^{-2.56} scaling, cited here for cross-consistency with α_Z = -2.613, is itself a GeoVac internal result — so the agreement is a within-corpus self-check, not external validation.
- **C9** (apparatus identity at machine precision) is a self-check of standard quantum statistical mechanics applied to GeoVac eigenstates; the identity itself is EXTERNAL, the verification is GEOVAC-only.
- **C10** (PSLQ negative on S_full) is GEOVAC-only, but is honestly framed as a negative result — no novelty / discovery claim is made.

**The only GEOVAC-ONLY claim that drives a substantive structural conclusion is C6's γ_∞ ≈ 1.96 extrapolation, which (a) does not actually survive my independent Richardson and (b) is connected to a 2nd-order RS expectation (γ=2) that is correct EXTERNAL physics. The interpretive question is whether γ → 2 (from above) or to a number slightly below 2 (the paper's claim) — see overstatement.**

### Overstatement findings

1. **"γ_∞ ≈ 1.96 in the joint Z, n_max → ∞ limit (Richardson on n_max=2,3,4,5)"** (abstract; §VI.B Table III caption "Richardson extrapolation on n_max=3,4,5 gives an asymptotic slope γ_∞ ≈ 1.96"). Richardson / Aitken Δ² on the actually computed sequence (1.9755, 1.9662, 1.9589) gives 1.93–1.94, not 1.96. The 1.96 figure is *essentially the n_max=5 value* (1.9589 rounded up). The directional claim — "below the non-degenerate second-order Rayleigh–Schrödinger value of 2" — is **supported by the data**; the magnitude of "below" is mis-stated.
   - **Suggested replacement:** "The local exponent decreases monotonically with n_max from 1.992 (n_max=2) to 1.959 (n_max=5); the sequence is below the non-degenerate second-order Rayleigh–Schrödinger value of 2 and is decreasing further with n_max, consistent with γ_∞ < 2. A 1/n_max-linear extrapolation on the n_max=3,4,5 sub-sequence yields γ_∞ ≈ 1.93. The residual gap is attributed to multi-shell aggregation in the basis-extension limit." This is honest about the extrapolation method and the actual extrapolated value.
   - **Severity: MEDIUM.** The structural conclusion ("γ_∞ < 2 by a basis-aggregation residue") is correct and load-bearing for the §VI.B reframing. The numerical figure is what's off.

2. **"this is a rigidity statement structurally parallel to the Fock and Bargmann–Segal theorems of Paper~24"** (abstract). Paper 24's two named theorems are the Fock projection rigidity theorem (S^3 unique to -Z/r) and the HO rigidity theorem (π-free graph certificate). The HO zero-entropy theorem in Paper 27 §VI.A is genuinely a third rigidity result; "structurally parallel to" is at most a heuristic statement. The actual structural mechanism (Moshinsky–Talmi N_tot conservation) is closer to a *consequence* of the Paper 24 HO π-free skeleton than a "parallel" to it.
   - **Suggested replacement:** "a rigidity statement on the Bargmann–Segal lattice that extends the Paper 24 HO π-free skeleton (rational matrix elements, Moshinsky–Talmi total-quanta conservation) from the operator level to the ground-state entanglement spectrum." Honest about what the parallel actually is.
   - **Severity: LOW.**

3. **"the local exponent extrapolates to γ_∞ ≈ 1.96 in the joint Z, n_max → ∞ limit"** (abstract first paragraph). Same overstatement as Finding 1, repeated in the abstract. Same suggested replacement.
   - **Severity: MEDIUM** (it's in the abstract — first thing a reader sees).

4. **"All ground-state entanglement is generated by the two-body interaction V_ee, with the non-degeneracy qualifier essential: open-shell and dissociated-bond states escape the floor"** (abstract). This is accurate — the open-shell qualifier IS in the abstract — but the qualifier "non-degenerate" is sandwiched in a way that the casual reader might miss it. The §II.C scope discussion is **excellent** (Li, Be, He carved out concretely); the abstract underplays how serious the qualifier is. Minor framing issue, not a real overstatement.
   - **Suggested replacement:** small reorder of the abstract to lead the qualifier — "of any non-degenerate closed-shell ground state."
   - **Severity: LOW.**

5. **"five results are presented"** (abstract opening). The five are: (i) one-body entanglement-inert (proven externally), (ii) area-law as pair counting (combinatorial), (iii) cusp hot-node localization (numerical observation), (iv) HO zero-entropy (theorem), (v) universal scaling exponent (numerical fit with extrapolated γ_∞). The mix is uneven — three are theorem-grade, two are numerical fits. The opening sentence treats all five at the same conceptual weight, which is generous to (iii) and (v).
   - **Suggested replacement:** "Three theorem-grade results and two numerical-fit observations are presented" — or absorb (iii) and (v) into the discussion without elevating to abstract-headline status. **Optional/minor.**

---

## Pass B — Citation and novelty

### Citation table

| \cite key | Claimed as | Verdict | What I found (URL / DOI) |
|-----------|------------|---------|--------------------------|
| Peschel2003 | "Calculation of reduced density matrices from correlation functions, J. Phys. A 36, L205 (2003)" — used to support the free-fermion / Slater-determinant correlation-matrix fact | **CITE-OK** | DOI 10.1088/0305-4470/36/14/101 confirmed; abstract matches use in §II.B and the §Equation-Verification paragraph closing. |
| LiHaldane2008 | "Entanglement spectrum as a generalization of the entanglement entropy, Phys. Rev. Lett. 101, 010504 (2008)" — cited in §I as folklore in quantum-info literature for the kinetic/potential asymmetry | **CITE-OK** | PRL 101, 010504 (2008) confirmed; this is the Li–Haldane entanglement-spectrum paper on Moore–Read ν=5/2; Paper 27 cites it for entanglement-spectrum folklore rather than for the specific topological content. The citation is in-scope. |
| GhirardiMarinatto2004 | "General criterion for the entanglement of two indistinguishable particles, PRA 70, 012109 (2004)" — supports the antisymmetrization (Slater) entropy framing in §II.B | **CITE-OK** | PRA 70, 012109 (2004) confirmed via arXiv quant-ph/0401065. The paper's content (Slater–Schmidt number, vN entropy of 1-particle RDO to detect entanglement vs antisymmetrization) directly supports Paper 27's use. |
| Moshinsky1996 | "M. Moshinsky and Yu. F. Smirnov, The Harmonic Oscillator in Modern Physics (Harwood Academic, 1996)" — supports the Moshinsky–Talmi total-quanta conservation result | **CITE-OK** | Confirmed: Harwood Academic, 1996, ISBN 9783718606214. Moshinsky–Talmi transformation is the canonical reference for N_tot conservation under central two-body potentials in the 3D HO basis. |
| GeoVac_Paper0, 5, 7, 14, 18, 23, 24, 26 | Internal GeoVac papers | **CITE-OK** (existence) | All exist in the repo at `papers/group{1..6}/...`. The §I-§VI uses are appropriate (Paper 5 explicitly flagged "conjectural" in the bibliography, Paper 26's Z^{-2.56} cross-referenced as Z-scaling consistency). |
| EP1_data, EP2_data, EP3_data, EnergyGraph | GeoVac project data files | **CITE-OK** | All files confirmed present in `debug/data/`. Note: EP3 cited via `debug/data/ep3_area_law_analysis.md` exists as `.md` not `.json` — fine, the cite explicitly says "memo". |

### Problems found

**None.** Every citation checks out. No CITE-MISATTRIBUTED, CITE-DOESNT-SUPPORT, CITE-CANT-FIND, or CITE-WRONG-METADATA findings on this paper. The reference list is unusually clean for a GeoVac paper — 4 external + 11 internal, all correct.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|------------------|----------|----------|------------------|----------------|
| "its sharp form is rarely written out, and its consequences on sparse, lattice-realized Hamiltonians are seldom exhibited explicitly" | §I, second paragraph | not formally searched (paper does NOT claim novelty here, only "rarely written out") | n/a | No action — paper does not assert "first." The framing is appropriately humble. |
| The HO zero-entropy theorem (§VI.A) is framed as: "is now established... strengthened to a structural statement by Moshinsky–Talmi total-quanta conservation" — paper presents it as a new GeoVac result | §VI.A | Conceptually: closed-shell two-fermion HO GS is well-known in nuclear structure literature to be (0s)² ⇒ single SD. The structural step "⇒ S=0 for any central V via M–T conservation" is straightforward once stated. | The full statement as a *theorem* may be novel to the GeoVac series, but "two HO closed-shell electrons stay in a single Slater determinant under any central interaction" is essentially folklore in the harmonic-oscillator / Talmi-bracket nuclear-structure community. | Paper does NOT claim "first in the literature." Framing is acceptable. **No edit needed.** A domain expert would likely recognize the result as a clean reformulation of a folklore fact, not a discovery. |
| "γ_∞ ≈ 1.96, slightly below the non-degenerate second-order Rayleigh–Schrödinger value of 2" (with "residual gap attributed to multi-shell aggregation") | abstract, §VI.B | This is a fit / extrapolation claim, not a novelty claim. | n/a | See Pass A overstatement-Findings 1+3. |

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---------|------|---------|----------|
| 1 | γ_∞ ≈ 1.96 "Richardson on n_max=2,3,4,5" — actual extrapolation gives ~1.93 | A | C (overstated) | **MEDIUM** |
| 2 | Same overstatement, repeated in abstract | A | C (overstated) | **MEDIUM** |
| 3 | "structurally parallel to the Fock and Bargmann–Segal theorems of Paper 24" — heuristic framing | A | C (mild) | **LOW** |
| 4 | Abstract underplays non-degeneracy qualifier (§II.C scope is good, abstract less so) | A | C (mild) | **LOW** |
| 5 | "Five results" treats theorem-grade and numerical-fit results at uniform weight | A | C (very mild / stylistic) | **LOW** |

No HIGH severity findings. No CITE-* problems. No math errors. No CITE-MISATTRIBUTED.

---

## A–E counts (Pass A)

- **A (externally verified):** C1, C2, C3, C4, C5, C8 = **6**
- **B (internally consistent only):** C7, C9 = **2**
- **C (overstated):** C6, C11 = **2** (with the abstract-level repetition of C6, severity-counted separately)
- **D (unverifiable here):** C10 = **1**
- **E (wrong):** **0**

## Citation verdict counts (Pass B)

- **CITE-OK:** 15
- **CITE-WRONG-METADATA / MISATTRIBUTED / DOESNT-SUPPORT / CANT-FIND:** 0

## Severity totals

- **HIGH:** 0
- **MEDIUM:** 2 (γ_∞ figure, repeated)
- **LOW:** 3 (parallel-to-Paper-24 framing, abstract qualifier ordering, "five results" weighting)

---

## Broadcast readiness: **GREEN-leaning YELLOW**

This is one of the cleaner GeoVac papers I have audited at the text-source level. The proposition's mathematical content (C1) is correctly attributed to external free-fermion literature; the area-law combinatorial reframe (C3) is symbolically exact; the HO zero-entropy theorem (C5) is structurally sound via Moshinsky–Talmi and numerically clean. **All four external citations check out without error**, and the EP-1 / energy-graph numerical claims (C2, C4) reproduce bit-exactly from the data files.

The single substantive issue is the γ_∞ ≈ 1.96 extrapolation claim (C6) which is repeated in the abstract: the local slopes are accurate, but "Richardson extrapolation" is not actually performed — the 1.96 figure is effectively the n_max=5 value rounded, while an honest 1/n_max extrapolation yields γ_∞ ≈ 1.93. The downstream structural conclusion (γ_∞ < 2, with residue from multi-shell aggregation) is **supported by the data either way**, but the numerical figure should be softened to either (a) "the n_max=5 local slope of 1.959, decreasing further with n_max" or (b) an honest extrapolation yielding ~1.93. This is a MEDIUM-severity errata-batch fix, not a broadcast blocker.

The "structurally parallel to Paper 24 theorems" framing (Finding 3) and the abstract's qualifier ordering (Finding 4) are LOW-severity polish.

**Recommendation:** YELLOW with two MEDIUM fixes (Findings 1+2, both about γ_∞) batched into an errata sprint; otherwise broadcast-ready. After the γ_∞ fix this paper would be a confident GREEN.

---

## What I could NOT verify (hand to a human expert)

1. **γ_∞ ≈ 1.96 extrapolation method.** My independent Richardson/Aitken extrapolation on the published (γ_3, γ_4, γ_5) sequence gives ~1.93. The paper claims 1.96 from "Richardson on n_max=3,4,5" but does not state the Richardson formula used. A domain expert / numerical-analysis specialist should review whether the paper's extrapolation is in fact a *different* extrapolation scheme than the standard Aitken Δ² or 1/n-linear (perhaps a 1/n² scheme, or weighted by something), and either justify the 1.96 figure or accept the 1.93 figure.

2. **C7's R² = 0.990 over a 22-point mixed dataset (γ = 2.228, A = 0.157).** Independent re-fitting on the combined `ep2c + ep2k LiH-bond R-sweep` dataset would be needed to confirm. Within scope of an errata-sprint diagnostic if the paper is going to broadcast.

3. **The "operational temperature decoding" §VI.C identity** is standard quantum statistical mechanics. The 12,312-form PSLQ negative result at 150 dps was not independently re-run. The framing (S_full belongs on the "outside-Mellin-ring" list with K, c, Wolfenstein) is consistent with CLAUDE.md §1.7 WH4 / Paper 18 §IV but only a specialist in PSLQ-on-spectral-zeta would recognize whether the 12,312-form basis is genuinely exhaustive of the master Mellin engine + adjacent rings.

4. **Whether the HO zero-entropy theorem is genuinely a new result.** The paper does not claim "first" — but a nuclear-structure expert familiar with the Moshinsky–Talmi literature might recognize it as folklore. Not a defect; just a flag.

---

## Summary line for dispatcher

**Paper 27 — "Entropy as a Projection Artifact: One-Body Operators are Entanglement-Inert on Sparse Lattices" — YELLOW (one MEDIUM, three LOW); A=6, B=2, C=2, D=1, E=0; CITE-OK=15, other CITE-*=0; HIGH=0, MEDIUM=2, LOW=3. Top finding: the abstract's "γ_∞ ≈ 1.96, Richardson on n_max=2,3,4,5" claim is the n_max=5 local slope (1.9589) rounded, not an actual extrapolation — an honest 1/n_max extrapolation on (γ_3, γ_4, γ_5) gives γ_∞ ≈ 1.93; downstream conclusion (γ_∞ < 2 by basis-aggregation) is supported either way, only the numerical figure needs softening.**
