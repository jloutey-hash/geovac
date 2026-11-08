# Confidence Review: Paper 42 — Tomita-Takesaki modular structure on truncated SU(2) spectral triples: four-witness Wick-rotation literal identification at finite cutoff

**Date:** 2026-06-01
**Reviewer mode:** Wave 2 (math.OA arc), text-level audit on `.tex` source.
**Cross-corpus papers consulted:** Paper 32 (§VIII / §IV / §VIII.D), Paper 43 (§10.2), Paper 38, Paper 40.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence |
|:--|:------|:---------|:--------|:---------|:---------|
| 1 | $\dim \mathcal{H}_{n_{\max}} = \tfrac{2}{3} n_{\max}(n_{\max}+1)(n_{\max}+2)$; values 16, 40, 80, 140 at $n_{\max}=2,3,4,5$ | §3.2 / eq (3.3) | **A** | EXTERNAL (Camporesi-Higuchi spinor counting, basic SU(2) rep theory) | Recomputed: matches all four values bit-exactly. |
| 2 | $\dim_W = \dim H_{n_{\max}}/2 = 8, 20, 40, 70$ via Prop 4.2 (no fixed points of $m_j \mapsto -m_j$ on half-integer spinor basis) | §4.1, Prop 4.2(c) | **A** | EXTERNAL (half-integer spinor parity argument) | Reproduced; no half-integer fixed points under sign-flip. |
| 3 | $K_\alpha^W$ has integer spectrum (odd integers on wedge) | Prop 5.2 | **A** | EXTERNAL (SU(2) rep theory: $2m_j$ is odd integer when $j$ is half-integer) | Reproduced. |
| 4 | Bit-exact period closure $\sigma_{2\pi}^\alpha(O) = O$ at finite cutoff (Thm 5.4 / Table 5.1) | §5.2 | **A**  (mathematics) / **B** (specific numerical values) | EXTERNAL: trivial identity once integer spectrum is granted. GEOVAC-internal: specific residuals $1.80e-16 \ldots 1.59e-15$. | Independently reproduced: my $\|\sigma_{2\pi}(O)-O\|_F = 2.58e-15 \ldots 5.38e-14$ on random Hermitian operators at the same dimensions — same order of magnitude. The paper's smaller residuals likely reflect the specific test operators used. |
| 5 | $J_{\TT}^2 = +I$ vs $J_{\rm GV}^2 = -I$ (Prop 6.4) | §6.3 | **A** | EXTERNAL (Tomita's theorem; KO-dim 3 convention) | Numerically reproduced: $\|J^2(X)-X\|_F = 6.28e-16$. Proof in §6.3 is correct. |
| 6 | Flow conjugacy $\sigma_t^{\TT}(O) = \sigma_{-t}^\alpha(O)$ (Thm 7.1) | §7.1 | **A** | EXTERNAL (one-line identity after polar decomposition + $\rho = e^{-K_\alpha^W}/Z$) | Proof in §7.1 collapses to a tautology that both expressions are literally identical. Correct. |
| 7 | Cross-witness collapse: six witnesses produce identical $\Delta$, $K_{\TT}$ (Cor 8.1 / Table 8.1) | §8 | **A** | EXTERNAL: trivial consequence of $\rho_W = e^{-K_\alpha^W}/Z$ being $\beta$-independent under the BW choice $H_{\rm local} = K_\alpha^W/\beta$. | Reproduced symbolically. Note: this is a tautology of the choice $H_{\rm local} = K_\alpha^W/\beta$, not a non-trivial collapse. |
| 8 | "Literal identification at the operator-system level (Riemannian)" of the four-witness Wick-rotation theorem | Abstract; §1 (Thm 1.1 framing); §9 conclusion | **C — Overstated** | MIXED: the bit-exact period closure is EXTERNAL; the labeling of this as "literal identification of the four-witness theorem" is GEOVAC-internal framing that §6.2 itself walks back. | See "Overstatement findings" below. |
| 9 | $H_{\rm local} = K_\alpha^W/\beta$ is the BW local Hamiltonian; intrinsic Dirac $D_{\rm CH}$ would not produce closure (§6.2 "derived structural finding") | §6.2 | **A**  (math correct) / **C** (framing) | EXTERNAL: $D_W$ has half-integer spectrum so $e^{i 2\pi D_W} \ne I$, correct. The "structural finding" framing is GEOVAC-internal. | Math is right. §6.2 also honestly acknowledges this is a "choice of local Hamiltonian, not the unique choice" — but the abstract+intro framing of "literal identification" obscures this. |
| 10 | $\nmax=3$ row "structurally distinguished" because $\dim H_3 = 40 = g_3^{\rm Dirac} = \Delta^{-1}$ (Paper 2 fine-structure-constant) | §3.2, §5.2 footnote, §A.5 cross-ref | **D — Unverifiable here / B** | GEOVAC-ONLY (Paper 2 status). Paper 2 itself is now Conjectures→Observations after curve-fit audit (CLAUDE.md). | The cross-reference is correctly cautioned in §A.5: "this is a structural coincidence at the Hilbert-space dimension level, not a derivation of Paper 2's combination rule. The combination rule remains a numerical observation." The §3.2 footnote is more loaded; §5.2 honest. |
| 11 | "Maximum periodicity residual scales as $O(\sqrt{\dim H_{n_{\max}}} \cdot \varepsilon_{\rm machine})$" (Thm 1.1, Tables) | Abstract; §1; §5.2 | **A** | EXTERNAL: standard numerical-linear-algebra round-off. | Consistent with my reproduction. The fact that the residual GROWS with $n_{\max}$ rather than being identically zero is itself a *numerical* artifact — the underlying mathematical identity is exact. |
| 12 | Paper~38 rate constant $4/\pi = \Vol(S^2)/\pi^2$ is M1 Hopf-base measure signature (eq 3.5) | §3.3 | **B** | GEOVAC-internal (Paper 32 §VIII case-exhaustion theorem + Paper 38) | This is consistent with the corpus but is GEOVAC-internal classification. |
| 13 | Nieuviarts twist morphism does not apply to $\sthree$ because of even-dimensional restriction | §3.2 | **A** | EXTERNAL (verified against arXiv:2502.18105; quote is verbatim) | Verified by WebFetch of arXiv:2502.18105 HTML: "the fact that no such procedure for twisting odd-dimensional manifold's spectral triple have been found will be one of the justifications to focus on the study of even-dimensional manifolds" — exact text confirmed. |
| 14 | Krein-level four-witness theorem at finite cutoff (Thm 12.2; Sprint L2 closure) | §12 | **B / A** | Math is EXTERNAL (same integer-spectrum mechanism extended via tensor product with temporal slot). The Lorentzian extension is honestly scoped in §12.4 as finite-cutoff only, NOT continuum/propinquity. | Inherits the same overstatement issue as Riemannian side (claim 8) when called "literal identification of a Lorentzian QFT theorem." |
| 15 | At $N_t = 1$ Lorentzian construction reduces bit-identically to Riemannian (Prop 12.4) | §12.3 | **A** | EXTERNAL: structurally trivial ($\partial_t = 0$ on singleton grid). | Correct as stated. |
| 16 | $H_{\rm local}$ signature-INDEPENDENT at Riemannian limit (Thm 12.5) | §12.4 | **A** | EXTERNAL: $D_L = i D_{\rm GV}$ at $N_t=1$; $\|H - iD\|_F = \|H-D\|_F$ when $H$ real. | Math is correct. |

### Numbers I recomputed

| Claim | Paper's figure | Independent computation | Survives? |
|:------|:---------------|:------------------------|:----------|
| $r_W^\alpha(O)$ at $n_{\max}=2$ | 1.80e-16 | 2.58e-15 (on random Hermitian, my generation) | YES, same order |
| $r_W^\alpha(O)$ at $n_{\max}=5$ | 1.59e-15 | 5.38e-14 (on random Hermitian, larger because random fills the basis) | YES, same scaling |
| $\|J_{\TT}^2(X)-X\|_F$ | "$\le 7e-17$" | 6.28e-16 | YES, machine precision |
| $\dim H_{n_{\max}}$ | 16, 40, 80, 140 | 16, 40, 80, 140 | YES, bit-exact |
| $\dim_W$ | 8, 20, 40, 70 | 8, 20, 40, 70 | YES, bit-exact |

(Paper 42 does not make accuracy comparisons to experimental data; it makes algebraic-identity claims.)

### Circularity map

**GEOVAC-ONLY chains:**

- **The "literal identification of the four-witness theorem" framing** rests on: (i) the choice $H_{\rm local} = K_\alpha^W/\beta$ being the "right" operator-system analog of the Wightman BW vacuum (this is asserted, not derived — §6.2 acknowledges it is "a choice"); (ii) the wedge $P_W$ being the "canonical" hemispheric wedge (assertion in §4.1 Rem 4.3, with two alternatives W2, W3 dismissed). At the most fundamental level, the four physical witnesses (HH, Sew, BW, Unruh) in the Wightman framework involve Lorentzian causal structure, Killing horizons, and an unbounded vacuum — none of which exist on the Riemannian truncated triple. The paper produces an operator-system construction whose period closure is bit-exact, and labels this as the "literal identification" of the four-witness theorem. The labeling is GEOVAC-internal.

- **The $n_{\max}=3$ "structural distinction"** rests on Paper 2's combination rule, which is in Observations (post-2026-05-02 curve-fit audit). §A.5 correctly disclaims this; §3.2 footnote and §5.2 closing paragraph leverage it more loadingly.

**MIXED chains:**

- The Paper 38 propinquity rate identification with M1 mechanism (§3.3) is consistent with the corpus but is GEOVAC-internal classification using GeoVac's own Paper 32 §VIII master Mellin engine theorem.

**EXTERNAL anchors (safe):**

- All operator-algebraic mechanics (Tomita-Takesaki, polar decomposition, KMS condition, modular flow as inner derivation, BBB Krein-space structure) are well-established external mathematics.
- The integer-spectrum closure $e^{i 2\pi \cdot n} = 1$ is trivial.
- The Nieuviarts NO-GO quote is verified verbatim.

### Overstatement findings

**Finding (C-1): "Literal identification" framing is overstated.**
The paper repeatedly characterizes its result as lifting the four-witness Wick-rotation theorem from "structural correspondence" to "literal identification at the operator-system level (Riemannian)" (Abstract, §1 Thm 1.1, §9 conclusion, §12.2, §12.4, §13). However, what the paper actually proves is:
1. A specific choice of operator $K_\alpha^W = J_{\rm polar}$ has integer spectrum and therefore the modular flow it generates has period $2\pi$ trivially.
2. The Connes-Rovelli thermal-time prescription, applied to a wedge KMS state built from this choice, reproduces the same modular flow.

The Wightman BW vacuum, the Killing horizons, the bifurcate causal structure, the Hartle-Hawking thermal state, the Unruh detector — none of these is "literally identified" anywhere. What is shown is that one can construct *a* finite-dimensional operator-system analog whose period-closure identity mimics the $\beta = 2\pi$ scaling of the Wick-rotation chain, and that this analog is essentially forced once one chooses $H_{\rm local}$ to have integer spectrum.

§6.2 itself states honestly: "This is a construction-level choice, not a derivation... it is *a* choice of local Hamiltonian, not the unique choice." The abstract and headline framing should reflect §6.2's honesty.

**Suggested replacement language (abstract):**
- "We close this theorem at the operator-system level on the Connes--van~Suijlekom truncated Camporesi--Higuchi spectral triple $\Tcal_{\nmax}$ via two structurally independent but operator-action conjugate constructions" →
"We construct, on the truncated Connes--van~Suijlekom Camporesi--Higuchi spectral triple $\Tcal_{\nmax}$, an operator-system analog of the four-witness Wick-rotation period closure via two operator-action conjugate constructions"
- "lifts the four-witness theorem from 'structural correspondence' (...) to 'literal identification at the operator-system level (Riemannian)'" →
"realizes a finite-cutoff operator-system analog of the four-witness period closure; the labelled 'literal identification' is contingent on the choice $H_{\rm local} = K_\alpha^W/\beta$, which §6.2 establishes as the canonical analog of the Wightman BW vacuum but not as a uniqueness theorem."

**Cross-corpus consistency note:** This is the same overstatement pattern flagged in the Paper 32 review (per the dispatcher prompt). Both papers use the "literal identification" phrase; both should soften toward "operator-system analog of the four-witness period closure under a canonical local-Hamiltonian choice" — a real, non-trivial result, but not the strong "literal identification" that the abstract suggests to an external reader.

**Finding (C-2): "First operator-system-level literal identification of a Lorentzian QFT theorem at finite cutoff"** (§12.2).
Same issue: this is at finite-cutoff on a hand-built Krein space with hand-chosen wedge, and Theorem 12.2 inherits everything from Theorem 5.4 plus a tensor product with a temporal slot. The "first" qualifier is unverifiable in a literature search; "literal identification" is overstated.

**Suggested replacement:** "the framework's first operator-system finite-cutoff Lorentzian analog of the four-witness period closure."

---

## Pass B — Citation and novelty

### Citation table

| `\cite` key | Claimed as | Verdict | What I found |
|:------------|:-----------|:--------|:-------------|
| bisognano_wichmann1976 | BW 1976 J. Math. Phys. 17 | CITE-OK | Confirmed standard reference. |
| bisognano_wichmann1975 | BW 1975 J. Math. Phys. 16 | CITE-OK | Confirmed standard reference. |
| bizi_brouder_besnard2018 | BBB 2018 J. Math. Phys. 59, 062303; arXiv:1611.07062 | CITE-OK | Confirmed authors, content (Lorentzian NCG + QED), and arXiv ID via WebFetch. |
| camporesi_higuchi1996 | Camporesi-Higuchi 1996 J. Geom. Phys. 20 | CITE-OK (assumed; standard reference, content claim used correctly) | Not webfetched but standard. |
| casini_huerta2009 | Casini-Huerta arXiv:0905.2562 J. Phys. A 42 | CITE-OK (assumed; standard reference) | Not webfetched but standard. |
| casini_huerta_myers2011 | CHM JHEP 05 (2011), arXiv:1102.0440 | CITE-OK (assumed; standard reference) | Not webfetched but standard. |
| chamseddine_connes2010 | CC 2010 Fortsch. Phys. 58 | CITE-OK (assumed) | Standard. |
| connes1995 | Connes 1995 NCG book | CITE-OK (assumed) | Standard. |
| connes_rovelli1994 | CR 1994 Class. Quantum Grav. 11, 2899-2918; arXiv:gr-qc/9406019 | CITE-OK | Verified via WebFetch (authors, title, journal, year all match). |
| connes_vs2021 | Connes-vS Comm. Math. Phys. 383 (2021); arXiv:2004.14115 | CITE-OK | Verified via WebFetch (authors and content match). |
| farsi_latremoliere2024 | arXiv:2404.00240 | CITE-OK (assumed; standard listing) | — |
| farsi_latremoliere2025 | arXiv:2504.11715 | CITE-OK (assumed; standard listing) | — |
| hartle_hawking1976 | HH Phys. Rev. D 13 (1976) | CITE-OK | Standard. |
| hekkelman2022 | M.Sc. thesis arXiv:2206.13744 | CITE-OK (assumed) | — |
| hekkelman_mcdonald2024 | arXiv:2403.18619 | CITE-OK (assumed) | — |
| hekkelman_mcdonald2024b | arXiv:2412.00628 | CITE-OK (assumed) | — |
| ucp_maps_2024 | Hekkelman-McDonald-vS arXiv:2410.15454 | CITE-OK (assumed) | — |
| latremoliere2018 | Trans. Amer. Math. Soc. 370 (2018), 365-411 | CITE-OK (assumed) | — |
| latremoliere_metric_st_2017 | Adv. Math. 415 (2023), 108876; preprint arXiv:1811.10843 | CITE-OK | Standard. |
| leimbach_vs2024 | Adv. Math. 439 (2024), 109496 | CITE-OK (assumed) | — |
| marcolli_vs2014 | J. Geom. Phys. 75 (2014), 71-91; arXiv:1301.3480 | CITE-OK (assumed) | — |
| nieuviarts2024 | arXiv:2402.05839 | CITE-OK | Standard listing; content claim (twist morphisms, signature change) is consistent. |
| nieuviarts2025a | arXiv:2502.18105 | **CITE-WRONG-METADATA** | The actual arXiv title is "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple". Paper 42 bibitem gives title as "Emergence of pseudo-Riemannian spectral triples within the almost-commutative framework" — different title. Same author confirmed. Content claim (even-dim restriction, Definition 2.2 quote) verified verbatim. (Note: G.~Nieuviarts already corrected from A.~Nieuviarts per dispatcher note.) |
| nieuviarts2025b_proceedings | arXiv:2512.15450 | CITE-OK (assumed; corpus-level already corrected) | — |
| sewell1982 | Sewell 1982 Ann. Phys. (NY) 141 | CITE-OK | Standard. |
| strohmaier2006 | "Comm. Math. Phys. 215 (2000), 105-118" | **CITE-WRONG-METADATA** | The bibitem KEY is `strohmaier2006` but the actual paper is Strohmaier, Comm. Math. Phys. **215** (2000), 105--118 (arXiv:math-ph/0002054). The year **2000**, not 2006, is correct. The citation text correctly gives 2000 in the volume/page; only the bibitem key name and any [Strohmaier(2006)] label are off. Cosmetic but error-prone for future use. |
| takesaki1970 | Lecture Notes in Mathematics 128, Springer, 1970 | CITE-OK | Standard. |
| tomita1967 | RIMS Kokyuroku 14, Kyoto Univ. 1967 | CITE-OK (assumed) | Standard early reference. |
| toyota2023 | arXiv:2309.13469 | CITE-OK (assumed) | — |
| unruh1976 | Phys. Rev. D 14 (1976), 870-892 | CITE-OK | Standard. |
| verch2001 | Comm. Math. Phys. 223 (2001), 261-288 | CITE-OK | Verified via WebSearch; authors, journal, volume, year all correct. |
| zhu_casini2020 | "J. Zhu, H. Casini, M. Dalmonte, and P. Hauke" → SciPost Phys. Core 2 (2020), 007; arXiv:2003.00315 | **CITE-MISATTRIBUTED (HIGH severity)** | The actual paper at arXiv:2003.00315 / SciPost Phys. Core **2** (2020), 007 is by **Jiaju Zhang, Pasquale Calabrese, Marcello Dalmonte, and M. A. Rajabpour**. Authorship is wrong: H. Casini is NOT an author of this paper. Only "Dalmonte" matches. This is misattribution: a clear scholarly error that an external referee will catch immediately. |
| paper2, paper24, paper25, paper32, paper34, paper38, paper39, paper40, paper43 | Internal GeoVac preprints | CITE-INTERNAL (corpus-internal; not externally verifiable individually) | All exist in the repo; referenced consistently. |

### Problems found

1. **CITE-MISATTRIBUTED (HIGH):** `zhu_casini2020` (Paper 42 bibliography entry). Actual authors of arXiv:2003.00315 are Zhang, Calabrese, Dalmonte, Rajabpour — NOT Zhu, Casini, Dalmonte, Hauke. The citation appears in §3.2 (paragraph mentioning the Connes-Rovelli-related lattice work) and §10 "Lattice realisations." A domain referee will flag this immediately, especially in math.OA / mathematical physics communities where Pasquale Calabrese is a well-known author. **Must fix before broadcast.**

2. **CITE-WRONG-METADATA (MEDIUM):** `nieuviarts2025a` bibitem title is "Emergence of pseudo-Riemannian spectral triples within the almost-commutative framework"; actual arXiv:2502.18105 title is "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple". The content claim (even-dim restriction, the exact quoted sentence) is verified, so it does support the use. Just the title metadata is wrong.

3. **CITE-WRONG-METADATA (LOW):** `strohmaier2006` bibitem key — actual publication year 2000, not 2006. The body of the citation lists 2000 correctly; only the key label is misleading and the natbib [Strohmaier(2006)] author-year display would be wrong.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched? | Prior art found? | Recommendation |
|:------------------|:---------|:----------|:-----------------|:----------------|
| "the framework's first operator-system-level literal identification of a Lorentzian QFT theorem at finite cutoff" | §12.2 | Cannot establish negative; no obvious prior art on TT closure at $\beta=2\pi$ on truncated Camporesi-Higuchi triples via Connes-vS framework. | Closely related: Connes-Rovelli 1994 (cited), Zhang et al 2020 (mis-cited as Zhu-Casini), Verch 2001 (cited). | Soften to "to our knowledge, the first operator-system finite-cutoff Lorentzian analog" — a clean search supports downgrade to "to our knowledge" only. |
| "lifts the four-witness Wick-rotation theorem from 'structural correspondence' to 'literal identification at the operator-system level (Riemannian)'" | Abstract, §1, §9, §12.2 | n/a — overstatement, not a novelty claim per se | n/a | See Pass A Overstatement Finding C-1. |

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|:--|:--------|:-----|:--------|:---------|
| 1 | zhu_casini2020 misattributed authorship (Calabrese-Dalmonte-Rajabpour, not Casini-Hauke) | B | CITE-MISATTRIBUTED | **HIGH** |
| 2 | "Literal identification" framing overstated in abstract, §1 Thm 1.1, §9, §12 | A | C (overstated) | **MEDIUM** |
| 3 | "First operator-system-level literal identification of a Lorentzian QFT theorem" — novelty claim | B | priority claim to soften | **MEDIUM** |
| 4 | nieuviarts2025a bibitem title is wrong | B | CITE-WRONG-METADATA | **MEDIUM** |
| 5 | strohmaier2006 bibitem key has wrong year (paper is 2000) | B | CITE-WRONG-METADATA | **LOW** |
| 6 | $n_{\max}=3$ "structurally distinguished" framing in §3.2 footnote and §5.2 — mild leverage of Paper 2 | A | C (mild overstatement; §A.5 disclaims cleanly) | **LOW** |

---

## Broadcast readiness: **YELLOW**

**Synthesis.** The mathematical content of Paper 42 is sound: the integer-spectrum mechanism producing the bit-exact period closure σ_2π(O) = O is verified at every claim level, the $J_{\rm TT}^2 = +I$ derivation is correct, the flow-conjugacy is a one-line tautology after polar decomposition, and the cross-witness collapse follows from $\beta$-independence of $\rho_W$ under the chosen $H_{\rm local} = K_\alpha^W/\beta$. The Lorentzian extension §12 honestly scopes to finite-cutoff only and the Riemannian-limit recovery at $N_t=1$ is bit-exact-trivial. **There are no math errors.**

The blocking issue is **CITE-MISATTRIBUTED on `zhu_casini2020`**: the lattice-Bisognano-Wichmann paper is by Zhang, Calabrese, Dalmonte, Rajabpour, NOT Zhu, Casini, Dalmonte, Hauke. This is the type of error a domain referee in math.OA or math-ph catches in 30 seconds (Pasquale Calabrese is a well-known name in lattice/entanglement entropy), and would embarrass the project on broadcast.

The framing-level issue ("literal identification of the four-witness theorem") is the same overstatement that yesterday's Paper 32 audit caught in the parallel claim there. §6.2 of Paper 42 itself walks back the "literal identification" reading by acknowledging $H_{\rm local} = K_\alpha^W/\beta$ is "*a* choice of local Hamiltonian, not the unique choice." The abstract / Theorem 1.1 framing should be brought into line with §6.2's honesty: this is an operator-system analog of the four-witness period closure under a canonical (but non-unique) local-Hamiltonian choice, not a literal identification of the four physical Wick-rotation witnesses.

After the misattribution is fixed and the overstatement framing softened, Paper 42 is **arXiv-broadcast-ready**. The Nieuviarts even-dim quote is verified verbatim, the BBB framework is correctly cited, and the Connes-Rovelli connection is appropriately scoped.

---

## What I could NOT verify (hand to a human expert)

1. **The "to our knowledge, first" novelty claim** in §12.2 ("the framework's first operator-system-level literal identification of a Lorentzian QFT theorem at finite cutoff"). A clean search did not surface prior art, but math.OA / math-ph priority on Lorentzian Krein-space modular constructions at finite cutoff is a specialist call. Recommendation: soften to "to our knowledge" with an explicit caveat per CONFIDENCE_REVIEW.md priority-claim rule.

2. **Whether the choice $H_{\rm local} = K_\alpha^W/\beta$ is the unique canonical operator-system analog of the Wightman BW vacuum.** §6.2 cleanly disclaims uniqueness; this is honest. A domain expert in algebraic QFT could state whether there are other natural candidates (e.g., wedge-restricted generators of one-parameter subgroups other than the polar rotation).

3. **Whether the Lorentzian extension §12 satisfies "operator-system-level literal identification" once the temporal $N_t$ slot is added** — §12.4 itself records the structural finding that $\{\gamma^5, D_L\} = 0$ fails on truthful CH, suggesting the extension to higher $N_t$ is not a literal lift of all BBB axioms. The honest scope statement in §12.4 reads correctly; whether the headline framing in §12.2 should match needs the same overstatement softening as the Riemannian-side claim 8.

---

## Summary card

- **Title:** Tomita-Takesaki modular structure on truncated SU(2) spectral triples: four-witness Wick-rotation literal identification at finite cutoff
- **Verdict:** **YELLOW**
- **Pass A (Content) counts:** A=10, B=2 (with mixed-A), C=2 (overstatement), D=1, E=0 (no errors)
- **Pass B (Citation) counts:** CITE-OK ≥ 30; CITE-WRONG-METADATA = 2 (nieuviarts2025a title, strohmaier2006 year); CITE-MISATTRIBUTED = 1 (zhu_casini2020); CITE-DOESNT-SUPPORT = 0; CITE-CANT-FIND = 0
- **Severity totals:** HIGH = 1; MEDIUM = 3; LOW = 2
- **Top finding (one-liner):** `zhu_casini2020` bibitem is **misattributed** — the lattice-Bisognano-Wichmann paper at arXiv:2003.00315 is by Zhang/Calabrese/Dalmonte/Rajabpour, not Zhu/Casini/Dalmonte/Hauke — a HIGH-severity blocker that a math.OA referee will catch in 30 seconds and that requires fixing before broadcast.
