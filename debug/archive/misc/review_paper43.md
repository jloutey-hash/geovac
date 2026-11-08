# Confidence Review: Paper 43 — Lorentzian extension of the four-witness Wick-rotation theorem on truncated $S^3 \times \mathbb{R}$ spectral triples at finite cutoff

## Calibration check
N/A — this is a substantive review, not a calibration pass.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | At every $(\nmax,\Nt)\in\{1,2,3\}\times\{1,11,21\}$ the BW-$\alpha$ period closure $\sigma_{2\pi}^{L,\alpha}(O)=O$ holds bit-exactly | Thm 1.1(i), Thm 3.1 | A | EXTERNAL math (integer spectrum $\Rightarrow$ $e^{2\pi i n}=1$) + MIXED (relies on GeoVac-specific integer spectrum of $K_\alpha^W$ inherited from Paper 42) | Proof reduces to a well-known algebraic identity; the spectrum claim is from Paper 42. Logic verified. |
| 2 | Tomita BW-$\gamma$ period closure $\sigma_{2\pi}^{L,\TT}(a)=a$ bit-exactly | Thm 1.1(ii), Thm 4.1 | A | EXTERNAL (Tomita modular theory $\Delta=(\rho^{-1})^T\otimes\rho$ on $M_n(\mathbb{C})$) + MIXED (integer spectrum of $\KL^{\alpha,W}$) | Eq (28) algebra correct; period closure follows from integer spectrum. |
| 3 | Flow conjugacy $\sigma_t^{L,\TT}(a)=\sigma_{-t}^{L,\alpha}(a)$ bit-exact at general $t$ | Thm 1.1(iii), Thm 5.1 | A | EXTERNAL (Tomita flow $\Delta^{it}\cdot\Delta^{-it}$ vs unitary conjugation $e^{itK}\cdot e^{-itK}$ for $\rho=e^{-K}/Z$) | Proof is one-line algebraic identity; correct. |
| 4 | Six-witness collapse on $\rho_W^L$: exactly zero cross-witness residual | Thm 1.1(iv), Cor 6.1 | A | EXTERNAL (state-dependence of modular operator $\Delta$ on $\rho$ alone) | Correct: under BW choice, $\rho_W^L$ is $\beta$-independent and all witnesses produce the same matrix. |
| 5 | $\Frob{H_{\mathrm{local}}-\DL^W}\big|_{\Nt=1}$ equals Riemannian residual bit-exactly (Theorem 7.1, the headline) | Thm 7.1 | **C** (proof underspecified; result correct) | MIXED — depends on later Thm 7.8 Pythagorean orthogonality | Direct computation at $\nmax=1$: r=2.1332 matches table bit-exact, and the bit-equality $\Frob{H-D}=\Frob{H-iD}$ requires $\mathrm{Tr}(HD)=0$, which is Thm 7.8. The proof as written cites only the trivial $\mathrm{Re}(\mathrm{Tr}(H\cdot iD))=0$ (Lorentzian side), leaving the Riemannian-side cross-term silently dependent on a not-yet-proved theorem. |
| 6 | BBB universal axiom $\{\chi,D\}=0$ fails on truthful $\DGV$; mutual-inconsistency of (chirality-as-$\gamma^5$, chirality-diagonal $\DGV$, BBB axiom) | §3.3.5 + Prop 3.4 + §7.5 | A | EXTERNAL (BBB §5(v)) + GEOVAC-ONLY ($\DGV$ chirality-diagonal by construction) | Directly verified: $[\gamma^5,\DGV]=0$ on the chirality-diagonal construction, so $\{\gamma^5,\DGV\}=2\gamma^5\DGV\ne 0$ unless $\DGV=0$. |
| 7 | Pythagorean orthogonality: $\langle H_{\mathrm{local}},\DL^W\rangle_{\mathrm{HS}}=0$ bit-exact; closed form r² = κ²·S/(4π²) + D | Cor 7.6, Thm 7.8 | A | EXTERNAL (Lemma 7.7 parity argument is standard linear algebra) + MIXED (parity inputs I1, I2 cite Papers 42, 43, 44) | Re-verified Lemma 7.7 computationally with random Π-even A and Π-odd B; Tr(A·B)=0 exactly. Closed form r²=S/(4π²)+D verified at n_max∈{1,2,3,4,5,6} using debug script: r at n_max=1 gives 2.1332 (matches table), n_max=2 gives 6.5275, n_max=3 gives 13.854. PSLQ relation `+1/π² - 4·1/(4π²) = 0` confirmed. |
| 8 | M3 sub-mechanism convention-dependent; sectional in $n_\text{fock}$-parity not signature | §7.10, §7.11 | A | MIXED — internal: depends on Paper 28 vertex-parity convention | Logically clean: the two readings (n_fock-parity vs chirality-pairing) give different verdicts, and the resolution is that M3 in Paper 28's convention is signature-blind. |
| 9 | Riemannian-limit recovery at $\Nt=1$ bit-identical (Frobenius residual exactly 0) | Thm 8.1, Thm 1.1(v) | A | EXTERNAL (Kronecker with $I_1$ identity) | Trivially correct. |
| 10 | "No Lorentzian propinquity construction is published in the spectral-truncation literature as of May 2026" (priority claim) | §9.2, §10.1 (O1) | D | EXTERNAL (literature search) | Cannot confirm novelty; the recently-emerged Nieuviarts 2025/2026 program is correctly placed out-of-scope (odd KO-dim 3); see priority claims table below. |
| 11 | Connes axiom audit at (4,6) — four BBB-predicted-sign relations bit-exact | Thm 3.3 | A | EXTERNAL (Clifford algebra) + MIXED (chiral-basis tensor decomposition of $U_4 = i\gamma^2$) | Proof structure correct. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value/error | Survives? |
|---|---|---|---|---|
| $\Frob{H_L-\DL^W}\big|_{\Nt=1, \nmax=1} = 2.1332$ | 2.1332 | Direct computation (4-dim wedge with $K_\alpha^W = \mathrm{diag}(1,1)$, $D_W^{\mathrm{GV}}=\mathrm{diag}(+3/2,-3/2)$, $H_{\mathrm{local}}=K_\alpha^W/(2\pi)$, $\beta=2\pi$) | r=2.133227740 (matches to ≥10 dps) | YES |
| Closed-form value at $\nmax=1$: $r^2 = S(1)/(4\pi^2)+D(1)$ where $S(1)=2, D(1)=4.5$ | $r=2.1332$ | Closed form $r^2 = 2/(4\pi^2)+4.5 = 4.5507$ | $r=2.13322774$ (matches) | YES |
| Closed-form at $\nmax\in\{1,...,6\}$ | "100 dps PSLQ verified" | Debug script `debug/h_local_residual_closed_form_verify.py` | All match to ≥15 dps | YES |
| Lemma 7.7 (parity ⇒ HS-orthogonality) | $\mathrm{Tr}(A^\dagger B)=0$ | Random Π-even A and Π-odd B (numpy n=6) | 0.0 exactly | YES |

### Circularity map

- **Claim 5 (headline Thm 7.1) — partial GEOVAC-ONLY chain via forward reference.** The proof of Thm 7.1 as written depends silently on Thm 7.8 (Pythagorean orthogonality, proved later in the same section). The result is correct, but the proof is logically downstream of 7.8, not the other way around. This is a presentation defect, not a math error: the bit-exact equality $\Frob{H-D}^2=\Frob{H}^2+\Frob{D}^2$ on the Riemannian side requires $\mathrm{Tr}(HD)=0$, which the proof asserts without justification.

- **Claim 7 (Pythagorean orthogonality) — GEOVAC-ONLY chain.** The parity inputs (I1) $[\Pi_W, H_{\mathrm{local}}]=0$ and (I2) $\{\Pi_W, \DL^W\}=0$ both rely on cited structural inputs from Papers 42, 43, 44 not re-derived here. Honestly flagged in §7.6's "Honest scope (post-formal-proof)" paragraph (caveat i).

- **Claim 11 (Connes axiom audit) — MIXED chain.** External Clifford algebra is fine; the GeoVac-specific tensor decomposition of $U_4 = i\gamma^2$ into chirality-swap × spin-charge-conjugation is a representation choice tied to Sprint L2-B convention.

- **Claim 6 (BBB universal axiom failure) — internal-only consequence.** The three-way mutual inconsistency holds for ANY chirality-diagonal Dirac (not just GeoVac's), so the result is structurally sound; only the path to it is internal.

### Overstatement findings

| Exact phrase | Verdict | Suggested honest replacement |
|---|---|---|
| Abstract: "the geometric BW-$\alpha$ modular Hamiltonian $\KL^\alpha = J_{\mathrm{polar}} \otimes I_{\Nt}$ has integer spectrum and produces the bit-exact period closure" | A | OK as stated |
| Abstract: "$H_{\mathrm{local}} \ne \DL^W$ scope finding [...] is *signature-independent* at the Riemannian reduction $\Nt = 1$" | A | OK — accurately reflects the body (Riemannian limit only) |
| Abstract & intro: "the four-witness theorem [...] has its operator-system analog literally satisfied at signature (3,1) at finite cutoff, not merely as a structural correspondence" | C-soft | "literally satisfied" is borderline strong. The closure depends on the integer spectrum of $\KL^{\alpha,W}$ (a GeoVac feature) and on the choice of *truthful* $\DGV$ (which fails BBB axiom). Suggest: "literally satisfied [...] under the R1 resolution of the BBB universal-axiom finding (truthful $\DGV$)" — already partially flagged in §3.3.5 but not in abstract. |
| §1: "no Lorentzian propinquity construction is published in the spectral-truncation literature as of May 2026" | D (priority claim) | Soft to: "to our knowledge, no Lorentzian propinquity construction is published [...]". §9.2 already does this softening; the abstract / intro should match. |

## Pass B — Citation and novelty

### Citation table

| \cite key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `bisognano_wichmann1976` | BW 1976, JMP 17, 303–321 | CITE-OK | Standard reference, confirmed |
| `bizi_brouder_besnard2018` | BBB 2018, JMP 59, 062303; arXiv:1611.07062 | CITE-OK (with minor singular/plural typo: paper title says "application" but arXiv says "applications") | Confirmed via arXiv |
| `camporesi_higuchi1996` | Camporesi-Higuchi 1996, JGP 20, 1–18 | CITE-OK | Standard reference |
| `chamseddine_connes2010` | Chamseddine-Connes 2010, Fortsch. Phys. 58, 553–600 | CITE-OK | Standard |
| `connes1995` | Connes, "Noncommutative Geometry", Academic Press 1995 | CITE-OK | Standard monograph |
| `connes_rovelli1994` | Connes-Rovelli 1994, CQG 11, 2899–2918; arXiv:gr-qc/9406019 | CITE-OK | Confirmed |
| `connes_vs2021` | Connes-vS 2021, CMP 383, 2021–2067; arXiv:2004.14115 | CITE-OK (minor: paper published 2020, but bibitem says 2021 — within the customary published-date range) | Confirmed |
| `farsi_latremoliere2024` | Farsi-Latrémolière 2024, arXiv:2404.00240 | CITE-OK | Confirmed: "Collapse in Noncommutative Geometry and Spectral Continuity" |
| `farsi_latremoliere2025` | Farsi-Latrémolière 2025, arXiv:2504.11715 | CITE-OK (assumed correct based on earlier corpus checks) | Not independently re-verified in this pass |
| `franco_eckstein2014` | Franco-Eckstein, CQG 30, 135007; arXiv:1212.5171 | CITE-WRONG-METADATA (LOW) | Confirmed authors+venue, but bibitem KEY says "franco_eckstein2014" while published 2013; cosmetic |
| `hartle_hawking1976` | Hartle-Hawking, PRD 13, 2188–2203 | CITE-OK | Standard reference (could not direct-fetch APS DOI due to 403 but standard) |
| `hekkelman2022` | Hekkelman MSc thesis 2022, arXiv:2206.13744 | CITE-OK (assumed; not re-verified) | — |
| **`hekkelman_mcdonald2024`** | **Hekkelman-McDonald, "Spectral truncations of $T^d$ and quantum metric geometry," J. Noncommut. Geom., to appear (2024); arXiv:2403.18619** | **CITE-MISATTRIBUTED (HIGH)** | **arXiv:2403.18619 is "Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures" by Calderón, Rucci, Chichizola — a parallel-computing paper. Eva-Maria Hekkelman's homepage (emhekkelman.nl) lists only TWO arXiv papers: 2412.00628 and 2512.14581. There is no Hekkelman-McDonald paper on T^d.** The actual Riemannian "Gromov–Hausdorff Convergence of Spectral Truncations for Tori" is **Leimbach-van Suijlekom, arXiv:2302.07877**, already cited correctly as `leimbach_vs2024`. The "Hekkelman-McDonald T^d" entry appears to be a fabricated or hallucinated citation. **Propagated to 13 other math.OA papers in the corpus.** |
| `hekkelman_mcdonald2024b` | HM 2024b, J. Funct. Anal. to appear; arXiv:2412.00628 | CITE-OK | Confirmed |
| `latremoliere2018` | Latrémolière, Trans. AMS 370 (2018), 365–411 | CITE-OK (assumed standard) | — |
| `latremoliere_metric_st_2017` | Latrémolière, Adv. Math. 415 (2023); preprint arXiv:1811.10843 | CITE-OK (assumed) | — |
| `leimbach_vs2024` | Leimbach-vS 2024, Adv. Math. 439, Paper No. 109496 | CITE-OK | Confirmed (arXiv:2302.07877) — title is "Gromov–Hausdorff Convergence of Spectral Truncations for Tori", venue and volume match |
| `marcolli_vs2014` | Marcolli-vS 2014, JGP 75, 71–91; arXiv:1301.3480 | CITE-OK (standard) | — |
| `nieuviarts2024` | G. Nieuviarts 2024, arXiv:2402.05839 | CITE-OK (assumed; first-name initial correctly updated to G.) | — |
| `nieuviarts2025a` | G. Nieuviarts 2025, arXiv:2502.18105 v3 | **CITE-WRONG-METADATA (MEDIUM)** | arXiv:2502.18105 actual title is **"Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple"** — the bibitem says "Emergence of pseudo-Riemannian spectral triples within the almost-commutative framework." Title mismatch. The single-author correction (G. not A.) is correct. |
| `nieuviarts2025b_proceedings` | G. Nieuviarts 2025b, arXiv:2512.15450 v2 | **CITE-WRONG-METADATA (MEDIUM)** | arXiv:2512.15450 actual title is **"Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry"** — bibitem title is paraphrased to "Emergence of pseudo-Riemannian structures from twisted spectral triples within the almost-commutative framework (proceedings synthesis)." |
| `sewell1982` | Sewell 1982, Ann. Phys. 141, 201–224 | CITE-OK (standard) | — |
| `strohmaier2006` | Strohmaier 2006, JGP 56, 175–195; arXiv:math-ph/0110001 | **CITE-WRONG-METADATA (LOW)** | actual title is "On Noncommutative and **semi**-Riemannian Geometry" (the published title uses "semi-Riemannian", a synonym for pseudo-Riemannian). The Paper 43 bibitem says "pseudo-Riemannian" — meaning preserved, title cosmetically wrong. |
| `takesaki1970` | Takesaki LNM 128, Springer 1970 | CITE-OK (standard) | — |
| `toyota2023` | M. Toyota 2023, arXiv:2309.13469 | **CITE-WRONG-METADATA (LOW)** | Author's first name is **Ryo**, not M. The bibitem says "M.~Toyota". |
| `unruh1976` | Unruh 1976, PRD 14, 870–892 | CITE-OK (standard) | — |
| `vandungen2016` | van den Dungen 2016, MPAG 19:4; arXiv:1505.01939 | CITE-OK | Confirmed |
| `paper24, paper32, paper38, paper40_unified, paper42, paper45, paper46, paper47` | GeoVac internal preprints | CITE-OK (internal) | All files exist in the corpus |

### Problems found

1. **`hekkelman_mcdonald2024` is CITE-MISATTRIBUTED (HIGH).** The bibitem references arXiv:2403.18619, which is "Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures" (Calderón, Rucci, Chichizola — a parallel-computing paper). There is no Hekkelman-McDonald paper at this ID; Hekkelman's homepage lists only 2412.00628 (already correctly cited as `hekkelman_mcdonald2024b`) and 2512.14581. The actual paper proving GH convergence of spectral truncations of T^d is **Leimbach-vS, arXiv:2302.07877** — already correctly cited as `leimbach_vs2024`. The `hekkelman_mcdonald2024` bibitem appears to be a duplicate / hallucination of the Leimbach-vS torus paper attributed to the wrong authors and with a wrong arXiv ID. Used at lines 1830-1835 (§9.2, "the existing propinquity-convergence literature […] all operates on Riemannian spaces") — a load-bearing citation supporting the priority claim that no Lorentzian propinquity construction exists. **The corpus-wide spread (13 math.OA papers cite this same wrong bibitem) means this is a systematic corpus issue, not paper-43-specific.**

2. **`nieuviarts2025a` and `nieuviarts2025b_proceedings` titles wrong (MEDIUM).** Both bibitems use paraphrased titles that don't match the actual arXiv titles. The dispatcher prompt already mentioned that the first-name initial fix (A. → G.) was applied corpus-wide today; the title-level fix has NOT been applied.

3. **`strohmaier2006` and `toyota2023` minor metadata issues (LOW).** Strohmaier title cosmetic; Toyota first-name initial wrong (M. → R.).

4. **Thm 7.1 proof is logically downstream of Thm 7.8 (MEDIUM presentation defect, not a math error).** The headline H_local-signature-independence proof asserts both $\Frob{H-iD}^2 = \Frob{H}^2+\Frob{D}^2$ AND $\Frob{H-D}^2 = \Frob{H}^2+\Frob{D}^2$, but only justifies the first equality (which is automatic from H real and iD imaginary). The second equality requires $\mathrm{Tr}(HD)=0$, which is the Pythagorean orthogonality proved later in §7.2. Either reorder, or have Thm 7.1's proof cite Thm 7.8 (which would create a citation-out-of-order). Result is correct; presentation defect.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "no Lorentzian propinquity construction is published in the spectral-truncation literature as of May 2026" | §1, §9.2, §10.1 (O1) | Verified the cited Riemannian-side authors (Latrémolière, Hekkelman-McDonald [misattributed; only 2412.00628 exists], Toyota, Farsi-Latrémolière, Leimbach-vS) all work in Riemannian framework | No prior art found for a Lorentzian construction; the cited Nieuviarts 2025/2026 program is correctly placed out-of-scope (odd KO-dim 3) | Already softened to "to our knowledge" in §9.2; abstract and §1 should also soften from absolute to "to our knowledge". §10.1 O1 says "remains open" which is fine. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `hekkelman_mcdonald2024` bibitem points to wrong arXiv ID 2403.18619 (actually parallel-computing paper); the cited paper does not exist as Hekkelman-McDonald work | B | CITE-MISATTRIBUTED / CITE-CANT-FIND | **HIGH** |
| `nieuviarts2025a` bibitem title wrong: actual arXiv title is "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple" | B | CITE-WRONG-METADATA | MEDIUM |
| `nieuviarts2025b_proceedings` bibitem title wrong: actual arXiv title is "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" | B | CITE-WRONG-METADATA | MEDIUM |
| Thm 7.1 proof relies on later Thm 7.8 (Pythagorean orthogonality) without citing it; the Riemannian-side cross-term vanishing is not justified by the stated argument | A | C (proof underspecified; result correct) | MEDIUM |
| "literally satisfied at signature (3,1)" in abstract elides the BBB universal-axiom failure (truthful $\DGV$ caveat) | A | C-soft | MEDIUM |
| `strohmaier2006` bibitem title says "pseudo-Riemannian" but actual arXiv title says "semi-Riemannian" (synonymous, cosmetic) | B | CITE-WRONG-METADATA | LOW |
| `toyota2023` bibitem author initial "M.~Toyota" but actual first name is Ryo | B | CITE-WRONG-METADATA | LOW |
| `bizi_brouder_besnard2018` bibitem title "application" but actual title "applications" (plural) | B | CITE-WRONG-METADATA | LOW |
| `connes_vs2021` bibitem year 2021 vs actual publication 2020 (close enough; CMP volume 383) | B | possibly CITE-WRONG-METADATA | LOW |
| Novelty claim "no Lorentzian propinquity construction is published" in abstract/§1 should soften to "to our knowledge" (already done in §9.2) | A/B | C / priority claim | LOW |

## Broadcast readiness: **YELLOW**

**Synthesis.** Paper 43 is mathematically substantive and the content claims are by-and-large well-supported: the BW-α and BW-γ period closures, the flow conjugacy, the six-witness collapse, and the Riemannian-limit recovery are all bit-exact algebraic identities that I independently re-verified. The Pythagorean orthogonality (Cor 7.6, Thm 7.8) is a genuine new structural result with a correct proof (Lemma 7.7 + parity inputs), and the closed-form residual $r^2 = \kappa_g^2 S(\nmax)/(4\pi^2) + D(\nmax)$ matches the table values to ≥15 dps across $\nmax\in\{1,\ldots,6\}$. The headline Theorem 7.1 (H_local signature-independence at the Riemannian reduction) is correct in result, but the proof as written has a forward-reference gap: it depends on the Pythagorean orthogonality proved later in the same section. This is a MEDIUM presentation defect that should be addressed before broadcast.

The blocking issue is a **CITE-MISATTRIBUTED bibitem** (`hekkelman_mcdonald2024`, arXiv:2403.18619) that points to a parallel-computing paper unrelated to spectral triples; the citation also appears in 12 other math.OA papers in the corpus, indicating this is a systematic corpus problem. Combined with two CITE-WRONG-METADATA Nieuviarts titles (already partially fixed today by the A. → G. corpus pass — title fixes pending) and a handful of LOW-severity cosmetic metadata errors, the citation discipline needs one more pass before broadcast.

The structural content of the paper is sound. The Riemannian-limit recovery falsifier is real (bit-exact agreement with Paper 42 at $\Nt=1$). The cross-references to Papers 42, 44, 45, 46, 47 are consistent with those papers' contents.

## What I could NOT verify (hand to a human expert)

- The exact bookkeeping at finite $\nmax$ for the BBB Connes axiom audit at (4,6) — the four signs ε, ε'', κ, κ'' are stated as coming from BBB Table 1, but I didn't independently retrieve BBB Table 1 to check. A domain expert in indefinite spectral triples should verify Eq (18) and (19) sign assignments.
- Whether the BBB universal axiom violation $\{\chi, D\} \ne 0$ disqualifies the construction from being called a "Krein spectral triple" in any of the field's standard sense (the paper says it "does not constitute a full BBB indefinite spectral triple in the strict §5(v) sense" but accepts the construction; whether this is acceptable is a community-norms question).
- The actual contents of the supporting sprint memos (`debug/l2_*_memo.md`) that the paper cites as the source of the numerical panel; I sampled `debug/h_local_residual_closed_form_verify.py` and it executes cleanly, but the full L2-D / L2-E pipeline was not re-run.
- Whether the temporal-derivative refinement at $\Nt > 1$ (Table 7.4) is the right physical content; the paper interprets it as "the temporal-derivative content $i\gamma^0 \otimes \partial_t$ in $\DL$ that $H_{\mathrm{local}}$ structurally does not capture," but other interpretations (e.g., the residual is an artifact of the centered-FD scheme at finite $T_{\max}$) cannot be ruled out from internal evidence alone.
