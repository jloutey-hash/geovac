# Confidence Audit: Paper 32 — The GeoVac Spectral Triple: A Synthesis Paper Locking the Framework into the Marcolli–van Suijlekom Gauge-Network Lineage

**Auditor:** External skeptical referee (CONFIDENCE_AUDITOR persona)
**Date:** 2026-06-01
**Scope:** Wave 1 text-level audit, cross-corpus with Papers 25, 28, 30, 31, 38

---

## Calibration check

Not a calibration paper.

---

## Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | A finite spectral triple $(\mathcal{A}_{\mathrm{GV}}, \mathcal{H}_{\mathrm{GV}}, D_{\mathrm{GV}})$ on the Fock-projected $S^3$ graph satisfying the Connes axioms at finite $n_{\max}$ | Abstract; Thm 1 §III.6 | **B** (Internally consistent) | MIXED: standard Connes definitions (EXTERNAL) + GeoVac-specific construction (GEOVAC-ONLY) | Definitions and theorem statements are coherent. Order-zero, bounded commutators, compact resolvent at finite dim are trivially satisfied. The non-trivial test (real structure $J$ exists at finite $n_{\max}$ with the predicted KO-dim-3 signs) rests on internal modules and is reported as bit-exact in `tests/test_real_structure.py` per §IV — verifiable but I did not re-run. |
| 2 | "Matches $\alpha^{-1}$ at $8.8\times 10^{-8}$" (combination $K = \pi(B+F-\Delta)$, Paper 2 carried into this paper) | Abstract; §VI Rem K_under_theorem | **C** (Overstated) | EXTERNAL (CODATA α) + EXTERNAL (mpmath) | Recomputed: $K = 137.03606441448\ldots$, CODATA $1/\alpha = 137.035999084\ldots$. **Raw relative error is $4.77\times10^{-7}$, ~5.4× larger than the abstract claim of $8.8\times10^{-8}$.** The $8.8\times10^{-8}$ figure is the *post-Sprint-A residual* $|K - 1/\alpha - \alpha^2|/(1/\alpha) = 8.81\times10^{-8}$, not the raw match. The abstract sentence does not flag this distinction. Cross-corpus: CLAUDE.md §1.7 WH5 and Paper 2 both note the Sprint-A residual context, but Paper 32's abstract presents the number as the direct $K$ vs $1/\alpha$ match. (See Numbers I recomputed.) |
| 3 | $B = 42$ from $\sum_{n=1}^{3}(2l+1)l(l+1)$ on Fock shells | §VI Eq.s after K formula | **A** | EXTERNAL (direct computation) | Verified: $0 + (0 + 6) + (0 + 6 + 30) = 42$ from $(n,l) \in \{(1,0),(2,0),(2,1),(3,0),(3,1),(3,2)\}$ with weights $(2l+1)l(l+1) \in \{0, 0, 6, 0, 6, 30\}$. |
| 4 | $F = \pi^2/6 = \zeta(2)$ identification | §VI | **A** | EXTERNAL (Euler 1735) | Standard; verified. |
| 5 | $\Delta^{-1} = 40 = g_3^{\mathrm{Dirac}}$ via Camporesi–Higuchi degeneracy $g_n = 2(n+1)(n+2)$ at $n=3$ | §VI | **A** | EXTERNAL (Camporesi–Higuchi 1996) | Verified: $g_3 = 2\cdot 4\cdot 5 = 40$. |
| 6 | $\pi^3\alpha^3$ matches Sprint-A residual within 0.25% | §VI Rem sprint_A | **A** | EXTERNAL (mpmath recompute) | Verified: $K - 1/\alpha - \alpha^2 = 1.2079\times10^{-5}$, $\pi^3\alpha^3 = 1.2049\times10^{-5}$, rel diff $0.251\%$. (Caveat: a coincidence at single-CODATA-point precision, not a derivation; the paper explicitly flags this as structural hint.) |
| 7 | $\mathrm{prop}(\mathcal{O}_{n_{\max}}) = 2$ at $n_{\max} \in \{2,3,4\}$, matching Toeplitz $C(S^1)^{(n)}$ | Prop 1, §III | **B** (Internally consistent) | MIXED: Connes–vS 2021 Prop 4.2 (EXTERNAL) defines the invariant; numerical verification (GEOVAC-ONLY) at the stated dim sequences | The dimension table $N=5,14,30; N^2=25,196,900$ checks. The middle column $\dim(\mathcal{O}) = 14, 55, 140$ requires the SO(4) multiplier construction. Per paper, verified by `tests/test_operator_system.py` (24 tests, ~20s). Not independently re-run. |
| 8 | $\pi$-source case-exhaustion theorem (M1/M2/M3 partition of $\pi$ in Paper 34 projections) | Thm 2 §VIII | **C** (Overstated in scope) | MIXED | The theorem **statement says "any finite composition of projections drawn from the Paper~34 list \S\,III.1--\S\,III.15"** but Paper 34 v3.2.x catalogues **twenty-eight projections** (verified by grep on `papers/group6_precision_observations/paper_34_projection_taxonomy.tex`: "We catalogue twenty-eight projection mechanisms" at line 45, repeated at lines 220, 2736, etc.). The proof at line 1675 also says "each of the fifteen Paper~34 projections is either $\pi$-free... or $\pi$-bearing...". The theorem is asserted for the 15-projection corpus from when it was first written; projections 16–28 (Sprint 1 and Sprint 2 dictionary-completion entries, per CLAUDE.md §6 Paper 34 entry) are NOT covered. This is a corpus-drift defect: the theorem holds for the 15 it covers, but its current statement claims coverage of "any finite composition of Paper 34 projections", which is now strictly more than 15. The remark about a "sixteenth projection" as falsifier (Rem at line 1734) is internally consistent with the 15-bound theorem but inconsistent with the corpus state. |
| 9 | GH-convergence theorem $\Lambda(\mathcal{T}_{n_{\max}}, \mathcal{T}_{S^3}) \le C_3\cdot\gamma_{n_{\max}}\to 0$ with $C_3 = 1$ | Thm 3 §VIII | **B** (Internally consistent at qualitative-rate level) | MIXED: Latrémolière propinquity framework (EXTERNAL) + five-lemma proof sketch in five distinct memos (GEOVAC-ONLY for L1'…L5 in their stated assembly) | The proof is a "sketch" per §VIII; the five lemmas are each documented in separate memos (`debug/r25_l{1,2,3,4,5}_proof_memo.md`). Numerical panel $\Lambda \in \{2.075, 1.610, 1.322\}$ at $n_{\max} \in \{2,3,4\}$ is reported and monotone-decreasing. The honest scope statement is clean ("qualitative-rate; asymptotic rate $O(\log n / n)$ consistent with but not rigorously proved"). Cross-corpus: Paper 38 is the standalone writeup. The paper is appropriately cautious that this is a sketch + computational verification, not a rigorous theorem in the published sense; the "PROVEN" wording in CLAUDE.md §1.7 is internal terminology. |
| 10 | Four BBB-predicted signs $\varepsilon = +1, \varepsilon'' = -1, \kappa = -1, \kappa'' = +1$ at $(m,n) = (4,6)$ hold bit-exact on the Krein lift | §VIII.D Tab 2 (axiom_audit_lorentzian) | **B** (Internally consistent) | MIXED: BBB 2018 Table 1 lookup (EXTERNAL — Citation Checker to verify) + Krein-space construction (GEOVAC-ONLY) | Signs lifted from \cite{bizi_brouder_besnard2018}. Bit-exact verification per `tests/test_connes_axiom_audit_31.py` (75 tests). Not re-run. |
| 11 | BBB universal $\{\chi, D\} = 0$ FAILS on truthful $D_{\mathrm{GV}}$; resolution: R1 (accept) + R2 (offdiag CH) | §VIII.D Tab 2; §sec:l2c_lorentzian_dirac | **B** (Honest finding, scope clear) | GEOVAC-ONLY | The honest scope of this finding is well-stated: a structural feature of GeoVac's chirality-diagonal $D_{\mathrm{GV}}$ in the Peskin–Schroeder convention, NOT a derivation error. Three resolutions are flagged with R1+R2 recommended. This is a model of how to handle a load-bearing failure honestly. |
| 12 | H1 verdict: Higgs admitted, Yukawa not selected by GeoVac | §sec:higgs_h1 | **B** | GEOVAC-ONLY | Internally consistent decomposition into gauge + Higgs sectors of inner fluctuations. The verdict "construction admits a Higgs given imposed Yukawa" is appropriately cautious. |
| 13 | G3 verdict: $\gamma_{\mathrm{GV}}$ and $\gamma_F$ are independent commuting $\mathbb{Z}_2$'s; $\|\Delta\|_{\mathrm{op}} = 2$ exactly at every $n_{\max} \in \{1,2,3\}$ | §sec:g3 Tab 3 | **A** (structural / linear algebra) | MIXED: standard tensor-product theory (EXTERNAL) + spec verification (computational) | The spectrum $\{-2:k, 0:2k, +2:k\}$ in 1:2:1 ratios is the standard spectrum of $A \otimes I - I \otimes B$ where $A, B$ both have spectrum $\{-1, +1\}$ in equal multiplicities. This is a tautology of tensor-product spectra and the verdict is sound. |
| 14 | G4a closure (Sprint G4a, 2026-05-31): full $\mathbb{C}\oplus\mathbb{H}\oplus M_3(\mathbb{C})$ on $\mathcal{T}_{S^3}$; six Connes axioms at machine-zero residual; gauge $U(1)\times SU(2)\times SU(3)$ recovered | §sec:g4 G4a closure paragraph | **B** | GEOVAC-ONLY | New module `geovac/standard_model_triple.py`, 45 tests. Verdict is appropriately "POSITIVE-THIN" (Yukawa free). |
| 15 | "Each structural ingredient has published precedent; what is not duplicated is the $\alpha$ prediction" (lineage placement) | Abstract; §IX Obs lineage | **B** (Honest) | EXTERNAL (lineage refs) + GEOVAC-ONLY (α prediction not duplicated) | Lineage citations (Marcolli–vS, Perez-Sanchez, Connes–vS, Chamseddine–Connes) are presented as published. Citation Checker should confirm these and the more recent Hekkelman–McDonald / Latrémolière 2026 refs. The "not duplicated" claim is structurally correct: no published framework predicts α from a discrete spectral action. |
| 16 | "Four-way $S^3$ coincidence" framing reduced from WH4 maximal claim to "deflated" reading | §I (intro), WH4 footnote-style | **B** (Internally honest, cross-corpus aligned) | GEOVAC-ONLY | CLAUDE.md §1.7 WH4 explicitly deflates to "single Fock-projection statement plus three forced consequences." Paper 32's body is consistent with the deflated framing — the abstract says "four projections of the single triple constructed here," which is the structural reading. |
| 17 | Bisognano–Wichmann reading: M1 $2\pi$ on $S^1_\tau$ IS the BW modular flow period, via published Wick rotation | §VIII Rem bisognano_wichmann_reading | **B** (Structural correspondence, clean scope) | EXTERNAL (HH 1976, Sewell 1982, BW 1976, Unruh 1976) + GEOVAC-ONLY (mapping to M1) | Honestly framed as "structural correspondence, not literal identification" initially; then upgraded to "literal identification at the operator-system level (Riemannian)" via Sprint L1. The upgrade rests on the Sprint L1 finding that $\sigma_{2\pi}(O) = O$ bit-exact at finite $n_{\max}$, which I take at face value (computational, internally verifiable). |
| 18 | Four-witness theorem: HH + Sew + BW + Unruh collapse to one operator-system test via integer-spectrum K_α | §VIII Rem bisognano_wichmann_reading "Sprint L1 closure" | **B** | GEOVAC-ONLY (the collapse to a single test) | The structural mechanism — $e^{i 2\pi n} = 1$ for integer-spectrum boost generator — is mathematically correct. The "four-witness collapse" framing is GeoVac-internal vocabulary; the underlying observation is that all four witnesses share a $2\pi$ KMS period, which is standard QFT. |

---

## Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value/error | Survives? |
|-------|----------------|------------------------|---------------------------|-----------|
| $K = \pi(B+F-\Delta)$ matches $1/\alpha$ at $8.8\times10^{-8}$ | $8.8\times10^{-8}$ | CODATA $1/\alpha = 137.035999084$ (mpmath, 50 dps) | **$4.77\times10^{-7}$** raw rel. err. The $8.8\times10^{-8}$ figure equals $|K - 1/\alpha - \alpha^2|/(1/\alpha)$ — the post-Sprint-A residual, not the raw match. | **NO** (as written in abstract). The body is honest in §VI Rem sprint_A; the abstract is misleading. |
| $B = 42$ | $42$ | $\sum_n \sum_l (2l+1)l(l+1)$ for $n \in [1,3], l \in [0, n-1]$ | $42$ exactly (sympy) | YES |
| $F = \pi^2/6$ | $\zeta(2)$ | $\pi^2/6 = 1.6449\ldots$ | $\pi^2/6 = 1.6449340668\ldots$ | YES |
| $\Delta^{-1} = g_3^{\mathrm{Dirac}} = 40$ | $40$ | CH degeneracy $2(n+1)(n+2)$ at $n=3$ | $2 \cdot 4 \cdot 5 = 40$ | YES |
| $\pi^3 \alpha^3$ residual match | within 0.25% | mpmath, 50 dps | $\frac{|R - \pi^3\alpha^3|}{\pi^3\alpha^3} = 0.2513\%$ | YES (just barely; would round to 0.3% at more precision) |
| $N_{\mathrm{Dirac}}(n_{\max}) = \frac{2}{3}n(n+1)(n+2)$ | $4, 16, 40, 80, 140$ at $n_{\max} = 1..5$ | sum $\sum_{n=0}^{n_{\max}-1} 2(n+1)(n+2)$ | $4, 16, 40, 80, 140$ exact | YES |
| $N_{\mathrm{Fock}}(n_{\max}) = \sum n^2$ | $5, 14, 30$ at $n_{\max} = 2,3,4$ | direct sum | $5, 14, 30$ | YES |

---

## Circularity map (GEOVAC-ONLY chains)

The following claims, if any single upstream GeoVac result is wrong, fall:

1. **GH-convergence theorem (Thm 3)**: rests on five separate lemmas (L1' = R3.5, L2 = central Fejér, L3 = Lichnerowicz, L4 = Berezin, L5 = Latrémolière assembly), each in its own memo. The Latrémolière propinquity framework is EXTERNAL but the application to the truthful Camporesi–Higuchi spectral triple, and the bound $C_3 = 1 \cdot \gamma_{n_{\max}}$, is GEOVAC-ONLY synthesis. Paper 38 is the standalone writeup; if Paper 38's proof has a gap, this theorem's "proof sketch" framing here is honest because §VIII calls it a sketch, but a domain expert review of Paper 38 is needed to verify the lemmas.

2. **prop = 2 at every $n_{\max}$**: rests on the SO(4) multiplier basis construction (`geovac/operator_system.py`) and the vec-stack rank algorithm. The Connes–vS framework defining prop is EXTERNAL; the application to the Fock-projected $S^3$ truncation is GEOVAC-ONLY. The robustness-under-placeholder claim (Avery–Wen–Avery upgrade preserved bit-identical) is reported but I cannot independently verify without running the code.

3. **Sprint L1 modular Hamiltonian closure ($\sigma_{2\pi}(O) = O$ bit-exact)**: rests on the BW-α realization $K = J_{\mathrm{polar}}$ on the wedge spinor basis having integer spectrum. The mechanism ($e^{i2\pi n} = 1$ for integer $n$) is EXTERNAL elementary algebra; the construction of $J_{\mathrm{polar}}$ on the truncated CH spinor bundle is GEOVAC-ONLY.

4. **Master Mellin engine partition M1/M2/M3** is the Sprint TS-E1 synthesis: each individual mechanism is documented elsewhere (Paper 18, Paper 28), but the unification under a single $\mathrm{Tr}(D^k e^{-tD^2})$ form is GeoVac-only synthesis. The case-exhaustion theorem rests on this synthesis being correct AND on a complete classification of Paper 34 projections.

5. **Case-exhaustion theorem (15 projections)**: GEOVAC-ONLY synthesis; flagged in the Overstatement findings as covering only 15 of the 28 currently-catalogued projections.

---

## Overstatement findings

1. **Abstract**: "the empirically interesting property that one of its spectral-action coefficient combinations matches $\alpha^{-1}$ at $8.8 \times 10^{-8}$."

   **Issue**: As written, a reader interprets this as $|K - 1/\alpha|/(1/\alpha) \approx 8.8\times10^{-8}$. The raw relative error is $4.77\times10^{-7}$; the $8.8\times10^{-8}$ figure is the *residual after subtracting $\alpha^2$* (Sprint A). The abstract does not make this distinction.

   **Suggested replacement**: "the empirically interesting property that the residual $K - 1/\alpha - \alpha^2$ between one of its spectral-action coefficient combinations and $\alpha^{-1}$ (after subtraction of the leading $\alpha^2$ Sprint-A correction) is $1.2\times10^{-5}$, a relative residual of $8.8\times10^{-8}$ that matches $\pi^3\alpha^3$ to 0.25%."

   Alternative softening: "matches $\alpha^{-1}$ at approximately $5\times10^{-7}$ raw, with a Sprint-A residual term $\pi^3\alpha^3$ that closes the gap to $\sim10^{-8}$."

2. **§VIII Thm 2 (Case-exhaustion theorem)**: "any finite composition of projections drawn from the Paper~34 list~\cite{paper34} \S\,III.1--\S\,III.15".

   **Issue**: Paper 34 v3.2.x now catalogues twenty-eight projection mechanisms (verified by grep). The theorem and its proof reference fifteen. Projections 16–28 (rest-mass refinement, observation/temporal-window split, nuclear charge-density / Foldy–Friar, nuclear magnetization-density / Zemach, nuclear tensor multipole, PK / core-valence orthogonality, multipole / Gaunt termination, bipolar harmonic / Drake combining, symmetry / Young tableau, adiabatic / BO, coupled-channel / adiabatic curve, gauge choice, Wick rotation, apparatus identity) are NOT covered by the proof body.

   **Suggested replacement**: Either (a) restrict the theorem statement to "the first fifteen Paper 34 projections (the corpus state as of [date])" and add an explicit note that projections 16–28 are not yet audited; or (b) extend the proof body to cover projections 16–28 (a real piece of additional work). Option (a) preserves the existing proof; option (b) is the substantive fix.

3. **Abstract**: "It made the open questions about $\alpha$, the fiber rank, and the four-way coincidence into precise spectral-triple questions instead of free-floating conjectures." (Conclusion §X)

   **Issue**: Fair, but the framing "four-way coincidence" was deflated in CLAUDE.md WH4 to "one Fock-projection statement plus three forced consequences." The conclusion language sits at the older "four roles coincide" framing; the abstract is more careful ("four projections of the single triple"). Minor inconsistency — abstract and §I are aligned with the deflated reading, conclusion uses older framing.

   **Suggested replacement**: In the conclusion, rephrase "four-way $S^3$ coincidence" → "the $S^3 = \mathrm{SU}(2)$ coincidence (manifold image, Hopf base, spin carrier, gauge manifold)" or use the §I "four projections of the single triple" language verbatim.

4. **§VIII Rem bisognano_wichmann_reading**: "the Wick-rotation theorem (Hawking + Sewell + BW + Unruh) **lifts** from structural correspondence at the metric-functional level (Sprint TD Track 4, Unruh-pendant) to **literal identification at the operator-system level** (Lorentzian, finite cutoff) at every Krein cutoff."

   **Issue**: "Literal identification" is a strong claim. The verification is: $\sigma_{2\pi}(O) = O$ holds bit-exactly because $K_α$ has integer spectrum, so $e^{i 2\pi n} = 1$. This is a *necessary condition* for BW identification but it is also satisfied by any operator with integer spectrum and period $2\pi$, regardless of physical interpretation. The "literal identification" framing implies more — that the GeoVac wedge KMS state IS the Wightman BW vacuum.

   **Suggested replacement**: "lifts from structural correspondence to operator-system-level *consistency* (period closure $\sigma_{2\pi}(O) = O$ holds bit-exact at every Krein cutoff) — a necessary condition for BW identification, not a full proof of identification with the continuum Wightman vacuum." The current language is internally honest in the surrounding scope discussion ("$H_{\mathrm{local}} \ne D_W$" finding); the headline phrasing oversells slightly.

---

## What I could NOT verify (hand to a domain expert)

1. **The five-lemma GH-convergence proof in detail**. L1'–L5 in five separate memos; the assembly is a "proof sketch" per §VIII. A math.OA referee on Latrémolière propinquity would need to check that the five lemmas compose to a genuine propinquity-bound proof, especially L3's Lipschitz constant $C_3 = 1$ on the natural panel and L5's tunneling-pair height/reach bookkeeping. This is the load-bearing claim that supports CLAUDE.md's "WH1 PROVEN" status. Paper 38 is the standalone target for that review.

2. **The Avery–Wen–Avery three-Y integral implementation** and its bit-identical agreement with the placeholder at the prop=2 level. Requires an Avery-school expert.

3. **The Marcolli–vS / Perez-Sanchez lineage placement**. Whether the gauge-network framework strictly contains the GeoVac construction is a math.OA / NCG-lineage judgment call. The paper's framing is appropriately humble ("structurally a specialization") but a referee on the gauge-network literature should confirm.

4. **The BBB sign table at $(m, n) = (4, 6)$**. Paper 32 reads $\varepsilon = +1, \varepsilon'' = -1, \kappa = -1, \kappa'' = +1$ directly from arXiv:1611.07062v2 Table 1. Citation Checker should confirm the v2 PDF and the translation $(s, t) = (3, 1) \leftrightarrow (m, n) = (4, 6)$ via BBB Table 3 ($m = t + s$, $n = t - s \pmod 8 = 4, -2 \equiv 6$).

5. **Connes–vS 2021 Definition 2.39 and Proposition 4.2 (propagation number)**. Citation Checker should confirm the lines.

6. **Latrémolière 2026 (arXiv:2603.19128) "Spectral continuity of almost commutative manifolds"**. Citation Checker should confirm this exists and that "the two-infinite case is explicitly not handled" claim (used in §VIII.D W2b) matches the paper.

7. **The matter-antimatter doubling convention at KO-dim 6 of $\mathcal{A}_F$**. The paper cites Connes–Marcolli 2008 Ch. 13. Citation Checker should verify.

---

## Broadcast readiness: YELLOW

The paper is a careful synthesis with honest scope statements throughout, an explicit construction backed by reproducible code modules, and a self-consistent axiom audit. The lineage placement is appropriately humble. The major math.OA-style theorems (Thm 3 GH-convergence, Thm 2 case-exhaustion) are presented with proof sketches and named open questions.

**Three issues block GREEN**:

1. **Abstract's "$8.8 \times 10^{-8}$" figure is a residual, not a raw match.** This is the single most likely number for a hostile reader to recompute and call wrong. The body explains the distinction; the abstract does not. Fix the abstract sentence and the corresponding §VI Rem K_under_theorem language.

2. **Case-exhaustion theorem references "fifteen" Paper 34 projections; Paper 34 now has twenty-eight.** Corpus drift; either restrict the theorem scope explicitly or extend the proof to the new projections. The theorem as stated is now strictly stronger than what is proved — a math reader will flag this immediately.

3. **Conclusion uses the older "four-way coincidence" framing** while abstract+intro use the deflated "four projections of one triple" framing. Minor inconsistency, fix in conclusion.

Once those three text-level fixes land, the paper moves to GREEN at the text-audit level. The deeper math.OA review of the five-lemma proof and the lineage placement remains to be done by domain experts (Paper 38 referee + NCG-lineage referee), but those are external dependencies, not Paper 32 defects.

**A–E counts**: A: 5 | B: 11 | C: 2 | D: 0 | E: 0
**Severity counts**: Critical: 0 | High: 2 (claim 2 abstract figure, claim 8 case-exhaustion scope) | Medium: 1 (Sprint L1 "literal identification" overstatement) | Low: 1 (conclusion framing)
**Top finding**: The abstract's headline number ($8.8 \times 10^{-8}$ match to $\alpha^{-1}$) is the Sprint-A residual relative-error $|K - 1/\alpha - \alpha^2|/(1/\alpha)$, not the raw $|K - 1/\alpha|/(1/\alpha) \approx 4.77 \times 10^{-7}$. A hostile reader will recompute and call this an overstatement. The §VI body is honest about the structure; the abstract needs the same explicitness.
