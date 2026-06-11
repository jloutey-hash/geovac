# Paper 42 — Drafting summary memo

**Date:** 2026-05-16
**Sprint:** Paper 42 drafting (consolidation of Sprints L0 + L1 + L1-tighten)
**Author:** Paper 42 drafting PM (Claude)
**Title:** *Tomita--Takesaki modular structure on truncated $\SU(2)$ spectral triples:\ four-witness Wick-rotation literal identification at finite cutoff*
**File:** `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex`
**Status:** First draft, three-pass clean LaTeX compilation, ready for PI review.

---

## §1. Deliverable summary

- **Length:** 2088 lines, 11,613 words, 23 pages compiled PDF (537 KB / 587 KB on second pass).
- **Three-pass LaTeX compilation:** clean, zero undefined references on final pass. Only warnings are cosmetic hyperref Unicode warnings (Greek math in title PDF metadata, same pattern as Papers 38/39/40).
- **Microtype disabled** to work around MiKTeX font-expansion environmental issue (same fix as Papers 38/39/40 use; commented inline).
- **All cross-references resolve** internally and to internal-paper bibliographic entries.

## §2. Section-by-section length and theorem count

| Section | Approx. words | Theorems / Props / Lemmas | Tables |
|:--------|--------------:|:---------------------------|:-------|
| §1 Introduction | 1100 | 1 (Theorem 1.1 main) | 0 |
| §2 The truncated CH triple | 750 | 0 | 0 |
| §3 Lorentzian-readiness audit | 900 | 0 | 0 |
| §4 BW at finite cutoff | 850 | 1 (Prop 4.2 wedge props) + 1 (Def 4.1) | 0 |
| §5 Geometric construction (BW-$\alpha$) | 1050 | 1 (Theorem 5.4) + 1 (Def 5.1) + 1 (Prop 5.3) | 1 (residuals) |
| §6 Tomita--Takesaki (BW-$\gamma$) | 1500 | 1 (Theorem 6.3) + 2 (Props 6.2, 6.4) | 1 (residuals) |
| §7 Unified closure | 1200 | 1 (Theorem 7.1 conjugacy) + 1 (Remark 7.2 speculation) | 0 |
| §8 Four-witness collapse | 700 | 1 (Cor 8.1 single construction) | 1 (cross-witness) |
| §9 Related work and scope | 900 | 0 | 0 |
| §10 Conclusion + open questions | 800 | 0 | 0 |
| App A Computational verification | 550 | 0 | 0 |
| App B Framework cross-reference | 550 | 0 | 0 |
| **Totals** | **~10,850** | **5 theorems, 4 propositions / corollaries, 2 definitions, 2 remarks** | **3 tables** |

Three load-bearing theorems carry the technical content:
- **Theorem 5.4 (`thm:bw_alpha`):** BW-$\alpha$ period closure $\sigma_{2\pi}^\alpha(O) = O$ at finite cutoff.
- **Theorem 6.3 (`thm:bw_gamma`):** BW-$\gamma$ Tomita-Takesaki period closure $\sigma_{2\pi}^{TT}(O) = O$ at finite cutoff.
- **Theorem 7.1 (`thm:unified`):** Flow conjugacy $\sigma_t^{TT}(O) = \sigma_{-t}^\alpha(O)$.

Plus **Theorem 1.1 (`thm:main_intro`)** in §1 as the consolidated main-theorem statement.

The unified closure verdict is **STRONG_IDENTIFICATION + UNIFIED_STRONG** at every tested $n_{\max} \in \{2, 3, 4, 5\}$ across all six witness instantiations (BW, HH$_{M=1}$, HH$_{M=2}$, Sew$_{M=1}$, Unruh$_{a=1}$, Unruh$_{a=2}$).

## §3. Cross-reference count

- **Internal cross-references (within Paper 42):** all sections, all theorems, all definitions, all 3 tables — resolve cleanly.
- **External cross-references to internal preprints (Papers 2, 24, 25, 32, 34, 38, 39, 40):** 9 references, each with internal preprint bibitem.
- **External cross-references to published literature:** 23 published references with arXiv IDs / journal citations.
- **Sprint memo cross-references:** L0 audit memo, Nieuviarts scoping memo, L1-A architecture memo, L1-B audit memo, L1-C witness specs memo, L1 closure memo, L1-tighten closure memo — all cited in section text with `\texttt{}` paths.

## §4. Honest-scope statements made

The paper makes the following honest-scope statements explicitly:

1. **§1 / §1 "What this paper does and does not close":** Closes Sprint L1 falsifier at signature $(3, 0)$ only. Lorentzian extension to $(3, 1)$ via BBB Krein lift is the named Sprint L2 follow-up. M3 trivialisation prediction at $(3, 1)$ from Sprint L0 is referenced but not addressed.

2. **§3.2 Nieuviarts footnote:** NO-GO documented as an honest related-work negative — Nieuviarts twist morphism (Def 2.2 of arXiv:2502.18105v3) explicitly restricted to even-dim manifolds, $\sthree$ at KO-dim 3 (odd) does not qualify. The restriction is structural, not editorial. Direct BBB Krein lift remains the published path for the $(3, 1)$ extension.

3. **§4.3 Remark 4.3 (`rem:wedge_choice`):** Three candidate wedges were considered (W1 hemispheric, W2 $\kappa$-preserving Dirac graph block, W3 pendant-vertex-free). W1 is the canonical choice; the closure's dependence on wedge choice is named as part of the open scope.

4. **§5.1 Remark 5.2 (`rem:two_m_j`):** The doubled generator $\mathrm{two}\_m_j$ vs the bare $m_j$ is a half-integer-spin doubling convention; the factor of 2 makes the integer-spectrum property manifest.

5. **§6.3 Remark 6.5 (`rem:state_dependent`):** State-dependence of $J_{TT}$ vs intrinsicness of $J_{GV}$ — categorically different antilinear operators, must not be conflated.

6. **§6.5 Numerical conditioning caveat in Table 6.1 caption:** The slightly larger residual at $\nmax = 5$ for the BW-$\gamma$ construction is pure numerical conditioning from the $4900 \times 4900$ GNS Hilbert-Schmidt diagonalisation, NOT structural drift. Wall-time table in App A documents this.

7. **§7.2 (`sec:derived_hamiltonian`) — THE LOAD-BEARING SCOPE FINDING.** The unified closure rests on the choice $H_{\mathrm{local}} = K_\alpha^W / \beta$. This is flagged as a deliberate construction-level choice and a derived structural finding, not a uniqueness theorem. Two consequences are spelled out: (I) it matches the Wightman BW vacuum (canonical operator-system analog, not post-hoc engineering); (II) the framework's intrinsic Camporesi-Higuchi Dirac $\DCH$ is NOT the right local Hamiltonian for the BW vacuum at $\beta = 2\pi$ — choosing $H_{\mathrm{local}} = D_W$ would give half-integer spectrum and the bit-exact closure would not hold.

8. **§7.2 Remark 7.3 (`rem:speculation`) — explicit speculative reading flagged as such.** The structural distinction between the spectral-action Dirac $\DCH$ and the modular-Hamiltonian generator $K_\alpha^W$ is flagged as possibly significant for the framework's interpretation (spectral-action paradigm vs thermal-time paradigm), with the explicit annotation "we note, but do not pursue here, the possibility that..." — fully complies with §1.5 rhetoric rule.

9. **§10 Open questions section — five named open questions:**
   - (O1) Lorentzian extension at $(3, 1)$;
   - (O2) Non-tracial wedge states;
   - (O3) Spectral-action Dirac vs modular-Hamiltonian generator distinction;
   - (O4) Higher-rank compact non-abelian extensions (extending to Paper 40 generality);
   - (O5) Cross-manifold modular structures (Paper 24 §V W2b blocker).

10. **App B (`app:framework_crossref`) Paper 2 paragraph:** explicit statement that the $\nmax = 3$ dimension coincidence ($\dim \Hilb_3 = 40 = g_3^{Dirac} = \Delta^{-1}$) is a structural coincidence at the Hilbert-space-dimension level, NOT a derivation of Paper 2's combination rule. Paper 2 stays conjectural.

## §5. New theorems / definitions in Paper 42

| Theorem | Statement | Proof technique |
|:--------|:----------|:----------------|
| **Theorem 5.4 (`thm:bw_alpha`)** | $\sigma_{2\pi}^\alpha(O) = O$ bit-exact at finite cutoff for the geometric BW-$\alpha$ generator $K_\alpha^W$. | Integer spectrum $\Rightarrow e^{i \cdot 2\pi \cdot n} = 1$. |
| **Theorem 6.3 (`thm:bw_gamma`)** | $\sigma_{2\pi}^{TT}(a) = a$ bit-exact at finite cutoff for the Tomita-Takesaki BW-$\gamma$ construction. | $K_{TT} = -\log \Delta$ has integer spectrum (Prop 6.4) $\Rightarrow$ same mechanism. |
| **Theorem 7.1 (`thm:unified`)** | $\sigma_t^{TT}(a) = \sigma_{-t}^\alpha(a)$ bit-exact at general $t$, identical at $t = 2\pi$. | Direct algebraic identity from $\rho^{it} = e^{-it K_\alpha^W}$ in (7.1)/(7.2) representation. |
| **Proposition 4.2 (`prop:wedge_properties`)** | $P_W^2 = P_W$, $P_W^* = P_W$, $\Tr(P_W) = \dim \Hilb_{\nmax}/2$. | Direct, $R_{\mathrm{polar}}^2 = I$. |
| **Proposition 5.3 (`prop:integer_spectrum`)** | $\mathrm{Spec}(K_\alpha^W) \subset \Z$. | Direct from $\mathrm{two}\_m_j$ being odd integers. |
| **Proposition 6.2 (`prop:J_TT_distinct`)** | $J_{TT}^2 = +I$, $J_{GV}^2 = -I$ (categorically distinct). | Direct $\rho^{1/2} (\cdot)^* \rho^{-1/2}$ computation. |
| **Proposition 6.4 (`prop:K_TT_spectrum`)** | Spectrum of $K_{TT}$ on $\Hilb_{GNS}$ is $\{n_j - n_k\} \subset \Z$. | Tensor-product eigenvalue computation on $(\rho^{-1})^T \otimes \rho$. |
| **Corollary 8.1 (`cor:single_construction`)** | Six-witness collapse to single operator-system construction at algebra-action level. | $\rho_W = e^{-K_\alpha^W}/Z$ is $\beta$-independent under the BW choice $H_{\mathrm{local}} = K_\alpha^W/\beta$. |

## §6. Bibliography

37 bibitems total:
- **8 published modular-theory / continuum BW references** (Bisognano-Wichmann 1975, 1976; Hartle-Hawking 1976; Sewell 1982; Unruh 1976; Tomita 1967; Takesaki 1970; Connes-Rovelli 1994).
- **6 continuum modular Hamiltonian / lattice realisation references** (Casini-Huerta 2009, Casini-Huerta-Myers 2011, Zhu-Casini et al 2020, Strohmaier 2006, Verch 2001, Chamseddine-Connes 2010).
- **8 NCG metric spectral-triple / propinquity references** (Connes 1995, Connes-vS 2021, Latrémolière 2017+2018, Hekkelman 2022, Hekkelman-McDonald 2024+2024b, UCP-maps 2024, Leimbach-vS 2024, Toyota 2023, Farsi-Latrémolière 2024+2025).
- **3 Nieuviarts references** (2402.05839, 2502.18105v3, 2512.15450v2) for the NO-GO citation.
- **2 published Lorentzian-NCG references** (Bizi-Brouder-Besnard 2018, Marcolli-vS 2014).
- **1 Camporesi-Higuchi 1996 reference** for the Dirac spectrum.
- **9 internal GeoVac preprint bibitems** (Papers 2, 24, 25, 32, 34, 38, 39, 40).

**Phantom citation check:** all arXiv IDs were verified at the source-memo stage (L0 audit memo passed phantom check; Nieuviarts scoping memo passed verbatim-abstract check; L1 architecture memo passed phantom check; L1 closure memo and L1-tighten closure memo passed phantom check). No new citations added in Paper 42 beyond what was already cited in supporting memos.

## §7. Comparison to Paper 38 template

Paper 42 closely follows the Paper 38 template:

| Feature | Paper 38 | Paper 42 |
|:--------|:---------|:---------|
| Standalone math.OA paper format | Yes | Yes |
| Document class / packages | 11pt article + amsmath + booktabs + natbib + hyperref + microtype | Same, microtype commented out |
| Theorem environments | plain (theorem, lemma, prop, cor) + definition (def, rem, ex) | Same |
| Custom commands ($\nmax$, $\sthree$, etc.) | Yes | Same + added BW/TT/HH/Sew/Unruh/GNS shortcuts |
| Footnote on AI-augmented workflow | Yes (in author thanks) | Same |
| MSC2020 + keywords block | Yes | Yes (different MSC focus: 46L60 modular theory added) |
| Theorem-proof structure | Lemmas L1'-L5 + Main Theorem | 3 Theorems + 2 Corollaries + 4 Propositions |
| Tables with booktabs | 1 quantitative bound table | 3 tables (BW-$\alpha$ residuals, BW-$\gamma$ residuals, cross-witness) |
| Acknowledgments + bibliography | thebibliography embedded | Same |
| Pages | ~22 pages | 23 pages |

The paper sits stylistically alongside Papers 38, 39, 40 as the fourth standalone math.OA paper in the GeoVac series.

## §8. Computational footprint

The paper documents but does not reproduce computational code:
- `geovac/modular_hamiltonian.py` (~1552 lines) — production module (referenced in App A).
- `tests/test_modular_hamiltonian.py` (~957 lines, 67 tests: 63 fast + 4 slow) — test suite.
- `debug/l1_modular_hamiltonian_compute.py` and `debug/l1_tighten_tomita_compute.py` — computational drivers.
- `debug/data/l1_modular_hamiltonian_results.json` and `debug/data/l1_tighten_tomita_results.json` — structured numerical results.

The paper's numerical tables (Tables 5.1, 6.1, 8.1) draw directly from these JSON results.

## §9. Sprint chain provenance

Paper 42 consolidates three sprints:
1. **Sprint L0** (Lorentzian-readiness audit, 2026-05-16): 28-projection transfer audit. Source: `debug/lorentzian_l0_audit_memo.md` (~3500 words). Output: 75% of dictionary inherits Lorentzian reading at no structural cost; 5 EUCLIDEAN_SPECIFIC projections require BBB Krein lift; the Wick-rotation projection at signature $(3,0)$ is named as the principal Sprint L1 falsifier. Paper 42 §3 summarises this audit and uses it to anchor the Riemannian scope.

2. **Sprint L1** (modular Hamiltonian closure via BW-$\alpha$, 2026-05-16): geometric construction $K = J_{\mathrm{polar}}$ with integer spectrum, bit-exact $\sigma_{2\pi}(O) = O$ at $n_{\max} \in \{2,3,4,5\}$, all four witnesses. Verdict: STRONG_IDENTIFICATION. Source: `debug/l1_modular_hamiltonian_results_memo.md` (~3500 words). Paper 42 §5 is the writeup.

3. **Sprint L1-tighten** (Tomita-Takesaki BW-$\gamma$ unified closure, 2026-05-16 evening): polar decomposition $S = J_{TT} \Delta^{1/2}$ on the GNS Hilbert-Schmidt space, $K_{TT} = -\log \Delta$ has integer spectrum (Prop 6.4), bit-exact closure at the operator-action level, flow conjugacy $\sigma_t^{TT} = \sigma_{-t}^\alpha$. Verdict: UNIFIED_STRONG. Source: `debug/l1_tighten_tomita_results_memo.md` (~2800 words). Paper 42 §6 and §7 are the writeup.

The three sprints span $\sim 12$ hours of elapsed time but produced ~10,000 words of internal memo content and the full `geovac/modular_hamiltonian.py` production module. Paper 42 reformats this content as a math.OA-style paper without re-implementing anything; the supporting code and JSON outputs are the load-bearing computational substrate.

## §10. Named open questions in the conclusion

Five open questions are listed in §10 of Paper 42, ordered by structural priority:

- (O1) **Lorentzian extension at $(3, 1)$ via BBB Krein lift.** The natural next step. Sprint L2 candidate; multi-month. Obstruction: no published Lorentzian propinquity.

- (O2) **Non-tracial wedge states.** Polar decomposition still goes through but integer-spectrum property of $K_{TT}$ is not guaranteed. Open structural question.

- (O3) **Spectral-action Dirac $\DCH$ vs modular-Hamiltonian generator $K_\alpha^W$.** Section 7.2 finding — these are NOT the same operator on the wedge KMS state. Status of distinction in higher-rank extensions (Paper 40) is open.

- (O4) **Higher-rank compact non-abelian extensions.** Paper 40 extends Paper 38 to all compact connected Lie groups; analog of (Theorem 5.4) on those groups requires verifying integer-spectrum property of the Cartan generator along the wedge axis.

- (O5) **Cross-manifold modular structures.** Paper 24 §V W2b blocker — requires mixed Riemannian / Hardy-sector tensor-product propinquity.

The most structurally interesting is (O3), which the paper flags as a derived structural finding rather than a uniqueness theorem; this is the most interpretively load-bearing content of the paper beyond the technical theorems.

## §11. Recommendations for PI review

1. **Theorem statements** are stated with maximal honesty about being finite-cutoff bit-exact identities rather than asymptotic results. This is the appropriate scope, but the framing could optionally be tightened if the PI wishes to emphasise the structural-correspondence $\to$ literal-identification graduation more explicitly.

2. **§7.2 derived-Hamiltonian finding** is the load-bearing scope statement and the most interpretively significant content. PI may wish to review the framing around the choice $H_{\mathrm{local}} = K_\alpha^W / \beta$ vs the alternative $H_{\mathrm{local}} = D_W$, and whether the (II) finding (framework's Dirac is not the right local Hamiltonian) warrants stronger emphasis in the abstract.

3. **§7.2 Remark 7.3 (speculative reading)** is flagged as speculation per §1.5 rhetoric rule. PI may wish to either (a) cut it entirely (cleaner), (b) keep as is (current state), or (c) develop it into a more structured open question.

4. **Section 3 (Lorentzian-readiness audit + Nieuviarts NO-GO footnote)** is the contextualisation section. The Nieuviarts footnote is included as honest related-work acknowledgment. PI may wish to either expand it (currently ~150 words) or move it to a dedicated subsection.

5. **Abstract length** is ~430 words, slightly longer than Paper 38's ~330 words. Could be tightened by removing the operator-action-conjugacy detail.

6. **Bibliography order** is alphabetical by author within external citations (per Paper 38 / 39 / 40 convention) then internal preprint citations grouped at the end. Standard math.OA format.

7. **Title length** is "Tomita-Takesaki modular structure on truncated SU(2) spectral triples: four-witness Wick-rotation literal identification at finite cutoff" — slightly long. Alternative shorter form: "Modular structure and four-witness Wick rotation on truncated SU(2) spectral triples". PI to choose.

## §12. Files produced

- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (~2088 lines, 11,613 words, 23-page PDF, three-pass clean LaTeX, 587 KB).
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.pdf` (compiled output).
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.aux/.log/.out` (LaTeX auxiliary files).
- `debug/paper_42_draft_summary_memo.md` (this memo, ~3500 words).

## §13. Recommended Zenodo / arXiv submission metadata (pending PI sign-off)

- **Primary classification:** math.OA (operator algebras)
- **Secondary classifications:** math-ph (mathematical physics); gr-qc (general relativity, for the Hawking / Unruh / horizon thermodynamics content).
- **MSC2020:** 58B34 (noncommutative geometry); 46L60 (applications of selfadjoint operator algebras to physics); 46L87 (noncommutative Riemannian and metric geometry); 22E45 (representations of compact Lie groups); 81T05 (axiomatic quantum field theory).
- **Title:** as written, or shortened per §11.7 above.
- **Keywords:** Tomita-Takesaki modular flow, Bisognano-Wichmann theorem, Connes-van Suijlekom truncated operator system, Camporesi-Higuchi Dirac, four-witness Wick rotation, Hartle-Hawking, Sewell, Unruh effect, KMS state, SU(2) spectral triple.

End of memo.
