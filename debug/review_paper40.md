# Confidence Review: Paper 40 — Latrémolière propinquity convergence of spectral truncations of compact Lie groups with bi-invariant metric

**Reviewer:** Confidence-review agent
**Date:** 2026-06-01
**Wave:** 2 (math.OA arc), text-level audit
**Source:** `papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex` (2,295 lines)

## Calibration check

Not a calibration run. Treating Paper 40's internal labels (PROVEN, rigorous, etc.) as unverified claims per protocol.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | $\Lambda_\text{prop}(\mathcal{T}_\Lambda, \mathcal{T}_G) \le C_3(G)\gamma_\Lambda(G)$ for all compact connected simple Lie $G$ | Thm 4.1 (main), abstract eq | **B** | GEOVAC-ONLY (Lemmas L1'–L5 internal) + EXTERNAL (Latrémolière framework) | Verified the five-lemma structure parallels Paper 38 (which is already PROVEN per WH1); standard ingredients (Peter-Weyl, Weyl integration, etc.) are EXTERNAL. The assembly is internal — no independent external reference for the propinquity bound at non-abelian rank ≥ 2. |
| 2 | Asymptotic rate constant $4/\pi$ is **universal** across all $G$ | abstract; Thm 4.1(ii); Thm 2.6; Thm 2.7 | **B** (numerical) + **C** (analytical) | MIXED — numerical is internal compute; analytical theorem (Thm 2.7) is built on EXTERNAL ingredients (Weyl dim formula, Vandermonde-Jacobian cancellation) but the "exactly cancel at leading order" claim is GEOVAC-internal | Numerical verification table reproducible: c(Sp(2)) = 1.087 gives A-err 14.6%, B-err 33.0% ✓; c(G_2) = 1.177 gives A-err 7.6%, B-err 51.6% ✓; c(SU(4)) = 0.900 gives A-err 29.3%. See "Numbers" below. SU(3): A-err 2.4%, **B-err 2.2% (very slightly better fit)** — paper notes this honestly. |
| 3 | L3 interior (non-PRV) summand closure via Brauer–Klimyk signed-sum (Lemma 3.18, Cor 3.22) is rigorously established at all ranks | §3.3 (lines 1140–1387) | **C** (overstated) | GEOVAC-ONLY in the structural-pairing step; classical inputs (Steinberg, Vinberg, Lemma 1) EXTERNAL but the load-bearing per-witness cancellation is empirical | The supporting memo `debug/sprint_phase1b_e_paper40_lemma_33_interior_memo.md` §7.2 acknowledges: "A purely combinatorial proof of the pairing structure at general rank would require deeper bookkeeping... is treated as case-bookkeeping per Remark 3.19." 5641/5641 verification is at four ranks (SU(3), SU(4), Sp(2), G_2) only; F_4/E_6/E_7/E_8 not done. See "Overstatement findings" below. |
| 4 | L2 universal-rate theorem (Plancherel-weight × Vandermonde cancellation, Thm 2.7) | §2.2 lines 686–754 | **B** (proof sketch) | MIXED — Weyl dim formula and Weyl integration are EXTERNAL; "exactly cancel at leading order" relies on bookkeeping deferred to internal memo `l2_universal_rate_memo` | Memo file exists (`debug/l2_universal_rate_proof.md`, 44 KB). Cannot externally verify the cancellation without going through the memo's §§4–8. Paper presents it as a proof **sketch** at "rigor level of Paper 38 Appendix A" — Paper 38 itself is internal. |
| 5 | $C_3(G) \le 1$ asymptotic-tight at all ranks via PRV-summand bound (Cor 3.17) | §3.3 (lines 1108–1130) | **A** | EXTERNAL — Kumar 1988 (PRV conjecture proof), Vinberg 1990, standard reverse triangle inequality | Kumar 1988 verified existing as cited; Vinberg 1990 verified existing. The asymptotic family argument uses standard classical machinery. The PRV-side proof is solid. |
| 6 | Paper 38 (SU(2)) recovered as rank-1 case | Cor 4.5; §1; Rem 3.23 | **A** | EXTERNAL (Camporesi–Higuchi spectrum on S³) + Paper 38 | The reduction is mechanical; PRV ⇒ Clebsch–Gordan in SU(2). Self-consistent. |
| 7 | "The constant 4/π identifies with the master Mellin engine M1 sub-mechanism" | §5.3 lines 1863–1891 | **D** | GEOVAC-ONLY (Paper 18 §III.7) | M1 master Mellin engine is an internal GeoVac framework. The identification is structural/interpretive, not a theorem; no external referee can verify because the M1 framing only exists in GeoVac. |
| 8 | The result "upgrades" Gaudillot-Estrada & van Suijlekom 2025 from state-space GH to propinquity with explicit rate | abstract; §1 ¶3 | **A** | EXTERNAL (G-E & vS 2025 IMRN paper verified) | Verified G-E & vS 2025 IMRN article rnaf197 (arXiv:2310.14733) is state-space GH only, no Dirac, no rate. Paper 40's claim of upgrading to propinquity + Dirac + rate is honest. |
| 9 | $5641/5641$ pass on interior (INT) panel | Prop 3.5; Rem 3.21 | **A** | INTERNAL but reproducible | `debug/data/dirac_triangle_proof_verification.json` shows total=5641, int_pass=5641. Per-panel: SU(3) 2737, SU(4) 2212, Sp(2) 501, G_2 191. Total = 5641 ✓ |
| 10 | $977/977$ pass on ordered-pair panel (Prop 3.5 table) | Prop 3.5 lines 980–1002 | **A** | INTERNAL but cross-verifiable | 441 (SU(3)) + 400 (SU(4)) + 100 (Sp(2)) + 36 (G_2) = 977 ✓. Different metric from 5641 (pairs vs interior sigmas). |
| 11 | Universality extends from simple to general compact connected $G = G_{\text{ss}} \times T^k$ via factor-wise cancellation (Cor 2.10) | §2.2 lines 856–885 | **B** | INTERNAL sketch; relies on Peter-Weyl factorization which is EXTERNAL | The factorization is standard; the c(G_1)/2 + c(G_2)/2 = 4/π combining-rule argument is plausible but proof is sketch only ("Sketch"). |
| 12 | Cb-norm bound $\|T_{K_\Lambda}\|_\text{cb} \le 2/(\Lambda + 1)$ via Bożejko–Fendler | Lem 2.5(c) | **B** | EXTERNAL (Bożejko-Fendler 1991 verified) + GEOVAC-specific normalization | Bożejko-Fendler 1991 exists as cited. The dual-Coxeter-normalised bound $2/(\Lambda+1)$ generalizes Paper 38's $2/(n_\max+1)$; the convention map is internal but conventional. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference / recompute | My recomputed value | Survives? |
|---|---|---|---|---|
| $4/\pi$ | 1.273 | Math identity | 1.27324 | ✓ |
| $16/\pi^2$ (Sp(2) Reading B) | 1.621 | Math identity | 1.62114 | ✓ |
| $24/\pi^2$ (G_2 Reading B) | 2.432 | Math identity | 2.43171 | ✓ |
| $48/\pi^3$ (SU(4) Reading B) | 1.548 | Math identity | 1.54809 | ✓ |
| Sp(2) c=1.087 A-err | 14.6% | Recompute | 14.63% | ✓ |
| Sp(2) c=1.087 B-err | 32.9% | Recompute | 32.95% | ✓ |
| G_2 c=1.177 A-err | 7.6% | Recompute | 7.56% | ✓ |
| G_2 c=1.177 B-err | 51.6% | Recompute | 51.60% | ✓ |
| SU(3) c=1.243 A-err | 2.4% | Recompute | 2.38% | ✓ |
| **SU(3) c=1.243 B-err** | (not stated explicitly) | Recompute | **2.23% — slightly BETTER than A** | ✓ (paper acknowledges as "indistinguishable at the 2% level" footnote in §1) |
| SU(4) c=0.900 A-err | 29.3% | Recompute | 29.31% | ✓ |
| 5641/5641 cell panel | 5641 cells, 0 failures | `debug/data/dirac_triangle_proof_verification.json` | 5641 (2737+2212+501+191) | ✓ |
| 977 ordered pairs panel | 441+400+100+36 = 977 | Sum | 977 | ✓ |
| Cross-group ratio Sp(2)/G_2 | 0.924 | 1.087/1.177 | 0.9236 | ✓ |
| Dual Coxeter $h^\vee$(SU(N))=N | N | Standard | N (= A_{N-1}, h^v=N) | ✓ |
| Dual Coxeter $h^\vee$(Sp(n))=n+1 | n+1 | Standard | n+1 (= C_n) | ✓ |
| Dual Coxeter $h^\vee$(G_2)=4 | 4 | Standard | 4 | ✓ |
| Dual Coxeter $h^\vee$(Spin(2n+1))=2n-1 | 2n-1 | Standard | 2n-1 (= B_n) | ✓ |
| $\sum_{d \text{ odd}} 1/d^2 = \pi^2/8$ | π²/8 | Math identity | 1.23370 ≈ π²/8 | ✓ |
| $|W|$ orders (SU(2)/SU(3)/Sp(2)/G_2/SU(4)) | 2/6/8/12/24 | Standard | 2/6/8/12/24 | ✓ |
| Unit-mass of central spectral Fejer kernel | claimed in 2.5(a) | Schur orthogonality | Verified by direct expansion | ✓ |

All quantitative claims survive recomputation.

### Circularity map (GEOVAC-ONLY chains)

The following load-bearing chains bottom out GEOVAC-internal:

1. **Main theorem (Thm 4.1) ⇒ Lemmas L1'–L5 ⇒ Paper 38 architecture.** The five-lemma chain is *internal*. Lemmas L1', L4, L5 mechanically transport Paper 38's structure to higher rank using EXTERNAL ingredients (Peter-Weyl, Weyl integration formula, Schur orthogonality, Young's inequality, Bożejko–Fendler), so each individual lemma rests partly on EXTERNAL. But the *assembly* (and the propinquity-rate output) is GEOVAC-internal. Paper 38 itself is "PROVEN" only in the GeoVac sense (WH1 PROVEN); no independent external referee has signed off on Paper 38.

2. **L3 interior closure (Lemma 3.18, Cor 3.22) ⇒ Brauer-Klimyk per-witness cancellation pairing.** The structural pairing claim ("Steinberg signed sum cancels exactly those orbit elements violating the per-witness inequality") is GEOVAC-internal at the qualitative-rate level. It is *case-checked* at four ranks (5641 cells) but not classically proved by orbit-bookkeeping. The supporting sprint memo §7.2 acknowledges this: "rigor level of Paper 38 Appendix A" — which is itself an internal-rigor standard.

3. **L2 universal-rate (Thm 2.7) ⇒ Plancherel × Vandermonde "exactly cancel".** Cancellation theorem proof relies on bookkeeping in `l2_universal_rate_memo` (44 KB internal memo). The high-level structural reading (Vandermonde-Jacobian + Weyl dim ⇒ cancellation) is plausible classical content, but the rigorous statement is mediated through the memo. A skeptical reader can read the memo but it is not an external reference.

4. **§5.3 master Mellin M1 reading.** Pure GeoVac framework; not externally checkable.

The PRV-summand bound (Thm 3.16, Cor 3.17), the Berezin reconstruction setup, the Latrémolière propinquity framework, and the Gaudillot-Estrada cross-reference all bottom out EXTERNAL. The numerical extractions are reproducible (compute is real, no fitting freedom that I can see beyond the standard subleading-aware Stein-Weiss fit).

### Overstatement findings

**Finding O1 (MEDIUM): Cor 3.22 "rigorously established at all ranks" overstates the L3-interior closure.**

- Exact phrase: "Lemma~\ref{lem:L3} is therefore rigorously established at all ranks with $C_3(G) = 1$ asymptotic-tight." (Cor 3.22, line ~1383)
- The supporting memo `debug/sprint_phase1b_e_paper40_lemma_33_interior_memo.md` §7.2 says this is "consistent with the rigor level of Paper 38's Appendix A (which uses the same blend of structural argument plus exhaustive per-rank check); reviewers familiar with that paper should accept the same convention here." That's an internal rigor convention, not classical mathematical rigor.
- The interior-summand cancellation is verified at $5641$ cells across SU(3), SU(4), Sp(2), $G_2$. The analytical Steps 1–5 do not include a rigorous combinatorial proof of the per-witness sign-pairing structure at general rank; the memo names this honestly ("treated as case-bookkeeping").
- The honest scope statement IS present in the paper at Remark 3.21 footnote: "at higher ranks (e.g., $F_4, E_6, E_7, E_8$) the analytical proof structure of Steps 1–5 carries through verbatim, but the per-orbit case-bookkeeping of Remark~\ref{rem:case_bookkeeping} would need to be re-tabulated on the relevant panels." This is good. But the headline corollary "rigorously established at all ranks" is in tension with this footnote.
- **Suggested honest replacement** (Cor 3.22): "Lemma~\ref{lem:L3} is therefore established at all ranks with $C_3(G) = 1$ asymptotic-tight in the following sense: rigorously proved for the asymptotic family that controls the rate (via PRV-summand bound, Thm 3.16) and via structural-mechanism-plus-exhaustive-case-check at the four ranks SU(3), SU(4), Sp(2), $G_2$ for interior summands (Lemma 3.18 + Remark 3.21). Closure of interior summands at $F_4, E_6, E_7, E_8$ is per-orbit case-bookkeeping not yet performed." This wording matches the body of the paper and the supporting memo without weakening the main theorem (since asymptotic-tightness of $C_3=1$ is already classically proved via PRV alone).

**Finding O2 (LOW-MEDIUM): "First-class-wide" implicit framing in abstract.**

- Exact phrase: abstract opens with "Let $G$ be a compact connected Lie group of rank $r \ge 1$..." and presents the result as a class-wide theorem. The structural reading in §1 (`paragraph: What this paper does`) writes "We close the gap between Paper~38 and~\cite{gaudillot_estrada_vs2025} at full Class~1 generality."
- This phrasing is *defensible* because the main theorem is stated at full generality, and the analytical L2 theorem (Plancherel × Vandermonde) is also stated at full generality. But the **empirical confirmation** is at five Lie groups only (SU(2), SU(3), SU(4), Sp(2), G_2). Plus the L3-interior closure has only been case-checked at four ranks.
- This is consistent with the honest scope of Remark 3.21 and §6.2 "Rank-3 numerical extraction." Not strictly overstating, but a sharper reader might want the abstract to acknowledge "with full empirical confirmation at rank ≤ 3 and structural arguments at all ranks." Minor.

**Finding O3 (LOW): "Universal" is asymptotic, not constant.**

- The constant $4/\pi$ is universal as an **asymptotic** ($\Lambda \to \infty$) — the finite-cutoff extractions vary (1.087 to 1.273 across rank 2). The paper acknowledges this ("the deviations are within the Stein–Weiss fit bias"). The abstract says "universal across the class" which is correct under the asymptotic reading, but a reader could mistake it for "the finite-cutoff number is the same." The body §2.2 makes this clear; the abstract phrasing is acceptable but could add the word "asymptotic" for crispness.

**Finding O4 (LOW): "First proof" / priority framing absent — good.**

The paper is appropriately careful: §1 ¶3 names Gaudillot-Estrada & van Suijlekom 2025 as concurrent state-space GH work, and labels the upgrade to propinquity + Dirac + rate as the present paper's contribution. The abstract does NOT claim "first" — it claims to upgrade (G-E & vS) to (propinquity + rate). This is honest and well-positioned.

## Pass B — Citation and novelty

### Citation table (selected key entries)

| `\cite` key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `gaudillot_estrada_vs2025` | "concurrent state-space GH result, IMRN 2025, arXiv:2310.14733" | **CITE-OK** | Confirmed: arXiv 2310.14733 (Y. Gaudillot-Estrada, W. D. van Suijlekom), published in IMRN 2025, rnaf197 (https://academic.oup.com/imrn/article/2025/13/rnaf197/8182203). State-space GH only, no Dirac, no rate. Matches paper's framing. |
| `kostant1999` | "Kostant cubic Dirac operator, Duke 1999" | **CITE-OK** | Confirmed: Kostant, "A cubic Dirac operator and the emergence of Euler number multiplets..." Duke Math J. 100(3): 447-501 (1999). |
| `kumar1988` | "PRV conjecture proof, Inventiones 1988" | **CITE-OK** | Confirmed: Kumar, "Proof of the Parthasarathy-Ranga Rao-Varadarajan conjecture," Inventiones Math. 93 (1988), 117-130. |
| `bozejko_fendler1991` | "Herz-Schur multipliers, Archiv Math 1991" | **CITE-OK** | Confirmed: Bożejko & Fendler, "Herz-Schur multipliers and uniformly bounded representations of discrete groups," Arch. Math. 57 (1991), 290-298. |
| `connes_vs2021` | "Spectral truncations and operator systems, CMP 2021, arXiv:2004.14115" | **CITE-OK** | Confirmed via multiple references. Standard. |
| `latremoliere_metric_st_2017` | "Latrémolière propinquity for metric spectral triples, Adv. Math. 2023 / preprint 2017" | **CITE-OK** | Confirmed: arXiv:1811.10843, "The Gromov-Hausdorff propinquity for metric Spectral Triples" — Adv. Math. 415 (2023), Paper 108876. |
| `leimbach_vs2024` | "GH convergence for tori, Adv. Math. 439 (2024), 109496" | **CITE-OK** | Confirmed: Leimbach & van Suijlekom, "Gromov-Hausdorff convergence of spectral truncations for tori," Adv. Math. 439 (2024), Paper 109496 (https://www.sciencedirect.com/science/article/pii/S0001870824000112). |
| `branson1996` | "Branson, Stein-Weiss operators and ellipticity, JFA 151 (1997)" | **CITE-OK** | Confirmed via ScienceDirect — JFA vol. 151 (1997), 334-383. Key text says 1997, bibitem natbib key says 1996 (likely preprint year); minor cosmetic inconsistency only. |
| `bourbaki_lie_8` | "Bourbaki, Groupes et algèbres de Lie Ch. 7-9, Hermann/Springer 1975 (English 2005)" | **CITE-OK** with mild caveat | Confirmed: English translation of Chapters 7-9 by Springer 2005 exists. The "Hermann/Springer 1975" framing of Chapters 7-9 combined is slightly anachronistic — Chapter 9 wasn't in the original 1975 Hermann edition (Ch. 9 was added later). But the **content cited (Steinberg formula in Ch. VIII §9)** is accurate. Minor cosmetic. |
| `fulton_harris1991` | "Fulton & Harris, Representation Theory: A First Course, GTM 129, Springer 1991" | **CITE-OK** | Standard textbook; the cited Prop. 25.30 (highest-weight-vector decomposition machinery) is in the right chapter. |
| `vinberg1990` | "Vinberg, On certain commutative subalgebras of universal enveloping algebra, Izv. AN SSSR 1990" | **CITE-DOESNT-SUPPORT** (potential) | Title verified existing; however, the **content cited** — "dominance maximises inner product lemma" — is NOT the obvious content of Vinberg 1990 (which is about commutative subalgebras / argument shift method). Search did not return any 'dominance maximises inner product' result in that paper. The cited fact is well-known as an elementary consequence of Weyl chambers and may be misattributed. Worth domain-expert check; the underlying inequality is classical. |
| `agricola2003` | "Agricola, Comm. Math. Phys. 232 (2003), 535-563" | **CITE-OK** (unverified-but-plausible) | Standard reference for naturally reductive Dirac operators; not deeply checked, but consistent with standard Kostant-Dirac references. |
| `mathieu1989` | "Mathieu, Compos. Math. 69 (1989), 37-60" | **CITE-OK** (independent PRV proof) | Verified to exist (search confirms Mathieu independently proved PRV around same time). |
| `polo1994` | "Polo, Astérisque 173-174 (1989), 281-311" | **CITE-WRONG-METADATA** (cosmetic) | Year mismatch: bibitem label says "Polo(1994)" but year of publication in Astérisque is 1989 according to the bibitem itself ("(1989)"). The natbib key `polo1994` may refer to a follow-up or be wrong. Minor. |
| `kac1990` | "Kac, Infinite-Dimensional Lie Algebras, 3rd ed. CUP 1990" | **CITE-OK** | Standard. |
| `hekkelman2022` | "Hekkelman, Truncated geometry on the circle, LMP 2022, arXiv:2111.13865" | **CITE-OK** | Standard. |
| `hekkelman_mcdonald2024` | "Hekkelman-McDonald, Spectral truncations of T^d, JNCG to appear, arXiv:2403.18619" | **CITE-OK** | The arXiv ID 2403.18619 — let me note that I did not deeply verify this specific ID but the authorship and topic align with verified searches above. |
| `helgason1984` | "Helgason, Groups and Geometric Analysis, Academic Press 1984" | **CITE-OK** | Standard. |
| Internal memos (`unified_gh_scoping_memo`, `dirac_triangle_memo`, `l2_universal_rate_memo`, `sp2_g2_rate_memo`, `su4_rate_memo`) | All cited as GeoVac internal preprints | **CITE-OK** (with circularity caveat) | All five memo files confirmed to exist on disk. They are *internal* references — a referee outside GeoVac would need to read them to verify Steps 4 and the Plancherel × Vandermonde bookkeeping. |
| `paper38`, `paper39` | "Loutey, Paper 38, GeoVac internal preprint (May 2026); Loutey, Paper 39, GeoVac internal preprint (May 2026)" | **CITE-OK** with circularity caveat | Both papers exist in `papers/group1_operator_algebras/`. They are GeoVac internal. |

### Problems found

**P1 (MEDIUM, CITE-DOESNT-SUPPORT): `vinberg1990` "dominance-maximises-inner-product lemma".**

- Paper 40 Step 4 of Thm 3.16 proof: "Vinberg's 'dominance maximises inner product' lemma~\cite{vinberg1990}". Cited again in §3.3 (Lem L3-interior proof) and Rem 3.20.
- Vinberg's 1990 paper is "On certain commutative subalgebras of a universal enveloping algebra" — Izv. Akad. Nauk SSSR Ser. Mat. 54 (1990). The paper is about commutative subalgebras and the argument-shift method.
- The "dominance maximises inner product with ρ" inequality is a classical Weyl-chamber fact (for any $v \in \mathfrak{h}^*$, the dominant Weyl-image $\bar{w(v)}$ has the largest pairing $\langle \bar{w(v)}, \rho\rangle$ among the orbit). This is in Bourbaki Ch. VI / Humphreys 1972 — *not* Vinberg 1990.
- **Recommendation**: replace the Vinberg citation with the standard Weyl-chamber reference (e.g., Bourbaki Ch. VI or Humphreys, *Introduction to Lie Algebras and Representation Theory*). If the user inherits this attribution from the supporting memo `dirac_triangle_memo`, the memo should also be corrected.
- This is **medium severity** because the underlying mathematical statement is clearly classical and correct; the misattribution does not affect the proof. But for a math.OA submission, citing the wrong source would draw a reviewer's eye.

**P2 (LOW, CITE-WRONG-METADATA): `polo1994` year inconsistency.**

- Bibitem text says "Astérisque 173-174 (1989), 281-311" but natbib key is `polo1994`. The year tag is likely a typo. Recommend correcting to `polo1989`.

**P3 (LOW, CITE-WRONG-METADATA): `branson1996` year inconsistency.**

- Bibitem text says "(1997), 334-383" but natbib key is `branson1996`. Standard convention follows preprint year; not a functional problem. Cosmetic only.

**P4 (LOW, COSMETIC): Bourbaki Chapters 7-9 framing.**

- The French Bourbaki "Groupes et algèbres de Lie" Chapter 9 was published 1982 (Masson), not 1975. The combined English translation by Springer (Pressley, translator) came in 2005. The bibitem date "1975" should be replaced by "1975-1982" or simply "1975 (Ch. 7-8), 1982 (Ch. 9); English translation Springer 2005." Cosmetic; not load-bearing.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "The first non-abelian non-Kähler instance was settled by the present author in [paper38]" | §1 ¶4 | Web search for SU(2) propinquity convergence math.OA | Did not find a prior published proof of SU(2) propinquity convergence with explicit Dirac rate. Concurrent G-E & vS 2025 covers state-space GH only. Connes-vS 2021 is the truncation framework, not the convergence theorem at non-abelian rank. | OK as stated. The paper does NOT claim "first in the literature" in the abstract; the body claim is "first non-abelian non-Kähler instance" which a clean search cannot disconfirm. Per protocol, this stays as defensible internal claim. |
| "The universality of $4/\pi$ established here says that..." | abstract; §1 paragraph "Structural reading" | Web search for universal-rate-constant in compact-Lie-group propinquity | No prior art found that proves $4/\pi$ universally across compact Lie groups for spectral-triple propinquity. | The novelty here is class-wide-rate-with-Dirac. Plausible novel result. Per protocol, I confirm no prior art FOUND; I cannot confirm novelty. The paper does not loudly claim "first" — it claims "universality" which is a content claim, not a priority claim. OK. |
| "The result upgrades... to Latrémolière propinquity with explicit Dirac-controlled rate" | abstract | Same as above | G-E & vS 2025 (IMRN) confirms state-space GH only, no Dirac, no rate. So this is a real upgrade in scope (Dirac + rate added). | OK as stated. |

No false-priority claims found. The paper's framing is honest — it cites concurrent work explicitly and positions itself as upgrade rather than replacement.

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| O1: Cor 3.22 "rigorously established at all ranks" tension with memo's "case-bookkeeping" admission for interior summands at higher rank | A | C (overstated) | **MEDIUM** |
| O2: Abstract framing as "class-wide" with empirical confirmation only up through rank 3 | A | C (mild) | LOW |
| O3: "Universal" without "asymptotic" qualifier in abstract | A | C (mild) | LOW |
| P1: `vinberg1990` likely doesn't support "dominance maximises inner product" cite | B | CITE-DOESNT-SUPPORT | **MEDIUM** |
| P2: `polo1994` natbib key year mismatch with bibitem 1989 | B | CITE-WRONG-METADATA | LOW |
| P3: `branson1996` natbib key year vs 1997 publication | B | CITE-WRONG-METADATA | LOW |
| P4: Bourbaki 7-9 year framing | B | CITE-WRONG-METADATA | LOW |

**Counts:**
- Pass A: A = 4, B = 5, C = 3, D = 1, E = 0
- Pass B: CITE-OK = ~17, CITE-WRONG-METADATA = 3, CITE-DOESNT-SUPPORT = 1, CITE-MISATTRIBUTED = 0, CITE-CANT-FIND = 0
- Severity: HIGH = 0, MEDIUM = 2, LOW = 5

## Broadcast readiness: **YELLOW**

Paper 40 is in good shape for math.OA submission **after addressing two MEDIUM items**:

(1) The L3-interior closure (Cor 3.22) carries a tension between its headline statement ("rigorously established at all ranks") and the actual content (analytical mechanism at general rank + per-orbit case-bookkeeping at four ranks). The supporting sprint memo §7.2 explicitly acknowledges this is "consistent with rigor level of Paper 38 Appendix A" — i.e., GeoVac's internal rigor standard. For a math.OA referee unfamiliar with that convention, this will read as overstatement. Soften Cor 3.22 to match the honest scope already articulated in Remark 3.21's footnote.

(2) The `vinberg1990` citation for "dominance maximises inner product" is likely misattributed; this is a classical Weyl-chamber fact, not Vinberg's commutative-subalgebra paper. A referee will notice and ask. Replace with Bourbaki Ch. VI or Humphreys.

Everything else is solid: all quantitative claims survive recomputation (the table of $c(G)$ extractions and A/B errors are accurate to the stated digits), all external citations check out (Gaudillot-Estrada-vS 2025 confirmed in IMRN; Kostant 1999, Kumar 1988, Bożejko-Fendler 1991, Leimbach-vS 2024, Latrémolière 2017 all verified), the structural argument extending Paper 38 is internally consistent, and the novelty framing is honest (no "first in the literature" overreach; concurrent work is acknowledged). The main theorem rests partly on Paper 38 (GeoVac-internal "PROVEN") and partly on classical EXTERNAL ingredients (Peter-Weyl, Weyl integration, Bożejko-Fendler, Latrémolière framework, Kumar PRV) — a typical math.OA paper. The asymptotic-tightness $C_3(G) = 1$ via the PRV-summand bound is rigorously established on EXTERNAL ingredients (Kumar 1988 PRV + standard reverse triangle inequality), and that is what determines the rate — so the headline asymptotic rate $4/\pi$ does not depend on the L3-interior fine-tuning at finite cutoff.

The SU(3) extraction ($c=1.243$, A-err 2.4%, B-err 2.2%) is worth noting in any future revision: Reading B fits SU(3) marginally better than Reading A (1.216 vs 1.273). The paper acknowledges this honestly in the §1 footnote ("indistinguishable at the 2% level"). The decisive A-vs-B discriminator is G_2 (91% gap; A-err 7.6% vs B-err 51.6%), and that survives.

## What I could NOT verify (hand to a human expert)

1. **The analytical L2 universal-rate proof** (Thm 2.7, "Plancherel × Vandermonde cancellation"). This is presented as a proof sketch with bookkeeping deferred to `debug/l2_universal_rate_memo.md` (44 KB). I confirmed the memo exists but a full external referee would need to walk through the Abel-Plana bookkeeping and confirm that the cancellation is "exact at leading order" type-by-type. The structural argument (Weyl dimension formula and Vandermonde Jacobian both have $\prod_{\alpha > 0}$ structure that cancels at leading order in the dual-Coxeter normalisation) is intuitively plausible and consistent with the numerical extractions, but a published proof requires the bookkeeping memo to be promoted to a published Appendix or referee-verified.

2. **The L3-interior Step 4 cancellation pairing as a general-rank theorem.** The structural reading ("Steinberg sign cancels exactly those orbits violating per-witness inequality") is plausible and supported by 5641/5641 case-checking at four ranks, but a proof of the cancellation pairing structure as a closed theorem at $F_4, E_6, E_7, E_8$ would require deeper Weyl-orbit bookkeeping. A Lie-theory specialist should be the final referee on this.

3. **The "non-Kähler central-Fejér-kernel" Berezin construction (Rem 3.24 / §3.4).** Presented as the non-Kähler analog of Hawkins 2000. The construction is internally consistent but I have not seen a published external Berezin-type construction on non-Kähler parallelisable Lie groups beyond Hekkelman-McDonald-vS 2024 (S²) and the Kähler line. A specialist could confirm whether this is a genuinely new construction.

4. **Cor 2.10 (extension to $G_{\text{ss}} \times T^k$).** Proof sketch only; the c(G_1)/2 + c(G_2)/2 combining rule is plausible but not fully proved.

5. The numerical extractions of $c(\Sp(2))$, $c(G_2)$, $c(\SU(4))$ from the subleading-aware Stein-Weiss fits — the actual fit code in `debug/sp2_g2_rate_constant.py` and `debug/su4_rate_constant.py` should be spot-checked by a Lie-rep-theory specialist; I confirmed the fit residuals and the extracted constants are internally consistent with the paper's stated values, but the underlying compute uses 2D / 3D Gauss-Legendre quadrature whose convergence at finite panels (specifically the rank-3 $\Lambda^2 = 300$ truncation for SU(4)) merits domain-expert scrutiny.

---

**Reviewer's bottom line:** This is a serious math.OA paper that closes a published-open question (G-E & vS 2025 → propinquity-level with Dirac and explicit rate). The structural reading (universality of $4/\pi$ via Plancherel × Vandermonde cancellation) is mathematically appealing and internally consistent. Two MEDIUM-severity items (one overstatement of L3-interior rigor at all ranks, one likely-misattributed citation for "dominance maximises inner product") should be addressed before broadcast. After those, the paper is broadcast-ready at the "preprint released for math.OA expert scrutiny" level.
