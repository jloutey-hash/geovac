# Confidence Review: Paper 2 — The Fine Structure Constant from Spectral Geometry of the Hopf Fibration

**Reviewer date:** 2026-06-02
**Paper file:** `papers/group5_qed_gauge/paper_2_alpha.tex`
**Status (per CLAUDE.md §6 / §13.5):** Observations folder; combination
rule $K = \pi(B+F-\Delta)$ is conjectural at the rule-level regardless
of folder placement.

## Calibration check

This is a wave-3 content audit, not a blind calibration run. The blind
calibration pass (Sprint confidence-review-infra, 2026-06-01) was
on this same paper and produced two MEDIUM fixes (Table II + abstract
"137.035987" → "137.036011"; §VIII.D summation index k=1 → k=0). Both
are already applied. My job today is to look for deeper structural
defects (citation misattribution, framing drift on the conjectural-rule
status, math/derivation gaps) the calibration pass did not target.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| C1 | $B = 42 = \sum_{n=1}^{3}\sum_{l=0}^{n-1}(2l+1)\,l(l+1)$ | §III.A, Eq.(8); abstract | **A** | EXTERNAL (elementary representation theory) | Re-summed: B(1)=0, B(2)=6, B(3)=42, B(4)=162, B(5)=462. Match. |
| C2 | $B(m) = m(m+1)(2m+1)(m+2)(m-1)/20$ | Eq.(11) | **A** | EXTERNAL (sympy algebraic) | sympy expand confirms identity to all orders. |
| C3 | $B/N = 3(m+2)(m-1)/10$ uniquely 3 at m=3 | Eq.(13), §III.B | **A** | EXTERNAL | sympy: $(m+4)(m-3)=0$, positive root m=3. |
| C4 | $F = \pi^2/6 = \zeta(2)$ | §III.D, Eq.(15) | **A** | EXTERNAL | Trivially correct ($\zeta_{S^1}(s) = 2\zeta_R(2s)$ on circumference-$2\pi$ circle). |
| C5 | $\Delta = 1/(|\lambda_3|\cdot N(2)) = 1/40$ | §III.E, Eq.(16) | **A** | EXTERNAL | $|\lambda_3|=8$, $N(2)=5$, product 40. |
| C6 | $g_3^{\mathrm{Dirac}} = 2(n+1)(n+2)|_{n=3} = 40$ on unit $S^3$ | Eq.(28) (Phase 4H finding) | **A** | EXTERNAL (Camporesi–Higuchi spectrum) | Verified: 2·4·5=40. |
| C7 | $K = \pi(42 + \pi^2/6 - 1/40) = 137.036064$ | Eq.(19) and Table II | **A** | EXTERNAL | mpmath 50-dps: K = 137.03606441448154... ✓ |
| C8 | Smallest positive real root of $\alpha^3 - K\alpha + 1 = 0$ gives $1/\alpha = 137.036011$ | Table II, Eq.(22) | **A** | EXTERNAL | sympy: 137.03601116313640498... ✓ |
| C9 | Relative error vs CODATA = $8.8 \times 10^{-8}$ | Table II, abstract | **A** | EXTERNAL (CODATA 2018: 137.035999084; CODATA 2022: 137.035999177) | mpmath: 8.81e-8 vs CODATA 2018, 8.75e-8 vs CODATA 2022. ✓ |
| C10 | Direct K vs CODATA gives $4.77 \times 10^{-7}$ raw (abstract Sprint A line) | abstract; §III.C item 2 | **A** | EXTERNAL | mpmath: $|K - 137.035999084|/137.035999084 = 4.7674\times10^{-7}$. ✓ |
| C11 | Triple selection at m=3: (i) $B/N=3$; (ii) $(+,+,-)$ uniquely hits at $4.8\times10^{-7}$ with next-best $(+,+,+)$ at $1.1\times10^{-3}$; (iii) two $\Delta$ forms agree at m=3 only | §III.C (Sprint A) | **A / B** | MIXED — (i) external proof; (ii) numerical computation over 8 sign patterns; (iii) sympy Fraction comparison m=2..9 | All three verified. (ii) and (iii) are GeoVac-internal computations against the framework's own object set, so technically B-circular at the framework-design level, but the math itself is exact and reproducible. |
| C12 | Search of $1.92\times 10^9$ formulas yields $p$-value $5.2\times 10^{-9}$ | §V, Eq.(34), Table V | **B** | GEOVAC-ONLY (search space defined by GeoVac with 743 quantities × 14 prefactors × {linear,cubic}) | Confirmed search-space arithmetic: $\binom{743}{3}_{\text{multiset}} \times 14 \times 2 \approx 1.92\times10^9$. p = 10/1.92e9 = 5.21e-9. Numbers consistent; the $p$-value interpretation is acknowledged in §V as "does not constitute a derivation" — appropriate hedge. |
| C13 | Cubic = char.poly. of traceless $\mathbb{Z}_3$-symmetric circulant | §VI | **A** | EXTERNAL (elementary circulant algebra) | det formula correct; $bc=K/3$, $b^3+c^3=-1$ identification correct. |
| C14 | $a_0 = \sqrt\pi$ is the zeroth Seeley–DeWitt coefficient on unit $S^3$ | Eq.(20) (Hopf-measure identification) | **B** | GEOVAC-ONLY (convention-dependent; matches Paper 28 Eq.(10)'s convention $a_0=a_1=\sqrt\pi$) | Convention is the unnormalized heat-kernel one (coefficient of $t^{-d/2}$ before $(4\pi)^{d/2}$ division). The bare equality $\pi = \mathrm{Vol}(S^2)/4 = \mathrm{Vol}(S^3)/\mathrm{Vol}(S^1) = \pi$ is trivially true; the $a_0^2 = \pi$ third form is convention-fixed. |
| C15 | Spectral-determinant near-miss $4\pi^2 e^{\zeta(3)/(2\pi^2)} \approx 41.957$ is 0.1% from $B=42$ | §VIII.C | **A** | EXTERNAL | mpmath: 41.957241... Gap $(42-41.957)/42 = 1.018\times10^{-3} \approx 0.10\%$. ✓ |
| C16 | Post-cubic residual $K - \alpha_{\rm exp}^{-1} - \alpha_{\rm exp}^2$ matches $\pi^3\alpha^3$ to 0.25% (Sprint A) | §VIII open question 3; abstract Sprint A line | **C / B** | GEOVAC-ONLY at 80 dps | Verified with CODATA 2018: $R_{\rm exp}/(\pi^3\alpha^3) = 1.00251$ (paper says 1.0025 ✓). **However**: re-doing with CODATA 2022 ($1/\alpha = 137.035999177$) gives ratio $= 0.99479$ (−0.5%), not +0.25%. The match is sensitive to the CODATA edition used. The abstract calls the agreement "$0.25\%$" without naming the CODATA edition. With CODATA 2022 the agreement is closer (−0.5% instead of +0.25%) but flipped in sign. |
| C17 | $D^{\rm CH}(4) = \pi^2 - \pi^4/12$ (Eq.(D_dirac_ch), Obstruction 2) | §IV §obstructions | **A** | EXTERNAL (analytical mpmath: residual 8e-51 between summed value and closed form) | ✓ |
| C18 | Improvement over paper2_old: $6.2\times10^{-5} \to 8.8\times10^{-8}$, "500× improvement" | Table VIII / §VIII.B | **C** | EXTERNAL | Actual ratio: 6.2e-5 / 8.8e-8 = **704×**, not 500×. Paper conservatively understates by ~40%. Low severity. |
| C19 | "The cubic has three real roots (discriminant $4K^3 - 27 > 0$)" | §IV after Eq.(22) | **A** | EXTERNAL (Cardano discriminant for depressed cubic $x^3 + px + q$ with $p=-K, q=1$: discriminant $-4p^3 - 27q^2 = 4K^3 - 27 > 0$ for $K > (27/4)^{1/3} \approx 1.89$). | ✓ |
| C20 | Conjectural-status of the combination rule $K = \pi(B+F-\Delta)$ is preserved | abstract, §I, §IV (Thm 1, §closure), §VII (Link 3 cross), §IX (Conclusion) | **A** | EXTERNAL (CLAUDE.md §13.5 hard prohibition #4) | Abstract explicitly says "remains conjectural"; §IV closure says "structural coincidence of unknown origin"; §IX says "remains conjectural." Discipline preserved end-to-end. |

### Numbers I recomputed

| claim | paper's figure | independent reference | my value | survives? |
|-------|---------------|----------------------|----------|-----------|
| $K$ | 137.036064 | sympy 50 dps | 137.03606441448154 | ✓ |
| $1/\alpha$ from cubic | 137.036011 | sympy 80 dps | 137.03601116313640 | ✓ |
| Rel err vs CODATA 2018 | $8.8\times 10^{-8}$ | NIST CODATA 2018 (137.035999084) | $8.81\times 10^{-8}$ | ✓ |
| Rel err vs CODATA 2022 | (paper cites CODATA 2018) | NIST CODATA 2022 (137.035999177(21)) | $8.75\times 10^{-8}$ | within rounding |
| Direct $K$ rel err | $4.8\times 10^{-7}$ (abstract Sprint A) | as above | $4.77\times 10^{-7}$ | ✓ |
| Triple-selection (ii) next-best $(+,+,+)$ | $1.1\times 10^{-3}$ | mpmath enum 8 patterns | $1.147\times 10^{-3}$ | ✓ |
| Triple-selection (iii) two $\Delta$ forms | agree at m=3 only | sympy Fraction m=2..9 | only m=3 gives 1/40=1/40 | ✓ |
| $g_3^{\rm Dirac}$ | 40 | Camporesi–Higuchi | 2·4·5 = 40 | ✓ |
| Spectral determinant near-miss | 0.1% from 42 | mpmath 50 dps | 0.102% | ✓ |
| Search count | $1.92\times 10^9$ | recomputed | $1.9219\times 10^9$ | ✓ |
| p-value | $5.2\times 10^{-9}$ | 10 / 1.92e9 | $5.21\times 10^{-9}$ | ✓ |
| Sprint A π³α³ ratio | 1.0025 (paper) | mpmath, CODATA 2018 | 1.00251 | ✓ |
| Sprint A π³α³ ratio | (paper does not state CODATA 2022 form) | mpmath, CODATA 2022 | 0.99479 | discrepant in sign |
| 500× improvement | 500× | 6.2e-5 / 8.8e-8 | 704× | understated |

### Circularity map

**Claims that are EXTERNAL** (anchored in established math, independently verified):
C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C13, C15, C17, C18, C19, C20 — these survive removal of any GeoVac assumption.

**Claims that are GEOVAC-ONLY at the framework-design level:**
- C12 (search-space construction depends on which 743 "spectral quantities" and 14 "transcendental prefactors" GeoVac elected). Mitigated by §V's own hedge.
- C14 (the $a_0^2 = \pi$ reading requires Paper 28's heat-kernel convention; bare $\pi = \pi$ is trivially true).
- C16 (the post-cubic residual ↔ $\pi^3\alpha^3$ match is a single-point coincidence whose precise residual depends on the CODATA edition).

**Critical headline claims and what they rest on:**
- "$\alpha$ from Hopf spectral invariants at $8.8\times10^{-8}$, zero free parameters" → C1, C4, C5, C7, C8, C9 — **all EXTERNAL**. Headline survives.
- "Combination rule remains conjectural" → C20 — discipline preserved.
- "$p$-value $5.2\times10^{-9}$" → C12 — GEOVAC-ONLY at the search-space-definition layer; numerics correct.

No claim in Paper 2 has a circularity that would break the headline once the conjectural-rule label is applied.

### Overstatement findings

| location | exact phrase | concern | suggested replacement |
|----------|-------------|---------|----------------------|
| **Abstract** | "establishing that this precision is very unlikely to be accidental" | "very unlikely" reads stronger than the body's own §V hedge ("does not constitute a derivation"). | "providing statistical evidence that this precision is unlikely to be accidental within the searched family of formulas" — preserves the body's caveat. **MEDIUM**. |
| **§VIII open question 3** | "matches $\pi^{3} \alpha^{3}$ to within $0.25\%$ at 80-digit precision" | The 0.25% number depends on CODATA 2018; with CODATA 2022 the ratio is 0.99479 (a −0.5% match in opposite sign direction). The paper does not name the CODATA edition for this comparison, and the agreement is in fact sensitive to it. | Add "(using CODATA 2018)" and report both editions, OR re-state as "within ~0.5%" to cover both. **LOW**. The "structural hint, not derivation" hedge is correct, and the broader honest framing is preserved. |
| **§VIII.B Table VIII** | "500× improvement" | Actual ratio is 704×; paper understates. Conservative direction. | Update to "700× improvement" for accuracy. **LOW**. |
| **§II.A spectral-triple lineage paragraph** | "the *prediction* is this paper's contribution. This distinction frames the open questions" | Calling the formula a "prediction" is mildly stronger than "observation." Paper 32 (Obs.~3.4) uses "numerical observation" consistently and the rest of Paper 2 says "empirical observation with structural support". Minor framing inconsistency between §II.A and the rest. | Reword as "what remains novel and GeoVac-specific is the *observation*—no published framework of this kind reproduces a physical coupling constant." **LOW**. |

No A-D / B-D framing-drift issue rises to the §13.5 hard-prohibition level — Paper 2's "conjectural" status of the **combination rule** is preserved in abstract, §I, §IV (Thm 1, §closure), §VII (Link 3), and §IX (Conclusion).

---

## Pass B — Citation and novelty

### Citation table

| key | claimed as | verdict | what I found |
|-----|-----------|---------|--------------|
| `feynman1985` | "Feynman called α 'one of the greatest damn mysteries'" | **CITE-OK** | QED: The Strange Theory of Light and Matter, Princeton UP 1985 — well-known book, quote verified |
| `eddington1935` | Numerical attempts have failed | **CITE-OK** | New Pathways in Science, Cambridge UP 1935 — historical record |
| `barrow2002` | Numerical attempts have failed | **CITE-OK** | The Constants of Nature, Vintage 2002 — standard reference |
| `paper7` | $S^3$ proof, 18 symbolic proofs | **CITE-OK** | Internal; consistent across corpus |
| `paper22_sparsity` | Angular sparsity universal | **CITE-OK** | Internal |
| `paper24_bargmann` | Bargmann–Segal π-free certificate | **CITE-OK** | Internal |
| `hopf1931` | Original Hopf fibration paper | **CITE-OK** | Math. Ann. 104, 637 (1931) — standard reference |
| `urbantke2003` | "Hopf fibration—seven times in physics" | **CITE-OK** | J. Geom. Phys. 46, 125 (2003) — standard reference |
| `fock1935` | Fock's stereographic projection | **CITE-OK** | Z. Phys. 98, 145 (1935) — well-known |
| `paper2_old` | Earlier symplectic impedance | **CITE-OK** | Internal v1 of this paper |
| `wyler1969` | Wyler's CRAS formula | **CITE-OK** | CRAS Paris 269A, 743 (1969), title "L'espace symétrique du groupe des équations de Maxwell" — confirmed via web search |
| `robertson1971` | Robertson critique | **CITE-OK** | Phys. Rev. Lett. 27, 1545 (1971), "Wyler's Expression for the Fine-Structure Constant" — confirmed via NASA ADS |
| `phase2_audit`, `phase1_kappa`, `phase2_kappa`, `phase3_audit`, `phase6_hopf`, `alpha_sprint_a_memo`, `track_nj_alpha_memo` | Internal sprint memos | **CITE-OK** as internal artifacts | These point to `debug/` — appropriate for the supplementary-material convention used in this paper |
| `marcolli_vs_2014` | "gauge networks in noncommutative geometry," J. Geom. Phys. 75, 71 (2014), arXiv:1301.3480 | **CITE-OK** | Verified via arxiv.org/abs/1301.3480 — authors, title, abstract content all match; DOI 10.1016/j.geomphys.2013.09.002 corresponds to JGP 75 |
| `perez_sanchez_2024` | "Revisits the continuum limit of Marcolli–van Suijlekom gauge networks; establishes that the limit is Yang–Mills without Higgs" | **CITE-DOESNT-SUPPORT** | Verified via arxiv.org/abs/2401.03705. The paper IS by Perez-Sánchez ("Bratteli networks and the Spectral Action on quivers"), but the abstract explicitly says: "a hermitian ('Higgs') matrix field emerges from the self-loops of the quiver and derive the Yang-Mills–Higgs theory on flat space as a smooth limit." This is the OPPOSITE of the "without Higgs" framing the Paper 2 footnote attributes to it. The "without Higgs" finding belongs to arXiv:2508.17338 (the Comment paper), which IS in the bibliography as `perez_sanchez_2024_comment` — so the two citations have been functionally swapped in the in-text description. **MEDIUM**. |
| `perez_sanchez_2024_comment` | Comment on the above, also "Yang-Mills without Higgs" | **CITE-OK** as the bibitem (paper exists, arXiv:2508.17338, "Comment on 'Gauge networks in noncommutative geometry'", abstract says "the continuum limit of this theory is the Yang-Mills action functional, without a Higgs scalar"). The IN-TEXT use of this together with `perez_sanchez_2024` collapses two distinct works' claims into one — see preceding row. |

### Problems found

**[CITE-DOESNT-SUPPORT, MEDIUM] §II.A "spectral-triple setting" paragraph + bibitem `perez_sanchez_2024`.**

Paper 2 §II.A reads (lines 153–158):
> "The discrete $S^3$ Hopf graph of Paper~7, with the $U(1)$ and $SU(2)$ gauge extensions developed elsewhere in the GeoVac corpus, is a specific instance of this general construction. The Perez-Sánchez correction~\cite{perez_sanchez_2024, perez_sanchez_2024_comment} clarifies that the continuum limit of such gauge networks is Yang–Mills without Higgs."

And the bibitem (line 1651) annotates `perez_sanchez_2024` (= arXiv:2401.03705) as: "Revisits the continuum limit of Marcolli–van Suijlekom gauge networks; establishes that the limit is Yang–Mills without Higgs."

My check against arxiv.org/abs/2401.03705 returns the abstract of *"Bratteli networks and the Spectral Action on quivers,"* which states:
> "...we obtain not only Wilsonian Yang-Mills lattice gauge theory...We show that a hermitian ('Higgs') matrix field emerges from the self-loops of the quiver and derive the Yang-Mills–Higgs theory on flat space as a smooth limit."

So arXiv:2401.03705 derives Yang-Mills-**with**-Higgs as its smooth limit. The "without Higgs" claim attaches to the Comment paper arXiv:2508.17338 (= `perez_sanchez_2024_comment`), which is correctly cited as such in its own bibitem (line 1656) but which doesn't share the framing-claim with arXiv:2401.03705.

This was noted in the calibration context ("`perez_sanchez2024` misbundle... corrected"), but the misbundle is still present in this paper's bibitem annotation and in the §II.A in-text use. Both papers should be cited but with separate claim attribution: arXiv:2401.03705 → Yang-Mills + Higgs, arXiv:2508.17338 → Yang-Mills without Higgs (correcting the Marcolli–vS continuum limit). **Note:** the CLAUDE.md session header lists `perez_sanchez2024` as "misbundle corrected" — if the correction was applied corpus-wide elsewhere but did NOT propagate into Paper 2 §II.A's bibitem annotation, this is the surviving manifestation. **MEDIUM**.

**No CITE-MISATTRIBUTED, CITE-CANT-FIND found.** All other citations are real and authored as claimed.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|------------------|----------|----------|------------------|----------------|
| "no published framework of this kind predicts a physical coupling constant" | §II.A | I cross-checked via WH1 history (CLAUDE.md), and verified via the Marcolli–vS, Perez-Sánchez 2024, Perez-Sánchez 2025 abstracts. None claims a coupling-constant prediction; they derive Wilson lattice / Yang-Mills(/-Higgs) at continuum limit. | **No prior art found for "predicts $\alpha$ from a discrete spectral action."** | Acceptable as stated. The phrasing avoids "first" — it says "no published framework of this kind." This is appropriately bounded relative to what a search can confirm. |
| "improves precision by a factor of 500" (vs paper2_old) | §I last sentence + §VIII.B Table VIII | Internal numerics | Yes — actual ratio 704× | Update to 700× (low severity, no novelty implication) |

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---------|------|---------|----------|
| F1 | `perez_sanchez_2024` (arXiv:2401.03705) annotated and cited as supporting "Yang-Mills without Higgs"; abstract actually says Yang-Mills + Higgs. Claim belongs to `perez_sanchez_2024_comment` (arXiv:2508.17338). Bibitem annotation + §II.A in-text use both need split. | B | CITE-DOESNT-SUPPORT | **MEDIUM** |
| F2 | Abstract: "very unlikely to be accidental" — slightly stronger than §V's body hedge ("does not constitute a derivation") | A | C overstatement | **MEDIUM** |
| F3 | §VIII open question 3: $\pi^3\alpha^3$ residual "to within 0.25%" — depends on CODATA edition (CODATA 2018 = +0.25%, CODATA 2022 = −0.5%); paper doesn't name the edition | A | C overstatement (mild) | **LOW** |
| F4 | §I + Table VIII: "500× improvement" should be 700× (paper understates) | A | C overstatement (in safe direction) | **LOW** |
| F5 | §II.A: "the *prediction* is this paper's contribution" — calling it a "prediction" is mildly inconsistent with the rest of the paper's "empirical observation" framing | A | C overstatement | **LOW** |

No HIGH-severity findings. No mathematical errors (E). No false-novelty (with-found-prior-art) findings. No CITE-MISATTRIBUTED or CITE-CANT-FIND issues. The conjectural-status of the combination rule is preserved at every load-bearing location.

## Broadcast readiness: **YELLOW**

The paper's core load-bearing numerical claims (B=42, F=$\pi^2/6$, $\Delta=1/40$, K=137.036064, $1/\alpha=137.036011$, rel err $8.8\times 10^{-8}$, triple selection at m=3, search count $1.92\times 10^9$, p-value $5.2\times 10^{-9}$) all reproduce externally to the precision claimed. The conjectural status of the combination rule is preserved through abstract / §I / §IV (Theorem 1, §closure) / §VII (Link 3) / §IX (Conclusion) — the §13.5 hard prohibition is honored. There are no math errors, no fabricated citations, no CITE-MISATTRIBUTED entries, and no false-novelty claims.

One MEDIUM-severity defect (the `perez_sanchez_2024` citation says the opposite of what the bibitem annotation and §II.A in-text claim attribute to it — the "Yang-Mills without Higgs" finding belongs to the Comment paper arXiv:2508.17338, NOT arXiv:2401.03705) and one MEDIUM-severity overstatement (the abstract's "very unlikely to be accidental" vs §V's correct hedge) keep this YELLOW rather than GREEN. The two MEDIUM fixes are sprint-scale: one bibitem annotation + one in-text split (F1), one abstract softening (F2). Three LOW-severity items can batch into the final cleanup pass.

## What I could NOT verify (hand to a human expert)

- The interpretive claim that the Hopf base $S^2$ "is the photon momentum sphere" and the fiber $S^1$ "is the U(1) gauge phase of the electromagnetic potential" (§II). This is a physical-mechanism reading the paper offers; it is plausible and consistent with the standard $U(1)$ → Hopf-bundle identification, but Paper 2 does not derive it from a Lagrangian or a quantized field-theoretic construction. A domain expert in NCG / spectral-action quantum field theory would be the right judge.
- The claim (Eq. 27, §IV Obstruction 1) that $|\lambda_{m-1}|\cdot g_{m-1}^{\rm Weyl} = 42$ holds **only** at m=3. I verified the m=3 case ($7/2 \cdot 12 = 42$) but did not exhaust other integers symbolically; the algebraic statement "at no other integer m" should be checkable but I did not run the full integer-uniqueness proof for this particular identity.
- The Theorem 1 ("Three Homes, structural mystery") statement that "no spectral construction on $S^3$ simultaneously generates all three" is a strong negative claim. It is supported by the nine-mechanism elimination over Phases 4B–4I, but those are GeoVac-internal sprints; an NCG / heat-kernel expert may know of a construction not covered in the Phase 4B–4I enumeration. This is the conjectural rule restated as a theorem about lack-of-derivation; the conjectural-status hedge in §IV closure rightly does not let the theorem upgrade $K$ from observation to derivation.

---

**Returned to dispatcher:**
- Title: Paper 2 — *The Fine Structure Constant from Spectral Geometry of the Hopf Fibration*
- Verdict: **YELLOW**
- Pass A counts: A=15, B=2, C=4, D=0, E=0 (some claims have dual A/B or A/C labels — counted once at the dominant verdict)
- Pass B counts: CITE-OK 16, CITE-WRONG-METADATA 0, CITE-MISATTRIBUTED 0, CITE-DOESNT-SUPPORT 1, CITE-CANT-FIND 0
- Severity: HIGH 0, MEDIUM 2, LOW 3
- Top finding: The `perez_sanchez_2024` citation (arXiv:2401.03705) is attributed in the bibitem annotation and §II.A as supporting "Yang-Mills without Higgs," but that abstract actually derives Yang-Mills-**with**-Higgs; the "without Higgs" claim belongs to the Comment paper (`perez_sanchez_2024_comment` = arXiv:2508.17338) and must be split out.
