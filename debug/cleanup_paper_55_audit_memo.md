# Paper 55 Pre-Outreach Audit Memo

**Date:** 2026-06-08
**Auditor scope:** Adversarial pre-outreach audit of `papers/group3_foundations/paper_55_periods_of_geovac.tex` for the Brown (Oxford/IHÉS) + Kleinschmidt (AEI Potsdam) outreach.
**Compile state at audit:** 38 pages, three-pass clean WITH 4 undefined cites + 1 undefined ref. Tests pass (45 passed, 1 skipped, 2.13s wall).

---

## 1. BIBITEM AUDIT

| Bibkey | Claimed source | Verification status | Correction needed |
|:-------|:---------------|:-------------------:|:------------------|
| `bost_connes1995` | Bost–Connes, *Selecta Math (N.S.)* 1, 411–457 (1995) | **VERIFIED** | None — title, authors, journal, volume, pages, year all match. |
| `connes_consani2014_scaling` | Connes–Consani, *C. R. Acad. Sci. Paris Ser. I* 354, 1–6 (2016); arXiv:1603.03191 (2016) | **WRONG (arXiv mismatch)** | Bibkey labels year 2014 but text says 2016 (minor). MAJOR: arXiv:1603.03191 is the LATER *Selecta Math.* 23 (2017) "Geometry of the scaling site" paper. The correct arXiv for the C. R. Math. Acad. Sci. Paris note is **arXiv:1507.05818**. Also journal name should be "C. R. Math. Acad. Sci. Paris" not "C. R. Acad. Sci. Paris". |
| `connes_consani2016_geometry` | Connes–Consani, "Geometry of the scaling site", *Selecta Math (NS)* 23, 1803–1850 (2017) | **VERIFIED** | None on metadata. Note that arXiv:1603.03191 belongs HERE, not in `connes_consani2014_scaling`. Add arXiv:1603.03191 to this bibitem. |
| `connes_consani_moscovici2025_zeta` | Connes–Consani–Moscovici, "Zeta spectral triples and the explicit formulas", arXiv:2511.22755 (2025) | **PARTIAL** | arXiv ID, year, authors confirmed. **Actual arXiv title is "Zeta Spectral Triples"** (no subtitle). Remove "and the explicit formulas" from bibitem title. The PI verification flag note is good as-is. |
| `greenfield_marcolli_teh2014` | Greenfield–Marcolli–Teh, "Twisted spectral triples and quantum statistical mechanical systems", *p-Adic Numbers, Ultrametric Analysis, and Applications* 6, 81–104 (2014); arXiv:1305.5492 (2013) | **VERIFIED** | None. Published title differs from arXiv title ("Type III σ-spectral triples…") but the published title is what's cited. Volume 6, pp 81–104, 2014 verified. |
| `loutey_paper2` | … `loutey_paper51` | **VERIFIED (internal)** | All internal corpus citations. Out of scope for adversarial-reader concern; consistent with paper-list. |
| `fock1935` | V. Fock, "Zur Theorie des Wasserstoffatoms", *Z. Physik* 98, 145–154 (1935) | **VERIFIED** | Classic. None. |
| `camporesi_higuchi1996` | Camporesi–Higuchi, *J. Geom. Phys.* 20, 1–18 (1996) | **VERIFIED** | All metadata confirmed. Add arXiv:gr-qc/9505009 for arXiv-side reader convenience (optional). |
| `kontsevich_zagier2001` | Kontsevich–Zagier, "Periods", in *Mathematics Unlimited — 2001 and Beyond*, Engquist & Schmid eds., Springer (2001) pp. 771–808 | **VERIFIED** | None. |
| `brown2012` | F. Brown, "Mixed Tate motives over Z", *Ann. of Math. (2)* 175, 949–976 (2012) | **VERIFIED** | None. arXiv:1102.1312 available; consider adding. |
| `deligne2010` | P. Deligne, "Le groupe fondamental unipotent motivique de $\mathbb{G}_m - \mu_N$, pour $N = 2, 3, 4, 6$ ou $8$", *Publ. Math. IHES* 112, 101–141 (2010); preprint arXiv:math/0302267 (2005) | **WRONG (arXiv mismatch)** | Journal/title/year/page metadata CORRECT. arXiv:math/0302267 is **Deligne–Goncharov 2003 "Groupes fondamentaux motiviques de Tate mixte"** — a DIFFERENT paper. Deligne's 2010 IHES paper has NO standalone arXiv preprint (DOI 10.1007/s10240-010-0027-6, accessible via NUMDAM/Springer). **Drop the arXiv ID OR move it to a separate bibitem citing Deligne–Goncharov 2005.** This is the single most embarrassing error for a Brown-lineage outreach — Brown cites both Deligne 2010 and Deligne–Goncharov constantly and will spot this immediately. |
| `glanois2015` | C. Glanois, "Motivic unipotent fundamental groupoid of $\mathbb{G}_m \setminus \mu_N$ for $N = 2, 3, 4, 6, 8$ and Galois descents", arXiv:1411.4947 v2 (2015); pub. *J. Number Theory* 182, 36–90 (2018) | **VERIFIED** | None. |
| `fathizadeh_marcolli2016` | F. Fathizadeh, M. Marcolli, "Periods and motives in the spectral action of Robertson–Walker spacetimes", *Comm. Math. Phys.* 356, 641–671 (2017); arXiv:1611.01815 (2016) | **VERIFIED** | All metadata confirmed including 2017 publication after 2016 submission. |
| `marcolli_tabuada2014a` | Marcolli–Tabuada, "Noncommutative motives, numerical equivalence, and semi-simplicity", *Amer. J. Math.* 136, 59–75 (2014); arXiv:1110.2438 (2011) | **WRONG (arXiv mismatch)** | Title/journal/pages/year CORRECT. But arXiv:1110.2438 is the DIFFERENT Marcolli–Tabuada paper "Noncommutative numerical motives, Tannakian structures, and motivic Galois groups". The correct arXiv ID for the *Amer. J. Math.* 136 paper is **arXiv:1105.2950** (submitted May 2011). |
| `karamata1962` | Cited at line 2218 (Tauberian rate uniformity); body claims "Karamata 1962" | **UNVERIFIABLE / MISSING BIBITEM** | **No bibitem in bibliography.** natbib warning at compile. Karamata's foundational Tauberian work is 1937/1938 (Acta Sci. Ind. 450; Hansischen Univ. 12). The "1962" date is suspicious — possibly a confused reference to a later proceedings or Bingham–Goldie–Teugels-style recompilation. **Either add the actual bibitem, replace with Korevaar 2004 reference (already in mind), or remove the citation in favour of Korevaar+Tenenbaum already cited.** |
| `korevaar2004` | Cited at line 2218 (Tauberian rate uniformity) | **MISSING BIBITEM** | natbib warning. The reference exists (Jacob Korevaar, *Tauberian Theory: A Century of Developments*, Grundlehren 329, Springer 2004) but **no bibitem in bibliography**. Add it. |
| `tenenbaum2015` | Cited at line 2218 (Tauberian rate uniformity) | **MISSING BIBITEM** | natbib warning. Reference exists (Gérald Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*; the 3rd ed. is 2015, AMS GSM 163) but **no bibitem**. Add it. |

**Bibitem audit summary:** 11 of 18 cited published-author bibitems are VERIFIED clean. **4 have arXiv-ID mismatches (`deligne2010`, `marcolli_tabuada2014a`, `connes_consani2014_scaling`, `connes_consani_moscovici2025_zeta`); 3 are entirely missing from the bibliography (`karamata1962`, `korevaar2004`, `tenenbaum2015`).** The `deligne2010` error is the most damaging for the targeted readers — Brown will recognise both Deligne 2010 IHES and Deligne–Goncharov 2003 instantly.

---

## 2. SCOPE-HONESTY AUDIT

| Claim | Stated label | Audit | Verdict |
|:------|:-------------|:------|:--------|
| **Theorem 3.1 (M1 pure-Tate, line 344)** | Theorem | Statement is precisely "$M_1 \subset \Q[\pi, \pi^{-1}]$, depth 0, Tate weight $\in \Z$" — verified at 5 corpus witnesses. The proof is short and rests on Kontsevich–Zagier. The "localisation at $\pi^{-1}$" is honest about the formal extension. | **GREEN — scope-honest.** |
| **Theorem 4.1 (M2 pure-Tate, line 418)** | Theorem | The proof sketch is via specialisation of Fathizadeh–Marcolli Thm 6.2 along $\{1\}\times\{0\}^{2n}$, with the Grothendieck class explicitly computed via F–M Thm 7.6 + Lemma 7.1.5 (factorisation of $[Z_{1,2n}]$ as $[\mathbb{P}^n](1+\mathbb{L}^n)$). The pure-Tate refinement over the F–M generic mixed-Tate result is the substantive new content. | **GREEN — proof sketch is honest, refinement is well-flagged.** |
| **Remark 4.2 (Distinction from M3 cyclotomic-level-4)** | Remark | The $\Q(i)$ Witt-splitting (M2) vs $\Q(\zeta_4)$ cyclotomic-conductor (M3) distinction is precise and a legitimate refinement — flagging the deeper unification question as Q5'. | **GREEN.** |
| **Theorem 5.1 (M3 cyclotomic-mixed-Tate level $\le 4$)** | Theorem | Three sub-sectors: classical Brown MZV at $s$ even/odd; vertex-restricted $\chi_{-4}$ via Hurwitz at $1/4,3/4$; depth-$k$ loop tower. Proof sketch is rigorous up to inheriting Deligne 2010 + Glanois 2015 results. | **GREEN — internally consistent;** see bibitem audit above on the `deligne2010` arXiv error. |
| **Proposition 5.2 (Vertex parity is Deligne–Glanois descent)** | Proposition | Claim is the structural identification of the GeoVac vertex-parity-as-$\chi_{-4}$ with the Glanois 2015 level-4-to-level-2 descent operation. The argument is well-laid-out but rests entirely on Glanois Cor 1.1–1.2. **No new motivic conjecture needed.** | **GREEN — scope-honest.** |
| **Theorems 5.3–5.5 (M3 on $S^5$)** | Theorems + Corollary | The three-term polynomial structure $Z(s)$ is derived; Bit-exact verification at $s \in \{6,7,8,9,10\}$ documented. Honest level-$\le 4$ closure with $5/4, 7/4$ shifts reducing to $1/4, 3/4$. | **GREEN.** Empirical $S^{(5)}_{\min}$ PSLQ identification correctly flagged as 2–3 week sprint follow-on. |
| **Theorems 5.6–5.8 (NA-1 trace-functional collapse)** | Theorems + Remarks | Off-diagonal-blind property of the trace, sharpening reading A/B disambiguation. Honest scope on what the JLO probe shows (Reading A POSITIVE) vs what it forecloses (substrate enrichment cannot recover A/B). | **GREEN — explicitly honest scope statements.** |
| **Proposition 7.3 (Thermal pure-Tate)** | Proposition | Closed-form coefficients $a_k^{4D,scalar} = 2\pi^2 \beta/k!$ etc.; factorisation $a_k^{4D} = \beta \cdot a_k^{S^3}$. Stefan–Boltzmann $\pi$-injection explicitly attributed to Mellin (not SD) — orthogonal mechanism flagged. | **GREEN — sub-mechanism distinction is honest.** |
| **Proposition 7.4 ($S^5$ M2 pure-Tate)** | Proposition | Three-term Dirac SD exact + scalar SD closed-form $a_k^\Delta = (6-k)\cdot 4^{k-1}\cdot 2/(3\cdot k!)\cdot\pi^3$. Even-weight ($S^3$, $\pi^{2k}\Q$) vs odd-weight slice ($S^5$, $\pi^3\Q$) is genuinely new dimension-parity sharpening. | **GREEN.** |
| **Inner-factor closure (§7.5)** | "closed in negative" | Combines Paper 18 §IV.6 $\eta$-trivialization + AC factorization + Sprint Yukawa-PSLQ negative. Conclusion: inner-factor sixth tier is categorically disjoint from M1/M2/M3. | **GREEN — closure verdict is honest;** explicitly says "no new sub-mechanism required" rather than overclaiming. |
| **Q5' Stage 1/Stage 2 narrative (§7.6, lines 1787–3012)** | "DEFLATED" / "VIABLE" / open multi-year | The narrative is approximately **1,225 lines of text** describing a sequence of ~30 sub-sprints (TC-1a..1f, PS-1..4, etc.), each with bit-exact zero-residual counts. The overall verdict is honest: Stage 1 closed at the verified-precondition level; Stage 2 (Tannakian closure proper, $U^*$ identification) explicitly multi-year. | **YELLOW — scope-honest but BLOATED.** The 1,225-line single subsection narrating 30 sub-sprints is the antithesis of Brown's compact-and-careful style. For outreach, this should be compressed to ~1.5–2 pages with the granular sprint narrative moved to a long appendix or factored into Paper 56 (which already exists). See readiness verdict. |
| **A1-Matsubara closure (§7.2 thermal sub-case)** | "closed in the negative" | Proposition 7.3 is the substantive claim. Body explicitly says Matsubara mode sum lives in IR / large-$t$, invisible to SD. | **GREEN.** |
| **A6 closure (M3 on $S^5$)** | "closed POSITIVE" | Corollary 5.10 + proof sketch + bit-exact verification at $s\in\{6..10\}$ matching $\le 10^{-12}$ down to $10^{-33}$. | **GREEN.** |
| **A7 verdict (Q5' deeper M2/M3 unification)** | "HALF-STRUCTURAL" | "Shared field at number-field level, no shared Tannakian symmetry known". This is appropriately cautious. | **GREEN.** |
| **Appendix A (BC / scaling-site)** | "structurally orthogonal" | Records GeoVac ≠ BC substrate at four-axis orthogonality; explicitly does NOT claim a BC-RH path for GeoVac. | **GREEN — exemplary scope honesty for the parallel arc.** |

**Scope-honesty summary:** The substantive theorems (3.1, 4.1, 5.1, 5.3–5.5, 5.6–5.8, 7.3, 7.4) are all scope-honest. The Q5' open-questions section (§7.6) is honest but bloated — see readiness verdict.

---

## 3. §13.4a EQUATION VERIFICATION

| Theorem / equation | Verification status | Test file (if any) | Adequacy |
|:-------------------|:--------------------|:-------------------|:---------|
| **Thm 3.1 (M1 pure-Tate)** | VERIFIED | `tests/test_paper55_m1_pure_tate.py` (316 lines) — PSLQ at 80 dps against $\{\pi^n\}_{n=-4..4}$ for $W_1 = \pi$, $W_2 = 4/\pi$, $W_4 = 1/\pi^2$. | **ADEQUATE.** But the test file is NOT referenced from Paper 55 — only `test_paper55_M2_mixed_tate.py` is cited at line 623. |
| **Thm 4.1 (M2 pure-Tate)** | VERIFIED | `tests/test_paper55_M2_mixed_tate.py` (563 lines) — referenced in Remark 4.3 at line 623. | **ADEQUATE.** |
| **Thm 5.1 (M3 cyclotomic mixed-Tate)** | VERIFIED | `tests/test_paper55_m3_cyclotomic_mixed_tate.py` (497 lines) | **ADEQUATE.** Not referenced from Paper 55. |
| **Thm 5.3 (Parity discriminant on $S^5$)** | DRIVER ONLY | `debug/sprint_a6_m3_s5_derivation.py` referenced at line 1012 | **WEAK** — driver in `debug/`, not in `tests/`. For pre-outreach, promote the $s\in\{6..10\}$ bit-exact verification to a regression test. |
| **Thm 5.4 ($\chi_{-4}$ identity on $S^5$)** | DRIVER ONLY | Same driver as 5.3 | Same as 5.3. |
| **Cor 5.5 (Cyclotomic mixed-Tate on $S^5$)** | INHERITS | Inherits Glanois 2015 + Hurwitz shift identity; no direct numerical test. | **STRUCTURAL ARGUMENT — adequate.** |
| **Thm 5.6 (NA-1 diagonal-collapse)** | VERIFIED | Sprint memo `debug/sprint_na1_depth2_mellin_memo.md` cites bit-exact at $s_{tot}\in\{3..7\}$, 200 dps, data file `debug/data/na1_depth2_mellin_results.json` | **ADEQUATE if data file exists.** Not promoted to regression test. |
| **Thm 5.7 (Closed form for $M_3^{\gamma_P}$)** | VERIFIED | PSLQ at 50/100/200 dps in `debug/sprint_na1_depth2_mellin_compute.py` | **ADEQUATE** as data; not in `tests/`. |
| **Thm 5.8 (NA-1-offdiag trace-functional collapse)** | PROOF-ONLY | Theorem with bit-exact proof in text; structural argument (proof at line 1201–1216 of Paper 55) | **ADEQUATE — proof is the verification.** |
| **Thm 5.9 (JLO depth-2 $S_3$-cocommutativity)** | VERIFIED | `debug/sprint_jlo_depth2_compute.py` bit-exact in sympy rational | **ADEQUATE.** Not promoted to regression test. |
| **Prop 7.3 (Thermal pure-Tate)** | STRUCTURAL ARGUMENT | Inherits Paper 51 Thm 3.1 + Cor 2.1 via $\beta$-factorisation. Algebraic identity, no separate numerical test required. | **ADEQUATE.** |
| **Prop 7.4 ($S^5$ M2 pure-Tate)** | DRIVER ONLY | `debug/sprint_a2_s5_sd_coefficients.py` cited at line 1711 | **WEAK** — promote to regression test (parallel to `test_paper55_M2_mixed_tate.py`). |
| **Eq. 4.4 (Grothendieck-class explicit form)** | DRIVER ONLY | `debug/sprint_a8_grothendieck_class_memo.md` "verified symbolically at $n=1,\ldots,5$" | **WEAK** — promote to regression test. |

**§13.4a verification summary:**
- **3 of 12 theorems/key propositions** are backed by tests in `tests/` (M1, M2, M3 Paper 55 tests).
- **6 of 12** rely on `debug/` drivers without regression promotion.
- **3 of 12** are structural arguments (proof IS the verification).

**Gap-closing recommendation for pre-outreach:** add 3 regression-test files:
1. `tests/test_paper55_s5_m2.py` — verifies Prop 7.4 ($S^5$ Dirac three-term + scalar closed form, with $a_6^\Delta = 0$ check).
2. `tests/test_paper55_s5_m3.py` — verifies Thms 5.3–5.4 at $s\in\{6..10\}$.
3. `tests/test_paper55_grothendieck.py` — verifies Eq. 4.4 symbolically at $n=1,\ldots,5$.

Effort estimate: ~3–4 hours total (drivers exist; promotion is mechanical).

---

## 4. LATEX CLEANLINESS

**Three-pass compile:** SUCCESS, 38 pages, 811050 bytes. **Five blocking issues remain:**

1. **Undefined reference `sec:open_m2_m3`** (lines 3038, 3069). Actual label is `subsec:open_m2_m3` at line 1788. **One-line patch:** change `\S\ref{sec:open_m2_m3}` → `\S\ref{subsec:open_m2_m3}` at both call sites.
2. **Missing bibitem `karamata1962`** (line 2218). See §1.
3. **Missing bibitem `korevaar2004`** (line 2218). See §1.
4. **Missing bibitem `tenenbaum2015`** (line 2218). See §1.
5. **arXiv ID errors** on `deligne2010`, `marcolli_tabuada2014a`, `connes_consani2014_scaling`, `connes_consani_moscovici2025_zeta` title. See §1.

**Other LaTeX issues (cosmetic, non-blocking):**
- 9 overfull `\hbox` warnings (lines 107–115, 166, 196–205, 241–244, 379–390, 760–773, 851–866, 1228–1248, 1253–1271, 1528–1537). Most are minor (< 50 pt); line 1528–1537 is 118 pt — substantial but visually tolerable. Mitigation: gentle rewording, particularly around the Mellin-engine formulas with hard line breaks.
- ~50 `Package hyperref Warning: Token not allowed in a PDF string (Unicode)` warnings — caused by `\mathbb{...}` and similar math inside `\section`/`\subsection` titles which hyperref tries to flatten into PDF bookmark strings. Cosmetic only, does not affect rendering. Mitigation: `\texorpdfstring{\mathbb{Z}[i,1/2]}{Z[i,1/2]}` etc. in affected section titles if a clean log is wanted.

**Bottom-line LaTeX:** compile is clean PDF output; warnings are surface-level fixes; only 5 items above are blocking pre-submission.

---

## 5. CROSS-PAPER CONSISTENCY

| Cross-paper claim | Aligned? | Notes |
|:------------------|:--------:|:------|
| Paper 32 §VIII `thm:pi_source_case_exhaustion` ↔ Paper 55 `thm:case_exhaustion` (line 304) | ✓ | Paper 55 uses a local label; attribution to Paper 32 §VIII is explicit. |
| Paper 32 §VIII Corollaries `cor:m1_pure_tate`, `cor:m2_mixed_tate`, `cor:m3_cyclotomic_mixed_tate` ↔ Paper 55 Theorems 3.1, 4.1, 5.1 | ✓ | Paper 32 has Corollaries; Paper 55 has full Theorems with proof sketches and witness panels. Standard pattern (Paper 55 is the "standalone with proofs" version of Paper 32 §VIII content). |
| Paper 18 §III.7 master Mellin engine M1/M2/M3 ↔ Paper 55 sections 3/4/5 | ✓ | Identical M1/M2/M3 trinity. |
| Paper 18 §IV.6 `thm:eta_trivialization`, `thm:ac_factorization` (lines 2073, 2097) ↔ Paper 55 §7.5 inner-factor closure | ✓ | Paper 55 cites both labels explicitly and uses them correctly. |
| Paper 28 Theorem 1 ($\zeta_{D^2}(s)$ closed forms) ↔ Paper 55 §4 M2 witnesses table (lines 636–644) | ✓ | Numerical agreement spot-checked: $\zeta_{D^2}(2) = \pi^2 - \pi^4/12$, $\zeta_{D^2}(3) = \pi^4/3 - \pi^6/30$, etc. |
| Paper 28 Theorem 3 ($D_{\text{even}}-D_{\text{odd}}$ identity) ↔ Paper 55 §5.4 lines 814–820 | ✓ | $D_{\text{even}}(4), D_{\text{odd}}(4)$ values $\frac{\pi^2}{2} - \frac{\pi^4}{24} \mp 4G \pm 4\beta(4)$ confirmed by Paper 28 references. |
| Paper 50 F-theorem values $F_s, F_D$ ↔ Paper 55 Eqs. 6.4–6.5 (lines 1461–1463) | ✓ | Numerical match: $F_s = \log 2/8 - 3\zeta(3)/(16\pi^2)$, $F_D = \log 2/4 + 3\zeta(3)/(8\pi^2)$. Algebra of $F_D \pm 2F_s$ checks out. |
| Paper 51 Cor 2.1 / Thm 3.1 (Dirac and scalar SD on $S^3$) ↔ Paper 55 §4 witness panel | ✓ | $a_0^{D^2} = 4\pi^2$, $a_1^{D^2} = -2\pi^2$, $a_k^\Delta = 2\pi^2/k!$ — consistent across papers. |
| Paper 38 universal $4/\pi$ rate ↔ Paper 55 Eq. 6.6 line 1507 | ✓ | $4/\pi = \text{Vol}(S^2)/\pi^2$ — same M1 reading. |
| Paper 56 `thm:injection_g4` ↔ Paper 55 Remark 5.13 (line 1316) | ✓ | Paper 56 contains TC-1a..1f (which paper 55 §7.6 narrates inline). **DUPLICATION RISK:** Paper 55 §7.6 lines 2797–2998 narrate TC-1a..1f essentially verbatim from Paper 56 §sec:tc1a..1f. For outreach, recommend Paper 55 §7.6 cite Paper 56 instead of in-line narration. |
| Paper 18 §III.7 sees M2 as "Mixed-Tate strict-subset pure-Tate" ↔ Paper 55 §4 same claim | ✓ | Consistent. |
| MEMORY rule "Discrete-for-skeleton — exact arithmetic NOT PSLQ + N_max sweeps" ↔ Paper 55 verification methodology | ✓ | Paper 55 verifications are bit-exact sympy.Rational for skeleton-side claims; PSLQ only applied to Layer-2 / period-output identification. Aligned with discipline. |

**Cross-paper consistency summary:** No contradictions. One DUPLICATION issue: Paper 55 §7.6 narrates ~200 lines of Paper 56's TC-1a..1f content inline. For outreach, this should cite Paper 56 instead.

---

## 6. READINESS VERDICT

**Overall: YELLOW** — paper is structurally sound and scope-honest, but five concrete issues must be patched before sending to Brown / Kleinschmidt. Estimated total patch effort: **4–6 hours**.

### Blocking issues (must fix)

1. **B1 (15 min).** `deligne2010` arXiv ID error: remove "preprint arXiv:math/0302267 (2005)" from this bibitem. The IHES paper has no arXiv preprint; that arXiv ID belongs to Deligne–Goncharov 2003. **This is the single most critical issue for the Brown outreach** — Brown will recognise both papers instantly and the error will dent credibility immediately.

2. **B2 (10 min).** `marcolli_tabuada2014a` arXiv ID error: change "arXiv:1110.2438 (2011)" → "arXiv:1105.2950 (2011)". (1110.2438 belongs to the different Marcolli–Tabuada paper on Tannakian structures.)

3. **B3 (10 min).** `connes_consani2014_scaling` mix-up: split the arXiv ID — change "preprint arXiv:1603.03191 (2016)" → "preprint arXiv:1507.05818 (2015)". Add arXiv:1603.03191 to `connes_consani2016_geometry` instead.

4. **B4 (5 min).** Trim title in `connes_consani_moscovici2025_zeta`: remove "and the explicit formulas" subtitle (actual arXiv title is just "Zeta Spectral Triples").

5. **B5 (15 min).** Add 3 missing bibitems for the Tauberian sprint paragraph:
   - `karamata1962`: replace with Karamata 1937 (Acta Sci. Ind. 450) OR remove this citation (Korevaar 2004 + Tenenbaum 2015 already redundant).
   - `korevaar2004`: J. Korevaar, *Tauberian Theory: A Century of Developments*, Grundlehren 329, Springer (2004).
   - `tenenbaum2015`: G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, 3rd ed., AMS GSM 163 (2015).

6. **B6 (1 min).** Fix undefined ref `\S\ref{sec:open_m2_m3}` → `\S\ref{subsec:open_m2_m3}` at lines 3038 and 3069.

### Recommended (non-blocking but advisable)

7. **R1 (1–2 hours).** **Compress §7.6 Q5' narrative** from ~1225 lines to ~2 pages. Move the granular sprint chronicle (TC-1a..1f, PS-1..4, ~30 sub-sprints) to either an appendix or a footnote pointer to Paper 56 (which already houses TC-1a..1f at theorem-grade). For Brown/Kleinschmidt outreach, the current narrative reads more like an internal lab notebook than a finished paper — it will be the single weakest impression made. **The substantive content (Q5' structural verdict, multi-year roadmap, viable-lineage statement) can fit in 2 pages without losing any honest scope.**

8. **R2 (3–4 hours).** Promote 3 `debug/` drivers to regression tests:
   - `tests/test_paper55_s5_m2.py` (Prop 7.4).
   - `tests/test_paper55_s5_m3.py` (Thms 5.3–5.4).
   - `tests/test_paper55_grothendieck.py` (Eq. 4.4 symbolic at $n=1,\ldots,5$).
   
   Closes §13.4a verification gap for the substantive new $S^5$ content.

9. **R3 (15 min).** Add cross-references from Paper 55 text to `tests/test_paper55_m1_pure_tate.py` and `tests/test_paper55_m3_cyclotomic_mixed_tate.py` (parallel to the M2 reference at line 623), so the verification audit trail is visible from the paper.

10. **R4 (cosmetic, 30 min).** Address 1–2 worst overfull `\hbox` warnings (especially line 1528–1537 at 118 pt) with mild rewording.

### Not blocking

- The 50 hyperref Unicode warnings (cosmetic).
- Most overfull \hbox warnings (< 50 pt each).
- The Q5' multi-year frontier framing — the honest scope is intact, just verbose.

**Recommendation:** Patches B1–B6 (~1 hour) are mandatory pre-outreach. R1 (Q5' compression, ~2 hours) is strongly recommended — the current §7.6 is the weakest part of the paper for the targeted audience. R2 (3–4 hours regression tests) and R3 (15 min cross-refs) close the §13.4a verification audit trail and are good hygiene. **Total effort to GREEN: ~6 hours.** Without R1, the paper is technically sound but stylistically unfit for Brown-tier outreach.

---

**End of memo.**
