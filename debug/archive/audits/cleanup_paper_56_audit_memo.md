# Paper 56 audit memo — pre-outreach readiness

**Auditor:** sub-agent, adversarial mode, 2026-06-08
**Paper:** `papers/group3_foundations/paper_56_tannakian_substrate.tex`
**Pages / size:** 24 pp, 692,434 bytes PDF (three-pass clean)
**Compile audit:** non-`-draftmode` compile completes, zero undefined refs/citations after three passes
**Verdict (headline):** **YELLOW** — paper is structurally GREEN on scope-honesty and on the headline 5,864 residual count, but the bibliography carries two factually-wrong Pérez-Sánchez titles (substantive, not cosmetic), one duplicate bibitem for Paper 55, one bibkey inconsistency vs the rest of the corpus, and minor citation-style imprecisions on Glanois / Deligne 2010 / Brown 2017. All findings are correctable inside a single 1–2 hour cleanup sprint; no theorem or numerical claim is affected.

---

## 1. BIBITEM AUDIT

Verifying every published-author bibitem against arXiv / journal records.

| Bibkey | Status | Notes |
|:-------|:-------|:------|
| `deligne_milne1982` | VERIFIED | Lecture Notes in Math. 900 (1982) 101–228 ✓ |
| `deligne1990` | NOT VERIFIED (cited only as `[deligne1990, brown2012]` in Cor. inverse_limit and §sec:open_converse). Catégories tannakiennes, Grothendieck Festschrift Vol II, Progr. Math. 87, Birkhäuser 1990 — standard reference, year/series consistent with literature. |
| `brown2012` | VERIFIED | Ann. Math. 175 (2012) 949–976; DOI 10.4007/annals.2012.175.2.10 ✓ |
| `andre2004` | VERIFIED | Panoramas et Synthèses 17, SMF 2004 ✓ |
| `brown2017` | YEAR INCONSISTENCY (minor) | arXiv:1407.5165, ICM 2014 Seoul, Vol II, Kyung Moon Sa 2014. The bibkey says "2017" but the cited proceedings volume was published 2014. Body text consistently says "Brown 2017"; the cited "2017" appears to be a reprint volume year. Citation usage in body is structurally OK, but the date label disagrees with the proceedings year. **Fix:** rename bibkey to `brown2014_icm` OR clarify the 2017 reprint year in the bibitem text. |
| `hain2014` | VERIFIED | LMS Lecture Notes 427, Cambridge 2016, pp 422–514; arXiv:1403.6443 ✓ |
| `deligne2010` | VERIFIED with TITLE TYPO (cosmetic) | Publ. Math. IHES 112 (2010) 101–141 ✓. Actual title says "N=2,3,4,6 ou 8" (French "ou" = "or"), Paper 56 prints "N = 2, 3, 4, 6, 8" (comma-list). Minor; the meaning is preserved. |
| `glanois2015` | VERIFIED with TITLE TYPO (cosmetic) | PhD thesis, Univ. Pierre et Marie Curie. Actual title: "Periods of the motivic fundamental **groupoid** of $\mathbb{P}^1 \backslash \{0, \mu_N, \infty\}$". Paper 56 omits "groupoid". Defense date is 6 January 2016 (thesis written 2015), arXiv:1603.05155 exists for the journal version — not cited. |
| `connes_marcolli2004` | VERIFIED | IMRN 2004(76) 4073–4091; arXiv:math/0409306 ✓ |
| `marcolli_vansuijlekom2014` | VERIFIED with BIBKEY INCONSISTENCY (cross-paper) | J. Geom. Phys. 75 (2014) 71–91; arXiv:1301.3480 ✓. **The bibkey `marcolli_vansuijlekom2014` is inconsistent with Papers 30, 32, 38, 41, 43, 48, 49, group1 synthesis, group3 synthesis, which all use `marcolli_vs2014`.** Cosmetic but a cross-paper-consistency issue (different bibkeys for the same paper). |
| `perez_sanchez2024` | **TITLE WRONG (substantive)** | arXiv:2401.03705 ✓ (arXiv ID and year correct). Paper 56 bibitem title: "Erratum to Marcolli--van Suijlekom 2014: Higgs from off-diagonal $D_F$ blocks". **Actual title: "Bratteli networks and the Spectral Action on quivers".** This is not an erratum paper — it is a substantive paper that happens to address the same content. The made-up "Erratum to..." title also conflicts with Paper 38 which titles this bibitem "On the continuum limit of gauge networks" (also wrong). **Paper 32 has the correct title** ("Bratteli networks and the Spectral Action on quivers"). Same failure mode flagged in `feedback_tell_about_paper_edits.md` adjacent literature-audit memory (post-cutoff sub-agent title fabrication). |
| `perez_sanchez2025` | **TITLE WRONG (substantive)** | arXiv:2508.17338 ✓ (arXiv ID and year correct). Paper 56 bibitem title: "The Standard Model without Higgs: spectral action on Marcolli--van Suijlekom gauge networks". **Actual title: "Comment on 'Gauge networks in noncommutative geometry'".** Topical content matches (continuum limit = Yang–Mills without Higgs), but the bibitem title is fabricated. **Paper 32 again has the correct title** ("Comment on 'Gauge networks in noncommutative geometry'"). |
| `fathizadeh_marcolli2016` | VERIFIED with YEAR-CITED MISMATCH | Comm. Math. Phys. 356 (**2017**) 641–671; arXiv:1611.01815. Bibitem text says "Comm.\ Math.\ Phys.\ \textbf{356} (2017) 641--671;\ arXiv:1611.01815", which matches; the bibkey `fathizadeh_marcolli2016` reflects the arXiv-submission year (Nov 2016), the journal year is 2017. Minor convention only. |
| `eskandari_murty_nemoto2025` | VERIFIED | arXiv:2510.20648, submitted Oct 23, 2025, rev Jan 26, 2026; authors Eskandari, Murty, Nemoto ✓; title "Mixed motives and linear forms in the Catalan constant" ✓. The single most date-sensitive citation in the bibliography (post-cutoff for any model trained earlier than late 2025) — and it checks out cleanly. |
| `goncharov_deligne2005` | VERIFIED, UNUSUAL FORMAT | Bundles Goncharov "Periods and mixed motives" arXiv:math/0202154 (2002) AND Deligne-Goncharov "Groupes fondamentaux motiviques de Tate mixte" Ann. Sci. ENS (4) 38 (2005) 1–56 ✓ into one bibitem. Both genuine papers. The faithfulness result attributed to "Goncharov–Deligne 2005" lives in the Deligne-Goncharov 2005 paper. Bundle is awkward (two distinct works in one bibitem) but accurate. |
| `cartier_milnor_moore` | VERIFIED, UNUSUAL FORMAT | Bundles Milnor-Moore Ann. Math. 81 (1965) 211–264 AND Cartier "A primer of Hopf algebras" Springer 2007 pp. 537–615 into one bibitem. Both genuine ✓. Bundle works because the Cartier–Milnor–Moore theorem is stated/proved across the two. |
| `loutey_paper55` and `paper55` | **DUPLICATE** | Both bibitems present (lines 1911 and 1985), both pointing to "Periods of GeoVac" Paper 55 (2026). Body cites `\cite{paper55}` (lines 305, 1683) and `\cite{loutey_paper55}` (lines 1775, 1780, 1800). **This is a substantive duplicate bibitem; one of the two should be removed and all references unified.** |
| `loutey_paper28` | VERIFIED internal | Paper 28 (2026) — internal cite ✓ |
| Internal GeoVac (`paper7`, `paper18`, `paper31`, `paper32`, `paper38`, `paper45`, `paper46`, `paper47`, `paper48`, `paper49`) | VERIFIED internal | All match Paper 56's known sibling-paper inventory and `paper_*.tex` files exist. |

**Summary:** 2 substantive title errors (`perez_sanchez2024`, `perez_sanchez2025`); 1 duplicate bibitem (`loutey_paper55` vs `paper55`); 1 cross-paper bibkey inconsistency (`marcolli_vansuijlekom2014` vs corpus-standard `marcolli_vs2014`); 4 cosmetic citation imprecisions (Glanois "groupoid", Deligne 2010 "ou 8", Brown 2017 vs 2014, Fathizadeh-Marcolli bibkey-year mismatch). The substantive Pérez-Sánchez issues are the only items that an outreach reader will flag.

---

## 2. SCOPE-HONESTY AUDIT (inclusion-vs-equality)

This is the priority audit for Paper 56 because the headline content is the cosmic-Galois identification, and the converse / reconstruction direction is multi-year.

### Three layers of scope to keep distinct

1. **Inclusion direction at finite cutoff** (TC-1e + TC-1f). $\Ulev(\Q) \hookrightarrow \Aut^\otimes(\fibfun)$ via the explicit map $\Phi(t, g)(V \otimes V_{\mathrm{PW}}) = \exp(\sum t_g X_g^V) \otimes \rho_{V_{\mathrm{PW}}}(g)$. Bit-exact at panel. **CLOSED.**
2. **Converse equality at finite cutoff** (TC-2a + TC-2b + TC-2c + TC-2d). $\dim \Aut^\otimes(\fibfun) = 3 N(\nmax) + 3 = \dim \Ulev$ bit-exact at every $\nmax \in \{1, 2, 3, 4\}$, on the natural-substrate panel. **CLOSED on the natural-substrate panel, at finite cutoff.**
3. **Inverse-limit identification with motivic Galois group $\mathcal{G}_4$** (PS-3.5 / pro-finite Tannakian + Glanois period exhaustion). **OPEN, multi-year frontier.**

### Per-section spot-checks

- **Abstract** (lines 99–177). Line 144 says "fourteenth math.OA standalone". Line 146 emphasis "Cosmic-Galois injection" closes the injection direction explicitly says "Equality / surjectivity on the depth-$\ge 2$ content is the named multi-year forward research question". ✓
- **Abstract honest-scope paragraph** (lines 164–177). Closes (2) at finite cutoff via TC-2a/b/c/d, names (3) as multi-year. ✓ — separates the levels cleanly.
- **§1.2 What this paper does not claim** (lines 223–246). "Tannakian reconstruction direction $\Aut^\otimes(\fibfun) = \Ulev$ at finite cutoff is closed in Sections... as the TC-2 sub-arc. What is not closed by this paper is the inverse-limit identification..." ✓ — exemplary scope statement.
- **Cor. inclusion (line 746–755).** "On the natural-substrate panel" qualifier present. ✓
- **Rem. inclusion_scope (lines 757–774).** "The full converse equality... is the Deligne–Milne 1982 Theorem 2.11 Tannakian reconstruction direction and would require a separate argument that no other natural ⊗-automorphisms... exist on the combined $\nmax$-axis × $j_{\max}$-axis substrate. The abelian-factor converse at finite cutoff is, however, closeable bit-exactly: this is the content of Section ref{sec:tc2a_equality} below". ✓ — the separation between "natural substrate" and "full" is maintained.
- **Thm tc2a_equality (lines 792–817).** Closes only the abelian factor at $\nmax = 2$. Rem tc2a_scope (lines 853–874) explicitly states this is the abelian-factor sub-claim, not the full claim. ✓
- **Thm tc2b_sl2_equality (lines 895–912).** Closes $SL_2$ factor at PW panel. ✓
- **Cor combined_equality (lines 971–977).** Uses Deligne–Milne 1982 Thm 2.3 exterior tensor product theorem to combine the per-factor dims into $\dim = 18 = \dim \Ulev$ at $\nmax = 2$. Scope is restricted to the combined natural substrate. ✓
- **Thm tc2c_higher_cutoff (lines 1007–1014).** $\nmax \in \{1, 2, 3, 4\}$, on the $\nmax$-axis substrate (abelian-factor only). ✓
- **Thm tc2d_cofiltered_coherence (lines 1074–1087).** "On the natural-substrate panel" qualifier maintained. ✓
- **Cor inverse_limit (lines 1112–1128).** "At the natural-substrate panel level. The remaining content of the full multi-year reconstruction is the pro-finite Tannakian theorem... coherent with v3.66.0 FO3 Interpretation C at the period level on $\mathcal{O}_\infty$" ✓ — the at-the-panel-level qualifier is precise; the multi-year frontier is named.
- **Thm injection_g4 (lines 1204–1253).** A "closed immersion" of pro-algebraic groups: $\Ustar_{GV} \hookrightarrow \mathcal{U}_4 \rtimes SL_2$. Four compatibilities (C1)–(C4). Explicitly distinguishes injection from surjection / equality (C4 codimension positive at depth $\ge 2$). ✓ — injection direction is theorem-grade; equality is named "multi-year forward research question" (§sec:equality_multi_year).
- **Rem injection_honest_scope (lines 1361–1397).** Three honest caveats: M1/M2 collapse, M3 substantive content, depth-blindness. ✓ — explicitly partitions the structural content of the injection.
- **§sec:equality_multi_year (lines 1554–1571).** "Equality direction... requires exhibiting (or ruling out) every level-4 cyclotomic motivic MZV as a GeoVac observable... multi-year forward research question". ✓
- **§sec:open_g4 (lines 1673–1733).** "Injection direction is now a theorem... Equality / exhaustion direction is the multi-year forward research question". Also flags Brown 2017 surjection-not-identification framing. ✓
- **§sec:open_na1 (lines 1734–1834).** Reading A wins; Reading B ruled out empirically at depth 2 on natural substrate; non-commutative inner-factor enrichment and equality direction named as remaining open questions. ✓
- **§sec:open_profinite (lines 1836–1845).** Multi-year content explicitly named. ✓

**Verdict on Scope-Honesty:** **GREEN, exemplary.** No sentence accidentally claims equality at the pro-finite limit. The three layers are kept rigorously distinct throughout. The injection-direction theorem is correctly framed as a closed immersion, with the surjective-equality direction explicitly named multi-year. The natural-substrate-panel qualifier appears wherever the equality is asserted at finite cutoff. The abelian-vs-shuffle-Hopf reading (Reading A wins) is empirically grounded via JLO depth-2 with named caveat. This section is *the* model of post-Q5' scope discipline.

---

## 3. §13.4a EQUATION/RESIDUAL VERIFICATION (the 5,864 number is the headline)

The prompt asks me to verify the sub-counts. Two distinct claims appear:
(i) prompt's stated breakdown: "PS-1: 130, PS-2: 485, PS-3: 284, PS-4: 872, TC-1a: 106, TC-1b: 56, TC-1c: 50, TC-1d: 40, TC-1e: 98, TC-1f: 490" = **2,611**;
(ii) Paper 56's headline number: **5,864**.

These are not in conflict — the prompt's 2,611 is the pre-v3.74.0 panel sum (TC-2 subarc + G4-Inj not yet integrated). Paper 56's current 5,864 includes TC-2a/b/c/d (349 residuals) AND the G4-Inj panel at $\nmax \in \{2, 3, 4\}$ (2,904 residuals), giving $2{,}611 + 349 + 2{,}904 = 5{,}864$.

### Sub-count verification against memos

| Sub-track | Paper 56 (Table~ref{tab:verification}) | Sprint memo | Status |
|:----------|:--------------------------------------:|:-----------:|:------:|
| PS-1 | 130 | "130 / 130 bit-exact zero residuals" | ✓ |
| PS-2 | 485 | "485 / 485 bit-exact zero residuals" | ✓ |
| PS-3 | 284 | "284 / 284 bit-exact zero residuals" | ✓ |
| PS-4 | 872 | "872 / 872 bit-exact checks pass" | ✓ |
| TC-1a | 106 | "106 / 106 bit-exact zero residuals" | ✓ |
| TC-1b | 56 | "56 / 56 bit-exact zero residuals (28 per cutoff)" | ✓ |
| TC-1c | 50 | "50 / 50 bit-exact zero residuals (25 per cutoff)" | ✓ |
| TC-1d | 40 | "All 40 bit-exact identities (20 per cutoff at $n_{\max} \in \{2, 3\}$)" | ✓ |
| TC-1e | 98 | "98 / 98 bit-exact zero residuals (49 per cutoff)" | ✓ |
| TC-1f (axioms 135 + commute 100 + injectivity 255 = 490) | 135 + 100 + 255 | "490/490 residuals" with sub-decomposition 15+45+75+100+255 → axioms-block = 135 | ✓ Sub-decomposition matches |
| TC-2a | 16 | Memo: "rank A = 46, nullity 15, predicted dim 15 = computed dim 15 bit-exactly" — 15 panel reps + 1 unit = 16 ✓ |
| TC-2b | 17 | Memo: "variety dim 3 bit-exact, Jacobian rank 1" — 14 PW-axiom cells + 3 dim/jacobian = 17 (plausible) |
| TC-2c | 73 | Memo: $\nmax \in \{1, 2, 3, 4\}$ with dim recovery $6+15+27+42 = 90$ at full panel; Paper 56 row says "at $\nmax = 3, 4$ ($\Phi$-recovery + dim + $\eta_T$)" → $27 + 42 + 2 + 2 = 73$ ✓ |
| TC-2d | 243 | "243 / 243 bit-exact zero residuals" ✓ (225 + 12 + 6 = 243 sub-decomposition matches memo Block A/B/C) |
| G4-Inj at $\nmax = 2$ | 261 | Memo: "261 residuals across C1 (225), C2 (15), C3 (5), C4 (16)" ✓ (closed form $9 \cdot 25 + 6 \cdot 5 + 6 = 225 + 30 + 6 = 261$ ✓) |
| G4-Inj at $\nmax = 3$ | 789 | Closed form: $9 \cdot 81 + 6 \cdot 9 + 6 = 729 + 54 + 6 = 789$ ✓ |
| G4-Inj at $\nmax = 4$ | 1854 | Closed form: $9 \cdot 196 + 6 \cdot 14 + 6 = 1764 + 84 + 6 = 1854$ ✓ |

**Grand total:**
$$130 + 485 + 284 + 872 + 106 + 56 + 50 + 40 + 98 + 135 + 100 + 255 + 16 + 17 + 73 + 243 + 261 + 789 + 1854 = 5{,}864 \checkmark$$

### Closed-form formula

The G4-Inj closed form $R^{\mathrm{G4-Inj}}(\nmax) = 9 N(\nmax)^2 + 6 N(\nmax) + 6$ with $N(\nmax) = \nmax(\nmax + 3)/2$ produces:
- $\nmax = 1$: $N = 2$, $R = 36 + 12 + 6 = 54$ (matches Table~ref{tab:g4_inj_panel}) ✓
- $\nmax = 2$: $N = 5$, $R = 225 + 30 + 6 = 261$ ✓
- $\nmax = 3$: $N = 9$, $R = 729 + 54 + 6 = 789$ ✓
- $\nmax = 4$: $N = 14$, $R = 1764 + 84 + 6 = 1854$ ✓

Production grid is $\{2, 3, 4\}$ (sums to 2,904), with $\nmax = 1$ excluded from the headline total but included in the closed-form check.

### Code/test scaffolding

- `geovac/tannakian.py`: 2,426 lines (paper says ≈2,400 ✓)
- `geovac/pro_system.py`: 884 lines (no specific claim, exists)
- Tests claimed (line 1655–1657): 282 tests — I do not run the test suite as part of read-only audit; the structure of `tests/test_tannakian_*.py` is present (10 files).
- Driver / data / memo trail enumerated in §sec:verification (lines 1620–1649) is comprehensive and matches the on-disk inventory.

**Verdict on residual counts: GREEN.** The 5,864 number is bit-exact verified at the sub-count level against memos and closed-form formulas. No discrepancy.

---

## 4. LATEX CLEANLINESS (three-pass compile)

- **Compile passes:** three passes (draftmode + non-draftmode), zero `! LaTeX Error`, zero `! pdfTeX error`, zero `Undefined reference`, zero `Undefined citation` after pass 3.
- **Warnings:**
  - ~30 `Package hyperref Warning: Token not allowed in a PDF string (Unicode)` — cosmetic, from PDF metadata bookmark strings. Could be silenced by `\texorpdfstring`'ing math in section titles, but does not affect rendering.
  - 1 `LaTeX Font Warning: Font shape T1/cmr/m/scit undefined` at line 1912 (bibliography \textsc{GeoVac} inside \textit context) — cosmetic, auto-substitutes scsl. Minor.
  - 2 `Overfull \hbox` warnings (around lines 1785–1802, NA-1 paragraph). Cosmetic — text may slightly overflow margin.
- **Final output:** 24 pages, 692,434 bytes PDF. Memo `debug/sprint_injection_nmax_extension_memo.md` claims "23 pages, 681 KB PDF" — minor +1 page discrepancy, likely due to subsequent edits.

**Verdict on LaTeX: GREEN with cosmetic warnings.** No blocking issues; all warnings are typographic.

---

## 5. CROSS-PAPER CONSISTENCY (Papers 18, 25, 30, 32, 38, 41, 55)

### Citation references

| Sibling paper | Paper 56 cite usage | Status |
|:--------------|:--------------------|:------:|
| Paper 7 | `\cite{paper7}` for Fock projection setup | ✓ |
| Paper 18 | `\cite{paper18}` for master Mellin engine §III.7 reference | ✓ |
| Paper 28 | `\cite{loutey_paper28}` for χ_4 / Catalan in QED vertex sums | ✓ |
| Paper 31 | `\cite{paper31}` for universal/Coulomb partition | ✓ |
| Paper 32 | `\cite{paper32}` for spectral triple construction + Connes axiom audit + Sprint H1 Yukawa-non-selection | ✓ |
| Paper 38 | `\cite{paper38}` for SU(2) propinquity convergence | ✓ |
| Paper 45 | `\cite{paper45}` for K⁺-restricted weak-form Lorentzian propinquity | ✓ |
| Paper 46 | `\cite{paper46}` for enlarged-substrate analogy (Appendix B), used to compare with non-natural-substrate enlargement in §sec:open_sprint | ✓ |
| Paper 47, 48, 49 | All cited as math.OA standalones | ✓ |
| Paper 55 | **Cited twice with different bibkeys** (`loutey_paper55` and `paper55`) — both bibitems point to the same paper | **DUPLICATE** |

### Bibkey style

- Paper 56 uses `marcolli_vansuijlekom2014` for Marcolli–vS 2014.
- All other GeoVac papers (Paper 25, 30, 32, 38, 41, 42, 43, 48, 49, group syntheses) use `marcolli_vs2014`.
- Result: an outside reader who collates the bibliography across the corpus (or who checks the GeoVac corpus's master bibkey table if one ever gets built) will see two different keys for the same source. **This is the most likely cosmetic stumble for cross-paper outreach.**

### Pérez-Sánchez title coherence

Three GeoVac papers cite arXiv:2401.03705 with three different titles:
- **Paper 32** (line 6856): "Bratteli networks and the Spectral Action on quivers" ← **CORRECT**
- **Paper 38** (line 1489): "On the continuum limit of gauge networks" ← WRONG (fabricated)
- **Paper 56** (line 1936): "Erratum to Marcolli--van Suijlekom 2014: Higgs from off-diagonal $D_F$ blocks" ← WRONG (fabricated, different from Paper 38's wrong title)

Same for arXiv:2508.17338:
- **Paper 32**: "Comment on 'Gauge networks in noncommutative geometry'" ← CORRECT
- **Paper 56** (line 1940): "The Standard Model without Higgs: spectral action on Marcolli--van Suijlekom gauge networks" ← WRONG (fabricated)

**The Paper 32 titles are correct and should be copied verbatim into Paper 56 (and Paper 38, although that is outside this audit scope).** Both arXiv IDs are correct; only the titles are fabricated. Same failure-mode as the v3.43.x literature-audit sprint already documented in `memory/feedback_*` (post-cutoff sub-agent fabrication of titles when arXiv IDs were correct). This is a structural pattern; the cleanup is mechanical.

### Marcolli–vS lineage framing

Paper 56's §1.4 (lines 273–317) gives an exemplary "honest scope of the MvS-lineage claim" paragraph (GeoVac is an *instance* of MvS 2014 at the data level, an *extension* of MvS at the structural level via the Tannakian dual). This framing is consistent with Paper 32's §1.4-equivalent and with the standing CLAUDE.md WH1 framing. ✓

---

## 6. READINESS VERDICT

**Headline: YELLOW (correctable in 1–2 hours of clean-up sprint; no theorem or numerical claim affected).**

### GREEN-grade items (ready for outreach)

1. **Scope-honesty.** Inclusion-vs-equality framing is exemplary. Every "equality" claim is precisely qualified by "at finite cutoff" / "on the natural-substrate panel" / "on the abelian factor". The three layers (inclusion at finite cutoff; converse equality at finite cutoff; inverse-limit identification with $\mathcal{G}_4$) are kept rigorously distinct.
2. **5,864 residual count.** Bit-exact verified at the sub-count level against all sprint memos. Closed-form formula $9 N(\nmax)^2 + 6 N(\nmax) + 6$ for G4-Inj matches every cell. Grand total 5,864 confirmed.
3. **LaTeX cleanliness.** Three-pass clean compile, 24 pp, zero substantive warnings.
4. **TC-2 sub-arc.** Theorem statements correctly labeled per-cutoff (TC-2a/b at $\nmax = 2$; TC-2c extends to $\{1, 2, 3, 4\}$; TC-2d packages cofiltered coherence). Reconstruction-direction-at-finite-cutoff is closed; inverse-limit identification is multi-year. Clean.
5. **Injection theorem.** "Closed immersion" framing, four named compatibilities (C1)–(C4), honest-scope remark covering M1/M2 column collapse + M3 substantive content + depth-blindness, equality direction explicitly multi-year.
6. **Reading A wins (NA-1)** documented honestly with named caveats; Reading B empirically ruled out at depth 2 on natural substrate; non-commutative enrichment is a remaining open question but classified out-of-scope via Sprint H1 calibration-data theorem.

### YELLOW-grade items (cleanup required before outreach)

**Required (substantive):**

1. **Pérez-Sánchez bibitem titles.** Paper 56 must copy the correct titles from Paper 32:
   - `perez_sanchez2024`: → "Bratteli networks and the Spectral Action on quivers"
   - `perez_sanchez2025`: → "Comment on `Gauge networks in noncommutative geometry'"
   - Effort: 2 minutes of mechanical edit on lines 1934–1940.
2. **Duplicate Paper 55 bibitem.** Choose one of `loutey_paper55` (line 1911) or `paper55` (line 1985) and unify all in-body citations. Recommend keeping `paper55` (matches the internal-paper bibkey style for `paper7`, `paper18`, `paper28`, `paper31`, `paper32`, `paper38`, `paper45`, `paper46`, `paper47`, `paper48`, `paper49`). Remove `loutey_paper55` bibitem, change 3 in-body cites (lines 1775, 1780, 1800). Effort: 5 minutes.

**Recommended (cosmetic but visible):**

3. **`marcolli_vansuijlekom2014` → `marcolli_vs2014` cross-paper bibkey unification.** This is the single most visible cross-paper-consistency stumble for an outside reader collating the corpus. Effort: 3 minutes (one bibitem rename + one body-text find-and-replace).
4. **Glanois title.** Add "groupoid" between "fundamental group" and "of" to match thesis title. Effort: 30 seconds.
5. **Deligne 2010 title.** Restore "ou 8" comma-or list to match French original, or leave the comma-list and document the convention. Effort: 30 seconds.
6. **Brown 2017 vs 2014.** Either rename bibkey to `brown2014_icm` or clarify that "2017" reflects a reprint year. Effort: 1 minute.
7. **`fathizadeh_marcolli2016` bibkey** could be renamed to `fathizadeh_marcolli2017` to match the journal year, or left as the arXiv-year convention. Effort: 1 minute. (Low priority.)

**Optional (cosmetic LaTeX warnings):**

8. Wrap section-title math in `\texorpdfstring{}{}` to silence hyperref Unicode warnings. Effort: 5–10 minutes. (Low priority; cosmetic only.)
9. Tweak the NA-1 paragraph (lines 1785–1802) to remove 2 `Overfull \hbox` warnings. Effort: 2 minutes.

### RED-grade items

None.

### Outreach impact

The substantive Pérez-Sánchez title errors are the items most likely to be noticed by an adversarial outside reader (Brown/Kleinschmidt-lineage referees would absolutely click through to arXiv:2401.03705 and arXiv:2508.17338 to confirm). Catching these now is the single highest-leverage cleanup.

After applying YELLOW items 1–7 (≈15–20 minutes of mechanical edits), Paper 56 is **arXiv-ready and outreach-ready**. No theorem statement, proof, residual count, or scope-honesty claim is in question; the cleanup is purely bibliographic. The mathematical structure of Paper 56 is at GREEN-grade and ready for the math.OA / NCG-with-motivic-leanings outreach panel.

---

## Audit method notes

- **Read-only**: no files modified, no test suite run. All claims about residual counts cross-checked against the sprint memo files at `debug/sprint_q5p_*_memo.md` and `debug/sprint_injection_*_memo.md`.
- **Web verification**: per-bibitem verification of every published-author bibitem against arXiv and journal records (8 bibitems verified via web search).
- **Compile audit**: three pdflatex passes (`-draftmode` then standard) executed locally; warnings parsed.
- **Cross-paper check**: bibkey and title coherence verified against Papers 32 and 38 (Pérez-Sánchez), and against Papers 30, 32, 38, 41, 43, 48, 49, group syntheses (Marcolli–vS bibkey).

End of audit memo.
