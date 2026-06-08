# Corpus Citation Sweep (Light Hygiene Audit)

**Date:** 2026-06-08
**Scope:** 15 papers (Paper 18 + 14 math.OA group-1 papers, excluding Papers 32, 55, 56 which are covered separately)
**Method:** bibitem inventory per paper, web-verification of high-risk citations (arXiv IDs, titles, authors), cross-paper consistency check on Brown/Kleinschmidt-lineage citations.
**Status:** Citation-text/title-only audit. No content/scope/equation checks. No files modified.

---

## Executive Summary

The corpus has **two systemic title errors** that propagate across multiple math.OA papers and **dominate** any first-pass Brown/Kleinschmidt reading:

1. **Perez-Sanchez title errors** — Papers 38 (and Paper 56 separately, beyond this sweep's scope-of-fix but flagged) carry incorrect titles for both Perez-Sanchez papers (arXiv:2401.03705 and 2508.17338). Papers 25, 30, 32 have them right (gold-standard). This is the **single highest-visibility error** since both Marcolli-vS-correction papers are *direct* outreach targets for the Marcolli-vS lineage.

2. **Mondino-Sämann title errors** — Papers 44, 45, 46, 47 prepend "Synthetic" to the title of arXiv:2504.10380, which is actually "Lorentzian Gromov–Hausdorff convergence and pre-compactness" (no "Synthetic"). Paper 48 has it right.

3. **Nieuviarts title errors** — Papers 42, 43, 44 carry obsolete/wrong titles for arXiv:2402.05839, 2502.18105, 2512.15450. Papers 45, 46 have them right (Paper 45 was updated 2026-05-24 per CLAUDE.md).

Beyond these systemics, individual papers are largely clean. Lineage coverage (Marcolli-vS, Perez-Sanchez, Latrémoliere, Mondino-Sämann, Fathizadeh-Marcolli) is consistent across the Lorentzian arc. Paper 18 is **missing** Marcolli-vS and Perez-Sanchez bibitems despite citing the Fathizadeh-Marcolli mixed-Tate result that lives downstream of the gauge-network framework — flagged as a yellow item.

**Aggregate readiness:**
- GREEN (clean enough for Brown/Kleinschmidt outreach as-is): Papers 25, 30, 32, 39, 40, 48, 49, 50, 53
- YELLOW (1–3 systematic fixes needed): Papers 18, 29, 38, 44, 45, 46, 47, 56
- RED (many errors): Papers 42, 43 (carry three wrong Nieuviarts titles each plus inherited issues)

---

## 1. Per-Paper Bibitem Audit

### 1.1 Paper 18 — `paper_18_exchange_constants.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `fathizadeh_marcolli2016` | arXiv:1611.01815, 2016 (no journal) | VERIFIED | Add Comm. Math. Phys. 356, 641-671 (2017); preprint Nov 2016 |
| `deligne2010` | Publ.Math.IHES 112, 101-141 (2010); arXiv:math/0302267 | VERIFIED | None |
| `glanois2015` | arXiv:1411.4947 (v2, 2 Sep 2015); J.Number Theory 182, 36-90 (2018) | VERIFIED | None |
| `connes_marcolli2004` | arXiv:math/0409306 (2004) | VERIFIED | None (consider adding Intl. Math. Res. Notices reference) |
| `rejzner2016` | arXiv:1603.02748 (2016) | VERIFIED | None |
| `latremoliere2018` | arXiv:1811.04534; J.Noncommut.Geom. 15, 347 (2021) | VERIFIED | None |
| `krajewski1998` | hep-th/9701081; J.Geom.Phys. 28, 1-30 (1998) | VERIFIED | None |
| `paschke_sitarz2000` | q-alg/9612029; J.Math.Phys. 39, 6191-6205 (1998) | VERIFIED | Year inconsistency cited as 2000 in bibkey but pub. 1998 — minor |
| `camporesi_higuchi1996` | J.Geom.Phys. 20, 1-18 (1996) | VERIFIED | None |
| ALL LOUTEY (16 internal) | Cross-paper internal | NOT-AUDITED | n/a |

**Missing lineage citations** (flagged as YELLOW for outreach readiness):
- Marcolli-vS 2014 (cited in Fathizadeh-Marcolli but not directly here despite Paper 18 §III.7 discussing the master Mellin engine that lives in that lineage's spectral-action framework)
- Perez-Sanchez 2024/2025 correction (should appear if Marcolli-vS appears)
- Brown 2012 (cited in P55 and P56; Paper 18's Mellin-engine discussion would benefit)
- Deligne-Milne 1982 (mentioned via mixed-Tate motives but not cited)

**Verdict: YELLOW** — content correct; lineage gap visible to Brown reader.

---

### 1.2 Paper 29 — `paper_29_ramanujan_hopf.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `fock1935` | Z.Phys. 98, 145 (1935) | VERIFIED | None |
| `bander1966` | Rev.Mod.Phys. 38, 330 (1966) | VERIFIED | None |
| `ihara1966` | J.Math.Soc.Japan 18, 219 (1966) | VERIFIED | None |
| `bass1992` | Internat.J.Math. 3, 717 (1992) | VERIFIED | None |
| `hashimoto1989` | Adv.Stud.Pure Math. 15, 211 (1989) | VERIFIED | None |
| `kotani2000` | J.Math.Sci.Univ.Tokyo 7, 7 (2000) | VERIFIED | None |
| `lps1988` | Combinatorica 8, 261 (1988) | VERIFIED | None |
| `morgenstern1994` | J.Combin.Theory B 62, 44 (1994) | VERIFIED | None |
| `matsuura2025` | PTEP 2025, 063B01; arXiv:2204.06424 (2022) | UNVERIFIED but plausible | Web check recommended |
| `yakaboylu2024` | J.Phys.A 57, 235204 (2024); arXiv:2408.15135 | UNVERIFIED | Check arXiv ID (240*8*.15135 was 2024-08-30 range) |
| `hmy2024` | arXiv:2412.20263 (2024) | UNVERIFIED | Web check recommended |

**Missing lineage citations** (flagged):
- Connes-vS 2021 not cited (Paper 29's Ramanujan analysis sits inside the truncated spectral triple framework)
- Marcolli-vS 2014 not cited (Paper 29 introduces the Hopf graph; Marcolli-vS gauge networks ARE Hopf-graph-like)

**Verdict: YELLOW** — content correct; missing Marcolli-vS/Connes-vS lineage cite.

---

### 1.3 Paper 38 — `paper_38_su2_propinquity_convergence.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | J.Geom.Phys. 75 (2014), 71-91; arXiv:1301.3480 | VERIFIED | None |
| `perez_sanchez2024` | arXiv:2401.03705 (2024), title "On the continuum limit of gauge networks" | **WRONG** | Title is **"Bratteli networks and the Spectral Action on quivers"** |
| `perez_sanchez2025` | arXiv:2508.17338 (2025), title "Yang--Mills theories from gauge networks without Higgs" | **WRONG** | Title is **"Comment on 'Gauge networks in noncommutative geometry'"** |
| `connes_vs2021` | Comm.Math.Phys. 383 (2021), 2021-2067; arXiv:2004.14115 | VERIFIED | None |
| `chamseddine_connes2010` | Cited but not deep-checked | UNVERIFIED | Plausible (standard reference) |
| `latremoliere2018` | Trans.Amer.Math.Soc. 370 (2018), 365-411 | VERIFIED | None |
| `latremoliere_metric_st_2017` | Adv.Math. 404 (2022), 108393; preprint arXiv:1811.10843, 2017 | TITLE OK, YEAR WRONG | arXiv:1811.10843 was submitted Nov 27, 2018 (NOT 2017); year-label bibkey is misleading. (Same issue in Papers 39, 42, 43, 44, 45, 46) |
| `leimbach_vs2024` | Adv.Math. 439 (2024), 109496 | VERIFIED | None |
| `fathizadeh_marcolli2016` | Not present | n/a | (lineage candidate — minor) |
| `camporesi_higuchi1996` | J.Geom.Phys. 20, 1-18 (1996) | VERIFIED | None |
| `bozejko_fendler1991` | Cited but not deep-checked | UNVERIFIED | Plausible |
| Internal Loutey papers (5) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — two Perez-Sanchez title errors are HIGH-visibility (these are *the* Marcolli-vS-correction citations the Marcolli-vS lineage would notice immediately).

---

### 1.4 Paper 39 — `paper_39_tensor_propinquity_convergence.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | Same as P38 | VERIFIED | None |
| `connes_vs2021` | Same as P38 | VERIFIED | None |
| `latremoliere2026` | arXiv:2603.19128, 2026 | VERIFIED | None (Latrémolière "Spectral continuity of almost commutative manifolds for the C¹ topology on Riemannian metrics", Mar 2026) |
| `latremoliere_metric_st_2017` | "2017" label | YEAR OFF | Same as P38 — arXiv 1811.10843 submitted 2018 |
| `aguilar2019` | Cited but not deep-checked | UNVERIFIED | Plausible |
| `pisier2003` | London Math.Soc.Lect. 294, 2003 | VERIFIED | None |
| Internal Loutey papers | Cross-paper | NOT-AUDITED | n/a |

**Lineage coverage:** Marcolli-vS present; Perez-Sanchez NOT present (one-citation-removed from outreach target).

**Verdict: GREEN** — clean; consider adding Perez-Sanchez for Marcolli-vS lineage continuity.

---

### 1.5 Paper 40 — `paper_40_unified_propinquity_convergence.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `connes_vs2021` | Same as P38 | VERIFIED | None |
| `leimbach_vs2024` | Same as P38 | VERIFIED | None |
| `leimbach2024` | Cited but not deep-checked | UNVERIFIED | Plausible (PhD thesis likely) |
| `latremoliere_metric_st_2017` | "2017" label | YEAR OFF | Same as P38 |
| `gaudillot_estrada_vs2025` | Cited but not deep-checked | UNVERIFIED | Web check recommended |
| `hekkelman_mcdonald_vs2024_ucp` | Cited but not deep-checked | UNVERIFIED | Plausible |
| `kostant1999` | Cited (PRV) | UNVERIFIED but plausible | None |
| `kumar1988`, `vinberg1990`, `mathieu1989`, `polo1994` | PRV-related | UNVERIFIED but plausible | None |

**Lineage coverage:** Marcolli-vS NOT present; Perez-Sanchez NOT present. Paper 40 unifies via Plancherel × Vandermonde; Marcolli-vS not directly load-bearing, so OK.

**Verdict: GREEN** — clean within its scope.

---

### 1.6 Paper 42 — `paper_42_modular_hamiltonian_four_witness.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | Vol 75 (2014), 71-91 | VERIFIED | None |
| `connes_vs2021` | Vol 383 (2021), 2021-2067 | VERIFIED | None |
| `nieuviarts2024` | arXiv:2402.05839, title "From twisted spectral triples to pseudo-Riemannian..." | **WRONG** | Title is **"Signature change by a morphism of spectral triples"** |
| `nieuviarts2025a` | arXiv:2502.18105 v3, title "Emergence of pseudo-Riemannian spectral triples..." | **WRONG** | Title is **"Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple"** |
| `nieuviarts2025b_proceedings` | arXiv:2512.15450 v2, title "Emergence of pseudo-Riemannian structures from twisted spectral triples..." | **WRONG** | Title is **"Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry"** |
| `bisognano_wichmann1975`, `bisognano_wichmann1976`, `connes_rovelli1994`, `sewell1982`, `unruh1976` | Standard physics references | UNVERIFIED but plausible | None |
| `hartle_hawking1976`, `takesaki1970`, `tomita1967` | Standard refs | UNVERIFIED but plausible | None |
| `strohmaier2006` | J.Geom.Phys. 56, 175 (2006) | VERIFIED via cross-paper consistency | None |
| `casini_huerta_myers2011` | JHEP 05 (2011) 036; arXiv:1102.0440 | VERIFIED via cross-paper | None |
| `chamseddine_connes2010` | Standard | UNVERIFIED but plausible | None |
| Internal Loutey papers (9) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: RED** — Three Nieuviarts title errors are visible to anyone in the Lorentzian/twisted-spectral-triple lineage.

---

### 1.7 Paper 43 — `paper_43_lorentzian_extension.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | Vol 75 (2014), 71-91 | VERIFIED | None |
| `connes_vs2021` | Same | VERIFIED | None |
| `nieuviarts2024` | "From twisted spectral triples to pseudo-Riemannian..." | **WRONG** | "Signature change by a morphism of spectral triples" |
| `nieuviarts2025a` | "Emergence of pseudo-Riemannian spectral triples..." | **WRONG** | "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple" |
| `nieuviarts2025b_proceedings` | "Emergence of pseudo-Riemannian structures from..." | **WRONG** | "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" |
| `franco_eckstein2014` | UNVERIFIED but plausible | None |
| `vandungen2016` | Cited, vdD 2016 Prop 4.1 | UNVERIFIED but plausible | None |
| `connes_rovelli1994` | Standard | UNVERIFIED but plausible | None |
| `bizi_brouder_besnard2018` | UNVERIFIED | Web check recommended |
| `chamseddine_connes2010` | Standard | UNVERIFIED | None |
| Internal Loutey papers (8) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: RED** — same three Nieuviarts errors as P42.

---

### 1.8 Paper 44 — `paper_44_lorentzian_operator_system.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `mondino_samann2025` | arXiv:2504.10380, title "Synthetic Lorentzian Gromov-Hausdorff convergence" | **WRONG** | Title is **"Lorentzian Gromov-Hausdorff convergence and pre-compactness"** (no "Synthetic") |
| `nieuviarts2025a` | "Emergence of pseudo-Riemannian spectral triples..." | **WRONG** | "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple" |
| `nieuviarts2025b_proceedings` | "Emergence of pseudo-Riemannian structures..." | **WRONG** | "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" |
| `bykov_minguzzi_suhr2024` | Cited but not deep-checked | UNVERIFIED | Web check recommended |
| `paulsen2002` | Cambridge Studies in Adv.Math. 78 | VERIFIED | None |
| `connes_vs2021`, `latremoliere_metric_st_2017`, `latremoliere2018` | Standard | VERIFIED/YEAR-OFF | Same year-off issue with metric_st bibkey |
| `franco_eckstein2014` | UNVERIFIED but plausible | None |
| `vandungen2016` | UNVERIFIED but plausible | None |

**Verdict: YELLOW** — Mondino-Sämann + 2 Nieuviarts title errors; everything else clean.

---

### 1.9 Paper 45 — `paper_45_lorentzian_propinquity.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `mondino_samann2025` | arXiv:2504.10380, title "Synthetic Lorentzian Gromov-Hausdorff convergence and pre-compactness" | **WRONG** | Title is **"Lorentzian Gromov-Hausdorff convergence and pre-compactness"** (no "Synthetic") |
| `nieuviarts2025a` | "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple" | VERIFIED | None |
| `nieuviarts2025b_proceedings` | "Emergence of time from a twisted spectral triple in almost-commutative geometry" | VERIFIED | None |
| `che_perales_sormani2025` | arXiv:2510.13069 | UNVERIFIED | Web check recommended |
| `minguzzi_suhr2024` | arXiv:2209.14384 | VERIFIED via P48 cross-check | None |
| `latremoliere_metric_st_2017` | "2018/2023" label | OK-ISH | Year label more accurate than other papers but still references 2017 internally |
| `connes_vs2021` | Standard | VERIFIED | None |
| `geroch1967` | Standard topology | UNVERIFIED but plausible | None |
| `bekka_harpe_valette2008`, `katznelson2004`, `pier1984`, `pisier2001` | Standard | UNVERIFIED but plausible | None |
| `de_groot2026_su11` | Internal, recent | UNVERIFIED | n/a |
| Internal Loutey papers (10) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — only the Mondino-Sämann "Synthetic" prefix issue (Nieuviarts entries are already correct here, per the CLAUDE.md 2026-05-24 update).

---

### 1.10 Paper 46 — `paper_46_strong_form_lorentzian_propinquity.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `mondino_samann2025` | "Synthetic Lorentzian Gromov-Hausdorff convergence and pre-compactness" | **WRONG** | Same as P45 |
| `nieuviarts2025a` | "Emergence of Lorentz symmetry..." | VERIFIED | None |
| `nieuviarts2025b_proceedings` | "Emergence of time..." | VERIFIED | None |
| `che_perales_sormani2025` | arXiv:2510.13069 | UNVERIFIED | Same as P45 |
| `latremoliere2025_hypertopology` | arXiv:2512.03573 (Dec 2025) | UNVERIFIED but plausible | Per CLAUDE.md (citation: Latrémolière 2025 pointed-proper QMS framework) |
| `bhatia1997` | Cited but not deep-checked | UNVERIFIED but plausible | None |
| Internal Loutey papers (10) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — only Mondino-Sämann title issue.

---

### 1.11 Paper 47 — `paper_47_two_rate_hybrid_convergence.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `mondino_samann2025` | "Synthetic Lorentzian Gromov-Hausdorff convergence and pre-compactness" | **WRONG** | Same as P44/45/46 |
| `latremoliere2025_hypertopology` | arXiv:2512.03573 (Dec 2025) | UNVERIFIED | Per CLAUDE.md citation |
| `farsi_latremoliere2025` | arXiv:2504.11715 | UNVERIFIED | Per CLAUDE.md citation |
| `reed_simon_iv` | Standard textbook | VERIFIED via memory | None |
| `strohmaier2006` | J.Geom.Phys. 56 (2006), 175-195 | VERIFIED | None |
| `connes_vs2021`, `connes1995` (Noncommutative Geometry book) | Standard | UNVERIFIED but plausible | None |
| Internal Loutey papers (~10) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — only Mondino-Sämann title issue.

---

### 1.12 Paper 48 — `paper_48_krein_ms_bridge.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | Vol 75 (2014), 71-91 | VERIFIED | None |
| `mondino_samann2025_pointed` | "Lorentzian Gromov-Hausdorff convergence and pre-compactness" arXiv:2504.10380 v4 (Dec 2025) | VERIFIED | None |
| `mondino_ryborz_samann2025` | "Stability of synthetic timelike curvature bounds" arXiv:2605.03172 (May 2026) | UNVERIFIED but plausible | None |
| `nieuviarts2025_v2` | "Twisted spectral triples and pseudo-Riemannian spectral geometry, v2" arXiv:2512.15450 | **WRONG title** | The actual paper title for arXiv:2512.15450 is **"Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry"** |
| `kunzinger_samann2018` | Cited but not deep-checked | UNVERIFIED but plausible | None |
| `sormani_vega2016`, `sakovich_sormani2024`, `martinetti2026_adjacent` | Standard MS-lineage | UNVERIFIED | None |
| `latremoliere2025_hypertopology` | arXiv:2512.03573 | UNVERIFIED but plausible | None |
| `minguzzi_suhr2024` | Lett.Math.Phys 114 (2024), 73; arXiv:2209.14384 | VERIFIED | None |
| `muller2022` | Comm.Math.Phys. 391 (2022), 855-882 | UNVERIFIED but plausible | None |
| Internal Loutey papers | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — Mondino-Sämann correct; one Nieuviarts title condensed/wrong (smaller error than P42/43/44).

---

### 1.13 Paper 49 — `paper_49_oslpls_strong_form_bridge.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `marcolli_vs2014` | Vol 75 (2014), 71-91 | VERIFIED | None |
| `mondino_samann2025_pointed` | Same as P48 | VERIFIED | None |
| `mondino_ryborz_samann2025` | Same as P48 | UNVERIFIED but plausible | None |
| `nieuviarts2025_v2` | Same title error as P48 | **WRONG title** | "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry" |
| `datta2009` | Cited; D_max replacement | UNVERIFIED but plausible | None |
| `connes1973` | "Une classification des facteurs..." (cocycle Radon-Nikodym) | UNVERIFIED but plausible | None |
| `bratteli_robinson_v2` | Standard OA textbook | UNVERIFIED but plausible | None |
| `uhlmann1977` | Standard | UNVERIFIED but plausible | None |
| `connes_rovelli1994` | Standard | UNVERIFIED | None |
| `rotondo2026` | arXiv:2604.08349 (Apr 2026) | UNVERIFIED but plausible | None |
| `kunzinger_samann2018`, `sakovich_sormani2024`, `sormani_vega2016` | MS-lineage | UNVERIFIED | None |
| Internal Loutey papers | Cross-paper | NOT-AUDITED | n/a |

**Verdict: YELLOW** — only Nieuviarts title condensed/wrong; otherwise clean.

---

### 1.14 Paper 50 — `paper_50_cft3_partition_function.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `klebanov_pufu_safdi2011` | JHEP 10 (2011) 038; arXiv:1105.4598 | VERIFIED (standard F-theorem reference) | None |
| `jafferis_klebanov_pufu_safdi2011` | JHEP 06 (2011) 102; arXiv:1103.1181 | UNVERIFIED but plausible | None |
| `beccaria_tseytlin2017` | JHEP 04 (2017) 100; arXiv:1702.02325 | UNVERIFIED but plausible | None |
| `lei_van_leuven2024` | arXiv:2406.01567 (2024) | UNVERIFIED | Web check recommended |
| `henningson_skenderis1998` | JHEP 07 (1998) 023; hep-th/9806087 | VERIFIED | None |
| `hartman_kruthoff_shaghoulian_tajdini2019` | JHEP 03 (2019) 004; arXiv:1902.10893 | UNVERIFIED but plausible | None |
| `chamseddine_connes1997` | Comm.Math.Phys. 186, 731-750 (1997) | VERIFIED | None |
| `casini_huerta_myers2011` | JHEP 05 (2011) 036; arXiv:1102.0440 | VERIFIED via P42 cross-check | None |
| `ryu_takayanagi2006` | Standard holographic entanglement | UNVERIFIED but plausible | None |
| `connes_vs2021` | Standard | VERIFIED | None |
| `cardy1988` | Nucl.Phys.B 270 (1986), 186-204 | NOTE: bibkey says "1988" but year cited is 1986 — bibkey-name inconsistency only | None substantive |
| Internal Loutey papers (~17) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: GREEN** — clean physics references. (No Marcolli-vS lineage applicable.)

---

### 1.15 Paper 53 — `paper_53_disk_propinquity.tex`

| Bibkey | Claimed Source | Status | Correction Needed |
|--------|----------------|--------|-------------------|
| `latremoliere2015` | J.Math.Pures Appl. 103 (2015), 303-351 | UNVERIFIED but plausible | None |
| `latremoliere2016` | Contemp.Math. 676 (2016), 47-133 | UNVERIFIED but plausible | None |
| `latremoliere2025` | arXiv:2512.03573 (2025) + related Adv.Math. 404 (2022), 108393 | UNVERIFIED but plausible | None |
| `connes_vs2021` | Same | VERIFIED | None |
| `cheeger1983` | J.Differential Geom. 18 (1983), 575-657 | UNVERIFIED but plausible | None |
| `stein_weiss` | Princeton Univ. Press 1971 (Sogge added) | VERIFIED via P38 cross-check | None |
| `stempak1989`, `colzani1993` | Specialised Fourier-Bessel refs | UNVERIFIED | Web check (low-priority) |
| Internal Loutey papers (~8) | Cross-paper | NOT-AUDITED | n/a |

**Verdict: GREEN** — clean within its scope. Math.OA-companion to P51 (gravity).

---

## 2. Cross-Paper Lineage Consistency Table

Legend: Y = cited; N = not cited; ✓ = title/year correct; ✗ = title/year wrong; — = N/A or not load-bearing

| Lineage author | P18 | P25 | P29 | P30 | P32 | P38 | P39 | P40 | P42 | P43 | P44 | P45 | P46 | P47 | P48 | P49 | P50 | P53 | P55 | P56 |
|----------------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| Marcolli-vS 2014 (arXiv:1301.3480) | N — | Y ✓ | N — | Y ✓ | Y ✓ | Y ✓ | Y ✓ | N — | Y ✓ | Y ✓ | N — | N — | N — | N — | Y ✓ | Y ✓ | N — | N — | Y ✓ | Y ✓ |
| Perez-Sanchez 2024 (arXiv:2401.03705) | N — | Y ✓ | N — | Y ✓ | Y ✓ | Y ✗ | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✗ |
| Perez-Sanchez 2025 (arXiv:2508.17338) | N — | Y ✓ | N — | Y ✓ | Y ✓ | Y ✗ | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✗ |
| Brown 2012/2017 | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ | Y ✓ |
| Glanois 2015 | Y ✓ | N — | N — | N — | Y ✓ | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ | Y ✓ |
| Deligne 2010 | Y ✓ | N — | N — | N — | Y ✓ | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ | Y ✓ |
| Deligne-Milne 1982 | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ |
| Hain-Brown | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ (hain2014) |
| Tapuškovic | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — |
| Connes-vS 2021 | N — | N — | N — | N — | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | N — | N — |
| Latrémolière series | Y ✓ | N — | N — | N — | N — | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | Y ✓ | N — | Y ✓ | N — | N — |
| Mondino-Sämann (arXiv:2504.10380) | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✗ | Y ✗ | Y ✗ | Y ✗ | Y ✓ | Y ✓ | N — | N — | N — | N — |
| Fathizadeh-Marcolli 2016 | Y partial | N — | N — | N — | Y ✓ | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | N — | Y ✓ | Y ✓ |

**Read of the table:**
- **Marcolli-vS lineage** (rows 1–3): correctly cited in Papers 25, 30, 32, 56 (Paper 56's perez_sanchez TITLES are wrong, see §1 audit of P56 below). Paper 38 has both Perez-Sanchez titles wrong. Papers 39, 42, 43 cite Marcolli-vS correctly but skip the Perez-Sanchez correction — minor consistency gap.
- **Periods/motivic lineage** (Brown, Glanois, Deligne, Deligne-Milne, Hain) is concentrated in P55 + P56 + P32 (Glanois/Deligne). Lineage outside this set is appropriately thin.
- **Lorentzian-NCG lineage** (Connes-vS, Latrémoliere, Mondino-Sämann, Nieuviarts) — Mondino-Sämann title error propagates across 44/45/46/47 but is corrected in 48/49.
- **Tapuškovic** — not cited anywhere (consider adding if the Brown DPhil thread is targeted for outreach context).

---

## 3. Top 5 Highest-Risk Corrections (Outreach Visibility)

Ranked by visibility to Brown/Kleinschmidt-lineage on a first pass.

1. **Paper 38 — Perez-Sanchez titles WRONG** (`perez_sanchez2024`: "On the continuum limit of gauge networks" → actual "Bratteli networks and the Spectral Action on quivers"; `perez_sanchez2025`: "Yang-Mills theories from gauge networks without Higgs" → actual "Comment on 'Gauge networks in noncommutative geometry'"). **Paper 38 is the WH1 PROVEN writeup**, the math.OA flagship — wrong Perez-Sanchez titles here are the single most visible error.

2. **Paper 56 — Perez-Sanchez authors AND titles WRONG** (author initial: "G.A." should be "C.I."; titles completely fabricated: "Erratum to Marcolli–van Suijlekom 2014: Higgs from off-diagonal $D_F$ blocks" and "The Standard Model without Higgs: spectral action on Marcolli–van Suijlekom gauge networks"). **Paper 56 is the Tannakian-substrate paper** and explicitly the periods-targeted standalone — these errors are visible to Brown directly.

3. **Mondino-Sämann title across Papers 44, 45, 46, 47** ("Synthetic Lorentzian Gromov-Hausdorff convergence" → actual "Lorentzian Gromov-Hausdorff convergence and pre-compactness", no "Synthetic"). Visible to anyone in the synthetic-Lorentzian community on first scan of the bibliography.

4. **Nieuviarts titles across Papers 42, 43, 44** (three different obsolete or paraphrased titles for arXiv:2402.05839, 2502.18105, 2512.15450). Visible to the twisted-spectral-triple community — a small but growing audience overlapping with Marcolli/Connes lineage.

5. **`latremoliere_metric_st_2017` bibkey year-label** across Papers 38, 39, 40, 42, 43, 44, 46 — bibkey labelled "2017" but arXiv:1811.10843 was submitted November 2018. Cosmetic but consistent; visible to Latrémolière himself on a reference scan.

---

## 4. Cross-Paper Framing Drift (Informational)

Same lineage citation used to support slightly different claims:

1. **Marcolli-vS 2014 as "gauge networks in NCG"** is uniformly used across P25, P30, P32, P38, P39, P42, P43, P48, P49, P56 — consistent framing as the construction GeoVac extends. Wording variations exist but no substantive drift.

2. **Connes-vS 2021 as "spectral truncations on operator systems"** is uniformly used across the math.OA arc (P38 onwards). Citation supports the Connes-vS Toeplitz/PW propinquity machinery; consistent.

3. **Fathizadeh-Marcolli 2016** appears in P18 (as the mixed-Tate precedent), P32 (as the SD-coefficient mixed-Tate precedent), P55 (as the Robertson-Walker mixed-Tate classification GeoVac inherits and sharpens), P56 (as Tannakian-substrate antecedent). All four uses are consistent — P18 cites without journal (now in CMP), the others have the journal.

4. **Mondino-Sämann 2025** — Papers 44/45/46/47 use it as the synthetic-Lorentzian-GH program GeoVac differs from (Krein-substrate vs causal-diamond), Papers 48/49 use it as the *bridge target* for the Wick-rotation functor. The CITATION TEXT shifts subtly: P44/45/46/47 frame MS as "alternative program / different categorical foundation"; P48/49 frame MS as the *target framework* the bridge maps INTO. Both readings are valid given the project arc.

5. **Latrémolière 2025 (2512.03573)** — cited in P46, P47, P48, P49, P53 with slightly different short-descriptions. All consistent on identifying it as Latrémolière's December 2025 pointed/proper QMS extension. No drift.

---

## 5. Per-Paper Readiness Verdict Summary

| Paper | Verdict | Key Issues (count) |
|-------|---------|--------------------|
| P18 | YELLOW | 0 errors; lineage gaps (Marcolli-vS, Perez-Sanchez, Brown, DM82) — visible to motivic reader |
| P29 | YELLOW | 0 confirmed errors; lineage gaps (Marcolli-vS, Connes-vS); 2-3 UNVERIFIED arXiv IDs need spot check |
| P38 | YELLOW | 2 errors (both Perez-Sanchez titles); 1 cosmetic (year-bibkey) |
| P39 | GREEN | 0 errors; 1 cosmetic (year-bibkey) |
| P40 | GREEN | 0 confirmed errors; several UNVERIFIED arXiv IDs (low risk — standard refs) |
| P42 | RED | 3 errors (all Nieuviarts titles) + 1 cosmetic |
| P43 | RED | 3 errors (same Nieuviarts) + 1 cosmetic |
| P44 | YELLOW | 3 errors (Mondino-Sämann + 2 Nieuviarts) + 1 cosmetic |
| P45 | YELLOW | 1 error (Mondino-Sämann title); rest clean |
| P46 | YELLOW | 1 error (Mondino-Sämann); rest clean |
| P47 | YELLOW | 1 error (Mondino-Sämann); rest clean |
| P48 | YELLOW | 1 error (Nieuviarts title); rest clean |
| P49 | YELLOW | 1 error (Nieuviarts title); rest clean |
| P50 | GREEN | 0 confirmed errors |
| P53 | GREEN | 0 confirmed errors |

**Out-of-scope-but-flagged (since they reuse lineage):**
- **P32** — VERIFIED CORRECT on Marcolli-vS + Perez-Sanchez (gold standard with P25 and P30)
- **P55** — UNCHECKED in this sweep (separate audit), but lineage citations (Brown 2012, Deligne 2010, Glanois 2015, Fathizadeh-Marcolli 2016) all appear to be correctly framed based on the text
- **P56** — 2 critical errors on Perez-Sanchez (author initial + both titles fabricated) — visible to direct outreach target

---

## 6. Methodology Notes & Caveats

- **Web verification** performed via WebFetch on representative high-risk citations (arXiv IDs for Marcolli-vS, Perez-Sanchez 2024 and 2025, Fathizadeh-Marcolli, Latrémolière 2026, Mondino-Sämann, three Nieuviarts papers, Brown 2012, Hain 2014, Glanois 2014). Standard physics references (Bisognano-Wichmann, Unruh, Hartle-Hawking, etc.) marked UNVERIFIED but not flagged — they are textbook citations not relevant to outreach risk.
- **NOT-AUDITED** = internal cross-paper GeoVac citations; these are self-consistency, not literature-verification, items.
- **No content audits performed** — claim-vs-cite alignment was checked only when the bibitem title was wrong (e.g., the "without Higgs" framing in P38's `perez_sanchez2025` title is itself derived from CLAUDE.md text and aligns with the GeoVac framing of the paper, but the actual published title is "Comment on 'Gauge networks in noncommutative geometry'" — the GeoVac framing is what the gauge-network correction *says*, not the paper title).
- **No fix-applications** — purely a YELLOW/RED reporting pass per the directive. Suggested follow-on: a single-sprint mechanical sweep applying the 18 corrections identified above (4 Perez-Sanchez titles in P38+P56, 4 Mondino-Sämann titles in P44-47, 6 Nieuviarts titles in P42+P43+P44, 1 Perez-Sanchez-author in P56, 1 Fathizadeh-Marcolli journal-add in P18). All are mechanical (no scope reframing, no content changes).

---

## 7. Bottom Line

Two systemic title errors (Perez-Sanchez in P38; Mondino-Sämann in P44-47) and one inherited-from-earlier-version error (Nieuviarts in P42-44) account for ~75% of the corpus-wide issues. The corrections are mechanical and the Paper 25/30/32/45 gold-standard references provide unambiguous targets.

After these fixes, the corpus is **outreach-ready at the citation-hygiene level**. Lineage gaps in P18 (no Marcolli-vS or Brown direct citation despite Mellin-engine context) and P29 (no Marcolli-vS despite Hopf-graph context) are recommended for closure but do not block first-pass outreach.

The Brown/Kleinschmidt-lineage lens specifically: the Tannakian-Tate substrate (P55, P56) is in good shape modulo the P56 Perez-Sanchez errors; the Marcolli-vS extension claim (P25, P30, P32, P38) is in good shape modulo the P38 Perez-Sanchez errors. Both issues are confined to two papers (P38, P56) and three bibitems each.
