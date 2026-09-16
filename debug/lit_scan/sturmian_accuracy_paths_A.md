# Lit scan A — the molecular Coulomb-Sturmian *accuracy* record

**Date:** 2026-09-13
**Axis:** A (assigned). Has anyone reached chemical accuracy (<= ~1.6 mHa) on a
MOLECULE (two or more centres) with a Sturmian basis, and what technical move
did it?
**Out of scope (other agents):** explicit correlation; STO/B-function integral
technology; two-centre *spheroidal* Sturmians.
**Budget:** 25/25 web calls spent. 9 sources cited, 5 VERIFIED.

---

## Bottom line

**No.** Across the whole one-centre-per-nucleus Sturmian lineage (Avery & Avery;
Herbst–Avery–Dreuw / molsturm; Aquilanti–Coletti) there is **no published
molecular calculation at chemical accuracy, and no published molecular error
budget at all.** The literature partitions cleanly into (a) *method proposals*,
(b) *pilot calculations* (the authors' own word), and (c) *integral technology*.
The most recent primary source in the lineage (Herbst–Avery–Dreuw 2019, full
text read) is **atoms only, Hartree–Fock only**, and lists "Molecular Coulomb
Sturmians?" as an open question.

The one genuinely useful find is a **named technical move the corpus has not
recorded**: Avery's molecular ERI construction resolves the *product* of two
Coulomb Sturmians onto an auxiliary CS set at the **doubled exponent 2k**, by
solving one linear system against the Shibuya–Wulfman matrix S(2k, S). Herbst's
own slide labels it "**Like exact density fitting**". It is integral technology,
not accuracy — but it is the only place in the canon where a *second scale*
appears, and the accuracy bottleneck in both canons is exactly the absence of a
second scale.

---

## Table

| Path | Best demonstrated accuracy + system + reference method | Enabling move | What GeoVac would forfeit | Verdict | Primary source |
|:--|:--|:--|:--|:--|:--|
| **A1. Herbst–Avery–Dreuw / molsturm** (shared, variationally optimized k) | **ATOMS ONLY.** Li–Ne, Na–Ar at **Hartree–Fock**, vs numerical-HF reference values (their Table I). n_max = 10 gives "4 to 5 digits of accuracy in the HF energy". With no l-cutoff: "reaching below absolute errors of **0.1 Hartree requires beyond 80 basis functions**, making calculations ... rather impractical". Convergence "noticeably sublinear and overall comparatively slow". Correlated MP2/FCI on Be shown in the 2019 talk (molsturm to pyscf), **no error budget published**. **Zero molecules.** | Separate cutoffs (n_max, l_max, m_max); one exponent k common to all functions, taken as a 4th *variational* parameter; k_opt tracks the Clementi–Raimondi average effective charge. The l_max cutoff drops basis growth from cubic to linear in n_max. | Nothing GeoVac has not already paid — this *is* GeoVac's free-scale posing. No molecular result to import. | **STOP** for molecules (nothing demonstrated); BORDERLINE as the atomic state-of-the-art datum | M. F. Herbst, J. E. Avery, A. Dreuw, "Quantum chemistry with Coulomb Sturmians: Construction and convergence of Coulomb Sturmian basis sets at Hartree-Fock level", *Phys. Rev. A* **99**, 012512 (2019), DOI 10.1103/PhysRevA.99.012512; arXiv:1811.05777. **VERIFIED (full text read)** |
| **A2. Avery & Avery isoenergetic many-centre CS / Generalized Sturmian Method** | **None reported.** JPCA 2009 is a method proposal; its abstract carries no energy. Mol. Phys. 2012 states "**pilot calculations** are performed on N-electron molecules". Secondary reports mention ~10 CS per nucleus and LiH(3+) — a **one-electron** two-centre ion — energy-vs-R curves. | Isoenergetic configurations from many-centre CS, analogous to Goscinskian atomic configurations; for **diatomics all relevant integrals are pure functions of s = kR**, so they can be computed once and stored (quadrature-free, parameter-free table). | Nothing — GeoVac already implements this build (Paper 60). No accuracy to import. | **STOP** (no demonstrated accuracy; "pilot" is the authors' own word) | J. Avery, J. Avery, "Can Coulomb Sturmians Be Used as a Basis for N-Electron Molecular Calculations?", *J. Phys. Chem. A* **113**(52), 14565–14572 (2009), DOI 10.1021/jp9040502. **VERIFIED (abstract read via EuropePMC)**. Companion: J. S. Avery, J. E. Avery, "Coulomb Sturmians as a basis for molecular calculations", *Mol. Phys.* **110**, 1593–1608 (2012), DOI 10.1080/00268976.2012.658876 — **UNVERIFIED** (T&F and ResearchGate both 403; abstract via search summary only) |
| **A3. Avery molecular-ERI technology: product resolution at doubled exponent 2k ("exact density fitting")** | **No molecular energies** — integral technology. Reported property: the resolution object is sparse (Herbst: "again sparse (220 MB for n = 20)"). | 1e multicentre integrals = **Shibuya–Wulfman integrals**, analytic in S = k(R_A − R_A'), times a sparse pre-computable coupling coefficient. ERIs: expand the *product* chi_1(r) chi_2(r) in an auxiliary CS set at exponent **2k** by solving one linear system against S(2k, S); then (12&#124;34) = sum C_21 I C_34. Herbst's slide: "**Like exact density fitting**". | Not the quadrature-free axis (GeoVac keeps that) and not the metric-free axis. It **does** forfeit l-selection: a density-fitting resolution mixes l in the auxiliary index, so atomic Gaunt/6j sparsity does not survive the product step. | **BORDERLINE** — nameable, not excluded by any GeoVac wall, but buys *decidability/speed*, not accuracy (exact != accurate). The 2k second scale is the interesting residue. | J. E. Avery, J. S. Avery, "Molecular Integrals for Exponential-Type Orbitals Using Hyperspherical Harmonics", *Adv. Quantum Chem.* **70**, 265–324 (2015), DOI 10.1016/bs.aiq.2014.07.004. **VERIFIED citation (Crossref); mechanism read off Herbst's 2019 Metz slides, which cite this chapter.** Earlier: J. E. Avery, "Fast Electron Repulsion Integrals for Molecular Coulomb Sturmians", *Adv. Quantum Chem.* **67**, 129–151 (2013), DOI 10.1016/b978-0-12-411544-6.00006-6 (VERIFIED citation) |
| **A4. Aquilanti–Coletti momentum-space Sturmians for molecules** | **H2 potential energy curve, "encouraging results already at the level of minimal basis set"** — tens of mHa at best; no error quoted, no correlated reference. Companion work on **H2+** (one electron). | Valence-bond-type expansion in momentum-space Sturmian / hyperspherical harmonics; integrals partly via angular-momentum algebra in momentum space, partly routed through existing **STO integral programs**. | The STO-routing step forfeits the pure-number s = kR table (it reintroduces an external integral engine) and does not restore l-sparsity. | **STOP** (minimal basis; no accuracy demonstrated) | J. Avery, D. Z. Ostrovsky (as indexed), C. Coletti, V. Aquilanti, "Sturmian orbitals and molecular structure", *J. Mol. Struct. THEOCHEM* (2004), ScienceDirect PII S0166128004004725. **UNVERIFIED** (ScienceDirect 403; abstract via two independent search summaries). Related: Aquilanti, Cavalli, Coletti, *Chem. Phys.* **214**, 1 (1997); *Phys. Rev. Lett.* **80**, 3209 (1998) — VERIFIED citations only |
| **A5. Goscinski / Shull–Löwdin / Rotenberg / Klahn–Bingel completeness line** | **No molecular calculation.** Content is negative theory: L2-completeness of a basis is **not sufficient** for energy convergence. A completeness *taxonomy* — overcompleteness / exact completeness / incompleteness — plus an "asymptotic dimension" measure. | n/a (theory of completeness vs convergence) | n/a | **STOP as an accuracy path; GO as a citation debt** for Paper 60's overcompleteness section | B. Klahn, W. A. Bingel, "The convergence of the Rayleigh-Ritz method in quantum chemistry", *Theor. Chim. Acta* **44**, 9–26 (I) and 27–43 (II) (1977), DOI 10.1007/BF00548027 — already in `frames_riesz_overcompleteness_memo.md` ref 42. **NEW, not in the corpus:** B. Klahn, W. A. Bingel, "Completeness and linear independence of basis sets used in quantum chemistry", *Int. J. Quantum Chem.* **11**(6), 943–957 (1977), DOI 10.1002/qua.560110607 — **UNVERIFIED** (Wiley 403; taxonomy vocabulary from search snippet only) |
| **A6. Kereselidze–Chkadua–Defrance spheroidal CS** (flagged, NOT my axis) | not assessed | CS re-expressed in spheroidal coordinates for diatomics | n/a | **hand to the spheroidal-Sturmian agent** | T. Kereselidze, G. Chkadua, P. Defrance, "Coulomb Sturmians in spheroidal coordinates and their application for diatomic molecular calculations", *Mol. Phys.* **113**, 3471–3479 (2015), DOI 10.1080/00268976.2015.1036146. **VERIFIED citation (Crossref); content not read** |

---

## Contradictions / sharpenings against corpus statements

1. **CLAUDE.md §3.5 guardrail row overstates Herbst–Avery–Dreuw.** The row reads
   "...it binds (Avery SW closed forms; Herbst-Avery-Dreuw)". HAD 2019 contains
   **no molecular calculation of any kind** — atoms only, HF only, molecules in
   the Outlook. Paper 8's own Remark (L1005) is correct and careful — "consistent
   with published Coulomb--Sturmian Hartree--Fock calculations (there for atoms;
   the molecular cross-$n$ binding is the escape this Remark establishes)" — so
   this is a **summary-surface defect in CLAUDE.md, not in the paper**, exactly
   the class the §9 Summary-Surface Reading Rule names. Recommend rewording so
   the binding claim is attributed to the corpus's own cross-n overlap
   computation, not to HAD.

2. **HAD 2019 is NOT an external test of `eq:scale_lock`.** Its k is shared and
   *variationally optimized*; the isoenergetic **locked** posing is never tested
   there. Its finding that "convergence is achieved regardless of the value of k"
   and that k affects the *rate* but not the *trend* is consistent with the
   corpus but says nothing about the locked branch. Do not cite it either way on
   the lock.

3. **A hard external number supporting "accuracy = basis size".** HAD: with no
   l-cutoff, second-half second-row atoms need **more than 80 basis functions to
   get below 0.1 Ha**, and convergence is "noticeably sublinear". An independent
   primary-source confirmation of the corpus's "exact != accurate; accuracy is
   max_n" line, at HF level and one centre. Usable in Paper 60.

4. **Klahn & Bingel IJQC 11, 943 (1977) is a citation debt.** The corpus cites
   their *Theor. Chim. Acta* pair but not the IJQC paper, whose reported subject
   is precisely a taxonomy of **overcompleteness / exact completeness /
   incompleteness** with an "asymptotic dimension" measure — the same vocabulary
   as the v5.11.2 "overcompleteness is ONE DIRECTION" finding. Owed: verify the
   abstract (Wiley 403 this pass) and cite or distinguish before Paper 60's
   overcompleteness section is certified.

5. **The doubled exponent 2k is not in the corpus.** The corpus formalizes the SW
   matrix as "the molecular metric V_0" at exponent k. Avery's molecular ERI
   route additionally requires the SW matrix at **2k** as a product-resolution
   metric. Repo-grepped: no occurrence of a 2k auxiliary object. A distinct
   structure, and the only second scale anywhere in the canon.

---

## Most promising untraveled path

**The two-scale question, posed at two centres.** Both canons independently
locate the accuracy floor in the *same* place, and neither has tested the fix in
a molecule. Avery's own explanation of the helium floor (relayed, recorded in
`memory/avery_method_and_prior_art_gaps.md` item 4) is that a Goscinskian
configuration puts every electron of a shell at one exponent, so in–out radial
correlation — which variational treatments buy with a **split-shell 1s1s'** pair
carrying two independent exponents — is unavailable for *any* weighting
potential. Herbst's 2019 talk lists "**Basis sets with multiple ks?**" as an open
question, unanswered. And Avery's molecular ERI machinery already carries a
second scale, **2k**, but only in the *product/auxiliary* space, never in the
orbital space. The field has the diagnosis, has the machinery for a doubled
exponent, and has never put the two together at two centres.

What makes this GeoVac-specific rather than generic: the corpus has already run
the atomic half of the experiment and **it worked** — "free per-shell lambda ...
is complete and DOES work", explicitly distinguished in the failed-approaches
ledger from the per-n hydrogenic scaling k_n = Z/n that plateaued near 60 mHa
because the bound hydrogenic set is incomplete. The atomic verdict is positive;
the molecular verdict is untested. The sharp, cheap question: **does a two-scale
(split-shell) Coulomb-Sturmian basis survive the molecular shared-scale posing,
and at what price?** The tension is explicit and decidable — Avery's diatomic
construction is quadrature-free *precisely because* every integral is a pure
function of the single ratio s = kR, so two orbital scales give two ratios
(k1·R, k2·R) and the pure-number table becomes a two-parameter table. A real
cost, but a bounded, nameable one, not a wall: a 2-D table is still
quadrature-free and still parameter-free once tabulated. Three things to
measure, in order: (i) does the split-shell lever transfer from the atomic
ladder to the diatomic build at all; (ii) what happens to the conditioning the
band-Toeplitz preconditioner just bounded (two scales change the symbol); and
(iii) does the metric-free locked posing survive — the derivation
(beta_mu − beta_nu)·<Phi_mu|V_0|Phi_nu> = 0 assumes one beta per configuration,
which is exactly what a split shell breaks.

**Honest cap.** It attacks the accuracy axis, which is the right axis, but it is
likely to cost the metric-free posing — and the corpus already knows (Paper 60)
that freeing the scale reaches chemical accuracy while reinstating the L2 metric.
The most probable outcome is a *quantified* restatement of the existing trade,
not an escape from it. Worth doing because the quantity — how much accuracy a
second scale buys per unit of metric reinstated, at two centres — is not measured
anywhere, by anyone.

---

## Coverage gaps

Nothing gathered was lost to the interruption; every finding above was already
in context when the run was resumed, and no claim here is reconstructed from
memory. The genuine gaps are sources never reached:

- **Mol. Phys. 110, 1593 (2012)** full text — the single most likely place a
  molecular Sturmian energy is actually printed. T&F and ResearchGate both 403.
- **The two Avery books** (1989 Kluwer; 2006 World Scientific) — not reachable by
  any web route; the standing unverified gap noted in the earlier conditioning
  scan applies here too.
- **THEOCHEM 2004** full text (ScienceDirect 403) — the H2 minimal-basis result
  rests on two independent search summaries, not the abstract at source.
- **Klahn & Bingel, IJQC 11, 943 (1977)** abstract (Wiley 403) — the
  overcompleteness-taxonomy claim is snippet-level only.
- **Herbst PhD thesis** (ref [95] of the PRA paper) — not fetched; may contain
  correlated CS results beyond the talk's Be MP2/FCI plot.

---

## Search trail

**Queries run (25 web calls):**

1. WS Herbst/Avery/Dreuw PRA 2019 + molecules → PRA/arXiv/molsturm hits.
2. WS molsturm + multicenter → surfaced the Avery "Fast ERI for Molecular
   Coulomb Sturmians" chapter.
3. WF arxiv.org/abs/1811.05777 → abstract only; confirmed "second and third row
   atoms".
4. WF michael-herbst.com/research/coulomb_sturmians/ → no numbers.
5. WF arxiv.org/pdf/1811.05777 → returned binary, **but saved locally**;
   extracted with `pdftotext -layout` and searched in full. The load-bearing read.
6. WF michael-herbst.com/talks/2019.05.03_coulomb_sturmians_metz.pdf → binary,
   saved locally, extracted. Source of the "exact density fitting" and "Molecular
   Coulomb Sturmians?" findings.
7. WS Avery many-centre CS H2/LiH energies → found JPCA 2009 and THEOCHEM 2004.
8. WS generalized Sturmian H2 PEC → "encouraging results at minimal basis set".
9. WF pubmed 19807119 → cookie wall, **dead end**.
10. WF sciencedirect S0166128004004725 → **403, dead end**.
11. WF EuropePMC REST `EXT_ID:19807119` → **HIT**: full JPCA 2009 citation +
    abstract. (Reusable: EuropePMC REST bypasses the PubMed cookie wall.)
12. WS "Can Coulomb Sturmians be used..." → corroborated abstract; "automatic
    scaling" advantage claim.
13–14. WF api.semanticscholar.org ×2 → **429 rate-limited, dead end** both times.
15. WS Aquilanti/Coletti THEOCHEM H2 → corroborated minimal-basis framing;
    surfaced the H2+ momentum-space paper.
16. WS Sturmian H2 total energy vs Kołos–Wolniewicz → no Sturmian *molecular*
    energy found; only a Sturmian **H2+** datum and the KW references.
17. WF api.crossref.org (Avery Mol. Phys. 2012) → **HIT**: 5 exact citations
    including the Kereselidze spheroidal paper.
18. WS Herbst thesis molsturm MP2/FCI → no molecular results; confirmed atoms.
19. WS "James Emil Avery" thesis → MSc 2008 "The generalised Sturmian method",
    PhD 2011 on nanostructures; surfaced "10 basis elements per nucleus" and
    **LiH(3+)** (one electron).
20. WS "Coulomb Sturmians as a basis for molecular calculations" → "**pilot
    calculations** are performed on N-electron molecules".
21. WF api.crossref.org (AQC vol. 70) → **HIT**: exact title "Molecular Integrals
    for Exponential-Type Orbitals Using Hyperspherical Harmonics", 265–324,
    DOI 10.1016/bs.aiq.2014.07.004.
22. WS Avery GSM pilot-calculation molecule identity → no neutral molecule named;
    corroborated LiH(3+).
23. WS Weniger / Klahn / Bingel completeness → **HIT**: Klahn & Bingel IJQC 11,
    943–957 (1977) and the overcomplete / exactly-complete / incomplete taxonomy.
24. WF researchgate 254265141 → **403, dead end**.
25. WF Wiley 10.1002/qua.560110607 → **403, dead end**.

**Local (free) work:** read the four pre-existing memos named in the task;
`pdftotext` on the two locally-saved PDFs (this produced every verbatim quote in
A1 and A3); grepped the repo for Klahn/Bingel, "asymptotic dimension", the 2k
object, the Herbst citation contexts (Paper 8 L1005, Paper 60 L288/L729), and the
per-n-scaling ledger row.

**Access notes for the next scanner.** ScienceDirect, Wiley, Taylor & Francis and
ResearchGate all return 403 to WebFetch; PubMed returns a cookie wall; Semantic
Scholar rate-limits aggressively. **Crossref REST and EuropePMC REST both work**
and should be the default for citation and abstract retrieval. arXiv and
personal-site PDFs come back as binary through WebFetch **but are saved to the
tool-results directory**, where `pdftotext -layout` recovers them fully — the
single highest-value trick of this pass, and what turned two "unreadable" fetches
into the two primary sources.
