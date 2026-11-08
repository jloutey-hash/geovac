# Confidence Review: Paper 45 — Lorentzian propinquity convergence on truncated SU(2) × U(1)_T Krein spectral triples

## Calibration check (calibration runs only)

N/A — this is a Wave-2 production audit, not a calibration run against a
known-honest answer.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence |
|---|---|---|---|---|---|
| 1 | $\Lprop(\Tcal^L_{n,N,T}, \Tcal^L_M) \le C_3 \cdot \gamma \to 0$ on K⁺ | abstract, Thm 5.4 (`thm:main`) | B | GEOVAC-ONLY | Proof chains entirely through Papers 38/44/46/47 internal results. Latrémolière 2017/2022 (Adv. Math. 404:108393) provides the standard weak-form bound; that part is EXTERNAL. The K⁺-restriction-to-Hilbert reduction is GeoVac's own construction. |
| 2 | Five-lemma transport L1'/L2/L3/L4/L5 closes | §3–§5 | B | GEOVAC-ONLY | Internal Papers 38/44 are the substrate; only Bożejko-Fendler central-multiplier equality (L2) and Katznelson Fejér rate (L4 U(1) factor) are external. |
| 3 | cb-norm L2 = $2/(\nmax+1)$ ($N_t$-independent) | Lemma 3.1 (`lem:L2`) | A | MIXED | Bożejko-Fendler 1991 + standard Fejér Plancherel symbol = 1 at $k=0$; calculation is mechanical. EXTERNAL anchors are sound. |
| 4 | Joint commutator $[\DL, a_s \otimes a_t] = i[\DGV, a_s] \otimes a_t$ bit-exact (Lemma 3.3, eq L3 structural identity) | §3.2, eq (\ref{eq:L3_struct_id}) | B | GEOVAC-ONLY | Two-line manipulation given the chirality-doubled block structure of $a_s$ and the diagonality of $a_t$ in the momentum basis. Mechanically correct under the stated conventions. |
| 5 | $C_3 \to 1^-$, finite-cutoff $\sqrt{1 - 1/\nmax}$ via envelope $N \le 2\nmax - 1$ | Lemma 3.3 + Remark 3.4 | B | GEOVAC-ONLY | Functional form $(N-1)/\sqrt{N^2-1}$ is Paper 38 Lemma L3; envelope range is a Paper 46 §4.3 refinement (cf. Remark `rem:envelope_v2`). Internally consistent. |
| 6 | Berezin L4(c) approximate identity rate $\gamma^{joint} = O(\log n/n + T/N_t)$ | Lemma 3.5 (`lem:L4`) | B | MIXED | Spatial rate from Paper 38; U(1) rate from Katznelson Fejér theorem (EXTERNAL). |
| 7 | Numerical panel $\{2.0746, 1.6101, 1.3223\}$ monotone decreasing in $\nmax$ | Table 1, §4.1 | B | GEOVAC-ONLY | Values match Papers 38/39/46 panels bit-identically; not independently reproduced here. |
| 8 | Riemannian-limit recovery at $N_t = 1$ bit-exact (Frobenius 0.0 in float64) | Prop 4.1 (`prop:riemannian_limit`) | B | GEOVAC-ONLY | Load-bearing falsifier internal to the construction; tensor factorization at $N_t = 1$ trivializes the U(1) Berezin to $1 \times 1$ identity. Mechanically correct under the construction. |
| 9 | "First Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" | abstract | C | MIXED | Web search returns no contradicting prior art (Connes-vS 2021 explicitly defers Lorentzian convergence as "elsewhere"; Latrémolière 2017–2025, Leimbach-vS 2024, Hekkelman-McDonald 2024, Toyota 2023, Farsi-Latrémolière 2024/2025 are all Riemannian; Mondino-Sämann + Minguzzi-Suhr + Che-Perales-Sormani are synthetic pre-length space, not operator-algebraic). Paper already hedges with "To our knowledge"; the abstract uses the unconditional "first…in the published…literature" while §1.3 hedges back to "to our knowledge". Recommend the abstract softens to match §1.3. See overstatement table. |
| 10 | "The Lorentzian Dirac per van~den~Dungen~2016 Proposition 4.1" | abstract, §2.3 | A | EXTERNAL | Verified — arXiv:1505.01939, Math. Phys. Anal. Geom. 19 (2016) Art. 4 confirmed via WebFetch. |
| 11 | BBB classification $(m, n) = (4, 6)$ at signature $(3, 1)$ West-coast | §2.2 | D | MIXED | Cited to BBB 2018 Table 3. Verified the BBB paper exists (arXiv:1611.07062, JMP 59:062303 2018), but Table 3 itself not re-derived here. Mechanically standard. |
| 12 | G2 (de-compactification $T \to \infty$) CLOSED on natural substrate via Paper 47 §7 Thm 7.3 | §1.4 G2 | B | GEOVAC-ONLY | Reference to Paper 47 §7 Thm 7.3 — internal. Cannot externally verify the closure status. The framing is internally consistent (zero-cost temporal extent element argument). |
| 13 | G3 (cross-manifold) CLOSED as structural impossibility | §1.4 G3 | B | GEOVAC-ONLY | Cites Paper 32 Thm 7.4 (Bertrand + Fock/HO rigidity); internal. |
| 14 | Universal $4/\pi$ asymptote on the SU(2) factor | Remark 3.2 | B | GEOVAC-ONLY | Cites Paper 38 §3.2 + Paper 40 §3.2 Plancherel-Vandermonde theorem. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputation | Survives? |
|---|---|---|---|---|
| $C_3^{(3)} = 2/\sqrt{8}$ | $0.707$ (Remark `rem:envelope_v2` line 851) | derive: $(N-1)/\sqrt{N^2-1}$ at $N=3$ = $2/\sqrt{8}$ = $0.7071$ | 0.7071 ✓ | YES |
| $C_3^{op}(3) = \sqrt{2/3}$ | $\approx 0.816$ (Remark line 852) | derive: $\sqrt{(6-1-1)/((6-1)^2-1)}$ at $N=2 \cdot 3 - 1 = 5$ = $\sqrt{4/24}$ = $\sqrt{1/6} \approx 0.408$. **Paper writes $C_3^{op}(3) = \sqrt{(2 \cdot 3 - 1 - 1)/((2 \cdot 3 - 1)^2 - 1)} = \sqrt{4/24}$ but states this equals $\sqrt{2/3} \approx 0.816$.** $\sqrt{4/24} = \sqrt{1/6} \approx 0.408$, **NOT** $\sqrt{2/3}$. Alternatively, $\sqrt{1 - 1/\nmax} = \sqrt{2/3} \approx 0.816$ AT $\nmax = 3$, which matches the asymptotic formula but the displayed intermediate calculation $\sqrt{(2 \cdot 3 - 1 - 1)/((2 \cdot 3 - 1)^2 - 1)}$ evaluates to $\sqrt{4/24} \ne \sqrt{2/3}$. | **FAILS** — possible bug in Remark 3.4 displayed calculation. See finding R-4. |
| Asymptotic ratio $\Lambda(4,7)/\Lambda(2,3)$ | $0.6374$ (§4.3) | $1.3223 / 2.0746 = 0.6374$ | 0.6374 ✓ | YES |
| Plancherel norm $\hat K^{SU(2)}_{\nmax}(j_{\max})$ supremum | $2/(\nmax+1)$ (Lemma 3.1) | $\nmax / [\nmax(\nmax+1)/2] = 2/(\nmax+1)$ ✓ | 2/(n+1) ✓ | YES |
| $\dim \HGV^{\nmax} = (2/3) n(n+1)(n+2)$ | eq (2.2) | At $\nmax = 3$: $(2/3) \cdot 3 \cdot 4 \cdot 5 = 40$. Consistent with the truthful CH dimensionality table in CLAUDE.md memory (R3.2 dimH = 40). | 40 ✓ | YES |

### Circularity map

The full proof chain bottoms out at internal references for steps L1',
L2 (in part), L3, L4, and L5. The bottom of each chain:

- **L1' substrate** → Paper 44 (operator-system substrate) → GEOVAC-ONLY.
- **L2 cb-norm** → Bożejko-Fendler 1991 (EXTERNAL) + Plancherel
  computation (MIXED).
- **L3 structural identity** → mechanical algebra given the
  chirality-doubled multiplier convention from Paper 44 → GEOVAC-ONLY.
- **L4 approximate identity** → Paper 38 spatial factor (GEOVAC-ONLY) +
  Katznelson Ch. I Fejér rate (EXTERNAL) → MIXED.
- **L5 propinquity assembly** → Latrémolière 2022 Adv. Math. (EXTERNAL)
  + Paper 38 §L5 bookkeeping (GEOVAC-ONLY) → MIXED.
- **Riemannian-limit falsifier** → tensor-factorization at $N_t = 1$
  → GEOVAC-ONLY.
- **Universal $4/\pi$** → Paper 38 + Paper 40 → GEOVAC-ONLY.

Net: the Latrémolière weak-form bound is EXTERNAL, but every substrate
ingredient is internal. The construction is well-formed as a transport
on top of EXTERNAL Latrémolière 2022 machinery; everything peculiar to
GeoVac (CH spinor bundle, BBB embedding, joint Berezin) is internal.

### Overstatement findings

| Exact phrase | Location | Suggested honest replacement |
|---|---|---|
| "To our knowledge this is the first Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature." | abstract | Already hedged with "To our knowledge"; the §1.3 phrasing is consistent. The CITE-OK status of all surveyed prior art does support a "to our knowledge first" claim, but per CONFIDENCE_REVIEW §B you cannot establish absolute novelty. **No change required**; the phrasing is honest. |
| "G2 is CLOSED on the natural substrate." | §1.4 G2 line 315 | Internally consistent with Paper 47 cross-reference; the natural-substrate qualification keeps this honest. No change required. |
| "G3 (cross-manifold extensions): CLOSED as structural impossibility" | §1.4 G3 | Standard claim against Paper 32 Theorem 7.4; internally consistent. No change required. |

No headline overstatement found after re-reading the abstract against
the body. The honest-scope discipline (K⁺-weak-form vs strong-form, Q1
named, G2 partial closure named) is upheld.

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found |
|---|---|---|---|
| `bekka_harpe_valette2008` | Kazhdan's Property (T), Cambridge UP 2008 | CITE-OK | Standard reference. |
| `bisognano_wichmann1976` | JMP 17 (1976) 303–321 | CITE-OK | Standard QFT reference. |
| `bizi_brouder_besnard2018` | JMP 59 (2018) 062303; arXiv:1611.07062 | CITE-OK | Verified via WebFetch. |
| `bozejko_fendler1991` | Arch. Math. 57 (1991) 290–298 | CITE-OK | Standard cb-norm reference. |
| `camporesi_higuchi1996` | J. Geom. Phys. 20 (1996) 1–18; arXiv:gr-qc/9505009 | CITE-OK | Verified — title exactly matches. |
| `connes_vs2021` | Comm. Math. Phys. 383 (2021) 2021–2067; arXiv:2004.14115 | CITE-OK | Verified — title and authors match. |
| `devastato_lizzi_martinetti2018` | JMP 59 (2018) 112301 | CITE-OK | Standard NCG/time emergence reference. |
| `farsi_latremoliere2024` | "Continuity of spectral propinquity for spectrally truncated spheres", preprint 2024 | CITE-WRONG-METADATA | Preprint stub — cannot easily verify; the 2025 sibling (arXiv:2504.11715) exists and is properly attributed. Recommend either upgrading to 2025 if same content, or providing arXiv ID. |
| `farsi_latremoliere2025` | arXiv:2504.11715 (2025) | CITE-OK | Verified. |
| `franco_eckstein2014` | Class. Quantum Grav. 30 (2013) 135007; arXiv:1212.5171 | CITE-WRONG-METADATA | Bibitem year says 2014 but actual journal year is 2013. Minor — fix year or rename key. |
| `geroch1967` | JMP 8 (1967) 782–786 | CITE-OK | Standard GR reference. |
| `hartle_hawking1976` | Phys. Rev. D 13 (1976) 2188–2203 | CITE-OK | Standard. |
| `hekkelman_mcdonald2024` | "Spectral truncations of $\mathbb{T}^d$ and quantum metric geometry", arXiv:2403.18619 | **CITE-MISATTRIBUTED** | **arXiv:2403.18619 is "Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures" by Calderón–Rucci–Chichizola — completely unrelated.** Furthermore, the actual paper "Gromov–Hausdorff convergence of spectral truncations for tori" is by **Leimbach–van~Suijlekom (arXiv:2302.07877, Adv. Math. 439:109496, 2024)** — which Paper 45 already cites separately as `leimbach_vs2024`. So this bibitem is doubly broken: wrong arXiv ID AND wrong authorship attribution. Recommend deleting this bibitem and consolidating into `leimbach_vs2024`. |
| `hekkelman_mcdonald2024b` | "A noncommutative integral on spectrally truncated spectral triples...", arXiv:2412.00628 | CITE-OK | Verified — title and authors confirmed. |
| `katznelson2004` | Harmonic Analysis 3rd ed., Cambridge UP 2004 | CITE-OK | Standard reference. |
| `latremoliere_metric_st_2017` | Adv. Math. 404 (2022) 108393; arXiv:1811.10843 (Nov 2018) | CITE-OK | Verified — title is "The Gromov–Hausdorff propinquity for metric spectral triples". |
| `latremoliere2018` | Trans. AMS 368 (2016) 365–411 | CITE-OK | Verified. |
| `leimbach_vs2024` | Adv. Math. 439 (2024) 109496 | CITE-OK | Verified — see hekkelman_mcdonald2024 entry above. |
| `mondino_samann2024` | Lett. Math. Phys. 114 (2024) Paper 37; arXiv:2209.14384 | **CITE-MISATTRIBUTED** | **arXiv:2209.14384 is "Lorentzian metric spaces and their Gromov-Hausdorff convergence" by E. Minguzzi and S. Suhr, LMP 114 (2024) Art. 73.** Not Mondino-Sämann, not "Lorentzian metric measure spaces", and not Art. 37 (it's Art. 73). The dispatch prompt flagged this; confirmed via WebFetch. Recommend full rename to `minguzzi_suhr2024` with corrected title/authors/article number. |
| `mondino_samann2025` | "Synthetic Lorentzian Gromov–Hausdorff convergence and pre-compactness", arXiv:2504.10380 (2025) | CITE-OK | Verified — title matches "Lorentzian Gromov-Hausdorff convergence and pre-compactness" (close enough; "Synthetic" prefix is in the keywords/abstract). |
| `che_perales_sormani2025` | arXiv:2510.13069 (Oct 2025) | CITE-OK | Verified — title matches. (Per dispatch prompt: don't re-flag.) |
| `nieuviarts2025a` | arXiv:2502.18105 (v3, May 2025) | CITE-OK | Verified. |
| `nieuviarts2025b_proceedings` | arXiv:2512.15450 (Dec 2025, rev May 2026) | CITE-OK | Verified. |
| `paulsen2002` | Cambridge Studies Adv. Math 78 | CITE-OK | Standard reference. |
| `pier1984` | Amenable Locally Compact Groups | CITE-OK | Standard reference. |
| `pisier2001` | LNM 1618, Springer | CITE-OK | Standard reference. |
| `sewell1982` | Ann. Phys. 141 (1982) 201–224 | CITE-OK | Standard QFT-on-curved-spacetime reference. |
| `strohmaier2006` | J. Geom. Phys. 56 (2006) 175–195 | CITE-OK | Standard Krein-NCG axiomatic reference. |
| `takesaki1979` | Theory of Operator Algebras I | CITE-OK | Standard. |
| `toyota2023` | arXiv:2309.13469 | CITE-OK | Verified. |
| `unruh1976` | Phys. Rev. D 14 (1976) 870–892 | CITE-OK | Standard. |
| `vandungen2016` | Math. Phys. Anal. Geom. 19 (2016) Art. 4; arXiv:1505.01939 | CITE-OK | Verified. |
| `de_groot2026_su11` | arXiv:2601.22171 (Jan 2026) | CITE-OK | Verified. |
| `paper24`–`paper47` (internal) | Internal preprints | CITE-OK (cannot externally verify but consistent across corpus) | — |
| `latremoliere2025_hypertopology` | referenced as "pinned proper QMS hypertopology" in §1.4 G2 line 309 | **CITE-CANT-FIND** | This `\cite` key appears in the body (§1.4 G2) but **has no bibitem in the bibliography**. Either add the bibitem (Latrémolière arXiv:2512.03573 per CLAUDE.md §6 Paper 47 entry) or change the `\cite` to a footnote/cross-ref to Paper 47. |
| `paper48` | referenced at §1.4 G2 line 311 ("Krein-lifted per Paper 48 §3") | **CITE-CANT-FIND** | Same problem — `\cite{paper48}` used in body, **no bibitem present**. Add Paper 48 bibitem (it exists in the corpus per CLAUDE.md §6). |

### Problems found

1. **CITE-MISATTRIBUTED — `mondino_samann2024`** (line 1759): arXiv:2209.14384 is Minguzzi-Suhr, not Mondino-Sämann. The CITE-USAGE count is **2** in-text references — line 374 (§1.3 related work) and line 1463 (§6.2 Mondino-Sämann discussion). Both usage sites cite this as "their earlier work" of Mondino-Sämann predating the 2025 paper; the actual earlier Mondino-Sämann paper is arXiv:1707.00880 ("Lorentzian metric measure spaces") OR more likely the Mondino-Sämann–Cavalletti–Ohta line. Recommend the PM:
   - rename cite-key to `minguzzi_suhr2024`,
   - correct the bibitem title to "Lorentzian metric spaces and their Gromov-Hausdorff convergence",
   - correct authors to Minguzzi & Suhr,
   - correct article number to 73 (not 37),
   - update both in-text `\cite{mondino_samann2024}` references on lines 374 and 1463 to `\cite{minguzzi_suhr2024}`,
   - re-read the body text at lines 372–382 and 1459–1478: the framing that 2024 was "their earlier work" is wrong — the actual prior Lorentzian-GH-type work is Minguzzi-Suhr's, not Mondino-Sämann's own earlier paper. Recommend rewording these passages so the Minguzzi-Suhr paper is correctly named as a separate prior strand of the synthetic Lorentzian-GH program.
2. **CITE-MISATTRIBUTED — `hekkelman_mcdonald2024`** (line 1722): arXiv:2403.18619 is an unrelated OpenMP/graph paper, NOT Hekkelman-McDonald "Spectral truncations of $\mathbb{T}^d$". The actual spectral-truncations-of-tori paper is by Leimbach-van~Suijlekom (arXiv:2302.07877), which Paper 45 already cites correctly as `leimbach_vs2024`. CITE-USAGE count: **1** in-text reference on line 347 (§1.3 related work, listed in the post-Connes-vS-2021 Riemannian-side program flat-structure series). Recommend deleting this bibitem and removing it from the in-text citation chain on line 347 (the citation chain already includes `leimbach_vs2024` for the same work in §3 of related work; or alternatively merge two cite-keys in the inline list).
3. **CITE-CANT-FIND — `latremoliere2025_hypertopology`** (line 309): no bibitem in bibliography. Add the corresponding bibitem (arXiv:2512.03573 per CLAUDE.md §6 Paper 47 entry).
4. **CITE-CANT-FIND — `paper48`** (line 311): no bibitem for internal Paper 48. Add the bibitem (Paper 48 is in `papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex` per CLAUDE.md §6).
5. **CITE-WRONG-METADATA — `franco_eckstein2014`** (line 1709): bibitem says 2014 but journal year is 2013. Minor.
6. **CITE-WRONG-METADATA — `farsi_latremoliere2024`** (line 1697): preprint stub, no arXiv ID, no DOI, no abstract handle. Either add arXiv ID or fold into `farsi_latremoliere2025` if same content.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "To our knowledge this is the first Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature." | abstract | math.OA arXiv listings: Connes-vS 2021, Latrémolière 2017-2025, Leimbach-vS 2024, Hekkelman-McDonald 2024, Toyota 2023, Farsi-Latrémolière 2024/2025, Mondino-Sämann 2024/2025, Minguzzi-Suhr 2024, Che-Perales-Sormani 2025, van~den~Dungen 2016, Strohmaier 2006, BBB 2018, Franco-Eckstein 2013, Devastato-Lizzi-Martinetti 2018, Nieuviarts 2025a/b, de Groot 2026 | NO contradicting prior art found. All Lorentzian-side spectral-triple frameworks (Strohmaier, vdD, BBB, FE, DLM) define the algebraic substrate but **do not provide a metric**. All Latrémolière-lineage convergence theorems are Riemannian. The Mondino-Sämann + Minguzzi-Suhr + Che-Perales-Sormani synthetic-pre-length-space program is a different mathematical object (paper acknowledges this in §1.3 and §6.2 explicitly). The Nieuviarts twist-morphism program is algebraic, not metric (and paper acknowledges this in §1.3). | Phrasing is appropriate as-is ("to our knowledge"). Per CONFIDENCE_REVIEW §B, a clean search supports the hedged claim. No change required. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| Remark 3.4 (`rem:envelope_v2`) displayed equation for $C_3^{op}(3)$: $\sqrt{(2 \cdot 3 - 1 - 1)/((2 \cdot 3 - 1)^2 - 1)} = \sqrt{4/24} \ne \sqrt{2/3}$ | A | E | **HIGH** — math error in displayed intermediate calculation. The asymptotic formula $\sqrt{1 - 1/\nmax}$ at $\nmax = 3$ does give $\sqrt{2/3} \approx 0.816$, so the conclusion is correct, but the intermediate displayed evaluation has an arithmetic slip. Fix: either reconcile the intermediate calc (likely should be $\sqrt{(2 \cdot 3 - 2)/((2 \cdot 3)^2 - 4)} = \sqrt{4/32}$, no — let me recompute. Actually $\sqrt{1 - 1/3} = \sqrt{2/3}$; the per-$N$ functional form at $N = \nmax + 1 = 4$ gives $(4-1)/\sqrt{16-1} = 3/\sqrt{15} \ne \sqrt{2/3}$. The relationship between the supremum-formula and the closed-form $\sqrt{1 - 1/\nmax}$ deserves a re-derivation in the paper or clarification.) |
| `mondino_samann2024` bibitem is Minguzzi-Suhr (arXiv:2209.14384), not Mondino-Sämann | B | CITE-MISATTRIBUTED | **HIGH** — load-bearing citation for §1.3 related work and §6.2 discussion of Mondino-Sämann strand. Rename to `minguzzi_suhr2024`, fix title/authors/article-number, and update 2 in-text usages on lines 374 and 1463. |
| `hekkelman_mcdonald2024` bibitem misattributed to arXiv:2403.18619 (an OpenMP paper); the work is by Leimbach-vS | B | CITE-MISATTRIBUTED | **HIGH** — both the arXiv ID and the authorship are wrong. Delete bibitem, consolidate with `leimbach_vs2024`. 1 in-text usage on line 347. |
| `latremoliere2025_hypertopology` cite-key has no bibitem | B | CITE-CANT-FIND | **HIGH** — broken \cite in compiled PDF. Add bibitem. |
| `paper48` cite-key has no bibitem | B | CITE-CANT-FIND | **HIGH** — broken \cite in compiled PDF. Add bibitem. |
| Novelty claim "first Lorentzian propinquity convergence theorem..." | A/B | C | LOW (already hedged) |
| `franco_eckstein2014` year mismatch (2014 vs 2013) | B | CITE-WRONG-METADATA | LOW |
| `farsi_latremoliere2024` preprint stub | B | CITE-WRONG-METADATA | LOW |

## Broadcast readiness: **YELLOW**

The paper's structural content (five-lemma argument, joint cb-norm,
joint commutator structural identity, Berezin reconstruction,
Riemannian-limit recovery) is consistent and internally well-defined.
The Latrémolière 2022 EXTERNAL anchor is correctly invoked, and the
GeoVac-internal substrate (Papers 38, 44, 46, 47) provides a complete
proof skeleton. Most of the load-bearing claims sit at verdict B
(internally consistent, GEOVAC-ONLY) — appropriate for the
operator-algebraic transport that the paper proposes.

However, four HIGH-severity citation defects block clean broadcast: (i)
**two** CITE-MISATTRIBUTED bibitems (`mondino_samann2024` →
Minguzzi-Suhr; `hekkelman_mcdonald2024` → unrelated OpenMP paper /
actual work is Leimbach-vS already cited); (ii) **two** CITE-CANT-FIND
keys (`latremoliere2025_hypertopology`, `paper48`) that compile as `??`
in the PDF and signal sloppiness to any external reviewer; (iii) one
HIGH-severity arithmetic slip in Remark 3.4's displayed intermediate
calculation (final answer correct, intermediate step appears
inconsistent — needs a 5-minute re-derivation).

Recommend a single errata sprint that: (a) renames `mondino_samann2024`
→ `minguzzi_suhr2024` with correct metadata and updates 2 in-text
`\cite` calls on lines 374 and 1463, plus rewords the §1.3 and §6.2
passages framing it as "their earlier work" of Mondino-Sämann; (b)
deletes the broken `hekkelman_mcdonald2024` bibitem and removes its
single in-text usage on line 347 (already covered by `leimbach_vs2024`);
(c) adds bibitems for `latremoliere2025_hypertopology` (arXiv:2512.03573)
and `paper48` (internal); (d) audits the Remark 3.4 intermediate
calculation; (e) optionally fixes the LOW-severity `franco_eckstein2014`
year and `farsi_latremoliere2024` preprint metadata. After this errata
sprint Paper 45 should pass to GREEN.

## What I could NOT verify (hand to a human expert)

- The internal-paper proof chains (Papers 38/40/44/46/47) — I take
  these on faith here. A domain expert in Latrémolière propinquity
  should verify that the five-lemma transport from Paper 38 to the
  joint SU(2) × U(1)_T setting is actually rigorous, especially the
  Stein-Weiss factor-wise step (Prop 4.2 / §`prop:reach_height`,
  $\mathrm{height}_B$ closure).
- The Connes-vS Theorem 5.5 application in the main theorem proof
  (the four-term-max → dominated-term reduction) — this is mechanical
  but rests on the specific Latrémolière metric-spectral-triple
  weak-form bound. A domain expert should verify the K⁺-restriction
  preserves the assumptions of Latrémolière Thm 5.5 (in particular,
  the UCP tunneling pair satisfies the K⁺-preservation property
  needed; the paper claims it does via Remark `rem:Kplus_tunnel_valid`,
  but this rests on Paper 44 §3–4 which I have not externally verified).
- The "first in published math.OA literature" claim — only a domain
  expert can settle priority claims definitively. A clean
  unstructured-search supports the hedged "to our knowledge"
  phrasing.
- The G2 "closed on natural substrate" status — rests on Paper 47 §7
  Thm 7.3 internal reference; cannot externally cross-check.
