# Confidence Review: Paper 50 — Boundary CFT$_3$-on-$S^3$ observables on the truncated GeoVac spectral triple

## Calibration check
Not a calibration run. Wave-3 audit; prior runs found Paper 50 cleaner than the 44–49 cluster (YELLOW in Sweep A).

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| A1 | $F_s = \tfrac{\log 2}{8} - \tfrac{3\zeta(3)}{16\pi^2}$ matches KPS bit-exactly for conformally coupled scalar on $S^3$ | Thm 3.1, abstract | **A** | EXTERNAL (KPS Eq. 41) | Computed both: bit-identical, $0.0638070547761982\ldots$. KPS PDF p. 11 visually confirmed. |
| A2 | $F_D = \tfrac{\log 2}{4} + \tfrac{3\zeta(3)}{8\pi^2}$ matches KPS for "free Weyl Camporesi-Higuchi Dirac on $S^3$" | Thm 3.2 | **A** (value) / **C** (terminology) | EXTERNAL (KPS Eq. 48) | KPS Eq. 48 gives this formula for the **complex Dirac** ("free massless complex Dirac fermion"). KPS d=3 dim_γ=2; total $F_D \approx 0.219$. Paper 50 labels this "Weyl (single chirality)" — in d=3 (odd) there is no Weyl chirality; the value is right but the label is non-standard. |
| A3 | Dual-basis projection: $F_D + 2 F_s = \log(2)/2$ (M2) and $F_D - 2 F_s = 3\zeta(3)/(4\pi^2)$ (M3) | Thm 3.5 | **A** | EXTERNAL (algebraic identity given A1+A2) | Symbolic check: both bit-identical at 50 dps. Trivial corollary. |
| A4 | "First verification of master Mellin engine on non-GeoVac-internal physics observables" | Rem 3.6 | **C** | GEOVAC-ONLY (Paper 18 §III.7 + corpus history) | Novelty claim — properly hedged in body ("we cannot verify novelty" caveat applies). The structural reading is plausible but the M1/M2/M3 classification IS the GeoVac framework; only the **target observable** (KPS values) is external. Soften "first verification" → "first verification known to us on independently published continuum CFT data". |
| A5 | $S(\rho_W) \sim 2\log n_{\max} + (\coth 1 - \log(2\sinh 1))$, constant $= 0.458448\ldots$ | Eq. (3.6) | **A** | EXTERNAL (direct computation) | Independently computed $S(\rho_W)$ at $n_{\max} \in \{2,5,10,50,100,500,2000\}$; asymptotic residual decays as expected (~$10^{-5}$ at $n_{\max}=2000$). Constant matches to 50 dps. PSLQ relation $[-1,1,-1,-1]$ checks symbolically. |
| A6 | $S^3$ Wedge boundary dim = 2 (slope of log scaling) | Sec 4 | **A** | EXTERNAL | Follows from $g(2k+1) \sim n_{\max}^2$ shell degeneracy; matches direct computation. |
| A7 | "F coefficient does NOT appear in $S(\rho_W)$" (Prop 4.3) | Prop 4.3 | **B** | GEOVAC-ONLY (PSLQ on the framework's own constants) | The PSLQ disjointness statement is internally consistent. The structural claim is stronger than what's proved — PSLQ at coeff ceiling $10^8$ rules out simple relations but not arbitrary ones. Wording "provably outside" is overstated → "outside the M2 ⊕ M3 ring at PSLQ ceiling $10^8$". |
| A8 | $K_\alpha = J_{\rm polar}$ = continuum CHM at finite cutoff; convergence rate $\Lambda \le C_3 \gamma_{n_{\max}}$ | Thm 5.1 | **B** | GEOVAC-ONLY (Paper 38 + Paper 42 internal results) | Convergence rate inheritance is internal. The CHM operator-level identification is structural / formal, not bit-exact verification; "identifies" wording is honest. |
| A9 | **$F_s^{S^5} = -\tfrac{\log 2}{32} - \tfrac{\zeta(3)}{32\pi^2} + \tfrac{15\zeta(5)}{64\pi^4} \approx -0.02297$** | Thm 7.1 | **E** (WRONG — 4× the KPS value) | claims EXTERNAL match to KPS Table 1 but FAILS | KPS Table 1 d=5: $F_S = -\tfrac{1}{256}\bigl(2\log 2 + \tfrac{2\zeta(3)}{\pi^2} - \tfrac{15\zeta(5)}{\pi^4}\bigr) \approx -5.74\times 10^{-3}$. Paper 50's value is **exactly 4×** larger. Root cause: Paper 50's stated multiplicity $(2n+4)(n+1)(n+2)(n+3)/6$ should be $/24$ (KPS Eq. 37). The /6 formula is 4× too large at every $n$. The PSLQ found a relation matching the wrong-spectrum value; the "bit-exact match" header claim of §7 (line 709: "The bit-exact match of Sec. 3 extends to round $S^5$") is **FALSE**. |
| A10 | $D'_{\rm Dirac}(0)^{S^5} = -\tfrac{3\log 2}{128} - \tfrac{5\zeta(3)}{128\pi^2} - \tfrac{15\zeta(5)}{256\pi^4} \approx -0.02163$ | Thm 7.2 | **A** (value) / **C** (header) | MIXED — matches KPS per-dim-γ, not total | Matches KPS Table 2 d=5 **divided by dim γ = 4** (i.e., the per-component value $\approx -2.16\times 10^{-2}$ in Table 2 listing). Total KPS $F_D \approx -0.0865$. So the Paper 50 value is one-quarter of the KPS total. The Dirac "Weyl/single chirality" framing in d=5 is closer to defensible (d=5 has Weyl spinors with chirality matrix) but the absolute factor-vs-total deserves explicit reconciliation. |
| A11 | Log 3 cancellation on $S^5$: coefficient $(81 - 90 + 9)/16 = 0$ | Thm 7.3 | **A** | EXTERNAL (algebraic identity) | Bit-exact zero. Confirmed. |
| A12 | Dual-basis projection does NOT extend to $S^5$: 3-dim M-ring vs 2-dim (scalar, Dirac) plane | Prop 7.4 | **A** (linear algebra is correct) / **B** (interpretation rests on A9–A10) | EXTERNAL algebra / GEOVAC interpretation | The over-determined linear system is correctly stated. But the propositions about which transcendentals appear depend on A9–A10; A9's 4× error doesn't change ring membership, so Prop 7.4 conclusion stands. |
| A13 | $S^7$ scalar PSLQ-fails-on-simple-ring negative result | Rem 7.5 | **D** | not verified here | Plausible (multi-basis sweep with named bases). I did not reproduce. |
| A14 | Bulk-side blocked (Prop 4.5 RT, Prop 6.1 AdS/H^4) | Sec 6 | **A** (honest scope statement) | EXTERNAL | Properly framed as "blocked at framework's existing infrastructure", with three named obstructions. |
| A15 | Henningson-Skenderis holographic Weyl anomaly $\leftrightarrow$ Connes-Chamseddine spectral action: both in M2 | Sec 8.4 paragraph | **C** | MIXED | "Structural correspondence" reading is reasonable; both quantities involve heat-kernel asymptotics. But "not a numerical coincidence but a structural inheritance" overstates — it's an architectural commonality, not a derived theorem. |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value/error | survives? |
|---|---|---|---|---|
| $F_s^{S^3}$ | $0.063807$ | KPS Eq. 41: $(1/16)(2\log 2 - 3\zeta(3)/\pi^2)$ | $0.063807054776198\ldots$ | YES — bit-identical |
| $F_D^{S^3}$ | $0.218959$ | KPS Eq. 48: $\log(2)/4 + 3\zeta(3)/(8\pi^2)$ | $0.218959480727576\ldots$ | YES — bit-identical |
| Dual basis $F_D + 2F_s$ | $\log(2)/2$ | algebraic | $0.346573590279972\ldots$ | YES |
| Dual basis $F_D - 2F_s$ | $3\zeta(3)/(4\pi^2)$ | algebraic | $0.091345371175179\ldots$ | YES |
| $C_\infty = \coth 1 - \log(2\sinh 1)$ | $0.458448$ | direct mpmath | $0.458448743368190\ldots$ | YES |
| Entropy asymptotic $S(\rho_W) - 2\log n_{\max} - C_\infty$ at $n_{\max}=2000$ | "$O(1/n_{\max})$" | direct computation | $-1.84\times 10^{-5}$ | YES (consistent with $O(1/n_{\max})$) |
| Log 3 cancellation coefficient | $0$ | $(81-90+9)/16$ | $0$ | YES |
| **$F_s^{S^5}$** | **$-0.02297$ (= KPS d=5 Table 1)** | KPS Table 1 d=5: $-1/256\bigl(2\log 2 + 2\zeta(3)/\pi^2 - 15\zeta(5)/\pi^4\bigr) \approx -5.74\times 10^{-3}$ | $-5.7430\times 10^{-3}$ | **NO — Paper 50 is exactly 4× KPS** |
| $F_D^{S^5}$ | $-0.02163$ (claimed "= KPS d=5") | KPS Table 2 d=5 total (dim_γ=4): $-0.08651$; per dim_γ: $-0.02163$ | matches per-dim-γ, NOT total | Partial: the formula matches the *per-component* KPS entry, not the *total* $F_D$ |

### Circularity map

**GEOVAC-ONLY chains (load-bearing):**
- **A4 / Rem 3.6 (master Mellin engine novelty):** M1/M2/M3 classification is Paper 18 §III.7; the categorization of $\log 2$ as M2 (via $\zeta'_R(0)$ heat-kernel asymptotic) and $\zeta(3)/\pi^2$ as M3 (via $\zeta'_R(-2)$ Hurwitz derivative) is GeoVac taxonomy. KPS values themselves are external; classification is internal.
- **A7 (Prop 4.3, F vs S decomposition):** PSLQ disjointness uses framework's own arithmetic conventions; "categorically different observables" is an interpretive label, not a theorem.
- **A8 (Thm 5.1, CHM identification):** Convergence rate inherits from Papers 38, 42, 45; all internal.
- **A15 (HS↔CC M2 inheritance):** Internal-to-Paper-18 taxonomic placement.

**MIXED chains:**
- **A1, A2, A3:** External (KPS) at the value level; M2/M3 transcendental tagging is GeoVac.
- **A9, A10:** Claim external match to KPS Table 1/2; A9 fails the match (E); A10 matches per-dim-γ only.

### Overstatement findings

| exact phrase | location | suggested honest replacement |
|---|---|---|
| "**bit-exactly** for both the conformally coupled scalar... and the free Camporesi--Higuchi Dirac" | Abstract | Keep for $S^3$ where it's true; add caveat that the *terminology* "Weyl" in d=3 corresponds to KPS's complex Dirac. |
| "the **first verification of the master Mellin engine on non-GeoVac-internal physics observables**" | Abstract, Rem 3.6 | "to our knowledge, the first verification..."; explicitly note that the M1/M2/M3 taxonomy is itself the GeoVac framework, only the **target value** (KPS partition function) is external. |
| "the bit-exact match of Sec. 3 **extends to round** $S^{5}$" | §7 line 709 | **MUST FIX** — does not extend. $F_s^{S^5}$ is **4× the KPS Table 1 value**. Either (i) the multiplicity formula $(2n+4)(n+1)(n+2)(n+3)/6$ has a $/24$ typo, in which case the value needs recomputation, or (ii) Paper 50's value is for a non-standard normalization (e.g., 4 real components or different Laplacian normalization) which must be stated. |
| "is **provably outside** the spectral-zeta ring" (entropy constant) | §4 | "is outside the M2 ⊕ M3 ring at PSLQ coefficient ceiling $10^8$" — PSLQ provides empirical search, not algebraic proof. |
| "the slope to be **exactly** $2$" | §4 | technically correct in the $n_{\max} \to \infty$ limit; one-line of justification (Cesàro on the multiplicity polynomial) would help. |
| HS↔CC "not a numerical coincidence but a **structural inheritance**" | §8.4 | "...a shared architectural property (both reduce to heat-kernel Mellin moments at operator order $k=2$)". |

## Pass B — Citation and novelty

### Citation table

| key | claimed as | verdict | what I found |
|---|---|---|---|
| `klebanov_pufu_safdi2011` | F-theorem values for scalar + Dirac on $S^3$ | **CITE-OK** | arXiv:1105.4598 confirmed; Eq. (41) and Eq. (48) match Paper 50's $S^3$ values bit-exactly. |
| `jafferis_klebanov_pufu_safdi2011` | "Towards the F-theorem: $N=2$ on $S^3$" | **CITE-OK** | arXiv:1103.1181 confirmed. |
| `beccaria_tseytlin2017` (arXiv:1702.02325) | "$\mathcal{N}=4$ conformal supergravity: the complete one-loop dilatation operator", JHEP 04 (2017) 100 | **CITE-MISATTRIBUTED — HIGH** | arXiv:1702.02325 is "$2^\infty$-Selmer groups, $2^\infty$-class groups, and Goldfeld's conjecture" by Alexander Smith (number theory). The correct Beccaria–Tseytlin paper relevant to free CFTs on $S^5$ is likely arXiv:1707.02456 ("$C_T$ for conformal higher spin fields from partition function on conically deformed sphere") or arXiv:1411.3585. **Fix the arXiv ID before broadcast.** This citation is load-bearing for Thms 7.1, 7.2, the §8.5 catalogue table, and §8.6 holographic-central-charge / Maxwell entries. |
| `hartman_kruthoff_shaghoulian_tajdini2019` (arXiv:1902.10893) | "Holography at finite cutoff with a $T^2$ deformation", JHEP 03 (2019) 004 | **CITE-MISATTRIBUTED — HIGH** | arXiv:1902.10893 is "Sphere partition functions and cut-off AdS" by Caputa, Datta, Shyam. The correct Hartman–Kruthoff–Shaghoulian–Tajdini arXiv ID is **1807.11401**. Title and JHEP reference are otherwise correct. |
| `henningson_skenderis1998` (arXiv:hep-th/9806087) | "The holographic Weyl anomaly", JHEP 07 (1998) 023 | **CITE-OK** | confirmed. |
| `chamseddine_connes1997` | "The spectral action principle", CMP 186 (1997) 731 | **CITE-OK** | standard reference, confirmed elsewhere in corpus. |
| `casini_huerta_myers2011` (arXiv:1102.0440) | hemisphere modular Hamiltonian formula $K = 2\pi\int \xi T^{00}$ | **CITE-OK** (title confirmed; specific eq not verified) | Paper confirmed as "Towards a derivation of holographic entanglement entropy". The specific CHM hemisphere formula is in §3 of that paper; I did not verify the exact equation. |
| `bisognano_wichmann1976` | duality condition modular flow | **CITE-OK** (standard external reference) | not web-verified; widely used. |
| `ryu_takayanagi2006` (hep-th/0603001) | minimum-surface entanglement | **CITE-OK** | standard. |
| `jafferis_lewkowycz_maldacena_suh2016` (1512.06431) | JLMS bulk-modular reconstruction | **CITE-OK** | standard. |
| `pastawski_yoshida_harlow_preskill2015` (1503.06237) | HaPPY perfect-tensor codes | **CITE-OK** | standard. |
| `maldacena2002` (astro-ph/0210603) | "Non-Gaussian features..." used for dS/CFT Euclidean sharpening | **CITE-OK** (title correct, use is one application of the paper, not its main subject — defensible) | confirmed. The paper does discuss dS/CFT $\Psi = Z_{\rm CFT}$ in passing. |
| `strominger2001` (hep-th/0106113) | dS/CFT proposal | **CITE-OK** | confirmed. |
| `cardy1988` | "Operator content of 2D CFT", NPB 270 (1986) 186 | **CITE-WRONG-METADATA — LOW** | Year in key (1988) vs cited year (1986) inconsistent; the actual NPB 270 paper is 1986. Cosmetic, fixable by renaming key. |
| `lei_van_leuven2024` (2406.01567) | "Modularity in $d>2$ free CFT" | **CITE-OK** (almost) | Title confirmed ("Modularity in $d>2$ free conformal field theory"), authors Yang Lei + Sam van Leuven. Paper 50 key `lei_van_leuven2024` has "van Leuven" — surname correct; spelling consistent in bibitem. **CITE-OK**. |
| `connes_vs2021` (arXiv:2004.14115) | Connes–van Suijlekom truncations | **CITE-OK** | standard CMP reference. |
| `paper7`, `paper18`, `paper22`, `paper23`, `paper28`, `paper29`, `paper32`, `paper35`, `paper38`, `paper40`, `paper42`–`paper49` | internal GeoVac preprints | **GEOVAC-INTERNAL** | not externally verifiable; reader should be aware that these are not yet peer-reviewed. |
| `debug:bh_phase0`, `debug:sprint_td_track5`, `debug:external_input_three_class_partition` | internal memos | **GEOVAC-INTERNAL** | unverifiable from outside corpus. |
| `fock1935`, `camporesi_higuchi1996` | foundational refs | **CITE-OK** (assumed; widely used in corpus) | not re-verified. |

### Problems found

1. **`beccaria_tseytlin2017` arXiv:1702.02325 — CITE-MISATTRIBUTED (HIGH).** The ID points to a number theory paper by Alexander Smith. Correct ID is probably arXiv:1707.02456 (need PI to confirm which Beccaria–Tseytlin paper was intended). Load-bearing: Thm 7.1 invokes this for the $S^5$ scalar continuum reference, and the §8.5 catalogue lists it for Maxwell on $S^3$/$S^5$ and higher-rank tensors.

2. **`hartman_kruthoff_shaghoulian_tajdini2019` arXiv:1902.10893 — CITE-MISATTRIBUTED (HIGH).** The correct ID is arXiv:1807.11401. The metadata (authors, title, JHEP reference) is otherwise consistent. Cited only in §8.4 (TT-bar discussion); not load-bearing for the main theorems.

3. **`cardy1988` year-key inconsistency — LOW.** Key says 1988, bibitem says 1986. Cosmetic.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|---|---|---|---|---|
| "the first verification of the master Mellin engine on non-GeoVac-internal physics observables" | Abstract, Rem 3.6 | Not deeply searched (claim is about internal framework's external first); strict novelty unverifiable | N/A | Soften to "to our knowledge, the first verification of the master Mellin engine on independently-published continuum CFT data". |
| "The present paper is the **first connecting that math.OA arc to a physics literature**" | §8.7 | — | N/A | Soften to "to our knowledge, the first connecting the GeoVac math.OA arc to the CFT-on-sphere partition function literature". Honest-ceiling: I cannot search "first to connect" claims. |
| "twelfth math.OA standalone in the GeoVac series" | Abstract, §8.7 | internal accounting | N/A | Verifiable from CLAUDE.md §6; not subject to external search. |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| §7.1 $F_s^{S^5}$ value is **4× the KPS Table 1 value** (multiplicity formula $/6$ should be $/24$); §7 abstract claim "bit-exact match extends to $S^5$" is false as stated | A | E (math error) | **HIGH** |
| `beccaria_tseytlin2017` arXiv:1702.02325 points to a number theory paper | B | CITE-MISATTRIBUTED | **HIGH** |
| `hartman_kruthoff_shaghoulian_tajdini2019` arXiv:1902.10893 should be 1807.11401 | B | CITE-MISATTRIBUTED | **HIGH** |
| §7.2 $D'_{\rm Dirac}(0)^{S^5}$ value matches KPS Table 2 **per-dim-γ entry**, not the total; "Weyl/single chirality" framing in d=3 is non-standard | A | C (overstatement) | **MEDIUM** |
| "first verification of master Mellin engine on non-GeoVac-internal data" — honest-ceiling novelty | A | C | **MEDIUM** |
| "provably outside the spectral-zeta ring" — PSLQ ≠ proof | A | C | **MEDIUM** |
| "the bit-exact match... extends to $S^5$" header sentence | A | C (also flagged in HIGH) | **MEDIUM** |
| HS ↔ CC "structural inheritance not numerical coincidence" — overstated | A | C | **MEDIUM** |
| §1 introduction quotes propinquity rate "$\sim (4/\pi)\log(n_{\max})/n_{\max}$" — internal-only inheritance from Paper 38 | A | B | **LOW** |
| `cardy1988` year-key vs bibitem-year mismatch (1988 key vs 1986 NPB) | B | CITE-WRONG-METADATA | **LOW** |

## Broadcast readiness: **RED**

Paper 50's $S^3$ content (Sec 3, 4, 5) is **GREEN** — the KPS bit-exact match is verified, the dual-basis projection is a correct algebraic identity, and the entropy asymptotic is reproducible. But §7 (the $S^5$ extension) contains a **load-bearing math error**: the stated scalar value $F_s^{S^5} \approx -0.02297$ is exactly **4× the KPS Table 1 d=5 entry** ($-5.74\times10^{-3}$), traceable to a multiplicity formula $(2n+4)(n+1)(n+2)(n+3)/6$ that should be $/24$ (KPS Eq. 37). The "bit-exact match extends to $S^5$" claim at line 709 is therefore false as written. The Dirac $S^5$ value matches the per-dim-γ entry of KPS Table 2 but not the total $F_D$, which deserves explicit reconciliation. Two arXiv IDs are misattributed (`beccaria_tseytlin2017` → number-theory paper; `hartman_kruthoff_shaghoulian_tajdini2019` → wrong author group). The $S^3$ side could broadcast green by itself; §7 needs either a multiplicity-correction sprint that re-derives the $S^5$ values cleanly, or §7 should be reframed as "framework's discrete machinery evaluated on a specific normalization convention" with explicit acknowledgement that the absolute value is 4× the KPS Table 1 convention. The novelty hedge ("first verification on non-GeoVac-internal data") needs the standard "to our knowledge" softening per the CONFIDENCE_REVIEW honest-ceiling rule.

## What I could NOT verify (hand to a human expert)

- Whether the "Weyl Camporesi–Higuchi Dirac" framing is consistent across the GeoVac corpus (Paper 22 §spinor sector, Paper 23, Paper 28) at the multiplicity-counting level. The d=3 "Weyl" label is unusual; d=5 has Weyl spinors with chirality but the per-dim-γ vs total ambiguity should be settled corpus-wide.
- Whether the $S^7$ negative result (Rem 7.5, 15+ basis-attempt PSLQ failure) is reproducible at the stated tolerances; I did not re-run.
- Whether the Maxwell/higher-rank-tensor catalogue entries (§8.5) carry through under the corrected $S^5$ scalar normalization.
- Whether the propagated $F_s^{S^5}$ error invalidates Thm 7.4 (log-3 cancellation): structurally Thm 7.4 is about the polynomial $(u^2 - 9/4)(u^2 - 1/4)$ coefficient and is independent of the overall multiplicity prefactor, so it should survive — but worth a clean re-derivation.
- Which precise Beccaria–Tseytlin paper was intended; CLAUDE.md §6 references should have a verifiable arXiv ID.
