# Confidence Review: Paper 44 — Operator-system extension of the BBB Krein spectral triple at finite cutoff: propagation number and envelope dependence on the Camporesi-Higuchi spinor bundle

## Calibration check

Not a calibration run. Wave 2 Lorentzian-arc audit; verdicts grounded in (a) external references (CITE checks) and (b) cross-corpus reads of Papers 32 / 45 / 43.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | $\Op_{\nmax,\Nt}^L$ is *-closed and contains identity (Prop 3.2 = `prop:axioms`) | §3, Prop. axioms | B | GEOVAC-ONLY (formal derivation from defs in §2.4-§2.5) | Direct symbolic check from selection rules and tensor-product structure; mechanical. |
| 2 | Riemannian-limit recovery at $\Nt=1$ bit-exact (Thm 4.1) | §4, Thm `thm:riemannian_limit` | B | GEOVAC-ONLY (rests on test file `test_operator_system_lorentzian.py`) | The proof argument that $M_0^{\temp}=1$ ⇒ tensor with identity is the matrix verbatim is correct symbolically; the bit-exactness is a code-vs-code statement, not external validation. |
| 3 | $\prop_{\achievable}(\Op^L)=2$ for $\nmax\ge 2$, $\Nt\ge 1$ (Thm 5.3 eq. 5.7) | §5, eq. `eq:prop_ach` | B | GEOVAC-ONLY (rests on Paper 32 §III prop=2, which is GeoVac-internal) | Argument is structurally correct: chirality-doubling $\oplus$ commutative temporal preserves prop=2 from Weyl sector. Paper 32 §III itself rests on its own test suite, not external published prop=2 on $S^3$. |
| 4 | $\prop_{\mathrm{full}}(\Op^L)=\infty$ at every finite cell | §5, eq. `eq:prop_full_inf` | A | EXTERNAL (basic linear algebra) | Proof in §5.3 is correct: $(\Op^L)^k \subset \mathcal{V}_\achievable$ for every $k$ by chirality-block + commutative-temporal closure; $\mathcal{V}_\achievable \subsetneq \Bcal(\Krein)$ strictly. Standard. |
| 5 | Krein-positive restriction is trivial: $\Op^{L,+} = \Op^L$ (Prop 6.2) | §6, `prop:krein_positivity_trivial` | A | EXTERNAL (direct 2×2 block computation) | The 2×2 matrix calculation $\gamma^0 (M\oplus M) = (M\oplus M)\gamma^0$ is correct and trivial. |
| 6 | Witness pair Frobenius residual: 0.3812 at $\nmax=2$, 0.3589 at $\nmax=3$ (Thm 7.1) | §7, Thm `thm:witness_pair` | B | GEOVAC-ONLY (JSON file `debug/data/l3a_1_lorentzian_operator_system.json`) | Magnitude is consistent with the Paper 32 §III Weyl-sector ~14.9% residual scaled by chirality-doubling (2× to ~30%, then further inflated by basis filling). No independent reference. |
| 7 | Connes-vS Toeplitz $\prop(C(S^1)^{(n)})=2$, cited as "verbatim" anchor (Prop 4.2) | Abstract; §1 lines 182-183; §5.2 `prop:paper32_prop` proof | A | EXTERNAL (Connes-vS arXiv:2004.14115) | The Toeplitz prop=2 result is genuinely in Connes-vS §4 (per web verification of paper structure and "Spectral Truncations" abstract). The specific proposition number ("Proposition 4.2") is not verified verbatim but the structural claim is sound. |
| 8 | BBB classification places $(3,1)$ West-coast at $(m,n)=(4,6)$ with primitive signs (5.16) | §2.3, eq. `eq:bbb_signs` | A | EXTERNAL (BBB 2018) | BBB 2018 (J. Math. Phys. 59, 062303; arXiv:1611.07062) verified by web search; the $(m,n)$ assignment is the central content of that paper. |
| 9 | van den Dungen Prop 4.1 lift gives $\DL = i(\gamma^0\otimes\partial_t + D_{GV}\otimes I)$ | §2.2, eq. `eq:DL_def` | B | MIXED (vdD framework external; specific GeoVac instantiation internal) | van den Dungen 2016 paper exists and is the right framework; the specific application to $S^3$ truncations is GeoVac-internal. |
| 10 | "No published Lorentzian propinquity construction exists as of May 2026" (abstract) | Abstract, lines 138-143 | D | Cannot be confirmed by search | This is a priority claim — see Pass B novelty section. Web search shows Latrémolière 2017/2022, Hekkelman-McDonald 2024, Farsi-Latrémolière 2024/2025 are strictly Riemannian; no Lorentzian propinquity found, but this is "honest absence of evidence" only. |
| 11 | Higher-rank Lie group extension (§9.5): "prop=2 should hold at every $G$" by Paper 40 L1' lemma | §9.5, lines 1378-1383 | D | MIXED, somewhat speculative | Conjectural; presented honestly as a "would be a natural ~4 week extension." Acceptable scope-honest framing. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value | Survives? |
|---|---|---|---|---|
| Dim of $\HGV^{\nmax}$ at $\nmax=2$ | $16 = (2/3)\cdot 2 \cdot 3 \cdot 4$ | direct formula | $16$ ✓ | yes |
| Dim of $\HGV^{\nmax}$ at $\nmax=3$ | $40 = (2/3)\cdot 3\cdot 4\cdot 5$ | direct formula | $40$ ✓ | yes |
| Dim of $\HGV^{\nmax}$ at $\nmax=5$ | $140 = (2/3)\cdot 5\cdot 6\cdot 7$ | direct formula | $140$ ✓ | yes |
| Achievable envelope dim at $(\nmax,\Nt)=(2,1)$ | $64 = 8^2\cdot 1$ | $\dim_\Weyl=8$ at $\nmax=2$ → $\dim_\Weyl^2=64$ | $64$ ✓ | yes |
| Achievable envelope dim at $(3,1)$ | $400 = 20^2$ | $\dim_\Weyl=20$ at $\nmax=3$ → $20^2=400$ | $400$ ✓ | yes |
| $\dim(\Op^L)$ at $(2,5)$ | $70$ = $14\cdot 5$ | $14$ spatial × $\Nt=5$ | $70$ ✓ | yes |
| Wedge dim at $(\nmax,\Nt)=(2,3)$ | $16$ | $(16/2)\cdot\lceil 3/2\rceil = 8\cdot 2 = 16$ ✓ | $16$ ✓ | yes |

All arithmetic survives.

### Circularity map (GEOVAC-ONLY chains)

- **Claim 2 (bit-exact Riem limit at $\Nt=1$)** rests on Paper 32 §III "FullDiracTruncatedOperatorSystem" construction + paper's own test file. There is no published external benchmark for "round-$S^3$ chirality-doubled truncated operator system" against which to validate. Acceptable as "internal-consistency falsifier," but not external truth.
- **Claim 3 ($\prop_{\achievable}=2$)** rests on Proposition `prop:paper32_prop` (Paper 32 §III), itself an internal computation. The chirality-doubling and temporal-tensoring arguments lifting prop=2 are correct algebraically, but the base case is GeoVac-internal.
- **Claim 6 (witness pair 0.3812 / 0.3589)** rests on Paper 32 §III Weyl-sector residual ~14.9% (internal), which is then scaled by chirality-doubling. No external benchmark.
- **Claim 1 (operator system axioms)** rests on the Wigner relation `eq:Y_conj` for complex hyperspherical harmonics, which is external (standard) — but the lift to the spinor bundle uses GeoVac's Clebsch-Gordan + 3-Y machinery (internal).

The strongest external anchor is **Claim 4 ($\prop_{\mathrm{full}}=\infty$)**: this is a clean linear-algebra ceiling argument independent of GeoVac internals.

### Overstatement findings

1. **Abstract line 116-117**: "We prove $\Op_{\nmax, \Nt}^L$ is a Lorentzian truncated operator system."
   - Body shows: it is a $*$-closed, unital, non-multiplicatively-closed linear subspace. This **is** the standard definition of an operator system (per Paulsen 2002, cited correctly).
   - Verdict: **OK**, the body supports the abstract; this is conventional terminology.

2. **Abstract lines 122-124**: "we obtain $\prop(\Op_{\nmax,\Nt}^L) = 2$ for every $\nmax \ge 2$, matching Paper 32 §III verbatim."
   - Theorem 5.3 actually proves only $\nmax \in \{1, 2, 3\}$ with $\Nt \in \{1, 3, 5\}$ numerically and a structural argument for general $(\nmax, \Nt)$ via chirality-doubling preservation.
   - The "for every $\nmax \ge 2$" claim is supported by the structural argument, not numerically.
   - Verdict: **B/borderline-C, LOW**. Suggest soft tightening: "we prove $\prop_{\achievable}=2$ structurally for every $\nmax \ge 2$ (numerically verified at $\nmax \in \{2, 3\}$)."

3. **Abstract line 138-143**: "No published Lorentzian propinquity construction exists as of May 2026: Latrémolière 2017/2023 and Hekkelman-McDonald 2024 are strictly Riemannian."
   - This is a priority/novelty claim. Web search corroborates "strictly Riemannian" for the cited works, but absence-of-evidence is not proof.
   - Verdict: **C, LOW**. Suggest "to our knowledge, no published Lorentzian propinquity construction exists." Paper does already hedge well in §1.2(b) and §9.1 ("To our knowledge, no construction..."), so abstract should match.

4. **§1.1 item 3 (line 245-249)**: "verbatim with Paper 32 §III" — accurate, since the Weyl-sector prop=2 transports verbatim under chirality-block doubling. **OK.**

5. **§5.4 / `rem:cvs_sharpening` (lines 934-991)**: The "envelope-relative informative content" framing is conceptual but well-hedged. **OK.**

Overall: paper is honest about scope. The "honest scope" paragraph in conclusion (lines 1447-1462) explicitly names continuum / GH-limit / propinquity convergence / Krein-positive state-space convergence / non-commutative-temporal multipliers / gamma-matrix multipliers as all open. **This is exemplary scope discipline.**

## Pass B — Citation and novelty

### Citation table

| `\cite` key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `avery_wen_avery2002` | "Avery, Wen, Avery (2002), J. Math. Chem. 32, 65-78" — 3-Y integrals on the 3-sphere via Gegenbauer recurrence (§2.4, line 498) | **CITE-CANT-FIND** | Could not locate this specific paper. The genuine Avery/Wen collaboration is Wen & Avery, "Some properties of hyperspherical harmonics," J. Math. Phys. 26 (1985) 396-403. The 2002 J. Math. Chem. issue exists but the article at vol 32 pages 65-78 by these authors could not be confirmed via web search. **This may be a fabricated citation parallel to the verified-fake `avery_wen_avery1986` flagged in CLAUDE.md §3 (Track v3.19.0).** Papers 38, 39 use a different fake variant; Paper 44 uses a third unverified variant. |
| `connes_vs2021` | "Connes-vS, Comm. Math. Phys. 383 (2021) 2021-2067; arXiv:2004.14115" — operator system + propagation number, Toeplitz prop=2 | **CITE-OK** | Verified. arXiv:2004.14115. Toeplitz prop=2 claim supported (Section 4 "The Fejér-Riesz operator system"). Proposition number "4.2" stated in Paper 44 not bit-verified but structurally correct. |
| `bizi_brouder_besnard2018` | "BBB, J. Math. Phys. 59 (2018) 062303; arXiv:1611.07062" — $(m,n)$ classification of indefinite spectral triples | **CITE-OK** | Verified. |
| `camporesi_higuchi1996` | "Camporesi-Higuchi, J. Geom. Phys. 20 (1996) 1-18" — Dirac on spheres | **CITE-OK** | Verified. arXiv:gr-qc/9505009. |
| `vandungen2016` | "van den Dungen, Math. Phys. Anal. Geom. 19 (2016) Art. 4; arXiv:1505.01939" — Krein spectral triples, Prop 4.1 lift | **CITE-OK** | Verified. |
| `latremoliere_metric_st_2017` | "Latrémolière, Adv. Math. 404 (2022), 108393, 56pp; preprint arXiv:1811.10843 (2017)" | **CITE-OK** | Verified per dispatcher context + own check. |
| `latremoliere2018` | "Latrémolière, Trans. AMS 368 (2016)" — Gromov-Hausdorff propinquity | **CITE-OK** | Per dispatcher context (Wave 1 fix). |
| `hekkelman2022` | "Hekkelman, Truncated Geometry on the Circle, arXiv:2111.13865, 2022" | **CITE-OK** | Per dispatcher context (Wave 1 fix). |
| `hekkelman_mcdonald2024` | "Hekkelman-McDonald, J. Noncommut. Geom., to appear (2024); arXiv:2403.18619" — spectral truncations of $T^d$ | **CITE-OK** (likely; standard load-bearing-Latrémolière-companion reference) | Not directly re-verified in this session but consistent with corpus. |
| `hekkelman_mcdonald2024b` | "Hekkelman-McDonald, J. Funct. Anal., to appear (2025); arXiv:2412.00628" | **CITE-OK** (likely) | Not re-verified. |
| `leimbach_vs2024` | "M. Leimbach and W.D. van Suijlekom, Adv. Math. 439 (2024), 109496" — GH convergence of torus spectral truncations | **CITE-OK** | Per dispatcher context (Wave 1 fix: M. Leimbach, not L. Leimbach). |
| `bykov_minguzzi_suhr2024` | "Bykov-Minguzzi-Suhr, arXiv:2412.04311, 2024" — Lorentzian metric spaces and GH convergence | **Plausible**, not re-verified | |
| `mondino_samann2025` | "Mondino-Sämann, arXiv:2504.10380, 2025" — synthetic Lorentzian GH | **Plausible**, not re-verified | Consistent with Paper 48 usage. |
| `nieuviarts2025a` | "Nieuviarts, arXiv:2502.18105 (v3, May 2025)" | **CITE-OK** | Consistent with CLAUDE.md §1 entries on Nieuviarts arXiv:2502.18105. |
| `nieuviarts2025b_proceedings` | "Nieuviarts, arXiv:2512.15450 (v2, May 2026)" — proceedings synthesis | **CITE-OK** | Consistent with Paper 43 §3.2 footnote and CLAUDE.md Paper 43 entry. |
| `paulsen2002` | "Paulsen, Cambridge Studies in Advanced Mathematics 78, 2002" | **CITE-OK** | Standard reference for operator systems. |
| `franco_eckstein2014` | "Franco-Eckstein, Class. Quantum Grav. 30 (2013), 135007; arXiv:1212.5171" | **CITE-OK** | Standard. Note: year tag in bibitem says (2014) but vol 30 (2013) is listed — venue-year mismatch is conventional (online 2013 / vol 30 / cite-year 2014). LOW. |
| `devastato_lizzi_martinetti2018` | "DLM, J. Math. Phys. 59 (2018), 092304" | **Plausible** | Not re-verified, but consistent. |
| `strohmaier2006` | "Strohmaier, J. Geom. Phys. 56 (2006), 175-195" | **CITE-OK** | Standard. |
| `bisognano_wichmann1976` | "Bisognano-Wichmann, J. Math. Phys. 17 (1976), 303-321" | **CITE-OK** | Classic. |
| `hartle_hawking1976` | "Hartle-Hawking, Phys. Rev. D 13 (1976), 2188-2203" | **CITE-OK** | Classic. |
| `sewell1982` | "Sewell, Ann. Phys. (NY) 141 (1982), 201-224" | **CITE-OK** | Classic. |
| `unruh1976` | "Unruh, Phys. Rev. D 14 (1976), 870-892" | **CITE-OK** | Classic. |
| `farsi_latremoliere2024` | "F-L, arXiv:2404.00240, 2024" | **CITE-OK** (likely) | Consistent with cross-corpus usage. |
| `farsi_latremoliere2025` | "F-L, arXiv:2504.11715, 2025" | **CITE-OK** | Consistent with Paper 47 bibliography update note. |
| `toyota2023` | "Toyota, arXiv:2309.13469, 2023" | **CITE-OK** (likely) | |
| `paper24, paper32, paper38, paper39, paper40_unified, paper42, paper43` | Internal | Internal cross-references | Not separately CITE-checked in Pass B (all confirmed in corpus). |

### Problems found

**CITE-CANT-FIND (HIGH-severity load-bearing):**
- **`avery_wen_avery2002`** ("Avery, Wen, Avery (2002), J. Math. Chem. 32, 65-78"): This is the 3-Y integral foundation cited at §2.4 line 498. The genuine collaboration is **Wen & Avery 1985, J. Math. Phys. 26(3), 396-403** ("Some properties of hyperspherical harmonics"). The 2002 J. Math. Chem. variant could not be located via web search; vol 32 (2002) exists but the article at pp. 65-78 by these authors was not found. This is the same pattern as the verified-fake `avery_wen_avery1986` flagged in CLAUDE.md §3 (Track v3.19.0). **Severity: HIGH** — load-bearing reference for the spinor multiplier construction.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "No published Lorentzian propinquity construction exists as of May 2026" | Abstract line 138-143 | "Lorentzian propinquity," Latrémolière + Lorentzian, Hekkelman + Lorentzian, Farsi-Latrémolière + Lorentzian | No (consistent with stated absence) | **Soften to** "to our knowledge, no published Lorentzian propinquity construction exists." The body language §1.2(b), §9.1 is already correctly hedged; align abstract. **MEDIUM**. |
| "We exhibit a witness pair lifting the Paper 32 §III witness verbatim" | Abstract; §7 | Not a novelty claim — restatement of internal construction | n/a | OK. |
| "Envelope-dependent propagation number" framing | Abstract; Thm 5.3 | "envelope-relative propagation number" + operator systems | No published prior use of "envelope-dependent" / "envelope-relative" terminology for propagation number found. Likely original phrasing for a real observation. | Acceptable. The structural distinction is real (chirality-doubling + commutative-temporal block scalar multipliers from full envelope). |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `avery_wen_avery2002` not findable / likely fabricated J. Math. Chem. citation | B | CITE-CANT-FIND | **HIGH** |
| Abstract "for every $\nmax \ge 2$" prop=2 numerically verified only at $\nmax \in \{2,3\}$ | A | C (mild overstatement) | LOW |
| Abstract "No published Lorentzian propinquity construction exists" stronger than body's "to our knowledge" | A | C | MEDIUM |
| Memo-JSON discrepancy at $\nmax=3$ (acknowledged in `rem:nmax3_residual`) | A | self-flagged | LOW (already managed) |
| Reliance on GeoVac-internal Paper 32 §III prop=2 as the structural anchor for Claim 3 | A | B (circularity acknowledged) | LOW (paper makes this explicit) |
| Franco-Eckstein bibitem year (2014) vs venue year (2013) | B | CITE-WRONG-METADATA (cosmetic) | LOW |

**HIGH: 1, MEDIUM: 1, LOW: 4.**

## Broadcast readiness: YELLOW

The paper is mathematically clean, scope-honest, and structurally well-organized. The novel content (envelope-dependent propagation, trivial Krein-positive restriction) is supported by correct linear-algebra arguments and consistent with the rest of the Lorentzian arc (Papers 42, 43, 45). All chains of inference are explicit, and the Honest Scope conclusion is exemplary. **The one blocking concern** is the `avery_wen_avery2002` citation (CITE-CANT-FIND, HIGH-severity, load-bearing for the spatial multiplier construction). This is the same fabrication-pattern that CLAUDE.md §3 records for the verified-fake `avery_wen_avery1986`. **Recommended action: replace with Wen & Avery, J. Math. Phys. 26 (1985), 396-403 — "Some properties of hyperspherical harmonics" (verified extant)** unless the 2002 paper can be substantively located. Two minor abstract softenings (MEDIUM/LOW) round out the errata pass. After these fixes the paper is GREEN.

## What I could NOT verify (hand to a human expert)

- Whether `avery_wen_avery2002` (J. Math. Chem. 32, pages 65-78) genuinely exists or is a fabricated/misattributed citation. The author should check personal records / Google Scholar profile.
- Whether `prop = 2` for the Toeplitz $C(S^1)^{(n)}$ is literally "Proposition 4.2" of Connes-vS 2021 or a different proposition number in the same Section 4 (structural claim is correct either way).
- Whether the witness-pair Frobenius residual scaling 0.3812 (chirality-doubled $\nmax=2$) is the expected lift from Paper 32 §III's ~14.9% (Weyl sector). The lift mechanism (doubling) seems right but the precise factor is not externally checkable.
- Conjecture in §9.5 that $\prop_\achievable = 2$ extends to all compact $G$ in the Paper 40 class — not verified at $\SU(3)$, $\SU(4)$, $\mathrm{Sp}(2)$, $G_2$ (paper is honest about this being an extension target, not closed).
