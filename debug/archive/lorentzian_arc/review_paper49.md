# Confidence Review: Paper 49 — Operator-system Lorentzian pre-length spaces and the strong-form Krein--Mondino--Sämann bridge

**Source:** `papers/group1_operator_algebras/paper_49_oslpls_strong_form_bridge.tex` (May 2026 draft, 2946 lines).
**Reviewer:** GeoVac Confidence Reviewer (CONFIDENCE_AUDITOR + CITATION_CHECKER unified pass), Wave 2 Lorentzian-arc audit.
**Date of review:** 2026-06-02.

---

## Calibration check

Not a calibration run; standard audit.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict (A–E) | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| A1 | The Q1′ open question of Paper 48 is closed at strict-strong-form on the full M≠0 enlarged substrate. | Abstract, §1.1, Conclusion | **B** | GEOVAC-ONLY | "Closed" is internally defined; depends on whether claims A2–A4 actually hold. |
| A2 | OSLPLS is a well-defined category containing commutative Mondino--Sämann as a faithful sub-category. | §3, Thm 3.7 (`thm:embedding_iota`) | **B** | MIXED | The construction is reasonable as stated; the faithful-sub-category claim rests on standard Gelfand-Naimark. Sketch only — proof-grade verification is not provided. |
| A3 | The bridge functor W^flip: KreinMetaMet_pp^{flip,full} → OSLPLS extends the Paper 48 W functor to the full M≠0 substrate. | §4, Def 4.4, Thm 4.10 | **B** | GEOVAC-ONLY | Construction is mechanically reasonable; verification depends on Krein-side bridge claims of Paper 48. |
| A4 | (B2′) Strict super-additivity of ℓ^{OS} on three-orbit triples via "Lindblad/Uhlmann chain inequality" (eq:lindblad_uhlmann): ΔS^{1→3} ≤ ΔS^{1→2} + ΔS^{2→3} (strict when pairwise distinct). | §5.4, Thm 5.10 (`thm:strict_super_additivity`) | **E** | MIXED | **LOAD-BEARING LOGICAL GAP.** Quantum relative entropy does NOT satisfy a triangle inequality. The proof states the chain inequality "follows from the triangle-type relation for relative entropies under composition of quantum channels" — but no such general triangle relation exists in the literature. The max-relative entropy does, but standard relative entropy does not (literature confirms this is a known fact; e.g. Wilde lecture notes 2015 explicitly say "the relative entropy does not satisfy the triangle inequality, so this simple proof does not work"). Without a triangle inequality on quantum relative entropy, the proof structure as written is incomplete. See **HIGH-severity finding F1** below. |
| A5 | The Triple-Intersection Cocycle Identity (TICI) is implicit in Connes 1973 (C4) but not previously named explicitly. | §5.3, Thm 5.7 (`thm:tici`) | **C/D** | MIXED | The cocycle composition rule (C4) is genuine Connes 1973. Whether TICI is "not previously named" is a priority claim I cannot confirm (Pass B). The mathematical content is plausible (it is a direct corollary of (C4) by composition), but the proof in the paper (§5.3) uses an unjustified rearrangement step: "Using u_{23}·u_{12} = σ_{t}^{ω_1}(u_{23})·u_{12}, or equivalently identifying via (C4)" — this rearrangement is the content of (C4) restated, not derived; it is circular in the proof as written. |
| A6 | Numerical panel: Λ joint propinquity bit-exact at (2,3), (3,5), (4,7) matching Paper 45. | §7.1, Tab 7.1 | **B** | GEOVAC-ONLY | Numbers match Paper 45 by construction; "bit-exact" is a statement of internal consistency, not external truth. The very fact that the strong-form Λ equals the K⁺-weak-form Λ from Paper 45 is described as "free upgrade" in the paper — this is internal-consistency at most, not novel verification. |
| A7 | Uhlmann monotonicity deficits substantially positive (66, 69, 81 nats) at three panel cells. | §7.2, Tab 7.2 | **B** | GEOVAC-ONLY | Verifies positivity on the SELECTED triples found by the "greedy K-commutator + operator-commutator algorithm" — i.e. these are triples chosen to make the deficit positive. This is consistent with the claim but does not verify the THEOREM (which states super-additivity holds for ALL three-orbit triples). |
| A8 | Riemannian-limit recovery at N_t=1 bit-exact. | §7.3, Tab 7.3 | **A** | EXTERNAL | Internal consistency check inherited from Paper 38/42 architecture; this is a well-defined and load-bearing falsifier. |
| A9 | Twin-paradox-as-quantum-information: relative-entropy monotonicity supplies the strict inequality of synthetic Lorentzian reverse triangle. | Abstract, §1.3, §5.5, §6.1 | **C** | MIXED | Headline claim is rhetorically strong. If A4 (the chain inequality) does not hold without further qualification, then this dual identification is also incomplete. The physical intuition is appealing but the operator-algebraic content is gated by the gap in Theorem 3.3. Suggested honest replacement: "The OSLPLS reading proposes quantum relative-entropy monotonicity as a candidate structural ingredient; whether the proposed chain inequality holds in full generality is open." |
| A10 | Q2′ (non-commutative MS) is CLOSED by the OSLPLS construction. | §9.1 | **C** | GEOVAC-ONLY | Section 9.1 reframes Q2′ as "OSLPLS *is* the non-commutative LPLS concept, built from the operator-algebraic side." This is a re-definition: the original Q2′ asked for a SYNTHETIC-side extension; the paper closes it by claiming the operator-algebraic-side OSLPLS replaces it. This is honest in §9.1 but the abstract's "Paper 49 *does not* close Q2′" is then contradicted by §9.1's "It also resolves Q2′." Two passages of the paper are inconsistent. |
| A11 | Bulk-gravity dual via Bousso--Casini--Fisher--Maldacena 2020 kink transform. | §8 | **E** (citation) → see B-pass | MIXED | The connection to BCFM/kink-transform is conjectural and clearly labeled as such; **but the citation has wrong authors** (see Pass B, **HIGH-severity finding F2**). |
| A12 | "First explicit identification of quantum relative-entropy monotonicity as load-bearing for the reverse triangle inequality of synthetic Lorentzian geometry." | Abstract, §1.3, §6.1, §6.2 | **D** | UNVERIFIABLE | Cannot confirm or refute novelty claim absent a domain expert; downgrade to "to our knowledge, the first explicit" (paper does include this hedge in places but not all). |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value/error | Survives? |
|-------|---------------|----------------------|--------------------------|-----------|
| Λ(2,3) joint propinquity | 2.0745510936998897 | Paper 45 reference panel (GeoVac-internal) | N/A — this is internal cross-check only | INTERNAL CONSISTENCY |
| Λ(3,5) | 1.6100599680657361 | Paper 45 | N/A | INTERNAL CONSISTENCY |
| Λ(4,7) | 1.3223327942828407 | Paper 45 | N/A | INTERNAL CONSISTENCY |
| Uhlmann deficit (2,3) | 66.998 nats | None (paper's own panel) | N/A — selected triple | INTERNAL ONLY |
| Riemannian-limit residual N_t=1 | 0.0 (float64) | Paper 38/42 architecture | N/A | INTERNAL CONSISTENCY |
| BCFM kink rapidity | 2π·s | arXiv:2007.00230 | The arXiv page confirms a kink transform of bulk initial data about the RT surface; precise rapidity-parameter value 2πs is consistent with the canonical β=2π modular period but I did not see the explicit "2πs" identification in the abstract metadata. Plausible. | LIKELY OK structurally |

**Note:** No external benchmark exists for any of the load-bearing numerical claims (they are all internal architectural consistency checks against earlier GeoVac papers). The claim "Λ matches Paper 45 reference exactly" is mathematically tautological under the "free upgrade" reading because the strong-form Λ is constructed to equal the weak-form Λ by Paper 46 Theorem 3.5 — the bit-exact agreement is a definitional consequence, not a verification.

### Circularity map

**The full claim chain bottoms out GEOVAC-ONLY at the load-bearing entry points:**

1. **Λ joint propinquity** → Paper 46 (strong-form) → Paper 45 (K⁺-weak-form) → Paper 38 (Riemannian SU(2) propinquity). All GeoVac papers. No external reproduction.
2. **Uhlmann panel** → Paper 49 internal driver `q1prime_phase2b3_panel_compute.py` selecting triples on three modular orbits → the construction depends on Paper 42 BW K_α = J_polar integer spectrum. GeoVac-internal.
3. **Strict super-additivity proof** → asserts chain inequality on relative entropies → claims this follows from Lindblad/Uhlmann monotonicity, but the actual chain inequality $S(\rho_1\|\rho_3) \le S(\rho_1\|\rho_2) + S(\rho_2\|\rho_3)$ is **NOT a theorem in the literature** (quantum relative entropy does not obey a triangle inequality). This is the **load-bearing gap**.
4. **TICI** → claimed as a corollary of Connes (C4); the proof in §5.3 is essentially restated (C4) rather than derived from it.
5. **OSLPLS axiom transport** → Theorem 3.1 walks through MS axioms 1-8. **Substantive observation in §3.1**: "axiom (d) [reverse triangle] is the only MS axiom whose proof uses commutativity structurally." This is honest about the proof difficulty.

**Externally verified ingredients:** Connes 1973 (cocycle exists, intertwines, satisfies (C4)) — verified at level of citation; Tomita–Takesaki modular theory — standard; Connes–Rovelli 1994 thermal-time hypothesis at single KMS state — confirmed.

### Overstatement findings

| Exact phrase | Suggested honest replacement |
|--------------|------------------------------|
| Abstract: "established via the Connes–Rovelli thermal-time stack across distinct KMS states, **a construction not previously published**" | Soften to: "to our knowledge, not previously isolated as a stack-consistency property explicitly in the published literature." |
| Abstract / §1.3 / §5.5: "**To the author's knowledge this is the first explicit identification of an information-theoretic underpinning for the reverse triangle inequality of synthetic Lorentzian geometry.**" | OK as stated (already hedged with "to the author's knowledge"); but note that the claim's substance is gated by the Theorem 3.3 gap. |
| §9.1 title: "**Q2′: Non-commutative MS extension --- CLOSED by the OSLPLS construction**" | This contradicts the Abstract claim "Paper~49 does *not* close Q2′" and the Honest Scope statement. Recommend: either retract §9.1 closure framing OR retract abstract's "does not close" framing. The two cannot both stand. |
| §6.1: "**The OSLPLS reading reveals that the same inequality at the operator-algebraic level rests on a fundamentally information-theoretic inequality (the data-processing inequality of Lindblad).**" | Conditional: this rests on Theorem 3.3 being valid; if the chain-inequality gap is closed (e.g., by replacing the general statement with a structural property of the specific cocycle entropy production on this substrate), the statement stands. As is, soften to "is conjectured to rest on..." or fix the proof. |
| Abstract: "the BCFM 2020 paper" attribution incorrectly lists Bousso, Casini, Fisher, Maldacena. | This is a Pass-B citation error, not a Pass-A overstatement. Fix at bibitem level. |
| Abstract: "(B3′) pre-compactness via the Paper 44 propagation-number bound prop ≤ 4" | The propagation = 2 was confirmed only at (2,3); at (3,5) and (4,7) the compute is "blocked by memory limitations." This is documented honestly in §7.4. No edit needed. |

---

## Pass B — Citation and novelty

### Citation table (selected load-bearing entries)

| \cite key | Claimed as | Verdict | What I found |
|-----------|------------|---------|--------------|
| `connes1973` | Source of cocycle Radon-Nikodym + (C4) composition rule. | **CITE-OK** | Connes, "Une classification des facteurs de type III," Ann. Sci. ENS (4) 6 (1973) 133–252. Confirmed via numdam.org/item/ASENS_1973_4_6_2_133_0/. DOI 10.24033/asens.1247. Page range 133–252 correct. Content (cocycle + Cor 3.2 gauge uniqueness) confirmed. |
| `connes_rovelli1994` | Thermal-time hypothesis at single KMS state. | **CITE-OK** (minor) | A. Connes, C. Rovelli, Class. Quantum Grav. 11 (1994) 2899–2917, also arXiv:gr-qc/9406019. ADS confirms (1994CQGra..11.2899C). Minor: the original journal title uses "General Covariant" not "Generally Covariant"; the GeoVac citation uses "generally" — cosmetic, no severity. Page range 2899–2917 vs 2899–2918 in ADS (off-by-one, cosmetic). |
| `bratteli_robinson1981` | KMS-state structure under modular automorphism groups, Thm 5.3.10. | **CITE-OK** | Springer 1981 edition confirmed; book exists with that title and authors. Theorem 5.3.10 reference cannot be verified without consulting the book physically but the cited section/topic is consistent with the book's known content. |
| `uhlmann1977` | Relative entropy and Wigner-Yanase-Dyson-Lieb concavity; "data-processing inequality" usage. | **CITE-OK** (modulo proof concern) | Uhlmann, Comm. Math. Phys. 54 (1977) 21–32. Confirmed via projecteuclid.org. The paper does prove monotonicity of relative entropy under identity-preserving completely positive maps. **But the GeoVac paper invokes a "chain inequality" not present in Uhlmann's paper** — Uhlmann gives monotonicity, not a triangle/chain inequality. **CITE-DOESNT-SUPPORT** for the specific use in §5.4 Theorem 3.3 proof. |
| `lindblad1975` | "The data-processing inequality used in Theorem 3.3." | **CITE-DOESNT-SUPPORT** | Lindblad, Comm. Math. Phys. 40 (1975) 147–151. Confirmed: the paper proves monotonicity of relative entropy under trace-preserving CP maps. **However, this is NOT a triangle/chain inequality on three states.** The GeoVac paper's Theorem 3.3 proof invokes "Lindblad/Uhlmann chain inequality" — neither Lindblad nor Uhlmann establishes the chain inequality $S(\rho_1\|\rho_3) \le S(\rho_1\|\rho_2) + S(\rho_2\|\rho_3)$. The data-processing inequality and a triangle inequality are different statements. |
| `bousso_casini_fisher_maldacena2020` | "Gravity dual of Connes cocycle flow," Phys. Rev. D 102 (2020) 066008, arXiv:2007.00230. **GeoVac bibitem lists "R. Bousso, H. Casini, Z. Fisher, J. Maldacena".** | **CITE-WRONG-METADATA → near-MISATTRIBUTED** | arXiv:2007.00230's actual authors are **Raphael Bousso, Venkatesa Chandrasekaran, Pratik Rath, Arvin Shahbazi-Moghaddam**. The cited "Casini, Fisher, Maldacena" are NOT the authors of this paper. The CITE key `bousso_casini_fisher_maldacena2020` and three of the four bibitem author names are wrong. The arXiv ID, title, journal, volume, and page are all correct, so the paper IS the intended reference, but the author list is fabricated. This is exactly the **Fursaev/Solodukhin failure mode** flagged in CLAUDE.md §3. |
| `mondino_samann2025_pointed` | Title "Pointed Lorentzian Gromov-Hausdorff convergence of covered pre-length spaces," arXiv:2504.10380~v4 (December 2025). | **CITE-WRONG-METADATA** | The actual title of arXiv:2504.10380 (including v4) is **"Lorentzian Gromov-Hausdorff convergence and pre-compactness."** The GeoVac bibitem title is a paraphrase/description, not the actual paper title. The paper does discuss pointed LGH and covered LPLS (per Paper 48 §5.1 cross-reference), but the bibitem's "title" is not the published title. Same error propagates from Paper 48 (also CITE-WRONG-METADATA). |
| `connes_vs2021` | Spectral truncations in NCG and operator systems. Comm. Math. Phys. 383 (2021) 2021–2067, arXiv:2004.14115. | **CITE-OK** | Standard reference, content matches the use (UCP morphisms, propagation number, operator systems). |
| `kunzinger_samann2018` | Lorentzian length spaces, Ann. Global Anal. Geom. 54 (2018) 399–447, arXiv:1711.08990. | **CITE-OK** | Standard reference, content matches. |
| `latremoliere2018` | Trans. AMS 368 (2016) 365–411, arXiv:1302.4058. | **CITE-OK** (per user note) | User confirmed bibitem fix applied. arXiv:1302.4058 confirmed exists. Note: journal year 2016 not 2018 — but the user-supplied note says fix already applied. |
| `che_perales_sormani2025` | DGA 103 (2026), arXiv:2510.13069. | **CITE-OK** (per user note) | User confirmed fix applied. |
| `minguzzi_suhr2024` | arXiv:2209.14384. | **CITE-OK** (per user note) | User confirmed fix applied. |
| `sakovich_sormani2024` | arXiv:2410.16800. | **CITE-OK** (per user note) | User confirmed fix applied. |
| `wilde2017` | Modern textbook on data-processing. | **CITE-OK** | Standard text; supports the framing of Lindblad/Uhlmann DPI but NOT the specific chain inequality the GeoVac paper invokes. |
| `paper48` (Paper 48) | K⁺-weak-form Krein-MS bridge. | **CITE-OK** | Internal cross-reference; Paper 48 exists in the corpus and §1.4 explicitly opens Q1′. |
| `paper44`, `paper45`, `paper46`, `paper47`, `paper42`, `paper50` | Various GeoVac internal references. | **CITE-OK** (internal) | All present in corpus. |
| `paper24` | Four-layer Coulomb/HO asymmetry, §V. | **CITE-OK** (internal) | Confirmed present. |

### Problems found (CITE-MISATTRIBUTED / DOESNT-SUPPORT / CANT-FIND)

#### F2 [HIGH] — `bousso_casini_fisher_maldacena2020`: wrong author list (Fursaev-failure-mode)

**Exact phrase affected (bibitem):**
> "R. Bousso, H. Casini, Z. Fisher, J. Maldacena, 'Gravity dual of Connes cocycle flow,' Phys. Rev. D 102 (2020), 066008. arXiv:2007.00230."

**What I verified (WebFetch on arXiv:2007.00230):**
> Title: "Gravity Dual of Connes Cocycle Flow"
> Authors: Raphael Bousso, **Venkatesa Chandrasekaran, Pratik Rath, Arvin Shahbazi-Moghaddam**.

**Affected use in paper 49:**
- §8.1 "The Bousso–Casini–Fisher–Maldacena conjecture" (subsection title)
- §8.2 "BCFM (continuum, AdS₄ / S³)" (acronym throughout §8)
- §8.4 "BCFM 2020 paper" (positioning)
- Abstract / Conclusion (BCFM acronym propagates)
- bibitem `bousso_casini_fisher_maldacena2020`

**Recommended fix:** rename `bousso_casini_fisher_maldacena2020` → `bousso_chandrasekaran_rath_shahbazi_moghaddam2020` (or `bousso_etal2020`). Update author list. Replace "BCFM" acronym with "BCRS-M" or simply "Bousso et al. 2020" throughout. This is exactly the failure mode CLAUDE.md §3 records as "Fursaev–Solodukhin / `hep-th/9512134` failure mode."

#### F3 [HIGH] — `mondino_samann2025_pointed`: bibitem title is a paraphrase, not the published title

**Bibitem reads:** "Pointed Lorentzian Gromov–Hausdorff convergence of covered pre-length spaces."

**Actual arXiv:2504.10380 v4 title:** "Lorentzian Gromov-Hausdorff convergence and pre-compactness." (Confirmed via WebFetch.)

This is propagated from Paper 48 (where the same fake title appears). The paper does discuss pointed LGH and covered LPLS as concepts, but the title in the bibitem is fabricated as a paraphrase. **This is the Fursaev pattern again at the title level.**

**Recommended fix:** correct the bibitem title to the published one across Papers 48 and 49 simultaneously.

#### F1 [HIGH] — `lindblad1975` and `uhlmann1977`: cited as support for a "chain inequality" that the cited papers do not prove

**Cited claim in §5.4 (Theorem 3.3 proof):** "The chain inequality follows from the triangle-type relation for relative entropies under composition of quantum channels" with the inequality stated as
$$ \Delta S^{1\to 3} \le \Delta S^{1\to 2} + \Delta S^{2\to 3}, $$
attributed to Lindblad/Uhlmann via the "data-processing inequality."

**What Lindblad 1975 and Uhlmann 1977 actually prove:** monotonicity of relative entropy under (trace-preserving) completely positive maps: $S(\Phi\rho \| \Phi\sigma) \le S(\rho\|\sigma)$.

**What is NOT in either reference:** a triangle/chain inequality $S(\rho_1\|\rho_3) \le S(\rho_1\|\rho_2) + S(\rho_2\|\rho_3)$. Literature confirms (Wilde's QIT lecture notes, arXiv:1909.05826, arXiv:2510.16918) that **quantum relative entropy does NOT satisfy a triangle inequality**. The max-relative entropy does; the standard relative entropy does not. Chain rules in the literature (1909.05826) are stated for channel relative entropies in multipartite settings, not for three arbitrary states.

This is **CITE-DOESNT-SUPPORT** combined with a real **logical gap in the proof** (Pass A finding E in claim A4). The chain inequality, as stated, is not a consequence of the cited papers and is not a known general fact about quantum relative entropy.

The paper does qualify with "the cocycle entropy production is **bounded above** by this relative entropy" (§5.4), but this upper bound by itself does not give a chain inequality on the entropy productions. If a → b, b → c, a → c each has its own ΔS bounded above by the corresponding relative entropy, the bound provides no direct comparison between ΔS^{1→3} and ΔS^{1→2} + ΔS^{2→3}.

**Recommended fix:** EITHER (a) restate Theorem 3.3 as conditional on a hypothesized chain property of the specific cocycle entropy production (with the chain inequality identified as a structural conjecture, not a theorem), OR (b) provide an actual proof of the chain inequality from a different ingredient (e.g., the explicit TICI cocycle composition with a Pinsker-style bound, or a specific construction on the BW substrate). The current proof structure is incomplete.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|------|------|------|------|------|
| "To the author's knowledge this is the first explicit identification of an information-theoretic underpinning for the reverse triangle inequality of synthetic Lorentzian geometry." | Abstract; §1.3; §5.5; §6.1 | "data processing inequality" + "reverse triangle" + "Lorentzian" + "synthetic"; "cocycle Radon-Nikodym + super-additivity"; arXiv search for relative entropy on causal pre-length spaces. | I found no direct prior art making this specific connection. However, the BCFM/kink-transform line (arXiv:2007.00230) does identify the Connes cocycle flow with bulk-gravity kink transforms, and the Bousso-Faulkner et al. QNEC line ties relative entropy to causal geometry in QFT. The exact GeoVac framing (synthetic Lorentzian reverse triangle ← Lindblad DPI) is not in those works. **Acceptable as "to the author's knowledge"** — but conditional on F1 being closed. | OK as hedged. |
| "OSLPLS as a category is genuinely new: no published 'operator-system Lorentzian pre-length space' concept exists as of May 2026." | §1.5 (related work) | Searched math.OA for "operator-system Lorentzian pre-length space"; checked Martinetti 2603.03216, Nieuviarts 2512.15450, Bizi-Brouder-Besnard 2018, van den Dungen 2016 (Krein spectral triples line). | No published concept with this exact framing found. **Acceptable as "to the author's knowledge"** — but the concept is essentially a relabeling of "operator system + KMS state + modular flow" which is well-known machinery; the novel content is the *categorical framing*, not the underlying objects. | Softening recommended: "we introduce the term OSLPLS for the operator-system-substrate-plus-modular-time-separation construction; while each ingredient is standard, the categorical framing as a non-commutative LPLS extension is, to our knowledge, new." |
| "Triple-intersection cocycle identity (TICI), implicit in Connes 1973 (C4) but not previously isolated as a stack-consistency property." | §5.3, §6 | Searched "triple intersection cocycle" + Connes + Radon-Nikodym; checked Petz-Hiai / Ohya-Petz literature on quantum information divergence + cocycle. | No direct match for "TICI as stack-consistency property." But this is essentially Connes (C4) applied iteratively — it's a direct corollary, not new content. **The novelty claim is weak** — it amounts to giving a name to a direct corollary of an established theorem. | Softer framing: "We observe that the cocycle composition rule (C4) of Connes (1973), applied to triple intersections, yields a stack-consistency property which we will refer to as TICI; this property may not have been emphasized in this form elsewhere, though it is a direct consequence of (C4)." |

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| F1 | Logical gap in Theorem 3.3 (`thm:strict_super_additivity`) proof: "Lindblad/Uhlmann chain inequality" $\Delta S^{1\to 3} \le \Delta S^{1\to 2} + \Delta S^{2\to 3}$ is asserted with citations Lindblad 1975 + Uhlmann 1977, but neither paper establishes such a triangle inequality and quantum relative entropy is known NOT to satisfy a triangle inequality. | A & B | E (proof gap) + CITE-DOESNT-SUPPORT | **HIGH** |
| F2 | `bousso_casini_fisher_maldacena2020` bibitem: wrong author list. Actual authors of arXiv:2007.00230 are Bousso, Chandrasekaran, Rath, Shahbazi-Moghaddam (not Casini, Fisher, Maldacena). Fursaev-failure-mode pattern; impacts §8 throughout (BCFM acronym, subsection title). | B | CITE-WRONG-METADATA (near-MISATTRIBUTED) | **HIGH** |
| F3 | `mondino_samann2025_pointed` bibitem title is a paraphrase, not the published title. Actual title is "Lorentzian Gromov-Hausdorff convergence and pre-compactness." Propagates from Paper 48. | B | CITE-WRONG-METADATA | **HIGH** |
| F4 | §9.1 framing inconsistent with abstract: section titled "Q2′ ... CLOSED" but abstract says "Paper~49 does *not* close Q2$^{\prime}$." Cannot have both. | A | C (overstatement / inconsistency) | **MEDIUM** |
| F5 | `connes_rovelli1994` minor metadata: page range 2899-2917 vs ADS 2899-2918; title "Generally Covariant" vs original "General Covariant." Cosmetic. | B | CITE-WRONG-METADATA (cosmetic) | **LOW** |
| F6 | TICI novelty framing: presented as new isolation, but is direct corollary of Connes (C4). | B (novelty) | Recommend soften | **MEDIUM** |
| F7 | Uhlmann panel "verification" of strict super-additivity is on triples selected by a greedy algorithm screening for three different modular orbits — does not verify the theorem (all triples), only that some positive-deficit triples exist. Honestly framed in §7.6, but the abstract states "Uhlmann monotonicity is verified numerically" without the selection caveat. | A | C (mild overstatement in abstract) | **MEDIUM** |
| F8 | "Bit-exact agreement with Paper 45 reference" Λ panel is tautological under "free-upgrade" reading (Paper 46 Thm 3.5 forces equality by construction). Internal-consistency only, not novel verification. | A | B (circularity to flag) | **LOW** (already documented as "free upgrade" in §7.1 honest scope; not over-stated in body) |
| F9 | `lindblad1975` and `uhlmann1977` citations at "data-processing inequality" framing are valid for monotonicity but NOT for the chain inequality. (Same as F1; listed here in citation table for completeness.) | B | CITE-DOESNT-SUPPORT (for chain inequality) | **HIGH** (counted under F1) |
| F10 | `bratteli_robinson1981` Theorem 5.3.10 cited but not verified physically against the book. Likely correct (the book does contain modular theory at this level). | B | CITE-OK (provisional) | **LOW** |

**Summary counts (Pass A):**
- A: 1 (Riemannian-limit recovery)
- B: 6 (internal-consistency results)
- C: 4 (overstatements / inconsistencies)
- D: 2 (unverifiable novelty)
- E: 1 (logical gap in Thm 3.3)

**Summary counts (Pass B):**
- CITE-OK: 11 (Connes 1973, Connes-Rovelli 1994, Bratteli-Robinson, Connes-vS, Kunzinger-Sämann, Latrémolière 2018, Che-Perales-Sormani, Minguzzi-Suhr, Sakovich-Sormani, Wilde, internal GeoVac papers)
- CITE-WRONG-METADATA: 3 (BCFM, Mondino-Sämann title, Connes-Rovelli minor)
- CITE-MISATTRIBUTED: 0 (BCFM is wrong-metadata, not MIS-attributed — the arXiv ID points to the right paper, just the authors are wrong)
- CITE-DOESNT-SUPPORT: 2 (Lindblad, Uhlmann at chain-inequality use)
- CITE-CANT-FIND: 0

**Severity totals:**
- HIGH: 3 (F1, F2, F3)
- MEDIUM: 3 (F4, F6, F7)
- LOW: 3 (F5, F8, F10)

---

## Broadcast readiness: **RED**

This paper has three HIGH-severity findings that block public broadcast as-is:

**F1 is the most consequential.** Theorem 3.3 (`thm:strict_super_additivity`) is the substantive new mathematical content of the paper — it is the load-bearing claim that makes Paper 49 a genuine extension of Paper 48 and the centerpiece of the "twin-paradox-as-quantum-information" headline. Its proof invokes a "chain inequality" on quantum relative entropy that **is not a known general theorem** and is **not in the cited references** (Lindblad 1975, Uhlmann 1977). The paper does say the cocycle entropy production is "bounded above" by relative entropy, but the proof structure of Theorem 3.3 (eq:lindblad_uhlmann + eq:deficit) requires an actual chain inequality on the ΔS's themselves, which is asserted without justification. A domain expert in quantum information theory would flag this immediately. The paper should NOT be broadcast under "Bridge Theorem 6.4′-Q1′ theorem-grade rigor" framing until either (a) the chain inequality is given a real proof from this substrate's specific structure (a candidate: explicit form of the cocycle on the GeoVac substrate via Paper 42's integer-spectrum K_α), or (b) Theorem 3.3 is downgraded to a conjecture conditional on a substrate-specific chain property.

**F2 is the Fursaev pattern.** The BCFM bibitem incorrectly lists Casini, Fisher, Maldacena as co-authors of arXiv:2007.00230, when the actual co-authors are Chandrasekaran, Rath, Shahbazi-Moghaddam. The §8 "Bousso-Casini-Fisher-Maldacena conjecture" subsection title and the BCFM acronym throughout §8 propagate this error. A domain expert in AdS/CFT will spot this in 30 seconds and the credibility cost is high. This is exactly the failure mode CLAUDE.md §3 records for `hep-th/9512134`.

**F3 propagates from Paper 48.** Both Paper 48 and Paper 49 cite arXiv:2504.10380 with a paraphrased "title" that is not the actual title. Fix in both papers simultaneously.

The MEDIUM-severity findings (F4 abstract/§9.1 inconsistency on Q2′, F6 TICI novelty softening, F7 numerical panel framing) are framing-level and should be addressed in an errata pass after the HIGH-severity fixes.

The bones of the architecture (OSLPLS category, axiom transport, embedding functor, propagation number, Riemannian-limit recovery, Λ free-upgrade panel) appear sound and internally consistent. The honest scope sections (§5.7, §7.6, §9) are well-written and openly acknowledge limitations. The Q1′-Light → Phase-1 → Phase-2 narrative is coherent. The paper has good GeoVac discipline overall — the failures are at the two interface points to external work: (i) the citation-author hygiene on BCFM and Mondino-Sämann, and (ii) the application of Lindblad/Uhlmann to a setting they were not designed for.

**Recommended path:** address F1 in a dedicated sub-sprint (close the chain inequality at the cocycle-entropy-production level on the GeoVac substrate; if unfixable, demote Theorem 3.3 to a conjecture); fix F2 + F3 bibitem entries (mechanical, 30 minutes); reconcile abstract / §9.1 on Q2′ (F4, 5 minutes). After these, the paper should re-pass to YELLOW or GREEN.

---

## What I could NOT verify (hand to a human expert)

1. **Whether the cocycle entropy production on the GeoVac substrate satisfies a chain inequality even though general quantum relative entropy does not.** This is the central question for closing F1. A domain expert in modular theory + quantum information theory should be consulted. (Possible angle: the substrate-specific cocycle on a single wedge KMS state with finite modular spectrum may admit a stronger inequality than the general case; this would be substrate-specific and would be the substantive new content.)

2. **The strict super-additivity claim on all three-orbit triples** (versus the selected triples in the numerical panel). The Uhlmann monotonicity panel is empirically positive on a few chosen triples; whether it holds for ALL three-orbit triples is the actual theorem and requires either a structural proof on the substrate or a much larger random panel.

3. **Whether the "TICI" property is genuinely new or trivially known to operator-algebra specialists.** Pass-B novelty audit cannot resolve this; a Connes-vS specialist should be asked.

4. **The bulk-gravity-dual reading (§8) at the operator-system level.** Pure conjecture in the paper, clearly labeled; consultation with the BCFM authors (Bousso et al.) would be the only way to assess.

5. **Bratteli-Robinson Theorem 5.3.10 exact statement** — should be physically verified against the book to confirm it supports the use in Lemma 3.1 (`lem:orbit_kms`).

6. **Whether `Mondino-Sämann` 2504.10380 v4 contains the Def 3.8 / Def 3.12 numbering used in Paper 49's §2.5 transcription.** If the numbering does not match the actual paper, this would be additional CITE-WRONG-METADATA. Worth confirming.

---

## Summary line for dispatcher

**Paper 49** — Operator-system Lorentzian pre-length spaces and the strong-form Krein–MS bridge — **RED** — A: 1 / B: 6 / C: 4 / D: 2 / E: 1 — CITE-OK: 11, CITE-WRONG-METADATA: 3, CITE-DOESNT-SUPPORT: 2 — HIGH: 3, MEDIUM: 3, LOW: 3 — **Top finding:** Theorem 3.3's proof invokes a "Lindblad/Uhlmann chain inequality" on quantum relative entropy that is not in the cited references and is known NOT to hold in general (quantum relative entropy does not satisfy a triangle inequality), making the load-bearing strict-super-additivity claim of the paper resting on an unproven step.
