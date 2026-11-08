# Confidence Review: Paper 16 — The Periodic Table Recast as $S_N$ Representation Theory on $S^{3N-1}$

**File audited:** `papers/group4_quantum_computing/paper_16_periodicity.tex`
**Date:** 2026-06-01, Wave 1 re-fire (text-level audit, no driver re-runs)
**Auditor:** Confidence Reviewer (combined Pass A + Pass B)

## Calibration check

Not a calibration run.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|:-:|:------|:---------|:--------|:---------|:--------------------|
| 1 | $\Lambda^2 Y_\nu = \nu(\nu+3N-2) Y_\nu \equiv 2\mu_{\rm free} Y_\nu$ (SO$(3N)$ Casimir) | Eq. (2), §III.A | **A** | EXTERNAL | Standard SO$(3N)$ representation theory; matches Paper 7 Eq. and `Avery2000` convention. |
| 2 | Theorem 1: for any $N$-electron atom with $S<N/2$, spatial Young diagram has $\lambda_1=2$ | §III.B Thm 1 | **A** | EXTERNAL | Derivation correct: spin diagram $[N/2+S, N/2-S]$ with $n_\downarrow\ge 1$ has conjugate-diagram first column of height 2. Schur–Weyl duality + Young-diagram conjugation is textbook. Verified by hand for $N=2,3,4,7,10,11,12$. |
| 3 | $\nu_{\rm min}=N-\lambda_1$ (minimum allowed $\nu$ for spatial irrep of given $\lambda_1$) | Eq. (5), §III.C | **B** | GEOVAC-ONLY / MIXED | Cited to Avery 2000 §18. Could not independently fetch Avery's book; the convention is plausible and consistent with Paper 7's helium and lithium cases ($N=2$, $\lambda_1=2$, $\nu=0$; $N=3$, $\lambda_1=2$, $\nu=1$, $\mu=4$). Flagged: this result is load-bearing for the universal $\nu=N-2$ headline; an explicit reference page would help. |
| 4 | **Universal result $\nu=N-2$ for $S<N/2$** (Eq. 6, boxed) | §III.C | **A (subject to claim 3)** | MIXED | Follows immediately from claims 2 and 3. If Avery's $\nu_{\rm min} = N - \lambda_1$ is correct, then $\nu = N-2$ for all $S<N/2$ is forced. Headline derivation is sound. |
| 5 | **Boxed formula $\mu_{\rm free} = 2(N-2)^2$** (Eq. 4) | §IV Eq. (4) | **E (WRONG)** | EXTERNAL contradicts | The algebra in Eq. (4) is wrong. Substituting $\nu=N-2$ into the (correct) Casimir formula in Eq. (2) gives $\mu_{\rm free} = (N-2)(4N-4)/2 = 2(N-2)(N-1)$, NOT $2(N-2)^2$. The author appears to have simplified $(4N-4)/2$ as $2(N-2)$ instead of $2(N-1)$. Verified algebraically with sympy and by hand. **The first algebra step in Eq. (4) is correct; the final boxed simplification is wrong.** |
| 6 | Table I numerics (all rows except H/He) | §IV Table I | **E (CASCADING)** | EXTERNAL contradicts | Every $N\ge 3$ row of Table I encodes $\mu = 2(N-2)^2$ instead of the correct $\mu = 2(N-2)(N-1) = \nu(\nu+3N-2)/2$. Verified table-by-table. Examples: Li (paper 2, correct 4), Be (paper 8, correct 12), Ne (paper 128, correct 144), Ar (paper 512, correct 544), …, Og (paper 26912, correct 27144). |
| 7 | Cr (N=24): $\nu=22$, $\mu_{\rm free}=1012$ | §VII.B (transition metals discussion) | **A (in itself), but INTERNALLY INCONSISTENT with claims 5/6** | EXTERNAL | $\mu = 22\cdot(22+70)/2 = 22\cdot 92/2 = 1012$. Verified. **This is the correct Casimir value.** Under the paper's own boxed Eq. (4), Cr would have $2(24-2)^2 = 968$, not 1012. The text uses the correct Casimir formula here while Eq. (4) and Table I use the wrong one — direct internal contradiction. |
| 8 | $\mu_{\rm free}/N^2 \to 2$, $\mu_{\rm free}/\binom{N}{2}\to 4$, $\nu/(3N-2) \to 1/3$ as $N\to\infty$ | §IV.A | **A (limits) / B (rate of approach)** | EXTERNAL | All three asymptotic limits are correct under EITHER the wrong or the correct $\mu$ formula. The limits survive the algebra bug. (Numerical convergence is slower under the correct formula but still $\to 4$.) |
| 9 | $\ell_{\rm eff}=(3N-3)/2$, $n_{\rm eff}=(3N-1)/2$ | Eq. (3) | **A** | EXTERNAL | Standard hyperspherical centrifugal extraction. Matches Paper 7 Table III. |
| 10 | Four-type classification A/B/C/D (+E for open shells) via $S_N$ irrep shape | §V | **C (overstated)** | GEOVAC-ONLY framing | The mathematical identification of spatial-irrep shape with chemical type is correct and a legitimate reformulation. The abstract's claim that it "recasts the periodic table's group structure in representation-theoretic language" is honest (matches the body). **BUT**: Types B and D are both $[2^{N/2}]$ — the spatial irrep shape does NOT distinguish them — the difference is "occupancy pattern" (an empirical input), as §V.A admits explicitly. The §I-B summary "the shape of the $S_N$ irrep classifies atoms into chemical types" overreaches what §V actually shows, and §V.A correctly admits the limitation. Already self-flagged in the body. |
| 11 | Periodic law as repetition of irrep sequence $C\to D \to E\cdots\to B$ | §VI | **C (already correctly hedged)** | EXTERNAL/MIXED | The body explicitly notes "This is a reformulation, not a derivation from first principles: the shell-filling order itself … is an empirical input." That hedge survives in the abstract too ("This classification maps onto, rather than predicts, the known periodic table structure"). Good honesty. |
| 12 | Dirac instability is a metric (not topological) singularity on $S^{3N-1}$; topology smooth through $Z=137$ | §VII | **A** | EXTERNAL | The topological smoothness through $N=137$ is trivially true (none of $\nu$, $\mu$, the irrep, the manifold dimension is singular at $N=137$). The relativistic conformal singularity at $Z\alpha\to 1$ is standard QED. Sound separation. |
| 13 | Cr $3d^5 4s^1$ has 141 exchange pairs; $3d^4 4s^2$ has 136; difference 5 | §VII.B | **A** | EXTERNAL | Verified: $\binom{15}{2}+\binom{9}{2} = 105+36 = 141$; $\binom{14}{2}+\binom{10}{2} = 91+45 = 136$. |
| 14 | $\delta_{1s}=(Z\alpha)^2/(1-(Z\alpha)^2)$ values (Table II) | §VII.A Table II | **A (with one minor cosmetic)** | EXTERNAL | Recomputed: $Z=29 \to 0.047$ ✓, $Z=47\to 0.133$ ✓, $Z=79\to 0.498$ ✓, $Z=118\to 2.87$ ✓. Z=1 row shows $\delta_{1s}=0.00005$ but recomputation gives $\sim 1.06\times 10^{-4}$ (paper conflated $(Z\alpha)^2$ with $\delta$ — both columns say $5\times 10^{-5}$). At $Z=137$ paper says $\delta \sim 10^3$; actual value 1903 — both are "diverging" so OK qualitatively. |
| 15 | Eq. (9) Dirac 1s formula $E = mc^2(\sqrt{1-(Z\alpha)^2}-1)$ | §VII.A Eq. (9) | **A** | EXTERNAL | Standard Dirac 1s point-nucleus result. |
| 16 | Connection-to-Paper-0 claim: "Paper 0 derived the $2n^2$ orbital degeneracy from the geometric packing of nodes on the discrete $S^3$ lattice" | §VII.A | **C (mild)** | GEOVAC-ONLY | Paper 0 §III.A explicitly states the $2n^2$ grouping is "not uniquely forced by the two axioms" — it is a structural correspondence, not a derivation. Paper 16 slightly overstates Paper 0's own qualification. Recommend "Paper 0 motivated" or "Paper 0 obtained the $2n^2$ count as a structural correspondence". |
| 17 | "periodic law itself … has not been derived from a single mathematical principle" | §I | **C (overstated)** | EXTERNAL | Prior group-theoretic frameworks (Rumer–Fet group $SO(4,2)\times SU(2)$; arXiv:2501.18272 SO(4,4); MDPI Symmetry 14:137) exist. Recommend acknowledging this lineage in §II Known Results. (Inherited from prior CITE report.) |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value | survives? |
|:------|:---------------|:----------------------|:--------------------|:----------|
| $\mu_{\rm free}$ for Li ($N=3, \nu=1$) | 2 | Casimir $\nu(\nu+3N-2)/2$ | 4 | **NO** (paper Eq. 4 + Table I row wrong) |
| $\mu_{\rm free}$ for Ne ($N=10, \nu=8$) | 128 | Casimir | 144 | **NO** |
| $\mu_{\rm free}$ for Ar ($N=18, \nu=16$) | 512 | Casimir | 544 | **NO** |
| $\mu_{\rm free}$ for Og ($N=118, \nu=116$) | 26 912 | Casimir | 27 144 | **NO** |
| $\mu_{\rm free}$ for Cr ($N=24, \nu=22$) | 1012 (in §VII.B) | Casimir | 1012 | **YES** — paper text is right; **inconsistent with paper Table I** which would give 968 |
| Cr exchange pair counts (141 / 136 / 5) | 141, 136, 5 | direct combinatorics | 141, 136, 5 | **YES** |
| Per-pair limit $\mu/\binom{N}{2} \to 4$ | 4 | algebra | 4 | **YES** (true under either formula) |
| Per-electron limit $\mu/N^2 \to 2$ | 2 | algebra | 2 | **YES** (true under either formula) |
| Table II $\delta_{1s}$ values | (see table) | $(Z\alpha)^2/(1-(Z\alpha)^2)$ | match within precision except $Z=1$ cosmetic | **YES** modulo Z=1 cosmetic |
| Dirac 1s at $Z=137$ | "$\sim 10^3$" | direct | 1.9×10³ | **YES** qualitatively |

### Circularity map

**GEOVAC-ONLY anchors:**

- Eq. (5) ($\nu_{\rm min} = N - \lambda_1$) is cited to `Avery2000` (an external book), but the specific page/equation is not given and I could not independently fetch it. The universal-$\nu$ headline of the paper rides entirely on this one citation being correct. If Avery actually proves a different form of $\nu_{\rm min}$ (e.g., dependent on the full Young diagram shape, not just $\lambda_1$), the entire universal-$\nu$ derivation would need revisiting. **Recommend:** add page reference (Avery 2000, ch. X, eq. Y) for the load-bearing $\nu_{\rm min} = N - \lambda_1$ identity.

- The "topology smooth through $Z=137$" structural claim is internally consistent — it is just the observation that none of $\nu$, $\mu$, the irrep, etc. is singular at any finite $N$. No external dependency. **EXTERNAL-trivial**, not GeoVac-only.

**MIXED chains:**

- The four-type A/B/C/D classification is built from (i) the standard Young-diagram conjugation theorem (EXTERNAL), (ii) the universal-$\nu$ identification of §III (MIXED, depends on $\nu_{\rm min}=N-\lambda_1$), (iii) the empirical assignment of "noble gas," "alkali," "alkaline earth" labels (EXTERNAL: periodic table). The classification is a re-organization of empirical fact, not a derivation, as §VI honestly admits.

### Overstatement findings

| exact phrase | location | suggested honest replacement |
|:-------------|:---------|:------------------------------|
| "the periodic law itself---the recurrence of chemical properties with increasing atomic number---has not been derived from a single mathematical principle." | §I, ¶2 | "the periodic law … has been approached from several group-theoretic perspectives (e.g., the Rumer–Fet $SO(4,2)\times SU(2)$ framework), each capturing different aspects. The present work develops the complementary $S_N$ permutation-symmetry approach." (Adds one acknowledgment sentence; does not undercut the claim's spirit.) |
| §VII.A "Paper 0 derived the $2n^2$ orbital degeneracy from the geometric packing of nodes on the discrete $S^3$ lattice" | §VII.A | "Paper 0 obtained the $2n^2$ orbital count as a structural correspondence with the discrete $S^3$ packing" (matches Paper 0's own caveat in its §III.A that the cumulative grouping is "not uniquely forced by the two axioms"). |

The abstract's framing is otherwise well-hedged: "This classification maps onto, rather than predicts, the known periodic table structure" is appropriately honest. The "appears to be new" hedge on the universal $\nu=N-2$ identification is at the correct honest-ceiling.

## Pass B — Citation and novelty

### Citation table

| key | claimed as | verdict | source |
|:----|:-----------|:--------|:-------|
| `Mendeleev1869` | Periodic table established empirically 1869 | **CITE-OK** | Mendeleev, *Zh. Russ. Khim. Obshch.* vol. 1, 60 (1869). Historical-record confirmed. |
| `Fock1935` | Hyperspherical / SO(4) hydrogen | **CITE-OK** | *Z. Phys.* **98**, 145 (1935). DOI 10.1007/BF01336904. |
| `Fock1954` | Hyperspherical $N$-body framework | **CITE-OK-UNVERIFIABLE-PRIMARY** | *Izv. Akad. Nauk SSSR Ser. Fiz.* **18**, 161 (1954). Standard, frequently cited reference; could not fetch the primary directly. |
| `Avery2000` | Hyperspherical formalism; $\nu_{\rm min}=N-\lambda_1$ | **CITE-OK (existence) / CITE-CANT-INDEPENDENTLY-VERIFY (specific formula)** | Book confirmed: Kluwer, Progress in Theoretical Chemistry and Physics vol. 4, ISBN 9780792360872. The book IS the standard hyperspherical-Sturmian reference; the specific formula's location in the text not independently confirmed. |
| `Weyl1946` | Schur–Weyl duality | **CITE-OK** | Princeton, 2nd ed. 1946. |
| `FultonHarris1991` | Schur–Weyl duality (paired) | **CITE-OK** | Springer GTM 129. |
| `loutey_paper0` | GeoVac internal | **N/A** | (PI internal cross-ref.) |
| `loutey_paper1` | GeoVac internal | **N/A** | |
| `loutey_paper7` | GeoVac internal | **N/A** | |
| `loutey_paper13` | GeoVac internal | **N/A** | |
| `KnowlesHandy1984` | (in bibliography but **no `\cite` in body**) | **ORPHAN** | *Chem. Phys. Lett.* **111**, 315 (1984). The citation itself is correct; it is just unused in the text. |
| `Dirac1928` | Dirac equation | **CITE-OK** | *Proc. R. Soc. London A* **117**, 610 (1928). |
| `Greiner1985` | $Z\sim 137$ and $Z\sim 170$ Dirac dive limits | **CITE-OK** | Standard reference for supercritical QED. |
| `Pyykkoe2012` | Periodic table extension to $Z\le 172$ | **CITE-WRONG-METADATA (bib key only, cosmetic)** | Pyykkö, *Phys. Chem. Chem. Phys.* **13**, 161 (2011), DOI 10.1039/C0CP01575J. The bib *text* correctly shows 2011; the bib *key* `Pyykkoe2012` should be `Pyykkoe2011`. |
| `Schwerdtfeger2015` | "Relativistic effects in superheavy elements" book chapter | **CITE-CANT-FIND / likely CITE-MISATTRIBUTED** | (See Problem 1 below.) Re-flag per dispatcher request — the PI choice between Pershina 2014 chapter or Schwerdtfeger NPA 944 (2015) was deferred and the .tex still has the suspect entry verbatim. |

### Problems found

#### Problem 1 (still present — re-flag) — `Schwerdtfeger2015` appears confabulated

**Citation as printed in `paper_16_periodicity.tex` lines 807–811:**
```
P. Schwerdtfeger, L. F. Pashtedsky, A. Pernpointner, and B. Fricke,
"Relativistic effects in superheavy elements,"
in The Chemistry of Superheavy Elements, 2nd ed.
(Springer, Berlin, 2015).
```

**Verifiable problems** (re-checked via WebSearch June 2026):
1. The 2nd edition of *The Chemistry of Superheavy Elements* (Springer) is edited by **Schädel & Shaughnessy**, **2014** (not 2015), DOI 10.1007/978-3-642-37466-1.
2. The chapter on theoretical chemistry of superheavy elements is **V. Pershina** alone (ch. 3), not a four-author Schwerdtfeger paper.
3. The "L.F. Pashtedsky" is a misspelling/confabulation of **L. F. Pašteka** (a known Schwerdtfeger collaborator).
4. The genuine 2015 Schwerdtfeger–Pašteka publication is in *Nuclear Physics A* **944**, 551 (2015) — a journal review, NOT a book chapter — with author list **Schwerdtfeger, Pašteka, Punnett, Bowman** (not Pernpointner/Fricke).
5. WebSearch June 2026 returns the *Nucl. Phys. A* paper as the canonical 2015 reference; no four-author book chapter with the cited title is found.

This is the same failure mode as `hep-th/9512134` (the Fursaev–Solodukhin confabulation in CLAUDE.md §3, dead-ends row). RED severity.

**Recommendation:** Replace with one of:
- V. Pershina, "Theoretical Chemistry of the Heaviest Elements," in Schädel & Shaughnessy (eds.), *The Chemistry of Superheavy Elements*, 2nd ed., Springer 2014, ch. 3 (DOI 10.1007/978-3-642-37466-1_3); or
- P. Schwerdtfeger, L. F. Pašteka, A. Punnett, P. O. Bowman, "Relativistic and quantum electrodynamic effects in superheavy elements," *Nucl. Phys. A* **944**, 551 (2015), DOI 10.1016/j.nuclphysa.2015.02.005.

The Pyykkö 2011 reference already supports the periodic-table-extension content. The replacement is needed primarily for the $Z\sim 170$ pair-creation limit, which can also be supported by the already-cited Greiner–Müller–Rafelski 1985.

#### Problem 2 — `Pyykkoe2012` bib key cosmetic

Cosmetic only; the bib text correctly shows 2011. Rename key to `Pyykkoe2011` at next pass.

#### Problem 3 — `KnowlesHandy1984` orphan

`\cite{KnowlesHandy1984}` does not appear in the body. Either delete the bibitem or insert a citation where it fits the narrative (e.g., FCI computational context).

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|:-----------------|:---------|:---------|:-----------------|:---------------|
| "This specific identification appears to be new." (universal $\nu=N-2$ for $S<N/2$) | §II.B item 1 | Web searches for "S_N representation theory" + "hyperspherical" + "periodic table"; for "minimum K hyperspherical Young diagram"; group-theoretic periodic-table literature | No direct prior art found for the $S_N$/Young-diagram $\nu=N-2$ identification. The neighboring literature uses $SO(4,2)\times SU(2)$ (Rumer–Fet) or $SO(4,4)$ rather than $S_N$. | KEEP "appears to be new" — at the honest ceiling. **Add** one sentence in §II Known Results acknowledging the Rumer–Fet / $SO(4,2)\times SU(2)$ / MDPI 14:137 / arXiv 2501.18272 lineage as complementary prior art. |
| "The four-type classification (A/B/C/D)…" | §II.B item 2 | Same searches + Young-diagram periodic table classifications | No prior Young-diagram-shape four-type classification found. | KEEP framing as "reformulation," already honest. |
| "$\mu_{\rm free}/\binom{N}{2}\to 4$" per-pair limit | §II.B item 3 | Searches for "Pauli centrifugal" + per-pair limits | No prior art found. Elementary consequence of universal-$\nu$ formula. | KEEP. |
| "The Dirac instability as metric vs.\ topological singularity" | §II.B item 4 | Prior framings of Dirac dive in terms of manifold topology | The physics is standard; the $S^{3N-1}$-topology-vs-metric framing is the novel contribution. | KEEP as appropriately scoped. |
| "the periodic law itself … has not been derived from a single mathematical principle" | §I, ¶2 | Group-theoretic derivations of periodic law | Rumer–Fet and SO(4,2)×SU(2) frameworks claim group-theoretic derivations | **MEDIUM/YELLOW** — recommend §II/§II.A acknowledgment of this lineage. |

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|:-:|:--------|:-----|:--------|:---------|
| 1 | **Eq. (4) boxed formula $\mu_{\rm free}=2(N-2)^2$ is algebraically wrong; correct value is $2(N-2)(N-1)$. Eq. (2) and Eq. (4) are inconsistent.** | A | E | **HIGH** |
| 2 | **Table I $\mu_{\rm free}$ column wrong from N=3 onward (cascades from Eq. (4))** | A | E | **HIGH** |
| 3 | **Internal inconsistency: §VII.B uses correct Casimir for Cr ($\mu=1012$), while §IV Eq. (4) + Table I use wrong formula (would give 968 for Cr)** | A | E | **HIGH** |
| 4 | `Schwerdtfeger2015` likely confabulated four-author book chapter; real 2015 Schwerdtfeger paper is *Nucl. Phys. A* with different authors | B | CITE-CANT-FIND / CITE-MISATTRIBUTED | **HIGH** |
| 5 | "Paper 0 derived $2n^2$ orbital degeneracy" overstates Paper 0's own qualified "structural correspondence" framing | A | C | LOW |
| 6 | "periodic law has not been derived from a single mathematical principle" ignores Rumer–Fet/SO(4,2)×SU(2) literature | A/B | C | MEDIUM |
| 7 | $\nu_{\rm min}=N-\lambda_1$ load-bearing identity cited to Avery 2000 without page/equation; could not independently verify | A/B | B (content) / CITE-CANT-INDEPENDENTLY-VERIFY | MEDIUM |
| 8 | Table II Z=1 row reports $\delta_{1s}=5\times 10^{-5}$ but $(Z\alpha)^2 = 5\times 10^{-5}$ — the two columns are conflated at Z=1; corrected value $\delta\sim 1.06\times 10^{-4}$ | A | E (cosmetic) | LOW |
| 9 | `Pyykkoe2012` bib key should be `Pyykkoe2011` (text already correct) | B | CITE-WRONG-METADATA | LOW |
| 10 | `KnowlesHandy1984` is an orphan bibitem (no `\cite` in body) | B | (orphan) | LOW |

**Totals — Pass A verdicts:** A: 11 · B: 2 · C: 3 · D: 0 · E: 4
**Totals — Pass B verdicts:** CITE-OK: 9 · CITE-WRONG-METADATA: 1 · CITE-CANT-FIND/MISATTRIBUTED: 1 · CITE-CANT-INDEPENDENTLY-VERIFY: 2 · ORPHAN: 1 · N/A (internal): 4
**Severity totals:** HIGH: 4 · MEDIUM: 2 · LOW: 4

## Broadcast readiness: **RED**

**RED — block until math errors fixed.**

Paper 16's content audit produced one structural finding that overrides the otherwise sound headline derivation: the boxed formula $\mu_{\rm free}=2(N-2)^2$ in Eq. (4) and every dependent row of Table I are mathematically wrong. The correct value is $\mu_{\rm free}=2(N-2)(N-1) = \nu(\nu+3N-2)/2$ at $\nu=N-2$. This cascades through the central numerical content of the paper. It is also internally inconsistent: §VII.B's Cr discussion uses the **correct** Casimir value 1012 (which is $22\cdot 92/2$, NOT $2\cdot 22^2 = 968$), so the paper is contradicting itself in print. The fix is a one-line algebra correction in Eq. (4), a one-column rewrite of Table I, and a sentence each in §II.B item 3 (per-pair limit) and §VII.A (the $N=137$ factoid which currently reads "$\mu_{\rm free}=36450$" but the Casimir actually gives 36720).

A domain expert reading Paper 16 would catch the boxed-formula error immediately by doing the substitution: $(4N-4)/2 = 2(N-1)$, not $2(N-2)$. This must be fixed before broadcast.

Independently, the `Schwerdtfeger2015` confabulation flagged in the prior CITE-only report is still present in the .tex and should be replaced (Pershina 2014 or Schwerdtfeger NPA 944 — the PI choice was deferred and not yet applied). This is the same failure mode as the Fursaev–Solodukhin / `hep-th/9512134` dead-end in CLAUDE.md §3.

Once the math fix and the citation replacement are applied, the paper is otherwise on solid ground: Theorem 1's proof of universal $\lambda_1=2$ is correct; the universal-$\nu$ headline survives unchanged (since $\nu=N-2$ is unaffected); the four-type A/B/C/D classification is honestly framed as reformulation; the per-pair → 4 and per-electron → 2 asymptotic limits survive under either formula; the Cr exchange-pair combinatorics are verified; and the Dirac-instability/topology-smoothness reading is standard physics. Two medium-severity framing items (Paper 0 overstatement; missing acknowledgment of Rumer–Fet lineage) can be batched into the errata pass alongside the math fix.

## What I could NOT verify (hand to a human expert)

- **Avery 2000 $\nu_{\rm min}=N-\lambda_1$ identity.** The load-bearing identity Eq. (5) is cited to Avery's 2000 book (Kluwer), which I could not fetch online. The whole universal-$\nu$ headline depends on this being the correct general statement. A reader with library access to Avery's book should confirm the page/equation. (Alternatively: re-derive the identity from Schur–Weyl + hyperspherical-harmonic transformation properties under $S_N$, independently of Avery's text.)
- **Fock 1954.** Soviet-era *Izv. Akad. Nauk SSSR* journal, not indexed online; cannot independently verify. The reference is standard in the few-body literature.
- **Domain priority claim.** Whether the universal $\nu=N-2$ / four-type classification is genuinely new in the $S_N$ framework (vs. covered in some unindexed prior book or paywalled book chapter) can only be settled by a domain expert in hyperspherical-Sturmian / atomic-structure group theory. The "appears to be new" hedge is at the honest ceiling for what a web search can support.

Sources consulted:
- [Schwerdtfeger et al. 2015, NASA-ADS](https://ui.adsabs.harvard.edu/abs/2015NuPhA.944..551S)
- [Schädel & Shaughnessy (eds.) 2014, Springer Link](https://link.springer.com/book/10.1007/978-3-642-37466-1)
- [Pershina 2014 chapter, Springer Link](https://link.springer.com/chapter/10.1007/978-3-642-37466-1_3)
- [Pyykkö 2011, RSC](https://pubs.rsc.org/en/content/articlehtml/2011/cp/c0cp01575j)
- [Periodic Table and SO(4,4), arXiv:2501.18272](https://arxiv.org/pdf/2501.18272)
- [SO(4,2)×SU(2) group-theoretic periodic system, MDPI Symmetry 14:137](https://www.mdpi.com/2073-8994/14/1/137)
