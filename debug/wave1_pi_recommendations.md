# Wave 1 — PI Recommendations

**Date:** 2026-06-01
**Scope:** Items from the always-load foundational papers (0, 1, 7, 14, 16, 22, 23, 24, 27, 31, 32) audit that need PI judgment.
**Companion:** `debug/wave1_salvage_memo.md` (partial-run record); `debug/review_paper{0,1,7,14,16,22,23,24,27,31,32}.md` (raw audit reports).

## 1. Mechanical fixes already applied (no decision needed)

20 edits across 8 files in this session:

| File | Edit | Source audit |
|:-----|:-----|:-------------|
| Paper 0 §V | H₂ recovery 94.1% → 96.0% (Paper 15 figure) | Paper 0 AUDITOR |
| Paper 7 §VI.B | S⁵ helium eigenvalues `0,3,8,15` → `0,6,16,30` (even-ν singlet) | Paper 7 AUDITOR |
| Paper 32 (×4) | Removed `perez_sanchez2024` from "Yang-Mills-without-Higgs correction" cite bundles | Paper 32 CITATION |
| Paper 32 bib | `bizi_brouder_besnard2018` title (right ID, wrong title) | Paper 32 CITATION |
| Paper 32 bib | `cacic2009` metadata (vol/pages/year) | Paper 32 CITATION |
| Paper 31 line 200 | `\cite{perez_sanchez2024}` → `\cite{perez_sanchez2025}` for "correction" claim | Paper 31 REVIEW |
| Paper 31 bib | Split `perez_sanchez2024` synthesized-title bibitem into two clean entries | Paper 31 REVIEW |
| Paper 31 bib | `bizi_brouder_besnard2018` title fix (same as Paper 32) | Paper 31 REVIEW |
| Paper 31 bib (×2) | Nieuviarts author `A.` → `G.` (Gaston) | Paper 31 REVIEW |
| Paper 30 line 139 | Removed `perez_sanchez2024` from "Y-M-without-Higgs" cite | corpus-wide grep |
| Paper 38 line 1049 | Removed `perez_sanchez2024` from "continuum-limit correction" cite | corpus-wide grep |
| Synthesis group1 | Removed `perez_sanchez2024` from "Y-M-without-Higgs clarification" cite | corpus-wide grep |
| Paper 32 bib | Nieuviarts `A.` → `G.` (1 instance, missed by Paper 32 audit) | corpus-wide grep |
| Paper 42 bib | Nieuviarts `A.` → `G.` (3 instances) | corpus-wide grep |
| Paper 43 bib | Nieuviarts `A.` → `G.` (3 instances) | corpus-wide grep |
| Paper 44 bib | Nieuviarts `A.` → `G.` (2 instances) | corpus-wide grep |

## 2. HIGH-severity items needing PI decision (the four REDs)

### R1. Paper 14 — `trenev2025` does not contain the headline Pauli counts

**Issue:** Paper 14's abstract claims "two or more orders of magnitude" advantage over Gaussian baselines. The load-bearing Gaussian Pauli counts (LiH/H₂O/BeH₂ at STO-3G/cc-pVDZ/cc-pVTZ: 276/5851/63519 and 551/8921/107382) are attributed to `trenev2025` (arXiv:2311.03719). Verified: arXiv:2311.03719 is about *vibrational structure of acetylene-like polyynes*, not LiH/H₂O electronic structure, and does not contain those Pauli counts.

**Decision needed:**
- Are these counts our own re-computation against re-derived Gaussian baselines?
- Or are they from a different published source (Lee/Berry/Gidney PRX Quantum 2, 030305 is the closest candidate from the bibliography)?
- Or a different Trenev paper that was misattributed?

**Recommendation:** Single targeted agent dispatch to identify the right source (or PI hand-supplies). Until resolved, the abstract's headline cannot be defended.

**Secondary:** Paper 14 §V (spinor section) is missing 7 bibitems entirely (`Sunaga2025`, `BJL`, `Szmytkowski2007`, `MartinezYRomero2004`, `Dyall`, `BreitPauli`, `ChildsBerry`). Mechanical to add once the right works are identified.

---

### R2. Paper 16 — Eq.(4) `μ_free = 2(N-2)²` is mathematically wrong

**Issue:** Substituting ν=N−2 into the S^(3N−1) Casimir formula gives `μ_free = 2(N-2)(N-1)`, not `2(N-2)²`. Verified at N=3 (correct 4, paper boxes 2) and N=4 (correct 12, paper boxes 8). The paper *contradicts itself*: §VII.B uses Cr's μ=1012 — incompatible with the boxed Eq.(4). Table I cascades from N=3 onward.

**Decision needed:**
- (a) **Drop boxed Eq.(4)** — the universal ν=N−2 derivation (Theorem 1) survives without it, since μ values can be quoted directly from the Casimir formula.
- (b) **Replace Eq.(4)** with `μ_free = 2(N-2)(N-1)` and recompute Table I.
- (c) **Re-derive** whatever motivated the `(N-2)²` form (was there a separate combinatorial argument?).

**Recommendation:** Option (b) — the formula is meaningful, just wrong, and recomputing Table I is mechanical once the formula is fixed. Then audit §VII.B for consistency.

**Secondary (held over from yesterday):** `Schwerdtfeger2015` book chapter does not appear to exist as described. Replace with **Pershina 2014** (the real chapter in the real *Chemistry of Superheavy Elements* 2nd ed.) or **Schwerdtfeger et al. *Nucl. Phys. A* 944, 551 (2015)** — depending on which content Paper 16 actually needs the cite for.

---

### R3. Paper 22 — 1.44% angular density uses stricter-than-Coulomb selection rule

**Issue:** The headline angular ERI density (1.44% at l_max=3) is computed under "pair-diagonal-m" (`m_a=m_c` AND `m_b=m_d` per pair), which is *strictly stronger* than the physical Coulomb rule `m_a+m_b=m_c+m_d` (total magnetic conservation). Under the physical rule, D(l_max=3) = 6.06% — ~4.2× higher. The paper's own §V Table II acknowledges the distinction, but the abstract / §III / §IV / §VII present 1.44% as *the* Coulomb angular density. Cascades into Papers 14, 23, 31 which cite 1.44% as load-bearing.

**Decision needed:**
- (a) **Defend pair-diagonal-m as the relevant sparsity** for the actual qubit-encoding case, because composed encoding fixes m per orbital block. Then the 1.44% is right for the qubit cost, and the abstract should say "pair-diagonal-m angular density" explicitly. Paper 14's downstream cost claim survives unchanged.
- (b) **Concede pair-diagonal-m is stricter** and update the headline to 6.06%. Paper 14's "two or more orders of magnitude" claim then needs to be re-checked under the looser rule.

**Recommendation:** This is a physics-judgment call only the PI can make. (My read: option (a) seems likely defensible — composed encoding does fix m per orbital — but the paper currently *neither* states this *nor* verifies the qubit cost depends on pair-diagonal-m and not on total-m conservation. Whichever option is right, the abstract needs to be unambiguous, and Paper 14's downstream count needs an audit under the chosen reading.)

**Secondary:** 2 CITE-CANT-FIND for Dyall and Grant in §V spinor section — same missing-bibitem failure mode as Paper 14.

---

### R4. Paper 23 — three CITE-MISATTRIBUTED precision-AMO citations

**Issue:** Three load-bearing citations in the +286 ppm chain breakdown (§X) point to real DOIs that resolve to different papers. Fursaev–Solodukhin failure mode in miniature, three times:

| Bibkey | What the paper says | What the DOI/PRL actually points to |
|:-------|:--------------------|:-------------------------------------|
| `Pachucki2023` | PRL 130, 023004 = "Three-Photon Exchange" by Pachucki–Patkóš–Yerokhin | PRL 130, 023004 = "Recoil Corrections" by Pachucki–Yerokhin |
| `PachuckiYerokhin2010` | PRL 104, 070403 = deuterium HFS paper | PRL 104, 070403 = "Fine Structure of Heliumlike Ions" |
| `Eides2024` | DOI 10.1016/j.physletb.2024.139049 = Eides hyperfine paper | Same DOI = black-hole physics by Hadi-Akbarieh |

The physics survives: 12 numerical claims reproduce bit-exact (Roothaan 5/8 + Taylor, λ_n=85.7, BS=2.72e-4, finite-size ΔE=1.01e-10, A_hf=2.16e-7, Eides Zemach −39.50 ppm, D BF HFS 327.3975 MHz at +40 ppm vs Wineland-Ramsey).

**Decision needed:** PI identifies (or directs a targeted agent to find) the correct citations for each.

**Secondary MEDIUM:** `FriarPayne2005` cites wrong journal (PRA→PRC) and wrong article number (014501→014002). Mechanical once verified.

---

## 3. Cross-paper structural items (HIGH/MEDIUM, deeper than mechanical)

### S1. Paper 18 §III vs §VIII Ω-normalization (blocks Paper 0 + Paper 7 broadcast)

**Issue:** Paper 18 §III defines $\Omega = (1+p^2/p_0^2)^{-1}$ (gives $\Omega(0)=1$). Paper 18 §VIII defines $\Omega = 2p_0/(p^2+p_0^2)$ (gives $\Omega(0)=2/p_0$). Paper 0 §VI.C and Paper 18 §VII assume $\Omega(0)=2$ (i.e., §VIII convention with $p_0=1$). The κ=−1/16 derivation in Paper 0 inherits this unstated double assumption.

**Recommendation:** Pick one canonical Ω convention in Paper 18 (or state both clearly with their relationship), then add a half-sentence in Paper 0 §VI.C noting the $p_0=1$ ground-state-shell convention. Mechanical fix once the convention is chosen.

### S2. Three-layer → five-layer Coulomb/HO asymmetry, corpus-wide

**Issue:** Paper 31 §VII describes the asymmetry as "three layers, not two." Paper 24 §V documents Layer 4 explicitly (modular-Hamiltonian Pythagorean orthogonality, Sprint L2-F.1). Paper 51 introduces Layer 5 (gravity termination via Bernoulli uniqueness, Sprint G4-5a, per memory `fifth-asymmetry-gravity-termination.md`). Paper 31 is two layers behind the corpus.

Plus Paper 24 itself has internal layer-count inconsistency:
- §IV footnote: "refined to five layers across Sprint 5 / Paper 30 / Sprint ST-SU3 / Paper 51"
- §V title: "Fourth layer: modular-Hamiltonian Pythagorean orthogonality" (this part is fine — it's the section that introduces Layer 4)
- §V conclusion line 699: "The four-layer Coulomb/HO asymmetry is consistent..."

**Decision needed:** What is Paper 24's scope?
- (a) Paper 24 covers Layers 1–4 only; Layer 5 belongs to Paper 51. Then §IV footnote should say "four layers in Paper 24, with Layer 5 added by Paper 51"; §V conclusion is correct as-is.
- (b) Paper 24 should cover all five (add a §V.6 "Fifth layer: gravity termination" subsection). Then §V conclusion needs "five-layer."

**Recommendation:** Option (a) — Paper 24's natural scope is the asymmetry up through modular structure; gravity is Paper 51's territory. After §IV footnote softening, Paper 31's "three-layer" section becomes the rewrite target.

**Paper 31 §VII rewrite scope:** Convert the "Three-Layer Sharpening (Sprint 5)" section into "Five-Layer Sharpening (Sprint 5 through Sprint G4-5a)" with two additional table rows for modular-Hamiltonian (Paper 24 §V) and gravity termination (Paper 51). Plus update the table caption, the §I forward reference, and the §IX summary that cite "three-layer."

### S3. Paper 32 abstract α figure + case-exhaustion theorem scope (carryover from yesterday)

**S3a — Abstract α figure:** "8.8×10⁻⁸ matches α⁻¹" is the post-Sprint-A residual $|K - 1/\alpha - \alpha^2|/(1/\alpha)$, not the raw match $|K-1/\alpha|/(1/\alpha) \approx 4.77\times 10^{-7}$. The body is honest in §VI Rem K_under_theorem; the abstract overshoots.

Suggested replacement options (auditor's text, lightly edited):
- (i) "the empirically interesting property that the residual $K - 1/\alpha - \alpha^2$ between one of its spectral-action coefficient combinations and $\alpha^{-1}$ (after subtraction of the leading $\alpha^2$ Sprint-A correction) is $1.2\times10^{-5}$, a relative residual of $8.8\times10^{-8}$ that matches $\pi^3\alpha^3$ to 0.25%."
- (ii) "matches $\alpha^{-1}$ at approximately $5\times10^{-7}$ raw, with a Sprint-A residual term $\pi^3\alpha^3$ that closes the gap to $\sim10^{-8}$."

**Recommendation:** Option (i) — the body uses the residual framing; the abstract should match.

**S3b — Case-exhaustion theorem scope:** §VIII Thm 2 reads "any finite composition of projections drawn from the Paper 34 list §III.1–§III.15". Paper 34 v3.2.x catalogues twenty-eight projections. Projections 16–28 (rest-mass, observation/temporal-window, nuclear charge-density, nuclear magnetization, nuclear tensor multipole, PK / core-valence, multipole/Gaunt termination, bipolar harmonic, Young tableau, BO, coupled-channel, gauge choice, Wick rotation, apparatus identity) are NOT covered by the proof body.

**Decision needed:**
- (a) Restrict theorem to "the first fifteen Paper 34 projections (the corpus state as of [date])" plus an explicit note that projections 16–28 are not yet audited.
- (b) Extend the proof body to cover 16–28 (substantive additional work).

**Recommendation:** Option (a) for now; option (b) as a named follow-up sprint.

## 4. MEDIUM errata batch (consolidated)

Eight framing/softening items to batch into a single errata sprint after the four REDs are resolved:

1. **Paper 0 abstract softenings (4 phrases):** "shell capacities 2k−1" → "per-shell angular capacities 2k−1"; "spectrally converge to the hydrogen eigenvalues" → "ground-state eigenvalue spectrally converges to the hydrogen ground-state energy"; "the lattice is the invariant; the continuum geometries are projections" → dual-reading statement; "foundational motif of the framework" → "a foundational motif."
2. **Paper 1 §VI / `companion_alpha` outlook:** Currently promises a "helical photon gauge fiber / symplectic impedance" derivation of α; Paper 2 actually implements the Hopf cubic $\alpha^3 - K\alpha + 1 = 0$. Plus Paper 2 is now in Observations status, not Core. Rewrite §VI Outlook + update the `companion_alpha` bibitem.
3. **Paper 7 framing lags-corpus items:** (i) §III.A "κ alone maps graph eigenvalues to Rydberg" — Paper 18 §κ has $H = \kappa\mathcal{L} + W$. (ii) §I "we demonstrate convergence" — rigorously proven only in Paper 38.
4. **Paper 7 `barut1967` SO(4,2) cite:** Bibitem points to PR 156, 1541 (the SO(4,1) paper) but is used to ground SO(4,2) claims at §I line 67 and §III line 99. PR 157, 1180 (1967) is the SO(4,2) sibling paper. Need PR 157's exact published title verified, then bibitem split or replacement.
5. **Paper 16 §I courtesy cite:** "the periodic law has not been derived from a single mathematical principle" glosses over substantial group-theoretic-periodicity literature (Rumer–Fet SO(4,2)×SU(2); arXiv:2501.18272 SO(4,4)). Add a one-sentence acknowledgment.
6. **Paper 27 γ_∞ figure:** Abstract says "γ_∞ ≈ 1.96 (Richardson on n_max=2,3,4,5)." Actually 1.9589 = n_max=5 local slope; honest Aitken/linear-in-1/n extrapolation gives γ_∞ ≈ 1.93. Soften "Richardson" word + swap number; downstream conclusion (γ_∞ < 2 by multi-shell aggregation) is unchanged.
7. **Paper 32 Sprint L1 "literal identification":** §VIII Rem bisognano_wichmann_reading says the Wick-rotation theorem lifts to "literal identification at the operator-system level." The actual finding ($\sigma_{2\pi}(O) = O$ holds bit-exactly) is a *necessary* condition for BW identification, not full identification with the continuum Wightman vacuum. Soften to "operator-system-level consistency (period closure holds bit-exactly at every Krein cutoff)."
8. **Paper 32 conclusion "four-way coincidence":** Conclusion uses older framing while abstract+§I use deflated "four projections of one triple" (CLAUDE.md WH4). Rewrite conclusion phrasing for consistency.

## 5. LOW cleanup batch

Twelve items across all 11 papers — typos, stale cross-refs, bibitem-key cosmetic mismatches, missing arXiv IDs, orphan bibitems. Full list in the per-paper review files; batch into a final cleanup pass after the errata sprint.

## 6. Next-step decisions for the PI

Three orthogonal threads, each can be picked up independently:

1. **Resolve the four REDs.** R1, R2, R4 are the most directly actionable (R1 and R4 need the right citations identified; R2 is a math fix once you pick option a/b/c). R3 needs the physics judgment on pair-diagonal-m vs total-m.
2. **Approve the Paper 18 §III/§VIII Ω-harmonization** (S1) and the Paper 31 three→five layer rewrite (S2). Once you pick the convention and the scope, the rewrites are mechanical.
3. **Approve the Paper 32 abstract α + case-exhaustion scope rewordings** (S3a, S3b). Yesterday's deferral; the suggested texts are in §3 above.

Once 1–3 are resolved, the MEDIUM errata batch (§4) and LOW cleanup batch (§5) are straightforward and can run as a single sprint.

After Wave 1 is closed (all REDs converted to GREEN), Wave 2 (math.OA arc, Papers 38–49 with REVIEWER added) is the recommended next dispatch. Some of yesterday's and today's discoveries already preview what we'll find there: Paper 38 already had a `perez_sanchez2024` correction-cite misbundle (fixed today corpus-wide); Papers 42, 43, 44 had the Nieuviarts initial bug (fixed today). The Wave 2 audit will surface the math.OA-specific issues those papers carry that we haven't seen yet.
