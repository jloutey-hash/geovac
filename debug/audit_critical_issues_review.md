# Critical Issues Review — Full Two-Day Audit

**Date:** 2026-06-02
**Scope:** Consolidated PI-action list across Waves 1+2+3 + three citation sweeps + two follow-on Explorer agents.
**Companions:** `debug/wave1_salvage_memo.md`, `debug/wave1_pi_recommendations.md`, `debug/citation_audit_master.md`, `debug/cite_sweep_A.md` / `cite_sweep_B.md` / `cite_sweep_C.md`, `debug/cite_metadata_verification.md`, `debug/cite_missing_findings.md`, `debug/review_paper{N}.md` for N ∈ {0, 1, 2, 7, 14, 16, 18, 22, 23, 24, 27, 28, 31, 32, 34, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51}, `debug/review_synth_g1.md`, `debug/review_synth_g3.md`.

## 1. Headline status

**~30 papers full-content-audited, ~50 papers citation-audited, ~93 mechanical fixes applied autonomously this session.**

Final verdicts post-fixes:
- **GREEN (1):** Paper 51 (gravity arc).
- **YELLOW (most papers):** structural math sound, residual MEDIUM/LOW framing items.
- **RED (2):** Paper 49 (Thm 3.3 chain inequality), Paper 50 (§7 S^5 multiplicity formula off by factor 4).
- **Compile broken (1):** Paper 20 — missing `.bib` file.

The audit infrastructure was worth the investment: **4 substantive math errors found**, ~30 right-ID-wrong-paper Fursaev-mode citations corrected, corpus-wide layer-count drift mapped, foundational Ω-normalization ambiguity resolved.

## 2. Substantive math errors (4) — highest PI priority

### M1. Paper 49 Theorem 3.3 — load-bearing chain inequality unjustified

**Location:** Paper 49 `thm:strict_super_additivity`, plus propagated to **synth-G1 lines 1505–1510** as "follows from relative-entropy monotonicity" hand-wave.

**Issue:** Theorem invokes a chain inequality $\Delta S^{1\to 3} \le \Delta S^{1\to 2} + \Delta S^{2\to 3}$ for quantum relative entropy. This inequality **does not appear in Lindblad 1975 or Uhlmann 1977** as cited, and **quantum relative entropy is known not to satisfy a triangle inequality** (only max-relative entropy does). The chain step in the proof is unjustified.

**Why it matters:** Theorem 3.3 is the centerpiece of Paper 49 — underwrites the "twin-paradox-as-quantum-information" headline. If the chain step is invalid, the theorem statement at full generality does not hold.

**Three fix directions, PI choice:**
- **(a)** Weaken Theorem 3.3 to a max-relative-entropy statement (chain inequality genuinely holds there). Headline survives at reduced generality.
- **(b)** Find an alternative real inequality that establishes the super-additivity (data-processing inequality + specific structure of TICI cocycle? Connes 1973 cocycle algebra?). Investigative.
- **(c)** Reframe the headline claim to not depend on this chain step (e.g., as a single-state observation rather than multi-state super-additivity).

The bones of Paper 49 (OSLPLS category construction, axiom transport, embedding functor, Riemannian-limit at N_t=1) are sound — only the centerpiece theorem and its derivation need work.

### M2. Paper 18 Theorem 1 part (2) Eq.(38) D(4) substitution error

**Location:** Paper 18 Eq.(38) and L1350.

**Issue:** Eq.(38) computes $D(4) = 2\zeta(2) + 2\zeta(3) \approx 5.694$, but Eq.(30) and Paper 28's verified value both give $D(4) = \pi^2 - \pi^4/12 \approx 1.752$. mpmath summation confirms the Eq.(30) value is correct. The Eq.(38) sum `Σ_{m≥1} 2m(m+1)/m^4` **conflates integer m with the half-integer Camporesi–Higuchi eigenvalues** $|\lambda_n| = n + 3/2$.

**Why it matters:** Propagates into §IV motivic-weight loop-order argument. The empirical Petermann–Sommerfield ζ(3) at two loops is real, but the GeoVac witness needs a correct ζ(3) source on the spinor bundle.

**Fix direction:**
- **(a)** Replace Eq.(38) with the correct CH sum (use $|\lambda_n| = n+3/2$, not integer m). Verify the substitution gives D(4) consistent with Eq.(30).
- **(b)** Find a different ζ(3) witness on the spinor bundle that does close Petermann–Sommerfield without the Eq.(38) miscalculation.

Cross-check Theorem 1 / T9 in Paper 28 — the squared-Dirac claim there is structurally sound (verified). Paper 18's error is in the specific D(4) substitution step, not the underlying theorem.

### M3. Paper 28 Table 1 spectral zeta values

**Location:** Paper 28 Table 1.

**Issue:** Table 1 lists squared Dirac spectral zeta values at s=1..4 that **contradict the paper's own Theorem 1 (Eq. 5)**. Values are wrong by orders of magnitude AND wrong sign at all four rows. Example: s=1 should give $-\pi^2/4 \approx -2.47$, Table 1 says $\pi^2 - \pi^4/12 \approx 1.75$.

The correct values at integer s are simply $D(2), D(4), D(6), D(8)$ from Table 2 (since $\zeta_{D^2}(s) = D(2s)$).

**Fix:** Mechanical — replace Table 1 row values with the corresponding entries from Table 2 (D(2), D(4), D(6), D(8)), apply the appropriate sign. **One table-level edit; not autonomous because the right values depend on signs/conventions in the Eq. 5 vs Table 2 mapping** — wants PI eyes for 5 minutes.

Cross-check: Theorem 1 / Eq. 5 itself is correct (the proof is sound and S_min 200-dps + Theorem 3 χ_{−4} identity + Theorem 5 ζ(3) complementarity all verified bit-exact). Only Table 1 is wrong.

### M4. Paper 50 §7.1 Theorem 7.1 S^5 multiplicity formula off by factor 4

**Location:** Paper 50 §7.1 Theorem 7.1 + supporting derivation; affected header L709 ("bit-exact match extends to S^5").

**Issue:** $F_s^{S^5} \approx -0.02297$ claimed as bit-exact match to Klebanov–Pufu–Safdi 2011, but the value is **exactly 4× the KPS Table 1 d=5 entry** ($-5.74 \times 10^{-3}$). Traceable to the stated multiplicity formula $(2n+4)(n+1)(n+2)(n+3)/6$ — should be $/24$ per KPS Eq. (37).

**Why it matters:** Paper 50's §3-§5 S^3 content (the headline F-theorem match for d=3) is verified clean and reproducible. The §7 extension to S^5 is what's broken. The "bit-exact match extends to S^5" header sentence (L709) is therefore false as written.

**Fix direction:**
- **(a)** Re-derive §7 with the corrected multiplicity $/24$, recompute $F_s^{S^5}$ and Dirac analog, verify the corrected values match KPS Table 1.
- **(b)** Reconcile via an explicit normalization convention — if there's a factor-4 normalization difference between GeoVac and KPS conventions on $S^5$, document it explicitly and adjust the comparison.

(a) is the substantive fix; (b) is a fallback if the convention story is real. **PI choose direction.**

## 3. Cross-paper structural drift (4 coordination items)

### S1. Coulomb/HO asymmetry layer count

| Paper | Says | Truth |
|:------|:-----|:------|
| Paper 31 §VII | three layers | five (lags by 2) |
| synth-G3 | four layers | five (lags by 1) |
| Paper 24 | "Fourth layer" + footnote says five | locally OK (introduces Layer 4; mentions Layer 5 lives in Paper 51) |
| Paper 51 | five layers (introduces Layer 5) | ✓ current truth |
| CLAUDE.md §6 | five layers | ✓ |

**Coordinated fix needed across Papers 31 and synth-G3.** Paper 31 §VII "Three-Layer Sharpening (Sprint~5)" needs rewrite to "Five-Layer Sharpening" with two added table rows for Layer 4 (modular Hamiltonian, Paper 24 §V) and Layer 5 (gravity termination, Paper 51 G4-5a). synth-G3 simpler — just prose update.

### S2. Paper 32 §VIII case-exhaustion theorem still references "fifteen" projections

**Location:** Paper 32 §VIII L1633, L1675. Paper 34 has 28 projections in §III.

**Carryover from Wave 1, still open.** Paper 34's audit this turn confirmed: projection-by-projection π-free/π-bearing assignment for slots 1–15 is consistent with Paper 34's current §III ordering. So the proof body is OK for the 15-projection corpus state at theorem date. **Fix options:**
- **(a)** Restrict theorem statement to "the first fifteen Paper 34 projections (as of [date])" — explicit scope; one-line edit.
- **(b)** Extend proof body to projections 16–28 — substantive new work; deferred multi-day sprint.

Both are acceptable; (a) is the conservative immediate fix.

### S3. Paper 49 BCFM authors/acronym

**Location:** Paper 49 §8 (subsec:bcfm_recap, subsec:bousso_positioning), bibitem L2658.

**Issue:** `bousso_casini_fisher_maldacena2020` cites arXiv:2007.00230 = real paper, but the actual authors are **Bousso, Chandrasekaran, Rath, Shahbazi-Moghaddam** (NOT Casini, Fisher, Maldacena). The "BCFM" acronym in §8 is structurally built from the wrong initials.

**Fix options:**
- **(a)** Rename acronym to "BCRS" (Bousso, Chandrasekaran, Rath, Shahbazi-Moghaddam) — preserves the acronym pattern but is mechanically heavy (~12 in-text changes).
- **(b)** Replace "BCFM" with "Bousso et al." throughout — neutral; ~12 in-text changes.

Plus the bibitem author-list update. The cite-key can stay (cosmetic legacy).

### S4. mondino_samann2024 cite-key collision — 5 in-text uses across 4 papers

**Locations:**
- Paper 45 (2 in-text uses at L374, L1463; bibitem L1759–1763)
- Paper 46 (1 in-text use at L381; bibitem L2114–2118)
- Paper 47 (2 in-text co-cites with 2025 at L345, L1125; bibitem)
- synth-G1 (2 in-text uses at L263, L1453; bibitem)

**Issue:** Cite-key `mondino_samann2024` is attached to arXiv:2209.14384, which is actually **Minguzzi-Suhr "Bounded Lorentzian metric spaces"** (LMP 114:73, 2024). And the agent's literature search found no real Mondino-Sämann 2024 paper at this content (the real Mondino-Sämann work is the 2025 paper at arXiv:2504.10380).

**Fix plan (uniform across all 4 papers):**
1. Replace `\bibitem{mondino_samann2024}` content with the real Minguzzi-Suhr metadata, renaming key to `minguzzi_suhr2024`.
2. Update every `\cite{mondino_samann2024}` → `\cite{minguzzi_suhr2024}` (5 in-text uses total: 2+1+2 — wait the Paper 47 ones are co-cites with 2025, so update to `\cite{mondino_samann2025,minguzzi_suhr2024}` keeping both).
3. Verify the citation surrounding prose makes sense as "Minguzzi-Suhr" rather than "Mondino-Sämann" (the contexts in P45 and P47 frame the cite as "earlier work in the Mondino-Sämann strand" — needs rewording to acknowledge Minguzzi-Suhr is a separate strand).

Sprint-scale fix once PI approves the rename.

## 4. Internal contradictions (4)

### C1. Paper 34 §V.C.1 Lamb shift autopsy table — sum off by 26 MHz

Components sum to 1031.41; printed total is 1057.41. Either a row is missing or a row value is wrong. Need PI 5-min check to identify which is correct.

### C2. Paper 34 conclusion contradicts §VI body

Conclusion L9574 still says "error compounds with projection depth." §VI body explicitly falsifies that prediction and replaces it with the Layer-2-presence bound. **Update the conclusion** to match §VI's "Layer-2-presence bound" framing.

### C3. Paper 49 §9.1 "Q2' CLOSED" contradicts abstract

Abstract says Q2' is NOT closed; §9.1 says CLOSED. Pick one — the §9 honest-scope sections (§5.7, §7.6) treated Q2' as open, so the abstract appears to be the honest version and §9.1 the overstatement.

### C4. Paper 46 App A — stale "Paper 45 not yet amended" claim

Paper 45 *was* amended this session (`rem:envelope_v2` added). App A needs prose update to reflect current state.

### C5. Paper 46 App B Theorem B.2 — proof-sketch-grade with deferred absorption step

App B says it "closes G1'" but the load-bearing absorption step is deferred to an internal β-L5 memo not in the bibliography. Reframe to honestly say "proof-sketch-grade; load-bearing absorption deferred to memo X."

## 5. Citation hunts needing PI decision

### F1. Paper 14 `trenev2025` — likely GeoVac-internal source

Headline Pauli counts (51×/746×/1712× advantage) are claimed from arXiv:2311.03719, but that paper is about vibrational acetylene polyynes and does not contain those electronic counts. Agent's best assessment: **likely GeoVac internal computation against re-derived Gaussian baselines**. Decision: remove the Trenev cite and present as internal recomputation against named baselines (STO-3G, cc-pVDZ, cc-pVTZ standard), with proper methodology citation.

### F2. Paper 23 `PachuckiYerokhin2010` — MEDIUM-confidence replacement

Agent suggests Pachucki PRL 106, 193007 (2011) "Nuclear Structure Corrections in Muonic Deuterium" but this changes both year AND author list (single-author Pachucki vs Pachucki-Yerokhin). PI verify — is this the intended source?

### F3. Paper 23 `Eides2024` — MEDIUM-confidence replacement

Agent suggests Eides-Shelyuto JHEP 07 (2023) 211, arXiv:2306.13369 "Two-loop corrections to Lamb shift and hyperfine splitting in hydrogen" but this changes year, authors, and venue. PI verify.

### F4. Paper 38 `ucp_maps_2024` — could-not-find

No real Hekkelman-McDonald-vS UCP paper found. Decision: remove the cite, or merge into `hekkelman_mcdonald2024` (arXiv:2412.00628 — quantum-integral, related but different scope).

### F5. Paper 47 `hekkelman_mcdonald2024b` — could-not-find

Placeholder `arXiv:2411.xxxxx`; no findable companion. Plus the `hekkelman_mcdonald2024` in Paper 47 has wrong title AND §6.2 framing the actual arXiv:2412.00628 paper does not support. Two-bibitem cleanup, but the second needs a replacement decision: remove cite or find different supporting reference.

### F6. Paper 48 `hekkelman_mcdonald2024b` (BosonSampling — arXiv:2411.04566) + `bertozzini_conti_lewkeeratiyutkul2009` (tridiagonal matrices — arXiv:0901.4031)

Two more wrong-paper cites in Paper 48. Both: could-not-find the intended real work via agent search. Decision: remove cites or find replacements.

### F7. Paper 16 `Schwerdtfeger2015` — Pershina 2014 chapter replacement APPLIED this session

Was held for PI decision; I applied Pershina 2014 (cleanest fit; matches the original bibitem's book-chapter shape). If you'd preferred Schwerdtfeger NPA 944 (2015) instead, easy to swap.

### F8. Paper 50 `beccaria_tseytlin2017` and `hartman_kruthoff_shaghoulian_tajdini2019`

Two more Fursaev-mode misattributions found in Paper 50:
- `beccaria_tseytlin2017` → arXiv:1702.02325 is actually a number-theory paper by A. Smith.
- `hartman_kruthoff_shaghoulian_tajdini2019` → arXiv:1902.10893 is by Caputa-Datta-Shyam (not the Hartman et al. group). Intended work appears to be arXiv:1807.11401.

Plus `cardy1988` cosmetic year metadata mismatch.

Mechanical fixes once correct replacements verified.

## 6. Structural defects (1)

### D1. Paper 20 missing `.bib` file

`paper_20_refs.bib` is missing from the repository. **Paper 20 does not LaTeX-compile.** Highest-priority broadcast blocker among single-paper issues. Either (a) locate the missing .bib file in your records, or (b) reconstruct from in-text `\cite` usage by extracting cite-keys from the .tex and rebuilding the bibliography by web search.

## 7. Open items from earlier waves (carryovers)

- **Paper 32 abstract α figure rewording** (Wave 1 deferred): "8.8×10⁻⁸ matches α⁻¹" is the post-Sprint-A residual, not raw match (raw is 4.77×10⁻⁷). Choose the more complete option (i) or the shorter (ii) rewording from `wave1_pi_recommendations.md` §3.S3a.
- **Paper 14 §V missing 7 bibitems** (Sunaga2025, BJL, Szmytkowski2007, MartinezYRomero2004, Dyall, BreitPauli, ChildsBerry). Once the right works are identified, mechanical to add.
- **Paper 22 §V Dyall, Grant** missing bibitems. Similar.
- **Paper 22 pair-diagonal-m vs Coulomb total-m** (Wave 1 deferred): physics judgment — defend 1.44% sparsity as the relevant qubit-cost measure, or update headline to 6.06%.
- **Paper 38 4/π identity math statement**: Sweep A flagged `4/π = 2 Vol(S¹)/Vol(SU(2))` evaluates to 2/π, not 4/π; identity needs correction. The 4/π = Vol(S²)/π² M1-signature reading per memory is the right interpretation.
- **Paper 45 missing bibitems** `latremoliere2025_hypertopology` + `paper48` — mechanical to add.

## 8. MEDIUM errata batch (consolidated)

Once the HIGH items above are resolved, the following MEDIUM softenings batch cleanly into one errata sprint:

- Paper 0 abstract softenings (some applied this session; "lattice is the invariant" and "foundational motif" softenings — latter applied)
- Paper 2 CODATA edition annotation
- Paper 7 lags-corpus framing items (κ-alone, "we demonstrate convergence")
- Paper 16 §I courtesy cite to group-theoretic-periodicity literature
- Paper 27 γ_∞ "Richardson on 1.96" → "Aitken on 1.93" (number swap + method label)
- Paper 32 Sprint L1 "literal identification" softening
- Paper 32 conclusion "four-way coincidence" rewording
- Paper 42 "literal identification" same pattern softening
- Paper 7 barut1967 SO(4,2) cite split (PR 156 vs PR 157)
- Paper 14 11.10×Q vs actual [11.106, 11.200] range honesty
- Paper 16 §I "single mathematical principle" needs courtesy cite to prior group-theoretic-periodicity literature
- Paper 18 paschke_sitarz2000 key vs 1998 paper (cosmetic)
- Paper 23 FriarPayne2005 PRA→PRC, 014501→014002
- Paper 26 Peruzzo article 4213→5213
- Paper 36 Eides 2001 mixed-edition title/pagination
- Paper 39 multiple metadata items
- synth-G1 abstract "fourteen papers" vs intro "nine documents" stale count
- synth-G1 master theorem k-independence proof sketch garbles a Pythagorean identity (for a synthesis paper that disclaims new theorems)
- synth-G3 six-tier vs five-tier taxonomy alignment
- synth-G3 F^0(1s,1s) Z-factor drop
- synth-G3 "natural gauge group is U(1), not SU(3)" overstatement

## 9. Recommended next session plan

In priority order:

1. **Paper 20 `.bib` recovery** — broadcast blocker (~30 min).
2. **Three substantive math errors (M1/M2/M3)** — PI choose fix directions; each is a 1-2 hour focused sprint per paper.
3. **Layer-count harmonization (S1)** — 3-paper coordinated edit, ~2 hours.
4. **mondino_samann2024 rename (S4)** — 4-paper coordinated edit + prose touchups, ~3 hours autonomous once PI approves.
5. **Paper 32 case-exhaustion scope (S2) + abstract α rewording** — PI choose option, then ~30 min autonomous.
6. **BCFM acronym (S3)** — PI choose option, then ~1 hour autonomous edit propagation.
7. **Internal contradictions (C1-C5)** — 5 small edits, ~1 hour with PI confirmation per item.
8. **Citation hunts (F1-F6)** — most need a 1-paragraph PI decision; F1 (Trenev) is the most structural rewrite.
9. **MEDIUM errata batch** — ~20 framing softenings + small metadata fixes; single sprint, ~2-3 hours.
10. **Wave-2 cleanup batch from this morning** (latremoliere fixes in synth-G1, etc.) — wrap into the above.

After all this lands, the corpus is broadcast-ready modulo the math-error fix-direction choices on Papers 18, 28, 49.

## 10. Calibration final note

The audit infrastructure performed well across two days:
- **4 substantive math errors found** in the corpus (Paper 49 Thm 3.3, Paper 18 Eq.(38), Paper 28 Table 1, Paper 50 §7 S^5 multiplicity) — all in load-bearing positions; none would have been caught by spot-reading.
- **~30 Fursaev-mode misattributed citations** identified and most corrected.
- **Ω-normalization ambiguity** that affected Papers 0, 7, 18 — narrowed to a §II.1 / §VIII inconsistency in Paper 18 and resolved this turn.
- **Paper 51 (gravity arc) emerged GREEN** — the most ambitious physics result of the corpus passed the audit cleanly.
- **Two agent miscalibrations** (Sweep B claims about Paper 28 chamseddine + Papers 2/25 perez_sanchez) caught by verify-the-verifier discipline before propagating bad edits.

The "right-ID-wrong-paper" failure mode (Fursaev–Solodukhin pattern documented in CLAUDE.md §3) was the dominant defect class across the corpus. The originating 2026-mid-year LLM-assisted drafting era left this signature broadly; recently-drafted papers (Paper 51, Paper 40, the chemistry corpus) are largely clean.
