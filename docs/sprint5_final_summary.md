# Sprint 5 Final Summary — Varshalovich + S⁵ Gauge + Li/Be Fine Structure

**Version:** v2.16.0 (April 16, 2026)
**Status:** Three parallel tracks; one partial, one mixed, one **major positive surprise**

---

## Executive Summary

| Track | Result |
|-------|--------|
| **DV** | Partial. Bipolar channel structure characterized (SS direct = unique (0,2), SS exchange = unique (1,1), SOO has ZERO bipolar channels — confirming it's not pure Y^(1)/r²). Drake M^K basis and bipolar basis are structurally mismatched — no closed-form reconciliation found. Mixing ratios (3/50, -2/5, 3/2, -1) from Sprint 4 BF-D's rational search remain the working result. |
| **S5** | MIXED verdict. U(1) abelian transfers verbatim. SU(3) non-abelian NOT natural. m_l-quotient NOT a CP² spectral discretization. Paper 24 gains §V.D corollary: Coulomb/HO asymmetry is now **three-layered** (spectrum, calibration π, Wilson physical content). |
| **CP** | **MAJOR POSITIVE:** Both Li 2²P (8.89%) and Be 2s2p ³P (2.76% span) closed to <20%. Sprint 3 BF-E's "Tier 3+ honest negatives" were actually a **convention bug** in BR-C that coincided with NIST only for He (Z=2). Fix: standard BP formula + Slater 1930 rules. |

---

## Track DV: Drake Bipolar Harmonics — Honest Partial

**Deliverables:** `debug/dv_drake_bipolar_memo.md`, `debug/dv_step*.py` (12 scripts).

**Positive findings:**
1. Bipolar channel structure for He (1s)(2p) ³P fully characterized in exact sympy:
   - SS (K=2): Direct has UNIQUE (k₁=0, k₂=2); Exchange has UNIQUE (k₁=1, k₂=1)
   - SOO (K=1): BOTH direct and exchange have ZERO allowed bipolar channels — confirming SOO is NOT a pure Y^(1)/r₁₂² tensor (involves p⃗₁ × r̂₁₂ structure)
2. Drake M^K integrals split symbolically into Region I + Region II via production module's `_t_kernel_region_I`
3. Numerical identity confirmed: Bipolar(k₁, k₂, K=2)_direct = M^{k₁}_dir,RegionI + M^{k₂}_dir,RegionII

**Honest negative:** Drake's M^K basis and the bipolar basis are structurally different — no choice of bipolar prefactor β(k₁, k₂, K) reconciles them. Numerical near-hit β=1/√(4π) hints at a Y^K↔C^K convention factor but no exact closed form closes the derivation.

**Recommendation for future work:** Use PSLQ on direct 6D numerical evaluation of ⟨³P_J | H_SS | ³P_J⟩ at multiple Z, identifying rationals in Drake's M^K basis — bypass the bipolar convention chase.

**Paper 14 §V updated** with Sprint 5 DV characterization (Drake J-pattern from Racah 6j, mixing ratios from BF-D rational search with DV structural characterization of why pure-bipolar doesn't close).

---

## Track S5: S⁵ Gauge Structure — MIXED Verdict

**Deliverables:** `debug/s5_bargmann_segal_graph.py`, `debug/s5_edge_laplacian_analysis.py`, `debug/s5_gauge_structure_memo.md`, `debug/data/s5_*.json`.

**Positive findings (matches Paper 24 certification):**
- N_max=5 graph: 56 nodes, 165 edges, π-free (all weights exact sympy Rational)
- Connected (β₀=1); shell degeneracies [1, 3, 6, 10, 15, 21] = (N+1)(N+2)/2
- β₁ = 110 independent cycle classes
- SVD identity between nonzero spectra of L₀ and L₁ verified exactly

**Mixed verdict on gauge structure:**
1. **U(1) abelian transfer: POSITIVE.** The edge Laplacian admits the same discrete Hodge-1 structure as Paper 25's S³ case.
2. **SU(3) non-abelian: NEGATIVE.** The (N,0) SU(3) symmetric reps have raising/lowering matrix elements, but these don't form a natural non-abelian connection in the lattice gauge sense.
3. **CP² spectral identification: NEGATIVE.** The m_l-quotient (12 sectors) gives a Laplacian whose spectrum matches CP²'s Fubini-Study eigenvalues only up to 25% (best linear rescaling). Ratios vary 0.08–0.19 across pairs; no uniform fit.

**New three-layer Coulomb/HO asymmetry (Paper 24 §V.D, NEW):**

| Layer | S³ Coulomb | S⁵ Bargmann |
|-------|------------|-------------|
| Spectrum from L₀ = D−A | Yes (via κ=-1/16) | No (HO spectrum is diagonal) |
| Calibration π | Yes | No (π-free) |
| Wilson physical content | Full (matter+photon) | Reduced (U(1) universal, no matter-propagator role) |

**Paper 25 §VII.1 OPEN question is now ANSWERED** with mixed verdict.

**No natural S⁵ analog of K = π(B+F−Δ):** F-candidate = (1/2)[ζ_R(4) + 3ζ_R(5) + 2ζ_R(6)] — mixed even/odd zetas, no clean π²/6 projection. B and Δ analogs rely on Coulomb projection (restricted by Paper 23 Fock rigidity to S³).

**Paper 24 updated:** new §V.D gauge-structure-reduction corollary (Corollary 4) after §V.C spinor-Lichnerowicz. Paper 25 §VII.1 rewritten from speculative open question to answered mixed-verdict subsection.

---

## Track CP: Li/Be Fine Structure — MAJOR POSITIVE SURPRISE

**Deliverables:** `debug/cp_li_core_polarization.py`, `debug/cp_be_multiconfig.py`, `debug/cp_fine_structure_memo.md`, `debug/data/cp_{li,be}_results.json`, 5 new tests in `tests/test_breit_integrals.py`.

### Results

| System | Sprint 3 BF-E | Sprint 5 CP | Target | Status |
|--------|--------------:|------------:|:------:|:------:|
| Li 2²P splitting | +118% or +553% | **+8.89%** | <20% | **MET** |
| Be 2s2p ³P span (P₀−P₂) | +88.4% | **+2.76%** | <20% | **MET** |
| Be 2s2p ³P P₀−P₁ | — | +18.89% | <20% | MET |
| Be 2s2p ³P P₁−P₂ | — | +6.28% | <20% | MET |

### The actual root cause

**Sprint 3 BF-E's "Tier 3+ honest negatives" were a convention bug, not missing physics.**

BR-C used `ζ = α²·Z_nuc·Z_eff³/[n³l(l+½)(l+1)]` which happens to coincide numerically with the standard formula only for He (where Z_nuc = 2 = 2·Z_val and the explicit (ζ/2) factor in E_SO(J) compensates). For Li (Z_nuc=3) and Be (Z_nuc=4) the factor-of-Z_nuc overcount is UNCOMPENSATED, producing the false 118%/553%/88% errors.

**Corrected physics (standard textbook BP):**
- ζ = (α²/2)·Z_val·⟨1/r³⟩_{2p} where Z_val = asymptotic valence charge (after full shielding)
- Z_eff from Slater's 1930 rules (0.85 per 1s, 0.35 per same-shell electron)
- Li 2p: Z_val=1, Z_eff=1 (full [He] shielding)
- Be 2p in 2s2p: Z_val=1, Z_eff=1.95 (Slater: 4 − 1.70 − 0.35)

### Structural findings (hypothesis NEGATIVES upgraded to positive insight)

1. **Core polarization WORSENS Li, not improves it.** Adding Migdalek-Bylicki V_cp (α_d=0.192, r_c=0.55 a.u., published Table I values) increases the error from +8.9% to +25%. CP is attractive and increases ζ further, while the baseline was already correctly predicting. Sprint 3's hypothesis that CP was the missing physics was INCORRECT.

2. **2-config 2s2p ↔ 2p² mixing is PARITY FORBIDDEN.** 2s2p has parity (−1)^(0+1) = −1 (odd); 2p² has parity (−1)^(1+1) = +1 (even). The 1/r₁₂ operator is parity-even, so ⟨odd|parity-even|even⟩ = 0. Verified symbolically via Slater-Condon — both direct and exchange two-electron integrals vanish. The proper same-parity mixing partner is 2s3p ³P (small ~1% admixture, not needed since convention+Slater fix already closes the gap).

**Paper 14 §V updated:** the Li/Be honest-negative paragraph replaced with Sprint 5 CP closure narrative. Paper 14 fine-structure table gains working rows for all three atoms (He, Li, Be) at <20% error.

**Tests:** 5 new tests in `tests/test_breit_integrals.py` assert the convention fix for Li and Be; 36 existing tests still pass.

---

## Structural Takeaways

1. **The "honest negatives" from earlier sprints have accumulated, and some can be revisited** — Sprint 5 CP showed that Sprint 3 BF-E was a convention bug, not Tier-3+ physics. This suggests a productive pattern: when an honest negative persists, re-examining the convention/setup before adding new physics can be more effective than assuming deeper theory is needed.

2. **The Coulomb/HO asymmetry is now three-layered** (Paper 24 §V.D). Coulomb has spectrum-from-Laplacian + calibration π + full Wilson structure; HO has none of the three. The asymmetry deepens with each sprint.

3. **Drake mixing ratios remain an open puzzle.** Sprint 4 DD closed the J-pattern via Racah 6j; Sprint 5 DV characterized the bipolar channel structure but couldn't close the mixing ratios. The result's practical usability (He 2³P at 0.20%) is unchanged; the structural derivation is partial.

---

## Files Modified (Production)

- `geovac/breit_integrals.py` — docstring updates (DV bipolar references + Drake references)
- `papers/core/paper_14_qubit_encoding.tex` — §V Breit-Pauli paragraph: DV partial + CP closure (Li/Be <20% via convention fix + Slater)
- `papers/core/paper_24_bargmann_segal.tex` — new §V.D gauge-structure-reduction corollary (Corollary 4)
- `papers/synthesis/paper_25_hopf_gauge_structure.tex` — §VII.1 rewritten from open question to answered mixed-verdict
- `tests/test_breit_integrals.py` — 5 new CP tests (Li/Be convention fix + parity verification)
- `CLAUDE.md` — version v2.15.0 → v2.16.0; §2 Sprint 5 bullet

## Files Created (Debug)

- DV: `debug/dv_drake_bipolar_memo.md`, `debug/dv_drake_bipolar.py`, `debug/dv_step*.py` (12 scripts)
- S5: `debug/s5_bargmann_segal_graph.py`, `debug/s5_edge_laplacian_analysis.py`, `debug/s5_gauge_structure_memo.md`, `debug/data/s5_{bargmann_graph,graph_spectrum}.json`
- CP: `debug/cp_li_core_polarization.py`, `debug/cp_be_multiconfig.py`, `debug/cp_fine_structure_memo.md`, `debug/data/cp_{li,be}_results.json`
- Plans/summaries: `docs/sprint5_tier4_plan.md`, `docs/sprint5_final_summary.md` (this file)

---

## Sprint 6 Candidates

Surfaced by Sprint 5:

1. **DV follow-up via PSLQ** — direct 6D numerical evaluation of ⟨³P_J|H_SS|³P_J⟩ at multiple Z, PSLQ-identify rationals in Drake's M^K basis. Bypasses the bipolar convention chase.

2. **Re-examine prior honest negatives for convention bugs** — CP showed this pattern works. Candidates: the α combination rule's B/F/Δ could have been mis-normalized somewhere; graph validity boundary at Z_c might have a setup issue.

3. **Paper 18 taxonomy consolidation** (carried over from Sprint 4): distributional embedding + log-embedding + Racah-intrinsic vs tensor-dependent.

4. **[Rn] frozen core for RaH** — direct Sunaga matched comparison.

5. **Tier-3 T7 Kramers-Pasternak** — closed-form Dirac-Coulomb ⟨r^{-2}⟩, ⟨r^{-3}⟩ for n_r ≥ 1.
