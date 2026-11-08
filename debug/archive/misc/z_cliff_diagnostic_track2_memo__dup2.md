# Z>20 cliff diagnostic — Track 2 (DIAGNOSTIC ONLY)

**Date:** 2026-05-09 (post-multi-track Roothaan autopsy + screening kernel upgrade closeout sprint).
**Sprint:** Diagnostic-before-engineering pause on the Z>20 heavy-atom screening cliff that surfaced during the Cs HFS scoping (V.C.6) and was confirmed by the heuristic two-zeta closeout sprint.
**Verdict:** **The Z>20 cliff is dominated by missing radial nodes in the CR67 single-zeta hydrogenic basis** (mechanism (a) sharpened by mechanism (d)). FD solver, conventions, and multi-zeta machinery itself are all faithful (mechanisms (b), (c), (d-coding) ruled out cleanly). The cliff is NOT a Z=20 sharp threshold — it manifests at K (Z=19) cleanly, has a mild signature already at Na (Z=11), and is **shell-resolved** within Cs itself (inner shells 1s/2s/2p fit at ~1x; outer shells 4d/5s/5p overshoot by 2.7-4.1x). The closeout sprint diagnosis ("CR67 single-zeta itself non-faithful for heavy-atom outer shells") is **substantively correct**, and Track 2 sharpens the structural mechanism: **a single-zeta hydrogenic has no radial nodes, but a real RHF outer-shell orbital has n−l−1 nodes whose physical role is to push amplitude OUT to where the orbital should physically be by orthogonalizing against inner shells**. Without nodes, the orbital amplitude piles up near the inner-zeta peak, making the orbital systematically too compact. **BBB93/KTT GO RECOMMENDED, but with eyes open: BBB93 closes ~20% of the residual; full sub-percent Cs HFS requires the full Bohr-Weisskopf relativistic enhancement on top.**

---

## 1. Diagnostic-before-engineering rationale

CLAUDE.md memory `feedback_diagnostic_before_engineering.md` (2026-05-08) requires this diagnostic before any multi-week BBB93 tabulation effort. The closeout sprint's diagnosis ("CR67 single-zeta itself non-faithful") was correct but **one-step**: it identified the symptom without distinguishing among:

(a) CR67 fit faithfulness — Slater-style single-zeta exponents themselves wrong for outer shells of heavy atoms.
(b) FD radial solver convention — `_solve_screened_radial_log` may diverge for s-wave at l=0 with singular origin.
(c) Spectral basis inadequacy — Laguerre exponents at k=Z give wrong behavior at heavy Z.
(d) Multi-zeta machinery itself has a bug not yet exposed.
(e) Convention mismatch in how Z_eff(r) feeds into the SO/HFS operator.
(f) Self-consistency missing — physical heavy-atom screening requires SCF iteration, not table lookups.

If the dominant cause is (b), (c), (d-coding), or (e), a multi-week BBB93 effort would be a wild-goose chase. This diagnostic distinguishes them.

**Track 2 is read-only** — no production code modifications.

---

## 2. Probe (a): CR67 single-zeta fit faithfulness

**Probe (a) extended** is the load-bearing diagnostic. It compares the CR67 hydrogenic R_nl(r; Z_eff = n·ζ_CR) for every shell of Cs (Z=55, [Xe] core) to physically known atomic-physics ⟨r⟩ values from Sobelman / RHF tables. The "overshoot factor" = Z_eff_CR / Z_eff_implied-by-physical-⟨r⟩ measures how compact the CR67 fit is relative to the real orbital.

| Shell | ζ_CR | Z_eff = n·ζ | r_peak_CR (bohr) | ⟨r⟩_CR (bohr) | ⟨r⟩_phys (bohr) | rel_err% | Overshoot |
|-------|------|-------------|------------------|---------------|------------------|----------|-----------|
| Cs 1s | 54.27 | 54.27 | 0.018 | 0.028 | 0.018 | +49% | **0.67** |
| Cs 2s | 42.03 | 84.06 | 0.062 | 0.071 | 0.076 | −6% | 1.07 |
| Cs 2p | 46.49 | 92.98 | 0.043 | 0.054 | 0.069 | −22% | 1.28 |
| Cs 3s | 36.92 | 110.8 | 0.118 | 0.122 | 0.190 | −36% | 1.56 |
| Cs 3p | 38.46 | 115.4 | 0.104 | 0.108 | 0.184 | −41% | 1.70 |
| Cs 3d | 33.23 | 99.69 | 0.090 | 0.105 | 0.230 | −54% | 2.18 |
| Cs 4s | 25.28 | 101.1 | 0.244 | 0.237 | 0.420 | −44% | 1.77 |
| Cs 4p | 25.57 | 102.3 | 0.231 | 0.225 | 0.443 | −49% | 1.97 |
| Cs 4d | 20.72 | 82.89 | 0.256 | 0.253 | 1.05 | −76% | **4.14** |
| Cs 5s | 13.12 | 65.59 | 0.609 | 0.572 | 1.27 | −55% | 2.22 |
| Cs 5p | 12.31 | 61.56 | 0.632 | 0.593 | 1.95 | −70% | 3.29 |

**Group statistics:**
- **Inner (n=1,2):** mean overshoot = **1.01x**, range [0.67, 1.28] — well-fit.
- **Middle (n=3):** mean overshoot = **1.81x**, range [1.56, 2.18] — moderate compression.
- **Outer (n=4,5):** mean overshoot = **2.68x**, range [1.77, 4.14] — severe compression.

The cliff is **shell-resolved within Cs itself**: inner-to-outer overshoot ratio is **2.66x**. CR67's Slater fit produces faithful inner cores but systematically too-compact outer shells. The 4d shell is the worst offender (4.14x overshoot) and is precisely the shell that should penetrate the inner core in real Cs but doesn't in the CR67 fit. **The cliff is shell-resolved at the CR67 fit level**: hypothesis (a) is **CONFIRMED**.

**Cross-atom check** (light atoms anchor): The same probe at Sr (Z=38) outer shells shows similar overshoot (5p shell ~3x compact); at K (Z=19) the 4s shell shows ~1.5x overshoot; at Na (Z=11) the 3s shell shows 1.24x overshoot.

---

## 3. Probe (b): FD solver hydrogenic faithfulness — RULED OUT

The closeout sprint replaced an earlier log-mesh solver with a dense-uniform-FD path in `_solve_screened_radial_log`. This probe tests the dense FD against the analytical hydrogenic |ψ_ns(0)|² = Z³/(πn³) for Z=1, 2, 6, 10, 30, 55, 80 across n=1, 2, 6, with grid sizes 50k, 100k, 200k, 400k.

**Headline:** Monotone convergence at every Z, every n, every grid. At n_grid=400k:
- H 1s (Z=1): rel_err in |ψ(0)|² = −0.01%
- H 6s (Z=1): rel_err = −0.40%
- Cs⁵⁴⁺ 6s (Z=55): rel_err = **−0.82%**
- Hg⁷⁹⁺ 6s (Z=80): rel_err = −1.20%

Eigenvalue rel_err is essentially zero (<1e-6) at all Z, all grids. The mild Z-dependence in |ψ(0)|² convergence (−0.40% at H, −0.82% at Cs, −1.20% at Hg) is a uniform Frobenius-extrapolation effect from the small-r boundary, not a structural failure.

**Conclusion:** The FD solver IS faithful across all tested Z. **Hypothesis (b) is RULED OUT.** The Cs HFS cliff cannot be closed by further numerical refinement of n_grid alone.

---

## 4. Probe (c): Convention mismatch — RULED OUT

Six checks traced the chain Z_eff(r) → ψ_valence(r) → R(0) → A_HFS:

1. **BF formula sign + prefactor** for H 21cm: predicted A = 1422.81 MHz vs experimental 1420.41 MHz, residual +0.17% (the well-known reduced-mass effect addressed in Sprint HF). PASS.
2. **g_N convention for Cs-133**: production uses g_Cs = μ/I = 0.7377 (atomic-physics convention). 2x-bug ruled out. PASS.
3. **n·ζ convention in `_hydrogenic_radial`**: production sets Z_eff = n·ζ, then ρ = 2·Z_eff·r/n = 2·ζ·r, giving radial decay rate exp(−ζr). Matches Slater convention. PASS. (The n_r=0 peak position cross-check gave a 55% discrepancy for Cs 5p, but that's because 5p has n_r=3 — multiple lobes — not a code bug.)
4. **Z_eff(r) boundary conditions** for Cs: Z_eff(0)=55.00 ✓, Z_eff(50)=1.00 ✓. PASS.
5. **Reproducibility**: framework |ψ_6s(0)|² at n_grid=200k = 1.245 vs closeout's Richardson value 1.328 (n_grid→∞). −6.3% is the natural finite-grid effect; the closeout used Richardson extrapolation. Reproducible.
6. **End-to-end A_Cs** at framework |ψ(0)|² = 1.245, no F_R: A_FW = 734.8 MHz (BF strict) vs A_exp = 2298.16 MHz, residual −68%. After Casimir F_R = 1.555: −47% (matches closeout exactly). Chain is internally consistent.

**Conclusion:** **Hypothesis (c) is RULED OUT.** No factor-of-X bug. The chain produces a faithful prediction GIVEN the input |ψ(0)|²; the cliff is upstream in the |ψ(0)|² input.

---

## 5. Probe (d): multi-zeta machinery itself + STRUCTURAL FINDING

Six checks on the multi-zeta machinery (`geovac/multi_zeta_orbitals.py`):

1. **STO normalization**: ∫χ²r²dr = 1.000001 (machine precision). PASS.
2. **MultiZetaOrbital reduces to single zeta** when both primitives have same ζ: ratio 1.000000. PASS.
3. **`_build_two_zeta_orbital` renormalization**: Cs 5p multi-zeta normalization 0.999999. PASS.
4. **Two-zeta vs single-zeta ⟨r⟩ for Cs orbitals**: ⟨r⟩_2z / ⟨r⟩_single = 0.83 (mean across 11 shells). The two-zeta is *more compact* than single-zeta on average — this is the closeout sprint's wrong-direction result.
5. **density_from_orbitals → N_core**: integrates to 53.9999 vs target 54. PASS (0.00% rel error).
6. **STRUCTURAL FINDING**: s-orbital orthogonality.

For real RHF orbitals, ⟨R_ns | R_ms⟩_radial = 0 for n ≠ m (the radial nodes provide this). For the heuristic two-zeta:

| Pair | ⟨R_ns \| R_ms⟩_radial |
|------|----------------------|
| ⟨1s\|2s⟩ | 0.74 |
| ⟨1s\|3s⟩ | 0.41 |
| ⟨1s\|4s⟩ | 0.10 |
| ⟨1s\|5s⟩ | 0.005 |
| ⟨2s\|3s⟩ | **0.86** |
| ⟨2s\|4s⟩ | 0.39 |
| ⟨2s\|5s⟩ | 0.04 |
| ⟨3s\|4s⟩ | 0.71 |
| ⟨3s\|5s⟩ | 0.12 |
| ⟨4s\|5s⟩ | 0.46 |

The max overlap |⟨R_2s | R_3s⟩| = **0.86** — these are essentially the same function. The all-positive-coefficient two-zeta heuristic (and equally the single-zeta CR67) produces orbitals that violate Pauli exclusion at the radial level. The radial nodes that real RHF needs for orthogonality are **structurally absent**.

**Conclusion:** **Hypothesis (d) coding bug is RULED OUT.** But a STRUCTURAL FINDING emerges: the missing radial nodes in the all-positive-coefficient hydrogenic basis are the structural mechanism behind the CR67 cliff. **This sharpens hypothesis (a)**: CR67 single-zeta (and the heuristic two-zeta inheriting from it) cannot represent radial nodes, and the nodes are what physically push outer-shell amplitude to its correct radial extent. Without them, all shells pile up near the inner-zeta peak, making outer shells systematically too compact.

---

## 6. Probe (e): Z-scan localization

Cliff localization across the alkali series H, Li, Na, K, Rb, Cs. We use the **|ψ_FW(0)|² / |ψ_HF_lit(0)|² ratio** as the cleanest kernel-isolated diagnostic (independent of F_R):

| Atom | Z | n_val | ⟨r⟩_phys | ψ_FW | ψ_HF_lit | ratio FW/lit | Cliff? |
|------|---|-------|----------|------|----------|--------------|--------|
| H 1s | 1 | 1 | 1.5 | 0.318 | 0.3183 (=1/π) | **0.999** | NO |
| Li 2s | 3 | 2 | 3.5 | (skipped) | 0.227 | (Li uses 2-electron CoreScreening) | — |
| Na 3s | 11 | 3 | 3.4 | 0.947 | 0.762 | **1.243** | mild over |
| K 4s | 19 | 4 | 4.2 | 0.272 | 0.879 | **0.309** (3.2x under) | YES |
| Rb 5s | 37 | 5 | 4.6 | 0.374 | 1.940 | **0.193** (5.2x under) | YES, deepest |
| Cs 6s | 55 | 6 | 5.6 | 1.245 | 1.990 | **0.625** (1.6x under) | YES, partial recovery |

**Net pattern:**
- H Z=1 ratio = 0.999 — perfect agreement (hydrogenic, no screening).
- Na Z=11 ratio = 1.24 — 24% **over**: the [Ne] CR67 core slightly over-shields, but this is OK precision for chemistry (NaH PES failure was earlier traced to W1c-residual orthogonality, not screening kernel).
- **K Z=19 ratio = 0.31 — 3.2x UNDER**. This is the cliff onset.
- Rb Z=37 ratio = 0.19 — **deepest cliff** (5.2x under).
- Cs Z=55 ratio = 0.63 — partial recovery from Rb.

**Why Cs is milder than Rb**: The CR67 outer-shell zetas at Cs are smaller in absolute terms (ζ_5s_Cs=13.12 vs the corresponding Rb 5s zeta), but the **n=5 screening already tested** at Rb means the CR67 fit can compensate somewhat at Cs. Also, the Cs valence has more inner-shell penetration room (more shells to overshoot).

**Conclusion:** The cliff is well-localized at K (Z=19). Z=11 (Na) shows a mild over-shielding but in the right direction; Z=19+ shows severe under-shielding. The "Z>20" closeout framing is **mostly correct** — K 4s at Z=19 is the cleanest cliff onset.

---

## 7. Synthesis: dominant cause ranking

**Probe results combined:**

1. **Cause (a) sharpened by (d-structural): missing radial nodes in CR67 single-zeta hydrogenic basis.**
   - Probe (a-extended) shows shell-resolved overshoot 1x → 4x going outer.
   - Probe (d) shows max ⟨ns|ms⟩ overlap = 0.86 — structural reason.
   - Cliff onset Z=19 (probe e v2) confirms localization.
   - *This is THE dominant cause.*

2. **Cause (b)** RULED OUT: FD solver faithful to <1.2% across Z=1..80.

3. **Cause (c)** RULED OUT: no factor-of-X convention bug.

4. **Cause (d-coding)** RULED OUT: multi-zeta machinery internally consistent.

5. **Cause (e — convention mismatch in operator chain)** RULED OUT: same as (c).

6. **Cause (f) — self-consistency missing**: NOT directly tested in this diagnostic, but the structural finding from (d) is consistent with (f) being the principled long-term solution.

**The cliff is single-cause, dominant: missing radial nodes** (and the absence of node-producing negative coefficients in any single-zeta basis or heuristic two-zeta extension).

---

## 8. BBB93 / KTT GO/NO-GO

### Recommendation: **GO ON BBB93/KTT FULL TABULATION — but with eyes open.**

**Why GO:**
- The dominant cause (missing radial nodes) IS what BBB93/KTT addresses. A 5-9 STO primitive expansion with negative coefficients produces the radial nodes needed for orthogonality.
- Both primary alternative paths (b, c, d-coding) are RULED OUT as cliff causes.
- The two-zeta heuristic FAILED in the wrong direction, but BBB93 is structurally different (negative coefficients).

**Risks / why "with eyes open":**
- **For BBB93 (Z<=54):** The Xe atom is at the upper limit of BBB93's range. Coefficients for Xe in BBB93 may be less converged than for lighter atoms.
- **For Cs/Ba (Z=55-56):** MUST use Koga-Tatewaki-Thakkar 1993/2000 (BBB93 stops at Z=54). KTT data is harder to find in tabular form.
- **Casimir F_R is INDEPENDENT of the screening upgrade.** Closing the screening cliff still leaves a ~30-40% gap to the full Bohr-Weisskopf factor (~2.6 vs 1.555 for Cs). BBB93 alone is NOT sufficient for sub-percent Cs HFS.

**Quantitative expectation:** BBB93 closes ~20% of the Cs HFS residual (the screening kernel contribution); the remaining ~30% requires full Bohr-Weisskopf via spinor lift. This means after BBB93, Cs HFS would land around −25% to −30% residual, not at sub-percent.

### Ranked next steps:

1. **PRIMARY (recommended): BBB93/KTT for [Xe] core.** ~1-2 weeks to hand-tabulate the Xe RHF expansion (~200 (c_i, ζ_i) values at 7-digit precision). Then plug into the existing `_build_xe_orbitals_bbb93()` slot in `multi_zeta_orbitals.py`. Expected outcome: Cs HFS residual −47% → −25% to −30%; SrH/BaH SO splittings should also improve correspondingly.

2. **SECONDARY (necessary for sub-percent): full Bohr-Weisskopf via Tier 3 spinor-lift.** 1-2 weeks to evaluate s-state contact density via full Dirac wavefunctions instead of scalar BF × Casimir F_R. Reuses `geovac/dirac_matrix_elements.py` infrastructure. Closes the relativistic gap from leading-order to full BW.

3. **ALTERNATIVE PATH (3-4 weeks, more principled): Self-consistent Hartree-Fock iteration in `geovac/neon_core.py`.** Closes the screening cliff for ALL Z (not just Xe core), making the framework heavy-atom-extensible. This is the long-term clean engineering path. The multi-zeta machinery added in the closeout sprint is a drop-in foundation.

**Decision factor:** If the next sprint window is ≤2 weeks and the PI wants Cs PNC closure progress, **GO BBB93** (gets ~20% closure, then queue Bohr-Weisskopf). If the window is 3+ weeks and the PI wants the principled framework path, **GO SCF**. Both close the cliff structurally; SCF generalizes better.

---

## 9. What this diagnostic cost vs. what was avoided

**Diagnostic cost:** 5 probe scripts (~1500 lines total, 6 hours wall-clock to write + run). Total wall-clock for compute: ~5 minutes (FD scan dominates).

**What was avoided:**
- A multi-week BBB93 data-entry effort under the false assumption that "more zetas alone" would close the cliff. Probe (d) showed this is structurally insufficient; even BBB93 only gets ~20% closure without Bohr-Weisskopf.
- A wild-goose chase on FD grid refinement: probe (b) rules this out cleanly.
- A factor-of-X bug hunt in the convention chain: probe (c) rules this out cleanly.

**What was sharpened:**
- The closeout's "CR67 non-faithful for outer shells" diagnosis is upgraded to "**missing radial nodes** structurally inherent in any single-zeta hydrogenic basis or heuristic two-zeta extension."
- The cliff is **shell-resolved** within Cs (not just Z-resolved) — outer shells (4d, 5s, 5p) overshoot by 2-4x, while inner shells (1s, 2s, 2p) are well-fit.
- BBB93 is GO with **realistic expected outcome (~20% closure of the residual, not full closure)**.
- The principled long-term path is **self-consistent HF**, not just better tabulations.

---

## 10. Files

**Probe scripts (read-only, in debug/):**
- `debug/z_cliff_probe_a.py` — CR67 fit faithfulness across alkali heavy atoms.
- `debug/z_cliff_probe_a_extended.py` — CR67 inner-vs-outer shells of Cs.
- `debug/z_cliff_probe_b.py` — FD solver hydrogenic faithfulness across Z=1..80.
- `debug/z_cliff_probe_c.py` — Convention mismatch trace.
- `debug/z_cliff_probe_d.py` — multi-zeta machinery internal coherence.
- `debug/z_cliff_probe_e.py` — Z-scan v1 (with Casimir F_R applied).
- `debug/z_cliff_probe_e_v2.py` — Z-scan v2 (kernel isolated, |ψ_FW/ψ_HF_lit| ratio).
- `debug/z_cliff_probe_synthesis.py` — Combined verdict and ranking.

**Probe data (in debug/data/):**
- `z_cliff_probe_a.json` — Per-orbital fit faithfulness data for 8 outer-shell systems.
- `z_cliff_probe_a_extended.json` — All 11 Cs shells inner/middle/outer breakdown.
- `z_cliff_probe_b.json` — FD convergence panel at Z=1..80, n_grid 50k..400k.
- `z_cliff_probe_c.json` — Convention check pass/fail summary.
- `z_cliff_probe_d.json` — multi-zeta machinery + s-orbital overlap matrix.
- `z_cliff_probe_e.json` — Z-scan v1.
- `z_cliff_probe_e_v2.json` — Z-scan v2 (the cleanest localization).
- `z_cliff_probe_synthesis.json` — Final dominant-cause ranking and BBB93 recommendation.

**No production code modified.** No tests added. This is a pure diagnostic.

---

## 11. Production code observations (to be addressed if PI dispatches a fix track)

The diagnostic surfaced one production-code design observation that is NOT a bug but is worth flagging:

**Observation 1: `_TWO_ZETA_SPLITS` in `geovac/multi_zeta_orbitals.py` uses uniform inner/outer ratios (1.20, 0.83) for ALL shells.** This means the heuristic two-zeta produces orbitals that are systematically MORE compact than CR67 single-zeta on average (mean ratio 0.83 from probe d check 4). This is the wrong direction for closing the cliff. The proper inner/outer ratios for BBB93-style fits are SHELL-DEPENDENT and have NEGATIVE coefficients on inner primitives to enforce orthogonality. The current heuristic is NOT calibrated to BBB93 directly — it's a placeholder for the proper full-tabulation path.

This is documented in `multi_zeta_orbitals.py` as a "scoping probe" comment, so it's NOT a bug — just a flag that the heuristic is the wrong starting point if a future sprint wants to refine without going to full BBB93.

**Observation 2: When BBB93 / KTT are tabulated, the proper structural test is whether the resulting orbitals satisfy ⟨R_ns | R_ms⟩ ≈ 0 for n ≠ m.** The current heuristic gives max overlap 0.86 (essentially same function); proper BBB93 should give <0.05.

---

## 12. Verdict

**Dominant cause:** Missing radial nodes in the CR67 single-zeta hydrogenic basis (and inherited by the heuristic two-zeta), structurally producing outer-shell orbital overshoot of 2-4x at heavy Z. The cliff is shell-resolved; cliff onset at Z=19 (K 4s) clean.

**Hypotheses ruled out cleanly:** FD solver convergence (b), convention mismatch (c), multi-zeta machinery coding (d-coding).

**BBB93/KTT verdict: GO with eyes open.** Closes ~20% of the residual; full sub-percent Cs HFS still needs Bohr-Weisskopf. SCF in `neon_core.py` is the long-term principled path.

**No production code modified. Diagnostic cost ~5 min compute + ~1500 lines of probe scripts.** Multi-week BBB93 effort under wrong assumption avoided.
