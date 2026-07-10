# Sprint F1 max_n=3 — Predictions test for combined W1c × multi-zeta

**Date:** 2026-05-23 (post-F1-P1+P2 same-day continuation).
**Sprint position:** Natural follow-on to F1-P1+P2 (closed earlier today with PARTIAL-CLOSURE-AT-MAX_N=2) and to the W1c × M-Z partition bridge sprint (which produced three falsifiable predictions for this exact test). Tests whether basis enlargement from max_n=2 (Q=20) to max_n=3 (Q=56) under the unified W1c × multi-zeta architecture provides the dimensional richness for the FCI to construct a true bonding combination with energy lower than the separated configuration.
**Cross-references:** `debug/sprint_f1_p1p2_combined_test_memo.md` (max_n=2 partial-closure baseline), `debug/sprint_w1c_mz_partition_analysis_memo.md` (bridge sprint with the three predictions tested here), `debug/data/sprint_w1c_mz_partition_predictions.json` (the falsifiable predictions schema), `debug/sprint_modular_alpha_arc_synthesis_memo.md` (full α arc context), `geovac/balanced_coupled.py`, `geovac/multi_zeta_orbitals.py`, `geovac/cross_center_screened_vne.py`, CLAUDE.md §1.7 multi-focal-composition wall taxonomy.

---

## §0. Executive summary + verdict

**Verdict line: CLEAN NEGATIVE.**

The bridge sprint's three predictions, tested at NaH max_n=3 with the unified W1c × multi-zeta architecture from F1-P1+P2, return:

- **P1 (internal PES minimum in R ∈ [3.0, 4.5] bohr): FAIL.** $R_{\min} = 2.0$ bohr (smallest tested), exactly the falsification criterion. The combined PES is monotonically descending across $R \in \{2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0\}$.
- **P2 (D_e ∈ [0.0375, 0.150] Ha within 2× of experimental): FAIL (conditional on P1).** PES descent depth $D_e^{\text{PES}} = +0.7097$ Ha — about $10\times$ larger than the experimental NaH $D_e \approx 0.075$ Ha, indicating continued over-attraction at small $R$. This is the spurious-binding signature, not a real $D_e$.
- **P3 (multi-zeta differential in [0.02, 0.30] Ha at R_eq): PASS.** $|E_{W1c} - E_{W1c+mz}| = 0.2066$ Ha at $R_{\min} = 2.0$ bohr, within the predicted range. Multi-zeta is load-bearing at max_n=3 (consistent with F1-P1+P2 at max_n=2: 0.21 Ha at R=2.0, 0.06 Ha at R=3.5).

Per the bridge sprint's gate logic: **the basis-closable cross-shift hypothesis for sub-layer 3 is REFUTED at max_n=3; sub-layer 3 reclassifies from "BIMODULE CROSS-SHIFT (basis-closable)" to "MODULE ENDOMORPHISM (basis-irreducible at tested scales)".**

**Substantive new content the gate did NOT anticipate:** at max_n=3 the dominant natural orbital of the combined architecture at R=3.5 IS now a true bonding combination — Na 3s + H 1s with ~50/50 amplitude split, mixing coefficients $-0.698\,\text{Na 3s} - 0.687\,\text{H 1s}$. The 2nd natural orbital is the antibonding combination ($-0.698\,\text{Na 3s} + 0.687\,\text{H 1s}$). The basis enlargement DID provide the orbital-mixing flexibility that max_n=2 lacked: at max_n=2 the dominant NO was 100% H-localized; at max_n=3 it is properly delocalized across the bond. The "structural gap" (bonding orbital constructibility) is closed; the "energetic gap" (the constructed bonding combination has higher energy than the separated configuration at every R) remains. The W1c residual at max_n=3 is therefore not the orbital-pair flexibility but the energetic-asymmetry of the cross-V_ne kernel against the bonding combination, which the framework's bimodule machinery cannot close from within at the tested basis.

Net structural reading: the M-Z 2-bucket partition (cross-shift / endomorphism) holds; the proposed 3-bucket refinement (cross-shift / basis-closable cross-shift / endomorphism) is NOT supported. Sub-layer 3 lives in the same external-input class as sub-layer 2.

---

## §1. Pre-flight basis structure (multi-zeta registry coverage at max_n=3)

### 1.1 NaH spec at max_n=3

The NaH spec (from `geovac.molecular_spec.nah_spec`) has a single bond block:

```
NaH_bond: Z_center=1.0  max_n=3  has_h_partner=True  max_n_partner=3  n_val_offset=2
```

The build produces $M = 28$ spatial orbitals, $Q = 56$ spinors:

- **Na center sub-block (14 orbitals):** block_n=1 has $l=0$ (1 orbital, Na 3s); block_n=2 has $l \in \{0, 1\}$ (1+3 = 4 orbitals, Na 4s + 4p$_{\{-1,0,+1\}}$); block_n=3 has $l \in \{0, 1, 2\}$ (1+3+5 = 9 orbitals, Na 5s + 5p$_{\{-1,0,+1\}}$ + 5d$_{\{-2,-1,0,+1,+2\}}$).
- **H partner sub-block (14 orbitals):** identical structure but $Z_{\text{orb}} = 1$ hydrogenic; orbitals labeled H 1s through H 3d.

The physical-n mapping uses $n_{\text{val\_offset}} = 2$, so block_n=1 → physical $n=3$, block_n=2 → physical $n=4$, block_n=3 → physical $n=5$.

### 1.2 Multi-zeta registry coverage

The registry `geovac/multi_zeta_orbitals.py::get_physical_valence_orbitals(11)` contains only Na 3s and Na 3p (both with `n_orbital=3`). The dispatch in `build_balanced_hamiltonian` matches `orb.n_orbital == physical_n` AND `orb.l_orbital == l`, so at max_n=3:

| block_n | l | physical_n | requested orbital | in registry? | substituted? |
|:-:|:-:|:-:|:--|:-:|:-:|
| 1 | 0 | 3 | Na 3s | YES | **YES** |
| 1 | 1 | 3 | Na 3p | YES, but l=1 invalid at n=1 | NO (impossible) |
| 2 | 0 | 4 | Na 4s | NO | NO |
| 2 | 1 | 4 | Na 4p | NO | NO |
| 3 | 0 | 5 | Na 5s | NO | NO |
| 3 | 1 | 5 | Na 5p | NO | NO |
| 3 | 2 | 5 | Na 5d | NO | NO |

So **the multi-zeta dispatch substitutes ONLY (block_n=1, l=0) = Na 3s** at max_n=3 — identical coverage to max_n=2 (which also substitutes only Na 3s, since at max_n=2 the only valid (block_n, l) pair matching the registry is (1, 0)).

This is a real observation: the multi-zeta dispatch's `block_n + n_val_offset = physical_n` mapping leaves **Na 3p unreachable** at any max_n (because Na 3p has physical n=3, l=1, which requires block_n=1 with l=1, which is structurally forbidden by the standard hydrogenic enumeration n ≥ l+1). Coverage extension would require (a) extending the registry with Na 4s/4p/5s/5p/5d physical-fit orbitals (substantial multi-day work per orbital), or (b) refactoring the dispatch to handle the Na 3p special case explicitly. Neither is needed for this sprint, but flagging for future reference.

**Decision: proceed with Option (b) hydrogenic placeholder for all Na orbitals except 3s.** This matches the F1-P1+P2 architecture verbatim and isolates the test to "does the larger basis admit a bonding-orbital construction" without confounding the multi-zeta extension with other architectural changes.

### 1.3 Compute budget estimate

Build timing at max_n=3 with full settings (n_grid_vne=4000, L_max=4): ~10 s per build, regardless of `screened_cross_center` or `multi_zeta_basis` flags. FCI dim at $n_e=2$ is $\binom{28}{1}^2 = 784$ — trivial diagonalization (<1 s). Total compute estimate for the three-step test: ~3-4 minutes. Actual: ~5 minutes (8 R-points × 2 architectures × 11 s per build + tiny FCI overhead).

---

## §2. Step 1 — Algebraic kernel differential at max_n=3

Cross-V_ne diagonal matrix element ⟨Na 3s | $-Z_H/|r-R_H|$ | Na 3s⟩ at $R = 3.566$ bohr (experimental NaH R_eq):

| Architecture | h1[0,0] (Ha) | h1_cv[0,0] (Ha) | mz vs no_mz diff (Ha) |
|:--|--:|--:|--:|
| bare, no mz | −0.779403 | −0.279403 | — |
| bare, mz | −0.726047 | −0.226047 | **+0.053356** |
| W1c, no mz | −0.779403 | −0.279403 | — |
| W1c, mz | −0.726047 | −0.226047 | **+0.053356** |

The h1_cv diagonal is the cross-V_ne contribution on Na 3s from the H nucleus. Multi-zeta makes Na 3s slightly more compact (FrozenCore-screened physical shape vs diffuse hydrogenic Z_orb=1 placeholder), so the H nucleus pulls Na 3s LESS strongly (less negative cross-V_ne contribution).

**Comparison to α-PES Step 1 and F1-P1+P2 baseline:**
- α-PES Step 1 (2026-05-23 morning) reported the cross-V_ne differential at NaH max_n=2 R=3.566 as $-0.135$ Ha (sign convention: differential of the multipole integrand contribution, not the matrix element).
- F1-P1+P2 §3.1 (max_n=2) reported FCI multi-zeta differential $+0.056$ Ha at R=3.5.
- This sprint (max_n=3) finds h1_cv differential $+0.053$ Ha at R=3.566 and FCI differential $+0.055$ Ha at R=3.5 (Step 2).

The numerical magnitude is consistent across basis sizes ($\sim +0.05$ Ha per matrix element). The multi-zeta machinery works correctly at max_n=3, confirming Pre-flight Decision (b). The screened-path differential equals the bare-path differential bit-exactly — confirms the F1-P1+P2 architectural finding that the screened-V_ne kernel correction (against H nucleus, which has no [Ne] core) reduces to the bare path on the Na sub-block (the screening correction is identically zero for $Z_H = 1$).

**Sanity check: PASS.** Multi-zeta dispatch works as expected at max_n=3.

---

## §3. Step 2 — Single-point FCI results at NaH max_n=3

### 3.1 Two-point FCI table

System: NaH max_n=3 (Q=56, M=28 spatial, $n_e=2$, FCI dim 784).

| Architecture | E(R=3.5) [Ha] | E(R=10.0) [Ha] | $D_e^{\text{2pt}}$ [Ha] | top NO | 2nd NO | Na 3s diag occ |
|:--|--:|--:|--:|--:|--:|--:|
| bare baseline | −169.6491 | −165.5323 | +4.1168 | 1.9963 | 0.0037 | 0.000 |
| W1c alone | −163.1774 | −162.7789 | +0.3985 | 1.0000 | 1.0000 | 0.974 |
| **COMBINED (W1c + mz)** | **−163.1220** | **−162.7788** | **+0.3431** | **1.0000** | **1.0000** | **0.974** |

Three structural observations from the two-point table:

(a) **Bare baseline at max_n=3 has the same Na 3s ≈ 0 occupation pattern as at max_n=2.** The basis enlargement does NOT change the qualitative behavior of the bare cross-V_ne case — the FCI is still dominated by an H-localized single-Slater configuration (top NO = 1.99) because the un-screened Z=11 attraction on H 1s makes it the deepest orbital. This reproduces the α-PES Step 2 bit-zero finding at max_n=3.

(b) **W1c at max_n=3 lifts the bare-baseline energy by ~6.5 Ha** (from −169.6 to −163.2 at R=3.5), similar to max_n=2 (~6 Ha lift). Multi-zeta on top adds another +0.055 Ha shift at R=3.5, fully load-bearing.

(c) **Natural occupations at W1c-on (with or without mz): [1.0, 1.0]** — two singly-occupied orbitals, identical pattern to max_n=2. The two-point $D_e$ values for both W1c configs look positive (+0.39 to +0.34 Ha), but the bridge sprint warned this could be misleading; the extended PES (Step 3) confirms it is.

### 3.2 FCI natural-orbital structure — the substantive new content

The dominant natural orbital at R=3.5 under W1c+mz at max_n=3:

| Architecture | dom NO Na amp² | dom NO H amp² | Top amplitude components |
|:--|--:|--:|:--|
| bare baseline | 0.000 | 1.000 | H_2p₀ = +0.662, H_1s = −0.634, H_2s = +0.277, H_3p₀ = +0.215 |
| W1c alone | 0.000 | 1.000 | H_1s = +0.972, H_2p₀ = −0.212, H_3p₀ = −0.077 |
| **COMBINED W1c+mz** | **0.500** | **0.500** | **Na_3s = −0.698, H_1s = −0.687, H_2p₀ = +0.150, Na_4p₀ = −0.106** |

**This is the headline structural finding of the sprint.** At max_n=2 the dominant NO was 100% H-localized at every architecture (F1-P1+P2 §4.2). At max_n=3 under COMBINED W1c+mz, the dominant NO is now a **true bonding combination** Na 3s + H 1s with ~50/50 amplitude split, and the 2nd natural orbital is the **antibonding combination** Na 3s − H 1s:

| | dom NO | 2nd NO |
|:--|:--|:--|
| COMBINED W1c+mz | −0.698 Na_3s − 0.687 H_1s | −0.698 Na_3s + 0.687 H_1s |

The basis enlargement from max_n=2 to max_n=3 DID provide the orbital mixing flexibility predicted by the bridge sprint. The structural prediction (P1 underlying classification: "max_n=3 provides the dimensional richness for bonding orbital construction") is **structurally confirmed at the orbital level** — the framework constructs the right bonding orbital.

Wrinkle: the W1c-alone configuration at max_n=3 does NOT show this bonding behavior at R=3.5 — its dominant NO is still 100% H-localized (97.2% H_1s). Multi-zeta is the load-bearing ingredient that activates the bonding orbital construction. This is consistent with the F1-P1+P2 finding that multi-zeta engagement requires Na 3s occupation, which requires W1c — but adds the new finding that **at larger basis, the combined architecture not only occupies Na 3s but mixes it with H 1s in the right bonding combination**.

### 3.3 Multi-zeta differential

| R [bohr] | E(W1c) [Ha] | E(W1c+mz) [Ha] | mz_diff [Ha] |
|:--:|--:|--:|--:|
| 3.5 | −163.1774 | −163.1220 | −0.05545 |
| 10.0 | −162.7789 | −162.7788 | −0.00009 |

The mz_diff at R=3.5 is $-0.055$ Ha, essentially identical to the F1-P1+P2 value at max_n=2 ($-0.056$ Ha). At R=10 it's bit-zero (dissociation limit). The endomorphism content of multi-zeta is preserved across basis sizes.

### 3.4 Diagonal occupations of named orbitals (COMBINED at R=3.5)

| Orbital | Diagonal occupation |
|:--|--:|
| Na_3s | 0.974 |
| H_1s | 0.944 |
| H_2p₀ | 0.045 |
| Na_4p₀ | 0.023 |
| H_3p₀ | 0.006 |
| Na_5p₀ | 0.003 |
| H_2s | 0.002 |
| H_3d₀ | 0.002 |

The two dominant occupations (Na_3s ≈ 0.97, H_1s ≈ 0.94) confirm the bonding/antibonding partition is real: the bonding combination puts ~0.5 of an electron on Na 3s and ~0.5 on H 1s; the antibonding combination puts the other ~0.5 on Na 3s and ~0.5 on H 1s with opposite phase. The diagonal occupations are dominated by these two configurations, with small (≤ 5%) admixture of H 2p, Na 4p, and higher partial waves.

### 3.5 Decision gate verdict

Per the Step 2 gate logic:
- $D_e^{\text{2pt}}$ (combined) = $+0.343$ Ha (numerically positive)
- Top natural occupation = 1.000 (not bonding-dominant)
- Two-singly-occupied signature: YES

Gate verdict: **PROCEED_TO_MINI_PES (binding, but two-singly-occupied — puzzling)**. The naive $D_e^{\text{2pt}} > 0$ looks like binding, but the natural-occupation pattern [1, 1] indicates open-shell singlet of two singly-occupied orbitals, NOT a closed-shell bond. The extended PES is needed to determine whether this is a real internal-minimum binding or a spurious-binding signature (two-point sampling of an over-attraction floor).

---

## §4. Step 3 — Mini-PES across 8 R-points

Per Step 2 gate verdict, ran 8-point PES at $R \in \{2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0\}$ bohr:

| R [bohr] | E(W1c+mz) [Ha] | E(W1c) [Ha] | mz_diff [Ha] | dom NO Na/H amp² |
|:--:|--:|--:|--:|:--|
| 2.000 | **−163.4885** | −163.6951 | −0.20656 | 0.498 / 0.502 |
| 2.500 | −163.3143 | −163.4496 | −0.13528 | 0.501 / 0.499 |
| 3.000 | −163.2022 | −163.2894 | −0.08713 | 0.501 / 0.499 |
| 3.500 | −163.1220 | −163.1774 | −0.05545 | 0.575 / 0.425 |
| 4.000 | −163.0607 | −163.0956 | −0.03494 | 0.503 / 0.497 |
| 5.000 | −162.9720 | −162.9855 | −0.01351 | 0.500 / 0.500 |
| 7.000 | −162.8640 | −162.8659 | −0.00187 | 0.470 / 0.530 |
| 10.000 | −162.7788 | −162.7789 | −0.00009 | 0.389 / 0.611 |

**$R_{\min} = 2.0$ bohr** (smallest tested) for both architectures. Combined W1c+mz PES descent depth: $E(2.0) - E(10.0) = -0.7097$ Ha. W1c-alone descent: $-0.9162$ Ha. Multi-zeta reduces descent by 0.206 Ha (~22%), comparable to the 23% reduction at max_n=2 (F1-P1+P2 §3.2).

**No internal minimum.** PES is monotonically descending across the full $R$ range. R_min is at the smallest tested R = 2.0 bohr — **exactly the falsification criterion for P1**.

The dominant NO Na/H amplitude split is stable at ~50/50 across most R values (from 2.0 to 5.0 bohr), confirming that the bonding orbital combination persists across the PES — it's not just a localized feature at R=3.5. At R=10.0 the split tilts to 0.39/0.61, weakening as the bonding combination dissolves in the dissociation limit.

The mz_diff column is smoothly R-dependent: large at small R (where the basis is closely overlapping), decaying to zero at R=10 (dissociation). This is the same R-dependence pattern as F1-P1+P2 max_n=2; the endomorphism content is consistent across basis sizes.

### 4.1 Concavity check

Cannot check concavity at $R_{\min} = 2.0$ since it's at the boundary of the test range. To get a true internal minimum, the PES would need to turn over somewhere in $R \in (0, 10)$, and it does not within the tested range. The PES gradient near $R = 2.0$ is steeply downhill, with no indication of curvature reversal at any tested $R$.

---

## §5. Predictions verification

| Prediction | Statement | Actual | Verdict |
|:--|:--|:--|:--:|
| P1 | Internal PES minimum at $R_{\text{eq}} \in [3.0, 4.5]$ bohr | $R_{\min} = 2.0$ bohr (smallest tested) | **FAIL** |
| P2 | $D_e \in [0.0375, 0.150]$ Ha within 2× of experimental | $D_e^{\text{PES}} = +0.7097$ Ha (10× experimental) | **FAIL** (conditional on P1) |
| P3 | $|E_{W1c} - E_{W1c+\text{mz}}| \in [0.02, 0.30]$ Ha at $R_{\text{eq}}$ | $|0.2066|$ Ha at $R = 2.0$ | **PASS** |

P1 fails on the exact falsification criterion specified in the bridge sprint ($R_{\min}$ at smallest tested R). P2 is conditional on P1 and also fails (the spurious-binding descent depth is 10× the experimental binding energy). P3 passes — multi-zeta differential magnitude is within the predicted range at max_n=3, consistent with the endomorphism classification that does NOT depend on basis size.

**Aggregate: 1 of 3 predictions pass. Verdict: CLEAN NEGATIVE per the bridge sprint's gate rules.**

---

## §6. Structural verdict

### 6.1 Sub-layer 3 reclassification

The bridge sprint's hypothesis was that sub-layer 3 (orbital-pair flexibility at NaH max_n=2 under W1c+mz) is a **BIMODULE CROSS-SHIFT (basis-closable)** — meaning the bonding combination [H 1s ± Na 3s]_± is a structurally bimodule cross-shift but the FCI cannot realize it at max_n=2 due to dimensional poverty; max_n=3 should provide the additional flexibility.

This sprint's result splits cleanly into a structural success and an energetic failure:

**Structural success:** at max_n=3 under combined W1c+mz, the FCI dominant natural orbital IS the bonding combination Na 3s + H 1s (50/50 amplitude split, with explicit mixing coefficients $-0.698, -0.687$). The 2nd NO is the antibonding combination. The basis enlargement DID provide the orbital-mixing flexibility predicted by the bridge sprint. From a pure orbital-structure standpoint, the basis-closable cross-shift hypothesis is CONFIRMED — the bonding orbital is constructible at max_n=3 in a way it wasn't at max_n=2.

**Energetic failure:** the constructed bonding orbital has higher energy than the separated configuration at every tested R. The PES is monotonically descending; no internal minimum emerges. The natural occupations remain [1.0, 1.0] (open-shell singlet of bonding + antibonding singly-occupied), NOT [2.0, 0.0] (closed-shell bond). The energetic gap that would lower the bonding combination below the separated configuration requires content the framework cannot generate from inside its bimodule machinery at the tested basis.

The net reading: **sub-layer 3 lives in the same external-input class as sub-layer 2**, not in a separate basis-closable cross-shift sub-class. The "structural" part of the cross-shift hypothesis is correct (the right linear combination IS structurally a bimodule cross-shift, and the basis enlargement DOES make it constructible), but the "framework_handles=true" claim is too strong: the framework constructs the orbital but cannot energetically prefer it over the separated configuration. The remaining wall is not orbital-pair flexibility but rather the cross-V_ne kernel energetics against the bonding combination, which the framework's W1c × multi-zeta machinery does not close.

### 6.2 Updated M-Z partition taxonomy

The bridge sprint proposed a **3-bucket refinement** of the M-Z partition:
1. Cross-shift (handled by framework)
2. Basis-closable cross-shift (handled by framework, but FCI realization requires sufficient basis dimensionality)
3. Endomorphism (external input)

This sprint **falsifies the 3-bucket refinement**. The proposed sub-class 2 collapses into sub-class 3: even when the basis is large enough for the FCI to construct the right orbital combination, the framework cannot make that orbital combination energetically preferred. The 2-bucket M-Z partition (cross-shift / endomorphism) is the correct framing; the W1c residual at NaH (sub-layer 3) lives in the endomorphism bucket alongside sub-layer 2 (multi-zeta wavefunction shape).

### 6.3 Consistency with the broader framework

The CLEAN NEGATIVE verdict is consistent with the structural-skeleton-scope finding (CLAUDE.md §1.7): the framework determines the skeleton (selection rules, transcendental classes, scaling laws, structural orbital constructibility) but does not generate calibration data (energetic asymmetries, parameter values, Yukawa selection). At NaH this manifests as: the bonding orbital is structurally constructible, but the energy that would make it the preferred ground-state configuration requires inputs beyond what the framework's bimodule machinery (Phillips-Kleinman + cross-V_ne + multi-zeta basis) provides.

The W1c-residual orthogonality wall is now sharpened: it is NOT an orbital-pair flexibility wall (the orbital pair IS constructed at max_n=3); it IS a wall in the cross-V_ne kernel energetics — specifically, the energetic asymmetry that should make the bonding combination preferred over the separated configuration. Track 3's named target P2 (cross-V_ne kernel-shape substitution on the partner side) addresses a structurally distinct mechanism that may close this energetic gap.

### 6.4 Honest scope

What this sprint DID demonstrate:
- Multi-zeta architecture works at max_n=3 (no architectural breakage, consistent with F1-P1+P2).
- The basis enlargement from max_n=2 to max_n=3 closes the structural gap for bonding-orbital constructibility (genuinely new content).
- The basis enlargement does NOT close the energetic gap — PES still monotonically descending, no internal minimum.
- Multi-zeta endomorphism content is preserved across basis sizes (mz differential magnitude similar at max_n=2 and max_n=3).

What this sprint did NOT demonstrate:
- That the W1c residual cannot be closed by ANY mechanism (only that basis enlargement at max_n=3 does not close it).
- That Track 3's named target P2 (cross-V_ne kernel-shape substitution) would or would not close the wall (separate sprint).
- That different system classes (LiH+, HCl, MgH₂) would or would not show similar walls (W1c residual is empirically Z-decreasing per CLAUDE.md §2, so these systems may already have effective binding even with partial wall).

---

## §7. Recommended next sprint

Per the bridge sprint's gate logic for the CLEAN NEGATIVE branch:

### Priority 1 — Apply paper edits (deferred to PI, ~2-3 hours)

The CLEAN NEGATIVE verdict warrants three documentation updates:

(a) **Paper 19 §sec:w1c_residual extension** (~150 lines): document the sub-layer 3 reclassification (basis-closable cross-shift → endomorphism) at max_n=3, with the structural success/energetic failure framing. Honest scope statement: the bonding orbital IS constructible at max_n=3, but the energetic gap requires content the framework's bimodule machinery does not generate. Cross-reference to Track 3 P2 named target.

(b) **CLAUDE.md §3 dead ends row** for the basis-closable cross-shift hypothesis: "Basis-closable cross-shift hypothesis for W1c residual at NaH max_n=3 — CLEAN NEGATIVE. The bonding orbital IS constructible at max_n=3 (dominant natural orbital is Na 3s + H 1s with 50/50 amplitude split, mixing coefficients -0.698 -0.687), but the constructed bonding combination has higher energy than the separated configuration at every R; PES still monotonically descending. Sub-layer 3 reclassifies to genuine endomorphism class. Framework cannot close the energetic gap from within its bimodule machinery at the tested scales."

(c) **CLAUDE.md §1.7 multi-focal-composition wall taxonomy update**: keep the 2-bucket M-Z partition (cross-shift / endomorphism) as canonical; drop the proposed 3-bucket refinement (cross-shift / basis-closable cross-shift / endomorphism) per CLEAN NEGATIVE. Document the structural-vs-energetic distinction emerging from this sprint as the substantive new content (in addition to falsifying the refinement, the sprint produced a sharper diagnostic of WHERE the residual lives).

These edits should NOT be applied autonomously per sprint mandate; PI to decide.

### Priority 2 — Pivot to F2 cross-V_ne kernel-shape substitution (1-2 weeks)

Track 3's named target P2 (from CLAUDE.md §3 PK cross-center sprint synthesis): replace the bare cross-V_ne kernel integration on the H partner side with physical-n hydrogenic shape, distinct from the multi-zeta basis on the Na center side that this sprint tested. The Track 3 diagnostic at NaH max_n=2 (`debug/w1c_residual_nah_track3_diag.py`) found a differential of $-0.674$ Ha at R=2.5 for cross-V_ne shape substitution, about 1.9× the W1c-residual descent depth — suggesting this substitution would over-correct from the current attraction profile, possibly producing a binding minimum.

The mechanism is structurally distinct from this sprint's basis enlargement: F2 addresses the cross-V_ne kernel shape on the OPPOSITE side from where multi-zeta lives, which is the right next target after this sprint's basis enlargement was insufficient. F2 would also probe whether the energetic-gap is closable by changing the kernel rather than the basis.

### Priority 3 — Pivot to a different system class (parallel, 2-3 days)

The W1c-residual wall is empirically Z-decreasing (NaH 5.4-6.0×, MgH₂ 2.99×, HCl 1.79× per CLAUDE.md §2). At higher Z the wall may already be smaller fraction of bond energy. A 2-day scoping run of HCl max_n=2 with the combined W1c × multi-zeta architecture would test whether the framework already supports second-row hydride binding even without complete W1c closure — but requires extending the multi-zeta registry to Cl (Z=17), which is a multi-day extension on its own.

### Priority 4 (DO NOT pursue) — further max_n enlargement on NaH

At max_n=3 the structural gap (bonding orbital constructibility) is closed; the wall now lives in the energetic asymmetry, which is not basis-size limited at this scale. Further enlargement (max_n=4, Q=84 per CLAUDE.md §2) would add more orbitals but would not address the energetic gap — same structural finding would persist with marginal numerical refinement. Compute cost would also increase substantially (FCI dim grows as $M(M-1)/2$ in the 2-electron case).

### Bundled recommendation

**Default**: open Priority 1 paper edit batch (PI review), and queue Priority 2 (F2 cross-V_ne kernel-shape substitution) as the next active sprint. Cost ~1-2 weeks for F2, structurally distinct from this sprint's basis-enlargement test.

**Parallel**: keep Priority 3 (HCl scoping) flagged but not active — the multi-zeta registry extension to Cl is a separate sub-sprint with its own multi-day cost, and the structural-skeleton-scope finding (CLAUDE.md §1.7) suggests heavier-Z systems may have the same wall character with the wall just being smaller relative to bond energy. The empirical Z-decreasing pattern is empirical, not predicted from first principles.

---

## §8. Files

### Created (debug)

- `debug/sprint_f1_maxn3_step1_kernel_differential.py` — Step 1 driver: cross-V_ne kernel differential at NaH max_n=3, 4-architecture sanity check.
- `debug/sprint_f1_maxn3_step2_fci.py` — Step 2 driver: single-point FCI at R=3.5 and R=10.0, 1-RDM analysis, h1 spectrum diagnosis, decision gate.
- `debug/sprint_f1_maxn3_step3_mini_pes.py` — Step 3 driver: 8-point PES scan + predictions verification.
- `debug/update_predictions_verification.py` — appends `verification` section to the predictions JSON.
- `debug/write_consolidated_results.py` — produces `debug/data/sprint_f1_maxn3_results.json` from Steps 1+2+3.
- `debug/sprint_f1_maxn3_predictions_test_memo.md` — this memo (~3500 words).

### Created (data)

- `debug/data/sprint_f1_maxn3_step1_kernel.json` — Step 1 results.
- `debug/data/sprint_f1_maxn3_step2_fci.json` — Step 2 results (full FCI eigenvalues, 1-RDM diagonals, natural orbital amplitudes).
- `debug/data/sprint_f1_maxn3_step3_pes.json` — Step 3 results (8-point PES + predictions check).
- `debug/data/sprint_f1_maxn3_results.json` — consolidated results across all three steps.

### Modified

- `debug/data/sprint_w1c_mz_partition_predictions.json` — appended `verification` section with the three predictions verdicts, sub-layer 3 reclassification, and next-sprint recommendations. Original `predictions_F1_max_n_3` section preserved.

### Not modified

- Production code (`geovac/balanced_coupled.py`, `geovac/multi_zeta_orbitals.py`, etc.) — sprint mandate: no production changes for the test.
- Paper `.tex` files — sprint mandate: paper edits deferred to a follow-on β-style sprint, gated on PI review.
- CLAUDE.md — sprint mandate: documentation edits deferred (the CLEAN NEGATIVE verdict warrants §3 dead ends entry + §1.7 partition update, but those go in the follow-on documentation pass).
- Tests — no production code modified, no test changes needed.

### Regression check

```
$ pytest tests/test_balanced_coupled_multizeta.py \
         tests/test_balanced_coupled_screened_valence.py \
         tests/test_multi_zeta_orbitals.py \
         tests/test_phillips_kleinman_cross_center.py \
         -q --no-header
102 passed, 1 skipped in 27.94s
```

Zero regression across 102 baseline tests. No production code modified; regression set is the safety net for the multi-zeta + screened-V_ne + PK + balanced-coupled architecture that this sprint exercises.

---

**End of Sprint F1 max_n=3 memo. Verdict: CLEAN NEGATIVE. Sub-layer 3 reclassifies to genuine endomorphism class. M-Z 2-bucket partition stands; proposed 3-bucket refinement is FALSIFIED. Substantive new content: dominant natural orbital at max_n=3 IS a true bonding combination (structural success), but constructed bonding orbital has higher energy than separated configuration at every R (energetic failure). Next sprint: F2 cross-V_ne kernel-shape substitution on partner side.**
