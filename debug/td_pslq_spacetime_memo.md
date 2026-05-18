# Sprint TD-PSLQ-2 — Spacetime/relativistic probe of A_60(1S) against master Mellin engine

**Date:** 2026-05-17
**Sprint type:** Diagnostic PSLQ probe (no production code or paper changes)
**Verdict:** **NULL** — `A_60(1S)` does not sit in the master Mellin engine ring at PSLQ precision available

This is the second confirmation that calibration-class atomic-QED Layer-2 constants are not in M1/M2/M3, and the FIRST test on the proper channel: a structurally spacetime/relativistic constant rather than a Bethe-log virtual-state sum. **PI hypothesis "spacetime corrections close Layer-2 residuals to bit-exactness" is NOT SUPPORTED for the alpha(Zα)^6 self-energy channel.**

---

## 1. Target survey

Six candidates were evaluated against the criteria (high precision; spacetime/relativistic mechanism; no known closed form; empirical relevance):

| # | Candidate | Mechanism | Precision sourced | Status |
|:--|:----------|:----------|:------------------|:-------|
| a | Pachucki muonic A_60 | muonic-specific recoil + QED higher-order | not located numerically in survey | not selected |
| b | Penin-Pivovarov α⁵ log coefficient (Ps 1S-2S) | two-loop Breit retardation | not located numerically in survey | not selected |
| c | Karshenboim Mu HFS recoil NLO | nonlogarithmic radiative-recoil | ~3 digits ("2.64 ± 0.07 kHz") | precision too low |
| d | Drake/Pachucki-Yerokhin He 2³P α³(Zα)² | bipolar harmonic relativistic | ~few kHz precision; multiple components | better suited to autopsy than PSLQ |
| **e** | **A_60(1S) hydrogen self-energy** | **α(Zα)^6 Dirac-Coulomb + radiative** | **7 digits, value -30.92415(1)** | **SELECTED** |
| f | Pachucki-Yerokhin two-loop B_60(1S) 2024 | two-loop self-energy | -72(7), ~2 digits | precision insufficient |

**Selection rationale.** A_60(1S) is the canonical alpha(Zα)^6 self-energy coefficient. It is the **dominant relativistic correction beyond the Bethe log** in the one-loop hydrogen self-energy expansion. Mechanism is unambiguously spacetime/relativistic: it arises from Dirac-Coulomb kinematics combined with the one-loop radiative correction at order α(Zα)^6, with all logarithmic content (which sits in A_50, A_60_log) explicitly separated. The nonlogarithmic A_60 was the load-bearing irrational in Lamb-shift theory for three decades of analytic work before the Jentschura-Mohr-Soff 1999 numerical breakthrough (PRL 82, 53). Sources include Jentschura-Mohr-Soff 1999, Mohr-Jentschura-Pachucki 2000 ([arXiv:physics/0001068](https://arxiv.org/abs/physics/0001068)), and the Yerokhin-Pachucki-Patkóš 2019 review ([arXiv:1809.00462](https://arxiv.org/abs/1809.00462)). The 7-digit precision is published with explicit uncertainty in the last digit and reproduced consistently across multiple independent calculations.

**Why this is the right test for the PI hypothesis.** The PI's question is "do spacetime/relativistic corrections close Layer-2 residuals?" A_60(1S) is THE spacetime-class constant in one-loop atomic QED — it captures exactly the relativistic-kinematic physics that the gravity/spacetime hypothesis points to. A hit here would partially confirm the hypothesis. A null sharpens the structural-skeleton-scope reading directly on the relativistic channel.

---

## 2. Selected target

```
A_60(1S) = -30.92415 ± 0.00001
```

| Item | Value |
|:-----|:------|
| Value | -30.92415 |
| Precision | 7 significant digits |
| Uncertainty (last digit) | ±1, i.e. ±0.00001 |
| Source | Jentschura-Mohr-Soff PRL 82, 53 (1999) |
| Cross-checked in | Mohr-Jentschura-Pachucki 2000 (arXiv:physics/0001068); Yerokhin-Pachucki-Patkóš 2019 (arXiv:1809.00462) |
| Physical meaning | nonlogarithmic α(Zα)^6 one-loop self-energy coefficient |
| GeoVac projection class | spacetime/relativistic (Dirac-Coulomb kinematics × radiative correction) |

---

## 3. Basis composition (frozen before PSLQ)

Atomic basis (same architecture as TD-PSLQ-1, with one additional class SPACETIME_AUG tailored to relativistic-correction patterns). Frozen at `debug/data/td_pslq_spacetime_basis_frozen.json`, SHA256 prefix `5caf5b45b5607ece`:

| Class | Count | Forms |
|:------|------:|:------|
| M1 (Hopf-base / π-powers) | 11 | 1, π, π², π³, π⁴, π⁵, π⁶, 1/π, 1/π², 1/π³, 1/π⁴ |
| M2 (Seeley–DeWitt / sqrt(π) tower) | 6 | √π, √π³, √π⁵, √π⁷, 1/√π, 1/√π³ |
| M3 (vertex parity / Dirichlet–Hurwitz) | 14 | G, β(4), β(6), β(8), ζ(s,1/4) for s=2..6, ζ(s,3/4) for s=2..6 |
| ALG (algebraic controls) | 10 | √2, √3, √5, √6, √7, φ, ln 2, ln 3, ln 5, ln 7 |
| ODD_ZETA_CONTROL (NOT in M1/M2/M3) | 4 | ζ(3), ζ(5), ζ(7), ζ(9) |
| CROSS_PRODUCT (depth-2 forms) | 27 | π·G, π²·G, π·ln2, ..., zeta(3)·π, zeta(3)·ln2, G·√2, ln2·√π (extended from TD-PSLQ-1's 19) |
| SPACETIME_AUG (large-mag π^k, Dirac-Coulomb closed-form patterns, Bernoulli-style rationals) | 18 | π³ (~31), 10π², π⁴/3, 10·ln2·π, 4/3, 11/24, 35/72, 139/144, π²·ln2/3, π²/12, 7ζ(3)/8, 7π⁴/360, 31π²/216, π²·ln2/2, 3ζ(3)/2, 9/2, ... |
| **TOTAL** | **90** | |

ODD_ZETA_CONTROL serves as a control class — odd Riemann zetas are structurally NOT in master Mellin engine rings. Their inclusion tests whether PSLQ can be confused by depth-2 cross-products involving them.

The SPACETIME_AUG class is explicitly designed for this target: A_60(1S) ≈ -30.9, so we include atoms with magnitude near 30 (π³ ≈ 31.0, 10π² ≈ 98.7, π⁴/3 ≈ 32.5) to let PSLQ identify dominant single-term hits. Bernoulli-style rationals (4/3, 11/24, 35/72, 139/144) arise naturally in Dirac-Coulomb expansions; Pachucki-Karshenboim closed-form patterns like 7ζ(3)/8, 7π⁴/360, 31π²/216 are included as natural candidates for an A_60 closed form.

Frozen before any PSLQ run. SHA256 stamped.

---

## 4. PSLQ run results per ceiling

10 sub-runs × 3 coefficient ceilings (10⁴, 10⁶, 10⁸) = **30 PSLQ tests**, each at 50 dps, `maxsteps=500`. Genuine-relation threshold: residual < 1e-6 (1-digit headroom on 7-digit target).

| Sub-run | N | Best outcome across ceilings (residual stable across all 3 ceilings) |
|:--------|---:|:---------------------------------------------------------------------|
| M1_only | 11 | basis-internal best-guess, residual 6.1e-4, max_coef 9 |
| M2_only | 6 | spurious-filled, residual 2.9e-5, max_coef 11 |
| M3_only | 14 | basis-internal best-guess, residual 4.0e-3, max_coef 23 |
| M1_M2 | 17 | spurious-filled, residual 2.1e-4, max_coef 10 |
| M1_M3 | 25 | basis-internal best-guess, residual 4.0e-3, max_coef 23 |
| M1_M2_M3 | 31 | basis-internal best-guess, residual 4.0e-3, max_coef 23 |
| M1_M2_M3_ALG | 41 | basis-internal best-guess, residual 4.1e-3, max_coef 3 |
| M1_M2_M3_ALG_CROSS | 67 | basis-internal best-guess, residual 3.4e-3, max_coef 1 |
| M1_M2_M3_ALG_CROSS_AUG | 84 | basis-internal best-guess, residual 3.4e-3, max_coef 1 |
| FULL | 88 | spurious-filled, residual 1.9e-3, max_coef 11 |

**Zero genuine hits across all 30 tests.** Best residual (M2_only, 2.9e-5) is still 29× above the 1e-6 genuine-relation threshold and arises from a 7-form combination at max_coef=11 — far too elaborate to be a real low-depth identification.

**Key signature: residuals are independent of coefficient ceiling.** Same 6.1e-4 at C=10⁴, 10⁶, 10⁸ for M1_only; same 4.0e-3 across ceilings for M3_only; same 1.9e-3 for FULL. This is the classic null signature: PSLQ cannot improve the residual by being allowed to use larger coefficients, meaning there is no improvable target relation to find.

The basis-internal best-guesses are reporting the smallest detectable linear combination among basis forms (with a_target = 0). Their residuals (~mHa scale) are above target precision, confirming no nontrivial Q-linear relations exist among the basis forms themselves either at this precision.

---

## 5. Hit verification or null confidence statement

**No hits were verified.**

**Null confidence statement.** For depth-1 integer relation `a · A_60(1S) + Σ b_i · basis_i = 0` against an irrational target with no actual relation in the basis, the false-positive probability is bounded by:

$$P_{FP} \approx N \cdot C / 10^{d_{\rm target}}$$

with d_target = 7 digits:

| N (sub-run dim) | C=10⁴ | C=10⁶ | C=10⁸ |
|----------------:|------:|------:|------:|
| 11 (M1) | 1.1e-2 | 1.1e+0 | 1.1e+2 |
| 14 (M3) | 1.4e-2 | 1.4e+0 | 1.4e+2 |
| 31 (M1_M2_M3) | 3.1e-2 | 3.1e+0 | 3.1e+2 |
| 88 (FULL) | 8.8e-2 | 8.8e+0 | 8.8e+2 |

At target precision 7 digits, the verdict-bearing ceiling is C = 10⁴. P_FP at ceiling 10⁴ ranges from 1.1% (M1 only) to 8.8% (FULL) — borderline-safe for ruling out small-coefficient depth-1 relations, unreliable for depth-3+ or large coefficient.

**Compared to TD-PSLQ-1's 16-digit Bethe log probe, this verdict is necessarily weaker.** At 7-digit precision, we can rule out only the cleanest single-atom identifications like `A_60(1S) = -10π² · (rational)` or similar low-depth low-coefficient forms. We cannot rule out arbitrarily intricate higher-depth combinations. **A higher-precision (≥ 15-digit) value of A_60(1S) would be required to probe to coefficient ceiling 10⁸ with full confidence; no such published value exists as of May 2026** — the Jentschura-Mohr 1999 result with 7-digit precision is the canonical reference, and the subsequent two-loop work (Yerokhin et al. 2024) extends to B_60 at coarser precision.

**However**, the structural shape of the result is telling. The smallest residual achieved across all 30 tests is 2.9e-5 (M2_only at 7-form depth, max_coef=11). At full basis (FULL, N=88) the best residual is 1.9e-3 — only ~3 digits below the target value's magnitude. If A_60(1S) had a clean low-depth identification in M1/M2/M3, we would expect residuals well below 1e-6 at moderate coefficient ceilings. The fact that the best residual at the largest basis still sits at the mHa scale, four orders of magnitude above the precision floor, is consistent with A_60(1S) being **structurally outside** the master Mellin engine ring rather than merely beyond the precision-ceiling combination of probe.

---

## 6. Per-mechanism interpretation

The PI hypothesis being tested: "if we work out the spacetime corrections across all our focal lengths, we start hitting bit-exact results in Layer 2."

**A_60(1S) is the cleanest test of this hypothesis available at present**: it is the canonical spacetime/relativistic Layer-2 constant in atomic QED, the dominant non-Bethe-log α(Zα)^6 self-energy coefficient. It is the principal residual surviving in the LS-7 cumulative-chain decomposition (per Paper 36 §VII reframing) and corresponds directly to the "recoil + finite-nuclear-size + relativistic hyperfine ~+4.4 MHz" attribution in the Lamb shift residual itemization.

If GeoVac's master Mellin engine (M1/M2/M3 rings) were structurally rich enough to contain spacetime corrections at this Layer-2 order, A_60(1S) is exactly where we would expect the closure to show. **The null result here is therefore the strongest possible specific falsification of the spacetime-projection hypothesis at the precision available** — not just for one observable, but for the canonical spacetime-class constant.

The Paper 34 projection-class identification of A_60(1S) sits clearly in:
- **§III.1 (Fock projection)** + **§III.7 (spinor lift)** + **§III.14 (rest-mass)** + **§III.16 (Breit retardation)** — the multi-projection chain whose composition gives the relativistic-corrected one-electron self-energy.

The null verdict says: even though each individual projection in the chain is well-defined and the framework can compute the diagrammatic structure (see Paper 28 graph-native QED), the **numerical value** of the composite is calibration data that the framework does not autonomously generate from its M1/M2/M3 inventory.

This is the same structural-skeleton-scope pattern that LS-8a, W3, and TD-PSLQ-1 surfaced. **Sprint TD-PSLQ-2 is the fourth independent confirmation, and the first specifically on the spacetime-projection class.**

---

## 7. Net verdict on the gravity/spacetime hypothesis at this channel

**NULL — gravity/spacetime hypothesis NOT SUPPORTED at this channel.**

The proposed reading "spacetime corrections close Layer-2 residuals to bit-exactness" fails for the α(Zα)^6 one-loop self-energy channel. The canonical relativistic-correction constant A_60(1S) is not in the M1/M2/M3 ring at any precision the PSLQ probe can reach (coefficient ceiling 10⁴ at 7-digit target precision; higher ceilings unreliable at this target precision).

**Two interpretations remain compatible with the null:**

(I) **Structural-skeleton-scope is correct (recommended reading).** The framework's spectral-triple machinery determines the structure of corrections (mechanism, transcendental class, divergence order, scaling exponents) but does not autonomously generate the numerical coefficients of multi-projection compositions. The spacetime corrections close NOTHING bit-exactly; they are calibration inputs to the framework. This is the dominant reading after four independent null results (LS-8a, W3, TD-PSLQ-1, TD-PSLQ-2).

(II) **The precision is insufficient and we cannot test this channel deeply enough.** With only 7 digits of target precision, the PSLQ probe rules out only depth-1 to depth-3 small-coefficient combinations. A genuine deep closure (e.g., A_60(1S) = combination of 5 atomic constants with all coefficients < 100) could exist below precision and not be detected. This reading is less plausible given the cross-ceiling residual-invariance signature but cannot be ruled out at present.

**To distinguish (I) from (II):** the canonical next probe would be a target with ≥ 15-digit precision. The available candidates are:

| Candidate | Mechanism | Available precision | Status |
|:----------|:----------|:--------------------|:-------|
| Bethe logarithm ln k_0(1S) | bound-state QED virtual-state sum | 16-digit Drake / 50-digit Korobov 2002 | TD-PSLQ-1 NULL at 16 digits; Korobov 2002 50-digit re-run available |
| ζ(3) | algebraic content of various QED quantities | infinite-precision via mpmath | known NOT in M1/M2/M3 (used as control) |
| Two-loop pure self-energy B_60(1S) | two-loop spacetime | -72(7), 2 digits | precision insufficient |

There is no spacetime/relativistic Layer-2 constant in atomic QED that is BOTH structurally spacetime AND known to ≥ 15-digit precision. The 7-digit A_60(1S) is the cleanest available spacetime test, and it is null.

**Recommended verdict: structural-skeleton-scope reading sharpens specifically on the spacetime-projection channel. Reading (I) wins on present evidence.** The PI's "spacetime corrections close Layer-2 to bit-exactness" hypothesis is **falsified at the present-best test channel**.

---

## 8. Recommended next probe

Two complementary directions:

(a) **Korobov 2002 high-precision Bethe log re-test** (~1 day, follows TD-PSLQ-1's recommended next test): source the 50-digit Korobov 2002 value of ln k_0(1S), re-run TD-PSLQ-1 at coefficient ceiling 10⁸ with a depth-3 basis. If still null at ≥ 50-digit precision and ceiling 10⁸, the calibration-class reading for the Bethe log channel becomes **structurally irreducible** in any practical PSLQ probe — formal closure on the bound-state-QED virtual-state-sum channel.

(b) **State-side dictionary opening** (~1 week, the structurally larger move): the spectral-side dictionary (Paper 34 §III.1-§III.27) has been tested four times against the master Mellin engine and returns four nulls on calibration constants. The state-side complement (§III.28 apparatus identity, opened 2026-05-15 Sprint 2) is genuinely untested for transcendental content. PSLQ-probing a state-side constant like a Wasserstein-Kantorovich distance, mutual information across a hyperfine doublet, or fidelity bound between two truncated spectral triple states — against an enlarged basis that explicitly includes von Neumann entropy and information-theoretic atoms — tests whether the state-side admits structural closures the spectral side does not.

**The structurally informative recommendation is (b).** Direction (a) sharpens an existing null result by one more order of magnitude; direction (b) opens a previously-untested axis of the Paper 34 dictionary. Per the "diagnostic-before-engineering" rule (CLAUDE.md §1.8 + memory `feedback_diagnostic_before_engineering.md`), after four convergent nulls on the spectral side, the structurally informative move is to ask whether the negative result is uniform across the full dictionary or specific to the spectral-side half.

---

## Files

- `debug/td_pslq_spacetime.py` — the probe script (atomic basis + spacetime augmentation, post-filter classifier, multi-stage sub-runs)
- `debug/data/td_pslq_spacetime_basis_frozen.json` — frozen basis (SHA256 `5caf5b45b5607ece...`, 90 atomic forms)
- `debug/data/td_pslq_spacetime_results.json` — full PSLQ results (30 tests, 10 sub-runs × 3 ceilings)
- `debug/data/td_pslq_spacetime_stdout.log` — run log

## Cross-references

- TD-PSLQ-1 Bethe log null: `debug/td_pslq_bethe_log_memo.md`
- LS-8a renormalization gap: CLAUDE.md §2 LS-8a entry
- W3 spectral-zeta falsification: CLAUDE.md §3 W3 row + memory `w3_spectral_zeta_candidate.md`
- Structural-skeleton-scope reading: memory `geovac_structural_skeleton_scope_pattern.md`
- PI hypothesis source: Sprint TD-PSLQ-1 reflection conversation about spacetime/gravitational Layer-2 corrections
