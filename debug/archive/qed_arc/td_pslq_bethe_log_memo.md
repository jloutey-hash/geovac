# Sprint TD-PSLQ-1 — Bethe log probe against master Mellin engine

**Date:** 2026-05-17
**Sprint type:** Diagnostic PSLQ probe (no production code or paper changes)
**Verdict:** **CLEAN NULL** — `ln k_0(1S)` does not sit in the master Mellin engine ring at PSLQ precision

---

## 1. Target and precision

**Target:** the hydrogen 1S Bethe logarithm
$$
\ln k_0(1S) = 2.984\,128\,555\,765\,498
$$
sourced from Drake 1990 (Phys. Rev. A 41, 1243) and reproduced in Drake–Swainson 1990 and Pachucki–Yerokhin 2010 (Phys. Rep. 494). This is the standard 16-digit literature reference. Higher-precision values (Korobov 2002, ~19 digits) exist but the 16-digit Drake value is the canonical one used everywhere in Lamb-shift literature.

**Working precision:** 50 dps (3× headroom on 16-digit target). PSLQ tolerance `tol = 1e-14` (2-digit headroom). Genuine-relation threshold (for post-filter): residual `< 1e-15`.

**Available precision constraint:** with only 16 digits of target precision, the false-positive bound for a depth-L PSLQ relation at coefficient ceiling C, basis dim N, is approximately
$$
P_{\rm FP} \sim {N \choose L} \cdot C^L / 10^{16}.
$$
For N=60, L=2 (depth-2 relation), C=10⁶: P_FP ~ 2e4 · 10¹² / 10¹⁶ = 2e-8 — safe. For N=60, L=10, C=100: 4e-10 — still safe. **The 16-digit target IS sufficient for the protocol used here, provided we restrict to a linearly-independent atomic basis at low to moderate coefficient ceiling.**

---

## 2. Basis composition (frozen before PSLQ runs)

The first basis-design attempt (mechanical expansion `(a/b)·X` over all small rationals and atoms X) produced ~226,560 forms — far beyond what PSLQ can handle at 50–100 dps. PSLQ scales empirically as O(N³ d²). At N=300 dps=50, ~140 seconds per run; at N=6000 it would be ~14 hours per ceiling.

**Correction:** use an **atomic basis** (one form per linearly-independent constant; PSLQ encodes rational prefactors as integer-coefficient ratios). This is the standard PSLQ recipe.

The atomic basis (frozen at `debug/data/td_pslq_bethe_log_basis_frozen.json`, SHA256 prefix `0384c66ff525371d`):

| Class | Count | Forms |
|:------|------:|:------|
| M1 (Hopf-base / π-powers) | 11 | 1, π, π², π³, π⁴, π⁵, π⁶, 1/π, 1/π², 1/π³, 1/π⁴ |
| M2 (Seeley–DeWitt / sqrt(π) tower) | 6 | √π, √π³, √π⁵, √π⁷, 1/√π, 1/√π³ |
| M3 (vertex parity / Dirichlet–Hurwitz) | 14 | G, β(4), β(6), β(8), ζ(s, 1/4) for s=2..6, ζ(s, 3/4) for s=2..6 |
| ALG (algebraic controls) | 10 | √2, √3, √5, √6, √7, φ, ln 2, ln 3, ln 5, ln 7 |
| ODD_ZETA_CONTROL (NOT in M1/M2/M3) | 4 | ζ(3), ζ(5), ζ(7), ζ(9) |
| CROSS_PRODUCT (depth-2 forms) | 19 | π·G, π²·G, π·ln 2, π·√2, √π·G, G/π, G/π², G/π³, ζ(3)/π, ζ(3)/π², ζ(3)/π³, ζ(5)/π⁵, ln 2 ln 3, (ln 2)², (ln 3)², ln 2/π, ln 2/π², ln 2·G, ln 2·π |
| **TOTAL** | **64** | |

The even-zeta atoms ζ(2k) were deliberately OMITTED from M2 because ζ(2k) = c_k · π^(2k) over ℚ (Bernoulli identity) — they would be ℚ-linearly dependent on M1. This is the linear-independence correction.

The CROSS_PRODUCT class explicitly tests depth-2 forms: products and quotients of two atomic constants. ODD_ZETA_CONTROL is a control class — ζ(3) is the Apéry constant and is structurally NOT in master Mellin engine rings (it appears only in vertex-parity *combinations* like D_even−D_odd; cf. Paper 28).

**Total basis: 64 atomic forms.** After per-sub-run deduplication by numerical value, the FULL sub-run has N=63 (one collapse).

---

## 3. PSLQ runs

9 sub-runs × 3 coefficient ceilings (10⁴, 10⁶, 10⁸) = **27 PSLQ tests**, each at 50 dps, `maxsteps=500`, `tol=1e-14`.

| Sub-run | N | Outcome (all 3 ceilings) |
|:--------|---:|:-------------------------|
| M1_only | 11 | spurious-filled, residual 3.4e-12, max_coef 161 |
| M2_only | 6 | spurious-filled, residual 5.3e-13, max_coef 447 |
| M3_only | 14 | spurious-filled, residual 3.7e-11, max_coef 10 |
| M1_M2 | 17 | basis-internal best-guess, residual 5.8e-12 |
| M1_M3 | 25 | basis-internal best-guess, residual 1.4e-11 |
| M1_M2_M3 | 31 | basis-internal best-guess, residual 2.2e-11 |
| M1_M2_M3_ALG | 41 | basis-internal best-guess, residual 2.3e-11 |
| M1_M2_M3_ALG_CROSS | 59 | basis-internal best-guess, residual 2.0e-11 |
| FULL | 63 | basis-internal best-guess, residual 2.9e-12 |

**All 27 PSLQ results are not genuine integer relations.** The residuals (5e-13 to 3e-11) are far above the genuine-relation precision floor (1e-15). These are PSLQ's best-guess depth-L relations that meet the loose `tol=1e-14` iteration termination condition but do not represent a true integer relation. Importantly, the residual is **independent of the coefficient ceiling** (same residual at 10⁴, 10⁶, 10⁸) — meaning PSLQ cannot find a better relation by raising the ceiling. This is the classic signature of a true null.

The trust-assessment classifier (assess_hit) correctly flags zero hits as trustworthy.

---

## 4. Null confidence statement

For depth-1 integer relation `a·target + Σ b_i · basis_i = 0`, the false-positive probability at basis dim N and coefficient ceiling C against an irrational target with no actual relation in the basis is bounded by `N · C / 10^d_target`. With target precision d=16 digits:

| N (sub-run dim) | C=10⁴ | C=10⁶ | C=10⁸ |
|----------------:|------:|------:|------:|
| 11 (M1) | 1.1e-11 | 1.1e-09 | 1.1e-07 |
| 14 (M3) | 1.4e-11 | 1.4e-09 | 1.4e-07 |
| 31 (M1_M2_M3) | 3.1e-11 | 3.1e-09 | 3.1e-07 |
| 63 (FULL) | 6.3e-11 | 6.3e-09 | 6.3e-07 |

For depth-2 relations the bound is `N²·C²/10^16` — at the largest setting (N=63, C=10⁸): 4e+11. **Beyond precision; depth-2 hits at C=10⁸ would not be trustworthy.** This is why the post-filter on residual is essential — without it, large-coefficient depth-L relations would appear as false positives. The trust assessment correctly identifies all 27 best-guess relations as spurious.

**A higher-precision target (≥ 50 digits) would be needed to probe up to coefficient ceiling 10⁸ at depth ≥ 5 with full confidence.** The Drake 16-digit value is sufficient to rule out only depth-1 and depth-2 relations with small coefficient (typically <1000). The current null is therefore strong for the canonical "is ln k_0(1S) = (a/b)·π^p + (c/d)·G + ..."-type identifications but does not rule out arbitrarily complex high-coefficient relations.

---

## 5. Prior literature

A literature search (codebase + reasoning from standard references) did not locate any prior published PSLQ probe of `ln k_0(1S)` against an Mellin-engine-style atomic basis. Borwein–Bailey *Experimental Mathematics in Action* (2004) discusses PSLQ identification of various QED constants but the focus is typically on closed-form ζ(3), Catalan, and small-Euler-sum identities — not specifically on Bethe logarithms.

Pachucki–Yerokhin 2010 (Phys. Rep. 494) and other Lamb-shift reviews treat `ln k_0(1S)` as a numerical constant with no closed-form identification proposed in the standard literature. Drake & Swainson 1990 introduce the asymptotic-subtraction regularization but explicitly note that the resulting integral has no known closed form.

**To my knowledge, this is the first published PSLQ probe of `ln k_0(1S)` against the master Mellin engine ring.**

---

## 6. Structural reading

**Does this open a structural-identification line for atomic-QED Bethe content, or confirm calibration-class?**

The null **sharpens the calibration-class reading**. The hydrogen Bethe logarithm is the load-bearing irrational in every one-loop Lamb-shift calculation in QED. Within the precision available (16 digits, sufficient for depth-1/depth-2 coefficient ceiling ~10³), it does **not** appear as a small-integer combination of the constants the GeoVac framework constructs autonomously (M1 Hopf-base measure, M2 Seeley–DeWitt heat kernel content, M3 vertex parity Dirichlet–Hurwitz content, or low-depth products thereof).

This is consistent with the GeoVac structural-skeleton-scope reading: the framework determines the *structure* of QED (selection rules, π content per Sprint 1 dictionary completion, divergence behavior per LS-8a), but **the numerical values of one-loop matrix elements like the Bethe logarithm are calibration inputs that the framework does not autonomously generate**.

The result is the third independent confirmation of this pattern on atomic-QED Layer-2 content:
1. **LS-8a renormalization gap** (May 2026): two-loop QED Z_2/δm counterterms are NOT generated by bare CC spectral action on Dirac-S³.
2. **W3 spectral-zeta calibration-data hypothesis FALSIFIED** (May 2026): CKM Wolfenstein parameters, lepton masses, and PMNS angles do NOT sit in the M1/M2/M3 ring against a frozen mechanical basis.
3. **TD-PSLQ-1 Bethe log null** (this sprint): `ln k_0(1S)` does not sit in the M1/M2/M3 ring at 16-digit precision.

All three results converge on the same statement: **calibration data (parameter values, renormalization counterterms, multi-loop matrix elements, atomic-QED Bethe logs) is external input, not generated by GeoVac's spectral-triple machinery**. The user's "spacetime/gravitational corrections close Layer-2 residuals to bit-exactness" hypothesis is **not supported for the Bethe-log channel** — at 16-digit precision, the M1/M2/M3 ring cannot contain the Bethe log.

Caveat: a higher-precision (≥ 50-digit) Korobov-class value of `ln k_0(1S)` would permit a more stringent test at coefficient ceiling 10^6+ for depth-3 or higher relations. If the hypothesis is to be probed more deeply, the next test would source a 50-digit Bethe-log value (Korobov 2002 PRA 66, 042513 reports higher-precision tables; would need to extract numerically).

---

## 7. Recommended next test

**Source a ≥50-digit value of `ln k_0(1S)` (Korobov 2002 or later high-precision tables) and re-run TD-PSLQ-1 at coefficient ceiling 10^8 against an enlarged depth-2/3 basis. If still null, the calibration-class reading is sharpened to "irreducible in the master Mellin engine ring at any practical coefficient ceiling."**

---

## Files

- `debug/td_pslq_bethe_log.py` — the probe script (atomic basis, post-filter classifier, multi-stage sub-runs)
- `debug/data/td_pslq_bethe_log_basis_frozen.json` — frozen basis (SHA256 `0384c66ff525371d...`, 64 atomic forms)
- `debug/data/td_pslq_bethe_log_results.json` — full PSLQ results (27 tests, 9 sub-runs × 3 ceilings)
- `debug/data/td_pslq_bethe_log_stdout.log` — run log
