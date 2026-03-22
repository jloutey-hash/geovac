# Phase 5: Residual Structure Analysis

**Date:** 2026-03-22
**Formula:** alpha^3 - K*alpha + 1 = 0, K = pi*(42 + pi^2/6 - 1/40)

## 1. Exact Residual (50-digit precision)

| Quantity | Value |
|:---------|:------|
| K | 137.03606441448154121... |
| alpha_formula^{-1} | 137.03601116313640854... |
| alpha_formula | 0.0072973519260534825... |

**Against CODATA 2018:** alpha^{-1} = 137.035999084(21)

| Quantity | Value |
|:---------|:------|
| delta(alpha^{-1}) | **-1.2079 x 10^{-5}** (formula overshoots) |
| delta(alpha) | +6.43 x 10^{-10} |
| relative error | **8.81 x 10^{-8}** |
| sigma distance | 575 sigma |

The formula predicts alpha^{-1} = 137.03601... while CODATA gives 137.03600...
The sign is negative: the formula's K is slightly too large.

## 2. QED Radiative Corrections — ALL NEGATIVE

No QED correction term matches the residual at useful precision:

| Candidate | Ratio to delta_inv | Assessment |
|:----------|:-------------------|:-----------|
| alpha/(2*pi) [Schwinger] | 96.15 | No match |
| alpha^2/(2*pi) | 0.702 | Nearest, but 30% off |
| alpha/pi | 192.3 | No match |
| alpha^2/pi | 1.403 | ~sqrt(2)? Numerology |
| alpha^3 | 0.0322 | ~1/31, no match |
| (alpha/pi)^2 | 0.447 | No match |

**Conclusion:** The residual does not correspond to any standard QED perturbative correction. This is expected: at q^2 = 0 (Thomson limit), alpha is already the fully renormalized physical coupling. The formula targets the infrared value, not a running coupling.

## 3. Spectral Invariant Corrections — ALL NEGATIVE

No combination of Hopf bundle spectral data (B, F, Delta, K) matches delta_inv:

| Candidate | Ratio to delta_inv |
|:----------|:-------------------|
| 1/(K*B) | 0.070 |
| Delta/K | 0.066 |
| Delta^2 | 0.019 |
| 1/K^2 | 0.227 |

All ratios are far from unity or simple fractions. The residual is not expressible as a simple function of the formula's own ingredients.

## 4. Mathematical Constants Search

**Best match:** pi^4 / 8064243 matches |delta_inv| to **2 x 10^{-9}** relative precision.

However, 8064243 = 3^2 x 11 x 81457 has no obvious spectral or physical meaning. With ~10^7 candidates tested, a match at 10^{-9} level is not statistically anomalous (expected best match ~10^{-7} from birthday-style statistics, so this is mildly notable but not compelling).

No match at <0.1% involves only recognized mathematical constants with small denominators.

## 5. Self-Consistency: Epsilon Correction to K

To match CODATA 2018 exactly:

| Quantity | Value |
|:---------|:------|
| epsilon | 1.2079 x 10^{-5} |
| K_corrected | 137.036076... |
| epsilon/K | 8.81 x 10^{-8} |
| 1/epsilon | 82787.4... |

**Does epsilon match any spectral quantity?**

| Candidate | epsilon/candidate | Assessment |
|:----------|:-----------------|:-----------|
| Delta^2 = 1/1600 | 0.0193 | No |
| F*Delta = pi^2/240 | 2.94 x 10^{-4} | No |
| 1/B = 1/42 | 5.07 x 10^{-4} | No |
| alpha^2/pi | 0.713 | Nearest, but 29% off |
| Schwinger alpha/(2*pi) | 0.0104 | No |
| spectral det gap 0.04276 | 2.82 x 10^{-4} | No |

1/epsilon ~ 82787 is not a recognizable integer in any spectral context.

**Conclusion:** The correction needed is ~8.8 x 10^{-8} of K itself, with no identifiable structure.

## 6. CODATA Vintage Comparison

| Vintage | alpha^{-1} | delta_inv (x 10^{-5}) | rel_err (x 10^{-8}) |
|:--------|:-----------|:---------------------|:--------------------|
| 2014 | 137.035999139(31) | -1.2024 | 8.77 |
| 2018 | 137.035999084(21) | -1.2079 | 8.81 |
| **2022** | **137.035999177(13)** | **-1.1986** | **8.75** |

**Trend: NON-MONOTONIC.** The 2018 value moved slightly away from the formula (delta grew from 1.2024 to 1.2079 x 10^{-5}), but the 2022 value bounced back closer (1.1986 x 10^{-5}). The movement is within ~1% of the residual itself — measurement noise at this precision.

All three vintages give relative error in the narrow band **[8.75, 8.81] x 10^{-8}**. The residual is stable across measurement improvements. This is consistent with a genuine ~10^{-7} level gap intrinsic to the formula.

**Note on 2022:** The CODATA 2022 value uses the TAMU Cs-133 recoil measurement which shifted alpha upward. The sigma distance actually *increases* from 575 to 922 because uncertainty shrank faster than the residual.

## Overall Assessment

**The residual is intrinsic to the formula and does not match any known physics correction.**

1. It is NOT a QED radiative correction (wrong order of magnitude for all terms)
2. It is NOT expressible in terms of the formula's own spectral ingredients
3. It is NOT a recognizable mathematical constant
4. It is NOT trending toward zero as measurements improve
5. The needed K-correction (epsilon ~ 1.21 x 10^{-5}) has no identifiable structure

The paper's current statement (Sec VIII.E, open question 3) that the residual is "intrinsic to the formula" is **confirmed and strengthened** by this analysis. No modification to the paper is warranted.

The residual likely reflects the limit of what a purely topological formula (using only spectral invariants of S^1, S^2, S^3) can capture. The missing ~10^{-7} may encode dynamics (vacuum fluctuations, higher-loop topology) that the static Hopf bundle geometry cannot access. This remains an open question.
