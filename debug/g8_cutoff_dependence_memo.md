# Sprint G8 — Cutoff-function dependence of $G_{\rm eff}$ and $\Lambda_{cc}$

**Date:** 2026-05-28
**Path:** Gravity arc consolidation. Generalizes G7's Gaussian-cutoff Newton constant to arbitrary cutoff function via Mellin transform $\phi(s) = \int_0^\infty f(x) x^{s-1} dx$.
**Verdict:** **POSITIVE-CALIBRATION-CLARIFICATION.** All gravity predictions ($G_{\rm eff}$, $\Lambda_{cc}$, $R_{\rm crit}$) are cutoff-dependent; the cutoff function is Class 1 calibration data per `memory/external_input_three_class_partition.md`. Consistent with standard CC.

## 1. Setup

G7 derived $G_{\rm eff} = 6\pi/\Lambda^2$ and $\Lambda_{cc} = 6\Lambda^2$ for the Gaussian cutoff $f(x) = e^{-x}$. The question: how do these change for arbitrary $f$?

The Mellin transform $\phi(s) = \int_0^\infty f(x) x^{s-1} dx$ captures the cutoff's spectral content. From the G2 spectral action on $S^3_R \times S^1_\beta$:
$$S(R, \beta, \Lambda) = \phi(2) \cdot \frac{\beta R^3}{4} \Lambda^4 - \phi(1) \cdot \frac{\beta R}{8} \Lambda^2$$

## 2. General formulas

Matching to Einstein-Hilbert plus cosmological constant:
$$S_{\rm EH} = -\frac{1}{16\pi G_{\rm eff}}\int R_{\rm scalar} \sqrt{g} + \frac{\Lambda_{cc}}{8\pi G_{\rm eff}}\int \sqrt{g}$$

with $\text{Vol}(S^3_R \times S^1_\beta) = 2\pi^2 R^3 \beta$ and $\int R_{\rm scalar} \sqrt{g} = 12\pi^2 R \beta$:

$$\boxed{\;G_{\rm eff}(f, \Lambda) = \frac{6\pi}{\phi(1)\,\Lambda^2}\;}$$

$$\boxed{\;\Lambda_{cc}(f, \Lambda) = \frac{6\,\phi(2)}{\phi(1)}\,\Lambda^2\;}$$

$$\boxed{\;R_{\rm crit}\,\Lambda = \sqrt{\frac{\phi(1)}{6\,\phi(2)}}\;}$$

## 3. Three specific cutoffs

| Cutoff $f(x)$ | $\phi(1)$ | $\phi(2)$ | $G_{\rm eff} \Lambda^2$ | $\Lambda_{cc}/\Lambda^2$ | $R_{\rm crit} \Lambda$ |
|---|---|---|---|---|---|
| Gaussian $e^{-x}$ | $1$ | $1$ | $6\pi$ | $6$ | $1/\sqrt{6} \approx 0.408$ |
| Polynomial $e^{-x^2}$ | $\sqrt{\pi}/2$ | $1/2$ | $12\sqrt{\pi}$ | $6/\sqrt{\pi}$ | $\sqrt{\sqrt{\pi}/6} \approx 0.544$ |
| Sharp $\Theta(1-x)$ | $1$ | $1/2$ | $6\pi$ | $3$ | $1/\sqrt{3} \approx 0.577$ |

All values verified symbolically.

## 4. Cutoff-independence analysis

Tested combinations for cutoff-independence:

| Quantity | Form | Cutoff-independent? |
|---|---|---|
| $G_{\rm eff} \cdot \Lambda^2$ | $6\pi/\phi(1)$ | No (depends on $\phi(1)$) |
| $\Lambda_{cc}/\Lambda^2$ | $6\phi(2)/\phi(1)$ | No (depends on ratio) |
| $G_{\rm eff} \cdot \Lambda_{cc}$ | $36\pi \phi(2)/(\phi(1)^2 \Lambda^2)$ | No |
| $(G_{\rm eff}\Lambda^2)/(\Lambda_{cc}/\Lambda^2)$ | $\pi/\phi(2)$ | No |
| $R_{\rm crit} \Lambda$ | $\sqrt{\phi(1)/(6\phi(2))}$ | No |

**Conclusion: NO cutoff-independent gravity prediction exists.** All depend on the cutoff function's $\phi$-moments.

## 5. Class-1 calibration data classification

Per `memory/external_input_three_class_partition.md`, GeoVac's external inputs partition into three structurally distinct classes:

1. **Class 1 — Calibration parameter values** (Yukawas, $Z_2$, $\delta m$, multi-loop QED, Born rule, **cutoff function**)
2. Class 2 — Multi-focal composition (HF-3/4/5)
3. Class 3 — Multi-determinant correlation (W1e)

**The cutoff function $f$ is Class 1 calibration data.** The framework gives proportionality relations but does NOT select $\phi(1), \phi(2)$ from first principles.

## 6. Consistency with standard CC

This is the standard situation in Chamseddine-Connes (1997, 2010) spectral action. The cutoff function $f$ is external input; the framework determines:
- The **structural shape** of $S[D]$ (spectrum, Seeley-DeWitt coefficients)
- The **proportionality** of $G_{\rm eff}$, $\Lambda_{cc}$ to specific $\phi$-moments
- The **functional form** of $G_{\rm eff}(f)$, $\Lambda_{cc}(f)$, $R_{\rm crit}(f)$

The framework does NOT determine:
- The numerical value of $G_{\rm eff}$ or $\Lambda_{cc}$ in absolute units
- The selection of which cutoff function to use

This is the standard structural-skeleton-scope reading sharpened from `memory/geovac_structural_skeleton_scope_pattern.md`.

## 7. Structural-skeleton scope reading

| Layer | Property | Cutoff-dependent? |
|---|---|---|
| **Structural** | Two-term exactness (Paper 28 G1) | No |
| **Structural** | $\zeta_{\rm unit}(-k) = 0$ identity | No |
| **Structural** | $S^3 \times S^1$ factorization (G2) | No |
| **Structural** | Bianchi cancellation in graviton sector (G3) | No |
| **Structural** | Bekenstein-Hawking $S = A/(4G)$ form (G4) | No |
| **Structural** | de Sitter vacuum solution (G7) | No |
| **Structural** | Extremality ratios $R_{\rm crit} \Lambda$ form | No |
| **Numerical** | $G_{\rm eff}$ value | Yes (depends on $\phi(1)$) |
| **Numerical** | $\Lambda_{cc}$ value | Yes (depends on $\phi(2)/\phi(1)$) |
| **Numerical** | $R_{\rm crit}$ value | Yes (depends on ratio) |

**Net: structural findings are cutoff-INDEPENDENT; numerical gravity predictions are cutoff-DEPENDENT.**

## 8. Honest scope

**Reached:**
- General $G_{\rm eff}(f, \Lambda)$, $\Lambda_{cc}(f, \Lambda)$, $R_{\rm crit}(f, \Lambda)$ formulas derived ✓
- Three specific cutoffs computed (Gaussian, polynomial, sharp) ✓
- Cutoff-independence analysis: no simple universal combination ✓
- Class 1 calibration data classification ✓
- Structural vs numerical scope clarification ✓

**Not reached:**
- First-principles selection of cutoff function (open question; structurally external)
- Connection of $\phi$-moments to other framework calibration data (e.g., Yukawas)

## 9. Verdict

**POSITIVE-CALIBRATION-CLARIFICATION.**

The gravity arc has produced two distinct classes of result:

1. **Structural findings** (G1, G2, G3, G4, G5, G6-Diag-Full, G7): cutoff-independent properties of the framework — closed forms, theorems, structural mechanisms. These are framework predictions in the strong sense.

2. **Numerical gravity predictions** (G7 specifically): $G_{\rm eff} = 6\pi/\Lambda^2$, $\Lambda_{cc} = 6\Lambda^2$ for the Gaussian cutoff. These depend on the cutoff function choice.

G8 clarifies that this is the standard situation in CC's spectral action: the framework determines the structural form, the cutoff is external calibration data.

This is consistent with the broader structural-skeleton-scope reading: GeoVac determines the structural skeleton (Class A); calibration data (cutoff function, Yukawas, etc.) is external (Class 1).

## 10. Gravity arc consolidation status

After G1 through G8, the gravity arc has produced:

| Sprint | Output | Type |
|---|---|---|
| G1 | Two-term exactness; $\zeta_{\rm unit}(-k) = 0$ identity | Structural theorem |
| G2 | Factorization $K_{S^3 \times S^1} = K_{S^3} \cdot K_{S^1}$ | Structural theorem |
| G3 | Closed-form graviton heat traces | Structural theorem |
| G4 first-pass | Continuum $S_E = A/(4G)$ | Structural derivation |
| G4-1 | $S^2$ Dirac heat trace $a_1 = -4\pi/3$ | Structural calculation |
| G4-2 | $S_{BH} = A\Lambda^2/(12\pi)$ via replica | Structural mechanism |
| G5 | $s_{\rm min} = -\Lambda/(12\sqrt{6})$ | Closed-form extremum |
| G6-Diag | Quadratic form on (1,1) irrep | Structural setup |
| G6-Diag-Full | Gauge classification via $V = i[X, D_0]$ image | Structural method |
| G7 | $G_{\rm eff} = 6\pi/\Lambda^2$ (Gaussian) | Numerical prediction |
| **G8** | **General $G_{\rm eff}(f), \Lambda_{cc}(f)$ formulas** | **Calibration clarification** |

The gravity arc has reached natural sprint-scale closure. Multi-month G4-3 onwards (discrete-substrate construction) is the structural target for future commitments.

## 11. Files

- `debug/g8_cutoff_dependence.py` — symbolic derivation driver
- `debug/data/g8_cutoff_dependence.json` — structured results
- `debug/g8_cutoff_dependence_memo.md` — this memo
- `papers/group5_qed_gauge/paper_28_qed_s3.tex` — §4.16 added

## 12. Cross-references

- **G7** (`debug/g7_extremality_newton_memo.md`): Gaussian-cutoff $G_{\rm eff} = 6\pi/\Lambda^2$; G8 generalizes
- **External-input three-class partition** (`memory/external_input_three_class_partition.md`): cutoff function is Class 1 calibration
- **Structural-skeleton-scope pattern** (`memory/geovac_structural_skeleton_scope_pattern.md`): G8 sharpens the structural vs calibration boundary
- **Chamseddine-Connes 1997, 2010**: standard CC spectral action formalism with cutoff function
