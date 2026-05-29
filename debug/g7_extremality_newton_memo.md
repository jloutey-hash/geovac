# Sprint G7 — Background extremality + Newton constant + cosmological constant

**Date:** 2026-05-28
**Path:** Gravity arc consolidation sprint. Combines two clean sprint-scale findings building on G2 and G6-Diag-Full.
**Verdict:** **POSITIVE.** (A) The Gaussian-weighted integral of the physical $(1,1)$ eigenvalue formula vanishes exactly — confirms background extremality. (B) Newton constant $G_{\rm eff} = 6\pi/\Lambda^2$ and cosmological constant $\Lambda_{cc} = 6\Lambda^2$ derived from G2 coefficients, identifying the framework's gravitational parameters in standard CC language with the well-known cosmological-constant-scale gap.

## 1. Part A: Background extremality consistency check

From G6-Diag-Full, the physical (1,1)-irrep graviton-candidate eigenvalues are
$$A_\lambda = a_\lambda \left(\frac{4\lambda^2}{\Lambda^4} - \frac{2}{\Lambda^2}\right), \qquad a_\lambda = e^{-\lambda^2/\Lambda^2}$$

The continuum-limit total contribution to the spectral action quadratic form, integrated over levels with the Gaussian weight, is (in dimensionless $u = \lambda/\Lambda$):
$$I = \int_0^\infty (4u^2 - 2)\,e^{-u^2}\,du$$

Component integrals (standard Gaussian moments):
- $\int_0^\infty u^2\,e^{-u^2}\,du = \sqrt{\pi}/4$
- $\int_0^\infty e^{-u^2}\,du = \sqrt{\pi}/2$

Therefore:
$$I = 4 \cdot \frac{\sqrt{\pi}}{4} - 2 \cdot \frac{\sqrt{\pi}}{2} = \sqrt{\pi} - \sqrt{\pi} = 0$$

**Exactly zero.** Symbolically verified by sympy.

**Interpretation.** This is the consistency condition that the spectral action is EXTREMIZED at the background. At an extremum, $\delta S = 0$ for arbitrary infinitesimal perturbations $V$. The integrated $A_\lambda$ formula is the kinetic-form trace along the (1,1) graviton-irrep direction; its vanishing reflects:
- **Spectral action is stationary at the background**, consistent with the variational principle defining the de Sitter vacuum (G1, G2, G5)
- **Sub-leading corrections** (finite-$n_{\max}$ deviations, non-Gaussian cutoffs) give the actual graviton kinetic structure
- **Diffeomorphism invariance** at the limit of the full CC framework is consistent with the framework's substrate-level structure

This is not a NEW dynamical claim — it's a STRUCTURAL CONSISTENCY CHECK that the framework's background and perturbation-theory analyses cohere.

## 2. Part B: Newton constant + cosmological constant from G2

G2's spectral action on $S^3_R \times S^1_\beta$ for Gaussian cutoff:
$$S(R, \beta, \Lambda) = \frac{\beta R^3}{4}\,\Lambda^4 - \frac{\beta R}{8}\,\Lambda^2$$

Standard Einstein-Hilbert + cosmological constant action:
$$S_{EH} = -\frac{1}{16\pi G_{\rm eff}}\int (R_{\rm scalar} - 2\Lambda_{cc})\sqrt{g}\,d^4x$$

### 2.1 Geometric integrals on $S^3_R \times S^1_\beta$

- Volume: $V = 2\pi^2 R^3 \beta$
- Scalar curvature: $R_{\rm scalar} = 6/R^2$ (only $S^3$ contributes; $S^1$ flat)
- Integrated scalar curvature: $\int R_{\rm scalar}\sqrt{g}\,d^4x = (6/R^2)(2\pi^2 R^3 \beta) = 12\pi^2 R\beta$

### 2.2 Matching coefficients

**$\Lambda^4$ term** (cosmological constant):
$$\frac{\Lambda_{cc}}{8\pi G_{\rm eff}}\cdot 2\pi^2 R^3 \beta = \frac{\beta R^3}{4}\Lambda^4$$
$$\frac{\Lambda_{cc}}{G_{\rm eff}} = \frac{\Lambda^4}{\pi}$$

**$\Lambda^2$ term** (Einstein-Hilbert):
$$-\frac{1}{16\pi G_{\rm eff}}\cdot 12\pi^2 R\beta = -\frac{\beta R}{8}\Lambda^2$$
$$\frac{1}{G_{\rm eff}} = \frac{\Lambda^2}{6\pi}$$

### 2.3 Solution

$$\boxed{\;G_{\rm eff} = \frac{6\pi}{\Lambda^2}, \qquad \Lambda_{cc} = 6\Lambda^2\;}$$

Symbolically verified in `debug/g7_extremality_newton.py` (sympy).

### 2.4 Planck-scale identification

Identifying $\Lambda$ as Planck mass $M_{\rm Planck}$:
- $G_{\rm eff} = 6\pi/M_{\rm Planck}^2 = 6\pi G_{\rm Newton}$ (Planck-scale Newton's constant)
- $\Lambda_{cc} = 6 M_{\rm Planck}^2 = 6/G_{\rm eff}$ (Planck-scale cosmological constant)

### 2.5 Cosmological constant scale problem

- Predicted: $\Lambda_{cc}/M_{\rm Planck}^2 = 6$
- Observed (cosmological dark energy): $\Lambda_{cc}/M_{\rm Planck}^2 \sim 10^{-120}$
- **Gap**: $\sim 120$ orders of magnitude

This is the **standard CC cosmological-constant-scale problem**, inherited from continuum Connes-Chamseddine, NOT specific to GeoVac. The framework's discrete substrate does not modify this leading-order prediction.

The standard interpretations:
- The spectral action's predicted $\Lambda_{cc}$ is the UV-cutoff-scale value; physical $\Lambda_{cc}$ would require renormalization or environmental selection (anthropic, multiverse)
- The cosmological-constant scale is calibration data (Class 1 per `memory/external_input_three_class_partition.md`), not derivable from the framework's first principles
- Matter contributions (Higgs sector in CC, scalar field on $S^3$ in GeoVac via G3's heat-kernel structure) shift $\Lambda_{cc}$ but don't close the 120-orders-of-magnitude gap

## 3. Connection to gravity arc

### 3.1 Sprint dependencies

- **G1** (`sprint_g1_path1_spectral_action_S3_radius.md`): two-term exactness on $S^3_R$, $\zeta_{\rm unit}(-k) = 0$ identity
- **G2** (`sprint_g2_spectral_action_S3_x_S1.md`): exact 4D CC form on $S^3 \times S^1_\beta$; G7's Newton constant + $\Lambda_{cc}$ derivation uses G2's coefficients directly
- **G3** (`sprint_g3_scalar_TT_S3.md`): spinor-bundle specificity; gravity sector inherits standard CC
- **G5** (`sprint_g5_decompactified_memo.md`): closed-form $s_{\rm min} = -\Lambda/(12\sqrt{6})$ on decompactified $S^3 \times \mathbb{R}_\tau$
- **G6-Diag-Full** (`sprint_g6_diag_full.md`): physical (1,1) eigenvalues $A_\lambda$; G7's extremality check uses this formula directly

### 3.2 Master Mellin engine classification

Per Paper 18 §III.7 (master Mellin engine), gravity content is M2 (Seeley-DeWitt heat-kernel coefficients $\sqrt{\pi} \cdot \mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$). G2's coefficients $1/(8\pi^2)$ and $-1/(96\pi^2)$ are M2 ring members.

GeoVac's gravity content is therefore:
- **M2 sector** (Seeley-DeWitt): captures Einstein-Hilbert + cosmological constant + (in continuum) higher curvature
- **M2 on Dirac substrate**: terminates at 2 terms (Paper 28 two-term exactness theorem) — clean structural property
- **M2 on scalar/tensor sectors**: full infinite series with $a_k = 2\pi^2/k!$ (G3) or similar — standard CC

This identifies gravity as living entirely within the master Mellin engine's M2 sub-mechanism.

## 4. Honest scope

**Reached:**
- Extremality consistency check verified ✓
- Newton constant $G_{\rm eff} = 6\pi/\Lambda^2$ derived ✓
- Cosmological constant $\Lambda_{cc} = 6\Lambda^2$ derived ✓
- Standard CC cosmological-constant-scale gap identified and acknowledged ✓
- Connection to master Mellin engine M2 sub-mechanism ✓

**Not reached:**
- Resolution of the cosmological-constant-scale problem (standard CC issue, not GeoVac-specific)
- Verification of $G_{\rm eff}$ against observed Newton's constant (requires identifying $\Lambda \sim M_{\rm Planck}$ from independent considerations)
- Matter contribution to $\Lambda_{cc}$ (would require G7-extension to compute scalar + Dirac + gauge contributions)
- Inner fluctuation analysis (Marcolli-vS gauge field structure)
- Discrete-substrate corrections to $G_{\rm eff}$ and $\Lambda_{cc}$ at finite $n_{\max}$

## 5. Verdict

Part A confirms structural consistency at the background; Part B identifies the framework's gravitational parameters in standard CC language. Together, these consolidate the gravity arc's connection to standard predictions:

- GeoVac reproduces CC gravity with $G_{\rm eff} \sim 1/M_{\rm Planck}^2$ and $\Lambda_{cc} \sim M_{\rm Planck}^2$
- Inherits the standard cosmological-constant-scale gap
- Does not provide a new mechanism to close this gap at the level analyzed

The gravity arc's net statement is now sharp:
- **CC gravity (Einstein-Hilbert + cosmological constant) emerges on the Dirac sector** with two-term exactness (G1, G2)
- **Gravitons exist at substrate level** but their physical interpretation requires multi-month full FP analysis (G6)
- **Bekenstein-Hawking entropy** standard derivation works at continuum (G4); discrete analog is multi-month
- **Cosmological constant scale gap** is standard CC issue, inherited

## 6. Files

- `debug/g7_extremality_newton.py` — symbolic derivation driver
- `debug/data/g7_extremality_newton.json` — structured results
- `debug/g7_extremality_newton_memo.md` — this memo

## 7. Cross-references

- **G2** — provides the spectral action coefficients used in Part B
- **G6-Diag-Full** — provides the (1,1) eigenvalue formula used in Part A
- **Paper 28 §4.7-4.12** — gravity arc structural sequence (G1 through G5)
- **Paper 18 §III.7** — master Mellin engine partition (M2 = gravity)
- **memory/external_input_three_class_partition** — cosmological constant scale as Class 1 calibration data
