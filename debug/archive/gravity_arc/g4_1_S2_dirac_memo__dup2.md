# Sprint G4-1 — Dirac on $S^2$ and the cigar's spatial spinor structure

**Date:** 2026-05-28
**Path:** Gravity arc, first sub-sprint of G4 full. Tests whether the framework's clean two-term exactness on $S^3$ Dirac propagates to the $S^2$ Dirac that appears in the cigar's spatial section.
**Verdict:** **POSITIVE-STRUCTURAL-FINDING.** $S^2$ Dirac has standard CC infinite SD series, NOT the two-term exact form of $S^3$ Dirac. The framework's clean two-term exactness is specific to the half-integer-shifted Dirac spectrum on $S^3$ (the GeoVac substrate). For G4 full: the cigar's spatial $S^2$ component inherits standard CC continuum behavior with higher-curvature corrections at all orders.

## 1. Why this matters for G4 full

G4 full requires constructing a discrete spectral triple on the cigar (Euclidean Schwarzschild) and computing $S_{BH} = A/4$ from the discrete substrate. The cigar's spatial section is the warped product $\mathbb{R}_+ \times S^2_r$ where $r$ is the radial coordinate.

The natural first question: does the framework's structural cleanness on $S^3$ Dirac (two-term exactness, Paper 28 theorem) propagate to the $S^2$ Dirac that appears in this spatial section? If YES, the discrete construction inherits the clean structure. If NO, it follows standard CC continuum behavior.

This sub-sprint answers this question (NO, see below), which sharpens the scope of G4 full.

## 2. Setup

**Dirac on $S^d$ (Camporesi-Higuchi 1996):**
$$|\lambda_n^{S^d}| = n + d/2$$
with multiplicity depending on dimension and chirality.

**For $S^2$ ($d = 2$):**
- Eigenvalues: $|\lambda_n^{S^2}| = n + 1$ (**INTEGER** shift, $n = 0, 1, 2, \ldots$)
- Multiplicity per chirality (Weyl): $g_n^{\rm Weyl} = 2(n+1)$
- Total (full Dirac): $g_n^{\rm Dirac} = 4(n+1)$

**For $S^3$ ($d = 3$, GeoVac substrate):**
- Eigenvalues: $|\lambda_n^{S^3}| = n + 3/2$ (**HALF-INTEGER** shift)
- Multiplicity: $g_n^{\rm Dirac} = 2(n+1)(n+2)$

## 3. Heat trace on $S^2$ Dirac

$$K_{S^2}^{\rm Dirac}(t) = \sum_n 4(n+1)\,e^{-(n+1)^2 t}$$

Setting $u = n + 1$:
$$K(t) = 4 \sum_{u \geq 1} u\,e^{-u^2 t}$$

### 3.1 Leading asymptotic

At small $t$, the sum is dominated by $u \sim 1/\sqrt{t}$, approximating by an integral:
$$K(t) \approx 4 \int_0^\infty u\,e^{-u^2 t}\,du = 4 \cdot \frac{1}{2t} = \frac{2}{t}$$

Numerical verification (mpmath, 50 dps, $n_{\max} = 2000$):

| $t$ | $K_{\rm direct}$ | $K_{\rm leading} = 2/t$ | rel diff |
|---|---|---|---|
| $0.1$ | $19.663$ | $20.000$ | $1.7 \times 10^{-2}$ |
| $0.01$ | $199.667$ | $200.000$ | $1.7 \times 10^{-3}$ |
| $0.001$ | $1999.667$ | $2000.000$ | $1.7 \times 10^{-4}$ |
| $0.0001$ | $19999.67$ | $20000.00$ | $1.7 \times 10^{-5}$ |

Matches expected $K(t) \sim 2/t = \dim_S \cdot \text{Vol}(S^2)/(4\pi t)$ with $\dim_S = 2$ and $\text{Vol} = 4\pi$.

### 3.2 First subleading (the $a_1$ coefficient)

Standard SD form: $K(t) = (4\pi t)^{-1}[a_0 + a_1 t + a_2 t^2 + \ldots]$

For unit $S^2$ with Dirac operator (Vassilevich 2003, Branson-Gilkey):
- $a_0 = \dim_S \cdot V = 2 \cdot 4\pi = 8\pi$ (giving leading $2/t$ ✓)
- $a_1 = \dim_S \cdot (R_{\rm scalar}/6 - E_{\rm Lich}) \cdot V$

On unit $S^2$: $R_{\rm scalar} = 2$, $E_{\rm Lich} = R_{\rm scalar}/4 = 1/2$. So
$$a_1 = 2 \cdot (1/3 - 1/2) \cdot 4\pi = 2 \cdot (-1/6) \cdot 4\pi = -\frac{4\pi}{3}$$

**$a_1 \neq 0$**, in stark contrast to $S^3$ Dirac where all higher SD coefficients vanish.

**Numerical extraction** of $a_1$ from the heat trace:
$$K(t) - \frac{2}{t} \to \frac{a_1}{4\pi} = -\frac{1}{3} \text{ as } t \to 0$$

| $t$ | $K(t) - 2/t$ | Expected $-1/3$ |
|---|---|---|
| $0.1$ | $-0.3367$ | $-0.3333\ldots$ |
| $0.05$ | $-0.3350$ | $-0.3333\ldots$ |
| $0.01$ | $-0.3337$ | $-0.3333\ldots$ |
| $0.001$ | $-0.33337$ | $-0.3333\ldots$ |

Confirmed: $a_1 = -4\pi/3$ is the first subleading SD coefficient. Higher $a_k$ ($k \geq 2$) are also nonzero (standard CC infinite series).

## 4. Comparison: $S^3$ Dirac vs $S^2$ Dirac

| | $S^3$ Dirac | $S^2$ Dirac |
|---|---|---|
| Spectrum shift | **half-integer** $n + 3/2$ | **integer** $n + 1$ |
| Heat trace | $\frac{\sqrt{\pi}}{2}t^{-3/2} - \frac{\sqrt{\pi}}{4}t^{-1/2}$ | $\frac{2}{t} - \frac{1}{3} + O(t)$ |
| Two-term exact? | YES (Paper 28 theorem) | NO |
| $a_2$, $a_3$, ... | $= 0$ | $\neq 0$ |
| Jacobi inversion | $\theta_2$ (half-integer) | $\theta_3$ (integer) + derivative |

**Structural distinction:** two-term exactness on $S^3$ Dirac comes from the half-integer shift, which makes the Jacobi $\theta_2$ inversion give exactly two power-law terms. $S^2$ Dirac uses integer shift with $\theta_3$ inversion + exponential prefactor, giving the standard CC infinite SD series.

## 5. Consistency with Sprint G3

G3 established that the SCALAR Laplacian on $S^3$ has full SD series with closed-form $a_k = 2\pi^2/k!$. The pattern is:
- **Half-integer Dirac on $S^3$**: TWO-TERM exact (Paper 28)
- **Integer scalar Laplacian on $S^3$**: full SD series ($a_k = 2\pi^2/k!$, G3)
- **Integer Dirac on $S^2$**: full SD series ($a_0 = 8\pi$, $a_1 = -4\pi/3$, etc., G4-1)

**Two-term exactness is SPECIFIC to half-integer-shifted Dirac.** Other dimensions / operators / sectors inherit standard CC continuum behavior.

## 6. Implications for G4 full

The cigar's spatial section $\mathbb{R}_+ \times S^2_r$ has:
- $S^2$ component: standard CC infinite SD series (G4-1)
- $\mathbb{R}_+$ radial direction: needs regularization (warped product, $r$-dependence)
- Horizon boundary at $r = 2M$: where the BH entropy contribution lives

**The discrete-substrate BH entropy derivation cannot exploit a clean two-term form on the spatial $S^2$.** It must follow the standard CC heat-kernel asymptotic structure with higher-curvature corrections at all orders.

**For G4 full this means:**
1. Define a discrete spectral triple on $\mathbb{R}_+ \times S^2 \times S^1_\beta$ with horizon boundary
2. The discrete $S^2$ component will have integer-shifted spectrum (Dirac on the framework's $S^2$ analog)
3. The heat-kernel asymptotic on the discrete substrate will follow standard CC structure (not two-term)
4. The horizon boundary contribution must give $A/(4G)$ via the boundary $a_2$ coefficient — this is the load-bearing piece

The framework's clean structure on the Dirac sector ($S^3$ specifically) does NOT extend to G4 full. The discrete-substrate construction inherits the standard CC complexity.

## 7. Honest scope of this sub-sprint

**Reached:**
- $S^2$ Dirac spectrum identified ($\lambda_n = n + 1$, $g_n = 4(n+1)$) ✓
- Heat trace leading asymptotic $K(t) \to 2/t$ verified numerically ✓
- $a_1 = -4\pi/3 \neq 0$ identified and verified ✓
- Structural distinction $S^3$ vs $S^2$ documented ✓
- Implications for G4 full identified ✓

**Not reached (subsequent G4 sub-sprints):**
- G4-2: define the warped product $\mathbb{R}_+ \times S^2_r$ structure on a discrete substrate
- G4-3: compute the discrete Dirac spectrum on the warped product
- G4-4: identify the boundary heat-kernel coefficient at the horizon
- G4-5: verify $A/(4G)$ emerges from the boundary contribution
- G4-6: thermodynamic entropy extraction

These are the remaining sub-sprints of G4 full.

## 8. Verdict

**POSITIVE-STRUCTURAL-FINDING.** The framework's clean two-term exactness on $S^3$ Dirac (Paper 28 theorem) does NOT propagate to the $S^2$ Dirac that appears in the cigar's spatial section. Two-term exactness is **specific to the half-integer-shifted $S^3$ Dirac spectrum (the GeoVac substrate)**.

For G4 full, this means the discrete-substrate BH entropy derivation must follow standard CC continuum behavior on the spatial $S^2$, with higher-curvature corrections at all orders. The framework's structural cleanness on the gravity arc is limited to G1/G2 (Dirac sector on $S^3 \times S^1_\beta$) and does not extend to G4 full's cigar geometry.

## 9. Files

- `debug/g4_1_S2_dirac.py` — driver (symbolic derivation + numerical verification)
- `debug/data/g4_1_S2_dirac.json` — structured results
- `debug/g4_1_S2_dirac_memo.md` — this memo

## 10. Cross-references

- **G3** (`sprint_g3_scalar_TT_S3.md`) — established spinor-bundle specificity of two-term exactness; G4-1 confirms this extends to other dimensions ($S^2$)
- **G4 first-pass** (`debug/g4_cigar_BH_entropy_memo.md`) — continuum-level $S_{BH} = A/(4G)$ derivation; G4-1 is the first sub-sprint of the discrete-substrate version
- **Paper 28 §4.9** (G3 addition) — three closed forms; G4-1 adds the $S^2$ Dirac case to this family
- **Camporesi-Higuchi 1996** — general $S^d$ Dirac spectrum formulas

## 11. Next sub-sprint

**G4-2: warped product structure for $\mathbb{R}_+ \times S^2_r$ on the discrete substrate.**

The natural next step: define the warped product structure for the cigar's spatial section on a discrete substrate. This involves:
- Choosing a radial discretization (Laguerre? finite difference? algebraic basis?)
- Specifying how the $S^2$ radius varies with the discrete radial coordinate
- Identifying the boundary at the horizon
- Computing the heat-trace on the warped-product substrate

Sprint-scale (probably 1-2 sessions).
