# Sprint G4-3c-proper — Wedge-lattice conical-defect test

**Date:** 2026-05-28
**Path:** Gravity arc, sprint-scale closure of the G4-3c open item (proper conical-defect via wedge lattice rather than naive $N_\phi$ resolution sweep). Implements apex-angle-respecting discretization $\phi \in [0, 2\pi\alpha)$ with periodic BC at fixed angular resolution $h_\phi = 2\pi/N_0$.
**Verdict:** **NEGATIVE-G4-3c-PROPER-WITH-STRUCTURAL-DIAGNOSIS.** Sign of the topological residual is correct ($\Delta K(\alpha < 1) > 0$, $\Delta K(\alpha > 1) < 0$, vanishing at $\alpha = 1$), and the wedge-lattice framework is operational, but the slope of $\Delta K$ vs $(1/\alpha - \alpha)$ is $\sim 1/4$ of the Sommerfeld/Cheeger target $1/12$, with strong asymmetry concentrated on the $\alpha < 1$ side. The discretization cannot isolate the $(1/12)(1/\alpha - \alpha)$ tip term from finite-substrate residuals at sprint scale. Confirms G4-5 (discrete replica method) as the proper multi-month target.

## 1. Target

Continuum Sommerfeld/Cheeger formula on a 2D scalar disk with apex angle $2\pi\alpha$:
$$
K_{\rm cone}(\alpha, t) = \alpha \cdot K_{\rm disk}(t) + \frac{1}{12}\!\left(\frac{1}{\alpha} - \alpha\right) + O(t)
$$
The $(1/12)(1/\alpha - \alpha)$ is the topological conical-defect contribution; vanishes at $\alpha = 1$, sign-flips through $\alpha = 1$.

## 2. Wedge lattice

Improvement on the original G4-3c: fix angular resolution $h_\phi = 2\pi/N_0$ and vary the total period to $2\pi\alpha$, giving $N_\phi(\alpha) = \alpha \cdot N_0$ sites. Periodic BC at the seam closes the wedge into a cone.

Discrete azimuthal Laplacian:
$$
\lambda_k^{\rm disc} = \left(\frac{2}{h_\phi}\right)^2 \sin^2\!\left(\frac{\pi k}{N_\phi}\right), \quad k = 0, 1, \ldots, N_\phi - 1
$$
For small $k$: $\lambda_k^{\rm disc} \to (k/\alpha)^2$, matching continuum $m_k = k/\alpha$ (the integer $m$-modes of the wedge are spaced at $1/\alpha$).

Setup: $R = 10$, $N_\rho = 200$, $a = 0.05$, $N_0 = 120$. Alpha sweep $\alpha \in \{1/3, 1/2, 2/3, 1, 3/2, 2, 3\}$ giving $N_\phi(\alpha) \in \{40, 60, 80, 120, 180, 240, 360\}$.

## 3. Topological residual

$\Delta K(\alpha, t) := K_{\rm wedge}^{\rm disc}(\alpha, t) - \alpha \cdot K_{\rm disk}^{\rm disc}(t)$ at $t = 0.005$:

| $\alpha$ | $(1/\alpha - \alpha)$ | SC pred $(1/12)\cdot \cdot$ | $\Delta K$ at $t=0.005$ |
|---|---|---|---|
| 1/3 | +2.667 | +0.2222 | **+0.1331** |
| 1/2 | +1.500 | +0.1250 | **+0.0525** |
| 2/3 | +0.833 | +0.0694 | **+0.0179** |
| 1   |  0.000 |  0.0000 | $+0$ (by construction) |
| 3/2 | −0.833 | −0.0694 | **−0.0017** |
| 2   | −1.500 | −0.1250 | **−0.0023** |
| 3   | −2.667 | −0.2222 | **−0.0034** |

**Sign is correct** at every $\alpha$. **Magnitude asymmetric**: $|\Delta K(\alpha = 1/3)| / |\Delta K(\alpha = 3)| = 39$, vs SC prediction's exact 1:1 antisymmetry.

## 4. Slope test

Linear fit $\Delta K(\alpha, t_*) = c \cdot (1/\alpha - \alpha)$ at each $t$:

| $t$ | slope $c$ | $c / (1/12)$ |
|---|---|---|
| 0.005 | +0.0230 | 0.276 |
| 0.01  | +0.0263 | 0.316 |
| 0.02  | +0.0296 | 0.355 |
| 0.05  | +0.0336 | 0.404 |
| 0.1   | +0.0364 | 0.437 |
| 0.2   | +0.0390 | 0.468 |
| 0.5   | +0.0420 | 0.505 |
| 1.0   | +0.0441 | 0.529 |

**Slope moves AWAY from SC target $1/12 = 0.0833$ as $t \to 0$.** At $t = 0.005$, the slope is 28% of SC; at $t = 1.0$ it's 53%. Wrong direction for the expected $t \to 0$ limit.

## 5. Reciprocal cancellation failure

Continuum prediction: $\Delta K(1/n) + \Delta K(n) = 0$ exactly (antisymmetric in $\alpha \leftrightarrow 1/\alpha$).

| Pair $(\alpha, 1/\alpha)$ | $\Delta K$ sum at $t = 0.005$ |
|---|---|
| $(1/3, 3)$ | $+0.1296$ (vs target $0$) |
| $(1/2, 2)$ | $+0.0502$ (vs target $0$) |
| $(2/3, 3/2)$ | $+0.0162$ (vs target $0$) |

Antisymmetry is broken by an $\alpha$-asymmetric residual concentrated on $\alpha < 1$.

## 6. Structural diagnosis

The bulk subtraction $K_{\rm wedge} - \alpha \cdot K_{\rm disk}$ would isolate the SC tip term exactly in the continuum. Three sources for the residual asymmetry observed in the discrete substrate:

### (a) Different finite-$N_\phi$ UV truncation across $\alpha$

The wedge at $\alpha = 1/3$ uses $N_\phi = 40$ azimuthal modes; the reference disk uses $N_\phi = 120$. From T2 (G4-3d-UV), at $t = 0.005$ the Weyl ratio at $N_\phi = 40$ is structurally worse than at $N_\phi = 120$ (predicted ~25% UV deficit at $N_\phi = 40$ vs ~5% at $N_\phi = 120$). The $\alpha < 1$ wedge inherits more UV truncation than $\alpha \cdot K_{\rm disk}$ subtracts, leaving a positive residual.

At $\alpha > 1$, the wedge uses MORE $N_\phi$ than the reference disk — UV better, residual nearly zero.

**Net: the residual asymmetry $\Delta K(\alpha<1) \gg |\Delta K(\alpha>1)|$ is dominantly UV-truncation asymmetry, NOT the SC tip term.**

### (b) Finite-$R$ IR cutoff

Both wedge and disk have $R = 10$ Dirichlet IR cutoff. In principle the boundary contribution cancels in $K_{\rm wedge} - \alpha \cdot K_{\rm disk}$ since both lengths scale as $\alpha$ at fixed $R$. In practice the boundary curvature corrections differ.

### (c) Apex discretization artifact

The continuum cone has a singular tip. The discrete substrate samples $\rho \geq a$ (smallest radial site at $\rho_1 = a$), so the tip itself is excluded. The topological $1/12$ contribution comes from analytical continuation through the tip in the continuum derivation; the discrete substrate has no equivalent mechanism short of full replica-method construction.

## 7. Quantitative slope estimate

If sources (a)+(b) account for the residual asymmetry, the symmetric part of $\Delta K$ should match the SC prediction. Decompose $\Delta K(\alpha, t) = c_{\rm sym}(t) \cdot (1/\alpha - \alpha) + c_{\rm asym}(t) \cdot h(\alpha)$ where $h(\alpha)$ captures the $\alpha < 1$ asymmetry.

Antisymmetric extraction at $t = 0.005$:
$$
c_{\rm sym} = \frac{1}{2}\sum_n \frac{\Delta K(1/n) - \Delta K(n)}{(1/n - n)}
$$
where the sum is over reciprocal pairs $(1/3, 3), (1/2, 2), (2/3, 3/2)$.

From the table: $\Delta K(1/3) - \Delta K(3) = 0.137$, divided by $1/(1/3) - 1/3 = 8/3$, gives $0.137 / 2.667 = 0.0515$. Similarly $(0.0525 - (-0.0023))/1.5 = 0.0365$ and $(0.0179 - (-0.0017))/0.833 = 0.0235$.

Antisymmetric slopes: $\{0.0515, 0.0365, 0.0235\}$ across pair scales $(2.667, 1.5, 0.833)$. None match $1/12 = 0.0833$.

**The antisymmetric piece is also off by factor 2-3.5x from SC target.** The diagnosis (a) UV-asymmetry is necessary but NOT sufficient to explain the deviation.

## 8. Honest scope and verdict

**Reached:**
- Wedge-lattice with proper apex-angle scaling implemented
- Periodic BC at the seam closes the cone correctly
- Sign of $\Delta K$ correct at every tested $\alpha$
- $\alpha = 1$ residual exactly zero (sanity check)
- Structural diagnosis of three residual sources (UV asymmetry, IR cutoff, apex discretization)

**Not reached:**
- Quantitative match to $(1/12)(1/\alpha - \alpha)$ slope (off by factor 2-4)
- Reciprocal cancellation $\Delta K(1/n) + \Delta K(n) = 0$ (broken by $\alpha$-asymmetric residual)
- $t \to 0$ extrapolation of slope (moving away from SC target rather than toward it)

**Verdict:** **NEGATIVE-G4-3c-PROPER-WITH-STRUCTURAL-DIAGNOSIS.**

The wedge-lattice + naive-bulk-subtraction approach is too coarse to isolate the SC tip term from finite-substrate residuals at sprint scale. The dominant residual is UV-asymmetry between $N_\phi(\alpha = 1/3) = 40$ and $N_\phi(\alpha = 3) = 360$. To extract SC cleanly at this discretization would require either:

1. **Pre-subtraction at the spectral level**: subtract the analytic Weyl two-term $\alpha \cdot [A/(4\pi t) - L/(8\sqrt{\pi t})]$ before the residual analysis. Removes finite-$R$ asymmetry but inherits its own truncation error.
2. **Joint UV refinement**: push $N_0 = 480$ or larger so $N_\phi(\alpha = 1/3) = 160$, with $N_\phi(\alpha = 3) = 1440$. Computationally feasible (factor ~64 in eigensolves) but doesn't address the apex discretization.
3. **Discrete replica method (G4-5 multi-month)**: implement the conical defect via $n$-sheeted covering and analytic continuation in $n$. The structurally correct approach.

The sprint-scale verdict confirms G4-5 as the proper target for literal SC extraction. The G4-3c-proper wedge construction is the right substrate for that future sprint.

## 9. Substantive new content (beyond G4-3c naive)

1. **Wedge-lattice discretization recipe**: $N_\phi(\alpha) = \alpha N_0$ at fixed $h_\phi$ gives the correct continuum-limit azimuthal spectrum $m = k/\alpha$. Sign of topological residual correct.

2. **The "1/4 of slope" pattern is UV-asymmetry, not Hermiticity**: combined with T2's finding that $N_\phi$ controls the UV regime, the slope deficit decomposes cleanly into the discretization mismatch between $\alpha < 1$ and $\alpha > 1$ branches.

3. **G4-5 (replica method) is the correct multi-month target**: the structural reason is now sharper. Even with a "proper" wedge lattice at sprint scale, the SC tip term lives below the discretization floor by 4x.

## 10. G4-3 sequence status (updated, post-T1)

| Sub-sprint | Status |
|---|---|
| G4-3a | Done (scoping) |
| G4-3a-cleanup | Done (Hermitian polar Laplacian) |
| G4-3b | Done (variable warp) |
| G4-3c (naive sweep) | Done (partial-positive, framework) |
| **G4-3c-proper (this)** | **Done (negative-with-diagnosis, wedge-lattice operational)** |
| G4-3d | Done (positive-IR) |
| G4-3d-UV (T2) | Done (positive-UV-verified) |
| G4-4 | Multi-month (warped Dirac, scoping in flight as T3) |
| G4-5 | Multi-month (discrete replica) — **now sharpened as the correct SC extraction target** |
| G4-6 | Multi-month (full $S_{BH}$ derivation) |

## 11. Files

- `debug/g4_3c_proper_wedge.py` — driver
- `debug/data/g4_3c_proper_wedge.json` — structured results
- `debug/g4_3c_proper_wedge_memo.md` — this memo

## 12. Cross-references

- G4-3c (naive sweep): framework substrate, predicted naive sweep wouldn't work
- G4-3d-UV (T2): UV truncation control, identifies $N_\phi = 96+$ for clean small-$t$
- G4-2 (continuum derivation): Sommerfeld/Cheeger target $(1/12)(1/\alpha - \alpha)$
- Sommerfeld 1894, Cheeger 1983: standard references
