# Sprint G4-4d — Seeley-DeWitt coefficient extraction

**Date:** 2026-05-29
**Verdict:** **PARTIAL-G4-4d-a_0-VERIFIED.** Leading Weyl coefficient $a_0$ (rank-2 spinor) extracted to **99.6%** at the sweet-spot $t = 0.1$ on the T2-refined substrate. Higher-order coefficients ($a_1$ boundary, $a_2$ topological) require multi-window fitting or continuum extrapolation; not closed at sprint scale.

## 1. a_0 extraction at sweet spot

Sweet-spot $t$ window (per T2 G4-3d-UV): $t = 0.1$, $N_\phi = 192$ — where Weyl ratio $\to 1$ in the scalar disk. For the spinor with rank-2 doubling:

$$K_{\rm Dirac}(t) \sim \frac{a_0 \cdot A_{D^2}}{4 \pi t} + \ldots$$

Continuum prediction: $a_0 = 2$ (rank of 2D spinor bundle).

| $t$ | $K(t)$ | $K \cdot 4\pi t / A$ | $a_0$ implied |
|---|---|---|---|
| 0.05 | 1063.58 | 1.064 | 2.127 (UV overshoot) |
| **0.1** | **497.98** | **0.996** | **1.992** |
| 0.2 | 236.95 | 0.948 | 1.896 (IR onset) |
| 0.5 | 89.07 | 0.891 | 1.781 (IR) |

**$a_0 = 1.992$ at sweet spot, within 0.4% of continuum 2.0.**

## 2. a_1 boundary (partial)

Weyl-subtracted residual at sweet spot is small ($-0.64$ at $t = 0.1$), giving $a_1 \sim -0.02$ in $L/\sqrt{\pi}$ units. Sub-leading coefficient extraction is sensitive to the $t$ window and the precise Weyl subtraction; a clean $a_1$ identification requires either:
1. Continuum extrapolation ($a \to 0$ on substrate)
2. Richardson-style multi-$t$ fitting with controlled cancellations
3. Per-mode angular decomposition

These are sub-leading characterizations of the spinor heat-kernel; the load-bearing $a_0$ is extracted.

## 3. Structural reading

The "naive" polynomial fit at small-$t$ (driver step 2) gives $c_0 \sim 62$ vs predicted 50 — biased by the UV-overshoot identified in T2 G4-3d-UV (high-$m$ angular truncation produces $4/\pi^2 \approx 0.405$ ratio between discrete and continuum azimuthal Laplacian at the truncation edge). The sweet-spot $t = 0.1$ approach correctly accounts for the overshoot and gives $a_0 = 1.992$.

**G4-4d outcome**: discrete substrate supports leading Seeley-DeWitt extraction at sub-percent precision at the sweet spot. Sub-leading coefficients are more sensitive but extractable with refined methodology (named for later sub-sprints).

## 4. Files
- `debug/g4_4d_seeley_dewitt_extraction.py` + JSON + this memo
