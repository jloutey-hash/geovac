# Sprint G4 first-pass — Bekenstein-Hawking entropy from cigar spectral action

**Date:** 2026-05-28
**Type:** Conceptual / structural first-pass. Verifies the CC continuum derivation of $S_{BH} = A/(4G)$ from the Euclidean Schwarzschild cigar and identifies the load-bearing structural mechanism (horizon boundary heat-kernel coefficient). Identifies what GeoVac discrete substrate would need for the full derivation.
**Verdict:** **POSITIVE-CONCEPTUAL.** Standard CC derivation reproduces $S_{BH} = A/(4G)$ at the continuum level via symbolic computation. The load-bearing mechanism is the horizon boundary heat-kernel coefficient (giving GHY $\sim A/(4G)$ at the smooth conical tip). GeoVac discrete-substrate analog requires constructing a discrete spectral triple on the cigar — multi-month per scoping memo.

## 1. Setup

Euclidean Schwarzschild metric:
$$ds^2 = \left(1 - \frac{2M}{r}\right) d\tau^2 + \left(1 - \frac{2M}{r}\right)^{-1} dr^2 + r^2 d\Omega_2^2$$

Outside the horizon $r > 2M$. Regularity at $r = 2M$ requires $\tau$ periodic with period
$$\beta = 8\pi M = \frac{1}{T_H}$$

The cigar (Euclidean Schwarzschild) has:

| Quantity | Value |
|---|---|
| Horizon radius | $r_h = 2M$ |
| Hawking temperature | $T_H = 1/(8\pi M)$ |
| Inverse temperature | $\beta = 8\pi M$ |
| Horizon area | $A = 4\pi r_h^2 = 16\pi M^2$ |
| Bekenstein-Hawking entropy | $S_{BH} = A/(4G) = 4\pi M^2/G$ |

Sprint TD Track 4 (2026-05-08) already reproduced $T_H = 1/(8\pi M)$ from the Matsubara structure on the $\tau$-circle. G4 first-pass derives the entropy $S_{BH} = A/(4G)$ structurally.

## 2. Standard CC continuum derivation

**Euclidean on-shell action.** Gibbons-Hawking 1977: after subtraction of asymptotic flat-space reference, the on-shell action of Euclidean Schwarzschild is
$$I_E^{Schw} = \frac{A}{4G} = \frac{\beta M}{2}$$

(two equivalent forms; same numerical value via $\beta = 8\pi M$ and $A = 16\pi M^2$). The asymptotic boundary's $\beta M$ contribution cancels against the flat-space subtraction; the bulk Ricci-flat $R = 0$ vanishes; the GHY boundary at the horizon (smooth conical tip with $2\pi$ opening angle) contributes the $A/(4G)$.

**Thermodynamic entropy.** Using $Z = e^{-I_E}$ and $F = -T\ln Z = T I_E$, the entropy is
$$S_{thermo} = -\frac{\partial F}{\partial T} = -I_E + \beta\,\frac{\partial I_E}{\partial \beta}$$

Treating $A = A(\beta) = \beta^2/(4\pi)$ (since $M = \beta/(8\pi)$, $r_h = 2M = \beta/(4\pi)$, $A = 4\pi r_h^2 = \beta^2/(4\pi)$):
- $I_E = A/(4G) = \beta^2/(16\pi G)$
- $\partial I_E/\partial \beta = \beta/(8\pi G)$
- $\beta\,\partial I_E/\partial \beta = \beta^2/(8\pi G) = A/(2G)$
- $S_{thermo} = -A/(4G) + A/(2G) = A/(4G)$ ✓

Symbolically verified in `debug/g4_cigar_BH_entropy.py` (sympy, `S_thermo - A/(4G) = 0` exact).

In terms of $M$: $S_{thermo} = 4\pi M^2/G$ — the standard Bekenstein-Hawking result.

## 3. Heat-kernel coefficient structure

The CC spectral action expansion on a 4-manifold is
$$S_{spec} = c_4 \Lambda^4 \int \sqrt{g}\,d^4x + c_2 \Lambda^2 \int R\sqrt{g}\,d^4x + (\text{boundary terms}) + (\text{higher curvature})$$

On Euclidean Schwarzschild (Ricci-flat, $R = 0$):

| Term | Value |
|---|---|
| Bulk $\int R \sqrt{g}$ | 0 (Ricci-flat) |
| Bulk $\Lambda^4 V$ | divergent (non-compact); regulated by background subtraction |
| Asymptotic boundary | $\beta M$ (cancels against flat-space reference) |
| Horizon boundary (smooth conical tip) | **$-A/(4G)$** |

**Load-bearing structural piece**: the horizon boundary heat-kernel coefficient at the smooth conical tip ($2\pi$ opening angle). The standard GHY term integrated over the horizon $\Sigma = S^2_{r_h}$ gives:
$$\frac{1}{8\pi G}\oint_\Sigma K \sqrt{h}\,d^3x = \frac{A}{4G}$$

This is what generates the Bekenstein-Hawking entropy via the thermodynamic relation.

## 4. The GeoVac discrete substrate question

GeoVac's substrate (Camporesi-Higuchi Dirac on truncated $S^3$ at finite $n_{\max}$, sourced by the Fock projection) does NOT natively describe the cigar geometry:

- The cigar's spatial section is $\Sigma \cong \mathbb{R} \times S^2$ (radial + 2-sphere), NOT $S^3$
- The cigar's temporal direction is compact $S^1_\beta$ but the radial direction is non-compact $\mathbb{R}_r$ (asymptotic to $\infty$)
- The horizon at $r = 2M$ is a smooth conical tip in the $(\tau, r)$ plane, where the metric coefficient $(1 - 2M/r)$ vanishes

**To reproduce $S_{BH} = A/(4G)$ on a GeoVac-style discrete substrate would require:**

1. **Define a discrete spectral triple on the cigar** (multi-month)
   - Non-compact radial $\mathbb{R}_r$ direction needs regularization (asymptotic boundary cutoff)
   - Smooth conical tip at the horizon needs careful discretization (preserve the $2\pi$ regularity)
   - Spatial $S^2$ sphere instead of $S^3$
2. Compute the discrete Dirac spectrum on the cigar
3. Apply the spectral action with appropriate cutoff function
4. Verify the heat-kernel asymptotic produces the boundary coefficient equivalent to GHY at the horizon
5. Extract $S_{BH} = A/(4G)$ via the thermodynamic relation $S_{thermo} = -I_E + \beta\partial_\beta I_E$

**This is the multi-month G4 full sprint.** Closing requires non-trivial new infrastructure (discrete spectral triple on $\mathbb{R} \times S^2 \times S^1$ with horizon boundary).

## 5. Connection to the rest of the GeoVac framework

- **Sprint TD Track 4** (2026-05-08, "Cigar Hawking temperature"): reproduced $T_H = 1/(8\pi M)$ from Matsubara on the $\tau$-circle. POSITIVE-LIMITED verdict. G4 deferred at the time.
- **Sprint G2** (2026-05-28): exact 4D CC Einstein-Hilbert + cosmological constant on $S^3 \times S^1_\beta$. Provides the CC spectral action's bulk structure but NOT the boundary heat-kernel contributions needed for horizon BH entropy.
- **Sprint Unruh-pendant + Bisognano-Wichmann** (2026-05-10): four-witness Wick-rotation theorem (Hawking, Sewell, BW, Unruh) at operator-system level. Provides the modular structure but operates on $S^3$ background, not cigar.
- **Paper 49 §11 BCFM connection** (Bousso-Casini-Fisher-Maldacena 2020): bulk-dual reading of Connes cocycle flow. Provides AdS/CFT-side anchor for the BH entropy question but doesn't simplify the discrete-substrate construction.

## 6. Honest scope

**Reached:**
- Symbolic verification of $S_{thermo} = A/(4G)$ from Euclidean Schwarzschild action via thermodynamic relation ✓
- Identification of horizon boundary heat-kernel coefficient as the load-bearing mechanism ✓
- Structural framework for CC continuum derivation in place ✓

**Not reached:**
- Discrete spectral triple on the cigar (multi-month, several non-trivial components)
- Direct computation of the boundary heat-kernel coefficient on a discrete cigar
- Verification that the discretization preserves the GHY-equivalent structure
- Numerical extraction of $A/4$ from the discrete spectral action

## 7. Verdict

**POSITIVE-CONCEPTUAL.** The standard CC derivation of $S_{BH} = A/(4G)$ is well-understood at the continuum level and verified symbolically here. The load-bearing structural piece is the horizon boundary heat-kernel coefficient, which generates the GHY term that thermodynamic processing turns into the Bekenstein-Hawking entropy.

For GeoVac, the discrete-substrate analog requires constructing a discrete spectral triple on the cigar geometry. This is multi-month per the gravity-arc scoping. The first-pass establishes the structural target and the load-bearing mechanism.

## 8. Files

- `debug/g4_cigar_BH_entropy.py` — symbolic derivation driver
- `debug/data/g4_cigar_BH_entropy.json` — structured results (sympy expressions as strings)
- `debug/g4_cigar_BH_entropy_memo.md` — this memo

## 9. Cross-references

- **Sprint TD Track 4** (CLAUDE.md §2) — cigar Hawking temperature, POSITIVE-LIMITED, M2 deferred
- **Sprint G1, G2, G3** (v3.4.0-v3.6.0) — Path 1 structural-test phase; CC on $S^3$ and $S^3 \times S^1_\beta$
- **G6 scoping memo** — multi-month G6 is the analog for graviton dynamics; G4 full is the analog for BH entropy
- **Gibbons-Hawking 1977** — standard Euclidean derivation
- **CC literature** (Chamseddine-Connes 1997 et seq., Solodukhin 1995, Frolov-Fursaev 1997) — boundary heat-kernel coefficient calculations
