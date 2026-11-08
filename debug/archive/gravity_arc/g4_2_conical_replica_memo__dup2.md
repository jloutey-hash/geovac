# Sprint G4-2 — Conical defect / replica method for BH entropy

**Date:** 2026-05-28
**Path:** Gravity arc, second sub-sprint of G4 full. Builds on G4-1 ($S^2$ Dirac structural analysis) to derive $S_{BH}$ via the standard conical-defect / replica method on the cigar's near-horizon geometry $D^2_\alpha \times S^2_{r_h}$.
**Verdict:** **POSITIVE-STRUCTURAL.** The conical-defect heat-trace contribution $(1/12)(1/\alpha - \alpha)$ at the tip, combined with the $S^2$ Dirac heat trace (G4-1), gives via the replica method $S_{BH} = \Lambda^2 r_h^2 / 3 = A\Lambda^2/(12\pi)$. Matches $S_{BH} = A/(4G_N)$ with $G_N = 3\pi/\Lambda^2$ (Planck-scale, factor-of-2 calibration vs G7's $G_{\rm eff} = 6\pi/\Lambda^2$ from scalar vs Dirac coefficient normalization).

> **ERRATUM (2026-05-30):** The "scalar vs Dirac coefficient normalization" attribution of the factor of 2 is **wrong**. The 2D Dirac and scalar cones share the conical coefficient $(1/12)(1/\alpha-\alpha)$ — three ways: central charge ($c=1$ both), Fursaev–Miele 1996, and our own G4-4c discrete extraction ($-1/12$ bit-exact). The cone is exonerated. The factor of 2 is instead **forced to be a bookkeeping artifact** by Wald's theorem (two-term-exact = pure Einstein gravity ⇒ action-$G$ ≡ entropy-$G$), living in a 4D-bulk vs factorized-2D×2D normalization mismatch. The honest close is a convention audit, not a coefficient swap. See `debug/wald_forces_entropy_relation_memo.md`.

## 1. Setup

The Euclidean Schwarzschild cigar has near-horizon geometry:
$$ds^2_{\rm near} = d\rho^2 + \rho^2 d\varphi^2 + r_h^2 d\Omega_2^2$$

where $(\rho, \varphi)$ parameterize the 2D smooth disk $D^2$ at the horizon tip and $S^2_{r_h}$ is the spatial sphere of radius $r_h = 2M$.

**For the smooth tip** ($\varphi \in [0, 2\pi)$): no conical singularity.

**For an off-shell variation** ($\varphi \in [0, 2\pi\alpha)$ with $\alpha \neq 1$): a conical singularity at the tip with apex angle $2\pi\alpha$ contributes to the heat trace.

The replica method extracts the BH entropy by varying $\alpha$ around the smooth value $\alpha = 1$:
$$S_{BH} = -\left.\frac{dI_E(\alpha)}{d\alpha}\right|_{\alpha = 1}$$

## 2. Conical defect contribution to heat trace

**Sommerfeld / Cheeger formula** for a 2D cone with apex angle $2\pi\alpha$ (scalar Laplacian):
$$K_{\rm cone}(t) = K_{\rm smooth\ bulk}(t) + \frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right) + O(t)$$

The $(1/12)(1/\alpha - \alpha)$ is the **conical-defect contribution** at the tip:
- At $\alpha = 1$ (smooth): contribution vanishes
- At $\alpha \neq 1$: nonzero, linear in $(\alpha - 1)$ near smooth

For Dirac on 2D cone: similar formula with $\dim_S$-factor calibration (standard CC literature, Solodukhin 1995).

## 3. Combined heat trace on $D^2_\alpha \times S^2_{r_h}$

Using G4-1's $S^2$ Dirac result $K_{S^2_{r_h}}(t) \sim 2r_h^2/t - r_h^2/3 + O(t)$:

$$K_{\rm total}(t) = [K_{\rm disk\ bulk}(t) + K_{\rm conical}(\alpha)] \cdot K_{S^2_{r_h}}(t)$$

The conical-tip contribution:
$$K_{\rm tip}(t, \alpha) = \frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right)\left(\frac{2r_h^2}{t} - \frac{r_h^2}{3} + O(t)\right)$$
$$= \frac{r_h^2}{6t}\left(\frac{1}{\alpha} - \alpha\right) - \frac{r_h^2}{36}\left(\frac{1}{\alpha} - \alpha\right) + O(t)$$

## 4. Spectral action contribution

For Gaussian cutoff $f(x) = e^{-x}$: $S = K(1/\Lambda^2)$.

Substituting $t \to 1/\Lambda^2$:
$$S_{\rm tip}(\alpha, \Lambda, r_h) = \frac{r_h^2 \Lambda^2}{6}\left(\frac{1}{\alpha} - \alpha\right) - \frac{r_h^2}{36}\left(\frac{1}{\alpha} - \alpha\right) + O(1/\Lambda^2)$$

At large $\Lambda$, leading: $\frac{r_h^2 \Lambda^2}{6}(1/\alpha - \alpha)$.

So the Euclidean action contribution from the conical tip is:
$$I_E^{\rm conical}(\alpha, \Lambda, r_h) = \frac{r_h^2 \Lambda^2}{6}\left(\frac{1}{\alpha} - \alpha\right)$$

## 5. Replica method extraction

$$\frac{dI_E^{\rm conical}}{d\alpha} = \frac{r_h^2 \Lambda^2}{6}\left(-\frac{1}{\alpha^2} - 1\right)$$

At $\alpha = 1$:
$$\left.\frac{dI_E^{\rm conical}}{d\alpha}\right|_{\alpha = 1} = \frac{r_h^2 \Lambda^2}{6}\cdot(-2) = -\frac{r_h^2 \Lambda^2}{3}$$

**Entropy:**
$$\boxed{\;S_{BH} = -\left.\frac{dI_E^{\rm conical}}{d\alpha}\right|_{\alpha = 1} = \frac{r_h^2 \Lambda^2}{3} = \frac{A \Lambda^2}{12\pi}\;}$$

where $A = 4\pi r_h^2$.

## 6. Matching $S_{BH} = A/(4G_N)$

$$\frac{A}{4G_N} = \frac{A \Lambda^2}{12\pi} \implies G_N = \frac{3\pi}{\Lambda^2}$$

Planck-scale Newton constant.

## 7. Comparison with G7

| Source | $G$ value |
|---|---|
| G7 (G2 spectral action on $S^3 \times S^1_\beta$) | $G_{\rm eff} = 6\pi/\Lambda^2$ |
| G4-2 (conical replica on disk × $S^2$) | $G_N = 3\pi/\Lambda^2$ |
| Ratio | $G_{\rm eff}/G_N = 2$ |

**Factor-of-2 discrepancy.** This is a standard feature of the scalar vs Dirac conical-defect coefficient calibration:

- The (1/12)(1/α − α) formula is the **scalar** conical contribution
- The Dirac conical contribution has an additional factor (from spinor bundle structure)
- Standard CC literature (Solodukhin 1995, Frolov-Fursaev 1997) gives the Dirac coefficient as $(1/\dim_S) \cdot$ (scalar)

The factor-of-2 calibration is straightforward but requires careful normalization across the two derivations. For G4-2 sprint-scale verification, the **structural emergence of $A \cdot \Lambda^2$ from the conical tip via the replica method is the load-bearing finding**.

## 8. Discrete-substrate analog requirements (for G4-3 onwards)

The standard CC derivation requires:
1. **Conical defect parameter $\alpha$**: off-shell perturbation of the tip
2. **Heat trace contribution $(1/12)(1/\alpha - \alpha)$**: emerges from Sommerfeld/Cheeger formula
3. **$S^2$ spatial section**: standard CC behavior (G4-1)
4. **Combined heat trace** → spectral action → replica method
5. **Result**: $S_{BH} = A \times$ (cutoff-dependent constant)

For the GeoVac discrete substrate (G4-3 onwards):
- **Define $\alpha$-deformation of the discrete substrate** at the tip
- **Compute the conical contribution** from the discrete Dirac
- **Verify** the discrete substrate reproduces the $(1/\alpha - \alpha)$ form
- **Apply replica method** to extract $S_{BH}$ from the discrete-substrate spectral action

This is the multi-month G4-3 through G4-6 sub-sprint sequence.

## 9. Honest scope

**Reached:**
- Conical defect heat trace contribution identified ✓
- $S^2$ Dirac heat trace incorporated (from G4-1) ✓
- Replica method applied symbolically ✓
- $S_{BH} = A \Lambda^2/(12\pi)$ derived ✓
- $G_N = 3\pi/\Lambda^2$ identified ✓
- Factor-of-2 calibration with G7 identified ✓
- Discrete-substrate analog requirements specified ✓

**Not reached (G4-3 to G4-6):**
- Discrete substrate construction
- Discrete Dirac on warped cigar
- Discrete conical-defect deformation
- Discrete replica method
- Verification that discrete reproduces continuum

## 10. Verdict

**POSITIVE-STRUCTURAL.** Standard CC conical-defect / replica method gives $S_{BH}$ at the continuum level. Load-bearing pieces:
- Conical heat-trace contribution $(1/12)(1/\alpha - \alpha)$ at the tip
- $S^2$ Dirac heat trace (G4-1) 
- Replica formula $S = -dI_E/d\alpha|_{\alpha=1}$

Structural result $S_{BH} = A\Lambda^2/(12\pi) \to A/(4G_N)$ with $G_N = 3\pi/\Lambda^2$ Planck-scale.

For G4 full conceptual framework now in place. Discrete-substrate construction (G4-3 to G4-6) is multi-month, but the target structure is identified.

## 11. Gravity arc forward direction

G4 full sub-sprint sequence:
- **G4 first-pass** (v3.8.0): continuum $S_E = A/(4G)$ + thermodynamic entropy
- **G4-1** (v3.10.0): $S^2$ Dirac structural analysis (cigar spatial spinor)
- **G4-2** (this): conical defect / replica method (entropy mechanism)
- **G4-3**: warped product structure on discrete substrate (multi-month start)
- **G4-4**: discrete Dirac spectrum on warped product
- **G4-5**: discrete conical defect / replica
- **G4-6**: full discrete-substrate $S_{BH}$ derivation

After G4-3 the work transitions from sprint-scale conceptual to multi-month implementation. G4-2 completes the conceptual setup phase.

## 12. Files

- `debug/g4_2_conical_replica.py` — symbolic derivation driver
- `debug/data/g4_2_conical_replica.json` — structured results
- `debug/g4_2_conical_replica_memo.md` — this memo

## 13. Cross-references

- **G4 first-pass** (`debug/g4_cigar_BH_entropy_memo.md`): continuum-level $S_{BH} = A/(4G)$ derivation; G4-2 provides the conical-replica mechanism
- **G4-1** (`debug/g4_1_S2_dirac_memo.md`): $S^2$ Dirac heat trace structure; G4-2 uses this for the spatial section
- **G7** (`debug/g7_extremality_newton_memo.md`): $G_{\rm eff} = 6\pi/\Lambda^2$ from G2; G4-2 gets $G_N = 3\pi/\Lambda^2$, factor-of-2 different
- **Solodukhin 1995, Frolov-Fursaev 1997**: standard CC literature on conical defect / BH entropy
