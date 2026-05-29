# Sprint G5 — Spectral action on the decompactified $S^3 \times \mathbb{R}_\tau$

**Date:** 2026-05-28
**Path:** Gravity arc, G5. Sprint-scale extension of G2 making the $\beta \to \infty$ minimum explicit.
**Verdict:** **POSITIVE.** Heat-kernel decompactifies cleanly: $K_{S^1}(t,\beta)/\beta \to 1/(2\sqrt{\pi t})$. Spectral action density on $S^3_R \times \mathbb{R}_\tau$ is EXACTLY two-term: $s(R, \Lambda) = (R^3/4)\Lambda^4 - (R/8)\Lambda^2$. Extremum at $R_{\rm crit} \Lambda = 1/\sqrt{6}$ — same as G1, G2. Minimum action density $s_{\rm min} = -\Lambda/(12\sqrt{6})$ — clean closed form, symbolically verified.

## 1. Setup

G2 found the joint minimum of the spectral action on $S^3_R \times S^1_\beta$ is at $R_{\rm crit} \Lambda = 1/\sqrt{6}$ and $\beta \to \infty$, corresponding to the zero-temperature de Sitter vacuum. G5 takes the $\beta \to \infty$ limit explicitly, decompactifying the temporal direction to $\mathbb{R}_\tau$.

## 2. Heat-kernel decompactification

On $S^3 \times \mathbb{R}_\tau$, the heat kernel factorizes as
$$K_{\rm full}(t, R) = K_{S^3}(t, R) \cdot K_\mathbb{R}(t)$$
per unit length of $\mathbb{R}_\tau$, where
$$K_\mathbb{R}(t) = \int \frac{dk}{2\pi}\,e^{-k^2 t} = \frac{1}{2\sqrt{\pi t}}$$

This is the heat-kernel **density** on the real line. The compact-time analog from G2 was $K_{S^1}(t, \beta) = \beta/(2\sqrt{\pi t})$ + exp small; the per-length density $K_{S^1}/\beta \to 1/(2\sqrt{\pi t})$ matches $K_\mathbb{R}$ exactly.

Using $K_{S^3}(t, R) = (\sqrt{\pi}/2) R^3 t^{-3/2} - (\sqrt{\pi}/4) R t^{-1/2}$ + exp small (Paper 28 two-term exactness):
$$\boxed{\;K_{\rm full}(t, R) = \frac{R^3}{4}\,t^{-2} - \frac{R}{8}\,t^{-1} + O(\text{exp small})\;}$$

EXACTLY two-term — same structure as G2's per-$\beta$ density. Two-term exactness propagates to the decompactified case.

## 3. Spectral action density

For Gaussian cutoff $f(x) = e^{-x}$ at $\Lambda$, the spectral action density factorizes:
$$s(R, \Lambda) = K_{\rm full}(1/\Lambda^2, R) = \frac{R^3}{4}\Lambda^4 - \frac{R}{8}\Lambda^2$$

Verified symbolically by sympy. Same form as G2's per-$\beta$ density.

## 4. Extremum at $R_{\rm crit}\Lambda = 1/\sqrt{6}$

$\partial s/\partial R = 0$:
$$\frac{3R^2}{4}\Lambda^4 - \frac{1}{8}\Lambda^2 = 0 \implies R^2 = \frac{1}{6\Lambda^2}$$

$$\boxed{\;R_{\rm crit}\Lambda = \frac{1}{\sqrt{6}}\;}$$

— same as G1 (one-parameter $S^3_R$) and G2 (joint extremum in $R$ at fixed $\beta$). Symbolically verified: $R_{\rm crit} = \sqrt{6}/(6\Lambda) = 1/(\sqrt{6}\Lambda)$.

## 5. Minimum action density

At the extremum:
$$s_{\rm min} = s(R_{\rm crit}, \Lambda) = \frac{(1/\sqrt{6}\Lambda)^3}{4}\Lambda^4 - \frac{1/\sqrt{6}\Lambda}{8}\Lambda^2 = \frac{\Lambda}{24\sqrt{6}} - \frac{\Lambda}{8\sqrt{6}} = -\frac{\Lambda}{12\sqrt{6}}$$

$$\boxed{\;s_{\rm min} = -\frac{\Lambda}{12\sqrt{6}}\;}$$

Symbolically verified ($-\sqrt{6}\Lambda/72 = -\Lambda/(12\sqrt{6})$, sympy residual 0).

## 6. Physical reading: zero-temperature de Sitter vacuum

G5 makes G2's $\beta \to \infty$ joint minimum **explicit**:

- The decompactified time direction $\mathbb{R}_\tau$ is the natural setting for the "vacuum" state (no thermal compactification = zero temperature).
- The spectral action density on $S^3 \times \mathbb{R}_\tau$ has the same two-term structure as G2's per-$\beta$ density.
- The extremum $R_{\rm crit} \sim 1/\Lambda$ is the de Sitter vacuum radius (Planck-scale when $\Lambda$ is identified as Planck mass).
- The minimum action density $s_{\rm min} = -\Lambda/(12\sqrt{6})$ is the "vacuum action density" of the de Sitter vacuum.

This is the standard Connes-Chamseddine spectral-action-principle prediction of a Planck-scale de Sitter vacuum, with the well-known cosmological-constant-scale gap (observed $\Lambda_{cc} \ll \Lambda^3$ ~ Planck-scale).

## 7. Verdict and comparison

| Sprint | Geometry | Spectral action | Extremum |
|---|---|---|---|
| G1 | $S^3_R$ | $\phi(3/2) (\Lambda R)^3 - \tfrac{1}{4}\phi(1/2)(\Lambda R)$ | $u_{\rm crit} = R_{\rm crit}\Lambda = 1/\sqrt{6}$ |
| G2 | $S^3_R \times S^1_\beta$ | $\beta R^3 \Lambda^4/4 - \beta R \Lambda^2/8$ | $R_{\rm crit}\Lambda = 1/\sqrt{6}$, $\beta \to \infty$ |
| G5 | $S^3_R \times \mathbb{R}_\tau$ | $s = R^3 \Lambda^4/4 - R \Lambda^2/8$ | $R_{\rm crit}\Lambda = 1/\sqrt{6}$, $s_{\rm min} = -\Lambda/(12\sqrt{6})$ |

The structural finding: **the framework's preferred zero-temperature de Sitter vacuum** is described identically across the three geometries. The radius $R_{\rm crit} \Lambda = 1/\sqrt{6}$ is robust to the temporal direction (compact $S^1_\beta$ or non-compact $\mathbb{R}_\tau$, or absent in the 3D G1 case). The vacuum action density $-\Lambda/(12\sqrt{6})$ is the clean closed-form expression of this minimum.

**G2's $\beta \to \infty$ limit is now manifest in G5.** The "decompactification" is cleanly handled by the heat-kernel density per unit time.

## 8. Honest scope

**Reached:**
- Heat-kernel density on $\mathbb{R}_\tau$ ($K_\mathbb{R} = 1/(2\sqrt{\pi t})$) connects cleanly to $K_{S^1}/\beta$ ✓
- Two-term exactness propagates to the non-compact case ✓
- Spectral action density on $S^3 \times \mathbb{R}_\tau$ has same form as G2 per-$\beta$ ✓
- $R_{\rm crit} \Lambda = 1/\sqrt{6}$ explicit ✓
- $s_{\rm min} = -\Lambda/(12\sqrt{6})$ closed-form ✓
- Zero-temperature de Sitter vacuum interpretation explicit ✓

**Not reached:**
- Stefan-Boltzmann thermal corrections at finite $\beta$ (these vanish as $\beta \to \infty$, but their behavior at large but finite $\beta$ is a sub-question)
- Correction to vacuum from finite-$n_{\max}$ truncation of the GeoVac substrate (i.e., applying G5 to the discrete CH spectrum instead of the continuum heat-kernel asymptotic)
- Connection to actual cosmological observations (the $\Lambda \sim$ Planck-scale prediction has the standard CC cosmological-constant-problem gap)

## 9. Cross-references

- **Sprint G1** (`sprint_g1_path1_spectral_action_S3_radius.md`) — one-parameter $S^3_R$, same $u_{\rm crit}$
- **Sprint G2** (`sprint_g2_spectral_action_S3_x_S1.md`) — two-parameter $S^3 \times S^1_\beta$, $\beta \to \infty$ minimum identified
- **Sprint G3** (`sprint_g3_scalar_TT_S3.md`) — spinor-bundle specificity of two-term exactness (the gravity sector inherits standard CC continuum, NOT the clean form derived here)
- **Sprint G4 first-pass** (`debug/g4_cigar_BH_entropy_memo.md`) — Bekenstein-Hawking on cigar, structurally complementary to G5's de Sitter vacuum
- **Paper 28 §4.7, §4.8** (G1, G2 additions) — same two-term structure, now in non-compact case

## 10. Files produced

- `debug/g5_decompactified_S3_x_R.py` — symbolic derivation driver
- `debug/data/g5_decompactified.json` — structured results
- `debug/g5_decompactified_memo.md` — this memo
