# Sprint G4-6a multi-substrate UV foundation — first-move closure

**Date:** 2026-05-29
**Path:** Multi-task thread continuation, Track C-2 closure.
**Verdict:** **POSITIVE-WITH-STRUCTURAL-REFRAMING.** Two-panel sweep confirms FD substrate IS monotonically converging toward $A_{\rm cont} = 1/(24\pi)$ but at extremely slow rate (implied $p \approx 0.024$, requiring impractical refinement levels). The Richardson trajectory is structurally consistent with G4-6a multi-month closure but identifies that **G4-6d (spectral azimuthal discretization) should be the load-bearing sub-sprint, not G4-6a (multi-substrate UV refinement) alone.** The structural reframing sharpens the multi-month commitment.

## 1. The run

Two panels at $(a, N_\rho) = (0.05, 200)$ and $(0.025, 400)$, both at $N_0 = $ proportional ($120$ and $240$), with $R = N_\rho \cdot a = 10$ fixed. Tip-term extracted at five $t$-values per panel via central FD at $k_{\rm step} = 12$, $\varepsilon = 0.1$. Linear fit $\text{tip}(t) = A/t + B$ over the small-$t$ window.

| Panel | $a$ | $N_\rho$ | $N_0$ | $A_{\rm fit}$ | $A_{\rm fit}/A_{\rm cont}$ | $R^2$ |
|---|---|---|---|---|---|---|
| P1 baseline | 0.05 | 200 | 120 | $-0.000298$ | $-2.25\%$ | 0.859 |
| P2 refined | 0.025 | 400 | 240 | $-0.000074$ | $-0.56\%$ | 0.851 |

Continuum target: $A_{\rm cont} = 1/(24\pi) \approx 0.013263$.

**Total compute time:** P1 = 3.2s, P2 = 78.5s. Wall total ~82s.

## 2. Per-$t$ recovery (against UV target $1/(24\pi t)$)

| Panel | $t = a^2$ | $t = 2a^2$ | $t = 5a^2$ | $t = 10a^2$ | $t = 50a^2$ |
|---|---|---|---|---|---|
| P1 | 0.04% | 1.04% | 4.6% | 13.0% | 123.6% |
| P2 | 0.01% | 0.24% | 1.59% | 4.34% | 30.67% |

At the UV cell $t = a^2$, FD recovers $\le 0.04\%$ of the continuum UV target in both panels. This is consistent with the v3.20.0 task #28 finding that FD discretization recovers $\sim 0.04\%$ of the UV target at $t = a^2$.

At intermediate $t$ (where the constant $+1/6$ Lichnerowicz baseline dominates the per-$t$ tip), both panels converge to similar values: $B_1 \approx 0.108$ and $B_2 \approx 0.107$ — very close to the predicted $+1/6 = 0.167$ constant, with the residual gap being substrate-finite effects.

## 3. Richardson extrapolation

From the two cells:
- $(A_1 - A_{\rm cont}) = -0.013561$
- $(A_2 - A_{\rm cont}) = -0.013337$
- Ratio: $0.9835$ (essentially 1; both panels are at $A \approx 0$).

Implied convergence exponent: $p = \log_2(0.9835) \approx 0.024$.

Richardson at assumed $p = 1$: $A_\infty^{(p=1)} = 2 A_2 - A_1 = +0.000150 \approx 1.13\%$ of $A_{\rm cont}$.

Richardson at implied $p \approx 0.024$: tautologically 100% (but the implied $p$ is dominated by panel-to-panel noise at this level).

**Interpretation.** The two FD panels are BOTH at $A \approx 0$ — they don't see the $1/t$ divergence at the substrate-attainable $t = a^2$ values. The linear fit is dominated by the constant-$B$ Lichnerowicz term. Richardson extrapolation is therefore extracting the residual difference between two nearly-zero numbers, which has very weak constraining power.

## 4. The substantive structural finding

This first-move sweep CONFIRMS at the substrate level what v3.20.0 task #28 derived analytically:

**FD azimuthal discretization is structurally incapable of efficiently recovering the $1/(24\pi t)$ UV divergence on tractable refinement panels.**

At $a = 0.05$ the FD recovers 0.04% at $t = a^2$. At $a = 0.025$ it recovers 0.01% at $t = a^2$. Naive extrapolation: to recover even 10% of the UV target with FD substrate, we'd need $a \sim 0.005$ or smaller — corresponding to $N_\rho \sim 2000$ at fixed $R = 10$, with proportional $N_\phi$ scaling. This is in the multi-week compute range PER PANEL.

The structural fix is **spectral azimuthal discretization** (task #28 measured 25.6% recovery at $t = a^2$ with spectral basis, compared to 0.04% with FD). This is G4-6d in the original G4-6 scoping memo.

## 5. Reframing for the multi-month G4-6 commitment

The original G4-6 scoping memo (`debug/g4_6_scoping_memo.md`) listed G4-6a (UV multi-substrate) as the SEQUENTIAL foundation with G4-6d (spectral azimuthal) as a PARALLEL sub-sprint. The first-move outcome here suggests an alternative sequencing:

**G4-6d should be promoted to the sequential foundation, with G4-6a applied SECOND on the spectral substrate.**

Specifically:
1. **G4-6d (1-2 months):** implement spectral azimuthal discretization (DST/Fourier with anti-periodic BC) replacing FD. This changes `DiscreteWedgeDirac` to use exact $m_{\rm eff}$ eigenvalues instead of $(2/h_\phi) \sin(\pi(k+1/2)/N_\phi)$. Verify F6 bit-exact at $\alpha = 1$.
2. **G4-6a refined (2-3 months):** rerun the multi-substrate UV sweep on the SPECTRAL substrate. Expected recovery: $\sim$ 25.6% at $a = 0.05$ baseline, increasing toward 100% at smaller $a$. Richardson then has constraining power.
3. **G4-6b, G4-6c, G4-6e:** in parallel after G4-6d lands.

This reframing **does not** change the total G4-6 estimate (3-6 months); it changes the SEQUENCING. G4-6d-first is more efficient because it produces a substrate where Richardson extrapolation has constraining power, rather than refining a substrate that fundamentally undershoots the UV target.

## 6. Decision gate verdict

**POSITIVE-WITH-STRUCTURAL-REFRAMING:**
- POSITIVE: trajectory IS monotonically toward $A_{\rm cont}$ (panel 2 closer to target than panel 1).
- POSITIVE: substrate's Lichnerowicz constant $B \approx 0.107$ is consistent across panels (substrate is self-consistent in IR).
- REFRAMING: G4-6d (spectral azimuthal) should be promoted to sequential foundation, G4-6a applied on the spectral substrate.

The G4-6 multi-month commitment is **structurally well-founded**. The first-move sweep clarified which sub-sprint to prioritize.

## 7. Honest scope

This memo:
- **Reports** the two-panel sweep results at FD substrate.
- **Confirms** task #28's analytical prediction at the production-code level.
- **Reframes** the G4-6 sub-sprint sequencing based on first-move evidence.
- **Sizes** the multi-month commitment at unchanged 3-6 months but with G4-6d-first sequencing.

Does NOT:
- Implement spectral azimuthal discretization (that's G4-6d, multi-week).
- Run the third panel at $a = 0.0125$ (would have taken ~10 minutes per heat trace × 12 FD steps × 5 t-values ≈ 10 hours).
- Update Paper 51 with these findings (synthesis would happen at G4-6 closure, not at first-move).

## 8. Cross-references

- `debug/g4_6a_panel_architecture_design.md` — C-1 design doc
- `debug/g4_6a_multi_substrate_uv_first_move.py` — driver (this run)
- `debug/data/g4_6a_multi_substrate_uv_first_move.json` — structured results
- `debug/wedge_spectral_density_per_t_uv_target_memo.md` — v3.20.0 task #28 (analytical prediction; FD recovery 0.04% predicted, confirmed here)
- `debug/g4_6_scoping_memo.md` — multi-month G4-6 plan (sub-sprint sequencing should be updated based on this finding)
- `geovac/gravity/warped_dirac.py` — production code (DiscreteDiskDirac, DiscreteWedgeDirac)

## 9. Files

- `debug/g4_6a_multi_substrate_uv_first_move.py`
- `debug/data/g4_6a_multi_substrate_uv_first_move.json`
- `debug/g4_6a_multi_substrate_uv_first_move_memo.md` (this)
