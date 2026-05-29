# Sprint G4-6a — Multi-substrate UV panel architecture design

**Date:** 2026-05-29
**Path:** Multi-task thread continuation, Track C-1. Design doc preceding C-2 implementation.

## 1. Goal

Per the continuation prompt, the G4-6a sub-sprint extracts the leading $1/(24\pi)$ coefficient of the tip term across multiple substrate refinements, demonstrating that Richardson extrapolation toward the continuum target is feasible.

## 2. Production code state (verified)

`geovac/gravity/warped_dirac.py` (1049 lines) exists with:
- `DiscreteDiskDirac(N_rho, a, N_phi)` — α=1 disk reference
- `DiscreteWedgeDirac(N_rho, a, N_phi, alpha)` — wedge at general α
- Both provide `heat_trace(t)` and `squared_eigenvalues()` methods
- Anti-periodic spinor BC, rank-2 spinor bundle doubling

Existing G4-4f driver (`debug/g4_4f_replica_dK_dalpha.py`) computes the tip term via central FD:
$$\text{tip}(t) = \frac{K_{\text{wedge}}(\alpha_+, t) - K_{\text{wedge}}(\alpha_-, t)}{\alpha_+ - \alpha_-} - K_{\text{disk}}(t)$$
with $\alpha_\pm = 1 \pm k/N_0$ for integer $k$ keeping $N_\phi = \alpha N_0$ integer.

## 3. C-2 first-move scope

Panel cells: two tractable cells (not three) to keep main-session compute under control.

| Panel | $a$ | $N_\rho$ | $N_0$ | Compute cost |
|---|---|---|---|---|
| 1 | 0.05 | 200 | 120 | ~few seconds |
| 2 | 0.025 | 400 | 240 | ~30-60 seconds |

The largest cell ($a = 0.0125$, $N_\rho = 800$, $N_0 = 480$) requires $\sim$ several minutes per heat trace × 12 = 1-2 hours; defer to follow-on if Panel 1 + Panel 2 trend is informative.

## 4. Algorithm per panel

For each $(a, N_\rho, N_0)$:

1. Compute `disk = DiscreteDiskDirac(N_rho, a, N_phi=N_0)`.
2. Compute tip term at $t = a^2$ via central FD with $k = 12$ (eps = 0.1, matching G4-4f stable point).
3. Compute tip term at $t \in \{a^2, 2a^2, 5a^2, 10a^2, 50a^2\}$ (small-$t$ regime where $1/(24\pi t)$ dominates).
4. Linear fit `tip(t) = A/t + B` over the small-$t$ panel.
5. Extract $A$ and compare to $A_{\rm cont} = 1/(24\pi) \approx 0.01326$.

## 5. Richardson extrapolation across panels

Given $A_1$ at $a_1 = 0.05$ and $A_2$ at $a_2 = 0.025 = a_1/2$:

Expected convergence: $A(a) = A_{\rm cont} + C \cdot a^p + O(a^{p+1})$ for some $p > 0$.

Richardson estimate at $a \to 0$:
$$A_\infty^{\rm Richardson} \approx \frac{2^p A_2 - A_1}{2^p - 1}$$
with $p$ either assumed (e.g., $p = 1$ for linear convergence) or fitted from two cells (one parameter from two values).

Two-panel cells are at the minimum needed to estimate $p$; a third panel would tighten the rate estimate. Verdict on G4-6a closure trajectory:
- If $A_1 \approx 0.25 \cdot A_{\rm cont}$ and $A_2 > A_1$ with clean Richardson trajectory toward $A_{\rm cont}$ → POSITIVE foundation for G4-6a multi-month.
- If $A_2 \le A_1$ or trajectory is non-monotonic → identify whether substrate has structural obstruction; reassess G4-6a.

## 6. Output

- `debug/g4_6a_multi_substrate_uv_first_move.py` — driver
- `debug/data/g4_6a_multi_substrate_uv_first_move.json` — structured panel results + Richardson extrapolation
- `debug/g4_6a_multi_substrate_uv_first_move_memo.md` — closure verdict memo

## 7. Honest scope

This is FIRST-MOVE, not full G4-6a closure. The continuation prompt sized G4-6a at 2-3 months. The first-move panel sweep demonstrates whether the trajectory IS converging cleanly; if yes, the multi-month commitment is well-founded. If no, we surface a structural obstruction early.

## 8. Cross-references

- `debug/g4_6_scoping_memo.md` — multi-month G4-6 architectural plan
- `debug/wedge_spectral_density_per_t_uv_target_memo.md` — task #28, per-t UV target $1/(24\pi t)$ derivation
- `debug/g4_4f_replica_dK_dalpha.py` — existing tip-term driver (template)
- `geovac/gravity/warped_dirac.py` — production code (DiscreteDiskDirac, DiscreteWedgeDirac)
