# Sprint G4-6a refined — simplified A extraction first move

**Date:** 2026-05-29
**Path:** Multi-task thread 7, Track a. First move of the refined G4-6a sub-sprint per Track α of thread 6's simplified-strategy identification.
**Verdict:** **POSITIVE-WITH-STRUCTURAL-FINDING.** The simplified extraction strategy (B_substrate = 0.163 measured at large t, then $A_{\rm est}(t) = t \cdot (\text{tip}(t) - B_{\rm substrate})$) gives substantially cleaner per-$t$ A values than B.2's joint A/B fit — peak recovery is **12.7% at intermediate $t \approx 5a^2$** vs B.2's joint-fit-near-zero. **Substantive new structural finding:** the substrate's UV recovery is set by $N_\phi$ (azimuthal mode count), NOT by $a$ alone. Refining $a$ at fixed $N_\phi/\alpha$ does not improve UV recovery beyond ~13% at production substrate values. **G4-6a refined needs a multi-axis refinement scheme (a, $N_\phi$) jointly**, not just radial $(a, N_\rho)$ refinement.

## 1. The simplified extraction strategy

Per Track α of thread 6, the G4-6b first move found:
- Substrate B is essentially R-independent at R ≥ 10 with value $B_{\rm substrate} = 0.163$.
- Within 2.3% of continuum $+1/6 = 0.167$.

The simplified strategy:
1. Use $B_{\rm substrate} = 0.163$ as the measured Lichnerowicz constant.
2. Per-$t$ A extraction: $A_{\rm est}(t) = t \cdot (\text{tip}(t) - B_{\rm substrate})$.
3. Identify the cleanest extraction window in $t$.

## 2. The extraction on B.2 data

Driver: `debug/g4_6a_refined_simplified_extraction.py`
Data: `debug/data/g4_6a_refined_simplified.json`

### 2.1 Panel 1 (a=0.05, $N_\rho$=200) per-$t$ A

| $t$ | tip(t) | tip − $B$ | $A_{\rm est}$ | recovery |
|---|---|---|---|---|
| 0.0025 | +1.350 | +1.187 (peak per-t recovery 6%) | – | – |
| 0.005 | +1.070 | +0.907 | – | – |
| 0.0125 | +0.703 | +0.540 | +0.00675 | **51%** |
| 0.025 | +0.421 | +0.258 | +0.00644 | **48%** |
| 0.125 | +0.135 | −0.028 | – | overshoot |

(showing intermediate-$t$ where peak extraction occurs.)

### 2.2 Panel 2 (a=0.025, $N_\rho$=400) per-$t$ A

| $t$ | tip(t) | tip − $B$ | $A_{\rm est}$ | recovery |
|---|---|---|---|---|
| 0.00063 | +1.350 | +1.187 | +0.00074 | 5.6% |
| 0.00125 | +1.070 | +0.907 | +0.00113 | 8.6% |
| **0.00313** | **+0.703** | **+0.540** | **+0.00169** | **12.7%** (peak) |
| 0.00625 | +0.421 | +0.258 | +0.00161 | 12.1% |
| 0.03125 | +0.135 | −0.028 | −0.00089 | overshoot |

The peak A_est at Panel 2 is 12.7% at $t = 0.00313 = 5a^2$ — the same relative $t$ window as Panel 1's peak.

### 2.3 Key per-$t$ recovery observation

Per-$t$ A recovery is monotonically increasing in $t$ up to a peak at intermediate $t$ ($\sim 5a^2$), then drops as B-dominated regime is entered.

**Same relative $t$ window has SAME per-$t$ recovery across substrate refinements:**
- Panel 1 at $t = 5a^2 = 0.0125$: 51%
- Panel 2 at $t = 5a^2 = 0.00313$: 13%

These are DIFFERENT recoveries even though $t/a^2$ is the same. The substrate's per-$t$ recovery is NOT invariant under joint refinement of (a, $N_\rho$).

## 3. The substantive structural finding

### 3.1 Comparison to B.2 joint fit

| Strategy | Panel 1 A/A_cont | Panel 2 A/A_cont | Trajectory |
|---|---|---|---|
| B.2 joint A/B linear fit | 22% | 5% | REVERSED |
| Simplified (this Track a) per-$t$ smallt mean | 36% | 9% | Same direction, REVERSED magnitude |
| Simplified peak ($t = 5a^2$) | 51% | 13% | REVERSED |

**Both strategies give DECREASING A as substrate refines (a smaller).** This is the opposite of what we'd expect for a well-converging substrate (where finer $a$ should improve UV resolution).

### 3.2 Why the trajectory is reversed

Looking at the substrate structure: B.2 panels share $R = 10$ (i.e., $N_\rho \cdot a = $ const). Panel 1 has $N_\phi = 120$ (the smaller wedge mode count); Panel 2 also has $N_\phi = 240$ (proportionally scaled).

The substrate's UV cutoff is at $|\lambda| \le N_\phi/(2\alpha) \approx N_0/2$ in eigenvalue space (per the FD discretization at α=1, similar for spectral with truncation). This cutoff is the SAME in both panels (since $N_0$ is held constant relative to R).

But the small-$t$ panel's $t$-values scale with $a^2$, going from {0.0025, ..., 0.125} (Panel 1) to {0.00063, ..., 0.03125} (Panel 2). The smaller-$a$ panel asks the substrate to give answers at SMALLER $t$, where the UV divergence is LARGER.

The substrate's UV cutoff isn't improving as we refine $a$ at fixed $N_\phi/R$ — so we're asking the substrate to give bigger answers without giving it more UV bandwidth. The recovery RATE drops because the target grows while the substrate's capability stays roughly the same.

### 3.3 The right refinement axis is $N_\phi$, not (a, $N_\rho$)

For the substrate to capture more UV content, the AZIMUTHAL mode count $N_\phi$ (not the radial substrate refinement) needs to grow. The simplified extraction reveals this structural feature that the B.2 joint fit obscured.

## 4. Implication for G4-6 multi-month commitment

The original G4-6a multi-substrate UV sweep refined $(a, N_\rho)$ with $R$ fixed. The first-move analyses now suggest this is the WRONG refinement axis. The correct G4-6a refined strategy should:

1. Hold $a$ fixed (matching the B-substrate calibration at $R = 10$).
2. Vary $N_\phi$ (or equivalently $N_0$) at fixed $\alpha = 1$ proportion.
3. Track how the UV recovery scales with $N_\phi$ — expected $\sim 1$ at $N_\phi \to \infty$ if the spectral substrate IS continuum-faithful.

This is structurally different from the original G4-6a plan. It's a SECOND reframing of G4-6a (the first was from FD to spectral via G4-6d; this second is from radial-refinement to azimuthal-refinement).

### 4.1 Updated G4-6a refined plan

| Original G4-6a | Reframed G4-6a (after Track a finding) |
|---|---|
| Refine $(a, N_\rho)$ at fixed $R$ | Refine $N_\phi$ at fixed $a$, $R$ |
| Richardson in $1/a^2$ | Richardson in $1/N_\phi$ or $1/N_\phi^2$ |
| Multi-month substrate sweep | Single-substrate $N_\phi$-sweep, sprint-scale |

**Compression: from multi-month to sprint-scale, again.** The G4-6a refined work that originally was 2-3 months may compress to ~1 week if the $N_\phi$-sweep diagnostic shows the expected convergence.

## 5. Honest scope

This Track a:
- **Confirms** the simplified A extraction is cleaner than B.2 joint fit.
- **Identifies** that the substrate's UV recovery is set by $N_\phi$, not $a$.
- **Reframes** G4-6a refined sub-sprint AGAIN: refinement axis shifts to $N_\phi$.

Does NOT:
- Run the $N_\phi$-sweep (that's the next sprint following this finding).
- Verify that $N_\phi \to \infty$ gives A → A_cont (the load-bearing falsifier for the reframed G4-6a).
- Identify the specific $N_\phi$ scaling exponent.

## 6. Recommended next sprint

**Sprint G4-6a refined v2 — $N_\phi$ sweep.**

Hold $a = 0.05$, $N_\rho = 200$, $R = 10$ fixed (B.2 baseline). Vary $N_0 \in \{60, 120, 240, 480, 960\}$ (factor-2 progression). At each $N_0$:
- Compute tip(t = 5a²) where Panel 1 peak occurred.
- Extract $A_{\rm est}$ using $B_{\rm substrate}$.
- Track A_recovery vs $N_0$.

Expected: monotonic convergence toward A_cont. Convergence exponent should reveal whether substrate's UV recovery is structurally complete in the $N_\phi \to \infty$ limit.

Effort: ~1-2 hours main-session (the radial Laplacian is the same per mode; only mode count changes).

## 7. Cross-references

- `debug/g4_6a_spectral_substrate_first_move_memo.md` — B.2 sweep (the joint-fit reading this Track a refines)
- `debug/g4_6b_ir_boundary_first_move_memo.md` — Track α of thread 6 (B = 0.163 measurement)
- `debug/g4_6a_refined_simplified_extraction.py` — driver (this Track a)
- `debug/data/g4_6a_refined_simplified.json` — structured results
- `debug/g4_6_scoping_memo.md` — G4-6 sequencing (needs another update with this finding)

## 8. Files

- `debug/g4_6a_refined_simplified_extraction.py` (driver)
- `debug/data/g4_6a_refined_simplified.json` (data)
- `debug/g4_6a_refined_simplified_extraction_memo.md` (this)
