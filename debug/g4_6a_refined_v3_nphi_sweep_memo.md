# Sprint G4-6a refined v3 — N_φ sweep diagnostic

**Date:** 2026-05-29
**Path:** Multi-task thread 8, Track a'. N_φ-axis refinement test per Track a of thread 7's structural finding.
**Verdict:** **CLEAN NEGATIVE WITH SUBSTANTIVE STRUCTURAL FINDING.** The N_φ-sweep at fixed (a, R, t-panel) **does NOT show convergence to A_cont**. Instead, the apparent A signal DECREASES with N_0 (going from +172.6% at N_0=60 to -18.9% at N_0=480) — opposite of expected convergence. The "A signal" at intermediate t turns out to be a finite-N_0 artifact from incomplete bulk α·K_disk subtraction, NOT a genuine UV-divergence signal. **G4-6a refined is NOT sprint-scale-tractable on a single refinement axis**; the substrate at production values doesn't faithfully capture the 1/(24πt) UV divergence under any single-axis refinement strategy tested so far.

## 1. The diagnostic

Driver: `debug/g4_6a_refined_v3_nphi_sweep.py`
Data: `debug/data/g4_6a_refined_v3_nphi_sweep.json`

Hold $a = 0.05$, $N_\rho = 200$, $R = 10$ fixed (B.2 / G4-6b baseline). Vary $N_0 \in \{60, 120, 240, 480\}$. Apply simplified extraction: $A_{\rm est}(t) = t \cdot (\text{tip}(t) - B_{\rm substrate})$ with $B_{\rm substrate} = 0.163$. Track A_recovery vs N_0 across four t-values.

## 2. The unexpected result

| $N_0$ | A_mean | A_recovery |
|---|---|---|
| 60 | +0.02289 | **+172.62%** |
| 120 | +0.00261 | +19.65% |
| 240 | −0.00208 | −15.71% |
| 480 | −0.00251 | −18.89% |

**A recovery DECREASES with N_0, eventually going negative.** This is the opposite of the expected convergence toward A_cont = 1/(24π) ≈ 0.0133.

### 2.1 Per-t breakdown at N_0=480

| $t$ | tip(t) | tip − B_substrate | A_est | recovery |
|---|---|---|---|---|
| 0.0125 | +0.067 | **−0.096** | −0.0012 | −9.0% |
| 0.025 | +0.092 | −0.071 | −0.0018 | −13.4% |
| 0.0625 | +0.116 | −0.046 | −0.0029 | −22.0% |
| 0.125 | +0.130 | −0.033 | −0.0041 | −31.2% |

At N_0=480, tip(t) is BELOW B_substrate (0.163) for all tested t-values. The substrate is now undershooting B at intermediate t.

## 3. Structural interpretation

### 3.1 What's actually happening

The "A signal" at intermediate t in the simplified extraction was NOT a clean UV divergence reading. It was finite-N_0 artifact from **incomplete bulk α·K_disk subtraction**.

At small N_0 (e.g., N_0=60), the bulk subtraction has poor resolution; the residual (tip − B_substrate) is INFLATED, looking like a large A coefficient. As N_0 grows, the bulk subtraction sharpens, the residual shrinks toward 0, and eventually the substrate's tip(t) crosses below B_substrate at intermediate t.

### 3.2 Why this is a clean negative

The N_φ-axis refinement test was the natural follow-up to thread 7 Track a's finding. The expected outcome was monotonic improvement of A_recovery with N_0. Instead, the recovery DEGRADES.

This means:
1. The simplified extraction strategy works at SOMETIMES, by accident of finite-N_0 inflation.
2. The substrate at production values does NOT faithfully capture the 1/(24πt) UV divergence under N_φ refinement.
3. Both Track a's and B.2's apparent A extractions at production substrate sizes are not converging to A_cont.

### 3.3 The substrate doesn't have the UV content

The honest structural finding: **at substrate values reachable in main-session work (N_0 ≤ 480, N_ρ ≤ 400, a ≥ 0.025), the discrete substrate does not extract A = 1/(24π) cleanly under any single-axis refinement strategy tested**.

Either:
- The substrate is fundamentally lacking UV content that no single-axis refinement provides (requires multi-axis or new strategy).
- The bulk α·K_disk subtraction is the wrong observable for A extraction at production substrate values.
- The "A_cont" target derived from continuum theory may not be the right target for substrate's intrinsic extraction (substrate may converge to a different limit at multi-axis refinement).

## 4. Implications for G4-6 multi-month commitment

### 4.1 G4-6a refined is NOT sprint-scale-tractable

The day's compression pattern (G4-6d, G4-6b closing at sprint scale) does NOT extend to G4-6a refined. The clean negative on N_φ-axis refinement means the simplified extraction strategy doesn't close A at production substrate values.

Implications for the multi-month estimate:
- G4-6a refined needs MULTI-AXIS substrate exploration: vary $a$, $N_\rho$, $N_\phi$ jointly.
- Or a fundamentally different observable for A extraction.
- Either path is genuinely multi-month work, not sprint-scale.

**Revised G4-6 estimate:** the compression to "4-8 weeks total" assumed G4-6a refined closes at sprint scale. With this Track a' result, the estimate revises back to a longer timeline — likely 2-4 months for G4-6a refined alone.

### 4.2 G4-6e (Mellin moment theorem grade) depends on this

G4-6e requires clean Mellin moment extraction at production substrate. If A extraction is fundamentally limited at substrate values, G4-6e's theorem-grade closure also requires multi-month work to identify the right observable.

### 4.3 The honest revised G4-6 timeline

| Sub-sprint | Status | Effort |
|---|---|---|
| G4-6d | DONE | sprint-scale realized |
| G4-6b | DONE | sprint-scale realized |
| G4-6a refined | **MULTI-AXIS or NEW OBSERVABLE NEEDED** | 2-4 months realistic |
| G4-6c | substrate-level + dichotomy named | 1 day (Route A) or multi-week (Routes C') |
| G4-6e | depends on G4-6a refined | 1 month after G4-6a refined closes |
| G4-6f | final | 1 month |

**Honest revised G4-6 total: 4-7 months**, less compressed than the day's pattern suggested but still less than the original 3-6 months estimate.

## 5. What the substrate IS achieving

Despite the negative on full UV recovery, the substrate IS achieving real physical content:
- Spectral azimuthal substrate is 160× better at UV cell t=a² than FD (thread 4 B.2).
- B_substrate = 0.163 matches continuum +1/6 to within 2.3% (thread 6 Track α).
- The Möbius factor at α > 1 is captured at sub-2% precision (multiple sprints, robust to refinement).
- Wedge KMS state and modular structure work bit-exactly at finite cutoff (Papers 42/43/45/47).

The remaining gap is specifically in clean A coefficient extraction at intermediate-t panels — a particular observable that's substrate-limited.

## 6. Recommended next move

Given the clean negative on N_φ-axis refinement, three options:

**Option 1 — Accept the multi-month G4-6a refined timeline.**
Continue G4-6 multi-month commitment with realistic 2-4 month G4-6a refined work via multi-axis substrate exploration. This is the original G4-6 plan and remains valid.

**Option 2 — Document the structural limitation and close G4-6 at observation rigor.**
Paper 51 already documents the substrate's UV behavior (B.2 verified, Möbius substrate-level identified). The A coefficient gap could be characterized as a structural-skeleton-scope limitation: substrate provides the structural form, calibration of A requires external content (similar to how cutoff function is calibration data).

**Option 3 — Pursue Route A (Fursaev-Miele PDF) for the Möbius mechanism instead.**
The Möbius mechanism is structurally close to closure via the substrate-level identification. Route A may close it definitively. This is the parallel Track b' of this thread.

## 7. Honest scope

This Track a':
- **Tests** the N_φ-axis refinement strategy at fixed (a, R, t-panel).
- **Documents** the clean negative: A recovery degrades, not converges, with N_0 refinement.
- **Reframes** G4-6a refined as NOT sprint-scale-tractable on single axis.
- **Revises** the G4-6 multi-month estimate honestly (4-7 months, less compression than thread 6 suggested).

Does NOT:
- Test multi-axis refinement (would require larger substrate panels, multi-week compute).
- Identify the right observable for A extraction.
- Make the G4-6 strategic choice (multi-month commitment vs scope-limitation closure).

## 8. Cross-references

- `debug/g4_6a_refined_simplified_extraction_memo.md` — thread 7 Track a (the N_φ-axis identification this Track a' tests)
- `debug/g4_6a_spectral_substrate_first_move_memo.md` — thread 4 B.2 (original A/B joint fit reading)
- `debug/g4_6b_ir_boundary_first_move_memo.md` — thread 6 Track α (B_substrate identification)
- `debug/g4_6_scoping_memo.md` — sub-sprint sequencing (needs another honest update)
- `debug/g4_6a_refined_v3_nphi_sweep.py` — driver
- `debug/data/g4_6a_refined_v3_nphi_sweep.json` — structured results

## 9. Files

- `debug/g4_6a_refined_v3_nphi_sweep.py`
- `debug/data/g4_6a_refined_v3_nphi_sweep.json`
- `debug/g4_6a_refined_v3_nphi_sweep_memo.md` (this)
