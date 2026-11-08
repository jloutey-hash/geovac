# Sprint GD-5 — helicity-2 modes are positive too: GD-4 concern resolved at necessary-condition level

**Date:** 2026-05-29. **Type:** ① follow-on. **Verdict:** POSITIVE — the helicity-2 (|Δj|=2, TT-graviton-class) within-sector modes carry the same positive A_λ as the (1,1), because A_λ is irrep-blind within a block. So G6's necessary condition extends to the actual TT-graviton irrep sector; the GD-4 "G6 tested the wrong sector" concern is resolved at the necessary-condition level. The TT-vs-trace *distinction* (sufficient condition) remains the multi-month FP task.

## 1. The question (from GD-4)

GD-4 found G6 tested the (1,1) = helicity-0 (scalar/trace-class) sector, while the TT graviton is helicity-2 (|Δj|=2). Do the helicity-2 modes — which the within-sector blocks contain but G6 didn't examine — also carry positive kinetic eigenvalue?

## 2. Result

G6's A_λ = a(4λ²/Λ⁴ − 2/Λ²) is a function of the sector eigenvalue λ alone (the spectral action sees only the degenerate-sector eigenvalue), so it is **irrep-blind**: every irrep in a within-sector block — (1,1), (2,0), (0,2), … — inherits the same A_λ. At Λ²=6:

| n | (j_L,j_R) | λ | A_λ | (1,1) h=0 | (2,0)/(0,2) h=2 |
|:-:|:---------:|:-:|:---:|:---------:|:---------------:|
| 1 | (1, ½) | 2.5 | +0.127 | +0.127 | +0.127 |
| 2 | (3/2, 1) | 3.5 | +0.133 | +0.133 | +0.133 |
| 3 | (2, 3/2) | 4.5 | +0.066 | +0.066 | +0.066 |

The helicity-2 (TT-graviton-class) modes are present (GD-4) and positive (same A_λ).

## 3. Resolution + honest scope

- **Necessary condition extends to the TT sector.** G6's positive result is NOT confined to the scalar/trace (helicity-0) sector — the actual TT-graviton irrep class (helicity-2) is also positive. GD-4's concern is resolved at this level.
- **Honest:** the irrep-blindness that makes this "free" IS the G6 approximation. The **sufficient** condition — that TT modes have a kinetic structure *distinct* from the trace (Fierz–Pauli) — requires lifting the irrep-blind A_λ to an irrep-resolved second variation. That is exactly the multi-month G6, and it is where the irrep-blindness must break.

## 4. ① arc summary (GD-1 + GD-4 + GD-5)

- GD-1: the 9→2 polarization reduction is a continuum-limit (propinquity) phenomenon (discrete inner gauge is cross-sector).
- GD-4: G6 tested helicity-0 (1,1); TT graviton is helicity-2 (|Δj|=2).
- GD-5: helicity-2 modes are also positive (A_λ irrep-blind) → necessary condition holds for the TT sector.

Net: the necessary condition (positive modes in the TT-graviton irrep) is confirmed; the multi-month G6 is sharply characterized — (i) irrep-resolved second variation to distinguish TT from trace, (ii) the δD↔h_μν dictionary, (iii) the 9→2 gauge reduction in the propinquity limit.

## 5. Documentation
- Paper 28 §4.10: GD-4 paragraph extended with the GD-5 resolution.

## Files
- `debug/sprint_gd5_helicity2_modes.py`
