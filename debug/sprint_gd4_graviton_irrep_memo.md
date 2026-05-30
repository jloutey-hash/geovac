# Sprint GD-4 — graviton irrep diagnostic: G6 tested helicity-0, the TT graviton is helicity-2

**Date:** 2026-05-29. **Type:** diagnostic (① follow-on, rep-theory). **Verdict:** SHARPENING (raises a concrete must-check for the multi-month G6). G6's "(1,1) graviton candidate" is the helicity-0 (|j_L−j_R|=0, scalar/trace-class) sector; the physical TT graviton is helicity ±2 (|j_L−j_R|=2), which the within-sector blocks contain but G6 did NOT examine. Whether G6's (1,1) maps to the physical graviton requires the CC δD↔h_μν dictionary — the unresolved multi-month piece.

## 1. The rep-theory fact

Helicity of an SO(4)=SU(2)_L×SU(2)_R irrep (j_L,j_R) is h = |j_L − j_R|. The TT spin-2 graviton is helicity ±2 → **|j_L−j_R| = 2**. The lowest TT tensor harmonics on S³ are (2,0)⊕(0,2); generally (k+2,k)⊕(k,k+2). The **(1,1)** irrep has |1−1| = 0 → **helicity 0** → scalar/trace-class, NOT the TT graviton.

## 2. What G6 tested

G6 searched the within-sector blocks sector_i⊗sector_i (CH spinor sectors (j_L,j_R)=((n+1)/2, n/2)) for the **(1,1)** irrep and found positive eigenvalues. But those blocks also contain the **(2,0)⊕(0,2)** (helicity-2) reps — verified present at n=1,2,3,4 — which G6 did NOT examine. So the G6 "necessary-condition" positive result was on the **helicity-0 (scalar/trace) sector**, not the helicity-2 TT-graviton sector.

## 3. The honest caveat (why this is a sharpening, not a refutation)

In Connes–Chamseddine, the physical graviton enters through the δD-perturbation via the δD↔h_μν dictionary, which can shift irrep assignments. The **local** graviton polarization space IS the traceless-symmetric 2-tensor = (1,1); the **global** TT harmonic modes on S³ are |Δj|=2. Which of these G6's bilinear (1,1) corresponds to depends on the dictionary — exactly the unresolved multi-month piece. So GD-4 does not refute G6; it identifies a specific gap: G6 cannot claim to have tested the TT graviton sector until either (a) the |Δj|=2 δD-modes are examined, or (b) the dictionary is established showing the (1,1) bilinear carries the physical graviton.

## 4. Combined ① picture (GD-1 + GD-4)

The multi-month G6 has two concrete must-dos, now sharply identified:
1. **(GD-4)** Examine the helicity-2 (|Δj|=2) δD-modes — not just (1,1) — and fix the δD↔h_μν dictionary so the test is on the physical graviton sector.
2. **(GD-1)** Realize the 9→2 polarization reduction in the continuum (propinquity) limit, since the discrete inner-automorphism gauge is cross-sector (vanishes within a block) while the continuum diffeomorphism gauge acts within the tensor.

Both are continuum-limit / dictionary tasks tied to the Paper 38 propinquity framework — consistent with the deferral of full Fierz–Pauli to multi-month G6.

## 5. Honest scope

The helicity/irrep facts are solid rep theory. The interpretation (whether G6 tested the "wrong" sector) is appropriately hedged by the δD-dictionary caveat. Net: a concrete, honest sharpening of what the multi-month G6 must verify, NOT a claim that G6's necessary-condition result is wrong.

## 6. Documentation
- Paper 28 §4.10: helicity/irrep paragraph added after the GD-1 paragraph.

## Files
- `debug/sprint_gd4_graviton_irrep.py`
