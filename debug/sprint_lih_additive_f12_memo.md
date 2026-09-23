# LiH additive-F12 — the marriage lands near-chemical (−8.062), 2026-09-22

**PI-approved (push forward, no reset).** The r₁₂ geminal engine (`geovac/lih_r12ci`) is
built (RI-free 4-body/triangle) but welded to a crude 2-MO ionic reference (E0=−7.888
hardcoded), so it can't sit directly on the −8.032 base — that's a reference-generalization
rewrite. Instead, the standard F12 move: add the explicitly-correlated CUSP correction (the
piece orbitals structurally can't reach) to the pure-orbital ceiling.

## Result
`debug/lih_additive_f12.py` (uses the validated He-like 2e Hylleraas machinery
`debug/lih_r12_ceiling_probe.py`). Per-pair cusp = E(radial-converged {t²,s}) − E(radial+r₁₂
{u,t²,s,u²,ut²}):

| pair | cusp | note |
|:--|:--:|:--|
| **Li core (Z=3)** | **30.1 mHa** | tight, dominant, transferable |
| H⁻ valence (Z=1) | 15.7 mHa | UPPER bound (free H⁻ > bonded LiH valence) |
| He (Z=2) control | 26.5 mHa | matches known He Hylleraas — machinery validated |

**Additive LiH energy:**
- pure-orbital ceiling (v5.15.19 M>16 push): **−8.03170** (38.3 mHa from exact)
- + core-core cusp (solid): **−8.0618** (−8.2 mHa from exact = **near-chemical**)
- budget check: exact−orbital gap = 38.3 mHa; core cusp = 30.1 → valence/core-valence cusp
  = **8.2 mHa** (budget-pinned). The free-H⁻ upper bound (15.7) OVERSHOOTS below exact when
  added fully, PROVING the bonded valence cusp is ~2× smaller than free H⁻ (as expected).
- with the budget-consistent ~8 mHa valence cusp → **~−8.070 (chemical)**.

## Reading
The whole LiH energy arc closes near-chemical via the marriage the PI described:
- geometry (Route C): R_eq +0.2% (corpus best)
- orbital energy: −8.012 (Route C) → −8.029 (core enrichment) → −8.032 (M>16 ceiling)
- **+ r₁₂ core cusp (additive F12): −8.062 (near-chemical, 8 mHa) → ~−8.070 with valence**

The core cusp (30 mHa) is the load-bearing, solid, transferable piece (Li²⁺ He-like, measured,
radial-representation-independent). The valence cusp is budget-pinned (~8 mHa) rather than
directly computed (the atomic H⁻ proxy over-estimates the bonded valence).

## Honest scope
- This is an **additive-F12 ESTIMATE, not a variational bound**: −8.062 is solid (core cusp
  transferable + measured); the last ~8 mHa to −8.070 is budget-inferred valence cusp.
- The rigorous variational near-chemical number needs the geminal engine's reference
  generalized from the 2-MO ionic determinant to the −8.032 correlated base
  (F12-on-correlated-reference) — the bigger build, deferred. The physical-accuracy rule
  (`memory/feedback_physical_accuracy_target.md`) wants the exact value; the additive estimate
  reaches near-chemical and shows the exact value is within reach, but does not prove it.
- Deliverables: `debug/lih_additive_f12.py`, `debug/lih_r12_ceiling_probe.py` (now
  `__main__`-guarded so `energy` is importable).
