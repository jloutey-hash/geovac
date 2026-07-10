# Sprint GD-6 — convention probe: the Möbius is a genuine R-independent tip, not a subtraction artifact

**Date:** 2026-05-29. **Type:** ② follow-on. **Verdict:** the Möbius slope is **R-independent (bit-identical across R = 5, 7.5, 10, 15)**. The Δ_K = K_wedge − α·K_disk subtraction cleanly isolates a pure tip constant (bulk R² and perimeter R cancel); the Möbius-vs-SC difference is a genuine tip effect, **not** an artifact of the subtraction convention.

## 1. The question (from GD-3)

GD-3 flagged: is the Möbius-vs-SC difference partly an artifact of the Δ_K = K_wedge − α·K_disk (bulk-normalized) subtraction convention? Test: R-independence. The bulk-Weyl term (~R²/t) cancels in Δ_K by construction; if the perimeter (~R/√t) also cancels geometrically, Δ_K is a pure R-independent tip → slope R-stable → genuine tip. If an edge residual survives, the slope drifts with R → convention artifact.

## 2. Result (R-sweep, FD substrate, a=0.05, N_0=80, t=1)

| α | slope/Möbius at R=5 | R=7.5 | R=10 | R=15 |
|:-:|:-------------------:|:-----:|:----:|:----:|
| 2 | 1.0121 | 1.0121 | 1.0121 | 1.0121 |
| 3 | 0.9766 | 0.9766 | 0.9766 | 0.9766 |
| 5 | 0.9773 | 0.9773 | 0.9773 | 0.9773 |

**Bit-identical across a 3× range of R** (spread 0.0000). The Δ_K subtraction isolates an R-independent tip constant — both the bulk (R²) and the perimeter (R) terms cancel exactly. The Möbius is a genuine **tip** quantity.

## 3. Reading

- **NOT a convention/edge artifact.** The R-independence rules out the hypothesis that Möbius-vs-SC is a residual of the α·K_disk subtraction. The subtraction is clean; the tip is real.
- **Refined coefficient picture.** The slope/Möbius ratio is R-independent but ~2% off unity at these params (1.012, 0.977, 0.977) — i.e. the genuine tip is *close to* α/(2α−1) but the exact-coefficient match is t≈1/parameter-sharp (GD-2; the validation reached 0.03% at α=10 with larger params). So: the tip is a clean R-independent quantity; the Möbius form α/(2α−1) is its (t≈1-sharp) description.

## 4. The full ② arc (GD-2 + GD-3 + GD-6)

The Möbius **form** is robust on **three independent axes**:
- **t** (GD-2): slope is Möbius-shaped, never SC-shaped, at every t.
- **discretization** (GD-3): FD-azimuthal and exact-spectral agree to ~4 sig figs (substrate-class-universal).
- **R** (GD-6): bit-identical across R = 5→15 (genuine tip, not edge).

Only the *exact coefficient's* match to α/(2α−1) is t≈1-tuned (GD-2); the soft-IR "mechanism" is a t≈1 coincidence (GD-2, demoted). So the standing is sharp: **the Möbius is a genuine, R-independent, discretization-independent substrate-class tip effect with no continuum analog (Route A); its α/(2α−1) form is robust; its exact coefficient is a t≈1-sharp approximation; the continuum mechanism remains open.**

## 5. Documentation
- Paper 51 Möbius caveat: R-independence / genuine-tip finding added (form robust on three axes).

## Files
- `debug/sprint_gd6_moebius_convention.py`
