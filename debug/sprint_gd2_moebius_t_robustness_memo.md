# Sprint GD-2 — t-robustness audit of the Möbius α>1 mechanism

**Date:** 2026-05-29. **Type:** diagnostic / audit (curve-fit-audit discipline applied to the gravity sector). **Verdict:** the Möbius **functional form** (α/(2α−1)-shape, categorically distinct from continuum SC) is **robust in t**; the **exact coefficient** (the "sub-2% / 0.03%" precision) and the **soft-IR-fraction "mechanism"** (Route C) are both **t≈1-tuned**, not t-robust. The prior precision claims were single-t (a chosen "sweet spot"); the mechanism remains open.

## 1. Motivation

The α>1 conical slope was empirically locked at −(1/12)·α/(2α−1) (sub-2% across 6 α, 0.03% at α=10), and Route C "identified the mechanism" as soft_IR_frac(α) → 1/(2α) at sub-percent. **Both used a single t = 1.0** (the mode-decomp driver explicitly calls it the "sweet-spot t window"). The audit discipline demands a robustness check: a continuum-Weyl estimate gives soft_IR_frac ∼ 2√(t/π)/α — t-dependent. So: sweep t.

## 2. Part 1 — soft-IR fraction (Route C "mechanism")

Driver `debug/sprint_gd2_moebius_t_robustness.py` (spectral substrate, N_ρ=100, a=0.05, N_0=60). The relation (1−soft_IR_frac)·F should equal 1/2 if the mechanism is structural:

| α | t=0.25 | t=0.5 | t=1.0 | t=2.0 | t=4.0 |
|:-:|:------:|:-----:|:-----:|:-----:|:-----:|
| 2 | 0.589 | 0.555 | **0.504** | 0.427 | 0.309 |
| 3 | 0.553 | 0.532 | **0.500** | 0.450 | 0.369 |
| 5 | 0.529 | 0.517 | **0.499** | 0.470 | 0.422 |

Hits 1/2 only at t≈1; spread 0.11–0.28 across t. **The soft_IR_frac → 1/(2α) "mechanism" is a t≈1 coincidence, not a structural identity.** DEMOTED.

## 3. Part 2 — the Möbius slope itself

Driver `debug/sprint_gd2_moebius_slope_t.py` (FD substrate, N_ρ=150, a=0.05, N_0=80; reproduces the validation's t=1 result at ~1–2%). slope_meas = [K_wedge − α·K_disk]/(1/α−α):

| α | quantity | t=0.25 | t=0.5 | t=1.0 | t=2.0 | t=4.0 |
|:-:|:---------|:------:|:-----:|:-----:|:-----:|:-----:|
| 2 | slope/Möbius | 0.819 | 0.923 | **1.012** | 1.089 | 1.154 |
| 3 | slope/Möbius | 0.780 | 0.884 | **0.977** | 1.058 | 1.130 |
| 5 | slope/Möbius | 0.780 | 0.885 | **0.977** | 1.059 | 1.132 |

**Two findings:**
1. **Form is robust:** at every t and α, slope_meas is far closer to Möbius (−(1/12)α/(2α−1)) than to continuum SC (−(1/12)(1/α−α)). E.g. α=2: measured ≈ −0.05, Möbius = −0.0556, SC = **+0.125** (wrong sign, wildly off). The α>1 slope is unambiguously Möbius-shaped, not SC-shaped — robustly.
2. **Coefficient is t≈1-tuned:** slope/Möbius crosses 1 near t≈1 but drifts to ~0.78 (t=0.25) and ~1.13–1.15 (t=4). Spread ~0.33–0.35. The "sub-2% / 0.03%" precision was a t≈1 crossing, NOT a t-robust coefficient.

## 4. Synthesis

- **Robust (stands):** the α>1 conical slope follows the Möbius functional form α/(2α−1), categorically distinct from the continuum Sommerfeld–Cheeger form. This is the real, robust content of the α>1 work.
- **t≈1-tuned (corrected):** the exact coefficient match (sub-2%, 0.03%) and the soft_IR_frac → 1/(2α) "mechanism" are both sweet-spot artifacts at t≈1; they drift ±15–35% over t ∈ [0.25, 4].
- **Consistent with v3.20.0 per-t UV finding:** substrate tip extraction is genuinely t-dependent (per-t UV target = 1/(24πt)); clean continuum coefficients appear only at a sweet-spot t. The Möbius coefficient's t≈1-tuning is the same phenomenon.
- **Mechanism status:** OPEN. The soft-IR explanation is demoted to a t=1 coincidence. The form's robustness (Möbius vs SC) is unexplained at the continuum level — consistent with Route A's finding that no published continuum Möbius exists.

## 5. Honest framing

The prior work chose t=1 as a deliberate "sweet spot," so this is a sharpening, not a caught error: the audit quantifies that the clean Möbius *coefficient* is sweet-spot-specific while the *form* is robust. The right presentation is: **Möbius form = robust gravity-sector finding; exact coefficient + soft-IR mechanism = sweet-spot (t≈1) artifacts.**

## 6. Documentation
- Paper 51 Möbius paragraph: t-robustness caveat added (form robust, coefficient/mechanism t≈1-tuned).
- CLAUDE.md §3 dead-ends: soft_IR_frac → 1/(2α) "mechanism" demoted to t=1 coincidence.

## Files
- `debug/sprint_gd2_moebius_t_robustness.py` (soft-IR fraction sweep)
- `debug/sprint_gd2_moebius_slope_t.py` (slope sweep)
