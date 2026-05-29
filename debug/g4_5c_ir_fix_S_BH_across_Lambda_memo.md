# Sprint G4-5c IR-fix — S_BH across Lambda with structural cutoff cures

**Date:** 2026-05-29
**Status:** NEGATIVE-ON-CUTOFF-CURE / SUBSTANTIVE-DIAGNOSTIC-POSITIVE
**Driver:** `debug/g4_5c_ir_fix_S_BH_across_Lambda.py`
**Data:** `debug/data/g4_5c_ir_fix_S_BH_across_Lambda.json`

## Summary

Three structural cutoff variants (Gaussian f(x) = e⁻ˣ, polynomial f(x) = e⁻ˣ², sharp
f(x) = Θ(1−x)) tested on the joint variable-warp + conical-defect substrate of
Sprint G4-5c, across six Lambda values {0.5, 1, 1.5, 2, 3, 5}, with extended t-grid
{0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50}.

**Headline finding (substantive, gate-unanticipated):**
The IR over-count at small Λ is **NOT** an IR tail-suppression failure.
The Mellin integral mass is dominated by the **UV peak of tip(t) at t ∈ [0.02, 0.5]**
(~88% of integral mass at Λ=0.5), where t·Λ² ≪ 1 for all tested Λ ≤ 1 — so the
cutoff function is essentially 1 in the dominant region regardless of which cutoff
family is used. Polynomial e⁻ˣ² and sharp Θ(1−x) do not improve over Gaussian e⁻ˣ
in the IR regime; they are slightly worse (because at small x the polynomial e⁻ˣ²
is closer to 1 than e⁻ˣ).

**Cutoff cure verdict: NEGATIVE.** None of the three cutoffs achieves S_BH within
factor 2 at all Λ. Best result is Gaussian with 2/6 cells passing.

**Per-Λ best cutoff:**

| Λ   | Gaussian | Polynomial | Sharp  | Best     | Continuum |
|----:|---------:|-----------:|-------:|:--------:|----------:|
| 0.5 |  48.44   |  50.46     | 50.98  | gaussian |  0.333    |
| 1.0 |  10.56   |  11.38     | 12.26  | gaussian |  1.333    |
| 1.5 |   4.04   |   4.44     |  4.50  | gaussian |  3.000    |
| 2.0 |   1.96   |   2.18     |  2.53  | gaussian |  5.333    |
| 3.0 |   0.64   |   0.71     |  0.89  | sharp    | 12.000    |
| 5.0 |   0.13   |   0.14     |  0.14  | sharp    | 33.333    |

Two cells pass the gate under Gaussian: Λ=2 (ratio 1.96, just inside factor 2) and
Λ=3 (ratio 0.64, comfortable). Λ=1.5 is borderline at 4.04. Λ=0.5 and Λ=1 are
massive over-counts; Λ=5 is a massive under-count (the t-grid does not reach small
enough t for Λ=5, so the integral misses UV contributions).

**Extended t-grid effect (vs G4-5c 8-pt grid):** The G4-5c 8-pt grid gave
ratio 29.08 / 5.81 / 0.85 at Λ ∈ {0.5, 1, 2}. The 13-pt extended grid (this sprint)
gives 48.44 / 10.56 / 1.96. The discrepancy is because the extended grid added
three small-t points {0.005, 0.01, 0.02} which contribute substantially to the
Mellin integral (the integrand at small t is large because tip(t) is rising and
1/log(t) has heavy weight). The G4-5c result was an artifact of the coarse t-grid
that happened to undersample the UV peak; the more honest extended-grid result
is the larger ratios reported here.

## Diagnostic — where does the integral mass live at Λ = 0.5?

Per-cell trapezoid contributions at Λ = 0.5, Gaussian cutoff:

| t_low  | t_high | integrand_lo | integrand_hi | cell contrib | cumul frac |
|-------:|-------:|-------------:|-------------:|-------------:|----------:|
| 0.005  | 0.010  |       2.17   |       4.70   |       2.38   |     0.074 |
| 0.010  | 0.020  |       4.70   |       6.55   |       3.90   |     0.195 |
| 0.020  | 0.050  |       6.55   |       7.88   |       6.61   |     0.400 |
| 0.050  | 0.100  |       7.88   |       7.84   |       5.45   |     0.569 |
| 0.100  | 0.200  |       7.84   |       6.74   |       5.05   |     0.726 |
| 0.200  | 0.500  |       6.74   |       4.09   |       4.96   |     0.879 |
| 0.500  | 1.000  |       4.09   |       2.18   |       2.17   |     0.947 |
| 1.000  | 2.000  |       2.18   |       0.98   |       1.10   |     0.981 |
| 2.000  | 5.000  |       0.98   |       0.23   |       0.56   |     0.998 |
| 5.000  | 10.0   |       0.23   |       0.04   |       0.10   |     1.001 |
| 10.0   | 20.0   |       0.04   |     2 × 10⁻³ |       0.015  |     1.001 |
| 20.0   | 50.0   |     2 × 10⁻³ |     2 × 10⁻⁵ |     9 × 10⁻⁴ |     1.001 |

**Mass partition at Λ = 0.5 (Gaussian):**
- Small-t (UV) regime t ≤ 0.05: 40% of integral
- Mid-t regime 0.05 ≤ t ≤ 0.5: 48% of integral
- Large-t (IR) regime t ≥ 0.5: 12% of integral

**The integral mass is concentrated where the Gaussian cutoff e⁻ᵗ·Λ² ≈ 1** — at
Λ = 0.5 and t = 0.05, x = t·Λ² = 0.0125, so e⁻ˣ ≈ 0.988. The cutoff is essentially
inactive across the bulk of the integral. The result is independent of cutoff
family for small Λ.

At Λ = 5 the situation reverses: the cutoff kicks in at t ≥ 1/25 = 0.04, suppressing
the entire UV peak that carries the mass. This is why the Λ=5 ratio drops to 0.13:
the substrate's tip-term mass is no longer counted.

## Structural reading — why the cutoff cure fails

The structural assumption underlying the Mellin integration S_BH = (1/2) ∫ (dt/t)
f(tΛ²) tip(t) is that the integrand has a **scale-separated structure**: the cutoff
f(tΛ²) selects the regime t ~ 1/Λ² and the tip(t) profile near that scale provides
the horizon contribution.

The discrete substrate violates this assumption. The tip(t) profile peaks at
t ~ 0.05 — a substrate-fixed UV scale (related to the lattice spacing a = 0.05 and
the warp-tip scale r_h = 2). It is **NOT** Λ-dependent. The Mellin integral for
Λ < 1 simply collects this fixed UV peak times log(t) measure, giving a
Λ-independent J ~ 30 + small corrections from the cutoff hitting the tail.

The continuum prediction S_BH = r_h²Λ²/3 has a clean Λ² scaling because in the
continuum the heat-kernel tip contribution itself scales as Λ². In the discrete
substrate, the tip-term magnitude is set by the substrate's own UV cutoff (a²)
and does not scale with Λ.

**This is a regime-of-validity finding, not a cutoff-engineering problem.**
The CC heat-trace Mellin integration framework requires Λ in a window where:

1. **UV bound:** Λ < 1/a (so the cutoff is inside the substrate's UV resolution
   — at a = 0.05, this gives Λ < 20).
2. **IR bound:** Λ > 1/R (so the cutoff is outside the substrate's IR cutoff
   — at R = 10, this gives Λ > 0.1).
3. **Tip-scale bound:** Λ should be comparable to 1/√(tip_peak_t).
   With tip-peak at t ≈ 0.05 → 1/√0.05 ≈ 4.5. So Λ in the **2 – 5 range**
   is where the substrate's UV peak in tip(t) is the regime the Gaussian
   cutoff is actually selecting.
4. **Sub-leading bound:** Λ < r_h⁻¹ × (substrate sub-leading scale) — at
   r_h = 2, this gives Λ < 0.5 × ... — and so on.

The G4-5c Λ=2 cell sits exactly in the sub-leading peak region (1/√0.05 ≈ 4.5
ish, Λ² = 4 → t·Λ² = 1 at t = 0.25). The Λ=3 cell of this sprint at ratio 0.64
sits even closer to the peak scale-matching. Smaller Λ misses the regime;
larger Λ over-cuts the substrate UV content.

## Three possible cures (NOT pursued here)

The diagnostic suggests three structurally-different cures, each beyond
sprint-scale:

### Cure 1: Sub-leading bulk subtraction

The bulk subtraction in G4-5c is K(α=1) — the linear-in-α leading bulk. The
extracted tip(t) still contains α-independent sub-leading bulk pieces (e.g., warp
corrections to Weyl coefficients) that the simple bulk subtraction misses. A
proper refinement would subtract a refined bulk including sub-leading Weyl terms:

  tip(t) → tip(t) - (sub-leading Weyl coefficient)·a²·K(α=1)

The diagnostic for this would be a reciprocal-cancellation test analogous to
G4-4f: rescaling Λ should produce predictable changes in the J integral if the
bulk subtraction is correct, and unpredictable changes if it is incomplete.

### Cure 2: Substrate UV refinement

Push a → 0 (smaller lattice spacing) so the tip(t) peak moves to t < a² and
the Λ-window for clean S_BH extraction expands. At a = 0.05, the tip peak is
at t ≈ 0.05 ≈ a; at a = 0.01, the tip peak would shift to t ~ 0.01, allowing
Λ up to ~10 with clean scaling. This is a multi-week compute (N_rho would
need to scale from 200 to 1000).

### Cure 3: Effective-horizon-area rescaling

The continuum prediction r_h²Λ²/3 may need rescaling if the discrete substrate's
"effective horizon area" differs from 4πr_h² due to discretization-of-conical-tip
effects. The G4-2 / G4-3 sprints documented that proper conical-defect
discretization (Sommerfeld-Cheeger (1/12)(1/α - α) factor) is not naively
realized by the wedge-azimuthal construction. A proper area-calibration sprint
would compute A_eff(N_φ, a) from the substrate independently and use that
instead of 4πr_h².

## What was learned

The cutoff variants are mathematically inadequate to address the IR-overcount;
the issue lives at a deeper structural level (tip-term profile vs Λ-scaling).
This is a **clean negative on the proposed cure** with **substantive new
diagnostic content**: the failure is localized to the regime-of-validity of
the Mellin integration framework on a fixed-UV-resolution substrate.

The G4-5c Λ=2 PARTIAL closure was already at the favorable edge of this regime;
the new t-grid shows Λ=2 actually sits closer to ratio 1.96 (not 0.85) when
the UV cells are properly sampled. The G4-5c PARTIAL verdict therefore stays
within the same PARTIAL band — the more honest extended-grid extraction gives
ratio 1.96 (still inside factor 2), with the Λ=3 cell at 0.64 as a new
favorable data point not in the G4-5c sprint.

## Per-Λ best window — recommendation

Recommended Λ window for clean S_BH extraction at fixed (N_rho = 200, a = 0.05,
r_h = 2): **Λ ∈ [2, 3]** under Gaussian cutoff. This gives the only two cells
within factor 2 of continuum:

  - Λ = 2: ratio 1.96 (within factor 2)
  - Λ = 3: ratio 0.64 (within factor 2)

The Λ = 0.5, 1.0, 1.5 cells overcount because the substrate's UV peak in tip(t)
is not cut by the small Λ; the Λ = 5 cell undercounts because the substrate's
tip-term UV content is cut by the strong Λ before it can contribute.

## Honest scope

This is a sprint-scale diagnostic (one driver, three cutoffs, six Λ). It rules out
the simple "switch cutoff family" cure with a clean structural reading. Substantive
cures (sub-leading bulk subtraction; smaller a; effective-horizon-area calibration)
are multi-week sprints flagged as named follow-ons in the G4-5 multi-month plan.

The construction (joint-warp + conical-defect at F11 closure level) is unchanged
from G4-5c; F6 falsifiers remain bit-exact. This sprint operates on the post-
processing Mellin integration only.

## Verdict by gate

| Gate criterion | Result |
|:--------------|:------:|
| POSITIVE (S_BH within factor 2 at all Λ in {0.5, 1, 2, 3, 5}) | FAIL |
| PARTIAL (improvement over G4-5c at Λ=0.5, 1.0) | FAIL — extended grid actually worsens |
| PARTIAL (ratio < factor 2 at 3 of 5 Λ values) | FAIL — best is 2/6 (Gaussian at Λ=2, 3) |
| NEGATIVE (cure does not improve IR cells) | MET |

**Verdict: NEGATIVE on cure; POSITIVE on diagnostic.**

Substantive new content: the IR-overcount is structural (substrate UV peak
versus Λ scaling), not a tail-suppression issue. Three structurally-different
cures named for follow-on sprints (sub-leading bulk subtraction; smaller a;
effective-horizon-area calibration).

## Files

- `debug/g4_5c_ir_fix_S_BH_across_Lambda.py` — driver + three cutoff variants + per-cell diagnostic
- `debug/data/g4_5c_ir_fix_S_BH_across_Lambda.json` — full S_BH table + diagnostic
- `debug/g4_5c_ir_fix_S_BH_across_Lambda_memo.md` — this memo

No `geovac/gravity/` modifications. No `tests/` modifications. Uses `np.trapezoid`
throughout. Imports `JointWarpConicalDirac` from the G4-5c driver via `sys.path`
insertion.
