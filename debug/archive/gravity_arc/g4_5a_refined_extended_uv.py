"""Sprint G4-5a refined -- Tip-only replica integration with extended UV t-grid.

Refines the G4-5a first move (debug/g4_5a_first_move_tip_replica.py) by
extending the t-grid down to the substrate UV cutoff t = a^2 = 0.0025,
capturing the contribution from the previously-unsampled UV window
[a^2, 0.1] where the Gaussian cutoff exp(-t Lambda^2) concentrates weight
at large Lambda.

Background
----------
- G4-5a first move sampled only t >= 0.1. At Lambda = 2.0, the Gaussian
  cutoff peaks at t ~ 1/Lambda^2 = 0.25, but exp(-t Lambda^2) decays
  slowly enough that significant weight extends below t = 0.1.
- The discrete/continuum ratio at G4-5a baseline:
    Lambda = 0.5: 0.37
    Lambda = 1.0: 0.26
    Lambda = 1.5: 0.18
    Lambda = 2.0: 0.13
  The ratio drops as Lambda grows because the missing UV window has more
  weight at large Lambda.
- T2 G4-3d-UV identified an angular truncation overshoot at small t
  (ratio 4/pi^2 ~ 0.405 between discrete and continuum azimuthal
  Laplacian at the truncation edge). Tip term may be affected, but
  G4-4f cross-t analysis already showed tip recovery 76% at t = 0.1
  declining toward UV.

Method
------
1. Extend t-grid: t_extended = {0.0025, 0.005, 0.01, 0.02, 0.05, 0.1,
   0.2, 0.5, 1, 2, 5, 10} (12 points spanning UV to IR).
2. Recompute K_wedge(alpha = 1 +/- eps) and K_disk at each t using the
   same panel as G4-5a: R = 10, a = 0.05, N_rho = 200, N_0 = 120,
   eps = 12/120 = 0.1.
3. Compute tip term Delta'(t) = dK/dalpha - K_disk at each t and
   tabulate vs +1/6 continuum prediction.
4. Integrate over t with Gaussian cutoff exp(-t Lambda^2) for
   Lambda in {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}.
5. Compare to:
   - Rough continuum Mellin moment M_0(Lambda) ~ -log(Lambda^2 * a^2)
   - Improved Mellin moment: exact Gaussian Mellin on [a^2, R^2]
     = E1(a^2 * Lambda^2) - E1(R^2 * Lambda^2)
     where E1 is the exponential integral E_1(x) = int_x^inf exp(-u)/u du
6. Report discrete/continuum ratio improvement vs G4-5a baseline.

Decision gate (per sprint plan)
-------------------------------
- POSITIVE: ratio > 0.5 at Lambda in [0.5, 2.0]
- PARTIAL: ratio improves over G4-5a baseline but does not reach 0.5
- NEGATIVE: extension does not improve recovery (UV overshoot dominates)
"""

import json
from pathlib import Path

import numpy as np
from scipy.special import exp1  # E_1(x) = int_x^inf exp(-u)/u du

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteWedgeDirac,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_refined_extended_uv.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def gaussian_mellin_E1(Lambda: float, t_min: float, t_max: float) -> float:
    """Exact Gaussian Mellin moment.

    M_0(Lambda; t_min, t_max) = int_{t_min}^{t_max} (dt/t) exp(-t Lambda^2)

    Substitute u = t Lambda^2:
    = int_{t_min Lambda^2}^{t_max Lambda^2} (du/u) exp(-u)
    = E1(t_min * Lambda^2) - E1(t_max * Lambda^2)
    """
    return float(exp1(t_min * Lambda**2) - exp1(t_max * Lambda**2))


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-5a refined -- Tip-only replica integration with UV extension")
    print("=" * 72)

    # Substrate (matches G4-5a / G4-4f sweet-spot panel)
    R = 10.0; a = 0.05; N_rho = 200; N_0 = 120
    print(f"\n[Setup] R={R}, a={a}, N_rho={N_rho}, N_0={N_0}")
    print(f"  Substrate UV cutoff: t_UV ~ a^2 = {a**2}")
    print(f"  Substrate IR cutoff: t_IR ~ R^2 = {R**2}")

    # Replica eps: matches G4-4f best window at recovery 96.69%
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps = (alpha_plus - alpha_minus) / 2
    print(f"\n[Replica parameters]")
    print(f"  alpha_+ = {alpha_plus}, alpha_- = {alpha_minus}, eps = {eps}")

    # Extended t-grid: spans substrate UV (a^2 = 0.0025) to IR (R^2 = 100)
    # 12 points, log-spaced enough to handle Gaussian Mellin
    t_grid = np.array([0.0025, 0.005, 0.01, 0.02, 0.05, 0.1,
                       0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    print(f"\n[Extended t grid] {t_grid.tolist()}")
    print(f"  Spans UV t = a^2 = {a**2} to t = 10")
    print(f"  Previous G4-5a grid started at t = 0.1; "
          f"this extends by 4 UV points below.")

    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    wedge_plus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
    )
    wedge_minus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
    )

    print(f"\n[Computing K(t) at each t in extended grid...]")
    K_disk_arr = np.array([disk.heat_trace(t) for t in t_grid])
    K_plus_arr = np.array([wedge_plus.heat_trace(t) for t in t_grid])
    K_minus_arr = np.array([wedge_minus.heat_trace(t) for t in t_grid])

    # Replica derivative at each t
    dK_dalpha = (K_plus_arr - K_minus_arr) / (alpha_plus - alpha_minus)
    # Tip-only contribution = dK/dalpha - K_disk (subtract bulk linear-in-alpha part)
    tip_term = dK_dalpha - K_disk_arr

    print(f"\n  {'t':>7}  {'K_disk':>13}  {'dK/dalpha':>13}  {'tip term':>12}  "
          f"{'tip/(1/6)':>10}")
    print("  " + "-" * 70)
    for i, t in enumerate(t_grid):
        rec = tip_term[i] / (1/6)
        marker = " *NEW UV*" if t < 0.1 else ""
        print(f"  {t:>7.4f}  {K_disk_arr[i]:>13.4f}  {dK_dalpha[i]:>13.4f}  "
              f"{tip_term[i]:>+12.6f}  {rec:>10.4f}{marker}")

    results["heat_trace_data"] = {
        "t_grid": t_grid.tolist(),
        "K_disk": K_disk_arr.tolist(),
        "dK_dalpha": dK_dalpha.tolist(),
        "tip_term": tip_term.tolist(),
        "tip_recovery_vs_1_over_6": (tip_term / (1/6)).tolist(),
    }

    # ------------------------------------------------------------------
    # Step 2: Integrate over t with Gaussian cutoff (extended grid)
    # ------------------------------------------------------------------
    print("\n[Step 2] Integrate tip term with Gaussian cutoff f(x) = exp(-x):")
    print("  Now with extended UV t-grid starting at t = a^2 = 0.0025.")
    print()

    Lambda_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    print(f"  Lambda values: {Lambda_values}")
    print(f"  Integral: J(Lambda) = int_{{t_min}}^{{t_max}} (dt/t) exp(-t Lambda^2) tip(t)")
    print(f"  S_tip(Lambda) = +J/2  (sign convention: dI_E/dalpha = -(1/2) S_CC; "
          f"S_tip = +J/2)")
    print()

    # Use trapezoidal rule on log(t) grid
    # integral (dt/t) g(t) = integral d(log t) g(t)
    log_t = np.log(t_grid)
    integrate_results = {}

    print(f"  {'Lambda':>8}  {'J integral':>14}  {'S_tip (disc)':>14}  "
          f"{'M_0 rough':>10}  {'S_tip rough':>12}  "
          f"{'M_0 exact':>12}  {'S_tip exact':>13}")
    print("  " + "-" * 95)

    for Lambda in Lambda_values:
        f_cutoff = np.exp(-t_grid * Lambda**2)
        # integral d(log t) of integrand
        integrand_log = f_cutoff * tip_term
        J = float(np.trapezoid(integrand_log, log_t))
        S_tip_disc = +J / 2

        # Rough Mellin estimate (G4-5a-style): M_0 ~ -log(Lambda^2 * a^2)
        M_0_rough = -np.log(Lambda**2 * a**2) if Lambda**2 * a**2 > 0 else 0
        S_tip_rough = +(1/12) * M_0_rough

        # Improved Mellin: exact Gaussian moment on [a^2, R^2]
        M_0_exact = gaussian_mellin_E1(Lambda, a**2, R**2)
        S_tip_exact = +(1/12) * M_0_exact

        ratio_rough = S_tip_disc / S_tip_rough if abs(S_tip_rough) > 1e-12 else float("inf")
        ratio_exact = S_tip_disc / S_tip_exact if abs(S_tip_exact) > 1e-12 else float("inf")

        integrate_results[Lambda] = {
            "J": J,
            "S_tip_disc": S_tip_disc,
            "M_0_rough": M_0_rough,
            "S_tip_rough": S_tip_rough,
            "ratio_rough": ratio_rough,
            "M_0_exact": M_0_exact,
            "S_tip_exact": S_tip_exact,
            "ratio_exact": ratio_exact,
        }
        print(f"  {Lambda:>8.2f}  {J:>14.6e}  {S_tip_disc:>+14.6e}  "
              f"{M_0_rough:>10.4f}  {S_tip_rough:>+12.6f}  "
              f"{M_0_exact:>12.6f}  {S_tip_exact:>+13.6f}")

    results["integrate_results"] = {
        str(k): v for k, v in integrate_results.items()
    }

    # ------------------------------------------------------------------
    # Step 3: Ratio improvement table vs G4-5a baseline
    # ------------------------------------------------------------------
    print("\n[Step 3] Ratio improvement vs G4-5a baseline:")
    print()
    # G4-5a reported these ratios using the rough Mellin (M_0 ~ -log(Lambda^2 a^2))
    g4_5a_baseline_rough = {0.5: 0.37, 1.0: 0.26, 1.5: 0.18, 2.0: 0.13}

    print(f"  Baseline (G4-5a, rough Mellin): ratios at Lambda = "
          f"0.5/1.0/1.5/2.0: 0.37/0.26/0.18/0.13")
    print()
    print(f"  {'Lambda':>8}  {'G4-5a rough':>12}  "
          f"{'G4-5a-ref rough':>16}  {'delta rough':>12}  "
          f"{'G4-5a-ref exact':>16}")
    print("  " + "-" * 75)

    ratio_improvement = {}
    for Lambda in Lambda_values:
        baseline = g4_5a_baseline_rough.get(Lambda, None)
        refined_rough = integrate_results[Lambda]["ratio_rough"]
        refined_exact = integrate_results[Lambda]["ratio_exact"]
        if baseline is not None:
            delta = refined_rough - baseline
            print(f"  {Lambda:>8.2f}  {baseline:>12.4f}  {refined_rough:>16.4f}  "
                  f"{delta:>+12.4f}  {refined_exact:>16.4f}")
        else:
            print(f"  {Lambda:>8.2f}  {'(new)':>12}  {refined_rough:>16.4f}  "
                  f"{'(new)':>12}  {refined_exact:>16.4f}")
        ratio_improvement[Lambda] = {
            "baseline_rough": baseline,
            "refined_rough": refined_rough,
            "refined_exact": refined_exact,
            "delta_rough": (refined_rough - baseline) if baseline is not None else None,
        }

    results["ratio_improvement_vs_g4_5a"] = {
        str(k): v for k, v in ratio_improvement.items()
    }

    # ------------------------------------------------------------------
    # Step 4: Decision gate
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("\n[Step 4] Decision gate:")
    print()
    print(f"  POSITIVE: ratio > 0.5 at Lambda in [0.5, 2.0]")
    print(f"  PARTIAL: ratio improves over G4-5a baseline but does not reach 0.5")
    print(f"  NEGATIVE: extension does not improve recovery")
    print()

    # Use rough-Mellin ratios for apples-to-apples comparison with G4-5a baseline
    inrange_Lambdas = [L for L in Lambda_values if 0.5 <= L <= 2.0]
    inrange_rough_ratios = [integrate_results[L]["ratio_rough"]
                             for L in inrange_Lambdas]
    inrange_exact_ratios = [integrate_results[L]["ratio_exact"]
                             for L in inrange_Lambdas]
    inrange_baseline_ratios = [g4_5a_baseline_rough[L] for L in inrange_Lambdas]

    positive_rough = all(r > 0.5 for r in inrange_rough_ratios)
    positive_exact = all(r > 0.5 for r in inrange_exact_ratios)
    improves_all = all(
        integrate_results[L]["ratio_rough"] > g4_5a_baseline_rough[L]
        for L in inrange_Lambdas
    )

    min_rough = min(inrange_rough_ratios)
    min_exact = min(inrange_exact_ratios)
    avg_improvement_rough = float(np.mean([
        integrate_results[L]["ratio_rough"] - g4_5a_baseline_rough[L]
        for L in inrange_Lambdas
    ]))

    print(f"  In-range Lambdas tested: {inrange_Lambdas}")
    print(f"  Refined rough ratios: {[f'{r:.4f}' for r in inrange_rough_ratios]}")
    print(f"  Refined exact ratios: {[f'{r:.4f}' for r in inrange_exact_ratios]}")
    print(f"  Baseline rough ratios: {[f'{r:.4f}' for r in inrange_baseline_ratios]}")
    print()
    print(f"  min(rough ratio) on [0.5, 2.0]: {min_rough:.4f}")
    print(f"  min(exact ratio) on [0.5, 2.0]: {min_exact:.4f}")
    print(f"  Avg improvement (rough) vs G4-5a: {avg_improvement_rough:+.4f}")
    print()
    print(f"  All in-range Lambdas rough > 0.5: {positive_rough}")
    print(f"  All in-range Lambdas exact > 0.5: {positive_exact}")
    print(f"  All in-range Lambdas improved vs G4-5a (rough): {improves_all}")

    if positive_rough or positive_exact:
        verdict = "POSITIVE-G4-5a-REFINED-EXTENDED-UV"
        msg = ("Extended UV t-grid lifts discrete/continuum ratio above 0.5 on "
               "at least one Mellin convention across Lambda in [0.5, 2.0]. "
               "Tip-only replica integration now captures substrate UV "
               "contribution previously missing.")
    elif improves_all:
        verdict = "PARTIAL-G4-5a-REFINED-EXTENDED-UV"
        msg = ("Extended UV t-grid improves discrete/continuum ratio vs G4-5a "
               "baseline at every in-range Lambda but does not reach 0.5. "
               "Residual gap likely tied to UV high-m angular overshoot "
               "(T2 G4-3d-UV finding: 4/pi^2 ratio at angular truncation edge).")
    else:
        verdict = "NEGATIVE-G4-5a-REFINED-EXTENDED-UV"
        msg = ("Extended UV t-grid does not improve recovery; UV overshoot or "
               "other substrate effect dominates. G4-5b should use a different "
               "regularization strategy (e.g. dedicated UV tip subtraction).")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["min_ratio_rough"] = min_rough
    results["min_ratio_exact"] = min_exact
    results["avg_improvement_rough"] = avg_improvement_rough
    results["positive_rough"] = positive_rough
    results["positive_exact"] = positive_exact
    results["improves_all_vs_g4_5a"] = improves_all
    results["g4_5a_baseline_rough"] = g4_5a_baseline_rough

    # ------------------------------------------------------------------
    # Step 5: UV recovery diagnostic
    # ------------------------------------------------------------------
    print("\n[Step 5] UV recovery diagnostic (tip term per t):")
    print()
    print(f"  How well does the discrete tip term track the continuum +1/6 "
          f"prediction in the UV?")
    print()
    uv_indices = [i for i, t in enumerate(t_grid) if t < 0.1]
    for i in uv_indices:
        rec = tip_term[i] / (1/6)
        print(f"    t = {t_grid[i]:.4f}: tip = {tip_term[i]:+.4f}, "
              f"recovery vs 1/6 = {rec*100:.1f}%")

    if len(uv_indices) > 0:
        uv_recoveries = [tip_term[i] / (1/6) for i in uv_indices]
        results["uv_diagnostic"] = {
            "uv_t_values": [float(t_grid[i]) for i in uv_indices],
            "uv_tip_terms": [float(tip_term[i]) for i in uv_indices],
            "uv_recoveries": [float(r) for r in uv_recoveries],
            "min_uv_recovery": float(min(uv_recoveries)),
            "max_uv_recovery": float(max(uv_recoveries)),
        }

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
