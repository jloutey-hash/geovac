"""Sprint G4-5c IR-fix -- S_BH extraction across multiple Lambda with structural cutoff cures.

Background: G4-5c (`debug/g4_5c_joint_warp_conical_memo.md`) extracted S_BH at
Lambda=2 within factor 0.85 of continuum r_h^2 Lambda^2 / 3, but Lambda=0.5 gave
29x over and Lambda=1.0 gave 5.8x over. The mechanism is incomplete suppression
of the IR tail of tip(t) by the Gaussian cutoff f(x) = exp(-x) at small Lambda.

Sprint goal: improve discrete/continuum ratio across multiple Lambda by testing
three structural cures:

  (A) Gaussian f(x) = exp(-x)           [baseline, matches G4-5c]
  (B) Polynomial f(x) = exp(-x^2)        [faster IR decay]
  (C) Sharp f(x) = Theta(1 - x)           [strict UV/IR window]

Method:
  1. Recompute dK/dalpha and tip(t) at extended t-grid spanning UV to IR.
  2. For each Lambda in {0.5, 1, 1.5, 2, 3, 5}, integrate three cutoff variants.
  3. Report S_BH^discrete / (r_h^2 Lambda^2 / 3) ratio for each (Lambda, f) cell.
  4. Diagnostic per-cell integrand contribution to identify the dominant t.

Convention: Cutoffs f(t Lambda^2) act on x = t * Lambda^2. The Mellin integral
is S_BH = (1/2) integral (dt/t) f(t Lambda^2) tip(t), in log(t) coordinates
as trapezoid sum over the t-grid.

DO NOT modify code in geovac/gravity/ or tests/. Use np.trapezoid (not np.trapz).
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np

# Make debug/ importable so we can pull in the G4-5c class
_DEBUG_DIR = Path(__file__).parent
if str(_DEBUG_DIR) not in sys.path:
    sys.path.insert(0, str(_DEBUG_DIR))

from g4_5c_joint_warp_conical import JointWarpConicalDirac
from geovac.gravity.warped_dirac import S2DiracSpectrum


OUT_JSON = Path(__file__).parent / "data" / "g4_5c_ir_fix_S_BH_across_Lambda.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# ============================================================================
# Cutoff variants
# ============================================================================


def cutoff_gaussian(x: np.ndarray) -> np.ndarray:
    """f(x) = exp(-x). Standard CC cutoff (G4-5c baseline)."""
    return np.exp(-x)


def cutoff_polynomial(x: np.ndarray) -> np.ndarray:
    """f(x) = exp(-x^2). Faster IR decay; matches G4-5d-refined polynomial regime."""
    return np.exp(-(x ** 2))


def cutoff_sharp(x: np.ndarray) -> np.ndarray:
    """f(x) = Theta(1 - x). Strict UV/IR window."""
    return np.where(x <= 1.0, 1.0, 0.0)


CUTOFFS = {
    "gaussian": cutoff_gaussian,
    "polynomial": cutoff_polynomial,
    "sharp": cutoff_sharp,
}


# ============================================================================
# Tip term computation
# ============================================================================


def compute_tip_term(
    t_grid: np.ndarray,
    N_rho: int,
    a: float,
    N_0: int,
    k_step: int,
    r_h: float,
    l_max: int,
) -> Dict[str, np.ndarray]:
    """Compute dK/dalpha and tip(t) at extended t-grid via central FD at alpha=1.

    Mirrors G4-5c replica setup: alpha_+ = (N_0 + k_step)/N_0, alpha_- = (N_0 - k_step)/N_0.
    """
    sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0

    joint_plus = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus, r_h=r_h, sphere=sphere,
    )
    joint_minus = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus, r_h=r_h, sphere=sphere,
    )
    joint_alpha1 = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0, r_h=r_h, sphere=sphere,
    )

    K_alpha1 = np.zeros(len(t_grid))
    K_plus = np.zeros(len(t_grid))
    K_minus = np.zeros(len(t_grid))

    for i, t in enumerate(t_grid):
        K_alpha1[i] = joint_alpha1.heat_trace(float(t))
        K_plus[i] = joint_plus.heat_trace(float(t))
        K_minus[i] = joint_minus.heat_trace(float(t))

    dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
    tip_term = dK_dalpha - K_alpha1

    return {
        "K_alpha1": K_alpha1,
        "K_plus": K_plus,
        "K_minus": K_minus,
        "dK_dalpha": dK_dalpha,
        "tip_term": tip_term,
        "alpha_plus": alpha_plus,
        "alpha_minus": alpha_minus,
        "eps_alpha": (alpha_plus - alpha_minus) / 2,
    }


# ============================================================================
# Mellin integration with arbitrary cutoff
# ============================================================================


def mellin_integral(
    t_grid: np.ndarray,
    tip_term: np.ndarray,
    Lambda: float,
    cutoff_fn: Callable[[np.ndarray], np.ndarray],
) -> Dict[str, float]:
    """Compute S_BH = (1/2) integral (dt/t) f(t Lambda^2) tip(t).

    Integrand in log(t) coordinates is f * tip. Trapezoid sum.

    Also returns per-cell integrand contributions for diagnostic.
    """
    x = t_grid * (Lambda ** 2)
    f_cutoff = cutoff_fn(x)
    integrand_log = f_cutoff * tip_term

    log_t = np.log(t_grid)
    J = float(np.trapezoid(integrand_log, log_t))
    S_BH = J / 2.0

    return {
        "J_integral": J,
        "S_BH_extracted": S_BH,
        "cutoff_values": f_cutoff,
        "integrand_log_t": integrand_log,
    }


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    results: Dict = {}
    print("=" * 76)
    print("Sprint G4-5c IR-fix -- S_BH across Lambda with structural cutoff cures")
    print("=" * 76)

    # ------------------------------------------------------------------
    # Substrate parameters (match G4-5c baseline)
    # ------------------------------------------------------------------
    N_rho = 200
    a = 0.05
    N_0 = 120
    r_h = 2.0
    l_max = 3
    k_step = 12
    R = N_rho * a

    print(f"\n[Setup -- matches G4-5c baseline]")
    print(f"  R = {R}, N_rho = {N_rho}, a = {a}")
    print(f"  N_0 = {N_0}, r_h = {r_h}, l_max = {l_max}")
    print(f"  k_step = {k_step}, eps = {k_step/N_0}")

    results["substrate"] = {
        "N_rho": N_rho, "a": a, "N_0": N_0, "r_h": r_h,
        "l_max": l_max, "k_step": k_step, "R": R,
    }

    # ------------------------------------------------------------------
    # Extended t-grid covering UV (a^2 = 0.0025) through IR (~R^2 = 100)
    # ------------------------------------------------------------------
    t_grid = np.array([
        0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
        1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
    ])
    print(f"\n[Extended t-grid] N_t = {len(t_grid)}")
    print(f"  t range: [{t_grid.min()}, {t_grid.max()}]")
    print(f"  UV bound: a^2 = {a**2}, IR bound: R^2 = {R**2}")

    results["t_grid"] = t_grid.tolist()

    # ------------------------------------------------------------------
    # Compute tip term at extended t-grid (this is the expensive step)
    # ------------------------------------------------------------------
    print(f"\n[Computing K_joint and tip(t) at extended t-grid...]")
    print(f"  Each cell needs 3 heat-trace evaluations (alpha=+, -, 1)")
    print(f"  Total: {3 * len(t_grid)} heat-trace evaluations")
    print()

    tip_data = compute_tip_term(
        t_grid=t_grid, N_rho=N_rho, a=a, N_0=N_0,
        k_step=k_step, r_h=r_h, l_max=l_max,
    )

    K_alpha1 = tip_data["K_alpha1"]
    dK_dalpha = tip_data["dK_dalpha"]
    tip_term = tip_data["tip_term"]

    print(f"  {'t':>8}  {'K(alpha=1)':>14}  {'dK/dalpha':>14}  {'tip(t)':>14}")
    print("  " + "-" * 60)
    for i, t in enumerate(t_grid):
        print(f"  {t:>8.4f}  {K_alpha1[i]:>14.4f}  {dK_dalpha[i]:>14.4f}  "
              f"{tip_term[i]:>14.4f}")

    results["tip_data"] = {
        "K_alpha1": K_alpha1.tolist(),
        "K_plus": tip_data["K_plus"].tolist(),
        "K_minus": tip_data["K_minus"].tolist(),
        "dK_dalpha": dK_dalpha.tolist(),
        "tip_term": tip_term.tolist(),
        "alpha_plus": tip_data["alpha_plus"],
        "alpha_minus": tip_data["alpha_minus"],
        "eps_alpha": tip_data["eps_alpha"],
    }

    # ------------------------------------------------------------------
    # Mellin integration with three cutoff variants, across Lambda
    # ------------------------------------------------------------------
    Lambda_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    cutoff_names = ["gaussian", "polynomial", "sharp"]

    print("\n" + "=" * 76)
    print("[S_BH extraction] Three cutoff variants x six Lambda values")
    print("=" * 76)

    S_BH_table: Dict[str, Dict[str, Dict]] = {}

    for cutoff_name in cutoff_names:
        cutoff_fn = CUTOFFS[cutoff_name]
        print(f"\n[Cutoff: {cutoff_name}]")
        if cutoff_name == "gaussian":
            print(f"  f(x) = exp(-x)  -- baseline (G4-5c convention)")
        elif cutoff_name == "polynomial":
            print(f"  f(x) = exp(-x^2)  -- faster IR decay")
        elif cutoff_name == "sharp":
            print(f"  f(x) = Theta(1 - x)  -- strict UV/IR window")

        print(f"\n  {'Lambda':>8}  {'J':>14}  {'S_BH':>14}  "
              f"{'continuum':>14}  {'ratio':>10}")
        print("  " + "-" * 68)

        S_BH_cutoff: Dict[str, Dict] = {}
        for Lambda in Lambda_values:
            integ = mellin_integral(
                t_grid=t_grid, tip_term=tip_term,
                Lambda=Lambda, cutoff_fn=cutoff_fn,
            )
            S_BH_continuum = r_h ** 2 * Lambda ** 2 / 3.0
            S_BH_extracted = integ["S_BH_extracted"]
            ratio = (
                S_BH_extracted / S_BH_continuum
                if abs(S_BH_continuum) > 1e-15 else float("inf")
            )

            S_BH_cutoff[str(Lambda)] = {
                "Lambda": Lambda,
                "J_integral": integ["J_integral"],
                "S_BH_extracted": S_BH_extracted,
                "S_BH_continuum": S_BH_continuum,
                "ratio_disc_to_cont": ratio,
                "cutoff_values": integ["cutoff_values"].tolist(),
                "integrand_log_t": integ["integrand_log_t"].tolist(),
            }
            print(f"  {Lambda:>8.2f}  {integ['J_integral']:>14.4f}  "
                  f"{S_BH_extracted:>14.4f}  {S_BH_continuum:>14.4f}  "
                  f"{ratio:>10.4f}")

        S_BH_table[cutoff_name] = S_BH_cutoff

    results["S_BH_table"] = S_BH_table

    # ------------------------------------------------------------------
    # Per-cell diagnostic: at Lambda=0.5, where does the integral mass live?
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Diagnostic] Per-t-cell integrand contribution at Lambda=0.5")
    print("=" * 76)

    Lambda_diag = 0.5
    log_t = np.log(t_grid)
    # Pair-wise trapezoid contributions ((t_i+t_{i+1})/2 * (log t_{i+1} - log t_i)) ~ sum dlog
    dlogt = np.diff(log_t)
    diag_results: Dict[str, Dict] = {}

    for cutoff_name in cutoff_names:
        cutoff_fn = CUTOFFS[cutoff_name]
        x = t_grid * Lambda_diag ** 2
        f_cutoff = cutoff_fn(x)
        integrand = f_cutoff * tip_term
        # trapezoid pair-cell contribution between t_i and t_{i+1}:
        # = (integrand_i + integrand_{i+1}) / 2 * (log t_{i+1} - log t_i)
        cell_contrib = 0.5 * (integrand[:-1] + integrand[1:]) * dlogt
        S_BH_total = float(np.sum(cell_contrib)) / 2.0
        # find biggest contributor cell
        biggest_idx = int(np.argmax(np.abs(cell_contrib)))
        biggest_t_low = float(t_grid[biggest_idx])
        biggest_t_high = float(t_grid[biggest_idx + 1])
        biggest_contrib = float(cell_contrib[biggest_idx])

        # cumulative fraction
        if abs(S_BH_total) > 1e-15:
            frac_per_cell = cell_contrib / (2.0 * S_BH_total)
        else:
            frac_per_cell = np.zeros_like(cell_contrib)

        print(f"\n  Cutoff: {cutoff_name}")
        print(f"  S_BH(Lambda={Lambda_diag}) total = {S_BH_total:.4f}")
        print(f"  Biggest cell: t in [{biggest_t_low}, {biggest_t_high}], "
              f"contribution = {biggest_contrib:.4f} (S_BH contrib "
              f"{biggest_contrib/2:.4f})")
        print(f"\n  {'t_low':>8}  {'t_high':>8}  {'integrand_low':>14}  "
              f"{'integrand_high':>14}  {'cell contrib':>14}  {'frac':>8}")
        print("  " + "-" * 78)
        for i in range(len(t_grid) - 1):
            print(f"  {t_grid[i]:>8.4f}  {t_grid[i+1]:>8.4f}  "
                  f"{integrand[i]:>14.4f}  {integrand[i+1]:>14.4f}  "
                  f"{cell_contrib[i]:>14.4f}  {frac_per_cell[i]:>+8.3f}")

        diag_results[cutoff_name] = {
            "Lambda": Lambda_diag,
            "S_BH_total": S_BH_total,
            "biggest_t_low": biggest_t_low,
            "biggest_t_high": biggest_t_high,
            "biggest_cell_contribution": biggest_contrib,
            "biggest_cell_S_BH_contribution": biggest_contrib / 2.0,
            "cell_contributions": cell_contrib.tolist(),
            "fraction_per_cell": frac_per_cell.tolist(),
        }

    results["diagnostic_Lambda_0p5"] = diag_results

    # ------------------------------------------------------------------
    # Summary -- best cutoff per Lambda
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Summary] Best cutoff per Lambda + factor-2 gate check")
    print("=" * 76)

    summary_table: Dict[str, Dict] = {}
    print(f"\n  {'Lambda':>8}  {'Gaussian':>10}  {'Polynomial':>12}  "
          f"{'Sharp':>10}  {'Best':>10}")
    print("  " + "-" * 62)
    for Lambda in Lambda_values:
        ratios = {
            cn: S_BH_table[cn][str(Lambda)]["ratio_disc_to_cont"]
            for cn in cutoff_names
        }
        # "best" = closest to 1.0
        best_cutoff = min(ratios.items(), key=lambda kv: abs(np.log(kv[1])) if kv[1] > 0 else float("inf"))
        summary_table[str(Lambda)] = {
            "Lambda": Lambda,
            "ratio_gaussian": ratios["gaussian"],
            "ratio_polynomial": ratios["polynomial"],
            "ratio_sharp": ratios["sharp"],
            "best_cutoff": best_cutoff[0],
            "best_ratio": best_cutoff[1],
        }
        print(f"  {Lambda:>8.2f}  {ratios['gaussian']:>10.4f}  "
              f"{ratios['polynomial']:>12.4f}  {ratios['sharp']:>10.4f}  "
              f"{best_cutoff[0]:>10}")

    results["summary_best_cutoff_per_Lambda"] = summary_table

    # ------------------------------------------------------------------
    # Factor-2 gate check per cutoff
    # ------------------------------------------------------------------
    print(f"\n[Factor-2 gate check]")
    gate_results: Dict[str, Dict] = {}
    for cutoff_name in cutoff_names:
        ratios = [
            S_BH_table[cutoff_name][str(L)]["ratio_disc_to_cont"]
            for L in Lambda_values
        ]
        n_within = sum(1 for r in ratios if 0.5 <= r <= 2.0)
        all_within = all(0.5 <= r <= 2.0 for r in ratios)
        gate_results[cutoff_name] = {
            "ratios_all_Lambda": ratios,
            "n_within_factor_2": n_within,
            "n_total": len(Lambda_values),
            "all_within_factor_2": all_within,
        }
        print(f"  {cutoff_name:>12}: {n_within}/{len(Lambda_values)} Lambda cells "
              f"within factor 2 (all_pass = {all_within})")

    results["factor_2_gate_per_cutoff"] = gate_results

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Verdict]")
    print("=" * 76)

    # G4-5c baseline gate: Gaussian, Lambda in {0.5, 1, 2}
    baseline_ratios = [
        S_BH_table["gaussian"][str(L)]["ratio_disc_to_cont"]
        for L in [0.5, 1.0, 2.0]
    ]
    g4_5c_baseline_n_within = sum(1 for r in baseline_ratios if 0.5 <= r <= 2.0)

    # POSITIVE: ALL 5 (effectively all 6) Lambda cells within factor 2 for any cutoff
    any_cutoff_positive = any(
        gate_results[cn]["all_within_factor_2"] for cn in cutoff_names
    )

    # PARTIAL: improvement over G4-5c at Lambda=0.5, 1.0; OR ratio < factor 2 at 3/5
    best_count = max(
        gate_results[cn]["n_within_factor_2"] for cn in cutoff_names
    )
    partial_3_of_5 = best_count >= 3

    # Check improvement at Lambda=0.5 / 1.0 vs baseline (G4-5c gaussian = 29.08, 5.81)
    g4_5c_ratio_lambda_0p5 = 29.08
    g4_5c_ratio_lambda_1p0 = 5.81

    best_lambda_0p5 = min(
        abs(np.log(S_BH_table[cn][str(0.5)]["ratio_disc_to_cont"]))
        for cn in cutoff_names
        if S_BH_table[cn][str(0.5)]["ratio_disc_to_cont"] > 0
    )
    best_lambda_1p0 = min(
        abs(np.log(S_BH_table[cn][str(1.0)]["ratio_disc_to_cont"]))
        for cn in cutoff_names
        if S_BH_table[cn][str(1.0)]["ratio_disc_to_cont"] > 0
    )
    improved_at_0p5 = np.exp(best_lambda_0p5) < g4_5c_ratio_lambda_0p5
    improved_at_1p0 = np.exp(best_lambda_1p0) < g4_5c_ratio_lambda_1p0
    partial_improvement = improved_at_0p5 and improved_at_1p0

    if any_cutoff_positive:
        verdict = "POSITIVE"
        msg = (
            f"At least one cutoff achieves S_BH within factor 2 of continuum at "
            f"ALL Lambda in {Lambda_values}."
        )
    elif partial_3_of_5 or partial_improvement:
        verdict = "PARTIAL"
        msg = (
            f"Best cutoff achieves {best_count}/{len(Lambda_values)} cells within factor 2. "
            f"Improvement vs G4-5c baseline at Lambda=0.5: {improved_at_0p5}, "
            f"Lambda=1.0: {improved_at_1p0}."
        )
    else:
        verdict = "NEGATIVE"
        msg = (
            "None of the three cutoffs improves the IR cells beyond G4-5c baseline."
        )

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")
    print()
    print(f"  G4-5c baseline (Gaussian):")
    for L in [0.5, 1.0, 2.0]:
        print(f"    Lambda={L}: ratio = {S_BH_table['gaussian'][str(L)]['ratio_disc_to_cont']:.4f}")
    print(f"  Best per-Lambda cutoff:")
    for L in Lambda_values:
        bc = summary_table[str(L)]["best_cutoff"]
        br = summary_table[str(L)]["best_ratio"]
        print(f"    Lambda={L}: best={bc} ratio={br:.4f}")

    results["verdict"] = verdict
    results["verdict_message"] = msg
    results["g4_5c_baseline_n_within_3_Lambda"] = g4_5c_baseline_n_within
    results["best_count_within_factor_2"] = best_count
    results["partial_improvement_at_Lambda_0p5_and_1p0"] = partial_improvement

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 76)


if __name__ == "__main__":
    main()
