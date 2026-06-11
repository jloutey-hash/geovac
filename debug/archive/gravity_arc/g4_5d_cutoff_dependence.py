"""Sprint G4-5d -- Cutoff-function dependence sweep for tip-only S_BH.

Closes F12 falsifier of G4-5: verify that the tip contribution to S_BH
on the discrete cigar substrate scales with the second Mellin moment
phi(2) of the cutoff function f, as predicted by G8 (Paper 51 sec.10).

Method
------
1. Recompute the tip term dDelta_K/dalpha - K_disk on the G4-5a sweet-spot
   panel (R=10, a=0.05, N_rho=200, N_0=120, eps=0.1) at the 8-point
   t-grid t in {0.1, 0.2, 0.5, 1, 2, 5, 10, 20}.

2. For each of three standard cutoff functions f, integrate

     S_tip(Lambda, f) = (1/2) * integral (dt/t) f(t Lambda^2) tip(t)

   on the t-grid via trapezoidal quadrature in log(t).

   The three cutoffs:
     - f_Gauss(x)    = exp(-x)       -> phi(s) = Gamma(s),       phi(2) = 1
     - f_sharp(x)    = Theta(1 - x)  -> phi(s) = 1/s,            phi(2) = 1/2
     - f_poly(x)     = exp(-x^2)     -> phi(s) = Gamma(s/2)/2,   phi(2) = 1/2

3. At fixed Lambda, compute pairwise ratios S_tip(Gauss)/S_tip(sharp)
   and S_tip(Gauss)/S_tip(poly).  Continuum G8 prediction (rigid
   Mellin-moment scaling) is

     S_tip(f1) / S_tip(f2)   =   phi_{f1}(2)  /  phi_{f2}(2)

   i.e. both ratios should equal 1 / (1/2) = 2.

4. Tabulate S_tip(Lambda, f) for Lambda in {0.5, 1, 2} and all three f.

5. Verify the Lambda-dependence pattern:
     - sharp cuts off all weight above t = 1/Lambda^2
     - Gaussian decays as exp(-t Lambda^2)
     - polynomial decays as exp(-(t Lambda^2)^2)

   Each should give S_tip decreasing with Lambda but at different rates.

Decision gate
-------------
- POSITIVE: empirical ratios within 20% of phi(2) prediction at fixed Lambda
- PARTIAL : correct ordering but quantitative gap >20%
- NEGATIVE: cutoff dependence does NOT follow Mellin-moment scaling
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteWedgeDirac,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5d_cutoff_dependence.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# ----------------------------------------------------------------------
#  Three standard cutoff functions and their second Mellin moments
# ----------------------------------------------------------------------

def f_gauss(x: np.ndarray) -> np.ndarray:
    """Gaussian cutoff f(x) = exp(-x)."""
    return np.exp(-x)


def f_sharp(x: np.ndarray) -> np.ndarray:
    """Sharp cutoff f(x) = Theta(1 - x)."""
    return (x <= 1.0).astype(float)


def f_poly(x: np.ndarray) -> np.ndarray:
    """Polynomial-decay cutoff f(x) = exp(-x^2)."""
    return np.exp(-(x ** 2))


CUTOFFS = {
    "gauss": {
        "label": "Gaussian f(x)=exp(-x)",
        "f": f_gauss,
        "phi_s": "Gamma(s)",
        "phi_2": 1.0,
    },
    "sharp": {
        "label": "Sharp f(x)=Theta(1-x)",
        "f": f_sharp,
        "phi_s": "1/s",
        "phi_2": 0.5,
    },
    "poly": {
        "label": "Polynomial f(x)=exp(-x^2)",
        "f": f_poly,
        "phi_s": "Gamma(s/2)/2",
        "phi_2": 0.5,
    },
}


def S_tip_from_grid(t_grid: np.ndarray, tip_term: np.ndarray,
                    Lambda: float, cutoff_key: str) -> dict:
    """Integrate (1/2) * integral (dt/t) f(t Lambda^2) tip(t) by
    log-trapezoidal on the substrate t-grid.
    """
    f = CUTOFFS[cutoff_key]["f"]
    x = t_grid * Lambda ** 2
    fx = f(x)
    # integrand for integral (dt/t) g(t) on log-t grid:
    #   integral d(log t) g(t)
    integrand_log = fx * tip_term
    log_t = np.log(t_grid)
    J = float(np.trapezoid(integrand_log, log_t))
    S_tip = +J / 2
    return {
        "J": J,
        "S_tip": S_tip,
        "log_S_tip": (float(np.log(abs(S_tip))) if abs(S_tip) > 1e-15
                      else float("-inf")),
        "x_grid": x.tolist(),
        "f_x": fx.tolist(),
        "integrand_log": integrand_log.tolist(),
    }


def main() -> None:
    results: dict = {}
    print("=" * 72)
    print("Sprint G4-5d -- Cutoff-function dependence sweep")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Step 0 -- recompute tip term on the G4-5a sweet-spot panel
    # ------------------------------------------------------------------
    R = 10.0
    a = 0.05
    N_rho = 200
    N_0 = 120
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps = (alpha_plus - alpha_minus) / 2
    print(f"\n[Substrate]")
    print(f"  R = {R}, a = {a}, N_rho = {N_rho}, N_0 = {N_0}")
    print(f"  alpha_+ = {alpha_plus:.4f}, alpha_- = {alpha_minus:.4f}, eps = {eps:.4f}")

    t_grid = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    print(f"  t-grid: {t_grid.tolist()}")

    print(f"\n[Step 0] Recompute tip term on substrate ...")
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    wedge_plus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
    )
    wedge_minus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
    )

    K_disk = np.array([disk.heat_trace(t) for t in t_grid])
    K_plus = np.array([wedge_plus.heat_trace(t) for t in t_grid])
    K_minus = np.array([wedge_minus.heat_trace(t) for t in t_grid])
    dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
    tip_term = dK_dalpha - K_disk

    print(f"\n  {'t':>5}  {'K_disk':>12}  {'dK/dalpha':>12}  {'tip term':>12}")
    print("  " + "-" * 50)
    for i, t in enumerate(t_grid):
        print(f"  {t:>5.1f}  {K_disk[i]:>12.4f}  {dK_dalpha[i]:>12.4f}  "
              f"{tip_term[i]:>+12.6f}")

    results["substrate"] = {
        "R": R, "a": a, "N_rho": N_rho, "N_0": N_0,
        "alpha_plus": alpha_plus, "alpha_minus": alpha_minus, "eps": eps,
    }
    results["heat_trace_data"] = {
        "t_grid": t_grid.tolist(),
        "K_disk": K_disk.tolist(),
        "K_plus": K_plus.tolist(),
        "K_minus": K_minus.tolist(),
        "dK_dalpha": dK_dalpha.tolist(),
        "tip_term": tip_term.tolist(),
    }

    # Sanity check vs G4-5a stored data
    expected_tip_at_t1 = 0.15376729381481624
    diff = abs(tip_term[3] - expected_tip_at_t1)
    print(f"\n[Sanity] tip_term(t=1.0) = {tip_term[3]:.10f}")
    print(f"         G4-5a value      = {expected_tip_at_t1:.10f}")
    print(f"         diff             = {diff:.3e}")
    assert diff < 1e-9, "tip term diverges from G4-5a sweet-spot panel"
    print(f"         (matches G4-5a within 1e-9 -- substrate reproducible)")

    # ------------------------------------------------------------------
    # Step 1 -- S_tip(Lambda, f) for three cutoffs, three Lambda
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[Step 1] S_tip(Lambda, f) for three cutoffs:")
    print("=" * 72)

    Lambda_values = [0.5, 1.0, 2.0]

    table_rows = []
    sweep_results: dict = {}

    for Lambda in Lambda_values:
        sweep_results[str(Lambda)] = {}
        print(f"\n[Lambda = {Lambda}]")
        print(f"  {'Cutoff':<10}  {'phi(2)':>8}  {'J':>14}  "
              f"{'S_tip':>14}  {'log|S_tip|':>10}")
        print("  " + "-" * 65)
        for key, info in CUTOFFS.items():
            res = S_tip_from_grid(t_grid, tip_term, Lambda, key)
            res["phi_2"] = info["phi_2"]
            res["label"] = info["label"]
            sweep_results[str(Lambda)][key] = res
            print(f"  {key:<10}  {info['phi_2']:>8.4f}  {res['J']:>14.6e}  "
                  f"{res['S_tip']:>+14.6e}  {res['log_S_tip']:>+10.4f}")
            table_rows.append((Lambda, key, res["S_tip"]))

    results["sweep"] = sweep_results

    # ------------------------------------------------------------------
    # Step 2 -- ratios at each Lambda vs phi(2) prediction
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[Step 2] Ratios at each Lambda vs phi(2) prediction:")
    print("=" * 72)
    print()
    print("  Rigid Mellin-moment scaling (G8):")
    print("    S_tip(f1)/S_tip(f2) = phi_{f1}(2)/phi_{f2}(2)")
    print(f"    S_tip(gauss)/S_tip(sharp) predicted = 1/(1/2) = 2.00")
    print(f"    S_tip(gauss)/S_tip(poly)  predicted = 1/(1/2) = 2.00")
    print(f"    S_tip(sharp)/S_tip(poly)  predicted = (1/2)/(1/2) = 1.00")
    print()

    ratio_table: dict = {}
    print(f"  {'Lambda':>8}  {'g/sh empir':>10}  {'g/sh pred':>10}  "
          f"{'dev %':>8}    "
          f"{'g/po empir':>10}  {'g/po pred':>10}  {'dev %':>8}    "
          f"{'sh/po empir':>10}  {'sh/po pred':>10}  {'dev %':>8}")
    print("  " + "-" * 112)
    for Lambda in Lambda_values:
        S_g = sweep_results[str(Lambda)]["gauss"]["S_tip"]
        S_s = sweep_results[str(Lambda)]["sharp"]["S_tip"]
        S_p = sweep_results[str(Lambda)]["poly"]["S_tip"]

        gs_emp = S_g / S_s if abs(S_s) > 1e-15 else float("inf")
        gs_pred = 1.0 / 0.5
        gs_dev = 100.0 * (gs_emp - gs_pred) / gs_pred

        gp_emp = S_g / S_p if abs(S_p) > 1e-15 else float("inf")
        gp_pred = 1.0 / 0.5
        gp_dev = 100.0 * (gp_emp - gp_pred) / gp_pred

        sp_emp = S_s / S_p if abs(S_p) > 1e-15 else float("inf")
        sp_pred = 0.5 / 0.5
        sp_dev = 100.0 * (sp_emp - sp_pred) / sp_pred

        ratio_table[str(Lambda)] = {
            "S_gauss": S_g,
            "S_sharp": S_s,
            "S_poly": S_p,
            "gauss_over_sharp_empirical": gs_emp,
            "gauss_over_sharp_predicted": gs_pred,
            "gauss_over_sharp_dev_pct": gs_dev,
            "gauss_over_poly_empirical": gp_emp,
            "gauss_over_poly_predicted": gp_pred,
            "gauss_over_poly_dev_pct": gp_dev,
            "sharp_over_poly_empirical": sp_emp,
            "sharp_over_poly_predicted": sp_pred,
            "sharp_over_poly_dev_pct": sp_dev,
        }

        print(f"  {Lambda:>8.2f}  {gs_emp:>10.4f}  {gs_pred:>10.4f}  "
              f"{gs_dev:>+7.2f}%    "
              f"{gp_emp:>10.4f}  {gp_pred:>10.4f}  {gp_dev:>+7.2f}%    "
              f"{sp_emp:>10.4f}  {sp_pred:>10.4f}  {sp_dev:>+7.2f}%")

    results["ratio_table"] = ratio_table

    # ------------------------------------------------------------------
    # Step 3 -- Lambda-dependence pattern check for each cutoff
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[Step 3] Lambda-dependence pattern for each cutoff:")
    print("=" * 72)
    print()
    print("  Expectation:")
    print("    sharp:   S_tip decreases rapidly (hard cutoff at t=1/Lambda^2).")
    print("    gauss:   S_tip decreases exponentially with t Lambda^2.")
    print("    poly:    S_tip decreases as exp(-(tLambda^2)^2), even sharper.")
    print()

    Lambda_dep: dict = {}
    for key in CUTOFFS:
        S_vals = [sweep_results[str(L)][key]["S_tip"] for L in Lambda_values]
        Lambda_dep[key] = {
            "Lambda": Lambda_values,
            "S_tip": S_vals,
            "S_tip_ratios": [
                S_vals[i+1] / S_vals[i] if abs(S_vals[i]) > 1e-15
                else float("inf")
                for i in range(len(S_vals) - 1)
            ],
            "decreasing": all(S_vals[i+1] < S_vals[i]
                              for i in range(len(S_vals) - 1)),
        }
        print(f"  {key:<10}: S_tip = {S_vals}")
        print(f"             ratios L(i+1)/L(i) = "
              f"{[f'{r:.4f}' for r in Lambda_dep[key]['S_tip_ratios']]}")
        print(f"             decreasing with Lambda: "
              f"{Lambda_dep[key]['decreasing']}")

    results["Lambda_dependence"] = Lambda_dep

    # ------------------------------------------------------------------
    # Step 4 -- Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[Step 4] Verdict")
    print("=" * 72)

    # Aggregate deviation across the three ratio classes at all Lambda
    all_devs = []
    for Lambda in Lambda_values:
        row = ratio_table[str(Lambda)]
        all_devs.append(abs(row["gauss_over_sharp_dev_pct"]))
        all_devs.append(abs(row["gauss_over_poly_dev_pct"]))
        all_devs.append(abs(row["sharp_over_poly_dev_pct"]))

    max_dev = max(all_devs)
    mean_dev = float(np.mean(all_devs))

    # Ordering check: S_gauss > S_sharp > S_poly at every Lambda
    # (Gaussian has the largest phi(2), polynomial the smallest tail-weight
    # after the second moment, but actually phi_sharp(2) = phi_poly(2),
    # so the empirical ordering is not strictly forced by phi(2) alone.)
    correct_ordering_at_each_L = all(
        sweep_results[str(L)]["gauss"]["S_tip"]
        >= sweep_results[str(L)]["sharp"]["S_tip"]
        and sweep_results[str(L)]["gauss"]["S_tip"]
        >= sweep_results[str(L)]["poly"]["S_tip"]
        for L in Lambda_values
    )

    print(f"\n  Max deviation from phi(2) prediction across all panel cells: "
          f"{max_dev:.2f}%")
    print(f"  Mean deviation                                              : "
          f"{mean_dev:.2f}%")
    print(f"  S_gauss >= S_sharp, S_poly at every Lambda                   : "
          f"{correct_ordering_at_each_L}")

    if max_dev <= 20.0 and correct_ordering_at_each_L:
        verdict = "POSITIVE-G4-5d-CUTOFF-MELLIN-SCALING-VERIFIED"
        msg = ("Empirical S_tip(Lambda, f) ratios match G8 phi(2) prediction "
               "within 20% at all three Lambda and all three cutoff functions. "
               "F12 falsifier closed: cutoff dependence is rigidly tied to the "
               "second Mellin moment, as predicted on the continuum side.")
    elif max_dev <= 50.0 and correct_ordering_at_each_L:
        verdict = "PARTIAL-G4-5d-CUTOFF-MELLIN-SCALING-ORDERING-ONLY"
        msg = ("Ordering correct (S_gauss >= S_sharp, S_poly), but quantitative "
               "ratios deviate >20% from phi(2). Discrete-substrate finite-size "
               "corrections enter at first order; continuum limit may close.")
    else:
        verdict = "NEGATIVE-G4-5d-CUTOFF-MELLIN-SCALING"
        msg = ("Ratios do not follow Mellin-moment scaling. G8 prediction "
               "does not transport to the discrete substrate at this resolution.")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["max_deviation_pct"] = max_dev
    results["mean_deviation_pct"] = mean_dev
    results["correct_ordering"] = correct_ordering_at_each_L
    results["msg"] = msg

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
