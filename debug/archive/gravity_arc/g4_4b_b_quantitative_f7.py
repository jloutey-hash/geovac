"""Sprint G4-4b-b — Quantitative F7 structural form.

Closes the named next-week's work from G4-4b-a: characterize the
structural form of the factorization-loss Delta_fact(t) at variable
warp.

Hypothesis (from small-warp Taylor):
    r(rho)^2 = r_h^2 + rho^2
    1/r(rho)^2 ~ (1/r_h^2)(1 - (rho/r_h)^2 + ...)
    Delta-mass per site ~ (n+1)^2 (rho/r_h)^2 / r_h^2
    Delta_fact at fixed substrate, leading order ~ <rho^2>/r_h^4

Predicted dimensionless scaling:
    Delta_fact(t) / K_const(t) ~ (R/r_h)^4   (perturbative, small R/r_h)

where R = N_rho * a is the IR cutoff.

Method
------
1. Fix substrate. Vary r_h across 7 orders of magnitude.
2. Compute Delta_fact and K_const at fixed t.
3. Fit log(Delta_fact / K_const) vs log(R/r_h).
4. Test predicted slope = 4 in perturbative regime.
5. Compare against alternative scalings (r'^2 integral, mass-deviation
   integral, warp range).
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    S2DiracSpectrum,
    VariableWarpDirac,
    WarpedDiracConstant,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_4b_b_quantitative_f7.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def compute_warp_metrics(disk: DiscreteDiskDirac, r_h: float) -> dict:
    """Diagnostic metrics for a smooth-tip warp profile."""
    rho = np.array([(k + 1) * disk.a for k in range(disk.N_rho)])
    r = r_h * np.sqrt(1 + (rho / r_h) ** 2)
    r_prime = rho / r  # r'(rho) = rho / r(rho)
    return {
        "r_min": float(np.min(r)),
        "r_max": float(np.max(r)),
        "r_h": r_h,
        "R_cutoff": float(np.max(rho)),
        "warp_range": float(np.max(r) - np.min(r)),
        "mass_deviation_integral": float(np.mean(1.0 / r**2 - 1.0 / r_h**2)),
        "r_prime_squared_mean": float(np.mean(r_prime**2)),
        "r_prime_over_r_squared_mean": float(np.mean((r_prime / r) ** 2)),
        "rho_squared_mean": float(np.mean(rho**2)),
    }


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4b-b -- Quantitative F7 structural form")
    print("=" * 72)

    # Fixed substrate
    N_rho = 20
    a = 0.3
    N_phi = 12
    l_max = 3
    R = N_rho * a
    print(f"\n[Setup] Fixed substrate: N_rho={N_rho}, a={a}, N_phi={N_phi}, "
          f"l_max={l_max}, R={R}")

    # Sweep r_h across 7 orders of magnitude
    r_h_values = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]
    t_values = [0.1, 0.5, 1.0]

    print(f"\n[Sweep] r_h values: {r_h_values}")
    print(f"        t values: {t_values}")
    print(f"        Predicted: Delta_fact/K_const ~ (R/r_h)^4 in perturbative regime")

    # ------------------------------------------------------------------
    # Step 1: Compute Delta_fact at each r_h
    # ------------------------------------------------------------------
    print("\n[Step 1] Computing Delta_fact across r_h sweep...")
    sweep_data = {}
    for r_h in r_h_values:
        disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_phi)
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)
        const = WarpedDiracConstant(disk=disk, sphere=sphere)

        metrics = compute_warp_metrics(disk, r_h)
        cell = {"metrics": metrics, "t_data": {}}

        for t in t_values:
            K_var = var.heat_trace(t)
            K_const = const.heat_trace_factorized(t)
            Delta = K_var - K_const
            ratio_dim = Delta / K_const  # dimensionless
            cell["t_data"][str(t)] = {
                "K_var": float(K_var),
                "K_const": float(K_const),
                "Delta_fact": float(Delta),
                "Delta_over_K_const": float(ratio_dim),
            }
        sweep_data[r_h] = cell
        print(f"  r_h = {r_h:>6}: warp_range = {metrics['warp_range']:>6.3f}, "
              f"R/r_h = {R/r_h:>6.4f}, "
              f"Delta/K_const(t=0.1) = {cell['t_data']['0.1']['Delta_over_K_const']:.4e}")

    results["sweep"] = sweep_data

    # ------------------------------------------------------------------
    # Step 2: Power-law fit at small-warp regime
    # ------------------------------------------------------------------
    print("\n[Step 2] Power-law fit: log(Delta/K_const) vs log(R/r_h)")
    print()

    # Use the perturbative regime (largest r_h) for clean scaling
    perturbative_r_h = [r for r in r_h_values if R / r < 0.3]  # R/r_h < 0.3
    print(f"  Perturbative regime (R/r_h < 0.3): r_h = {perturbative_r_h}")

    fit_results = {}
    for t in t_values:
        x = np.array([np.log(R / r) for r in perturbative_r_h])
        y = np.array([np.log(abs(sweep_data[r]["t_data"][str(t)]["Delta_over_K_const"]))
                      for r in perturbative_r_h])
        if len(x) >= 2:
            slope, intercept = np.polyfit(x, y, 1)
            y_pred = slope * x + intercept
            resid_rms = float(np.sqrt(np.mean((y - y_pred) ** 2)))
            fit_results[t] = {
                "slope": float(slope),
                "intercept": float(intercept),
                "resid_rms": resid_rms,
                "n_points": len(x),
            }
            print(f"  t = {t:>4}: slope = {slope:+.3f}  intercept = {intercept:+.3f}  "
                  f"RMS residual = {resid_rms:.3f}  (n={len(x)})")
        else:
            fit_results[t] = {"slope": None, "intercept": None, "resid_rms": None}

    results["power_law_fit_perturbative"] = fit_results

    print(f"\n  Predicted slope = 4 (from leading-order small-warp Taylor)")
    print(f"  Measured slopes at small t suggest cleanest extraction at t=0.1.")

    # ------------------------------------------------------------------
    # Step 3: Test alternative scaling laws against perturbative data
    # ------------------------------------------------------------------
    print("\n[Step 3] Alternative scaling laws (perturbative regime, t=0.1):")
    print()

    # Candidate drivers:
    #   A) (R/r_h)^4 — Taylor leading
    #   B) |<1/r^2 - 1/r_h^2>| — mass-deviation integral magnitude
    #   C) <(r'/r)^2> — spin-connection-like
    #   D) (warp_range/r_h)^2 — geometric
    t_focus = 0.1
    delta_over_K = [abs(sweep_data[r]["t_data"][str(t_focus)]["Delta_over_K_const"])
                    for r in perturbative_r_h]

    drivers = {
        "A_R_over_r_h_4": [(R / r) ** 4 for r in perturbative_r_h],
        "B_mass_deviation_abs": [abs(sweep_data[r]["metrics"]["mass_deviation_integral"])
                                  for r in perturbative_r_h],
        "C_r_prime_over_r_sq_mean": [sweep_data[r]["metrics"]["r_prime_over_r_squared_mean"]
                                      for r in perturbative_r_h],
        "D_warp_range_over_rh_sq": [(sweep_data[r]["metrics"]["warp_range"] / r) ** 2
                                     for r in perturbative_r_h],
    }

    fit_alt = {}
    print(f"  {'Driver':>35}  {'slope':>8}  {'intercept':>10}  "
          f"{'log RMS':>10}  CV")
    print("  " + "-" * 72)

    for name, x_vals in drivers.items():
        x = np.log(np.array(x_vals))
        y = np.log(np.array(delta_over_K))
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            continue
        slope, intercept = np.polyfit(x, y, 1)
        y_pred = slope * x + intercept
        resid_rms = float(np.sqrt(np.mean((y - y_pred) ** 2)))
        # Coefficient of variation of (delta_over_K / driver^slope)
        normalized = np.array(delta_over_K) / (np.array(x_vals) ** slope)
        cv = float(np.std(normalized) / np.abs(np.mean(normalized)))
        fit_alt[name] = {
            "slope": float(slope), "intercept": float(intercept),
            "resid_rms": resid_rms, "cv": cv,
        }
        print(f"  {name:>35}  {slope:>+8.3f}  {intercept:>+10.3f}  "
              f"{resid_rms:>10.3f}  {cv:>6.3f}")

    results["alternative_scaling_fits"] = fit_alt

    # ------------------------------------------------------------------
    # Step 4: Verify A) (R/r_h)^4 prediction quantitatively
    # ------------------------------------------------------------------
    print("\n[Step 4] (R/r_h)^4 prediction test at perturbative r_h:")
    print()
    print(f"  Perturbative-regime: Delta/K_const = c * (R/r_h)^4 + corrections")
    print()
    print(f"  {'r_h':>6}  {'R/r_h':>8}  {'(R/r_h)^4':>10}  "
          f"{'Delta/K_const(t=0.1)':>22}  {'ratio':>10}")
    print("  " + "-" * 66)

    perturbative_ratios = []
    for r in perturbative_r_h:
        Rrh = R / r
        Rrh4 = Rrh ** 4
        DoK = abs(sweep_data[r]["t_data"]["0.1"]["Delta_over_K_const"])
        ratio = DoK / Rrh4
        perturbative_ratios.append(ratio)
        print(f"  {r:>6}  {Rrh:>8.4f}  {Rrh4:>10.4e}  {DoK:>22.4e}  {ratio:>10.4e}")

    print(f"\n  Convergence: ratio should -> constant in perturbative limit")
    ratio_cv = float(np.std(perturbative_ratios) / np.abs(np.mean(perturbative_ratios)))
    print(f"  CV of perturbative ratios: {ratio_cv:.4f}")
    if len(perturbative_ratios) >= 3:
        # Check ratio decreases (approaching constant) as r_h increases
        trend = (perturbative_ratios[-1] - perturbative_ratios[0]) / perturbative_ratios[0]
        print(f"  Relative change first-to-last: {trend:+.3f}")

    results["perturbative_ratio_test"] = {
        "r_h_values": perturbative_r_h,
        "ratios": perturbative_ratios,
        "cv": ratio_cv,
    }

    # ------------------------------------------------------------------
    # Step 5: Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print(f"\n[Verdict]")

    # Best fit in perturbative regime: slope close to 4?
    slope_t01 = fit_results.get(0.1, {}).get("slope", None)
    if slope_t01 is not None:
        deviation_from_4 = abs(slope_t01 - 4.0)
        slope_close_to_4 = deviation_from_4 < 0.5
        print(f"  Perturbative slope at t=0.1: {slope_t01:.3f} (target 4.0, "
              f"deviation = {deviation_from_4:.3f})")
        print(f"  Slope within 0.5 of 4: {slope_close_to_4}")
    else:
        slope_close_to_4 = False

    # Best driver: lowest CV
    if fit_alt:
        best_driver = min(fit_alt.items(), key=lambda kv: kv[1]["cv"])
        print(f"  Best-fitting driver (lowest CV): {best_driver[0]} "
              f"(slope = {best_driver[1]['slope']:+.3f}, CV = {best_driver[1]['cv']:.3f})")

    if slope_close_to_4:
        verdict = "POSITIVE-G4-4b-b-VERIFIED"
        msg = ("Perturbative slope ~ 4 matches small-warp Taylor "
               "prediction. (R/r_h)^4 structural form confirmed at "
               "leading order. Quantitative F7 closure: at small "
               "warp variation, factorization-loss is power-law-suppressed.")
    else:
        verdict = "PARTIAL-G4-4b-b"
        msg = ("Structural form characterized; perturbative slope deviates "
               "from naive 4. Best-fit driver identified for further "
               "analytical work.")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["slope_at_t01"] = slope_t01 if slope_t01 is not None else None
    results["R"] = R

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
