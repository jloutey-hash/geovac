"""Sprint G4-4b-c — F5 asymptotic-free recovery saturation.

Per the G4-4b scoping memo §5 F5: at the asymptotic-free regime
(rho >> r_h), the variable warp Dirac approaches the form where
the S^2 mass (n+1)^2/r(rho)^2 -> 0. The per-radial-site contribution
saturates as

    K_S^2_per_site(t, rho -> infty) -> sum_{n=0..l_max} 8(n+1)
                                      = 4 (l_max+1)(l_max+2)

In the strict r_h -> 0 limit, ALL radial sites become asymptotic-free,
and the full heat trace saturates as

    K_var(t; r_h -> 0) / K_disk(t) -> 4 (l_max+1)(l_max+2)

(independent of t and the substrate beyond the disk part).

For l_max = 3: saturation value = 4 * 4 * 5 = 80.

Method
------
Fix substrate (N_rho, a, N_phi, l_max). Sweep r_h to small values.
Compute K_var(t) and K_disk(t). Verify ratio approaches 80 as r_h -> 0.

Cross-check: at finite l_max, saturation depends on l_max:
    saturation(l_max) = 4 (l_max+1)(l_max+2)
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    S2DiracSpectrum,
    VariableWarpDirac,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_4b_c_asymptotic_free.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4b-c -- F5 asymptotic-free recovery saturation")
    print("=" * 72)

    # Fixed substrate
    N_rho = 20
    a = 0.3
    N_phi = 12
    l_max = 3
    R = N_rho * a
    print(f"\n[Setup] Fixed substrate: N_rho={N_rho}, a={a}, N_phi={N_phi}, "
          f"l_max={l_max}, R={R}")

    # Saturation prediction
    saturation_predicted = 4 * (l_max + 1) * (l_max + 2)
    print(f"\n[Prediction] At r_h -> 0, K_var/K_disk -> 4(l_max+1)(l_max+2) "
          f"= {saturation_predicted}")

    # Disk-Dirac heat trace (r_h-independent)
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_phi)
    t_values = [0.1, 0.5, 1.0]
    K_disk_at_t = {t: disk.heat_trace(t) for t in t_values}
    print(f"\n[Reference] K_disk(t):")
    for t in t_values:
        print(f"  K_disk({t}) = {K_disk_at_t[t]:.4f}")
    results["K_disk"] = K_disk_at_t

    # ------------------------------------------------------------------
    # Step 1: r_h sweep, compute K_var/K_disk
    # ------------------------------------------------------------------
    print("\n[Step 1] r_h sweep -- K_var(t) / K_disk(t):")
    print()

    r_h_values = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    print(f"  r_h sweep: {r_h_values}")
    print(f"  R/r_h: {[f'{R/r:.1f}' for r in r_h_values]}")
    print()

    header = f"  {'r_h':>6}  {'R/r_h':>8}  " + "  ".join(
        [f"K_var/K_disk(t={t})".rjust(20) for t in t_values]
    )
    print(header)
    print("  " + "-" * (16 + 22 * len(t_values)))

    sweep_data = {}
    for r_h in r_h_values:
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)

        ratios = {}
        for t in t_values:
            K_var = var.heat_trace(t)
            ratio = K_var / K_disk_at_t[t]
            ratios[t] = ratio

        sweep_data[r_h] = {"R_over_r_h": R / r_h, "ratios": ratios}
        row = f"  {r_h:>6.2f}  {R/r_h:>8.1f}  " + "  ".join(
            [f"{ratios[t]:>20.4f}" for t in t_values]
        )
        print(row)

    results["sweep"] = {str(r): v for r, v in sweep_data.items()}

    # ------------------------------------------------------------------
    # Step 2: Check t-independence of saturation
    # ------------------------------------------------------------------
    print(f"\n[Step 2] t-independence of saturation (at smallest r_h):")
    smallest_r_h = r_h_values[-1]
    smallest_ratios = sweep_data[smallest_r_h]["ratios"]
    print(f"  At r_h = {smallest_r_h}:")
    for t in t_values:
        rel_err_from_pred = (smallest_ratios[t] - saturation_predicted) / saturation_predicted
        print(f"    t = {t}: ratio = {smallest_ratios[t]:.4f}  "
              f"vs prediction {saturation_predicted}  "
              f"(rel_err = {rel_err_from_pred:+.4f})")

    # ------------------------------------------------------------------
    # Step 3: Convergence rate as r_h -> 0
    # ------------------------------------------------------------------
    print(f"\n[Step 3] Convergence to saturation as r_h -> 0:")
    print()
    t_focus = 0.5
    print(f"  At t = {t_focus}:")
    print(f"  {'r_h':>6}  {'R/r_h':>8}  {'ratio':>12}  {'saturation':>12}  {'deficit':>12}")
    print("  " + "-" * 60)

    deficits = []
    for r_h in r_h_values:
        ratio = sweep_data[r_h]["ratios"][t_focus]
        deficit = (saturation_predicted - ratio) / saturation_predicted
        deficits.append((r_h, R / r_h, ratio, deficit))
        print(f"  {r_h:>6.3f}  {R/r_h:>8.1f}  {ratio:>12.4f}  "
              f"{saturation_predicted:>12}  {deficit:>+12.4f}")

    # Fit log(deficit) vs log(r_h) to extract power-law convergence
    valid = [(r, rrr, r_a, d) for r, rrr, r_a, d in deficits if d > 0]
    if len(valid) >= 3:
        r_h_arr = np.array([v[0] for v in valid])
        deficits_arr = np.array([v[3] for v in valid])
        # Take last few points (smallest deficit, deepest convergence)
        x = np.log(r_h_arr)
        y = np.log(deficits_arr)
        slope, intercept = np.polyfit(x, y, 1)
        print(f"\n  Convergence fit: log(deficit) ~ {slope:+.3f} * log(r_h) + {intercept:+.3f}")
        print(f"  Power-law deficit vs r_h: deficit ~ r_h^{slope:+.3f}")

    results["convergence_fit"] = {"slope": float(slope), "intercept": float(intercept)} if len(valid) >= 3 else None

    # ------------------------------------------------------------------
    # Step 4: l_max scaling cross-check
    # ------------------------------------------------------------------
    print(f"\n[Step 4] l_max scaling cross-check (at r_h = {r_h_values[-2]}):")
    print(f"  Saturation prediction: 4(l_max+1)(l_max+2)")
    print()
    print(f"  {'l_max':>7}  {'prediction':>12}  {'measured':>12}  {'rel_err':>10}")
    print("  " + "-" * 50)

    l_max_test_values = [1, 2, 3, 4]
    r_h_test = r_h_values[-2]
    t_test = t_values[1]
    l_max_data = {}
    for l_m in l_max_test_values:
        sphere = S2DiracSpectrum(l_max=l_m, r_h=r_h_test)
        var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h_test)
        K_var = var.heat_trace(t_test)
        ratio = K_var / K_disk_at_t[t_test]
        pred = 4 * (l_m + 1) * (l_m + 2)
        rel_err = (ratio - pred) / pred
        l_max_data[l_m] = {"ratio": float(ratio), "prediction": pred,
                            "rel_err": float(rel_err)}
        print(f"  {l_m:>7}  {pred:>12}  {ratio:>12.4f}  {rel_err:>+10.4f}")

    results["l_max_check"] = {"r_h": r_h_test, "t": t_test, "data": l_max_data}

    # ------------------------------------------------------------------
    # Step 5: Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    # Convergence verdict: smallest r_h ratio within 5% of prediction
    smallest_rel_err = abs((smallest_ratios[t_values[1]] - saturation_predicted)
                          / saturation_predicted)
    F5_close = smallest_rel_err < 0.05

    # l_max scaling: all l_max values close to prediction
    l_max_max_err = max(abs(d["rel_err"]) for d in l_max_data.values())
    l_max_close = l_max_max_err < 0.10

    # t-independence: spread in ratio across t small at saturation
    t_ratios = [smallest_ratios[t] for t in t_values]
    t_spread = (max(t_ratios) - min(t_ratios)) / np.mean(t_ratios)
    t_independent = t_spread < 0.05

    print(f"\n[Verdict]")
    print(f"  Saturation within 5% of prediction at smallest r_h: {F5_close} "
          f"(rel_err = {smallest_rel_err:+.4f})")
    print(f"  l_max scaling matches 4(l_max+1)(l_max+2): {l_max_close} "
          f"(max rel_err = {l_max_max_err:+.4f})")
    print(f"  t-independence at saturation: {t_independent} "
          f"(spread = {t_spread:+.4f})")

    if F5_close and l_max_close and t_independent:
        verdict = "POSITIVE-G4-4b-c-VERIFIED"
        msg = ("F5 asymptotic-free saturation verified: K_var/K_disk -> "
               "4(l_max+1)(l_max+2) as r_h -> 0, t-independent, "
               "matches l_max scaling. Sprint-scale closure.")
    elif F5_close:
        verdict = "POSITIVE-G4-4b-c-PARTIAL"
        msg = "Saturation at small r_h verified; l_max or t cross-checks need work."
    else:
        verdict = "PARTIAL-G4-4b-c"
        msg = "Saturation trend visible; not yet within 5% at sprint-scale r_h."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["saturation_predicted"] = saturation_predicted

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
