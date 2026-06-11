"""Sprint G4-4b-d first move — Level 1.5 spin-connection scalar correction.

Opens the G4-4b-d sub-sprint (Level 2 operator-level Dirac with spin
connection). The first move adds the leading SCALAR contribution
(r'/r)^2 to D^2 via the include_spin_connection flag on
VariableWarpDirac.H_block, and characterizes its magnitude vs Level 1.

Continuum reduction (Camporesi 1996):
    D = D_disk + gamma^5 D_S^2/r + (r'/r) gamma^rho
    D^2 = D_disk^2 + D_S^2^2/r^2 + (r'/r)^2 + {gamma^rho cross terms}

The (r'/r)^2 piece is a universal SCALAR shift on the radial diagonal
(n-independent), distinct from the (n+1)^2/r^2 mass term. Adding it
yields a "Level 1.5" approximation between Level 1 (G4-4b-a) and full
Level 2 (deferred multi-week sprint).

Falsifiers
----------
F6 EXTENSION: at constant warp r(rho) = r_h identically, (r'/r) = 0,
   so Level 1.5 = Level 1 = G4-4a constant warp BIT-EXACT.
   LOAD-BEARING for the spin-connection construction correctness.

QUANT: at variable warp, Level 1.5 gives a smaller K_var than Level 1
   (spin connection raises eigenvalues, lowers heat trace). Quantify
   the magnitude.
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

OUT_JSON = Path(__file__).parent / "data" / "g4_4b_d_spin_connection_scalar.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4b-d first move -- Level 1.5 spin-connection scalar")
    print("=" * 72)

    # ------------------------------------------------------------------
    # F6 EXTENSION: bit-exact reduction at constant warp
    # ------------------------------------------------------------------
    print("\n[F6 EXTENSION] Constant warp: Level 1.5 = Level 1 = constant warp")
    print()

    panels = [
        ("small",  10, 0.5, 6,  1.5, 2),
        ("medium", 20, 0.3, 12, 2.0, 3),
    ]
    t_panel = [0.1, 0.5, 1.0]
    F6_ext = {}
    for name, Nr, a, Np, r_h, l_max in panels:
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        var_const = VariableWarpDirac.constant(
            disk=disk, sphere=sphere, r_h=r_h,
        )
        const = WarpedDiracConstant(disk=disk, sphere=sphere)
        cell = {}
        all_passed = True
        for t in t_panel:
            K_lvl15_const = var_const.heat_trace(t, include_spin_connection=True)
            K_lvl1_const = var_const.heat_trace(t, include_spin_connection=False)
            K_g4_4a = const.heat_trace_factorized(t)

            rel_err_lvl15 = abs(K_lvl15_const - K_g4_4a) / abs(K_g4_4a)
            rel_err_lvl1 = abs(K_lvl1_const - K_g4_4a) / abs(K_g4_4a)
            passed = bool(rel_err_lvl15 < 1e-10 and rel_err_lvl1 < 1e-10)
            all_passed = all_passed and passed
            cell[str(t)] = {
                "K_lvl15_const": K_lvl15_const,
                "K_lvl1_const": K_lvl1_const,
                "K_g4_4a": K_g4_4a,
                "rel_err_lvl15": rel_err_lvl15,
                "rel_err_lvl1": rel_err_lvl1,
                "passed": passed,
            }
        F6_ext[name] = cell
        print(f"  Panel '{name}' (Nr={Nr}, Np={Np}, l_max={l_max}, r_h={r_h}):")
        for t in t_panel:
            entry = cell[str(t)]
            mark = "PASS" if entry["passed"] else "FAIL"
            print(f"    t={t:>4}: rel_err_lvl1.5 = {entry['rel_err_lvl15']:.2e}, "
                  f"rel_err_lvl1 = {entry['rel_err_lvl1']:.2e} [{mark}]")
        print(f"  Panel '{name}' all_passed: {all_passed}")
    results["F6_extension_constant_warp"] = F6_ext

    # ------------------------------------------------------------------
    # QUANT: Level 1.5 vs Level 1 at variable warp
    # ------------------------------------------------------------------
    print("\n[QUANT] Level 1.5 vs Level 1 at variable warp (smooth-tip):")
    quant_panels = [
        ("rh_5",  20, 0.3, 12, 5.0, 3),
        ("rh_2",  20, 0.3, 12, 2.0, 3),
        ("rh_1",  20, 0.3, 12, 1.0, 3),
        ("rh_05", 20, 0.3, 12, 0.5, 3),
    ]
    quant_results = {}
    print(f"\n  {'Panel':>8}  {'r_h':>5}  {'(r prime/r)^2 mean':>18}  "
          f"{'t':>4}  {'K_lvl15':>10}  {'K_lvl1':>10}  "
          f"{'K_var-K_lvl15':>14}  {'K_lvl15/K_lvl1':>14}")
    print("  " + "-" * 110)
    for name, Nr, a, Np, r_h, l_max in quant_panels:
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        var = VariableWarpDirac.smooth_tip(
            disk=disk, sphere=sphere, r_h=r_h,
        )
        spin_conn_sq_mean = float(np.mean(var.warp_derivative_over_warp() ** 2))
        cell = {"spin_conn_sq_mean": spin_conn_sq_mean}
        for t in t_panel:
            K_lvl15 = var.heat_trace(t, include_spin_connection=True)
            K_lvl1 = var.heat_trace(t, include_spin_connection=False)
            K_lvl1_minus_K_lvl15 = K_lvl1 - K_lvl15
            ratio = K_lvl15 / K_lvl1
            cell[str(t)] = {
                "K_lvl15": K_lvl15,
                "K_lvl1": K_lvl1,
                "Delta_lvl1_minus_lvl15": K_lvl1_minus_K_lvl15,
                "ratio_lvl15_over_lvl1": ratio,
            }
            print(f"  {name:>8}  {r_h:>5}  {spin_conn_sq_mean:>18.4f}  "
                  f"{t:>4}  {K_lvl15:>10.2f}  {K_lvl1:>10.2f}  "
                  f"{K_lvl1_minus_K_lvl15:>14.4e}  {ratio:>14.4f}")
        quant_results[name] = cell

    results["quant_variable_warp"] = quant_results

    # ------------------------------------------------------------------
    # Trend: spin-connection correction scaling
    # ------------------------------------------------------------------
    print("\n[Trend] Spin-connection correction (K_lvl1 - K_lvl15) vs (r'/r)^2:")
    print()
    t_focus = 0.5
    print(f"  At t = {t_focus}:")
    print(f"  {'Panel':>8}  {'r_h':>5}  {'(r prime/r)^2 mean':>18}  "
          f"{'Delta_corr':>14}  {'Delta_corr / K_lvl1':>20}")
    print("  " + "-" * 80)
    spin_data = []
    for name, _, _, _, r_h, _ in quant_panels:
        entry = quant_results[name][str(t_focus)]
        rps_mean = quant_results[name]["spin_conn_sq_mean"]
        delta_corr = entry["Delta_lvl1_minus_lvl15"]
        normalized = delta_corr / entry["K_lvl1"]
        spin_data.append((r_h, rps_mean, delta_corr, normalized))
        print(f"  {name:>8}  {r_h:>5}  {rps_mean:>18.4f}  "
              f"{delta_corr:>14.4e}  {normalized:>20.4e}")

    # Fit log(Delta_corr) vs log((r'/r)^2 mean)
    if len(spin_data) >= 3:
        x = np.log(np.array([d[1] for d in spin_data]))
        y = np.log(np.array([d[2] for d in spin_data]))
        slope, intercept = np.polyfit(x, y, 1)
        print(f"\n  Linear fit log(Delta_corr) ~ {slope:+.3f} * log((r'/r)^2) + {intercept:+.3f}")
        print(f"  Power-law: Delta_corr ~ ((r'/r)^2)^{slope:.3f}")
        results["scaling_fit"] = {
            "slope": float(slope),
            "intercept": float(intercept),
        }

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    F6_all = all(
        F6_ext[name][str(t)]["passed"]
        for name, _, _, _, _, _ in panels
        for t in t_panel
    )

    # Sign check: K_lvl15 < K_lvl1 at variable warp (correction is positive)
    sign_correct = all(
        quant_results[name][str(t)]["K_lvl15"] < quant_results[name][str(t)]["K_lvl1"]
        for name, _, _, _, _, _ in quant_panels
        for t in t_panel
    )

    print(f"\n[Verdict]")
    print(f"  F6 extension bit-exact at constant warp:           {F6_all}")
    print(f"  K_lvl15 < K_lvl1 at variable warp (sign correct):  {sign_correct}")

    if F6_all and sign_correct:
        verdict = "POSITIVE-G4-4b-d-FIRST-MOVE-VERIFIED"
        msg = ("F6 extension bit-exact at constant warp; Level 1.5 scalar "
               "spin-connection correction is operational. Quantitative "
               "magnitude characterized. Full Level 2 (with gamma^rho "
               "mixing cross terms) deferred to G4-4b-d multi-week.")
    elif F6_all:
        verdict = "POSITIVE-G4-4b-d-FIRST-MOVE-PARTIAL"
        msg = "F6 reduction OK; correction sign at variable warp needs check."
    else:
        verdict = "NEGATIVE-G4-4b-d-FIRST-MOVE"
        msg = "F6 reduction failed at constant warp; structural issue."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["F6_all"] = F6_all
    results["sign_correct"] = sign_correct

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
