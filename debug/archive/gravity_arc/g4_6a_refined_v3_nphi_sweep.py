"""Sprint G4-6a refined v3 -- N_phi sweep diagnostic.

Per Track a of thread 7 (debug/g4_6a_refined_simplified_extraction_memo.md):
substrate UV recovery is set by N_phi (azimuthal mode count), not by a
alone. Hold (a, R) fixed; vary N_0; track A_recovery.

Expected: monotonic convergence of A_recovery toward A_cont as N_0 grows.

Setup:
- a = 0.05 fixed (matches B.2 / G4-6b baseline)
- N_rho = 200 fixed (R = 10)
- N_0 in {60, 120, 240, 480}
- t at peak extraction window (Panel 1's intermediate-t peak)
- B_substrate = 0.163 (from G4-6b)
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_6a_refined_v3_nphi_sweep.json"
OUT_JSON.parent.mkdir(exist_ok=True)


B_SUBSTRATE = 0.163  # From G4-6b at R = 10
A_CONT = 1.0 / (24.0 * np.pi)


def extract_A_at_t(N_rho: int, a: float, N_0: int, t: float,
                   k_step: int = 12) -> dict:
    """Extract A coefficient at given t via simplified strategy.

    A_est = t * (tip(t) - B_substrate)
    """
    disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_0)
    K_disk = disk.heat_trace(t)

    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0

    wedge_plus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
    )
    wedge_minus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
    )
    K_plus = wedge_plus.heat_trace(t)
    K_minus = wedge_minus.heat_trace(t)
    dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
    tip = dK_dalpha - K_disk

    A_est = t * (tip - B_SUBSTRATE)
    return {
        "t": t,
        "tip": tip,
        "tip_minus_B": tip - B_SUBSTRATE,
        "A_est": A_est,
        "A_recovery": A_est / A_CONT,
    }


def main() -> None:
    print("=" * 72)
    print("Sprint G4-6a refined v3 -- N_phi SWEEP")
    print("=" * 72)

    a = 0.05
    N_rho = 200
    R = N_rho * a
    B_subs = B_SUBSTRATE

    print(f"\n  Fixed: a={a}, N_rho={N_rho}, R={R}, B_substrate={B_subs}")
    print(f"  Continuum A target = 1/(24*pi) = {A_CONT:.6f}")
    print(f"  Strategy: A_est = t * (tip(t) - B_substrate)")

    # Multiple t values to track per-N_0 trajectory
    # Use Panel 1's peak window (intermediate t ~ 5a^2 to 10a^2)
    t_values = [5 * a**2, 10 * a**2, 25 * a**2, 50 * a**2]
    print(f"  t values: {[f'{t:.5f}' for t in t_values]}")

    N_0_values = [60, 120, 240, 480]
    print(f"  N_0 sweep: {N_0_values}")

    sweep_data = []
    for N_0 in N_0_values:
        print(f"\n{'=' * 72}")
        print(f"N_0 = {N_0}")
        print(f"{'=' * 72}")
        t0 = time.time()
        cell_data = {"N_0": N_0, "extractions": []}
        for t in t_values:
            result = extract_A_at_t(N_rho, a, N_0, t)
            cell_data["extractions"].append(result)
            elapsed = time.time() - t0
            print(f"  t={t:.5f}: tip={result['tip']:+.4f}, "
                  f"tip-B={result['tip_minus_B']:+.4f}, "
                  f"A_est={result['A_est']:+.5f}, "
                  f"recovery={result['A_recovery']*100:+.2f}%  [{elapsed:.1f}s]")
        # Take mean across the t values
        A_means = [e["A_est"] for e in cell_data["extractions"]]
        cell_data["A_mean"] = float(np.mean(A_means))
        cell_data["A_recovery_mean"] = float(np.mean(A_means) / A_CONT)
        sweep_data.append(cell_data)
        print(f"\n  Cell summary: A_mean = {cell_data['A_mean']:+.5f} "
              f"({cell_data['A_recovery_mean']*100:+.2f}% of A_cont)")

    # ---------------------------------------------------------------------
    # Analysis
    # ---------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Trajectory analysis")
    print(f"{'=' * 72}")
    print(f"\n  {'N_0':>6}  {'A_mean':>10}  {'A_recovery':>12}")
    for cell in sweep_data:
        print(f"  {cell['N_0']:>6d}  {cell['A_mean']:>+10.5f}  "
              f"{cell['A_recovery_mean']*100:>+11.2f}%")

    # Check monotonic improvement
    recoveries = [cell["A_recovery_mean"] for cell in sweep_data]
    monotonic = all(recoveries[i+1] > recoveries[i] for i in range(len(recoveries)-1))
    print(f"\n  Monotonic improvement with N_0: {monotonic}")

    if recoveries[-1] > recoveries[0]:
        improvement = recoveries[-1] / max(abs(recoveries[0]), 0.001)
        print(f"  Improvement N_0=60 -> N_0=480: "
              f"{recoveries[0]*100:.2f}% -> {recoveries[-1]*100:.2f}% "
              f"({improvement:.2f}x ratio)")

    # Richardson extrapolation in 1/N_0
    print(f"\n  Richardson extrapolation in 1/N_0:")
    N_arr = np.array([cell["N_0"] for cell in sweep_data])
    A_arr = np.array([cell["A_mean"] for cell in sweep_data])
    # Fit A(N) = A_inf + C/N^p
    # Two-parameter fit: A_inf + C/N (assume p=1)
    X = np.vstack([np.ones_like(N_arr, dtype=float), 1.0 / N_arr]).T
    A_inf_p1, C_p1 = np.linalg.lstsq(X, A_arr, rcond=None)[0]
    print(f"    Assuming 1/N convergence:")
    print(f"      A_inf = {A_inf_p1:+.5f} ({A_inf_p1/A_CONT*100:.2f}% of A_cont)")
    print(f"      C     = {C_p1:.5f}")

    # 1/N^2 fit
    X = np.vstack([np.ones_like(N_arr, dtype=float), 1.0 / N_arr**2]).T
    A_inf_p2, C_p2 = np.linalg.lstsq(X, A_arr, rcond=None)[0]
    print(f"    Assuming 1/N^2 convergence:")
    print(f"      A_inf = {A_inf_p2:+.5f} ({A_inf_p2/A_CONT*100:.2f}% of A_cont)")
    print(f"      C     = {C_p2:.5f}")

    # Verdict
    print(f"\n{'=' * 72}")
    print("Verdict")
    print(f"{'=' * 72}")
    if monotonic and recoveries[-1] > 0.5:
        verdict = "POSITIVE-CONVERGING"
        print(f"\n  [{verdict}] A trajectory monotonically improves with N_0 and reaches >50%.")
        print(f"  N_phi-sweep IS the correct refinement axis (confirming Track a finding).")
        print(f"  G4-6a refined can close to high precision via N_phi refinement.")
    elif monotonic:
        verdict = "POSITIVE-MONOTONIC-BUT-SLOW"
        print(f"\n  [{verdict}] A trajectory improves monotonically but not yet >50%.")
        print(f"  Larger N_0 values may be needed for full convergence.")
    else:
        verdict = "NEUTRAL-NON-MONOTONIC"
        print(f"\n  [{verdict}] A trajectory not cleanly monotonic in N_0.")
        print(f"  Additional structural factor at play; investigate other axes.")

    results = {
        "B_substrate": B_subs,
        "A_cont": A_CONT,
        "a": a,
        "N_rho": N_rho,
        "R": R,
        "t_values": t_values,
        "N_0_sweep": sweep_data,
        "monotonic_improvement": monotonic,
        "Richardson": {
            "A_inf_1_over_N": float(A_inf_p1),
            "C_1_over_N": float(C_p1),
            "A_inf_1_over_N2": float(A_inf_p2),
            "C_1_over_N2": float(C_p2),
        },
        "verdict": verdict,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
