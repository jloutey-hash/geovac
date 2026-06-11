"""Sprint G4-4c week 2 -- alpha > 1 branch refinement.

Per G4-4c first move closure memo: spinor SC slope -1/12 extracted at
4-digit precision at alpha < 1, but only 60-80% recovery at alpha > 1.
This sprint tests whether UV refinement (larger N_0) closes the
asymmetry.

Method
------
Sweep N_0 in {120, 240, 480} at fixed R = 10, a = 0.05. For each N_0:
  - alpha < 1: N_phi = alpha * N_0
  - alpha > 1: N_phi = alpha * N_0
Compute slope Delta_K^Dirac / (1/alpha - alpha) at fixed t = 1.0.
Check convergence to -1/12.
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

OUT_JSON = Path(__file__).parent / "data" / "g4_4c_week2_alpha_gt_1_refinement.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4c week 2 -- alpha > 1 branch refinement")
    print("=" * 72)

    R = 10.0
    a = 0.05
    N_rho = int(R / a)
    print(f"\n[Setup] R = {R}, a = {a}, N_rho = {N_rho}")

    N_0_values = [120, 240, 480]
    alpha_pairs = [
        ("1/3", 1.0/3.0),
        ("1/2", 0.5),
        ("2/3", 2.0/3.0),
        ("3/2", 1.5),
        ("2",   2.0),
        ("3",   3.0),
    ]
    t_focus = 1.0
    print(f"  N_0 sweep: {N_0_values}")
    print(f"  Alpha sweep: {[n for n, _ in alpha_pairs]}")
    print(f"  Target slope: -1/12 = {-1/12:.6f}")
    print(f"  t = {t_focus}")

    results["setup"] = {
        "R": R, "a": a, "N_rho": N_rho,
        "N_0_values": N_0_values, "t_focus": t_focus,
    }

    sweep_data = {}
    for N_0 in N_0_values:
        print(f"\n[N_0 = {N_0}] Computing slopes...")
        disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
        K_disk_at_t = disk.heat_trace(t_focus)

        slopes = {}
        for name, alpha in alpha_pairs:
            N_phi = int(round(alpha * N_0))
            wedge = DiscreteWedgeDirac(
                N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha,
            )
            K_wedge = wedge.heat_trace(t_focus)
            delta = K_wedge - alpha * K_disk_at_t
            sc = 1.0 / alpha - alpha
            slope = delta / sc
            slopes[name] = {
                "alpha": alpha,
                "N_phi": N_phi,
                "delta": float(delta),
                "1/a-a": float(sc),
                "slope": float(slope),
                "deviation_from_minus_1_over_12": float(slope - (-1/12)),
            }
            print(f"  alpha = {name:>4}  N_phi = {N_phi:>4}  "
                  f"Delta = {delta:+10.4f}  slope = {slope:+.6f}  "
                  f"deviation = {slope - (-1/12):+.6f}")

        sweep_data[N_0] = slopes

    results["sweep"] = {str(k): v for k, v in sweep_data.items()}

    # ------------------------------------------------------------------
    # Comparison: alpha > 1 branch vs alpha < 1 branch
    # ------------------------------------------------------------------
    print("\n[Convergence] alpha > 1 branch vs alpha < 1 branch:")
    print()
    print(f"  Target slope: -1/12 = {-1/12:+.6f}")
    print()
    print(f"  Mean recovery (slope / target) at each N_0:")
    print(f"  {'N_0':>5}  {'< 1 mean':>12}  {'> 1 mean':>12}  {'< 1/-1/12':>10}  {'> 1/-1/12':>10}")
    print("  " + "-" * 64)

    for N_0 in N_0_values:
        slopes_lt1 = [sweep_data[N_0][n]["slope"] for n in ["1/3", "1/2", "2/3"]]
        slopes_gt1 = [sweep_data[N_0][n]["slope"] for n in ["3/2", "2", "3"]]
        mean_lt1 = np.mean(slopes_lt1)
        mean_gt1 = np.mean(slopes_gt1)
        ratio_lt1 = mean_lt1 / (-1/12)
        ratio_gt1 = mean_gt1 / (-1/12)
        print(f"  {N_0:>5}  {mean_lt1:>+12.6f}  {mean_gt1:>+12.6f}  "
              f"{ratio_lt1:>10.4f}  {ratio_gt1:>10.4f}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    # Best alpha > 1 recovery
    best_N_0 = N_0_values[-1]
    slopes_gt1_best = [sweep_data[best_N_0][n]["slope"] for n in ["3/2", "2", "3"]]
    mean_gt1_best = np.mean(slopes_gt1_best)
    recovery_gt1 = mean_gt1_best / (-1/12)
    improved = recovery_gt1 > 0.85  # > 85% recovery at fine substrate

    print(f"\n[Verdict]")
    print(f"  At N_0 = {best_N_0} (finest), alpha > 1 mean slope = "
          f"{mean_gt1_best:+.6f}")
    print(f"  Recovery = {recovery_gt1:.4f} ({recovery_gt1*100:.1f}% of -1/12)")
    print(f"  Improvement over G4-4c first move (N_0=120 was ~ 0.65): {improved}")

    if recovery_gt1 > 0.95:
        verdict = "POSITIVE-G4-4c-WEEK2-CONVERGED"
        msg = (f"alpha > 1 branch converges to -1/12 at fine substrate "
               f"(N_0 = {best_N_0}). UV asymmetry is structurally a finite-N_0 "
               f"effect that closes with refinement.")
    elif improved:
        verdict = "POSITIVE-G4-4c-WEEK2-IMPROVED"
        msg = (f"alpha > 1 branch improves substantially (>85% recovery at "
               f"N_0={best_N_0}). Full convergence needs even larger N_0.")
    else:
        verdict = "PARTIAL-G4-4c-WEEK2"
        msg = "alpha > 1 branch shows weak UV-refinement response; structural."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["mean_slope_gt1_finest"] = float(mean_gt1_best)
    results["recovery_gt1_finest"] = float(recovery_gt1)

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
