"""Sprint G4-6b -- IR boundary regularization first-move.

Per the reframed G4-6 sequencing (thread 4 Track C of 2026-05-29):
G4-6b is the sequential prerequisite to G4-6a refined. The B.2 spectral
substrate sweep found B_fit ~= 0.30 in both panels (with R = 10 fixed),
overshooting the continuum Lichnerowicz constant +1/6 ~= 0.167.

This driver probes whether the IR cutoff R = N_rho * a is the dominant
source of B inflation. Procedure:

1. Hold a = 0.05 fixed (matching B.2 baseline).
2. Vary N_rho in {100, 200, 400, 800} -> R in {5, 10, 20, 40}.
3. Compute tip(t = 1.0) at each R; at this t, A/t ~ 0.013 is small
   compared to B ~ 0.3, so tip(1.0) is approximately B.
4. Check whether B converges to +1/6 as R increases.

If R -> infinity gives B -> +1/6, the substrate-finite contribution is
analytically subtractable as a 1/R^2 correction, enabling G4-6a refined
clean A extraction.

If not, the B inflation comes from elsewhere (UV truncation,
discretization), and the B-subtraction strategy needs alternative form.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_6b_ir_boundary_first_move.json"
OUT_JSON.parent.mkdir(exist_ok=True)


B_TARGET = 1.0 / 6.0  # Continuum Lichnerowicz constant


def compute_tip_at_t(N_rho: int, a: float, N_0: int, t: float,
                     k_step: int = 12) -> float:
    """Central FD: tip(t) = dK_wedge/d_alpha - K_disk at alpha = 1."""
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
    return tip


def main() -> None:
    print("=" * 72)
    print("Sprint G4-6b -- IR boundary regularization FIRST MOVE")
    print("=" * 72)

    a = 0.05  # fixed (B.2 baseline)
    N_0 = 120
    t_panel = [1.0, 5.0, 10.0]  # large-t, B-dominated regime

    print(f"\n  Setup:")
    print(f"    a = {a} (fixed)")
    print(f"    N_0 = {N_0} (fixed azimuthal mode count)")
    print(f"    t values = {t_panel} (B-dominated regime where A/t << B)")
    print(f"    Continuum B target = +1/6 = {B_TARGET:.6f}")

    # ------------------------------------------------------------------
    # Step 1: Vary R via N_rho
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 1: Vary R = N_rho * a")
    print(f"{'=' * 72}")

    N_rho_values = [100, 200, 400, 600]  # R in {5, 10, 20, 30}
    panel_data = []

    for N_rho in N_rho_values:
        R = N_rho * a
        print(f"\n  N_rho = {N_rho}, R = {R}")
        t0 = time.time()
        panel = {"N_rho": N_rho, "a": a, "R": R, "tip_values": []}
        for t in t_panel:
            tip = compute_tip_at_t(N_rho=N_rho, a=a, N_0=N_0, t=t)
            elapsed = time.time() - t0
            panel["tip_values"].append({"t": t, "tip": tip})
            A_over_t = (tip - B_TARGET) * 0  # placeholder, real fit later
            print(f"    t={t}: tip = {tip:+.6f}, "
                  f"tip - B_target = {tip - B_TARGET:+.6f}  [{elapsed:.1f}s]")
        panel_data.append(panel)

    # ------------------------------------------------------------------
    # Step 2: B extrapolation
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 2: B extrapolation as R -> infinity")
    print(f"{'=' * 72}")

    print(f"\n  Per-R B estimates (using tip(t=10), where A/t ~ 0.001 << B):")
    print(f"    {'N_rho':>6}  {'R':>5}  {'B_est':>10}  {'B - B_target':>14}  {'rel':>8}")
    R_values = []
    B_estimates = []
    for panel in panel_data:
        B_est = next(tv["tip"] for tv in panel["tip_values"] if tv["t"] == 10.0)
        diff = B_est - B_TARGET
        rel = diff / B_TARGET * 100
        R_values.append(panel["R"])
        B_estimates.append(B_est)
        print(f"    {panel['N_rho']:>6d}  {panel['R']:>5.1f}  "
              f"{B_est:>10.6f}  {diff:>+14.6f}  {rel:>+7.2f}%")

    # Richardson-style extrapolation
    print(f"\n  Richardson extrapolation in 1/R^2:")
    if len(R_values) >= 2:
        R_arr = np.array(R_values)
        B_arr = np.array(B_estimates)
        # Fit B = B_inf + C/R^2
        X = np.vstack([np.ones_like(R_arr), 1.0 / R_arr ** 2]).T
        B_inf, C = np.linalg.lstsq(X, B_arr, rcond=None)[0]
        print(f"    Fit: B(R) = B_inf + C/R^2")
        print(f"    B_inf = {B_inf:.6f}")
        print(f"    C     = {C:.6f}")
        print(f"    B_target (continuum) = {B_TARGET:.6f}")
        print(f"    B_inf - B_target = {B_inf - B_TARGET:+.6f} "
              f"({(B_inf - B_TARGET)/B_TARGET*100:+.2f}% rel)")
        print(f"    R^2 of fit: {np.corrcoef(R_arr**(-2), B_arr)[0,1]**2:.4f}")

        # 1/R^2 trend test
        rel_errs = (np.array(B_estimates) - B_TARGET) / B_TARGET
        print(f"\n  Per-panel relative error vs target:")
        print(f"    R=5:  rel_err = {rel_errs[0]*100:+.2f}%")
        if len(rel_errs) >= 2:
            print(f"    R=10: rel_err = {rel_errs[1]*100:+.2f}%")
        if len(rel_errs) >= 3:
            print(f"    R=20: rel_err = {rel_errs[2]*100:+.2f}%")

        # If error decreases monotonically with R, IR boundary is dominant
        monotonic_decrease = all(
            abs(rel_errs[i+1]) < abs(rel_errs[i])
            for i in range(len(rel_errs) - 1)
        )
        print(f"    Monotonic decrease |error|: {monotonic_decrease}")
    else:
        B_inf = C = float("nan")
        monotonic_decrease = False

    # ------------------------------------------------------------------
    # Step 3: Verdict
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 3: Verdict")
    print(f"{'=' * 72}")

    if not np.isnan(B_inf):
        if abs(B_inf - B_TARGET) / B_TARGET < 0.1:
            verdict = "POSITIVE-IR-BOUNDARY-DOMINANT"
            print(f"\n  [{verdict}] B_inf converges to within 10% of B_target.")
            print(f"  IR boundary is the dominant source of B inflation.")
            print(f"  G4-6b can subtract 1/R^2 contribution analytically;")
            print(f"  G4-6a refined will then have clean A extraction.")
        elif monotonic_decrease:
            verdict = "POSITIVE-WITH-LARGER-R-NEEDED"
            print(f"\n  [{verdict}] Error decreases monotonically with R")
            print(f"  but B_inf doesn't converge within 10% at tested R values.")
            print(f"  Larger R panels needed to identify true B_inf;")
            print(f"  IR boundary remains a candidate dominant source.")
        else:
            verdict = "NEUTRAL-NON-MONOTONIC"
            print(f"\n  [{verdict}] Error does not decrease monotonically with R.")
            print(f"  IR boundary is NOT the dominant source of B inflation.")
            print(f"  G4-6b needs alternative B-subtraction strategy.")
    else:
        verdict = "INCONCLUSIVE"
        print(f"\n  [{verdict}] Need more data points for extrapolation.")

    # Save
    results = {
        "B_target": B_TARGET,
        "B_target_description": "Continuum Lichnerowicz constant +1/6",
        "a": a,
        "N_0": N_0,
        "panels": panel_data,
        "B_extrapolation": {
            "B_inf": float(B_inf) if not np.isnan(B_inf) else None,
            "C_coefficient": float(C) if not np.isnan(C) else None,
            "monotonic_decrease": monotonic_decrease,
            "verdict": verdict,
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
