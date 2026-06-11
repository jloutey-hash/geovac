"""Sprint G4-6a -- Multi-substrate UV foundation, first-move panel sweep.

Per the v3.20.0 continuation prompt and Task 11 mode-counting diagnostic:
- Per-t UV target: tip(t) = 1/(24*pi*t) + const + O(t)
- A_cont = 1/(24*pi) ~ 0.013263
- Spectral discretization at a=0.05 recovers ~23% of A_cont

This driver sweeps two substrate refinements (a=0.05, N_rho=200) and
(a=0.025, N_rho=400), extracts the leading 1/t coefficient A from a
linear fit over small-t, and performs a two-cell Richardson extrapolation
toward A_cont.

The decision gate:
- If A_2 > A_1 with clean Richardson trajectory toward A_cont -> positive
  foundation for G4-6a multi-month commitment.
- Otherwise -> identify structural obstruction; reassess.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

OUT_JSON = Path(__file__).parent / "data" / "g4_6a_multi_substrate_uv_first_move.json"
OUT_JSON.parent.mkdir(exist_ok=True)


A_CONT = 1.0 / (24.0 * np.pi)


def compute_tip_term_at_t(N_rho: int, a: float, N_0: int, t: float,
                          k_step: int = 12) -> dict:
    """Compute tip term at a single t value via central FD.

    tip(t) = dK_wedge/dalpha|_{alpha=1} - K_disk(t)

    Uses alpha_pm = (N_0 +/- k_step)/N_0 so that N_phi_pm = N_0 +/- k_step
    is integer; spinor anti-periodic BC.
    """
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    K_disk = disk.heat_trace(t)

    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0

    wedge_plus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
    )
    wedge_minus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
    )
    K_plus = wedge_plus.heat_trace(t)
    K_minus = wedge_minus.heat_trace(t)
    dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
    tip_term = dK_dalpha - K_disk

    return {
        "t": t,
        "K_disk": K_disk,
        "K_plus": K_plus,
        "K_minus": K_minus,
        "dK_dalpha": dK_dalpha,
        "tip_term": tip_term,
    }


def fit_A_over_t(t_values, tip_values):
    """Linear fit tip(t) = A/t + B over small-t panel."""
    t_arr = np.array(t_values)
    tip_arr = np.array(tip_values)
    # Design matrix for tip = A * (1/t) + B
    X = np.vstack([1.0 / t_arr, np.ones_like(t_arr)]).T
    A, B = np.linalg.lstsq(X, tip_arr, rcond=None)[0]
    # Residuals
    pred = A / t_arr + B
    rss = float(np.sum((tip_arr - pred) ** 2))
    tss = float(np.sum((tip_arr - tip_arr.mean()) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else float("nan")
    return float(A), float(B), r2


def run_panel(a: float, N_rho: int, N_0: int, label: str) -> dict:
    """Run one substrate panel: compute tip(t) at multiple small-t values
    and extract leading A coefficient.
    """
    print(f"\n{'=' * 72}")
    print(f"Panel {label}: a={a}, N_rho={N_rho}, N_0={N_0}, R={N_rho * a}")
    print(f"{'=' * 72}")

    t_values = [a ** 2, 2 * a ** 2, 5 * a ** 2, 10 * a ** 2, 50 * a ** 2]
    print(f"  t values: {[f'{t:.5f}' for t in t_values]}")
    print(f"  1/(24*pi*t) targets:")
    for t in t_values:
        target = 1.0 / (24.0 * np.pi * t)
        print(f"    t={t:.5f}: 1/(24*pi*t) = {target:.4f}")

    print(f"\n  Computing tip term at each t (k_step = 12 -> eps = 0.1):")
    t0 = time.time()
    raw_data = []
    for t in t_values:
        result = compute_tip_term_at_t(N_rho=N_rho, a=a, N_0=N_0, t=t, k_step=12)
        raw_data.append(result)
        target = 1.0 / (24.0 * np.pi * t)
        recovery = result["tip_term"] / target if target > 0 else float("nan")
        elapsed = time.time() - t0
        print(f"    t={t:.5f}: tip = {result['tip_term']:+.6f}, "
              f"target = {target:.4f}, recovery = {recovery:.4f} "
              f"({recovery * 100:.2f}%)  [{elapsed:.1f}s]")

    # Fit A/t + B
    tip_values = [r["tip_term"] for r in raw_data]
    A_fit, B_fit, r2 = fit_A_over_t(t_values, tip_values)
    A_recovery = A_fit / A_CONT

    print(f"\n  Linear fit tip(t) = A/t + B:")
    print(f"    A_fit  = {A_fit:+.6f}")
    print(f"    B_fit  = {B_fit:+.6f}")
    print(f"    A_cont = {A_CONT:.6f}")
    print(f"    A_fit / A_cont = {A_recovery:.4f} ({A_recovery * 100:.2f}%)")
    print(f"    R^2 = {r2:.6f}")
    print(f"  Total compute time: {time.time() - t0:.1f}s")

    return {
        "label": label,
        "a": a,
        "N_rho": N_rho,
        "N_0": N_0,
        "R": N_rho * a,
        "raw_data": raw_data,
        "A_fit": A_fit,
        "B_fit": B_fit,
        "A_recovery_vs_cont": A_recovery,
        "r2": r2,
    }


def main() -> None:
    results = {}
    results["A_cont"] = float(A_CONT)
    results["A_cont_description"] = "1/(24*pi), continuum UV target"

    print("=" * 72)
    print("Sprint G4-6a -- Multi-substrate UV foundation FIRST MOVE")
    print("=" * 72)
    print(f"\n  A_cont = 1/(24*pi) = {A_CONT:.6f}")

    # Panel 1: a = 0.05, N_rho = 200
    panel_1 = run_panel(a=0.05, N_rho=200, N_0=120, label="P1_baseline")

    # Panel 2: a = 0.025, N_rho = 400 (refinement)
    panel_2 = run_panel(a=0.025, N_rho=400, N_0=240, label="P2_refined")

    results["panels"] = [panel_1, panel_2]

    # ---------------------------------------------------------------------
    # Richardson extrapolation
    # ---------------------------------------------------------------------
    A1 = panel_1["A_fit"]
    A2 = panel_2["A_fit"]
    a1 = panel_1["a"]
    a2 = panel_2["a"]
    ratio = a1 / a2  # = 2 by construction

    print(f"\n{'=' * 72}")
    print("Richardson extrapolation: two-cell first-move")
    print(f"{'=' * 72}")
    print(f"\n  Panel 1 (a={a1}):  A = {A1:+.6f} ({A1/A_CONT*100:.2f}% of A_cont)")
    print(f"  Panel 2 (a={a2}):  A = {A2:+.6f} ({A2/A_CONT*100:.2f}% of A_cont)")

    # Assume A(a) = A_cont + C * a^p; from two cells with a1 = 2*a2:
    #   A1 - A_cont = C * a1^p
    #   A2 - A_cont = C * a2^p
    # Ratio (A1 - A_cont)/(A2 - A_cont) = (a1/a2)^p = 2^p
    # 2^p = (A1 - A_cont)/(A2 - A_cont)
    diff_1 = A1 - A_CONT
    diff_2 = A2 - A_CONT

    if diff_1 != 0 and diff_2 != 0 and diff_1 * diff_2 > 0:
        # Both same sign; can extract p
        p_implied = np.log2(diff_1 / diff_2) if diff_1 / diff_2 > 0 else float("nan")
        # Richardson at a -> 0: A_inf = (2^p * A_2 - A_1) / (2^p - 1)
        # For an assumed p=1 (linear convergence):
        A_richardson_p1 = (2.0 * A2 - A1) / 1.0  # = 2*A2 - A1
        # For implied p:
        if not np.isnan(p_implied):
            try:
                two_p = 2.0 ** p_implied
                A_richardson_implied = (two_p * A2 - A1) / (two_p - 1.0)
            except Exception:
                A_richardson_implied = float("nan")
        else:
            A_richardson_implied = float("nan")
    else:
        p_implied = float("nan")
        A_richardson_p1 = (2.0 * A2 - A1)
        A_richardson_implied = float("nan")

    print(f"\n  Implied convergence exponent p: {p_implied:.4f}")
    print(f"  (A_1 - A_cont) = {diff_1:+.6f}")
    print(f"  (A_2 - A_cont) = {diff_2:+.6f}")
    print(f"  Ratio (a1/a2)^p = {ratio**p_implied if not np.isnan(p_implied) else float('nan'):.4f}")

    print(f"\n  Richardson estimates of A(a -> 0):")
    print(f"    Assuming p=1: A_inf = 2*A_2 - A_1 = {A_richardson_p1:+.6f} "
          f"({A_richardson_p1/A_CONT*100:.2f}% of A_cont)")
    if not np.isnan(A_richardson_implied):
        print(f"    Using implied p={p_implied:.3f}: A_inf = {A_richardson_implied:+.6f} "
              f"({A_richardson_implied/A_CONT*100:.2f}% of A_cont)")

    # Decision gate
    print(f"\n{'=' * 72}")
    print("DECISION GATE")
    print(f"{'=' * 72}")
    monotonic_toward_cont = (diff_1 < 0 and diff_2 < 0 and diff_2 > diff_1) or \
                            (diff_1 > 0 and diff_2 > 0 and diff_2 < diff_1)
    if monotonic_toward_cont:
        print("  [POSITIVE] A_2 is closer to A_cont than A_1; trajectory")
        print("  is converging toward continuum target. G4-6a foundation")
        print("  CONSISTENT with multi-month commitment.")
    else:
        print("  [CAUTION] Trajectory is NOT cleanly monotonic toward A_cont.")
        print("  Possible explanations: small-t fit window includes intermediate")
        print("  regime where B constant dominates; substrate has different")
        print("  convergence behavior than assumed; additional panel needed.")

    results["richardson"] = {
        "A_1": float(A1),
        "A_2": float(A2),
        "a_1": float(a1),
        "a_2": float(a2),
        "diff_1": float(diff_1),
        "diff_2": float(diff_2),
        "p_implied": float(p_implied),
        "A_richardson_p1": float(A_richardson_p1),
        "A_richardson_implied": float(A_richardson_implied),
        "monotonic_toward_cont": bool(monotonic_toward_cont),
    }

    # Save
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
