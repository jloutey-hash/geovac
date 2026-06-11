"""Sprint G4-6a -- Multi-substrate UV foundation on SPECTRAL substrate.

Direct follow-on to debug/g4_6a_multi_substrate_uv_first_move.py which
ran the same panel sweep on the FD substrate and found 0.04% recovery
at t = a^2 (consistent with task #28 analytical prediction).

Per the G4-6d reframing (debug/g4_6a_multi_substrate_uv_first_move_memo.md):
spectral azimuthal discretization is the structural fix. This driver
reuses the panel architecture but with DiscreteWedgeDiracSpectral.

Expected per task #28: ~25.6% recovery at t = a^2 for spectral
substrate (vs 0.04% for FD).
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_6a_multi_substrate_uv_first_move_spectral.json"
OUT_JSON.parent.mkdir(exist_ok=True)


A_CONT = 1.0 / (24.0 * np.pi)


def compute_tip_term_at_t(N_rho: int, a: float, N_0: int, t: float,
                          k_step: int = 12) -> dict:
    """Compute tip term at a single t value via central FD.

    SPECTRAL version: uses DiscreteWedgeDiracSpectral instead of FD.
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
    X = np.vstack([1.0 / t_arr, np.ones_like(t_arr)]).T
    A, B = np.linalg.lstsq(X, tip_arr, rcond=None)[0]
    pred = A / t_arr + B
    rss = float(np.sum((tip_arr - pred) ** 2))
    tss = float(np.sum((tip_arr - tip_arr.mean()) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else float("nan")
    return float(A), float(B), r2


def run_panel(a: float, N_rho: int, N_0: int, label: str) -> dict:
    """Run one substrate panel: compute tip(t) at multiple small-t values."""
    print(f"\n{'=' * 72}")
    print(f"Panel {label} (SPECTRAL): a={a}, N_rho={N_rho}, N_0={N_0}, R={N_rho * a}")
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
    results["substrate"] = "spectral (exact m_eff)"

    print("=" * 72)
    print("Sprint G4-6a -- SPECTRAL substrate UV foundation FIRST MOVE")
    print("=" * 72)
    print(f"\n  A_cont = 1/(24*pi) = {A_CONT:.6f}")
    print(f"  Substrate: DiscreteWedgeDiracSpectral (exact m_eff = (k+1/2)/alpha)")
    print(f"  Compare to FD result: 0.04% recovery at t=a^2")

    panel_1 = run_panel(a=0.05, N_rho=200, N_0=120, label="P1_baseline_spec")
    panel_2 = run_panel(a=0.025, N_rho=400, N_0=240, label="P2_refined_spec")

    results["panels"] = [panel_1, panel_2]

    A1 = panel_1["A_fit"]
    A2 = panel_2["A_fit"]
    a1 = panel_1["a"]
    a2 = panel_2["a"]

    print(f"\n{'=' * 72}")
    print("Richardson extrapolation (SPECTRAL)")
    print(f"{'=' * 72}")
    print(f"\n  Panel 1 (a={a1}):  A = {A1:+.6f} ({A1/A_CONT*100:.2f}% of A_cont)")
    print(f"  Panel 2 (a={a2}):  A = {A2:+.6f} ({A2/A_CONT*100:.2f}% of A_cont)")

    diff_1 = A1 - A_CONT
    diff_2 = A2 - A_CONT

    if diff_1 != 0 and diff_2 != 0 and diff_1 * diff_2 > 0:
        p_implied = np.log2(diff_1 / diff_2) if diff_1 / diff_2 > 0 else float("nan")
        A_richardson_p1 = (2.0 * A2 - A1)
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

    print(f"\n  Richardson estimates of A(a -> 0):")
    print(f"    Assuming p=1: A_inf = 2*A_2 - A_1 = {A_richardson_p1:+.6f} "
          f"({A_richardson_p1/A_CONT*100:.2f}% of A_cont)")
    if not np.isnan(A_richardson_implied):
        print(f"    Using implied p={p_implied:.3f}: A_inf = {A_richardson_implied:+.6f} "
              f"({A_richardson_implied/A_CONT*100:.2f}% of A_cont)")

    print(f"\n{'=' * 72}")
    print("DECISION GATE (spectral vs FD)")
    print(f"{'=' * 72}")
    # FD result for comparison
    fd_path = Path(__file__).parent / "data" / "g4_6a_multi_substrate_uv_first_move.json"
    if fd_path.exists():
        with open(fd_path) as f:
            fd_results = json.load(f)
        fd_A2_recovery = fd_results["panels"][1]["A_recovery_vs_cont"]
        improvement = A2 / A_CONT - fd_A2_recovery
        print(f"  Spectral P2 recovery: {A2/A_CONT*100:.2f}% of A_cont")
        print(f"  FD       P2 recovery: {fd_A2_recovery*100:.2f}% of A_cont")
        print(f"  Improvement: {improvement*100:+.2f} percentage points")
        if A2 / A_CONT > fd_A2_recovery + 0.05:  # at least 5pp improvement
            print(f"\n  [POSITIVE] Spectral substrate substantially outperforms FD.")
            print(f"  G4-6d reframing CONFIRMED at production-code level.")
        else:
            print(f"\n  [CAUTION] Spectral does NOT substantially improve over FD.")
            print(f"  Need to investigate whether sub-percent panel range is in")
            print(f"  UV-divergence regime or intermediate-t regime.")

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
    }

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
