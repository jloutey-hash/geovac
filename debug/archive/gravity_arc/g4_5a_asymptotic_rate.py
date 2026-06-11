"""G4-5a rate sharpening: Does tip convergence reach O(a^2) asymptotically?

The Part 1 panel showed p accelerating: 1.08 -> 1.23 -> 1.70.
This suggests approaching O(a^2) from below. Adding a=0.00625 (N_rho=1600)
should confirm whether the rate stabilizes at p=2.

The structural prediction: the leading O(a) correction is alpha-independent
(cancels in the tip), and the next correction should be O(a^2) from standard
centered-difference FD. So the tip rate should be asymptotically O(a^2).

If confirmed, the R2 theorem rate upgrades from O(a^1.7) (observed at current
scale) to O(a^2) (asymptotic, stronger statement for the paper).
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_asymptotic_rate.json"

B_TARGET = 1.0 / 6.0


def compute_tip(N_rho: int, a: float, N_phi: int, t: float,
                k_step: int = 12) -> float:
    """tip(t) = dK/dalpha|_{alpha=1} - K_disk via central FD."""
    N_0 = N_phi
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
    return dK_dalpha - K_disk


def main() -> None:
    print("=" * 72)
    print("G4-5a: Asymptotic rate — does tip convergence reach O(a^2)?")
    print("=" * 72)

    R = 10.0
    N_phi = 120
    t = 10.0
    k_step = 12

    # Extended panel including finer grid
    a_panel = [0.10, 0.05, 0.025, 0.0125, 0.00625]
    N_rho_panel = [int(R / a) for a in a_panel]

    print(f"\n  R={R}, N_phi={N_phi}, t={t}")
    print(f"  Panel: a = {a_panel}")
    print(f"  N_rho = {N_rho_panel}")
    print(f"  Target: B = 1/6 = {B_TARGET:.10f}")
    print()

    results = []
    for a in a_panel:
        N_rho = int(R / a)
        t0 = time.time()
        B_est = compute_tip(N_rho=N_rho, a=a, N_phi=N_phi, t=t, k_step=k_step)
        elapsed = time.time() - t0
        err_abs = B_est - B_TARGET
        err_pct = err_abs / B_TARGET * 100
        results.append({
            "a": a, "N_rho": N_rho, "B_est": B_est,
            "err_abs": err_abs, "err_pct": err_pct, "elapsed_s": elapsed,
        })
        print(f"  a={a:.5f}  N_rho={N_rho:5d}  B={B_est:.10f}  "
              f"err={err_pct:+.5f}%  ({elapsed:.1f}s)")

    # Convergence order from consecutive error ratios
    print(f"\n  Convergence order (from consecutive absolute errors):")
    errs = [abs(r["err_abs"]) for r in results]
    for i in range(len(errs) - 1):
        ratio = errs[i] / errs[i+1] if errs[i+1] > 0 else float('inf')
        p = np.log2(ratio) if ratio > 0 else 0
        print(f"    a={a_panel[i]:.5f} -> {a_panel[i+1]:.5f}: "
              f"err_ratio={ratio:.4f}, p={p:.4f}")

    # Richardson extrapolation at order 2 (all consecutive pairs)
    print(f"\n  Richardson (a^2):")
    for i in range(len(results) - 1):
        a1, B1 = results[i]["a"], results[i]["B_est"]
        a2, B2 = results[i+1]["a"], results[i+1]["B_est"]
        r = (a1 / a2) ** 2
        B_rich = (r * B2 - B1) / (r - 1)
        err_rich = (B_rich - B_TARGET) / B_TARGET * 100
        print(f"    ({a1:.5f}, {a2:.5f}): B_rich={B_rich:.10f}  err={err_rich:+.6f}%")

    # Also try Richardson at order 1.7 (empirical from earlier)
    print(f"\n  Richardson (a^1.7, empirical):")
    for i in range(len(results) - 1):
        a1, B1 = results[i]["a"], results[i]["B_est"]
        a2, B2 = results[i+1]["a"], results[i+1]["B_est"]
        r = (a1 / a2) ** 1.7
        B_rich = (r * B2 - B1) / (r - 1)
        err_rich = (B_rich - B_TARGET) / B_TARGET * 100
        print(f"    ({a1:.5f}, {a2:.5f}): B_rich={B_rich:.10f}  err={err_rich:+.6f}%")

    # Save
    output = {
        "description": "G4-5a asymptotic rate test — extended to a=0.00625",
        "parameters": {"R": R, "N_phi": N_phi, "t": t},
        "panel": results,
    }
    OUT_JSON.write_text(json.dumps(output, indent=2))
    print(f"\n  Data saved: {OUT_JSON}")

    # Verdict
    print("\n" + "=" * 72)
    if len(errs) >= 5:
        final_ratio = errs[-2] / errs[-1] if errs[-1] > 0 else 0
        final_p = np.log2(final_ratio) if final_ratio > 0 else 0
        if final_p >= 1.9:
            print(f"VERDICT: O(a^2) CONFIRMED — finest-pair p = {final_p:.3f}")
        elif final_p >= 1.5:
            print(f"VERDICT: APPROACHING O(a^2) — finest-pair p = {final_p:.3f}")
            print(f"  (still accelerating; need one more halving to confirm)")
        else:
            print(f"VERDICT: SUB-QUADRATIC — finest-pair p = {final_p:.3f}")
    print("=" * 72)


if __name__ == "__main__":
    main()
