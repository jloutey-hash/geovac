"""G4-5a: Does the entropy coefficient B converge to 1/6 under UV refinement?

Existing data (G4-6b): B = 0.163 at a=0.05, N_rho=200, R=10. That's -2.3%
from continuum 1/6 = 0.16667. Richardson in 1/R^2 gave B_inf = 0.173 (+3.96%).

This driver tests UV convergence: decrease a (increase N_rho) at FIXED R=10,
measuring B at large t=10 where A/t << B. If B -> 1/6 as a -> 0, the gap is
purely discretization error and the substrate faithfully captures the entropy.

Panel:
  a = 0.10, 0.05, 0.025, 0.0125   (halving each step)
  N_rho = R/a = 100, 200, 400, 800
  N_phi = 120 (fixed azimuthal mode count)
  t = 10.0 (B-dominated regime)
  k_step = 12 (central FD for d/dalpha)

Expected: B(a) = 1/6 + C * a^p + ... with p ~ 2 (standard FD discretization).
Richardson extrapolation in a^2 should give B -> 1/6 if this is clean.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_tip_uv_convergence.json"
OUT_JSON.parent.mkdir(exist_ok=True)

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
    tip = dK_dalpha - K_disk
    return tip


def richardson_extrapolate(a_vals, B_vals, order=2):
    """Richardson extrapolation assuming B(a) = B_inf + C*a^order."""
    results = []
    for i in range(len(a_vals) - 1):
        a1, B1 = a_vals[i], B_vals[i]
        a2, B2 = a_vals[i + 1], B_vals[i + 1]
        r = (a1 / a2) ** order
        B_rich = (r * B2 - B1) / (r - 1)
        results.append(B_rich)
    return results


def main() -> None:
    print("=" * 72)
    print("G4-5a: Entropy coefficient B — UV convergence under a -> 0")
    print("=" * 72)

    R = 10.0
    N_phi = 120
    t = 10.0
    k_step = 12

    a_panel = [0.10, 0.05, 0.025, 0.0125]

    print(f"\n  Fixed: R = {R}, N_phi = {N_phi}, t = {t}, k_step = {k_step}")
    print(f"  Target: B = 1/6 = {B_TARGET:.8f}")
    print(f"  Panel: a = {a_panel}")
    print(f"  N_rho = {[int(R / a) for a in a_panel]}")
    print()

    results = []
    for a in a_panel:
        N_rho = int(R / a)
        t0 = time.time()
        B_est = compute_tip(N_rho=N_rho, a=a, N_phi=N_phi, t=t, k_step=k_step)
        elapsed = time.time() - t0
        err = (B_est - B_TARGET) / B_TARGET * 100
        results.append({
            "a": a,
            "N_rho": N_rho,
            "B_est": B_est,
            "B_target": B_TARGET,
            "rel_err_pct": err,
            "elapsed_s": elapsed,
        })
        print(f"  a={a:.4f}  N_rho={N_rho:4d}  B={B_est:.8f}  "
              f"err={err:+.4f}%  ({elapsed:.1f}s)")

    # Richardson extrapolation in a^2
    a_vals = [r["a"] for r in results]
    B_vals = [r["B_est"] for r in results]
    rich_a2 = richardson_extrapolate(a_vals, B_vals, order=2)

    print(f"\n  Richardson extrapolation (a^2 correction):")
    for i, B_r in enumerate(rich_a2):
        err = (B_r - B_TARGET) / B_TARGET * 100
        print(f"    pair ({a_vals[i]:.4f}, {a_vals[i+1]:.4f}): "
              f"B_rich = {B_r:.8f}  err = {err:+.4f}%")

    # Also try a^1 Richardson (in case the leading correction is O(a))
    rich_a1 = richardson_extrapolate(a_vals, B_vals, order=1)
    print(f"\n  Richardson extrapolation (a^1 correction):")
    for i, B_r in enumerate(rich_a1):
        err = (B_r - B_TARGET) / B_TARGET * 100
        print(f"    pair ({a_vals[i]:.4f}, {a_vals[i+1]:.4f}): "
              f"B_rich = {B_r:.8f}  err = {err:+.4f}%")

    # Convergence order estimate from consecutive ratios
    if len(B_vals) >= 3:
        print(f"\n  Convergence order estimate (from consecutive error ratios):")
        for i in range(len(B_vals) - 2):
            e1 = abs(B_vals[i] - B_TARGET)
            e2 = abs(B_vals[i + 1] - B_TARGET)
            e3 = abs(B_vals[i + 2] - B_TARGET)
            if e2 > 0 and e3 > 0:
                # e ~ C * a^p => log(e1/e2) / log(a1/a2) ~ p
                p12 = np.log(e1 / e2) / np.log(a_vals[i] / a_vals[i + 1])
                p23 = np.log(e2 / e3) / np.log(a_vals[i + 1] / a_vals[i + 2])
                print(f"    p({a_vals[i]:.4f}->{a_vals[i+1]:.4f}) = {p12:.3f}")
                print(f"    p({a_vals[i+1]:.4f}->{a_vals[i+2]:.4f}) = {p23:.3f}")

    # Also measure at t=5 and t=20 to check t-independence of B
    print(f"\n  B vs t check (at a=0.025, N_rho=400):")
    for t_check in [5.0, 10.0, 20.0, 50.0]:
        B_t = compute_tip(N_rho=400, a=0.025, N_phi=N_phi, t=t_check, k_step=k_step)
        err = (B_t - B_TARGET) / B_TARGET * 100
        print(f"    t={t_check:5.1f}  B={B_t:.8f}  err={err:+.4f}%")

    # Save
    output = {
        "description": "G4-5a UV convergence of entropy coefficient B under a -> 0",
        "parameters": {"R": R, "N_phi": N_phi, "t": t, "k_step": k_step},
        "B_target": B_TARGET,
        "panel": results,
        "richardson_a2": [{"pair": f"({a_vals[i]:.4f},{a_vals[i+1]:.4f})",
                           "B_rich": B_r}
                          for i, B_r in enumerate(rich_a2)],
        "richardson_a1": [{"pair": f"({a_vals[i]:.4f},{a_vals[i+1]:.4f})",
                           "B_rich": B_r}
                          for i, B_r in enumerate(rich_a1)],
    }
    OUT_JSON.write_text(json.dumps(output, indent=2))
    print(f"\n  Data saved: {OUT_JSON}")

    # Verdict
    print("\n" + "=" * 72)
    if len(rich_a2) >= 2:
        best_rich = rich_a2[-1]
        best_err = abs((best_rich - B_TARGET) / B_TARGET * 100)
        if best_err < 1.0:
            print("VERDICT: POSITIVE — B converges to 1/6 under UV refinement")
            print(f"  Best Richardson: {best_rich:.8f} ({best_err:.3f}% from 1/6)")
        elif best_err < 5.0:
            print("VERDICT: PARTIAL — B approaches 1/6 but slowly")
            print(f"  Best Richardson: {best_rich:.8f} ({best_err:.3f}% from 1/6)")
        else:
            print("VERDICT: NEGATIVE — B does not converge to 1/6 cleanly")
            print(f"  Best Richardson: {best_rich:.8f} ({best_err:.3f}% from 1/6)")
    print("=" * 72)


if __name__ == "__main__":
    main()
