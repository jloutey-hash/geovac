"""G4-5a part 2: Is the UV convergence uniform in alpha near alpha=1?

R2 Layer 3 (L6) requires: lim_n and d/dalpha commute.
This means the heat trace convergence must be uniform in alpha on
[1-delta, 1+delta]. Concretely: does B(a, alpha) converge to B_cont(alpha)
at the SAME rate for alpha near 1?

We compute tip(t, alpha) = K_wedge(t, alpha) - K_disk(t) at several
alpha values and check that the a-dependence is alpha-independent.

For the replica derivative: S_tip = d/dalpha|_{alpha=1} [K_wedge - K_disk]
The L6 question is whether this derivative commutes with a -> 0.
We test numerically by computing the finite-difference derivative at
several alpha offsets and checking stability.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_alpha_uniformity.json"

B_TARGET = 1.0 / 6.0


def compute_dK_dalpha(N_rho: int, a: float, N_phi: int, t: float,
                      alpha_center: float = 1.0,
                      dalpha: float = 0.1) -> float:
    """Compute dK/dalpha at alpha_center via central FD with step dalpha."""
    alpha_plus = alpha_center + dalpha
    alpha_minus = alpha_center - dalpha

    N_plus = int(round(N_phi * alpha_plus))
    N_minus = int(round(N_phi * alpha_minus))
    alpha_plus_actual = N_plus / N_phi
    alpha_minus_actual = N_minus / N_phi

    wedge_plus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus_actual,
    )
    wedge_minus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus_actual,
    )
    K_plus = wedge_plus.heat_trace(t)
    K_minus = wedge_minus.heat_trace(t)
    dK = (K_plus - K_minus) / (alpha_plus_actual - alpha_minus_actual)
    return dK


def compute_tip(N_rho: int, a: float, N_phi: int, t: float,
                dalpha: float = 0.1) -> float:
    """tip(t) = dK/dalpha|_{alpha=1} - K_disk."""
    disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi)
    K_disk = disk.heat_trace(t)
    dK = compute_dK_dalpha(N_rho, a, N_phi, t, alpha_center=1.0, dalpha=dalpha)
    return dK - K_disk


def main() -> None:
    print("=" * 72)
    print("G4-5a part 2: Alpha-uniformity of UV convergence (L6 test)")
    print("=" * 72)

    R = 10.0
    N_phi = 120
    t = 10.0

    # Test 1: Does the replica derivative converge as a -> 0?
    # Use different dalpha values to check FD step stability
    print("\n  === Test 1: FD step convergence at fixed a=0.025 ===")
    a = 0.025
    N_rho = int(R / a)
    dalpha_panel = [0.20, 0.15, 0.10, 0.05]
    print(f"  a={a}, N_rho={N_rho}, t={t}")
    for da in dalpha_panel:
        tip = compute_tip(N_rho, a, N_phi, t, dalpha=da)
        err = (tip - B_TARGET) / B_TARGET * 100
        print(f"    dalpha={da:.2f}  tip={tip:.8f}  err={err:+.4f}%")

    # Test 2: UV convergence of the derivative at fixed dalpha
    print("\n  === Test 2: UV convergence of tip at dalpha=0.10 ===")
    a_panel = [0.05, 0.025, 0.0125]
    dalpha = 0.10
    results_main = []
    for a in a_panel:
        N_rho = int(R / a)
        t0 = time.time()
        tip = compute_tip(N_rho, a, N_phi, t, dalpha=dalpha)
        elapsed = time.time() - t0
        err = (tip - B_TARGET) / B_TARGET * 100
        results_main.append({"a": a, "N_rho": N_rho, "tip": tip, "err_pct": err})
        print(f"    a={a:.4f}  N_rho={N_rho:4d}  tip={tip:.8f}  "
              f"err={err:+.4f}%  ({elapsed:.1f}s)")

    # Test 3: Compute K_wedge at several alpha values, check that
    # (K_wedge(a, alpha) - K_wedge_cont(alpha)) / a^p is alpha-independent
    print("\n  === Test 3: K_wedge convergence at different alpha values ===")
    alpha_panel = [0.90, 0.95, 1.00, 1.05, 1.10]
    a_panel_3 = [0.05, 0.025]

    print(f"  (Checking heat trace K_wedge at each alpha, at two a values)")
    print(f"  {'alpha':>6}  {'K(a=0.05)':>12}  {'K(a=0.025)':>12}  {'ratio':>8}  "
          f"{'diff':>12}")
    results_alpha = []
    for alpha in alpha_panel:
        Ks = []
        for a in a_panel_3:
            N_rho = int(R / a)
            N_w = int(round(N_phi * alpha))
            alpha_actual = N_w / N_phi
            wedge = DiscreteWedgeDiracSpectral(
                N_rho=N_rho, a=a, N_phi=N_w, alpha=alpha_actual,
            )
            K = wedge.heat_trace(t)
            Ks.append(K)
        diff = Ks[1] - Ks[0]
        ratio = diff / Ks[0] if Ks[0] != 0 else float('inf')
        results_alpha.append({
            "alpha": alpha,
            "K_coarse": Ks[0], "K_fine": Ks[1],
            "diff": diff, "ratio": ratio
        })
        print(f"  {alpha:6.2f}  {Ks[0]:12.6f}  {Ks[1]:12.6f}  "
              f"{ratio:8.5f}  {diff:12.6f}")

    # Check: is the ratio (K_fine - K_coarse)/K roughly alpha-independent?
    ratios = [r["ratio"] for r in results_alpha]
    ratio_cv = np.std(ratios) / np.mean(ratios) if np.mean(ratios) != 0 else 0
    print(f"\n  Ratio CV across alpha values: {ratio_cv:.4f}")
    print(f"  (CV < 0.1 = strongly uniform; CV > 0.5 = alpha-dependent)")

    # Test 4: Second derivative d²K/dalpha² at alpha=1 (convergence check)
    # This tests whether the derivative has convergent structure
    print("\n  === Test 4: d²K/dalpha² at alpha=1 (second replica derivative) ===")
    a_panel_4 = [0.05, 0.025]
    dalpha = 0.10
    for a in a_panel_4:
        N_rho = int(R / a)
        disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi)
        K_disk = disk.heat_trace(t)

        # Three-point stencil for d²K/dalpha²
        N_plus = int(round(N_phi * (1.0 + dalpha)))
        N_minus = int(round(N_phi * (1.0 - dalpha)))
        alpha_plus = N_plus / N_phi
        alpha_minus = N_minus / N_phi
        N_center = N_phi
        alpha_center = 1.0

        wedge_plus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus)
        wedge_minus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus)
        wedge_center = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_center, alpha=alpha_center)

        K_p = wedge_plus.heat_trace(t)
        K_m = wedge_minus.heat_trace(t)
        K_c = wedge_center.heat_trace(t)

        d2K = (K_p - 2*K_c + K_m) / ((alpha_plus - alpha_center)**2)
        print(f"    a={a:.4f}  d²K/dalpha² = {d2K:.6f}")

    # Save
    output = {
        "description": "G4-5a alpha-uniformity test for L6 (lim and d/dalpha commute)",
        "test3_alpha_ratios": results_alpha,
        "test3_ratio_CV": ratio_cv,
        "test2_uv_convergence": results_main,
    }
    OUT_JSON.write_text(json.dumps(output, indent=2, default=str))
    print(f"\n  Data saved: {OUT_JSON}")

    # Verdict
    print("\n" + "=" * 72)
    if ratio_cv < 0.1:
        print("VERDICT: STRONG POSITIVE — UV convergence is uniform in alpha")
        print(f"  K(fine)-K(coarse) ratio CV = {ratio_cv:.4f} across alpha panel")
        print("  L6 (lim and d/dalpha commute) is numerically supported")
    elif ratio_cv < 0.3:
        print("VERDICT: MODERATE POSITIVE — UV convergence weakly alpha-dependent")
        print(f"  Ratio CV = {ratio_cv:.4f}")
    else:
        print("VERDICT: NEGATIVE — UV convergence is alpha-dependent")
        print(f"  Ratio CV = {ratio_cv:.4f}")
    print("=" * 72)


if __name__ == "__main__":
    main()
