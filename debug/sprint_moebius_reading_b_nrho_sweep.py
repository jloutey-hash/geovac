"""Reading B verification -- N_rho sweep for Moebius factor.

The Moebius factor F(alpha) = alpha/(2alpha-1) was empirically verified:
- N_0-independent at N_0 in {120, 240, 480} (G4-4c week 2)
- alpha-extrapolatable to alpha in {4, 5, 10} (v3.20.0 task #25)

What is NOT yet verified is N_rho-independence (the radial refinement).

Reading B (substrate-universal) prediction: the Moebius slope should
be invariant under N_rho refinement, since the factor is tied to the
azimuthal mode structure (anti-periodic spinor BC, lowest eigenvalue
1/(2alpha)), NOT to the radial discretization.

Setup:
- Fixed alpha = 2 (the Moebius prediction is -1/18 ~= -0.0556)
- Fixed a = 0.05 (matching baseline)
- Fixed N_0 = 120 (matching G4-4c week 2 baseline)
- Vary N_rho in {100, 200, 400, 600}
- Measure the slope d(Delta_K)/d(alpha) at alpha = 2 via central FD

If slope is N_rho-independent at -1/18, Reading B is supported across
radial refinement too.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDiracSpectral,
)

OUT_JSON = Path(__file__).parent / "data" / "sprint_moebius_reading_b_nrho_sweep.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def measure_moebius_slope(N_rho: int, a: float, N_0: int, alpha: float,
                          t: float, eps: float = 0.05) -> dict:
    """Measure d(tip)/d(alpha) at given alpha via central FD over alpha range.

    The Moebius prediction is: slope at alpha = -1/12 * alpha/(2*alpha-1).
    """
    alpha_plus = alpha + eps
    alpha_minus = alpha - eps

    # For each alpha, N_phi must be alpha * N_0 rounded to integer
    N_plus = max(2, int(round(alpha_plus * N_0)))
    N_minus = max(2, int(round(alpha_minus * N_0)))
    # The actual alphas after rounding
    actual_alpha_plus = N_plus / N_0
    actual_alpha_minus = N_minus / N_0

    # Compute tip at each alpha via inner central FD over wedge mode count
    # but use small inner k_step for the tip at each alpha
    # Actually let me use direct: at this alpha, tip is just the value
    # of the wedge minus the bulk piece alpha*K_disk

    disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_0)
    K_disk = disk.heat_trace(t)

    # tip at alpha_plus
    wedge_plus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=actual_alpha_plus,
    )
    K_wedge_plus = wedge_plus.heat_trace(t)
    tip_plus = K_wedge_plus - actual_alpha_plus * K_disk

    # tip at alpha_minus
    wedge_minus = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=actual_alpha_minus,
    )
    K_wedge_minus = wedge_minus.heat_trace(t)
    tip_minus = K_wedge_minus - actual_alpha_minus * K_disk

    # Central FD slope: d(tip)/d(alpha)
    slope = (tip_plus - tip_minus) / (actual_alpha_plus - actual_alpha_minus)

    # Moebius prediction at this alpha
    moebius_pred = -(1.0 / 12.0) * alpha / (2.0 * alpha - 1.0)

    return {
        "alpha": alpha,
        "actual_alpha_plus": actual_alpha_plus,
        "actual_alpha_minus": actual_alpha_minus,
        "tip_plus": tip_plus,
        "tip_minus": tip_minus,
        "slope_measured": slope,
        "moebius_pred": moebius_pred,
        "rel_error": (slope - moebius_pred) / moebius_pred if moebius_pred != 0 else float("nan"),
    }


def main() -> None:
    print("=" * 72)
    print("Reading B verification -- N_rho sweep for Moebius factor")
    print("=" * 72)

    a = 0.05
    N_0 = 120
    alpha = 2.0  # Excess angle where Moebius is non-trivial
    t = 1.0  # Sweet-spot t window matching G4-4c week 2

    moebius_pred = -(1.0 / 12.0) * alpha / (2.0 * alpha - 1.0)
    print(f"\n  Fixed: a={a}, N_0={N_0}, alpha={alpha}, t={t}")
    print(f"  Moebius prediction at alpha=2: -1/12 * 2/3 = -1/18 = {moebius_pred:.6f}")

    N_rho_values = [100, 200, 400, 600]
    print(f"  N_rho sweep: {N_rho_values}")

    sweep_data = []
    for N_rho in N_rho_values:
        print(f"\n{'=' * 72}")
        print(f"N_rho = {N_rho}, R = {N_rho * a}")
        print(f"{'=' * 72}")
        t0 = time.time()
        result = measure_moebius_slope(N_rho, a, N_0, alpha, t)
        sweep_data.append({"N_rho": N_rho, "R": N_rho * a, **result})
        elapsed = time.time() - t0
        print(f"  alpha_plus={result['actual_alpha_plus']:.4f}, "
              f"alpha_minus={result['actual_alpha_minus']:.4f}")
        print(f"  tip(alpha_plus)  = {result['tip_plus']:+.6f}")
        print(f"  tip(alpha_minus) = {result['tip_minus']:+.6f}")
        print(f"  slope measured = {result['slope_measured']:+.6f}")
        print(f"  Moebius pred   = {result['moebius_pred']:+.6f}")
        print(f"  rel error      = {result['rel_error']*100:+.3f}%")
        print(f"  compute time   = {elapsed:.1f}s")

    # ---------------------------------------------------------------------
    # Analysis
    # ---------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("N_rho invariance check (Reading B verification)")
    print(f"{'=' * 72}")
    slopes = [s["slope_measured"] for s in sweep_data]
    rel_errs = [s["rel_error"] for s in sweep_data]

    print(f"\n  {'N_rho':>6}  {'slope':>10}  {'Moebius pred':>14}  {'rel err':>10}")
    for s in sweep_data:
        print(f"  {s['N_rho']:>6d}  {s['slope_measured']:>+10.6f}  "
              f"{s['moebius_pred']:>+14.6f}  {s['rel_error']*100:>+9.3f}%")

    # Variation in slope across N_rho
    slope_arr = np.array(slopes)
    slope_mean = float(np.mean(slope_arr))
    slope_std = float(np.std(slope_arr))
    slope_cv = slope_std / abs(slope_mean) if slope_mean != 0 else float("nan")
    print(f"\n  Slope mean across N_rho: {slope_mean:+.6f}")
    print(f"  Slope std:               {slope_std:.6f}")
    print(f"  Coefficient of variation: {slope_cv*100:.3f}%")

    # Mean rel err vs Moebius
    mean_rel_err = float(np.mean(np.abs(rel_errs)))
    print(f"\n  Mean |rel err| vs Moebius: {mean_rel_err*100:.3f}%")

    # Verdict
    print(f"\n{'=' * 72}")
    print("Verdict")
    print(f"{'=' * 72}")
    # Reading B passes if slope is N_rho-stable AND matches Moebius
    n_rho_stable = slope_cv < 0.05  # Less than 5% coefficient of variation
    moebius_match = mean_rel_err < 0.05  # Within 5% of Moebius prediction

    if n_rho_stable and moebius_match:
        verdict = "POSITIVE-READING-B-SUPPORTED"
        print(f"\n  [{verdict}] Moebius slope is N_rho-stable (CV < 5%) AND")
        print(f"  matches Moebius prediction (mean rel err < 5%) across N_rho.")
        print(f"  Reading B (substrate-universal) is strongly supported:")
        print(f"  Moebius factor is invariant under radial refinement, consistent")
        print(f"  with substrate-level mechanism tied to azimuthal lowest eigenvalue.")
    elif n_rho_stable:
        verdict = "PARTIAL-N_RHO-STABLE-MOEBIUS-MISS"
        print(f"\n  [{verdict}] Slope is N_rho-stable but doesn't match Moebius.")
        print(f"  May indicate the substrate has a different stable mechanism.")
    elif moebius_match:
        verdict = "PARTIAL-MOEBIUS-MATCH-VARIABLE"
        print(f"\n  [{verdict}] Slope matches Moebius on average but varies with N_rho.")
    else:
        verdict = "NEUTRAL"
        print(f"\n  [{verdict}] Neither N_rho stability nor Moebius match.")

    results = {
        "alpha": alpha,
        "t": t,
        "a": a,
        "N_0": N_0,
        "moebius_pred": moebius_pred,
        "N_rho_sweep": sweep_data,
        "slope_mean": slope_mean,
        "slope_std": slope_std,
        "slope_cv": slope_cv,
        "mean_rel_err_vs_moebius": mean_rel_err,
        "n_rho_stable": n_rho_stable,
        "moebius_match": moebius_match,
        "verdict": verdict,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
