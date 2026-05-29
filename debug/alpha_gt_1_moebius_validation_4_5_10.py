"""Task #25 -- alpha > 1 Moebius ansatz validation at alpha in {4, 5, 10}.

Closes the empirical validation gate for the v3.19.0 Track 5 Moebius
closed form against new alpha values NOT in the original fit set
(which was alpha in {1.5, 2.0, 3.0}).

Closed form (v3.19.0 Track 5):
    slope^{Dirac, ex}(alpha) = -(1/12) * alpha / (2*alpha - 1)   for alpha > 1
    Delta_K^{Dirac, ex}(alpha) = (alpha^2 - 1) / (24 * (alpha - 1/2))

Predicted slopes and recoveries at the new alpha values:
    alpha=4:   slope = -1/12 * 4/7  = -0.047619    recovery = 4/7  = 0.5714  (57.14%)
    alpha=5:   slope = -1/12 * 5/9  = -0.046296    recovery = 5/9  = 0.5556  (55.56%)
    alpha=10:  slope = -1/12 * 10/19 = -0.043860   recovery = 10/19 = 0.5263  (52.63%)

Decision gate
-------------
- POSITIVE: average relative error < 3% across {4, 5, 10}
- POSITIVE-with-caveat: avg rel err < 5%
- NEGATIVE: avg rel err >= 5% (three-point fit doesn't extrapolate)

Substrate (matches G4-4c week 2, N_0-independence verified there):
- R = 10, a = 0.05, N_rho = 200, N_0 = 120, t = 1.0
- alpha in {4, 5, 10}
- N_phi = alpha * N_0 = {480, 600, 1200}

Files:
- debug/data/alpha_gt_1_moebius_validation_4_5_10.json
- debug/alpha_gt_1_moebius_validation_4_5_10_memo.md (separate)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

OUT_JSON = Path(__file__).parent / "data" / "alpha_gt_1_moebius_validation_4_5_10.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main():
    print("=" * 76)
    print("Task #25 -- alpha > 1 Moebius validation at alpha in {4, 5, 10}")
    print("=" * 76)

    # Substrate
    R = 10.0
    a = 0.05
    N_rho = 200
    N_0 = 120
    t_eval = 1.0
    alphas = [4.0, 5.0, 10.0]

    # Continuum target
    SC_continuum = -1.0 / 12.0  # standard SC coefficient at alpha < 1

    print(f"\n[Substrate]")
    print(f"  R = {R}, a = {a}, N_rho = {N_rho}, N_0 = {N_0}, t = {t_eval}")
    print(f"  alpha values: {alphas}")
    print(f"  N_phi values: {[int(alpha * N_0) for alpha in alphas]}")
    print(f"  Continuum SC target: -1/12 = {SC_continuum:.10f}")

    # ------------------------------------------------------------------
    # Step 0 -- disk baseline
    # ------------------------------------------------------------------
    print(f"\n[Step 0] Disk-Dirac baseline at N_0 = {N_0} ...")
    t_start = time.time()
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    K_disk = float(disk.heat_trace(t_eval))
    t_disk = time.time() - t_start
    print(f"  K_disk(t={t_eval}) = {K_disk:.6f}")
    print(f"  elapsed: {t_disk:.1f} s")

    # ------------------------------------------------------------------
    # Step 1 -- measure slopes at alpha in {4, 5, 10}
    # ------------------------------------------------------------------
    print(f"\n[Step 1] Measure slopes at alpha in {alphas}")

    results = {
        "substrate": {
            "R": R, "a": a, "N_rho": N_rho, "N_0": N_0, "t_eval": t_eval,
            "alphas": alphas,
            "K_disk": K_disk,
            "SC_continuum_target": SC_continuum,
        },
        "ansatz": {
            "formula_slope": "-(1/12) * alpha / (2*alpha - 1)",
            "formula_Delta_K": "(alpha^2 - 1) / (24 * (alpha - 1/2))",
            "F_moebius": "F(alpha) = alpha / (2*alpha - 1)",
            "asymptote": "F(alpha) -> 1/2 as alpha -> infinity",
        },
        "data": [],
    }

    for alpha in alphas:
        N_phi = int(round(alpha * N_0))
        print(f"\n  alpha = {alpha:.1f}, N_phi = {N_phi}")
        t_start = time.time()
        wedge = DiscreteWedgeDirac(N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha)
        K_wedge = float(wedge.heat_trace(t_eval))
        t_w = time.time() - t_start
        # Slope = Delta_K / (1/alpha - alpha) where Delta_K = K_wedge - alpha*K_disk
        # Equivalent at the SC level to Delta_K = K_wedge - K_disk_at_alpha_x_bulk
        # but the standard convention is K_wedge - K_disk*alpha (bulk normalization
        # so the bulk piece scales linearly with alpha).
        # At alpha=1: K_wedge = K_disk, Delta_K = 0, slope undefined (1/alpha - alpha = 0).
        # G4-4c week 2 used Delta_K = K_wedge - alpha * K_disk.
        Delta_K = K_wedge - alpha * K_disk
        denom = 1.0 / alpha - alpha
        slope_meas = Delta_K / denom

        # Moebius prediction
        slope_pred = -1.0 / 12.0 * alpha / (2.0 * alpha - 1.0)
        Delta_K_pred = (alpha ** 2 - 1.0) / (24.0 * (alpha - 0.5))

        # Recovery
        recovery_meas = slope_meas / SC_continuum
        recovery_pred = alpha / (2.0 * alpha - 1.0)

        # Relative errors
        rel_err_slope = (slope_pred - slope_meas) / slope_meas
        rel_err_Delta = (Delta_K_pred - Delta_K) / Delta_K
        rel_err_recovery = (recovery_pred - recovery_meas) / recovery_meas

        print(f"    K_wedge({t_eval}) = {K_wedge:.6f}  ({t_w:.1f} s)")
        print(f"    Delta_K = {Delta_K:.6f}, denom (1/a - a) = {denom:.6f}")
        print(f"    slope measured  = {slope_meas:.6f}")
        print(f"    slope predicted = {slope_pred:.6f}")
        print(f"    rel err slope   = {100 * rel_err_slope:+.3f}%")
        print(f"    recovery measured = {recovery_meas:.4f} ({100*recovery_meas:.2f}%)")
        print(f"    recovery predicted = {recovery_pred:.4f} ({100*recovery_pred:.2f}%)")

        results["data"].append({
            "alpha": alpha,
            "N_phi": N_phi,
            "K_wedge": K_wedge,
            "Delta_K_measured": Delta_K,
            "Delta_K_predicted": Delta_K_pred,
            "slope_measured": slope_meas,
            "slope_predicted": slope_pred,
            "recovery_measured": recovery_meas,
            "recovery_predicted": recovery_pred,
            "rel_err_slope_pct": 100 * rel_err_slope,
            "rel_err_Delta_K_pct": 100 * rel_err_Delta,
            "rel_err_recovery_pct": 100 * rel_err_recovery,
            "elapsed_s": t_w,
        })

    # ------------------------------------------------------------------
    # Step 2 -- verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 2] Summary and verdict")
    print("=" * 76)
    print()
    print(f"  {'alpha':>6}  {'slope_meas':>12}  {'slope_pred':>12}  "
          f"{'rel err':>10}  {'rec meas':>10}  {'rec pred':>10}")
    print("  " + "-" * 72)
    rel_errs = []
    for d in results["data"]:
        print(f"  {d['alpha']:>6.1f}  {d['slope_measured']:>+12.6f}  "
              f"{d['slope_predicted']:>+12.6f}  "
              f"{d['rel_err_slope_pct']:>+7.3f}%   "
              f"{d['recovery_measured']:>10.4f}  "
              f"{d['recovery_predicted']:>10.4f}")
        rel_errs.append(abs(d['rel_err_slope_pct']))

    mean_rel_err = float(np.mean(rel_errs))
    max_rel_err = float(np.max(rel_errs))

    print(f"\n  Mean |rel err|: {mean_rel_err:.3f}%")
    print(f"  Max  |rel err|: {max_rel_err:.3f}%")
    print(f"\n  v3.19.0 original fit set (alpha in {{1.5, 2, 3}}): avg 2.3%")
    print(f"  v3.19.0 decision gate: 10%")
    print(f"  Task #25 decision gate: 3% (POSITIVE), 5% (POSITIVE-with-caveat)")

    if mean_rel_err < 3.0:
        verdict = "POSITIVE-EMPIRICAL-LOCK"
        msg = (f"Moebius ansatz extrapolates cleanly from the original fit set "
               f"{{1.5, 2, 3}} to the new alpha values {{4, 5, 10}}. Mean rel err "
               f"{mean_rel_err:.2f}% < 3% gate. v3.19.0 Track 5 headline result "
               f"empirically locked as a structural identification.")
    elif mean_rel_err < 5.0:
        verdict = "POSITIVE-WITH-CAVEAT"
        msg = (f"Mean rel err {mean_rel_err:.2f}% < 5% gate but exceeds strict "
               f"3%. Moebius ansatz extrapolates well but with finite-t "
               f"corrections growing at large alpha. Structural identification "
               f"holds; quantitative match wider than at the fit points.")
    else:
        verdict = "NEGATIVE-FIT-DOESNT-EXTRAPOLATE"
        msg = (f"Mean rel err {mean_rel_err:.2f}% exceeds 5% gate. v3.19.0 "
               f"three-point fit may be Diophantine coincidence rather than "
               f"structural identification. Needs further investigation.")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["mean_rel_err_pct"] = mean_rel_err
    results["max_rel_err_pct"] = max_rel_err
    results["verdict"] = verdict
    results["msg"] = msg

    # Save
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\n  Results saved to {OUT_JSON}")
    print("=" * 76)


if __name__ == "__main__":
    main()
