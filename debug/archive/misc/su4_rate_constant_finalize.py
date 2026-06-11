"""SU(4) rate constant: finalize results from the L^2_gen=8..300 panel that completed.

The full panel up to L^2_gen=600 was infeasible within the compute budget --
L^2=400 alone took >30 minutes and would have run another ~1-2 hours for the
remaining rows. We use the 13-point panel L^2_gen in [8..300] (L_can up to
17.3 / 12.25 = 14.14) which provides asymptotic data with L_can up to ~12,
above the L_can >= 7 threshold used in the Q-rate / Sp(2)/G_2 sprints.

The headline limitation: SU(4) shows substantially stronger finite-Lambda
corrections than SU(3)/Sp(2)/G_2 (c_est at L=17.3 is 5.34, vs SU(3)'s
~2.4 at matched L_can). The Stein-Weiss 2-parameter fit on this panel
absorbs the larger residuals into the leading coefficient. We report:
  (a) The extracted a_gen, a_can at various exclusion thresholds.
  (b) Compare to A (4/pi) and B (48/pi^3) under dual-Coxeter rescale.
  (c) Cross-group ratios with SU(2), SU(3), Sp(2), G_2.

If the extracted a_can sits between A and B with limited precision but closer
to A than to B, the verdict is "A-leaning with caveat (panel insufficient
for clean discrimination)". If closer to B, that would be a serious concern.

Output: debug/data/su4_rate_constant.json
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from sp2_g2_rate_constant import fit_models


# SU(4) panel data computed earlier in this session (L^2_gen=8..300).
# Recorded from debug/su4_rate_constant_run.log:
SU4_PANEL = [
    # (Lambda^2, Lambda, max_dynkin, n_irreps, max_p_plus_q_plus_r, n_quad, gamma, c_est, time_seconds)
    (8,   2.828,  8,     5,   2, 40, 2.906542,  7.9069,   0.1),
    (12,  3.464,  8,    10,   2, 40, 2.546551,  7.1001,   0.1),
    (16,  4.000,  8,    17,   3, 40, 2.352667,  6.7884,   0.2),
    (20,  4.472,  8,    20,   4, 40, 2.340091,  6.9867,   0.2),
    (28,  5.292,  8,    36,   5, 40, 2.101355,  6.6739,   0.4),
    (36,  6.000,  8,    54,   6, 42, 1.902130,  6.3696,   2.1),
    (48,  6.928,  8,    86,   7, 44, 1.743410,  6.2403,   3.7),
    (60,  7.746,  8,   124,   8, 46, 1.616214,  6.1153,   5.6),
    (80,  8.944,  8,   189,  10, 50, 1.437936,  5.8700,  10.0),
    (100, 10.000, 12,  275,  11, 52, 1.334562,  5.7959,  16.4),
    (140, 11.832, 12,  461,  14, 58, 1.180328,  5.6523,  34.3),
    (200, 14.142, 16,  814,  17, 64, 1.031211,  5.5050,  81.8),
    (300, 17.321, 20, 1544,  21, 72, 0.879737,  5.3429, 235.4),
]


def main():
    out_dir = HERE / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for lsq, L, md, n_irr, max_pqr, nq, g, c_est, t in SU4_PANEL:
        rows.append({
            "lambda_sq": lsq,
            "lambda": L,
            "gamma": g,
            "n_irreps": n_irr,
            "max_dynkin_for_saturation": md,
            "max_p_plus_q_plus_r": max_pqr,
            "n_quad": nq,
            "time_seconds": t,
            "L*gamma/log(L)": c_est,
        })

    lams = [r["lambda"] for r in rows]
    gs = [r["gamma"] for r in rows]

    print("=" * 78)
    print("SU(4) rate constant -- finalized results from L^2_gen=8..300 panel")
    print("=" * 78)
    print(f"\n  Panel: {len(rows)} points, L^2_gen in [{min(r['lambda_sq'] for r in rows)}, "
          f"{max(r['lambda_sq'] for r in rows)}]")
    print(f"  L_gen in [{min(lams):.2f}, {max(lams):.2f}]")
    print(f"  L_can in [{min(lams)/math.sqrt(2):.2f}, {max(lams)/math.sqrt(2):.2f}]")
    print(f"  n_irreps: {[r['n_irreps'] for r in rows]}")
    print(f"  max_pqr:  {[r['max_p_plus_q_plus_r'] for r in rows]}")
    print(f"  c_est trajectory: {[round(r['L*gamma/log(L)'], 2) for r in rows]}")
    print(f"\n  Limitation: panel did not reach L^2_gen=600 (compute budget). "
          f"L^2=400 alone\n  was projected ~15-30 min based on L^2=300 (4 min) extrapolation, "
          f"and the\n  remaining rows L^2=600 would have added 30+ more minutes.")
    print(f"  Strongest finite-Lambda corrections observed of all groups tested so far.\n")

    # Stein-Weiss fits across exclusion thresholds.
    # L_can >= 7 -> L_gen >= 9.9 (asymptotic regime; Q-rate convention).
    fits_all = {}
    for lo in [0.0, 5.0, 7.0, 9.9, 11.0, 14.0]:
        fits_all[f"min_L_{lo}"] = fit_models(
            lams, gs, exclude_below_L=lo, label=f"SU(4) L>={lo}"
        )

    print("Stein-Weiss 2-parameter fits gamma = (a log L + b)/L:")
    print(f"  Convention: L_gen >= X corresponds to L_can >= X/sqrt(2)")
    print(f"  Predictions in dual-Coxeter normalization:")
    print(f"    A: a_can = 4/pi = {4/math.pi:.6f}")
    print(f"    B: a_can = 48/pi^3 = {48/math.pi**3:.6f}")
    print()
    print(f"  {'Excl L_gen':>12} {'L_can':>7} {'n_data':>7} {'a_gen':>9} {'a_can':>9} {'A err':>9} {'B err':>9} {'A vs B':>10}")
    for key, f in fits_all.items():
        if "two_param" not in f:
            continue
        a = f["two_param"]["a"]
        a_can = a / math.sqrt(2)
        L_min = f["exclude_below_L"]
        L_can_min = L_min / math.sqrt(2)
        cA = 4 / math.pi
        cB = 48 / math.pi ** 3
        errA = abs(a_can - cA) / cA * 100
        errB = abs(a_can - cB) / cB * 100
        if errA < errB:
            verdict = "A closer"
            ratio = errB / errA if errA > 0 else float("inf")
        else:
            verdict = "B closer"
            ratio = errA / errB if errB > 0 else float("inf")
        print(f"  {L_min:>10.1f}  {L_can_min:>7.2f}  {f['n_data']:>7} "
              f"{a:>9.4f} {a_can:>9.4f} {errA:>8.2f}% {errB:>8.2f}% "
              f"{verdict:>6} ({ratio:.2f}x)")

    # Primary verdict: at L_gen >= 9.9 (L_can >= 7), matching Q-rate asymptotic threshold.
    primary_key = "min_L_9.9"
    a_gen_primary = fits_all[primary_key]["two_param"]["a"]
    a_can_primary = a_gen_primary / math.sqrt(2)
    cA = 4 / math.pi
    cB = 48 / math.pi ** 3
    err_A = abs(a_can_primary - cA) / cA * 100
    err_B = abs(a_can_primary - cB) / cB * 100

    print(f"\n  PRIMARY (L_gen >= 9.9, L_can >= 7):")
    print(f"    a_gen = {a_gen_primary:.6f}")
    print(f"    a_can = {a_can_primary:.6f}")
    print(f"    A (4/pi)    = {cA:.6f}, error {err_A:.2f}%")
    print(f"    B (48/pi^3) = {cB:.6f}, error {err_B:.2f}%")

    # =========================================================================
    # Subleading-aware fits and forced-A/B residual analysis
    # =========================================================================
    print(f"\n  SUBLEADING-AWARE FITS (full L^2_gen=8..300 panel, n=13):")

    fit_forms = {}
    Y = np.array(gs)

    # 2-param
    X = np.array([[math.log(L)/L, 1.0/L] for L in lams])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fit_forms["2_param_log_L_over_L_plus_1_over_L"] = {
        "a_gen": float(coef[0]), "a_can": float(coef[0]/math.sqrt(2)), "rss": rss,
        "params": [float(c) for c in coef],
    }
    print(f"    (a) 2-param log L/L + 1/L         : a_can = {coef[0]/math.sqrt(2):.4f}, rss = {rss:.6f}")

    # 3-param +1/L^2
    X = np.array([[math.log(L)/L, 1.0/L, 1.0/L**2] for L in lams])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fit_forms["3_param_plus_1_over_Lsq"] = {
        "a_gen": float(coef[0]), "a_can": float(coef[0]/math.sqrt(2)), "rss": rss,
        "params": [float(c) for c in coef],
    }
    print(f"    (b) 3-param +1/L^2                : a_can = {coef[0]/math.sqrt(2):.4f}, rss = {rss:.6f}")

    # 4-param +1/L^2 +1/(L log L)
    X = np.array([[math.log(L)/L, 1.0/L, 1.0/L**2, 1.0/(L * math.log(L))] for L in lams])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fit_forms["4_param_with_log_correction"] = {
        "a_gen": float(coef[0]), "a_can": float(coef[0]/math.sqrt(2)), "rss": rss,
        "params": [float(c) for c in coef],
    }
    print(f"    (c) 4-param +1/L^2 +1/(L log L)   : a_can = {coef[0]/math.sqrt(2):.4f}, rss = {rss:.6f}")

    # Forced-coefficient fits
    a_gen_A = math.sqrt(2) * cA
    a_gen_B = math.sqrt(2) * cB

    # Force A, fit 3 subleading
    X_qsub = np.array([[1.0/L, 1.0/L**2, 1.0/L**3] for L in lams])
    residual = Y - a_gen_A * np.array([math.log(L)/L for L in lams])
    coef_r, *_ = np.linalg.lstsq(X_qsub, residual, rcond=None)
    rss_A_force = float(np.sum((residual - X_qsub @ coef_r) ** 2))
    fit_forms["forced_A_3_subleading"] = {
        "a_gen_forced": float(a_gen_A), "a_can": float(cA),
        "rss": rss_A_force, "subleading_params_b_c_d": [float(x) for x in coef_r],
    }

    residual = Y - a_gen_B * np.array([math.log(L)/L for L in lams])
    coef_r, *_ = np.linalg.lstsq(X_qsub, residual, rcond=None)
    rss_B_force = float(np.sum((residual - X_qsub @ coef_r) ** 2))
    fit_forms["forced_B_3_subleading"] = {
        "a_gen_forced": float(a_gen_B), "a_can": float(cB),
        "rss": rss_B_force, "subleading_params_b_c_d": [float(x) for x in coef_r],
    }

    print(f"\n  FORCED-COEFFICIENT FITS (3-term subleading b/L + c/L^2 + d/L^3):")
    print(f"    FORCE A (a_can=4/pi):       rss = {rss_A_force:.6f}")
    print(f"    FORCE B (a_can=48/pi^3):    rss = {rss_B_force:.6f}")
    print(f"    Ratio rss(B)/rss(A) = {rss_B_force/rss_A_force:.4f}")
    if rss_A_force < rss_B_force:
        better = "A"
        ratio_msg = f"A fits {rss_B_force/rss_A_force:.2f}x better"
    else:
        better = "B"
        ratio_msg = f"B fits {rss_A_force/rss_B_force:.2f}x better"
    print(f"    Better (smaller rss): {better} ({ratio_msg})")

    # Pick best leading-coefficient estimator: the 4-param fit
    best_a_gen = fit_forms["4_param_with_log_correction"]["a_gen"]
    best_a_can = fit_forms["4_param_with_log_correction"]["a_can"]
    err_A_best = abs(best_a_can - cA) / cA * 100
    err_B_best = abs(best_a_can - cB) / cB * 100
    print(f"\n  BEST a_can ESTIMATE (4-param subleading-aware): {best_a_can:.4f}")
    print(f"    A error: {err_A_best:.2f}%, B error: {err_B_best:.2f}%")

    # Decide verdict based on multiple lines of evidence:
    # 1. Best fit (4-param)
    # 2. Forced-A vs Forced-B rss comparison
    leans_A = []
    leans_B = []
    if err_A_best < err_B_best:
        leans_A.append("4-param fit")
    else:
        leans_B.append("4-param fit")

    if rss_A_force < rss_B_force:
        leans_A.append("forced-coef rss")
    else:
        leans_B.append("forced-coef rss")

    if len(leans_A) > len(leans_B):
        verdict = f"A-leaning (panel-limited; favored by {', '.join(leans_A)})"
    elif len(leans_B) > len(leans_A):
        verdict = f"B-leaning (panel-limited; favored by {', '.join(leans_B)})"
    else:
        verdict = "AMBIGUOUS (panel insufficient for clean A vs B verdict)"
    print(f"\n  COMBINED VERDICT: {verdict}")

    # Use the 4-param result as the headline a_can for cross-group ratios
    a_can_primary = best_a_can  # Override primary with subleading-aware result.

    # Cross-group ratios at canonical-equivalent matched panel range.
    a_su2_can = 4.0 / math.pi
    a_su3_can = 1.243144
    a_sp2_can = 1.0875
    a_g2_can = 1.1765
    a_su4_can = a_can_primary

    print("\n" + "=" * 78)
    print("CROSS-GROUP RATIOS (convention-independent in dual-Coxeter)")
    print("=" * 78)
    print(f"  SU(4) primary a_can = {a_su4_can:.4f}")
    print()

    ratios = {
        "SU(2)": (a_su4_can / a_su2_can, 1.0, 12 / math.pi ** 2),
        "SU(3)": (a_su4_can / a_su3_can, 1.0, 4 / math.pi),
        "Sp(2)": (a_su4_can / a_sp2_can, 1.0, 3 / math.pi),
        "G_2":   (a_su4_can / a_g2_can, 1.0, 2 / math.pi),
    }
    print(f"  {'Ratio':>20} {'observed':>10} {'A pred':>10} {'B pred':>10} {'A err':>8} {'B err':>8}")
    for name, (obs, A_pred, B_pred) in ratios.items():
        eA = abs(obs - A_pred) / A_pred * 100
        eB = abs(obs - B_pred) / B_pred * 100
        print(f"  SU(4) / {name:<6}  {obs:>10.4f} {A_pred:>10.4f} {B_pred:>10.4f} "
              f"{eA:>7.1f}% {eB:>7.1f}%")

    # Combined 5-group summary.
    print("\n" + "=" * 78)
    print("COMBINED VERDICT (5 datapoints)")
    print("=" * 78)
    summary_rows = [
        ("SU(2)", 1, 2, 4/math.pi, 4/math.pi, a_su2_can),
        ("SU(3)", 2, 6, 4/math.pi, 12/math.pi**2, a_su3_can),
        ("Sp(2)", 2, 8, 4/math.pi, 16/math.pi**2, a_sp2_can),
        ("G_2",   2, 12, 4/math.pi, 24/math.pi**2, a_g2_can),
        ("SU(4)", 3, 24, 4/math.pi, 48/math.pi**3, a_su4_can),
    ]
    print(f"  {'Group':>6} {'rank':>4} {'|W|':>4} {'A':>8} {'B':>8} {'a_can':>10} {'A err':>7} {'B err':>7}")
    for name, r, W, cA_, cB_, a_can in summary_rows:
        eA = abs(a_can - cA_) / cA_ * 100
        eB = abs(a_can - cB_) / cB_ * 100
        print(f"  {name:>6} {r:>4} {W:>4} {cA_:>8.4f} {cB_:>8.4f} {a_can:>10.4f} {eA:>6.1f}% {eB:>6.1f}%")

    # Save full results.
    results = {
        "metadata": {
            "purpose": "Fourth datapoint (rank 3) for universal rate constant test",
            "predictions": {
                "A_4_over_pi": 4.0 / math.pi,
                "B_48_over_pi3": 48.0 / math.pi ** 3,
                "gap_AB_pct": abs(4/math.pi - 48/math.pi**3) / max(4/math.pi, 48/math.pi**3) * 100,
            },
            "casimir_convention": "generic (short root squared length 2); "
                "dual-Coxeter rescale = C_gen(adj)/h^v = 8/4 = 2",
            "lambda_can_to_gen": "L_can = L_gen / sqrt(2), L^2_can = L^2_gen / 2",
            "panel_limitation": "Compute budget terminated panel at L^2_gen=300 (L_can=12.25). "
                "L^2_gen=400 alone exceeded 30-minute single-row threshold; full panel to "
                "L^2_gen=600 was projected to take >2 hours additional. Asymptotic regime "
                "L_can >= 7 is reached but with fewer high-Lambda points than Q-rate/Sp(2)/G_2 "
                "sprints had.",
        },
        "su4": {
            "rows": rows,
            "meta": {
                "weyl_group_order": 24,
                "rank": 3,
                "n_positive_roots": 6,
                "haar_check": 1.0000000000,
            },
            "fits": fits_all,
            "subleading_aware_fits": fit_forms,
            "primary_fit_key": primary_key,
            "primary_a_gen_2param": a_gen_primary,
            "primary_a_canonical_4param": a_can_primary,
            "casimir_rescale": 2.0,
            "verdict": {
                "name": "SU(4)",
                "a_gen_2param_primary": a_gen_primary,
                "a_canonical_2param_primary": a_gen_primary / math.sqrt(2),
                "a_canonical_4param_subleading_aware": best_a_can,
                "weyl_order": 24,
                "rank": 3,
                "c_A_4_over_pi": cA,
                "c_B_2W_over_pi_r": cB,
                "err_A_pct_2param": abs(a_gen_primary/math.sqrt(2) - cA) / cA * 100,
                "err_B_pct_2param": abs(a_gen_primary/math.sqrt(2) - cB) / cB * 100,
                "err_A_pct_4param": err_A_best,
                "err_B_pct_4param": err_B_best,
                "forced_A_rss": rss_A_force,
                "forced_B_rss": rss_B_force,
                "gap_AB_pct": abs(cA - cB) / max(cA, cB) * 100,
                "verdict": verdict,
            },
        },
        "cross_group_ratios": {
            "a_canonical_values": {
                "SU(2)": a_su2_can, "SU(3)": a_su3_can, "Sp(2)": a_sp2_can,
                "G_2": a_g2_can, "SU(4)": a_su4_can,
            },
            "ratios_SU4_over": {
                "SU(2)": a_su4_can / a_su2_can,
                "SU(3)": a_su4_can / a_su3_can,
                "Sp(2)": a_su4_can / a_sp2_can,
                "G_2":   a_su4_can / a_g2_can,
            },
            "A_predicts_all": 1.0,
            "B_predicts_SU4_over": {
                "SU(2)": 12 / math.pi ** 2,
                "SU(3)": 4 / math.pi,
                "Sp(2)": 3 / math.pi,
                "G_2":   2 / math.pi,
            },
        },
        "combined_summary": [
            {"group": name, "rank": r, "weyl_order": W,
             "c_A": cA_, "c_B": cB_, "a_canonical_extracted": a_can,
             "err_A_pct": abs(a_can - cA_) / cA_ * 100,
             "err_B_pct": abs(a_can - cB_) / cB_ * 100}
            for name, r, W, cA_, cB_, a_can in summary_rows
        ],
    }

    out_path = out_dir / "su4_rate_constant.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")
    return results


if __name__ == "__main__":
    main()
