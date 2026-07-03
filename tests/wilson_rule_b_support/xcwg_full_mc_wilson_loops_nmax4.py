"""
XCWG-G extension to n_max=4 on Rule B (2026-05-16).

Background
==========
The XCWG-G sprint (debug/xcwg_full_mc_wilson_loops.py) measured the Wilson
loop string tension sigma(beta) on Rule B at n_max in {2, 3} via full
Metropolis-Hastings MC of the compact U(1) action. The ensemble-mean
per-area fit gave sigma_ens(beta) > 0 monotonically decreasing — the right
3D U(1) shape — BUT the joint sigma*A + mu*L fit and the fixed-perimeter
L=4 fit both showed the genuine area-law coefficient sigma_comb is
statistically zero. The apparent confinement signal was carried by the
perimeter term mu_comb. This was characterized as a finite-volume
perimeter-dominance artifact (Paper 41 v4, WEAK PASS with honest scope).

The XCWG-G memo recommended n_max=4 (V=60, E=312, beta_1=253) as the
cleanest discriminator. THAT IS THIS SCRIPT.

Question
========
At n_max=4 with V=60, E=312, does the joint sigma*A + mu*L fit show
statistically nonzero sigma_comb (area-law dominance)? Three possible
outcomes:
  (a) sigma_comb > 0 statistically at n_max=4: finite-volume perimeter
      dominance recedes; WEAK PASS strengthens toward CLEAN PASS.
  (b) sigma_comb still ~ 0 at n_max=4: perimeter dominance persists;
      need even larger graphs (or a different observable).
  (c) sigma_comb negative: structural surprise (would falsify the
      sixth witness as a confinement diagnostic).

Budget
======
n_max=4 has E=312 edges vs n_max=3's E=106 (~3x slower per sweep).
Cluster enumeration at A=5,6 is more expensive too.
Settings tuned for ~30-60 min single-agent runtime:
  - 3 beta values: {0.3, 1.0, 3.0}
  - n_therm = 600 sweeps
  - n_sample = 3000 samples, sample_interval = 10 sweeps
  - A_targets = [1, 2, 3, 4, 5, 6]
  - max_clusters_per_area = 30
  - n_loops_per_area = 8 (mild reduction from XCWG-G's 10)
  - plaq_cap = 1500 (a subset of the 5047 4-plaquettes; preserves the
    cluster-growth machinery without inflating the d_1 matrix)

Outputs
=======
  - tests/wilson_rule_b_support/data/xcwg_full_mc_wilson_loops_nmax4.json
  - tests/wilson_rule_b_support/data/xcwg_full_mc_wilson_loops_nmax4.png (if matplotlib)
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Reuse all the XCWG-G machinery
from xcwg_full_mc_wilson_loops import (  # noqa: E402
    run_at_nmax,
    polyakov_cross_check,
    fit_sigma_area_law,
    fit_sigma_perimeter_combined,
    fit_sigma_fixed_perimeter,
    fit_sigma_jackknife,
)

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    sys.stdout.reconfigure(line_buffering=True)


def main():
    out: Dict = {
        "sprint": "XCWG-G n_max=4 extension on Rule B compact U(1) (2026-05-16)",
        "motivation": (
            "XCWG-G at n_max in {2,3} found sigma_ens(beta) > 0 monotone (good "
            "shape) BUT joint sigma*A + mu*L fit gave sigma_comb statistically "
            "zero — apparent confinement signal was carried by perimeter term. "
            "Paper 41 v4 finite-volume scope. Here we extend to n_max=4 "
            "(V=60, E=312, beta_1=253) to test whether perimeter dominance "
            "recedes at the larger graph size."
        ),
        "question": (
            "Does the joint sigma_comb fit show statistically nonzero "
            "area-law coefficient at n_max=4?"
        ),
    }

    # XCWG-F Polyakov rho_M fit parameters (from XCWG-F memo, n_max=2)
    # Used for the Polyakov rate cross-check
    rho_M_A_xcwgF = 0.7471301473307033
    rho_M_c_xcwgF = 9.40152906810695
    out["xcwg_f_rho_M_fit_used_for_polyakov_cross_check"] = {
        "rho_M_A": rho_M_A_xcwgF,
        "rho_M_c": rho_M_c_xcwgF,
        "r_squared": 0.9815,
        "fit_beta_range": [0.05, 1.0],
        "note": (
            "Polyakov 3D compact U(1): sigma(beta) ~ exp(-(c_rho/2) beta) "
            "in the dilute regime; predicted slope c_sigma_predicted = "
            f"c_rho/2 = {0.5 * rho_M_c_xcwgF:.4f}."
        ),
    }

    # ---------- Budget-tuned settings for n_max=4 ----------
    # Beta grid focused on the area-vs-perimeter separation question.
    # 0.3: strong-coupling regime where confinement should be most visible
    # 1.0: intermediate
    # 3.0: dilute-monopole regime where Polyakov rate prediction is tested
    beta_grid_n4 = [0.3, 1.0, 3.0]
    out["beta_grid"] = beta_grid_n4

    # ---------- n_max = 4 run ----------
    print("\n" + "#" * 72)
    print("# n_max = 4 (V=60, E=312, beta_1=253) -- XCWG-G area-vs-perimeter test")
    print("#" * 72)
    t_start = time.time()
    try:
        out["n_max_4"] = run_at_nmax(
            n_max=4,
            beta_grid=beta_grid_n4,
            n_therm=600,
            n_sample=3000,
            sample_interval=10,
            n_loops_per_area=8,
            rng_seed=44,
            A_targets=[1, 2, 3, 4, 5, 6],
            max_clusters_per_area=30,
            plaq_cap=1500,
        )
    except Exception as e:
        out["n_max_4"] = {"error": repr(e)}
        print(f"  [n_max=4 FAILED]: {e}")
        import traceback
        traceback.print_exc()
    elapsed_n4 = time.time() - t_start
    print(f"\n  n_max=4 total elapsed: {elapsed_n4 / 60:.1f} min")
    out["n_max_4_wall_time_seconds"] = float(elapsed_n4)

    # ---------- Polyakov cross-check at n_max=4 ----------
    print("\n" + "=" * 72)
    print("POLYAKOV CROSS-CHECK (n_max = 4)")
    print("=" * 72)
    if "sigma_results" in out["n_max_4"]:
        pol4 = polyakov_cross_check(
            out["n_max_4"]["sigma_results"],
            rho_M_A=rho_M_A_xcwgF, rho_M_c=rho_M_c_xcwgF,
            beta_min_fit=0.3, beta_max_fit=3.0,
        )
        out["polyakov_cross_check_n_max_4"] = pol4
        print(f"  fit beta range: {pol4.get('beta_range')}")
        print(f"  n points fit: {pol4.get('n_points')}")
        if "c_sigma_MC" in pol4:
            print(f"  c_sigma_MC          = {pol4['c_sigma_MC']:.4f}")
            print(f"  c_sigma_predicted   = "
                  f"{pol4['c_sigma_predicted_polyakov_half_of_c_rho']:.4f}  "
                  f"(= c_rho/2 from XCWG-F)")
            print(f"  relative error      = "
                  f"{pol4['c_sigma_MC_vs_predicted_rel_error']:.4f}")
            print(f"  R^2 of exp fit      = {pol4['r_squared']:.4f}")
            print(f"  qualitative agree   = {pol4['qualitative_agreement']}")

    # ---------- Trend comparison across n_max = 2, 3, 4 ----------
    print("\n" + "=" * 72)
    print("TREND ACROSS n_max = 2, 3, 4 (area-law vs perimeter separation)")
    print("=" * 72)

    # n_max=2 reference values from XCWG-G memo (Table in §3)
    nmax2_table = [
        # (beta, sigma_ens, sigma_ens_se, sigma_comb, sigma_comb_se, mu_comb)
        (0.3, 0.174, 0.009, -0.006, 0.003, +0.217),
        (1.0, 0.033, 0.001, -0.002, 0.000, +0.034),
        (3.0, 0.011, 0.000, -0.000, 0.000, +0.010),
    ]
    # n_max=3 reference values from XCWG-G memo (Table in §6)
    nmax3_table = [
        (0.3, 0.131, 0.010, -0.024, 0.005, +0.178),
        (1.0, 0.021, 0.002, -0.016, 0.001, +0.048),
        # 3.0 not in original table; XCWG-G ran 0.1,0.3,0.5,1,2,5 at n_max=3
        # use beta=2 as proxy entry
        (2.0, 0.009, 0.001, -0.009, 0.001, +0.022),
    ]
    print(f"\n  At beta=0.3 (cleanest area-law probe regime):")
    print(f"  {'n_max':>6} {'sigma_ens':>13} {'sigma_comb':>16} "
          f"{'mu_comb':>10} {'sigma_comb_zsig':>16}")

    # n_max=2,3 from reference tables (literature)
    for n, tbl in [(2, nmax2_table), (3, nmax3_table)]:
        row = tbl[0]
        b, sens, sens_se, scom, scom_se, mu = row
        zsig = scom / scom_se if scom_se > 0 else 0
        print(f"  {n:>6} {sens:+7.4f}+/-{sens_se:.4f} "
              f"{scom:+7.4f}+/-{scom_se:.4f} {mu:+8.4f} {zsig:+12.2f}")

    nmax_trend = {
        "beta_grid_used": beta_grid_n4,
        "n_max_2_reference_xcwg_g": nmax2_table,
        "n_max_3_reference_xcwg_g": nmax3_table,
        "n_max_4_this_sprint": [],
    }

    if "sigma_results" in out["n_max_4"]:
        for r in out["n_max_4"]["sigma_results"]:
            b = r["beta"]
            sens = r.get("sigma_MC")
            sens_se = r.get("sigma_MC_se_wls", 0.0) or 0.0
            scom = r.get("sigma_MC_combined")
            scom_se = r.get("sigma_MC_combined_se", 0.0) or 0.0
            mu = r.get("mu_MC_combined")
            mu_se = r.get("mu_MC_combined_se", 0.0) or 0.0
            zsig = scom / scom_se if (scom is not None and scom_se > 0) else 0
            zmu = mu / mu_se if (mu is not None and mu_se > 0) else 0
            print(f"  {4:>6} "
                  f"{sens if sens is not None else 0:+7.4f}+/-{sens_se:.4f} "
                  f"{scom if scom is not None else 0:+7.4f}+/-{scom_se:.4f} "
                  f"{mu if mu is not None else 0:+8.4f} {zsig:+12.2f}"
                  f"  (b={b:.2f})")
            nmax_trend["n_max_4_this_sprint"].append({
                "beta": float(b),
                "sigma_ens": float(sens) if sens is not None else None,
                "sigma_ens_se": float(sens_se),
                "sigma_comb": float(scom) if scom is not None else None,
                "sigma_comb_se": float(scom_se),
                "mu_comb": float(mu) if mu is not None else None,
                "mu_comb_se": float(mu_se),
                "sigma_comb_z_significance": float(zsig),
                "mu_comb_z_significance": float(zmu),
            })

    out["trend_n_max_2_3_4"] = nmax_trend

    # ---------- Verdict ----------
    print("\n" + "=" * 72)
    print("VERDICT -- n_max=4 area-vs-perimeter separation")
    print("=" * 72)

    if "sigma_results" not in out["n_max_4"]:
        verdict = {"comment": "n_max=4 run failed; no verdict possible"}
        print("  n_max=4 run failed; cannot verdict")
    else:
        # Check sigma_comb sign and significance at each beta
        sig_comb_pos_signif: List[bool] = []
        sig_comb_neg_signif: List[bool] = []
        sig_comb_near_zero: List[bool] = []
        per_beta_status: List[Dict] = []
        for r in out["n_max_4"]["sigma_results"]:
            sc = r.get("sigma_MC_combined")
            scse = r.get("sigma_MC_combined_se", 0.0) or 0.0
            if sc is None:
                per_beta_status.append({
                    "beta": r["beta"], "status": "no_fit",
                })
                continue
            # >= 2 sigma from zero
            is_pos = (sc - 2 * scse) > 0
            is_neg = (sc + 2 * scse) < 0
            is_near_zero = (not is_pos) and (not is_neg)
            sig_comb_pos_signif.append(is_pos)
            sig_comb_neg_signif.append(is_neg)
            sig_comb_near_zero.append(is_near_zero)
            per_beta_status.append({
                "beta": float(r["beta"]),
                "sigma_comb": float(sc),
                "sigma_comb_se": float(scse),
                "z_score": float(sc / scse) if scse > 0 else 0.0,
                "status": "positive_2sigma" if is_pos
                          else ("negative_2sigma" if is_neg
                                else "near_zero_within_2sigma"),
            })

        # Aggregate verdict
        any_pos = any(sig_comb_pos_signif)
        all_pos = all(sig_comb_pos_signif) if sig_comb_pos_signif else False
        any_neg = any(sig_comb_neg_signif)
        all_near_zero = all(sig_comb_near_zero) if sig_comb_near_zero else False

        if all_pos:
            outcome = "AREA_LAW_DOMINANT_AT_NMAX4"
            text = ("sigma_comb is statistically positive (≥2sigma) at ALL "
                    "tested beta. Finite-volume perimeter dominance recedes; "
                    "WEAK PASS strengthens toward CLEAN PASS.")
        elif any_pos and not any_neg:
            outcome = "AREA_LAW_PARTIAL_AT_NMAX4"
            text = ("sigma_comb is statistically positive (≥2sigma) at SOME "
                    "but not all beta. Perimeter dominance partially recedes; "
                    "WEAK PASS improves with explicit regime-dependent scope.")
        elif all_near_zero:
            outcome = "PERIMETER_DOMINANT_PERSISTS_AT_NMAX4"
            text = ("sigma_comb still statistically zero at n_max=4. "
                    "Finite-volume perimeter dominance persists at this graph "
                    "size; need even larger graphs (or a different observable) "
                    "to cleanly separate area law.")
        elif any_neg:
            outcome = "NEGATIVE_SIGMA_COMB_STRUCTURAL_SURPRISE"
            text = ("sigma_comb is statistically NEGATIVE at some beta. This "
                    "is a structural surprise — area-law coefficient cannot "
                    "be intrinsically negative in a confining theory. Likely "
                    "indicates a non-trivial finite-volume correlation between "
                    "perimeter and area in the cluster ensemble that the "
                    "linear joint fit cannot disentangle.")
        else:
            outcome = "MIXED_INCONCLUSIVE"
            text = "Mixed results across beta; need more statistics or larger n_max."

        verdict = {
            "outcome": outcome,
            "verdict_text": text,
            "per_beta_status": per_beta_status,
        }
        print(f"\n  Outcome: {outcome}")
        print(f"  {text}")
        print(f"\n  Per-beta sigma_comb status:")
        for ps in per_beta_status:
            if "z_score" in ps:
                print(f"    beta={ps['beta']:.2f}: sigma_comb = "
                      f"{ps['sigma_comb']:+.4f} +/- {ps['sigma_comb_se']:.4f} "
                      f"(z = {ps['z_score']:+.2f}) -> {ps['status']}")
            else:
                print(f"    beta={ps['beta']:.2f}: {ps['status']}")

    out["verdict"] = verdict

    # ---------- Save ----------
    def _sanitize(obj):
        if isinstance(obj, dict):
            return {str(k): _sanitize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_sanitize(x) for x in obj]
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    out_path = os.path.join(_HERE, "data", "xcwg_full_mc_wilson_loops_nmax4.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(_sanitize(out), f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # ---------- Plot ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if "sigma_results" not in out["n_max_4"]:
            print("  [plot skipped]: n_max=4 run failed")
        else:
            fig, ax = plt.subplots(1, 3, figsize=(16, 5))

            # Panel 1: <W(A)> vs A on log-y axis at the 3 betas
            n4 = out["n_max_4"]
            cmap = plt.get_cmap("viridis")
            for i, r in enumerate(n4["beta_results"]):
                beta = r["beta"]
                A_arr = sorted(r["per_area"].keys())
                W_arr = [r["per_area"][A]["ens_mean"] for A in A_arr]
                W_se = [r["per_area"][A]["ens_se"] for A in A_arr]
                color = cmap(i / max(len(n4["beta_results"]) - 1, 1))
                ax[0].errorbar(A_arr, W_arr, yerr=W_se, marker="o", linestyle="-",
                               color=color, label=f"β={beta:.2g}", capsize=2)
            ax[0].set_yscale("symlog", linthresh=1e-3)
            ax[0].set_xlabel(r"Area $A$ (plaquettes)")
            ax[0].set_ylabel(r"$\langle W(C) \rangle_\beta$")
            ax[0].set_title("Wilson loop ⟨W⟩ vs area (n_max=4)")
            ax[0].grid(True, alpha=0.3)
            ax[0].legend(loc="best", fontsize=8)
            ax[0].axhline(0, color="k", linewidth=0.5, alpha=0.3)

            # Panel 2: sigma_ens vs sigma_comb at each beta
            betas = [r["beta"] for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            s_ens = [r["sigma_MC"] for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            s_ens_se = [r.get("sigma_MC_se_wls", 0.0) or 0.0
                        for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            s_comb = [r.get("sigma_MC_combined") or 0.0
                      for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            s_comb_se = [r.get("sigma_MC_combined_se") or 0.0
                         for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            mu_comb = [r.get("mu_MC_combined") or 0.0
                       for r in n4["sigma_results"] if r.get("sigma_MC") is not None]
            ax[1].errorbar(betas, s_ens, yerr=s_ens_se, marker="o", linestyle="-",
                           label=r"$\sigma_{\rm ens}$ (perim-contaminated)", capsize=3)
            ax[1].errorbar(betas, s_comb, yerr=s_comb_se, marker="s", linestyle="-",
                           label=r"$\sigma_{\rm comb}$ (area-only)", capsize=3)
            ax[1].plot(betas, mu_comb, marker="d", linestyle="--",
                       label=r"$\mu_{\rm comb}$ (perim coeff)")
            ax[1].set_xlabel(r"$\beta$")
            ax[1].set_ylabel(r"$\sigma$ or $\mu$")
            ax[1].set_title("Area vs perimeter at n_max=4")
            ax[1].grid(True, alpha=0.3)
            ax[1].axhline(0, color="k", linewidth=0.5, alpha=0.3)
            ax[1].legend(loc="best", fontsize=8)

            # Panel 3: trend n_max=2,3,4 at beta=0.3 (or smallest common beta)
            n_max_list = [2, 3, 4]
            sigma_ens_at_03 = [0.174, 0.131, None]
            sigma_ens_se_at_03 = [0.009, 0.010, None]
            sigma_comb_at_03 = [-0.006, -0.024, None]
            sigma_comb_se_at_03 = [0.003, 0.005, None]

            # Fill n_max=4 value from data
            for r in n4["sigma_results"]:
                if abs(r["beta"] - 0.3) < 0.01:
                    sigma_ens_at_03[2] = r.get("sigma_MC")
                    sigma_ens_se_at_03[2] = r.get("sigma_MC_se_wls", 0.0) or 0.0
                    sigma_comb_at_03[2] = r.get("sigma_MC_combined")
                    sigma_comb_se_at_03[2] = r.get("sigma_MC_combined_se", 0.0) or 0.0
                    break

            x_arr = n_max_list
            y_ens = sigma_ens_at_03
            ye_ens = sigma_ens_se_at_03
            y_com = sigma_comb_at_03
            ye_com = sigma_comb_se_at_03
            if all(v is not None for v in y_ens):
                ax[2].errorbar(x_arr, y_ens, yerr=ye_ens, marker="o", linestyle="-",
                               label=r"$\sigma_{\rm ens}$ (perim-contaminated)", capsize=3)
            if all(v is not None for v in y_com):
                ax[2].errorbar(x_arr, y_com, yerr=ye_com, marker="s", linestyle="-",
                               label=r"$\sigma_{\rm comb}$ (area-only)", capsize=3)
            ax[2].axhline(0, color="k", linewidth=0.5, alpha=0.3)
            ax[2].set_xticks(n_max_list)
            ax[2].set_xlabel(r"$n_{\rm max}$")
            ax[2].set_ylabel(r"$\sigma$ at $\beta=0.3$")
            ax[2].set_title("Area-vs-perimeter trend with n_max")
            ax[2].grid(True, alpha=0.3)
            ax[2].legend(loc="best", fontsize=8)

            fig.tight_layout()
            plot_path = os.path.join(_HERE, "data", "xcwg_full_mc_wilson_loops_nmax4.png")
            os.makedirs(os.path.dirname(plot_path), exist_ok=True)
            fig.savefig(plot_path, dpi=120)
            plt.close(fig)
            print(f"Wrote {plot_path}")
    except Exception as e:
        print(f"  [plot skipped]: {e}")


if __name__ == "__main__":
    main()
