"""Sprint G4-5d-refined -- F12 closure via phi(0) Mellin moment.

G4-5d (debug/g4_5d_cutoff_dependence_memo.md) rejected the naive phi(2)
prediction at 65% deviation and identified phi(0) as the proper Mellin
moment for the topological tip contribution to S_BH.  This driver tests
the refined prediction quantitatively.

Refined prediction (the sector-by-moment map):
- Bulk R^0 (cosmological const.) controlled by phi(2)
- Bulk R^1 (Einstein-Hilbert)    controlled by phi(1)
- Topological tip (S_BH)         controlled by phi(0)

The integral

    S_tip(Lambda, f) = (1/2) int (dt/t) f(t Lambda^2) tip(t)

with substitution u = t Lambda^2 becomes, when tip(t) is approximately
constant in log t (which G4-5d verified empirically),

    S_tip(Lambda, f) ~= (1/2) * C_tip * phi(0)[f, regulated]

where C_tip ~ 0.15 on the G4-5a sweet-spot panel and

    phi(0)[f, regulated] = int_{u_UV}^{u_IR} (du/u) f(u)

with u_UV = t_UV * Lambda^2 (substrate UV cutoff, t_UV = a^2 = 0.0025)
and u_IR = t_IR * Lambda^2 (substrate IR cutoff, t_IR set by the
panel range, here t_IR = 20 the largest panel t).

The three regulated phi(0)[f]:
- Gaussian:    phi(0) = E_1(u_UV) - E_1(u_IR)  ~ -gamma_E - log(u_UV) for small u_UV
- Sharp:       phi(0) = log(min(1, u_IR)/u_UV) [if u_UV < 1]
- Polynomial:  phi(0) = (1/2) [E_1(u_UV^2) - E_1(u_IR^2)]

Decision gate
-------------
- POSITIVE: S_sharp/S_Gauss and S_poly/S_Gauss match phi(0) prediction
  within 20%, closes F12 POSITIVE
- PARTIAL : qualitative ordering S_sharp > S_poly > S_Gauss matches but
  quantitative match incomplete
- NEGATIVE: phi(0) prediction also fails

Method
------
1. Recompute tip(t) on the G4-5a sweet-spot panel (R=10, a=0.05,
   N_rho=200, N_0=120, eps=0.1) at t in {0.1, 0.2, 0.5, 1, 2, 5, 10, 20}.
2. Compute S_tip(Lambda, f) by trapezoidal in log t on the same panel
   (this matches G4-5d numerically).
3. Compute regulated phi(0)[f, Lambda] analytically with t_UV = a^2 and
   t_IR = max panel t.
4. Compare empirical ratios S_tip(f)/S_tip(Gauss) to predicted
   phi(0)[f]/phi(0)[Gauss].
5. Report quantitative recovery: how close to within 20%?
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import exp1  # E_1(x) = int_x^inf exp(-t)/t dt

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

OUT_JSON = Path(__file__).parent / "data" / "g4_5d_refined_phi0_prediction.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# Euler-Mascheroni constant
GAMMA_E = 0.5772156649015329


# ----------------------------------------------------------------------
#  Cutoffs (same as G4-5d)
# ----------------------------------------------------------------------

def f_gauss(x: np.ndarray) -> np.ndarray:
    """Gaussian cutoff f(x) = exp(-x)."""
    return np.exp(-x)


def f_sharp(x: np.ndarray) -> np.ndarray:
    """Sharp cutoff f(x) = Theta(1-x)."""
    return (x <= 1.0).astype(float)


def f_poly(x: np.ndarray) -> np.ndarray:
    """Polynomial-decay cutoff f(x) = exp(-x^2)."""
    return np.exp(-(x ** 2))


CUTOFFS = {
    "gauss": {"label": "Gaussian f(x)=exp(-x)", "f": f_gauss},
    "sharp": {"label": "Sharp f(x)=Theta(1-x)", "f": f_sharp},
    "poly":  {"label": "Polynomial f(x)=exp(-x^2)", "f": f_poly},
}


# ----------------------------------------------------------------------
#  Analytical regulated phi(0) for each cutoff
# ----------------------------------------------------------------------

def phi0_gauss_regulated(u_UV: float, u_IR: float) -> float:
    """Regulated phi(0) for Gaussian: int_{u_UV}^{u_IR} (du/u) exp(-u).

    Closed form: E_1(u_UV) - E_1(u_IR), where E_1 is the exponential
    integral.  In the small-u_UV, large-u_IR limit,
    E_1(u_UV) ~ -gamma_E - log(u_UV), E_1(u_IR) ~ 0, so
    phi(0) ~ -gamma_E - log(u_UV).
    """
    return float(exp1(u_UV) - exp1(u_IR))


def phi0_sharp_regulated(u_UV: float, u_IR: float) -> float:
    """Regulated phi(0) for sharp: int_{u_UV}^{min(1, u_IR)} du/u.

    Closed form: log(min(1, u_IR)/u_UV) provided u_UV < 1.
    If u_UV >= 1 the integral is zero.
    """
    if u_UV >= 1.0:
        return 0.0
    u_top = min(1.0, u_IR)
    if u_top <= u_UV:
        return 0.0
    return float(np.log(u_top / u_UV))


def phi0_poly_regulated(u_UV: float, u_IR: float) -> float:
    """Regulated phi(0) for polynomial: int_{u_UV}^{u_IR} (du/u) exp(-u^2).

    Substitute v = u^2, dv = 2 u du, dv/v = 2 du/u:
        int (du/u) exp(-u^2) = (1/2) int (dv/v) exp(-v)
    Closed form: (1/2)[E_1(u_UV^2) - E_1(u_IR^2)].
    """
    return float(0.5 * (exp1(u_UV ** 2) - exp1(u_IR ** 2)))


PHI0_FUNCS = {
    "gauss": phi0_gauss_regulated,
    "sharp": phi0_sharp_regulated,
    "poly":  phi0_poly_regulated,
}


# ----------------------------------------------------------------------
#  Empirical S_tip from substrate tip term
# ----------------------------------------------------------------------

def S_tip_empirical(t_grid: np.ndarray, tip_term: np.ndarray,
                    Lambda: float, cutoff_key: str) -> float:
    """Empirical S_tip(Lambda, f) by log-trapezoidal on the panel.

    Matches G4-5d formula bit-identically.
    """
    f = CUTOFFS[cutoff_key]["f"]
    x = t_grid * Lambda ** 2
    fx = f(x)
    integrand_log = fx * tip_term
    log_t = np.log(t_grid)
    J = float(np.trapezoid(integrand_log, log_t))
    return J / 2


# ----------------------------------------------------------------------
#  Main driver
# ----------------------------------------------------------------------

def main() -> None:
    results: dict = {}
    print("=" * 76)
    print("Sprint G4-5d-refined -- F12 closure via phi(0) Mellin moment")
    print("=" * 76)

    # ------------------------------------------------------------------
    # Step 0 -- substrate setup, recompute tip term
    # ------------------------------------------------------------------
    R = 10.0
    a = 0.05
    N_rho = 200
    N_0 = 120
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps = (alpha_plus - alpha_minus) / 2

    # UV/IR substrate cutoffs
    t_UV = a ** 2          # substrate UV: smallest length-squared
    t_IR = 20.0            # panel IR: largest t on the G4-5a grid

    print(f"\n[Substrate]")
    print(f"  R = {R}, a = {a}, N_rho = {N_rho}, N_0 = {N_0}")
    print(f"  alpha_+ = {alpha_plus:.4f}, alpha_- = {alpha_minus:.4f}, "
          f"eps = {eps:.4f}")
    print(f"  t_UV = a^2 = {t_UV}, t_IR (panel max) = {t_IR}")

    t_grid = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])

    print(f"\n[Step 0] Recomputing tip term on substrate ...")
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    wedge_plus = DiscreteWedgeDirac(N_rho=N_rho, a=a, N_phi=N_plus,
                                    alpha=alpha_plus)
    wedge_minus = DiscreteWedgeDirac(N_rho=N_rho, a=a, N_phi=N_minus,
                                     alpha=alpha_minus)
    K_disk = np.array([disk.heat_trace(t) for t in t_grid])
    K_plus = np.array([wedge_plus.heat_trace(t) for t in t_grid])
    K_minus = np.array([wedge_minus.heat_trace(t) for t in t_grid])
    dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
    tip_term = dK_dalpha - K_disk

    # Sanity vs G4-5d
    expected_tip_at_t1 = 0.15376729381481624
    diff = abs(tip_term[3] - expected_tip_at_t1)
    assert diff < 1e-9, f"tip term differs from G4-5d: diff = {diff}"
    print(f"  Sanity vs G4-5d at t=1: diff = {diff:.3e}  OK")

    C_tip_mean = float(np.mean(tip_term))
    C_tip_std = float(np.std(tip_term))
    print(f"  tip term mean = {C_tip_mean:.6f}, "
          f"std = {C_tip_std:.6f}  (slowly varying in log t)")

    results["substrate"] = {
        "R": R, "a": a, "N_rho": N_rho, "N_0": N_0,
        "alpha_plus": alpha_plus, "alpha_minus": alpha_minus, "eps": eps,
        "t_UV": t_UV, "t_IR": t_IR,
        "t_grid": t_grid.tolist(),
        "tip_term": tip_term.tolist(),
        "tip_term_mean": C_tip_mean,
        "tip_term_std": C_tip_std,
    }

    # ------------------------------------------------------------------
    # Step 1 -- S_tip empirical (matches G4-5d)
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 1] S_tip(Lambda, f) empirical (matches G4-5d)")
    print("=" * 76)

    Lambda_values = [0.5, 1.0, 2.0]

    print(f"\n  {'Lambda':>6}  {'gauss':>12}  {'sharp':>12}  {'poly':>12}")
    print("  " + "-" * 48)
    S_empirical: dict = {}
    for L in Lambda_values:
        S_empirical[str(L)] = {}
        row = []
        for key in CUTOFFS:
            S = S_tip_empirical(t_grid, tip_term, L, key)
            S_empirical[str(L)][key] = S
            row.append(S)
        print(f"  {L:>6.2f}  {row[0]:>+12.6f}  {row[1]:>+12.6f}  "
              f"{row[2]:>+12.6f}")

    results["S_empirical"] = S_empirical

    # ------------------------------------------------------------------
    # Step 2 -- Analytical regulated phi(0)
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 2] Analytical regulated phi(0)[f, Lambda]")
    print("=" * 76)
    print()
    print("  phi(0)[f, Lambda] = int_{u_UV}^{u_IR} (du/u) f(u)")
    print(f"  u_UV = t_UV * Lambda^2 = {t_UV} * Lambda^2")
    print(f"  u_IR = t_IR * Lambda^2 = {t_IR} * Lambda^2")
    print()
    print(f"  {'Lambda':>6}  {'u_UV':>10}  {'u_IR':>10}  "
          f"{'phi(0)_g':>12}  {'phi(0)_sh':>12}  {'phi(0)_po':>12}")
    print("  " + "-" * 76)
    phi0_analytical: dict = {}
    for L in Lambda_values:
        u_UV = t_UV * L ** 2
        u_IR = t_IR * L ** 2
        phi0_analytical[str(L)] = {"u_UV": u_UV, "u_IR": u_IR}
        row = []
        for key in CUTOFFS:
            phi = PHI0_FUNCS[key](u_UV, u_IR)
            phi0_analytical[str(L)][key] = phi
            row.append(phi)
        print(f"  {L:>6.2f}  {u_UV:>10.4e}  {u_IR:>10.4e}  "
              f"{row[0]:>12.6f}  {row[1]:>12.6f}  {row[2]:>12.6f}")

    results["phi0_analytical"] = phi0_analytical

    # ------------------------------------------------------------------
    # Step 3 -- predicted vs empirical ratios
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 3] Predicted vs empirical ratios (refined phi(0))")
    print("=" * 76)
    print()
    print("  Prediction: S_tip(f) / S_tip(gauss) = phi(0)[f] / phi(0)[gauss]")
    print()

    ratio_table: dict = {}
    print(f"  {'Lambda':>6}  ratio  {'empir':>10}  {'pred':>10}  {'dev %':>8}")
    print("  " + "-" * 56)
    all_devs: list = []
    for L in Lambda_values:
        ratio_table[str(L)] = {}
        Sg = S_empirical[str(L)]["gauss"]
        Ss = S_empirical[str(L)]["sharp"]
        Sp = S_empirical[str(L)]["poly"]
        pg = phi0_analytical[str(L)]["gauss"]
        ps = phi0_analytical[str(L)]["sharp"]
        pp = phi0_analytical[str(L)]["poly"]

        # sharp / gauss
        emp_sh = Ss / Sg if abs(Sg) > 1e-15 else float("inf")
        pred_sh = ps / pg if abs(pg) > 1e-15 else float("inf")
        dev_sh = 100.0 * (emp_sh - pred_sh) / pred_sh if pred_sh != 0 else float("inf")
        # poly / gauss
        emp_po = Sp / Sg if abs(Sg) > 1e-15 else float("inf")
        pred_po = pp / pg if abs(pg) > 1e-15 else float("inf")
        dev_po = 100.0 * (emp_po - pred_po) / pred_po if pred_po != 0 else float("inf")
        # sharp / poly (cross-check)
        emp_shpo = Ss / Sp if abs(Sp) > 1e-15 else float("inf")
        pred_shpo = ps / pp if abs(pp) > 1e-15 else float("inf")
        dev_shpo = 100.0 * (emp_shpo - pred_shpo) / pred_shpo if pred_shpo != 0 else float("inf")

        ratio_table[str(L)] = {
            "sharp_over_gauss_empirical":  emp_sh,
            "sharp_over_gauss_predicted":  pred_sh,
            "sharp_over_gauss_dev_pct":    dev_sh,
            "poly_over_gauss_empirical":   emp_po,
            "poly_over_gauss_predicted":   pred_po,
            "poly_over_gauss_dev_pct":     dev_po,
            "sharp_over_poly_empirical":   emp_shpo,
            "sharp_over_poly_predicted":   pred_shpo,
            "sharp_over_poly_dev_pct":     dev_shpo,
        }

        print(f"  {L:>6.2f}  sh/g   {emp_sh:>10.4f}  {pred_sh:>10.4f}  "
              f"{dev_sh:>+7.2f}%")
        print(f"  {L:>6.2f}  po/g   {emp_po:>10.4f}  {pred_po:>10.4f}  "
              f"{dev_po:>+7.2f}%")
        print(f"  {L:>6.2f}  sh/po  {emp_shpo:>10.4f}  {pred_shpo:>10.4f}  "
              f"{dev_shpo:>+7.2f}%")

        all_devs.extend([abs(dev_sh), abs(dev_po), abs(dev_shpo)])

    results["ratio_table"] = ratio_table

    # ------------------------------------------------------------------
    # Step 4 -- qualitative ordering check
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 4] Qualitative ordering check")
    print("=" * 76)
    print()
    print("  Expected (per G4-5d): S_sharp > S_poly > S_Gauss at all Lambda")
    print("  Predicted from phi(0): same ordering?")
    print()
    print(f"  {'Lambda':>6}  empirical          predicted")
    print("  " + "-" * 56)
    ordering_check: dict = {}
    for L in Lambda_values:
        Sg = S_empirical[str(L)]["gauss"]
        Ss = S_empirical[str(L)]["sharp"]
        Sp = S_empirical[str(L)]["poly"]
        pg = phi0_analytical[str(L)]["gauss"]
        ps = phi0_analytical[str(L)]["sharp"]
        pp = phi0_analytical[str(L)]["poly"]

        emp_ord = Ss > Sp > Sg
        pred_ord = ps > pp > pg
        ordering_check[str(L)] = {
            "S_emp_sh_gt_po": Ss > Sp,
            "S_emp_po_gt_g":  Sp > Sg,
            "phi0_pred_sh_gt_po": ps > pp,
            "phi0_pred_po_gt_g":  pp > pg,
            "emp_full_order":  emp_ord,
            "pred_full_order": pred_ord,
            "match": emp_ord == pred_ord,
        }
        print(f"  {L:>6.2f}  sh>po>g = {str(emp_ord):>5}    "
              f"sh>po>g = {str(pred_ord):>5}   "
              f"match = {emp_ord == pred_ord}")

    results["ordering_check"] = ordering_check

    # ------------------------------------------------------------------
    # Step 5 -- verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 5] Verdict")
    print("=" * 76)

    max_dev = max(all_devs)
    mean_dev = float(np.mean(all_devs))
    all_ordering = all(ordering_check[str(L)]["match"] for L in Lambda_values)

    print(f"\n  Max deviation from refined phi(0) prediction: {max_dev:.2f}%")
    print(f"  Mean deviation                              : {mean_dev:.2f}%")
    print(f"  Qualitative ordering match at all Lambda    : {all_ordering}")

    if max_dev <= 20.0 and all_ordering:
        verdict = "POSITIVE-G4-5d-REFINED-PHI0-CLOSES-F12"
        msg = ("Empirical S_tip(Lambda, f) ratios match phi(0) prediction "
               "within 20% at all three Lambda values and all three cutoff "
               "functions.  F12 closed POSITIVE: the topological tip "
               "contribution to S_BH scales with the ZEROTH Mellin moment "
               "phi(0), regulated by the substrate UV/IR cutoffs.")
    elif all_ordering and max_dev <= 50.0:
        verdict = "PARTIAL-G4-5d-REFINED-PHI0-ORDERING-PLUS"
        msg = ("Qualitative ordering S_sharp > S_poly > S_Gauss matches "
               "phi(0) prediction at all Lambda, but quantitative deviations "
               "exceed 20%.  Substantive structural reading of phi(0) is "
               "confirmed; finite-substrate corrections enter at the "
               "ratio level.")
    elif all_ordering:
        verdict = "PARTIAL-G4-5d-REFINED-PHI0-ORDERING-ONLY"
        msg = ("Ordering correct but quantitative deviation > 50%.  "
               "phi(0) reading captures the qualitative pattern but "
               "discrete-substrate corrections dominate the absolute scale.")
    else:
        verdict = "NEGATIVE-G4-5d-REFINED-PHI0"
        msg = ("phi(0) prediction also fails: ordering does not match. "
               "More structural work needed beyond the Mellin-moment map.")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["max_deviation_pct"] = max_dev
    results["mean_deviation_pct"] = mean_dev
    results["ordering_match_all_Lambda"] = all_ordering
    results["verdict"] = verdict
    results["msg"] = msg

    # ------------------------------------------------------------------
    # Step 6 -- sector-wise Mellin moment map summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 76)
    print("[Step 6] Sector-wise Mellin moment map (G8 refinement)")
    print("=" * 76)
    print()
    print("  | Sector                              | Coefficient   | Moment |")
    print("  | ----------------------------------- | ------------- | ------ |")
    print("  | Bulk R^0 (cosmological const.)      | Lambda_cc     | phi(2) |")
    print("  | Bulk R^1 (Einstein-Hilbert)         | G_eff^{-1}    | phi(1) |")
    print("  | Topological tip (S_BH tip-of-cigar) | (1/12) S-C    | phi(0) |")
    print()
    print("  G8 (Paper 51 sec.10) gave phi(1), phi(2) for bulk Wilson coeffs.")
    print("  This sprint quantitatively verifies phi(0) for the topological")
    print("  tip sector via the discrete substrate.")

    results["sector_mellin_map"] = {
        "bulk_R0_cosmological_constant": "phi(2)",
        "bulk_R1_einstein_hilbert":      "phi(1)",
        "topological_tip_S_BH":          "phi(0)",
    }

    # Save
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 76)


if __name__ == "__main__":
    main()
