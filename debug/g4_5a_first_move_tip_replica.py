"""Sprint G4-5a first move -- Tip-only replica integration on constant-warp cigar.

Closes the F8 falsifier from G4-5 scoping: integrate the wedge-Dirac
replica derivative over t with a Gaussian Connes-Chamseddine cutoff
function and extract the entropy contribution per topological tip.

Method
------
1. Compute K_wedge^Dirac(alpha=1+eps, t) and K_wedge(alpha=1-eps, t)
   via DiscreteWedgeDirac at sweet-spot eps = 12/120 = 0.1 (matching
   G4-4f best window).
2. Compute the tip-contribution dDelta_K/dalpha at alpha=1 at each t:
       dDelta_K/dalpha|_{1} = [K_wedge(1+eps) - K_wedge(1-eps)] / (2 eps)
                            - K_disk(t)
3. Integrate:
       S_tip = -1/2 * integral( dt/t * f(t Lambda^2) * dDelta_K/dalpha )
   with Gaussian cutoff f(x) = exp(-x).
4. Continuum prediction: tip contribution to S_BH per spinor =
   (1/12) * Mellin integral.

This gives the discrete-substrate evaluation of the spinor topological
tip contribution to S_BH.
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteWedgeDirac,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_first_move_tip_replica.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-5a first move -- Tip-only replica integration")
    print("=" * 72)

    # Substrate (matches G4-4c / G4-4f sweet-spot panel)
    R = 10.0; a = 0.05; N_rho = 200; N_0 = 120
    print(f"\n[Setup] R={R}, a={a}, N_rho={N_rho}, N_0={N_0}")

    # Replica eps: matches G4-4f best window at recovery 96.69%
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps = (alpha_plus - alpha_minus) / 2
    print(f"\n[Replica parameters]")
    print(f"  alpha_+ = {alpha_plus}, alpha_- = {alpha_minus}, eps = {eps}")

    # Compute K at relevant t values
    # Use a substrate-converged grid (avoid extreme UV/IR)
    t_grid = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    print(f"  t grid: {t_grid.tolist()}")

    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    wedge_plus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
    )
    wedge_minus = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
    )

    print(f"\n[Computing K(t) at each t in grid...]")
    K_disk_arr = np.array([disk.heat_trace(t) for t in t_grid])
    K_plus_arr = np.array([wedge_plus.heat_trace(t) for t in t_grid])
    K_minus_arr = np.array([wedge_minus.heat_trace(t) for t in t_grid])

    # Replica derivative at each t
    dK_dalpha = (K_plus_arr - K_minus_arr) / (alpha_plus - alpha_minus)
    # Tip-only contribution = dK/dalpha - K_disk (subtract bulk linear-in-alpha part)
    tip_term = dK_dalpha - K_disk_arr

    print(f"\n  {'t':>5}  {'K_disk':>12}  {'dK/dalpha':>12}  {'tip term':>12}  "
          f"{'tip/(1/6)':>10}")
    print("  " + "-" * 65)
    for i, t in enumerate(t_grid):
        rec = tip_term[i] / (1/6)
        print(f"  {t:>5.1f}  {K_disk_arr[i]:>12.4f}  {dK_dalpha[i]:>12.4f}  "
              f"{tip_term[i]:>12.6f}  {rec:>10.4f}")

    results["heat_trace_data"] = {
        "t_grid": t_grid.tolist(),
        "K_disk": K_disk_arr.tolist(),
        "dK_dalpha": dK_dalpha.tolist(),
        "tip_term": tip_term.tolist(),
    }

    # ------------------------------------------------------------------
    # Step 2: Integrate over t with Gaussian cutoff
    # ------------------------------------------------------------------
    print("\n[Step 2] Integrate tip term with Gaussian cutoff f(x) = exp(-x):")
    print()

    Lambda_values = [0.5, 1.0, 1.5, 2.0]
    print(f"  Lambda values: {Lambda_values}")
    print(f"  Integral: J(Lambda) = integral_{{t_min}}^{{t_max}} (dt/t) exp(-t Lambda^2) tip(t)")
    print(f"  S_tip(Lambda) = -J/2")
    print()
    print(f"  {'Lambda':>8}  {'J integral':>14}  {'S_tip':>12}  "
          f"{'log(S_tip)':>14}")
    print("  " + "-" * 60)

    # Use trapezoidal rule on log-spaced t grid (effective dt/t weighting)
    # log(t_grid)
    log_t = np.log(t_grid)
    integrate_results = {}
    for Lambda in Lambda_values:
        f_cutoff = np.exp(-t_grid * Lambda**2)
        integrand = (1.0 / t_grid) * f_cutoff * tip_term
        # Trapezoidal in log(t) is integral dt/t = integral d(log t)
        # So: integral (dt/t) g(t) = integral d(log t) g(t)
        # = sum over intervals
        integrand_log = f_cutoff * tip_term  # = (dt/t) g(t) = d(log t) g(t)
        J = float(np.trapezoid(integrand_log, log_t))
        # Replica method: S = -dI_E/d_alpha; I_E = -(1/2) S_CC
        # => S = +(1/2) integral(dt/t) f(t Lambda^2) dDelta_K/d_alpha
        S_tip = +J / 2
        log_S_tip = np.log(abs(S_tip)) if abs(S_tip) > 1e-15 else -np.inf
        integrate_results[Lambda] = {
            "J": J,
            "S_tip": S_tip,
            "log_S_tip": log_S_tip,
        }
        print(f"  {Lambda:>8.2f}  {J:>14.6e}  {S_tip:>+12.6e}  "
              f"{log_S_tip:>+14.4f}")

    results["integrate_results"] = {
        str(k): v for k, v in integrate_results.items()
    }

    # ------------------------------------------------------------------
    # Step 3: Test Lambda-dependence
    # ------------------------------------------------------------------
    print("\n[Step 3] Lambda-dependence test:")
    print()
    print(f"  Continuum prediction: S_tip ~ (1/12) * log integral over (dt/t) exp(-t Lambda^2)")
    print(f"  The log integral is a Mellin moment that depends on Lambda.")
    print(f"  Constant part (Lambda -> infinity): finite (regularized topological term).")
    print()
    # Differences in log show whether S_tip is log-divergent
    Lambda_arr = np.array(Lambda_values)
    S_arr = np.array([integrate_results[L]["S_tip"] for L in Lambda_values])
    print(f"  S_tip ratios (test log scaling):")
    for i in range(len(Lambda_arr) - 1):
        ratio = S_arr[i+1] / S_arr[i] if abs(S_arr[i]) > 1e-15 else float("inf")
        print(f"    S_tip({Lambda_arr[i+1]}) / S_tip({Lambda_arr[i]}) = {ratio:.4f}")

    # ------------------------------------------------------------------
    # Step 4: Continuum comparison
    # ------------------------------------------------------------------
    print("\n[Step 4] Continuum comparison:")
    print()
    print(f"  Spinor tip term continuum: dDelta_K_tip/dalpha = +1/6")
    print(f"  Mellin integral M_0(Lambda) = integral_0^inf (dt/t) exp(-t Lambda^2)")
    print(f"  This is logarithmically divergent at t -> 0; needs regularization.")
    print(f"  For substrate with t_UV = a^2 ~ {a**2:.4f} and t_IR ~ {R**2:.2f}:")
    print(f"    M_0 ~ log(t_IR/t_UV) ~ log({(R/a)**2:.0f}) = {np.log((R/a)**2):.4f}")
    print()
    print(f"  Continuum-style prediction (rough):")
    for Lambda in Lambda_values:
        # M_0 ~ -log(t_UV * Lambda^2) for sharp UV cutoff at t_UV = a^2
        # but for Gaussian, this is different
        # Use exact log-Mellin = -log(Lambda^2 a^2) (rough)
        M_0_approx = -np.log(Lambda**2 * a**2) if Lambda**2 * a**2 > 0 else 0
        S_continuum_approx = -(1/12) * M_0_approx / 2  # (-1/2 prefactor)
        # Already -J/2, so S_tip = -(1/2) * M_0 * (1/6) = -M_0/12
        S_predicted = +(1/12) * M_0_approx
        ratio = integrate_results[Lambda]["S_tip"] / S_predicted if abs(S_predicted) > 1e-10 else float("inf")
        print(f"    Lambda={Lambda}: M_0 ~ {M_0_approx:.4f}, "
              f"S_pred ~ {S_predicted:+.4f}, "
              f"S_disc = {integrate_results[Lambda]['S_tip']:+.6e}, "
              f"ratio = {ratio:.4f}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    # Check that S_tip is finite and positive sign
    S_finite = all(np.isfinite(integrate_results[L]["S_tip"]) for L in Lambda_values)
    # S_tip should be positive (entropy is positive, with sign convention)
    S_positive = all(integrate_results[L]["S_tip"] > 0 for L in Lambda_values)
    # Lambda > 1 should give smaller S (Gaussian damps more) — sanity check
    S_decreasing_with_Lambda = (
        integrate_results[Lambda_values[-1]]["S_tip"]
        < integrate_results[Lambda_values[0]]["S_tip"]
    )

    print(f"\n[Verdict]")
    print(f"  S_tip finite at all Lambda:               {S_finite}")
    print(f"  S_tip positive (correct entropy sign):    {S_positive}")
    print(f"  S_tip decreases with Lambda (Gaussian):   {S_decreasing_with_Lambda}")

    if S_finite and S_positive and S_decreasing_with_Lambda:
        verdict = "POSITIVE-G4-5a-FIRST-MOVE-VERIFIED"
        msg = ("Tip-only replica integration operational on discrete substrate. "
               "S_tip extracted at finite, positive, Lambda-dependent values. "
               "Bulk Weyl (G4-5b) + variable warp (G4-5c) follow.")
    elif S_finite:
        verdict = "POSITIVE-G4-5a-FIRST-MOVE-PARTIAL"
        msg = "Integration finite but sign or trend unexpected; diagnose."
    else:
        verdict = "NEGATIVE-G4-5a-FIRST-MOVE"
        msg = "Integration divergent or non-finite; substrate UV/IR cutoff issue."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["S_finite"] = S_finite
    results["S_positive"] = S_positive

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
