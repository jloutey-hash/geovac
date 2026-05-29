"""Sprint G4-4b-a first move (variable warp Dirac, Level 1).

Closes the named first-move work from the G4-4b scoping memo:
  - VariableWarpDirac at Level 1 (leading position-dependent S^2 mass)
  - F4 tip-regularity at the warp-derivative array
  - F6 Riemannian-limit bit-exact reduction at constant r(rho) = r_h
  - F7 first-pass factorization-loss quantification

Setup
-----
Smooth-tip warp r(rho) = r_h * sqrt(1 + (rho/r_h)^2) on the G4-3a-cleanup
Hermitian polar substrate. The cigar squared-Dirac at variable warp is
diagonalized per (S^2 mode n, azimuthal Fourier mode k_phi) block:

    H_{n, k_phi} = L_disk(k_phi) + diag((n+1)^2 / r(rho)^2)

Each radial eigenvalue contributes with multiplicity 16(n+1).

Exit gates
----------
- F4: warp-derivative array finite + tip value within 50% of rho_1/r_h^2
- F6: K_var(t) - K_const(t) bit-exact (rel_err < 1e-10) at constant r
- F7: K_var > K_const at variable warp (positive Delta_fact)
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    S2DiracSpectrum,
    VariableWarpDirac,
    WarpedDiracConstant,
    verify_F4_tip_regular,
    verify_F6_riemannian_limit,
    verify_F7_factorization_loss,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_4b_a_first_move.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4b-a first move -- Variable warp Dirac (Level 1)")
    print("=" * 72)

    # ------------------------------------------------------------------
    # F4: tip-regularity of warp-derivative term
    # ------------------------------------------------------------------
    print("\n[F4] Warp-derivative tip-regularity r'(rho)/r(rho) at rho -> 0:")
    f4_panels = [
        ("small",  10, 0.5, 8,  2.0),  # (Nr, a, Nphi, r_h)
        ("medium", 20, 0.3, 12, 2.0),
        ("larger", 30, 0.2, 16, 2.0),
    ]
    F4_results = {}
    for name, Nr, a, Np, r_h in f4_panels:
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        res = verify_F4_tip_regular(disk=disk, r_h=r_h, factor_tol=0.50)
        F4_results[name] = res
        marker = "PASS" if res["passed"] else "FAIL"
        print(f"  Panel '{name}' (a={a}, r_h={r_h}):")
        print(f"    r'/r at rho_1 = {res['deriv_at_rho_1']:.4e}  "
              f"(target rho_1/r_h^2 = {res['expected_at_rho_1']:.4e}, "
              f"rel_err = {res['rel_err_tip']:.3f}) [{marker}]")
        print(f"    max r'/r over substrate = {res['max_deriv_value']:.4e} "
              f"(finite: {res['deriv_array_finite']})")

    results["F4"] = F4_results

    # ------------------------------------------------------------------
    # F6: load-bearing Riemannian-limit at constant r(rho) = r_h
    # ------------------------------------------------------------------
    print("\n[F6 LOAD-BEARING] Constant-warp reduction at r(rho) = r_h:")
    f6_panels = [
        ("small",  10, 0.5, 6,  1.5, 2),
        ("medium", 20, 0.3, 12, 2.0, 3),
    ]
    t_panel = [0.05, 0.1, 0.5, 1.0]
    F6_results = {}
    for name, Nr, a, Np, r_h, l_max in f6_panels:
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        res = verify_F6_riemannian_limit(
            disk=disk, sphere=sphere, r_h=r_h, t_values=t_panel, tol=1e-10,
        )
        F6_results[name] = res
        print(f"\n  Panel '{name}' (Nr={Nr}, Np={Np}, l_max={l_max}, r_h={r_h}):")
        for t in t_panel:
            entry = res[str(t)]
            marker = "PASS" if entry["passed"] else "FAIL"
            print(f"    t={t:>4}: K_var={entry['K_var']:.4e}  "
                  f"K_const={entry['K_const_factorized']:.4e}  "
                  f"rel_err={entry['rel_err']:.2e} [{marker}]")
        print(f"  Panel '{name}' all_passed: {res['all_passed']}")

    results["F6"] = F6_results

    # ------------------------------------------------------------------
    # F7: first-pass factorization-loss at smooth-tip variable warp
    # ------------------------------------------------------------------
    print("\n[F7 first-pass] Factorization-loss at smooth-tip variable warp:")
    f7_panels = [
        ("small_rh_2", 20, 0.3, 12, 2.0, 3),  # tip-dominated regime
        ("medium_rh_5", 20, 0.3, 12, 5.0, 3),  # closer-to-constant regime
        ("medium_rh_1", 20, 0.3, 12, 1.0, 3),  # asymptotic-dominated regime
    ]
    F7_results = {}
    for name, Nr, a, Np, r_h, l_max in f7_panels:
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
        # Show warp profile statistics
        var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)
        warp_min = float(np.min(var.warp_profile))
        warp_max = float(np.max(var.warp_profile))
        warp_mean_relative = float(np.mean(var.warp_profile) / r_h)
        print(f"\n  Panel '{name}' (Nr={Nr}, Np={Np}, l_max={l_max}, r_h={r_h}):")
        print(f"    Warp profile: r in [{warp_min:.3f}, {warp_max:.3f}], "
              f"mean/r_h = {warp_mean_relative:.3f}")
        res = verify_F7_factorization_loss(
            disk=disk, sphere=sphere, r_h=r_h, t_values=t_panel,
        )
        F7_results[name] = res
        F7_results[name]["warp_min"] = warp_min
        F7_results[name]["warp_max"] = warp_max
        for t in t_panel:
            entry = res[str(t)]
            sign_marker = "(+)" if entry["sign_positive"] else "(-)"
            print(f"    t={t:>4}: K_var={entry['K_var']:.4e}  "
                  f"K_const={entry['K_const_factorized']:.4e}  "
                  f"Delta={entry['Delta_fact']:+.4e}  "
                  f"ratio={entry['ratio']:.4f} {sign_marker}")

    results["F7"] = F7_results

    # ------------------------------------------------------------------
    # F7 monotonicity in warp variation
    # ------------------------------------------------------------------
    print("\n[F7 monotonicity] Delta_fact should grow with warp variation:")
    print(f"  Ordering by tip parameter: small r_h = strong variation, "
          f"large r_h = weak variation")
    delta_at_t01 = {}
    for name, _, _, _, r_h, _ in f7_panels:
        Delta = F7_results[name]["0.1"]["Delta_fact"]
        warp_range = F7_results[name]["warp_max"] - F7_results[name]["warp_min"]
        delta_at_t01[name] = (r_h, warp_range, Delta)
        print(f"  r_h={r_h:>4}: warp range = {warp_range:.3f},  "
              f"Delta_fact(t=0.1) = {Delta:+.4e}")

    # Trend check
    r_h_sorted = sorted(delta_at_t01.items(), key=lambda kv: -kv[1][0])
    print(f"\n  Trend: as r_h decreases (more variation), Delta should grow.")
    for name, (r_h, var_range, delta) in r_h_sorted:
        print(f"    r_h={r_h:>4}: Delta_fact = {delta:+.4e}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    F4_all = all(F4_results[n]["passed"] for n, _, _, _, _ in f4_panels)
    F6_all = all(F6_results[n]["all_passed"] for n, _, _, _, _, _ in f6_panels)
    F7_all_positive = all(
        F7_results[n][str(t)]["sign_positive"]
        for n, _, _, _, _, _ in f7_panels
        for t in t_panel
    )

    print(f"\n[Verdict]")
    print(f"  F4 (tip-regularity):                       {F4_all}")
    print(f"  F6 (Riemannian-limit bit-exact, LOAD-BEARING): {F6_all}")
    print(f"  F7 (Delta_fact > 0 at variable warp):      {F7_all_positive}")

    if F4_all and F6_all and F7_all_positive:
        verdict = "POSITIVE-G4-4b-a-FIRST-MOVE-VERIFIED"
        msg = ("F4 tip-regularity + F6 Riemannian-limit bit-exact + F7 "
               "first-pass factorization-loss positive. Variable warp "
               "architecture at Level 1 closes at sprint scale. Next "
               "G4-4b weeks: F7 structural form quantification (b), F5 "
               "asymptotic-free recovery (c), explicit Level 2 with "
               "spin-connection (d).")
    elif F6_all:
        verdict = "POSITIVE-G4-4b-a-FIRST-MOVE-PARTIAL"
        msg = "F6 load-bearing reduction bit-exact; F4 or F7 needs refinement."
    else:
        verdict = "NEGATIVE-G4-4b-a-FIRST-MOVE"
        msg = "F6 load-bearing reduction failed; structural issue."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["F4_all"] = F4_all
    results["F6_all"] = F6_all
    results["F7_all_positive"] = F7_all_positive

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
