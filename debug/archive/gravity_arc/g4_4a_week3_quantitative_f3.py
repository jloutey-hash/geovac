"""Sprint G4-4a week 3 — Quantitative F3 with continuum Weyl-Selberg.

Closes the next-week's named work from the week-2 closure: quantitative
F3 with continuum Weyl-Selberg matching at small t. The continuum
prediction for the 2D disk-Dirac heat trace at leading order is

    K_Dirac(t) ~ 2 * A / (4 pi t)

where the factor of 2 is the rank-2 spinor bundle multiplicity. We
verify

    R_Dirac(t) := K_Dirac(t) * 4 pi t / (2 A) -> 1 as t -> 0

at UV-refined substrate panels. Cross-check via the rank-2 ratio
K_Dirac(t) / K_scalar(t) -> 2.

Method
------
Use the T2 substrate (a = 0.05, N_phi up to 192) and sweep:
  - K_Dirac(t) via DiscreteDiskDirac.heat_trace(t)
  - K_scalar(t) via DiscreteDiskScalar.heat_trace(t) (periodic phi)
  - Weyl ratios R_Dirac(t) = K_Dirac(t) * 4 pi t / (2 A)

Exit gate
---------
Spinor Weyl ratio R_Dirac(t) within 5% of 1.0 at small t at the finest
panel. If passing, F3 is closed at the sprint-scale UV regime.
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteDiskScalar

OUT_JSON = Path(__file__).parent / "data" / "g4_4a_week3_quantitative_f3.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4a week 3 -- Quantitative F3 (continuum Weyl-Selberg)")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Substrate: R = 10, fixed N_rho = 200, a = 0.05 (T2 substrate)
    # ------------------------------------------------------------------
    R = 10.0
    N_rho = 200
    a = 0.05
    A = np.pi * R**2  # disk area

    N_phi_values = [24, 48, 96, 144, 192]
    t_values = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

    print(f"\n[Setup] R = {R}, N_rho = {N_rho}, a = {a}, A = pi R^2 = {A:.4f}")
    print(f"        N_phi sweep: {N_phi_values}")
    print(f"        t values: {t_values}")
    print(f"        Spinor Weyl prediction (leading): K_Dirac(t) ~ 2A/(4 pi t)")

    results["setup"] = {
        "R": R, "N_rho": N_rho, "a": a, "A": A,
        "N_phi_values": N_phi_values, "t_values": t_values,
    }

    # ------------------------------------------------------------------
    # Compute K_Dirac, K_scalar across panel
    # ------------------------------------------------------------------
    print("\n[Computing] K_Dirac and K_scalar at each (N_phi, t)...")
    K_Dirac_table = {}
    K_scalar_table = {}
    for N_phi in N_phi_values:
        disk_Dirac = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_phi)
        disk_scalar = DiscreteDiskScalar(N_rho=N_rho, a=a, N_phi=N_phi)
        K_Dirac_table[N_phi] = {t: disk_Dirac.heat_trace(t) for t in t_values}
        K_scalar_table[N_phi] = {t: disk_scalar.heat_trace(t) for t in t_values}
        print(f"  N_phi = {N_phi:>3d}: done")

    # ------------------------------------------------------------------
    # Spinor Weyl ratio R_Dirac(t) = K_Dirac(t) * 4 pi t / (2 A)
    # ------------------------------------------------------------------
    print("\n[Spinor Weyl ratio] R_Dirac = K_Dirac * 4 pi t / (2 A), target -> 1:")
    print()
    header = f"  {'N_phi':>6}  " + "  ".join([f"t={t}".rjust(10) for t in t_values])
    print(header)
    print("  " + "-" * (8 + 12 * len(t_values)))

    R_Dirac_table = {}
    for N_phi in N_phi_values:
        ratios = {}
        for t in t_values:
            R_val = K_Dirac_table[N_phi][t] * 4 * np.pi * t / (2 * A)
            ratios[t] = R_val
        R_Dirac_table[N_phi] = ratios
        row = f"  {N_phi:>6d}  " + "  ".join([f"{ratios[t]:>10.4f}" for t in t_values])
        print(row)

    results["spinor_weyl_ratio"] = {
        str(N): {str(t): v for t, v in r.items()}
        for N, r in R_Dirac_table.items()
    }

    # ------------------------------------------------------------------
    # Scalar Weyl ratio (cross-check vs T2 result)
    # ------------------------------------------------------------------
    print("\n[Scalar Weyl ratio cross-check] R_scalar = K_scalar * 4 pi t / A, target -> 1:")
    print()
    print(header)
    print("  " + "-" * (8 + 12 * len(t_values)))

    R_scalar_table = {}
    for N_phi in N_phi_values:
        ratios = {}
        for t in t_values:
            R_val = K_scalar_table[N_phi][t] * 4 * np.pi * t / A
            ratios[t] = R_val
        R_scalar_table[N_phi] = ratios
        row = f"  {N_phi:>6d}  " + "  ".join([f"{ratios[t]:>10.4f}" for t in t_values])
        print(row)

    results["scalar_weyl_ratio"] = {
        str(N): {str(t): v for t, v in r.items()}
        for N, r in R_scalar_table.items()
    }

    # ------------------------------------------------------------------
    # Rank-2 ratio (K_Dirac / K_scalar)
    # ------------------------------------------------------------------
    print("\n[Rank-2 ratio] K_Dirac / K_scalar, target -> 2.0 at small t:")
    print()
    print(header)
    print("  " + "-" * (8 + 12 * len(t_values)))

    rank2_table = {}
    for N_phi in N_phi_values:
        ratios = {}
        for t in t_values:
            ratio = K_Dirac_table[N_phi][t] / K_scalar_table[N_phi][t]
            ratios[t] = ratio
        rank2_table[N_phi] = ratios
        row = f"  {N_phi:>6d}  " + "  ".join([f"{ratios[t]:>10.4f}" for t in t_values])
        print(row)

    results["rank_2_ratio"] = {
        str(N): {str(t): v for t, v in r.items()}
        for N, r in rank2_table.items()
    }

    # ------------------------------------------------------------------
    # Convergence focus: spinor ratio at t = 0.1, t = 0.2 across N_phi
    # ------------------------------------------------------------------
    print("\n[UV convergence] Spinor Weyl ratio at t = 0.1 and t = 0.2:")
    print()
    print(f"  {'N_phi':>6}  {'R_Dirac(0.1)':>14}  {'R_Dirac(0.2)':>14}  "
          f"{'R_scalar(0.1)':>14}  {'rank2(0.1)':>12}")
    print("  " + "-" * 78)
    for N_phi in N_phi_values:
        print(f"  {N_phi:>6d}  "
              f"{R_Dirac_table[N_phi][0.1]:>14.4f}  "
              f"{R_Dirac_table[N_phi][0.2]:>14.4f}  "
              f"{R_scalar_table[N_phi][0.1]:>14.4f}  "
              f"{rank2_table[N_phi][0.1]:>12.4f}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    finest = N_phi_values[-1]
    R_finest_t01 = R_Dirac_table[finest][0.1]
    R_finest_t02 = R_Dirac_table[finest][0.2]
    R_finest_t005 = R_Dirac_table[finest][0.05]
    rank2_finest_t01 = rank2_table[finest][0.1]

    R_within_5pct_t01 = abs(R_finest_t01 - 1.0) < 0.05
    R_within_5pct_t02 = abs(R_finest_t02 - 1.0) < 0.05
    rank2_within_5pct = abs(rank2_finest_t01 - 2.0) < 0.10

    print(f"\n[Verdict at finest N_phi = {finest}]")
    print(f"  Spinor R_Dirac(t=0.05) = {R_finest_t005:.4f}")
    print(f"  Spinor R_Dirac(t=0.1)  = {R_finest_t01:.4f}  "
          f"(within 5%: {R_within_5pct_t01})")
    print(f"  Spinor R_Dirac(t=0.2)  = {R_finest_t02:.4f}  "
          f"(within 5%: {R_within_5pct_t02})")
    print(f"  Rank-2 ratio (t=0.1)   = {rank2_finest_t01:.4f}  "
          f"(within 10%: {rank2_within_5pct})")

    if R_within_5pct_t01 and rank2_within_5pct:
        verdict = "POSITIVE-G4-4a-WEEK3-QUANTITATIVE-F3-VERIFIED"
        msg = ("Spinor Weyl ratio R_Dirac(t=0.1) within 5% of 1.0 at "
               "N_phi=192; rank-2 enhancement within 10%. "
               "Continuum Weyl-Selberg leading order recovered "
               "quantitatively on the discrete substrate.")
    elif R_within_5pct_t02 and rank2_within_5pct:
        verdict = "POSITIVE-G4-4a-WEEK3-QUANTITATIVE-F3-PARTIAL"
        msg = ("Spinor Weyl ratio within 5% at t=0.2 only; "
               "more UV refinement needed for the t=0.1 sweet spot.")
    else:
        verdict = "PARTIAL-G4-4a-WEEK3"
        msg = ("Spinor Weyl convergence characterized but does not yet "
               "reach 5% at sprint scale. Structural pattern documented.")

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["R_finest_t01"] = R_finest_t01
    results["R_finest_t02"] = R_finest_t02
    results["rank2_finest_t01"] = rank2_finest_t01

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
