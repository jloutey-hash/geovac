"""Sprint G4-6a refined first move -- simplified A extraction.

Per Track alpha of thread 6 (2026-05-29): G4-6b first move found
B_substrate ~ 0.163 at R = 10 (matching continuum +1/6 to within 2.3%).
The B.2 small-t-panel linear fit produced B_fit = 0.290-0.318 because
joint A/B fit was poorly conditioned by substrate's UV undershoot in A.

This driver reprocesses the B.2 spectral panel data with the simplified
strategy:
1. Use B_substrate = 0.163 (measured at large t).
2. Compute A_estimate(t) = t * (tip(t) - B_substrate) at each small-t.
3. Richardson-extrapolate A across substrate refinements.

Goal: verify that the A extraction now has constraining power
(vs B.2's joint fit which gave A ~ 0 at both panels).
"""

import json
from pathlib import Path

import numpy as np

OUT_JSON = Path(__file__).parent / "data" / "g4_6a_refined_simplified.json"
B_2_DATA = Path(__file__).parent / "data" / "g4_6a_multi_substrate_uv_first_move_spectral.json"

# From G4-6b: B at large t (asymptotic in R) for R >= 10
B_SUBSTRATE = 0.163  # value at R = 10 panel; R-independent for R >= 10
A_CONT = 1.0 / (24.0 * np.pi)


def reprocess_panel(panel: dict, B_subs: float) -> dict:
    """Apply simplified A extraction to one substrate panel.

    For each (t, tip) datapoint, compute:
        A_estimate(t) = t * (tip(t) - B_subs)

    Average A across the small-t window.
    """
    raw = panel["raw_data"]
    per_t_A = []
    for entry in raw:
        t = entry["t"]
        tip = entry["tip_term"]
        # Simplified extraction: A/t = tip - B, so A = t * (tip - B)
        A_est = t * (tip - B_subs)
        per_t_A.append({
            "t": t,
            "tip": tip,
            "tip_minus_B": tip - B_subs,
            "A_est": A_est,
            "A_recovery_vs_cont": A_est / A_CONT,
        })

    # Average and standard deviation
    A_arr = np.array([entry["A_est"] for entry in per_t_A])
    return {
        "label": panel["label"],
        "a": panel["a"],
        "N_rho": panel["N_rho"],
        "N_0": panel["N_0"],
        "per_t_extraction": per_t_A,
        "A_mean": float(np.mean(A_arr)),
        "A_std": float(np.std(A_arr)),
        "A_median": float(np.median(A_arr)),
        # Weighted by 1/t (downweight intermediate t where B dominates)
        # Use small-t points (t <= 5*a^2) for cleaner UV signal
        "A_smallt_mean": float(np.mean([
            entry["A_est"] for entry in per_t_A if entry["t"] <= 5 * panel["a"]**2
        ])),
    }


def main() -> None:
    print("=" * 72)
    print("Sprint G4-6a refined -- Simplified A extraction first move")
    print("=" * 72)

    # Load B.2 data
    with open(B_2_DATA) as f:
        b2 = json.load(f)

    print(f"\n  B_substrate = {B_SUBSTRATE:.6f} (from G4-6b at R = 10 panel)")
    print(f"  Continuum A target = 1/(24*pi) = {A_CONT:.6f}")
    print(f"  Strategy: A_est(t) = t * (tip(t) - B_substrate)")

    panels = b2["panels"]
    refined_panels = []

    for panel in panels:
        print(f"\n{'=' * 72}")
        print(f"Panel {panel['label']}: a={panel['a']}, N_rho={panel['N_rho']}")
        print(f"{'=' * 72}")
        refined = reprocess_panel(panel, B_SUBSTRATE)
        refined_panels.append(refined)

        print(f"\n  Per-t A extraction:")
        print(f"  {'t':>10}  {'tip(t)':>12}  {'tip - B_subs':>14}  "
              f"{'A_est':>12}  {'recovery':>10}")
        for entry in refined["per_t_extraction"]:
            print(f"  {entry['t']:>10.5f}  {entry['tip']:>+12.6f}  "
                  f"{entry['tip_minus_B']:>+14.6f}  {entry['A_est']:>+12.6f}  "
                  f"{entry['A_recovery_vs_cont']*100:>+9.2f}%")

        print(f"\n  Statistical summary:")
        print(f"    A_mean         = {refined['A_mean']:+.6f} "
              f"({refined['A_mean']/A_CONT*100:.2f}% of A_cont)")
        print(f"    A_median       = {refined['A_median']:+.6f} "
              f"({refined['A_median']/A_CONT*100:.2f}% of A_cont)")
        print(f"    A_std          = {refined['A_std']:.6f}")
        print(f"    A_smallt_mean  = {refined['A_smallt_mean']:+.6f} "
              f"({refined['A_smallt_mean']/A_CONT*100:.2f}% of A_cont)")
        print(f"    [smallt: only t <= 5*a^2 = {5*panel['a']**2:.5f}]")

    # ---------------------------------------------------------------------
    # Richardson extrapolation
    # ---------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Richardson extrapolation: refined simplified strategy")
    print(f"{'=' * 72}")

    # Use A_smallt_mean (best UV signal) for Richardson
    A1 = refined_panels[0]["A_smallt_mean"]
    A2 = refined_panels[1]["A_smallt_mean"]
    a1 = refined_panels[0]["a"]
    a2 = refined_panels[1]["a"]

    print(f"\n  A_smallt at each panel:")
    print(f"    Panel 1 (a={a1}): A = {A1:+.6f} ({A1/A_CONT*100:.2f}% of A_cont)")
    print(f"    Panel 2 (a={a2}): A = {A2:+.6f} ({A2/A_CONT*100:.2f}% of A_cont)")

    diff_1 = A1 - A_CONT
    diff_2 = A2 - A_CONT
    print(f"\n  Diff from A_cont:")
    print(f"    Panel 1: A1 - A_cont = {diff_1:+.6f}")
    print(f"    Panel 2: A2 - A_cont = {diff_2:+.6f}")

    if diff_1 != 0 and diff_2 != 0 and diff_1 * diff_2 > 0:
        ratio = diff_1 / diff_2
        if ratio > 0:
            p_implied = np.log2(ratio)
            print(f"    Implied convergence exponent p = log2(ratio) = {p_implied:.4f}")
            # Richardson: A_inf = (2^p * A2 - A1) / (2^p - 1)
            two_p = 2.0 ** p_implied
            if abs(two_p - 1) > 1e-10:
                A_richardson = (two_p * A2 - A1) / (two_p - 1.0)
                print(f"    Richardson A_inf = {A_richardson:+.6f} "
                      f"({A_richardson/A_CONT*100:.2f}% of A_cont)")
            else:
                A_richardson = float("nan")
        else:
            p_implied = float("nan")
            A_richardson = float("nan")
    else:
        p_implied = float("nan")
        A_richardson = float("nan")

    # Richardson at p = 1 (linear convergence)
    A_richardson_p1 = 2 * A2 - A1
    print(f"\n  Assuming p = 1 (linear convergence in a):")
    print(f"    A_inf = 2*A2 - A1 = {A_richardson_p1:+.6f} "
          f"({A_richardson_p1/A_CONT*100:.2f}% of A_cont)")

    # Richardson at p = 2 (quadratic convergence in a)
    A_richardson_p2 = (4 * A2 - A1) / 3.0
    print(f"\n  Assuming p = 2 (quadratic convergence in a):")
    print(f"    A_inf = (4*A2 - A1)/3 = {A_richardson_p2:+.6f} "
          f"({A_richardson_p2/A_CONT*100:.2f}% of A_cont)")

    # ---------------------------------------------------------------------
    # Verdict
    # ---------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Verdict")
    print(f"{'=' * 72}")

    # Compare to B.2 joint-fit result
    print(f"\n  B.2 joint A/B fit gave:")
    print(f"    A1 (a=0.05)  = +0.002883 (21.74% of A_cont)")
    print(f"    A2 (a=0.025) = +0.000716 (5.40% of A_cont)")
    print(f"    Trajectory: REVERSED (A2 < A1; Richardson trajectory bad)")

    print(f"\n  G4-6a refined simplified extraction:")
    print(f"    A1 = {A1:+.6f} ({A1/A_CONT*100:.2f}% of A_cont)")
    print(f"    A2 = {A2:+.6f} ({A2/A_CONT*100:.2f}% of A_cont)")

    if abs(A2) > abs(A1) and A2 / A_CONT > A1 / A_CONT:
        print(f"    Trajectory: TOWARD A_cont (A2 > A1 in recovery)")
        print(f"    Status: POSITIVE -- simplified strategy has constraining power.")
    elif abs(A2) > abs(A1) * 1.1:
        print(f"    Trajectory: monotonic improvement.")
        print(f"    Status: POSITIVE -- simplified strategy provides cleaner A extraction.")
    else:
        print(f"    Trajectory: weak; substrate may not extend to small-enough a for")
        print(f"    sub-percent A extraction at current resolution.")

    results = {
        "B_substrate": B_SUBSTRATE,
        "A_cont": A_CONT,
        "panels_refined": refined_panels,
        "Richardson": {
            "A1": float(A1),
            "A2": float(A2),
            "a1": a1,
            "a2": a2,
            "p_implied": float(p_implied) if not np.isnan(p_implied) else None,
            "A_richardson_implied": float(A_richardson) if not np.isnan(A_richardson) else None,
            "A_richardson_p1": float(A_richardson_p1),
            "A_richardson_p2": float(A_richardson_p2),
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
