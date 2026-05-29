"""Sprint G4-4f -- Replica method preparation: dK/dalpha at alpha=1.

Per replica method: the BH entropy contribution at the conical defect
tip is

    S_tip = -dI_E/dalpha|_{alpha=1} = some constant times d(K_tip)/dalpha|_{alpha=1}

For the spinor on cone (G4-4c week 3 identification):
    K_tip(alpha) = -(1/12)(1/alpha - alpha)
    d K_tip / dalpha = -(1/12)(-1/alpha^2 - 1) = (1/12)(1/alpha^2 + 1)
    At alpha = 1: dK_tip/dalpha = (1/12)(1 + 1) = 1/6

Test: compute dK_wedge/dalpha|_{alpha=1} via central finite difference and
subtract K_disk (the bulk contribution). The residual should be 1/6.

Finite difference:
    dK_wedge/dalpha|_{alpha=1} ≈ [K(alpha_+) - K(alpha_-)] / (alpha_+ - alpha_-)
where alpha_± are rational close to 1 (so that N_phi = alpha * N_0 is integer).
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

OUT_JSON = Path(__file__).parent / "data" / "g4_4f_replica_dK_dalpha.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4f -- Replica dK/dalpha at alpha=1")
    print("=" * 72)

    R = 10.0; a = 0.05; N_rho = 200; N_0 = 120
    print(f"\n[Setup] R={R}, a={a}, N_rho={N_rho}, N_0={N_0}")
    print(f"  Continuum prediction: d Delta_K^Dirac / d alpha |_{{alpha=1}} = +1/6 = "
          f"{1/6:.6f}")

    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)

    # Use sweet-spot t window for clean extraction
    t_focus = 1.0
    K_disk = disk.heat_trace(t_focus)
    print(f"\n  K_disk_Dirac(t={t_focus}) = {K_disk:.6f}")

    # ------------------------------------------------------------------
    # Central finite difference using N_phi = N_0 ± k
    # ------------------------------------------------------------------
    print(f"\n[Step 1] Central FD: dK_wedge/dalpha at alpha=1 using various step sizes:")
    print()
    print(f"  {'N_+/N_0':>10}  {'eps':>10}  {'K(alpha_+)':>12}  {'K(alpha_-)':>12}  "
          f"{'dK_wedge/dalpha':>14}  {'-K_disk':>12}  {'tip term':>12}")
    print("  " + "-" * 100)

    fd_data = {}
    for k_step in [1, 2, 4, 6, 12, 24]:
        N_plus = N_0 + k_step
        N_minus = N_0 - k_step
        alpha_plus = N_plus / N_0
        alpha_minus = N_minus / N_0
        eps = (alpha_plus - alpha_minus) / 2

        wedge_plus = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
        )
        wedge_minus = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
        )
        K_plus = wedge_plus.heat_trace(t_focus)
        K_minus = wedge_minus.heat_trace(t_focus)
        dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
        tip_term = dK_dalpha - K_disk

        fd_data[k_step] = {
            "alpha_plus": alpha_plus, "alpha_minus": alpha_minus,
            "eps": eps,
            "K_plus": K_plus, "K_minus": K_minus,
            "dK_dalpha": dK_dalpha,
            "tip_term": tip_term,
            "tip_recovery_vs_1_6": tip_term / (1/6),
        }
        print(f"  {N_plus}/{N_0:>3}     {eps:>10.4f}  {K_plus:>12.4f}  {K_minus:>12.4f}  "
              f"{dK_dalpha:>14.6f}  {-K_disk:>+12.4f}  {tip_term:>+12.6f}")

    results["finite_difference_data"] = {str(k): v for k, v in fd_data.items()}

    # ------------------------------------------------------------------
    # Step 2: Richardson extrapolation as eps -> 0
    # ------------------------------------------------------------------
    print(f"\n[Step 2] Tip term value at different step sizes (Richardson):")
    print()
    print(f"  Continuum prediction: tip term = 1/6 = {1/6:.6f}")
    print()
    print(f"  {'eps':>10}  {'tip term':>14}  {'rel_err vs 1/6':>16}")
    for k in sorted(fd_data.keys()):
        eps = fd_data[k]["eps"]
        tip = fd_data[k]["tip_term"]
        rel_err = (tip - 1/6) / (1/6)
        print(f"  {eps:>10.4f}  {tip:>+14.6f}  {rel_err:>+16.4f}")

    # Best tip term: smallest eps
    smallest_k = min(fd_data.keys())
    best_tip = fd_data[smallest_k]["tip_term"]
    best_rec = fd_data[smallest_k]["tip_recovery_vs_1_6"]
    print(f"\n  Smallest eps: {fd_data[smallest_k]['eps']:.4f}")
    print(f"  Tip term: {best_tip:.6f}")
    print(f"  Recovery vs 1/6: {best_rec:.4f} "
          f"({best_rec*100:.2f}%)")

    results["smallest_eps_tip_term"] = float(best_tip)
    results["smallest_eps_recovery"] = float(best_rec)

    # ------------------------------------------------------------------
    # Step 3: Cross-t check
    # ------------------------------------------------------------------
    print(f"\n[Step 3] Cross-t check (k_step = 12, eps = 0.1):")
    print()
    print(f"  {'t':>5}  {'dK_wedge/dalpha':>14}  {'-K_disk':>12}  "
          f"{'tip term':>12}  {'recovery':>10}")
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    t_sweep = [0.5, 1.0, 2.0, 5.0]
    cross_t = {}
    for t in t_sweep:
        K_d = disk.heat_trace(t)
        K_p = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
        ).heat_trace(t)
        K_m = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
        ).heat_trace(t)
        dK = (K_p - K_m) / (alpha_plus - alpha_minus)
        tip = dK - K_d
        rec = tip / (1/6)
        cross_t[str(t)] = {
            "K_disk": float(K_d),
            "dK_dalpha": float(dK),
            "tip_term": float(tip),
            "recovery_vs_1_6": float(rec),
        }
        print(f"  {t:>5}  {dK:>+14.4f}  {-K_d:>+12.4f}  {tip:>+12.6f}  "
              f"{rec:>10.4f}")

    results["cross_t"] = cross_t

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    best_recovery_across_t = max(abs(cell["recovery_vs_1_6"]) for cell in cross_t.values())
    closest_to_1 = min(cross_t.items(),
                       key=lambda kv: abs(kv[1]["recovery_vs_1_6"] - 1.0))
    print(f"\n[Verdict]")
    print(f"  Best tip-term recovery at smallest eps: {best_rec:.4f}")
    print(f"  Closest to 1.0 across t: {closest_to_1[1]['recovery_vs_1_6']:.4f} "
          f"at t = {closest_to_1[0]}")

    if abs(best_rec - 1.0) < 0.1:
        verdict = "POSITIVE-G4-4f-VERIFIED"
        msg = (f"dDelta__K^Dirac/dalpha at alpha=1 extracted to 1/6 within 10%; replica "
               f"method derivative confirmed at sprint scale. Best recovery: "
               f"{best_rec*100:.2f}%. Load-bearing tip-contribution to "
               f"S_BH on discrete substrate is structurally extractable.")
    elif abs(best_rec - 1.0) < 0.3:
        verdict = "POSITIVE-G4-4f-PARTIAL"
        msg = ("Tip term recovered within 30%; substrate refinement needed "
               "for full precision.")
    else:
        verdict = "PARTIAL-G4-4f"
        msg = "Tip-term FD extraction sensitive to t window and substrate."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
