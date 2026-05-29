"""Sprint Moebius mechanism diagnostic via discrete mode decomposition.

Task 11 of the prior thread identified the harmonic-conjugate algebraic
structure: F(alpha) = alpha/(2alpha-1) satisfies 1/alpha + 1/F = 2. Task 14 sharpened
with a structural-derivation attempt. Task 25 (this driver) tests
whether the soft-IR mode m=0 (lowest spinor angular eigenvalue
1/(2alpha)) provides the dominant Moebius contribution.

Approach:
1. Decompose K_wedge(alpha, t) by angular mode m on the spectral substrate.
2. Compute per-mode contributions to dK/dalpha at multiple alpha values.
3. Check whether (a) the soft-IR mode dominates the variation, (b) the
   variation across alpha follows the Moebius alpha/(2alpha-1) pattern.

The mode-resolved heat trace contribution at mode k_idx is:
  K_m(alpha, t) = 2 * sum_radial exp(-t * eigenvalue_m_radial)
where eigenvalue_m_radial uses m_eff = (k+1/2)/alpha and the 2 factor is
spinor doubling.
"""

import json
import time
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import DiscreteWedgeDiracSpectral

OUT_JSON = Path(__file__).parent / "data" / "sprint_moebius_mode_decomposition.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def mode_resolved_heat_trace(N_rho: int, a: float, N_phi: int, alpha: float,
                              t: float):
    """Compute per-mode contributions to the wedge spectral heat trace.

    Returns list of dicts with mode index k, m_eff, and per-mode K_m(t).
    """
    wedge = DiscreteWedgeDiracSpectral(
        N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha,
    )
    mode_contribs = []
    for k_idx in range(N_phi):
        if k_idx <= N_phi // 2:
            k = k_idx
        else:
            k = k_idx - N_phi
        m_eff = float((k + 0.5) / alpha)
        H_rad = wedge._hermitian_radial_laplacian(m_eff)
        evals = np.linalg.eigvalsh(H_rad)
        # Spinor doubling: each eigenvalue counted twice
        K_m = float(2.0 * np.sum(np.exp(-evals * t)))
        mode_contribs.append({
            "k_idx": k_idx,
            "k_signed": k,
            "m_eff": m_eff,
            "K_m": K_m,
        })
    return mode_contribs


def compute_per_mode_dK_dalpha(N_rho: int, a: float, N_0: int, t: float,
                                k_step: int = 12) -> list:
    """Compute per-mode dK/dalpha at alpha=1 via central FD.

    The wedge with N_phi = N_0 ± k_step modes has alpha = (N_0 ± k_step)/N_0.
    Each mode-pair gives dK_m/dalpha contribution.

    Note: at alpha_+ and alpha_-, the N_phi values differ, so we have different
    mode-index ranges. We match by SIGNED k (the angular momentum
    label) to compute the per-mode contribution properly.
    """
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0

    modes_plus = mode_resolved_heat_trace(N_rho, a, N_plus, alpha_plus, t)
    modes_minus = mode_resolved_heat_trace(N_rho, a, N_minus, alpha_minus, t)

    # Build dict by signed k for matching
    plus_by_k = {m["k_signed"]: m for m in modes_plus}
    minus_by_k = {m["k_signed"]: m for m in modes_minus}

    # Intersection: modes available at both
    common_k = sorted(set(plus_by_k.keys()) & set(minus_by_k.keys()))

    per_mode_dK = []
    for k in common_k:
        K_plus = plus_by_k[k]["K_m"]
        K_minus = minus_by_k[k]["K_m"]
        dK_dalpha = (K_plus - K_minus) / (alpha_plus - alpha_minus)
        per_mode_dK.append({
            "k_signed": k,
            "K_plus": K_plus,
            "K_minus": K_minus,
            "dK_m_dalpha": dK_dalpha,
            "m_eff_at_1": (k + 0.5),  # m_eff at alpha=1
        })

    return per_mode_dK


def main() -> None:
    print("=" * 72)
    print("Sprint Moebius mode decomposition diagnostic")
    print("=" * 72)

    results = {}

    # Setup: small substrate for fast iteration
    N_rho, a, N_0 = 100, 0.05, 60
    t = 1.0  # Sweet-spot t window for clean dK/dalpha at alpha=1
    print(f"\n  Substrate: N_rho={N_rho}, a={a}, N_0={N_0}, t={t}")
    print(f"  Approach: central FD dK/dalpha at alpha=1, decomposed per angular mode")
    print(f"  Continuum dK/dalpha at alpha=1: +1/6 = {1/6:.6f}")

    # ------------------------------------------------------------------
    # Step 1: Per-mode dK/dalpha at alpha=1
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 1: Per-mode dK/dalpha at alpha=1")
    print(f"{'=' * 72}")

    t0 = time.time()
    per_mode = compute_per_mode_dK_dalpha(N_rho, a, N_0, t, k_step=12)
    print(f"  Compute time: {time.time() - t0:.1f}s")
    print(f"  {len(per_mode)} modes contributing (intersection of N_+ and N_- mode sets)")

    # Sort by |m_eff_at_1| to identify dominant modes
    per_mode_sorted = sorted(per_mode, key=lambda m: abs(m["m_eff_at_1"]))

    print(f"\n  Top 10 modes by |m_eff| (most IR):")
    print(f"  {'k':>5}  {'m_eff(alpha=1)':>12}  {'dK_m/dalpha':>14}  {'fraction':>10}")
    total_dK = sum(m["dK_m_dalpha"] for m in per_mode)
    for m in per_mode_sorted[:10]:
        frac = m["dK_m_dalpha"] / total_dK if total_dK != 0 else 0
        print(f"  {m['k_signed']:>5d}  {m['m_eff_at_1']:>12.3f}  "
              f"{m['dK_m_dalpha']:>+14.6f}  {frac:>10.4f}")

    # Total dK/dalpha across all modes
    print(f"\n  Total dK/dalpha (sum across all modes): {total_dK:.6f}")
    print(f"  Continuum target: 1/6 = {1/6:.6f}")
    print(f"  Recovery: {total_dK / (1/6):.4f} ({total_dK / (1/6) * 100:.2f}%)")

    # Soft-IR contribution (m_eff in [-1/2, 1/2])
    soft_ir = [m for m in per_mode if abs(m["m_eff_at_1"]) <= 1.0]
    soft_ir_dK = sum(m["dK_m_dalpha"] for m in soft_ir)
    soft_ir_frac = soft_ir_dK / total_dK if total_dK != 0 else 0
    print(f"\n  Soft-IR modes (|m_eff(alpha=1)| <= 1.0): "
          f"sum dK/dalpha = {soft_ir_dK:.6f} ({soft_ir_frac*100:.2f}% of total)")

    results["t"] = t
    results["substrate"] = {"N_rho": N_rho, "a": a, "N_0": N_0}
    results["total_dK_dalpha"] = float(total_dK)
    results["continuum_target"] = float(1/6)
    results["recovery"] = float(total_dK / (1/6))
    results["soft_IR_sum"] = float(soft_ir_dK)
    results["soft_IR_fraction"] = float(soft_ir_frac)
    results["per_mode"] = per_mode

    # ------------------------------------------------------------------
    # Step 2: Vary alpha and check Moebius pattern
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 2: Per-mode contributions at multiple alpha (excess/deficit)")
    print(f"{'=' * 72}")

    # For each alpha, compute the total heat trace tip term
    # and decompose by soft-IR vs other.
    alpha_values = [0.5, 1.0, 1.5, 2.0, 3.0]
    alpha_data = []

    for alpha in alpha_values:
        # Compute K_wedge(alpha, t) at this alpha with N_phi = alpha * N_0 (integer)
        N_phi_alpha = int(round(alpha * N_0))
        actual_alpha = N_phi_alpha / N_0
        modes_alpha = mode_resolved_heat_trace(
            N_rho, a, N_phi_alpha, actual_alpha, t,
        )
        K_total = sum(m["K_m"] for m in modes_alpha)
        soft_ir_K = sum(
            m["K_m"] for m in modes_alpha
            if abs(m["m_eff"]) <= 1.0 / actual_alpha + 1e-6
        )
        # Moebius factor at this alpha (for alpha > 1)
        if actual_alpha > 1:
            F_moebius = actual_alpha / (2 * actual_alpha - 1)
        elif actual_alpha < 1:
            F_moebius = 1.0  # Standard SC, no Moebius
        else:
            F_moebius = 1.0
        # SC continuum prediction
        SC_pred = -(1/12) * (1/actual_alpha - actual_alpha)
        Moebius_pred = SC_pred * F_moebius

        alpha_data.append({
            "alpha": actual_alpha,
            "N_phi": N_phi_alpha,
            "K_total": K_total,
            "soft_IR_K": soft_ir_K,
            "soft_IR_fraction": soft_ir_K / K_total if K_total != 0 else 0,
            "SC_pred": SC_pred,
            "Moebius_pred": Moebius_pred,
            "F_moebius": F_moebius,
        })

        print(f"\n  alpha = {actual_alpha:.4f} (N_phi = {N_phi_alpha})")
        print(f"    K_total                  = {K_total:.4f}")
        print(f"    Soft-IR K contribution   = {soft_ir_K:.4f} "
              f"({soft_ir_K/K_total*100:.2f}%)")
        print(f"    Standard SC prediction   = {SC_pred:+.6f}")
        print(f"    Moebius F(alpha)              = {F_moebius:.4f}")
        print(f"    Moebius-modified SC pred  = {Moebius_pred:+.6f}")

    results["alpha_sweep"] = alpha_data

    # ------------------------------------------------------------------
    # Step 3: Verdict
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Step 3: Verdict")
    print(f"{'=' * 72}")

    # Check whether soft-IR fraction varies with alpha in a way that
    # mirrors the Moebius factor
    print(f"\n  alpha-dependence of soft-IR fraction:")
    print(f"  {'alpha':>8}  {'soft-IR frac':>14}  {'F_moebius':>10}")
    for d in alpha_data:
        print(f"  {d['alpha']:>8.4f}  {d['soft_IR_fraction']:>14.4f}  {d['F_moebius']:>10.4f}")

    # Check if soft-IR fraction at alpha > 1 trends as Moebius
    if total_dK > 0:
        if soft_ir_frac > 0.5:
            print(f"\n  [PARTIAL POSITIVE] Soft-IR modes carry > 50% of dK/dalpha.")
            print(f"  Consistent with soft-IR mode being load-bearing for tip term.")
            print(f"  Mechanism sketch (Task 11 §3.2) gains supporting evidence.")
        else:
            print(f"\n  [NEUTRAL] Soft-IR modes carry < 50% of dK/dalpha.")
            print(f"  Mechanism is distributed across many modes; soft-IR alone")
            print(f"  not the dominant fingerprint.")
    else:
        print(f"\n  [INCONCLUSIVE] Sign/magnitude unusual; need follow-up.")

    # Save
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved results to {OUT_JSON}")


if __name__ == "__main__":
    main()
