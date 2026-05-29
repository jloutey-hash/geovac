"""Sprint G4-3d-UV extension — UV regime via N_phi refinement.

Tests the G4-3d prediction (memo Sec. 5): the small-t saturation of the
Weyl ratio K(t) * 4 pi t / A_disk at N_phi = 24 is angular truncation,
not a Hermiticity issue. Predicted threshold: N_phi >~ 96 to reach
Weyl law at t = 0.01.

Method
------
Fix R = N_rho * a = 10, fix radial discretization at N_rho = 200,
a = 0.05 (sufficient radial mode count). Sweep N_phi in {24, 48, 96,
144, 192}. Compute K(t) at fine small-t grid; verify ratio approaches
1 as N_phi grows.

Exit gate
---------
Weyl ratio at t = 0.05 within 5% of 1.0 at fine N_phi: VERIFIED-UV.
Otherwise: characterize residual saturation, sharpen the diagnostic.
"""

import json
import time
from pathlib import Path

import numpy as np

OUT_JSON = Path(__file__).parent / "data" / "g4_3d_uv_extension.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def hermitian_radial_laplacian(N_rho, a, m):
    """Symmetric tridiagonal radial Hamiltonian from G4-3a-cleanup."""
    H = np.zeros((N_rho, N_rho))
    for i in range(N_rho):
        k = i + 1
        rho_k = k * a
        H[i, i] = 2.0 / a**2 + (m * m - 0.25) / rho_k**2
        if i > 0:
            H[i, i - 1] = -1.0 / a**2
        if i < N_rho - 1:
            H[i, i + 1] = -1.0 / a**2
    return H


def disk_eigenvalues(N_rho, a, N_phi):
    """Disk Laplacian eigenvalues with N_phi azimuthal modes."""
    eigenvalues = []
    h_phi = 2 * np.pi / N_phi
    for m_idx in range(N_phi):
        if m_idx <= N_phi // 2:
            m = m_idx
        else:
            m = m_idx - N_phi
        m_eff_sq = (2.0 / h_phi)**2 * np.sin(np.pi * m / N_phi)**2
        m_eff = np.sqrt(m_eff_sq)
        H_rad = hermitian_radial_laplacian(N_rho, a, m_eff)
        evals = np.linalg.eigvalsh(H_rad)
        eigenvalues.extend(evals.tolist())
    return np.array(sorted(eigenvalues))


def main():
    results = {}
    print("=" * 72)
    print("Sprint G4-3d-UV: N_phi refinement test")
    print("=" * 72)

    # Fixed substrate: R = 10, N_rho = 200, a = 0.05
    R_target = 10.0
    N_rho = 200
    a = 0.05
    A_disk = np.pi * R_target**2

    N_phi_values = [24, 48, 96, 144, 192]
    t_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    print(f"\n[Setup] R = {R_target}, N_rho = {N_rho}, a = {a}")
    print(f"        A_disk = pi R^2 = {A_disk:.4f}")
    print(f"        N_phi sweep: {N_phi_values}")
    print(f"        t values: {t_values}")

    results["setup"] = {
        "R": R_target, "N_rho": N_rho, "a": a, "A_disk": A_disk,
        "N_phi_values": N_phi_values, "t_values": t_values,
    }

    # Prediction from G4-3d memo Sec 5:
    # At t = 0.01, ratio at N_phi=24 plateaus at ~0.25; need N_phi ~ 96+
    print("\n[Prediction] G4-3d memo Sec 5: ratio plateau ~ N_phi^used / N_phi^full at small t")
    print(f"  At t=0.01, predicted N_phi >~ 96 to reach Weyl law")

    # -----------------------------------------------------------------------
    # Step 1: Compute eigenvalues at each N_phi
    # -----------------------------------------------------------------------
    print("\n[Step 1] Computing disk eigenvalues at each N_phi...")

    K_table = {}
    eigenvalue_stats = {}
    for N_phi in N_phi_values:
        t_start = time.time()
        evals = disk_eigenvalues(N_rho, a, N_phi)
        elapsed = time.time() - t_start
        n_modes = len(evals)
        lambda_max = float(evals[-1])
        lambda_min = float(evals[0])
        print(f"  N_phi = {N_phi:>3d}: {n_modes} modes, "
              f"lambda in [{lambda_min:.4f}, {lambda_max:.2f}], "
              f"{elapsed:.2f}s")

        K_vals = {t: float(np.sum(np.exp(-evals * t))) for t in t_values}
        K_table[N_phi] = K_vals
        eigenvalue_stats[N_phi] = {
            "n_modes": n_modes,
            "lambda_min": lambda_min,
            "lambda_max": lambda_max,
            "elapsed_s": float(elapsed),
        }

    results["eigenvalue_stats"] = eigenvalue_stats

    # -----------------------------------------------------------------------
    # Step 2: Weyl ratio table
    # -----------------------------------------------------------------------
    print("\n[Step 2] Weyl ratio K(t) * 4 pi t / A_disk:")
    print()
    header = f"  {'N_phi':>6}  " + "  ".join([f"t={t}".rjust(10) for t in t_values])
    print(header)
    print("  " + "-" * (10 + 12 * len(t_values)))

    ratio_table = {}
    for N_phi in N_phi_values:
        ratios = {t: K_table[N_phi][t] * 4 * np.pi * t / A_disk for t in t_values}
        ratio_table[N_phi] = ratios
        row = f"  {N_phi:>6d}  " + "  ".join([f"{ratios[t]:>10.4f}" for t in t_values])
        print(row)

    results["weyl_ratios"] = {str(k): v for k, v in ratio_table.items()}

    # -----------------------------------------------------------------------
    # Step 3: UV convergence — focus on small t
    # -----------------------------------------------------------------------
    print("\n[Step 3] UV convergence as N_phi grows:")
    uv_focus_t = [0.01, 0.02, 0.05]
    for t in uv_focus_t:
        print(f"\n  t = {t}:")
        print(f"    {'N_phi':>6}  {'K':>12}  {'K_Weyl':>12}  {'Ratio':>10}  {'Deficit':>10}")
        print("    " + "-" * 58)
        K_weyl = A_disk / (4 * np.pi * t)
        for N_phi in N_phi_values:
            K_disc = K_table[N_phi][t]
            ratio = K_disc / K_weyl
            deficit = 1.0 - ratio
            print(f"    {N_phi:>6d}  {K_disc:>12.4f}  {K_weyl:>12.4f}  "
                  f"{ratio:>10.4f}  {deficit:>10.4f}")

    results["uv_convergence_focus"] = {
        str(t): {
            "K_Weyl": A_disk / (4 * np.pi * t),
            "K_discrete": {str(N): K_table[N][t] for N in N_phi_values},
            "ratio": {str(N): K_table[N][t] / (A_disk / (4 * np.pi * t))
                      for N in N_phi_values},
        }
        for t in uv_focus_t
    }

    # -----------------------------------------------------------------------
    # Step 4: Convergence rate at t = 0.05
    # -----------------------------------------------------------------------
    print("\n[Step 4] Convergence rate at t = 0.05 (Weyl ratio vs N_phi):")
    t_anchor = 0.05
    ratios_at_anchor = np.array([ratio_table[N][t_anchor] for N in N_phi_values])
    N_phi_arr = np.array(N_phi_values)

    # Deficit (1 - ratio) should decay; fit log|deficit| vs log(N_phi)
    deficits = 1.0 - ratios_at_anchor
    if np.all(deficits > 0):
        log_N = np.log(N_phi_arr)
        log_def = np.log(deficits)
        rate = np.polyfit(log_N, log_def, 1)[0]
        print(f"  Deficit at t={t_anchor}: {deficits.tolist()}")
        print(f"  Linear fit log(deficit) vs log(N_phi): slope = {rate:.3f}")
        print(f"  (slope ~ -1: O(1/N_phi); slope ~ -2: O(1/N_phi^2))")
        results["uv_convergence_rate"] = float(rate)
    else:
        print(f"  Deficit at t={t_anchor}: {deficits.tolist()}")
        print(f"  Some negative deficits (overshoot); skipping fit.")
        results["uv_convergence_rate"] = None

    # -----------------------------------------------------------------------
    # Step 5: Verdict
    # -----------------------------------------------------------------------
    print("\n[Step 5] G4-3d-UV verdict:")
    print()
    finest_N_phi = N_phi_values[-1]
    ratio_t005 = ratio_table[finest_N_phi][0.05]
    ratio_t001 = ratio_table[finest_N_phi][0.01]

    if abs(ratio_t005 - 1.0) < 0.05:
        verdict = "POSITIVE-G4-3d-UV-VERIFIED"
        msg = "Weyl ratio at t=0.05 within 5% at fine N_phi."
    elif abs(ratio_t005 - 1.0) < 0.10:
        verdict = "POSITIVE-G4-3d-UV-CONVERGING-10PCT"
        msg = "Weyl ratio at t=0.05 within 10% at fine N_phi; converging."
    else:
        verdict = "PARTIAL-G4-3d-UV-CHARACTERIZED"
        msg = "Weyl ratio at t=0.05 not within 10%; saturation characterized."

    print(f"  At N_phi = {finest_N_phi}:")
    print(f"    Weyl ratio at t = 0.05: {ratio_t005:.4f}")
    print(f"    Weyl ratio at t = 0.01: {ratio_t001:.4f}")
    print()
    print(f"  Verdict: {verdict}")
    print(f"  {msg}")
    print()

    # G4-3d original prediction check
    if N_phi_values[2] == 96:
        ratio_96_t001 = ratio_table[96][0.01]
        ratio_24_t001 = ratio_table[24][0.01]
        gain = ratio_96_t001 / ratio_24_t001
        print(f"  Original G4-3d Sec 5 prediction check:")
        print(f"    Ratio at N_phi=24, t=0.01: {ratio_24_t001:.4f}")
        print(f"    Ratio at N_phi=96, t=0.01: {ratio_96_t001:.4f}")
        print(f"    Gain (N_phi=96 / N_phi=24): {gain:.2f}x")
        print(f"    (Memo predicted ~4x improvement with N_phi=96)")
        results["prediction_check"] = {
            "ratio_N24_t001": ratio_24_t001,
            "ratio_N96_t001": ratio_96_t001,
            "gain_96_over_24": gain,
        }

    results["verdict"] = verdict

    # Save
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
