"""Sprint G4-3c-proper — Wedge-lattice conical-defect test.

Closes the G4-3c open item: proper wedge lattice with apex angle
2*pi*alpha at fixed angular resolution h_phi = 2*pi/N_0 (so that
N_phi = alpha * N_0). Periodic BC at the seam (closes the wedge
into a cone). The original G4-3c naive sweep varied N_phi at fixed
period 2*pi (resolution sweep, not apex-angle sweep) and could not
extract the Sommerfeld/Cheeger (1/12)(1/alpha - alpha) coefficient.

Sommerfeld/Cheeger formula (continuum):
    K_cone(alpha, t) = alpha * K_disk(t) + (1/12)(1/alpha - alpha) + O(t)

Discrete test:
    Delta K(alpha, t) := K_wedge_disc(alpha, t) - alpha * K_disk_disc(t)
    -> (1/12)(1/alpha - alpha) as t -> 0 (continuum limit)

Setup
-----
N_0 = 120 (divisible by 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, ...).
Alpha sweep: {1/3, 1/2, 2/3, 1, 3/2, 2, 3}. Gives N_phi values
{40, 60, 80, 120, 180, 240, 360}.
Radial: N_rho = 200, a = 0.05 (matches T2 substrate).
"""

import json
from pathlib import Path

import numpy as np

OUT_JSON = Path(__file__).parent / "data" / "g4_3c_proper_wedge.json"
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


def wedge_eigenvalues(N_rho, a, N_phi, alpha):
    """Wedge lattice with apex angle 2*pi*alpha, periodic BC.

    The radial Hamiltonian uses centrifugal m_eff from the discrete
    azimuthal Laplacian. With period 2*pi*alpha and N_phi sites:
        h_phi = 2*pi*alpha / N_phi
        m_eff_sq[k] = (2/h_phi)^2 * sin^2(pi*k/N_phi)

    In the small-k continuum limit: m_eff_sq -> (k/alpha)^2 matching
    the continuum wedge prediction m = k/alpha for k integer.
    """
    eigenvalues = []
    h_phi = 2 * np.pi * alpha / N_phi
    for k_idx in range(N_phi):
        if k_idx <= N_phi // 2:
            k = k_idx
        else:
            k = k_idx - N_phi
        m_eff_sq = (2.0 / h_phi)**2 * np.sin(np.pi * k / N_phi)**2
        m_eff = np.sqrt(m_eff_sq)
        H_rad = hermitian_radial_laplacian(N_rho, a, m_eff)
        evals = np.linalg.eigvalsh(H_rad)
        eigenvalues.extend(evals.tolist())
    return np.array(sorted(eigenvalues))


def main():
    results = {}
    print("=" * 72)
    print("Sprint G4-3c-proper: Wedge-lattice conical-defect test")
    print("=" * 72)

    N_rho = 200
    a = 0.05
    R = N_rho * a
    N_0 = 120

    alpha_pairs = [
        ("1/3", 1.0/3.0, 40),
        ("1/2", 0.5,      60),
        ("2/3", 2.0/3.0, 80),
        ("1",   1.0,      120),
        ("3/2", 1.5,      180),
        ("2",   2.0,      240),
        ("3",   3.0,      360),
    ]

    t_values = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    print(f"\n[Setup] R = {R}, N_rho = {N_rho}, a = {a}, N_0 = {N_0}")
    print(f"        h_phi (fixed across alpha) = 2*pi/N_0 = {2*np.pi/N_0:.5f}")
    print(f"        Alpha sweep: {[name for name, _, _ in alpha_pairs]}")
    print(f"        N_phi(alpha) = alpha * N_0: {[Nphi for _, _, Nphi in alpha_pairs]}")
    print(f"        t values: {t_values}")
    results["setup"] = {
        "R": R, "N_rho": N_rho, "a": a, "N_0": N_0,
        "alpha_values": [(name, val, Nphi) for name, val, Nphi in alpha_pairs],
        "t_values": t_values,
    }

    # ------------------------------------------------------------------
    # Step 1: Compute K_wedge(alpha, t) at each alpha
    # ------------------------------------------------------------------
    print("\n[Step 1] Computing wedge heat traces...")
    K_table = {}
    for name, alpha, N_phi in alpha_pairs:
        evals = wedge_eigenvalues(N_rho, a, N_phi, alpha)
        K_vals = {t: float(np.sum(np.exp(-evals * t))) for t in t_values}
        K_table[name] = (alpha, N_phi, K_vals)
        print(f"  alpha = {name:>4}  N_phi = {N_phi:>3}  "
              f"{len(evals)} modes, lambda_max = {evals[-1]:.2f}")

    # Sanity: alpha = 1 should give the standard disk heat trace
    K_disk = K_table["1"][2]
    print(f"\n  Reference (alpha = 1, standard disk):")
    print(f"  " + "  ".join([f"K(t={t}) = {K_disk[t]:.4f}" for t in t_values]))

    # ------------------------------------------------------------------
    # Step 2: Residual Delta K(alpha, t) = K_wedge(alpha) - alpha * K_disk
    # ------------------------------------------------------------------
    print("\n[Step 2] Topological residual Delta K(alpha, t) = K_wedge - alpha * K_disk:")
    print()
    print(f"  Sommerfeld/Cheeger prediction at small t: (1/12)(1/alpha - alpha)")
    print()
    sc_predictions = {name: (1.0/12.0) * (1.0/alpha - alpha)
                     for name, alpha, _ in alpha_pairs}
    print(f"  Predicted residual:")
    for name, _, _ in alpha_pairs:
        print(f"    alpha = {name:>4}: (1/12)(1/alpha - alpha) = {sc_predictions[name]:+.6f}")

    print()
    print(f"  Computed residual Delta K(alpha, t):")
    header = f"  {'alpha':>6}  {'SC pred':>10}  " + "  ".join([f"t={t}".rjust(11) for t in t_values])
    print(header)
    print("  " + "-" * (20 + 13 * len(t_values)))

    delta_K_table = {}
    for name, alpha, N_phi in alpha_pairs:
        K_vals = K_table[name][2]
        delta = {t: K_vals[t] - alpha * K_disk[t] for t in t_values}
        delta_K_table[name] = (alpha, delta)
        row_data = "  ".join([f"{delta[t]:+11.4f}" for t in t_values])
        print(f"  {name:>6}  {sc_predictions[name]:+10.5f}  {row_data}")

    results["K_wedge"] = {
        name: {"alpha": alpha, "N_phi": N_phi, "K": K_vals}
        for name, (alpha, N_phi, K_vals) in K_table.items()
    }
    results["delta_K"] = {
        name: {"alpha": alpha, "SC_prediction": sc_predictions[name], "delta_K": delta}
        for name, (alpha, delta) in delta_K_table.items()
    }

    # ------------------------------------------------------------------
    # Step 3: Slope test - fit Delta K(alpha, t_small) vs (1/alpha - alpha)
    # ------------------------------------------------------------------
    print("\n[Step 3] Slope test: fit Delta K vs (1/alpha - alpha) at each t.")
    print(f"         Predicted slope = 1/12 = {1.0/12.0:.6f}")
    print()
    sc_x = np.array([(1.0/alpha - alpha) for _, alpha, _ in alpha_pairs])
    print(f"  (1/alpha - alpha) values: {sc_x.tolist()}")
    print()

    slope_table = {}
    for t in t_values:
        y = np.array([delta_K_table[name][1][t] for name, _, _ in alpha_pairs])
        # Linear fit through origin (since at (1/alpha-alpha)=0 i.e. alpha=1, delta_K=0)
        slope = np.sum(sc_x * y) / np.sum(sc_x * sc_x)
        # Also do unconstrained linear fit
        slope_full, intercept = np.polyfit(sc_x, y, 1)
        # Residuals
        y_pred = slope * sc_x
        resid_rms = float(np.sqrt(np.mean((y - y_pred)**2)))
        slope_table[t] = {
            "slope_through_origin": float(slope),
            "slope_unconstrained": float(slope_full),
            "intercept_unconstrained": float(intercept),
            "resid_rms_origin": resid_rms,
        }
        print(f"  t = {t:>5}:  slope (origin) = {slope:+.6f}  "
              f"slope (free) = {slope_full:+.6f}  intercept = {intercept:+.4f}  "
              f"resid RMS = {resid_rms:.4f}")

    results["slope_test"] = slope_table

    # ------------------------------------------------------------------
    # Step 4: t-extrapolation of slope
    # ------------------------------------------------------------------
    print("\n[Step 4] t-extrapolation: how does the slope behave as t -> 0?")
    print()
    print(f"  Target (Sommerfeld/Cheeger): slope -> 1/12 = {1.0/12.0:.6f} as t -> 0")
    print(f"  At t -> infinity the slope reflects spectral structure, not the tip.")

    # Show slope trajectory
    print(f"\n  Slope trajectory:")
    for t in t_values:
        print(f"    t = {t:>5}: slope (origin) = {slope_table[t]['slope_through_origin']:+.6f}")

    # Test ratio: alpha = 1/n + alpha = n should give zero delta_K sum
    print(f"\n[Step 5] Reciprocal cancellation test:")
    print(f"  At each t: Delta K(1/n) + Delta K(n) should = 0 (topological cancellation)")
    print()
    reciprocal_pairs = [("1/3", "3"), ("1/2", "2"), ("2/3", "3/2")]
    recip_table = {}
    for name_inv, name_n in reciprocal_pairs:
        print(f"  Pair (alpha = {name_inv}, alpha = {name_n}):")
        recip_table[f"{name_inv}+{name_n}"] = {}
        for t in t_values:
            d_inv = delta_K_table[name_inv][1][t]
            d_n = delta_K_table[name_n][1][t]
            total = d_inv + d_n
            recip_table[f"{name_inv}+{name_n}"][str(t)] = float(total)
            print(f"    t = {t:>5}: Delta K({name_inv}) + Delta K({name_n}) = "
                  f"{d_inv:+8.4f} + {d_n:+8.4f} = {total:+8.4f}")
        print()
    results["reciprocal_cancellation"] = recip_table

    # ------------------------------------------------------------------
    # Step 6: Verdict
    # ------------------------------------------------------------------
    print("\n[Step 6] Verdict:")
    print()
    # Take slope at smallest t as the conical defect signature
    slope_smallest_t = slope_table[t_values[0]]["slope_through_origin"]
    rel_err = (slope_smallest_t - 1.0/12.0) / (1.0/12.0)
    print(f"  Slope at t = {t_values[0]}: {slope_smallest_t:+.6f}")
    print(f"  Sommerfeld/Cheeger target:  {1.0/12.0:+.6f}")
    print(f"  Relative error: {rel_err*100:+.2f}%")
    print()

    if abs(rel_err) < 0.20:
        verdict = "POSITIVE-G4-3c-PROPER-VERIFIED"
        msg = "Sommerfeld/Cheeger 1/12 coefficient extracted within 20%."
    elif abs(rel_err) < 0.50:
        verdict = "PARTIAL-G4-3c-PROPER-SIGN-CORRECT"
        msg = "Slope sign and OoM correct; quantitative match needs finer extraction."
    else:
        verdict = "NEGATIVE-G4-3c-PROPER"
        msg = "Wedge-lattice + naive K_disk subtraction does not give Sommerfeld/Cheeger at this scale."

    print(f"  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["slope_at_smallest_t"] = float(slope_smallest_t)
    results["sc_target"] = 1.0/12.0
    results["rel_err"] = float(rel_err)

    # Save
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
