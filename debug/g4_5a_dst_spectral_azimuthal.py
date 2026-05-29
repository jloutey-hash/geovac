"""Sprint G4-5a DST -- Spectral azimuthal discretization for UV tip recovery.

Replaces the FD azimuthal Laplacian discretization in DiscreteDiskDirac /
DiscreteWedgeDirac with an EXACT spectral discretization using the actual
azimuthal eigenvalues (integer for periodic, half-integer for fermionic
anti-periodic) instead of the FD form

    lambda_k^FD = (2/h_phi)^2 * sin^2(pi (k + 0.5) / N_phi).

The FD form has a UV overshoot at the truncation edge: at k ~ N_phi/2,
lambda_k^FD / lambda_k^continuum -> 4/pi^2 ~ 0.405 (T2 G4-3d-UV finding).
This produces undershoot of the high-m azimuthal Laplacian eigenvalues,
which dominates the small-t heat-trace contribution where the tip term
lives.

The cure: use exact azimuthal eigenvalues m_eff = m where
  - disk anti-periodic: m_eff(k) = k + 0.5, k in symmetric range
  - wedge anti-periodic (period 2*pi*alpha): m_eff(k) = (k + 0.5) / alpha
  - periodic scalar (for sanity check): m_eff(k) = k

Method
------
1. Implement spectral driver functions disk_dirac_heat_trace_spectral and
   wedge_dirac_heat_trace_spectral that use exact m_eff (no FD sin^2).
2. F6 sanity: at alpha = 1, spectral wedge equals spectral disk
   eigenvalues bit-exactly.
3. Per-t tip recovery diagnostic: at substrate panel (R=10, a=0.05,
   N_rho=200, N_0=120, eps=12/120=0.1), compute Delta'_spectral(t) at the
   G4-5a-refined t-grid {0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
   1, 2, 5, 10}. At t = 0.0025, G4-5a-refined FD gave 1.3% recovery;
   spectral should give much higher (>50%, ideally >90%).
4. Integrate with Gaussian cutoff and compute S_tip(Lambda) at
   Lambda in {0.5, 1, 1.5, 2}.

Decision gate
-------------
- POSITIVE: per-t tip recovery at t=0.0025 > 90%; integrated S_tip ratio
  > 0.9 at all Lambda in {0.5, 1, 2}.
- PARTIAL: per-t recovery > 50% at t=0.0025; or integration ratio > 0.7.
- NEGATIVE: spectral azimuthal does not significantly improve UV
  recovery.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

import numpy as np
from scipy.special import exp1

OUT_JSON = Path(__file__).parent / "data" / "g4_5a_dst_spectral_azimuthal.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# =================================================================
# Spectral azimuthal m_eff arrays
# =================================================================


def _m_eff_spectral_disk_antiperiodic(N_phi: int) -> np.ndarray:
    """Exact m_eff array for disk anti-periodic BC.

    Fermionic anti-periodic on [0, 2*pi): m takes half-integer values
    m in {1/2, -1/2, 3/2, -3/2, ...}. Truncated at N_phi modes (the
    same count as the FD discretization for matched comparison),
    indexed symmetrically:
      k_idx in {0, 1, ..., N_phi - 1}
      k = k_idx if k_idx <= N_phi // 2 else k_idx - N_phi
      m_eff(k) = k + 0.5

    Note: this uses the SAME symmetric k mapping as DiscreteDiskDirac
    in geovac/gravity/warped_dirac.py.
    """
    k_arr = np.array(
        [k_idx if k_idx <= N_phi // 2 else k_idx - N_phi
         for k_idx in range(N_phi)],
        dtype=float,
    )
    return k_arr + 0.5


def _m_eff_spectral_disk_periodic(N_phi: int) -> np.ndarray:
    """Exact m_eff array for disk periodic BC (scalar sanity check).

    Bosonic periodic: m takes integer values m in {0, 1, -1, 2, -2, ...}.
    Truncated at N_phi modes with the same symmetric k mapping.
    """
    k_arr = np.array(
        [k_idx if k_idx <= N_phi // 2 else k_idx - N_phi
         for k_idx in range(N_phi)],
        dtype=float,
    )
    return k_arr  # integer m


def _m_eff_spectral_wedge_antiperiodic(N_phi: int, alpha: float) -> np.ndarray:
    """Exact m_eff array for wedge anti-periodic BC with apex 2*pi*alpha.

    Fermionic anti-periodic on [0, 2*pi*alpha):
        phi-eigenvalues in continuum:  m = (k + 1/2) / alpha for integer k.

    Same symmetric k mapping as DiscreteWedgeDirac.

    At alpha = 1, this reduces bit-exactly to
    _m_eff_spectral_disk_antiperiodic(N_phi).
    """
    k_arr = np.array(
        [k_idx if k_idx <= N_phi // 2 else k_idx - N_phi
         for k_idx in range(N_phi)],
        dtype=float,
    )
    return (k_arr + 0.5) / alpha


# =================================================================
# Hermitian polar radial Laplacian (same as geovac.gravity.warped_dirac)
# =================================================================


def _hermitian_radial_laplacian(N_rho: int, a: float, m_eff: float) -> np.ndarray:
    """Symmetric tridiagonal radial Laplacian on L^2(rho drho).

    Uses u = sqrt(rho) f substitution (G4-3a-cleanup convention). The
    centrifugal term becomes (m_eff^2 - 1/4) / rho^2 in the
    u-representation.

    This MATCHES geovac/gravity/warped_dirac.py
    DiscreteDiskDirac._hermitian_radial_laplacian and
    DiscreteWedgeDirac._hermitian_radial_laplacian exactly.
    """
    H = np.zeros((N_rho, N_rho))
    for i in range(N_rho):
        k = i + 1
        rho = k * a
        H[i, i] = 2.0 / a ** 2 + (m_eff ** 2 - 0.25) / rho ** 2
        if i > 0:
            H[i, i - 1] = -1.0 / a ** 2
        if i < N_rho - 1:
            H[i, i + 1] = -1.0 / a ** 2
    return H


# =================================================================
# Spectral heat-trace drivers
# =================================================================


def disk_dirac_eigenvalues_spectral(
    N_rho: int, a: float, N_phi: int,
) -> np.ndarray:
    """All eigenvalues of D_{D^2}^2 with spectral azimuthal discretization.

    Replaces the FD eigenvalues
        m_eff_sq = (2/h_phi)^2 * sin^2(pi (k+0.5) / N_phi)
    with the EXACT spectral values m_eff = k + 0.5 for the anti-periodic
    bosonic case.

    The radial Laplacian is the SAME hermitian polar form. Rank-2 spinor
    bundle: each scalar eigenvalue is doubled.

    Returns sorted array of length 2 * N_rho * N_phi.
    """
    m_eff_arr = _m_eff_spectral_disk_antiperiodic(N_phi)
    scalar_eigs: List[float] = []
    for m_eff in m_eff_arr:
        H = _hermitian_radial_laplacian(N_rho, a, float(m_eff))
        evals = np.linalg.eigvalsh(H)
        scalar_eigs.extend(evals.tolist())
    scalar_arr = np.array(scalar_eigs)
    return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))


def disk_dirac_heat_trace_spectral(
    N_rho: int, a: float, N_phi: int, t: float,
) -> float:
    """Tr exp(-t D^2) for disk-Dirac with spectral azimuthal."""
    eigs = disk_dirac_eigenvalues_spectral(N_rho, a, N_phi)
    return float(np.sum(np.exp(-eigs * t)))


def wedge_dirac_eigenvalues_spectral(
    N_rho: int, a: float, N_phi: int, alpha: float,
) -> np.ndarray:
    """All eigenvalues of D^2 on wedge with spectral azimuthal.

    Uses m_eff(k) = (k + 0.5) / alpha (exact spectral) instead of the
    FD sin^2 form. At alpha = 1, this reduces bit-exactly to
    disk_dirac_eigenvalues_spectral.

    Returns sorted array of length 2 * N_rho * N_phi.
    """
    m_eff_arr = _m_eff_spectral_wedge_antiperiodic(N_phi, alpha)
    scalar_eigs: List[float] = []
    for m_eff in m_eff_arr:
        H = _hermitian_radial_laplacian(N_rho, a, float(m_eff))
        evals = np.linalg.eigvalsh(H)
        scalar_eigs.extend(evals.tolist())
    scalar_arr = np.array(scalar_eigs)
    return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))


def wedge_dirac_heat_trace_spectral(
    N_rho: int, a: float, N_phi: int, alpha: float, t: float,
) -> float:
    """Tr exp(-t D^2) for wedge-Dirac with spectral azimuthal."""
    eigs = wedge_dirac_eigenvalues_spectral(N_rho, a, N_phi, alpha)
    return float(np.sum(np.exp(-eigs * t)))


# =================================================================
# F6 sanity: alpha = 1 wedge = disk
# =================================================================


def verify_F6_spectral_alpha1_match(
    N_rho: int, a: float, N_phi: int, tol: float = 1e-13,
) -> dict:
    """At alpha = 1, spectral wedge eigenvalues = spectral disk eigenvalues."""
    disk_eigs = disk_dirac_eigenvalues_spectral(N_rho, a, N_phi)
    wedge_eigs = wedge_dirac_eigenvalues_spectral(N_rho, a, N_phi, alpha=1.0)
    max_diff = float(np.max(np.abs(disk_eigs - wedge_eigs)))
    return {
        "max_diff": max_diff,
        "passed": bool(max_diff < tol),
        "N_eigenvalues": len(disk_eigs),
        "tolerance": tol,
    }


# =================================================================
# Gaussian Mellin moment
# =================================================================


def gaussian_mellin_E1(Lambda: float, t_min: float, t_max: float) -> float:
    """Exact Gaussian Mellin moment on [t_min, t_max].

    M_0(Lambda; t_min, t_max) = int_{t_min}^{t_max} (dt/t) exp(-t Lambda^2)
                              = E_1(t_min Lambda^2) - E_1(t_max Lambda^2)
    """
    return float(exp1(t_min * Lambda ** 2) - exp1(t_max * Lambda ** 2))


# =================================================================
# Main
# =================================================================


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-5a DST -- Spectral azimuthal discretization (UV cure)")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Step 1: Substrate setup (matches G4-5a-refined panel)
    # ------------------------------------------------------------------
    R = 10.0
    a = 0.05
    N_rho = 200
    N_0 = 120
    print(f"\n[Setup] R={R}, a={a}, N_rho={N_rho}, N_0={N_0}")
    print(f"  Substrate UV cutoff: t_UV ~ a^2 = {a ** 2}")
    print(f"  Substrate IR cutoff: t_IR ~ R^2 = {R ** 2}")

    # Replica eps: matches G4-4f best window
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps = (alpha_plus - alpha_minus) / 2
    print("\n[Replica parameters]")
    print(f"  alpha_+ = {alpha_plus}, alpha_- = {alpha_minus}, eps = {eps}")
    print(f"  N_plus = {N_plus}, N_minus = {N_minus}")

    results["setup"] = {
        "R": R, "a": a, "N_rho": N_rho, "N_0": N_0,
        "k_step": k_step,
        "N_plus": N_plus, "N_minus": N_minus,
        "alpha_plus": alpha_plus, "alpha_minus": alpha_minus,
        "eps": eps,
    }

    # ------------------------------------------------------------------
    # Step 2: F6 sanity at sprint-scale params
    # ------------------------------------------------------------------
    print("\n[Step 2] F6 sanity: at alpha = 1, spectral wedge = spectral disk")
    # Use moderate panel for F6 sanity
    f6_N_rho = 50; f6_a = 0.1; f6_N_phi = 24
    f6 = verify_F6_spectral_alpha1_match(
        N_rho=f6_N_rho, a=f6_a, N_phi=f6_N_phi,
    )
    print(f"  Test panel: N_rho={f6_N_rho}, a={f6_a}, N_phi={f6_N_phi}")
    print(f"  max |disk_spectral - wedge_spectral(alpha=1)|: "
          f"{f6['max_diff']:.3e}")
    print(f"  Eigenvalues: {f6['N_eigenvalues']}")
    print(f"  F6 PASSED: {f6['passed']}")
    results["F6_sanity"] = f6

    # Print spectral m_eff range vs FD
    m_eff_spec = _m_eff_spectral_disk_antiperiodic(N_0)
    h_phi = 2 * np.pi / N_0
    m_eff_FD_sq = (2.0 / h_phi) ** 2 * np.sin(
        np.pi * (np.arange(N_0) + 0.5) / N_0
    ) ** 2
    # Use symmetric k mapping for FD
    k_arr = np.array(
        [k_idx if k_idx <= N_0 // 2 else k_idx - N_0 for k_idx in range(N_0)],
        dtype=float,
    )
    m_eff_FD_sq_sym = (2.0 / h_phi) ** 2 * np.sin(
        np.pi * (k_arr + 0.5) / N_0
    ) ** 2

    print("\n  Diagnostic: spectral vs FD m_eff^2 at N_0 = 120")
    print(f"    Lowest |m_eff^2|: spec={(m_eff_spec[0]) ** 2:.6f}, "
          f"FD={m_eff_FD_sq_sym[0]:.6f}, "
          f"ratio={m_eff_FD_sq_sym[0] / (m_eff_spec[0]) ** 2:.6f}")
    # find max-|k| modes
    idx_high = int(np.argmax(np.abs(m_eff_spec)))
    print(f"    Highest |m_eff^2|: spec={(m_eff_spec[idx_high]) ** 2:.6f}, "
          f"FD={m_eff_FD_sq_sym[idx_high]:.6f}, "
          f"ratio={m_eff_FD_sq_sym[idx_high] / (m_eff_spec[idx_high]) ** 2:.6f}")
    print(f"    (Predicted ratio at edge: 4/pi^2 = {4 / np.pi ** 2:.6f})")

    results["spectral_vs_FD_meff"] = {
        "N_phi": N_0,
        "lowest_spec_sq": float(m_eff_spec[0] ** 2),
        "lowest_FD_sq": float(m_eff_FD_sq_sym[0]),
        "lowest_ratio_FD_over_spec": float(
            m_eff_FD_sq_sym[0] / m_eff_spec[0] ** 2
        ),
        "highest_spec_sq": float(m_eff_spec[idx_high] ** 2),
        "highest_FD_sq": float(m_eff_FD_sq_sym[idx_high]),
        "highest_ratio_FD_over_spec": float(
            m_eff_FD_sq_sym[idx_high] / m_eff_spec[idx_high] ** 2
        ),
        "predicted_edge_ratio_4_over_pi2": float(4 / np.pi ** 2),
    }

    # ------------------------------------------------------------------
    # Step 3: Per-t tip recovery diagnostic (the headline)
    # ------------------------------------------------------------------
    t_grid = np.array([0.0025, 0.005, 0.01, 0.02, 0.05, 0.1,
                       0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    print(f"\n[Step 3] Per-t tip recovery diagnostic (12-point t-grid)")
    print(f"  t-grid: {t_grid.tolist()}")
    print(f"  Computing K(t) at each t via spectral azimuthal...")
    print(f"  (Total: 3 panels x 12 t-points = 36 heat traces)")
    print()

    K_disk_arr = np.array([
        disk_dirac_heat_trace_spectral(N_rho, a, N_0, float(t))
        for t in t_grid
    ])
    K_plus_arr = np.array([
        wedge_dirac_heat_trace_spectral(N_rho, a, N_plus, alpha_plus, float(t))
        for t in t_grid
    ])
    K_minus_arr = np.array([
        wedge_dirac_heat_trace_spectral(N_rho, a, N_minus, alpha_minus, float(t))
        for t in t_grid
    ])

    dK_dalpha = (K_plus_arr - K_minus_arr) / (alpha_plus - alpha_minus)
    tip_term = dK_dalpha - K_disk_arr

    # G4-5a-refined FD baseline (taken from existing JSON)
    g4_5a_refined_FD = {
        0.0025: 0.013142053245246643,
        0.005:  0.1653271082086576,
        0.01:   0.362258103385102,
        0.02:   0.5187012123651584,
        0.05:   0.6759082480057259,
        0.1:    0.7634,
        0.2:    0.8293,
        0.5:    0.8907,
        1.0:    0.9229,
        2.0:    0.9456,
        5.0:    0.9670,
        10.0:   0.9770,
    }

    print(f"  {'t':>8}  {'K_disk':>13}  {'dK/dalpha':>13}  {'tip term':>12}  "
          f"{'spec rec':>10}  {'FD rec (G4-5a-ref)':>20}")
    print("  " + "-" * 86)
    per_t_data = []
    for i, t in enumerate(t_grid):
        spec_rec = tip_term[i] / (1.0 / 6.0)
        FD_rec = g4_5a_refined_FD.get(float(t), float("nan"))
        marker = " *UV*" if t < 0.1 else ""
        print(f"  {t:>8.4f}  {K_disk_arr[i]:>13.4f}  {dK_dalpha[i]:>13.4f}  "
              f"{tip_term[i]:>+12.6f}  {spec_rec * 100:>9.2f}%  "
              f"{FD_rec * 100:>18.2f}%{marker}")
        per_t_data.append({
            "t": float(t),
            "K_disk": float(K_disk_arr[i]),
            "dK_dalpha": float(dK_dalpha[i]),
            "tip_term": float(tip_term[i]),
            "recovery_spec": float(spec_rec),
            "recovery_FD_baseline": float(FD_rec),
        })

    results["per_t_diagnostic"] = per_t_data
    results["t_grid"] = t_grid.tolist()
    results["K_disk_spec"] = K_disk_arr.tolist()
    results["K_plus_spec"] = K_plus_arr.tolist()
    results["K_minus_spec"] = K_minus_arr.tolist()
    results["dK_dalpha_spec"] = dK_dalpha.tolist()
    results["tip_term_spec"] = tip_term.tolist()
    results["recovery_spec"] = (tip_term / (1.0 / 6.0)).tolist()
    results["recovery_FD_baseline"] = [
        g4_5a_refined_FD.get(float(t), float("nan")) for t in t_grid
    ]

    # Key UV diagnostic value
    rec_at_uv = tip_term[0] / (1.0 / 6.0)
    print()
    print(f"  HEADLINE: spectral tip recovery at t = a^2 = 0.0025: "
          f"{rec_at_uv * 100:.2f}%")
    print(f"  G4-5a-refined FD baseline at same t: "
          f"{g4_5a_refined_FD[0.0025] * 100:.2f}%")
    print(f"  Improvement: x{rec_at_uv / g4_5a_refined_FD[0.0025]:.1f}")

    results["uv_headline"] = {
        "t_uv": float(t_grid[0]),
        "spec_recovery": float(rec_at_uv),
        "FD_recovery_baseline": float(g4_5a_refined_FD[0.0025]),
        "improvement_factor": float(
            rec_at_uv / g4_5a_refined_FD[0.0025]
        ),
    }

    # ------------------------------------------------------------------
    # Step 4: Integrate over t with Gaussian cutoff
    # ------------------------------------------------------------------
    print("\n[Step 4] Integrate tip term with Gaussian cutoff")
    print()
    Lambda_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    log_t = np.log(t_grid)

    # G4-5a-refined baseline (exact Mellin) for comparison
    g4_5a_refined_exact = {
        0.5: 0.6345, 1.0: 0.5717, 1.5: 0.5218, 2.0: 0.4835,
        3.0: 0.4237, 5.0: 0.3358,
    }

    print(f"  {'Lambda':>8}  {'J integral':>14}  {'S_tip (spec)':>14}  "
          f"{'M_0 exact':>11}  {'S_tip pred':>11}  "
          f"{'spec/pred ratio':>15}  {'FD ratio (G4-5a-ref)':>21}")
    print("  " + "-" * 100)

    integrate_results = {}
    for Lambda in Lambda_values:
        f_cutoff = np.exp(-t_grid * Lambda ** 2)
        integrand_log = f_cutoff * tip_term
        J = float(np.trapezoid(integrand_log, log_t))
        S_tip_spec = +J / 2

        M_0_exact = gaussian_mellin_E1(Lambda, a ** 2, R ** 2)
        S_tip_pred = +(1 / 12) * M_0_exact

        ratio = S_tip_spec / S_tip_pred if abs(S_tip_pred) > 1e-12 else float("inf")
        ratio_FD = g4_5a_refined_exact.get(Lambda, float("nan"))

        integrate_results[Lambda] = {
            "J": J,
            "S_tip_spec": S_tip_spec,
            "M_0_exact": M_0_exact,
            "S_tip_pred": S_tip_pred,
            "ratio_spec": ratio,
            "ratio_FD_baseline": ratio_FD,
        }
        print(f"  {Lambda:>8.2f}  {J:>14.6e}  {S_tip_spec:>+14.6e}  "
              f"{M_0_exact:>11.6f}  {S_tip_pred:>+11.6f}  "
              f"{ratio:>15.4f}  {ratio_FD:>21.4f}")

    results["integrate_results"] = {
        str(k): v for k, v in integrate_results.items()
    }

    # ------------------------------------------------------------------
    # Step 5: Decision gate
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("\n[Step 5] Decision gate:")
    print()
    print("  POSITIVE: per-t recovery at t=0.0025 > 90% AND integrated ratio")
    print("           > 0.9 at all Lambda in {0.5, 1, 2}")
    print("  PARTIAL: per-t recovery at t=0.0025 > 50%, OR integrated ratio")
    print("           > 0.7 at all Lambda in {0.5, 1, 2}")
    print("  NEGATIVE: spectral does not significantly improve UV recovery")
    print()

    rec_at_uv = float(tip_term[0] / (1.0 / 6.0))
    target_Lambdas = [0.5, 1.0, 2.0]
    integ_ratios = [integrate_results[L]["ratio_spec"] for L in target_Lambdas]
    fd_ratios = [integrate_results[L]["ratio_FD_baseline"] for L in target_Lambdas]

    min_integ = float(min(integ_ratios))
    avg_integ_improvement = float(np.mean(
        [r - fd for r, fd in zip(integ_ratios, fd_ratios)]
    ))

    print(f"  t = 0.0025 spectral recovery: {rec_at_uv * 100:.2f}%")
    print(f"  Target Lambdas: {target_Lambdas}")
    print(f"  Spectral integrated ratios: {[f'{r:.4f}' for r in integ_ratios]}")
    print(f"  FD baseline (exact) ratios: {[f'{r:.4f}' for r in fd_ratios]}")
    print(f"  min(spec integrated ratio): {min_integ:.4f}")
    print(f"  avg improvement (spec - FD): {avg_integ_improvement:+.4f}")
    print()

    cond_positive_per_t = rec_at_uv > 0.9
    cond_positive_integ = all(r > 0.9 for r in integ_ratios)
    cond_partial_per_t = rec_at_uv > 0.5
    cond_partial_integ = all(r > 0.7 for r in integ_ratios)
    improves_all = all(
        integrate_results[L]["ratio_spec"] > integrate_results[L]["ratio_FD_baseline"]
        for L in target_Lambdas
    )

    print(f"  per-t at t=0.0025 > 90%? {cond_positive_per_t}")
    print(f"  per-t at t=0.0025 > 50%? {cond_partial_per_t}")
    print(f"  integrated ratio > 0.9 at all targets? {cond_positive_integ}")
    print(f"  integrated ratio > 0.7 at all targets? {cond_partial_integ}")
    print(f"  spectral improves over FD at all targets? {improves_all}")
    print()

    if cond_positive_per_t and cond_positive_integ:
        verdict = "POSITIVE-G4-5a-DST-SPECTRAL"
        msg = (
            "Spectral azimuthal cure works: UV per-t recovery > 90% AND "
            "integrated S_tip ratio > 0.9 at all in-range Lambdas. The "
            "T2 G4-3d-UV high-m angular truncation overshoot was indeed "
            "the residual barrier to full UV closure."
        )
    elif cond_partial_per_t or cond_partial_integ:
        verdict = "PARTIAL-G4-5a-DST-SPECTRAL"
        if cond_partial_per_t and cond_partial_integ:
            sub = "BOTH"
        elif cond_partial_per_t:
            sub = "PER-T-ONLY"
        else:
            sub = "INTEGRATED-ONLY"
        msg = (
            f"Spectral azimuthal substantially improves UV recovery vs the "
            f"FD G4-5a-refined baseline ({sub}). Per-t recovery at t=a^2 "
            f"climbs from 1.3% (FD) to {rec_at_uv * 100:.1f}% (spec). "
            f"min integrated ratio across {target_Lambdas} is "
            f"{min_integ:.4f}. Confirms T2 G4-3d-UV finding: the angular "
            f"truncation was the dominant UV barrier; residual gap likely "
            f"tied to remaining substrate effects (finite N_rho, radial "
            f"discretization, IR cutoff)."
        )
    elif improves_all:
        verdict = "WEAK-PARTIAL-G4-5a-DST-SPECTRAL"
        msg = (
            "Spectral azimuthal improves over FD at every in-range Lambda "
            "but does not meet PARTIAL gate (50% per-t or 70% integrated). "
            "T2 angular overshoot was a contributing barrier but not the "
            "dominant one; other substrate effects must be at work."
        )
    else:
        verdict = "NEGATIVE-G4-5a-DST-SPECTRAL"
        msg = (
            "Spectral azimuthal does not significantly improve UV recovery. "
            "T2 G4-3d-UV angular truncation was NOT the dominant barrier; "
            "the UV gap lies elsewhere (radial discretization, IR cutoff, "
            "or replica-method-intrinsic effect)."
        )

    print(f"  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["min_integrated_ratio_targets"] = min_integ
    results["avg_integ_improvement_vs_FD"] = avg_integ_improvement
    results["per_t_uv_recovery"] = rec_at_uv
    results["FD_per_t_uv_recovery"] = float(g4_5a_refined_FD[0.0025])
    results["cond_positive_per_t"] = cond_positive_per_t
    results["cond_positive_integ"] = cond_positive_integ
    results["cond_partial_per_t"] = cond_partial_per_t
    results["cond_partial_integ"] = cond_partial_integ
    results["improves_all_targets"] = improves_all

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
