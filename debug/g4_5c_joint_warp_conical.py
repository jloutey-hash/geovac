"""Sprint G4-5c -- Joint variable-warp + conical-defect Dirac (F11 closure of G4-5).

Driver-level construction (no production code mods). Combines:
  - Wedge azimuthal sector: anti-periodic, period 2*pi*alpha, h_phi = 2*pi*alpha/N_phi
  - Variable warp radial: smooth-tip r(rho) = r_h * sqrt(1 + (rho/r_h)^2)
  - S^2 sector with position-dependent mass (n+1)^2 / r(rho)^2

The squared Dirac on the cigar at variable warp r(rho) + conical defect alpha
(Level 1, per G4-4b convention -- spin-connection cross term deferred):

    D_cigar^2 = (D_{wedge_alpha})^2 (x) I_{S^2} + I_{wedge} (x) D_{S^2}^2 / r(rho)^2

Per (S^2 mode n, azimuthal Fourier mode k_phi):

    m_eff_sq = (2/h_phi)^2 * sin^2(pi (k+0.5)/N_phi)   with h_phi = 2*pi*alpha/N_phi
    H_{n, k_phi}(rho) = L_disk(m_eff) + diag((n+1)^2 / r(rho)^2)

Multiplicity per (n, k_phi) radial eigenvalue: 16(n+1)
    [2 (rank-2 disk spinor) * 8(n+1) (S^2 Dirac: both signs * angular deg 4(n+1))]

F6 LOAD-BEARING extensions:
    F6-A: at alpha=1, JointWarpConicalDirac = G4-4b's VariableWarpDirac.smooth_tip
          bit-exact (rel_err < 1e-10 at every t).
    F6-B: at constant warp r(rho) = r_h, JointWarpConicalDirac =
          G4-4c's DiscreteWedgeDirac x S^2 spectrum bit-exact.

S_BH extraction via replica derivative at alpha=1:
    dK_joint/dalpha|_{1} computed by central FD at eps = k_step/N_0 = 12/120 = 0.1
    Bulk subtraction: tip_term(t) = dK/dalpha - K_alpha1(t)  [bulk linear-in-alpha]
    Mellin integration: S_BH = +(1/2) integral (dt/t) f(t Lambda^2) tip_term

Continuum prediction (full cigar with horizon at conical tip):
    S_BH = A_horizon * Lambda^2 / (12 pi) = (4 pi r_h^2) Lambda^2 / (12 pi)
         = r_h^2 Lambda^2 / 3
    For r_h = 2, Lambda = 1: S_BH = 4/3 ~ 1.33
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    VariableWarpDirac,
    DiscreteWedgeDirac,
    S2DiracSpectrum,
    WarpedDiracConstant,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_5c_joint_warp_conical.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# ============================================================================
# JointWarpConicalDirac class -- driver-level (NOT in geovac/)
# ============================================================================


@dataclass
class JointWarpConicalDirac:
    """Squared Dirac on cigar with variable warp r(rho) + conical defect alpha.

    Driver-level construction combining G4-4b's VariableWarpDirac smooth-tip
    profile with G4-4c's DiscreteWedgeDirac azimuthal wedge structure, at
    Level 1 (position-dependent S^2 mass only; spin-connection cross term
    deferred per G4-4b convention).

    Parameters
    ----------
    N_rho : int
        Number of radial sites.
    a : float
        Radial lattice spacing.
    N_phi : int
        Number of azimuthal sites on the wedge.
    alpha : float
        Conical-defect parameter; apex angle = 2 pi alpha.
    r_h : float
        Horizon radius (and S^2 reference radius).
    sphere : S2DiracSpectrum
        Camporesi-Higuchi S^2 spectrum (l_max = lmax).

    Convention: proper-wedge-lattice (fix h_phi_ref = 2 pi / N_0, choose
    N_phi(alpha) = alpha * N_0). This is enforced at the caller level.
    """

    N_rho: int
    a: float
    N_phi: int
    alpha: float
    r_h: float
    sphere: S2DiracSpectrum

    def __post_init__(self) -> None:
        if self.N_rho < 1:
            raise ValueError(f"N_rho must be >= 1, got {self.N_rho}")
        if self.a <= 0:
            raise ValueError(f"a must be > 0, got {self.a}")
        if self.N_phi < 2:
            raise ValueError(f"N_phi must be >= 2, got {self.N_phi}")
        if self.alpha <= 0:
            raise ValueError(f"alpha must be > 0, got {self.alpha}")
        if self.r_h <= 0:
            raise ValueError(f"r_h must be > 0, got {self.r_h}")

    @property
    def h_phi(self) -> float:
        """h_phi = 2 pi alpha / N_phi (wedge convention)."""
        return 2 * np.pi * self.alpha / self.N_phi

    @property
    def rho_array(self) -> np.ndarray:
        """rho_k = k * a for k = 1, ..., N_rho."""
        return np.array([(k + 1) * self.a for k in range(self.N_rho)])

    @property
    def warp_profile(self) -> np.ndarray:
        """Smooth-tip warp r(rho) = r_h * sqrt(1 + (rho/r_h)^2)."""
        rho = self.rho_array
        return self.r_h * np.sqrt(1.0 + (rho / self.r_h) ** 2)

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Symmetric tridiagonal radial Laplacian (G4-3a-cleanup convention).

        u = sqrt(rho) f substitution; centrifugal -> (m_eff^2 - 1/4)/rho^2.
        """
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def H_block(self, n: int, k_phi_idx: int) -> np.ndarray:
        """Radial Hamiltonian L_wedge(k_phi) + diag((n+1)^2 / r(rho)^2).

        Combines the wedge-azimuthal m_eff (anti-periodic, period 2*pi*alpha)
        with the variable-warp position-dependent mass.
        """
        # Wedge anti-periodic Fourier mode k -> m_eff
        if k_phi_idx <= self.N_phi // 2:
            k = k_phi_idx
        else:
            k = k_phi_idx - self.N_phi
        m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
            np.pi * (k + 0.5) / self.N_phi
        ) ** 2
        m_eff = float(np.sqrt(m_eff_sq))
        L_k = self._hermitian_radial_laplacian(m_eff)
        # Variable-warp position-dependent S^2 mass diagonal
        mass_diag = (n + 1) ** 2 / self.warp_profile ** 2
        return L_k + np.diag(mass_diag)

    def heat_trace(self, t: float) -> float:
        """Tr exp(-t D_cigar^2) on the joint substrate.

        Direct sum over (n, k_phi) blocks with multiplicity 16(n+1)
        per radial eigenvalue.
        """
        K = 0.0
        for n in range(self.sphere.l_max + 1):
            mult = 16 * (n + 1)
            for k_phi_idx in range(self.N_phi):
                H = self.H_block(n, k_phi_idx)
                evals = np.linalg.eigvalsh(H)
                K += mult * float(np.sum(np.exp(-evals * t)))
        return K


# ============================================================================
# F6 extension checks (LOAD-BEARING)
# ============================================================================


def verify_F6_A_alpha1_reduces_to_variable_warp(
    N_rho: int,
    a: float,
    N_0: int,
    r_h: float,
    l_max: int,
    t_values: List[float],
    tol: float = 1e-10,
) -> dict:
    """At alpha=1, joint reduces to G4-4b VariableWarpDirac.smooth_tip bit-exact.

    Both objects share:
      - Same anti-periodic azimuthal BC, period 2*pi (alpha=1), N_phi = N_0
      - Same smooth-tip warp profile r(rho)
      - Same S^2 sphere spectrum
    Therefore heat traces must match to floating-point precision.
    """
    sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
    # Joint at alpha=1
    joint_a1 = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0, r_h=r_h, sphere=sphere
    )
    # Variable-warp smooth-tip reference
    disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    var_ref = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)

    results = {}
    all_passed = True
    for t in t_values:
        K_joint = joint_a1.heat_trace(t)
        K_ref = var_ref.heat_trace(t)
        denom = abs(K_ref) if K_ref != 0 else 1.0
        rel_err = abs(K_joint - K_ref) / denom
        passed = bool(rel_err < tol)
        all_passed = all_passed and passed
        results[str(t)] = {
            "K_joint_alpha1": float(K_joint),
            "K_var_smooth_tip": float(K_ref),
            "rel_err": float(rel_err),
            "passed": passed,
        }
    results["all_passed"] = all_passed
    results["tolerance"] = tol
    return results


def verify_F6_B_constant_warp_reduces_to_wedge_x_S2(
    N_rho: int,
    a: float,
    N_phi: int,
    alpha: float,
    r_h: float,
    l_max: int,
    t_values: List[float],
    tol: float = 1e-10,
) -> dict:
    """At constant warp r(rho) = r_h, joint reduces to wedge x S^2 spectrum.

    With r(rho) = r_h (constant warp), the position-dependent mass collapses
    to a constant (n+1)^2 / r_h^2, and the joint factorizes:
        D^2 = D_wedge^2 (x) I + I (x) (D_{S^2}/r_h)^2
        K_joint = K_wedge * K_{S^2}
    where K_wedge is DiscreteWedgeDirac.heat_trace and K_{S^2} is the
    Camporesi-Higuchi S^2 spinor heat trace at radius r_h.

    Constructs an auxiliary JointWarpConicalDirac with override:
    we instantiate via private kwarg by patching warp_profile after init.
    """
    sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
    joint = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha, r_h=r_h, sphere=sphere
    )
    # Patch warp profile to constant r_h (override smooth-tip)
    # Use a derived class with the constant-warp profile
    joint_const = _JointConstantWarp(
        N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha, r_h=r_h, sphere=sphere
    )

    wedge = DiscreteWedgeDirac(N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha)

    # K_{S^2} at radius r_h
    s2_eigs = sphere.squared_eigenvalues()

    results = {}
    all_passed = True
    for t in t_values:
        K_joint_const = joint_const.heat_trace(t)
        # Factorized: K_wedge * K_{S^2}
        K_wedge = wedge.heat_trace(t)
        K_S2 = float(np.sum(np.exp(-s2_eigs * t)))
        K_factorized = K_wedge * K_S2
        denom = abs(K_factorized) if K_factorized != 0 else 1.0
        rel_err = abs(K_joint_const - K_factorized) / denom
        passed = bool(rel_err < tol)
        all_passed = all_passed and passed
        results[str(t)] = {
            "K_joint_const_warp": float(K_joint_const),
            "K_wedge_times_K_S2": float(K_factorized),
            "rel_err": float(rel_err),
            "passed": passed,
        }
    results["all_passed"] = all_passed
    results["tolerance"] = tol
    return results


@dataclass
class _JointConstantWarp(JointWarpConicalDirac):
    """Override warp_profile to constant r_h for F6-B check."""

    @property
    def warp_profile(self) -> np.ndarray:
        return np.full(self.N_rho, self.r_h)


# ============================================================================
# Main: F11 closure driver
# ============================================================================


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-5c -- F11 closure: joint warp + conical-defect Dirac")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Substrate parameters (match G4-4 panel)
    # ------------------------------------------------------------------
    N_rho = 200
    a = 0.05
    N_0 = 120
    r_h = 2.0
    l_max = 3
    R = N_rho * a

    print(f"\n[Setup]")
    print(f"  R = {R}, N_rho = {N_rho}, a = {a}")
    print(f"  N_0 = {N_0} (reference smooth-disk azimuth count)")
    print(f"  r_h = {r_h} (horizon radius)")
    print(f"  l_max = {l_max} (S^2 truncation)")

    results["substrate"] = {
        "N_rho": N_rho, "a": a, "N_0": N_0, "r_h": r_h,
        "l_max": l_max, "R": R,
    }

    # ------------------------------------------------------------------
    # F6 extension checks: load-bearing
    # ------------------------------------------------------------------
    t_F6 = [0.05, 0.1, 0.5, 1.0]

    print("\n" + "=" * 72)
    print("[F6-A LOAD-BEARING] alpha=1 reduces to G4-4b VariableWarpDirac.smooth_tip")
    print("=" * 72)
    F6_A_results = verify_F6_A_alpha1_reduces_to_variable_warp(
        N_rho=N_rho, a=a, N_0=N_0, r_h=r_h, l_max=l_max, t_values=t_F6,
    )
    print(f"\n  {'t':>5}  {'K_joint(alpha=1)':>18}  {'K_var_smooth_tip':>18}  "
          f"{'rel_err':>12}  {'pass':>5}")
    print("  " + "-" * 70)
    for t in t_F6:
        r = F6_A_results[str(t)]
        marker = "PASS" if r["passed"] else "FAIL"
        print(f"  {t:>5.2f}  {r['K_joint_alpha1']:>18.4f}  "
              f"{r['K_var_smooth_tip']:>18.4f}  {r['rel_err']:>12.3e}  {marker:>5}")
    print(f"\n  F6-A all_passed: {F6_A_results['all_passed']}")
    results["F6_A_variable_warp_recovery"] = F6_A_results

    print("\n" + "=" * 72)
    print("[F6-B LOAD-BEARING] constant warp reduces to DiscreteWedgeDirac x S^2")
    print("=" * 72)
    # Test F6-B at multiple alpha values
    alpha_F6_B = [0.5, 1.0, 1.5]
    F6_B_all = {}
    for alpha_test in alpha_F6_B:
        N_phi_test = int(round(alpha_test * N_0))
        F6_B_results = verify_F6_B_constant_warp_reduces_to_wedge_x_S2(
            N_rho=N_rho, a=a, N_phi=N_phi_test, alpha=alpha_test,
            r_h=r_h, l_max=l_max, t_values=t_F6,
        )
        print(f"\n  Alpha = {alpha_test} (N_phi = {N_phi_test}):")
        print(f"  {'t':>5}  {'K_joint_const':>16}  {'K_wedge*K_{S^2}':>16}  "
              f"{'rel_err':>12}  {'pass':>5}")
        print("  " + "-" * 65)
        for t in t_F6:
            r = F6_B_results[str(t)]
            marker = "PASS" if r["passed"] else "FAIL"
            print(f"  {t:>5.2f}  {r['K_joint_const_warp']:>16.4f}  "
                  f"{r['K_wedge_times_K_S2']:>16.4f}  {r['rel_err']:>12.3e}  "
                  f"{marker:>5}")
        print(f"  F6-B at alpha={alpha_test} all_passed: "
              f"{F6_B_results['all_passed']}")
        F6_B_all[str(alpha_test)] = F6_B_results

    results["F6_B_constant_warp_factorization"] = F6_B_all

    F6_A_all_passed = F6_A_results["all_passed"]
    F6_B_all_passed = all(
        F6_B_all[str(a_test)]["all_passed"] for a_test in alpha_F6_B
    )
    F6_extensions_passed = F6_A_all_passed and F6_B_all_passed
    print(f"\n[F6 extension overall] passed: {F6_extensions_passed}")
    results["F6_extensions_passed"] = F6_extensions_passed

    # ------------------------------------------------------------------
    # S_BH extraction via replica derivative
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[S_BH extraction] Replica derivative at alpha=1, joint smooth-tip warp")
    print("=" * 72)

    # Replica eps: matches G4-4f best window (k_step = 12)
    k_step = 12
    N_plus = N_0 + k_step
    N_minus = N_0 - k_step
    alpha_plus = N_plus / N_0
    alpha_minus = N_minus / N_0
    eps_alpha = (alpha_plus - alpha_minus) / 2

    print(f"\n[Replica parameters]")
    print(f"  k_step = {k_step}, alpha_+ = {alpha_plus}, alpha_- = {alpha_minus}")
    print(f"  eps = (alpha_+ - alpha_-)/2 = {eps_alpha}")

    # t grid for Mellin integration
    t_grid = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    print(f"  t grid: {t_grid.tolist()}")

    sphere = S2DiracSpectrum(l_max=l_max, r_h=r_h)
    joint_plus = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_plus, alpha=alpha_plus,
        r_h=r_h, sphere=sphere,
    )
    joint_minus = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_minus, alpha=alpha_minus,
        r_h=r_h, sphere=sphere,
    )
    joint_alpha1 = JointWarpConicalDirac(
        N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0,
        r_h=r_h, sphere=sphere,
    )

    print(f"\n[Computing K_joint(alpha, t) at each t in grid...]")
    K_alpha1_arr = np.zeros(len(t_grid))
    K_plus_arr = np.zeros(len(t_grid))
    K_minus_arr = np.zeros(len(t_grid))

    for i, t in enumerate(t_grid):
        K_alpha1_arr[i] = joint_alpha1.heat_trace(t)
        K_plus_arr[i] = joint_plus.heat_trace(t)
        K_minus_arr[i] = joint_minus.heat_trace(t)
        print(f"  t = {t:>5.2f}: K(alpha=1) = {K_alpha1_arr[i]:>14.4f}, "
              f"K(+) = {K_plus_arr[i]:>14.4f}, K(-) = {K_minus_arr[i]:>14.4f}")

    # Central FD derivative
    dK_dalpha = (K_plus_arr - K_minus_arr) / (alpha_plus - alpha_minus)
    # Bulk subtraction: linear-in-alpha part is K_alpha1 (Weyl bulk scales with vol ~ alpha)
    tip_term = dK_dalpha - K_alpha1_arr

    print(f"\n  {'t':>5}  {'K(alpha=1)':>14}  {'dK/dalpha':>14}  "
          f"{'tip_term':>14}")
    print("  " + "-" * 60)
    for i, t in enumerate(t_grid):
        print(f"  {t:>5.2f}  {K_alpha1_arr[i]:>14.4f}  "
              f"{dK_dalpha[i]:>14.4f}  {tip_term[i]:>14.4f}")

    results["replica_data"] = {
        "k_step": k_step, "eps_alpha": eps_alpha,
        "alpha_plus": alpha_plus, "alpha_minus": alpha_minus,
        "t_grid": t_grid.tolist(),
        "K_alpha1": K_alpha1_arr.tolist(),
        "K_plus": K_plus_arr.tolist(),
        "K_minus": K_minus_arr.tolist(),
        "dK_dalpha": dK_dalpha.tolist(),
        "tip_term": tip_term.tolist(),
    }

    # ------------------------------------------------------------------
    # Mellin integration for S_BH at Lambda values
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("[Mellin integration] S_BH = +(1/2) integral (dt/t) f(t Lambda^2) tip(t)")
    print("=" * 72)

    Lambda_values = [0.5, 1.0, 2.0]
    print(f"\n  Lambda values: {Lambda_values}")
    print(f"  Gaussian cutoff: f(x) = exp(-x)")
    print(f"  Integration: trapezoid in log(t)")

    log_t = np.log(t_grid)

    S_BH_results = {}
    print(f"\n  {'Lambda':>8}  {'J integral':>14}  {'S_BH (extr)':>14}  "
          f"{'S_BH (cont)':>14}  {'ratio':>10}")
    print("  " + "-" * 70)

    for Lambda in Lambda_values:
        f_cutoff = np.exp(-t_grid * Lambda**2)
        # integral (dt/t) g(t) = integral d(log t) g(t)
        # So integrand in log space is just f_cutoff * tip_term
        integrand_log = f_cutoff * tip_term
        J = float(np.trapezoid(integrand_log, log_t))
        S_BH_extracted = +J / 2

        # Continuum prediction: S_BH = r_h^2 Lambda^2 / 3
        S_BH_continuum = r_h**2 * Lambda**2 / 3.0
        ratio = (
            S_BH_extracted / S_BH_continuum
            if abs(S_BH_continuum) > 1e-15 else float("inf")
        )

        S_BH_results[str(Lambda)] = {
            "Lambda": Lambda,
            "J_integral": J,
            "S_BH_extracted": S_BH_extracted,
            "S_BH_continuum_prediction": S_BH_continuum,
            "ratio_disc_to_cont": ratio,
        }
        print(f"  {Lambda:>8.2f}  {J:>14.6e}  {S_BH_extracted:>+14.6e}  "
              f"{S_BH_continuum:>+14.6e}  {ratio:>10.4f}")

    results["S_BH_extraction"] = S_BH_results

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    # F6 extension verified
    F6_verified = F6_extensions_passed

    # S_BH magnitude check: within factor 2 of continuum
    S_within_factor_2 = True
    for L in Lambda_values:
        ratio = S_BH_results[str(L)]["ratio_disc_to_cont"]
        if not (0.5 <= ratio <= 2.0):
            S_within_factor_2 = False
            break

    # S_BH sign check: positive
    S_positive = all(
        S_BH_results[str(L)]["S_BH_extracted"] > 0 for L in Lambda_values
    )

    print(f"\n[Verdict]")
    print(f"  F6-A alpha=1 -> VariableWarpDirac:           {F6_A_all_passed}")
    print(f"  F6-B const warp -> wedge x S^2:              {F6_B_all_passed}")
    print(f"  S_BH positive (correct entropy sign):        {S_positive}")
    print(f"  S_BH within factor 2 of continuum:           {S_within_factor_2}")

    if F6_verified and S_positive and S_within_factor_2:
        verdict = "POSITIVE-G4-5c-F11-CLOSED"
        msg = (
            "Joint variable-warp + conical-defect Dirac F11 closed: F6 extensions "
            "verified bit-exact (alpha=1 -> VariableWarpDirac; const warp -> "
            "wedge x S^2); S_BH extracted at positive value within factor 2 of "
            "continuum r_h^2 Lambda^2/3 prediction."
        )
    elif F6_verified and S_positive:
        verdict = "PARTIAL-G4-5c-F11-MAGNITUDE-OFF"
        msg = (
            "F6 extensions verified bit-exact and S_BH sign positive, but "
            "magnitude off continuum prediction by > factor 2; substrate UV/IR "
            "or sub-leading bulk subtraction may need refinement."
        )
    elif F6_verified:
        verdict = "PARTIAL-G4-5c-F11-SIGN-ISSUE"
        msg = (
            "F6 extensions verified but S_BH sign or magnitude unexpected; "
            "diagnose bulk subtraction or Mellin convention."
        )
    else:
        verdict = "NEGATIVE-G4-5c-F11"
        msg = (
            "F6 extension failed: joint construction does not reduce to "
            "G4-4b/4c references; construction bug."
        )

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["F6_verified"] = F6_verified
    results["S_BH_positive"] = S_positive
    results["S_BH_within_factor_2"] = S_within_factor_2

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
