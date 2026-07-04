"""Pins for the three Paper 51 §G4-4 headline numbers (§C8 backing).

Paper 51 (papers/group5_qed_gauge/paper_51_gravity_arc.tex) §G4-4 states
three MEASURED sprint-scale headline numbers that previously had no live
test (group5 1st-cert finding E11):

  (1) Cone-Dirac saturation (§G4-4b, "Cone-Dirac saturation" paragraph):
        K^Dirac_var(t; r_h -> 0) = K^Dirac_cone(t)
      "verified bit-exact to 6+ digits at r_h = 10^-3", where K_cone is
      the singular-tip warp r(rho) = rho, NOT the naive ceiling
      4(l_max+1)(l_max+2) * K_disk.

  (2) SC slope (§G4-4c + BC-sector table): the spinor conical-defect tip
      coefficient Delta_K^{Dirac,tip}(alpha) = -(1/12)(1/alpha - alpha),
      mean recovery 99.4% across alpha in {2/5, 1/2, 3/5, 2/3} at t = 1.0
      (the G4-4e side-by-side panel), and best single extraction
      recovery = 1.000021 (rel err 2.1e-5) at alpha = 2/5, t = 2.0
      (G4-4c week 3).

  (3) Seeley-DeWitt a_0 (§G4-4d): at the sweet spot t = 0.1, N_phi = 192,
        a_0^Dirac = 1.992   (continuum 2.0, recovery 99.6%).

All three are recomputed here AT THE PAPER'S OWN SUBSTRATE from the
tested production module geovac/gravity/warped_dirac.py (no archived
debug/ driver needed):

  (1) substrate N_rho=20, a=0.3, N_phi=12, l_max=3 (the G4-4b-c sprint
      substrate; same substrate as TestF5ConeDiracSaturation in
      tests/test_paper51_L6_full.py, which checks only ratio bounds +
      plateau stability -- this file adds the INDEPENDENT cone-Dirac
      reference the paper's claim actually compares against);
  (2) substrate R=10, a=0.05, N_rho=200, N_0=120 (production);
  (3) substrate N_rho=200, a=0.05, N_phi=192 (production UV-refined).

Validated-before-pinned 2026-07-03 (values measured, then asserted):
  (1) rel diff |K_var(r_h=1e-3) - K_cone| / K_cone = 5.4e-8 .. 6.4e-8
      across t in {0.1, 0.5, 1.0} (7+ digits; paper claims 6+), with
      quadratic r_h^2 convergence (6e-6 at r_h=1e-2, 6e-10 at r_h=1e-4);
      plateau ratio K_cone/K_disk(t=0.5) = 50.847785.
  (2) recoveries 0.999902 / 0.998057 / 0.992821 / 0.986470,
      mean = 0.994312; best cell 1.0000212.
  (3) a_0 = 1.991913.

Total runtime ~ 3 s (no slow marks needed).

Per CLAUDE.md §13.4a these are numerical cross-checks of MEASURED
(sprint-scale) claims: the asserts pin the published values with bands
that the measured values pass with margin but that would catch any
structural regression (wrong BC, wrong multiplicity, wrong warp).
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteWedgeDirac,
    S2DiracSpectrum,
    VariableWarpDirac,
)


# =====================================================================
# (1) Cone-Dirac saturation: K_var(r_h -> 0) = K_cone
# =====================================================================
#
# The independent reference is built HERE, from scratch (plain numpy
# tridiagonals), not via VariableWarpDirac: singular-tip warp
# r(rho) = rho means the S^2 mass term is (n+1)^2 / rho_k^2 exactly.
# Per (n, k_phi) block the radial Hamiltonian is
#
#   H = tridiag(-1/a^2,
#               2/a^2 + (m_eff^2 - 1/4)/rho^2 + (n+1)^2/rho^2,
#               -1/a^2)
#
# (Hermitian polar Laplacian, u = sqrt(rho) f convention, anti-periodic
# FD azimuthal m_eff) with multiplicity 16(n+1) = 2 (disk spinor rank)
# x 8(n+1) (S^2 Dirac, both signs).

SUBSTRATE = dict(N_rho=20, a=0.3, N_phi=12)
L_MAX = 3
NAIVE_CEILING = 4 * (L_MAX + 1) * (L_MAX + 2)  # = 80


def _m_eff_antiperiodic(k_idx: int, N_phi: int) -> float:
    """Anti-periodic FD azimuthal eigenvalue (disk-Dirac convention)."""
    if k_idx <= N_phi // 2:
        k = k_idx
    else:
        k = k_idx - N_phi
    h_phi = 2 * np.pi / N_phi
    return abs(2.0 / h_phi * np.sin(np.pi * (k + 0.5) / N_phi))


def _K_cone_independent(t: float) -> float:
    """Cone-Dirac heat trace (r(rho) = rho), built independently."""
    N_rho = SUBSTRATE["N_rho"]
    a = SUBSTRATE["a"]
    N_phi = SUBSTRATE["N_phi"]
    rho = np.arange(1, N_rho + 1) * a
    K = 0.0
    for n in range(L_MAX + 1):
        mult = 16 * (n + 1)
        for k_idx in range(N_phi):
            m = _m_eff_antiperiodic(k_idx, N_phi)
            diag = (2.0 / a**2 + (m**2 - 0.25) / rho**2
                    + (n + 1) ** 2 / rho**2)
            H = np.diag(diag)
            off = -np.ones(N_rho - 1) / a**2
            H += np.diag(off, 1) + np.diag(off, -1)
            evals = np.linalg.eigvalsh(H)
            K += mult * float(np.sum(np.exp(-evals * t)))
    return K


def _K_var(r_h: float, t: float, disk: DiscreteDiskDirac) -> float:
    sphere = S2DiracSpectrum(l_max=L_MAX, r_h=r_h)
    var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)
    return var.heat_trace(t)


class TestConeDiracSaturationReference:
    """§G4-4b headline: K_var(r_h -> 0) = K_cone to 6+ digits at r_h=1e-3.

    Substrate: N_rho=20, a=0.3, N_phi=12, l_max=3 (the G4-4b-c sprint
    substrate on which the paper's claim was recorded). At this substrate
    the measured agreement at r_h = 10^-3 is 5.4e-8..6.4e-8 relative
    (7+ digits), comfortably inside the paper's "6+ digits".

    Honest-scope note (2026-07-03 delta-2 review): at this coarse
    substrate the smallest lattice site rho_1 = a = 0.3 >> r_h = 1e-3,
    so sqrt(rho^2 + r_h^2) ~ rho to ~1.7e-6 pointwise and much of the
    single-r_h agreement is forced by scale separation; the load-bearing
    falsifier is test_kvar_converges_to_cone_quadratically_in_rh below,
    whose O(r_h^2) scaling law a wrong operator would break.
    """

    def test_kvar_matches_independent_cone_to_6_digits_at_rh_1e3(self):
        """|K_var(1e-3) - K_cone| / K_cone < 1e-6 at t in {0.1, 0.5, 1.0}.

        1e-6 relative IS the paper's "6+ digits" gate; the measured
        values (5.4e-8..6.4e-8) pass it with ~15x margin.
        """
        disk = DiscreteDiskDirac(**SUBSTRATE)
        for t in (0.1, 0.5, 1.0):
            K_cone = _K_cone_independent(t)
            K_var = _K_var(1e-3, t, disk)
            rel = abs(K_var - K_cone) / K_cone
            assert rel < 1e-6, (
                f"t={t}: K_var(r_h=1e-3)={K_var!r} vs independent "
                f"K_cone={K_cone!r}, rel diff {rel:.3e} >= 1e-6"
            )
            # Tighter regression band around the measured 5.4e-8..6.4e-8:
            assert rel < 2e-7, (
                f"t={t}: rel diff {rel:.3e} regressed above the measured "
                f"~6e-8 band (validated 2026-07-03)"
            )

    def test_kvar_converges_to_cone_quadratically_in_rh(self):
        """rel diff drops ~100x per decade of r_h (smooth-tip r(rho) =
        sqrt(rho^2 + r_h^2) deviates from the cone at O(r_h^2)).

        Measured at t=0.5: 5.43e-6 (r_h=1e-2) -> 5.43e-8 (1e-3)
        -> 5.43e-10 (1e-4)."""
        disk = DiscreteDiskDirac(**SUBSTRATE)
        t = 0.5
        K_cone = _K_cone_independent(t)
        rels = [abs(_K_var(r_h, t, disk) - K_cone) / K_cone
                for r_h in (1e-2, 1e-3, 1e-4)]
        for i in range(1, len(rels)):
            factor = rels[i - 1] / rels[i]
            assert 50.0 < factor < 200.0, (
                f"Expected ~100x per r_h decade (quadratic); step {i} "
                f"gave {factor:.1f}x (rels={rels})"
            )

    def test_cone_plateau_ratio_value_and_below_naive_ceiling(self):
        """The saturation plateau K_cone/K_disk at t=0.5 equals 50.848
        (validated 50.847785), NOT the naive ceiling
        4(l_max+1)(l_max+2) = 80. This pins the value the prior F5 test
        (tests/test_paper51_L6_full.py) only bounds."""
        disk = DiscreteDiskDirac(**SUBSTRATE)
        t = 0.5
        ratio = _K_cone_independent(t) / disk.heat_trace(t)
        assert abs(ratio - 50.848) < 0.01, (
            f"K_cone/K_disk(t=0.5) = {ratio:.6f}, expected 50.848 +/- 0.01"
        )
        assert ratio < NAIVE_CEILING


# =====================================================================
# (2) SC slope: spinor -1/12 recovery 99.4% (mean) / 1.000021 (best)
# =====================================================================


class TestSpinorSCSlopeRecovery:
    """§G4-4c/§G4-4e headline: spinor conical-defect tip coefficient
    -(1/12)(1/alpha - alpha) recovered at 99.4% (mean over the G4-4e
    panel) on the production substrate R=10, a=0.05, N_rho=200, N_0=120.

    Extraction (the g4_4c/g4_4e protocol): at fixed t,
        slope = [K_wedge(alpha) - alpha * K_disk] / (1/alpha - alpha),
        recovery = slope / (-1/12).
    Proper wedge lattice: h_phi fixed by N_0 = 120, N_phi = alpha * N_0.
    """

    R = 10.0
    A_LAT = 0.05
    N_RHO = 200
    N_0 = 120
    # (label, alpha, N_phi = alpha * N_0)
    PANEL = [("2/5", 2.0 / 5.0, 48), ("1/2", 0.5, 60),
             ("3/5", 3.0 / 5.0, 72), ("2/3", 2.0 / 3.0, 80)]

    def _recovery(self, alpha: float, N_phi: int, t: float,
                  K_disk: float) -> float:
        wedge = DiscreteWedgeDirac(
            N_rho=self.N_RHO, a=self.A_LAT, N_phi=N_phi, alpha=alpha,
        )
        slope = (wedge.heat_trace(t) - alpha * K_disk) / (1 / alpha - alpha)
        return slope / (-1.0 / 12.0)

    def test_mean_recovery_99_4_percent_at_t1(self):
        """Mean recovery over the G4-4e panel at t = 1.0 equals 0.994
        (validated 0.994312; band +/- 0.004). Each individual alpha
        recovers in [0.985, 1.001] (validated 0.98647..0.99990)."""
        disk = DiscreteDiskDirac(
            N_rho=self.N_RHO, a=self.A_LAT, N_phi=self.N_0,
        )
        K_disk = disk.heat_trace(1.0)
        recs = []
        for label, alpha, N_phi in self.PANEL:
            rec = self._recovery(alpha, N_phi, 1.0, K_disk)
            assert 0.985 < rec < 1.001, (
                f"alpha={label}: recovery {rec:.6f} outside [0.985, 1.001]"
            )
            recs.append(rec)
        mean = float(np.mean(recs))
        assert abs(mean - 0.994) < 0.004, (
            f"Mean spinor SC recovery {mean:.6f}, paper headline 0.994"
        )

    def test_best_extraction_alpha_2_5_t2(self):
        """G4-4c week-3 best cell: alpha = 2/5, t = 2.0 gives recovery
        1.000021 (rel err 2.1e-5). Validated 1.0000212; band 1e-4."""
        disk = DiscreteDiskDirac(
            N_rho=self.N_RHO, a=self.A_LAT, N_phi=self.N_0,
        )
        K_disk = disk.heat_trace(2.0)
        rec = self._recovery(2.0 / 5.0, 48, 2.0, K_disk)
        assert abs(rec - 1.0) < 1e-4, (
            f"Best-cell recovery {rec:.7f}, paper claims rel err 2.1e-5"
        )


# =====================================================================
# (3) Seeley-DeWitt a_0 sweet spot: a_0^Dirac = 1.992
# =====================================================================


class TestA0DiracSweetSpot:
    """§G4-4d headline: a_0^Dirac = 1.992 (continuum 2, recovery 99.6%)
    at the sweet spot t = 0.1 on the UV-refined production substrate
    N_rho=200, a=0.05, N_phi=192 (R = 10, A = pi R^2).

    Extraction: continuum Weyl leading order K_Dirac(t) ~ a_0 A/(4 pi t)
    with a_0 = 2 (rank-2 spinor bundle), so
        a_0(measured) = K(0.1) * 4 pi t / A.
    Validated 2026-07-03: K(0.1) = 497.9782, a_0 = 1.991913.
    """

    def test_a0_dirac_equals_1_992_at_sweet_spot(self):
        disk = DiscreteDiskDirac(N_rho=200, a=0.05, N_phi=192)
        R = disk.R
        assert R == pytest.approx(10.0)
        A = np.pi * R**2
        t = 0.1
        a0 = disk.heat_trace(t) * 4 * np.pi * t / A
        # Pin the published 3-decimal value with a +/- 0.002 band
        # (measured 1.991913).
        assert abs(a0 - 1.992) < 2e-3, (
            f"a_0^Dirac(t=0.1, N_phi=192) = {a0:.6f}, paper pins 1.992"
        )
        # And the recovery statement (99.6% of continuum 2.0).
        assert abs(a0 / 2.0 - 0.996) < 1e-3
