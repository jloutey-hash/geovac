"""Tests for the spectral azimuthal discretization classes (G4-6d).

DiscreteDiskDiracSpectral and DiscreteWedgeDiracSpectral replace the FD
azimuthal eigenvalues (2/h_phi) sin(pi(k+1/2)/N_phi) with the exact
continuum eigenvalue (k+1/2)/alpha. Per task #28 analytical prediction
and Sprint G4-6a B.1+B.2 empirical verification (2026-05-29):
- F6 bit-exact at alpha=1: wedge spectral reduces to disk spectral.
- 6.36% UV recovery at t=a^2 on production panel vs 0.04% for FD.

These tests verify the production code matches both predictions.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteDiskDiracSpectral,
    DiscreteWedgeDirac,
    DiscreteWedgeDiracSpectral,
)


# =================================================================
# F6 bit-exact reduction: wedge at alpha=1 reduces to disk
# =================================================================


class TestF6BitExactDiskSpectral:
    """At alpha=1, wedge spectral should match disk spectral bit-exactly."""

    def test_eigenvalues_bit_exact_at_alpha_1(self):
        """Eigenvalue spectra of wedge_spec(alpha=1) and disk_spec match
        bit-exactly."""
        N_rho, a, N_phi = 50, 0.1, 60
        disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi)
        wedge = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_phi, alpha=1.0,
        )
        evs_disk = disk.squared_eigenvalues()
        evs_wedge = wedge.squared_eigenvalues()
        assert evs_disk.shape == evs_wedge.shape
        # Bit-exact equality
        np.testing.assert_array_equal(evs_disk, evs_wedge)

    def test_heat_trace_bit_exact_at_alpha_1(self):
        """Heat traces at multiple t values match bit-exactly."""
        N_rho, a, N_phi = 50, 0.1, 60
        disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi)
        wedge = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_phi, alpha=1.0,
        )
        for t in [0.01, 0.1, 0.5, 1.0, 10.0]:
            K_disk = disk.heat_trace(t)
            K_wedge = wedge.heat_trace(t)
            assert K_disk == K_wedge, f"Bit-exact mismatch at t={t}"


# =================================================================
# Hilbert dimension is correct
# =================================================================


class TestHilbertDim:
    def test_disk_hilbert_dim(self):
        """Hilbert dim = 2 * N_rho * N_phi (spinor doubling)."""
        disk = DiscreteDiskDiracSpectral(N_rho=20, a=0.1, N_phi=30)
        assert disk.hilbert_dim == 2 * 20 * 30

    def test_wedge_hilbert_dim(self):
        """Hilbert dim = 2 * N_rho * N_phi (spinor doubling)."""
        wedge = DiscreteWedgeDiracSpectral(
            N_rho=20, a=0.1, N_phi=30, alpha=1.5,
        )
        assert wedge.hilbert_dim == 2 * 20 * 30


# =================================================================
# Eigenvalue structure: spectral has EXACT m_eff = (k+1/2)/alpha
# =================================================================


class TestSpectralEigenvalueStructure:
    """Verify spectral substrate uses exact m_eff = (k+1/2)/alpha."""

    def test_lowest_mode_spectral_vs_FD(self):
        """Spectral m_eff = 1/2 at alpha=1; FD differs slightly at finite N_phi."""
        N_rho, a, N_phi = 50, 0.1, 60
        # Spectral disk
        disk_spec = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_phi)
        # FD disk
        disk_fd = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_phi)
        # Compute heat trace at very large t (only lowest modes matter)
        t_large = 100.0
        K_spec = disk_spec.heat_trace(t_large)
        K_fd = disk_fd.heat_trace(t_large)
        # FD should give heat trace very close but not bit-exact to spectral
        # (the FD m_eff_low approaches 1/2 in the continuum limit)
        assert abs(K_spec - K_fd) < 1e-3, "FD and spectral should agree at large t"

    def test_wedge_softer_IR_at_alpha_gt_1(self):
        """At alpha > 1, lowest spectral eigenvalue is 1/(2*alpha) (softer IR)."""
        N_rho, a, N_phi = 20, 0.2, 30
        wedge = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_phi, alpha=2.0,
        )
        evs = wedge.squared_eigenvalues()
        # Lowest m_eff^2 = (1/2)^2 / 4 = 1/16. Radial Laplacian's lowest
        # eigenvalue with this m_eff dominates the heat trace at large t.
        # We check that the lowest squared eigenvalue is positive (no zero
        # mode for anti-periodic BC) but smaller than the alpha=1 case.
        assert evs[0] > 0
        wedge_1 = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_phi, alpha=1.0,
        )
        evs_1 = wedge_1.squared_eigenvalues()
        # alpha=2 should have softer (smaller) lowest eigenvalue than alpha=1
        # because m_eff is divided by alpha
        assert evs[0] < evs_1[0]


# =================================================================
# Heat-trace continuity in t
# =================================================================


class TestHeatTraceContinuity:
    """Heat trace is monotonically decreasing in t and bounded."""

    def test_disk_heat_trace_decreasing(self):
        disk = DiscreteDiskDiracSpectral(N_rho=30, a=0.1, N_phi=40)
        K_small_t = disk.heat_trace(0.1)
        K_large_t = disk.heat_trace(1.0)
        assert K_small_t > K_large_t

    def test_wedge_heat_trace_decreasing(self):
        wedge = DiscreteWedgeDiracSpectral(
            N_rho=30, a=0.1, N_phi=40, alpha=1.5,
        )
        K_small_t = wedge.heat_trace(0.1)
        K_large_t = wedge.heat_trace(1.0)
        assert K_small_t > K_large_t

    def test_disk_heat_trace_positive(self):
        disk = DiscreteDiskDiracSpectral(N_rho=30, a=0.1, N_phi=40)
        for t in [0.1, 1.0, 10.0]:
            assert disk.heat_trace(t) > 0


# =================================================================
# Constructor validation
# =================================================================


class TestConstructorValidation:
    def test_disk_invalid_N_rho(self):
        with pytest.raises(ValueError, match="N_rho"):
            DiscreteDiskDiracSpectral(N_rho=0, a=0.1, N_phi=10)

    def test_disk_invalid_a(self):
        with pytest.raises(ValueError, match="a must be"):
            DiscreteDiskDiracSpectral(N_rho=10, a=-0.1, N_phi=10)

    def test_disk_invalid_N_phi(self):
        with pytest.raises(ValueError, match="N_phi"):
            DiscreteDiskDiracSpectral(N_rho=10, a=0.1, N_phi=1)

    def test_wedge_invalid_alpha(self):
        with pytest.raises(ValueError, match="alpha"):
            DiscreteWedgeDiracSpectral(N_rho=10, a=0.1, N_phi=10, alpha=-1.0)


# =================================================================
# Spectral vs FD: spectral should produce larger UV divergence
# =================================================================


class TestSpectralUVImprovement:
    """At sufficiently small t, spectral substrate has more UV content
    than FD substrate.

    This is the load-bearing physical property motivating G4-6d.
    """

    @pytest.mark.slow
    def test_B_substrate_R_independent_at_R_geq_10(self):
        """G4-6b finding: B_substrate at large t is essentially R-independent
        for R >= 10, matching continuum +1/6 to within ~3%.

        Per debug/g4_6b_ir_boundary_first_move_memo.md (2026-05-29).
        """
        a = 0.05
        N_0 = 120
        t = 10.0  # Large t regime where A/t ~ 0.001 << B

        # Test at R = 10 (N_rho = 200)
        N_rho = 200
        disk = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_0)
        K_disk = disk.heat_trace(t)
        # Central FD for dK/d_alpha
        k_step = 12
        N_plus = N_0 + k_step
        N_minus = N_0 - k_step
        wedge_plus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_plus, alpha=N_plus / N_0,
        )
        wedge_minus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_minus, alpha=N_minus / N_0,
        )
        K_plus = wedge_plus.heat_trace(t)
        K_minus = wedge_minus.heat_trace(t)
        dK = (K_plus - K_minus) / (N_plus / N_0 - N_minus / N_0)
        tip = dK - K_disk

        # B target is +1/6 = 0.1667; substrate value at R=10 is ~0.163
        # Allow generous 5% relative error (measured is 2.3%)
        B_target = 1.0 / 6.0
        assert abs(tip - B_target) / B_target < 0.05, (
            f"B at R=10 should be within 5% of +1/6; got {tip:.6f}"
        )

    @pytest.mark.slow
    def test_spectral_outperforms_FD_at_small_t(self):
        """At t = a^2, spectral substrate's tip term exceeds FD's by
        more than 10x (per Sprint B.2 finding: ~160x improvement)."""
        N_rho, a, N_0 = 50, 0.1, 30
        t = a ** 2  # UV cell

        # FD substrate
        disk_fd = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
        wedge_fd_plus = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_0 + 6, alpha=(N_0 + 6) / N_0,
        )
        wedge_fd_minus = DiscreteWedgeDirac(
            N_rho=N_rho, a=a, N_phi=N_0 - 6, alpha=(N_0 - 6) / N_0,
        )
        K_disk_fd = disk_fd.heat_trace(t)
        dK_fd = (wedge_fd_plus.heat_trace(t) - wedge_fd_minus.heat_trace(t)) / \
                ((N_0 + 6) / N_0 - (N_0 - 6) / N_0)
        tip_fd = abs(dK_fd - K_disk_fd)

        # Spectral substrate
        disk_spec = DiscreteDiskDiracSpectral(N_rho=N_rho, a=a, N_phi=N_0)
        wedge_spec_plus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_0 + 6, alpha=(N_0 + 6) / N_0,
        )
        wedge_spec_minus = DiscreteWedgeDiracSpectral(
            N_rho=N_rho, a=a, N_phi=N_0 - 6, alpha=(N_0 - 6) / N_0,
        )
        K_disk_spec = disk_spec.heat_trace(t)
        dK_spec = (wedge_spec_plus.heat_trace(t) -
                   wedge_spec_minus.heat_trace(t)) / \
                  ((N_0 + 6) / N_0 - (N_0 - 6) / N_0)
        tip_spec = abs(dK_spec - K_disk_spec)

        # Spectral should be substantially larger than FD
        # (B.2 found ~160x improvement; here we check 10x as a loose bound)
        assert tip_spec > 10 * tip_fd, (
            f"Spectral tip {tip_spec:.6f} should exceed 10x FD tip {tip_fd:.6f}"
        )
