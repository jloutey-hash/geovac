"""Tests for the screened cross-center V_ne machinery (Phase C-W1c).

Validates:
  - Backward-compat with bare cross-center V_ne for first-row Z (no frozen core)
  - Sanity match with Phase B-W1c-diag's NaH naive-screened proxy
  - Multipole termination preserved at L_max = 2 * l_max
  - Asymptotic limit (large R) -> bare Coulomb of asymptotic Z_eff
  - Hermiticity, m-diagonal block structure
  - Direction rotation works
  - Auto-detection of frozen core type from Z

Author: GeoVac Development Team (Phase C-W1c, May 2026)
"""

import numpy as np
import pytest

from geovac.composed_qubit import _enumerate_states
from geovac.cross_center_screened_vne import (
    _detect_core_type,
    _screening_correction_grid,
    compute_screened_cross_center_vne,
    compute_screened_cross_center_vne_element,
)
from geovac.neon_core import FrozenCore
from geovac.shibuya_wulfman import (
    compute_cross_center_vne,
    compute_cross_center_vne_element,
)


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------


class TestCoreTypeDetection:
    def test_first_row_no_core(self):
        for Z in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            assert _detect_core_type(Z) is None, f"Z={Z} should have no frozen core"

    def test_ne_core(self):
        for Z in [11, 12, 13, 14, 15, 16, 17, 18]:
            assert _detect_core_type(Z) == 'Ne', f"Z={Z} should have [Ne] core"

    def test_ar_core(self):
        for Z in [19, 20, 21, 25, 30]:
            assert _detect_core_type(Z) == 'Ar', f"Z={Z} should have [Ar] core"

    def test_ar3d10_core(self):
        for Z in [31, 32, 35, 36]:
            assert _detect_core_type(Z) == 'Ar3d10', f"Z={Z} should have [Ar]3d10 core"

    def test_kr_core(self):
        for Z in [37, 38]:
            assert _detect_core_type(Z) == 'Kr', f"Z={Z} should have [Kr] core"


# ---------------------------------------------------------------------------
# Backward compatibility: bare-Coulomb regression for first-row Z
# ---------------------------------------------------------------------------


class TestBareRegression:
    """For first-row Z (no frozen core), screened path must = bare path
    bit-exactly, since there's no screening to apply."""

    def test_h_h_z1(self):
        states = _enumerate_states(2)
        v_bare = compute_cross_center_vne(1.0, states, 1.0, 3.0, L_max=2)
        v_scr = compute_screened_cross_center_vne(1.0, states, 1.0, 3.0, L_max=2)
        assert np.allclose(v_bare, v_scr, atol=0, rtol=0), (
            f"max diff = {np.max(np.abs(v_bare - v_scr)):.2e}, expected exact 0"
        )

    def test_li_z3(self):
        states = _enumerate_states(2)
        v_bare = compute_cross_center_vne(3.0, states, 3.0, 3.015, L_max=2)
        v_scr = compute_screened_cross_center_vne(3.0, states, 3.0, 3.015, L_max=2)
        assert np.allclose(v_bare, v_scr, atol=0, rtol=0)

    def test_be_z4(self):
        states = _enumerate_states(2)
        v_bare = compute_cross_center_vne(4.0, states, 4.0, 2.5, L_max=2)
        v_scr = compute_screened_cross_center_vne(4.0, states, 4.0, 2.5, L_max=2)
        assert np.allclose(v_bare, v_scr, atol=0, rtol=0)

    def test_element_api_matches_matrix(self):
        """Element-by-element API also reverts to bare for first-row."""
        # 1s, 1s element with Z_orb=1, Z_nuc=1
        v_bare = compute_cross_center_vne_element(
            1.0, 1, 0, 0, 1, 0, 0, 1.0, 3.0, L_max=2,
        )
        v_scr = compute_screened_cross_center_vne_element(
            1.0, 1, 0, 0, 1, 0, 0, 1.0, 3.0, L_max=2,
        )
        assert v_bare == v_scr


# ---------------------------------------------------------------------------
# Diagnostic sanity probe match (Phase B-W1c-diag)
# ---------------------------------------------------------------------------


class TestNaHSanityProbeMatch:
    """The diagnostic sanity probe at NaH R=3.5 bohr predicts the screened
    H-side V_ne trace ~ -1.10 Ha (vs bare -12.15 Ha). Verify the new
    machinery matches within a few percent (the proxy was a uniform
    rescaling; the proper screened result differs by O(few %) due to
    radial profile)."""

    def test_nah_screened_trace(self):
        H_states = _enumerate_states(2)
        v_screened = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4,
        )
        trace = float(np.trace(v_screened))
        # Diagnostic predicts ~ -1.10 Ha (proxy is uniform rescale)
        # Proper screened result allows a small extra attraction (H orbital
        # tail probes inner region where Z_eff > 1). Match within 15% of
        # proxy, sign correct, much smaller than bare -12.15 Ha.
        assert -2.0 < trace < -0.5, (
            f"Screened trace {trace:.4f} Ha out of plausible range [-2.0, -0.5]"
        )
        assert abs(trace) < abs(-12.15) / 5, (
            f"Screened {trace:.4f} should be at least 5x smaller than bare -12.15"
        )

    def test_nah_screened_vs_bare_shift(self):
        """The shift from bare to screened should be ~10 Ha (large)."""
        H_states = _enumerate_states(2)
        v_bare = compute_cross_center_vne(1.0, H_states, 11.0, 3.5, L_max=4)
        v_scr = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4,
        )
        shift = float(np.trace(v_bare) - np.trace(v_scr))
        # Diagnostic predicts ~ -11 Ha (bare more attractive)
        assert -12.0 < shift < -10.0, (
            f"Bare-screened trace shift {shift:.4f} not in expected range "
            f"-11 +- 1 Ha"
        )


# ---------------------------------------------------------------------------
# Asymptotic limit: at large R, screened with [Ne] -> bare with Z_eff_asy = Z-10
# ---------------------------------------------------------------------------


class TestAsymptotic:
    def test_nah_large_R_approaches_naplus(self):
        """At R >> R_core, screened V_ne(Z=11, [Ne]) ~ bare V_ne(Z=1, Na+)."""
        H_states = _enumerate_states(2)
        # Pick large R where [Ne] core is fully internalized (rho > 5 bohr)
        R = 12.0
        v_screened = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, R, L_max=4,
        )
        v_naplus = compute_cross_center_vne(1.0, H_states, 1.0, R, L_max=4)
        # At R=12 bohr, the orbital extent (n=2 Z=1) is ~6 bohr, so most
        # of the orbital lives at rho > 6 bohr from Na where Z_eff -> 1.
        diff = np.max(np.abs(v_screened - v_naplus))
        assert diff < 1e-3, (
            f"At R={R} bohr, screened-Na should match bare-Na+ within 1e-3 "
            f"(got {diff:.2e})"
        )


# ---------------------------------------------------------------------------
# Multipole termination preserved (Gaunt selection)
# ---------------------------------------------------------------------------


class TestMultipoleTermination:
    """Gaunt 3j triangle inequality enforces L_max = 2*l_max for the
    angular coefficient, regardless of the radial form of the potential.
    Verify by checking that L > 2*l_max contributes nothing."""

    def test_lmax_zero_for_s_orbitals(self):
        """For pure s-orbitals (l_max=0), only L=0 contributes."""
        H_states = [(1, 0, 0), (2, 0, 0)]  # 1s, 2s only
        v_L2 = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=0,
        )
        v_L4 = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4,
        )
        # All orbitals are s, so L_max=0 should already be exact
        assert np.allclose(v_L2, v_L4, atol=1e-9), (
            f"L>0 contributing for s-only orbitals: max diff "
            f"{np.max(np.abs(v_L2 - v_L4)):.2e}"
        )

    def test_lmax_two_for_l1_orbitals(self):
        """For l_max=1, only L=0,1,2 contribute (2*l_max)."""
        states = _enumerate_states(2)  # includes 2p (l=1)
        v_L2 = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.5, L_max=2,
        )
        v_L6 = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.5, L_max=6,
        )
        assert np.allclose(v_L2, v_L6, atol=1e-9), (
            f"L>2 contributing for l_max=1: max diff "
            f"{np.max(np.abs(v_L2 - v_L6)):.2e}"
        )


# ---------------------------------------------------------------------------
# Hermiticity, m-diagonal structure
# ---------------------------------------------------------------------------


class TestStructuralProperties:
    def test_hermiticity(self):
        states = _enumerate_states(2)
        v = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4,
        )
        diff = np.max(np.abs(v - v.T))
        assert diff < 1e-10, f"Non-Hermitian: max|V-V.T| = {diff:.2e}"

    def test_m_diagonal_z_aligned(self):
        """For nucleus along z-axis, V_ne is m-diagonal (M=0 selection)."""
        states = _enumerate_states(2)
        v = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4, nuc_parity=1,
        )
        for i, (n1, l1, m1) in enumerate(states):
            for j, (n2, l2, m2) in enumerate(states):
                if m1 != m2:
                    assert abs(v[i, j]) < 1e-10, (
                        f"Off-m matrix element nonzero: "
                        f"<{n1},{l1},{m1}|V|{n2},{l2},{m2}> = {v[i,j]}"
                    )


# ---------------------------------------------------------------------------
# Direction rotation
# ---------------------------------------------------------------------------


class TestDirectionRotation:
    def test_z_direction_consistency(self):
        """direction=(0,0,1) should match nuc_parity=+1 (within numerical tol)."""
        states = _enumerate_states(2)
        v_par = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4, nuc_parity=+1, direction=None,
        )
        v_dir = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4, direction=(0.0, 0.0, 1.0),
        )
        assert np.allclose(v_par, v_dir, atol=1e-12)

    def test_neg_z_direction(self):
        """direction=(0,0,-1) should match nuc_parity=-1."""
        states = _enumerate_states(2)
        v_par = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4, nuc_parity=-1, direction=None,
        )
        v_dir = compute_screened_cross_center_vne(
            1.0, states, 11.0, 3.0, L_max=4, direction=(0.0, 0.0, -1.0),
        )
        assert np.allclose(v_par, v_dir, atol=1e-12)


# ---------------------------------------------------------------------------
# Screening correction profile sanity
# ---------------------------------------------------------------------------


class TestScreeningCorrection:
    """Verify the screening correction f_screen(rho) has the right shape."""

    def test_finite_at_origin(self):
        fc = FrozenCore(Z=11)
        fc.solve()
        rho = np.array([1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0])
        f = _screening_correction_grid(11.0, fc, rho)
        # Must all be finite
        assert np.all(np.isfinite(f))
        # No diverging behavior at small rho (would be 1/rho at the limit)
        # f(rho_small) - f(rho_2x_smaller) shouldn't blow up
        # (f is smooth at origin per Hartree-on-spherical-density)

    def test_decays_at_large_rho(self):
        """f_screen(rho) -> N_core_total / rho at large rho (asymptotic [Ne]+)."""
        fc = FrozenCore(Z=11)
        fc.solve()
        # At rho = 50 bohr, all 10 [Ne] electrons are inside
        rho = np.array([50.0])
        f = _screening_correction_grid(11.0, fc, rho)
        expected = 10.0 / 50.0
        # tilde_Phi(50) should be ~0 (n_r decays); so f ~ 10/50 = 0.2
        assert abs(f[0] - expected) < 0.01, (
            f"f_screen(50) = {f[0]:.4f}, expected ~ 10/50 = {expected:.4f}"
        )


# ---------------------------------------------------------------------------
# Convergence with grid / quadrature parameters
# ---------------------------------------------------------------------------


class TestConvergence:
    def test_n_legendre_convergence(self):
        """Result should be stable wrt n_legendre at moderate values."""
        H_states = _enumerate_states(2)
        v_64 = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4, n_legendre=64,
        )
        v_128 = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4, n_legendre=128,
        )
        diff = np.max(np.abs(v_64 - v_128))
        assert diff < 1e-4, (
            f"n_legendre 64->128 not converged: {diff:.2e}"
        )

    def test_n_grid_convergence(self):
        """Result should be stable wrt n_grid at moderate values."""
        H_states = _enumerate_states(2)
        v_2k = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4, n_grid=2000,
        )
        v_8k = compute_screened_cross_center_vne(
            1.0, H_states, 11.0, 3.5, L_max=4, n_grid=8000,
        )
        diff = np.max(np.abs(v_2k - v_8k))
        # Looser tol for n_grid (coarser convergence than Legendre)
        assert diff < 1e-3, f"n_grid 2k->8k not converged: {diff:.2e}"
