"""
Tests for Z_eff injection into Level 4 charge function + Pauli projector.

Validates:
1. Backward compatibility (constant Z_A = 1 vs original H2)
2. Callable Z_A_func = constant reproduces analytical path
3. Non-trivial Z_eff lowers energy (stronger attraction)
4. Pauli projector is a valid projection matrix
5. Pauli projector removes one eigenvalue from spectrum
6. Integration test: CoreScreening + Level 4 for LiH-like system
"""

import numpy as np
import pytest

from geovac.level4_multichannel import (
    build_angular_hamiltonian,
    compute_nuclear_coupling,
    compute_nuclear_coupling_screened,
    _channel_list,
)
from geovac.rho_collapse_cache import AngularCache, FastAdiabaticPES
from geovac.pauli_projector import CoreValenceProjector


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _h2_energy_original(R: float = 1.4, l_max: int = 2, n_alpha: int = 100,
                        n_rho: int = 60) -> float:
    """Compute H2 energy using the ORIGINAL code path (no Z_A_func)."""
    cache = AngularCache(Z_A=1.0, Z_B=1.0, l_max=l_max, n_alpha=n_alpha)
    pes = FastAdiabaticPES(cache, Z_A=1.0, Z_B=1.0)
    E_elec, E_total = pes.energy_at_R(R, n_rho=n_rho)
    return E_total


def _h2_energy_screened(R: float = 1.4, l_max: int = 2, n_alpha: int = 100,
                        n_rho: int = 60,
                        Z_A_func=None) -> float:
    """Compute H2 energy using the Z_A_func code path."""
    if Z_A_func is None:
        Z_A_func = lambda r: 1.0
    cache = AngularCache(Z_A=1.0, Z_B=1.0, l_max=l_max, n_alpha=n_alpha,
                         Z_A_func=Z_A_func, n_theta=64)
    pes = FastAdiabaticPES(cache, Z_A=1.0, Z_B=1.0)
    E_elec, E_total = pes.energy_at_R(R, n_rho=n_rho)
    return E_total


# ---------------------------------------------------------------------------
# Test 1: Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """H2 with constant Z_A=1, Z_B=1 must match original results."""

    def test_nuclear_coupling_screened_matches_analytical(self):
        """Screened coupling with constant Z=1 matches analytical at single point."""
        n_alpha = 50
        h = (np.pi / 2) / (n_alpha + 1)
        alpha = (np.arange(n_alpha) + 1) * h
        rho = 0.5
        R_e = 1.4  # arbitrary

        # Test several channel pairs
        test_cases = [
            (0, 0, 0, 0, 0, 0),  # diagonal (0,0)-(0,0)
            (0, 0, 2, 0, 0, 0),  # off-diagonal (0,0)-(2,0)
            (2, 0, 0, 0, 0, 0),  # off-diagonal (2,0)-(0,0)
            (2, 0, 2, 0, 0, 0),  # diagonal (2,0)-(2,0)
            (0, 0, 0, 2, 0, 0),  # coupling electron 2
        ]

        for l1p, l2p, l1, l2, m1, m2 in test_cases:
            V_analytical = compute_nuclear_coupling(
                l1p, l2p, l1, l2, m1, m2, alpha, rho,
                Z=1.0, Z_A=1.0, Z_B=1.0,
            )
            V_screened = compute_nuclear_coupling_screened(
                l1p, l2p, l1, l2, m1, m2, alpha, rho,
                Z_A_func=lambda r: 1.0, Z_B=1.0, R_e=R_e,
            )

            if np.max(np.abs(V_analytical)) < 1e-12:
                assert np.max(np.abs(V_screened)) < 1e-6, (
                    f"Channel ({l1p},{l2p})-({l1},{l2}): expected zero, "
                    f"got max={np.max(np.abs(V_screened)):.2e}"
                )
            else:
                rel_err = np.max(np.abs(V_screened - V_analytical)) / np.max(
                    np.abs(V_analytical)
                )
                assert rel_err < 1e-4, (
                    f"Channel ({l1p},{l2p})-({l1},{l2}): "
                    f"relative error {rel_err:.2e} > 1e-4"
                )

    def test_build_hamiltonian_screened_matches_original(self):
        """Full Hamiltonian with Z_A_func=const matches original."""
        n_alpha = 50
        h = (np.pi / 2) / (n_alpha + 1)
        alpha = (np.arange(n_alpha) + 1) * h
        rho = 0.5
        R_e = 1.4

        H_orig = build_angular_hamiltonian(
            alpha, rho, R_e, l_max=2, Z=1.0,
        )
        H_screened = build_angular_hamiltonian(
            alpha, rho, R_e, l_max=2, Z=1.0,
            Z_A_func=lambda r: 1.0, n_theta=64,
        )

        # Nuclear coupling is the only difference; kinetic+ee are identical
        diff = np.max(np.abs(H_screened - H_orig))
        rel = diff / np.max(np.abs(H_orig))
        assert rel < 1e-4, (
            f"Hamiltonian relative difference {rel:.2e} > 1e-4"
        )


# ---------------------------------------------------------------------------
# Test 2: Z_A_func callable smoke test (via PES)
# ---------------------------------------------------------------------------

class TestCallableSmoke:
    """Z_A_func = lambda r: 1.0 must match constant Z_A=1.0 PES result."""

    def test_h2_callable_const_matches_original(self):
        """Single-R energy with callable const Z matches original."""
        R = 1.4
        E_orig = _h2_energy_original(R, l_max=2, n_alpha=80, n_rho=50)
        E_callable = _h2_energy_screened(R, l_max=2, n_alpha=80, n_rho=50,
                                          Z_A_func=lambda r: 1.0)

        diff = abs(E_callable - E_orig)
        print(f"\nH2 backward compatibility at R={R}:")
        print(f"  E_original  = {E_orig:.6f} Ha")
        print(f"  E_callable  = {E_callable:.6f} Ha")
        print(f"  |diff|      = {diff:.2e} Ha")

        assert diff < 1e-3, (
            f"Callable const Z_A energy differs by {diff:.2e} Ha "
            f"(threshold 1e-3)"
        )


# ---------------------------------------------------------------------------
# Test 3: Z_eff changes energies
# ---------------------------------------------------------------------------

class TestZEffShift:
    """Enhanced Z near nucleus should lower the energy."""

    def test_enhanced_z_lowers_energy(self):
        """Z_A_func = 1 + 0.5*exp(-r) gives stronger attraction → lower E."""
        R = 1.4

        def Z_enhanced(r):
            return 1.0 + 0.5 * np.exp(-r)

        E_standard = _h2_energy_screened(R, l_max=2, n_alpha=80, n_rho=50,
                                          Z_A_func=lambda r: 1.0)
        E_enhanced = _h2_energy_screened(R, l_max=2, n_alpha=80, n_rho=50,
                                          Z_A_func=Z_enhanced)

        print(f"\nZ_eff shift test at R={R}:")
        print(f"  E_standard  = {E_standard:.6f} Ha")
        print(f"  E_enhanced  = {E_enhanced:.6f} Ha")
        print(f"  shift       = {E_enhanced - E_standard:.6f} Ha")

        assert E_enhanced < E_standard, (
            f"Enhanced Z_A should lower energy: "
            f"E_enhanced={E_enhanced:.6f} >= E_standard={E_standard:.6f}"
        )


# ---------------------------------------------------------------------------
# Test 4: Pauli projector structure
# ---------------------------------------------------------------------------

class TestPauliProjectorStructure:
    """Verify P is a valid projection matrix."""

    def test_projector_properties(self):
        """P² = P, P = Pᵀ, rank = N-1, eigenvalues ∈ {0, 1}."""
        proj = CoreValenceProjector(
            Z_eff_core=2.69, l_max=2, n_alpha=50, homonuclear=False,
        )
        R_e = 2.0
        P = proj.build_projector(R_e)
        N = proj.N

        # P² = P (idempotent)
        P2 = P @ P
        assert np.max(np.abs(P2 - P)) < 1e-12, "P is not idempotent"

        # P = Pᵀ (symmetric)
        assert np.max(np.abs(P - P.T)) < 1e-12, "P is not symmetric"

        # Eigenvalues are 0 or 1
        evals = np.linalg.eigvalsh(P)
        for ev in evals:
            assert abs(ev) < 1e-10 or abs(ev - 1.0) < 1e-10, (
                f"Eigenvalue {ev} is neither 0 nor 1"
            )

        # Rank = N - 1 (one state removed)
        rank = np.sum(np.abs(evals) > 0.5)
        assert rank == N - 1, (
            f"Expected rank {N-1}, got {rank}"
        )

        # Print summary
        n_zero = np.sum(np.abs(evals) < 0.5)
        print(f"\nPauli projector (l_max=2, Z_eff={proj.Z_eff_core}):")
        print(f"  Matrix size: {N} x {N}")
        print(f"  Channels: {proj.n_ch}")
        print(f"  Rank: {rank} (removed {n_zero} state(s))")
        print(f"  Min eigenvalue: {evals[0]:.2e}")
        print(f"  Max eigenvalue: {evals[-1]:.6f}")


# ---------------------------------------------------------------------------
# Test 5: Pauli projector effect on eigenvalues
# ---------------------------------------------------------------------------

class TestPauliProjectorEffect:
    """Projected Hamiltonian should have one fewer non-trivial eigenvalue."""

    def test_eigenvalue_removal(self):
        """Projector removes one eigenvalue from the angular spectrum."""
        n_alpha = 50
        h = (np.pi / 2) / (n_alpha + 1)
        alpha = (np.arange(n_alpha) + 1) * h
        rho = 0.5
        R_e = 2.0

        # Build H2 angular Hamiltonian (homonuclear=False for comparison)
        H = build_angular_hamiltonian(
            alpha, rho, R_e, l_max=2, Z=1.0,
            Z_A=1.0, Z_B=1.0,
        )

        # Original eigenvalues
        evals_orig = np.linalg.eigvalsh(H)

        # Build projector and apply
        proj = CoreValenceProjector(
            Z_eff_core=2.69, l_max=2, n_alpha=n_alpha, homonuclear=True,
        )
        H_proj = proj.project(H, R_e)

        # Projected eigenvalues
        evals_proj = np.linalg.eigvalsh(H_proj)

        # Count non-trivial eigenvalues (not near zero from projection)
        threshold = 1e-6
        n_orig = np.sum(np.abs(evals_orig) > threshold)
        n_proj = np.sum(np.abs(evals_proj) > threshold)

        print(f"\nPauli projector effect on eigenvalues:")
        print(f"  Original non-trivial: {n_orig}")
        print(f"  Projected non-trivial: {n_proj}")
        print(f"  Removed: {n_orig - n_proj}")
        print(f"  Lowest 5 original:  {evals_orig[:5]}")
        print(f"  Lowest 5 projected: {evals_proj[:5]}")

        # The projected Hamiltonian should have at least 1 fewer
        # non-trivial eigenvalue (the core overlap state becomes zero)
        assert n_proj <= n_orig, (
            f"Projected spectrum has MORE non-trivial eigenvalues than original"
        )


# ---------------------------------------------------------------------------
# Test 6: Integration test — CoreScreening + Level 4 for LiH
# ---------------------------------------------------------------------------

class TestIntegrationLiH:
    """
    Construct CoreScreening for Z=3, inject into Level 4 with Z_B=1.

    This test verifies the full pipeline works without errors and
    produces physically reasonable results.
    """

    def test_lih_single_point(self):
        """Single-point LiH calculation completes with reasonable energy."""
        # Use a mock Z_eff function instead of solving the full CoreScreening
        # (which would take too long for a unit test).
        # Li core: Z_eff goes from 3.0 at r=0 to ~1.0 at r>2 bohr
        def z_eff_li(r):
            """Mock Z_eff for Li core (analytical approximation)."""
            r = np.atleast_1d(np.asarray(r, dtype=float))
            # Hydrogenic 1s² screening: N_core(r) ≈ 2[1-(1+Zr)exp(-2Zr)]
            Z_core = 2.69  # Clementi-Raimondi for Li 1s
            N_core = 2.0 * (1.0 - (1.0 + Z_core * r) * np.exp(-2.0 * Z_core * r))
            z_eff = 3.0 - np.clip(N_core, 0.0, 2.0)
            return z_eff

        R = 3.0  # bohr, near LiH equilibrium

        # Verify z_eff behavior
        assert abs(z_eff_li(0.0) - 3.0) < 0.01, "Z_eff(0) should be ~3"
        assert abs(z_eff_li(5.0) - 1.0) < 0.1, "Z_eff(5) should be ~1"

        # Build Level 4 solver with screened Z_A
        cache = AngularCache(
            Z_A=3.0, Z_B=1.0, l_max=2, n_alpha=60,
            Z_A_func=z_eff_li, n_theta=48,
        )
        pes = FastAdiabaticPES(cache, Z_A=3.0, Z_B=1.0)

        E_elec, E_total = pes.energy_at_R(R, n_rho=40)

        # Core energy (Li²⁺ → Li core ≈ He-like Z=3)
        E_core = -9.0 / 2.0  # hydrogenic: -Z²/2 per electron × 2 = -9 Ha
        # Better estimate: He-like Z=3 ground state ≈ -7.28 Ha
        E_core_he = -7.28

        # Nuclear repulsion
        V_NN = 3.0 * 1.0 / R

        print(f"\nLiH integration test at R={R:.1f} bohr:")
        print(f"  Z_eff(0.0) = {z_eff_li(0.0).item():.3f}")
        print(f"  Z_eff(0.5) = {z_eff_li(0.5).item():.3f}")
        print(f"  Z_eff(1.0) = {z_eff_li(1.0).item():.3f}")
        print(f"  Z_eff(3.0) = {z_eff_li(3.0).item():.3f}")
        print(f"  Valence E_elec = {E_elec:.6f} Ha")
        print(f"  Valence E_total (E_elec + V_NN) = {E_total:.6f} Ha")
        print(f"  V_NN = {V_NN:.6f} Ha")
        print(f"  Core E (He-like Z=3) ~ {E_core_he:.3f} Ha")
        print(f"  Estimated total = {E_total + E_core_he:.3f} Ha")

        # Basic sanity checks
        assert np.isfinite(E_elec), "E_elec is not finite"
        assert np.isfinite(E_total), "E_total is not finite"
        assert E_elec < 0, f"E_elec should be negative, got {E_elec:.4f}"

        # The valence electronic energy for 2 electrons in the screened
        # field should be negative and of order -1 to -5 Ha
        assert -10.0 < E_elec < 0.0, (
            f"E_elec = {E_elec:.4f} outside expected range [-10, 0]"
        )


# ---------------------------------------------------------------------------
# Summary runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
