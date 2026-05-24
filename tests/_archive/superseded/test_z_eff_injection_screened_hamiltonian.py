"""
ARCHIVED 2026-05-23 (Cleanup Track B): TestBackwardCompatibility extracted from
tests/test_z_eff_injection.py. Both tests depend on
compute_nuclear_coupling_screened (gone from geovac/level4_multichannel.py) and
on the Z_A_func kwarg of build_angular_hamiltonian (also gone). Both were
removed in v2.7.0 (commit 8d692a0). The Z_A_func code path now lives only in
geovac/rho_collapse_cache.py::AngularCache, which is exercised by the live
TestCallableSmoke, TestZEffShift, and TestIntegrationLiH classes that remain in
tests/test_z_eff_injection.py.
"""

import numpy as np
import pytest

from geovac.level4_multichannel import (
    build_angular_hamiltonian,
    compute_nuclear_coupling,
    compute_nuclear_coupling_screened,
    _channel_list,
)


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
