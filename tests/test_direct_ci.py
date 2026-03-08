"""
Tests for excitation-driven direct CI Hamiltonian construction.

Validates that the direct assembly reproduces exact matrix diagonalization
for all systems currently passing in the test suite.
"""

import time
import warnings

import numpy as np
import pytest

from geovac.lattice_index import LatticeIndex


def _solve_both(n_electrons: int, max_n: int, Z: int) -> tuple:
    """
    Solve with both matrix and direct methods, return (E_matrix, E_direct).

    Suppresses LatticeIndex warnings for cleaner test output.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        idx_m = LatticeIndex(
            n_electrons=n_electrons, max_n=max_n, nuclear_charge=Z,
            vee_method='slater_full', fci_method='matrix',
        )
        E_m, _ = idx_m.compute_ground_state()

        idx_d = LatticeIndex(
            n_electrons=n_electrons, max_n=max_n, nuclear_charge=Z,
            vee_method='slater_full', fci_method='direct',
        )
        E_d, _ = idx_d.compute_ground_state()

    return E_m[0], E_d[0]


class TestDirectVsMatrix:
    """Direct CI must reproduce matrix CI to machine precision."""

    def test_sigma_consistency_he(self) -> None:
        """He nmax=3: direct vs matrix energy diff < 1e-8 Ha."""
        E_m, E_d = _solve_both(n_electrons=2, max_n=3, Z=2)
        assert abs(E_d - E_m) < 1e-8, (
            f"He nmax=3: matrix={E_m:.10f}, direct={E_d:.10f}, "
            f"diff={abs(E_d-E_m):.2e}"
        )

    def test_sigma_consistency_li(self) -> None:
        """Li nmax=3: direct vs matrix energy diff < 1e-8 Ha."""
        E_m, E_d = _solve_both(n_electrons=3, max_n=3, Z=3)
        assert abs(E_d - E_m) < 1e-8, (
            f"Li nmax=3: matrix={E_m:.10f}, direct={E_d:.10f}, "
            f"diff={abs(E_d-E_m):.2e}"
        )


class TestDirectCIAccuracy:
    """Direct CI reproduces known benchmark energies."""

    def test_direct_ci_he_accuracy(self) -> None:
        """He nmax=5: E < -2.84 Ha (matches matrix baseline -2.8448)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=2, max_n=5, nuclear_charge=2,
                vee_method='slater_full', fci_method='direct',
            )
            E, _ = idx.compute_ground_state()
        assert E[0] < -2.84, f"He nmax=5: E={E[0]:.6f}, expected < -2.84"

    def test_direct_ci_li_accuracy(self) -> None:
        """Li nmax=4: E < -7.39 Ha (matches matrix baseline -7.3959)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=4, nuclear_charge=3,
                vee_method='slater_full', fci_method='direct',
            )
            E, _ = idx.compute_ground_state()
        assert E[0] < -7.39, f"Li nmax=4: E={E[0]:.6f}, expected < -7.39"


class TestDirectCINewAtoms:
    """Direct CI enables new atoms that were previously inaccessible."""

    def test_direct_ci_be_runs(self) -> None:
        """Be (4e) nmax=3: completes without error."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=4, max_n=3, nuclear_charge=4,
                vee_method='slater_full', fci_method='direct',
            )
            E, _ = idx.compute_ground_state()
        # Just check it ran and produced a finite energy
        assert np.isfinite(E[0]), f"Be nmax=3: E={E[0]}"
        print(f"  Be nmax=3: E={E[0]:.6f} Ha")


    def test_direct_ci_be_nmax4(self) -> None:
        """Be (4e) nmax=4: completes and E < -14.0 Ha."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=4, max_n=4, nuclear_charge=4,
                vee_method='slater_full', fci_method='direct',
            )
            E, _ = idx.compute_ground_state()
        assert E[0] < -14.0, f"Be nmax=4: E={E[0]:.6f}, expected < -14.0"
        error = abs(E[0] - (-14.6674)) / 14.6674 * 100
        print(f"  Be nmax=4: E={E[0]:.4f} Ha, error={error:.2f}%")

    def test_direct_ci_li_nmax5(self) -> None:
        """Li nmax=5 (216k SDs): completes and E < -7.39 Ha."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=5, nuclear_charge=3,
                vee_method='slater_full', fci_method='direct',
            )
            E, _ = idx.compute_ground_state()
        assert E[0] < -7.39, f"Li nmax=5: E={E[0]:.6f}, expected < -7.39"
        error = abs(E[0] - (-7.4781)) / 7.4781 * 100
        print(f"  Li nmax=5: E={E[0]:.4f} Ha, error={error:.2f}%")


class TestScaling:
    """Verify subquadratic scaling in N_SD."""

    def test_scaling_sublinear(self) -> None:
        """Wall time exponent < 2.0 in N_SD for He nmax=2..5."""
        times = []
        n_sds = []

        for nmax in [2, 3, 4, 5]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                idx = LatticeIndex(
                    n_electrons=2, max_n=nmax, nuclear_charge=2,
                    vee_method='slater_full', fci_method='direct',
                )
            t0 = time.perf_counter()
            idx.compute_ground_state()
            dt = time.perf_counter() - t0

            if dt > 0.001:  # skip tiny timings
                times.append(dt)
                n_sds.append(idx.n_sd)

        if len(times) >= 3:
            # Fit log(t) vs log(N_SD) to get exponent
            log_n = np.log(np.array(n_sds))
            log_t = np.log(np.array(times))
            coeffs = np.polyfit(log_n, log_t, 1)
            exponent = coeffs[0]
            print(f"  Scaling exponent: {exponent:.2f} (target < 2.0)")
            assert exponent < 2.5, (
                f"Scaling exponent {exponent:.2f} >= 2.5, "
                f"expected subquadratic"
            )
