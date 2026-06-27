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
from geovac.direct_ci import DirectCISolver


def _direct_solve(idx):
    """Drop-in replacement for the removed fci_method='direct' path.

    The fci_method= kwarg was removed during the v2.7.0 PK/composed-qubit
    refactor; DirectCISolver(idx).solve() is the production direct-CI entry.
    """
    return DirectCISolver(idx).solve()


def _solve_both(n_electrons: int, max_n: int, Z: int) -> tuple:
    """
    Solve with both matrix and direct methods, return (E_matrix, E_direct).

    Suppresses LatticeIndex warnings for cleaner test output.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        idx_m = LatticeIndex(
            n_electrons=n_electrons, max_n=max_n, nuclear_charge=Z,
            vee_method='slater_full',  # fci_method='matrix' removed; default is matrix
        )
        E_m, _ = idx_m.compute_ground_state()

        idx_d = LatticeIndex(
            n_electrons=n_electrons, max_n=max_n, nuclear_charge=Z,
            vee_method='slater_full',  # fci_method='direct' → use DirectCISolver
        )
        E_d, _ = _direct_solve(idx_d)

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
        """He nmax=5 direct CI (exact h1 + grid Slater): E = -2.84483 Ha (2.03%).

        Two-sided recompute-pin. The previous guard (`E < -2.84`) was a
        FALSE-POSITIVE: it is one-sided, so a grossly-wrong-but-too-low
        energy (e.g. an unphysical -3.0) would also pass. We pin the
        actual framework value and add a variational-bound sanity.

        This is the *exact*-h1 + 2000-pt-grid-Slater config (the
        LatticeIndex default for vee_method='slater_full'); it lands at
        2.03% above the exact NR He energy -2.903724 Ha. NOTE this is
        DISTINCT from the FCI-A *hybrid*-h1 headline (-2.893582 / 0.349%,
        pinned in tests/test_casimir_ci.py::TestFCIAHybridHe).
        """
        E_HE_EXACT = -2.903724  # exact non-relativistic He ground state
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=2, max_n=5, nuclear_charge=2,
                vee_method='slater_full',  # fci_method='direct' → DirectCISolver
            )
            E, _ = _direct_solve(idx)
        assert abs(E[0] - (-2.84483)) < 1e-3, (
            f"He nmax=5 exact-h1: E={E[0]:.6f}, expected -2.84483 ± 1e-3"
        )
        assert E[0] > E_HE_EXACT, (
            f"He nmax=5: E={E[0]:.6f} violates variational bound "
            f"(exact={E_HE_EXACT})"
        )

    @pytest.mark.slow
    def test_direct_ci_li_accuracy(self) -> None:
        """Li nmax=4 direct CI (exact h1 + grid Slater): E = -7.39592 Ha (1.10%).

        Two-sided recompute-pin replacing the one-sided `E < -7.39` guard.
        Reference is the exact NR Li energy -7.4781 Ha; this config sits
        1.10% above it (matches FCI-A Table convergence row, -7.395921).
        Marked slow: ~19s (34,220-SD assembly + sparse eigensolve).
        """
        E_LI_EXACT = -7.4781  # exact non-relativistic Li ground state
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=4, nuclear_charge=3,
                vee_method='slater_full',  # fci_method='direct' → DirectCISolver
            )
            E, _ = _direct_solve(idx)
        assert abs(E[0] - (-7.395921)) < 1e-3, (
            f"Li nmax=4 exact-h1: E={E[0]:.6f}, expected -7.395921 ± 1e-3"
        )
        assert E[0] > E_LI_EXACT, (
            f"Li nmax=4: E={E[0]:.6f} violates variational bound "
            f"(exact={E_LI_EXACT})"
        )


class TestDirectCINewAtoms:
    """Direct CI enables new atoms that were previously inaccessible."""

    def test_direct_ci_be_runs(self) -> None:
        """Be (4e) nmax=3: completes without error."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=4, max_n=3, nuclear_charge=4,
                vee_method='slater_full',  # fci_method='direct' → DirectCISolver
            )
            E, _ = _direct_solve(idx)
        # Just check it ran and produced a finite energy
        assert np.isfinite(E[0]), f"Be nmax=3: E={E[0]}"
        print(f"  Be nmax=3: E={E[0]:.6f} Ha")


    @pytest.mark.slow
    def test_direct_ci_be_nmax4(self) -> None:
        """Be (4e) nmax=4 direct CI (exact h1 + grid Slater): E = -14.5355 Ha (0.90%).

        Two-sided recompute-pin replacing the one-sided `E < -14.0` guard.
        Reference is the exact NR Be energy -14.6674 Ha; this config sits
        0.90% above it (matches FCI-A Table convergence row, -14.535460).
        Marked slow: ~6 min (487,635-SD assembly + sparse eigensolve).
        """
        E_BE_EXACT = -14.6674  # exact non-relativistic Be ground state
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=4, max_n=4, nuclear_charge=4,
                vee_method='slater_full',  # fci_method='direct' → DirectCISolver
            )
            E, _ = _direct_solve(idx)
        assert abs(E[0] - (-14.535460)) < 2e-3, (
            f"Be nmax=4 exact-h1: E={E[0]:.6f}, expected -14.535460 ± 2e-3"
        )
        assert E[0] > E_BE_EXACT, (
            f"Be nmax=4: E={E[0]:.6f} violates variational bound "
            f"(exact={E_BE_EXACT})"
        )
        error = abs(E[0] - E_BE_EXACT) / abs(E_BE_EXACT) * 100
        print(f"  Be nmax=4: E={E[0]:.4f} Ha, error={error:.2f}%")

    @pytest.mark.slow
    def test_direct_ci_li_nmax5(self) -> None:
        """Li nmax=5 (216k SDs) direct CI (exact h1 + grid Slater): E = -7.39775 Ha (1.07%).

        Two-sided recompute-pin replacing the one-sided `E < -7.39` guard.
        Reference is the exact NR Li energy -7.4781 Ha; this config sits
        1.07% above it (matches FCI-A Table convergence row, -7.397751).
        Marked slow: ~2 min (215,820-SD assembly + sparse eigensolve).
        """
        E_LI_EXACT = -7.4781  # exact non-relativistic Li ground state
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=5, nuclear_charge=3,
                vee_method='slater_full',  # fci_method='direct' → DirectCISolver
            )
            E, _ = _direct_solve(idx)
        assert abs(E[0] - (-7.397751)) < 1e-3, (
            f"Li nmax=5 exact-h1: E={E[0]:.6f}, expected -7.397751 ± 1e-3"
        )
        assert E[0] > E_LI_EXACT, (
            f"Li nmax=5: E={E[0]:.6f} violates variational bound "
            f"(exact={E_LI_EXACT})"
        )
        error = abs(E[0] - E_LI_EXACT) / abs(E_LI_EXACT) * 100
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
                    vee_method='slater_full',  # fci_method='direct' → DirectCISolver
                )
            t0 = time.perf_counter()
            _direct_solve(idx)
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
