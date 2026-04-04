"""
Tests for the IBM Quantum VQE Demo
===================================

All tests use local simulator only -- no IBM account or network required.
Marked slow because VQE optimization on 10-qubit systems takes 10-30s.

Author: GeoVac Development Team
Date: April 2026
"""

import subprocess
import sys
import warnings

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _import_demo():
    """Import demo module functions."""
    sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent / 'demo'))
    import ibm_quantum_demo as demo
    return demo


# ---------------------------------------------------------------------------
# Unit tests (fast)
# ---------------------------------------------------------------------------

class TestHamiltonianBuilders:
    """Test Hamiltonian construction (fast, no VQE)."""

    def test_build_geovac_h2(self) -> None:
        """GeoVac H2 bond-pair Hamiltonian has expected structure."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, info = demo.build_geovac_h2()

        assert info['n_qubits'] == 10
        assert info['n_terms'] == 112
        assert abs(info['one_norm'] - 8.17) < 0.5
        assert spo.num_qubits == 10

    def test_build_gaussian_h2(self) -> None:
        """Gaussian STO-3G H2 has expected structure."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, info = demo.build_gaussian_h2()

        assert info['n_qubits'] == 4
        assert info['n_terms'] == 15
        assert spo.num_qubits == 4

    def test_exact_ground_state_gaussian(self) -> None:
        """Exact diag of Gaussian STO-3G matches literature."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, _ = demo.build_gaussian_h2()
        exact = demo.exact_ground_state(spo)

        # Literature: -1.1373 Ha for H2 STO-3G at R=1.4
        assert abs(exact - (-1.1373)) < 0.001

    def test_exact_ground_state_geovac(self) -> None:
        """Exact diag of GeoVac H2 returns a finite energy."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, _ = demo.build_geovac_h2()
        exact = demo.exact_ground_state(spo)

        assert np.isfinite(exact)
        # GeoVac bond-pair encoding uses a single-center reference;
        # the absolute energy can be positive. Just verify it is finite.
        assert abs(exact) < 100  # sanity bound

    def test_qwc_groups_geovac(self) -> None:
        """GeoVac H2 has ~21 QWC groups."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, info = demo.build_geovac_h2()

        assert info['n_qwc_groups'] is not None
        assert info['n_qwc_groups'] == 21


# ---------------------------------------------------------------------------
# VQE tests (slow -- require optimization)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestVQEStatevector:
    """Test statevector VQE convergence."""

    def test_gaussian_vqe_converges(self) -> None:
        """Gaussian STO-3G VQE converges within 1 mHa (4 qubits, easy)."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, info = demo.build_gaussian_h2()

        exact = demo.exact_ground_state(spo)
        result = demo.run_vqe_statevector(
            spo, n_restarts=5, maxiter=1000, verbose=False,
        )

        error_mha = abs(result['energy'] - exact) * 1000
        assert error_mha < 1.0, f"Gaussian VQE error {error_mha:.3f} mHa > 1.0 mHa"

    def test_geovac_vqe_reasonable(self) -> None:
        """GeoVac H2 VQE produces a reasonable result (10 qubits, 80 params).

        The 10-qubit bond-pair system with 80 variational parameters is
        a harder optimization landscape than 4-qubit STO-3G. With 5 restarts
        and 1000 iterations, we expect < 20 mHa error. Sub-mHa convergence
        would require many more restarts or a better optimizer.
        """
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, info = demo.build_geovac_h2()

        exact = demo.exact_ground_state(spo)
        result = demo.run_vqe_statevector(
            spo, n_restarts=5, maxiter=1000, verbose=False,
        )

        error_mha = abs(result['energy'] - exact) * 1000
        assert error_mha < 20.0, f"GeoVac VQE error {error_mha:.3f} mHa > 20 mHa"

    def test_vqe_returns_valid_result(self) -> None:
        """VQE result dict has all expected keys."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, _ = demo.build_gaussian_h2()

        result = demo.run_vqe_statevector(
            spo, n_restarts=1, maxiter=100, verbose=False,
        )

        assert 'energy' in result
        assert 'wall_time' in result
        assert 'method' in result
        assert result['method'] == 'statevector'
        assert np.isfinite(result['energy'])
        assert result['wall_time'] > 0


@pytest.mark.slow
class TestVQEShotBased:
    """Test shot-based VQE."""

    def test_shot_based_reasonable(self) -> None:
        """Shot-based VQE on Gaussian STO-3G within 10 mHa."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spo, _ = demo.build_gaussian_h2()

        exact = demo.exact_ground_state(spo)
        result = demo.run_vqe_shots(
            spo, shots=10000, n_restarts=3, maxiter=500, verbose=False,
        )

        error_mha = abs(result['energy'] - exact) * 1000
        # Shot noise makes convergence harder; 10 mHa is a reasonable bar
        assert error_mha < 10.0, f"Shot-based VQE error {error_mha:.3f} mHa > 10 mHa"
        assert result['shots'] == 10000


@pytest.mark.slow
class TestCompareGaussian:
    """Test the comparison workflow."""

    def test_comparison_runs(self) -> None:
        """Both Hamiltonians can be built and compared."""
        demo = _import_demo()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geo_spo, geo_info = demo.build_geovac_h2()
            gau_spo, gau_info = demo.build_gaussian_h2()

        geo_exact = demo.exact_ground_state(geo_spo)
        gau_exact = demo.exact_ground_state(gau_spo)

        geo_vqe = demo.run_vqe_statevector(
            geo_spo, n_restarts=3, maxiter=500, verbose=False,
        )
        gau_vqe = demo.run_vqe_statevector(
            gau_spo, n_restarts=3, maxiter=500, verbose=False,
        )

        # Both should produce finite energies
        assert np.isfinite(geo_vqe['energy'])
        assert np.isfinite(gau_vqe['energy'])

        # GeoVac has more qubits and terms
        assert geo_info['n_qubits'] > gau_info['n_qubits']
        assert geo_info['n_terms'] > gau_info['n_terms']


# ---------------------------------------------------------------------------
# CLI integration test
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestCLI:
    """Test the command-line interface."""

    def test_cli_simulator_help(self) -> None:
        """Script --help runs without error."""
        result = subprocess.run(
            [sys.executable, 'demo/ibm_quantum_demo.py', '--help'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert 'GeoVac' in result.stdout
