"""
Tests for ecosystem_export.py — full library validation.

Verifies:
  1. All 30 systems build successfully via hamiltonian()
  2. Known Pauli counts match expected values
  3. Universal Pauli/Q ratios hold
  4. Isostructural invariance across rows

Author: GeoVac Development Team
Date: April 2026
"""

import sys
from pathlib import Path

import pytest

# Ensure local project root is on sys.path
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from geovac.ecosystem_export import hamiltonian


# ---------------------------------------------------------------------------
# Expected Pauli counts (non-identity terms) at max_n=2
# ---------------------------------------------------------------------------
EXPECTED_PAULI = {
    # Atomic / diatomic
    'He':   119,
    'H2':   111,
    # First-row main-group
    'LiH':  333,
    'BeH2': 555,
    'CH4':  999,
    'NH3':  888,
    'H2O':  777,
    'HF':   666,
    # Second-row
    'NaH':  222,
    'MgH2': 444,
    'SiH4': 888,
    'PH3':  777,
    'H2S':  666,
    'HCl':  555,
    # Third-row s-block
    'KH':   222,
    'CaH2': 444,
    # Third-row p-block
    'GeH4': 888,
    'AsH3': 777,
    'H2Se': 666,
    'HBr':  555,
    # TM hydrides (all identical)
    'ScH':  277, 'TiH':  277, 'VH':   277, 'CrH':  277, 'MnH':  277,
    'FeH':  277, 'CoH':  277, 'NiH':  277, 'CuH':  277, 'ZnH':  277,
}


# ---------------------------------------------------------------------------
# Test: all systems build successfully
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("system", list(EXPECTED_PAULI.keys()))
def test_system_builds(system: str) -> None:
    """Each system should build without error."""
    H = hamiltonian(system, verbose=False)
    assert H.n_qubits > 0
    assert H.n_terms > 0
    assert H.one_norm > 0


# ---------------------------------------------------------------------------
# Test: Pauli counts match expected values
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("system,expected", list(EXPECTED_PAULI.items()))
def test_pauli_counts(system: str, expected: int) -> None:
    """Non-identity Pauli count must match expected value exactly."""
    H = hamiltonian(system, verbose=False)
    # n_terms includes identity; non-identity = n_terms - 1
    n_pauli = H.n_terms - 1
    assert n_pauli == expected, (
        f"{system}: expected {expected} non-identity Pauli terms, got {n_pauli}"
    )


# ---------------------------------------------------------------------------
# Test: universal Pauli/Q ratio for main-group hydrides
# ---------------------------------------------------------------------------
MAIN_GROUP = [
    'LiH', 'BeH2', 'CH4', 'NH3', 'H2O', 'HF',
    'NaH', 'MgH2', 'SiH4', 'PH3', 'H2S', 'HCl',
    'KH', 'CaH2', 'GeH4', 'AsH3', 'H2Se', 'HBr',
    'H2',
]

@pytest.mark.parametrize("system", MAIN_GROUP)
def test_main_group_ratio(system: str) -> None:
    """All main-group hydrides should have Pauli/Q = 11.10."""
    H = hamiltonian(system, verbose=False)
    n_pauli = H.n_terms - 1
    ratio = n_pauli / H.n_qubits
    assert abs(ratio - 11.10) < 0.01, (
        f"{system}: Pauli/Q = {ratio:.2f}, expected 11.10"
    )


# ---------------------------------------------------------------------------
# Test: TM hydride ratio
# ---------------------------------------------------------------------------
TM_HYDRIDES = ['ScH', 'TiH', 'VH', 'CrH', 'MnH',
               'FeH', 'CoH', 'NiH', 'CuH', 'ZnH']

@pytest.mark.parametrize("system", TM_HYDRIDES)
def test_tm_ratio(system: str) -> None:
    """All TM hydrides should have Pauli/Q = 9.23."""
    H = hamiltonian(system, verbose=False)
    n_pauli = H.n_terms - 1
    ratio = n_pauli / H.n_qubits
    assert abs(ratio - 9.23) < 0.01, (
        f"{system}: Pauli/Q = {ratio:.2f}, expected 9.23"
    )


# ---------------------------------------------------------------------------
# Test: isostructural invariance
# ---------------------------------------------------------------------------
ISOSTRUCTURAL_GROUPS = [
    # Same block topology → same Pauli count (and same Q)
    (['NaH', 'KH'], 222, 20),           # alkali hydrides, no core block
    (['MgH2', 'CaH2'], 444, 40),        # alkaline earth hydrides
    (['SiH4', 'GeH4'], 888, 80),        # group 14 tetrahydrides
    (['PH3', 'AsH3'], 777, 70),         # group 15 trihydrides
    (['H2S', 'H2Se'], 666, 60),         # group 16 dihydrides
    (['HCl', 'HBr'], 555, 50),          # group 17 monohydrides
]

@pytest.mark.parametrize("systems,expected_pauli,expected_Q", ISOSTRUCTURAL_GROUPS)
def test_isostructural_invariance(
    systems: list, expected_pauli: int, expected_Q: int,
) -> None:
    """Frozen-core molecules with same block topology must have identical
    Pauli counts and qubit counts across rows."""
    for system in systems:
        H = hamiltonian(system, verbose=False)
        n_pauli = H.n_terms - 1
        assert n_pauli == expected_pauli, (
            f"{system}: N_pauli={n_pauli}, expected {expected_pauli}"
        )
        assert H.n_qubits == expected_Q, (
            f"{system}: Q={H.n_qubits}, expected {expected_Q}"
        )


# ---------------------------------------------------------------------------
# Test: TM hydride isostructural invariance (all 10 identical)
# ---------------------------------------------------------------------------
def test_tm_isostructural() -> None:
    """All 10 TM hydrides must have identical Q=30 and N_pauli=277."""
    for system in TM_HYDRIDES:
        H = hamiltonian(system, verbose=False)
        n_pauli = H.n_terms - 1
        assert H.n_qubits == 30, f"{system}: Q={H.n_qubits}"
        assert n_pauli == 277, f"{system}: N_pauli={n_pauli}"
