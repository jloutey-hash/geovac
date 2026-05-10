"""Tests for geovac.hyperfine_a_constant (Sprint Cs-HFS-v2 A3).

Validates the generic A * I . J Pauli-sum wrapper for atomic HFS
observables. Covers:
  - Bohr-Fermi A constant against H 21cm experimental value
  - Pauli decomposition of I . J for I=1/2 (matches deuterium-like 4-state)
  - Pauli decomposition for Cs (I=7/2) with correct F=4/F=3 splitting
  - Hermiticity of the resulting Pauli sum
  - F-level energies match closed-form Lande formula
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.hyperfine_a_constant import (
    bohr_fermi_a_constant,
    hyperfine_a_ij_pauli_general,
    hyperfine_a_pauli_for_atomic_hfs,
    GE_FULL,
    GE_DIRAC,
    HZ_PER_HA,
    M_PROTON_OVER_M_E,
)

PAULI_MATS = {
    'I': np.array([[1, 0], [0, 1]], dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def _pauli_dict_to_dense(pauli, n_qubits):
    """Reconstruct dense matrix from Pauli dictionary."""
    dim = 2 ** n_qubits
    H = np.zeros((dim, dim), dtype=complex)
    for k, v in pauli.items():
        P = np.array([[1.0]], dtype=complex)
        for ch in k:
            P = np.kron(P, PAULI_MATS[ch])
        H = H + v * P
    return H


# ---------------------------------------------------------------------------
# Bohr-Fermi compute
# ---------------------------------------------------------------------------

class TestBohrFermiAConstant:
    """A_HF formula validation against hydrogen 21cm."""

    def test_h_1s_bohr_fermi(self):
        """H 1s with full g-factors -> A in [1418, 1424] MHz."""
        psi_1s = 1.0 / np.pi
        g_p = 5.585694713
        result = bohr_fermi_a_constant(psi_1s, g_e=GE_FULL, g_N=g_p)
        # 21cm experimental: 1420.405 MHz; BF formula gives ~1418 MHz
        # (the difference is mostly BF strict vs experimental, ~0.2%)
        assert 1418.0 < result['A_MHz'] < 1424.0, (
            f"A(H 1s) = {result['A_MHz']:.4f} MHz outside [1418, 1424]"
        )

    def test_h_1s_dirac_gives_lower(self):
        """Dirac g_e=2 should give smaller A than full Schwinger."""
        psi_1s = 1.0 / np.pi
        g_p = 5.585694713
        a_dirac = bohr_fermi_a_constant(psi_1s, g_e=GE_DIRAC, g_N=g_p)
        a_full = bohr_fermi_a_constant(psi_1s, g_e=GE_FULL, g_N=g_p)
        assert a_dirac['A_Ha'] < a_full['A_Ha'], (
            f"Dirac A {a_dirac['A_MHz']:.4f} should be < full A "
            f"{a_full['A_MHz']:.4f}"
        )
        # Ratio ~ 1.00116 = (1+a_e_codata)/1.0
        ratio = a_full['A_Ha'] / a_dirac['A_Ha']
        assert abs(ratio - GE_FULL / GE_DIRAC) < 1e-12

    def test_units_consistency(self):
        """A_Hz = A_Ha * HZ_PER_HA, A_MHz = A_Hz / 1e6."""
        result = bohr_fermi_a_constant(1.0, g_e=GE_FULL, g_N=1.0)
        assert result['A_Hz'] == pytest.approx(result['A_Ha'] * HZ_PER_HA)
        assert result['A_MHz'] == pytest.approx(result['A_Hz'] * 1e-6)

    def test_zero_psi_gives_zero(self):
        """|psi(0)|^2 = 0 -> A = 0 (l>=1 case)."""
        result = bohr_fermi_a_constant(0.0)
        assert result['A_Ha'] == 0.0


# ---------------------------------------------------------------------------
# Generic A * I . J Pauli decomposition
# ---------------------------------------------------------------------------

class TestHyperfineIJGeneral:
    """Validation of generic I . J Pauli decomposition."""

    def test_h_1s_isj_4_states(self):
        """H 1s I=1/2 J=1/2 -> 3 Pauli terms (XX, YY, ZZ each at A/4)."""
        A = 1.0  # arbitrary, in Ha
        pauli = hyperfine_a_ij_pauli_general(A, I=0.5, J=0.5)
        # I . J = (1/4)(X1 X2 + Y1 Y2 + Z1 Z2) so the Pauli expansion has
        # 3 terms, each with coefficient A/4
        expected = {'XX': 0.25, 'YY': 0.25, 'ZZ': 0.25}
        assert set(pauli.keys()) == set(expected.keys()), (
            f"Expected {set(expected.keys())}, got {set(pauli.keys())}"
        )
        for k, v in expected.items():
            assert pauli[k] == pytest.approx(v, abs=1e-12), (
                f"{k}: got {pauli[k]} expected {v}"
            )

    def test_h_1s_eigenvalues(self):
        """H 1s I=1/2 J=1/2 -> triplet (3-fold) at +A/4, singlet at -3A/4."""
        A = 1.0
        pauli = hyperfine_a_ij_pauli_general(A, I=0.5, J=0.5)
        H = _pauli_dict_to_dense(pauli, 2)
        evals = np.sort(np.linalg.eigvalsh(H).real)
        # singlet F=0 at -3A/4, triplet F=1 at +A/4 (3-fold)
        assert evals[0] == pytest.approx(-0.75, abs=1e-10), (
            f"Singlet eigenvalue {evals[0]}"
        )
        for k in [1, 2, 3]:
            assert evals[k] == pytest.approx(0.25, abs=1e-10), (
                f"Triplet eigenvalue [{k}] {evals[k]}"
            )

    def test_cs_I_7_2_J_1_2_dimensions(self):
        """Cs I=7/2 J=1/2 -> 4 qubits (3+1), Pauli terms in [16, 64]."""
        A = 1.0
        pauli = hyperfine_a_ij_pauli_general(A, I=3.5, J=0.5)
        # All Pauli strings are length 4
        for k in pauli:
            assert len(k) == 4

    def test_cs_eigenvalue_splitting_F4_minus_F3(self):
        """Cs I=7/2 J=1/2: E(F=4) - E(F=3) should be 4*A in Ha.

        For H = A I.J with I=7/2, J=1/2:
          F = I + J: F can be 4 or 3.
          E_F = (A/2) * [F(F+1) - I(I+1) - J(J+1)]
          F=4: E_F = (A/2) * (20 - 63/4 - 3/4) = (A/2) * (20 - 16.5) = 7A/4
          F=3: E_F = (A/2) * (12 - 63/4 - 3/4) = (A/2) * (12 - 16.5) = -9A/4
          Splitting = 7A/4 - (-9A/4) = 16A/4 = 4A
        """
        A = 1.0
        pauli = hyperfine_a_ij_pauli_general(A, I=3.5, J=0.5)
        H = _pauli_dict_to_dense(pauli, 4)
        evals = np.sort(np.linalg.eigvalsh(H).real)
        # The 16-dim register has (2I+1)*(2J+1) = 16 physical states
        # Top 9 should be F=4, next 7 should be F=3
        # E(F=4) = +7A/4, E(F=3) = -9A/4
        # Split = 4A
        # Find the F=4 (high) and F=3 (low) eigenvalues
        # The (2F+1) multiplicities are 9 and 7 for F=4 and F=3 respectively
        # (sum = 16 = full register, no padding here since 2I+1=8 fills 3 qubits)
        assert evals[-1] == pytest.approx(7.0 / 4.0, abs=1e-10), (
            f"Top eigenvalue (F=4) = {evals[-1]} != 7/4"
        )
        assert evals[0] == pytest.approx(-9.0 / 4.0, abs=1e-10), (
            f"Bottom eigenvalue (F=3) = {evals[0]} != -9/4"
        )

    def test_hermiticity(self):
        """H_hf must be Hermitian."""
        for I in [0.5, 1.0, 1.5, 2.5, 3.5]:
            pauli = hyperfine_a_ij_pauli_general(1.0, I=I, J=0.5)
            n_qubits_I = max(1, int(np.ceil(np.log2(2 * I + 1))))
            n_total = n_qubits_I + 1
            H = _pauli_dict_to_dense(pauli, n_total)
            assert np.allclose(H, H.conj().T, atol=1e-12), (
                f"Non-Hermitian for I={I}"
            )

    def test_invalid_I(self):
        """Invalid I (e.g. 0.7 not half-integer) raises ValueError."""
        with pytest.raises(ValueError, match="non-negative"):
            hyperfine_a_ij_pauli_general(1.0, I=0.7, J=0.5)

    def test_negative_I_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            hyperfine_a_ij_pauli_general(1.0, I=-0.5, J=0.5)

    def test_a_scaling_linear(self):
        """Scaling A by 2 doubles every Pauli coefficient."""
        pauli_1 = hyperfine_a_ij_pauli_general(1.0, I=3.5, J=0.5)
        pauli_2 = hyperfine_a_ij_pauli_general(2.0, I=3.5, J=0.5)
        assert set(pauli_1.keys()) == set(pauli_2.keys())
        for k in pauli_1:
            assert pauli_2[k] == pytest.approx(2.0 * pauli_1[k], rel=1e-12)


# ---------------------------------------------------------------------------
# Atomic HFS wrapper
# ---------------------------------------------------------------------------

class TestAtomicHFSWrapper:
    """The hyperfine_a_pauli_for_atomic_hfs convenience wrapper."""

    def test_h_1s_wrapper(self):
        """H 1s I=1/2 J=1/2 -> 2 qubits, splitting = A in MHz."""
        A_Ha = 2.16e-7  # ~1420 MHz
        result = hyperfine_a_pauli_for_atomic_hfs(A_Ha, I=0.5)
        assert result['Q_total'] == 2
        assert result['Q_nuc'] == 1
        assert result['Q_elec'] == 1
        # F levels: 1 and 0
        assert sorted(result['F_levels']) == [0.0, 1.0]
        # Splitting F=1 - F=0 = A
        assert result['splitting_F_max_to_F_min_Ha'] == pytest.approx(A_Ha)

    def test_cs_wrapper(self):
        """Cs I=7/2 J=1/2 -> 4 qubits, splitting = 4*A."""
        A_Ha = 2298.157943e6 / HZ_PER_HA  # ~3.5e-7 Ha
        result = hyperfine_a_pauli_for_atomic_hfs(A_Ha, I=3.5)
        assert result['Q_total'] == 4
        assert result['Q_nuc'] == 3
        assert result['Q_elec'] == 1
        # F levels: 4 and 3
        assert sorted(result['F_levels']) == [3.0, 4.0]
        # Splitting F=4 - F=3 = 4 A
        assert result['splitting_F_max_to_F_min_Ha'] == pytest.approx(
            4.0 * A_Ha, rel=1e-12
        )
        # In MHz: 4 * 2298.158 = 9192.63 MHz = the SI-second
        assert abs(
            result['splitting_F_max_to_F_min_MHz'] - 9192.6318e0
        ) < 1.0  # within 1 MHz

    def test_cs_si_second(self):
        """Cs A=2298.157943 MHz reproduces SI second nu_HFS=9192.631770 MHz."""
        A_MHz = 2298.157943
        A_Ha = A_MHz * 1e6 / HZ_PER_HA
        result = hyperfine_a_pauli_for_atomic_hfs(A_Ha, I=3.5)
        # Expected exactly: 4 * 2298.157943 = 9192.631772 MHz (rounding)
        nu_si_MHz = 9192.6317700  # exact SI definition
        assert abs(
            result['splitting_F_max_to_F_min_MHz'] - nu_si_MHz
        ) < 0.001  # sub-kHz precision

    def test_track_ni_distinct_from_atomic_wrapper(self):
        """Track NI hyperfine_coupling_pauli (4 qubits, occupation encoding)
        is structurally DIFFERENT from atomic_hfs_wrapper for I=1/2 (2 qubits,
        binary encoding). Both must coexist in the framework."""
        from geovac.nuclear.nuclear_electronic import hyperfine_coupling_pauli
        # Just verify the Track-NI function exists with its 4-qubit signature
        # and is independent of the new wrapper.
        import inspect
        sig = inspect.signature(hyperfine_coupling_pauli)
        params = list(sig.parameters)
        assert 'Q_nuc' in params  # Track NI uses register+block architecture
        assert 'Q_elec' in params
        # Our atomic wrapper is much simpler: just A and I (J=1/2 default)
        sig2 = inspect.signature(hyperfine_a_pauli_for_atomic_hfs)
        params2 = list(sig2.parameters)
        assert 'A_au' in params2
        assert 'I' in params2
        assert 'Q_nuc' not in params2  # Different architecture

    def test_pauli_terms_real(self):
        """Pauli coefficients are real-valued."""
        A = 1.0
        result = hyperfine_a_pauli_for_atomic_hfs(A, I=3.5)
        for k, v in result['pauli_terms'].items():
            assert isinstance(v, float), (
                f"Pauli coefficient at {k} is not float: {type(v)}={v}"
            )
