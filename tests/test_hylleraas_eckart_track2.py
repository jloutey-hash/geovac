"""Track 2 tests for Hylleraas-Eckart double-alpha integration into hylleraas_r12.py.

Verifies:
1. beta=0 bit-identical regression: mode='eckart_double_alpha' with beta=0
   reproduces mode='single_alpha' for S, V_ne, V_ee and (with matching
   quadrature defaults) for kinetic T, hence for total energy.
2. Hermiticity of S and H at beta > 0.
3. Variational bound: E_var >= E_exact for He 1^1S at all tested beta.
4. Eckart matrix elements: per-pair B->0 reduces to single-alpha.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.hylleraas_r12 import (
    HylleraasBasisFn,
    assemble_matrices,
    hylleraas_basis_3p,
    hylleraas_basis_6p,
    hylleraas_basis_total_degree,
    overlap_element,
    overlap_element_eckart,
    potential_vne_element,
    potential_vne_element_eckart,
    potential_vee_element,
    potential_vee_element_eckart,
    solve_hylleraas_state,
)


# ---------------------------------------------------------------------------
# Per-element beta=0 reduction to single-alpha
# ---------------------------------------------------------------------------

class TestEckartElementsReduceAtBetaZero:
    BASIS = hylleraas_basis_total_degree(3)
    ALPHA = 1.6875
    Z = 2.0

    def test_overlap_per_pair(self):
        worst = 0.0
        for bf_p in self.BASIS:
            for bf_q in self.BASIS:
                S_single = overlap_element(bf_p, bf_q, self.ALPHA)
                S_eckart = overlap_element_eckart(
                    bf_p, bf_q, self.ALPHA, beta=0.0, spin='singlet'
                )
                rel = abs(S_eckart - S_single) / max(abs(S_single), 1e-30)
                worst = max(worst, rel)
        assert worst < 1e-12, f"Worst overlap rel err: {worst:.2e}"

    def test_vne_per_pair(self):
        worst = 0.0
        for bf_p in self.BASIS:
            for bf_q in self.BASIS:
                V_single = potential_vne_element(bf_p, bf_q, self.ALPHA, self.Z)
                V_eckart = potential_vne_element_eckart(
                    bf_p, bf_q, self.ALPHA, beta=0.0, Z=self.Z, spin='singlet'
                )
                rel = abs(V_eckart - V_single) / max(abs(V_single), 1e-30)
                worst = max(worst, rel)
        assert worst < 1e-12, f"Worst V_ne rel err: {worst:.2e}"

    def test_vee_per_pair(self):
        worst = 0.0
        for bf_p in self.BASIS:
            for bf_q in self.BASIS:
                V_single = potential_vee_element(bf_p, bf_q, self.ALPHA)
                V_eckart = potential_vee_element_eckart(
                    bf_p, bf_q, self.ALPHA, beta=0.0, spin='singlet'
                )
                rel = abs(V_eckart - V_single) / max(abs(V_single), 1e-30)
                worst = max(worst, rel)
        assert worst < 1e-12, f"Worst V_ee rel err: {worst:.2e}"


# ---------------------------------------------------------------------------
# Full-matrix beta=0 regression
# ---------------------------------------------------------------------------

class TestFullMatrixBetaZeroRegression:
    """At beta=0, mode='eckart_double_alpha' must reproduce mode='single_alpha'
    bit-identically across H and S matrices for He 1^1S.
    """

    @pytest.mark.parametrize("basis_size", ["3p", "6p", "omega_3"])
    def test_matrix_bit_identical(self, basis_size):
        if basis_size == "3p":
            basis = hylleraas_basis_3p()
        elif basis_size == "6p":
            basis = hylleraas_basis_6p()
        else:
            basis = hylleraas_basis_total_degree(int(basis_size.split("_")[1]))

        alpha = 1.6875
        Z = 2.0
        H1, S1 = assemble_matrices(basis, alpha, Z, mode='single_alpha')
        H2, S2 = assemble_matrices(basis, alpha, Z,
                                    mode='eckart_double_alpha', beta=0.0)
        # S must be machine-precision bit-identical.
        assert np.max(np.abs(S1 - S2)) < 1e-13, (
            f"{basis_size}: max |S_single - S_eckart| = {np.max(np.abs(S1 - S2)):.2e}"
        )
        # H must be bit-identical (kinetic uses same quad defaults in both paths).
        assert np.max(np.abs(H1 - H2)) < 1e-13, (
            f"{basis_size}: max |H_single - H_eckart| = {np.max(np.abs(H1 - H2)):.2e}"
        )

    def test_energy_bit_identical_3p(self):
        """He 1^1S energy at beta=0 must be bit-identical between modes."""
        basis = hylleraas_basis_3p()
        alpha = 1.6875
        Z = 2.0
        res1 = solve_hylleraas_state(basis, alpha, Z, mode='single_alpha')
        res2 = solve_hylleraas_state(basis, alpha, Z,
                                      mode='eckart_double_alpha', beta=0.0)
        assert abs(res1.energy - res2.energy) < 1e-13, (
            f"E_single={res1.energy}, E_eckart_b0={res2.energy}"
        )


# ---------------------------------------------------------------------------
# Hermiticity at beta > 0
# ---------------------------------------------------------------------------

class TestHermiticityAtBetaNonzero:
    @pytest.mark.parametrize("beta", [0.1, 0.3, 0.5])
    def test_overlap_hermitian(self, beta):
        basis = hylleraas_basis_total_degree(3)
        alpha = 1.6875
        for bf_p in basis:
            for bf_q in basis:
                v_pq = overlap_element_eckart(bf_p, bf_q, alpha, beta)
                v_qp = overlap_element_eckart(bf_q, bf_p, alpha, beta)
                # These are functions of (L, M, N) summed in the same way,
                # so should be bit-identical.
                assert v_pq == v_qp

    @pytest.mark.parametrize("beta", [0.1, 0.3])
    def test_matrix_hermitian(self, beta):
        basis = hylleraas_basis_3p()
        H, S = assemble_matrices(
            basis, 1.6875, 2.0,
            mode='eckart_double_alpha', beta=beta,
        )
        # Symmetric to machine precision (matrix is explicitly symmetrized).
        assert np.max(np.abs(H - H.T)) < 1e-13
        assert np.max(np.abs(S - S.T)) < 1e-13


# ---------------------------------------------------------------------------
# Variational bound at beta > 0
# ---------------------------------------------------------------------------

class TestVariationalBoundEckart:
    """At beta > 0, the He 1^1S variational energy must remain above the
    exact non-relativistic value (-2.903724 Ha)."""

    EXACT_HE_1S = -2.903724
    ALPHA = 1.6875

    @pytest.mark.parametrize("beta", [0.0, 0.1, 0.3, 0.5, 1.0])
    def test_he_1s_variational_3p(self, beta):
        basis = hylleraas_basis_3p()
        res = solve_hylleraas_state(
            basis, self.ALPHA, 2.0,
            mode='eckart_double_alpha', beta=beta,
        )
        # Variational principle: E_var >= E_exact.
        assert res.energy >= self.EXACT_HE_1S, (
            f"beta={beta}: E={res.energy} below exact {self.EXACT_HE_1S}"
        )

    def test_he_1s_beta_zero_optimal_at_3p(self):
        """For He 1^1S in the 3p basis, beta=0 should give the lowest energy
        (since the ground state is symmetric in r_1, r_2 -- no asymmetry
        gain). Asymmetric trial only helps excited / antisymmetric states.
        """
        basis = hylleraas_basis_3p()
        E_at_zero = solve_hylleraas_state(
            basis, self.ALPHA, 2.0, mode='eckart_double_alpha', beta=0.0,
        ).energy
        for beta in [0.1, 0.3, 0.5]:
            E = solve_hylleraas_state(
                basis, self.ALPHA, 2.0, mode='eckart_double_alpha', beta=beta,
            ).energy
            assert E >= E_at_zero - 1e-10, (
                f"beta={beta} gives E={E} below E(beta=0)={E_at_zero}; "
                "ground state should not benefit from asymmetry."
            )


# ---------------------------------------------------------------------------
# Mode kwarg validation
# ---------------------------------------------------------------------------

class TestModeKwarg:
    def test_invalid_mode_rejected(self):
        with pytest.raises(ValueError):
            assemble_matrices(
                hylleraas_basis_3p(), 1.6875, 2.0, mode='nonsense'
            )

    def test_invalid_spin_rejected(self):
        with pytest.raises(ValueError):
            assemble_matrices(
                hylleraas_basis_3p(), 1.6875, 2.0,
                mode='eckart_double_alpha', beta=0.3, spin='quintet',
            )
