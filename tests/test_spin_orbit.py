"""Tests for geovac.spin_orbit (Track T2 deliverable).

Verifies the closed-form Breit-Pauli spin-orbit matrix element in
the (n, κ, m_j) basis, the Kramers l=0 cancellation, and the
exact Z⁴ scaling.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Integer, Rational, Symbol, simplify, expand

from geovac.dirac_matrix_elements import (
    DiracLabel,
    alpha_sym,
    Z_sym,
    iter_dirac_labels,
)
from geovac.spin_orbit import (
    so_diagonal_matrix_element,
    build_so_hamiltonian_block,
    verify_z4_scaling,
)


# ---------------------------------------------------------------------------
# Kramers: l = 0 cancellation
# ---------------------------------------------------------------------------


class TestKramersL0:
    """For l=0 (κ=−1) the SO matrix element must be exactly 0 regardless
    of n, Z, α — the factor (κ+1) in ⟨L·S⟩ vanishes before the divergent
    ⟨1/r³⟩ is ever touched.
    """

    def test_1s_hydrogen(self):
        # n=1, κ=−1 (l=0, j=1/2). Only allowed κ for n=1.
        val = so_diagonal_matrix_element(1, -1, Z=1)
        assert val == 0

    def test_2s_hydrogen(self):
        val = so_diagonal_matrix_element(2, -1, Z=1)
        assert val == 0

    def test_2s_symbolic_Z(self):
        val = so_diagonal_matrix_element(2, -1)  # symbolic Z
        assert val == 0

    def test_3s_heavy(self):
        # Even at Z=38 (Sr): SO contribution for s-states is zero.
        val = so_diagonal_matrix_element(3, -1, Z=38)
        assert val == 0

    def test_l0_never_evaluates_divergent_r3(self):
        # The routine must short-circuit BEFORE calling
        # inverse_r_cubed_hydrogenic(n, 0). Relying on ValueError from
        # that function being NOT raised is the test.
        # (If it were called, we'd get a ValueError — we don't.)
        val = so_diagonal_matrix_element(5, -1, Z=Integer(1))
        assert val == 0


# ---------------------------------------------------------------------------
# Hydrogen 2p matrix elements (benchmark)
# ---------------------------------------------------------------------------


class TestHydrogen2pBreitPauli:
    """The matrix elements ⟨2p_j | H_SO | 2p_j⟩ for hydrogen are:

        H_SO(2p_{3/2}) = (Zα²/2) · ⟨L·S⟩ · ⟨1/r³⟩
                       = (α²/2) · (1/2)  · (1/24)   =  α² / 96
        H_SO(2p_{1/2}) = (α²/2) · (−1)   · (1/24)   = −α² / 48

    These are the Breit-Pauli leading-order (non-relativistic-limit)
    matrix elements. The combined fine-structure splitting including
    Darwin + relativistic KE terms is the familiar α²/16 per Hartree
    (at fixed n), but *the spin-orbit contribution alone* is α²/32 —
    derived below as the explicit test.
    """

    def test_2p_three_halves(self):
        # κ = −2  (j = l + 1/2 = 3/2, l = 1)
        val = so_diagonal_matrix_element(2, -2, Z=1)
        expected = alpha_sym**2 / Integer(96)
        assert simplify(val - expected) == 0

    def test_2p_one_half(self):
        # κ = +1  (j = l − 1/2 = 1/2, l = 1)
        val = so_diagonal_matrix_element(2, +1, Z=1)
        expected = -alpha_sym**2 / Integer(48)
        assert simplify(val - expected) == 0

    def test_2p_SO_splitting(self):
        """The SO-only piece of the 2p fine-structure splitting:

            Δ_SO = H_SO(2p_{3/2}) − H_SO(2p_{1/2})
                 =  α²/96 − (−α²/48)
                 =  α²/96 + α²/48
                 =  α²/96 + 2α²/96
                 =  3α²/96  =  α²/32 .

        This is the *spin-orbit* contribution. The full hydrogenic fine-
        structure splitting (Dirac-equation) for the 2p pair is
        α²·m_e·c²·[1/4 − 1/2]·(1/2 − 1) / (different normalization), and
        in atomic-unit Ha is commonly quoted as α⁴/32 (i.e. α²/32 in the
        α²-reduced form when Z⁴ is absorbed). The intermediate H_SO only
        captures the Breit-Pauli spin-orbit term, not the Darwin nor
        mass-velocity corrections. α²/32 is the number the closed-form
        ⟨H_SO⟩ delivers and is what the test checks.
        """
        top = so_diagonal_matrix_element(2, -2, Z=1)
        bot = so_diagonal_matrix_element(2, +1, Z=1)
        splitting = simplify(top - bot)
        assert simplify(splitting - alpha_sym**2 / Integer(32)) == 0


# ---------------------------------------------------------------------------
# Z⁴ scaling
# ---------------------------------------------------------------------------


class TestZ4Scaling:
    """For any (n, κ) with l ≥ 1 the SO matrix element is exactly
    proportional to Z⁴ in the Coulomb potential. Verified symbolically
    across Tier-2 target atoms.
    """

    def test_symbolic_Z4_structure(self):
        # With symbolic Z, the expression must factor as Z⁴ × (stuff in α, n, κ).
        val = so_diagonal_matrix_element(2, -2, Z=Z_sym)
        # Divide by Z⁴: remainder must be Z-free.
        ratio = simplify(val / Z_sym**4)
        assert Z_sym not in ratio.free_symbols

    def test_Z4_ratios_2p32(self):
        """H_SO(2p_{3/2}, Z) / H_SO(2p_{3/2}, Z_ref) = (Z/Z_ref)⁴ exactly."""
        results = verify_z4_scaling(2, -2, Z_values=[1, 3, 4, 38])
        ref = results[1]
        for Z, val in results.items():
            expected_ratio = Integer(Z) ** 4
            ratio = simplify(val / ref)
            assert simplify(ratio - expected_ratio) == 0, (
                f"Z={Z}: ratio {ratio} != {expected_ratio}")

    def test_Z4_ratios_2p12(self):
        results = verify_z4_scaling(2, +1, Z_values=[1, 3, 4, 38])
        ref = results[1]
        for Z, val in results.items():
            expected_ratio = Integer(Z) ** 4
            ratio = simplify(val / ref)
            assert simplify(ratio - expected_ratio) == 0

    def test_Sr_over_hydrogen_factor(self):
        """SrH implication: at Z=38, SO is 38⁴ = 2,085,136 times hydrogen.
        This is why heavy atoms need relativistic treatment.
        """
        results = verify_z4_scaling(2, -2, Z_values=[1, 38])
        ratio = simplify(results[38] / results[1])
        assert ratio == Integer(38) ** 4
        assert int(ratio) == 2085136


# ---------------------------------------------------------------------------
# Non-relativistic limit α → 0
# ---------------------------------------------------------------------------


class TestNonRelLimit:
    """In the limit α → 0 every SO matrix element must vanish (the SO
    operator is explicitly O(α²)).
    """

    def test_2p32_hydrogen_alpha_zero(self):
        val = so_diagonal_matrix_element(2, -2, Z=1)
        assert val.subs(alpha_sym, 0) == 0

    def test_2p12_hydrogen_alpha_zero(self):
        val = so_diagonal_matrix_element(2, +1, Z=1)
        assert val.subs(alpha_sym, 0) == 0

    def test_heavy_alpha_zero(self):
        val = so_diagonal_matrix_element(3, -2, Z=38)
        assert val.subs(alpha_sym, 0) == 0


# ---------------------------------------------------------------------------
# Block assembly
# ---------------------------------------------------------------------------


class TestBlockAssembly:
    """build_so_hamiltonian_block(n_max, Z) must produce one entry per
    DiracLabel in iter_dirac_labels(n_max), with l=0 entries = 0 and
    all l≥1 entries equal to the closed-form so_diagonal_matrix_element.
    """

    def test_n_max_1(self):
        # Only n=1, l=0, κ=−1; two m_j values (±1/2) → 2 labels.
        block = build_so_hamiltonian_block(1, Z=1)
        labels = list(block.diag.keys())
        assert len(labels) == 2
        for lab in labels:
            assert lab.n_fock == 1 and lab.kappa == -1
            assert block.diag[lab] == 0

    def test_n_max_2_cardinality(self):
        # n=1: 2 labels (1s, m=±1/2)
        # n=2: 2s (κ=−1, 2 m_j) + 2p_{3/2} (κ=−2, 4 m_j) + 2p_{1/2} (κ=+1, 2 m_j)
        #    = 2 + 4 + 2 = 8 labels
        # Total: 2 + 8 = 10 labels
        block = build_so_hamiltonian_block(2, Z=1)
        assert len(block) == 10

    def test_m_j_independence(self):
        """H_SO is m_j-independent: all four m_j values of 2p_{3/2} must
        have the same eigenvalue.
        """
        block = build_so_hamiltonian_block(2, Z=1)
        p32_vals = [
            block.diag[DiracLabel(n_fock=2, kappa=-2, two_m_j=mj)]
            for mj in (-3, -1, +1, +3)
        ]
        for v in p32_vals:
            assert simplify(v - p32_vals[0]) == 0
        assert simplify(p32_vals[0] - alpha_sym**2 / Integer(96)) == 0

    def test_l0_all_zero(self):
        """Every l=0 label (κ=−1, any n, any m_j) must have 0 eigenvalue."""
        block = build_so_hamiltonian_block(3, Z=38)
        for lab, val in block.diag.items():
            if lab.l == 0:
                assert val == 0

    def test_consistency_with_direct_call(self):
        """Block values must match direct so_diagonal_matrix_element calls."""
        block = build_so_hamiltonian_block(2, Z=Integer(3))
        for lab, val in block.diag.items():
            direct = so_diagonal_matrix_element(lab.n_fock, lab.kappa, Z=Integer(3))
            assert simplify(val - direct) == 0


# ---------------------------------------------------------------------------
# Regression: imports and T1/D1 infrastructure still healthy
# ---------------------------------------------------------------------------


class TestRegressionTrack1():
    """Sanity: T2 imports of T1 API must still resolve and produce the
    T1-advertised closed forms. If anyone edits T1 and breaks these,
    T2 CI catches it.
    """

    def test_T1_LS_identity(self):
        from geovac.dirac_matrix_elements import angular_matrix_L_dot_S
        # κ=−2: expected (−(−2)−1)/2 = 1/2
        assert angular_matrix_L_dot_S(-2, 1, -2, 1) == Rational(1, 2)
        # κ=+1: expected (−1−1)/2 = −1
        assert angular_matrix_L_dot_S(+1, 1, +1, 1) == Rational(-1, 1)

    def test_T1_r3_identity(self):
        from geovac.dirac_matrix_elements import inverse_r_cubed_hydrogenic
        # ⟨2p| 1/r³ |2p⟩ for Z=1: 1 / (8 · 1 · 3/2 · 2) = 1/24
        val = inverse_r_cubed_hydrogenic(2, 1, Z=Integer(1))
        assert simplify(val - Rational(1, 24)) == 0
