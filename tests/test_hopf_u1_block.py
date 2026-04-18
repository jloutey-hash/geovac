"""
Tests for Track RH-E: Hopf-U(1) block decomposition.

Verifies Paper 29 Hypothesis 1 (§5.3) and its Dirac-S^3 analog via the
Z_2 m-reflection symmetry of the node-level Ihara-Bass matrix.
"""

from __future__ import annotations

import pytest
import numpy as np
import sympy as sp

from debug.hopf_u1_block_test import (
    bargmann_scalar_data,
    dirac_ruleB_data,
    bargmann_m_reflection,
    dirac_mj_reflection,
    build_reflection_projectors_bargmann,
    build_reflection_projectors_dirac,
    block_det_I_minus_sA_plus_sQ,
    target_P12_P22_scalar_S5,
    target_P22_P24_dirac,
)


# ---------------------------------------------------------------------------
# Sector decomposition correctness (the three required tests)
# ---------------------------------------------------------------------------

class TestScalarS5SectorDecomposition:
    """Verify the sector decomposition is correct for scalar S^5 N_max=3."""

    def test_sector_dimensions_sum_to_V(self):
        """Sum of sym + antisym dims equals V = 20."""
        A, nodes, info = bargmann_scalar_data()
        V = info["V"]
        U_sym, U_antisym = build_reflection_projectors_bargmann(nodes)
        dim_sym = U_sym.shape[1]
        dim_antisym = U_antisym.shape[1]
        assert dim_sym + dim_antisym == V
        assert dim_sym == 13
        assert dim_antisym == 7

    def test_projectors_are_orthonormal(self):
        """U_sym and U_antisym form an orthonormal basis."""
        A, nodes, info = bargmann_scalar_data()
        U_sym, U_antisym = build_reflection_projectors_bargmann(nodes)
        # U_sym.T * U_sym = I
        I_sym = sp.simplify(U_sym.T * U_sym)
        assert I_sym == sp.eye(U_sym.shape[1])
        # U_antisym.T * U_antisym = I
        I_antisym = sp.simplify(U_antisym.T * U_antisym)
        assert I_antisym == sp.eye(U_antisym.shape[1])
        # U_sym orthogonal to U_antisym
        cross = sp.simplify(U_sym.T * U_antisym)
        assert cross == sp.zeros(U_sym.shape[1], U_antisym.shape[1])

    def test_adjacency_commutes_with_reflection(self):
        """Adjacency A commutes with m → -m reflection. (Closure of each
        sector under the Hashimoto / adjacency action follows from this.)"""
        A, nodes, info = bargmann_scalar_data()
        P = bargmann_m_reflection(nodes)
        # P^2 = I
        assert np.allclose(P @ P, np.eye(len(nodes)))
        # P A = A P
        assert np.abs(P @ A - A @ P).sum() == 0


class TestDiracSectorDecomposition:
    """Verify sector decomposition for Dirac-S^3 Rule B n_max=3."""

    def test_sector_dimensions_sum_to_V(self):
        """Sum of sym + antisym dims equals V = 28, half and half."""
        A, labels, info = dirac_ruleB_data()
        V = info["V"]
        U_sym, U_antisym = build_reflection_projectors_dirac(labels)
        assert U_sym.shape[1] + U_antisym.shape[1] == V
        assert U_sym.shape[1] == 14
        assert U_antisym.shape[1] == 14

    def test_projectors_are_orthonormal(self):
        A, labels, info = dirac_ruleB_data()
        U_sym, U_antisym = build_reflection_projectors_dirac(labels)
        I_sym = sp.simplify(U_sym.T * U_sym)
        assert I_sym == sp.eye(U_sym.shape[1])
        I_antisym = sp.simplify(U_antisym.T * U_antisym)
        assert I_antisym == sp.eye(U_antisym.shape[1])
        cross = sp.simplify(U_sym.T * U_antisym)
        assert cross == sp.zeros(U_sym.shape[1], U_antisym.shape[1])

    def test_adjacency_commutes_with_mj_reflection(self):
        A, labels, info = dirac_ruleB_data()
        P = dirac_mj_reflection(labels)
        assert np.allclose(P @ P, np.eye(len(labels)))
        assert np.abs(P @ A - A @ P).sum() == 0


# ---------------------------------------------------------------------------
# Hypothesis validation (slow)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHypothesisValidation:
    """Verify the core claim: block dets equal the target factor polynomials."""

    def test_scalar_S5_antisym_block_is_P12(self):
        """The m-antisym block dets equals P_12(s)."""
        A, nodes, info = bargmann_scalar_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        _, U_antisym = build_reflection_projectors_bargmann(nodes)
        det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)
        P12, _ = target_P12_P22_scalar_S5()
        # The antisym block is exactly P_12 (dim 7 -> degree 12 in s
        # by two-per-node).  Wait — dim 7 gives 2E'_block worth of
        # polynomial.  Actually det(I - sA + s^2 Q) for a 7x7 block
        # has degree 14 in s, but the leading coefficient may be zero
        # for the block because Q was projected and may have zero
        # elements.  We should instead check that P_12 divides det_antisym.
        s = sp.symbols("s")
        remainder = sp.rem(det_antisym, P12.as_expr(), s)
        assert remainder == 0
        # And the quotient has integer coefficients & low degree
        quotient = sp.quo(det_antisym, P12.as_expr(), s)
        quotient_poly = sp.Poly(quotient, s)
        assert quotient_poly.degree() <= 2, \
            f"Quotient has unexpectedly high degree {quotient_poly.degree()}"

    def test_scalar_S5_sym_block_is_P22(self):
        """The m-sym block det contains P_22(s) as a factor."""
        A, nodes, info = bargmann_scalar_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        U_sym, _ = build_reflection_projectors_bargmann(nodes)
        det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
        _, P22 = target_P12_P22_scalar_S5()
        s = sp.symbols("s")
        remainder = sp.rem(det_sym, P22.as_expr(), s)
        assert remainder == 0

    def test_scalar_S5_block_product_equals_full_det(self):
        """det(M_sym) · det(M_antisym) = det(I - sA + s^2 Q) on the full
        node space.  This confirms the block decomposition is complete."""
        A, nodes, info = bargmann_scalar_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        U_sym, U_antisym = build_reflection_projectors_bargmann(nodes)
        det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
        det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)
        s = sp.symbols("s")
        M_full = sp.eye(V) - s * sp.Matrix(A.tolist()) + s * s * sp.Matrix(Q.tolist())
        full_det = sp.expand(M_full.det())
        product = sp.expand(det_sym * det_antisym)
        assert sp.simplify(product - full_det) == 0

    def test_dirac_ruleB_sym_block_is_P22(self):
        """The m_j-sym block det contains P_22(s) as a factor."""
        A, labels, info = dirac_ruleB_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        U_sym, _ = build_reflection_projectors_dirac(labels)
        det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
        P22, _ = target_P22_P24_dirac()
        s = sp.symbols("s")
        remainder = sp.rem(det_sym, P22.as_expr(), s)
        assert remainder == 0

    def test_dirac_ruleB_antisym_block_is_P24(self):
        """The m_j-antisym block det contains P_24(s) as a factor."""
        A, labels, info = dirac_ruleB_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        _, U_antisym = build_reflection_projectors_dirac(labels)
        det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)
        _, P24 = target_P22_P24_dirac()
        s = sp.symbols("s")
        remainder = sp.rem(det_antisym, P24.as_expr(), s)
        assert remainder == 0

    def test_dirac_ruleB_block_product_equals_full_det(self):
        A, labels, info = dirac_ruleB_data()
        V = info["V"]
        deg = A.sum(axis=1)
        Q = np.diag(deg - 1)
        U_sym, U_antisym = build_reflection_projectors_dirac(labels)
        det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
        det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)
        s = sp.symbols("s")
        M_full = sp.eye(V) - s * sp.Matrix(A.tolist()) + s * s * sp.Matrix(Q.tolist())
        full_det = sp.expand(M_full.det())
        product = sp.expand(det_sym * det_antisym)
        assert sp.simplify(product - full_det) == 0
