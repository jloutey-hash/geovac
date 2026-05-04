"""Tests for geovac.spinor_operator_system: Weyl-spinor lift on S^3.

Sprint WH1-R3.2 equation verification per CLAUDE.md Sec 13.4a.

Verifies:

  - Basis dimension formula: dim_H(n_max) = sum_{n=1..n_max} n*(n+1).
  - SpinorLabel validation (out-of-range two_m_j, l > n - 1, etc.).
  - Clebsch-Gordan coefficients sum-to-one identity (alpha_+^2 + alpha_-^2 = 1).
  - CG edge cases: m_j = +/-(l + 1/2) gives one alpha = 1, the other = 0.
  - Camporesi-Higuchi Dirac eigenvalues = n_fock + 1/2 exactly (sympy Rational).
  - Identity matrix is in O (constant-multiplier branch).
  - O is *-closed (each generator's conjugate transpose is in O).
  - Strict inequality dim(O) < N^2 (O is genuinely an operator system).
  - dim(O) at small n_max (regression baseline).
  - Performance smoke: build at n_max=3 < 30s.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import sympy as sp
from sympy import Rational, simplify, sqrt

from geovac.spinor_operator_system import (
    SpinorLabel,
    SpinorTruncatedOperatorSystem,
    build_spinor_multiplier_matrix,
    camporesi_higuchi_dirac_matrix,
    cg_coefficients,
    spinor_basis,
    spinor_dim,
    spinor_graph_distance,
    spinor_label_strings,
)


# ---------------------------------------------------------------------------
# Basis dimension and label validation
# ---------------------------------------------------------------------------


def test_spinor_dim_table():
    """dim_H(n_max) = sum_{n=1..n_max} n*(n+1)."""
    assert spinor_dim(1) == 2          # 1*2
    assert spinor_dim(2) == 8          # 1*2 + 2*3
    assert spinor_dim(3) == 20         # + 3*4
    assert spinor_dim(4) == 40         # + 4*5
    assert spinor_dim(5) == 70         # + 5*6


def test_spinor_basis_consistency():
    """basis size matches spinor_dim."""
    for n_max in (1, 2, 3, 4):
        b = spinor_basis(n_max)
        assert len(b) == spinor_dim(n_max)


def test_spinor_label_validation():
    """SpinorLabel rejects out-of-range quantum numbers."""
    SpinorLabel(n_fock=1, l=0, two_m_j=+1)   # OK
    SpinorLabel(n_fock=1, l=0, two_m_j=-1)   # OK
    SpinorLabel(n_fock=2, l=1, two_m_j=+3)   # OK (m_j = +3/2)
    with pytest.raises(ValueError):
        SpinorLabel(n_fock=0, l=0, two_m_j=+1)
    with pytest.raises(ValueError):
        SpinorLabel(n_fock=1, l=1, two_m_j=+1)  # l = n_fock violates l <= n-1
    with pytest.raises(ValueError):
        SpinorLabel(n_fock=1, l=0, two_m_j=+2)  # even (must be odd)
    with pytest.raises(ValueError):
        SpinorLabel(n_fock=2, l=1, two_m_j=+5)  # |two_m_j| > 2l + 1


def test_label_properties():
    """SpinorLabel m_j and j properties are sympy Rationals."""
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=+1)
    assert lab.m_j == Rational(1, 2)
    assert lab.j == Rational(3, 2)
    lab2 = SpinorLabel(n_fock=3, l=2, two_m_j=-3)
    assert lab2.m_j == Rational(-3, 2)
    assert lab2.j == Rational(5, 2)


# ---------------------------------------------------------------------------
# Clebsch-Gordan coefficients
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_cg_sum_to_one(n_max):
    """alpha_+^2 + alpha_-^2 = 1 EXACTLY (sympy) for every label."""
    for lab in spinor_basis(n_max):
        ap, am, _, _ = cg_coefficients(lab)
        s = ap ** 2 + am ** 2
        assert simplify(s) == 1, (
            f"CG sum failed for {lab}: alpha_+^2 + alpha_-^2 = {simplify(s)}"
        )


def test_cg_edge_cases_max_m_j():
    """m_j = +(l + 1/2) gives alpha_+ = 1, alpha_- = 0."""
    # n_fock=2, l=1, m_j = +3/2 (max)
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=+3)
    ap, am, _, _ = cg_coefficients(lab)
    assert ap == 1, f"alpha_+ = {ap}, expected 1"
    assert am == 0, f"alpha_- = {am}, expected 0"


def test_cg_edge_cases_min_m_j():
    """m_j = -(l + 1/2) gives alpha_+ = 0, alpha_- = 1."""
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=-3)
    ap, am, _, _ = cg_coefficients(lab)
    assert ap == 0, f"alpha_+ = {ap}, expected 0"
    assert am == 1, f"alpha_- = {am}, expected 1"


def test_cg_specific_values():
    """Spot check: (n=2, l=1, m_j=+1/2) gives alpha_+^2 = 2/3, alpha_-^2 = 1/3.

    From CG formula: alpha_+ = sqrt((l + m_j + 1/2)/(2l+1)) = sqrt((1 + 1/2 + 1/2)/3)
                            = sqrt(2/3).
                     alpha_- = sqrt((l - m_j + 1/2)/(2l+1)) = sqrt((1 - 1/2 + 1/2)/3)
                            = sqrt(1/3).
    """
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=+1)
    ap, am, _, _ = cg_coefficients(lab)
    assert simplify(ap ** 2 - Rational(2, 3)) == 0
    assert simplify(am ** 2 - Rational(1, 3)) == 0


def test_cg_m_l_indices():
    """Verify m_l_plus = m_j - 1/2 and m_l_minus = m_j + 1/2 are integers."""
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=+1)  # m_j = +1/2
    _, _, mlp, mlm = cg_coefficients(lab)
    assert mlp == 0, f"m_l_plus = {mlp}, expected 0"   # 1/2 - 1/2 = 0
    assert mlm == 1, f"m_l_minus = {mlm}, expected 1"  # 1/2 + 1/2 = 1


# ---------------------------------------------------------------------------
# Camporesi-Higuchi Dirac matrix
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_dirac_diagonal_eigenvalues(n_max):
    """Native CH Dirac is diagonal with eigenvalue n_fock + 1/2."""
    basis = spinor_basis(n_max)
    D = camporesi_higuchi_dirac_matrix(basis)
    # Off-diagonal entries are zero
    off = D - np.diag(np.diag(D))
    assert np.allclose(off, 0, atol=1e-15)
    # Diagonal entries match n_fock + 1/2
    for i, b in enumerate(basis):
        expected = float(b.n_fock) + 0.5
        assert abs(D[i, i].real - expected) < 1e-15
        assert abs(D[i, i].imag) < 1e-15


def test_dirac_eigenvalues_sympy_rational():
    """The CH eigenvalues are sympy Rationals (exact half-integers)."""
    from geovac.dirac_s3 import dirac_eigenvalue_abs
    for n_ch in range(5):
        lam = dirac_eigenvalue_abs(n_ch, convention="ch")
        assert isinstance(lam, Rational)
        assert lam == Rational(2 * n_ch + 3, 2)


# ---------------------------------------------------------------------------
# SpinorTruncatedOperatorSystem
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_identity_in_O(n_max):
    """Identity matrix is in spinor O (constant-multiplier branch)."""
    O = SpinorTruncatedOperatorSystem(n_max)
    in_O, residual = O.identity_in_O()
    assert in_O, f"identity not in spinor O at n_max={n_max}, residual={residual:.2e}"
    assert residual < 1e-10


@pytest.mark.parametrize("n_max", [2, 3])
def test_star_closed(n_max):
    """Spinor O is *-closed."""
    O = SpinorTruncatedOperatorSystem(n_max)
    star_ok, failures = O.is_star_closed()
    assert star_ok, (
        f"*-closure violated at n_max={n_max}, {len(failures)} failures: {failures}"
    )


@pytest.mark.parametrize("n_max", [2, 3])
def test_O_strictly_smaller_than_envelope(n_max):
    """dim(O) < dim_H^2 — O is a proper operator system, not the C*-envelope."""
    O = SpinorTruncatedOperatorSystem(n_max)
    assert O.dim < O.envelope_dim, (
        f"dim(O)={O.dim} = N^2={O.envelope_dim}; would be the full algebra"
    )
    assert O.dim > 1


def test_dim_O_table():
    """Regression baseline: dim(O) at small n_max."""
    expected = {
        1: 1,
        2: 14,
        3: 55,
    }
    for n_max, dim_expected in expected.items():
        O = SpinorTruncatedOperatorSystem(n_max)
        assert O.dim == dim_expected, (
            f"dim(spinor O_{n_max}) = {O.dim}, expected {dim_expected}"
        )


def test_dim_H_table():
    """dim_H matches spinor_dim formula."""
    for n_max, expected in [(1, 2), (2, 8), (3, 20), (4, 40)]:
        O = SpinorTruncatedOperatorSystem(n_max)
        assert O.dim_H == expected
        assert O.envelope_dim == expected ** 2


# ---------------------------------------------------------------------------
# Multiplier matrix structure
# ---------------------------------------------------------------------------


def test_constant_multiplier_is_proportional_to_identity():
    """M_{N=1, L=0, M=0} (constant function) lifts to a multiple of I.

    This is the key check that the CG-decomposition is correctly normalized
    so that scalar identities lift to spinor identities.
    """
    O = SpinorTruncatedOperatorSystem(2)
    # Find the (1, 0, 0) multiplier
    for label, M in zip(O.multiplier_labels, O.multiplier_matrices):
        if label == (1, 0, 0):
            # M should be c * I for some constant c (real, > 0)
            diag = np.diag(M)
            assert np.allclose(M - np.diag(diag), 0, atol=1e-12), (
                f"M_{{1,0,0}} has nonzero off-diagonal: max = {np.abs(M - np.diag(diag)).max()}"
            )
            assert np.allclose(diag, diag[0], atol=1e-12), (
                f"M_{{1,0,0}} diagonal not constant: {diag}"
            )
            return
    pytest.fail("M_{1,0,0} multiplier not found")


def test_multiplier_matrix_hermiticity_pairs():
    """For each (N, L, M) in O, there exists (N, L, -M) such that
    M^*_{NLM} is in span of {M_{NL, M'}}. (Result: O is *-closed.)
    """
    # This is implicit in test_star_closed; here we just spot-check that
    # the conjugate-transpose pattern works for one specific case.
    O = SpinorTruncatedOperatorSystem(2)
    # Find M_{2,1,+1}; its conjugate transpose should be in span of
    # {M_{2,1,-1}} (with a phase).
    for label, M in zip(O.multiplier_labels, O.multiplier_matrices):
        if label == (2, 1, 1):
            M_dag = M.conj().T
            in_O, residual = O.contains(M_dag)
            assert in_O, f"M_{{2,1,1}}^dag not in O: residual {residual:.2e}"
            return
    pytest.fail("M_{2,1,1} multiplier not found")


# ---------------------------------------------------------------------------
# Graph distance
# ---------------------------------------------------------------------------


def test_graph_distance_self_zero():
    lab = SpinorLabel(n_fock=2, l=1, two_m_j=+1)
    assert spinor_graph_distance(lab, lab) == 0


def test_graph_distance_symmetric():
    a = SpinorLabel(n_fock=1, l=0, two_m_j=+1)
    b = SpinorLabel(n_fock=2, l=1, two_m_j=-1)
    assert spinor_graph_distance(a, b) == spinor_graph_distance(b, a)


def test_graph_distance_specific():
    # |1, 0, +1/2> and |2, 0, +1/2>: differ in n only by 1
    a = SpinorLabel(n_fock=1, l=0, two_m_j=+1)
    b = SpinorLabel(n_fock=2, l=0, two_m_j=+1)
    assert spinor_graph_distance(a, b) == 1
    # |2, 1, +1/2> and |2, 1, -1/2>: differ in m_j by 1
    c = SpinorLabel(n_fock=2, l=1, two_m_j=+1)
    d = SpinorLabel(n_fock=2, l=1, two_m_j=-1)
    assert spinor_graph_distance(c, d) == 1
    # |2, 1, +3/2> and |2, 1, -3/2>: differ by 3
    e = SpinorLabel(n_fock=2, l=1, two_m_j=+3)
    f = SpinorLabel(n_fock=2, l=1, two_m_j=-3)
    assert spinor_graph_distance(e, f) == 3


# ---------------------------------------------------------------------------
# Performance smoke
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_build_n_max_3_under_30s():
    """SpinorTruncatedOperatorSystem(3) must build in < 30s."""
    t0 = time.time()
    O = SpinorTruncatedOperatorSystem(3)
    elapsed = time.time() - t0
    assert elapsed < 30, f"build took {elapsed:.1f}s, exceeds 30s budget"
    assert O.dim_H == 20
    assert O.dim == 55
