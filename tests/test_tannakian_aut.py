"""Tests for the TC-1e Aut^⊗(ω) inclusion of U^*_Levi.

See: Sprint Q5'-Tannakian-Closure TC-1e
(debug/sprint_q5p_tc1e_aut_inclusion_memo.md).
"""

import pytest
from sympy import Matrix, Rational, eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    levi_unipotent_action,
    trivial_rep,
    unit_object,
    verify_natural_auto_group_law,
    verify_natural_auto_invertibility,
    verify_natural_auto_naturality,
    verify_natural_auto_tensor,
    verify_natural_auto_unit,
)


def _J2(n_max):
    return FinDimRep(
        n_max=n_max, dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2",
    )


def _J3(n_max):
    return FinDimRep(
        n_max=n_max, dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3",
    )


def _K2(n_max):
    return FinDimRep(
        n_max=n_max, dim=2,
        endos={(2, 0, 1): Matrix([[0, 1], [0, 0]])},
        label="K2",
    )


# ---------------------------------------------------------------------
# levi_unipotent_action: closed-form sanity
# ---------------------------------------------------------------------


def test_levi_action_zero_param_returns_identity():
    J = _J2(2)
    eta = levi_unipotent_action({}, J)
    assert eta == eye(2)


def test_levi_action_jordan_one_param():
    J = _J2(2)
    t = {(1, 0, 0): Rational(1)}
    eta = levi_unipotent_action(t, J)
    # exp([[0,1],[0,0]]) = [[1,1],[0,1]]
    assert eta == Matrix([[1, 1], [0, 1]])


def test_levi_action_jordan_3d_rational_param():
    J = _J3(2)
    t = {(1, 0, 0): Rational(2, 3)}
    eta = levi_unipotent_action(t, J)
    # exp((2/3)*N_3) = I + (2/3) N + (2/3)^2/2 N^2 = I + (2/3) N + (2/9) N^2
    # N = [[0,1,0],[0,0,1],[0,0,0]], N^2 = [[0,0,1],[0,0,0],[0,0,0]]
    expected = Matrix([
        [1, Rational(2, 3), Rational(2, 9)],
        [0, 1, Rational(2, 3)],
        [0, 0, 1],
    ])
    assert eta == expected


def test_levi_action_trivial_rep_is_identity():
    T = trivial_rep(2, dim=2)
    t = {(1, 0, 0): Rational(5, 7)}
    eta = levi_unipotent_action(t, T)
    assert eta == eye(2)


def test_levi_action_unit_is_identity():
    one = unit_object(2)
    t = {(1, 0, 0): Rational(3)}
    eta = levi_unipotent_action(t, one)
    assert eta == eye(1)


# ---------------------------------------------------------------------
# Invertibility
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_invertibility_panel(n_max):
    panel_reps = [unit_object(n_max), trivial_rep(n_max, 2), _J2(n_max), _J3(n_max), _K2(n_max)]
    t = {(1, 0, 0): Rational(2, 5), (2, 0, 1): Rational(-3, 4)}
    for V in panel_reps:
        v = verify_natural_auto_invertibility(t, V)
        assert v["invertible"]
        assert v["bit_exact"]


# ---------------------------------------------------------------------
# Unit
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_unit_axiom(n_max):
    t = {(1, 0, 0): Rational(1), (2, 0, 1): Rational(-1)}
    v = verify_natural_auto_unit(t, n_max)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Naturality
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_naturality_mono(n_max):
    T1 = unit_object(n_max)
    J = _J2(n_max)
    f = RepMorphism(T1, J, Matrix([[1], [0]]), label="f")
    t = {(1, 0, 0): Rational(2, 3)}
    v = verify_natural_auto_naturality(t, f)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_naturality_epi(n_max):
    T1 = unit_object(n_max)
    J = _J2(n_max)
    f = RepMorphism(J, T1, Matrix([[0, 1]]), label="f")
    t = {(1, 0, 0): Rational(2, 3)}
    v = verify_natural_auto_naturality(t, f)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Tensor compatibility
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_tensor_compatibility_J2_J2(n_max):
    J = _J2(n_max)
    t = {(1, 0, 0): Rational(1, 3)}
    v = verify_natural_auto_tensor(t, J, J)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_tensor_compatibility_J2_K2(n_max):
    J = _J2(n_max)
    K = _K2(n_max)
    t = {(1, 0, 0): Rational(1, 2), (2, 0, 1): Rational(-1, 5)}
    v = verify_natural_auto_tensor(t, J, K)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_tensor_compatibility_unit_jordan(n_max):
    T1 = unit_object(n_max)
    J = _J3(n_max)
    t = {(1, 0, 0): Rational(2, 7)}
    v = verify_natural_auto_tensor(t, T1, J)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Group law
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_group_law_distinct_generators(n_max):
    J = _J2(n_max)
    K = _K2(n_max)
    t1 = {(1, 0, 0): Rational(1)}
    t2 = {(2, 0, 1): Rational(1)}
    v = verify_natural_auto_group_law(t1, t2, J)
    assert v["bit_exact"]
    v2 = verify_natural_auto_group_law(t1, t2, K)
    assert v2["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_group_law_doubling_same_param(n_max):
    J = _J3(n_max)
    t = {(1, 0, 0): Rational(3, 5)}
    v = verify_natural_auto_group_law(t, t, J)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_group_law_with_zero_is_identity(n_max):
    J = _J2(n_max)
    t1 = {(1, 0, 0): Rational(2, 3)}
    t2 = {}
    v = verify_natural_auto_group_law(t1, t2, J)
    assert v["bit_exact"]
