"""Tests for the TC-1f SL_2 / Peter-Weyl substrate and full-panel injectivity.

See debug/sprint_q5p_tc1f_sl2_inclusion_memo.md.
"""

import pytest
from sympy import Integer, Matrix, Rational, eye

from geovac.pro_system import primitive_generators
from geovac.tannakian import (
    FinDimRep,
    PWMorphism,
    PWRep,
    _pw_standard_rep,
    _pw_sym2_rep,
    _pw_trivial_rep,
    _sl2_sym2_action,
    primitive_generator_rep,
    sl2_standard_action,
    verify_ga_sl2_commute,
    verify_injectivity_at_generator,
    verify_sl2_group_homomorphism,
    verify_sl2_invertibility,
    verify_sl2_tensor,
)


def _g_identity():
    return Matrix([[1, 0], [0, 1]])


def _g_unipotent():
    return Matrix([[1, 1], [0, 1]])


def _g_torus():
    return Matrix([[2, 0], [0, Rational(1, 2)]])


def _g_generic():
    # det = 5*3 - 2*7 = 15 - 14 = 1
    return Matrix([[5, 2], [7, 3]])


def _g_weyl():
    return Matrix([[0, 1], [-1, 0]])


# ---------------------------------------------------------------------
# PWRep / PWMorphism plumbing
# ---------------------------------------------------------------------


def test_pwrep_standard_dim():
    V = _pw_standard_rep()
    assert V.dim == 2
    assert V.label == "V_fund"


def test_pwrep_sym2_dim():
    V = _pw_sym2_rep()
    assert V.dim == 3


def test_pwrep_trivial_dim():
    V = _pw_trivial_rep(dim=5)
    assert V.dim == 5
    for g in [_g_identity(), _g_unipotent(), _g_torus()]:
        assert V.rho(g) == eye(5)


def test_pwrep_standard_rho_is_identity_action():
    V = _pw_standard_rep()
    g = _g_generic()
    assert V.rho(g) == g


def test_pwrep_repr_contains_dim_and_label():
    V = _pw_standard_rep()
    s = repr(V)
    assert "PWRep" in s
    assert "2" in s
    assert "V_fund" in s


def test_pwmorphism_construction():
    V = _pw_standard_rep()
    W = _pw_standard_rep("W")
    F = PWMorphism(V, W, eye(2), label="id")
    assert F.source is V
    assert F.target is W
    assert F.matrix == eye(2)


# ---------------------------------------------------------------------
# Sym^2 action explicit formula
# ---------------------------------------------------------------------


def test_sym2_action_identity():
    g = _g_identity()
    M = _sl2_sym2_action(g)
    assert M == eye(3)


def test_sym2_action_unipotent_closed_form():
    g = _g_unipotent()  # [[1,1],[0,1]]
    # Sym^2 of g on basis {e1^2, e1 e2, e2^2}:
    # [[1, 2, 1], [0, 1, 1], [0, 0, 1]]
    expected = Matrix([
        [1, 2, 1],
        [0, 1, 1],
        [0, 0, 1],
    ])
    assert _sl2_sym2_action(g) == expected


def test_sym2_action_torus_diag():
    g = _g_torus()  # diag(2, 1/2)
    # Sym^2: diag(4, 1, 1/4)
    expected = Matrix([
        [4, 0, 0],
        [0, 1, 0],
        [0, 0, Rational(1, 4)],
    ])
    assert _sl2_sym2_action(g) == expected


# ---------------------------------------------------------------------
# SL_2 invertibility
# ---------------------------------------------------------------------


@pytest.mark.parametrize("g_label,g_fn", [
    ("identity", _g_identity),
    ("unipotent", _g_unipotent),
    ("torus", _g_torus),
    ("generic", _g_generic),
    ("weyl", _g_weyl),
])
@pytest.mark.parametrize("V_fn", [_pw_standard_rep, _pw_sym2_rep])
def test_sl2_invertibility(g_label, g_fn, V_fn):
    g = g_fn()
    V = V_fn()
    v = verify_sl2_invertibility(g, V)
    assert v["g_in_SL2"]
    assert v["rho_invertible"]
    assert v["bit_exact"]


def test_sl2_invertibility_trivial_rep():
    g = _g_generic()
    V = _pw_trivial_rep(dim=4)
    v = verify_sl2_invertibility(g, V)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Tensor compatibility
# ---------------------------------------------------------------------


@pytest.mark.parametrize("g_label,g_fn", [
    ("identity", _g_identity),
    ("unipotent", _g_unipotent),
    ("torus", _g_torus),
    ("generic", _g_generic),
    ("weyl", _g_weyl),
])
@pytest.mark.parametrize("V_fn,W_fn", [
    (_pw_standard_rep, _pw_standard_rep),
    (_pw_standard_rep, _pw_sym2_rep),
    (_pw_sym2_rep, _pw_sym2_rep),
])
def test_sl2_tensor_kronecker(g_label, g_fn, V_fn, W_fn):
    g = g_fn()
    V = V_fn()
    W = W_fn()
    v = verify_sl2_tensor(g, V, W)
    assert v["bit_exact"]


def test_sl2_tensor_with_trivial_is_kronecker_with_identity():
    g = _g_generic()
    V = _pw_standard_rep()
    triv = _pw_trivial_rep(dim=1)
    v = verify_sl2_tensor(g, V, triv)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Group homomorphism
# ---------------------------------------------------------------------


@pytest.mark.parametrize("g1_fn,g2_fn", [
    (_g_unipotent, _g_torus),
    (_g_generic, _g_weyl),
    (_g_torus, _g_unipotent),
    (_g_weyl, _g_weyl),
])
@pytest.mark.parametrize("V_fn", [_pw_standard_rep, _pw_sym2_rep])
def test_sl2_group_homomorphism(g1_fn, g2_fn, V_fn):
    g1 = g1_fn()
    g2 = g2_fn()
    V = V_fn()
    v = verify_sl2_group_homomorphism(g1, g2, V)
    assert v["bit_exact"]


def test_sl2_group_hom_with_identity_left():
    e = _g_identity()
    g = _g_generic()
    v = verify_sl2_group_homomorphism(e, g, _pw_standard_rep())
    assert v["bit_exact"]


def test_sl2_group_hom_with_identity_right():
    e = _g_identity()
    g = _g_generic()
    v = verify_sl2_group_homomorphism(g, e, _pw_sym2_rep())
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# G_a x SL_2 commute
# ---------------------------------------------------------------------


@pytest.mark.parametrize("g_label,g_fn", [
    ("identity", _g_identity),
    ("unipotent", _g_unipotent),
    ("torus", _g_torus),
    ("generic", _g_generic),
])
def test_ga_sl2_commute_J2_standard(g_label, g_fn):
    J = FinDimRep(
        n_max=2, dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2",
    )
    V_pw = _pw_standard_rep()
    t = {(1, 0, 0): Rational(2, 3), (2, 0, 1): Rational(-1, 5)}
    v = verify_ga_sl2_commute(t, g_fn(), J, V_pw)
    assert v["bit_exact"]


@pytest.mark.parametrize("g_label,g_fn", [
    ("torus", _g_torus),
    ("generic", _g_generic),
])
def test_ga_sl2_commute_J3_sym2(g_label, g_fn):
    J = FinDimRep(
        n_max=2, dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3",
    )
    V_pw = _pw_sym2_rep()
    t = {(1, 0, 0): Rational(1, 4)}
    v = verify_ga_sl2_commute(t, g_fn(), J, V_pw)
    assert v["bit_exact"]


def test_ga_sl2_commute_zero_t():
    """When t = 0, G_a side is identity, so commutativity is trivial."""
    J = FinDimRep(
        n_max=2, dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2",
    )
    V_pw = _pw_sym2_rep()
    v = verify_ga_sl2_commute({}, _g_generic(), J, V_pw)
    assert v["bit_exact"]


def test_ga_sl2_commute_identity_g():
    """When g = identity, SL_2 side is identity, so commutativity is trivial."""
    J = FinDimRep(
        n_max=2, dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2",
    )
    V_pw = _pw_standard_rep()
    t = {(1, 0, 0): Rational(7, 3)}
    v = verify_ga_sl2_commute(t, _g_identity(), J, V_pw)
    assert v["bit_exact"]


# ---------------------------------------------------------------------
# Generator-rep construction
# ---------------------------------------------------------------------


def test_primitive_generator_rep_dim_and_label():
    g = (1, 0, 0)
    V = primitive_generator_rep(2, g)
    assert V.dim == 2
    assert V.label == "V_1_0_0"
    assert V.n_max == 2


def test_primitive_generator_rep_endo_structure():
    g = (2, 1, 1)
    V = primitive_generator_rep(2, g)
    assert V.X(g) == Matrix([[0, 1], [0, 0]])
    other = (1, 0, 0) if g != (1, 0, 0) else (1, 1, 0)
    assert V.X(other) == Matrix.zeros(2, 2)


# ---------------------------------------------------------------------
# Full-panel injectivity at n_max = 2
# ---------------------------------------------------------------------


def test_injectivity_zero_t_gives_identity():
    """eta_{V_g}(0) = I for every generator g."""
    gens = primitive_generators(2)
    assert len(gens) == 15
    t = {}
    for g in gens:
        v = verify_injectivity_at_generator(t, g, 2)
        assert v["eta_is_identity"]
        assert v["bit_exact"]


def test_injectivity_single_t_detects_only_matching_generator():
    """t_g != 0 -> eta_{V_g}(t) != I and eta_{V_h}(t) = I for h != g."""
    gens = primitive_generators(2)
    for g_active in gens:
        t = {g_active: Rational(1)}
        for g_test in gens:
            v = verify_injectivity_at_generator(t, g_test, 2)
            if g_test == g_active:
                assert not v["eta_is_identity"]
            else:
                assert v["eta_is_identity"]
            assert v["bit_exact"]


def test_injectivity_generic_t_detects_active_generators():
    gens = primitive_generators(2)
    active = {gens[0], gens[5], gens[10]}
    t = {
        gens[0]: Rational(1, 2),
        gens[5]: Rational(-1, 3),
        gens[10]: Rational(7, 4),
    }
    for g in gens:
        v = verify_injectivity_at_generator(t, g, 2)
        if g in active:
            assert not v["eta_is_identity"]
        else:
            assert v["eta_is_identity"]
        assert v["bit_exact"]


def test_panel_size_at_n_max_2_is_15():
    """3 * N(2) = 3 * 5 = 15 generators."""
    gens = primitive_generators(2)
    assert len(gens) == 15
