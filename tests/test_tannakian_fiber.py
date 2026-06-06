"""Tests for the TC-1d fiber functor omega: Rep_fin(H_GV) -> Vec_Q.

See: Sprint Q5'-Tannakian-Closure TC-1d
(debug/sprint_q5p_tc1d_fiber_functor_memo.md).
"""

import pytest
from sympy import Matrix, eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    direct_sum,
    fiber_functor_morphism,
    fiber_functor_object,
    kernel,
    cokernel,
    tensor_rep,
    trivial_rep,
    unit_object,
    verify_omega_faithful,
    verify_omega_preserves_cokernel,
    verify_omega_preserves_direct_sum,
    verify_omega_preserves_kernel,
    verify_omega_tensor_preservation,
    verify_omega_unit,
    zero_rep,
)


# ---------------------------------------------------------------------
# Fiber functor on objects and morphisms
# ---------------------------------------------------------------------


def test_fiber_functor_object_returns_dim():
    T = trivial_rep(2, dim=3)
    assert fiber_functor_object(T) == 3


def test_fiber_functor_object_zero_rep():
    Z = zero_rep(2)
    assert fiber_functor_object(Z) == 0


def test_fiber_functor_morphism_returns_matrix():
    T = trivial_rep(2, dim=2)
    f = RepMorphism(T, T, eye(2))
    assert fiber_functor_morphism(f) == eye(2)


# ---------------------------------------------------------------------
# omega(1) = Q
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_unit(n_max):
    v = verify_omega_unit(n_max)
    assert v["bit_exact"]
    assert v["omega_unit_dim"] == 1


# ---------------------------------------------------------------------
# Tensor preservation
# ---------------------------------------------------------------------


def _jordan(n_max, dim=2):
    if dim == 2:
        M = Matrix([[0, 1], [0, 0]])
    elif dim == 3:
        M = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    else:
        raise ValueError("unsupported dim")
    return FinDimRep(n_max=n_max, dim=dim, endos={(1, 0, 0): M}, label=f"J{dim}")


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_tensor_preservation_trivial_jordan(n_max):
    T1 = unit_object(n_max)
    J = _jordan(n_max, 2)
    v = verify_omega_tensor_preservation(T1, J)
    assert v["bit_exact"]
    assert v["omega_M_tensor_N"] == 2


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_tensor_preservation_jordan_jordan(n_max):
    J = _jordan(n_max, 2)
    v = verify_omega_tensor_preservation(J, J)
    assert v["bit_exact"]
    assert v["omega_M_tensor_N"] == 4


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_tensor_preservation_jordan_jordan3(n_max):
    J = _jordan(n_max, 2)
    J3 = _jordan(n_max, 3)
    v = verify_omega_tensor_preservation(J, J3)
    assert v["bit_exact"]
    assert v["omega_M_tensor_N"] == 6


# ---------------------------------------------------------------------
# Kernel and cokernel preservation
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_preserves_kernel_of_projection(n_max):
    J = _jordan(n_max, 2)
    T = trivial_rep(n_max, dim=1)
    f = RepMorphism(J, T, Matrix([[0, 1]]))  # e_2 -> 1, e_1 -> 0
    v = verify_omega_preserves_kernel(f)
    assert v["bit_exact"]
    assert v["omega_ker_dim"] == 1


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_preserves_cokernel_of_injection(n_max):
    J = _jordan(n_max, 2)
    T = trivial_rep(n_max, dim=1)
    f = RepMorphism(T, J, Matrix([[1], [0]]))  # 1 -> e_1
    v = verify_omega_preserves_cokernel(f)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_preserves_kernel_of_zero(n_max):
    T = trivial_rep(n_max, dim=2)
    f = RepMorphism(T, T, Matrix.zeros(2, 2))
    v = verify_omega_preserves_kernel(f)
    assert v["bit_exact"]
    # ker of zero morphism is whole source
    assert v["omega_ker_dim"] == 2


# ---------------------------------------------------------------------
# Direct sum preservation
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_preserves_direct_sum_triv_triv(n_max):
    T1 = trivial_rep(n_max, dim=1)
    T1b = trivial_rep(n_max, dim=1)
    v = verify_omega_preserves_direct_sum(T1, T1b)
    assert v["bit_exact"]
    assert v["omega_R1_plus_R2"] == 2


@pytest.mark.parametrize("n_max", [2, 3])
def test_omega_preserves_direct_sum_jordan_pair(n_max):
    J = _jordan(n_max, 2)
    J3 = _jordan(n_max, 3)
    v = verify_omega_preserves_direct_sum(J, J3)
    assert v["bit_exact"]
    assert v["omega_R1_plus_R2"] == 5


# ---------------------------------------------------------------------
# Faithfulness
# ---------------------------------------------------------------------


def test_omega_faithful_identity_vs_zero():
    T = trivial_rep(2, dim=2)
    id_T = RepMorphism(T, T, eye(2))
    zero_T = RepMorphism(T, T, Matrix.zeros(2, 2))
    v = verify_omega_faithful(id_T, zero_T)
    assert v["bit_exact"]
    assert v["morphism_equal"] is False
    assert v["matrix_equal"] is False
    assert v["faithful_consistent"]


def test_omega_faithful_identity_vs_identity():
    T = trivial_rep(2, dim=2)
    id_T = RepMorphism(T, T, eye(2))
    id_T_b = RepMorphism(T, T, eye(2))
    v = verify_omega_faithful(id_T, id_T_b)
    assert v["bit_exact"]
    assert v["matrix_equal"] is True


def test_omega_faithful_shape_mismatch():
    T2 = trivial_rep(2, dim=2)
    T3 = trivial_rep(2, dim=3)
    f = RepMorphism(T2, T2, eye(2))
    g = RepMorphism(T3, T3, eye(3))
    v = verify_omega_faithful(f, g)
    assert v["shapes_match"] is False
