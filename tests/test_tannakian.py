"""Tests for geovac.tannakian — FinDimRep, RepMorphism, abelian axioms.

Verifies the rep category Rep_fin(H_GV(n_max)) is abelian via the
standard Deligne-Milne 1982 universal property checks at small reps.

See: Sprint Q5'-Tannakian-Closure TC-1a
(debug/sprint_q5p_tc1a_abelian_memo.md).
"""

import pytest
from sympy import Matrix, Rational, eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    cokernel,
    compose,
    direct_sum,
    kernel,
    trivial_rep,
    verify_cokernel_universal_property,
    verify_direct_sum_universal_property,
    verify_epi_eq_coker_ker,
    verify_kernel_universal_property,
    verify_mono_eq_ker_coker,
    verify_zero_object_axiom,
    zero_rep,
)


# ---------------------------------------------------------------------
# FinDimRep construction and validation
# ---------------------------------------------------------------------


def test_zero_rep_basics():
    Z = zero_rep(2)
    assert Z.dim == 0
    assert Z.is_zero_object()


def test_trivial_rep_basics():
    T = trivial_rep(2, dim=3)
    assert T.dim == 3
    assert not T.is_zero_object()
    assert len(T.non_zero_endos()) == 0


def test_FinDimRep_validates_nilpotency():
    # Identity matrix is NOT nilpotent
    with pytest.raises(ValueError, match="not nilpotent"):
        FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): eye(2)})


def test_FinDimRep_validates_commutativity():
    # Two non-commuting nilpotents (e.g., distinct upper triangulars)
    A = Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    B = Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])
    # A and B do commute (both nilpotent, but A*B and B*A = 0 since strict upper).
    # Actually for this AB = [[0,0,1],[0,0,0],[0,0,0]], BA = 0. Don't commute.
    # Verify the rejection:
    with pytest.raises(ValueError, match="do not commute"):
        FinDimRep(n_max=2, dim=3, endos={(1, 0, 0): A, (1, 1, 0): B})


def test_FinDimRep_rejects_unknown_generator():
    # (10, 5, 0) is not a primitive generator at n_max=2
    with pytest.raises(ValueError, match="not in primitive_generators"):
        FinDimRep(n_max=2, dim=2, endos={(10, 5, 0): Matrix.zeros(2, 2)})


def test_FinDimRep_rejects_bad_shape():
    with pytest.raises(ValueError, match="must be"):
        FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): Matrix.zeros(3, 3)})


# ---------------------------------------------------------------------
# RepMorphism
# ---------------------------------------------------------------------


def test_RepMorphism_validates_intertwining():
    # f: J2 -> J2 where J2 has X = [[0,1],[0,0]].
    # f = [[1, 0], [0, 2]] should NOT intertwine: f X = [[0, 1], [0, 0]],
    # X f = [[0, 2], [0, 0]]. Different.
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    f = Matrix([[1, 0], [0, 2]])
    with pytest.raises(ValueError, match="intertwining"):
        RepMorphism(R, R, f)


def test_RepMorphism_identity_intertwines():
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    f = RepMorphism(R, R, eye(2))  # identity intertwines
    assert f.is_injective()
    assert f.is_surjective()
    assert not f.is_zero()


def test_RepMorphism_zero_intertwines_trivially():
    R = trivial_rep(2, dim=2)
    f = RepMorphism(R, R, Matrix.zeros(2, 2))
    assert f.is_zero()
    assert not f.is_injective()


def test_compose_intertwines():
    T = trivial_rep(2, dim=1)
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    f = RepMorphism(T, R, Matrix([[1], [0]]))  # T -> R, sends 1 -> e_1
    g = RepMorphism(R, T, Matrix([[0, 1]]))   # R -> T, sends e_1 -> 0, e_2 -> 1
    h = compose(g, f)
    assert h.is_zero()  # 0 -> e_1 -> 0


# ---------------------------------------------------------------------
# Kernel
# ---------------------------------------------------------------------


def test_kernel_of_zero_morphism_is_source():
    T = trivial_rep(2, dim=2)
    f = RepMorphism(T, T, Matrix.zeros(2, 2))
    K, iota = kernel(f)
    assert K.dim == 2  # whole source


def test_kernel_of_identity_is_zero():
    T = trivial_rep(2, dim=2)
    f = RepMorphism(T, T, eye(2))
    K, iota = kernel(f)
    assert K.dim == 0  # only zero vector


def test_kernel_of_jordan_projection():
    # J2 = (e_1, e_2) with X(e_2) = e_1. Project onto second coord.
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    T = trivial_rep(2, dim=1)
    pi = RepMorphism(R, T, Matrix([[0, 1]]))  # e_2 -> 1, e_1 -> 0
    K, iota = kernel(pi)
    # ker = span(e_1), dim 1
    assert K.dim == 1
    # X restricted to span(e_1) is zero (since X(e_1) = 0)
    assert len(K.non_zero_endos()) == 0
    # iota is the inclusion e_1 -> e_1 in R
    assert iota.matrix == Matrix([[1], [0]])


# ---------------------------------------------------------------------
# Cokernel
# ---------------------------------------------------------------------


def test_cokernel_of_zero_morphism_is_target():
    T = trivial_rep(2, dim=2)
    f = RepMorphism(T, T, Matrix.zeros(2, 2))
    C, pi = cokernel(f)
    assert C.dim == 2


def test_cokernel_of_identity_is_zero():
    T = trivial_rep(2, dim=2)
    f = RepMorphism(T, T, eye(2))
    C, pi = cokernel(f)
    assert C.dim == 0


# ---------------------------------------------------------------------
# Direct sum
# ---------------------------------------------------------------------


def test_direct_sum_dim():
    T = trivial_rep(2, dim=1)
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    S, i1, i2, p1, p2 = direct_sum(T, R)
    assert S.dim == 3


def test_direct_sum_universal_property_T1_T1():
    T1 = trivial_rep(2, dim=1)
    T1b = trivial_rep(2, dim=1)
    S, i1, i2, p1, p2 = direct_sum(T1, T1b)
    v = verify_direct_sum_universal_property(T1, T1b, S, i1, i2, p1, p2)
    assert v["bit_exact"]


def test_direct_sum_block_diagonal_endos():
    T = trivial_rep(2, dim=1)
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    S, _, _, _, _ = direct_sum(R, T)
    # S has dim 3, X^S = block diag(J, 0) = [[0,1,0],[0,0,0],[0,0,0]]
    expected = Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    assert S.X((1, 0, 0)) == expected


# ---------------------------------------------------------------------
# Abelian axioms (end-to-end via the verify_* helpers)
# ---------------------------------------------------------------------


def _build_jordan_test_pair(n_max: int):
    J = Matrix([[0, 1], [0, 0]])
    R = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): J}, label="J2")
    T = trivial_rep(n_max, dim=1)
    return T, R


@pytest.mark.parametrize("n_max", [2, 3])
def test_kernel_universal_property_passes(n_max):
    T, R = _build_jordan_test_pair(n_max)
    pi = RepMorphism(R, T, Matrix([[0, 1]]))
    K, iota = kernel(pi)
    v = verify_kernel_universal_property(pi, K, iota)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_cokernel_universal_property_passes(n_max):
    T, R = _build_jordan_test_pair(n_max)
    pi = RepMorphism(R, T, Matrix([[0, 1]]))
    C, pi_q = cokernel(pi)
    v = verify_cokernel_universal_property(pi, C, pi_q)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_mono_eq_ker_coker(n_max):
    T, R = _build_jordan_test_pair(n_max)
    # iota: T -> R sending 1 -> e_1
    iota = RepMorphism(T, R, Matrix([[1], [0]]))
    v = verify_mono_eq_ker_coker(iota)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_epi_eq_coker_ker(n_max):
    T, R = _build_jordan_test_pair(n_max)
    pi = RepMorphism(R, T, Matrix([[0, 1]]))
    v = verify_epi_eq_coker_ker(pi)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_zero_object_axiom_passes(n_max):
    T, R = _build_jordan_test_pair(n_max)
    v_T = verify_zero_object_axiom(T)
    v_R = verify_zero_object_axiom(R)
    assert v_T["bit_exact"]
    assert v_R["bit_exact"]


# =====================================================================
# TC-1b --- Symmetric monoidal structure
# =====================================================================


from geovac.tannakian import (
    associator,
    braiding,
    left_unitor,
    right_unitor,
    tensor_morphism,
    tensor_rep,
    unit_object,
    verify_associator_intertwines,
    verify_braiding_intertwines,
    verify_braiding_symmetric,
    verify_hexagon_coherence,
    verify_pentagon_coherence,
    verify_tensor_diagonal_action,
    verify_tensor_functoriality,
    verify_triangle_coherence,
    verify_unitor_intertwines,
)


def test_tensor_rep_dim_is_product():
    T = trivial_rep(2, dim=2)
    J = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    R = tensor_rep(T, J)
    assert R.dim == 4


def test_tensor_rep_diagonal_action_T1_otimes_M_is_M():
    # 1 ⊗ M has the same endomorphisms as M (because 1 has zero action).
    T1 = unit_object(2)
    J = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    R = tensor_rep(T1, J)
    # X^R = 0 ⊗ I + I_1 ⊗ J = J
    assert R.X((1, 0, 0)) == Matrix([[0, 1], [0, 0]])


def test_tensor_rep_diagonal_action_J_otimes_J():
    J_mat = Matrix([[0, 1], [0, 0]])
    J = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J_mat})
    R = tensor_rep(J, J)
    # X^R = J ⊗ I_2 + I_2 ⊗ J (4x4 Kronecker sum).
    expected_J_kron_I = Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    expected_I_kron_J = Matrix([
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
    ])
    assert R.X((1, 0, 0)) == expected_J_kron_I + expected_I_kron_J


@pytest.mark.parametrize("n_max", [2, 3])
def test_tensor_diagonal_action_panel(n_max):
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    T = trivial_rep(n_max, dim=1)
    v = verify_tensor_diagonal_action(J, T)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_unitor_intertwines_for_jordan(n_max):
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    v = verify_unitor_intertwines(J)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_associator_intertwines(n_max):
    T = trivial_rep(n_max, dim=1)
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    v = verify_associator_intertwines(T, J, J)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_braiding_intertwines(n_max):
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    T = trivial_rep(n_max, dim=2)
    v = verify_braiding_intertwines(J, T)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_braiding_symmetric(n_max):
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    J3 = FinDimRep(n_max=n_max, dim=3, endos={
        (1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    })
    v = verify_braiding_symmetric(J, J3)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_pentagon_coherence(n_max):
    T = trivial_rep(n_max, dim=1)
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    v = verify_pentagon_coherence(T, J, J, T)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_triangle_coherence(n_max):
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    J3 = FinDimRep(n_max=n_max, dim=3, endos={
        (1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    })
    v = verify_triangle_coherence(J, J3)
    assert v["bit_exact"]


@pytest.mark.parametrize("n_max", [2, 3])
def test_hexagon_coherence(n_max):
    T = trivial_rep(n_max, dim=1)
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    J3 = FinDimRep(n_max=n_max, dim=3, endos={
        (1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    })
    v = verify_hexagon_coherence(J, J3, T)
    assert v["bit_exact"]


def test_tensor_functoriality_with_identities():
    n_max = 2
    T = trivial_rep(n_max, dim=1)
    J = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])})
    id_T = RepMorphism(T, T, eye(1), validate=False)
    id_J = RepMorphism(J, J, eye(2), validate=False)
    v = verify_tensor_functoriality(id_T, id_J, id_T, id_J)
    assert v["bit_exact"]


def test_unit_object_is_trivial_dim_1():
    one = unit_object(2)
    assert one.dim == 1
    assert len(one.non_zero_endos()) == 0
