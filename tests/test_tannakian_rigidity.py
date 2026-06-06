r"""Tests for geovac.tannakian TC-1c additions --- rigidity.

Verifies the rep category Rep_fin(H_GV(n_max)) is rigid: every finite-dim
rep V admits a dual V^vee = Hom_Q(V, Q) with the contragredient Hopf action
X_g^{V^vee} = -(X_g^V)^T, evaluation ev_V: V^vee otimes V -> 1, coevaluation
coev_V: 1 -> V otimes V^vee, and the two snake (zigzag) identities.

See: Sprint Q5'-Tannakian-Closure TC-1c
(debug/sprint_q5p_tc1c_rigidity_memo.md).
"""

import pytest
from sympy import Matrix, eye

from geovac.tannakian import (
    FinDimRep,
    coevaluation_morphism,
    dual_rep,
    evaluation_morphism,
    trivial_rep,
    unit_object,
    verify_coevaluation_intertwines,
    verify_double_dual_iso,
    verify_dual_action,
    verify_evaluation_intertwines,
    verify_snake_identity_first,
    verify_snake_identity_second,
    verify_unit_self_dual,
)


# ---------------------------------------------------------------------
# dual_rep construction
# ---------------------------------------------------------------------


def test_dual_rep_trivial_is_trivial():
    """Dual of a trivial rep is trivial (no non-zero endos)."""
    T = trivial_rep(2, dim=3)
    T.label = "T3"
    Td = dual_rep(T)
    assert Td.dim == 3
    assert len(Td.non_zero_endos()) == 0


def test_dual_rep_jordan2_at_nmax2():
    """Dual of the 2-dim Jordan block has the negative transpose endo."""
    J = Matrix([[0, 1], [0, 0]])
    V = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J}, label="J2")
    Vd = dual_rep(V)
    assert Vd.dim == 2
    # -(J)^T = -[[0,1],[0,0]]^T = -[[0,0],[1,0]] = [[0,0],[-1,0]]
    expected = Matrix([[0, 0], [-1, 0]])
    assert Vd.X((1, 0, 0)) == expected


def test_dual_rep_preserves_nilpotency():
    """The contragredient action of a nilpotent endo is nilpotent."""
    J = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    V = FinDimRep(n_max=3, dim=3, endos={(1, 0, 0): J}, label="J3")
    Vd = dual_rep(V)
    Xd = Vd.X((1, 0, 0))
    # Xd^3 = 0
    assert (Xd ** 3) == Matrix.zeros(3, 3)


def test_dual_rep_preserves_commutativity():
    """The contragredient action preserves pairwise commutativity."""
    # Two commuting nilpotent matrices: strict-upper-triangular E_{12} and E_{13}.
    A = Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    B = Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]])
    # AB = 0, BA = 0 (both strict upper, "to row 1 only")
    assert A * B == B * A
    V = FinDimRep(n_max=2, dim=3, endos={(1, 0, 0): A, (2, 0, 0): B}, label="V")
    Vd = dual_rep(V)
    Ad = Vd.X((1, 0, 0))
    Bd = Vd.X((2, 0, 0))
    assert Ad * Bd == Bd * Ad


# ---------------------------------------------------------------------
# Contragredient action verifier
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_verify_dual_action_jordan(n_max):
    J = Matrix([[0, 1], [0, 0]])
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): J}, label="J2")
    result = verify_dual_action(V)
    assert result["bit_exact"] is True
    assert result["mismatches"] == []
    assert result["nilpotency_ok"] is True
    assert result["commutativity_ok"] is True


@pytest.mark.parametrize("n_max", [2, 3])
def test_verify_dual_action_trivial(n_max):
    """Trivial rep has no non-zero endos; the audit is vacuously bit-exact."""
    T = trivial_rep(n_max, dim=2)
    T.label = "T2"
    result = verify_dual_action(T)
    assert result["bit_exact"] is True


# ---------------------------------------------------------------------
# Evaluation and coevaluation shapes
# ---------------------------------------------------------------------


def test_evaluation_morphism_shape_and_entries():
    """ev: V^vee ⊗ V → 1 is a 1 x dim(V)^2 matrix with deltas on the diagonal positions."""
    V = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    ev = evaluation_morphism(V)
    assert ev.matrix.rows == 1
    assert ev.matrix.cols == 4
    assert ev.matrix[0, 0] == 1  # phi_0 ⊗ e_0
    assert ev.matrix[0, 1] == 0  # phi_0 ⊗ e_1
    assert ev.matrix[0, 2] == 0  # phi_1 ⊗ e_0
    assert ev.matrix[0, 3] == 1  # phi_1 ⊗ e_1


def test_coevaluation_morphism_shape_and_entries():
    """coev: 1 → V ⊗ V^vee is a dim(V)^2 x 1 matrix with deltas on diagonal positions."""
    V = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    coev = coevaluation_morphism(V)
    assert coev.matrix.rows == 4
    assert coev.matrix.cols == 1
    assert coev.matrix[0, 0] == 1  # e_0 ⊗ phi_0
    assert coev.matrix[1, 0] == 0  # e_0 ⊗ phi_1
    assert coev.matrix[2, 0] == 0  # e_1 ⊗ phi_0
    assert coev.matrix[3, 0] == 1  # e_1 ⊗ phi_1


# ---------------------------------------------------------------------
# Intertwining verifiers
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_verify_evaluation_intertwines_jordan(n_max):
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    result = verify_evaluation_intertwines(V)
    assert result["bit_exact"] is True
    assert result["mismatches"] == []


@pytest.mark.parametrize("n_max", [2, 3])
def test_verify_coevaluation_intertwines_jordan(n_max):
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    result = verify_coevaluation_intertwines(V)
    assert result["bit_exact"] is True
    assert result["mismatches"] == []


@pytest.mark.parametrize("n_max", [2, 3])
def test_verify_evaluation_intertwines_trivial(n_max):
    T = trivial_rep(n_max, dim=3)
    T.label = "T3"
    result = verify_evaluation_intertwines(T)
    assert result["bit_exact"] is True


# ---------------------------------------------------------------------
# Snake identities
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_snake_first_identity_jordan2(n_max):
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    result = verify_snake_identity_first(V)
    assert result["bit_exact"] is True
    assert result["composite_shape"] == [2, 2]


@pytest.mark.parametrize("n_max", [2, 3])
def test_snake_second_identity_jordan2(n_max):
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    result = verify_snake_identity_second(V)
    assert result["bit_exact"] is True
    assert result["composite_shape"] == [2, 2]


def test_snake_first_identity_jordan3():
    V = FinDimRep(
        n_max=3, dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3",
    )
    assert verify_snake_identity_first(V)["bit_exact"] is True


def test_snake_second_identity_jordan3():
    V = FinDimRep(
        n_max=3, dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3",
    )
    assert verify_snake_identity_second(V)["bit_exact"] is True


@pytest.mark.parametrize("n_max", [2, 3])
def test_snake_identities_trivial_rep(n_max):
    """Snake identities hold for the trivial 2-dim rep (zero action)."""
    T = trivial_rep(n_max, dim=2)
    T.label = "T2"
    assert verify_snake_identity_first(T)["bit_exact"] is True
    assert verify_snake_identity_second(T)["bit_exact"] is True


@pytest.mark.parametrize("n_max", [2, 3])
def test_snake_identities_unit(n_max):
    """Snake identities hold for the tensor unit T1."""
    T1 = unit_object(n_max)
    T1.label = "T1"
    assert verify_snake_identity_first(T1)["bit_exact"] is True
    assert verify_snake_identity_second(T1)["bit_exact"] is True


# ---------------------------------------------------------------------
# Double dual and unit self-dual
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_double_dual_jordan2(n_max):
    V = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])}, label="J2")
    result = verify_double_dual_iso(V)
    assert result["bit_exact"] is True
    # (V^vee)^vee should have the same endo as V
    Vdd = dual_rep(dual_rep(V))
    assert Vdd.X((1, 0, 0)) == V.X((1, 0, 0))


@pytest.mark.parametrize("n_max", [2, 3])
def test_unit_self_dual(n_max):
    result = verify_unit_self_dual(n_max)
    assert result["bit_exact"] is True
    assert result["dim_match"] is True
    assert result["endos_match"] is True
