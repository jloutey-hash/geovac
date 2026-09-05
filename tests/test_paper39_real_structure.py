"""Backing test for Paper 39 lem:real_structure-T (the KO-6 tensor triple).

Written as its own pass (CLAUDE.md §9), after the discharge added the lemma
identifying the limit as the Dąbrowski–Dossena tensor product of two odd
(KO-3) Camporesi–Higuchi real spectral triples, with real structure
J_ab = J_a ⊗ σ₁K ⊗ J_b.

Each assertion names, in its docstring, THE WRONG ANSWER IT REJECTS.  The
load-bearing, non-tautological test is `test_sigma1_is_forced`: the KO-6
signs come out right ONLY for the σ₁ middle factor, and the single factor
must genuinely satisfy KO-3 for the product to satisfy KO-6 at all.

Convention: an antilinear J = U·K (K = complex conjugation) satisfies
`J X = s · X J` as operators iff `U conj(X) = s · X U`.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.real_structure import build_J_full_dirac

SIGMA = {
    "1": np.eye(2, dtype=complex),
    "x": np.array([[0, 1], [1, 0]], dtype=complex),   # σ₁
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),  # σ₂
    "z": np.array([[1, 0], [0, -1]], dtype=complex),   # σ₃
}


def _single_factor(n_max: int):
    """(D, U) for one CH factor: D self-adjoint, J = U·K at KO-3."""
    basis = full_dirac_basis(n_max)
    D = np.asarray(camporesi_higuchi_full_dirac_matrix(basis), dtype=complex)
    U = np.asarray(build_J_full_dirac(n_max).U, dtype=complex)
    return D, U


def _kron3(A, B, C):
    return np.kron(np.kron(A, B), C)


def _build_product(n_a, n_b, W):
    """The Clifford-doubled product with middle antiunitary factor W·K.

    Returns (D_ab, U_ab, gamma_ab).  D_ab = D_a⊗σ₁⊗1 + 1⊗σ₂⊗D_b,
    gamma_ab = 1⊗σ₃⊗1, J_ab = (U_a ⊗ W ⊗ U_b)·K.
    """
    D_a, U_a = _single_factor(n_a)
    D_b, U_b = _single_factor(n_b)
    Ia, Ib = np.eye(D_a.shape[0]), np.eye(D_b.shape[0])
    D_ab = (_kron3(D_a, SIGMA["x"], Ib) + _kron3(Ia, SIGMA["y"], D_b))
    U_ab = _kron3(U_a, W, U_b)
    gamma_ab = _kron3(Ia, SIGMA["z"], Ib)
    return D_ab, U_ab, gamma_ab


def _signs(D_ab, U_ab, gamma_ab):
    """Measured KO signs (J², JD/DJ, Jγ/γJ) as +1/-1 or nan if inconsistent."""
    # J² = U conj(U)
    J2 = U_ab @ np.conj(U_ab)
    s_sq = _scalar_sign(J2)
    # J D = s' D J  <=>  U conj(D) = s' D U
    s_d = _relation_sign(U_ab @ np.conj(D_ab), D_ab @ U_ab)
    # J γ = s'' γ J  <=>  U conj(γ) = s'' γ U
    s_g = _relation_sign(U_ab @ np.conj(gamma_ab), gamma_ab @ U_ab)
    return s_sq, s_d, s_g


def _scalar_sign(M):
    for s in (+1.0, -1.0):
        if np.allclose(M, s * np.eye(M.shape[0]), atol=1e-9):
            return int(s)
    return float("nan")


def _relation_sign(lhs, rhs_op):
    for s in (+1.0, -1.0):
        if np.allclose(lhs, s * rhs_op, atol=1e-9):
            return int(s)
    return float("nan")


# --------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 2])
def test_single_factor_is_KO3(n: int) -> None:
    """REJECTS: a single factor that does not satisfy KO-3.

    The product's KO-6 is not automatic — it needs J_a² = -1 and
    J_a D_a = +D_a J_a on each factor.  If this failed, every KO-6 check
    below would be vacuous, so it is asserted first.
    """
    D, U = _single_factor(n)
    assert np.allclose(D, D.conj().T, atol=1e-9)                 # D self-adjoint
    assert np.allclose(U @ np.conj(U), -np.eye(U.shape[0]), atol=1e-9)   # J² = -1
    assert np.allclose(U @ np.conj(D), D @ U, atol=1e-9)         # J D = +D J


@pytest.mark.parametrize("n_a,n_b", [(1, 1), (1, 2), (2, 2)])
def test_ko6_signs(n_a: int, n_b: int) -> None:
    """REJECTS: any KO-sign triple other than (+,+,-).

    The whole point of the discharge's limit-object claim: the doubled
    product is a KO-dimension-6 real spectral triple, signs (J², JD, Jγ) =
    (+1, +1, -1).  A wrong grading or Dirac assembly moves one of these.
    """
    s_sq, s_d, s_g = _signs(*_build_product(n_a, n_b, SIGMA["x"]))
    assert (s_sq, s_d, s_g) == (+1, +1, -1), (s_sq, s_d, s_g)


@pytest.mark.parametrize("n_a,n_b", [(1, 1), (2, 2)])
def test_sigma1_is_forced(n_a: int, n_b: int) -> None:
    """REJECTS: any middle factor other than σ₁ passing as the real structure.

    This is the non-tautological core.  Of the four antiunitaries W·K,
    W ∈ {1, σ₁, σ₂, σ₃}, ONLY σ₁ yields the KO-6 row against the ACTUAL CH
    Dirac: σ₂ gives J² = -1 (wrong parity); 1 and σ₃ fail J D = +D J on the
    second Dirac term (the D-dependent axis the suspiciously-clean heuristic
    targets).  A construction that got KO-6 for more than one W would mean the
    D-dependent relations were not really being tested.
    """
    passing = []
    for name, W in SIGMA.items():
        if (s := _signs(*_build_product(n_a, n_b, W))) == (+1, +1, -1):
            passing.append(name)
    assert passing == ["x"], f"expected only sigma_1, got {passing}"
    # and confirm the specific failure modes the lemma names
    assert _signs(*_build_product(n_a, n_b, SIGMA["y"]))[0] == -1   # σ₂: J²=-1
    assert _signs(*_build_product(n_a, n_b, SIGMA["1"]))[1] != +1   # 1: JD fails
    assert _signs(*_build_product(n_a, n_b, SIGMA["z"]))[1] != +1   # σ₃: JD fails


def test_grading_squares_to_one_and_anticommutes_with_dirac() -> None:
    """REJECTS: a grading that is not a Z2 grading of the product Dirac.

    γ_ab = 1⊗σ₃⊗1 must satisfy γ² = 1 and {γ, D_ab} = 0 (the two Dirac terms
    carry σ₁, σ₂, both anticommuting with σ₃).  A wrong Pauli slot breaks the
    anticommutation.
    """
    D_ab, _, gamma_ab = _build_product(2, 2, SIGMA["x"])
    assert np.allclose(gamma_ab @ gamma_ab, np.eye(gamma_ab.shape[0]), atol=1e-9)
    assert np.allclose(gamma_ab @ D_ab + D_ab @ gamma_ab,
                       np.zeros_like(D_ab), atol=1e-9)


def test_dirac_squares_block_diagonal() -> None:
    """REJECTS: cross terms in D²; the two Dirac summands must anticommute.

    D_ab² = D_a²⊗1⊗1 + 1⊗1⊗D_b² (no cross term) iff
    {D_a⊗σ₁⊗1, 1⊗σ₂⊗D_b} = 0, i.e. iff σ₁, σ₂ anticommute — the doubling's
    reason for existing.
    """
    D_a, _ = _single_factor(2)
    D_b, _ = _single_factor(2)
    Ia, Ib = np.eye(D_a.shape[0]), np.eye(D_b.shape[0])
    A = _kron3(D_a, SIGMA["x"], Ib)
    B = _kron3(Ia, SIGMA["y"], D_b)
    assert np.allclose(A @ B + B @ A, np.zeros_like(A), atol=1e-9)
    expected = _kron3(D_a @ D_a, np.eye(2), Ib) + _kron3(Ia, np.eye(2), D_b @ D_b)
    assert np.allclose((A + B) @ (A + B), expected, atol=1e-9)


# --------------------------------------------------------------------------
# Order-zero / order-one for the tensor triple (the reviewer's soft spot).
# Added 2026-09-05 after measuring that the finite-cutoff residual GROWS, not
# vanishes -- so the continuum claim rests on commutativity of the limit
# algebra (Paper 32), and these tests guard the corrected facts.

from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem
from geovac.connes_axiom_audit_31 import verify_order_zero, verify_order_one


def test_order01_residual_grows_not_vanishes() -> None:
    """REJECTS: the claim that the finite-cutoff order-0/1 residual vanishes.

    The truncated operator system is NOT an algebra, and the max order-0/1
    residual GROWS with n_max (measured 0.055->0.078 order-0, 0.101->0.203
    order-1 from n=2 to n=3).  The continuum object satisfies order-0/1 for a
    different reason (its limit algebra is commutative, Paper 32), NOT because
    this residual shrinks.  A future edit asserting the residual vanishes would
    fail here.
    """
    def resid(n):
        ops = FullDiracTruncatedOperatorSystem(n)
        A = ops.multiplier_matrices
        U = np.asarray(build_J_full_dirac(n).U, dtype=complex)
        D, _ = _single_factor(n)
        _, _, m0 = verify_order_zero(U, A, tol=1e-9)
        _, _, m1 = verify_order_one(U, A, D, tol=1e-9)
        return m0, m1
    m0_2, m1_2 = resid(2)
    m0_3, m1_3 = resid(3)
    assert m0_2 > 1e-3 and m1_2 > 1e-3          # genuinely nonzero at finite cutoff
    assert m0_3 > m0_2 and m1_3 > m1_2          # GROWS, does not vanish


def test_cross_factor_order_conditions_vanish() -> None:
    """REJECTS: a cross-factor order-0/1 obstruction in the product.

    The proved factor-reduction (memo 4.1) says a cross pair -- a on factor a,
    b on factor b -- satisfies both conditions exactly, because the middle 1_2
    conjugates to 1_2 and the legs decouple.  A wrong J_ab middle factor or a
    Dirac cross term would break this.  Control: a same-factor pair is NOT
    forced to vanish (so the test is not vacuously zero).
    """
    na = nb = 2
    D_a, U_a = _single_factor(na)
    D_b, U_b = _single_factor(nb)
    Ma = FullDiracTruncatedOperatorSystem(na).multiplier_matrices
    Mb = FullDiracTruncatedOperatorSystem(nb).multiplier_matrices
    Ia, Ib = np.eye(D_a.shape[0]), np.eye(D_b.shape[0])
    s1, s2 = SIGMA["x"], SIGMA["y"]
    U_ab = _kron3(U_a, s1, U_b)
    D_ab = _kron3(D_a, s1, Ib) + _kron3(Ia, s2, D_b)

    def order_res(a, b):
        JbJ = U_ab @ b.conj() @ U_ab.T            # verifier's convention U conj(b) U^T
        o0 = np.max(np.abs(a @ JbJ - JbJ @ a))
        Da = D_ab @ a - a @ D_ab
        o1 = np.max(np.abs(Da @ JbJ - JbJ @ Da))
        return o0, o1

    # cross pairs: a on factor a, b on factor b -- must vanish
    worst0 = worst1 = 0.0
    for A in Ma[:4]:
        for B in Mb[:4]:
            a = _kron3(A, np.eye(2), Ib)
            b = _kron3(Ia, np.eye(2), B)
            o0, o1 = order_res(a, b)
            worst0, worst1 = max(worst0, o0), max(worst1, o1)
    assert worst0 < 1e-9, worst0
    assert worst1 < 1e-9, worst1

    # control: same-factor (both on a) is allowed to be nonzero for some pair
    same = 0.0
    for A in Ma[:6]:
        for B in Ma[:6]:
            a = _kron3(A, np.eye(2), Ib)
            b = _kron3(B, np.eye(2), Ib)
            same = max(same, order_res(a, b)[0])
    assert same > 1e-3, "same-factor order-0 never fails; cross-test may be vacuous"
