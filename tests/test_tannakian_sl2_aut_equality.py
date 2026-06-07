"""Tests for the TC-2b SL_2 reconstruction on the PW panel.

See: Sprint Q5'-Tannakian-Closure TC-2b
(debug/sprint_q5p_tc2b_sl2_aut_equality_memo.md).

Verifies that dim Aut^⊗(ω) on the PW panel = 3 = dim SL_2 bit-exactly,
via the witness-variety method on the tensor decomposition
V_fund ⊗ V_fund = Sym^2 ⊕ V_triv.
"""

import importlib.util
from pathlib import Path

from sympy import Matrix, Symbol, Rational, Integer, eye, zeros

from geovac.tannakian import _sl2_sym2_action


def _load_driver():
    """Import the TC-2b driver as a module so we can test its helpers."""
    repo_root = Path(__file__).parent.parent
    driver_path = repo_root / "debug" / "compute_q5p_tc2b_sl2_aut_equality.py"
    spec = importlib.util.spec_from_file_location("tc2b_driver", driver_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


driver = _load_driver()


# ---------------------------------------------------------------------
# Phi decomposition isomorphism
# ---------------------------------------------------------------------


def test_phi_inverts_correctly():
    """Φ · Φ^{-1} = I_4 bit-exact."""
    Phi, Phi_inv = driver.build_phi_decomposition()
    assert Phi * Phi_inv == eye(4)
    assert Phi_inv * Phi == eye(4)


def test_phi_decomposition_identity_invariant():
    """For η = I, Φ · (I ⊗ I) · Φ^{-1} = I_4."""
    Phi, Phi_inv = driver.build_phi_decomposition()
    I = eye(2)
    I_tensor = driver._kron_q(I, I)
    result = (Phi * I_tensor * Phi_inv).expand()
    assert result == eye(4)


# ---------------------------------------------------------------------
# Tensor compatibility forces det = 1
# ---------------------------------------------------------------------


def test_top_left_block_equals_sym2_eta():
    """The 3x3 top-left block of Φ · (η ⊗ η) · Φ^{-1} equals Sym^2(η)."""
    result = driver.run_sl2_panel_test()
    assert result["top_left_eq_sym2_eta"] is True


def test_bottom_right_block_equals_det_eta():
    """The 1x1 bottom-right block of Φ · (η ⊗ η) · Φ^{-1} equals det(η)."""
    result = driver.run_sl2_panel_test()
    assert result["bottom_right_eq_det_eta"] is True


def test_off_diagonal_blocks_are_zero():
    """The Sym^2 ↔ V_triv off-diagonal blocks vanish bit-exactly."""
    result = driver.run_sl2_panel_test()
    assert result["off_diagonal_blocks_zero"] is True


def test_det_constraint_is_ad_minus_bc_minus_one():
    """The naturality-derived constraint is exactly ad - bc - 1 = 0."""
    result = driver.run_sl2_panel_test()
    assert result["det_constraint"] == "a*d - b*c - 1"


# ---------------------------------------------------------------------
# Variety dim = 3
# ---------------------------------------------------------------------


def test_jacobian_rank_at_identity_is_1():
    result = driver.run_sl2_panel_test()
    assert result["jacobian_rank_at_identity"] == 1


def test_jacobian_rank_at_generic_is_1():
    """At any SL_2-point, the det = 1 constraint has codim 1."""
    result = driver.run_sl2_panel_test()
    assert result["jacobian_rank_at_generic"] == 1


def test_sl2_variety_dim_equals_3():
    """Headline: dim Aut^⊗(ω) on PW panel = 3 = dim SL_2."""
    result = driver.run_sl2_panel_test()
    assert result["variety_dim"] == 3
    assert result["match"] is True


# ---------------------------------------------------------------------
# Combined sanity check
# ---------------------------------------------------------------------


def test_combined_eta_invertibility():
    """A generic combined η = exp(q E_{12}) ⊗ SL_2 element is invertible."""
    result = driver.run_combined_sanity_check(n_max=2)
    assert result["invertibility_holds"] is True


def test_combined_h_gv_naturality():
    """Combined η commutes with the H_GV action (X_g ⊗ I) bit-exact.

    Follows from η_{V_g} commuting with X_g (since η_{V_g} = exp(q_g X_g)
    and X_g commutes with itself).
    """
    result = driver.run_combined_sanity_check(n_max=2)
    assert result["h_gv_naturality_bit_exact"] is True


def test_combined_dim_18():
    """At n_max = 2 the combined Aut^⊗ has predicted dim 15 + 3 = 18.

    By Deligne-Milne 1982 Theorem 2.3 (exterior tensor product of neutral
    Tannakian categories), Aut^⊗(C_1 ⊗ C_2) = Aut^⊗(C_1) × Aut^⊗(C_2).
    Combining TC-2a (Aut^⊗(C_{H_GV}) = G_a^{15}) with TC-2b Part A
    (Aut^⊗(C_{SL_2}) = SL_2), we get dim = 15 + 3 = 18.
    """
    result = driver.run_combined_sanity_check(n_max=2)
    assert result["n_max_axis_dim_from_tc2a"] == 15
    assert result["sl2_axis_dim"] == 3
    assert result["combined_dim_predicted"] == 18


# ---------------------------------------------------------------------
# Cross-check: _sl2_sym2_action is consistent with the symmetric square
# ---------------------------------------------------------------------


def test_sl2_sym2_action_identity():
    """_sl2_sym2_action(I) = I_3."""
    assert _sl2_sym2_action(eye(2)) == eye(3)


def test_sl2_sym2_action_unipotent():
    """For g = [[1, 1], [0, 1]], _sl2_sym2_action(g) is upper triangular
    with diag (1, 1, 1) and explicit off-diagonal pattern."""
    g = Matrix([[1, 1], [0, 1]])
    expected = Matrix([
        [1, 2, 1],
        [0, 1, 1],
        [0, 0, 1],
    ])
    assert _sl2_sym2_action(g) == expected
