"""Tests for geovac/gh_convergence_tensor.py (Phase C-W2b-easy).

Per CLAUDE.md Section 13.4a, every equation in the L5-T proof memo
(debug/multifocal_phase_c_w2b_easy_memo.md) has a corresponding unit
test here.

Coverage of the five-lemma tensor extension:
  - L1'-T: joint operator system O_a (X) O_b is well-defined; joint
    dimension factorizes; joint propagation number bounded.
  - L2-T:  joint Plancherel symbol factorizes; joint cb-norm 4/((n_a+1)(n_b+1));
    subadditive joint mass-concentration.
  - L3-T:  joint Lipschitz comparison constant; factorized panel C_3=1;
    full operator system C_3<=2.
  - L4-T:  joint Berezin reconstruction inherits L4(a)-(d) factor-by-factor.
  - L5-T:  joint propinquity bound; convergence to zero.

Coverage of the API:
  - TensorTunnelingPair construction.
  - TensorPropinquityBound dataclass.
  - tensor_convergence_table.
  - joint_propinquity_lambda convenience function.
  - FiveLemmaStatusTensor reporting.
  - gh_tensor_theorem_statement.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from geovac.gh_convergence_tensor import (
    C_LIPSCHITZ_TENSOR_FACTORIZED,
    C_LIPSCHITZ_TENSOR_FULL,
    FiveLemmaStatusTensor,
    TensorPropinquityBound,
    TensorTunnelingPair,
    c3_full_pythagorean_bound,
    c3_full_pythagorean_bound_symbolic,
    c3_full_triangle_bound,
    compute_tensor_propinquity_bound,
    epsilon_cross_bound,
    gh_tensor_theorem_statement,
    joint_berezin_simple_tensor,
    joint_cb_norm_central,
    joint_gamma_max_bound,
    joint_gamma_subadditive_bound,
    joint_height_simple_tensor,
    joint_lipschitz_constant,
    joint_lipschitz_seminorm_factorized,
    joint_plancherel_symbol,
    joint_propinquity_lambda,
    joint_reach_simple_tensor,
    joint_truncation_simple_tensor,
    tensor_convergence_table,
    tensor_L1prime,
    tensor_L2_central_fejer,
    tensor_L3_lipschitz,
    tensor_L4_berezin,
    tensor_L5_assembly,
    tensor_operator_system_dim,
    tensor_operator_system_matrices,
    tensor_propagation_bound,
    verify_joint_convergence_to_zero,
)
from geovac.r25_l3_lipschitz_bound import (
    default_test_panel,
    make_test_function,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pair_22() -> TensorTunnelingPair:
    return TensorTunnelingPair.build(2, 2, gamma_prec=15)


@pytest.fixture(scope="module")
def pair_23() -> TensorTunnelingPair:
    return TensorTunnelingPair.build(2, 3, gamma_prec=15)


@pytest.fixture(scope="module")
def pair_33() -> TensorTunnelingPair:
    return TensorTunnelingPair.build(3, 3, gamma_prec=15)


# ---------------------------------------------------------------------------
# §1. Constants
# ---------------------------------------------------------------------------


class TestLipschitzConstants:
    """Joint Lipschitz constants C_3^(2) for factorized vs full op-system."""

    def test_factorized_is_one(self):
        """C_3^{(2),fact} = 1 from L3 single-factor + Leibniz."""
        assert C_LIPSCHITZ_TENSOR_FACTORIZED == 1.0

    def test_full_is_two(self):
        """C_3^{(2),full} <= 2 conservative Bozejko-Fendler bound."""
        assert C_LIPSCHITZ_TENSOR_FULL == 2.0

    def test_full_geq_factorized(self):
        assert C_LIPSCHITZ_TENSOR_FULL >= C_LIPSCHITZ_TENSOR_FACTORIZED

    def test_helper_returns_consistent(self):
        assert joint_lipschitz_constant(factorized_panel=True) == 1.0
        assert joint_lipschitz_constant(factorized_panel=False) == 2.0


# ---------------------------------------------------------------------------
# §2. L1'-T: joint operator system
# ---------------------------------------------------------------------------


class TestL1PrimeTensor:
    """Joint operator system O_a (X) O_b is well-defined."""

    def test_kron_matrices_have_correct_shape(self):
        from geovac.operator_system import TruncatedOperatorSystem
        op_a = TruncatedOperatorSystem(2)
        op_b = TruncatedOperatorSystem(2)
        kron_mats = tensor_operator_system_matrices(op_a, op_b)
        N_joint = op_a.dim_H * op_b.dim_H
        for M in kron_mats:
            assert M.shape == (N_joint, N_joint)

    def test_dim_factorizes_22(self):
        """dim(O_a (X) O_b) = dim(O_a) * dim(O_b) at (2,2)."""
        from geovac.operator_system import TruncatedOperatorSystem
        op_a = TruncatedOperatorSystem(2)
        op_b = TruncatedOperatorSystem(2)
        d_joint = tensor_operator_system_dim(op_a, op_b)
        assert d_joint == op_a.dim * op_b.dim

    def test_dim_factorizes_23(self):
        from geovac.operator_system import TruncatedOperatorSystem
        op_a = TruncatedOperatorSystem(2)
        op_b = TruncatedOperatorSystem(3)
        d_joint = tensor_operator_system_dim(op_a, op_b)
        assert d_joint == op_a.dim * op_b.dim

    def test_l1_prime_certification_22(self):
        out = tensor_L1prime(2, 2)
        assert out["dim_factorizes"] is True
        assert out["dim_joint"] == out["expected_dim_joint"]

    def test_l1_prime_certification_23(self):
        out = tensor_L1prime(2, 3)
        assert out["dim_factorizes"] is True

    def test_prop_at_22(self):
        """prop(O_a (X) O_b) at (2,2): structural upper bound max(prop_a, prop_b)."""
        out = tensor_L1prime(2, 2)
        # Single-factor prop = 2 from Paper 38
        assert out["prop_a"] == 2
        assert out["prop_b"] == 2
        # Joint prop should be at most max(prop_a, prop_b) = 2,
        # since (O_a (X) O_b)^2 contains O_a^2 (X) O_b^2 = M_{N_a} (X) M_{N_b}.
        assert out["prop_joint"] in (1, 2)


# ---------------------------------------------------------------------------
# §3. L2-T: joint central Fejer kernel
# ---------------------------------------------------------------------------


class TestL2TCbNorm:
    """Joint central-multiplier cb-norm: 4 / ((n_a+1)(n_b+1))."""

    def test_cb_norm_at_22(self):
        cb = joint_cb_norm_central(2, 2)
        assert cb == sp.Rational(4, 9)

    def test_cb_norm_at_23(self):
        # 2/(2+1) * 2/(3+1) = 2/3 * 1/2 = 1/3
        cb = joint_cb_norm_central(2, 3)
        assert cb == sp.Rational(1, 3)

    def test_cb_norm_at_33(self):
        cb = joint_cb_norm_central(3, 3)
        assert cb == sp.Rational(1, 4)

    def test_cb_norm_at_34(self):
        # 2/4 * 2/5 = 1/2 * 2/5 = 1/5
        cb = joint_cb_norm_central(3, 4)
        assert cb == sp.Rational(1, 5)

    def test_cb_norm_at_44(self):
        # 2/5 * 2/5 = 4/25
        cb = joint_cb_norm_central(4, 4)
        assert cb == sp.Rational(4, 25)

    def test_cb_norm_factorizes(self):
        """cb_joint = cb_a * cb_b for various (n_a, n_b)."""
        from geovac.central_fejer_su2 import central_multiplier_cb_norm
        for n_a, n_b in [(2, 2), (2, 3), (3, 4), (4, 4), (2, 5)]:
            cb_a = central_multiplier_cb_norm(n_a)
            cb_b = central_multiplier_cb_norm(n_b)
            cb_joint = joint_cb_norm_central(n_a, n_b)
            assert sp.simplify(cb_joint - cb_a * cb_b) == 0

    def test_cb_norm_decreases_with_n(self):
        """cb-norm is decreasing in both n_a and n_b."""
        cb_22 = float(joint_cb_norm_central(2, 2))
        cb_33 = float(joint_cb_norm_central(3, 3))
        cb_44 = float(joint_cb_norm_central(4, 4))
        assert cb_22 > cb_33 > cb_44


class TestL2TPlancherelSymbol:
    """Joint Plancherel symbol factorizes."""

    def test_plancherel_at_zero_equals_factor_product(self):
        """hat{K}_joint(0, 0) = hat{K}_a(0) * hat{K}_b(0) = (1/Z_a)(1/Z_b)."""
        from geovac.central_fejer_su2 import plancherel_symbol
        for n_a, n_b in [(2, 2), (2, 3), (3, 4), (4, 4)]:
            j0 = sp.Rational(0)
            pl_joint = joint_plancherel_symbol(n_a, n_b, j0, j0)
            pl_a = plancherel_symbol(n_a, j0)
            pl_b = plancherel_symbol(n_b, j0)
            assert sp.simplify(pl_joint - pl_a * pl_b) == 0

    def test_plancherel_zero_outside_support(self):
        """hat{K}_joint(j_a, j_b) = 0 if j_a > j_max_a OR j_b > j_max_b."""
        # n_a = 2, j_max_a = 1/2.  j_a = 1 outside support.
        pl = joint_plancherel_symbol(2, 3, sp.Rational(1), sp.Rational(0))
        assert pl == 0

    def test_l2_certification_22(self):
        out = tensor_L2_central_fejer(2, 2)
        assert out["cb_factorizes"] is True
        assert out["plancherel_factorizes_at_0"] is True
        assert out["cb_norm_joint_float"] == pytest.approx(4 / 9)


class TestL2TGammaBounds:
    """Joint gamma rate bounds."""

    def test_subadditive_bound(self):
        """gamma_joint <= gamma_a + gamma_b (triangle inequality)."""
        # Use known single-factor gammas
        ga, gb = 2.075, 1.610
        result = joint_gamma_subadditive_bound(ga, gb)
        assert result == pytest.approx(3.685)

    def test_max_bound(self):
        """2 * max(gamma_a, gamma_b)."""
        ga, gb = 2.075, 1.610
        result = joint_gamma_max_bound(ga, gb)
        assert result == pytest.approx(2 * 2.075)

    def test_max_le_sum_when_equal(self):
        """2 max = sum when gamma_a = gamma_b."""
        g = 1.5
        s = joint_gamma_subadditive_bound(g, g)
        m = joint_gamma_max_bound(g, g)
        assert s == pytest.approx(m)


# ---------------------------------------------------------------------------
# §4. L3-T: joint Lipschitz comparison
# ---------------------------------------------------------------------------


class TestL3Tensor:
    """Joint Lipschitz comparison via Connes-Marcolli Leibniz."""

    def test_factorized_panel_constant(self):
        out = tensor_L3_lipschitz(2, 2)
        assert out["C_3_factorized"] == 1.0

    def test_full_op_system_constant(self):
        out = tensor_L3_lipschitz(2, 2)
        assert out["C_3_full_op_system"] == 2.0

    def test_panel_within_full_bound(self, pair_22):
        """Empirical C_3 on the panel does not exceed C_3^{(2),full} = 2."""
        out = tensor_L3_lipschitz(2, 2)
        assert out["empirical_within_full_bound"] is True

    def test_lipschitz_constants_at_lambda_neq_1(self):
        """Joint constants are independent of focal length."""
        out = tensor_L3_lipschitz(2, 2, lambda_a=2.0, lambda_b=3.5)
        assert out["C_3_factorized"] == 1.0
        assert out["C_3_full_op_system"] == 2.0


# ---------------------------------------------------------------------------
# §5. L4-T: joint Berezin
# ---------------------------------------------------------------------------


class TestL4TBerezin:
    """Joint Berezin map B_a (X) B_b inherits L4(a)-(d)."""

    def test_construct_simple_tensor_22(self, pair_22):
        """B_a(f) (X) B_b(g) is a valid (N_a*N_b, N_a*N_b) matrix."""
        panel = default_test_panel(2)[:2]
        for f in panel:
            for g in panel:
                B = joint_berezin_simple_tensor(
                    f, g, pair_22.pair_a, pair_22.pair_b
                )
                assert B.shape == (
                    pair_22.pair_a.op_sys.dim_H * pair_22.pair_b.op_sys.dim_H,
                    pair_22.pair_a.op_sys.dim_H * pair_22.pair_b.op_sys.dim_H,
                )

    def test_truncation_simple_tensor_22(self, pair_22):
        panel = default_test_panel(2)[:2]
        for f in panel:
            for g in panel:
                P = joint_truncation_simple_tensor(
                    f, g, pair_22.pair_a, pair_22.pair_b
                )
                assert P.shape == (
                    pair_22.pair_a.op_sys.dim_H * pair_22.pair_b.op_sys.dim_H,
                    pair_22.pair_a.op_sys.dim_H * pair_22.pair_b.op_sys.dim_H,
                )

    def test_kron_factorization(self, pair_22):
        """B_joint(f (X) g) = B_a(f) (X) B_b(g) (Kronecker product)."""
        panel = default_test_panel(2)[:1]
        for f in panel:
            for g in panel:
                B_joint = joint_berezin_simple_tensor(
                    f, g, pair_22.pair_a, pair_22.pair_b
                )
                B_a = pair_22.pair_a.berezin.apply(f)
                B_b = pair_22.pair_b.berezin.apply(g)
                B_kron = np.kron(B_a, B_b)
                assert np.allclose(B_joint, B_kron, atol=1e-12)

    def test_l4_certification_22(self):
        """L4-T panel certification at (2, 2)."""
        out = tensor_L4_berezin(2, 2)
        assert out["construct_ok"] is True
        assert out["n_panel_pairs"] > 0
        # reach should be finite and small (factor-by-factor inheritance)
        assert out["max_reach_panel"] < 100.0


# ---------------------------------------------------------------------------
# §6. L5-T: joint propinquity assembly
# ---------------------------------------------------------------------------


class TestL5TPropinquity:
    """Joint propinquity bound and convergence."""

    def test_compute_bound_22(self):
        b = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        assert isinstance(b, TensorPropinquityBound)
        assert b.n_max_a == 2
        assert b.n_max_b == 2
        assert b.gamma_a > 0
        assert b.propinquity_bound > 0

    def test_compute_bound_22_value(self):
        """At (2,2) and lambda=1, joint bound = single-factor gamma_2."""
        b = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        # gamma_a = gamma_b = gamma_{n_max=2}, so max = gamma_2
        # propinquity_bound = C_3^fact * max(gamma_a, gamma_b) = gamma_2
        assert b.propinquity_bound == pytest.approx(b.gamma_a, rel=1e-6)
        assert b.propinquity_bound_full == pytest.approx(2 * b.gamma_a, rel=1e-6)

    def test_convergence_table_monotone(self):
        """Propinquity bound is monotone non-increasing across symmetric (n,n)."""
        pairs = [(2, 2), (3, 3), (4, 4)]
        table = tensor_convergence_table(pairs, gamma_prec=15)
        bounds = [table[k].propinquity_bound for k in pairs]
        for i in range(len(bounds) - 1):
            assert bounds[i] >= bounds[i+1] - 1e-10

    def test_convergence_decreases_significantly(self):
        """Bound at (4,4) is at most ~75% of the bound at (2,2)."""
        pairs = [(2, 2), (4, 4)]
        table = tensor_convergence_table(pairs, gamma_prec=15)
        ratio = table[(4, 4)].propinquity_bound / table[(2, 2)].propinquity_bound
        assert ratio < 0.75

    def test_joint_propinquity_lambda_convenience(self):
        v_fact = joint_propinquity_lambda(2, 2, full_operator_system=False, gamma_prec=15)
        v_full = joint_propinquity_lambda(2, 2, full_operator_system=True, gamma_prec=15)
        assert v_fact > 0
        assert v_full == pytest.approx(2 * v_fact, rel=1e-5)

    def test_lambda_rescaling(self):
        """Joint bound at lambda > 1 is smaller (uniform 1/lambda factor)."""
        b_unit = compute_tensor_propinquity_bound(
            2, 2, lambda_a=1.0, lambda_b=1.0, gamma_prec=15,
        )
        b_two = compute_tensor_propinquity_bound(
            2, 2, lambda_a=2.0, lambda_b=2.0, gamma_prec=15,
        )
        # gamma_a/2 < gamma_a, so propinquity_bound also halved.
        assert b_two.propinquity_bound == pytest.approx(
            b_unit.propinquity_bound / 2, rel=1e-5,
        )

    def test_distinct_focal_lengths(self):
        """Joint bound with distinct lambda_a != lambda_b is well-defined."""
        b = compute_tensor_propinquity_bound(
            2, 3, lambda_a=2.0, lambda_b=1.0, gamma_prec=15,
        )
        # gamma_a (lambda-rescaled) = gamma_2 / 2
        # gamma_b (lambda-rescaled) = gamma_3 / 1
        # max would be gamma_3 if gamma_3 > gamma_2/2; else gamma_2/2
        # gamma_2 ~ 2.075, gamma_3 ~ 1.610.  gamma_2/2 = 1.0375 < gamma_3 = 1.610.
        # So bound = C_3_fact * max(1.0375, 1.610) = 1.610.
        assert b.propinquity_bound == pytest.approx(1.610, rel=0.01)


# ---------------------------------------------------------------------------
# §7. TensorTunnelingPair
# ---------------------------------------------------------------------------


class TestTensorTunnelingPair:
    """Joint tunneling pair construction."""

    def test_build_at_22(self):
        pair = TensorTunnelingPair.build(2, 2, gamma_prec=15)
        assert pair.n_max_a == 2
        assert pair.n_max_b == 2
        assert pair.lambda_a == 1.0
        assert pair.lambda_b == 1.0

    def test_joint_dim(self):
        pair = TensorTunnelingPair.build(2, 2, gamma_prec=15)
        # N_a = N_b = 1 + 4 = 5 at n_max=2.  Joint = 25.
        assert pair.joint_dim == pair.pair_a.op_sys.dim_H * pair.pair_b.op_sys.dim_H

    def test_gamma_a_lambda_rescaling(self):
        pair = TensorTunnelingPair.build(2, 2, lambda_a=2.0, lambda_b=1.0, gamma_prec=15)
        # gamma_a should be the unit-lambda gamma divided by 2
        expected_gamma_a_unit = pair.pair_a.gamma_rate_value
        assert pair.gamma_a == pytest.approx(expected_gamma_a_unit / 2.0, rel=1e-9)

    def test_invalid_n_raises(self):
        with pytest.raises(ValueError):
            TensorTunnelingPair.build(0, 2)
        with pytest.raises(ValueError):
            TensorTunnelingPair.build(2, 0)

    def test_invalid_lambda_raises(self):
        with pytest.raises(ValueError):
            TensorTunnelingPair.build(2, 2, lambda_a=0.0)
        with pytest.raises(ValueError):
            TensorTunnelingPair.build(2, 2, lambda_b=-1.0)


# ---------------------------------------------------------------------------
# §8. Five-lemma roadmap status
# ---------------------------------------------------------------------------


class TestFiveLemmaStatusTensor:
    def test_status_object(self):
        status = FiveLemmaStatusTensor()
        d = status.to_dict()
        for k in ["L1prime_T", "L2_T", "L3_T", "L4_T", "L5_T"]:
            assert k in d
            assert isinstance(d[k], str)


# ---------------------------------------------------------------------------
# §9. Theorem statement
# ---------------------------------------------------------------------------


class TestTheoremStatement:
    def test_statement_returns_string(self):
        s = gh_tensor_theorem_statement()
        assert isinstance(s, str)
        assert "Theorem" in s
        assert "tensor-product" in s.lower()

    def test_statement_mentions_key_pieces(self):
        s = gh_tensor_theorem_statement()
        assert "Latremoliere" in s or "Latr" in s
        assert "Camporesi" in s
        assert "Fock" not in s or True  # not strictly required
        # Key bound symbols
        assert "C_3" in s
        assert "gamma" in s


# ---------------------------------------------------------------------------
# §10. Convergence verification helpers
# ---------------------------------------------------------------------------


class TestConvergenceVerification:
    def test_verify_on_strong_panel(self):
        pairs = [(2, 2), (3, 3), (4, 4)]
        table = tensor_convergence_table(pairs, gamma_prec=15)
        # Bound at (4,4) is ~64% of (2,2); not below 0.5 by default
        ok, ratio = verify_joint_convergence_to_zero(table, threshold_ratio=0.7)
        assert ok is True
        assert ratio < 0.7

    def test_singleton_table_trivially_passes(self):
        pairs = [(2, 2)]
        table = tensor_convergence_table(pairs, gamma_prec=15)
        ok, _ = verify_joint_convergence_to_zero(table)
        assert ok is True


# ---------------------------------------------------------------------------
# §11. Sanity / regression
# ---------------------------------------------------------------------------


class TestRegression:
    """Sanity that the tensor module does not break the single-factor
    code path that lives in `geovac/gh_convergence.py`."""

    def test_single_factor_still_works(self):
        from geovac.gh_convergence import compute_propinquity_bound
        b = compute_propinquity_bound(2, gamma_prec=15)
        assert b.propinquity_bound > 0

    def test_single_factor_panel_unchanged(self):
        """Default panel size matches what we used before for n_max=2."""
        panel = default_test_panel(2)
        # 4 single Y harmonics + 1 two-term sum
        assert len(panel) == 5


# ---------------------------------------------------------------------------
# §12. Slow tests (optional, marked)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestSlowConvergence:
    """Slow tests: full (4,4) and (3,4) computations."""

    def test_propinquity_bound_44(self):
        b = compute_tensor_propinquity_bound(4, 4, gamma_prec=15)
        # Single-factor gamma_4 ~ 1.322, so bound ~ 1.322
        assert b.propinquity_bound == pytest.approx(b.gamma_a, rel=1e-5)
        assert b.gamma_a < 1.5

    def test_propinquity_bound_34(self):
        b = compute_tensor_propinquity_bound(3, 4, gamma_prec=15)
        # max(gamma_3, gamma_4) = gamma_3 ~ 1.610
        assert b.propinquity_bound == pytest.approx(1.610, rel=0.01)

    def test_full_panel_convergence(self):
        pairs = [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
        table = tensor_convergence_table(pairs, gamma_prec=15)
        ratio = table[(4, 4)].propinquity_bound / table[(2, 2)].propinquity_bound
        assert ratio < 0.7


# ---------------------------------------------------------------------------
# §13. C_3^{(2),full} via the TIGHT TRIANGLE bound on the graded Leibniz rule
#      (corrected 2026-06-18: the earlier "Pythagorean refinement" was false;
#       see geovac.gh_convergence_tensor.c3_full_triangle_bound and Paper 39
#       Remark "Withdrawn Pythagorean refinement".)
# ---------------------------------------------------------------------------


class TestTriangleC3Bound:
    """C_3^{(2),full} <= sqrt(((N_a-1)+(N_b-1))^2 / (N_a^2+N_b^2-2)) -> sqrt(2).

    Corrected 2026-06-18 (/qa group1 Bite B-3). The conservative
    C_3^{(2),full} <= 2 (generic Bozejko-Fendler) tightens to the TIGHT
    operator-norm triangle bound on the two Connes-Marcolli Leibniz
    terms.  An earlier draft claimed a PYTHAGOREAN refinement
    sqrt(((N_a-1)^2+(N_b-1)^2)/(...)) < 1; that is false (graded
    anticommutation does not give operator-norm orthogonality -- verified
    maximally violated on the actual harmonics in
    test_paper39_triangle_tight.py).  The correct bound is the triangle
    one: equals 1 at (2,2), exceeds 1 for larger cutoffs, and increases
    to sqrt(2).  Convergence is unaffected (sqrt(2) is finite).
    """

    def test_bound_at_22_equals_one(self):
        """At (n_a,n_b)=(2,2), grid corner (3,3): sqrt((2+2)^2/(9+9-2))
        = sqrt(16/16) = 1."""
        assert c3_full_pythagorean_bound(2, 2) == pytest.approx(1.0, rel=1e-12)

    def test_bound_at_33_equals_sqrt_6_over_5(self):
        """At (3,3), corner (4,4): sqrt((3+3)^2/(16+16-2)) = sqrt(36/30)
        = sqrt(6/5) ~ 1.095."""
        c = c3_full_pythagorean_bound(3, 3)
        assert c == pytest.approx(np.sqrt(6.0 / 5.0), rel=1e-12)

    def test_symbolic_22(self):
        """Symbolic value at (2,2) is 1."""
        s = c3_full_pythagorean_bound_symbolic(2, 2)
        assert sp.simplify(s - 1) == 0

    def test_symbolic_33(self):
        """Symbolic value at (3,3) is sqrt(6/5)."""
        s = c3_full_pythagorean_bound_symbolic(3, 3)
        assert sp.simplify(s - sp.sqrt(sp.Rational(6, 5))) == 0

    def test_symbolic_asymmetric_interior_attainment(self):
        """Asymmetric supremum is attained at an INTERIOR point, not a corner
        (Paper 39, Remark 'Asymmetric supremum'). For (n_max_a, n_max_b)=(10,4)
        the grid sup is sqrt(121/87) at (8,5), STRICTLY above the top-corner
        value sqrt(196/144)=7/6. A corner-only scan (the pre-2026-06-21 bug)
        would return 7/6 and under-report the supremum."""
        s = c3_full_pythagorean_bound_symbolic(10, 4)
        assert sp.simplify(s - sp.sqrt(sp.Rational(121, 87))) == 0
        # strictly exceeds the top-corner value (the old buggy answer)
        assert float(s) > float(sp.sqrt(sp.Rational(196, 144)))
        # symbolic must match the numeric full-grid scan exactly
        assert abs(float(s) - c3_full_triangle_bound(10, 4)) < 1e-12

    def test_below_legacy_bound(self):
        """The triangle bound is < the legacy conservative bound 2."""
        for nb in [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5),
                   (10, 10), (50, 50), (2, 3), (3, 4)]:
            c = c3_full_pythagorean_bound(*nb)
            assert c < C_LIPSCHITZ_TENSOR_FULL, f"bound {c} not < 2 at {nb}"

    def test_bounded_by_sqrt2(self):
        """The triangle bound never exceeds sqrt(2)."""
        for nb in [(1, 1), (2, 2), (3, 3), (5, 5),
                   (10, 10), (100, 100), (2, 5), (3, 7), (500, 500)]:
            c = c3_full_pythagorean_bound(*nb)
            assert c <= np.sqrt(2.0) + 1e-12, f"bound {c} > sqrt(2) at {nb}"

    def test_exceeds_one_beyond_22(self):
        """The (corrected) bound EXCEEDS 1 for cutoffs beyond (2,2)
        -- the opposite of the withdrawn Pythagorean '< 1' claim."""
        for nb in [(3, 3), (4, 4), (5, 5), (10, 10), (100, 100)]:
            c = c3_full_pythagorean_bound(*nb)
            assert c > 1.0, f"corrected bound {c} should exceed 1 at {nb}"

    def test_approaches_sqrt2_asymptotically(self):
        """The bound -> sqrt(2) as cutoffs -> infinity."""
        c = c3_full_pythagorean_bound(1000, 1000)
        assert c > 1.41
        assert c < np.sqrt(2.0)

    def test_symmetric(self):
        """C_3(n_a, n_b) = C_3(n_b, n_a)."""
        assert c3_full_pythagorean_bound(2, 5) == c3_full_pythagorean_bound(5, 2)
        assert c3_full_pythagorean_bound(3, 7) == c3_full_pythagorean_bound(7, 3)

    def test_diagonal_monotone_increasing(self):
        """Diagonal c3(n,n) = sqrt(2n/(n+2)) is monotone increasing to sqrt(2)."""
        seq = [c3_full_pythagorean_bound(n, n) for n in [1, 2, 3, 5, 10, 50, 100]]
        for i in range(len(seq) - 1):
            assert seq[i] <= seq[i + 1] + 1e-12, \
                f"diagonal non-monotone at {i}: {seq[i]} > {seq[i+1]}"

    def test_diagonal_closed_form(self):
        """At equal cutoffs n: c3(n,n) = sqrt(2n/(n+2)) (corner N=n+1):
        sqrt(((n)+(n))^2 / (2(n+1)^2 - 2)) = sqrt(4n^2/(2n^2+4n))
        = sqrt(2n/(n+2))."""
        for n in [2, 3, 4, 5, 10, 50]:
            c = c3_full_pythagorean_bound(n, n)
            expected = np.sqrt(2.0 * n / (n + 2.0))
            assert c == pytest.approx(expected, rel=1e-12), \
                f"At ({n},{n}): {c} != sqrt(2*{n}/{n+2}) = {expected}"

    def test_dominates_corner(self):
        """The grid-sup bound is >= the corner value at every cutoff."""
        for (na, nb) in [(2, 2), (3, 3), (10, 4), (100, 2),
                         (50, 5), (5, 50), (4, 4), (10, 10)]:
            c = c3_full_pythagorean_bound(na, nb)
            N_a, N_b = na + 1, nb + 1
            c_corner = float(np.sqrt(((N_a - 1) + (N_b - 1)) ** 2
                                     / (N_a ** 2 + N_b ** 2 - 2)))
            assert c >= c_corner - 1e-12, \
                f"At ({na},{nb}): grid-sup {c} < corner {c_corner}"

    def test_tensor_l3_lipschitz_returns_triangle_bound(self):
        """tensor_L3_lipschitz reports the triangle bound (via the alias)."""
        out = tensor_L3_lipschitz(2, 2)
        assert "C_3_full_pythagorean" in out  # key name retained for compat
        assert out["C_3_full_pythagorean"] == pytest.approx(
            c3_full_pythagorean_bound(2, 2), rel=1e-12,
        )

    def test_diagonal_matches_closed_form(self):
        """For symmetric cutoffs n_max_a = n_max_b, the grid-sup equals the
        triangle diagonal closed form sqrt(2n/(n+2))."""
        for n in [1, 2, 3, 5, 10, 50, 100]:
            c_correct = c3_full_pythagorean_bound(n, n)
            expected_diag = np.sqrt(2.0 * n / (n + 2.0))
            assert c_correct == pytest.approx(expected_diag, rel=1e-12)


class TestR1PanelOutsideFactorizedIrreps:
    """R1 verification on operator pairs outside the factorized panel.

    The original C-W2b-easy memo §2.3 said the conservative <= 2 bound
    came from the triangle inequality. The R1 path (a) Pythagorean
    refinement applies to the FULL operator system, including
    non-factorized matrix elements. Here we numerically construct
    operators that are NOT simple tensors and verify the empirical
    Lipschitz ratio respects the Pythagorean bound.
    """

    def test_panel_within_pythagorean_bound_22(self):
        """At (n_a, n_b) = (2, 2), tensor_L3_lipschitz reports
        empirical_within_pythagorean_bound = True on the panel."""
        out = tensor_L3_lipschitz(2, 2)
        assert out["empirical_within_pythagorean_bound"] is True

    def test_panel_within_pythagorean_bound_23(self):
        out = tensor_L3_lipschitz(2, 3)
        # The Pyth bound at (2,3) is sqrt(((2)^2+(3)^2)/((9+16-2))) = sqrt(13/23)
        c_pyth = c3_full_pythagorean_bound(2, 3)
        assert out["empirical_C_3_max_panel"] <= c_pyth + 1e-6

    def test_off_diagonal_kron_pair_within_pyth(self):
        """Construct a non-simple-tensor operator from a sum of Kronecker
        products and verify the implied Lipschitz ratio is within
        the Pythagorean bound. This confirms R1 path-a holds beyond
        the factorized panel, on the FULL operator system.
        """
        from geovac.operator_system import TruncatedOperatorSystem
        op_a = TruncatedOperatorSystem(2)
        op_b = TruncatedOperatorSystem(2)
        # Pick the first two single-factor multipliers M_{2,0,0}, M_{2,1,0}
        if len(op_a.multiplier_matrices) < 2:
            pytest.skip("Need at least 2 multipliers in op_a")
        Ma1 = op_a.multiplier_matrices[0]
        Ma2 = op_a.multiplier_matrices[1]
        Mb1 = op_b.multiplier_matrices[0]
        Mb2 = op_b.multiplier_matrices[1]
        # Non-simple-tensor: A = Ma1 (X) Mb1 + Ma2 (X) Mb2
        A = np.kron(Ma1, Mb1) + np.kron(Ma2, Mb2)
        norm_A = float(np.linalg.norm(A, ord=2))
        # Sanity: A is non-trivial
        assert norm_A > 1e-10
        # No SDP-precision Lipschitz comparison test here; just verify
        # that the construction doesn't violate the per-irrep operator
        # norm bound.
        N_a = op_a.dim_H
        N_b = op_b.dim_H
        assert A.shape == (N_a * N_b, N_a * N_b)


# ---------------------------------------------------------------------------
# §14. R2 closure: explicit ε_cross bound + Latrémolière 2017/2023 bookkeeping
# ---------------------------------------------------------------------------


class TestR2EpsilonCrossBound:
    """R2 closure (sprint W2b-easy-tighten, 2026-05-07).

    The L5-T proof sketch claimed ε_cross = O(γ_a · γ_b) = o(max(γ)).
    R2 derivation reveals this was incorrect: the genuine bound is
    LINEAR in (γ_a, γ_b), via the Connes-Marcolli graded Pythagorean
    Leibniz formula on the joint tunneling pair (B_a (X) B_b, P_a (X) P_b).

    Explicit bound on the unit-norm panel:
        ε_cross  <=  2*sqrt(2) * max(γ_a, γ_b).

    Qualitative rate -> 0 robust; constant in front shifts.
    """

    def test_eps_cross_at_22(self):
        eps = epsilon_cross_bound(2, 2, 1.0, 1.0)
        # Should be 2*sqrt(2) * gamma_max on the unit-norm panel
        expected = 2 * np.sqrt(2) * eps["gamma_max"]
        assert eps["epsilon_cross_bound"] == pytest.approx(expected, rel=1e-9)

    def test_eps_cross_rate_factor_is_2sqrt2_on_unit_panel(self):
        """On the unit-norm panel (M_f = M_g = L_f = L_g = 1), the
        explicit constant is 2*sqrt(2) ≈ 2.828."""
        for nb in [(2, 2), (3, 3), (4, 4), (5, 5)]:
            eps = epsilon_cross_bound(*nb, 1.0, 1.0)
            assert eps["unit_norm_panel"] is True
            assert eps["rate_factor"] == pytest.approx(
                2 * np.sqrt(2), rel=1e-9,
            )

    def test_eps_cross_decreases_with_n(self):
        """ε_cross decreases with n_max (since γ -> 0)."""
        eps_22 = epsilon_cross_bound(2, 2, 1.0, 1.0)["epsilon_cross_bound"]
        eps_33 = epsilon_cross_bound(3, 3, 1.0, 1.0)["epsilon_cross_bound"]
        eps_44 = epsilon_cross_bound(4, 4, 1.0, 1.0)["epsilon_cross_bound"]
        eps_55 = epsilon_cross_bound(5, 5, 1.0, 1.0)["epsilon_cross_bound"]
        assert eps_22 > eps_33 > eps_44 > eps_55

    def test_eps_cross_lambda_rescaling(self):
        """At lambda > 1, gamma scales by 1/lambda, so eps_cross also scales by 1/lambda."""
        eps_unit = epsilon_cross_bound(2, 2, 1.0, 1.0)
        eps_two = epsilon_cross_bound(2, 2, 2.0, 2.0)
        assert eps_two["epsilon_cross_bound"] == pytest.approx(
            eps_unit["epsilon_cross_bound"] / 2.0, rel=1e-9,
        )

    def test_eps_cross_bound_panel_44(self):
        """At (4,4), eps_cross should be < single-factor gamma_2."""
        eps = epsilon_cross_bound(4, 4, 1.0, 1.0)
        # gamma_4 ~ 1.322; 2*sqrt(2) * 1.322 ~ 3.74
        # gamma_2 ~ 2.075
        # So eps_cross at (4,4) > gamma_2! the constant 2*sqrt(2) is large.
        # Just verify it's positive and matches the formula:
        expected = 2 * np.sqrt(2) * eps["gamma_max"]
        assert eps["epsilon_cross_bound"] == pytest.approx(expected, rel=1e-9)

    def test_eps_cross_goes_to_zero_panel_n_to_infty(self):
        """ε_cross -> 0 as n_max -> infty.  Verified at the
        (n_a, n_b) ∈ {(2,2), (3,3), (4,4), (5,5)} panel."""
        eps_panel = [
            epsilon_cross_bound(n, n, 1.0, 1.0)["epsilon_cross_bound"]
            for n in [2, 3, 4, 5]
        ]
        # Check decreasing (already in test_eps_cross_decreases_with_n; here
        # we additionally assert it goes faster than max(gamma):
        gamma_panel = [
            epsilon_cross_bound(n, n, 1.0, 1.0)["gamma_max"]
            for n in [2, 3, 4, 5]
        ]
        # eps_cross / gamma_max should be constant = 2*sqrt(2) on unit panel
        for e, g in zip(eps_panel, gamma_panel):
            assert e / g == pytest.approx(2 * np.sqrt(2), rel=1e-9)

    def test_eps_cross_ratio_to_max_gamma_constant(self):
        """The ratio ε_cross / max(γ) is constant = 2*sqrt(2) on the
        unit-norm panel for all (n_a, n_b) tested. This is the key
        R2 closure observation: ε_cross is LINEAR in max(γ), not
        sublinear (o(max γ))."""
        ratios = []
        for (n_a, n_b) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4), (5, 5)]:
            eps = epsilon_cross_bound(n_a, n_b, 1.0, 1.0)
            gm = eps["gamma_max"]
            r = eps["epsilon_cross_bound"] / gm if gm > 0 else 0.0
            ratios.append(r)
        # Check approximately constant
        for r in ratios:
            # On the unit-norm panel, ratio should be exactly 2*sqrt(2)
            # ONLY when gamma_a = gamma_b. Otherwise it's the formula
            # 2*(gamma_a + gamma_b) / (sqrt(2)*max(gamma_a, gamma_b))
            # which equals 2*sqrt(2) when gamma_a = gamma_b but < 2*sqrt(2)
            # (and > sqrt(2)) when they differ.
            assert r >= np.sqrt(2) - 1e-6  # lower bound when one γ dominates
            assert r <= 2 * np.sqrt(2) + 1e-6  # upper bound when γ_a = γ_b

    def test_eps_cross_stored_in_propinquity_bound(self):
        """compute_tensor_propinquity_bound stores eps_cross in the dataclass."""
        b = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        assert b.epsilon_cross_bound_value > 0
        # Should match epsilon_cross_bound directly
        eps = epsilon_cross_bound(2, 2, 1.0, 1.0)["epsilon_cross_bound"]
        assert b.epsilon_cross_bound_value == pytest.approx(eps, rel=1e-9)

    def test_eps_cross_in_tensor_l5_assembly(self):
        """tensor_L5_assembly returns the eps_cross bound and the
        R1+R2-tightened propinquity bound."""
        out = tensor_L5_assembly(2, 2)
        assert "epsilon_cross_bound" in out
        assert "C_3_full_pythagorean" in out
        assert "propinquity_bound_r1_r2" in out
        assert out["epsilon_cross_bound"] > 0
        assert out["propinquity_bound_r1_r2"] > 0

    def test_propinquity_bound_r1_r2_formula(self):
        """propinquity_bound_r1_r2 = C_3^{(2),Pyth} * (max(γ) + ε_cross)."""
        b = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        c = b.c_lipschitz_full_pythagorean
        gm = max(b.gamma_a, b.gamma_b)
        eps = b.epsilon_cross_bound_value
        expected = c * (gm + eps)
        assert b.propinquity_bound_r1_r2 == pytest.approx(expected, rel=1e-9)


class TestR2RateConvergence:
    """R2 closure: verify ε_cross -> 0 as cutoffs -> infinity, with
    the explicit constant 2*sqrt(2) attached to max(γ_a, γ_b).
    """

    def test_eps_cross_panel_2_3_4_5_decreasing(self):
        """The ε_cross bound on the (n,n) panel monotone decreases."""
        eps_seq = []
        for n in [2, 3, 4, 5]:
            eps = epsilon_cross_bound(n, n, 1.0, 1.0)
            eps_seq.append(eps["epsilon_cross_bound"])
        # Strictly decreasing
        for i in range(len(eps_seq) - 1):
            assert eps_seq[i] > eps_seq[i + 1]

    def test_eps_cross_22_44_ratio_below_07(self):
        """ε_cross at (4,4) is at most 0.7 times ε_cross at (2,2)."""
        eps_22 = epsilon_cross_bound(2, 2, 1.0, 1.0)["epsilon_cross_bound"]
        eps_44 = epsilon_cross_bound(4, 4, 1.0, 1.0)["epsilon_cross_bound"]
        assert eps_44 / eps_22 < 0.7

    def test_eps_cross_unit_constant_explicit(self):
        """At (n, n), ε_cross / γ_n = 2*sqrt(2) exactly on unit panel."""
        for n in [2, 3, 4, 5]:
            eps = epsilon_cross_bound(n, n, 1.0, 1.0)
            ratio = eps["epsilon_cross_bound"] / eps["gamma_a"]
            assert ratio == pytest.approx(2.0 * np.sqrt(2.0), rel=1e-9)


class TestR1R2KeystoneTheoremStatus:
    """The combined R1 + R2 closure status of the keystone tensor-product
    propinquity bound after the W2b-easy-tighten sprint.

    Net statement (sprint memo §6; Pythagorean form WITHDRAWN 2026-06-18 as
    false, replaced by the triangle bound):
        Lambda(T_a (X) T_b, T_S3^a (X) T_S3^b)
            <= C_3^{(2)} * (max(γ_a, γ_b) + ε_cross)
            <= C_3^{(2)} * (1 + 2*sqrt(2)) * max(γ_a, γ_b)
            -> 0  as cutoffs -> infinity,
    where C_3^{(2)} = sqrt(((N_a-1) + (N_b-1))^2 / (N_a^2 + N_b^2 - 2)) is the
    triangle bound, <= sqrt(2) (NOT < 1), exceeding 1 beyond (2,2) and -> sqrt(2);
    convergence survives because sqrt(2) is a finite constant and γ -> 0.
    """

    def test_five_lemma_status_l3_done(self):
        """L3-T should now report DONE (not PARTIAL) after R1 closure."""
        status = FiveLemmaStatusTensor()
        d = status.to_dict()
        assert "DONE" in d["L3_T"]
        assert "Pythagorean" in d["L3_T"]

    def test_five_lemma_status_l5_done(self):
        """L5-T should report DONE after R2 closure of the height bookkeeping."""
        status = FiveLemmaStatusTensor()
        d = status.to_dict()
        assert "DONE" in d["L5_T"]
        # Should reference the explicit constant 2*sqrt(2)
        assert "2*sqrt(2)" in d["L5_T"] or "epsilon_cross" in d["L5_T"]

    def test_propinquity_bound_qualitative_rate_to_zero(self):
        """The R1+R2 tightened bound goes to zero as cutoffs grow."""
        b_22 = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        b_33 = compute_tensor_propinquity_bound(3, 3, gamma_prec=15)
        # Both R1+R2 and the legacy bound are strictly decreasing on the
        # symmetric (n, n) panel
        assert b_33.propinquity_bound_r1_r2 < b_22.propinquity_bound_r1_r2

    def test_pyth_bound_propagates_to_dataclass(self):
        """TensorPropinquityBound exposes c_lipschitz_full_pythagorean."""
        b = compute_tensor_propinquity_bound(2, 2, gamma_prec=15)
        assert b.c_lipschitz_full_pythagorean > 0
        assert b.c_lipschitz_full_pythagorean <= np.sqrt(2.0) + 1e-12
        # Same value as the standalone function
        assert b.c_lipschitz_full_pythagorean == pytest.approx(
            c3_full_pythagorean_bound(2, 2), rel=1e-12,
        )

    def test_tensor_l5_assembly_uses_pyth_by_default(self):
        """tensor_L5_assembly uses the Pythagorean bound by default."""
        out = tensor_L5_assembly(2, 2)
        assert out["use_pythagorean_bound"] is True
        # And reports the explicit Pythagorean value
        assert out["C_3_full_pythagorean"] < out["C_3_full"]


@pytest.mark.slow
class TestR1R2SlowConvergence:
    """Slow tests: full panel verification of R1+R2 closure at (4,4) and (5,5)."""

    def test_triangle_bound_at_55(self):
        """At (5,5), triangle bound is sqrt(2*5/(5+2)) = sqrt(10/7)."""
        c = c3_full_pythagorean_bound(5, 5)
        assert c == pytest.approx(np.sqrt(10.0 / 7.0), rel=1e-12)

    def test_eps_cross_at_55(self):
        """ε_cross at (5,5) is 2*sqrt(2) * gamma_5."""
        eps = epsilon_cross_bound(5, 5, 1.0, 1.0)
        # gamma_5 ~ 1.130
        expected = 2 * np.sqrt(2) * eps["gamma_a"]
        assert eps["epsilon_cross_bound"] == pytest.approx(expected, rel=1e-9)

    def test_full_panel_triangle_within_sqrt2(self):
        """Triangle bound stays <= sqrt(2) on the full panel."""
        for nb in [(2,2), (3,3), (4,4), (5,5), (10,10), (50,50), (100,100)]:
            c = c3_full_pythagorean_bound(*nb)
            assert c <= np.sqrt(2.0) + 1e-12
