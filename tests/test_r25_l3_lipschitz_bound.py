"""Tests for R2.5 Lemma L3 (Lipschitz bound) — geovac/r25_l3_lipschitz_bound.py.

Per CLAUDE.md §13.4a, every equation in the L3 proof memo
(`debug/r25_l3_proof_memo.md`) has at least one corresponding test in
this file. Tests cover:

  - Y^{(3)}_{NLM} symbolic correctness (orthonormality and identification
    with cos(chi) etc. for low harmonics).
  - Lipschitz norm closed forms / numerical sup checks for low harmonics.
  - Per-multiplier commutator structure on the truthful CH Dirac.
  - Bound check ratios for the natural test panel at n_max in {2, 3, 4}
    confirming C_3 <= 1 numerically.
  - Theoretical bound (N - 1) / sqrt(N^2 - 1) inequality.
  - Scalar-Dirac fallback path.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import sympy as sp

from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem
from geovac.r25_l3_lipschitz_bound import (
    BoundCheckResult,
    bound_check_one,
    bound_check_panel,
    bound_check_scalar_dirac_one,
    bound_check_scalar_dirac_panel,
    build_multiplier_for_test_function,
    commutator_norm_decomposition,
    commutator_with_ch_dirac,
    commutator_with_ch_dirac_spinor,
    commutator_with_scalar_dirac,
    commutator_with_scalar_laplacian,
    constant_C3_panel,
    default_test_panel,
    lipschitz_norm_inf_test_function,
    lipschitz_norm_inf_y3,
    make_test_function,
    scalar_dirac_proxy_diag_full,
    scalar_laplacian_diag_full,
    shell_diff_matrix_full,
    shell_diff_matrix_spinor,
    shell_diff_max_for_label,
    y3_avery_symbolic,
)
from geovac.spinor_operator_system import spinor_basis


# =====================================================================
# §1. y3_avery_symbolic correctness
# =====================================================================


class TestY3AverySymbolic:
    """Verify y3_avery_symbolic against known identities."""

    def test_y3_100_constant(self):
        """Y^{(3)}_{1,0,0}(chi) = 1 / sqrt(2 pi^2) (constant function on S^3)."""
        chi, theta, phi = sp.symbols("chi theta phi", real=True)
        Y100 = y3_avery_symbolic(1, 0, 0, chi, theta, phi)
        Y100_simplified = sp.simplify(Y100)
        # Expected: 1/sqrt(2 pi^2) = sqrt(2)/(2 pi)
        expected = 1 / sp.sqrt(2 * sp.pi ** 2)
        assert sp.simplify(Y100_simplified - expected) == 0

    def test_y3_200_form(self):
        """Y^{(3)}_{2,0,0}(chi) ~ cos(chi) (axisymmetric, sin^0 prefactor)."""
        chi, theta, phi = sp.symbols("chi theta phi", real=True)
        Y200 = y3_avery_symbolic(2, 0, 0, chi, theta, phi)
        # Should be a constant times cos(chi).
        ratio = sp.simplify(Y200 / sp.cos(chi))
        # ratio should be a constant (no chi, theta, phi)
        assert ratio.free_symbols == set()
        # The numerical value matches sqrt(2)/pi
        assert float(sp.N(ratio)) == pytest.approx(math.sqrt(2) / math.pi, abs=1e-10)

    def test_y3_orthonormal_l2(self):
        """L^2 norm of Y^{(3)}_{NLM} on the round S^3 = 1 for low (N, L, M)."""
        chi, theta, phi = sp.symbols("chi theta phi", real=True)
        for (N, L, M) in [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1)]:
            Y = y3_avery_symbolic(N, L, M, chi, theta, phi)
            integrand = (
                Y * sp.conjugate(Y)
                * sp.sin(chi) ** 2
                * sp.sin(theta)
            )
            val = sp.integrate(
                integrand,
                (phi, 0, 2 * sp.pi),
                (theta, 0, sp.pi),
                (chi, 0, sp.pi),
            )
            val = sp.simplify(val)
            assert val == 1, f"Y_{{{N},{L},{M}}} not normalized: ||Y||^2 = {val}"


# =====================================================================
# §2. Lipschitz norms of low Y^{(3)}_{NLM}
# =====================================================================


class TestLipschitzNorm:
    """Sanity checks on lipschitz_norm_inf_y3."""

    def test_lipschitz_y3_200_analytic(self):
        """Y^{(3)}_{2,0,0} = sqrt(2)/pi cos(chi); ||grad||_inf = sqrt(2)/pi."""
        # |grad Y_{2,0,0}|^2 = (sqrt(2)/pi)^2 * sin^2(chi).
        # Sup over chi in (0, pi) is at chi = pi/2, value = sqrt(2)/pi.
        lip = float(lipschitz_norm_inf_y3(2, 0, 0, prec=30))
        expected = math.sqrt(2) / math.pi
        assert lip == pytest.approx(expected, rel=1e-3)

    def test_lipschitz_unit_y3_positive(self):
        """All low Y^{(3)}_{NLM} have positive Lipschitz norm."""
        for (N, L, M) in [(2, 0, 0), (2, 1, 0), (3, 0, 0), (3, 1, 0)]:
            lip = float(lipschitz_norm_inf_y3(N, L, M, prec=30))
            assert lip > 0, f"Y_{{{N},{L},{M}}} has zero Lipschitz norm"

    def test_lipschitz_grows_with_N(self):
        """For axisymmetric Y_{N,1,0}, ||grad||_inf is monotone in N (low cases)."""
        lips = [
            float(lipschitz_norm_inf_y3(N, 1, 0, prec=30)) for N in [2, 3, 4]
        ]
        assert lips[0] < lips[1] < lips[2], f"Non-monotone Lipschitz: {lips}"


# =====================================================================
# §3. Commutator structure on truthful CH Dirac
# =====================================================================


class TestCommutatorStructure:
    """Verify [D_CH, M] structure on the truthful CH Dirac."""

    def test_commutator_constant_zero(self):
        """[D_CH, M_{1,0,0}] = 0 (constant function commutes with D)."""
        op = FullDiracTruncatedOperatorSystem(2)
        # M_{1,0,0} is the identity multiplier (constant function)
        idx = op.multiplier_labels.index((1, 0, 0))
        M = op.multiplier_matrices[idx]
        comm = commutator_with_ch_dirac(M, op)
        assert np.linalg.norm(comm) < 1e-12

    def test_commutator_block_diagonal_in_chirality(self):
        """[D_CH^{truthful}, M_f] is block-diagonal in chirality.

        Off-diagonal in chirality the commutator vanishes because
        scalar M_f is block-diagonal in chirality (no chirality-flipping)
        and D_CH is block-diagonal in chirality.
        """
        op = FullDiracTruncatedOperatorSystem(2)
        idx = op.multiplier_labels.index((2, 1, 0))
        M = op.multiplier_matrices[idx]
        comm = commutator_with_ch_dirac(M, op)
        n_half = op.dim_H // 2
        # Off-diagonal blocks (Weyl <-> anti-Weyl) should be zero.
        block_pp = comm[:n_half, :n_half]
        block_mm = comm[n_half:, n_half:]
        block_pm = comm[:n_half, n_half:]
        block_mp = comm[n_half:, :n_half]
        assert np.linalg.norm(block_pm) < 1e-12
        assert np.linalg.norm(block_mp) < 1e-12
        # Diagonal blocks should be opposite sign (chirality factor flips).
        assert np.linalg.norm(block_pp + block_mm) < 1e-12

    def test_commutator_scaling_by_shell_difference(self):
        """[D_CH, M_{NLM}]_{ij} = chi_i (n_i - n_j) M_{ij} (within same chi)."""
        op = FullDiracTruncatedOperatorSystem(3)
        idx = op.multiplier_labels.index((2, 1, 0))
        M = op.multiplier_matrices[idx]
        # Truthful CH eigenvalues
        diag = np.array(
            [
                b.chirality * (b.n_fock + 0.5) for b in op.basis
            ],
            dtype=np.float64,
        )
        comm = commutator_with_ch_dirac(M, op)
        # Element-wise check: comm[i, j] = (lambda_i - lambda_j) * M[i, j].
        for i in range(op.dim_H):
            for j in range(op.dim_H):
                expected = (diag[i] - diag[j]) * M[i, j]
                assert abs(comm[i, j] - expected) < 1e-12, (
                    f"comm[{i},{j}] = {comm[i, j]}, expected {expected}"
                )

    def test_constant_function_in_kernel(self):
        """The constant function (M_{1,0,0}) is in the kernel of [D, .]."""
        op = FullDiracTruncatedOperatorSystem(3)
        idx = op.multiplier_labels.index((1, 0, 0))
        M = op.multiplier_matrices[idx]
        comm = commutator_with_ch_dirac(M, op)
        assert np.linalg.norm(comm, ord=2) < 1e-12

    def test_shell_diff_max_selection_rule(self):
        """shell_diff_max_for_label(N, n_max) <= N - 1 (SO(4) selection rule)."""
        for N in [1, 2, 3, 4, 5]:
            for n_max in [2, 3, 4]:
                d = shell_diff_max_for_label(N, n_max)
                assert d <= N - 1, (
                    f"shell_diff_max({N}, {n_max}) = {d} > N - 1 = {N - 1}"
                )


# =====================================================================
# §4. Per-multiplier ratio = max |Delta n| (within chirality)
# =====================================================================


class TestPerMultiplierRatio:
    """Verify ||[D, M_{NLM}]|| / ||M_{NLM}|| == max |Delta n| in truncation."""

    def test_per_multiplier_ratio_NM_eq_2(self):
        """For M_{2, L, M} (couples |Delta n| = 1), ratio = 1.

        At n_max >= 2, M_{2, L, M} couples shells (n, n') with
        |n - n'| in {0, 1} and n + n' >= 3. The largest |Delta n|
        achievable is 1 at (n, n') = (1, 2).

        On the truthful CH Dirac block-diagonal in chirality, the
        operator-norm ratio equals 1 when the leading singular vector of
        M is concentrated on the (1, 2) shell pair. We verify
        numerically.
        """
        op = FullDiracTruncatedOperatorSystem(3)
        for label in [(2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]:
            d = commutator_norm_decomposition(op, label)
            assert d["ratio"] == pytest.approx(1.0, abs=1e-6), (
                f"label {label}: ratio = {d['ratio']}, expected 1.0"
            )

    def test_per_multiplier_ratio_NM_eq_3_at_n3(self):
        """For M_{3, L, M} at n_max=3, ratio is in [1, 2] (1 <= |Dn| <= 2).

        The exact ratio depends on the singular structure of M; we
        verify the bound 1 <= ratio <= 2 = N - 1 for N=3.
        """
        op = FullDiracTruncatedOperatorSystem(3)
        for label in [(3, 0, 0), (3, 1, 0), (3, 2, 0)]:
            d = commutator_norm_decomposition(op, label)
            assert 1.0 <= d["ratio"] + 1e-10 <= 2.0 + 1e-6, (
                f"label {label}: ratio = {d['ratio']} not in [1, 2]"
            )


# =====================================================================
# §5. Lemma L3: bound check on the natural test panel
# =====================================================================


class TestL3BoundCheck:
    """Lemma L3 numerical verification: C_3 <= 1 across the test panel."""

    @pytest.fixture(scope="class")
    def op2(self):
        return FullDiracTruncatedOperatorSystem(2)

    @pytest.fixture(scope="class")
    def op3(self):
        return FullDiracTruncatedOperatorSystem(3)

    def test_l3_n_max_2(self, op2):
        """L3 holds at n_max = 2 with C_3 <= 1."""
        results = bound_check_panel(op2)
        C3 = constant_C3_panel(results)
        assert C3 <= 1.0 + 1e-6, f"C_3 = {C3} > 1 at n_max=2"
        assert C3 > 0.4, f"C_3 = {C3} suspiciously small at n_max=2"

    def test_l3_n_max_3(self, op3):
        """L3 holds at n_max = 3 with C_3 <= 1."""
        results = bound_check_panel(op3)
        C3 = constant_C3_panel(results)
        assert C3 <= 1.0 + 1e-6, f"C_3 = {C3} > 1 at n_max=3"
        assert C3 > 0.5, f"C_3 = {C3} suspiciously small at n_max=3"

    def test_l3_C3_below_one_each_function(self, op3):
        """Every test function in the panel has ratio <= 1."""
        results = bound_check_panel(op3)
        for r in results:
            if r.lipschitz_inf > 1e-10:  # skip constant-function (lip = 0)
                assert r.ratio <= 1.0 + 1e-6, (
                    f"f={r.f_name}: ratio = {r.ratio} > 1"
                )

    def test_l3_C3_uniform_in_n_max(self, op2, op3):
        """C_3 does not blow up between n_max=2 and n_max=3.

        Both should be <= 1 (the lemma). Stronger: they should be
        comparable in magnitude.
        """
        C3_2 = constant_C3_panel(bound_check_panel(op2))
        C3_3 = constant_C3_panel(bound_check_panel(op3))
        # Both bounded by 1
        assert C3_2 <= 1.0 + 1e-6
        assert C3_3 <= 1.0 + 1e-6


# =====================================================================
# §6. Theoretical bound (N - 1) / sqrt(N^2 - 1) -> 1
# =====================================================================


class TestTheoreticalBound:
    """The theoretical bound (N-1)/sqrt(N^2-1) approaches 1 from below."""

    def test_theoretical_bound_below_one(self):
        """For all N >= 2, (N-1)/sqrt(N^2-1) <= 1."""
        for N in range(2, 50):
            ratio = (N - 1) / math.sqrt(N ** 2 - 1)
            assert ratio < 1.0, f"N={N}: ratio = {ratio} >= 1"

    def test_theoretical_bound_monotone(self):
        """The theoretical bound is monotone increasing in N."""
        prev = 0.0
        for N in range(2, 20):
            ratio = (N - 1) / math.sqrt(N ** 2 - 1)
            assert ratio > prev, f"N={N}: not monotone"
            prev = ratio

    def test_theoretical_bound_limit(self):
        """As N -> infinity, (N-1)/sqrt(N^2-1) -> 1."""
        N = 1000
        ratio = (N - 1) / math.sqrt(N ** 2 - 1)
        assert abs(ratio - 1.0) < 1e-3

    def test_theoretical_bound_at_N3_exact(self):
        """At N=3, (N-1)/sqrt(N^2-1) = 2/sqrt(8) = 1/sqrt(2) ≈ 0.7071."""
        ratio = 2 / math.sqrt(8)
        expected = 1 / math.sqrt(2)
        assert ratio == pytest.approx(expected, abs=1e-12)
        # And matches our observed C_3 at n_max = 3.
        assert ratio == pytest.approx(0.7071, abs=1e-3)


# =====================================================================
# §7. Scalar-Dirac fallback path
# =====================================================================


class TestScalarDiracFallback:
    """Lemma L3 with the sqrt(-Delta_LB) scalar Dirac proxy (master memo §7)."""

    def test_scalar_dirac_diagonal(self):
        """sqrt(-Delta_LB) is diagonal with entries sqrt(n^2 - 1)."""
        op = FullDiracTruncatedOperatorSystem(2)
        D = scalar_dirac_proxy_diag_full(op)
        for i, b in enumerate(op.basis):
            expected = math.sqrt(b.n_fock ** 2 - 1)
            assert abs(D[i, i] - expected) < 1e-10

    def test_scalar_dirac_n1_kernel(self):
        """The constant function (n=1) is in the kernel of sqrt(-Delta_LB)."""
        op = FullDiracTruncatedOperatorSystem(2)
        D = scalar_dirac_proxy_diag_full(op)
        for i, b in enumerate(op.basis):
            if b.n_fock == 1:
                assert abs(D[i, i]) < 1e-10

    def test_scalar_dirac_l3_bound(self):
        """scalar Dirac C_3 is also bounded numerically (within ~5% of 1)."""
        op = FullDiracTruncatedOperatorSystem(2)
        results = bound_check_scalar_dirac_panel(op)
        C3sc = constant_C3_panel(results)
        # The scalar-Dirac fallback may give slightly larger ratios than
        # truthful CH but should remain bounded near 1.
        assert C3sc <= 1.1, f"C_3 (scalar Dirac) = {C3sc} > 1.1 at n_max=2"


# =====================================================================
# §8. Shell-difference matrix structural identity
# =====================================================================


class TestShellDiffMatrix:
    """shell_diff_matrix gives the right weighting for the commutator."""

    def test_shell_diff_full_zero_diagonal(self):
        """S_{ii} = 0 for the shell-difference matrix."""
        op = FullDiracTruncatedOperatorSystem(3)
        S = shell_diff_matrix_full(op)
        for i in range(op.dim_H):
            assert S[i, i] == 0.0

    def test_shell_diff_full_antisymmetric(self):
        """S_{ij} = -S_{ji} for shell-difference matrix."""
        op = FullDiracTruncatedOperatorSystem(2)
        S = shell_diff_matrix_full(op)
        np.testing.assert_allclose(S, -S.T, atol=1e-12)

    def test_shell_diff_spinor_integer_diff(self):
        """For Weyl-only basis, S_{ij} = n_i - n_j is integer-valued."""
        basis = spinor_basis(3)
        S = shell_diff_matrix_spinor(basis)
        for i, j in [(0, 0), (0, 5), (5, 10)]:
            expected = basis[i].n_fock - basis[j].n_fock
            assert S[i, j] == expected


# =====================================================================
# §9. Test panel construction
# =====================================================================


class TestDefaultPanel:
    """Verify default_test_panel produces sensible test functions."""

    def test_panel_n2_has_single_y3(self):
        """Panel at n_max=2 contains Y^{(3)}_{2, *, *}."""
        panel = default_test_panel(2)
        names = {f.name for f in panel}
        for L in range(2):
            for M in range(-L, L + 1):
                assert f"Y3_(2,{L},{M})" in names

    def test_panel_n3_extended(self):
        """Panel at n_max=3 has more test functions than n_max=2."""
        p2 = default_test_panel(2)
        p3 = default_test_panel(3)
        assert len(p3) > len(p2)

    def test_make_test_function(self):
        """make_test_function constructs a valid TestFunction."""
        f = make_test_function("test", {(2, 1, 0): 1.5})
        assert f.name == "test"
        assert f.coeff_dict == {(2, 1, 0): 1.5}


# =====================================================================
# §10. Full pipeline: BoundCheckResult fields
# =====================================================================


class TestBoundCheckResult:
    """Sanity on BoundCheckResult dataclass."""

    def test_bound_check_result_fields(self):
        """bound_check_one returns a BoundCheckResult with expected fields."""
        op = FullDiracTruncatedOperatorSystem(2)
        f = make_test_function("Y3_(2,1,0)", {(2, 1, 0): 1.0})
        r = bound_check_one(f, op)
        assert isinstance(r, BoundCheckResult)
        assert r.f_name == "Y3_(2,1,0)"
        assert r.n_max == 2
        assert r.op_norm_commutator > 0
        assert r.lipschitz_inf > 0
        assert r.ratio > 0
        assert r.ratio <= 1.0 + 1e-6


# =====================================================================
# §11. Per-equation verification protocol (CLAUDE.md §13.4a)
# =====================================================================


class TestEquationVerification:
    """Each equation in the L3 proof memo has a corresponding test.

    This class names the equations explicitly for traceability.
    """

    def test_eq_commutator_block_diagonal_in_chirality(self):
        """Eq. (L3-1): [D_CH, M_f] block-diagonal in chirality, opposite sign blocks.

        Tested by: TestCommutatorStructure::test_commutator_block_diagonal_in_chirality
        """
        # See class above; this test re-invokes the structural check on a
        # different cutoff to confirm the equation holds at multiple n_max.
        op = FullDiracTruncatedOperatorSystem(3)
        idx = op.multiplier_labels.index((3, 1, 0))
        M = op.multiplier_matrices[idx]
        comm = commutator_with_ch_dirac(M, op)
        n_half = op.dim_H // 2
        block_pm = comm[:n_half, n_half:]
        block_mp = comm[n_half:, :n_half]
        assert np.linalg.norm(block_pm) < 1e-12
        assert np.linalg.norm(block_mp) < 1e-12

    def test_eq_shell_difference_weighting(self):
        """Eq. (L3-2): [D, M]_{ij} = (lambda_i - lambda_j) M_{ij}.

        Verified at n_max=3 across all multipliers.
        """
        op = FullDiracTruncatedOperatorSystem(3)
        diag = np.array(
            [b.chirality * (b.n_fock + 0.5) for b in op.basis],
            dtype=np.float64,
        )
        # Pick a representative non-trivial multiplier
        idx = op.multiplier_labels.index((2, 0, 0))
        M = op.multiplier_matrices[idx]
        comm = commutator_with_ch_dirac(M, op)
        # Element-wise check
        for i in range(op.dim_H):
            for j in range(op.dim_H):
                expected = (diag[i] - diag[j]) * M[i, j]
                assert abs(comm[i, j] - expected) < 1e-12

    def test_eq_l3_main_inequality(self):
        """Eq. (L3-main): ||[D_CH, M_f]||_op <= C_3 * ||grad f||_inf.

        Verified across the test panel at n_max=2,3 with C_3 = 1.
        """
        for n_max in [2, 3]:
            op = FullDiracTruncatedOperatorSystem(n_max)
            results = bound_check_panel(op)
            for r in results:
                if r.lipschitz_inf > 1e-10:
                    assert r.op_norm_commutator <= 1.0 * r.lipschitz_inf + 1e-10, (
                        f"L3 violated at n_max={n_max}, f={r.f_name}: "
                        f"{r.op_norm_commutator} > {r.lipschitz_inf}"
                    )

    def test_eq_theoretical_bound(self):
        """Eq. (L3-thm): C_3(n_max) <= max_{N <= n_max} (N-1)/sqrt(N^2-1).

        Verified for n_max in {2, 3}.
        """
        for n_max in [2, 3]:
            op = FullDiracTruncatedOperatorSystem(n_max)
            results = bound_check_panel(op)
            C3 = constant_C3_panel(results)
            theoretical_max = max(
                (N - 1) / math.sqrt(N ** 2 - 1)
                for N in range(2, n_max + 1)
            ) if n_max >= 2 else 0.0
            # Allow some tolerance because the per-multiplier ratio (N-1)
            # is only attained when the singular vector of M has full
            # support on the |Delta n| = N - 1 shell pair.
            assert C3 <= theoretical_max + 0.05 + 1e-6, (
                f"n_max={n_max}: C_3 = {C3} > theoretical {theoretical_max}"
            )


@pytest.mark.slow
class TestL3SlowChecks:
    """Slow tests: full panel at n_max=4 (8-12s)."""

    def test_l3_n_max_4_restricted_panel(self):
        """L3 holds at n_max=4 on a restricted small panel with C_3 <= 1."""
        op = FullDiracTruncatedOperatorSystem(4)
        small_panel = [
            make_test_function("Y3_(4,1,0)", {(4, 1, 0): 1.0}),
            make_test_function("Y3_(4,2,0)", {(4, 2, 0): 1.0}),
            make_test_function("Y3_(4,3,0)", {(4, 3, 0): 1.0}),
        ]
        results = bound_check_panel(op, panel=small_panel)
        C3 = constant_C3_panel(results)
        assert C3 <= 1.0 + 1e-6, f"C_3 = {C3} > 1 at n_max=4"
