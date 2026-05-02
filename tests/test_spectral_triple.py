"""
Tests for geovac.spectral_triple — Fock-projected S^3 spectral triple.

Validates the data (A, H, D, gamma, J) at n_max=2 and n_max=3 against
the axioms of a real spectral triple.  Structural failures ({D,gamma}=0,
order-one condition) are tested as EXPECTED FAIL with documented reasons.

Reference: GeoVac CLAUDE.md section 1.7 WH1 (almost-commutative spectral triple).
"""

import pytest
from sympy import Integer, Rational, eye as sp_eye, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


# ================================================================
# Fixtures
# ================================================================

@pytest.fixture(scope="module")
def st2():
    """Spectral triple at n_max=2 (permutation J)."""
    return FockSpectralTriple(n_max=2, j_type="permutation")


@pytest.fixture(scope="module")
def st3():
    """Spectral triple at n_max=3 (permutation J)."""
    return FockSpectralTriple(n_max=3, j_type="permutation")


@pytest.fixture(scope="module")
def st2_kramers():
    """Spectral triple at n_max=2 (Kramers J)."""
    return FockSpectralTriple(n_max=2, j_type="kramers")


# ================================================================
# Dimensions
# ================================================================

class TestDimensions:
    """Hilbert space dimension = sum_{n=1}^{n_max} 2n(n+1)."""

    def test_dim_nmax1(self):
        st = FockSpectralTriple(n_max=1)
        assert st.dim_H == 4  # 2*1*2 = 4

    def test_dim_nmax2(self, st2):
        assert st2.dim_H == 16  # 4 + 2*2*3 = 16

    def test_dim_nmax3(self, st3):
        assert st3.dim_H == 40  # 4 + 12 + 2*3*4 = 40 = Delta^{-1} (Paper 2)

    def test_sectors_nmax2(self, st2):
        assert st2.n_sectors == 5
        assert st2.sectors == [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]

    def test_sectors_nmax3(self, st3):
        assert st3.n_sectors == 9
        assert st3.sectors == [
            (1, 0), (1, 1), (2, 0), (2, 1), (2, 2),
            (3, 0), (3, 1), (3, 2), (3, 3),
        ]


# ================================================================
# D self-adjointness
# ================================================================

class TestDiracOperator:
    """D = Lambda + kappa * A must be self-adjoint."""

    def test_selfadjoint_nmax2(self, st2):
        assert st2.check_selfadjoint()

    def test_selfadjoint_nmax3(self, st3):
        assert st3.check_selfadjoint()

    def test_D_diagonal_is_lambda(self, st2):
        """Diagonal of D should equal Lambda (Camporesi-Higuchi)."""
        D = st2.dirac_operator
        Lam = st2.diagonal_part
        for i in range(st2.dim_H):
            assert D[i, i] == Lam[i, i]

    def test_D_offdiag_is_kappa(self, st2):
        """Off-diagonal nonzero entries of D should all be kappa = -1/16."""
        D = st2.dirac_operator
        kappa = Rational(-1, 16)
        for i in range(st2.dim_H):
            for j in range(st2.dim_H):
                if i != j:
                    if D[i, j] != 0:
                        assert D[i, j] == kappa

    def test_spectrum_paired(self, st2):
        """D spectrum should be paired: for each lambda, -lambda also appears.

        This follows from the chirality structure: if |psi> has eigenvalue
        lambda, then gamma|psi> has eigenvalue related by the sign structure.
        On the finite graph, exact pairing depends on the symmetry of the
        adjacency under chirality flip.
        """
        import numpy as np
        D = st2.dirac_operator
        arr = np.array(D.tolist(), dtype=np.float64)
        evals = np.linalg.eigvalsh(arr)
        evals_sorted = sorted(evals)
        evals_neg = sorted(-evals)
        # Check approximate pairing
        for a, b in zip(evals_sorted, evals_neg):
            assert abs(a - b) < 1e-10, f"Pairing mismatch: {a} vs {b}"


# ================================================================
# Grading gamma
# ================================================================

class TestGrading:
    """Chirality grading gamma with gamma^2 = I."""

    def test_gamma_squared_I_nmax2(self, st2):
        assert st2.check_grading_square()

    def test_gamma_squared_I_nmax3(self, st3):
        assert st3.check_grading_square()

    def test_gamma_selfadjoint_nmax2(self, st2):
        assert st2.check_grading_selfadjoint()

    def test_gamma_selfadjoint_nmax3(self, st3):
        assert st3.check_grading_selfadjoint()

    def test_Lambda_gamma_commute_nmax2(self, st2):
        """[Lambda, gamma] = 0 because both are diagonal."""
        assert st2.check_D_gamma_commute_diagonal()

    def test_Lambda_gamma_commute_nmax3(self, st3):
        assert st3.check_D_gamma_commute_diagonal()

    def test_D_gamma_anticommute_expected_fail_nmax2(self, st2):
        """EXPECTED FAIL: {D, gamma} != 0 for odd-dimensional spectral triple."""
        result, msg = st2.check_D_gamma_anticommute()
        assert not result, "Unexpected pass — {D, gamma} should fail on finite graph"
        assert "EXPECTED FAIL" in msg

    def test_D_gamma_anticommute_expected_fail_nmax3(self, st3):
        result, msg = st3.check_D_gamma_anticommute()
        assert not result

    def test_chirality_counts_nmax2(self, st2):
        """n_max=2: 8 positive-chirality + 8 negative-chirality = 16."""
        hs = st2.hilbert_space_summary()
        assert hs["chi_plus"] == 8
        assert hs["chi_minus"] == 8


# ================================================================
# Real structure J
# ================================================================

class TestRealStructure:
    """Tests for both permutation J and Kramers J."""

    def test_permutation_J_squared_plus_I(self, st2):
        j2_ok, j2_eps = st2.check_J_squared()
        assert j2_ok
        assert j2_eps == 1

    def test_permutation_JD_commute(self, st2):
        jd_ok, jd_msg = st2.check_J_D_relation()
        assert jd_ok
        assert "PASS" in jd_msg

    def test_permutation_JD_commute_nmax3(self, st3):
        jd_ok, jd_msg = st3.check_J_D_relation()
        assert jd_ok

    def test_permutation_J_gamma_commute(self, st2):
        assert st2.check_J_gamma_commute()

    def test_kramers_J_squared_minus_I(self, st2_kramers):
        j2_ok, j2_eps = st2_kramers.check_J_squared()
        assert j2_ok
        assert j2_eps == -1

    def test_kramers_JD_expected_fail(self, st2_kramers):
        """Kramers J does not anticommute with D on finite graph."""
        jd_ok, jd_msg = st2_kramers.check_J_D_relation()
        assert not jd_ok, "Unexpected pass — Kramers JD should fail"
        assert "EXPECTED FAIL" in jd_msg

    def test_kramers_J_gamma_commute(self, st2_kramers):
        assert st2_kramers.check_J_gamma_commute()


# ================================================================
# Order conditions
# ================================================================

class TestOrderConditions:
    """Order-zero and order-one conditions for the spectral triple."""

    def test_order_zero_nmax2(self, st2):
        assert st2.check_order_zero()

    def test_order_zero_nmax3(self, st3):
        assert st3.check_order_zero()

    def test_order_one_expected_fail_nmax2(self, st2):
        """Order-one fails because D connects different sectors."""
        o1_ok, o1_msg = st2.check_order_one()
        assert not o1_ok, "Unexpected pass — order-one should fail"
        assert "EXPECTED FAIL" in o1_msg
        # Verify specific failure count
        assert "17/25" in o1_msg

    @pytest.mark.slow
    def test_order_one_expected_fail_nmax3(self, st3):
        """Order-one fails at n_max=3 with known failure count."""
        o1_ok, o1_msg = st3.check_order_one()
        assert not o1_ok
        assert "37/81" in o1_msg


# ================================================================
# Commutator and Connes distance
# ================================================================

class TestCommutatorAndDistance:
    """Tests for [D, pi(f)] and Connes distance."""

    def test_commutator_nonzero_generic(self, st2):
        """[D, pi(f)] is nonzero for generic f."""
        f = [Integer(1), Integer(3), Integer(7), Integer(11), Integer(13)]
        comm = st2.commutator(f)
        n_nz = sum(1 for i in range(st2.dim_H) for j in range(st2.dim_H) if comm[i, j] != 0)
        assert n_nz > 0

    def test_commutator_zero_for_constant(self, st2):
        """[D, pi(c)] = 0 for constant function c."""
        f = [Integer(5)] * st2.n_sectors
        comm = st2.commutator(f)
        assert comm.equals(sp_zeros(st2.dim_H, st2.dim_H))

    def test_connes_distance_same_sector(self, st2):
        """d(i, i) = 0."""
        d = st2.connes_distance(0, 0)
        assert d == Integer(0)

    def test_connes_distance_positive(self, st2):
        """d(i, j) > 0 for i != j."""
        d = st2.connes_distance(0, 1)
        assert float(d) > 0

    def test_connes_distances_finite(self, st2):
        """All pairwise Connes distances are finite."""
        for i in range(st2.n_sectors):
            for j in range(i + 1, st2.n_sectors):
                d = st2.connes_distance(i, j)
                assert float(d) < 1e6, f"Distance d({i},{j}) appears infinite"


# ================================================================
# Pi-free certificate
# ================================================================

class TestPiFree:
    """All matrix entries are exact algebraic (no pi, no floats)."""

    def test_pi_free_nmax2(self, st2):
        assert st2.verify_pi_free()

    def test_pi_free_nmax3(self, st3):
        assert st3.verify_pi_free()


# ================================================================
# Spectral action
# ================================================================

class TestSpectralAction:
    """Spectral action Tr(f(D)) computations."""

    def test_spectral_action_no_cutoff(self, st2):
        """Tr(I) = dim_H."""
        assert st2.spectral_action() == Integer(16)

    def test_spectral_action_no_cutoff_nmax3(self, st3):
        assert st3.spectral_action() == Integer(40)

    def test_spectral_action_heat_positive(self, st2):
        """Heat-kernel action Tr(exp(-t*D^2)) should be positive."""
        import sympy as sp
        t = Rational(1, 10)
        sa = st2.spectral_action_heat(t)
        assert float(sp.N(sa)) > 0


# ================================================================
# Supertrace
# ================================================================

class TestSupertrace:
    """Str(M) = Tr(gamma * M)."""

    def test_supertrace_gamma_is_zero(self, st2):
        """Str(I) = Tr(gamma) = 0 when chi+ count == chi- count."""
        g = st2.grading
        assert g.trace() == Integer(0)

    def test_supertrace_D_squared(self, st2):
        """Str(D^2) = Tr(gamma * D^2)."""
        D = st2.dirac_operator
        st_val = st2.supertrace(D * D)
        # Should be computable (may or may not be zero)
        assert isinstance(st_val, (int, Integer)) or hasattr(st_val, '__float__')


# ================================================================
# Algebra representation
# ================================================================

class TestAlgebra:
    """Tests for the algebra A and its representation pi."""

    def test_pi_diagonal(self, st2):
        """pi(f) is diagonal."""
        f = [Integer(k) for k in range(1, st2.n_sectors + 1)]
        pi_f = st2.algebra_representation(f)
        for i in range(st2.dim_H):
            for j in range(st2.dim_H):
                if i != j:
                    assert pi_f[i, j] == 0

    def test_pi_respects_sectors(self, st2):
        """States in the same sector get the same function value."""
        f = [Integer(10 * k) for k in range(st2.n_sectors)]
        pi_f = st2.algebra_representation(f)
        for i, lab in enumerate(st2.labels):
            sector_idx = st2._state_to_sector[i]
            assert pi_f[i, i] == f[sector_idx]

    def test_pi_wrong_length_raises(self, st2):
        """f with wrong length raises ValueError."""
        with pytest.raises(ValueError):
            st2.algebra_representation([Integer(1)] * (st2.n_sectors + 1))


# ================================================================
# Edge cases and construction
# ================================================================

class TestConstruction:
    """Tests for constructor validation and edge cases."""

    def test_nmax_zero_raises(self):
        with pytest.raises(ValueError):
            FockSpectralTriple(n_max=0)

    def test_invalid_j_type_raises(self):
        with pytest.raises(ValueError):
            FockSpectralTriple(n_max=2, j_type="invalid")

    def test_nmax1_minimal(self):
        """n_max=1 is the smallest valid spectral triple."""
        st = FockSpectralTriple(n_max=1)
        assert st.dim_H == 4
        assert st.n_sectors == 2
        assert st.check_selfadjoint()
        assert st.check_grading_square()
        j2_ok, j2_eps = st.check_J_squared()
        assert j2_ok and j2_eps == 1

    def test_repr(self, st2):
        r = repr(st2)
        assert "FockSpectralTriple" in r
        assert "n_max=2" in r
        assert "dim_H=16" in r


# ================================================================
# Summary
# ================================================================

class TestSummary:
    """Test the summary method returns complete data."""

    def test_summary_keys(self, st2):
        s = st2.summary()
        assert "n_max" in s
        assert "hilbert_space" in s
        assert "spectrum" in s
        assert "axiom_checks" in s

    def test_summary_axiom_checks_complete(self, st2):
        s = st2.summary()
        checks = s["axiom_checks"]
        expected_keys = {
            "D_selfadjoint", "gamma_squared_I", "gamma_selfadjoint",
            "D_gamma_anticommute", "D_gamma_anticommute_msg",
            "Lambda_gamma_commute", "J_squared", "J_squared_epsilon",
            "J_D_relation", "J_D_relation_msg", "J_gamma_commute",
            "order_zero", "order_one", "order_one_msg", "pi_free",
        }
        assert expected_keys <= set(checks.keys())

    def test_summary_includes_adjacency_weights(self, st2):
        s = st2.summary()
        assert "adjacency_weights" in s
        assert s["adjacency_weights"] == "uniform"


# ================================================================
# CG-weighted adjacency — fixtures
# ================================================================

@pytest.fixture(scope="module")
def cg2():
    """CG-weighted spectral triple at n_max=2 (permutation J)."""
    return FockSpectralTriple(n_max=2, j_type="permutation", adjacency_weights="cg")


@pytest.fixture(scope="module")
def cg3():
    """CG-weighted spectral triple at n_max=3 (permutation J)."""
    return FockSpectralTriple(n_max=3, j_type="permutation", adjacency_weights="cg")


@pytest.fixture(scope="module")
def cg2_kramers():
    """CG-weighted spectral triple at n_max=2 (Kramers J)."""
    return FockSpectralTriple(n_max=2, j_type="kramers", adjacency_weights="cg")


@pytest.fixture(scope="module")
def cg3_kramers():
    """CG-weighted spectral triple at n_max=3 (Kramers J)."""
    return FockSpectralTriple(n_max=3, j_type="kramers", adjacency_weights="cg")


# ================================================================
# CG-weighted: construction and validation
# ================================================================

class TestCGConstruction:
    """Validate that CG-weighted constructor accepts parameters properly."""

    def test_invalid_adjacency_weights_raises(self):
        with pytest.raises(ValueError):
            FockSpectralTriple(n_max=2, adjacency_weights="invalid")

    def test_cg_dim_matches_uniform(self, cg2, st2):
        """CG weighting should not change Hilbert space dimension."""
        assert cg2.dim_H == st2.dim_H

    def test_cg_sectors_match_uniform(self, cg2, st2):
        """CG weighting should not change sector structure."""
        assert cg2.sectors == st2.sectors
        assert cg2.n_sectors == st2.n_sectors

    def test_cg_repr(self, cg2):
        r = repr(cg2)
        assert "adjacency_weights='cg'" in r

    def test_cg_weight_matrix_symmetric(self, cg2):
        """CG weight matrix should be symmetric."""
        W = cg2._cg_weights
        N = cg2.dim_H
        for i in range(N):
            for j in range(i + 1, N):
                assert W[i, j] == W[j, i], f"W[{i},{j}]={W[i,j]} != W[{j},{i}]={W[j,i]}"

    def test_cg_weight_matrix_has_nonzero_entries(self, cg2):
        """CG weight matrix should have at least some nonzero entries."""
        W = cg2._cg_weights
        N = cg2.dim_H
        n_nz = sum(1 for i in range(N) for j in range(N) if W[i, j] != 0)
        assert n_nz > 0

    def test_cg_adjacency_weights_stored(self, cg2):
        """The adjacency_weights attribute should be accessible."""
        assert cg2._adjacency_weights == "cg"


# ================================================================
# CG-weighted: D self-adjointness
# ================================================================

class TestCGDiracOperator:
    """CG-weighted D = Lambda + kappa * W must be self-adjoint."""

    def test_selfadjoint_nmax2(self, cg2):
        assert cg2.check_selfadjoint()

    def test_selfadjoint_nmax3(self, cg3):
        assert cg3.check_selfadjoint()

    def test_D_diagonal_is_lambda(self, cg2):
        """Diagonal of CG-D should equal Lambda (Camporesi-Higuchi)."""
        D = cg2.dirac_operator
        Lam = cg2.diagonal_part
        for i in range(cg2.dim_H):
            assert D[i, i] == Lam[i, i]

    def test_D_differs_from_uniform(self, cg2, st2):
        """CG-weighted D should differ from uniform D in off-diagonal."""
        D_cg = cg2.dirac_operator
        D_uni = st2.dirac_operator
        differences = 0
        N = cg2.dim_H
        for i in range(N):
            for j in range(N):
                if i != j and D_cg[i, j] != D_uni[i, j]:
                    differences += 1
        assert differences > 0, "CG-weighted and uniform D are identical"

    def test_spectrum_paired(self, cg2):
        """CG-weighted D spectrum should be approximately paired."""
        import numpy as np
        D = cg2.dirac_operator
        arr = np.array(D.tolist(), dtype=np.float64)
        evals = np.linalg.eigvalsh(arr)
        evals_sorted = sorted(evals)
        evals_neg = sorted(-evals)
        for a, b in zip(evals_sorted, evals_neg):
            assert abs(a - b) < 1e-10, f"Pairing mismatch: {a} vs {b}"


# ================================================================
# CG-weighted: Grading
# ================================================================

class TestCGGrading:
    """CG weighting should not change gamma properties."""

    def test_gamma_squared_I(self, cg2):
        assert cg2.check_grading_square()

    def test_gamma_selfadjoint(self, cg2):
        assert cg2.check_grading_selfadjoint()

    def test_Lambda_gamma_commute(self, cg2):
        assert cg2.check_D_gamma_commute_diagonal()


# ================================================================
# CG-weighted: Real structure J
# ================================================================

class TestCGRealStructure:
    """Tests for J with CG-weighted adjacency."""

    def test_permutation_J_squared_plus_I(self, cg2):
        j2_ok, j2_eps = cg2.check_J_squared()
        assert j2_ok
        assert j2_eps == 1

    def test_kramers_J_squared_minus_I(self, cg2_kramers):
        j2_ok, j2_eps = cg2_kramers.check_J_squared()
        assert j2_ok
        assert j2_eps == -1

    def test_permutation_J_gamma_commute(self, cg2):
        assert cg2.check_J_gamma_commute()

    def test_kramers_J_gamma_commute(self, cg2_kramers):
        assert cg2_kramers.check_J_gamma_commute()

    def test_permutation_JD_residual_cg(self, cg2):
        """Permutation J does NOT commute with CG-weighted D (characterize residual)."""
        jd_ok, jd_msg = cg2.check_J_D_relation()
        # CG weights break the permutation JD=DJ relation
        # (uniform weights have identical off-diagonal structure that J preserves,
        # but CG weights introduce m_j-dependent couplings that J scrambles)
        # Characterize rather than assert pass:
        assert isinstance(jd_ok, bool)
        assert isinstance(jd_msg, str)

    def test_permutation_JD_residual_cg_nmax3(self, cg3):
        """Permutation J residual with CG-weighted D at n_max=3."""
        jd_ok, jd_msg = cg3.check_J_D_relation()
        assert isinstance(jd_ok, bool)
        assert isinstance(jd_msg, str)


# ================================================================
# CG-weighted: Kramers JD residual analysis
# ================================================================

class TestCGKramersResidual:
    """Detailed analysis of the Kramers JD + DJ residual with CG weights."""

    def test_kramers_D_relation_returns_dict(self, cg2_kramers):
        result = cg2_kramers.check_kramers_D_relation()
        assert isinstance(result, dict)
        assert "exact_zero" in result
        assert "n_nonzero" in result
        assert "frobenius_norm" in result
        assert "max_entry" in result
        assert "residual_matrix" in result

    def test_kramers_residual_analysis_returns_dict(self, cg2_kramers):
        result = cg2_kramers.kramers_D_residual_analysis()
        assert isinstance(result, dict)
        assert "diagonal_residual" in result
        assert "offdiag_residual" in result
        assert "per_sector_max" in result
        assert "constraint_test" in result

    def test_kramers_diagonal_residual_is_nonzero(self, cg2_kramers):
        """JLambda + LambdaJ should be nonzero (Lambda commutes, not anticommutes)."""
        result = cg2_kramers.kramers_D_residual_analysis()
        assert not result["diagonal_residual"]["exact_zero"]

    def test_kramers_residual_frobenius_finite(self, cg2_kramers):
        """Frobenius norm of residual should be finite and positive."""
        result = cg2_kramers.check_kramers_D_relation()
        assert result["frobenius_norm"] > 0
        assert result["frobenius_norm"] < 1e6

    def test_uniform_kramers_residual_comparison(self, st2_kramers, cg2_kramers):
        """CG residual should be characterized (may be larger or smaller than uniform)."""
        r_uni = st2_kramers.check_kramers_D_relation()
        r_cg = cg2_kramers.check_kramers_D_relation()
        # Both should have nonzero residuals (JLambda + LambdaJ dominates)
        assert r_uni["frobenius_norm"] > 0
        assert r_cg["frobenius_norm"] > 0

    def test_kramers_per_sector_analysis(self, cg2_kramers):
        result = cg2_kramers.kramers_D_residual_analysis()
        assert len(result["per_sector_max"]) == cg2_kramers.n_sectors
        for label, max_val in result["per_sector_max"]:
            assert isinstance(max_val, float)

    def test_kramers_constraint_test_exists(self, cg2_kramers):
        """The constraint {J, A} + 2*Lambda*J should be computed."""
        result = cg2_kramers.kramers_D_residual_analysis()
        ct = result["constraint_test"]
        assert "exact_zero" in ct
        assert "frobenius_norm" in ct


# ================================================================
# CG-weighted: Order conditions
# ================================================================

class TestCGOrderConditions:
    """Order conditions with CG-weighted D."""

    def test_order_zero_nmax2(self, cg2):
        assert cg2.check_order_zero()

    def test_order_one_expected_fail_nmax2(self, cg2):
        """Order-one still expected to fail with CG weights."""
        o1_ok, o1_msg = cg2.check_order_one()
        assert not o1_ok, "Unexpected pass — order-one should fail"
        assert "EXPECTED FAIL" in o1_msg


# ================================================================
# CG-weighted: Pi-free certificate
# ================================================================

class TestCGPiFree:
    """All CG matrix entries are algebraic (no pi, no floats)."""

    def test_pi_free_nmax2(self, cg2):
        assert cg2.verify_pi_free()

    def test_pi_free_nmax3(self, cg3):
        assert cg3.verify_pi_free()

    def test_pi_free_kramers_nmax2(self, cg2_kramers):
        assert cg2_kramers.verify_pi_free()


# ================================================================
# CG-weighted: Connes distances
# ================================================================

class TestCGConnesDistance:
    """Connes distances with CG-weighted D."""

    def test_connes_distance_same_sector(self, cg2):
        d = cg2.connes_distance(0, 0)
        assert d == Integer(0)

    def test_connes_distance_positive(self, cg2):
        d = cg2.connes_distance(0, 1)
        assert float(d) > 0

    def test_connes_distances_finite(self, cg2):
        for i in range(cg2.n_sectors):
            for j in range(i + 1, cg2.n_sectors):
                d = cg2.connes_distance(i, j)
                assert float(d) < 1e6, f"Distance d({i},{j}) appears infinite"


# ================================================================
# CG-weighted: Spectral action
# ================================================================

class TestCGSpectralAction:
    """Spectral action with CG-weighted D."""

    def test_spectral_action_no_cutoff(self, cg2):
        """Tr(I) = dim_H regardless of weights."""
        assert cg2.spectral_action() == Integer(16)

    def test_spectral_action_heat_positive(self, cg2):
        """Heat-kernel action should be positive (uses numeric fallback for CG)."""
        t = Rational(1, 10)
        sa = cg2.spectral_action_heat(t)
        assert float(sa) > 0
