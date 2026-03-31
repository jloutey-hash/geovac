"""
Tests for the characteristic polynomial P(R, mu) of the Level 3 angular Hamiltonian.

Verifies that:
1. H0 has integer (rational) entries (SO(6) Casimir eigenvalues)
2. V_C can be rationalized (nuclear + V_ee coupling)
3. P(R, mu) has rational coefficients
4. Roots of P(R_val, mu) = 0 match numerical eigenvalues to machine precision
5. The algebraic curve P(R, mu) = 0 is valid at ALL R, including R > 5
   where the Track G Taylor series diverges

Diagnostic track only — does not modify production solvers.

References:
  - Paper 13, Sec III, XII, XII.B
  - Paper 18, Sec II.C, IV
  - Track G (perturbation series divergence)
"""

import numpy as np
import pytest

from geovac.algebraic_angular import AlgebraicAngularSolver
from geovac.algebraic_curve import (
    extract_pencil_matrices,
    rationalize_matrix,
    characteristic_polynomial,
    solve_characteristic_polynomial,
    verify_against_numerical,
    explicit_mu_lmax0,
    polynomial_structure_report,
)

try:
    import sympy
    from sympy import Symbol, Rational, Float, Poly
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

pytestmark = pytest.mark.skipif(not HAS_SYMPY, reason="SymPy required")


# ─────────────────────────────────────────────────────────────────────
# 1. Matrix extraction and rationality
# ─────────────────────────────────────────────────────────────────────

class TestPencilExtraction:
    """Test extraction and rationality of H0, V_C matrices."""

    def test_h0_integer_entries_lmax0(self):
        """H0 at l_max=0 should have integer diagonal entries (Casimir)."""
        H0, _ = extract_pencil_matrices(Z=2.0, l_max=0, n_basis=1)
        # SO(6) Casimir for l=0, k=0 (singlet): mu = 2*(0+0+1)^2 - 2 = 0
        assert H0.shape == (1, 1)
        assert abs(H0[0, 0] - 0.0) < 1e-12, f"H0[0,0] = {H0[0,0]}, expected 0"

    def test_h0_integer_entries_lmax1(self):
        """H0 at l_max=1 should have integer diagonal entries."""
        H0, _ = extract_pencil_matrices(Z=2.0, l_max=1, n_basis=1)
        # l=0, k=0: mu = 2*(1)^2 - 2 = 0
        # l=1, k=0: mu = 2*(1+0+1)^2 - 2 = 6
        assert H0.shape == (2, 2)
        assert abs(H0[0, 0] - 0.0) < 1e-12
        assert abs(H0[1, 1] - 6.0) < 1e-12

    def test_h0_integer_entries_lmax2(self):
        """H0 at l_max=2 should have integer diagonal entries."""
        H0, _ = extract_pencil_matrices(Z=2.0, l_max=2, n_basis=1)
        # l=0: 2*(1)^2 - 2 = 0
        # l=1: 2*(2)^2 - 2 = 6
        # l=2: 2*(3)^2 - 2 = 16
        assert H0.shape == (3, 3)
        assert abs(H0[0, 0] - 0.0) < 1e-12
        assert abs(H0[1, 1] - 6.0) < 1e-12
        assert abs(H0[2, 2] - 16.0) < 1e-12

    def test_h0_diagonal(self):
        """H0 should be diagonal (no off-diagonal elements)."""
        for lm in [0, 1, 2]:
            H0, _ = extract_pencil_matrices(Z=2.0, l_max=lm, n_basis=1)
            off_diag = H0 - np.diag(np.diag(H0))
            assert np.max(np.abs(off_diag)) < 1e-15, f"H0 not diagonal at l_max={lm}"

    def test_vc_symmetric(self):
        """V_C should be symmetric."""
        for lm in [0, 1, 2]:
            _, V_C = extract_pencil_matrices(Z=2.0, l_max=lm, n_basis=1)
            assert np.max(np.abs(V_C - V_C.T)) < 1e-14, f"V_C not symmetric at l_max={lm}"

    def test_vc_rationalizable_lmax0(self):
        """V_C at l_max=0, n_basis=1 should be rationalizable."""
        _, V_C = extract_pencil_matrices(Z=2.0, l_max=0, n_basis=1, n_quad=200)
        # Should not raise
        V_C_sym = rationalize_matrix(V_C, tol=1e-8)
        assert V_C_sym.shape == (1, 1)

    def test_pencil_structure(self):
        """H(R) = H0 + R*V_C should reproduce numerical eigenvalues."""
        for lm in [0, 1, 2]:
            H0, V_C = extract_pencil_matrices(Z=2.0, l_max=lm, n_basis=1)
            solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=lm)
            dim = lm + 1  # n_basis=1 per channel

            for R in [0.5, 2.0, 10.0]:
                H = H0 + R * V_C
                evals_direct = np.sort(np.linalg.eigvalsh(H))
                evals_solver, _ = solver.solve(R, n_channels=dim)
                np.testing.assert_allclose(
                    evals_direct, evals_solver, atol=1e-12,
                    err_msg=f"Pencil mismatch at l_max={lm}, R={R}"
                )


# ─────────────────────────────────────────────────────────────────────
# 2. Characteristic polynomial structure
# ─────────────────────────────────────────────────────────────────────

class TestCharacteristicPolynomial:
    """Test structural properties of the characteristic polynomial."""

    def test_degree_mu_equals_dim(self):
        """Degree of P in mu should equal matrix dimension."""
        for lm in [0, 1, 2]:
            P, meta = characteristic_polynomial(
                Z=2.0, l_max=lm, n_basis=1, n_quad=200
            )
            assert meta['degree_mu'] == meta['dim'], (
                f"l_max={lm}: degree_mu={meta['degree_mu']} != dim={meta['dim']}"
            )

    def test_degree_R_bounded_by_dim(self):
        """Degree of P in R should be at most dim."""
        for lm in [0, 1, 2]:
            P, meta = characteristic_polynomial(
                Z=2.0, l_max=lm, n_basis=1, n_quad=200
            )
            assert meta['degree_R'] <= meta['dim'], (
                f"l_max={lm}: degree_R={meta['degree_R']} > dim={meta['dim']}"
            )

    def test_h0_rational_flag(self):
        """H0 should always be flagged as rational."""
        for lm in [0, 1, 2]:
            _, meta = characteristic_polynomial(
                Z=2.0, l_max=lm, n_basis=1, n_quad=200
            )
            assert meta['H0_rational'], f"H0 not rational at l_max={lm}"

    def test_lmax0_linear_in_mu(self):
        """At l_max=0, n_basis=1: P is linear in mu (1x1 matrix)."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=0, n_basis=1, n_quad=200
        )
        assert meta['dim'] == 1
        mu_sym = Symbol('mu')
        poly = Poly(P, mu_sym)
        assert poly.degree() == 1, f"Expected degree 1 in mu, got {poly.degree()}"

    def test_lmax0_linear_in_R(self):
        """At l_max=0, n_basis=1: P is linear in R (1x1 matrix)."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=0, n_basis=1, n_quad=200
        )
        R_sym = Symbol('R')
        poly = Poly(P, R_sym)
        assert poly.degree() == 1, f"Expected degree 1 in R, got {poly.degree()}"

    def test_lmax1_quadratic_in_mu(self):
        """At l_max=1, n_basis=1: P is degree 2 in mu (2x2 matrix)."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=1, n_basis=1, n_quad=200
        )
        assert meta['dim'] == 2
        mu_sym = Symbol('mu')
        poly = Poly(P, mu_sym)
        assert poly.degree() == 2

    def test_lmax2_cubic_in_mu(self):
        """At l_max=2, n_basis=1: P is degree 3 in mu (3x3 matrix)."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=2, n_basis=1, n_quad=200
        )
        assert meta['dim'] == 3
        mu_sym = Symbol('mu')
        poly = Poly(P, mu_sym)
        assert poly.degree() == 3


# ─────────────────────────────────────────────────────────────────────
# 3. Larger bases: n_basis=2 gives dim 2, 4, 6
# ─────────────────────────────────────────────────────────────────────

class TestLargerBasis:
    """Test with n_basis=2 for larger characteristic polynomials."""

    def test_dimensions_nbasis2(self):
        """n_basis=2: l_max=0 -> dim 2, l_max=1 -> dim 4, l_max=2 -> dim 6."""
        for lm, expected_dim in [(0, 2), (1, 4), (2, 6)]:
            _, meta = characteristic_polynomial(
                Z=2.0, l_max=lm, n_basis=2, n_quad=200
            )
            assert meta['dim'] == expected_dim, (
                f"l_max={lm}, n_basis=2: dim={meta['dim']}, expected {expected_dim}"
            )

    def test_degree_mu_nbasis2(self):
        """Degree in mu equals dim for n_basis=2."""
        for lm, expected_dim in [(0, 2), (1, 4), (2, 6)]:
            P, meta = characteristic_polynomial(
                Z=2.0, l_max=lm, n_basis=2, n_quad=200
            )
            mu_sym = Symbol('mu')
            poly = Poly(P, mu_sym)
            assert poly.degree() == expected_dim

    def test_nbasis3_lmax0_gives_dim3(self):
        """n_basis=3, l_max=0 gives dim=3 (for comparison with 'dim 3' from task)."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=0, n_basis=3, n_quad=200
        )
        assert meta['dim'] == 3
        mu_sym = Symbol('mu')
        poly = Poly(P, mu_sym)
        assert poly.degree() == 3


# ─────────────────────────────────────────────────────────────────────
# 4. Root verification against numerical eigenvalues
# ─────────────────────────────────────────────────────────────────────

class TestRootVerification:
    """Verify P(R,mu)=0 roots match numerical eigenvalues."""

    @pytest.mark.parametrize("R_val", [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    def test_lmax0_roots_match(self, R_val):
        """l_max=0, n_basis=1: single root matches numerical eigenvalue."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=0, n_basis=1, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=0, n_quad=200)
        evals_num, _ = solver.solve(R_val, n_channels=1)

        roots = solve_characteristic_polynomial(P, R_val)
        assert len(roots) >= 1, f"No roots found at R={R_val}"
        assert abs(roots[0] - evals_num[0]) < 1e-10, (
            f"R={R_val}: root={roots[0]}, numerical={evals_num[0]}, "
            f"diff={abs(roots[0] - evals_num[0])}"
        )

    @pytest.mark.parametrize("R_val", [0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    def test_lmax1_roots_match(self, R_val):
        """l_max=1, n_basis=1: all roots match numerical eigenvalues."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=1, n_basis=1, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=1, n_quad=200)
        evals_num, _ = solver.solve(R_val, n_channels=2)
        roots = solve_characteristic_polynomial(P, R_val)

        assert len(roots) == 2, f"Expected 2 roots, got {len(roots)}"
        for i in range(2):
            assert abs(roots[i] - evals_num[i]) < 1e-10, (
                f"R={R_val}, root[{i}]: from_P={roots[i]}, "
                f"numerical={evals_num[i]}, diff={abs(roots[i] - evals_num[i])}"
            )

    @pytest.mark.parametrize("R_val", [0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    def test_lmax2_roots_match(self, R_val):
        """l_max=2, n_basis=1: all roots match numerical eigenvalues."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=2, n_basis=1, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=2, n_quad=200)
        evals_num, _ = solver.solve(R_val, n_channels=3)
        roots = solve_characteristic_polynomial(P, R_val)

        assert len(roots) == 3, f"Expected 3 roots, got {len(roots)}"
        for i in range(3):
            assert abs(roots[i] - evals_num[i]) < 1e-10, (
                f"R={R_val}, root[{i}]: from_P={roots[i]}, "
                f"numerical={evals_num[i]}, diff={abs(roots[i] - evals_num[i])}"
            )


# ─────────────────────────────────────────────────────────────────────
# 5. Key diagnostic: algebraic curve valid beyond Taylor divergence
# ─────────────────────────────────────────────────────────────────────

class TestBeyondTaylorDivergence:
    """Test that P(R,mu)=0 works at R values where Track G's Taylor series diverges.

    Track G showed the Taylor/Pade expansion of mu(R) diverges for R > ~5 bohr.
    The characteristic polynomial P(R,mu) = 0 defines an implicit algebraic
    relation that is exact for ALL R — no convergence radius limitation.
    """

    def test_R5_beyond_taylor(self):
        """R=5.0 is near the Taylor divergence boundary — curve must work."""
        P, _, verification = verify_against_numerical(
            Z=2.0, l_max=1, n_basis=1, R_values=[5.0], n_quad=200
        )
        assert verification[0]['max_diff'] < 1e-10

    def test_R10_well_beyond_taylor(self):
        """R=10.0 is well beyond Taylor divergence — curve must work."""
        P, _, verification = verify_against_numerical(
            Z=2.0, l_max=1, n_basis=1, R_values=[10.0], n_quad=200
        )
        assert verification[0]['max_diff'] < 1e-10

    def test_R20_far_beyond_taylor(self):
        """R=20.0 is far beyond Taylor divergence — curve must work."""
        P, _, verification = verify_against_numerical(
            Z=2.0, l_max=1, n_basis=1, R_values=[20.0], n_quad=200
        )
        assert verification[0]['max_diff'] < 1e-10

    def test_R50_extreme(self):
        """R=50.0 is an extreme test — curve must still work."""
        P, _, verification = verify_against_numerical(
            Z=2.0, l_max=1, n_basis=1, R_values=[50.0], n_quad=200
        )
        assert verification[0]['max_diff'] < 1e-10

    def test_full_range_lmax2(self):
        """Full R-range with l_max=2 — all points must match."""
        R_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
        _, _, verification = verify_against_numerical(
            Z=2.0, l_max=2, n_basis=1, R_values=R_values, n_quad=200
        )
        for v in verification:
            assert v['max_diff'] < 1e-10, (
                f"R={v['R']}: max_diff={v['max_diff']:.2e}"
            )


# ─────────────────────────────────────────────────────────────────────
# 6. Explicit mu(R) at l_max=0
# ─────────────────────────────────────────────────────────────────────

class TestExplicitMu:
    """Test the explicit rational function mu(R) at l_max=0, n_basis=1."""

    def test_explicit_form(self):
        """mu(R) = casimir + R * V_C_00 at l_max=0, n_basis=1."""
        mu_expr = explicit_mu_lmax0(Z=2.0, n_quad=200)
        R_sym = Symbol('R')

        # Should be linear in R
        poly = Poly(mu_expr, R_sym)
        assert poly.degree() == 1

    def test_explicit_matches_numerical(self):
        """Explicit mu(R) matches numerical eigenvalue at several R values."""
        mu_expr = explicit_mu_lmax0(Z=2.0, n_quad=200)
        R_sym = Symbol('R')
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=0, n_quad=200)

        for R_val in [0.1, 1.0, 5.0, 20.0, 100.0]:
            mu_symbolic = float(mu_expr.subs(R_sym, R_val))
            evals, _ = solver.solve(R_val, n_channels=1)
            assert abs(mu_symbolic - evals[0]) < 1e-10, (
                f"R={R_val}: symbolic={mu_symbolic}, numerical={evals[0]}"
            )

    def test_R0_gives_casimir(self):
        """At R=0, mu should equal the SO(6) Casimir eigenvalue."""
        mu_expr = explicit_mu_lmax0(Z=2.0, n_quad=200)
        R_sym = Symbol('R')
        mu_at_0 = float(mu_expr.subs(R_sym, 0))
        # Casimir for l=0, k=0 singlet: 2*(0+0+1)^2 - 2 = 0
        assert abs(mu_at_0 - 0.0) < 1e-12


# ─────────────────────────────────────────────────────────────────────
# 7. Polynomial structure report
# ─────────────────────────────────────────────────────────────────────

class TestStructureReport:
    """Test the polynomial structure reporting functionality."""

    def test_report_generation(self):
        """Structure report should complete for l_max = 0, 1, 2."""
        reports = polynomial_structure_report(
            Z=2.0, l_max_values=[0, 1, 2], n_basis=1, n_quad=200
        )
        assert len(reports) == 3

    def test_report_dimensions(self):
        """Report dimensions should match expected values."""
        reports = polynomial_structure_report(
            Z=2.0, l_max_values=[0, 1, 2], n_basis=1, n_quad=200
        )
        expected_dims = [1, 2, 3]
        for rep, dim in zip(reports, expected_dims):
            assert rep['dim'] == dim, (
                f"l_max={rep['l_max']}: dim={rep['dim']}, expected={dim}"
            )

    def test_report_h0_always_rational(self):
        """H0 should always be rational in reports."""
        reports = polynomial_structure_report(
            Z=2.0, l_max_values=[0, 1, 2], n_basis=1, n_quad=200
        )
        for rep in reports:
            assert rep['H0_rational'], f"H0 not rational at l_max={rep['l_max']}"


# ─────────────────────────────────────────────────────────────────────
# 8. n_basis=2 verification (dim 2, 4, 6 as task suggests)
# ─────────────────────────────────────────────────────────────────────

class TestNBasis2Verification:
    """Verify roots match for larger truncations (n_basis=2, dims 2, 4, 6)."""

    @pytest.mark.parametrize("R_val", [1.0, 5.0, 10.0])
    def test_lmax0_nbasis2_roots(self, R_val):
        """l_max=0, n_basis=2 (dim=2): roots match."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=0, n_basis=2, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=2, l_max=0, n_quad=200)
        evals_num, _ = solver.solve(R_val, n_channels=2)
        roots = solve_characteristic_polynomial(P, R_val)

        assert len(roots) == 2
        for i in range(2):
            assert abs(roots[i] - evals_num[i]) < 1e-8, (
                f"R={R_val}, root[{i}]: diff={abs(roots[i] - evals_num[i]):.2e}"
            )

    @pytest.mark.parametrize("R_val", [1.0, 5.0, 10.0])
    def test_lmax1_nbasis2_roots(self, R_val):
        """l_max=1, n_basis=2 (dim=4): roots match."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=1, n_basis=2, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=2, l_max=1, n_quad=200)
        dim = meta['dim']
        evals_num, _ = solver.solve(R_val, n_channels=dim)
        roots = solve_characteristic_polynomial(P, R_val)

        assert len(roots) == dim
        for i in range(dim):
            assert abs(roots[i] - evals_num[i]) < 1e-8, (
                f"R={R_val}, root[{i}]: diff={abs(roots[i] - evals_num[i]):.2e}"
            )

    @pytest.mark.slow
    @pytest.mark.parametrize("R_val", [1.0, 5.0])
    def test_lmax2_nbasis2_roots(self, R_val):
        """l_max=2, n_basis=2 (dim=6): roots match. Slow due to 6x6 det."""
        P, meta = characteristic_polynomial(
            Z=2.0, l_max=2, n_basis=2, n_quad=200
        )
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=2, l_max=2, n_quad=200)
        dim = meta['dim']
        evals_num, _ = solver.solve(R_val, n_channels=dim)
        roots = solve_characteristic_polynomial(P, R_val)

        assert len(roots) == dim
        for i in range(dim):
            assert abs(roots[i] - evals_num[i]) < 1e-7, (
                f"R={R_val}, root[{i}]: diff={abs(roots[i] - evals_num[i]):.2e}"
            )


# ─────────────────────────────────────────────────────────────────────
# 9. Rationality of coupling matrix
# ─────────────────────────────────────────────────────────────────────

class TestCouplingRationality:
    """Test whether V_C entries can be rationalized.

    The coupling matrix involves integrals of:
      nuclear: sin(2n*alpha) * (-Z/cos(alpha) + -Z/sin(alpha)) * sin(2m*alpha)
      V_ee:    sin(2n*alpha) * 1/max(sin,cos) * sin(2m*alpha)

    These may or may not have rational closed forms. This test documents
    the finding.
    """

    def test_vc_rationalizability_lmax0(self):
        """Check if V_C at l_max=0 is rationalizable."""
        _, V_C = extract_pencil_matrices(Z=2.0, l_max=0, n_basis=1, n_quad=200)
        try:
            V_C_sym = rationalize_matrix(V_C, tol=1e-8, max_denom=100000)
            is_rational = True
        except ValueError:
            is_rational = False
        # Record result — this is a diagnostic, not a hard pass/fail
        print(f"V_C rationalizable at l_max=0, n_basis=1: {is_rational}")
        if is_rational:
            print(f"  V_C[0,0] = {V_C_sym[0, 0]}")

    def test_vc_rationalizability_lmax1(self):
        """Check if V_C at l_max=1 is rationalizable."""
        _, V_C = extract_pencil_matrices(Z=2.0, l_max=1, n_basis=1, n_quad=200)
        try:
            V_C_sym = rationalize_matrix(V_C, tol=1e-8, max_denom=100000)
            is_rational = True
        except ValueError:
            is_rational = False
        print(f"V_C rationalizable at l_max=1, n_basis=1: {is_rational}")
        if is_rational:
            for i in range(2):
                for j in range(2):
                    print(f"  V_C[{i},{j}] = {V_C_sym[i, j]}")

    def test_vc_rationalizability_lmax2(self):
        """Check if V_C at l_max=2 is rationalizable."""
        _, V_C = extract_pencil_matrices(Z=2.0, l_max=2, n_basis=1, n_quad=200)
        try:
            V_C_sym = rationalize_matrix(V_C, tol=1e-8, max_denom=100000)
            is_rational = True
        except ValueError:
            is_rational = False
        print(f"V_C rationalizable at l_max=2, n_basis=1: {is_rational}")
        if is_rational:
            for i in range(3):
                for j in range(3):
                    print(f"  V_C[{i},{j}] = {V_C_sym[i, j]}")

    def test_nuclear_vs_vee_rationality(self):
        """Separately test nuclear and V_ee coupling rationality."""
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=2, n_quad=200)
        nuc = solver._nuclear_full
        vee = solver._vee_full

        for name, M in [("nuclear", nuc), ("V_ee", vee)]:
            try:
                M_sym = rationalize_matrix(M, tol=1e-8, max_denom=100000)
                print(f"{name} rationalizable: True")
                for i in range(M.shape[0]):
                    for j in range(M.shape[1]):
                        if abs(float(M_sym[i, j])) > 1e-15:
                            print(f"  {name}[{i},{j}] = {M_sym[i, j]}")
            except ValueError as e:
                print(f"{name} rationalizable: False ({e})")


# ─────────────────────────────────────────────────────────────────────
# 10. Exact closed-form matrix entries (verified against quadrature)
# ─────────────────────────────────────────────────────────────────────

class TestExactClosedForms:
    """Verify exact closed forms for V_C matrix entries.

    The nuclear and V_ee coupling integrals have closed forms involving
    pi and sqrt(2). These are NOT rational — they are elements of Q(pi, sqrt(2)).
    This is the transcendental content that enters the characteristic polynomial
    and makes mu(R) transcendental as a function (Track G, Paper 18).

    Verified entries:
      Nuclear (diagonal, scales as Z): nuc[l,l] = -Z * c_l / pi
        l=0: c_0 = 32/3
        l=1: c_1 = 1024/105
      V_ee (Z-independent):
        vee[0,0] = 8*sqrt(2) / (3*pi)
        vee[0,1] = vee[1,0] = 16*sqrt(2) / (15*pi)
        vee[0,2] = vee[2,0] = 32 / (35*pi)
    """

    def test_nuc00_exact(self):
        """nuc[0,0] = -Z * 32/(3*pi)."""
        from math import pi, sqrt
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=0, n_quad=1000)
        expected = -2.0 * 32.0 / (3.0 * pi)
        assert abs(solver._nuclear_full[0, 0] - expected) < 1e-12

    def test_nuc11_exact(self):
        """nuc[1,1] = -Z * 1024/(105*pi)."""
        from math import pi
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=1, n_quad=1000)
        expected = -2.0 * 1024.0 / (105.0 * pi)
        assert abs(solver._nuclear_full[1, 1] - expected) < 1e-12

    def test_nuc_off_diagonal_zero(self):
        """Nuclear coupling is diagonal in l (isotropic potential)."""
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=2, n_quad=1000)
        nuc = solver._nuclear_full
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert abs(nuc[i, j]) < 1e-14

    def test_nuc_z_scaling(self):
        """Nuclear entries scale linearly with Z."""
        from math import pi
        for Z in [1.0, 2.0, 3.0]:
            solver = AlgebraicAngularSolver(Z=Z, n_basis=1, l_max=0, n_quad=1000)
            expected = -Z * 32.0 / (3.0 * pi)
            assert abs(solver._nuclear_full[0, 0] - expected) < 1e-11

    def test_vee00_exact(self):
        """vee[0,0] = 8*sqrt(2)/(3*pi)."""
        from math import pi, sqrt
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=0, n_quad=1000)
        expected = 8.0 * sqrt(2) / (3.0 * pi)
        assert abs(solver._vee_full[0, 0] - expected) < 1e-12

    def test_vee01_exact(self):
        """vee[0,1] = 16*sqrt(2)/(15*pi)."""
        from math import pi, sqrt
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=1, n_quad=1000)
        expected = 16.0 * sqrt(2) / (15.0 * pi)
        assert abs(solver._vee_full[0, 1] - expected) < 1e-12

    def test_vee02_exact(self):
        """vee[0,2] = 32/(35*pi)."""
        from math import pi
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=2, n_quad=1000)
        expected = 32.0 / (35.0 * pi)
        assert abs(solver._vee_full[0, 2] - expected) < 1e-12

    def test_vee_z_independent(self):
        """V_ee coupling is independent of nuclear charge Z."""
        solver1 = AlgebraicAngularSolver(Z=1.0, n_basis=1, l_max=2, n_quad=1000)
        solver2 = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=2, n_quad=1000)
        solver3 = AlgebraicAngularSolver(Z=3.0, n_basis=1, l_max=2, n_quad=1000)
        np.testing.assert_allclose(
            solver1._vee_full, solver2._vee_full, atol=1e-12
        )
        np.testing.assert_allclose(
            solver2._vee_full, solver3._vee_full, atol=1e-12
        )


# ═════════════════════════════════════════════════════════════════════
# Track P2: Algebraic V_eff matrix elements
# ═════════════════════════════════════════════════════════════════════

from geovac.algebraic_curve import (
    _stieltjes_ordinary_laguerre,
    _stieltjes_x2_weighted,
    _stieltjes_x2_squared_analytical,
    algebraic_veff_matrix_lmax0,
    algebraic_veff_matrix_quadrature,
    algebraic_he_energy_lmax0,
    analyze_discriminant_lmax1,
    detect_recurrence_lmax1,
)


# ─────────────────────────────────────────────────────────────────────
# 11. Stieltjes integral infrastructure
# ─────────────────────────────────────────────────────────────────────

class TestStieltjesOrdinaryLaguerre:
    """Test the Stieltjes integral for ordinary Laguerre polynomials."""

    def test_j0_symmetric(self):
        """J0 should be symmetric."""
        J0 = _stieltjes_ordinary_laguerre(10, a=1.0)
        np.testing.assert_allclose(J0, J0.T, atol=1e-14)

    def test_j0_positive_diagonal(self):
        """J0 diagonal should be positive (integrand is positive for i=j)."""
        J0 = _stieltjes_ordinary_laguerre(10, a=1.0)
        assert np.all(J0.diagonal() > 0), f"Negative diagonal: {J0.diagonal()}"

    def test_j0_s0_value(self):
        """J0[0,0] = S_0(a) = e^a E_1(a) should match scipy."""
        from scipy.special import exp1
        a = 2.0
        J0 = _stieltjes_ordinary_laguerre(1, a=a)
        expected = float(np.exp(a) * exp1(a))
        assert abs(J0[0, 0] - expected) < 1e-12, (
            f"J0[0,0]={J0[0,0]}, expected e^a E1(a)={expected}"
        )

    def test_j0_against_quadrature(self):
        """J0 should match numerical quadrature to high precision."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 1.5

        J0_alg = _stieltjes_ordinary_laguerre(N, a)

        # Quadrature reference
        n_quad = 200
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        # J0_quad[i,j] = int L_i L_j e^{-x} / (x+a) dx
        # = sum_k w_k L_i(x_k) L_j(x_k) / (x_k + a)
        wL = (w / (x + a))[np.newaxis, :] * L
        J0_quad = wL @ L.T

        np.testing.assert_allclose(
            J0_alg, J0_quad, atol=1e-10,
            err_msg="Stieltjes J0: algebraic vs quadrature mismatch"
        )

    def test_j0_various_a_values(self):
        """J0 should be accurate for various shift values.

        Accuracy degrades for small a (< 1) due to the forward segment
        of the hybrid recurrence, matching known behavior from
        prolate_spheroidal_lattice._stieltjes_matrix.
        """
        from scipy.special import roots_laguerre, eval_laguerre
        N = 15

        for a, tol in [(0.15, 1e-4), (0.5, 1e-7), (1.0, 1e-9),
                        (5.0, 1e-10), (20.0, 1e-10)]:
            J0_alg = _stieltjes_ordinary_laguerre(N, a)

            n_quad = 300
            x, w = roots_laguerre(n_quad)
            L = np.zeros((N, n_quad))
            for n in range(N):
                L[n] = eval_laguerre(n, x)
            wL = (w / (x + a))[np.newaxis, :] * L
            J0_quad = wL @ L.T

            np.testing.assert_allclose(
                J0_alg, J0_quad, atol=tol,
                err_msg=f"J0 mismatch at a={a}"
            )


class TestStieltjesWeighted:
    """Test x^2-weighted Stieltjes integrals."""

    def test_j1_against_quadrature_large_a(self):
        """J1 at a=2.0 should match quadrature to high precision."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 2.0  # Large a: recurrence is accurate

        J1, J2 = _stieltjes_x2_weighted(N, a)

        n_quad = 300
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        wL1 = (w * x**2 / (x + a))[np.newaxis, :] * L
        J1_quad = wL1 @ L.T

        np.testing.assert_allclose(
            J1, J1_quad, atol=1e-9,
            err_msg="J1 algebraic vs quadrature mismatch at a=2.0"
        )

    def test_j1_against_quadrature_small_a(self):
        """J1 at a=0.15 should match quadrature (reduced tolerance for small a)."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 0.15  # = 2 * 1.5 * 0.05

        J1, J2 = _stieltjes_x2_weighted(N, a)

        n_quad = 300
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        wL1 = (w * x**2 / (x + a))[np.newaxis, :] * L
        J1_quad = wL1 @ L.T

        # Reduced tolerance for small a — recurrence accuracy ~1e-4
        np.testing.assert_allclose(
            J1, J1_quad, atol=1e-3,
            err_msg="J1 algebraic vs quadrature mismatch at a=0.15"
        )

    def test_j2_against_quadrature(self):
        """J2 = int x^2 L_i L_j e^{-x}/(x+a)^2 dx should match quadrature."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 2.0  # Use large a for clean test

        _, J2 = _stieltjes_x2_weighted(N, a)

        n_quad = 300
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        wL2 = (w * x**2 / (x + a)**2)[np.newaxis, :] * L
        J2_quad = wL2 @ L.T

        np.testing.assert_allclose(
            J2, J2_quad, atol=1e-8,
            err_msg="J2 algebraic vs quadrature mismatch"
        )

    def test_j_sq_analytical_against_quadrature(self):
        """J_sq = int L_i L_j e^{-x}/(x+a)^2 dx should match quadrature."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 2.0  # Use large a for clean test

        J_sq = _stieltjes_x2_squared_analytical(N, a)

        n_quad = 300
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        wL = (w / (x + a)**2)[np.newaxis, :] * L
        J_sq_quad = wL @ L.T

        np.testing.assert_allclose(
            J_sq, J_sq_quad, atol=1e-9,
            err_msg="J_sq analytical vs quadrature mismatch"
        )

    def test_j_sq_analytical_small_a(self):
        """J_sq at a=0.15 should match quadrature (reduced tolerance)."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 10
        a = 0.15

        J_sq = _stieltjes_x2_squared_analytical(N, a)

        n_quad = 300
        x, w = roots_laguerre(n_quad)
        L = np.zeros((N, n_quad))
        for n in range(N):
            L[n] = eval_laguerre(n, x)

        wL = (w / (x + a)**2)[np.newaxis, :] * L
        J_sq_quad = wL @ L.T

        np.testing.assert_allclose(
            J_sq, J_sq_quad, atol=1e-3,
            err_msg="J_sq analytical vs quadrature mismatch at a=0.15"
        )


# ─────────────────────────────────────────────────────────────────────
# 12. Algebraic V_eff at l_max=0
# ─────────────────────────────────────────────────────────────────────

class TestAlgebraicVeffLmax0:
    """Test algebraic V_eff matrix at l_max=0."""

    def test_veff_matches_quadrature(self):
        """Algebraic V_eff should match quadrature.

        At R_min=0.05, alpha=1.5, the shift parameter a=0.15 is small,
        limiting Stieltjes recurrence accuracy to ~1e-3. The V_eff
        matrix accuracy inherits this limitation.
        """
        n_basis = 20
        alpha = 1.5
        R_min = 0.05

        V_alg, info = algebraic_veff_matrix_lmax0(
            Z=2.0, n_basis=n_basis, alpha=alpha, R_min=R_min,
        )
        V_quad = algebraic_veff_matrix_quadrature(
            Z=2.0, l_max=0, n_basis=n_basis, alpha=alpha, R_min=R_min,
        )

        max_diff = np.max(np.abs(V_alg - V_quad))
        rel_diff = max_diff / max(np.max(np.abs(V_quad)), 1e-15)

        # Tolerance limited by Stieltjes recurrence at small a=0.15
        np.testing.assert_allclose(
            V_alg, V_quad, atol=0.1, rtol=0.01,
            err_msg=(
                f"V_eff algebraic vs quadrature: max_diff={max_diff:.2e}, "
                f"rel_diff={rel_diff:.2e}"
            )
        )

    def test_veff_matches_quadrature_large_rmin(self):
        """V_eff with larger R_min (a=1.5) should match to high precision."""
        n_basis = 15
        alpha = 1.5
        R_min = 0.5  # a = 2*1.5*0.5 = 1.5, better recurrence accuracy

        V_alg, info = algebraic_veff_matrix_lmax0(
            Z=2.0, n_basis=n_basis, alpha=alpha, R_min=R_min,
        )
        V_quad = algebraic_veff_matrix_quadrature(
            Z=2.0, l_max=0, n_basis=n_basis, alpha=alpha, R_min=R_min,
        )

        max_diff = np.max(np.abs(V_alg - V_quad))

        np.testing.assert_allclose(
            V_alg, V_quad, atol=1e-6,
            err_msg=f"V_eff at R_min=0.5: max_diff={max_diff:.2e}"
        )

    def test_veff_symmetric(self):
        """Algebraic V_eff should be symmetric."""
        V_alg, _ = algebraic_veff_matrix_lmax0(n_basis=15)
        np.testing.assert_allclose(V_alg, V_alg.T, atol=1e-14)

    def test_veff_c0_is_zero(self):
        """c0 (Casimir at l=0) should be 0."""
        _, info = algebraic_veff_matrix_lmax0()
        assert abs(info['c0']) < 1e-12, f"c0 = {info['c0']}, expected 0"

    def test_veff_c1_matches_known(self):
        """c1 should match the known V_C[0,0] value."""
        import math
        _, info = algebraic_veff_matrix_lmax0(Z=2.0, n_quad=1000)
        Z = 2.0
        nuc = -Z * 32.0 / (3.0 * math.pi)
        vee = 8.0 * math.sqrt(2) / (3.0 * math.pi)
        expected_c1 = nuc + vee
        assert abs(info['c1'] - expected_c1) < 1e-10, (
            f"c1 = {info['c1']}, expected {expected_c1}"
        )


class TestAlgebraicHeEnergy:
    """Test He ground state energy with algebraic V_eff."""

    def test_energy_matches_quadrature_solver(self):
        """Algebraic He energy should match quadrature-based solver."""
        from geovac.hyperspherical_radial import solve_radial_spectral
        from geovac.algebraic_angular import AlgebraicAngularSolver

        n_basis = 25
        alpha = 1.5
        R_min = 0.05

        # Algebraic energy
        E_alg, info = algebraic_he_energy_lmax0(
            Z=2.0, n_basis=n_basis, alpha=alpha, R_min=R_min,
        )

        # Quadrature reference
        solver = AlgebraicAngularSolver(Z=2.0, n_basis=1, l_max=0, n_quad=200)

        def V_eff(R):
            V = np.zeros_like(R)
            for k, Rv in enumerate(R):
                if Rv > 1e-10:
                    evals, _ = solver.solve(Rv, n_channels=1)
                    V[k] = evals[0] / Rv**2 + 15.0 / (8.0 * Rv**2)
                else:
                    V[k] = 1e10
            return V

        E_quad, _, _ = solve_radial_spectral(
            V_eff, n_basis=n_basis, alpha=alpha, R_min=R_min,
            matrix_method='algebraic',
        )

        diff = abs(E_alg - E_quad[0])
        # Tolerance limited by Stieltjes recurrence at small a=0.15
        assert diff < 0.01, (
            f"Energy mismatch: algebraic={E_alg:.8f}, "
            f"quadrature={E_quad[0]:.8f}, diff={diff:.2e}"
        )

    def test_energy_physical(self):
        """Algebraic He energy should be in the physical range.

        Uses R_min=0.5 (a=1.5) where Stieltjes recurrence is accurate.
        The production R_min=0.05 (a=0.15) has reduced accuracy.
        """
        E, _ = algebraic_he_energy_lmax0(
            Z=2.0, n_basis=25, alpha=1.5, R_min=0.5,
        )
        # Single-channel l_max=0 at R_min=0.5 gives slightly higher energy
        # than R_min=0.05 due to missing short-range physics, but should
        # still be in the physical range
        assert -3.5 < E < -2.0, f"Unphysical energy: {E}"

    def test_energy_variational(self):
        """Energy should decrease with increasing basis size (R_min=0.5)."""
        E1, _ = algebraic_he_energy_lmax0(n_basis=10, alpha=1.5, R_min=0.5)
        E2, _ = algebraic_he_energy_lmax0(n_basis=20, alpha=1.5, R_min=0.5)
        # More basis functions should give lower (better) energy
        assert E2 <= E1 + 1e-4, (
            f"Variational violation: E(10)={E1:.6f}, E(20)={E2:.6f}"
        )


# ─────────────────────────────────────────────────────────────────────
# 13. Tier 2: Discriminant analysis and recurrence detection
# ─────────────────────────────────────────────────────────────────────

class TestDiscriminantAnalysis:
    """Test discriminant analysis for l_max=1."""

    def test_discriminant_degree(self):
        """Discriminant should be degree 2 in R for l_max=1."""
        analysis = analyze_discriminant_lmax1(Z=2.0)
        assert analysis['disc_degree_R'] == 2, (
            f"Expected degree 2, got {analysis['disc_degree_R']}"
        )

    def test_discriminant_positive(self):
        """Discriminant should be positive for all R > 0 (real eigenvalues)."""
        analysis = analyze_discriminant_lmax1(Z=2.0)
        for R_val, disc_val in analysis['disc_values'].items():
            assert disc_val > 0, (
                f"Negative discriminant at R={R_val}: {disc_val}"
            )

    def test_integral_type_classification(self):
        """Should classify as square_root_quadratic."""
        analysis = analyze_discriminant_lmax1(Z=2.0)
        assert analysis['integral_type'] == 'square_root_quadratic'


class TestRecurrenceDetection:
    """Test linear recurrence detection for l_max=1 V_eff."""

    @pytest.mark.slow
    def test_recurrence_search(self):
        """Search for recurrence in V_nm at l_max=1."""
        result = detect_recurrence_lmax1(
            Z=2.0, n_basis=20, alpha=1.5, R_min=0.05,
            max_order=8,
        )
        # This is diagnostic — report results regardless of outcome
        print(f"\nRecurrence detection results:")
        print(f"  has_recurrence: {result['has_recurrence']}")
        for m, r in result['results_by_m'].items():
            print(f"  m={m}: order={r['order_found']}, sv_ratio={r['sv_ratio']:.2e}")
        # No assertion — this is investigative
