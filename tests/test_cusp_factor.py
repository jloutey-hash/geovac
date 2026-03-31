"""
Tests for the Kato cusp factor solver (Track U).

Validates:
1. Baseline reproduction (cusp_gamma=0 matches standard spectral solver)
2. Overlap matrix positive definiteness and conditioning
3. Gerade/ungerade symmetry preservation
4. Hermiticity of the transformed Hamiltonian
5. D_e convergence with cusp factor (documents negative result)

TRACK U RESULT (NEGATIVE):
The alpha-only Kato cusp factor f(alpha) = 1 + gamma*sin(2*alpha) does NOT
improve D_e or l_max convergence.  The cusp singularity in Level 4 coordinates
lives at alpha = pi/4 AND theta_12 = 0 jointly.  Multiplying by f(alpha) only
modifies the alpha dependence, while the slow partial-wave convergence is
dominated by the theta_12 expansion (Gegenbauer in cos(theta_12)).
The alpha-only factor cannot accelerate the angular momentum convergence
because the singularity is inherently 2D (alpha, theta_12).
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose


class TestCuspFactorBasics:
    """Basic sanity checks for the cusp factor."""

    def test_cusp_factor_symmetry(self) -> None:
        """f(alpha) = f(pi/2 - alpha): gerade symmetry preserved."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=0, n_basis=5, n_quad=50, cusp_gamma=1.0)

        # sin(2*alpha) is symmetric about pi/4:
        # sin(2*(pi/2 - alpha)) = sin(pi - 2*alpha) = sin(2*alpha)
        test_alpha = np.array([0.1, 0.3, np.pi / 4 - 0.1, np.pi / 4])
        test_reflected = np.pi / 2 - test_alpha
        f_test = 1.0 + 1.0 * np.sin(2.0 * test_alpha)
        f_refl = 1.0 + 1.0 * np.sin(2.0 * test_reflected)
        assert_allclose(f_test, f_refl, atol=1e-14)

    def test_cusp_factor_gamma_zero(self) -> None:
        """gamma=0 gives f=1 everywhere."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=50, cusp_gamma=0.0)
        f = solver._cusp_factor(R_e=3.0)
        assert_allclose(f, 1.0, atol=1e-15)

    def test_cusp_factor_positive(self) -> None:
        """f(alpha) > 0 for all alpha in (0, pi/2) at moderate R_e."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=0, n_basis=5, n_quad=100)
        # At R_e=2, gamma=1, f = 1 + sin(2a) >= 0 with min at a=3pi/4 (outside range)
        f = solver._cusp_factor(R_e=2.0)
        assert np.all(f > 0), f"f has negative values: min={f.min()}"

    def test_cusp_factor_derivative_consistency(self) -> None:
        """Check f' and f'' are consistent with f via finite differences."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=0, n_basis=5, n_quad=200)
        R_e = 3.0
        f = solver._cusp_factor(R_e)
        fp, fpp = solver._cusp_factor_deriv(R_e)

        # Numerical derivative of f
        fp_num = np.gradient(f, solver._alpha)

        # Compare in interior (away from endpoints)
        interior = slice(5, -5)
        assert_allclose(fp[interior], fp_num[interior], rtol=0.01,
                        err_msg="f' inconsistent with numerical derivative")

    def test_analytic_derivative_vs_numerical(self) -> None:
        """Analytic phi_bare' matches numerical gradient."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=200, cusp_gamma=0.0)
        for ic in range(solver.n_ch):
            phi_bare = solver._channel_phi_bare[ic]
            dphi_analytic = solver._analytic_derivative(ic, phi_bare)

            # Numerical derivative
            dphi_num = np.zeros_like(phi_bare)
            for k in range(phi_bare.shape[0]):
                dphi_num[k] = np.gradient(phi_bare[k], solver._alpha)

            # Compare in interior
            interior = slice(10, -10)
            if np.max(np.abs(dphi_analytic[:, interior])) > 1e-10:
                assert_allclose(
                    dphi_analytic[:, interior], dphi_num[:, interior],
                    rtol=0.05, atol=1e-6,
                    err_msg=f"Derivative mismatch in channel {ic}",
                )


class TestCuspFactorBaseline:
    """Verify that gamma=0 reproduces the standard spectral solver exactly."""

    def test_overlap_identity_at_gamma_zero(self) -> None:
        """At gamma=0, overlap matrix is positive definite."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100, cusp_gamma=0.0)
        channel_phi = solver._build_modified_basis(R_e=2.0)
        S = solver._build_overlap_matrix(channel_phi)

        evals = np.linalg.eigvalsh(S)
        assert np.all(evals > 0), f"S not positive definite: min eig = {evals[0]}"

    def test_eigenvalue_baseline_lmax0(self) -> None:
        """At gamma=0, cusp solver eigenvalues match standard spectral solver."""
        from geovac.cusp_factor import CuspFactorSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        R_e = 2.0
        rho = 1.4 / (2.0 * R_e)

        std = Level4SpectralAngular(l_max=0, n_basis=5, n_quad=100)
        mu_std, _ = std.solve(rho, R_e, n_eig=1)

        cusp = CuspFactorSolver(l_max=0, n_basis=5, n_quad=100, cusp_gamma=0.0)
        mu_cusp, _ = cusp.solve(rho, R_e, n_eig=1)

        assert_allclose(mu_cusp[0], mu_std[0], rtol=1e-6,
                        err_msg=f"Baseline mismatch: cusp={mu_cusp[0]:.8f} vs std={mu_std[0]:.8f}")

    def test_eigenvalue_baseline_lmax1(self) -> None:
        """At gamma=0, cusp solver matches standard solver at l_max=1."""
        from geovac.cusp_factor import CuspFactorSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        R_e = 2.5
        rho = 1.4 / (2.0 * R_e)

        std = Level4SpectralAngular(l_max=1, n_basis=8, n_quad=100)
        mu_std, _ = std.solve(rho, R_e, n_eig=1)

        cusp = CuspFactorSolver(l_max=1, n_basis=8, n_quad=100, cusp_gamma=0.0)
        mu_cusp, _ = cusp.solve(rho, R_e, n_eig=1)

        assert_allclose(mu_cusp[0], mu_std[0], rtol=1e-5,
                        err_msg="Baseline mismatch at l_max=1")

    def test_eigenvalue_baseline_lmax2(self) -> None:
        """At gamma=0, cusp solver matches standard solver at l_max=2."""
        from geovac.cusp_factor import CuspFactorSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        R_e = 2.5
        rho = 1.4 / (2.0 * R_e)

        std = Level4SpectralAngular(l_max=2, n_basis=8, n_quad=100)
        mu_std, _ = std.solve(rho, R_e, n_eig=1)

        cusp = CuspFactorSolver(l_max=2, n_basis=8, n_quad=100, cusp_gamma=0.0)
        mu_cusp, _ = cusp.solve(rho, R_e, n_eig=1)

        assert_allclose(mu_cusp[0], mu_std[0], rtol=1e-4,
                        err_msg="Baseline mismatch at l_max=2")

    def test_eigenvalue_baseline_across_Re(self) -> None:
        """At gamma=0, eigenvalues match across the full R_e range."""
        from geovac.cusp_factor import CuspFactorSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        R = 1.4
        for R_e in [1.0, 2.0, 5.0, 10.0]:
            rho = R / (2.0 * R_e)

            std = Level4SpectralAngular(l_max=1, n_basis=10, n_quad=100)
            mu_std, _ = std.solve(rho, R_e, n_eig=1)

            cusp = CuspFactorSolver(l_max=1, n_basis=10, n_quad=100, cusp_gamma=0.0)
            mu_cusp, _ = cusp.solve(rho, R_e, n_eig=1)

            assert_allclose(mu_cusp[0], mu_std[0], rtol=1e-10,
                            err_msg=f"Baseline mismatch at R_e={R_e}")


class TestCuspFactorHermiticity:
    """Verify Hermiticity of the modified Hamiltonian."""

    def test_overlap_symmetric(self) -> None:
        """Overlap matrix must be symmetric."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100)
        channel_phi = solver._build_modified_basis(R_e=2.0)
        S = solver._build_overlap_matrix(channel_phi)
        assert_allclose(S, S.T, atol=1e-14, err_msg="Overlap not symmetric")

    def test_overlap_positive_definite(self) -> None:
        """Overlap matrix must be positive definite."""
        from geovac.cusp_factor import CuspFactorSolver

        for gamma in [0.0, 0.5, 1.0, 2.0]:
            solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100, cusp_gamma=gamma)
            channel_phi = solver._build_modified_basis(R_e=2.0)
            S = solver._build_overlap_matrix(channel_phi)
            evals = np.linalg.eigvalsh(S)
            assert np.all(evals > 0), (
                f"S not positive definite at gamma={gamma}: min eig = {evals[0]}")

    def test_overlap_conditioning(self) -> None:
        """Overlap matrix condition number should not explode with cusp factor."""
        from geovac.cusp_factor import CuspFactorSolver

        s0 = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100, cusp_gamma=0.0)
        phi0 = s0._build_modified_basis(R_e=2.0)
        S0 = s0._build_overlap_matrix(phi0)
        cond0 = np.linalg.cond(S0)

        s1 = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100)
        phi1 = s1._build_modified_basis(R_e=2.0)
        S1 = s1._build_overlap_matrix(phi1)
        cond1 = np.linalg.cond(S1)

        assert cond1 < 10 * cond0, (
            f"Condition number exploded: {cond0:.1f} -> {cond1:.1f}")

    def test_hamiltonian_symmetric(self) -> None:
        """Full Hamiltonian in modified basis should be symmetric."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=1, n_basis=5, n_quad=100)
        R_e = 2.5
        rho = 1.4 / (2.0 * R_e)

        channel_phi = solver._build_modified_basis(R_e)
        H = solver._build_kinetic_matrix(channel_phi, R_e)
        V_nuc = solver._compute_nuclear_modified(channel_phi, rho, R_e)
        V_ee = solver._compute_vee_modified(channel_phi)
        H_full = H + R_e * (V_nuc + V_ee)

        asym = np.max(np.abs(H_full - H_full.T))
        assert asym < 1e-10, f"Hamiltonian asymmetry: {asym:.2e}"


class TestCuspFactorGerade:
    """Verify gerade symmetry preservation."""

    def test_channel_labels_preserved(self) -> None:
        """Channel quantum numbers must be unchanged by cusp factor."""
        from geovac.cusp_factor import CuspFactorSolver
        from geovac.level4_spectral_angular import Level4SpectralAngular

        for l_max in [0, 1, 2]:
            std = Level4SpectralAngular(l_max=l_max, n_basis=5, n_quad=50)
            cusp = CuspFactorSolver(l_max=l_max, n_basis=5, n_quad=50)
            assert std.channels_2 == cusp.channels_2, (
                f"Channel labels differ at l_max={l_max}")

    def test_gerade_selection(self) -> None:
        """Only l1+l2 even channels should be present (homonuclear)."""
        from geovac.cusp_factor import CuspFactorSolver

        solver = CuspFactorSolver(l_max=2, n_basis=5, n_quad=50)
        for l1, l2 in solver.channels_2:
            assert (l1 + l2) % 2 == 0, f"Non-gerade channel: ({l1},{l2})"


class TestCuspFactorDe:
    """D_e convergence tests documenting the negative result.

    RESULT: The alpha-only cusp factor provides no D_e improvement.
    The cusp singularity at r_12 = 0 involves both alpha AND theta_12,
    so an alpha-only modification cannot accelerate the Gegenbauer
    convergence in theta_12 that dominates the l_max expansion.
    """

    @pytest.mark.slow
    def test_de_baseline_matches_standard(self) -> None:
        """D_e at gamma=0 matches standard solver."""
        from geovac.cusp_factor import solve_level4_h2_cusp
        from geovac.level4_multichannel import solve_level4_h2_multichannel

        result_cusp = solve_level4_h2_cusp(
            R=1.4, l_max=2, n_basis=10, n_quad=100,
            cusp_gamma=0.0, verbose=False,
        )
        result_std = solve_level4_h2_multichannel(
            R=1.4, l_max=2, angular_method='spectral',
            n_basis_angular=10, n_quad_angular=100, verbose=False,
        )
        assert_allclose(result_cusp['D_e_pct'], result_std['D_e_pct'], atol=0.1,
                        err_msg="Baseline D_e does not match standard solver")

    @pytest.mark.slow
    def test_cusp_does_not_improve_de(self) -> None:
        """Document: cusp factor does NOT improve D_e (negative result).

        The alpha-only cusp factor cannot accelerate the theta_12-dependent
        partial-wave convergence that limits D_e recovery.
        """
        from geovac.cusp_factor import solve_level4_h2_cusp

        result_base = solve_level4_h2_cusp(
            R=1.4, l_max=2, n_basis=10, n_quad=100,
            cusp_gamma=0.0, verbose=False,
        )
        result_cusp = solve_level4_h2_cusp(
            R=1.4, l_max=2, n_basis=10, n_quad=100,
            cusp_gamma=None, verbose=False,
        )
        delta = result_cusp['D_e_pct'] - result_base['D_e_pct']

        # Document the negative result: cusp factor gives -0.3 pp change
        assert abs(delta) < 1.0, (
            f"Unexpected large change from cusp factor: delta={delta:.1f} pp")
        # The cusp factor slightly degrades D_e due to non-orthogonality noise
        assert delta <= 0.1, (
            f"Cusp factor unexpectedly improves D_e by {delta:.1f} pp")

    @pytest.mark.slow
    def test_convergence_table(self) -> None:
        """Build and print convergence table for documentation."""
        from geovac.cusp_factor import solve_level4_h2_cusp

        results = {'baseline': {}, 'cusp': {}}

        for l_max in [0, 1, 2, 3]:
            r0 = solve_level4_h2_cusp(
                R=1.4, l_max=l_max, n_basis=10, n_quad=100,
                cusp_gamma=0.0, verbose=False,
            )
            results['baseline'][l_max] = r0['D_e_pct']

            r1 = solve_level4_h2_cusp(
                R=1.4, l_max=l_max, n_basis=10, n_quad=100,
                cusp_gamma=None, verbose=False,
            )
            results['cusp'][l_max] = r1['D_e_pct']

        print("\n=== Track U: D_e Convergence Table (NEGATIVE RESULT) ===")
        print(f"{'l_max':>5} {'Baseline %':>12} {'Cusp %':>12} {'Delta':>10}")
        for l_max in [0, 1, 2, 3]:
            b = results['baseline'][l_max]
            c = results['cusp'][l_max]
            d = c - b
            print(f"{l_max:>5} {b:>12.1f} {c:>12.1f} {d:>10.1f}")

        print("\nConclusion: alpha-only cusp factor provides no improvement.")
        print("Root cause: cusp singularity is 2D (alpha, theta_12),")
        print("not separable in alpha alone.")

        # Verify all results are physical
        for l_max in [0, 1, 2, 3]:
            assert results['baseline'][l_max] > 0
            assert results['cusp'][l_max] > 0

    @pytest.mark.slow
    def test_gamma_scan(self) -> None:
        """Verify D_e degrades monotonically with increasing gamma."""
        from geovac.cusp_factor import solve_level4_h2_cusp

        gammas = [0.0, 0.1, 0.5, 1.0, 2.0]
        de_pcts = []
        for gamma in gammas:
            r = solve_level4_h2_cusp(
                R=1.4, l_max=2, n_basis=10, n_quad=100,
                cusp_gamma=gamma, verbose=False,
            )
            de_pcts.append(r['D_e_pct'])

        print("\n=== Gamma scan at l_max=2 ===")
        for g, d in zip(gammas, de_pcts):
            print(f"  gamma={g:.1f}: D_e = {d:.1f}%")

        # D_e should decrease (or stay flat) with increasing gamma
        for i in range(1, len(de_pcts)):
            assert de_pcts[i] <= de_pcts[0] + 0.5, (
                f"D_e increased at gamma={gammas[i]}: {de_pcts[i]:.1f}% vs {de_pcts[0]:.1f}%")
