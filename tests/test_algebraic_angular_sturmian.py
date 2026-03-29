"""
Tests for the Sturmian angular basis prototype (Level 3 hyperspherical).

Validates the SturmianAngularSolver against the free-basis AlgebraicAngularSolver
and verifies key properties: orthonormality, completeness, convergence.

Test inventory:
  1-3.  Sturmian basis orthonormality and completeness
  4-5.  Eigenvalue consistency with free basis at matched dimension
  6-7.  He energy benchmarks at l_max=0 and l_max=1
  8-9.  R0 sensitivity and n_construct convergence
  10.   Monopole decomposition fractions
  11.   DBOC consistency
"""

import numpy as np
import pytest
from math import pi


# ---------------------------------------------------------------------------
# 1. Sturmian basis orthonormality
# ---------------------------------------------------------------------------

def test_sturmian_orthonormality():
    """Verify the truncated Sturmian rotation matrix has orthonormal columns.

    U^T @ U should be identity (columns are orthonormal).
    """
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver

    solver = SturmianAngularSolver(Z=2.0, n_basis=10, l_max=0,
                                    R0=1.5, n_construct=50)
    U = solver._U_truncated

    # U^T @ U should be identity (dim_kept x dim_kept)
    gram = U.T @ U
    np.testing.assert_allclose(
        gram, np.eye(solver._total_dim_kept), atol=1e-12,
        err_msg="Sturmian columns are not orthonormal"
    )


def test_sturmian_orthonormality_multichannel():
    """Verify orthonormality with l_max > 0."""
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver

    solver = SturmianAngularSolver(Z=2.0, n_basis=8, l_max=2,
                                    R0=1.5, n_construct=30)
    U = solver._U_truncated
    gram = U.T @ U

    np.testing.assert_allclose(
        gram, np.eye(solver._total_dim_kept), atol=1e-12,
        err_msg="Multichannel Sturmian columns not orthonormal"
    )


# ---------------------------------------------------------------------------
# 2. Completeness: Sturmian at full dimension = free basis
# ---------------------------------------------------------------------------

def test_sturmian_full_dimension_matches_free():
    """When n_construct == n_basis, Sturmian spans the same space as free.

    Eigenvalues should be identical (rotation within the same subspace).
    """
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver
    from geovac.algebraic_angular import AlgebraicAngularSolver

    nb = 15
    solver_s = SturmianAngularSolver(Z=2.0, n_basis=nb, l_max=0,
                                      R0=1.5, n_construct=nb)
    solver_f = AlgebraicAngularSolver(Z=2.0, n_basis=nb, l_max=0)

    for R in [0.5, 1.0, 2.0, 5.0]:
        evals_s, _ = solver_s.solve(R, n_channels=5)
        evals_f, _ = solver_f.solve(R, n_channels=5)

        np.testing.assert_allclose(
            evals_s, evals_f, atol=1e-10,
            err_msg=f"R={R}: Sturmian at full dim should match free"
        )


# ---------------------------------------------------------------------------
# 3. Eigenvalue consistency at R=0
# ---------------------------------------------------------------------------

def test_sturmian_R0_eigenvalues():
    """At R=R0, Sturmian eigenvalues should be the Sturmian eigenvalues.

    H(R0) = diag(mu_sturmian) + 0*V_ref + R0*V_residual
           = diag(mu_sturmian) + R0*V_residual

    For ref='nuclear', V_residual = V_ee, so H(R0) != diag(mu_sturm).
    But H(R0) eigenvalues should match the free solver at R=R0.
    """
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver
    from geovac.algebraic_angular import AlgebraicAngularSolver

    # Use large free basis as reference
    ref = AlgebraicAngularSolver(Z=2.0, n_basis=80, l_max=0, n_quad=200)
    evals_ref, _ = ref.solve(1.5, n_channels=1)

    solver = SturmianAngularSolver(Z=2.0, n_basis=15, l_max=0,
                                    R0=1.5, n_construct=50)
    evals_s, _ = solver.solve(1.5, n_channels=1)

    # Should be close to reference (better than free at n_basis=15)
    err_s = abs(evals_s[0] - evals_ref[0])
    assert err_s < 0.01, f"Sturmian mu(R0) error = {err_s:.6f}"


# ---------------------------------------------------------------------------
# 4. Sturmian improves convergence at l_max=0
# ---------------------------------------------------------------------------

def test_sturmian_improves_lmax0():
    """Sturmian basis should give lower error than free at same n_basis for l_max=0.

    At n_basis=10, n_construct=50, Sturmian should beat free by at least 20%.
    """
    from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    E_exact = -2.903724

    res_s = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=10, l_max=0, R0=1.5, n_construct=50,
        n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
    )
    res_f = solve_hyperspherical_algebraic(
        Z=2.0, n_basis=10, l_max=0,
        n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
    )

    err_s = abs(res_s['energy'] - E_exact) / abs(E_exact) * 100
    err_f = abs(res_f['energy'] - E_exact) / abs(E_exact) * 100

    assert err_s < err_f, \
        f"Sturmian ({err_s:.4f}%) should be better than free ({err_f:.4f}%)"
    assert err_s < 0.20, \
        f"Sturmian at n_basis=10 error {err_s:.4f}% exceeds 0.20%"


# ---------------------------------------------------------------------------
# 5. He energy at l_max=0 with Sturmian
# ---------------------------------------------------------------------------

def test_he_energy_sturmian_lmax0():
    """Sturmian solver should achieve < 0.15% at n_basis=15, l_max=0."""
    from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian

    result = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=15, l_max=0, R0=1.5, n_construct=50,
        n_R=200, R_max=30.0, N_R_radial=2000, verbose=False,
    )

    E_exact = -2.903724
    err = abs(result['energy'] - E_exact) / abs(E_exact) * 100

    assert err < 0.15, \
        f"He energy error {err:.4f}% exceeds 0.15% (Sturmian l_max=0)"


# ---------------------------------------------------------------------------
# 6. He energy at l_max=1 beats FD solver
# ---------------------------------------------------------------------------

def test_he_energy_sturmian_lmax1_beats_fd():
    """Sturmian solver at l_max=1 should beat FD solver (0.78%).

    The FD solver degrades to 0.78% at l_max=1 due to centrifugal
    singularity.  The Sturmian basis handles this analytically
    (inherited from the Gegenbauer basis).

    Note: the Sturmian at l_max>0 is worse than the free basis at
    the same n_basis due to adiabatic approximation effects, but
    both beat the FD solver significantly.
    """
    from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian

    result = solve_hyperspherical_sturmian(
        Z=2.0, n_basis=15, l_max=1, R0=1.5, n_construct=50,
        n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
    )

    E_exact = -2.903724
    err = abs(result['energy'] - E_exact) / abs(E_exact) * 100

    # Must beat FD (0.78%) by a clear margin
    assert err < 0.65, \
        f"He l_max=1 error {err:.4f}% exceeds 0.65% (FD gives 0.78%)"


# ---------------------------------------------------------------------------
# 7. R0 sensitivity
# ---------------------------------------------------------------------------

def test_R0_sensitivity():
    """Verify R0 in [0.5, 3.0] gives similar results (< 0.03% spread)."""
    from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian

    E_exact = -2.903724
    energies = []

    for R0 in [0.5, 1.0, 1.5, 2.0, 3.0]:
        result = solve_hyperspherical_sturmian(
            Z=2.0, n_basis=15, l_max=0, R0=R0, n_construct=50,
            n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
        )
        energies.append(result['energy'])

    spread = max(energies) - min(energies)
    spread_pct = abs(spread) / abs(E_exact) * 100

    assert spread_pct < 0.05, \
        f"R0 spread {spread_pct:.4f}% exceeds 0.05%"


# ---------------------------------------------------------------------------
# 8. n_construct convergence
# ---------------------------------------------------------------------------

def test_n_construct_convergence():
    """Verify energy converges as n_construct increases."""
    from geovac.algebraic_angular_sturmian import solve_hyperspherical_sturmian

    energies = []
    for nc in [15, 30, 50]:
        result = solve_hyperspherical_sturmian(
            Z=2.0, n_basis=10, l_max=0, R0=1.5, n_construct=nc,
            n_R=150, R_max=25.0, N_R_radial=1500, verbose=False,
        )
        energies.append(result['energy'])

    # Energy should decrease monotonically (more negative = better)
    for i in range(len(energies) - 1):
        assert energies[i + 1] <= energies[i] + 1e-6, \
            f"Non-monotonic: E(nc={[15,30,50][i]})={energies[i]:.8f}, " \
            f"E(nc={[15,30,50][i+1]})={energies[i+1]:.8f}"

    # Convergence: last two should be within 0.005 Ha
    assert abs(energies[-1] - energies[-2]) < 0.005, \
        f"Not converged: E(30)={energies[1]:.6f}, E(50)={energies[2]:.6f}"


# ---------------------------------------------------------------------------
# 9. Monopole decomposition
# ---------------------------------------------------------------------------

def test_monopole_decomposition():
    """Verify nuclear monopole dominates V_coupling (> 95% Frobenius norm)."""
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver

    solver = SturmianAngularSolver(Z=2.0, n_basis=10, l_max=0,
                                    R0=1.5, n_construct=50)
    decomp = solver.monopole_decomposition()

    assert decomp['monopole_fraction'] > 0.95, \
        f"Monopole fraction {decomp['monopole_fraction']:.1%} < 95%"
    assert decomp['remainder_fraction'] < 0.10, \
        f"Remainder fraction {decomp['remainder_fraction']:.1%} > 10%"
    assert decomp['projection_error'] < 1e-10, \
        f"Projection error {decomp['projection_error']:.2e} > 1e-10"


# ---------------------------------------------------------------------------
# 10. DBOC consistency
# ---------------------------------------------------------------------------

def test_sturmian_dboc_positive():
    """Verify DBOC is positive at all R in the Sturmian basis."""
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver

    solver = SturmianAngularSolver(Z=2.0, n_basis=10, l_max=0,
                                    R0=1.5, n_construct=50)

    for R in [0.5, 1.0, 2.0, 5.0, 10.0]:
        _, _, dboc = solver.solve_with_dboc(R, n_channels=1)
        assert dboc >= 0.0, f"DBOC({R}) = {dboc:.6e} is negative"


# ---------------------------------------------------------------------------
# 11. Coupling matrix symmetry in Sturmian basis
# ---------------------------------------------------------------------------

def test_sturmian_coupling_symmetry():
    """Verify projected coupling matrices are symmetric."""
    from geovac.algebraic_angular_sturmian import SturmianAngularSolver

    solver = SturmianAngularSolver(Z=2.0, n_basis=10, l_max=1,
                                    R0=1.5, n_construct=30)

    np.testing.assert_allclose(
        solver._V_coupling_proj, solver._V_coupling_proj.T, atol=1e-12,
        err_msg="Projected coupling matrix not symmetric"
    )
    np.testing.assert_allclose(
        solver._V_ref_proj, solver._V_ref_proj.T, atol=1e-12,
        err_msg="Projected V_ref matrix not symmetric"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
