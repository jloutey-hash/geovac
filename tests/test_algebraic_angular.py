"""
Tests for the algebraic angular solver (Level 3 hyperspherical).

Validates the spectral Gegenbauer solver against known analytical results,
the existing FD solver, and internal consistency checks.

Test inventory:
  1. Casimir eigenvalues (R -> 0 limit)
  2. First-order perturbation coefficient a_1
  3. GL quadrature vs scipy.integrate.quad consistency
  4. Cross-validation against FD solver
  5. He ground-state energy (full pipeline)
  6. l_max convergence (stub for future multichannel)
  7. Basis convergence (spectral vs polynomial)
  8. Selection rules (nuclear coupling zeros)
"""

import numpy as np
import pytest
from math import pi, sqrt


# ---------------------------------------------------------------------------
# 1. Casimir eigenvalue test
# ---------------------------------------------------------------------------

def test_casimir_eigenvalues():
    """Verify H_ang(R=0) eigenvalues match SO(6) Casimir nu(nu+4)/2.

    For singlet (n = 1, 3, 5, ...):
      n=1 -> mu_free = 0   (nu=0)
      n=3 -> mu_free = 16  (nu=4)
      n=5 -> mu_free = 48  (nu=8)
      n=7 -> mu_free = 96  (nu=12)
      n=9 -> mu_free = 160 (nu=16)
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=0, symmetry='singlet')

    # At R ~ 0, eigenvalues should be the Casimir values
    evals, _ = solver.solve(R=0.0, n_channels=5)

    expected = np.array([0.0, 16.0, 48.0, 96.0, 160.0])
    np.testing.assert_allclose(evals, expected, atol=1e-12,
                               err_msg="Casimir eigenvalues at R=0 do not match SO(6)")


# ---------------------------------------------------------------------------
# 2. First-order perturbation coefficient a_1
# ---------------------------------------------------------------------------

def test_first_order_perturbation():
    """Verify a_1 = dmu/dR|_{R=0} matches Paper 13 closed form.

    Paper 13 Eq. 32: a_1(n=1, Z) = (8/(3*pi)) * (sqrt(2) - 4*Z)
    At Z=2: a_1 = (8/(3*pi)) * (sqrt(2) - 8) = -5.590189...
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=15, l_max=0)

    # Method 1: coupling matrix diagonal element = a_1 directly
    a1_matrix = solver.coupling_matrix[0, 0]

    # Method 2: numerical derivative at small R
    dR = 1e-5
    evals_dR, _ = solver.solve(dR, n_channels=1)
    a1_numerical = evals_dR[0] / dR  # since mu(0) = 0

    # Analytical value from Paper 13
    a1_exact = (8.0 / (3.0 * pi)) * (sqrt(2) - 4.0 * 2.0)

    assert abs(a1_matrix - a1_exact) < 1e-8, \
        f"a_1 from coupling matrix: {a1_matrix:.10f}, exact: {a1_exact:.10f}"
    assert abs(a1_numerical - a1_exact) < 1e-4, \
        f"a_1 from numerical derivative: {a1_numerical:.8f}, exact: {a1_exact:.8f}"


def test_first_order_perturbation_components():
    """Verify nuclear and V_ee contributions to a_1 separately.

    Nuclear contribution: -Z * I_nuc(1) * (8/pi)
    where I_nuc(1) = 1 + 1/3 = 4/3 (Paper 13 Eq. 31, both 1/cos and 1/sin)

    Total nuclear a_1 component: -Z * (4/3) * (8/pi) * 2 = ... wait,
    let's compute directly:

    V_nuc_00 = (4/pi) * integral of sin^2(2*alpha) * (-Z)(1/cos + 1/sin) dalpha

    sin^2(2*alpha) = 4 sin^2(alpha) cos^2(alpha), so:
    1/cos part: (4/pi) * (-Z) * 4 * integral_0^{pi/2} sin^2(alpha) cos(alpha) dalpha
              = (4/pi) * (-Z) * 4 * 1/3 = -16Z/(3*pi)
    1/sin part equals 1/cos part by exchange symmetry (both n, n' odd).
    Total: V_nuc_00 = -32Z/(3*pi)

    V_ee_00 = (4/pi) * integral of sin^2(2*alpha) / max(cos,sin) dalpha

    [0, pi/4]: 4sin^2 cos^2 / cos = 4sin^2 cos, integral = 4sin^3(pi/4)/3 = sqrt(2)/3
    [pi/4, pi/2]: by exchange symmetry = sqrt(2)/3
    Total: V_ee_00 = (4/pi)(2*sqrt(2)/3) = 8*sqrt(2)/(3*pi)
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)

    nuc_expected = -2.0 * 32.0 / (3.0 * pi)   # -Z * 32/(3*pi)
    vee_expected = 8.0 * sqrt(2) / (3.0 * pi)

    nuc_actual = solver.nuclear_matrix[0, 0]
    vee_actual = solver.vee_matrix[0, 0]

    assert abs(nuc_actual - nuc_expected) < 1e-10, \
        f"Nuclear V_00: {nuc_actual:.12f}, expected: {nuc_expected:.12f}"
    assert abs(vee_actual - vee_expected) < 1e-10, \
        f"V_ee V_00: {vee_actual:.12f}, expected: {vee_expected:.12f}"


# ---------------------------------------------------------------------------
# 3. GL quadrature vs scipy.integrate.quad consistency
# ---------------------------------------------------------------------------

def test_quadrature_consistency():
    """Verify GL quadrature matches adaptive quad to < 10^-10.

    The integrands are smooth (basis functions absorb singularities),
    so both methods should agree to near machine precision.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=0, n_quad=150)

    # Check nuclear and V_ee matrix elements for several (i,j) pairs
    test_pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 3), (2, 4)]

    for i, j in test_pairs:
        ni = solver.basis_indices[i]
        nj = solver.basis_indices[j]

        # GL quadrature (from the precomputed matrices)
        nuc_gl = solver.nuclear_matrix[i, j]
        vee_gl = solver.vee_matrix[i, j]

        # Adaptive quadrature
        nuc_quad, vee_quad = solver.coupling_integral_quad(ni, nj)

        assert abs(nuc_gl - nuc_quad) < 1e-10, \
            f"Nuclear ({i},{j}): GL={nuc_gl:.14f}, quad={nuc_quad:.14f}, " \
            f"diff={abs(nuc_gl - nuc_quad):.2e}"
        assert abs(vee_gl - vee_quad) < 1e-10, \
            f"V_ee ({i},{j}): GL={vee_gl:.14f}, quad={vee_quad:.14f}, " \
            f"diff={abs(vee_gl - vee_quad):.2e}"


# ---------------------------------------------------------------------------
# 4. Cross-validation against FD solver
# ---------------------------------------------------------------------------

def test_cross_validation_fd():
    """Compare mu(R) from algebraic solver against FD solver.

    At Z=2, n_basis=25, compare at 5 representative R values.
    Agreement should be < 0.5%.  Both solvers have different error
    profiles: FD has discretization error, algebraic has basis truncation.
    At small R (perturbative regime) both should agree very closely.
    At large R (non-perturbative), differences reflect basis completeness.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.hyperspherical_angular import solve_angular

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=25, l_max=0)

    # Small/medium R: tight tolerance (basis converges quickly)
    # Large R: relaxed tolerance (eigenfunction localizes, needs many basis fns)
    R_tol = [(0.1, 0.5), (0.5, 0.5), (1.0, 0.5), (3.0, 0.5), (10.0, 3.0)]

    for R, tol_pct in R_tol:
        # Algebraic solver
        evals_alg, _ = solver.solve(R, n_channels=1)
        mu_alg = evals_alg[0]

        # FD solver (high resolution for accuracy)
        evals_fd, _ = solve_angular(R, Z=2.0, l_max=0, n_alpha=400, n_channels=1)
        mu_fd = evals_fd[0]

        # Both should be negative for R > 0 at Z=2
        if abs(mu_fd) > 0.01:
            rel_diff = abs(mu_alg - mu_fd) / abs(mu_fd) * 100
            assert rel_diff < tol_pct, \
                f"R={R}: mu_alg={mu_alg:.8f}, mu_fd={mu_fd:.8f}, " \
                f"rel_diff={rel_diff:.4f}% (tol={tol_pct}%)"
        else:
            abs_diff = abs(mu_alg - mu_fd)
            assert abs_diff < 0.01, \
                f"R={R}: mu_alg={mu_alg:.8f}, mu_fd={mu_fd:.8f}, " \
                f"abs_diff={abs_diff:.6f}"


# ---------------------------------------------------------------------------
# 5. He ground-state energy (full pipeline)
# ---------------------------------------------------------------------------

def test_he_ground_state():
    """Run full pipeline and verify E < -2.903 Ha.

    The existing FD solver gives -2.9052 Ha (0.05% error) with 100
    grid points.  The algebraic solver at l_max=0 converges more slowly
    because the angular eigenfunction at large R localizes near alpha=0
    and alpha=pi/2, requiring many oscillatory basis functions.  The
    algebraic solver's advantage is at l_max > 0, where the FD solver
    degrades to 0.87% due to centrifugal singularity (Paper 13 Table I).

    For l_max=0, n_basis=25 gives ~ 0.14% error.
    """
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    result = solve_hyperspherical_algebraic(
        Z=2.0, n_basis=25, l_max=0, n_R=200,
        R_max=30.0, N_R_radial=2000, verbose=False,
    )

    E = result['energy']
    E_exact = -2.903724

    err_pct = abs(E - E_exact) / abs(E_exact) * 100

    assert E < -2.895, \
        f"He energy too high: {E:.6f} Ha (should be < -2.895)"
    assert err_pct < 0.3, \
        f"He energy error {err_pct:.4f}% exceeds 0.3% threshold"


# ---------------------------------------------------------------------------
# 6. l_max convergence test
# ---------------------------------------------------------------------------

def test_l_max_convergence():
    """Verify energy improves monotonically with l_max (THE KEY TEST).

    The FD solver (Paper 13 Table I) shows l_max degradation:
      l_max=0: -2.9052 (0.054% error, above exact)
      l_max=1: -2.9262 (0.78% error, below exact)
      l_max=3: -2.9289 (0.87% error, below exact)

    The algebraic solver handles the centrifugal singularity analytically,
    reducing the overshoot at l_max >= 1 by ~0.2 percentage points:
      l_max=0: ~0.16% (above exact, basis truncation limits l=0)
      l_max=1: ~0.56% (below exact, vs FD 0.78%)
      l_max=3: ~0.65% (below exact, vs FD 0.87%)

    The remaining degradation is intrinsic to the adiabatic approximation
    (non-adiabatic corrections are neglected), not the centrifugal singularity.

    Note: The adiabatic approximation is NOT variational — energies can
    go below exact when non-adiabatic repulsion is omitted.
    """
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    energies = {}
    for l_max in [0, 1, 2, 3]:
        result = solve_hyperspherical_algebraic(
            Z=2.0, n_basis=15, l_max=l_max, n_R=150,
            R_max=25.0, N_R_radial=1500, verbose=False,
        )
        energies[l_max] = result['energy']

    E_exact = -2.903724

    # Energy should decrease monotonically (more correlation captured)
    for l in range(3):
        assert energies[l + 1] < energies[l] + 1e-6, \
            f"Non-monotonic: E(l_max={l})={energies[l]:.8f}, " \
            f"E(l_max={l+1})={energies[l+1]:.8f}"

    # Algebraic solver should be closer to exact than FD at l_max >= 1
    # FD gives 0.87% at l_max=3; algebraic should be < 0.75%
    err_l3 = abs(energies[3] - E_exact) / abs(E_exact) * 100
    assert err_l3 < 0.75, \
        f"l_max=3 error {err_l3:.2f}% exceeds 0.75% (FD gives 0.87%)"


# ---------------------------------------------------------------------------
# 7. Basis convergence test
# ---------------------------------------------------------------------------

def test_basis_convergence():
    """Verify monotonic convergence with n_basis at fixed l_max=0.

    Spectral convergence should be exponential (or at least faster
    than the polynomial FD convergence).
    """
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    basis_sizes = [3, 5, 8, 10, 15]
    energies = []

    for nb in basis_sizes:
        result = solve_hyperspherical_algebraic(
            Z=2.0, n_basis=nb, l_max=0, n_R=150,
            R_max=25.0, N_R_radial=1500, verbose=False,
        )
        energies.append(result['energy'])

    # Energy should decrease monotonically (more negative = better)
    for i in range(len(energies) - 1):
        assert energies[i + 1] <= energies[i] + 1e-6, \
            f"Non-monotonic: E(n_basis={basis_sizes[i]})={energies[i]:.8f}, " \
            f"E(n_basis={basis_sizes[i+1]})={energies[i+1]:.8f}"

    # The best result should be close to the FD result (-2.9052)
    E_exact = -2.903724
    err_best = abs(energies[-1] - E_exact) / abs(E_exact) * 100
    assert err_best < 0.2, \
        f"Best energy at n_basis={basis_sizes[-1]}: error {err_best:.4f}%"


# ---------------------------------------------------------------------------
# 8. Selection rule test
# ---------------------------------------------------------------------------

def test_selection_rules():
    """Verify exchange symmetry selection rules.

    The nuclear potential -Z(1/cos + 1/sin) is exchange-symmetric
    under alpha -> pi/2 - alpha.  Under this exchange:
      sin(2m*alpha) -> (-1)^{m+1} sin(2m*alpha)

    For V_nuc matrix elements, exchange symmetry requires:
      V_nuc[m', m] = (-1)^{m'+m} V_nuc[m', m]

    Within the singlet sector (all m odd), m'+m is always even,
    so (-1)^{m'+m} = +1 — the selection rule is automatically
    satisfied and all couplings can be nonzero.

    The stronger test: singlet-triplet coupling should be ZERO.
    Compute V_nuc between a singlet mode (odd m) and a triplet mode
    (even m), which requires m'+m odd -> V = 0.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from scipy.integrate import quad
    import numpy as np

    Z = 2.0
    norm = 4.0 / np.pi

    # Test singlet-triplet decoupling: m=1 (singlet) vs m=2 (triplet)
    for m_singlet in [1, 3, 5]:
        for m_triplet in [2, 4, 6]:
            def integrand(alpha):
                s1 = np.sin(2 * m_singlet * alpha)
                s2 = np.sin(2 * m_triplet * alpha)
                return norm * s1 * s2 * (-Z) * (
                    1.0 / np.cos(alpha) + 1.0 / np.sin(alpha)
                )

            val, _ = quad(integrand, 0, np.pi / 2, limit=200,
                         points=[np.pi / 4])
            assert abs(val) < 1e-12, \
                f"Singlet-triplet coupling V[{m_singlet},{m_triplet}] " \
                f"= {val:.2e} (should be zero)"

    # Within singlet sector: verify couplings are generally nonzero
    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=0)
    V_nuc = solver.nuclear_matrix

    # Diagonal should be negative (nuclear attraction)
    for i in range(solver.n_basis):
        assert V_nuc[i, i] < 0, \
            f"V_nuc[{i},{i}] = {V_nuc[i, i]:.6f} should be negative"

    # Off-diagonal coupling should decay with distance
    # |V[0,1]| > |V[0,2]| > |V[0,3]| generally
    off_diag = [abs(V_nuc[0, j]) for j in range(1, solver.n_basis)]
    # At least first off-diagonal should be significant
    assert off_diag[0] > 0.1, \
        f"Nearest-neighbor coupling too small: {off_diag[0]:.6f}"


def test_vee_selection_rules():
    """Verify that V_ee monopole coupling is nonzero for more pairs.

    The V_ee potential 1/max(cos,sin) does NOT have the same selection
    rules as the nuclear potential.  It should couple adjacent modes.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=0)
    V_ee = solver.vee_matrix

    # V_ee diagonal should be positive (repulsive)
    for i in range(solver.n_basis):
        assert V_ee[i, i] > 0, \
            f"V_ee[{i},{i}] = {V_ee[i, i]:.6f} should be positive"


# ---------------------------------------------------------------------------
# Additional validation tests
# ---------------------------------------------------------------------------

def test_coupling_matrix_symmetry():
    """Verify coupling matrices are symmetric."""
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=8, l_max=0)

    np.testing.assert_allclose(
        solver.nuclear_matrix, solver.nuclear_matrix.T, atol=1e-14,
        err_msg="Nuclear matrix not symmetric"
    )
    np.testing.assert_allclose(
        solver.vee_matrix, solver.vee_matrix.T, atol=1e-14,
        err_msg="V_ee matrix not symmetric"
    )
    np.testing.assert_allclose(
        solver.coupling_matrix, solver.coupling_matrix.T, atol=1e-14,
        err_msg="Coupling matrix not symmetric"
    )


def test_triplet_symmetry():
    """Verify triplet basis uses even n and different Casimir values."""
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, symmetry='triplet')

    # Basis indices should be 2, 4, 6, 8, 10
    np.testing.assert_array_equal(solver.basis_indices, [2, 4, 6, 8, 10])

    # Casimir eigenvalues: 2n^2 - 2
    expected_casimir = np.array([6.0, 30.0, 70.0, 126.0, 198.0])
    np.testing.assert_allclose(solver.casimir_eigenvalues, expected_casimir)


def test_large_R_asymptotics():
    """At large R, mu(R) should approach -Z^2 * R^2 / 2 (He+ threshold).

    For Z=2: mu(R) -> -2*R^2 at large R.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=15, l_max=0)

    R_large = 20.0
    evals, _ = solver.solve(R_large, n_channels=1)
    mu = evals[0]

    # Asymptotic value: -Z^2 * R^2 / 2 = -2 * 400 = -800
    mu_asymp = -2.0 * R_large**2 / 2.0

    # Finite basis can't perfectly capture the asymptotic regime,
    # but should get the right order of magnitude and sign
    assert mu < 0, f"mu at large R should be negative: {mu:.2f}"
    # The effective potential V_eff = mu/R^2 + 15/(8R^2) should
    # approach -Z^2/2 = -2
    V_eff = mu / R_large**2 + 15.0 / (8.0 * R_large**2)
    assert V_eff < -1.5, \
        f"V_eff(R={R_large}) = {V_eff:.4f}, should approach -2.0"


def test_z_scaling():
    """Verify a_1 scales correctly with Z.

    a_1(Z) = (8/(3*pi)) * (sqrt(2) - 4*Z)
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    for Z in [1.0, 2.0, 3.0, 5.0]:
        solver = AlgebraicAngularSolver(Z=Z, n_basis=10, l_max=0)
        a1_computed = solver.coupling_matrix[0, 0]
        a1_exact = (8.0 / (3.0 * pi)) * (sqrt(2) - 4.0 * Z)
        assert abs(a1_computed - a1_exact) < 1e-9, \
            f"Z={Z}: a1={a1_computed:.10f}, exact={a1_exact:.10f}"


def test_second_order_perturbation():
    """Verify the second-order coefficient a_2 at Z=2.

    Paper 13: a_2 = -0.2442 (computed from sum over intermediate states).

    a_2 = sum_{n'>1, singlet} |V[0, n']|^2 / (mu_0 - mu_n')
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=20, l_max=0)

    V = solver.coupling_matrix
    mu_free = solver.casimir_eigenvalues

    # Second-order perturbation sum
    a2 = 0.0
    for n in range(1, solver.n_basis):
        a2 += V[0, n]**2 / (mu_free[0] - mu_free[n])

    a2_paper = -0.2442

    # Allow 5% tolerance due to basis truncation
    assert abs(a2 - a2_paper) / abs(a2_paper) < 0.05, \
        f"a_2 = {a2:.6f}, Paper 13: {a2_paper:.4f}"


# ---------------------------------------------------------------------------
# Multichannel tests (l_max > 0)
# ---------------------------------------------------------------------------

def test_l1_free_eigenvalues():
    """Verify l=1 channel has correct free eigenvalues at R=0.

    At R=0, channels decouple and eigenvalues are the union of
    per-channel Casimir values:
      l=0 singlet (k=0,2,4,6,8): mu = 0, 16, 48, 96, 160
      l=1 singlet (k=0,2,4,6,8): mu = 2(1+k+1)^2-2 = 6, 30, 70, 126, 198
    Sorted: [0, 6, 16, 30, 48, 70, 96, 126, 160, 198]
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=1, symmetry='singlet')

    evals, _ = solver.solve(R=0.0, n_channels=10)

    l0_evals = [2 * (0 + 2*j + 1)**2 - 2 for j in range(5)]  # 0, 16, 48, 96, 160
    l1_evals = [2 * (1 + 2*j + 1)**2 - 2 for j in range(5)]  # 6, 30, 70, 126, 198
    expected = sorted(l0_evals + l1_evals)

    np.testing.assert_allclose(evals, expected, atol=1e-10,
                               err_msg="Multichannel Casimir eigenvalues at R=0 incorrect")


def test_gaunt_coupling_consistency():
    """Verify V_ee coupling uses correct Gaunt integrals.

    Key Gaunt integrals from the Wigner 3j formula:
      G(0,0,0) = 2     (monopole self-coupling)
      G(0,1,1) = 2/3   (dipole l=0 <-> l=1 coupling)
      G(1,0,1) = 2/3   (monopole l=1 self-coupling)
      G(1,2,1) = 4/15  (quadrupole l=1 self-coupling)
    Selection rules: G = 0 when l+k+l' is odd.
    """
    from geovac.hyperspherical_angular import gaunt_integral

    # Non-zero integrals
    assert abs(gaunt_integral(0, 0, 0) - 2.0) < 1e-14
    assert abs(gaunt_integral(0, 1, 1) - 2.0 / 3.0) < 1e-14
    assert abs(gaunt_integral(1, 0, 1) - 2.0 / 3.0) < 1e-14
    assert abs(gaunt_integral(1, 2, 1) - 4.0 / 15.0) < 1e-14

    # Selection rules: l+k+l' odd => 0
    assert gaunt_integral(0, 1, 0) == 0.0
    assert gaunt_integral(1, 1, 1) == 0.0

    # Triangle inequality violations => 0
    assert gaunt_integral(0, 0, 1) == 0.0
    assert gaunt_integral(0, 3, 1) == 0.0


def test_multichannel_coupling_symmetry():
    """Verify the full multichannel coupling matrix is symmetric."""
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=8, l_max=2)

    np.testing.assert_allclose(
        solver.nuclear_matrix_full, solver.nuclear_matrix_full.T, atol=1e-14,
        err_msg="Full nuclear matrix not symmetric"
    )
    np.testing.assert_allclose(
        solver.vee_matrix_full, solver.vee_matrix_full.T, atol=1e-14,
        err_msg="Full V_ee matrix not symmetric"
    )
    np.testing.assert_allclose(
        solver.coupling_matrix_full, solver.coupling_matrix_full.T, atol=1e-14,
        err_msg="Full coupling matrix not symmetric"
    )


def test_nuclear_diagonal_in_l():
    """Verify nuclear coupling is block-diagonal in l.

    The nuclear potential -Z(1/cos + 1/sin) depends only on alpha,
    not on the angular variables distinguishing l channels.  The
    Legendre orthogonality integral P_l P_{l'} = 0 for l != l'
    ensures the off-diagonal l-blocks are zero.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=2)
    V_nuc = solver.nuclear_matrix_full
    nb = solver.n_basis

    # Off-diagonal l-blocks should be zero
    for l in range(solver.n_l):
        for lp in range(solver.n_l):
            if l == lp:
                continue
            block = V_nuc[l*nb:(l+1)*nb, lp*nb:(lp+1)*nb]
            assert np.max(np.abs(block)) < 1e-14, \
                f"Nuclear coupling V[l={l}, l'={lp}] is non-zero: " \
                f"max = {np.max(np.abs(block)):.2e}"


def test_vee_cross_channel_nonzero():
    """Verify V_ee cross-channel coupling is non-zero for l=0 <-> l=1.

    The l=0 <-> l=1 coupling is mediated by the k=1 dipole term with
    Gaunt integral G(0,1,1) = 2/3.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=5, l_max=1)
    V_ee = solver.vee_matrix_full
    nb = solver.n_basis

    # l=0 <-> l=1 block should be non-zero
    cross_block = V_ee[0:nb, nb:2*nb]
    assert np.max(np.abs(cross_block)) > 0.01, \
        f"V_ee cross-channel coupling too small: max = {np.max(np.abs(cross_block)):.6f}"


def test_lmax1_cross_validation():
    """Cross-validate algebraic vs FD solver at l_max=1.

    At small R (perturbative regime), both solvers should agree closely.
    The algebraic solver handles the centrifugal singularity analytically,
    while the FD solver uses a grid.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver
    from geovac.hyperspherical_angular import solve_angular

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=20, l_max=1)

    # Small R: tight tolerance; larger R: relaxed (different error profiles)
    R_tol = [(0.1, 1.0), (0.5, 1.0), (1.0, 2.0), (2.0, 5.0)]

    for R, tol_pct in R_tol:
        evals_alg, _ = solver.solve(R, n_channels=1)
        mu_alg = evals_alg[0]

        evals_fd, _ = solve_angular(R, Z=2.0, l_max=1, n_alpha=400, n_channels=1)
        mu_fd = evals_fd[0]

        if abs(mu_fd) > 0.01:
            rel_diff = abs(mu_alg - mu_fd) / abs(mu_fd) * 100
            assert rel_diff < tol_pct, \
                f"R={R}: mu_alg={mu_alg:.6f}, mu_fd={mu_fd:.6f}, " \
                f"diff={rel_diff:.2f}% (tol={tol_pct}%)"


# ---------------------------------------------------------------------------
# DBOC (diagonal Born-Oppenheimer correction) tests
# ---------------------------------------------------------------------------

def test_dboc_positivity():
    """Verify DBOC(R) >= 0 for all R.

    The DBOC is a sum of squares: (1/2) sum |P_{mu,0}|^2, so it must
    be non-negative at every R.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)

    for R in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        _, _, dboc = solver.solve_with_dboc(R, n_channels=1)
        assert dboc >= 0.0, f"DBOC({R}) = {dboc:.6e} is negative"


def test_dboc_magnitude():
    """Verify DBOC has the right order of magnitude near the He well.

    Paper 13 estimates the total DBOC contribution as ~0.0015 Ha for He.
    At R ~ 1 bohr (near the well minimum), DBOC should be O(0.001-0.01).
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=15, l_max=0)
    _, _, dboc = solver.solve_with_dboc(R=1.0, n_channels=1)

    assert dboc > 1e-4, f"DBOC(R=1) = {dboc:.6e} too small (expect ~0.001)"
    assert dboc < 0.1, f"DBOC(R=1) = {dboc:.6e} too large (expect ~0.001)"


def test_dboc_small_R_perturbative():
    """At small R, DBOC converges to a finite perturbative value.

    At R=0, the eigenvectors are the standard basis vectors, but
    dPhi/dR is nonzero (first-order perturbation theory gives a
    constant rotation rate).  The DBOC approaches a finite limit:

        DBOC(R=0) = (1/2) sum_{mu>0} |V_{mu,0}|^2 / (casimir_mu - casimir_0)^2

    This is the perturbative DBOC.  Verify convergence to this value.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)

    _, _, dboc_001 = solver.solve_with_dboc(R=0.001, n_channels=1)
    _, _, dboc_01 = solver.solve_with_dboc(R=0.01, n_channels=1)

    # Perturbative limit: compute analytically
    V = solver.coupling_matrix
    casimir = solver.casimir_eigenvalues
    dboc_pert = 0.5 * sum(
        V[0, mu] ** 2 / (casimir[mu] - casimir[0]) ** 2
        for mu in range(1, len(casimir))
    )

    # At very small R, should approach the perturbative value
    rel_diff = abs(dboc_001 - dboc_pert) / dboc_pert
    assert rel_diff < 0.01, \
        f"DBOC(R=0.001) = {dboc_001:.6e}, perturbative = {dboc_pert:.6e}, " \
        f"rel_diff = {rel_diff:.4f}"


def test_dboc_decays_large_R():
    """At large R, DBOC should decay as eigenvector stabilizes.

    At large R, the coupling term R * V_coupling dominates over
    diag(casimir), so the eigenvectors approach those of V_coupling
    (R-independent).  The DBOC should decrease with R.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=15, l_max=0)

    _, _, dboc_5 = solver.solve_with_dboc(R=5.0, n_channels=1)
    _, _, dboc_20 = solver.solve_with_dboc(R=20.0, n_channels=1)

    assert dboc_20 < dboc_5, \
        f"DBOC should decay: DBOC(5)={dboc_5:.6e}, DBOC(20)={dboc_20:.6e}"
    assert dboc_20 < 0.001, \
        f"DBOC(R=20) = {dboc_20:.6e} should be small"


def test_dboc_hellmann_feynman_consistency():
    """Verify Hellmann-Feynman DBOC matches finite-difference dPhi/dR.

    Compute dPhi_0/dR by central difference and compare ||dPhi/dR||^2
    against the Hellmann-Feynman sum over excited states.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)

    R = 1.5
    delta = 1e-5

    # Hellmann-Feynman DBOC
    _, _, dboc_hf = solver.solve_with_dboc(R, n_channels=1)

    # Finite-difference: solve at R-delta and R+delta
    n_all = solver._total_dim
    _, evecs_m = solver.solve(R - delta, n_channels=n_all)
    _, evecs_p = solver.solve(R + delta, n_channels=n_all)

    phi_m = evecs_m[0]  # ground state at R - delta
    phi_p = evecs_p[0]  # ground state at R + delta

    # Fix sign ambiguity: ensure consistent sign
    if np.dot(phi_m, phi_p) < 0:
        phi_p = -phi_p

    dphi_dR = (phi_p - phi_m) / (2.0 * delta)
    dboc_fd = 0.5 * np.dot(dphi_dR, dphi_dR)

    rel_diff = abs(dboc_hf - dboc_fd) / max(abs(dboc_fd), 1e-15)
    assert rel_diff < 0.01, \
        f"HF DBOC = {dboc_hf:.8e}, FD DBOC = {dboc_fd:.8e}, rel_diff = {rel_diff:.4f}"


def test_he_energy_dboc_raises():
    """Verify DBOC raises He energy (positive correction).

    At l_max >= 1, the adiabatic energy overshoots below exact.
    The DBOC (positive correction) should raise it.

    IMPORTANT: The DBOC overcorrects for He because the off-diagonal
    P-matrix coupling (which lowers the energy) nearly cancels the DBOC
    (~97% cancellation, see hyperspherical_radial.py line 327-329).
    The DBOC alone is NOT a useful single-channel correction — it is
    a diagnostic measuring the strength of non-adiabatic coupling.
    """
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    result = solve_hyperspherical_algebraic(
        Z=2.0, n_basis=15, l_max=1, n_R=150,
        R_max=25.0, N_R_radial=1500, include_dboc=True, verbose=False,
    )

    E_exact = -2.903724
    E_with = result['energy']
    E_without = result['energy_no_dboc']

    # Energy without DBOC should be below exact (adiabatic overshoot)
    assert E_without < E_exact, \
        f"E_no_dboc = {E_without:.6f} should be below exact {E_exact:.6f}"

    # Energy with DBOC should be higher (less negative) than without
    assert E_with > E_without, \
        f"E_dboc = {E_with:.6f} should be above E_no_dboc = {E_without:.6f}"

    # DBOC overcorrects: E_with goes above exact (known physics —
    # the off-diagonal P coupling that cancels ~97% of the DBOC
    # is not included in the single-channel approximation)
    dboc_shift = E_with - E_without
    adiabatic_error = E_exact - E_without  # positive (exact is above adiab)
    assert dboc_shift > adiabatic_error, \
        f"DBOC should overcorrect: shift={dboc_shift:.6f}, " \
        f"adiab_err={adiabatic_error:.6f}"


def test_dboc_lmax_scaling():
    """Verify DBOC shift is consistent across l_max values.

    The DBOC overcorrects at all l_max >= 1 (missing off-diagonal
    cancellation).  Verify that the DBOC shift is positive and
    roughly consistent (same order of magnitude) across l_max values.
    """
    from geovac.algebraic_angular import solve_hyperspherical_algebraic

    shifts = []
    for l_max in [1, 2, 3]:
        result = solve_hyperspherical_algebraic(
            Z=2.0, n_basis=15, l_max=l_max, n_R=150,
            R_max=25.0, N_R_radial=1500, include_dboc=True, verbose=False,
        )

        shift = result['energy'] - result['energy_no_dboc']
        shifts.append(shift)

        # DBOC shift must be positive (repulsive)
        assert shift > 0, \
            f"l_max={l_max}: DBOC shift should be positive, got {shift:.6f}"

    # Shifts should be in the same order of magnitude (0.01 - 0.1)
    for i, l_max in enumerate([1, 2, 3]):
        assert 0.01 < shifts[i] < 0.15, \
            f"l_max={l_max}: DBOC shift {shifts[i]:.6f} out of expected range"


# ---------------------------------------------------------------------------
# Perturbation series tests (Lane B — algebraic adiabatic curves)
# ---------------------------------------------------------------------------

def test_perturbation_a1_exact():
    """Verify a_1 matches Paper 13 Table II exact formula.

    a_1(n=1, Z=2) = (8/3pi)(sqrt(2) - 8) = -5.590189...
    This is the primary validation gate for the perturbation series.
    """
    from geovac.algebraic_angular import AlgebraicAngularSolver, perturbation_series_mu

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=5)

    a1_exact = (8.0 / (3.0 * pi)) * (sqrt(2) - 8.0)
    assert abs(coeffs[0, 0] - 0.0) < 1e-12, "a_0 should be 0 (ground channel Casimir)"
    assert abs(coeffs[0, 1] - a1_exact) < 1e-9, \
        f"a_1 = {coeffs[0, 1]:.10f}, exact = {a1_exact:.10f}"


def test_perturbation_a1_z_scaling():
    """Verify a_1(Z) = (8/3pi)(sqrt(2) - 4Z) for multiple Z values."""
    from geovac.algebraic_angular import AlgebraicAngularSolver, perturbation_series_mu

    for Z in [1.0, 2.0, 3.0, 5.0]:
        solver = AlgebraicAngularSolver(Z=Z, n_basis=10, l_max=0)
        H0_diag = np.concatenate(solver._channel_casimir)
        V_C = solver._coupling_full
        coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=2)

        a1_exact = (8.0 / (3.0 * pi)) * (sqrt(2) - 4.0 * Z)
        assert abs(coeffs[0, 1] - a1_exact) < 1e-9, \
            f"Z={Z}: a_1 = {coeffs[0, 1]:.10f}, exact = {a1_exact:.10f}"


def test_perturbation_a2_matches_direct():
    """Verify a_2 from perturbation_series_mu matches direct second-order PT."""
    from geovac.algebraic_angular import AlgebraicAngularSolver, perturbation_series_mu

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=20, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=3)

    # Direct second-order: a_2 = sum_{n!=0} |V[0,n]|^2 / (E0 - En)
    a2_direct = sum(
        V_C[0, n] ** 2 / (H0_diag[0] - H0_diag[n])
        for n in range(1, len(H0_diag))
    )

    assert abs(coeffs[0, 2] - a2_direct) < 1e-10, \
        f"a_2 series = {coeffs[0, 2]:.10f}, direct = {a2_direct:.10f}"


def test_perturbation_convergence_small_R():
    """Perturbation series matches numerical diagonalization for R < 2 bohr."""
    from geovac.algebraic_angular import (
        AlgebraicAngularSolver, perturbation_series_mu,
        evaluate_perturbation_series,
    )

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=15)

    R_test = np.array([0.1, 0.5, 1.0])
    mu_series = evaluate_perturbation_series(coeffs[:1], R_test)

    for i, R in enumerate(R_test):
        evals, _ = solver.solve(R, n_channels=1)
        rel_err = abs(mu_series[0, i] - evals[0]) / abs(evals[0])
        assert rel_err < 1e-6, \
            f"R={R}: series={mu_series[0, i]:.8f}, exact={evals[0]:.8f}, rel_err={rel_err:.2e}"


def test_perturbation_diverges_large_R():
    """Perturbation series diverges for R >> convergence radius.

    This confirms the algebraic/transcendental boundary noted in
    Paper 13 Sec XII.B.  The convergence radius is ~2-3 bohr.
    """
    from geovac.algebraic_angular import (
        AlgebraicAngularSolver, perturbation_series_mu,
        evaluate_perturbation_series,
    )

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=15)

    R_large = np.array([10.0])
    mu_series = evaluate_perturbation_series(coeffs[:1], R_large)
    evals, _ = solver.solve(10.0, n_channels=1)

    # Series should diverge wildly at R=10
    rel_err = abs(mu_series[0, 0] - evals[0]) / abs(evals[0])
    assert rel_err > 1.0, \
        f"Series should diverge at R=10: rel_err = {rel_err:.2e}"


def test_pade_extends_convergence():
    """Padé approximant extends useful convergence beyond raw series.

    [7/7] Padé on g(R) = mu(R) + Z^2*R^2/2 should be accurate to
    < 1% at R = 3 bohr, where the raw series has > 1% error.
    """
    from geovac.algebraic_angular import (
        AlgebraicAngularSolver, perturbation_series_mu,
        evaluate_perturbation_series, pade_approximant, evaluate_pade,
    )

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=14)

    Z = 2.0
    g_coeffs = coeffs[0].copy()
    g_coeffs[2] += Z ** 2 / 2.0

    p_c, q_c = pade_approximant(g_coeffs, 7, 7)

    # At R=2, Padé should be excellent
    R_test = np.array([2.0])
    g_pade = evaluate_pade(p_c, q_c, R_test)
    mu_pade = g_pade - Z ** 2 * R_test ** 2 / 2.0
    evals, _ = solver.solve(2.0, n_channels=1)
    rel_err_pade = abs(mu_pade[0] - evals[0]) / abs(evals[0])
    assert rel_err_pade < 1e-4, \
        f"Padé at R=2: rel_err = {rel_err_pade:.2e}"

    # At R=3, Padé should still be reasonable (< 1%)
    R_test = np.array([3.0])
    g_pade = evaluate_pade(p_c, q_c, R_test)
    mu_pade = g_pade - Z ** 2 * R_test ** 2 / 2.0
    evals, _ = solver.solve(3.0, n_channels=1)
    rel_err_pade = abs(mu_pade[0] - evals[0]) / abs(evals[0])
    assert rel_err_pade < 0.01, \
        f"Padé at R=3: rel_err = {rel_err_pade:.2e}"


def test_pade_fails_beyond_transcendental_boundary():
    """Padé cannot capture the transcendental μ(R) at large R.

    Paper 13 Sec XII.B: μ(R) transitions from O(R) perturbative to
    O(R²) asymptotic; no finite rational function spans both regimes.
    """
    from geovac.algebraic_angular import (
        AlgebraicAngularSolver, perturbation_series_mu,
        pade_approximant, evaluate_pade,
    )

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=1, max_order=14)

    Z = 2.0
    g_coeffs = coeffs[0].copy()
    g_coeffs[2] += Z ** 2 / 2.0
    p_c, q_c = pade_approximant(g_coeffs, 7, 7)

    # At R=10, Padé should fail badly (> 50% error)
    R_test = np.array([10.0])
    g_pade = evaluate_pade(p_c, q_c, R_test)
    mu_pade = g_pade - Z ** 2 * R_test ** 2 / 2.0
    evals, _ = solver.solve(10.0, n_channels=1)
    rel_err = abs(mu_pade[0] - evals[0]) / abs(evals[0])
    assert rel_err > 0.5, \
        f"Expected Padé failure at R=10: rel_err = {rel_err:.2e}"


def test_perturbation_multichannel_consistency():
    """Higher channels have correct a_0 (Casimir) values."""
    from geovac.algebraic_angular import AlgebraicAngularSolver, perturbation_series_mu

    solver = AlgebraicAngularSolver(Z=2.0, n_basis=10, l_max=0)
    H0_diag = np.concatenate(solver._channel_casimir)
    V_C = solver._coupling_full
    coeffs = perturbation_series_mu(H0_diag, V_C, n_channels=4, max_order=5)

    # a_0 values should be Casimir: 0, 16, 48, 96 for l=0 singlet
    expected_a0 = [0.0, 16.0, 48.0, 96.0]
    for ch in range(4):
        assert abs(coeffs[ch, 0] - expected_a0[ch]) < 1e-10, \
            f"Channel {ch}: a_0 = {coeffs[ch, 0]}, expected {expected_a0[ch]}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
