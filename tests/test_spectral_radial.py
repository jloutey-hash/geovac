"""
Tests for the spectral Laguerre radial solver in prolate spheroidal coordinates.

Validates that the spectral solver (radial_method='spectral') reproduces
H2+ energies with dramatically fewer degrees of freedom than the FD solver
(20 basis functions vs 5000 grid points), at equal or better accuracy.

Exact references (Bates, Ledsham & Stewart 1953):
  R_eq = 1.997 bohr
  E_total(R_eq) = -0.6026342 Ha
  D_e = 0.1026 Ha
"""
import numpy as np
import pytest
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    scan_h2plus_pes,
    fit_spectroscopic_constants,
)


class TestSpectralBasic:
    """Basic functionality of the spectral radial solver."""

    def test_solve_r2_bound(self):
        """Spectral solver at R=2.0 returns bound energy."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_el, c2, A = lat.solve()
        E_tot = E_el + 0.5
        assert E_tot < -0.5, f"Should be bound, got E_tot={E_tot:.6f}"

    def test_c2_positive(self):
        """c^2 should be positive for bound state."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        _, c2, _ = lat.solve()
        assert c2 > 0, f"c^2={c2} should be positive"

    def test_separation_constant_positive(self):
        """Separation constant A should be positive for sigma_g ground state."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        _, _, A = lat.solve()
        assert A > 0, f"A={A} should be positive"

    def test_total_energy_includes_vnn(self):
        """total_energy() = E_elec + Z_A*Z_B/R."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_tot = lat.total_energy()
        E_el, _, _ = lat.solve()
        assert abs(E_tot - (E_el + 0.5)) < 1e-10


class TestSpectralAccuracy:
    """Accuracy benchmarks for the spectral solver."""

    def test_accuracy_vs_exact_r2(self):
        """n_basis=20 should give < 0.01% error at R=2.0."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        E_el, _, _ = lat.solve()
        E_tot = E_el + 0.5
        # Exact: -0.6026342 Ha
        err = abs(E_tot - (-0.6026342)) / 0.6026342 * 100
        assert err < 0.01, f"Error {err:.4f}% > 0.01% at n_basis=20"

    def test_accuracy_vs_exact_r1(self):
        """Check accuracy at R=1.0 (short range)."""
        lat = ProlateSpheroidalLattice(R=1.0, radial_method='spectral', n_basis=20)
        E_el, _, _ = lat.solve()
        # E_elec at R=1.0 should be well below -1.0 (united atom He+ limit -2.0)
        assert E_el < -1.0, f"E_elec={E_el:.4f} should be < -1.0 at R=1.0"

    def test_accuracy_vs_fd(self):
        """Spectral n_basis=20 should match or beat FD N_xi=5000."""
        lat_sp = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
        lat_fd = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        E_sp, _, _ = lat_sp.solve()
        E_fd, _, _ = lat_fd.solve()
        # Spectral should be lower (more accurate, variational)
        # or at least within 0.01 Ha
        assert E_sp <= E_fd + 0.01, \
            f"Spectral E={E_sp:.6f} much higher than FD E={E_fd:.6f}"

    def test_convergence_with_basis_size(self):
        """Energy should converge monotonically with increasing n_basis."""
        energies = []
        for nb in [5, 10, 15, 20]:
            lat = ProlateSpheroidalLattice(
                R=2.0, radial_method='spectral', n_basis=nb
            )
            E_el, _, _ = lat.solve()
            energies.append(E_el + 0.5)
        # Should be converging (becoming more negative or stable)
        for i in range(1, len(energies)):
            # Allow small tolerance for non-monotonicity due to alpha adaptation
            assert energies[i] <= energies[i - 1] + 1e-4, \
                f"E[nb={[5,10,15,20][i]}]={energies[i]:.8f} > " \
                f"E[nb={[5,10,15,20][i-1]}]={energies[i-1]:.8f}"

    def test_small_basis_still_reasonable(self):
        """Even n_basis=5 should give < 5% error."""
        lat = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=5)
        E_el, _, _ = lat.solve()
        E_tot = E_el + 0.5
        err = abs(E_tot - (-0.6026)) / 0.6026 * 100
        assert err < 5.0, f"Error {err:.2f}% > 5% even with n_basis=5"


class TestSpectralPES:
    """PES scan tests with the spectral solver."""

    def test_pes_has_minimum(self):
        """Spectral PES should have a minimum between R=1 and R=5."""
        R_values = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
        E_tot = []
        for R in R_values:
            lat = ProlateSpheroidalLattice(
                R=R, radial_method='spectral', n_basis=20
            )
            E_tot.append(lat.total_energy())
        E_tot = np.array(E_tot)
        idx_min = np.argmin(E_tot)
        assert 0 < idx_min < len(R_values) - 1, \
            f"Minimum at boundary (idx={idx_min})"

    def test_spectroscopic_constants(self):
        """Spectral PES should give accurate R_eq and D_e."""
        # Fine grid around equilibrium for accurate polynomial fit
        R_values = np.arange(1.5, 3.01, 0.05)
        E_tot = []
        for R in R_values:
            lat = ProlateSpheroidalLattice(
                R=R, radial_method='spectral', n_basis=20
            )
            E_tot.append(lat.total_energy())
        E_tot = np.array(E_tot)
        fit = fit_spectroscopic_constants(R_values, E_tot)

        assert not fit['boundary'], "Minimum should not be at boundary"
        # R_eq should be within 1% of exact 1.997
        assert abs(fit['R_eq'] - 1.997) / 1.997 < 0.01, \
            f"R_eq={fit['R_eq']:.4f}, expected ~1.997 (> 1% err)"
        # E_min should be within 0.01% of exact
        assert abs(fit['E_min'] - (-0.6026342)) / 0.6026342 < 0.0001, \
            f"E_min={fit['E_min']:.8f}, expected ~-0.6026342 (> 0.01% err)"
        # D_e > 0
        assert fit['D_e'] > 0, f"D_e={fit['D_e']:.6f} should be positive"


class TestSpectralDissociation:
    """Dissociation limit tests."""

    def test_large_r_approaches_h_atom(self):
        """At large R, E_total should approach -0.5 Ha."""
        lat = ProlateSpheroidalLattice(
            R=10.0, radial_method='spectral', n_basis=25
        )
        E_tot = lat.total_energy()
        assert abs(E_tot - (-0.5)) < 0.01, \
            f"E(R=10)={E_tot:.6f}, expected ~-0.5"


class TestSpectralHeteronuclear:
    """Heteronuclear tests (HeH^2+)."""

    def test_heh2plus_bound(self):
        """HeH^2+ should be well bound."""
        lat = ProlateSpheroidalLattice(
            R=2.0, Z_A=2, Z_B=1,
            radial_method='spectral', n_basis=20
        )
        E_tot = lat.total_energy()
        assert E_tot < -1.0, f"HeH2+ E={E_tot:.4f} should be well bound"


class TestSpectralPerformance:
    """Performance benchmarks."""

    def test_speedup_over_fd(self):
        """Spectral solver should be at least 10x faster than FD N_xi=5000."""
        # Spectral
        t0 = time.time()
        lat_sp = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20
        )
        lat_sp.solve()
        dt_sp = time.time() - t0

        # FD
        t0 = time.time()
        lat_fd = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        lat_fd.solve()
        dt_fd = time.time() - t0

        speedup = dt_fd / max(dt_sp, 1e-6)
        assert speedup > 10, \
            f"Speedup {speedup:.1f}x < 10x (spectral={dt_sp:.3f}s, fd={dt_fd:.3f}s)"

    def test_dimension_reduction(self):
        """Spectral should use <=25 basis functions vs 5000 FD grid points."""
        lat = ProlateSpheroidalLattice(
            R=2.0, radial_method='spectral', n_basis=20
        )
        assert lat.n_basis <= 25
        assert lat.n_basis < lat.N_xi / 100  # >100x reduction


class TestFDUnchanged:
    """Verify that adding the spectral option didn't break FD."""

    def test_fd_default_unchanged(self):
        """Default radial_method should be 'fd'."""
        lat = ProlateSpheroidalLattice(R=2.0)
        assert lat.radial_method == 'fd'

    def test_fd_results_unchanged(self):
        """FD results should be identical to before."""
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=5000, radial_method='fd')
        E_el, c2, A = lat.solve()
        E_tot = E_el + 0.5
        # Should match the known FD result
        err = abs(E_tot - (-0.6026)) / 0.6026 * 100
        assert err < 1.5, f"FD error {err:.2f}% > 1.5% (regression)"


from geovac.prolate_spheroidal_lattice import (
    _laguerre_moment_matrices,
    _build_laguerre_matrices_algebraic,
)


class TestLaguerreMomentMatrices:
    """Validate the algebraic Laguerre moment matrices M0, M1, M2."""

    def test_m0_is_identity(self):
        """M0 = delta_{mn} (Laguerre orthonormality)."""
        M0, _, _ = _laguerre_moment_matrices(10)
        np.testing.assert_allclose(M0, np.eye(10), atol=1e-15)

    def test_m1_tridiagonal_symmetry(self):
        """M1 should be symmetric and tridiagonal."""
        _, M1, _ = _laguerre_moment_matrices(10)
        np.testing.assert_allclose(M1, M1.T, atol=1e-15)
        # Check no entries beyond tridiagonal band
        for i in range(10):
            for j in range(10):
                if abs(i - j) > 1:
                    assert M1[i, j] == 0.0

    def test_m1_known_values(self):
        """Spot-check M1 entries against the recurrence formula."""
        _, M1, _ = _laguerre_moment_matrices(5)
        # M1[i,i] = 2i+1
        for i in range(5):
            assert M1[i, i] == 2 * i + 1
        # M1[i,i+1] = -(i+1)
        for i in range(4):
            assert M1[i, i + 1] == -(i + 1)

    def test_m2_pentadiagonal_symmetry(self):
        """M2 should be symmetric and pentadiagonal."""
        _, _, M2 = _laguerre_moment_matrices(10)
        np.testing.assert_allclose(M2, M2.T, atol=1e-15)
        for i in range(10):
            for j in range(10):
                if abs(i - j) > 2:
                    assert M2[i, j] == 0.0

    def test_m2_diagonal_values(self):
        """M2[j,j] = 6j^2 + 6j + 2."""
        _, _, M2 = _laguerre_moment_matrices(8)
        for j in range(8):
            expected = 6 * j * j + 6 * j + 2
            assert M2[j, j] == expected, f"M2[{j},{j}]={M2[j,j]}, expected {expected}"

    def test_moments_vs_quadrature(self):
        """Moment matrices should match Gauss-Laguerre quadrature reference."""
        from scipy.special import roots_laguerre, eval_laguerre
        N = 8
        M0, M1, M2 = _laguerre_moment_matrices(N)
        n_quad = 200
        x, w = roots_laguerre(n_quad)
        L = np.array([eval_laguerre(n, x) for n in range(N)])
        # Build reference moment matrices via quadrature
        M0_ref = (w[np.newaxis, :] * L) @ L.T
        M1_ref = (w[np.newaxis, :] * L * x[np.newaxis, :]) @ L.T
        M2_ref = (w[np.newaxis, :] * L * (x**2)[np.newaxis, :]) @ L.T
        np.testing.assert_allclose(M0, M0_ref, atol=1e-10)
        np.testing.assert_allclose(M1, M1_ref, atol=1e-10)
        np.testing.assert_allclose(M2, M2_ref, atol=1e-9)


class TestAlgebraicVsQuadratureMatrices:
    """Compare algebraic and quadrature matrix construction element-by-element."""

    def _get_both_matrices(self, R: float = 2.0, n_basis: int = 20):
        """Build H,S matrices both ways at converged c2, A."""
        sq = ProlateSpheroidalLattice(
            R=R, radial_method='spectral', n_basis=n_basis,
            matrix_method='quadrature'
        )
        E, c2, A = sq.solve()
        alpha = max(np.sqrt(max(c2, 0.01)), 0.5)
        H_q, S_q = sq._build_laguerre_matrices_quadrature(n_basis, alpha, A, c2)
        a_param = R * 2  # Z_A + Z_B = 2 for H2+
        H_a, S_a = _build_laguerre_matrices_algebraic(n_basis, alpha, A, a_param, c2, m=0)
        return H_q, S_q, H_a, S_a

    def test_overlap_agreement_r2(self):
        """Overlap matrices match to machine precision at R=2."""
        _, S_q, _, S_a = self._get_both_matrices(R=2.0)
        np.testing.assert_allclose(S_a, S_q, atol=1e-14)

    def test_hamiltonian_agreement_r2(self):
        """Hamiltonian matrices match to < 1e-11 at R=2."""
        H_q, _, H_a, _ = self._get_both_matrices(R=2.0)
        np.testing.assert_allclose(H_a, H_q, atol=1e-11)

    def test_overlap_agreement_r05(self):
        """Overlap matrices match at R=0.5 (strong binding).

        The algebraic overlap is exactly diagonal; quadrature introduces
        ~1e-14 off-diagonal noise from finite precision.
        """
        _, S_q, _, S_a = self._get_both_matrices(R=0.5)
        np.testing.assert_allclose(S_a, S_q, atol=2e-14)

    def test_hamiltonian_agreement_r05(self):
        """Hamiltonian matrices match at R=0.5."""
        H_q, _, H_a, _ = self._get_both_matrices(R=0.5)
        np.testing.assert_allclose(H_a, H_q, atol=1e-11)

    def test_hamiltonian_agreement_r8(self):
        """Hamiltonian matrices match at R=8.0 (near dissociation)."""
        H_q, _, H_a, _ = self._get_both_matrices(R=8.0)
        np.testing.assert_allclose(H_a, H_q, atol=1e-11)


class TestAlgebraicSolverEnergies:
    """Verify algebraic solver produces identical energies to quadrature."""

    def test_energy_r2(self):
        """Algebraic and quadrature give same energy at R=2."""
        sq = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', matrix_method='quadrature')
        sa = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', matrix_method='algebraic')
        Eq, _, _ = sq.solve()
        Ea, _, _ = sa.solve()
        assert abs(Ea - Eq) < 1e-13, f"Energy diff {abs(Ea-Eq):.2e} > 1e-13"

    def test_energy_multiple_r(self):
        """Algebraic matches quadrature across full PES range."""
        for R in [0.5, 1.0, 2.0, 4.0, 8.0, 15.0]:
            sq = ProlateSpheroidalLattice(R=R, radial_method='spectral', matrix_method='quadrature')
            sa = ProlateSpheroidalLattice(R=R, radial_method='spectral', matrix_method='algebraic')
            Eq, _, _ = sq.solve()
            Ea, _, _ = sa.solve()
            assert abs(Ea - Eq) < 1e-12, f"R={R}: diff {abs(Ea-Eq):.2e}"

    def test_pes_scan_algebraic(self):
        """Full PES scan with algebraic method matches quadrature."""
        R_vals = np.linspace(1.0, 6.0, 10)
        res_q = scan_h2plus_pes(R_vals, radial_method='spectral', matrix_method='quadrature', verbose=False)
        res_a = scan_h2plus_pes(R_vals, radial_method='spectral', matrix_method='algebraic', verbose=False)
        np.testing.assert_allclose(res_a['E_total'], res_q['E_total'], atol=1e-13)

    def test_algebraic_accuracy_vs_exact(self):
        """Algebraic solver achieves < 0.01% error vs exact H2+ at R=2."""
        sa = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', matrix_method='algebraic')
        E_tot = sa.total_energy()
        exact = -0.6026342
        err_pct = abs(E_tot - exact) / abs(exact) * 100
        assert err_pct < 0.01, f"Error {err_pct:.4f}% > 0.01%"

    def test_heteronuclear_algebraic(self):
        """Algebraic solver works for HeH2+ (heteronuclear)."""
        sq = ProlateSpheroidalLattice(R=2.0, Z_A=2, Z_B=1, radial_method='spectral', matrix_method='quadrature')
        sa = ProlateSpheroidalLattice(R=2.0, Z_A=2, Z_B=1, radial_method='spectral', matrix_method='algebraic')
        Eq, _, _ = sq.solve()
        Ea, _, _ = sa.solve()
        assert abs(Ea - Eq) < 1e-12, f"HeH2+ diff {abs(Ea-Eq):.2e}"

    def test_m_nonzero_raises(self):
        """Algebraic method should raise NotImplementedError for m!=0."""
        with pytest.raises(NotImplementedError, match="m=0"):
            _build_laguerre_matrices_algebraic(20, 1.0, 0.8, 4.0, 2.0, m=1)

    def test_m_nonzero_falls_back_to_quadrature(self):
        """For m!=0, spectral solver with algebraic falls back to quadrature."""
        sq = ProlateSpheroidalLattice(
            R=2.0, m=1, n_angular=0, n_radial=0,
            radial_method='spectral', matrix_method='quadrature'
        )
        sa = ProlateSpheroidalLattice(
            R=2.0, m=1, n_angular=0, n_radial=0,
            radial_method='spectral', matrix_method='algebraic'
        )
        Eq, _, _ = sq.solve()
        Ea, _, _ = sa.solve()
        assert abs(Ea - Eq) < 1e-13, f"m=1 fallback diff {abs(Ea-Eq):.2e}"

    def test_quadrature_default_unchanged(self):
        """Default matrix_method='quadrature' should not change results."""
        s_default = ProlateSpheroidalLattice(R=2.0, radial_method='spectral')
        s_explicit = ProlateSpheroidalLattice(R=2.0, radial_method='spectral', matrix_method='quadrature')
        Ed, _, _ = s_default.solve()
        Ee, _, _ = s_explicit.solve()
        assert Ed == Ee, "Default should be quadrature"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
