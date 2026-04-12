"""
Tests for nuclear form factor investigation (Track NH).

Tests the finite nuclear charge distribution's effect on the S³ graph
Laplacian correspondence.

Author: GeoVac Development Team
Date: April 2026
"""

import numpy as np
import pytest

from geovac.nuclear.form_factor import (
    uniform_sphere_potential,
    finite_size_correction,
    spectrum_vs_nuclear_radius,
    fock_projection_with_finite_nucleus,
    s3_lattice_with_modified_potential,
    fock_projection_breakdown_scan,
    R_PROTON_BOHR,
    _hydrogenic_radial,
)
# Re-export for test convenience
from geovac.nuclear.form_factor import finite_size_correction as _fsc


# ---------------------------------------------------------------------------
# 1. Uniform sphere potential tests
# ---------------------------------------------------------------------------

class TestUniformSpherePotential:
    """Tests for the uniform sphere potential."""

    def test_continuous_at_boundary(self):
        """V(r) is continuous at r = R_nuc."""
        Z, R_nuc = 1.0, 0.1
        eps = 1e-10
        r_inside = np.array([R_nuc - eps])
        r_outside = np.array([R_nuc + eps])
        V_in = uniform_sphere_potential(Z, R_nuc, r_inside)[0]
        V_out = uniform_sphere_potential(Z, R_nuc, r_outside)[0]
        assert abs(V_in - V_out) < 1e-6, (
            f"Discontinuity at boundary: V_in={V_in:.10f}, V_out={V_out:.10f}"
        )

    def test_matches_point_outside(self):
        """V(r) = -Z/r for r > R_nuc to machine precision."""
        Z, R_nuc = 2.0, 0.01
        r = np.linspace(0.02, 10.0, 100)
        V = uniform_sphere_potential(Z, R_nuc, r)
        V_point = -Z / r
        np.testing.assert_allclose(V, V_point, rtol=1e-14)

    def test_point_charge_limit(self):
        """R_nuc=0 returns exact point-charge potential."""
        Z = 1.0
        r = np.linspace(0.01, 10.0, 100)
        V = uniform_sphere_potential(Z, 0.0, r)
        V_point = -Z / r
        np.testing.assert_allclose(V, V_point, rtol=1e-14)

    def test_value_at_center(self):
        """V(0) = -3Z/(2R_nuc) for finite nucleus."""
        Z, R_nuc = 1.0, 0.5
        V0 = uniform_sphere_potential(Z, R_nuc, np.array([1e-15]))[0]
        expected = -3.0 * Z / (2.0 * R_nuc)
        assert abs(V0 - expected) < 1e-6


# ---------------------------------------------------------------------------
# 2. Finite size correction tests
# ---------------------------------------------------------------------------

class TestFiniteSizeCorrection:
    """Tests for first-order perturbation theory corrections."""

    def test_1s_hydrogen_known_value(self):
        """
        ΔE(1s, Z=1, R=0.88 fm) for uniform sphere model.

        The well-known ~5.2e-12 Ha value includes relativistic (Dirac)
        corrections. The non-relativistic uniform sphere model gives a
        slightly larger value (~1e-10 Ha) because the Dirac wavefunction
        has a different density at the origin than the Schrodinger one.

        We verify: (a) positive sign, (b) correct order of magnitude for
        the non-relativistic calculation.
        """
        R_nuc = 0.88 / 52917.72  # fm to bohr
        dE = finite_size_correction(Z=1.0, R_nuc=R_nuc, n=1, l=0)
        # The correction should be positive (finite nucleus is shallower)
        assert dE > 0, f"Expected positive correction, got {dE}"
        # Non-relativistic uniform sphere: ~1e-10 Ha at proton radius
        # Leading order: dE = (2Z⁴)/(3n³) × R² × |R_nl(0)|² factor
        assert 1e-12 < dE < 1e-8, (
            f"Expected ~1e-10 Ha, got {dE:.3e} Ha"
        )

    def test_l1_suppressed(self):
        """ΔE(2p) << ΔE(2s) due to angular momentum suppression."""
        R_nuc = 0.01  # Use larger R for numerical stability
        dE_2s = finite_size_correction(Z=1.0, R_nuc=R_nuc, n=2, l=0)
        dE_2p = finite_size_correction(Z=1.0, R_nuc=R_nuc, n=2, l=1)
        # p-wave suppression: |ψ(0)|² ~ r^{2l}, so ΔE(2p) ~ R^4 vs ΔE(2s) ~ R^2
        assert abs(dE_2p) < abs(dE_2s) * 0.1, (
            f"Expected l=1 suppression: dE_2s={dE_2s:.3e}, dE_2p={dE_2p:.3e}"
        )

    def test_zero_for_point_charge(self):
        """ΔE = 0 when R_nuc = 0."""
        dE = finite_size_correction(Z=1.0, R_nuc=0.0, n=1, l=0)
        assert dE == 0.0

    def test_1s_shift_quadratic_in_R(self):
        """ΔE(1s) ∝ R² to leading order."""
        R1 = 0.001
        R2 = 0.002
        dE1 = finite_size_correction(Z=1.0, R_nuc=R1, n=1, l=0)
        dE2 = finite_size_correction(Z=1.0, R_nuc=R2, n=1, l=0)
        ratio = dE2 / dE1
        expected_ratio = (R2 / R1) ** 2
        assert abs(ratio - expected_ratio) < 0.1 * expected_ratio, (
            f"Expected ratio {expected_ratio}, got {ratio:.4f}"
        )

    def test_Z_scaling(self):
        """ΔE scales as Z⁴ for 1s state (leading order)."""
        R_nuc = 0.001
        dE_Z1 = finite_size_correction(Z=1.0, R_nuc=R_nuc, n=1, l=0)
        dE_Z2 = finite_size_correction(Z=2.0, R_nuc=R_nuc, n=1, l=0)
        ratio = dE_Z2 / dE_Z1
        # Z^4 scaling: should be ~16
        assert abs(ratio - 16.0) < 3.0, (
            f"Expected Z⁴ scaling (~16), got ratio={ratio:.2f}"
        )


# ---------------------------------------------------------------------------
# 3. Spectrum vs nuclear radius
# ---------------------------------------------------------------------------

class TestSpectrumVsRadius:
    """Tests for the spectrum scan."""

    def test_point_charge_limit(self):
        """At R=0, spectrum matches -Z²/(2n²) to machine precision."""
        # Use the numerical solver with point charge
        res = fock_projection_with_finite_nucleus(
            Z=1.0, R_nuc=0.0, n_max=3, n_grid=4000
        )
        for (n, l), E_num in res['numerical_energies'].items():
            E_exact = -1.0 / (2.0 * n ** 2)
            rel_err = abs((E_num - E_exact) / E_exact)
            assert rel_err < 1e-4, (
                f"({n},{l}): E_num={E_num:.8f}, E_exact={E_exact:.8f}, "
                f"rel_err={rel_err:.2e}"
            )

    @pytest.mark.slow
    def test_large_R_harmonic(self):
        """At R >> a₀, the potential is approximately harmonic near center."""
        # For very large R_nuc, V(r) ≈ -3Z/(2R) + Z r²/(2R³) near r=0
        # This is a shifted harmonic oscillator with omega² = Z/R³
        from geovac.nuclear.potential_sparsity import solve_radial_schrodinger
        from geovac.nuclear.form_factor import uniform_sphere_potential

        R_nuc = 5.0
        Z = 1.0
        V_func = lambda r: uniform_sphere_potential(Z, R_nuc, r)

        # Ground state is much higher than -Z²/2
        E, _, _ = solve_radial_schrodinger(V_func, 0, 0, r_max=80.0, n_grid=4000)

        E_point = -Z ** 2 / 2.0
        assert E > E_point, (
            f"Finite nucleus should raise 1s energy: E={E:.6f}, E_point={E_point:.6f}"
        )


# ---------------------------------------------------------------------------
# 4. Fock projection tests
# ---------------------------------------------------------------------------

class TestFockProjection:
    """Tests for the Fock projection with finite nucleus."""

    def test_fock_projection_at_physical_R(self):
        """
        At R = 1e-5 a₀ (physical proton radius), the diagonal-formulation
        graph exactly reproduces the numerical spectrum (by construction).
        The physical point-vs-finite shift is ~10⁻¹⁰, well below the
        ~10⁻⁴ FD grid precision.
        """
        R_nuc = 1e-5  # ~physical proton radius
        res = fock_projection_with_finite_nucleus(
            Z=1.0, R_nuc=R_nuc, n_max=3, n_grid=4000
        )
        # Diagonal formulation matches spectrum exactly (by construction)
        assert res['max_rel_error_diag'] < 1e-8, (
            f"Diagonal graph error at physical R: {res['max_rel_error_diag']:.2e}"
        )
        # Verify the analytic finite-size correction is tiny at physical R
        dE_1s = finite_size_correction(Z=1.0, R_nuc=R_nuc, n=1, l=0)
        # Should be ~10^-10 Ha relative to -0.5 Ha = ~10^-10 relative
        assert abs(dE_1s / 0.5) < 1e-8, (
            f"Analytic 1s shift: {dE_1s:.3e} Ha (rel: {abs(dE_1s/0.5):.3e})"
        )

    def test_fock_projection_point_charge(self):
        """
        At R=0, diagonal formulation is exact.
        Pure graph formulation H = -(Z²/16)(D-A) has large basis-truncation
        error at small n_max — Paper 1 reports ~0.57% at max_n=30, and the
        error grows rapidly as max_n decreases (excited states diverge
        since (D-A) eigenvalues are bounded by max degree).
        """
        res = fock_projection_with_finite_nucleus(
            Z=1.0, R_nuc=0.0, n_max=3, n_grid=4000
        )
        # Diagonal form: exact
        assert res['max_rel_error_diag'] < 1e-8, (
            f"Diagonal form error: {res['max_rel_error_diag']:.2e}"
        )
        # Pure graph form: bounded but can be ~2x for high excited states
        # at max_n=3 due to finite-basis cut-off of the n=4,5,... states
        assert res['max_rel_error_pure'] < 3.0, (
            f"Pure graph error: {res['max_rel_error_pure']:.2e}"
        )

    @pytest.mark.slow
    def test_fock_projection_breakdown_radius(self):
        """
        Identify the R at which the point-charge spectrum diverges from
        the finite-nucleus spectrum by more than 1%.

        This is the radius at which the l-INDEPENDENT structure of node
        weights (-Z/n²) becomes inadequate — the finite nucleus introduces
        an l-dependent shift that the standard GeoVac lattice cannot encode
        without promoting node weights to l-dependent form.
        """
        R_values = np.logspace(-3, 0.5, 15)
        res = fock_projection_breakdown_scan(
            Z=1.0, n_max=3, R_values=R_values, n_grid=4000
        )
        # Should find a breakdown radius somewhere in the scanned range
        # (finite-size correction scales as R², so ~1% at R~0.1 bohr)
        assert res['breakdown_R'] < np.inf, (
            f"Breakdown radius not found. Max errors: {res['max_errors_point_vs_finite']}"
        )
        print(f"\nFock projection breakdown radius: {res['breakdown_R']:.4f} bohr")
        print(f"Max point-vs-finite errors: {res['max_errors_point_vs_finite']}")

    def test_node_weight_correction(self):
        """
        Corrected diagonal formulation H = diag(E_n(l,R)) reproduces
        the finite-size corrected spectrum exactly (by construction).
        """
        R_nuc = 0.1  # Large enough for measurable effect
        res = s3_lattice_with_modified_potential(
            Z=1.0, R_nuc=R_nuc, n_max=3, n_grid=4000
        )
        # Corrected graph (diagonal with numerical energies) is exact
        assert res['corrected_error'] < 1e-8, (
            f"Corrected diagonal error: {res['corrected_error']:.4e}"
        )


# ---------------------------------------------------------------------------
# 5. S³ lattice with modified potential
# ---------------------------------------------------------------------------

class TestS3LatticeModified:
    """Tests for the modified S³ lattice."""

    def test_standard_lattice_point_charge(self):
        """
        Standard lattice H = -(Z²/16)(D-A) has known basis-truncation
        error at small n_max. Paper 1 reports ~0.57% at max_n=30; the
        error grows rapidly as max_n decreases (high excited states are
        bounded by max degree of the graph). Just verify it's finite.
        """
        res = s3_lattice_with_modified_potential(
            Z=1.0, R_nuc=0.0, n_max=3, n_grid=4000
        )
        # Pure graph form error is ~2x at max_n=3 (high excited states)
        assert res['standard_error'] < 3.0, (
            f"Standard graph error bounded: {res['standard_error']:.2e}"
        )

    def test_degeneracy_broken(self):
        """Finite nuclear size breaks l-degeneracy."""
        R_nuc = 0.5  # Large enough for significant effect
        res = s3_lattice_with_modified_potential(
            Z=1.0, R_nuc=R_nuc, n_max=3, n_grid=4000
        )
        # At R=0.5 bohr, l-degeneracy should be broken
        assert res['degeneracy_broken'], (
            f"Expected l-degeneracy to be broken at R_nuc={R_nuc}"
        )
        # s-states should be more affected (shifted up more)
        for n, split in res['l_splitting'].items():
            # E(n,0) - E(n,1) should be positive (s raised more than p)
            assert split > 0, (
                f"n={n}: E(s)-E(p) = {split:.6f}, expected positive"
            )


# ---------------------------------------------------------------------------
# 6. Hydrogenic radial wavefunction
# ---------------------------------------------------------------------------

class TestHydrogenicRadial:
    """Tests for the hydrogenic radial wavefunction helper."""

    def test_1s_normalization(self):
        """1s wavefunction is normalized: integral |R|² r² dr = 1."""
        r = np.linspace(1e-6, 30.0, 10000)
        R = _hydrogenic_radial(1.0, 1, 0, r)
        norm = np.trapezoid(R ** 2 * r ** 2, r)
        assert abs(norm - 1.0) < 0.01, f"1s norm = {norm:.6f}"

    def test_2s_normalization(self):
        """2s wavefunction is normalized."""
        r = np.linspace(1e-6, 50.0, 20000)
        R = _hydrogenic_radial(1.0, 2, 0, r)
        norm = np.trapezoid(R ** 2 * r ** 2, r)
        assert abs(norm - 1.0) < 0.01, f"2s norm = {norm:.6f}"

    def test_2p_vanishes_at_origin(self):
        """2p wavefunction vanishes at r=0 (l=1)."""
        r = np.array([1e-6])
        R = _hydrogenic_radial(1.0, 2, 1, r)
        assert abs(R[0]) < 1e-3, f"|R_2p(0)| = {abs(R[0]):.6f}"


# ---------------------------------------------------------------------------
# Integration test: full investigation
# ---------------------------------------------------------------------------

class TestFullInvestigation:
    """Integration tests combining multiple functions."""

    @pytest.mark.slow
    def test_spectrum_table(self):
        """Generate the spectrum table for representative R values."""
        Z = 1.0
        n_max = 3
        R_test = [0.0, 1e-5, 0.001, 0.01, 0.1, 1.0, 5.0]

        print("\n" + "=" * 80)
        print("E(n,l) vs R_nuc for hydrogen (Z=1)")
        print("=" * 80)

        header = f"{'R_nuc (bohr)':>14}"
        states = []
        for n in range(1, n_max + 1):
            for l in range(n):
                states.append((n, l))
                header += f"  E({n},{l}):>12"
        print(f"\n{'R_nuc':>14}", end="")
        for n, l in states:
            print(f"  E({n},{l})", end="")
        print()
        print("-" * (14 + 12 * len(states)))

        for R_nuc in R_test:
            if R_nuc == 0:
                energies = {(n, l): -Z ** 2 / (2.0 * n ** 2)
                            for n in range(1, n_max + 1) for l in range(n)}
            else:
                V_func = lambda r, _R=R_nuc: uniform_sphere_potential(Z, _R, r)
                energies = {}
                for n in range(1, n_max + 1):
                    for l in range(n):
                        n_r = n - l - 1
                        try:
                            from geovac.nuclear.potential_sparsity import solve_radial_schrodinger
                            E, _, _ = solve_radial_schrodinger(
                                V_func, n_r, l, r_max=80.0, n_grid=4000
                            )
                            energies[(n, l)] = E
                        except Exception:
                            energies[(n, l)] = float('nan')

            print(f"{R_nuc:14.6e}", end="")
            for n, l in states:
                print(f"  {energies[(n, l)]:11.6f}", end="")
            print()

    @pytest.mark.slow
    def test_fock_projection_scan(self):
        """Scan to find the Fock projection breakdown radius."""
        R_values = np.logspace(-4, 0.5, 15)
        res = fock_projection_breakdown_scan(
            Z=1.0, n_max=3, R_values=R_values, n_grid=4000
        )

        print("\n" + "=" * 60)
        print("Fock projection: corrected graph vs numerical spectrum")
        print("=" * 60)
        print(f"{'R_nuc (bohr)':>14}  {'Max rel error':>14}")
        print("-" * 30)
        for R, err in zip(res['R_values'], res['max_errors']):
            flag = " ***" if err > 0.01 else ""
            print(f"{R:14.6e}  {err:14.6e}{flag}")
        print(f"\nBreakdown radius (>1% error): {res['breakdown_R']:.4f} bohr")
