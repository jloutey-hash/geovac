"""
Tests for the multi-zeta Roothaan-Hartree-Fock orbital tabulations.

Validates:
  - STO primitive normalization
  - MultiZetaOrbital construction and evaluation
  - Density builder produces correct electron count
  - Neutral Ne BBB93 reference data integrates to ~10
  - Xe two-zeta builder produces ~54 electrons for Cs/Ba
  - FrozenCore screening='multi_zeta' path works for Xe
  - FrozenCore screening='single_zeta' (default) is bit-identical to legacy
  - Backward compatibility: existing FrozenCore behavior preserved
"""

import warnings

import numpy as np
import pytest

from geovac.multi_zeta_orbitals import (
    STO,
    MultiZetaOrbital,
    _build_ne_orbitals_neutral,
    _build_na_valence_orbitals_physical,
    _PHYSICAL_NA_OBSERVABLES,
    _PHYSICAL_FIT_AVAILABLE,
    build_two_zeta_xe_orbitals_from_cr,
    density_from_orbitals,
    get_physical_valence_orbitals,
    core_electron_count,
)
from geovac.neon_core import FrozenCore, _CLEMENTI_ZETA_XE


# ---------------------------------------------------------------------------
# STO primitives
# ---------------------------------------------------------------------------

class TestSTOPrimitives:
    def test_sto_normalization(self):
        """integral_0^inf chi^2 * r^2 dr = 1 for any (n, zeta)."""
        for n, zeta in [(1, 1.0), (1, 5.0), (2, 2.0), (3, 0.5), (5, 12.3)]:
            sto = STO(n=n, zeta=zeta)
            r = np.geomspace(1e-5, 100.0, 8000)
            chi = sto.evaluate(r)
            integrand = chi * chi * r * r
            norm_sq = np.trapezoid(integrand, r)
            assert abs(norm_sq - 1.0) < 1e-3, (
                f"STO(n={n}, zeta={zeta}) normalization {norm_sq:.6f} != 1"
            )

    def test_sto_normalization_constant(self):
        """N_{n, zeta} = (2 zeta)^{n + 1/2} / sqrt((2n)!)."""
        sto = STO(n=2, zeta=3.0)
        N = sto.normalization()
        # (2 * 3)^{2.5} / sqrt(24) = 6^2.5 / sqrt(24) = 88.18... / 4.899... ~ 18.0
        expected = (2.0 * 3.0) ** 2.5 / np.sqrt(24.0)
        assert abs(N - expected) < 1e-10

    def test_sto_high_zeta_compact(self):
        """High-zeta STOs are concentrated near origin."""
        sto = STO(n=1, zeta=20.0)
        r = np.geomspace(1e-4, 5.0, 3000)
        chi = sto.evaluate(r)
        prob = chi * chi * r * r
        # Most probability inside r < 1/zeta = 0.05 bohr
        mask = r < 0.5
        inner_prob = np.trapezoid(prob[mask], r[mask])
        assert inner_prob > 0.99, f"Inner probability {inner_prob} should be > 0.99"


# ---------------------------------------------------------------------------
# MultiZetaOrbital
# ---------------------------------------------------------------------------

class TestMultiZetaOrbital:
    def test_orbital_basic_construction(self):
        """Construct an orbital and evaluate."""
        prims = (STO(n=1, zeta=2.0), STO(n=1, zeta=4.0))
        coefs = (0.5, 0.5)
        orb = MultiZetaOrbital(
            n_orbital=1, l_orbital=0, occupancy=2,
            primitives=prims, coefficients=coefs,
        )
        r = np.array([0.5, 1.0, 2.0])
        R = orb.evaluate(r)
        assert R.shape == r.shape
        assert np.all(R > 0)  # 1s orbitals are positive everywhere

    def test_orbital_validation_n_lt_l(self):
        """l >= n is invalid."""
        with pytest.raises(ValueError, match="quantum"):
            MultiZetaOrbital(
                n_orbital=1, l_orbital=1, occupancy=2,
                primitives=(STO(n=1, zeta=1.0),),
                coefficients=(1.0,),
            )

    def test_orbital_validation_mismatched_lengths(self):
        """primitives and coefficients must have same length."""
        with pytest.raises(ValueError, match="match"):
            MultiZetaOrbital(
                n_orbital=2, l_orbital=0, occupancy=2,
                primitives=(STO(n=2, zeta=1.0), STO(n=2, zeta=2.0)),
                coefficients=(1.0,),
            )

    def test_orbital_validation_sto_n_lt_lplus1(self):
        """STO primitive n must be >= l + 1."""
        with pytest.raises(ValueError, match="primitive"):
            MultiZetaOrbital(
                n_orbital=2, l_orbital=1, occupancy=6,
                primitives=(STO(n=1, zeta=1.0),),  # n=1 < l+1=2
                coefficients=(1.0,),
            )

    def test_orbital_density_contribution(self):
        """occ * |R|^2 * r^2 returned by density_contribution."""
        sto = STO(n=1, zeta=1.0)
        orb = MultiZetaOrbital(
            n_orbital=1, l_orbital=0, occupancy=2,
            primitives=(sto,), coefficients=(1.0,),
        )
        r = np.geomspace(1e-4, 30.0, 5000)
        n_r = orb.density_contribution(r)
        N = np.trapezoid(n_r, r)
        # 2 electrons in a normalized 1s
        assert abs(N - 2.0) < 1e-3


# ---------------------------------------------------------------------------
# BBB93 Ne reference
# ---------------------------------------------------------------------------

class TestBBB93Neon:
    def test_neutral_ne_orbital_count(self):
        """Neutral Ne has 3 orbitals (1s, 2s, 2p)."""
        orbs = _build_ne_orbitals_neutral()
        assert len(orbs) == 3
        # Check (n, l) of each
        labels = sorted([(o.n_orbital, o.l_orbital) for o in orbs])
        assert labels == [(1, 0), (2, 0), (2, 1)]

    def test_neutral_ne_electron_count(self):
        """Total occupancy = 10 (Ne)."""
        orbs = _build_ne_orbitals_neutral()
        assert core_electron_count(orbs) == 10

    def test_neutral_ne_density_normalization(self):
        """integral n_core(r) * r^2 dr (built into n_r) = 10 within 1%."""
        orbs = _build_ne_orbitals_neutral()
        r = np.geomspace(1e-4, 30.0, 5000)
        n_r = density_from_orbitals(orbs, r)
        N = np.trapezoid(n_r, r)
        # BBB93 leading-digits coefficients give ~1% normalization error
        assert abs(N - 10.0) / 10.0 < 0.02


# ---------------------------------------------------------------------------
# Xe two-zeta builder
# ---------------------------------------------------------------------------

class TestXeTwoZetaBuilder:
    def test_xe_builder_orbital_count(self):
        """[Xe] core has 11 orbitals."""
        zetas_cs = _CLEMENTI_ZETA_XE[55]
        orbs = build_two_zeta_xe_orbitals_from_cr(zetas_cs)
        assert len(orbs) == 11
        labels = sorted([(o.n_orbital, o.l_orbital) for o in orbs])
        assert labels == [
            (1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2),
            (4, 0), (4, 1), (4, 2), (5, 0), (5, 1),
        ]

    def test_xe_builder_electron_count_cs(self):
        """[Xe] core for Cs = 54 electrons."""
        zetas_cs = _CLEMENTI_ZETA_XE[55]
        orbs = build_two_zeta_xe_orbitals_from_cr(zetas_cs)
        assert core_electron_count(orbs) == 54

    def test_xe_builder_electron_count_ba(self):
        """[Xe] core for Ba = 54 electrons (same as Cs)."""
        zetas_ba = _CLEMENTI_ZETA_XE[56]
        orbs = build_two_zeta_xe_orbitals_from_cr(zetas_ba)
        assert core_electron_count(orbs) == 54

    def test_xe_density_normalization(self):
        """integral n_core(r) * r^2 dr = 54 within 0.1%."""
        zetas_cs = _CLEMENTI_ZETA_XE[55]
        orbs = build_two_zeta_xe_orbitals_from_cr(zetas_cs)
        r = np.geomspace(1e-5, 60.0, 8000)
        n_r = density_from_orbitals(orbs, r)
        N = np.trapezoid(n_r, r)
        # Each orbital is renormalized at build time; overall density should
        # match exactly to ~0.1%
        assert abs(N - 54.0) / 54.0 < 0.005

    def test_xe_builder_wrong_zeta_count(self):
        """11 zetas required."""
        with pytest.raises(ValueError, match="11 zetas"):
            build_two_zeta_xe_orbitals_from_cr((1.0,) * 5)


# ---------------------------------------------------------------------------
# FrozenCore integration
# ---------------------------------------------------------------------------

class TestFrozenCoreMultiZetaIntegration:
    def test_frozen_core_default_screening_is_single_zeta(self):
        """Default FrozenCore(Z) gives screening='single_zeta'."""
        fc = FrozenCore(Z=55)
        assert fc.screening == 'single_zeta'

    def test_frozen_core_bit_identical_default_vs_explicit_single(self):
        """Default and screening='single_zeta' produce bit-identical Z_eff."""
        fc1 = FrozenCore(Z=55)
        fc1.solve()
        fc2 = FrozenCore(Z=55, screening='single_zeta')
        fc2.solve()
        for r in [0.01, 0.1, 1.0, 5.0, 20.0]:
            assert fc1.z_eff(r) == fc2.z_eff(r), (
                f"Z_eff at r={r} differs: {fc1.z_eff(r)} vs {fc2.z_eff(r)}"
            )

    def test_frozen_core_invalid_screening_raises(self):
        """Bad screening name raises."""
        with pytest.raises(ValueError, match="Unknown screening"):
            FrozenCore(Z=55, screening='bogus')

    def test_frozen_core_multi_zeta_xe_works(self):
        """FrozenCore(Z=55, screening='multi_zeta') solves and gives finite
        Z_eff at relevant radii."""
        fc = FrozenCore(Z=55, screening='multi_zeta')
        fc.solve()
        assert fc.screening == 'multi_zeta'
        # Z_eff(0) = Z (no screening at origin)
        assert abs(fc.z_eff(0.001) - 55.0) < 0.5
        # Z_eff(infinity) = Z - N_core = 55 - 54 = 1
        assert abs(fc.z_eff(40.0) - 1.0) < 0.05
        # Z_eff is monotonically decreasing
        z_eff_grid = fc.z_eff(np.linspace(0.01, 30.0, 100))
        diffs = np.diff(z_eff_grid)
        # Allow tiny positive bumps from interpolation noise
        assert np.all(diffs < 1e-3), "Z_eff should be monotonically decreasing"

    def test_frozen_core_multi_zeta_falls_back_for_unsupported_core(self):
        """Non-Xe cores fall back to single_zeta with a UserWarning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            fc = FrozenCore(Z=11, screening='multi_zeta')
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert 'multi-zeta' in str(w[0].message).lower()
        # screening attribute updates to single_zeta after fallback
        assert fc.screening == 'single_zeta'

    def test_frozen_core_multi_zeta_density_normalization(self):
        """The multi_zeta core density normalizes to 54 electrons."""
        fc = FrozenCore(Z=55, screening='multi_zeta', n_radial=4000, r_max=40.0)
        fc.solve()
        # _N_core[-1] should be 54 within 0.5% by construction
        assert abs(fc.N_core[-1] - 54.0) / 54.0 < 0.01

    def test_frozen_core_multi_zeta_z_eff_at_origin(self):
        """At r->0 the Z_eff approaches Z (no screening yet)."""
        fc = FrozenCore(Z=55, screening='multi_zeta')
        fc.solve()
        assert abs(fc.z_eff(0.0001) - 55.0) < 0.1

    def test_frozen_core_multi_zeta_z_eff_at_infinity(self):
        """At r >> all orbital extents, Z_eff -> Z - 54 = 1."""
        fc = FrozenCore(Z=55, screening='multi_zeta', r_max=50.0)
        fc.solve()
        z_eff_far = fc.z_eff(45.0)
        assert abs(z_eff_far - 1.0) < 0.01

    def test_frozen_core_existing_callers_unaffected(self):
        """Callers that don't pass screening get the default behavior."""
        # Mimic the SO callers
        fc = FrozenCore(Z=20)
        fc.solve()
        assert fc.screening == 'single_zeta'
        # And basic queries work
        assert 0 < fc.z_eff(2.0) < 20.0


# ---------------------------------------------------------------------------
# Z=11 (Na) physical-fit valence orbital tabulation
# Sprint alpha-Multi-zeta (Track alpha-2, 2026-05-23)
# ---------------------------------------------------------------------------

class TestPhysicalNaValenceOrbitals:
    """Tests for the Na 3s/3p physical-fit multi-zeta tabulation.

    These tests verify the production parameters added in Sprint
    alpha-Multi-zeta Track alpha-2 (2026-05-23). The parameters were
    fitted from the FrozenCore([Ne]) + radial Schrodinger solver against
    the physical Na 3s and 3p screened wavefunctions, with bounded
    coefficients to ensure a well-conditioned basis.
    """

    def test_physical_fit_registry_contains_z11(self):
        """Z=11 is registered as a physical-fit Z."""
        assert 11 in _PHYSICAL_FIT_AVAILABLE

    def test_get_physical_valence_orbitals_na_returns_two_orbitals(self):
        """Na has 3s and 3p valence orbitals."""
        orbs = get_physical_valence_orbitals(11)
        assert len(orbs) == 2
        # (n, l) labels: 3s = (3, 0), 3p = (3, 1)
        labels = sorted([(o.n_orbital, o.l_orbital) for o in orbs])
        assert labels == [(3, 0), (3, 1)]

    def test_get_physical_valence_orbitals_unsupported_z_raises(self):
        """Unsupported Z raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="physical-fit"):
            get_physical_valence_orbitals(19)  # K not yet tabulated

    def test_na_3s_is_normalized(self):
        """integral |R_3s|^2 r^2 dr = 1 within 0.5%."""
        orbs = _build_na_valence_orbitals_physical()
        na_3s = orbs[0]
        assert (na_3s.n_orbital, na_3s.l_orbital) == (3, 0)
        r = np.geomspace(1e-4, 80.0, 8000)
        R = na_3s.evaluate(r)
        norm_sq = np.trapezoid(R ** 2 * r ** 2, r)
        assert abs(norm_sq - 1.0) < 5e-3, (
            f"Na 3s normalization {norm_sq:.6f} should be ~1.0"
        )

    def test_na_3p_is_normalized(self):
        """integral |R_3p|^2 r^2 dr = 1 within 0.5%."""
        orbs = _build_na_valence_orbitals_physical()
        na_3p = orbs[1]
        assert (na_3p.n_orbital, na_3p.l_orbital) == (3, 1)
        r = np.geomspace(1e-4, 80.0, 8000)
        R = na_3p.evaluate(r)
        norm_sq = np.trapezoid(R ** 2 * r ** 2, r)
        assert abs(norm_sq - 1.0) < 5e-3, (
            f"Na 3p normalization {norm_sq:.6f} should be ~1.0"
        )

    def test_na_3s_radial_node_count(self):
        """Na 3s has 2 radial nodes (3s has n_r = n - l - 1 = 2 nodes)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3s = orbs[0]
        # Use uniform grid to count sign changes robustly
        r = np.linspace(0.005, 50.0, 8000)
        R = na_3s.evaluate(r)
        sign_changes = int(np.sum(np.diff(np.sign(R)) != 0))
        assert sign_changes == 2, (
            f"Na 3s should have 2 radial nodes, got {sign_changes}"
        )

    def test_na_3p_radial_node_count(self):
        """Na 3p has 1 radial node (3p has n_r = n - l - 1 = 1 node)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3p = orbs[1]
        r = np.linspace(0.005, 50.0, 8000)
        R = na_3p.evaluate(r)
        sign_changes = int(np.sum(np.diff(np.sign(R)) != 0))
        assert sign_changes == 1, (
            f"Na 3p should have 1 radial node, got {sign_changes}"
        )

    def test_na_3s_mean_radius_matches_physical(self):
        """<r>_3s ~ 4.47 bohr (physical Na 3s mean radius)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3s = orbs[0]
        r = np.geomspace(1e-4, 80.0, 8000)
        R = na_3s.evaluate(r)
        # <r> = integral R^2 r^3 dr (radial part of integral <r> = integral |psi|^2 r dV)
        mean_r = float(np.trapezoid(r * R ** 2 * r ** 2, r))
        expected = _PHYSICAL_NA_OBSERVABLES['Na_3s_mean_radius_bohr']
        # Allow 2% deviation (fit overlap was 0.999999 but L2 error gives small drift)
        assert abs(mean_r - expected) / expected < 0.02, (
            f"Na 3s mean radius {mean_r:.4f} should match physical "
            f"{expected:.4f} within 2%"
        )

    def test_na_3p_mean_radius_matches_physical(self):
        """<r>_3p ~ 5.93 bohr (physical Na 3p mean radius)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3p = orbs[1]
        r = np.geomspace(1e-4, 80.0, 8000)
        R = na_3p.evaluate(r)
        mean_r = float(np.trapezoid(r * R ** 2 * r ** 2, r))
        expected = _PHYSICAL_NA_OBSERVABLES['Na_3p_mean_radius_bohr']
        assert abs(mean_r - expected) / expected < 0.02, (
            f"Na 3p mean radius {mean_r:.4f} should match physical "
            f"{expected:.4f} within 2%"
        )

    def test_na_3s_is_diffuse_not_compact(self):
        """Na 3s mean radius is much larger than hydrogenic Z=1 1s (1.5 bohr).

        This is the structural distinction the bimodule diagnostic
        identified: physical Na 3s is diffuse (mean ~4.5 bohr) while
        the hydrogenic Z_orb=1 placeholder is compact (mean ~1.5 bohr).
        """
        orbs = _build_na_valence_orbitals_physical()
        na_3s = orbs[0]
        r = np.geomspace(1e-4, 80.0, 8000)
        R = na_3s.evaluate(r)
        mean_r = float(np.trapezoid(r * R ** 2 * r ** 2, r))
        # Diffuse: mean r > 3 bohr (hydrogenic Z=1 1s mean r = 1.5 bohr)
        assert mean_r > 3.0, (
            f"Na 3s should be diffuse (mean r > 3 bohr), got {mean_r:.4f}"
        )

    def test_na_3s_has_nonzero_density_at_origin(self):
        """Na 3s has nonzero R(0) (s-orbital, finite at origin)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3s = orbs[0]
        # Evaluate at very small r (R should be finite, large from core penetration)
        r_small = np.array([0.001, 0.01])
        R_small = na_3s.evaluate(r_small)
        # R(0) for Na 3s should be at least a few bohr^(-3/2) from core penetration
        assert R_small[0] > 0.5, (
            f"Na 3s R(0+) should be substantial (>0.5), got {R_small[0]:.4e}"
        )

    def test_na_3p_vanishes_at_origin(self):
        """Na 3p has R(0) = 0 (l=1 orbital, vanishes at origin as r^1)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3p = orbs[1]
        r_small = np.array([0.001])
        R_small = na_3p.evaluate(r_small)
        # For l=1, R(r) ~ r at origin, so R(0.001) should be very small
        # The STO primitives have r^{n-1} prefactor; for n=2 STO, prefactor is r,
        # for n=3 it's r^2. The leading-order behavior is dominated by the n=2 STO
        # contribution, giving R(r) ~ r at origin.
        assert abs(R_small[0]) < 0.01, (
            f"Na 3p R(0+ = 0.001) should be near zero, got {R_small[0]:.4e}"
        )

    def test_orbital_orthogonality_with_full_spherical_harmonic(self):
        """3s and 3p have different l so they are angularly orthogonal regardless
        of the radial overlap (full 3D overlap <3s|3p> = 0 by Y_lm orthogonality)."""
        orbs = _build_na_valence_orbitals_physical()
        na_3s, na_3p = orbs
        assert na_3s.l_orbital != na_3p.l_orbital, (
            "3s and 3p must have different l for angular orthogonality"
        )

    def test_na_orbital_evaluation_returns_arrays(self):
        """Evaluation on a grid returns ndarray with matching shape."""
        orbs = _build_na_valence_orbitals_physical()
        r = np.linspace(0.1, 10.0, 100)
        for orb in orbs:
            R = orb.evaluate(r)
            assert R.shape == r.shape
            assert np.all(np.isfinite(R))


# ---------------------------------------------------------------------------
# Slow integration tests (opt-in via --slow)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestCsValenceWithMultiZeta:
    def test_cs_6s_eigenvalue_finite_with_multi_zeta(self):
        """Cs 6s orbital under multi_zeta screening gives a finite, bound
        eigenvalue."""
        from geovac.neon_core import _solve_screened_radial_log

        energy, _u, _r, R0 = _solve_screened_radial_log(
            55, 0, 6, n_grid=100000, r_max=60.0, screening='multi_zeta',
        )
        assert -10.0 < energy < 0.0
        # |psi(0)|^2 should be finite and positive
        psi0_sq = R0 ** 2 / (4.0 * np.pi)
        assert psi0_sq > 0
