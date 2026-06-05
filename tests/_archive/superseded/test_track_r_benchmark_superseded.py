"""
Track R benchmark: Compare algebraic vs quadrature defaults for Z_eff and exchange.

Tests that switching the production defaults produces negligible differences
in composed-geometry PES results for LiH, BeH₂, and H₂O.
"""

import time
import numpy as np
import pytest
import json
from pathlib import Path

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.composed_diatomic import _v_cross_nuc_1s


# ======================================================================
# Z_eff consistency tests
# ======================================================================

class TestZeffConsistency:
    """Verify algebraic Z_eff matches spline Z_eff for both Li and Be cores."""

    @pytest.fixture(scope='class')
    def li_spline(self):
        cs = CoreScreening(Z=3, l_max=2, n_alpha=200, zeff_method='spline')
        cs.solve()
        return cs

    @pytest.fixture(scope='class')
    def li_spectral(self):
        cs = CoreScreening(Z=3, l_max=2, n_alpha=200, zeff_method='spectral_laguerre')
        cs.solve()
        return cs

    @pytest.fixture(scope='class')
    def be_spline(self):
        cs = CoreScreening(Z=4, l_max=2, n_alpha=200, zeff_method='spline')
        cs.solve()
        return cs

    @pytest.fixture(scope='class')
    def be_spectral(self):
        cs = CoreScreening(Z=4, l_max=2, n_alpha=200, zeff_method='spectral_laguerre')
        cs.solve()
        return cs

    def test_li_z_eff_consistency(self, li_spline, li_spectral):
        """Li Z_eff: spectral and spline agree to < 0.01 at key r-points."""
        r_pts = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
        z_spline = li_spline.z_eff(r_pts)
        z_spectral = li_spectral.z_eff(r_pts)
        max_diff = np.max(np.abs(z_spline - z_spectral))
        print(f"  Li Z_eff max diff: {max_diff:.6f}")
        assert max_diff < 0.01, f"Li Z_eff divergence: max diff = {max_diff}"

    def test_be_z_eff_consistency(self, be_spline, be_spectral):
        """Be Z_eff: spectral falls back to spline if fit fails, giving exact agreement."""
        r_pts = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
        z_spline = be_spline.z_eff(r_pts)
        z_spectral = be_spectral.z_eff(r_pts)
        max_diff = np.max(np.abs(z_spline - z_spectral))
        print(f"  Be Z_eff max diff: {max_diff:.6f}")
        # Be (Z=4) has a tight density that causes polynomial cancellation
        # in the spectral form. The CoreScreening validation catches this
        # and falls back to spline, so the results should agree exactly.
        assert max_diff < 0.01, f"Be Z_eff divergence: max diff = {max_diff}"

    def test_be_spectral_falls_back(self, be_spectral):
        """Be spectral should fall back to spline due to tight density."""
        # The spectral fit for Be Z=4 fails validation and falls back
        assert be_spectral.zeff_method == 'spline', (
            f"Expected spline fallback for Be, got '{be_spectral.zeff_method}'"
        )

    def test_li_energy_same(self, li_spline, li_spectral):
        """Core energy should be identical (Z_eff doesn't affect core solve)."""
        assert abs(li_spline.energy - li_spectral.energy) < 1e-10

    def test_li_pk_same(self, li_spline, li_spectral):
        """PK parameters derived from spline vs spectral cores should match."""
        pk_s = AbInitioPK(li_spline, n_core=2)
        pk_a = AbInitioPK(li_spectral, n_core=2)
        assert abs(pk_s.A - pk_a.A) < 0.01
        assert abs(pk_s.B - pk_a.B) < 0.01


# ======================================================================
# LiH single-point energy comparison
# ======================================================================

class TestLiHSinglePoint:
    """Compare LiH composed energy at R_eq with both Z_eff methods."""

    @pytest.fixture(scope='class')
    def lih_energies(self):
        """Compute LiH composed energy at R=3.015 with both methods."""
        R = 3.015
        results = {}

        for method in ['spline', 'spectral_laguerre']:
            t0 = time.time()
            core = CoreScreening(Z=3, l_max=2, n_alpha=200,
                                zeff_method=method)
            core.solve()
            E_core = core.energy

            pk = AbInitioPK(core, n_core=2)
            pk_d = pk.pk_dict(atom='A')
            pk_d['channel_mode'] = 'l_dependent'

            result = solve_level4_h2_multichannel(
                R=R, Z_A=1.0, Z_B=1.0, l_max=2,
                n_alpha=100, n_Re=300, verbose=False,
                pk_potentials=[pk_d],
            )

            V_NN = 3.0 / R
            V_cross = _v_cross_nuc_1s(3.0, 2, 1.0, R)
            E_composed = E_core + V_cross + result['E_elec'] + V_NN
            dt = time.time() - t0

            results[method] = {
                'E_core': E_core,
                'E_elec': result['E_elec'],
                'E_composed': E_composed,
                'pk_A': pk.A,
                'pk_B': pk.B,
                'time': dt,
            }

        return results

    def test_energy_agreement(self, lih_energies):
        """LiH composed energy at R_eq should agree to < 0.001 Ha."""
        e_spline = lih_energies['spline']['E_composed']
        e_spectral = lih_energies['spectral_laguerre']['E_composed']
        diff = abs(e_spline - e_spectral)
        print(f"\n  LiH R=3.015:")
        print(f"    Spline:   E = {e_spline:.6f} Ha")
        print(f"    Spectral: E = {e_spectral:.6f} Ha")
        print(f"    |Delta|:  {diff:.6f} Ha")
        assert diff < 0.001, f"LiH energy difference too large: {diff}"

    def test_core_energy_identical(self, lih_energies):
        """E_core must be identical regardless of zeff_method."""
        diff = abs(lih_energies['spline']['E_core'] -
                   lih_energies['spectral_laguerre']['E_core'])
        assert diff < 1e-10

    def test_pk_parameters_identical(self, lih_energies):
        """PK parameters should be the same."""
        diff_A = abs(lih_energies['spline']['pk_A'] -
                     lih_energies['spectral_laguerre']['pk_A'])
        diff_B = abs(lih_energies['spline']['pk_B'] -
                     lih_energies['spectral_laguerre']['pk_B'])
        print(f"  PK A diff: {diff_A:.6f}, B diff: {diff_B:.6f}")
        assert diff_A < 0.01
        assert diff_B < 0.01

    def test_print_comparison(self, lih_energies):
        """Print formatted comparison table."""
        print("\n  LiH Single-Point Comparison (R=3.015):")
        print(f"  {'Metric':<20} {'Spline':>12} {'Spectral':>12} {'Delta':>12}")
        print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12}")
        for key in ['E_core', 'E_elec', 'E_composed', 'pk_A', 'pk_B', 'time']:
            sv = lih_energies['spline'][key]
            av = lih_energies['spectral_laguerre'][key]
            delta = av - sv
            print(f"  {key:<20} {sv:>12.6f} {av:>12.6f} {delta:>+12.6f}")


# ======================================================================
# LiH PES comparison (slow but comprehensive)
# ======================================================================

@pytest.mark.slow
class TestLiHPES:
    """Full PES comparison for LiH."""

    @pytest.fixture(scope='class')
    def lih_pes_results(self):
        """Run full PES with both methods."""
        from geovac.composed_diatomic import ComposedDiatomicSolver

        results = {}
        R_grid = np.concatenate([
            np.linspace(2.0, 2.5, 3),
            np.linspace(2.7, 4.0, 10),
            np.linspace(4.5, 7.0, 5),
        ])

        for label, zeff_method in [('quadrature', 'spline'),
                                     ('algebraic', 'spectral_laguerre')]:
            t0 = time.time()
            solver = ComposedDiatomicSolver.LiH_ab_initio(
                l_max=2, pk_channel_mode='l_dependent', verbose=False,
            )

            # Override core with specified zeff_method
            solver.core = CoreScreening(Z=3, l_max=2, n_alpha=200,
                                        zeff_method=zeff_method)
            solver.core.solve()
            solver.E_core = solver.core.energy
            solver.ab_initio_pk = AbInitioPK(solver.core, n_core=2)
            solver.pk_A = solver.ab_initio_pk.A
            solver.pk_B = solver.ab_initio_pk.B
            pk_d = solver.ab_initio_pk.pk_dict(atom='A')
            pk_d['channel_mode'] = 'l_dependent'
            solver.pk_potentials = [pk_d]

            pes = solver.scan_pes(R_grid=R_grid)
            spectro = solver.fit_spectroscopic_constants()

            results[label] = {
                'R_eq': spectro['R_eq'],
                'E_min': spectro['E_min'],
                'D_e': spectro['D_e'],
                'omega_e': spectro['omega_e'],
                'total_time': time.time() - t0,
            }

        return results

    def test_r_eq_agreement(self, lih_pes_results):
        """R_eq should agree to < 0.5%."""
        q = lih_pes_results['quadrature']['R_eq']
        a = lih_pes_results['algebraic']['R_eq']
        pct = abs(a - q) / q * 100
        print(f"\n  LiH R_eq: quad={q:.4f}, alg={a:.4f}, delta={pct:.3f}%")
        assert pct < 0.5, f"R_eq shift too large: {pct:.3f}%"

    def test_d_e_agreement(self, lih_pes_results):
        """D_e should agree to < 1%."""
        q = lih_pes_results['quadrature']['D_e']
        a = lih_pes_results['algebraic']['D_e']
        pct = abs(a - q) / q * 100
        print(f"  LiH D_e: quad={q:.6f}, alg={a:.6f}, delta={pct:.3f}%")
        assert pct < 1.0, f"D_e shift too large: {pct:.3f}%"

    def test_print_comparison(self, lih_pes_results):
        """Print formatted PES comparison table."""
        print("\n  LiH PES Comparison:")
        print(f"  {'Metric':<15} {'Quadrature':>12} {'Algebraic':>12} {'Delta':>10} {'Delta%':>8}")
        print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*10} {'-'*8}")
        for key in ['R_eq', 'E_min', 'D_e', 'omega_e', 'total_time']:
            qv = lih_pes_results['quadrature'][key]
            av = lih_pes_results['algebraic'][key]
            delta = av - qv
            pct = delta / abs(qv) * 100 if abs(qv) > 1e-15 else 0.0
            print(f"  {key:<15} {qv:>12.6f} {av:>12.6f} {delta:>+10.6f} {pct:>+7.3f}%")


# ======================================================================
# Exchange method consistency (algebraic_laguerre vs numerical)
# ======================================================================

class TestExchangeMethodConsistency:
    """Verify algebraic Slater F^0 matches numerical for inter-fiber coupling."""

    def test_slater_f0_consistency(self):
        """F^0 from algebraic_laguerre should match numerical to < 5%."""
        from geovac.inter_fiber_coupling import slater_fk_integral
        from geovac.algebraic_slater import slater_fk_integral_algebraic

        # Create a simple test density (1s hydrogenic)
        Z = 2.0
        n_r = 500
        r_max = 10.0
        dr = r_max / n_r
        r = (np.arange(n_r) + 0.5) * dr
        P = 4.0 * Z**3 * r**2 * np.exp(-2 * Z * r)

        f0_num = slater_fk_integral(r, P, P, k=0)
        f0_alg = slater_fk_integral_algebraic(r, P, P, k=0)

        rel_diff = abs(f0_alg - f0_num) / abs(f0_num) * 100
        print(f"\n  F^0 (1s Z={Z}): numerical={f0_num:.6f}, algebraic={f0_alg:.6f}, "
              f"diff={rel_diff:.2f}%")
        # Known F^0(1s,1s) = 5Z/8 = 1.25
        assert rel_diff < 5.0, f"F^0 disagreement too large: {rel_diff:.2f}%"


# ======================================================================
# Default verification
# ======================================================================

class TestNewDefaults:
    """Verify the new default settings are active."""

    def test_core_screening_default(self):
        """CoreScreening default should now be 'spectral_laguerre'."""
        cs = CoreScreening(Z=3)
        assert cs.zeff_method == 'spectral_laguerre'

    def test_slater_method_default(self):
        """full_exchange_inter_fiber_energy default should be 'algebraic_laguerre'."""
        import inspect
        from geovac.inter_fiber_coupling import full_exchange_inter_fiber_energy
        sig = inspect.signature(full_exchange_inter_fiber_energy)
        default = sig.parameters['slater_method'].default
        assert default == 'algebraic_laguerre', f"Expected 'algebraic_laguerre', got '{default}'"

    def test_monopole_method_default(self):
        """monopole_inter_fiber_energy default should be 'algebraic_laguerre'."""
        import inspect
        from geovac.inter_fiber_coupling import monopole_inter_fiber_energy
        sig = inspect.signature(monopole_inter_fiber_energy)
        default = sig.parameters['slater_method'].default
        assert default == 'algebraic_laguerre', f"Expected 'algebraic_laguerre', got '{default}'"
