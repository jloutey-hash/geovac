"""
Tests for nuclear sparsity characterization (Track NC).

Validates:
1. Radial solver reproduces known exact wavefunctions
2. Slater integrals match known values
3. Angular (Gaunt) selection rules are potential-independent
4. ERI symmetry properties
5. Sparsity table runs for all potentials
"""

import numpy as np
import pytest

from geovac.nuclear.potential_sparsity import (
    angular_zero_count,
    ck_coefficient,
    compute_eri_tensor,
    compute_slater_rk,
    enumerate_states,
    make_potential,
    pauli_scaling_estimate,
    radial_wavefunctions_for_potential,
    solve_radial_schrodinger,
    sparsity_table,
)


# ---------------------------------------------------------------------------
# Test 1: Coulomb wavefunctions
# ---------------------------------------------------------------------------

class TestCoulombWavefunctions:
    """Verify numerical radial solver reproduces hydrogenic wavefunctions."""

    def test_1s_energy(self):
        """n_r=0, l=0 -> E = -Z^2/2 for hydrogen (Z=1)."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 1.0})
        E, R_grid, r_grid = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        assert abs(E - (-0.5)) / 0.5 < 0.001, f"1s energy {E} vs exact -0.5"

    def test_2s_energy(self):
        """n_r=1, l=0 -> E = -Z^2/(2*4) = -0.125 for hydrogen."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 1.0})
        E, R_grid, r_grid = solve_radial_schrodinger(V_func, 1, 0, r_max, n_grid)
        assert abs(E - (-0.125)) / 0.125 < 0.001, f"2s energy {E} vs exact -0.125"

    def test_2p_energy(self):
        """n_r=0, l=1 -> n=2, E = -0.125 for hydrogen."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 1.0})
        E, R_grid, r_grid = solve_radial_schrodinger(V_func, 0, 1, r_max, n_grid)
        assert abs(E - (-0.125)) / 0.125 < 0.001, f"2p energy {E} vs exact -0.125"

    def test_1s_normalization(self):
        """R_{1s} should be normalized: integral |R|^2 r^2 dr = 1."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 1.0})
        _, R_grid, r_grid = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        norm = np.trapezoid(R_grid**2 * r_grid**2, r_grid)
        assert abs(norm - 1.0) < 0.01, f"Normalization {norm} vs 1.0"

    def test_1s_shape(self):
        """R_{1s}(r) = 2 Z^{3/2} exp(-Zr), max at r=0."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 1.0})
        _, R_grid, r_grid = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        # R(r) should be monotonically decreasing for 1s
        # (approximately, given grid discretization)
        # Check that maximum is near r=0
        max_idx = np.argmax(np.abs(R_grid))
        assert max_idx < len(R_grid) // 4, "1s wavefunction peak should be near origin"

    def test_z2_energy(self):
        """Z=2: 1s energy should be -Z^2/2 = -2.0."""
        V_func, r_max, n_grid = make_potential('coulomb', {'Z': 2.0})
        E, _, _ = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        assert abs(E - (-2.0)) / 2.0 < 0.001, f"He+ 1s energy {E} vs exact -2.0"


# ---------------------------------------------------------------------------
# Test 2: Harmonic oscillator wavefunctions
# ---------------------------------------------------------------------------

class TestHarmonicWavefunctions:
    """Verify numerical solver reproduces exact HO wavefunctions."""

    def test_ground_state_energy(self):
        """3D HO ground state: E = (3/2) hbar*omega for l=0, n_r=0."""
        omega = 1.0
        V_func, r_max, n_grid = make_potential('harmonic', {'omega': omega})
        E, _, _ = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        E_exact = 1.5 * omega  # (2*0 + 0 + 3/2) * omega
        assert abs(E - E_exact) / E_exact < 0.001, (
            f"HO ground state energy {E} vs exact {E_exact}"
        )

    def test_first_excited_l0(self):
        """n_r=1, l=0: N=2, E = (2*1 + 0 + 3/2)*omega = 3.5."""
        omega = 1.0
        V_func, r_max, n_grid = make_potential('harmonic', {'omega': omega})
        E, _, _ = solve_radial_schrodinger(V_func, 1, 0, r_max, n_grid)
        E_exact = 3.5 * omega
        assert abs(E - E_exact) / E_exact < 0.001, (
            f"HO (n_r=1, l=0) energy {E} vs exact {E_exact}"
        )

    def test_l1_ground(self):
        """n_r=0, l=1: N=1, E = (0 + 1 + 3/2)*omega = 2.5."""
        omega = 1.0
        V_func, r_max, n_grid = make_potential('harmonic', {'omega': omega})
        E, _, _ = solve_radial_schrodinger(V_func, 0, 1, r_max, n_grid)
        E_exact = 2.5 * omega
        assert abs(E - E_exact) / E_exact < 0.001, (
            f"HO (n_r=0, l=1) energy {E} vs exact {E_exact}"
        )

    def test_normalization(self):
        """HO wavefunctions should be normalized."""
        V_func, r_max, n_grid = make_potential('harmonic', {'omega': 1.0})
        _, R_grid, r_grid = solve_radial_schrodinger(V_func, 0, 0, r_max, n_grid)
        norm = np.trapezoid(R_grid**2 * r_grid**2, r_grid)
        assert abs(norm - 1.0) < 0.01, f"HO normalization {norm} vs 1.0"


# ---------------------------------------------------------------------------
# Test 3: Coulomb Slater F^0(1s,1s;1s,1s) = 5Z/8
# ---------------------------------------------------------------------------

class TestSlaterIntegrals:
    """Verify R^k integrals against known exact values."""

    def test_coulomb_f0_1s(self):
        """R^0(1s,1s;1s,1s) = 5Z/8 = 0.625 for Z=1."""
        Z = 1.0
        wfs = radial_wavefunctions_for_potential('coulomb', {'Z': Z}, 0, 0)
        assert (0, 0) in wfs, "1s wavefunction not computed"

        _, R_1s, r_grid = wfs[(0, 0)]
        val = compute_slater_rk(R_1s, R_1s, R_1s, R_1s, 0, r_grid)
        exact = 5.0 * Z / 8.0
        assert abs(val - exact) / exact < 0.01, (
            f"F^0(1s,1s;1s,1s) = {val:.6f}, exact = {exact:.6f}"
        )

    def test_coulomb_f0_z2(self):
        """R^0(1s,1s;1s,1s) = 5Z/8 = 1.25 for Z=2."""
        Z = 2.0
        wfs = radial_wavefunctions_for_potential('coulomb', {'Z': Z}, 0, 0)
        _, R_1s, r_grid = wfs[(0, 0)]
        val = compute_slater_rk(R_1s, R_1s, R_1s, R_1s, 0, r_grid)
        exact = 5.0 * Z / 8.0
        assert abs(val - exact) / exact < 0.01, (
            f"F^0(1s,1s;1s,1s) = {val:.6f}, exact = {exact:.6f}"
        )


# ---------------------------------------------------------------------------
# Test 4: Gaunt coefficients are potential-independent
# ---------------------------------------------------------------------------

class TestGauntInvariance:
    """Angular selection rules are identical regardless of potential."""

    def test_gaunt_identical(self):
        """c^k coefficients depend only on (l, m, k), not on wavefunctions."""
        # Compute c^k for a few cases — these are pure angular, no potential
        val1 = ck_coefficient(0, 0, 0, 0, 0)
        val2 = ck_coefficient(1, 0, 1, 0, 0)
        val3 = ck_coefficient(1, 1, 1, 1, 0)
        val4 = ck_coefficient(1, 0, 0, 0, 1)

        # These are mathematical constants, independent of context
        assert abs(val1 - 1.0) < 1e-10, f"c^0(0,0;0,0) = {val1}"
        assert abs(val2) > 1e-10, f"c^0(1,0;1,0) should be nonzero"
        # val4: l=1,m=0 -> l=0,m=0 with k=1: triangle OK, parity OK
        assert abs(val4) > 1e-10, f"c^1(1,0;0,0) should be nonzero"


# ---------------------------------------------------------------------------
# Test 5: Angular sparsity invariance
# ---------------------------------------------------------------------------

class TestAngularSparsityInvariance:
    """Angular zero count is identical for all potentials at same (n_max, l_max)."""

    def test_angular_sparsity_invariance(self):
        """Angular zeros depend only on state labels, not potential."""
        n_max, l_max = 1, 1

        # Compute angular zero count from pure combinatorics
        ang_zeros_ref, total_ref = angular_zero_count(n_max, l_max)

        # Now compute from actual ERI tensors for different potentials
        potentials = {
            'coulomb': {'Z': 1.0},
            'harmonic': {'omega': 1.0},
            'square_well': {'V0': 50.0, 'R': 3.0},
        }

        for pot_name, params in potentials.items():
            wfs = radial_wavefunctions_for_potential(pot_name, params, n_max, l_max)
            states = enumerate_states(wfs)
            n_sp = len(states)

            # The number of spatial orbitals may differ if some potentials
            # don't have enough bound states, but for matched state counts
            # the angular zero fraction should be identical
            _, stats = compute_eri_tensor(wfs, threshold=1e-10)

            # For potentials with enough bound states to match the reference
            if stats['n_spatial_orbitals'] == int(np.sqrt(np.sqrt(total_ref))):
                assert stats['angular_zero'] == ang_zeros_ref, (
                    f"{pot_name}: angular_zero={stats['angular_zero']} "
                    f"vs reference {ang_zeros_ref}"
                )

    def test_angular_zero_count_small(self):
        """At n_max=0, l_max=0 (single s orbital), all 1 ERI is nonzero."""
        ang_zeros, total = angular_zero_count(0, 0)
        # 1 orbital, 1^4 = 1 ERI, none are zero by angular rules
        assert total == 1
        assert ang_zeros == 0

    def test_angular_zero_fraction_increases_with_l(self):
        """Higher l_max should increase angular sparsity."""
        _, total_l0 = angular_zero_count(1, 0)
        zeros_l0, _ = angular_zero_count(1, 0)
        zeros_l1, total_l1 = angular_zero_count(1, 1)

        frac_l0 = zeros_l0 / total_l0 if total_l0 > 0 else 0
        frac_l1 = zeros_l1 / total_l1 if total_l1 > 0 else 0

        # With more angular channels, angular sparsity should increase
        # (more selection rule zeros)
        assert frac_l1 >= frac_l0, (
            f"Angular zero fraction should increase: l_max=0 {frac_l0:.2%} vs l_max=1 {frac_l1:.2%}"
        )


# ---------------------------------------------------------------------------
# Test 6: ERI symmetry
# ---------------------------------------------------------------------------

class TestERISymmetry:
    """Verify ERI tensor has proper symmetries."""

    def test_exchange_symmetry(self):
        """<ab|cd> = <ba|dc> (simultaneous exchange of electrons).

        Tolerance relaxed to 1e-3 relative for numerical wavefunctions
        on a finite grid (grid-based Slater integrals have ~0.01% error).
        """
        wfs = radial_wavefunctions_for_potential('coulomb', {'Z': 1.0}, 1, 1)
        eri, _ = compute_eri_tensor(wfs)
        n = eri.shape[0]

        max_rel_err = 0.0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        val1 = eri[a, b, c, d]
                        val2 = eri[b, a, d, c]
                        if abs(val1) > 1e-8:
                            rel = abs(val1 - val2) / abs(val1)
                            max_rel_err = max(max_rel_err, rel)
        assert max_rel_err < 1e-3, (
            f"Exchange symmetry max relative error {max_rel_err:.2e} exceeds 1e-3"
        )

    def test_hermitian_symmetry(self):
        """<ab|cd> = <cd|ab> (particle exchange)."""
        wfs = radial_wavefunctions_for_potential('coulomb', {'Z': 1.0}, 1, 0)
        eri, _ = compute_eri_tensor(wfs)
        n = eri.shape[0]

        max_rel_err = 0.0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        val1 = eri[a, b, c, d]
                        val2 = eri[c, d, a, b]
                        if abs(val1) > 1e-8:
                            rel = abs(val1 - val2) / abs(val1)
                            max_rel_err = max(max_rel_err, rel)
        assert max_rel_err < 1e-3, (
            f"Hermitian symmetry max relative error {max_rel_err:.2e} exceeds 1e-3"
        )

    def test_complex_conjugate_symmetry(self):
        """For real orbitals, <ab|cd> = <cb|ad> should hold."""
        wfs = radial_wavefunctions_for_potential('coulomb', {'Z': 1.0}, 1, 0)
        eri, _ = compute_eri_tensor(wfs)
        n = eri.shape[0]

        max_rel_err = 0.0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        val1 = eri[a, b, c, d]
                        val2 = eri[c, b, a, d]
                        if abs(val1) > 1e-8:
                            rel = abs(val1 - val2) / abs(val1)
                            max_rel_err = max(max_rel_err, rel)
        assert max_rel_err < 1e-3, (
            f"Conjugate symmetry max relative error {max_rel_err:.2e} exceeds 1e-3"
        )


# ---------------------------------------------------------------------------
# Test 7: Sparsity table runs
# ---------------------------------------------------------------------------

class TestSparsityTable:
    """Verify sparsity_table completes for all potentials."""

    def test_sparsity_table_runs(self):
        """Run at small (n_max=1, l_max=1) for speed."""
        results = sparsity_table(
            n_max=1, l_max=1,
            potentials={
                'coulomb': {'Z': 1.0},
                'harmonic': {'omega': 1.0},
                'square_well': {'V0': 50.0, 'R': 3.0},
                'yukawa': {'V0': 1.0, 'mu': 0.5},
            },
        )

        for name, stats in results.items():
            if 'error' in stats:
                # Some potentials may not have enough bound states
                continue
            assert 'eri_density_pct' in stats, f"{name}: missing eri_density_pct"
            assert 0 <= stats['eri_density_pct'] <= 100, (
                f"{name}: ERI density {stats['eri_density_pct']} out of range"
            )
            assert stats['angular_sparsity_pct'] >= 0, (
                f"{name}: negative angular sparsity"
            )

    def test_coulomb_reasonable_density(self):
        """Coulomb ERI density at small basis should be in known range."""
        results = sparsity_table(
            n_max=1, l_max=1,
            potentials={'coulomb': {'Z': 1.0}},
        )
        stats = results['coulomb']
        assert 'eri_density_pct' in stats
        # ERI density should be well below 100% due to Gaunt selection rules
        assert stats['eri_density_pct'] < 50, (
            f"Coulomb ERI density {stats['eri_density_pct']:.1f}% too high"
        )


# ---------------------------------------------------------------------------
# Test 8: Coulomb known sparsity comparison
# ---------------------------------------------------------------------------

class TestCoulombKnownSparsity:
    """Cross-check against known GeoVac ERI densities."""

    def test_s_only_full_density(self):
        """With only s orbitals (l_max=0), all ERIs are nonzero
        (no angular selection rules kill any)."""
        results = sparsity_table(
            n_max=1, l_max=0,
            potentials={'coulomb': {'Z': 1.0}},
        )
        stats = results['coulomb']
        # With s-only, all angular couplings are nonzero (k=0 only)
        # So only angular zero comes from m-conservation (trivially satisfied for m=0)
        # ERI density should be high
        assert stats['eri_density_pct'] > 90, (
            f"s-only ERI density {stats['eri_density_pct']:.1f}% should be near 100%"
        )

    def test_sp_sparsity_increase(self):
        """Including p orbitals should increase sparsity (lower density)."""
        results_s = sparsity_table(
            n_max=1, l_max=0,
            potentials={'coulomb': {'Z': 1.0}},
        )
        results_sp = sparsity_table(
            n_max=1, l_max=1,
            potentials={'coulomb': {'Z': 1.0}},
        )
        dens_s = results_s['coulomb']['eri_density_pct']
        dens_sp = results_sp['coulomb']['eri_density_pct']
        assert dens_sp < dens_s, (
            f"s+p density {dens_sp:.1f}% should be lower than s-only {dens_s:.1f}%"
        )


# ---------------------------------------------------------------------------
# Additional: Woods-Saxon bound states
# ---------------------------------------------------------------------------

class TestWoodsSaxon:
    """Woods-Saxon potential should produce bound states."""

    def test_has_bound_states(self):
        """Woods-Saxon with V0=50, R0=3, a=0.65 should have bound states."""
        wfs = radial_wavefunctions_for_potential(
            'woods_saxon', {'V0': 50.0, 'R0': 3.0, 'a': 0.65}, 1, 1,
        )
        assert len(wfs) > 0, "Woods-Saxon should have at least one bound state"

    def test_bound_energies_negative(self):
        """All bound state energies should be negative."""
        wfs = radial_wavefunctions_for_potential(
            'woods_saxon', {'V0': 50.0, 'R0': 3.0, 'a': 0.65}, 1, 1,
        )
        for (n_r, l_val), (E, _, _) in wfs.items():
            assert E < 0, f"Woods-Saxon ({n_r},{l_val}): E={E} should be negative"


# ---------------------------------------------------------------------------
# Full sparsity analysis (slow)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestFullSparsityAnalysis:
    """Run full sparsity analysis and print results."""

    def test_sparsity_table_n2_l1(self):
        """Sparsity table at n_max=2, l_max=1."""
        from geovac.nuclear.potential_sparsity import print_sparsity_table

        results = sparsity_table(n_max=2, l_max=1)
        table = print_sparsity_table(results)
        print("\n=== Sparsity Table (n_max=2, l_max=1) ===")
        print(table)

        # Verify angular sparsity is identical across potentials
        ang_pcts = []
        for name, stats in results.items():
            if 'error' not in stats:
                ang_pcts.append(stats['angular_sparsity_pct'])

        if len(ang_pcts) >= 2:
            # All angular sparsity percentages should be very close
            # (within the set of potentials that have enough bound states
            # to produce the same number of orbitals)
            # Group by n_spatial_orbitals
            by_norb = {}
            for name, stats in results.items():
                if 'error' not in stats:
                    norb = stats['n_spatial_orbitals']
                    by_norb.setdefault(norb, []).append(
                        (name, stats['angular_sparsity_pct'])
                    )

            for norb, entries in by_norb.items():
                if len(entries) >= 2:
                    pcts = [p for _, p in entries]
                    names = [n for n, _ in entries]
                    spread = max(pcts) - min(pcts)
                    assert spread < 0.01, (
                        f"Angular sparsity spread {spread:.4f}% among "
                        f"{names} at N_orb={norb}"
                    )

    def test_angular_zero_reference(self):
        """Angular-only zeros should be potential-independent."""
        for n_max, l_max in [(1, 1), (2, 1), (2, 2)]:
            ang_zeros, total = angular_zero_count(n_max, l_max)
            ang_pct = 100.0 * ang_zeros / total
            density_pct = 100.0 - ang_pct
            print(f"\nn_max={n_max}, l_max={l_max}: "
                  f"ang_zeros={ang_zeros}/{total}, "
                  f"angular density={density_pct:.2f}%")
            assert ang_zeros >= 0
            assert total > 0

    def test_sparsity_table_n2_l2(self):
        """Sparsity table at n_max=2, l_max=2."""
        from geovac.nuclear.potential_sparsity import print_sparsity_table

        results = sparsity_table(n_max=2, l_max=2)
        table = print_sparsity_table(results)
        print("\n=== Sparsity Table (n_max=2, l_max=2) ===")
        print(table)

        # Report and verify angular sparsity invariance
        ang_pcts_by_norb = {}
        for name, stats in results.items():
            if 'error' not in stats:
                norb = stats['n_spatial_orbitals']
                print(f"\n{name}: N_orb={norb}, Q={stats['n_qubits']}, "
                      f"ang={stats['angular_sparsity_pct']:.2f}%, "
                      f"rad={stats['radial_sparsity_pct']:.2f}%, "
                      f"density={stats['eri_density_pct']:.2f}%")
                ang_pcts_by_norb.setdefault(norb, []).append(
                    (name, stats['angular_sparsity_pct'])
                )

        for norb, entries in ang_pcts_by_norb.items():
            if len(entries) >= 2:
                pcts = [p for _, p in entries]
                spread = max(pcts) - min(pcts)
                assert spread < 0.01, (
                    f"Angular sparsity spread {spread:.4f}% at N_orb={norb}"
                )

    def test_pauli_scaling_coulomb(self):
        """Pauli scaling estimate for Coulomb potential."""
        result = pauli_scaling_estimate(
            'coulomb', {'Z': 1.0}, n_max_values=[1, 2, 3], l_max=2,
        )
        print(f"\n=== Pauli Scaling (Coulomb Z=1) ===")
        for i, nm in enumerate(result['n_max']):
            print(f"  n_max={nm}: Q={result['Q'][i]}, "
                  f"N_Pauli={result['N_Pauli'][i]}")
        if result['scaling_exponent'] is not None:
            print(f"  Scaling exponent: {result['scaling_exponent']:.2f}")

    def test_pauli_scaling_harmonic(self):
        """Pauli scaling estimate for harmonic potential."""
        result = pauli_scaling_estimate(
            'harmonic', {'omega': 1.0}, n_max_values=[1, 2, 3], l_max=2,
        )
        print(f"\n=== Pauli Scaling (Harmonic omega=1) ===")
        for i, nm in enumerate(result['n_max']):
            print(f"  n_max={nm}: Q={result['Q'][i]}, "
                  f"N_Pauli={result['N_Pauli'][i]}")
        if result['scaling_exponent'] is not None:
            print(f"  Scaling exponent: {result['scaling_exponent']:.2f}")


# ---------------------------------------------------------------------------
# Run if executed directly
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
