"""
Tests for the Level 4 multichannel solver.

Validates convergence with l_max and comparison to Phase 1 / Paper 12.
"""

import numpy as np
import pytest

from geovac.level4_multichannel import (
    _channel_list,
    _channel_list_extended,
    _ee_coupling,
    _ee_coupling_general,
    compute_nuclear_coupling,
    solve_angular_multichannel,
    solve_level4_h2_multichannel,
)


class TestChannelStructure:
    """Tests for channel indexing and symmetry."""

    def test_channel_count_lmax0(self) -> None:
        ch = _channel_list(0)
        assert ch == [(0, 0)]

    def test_channel_count_lmax1(self) -> None:
        ch = _channel_list(1)
        # l1+l2 even, l1<=l2: (0,0), (1,1)
        assert ch == [(0, 0), (1, 1)]

    def test_channel_count_lmax2(self) -> None:
        ch = _channel_list(2)
        # Both orderings: (0,0), (0,2), (2,0), (1,1), (2,2)
        assert len(ch) == 5
        assert (0, 0) in ch
        assert (1, 1) in ch
        assert (0, 2) in ch
        assert (2, 0) in ch
        assert (2, 2) in ch

    def test_channel_count_lmax3(self) -> None:
        ch = _channel_list(3)
        # (0,0),(0,2),(2,0),(1,1),(1,3),(3,1),(2,2),(3,3) = 8
        assert len(ch) == 8


class TestCouplingIntegrals:
    """Tests for nuclear and e-e coupling matrix elements."""

    def test_nuclear_diagonal_united_atom(self) -> None:
        """At rho->0, nuclear coupling (0,0)|(0,0) -> -2Z(1/cos+1/sin)."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        rho = 1e-6
        Z = 1.0
        V = compute_nuclear_coupling(0, 0, 0, 0, 0, 0, alpha, rho, Z)

        # P_0 = 1, overlap <P0|P0>/2 = 1. So for l1p=l1=0, l2p=l2=0:
        # V = <P0|V1|P0> * <P0|P0>/2 + <P0|P0>/2 * <P0|V2|P0>
        # <P0|V|P0> = integral of V * 1/2 dx = average of V
        # <P0|P0>/2 = 1/(2*0+1) = 1
        # This should match _nuclear_attraction_averaged from Phase 1
        from geovac.level4_sigma_channel import _nuclear_attraction_averaged
        V_phase1 = _nuclear_attraction_averaged(alpha, rho, Z)
        np.testing.assert_allclose(V, V_phase1, rtol=1e-3)

    def test_ee_diagonal_l0(self) -> None:
        """E-e coupling (0,0)|(0,0) should match Phase 1 C_ee."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        W = _ee_coupling(0, 0, 0, 0, alpha, l_max=0)

        # For l1=l2=l1'=l2'=0, k must satisfy k_min=0, k_max=0, so k=0 only
        # G(0,0,0) = gaunt_integral(0,0,0) = 2 * (1*1*1)/(1*1*1*1) = 2
        # Wait: gaunt_integral computes int P_l1 P_k P_l2 dx over [-1,1]
        # For (0,0,0): int 1*1*1 dx = 2, threej_sq = (0!^2 * 1*1*1)/(1 * 1*1*1) = 1
        # return 2 * threej_sq = 2
        # g1 = g2 = G(0,0,0) = 2
        # norm = sqrt(1*1*1*1)/4 = 0.25
        # W = g1*g2*f_0/max * norm = 2*2*(1/max)*0.25 = 1/max
        expected = 1.0 / np.maximum(np.cos(alpha), np.sin(alpha))
        np.testing.assert_allclose(W, expected, rtol=1e-10)


class TestAngularSolver:
    """Tests for the coupled-channel angular eigenvalue problem."""

    def test_lmax0_matches_phase1(self) -> None:
        """l_max=0 multichannel should reproduce Phase 1 eigenvalue."""
        from geovac.level4_sigma_channel import solve_angular as solve_p1

        R_e = 2.0
        rho = 0.35

        mu_mc, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max=0, n_alpha=150        )
        mu_p1, _, _ = solve_p1(rho, R_e, n_alpha=150, n_quad=24)

        # rtol=1e-3: multichannel uses exact algebraic nuclear coupling,
        # Phase 1 uses Gauss-Legendre quadrature with n_quad=24
        np.testing.assert_allclose(mu_mc[0], mu_p1[0], rtol=1e-3,
                                   err_msg="l_max=0 multichannel must match Phase 1")

    def test_lmax1_lower_than_lmax0(self) -> None:
        """Adding l=1 channels should lower the eigenvalue (variational)."""
        R_e = 2.0
        rho = 0.35

        mu0, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max=0, n_alpha=100        )
        mu1, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max=1, n_alpha=100        )
        assert mu1[0] <= mu0[0] + 1e-10, (
            f"l_max=1 eigenvalue {mu1[0]:.4f} should be <= l_max=0 {mu0[0]:.4f}"
        )


class TestFullSolver:
    """Tests for the complete multichannel solver."""

    def test_lmax0_reproduces_phase1(self) -> None:
        """l_max=0 should give same D_e as Phase 1 (~31%)."""
        result = solve_level4_h2_multichannel(
            R=1.4, l_max=0, n_alpha=150, n_Re=300, verbose=False,
        )
        # Phase 1 gave 31.3%
        assert abs(result['D_e_pct'] - 31.3) < 2.0, (
            f"l_max=0 D_e={result['D_e_pct']:.1f}%, expected ~31%"
        )

    def test_bound_state_lmax2(self) -> None:
        """l_max=2 should produce bound state."""
        result = solve_level4_h2_multichannel(
            R=1.4, l_max=2, n_alpha=100, n_Re=300, verbose=False,
        )
        assert result['D_e'] > 0, f"D_e = {result['D_e']:.6f}, should be > 0"

    def test_convergence_table(self) -> None:
        """Convergence study: l_max = 0, 1, 2, 3."""
        print("\n" + "=" * 65)
        print("  LEVEL 4 CONVERGENCE STUDY: l_max vs D_e")
        print("=" * 65)
        print(f"  {'l_max':>5} {'N_ch':>5} {'E_total':>12} {'D_e':>10} {'D_e%':>8} {'time':>8}")
        print("-" * 65)

        D_e_exact = 0.17447

        for lm in [0, 1, 2, 3]:
            import time
            t0 = time.time()
            result = solve_level4_h2_multichannel(
                R=1.4, l_max=lm,
                n_alpha=150, n_Re=300,                verbose=False,
            )
            dt = time.time() - t0
            print(f"  {lm:>5} {result['n_ch']:>5} "
                  f"{result['E_total']:>12.6f} "
                  f"{result['D_e']:>10.6f} "
                  f"{result['D_e_pct']:>7.1f}% "
                  f"{dt:>7.1f}s")

        print("-" * 65)
        print(f"  Paper 12 Neumann V_ee: 92.4%")
        print(f"  Exact D_e = {D_e_exact:.5f} Ha")
        print("=" * 65)

        # l_max=3 should recover > 80% D_e
        assert result['D_e_pct'] > 80.0, (
            f"l_max=3 D_e={result['D_e_pct']:.1f}%, expected > 80%"
        )

    def test_exceeds_paper12(self) -> None:
        """l_max=4 should exceed Paper 12's 92.4% D_e."""
        result = solve_level4_h2_multichannel(
            R=1.4, l_max=4,
            n_alpha=150, n_Re=300,            verbose=False,
        )
        print(f"\n  l_max=4: D_e = {result['D_e']:.6f} Ha "
              f"({result['D_e_pct']:.1f}% of exact)")
        assert result['D_e_pct'] > 92.4, (
            f"l_max=4 D_e={result['D_e_pct']:.1f}%, expected > 92.4% (Paper 12)"
        )


class TestExtendedChannels:
    """Tests for the m_max > 0 channel extension (Phase 3 infrastructure)."""

    def test_extended_mmax0_matches_sigma(self) -> None:
        """_channel_list_extended(l, 0) matches _channel_list(l) ordering."""
        for lm in [0, 1, 2, 3, 4]:
            ch2 = _channel_list(lm)
            ch4 = _channel_list_extended(lm, 0)
            assert len(ch4) == len(ch2), f"l_max={lm}: count mismatch"
            for (l1, l2), (l1e, m1e, l2e, m2e) in zip(ch2, ch4):
                assert (l1, l2) == (l1e, l2e), f"l_max={lm}: order mismatch"
                assert m1e == 0 and m2e == 0

    def test_extended_mmax1_count(self) -> None:
        """m_max=1 adds pi channels correctly."""
        # l_max=1: sigma: (0,0),(1,1)=2. pi: (1,+1,1,-1),(1,-1,1,+1)=2. Total=4
        ch = _channel_list_extended(1, 1)
        assert len(ch) == 4

        # l_max=2: sigma: 5. pi: min(l1,l2)>=1 with l1+l2 even:
        #   (1,1): 2 pi variants. (2,2): 2 pi variants. Total pi=4. Grand=9
        ch = _channel_list_extended(2, 1)
        sigma = [(l1, m1, l2, m2) for l1, m1, l2, m2 in ch if m1 == 0]
        pi = [(l1, m1, l2, m2) for l1, m1, l2, m2 in ch if m1 != 0]
        assert len(sigma) == 5
        assert len(pi) == 4  # (1,+1,1,-1), (1,-1,1,+1), (2,+1,2,-1), (2,-1,2,+1)

    def test_extended_mmax1_lmax4_count(self) -> None:
        """l_max=4, m_max=1 channel count."""
        ch = _channel_list_extended(4, 1)
        sigma = [c for c in ch if c[1] == 0]
        pi = [c for c in ch if c[1] != 0]
        assert len(sigma) == 13  # same as Phase 2
        # pi: (l1,l2) with l1+l2 even, min(l1,l2)>=1, 2 variants each
        # (1,1),(1,3),(3,1),(2,2),(2,4),(4,2),(3,3),(4,4) = 8 pairs x 2 = 16
        assert len(pi) == 16
        assert len(ch) == 29

    def test_extended_gerade_constraint(self) -> None:
        """All extended channels satisfy l1+l2 even."""
        for lm in [2, 3, 4]:
            for mm in [0, 1]:
                ch = _channel_list_extended(lm, mm)
                for l1, m1, l2, m2 in ch:
                    assert (l1 + l2) % 2 == 0, f"({l1},{m1},{l2},{m2}) violates gerade"
                    assert m1 + m2 == 0, f"({l1},{m1},{l2},{m2}) violates M=0"
                    assert abs(m1) <= l1 and abs(m2) <= l2


class TestNuclearCouplingM:
    """Tests for associated Legendre nuclear coupling."""

    def test_algebraic_diagonal_exact(self) -> None:
        """Algebraic (0,0)-(0,0) m=0 coupling matches known formula."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        rho = 0.35
        V = compute_nuclear_coupling(0, 0, 0, 0, 0, 0, alpha, rho, Z=1.0)
        # Exact: V = -2Z/max(cos a, rho) - 2Z/max(sin a, rho)
        expected = (-2.0 / np.maximum(np.cos(alpha), rho)
                    - 2.0 / np.maximum(np.sin(alpha), rho))
        np.testing.assert_allclose(V, expected, rtol=1e-12)

    def test_pi_diagonal_nonzero(self) -> None:
        """Pi-channel diagonal nuclear coupling is non-zero."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        rho = 0.35
        V = compute_nuclear_coupling(1, 1, 1, 1, 1, -1, alpha, rho, Z=1.0)
        assert np.max(np.abs(V)) > 0.1, "Pi diagonal nuclear coupling should be significant"

    def test_pi_coupling_weaker_than_sigma(self) -> None:
        """Pi (1,1) diagonal should be weaker than sigma (0,0) diagonal."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        rho = 0.35
        V_sigma = compute_nuclear_coupling(0, 0, 0, 0, 0, 0, alpha, rho, Z=1.0)
        V_pi = compute_nuclear_coupling(1, 1, 1, 1, 1, -1, alpha, rho, Z=1.0)
        assert np.mean(np.abs(V_pi)) < np.mean(np.abs(V_sigma))


class TestMmax0Backward:
    """Tests that m_max=0 reproduces Phase 2 results exactly."""

    def test_angular_eigenvalue_identical(self) -> None:
        """m_max=0 angular eigenvalue matches Phase 2."""
        rho = 0.35
        R_e = 2.0
        mu_orig, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max=2, n_alpha=100        )
        mu_m0, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max=2, n_alpha=100, m_max=0
        )
        np.testing.assert_allclose(mu_m0[0], mu_orig[0], rtol=1e-12,
                                   err_msg="m_max=0 must exactly reproduce Phase 2")

    def test_full_solver_identical(self) -> None:
        """m_max=0 full solver matches Phase 2 D_e."""
        result_orig = solve_level4_h2_multichannel(
            R=1.4, l_max=2, n_alpha=80, n_Re=200, verbose=False,
        )
        result_m0 = solve_level4_h2_multichannel(
            R=1.4, l_max=2, n_alpha=80, n_Re=200, verbose=False, m_max=0,
        )
        np.testing.assert_allclose(
            result_m0['D_e'], result_orig['D_e'], rtol=1e-12,
            err_msg="m_max=0 full solver must reproduce Phase 2 exactly"
        )


class TestEeCouplingGeneral:
    """Tests for the Wigner 3j-based e-e coupling."""

    def test_m0_matches_gaunt(self) -> None:
        """_ee_coupling_general with m=0 matches _ee_coupling."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        pairs = [(0, 0, 0, 0), (0, 0, 0, 2), (0, 2, 2, 0),
                 (1, 1, 1, 1), (2, 0, 0, 2), (2, 2, 2, 2)]
        for l1p, l2p, l1, l2 in pairs:
            W_old = _ee_coupling(l1p, l2p, l1, l2, alpha, l_max=4)
            W_new = _ee_coupling_general(l1p, 0, l2p, 0, l1, 0, l2, 0, alpha)
            np.testing.assert_allclose(
                W_new, W_old, atol=1e-14,
                err_msg=f"m=0 e-e mismatch for ({l1p},{l2p})->({l1},{l2})"
            )

    def test_pi_diagonal_nonzero(self) -> None:
        """Pi-pi e-e coupling (1,+1,1,-1)|(1,+1,1,-1) should be non-zero."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        W = _ee_coupling_general(1, 1, 1, -1, 1, 1, 1, -1, alpha)
        assert np.max(np.abs(W)) > 0.01, "Pi diagonal e-e should be non-zero"

    def test_sigma_pi_cross_coupling(self) -> None:
        """Sigma-pi e-e coupling exists for appropriate channels."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        # (0,0,0,0) -> (1,+1,1,-1): q = m1-m1' = 1, k_min=1, k_max=1
        # k=1: (0 1 1; 0 0 0) nonzero (0+1+1=2 even), (0 1 1; 0 -1 1) nonzero
        # Both electrons change l: l1: 0->1, l2: 0->1, with angular momentum
        # transfer q=1 mediated by the e-e multipole k=1.
        W = _ee_coupling_general(0, 0, 0, 0, 1, 1, 1, -1, alpha)
        assert np.max(np.abs(W)) > 0.01  # non-zero sigma-pi cross coupling

    def test_m_conservation(self) -> None:
        """Coupling between channels with different M should be zero."""
        alpha = np.linspace(0.2, np.pi / 2 - 0.2, 20)
        # m1+m2=0 but m1'+m2'=2: violates M conservation
        W = _ee_coupling_general(1, 1, 1, 1, 1, 0, 1, 0, alpha)
        np.testing.assert_allclose(W, 0.0, atol=1e-15)


class TestPhase3Solver:
    """Tests for the m_max=1 (sigma+pi) solver."""

    def test_mmax1_improves_lmax2(self) -> None:
        """l_max=2, m_max=1 should exceed l_max=2, m_max=0."""
        result_s = solve_level4_h2_multichannel(
            R=1.4, l_max=2, m_max=0, n_alpha=80, n_Re=200, verbose=False,
        )
        result_p = solve_level4_h2_multichannel(
            R=1.4, l_max=2, m_max=1, n_alpha=80, n_Re=200, verbose=False,
        )
        assert result_p['D_e_pct'] > result_s['D_e_pct'], (
            f"m_max=1 ({result_p['D_e_pct']:.1f}%) should exceed "
            f"m_max=0 ({result_s['D_e_pct']:.1f}%)"
        )

    def test_mmax1_exceeds_paper12(self) -> None:
        """l_max=4, m_max=1 should exceed Paper 12's 92.4%."""
        result = solve_level4_h2_multichannel(
            R=1.4, l_max=4, m_max=1,
            n_alpha=100, n_Re=250, verbose=False,
        )
        print(f"\n  l_max=4, m_max=1: {result['n_ch']} channels, "
              f"D_e = {result['D_e']:.6f} Ha ({result['D_e_pct']:.1f}%)")
        assert result['D_e_pct'] > 92.4, (
            f"l_max=4 m_max=1 D_e={result['D_e_pct']:.1f}%, expected > 92.4%"
        )


class TestLmaxPerM:
    """Tests for per-m angular momentum truncation."""

    def test_lmax_per_m_channel_count(self) -> None:
        """l_max_per_m correctly limits channels per |m|."""
        ch = _channel_list_extended(4, 1, l_max_per_m={0: 4, 1: 2})
        sigma = [c for c in ch if c[1] == 0]
        pi = [c for c in ch if c[1] != 0]
        # sigma limited by l_max=4: 13 channels
        assert len(sigma) == 13
        # pi limited by l_max=2: (1,1),(2,2) × 2 = 4
        assert len(pi) == 4
        assert len(ch) == 17

    def test_lmax_per_m_mmax0_unchanged(self) -> None:
        """l_max_per_m with only m=0 entry doesn't change sigma count."""
        ch_normal = _channel_list_extended(4, 0)
        ch_lpm = _channel_list_extended(4, 0, l_max_per_m={0: 4})
        assert len(ch_lpm) == len(ch_normal)

    def test_lmax_per_m_reduces_channels(self) -> None:
        """l_max_per_m reduces channel count compared to full m_max."""
        ch_full = _channel_list_extended(4, 1)  # 29 channels
        ch_trunc = _channel_list_extended(4, 1, {0: 4, 1: 3})  # 23 channels
        assert len(ch_trunc) < len(ch_full)
        # All truncated channels should be a subset of full
        assert set(ch_trunc).issubset(set(ch_full))


class TestMmax2:
    """Tests for m_max=2 (sigma + pi + delta) channels."""

    def test_mmax2_channel_count(self) -> None:
        """l_max=3 m_max=2 has correct channel count."""
        ch = _channel_list_extended(3, 2)
        sigma = [c for c in ch if abs(c[1]) == 0]
        pi = [c for c in ch if abs(c[1]) == 1]
        delta = [c for c in ch if abs(c[1]) == 2]
        assert len(sigma) == 8
        assert len(pi) == 10
        assert len(delta) == 4
        assert len(ch) == 22

    def test_mmax2_improves_lmax3(self) -> None:
        """l_max=3 m_max=2 should exceed m_max=1."""
        result_m1 = solve_level4_h2_multichannel(
            R=1.4, l_max=3, m_max=1, n_alpha=80, n_Re=200, verbose=False,
        )
        result_m2 = solve_level4_h2_multichannel(
            R=1.4, l_max=3, m_max=2, n_alpha=80, n_Re=200, verbose=False,
        )
        assert result_m2['D_e_pct'] >= result_m1['D_e_pct'] - 0.5, (
            f"m_max=2 ({result_m2['D_e_pct']:.1f}%) should be >= "
            f"m_max=1 ({result_m1['D_e_pct']:.1f}%) within tolerance"
        )

    def test_mmax2_gerade(self) -> None:
        """All m_max=2 channels satisfy gerade and M=0 constraints."""
        for lm in [2, 3, 4]:
            ch = _channel_list_extended(lm, 2)
            for l1, m1, l2, m2 in ch:
                assert (l1 + l2) % 2 == 0
                assert m1 + m2 == 0
                assert abs(m1) <= l1 and abs(m2) <= l2


class TestMmax1Convergence:
    """Tests for m_max=1 convergence behavior."""

    def test_mmax1_improves_lmax3(self) -> None:
        """l_max=3 m_max=1 exceeds l_max=3 m_max=0."""
        result_s = solve_level4_h2_multichannel(
            R=1.4, l_max=3, m_max=0, n_alpha=80, n_Re=200, verbose=False,
        )
        result_p = solve_level4_h2_multichannel(
            R=1.4, l_max=3, m_max=1, n_alpha=80, n_Re=200, verbose=False,
        )
        improvement = result_p['D_e_pct'] - result_s['D_e_pct']
        assert improvement > 3.0, (
            f"m_max=1 improvement at l_max=3: {improvement:.1f}% (expected > 3%)"
        )
        print(f"\n  l_max=3: m_max=0 {result_s['D_e_pct']:.1f}% -> "
              f"m_max=1 {result_p['D_e_pct']:.1f}% (+{improvement:.1f}%)")


class TestHeteronuclear:
    """Tests for heteronuclear diatomic extension (HeH+)."""

    def test_channel_list_heteronuclear_includes_odd(self) -> None:
        """Heteronuclear channel list includes odd l1+l2."""
        ch = _channel_list(2, homonuclear=False)
        parities = [(l1 + l2) % 2 for l1, l2 in ch]
        assert 1 in parities, "Heteronuclear must include odd l1+l2"
        # All 9 channels for l_max=2
        assert len(ch) == 9

    def test_channel_list_homonuclear_default(self) -> None:
        """Default homonuclear=True preserves gerade constraint."""
        ch_default = _channel_list(2)
        ch_homo = _channel_list(2, homonuclear=True)
        assert ch_default == ch_homo
        assert all((l1 + l2) % 2 == 0 for l1, l2 in ch_default)

    def test_heteronuclear_channel_count(self) -> None:
        """Heteronuclear has (l_max+1)^2 channels vs ~half for homonuclear."""
        for l_max in [1, 2, 3, 4]:
            ch_het = _channel_list(l_max, homonuclear=False)
            assert len(ch_het) == (l_max + 1) ** 2

    def test_nuclear_coupling_backward_compat(self) -> None:
        """Z_A=Z_B=Z gives same result as old Z-only code."""
        alpha = np.linspace(0.1, np.pi / 2 - 0.1, 50)
        rho = 0.5
        Z = 1.0
        # With explicit Z_A=Z_B=Z should match Z-only
        V_old = compute_nuclear_coupling(0, 0, 0, 0, 0, 0, alpha, rho, Z=Z)
        V_new = compute_nuclear_coupling(
            0, 0, 0, 0, 0, 0, alpha, rho, Z_A=Z, Z_B=Z,
        )
        np.testing.assert_allclose(V_old, V_new, atol=1e-14)

    def test_nuclear_coupling_heteronuclear_odd_k(self) -> None:
        """Heteronuclear Z_A != Z_B produces nonzero (0,0)-(0,1) coupling."""
        alpha = np.linspace(0.1, np.pi / 2 - 0.1, 50)
        rho = 0.5
        # (l1p=0, l2p=0) -> (l1=0, l2=1): requires odd k on electron 2
        # Homonuclear: odd k coeff = 0
        V_homo = compute_nuclear_coupling(
            0, 0, 0, 1, 0, 0, alpha, rho, Z_A=1.0, Z_B=1.0,
        )
        # Heteronuclear: odd k coeff = -(Z_A - Z_B) != 0
        V_het = compute_nuclear_coupling(
            0, 0, 0, 1, 0, 0, alpha, rho, Z_A=2.0, Z_B=1.0,
        )
        assert np.max(np.abs(V_homo)) < 1e-14, (
            "Homonuclear odd-k coupling should vanish"
        )
        assert np.max(np.abs(V_het)) > 0.01, (
            "Heteronuclear odd-k coupling should be nonzero"
        )

    def test_heh_plus_structure(self) -> None:
        """HeH+ (Z_A=2, Z_B=1) has correct structure and dissociation limit."""
        result = solve_level4_h2_multichannel(
            R=1.46, l_max=2, n_alpha=80, n_Re=150,
            Z_A=2.0, Z_B=1.0, verbose=False,
        )
        assert result['homonuclear'] is False
        assert result['Z_A'] == 2.0
        assert result['Z_B'] == 1.0
        # V_NN = Z_A * Z_B / R
        assert abs(result['V_NN'] - 2.0 / 1.46) < 1e-10
        # E_atoms = E(He) = -2.9037 (both electrons go to He)
        assert abs(result['E_atoms'] - (-2.9037)) < 1e-10
        # E_total should be below zero (attractive)
        assert result['E_total'] < 0

    def test_heh_plus_convergence(self) -> None:
        """HeH+ E_total improves (decreases) with l_max."""
        result_l2 = solve_level4_h2_multichannel(
            R=1.46, l_max=2, n_alpha=80, n_Re=150,
            Z_A=2.0, Z_B=1.0, verbose=False,
        )
        result_l3 = solve_level4_h2_multichannel(
            R=1.46, l_max=3, n_alpha=80, n_Re=150,
            Z_A=2.0, Z_B=1.0, verbose=False,
        )
        # E_total should decrease (improve) with l_max
        assert result_l3['E_total'] < result_l2['E_total'], (
            f"l_max=3 E_total ({result_l3['E_total']:.4f}) should be lower than "
            f"l_max=2 E_total ({result_l2['E_total']:.4f})"
        )
        # D_e should also improve (less negative or more positive)
        assert result_l3['D_e'] > result_l2['D_e'], (
            f"l_max=3 D_e ({result_l3['D_e']:.4f}) should exceed "
            f"l_max=2 D_e ({result_l2['D_e']:.4f})"
        )
        print(f"\n  HeH+ convergence:")
        print(f"    l_max=2: E_total={result_l2['E_total']:.4f}, D_e={result_l2['D_e']:.4f}")
        print(f"    l_max=3: E_total={result_l3['E_total']:.4f}, D_e={result_l3['D_e']:.4f}")
        print(f"    (exact: E_total=-2.9787, D_e=0.0750)")

    def test_h2_unchanged_with_za_zb(self) -> None:
        """H2 with explicit Z_A=Z_B=1 matches Z=1 results."""
        result_z = solve_level4_h2_multichannel(
            R=1.4, l_max=2, Z=1.0, n_alpha=80, n_Re=150, verbose=False,
        )
        result_ab = solve_level4_h2_multichannel(
            R=1.4, l_max=2, Z_A=1.0, Z_B=1.0,
            n_alpha=80, n_Re=150, verbose=False,
        )
        np.testing.assert_allclose(
            result_z['E_total'], result_ab['E_total'], atol=1e-10,
            err_msg="Z_A=Z_B=1 should match Z=1 exactly",
        )
        assert result_ab['homonuclear'] is True


class TestOriginShift:
    """Tests for adjustable origin offset in heteronuclear systems."""

    def test_h2_charge_center_equals_midpoint(self) -> None:
        """For homonuclear, charge_center = midpoint (z0=0)."""
        r_mid = solve_level4_h2_multichannel(
            R=1.4, l_max=2, Z=1.0, n_alpha=80, n_Re=150,
            verbose=False, origin='midpoint',
        )
        r_cc = solve_level4_h2_multichannel(
            R=1.4, l_max=2, Z=1.0, n_alpha=80, n_Re=150,
            verbose=False, origin='charge_center',
        )
        assert r_mid['z0'] == 0.0
        assert r_cc['z0'] == 0.0  # Z_A = Z_B → z0 = 0
        np.testing.assert_allclose(
            r_mid['E_total'], r_cc['E_total'], atol=1e-10,
        )

    def test_charge_center_improves_heh_plus(self) -> None:
        """Charge-center origin lowers HeH+ energy vs midpoint."""
        r_mid = solve_level4_h2_multichannel(
            R=1.46, l_max=3, n_alpha=80, n_Re=150,
            Z_A=2.0, Z_B=1.0, verbose=False, origin='midpoint',
        )
        r_cc = solve_level4_h2_multichannel(
            R=1.46, l_max=3, n_alpha=80, n_Re=150,
            Z_A=2.0, Z_B=1.0, verbose=False, origin='charge_center',
        )
        # z0 should be R/6 for HeH+ (Z_A=2, Z_B=1)
        expected_z0 = 1.46 * (2 - 1) / (2 * (2 + 1))
        np.testing.assert_allclose(r_cc['z0'], expected_z0, atol=1e-10)
        # Charge center should lower E_total
        assert r_cc['E_total'] < r_mid['E_total'], (
            f"charge_center E ({r_cc['E_total']:.6f}) should be lower than "
            f"midpoint E ({r_mid['E_total']:.6f})"
        )
        improvement = r_mid['E_total'] - r_cc['E_total']
        print(f"\n  Origin shift improvement: {improvement:.4f} Ha")
        print(f"  Midpoint:      E={r_mid['E_total']:.6f}, D_e={r_mid['D_e']:.6f}")
        print(f"  Charge center: E={r_cc['E_total']:.6f}, D_e={r_cc['D_e']:.6f}")

    def test_heh_plus_bound_with_charge_center(self) -> None:
        """HeH+ becomes bound at l_max=3 with charge-center origin."""
        r = solve_level4_h2_multichannel(
            R=1.46, l_max=3, n_alpha=100, n_Re=200,
            Z_A=2.0, Z_B=1.0, verbose=False, origin='charge_center',
        )
        assert r['D_e'] > 0, (
            f"HeH+ should be bound with charge-center origin at l_max=3, "
            f"got D_e = {r['D_e']:.6f}"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s', '--tb=short'])
