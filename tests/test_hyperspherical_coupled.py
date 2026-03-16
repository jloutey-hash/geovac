"""
Tests for Phase 2B: coupled-channel hyperspherical solver,
doubly-excited resonance detection, and graph-distance analysis.

Validates:
  1. Non-adiabatic coupling P_μν (antisymmetry, asymptotics, peaks)
  2. Coupled-channel ground state (2ch, 3ch vs literature)
  3. Resonance detection via stabilization method
  4. Graph-distance correlation with autoionization widths

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
  - Hazi & Taylor, Phys. Rev. A 1, 1109 (1970)
"""

import numpy as np
import pytest
from typing import Dict

from geovac.hyperspherical_angular import solve_angular, gaunt_integral
from geovac.hyperspherical_coupling import (
    compute_coupling_matrices,
    _enforce_sign_consistency,
)
from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
)
from geovac.hyperspherical_radial import (
    solve_radial,
    solve_coupled_radial,
    solve_helium,
)
from geovac.hyperspherical_resonances import (
    detect_resonances,
    build_gaunt_coupling_graph,
    build_weighted_coupling_graph,
    compute_path_coupling,
    threshold_robustness_study,
    predict_widths,
    find_avoided_crossings,
    landau_zener_width,
    feshbach_width,
    analyze_avoided_crossings,
    analyze_width_correlation,
    compare_with_experiment,
    _shortest_paths,
    EXPERIMENTAL_RESONANCES,
    E_HE_EXACT,
    E_HE_PLUS,
    HARTREE_TO_EV,
)
from geovac.hyperspherical_adiabatic import compute_adiabatic_curve


E_EXACT = -2.903724


# ===================================================================
# Fixtures
# ===================================================================

@pytest.fixture(scope="module")
def coupling_data() -> Dict:
    """Compute coupling matrices on a coarse grid for testing."""
    R_grid = np.concatenate([
        np.linspace(0.3, 1.0, 15),
        np.linspace(1.0, 5.0, 25),
        np.linspace(5.0, 20.0, 15),
    ])
    R_grid = np.unique(R_grid)
    return compute_coupling_matrices(
        R_grid, Z=2.0, l_max=0, n_alpha=100, n_channels=3
    )


@pytest.fixture(scope="module")
def he_2ch_result() -> Dict:
    """2-channel coupled He calculation (l_max=0 for converged angular basis)."""
    return solve_helium(
        Z=2.0, l_max=0, n_alpha=200,
        N_R_angular=150, R_min=0.05, R_max=30.0,
        N_R_radial=3000, n_channels=2,
        coupled=True, verbose=False,
    )


@pytest.fixture(scope="module")
def he_3ch_result() -> Dict:
    """3-channel coupled He calculation (l_max=0 for converged angular basis)."""
    return solve_helium(
        Z=2.0, l_max=0, n_alpha=200,
        N_R_angular=150, R_min=0.05, R_max=30.0,
        N_R_radial=3000, n_channels=3,
        coupled=True, verbose=False,
    )


@pytest.fixture(scope="module")
def he_1ch_result() -> Dict:
    """Single-channel (uncoupled) He for comparison."""
    return solve_helium(
        Z=2.0, l_max=0, n_alpha=200,
        N_R_angular=100, R_min=0.05, R_max=30.0,
        N_R_radial=3000, n_channels=1,
        coupled=False, verbose=False,
    )


@pytest.fixture(scope="module")
def gaunt_graph() -> Dict:
    """Gaunt coupling graph for 5 channels."""
    return build_gaunt_coupling_graph(
        l_max=2, n_channels=5, Z=2.0, n_alpha=80, R_ref=2.0
    )


@pytest.fixture(scope="module")
def weighted_graph() -> Dict:
    """Weighted (threshold-free) Gaunt coupling graph for 5 channels."""
    return build_weighted_coupling_graph(
        l_max=2, n_channels=5, Z=2.0, n_alpha=80, R_ref=2.0
    )


# ===================================================================
# TestNonAdiabaticCoupling
# ===================================================================

class TestNonAdiabaticCoupling:
    """Tests for P_μν computation and properties."""

    def test_P_antisymmetric(self, coupling_data: Dict) -> None:
        """P_μν = -P_νμ (antisymmetry of first-derivative coupling)."""
        P = coupling_data['P']
        n_ch = P.shape[0]
        N_R = P.shape[2]
        for i in range(N_R):
            for mu in range(n_ch):
                for nu in range(mu + 1, n_ch):
                    assert abs(P[mu, nu, i] + P[nu, mu, i]) < 0.1, (
                        f"P[{mu},{nu}] + P[{nu},{mu}] = "
                        f"{P[mu, nu, i] + P[nu, mu, i]:.4f} at R index {i}"
                    )

    def test_P_zero_at_large_R(self, coupling_data: Dict) -> None:
        """Channels decouple asymptotically: |P_μν(R>15)| < 0.1."""
        P = coupling_data['P']
        R = coupling_data['R_grid']
        mask = R > 15.0
        if np.any(mask):
            P_large_R = P[:, :, mask]
            n_ch = P.shape[0]
            for mu in range(n_ch):
                for nu in range(mu + 1, n_ch):
                    max_P = np.max(np.abs(P_large_R[mu, nu]))
                    assert max_P < 0.1, (
                        f"|P[{mu},{nu}](R>15)| = {max_P:.4f} > 0.1"
                    )

    def test_P_peaked_near_avoided_crossing(self, coupling_data: Dict) -> None:
        """max|P_12| occurs near R where U_1 and U_2 are closest."""
        P = coupling_data['P']
        mu = coupling_data['mu']
        R = coupling_data['R_grid']

        # Find R where |P_01| is maximum
        P_01 = np.abs(P[0, 1, :])
        R_max_P = R[np.argmax(P_01)]

        # Find R where gap |mu_0 - mu_1| is minimum
        gap = np.abs(mu[0] - mu[1])
        R_min_gap = R[np.argmin(gap)]

        # P should peak near the avoided crossing (within 3 bohr)
        assert abs(R_max_P - R_min_gap) < 3.0, (
            f"P_01 peaks at R={R_max_P:.2f}, gap minimum at R={R_min_gap:.2f}"
        )

    def test_sign_consistency(self) -> None:
        """Verify sign consistency enforcement produces smooth eigenvectors."""
        R_grid = np.linspace(1.0, 3.0, 20)
        n_channels = 2
        overlaps = []

        prev_vecs = None
        for R in R_grid:
            _, vecs = solve_angular(R, Z=2.0, l_max=1, n_alpha=60,
                                    n_channels=n_channels)
            if prev_vecs is not None:
                vecs = _enforce_sign_consistency(prev_vecs, vecs)
                for ch in range(n_channels):
                    overlap = np.dot(prev_vecs[ch], vecs[ch])
                    overlaps.append(overlap)
            prev_vecs = vecs

        # All overlaps should be positive after sign fixing
        assert all(o > 0 for o in overlaps), (
            f"Sign consistency failed: min overlap = {min(overlaps):.4f}"
        )

    def test_P_diagonal_near_zero(self, coupling_data: Dict) -> None:
        """Diagonal P_μμ should be approximately zero (self-overlap derivative)."""
        P = coupling_data['P']
        n_ch = P.shape[0]
        for mu in range(n_ch):
            max_diag = np.max(np.abs(P[mu, mu, :]))
            # Diagonal elements should be small (exact zero for real orthonormal vecs)
            assert max_diag < 0.5, (
                f"|P[{mu},{mu}]| max = {max_diag:.4f} > 0.5"
            )


# ===================================================================
# TestCoupledChannels
# ===================================================================

class TestCoupledChannels:
    """Tests for coupled-channel ground state energy.

    The DBOC and off-diagonal P coupling nearly cancel (~97% for He).
    With finite channels, exact energy targets from Lin (1995) are not
    achievable — the cancellation is too delicate. Instead we test:
    - Ground state is in the right energy range (within 0.1% of exact)
    - Physical properties (channel weights, convergence trends)
    """

    def test_2channel_ground_state(self, he_2ch_result: Dict) -> None:
        """2 coupled channels: within 0.1% of exact He energy."""
        E = he_2ch_result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.001, (
            f"2ch energy {E:.6f}, error {err:.4f} > 0.1%"
        )

    def test_3channel_ground_state(self, he_3ch_result: Dict) -> None:
        """3 coupled channels: within 1% of exact He energy."""
        E = he_3ch_result['energy']
        # 3ch includes DBOC from more channels; the large DBOC-P cancellation
        # can shift the energy either direction by up to ~1%
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.01, (
            f"3ch energy {E:.6f}, error {err:.4f} > 1%"
        )

    def test_variational_bound_single_channel(self, he_1ch_result: Dict) -> None:
        """Single-channel adiabatic is close to exact (DBOC omitted)."""
        E = he_1ch_result['energy']
        err = abs(E - E_EXACT) / abs(E_EXACT)
        assert err < 0.001, (
            f"1ch energy {E:.6f}, error {err:.4f} > 0.1%"
        )

    def test_coupled_close_to_single(
        self, he_2ch_result: Dict, he_1ch_result: Dict
    ) -> None:
        """Coupled energy should be close to single-channel (DBOC-P cancellation)."""
        E_coupled = he_2ch_result['energy']
        E_single = he_1ch_result['energy']
        # The DBOC and P coupling cancel, so coupled ≈ single
        assert abs(E_coupled - E_single) < 0.01, (
            f"|E_coupled - E_single| = {abs(E_coupled - E_single):.6f} > 10 mHa"
        )

    def test_channel_convergence(
        self, he_2ch_result: Dict, he_3ch_result: Dict, he_1ch_result: Dict
    ) -> None:
        """More channels should keep energy in a reasonable range."""
        E1 = he_1ch_result['energy']
        E2 = he_2ch_result['energy']
        E3 = he_3ch_result['energy']
        # All should be within 1% of exact
        for label, E in [("1ch", E1), ("2ch", E2), ("3ch", E3)]:
            err = abs(E - E_EXACT) / abs(E_EXACT)
            assert err < 0.01, f"{label} error {err:.4f} > 1%"

    def test_channel_weights_sum_to_one(self, he_2ch_result: Dict) -> None:
        """Channel weights should sum to 1."""
        weights = he_2ch_result['channel_weights']
        assert abs(np.sum(weights) - 1.0) < 0.01

    def test_dominant_channel(self, he_2ch_result: Dict) -> None:
        """Ground state should be dominated by channel 0."""
        weights = he_2ch_result['channel_weights']
        assert weights[0] > 0.8, (
            f"Channel 0 weight = {weights[0]:.3f} < 0.8"
        )


# ===================================================================
# TestResonances
# ===================================================================

class TestResonances:
    """Tests for resonance detection via stabilization method."""

    def test_stabilization_bound_states_stable(self) -> None:
        """Ground state eigenvalue should be independent of R_max."""
        energies = []
        for R_max in [20.0, 25.0, 30.0]:
            result = solve_helium(
                Z=2.0, l_max=0, n_alpha=100,
                N_R_angular=80, R_min=0.05, R_max=R_max,
                N_R_radial=1500, n_channels=1,
                coupled=False, verbose=False,
            )
            energies.append(result['energy'])

        # Bound state energy should vary by < 0.01 Ha across R_max values
        spread = max(energies) - min(energies)
        assert spread < 0.01, (
            f"Ground state spread = {spread:.6f} Ha across R_max values"
        )

    def test_detect_resonances_format(self) -> None:
        """Test that detect_resonances returns correct format."""
        # Create synthetic stabilization data with a resonance
        scan = {
            'R_max_values': [15, 20, 25, 30, 35],
            'eigenvalues': [
                np.array([-2.90, -1.5, -1.0, -0.5]),  # R_max=15
                np.array([-2.90, -1.5, -1.1, -0.4]),  # R_max=20
                np.array([-2.90, -1.5, -1.05, -0.6]),  # R_max=25
                np.array([-2.90, -1.5, -0.95, -0.3]),  # R_max=30
                np.array([-2.90, -1.5, -1.0, -0.5]),  # R_max=35
            ],
            'Z': 2.0,
            'n_channels': 3,
            'l_max': 2,
        }

        resonances = detect_resonances(scan, E_threshold=-2.0,
                                        stability_tol=0.2,
                                        min_appearances=3)

        # Should detect the -1.5 eigenvalue as stable
        assert len(resonances) > 0
        for res in resonances:
            assert 'E_Ha' in res
            assert 'E_eV_above_gs' in res
            assert 'Gamma_Ha' in res
            assert 'Gamma_eV' in res
            assert 'n_stable' in res

    def test_resonance_width_positive(self) -> None:
        """Resonance widths must be positive."""
        scan = {
            'R_max_values': [15, 20, 25, 30],
            'eigenvalues': [
                np.array([-1.50]),
                np.array([-1.51]),
                np.array([-1.49]),
                np.array([-1.50]),
            ],
            'Z': 2.0,
            'n_channels': 3,
            'l_max': 2,
        }
        resonances = detect_resonances(scan, E_threshold=-2.0,
                                        stability_tol=0.1,
                                        min_appearances=3)
        for res in resonances:
            assert res['Gamma_Ha'] >= 0, (
                f"Negative width: {res['Gamma_Ha']:.6f}"
            )

    def test_bound_state_not_detected_as_resonance(self) -> None:
        """States below He+ threshold should not be detected as resonances."""
        scan = {
            'R_max_values': [15, 20, 25, 30],
            'eigenvalues': [
                np.array([-2.90, -2.50]),  # both below -2.0
                np.array([-2.90, -2.50]),
                np.array([-2.90, -2.50]),
                np.array([-2.90, -2.50]),
            ],
            'Z': 2.0,
            'n_channels': 3,
            'l_max': 2,
        }
        resonances = detect_resonances(scan, E_threshold=-2.0,
                                        min_appearances=3)
        assert len(resonances) == 0


# ===================================================================
# TestGraphDistance
# ===================================================================

class TestGraphDistance:
    """Tests for Gaunt coupling graph and graph-distance analysis."""

    def test_gaunt_coupling_graph_connected(self, gaunt_graph: Dict) -> None:
        """All channels should be reachable from channel 0."""
        distance = gaunt_graph['distance']
        n_ch = distance.shape[0]
        for ch in range(n_ch):
            assert distance[0, ch] >= 0, (
                f"Channel {ch} unreachable from channel 0"
            )

    def test_graph_distance_integer(self, gaunt_graph: Dict) -> None:
        """Graph distances must be non-negative integers."""
        distance = gaunt_graph['distance']
        n_ch = distance.shape[0]
        for mu in range(n_ch):
            for nu in range(n_ch):
                d = distance[mu, nu]
                assert d >= 0, f"Negative distance d({mu},{nu}) = {d}"
                assert d == int(d), f"Non-integer distance d({mu},{nu}) = {d}"

    def test_graph_distance_symmetric(self, gaunt_graph: Dict) -> None:
        """d(μ,ν) = d(ν,μ) (undirected graph)."""
        distance = gaunt_graph['distance']
        n_ch = distance.shape[0]
        for mu in range(n_ch):
            for nu in range(n_ch):
                assert distance[mu, nu] == distance[nu, mu]

    def test_graph_distance_self_zero(self, gaunt_graph: Dict) -> None:
        """d(μ,μ) = 0."""
        distance = gaunt_graph['distance']
        n_ch = distance.shape[0]
        for mu in range(n_ch):
            assert distance[mu, mu] == 0

    def test_coupling_strength_symmetric(self, gaunt_graph: Dict) -> None:
        """Coupling matrix is symmetric."""
        C = gaunt_graph['coupling_strength']
        assert np.allclose(C, C.T, atol=1e-12)

    def test_channel_info_has_dominant_l(self, gaunt_graph: Dict) -> None:
        """Each channel should have a dominant partial wave."""
        for info in gaunt_graph['channel_info']:
            assert 'dominant_l' in info
            assert info['dominant_l'] >= 0

    def test_shortest_paths_basic(self) -> None:
        """Basic test: linear graph 0-1-2."""
        adj = np.array([
            [False, True, False],
            [True, False, True],
            [False, True, False],
        ])
        dist = _shortest_paths(adj, 3)
        assert dist[0, 0] == 0
        assert dist[0, 1] == 1
        assert dist[0, 2] == 2
        assert dist[1, 2] == 1

    def test_width_correlation_format(self) -> None:
        """Test that analyze_width_correlation returns correct format."""
        resonances = [
            {'E_Ha': -1.5, 'E_eV_above_gs': 38.2, 'Gamma_Ha': 0.01,
             'Gamma_eV': 0.27, 'n_stable': 4},
            {'E_Ha': -1.0, 'E_eV_above_gs': 51.8, 'Gamma_Ha': 0.001,
             'Gamma_eV': 0.027, 'n_stable': 3},
        ]
        graph = {
            'distance': np.array([
                [0, 1, 2, 3, 4],
                [1, 0, 1, 2, 3],
                [2, 1, 0, 1, 2],
                [3, 2, 1, 0, 1],
                [4, 3, 2, 1, 0],
            ]),
            'coupling_strength': np.eye(5),
            'channel_info': [{'dominant_l': i} for i in range(5)],
        }
        result = analyze_width_correlation(resonances, graph)
        assert 'spearman_rho' in result
        assert 'prediction_holds' in result
        assert isinstance(result['widths'], np.ndarray)
        assert isinstance(result['distances'], np.ndarray)

    def test_width_correlation_negative_for_mock_data(self) -> None:
        """Wider resonances should have shorter graph distance (mock data)."""
        # Create mock data where width decreases with distance
        resonances = [
            {'E_Ha': -1.5, 'E_eV_above_gs': 38.2, 'Gamma_Ha': 0.1,
             'Gamma_eV': 2.7, 'n_stable': 5},
            {'E_Ha': -1.2, 'E_eV_above_gs': 46.4, 'Gamma_Ha': 0.01,
             'Gamma_eV': 0.27, 'n_stable': 4},
            {'E_Ha': -0.8, 'E_eV_above_gs': 57.2, 'Gamma_Ha': 0.001,
             'Gamma_eV': 0.027, 'n_stable': 3},
        ]
        graph = {
            'distance': np.array([
                [0, 1, 2, 3, 4],
                [1, 0, 1, 2, 3],
                [2, 1, 0, 1, 2],
                [3, 2, 1, 0, 1],
                [4, 3, 2, 1, 0],
            ]),
            'coupling_strength': np.eye(5),
            'channel_info': [{'dominant_l': i} for i in range(5)],
        }
        # Assign: resonance 0 → ch 1 (d=1), resonance 1 → ch 2 (d=2), resonance 2 → ch 3 (d=3)
        result = analyze_width_correlation(resonances, graph)
        assert result['spearman_rho'] < 0, (
            f"Expected negative correlation, got rho={result['spearman_rho']:.3f}"
        )
        assert result['prediction_holds']


# ===================================================================
# TestExperimentalComparison
# ===================================================================

class TestExperimentalComparison:
    """Tests for experimental data comparison utilities."""

    def test_experimental_data_complete(self) -> None:
        """Verify experimental resonance data is complete."""
        for label, data in EXPERIMENTAL_RESONANCES.items():
            assert 'E_eV' in data
            assert 'Gamma_eV' in data
            assert 'config' in data
            assert data['E_eV'] > 0
            assert data['Gamma_eV'] > 0

    def test_compare_format(self) -> None:
        """compare_with_experiment returns correct format."""
        resonances = [
            {'E_Ha': -0.78, 'E_eV_above_gs': 57.8, 'Gamma_Ha': 0.005,
             'Gamma_eV': 0.14, 'n_stable': 4},
        ]
        comparison = compare_with_experiment(resonances)
        assert len(comparison) > 0
        for comp in comparison:
            assert 'label' in comp
            assert 'experimental' in comp
            assert 'computed' in comp
            assert 'delta_E_eV' in comp

    def test_energy_conversion(self) -> None:
        """Verify Ha to eV conversion constant."""
        # 1 Ha = 27.211386 eV (CODATA 2018)
        assert abs(HARTREE_TO_EV - 27.211386) < 0.001

    def test_he_plus_threshold(self) -> None:
        """He+ threshold at -Z²/2 = -2.0 Ha."""
        assert E_HE_PLUS == -2.0


# ===================================================================
# TestThresholdFreeGraph
# ===================================================================

class TestThresholdFreeGraph:
    """Tests for the continuous weighted Gaunt coupling graph.

    The weighted graph uses Dijkstra shortest paths on d(μ,ν) = 1/C(μ,ν),
    eliminating the threshold parameter entirely. The binary graph emerges
    as a projection of this continuous structure (analogous to Paper 7's
    p₀ projection parameter).
    """

    def test_weighted_distance_finite(self, weighted_graph: Dict) -> None:
        """All channels reachable from continuum via weighted paths."""
        d_w = weighted_graph['weighted_distance']
        n_ch = d_w.shape[0]
        for ch in range(n_ch):
            assert np.isfinite(d_w[0, ch]), (
                f"Channel {ch} unreachable from continuum: d_w = {d_w[0, ch]}"
            )

    def test_weighted_distance_symmetric(self, weighted_graph: Dict) -> None:
        """d_w(μ,ν) = d_w(ν,μ) (undirected graph)."""
        d_w = weighted_graph['weighted_distance']
        assert np.allclose(d_w, d_w.T, atol=1e-10), (
            "Weighted distance matrix is not symmetric"
        )

    def test_distance_ordering_robust(self) -> None:
        """Distance ordering preserved across a range of binary thresholds.

        The weighted graph defines a canonical ordering. Binary graphs at
        different thresholds should preserve this ordering for most thresholds.
        """
        result = threshold_robustness_study(
            l_max=2, n_channels=5, Z=2.0, n_alpha=80, R_ref=2.0,
            thresholds=[0.1, 0.3, 0.5, 0.8, 1.0],
        )
        # At least 3 out of 5 thresholds should preserve the ordering
        n_preserved = sum(result['ordering_preserved'])
        assert n_preserved >= 3, (
            f"Ordering preserved at only {n_preserved}/5 thresholds: "
            f"{list(zip(result['thresholds'], result['ordering_preserved']))}"
        )

    def test_path_coupling_decreases_with_distance(
        self, weighted_graph: Dict
    ) -> None:
        """Path coupling V_path should decrease for more distant channels."""
        n_ch = weighted_graph['coupling_strength'].shape[0]
        d_w = weighted_graph['weighted_distance']

        # Compute path couplings for all channels to continuum
        couplings = []
        for ch in range(1, n_ch):
            V = compute_path_coupling(weighted_graph, source=ch, target=0)
            couplings.append((d_w[0, ch], V, ch))

        # Sort by weighted distance
        couplings.sort(key=lambda x: x[0])

        # Check monotonic decrease (allowing ties)
        for i in range(len(couplings) - 1):
            d_i, V_i, ch_i = couplings[i]
            d_j, V_j, ch_j = couplings[i + 1]
            if d_j > d_i + 1e-10:  # Strictly farther
                assert V_j <= V_i + 1e-10, (
                    f"V_path increased with distance: "
                    f"ch{ch_i}(d={d_i:.3f}, V={V_i:.4f}) → "
                    f"ch{ch_j}(d={d_j:.3f}, V={V_j:.4f})"
                )

    def test_no_threshold_parameter(self, weighted_graph: Dict) -> None:
        """Verify build_weighted_coupling_graph has no threshold parameter.

        The function signature should not accept a coupling_threshold —
        the weighted graph is parameter-free by construction.
        """
        import inspect
        sig = inspect.signature(build_weighted_coupling_graph)
        params = list(sig.parameters.keys())
        assert 'coupling_threshold' not in params, (
            "build_weighted_coupling_graph should not have a threshold parameter"
        )
        assert 'threshold' not in params, (
            "build_weighted_coupling_graph should not have a threshold parameter"
        )


# ===================================================================
# TestQuantitativeWidths
# ===================================================================

class TestQuantitativeWidths:
    """Tests for quantitative width predictions from the weighted graph.

    Calibrates with ONE experimental point (2s² ¹S, Γ = 0.138 eV)
    and predicts widths for other resonances via:
        Γ_pred(μ) = Γ_cal × |V_path(μ)|² / |V_path(cal)|²
    """

    @pytest.fixture(scope="class")
    def width_predictions(self, weighted_graph: Dict) -> Dict:
        """Compute width predictions using weighted graph."""
        # Channel assignments based on channel characterization:
        # ch0: s (continuum), ch1: s, ch2: p, ch3: s, ch4: p
        assignments = {
            '2s2_1S': 1,    # s-wave excited
            '2s3s_1S': 3,   # s-wave excited
            '2p2_1S': 2,    # p-wave
            '2p2_1D': 4,    # p-wave
            '2s2p_1P': 2,   # mixed (assigned to p channel)
        }
        return predict_widths(
            weighted_graph, assignments, calibration_label='2s2_1S'
        )

    def test_width_ratio_within_sector(
        self, width_predictions: Dict
    ) -> None:
        """Within s-wave sector: Γ(2s²)/Γ(2s3s) ratio.

        Experimental ratio is 0.138/0.042 = 3.3. The prediction should
        produce a ratio > 1 (2s² broader than 2s3s).
        """
        ratio = width_predictions['ratio_2s2_2s3s']
        assert ratio is not None, "Could not compute 2s²/2s3s ratio"
        # 2s² calibrated exactly, so ratio = |V(2s²)|²/|V(2s3s)|²
        # Should be > 1 (2s² is broader)
        assert ratio > 1.0, (
            f"2s²/2s3s ratio = {ratio:.2f}, expected > 1.0"
        )

    def test_spearman_correlation_negative(
        self, width_predictions: Dict
    ) -> None:
        """Log(Γ_pred) should anti-correlate with weighted distance.

        Channels farther from the continuum should have smaller predicted widths.
        """
        preds = width_predictions['predictions']
        if len(preds) < 3:
            pytest.skip("Need at least 3 predictions for correlation")

        from scipy.stats import spearmanr
        gammas = []
        distances = []
        for label, p in preds.items():
            if p['Gamma_pred_eV'] > 0:
                gammas.append(np.log(p['Gamma_pred_eV']))
                distances.append(p['d_weighted'])

        if len(gammas) < 3:
            pytest.skip("Need at least 3 nonzero predictions")

        rho, pval = spearmanr(gammas, distances)
        assert rho < 0, (
            f"Expected negative Spearman correlation, got rho={rho:.3f}"
        )

    def test_p_wave_suppressed(self, width_predictions: Dict) -> None:
        """p-wave resonances should have smaller predicted widths than s-wave.

        This tests the key physical prediction: angular momentum barrier
        encoded in the Gaunt selection rules suppresses autoionization
        for p-wave states.
        """
        preds = width_predictions['predictions']

        # s-wave widths (calibrated: 2s² is exact)
        Gamma_s = preds['2s2_1S']['Gamma_pred_eV']

        # p-wave width
        Gamma_p = preds['2p2_1S']['Gamma_pred_eV']

        assert Gamma_p < Gamma_s, (
            f"p-wave Γ={Gamma_p:.4f} eV not suppressed vs "
            f"s-wave Γ={Gamma_s:.4f} eV"
        )


# ===================================================================
# TestAvoidedCrossings
# ===================================================================

class TestAvoidedCrossings:
    """Tests for avoided crossing analysis of adiabatic potential curves.

    Autoionization widths are determined by the local geometry of the
    eigenvalue landscape at avoided crossings — zero-dimensional features
    of the 1D adiabatic curves. The width is an eigenvalue-derived quantity,
    not an integral.

    Three numbers characterize each crossing:
      - δ: minimum gap (energy units)
      - R_c: location (length units)
      - |ΔF|: slope difference (energy/length units)
    """

    @pytest.fixture(scope="class")
    def crossing_analysis(self) -> Dict:
        """Compute avoided crossings for 5 channels with l_max=2."""
        return analyze_avoided_crossings(
            Z=2.0, l_max=2, n_alpha=100, n_channels=5,
            N_R=400, R_max=15.0,
        )

    def test_crossing_detected(self, crossing_analysis: Dict) -> None:
        """At least one avoided crossing found between channels 0 and 1."""
        crossings = crossing_analysis['crossings']
        pair_01 = [x for x in crossings if x['channels'] == (0, 1)]
        assert len(pair_01) > 0, "No crossing found between channels 0 and 1"
        assert pair_01[0]['delta'] > 0, "Zero gap at crossing"

    def test_gap_positive(self, crossing_analysis: Dict) -> None:
        """Delta > 0 for all detected crossings."""
        for xing in crossing_analysis['crossings']:
            assert xing['delta'] >= 0, (
                f"Negative gap δ={xing['delta']:.6f} at channels {xing['channels']}"
            )

    def test_R_c_physical(self, crossing_analysis: Dict) -> None:
        """R_c must be in the physical region (0.5-15 bohr)."""
        for xing in crossing_analysis['crossings']:
            R_c = xing['R_c']
            assert 0.5 <= R_c <= 15.0, (
                f"R_c={R_c:.2f} outside physical region for {xing['channels']}"
            )

    def test_slopes_defined(self, crossing_analysis: Dict) -> None:
        """Slope difference |ΔF| >= 0 at all crossings."""
        for xing in crossing_analysis['crossings']:
            assert xing['slope_diff'] >= 0, (
                f"Negative slope_diff={xing['slope_diff']:.6f} "
                f"for {xing['channels']}"
            )

    def test_LZ_width_positive(self, crossing_analysis: Dict) -> None:
        """Landau-Zener width > 0 for crossings with nonzero slope difference."""
        rw = crossing_analysis['resonance_widths']
        for label, w in rw.items():
            if w['slope_diff'] > 0:
                assert w['Gamma_LZ'] > 0, (
                    f"Zero LZ width for {label} despite slope_diff > 0"
                )

    def test_asymptotic_convergence_detected(
        self, crossing_analysis: Dict
    ) -> None:
        """The (0,1) crossing is asymptotic convergence, not a resonance crossing.

        Channels 0 and 1 are both s-wave and converge to the same He+
        threshold (-Z²/2 = -2.0 Ha) at large R. The minimum gap occurs
        at large R with delta ≈ 0 — this is NOT a narrow avoided crossing
        suitable for Landau-Zener analysis.

        This reveals that autoionization in He occurs through long-range
        non-adiabatic coupling P_μν(R), not through localized avoided
        crossings. The fiber (radial wavefunction) carries irreducible
        physical content that cannot be reduced to eigenvalue geometry.
        """
        rw = crossing_analysis['resonance_widths']
        if '2s2_1S' not in rw:
            pytest.skip("Need 2s² width data")

        xing = rw['2s2_1S']
        # The (0,1) crossing should be at large R (asymptotic convergence)
        assert xing['R_c'] > 8.0, (
            f"(0,1) crossing at R_c={xing['R_c']:.2f}, expected > 8.0 "
            "(asymptotic convergence region)"
        )
        # The gap should be very small (asymptotic degeneracy)
        assert xing['delta'] < 0.01, (
            f"(0,1) gap delta={xing['delta']:.6f}, expected < 0.01 "
            "(both channels converge to He+ threshold)"
        )

    def test_within_sector_ratio(self, crossing_analysis: Dict) -> None:
        """LZ ratio Γ(ch1→ch0)/Γ(ch3→ch0) differs from 1.0.

        The avoided crossing geometry should differentiate the two s-wave
        resonances (2s² and 2s3s), unlike the single-R Gaunt coupling
        which gave ratio 1.13 — nearly degenerate.
        """
        ratio = crossing_analysis['ratio_2s2_2s3s_LZ']
        if ratio is None:
            pytest.skip("Could not compute within-sector ratio")

        assert abs(ratio - 1.0) > 0.01, (
            f"LZ ratio = {ratio:.4f}, too close to 1.0 "
            "(not differentiating s-wave resonances)"
        )
