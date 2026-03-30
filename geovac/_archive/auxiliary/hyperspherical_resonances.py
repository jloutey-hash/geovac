"""
Doubly-excited resonance detection and graph-distance analysis.

Uses the stabilization method (Hazi & Taylor 1970) to identify
autoionizing resonances of helium from the coupled-channel
hyperspherical solver.

Tests the hypothesis that autoionization widths correlate with
graph distance on the hyperangular Gaunt coupling graph.

References:
  - Hazi & Taylor, Phys. Rev. A 1, 1109 (1970) — stabilization method
  - Fano, Phys. Rev. 124, 1866 (1961) — autoionization theory
  - Lin, Phys. Rep. 257, 1 (1995) — hyperspherical classification
  - Madden & Codling, Phys. Rev. Lett. 10, 516 (1963) — experimental data
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import shortest_path as dijkstra_shortest_path
from scipy.sparse import csr_matrix
from typing import Tuple, List, Dict, Optional

from geovac.hyperspherical_angular import solve_angular, gaunt_integral
from geovac.hyperspherical_coupling import compute_coupling_matrices
from geovac.hyperspherical_adiabatic import compute_adiabatic_curve, effective_potential
from geovac.hyperspherical_radial import solve_coupled_radial


# He ground state energy and ionization thresholds
E_HE_EXACT = -2.903724  # Ha
E_HE_PLUS = -2.0  # Ha (He+ threshold, 1s)
HARTREE_TO_EV = 27.211386  # eV per Ha

# Experimental doubly-excited resonances (eV above He ground state)
EXPERIMENTAL_RESONANCES = {
    '2s2_1S': {'E_eV': 57.82, 'Gamma_eV': 0.138, 'config': (0, 0)},
    '2s2p_1P': {'E_eV': 60.15, 'Gamma_eV': 0.037, 'config': (0, 1)},
    '2p2_1S': {'E_eV': 62.06, 'Gamma_eV': 0.0014, 'config': (1, 1)},
    '2s3s_1S': {'E_eV': 62.94, 'Gamma_eV': 0.042, 'config': (0, 0)},
    '2p2_1D': {'E_eV': 59.90, 'Gamma_eV': 0.070, 'config': (1, 1)},
}


def stabilization_scan(
    Z: float = 2.0,
    l_max: int = 2,
    n_alpha: int = 100,
    n_channels: int = 5,
    N_R_angular: int = 200,
    N_R_radial: int = 3000,
    R_max_values: Optional[List[float]] = None,
    n_states: int = 25,
    E_min: float = -3.5,
    E_max: float = 3.0,
    verbose: bool = True,
) -> Dict:
    """
    Perform a stabilization scan: solve coupled-channel equation at
    several R_max values and track eigenvalue stability.

    Bound states are stable (independent of R_max). Resonances appear
    as approximately stationary eigenvalues above the He+ threshold.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_alpha : int
        Angular grid points.
    n_channels : int
        Number of adiabatic channels.
    N_R_angular : int
        R grid points for angular solve.
    N_R_radial : int
        R grid points for radial solve.
    R_max_values : list of float, optional
        R_max values for the scan. Default: [15, 20, 25, 30, 35, 40].
    n_states : int
        Number of eigenvalues to find at each R_max.
    E_min, E_max : float
        Energy window (Ha) for eigenvalue search.
    verbose : bool
        Print progress.

    Returns
    -------
    result : dict
        Keys:
        - 'R_max_values': list of R_max used
        - 'eigenvalues': list of ndarrays, eigenvalues at each R_max
        - 'channel_data': list of coupling data dicts
        - 'Z', 'n_channels', 'l_max'
    """
    if R_max_values is None:
        R_max_values = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0]

    all_eigenvalues = []

    for R_max in R_max_values:
        if verbose:
            print(f"\n--- Stabilization scan: R_max = {R_max:.1f} bohr ---")

        # Compute coupling on adaptive grid
        R_grid_ang = np.concatenate([
            np.linspace(0.1, 1.0, N_R_angular // 3),
            np.linspace(1.0, 5.0, N_R_angular // 3),
            np.linspace(5.0, R_max, N_R_angular // 3 + 1),
        ])
        R_grid_ang = np.unique(R_grid_ang)

        coupling = compute_coupling_matrices(
            R_grid_ang, Z, l_max, n_alpha, n_channels
        )
        mu_arr = coupling['mu']
        P = coupling['P']

        # Build splines
        V_eff_splines = []
        for ch in range(n_channels):
            V_eff_ch = effective_potential(R_grid_ang, mu_arr[ch])
            V_eff_splines.append(
                CubicSpline(R_grid_ang, V_eff_ch, extrapolate=True)
            )

        P_splines = [[None] * n_channels for _ in range(n_channels)]
        for mu_idx in range(n_channels):
            for nu_idx in range(n_channels):
                P_splines[mu_idx][nu_idx] = CubicSpline(
                    R_grid_ang, P[mu_idx, nu_idx], extrapolate=True
                )

        # Solve coupled-channel equation
        sigma = (E_min + E_max) / 2.0
        E_all, F_all, R_grid_rad = solve_coupled_radial(
            V_eff_splines, P_splines,
            n_channels=n_channels,
            R_min=0.05, R_max=R_max,
            N_R=N_R_radial,
            n_states=n_states,
            sigma=sigma,
        )

        # Filter to energy window
        mask = (E_all >= E_min) & (E_all <= E_max)
        E_window = E_all[mask]

        all_eigenvalues.append(E_window)

        if verbose:
            n_bound = np.sum(E_window < E_HE_PLUS)
            n_above = np.sum(E_window >= E_HE_PLUS)
            print(f"  Found {len(E_window)} eigenvalues "
                  f"({n_bound} bound, {n_above} above threshold)")
            if len(E_window) > 0:
                print(f"  E range: [{E_window[0]:.4f}, {E_window[-1]:.4f}] Ha")

    return {
        'R_max_values': R_max_values,
        'eigenvalues': all_eigenvalues,
        'Z': Z,
        'n_channels': n_channels,
        'l_max': l_max,
    }


def detect_resonances(
    scan_result: Dict,
    E_threshold: float = -2.0,
    stability_tol: float = 0.01,
    min_appearances: int = 3,
) -> List[Dict]:
    """
    Detect resonances from a stabilization scan by finding eigenvalues
    that are approximately stationary across R_max values.

    Parameters
    ----------
    scan_result : dict
        Output from stabilization_scan.
    E_threshold : float
        Energy threshold (Ha). Only consider E > E_threshold (above He+).
    stability_tol : float
        Maximum variation (Ha) for an eigenvalue to be considered stable.
    min_appearances : int
        Minimum number of R_max values where the eigenvalue must appear.

    Returns
    -------
    resonances : list of dict
        Each dict has keys:
        - 'E_Ha': resonance energy (Ha)
        - 'E_eV_above_gs': energy above He ground state (eV)
        - 'Gamma_Ha': estimated width (Ha)
        - 'Gamma_eV': estimated width (eV)
        - 'n_stable': number of R_max values where eigenvalue is stable
    """
    eigenvalues = scan_result['eigenvalues']
    R_max_values = scan_result['R_max_values']

    # Collect all eigenvalues above threshold
    all_above = []
    for i, E_arr in enumerate(eigenvalues):
        for E in E_arr:
            if E > E_threshold:
                all_above.append((E, i))

    if not all_above:
        return []

    # Cluster eigenvalues that appear at similar energies across R_max
    all_above.sort(key=lambda x: x[0])
    resonances = []
    used = set()

    for i, (E_ref, rmax_idx) in enumerate(all_above):
        if i in used:
            continue

        # Find all eigenvalues within stability_tol of E_ref
        cluster = [(E_ref, rmax_idx)]
        cluster_indices = {i}
        rmax_seen = {rmax_idx}

        for j, (E_j, rmax_j) in enumerate(all_above):
            if j in used or j in cluster_indices:
                continue
            if rmax_j in rmax_seen:
                continue
            if abs(E_j - E_ref) < stability_tol:
                cluster.append((E_j, rmax_j))
                cluster_indices.add(j)
                rmax_seen.add(rmax_j)

        if len(cluster) >= min_appearances:
            Es = np.array([c[0] for c in cluster])
            E_mean = np.mean(Es)
            E_spread = np.max(Es) - np.min(Es)

            resonances.append({
                'E_Ha': E_mean,
                'E_eV_above_gs': (E_mean - E_HE_EXACT) * HARTREE_TO_EV,
                'Gamma_Ha': E_spread,
                'Gamma_eV': E_spread * HARTREE_TO_EV,
                'n_stable': len(cluster),
            })
            used |= cluster_indices

    # Sort by energy
    resonances.sort(key=lambda r: r['E_Ha'])
    return resonances


def build_gaunt_coupling_graph(
    l_max: int = 2,
    n_channels: int = 5,
    Z: float = 2.0,
    n_alpha: int = 100,
    R_ref: float = 2.0,
    coupling_threshold: float = 0.01,
) -> Dict:
    """
    Build the Gaunt coupling graph on the adiabatic channels.

    Channels are nodes. An edge connects channels μ and ν if their
    Gaunt coupling |G_μν| exceeds a threshold. The graph distance
    d(μ, ν) is the shortest path length.

    Parameters
    ----------
    l_max : int
        Maximum partial wave.
    n_channels : int
        Number of channels.
    Z : float
        Nuclear charge.
    n_alpha : int
        Angular grid points.
    R_ref : float
        Reference hyperradius for computing coupling strengths.
    coupling_threshold : float
        Minimum |coupling| to create an edge.

    Returns
    -------
    result : dict
        Keys:
        - 'adjacency': ndarray (n_channels, n_channels) boolean
        - 'coupling_strength': ndarray (n_channels, n_channels)
        - 'distance': ndarray (n_channels, n_channels) int (shortest path)
        - 'channel_info': list of dicts with dominant l for each channel
    """
    # Solve angular at R_ref to get eigenvectors
    mu, vecs = solve_angular(R_ref, Z, l_max, n_alpha, n_channels)

    n_l = l_max + 1

    # Compute coupling strength between channels via Gaunt integrals
    # Each eigenvector is in the (l, alpha) basis. The channel's l-composition
    # tells us which partial waves dominate.
    coupling = np.zeros((n_channels, n_channels))
    channel_info = []

    for ch in range(n_channels):
        # Compute l-composition of each channel
        l_weights = np.zeros(n_l)
        for l in range(n_l):
            l_weights[l] = np.sum(vecs[ch, l * n_alpha:(l + 1) * n_alpha]**2)
        dominant_l = np.argmax(l_weights)
        channel_info.append({
            'channel': ch,
            'dominant_l': int(dominant_l),
            'l_weights': l_weights,
            'mu': mu[ch],
        })

    # Coupling between channels via V_ee Gaunt integrals
    # Use the overlap of channel eigenvectors mediated by V_ee
    for mu_idx in range(n_channels):
        for nu_idx in range(mu_idx + 1, n_channels):
            # Compute direct coupling via Gaunt-weighted overlap
            l_mu = channel_info[mu_idx]['dominant_l']
            l_nu = channel_info[nu_idx]['dominant_l']

            # Sum over all Gaunt couplings weighted by l-composition
            strength = 0.0
            for l1 in range(n_l):
                w1 = channel_info[mu_idx]['l_weights'][l1]
                for l2 in range(n_l):
                    w2 = channel_info[nu_idx]['l_weights'][l2]
                    for k in range(abs(l1 - l2), l1 + l2 + 1):
                        g = gaunt_integral(l1, k, l2)
                        strength += w1 * w2 * abs(g)

            coupling[mu_idx, nu_idx] = strength
            coupling[nu_idx, mu_idx] = strength

    # Build adjacency matrix
    adjacency = coupling > coupling_threshold

    # Compute shortest path distances (BFS)
    distance = _shortest_paths(adjacency, n_channels)

    return {
        'adjacency': adjacency,
        'coupling_strength': coupling,
        'distance': distance,
        'channel_info': channel_info,
    }


def _shortest_paths(adjacency: np.ndarray, n: int) -> np.ndarray:
    """Compute all-pairs shortest paths on an undirected graph via BFS."""
    dist = np.full((n, n), -1, dtype=int)
    for start in range(n):
        dist[start, start] = 0
        queue = [start]
        head = 0
        while head < len(queue):
            u = queue[head]
            head += 1
            for v in range(n):
                if adjacency[u, v] and dist[start, v] == -1:
                    dist[start, v] = dist[start, u] + 1
                    queue.append(v)
    return dist


def analyze_width_correlation(
    resonances: List[Dict],
    graph_result: Dict,
    channel_assignments: Optional[Dict[int, int]] = None,
) -> Dict:
    """
    Analyze correlation between resonance widths and graph distance.

    Parameters
    ----------
    resonances : list of dict
        Detected resonances from detect_resonances.
    graph_result : dict
        Output from build_gaunt_coupling_graph.
    channel_assignments : dict, optional
        Map from resonance index to dominant channel index.
        If None, assigns based on energy ordering.

    Returns
    -------
    result : dict
        Keys:
        - 'widths': ndarray of resonance widths
        - 'distances': ndarray of graph distances to continuum (channel 0)
        - 'spearman_rho': Spearman rank correlation coefficient
        - 'spearman_p': p-value for correlation
        - 'prediction_holds': bool, True if negative correlation found
    """
    if len(resonances) < 2:
        return {
            'widths': np.array([]),
            'distances': np.array([]),
            'spearman_rho': 0.0,
            'spearman_p': 1.0,
            'prediction_holds': False,
        }

    distance = graph_result['distance']
    n_ch = distance.shape[0]

    widths = []
    distances = []

    for i, res in enumerate(resonances):
        if channel_assignments and i in channel_assignments:
            ch = channel_assignments[i]
        else:
            # Default: assign resonances to channels 1, 2, ... in order
            ch = min(i + 1, n_ch - 1)

        widths.append(res['Gamma_eV'])
        # Distance to channel 0 (lowest channel, connects to continuum)
        d = distance[ch, 0] if distance[ch, 0] >= 0 else n_ch
        distances.append(d)

    widths = np.array(widths)
    distances = np.array(distances)

    # Spearman rank correlation between log(Gamma) and d_graph
    if len(widths) >= 2 and np.all(widths > 0):
        from scipy.stats import spearmanr
        log_widths = np.log(widths)
        rho, p_value = spearmanr(log_widths, distances)
    else:
        rho, p_value = 0.0, 1.0

    return {
        'widths': widths,
        'distances': distances,
        'spearman_rho': rho,
        'spearman_p': p_value,
        'prediction_holds': rho < 0,
    }


def compare_with_experiment(
    resonances: List[Dict],
) -> List[Dict]:
    """
    Compare detected resonances with experimental data.

    Parameters
    ----------
    resonances : list of dict
        Detected resonances.

    Returns
    -------
    comparison : list of dict
        Each dict has keys: 'computed', 'experimental', 'delta_E_eV',
        'label'.
    """
    comparison = []

    for label, expt in EXPERIMENTAL_RESONANCES.items():
        E_expt = expt['E_eV']

        # Find closest computed resonance
        best_match = None
        best_delta = float('inf')
        for res in resonances:
            delta = abs(res['E_eV_above_gs'] - E_expt)
            if delta < best_delta:
                best_delta = delta
                best_match = res

        if best_match is not None:
            comparison.append({
                'label': label,
                'experimental': expt,
                'computed': best_match,
                'delta_E_eV': best_match['E_eV_above_gs'] - E_expt,
            })

    return comparison


def build_weighted_coupling_graph(
    l_max: int = 2,
    n_channels: int = 5,
    Z: float = 2.0,
    n_alpha: int = 100,
    R_ref: float = 2.0,
) -> Dict:
    """
    Build a continuous weighted Gaunt coupling graph — no threshold parameter.

    Edge weights are the coupling strengths C(μ,ν). Edge distances are
    d(μ,ν) = 1/C(μ,ν) for C > 0, infinity otherwise. Weighted shortest
    paths are computed via Dijkstra's algorithm.

    This is the threshold-free analogue of build_gaunt_coupling_graph.
    The binary graph emerges as a projection of this continuous structure,
    analogous to how the Schrödinger equation emerges from projecting the
    dimensionless S³ topology (Paper 7).

    Parameters
    ----------
    l_max : int
        Maximum partial wave.
    n_channels : int
        Number of channels.
    Z : float
        Nuclear charge.
    n_alpha : int
        Angular grid points.
    R_ref : float
        Reference hyperradius for computing coupling strengths.

    Returns
    -------
    result : dict
        Keys:
        - 'coupling_strength': ndarray (n_channels, n_channels)
        - 'weighted_distance': ndarray (n_channels, n_channels) — Dijkstra distances
        - 'predecessors': ndarray (n_channels, n_channels) — shortest-path predecessors
        - 'channel_info': list of dicts with dominant l for each channel
    """
    # Solve angular at R_ref to get eigenvectors
    mu, vecs = solve_angular(R_ref, Z, l_max, n_alpha, n_channels)

    n_l = l_max + 1

    # Compute coupling strength and channel info (same as build_gaunt_coupling_graph)
    coupling = np.zeros((n_channels, n_channels))
    channel_info = []

    for ch in range(n_channels):
        l_weights = np.zeros(n_l)
        for l in range(n_l):
            l_weights[l] = np.sum(vecs[ch, l * n_alpha:(l + 1) * n_alpha]**2)
        dominant_l = np.argmax(l_weights)
        channel_info.append({
            'channel': ch,
            'dominant_l': int(dominant_l),
            'l_weights': l_weights,
            'mu': mu[ch],
        })

    for mu_idx in range(n_channels):
        for nu_idx in range(mu_idx + 1, n_channels):
            strength = 0.0
            for l1 in range(n_l):
                w1 = channel_info[mu_idx]['l_weights'][l1]
                for l2 in range(n_l):
                    w2 = channel_info[nu_idx]['l_weights'][l2]
                    for k in range(abs(l1 - l2), l1 + l2 + 1):
                        g = gaunt_integral(l1, k, l2)
                        strength += w1 * w2 * abs(g)

            coupling[mu_idx, nu_idx] = strength
            coupling[nu_idx, mu_idx] = strength

    # Build distance matrix: d(μ,ν) = 1/C(μ,ν) for C > 0
    distance_matrix = np.zeros((n_channels, n_channels))
    for i in range(n_channels):
        for j in range(n_channels):
            if i == j:
                distance_matrix[i, j] = 0.0
            elif coupling[i, j] > 0:
                distance_matrix[i, j] = 1.0 / coupling[i, j]
            # Zero entries stay zero — will be treated as no edge by csr_matrix

    # Use Dijkstra's algorithm via scipy.sparse.csgraph
    graph = csr_matrix(distance_matrix)
    weighted_dist, predecessors = dijkstra_shortest_path(
        graph, directed=False, return_predecessors=True
    )

    return {
        'coupling_strength': coupling,
        'weighted_distance': weighted_dist,
        'predecessors': predecessors,
        'channel_info': channel_info,
    }


def compute_path_coupling(
    weighted_graph: Dict,
    source: int,
    target: int = 0,
) -> float:
    """
    Compute the path coupling V_path from source to target (continuum).

    V_path = product of edge coupling strengths along the Dijkstra shortest
    weighted path from source to target. This is a continuous analogue of
    the binary graph distance: instead of counting hops, we multiply the
    coupling strengths along the optimal route.

    For Fermi's golden rule: Γ ∝ |V_path|².

    Parameters
    ----------
    weighted_graph : dict
        Output from build_weighted_coupling_graph.
    source : int
        Source channel index.
    target : int
        Target channel index (default 0, the continuum).

    Returns
    -------
    V_path : float
        Product of edge weights along shortest weighted path.
        Returns 0.0 if no path exists.
    """
    predecessors = weighted_graph['predecessors']
    coupling = weighted_graph['coupling_strength']
    n_ch = coupling.shape[0]

    if source == target:
        return 1.0

    # Reconstruct path from predecessors
    path = []
    current = source
    visited = set()
    while current != target:
        if current in visited or current < 0 or current >= n_ch:
            return 0.0  # No path exists
        visited.add(current)
        path.append(current)
        current = predecessors[target, current]

    path.append(target)

    # Compute product of edge couplings along path
    V_path = 1.0
    for i in range(len(path) - 1):
        edge_coupling = coupling[path[i], path[i + 1]]
        if edge_coupling <= 0:
            return 0.0
        V_path *= edge_coupling

    return V_path


def threshold_robustness_study(
    l_max: int = 2,
    n_channels: int = 5,
    Z: float = 2.0,
    n_alpha: int = 100,
    R_ref: float = 2.0,
    thresholds: Optional[List[float]] = None,
) -> Dict:
    """
    Study how binary graph distances depend on the coupling threshold.

    Computes the binary graph (build_gaunt_coupling_graph) at multiple
    thresholds and checks whether the distance ordering is preserved.
    The weighted graph (build_weighted_coupling_graph) is threshold-free
    and serves as the reference.

    Parameters
    ----------
    l_max, n_channels, Z, n_alpha, R_ref : solver parameters
    thresholds : list of float, optional
        Coupling thresholds to scan. Default: [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0].

    Returns
    -------
    result : dict
        Keys:
        - 'thresholds': list of float
        - 'binary_distances': list of ndarray — distance matrices at each threshold
        - 'weighted_distances': ndarray — Dijkstra weighted distances (reference)
        - 'ordering_preserved': list of bool — whether rank ordering is preserved
        - 'coupling_strength': ndarray — raw coupling matrix
        - 'channel_info': list of dicts
    """
    if thresholds is None:
        thresholds = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]

    # Compute the weighted (threshold-free) graph
    wg = build_weighted_coupling_graph(l_max, n_channels, Z, n_alpha, R_ref)
    weighted_dist = wg['weighted_distance']
    coupling = wg['coupling_strength']
    channel_info = wg['channel_info']

    # Reference ordering: rank channels by weighted distance to channel 0
    ref_distances = weighted_dist[0, :]
    ref_order = np.argsort(ref_distances)

    binary_distances_list = []
    ordering_preserved = []

    for thresh in thresholds:
        adjacency = coupling > thresh
        dist = _shortest_paths(adjacency, n_channels)
        binary_distances_list.append(dist)

        # Check if distance ordering to channel 0 is preserved
        bin_dist_to_0 = dist[0, :].astype(float)
        # Replace -1 (unreachable) with infinity
        bin_dist_to_0[bin_dist_to_0 < 0] = np.inf

        bin_order = np.argsort(bin_dist_to_0)

        # Ordering preserved if rank correlation is monotonic
        # (channels closer in weighted distance are also closer in binary distance)
        # Use a simple check: for all i < j in ref_order,
        # bin_dist[ref_order[i]] <= bin_dist[ref_order[j]]
        preserved = True
        for i in range(len(ref_order)):
            for j in range(i + 1, len(ref_order)):
                ci, cj = ref_order[i], ref_order[j]
                if ref_distances[ci] < ref_distances[cj]:
                    if bin_dist_to_0[ci] > bin_dist_to_0[cj]:
                        preserved = False
                        break
            if not preserved:
                break

        ordering_preserved.append(preserved)

    return {
        'thresholds': thresholds,
        'binary_distances': binary_distances_list,
        'weighted_distances': weighted_dist,
        'ordering_preserved': ordering_preserved,
        'coupling_strength': coupling,
        'channel_info': channel_info,
    }


def predict_widths(
    weighted_graph: Dict,
    channel_assignments: Dict[str, int],
    calibration_label: str = '2s2_1S',
) -> Dict:
    """
    Predict autoionization widths from path coupling on the weighted graph.

    Uses ONE experimental calibration point to set the overall scale:
        Γ_pred(μ) = Γ_expt(cal) × |V_path(μ)|² / |V_path(cal)|²

    This is the graph analogue of Fermi's golden rule: the decay rate is
    proportional to the square of the effective coupling to the continuum,
    computed as the product of edge weights along the shortest weighted path.

    Parameters
    ----------
    weighted_graph : dict
        Output from build_weighted_coupling_graph.
    channel_assignments : dict
        Map from EXPERIMENTAL_RESONANCES label to channel index.
        E.g., {'2s2_1S': 1, '2s3s_1S': 3, '2p2_1S': 2, ...}
    calibration_label : str
        Which experimental state to use as the calibration point.

    Returns
    -------
    result : dict
        Keys:
        - 'predictions': dict mapping label → {Gamma_pred_eV, V_path, d_weighted, ...}
        - 'calibration': dict with calibration details
        - 'ratio_2s2_2s3s': predicted ratio Γ(2s²)/Γ(2s3s) (expt: 3.3)
    """
    expt = EXPERIMENTAL_RESONANCES

    # Compute path couplings for all assigned channels
    path_couplings = {}
    for label, ch in channel_assignments.items():
        V_path = compute_path_coupling(weighted_graph, source=ch, target=0)
        d_w = weighted_graph['weighted_distance'][0, ch]
        path_couplings[label] = {
            'V_path': V_path,
            'd_weighted': d_w,
            'channel': ch,
        }

    # Calibrate with one point
    cal = path_couplings[calibration_label]
    V_cal_sq = cal['V_path']**2
    Gamma_cal = expt[calibration_label]['Gamma_eV']

    if V_cal_sq < 1e-30:
        raise ValueError(
            f"Calibration channel {calibration_label} has zero path coupling"
        )

    # Predict all widths
    predictions = {}
    for label, pc in path_couplings.items():
        V_sq = pc['V_path']**2
        Gamma_pred = Gamma_cal * V_sq / V_cal_sq if V_cal_sq > 0 else 0.0
        predictions[label] = {
            'Gamma_pred_eV': Gamma_pred,
            'Gamma_expt_eV': expt[label]['Gamma_eV'] if label in expt else None,
            'V_path': pc['V_path'],
            'd_weighted': pc['d_weighted'],
            'channel': pc['channel'],
        }

    # Within-sector ratio: 2s²/2s3s (both s-wave, expt ratio = 3.3)
    ratio = None
    if '2s2_1S' in predictions and '2s3s_1S' in predictions:
        G_2s2 = predictions['2s2_1S']['Gamma_pred_eV']
        G_2s3s = predictions['2s3s_1S']['Gamma_pred_eV']
        if G_2s3s > 0:
            ratio = G_2s2 / G_2s3s

    return {
        'predictions': predictions,
        'calibration': {
            'label': calibration_label,
            'Gamma_expt_eV': Gamma_cal,
            'V_path': cal['V_path'],
            'V_path_sq': V_cal_sq,
        },
        'ratio_2s2_2s3s': ratio,
    }


def find_avoided_crossings(
    R_grid: np.ndarray,
    mu_curves: np.ndarray,
    Z: float = 2.0,
    R_search_min: float = 0.5,
    R_search_max: float = 15.0,
) -> List[Dict]:
    """
    Find avoided crossings between adiabatic potential curves.

    For each pair (μ, ν) with μ < ν, locates the hyperradius R_c where
    the effective potential gap |V_eff_μ(R) - V_eff_ν(R)| is minimized.
    Uses cubic spline interpolation for sub-grid resolution.

    The avoided crossing geometry — characterized by (R_c, δ, |ΔF|) —
    determines autoionization widths via the Landau-Zener formula,
    without any integrals over radial wavefunctions.

    Parameters
    ----------
    R_grid : ndarray of shape (N_R,)
        Hyperradius grid.
    mu_curves : ndarray of shape (n_channels, N_R)
        Angular eigenvalues mu(R) for each channel.
    Z : float
        Nuclear charge (for threshold information).
    R_search_min : float
        Minimum R for crossing search (bohr).
    R_search_max : float
        Maximum R for crossing search (bohr).

    Returns
    -------
    crossings : list of dict
        One entry per pair (μ, ν), each with keys:
        - 'channels': (mu, nu) pair indices
        - 'R_c': location of minimum gap (bohr)
        - 'delta': minimum gap in V_eff (Ha)
        - 'V_coupling': estimated coupling = delta/2 (Ha)
        - 'slope_diff': |dV_eff_mu/dR - dV_eff_nu/dR| at R_c (Ha/bohr)
        - 'V_eff_mu_at_Rc': V_eff of channel mu at crossing
        - 'V_eff_nu_at_Rc': V_eff of channel nu at crossing
    """
    n_channels = mu_curves.shape[0]

    # Build V_eff splines for each channel
    V_eff_data = np.zeros_like(mu_curves)
    for ch in range(n_channels):
        V_eff_data[ch] = effective_potential(R_grid, mu_curves[ch])

    V_eff_splines = []
    for ch in range(n_channels):
        V_eff_splines.append(CubicSpline(R_grid, V_eff_data[ch]))

    # Restrict search range
    R_lo = max(R_search_min, R_grid[0])
    R_hi = min(R_search_max, R_grid[-1])

    crossings = []

    for mu_idx in range(n_channels):
        for nu_idx in range(mu_idx + 1, n_channels):
            spl_mu = V_eff_splines[mu_idx]
            spl_nu = V_eff_splines[nu_idx]

            # Gap function: |V_eff_mu(R) - V_eff_nu(R)|
            def gap_func(R: float) -> float:
                return abs(float(spl_mu(R)) - float(spl_nu(R)))

            # Find minimum gap via minimize_scalar on the search interval
            result = minimize_scalar(
                gap_func, bounds=(R_lo, R_hi), method='bounded',
                options={'xatol': 1e-6}
            )

            R_c = float(result.x)
            delta = float(result.fun)

            # Slopes at R_c via spline derivatives
            dV_mu = float(spl_mu(R_c, 1))  # first derivative
            dV_nu = float(spl_nu(R_c, 1))
            slope_diff = abs(dV_mu - dV_nu)

            crossings.append({
                'channels': (mu_idx, nu_idx),
                'R_c': R_c,
                'delta': delta,
                'V_coupling': delta / 2.0,
                'slope_diff': slope_diff,
                'V_eff_mu_at_Rc': float(spl_mu(R_c)),
                'V_eff_nu_at_Rc': float(spl_nu(R_c)),
            })

    return crossings


def landau_zener_width(
    delta: float,
    slope_diff: float,
) -> float:
    """
    Landau-Zener estimate of autoionization width from avoided crossing geometry.

    In the Landau-Zener model, the non-adiabatic transition probability through
    an avoided crossing depends on the gap δ and the difference of diabatic
    slopes |ΔF|. The width (inverse lifetime) scales as:

        Γ ~ δ² / |ΔF|

    This is parameter-free for width *ratios*: the unknown prefactor cancels.

    Parameters
    ----------
    delta : float
        Minimum gap between adiabatic curves (Ha).
    slope_diff : float
        |dV_mu/dR - dV_nu/dR| at the crossing point (Ha/bohr).

    Returns
    -------
    width : float
        Relative width in Ha (arbitrary overall scale, meaningful for ratios).
        Returns 0.0 if slope_diff is zero.
    """
    if slope_diff < 1e-30:
        return 0.0
    return delta**2 / slope_diff


def feshbach_width(
    delta: float,
    E_res: float,
    V_eff_continuum: float,
) -> float:
    """
    Feshbach projection estimate of autoionization width.

    The Feshbach half-width for a resonance coupling to a continuum channel is:
        Γ/2 = π |V_coupling|² × ρ(E)

    where V_coupling = δ/2 (half the avoided crossing gap) and ρ(E) is the
    density of states in the continuum channel:
        ρ(E) = 1/(π × v_group) ≈ 1/(π × sqrt(2(E_res - V_eff_continuum)))

    So:
        Γ = δ² / (2 × sqrt(2(E_res - V_eff_continuum(R_c))))

    All quantities are eigenvalue-derived — no radial wavefunction integrals.

    Parameters
    ----------
    delta : float
        Minimum gap between adiabatic curves (Ha).
    E_res : float
        Resonance energy (Ha).
    V_eff_continuum : float
        Effective potential of the continuum channel at R_c (Ha).

    Returns
    -------
    width : float
        Estimated width (Ha). Returns 0.0 if E_res <= V_eff_continuum.
    """
    kinetic_energy = E_res - V_eff_continuum
    if kinetic_energy <= 0:
        return 0.0
    v_group = np.sqrt(2.0 * kinetic_energy)
    return delta**2 / (2.0 * v_group)


def analyze_avoided_crossings(
    Z: float = 2.0,
    l_max: int = 2,
    n_alpha: int = 100,
    n_channels: int = 5,
    N_R: int = 400,
    R_max: float = 15.0,
    channel_assignments: Optional[Dict[str, Tuple[int, int]]] = None,
) -> Dict:
    """
    Full avoided crossing analysis: compute adiabatic curves, find crossings,
    predict widths via Landau-Zener and Feshbach formulas.

    Parameters
    ----------
    Z, l_max, n_alpha, n_channels : angular solver parameters
    N_R : int
        Number of R grid points (dense grid for accurate crossing detection).
    R_max : float
        Maximum hyperradius.
    channel_assignments : dict, optional
        Map from experimental resonance label to (resonant_channel, continuum_channel).
        Default: {'2s2_1S': (1, 0), '2s3s_1S': (3, 0), '2p2_1S': (2, 0),
                  '2p2_1D': (4, 0), '2s2p_1P': (2, 0)}

    Returns
    -------
    result : dict
        Keys:
        - 'R_grid': ndarray
        - 'mu_curves': ndarray (n_channels, N_R)
        - 'V_eff': ndarray (n_channels, N_R)
        - 'crossings': list of all pairwise crossing dicts
        - 'resonance_widths': dict mapping label to width analysis
        - 'ratio_2s2_2s3s_LZ': Landau-Zener width ratio
        - 'ratio_2s2_2s3s_Feshbach': Feshbach width ratio
    """
    if channel_assignments is None:
        channel_assignments = {
            '2s2_1S': (1, 0),
            '2s3s_1S': (3, 0),
            '2p2_1S': (2, 0),
            '2p2_1D': (4, 0),
            '2s2p_1P': (2, 0),
        }

    # Dense R grid concentrated where crossings occur
    R_grid = np.concatenate([
        np.linspace(0.3, 1.0, N_R // 4),
        np.linspace(1.0, 4.0, N_R // 2),
        np.linspace(4.0, R_max, N_R // 4 + 1),
    ])
    R_grid = np.unique(R_grid)

    # Compute adiabatic curves
    mu_curves = compute_adiabatic_curve(R_grid, Z, l_max, n_alpha, n_channels)

    # Compute V_eff for each channel
    V_eff = np.zeros_like(mu_curves)
    for ch in range(n_channels):
        V_eff[ch] = effective_potential(R_grid, mu_curves[ch])

    # Find all avoided crossings
    crossings = find_avoided_crossings(R_grid, mu_curves, Z)

    # Build lookup: (mu, nu) -> crossing data
    crossing_lookup = {}
    for xing in crossings:
        pair = xing['channels']
        crossing_lookup[pair] = xing

    # Analyze widths for each experimental resonance
    resonance_widths = {}
    for label, (res_ch, cont_ch) in channel_assignments.items():
        pair = (min(res_ch, cont_ch), max(res_ch, cont_ch))
        if pair not in crossing_lookup:
            continue

        xing = crossing_lookup[pair]
        delta = xing['delta']
        slope_diff = xing['slope_diff']

        # Landau-Zener width
        Gamma_LZ = landau_zener_width(delta, slope_diff)

        # Feshbach width — need resonance energy
        # Use experimental energy converted to Ha (relative to He ground state)
        if label in EXPERIMENTAL_RESONANCES:
            E_res_eV = EXPERIMENTAL_RESONANCES[label]['E_eV']
            E_res_Ha = E_HE_EXACT + E_res_eV / HARTREE_TO_EV
        else:
            E_res_Ha = xing['V_eff_mu_at_Rc']  # fallback

        V_cont = xing['V_eff_nu_at_Rc'] if cont_ch == pair[1] else xing['V_eff_mu_at_Rc']
        Gamma_Feshbach = feshbach_width(delta, E_res_Ha, V_cont)

        resonance_widths[label] = {
            'channels': (res_ch, cont_ch),
            'R_c': xing['R_c'],
            'delta': delta,
            'V_coupling': xing['V_coupling'],
            'slope_diff': slope_diff,
            'Gamma_LZ': Gamma_LZ,
            'Gamma_Feshbach': Gamma_Feshbach,
            'Gamma_expt_eV': EXPERIMENTAL_RESONANCES.get(label, {}).get('Gamma_eV'),
        }

    # Within-sector ratios
    ratio_LZ = None
    ratio_Feshbach = None
    if '2s2_1S' in resonance_widths and '2s3s_1S' in resonance_widths:
        w1 = resonance_widths['2s2_1S']
        w2 = resonance_widths['2s3s_1S']
        if w2['Gamma_LZ'] > 0:
            ratio_LZ = w1['Gamma_LZ'] / w2['Gamma_LZ']
        if w2['Gamma_Feshbach'] > 0:
            ratio_Feshbach = w1['Gamma_Feshbach'] / w2['Gamma_Feshbach']

    return {
        'R_grid': R_grid,
        'mu_curves': mu_curves,
        'V_eff': V_eff,
        'crossings': crossings,
        'resonance_widths': resonance_widths,
        'ratio_2s2_2s3s_LZ': ratio_LZ,
        'ratio_2s2_2s3s_Feshbach': ratio_Feshbach,
    }


def plot_avoided_crossings(
    analysis: Dict,
    save_path: Optional[str] = None,
) -> None:
    """
    Plot adiabatic potential curves with avoided crossings annotated.

    Parameters
    ----------
    analysis : dict
        Output from analyze_avoided_crossings.
    save_path : str, optional
        Path to save figure.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    R_grid = analysis['R_grid']
    V_eff = analysis['V_eff']
    n_channels = V_eff.shape[0]
    crossings = analysis['crossings']
    resonance_widths = analysis.get('resonance_widths', {})

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    labels = ['ch0 (s, continuum)', 'ch1 (s)', 'ch2 (p)', 'ch3 (s)', 'ch4 (p)']

    for ch in range(n_channels):
        label = labels[ch] if ch < len(labels) else f'ch{ch}'
        color = colors[ch % len(colors)]
        ax.plot(R_grid, V_eff[ch], color=color, label=label, linewidth=1.5)

    # Annotate avoided crossings relevant to resonances
    for label, rw in resonance_widths.items():
        R_c = rw['R_c']
        delta = rw['delta']
        ch_res, ch_cont = rw['channels']
        ax.plot(R_c, (rw.get('V_eff_mu_at_Rc', 0) + rw.get('V_eff_nu_at_Rc', 0)) / 2
                if 'V_eff_mu_at_Rc' in rw else 0,
                'ko', markersize=8, zorder=5)

        # Find crossing data for annotation
        pair = (min(ch_res, ch_cont), max(ch_res, ch_cont))
        for xing in crossings:
            if xing['channels'] == pair:
                mid_V = (xing['V_eff_mu_at_Rc'] + xing['V_eff_nu_at_Rc']) / 2
                ax.plot(R_c, mid_V, 'ko', markersize=8, zorder=5)
                ax.annotate(
                    f'{label}\nδ={delta:.4f} Ha\nR_c={R_c:.2f}',
                    xy=(R_c, mid_V),
                    xytext=(R_c + 0.5, mid_V + 0.3),
                    fontsize=7,
                    arrowprops=dict(arrowstyle='->', color='gray'),
                )
                break

    # He+ threshold
    ax.axhline(E_HE_PLUS, color='gray', linestyle='--', alpha=0.5,
               label='He+ threshold (-2.0 Ha)')

    ax.set_xlabel('Hyperradius R (bohr)')
    ax.set_ylabel('V_eff(R) (Ha)')
    ax.set_title('Adiabatic Potential Curves — Avoided Crossings')
    ax.set_xlim(0.3, 12)
    ax.set_ylim(-5, 3)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved avoided crossing plot to {save_path}")
    plt.close()


def plot_stabilization_diagram(
    scan_result: Dict,
    resonances: Optional[List[Dict]] = None,
    save_path: Optional[str] = None,
) -> None:
    """
    Plot the stabilization diagram (eigenvalues vs R_max).

    Parameters
    ----------
    scan_result : dict
        Output from stabilization_scan.
    resonances : list of dict, optional
        Detected resonances to highlight.
    save_path : str, optional
        Path to save figure.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    R_max_values = scan_result['R_max_values']

    for i, (R_max, E_arr) in enumerate(zip(R_max_values, scan_result['eigenvalues'])):
        ax.plot([R_max] * len(E_arr), E_arr, 'k.', markersize=3, alpha=0.5)

    # He+ threshold
    ax.axhline(E_HE_PLUS, color='red', linestyle='--', alpha=0.7,
               label=f'He$^+$ threshold ({E_HE_PLUS:.1f} Ha)')

    # He ground state
    ax.axhline(E_HE_EXACT, color='blue', linestyle='--', alpha=0.7,
               label=f'He g.s. ({E_HE_EXACT:.4f} Ha)')

    if resonances:
        for res in resonances:
            ax.axhline(res['E_Ha'], color='green', linestyle=':',
                       alpha=0.5)
            ax.annotate(
                f"E={res['E_eV_above_gs']:.1f} eV",
                xy=(R_max_values[-1], res['E_Ha']),
                fontsize=8, color='green',
            )

    ax.set_xlabel('R_max (bohr)')
    ax.set_ylabel('Energy (Ha)')
    ax.set_title('Stabilization Diagram — Doubly-Excited He Resonances')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved stabilization diagram to {save_path}")
    plt.close()
