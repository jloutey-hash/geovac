"""
Wilson-Dirac operator on the Dirac packing graph with S3 spin connection.

Combines two ingredients:
1. Dirac packing graph topology (nodes = (n, kappa, m_j) with correct
   degeneracy pattern 2(n+1)(n+2) per Fock shell)
2. SU(2) parallel transport as spin connection on edges

The Fock graph test (wilson_dirac_s3_test.py) showed the scalar graph is
too sparse and has wrong state-counting. This test uses the spinor quantum
number graph which has the right state-counting for Dirac.

S3 embedding for Dirac states: each (n, kappa, m_j) is mapped to a point
on S3 via a generalization of the Fock projection that accounts for the
spinor structure.
"""

import numpy as np
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# === SU(2) machinery ===

sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
SIGMA = [sigma_x, sigma_y, sigma_z]


def s3_to_su2(x):
    return x[0] * I2 + 1j * (x[1] * sigma_x + x[2] * sigma_y + x[3] * sigma_z)


def su2_parallel_transport(g_x, g_y):
    return g_y @ np.conj(g_x.T)


def geodesic_distance_s3(x, y):
    dot = np.clip(np.dot(x, y), -1.0, 1.0)
    return np.arccos(dot)


def tangent_vector_s3(x, y):
    proj = y - np.dot(x, y) * x
    norm = np.linalg.norm(proj)
    if norm < 1e-14:
        return np.zeros(4)
    return proj / norm


def tangent_to_gamma(x, t):
    """Map tangent vector at x on S3 to 2x2 gamma matrix via left-invariant frame."""
    g = s3_to_su2(x)
    frame_r4 = []
    for k in range(3):
        mat = g @ (1j * SIGMA[k])
        v0 = mat[0, 0].real
        v1 = mat[0, 1].imag
        v2 = mat[1, 0].real
        v3 = mat[0, 0].imag
        frame_r4.append(np.array([v0, v1, v2, v3]))
    components = np.array([np.dot(t, e) for e in frame_r4])
    gamma = sum(components[k] * SIGMA[k] for k in range(3))
    return gamma


# === Dirac state generation ===

def dirac_states(n_max_fock):
    """Generate all (n, kappa, m_j) Dirac states for n=1..n_max_fock."""
    states = []
    for n in range(1, n_max_fock + 1):
        for l in range(n):
            # kappa = -(l+1), j = l + 1/2
            kappa = -(l + 1)
            j = l + 0.5
            for m_j_2 in range(int(-2*j), int(2*j) + 1, 2):
                states.append((n, kappa, m_j_2 / 2.0))
            # kappa = +l, j = l - 1/2 (only if l >= 1)
            if l >= 1:
                kappa = l
                j = l - 0.5
                for m_j_2 in range(int(-2*j), int(2*j) + 1, 2):
                    states.append((n, kappa, m_j_2 / 2.0))
    return states


def dirac_state_to_s3(n, kappa, m_j):
    """
    Map Dirac state (n, kappa, m_j) to a point on S3.

    Generalization of the Fock map:
    - chi_n = 2*arctan(1/n) from Fock projection (radial -> S3 latitude)
    - psi_1 encodes kappa (orbital/spin structure)
    - psi_2 encodes m_j (magnetic projection)

    For the spinor embedding, we use:
    - chi from the principal quantum number (same as scalar)
    - psi_1 from kappa: maps the range of kappa values at each n onto [0, pi]
    - psi_2 from m_j: maps the range of m_j values onto [0, 2pi)
    """
    chi = 2.0 * np.arctan(1.0 / n)

    # kappa ranges from -(n) to +(n-1), excluding 0
    # Map kappa to psi_1 in [0, pi]
    # Total kappa count at shell n: 2n-1 values
    # Use a uniform spacing
    all_kappas = []
    for l in range(n):
        all_kappas.append(-(l+1))
        if l >= 1:
            all_kappas.append(l)
    all_kappas.sort()
    n_kappas = len(all_kappas)
    if n_kappas > 1:
        k_idx = all_kappas.index(kappa)
        psi_1 = np.pi * k_idx / (n_kappas - 1)
    else:
        psi_1 = 0.0

    # m_j ranges from -j to +j where j = |kappa| - 1/2
    j = abs(kappa) - 0.5
    if j > 0:
        psi_2 = 2.0 * np.pi * (m_j + j) / (2 * j + 1)
    else:
        psi_2 = 0.0

    cos_half = np.cos(chi / 2.0)
    sin_half = np.sin(chi / 2.0)
    x = np.array([
        cos_half * np.cos(psi_1),
        cos_half * np.sin(psi_1),
        sin_half * np.cos(psi_2),
        sin_half * np.sin(psi_2),
    ])
    # Normalize (should already be ~1, but be safe)
    return x / np.linalg.norm(x)


# === Graph construction ===

def build_dirac_graph(n_max, edge_rule="standard"):
    """
    Build graph on Dirac quantum numbers (n, kappa, m_j).

    Edge rules:
    - "standard": m_j +/- 1 angular, n +/- 1 radial (same kappa, m_j)
    - "dipole": E1 selection rules (Delta l = +/-1, |Delta m_j| <= 1)
    - "full": standard + dipole combined
    """
    states = dirac_states(n_max)
    state_idx = {s: i for i, s in enumerate(states)}
    V = len(states)
    adj = np.zeros((V, V), dtype=float)
    edges = []

    for s in states:
        n, kappa, m_j = s
        i = state_idx[s]
        j_val = abs(kappa) - 0.5

        if edge_rule in ["standard", "full"]:
            # Angular: m_j +/- 1
            if m_j + 1 <= j_val:
                t = (n, kappa, m_j + 1)
                if t in state_idx:
                    jj = state_idx[t]
                    if adj[i, jj] == 0:
                        adj[i, jj] = 1.0
                        adj[jj, i] = 1.0
            # Radial: n +/- 1
            if n < n_max:
                t = (n + 1, kappa, m_j)
                if t in state_idx:
                    jj = state_idx[t]
                    if adj[i, jj] == 0:
                        adj[i, jj] = 1.0
                        adj[jj, i] = 1.0

        if edge_rule in ["dipole", "full"]:
            # E1 dipole: Delta l = +/-1, |Delta m_j| <= 1
            if kappa < 0:
                l_val = -kappa - 1
            else:
                l_val = kappa
            for dn in [0, 1]:
                n_new = n + dn
                if n_new < 1 or n_new > n_max:
                    continue
                for dl in [-1, 1]:
                    l_new = l_val + dl
                    if l_new < 0 or l_new >= n_new:
                        continue
                    for dm in [-1, 0, 1]:
                        new_mj = m_j + dm
                        for kp in [-(l_new + 1), l_new]:
                            if kp == 0:
                                continue
                            j_new = abs(kp) - 0.5
                            if abs(new_mj) <= j_new:
                                t = (n_new, kp, new_mj)
                                if t in state_idx and t != s:
                                    jj = state_idx[t]
                                    if adj[i, jj] == 0:
                                        adj[i, jj] = 1.0
                                        adj[jj, i] = 1.0

    # Extract edge list
    for i in range(V):
        for j in range(i+1, V):
            if adj[i, j] > 0.5:
                edges.append((i, j))

    # S3 embedding
    s3_coords = np.array([dirac_state_to_s3(*s) for s in states])
    su2_mats = np.array([s3_to_su2(s3_coords[i]) for i in range(V)])

    return states, edges, adj, s3_coords, su2_mats, state_idx


# === Wilson-Dirac operators ===

def build_wilson_dirac_pure_transport(V, edges, su2_mats):
    """Pure SU(2) parallel transport, Hermitian."""
    dim = 2 * V
    D = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        U_ij = su2_parallel_transport(su2_mats[i], su2_mats[j])
        U_ji = su2_parallel_transport(su2_mats[j], su2_mats[i])
        D[2*i:2*i+2, 2*j:2*j+2] += U_ij
        D[2*j:2*j+2, 2*i:2*i+2] += U_ji
    return D


def build_wilson_dirac_gamma_U(V, edges, s3_coords, su2_mats):
    """gamma * U Wilson-Dirac, anti-Hermitian contribution."""
    dim = 2 * V
    D = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        x_i, x_j = s3_coords[i], s3_coords[j]
        t_ij = tangent_vector_s3(x_i, x_j)
        t_ji = tangent_vector_s3(x_j, x_i)
        if np.linalg.norm(t_ij) < 1e-14:
            continue
        gamma_ij = tangent_to_gamma(x_i, t_ij)
        gamma_ji = tangent_to_gamma(x_j, t_ji)
        U_ij = su2_parallel_transport(su2_mats[i], su2_mats[j])
        U_ji = su2_parallel_transport(su2_mats[j], su2_mats[i])
        D[2*i:2*i+2, 2*j:2*j+2] += gamma_ij @ U_ij
        D[2*j:2*j+2, 2*i:2*i+2] += gamma_ji @ U_ji
    return D


def build_wilson_dirac_dist_norm(V, edges, s3_coords, su2_mats):
    """(1/d) gamma * U, distance normalized."""
    dim = 2 * V
    D = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        x_i, x_j = s3_coords[i], s3_coords[j]
        d = geodesic_distance_s3(x_i, x_j)
        if d < 1e-14:
            continue
        t_ij = tangent_vector_s3(x_i, x_j)
        t_ji = tangent_vector_s3(x_j, x_i)
        if np.linalg.norm(t_ij) < 1e-14:
            continue
        gamma_ij = tangent_to_gamma(x_i, t_ij)
        gamma_ji = tangent_to_gamma(x_j, t_ji)
        U_ij = su2_parallel_transport(su2_mats[i], su2_mats[j])
        U_ji = su2_parallel_transport(su2_mats[j], su2_mats[i])
        D[2*i:2*i+2, 2*j:2*j+2] += gamma_ij @ U_ij / d
        D[2*j:2*j+2, 2*i:2*i+2] += gamma_ji @ U_ji / d
    return D


def build_scalar_laplacian(V, adj):
    """Standard graph Laplacian for comparison (no spin)."""
    degree = np.sum(adj, axis=1)
    L = np.diag(degree) - adj
    return L


# === Analysis ===

def camporesi_higuchi_spectrum(n_max):
    levels = []
    for n in range(0, n_max + 1):
        lam = n + 1.5
        deg = 2 * (n + 1) * (n + 2)
        levels.append((n, lam, deg))
    return levels


def analyze_operator(D, label, ch_targets):
    """Diagonalize and compare to Camporesi-Higuchi targets."""
    herm_err = float(np.max(np.abs(D - D.conj().T)))
    if herm_err > 1e-8:
        evals_c = np.linalg.eigvals(D)
        evals = np.sort(evals_c.real)
        imag_max = float(np.max(np.abs(evals_c.imag)))
    else:
        evals = np.sort(np.linalg.eigvalsh(D))
        imag_max = 0.0

    # Cluster
    tol = 0.02
    clusters = []
    rem = list(evals)
    while rem:
        val = rem[0]
        cl = [v for v in rem if abs(v - val) < tol]
        rem = [v for v in rem if abs(v - val) >= tol]
        clusters.append((float(np.mean(cl)), len(cl)))

    # Match to CH
    matches = []
    for n_ch, lam, deg in ch_targets:
        for sign in [+1, -1]:
            tv = sign * lam
            best_c = min(clusters, key=lambda c: abs(c[0] - tv))
            matches.append({
                'n_CH': n_ch, 'sign': sign,
                'target': float(tv), 'target_deg': deg,
                'found': float(best_c[0]), 'found_deg': best_c[1],
                'error': float(abs(best_c[0] - tv))
            })

    return {
        'label': label,
        'dim': int(D.shape[0]),
        'hermiticity_error': herm_err,
        'max_imag': imag_max,
        'eigenvalues': [float(e) for e in evals],
        'clusters': clusters,
        'matches': matches,
        'spectral_range': [float(evals[0]), float(evals[-1])],
    }


def print_analysis(result):
    print(f"\n{'='*70}")
    print(f"  {result['label']}")
    print(f"  dim={result['dim']}, Hermiticity={result['hermiticity_error']:.2e}")
    if result['max_imag'] > 0:
        print(f"  Max imaginary: {result['max_imag']:.2e}")
    print(f"  Spectral range: [{result['spectral_range'][0]:.4f}, {result['spectral_range'][1]:.4f}]")
    print(f"{'='*70}")
    print(f"\n  Clusters ({len(result['clusters'])}):")
    for val, mult in result['clusters']:
        print(f"    lam = {val:+8.4f}  (mult {mult})")
    print(f"\n  CH comparison:")
    print(f"    {'n':>3s} {'sign':>5s} {'target':>8s} {'found':>8s} {'error':>8s} {'deg_t':>6s} {'deg_f':>6s}")
    for m in result['matches']:
        ok = "ok" if m['error'] < 0.1 and m['target_deg'] == m['found_deg'] else "X"
        print(f"    {m['n_CH']:3d} {m['sign']:+5d} {m['target']:+8.4f} {m['found']:+8.4f} "
              f"{m['error']:8.4f} {m['target_deg']:6d} {m['found_deg']:6d}  {ok}")


def main():
    all_results = {}

    for n_max in [2, 3]:
        print(f"\n{'#'*70}")
        print(f"#  n_max = {n_max} (Fock convention)")
        print(f"{'#'*70}")

        ch_targets = camporesi_higuchi_spectrum(n_max + 1)
        results = {}

        for edge_rule in ["standard", "dipole", "full"]:
            print(f"\n--- Edge rule: {edge_rule} ---")
            states, edges, adj, s3_coords, su2_mats, state_idx = build_dirac_graph(n_max, edge_rule)
            V = len(states)
            E = len(edges)

            # Shell decomposition
            shell_counts = {}
            for n, kappa, m_j in states:
                shell_counts[n] = shell_counts.get(n, 0) + 1
            print(f"  V={V}, E={E}")
            print(f"  Shell sizes: {dict(sorted(shell_counts.items()))}")
            ch_expected = {n: 2*n*(n+1) for n in range(1, n_max+1)}
            print(f"  CH expected (2n(n+1)): {ch_expected}")

            # Check S3 coordinates
            norms = np.linalg.norm(s3_coords, axis=1)
            print(f"  S3 coord norms: min={norms.min():.6f}, max={norms.max():.6f}")

            # Check for coincident points
            min_dist = float('inf')
            for ei, ej in edges[:20]:
                d = geodesic_distance_s3(s3_coords[ei], s3_coords[ej])
                min_dist = min(min_dist, d)
            print(f"  Min edge geodesic distance: {min_dist:.6f}")

            # Scalar Laplacian (for reference)
            L = build_scalar_laplacian(V, adj)
            evals_L = np.sort(np.linalg.eigvalsh(L))
            print(f"\n  Scalar Laplacian spectrum (first 10): {[f'{e:.4f}' for e in evals_L[:10]]}")

            # Wilson-Dirac variants
            variants = [
                ("Pure SU(2) transport", lambda: build_wilson_dirac_pure_transport(V, edges, su2_mats)),
                ("gamma*U Wilson-Dirac", lambda: build_wilson_dirac_gamma_U(V, edges, s3_coords, su2_mats)),
                ("(1/d) gamma*U dist-norm", lambda: build_wilson_dirac_dist_norm(V, edges, s3_coords, su2_mats)),
            ]

            results[edge_rule] = {}
            for vlabel, builder in variants:
                full_label = f"{vlabel} [{edge_rule}]"
                D = builder()
                result = analyze_operator(D, full_label, ch_targets)
                results[edge_rule][vlabel] = result
                print_analysis(result)

        all_results[str(n_max)] = {}
        for rule, variants in results.items():
            all_results[str(n_max)][rule] = {
                k: {kk: vv for kk, vv in v.items() if kk != 'eigenvalues'}
                for k, v in variants.items()
            }

    # Also store eigenvalues for n_max=2 only (compact)
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data', 'wilson_dirac_packing_test.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")


if __name__ == '__main__':
    main()
