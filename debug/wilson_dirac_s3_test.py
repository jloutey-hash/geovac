"""
Wilson-Dirac operator on the Fock graph with S3 spin connection.

The hypothesis: the Dirac spectrum on S3 (Camporesi-Higuchi: |lambda_n| = n+3/2,
degeneracy 2(n+1)(n+2)) can emerge from a graph operator that includes
the spin connection as SU(2) parallel transport on edges.

Key insight: the scalar graph Laplacian only encodes the metric (distances).
The Dirac operator requires BOTH metric AND spin connection (frame rotations).
On S3 = SU(2), parallel transport from point g_x to g_y is U = g_y * g_x^{-1}.
"""

import numpy as np
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hopf_bundle import state_to_s3, fock_chi, state_to_hopf_angles
from geovac.lattice import GeometricLattice


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


# === Graph construction ===

def build_fock_graph(n_max):
    lattice = GeometricLattice(max_n=n_max, nuclear_charge=1)
    states = list(lattice.states)
    state_idx = {s: i for i, s in enumerate(states)}
    s3_coords = np.array([state_to_s3(*s) for s in states])
    su2_mats = np.array([s3_to_su2(s3_coords[i]) for i in range(len(states))])
    edges = []
    adj = lattice.adjacency.toarray()
    V = len(states)
    for i in range(V):
        for j in range(i + 1, V):
            if abs(adj[i, j]) > 1e-12:
                edges.append((i, j))
    return states, edges, s3_coords, su2_mats, state_idx


# === Wilson-Dirac operator variants ===

def build_wilson_dirac_v1(V, edges, s3_coords, su2_mats):
    """V1: Pure SU(2) parallel transport, no gamma matrices."""
    dim = 2 * V
    D = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        U_ij = su2_parallel_transport(su2_mats[i], su2_mats[j])
        U_ji = su2_parallel_transport(su2_mats[j], su2_mats[i])
        D[2*i:2*i+2, 2*j:2*j+2] += U_ij
        D[2*j:2*j+2, 2*i:2*i+2] += U_ji
    return D


def build_wilson_dirac_v2(V, edges, s3_coords, su2_mats):
    """V2: Wilson-Dirac with gamma (tangent direction) * connection."""
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


def build_wilson_dirac_v3(V, edges, s3_coords, su2_mats):
    """V3: Distance-normalized Wilson-Dirac: (1/d) gamma * U."""
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


def build_wilson_dirac_v4(V, edges, s3_coords, su2_mats):
    """V4: Gamma only, no connection (control experiment)."""
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
        D[2*i:2*i+2, 2*j:2*j+2] += gamma_ij
        D[2*j:2*j+2, 2*i:2*i+2] += gamma_ji
    return D


def build_wilson_dirac_v5(V, edges, s3_coords, su2_mats):
    """V5: Lattice QCD convention (1-gamma)U/2, symmetrized."""
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
        D[2*i:2*i+2, 2*j:2*j+2] += 0.5 * (I2 - gamma_ij) @ U_ij
        D[2*j:2*j+2, 2*i:2*i+2] += 0.5 * (I2 - gamma_ji) @ U_ji
    D = (D + D.conj().T) / 2
    return D


# === Analysis ===

def camporesi_higuchi_spectrum(n_max):
    levels = []
    for n in range(0, n_max + 1):
        lam = n + 1.5
        deg = 2 * (n + 1) * (n + 2)
        levels.append({'n_CH': n, 'eigenvalue': lam, 'degeneracy': deg})
    return levels


def analyze_spectrum(D, label, n_max, target_levels):
    herm_err = float(np.max(np.abs(D - D.conj().T)))
    if herm_err > 1e-8:
        evals_complex = np.linalg.eigvals(D)
        evals_sorted = np.sort(evals_complex.real)
        evals_imag_max = float(np.max(np.abs(evals_complex.imag)))
    else:
        evals_sorted = np.sort(np.linalg.eigvalsh(D))
        evals_imag_max = 0.0
    tol = 0.05
    clusters = []
    remaining = list(evals_sorted)
    while remaining:
        val = remaining[0]
        cluster = [v for v in remaining if abs(v - val) < tol]
        remaining = [v for v in remaining if abs(v - val) >= tol]
        clusters.append({'value': float(np.mean(cluster)), 'multiplicity': len(cluster)})
    matches = []
    for target in target_levels:
        lam_target = target['eigenvalue']
        for sign in [+1, -1]:
            target_val = sign * lam_target
            best_match = None
            best_dist = float('inf')
            for c in clusters:
                dist = abs(c['value'] - target_val)
                if dist < best_dist:
                    best_dist = dist
                    best_match = c
            if best_match is not None:
                matches.append({
                    'target': float(target_val),
                    'found': float(best_match['value']),
                    'error': float(best_dist),
                    'target_deg': int(target['degeneracy']),
                    'found_deg': int(best_match['multiplicity'])
                })
    D2 = D @ D
    if herm_err < 1e-8:
        evals_D2 = np.sort(np.linalg.eigvalsh(D2))
    else:
        evals_D2 = np.sort(np.linalg.eigvals(D2).real)
    D2_targets = []
    for t in target_levels:
        D2_targets.append({
            'target': float(t['eigenvalue']**2),
            'degeneracy': 2 * int(t['degeneracy'])
        })
    return {
        'label': label,
        'n_max': n_max,
        'dim': int(D.shape[0]),
        'hermiticity_error': herm_err,
        'max_imag': evals_imag_max,
        'eigenvalues': [float(e) for e in evals_sorted],
        'clusters': [{'value': float(c['value']), 'mult': c['multiplicity']}
                     for c in clusters],
        'matches': matches,
        'D2_eigenvalues': [float(e) for e in evals_D2[:20]],
        'D2_targets': D2_targets,
        'spectral_range': [float(evals_sorted[0]), float(evals_sorted[-1])],
    }


def print_result(result):
    print(f"\n{'='*70}")
    print(f"  {result['label']}")
    print(f"  n_max={result['n_max']}, dim={result['dim']}")
    print(f"  Hermiticity error: {result['hermiticity_error']:.2e}")
    if result['max_imag'] > 0:
        print(f"  Max imaginary part: {result['max_imag']:.2e}")
    sr = result['spectral_range']
    print(f"  Spectral range: [{sr[0]:.4f}, {sr[1]:.4f}]")
    print(f"{'='*70}")
    print(f"\n  Eigenvalue clusters:")
    for c in result['clusters']:
        print(f"    lam = {c['value']:+8.4f}  (multiplicity {c['mult']})")
    print(f"\n  Camporesi-Higuchi comparison:")
    # header
    print("    %10s  %10s  %10s  %8s  %8s" % ("Target", "Found", "Error", "Deg_t", "Deg_f"))
    for m in result['matches']:
        marker = ' ok' if m['error'] < 0.1 and m['target_deg'] == m['found_deg'] else ' X'
        print(f"    {m['target']:+10.4f}  {m['found']:+10.4f}  {m['error']:10.4f}  {m['target_deg']:8d}  {m['found_deg']:8d}{marker}")
    print(f"\n  D^2 lowest eigenvalues:")
    for i, e in enumerate(result['D2_eigenvalues'][:10]):
        print(f"    D^2[{i}] = {e:.6f}")
    d2t = result['D2_targets'][:5]
    print(f"\n  D^2 targets (CH): ", end="")
    for t in d2t:
        print(f"{t['target']:.4f} (deg {t['degeneracy']}), ", end="")
    print()


def main():
    results = {}
    for n_max in [2, 3]:
        print(f"\n{'#'*70}")
        print(f"#  n_max = {n_max}")
        print(f"{'#'*70}")
        states, edges, s3_coords, su2_mats, state_idx = build_fock_graph(n_max)
        V = len(states)
        E = len(edges)
        print(f"\n  Graph: V={V}, E={E}")
        for i in range(V):
            g = su2_mats[i]
            det = np.linalg.det(g)
            unitarity = np.max(np.abs(g @ np.conj(g.T) - I2))
            if abs(det - 1.0) > 1e-10 or unitarity > 1e-10:
                print(f"  WARNING: state {states[i]} det={det:.6f}")
        print(f"\n  Edge geodesic distances (first 10):")
        for (i, j) in edges[:10]:
            d = geodesic_distance_s3(s3_coords[i], s3_coords[j])
            print(f"    {states[i]} -> {states[j]}: d = {d:.4f}")
        target = camporesi_higuchi_spectrum(n_max + 1)
        variants = [
            ("V1: Pure SU(2) transport (no gamma)", build_wilson_dirac_v1),
            ("V2: gamma*U (Wilson-Dirac)", build_wilson_dirac_v2),
            ("V3: (1/d) gamma*U (distance-normalized)", build_wilson_dirac_v3),
            ("V4: gamma only (no connection, control)", build_wilson_dirac_v4),
            ("V5: Lattice QCD (1-gamma)U/2 convention", build_wilson_dirac_v5),
        ]
        results[n_max] = {}
        for label, builder in variants:
            D = builder(V, edges, s3_coords, su2_mats)
            result = analyze_spectrum(D, label, n_max, target)
            results[n_max][label] = result
            print_result(result)
    serializable = {}
    for n_max, variants in results.items():
        serializable[str(n_max)] = {}
        for label, result in variants.items():
            serializable[str(n_max)][label] = result
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'wilson_dirac_s3_test.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(serializable, f, indent=2, default=str)
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
