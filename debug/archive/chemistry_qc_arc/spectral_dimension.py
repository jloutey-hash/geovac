"""
Spectral Dimension of the GeoVac Graph

Computes the spectral dimension d_s(sigma) from the heat kernel return probability
on the GeoVac graph Laplacian and compares with the smooth S^3 continuum.

Key structural finding: the GeoVac graph as defined in lattice.py decomposes into
n_max disconnected l-sectors (no edges connect different l values). Each l-sector
is a 2D grid: path graph in n (length n_max - l) x path graph in m (length 2l+1).
This gives the graph a fundamentally different topology from the smooth S^3.

We compute:
1. GeoVac graph Laplacian spectrum and spectral dimension
2. Smooth S^3 continuum (eigenvalues n^2-1, degeneracy n^2)
3. S^3 truncated at n_max (same spectrum, finite cutoff)
4. Cubic lattice sanity checks
5. Comparison: graph vs continuum, with analysis of dimensional crossover

Author: GeoVac Project
Date: April 2026
"""

import numpy as np
from scipy import sparse
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


###############################################################################
# Graph construction
###############################################################################

def build_geovac_graph_laplacian(n_max: int, weighted: bool = False):
    """
    Build the GeoVac graph Laplacian L = D - A.

    Nodes: (n, l, m) with n=1..n_max, l=0..n-1, m=-l..l
    Edges (uniform weights unless weighted=True):
        L+: (n,l,m) -> (n,l,m+1) if m < l
        L-: (n,l,m) -> (n,l,m-1) if m > -l
        T+: (n,l,m) -> (n+1,l,m) if n < n_max and l < n+1
        T-: (n,l,m) -> (n-1,l,m) if n > 1

    IMPORTANT: l never changes. The graph decomposes into n_max disconnected
    l-sectors. Each l-sector is a 2D grid: path(n_max-l) x path(2l+1).
    """
    states = []
    state_index = {}
    idx = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
                state_index[(n, l, m)] = idx
                idx += 1

    N = len(states)
    A = np.zeros((N, N), dtype=np.float64)

    for n, l, m in states:
        i = state_index[(n, l, m)]
        if m < l:
            j = state_index[(n, l, m + 1)]
            w = 1.0 / (n * n) if weighted else 1.0
            A[i, j] = w
            A[j, i] = w
        if n < n_max:
            target = (n + 1, l, m)
            if target in state_index:
                j = state_index[target]
                w = 1.0 / (n * (n + 1)) if weighted else 1.0
                A[i, j] = w
                A[j, i] = w

    D = np.diag(A.sum(axis=1))
    L = D - A
    return L, states


def count_states(n_max: int) -> int:
    return sum(n**2 for n in range(1, n_max + 1))


###############################################################################
# Heat kernel and spectral dimension
###############################################################################

def compute_heat_kernel(eigenvalues: np.ndarray, sigma_values: np.ndarray) -> np.ndarray:
    """P(sigma) = (1/N) sum_i exp(-sigma * lambda_i)"""
    N = len(eigenvalues)
    exponent = -np.outer(sigma_values, eigenvalues)
    # Clip for numerical stability
    exponent = np.clip(exponent, -500, 500)
    P = np.exp(exponent).sum(axis=1) / N
    return P


def spectral_dimension(sigma_values: np.ndarray, P_values: np.ndarray) -> np.ndarray:
    """d_s(sigma) = -2 d(log P) / d(log sigma) via central differences."""
    log_sigma = np.log(sigma_values)
    log_P = np.log(np.maximum(P_values, 1e-300))
    d_log_P = np.gradient(log_P, log_sigma)
    return -2.0 * d_log_P


###############################################################################
# Exact S^3 spectrum
###############################################################################

def s3_heat_kernel_truncated(sigma_values: np.ndarray, n_max: int) -> np.ndarray:
    """
    S^3 heat kernel truncated at n_max.
    Eigenvalues: n^2 - 1, degeneracy: n^2, for n = 1, ..., n_max.
    Identical to the smooth S^3 with l_max = n_max - 1.
    """
    ns = np.arange(1, n_max + 1, dtype=np.float64)
    degs = ns**2
    eigs = ns**2 - 1
    N = degs.sum()
    P = np.zeros_like(sigma_values)
    for i, s in enumerate(sigma_values):
        P[i] = np.sum(degs * np.exp(-s * eigs)) / N
    return P


def s3_heat_kernel_continuum(sigma_values: np.ndarray, l_max: int = 500) -> np.ndarray:
    """
    S^3 heat kernel on the smooth 3-sphere with large truncation.
    Eigenvalues: l(l+2) = (n-1)(n+1) = n^2 - 1, degeneracy: (l+1)^2 = n^2.
    """
    return s3_heat_kernel_truncated(sigma_values, n_max=l_max + 1)


###############################################################################
# Cubic lattice for sanity check
###############################################################################

def cubic_lattice_heat_kernel(d: int, sigma_values: np.ndarray, L: int = 50) -> np.ndarray:
    """
    Heat kernel on d-dimensional periodic cubic lattice of side L.
    Uses product structure for efficiency.
    """
    ks = np.arange(L)
    eigenvalues_1d = 2.0 * (1.0 - np.cos(2 * np.pi * ks / L))
    P = np.ones_like(sigma_values)
    for _ in range(d):
        factor = np.array([np.mean(np.exp(-sigma * eigenvalues_1d)) for sigma in sigma_values])
        P *= factor
    return P


###############################################################################
# Walk dimension via thermal variance
###############################################################################

def compute_msd(eigvals, eigvecs, states, sigma_values):
    """
    Compute thermal variance of n (proxy for MSD).
    Var_sigma(n) = <n^2>_sigma - <n>_sigma^2
    where <O>_sigma = Tr(O * diag([e^{-sigma L}]_{jj})) / Tr(diag([e^{-sigma L}]_{jj})).
    """
    N = len(states)
    n_values = np.array([s[0] for s in states], dtype=np.float64)
    msd = np.zeros(len(sigma_values))

    for i, sigma in enumerate(sigma_values):
        # Diagonal of the heat kernel matrix
        diag_exp = np.sum(eigvecs**2 * np.exp(-sigma * eigvals)[np.newaxis, :], axis=1)
        Z = np.sum(diag_exp)
        if Z < 1e-300:
            msd[i] = 0
            continue
        mean_n = np.sum(n_values * diag_exp) / Z
        mean_n2 = np.sum(n_values**2 * diag_exp) / Z
        msd[i] = mean_n2 - mean_n**2

    return msd


###############################################################################
# Per-l-sector analysis
###############################################################################

def analyze_per_sector(n_max: int, sigma_values: np.ndarray):
    """
    Analyze each l-sector separately. Each is a 2D grid graph:
    path(n_max - l) x path(2l+1).
    """
    sector_data = []

    for l in range(n_max):
        n_range = n_max - l  # number of n-values (n = l+1, ..., n_max)
        m_range = 2 * l + 1  # number of m-values
        N_sector = n_range * m_range

        if N_sector == 0:
            continue

        # Path graph eigenvalues
        if n_range > 1:
            eigs_n = np.array([2 - 2 * np.cos(np.pi * k / n_range) for k in range(n_range)])
        else:
            eigs_n = np.array([0.0])

        if m_range > 1:
            eigs_m = np.array([2 - 2 * np.cos(np.pi * k / m_range) for k in range(m_range)])
        else:
            eigs_m = np.array([0.0])

        # 2D grid eigenvalues = tensor sum
        eigs_sector = np.add.outer(eigs_n, eigs_m).ravel()

        P_sector = compute_heat_kernel(eigs_sector, sigma_values)
        ds_sector = spectral_dimension(sigma_values, P_sector)

        sector_data.append({
            'l': l,
            'n_range': n_range,
            'm_range': m_range,
            'N': N_sector,
            'eigs': eigs_sector,
            'P': P_sector,
            'ds': ds_sector,
            'gap': eigs_sector[eigs_sector > 1e-10].min() if np.any(eigs_sector > 1e-10) else 0,
            'max_eig': eigs_sector.max(),
        })

    return sector_data


###############################################################################
# Main analysis
###############################################################################

def run_spectral_dimension_analysis():
    print("=" * 70)
    print("GeoVac Graph Spectral Dimension Analysis")
    print("=" * 70)

    n_sigma = 400
    sigma_values = np.logspace(-3, 3, n_sigma)
    n_max_values = [5, 10, 15, 20, 30]

    results = {
        'sigma_values': sigma_values.tolist(),
        'n_max_values': n_max_values,
    }

    # =====================================================================
    # 1. GeoVac graph Laplacian: diagonalize and compute d_s
    # =====================================================================
    print("\n--- 1. GeoVac Graph Laplacian ---")

    graph_ds = {}
    graph_P = {}
    graph_eigvals = {}
    graph_info = {}

    for n_max in n_max_values:
        N = count_states(n_max)
        print(f"\nn_max = {n_max}: N = {N} nodes")

        t0 = time.time()
        L, states = build_geovac_graph_laplacian(n_max)
        eigvals = np.linalg.eigvalsh(L)
        t1 = time.time()

        n_zero = np.sum(np.abs(eigvals) < 1e-10)
        fiedler = eigvals[eigvals > 1e-10].min() if np.any(eigvals > 1e-10) else 0

        # Edge counts
        n_ang = sum(1 for n, l, m in states if m < l)
        n_rad = sum(1 for n, l, m in states if n < n_max and l < n + 1)

        print(f"  Time: {t1-t0:.2f}s | Zero eigs: {n_zero} (= n_max = {n_max} l-sectors)")
        print(f"  Fiedler value: {fiedler:.6f} | Max eigenvalue: {eigvals[-1]:.4f}")
        print(f"  Edges: {n_ang} angular + {n_rad} radial = {n_ang+n_rad} total")

        P = compute_heat_kernel(eigvals, sigma_values)
        ds = spectral_dimension(sigma_values, P)

        graph_ds[n_max] = ds
        graph_P[n_max] = P
        graph_eigvals[n_max] = eigvals
        graph_info[n_max] = {
            'N': N, 'n_zero': int(n_zero), 'fiedler': float(fiedler),
            'max_eig': float(eigvals[-1]), 'n_angular': n_ang, 'n_radial': n_rad,
        }

    # =====================================================================
    # 2. Smooth S^3 continuum (identical spectrum, different truncation)
    # =====================================================================
    print("\n--- 2. S^3 Continuum Spectrum ---")

    s3_ds = {}
    s3_P = {}

    # Truncated S^3 at same n_max as the graph
    for n_max in n_max_values:
        P_s3 = s3_heat_kernel_truncated(sigma_values, n_max)
        ds_s3 = spectral_dimension(sigma_values, P_s3)
        s3_ds[n_max] = ds_s3
        s3_P[n_max] = P_s3

    # Large truncation (effectively continuum)
    P_s3_cont = s3_heat_kernel_continuum(sigma_values, l_max=500)
    ds_s3_cont = spectral_dimension(sigma_values, P_s3_cont)
    s3_ds['continuum'] = ds_s3_cont
    s3_P['continuum'] = P_s3_cont

    # Print comparison
    print(f"\n{'':>12} {'ds_graph':>10} {'ds_S3_trunc':>12} {'ds_S3_cont':>12}")
    for sigma_probe in [0.001, 0.01, 0.1, 1.0, 10.0]:
        idx = np.argmin(np.abs(sigma_values - sigma_probe))
        print(f"  sigma={sigma_probe:<6} {graph_ds[30][idx]:>10.4f} {s3_ds[30][idx]:>12.4f} {ds_s3_cont[idx]:>12.4f}")

    # =====================================================================
    # 3. Per-l-sector analysis (understand dimensional structure)
    # =====================================================================
    print("\n--- 3. Per-l-Sector Analysis ---")

    sector_data_20 = analyze_per_sector(20, sigma_values)
    print(f"  n_max=20: {len(sector_data_20)} sectors")
    for sd in sector_data_20[:5]:
        print(f"    l={sd['l']}: {sd['n_range']}x{sd['m_range']} grid, N={sd['N']}, "
              f"gap={sd['gap']:.4f}, max={sd['max_eig']:.4f}")

    # Effective dimension of each sector at sigma=1
    idx_1 = np.argmin(np.abs(sigma_values - 1.0))
    print(f"\n  d_s at sigma=1.0:")
    for sd in sector_data_20[:6]:
        print(f"    l={sd['l']}: d_s = {sd['ds'][idx_1]:.4f} "
              f"({'1D path' if sd['m_range']==1 else f'{sd['n_range']}x{sd['m_range']} grid'})")

    # =====================================================================
    # 4. Cubic lattice sanity check
    # =====================================================================
    print("\n--- 4. Cubic Lattice Sanity Check ---")

    cubic_ds = {}
    for d_dim in [1, 2, 3]:
        P_c = cubic_lattice_heat_kernel(d_dim, sigma_values, L=40)
        ds_c = spectral_dimension(sigma_values, P_c)
        cubic_ds[d_dim] = ds_c
        idx_mid = n_sigma // 2
        print(f"  d={d_dim}: d_s at sigma=1 is {ds_c[idx_mid]:.4f} (expected {d_dim})")

    # =====================================================================
    # 5. Walk dimension
    # =====================================================================
    print("\n--- 5. Walk Dimension ---")

    walk_results = {}
    for n_max in [5, 10, 15]:
        N = count_states(n_max)
        L_mat, states = build_geovac_graph_laplacian(n_max)
        eigvals, eigvecs = np.linalg.eigh(L_mat)

        sigma_walk = np.logspace(-2, 2, 100)
        msd = compute_msd(eigvals, eigvecs, states, sigma_walk)

        # Fit in the scaling regime
        valid = (msd > 1e-6) & (sigma_walk > 0.1) & (sigma_walk < 5)
        if np.sum(valid) > 5:
            coeffs = np.polyfit(np.log(sigma_walk[valid]), np.log(msd[valid]), 1)
            exp = coeffs[0]
            d_w = 2.0 / exp if exp > 0 else float('inf')
            print(f"  n_max={n_max}: MSD exponent = {exp:.4f}, d_w = {d_w:.4f}")
        else:
            exp, d_w = None, None
            print(f"  n_max={n_max}: insufficient data for fit")

        walk_results[n_max] = {
            'sigma': sigma_walk.tolist(),
            'msd': msd.tolist(),
            'exponent': float(exp) if exp else None,
            'd_w': float(d_w) if d_w and d_w < 1e6 else None,
        }

    # =====================================================================
    # 6. Topologically-weighted graph
    # =====================================================================
    print("\n--- 6. Topologically-Weighted Graph ---")

    weighted_ds = {}
    for n_max in [10, 20]:
        L_w, _ = build_geovac_graph_laplacian(n_max, weighted=True)
        eigvals_w = np.linalg.eigvalsh(L_w)
        P_w = compute_heat_kernel(eigvals_w, sigma_values)
        ds_w = spectral_dimension(sigma_values, P_w)
        weighted_ds[n_max] = ds_w
        print(f"  n_max={n_max}: eigenvalue range [{eigvals_w[0]:.8f}, {eigvals_w[-1]:.6f}]")

    # =====================================================================
    # 7. Degree distribution
    # =====================================================================
    print("\n--- 7. Degree Distribution ---")
    for n_max in [10, 20]:
        L_mat, states = build_geovac_graph_laplacian(n_max)
        degrees = np.diag(L_mat).astype(int)
        unique_degrees, counts = np.unique(degrees, return_counts=True)
        print(f"  n_max={n_max}:")
        for deg, cnt in zip(unique_degrees, counts):
            frac = cnt / len(states) * 100
            print(f"    degree {deg}: {cnt:>5} nodes ({frac:.1f}%)")

    # =====================================================================
    # Save results
    # =====================================================================
    output_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(output_dir, 'data')
    plot_dir = os.path.join(output_dir, 'plots')

    save_data = {
        'sigma_values': sigma_values.tolist(),
        'n_max_values': n_max_values,
        'graph_info': graph_info,
        'graph_ds': {str(k): v.tolist() for k, v in graph_ds.items()},
        's3_ds': {str(k): v.tolist() for k, v in s3_ds.items()},
        'walk_dimension': walk_results,
        'key_findings': {},
    }

    # Key findings
    n_max_ref = 30
    ds_ref = graph_ds[n_max_ref]
    ds_s3_ref = s3_ds[n_max_ref]

    sigma_probes = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]
    findings = {}
    for sp in sigma_probes:
        idx = np.argmin(np.abs(sigma_values - sp))
        findings[str(sp)] = {
            'graph_ds': float(ds_ref[idx]),
            's3_trunc_ds': float(ds_s3_ref[idx]),
            's3_cont_ds': float(ds_s3_cont[idx]),
        }
    save_data['key_findings'] = findings

    peak_idx = np.argmax(ds_ref)
    save_data['peak_ds'] = {
        'value': float(ds_ref[peak_idx]),
        'sigma': float(sigma_values[peak_idx]),
    }

    with open(os.path.join(data_dir, 'spectral_dimension_results.json'), 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nResults saved to debug/data/spectral_dimension_results.json")

    # =====================================================================
    # PLOTS
    # =====================================================================

    # --- Plot 1: Main result: d_s(sigma) for graph vs S^3 ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(n_max_values)))
    for i, n_max in enumerate(n_max_values):
        N = count_states(n_max)
        ax.semilogx(sigma_values, graph_ds[n_max], color=colors[i],
                     linewidth=1.5, label=f'Graph n_max={n_max} (N={N})')

    ax.semilogx(sigma_values, ds_s3_cont, 'k--', linewidth=2.5, label='Smooth S$^3$ continuum')
    ax.semilogx(sigma_values, s3_ds[30], 'r:', linewidth=2, label='S$^3$ truncated (n_max=30)')

    ax.axhline(y=3, color='red', linestyle=':', alpha=0.3, linewidth=0.8)
    ax.axhline(y=2, color='blue', linestyle=':', alpha=0.3, linewidth=0.8)
    ax.axhline(y=1, color='green', linestyle=':', alpha=0.3, linewidth=0.8)
    ax.text(1e-3, 3.1, '$d_s=3$', fontsize=9, color='red', alpha=0.5)
    ax.text(1e-3, 2.1, '$d_s=2$', fontsize=9, color='blue', alpha=0.5)
    ax.text(1e-3, 1.1, '$d_s=1$', fontsize=9, color='green', alpha=0.5)

    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', fontsize=12)
    ax.set_title('GeoVac Graph vs Smooth S$^3$', fontsize=13)
    ax.legend(fontsize=7.5, loc='upper right')
    ax.set_ylim(-0.5, 5)
    ax.grid(True, alpha=0.2)

    # --- Plot 1b: Comparison at matched truncation ---
    ax = axes[1]
    for n_max in [10, 20, 30]:
        ax.semilogx(sigma_values, graph_ds[n_max], linewidth=1.5,
                     label=f'Graph n_max={n_max}')
        ax.semilogx(sigma_values, s3_ds[n_max], '--', linewidth=1.5,
                     label=f'S$^3$ trunc n_max={n_max}')

    ax.axhline(y=3, color='red', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', fontsize=12)
    ax.set_title('Matched Truncation: Graph vs S$^3$', fontsize=13)
    ax.legend(fontsize=8)
    ax.set_ylim(-0.5, 5)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_geovac.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_geovac.png")

    # --- Plot 2: Per-sector spectral dimension ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax = axes[0]
    sector_colors = plt.cm.plasma(np.linspace(0.1, 0.9, min(8, len(sector_data_20))))
    for i, sd in enumerate(sector_data_20[:8]):
        l = sd['l']
        ax.semilogx(sigma_values, sd['ds'], color=sector_colors[i], linewidth=1.2,
                     label=f"l={l} ({sd['n_range']}x{sd['m_range']})")

    ax.semilogx(sigma_values, graph_ds[20], 'k-', linewidth=2.5, label='Full graph (n_max=20)')
    ax.axhline(y=2, color='blue', linestyle=':', alpha=0.3)
    ax.axhline(y=1, color='green', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', fontsize=12)
    ax.set_title('Per-l-Sector Spectral Dimension (n_max=20)', fontsize=13)
    ax.legend(fontsize=7.5, ncol=2)
    ax.set_ylim(-0.5, 4)
    ax.grid(True, alpha=0.2)

    # --- Sector dimension diagram ---
    ax = axes[1]
    ls = np.arange(20)
    n_ranges = [20 - l for l in ls]
    m_ranges = [2*l + 1 for l in ls]
    aspects = [n_ranges[i] / m_ranges[i] if m_ranges[i] > 0 else float('inf') for i in range(len(ls))]

    ax.bar(ls - 0.15, n_ranges, width=0.3, alpha=0.7, label='n-range (radial)', color='steelblue')
    ax.bar(ls + 0.15, m_ranges, width=0.3, alpha=0.7, label='m-range (angular)', color='coral')
    ax.set_xlabel('Angular momentum $l$', fontsize=12)
    ax.set_ylabel('Grid dimension', fontsize=12)
    ax.set_title('Sector Grid Sizes (n_max=20)', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_sectors.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_sectors.png")

    # --- Plot 3: Heat kernel comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax = axes[0]
    for i, n_max in enumerate(n_max_values):
        N = count_states(n_max)
        ax.loglog(sigma_values, graph_P[n_max], color=colors[i],
                  linewidth=1.5, label=f'Graph n_max={n_max}')

    ax.loglog(sigma_values, P_s3_cont, 'k--', linewidth=2.5, label='Smooth S$^3$')

    sigma_ref = np.logspace(-2, 0, 50)
    ax.loglog(sigma_ref, 0.3 * sigma_ref**(-1.5), 'r:', alpha=0.4, label=r'$\propto\sigma^{-3/2}$ ($d_s$=3)')
    ax.loglog(sigma_ref, 0.3 * sigma_ref**(-1.0), 'b:', alpha=0.4, label=r'$\propto\sigma^{-1}$ ($d_s$=2)')
    ax.loglog(sigma_ref, 0.3 * sigma_ref**(-0.5), 'g:', alpha=0.4, label=r'$\propto\sigma^{-1/2}$ ($d_s$=1)')

    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Return probability $P(\sigma)$', fontsize=12)
    ax.set_title('Heat Kernel Return Probability', fontsize=13)
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.2)

    # --- Ratio plot ---
    ax = axes[1]
    for n_max in [10, 20, 30]:
        ratio = graph_P[n_max] / s3_heat_kernel_truncated(sigma_values, n_max)
        ax.semilogx(sigma_values, ratio, linewidth=1.5, label=f'n_max={n_max}')

    ax.axhline(y=1, color='k', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'$P_{\rm graph}(\sigma) / P_{S^3}(\sigma)$', fontsize=12)
    ax.set_title('Heat Kernel Ratio: Graph / S$^3$ (matched truncation)', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_heat_kernel.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_heat_kernel.png")

    # --- Plot 4: Comparison with cubic lattices ---
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.semilogx(sigma_values, graph_ds[20], 'b-', linewidth=2, label='GeoVac graph (n_max=20)')
    ax.semilogx(sigma_values, s3_ds[20], 'b--', linewidth=2, label='S$^3$ trunc (n_max=20)')
    ax.semilogx(sigma_values, ds_s3_cont, 'k--', linewidth=2, label='S$^3$ continuum')

    for d_dim, ls_style, col in [(1, ':', 'green'), (2, '-.', 'orange'), (3, '-.', 'red')]:
        ax.semilogx(sigma_values, cubic_ds[d_dim], linestyle=ls_style, color=col,
                     linewidth=1.5, alpha=0.7, label=f'{d_dim}D cubic (L=40)')

    ax.axhline(y=3, color='red', linestyle=':', alpha=0.3)
    ax.axhline(y=2, color='orange', linestyle=':', alpha=0.3)
    ax.axhline(y=1, color='green', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', fontsize=12)
    ax.set_title('Spectral Dimension: GeoVac vs Cubic Lattices vs S$^3$', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_ylim(-0.5, 5)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_comparison.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_comparison.png")

    # --- Plot 5: Uniform vs weighted ---
    fig, ax = plt.subplots(figsize=(10, 6))

    for n_max in [10, 20]:
        ax.semilogx(sigma_values, graph_ds[n_max], linewidth=1.5,
                     label=f'Uniform (n_max={n_max})')
        ax.semilogx(sigma_values, weighted_ds[n_max], '--', linewidth=1.5,
                     label=f'Weighted 1/(n$_1$n$_2$) (n_max={n_max})')

    ax.semilogx(sigma_values, ds_s3_cont, 'k--', linewidth=2, label='S$^3$ continuum')
    ax.axhline(y=3, color='red', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Spectral dimension $d_s(\sigma)$', fontsize=12)
    ax.set_title('Uniform vs Topologically-Weighted Graph', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_ylim(-0.5, 5)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_weighted.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_weighted.png")

    # --- Plot 6: Eigenvalue density comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    ax = axes[0]
    for n_max in [10, 20, 30]:
        eigv = graph_eigvals[n_max]
        ax.hist(eigv, bins=80, density=True, alpha=0.4, label=f'Graph n_max={n_max}')
    ax.set_xlabel('Eigenvalue $\\lambda$', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Graph Laplacian Eigenvalue Distribution', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.2)

    ax = axes[1]
    # S^3 continuum: eigenvalues n^2-1 with degeneracy n^2
    ns = np.arange(1, 31)
    eigs_s3 = ns**2 - 1
    degs_s3 = ns**2
    ax.bar(eigs_s3, degs_s3 / degs_s3.sum(), width=3, alpha=0.7, color='black',
           label='S$^3$ spectrum ($n^2-1$, deg $n^2$)')
    ax.set_xlabel('Eigenvalue', fontsize=12)
    ax.set_ylabel('Relative weight', fontsize=12)
    ax.set_title('S$^3$ Spectrum (vs Graph)', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 200)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_eigenvalues.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_eigenvalues.png")

    # --- Plot 7: Walk dimension ---
    fig, ax = plt.subplots(figsize=(10, 6))

    for n_max in [5, 10, 15]:
        wd = walk_results[n_max]
        sigma_w = np.array(wd['sigma'])
        msd_w = np.array(wd['msd'])
        valid = msd_w > 1e-8
        if np.any(valid):
            d_w_str = f", d_w={wd['d_w']:.1f}" if wd['d_w'] else ""
            ax.loglog(sigma_w[valid], msd_w[valid], linewidth=1.5,
                      label=f'n_max={n_max}{d_w_str}')

    sigma_ref2 = np.logspace(-1, 1, 50)
    ax.loglog(sigma_ref2, 0.1 * sigma_ref2**1.0, 'k:', alpha=0.4, label=r'$\propto\sigma^1$ (d_w=2)')
    ax.loglog(sigma_ref2, 0.1 * sigma_ref2**0.5, 'r:', alpha=0.4, label=r'$\propto\sigma^{0.5}$ (d_w=4)')

    ax.set_xlabel(r'Diffusion time $\sigma$', fontsize=12)
    ax.set_ylabel(r'Variance $\langle\Delta n^2\rangle$', fontsize=12)
    ax.set_title('Walk Dimension from Thermal Variance of $n$', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectral_dimension_walk.png'), dpi=150)
    plt.close()
    print("Saved: spectral_dimension_walk.png")

    # =====================================================================
    # Summary
    # =====================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\n1. STRUCTURAL FINDING: The GeoVac graph decomposes into n_max disconnected")
    print("   l-sectors. No edges connect different l values. Each sector is a 2D grid:")
    print("   path(n_max - l) in n  x  path(2l+1) in m.")
    print(f"   n_max=30: {30} disconnected components, {count_states(30)} total nodes.")

    print("\n2. Spectral dimension d_s(sigma) at selected probe values:")
    print(f"{'sigma':>8} {'Graph(30)':>10} {'S3_trunc(30)':>13} {'S3_cont':>10} {'Ratio G/S3':>11}")
    for sp in [0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0]:
        idx = np.argmin(np.abs(sigma_values - sp))
        g = graph_ds[30][idx]
        s_t = s3_ds[30][idx]
        s_c = ds_s3_cont[idx]
        r = g / s_t if abs(s_t) > 0.01 else float('nan')
        print(f"{sp:>8.3f} {g:>10.4f} {s_t:>13.4f} {s_c:>10.4f} {r:>11.4f}")

    peak_val = ds_ref[np.argmax(ds_ref)]
    peak_sigma = sigma_values[np.argmax(ds_ref)]
    print(f"\n   Peak d_s = {peak_val:.4f} at sigma = {peak_sigma:.4f}")

    s3_peak_val = ds_s3_cont[np.argmax(ds_s3_cont)]
    s3_peak_sigma = sigma_values[np.argmax(ds_s3_cont)]
    print(f"   S^3 continuum peak d_s = {s3_peak_val:.4f} at sigma = {s3_peak_sigma:.6f}")

    print("\n3. KEY PHYSICS:")
    print("   a) The smooth S^3 has d_s -> 3 as sigma -> 0 (UV), as expected for a")
    print("      3-dimensional manifold. The finite truncation reduces this.")
    print(f"   b) The graph peaks at d_s ~ {peak_val:.2f}, well below the S^3 value of 3.")
    print("   c) The graph's disconnected l-sectors prevent diffusion from accessing")
    print("      all angular momentum channels. Each sector behaves as a ~2D grid,")
    print("      so the graph looks effectively 2-dimensional.")
    print("   d) At large sigma (IR), both graph and S^3 show d_s -> 0 due to")
    print("      finite-size effects (compact manifold / finite graph).")
    print("   e) The S^3 continuum shows DIMENSIONAL REDUCTION: d_s -> 3 at UV,")
    print("      d_s -> 0 at IR. This is just the standard finite-size effect on S^3,")
    print("      NOT the CDT-type d_s -> 2 reduction.")

    print("\n4. DIMENSIONAL REDUCTION VERDICT:")
    print("   The GeoVac graph does NOT show CDT-type dimensional reduction (d_s: 4 -> 2).")
    print("   It shows a DIFFERENT phenomenon: the disconnected l-sector structure")
    print("   caps the effective spectral dimension at ~2, which is the dimension of")
    print("   the largest individual sectors (2D grids). The peak d_s ~ 2 is a")
    print("   consequence of the graph's topology (no l-changing edges), not a")
    print("   UV/IR crossover in the quantum-gravity sense.")
    print("   The smooth S^3 itself shows d_s -> 3 (its manifold dimension), with")
    print("   no sub-dimensional UV behavior.")

    print("\n5. WALK DIMENSION:")
    for n_max in [5, 10, 15]:
        wd = walk_results[n_max]
        if wd['d_w']:
            print(f"   n_max={n_max}: d_w = {wd['d_w']:.1f} (anomalous: d_w >> 2)")
    print("   The extremely high d_w reflects subdiffusive transport on the graph,")
    print("   caused by trapping in disconnected l-sectors.")

    return results


if __name__ == '__main__':
    results = run_spectral_dimension_analysis()
