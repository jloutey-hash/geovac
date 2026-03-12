"""
LiH FCI Performance Profiler
=============================
Profiles a single LiH FCI calculation at R=3.015, nmax=3.
Breaks down time by phase, measures sparsity, and identifies bottlenecks.

Output: debug/PERFORMANCE_PROFILE.md
"""
import sys
import os
import time
import tracemalloc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.lattice_index import MolecularLatticeIndex


def profile_construction(Z_A, Z_B, nmax_A, nmax_B, R, n_electrons):
    """Profile MolecularLatticeIndex construction phases."""
    print("=" * 70)
    print(f"PHASE 1: Construction (Z_A={Z_A}, Z_B={Z_B}, nmax={nmax_A}, R={R})")
    print("=" * 70)

    tracemalloc.start()
    t_total = time.perf_counter()

    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax_A, nmax_B=nmax_B,
        R=R, n_electrons=n_electrons,
        vee_method='slater_full',
        fci_method='direct',
    )

    t_construct = time.perf_counter() - t_total
    mem_construct = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print(f"\n  Construction time: {t_construct:.3f}s")
    print(f"  Peak memory: {mem_construct[1] / 1024**2:.1f} MB")
    print(f"  N_spatial: {mol._n_spatial}")
    print(f"  N_spinorb: {mol.n_sp}")
    print(f"  N_SD: {mol.n_sd:,}")
    print(f"  N_electrons: {mol.n_electrons}")

    return mol, t_construct, mem_construct[1]


def profile_sparsity(mol):
    """Analyze sparsity of all key matrices."""
    print("\n" + "=" * 70)
    print("PHASE 2: Sparsity Analysis")
    print("=" * 70)

    results = {}

    # H1 spatial matrix
    H1 = np.asarray(mol._H1_spatial.todense())
    n_spatial = H1.shape[0]
    h1_nnz = np.count_nonzero(H1)
    h1_total = H1.size
    h1_diag_nnz = np.count_nonzero(np.diag(H1))
    h1_offdiag_nnz = h1_nnz - h1_diag_nnz
    results['h1'] = {
        'shape': H1.shape,
        'nnz': h1_nnz,
        'total': h1_total,
        'density': h1_nnz / h1_total,
        'diag_nnz': h1_diag_nnz,
        'offdiag_nnz': h1_offdiag_nnz,
        'mem_dense_MB': H1.nbytes / 1024**2,
    }
    print(f"\n  H1_spatial ({n_spatial}x{n_spatial}):")
    print(f"    NNZ: {h1_nnz}/{h1_total} = {100*h1_nnz/h1_total:.1f}% nonzero")
    print(f"    Diagonal NNZ: {h1_diag_nnz}, Off-diagonal NNZ: {h1_offdiag_nnz}")
    print(f"    Dense memory: {H1.nbytes / 1024:.1f} KB")

    # ERI dict analysis
    eri = mol._eri
    n_eri = len(eri)
    eri_vals = np.array(list(eri.values()))
    eri_max_possible = n_spatial**4
    eri_abs = np.abs(eri_vals)

    # Threshold analysis
    thresholds = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    thresh_counts = {}
    for th in thresholds:
        thresh_counts[th] = np.sum(eri_abs > th)

    results['eri'] = {
        'n_stored': n_eri,
        'max_possible': eri_max_possible,
        'density': n_eri / eri_max_possible,
        'max_val': float(np.max(eri_abs)),
        'min_val': float(np.min(eri_abs)),
        'mean_val': float(np.mean(eri_abs)),
        'thresh_counts': thresh_counts,
        'mem_dict_est_MB': n_eri * (4*8 + 8 + 56) / 1024**2,  # 4 ints + float + overhead
        'mem_dense_4d_MB': eri_max_possible * 8 / 1024**2,
    }

    print(f"\n  ERI (n_spatial={n_spatial}, max 4D = {n_spatial}^4 = {eri_max_possible:,}):")
    print(f"    Stored entries: {n_eri:,} ({100*n_eri/eri_max_possible:.2f}% of dense)")
    print(f"    Value range: [{float(np.min(eri_vals)):.6f}, {float(np.max(eri_vals)):.6f}]")
    print(f"    |ERI| range: [{float(np.min(eri_abs)):.2e}, {float(np.max(eri_abs)):.2e}]")
    print(f"    Dict memory (est): {results['eri']['mem_dict_est_MB']:.2f} MB")
    print(f"    Dense 4D memory: {results['eri']['mem_dense_4d_MB']:.2f} MB")
    print(f"    Threshold analysis:")
    for th, cnt in thresh_counts.items():
        print(f"      |ERI| > {th:.0e}: {cnt:,} entries ({100*cnt/n_eri:.1f}% of stored)")

    # Categorize ERI by type (same-atom vs cross-atom)
    n_A = mol._n_spatial_A
    n_B = mol._n_spatial_B
    same_A, same_B, cross = 0, 0, 0
    for (a, b, c, d) in eri.keys():
        a_on_A = a < n_A
        b_on_A = b < n_A
        c_on_A = c < n_A
        d_on_A = d < n_A
        if a_on_A and b_on_A and c_on_A and d_on_A:
            same_A += 1
        elif not a_on_A and not b_on_A and not c_on_A and not d_on_A:
            same_B += 1
        else:
            cross += 1
    results['eri_types'] = {'same_A': same_A, 'same_B': same_B, 'cross': cross}
    print(f"    ERI by type: same-A={same_A}, same-B={same_B}, cross-atom={cross}")

    return results


def profile_direct_ci(mol):
    """Profile the DirectCI assembly and eigsh solve separately."""
    print("\n" + "=" * 70)
    print("PHASE 3: DirectCI Assembly + Eigensolve")
    print("=" * 70)

    from geovac.direct_ci import DirectCISolver

    # Phase 3a: DirectCI precomputation
    tracemalloc.start()
    t0 = time.perf_counter()
    solver = DirectCISolver(mol)
    t_precomp = time.perf_counter() - t0
    mem_after_precomp = tracemalloc.get_traced_memory()

    # Measure dense array sizes
    h1_dense_mem = solver._H1_dense.nbytes / 1024**2
    eri_4d_mem = solver._eri_4d.nbytes / 1024**2
    print(f"\n  DirectCI precompute: {t_precomp:.3f}s")
    print(f"    H1_dense: {solver._H1_dense.shape} = {h1_dense_mem:.3f} MB")
    print(f"    ERI_4d: {solver._eri_4d.shape} = {eri_4d_mem:.2f} MB")

    # ERI 4D sparsity
    eri_4d = solver._eri_4d
    eri_4d_nnz = np.count_nonzero(eri_4d)
    eri_4d_total = eri_4d.size
    print(f"    ERI_4d NNZ: {eri_4d_nnz:,}/{eri_4d_total:,} = {100*eri_4d_nnz/eri_4d_total:.2f}%")

    # Phase 3b: Hamiltonian assembly
    t0 = time.perf_counter()
    H = solver.assemble_hamiltonian()
    t_assembly = time.perf_counter() - t0
    mem_after_assembly = tracemalloc.get_traced_memory()

    # CI Hamiltonian sparsity
    ci_nnz = H.nnz
    ci_total = H.shape[0] * H.shape[1]
    ci_mem = (H.data.nbytes + H.indices.nbytes + H.indptr.nbytes) / 1024**2

    print(f"\n  Hamiltonian assembly: {t_assembly:.3f}s")
    print(f"    Shape: {H.shape}")
    print(f"    NNZ: {ci_nnz:,} / {ci_total:,} = {100*ci_nnz/ci_total:.4f}%")
    print(f"    CSR memory: {ci_mem:.1f} MB")
    print(f"    Avg NNZ per row: {ci_nnz / H.shape[0]:.1f}")

    # Phase 3c: Eigensolve
    t0 = time.perf_counter()
    k = 1
    rng = np.random.RandomState(42)
    v0 = rng.randn(H.shape[0])
    eigvals, eigvecs = eigsh_profiled(H, k, v0)
    t_eigsh = time.perf_counter() - t0
    mem_after_eigsh = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print(f"\n  Eigensolve (eigsh): {t_eigsh:.3f}s")
    print(f"    E_elec = {eigvals[0]:.6f} Ha")
    print(f"    Peak memory: {mem_after_eigsh[1] / 1024**2:.1f} MB")

    results = {
        't_precomp': t_precomp,
        't_assembly': t_assembly,
        't_eigsh': t_eigsh,
        'ci_nnz': ci_nnz,
        'ci_total': ci_total,
        'ci_density': ci_nnz / ci_total,
        'ci_mem_MB': ci_mem,
        'eri_4d_nnz': eri_4d_nnz,
        'eri_4d_total': eri_4d_total,
        'eri_4d_density': eri_4d_nnz / eri_4d_total,
        'h1_dense_MB': h1_dense_mem,
        'eri_4d_MB': eri_4d_mem,
        'eigval': float(eigvals[0]),
        'mem_peak_MB': mem_after_eigsh[1] / 1024**2,
    }

    return results


def eigsh_profiled(H, k, v0):
    """Wrapper around eigsh that counts matvec operations."""
    from scipy.sparse.linalg import eigsh, LinearOperator

    matvec_count = [0]
    matvec_time = [0.0]

    def counting_matvec(x):
        t0 = time.perf_counter()
        result = H @ x
        matvec_time[0] += time.perf_counter() - t0
        matvec_count[0] += 1
        return result

    A = LinearOperator(H.shape, matvec=counting_matvec, dtype=H.dtype)
    eigvals, eigvecs = eigsh(A, k=k, which="SA", v0=v0)

    print(f"    Matvec count: {matvec_count[0]}")
    print(f"    Total matvec time: {matvec_time[0]:.3f}s")
    if matvec_count[0] > 0:
        print(f"    Avg matvec time: {1000*matvec_time[0]/matvec_count[0]:.2f} ms")

    return eigvals, eigvecs


def profile_nmax2(Z_A, Z_B, R, n_electrons):
    """Quick profile at nmax=2 for scaling comparison."""
    print("\n" + "=" * 70)
    print("PHASE 4: nmax=2 Scaling Comparison")
    print("=" * 70)

    t0 = time.perf_counter()
    mol2 = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=2, nmax_B=2,
        R=R, n_electrons=n_electrons,
        vee_method='slater_full',
        fci_method='direct',
    )
    t_construct = time.perf_counter() - t0

    from geovac.direct_ci import DirectCISolver

    t0 = time.perf_counter()
    solver2 = DirectCISolver(mol2)
    t_precomp = time.perf_counter() - t0

    t0 = time.perf_counter()
    H2 = solver2.assemble_hamiltonian()
    t_assembly = time.perf_counter() - t0

    t0 = time.perf_counter()
    rng = np.random.RandomState(42)
    v0 = rng.randn(H2.shape[0])
    from scipy.sparse.linalg import eigsh
    eigvals, _ = eigsh(H2, k=1, which="SA", v0=v0)
    t_eigsh = time.perf_counter() - t0

    results = {
        'n_spatial': mol2._n_spatial,
        'n_sp': mol2.n_sp,
        'n_sd': mol2.n_sd,
        't_construct': t_construct,
        't_precomp': t_precomp,
        't_assembly': t_assembly,
        't_eigsh': t_eigsh,
        't_total': t_construct + t_precomp + t_assembly + t_eigsh,
        'ci_nnz': H2.nnz,
        'eigval': float(eigvals[0]),
    }

    print(f"\n  nmax=2: {mol2._n_spatial} spatial, {mol2.n_sp} spinorb, {mol2.n_sd:,} SDs")
    print(f"  Construction: {t_construct:.3f}s")
    print(f"  DirectCI precomp: {t_precomp:.3f}s")
    print(f"  Assembly: {t_assembly:.3f}s")
    print(f"  Eigensolve: {t_eigsh:.3f}s")
    print(f"  Total: {results['t_total']:.3f}s")
    print(f"  CI NNZ: {H2.nnz:,}")
    print(f"  E_elec = {eigvals[0]:.6f} Ha")

    return results


def main():
    R = 3.015
    Z_A, Z_B = 3, 1
    nmax = 3
    n_electrons = 4

    print("LiH FCI Performance Profile")
    print("=" * 70)
    print(f"System: LiH, R={R}, Z_A={Z_A}, Z_B={Z_B}, nmax={nmax}, N_el={n_electrons}")
    print(f"Expected: 28 spatial, 56 spinorb, C(56,4)=367,290 SDs")
    print()

    # Phase 1: Construction
    mol, t_construct, mem_construct = profile_construction(
        Z_A, Z_B, nmax, nmax, R, n_electrons
    )

    # Phase 2: Sparsity
    sparsity = profile_sparsity(mol)

    # Phase 3: DirectCI
    ci_results = profile_direct_ci(mol)

    # Phase 4: nmax=2 comparison
    nmax2_results = profile_nmax2(Z_A, Z_B, R, n_electrons)

    # Summary
    t_total_nmax3 = t_construct + ci_results['t_precomp'] + ci_results['t_assembly'] + ci_results['t_eigsh']

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  nmax=3 Total: {t_total_nmax3:.1f}s")
    print(f"    Construction:     {t_construct:.3f}s ({100*t_construct/t_total_nmax3:.1f}%)")
    print(f"    DirectCI precomp: {ci_results['t_precomp']:.3f}s ({100*ci_results['t_precomp']/t_total_nmax3:.1f}%)")
    print(f"    H assembly:       {ci_results['t_assembly']:.3f}s ({100*ci_results['t_assembly']/t_total_nmax3:.1f}%)")
    print(f"    Eigensolve:       {ci_results['t_eigsh']:.3f}s ({100*ci_results['t_eigsh']/t_total_nmax3:.1f}%)")
    print(f"\n  nmax=2 Total: {nmax2_results['t_total']:.1f}s")
    print(f"  Scaling factor (nmax3/nmax2):")
    if nmax2_results['t_assembly'] > 0.001:
        print(f"    Assembly: {ci_results['t_assembly']/nmax2_results['t_assembly']:.1f}x")
    print(f"    SDs: {mol.n_sd/nmax2_results['n_sd']:.1f}x ({nmax2_results['n_sd']:,} -> {mol.n_sd:,})")
    print(f"    CI NNZ: {ci_results['ci_nnz']/nmax2_results['ci_nnz']:.1f}x")

    print(f"\n  Memory breakdown (nmax=3):")
    print(f"    H1_dense: {ci_results['h1_dense_MB']:.3f} MB")
    print(f"    ERI_4d:   {ci_results['eri_4d_MB']:.2f} MB")
    print(f"    CI matrix: {ci_results['ci_mem_MB']:.1f} MB")
    print(f"    Peak tracked: {ci_results['mem_peak_MB']:.1f} MB")

    print(f"\n  Sparsity:")
    print(f"    H1: {100*sparsity['h1']['density']:.1f}% dense")
    print(f"    ERI dict: {sparsity['eri']['n_stored']:,} entries ({100*sparsity['eri']['density']:.2f}% of {sparsity['eri']['max_possible']:,})")
    print(f"    ERI_4d: {100*ci_results['eri_4d_density']:.2f}% dense")
    print(f"    CI Hamiltonian: {100*ci_results['ci_density']:.4f}% dense (avg {ci_results['ci_nnz']/mol.n_sd:.1f} nnz/row)")

    # Write results to file for the report
    with open(os.path.join(os.path.dirname(__file__), 'data', 'lih_fci_profile.txt'), 'w') as f:
        f.write("LiH FCI Performance Profile\n")
        f.write(f"Date: 2026-03-12\n")
        f.write(f"System: LiH R={R} nmax={nmax}\n\n")

        f.write("=== TIMING (seconds) ===\n")
        f.write(f"nmax3_construction = {t_construct:.3f}\n")
        f.write(f"nmax3_directci_precomp = {ci_results['t_precomp']:.3f}\n")
        f.write(f"nmax3_assembly = {ci_results['t_assembly']:.3f}\n")
        f.write(f"nmax3_eigsh = {ci_results['t_eigsh']:.3f}\n")
        f.write(f"nmax3_total = {t_total_nmax3:.3f}\n")
        f.write(f"nmax2_construction = {nmax2_results['t_construct']:.3f}\n")
        f.write(f"nmax2_assembly = {nmax2_results['t_assembly']:.3f}\n")
        f.write(f"nmax2_eigsh = {nmax2_results['t_eigsh']:.3f}\n")
        f.write(f"nmax2_total = {nmax2_results['t_total']:.3f}\n\n")

        f.write("=== SIZES ===\n")
        f.write(f"nmax3_n_spatial = {mol._n_spatial}\n")
        f.write(f"nmax3_n_sp = {mol.n_sp}\n")
        f.write(f"nmax3_n_sd = {mol.n_sd}\n")
        f.write(f"nmax2_n_spatial = {nmax2_results['n_spatial']}\n")
        f.write(f"nmax2_n_sp = {nmax2_results['n_sp']}\n")
        f.write(f"nmax2_n_sd = {nmax2_results['n_sd']}\n\n")

        f.write("=== SPARSITY ===\n")
        f.write(f"h1_nnz = {sparsity['h1']['nnz']}\n")
        f.write(f"h1_total = {sparsity['h1']['total']}\n")
        f.write(f"h1_density = {sparsity['h1']['density']:.4f}\n")
        f.write(f"eri_dict_entries = {sparsity['eri']['n_stored']}\n")
        f.write(f"eri_max_possible = {sparsity['eri']['max_possible']}\n")
        f.write(f"eri_dict_density = {sparsity['eri']['density']:.6f}\n")
        f.write(f"eri_4d_nnz = {ci_results['eri_4d_nnz']}\n")
        f.write(f"eri_4d_total = {ci_results['eri_4d_total']}\n")
        f.write(f"eri_4d_density = {ci_results['eri_4d_density']:.6f}\n")
        f.write(f"ci_nnz = {ci_results['ci_nnz']}\n")
        f.write(f"ci_total = {ci_results['ci_total']}\n")
        f.write(f"ci_density = {ci_results['ci_density']:.8f}\n\n")

        f.write("=== MEMORY (MB) ===\n")
        f.write(f"h1_dense = {ci_results['h1_dense_MB']:.4f}\n")
        f.write(f"eri_4d = {ci_results['eri_4d_MB']:.4f}\n")
        f.write(f"ci_csr = {ci_results['ci_mem_MB']:.1f}\n")
        f.write(f"peak_tracked = {ci_results['mem_peak_MB']:.1f}\n\n")

        f.write("=== ERI TYPES ===\n")
        f.write(f"same_A = {sparsity['eri_types']['same_A']}\n")
        f.write(f"same_B = {sparsity['eri_types']['same_B']}\n")
        f.write(f"cross = {sparsity['eri_types']['cross']}\n\n")

        f.write("=== ERI THRESHOLD ANALYSIS ===\n")
        for th, cnt in sparsity['eri']['thresh_counts'].items():
            f.write(f"|ERI|>{th:.0e} = {cnt}\n")

        f.write(f"\n=== ENERGIES (Ha) ===\n")
        f.write(f"nmax3_E_elec = {ci_results['eigval']:.6f}\n")
        f.write(f"nmax2_E_elec = {nmax2_results['eigval']:.6f}\n")
        f.write(f"V_NN = {mol.V_NN:.6f}\n")

    print(f"\n  Results saved to debug/data/lih_fci_profile.txt")


if __name__ == '__main__':
    main()
