"""
Benchmark: Molecular Bridge Construction Performance
=====================================================

v0.9.2: Vectorized COO assembly + adaptive sparsity masking.

Profiles the molecular adjacency matrix construction across
multiple molecule types and lattice sizes. Measures:
  1. Lattice construction time (GeometricLattice init)
  2. Molecular assembly time (_build_molecular_adjacency)
  3. Bridge construction breakdown (conformal factors, pruning)
  4. Total Hamiltonian build time

Date: February 23, 2026
"""

import numpy as np
import sys
import time

sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian, GeometricLattice


def benchmark_molecule(name: str, lattices: list, connectivity: list,
                       n_trials: int = 5, **kwargs) -> dict:
    """Time molecular Hamiltonian construction over multiple trials."""
    times_lattice = []
    times_assembly = []

    for _ in range(n_trials):
        # Rebuild lattices each trial for fair timing
        t0 = time.perf_counter()
        lats = []
        for spec in lattices:
            lats.append(GeometricLattice(
                max_n=spec['max_n'],
                nucleus_position=spec['position'],
                nuclear_charge=spec['Z'],
            ))
        t1 = time.perf_counter()
        times_lattice.append(t1 - t0)

        # Molecular assembly
        t2 = time.perf_counter()
        mol = MoleculeHamiltonian(
            lattices=lats,
            connectivity=connectivity,
            bridge_amplitude=kwargs.get('bridge_amplitude', 1.0),
            bridge_decay_rate=kwargs.get('bridge_decay_rate', 1.0),
        )
        t3 = time.perf_counter()
        times_assembly.append(t3 - t2)

    lat_ms = np.median(times_lattice) * 1000
    asm_ms = np.median(times_assembly) * 1000
    total_ms = lat_ms + asm_ms

    n_states = mol.n_total_states
    nnz = mol.adjacency.nnz
    n_bridges = sum(info['n_bridges_actual'] for info in mol.bridge_info)

    return {
        'name': name,
        'n_states': n_states,
        'nnz': nnz,
        'n_bridges': n_bridges,
        'lattice_ms': lat_ms,
        'assembly_ms': asm_ms,
        'total_ms': total_ms,
    }


if __name__ == '__main__':
    print("=" * 70)
    print("BENCHMARK: Molecular Bridge Construction")
    print("=" * 70)
    print(f"\nv0.9.2: Vectorized COO + adaptive sparsity masking")
    print(f"Method: scipy.sparse COO concatenation (no lil_matrix loops)")

    results = []

    # --- H2 (homonuclear, simplest case) ---
    for max_n in [5, 10, 15]:
        n_bridges = max_n * max_n
        r = benchmark_molecule(
            f'H2 (max_n={max_n})',
            lattices=[
                {'Z': 1, 'max_n': max_n, 'position': (0, 0, 0)},
                {'Z': 1, 'max_n': max_n, 'position': (0, 0, 1.4)},
            ],
            connectivity=[(0, 1, n_bridges)],
        )
        results.append(r)

    # --- LiH (heteronuclear, asymmetric conformal factors) ---
    for max_n in [5, 10, 15]:
        n_bridges = max_n * max_n
        r = benchmark_molecule(
            f'LiH (max_n={max_n})',
            lattices=[
                {'Z': 3, 'max_n': max_n, 'position': (0, 0, 0)},
                {'Z': 1, 'max_n': max_n, 'position': (0, 0, 3.015)},
            ],
            connectivity=[(0, 1, n_bridges)],
        )
        results.append(r)

    # --- H2O-like (3 atoms, 2 bonds) ---
    for max_n in [5, 10]:
        n_bridges = max_n * max_n
        r = benchmark_molecule(
            f'H2O (max_n={max_n})',
            lattices=[
                {'Z': 8, 'max_n': max_n, 'position': (0, 0, 0)},
                {'Z': 1, 'max_n': max_n, 'position': (0, 0.757, 0.587)},
                {'Z': 1, 'max_n': max_n, 'position': (0, -0.757, 0.587)},
            ],
            connectivity=[(0, 1, n_bridges), (0, 2, n_bridges)],
        )
        results.append(r)

    # --- Print results ---
    print(f"\n  {'System':<20}  {'States':>7}  {'NNZ':>7}  {'Bridges':>7}  "
          f"{'Lattice':>9}  {'Assembly':>9}  {'Total':>9}")
    print(f"  {'-'*20}  {'-'*7}  {'-'*7}  {'-'*7}  "
          f"{'-'*9}  {'-'*9}  {'-'*9}")

    for r in results:
        print(f"  {r['name']:<20}  {r['n_states']:>7}  {r['nnz']:>7}  "
              f"{r['n_bridges']:>7}  {r['lattice_ms']:>8.1f}ms  "
              f"{r['assembly_ms']:>8.1f}ms  {r['total_ms']:>8.1f}ms")

    # Verify O(N) scaling
    h2_results = [r for r in results if r['name'].startswith('H2 ')]
    if len(h2_results) >= 2:
        n1, t1 = h2_results[0]['n_states'], h2_results[0]['assembly_ms']
        n2, t2 = h2_results[-1]['n_states'], h2_results[-1]['assembly_ms']
        if t1 > 0 and n1 > 0:
            scaling = np.log(t2 / t1) / np.log(n2 / n1) if t2 > t1 else 0.0
            print(f"\n  H2 scaling exponent: O(N^{scaling:.2f})")

    print(f"\n{'='*70}")
