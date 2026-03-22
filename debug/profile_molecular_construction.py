"""Profile MolecularLatticeIndex construction to find bottlenecks."""
import time
import sys
sys.path.insert(0, '.')

# Profile BH (Z=5, Z=1) at nmax=3 — moderate size
from geovac.lattice_index import MolecularLatticeIndex, compute_cross_atom_J, _form_factor_nl

print("=" * 60)
print("PROFILING: BH (Z=5, Z=1) at nmax=3")
print("=" * 60)

# Count how many unique (na, la, nb, lb) pairs exist
from geovac.lattice_index import LatticeIndex
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    li_A = LatticeIndex(n_electrons=1, max_n=3, nuclear_charge=5,
                        vee_method='slater_full', h1_method='exact',
                        fci_method='matrix')
    li_B = LatticeIndex(n_electrons=1, max_n=3, nuclear_charge=1,
                        vee_method='slater_full', h1_method='exact',
                        fci_method='matrix')

states_A = li_A.lattice.states
states_B = li_B.lattice.states
print(f"\nAtom A states ({len(states_A)}):")
for i, s in enumerate(states_A):
    print(f"  {i}: n={s[0]}, l={s[1]}, m={s[2]}")
print(f"\nAtom B states ({len(states_B)}):")
for i, s in enumerate(states_B):
    print(f"  {i}: n={s[0]}, l={s[1]}, m={s[2]}")

# Count unique (na, la, nb, lb) pairs
pairs = set()
for na, la, ma in states_A:
    for nb, lb, mb in states_B:
        pairs.add((na, la, nb, lb))
print(f"\nUnique (na,la,nb,lb) pairs for cross-atom J: {len(pairs)}")
for p in sorted(pairs):
    print(f"  {p}")

# Time individual form factor calls
print("\n--- Form factor timing ---")
import numpy as np
q_test = 1.0

for n, l in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]:
    t0 = time.perf_counter()
    for _ in range(10):
        _form_factor_nl(n, l, 5.0, q_test)
    elapsed = (time.perf_counter() - t0) / 10
    print(f"  _form_factor_nl(n={n}, l={l}): {elapsed*1000:.3f} ms/call")

# Time a single cross-atom J integral (uncached)
print("\n--- Cross-atom J timing (uncached) ---")
import os
cache_dir = os.path.join(os.path.dirname(__file__), '..', 'geovac', 'cache')

# Clear cache for this test
for f in os.listdir(cache_dir):
    if f.startswith('cross_atom_J_') and 'Z5_1' in f and 'R3.0' in f:
        os.remove(os.path.join(cache_dir, f))

for na, la, nb, lb in [(1, 0, 1, 0), (2, 1, 1, 0), (3, 2, 1, 0), (3, 2, 3, 2)]:
    # Clear specific cache
    cache_file = os.path.join(
        cache_dir,
        f"cross_atom_J_{na}{la}_{nb}{lb}_R3.0000_Z5_1.npy"
    )
    if os.path.exists(cache_file):
        os.remove(cache_file)

    t0 = time.perf_counter()
    val = compute_cross_atom_J(na, la, nb, lb, 3.0, 5.0, 1.0)
    elapsed = time.perf_counter() - t0
    print(f"  J({na}{la},{nb}{lb}): {val:.6f} Ha, {elapsed:.3f}s")

# Now time full MolecularLatticeIndex construction
print("\n--- Full MolecularLatticeIndex construction ---")
# Clear all BH cache files
for f in os.listdir(cache_dir):
    if 'Z5_1' in f and f.startswith('cross_atom_J_'):
        os.remove(os.path.join(cache_dir, f))

t0 = time.perf_counter()
mol = MolecularLatticeIndex(
    Z_A=5, Z_B=1, nmax_A=3, nmax_B=3,
    R=3.0, n_electrons=6,
    vee_method='slater_full',
    enumerate_sds=False,
)
elapsed = time.perf_counter() - t0
print(f"\nTotal construction time: {elapsed:.3f}s")

# Now with cache (second run)
print("\n--- Second run (cached) ---")
t0 = time.perf_counter()
mol2 = MolecularLatticeIndex(
    Z_A=5, Z_B=1, nmax_A=3, nmax_B=3,
    R=3.0, n_electrons=6,
    vee_method='slater_full',
    enumerate_sds=False,
)
elapsed = time.perf_counter() - t0
print(f"\nCached construction time: {elapsed:.3f}s")
