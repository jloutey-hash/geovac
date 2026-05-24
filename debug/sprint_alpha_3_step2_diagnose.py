"""Diagnostic: trace why multi-zeta substitution gives bit-identical FCI."""

import numpy as np

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec

R = 3.5
spec = nah_spec(R=R, max_n=2)
n_e = sum(b.n_electrons for b in spec.blocks)

print(f"NaH at R={R}, max_n=2, n_e={n_e}")
print(f"Spec blocks: {len(spec.blocks)}")
for b in spec.blocks:
    print(f"  {b.label}: Z_center={b.Z_center}, n_electrons={b.n_electrons}, "
          f"Z_nuc_center={b.Z_nuc_center}, n_val_offset={b.n_val_offset}, "
          f"max_n={b.max_n}")

# Run A: baseline
r_a = build_balanced_hamiltonian(spec, R=R, multi_zeta_basis=False, verbose=False)
# Run B: multi-zeta
r_b = build_balanced_hamiltonian(spec, R=R, multi_zeta_basis=True, verbose=False)

print(f"\nh1 diff max: {np.max(np.abs(r_a['h1'] - r_b['h1'])):.6e}")
print(f"h1_cross_vne diff max: {np.max(np.abs(r_a['h1_cross_vne'] - r_b['h1_cross_vne'])):.6e}")
print(f"eri diff max: {np.max(np.abs(r_a['eri'] - r_b['eri'])):.6e}")
print(f"nuclear_repulsion diff: {r_a['nuclear_repulsion'] - r_b['nuclear_repulsion']:.6e}")

# The Na center sub-block lives at orbital indices 0..4 (5 orbitals,
# block_n=1 and block_n=2). Slot 0 is (block_n=1, l=0, m=0) = the Na 3s.
print(f"\nh1[0,0] (Na 3s on-site): baseline={r_a['h1'][0,0]:.6f}, multizeta={r_b['h1'][0,0]:.6f}, diff={r_b['h1'][0,0]-r_a['h1'][0,0]:+.6e}")
print(f"h1_cross_vne[0,0]: baseline={r_a['h1_cross_vne'][0,0]:.6f}, multizeta={r_b['h1_cross_vne'][0,0]:.6f}, diff={r_b['h1_cross_vne'][0,0]-r_a['h1_cross_vne'][0,0]:+.6e}")
print(f"h1_cross_vne[1,1] (block_n=2,l=0): baseline={r_a['h1_cross_vne'][1,1]:.6f}, multizeta={r_b['h1_cross_vne'][1,1]:.6f}")
print(f"h1_cross_vne[5,5] (NaH_bond_partner block_n=1,l=0 = H 1s): baseline={r_a['h1_cross_vne'][5,5]:.6f}, multizeta={r_b['h1_cross_vne'][5,5]:.6f}")
print(f"h1_no_pk[0,0] (kinetic+nucleus on Na 3s): {r_a['h1_no_pk'][0,0]:.6f}")
print(f"h1 final[0,0] = h1_no_pk + h1_cross_vne + ... : {r_a['h1'][0,0]:.6f}")

# Run the FCI separately to confirm
print("\nRunning coupled_fci...", flush=True)
fci_a = coupled_fci_energy(r_a, n_electrons=n_e, verbose=False)
fci_b = coupled_fci_energy(r_b, n_electrons=n_e, verbose=False)
print(f"FCI baseline: {fci_a['E_coupled']:.10f} Ha")
print(f"FCI multizeta: {fci_b['E_coupled']:.10f} Ha")
print(f"FCI diff: {fci_b['E_coupled']-fci_a['E_coupled']:+.6e} Ha")

# Eigenvalue comparison on raw h1
e_a, v_a = np.linalg.eigh(r_a['h1'])
e_b, v_b = np.linalg.eigh(r_b['h1'])
print(f"\nh1 eigenvalue diff (sorted): {np.max(np.abs(e_a - e_b)):.6e}")
print(f"h1 lowest eigenvalue baseline: {e_a[0]:.6f}")
print(f"h1 lowest eigenvalue multizeta: {e_b[0]:.6f}")
print(f"h1 lowest 5 eigenvalues:")
for i in range(min(5, len(e_a))):
    print(f"  i={i}: base={e_a[i]:.6f}  mz={e_b[i]:.6f}  diff={e_b[i]-e_a[i]:+.6e}")
