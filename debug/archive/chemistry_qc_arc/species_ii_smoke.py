"""Smoke test: hand-rolled FCI must match library E_coupled before trusting eigvec."""
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec
import debug.species_ii_ordering as s

spec = nah_spec()
res = build_balanced_hamiltonian(spec, R=3.566, screened_cross_center=True,
                                 multi_zeta_basis=True, cross_block_h1=True)
ref = coupled_fci_energy(res, 2)['E_coupled']
E, vec, a_s, b_s, nb = s.fci_ground(res, 2)
print(f"library E_coupled = {ref:.8f}")
print(f"my     E_gs       = {E:.8f}")
print(f"diff              = {abs(E-ref):.2e}   {'OK' if abs(E-ref) < 1e-6 else 'MISMATCH'}")
site_spatial, labels = s.site_of_spatial(res)
print("blocks:", labels)
print("site_spatial:", site_spatial.tolist())
print("norm(vec):", float((vec**2).sum()))
S = s.site_entanglement(vec, a_s, b_s, nb, site_spatial, res['M'])
print(f"S_site at R=3.566 = {S:.5f}")
