"""Debug Coulomb sign in He-4."""
import numpy as np
from geovac.nuclear.nuclear_hamiltonian import (
    compute_same_species_tbme, enumerate_sp_states, DeuteronSpec,
    build_he4_hamiltonian,
)

spec = DeuteronSpec(N_shells=2, hw=10.0)
states_p, _ = enumerate_sp_states(spec)

# Compute nuclear-only and nuclear+Coulomb TBME
print("Computing pp TBME without Coulomb...")
tbme_nuc = compute_same_species_tbme(states_p, spec, include_coulomb=False)
print("Computing pp TBME with Coulomb...")
tbme_coul = compute_same_species_tbme(states_p, spec, include_coulomb=True)

print(f"\nNum nuclear MEs: {len(tbme_nuc)}")
print(f"Num Coulomb MEs: {len(tbme_coul)}")

# Compare MEs for diagonal keys (0s^2)
# States 0 = (0s, m_l=0, +1/2), 1 = (0s, m_l=0, -1/2)
diag_keys = [(0,1,0,1), (1,0,1,0), (0,1,1,0), (1,0,0,1)]
print("\n0s^2 diagonal MEs:")
print(f"{'key':<15} {'nuc':>10} {'nuc+coul':>12} {'dCoul':>12}")
total_diff = 0.0
for k in diag_keys:
    v_nuc = tbme_nuc.get(k, 0.0)
    v_coul = tbme_coul.get(k, 0.0)
    dv = v_coul - v_nuc
    total_diff += dv
    print(f"{str(k):<15} {v_nuc:>10.4f} {v_coul:>12.4f} {dv:>12.4f}")
print(f"Sum of dCoul: {total_diff:.4f} (should be 4*0.69 = 2.76 MeV for antisym ME sum)")

# Direct diagonal evaluation: expected <ψ_0s^2|V_coul|ψ_0s^2> = 0.69 MeV
# In second quantization: (1/4) * sum_{ijkl} v_{ijkl}_AS <SD|a^+_i a^+_j a_l a_k|SD>
# For SD = (0,1), the only nonzero term is (i,j,k,l) with {i,j}={0,1} and {k,l}={0,1}
# There are 4 such combinations, each contributing phase*v_as/4 = 0.69 (shown above)

# Compute Coulomb energy as first-order perturbation
print("\nSum over all MEs, nuclear vs with Coulomb:")
nuc_sum = sum(abs(v) for v in tbme_nuc.values())
coul_sum = sum(abs(v) for v in tbme_coul.values())
print(f"  |V_nuc|: {nuc_sum:.3f}")
print(f"  |V_coul|: {coul_sum:.3f}")
print(f"  diff: {coul_sum - nuc_sum:.3f}")
