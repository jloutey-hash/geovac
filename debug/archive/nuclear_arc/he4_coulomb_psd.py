"""Check if Coulomb operator is PSD."""
import numpy as np
from geovac.nuclear.nuclear_hamiltonian import build_he4_hamiltonian

data_nuc = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=False)
data_coul = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=True)

V_coul = data_coul['H_matrix'] - data_nuc['H_matrix']

# Is V_coul symmetric?
sym = np.allclose(V_coul, V_coul.T)
print(f"V_coul symmetric: {sym}")

# Eigenvalues of V_coul
eigs = np.linalg.eigvalsh(V_coul)
print(f"V_coul eigenvalues: min={eigs.min():.4f}, max={eigs.max():.4f}")
print(f"Num negative eigenvalues: {np.sum(eigs < -1e-10)}")
print(f"Most negative eigenvalues: {eigs[:5]}")
print(f"Most positive eigenvalues: {eigs[-5:]}")

# Diagonal of V_coul
print(f"\nDiagonal of V_coul: min={V_coul.diagonal().min():.4f}, "
      f"max={V_coul.diagonal().max():.4f}")
print(f"Num negative diagonal: {np.sum(V_coul.diagonal() < -1e-10)}")
print(f"Sum of diagonal: {V_coul.diagonal().sum():.4f}")
