"""Benchmark script for He-4 Hamiltonian (Track NF)."""
from geovac.nuclear.nuclear_hamiltonian import (
    diagonalize_he4, build_he4_hamiltonian, analyze_he4_hamiltonian,
    build_deuteron_hamiltonian, analyze_deuteron_hamiltonian,
    he4_hw_scan, _compute_coulomb_me, enumerate_sp_states, DeuteronSpec,
    compute_same_species_tbme,
)
import numpy as np

print("=" * 70)
print("He-4 Benchmark (Track NF)")
print("=" * 70)

# 1. hw scan with and without Coulomb
print("\n[1] hw scan at N_shells=2 (nuclear only)")
print(f"{'hw (MeV)':>10} {'E_gs (MeV)':>12} {'E-E_free':>12}")
hw_list = [8, 10, 12, 14, 16, 18, 20, 22, 25]
for hw in hw_list:
    r = diagonalize_he4(N_shells=2, hw=hw, include_coulomb=False)
    E_free = 4 * 1.5 * hw
    print(f"{hw:10.1f} {r['E_gs']:12.4f} {r['E_gs']-E_free:12.4f}")

print("\n[2] hw scan at N_shells=2 (with Coulomb)")
print(f"{'hw (MeV)':>10} {'E_gs (MeV)':>12} {'dE_coul':>12}")
for hw in hw_list:
    r_wo = diagonalize_he4(N_shells=2, hw=hw, include_coulomb=False)
    r_w = diagonalize_he4(N_shells=2, hw=hw, include_coulomb=True)
    print(f"{hw:10.1f} {r_w['E_gs']:12.4f} {r_w['E_gs']-r_wo['E_gs']:12.4f}")

# 2. Diagonal Coulomb ME check
print("\n[3] Diagonal Coulomb ME check")
spec = DeuteronSpec(N_shells=2, hw=15.0)
states_p, _ = enumerate_sp_states(spec)
print(f"   b = {spec.b:.4f} fm")
# <0s+ 0s- | V_coul | 0s+ 0s->
cache = {}
s0 = states_p[0]  # 0s, m_l=0, m_s=+0.5
s1 = states_p[1]  # 0s, m_l=0, m_s=-0.5
v_direct = _compute_coulomb_me(s0, s1, s0, s1, spec.b, cache)
v_exchange = _compute_coulomb_me(s0, s1, s1, s0, spec.b, cache)
print(f"   <{s0} {s1}| V |{s0} {s1}> direct    = {v_direct:.6f} MeV")
print(f"   <{s0} {s1}| V |{s1} {s0}> exchange  = {v_exchange:.6f} MeV")
print(f"   Antisymmetrized = {v_direct - v_exchange:.6f} MeV")

# 3. Pauli analysis
print("\n[4] Pauli Hamiltonian Analysis (N_shells=2, hw=15)")
data_d = build_deuteron_hamiltonian(N_shells=2, hw=15.0)
data_he4 = build_he4_hamiltonian(N_shells=2, hw=15.0, include_coulomb=True)
data_he4_nc = build_he4_hamiltonian(N_shells=2, hw=15.0, include_coulomb=False)

analysis_d = analyze_deuteron_hamiltonian(data_d['H_pauli'], data_d['Q'])
analysis_he4 = analyze_he4_hamiltonian(data_he4['H_pauli'], data_he4['Q'])
analysis_he4_nc = analyze_he4_hamiltonian(data_he4_nc['H_pauli'], data_he4_nc['Q'])

print(f"\n{'Metric':<30} {'Deuteron':>12} {'He-4 (nuc)':>12} {'He-4 (+Coul)':>14}")
print(f"{'-'*70}")
print(f"{'Qubits (Q)':<30} {data_d['Q']:>12} {data_he4['Q']:>12} {data_he4['Q']:>14}")
print(f"{'Particles':<30} {'1p+1n':>12} {'2p+2n':>12} {'2p+2n':>14}")
print(f"{'Hilbert dim':<30} {data_d['H_matrix'].shape[0]:>12} "
      f"{data_he4['H_matrix'].shape[0]:>12} {data_he4['H_matrix'].shape[0]:>14}")
print(f"{'Pauli terms (non-I)':<30} {analysis_d['n_pauli_terms']:>12} "
      f"{analysis_he4_nc['n_pauli_terms']:>12} {analysis_he4['n_pauli_terms']:>14}")
print(f"{'Pauli terms (total)':<30} {analysis_d['n_pauli_terms_total']:>12} "
      f"{analysis_he4_nc['n_pauli_terms_total']:>12} {analysis_he4['n_pauli_terms_total']:>14}")
print(f"{'1-norm (non-I)':<30} {analysis_d['one_norm_ni']:>12.3f} "
      f"{analysis_he4_nc['one_norm_ni']:>12.3f} {analysis_he4['one_norm_ni']:>14.3f}")
print(f"{'1-norm (total)':<30} {analysis_d['one_norm_total']:>12.3f} "
      f"{analysis_he4_nc['one_norm_total']:>12.3f} {analysis_he4['one_norm_total']:>14.3f}")
print(f"{'QWC groups':<30} {analysis_d['n_qwc_groups']:>12} "
      f"{analysis_he4_nc['n_qwc_groups']:>12} {analysis_he4['n_qwc_groups']:>14}")
print(f"{'Z-only terms':<30} {analysis_d['n_z_only_terms']:>12} "
      f"{analysis_he4_nc['n_z_only_terms']:>12} {analysis_he4['n_z_only_terms']:>14}")
print(f"{'Non-Z terms (XY/mixed)':<30} {analysis_d['n_xy_terms']:>12} "
      f"{analysis_he4_nc['n_xy_terms']:>12} {analysis_he4['n_xy_terms']:>14}")

# Final scaling
print("\n[5] Scaling Comparison (same Q=16)")
print(f"   Deuteron (1p+1n) Pauli: {analysis_d['n_pauli_terms']}")
print(f"   He-4 (2p+2n) Pauli: {analysis_he4['n_pauli_terms']}")
if analysis_d['n_pauli_terms'] > 0:
    ratio = analysis_he4['n_pauli_terms'] / analysis_d['n_pauli_terms']
    print(f"   Ratio: {ratio:.2f}x")

print("\n[6] Ground state energies at hw=15 MeV")
print(f"   He-4 (no Coul): {diagonalize_he4(N_shells=2, hw=15, include_coulomb=False)['E_gs']:.4f} MeV")
print(f"   He-4 (Coul):    {diagonalize_he4(N_shells=2, hw=15, include_coulomb=True)['E_gs']:.4f} MeV")
print(f"   Exp. binding:   -28.296 MeV (reference)")
print(f"   Note: N_shells=2 is too small for convergence to experimental binding.")
