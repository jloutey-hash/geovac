#!/usr/bin/env python3
"""HeH+ 2D-variational recompute (QA group2 re-cert remediation, 2026-06-27).

Paper 15 reports 93.1% of D_e for HeH+ at l_max=4 (sigma+pi) using the ADIABATIC
solver (tab:heh_convergence caption: "The adiabatic solver is used for all
entries"), then compares it directly to H2's *variational* 94.1% as "parity."
The adiabatic solver overestimates D_e by ~11% (variational-bound violator).

This recomputes HeH+ at the SAME (R=1.46, l_max=4, sigma+pi, charge-center
origin) with the *2D variational* solver (n_coupled=-1) -- the same solver that
produced H2's 94.1% -- so the comparison is apples-to-apples.

Setup (paper_15 sec:heh_results + tab:heh_convergence):
  R = 1.46 bohr, charge-center origin z0 = R/6 = 0.2433
  E_atoms = E(He) = -2.9037 Ha   (dissociation HeH+ -> He + H+)
  D_e^exact = 0.0750 Ha (Bishop 1977), E_total^exact = -2.9787 Ha
  D_e = E_atoms - E_total ;  pct = D_e / 0.0750 * 100
"""
import time
from geovac.level4_multichannel import solve_level4_h2_multichannel

R = 1.46
E_ATOMS_HE = -2.9037     # He ground state (exact nonrel), the HeH+ dissociation limit
DE_EXACT = 0.0750        # Bishop 1977
ETOT_EXACT = -2.9787

t0 = time.time()
result = solve_level4_h2_multichannel(
    R=R, l_max=4, m_max=1, n_coupled=-1,        # 2D variational (NOT adiabatic)
    Z_A=2.0, Z_B=1.0, origin='charge_center',   # z0 = R/6 = 0.2433
    n_alpha=60, n_Re=120, verbose=True,          # paper's grid (same as H2 fixture)
)
elapsed = time.time() - t0

E_elec = result['E_elec']
E_total = result['E_total']
V_NN = 2.0 * 1.0 / R
D_e = E_ATOMS_HE - E_total
pct = D_e / DE_EXACT * 100.0

print("\n" + "=" * 64)
print("HeH+ 2D-VARIATIONAL RECOMPUTE (l_max=4, sigma+pi, charge-center)")
print("=" * 64)
print(f"  result keys     : {sorted(result.keys())}")
print(f"  N_channels      : {result.get('channels') and len(result['channels'])}")
print(f"  z0 (origin)     : {R/6:.4f} bohr (charge-center)")
print(f"  V_NN = Z_A Z_B/R: {V_NN:.6f} Ha")
print(f"  E_elec          : {E_elec:.6f} Ha")
print(f"  E_total         : {E_total:.6f} Ha   (exact {ETOT_EXACT})")
print(f"  E_atoms (He)    : {E_ATOMS_HE:.6f} Ha")
print(f"  D_e = E_atoms-E_total = {D_e:.6f} Ha   (exact {DE_EXACT})")
print(f"  D_e / D_e^exact : {pct:.1f}%   (adiabatic paper value: 93.1%)")
print(f"  variational?    : E_total {'>' if E_total > ETOT_EXACT else '<='} exact "
      f"({'OK above' if E_total > ETOT_EXACT else 'VIOLATION below'})")
print(f"  wall time       : {elapsed:.1f}s")
print("=" * 64)
