"""
Track BK: Investigate whether l-dependent PK produces correct s-p splitting.

Computes <n,l|V_PK|n,l> for hydrogenic orbitals at various Z_eff values
and compares to known HF orbital energy gaps.

PK parameters use Z²-scaled He-like defaults from composed_qubit.py:
  Li (Z=3): A=6.93, B=7.00 (Paper 17 Table 1)
  C  (Z=6): A=6.93*(6/3)²=27.72, B=7.00*(6/3)²=28.00
  N  (Z=7): A=6.93*(7/3)²=37.72, B=7.00*(7/3)²=38.11  (corrected: 7.00*(49/9)=38.11)
  O  (Z=8): A=49.28, B=49.78 (from _PK_HELIKE_DEFAULTS)

The composed architecture applies PK with delta_{l,0} — only l=0 orbitals
get the PK repulsive barrier. This pushes s-orbital energies UP.

In real atoms, core penetration makes s-orbitals MORE bound than p.
So the question is: does PK go the right way?
"""

import sys
sys.path.insert(0, "c:/Users/jlout/OneDrive/Desktop/Project_Geometric")

import numpy as np
from geovac.composed_qubit import _compute_pk_matrix_elements, _radial_wf_grid

# -----------------------------------------------------------------------
# PK parameters (Z²-scaled from Li²⁺)
# -----------------------------------------------------------------------
A_Li = 6.93
B_Li = 7.00

atoms = {
    'Li': {'Z': 3, 'n_core': 2, 'Z_eff_val': 1, 'n_val': 2,
           'A_pk': A_Li, 'B_pk': B_Li},
    'C':  {'Z': 6, 'n_core': 2, 'Z_eff_val': 4, 'n_val': 2,
           'A_pk': A_Li * (6/3)**2, 'B_pk': B_Li * (6/3)**2},
    'N':  {'Z': 7, 'n_core': 2, 'Z_eff_val': 5, 'n_val': 2,
           'A_pk': A_Li * (7/3)**2, 'B_pk': B_Li * (7/3)**2},
    'O':  {'Z': 8, 'n_core': 2, 'Z_eff_val': 6, 'n_val': 2,
           'A_pk': 49.28, 'B_pk': 49.78},
}

# Known HF orbital energies (Koopmans' theorem, from NIST/Clementi-Roetti)
# These are approximate — exact values depend on basis set
hf_orbital_energies = {
    'Li': {'2s': -0.196, '2p': -0.129},  # Li: 2s is more bound
    'C':  {'2s': -0.706, '2p': -0.434},  # C: 2s is more bound
    'N':  {'2s': -0.945, '2p': -0.573},  # N: 2s is more bound
    'O':  {'2s': -1.244, '2p': -0.632},  # O: 2s is more bound
}

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

print("=" * 78)
print("Track BK: PK s-p Splitting Investigation")
print("=" * 78)

results = {}

for name, params in atoms.items():
    Z_eff = params['Z_eff_val']
    A_pk = params['A_pk']
    B_pk = params['B_pk']
    n = params['n_val']

    print(f"\n{'-' * 78}")
    print(f"  {name} (Z={params['Z']}, Z_eff={Z_eff}, A_pk={A_pk:.2f}, B_pk={B_pk:.2f})")
    print(f"{'-' * 78}")

    # States: n=2 orbitals (2s, 2p_m=0, 2p_m=+1, 2p_m=-1)
    states_s = [(2, 0, 0)]  # 2s
    states_p = [(2, 1, 0)]  # 2p (m=0, representative)

    # Hydrogenic energy (l-degenerate)
    E_hydro = -Z_eff**2 / (2 * n**2)

    # PK matrix element for 2s
    h1_pk_s = _compute_pk_matrix_elements(Z_eff, states_s, A_pk, B_pk)
    pk_2s = h1_pk_s[0, 0]

    # PK matrix element for 2p (should be nonzero too, unless we apply delta_{l,0})
    h1_pk_p = _compute_pk_matrix_elements(Z_eff, states_p, A_pk, B_pk)
    pk_2p = h1_pk_p[0, 0]

    # With delta_{l,0} prescription: only l=0 gets PK
    E_2s_with_pk = E_hydro + pk_2s  # PK applied
    E_2p_with_pk = E_hydro          # PK NOT applied (delta_{l,0})

    # Without delta_{l,0}: all l get PK (for comparison)
    E_2s_all_pk = E_hydro + pk_2s
    E_2p_all_pk = E_hydro + pk_2p

    # HF reference
    hf_2s = hf_orbital_energies[name]['2s']
    hf_2p = hf_orbital_energies[name]['2p']
    hf_gap = hf_2p - hf_2s  # positive: 2p is higher (less bound)

    # GeoVac gap with delta_{l,0}
    gv_gap = E_2p_with_pk - E_2s_with_pk  # = -pk_2s (since pk_2s > 0)

    print(f"  Hydrogenic E(n=2) = {E_hydro:.4f} Ha")
    print(f"  <2s|V_PK|2s> = {pk_2s:.4f} Ha  (pushes 2s UP)")
    print(f"  <2p|V_PK|2p> = {pk_2p:.4f} Ha  (raw, before delta_{{l,0}})")
    print()
    print(f"  With delta_{{l,0}} PK:")
    print(f"    E(2s) = {E_2s_with_pk:.4f} Ha  (hydro + PK)")
    print(f"    E(2p) = {E_2p_with_pk:.4f} Ha  (hydro, no PK)")
    print(f"    Delta(2p-2s) = {gv_gap:.4f} Ha  ({'2p above 2s' if gv_gap > 0 else '2s above 2p'})")
    print()
    print(f"  HF reference (Koopmans):")
    print(f"    E(2s) = {hf_2s:.4f} Ha")
    print(f"    E(2p) = {hf_2p:.4f} Ha")
    print(f"    Delta(2p-2s) = {hf_gap:.4f} Ha  (2p above 2s = less bound)")
    print()

    # The key question: does PK produce the right SIGN of splitting?
    # In HF: 2s is more bound (lower) than 2p → Delta(2p-2s) > 0
    # PK pushes 2s UP → makes 2s LESS bound → wrong direction!
    if gv_gap < 0:
        direction = "WRONG — PK pushes 2s ABOVE 2p (opposite to real physics)"
    elif gv_gap > 0:
        direction = "RIGHT sign — 2p above 2s"
    else:
        direction = "No splitting"

    print(f"  Verdict: {direction}")
    print(f"  |GV gap| = {abs(gv_gap):.4f} Ha vs |HF gap| = {abs(hf_gap):.4f} Ha")

    results[name] = {
        'E_hydro': E_hydro,
        'pk_2s': pk_2s,
        'pk_2p': pk_2p,
        'E_2s_pk': E_2s_with_pk,
        'E_2p_pk': E_2p_with_pk,
        'gv_gap': gv_gap,
        'hf_gap': hf_gap,
    }

# -----------------------------------------------------------------------
# Also compute for n=3 states (3s, 3p, 3d) — relevant for Madelung rule
# -----------------------------------------------------------------------
print(f"\n{'=' * 78}")
print("  n=3 PK splittings (relevant for Madelung n+l ordering)")
print(f"{'=' * 78}")

for name in ['N', 'O']:
    params = atoms[name]
    Z_eff = params['Z_eff_val']
    A_pk = params['A_pk']
    B_pk = params['B_pk']

    E_hydro_3 = -Z_eff**2 / 18.0  # n=3

    states_3s = [(3, 0, 0)]
    states_3p = [(3, 1, 0)]
    states_3d = [(3, 2, 0)]

    pk_3s = _compute_pk_matrix_elements(Z_eff, states_3s, A_pk, B_pk)[0, 0]
    pk_3p = _compute_pk_matrix_elements(Z_eff, states_3p, A_pk, B_pk)[0, 0]
    pk_3d = _compute_pk_matrix_elements(Z_eff, states_3d, A_pk, B_pk)[0, 0]

    print(f"\n  {name}: n=3, Z_eff={Z_eff}")
    print(f"    E_hydro(n=3) = {E_hydro_3:.4f} Ha")
    print(f"    <3s|V_PK|3s> = {pk_3s:.6f} Ha")
    print(f"    <3p|V_PK|3p> = {pk_3p:.6f} Ha")
    print(f"    <3d|V_PK|3d> = {pk_3d:.6f} Ha")
    print(f"    With delta_{{l,0}}: E(3s)={E_hydro_3+pk_3s:.4f}, E(3p)=E(3d)={E_hydro_3:.4f}")
    print(f"    PK 3s-3p gap: {-pk_3s:.6f} Ha (PK pushes 3s up by this amount)")

# -----------------------------------------------------------------------
# Cross-n PK coupling: <2s|V_PK|3s> — does PK couple different n?
# -----------------------------------------------------------------------
print(f"\n{'=' * 78}")
print("  Cross-n PK coupling (off-diagonal)")
print(f"{'=' * 78}")

for name in ['Li', 'N']:
    params = atoms[name]
    Z_eff = params['Z_eff_val']
    A_pk = params['A_pk']
    B_pk = params['B_pk']

    states_ns = [(2, 0, 0), (3, 0, 0)]
    h1_pk = _compute_pk_matrix_elements(Z_eff, states_ns, A_pk, B_pk)
    print(f"\n  {name}: <2s|V_PK|2s>={h1_pk[0,0]:.6f}, <3s|V_PK|3s>={h1_pk[1,1]:.6f}, <2s|V_PK|3s>={h1_pk[0,1]:.6f}")

# -----------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------
print(f"\n{'=' * 78}")
print("  SUMMARY TABLE: PK s-p splitting vs HF")
print(f"{'=' * 78}")
print(f"  {'Atom':>4s}  {'Z_eff':>5s}  {'<2s|PK|2s>':>10s}  {'GV Delta(2p-2s)':>12s}  {'HF Delta(2p-2s)':>12s}  {'Sign?':>8s}")
print(f"  {'-'*4}  {'-'*5}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*8}")
for name in ['Li', 'C', 'N', 'O']:
    r = results[name]
    Z_eff = atoms[name]['Z_eff_val']
    sign = "WRONG" if r['gv_gap'] < 0 else "OK"
    print(f"  {name:>4s}  {Z_eff:>5.0f}  {r['pk_2s']:>10.4f}  {r['gv_gap']:>12.4f}  {r['hf_gap']:>12.4f}  {sign:>8s}")

print(f"\n  Key finding: PK pushes s-orbitals UP (less bound), but in real atoms")
print(f"  s-orbitals are MORE bound (below p). The PK effect has the WRONG SIGN")
print(f"  for reproducing the s-p splitting.")
print(f"\n  This is expected: PK creates a repulsive barrier to enforce core-valence")
print(f"  orthogonality. The core penetration effect that makes s < p in real atoms")
print(f"  comes from the r-dependent Z_eff(r) — s-orbitals 'see' more nuclear charge")
print(f"  near the nucleus. The constant Z_eff approximation in the composed")
print(f"  architecture misses this physics entirely.")
