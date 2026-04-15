"""
Track CUSP-3: TC at high l_max with proper particle-number-projected FCI.

The original BX-3 TC benchmark used qubit-space diagonalization (full 2^Q
matrix including wrong-particle-number sectors), giving a false positive
where TC appeared to beat standard FCI. The v2.9.0 corrected protocol uses
particle-number-projected FCI (Slater determinants in spin-orbital basis)
and showed:
  - n_max=1: TC 5.3% vs Std 1.93% (Std wins)
  - n_max=2: TC 3.6% vs Std 1.6% (Std wins)
  - n_max=3: TC 3.4% vs Std 0.36% (Std wins by 9x)

The TC gain (or lack thereof) at high n_max was never tested. If TC's
plateau at ~3.3% is structural (cusp-tail unchanged) and Std improves
~l^-2, then Std crosses TC at small n_max and the gap grows. But if TC's
plateau is artificial (n_max=3 being basis-limited, not TC-limited),
maybe at n_max=4 or 5 TC converges further.

Run: He composed at n_max=2,3,4,5 with Std (Hermitian) and TC (non-Herm)
PNP-FCI. Report energy, error, Pauli count, ratio.

Output: debug/data/cusp3_tc_high_lmax.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Use the corrected pipeline from tc_verification.
from debug.tc_verification import (
    fci_ground_energy, count_pauli_nonidentity, _he_spec, E_EXACT_HE,
)
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.tc_integrals import build_tc_composed_hamiltonian


def main():
    print("CUSP-3: TC vs Standard, particle-number-projected FCI, He")
    print("=" * 70)
    print(f"E_exact = {E_EXACT_HE:.6f} Ha")
    print()

    rows = []
    for max_n in (2, 3, 4, 5):
        print(f"--- n_max = {max_n} ---")
        spec = _he_spec(max_n)

        # Standard
        t0 = time.time()
        res_std = build_composed_hamiltonian(spec, verbose=False)
        E_std, _, n_det, _ = fci_ground_energy(
            res_std['h1'], res_std['eri'], n_electrons=2,
            nuclear_repulsion=0.0, hermitian=True)
        t_std = time.time() - t0
        N_pauli_std = count_pauli_nonidentity(res_std['qubit_op'])
        err_std_pct = 100 * abs(E_std - E_EXACT_HE) / abs(E_EXACT_HE)
        print(f"  Std: Q={res_std['Q']}, n_det={n_det}, "
              f"N_pauli={N_pauli_std}")
        print(f"       E={E_std:.6f}  err={err_std_pct:.4f}%  ({t_std:.1f}s)")

        # TC
        t0 = time.time()
        res_tc = build_tc_composed_hamiltonian(spec, pk_in_hamiltonian=False,
                                               verbose=False)
        E_tc, max_imag, _, _ = fci_ground_energy(
            res_tc['h1'], res_tc['eri'], n_electrons=2,
            nuclear_repulsion=res_tc['nuclear_repulsion'],
            hermitian=False)
        t_tc = time.time() - t0
        N_pauli_tc = count_pauli_nonidentity(res_tc['qubit_op'])
        err_tc_pct = 100 * abs(E_tc - E_EXACT_HE) / abs(E_EXACT_HE)
        print(f"  TC:  Q={res_tc['Q']}, "
              f"N_pauli={N_pauli_tc} ({N_pauli_tc/max(N_pauli_std,1):.2f}x)")
        print(f"       E={E_tc:.6f}  err={err_tc_pct:.4f}%  "
              f"max|imag|={max_imag:.2e}  ({t_tc:.1f}s)")

        winner = 'STD' if err_std_pct < err_tc_pct else 'TC'
        margin = abs(err_std_pct - err_tc_pct)
        print(f"  Winner: {winner} (margin {margin:.4f}pp)")

        rows.append({
            'n_max': max_n, 'Q_std': res_std['Q'], 'Q_tc': res_tc['Q'],
            'n_det': n_det,
            'N_pauli_std': N_pauli_std, 'N_pauli_tc': N_pauli_tc,
            'pauli_ratio': N_pauli_tc / max(N_pauli_std, 1),
            'E_std': float(E_std), 'E_tc': float(E_tc),
            'err_std_pct': float(err_std_pct),
            'err_tc_pct': float(err_tc_pct),
            'max_imag_tc': float(max_imag),
            'winner': winner,
            'margin_pp': float(margin),
            't_std_s': t_std, 't_tc_s': t_tc,
        })

    # Convergence table
    print("\n" + "=" * 70)
    print("Convergence table:")
    print(f"  {'n_max':>5} {'err_std%':>10} {'err_tc%':>10} "
          f"{'TC/Std Pauli':>13} {'winner':>7}")
    for r in rows:
        print(f"  {r['n_max']:>5} {r['err_std_pct']:>10.4f} "
              f"{r['err_tc_pct']:>10.4f} {r['pauli_ratio']:>13.2f} "
              f"{r['winner']:>7}")

    # TC plateau check: does err_tc continue to decrease?
    if len(rows) >= 3:
        ratios = [rows[i]['err_tc_pct'] / rows[i+1]['err_tc_pct']
                  for i in range(len(rows) - 1)]
        print(f"\n  TC error ratios (n_max -> n_max+1): "
              f"{[f'{r:.3f}' for r in ratios]}")
        # If ratios -> 1.0, TC has plateaued
        if max(ratios) < 1.1 and min(ratios) > 0.95:
            print(f"  TC has clearly PLATEAUED -- not a basis convergence")
        else:
            print(f"  TC is still converging (or oscillating)")

    out = {'rows': rows, 'E_exact': E_EXACT_HE}
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'cusp3_tc_high_lmax.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
