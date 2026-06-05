"""LiH PES scan via balanced_coupled + number-preserving FCI.

Tests whether the balanced_coupled path (which DOES include cross-center V_ne
via multipole — the architectural difference from composed_qubit) reproduces
Paper 17's R_eq ~ 3 bohr.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

from geovac.molecular_spec import lih_spec  # noqa: E402
from geovac.balanced_coupled import build_balanced_hamiltonian  # noqa: E402
from geovac.qubit_encoding import build_fermion_op_from_integrals  # noqa: E402
from openfermion.linalg.sparse_tools import (  # noqa: E402
    get_number_preserving_sparse_operator,
)


def lih_balanced_e_at_r(R: float, max_n: int = 2) -> dict:
    spec = lih_spec(max_n=max_n)  # spec is R-agnostic; R goes to builder
    t0 = time.time()
    result = build_balanced_hamiltonian(spec, R=R, verbose=False)
    t_build_int = time.time() - t0

    fop = build_fermion_op_from_integrals(
        result['h1'], result['eri'], result['nuclear_repulsion'],
    )
    Q = result['Q']
    M = result['M']
    n_total_electrons = sum(blk.n_electrons for blk in spec.blocks)

    t0 = time.time()
    H_sparse = get_number_preserving_sparse_operator(
        fop, num_qubits=Q,
        num_electrons=n_total_electrons,
        spin_preserving=True,
    )
    t_build_fci = time.time() - t0
    dim = H_sparse.shape[0]

    t0 = time.time()
    if dim <= 1024:
        E = float(np.linalg.eigvalsh(H_sparse.toarray())[0])
    else:
        E = float(spla.eigsh(H_sparse, k=1, which='SA', tol=1e-9)[0][0])
    t_diag = time.time() - t0

    return {
        'R': R, 'M': M, 'Q': Q, 'N': n_total_electrons,
        'dim': dim, 'E': E,
        't_build_int_s': t_build_int,
        't_build_fci_s': t_build_fci,
        't_diag_s': t_diag,
    }


def main() -> int:
    print('=== LiH balanced_coupled FCI PES scan (with cross-center V_ne) ===')
    print(f'{"R[bohr]":>8s}  {"E[Ha]":>12s}  {"dim":>8s}  '
          f'{"int[s]":>7s}  {"fci[s]":>7s}  {"diag[s]":>7s}')
    R_values = [1.5, 2.0, 2.5, 3.015, 3.5, 4.0, 4.5, 5.5, 7.0]
    rows = []
    for R in R_values:
        try:
            r = lih_balanced_e_at_r(R)
            print(f'{r["R"]:8.3f}  {r["E"]:12.6f}  {r["dim"]:8d}  '
                  f'{r["t_build_int_s"]:7.1f}  {r["t_build_fci_s"]:7.1f}  '
                  f'{r["t_diag_s"]:7.2f}')
            rows.append(r)
        except Exception as e:
            print(f'{R:8.3f}  ERROR: {e}')
            import traceback
            traceback.print_exc()

    if rows:
        Rs = np.array([r['R'] for r in rows])
        Es = np.array([r['E'] for r in rows])
        i_min = int(np.argmin(Es))
        print()
        print(f'min at R={Rs[i_min]:.3f}, E={Es[i_min]:.6f} Ha')
        print(f'expt R_eq = 3.015 bohr; paper 17 ab-initio composed = 3.21 (6.4%)')
        print(f'             paper 17 l-dep PK = 5.3% error')
        if 0 < i_min < len(Rs) - 1:
            i0, i1, i2 = i_min - 1, i_min, i_min + 1
            x = Rs[[i0, i1, i2]]; y = Es[[i0, i1, i2]]
            a, b, c = np.polyfit(x, y, 2)
            R_min_fit = -b / (2 * a)
            err = abs(R_min_fit - 3.015) / 3.015 * 100
            print(f'parabolic-fit R_eq = {R_min_fit:.4f} bohr '
                  f'({err:.2f}% from expt)')

    import json
    out_path = PROJECT_ROOT / 'debug/data/diag_lih_balanced_pes.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(rows, f, indent=2)
    print(f'wrote {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
