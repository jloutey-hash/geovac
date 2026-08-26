"""
Validation leg 2: real balanced-LiH integrals.

  n_max=2 : DirectCI4e(faithful=True) vs live coupled_fci_energy  (must be <=1e-9 Ha)
            DirectCI4e(faithful=False) = the sign-corrected physics
  n_max=3 : same comparison if the library solver is affordable (it is ~2.3 h/pt,
            so by default only the new solver runs and we report its numbers).

Also audits the 8-fold permutational symmetry of the balanced `eri` tensor,
which the closed-form same-spin block assumes.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec
from davidson_ci import DirectCI4e

R_TRUE = 3.015


def build(R: float, max_n: int, n_grid_vne: int = 8000, L_max: int = 4):
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=n_grid_vne, L_max=L_max,
                                     screened_cross_center=False, verbose=False)
    return ham, n_e, time.perf_counter() - t0


def eri_symmetry_audit(eri):
    e = eri
    return {
        'pq_swap  (pq|rs)=(qp|rs)': float(np.abs(e - e.transpose(1, 0, 2, 3)).max()),
        'rs_swap  (pq|rs)=(pq|sr)': float(np.abs(e - e.transpose(0, 1, 3, 2)).max()),
        'bra-ket  (pq|rs)=(rs|pq)': float(np.abs(e - e.transpose(2, 3, 0, 1)).max()),
        'max|eri|': float(np.abs(e).max()),
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--max_n', type=int, default=2)
    ap.add_argument('--lib', action='store_true', help='also run the library coupled_fci_energy')
    ap.add_argument('--R', type=float, default=R_TRUE)
    args = ap.parse_args()

    print("=" * 92)
    print(f"Balanced LiH  R={args.R}  max_n={args.max_n}")
    print("=" * 92)
    ham, n_e, t_build = build(args.R, args.max_n)
    M = ham['M']
    print(f"  build_balanced_hamiltonian: {t_build:.1f}s   M={M}  n_e={n_e}  "
          f"eri {ham['eri'].nbytes/1e6:.0f} MB  E_nuc field={ham['nuclear_repulsion']:+.6f}",
          flush=True)
    sym = eri_symmetry_audit(ham['eri'])
    print("  ERI permutational symmetry audit:")
    for k, v in sym.items():
        print(f"    {k:28s} {v:.3e}")

    res = {}
    for faithful in (True, False):
        tag = 'faithful(lib-sign)' if faithful else 'corrected'
        t0 = time.perf_counter()
        ci = DirectCI4e(ham['h1'], ham['eri'], ham['nuclear_repulsion'],
                        faithful=faithful, verbose=True)
        out = ci.ground_state(tol=1e-8, verbose=False, max_sub=10)
        out['t_total'] = time.perf_counter() - t0
        res[tag] = out
        print(f"  [{tag:19s}] E = {out['E']:+.12f} Ha   |r|={out['residual']:.2e}  "
              f"iters={out['n_iter']}  n_sigma={out['n_sigma']}  "
              f"t_sigma_avg={out['t_sigma_avg']:.2f}s  wall={out['t_total']:.1f}s", flush=True)
        del ci

    if args.lib:
        from geovac.coupled_composition import coupled_fci_energy
        t0 = time.perf_counter()
        lib = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
        t_lib = time.perf_counter() - t0
        d = res['faithful(lib-sign)']['E'] - lib['E_coupled']
        print(f"\n  LIBRARY coupled_fci_energy: E = {lib['E_coupled']:+.12f} Ha  "
              f"n_det={lib['n_det']:,}  wall={t_lib:.1f}s")
        print(f"  DELTA (faithful - library) = {d:+.3e} Ha   "
              f"{'PASS (<=1e-9)' if abs(d) <= 1e-9 else 'FAIL'}")
        print(f"  sign-bug energy shift (corrected - faithful) = "
              f"{res['corrected']['E'] - res['faithful(lib-sign)']['E']:+.9f} Ha")
        res['library'] = {'E': float(lib['E_coupled']), 'n_det': int(lib['n_det']),
                          'wall_s': t_lib, 'delta_vs_faithful': float(d)}

    os.makedirs('debug/data', exist_ok=True)
    fn = f'debug/data/davidson_validate_lih_n{args.max_n}.json'
    json.dump({'R': args.R, 'max_n': args.max_n, 'M': M, 'eri_sym': sym,
               'results': {k: {kk: vv for kk, vv in v.items() if kk != 'history'}
                           for k, v in res.items()}},
              open(fn, 'w'), indent=2)
    print(f"\n[saved] {fn}")
