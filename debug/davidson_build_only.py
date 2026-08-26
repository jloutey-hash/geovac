"""
Phase 1 of the two-phase n_max=4 sweep: build the balanced-LiH integrals at ONE R
and cache them to .npz.  The integral build is single-threaded, so several R points
can be built concurrently on separate cores while the (BLAS-threaded, memory-hungry)
Davidson solves run serially afterwards in debug/davidson_solve_cached.py.

Usage:  python -u debug/davidson_build_only.py --R 3.015 --max_n 4
"""
from __future__ import annotations
import argparse
import os
import time

import numpy as np
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--R', type=float, required=True)
    ap.add_argument('--max_n', type=int, required=True)
    ap.add_argument('--outdir', default='debug/data/ints')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fn = os.path.join(args.outdir, f'lih_bal_n{args.max_n}_R{args.R:.4f}.npz')
    if os.path.exists(fn):
        print(f"[skip] {fn} exists", flush=True)
        raise SystemExit(0)
    t0 = time.perf_counter()
    spec = lih_spec(R=args.R, max_n=args.max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=args.R, n_grid_vne=8000, L_max=4,
                                     screened_cross_center=False, verbose=False)
    dt = time.perf_counter() - t0
    tmp = fn + '.tmp.npz'
    np.savez(tmp, h1=ham['h1'], eri=ham['eri'],
             nuclear_repulsion=np.array(ham['nuclear_repulsion']),
             M=np.array(ham['M']), n_e=np.array(n_e), R=np.array(args.R),
             build_s=np.array(dt))
    os.replace(tmp, fn)
    print(f"[built] R={args.R:.4f} max_n={args.max_n} M={ham['M']} "
          f"in {dt:.1f}s -> {fn} ({os.path.getsize(fn)/1e6:.0f} MB)", flush=True)
