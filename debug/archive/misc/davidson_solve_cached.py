"""
Phase 2 of the two-phase sweep: Davidson-solve the cached balanced-LiH integrals.

Usage: python -u debug/davidson_solve_cached.py --max_n 4 --grid decider
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from davidson_ci import DirectCI4e
from davidson_pes import GRIDS, fit_shape, K_TRUE, R_TRUE

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--max_n', type=int, required=True)
    ap.add_argument('--grid', default='decider')
    ap.add_argument('--indir', default='debug/data/ints')
    ap.add_argument('--tol', type=float, default=1e-7)
    ap.add_argument('--budget', type=float, default=None)
    ap.add_argument('--max_sub', type=int, default=10)
    args = ap.parse_args()

    Rs = GRIDS[args.grid]
    fn_out = f'debug/data/davidson_pes_n{args.max_n}_{args.grid}.json'
    rows = []
    if os.path.exists(fn_out):                       # resume: keep finished points
        rows = json.load(open(fn_out))['rows']
        done = {round(r['R'], 4) for r in rows
                if 'faithful' in r and 'corrected' in r}
        print(f"  [resume] {len(done)} points already solved: {sorted(done)}", flush=True)
        Rs = [R for R in Rs if round(R, 4) not in done]
    t00 = time.perf_counter()
    for R in Rs:
        fn = os.path.join(args.indir, f'lih_bal_n{args.max_n}_R{R:.4f}.npz')
        if not os.path.exists(fn):
            print(f"  [missing] {fn} -- skipped", flush=True)
            continue
        t0 = time.perf_counter()
        z = np.load(fn)
        h1 = z['h1']; eri = z['eri']
        enuc = float(z['nuclear_repulsion']); M = int(z['M'])
        t_load = time.perf_counter() - t0
        print(f"  R={R:.4f}  loaded M={M} eri {eri.nbytes/1e6:.0f} MB in {t_load:.1f}s "
              f"(build was {float(z['build_s']):.0f}s)", flush=True)
        row = {'R': R, 'M': M, 'nuclear_repulsion': enuc,
               't_build_s': float(z['build_s']), 't_load_s': t_load}
        for faithful in (True, False):
            tag = 'faithful' if faithful else 'corrected'
            t1 = time.perf_counter()
            ci = DirectCI4e(h1, eri, enuc, faithful=faithful, verbose=True)
            t_setup = time.perf_counter() - t1
            r = ci.ground_state(tol=args.tol, verbose=True, max_sub=args.max_sub,
                                time_budget_s=args.budget)
            row[tag] = {'E': r['E'], 'residual': r['residual'], 'n_iter': r['n_iter'],
                        'n_sigma': r['n_sigma'], 't_sigma_avg': r['t_sigma_avg'],
                        'wall_s': r['wall_s'], 'converged': r['converged'],
                        'setup_s': t_setup}
            print(f"    [{tag:9s}] E={r['E']:+.12f}  |r|={r['residual']:.2e}  "
                  f"iters={r['n_iter']}  sigma_avg={r['t_sigma_avg']:.2f}s  "
                  f"setup={t_setup:.0f}s  dav={r['wall_s']:.0f}s  conv={r['converged']}",
                  flush=True)
            del ci
        del h1, eri, z
        rows.append(row)
        json.dump({'max_n': args.max_n, 'grid': args.grid, 'R_true': R_TRUE,
                   'k_true': K_TRUE, 'rows': rows,
                   'wall_s': time.perf_counter() - t00},
                  open(fn_out, 'w'), indent=2)
        print(f"  [partial saved -> {fn_out}]  cumulative {time.perf_counter()-t00:.0f}s\n",
              flush=True)

    for tag in ('faithful', 'corrected'):
        pts = [(r['R'], r[tag]['E']) for r in rows if tag in r]
        if len(pts) < 4:
            continue
        print(f"\n{tag}:")
        for R, E in pts:
            print(f"   R={R:.4f}  E={E:+.12f}")
        for order in (2, 3, 4):
            if len(pts) > order:
                f = fit_shape([p[0] for p in pts], [p[1] for p in pts], order=order)
                om = f['own_min']
                print(f"   fit order={order}: tilt={f['tilt_at_Rtrue']:+.6f}  "
                      f"curv={f['curv_at_Rtrue']:+.6f} "
                      f"({f['curv_over_ktrue_at_Rtrue']:.3f}x k)"
                      + (f"  R_eq={om['R_eq']:.4f} ({om['R_eq_err_pct']:+.2f}%) "
                         f"omega_e={om['omega_e_cm1']:.0f} ({om['omega_e_err_pct']:+.1f}%)"
                         if om else "  [no interior min]"))
    print(f"\n[saved] {fn_out}   total {time.perf_counter()-t00:.0f}s")
