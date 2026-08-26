"""
Balanced-LiH PES sweep with the matrix-free Davidson CI.

Usage:
  python -u debug/davidson_pes.py --max_n 3 --grid banked3
  python -u debug/davidson_pes.py --max_n 4 --grid stencil --budget 5400

Grids
  banked3 : the 5 R points of the cached n_max=3 curve in
            debug/sprint_pk_amplification_lih.py (validation grid)
  stencil : the 5-point h=0.05 stencil about R_true=3.015 used by
            debug/sprint_abc_tilt_sensitivity.py (tilt/curvature grid)
  own     : a wider grid that brackets each curve's OWN minimum (omega_e grid)
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
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004 * 1.007825) / (7.016004 + 1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65 * CM1_TO_HA) ** 2         # 0.065893 Ha/bohr^2
W_E_TRUE = 1405.65

GRIDS = {
    'banked3': [2.9, 3.015, 3.1, 3.3, 3.5],
    'stencil': [2.915, 2.965, 3.015, 3.065, 3.115],
    'own':     [2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5],
    'wide':    [2.9, 3.015, 3.1, 3.2, 3.3, 3.4, 3.5, 3.7],
    # decider grid: uniform h=0.1 from R_true, 6 pts -- resolves BOTH the
    # tilt/curvature at R_true AND the curve's own minimum (R_eq ~ 3.2-3.35)
    'decider': [3.015, 2.915, 3.115, 3.215, 3.315, 3.415],
    'one':     [3.015],
}


def run_point(R: float, max_n: int, faithful: bool, tol: float,
              budget: float | None, verbose_dav: bool) -> dict:
    t0 = time.perf_counter()
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=8000, L_max=4,
                                     screened_cross_center=False, verbose=False)
    t_build = time.perf_counter() - t0
    print(f"  R={R:.4f}  build {t_build:.1f}s (M={ham['M']}, "
          f"eri {ham['eri'].nbytes/1e6:.0f} MB)", flush=True)
    out = {}
    for f in ([faithful] if faithful is not None else [True, False]):
        t1 = time.perf_counter()
        ci = DirectCI4e(ham['h1'], ham['eri'], ham['nuclear_repulsion'],
                        faithful=f, verbose=True)
        r = ci.ground_state(tol=tol, verbose=verbose_dav, max_sub=10,
                            time_budget_s=budget)
        tag = 'faithful' if f else 'corrected'
        out[tag] = {'E': r['E'], 'residual': r['residual'], 'n_iter': r['n_iter'],
                    'n_sigma': r['n_sigma'], 't_sigma_avg': r['t_sigma_avg'],
                    'wall_s': r['wall_s'], 'converged': r['converged'],
                    'setup_s': time.perf_counter() - t1 - r['wall_s']}
        print(f"    [{tag:9s}] E={r['E']:+.12f}  |r|={r['residual']:.2e}  "
              f"iters={r['n_iter']}  sigma_avg={r['t_sigma_avg']:.2f}s  "
              f"dav_wall={r['wall_s']:.0f}s  conv={r['converged']}", flush=True)
        del ci
    out['R'] = R
    out['t_build_s'] = t_build
    out['M'] = int(ham['M'])
    out['nuclear_repulsion'] = float(ham['nuclear_repulsion'])
    del ham
    return out


def fit_shape(Rs, Es, R0=R_TRUE, order=3):
    """Tilt / curvature at R0 and the own-minimum omega_e, from a polynomial fit."""
    Rs = np.asarray(Rs, float); Es = np.asarray(Es, float)
    c = np.polyfit(Rs - R0, Es, order)
    p = np.poly1d(c)
    tilt = float(p.deriv(1)(0.0))
    curv = float(p.deriv(2)(0.0))
    # own minimum
    roots = p.deriv(1).roots
    real = [float(r.real) for r in roots
            if abs(r.imag) < 1e-9 and (Rs - R0).min() - 0.35 <= r.real <= (Rs - R0).max() + 0.35]
    own = {}
    if real:
        # pick the minimum with positive curvature closest to the data centre
        cands = [r for r in real if p.deriv(2)(r) > 0]
        if cands:
            x = min(cands, key=lambda r: abs(r))
            k = float(p.deriv(2)(x))
            own = {'R_eq': x + R0, 'curv_at_min': k,
                   'curv_over_ktrue': k / K_TRUE,
                   'omega_e_cm1': float(np.sqrt(k / MU_ME) / CM1_TO_HA),
                   'R_eq_err_pct': (x + R0 - R_TRUE) / R_TRUE * 100.0}
            own['omega_e_err_pct'] = (own['omega_e_cm1'] - W_E_TRUE) / W_E_TRUE * 100.0
    return {'tilt_at_Rtrue': tilt, 'curv_at_Rtrue': curv,
            'curv_over_ktrue_at_Rtrue': curv / K_TRUE, 'order': order,
            'own_min': own}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--max_n', type=int, required=True)
    ap.add_argument('--grid', default='stencil')
    ap.add_argument('--mode', default='both', choices=['both', 'faithful', 'corrected'])
    ap.add_argument('--tol', type=float, default=1e-7)
    ap.add_argument('--budget', type=float, default=None,
                    help='per-Davidson-solve wall budget in seconds')
    ap.add_argument('--verbose_dav', action='store_true')
    ap.add_argument('--tag', default='')
    args = ap.parse_args()

    faithful = {'both': None, 'faithful': True, 'corrected': False}[args.mode]
    Rs = GRIDS[args.grid]
    print("=" * 92)
    print(f"balanced LiH PES  max_n={args.max_n}  grid={args.grid} {Rs}  mode={args.mode}")
    print("=" * 92, flush=True)
    rows = []
    t0 = time.perf_counter()
    fn = (f"debug/data/davidson_pes_n{args.max_n}_{args.grid}"
          f"{('_' + args.tag) if args.tag else ''}.json")
    os.makedirs('debug/data', exist_ok=True)
    for R in Rs:
        rows.append(run_point(R, args.max_n, faithful, args.tol,
                              args.budget, args.verbose_dav))
        json.dump({'max_n': args.max_n, 'grid': args.grid, 'R_true': R_TRUE,
                   'k_true': K_TRUE, 'rows': rows,
                   'wall_s': time.perf_counter() - t0},
                  open(fn, 'w'), indent=2)
        print(f"  [partial saved -> {fn}]  cumulative {time.perf_counter()-t0:.0f}s\n",
              flush=True)

    print("=" * 92)
    for tag in ('faithful', 'corrected'):
        if tag not in rows[0]:
            continue
        Rl = [r['R'] for r in rows]
        El = [r[tag]['E'] for r in rows]
        print(f"\n{tag}:")
        for R, E in zip(Rl, El):
            print(f"   R={R:.4f}  E={E:+.12f}")
        if len(Rl) >= 4:
            for order in (2, 3):
                if len(Rl) > order:
                    f = fit_shape(Rl, El, order=order)
                    print(f"   fit order={order}: tilt(R_true)={f['tilt_at_Rtrue']:+.6f}  "
                          f"curv(R_true)={f['curv_at_Rtrue']:+.6f} "
                          f"({f['curv_over_ktrue_at_Rtrue']:.3f}x k_true)"
                          + (f"  own-min R_eq={f['own_min']['R_eq']:.4f} "
                             f"({f['own_min']['R_eq_err_pct']:+.2f}%) "
                             f"omega_e={f['own_min']['omega_e_cm1']:.0f} cm-1 "
                             f"({f['own_min']['omega_e_err_pct']:+.1f}%)"
                             if f['own_min'] else "  [no interior min in window]"))
    json.dump({'max_n': args.max_n, 'grid': args.grid, 'R_true': R_TRUE,
               'k_true': K_TRUE, 'rows': rows, 'wall_s': time.perf_counter() - t0},
              open(fn, 'w'), indent=2)
    print(f"\n[saved] {fn}   total wall {time.perf_counter()-t0:.0f}s")
