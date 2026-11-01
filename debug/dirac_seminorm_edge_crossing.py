"""Is the Dirac-seminorm item on the true CH window the SPECTRAL-WINDOW-COMPRESSION edge-crossing term?

Question (walls register, notes 2026-10-30/31): the cluster has two walls, the gauge defect and the order-one
defect, both the one sum  sum_{k not in W} (Dirac-gap weight) x (matrix elements crossing the window edge).
The third candidate is "convergence in the Dirac seminorm on the true window".

What this driver measures, at n_max = 1, 2, 3 (true window |lambda| <= n_max + 1/2, complete Dirac eigenspaces):

  1. The compressed Dirac seminorm's defect.  L_W(a_W) = ||[D_W, a_W]||.  Since P commutes with D,
     [D_W, a_W] = P[D, a]P exactly: the defect P[D,a]P - [D_W,a_W] is identically zero, for generators AND for
     products (band-1 functions).  A single commutator with diagonal D has no intermediate index, so no
     sum over k not in W can occur.  Contrast: the order-one defect (double commutator) is nonzero and equals the
     edge-crossing sum (reproduced here as a positive control of the setup).
  2. Control (non-tautology): a window P' = U P U^dag that does NOT commute with D (U mixes one window state
     with one outside state) gives a nonzero linear defect.  So the zero in (1) is the spectral-projection
     property, not an identity of any compression.
  3. The direction of the norm comparison.  ||P[D,a]P|| <= ||[D,a]P|| (the continuum column norm, exact here
     because a raises the level by at most one): compression is a contraction in the Dirac seminorm.  The
     discarded part is the crossing block Q[D,a]P with entries (lambda_k - lambda_i) a_ki, i in W, k not in W:
     the same matrix elements with the Dirac-gap weight, but FIRST order, a norm loss, not a broken identity.
  4. The kernel condition (the only measured Dirac-seminorm fact on the true window) holds: rerun of
     tests/test_p38_true_ch_window.py::test_dirac_commutator_kernel_is_the_scalars is done separately by pytest.

Run with OMP/OPENBLAS/MKL_NUM_THREADS=1.
"""
import importlib.util
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    'gb', os.path.join(ROOT, 'tests', 'test_paper32_window_gauge_breaking.py'))
gb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gb)

PAIRS = [(-1, 1), (1, 1), (1, -1), (-1, -1)]


def comm(lam, X):
    return lam[:, None] * X - X * lam[None, :]


def order_one_edge(lam, a, b, inW):
    w, out = np.nonzero(inW)[0], np.nonzero(~inW)[0]
    li, lk = lam[w], lam[out]
    t1 = -((li[:, None] - lk[None, :]) * a[np.ix_(w, out)]) @ b[np.ix_(out, w)]
    t2 = b[np.ix_(w, out)] @ ((lk[:, None] - li[None, :]) * a[np.ix_(out, w)])
    return t1 + t2


def run(n_max):
    S = gb.Sphere(2 * n_max + 4)            # two levels of margin: products exact on window columns
    lam = S.lam
    gens = [S.multiplier(M, Mp).toarray() for M, Mp in PAIRS]
    prods = [gens[i] @ gens[j] for i in range(4) for j in range(i, 4)]
    funcs = gens + prods
    inW = np.abs(lam) <= n_max + 0.5
    w = np.nonzero(inW)[0]
    out = {'n_max': n_max, 'dim_window': int(inW.sum())}

    # 1. linear defect (Dirac seminorm) and order-one defect (cluster wall) on the true window
    lin, o1, o1edge = 0.0, 0.0, 0.0
    for a in funcs:
        Da = comm(lam, a)
        aw = a[np.ix_(w, w)]
        lin = max(lin, np.abs(Da[np.ix_(w, w)] - comm(lam[w], aw)).max())
    for a in gens:
        for b in gens:
            aw, bw = a[np.ix_(w, w)], b[np.ix_(w, w)]
            Dw = comm(lam[w], aw)
            C = Dw @ bw - bw @ Dw
            o1 = max(o1, np.abs(C).max())
            o1edge = max(o1edge, np.abs(C - order_one_edge(lam, a, b, inW)).max())
    out['linear_defect_max'] = float(lin)
    out['order_one_defect_max'] = float(o1)
    out['order_one_minus_edge_sum_max'] = float(o1edge)

    # 2. control: a window that does not commute with D
    n = len(lam)
    Pm = np.diag(inW.astype(float))
    i0 = w[np.argmax(np.abs(lam[w]))]                 # top window state
    k0 = np.nonzero(~inW)[0][np.argmin(np.abs(lam[~inW]))]  # lowest outside state
    th = 0.3
    U = np.eye(n)
    U[i0, i0] = U[k0, k0] = np.cos(th)
    U[i0, k0], U[k0, i0] = -np.sin(th), np.sin(th)
    Pp = U @ Pm @ U.T
    Dm = np.diag(lam)
    ctrl = 0.0
    for a in gens:
        lhs = Pp @ (Dm @ a - a @ Dm) @ Pp
        DP, aP = Pp @ Dm @ Pp, Pp @ a @ Pp
        ctrl = max(ctrl, np.abs(lhs - (DP @ aP - aP @ DP)).max())
    out['control_rotated_window_linear_defect_max'] = float(ctrl)
    out['control_DP_commutator_norm'] = float(np.linalg.norm(Dm @ Pp - Pp @ Dm, 2))

    # 3. norm comparison: compressed vs continuum column norm; crossing block
    ratios, crossing = [], []
    for a in gens:
        Da = comm(lam, a)
        full_cols = np.linalg.norm(Da[:, w], 2)                 # ||[D,a] P||, exact
        comp = np.linalg.norm(Da[np.ix_(w, w)], 2)              # ||[D_W, a_W]||
        cross = np.linalg.norm(Da[np.ix_(np.nonzero(~inW)[0], w)], 2)   # ||Q [D,a] P||
        ratios.append(comp / full_cols)
        crossing.append(cross / full_cols)
    out['compressed_over_continuum_min'] = float(min(ratios))
    out['compressed_over_continuum_max'] = float(max(ratios))
    out['crossing_block_over_continuum_max'] = float(max(crossing))
    return out


def main():
    res = [run(n) for n in (1, 2, 3)]
    for r in res:
        print(json.dumps(r))
    ok = all(r['linear_defect_max'] < 1e-12 and r['order_one_defect_max'] > 0.1
             and r['order_one_minus_edge_sum_max'] < 1e-12
             and r['control_rotated_window_linear_defect_max'] > 1e-2
             and r['compressed_over_continuum_max'] <= 1 + 1e-12 for r in res)
    print('ALL CHECKS', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
