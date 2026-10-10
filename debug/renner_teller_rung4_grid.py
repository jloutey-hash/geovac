"""m-seam rung 4 driver (engine + tests in tests/test_renner_teller_rung4.py).

Usage: python debug/renner_teller_rung4_grid.py check | sweep LMAX... | full LMAX
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    'rung4', os.path.join(ROOT, 'tests', 'test_renner_teller_rung4.py'))
g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(g)


def check():
    r, T = g.fedvr(g.BOUNDS, 10)
    for l, ex in ((0, -0.5), (1, -0.125)):
        H = T + np.diag(l * (l + 1) / (2 * r ** 2) - 1 / r)
        print('hydrogen l=%d  E=%.12f  (exact %.6f)' % (l, np.linalg.eigvalsh(H)[0], ex))
    for order in (8, 10, 12):
        rr, TT = g.fedvr(g.BOUNDS, order)
        H1 = g.h0_block(rr, TT, g.harmonics(1, 'c', 9, +1), 18)
        print('order %d Nr %d  E_pi(lmax=9) = %.12f' % (order, len(rr), np.linalg.eigvalsh(H1)[0]))


def sweep(lmaxes):
    path = os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung4_grid.jsonl')
    for lm in lmaxes:
        t0 = time.time()
        o = g.a2_linear(lm)
        o['seconds'] = round(time.time() - t0, 1)
        print(json.dumps({k: o[k] for k in ('lmax', 'Nr', 'E_pi', 'a2', 'seconds')}),
              json.dumps(o['pieces']))
        sys.stdout.flush()
        with open(path, 'a') as f:
            f.write(json.dumps(o) + '\n')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'check':
        check()
    elif cmd == 'full':
        print('lmax %s: Richardson %.8f spread %.1e perturbative %.8f'
              % ((sys.argv[2],) + g.full_bent(int(sys.argv[2]))))
    else:
        sweep([int(v) for v in sys.argv[2:]])
