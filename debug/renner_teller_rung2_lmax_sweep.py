"""m-seam rung 2S, angular push: Sturmian fixed-R (R = 2) |m|=1 bending
coefficient a2 and pi energy vs (l_max, n_max), l_max = 2, 3, 4.
Machinery = tests/test_renner_teller_rung2{,_sturmian}.py.  Appends one JSON
line per point to debug/data/renner_teller_rung2_lmax_sweep.jsonl."""
import importlib.util
import json
import os
import sys
import time

import mpmath as mp
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    'rs', os.path.join(ROOT, 'tests', 'test_renner_teller_rung2_sturmian.py'))
rs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rs)
mp.mp.dps = 40
R = 2
OUT = os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung2_lmax_sweep.jsonl')
PI_X, PI_Y = ('x', 'l3c1'), ('y', 'l3s1')

points = [(2, n) for n in (3, 4, 5, 6)] + [(3, n) for n in (4, 5, 6)] + [(4, n) for n in (5, 6)]
if len(sys.argv) > 1:
    points = [tuple(int(v) for v in a.split(',')) for a in sys.argv[1:]]
for lmax, nm in points:
    t0 = time.time()
    B, T = rs.t_series(nm, l_max=lmax)
    fns = [sp.lambdify(rs.Ss, M, 'mpmath') for M in T]
    t1 = time.time()
    p0 = rs.fixed_r_p0(B, fns, R, 0, rs.IN, PI_X, mp.mpf('1.2'))
    def dE(sv):
        px = rs.fixed_r_p0(B, fns, R, sv, rs.IN, PI_X, p0)
        py = rs.fixed_r_p0(B, fns, R, sv, rs.OUT, PI_Y, p0)
        return -(px**2 - py**2) / 2
    h = [mp.mpf('2e-3'), mp.mpf('1e-3'), mp.mpf('5e-4')]
    q = [dE(v) / v**2 for v in h]
    r1, r2 = (4 * q[1] - q[0]) / 3, (4 * q[2] - q[1]) / 3
    row = {'l_max': lmax, 'n_max': nm, 'N': len(B), 'p0_pi': float(p0),
           'E_pi': float(-p0**2 / 2), 'a2': float(r2), 'a2_richardson_spread': float(abs(r1 - r2)),
           'build_s': round(t1 - t0, 1), 'solve_s': round(time.time() - t1, 1)}
    print(json.dumps(row))
    sys.stdout.flush()
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')
