"""m-seam rung 2S: Sturmian (Avery posing) vs hydrogenic (Galerkin) |m|=1
bending coefficient at fixed R = 2, across basis size (l <= 2).  Also the
self-consistent p0 / energy of the pi state.  Machinery =
tests/test_renner_teller_rung2{,_sturmian}.py.  Writes
debug/data/renner_teller_rung2_sturmian_sweep.json."""
import importlib.util
import json
import os
import sys

import mpmath as mp
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, 'tests', fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rs = load('rs', 'test_renner_teller_rung2_sturmian.py')
rt2 = rs.rt2
mp.mp.dps = 40
R = 2
rows = []
for nm in (2, 3, 4, 5, 6):
    B, T = rs.t_series(nm)
    fns = [sp.lambdify(rs.Ss, M, 'mpmath') for M in T]
    p0 = rs.fixed_r_p0(B, fns, R, 0, rs.IN, ('x',), mp.mpf('0.6'))
    def dE(sv):
        px = rs.fixed_r_p0(B, fns, R, sv, rs.IN, ('x',), p0)
        py = rs.fixed_r_p0(B, fns, R, sv, rs.OUT, ('y',), p0)
        return -(px**2 - py**2) / 2
    h = [mp.mpf('2e-3'), mp.mpf('1e-3')]
    q = [dE(v) / v**2 for v in h]
    a2 = (4 * q[1] - q[0]) / 3
    # fixed-S decomposition at the self-consistent S, in p0 units
    Tn = rs.neg_numeric(T, p0 * R)
    d_in, r_in = rt2.second_order(B, Tn, rs.IN, ('x',))
    d_out, r_out = rt2.second_order(B, Tn, rs.OUT, ('y',))
    row = {'n_max': nm, 'N': len(B), 'p0_pi': float(p0), 'E_pi': float(-p0**2 / 2),
           'a2_fixedR_energy': float(a2),
           'fixedS_direct_p0': float(-(d_in - d_out)), 'fixedS_relay_p0': float(-(r_in - r_out))}
    rows.append(row)
    print(row)
    sys.stdout.flush()

# hydrogenic Galerkin pi energy at R = 2 for comparison
for nm in (2, 3, 4, 5):
    B, H = rt2.hamiltonian_series(nm)
    Hn = rt2.numeric(H, R)
    e, _ = rt2.pair_state(B, Hn, rt2.IN, ('x',))
    print('hydrogenic n_max=%d  E_pi=% .8f' % (nm, float(e)))
    rows.append({'hydrogenic_n_max': nm, 'E_pi': float(e)})

with open(os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung2_sturmian_sweep.json'), 'w') as f:
    json.dump(rows, f, indent=1)
