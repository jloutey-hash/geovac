"""m-seam rung 2 measurements: direct vs relay share of the |m| = 1 bending
coefficient a2 (E_pi_x - E_pi_y = a2 s^2 + ...), across R and basis size,
plus the |m| = 2 s^4 coefficient.  Machinery = tests/test_renner_teller_rung2.py.
Writes debug/data/renner_teller_rung2_sweep.json."""
import importlib.util
import json
import os
import sys

import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    'rt2', os.path.join(ROOT, 'tests', 'test_renner_teller_rung2.py'))
rt2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rt2)

mp.mp.dps = 50
out = {'R_sweep': [], 'basis_sweep': []}


def m1(B, Hn):
    d_in, r_in = rt2.second_order(B, Hn, rt2.IN, ('x',))
    d_out, r_out = rt2.second_order(B, Hn, rt2.OUT, ('y',))
    return float(d_in - d_out), float(r_in - r_out)


def m2_s4(B, Hn):
    e0a, _ = rt2.pair_state(B, Hn, rt2.IN, ('x2-y2',))
    e0b, _ = rt2.pair_state(B, Hn, rt2.OUT, ('xy',))
    def split(sv):
        return rt2.tracked_level(B, Hn, sv, rt2.IN, e0a) - rt2.tracked_level(B, Hn, sv, rt2.OUT, e0b)
    h = mp.mpf('1e-3')
    q1, q2 = split(h) / h**4, split(h / 2) / (h / 2) ** 4
    return float((4 * q2 - q1) / 3)


B3, H3 = rt2.hamiltonian_series(3)
for Rv in (1, 2, 4, 8, 16):
    Hn = rt2.numeric(H3, Rv)
    d, rl = m1(B3, Hn)
    pc = -1.2 * 30 / Rv**3 * 2 / 2  # rung-1 point-charge direct term, 2p alone
    row = {'R': Rv, 'direct': d, 'relay': rl, 'a2': d + rl,
           'point_charge_direct_2p': -6 / 5 * 30 / Rv**3, 'delta_s4': m2_s4(B3, Hn)}
    out['R_sweep'].append(row)
    print('R=%-3s direct=% .6e relay=% .6e a2=% .6e  pc(2p)=% .3e  delta_s4=% .4e'
          % (Rv, d, rl, d + rl, row['point_charge_direct_2p'], row['delta_s4']))
    sys.stdout.flush()

for nm in (2, 3, 4, 5):
    B, H = rt2.hamiltonian_series(nm)
    Hn = rt2.numeric(H, 2)
    d, rl = m1(B, Hn)
    out['basis_sweep'].append({'n_max': nm, 'l_max': 2, 'R': 2, 'N': len(B),
                               'direct': d, 'relay': rl, 'a2': d + rl})
    print('n_max=%d (N=%d)  direct=% .6e relay=% .6e a2=% .6e' % (nm, len(B), d, rl, d + rl))
    sys.stdout.flush()

with open(os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung2_sweep.json'), 'w') as f:
    json.dump(out, f, indent=1)
