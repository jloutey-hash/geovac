"""Diagnostic: is the large-R relay growth (~1/R, sign flip) the hydrogenic
2s-2p degeneracy (accidental SO(4)) amplifying the dipole relay?  Break the
l-degeneracy with a quantum-defect-like shift on the s orbitals and watch the
relay.  Machinery = tests/test_renner_teller_rung2.py."""
import importlib.util
import json
import os

import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    'rt2', os.path.join(ROOT, 'tests', 'test_renner_teller_rung2.py'))
rt2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rt2)
mp.mp.dps = 50

B, H = rt2.hamiltonian_series(3)
s_idx = [i for i, b in enumerate(B) if b[1] == 0]
rows = []
for shift in (0.0, -0.02, -0.1):
    for Rv in (2, 8, 16):
        Hn = rt2.numeric(H, Rv)
        for i in s_idx:
            Hn[0][i, i] += mp.mpf(shift)
        d_in, r_in = rt2.second_order(B, Hn, rt2.IN, ('x',))
        d_out, r_out = rt2.second_order(B, Hn, rt2.OUT, ('y',))
        row = {'s_shift': shift, 'R': Rv, 'direct': float(d_in - d_out),
               'relay': float(r_in - r_out), 'a2': float(d_in + r_in - d_out - r_out)}
        rows.append(row)
        print('shift=% .2f R=%-3d direct=% .4e relay=% .4e a2=% .4e'
              % (shift, Rv, row['direct'], row['relay'], row['a2']))
with open(os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung2_defect.json'), 'w') as f:
    json.dump(rows, f, indent=1)
