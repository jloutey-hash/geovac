"""m-seam rung 5 driver; engine + tests in tests/test_walsh_ah2.py.

Usage: python debug/walsh_ah2.py slopes LMAX | curve LMAX
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('walsh', os.path.join(ROOT, 'tests', 'test_walsh_ah2.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'slopes':
        lmax = int(sys.argv[2])
        for name in w.MOLECULES:
            t0 = time.time()
            o = w.slopes(w.molecule(name), lmax)
            print(name, 'lmax', lmax, json.dumps({k: {q: round(v[q], 5) for q in v} for k, v in o.items()}),
                  '%.0fs' % (time.time() - t0))
    elif cmd == 'curve':
        lmax = int(sys.argv[2])
        rows = []
        for name in w.MOLECULES:
            mol = w.molecule(name, order=8)
            for ang in (180, 165, 150, 135, 120, 104.5, 90):
                sv = np.sqrt((1 + np.cos(np.radians(ang))) / 2)
                lv = w.walsh_levels(mol, lmax, sv)
                vnn = 1 / (2 * mol['R'] * np.sqrt(1 - sv ** 2))
                row = {'mol': name, 'angle': ang, 's': float(sv), 'lmax': lmax,
                       **{k: float(v) for k, v in lv.items()}, 'V_HH': float(vnn)}
                print(json.dumps(row)); sys.stdout.flush()
                rows.append(row)
        with open(os.path.join(ROOT, 'debug', 'data', 'walsh_ah2_curve_l%d.json' % lmax), 'w') as f:
            json.dump(rows, f, indent=1)
