"""Driver: water RHF / CCSD / CCSD(T) in Gaussian (cc-pVDZ) and Slater bases.
Engine + tests: tests/test_water_ccsd.py.  Memo: debug/water_barrier_memo.md.

Usage: python debug/water_ccsd.py ccpvdz | ccpvtz | slater_min   (ccpvtz: frozen-core only)
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
OUT = os.path.join(ROOT, 'debug', 'data', 'water_ccsd.jsonl')


def record(row):
    print(json.dumps(row)); sys.stdout.flush()
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


if __name__ == '__main__':
    nuc = w.water_nuclei()
    en = w.vnn(nuc)
    kind = sys.argv[1]
    t0 = time.time()
    if kind in ('ccpvdz', 'ccpvtz'):
        fns, T = w.gaussian_basis(os.path.join(ROOT, 'debug', 'data', 'basis', kind[:2] + '-' + kind[2:] + '_H_O.json'), nuc)
        S, h, g = w.ao_integrals(fns, nuc, T)
        nbf = T.shape[1]
    else:
        fns = w.slater_minimal(nuc, 12)
        S, h, g = w.ao_integrals(fns, nuc)
        nbf = len(fns)
    t1 = time.time()
    print('integrals %.0fs nbf %d' % (t1 - t0, nbf)); sys.stdout.flush()
    fc = 1 if kind == 'ccpvtz' else 0
    r = {k: float(v) for k, v in w.run(S, h, g, 10, en, frozen_core=fc).items()}
    record({'basis': kind, 'nbf': nbf, 'frozen_core': fc, **r, 'E_nuc': en,
            't_integrals': round(t1 - t0, 1), 't_solve': round(time.time() - t1, 1)})
