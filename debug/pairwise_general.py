"""Driver: pairwise three-centre ERIs on water (O 1s/2s/2p, H 1s), L-convergence
and the Gaussian-fit referee sequence.  Engine: tests/test_pairwise_general.py.

Usage: python debug/pairwise_general.py [L]
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('pg', os.path.join(ROOT, 'tests', 'test_pairwise_general.py'))
g = importlib.util.module_from_spec(_s)
_s.loader.exec_module(g)

ZO1, ZO2, ZH = 7.66, 2.275, 1.0   # Slater-rule exponents (O), hydrogen 1s


def cases():
    O, H1, H2 = g.O, g.H1, g.H2
    o = lambda n, z, k: g.orbital(O, n, z, k)
    h1, h2 = g.orbital(H1, 1, ZH, 's'), g.orbital(H2, 1, ZH, 's')
    out = {
        '(O2s H1|O2s H2)': (o(2, ZO2, 's'), h1, o(2, ZO2, 's'), h2),
        '(O2pz H1|O2pz H2)': (o(2, ZO2, 'pz'), h1, o(2, ZO2, 'pz'), h2),
        '(O2px H1|O2px H2)': (o(2, ZO2, 'px'), h1, o(2, ZO2, 'px'), h2),
        '(O2px H1|O2pz H2)': (o(2, ZO2, 'px'), h1, o(2, ZO2, 'pz'), h2),
        '(O1s H1|O2pz H2)': (o(1, ZO1, 's'), h1, o(2, ZO2, 'pz'), h2),
    }
    # shared hydrogen: translate so H1 sits at the origin
    sh = lambda orb: dict(orb, center=orb['center'] - H1)
    out['(H1 O2s|H1 H2)'] = tuple(sh(x) for x in (h1, o(2, ZO2, 's'), h1, h2))
    out['(H1 O2px|H1 H2)'] = tuple(sh(x) for x in (h1, o(2, ZO2, 'px'), h1, h2))
    return out


if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rows = []
    for name, (a, b, c, d) in cases().items():
        t0 = time.time()
        part = g.eri_pairwise(a, b, c, d, L)
        tp = time.time() - t0
        ref = {ng: g.eri_gauss(a, b, c, d, ng) for ng in (8, 12, 16)}
        row = {'case': name, 'L': L, 'value': part[-1],
               'partial': {l: part[l] for l in (0, 2, 4, 6, 8, 12, 16, L) if l <= L},
               'gauss': ref, 'seconds': round(tp, 1)}
        rows.append(row)
        print('%-20s pairwise %.12f | L<=4 %+.1e L<=8 %+.1e L<=12 %+.1e L<=16 %+.1e | gauss-pair: ng8 %+.1e ng12 %+.1e ng16 %+.1e'
              % (name, part[-1], part[4] - part[-1], part[8] - part[-1], part[12] - part[-1], part[16] - part[-1],
                 ref[8] - part[-1], ref[12] - part[-1], ref[16] - part[-1]))
        sys.stdout.flush()
    with open(os.path.join(ROOT, 'debug', 'data', 'pairwise_general_water.json'), 'w') as f:
        json.dump(rows, f, indent=1)
