"""Diagnostic for rung 3: at linear geometry and the pi state's self-consistent
S, list the y-even block eigenvalues nearest k_pi with their |m| character,
and the bend coupling <phi|dK/ds|pi> (finite difference), to see whether a
near-degenerate sigma/delta state drives the a2 swings."""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'debug'))
import renner_teller_rung3 as r3  # noqa: E402
from geovac import sturmian_general_sw as g  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

R = 2.0
rows = []
for spec in sys.argv[1:]:
    cpart, hpart = spec.split(':')
    Nc, Lc = map(int, cpart.split(','))
    Nh, Lh = map(int, hpart.split(','))
    cb, hb = g.labels(Nc, Lc), g.labels(Nh, Lh)
    def kpi(S):
        K, index = r3.kernel(0.0, S, cb, hb)
        c, _ = r3.blocks(index)
        return r3.pi_branch(K, index, c)
    S = brentq(lambda SS: SS / kpi(SS) - R, 0.6, 6.0, xtol=1e-13)
    K0, index = r3.kernel(0.0, S, cb, hb)
    h = 1e-4
    K1, _ = r3.kernel(h, S, cb, hb)
    c, _ = r3.blocks(index)
    ev, V = np.linalg.eigh(K0[np.ix_(c, c)])
    dK = (K1 - K0)[np.ix_(c, c)] / h
    kp = kpi(S)
    ip = int(np.argmin(np.abs(ev - kp)))
    mchar = [sorted({index[c[i]][1][2] for i in range(len(c)) if abs(V[i, t]) > 0.1}) for t in range(len(ev))]
    near = sorted(range(len(ev)), key=lambda t: abs(ev[t] - kp))[:6]
    print('basis', spec, 'S=%.6f k_pi=%.8f' % (S, kp))
    for t in near:
        cpl = float(V[:, t] @ dK @ V[:, ip])
        print('   k=%.8f  gap=% .2e  |m| in %s  dK/ds coupling=% .3e' % (ev[t], ev[t] - kp, mchar[t], cpl))
    rows.append({'basis': spec, 'S': S, 'k_pi': kp,
                 'neighbours': [{'k': float(ev[t]), 'gap': float(ev[t] - kp), 'm': mchar[t],
                                 'coupling': float(V[:, t] @ dK @ V[:, ip])} for t in near]})
with open(os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung3_neighbours.json'), 'w') as f:
    json.dump(rows, f, indent=1)
