"""Scan the AO operator budget on one coefficient file (record docs/predictions/2026-10-28_angle_operator.md).
Candidates (R, Q), R = 1..Rmax, Q = 1..Qmax, (R+1)(Q+1) <= n_fit - 10, in order of (R+1)(Q+1), then R.
For each: smallest singular value, gap s2/smin, worst held-out residual (last `held` rows).  FOUND = gap >= 1e15,
smin <= 10 eps, held <= 1e3 eps.  Stops at the first FOUND unless --all.
    python debug/ao_scan.py <coeffs.json> <dps> <Rmax> <Qmax> <held> <eps> <out.json> [--all]
"""
import json
import sys
import time

import mpmath as mp

import ao_guess as G

path, dps, Rmax, Qmax, held, eps, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), \
    int(sys.argv[5]), sys.argv[6], sys.argv[7]
scan_all = '--all' in sys.argv
mp.mp.dps = dps
eps = mp.mpf(eps)
F = [mp.mpf(s.strip('[').split(' +/-')[0]) for s in json.load(open(path))['F']]
L = len(F) - 1
cands = []
for R in range(1, Rmax + 1):
    for Q in range(1, Qmax + 1):
        nrow = L + 1 - R - Q
        nfit = nrow - held
        if (R + 1) * (Q + 1) <= nfit - 10:
            cands.append((R, Q, nfit))
cands.sort(key=lambda c: ((c[0] + 1) * (c[1] + 1), c[0]))
log, found = [], None
t0 = time.time()
for R, Q, nfit in cands:
    g = G.guess(F, R, Q, nfit)
    gap = g['s2'] / g['smin'] if g['smin'] else mp.inf
    ok = gap >= mp.mpf(10) ** 15 and g['smin'] <= 10 * eps and g['held_max'] <= 1000 * eps
    log.append(dict(R=R, Q=Q, smin=mp.nstr(g['smin'], 4), s2=mp.nstr(g['s2'], 4), gap=mp.nstr(gap, 4),
                    held=mp.nstr(g['held_max'], 4), found=ok))
    print(R, Q, log[-1]['smin'], log[-1]['gap'], log[-1]['held'], 'FOUND' if ok else '', round(time.time() - t0),
          flush=True)
    if ok and found is None:
        found = dict(R=R, Q=Q, keys=[list(k) for k in g['keys']], v=[mp.nstr(x, dps - 5) for x in g['v']])
        if not scan_all:
            break
json.dump(dict(file=path, log=log, found=found, seconds=time.time() - t0), open(out, 'w'), indent=0)
