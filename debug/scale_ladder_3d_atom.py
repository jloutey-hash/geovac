"""SL-3d ladder, atom side (record docs/predictions/2026-11-01_scale_ladder_3d.md).

On every rung of the registered scale ladder (record SL, docs/predictions/2026-10-30_scale_ladder.md) the
O s and p exponents stay at their registered values (debug/data/scale_ladder_atom.json).  Only the 3d
exponent changes: instead of being tied to the 3p, it minimizes the neutral O atom's frozen-core CCSD(T)
energy, UHF reference for the M_S = 1 triplet determinant, in the same Gaussian-fitted kit (the R3
diagnosis' one-centre probe, debug/scale_ladder_r3_diag_atom.atom).  Nothing is tuned on water.

Search: a coarse grid in k3d, then a bounded Brent search in ln k3d around the best grid point.
Every evaluation is appended to debug/data/scale_ladder_3d_atom_evals_<rung>.jsonl (a killed run loses
nothing, and a re-run reuses cached points); the optimum goes to debug/data/scale_ladder_3d_atom_<rung>.json.

    python debug/scale_ladder_3d_atom.py <rung> [<rung> ...]     rungs R2 R3 R4 R5
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('a', os.path.join(HERE, 'scale_ladder_r3_diag_atom.py'))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)           # loads scale_ladder (np.diag patch) and test_water_ccsd
sl = a.d.sl
DATA = sl.DATA
ATOM = json.load(open(os.path.join(DATA, 'scale_ladder_atom.json')))
GRID = [1.5, 2.0, 2.5, 3.0, 3.5, 4.5]


def sp_exponents(rung):
    """The registered (1s, 1s', 2s, 2p, 3p) exponents of a rung, and its old (tied) 3d."""
    ks, groups = ATOM[rung]['ks'], ATOM[rung]['groups']
    five = [float(ks[g]) for g in groups]
    return five, five[4]


def _evals(rung):
    return os.path.join(DATA, 'scale_ladder_3d_atom_evals_%s.jsonl' % rung)


def _cache(rung):
    out = {}
    if os.path.exists(_evals(rung)):
        for line in open(_evals(rung)):
            r = json.loads(line)
            out[(r['rung'], round(r['k3d'], 10))] = r
    return out


def evaluate(rung, k3d, cache):
    key = (rung, round(float(k3d), 10))
    if key in cache:
        return cache[key]
    five, _ = sp_exponents(rung)
    t0 = time.time()
    e = a.atom(tuple(five) + (float(k3d),))
    row = {'rung': rung, 'k3d': float(k3d), 'UHF': e['UHF'], 'CCSD_fc': e['CCSD_fc'],
           'CCSD(T)_fc': e['CCSD(T)_fc'], 't': round(time.time() - t0, 1)}
    with open(_evals(rung), 'a') as f:
        f.write(json.dumps(row) + '\n')
    cache[key] = row
    print('%s k3d %.6f  UHF %.6f  fcCCSD %.6f  fcCCSD(T) %.6f  %.0f s' % (
        rung, k3d, row['UHF'], row['CCSD_fc'], row['CCSD(T)_fc'], row['t']), flush=True)
    return row


def optimize(rung, key='CCSD(T)_fc'):
    cache = _cache(rung)
    five, k_old = sp_exponents(rung)
    old = evaluate(rung, k_old, cache)                     # the registered (tied) 3d, for the atom change
    grid = [evaluate(rung, k, cache) for k in GRID]
    i = int(np.argmin([r[key] for r in grid]))
    if i in (0, len(GRID) - 1):
        raise RuntimeError('%s: grid minimum at the edge (k3d = %.2f); widen GRID' % (rung, GRID[i]))
    lo, hi = math.log(GRID[i - 1]), math.log(GRID[i + 1])
    res = minimize_scalar(lambda x: evaluate(rung, math.exp(x), cache)[key], bounds=(lo, hi),
                          method='bounded', options=dict(xatol=1e-3))
    best = evaluate(rung, math.exp(res.x), cache)
    return {'s_p_exponents': five, 'k3d_old': k_old, 'k3d_new': best['k3d'], 'criterion': key,
            'old': old, 'new': best,
            'atom_change_mHa': {k: 1e3 * (best[k] - old[k]) for k in ('UHF', 'CCSD_fc', 'CCSD(T)_fc')},
            'n_evals': sum(1 for kk in cache if kk[0] == rung)}


if __name__ == '__main__':
    for rung in sys.argv[1:]:
        r = optimize(rung)
        out = os.path.join(DATA, 'scale_ladder_3d_atom_%s.json' % rung)
        with open(out + '.tmp', 'w') as f:
            json.dump(r, f, indent=1)
        os.replace(out + '.tmp', out)
        print('%s: k3d %.4f -> %.4f, atom fc-CCSD(T) change %.3f mHa' % (
            rung, r['k3d_old'], r['k3d_new'], r['atom_change_mHa']['CCSD(T)_fc']), flush=True)
