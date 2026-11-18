"""op52 diagnostic (no record), part 3: exact Legendre coefficients of f split by constant group, evaluated with the
constants' numerical values, and the K6 decay instrument beta_eff on each group.

Groups (A1/A2 alike): trio = D1-endpoint exponentials e^{-zB D1}{1, e^{-mu(D2-D1)}, e^{-mu(D2+D1)}};
pair = D2-endpoint exponentials e^{+-zC D1 - (zB+zC) D2}; rest = the other 16 constants.
Hypothesis H: the pair carries all of f's x* branch (G_L(pair)/G_L(f) -> 1, beta -> -4); the trio's constant-weighted
sum carries no x* branch even though each trio component's operator has 3/2.

    python debug/op52_groups.py <case> [Lmax] [dps]   -> debug/data/op52_groups_<case>.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X           # noqa: E402
import na_exact as E             # noqa: E402

DATA = HERE / 'data'


def mpq(v):
    return mp.mpf(int(v.p)) / int(v.q)


def main(case, Lmax=200, dps=1500):
    D1, D2, zB, zC, zA = X.params(case)
    mu = 2 * zA
    t0 = time.time()
    K = E.Field(None)
    pc = X.pieces(case, Lmax, K)
    f = X.combine(pc, zA, K)
    t_split = round(time.time() - t0, 1)
    mp.mp.dps = dps
    keys = set(k for q in f for k in q)
    consts = {k: E.value(k, mp) for k in keys}
    one = lambda e: (Fraction(e), '1', None)                     # noqa: E731
    trio = [one(-zB * D1), one(-zB * D1 - mu * (D2 - D1)), one(-zB * D1 - mu * (D2 + D1))]
    pair = [one(zC * D1 - (zB + zC) * D2), one(-zC * D1 - (zB + zC) * D2)]
    groups = {'f': None, 'trio': trio, 'pair': pair, 'rest': [k for k in keys if k not in trio + pair]}
    groups.update({f'trio{i}': [k] for i, k in enumerate(trio)})
    groups.update({f'pair{i}': [k] for i, k in enumerate(pair)})
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    xsm = mp.mpf(xs.numerator) / xs.denominator
    rho = 1 / (xsm + mp.sqrt(xsm * xsm - 1))
    Ls = [20, 40, 80, 120, 160, Lmax - 20]
    res = dict(case=case, Lmax=Lmax, dps=dps, t_split=t_split, groups={})
    G = {}
    for name, ks in groups.items():
        F = [sum((mpq(v) * consts[k] for k, v in q.items() if ks is None or k in ks), mp.mpf(0)) for q in f]
        G[name] = [F[L] / rho ** L for L in range(Lmax + 1)]
    for name, g in G.items():
        beta = {L: float(mp.log(abs(g[L + 1] / g[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))) for L in Ls}
        res['groups'][name] = dict(
            beta=beta,
            G=[mp.nstr(g[L], 8) for L in Ls],
            ratio_to_f=[mp.nstr(g[L] / G['f'][L], 10) for L in Ls],
            log10_absq=[float(mp.log10(abs(g[L] * rho ** L))) if g[L] != 0 else None for L in Ls])
        print(name, json.dumps(res['groups'][name]), flush=True)
    res['time'] = round(time.time() - t0, 1)
    (DATA / f'op52_groups_{case}.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], *(int(x) for x in a[1:]))
