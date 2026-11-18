"""op52 diagnostic (no record), part 2.
(a) x* exponents of the generic sum over named subsets of the x*-carrying pieces (cached MO sequences);
(b) which kernel (NA = 1/s, Y = e^{-mu s}/s, dY = d/dmu Y) feeds each x*-carrying constant (fresh engine run, L <= 12).

    python debug/op52_subsets.py <case> <p>   -> debug/data/op52_subsets_<case>_<p>.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import s1_controls as C          # noqa: E402
import aabc_exact as X           # noqa: E402
import na_exact as E             # noqa: E402
from op52_pieces import xstar_exps   # noqa: E402

DATA = HERE / 'data'


def main(case, p):
    seqs, _ = C.a_coeffs(case, p)
    D1, D2, zB, zC, zA = X.params(case)
    mu = 2 * zA
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    K = lambda e: (Fraction(e), '1', None)                     # noqa: E731
    d1 = [K(-zB * D1), K(-zB * D1 - mu * (D2 - D1)), K(-zB * D1 - mu * (D2 + D1))]     # D1-endpoint trio
    d2 = [K(zC * D1 - (zB + zC) * D2), K(-zC * D1 - (zB + zC) * D2)]                  # D2-endpoint pair
    assert all(k in seqs for k in d1 + d2), (d1, d2)
    subsets = dict(trio_D1=d1, pair_D2=d2, trio_minus_NA=d1[1:], NA_D1_plus_pair=[d1[0]] + d2,
                   Ytrio_plus_pair=d1[1:] + d2, all5=d1 + d2)
    out = dict(case=case, prime=p, xstar=str(xs), subsets={}, kernels={})
    for name, keys in subsets.items():
        g = C.generic_operator(C.generic_seq({k: seqs[k] for k in keys}, p, C.SEEDS[0]), p)
        e = xstar_exps(g['c'], p, xs)
        out['subsets'][name] = dict(keys=[str(k[0]) for k in keys], order=g['order'], degree=g['degree'],
                                    xstar=e, halves=[x for x in (e or []) if Fraction(x).denominator == 2])
        print(name, g['order'], out['subsets'][name]['halves'], flush=True)
    pc = X.pieces(case, 12, E.Field(p))
    for k in d1 + d2:
        out['kernels'][str(k[0])] = {w: any(k in q and int(q[k]) != 0 for q in pc[w][2:]) for w in pc}
    print(json.dumps(out['kernels']), flush=True)
    (DATA / f'op52_subsets_{case}_{p}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
