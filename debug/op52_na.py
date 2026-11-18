"""op52 diagnostic (no record), part 4: the one-electron integral (record NA's object, MO's control geometry N3).
Mechanism prediction, written before running: NA's D1-endpoint component e^{-zB D1} and its D2-endpoint pair
e^{+-zC D1 - (zB+zC) D2} each carry the point-charge branch (one half-integer exponent, 3/2); the pair sees the
kernel's point-charge jump at D2, so pair and D1 component differ first at the next order, giving 5/2 only in the
union (NA's {3/2, 5/2}).

    python debug/op52_na.py <p>   -> debug/data/op52_na_<p>.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import s1_controls as C          # noqa: E402
import na_exact as E             # noqa: E402
from op52_pieces import xstar_exps   # noqa: E402

DATA = HERE / 'data'


def main(p):
    t0 = time.time()
    D1, D2, zB, zC = (Fraction(v) for v in C.N3[:4])
    q = E.coefficients(*C.N3, 1000, E.Field(p))
    keys = sorted(set(k for d in q for k in d), key=str)
    seqs = {k: [int(d.get(k, flint.nmod(0, p))) for d in q] for k in keys}
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    one = lambda e: (Fraction(e), '1', None)                     # noqa: E731
    d1 = [one(-zB * D1)]
    d2 = [one(zC * D1 - (zB + zC) * D2), one(-zC * D1 - (zB + zC) * D2)]
    assert all(k in seqs for k in d1 + d2), keys
    out = dict(prime=p, xstar=str(xs), coeff_time=round(time.time() - t0, 1), subsets={})
    for name, ks in dict(D1=d1, pair0=d2[:1], pair1=d2[1:], pair=d2, D1_plus_pair=d1 + d2, all=keys).items():
        g = C.generic_operator(C.generic_seq({k: seqs[k] for k in ks}, p, C.SEEDS[0]), p)
        e = xstar_exps(g['c'], p, xs)
        out['subsets'][name] = dict(order=g['order'], xstar=e,
                                    halves=[x for x in (e or []) if Fraction(x).denominator == 2])
        print(name, g['order'], out['subsets'][name]['halves'], flush=True)
    (DATA / f'op52_na_{p}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(int(sys.argv[1]))
