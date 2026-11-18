"""op52 diagnostic (no record): full indicial exponents at x* of each constant piece's minimal operator, and of
generic sums over subsets of pieces, from MO's cached sequences (debug/data/s1_coeffs_*.json; nothing recomputed).

    python debug/op52_pieces.py A1 1073741441   -> debug/data/op52_pieces_A1_<p>.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import s1_controls as C          # noqa: E402
from aabc_exact import params    # noqa: E402

DATA = HERE / 'data'


def xstar_exps(c, p, xs):
    kx = ((-C.I.to_fp(xs, p)) % p, 1)
    pts = C.points(c, p)
    if kx not in pts:
        return None
    return sorted(pts[kx]['exps'], key=Fraction)


def main(case, p):
    seqs, _ = C.a_coeffs(case, p)            # cached
    D1, D2 = params(case)[:2]
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    out = dict(case=case, prime=p, xstar=str(xs), pieces={}, subsets={})
    t = time.time()
    carriers = []
    for k, sq in seqs.items():
        po = C.piece_operator(sq, p)
        e = xstar_exps(po['c'], p, xs) if po else None
        out['pieces'][repr(k)] = dict(order=po and po['order'], degree=po and po['degree'], xstar=e)
        if e is not None:
            carriers.append(k)
        print(repr(k), po and po['order'], e, flush=True)
    # generic sums: all x*-carriers; all pieces
    for name, keys in (('xstar_carriers', carriers), ('all', list(seqs))):
        sub = {k: seqs[k] for k in keys}
        g = C.generic_operator(C.generic_seq(sub, p, C.SEEDS[0]), p)
        out['subsets'][name] = dict(n=len(keys), order=g['order'], degree=g['degree'],
                                    xstar=xstar_exps(g['c'], p, xs))
        print(name, g['order'], out['subsets'][name]['xstar'], flush=True)
    out['time'] = round(time.time() - t, 1)
    (DATA / f'op52_pieces_{case}_{p}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
