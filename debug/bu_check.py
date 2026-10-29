"""Record BU, CHECKS ONLY (data enters here and nowhere upstream).
    python bu_check.py rec          : the derived order-11 recurrence (data/bu_recurrence_<p>.json) on the generic
                                      sequence and on each component, L = 0..1000, both primes
    python bu_check.py angle R Q    : a derived angle operator (data/bu_angle_<p>_<R>_<Q>.json) applied to the data
"""
import json
import sys

import flint

import na_modguess as G

PRIMES = [(1 << 61) - 1, 4611686018427387847]
VER = 'v2' if '--v2' in sys.argv else 'v1'
if '--v2' in sys.argv:
    sys.argv.remove('--v2')


def seqs(p):
    d = json.load(open(f'data/bu_check_seq_{p}.json'))
    return d['gen'], d['components']


def check_rec(p):
    rec = json.load(open(f'data/bu_recurrence_{VER}_{p}.json'))
    polys = [flint.nmod_poly(c, p) for c in rec['polys']]
    m = len(polys) - 1
    gen, comps = seqs(p)
    out = {}
    for name, s in [('gen', gen)] + list(comps.items()):
        bad = 0
        for l in range(len(s) - m):
            lv = flint.nmod(l, p)
            v = sum((polys[k](lv) * flint.nmod(s[l + k], p) for k in range(m + 1)), flint.nmod(0, p))
            bad += int(v) != 0
        out[name] = (bad, len(s) - m)
    return out


def apply_op(keys, vec, s, p):
    B = G.Basis([flint.nmod(v, p) for v in s], p)
    cols = {k: B.col(*k) for k in keys}
    n = min(len(c) for c in cols.values())
    res = [sum((flint.nmod(a, p) * cols[k][l] for k, a in zip(keys, vec)), flint.nmod(0, p)) for l in range(n)]
    return res


if __name__ == '__main__':
    if sys.argv[1] == 'rec':
        for p in PRIMES:
            print(p, check_rec(p), flush=True)
    elif sys.argv[1] == 'angle':
        R, Q = int(sys.argv[2]), int(sys.argv[3])
        for p in PRIMES:
            d = json.load(open(f'data/bu_angle_{VER}_{p}_{R}_{Q}.json'))
            keys = [tuple(k) for k in d['keys']]
            gen, comps = seqs(p)
            for vi, vec in enumerate(d['vecs']):
                r = apply_op(keys, vec, gen, p)
                nz = sum(int(v) != 0 for v in r)
                print(p, 'vec', vi, 'rows', len(r), 'nonzero on gen', nz,
                      {n: sum(int(v) != 0 for v in apply_op(keys, vec, s, p)) for n, s in comps.items()}, flush=True)
