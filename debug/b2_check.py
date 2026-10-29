"""Record B2, CHECKS ONLY (data enters here and nowhere upstream).
    python b2_check.py point <p> L0 [L0 ...] : the recurrence evaluated pointwise (kernel of the cyclic rows at the
                                              integer L0, no reconstruction) applied to the data at L0..L0+m
    python b2_check.py rec <file>            : a reconstructed recurrence on the generic sequence and on every
                                              component, all L with L + m <= 300
"""
import json
import sys

import flint

PRIMES = [(1 << 61) - 1, 4611686018427387847]


def seqs(p):
    d = json.load(open(f'data/b2_check_seq_{p}.json'))
    return d['gen'], d['components']


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'point':
        import b2_derive as Dv
        from bu_derive import stack
        p = int(sys.argv[2])
        S = Dv.B2Sys(p)
        m = Dv.cyclic_order(S, 123456789)
        gen, comps = seqs(p)
        for L0 in map(int, sys.argv[3:]):
            try:
                M = stack(S.rows(L0, m), p)
            except KeyError:
                print('p', p, 'L0', L0, 'IBP structure differs at this integer point (skipped)', flush=True)
                continue
            X, nul = M.transpose().nullspace()
            res = []
            for b in range(nul):
                v = [X[i, b] for i in range(m + 1)]
                s = sum((v[k] * flint.nmod(gen[L0 + k], p) for k in range(m + 1)), flint.nmod(0, p))
                bad = [name for name, sq in comps.items()
                       if int(sum((v[k] * flint.nmod(sq[L0 + k], p) for k in range(m + 1)), flint.nmod(0, p))) != 0]
                res.append((int(s) == 0, len(bad)))
            print('p', p, 'L0', L0, 'order', m, 'kernel dim', nul, 'gen zero / #components nonzero', res, flush=True)
    elif mode == 'rec':
        d = json.load(open(sys.argv[2]))
        p = d['p']
        polys = [flint.nmod_poly(c, p) for c in d['polys']]
        m = len(polys) - 1
        gen, comps = seqs(p)
        out = {}
        for name, s in [('gen', gen)] + list(comps.items()):
            bad = [l for l in range(len(s) - m)
                   if int(sum((polys[k](flint.nmod(l, p)) * flint.nmod(s[l + k], p) for k in range(m + 1)),
                              flint.nmod(0, p))) != 0]
            out[name] = (len(bad), len(s) - m, bad[:5])
        nz = {k: v for k, v in out.items() if v[0]}
        print('p', p, 'gen', out['gen'], 'components with failures', len(nz), nz if nz else '', flush=True)
