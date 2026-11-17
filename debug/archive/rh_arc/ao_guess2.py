"""Large-budget operator search for the bond-angle ODE (record AO2).  Same ansatz and Legendre-local operators as
debug/ao_guess.py:  sum_{i<=R, j<=Q} a_ij x^j T^i f = 0,  T = (1 - x^2) d/dx.

At one (R, Q) it returns ALL singular values of the row- and column-scaled matrix (eigenvalues of the Gram matrix in
ball arithmetic at high precision).  If an operator of order R and degree Q0 <= Q exists, then x^m L (m = 0..Q-Q0)
are null vectors too, so a cluster of Q - Q0 + 1 singular values sits at the data-noise floor below a cliff;
approximate (fake) operators give a graded, cliff-free spectrum.

    python debug/ao_guess2.py <coeffs.json> <R> <Q> <prec_bits> <out.json>
"""
import json
import sys
import time

import flint
from flint import arb, arb_mat, acb_mat


def X(s):
    out = []
    for l in range(len(s) - 1):
        a = arb(l) / (2 * l - 1) * s[l - 1] if l >= 1 else arb(0)
        out.append(a + arb(l + 1) / (2 * l + 3) * s[l + 1])
    return out


def T(s):
    out = []
    for l in range(len(s) - 1):
        a = arb((l - 1) * l) / (2 * l - 1) * s[l - 1] if l >= 1 else arb(0)
        out.append(arb((l + 1) * (l + 2)) / (2 * l + 3) * s[l + 1] - a)
    return out


def columns(F, R, Q):
    cols = {}
    ti = list(F)
    for i in range(R + 1):
        xj = ti
        for j in range(Q + 1):
            cols[(i, j)] = xj
            xj = X(xj)
        ti = T(ti)
    return cols


def spectrum(F, R, Q):
    cols = columns(F, R, Q)
    keys = sorted(cols)
    nrow = min(len(c) for c in cols.values())
    rows = []
    for l in range(nrow):
        row = [cols[k][l].mid() for k in keys]
        m = max(abs(v) for v in row)
        rows.append([v / m for v in row])
    ncol = len(keys)
    cs = []
    for j in range(ncol):
        cs.append(max(abs(rows[i][j]) for i in range(nrow)))
    A = arb_mat(nrow, ncol, [rows[i][j] / cs[j] for i in range(nrow) for j in range(ncol)])
    G = A.transpose() * A
    ev = acb_mat(G).eig(algorithm='approx')
    sv = sorted(abs(e.real.mid()).sqrt() for e in ev)
    return keys, nrow, sv


if __name__ == '__main__':
    path, R, Q, prec, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    flint.ctx.prec = prec
    F = [arb(s.strip('[').split(' +/-')[0]) if s.startswith('[') else arb(s) for s in json.load(open(path))['F']]
    t0 = time.time()
    keys, nrow, sv = spectrum(F, R, Q)
    lg = [float(s.mid().log()) / 2.302585092994046 if float(s.mid()) > 0 else -9999 for s in sv]
    gaps = [lg[i + 1] - lg[i] for i in range(len(lg) - 1)]
    k = max(range(len(gaps)), key=lambda i: gaps[i])
    json.dump(dict(file=path, R=R, Q=Q, nrow=nrow, ncol=len(keys), seconds=time.time() - t0, log10_sv=lg,
                   biggest_gap_after=k + 1, biggest_gap_decades=gaps[k]), open(out, 'w'), indent=0)
    print('R', R, 'Q', Q, 'rows', nrow, 'cols', len(keys), 'sec', round(time.time() - t0),
          'smallest log10', [round(v, 1) for v in lg[:6]], 'biggest gap', round(gaps[k], 1), 'after', k + 1, flush=True)
