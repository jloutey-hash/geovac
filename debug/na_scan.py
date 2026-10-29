"""Cluster scan for the bond-angle operator (record NA; ansatz and spectrum as debug/ao_guess2.py).
For R = 1..Rmax, Q = Qstart, Qstart+step, ... while (R+1)(Q+1) <= rows/4: the singular spectrum at (R, Q).
DETECT = the biggest gap between consecutive log10 singular values is >= 15 decades, with k values below it;
then Q0 = Q - k + 1 is the minimal degree at that order.  The first detection (smallest (R+1)(Q0+1)) is FOUND; its
null vector at (R, Q0) by fixing the largest component to 1 and solving the Gram system.
    python debug/na_scan.py <coeffs.json> <Rmax> <Qstart> <step> <prec_bits> <out.json>
"""
import json
import sys
import time

import flint
from flint import arb, arb_mat, acb_mat

import ao_guess2 as G2


def scaled_matrix(F, R, Q):
    cols = G2.columns(F, R, Q)
    keys = sorted(cols)
    nrow = min(len(c) for c in cols.values())
    rows = []
    for l in range(nrow):
        row = [cols[k][l].mid() for k in keys]
        m = max(abs(v) for v in row)
        rows.append([v / m for v in row])
    cs = [max(abs(rows[i][j]) for i in range(nrow)) for j in range(len(keys))]
    A = arb_mat(nrow, len(keys), [rows[i][j] / cs[j] for i in range(nrow) for j in range(len(keys))])
    return keys, A, cs


def log10(a):
    v = float(a.mid())
    return float(a.mid().log()) / 2.302585092994046 if v > 0 else -9999.0


def spectrum(A):
    Gm = A.transpose() * A
    ev = acb_mat(Gm).eig(algorithm='approx')
    return sorted(log10(abs(e.real.mid()).sqrt()) for e in ev), Gm


def null_vector(A, cs):
    Gm = A.transpose() * A
    n = Gm.ncols()
    best = None
    for piv in range(n):
        idx = [j for j in range(n) if j != piv]
        M = arb_mat(n - 1, n - 1, [Gm[i, j].mid() for i in idx for j in idx])
        b = arb_mat(n - 1, 1, [-Gm[i, piv].mid() for i in idx])
        try:
            y = M.solve(b, nonstop=True)
        except Exception:
            continue
        v = [arb(0)] * n
        v[piv] = arb(1)
        for t, i in enumerate(idx):
            v[i] = y[t, 0].mid()
        r = A * arb_mat(n, 1, v)
        res = max(abs(float(r[i, 0].mid())) for i in range(r.nrows()))
        if best is None or res < best[0]:
            best = (res, v)
        if res < 1e-100:
            break
    res, v = best
    return res, [vi / c for vi, c in zip(v, cs)]


if __name__ == '__main__':
    path, Rmax, Qstart, step, prec, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), \
        int(sys.argv[5]), sys.argv[6]
    flint.ctx.prec = prec
    F = [arb(s.strip('[').split(' +/-')[0]) for s in json.load(open(path))['F']]
    L = len(F) - 1
    t0 = time.time()
    log, det = [], []
    for R in range(1, Rmax + 1):
        Q = Qstart
        while (R + 1) * (Q + 1) <= (L + 1 - R - Q) // 4:
            keys, A, cs = scaled_matrix(F, R, Q)
            lg, _ = spectrum(A)
            gaps = [lg[i + 1] - lg[i] for i in range(len(lg) - 1)]
            k = max(range(len(gaps)), key=lambda i: gaps[i])
            rec = dict(R=R, Q=Q, rows=A.nrows(), cols=A.ncols(), smallest=[round(v, 1) for v in lg[:5]],
                       gap=round(gaps[k], 1), k=k + 1, sec=round(time.time() - t0))
            log.append(rec)
            print(rec, flush=True)
            if gaps[k] >= 15:
                det.append(dict(R=R, Q0=Q - (k + 1) + 1, at_Q=Q, k=k + 1, gap=gaps[k], floor=lg[0]))
                break
            Q += step
    found = None
    if det:
        d = min(det, key=lambda d: ((d['R'] + 1) * (d['Q0'] + 1), d['R']))
        keys, A, cs = scaled_matrix(F, d['R'], d['Q0'])
        lg, _ = spectrum(A)
        res, v = null_vector(A, cs)
        found = dict(d, check_smallest=[round(x, 1) for x in lg[:3]], null_residual=res,
                     keys=[list(k) for k in keys], v=[x.str(int(prec * 0.3) - 10) for x in v])
        print('FOUND', {k: found[k] for k in ('R', 'Q0', 'k', 'gap', 'check_smallest', 'null_residual')}, flush=True)
    json.dump(dict(file=path, log=log, detections=det, found=found, seconds=time.time() - t0), open(out, 'w'), indent=0)
