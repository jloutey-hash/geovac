"""KS route Q (record docs/predictions/2026-10-28_kite_speed_test.md): float64 graded tanh-sinh quadrature of the
route-B double integral (pf_common's kernel), numpy-vectorized, at each node of the shared outer grid.
assign 'correct': k2 = (1-u-v)^2 z, v split at the ridge v = 1-u; 'wrong': k2 = (u-v)^2 z (negative control C5),
same rule.  usage (child mode, one timed run): python ks_route_q.py <T12|T9> <h> <correct|wrong>"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json, sys, time
import numpy as np
from ks_check import kernel, A, B, G1, G2
from ks_check2 import ts01
import ks_common as C


def phi_q(wn, h, kmax, assign="correct"):
    z = wn * wn
    u, um, wu = ts01(h, kmax)
    U = u * um
    M = np.sqrt(A * A * u + B * B * um + U * z)
    cu = wu * np.sqrt(U) / M ** 3
    r, rm, wr = ts01(h, kmax)
    tot = 0.0
    for side in (0, 1):
        if side == 0:
            L = um[:, None]; v = L * r; vm = 1 - v; d = L * rm
        else:
            L = u[:, None]; vm = L * rm; v = 1 - vm; d = -L * r
        if assign == "wrong":
            d = u[:, None] - v
        wv = L * wr
        V = v * vm
        N = np.sqrt(G1 * G1 * vm + G2 * G2 * v + V * z)
        cv = wv * V / N ** 5
        K = kernel(M[:, None], N, d * d * z)
        tot += np.sum(cu[:, None] * cv * K)
    return 2 * G1 * G2 * A ** 3 * B ** 3 * tot


def run(target, h, assign="correct"):
    nd = C.outer_nodes(target)
    w = np.array([float(x) for x, _ in nd]); dw = np.array([float(x) for _, x in nd])
    kmax = int(round(5 / h))
    t0 = time.process_time()
    pre = t0
    phi = np.array([phi_q(x, h, kmax, assign) for x in w])
    t1 = time.process_time()
    Y = C.stage2(w, dw, phi)
    t2 = time.process_time()
    return Y, {"stage1": t1 - t0, "stage2": t2 - t1, "before_run": pre}


if __name__ == "__main__":
    tg, h, assign = sys.argv[1], eval(sys.argv[2]), sys.argv[3]
    Y, tm = run(tg, h, assign)
    print(json.dumps({"route": "Q", "target": tg, "h": h, "assign": assign, "Y": list(map(float, Y)),
                      "rel": C.rel_errs(Y), "cpu": tm}))
