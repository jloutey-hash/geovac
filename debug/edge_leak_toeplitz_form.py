"""Order-one defect on the true CH window of S^3 in Toeplitz/Hankel (Widom) form.

Notation: W = true window |lambda| <= n_max + 1/2, P its projection, Q = 1 - P, D diagonal.
  T(x)  = P x P   (compression, "Toeplitz")
  H(x)  = Q x P   (window -> outside, "Hankel")
  Ht(x) = P x Q   (outside -> window)
Widom-type identity (elementary, P + Q = 1):  T(xy) = T(x) T(y) + Ht(x) H(y).
Order-one defect, using [D,a]b = b[D,a] for functions (continuum order one):
  C(a,b) = [[D_W, a_W], b_W] = Ht(b) H([D,a]) - Ht([D,a]) H(b)
         = T(b[D,a]) - T(b)T([D,a]) - (T([D,a]b) - T([D,a])T(b)).
Checks (1) block form vs the guarded edge-crossing sum (tests/test_paper32_window_order_one.py)
and vs the direct commutator, (2) the Widom identity itself, (3) operator / trace norm scaling.

Run with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(ROOT, 'tests')


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(TESTS, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gb = _load('gb_tw', 'test_paper32_window_gauge_breaking.py')
o1 = _load('o1_tw', 'test_paper32_window_order_one.py')

PAIRS = o1.PAIRS  # the four D^{1/2}_{M M'}


def build(n_max, J2max):
    S = gb.Sphere(J2max)
    funcs = [S.multiplier(M, Mp).toarray() for M, Mp in PAIRS]
    inW = np.abs(S.lam) <= n_max + 0.5
    return S, funcs, inW


def blocks(x, w, out):
    return x[np.ix_(w, w)], x[np.ix_(out, w)], x[np.ix_(w, out)]   # T, H, Ht


def defect_block(lam, a, b, inW):
    w, out = np.nonzero(inW)[0], np.nonzero(~inW)[0]
    Da = lam[:, None] * a - a * lam[None, :]
    _, H_Da, Ht_Da = blocks(Da, w, out)
    _, H_b, Ht_b = blocks(b, w, out)
    return Ht_b @ H_Da - Ht_Da @ H_b


def check_agreement(n_list=(1, 2, 3, 4)):
    rows = []
    for n in n_list:
        S, funcs, inW = build(n, 2 * n + 2)            # the guarded construction
        w, out = np.nonzero(inW)[0], np.nonzero(~inW)[0]
        worst_be, worst_bd, worst_widom, worst_small, scale = 0., 0., 0., 0., 0.
        S2, funcs2, inW2 = build(n, n + 3)            # smaller ambient, used for the norm sweep
        for ia, a in enumerate(funcs):
            for ib, b in enumerate(funcs):
                Cd = o1.defect_direct(S.lam, a, b, inW)
                Ce = o1.defect_edge(S.lam, a, b, inW)
                Cb = defect_block(S.lam, a, b, inW)
                s = max(np.abs(Cd).max(), 1e-300)
                scale = max(scale, s)
                worst_be = max(worst_be, np.abs(Cb - Ce).max() / s)
                worst_bd = max(worst_bd, np.abs(Cb - Cd).max() / s)
                # Widom identity T(ab) = T(a)T(b) + Ht(a)H(b)
                Ta, Ha, Hta = blocks(a, w, out)
                Tb, Hb, Htb = blocks(b, w, out)
                Tab = (a @ b)[np.ix_(w, w)]
                r = Tab - Ta @ Tb - Hta @ Hb
                worst_widom = max(worst_widom, np.abs(r).max() / max(np.abs(Hta @ Hb).max(), 1e-300))
                Cs = defect_block(S2.lam, funcs2[ia], funcs2[ib], inW2)
                worst_small = max(worst_small, np.abs(Cs - Cd).max() / s)
        rows.append(dict(n_max=n, dimW=int(inW.sum()), max_entry=scale,
                         rel_block_vs_edge=worst_be, rel_block_vs_direct=worst_bd,
                         rel_widom_identity=worst_widom, rel_small_ambient=worst_small))
        print(rows[-1], flush=True)
    return rows


def norms(C):
    s = np.linalg.svd(C, compute_uv=False)
    return float(s[0]), float(s.sum()), float(np.sqrt((s ** 2).sum()))


def hermitian_coords(funcs):
    """Four real coordinates of SU(2) from D^{1/2}_{-+} and D^{1/2}_{++}: Re, Im parts."""
    a, b = funcs[0], funcs[1]
    return [(a + a.conj().T) / 2, (a - a.conj().T) / 2j, (b + b.conj().T) / 2, (b - b.conj().T) / 2j]


def norm_sweep(n_list):
    out = []
    for n in n_list:
        t0 = time.time()
        S, funcs, inW = build(n, n + 3)
        diag = [norms(defect_block(S.lam, f, f, inW)) for f in funcs]
        allp = [norms(defect_block(S.lam, a, b, inW)) for a in funcs for b in funcs]
        X = hermitian_coords(funcs)
        sumX = sum(defect_block(S.lam, x, x, inW) for x in X)
        rec = dict(n_max=n, dimW=int(inW.sum()),
                   diag_op=[d[0] for d in diag], diag_tr=[d[1] for d in diag], diag_hs=[d[2] for d in diag],
                   max16_op=max(p[0] for p in allp), max16_tr=max(p[1] for p in allp),
                   coordsum=norms(sumX),
                   coordsum_is_scalar_offdiag=float(np.abs(sumX - np.diag(np.diag(sumX))).max()),
                   secs=time.time() - t0)
        out.append(rec)
        print(n, rec['dimW'], 'op', np.round(rec['diag_op'], 6), 'tr', np.round(rec['diag_tr'], 4),
              'max16', round(rec['max16_op'], 6), round(rec['max16_tr'], 4),
              'coordsum', np.round(rec['coordsum'], 5), f"{rec['secs']:.1f}s", flush=True)
    return out


def slopes(ns, ys):
    ln, ly = np.log(ns), np.log(ys)
    local = [(ly[i + 1] - ly[i]) / (ln[i + 1] - ln[i]) for i in range(len(ns) - 1)]
    p = np.polyfit(ln, ly, 1)
    return local, float(p[0])


if __name__ == '__main__':
    nmax_sweep = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    agree = check_agreement()
    sweep = norm_sweep(range(1, nmax_sweep + 1))
    ns = np.array([r['n_max'] for r in sweep], float)
    fits = {}
    for key, ys in (('diag_op_f0', [r['diag_op'][0] for r in sweep]),
                    ('diag_tr_f0', [r['diag_tr'][0] for r in sweep]),
                    ('diag_hs_f0', [r['diag_hs'][0] for r in sweep]),
                    ('max16_op', [r['max16_op'] for r in sweep]),
                    ('max16_tr', [r['max16_tr'] for r in sweep]),
                    ('coordsum_op', [r['coordsum'][0] for r in sweep]),
                    ('coordsum_tr', [r['coordsum'][1] for r in sweep]),
                    ('dimW', [r['dimW'] for r in sweep])):
        ys = np.array(ys)
        loc, _ = slopes(ns, ys)
        m = ns >= 4
        _, tail = slopes(ns[m], ys[m])
        fits[key] = dict(local=loc, tail_slope_n_ge_4=tail, values=ys.tolist())
        print(key, 'local slopes', np.round(loc, 3), 'tail fit', round(tail, 3))
    os.makedirs(os.path.join(ROOT, 'debug', 'data'), exist_ok=True)
    with open(os.path.join(ROOT, 'debug', 'data', 'edge_leak_toeplitz_form.json'), 'w') as fh:
        json.dump(dict(agreement=agree, sweep=sweep, fits=fits), fh, indent=1)
