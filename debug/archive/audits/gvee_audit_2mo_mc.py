r"""g_Vee defect audit, step B: EXACT i.i.d. Monte-Carlo of every <f_p f_q / r_r> on the 2-MO
sigma reference |Phi0|^2 = D_p(1,2) D_p(3,4) / 4 (Loewdin 1s_A(2.70), 1s_B(1.00); geovac.lih_r12ci).

Independent of the analytic path: no grid, no kernel, no Neumann.  Each 2-electron block is sampled
EXACTLY (no Metropolis): r1 ~ (m0^2+m1^2)/2, then r2 ~ chi^2 with chi = (m0(r1) m1 - m1(r1) m0)/norm
(the exact conditional of a 2-orbital determinant); every orbital^2 draw is a rejection sample from
the 2(u^2 oA^2 + v^2 oB^2) envelope (acceptance >= 1/2).  Chunks are i.i.d. -> plain error bars.

    python debug/gvee_audit_2mo_mc.py --geminal linexp --gam 0.25 --nchunk 400 --chunk 250000 --seed 1
Writes debug/data/gvee_audit_2mo_mc_<gem>_g<gam>_s<seed>.npz  (per-chunk means of f, v, f v, f f v).
"""
from __future__ import annotations

import os
import sys
import argparse
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")

from geovac.lih_r12ci.energy import ZA, ZB, N_A, N_B, S_AB
from geovac.lih_r12ci.fourbody import CENTER_A, CENTER_B

_sv, _U = np.linalg.eigh(np.array([[1.0, S_AB], [S_AB, 1.0]]))
X = _U @ np.diag(1.0 / np.sqrt(_sv)) @ _U.T
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def oA(r):
    return N_A * np.exp(-ZA * np.linalg.norm(r - CENTER_A, axis=-1))


def oB(r):
    return N_B * np.exp(-ZB * np.linalg.norm(r - CENTER_B, axis=-1))


def sample_1s(rng, n, Z, C):
    rr = rng.gamma(3.0, 1.0 / (2.0 * Z), size=n)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    return C[None, :] + rr[:, None] * v


def sample_chi2(rng, u, v):
    """one sample per row from chi^2, chi = u oA + v oB (u, v arrays (n,)), normalised chi."""
    n = len(u)
    out = np.empty((n, 3))
    todo = np.arange(n)
    while len(todo):
        uu = u[todo]; vv = v[todo]
        pA = uu ** 2 / (uu ** 2 + vv ** 2)
        isA = rng.uniform(size=len(todo)) < pA
        r = np.where(isA[:, None], sample_1s(rng, len(todo), ZA, CENTER_A),
                     sample_1s(rng, len(todo), ZB, CENTER_B))
        a = oA(r); b = oB(r)
        acc = rng.uniform(size=len(todo)) * 2.0 * (uu ** 2 * a ** 2 + vv ** 2 * b ** 2) < (uu * a + vv * b) ** 2
        out[todo[acc]] = r[acc]
        todo = todo[~acc]
    return out


def sample_block(rng, n):
    p = rng.integers(0, 2, size=n)
    u1 = X[0, p]; v1 = X[1, p]
    r1 = sample_chi2(rng, u1, v1)
    a1 = oA(r1); b1 = oB(r1)
    m0 = X[0, 0] * a1 + X[1, 0] * b1
    m1 = X[0, 1] * a1 + X[1, 1] * b1
    # chi = m0(r1) m1 - m1(r1) m0  in AO coefficients
    u2 = m0 * X[0, 1] - m1 * X[0, 0]
    v2 = m0 * X[1, 1] - m1 * X[1, 0]
    r2 = sample_chi2(rng, u2, v2)
    return r1, r2


def fgem(kind, g, d):
    return d * np.exp(-g * d) if kind == 'linexp' else np.exp(-g * d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp')
    ap.add_argument('--gam', type=float, nargs='+', default=[0.25])
    ap.add_argument('--nchunk', type=int, default=100)
    ap.add_argument('--chunk', type=int, default=200000)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(1000 + a.seed)
    ng = len(a.gam)
    Fm = np.zeros((a.nchunk, ng, 6)); Vm = np.zeros((a.nchunk, 6))
    FVm = np.zeros((a.nchunk, ng, 6, 6)); FFVm = np.zeros((a.nchunk, ng, 6, 6, 6))
    t0 = time.time()
    for c in range(a.nchunk):
        r1, r2 = sample_block(rng, a.chunk)
        r3, r4 = sample_block(rng, a.chunk)
        R = {1: r1, 2: r2, 3: r3, 4: r4}
        d = np.stack([np.linalg.norm(R[i] - R[j], axis=1) for i, j in PAIRS6], axis=1)   # (n,6)
        v = 1.0 / d
        Vm[c] = v.mean(0)
        for ig, g in enumerate(a.gam):
            f = fgem(a.geminal, g, d)
            Fm[c, ig] = f.mean(0)
            FVm[c, ig] = np.einsum('ni,nk->ik', f, v) / a.chunk
            FFVm[c, ig] = np.einsum('ni,nj,nk->ijk', f, f, v, optimize=True) / a.chunk
        if c % 20 == 0:
            print(f"  chunk {c}/{a.nchunk}  {time.time()-t0:.0f}s", flush=True)
    tag = f"{a.geminal}_g{'-'.join(f'{g:g}' for g in a.gam)}_s{a.seed}"
    np.savez(os.path.join(DATA, f"gvee_audit_2mo_mc_{tag}.npz"), gam=np.array(a.gam), Fm=Fm, Vm=Vm,
             FVm=FVm, FFVm=FFVm, chunk=a.chunk)
    print(f"done {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
