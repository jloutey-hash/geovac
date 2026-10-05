r"""r12 t^2 STEP 1, TASK B (part 2) -- the two-f-edge Coulomb TRIANGLE for f = u t^2 e^{-gam u}.

T1 = <f13 f23 / r12>, T2 = <f12 f23 / r13> over |Phi0|^2 (2-MO reference), the one non-separable
term class of g_Vee.  With f(i,k) = k(u_ik) sum_c alpha_c(r_i) beta_c(r_k) (rank 3, k = u e^{-gam u}),
    f13 f23 = k13 k23 sum_{c,d} alpha_c(1) alpha_d(2) [beta_c beta_d](3)
so each (c,d) is the EXISTING triangle reduction (geovac.lih_r12ci.triangle.triangle_raw) with the
shared-vertex density mu -> mu beta_c beta_d and the Coulomb-edge densities reweighted by alpha_c,
alpha_d; the azimuthal-Fourier kernel F_m is built for k(u) instead of e^{-gam u}.
(T2: shared vertex 2 carries a D_p component; electron 1 its partner; electron 3 the marginal.)

Checked against brute-force MC with exact i.i.d. |Phi0|^2 sampling (r12t2_step1_separable.sample_pairs).
Positive control: f = e^{-0.5u} through the same wrapper (ONE weight) vs the same MC.
Run:  python debug/r12t2_step1_triangle.py [--nmc 4e7]
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'debug'))
from geovac.lih_r12ci import kernels as KG                                   # noqa: E402
from geovac.lih_r12ci import triangle as TR                                  # noqa: E402
from geovac.lih_r12ci.hVee import P00, P11, P01, rho_g                       # noqa: E402
from geovac.lih_r12ci.fourbody import CENTER_A                               # noqa: E402
from r12t2_step1_separable import sample_pairs                               # noqa: E402

rA_f = KG.rA.reshape(-1)
ONES = np.ones_like(rA_f)
COMPS = [(P00, P11, 1.0, 'P00', 'P11'), (P11, P00, 1.0, 'P11', 'P00'), (P01, P01, -2.0, 'P01', 'P01')]


def terms(kind):
    if kind == 'exp':
        return [(1.0, ONES, ONES)]
    return [(1.0, rA_f ** 2, ONES), (1.0, ONES, rA_f ** 2), (-2.0, rA_f, rA_f)]


def analytic(gam, kind, rank=30):
    func = (lambda d: np.exp(-gam * d)) if kind == 'exp' else (lambda d: d * np.exp(-gam * d))
    Fm = [KG.build_kernel_m_gen(func, m, nphi=48) for m in range(TR.MMAX + 1)]
    TR._MODE_CACHE.clear()
    T = terms(kind)
    T1 = T2 = 0.0
    for ic, (cc, ac, bc) in enumerate(T):           # edge (1,3): alpha_c on 1, beta_c on 3 (T1)
        for idd, (cd, ad, bd) in enumerate(T):      # edge (2,3): alpha_d on 2, beta_d on 3
            w = cc * cd
            # T1 = <f13 f23 /r12>: mu = rho * beta_c beta_d (vertex 3), left = g alpha_c, right = h alpha_d
            for g, h, k, gk, hk in COMPS:
                v, _ = TR.triangle_raw(Fm, rho_g * bc * bd, g * ac, h * ad,
                                       f"T1|{kind}|{gam}|{ic}{idd}", rank=rank)
                T1 += w * k * v
            # T2 = <f12 f23 /r13>: f12 = k12 alpha_c(1)beta_c(2); f23 = k23 alpha_d(2) beta_d(3)
            #   shared vertex 2 carries h * beta_c * alpha_d; e1 carries g alpha_c; e3 carries rho beta_d
            for g, h, k, gk, hk in COMPS:
                v, _ = TR.triangle_raw(Fm, h * bc * ad, g * ac, rho_g * bd,
                                       f"T2|{kind}|{gam}|{ic}{idd}|{hk}", rank=rank)
                T2 += w * k * v
    return 0.25 * T1, 0.25 * T2


def fval(kind, gam, d, rAi, rAj):
    return np.exp(-gam * d) if kind == 'exp' else d * (rAi - rAj) ** 2 * np.exp(-gam * d)


def mc(cases, n_total, batch=1_000_000, seed=777):
    rng = np.random.default_rng(seed)
    acc = {c: [] for c in cases}
    t0 = time.time()
    for b in range(int(n_total // batch)):
        u1, u2 = sample_pairs(batch, rng); d1, _ = sample_pairs(batch, rng)
        r = [u1, u2, d1]
        rA = [np.linalg.norm(x - CENTER_A, axis=-1) for x in r]
        D = lambda i, j: np.linalg.norm(r[i] - r[j], axis=-1)
        d12, d13, d23 = D(0, 1), D(0, 2), D(1, 2)
        for (gam, kind) in cases:
            f13 = fval(kind, gam, d13, rA[0], rA[2]); f23 = fval(kind, gam, d23, rA[1], rA[2])
            f12 = fval(kind, gam, d12, rA[0], rA[1])
            acc[(gam, kind)].append([(f13 * f23 / d12).mean(), (f12 * f23 / d13).mean()])
        if b % 10 == 0:
            print(f"   MC batch {b+1}  ({time.time()-t0:.0f}s)", flush=True)
    out = {}
    for c, rows in acc.items():
        R_ = np.array(rows)
        out[c] = (R_.mean(0), R_.std(0, ddof=1) / np.sqrt(len(R_)))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument('--nmc', type=float, default=4e7)
    ap.add_argument('--exact-neumann', action='store_true')
    args = ap.parse_args()
    KG.USE_EXACT_NEUMANN = bool(args.exact_neumann)
    cases = [(0.5, 'exp'), (1.5, 'ut2')]
    res = {}
    for gam, kind in cases:
        t0 = time.time()
        res[(gam, kind)] = analytic(gam, kind)
        print(f"[reduced] {kind} gam={gam}: T1={res[(gam,kind)][0]:.7f} T2={res[(gam,kind)][1]:.7f}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
    m = mc(cases, args.nmc)
    rows = []
    for c in cases:
        for i, nm in enumerate(('T1=<f13f23/r12>', 'T2=<f12f23/r13>')):
            av = res[c][i]; mv, me = m[c][0][i], m[c][1][i]
            rows.append(dict(case=f"{c[1]} g={c[0]}", term=nm, reduced=av, mc=mv, mc_err=me,
                             rel=abs(av - mv) / abs(mv), nsig=abs(av - mv) / me))
            print(f"{c[1]+' g='+str(c[0]):12s} {nm:17s} reduced {av:.7f}  MC {mv:.7f} +/- {me:.1e}"
                  f"  rel {abs(av-mv)/abs(mv):.1e}  {abs(av-mv)/me:.1f} sigma")
    with open(os.path.join(ROOT, 'debug', 'data', 'r12t2_step1_triangle' + ('_exactN' if KG.USE_EXACT_NEUMANN else '') + '.json'), 'w') as fh:
        json.dump(dict(nmc=args.nmc, exact_neumann=KG.USE_EXACT_NEUMANN, rows=rows), fh, indent=1, default=float)
