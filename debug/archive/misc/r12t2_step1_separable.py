r"""r12 t^2 STEP 1, TASK B -- the u t^2 g(u) geminal through the EXISTING RI-free reductions.

Separability (exact, finite rank 3):
    f(1,2) = u t^2 g(u),  t = r_1A - r_2A,  t^2 = r_1A^2 + r_2A^2 - 2 r_1A r_2A
           = k(u) [ r_1A^2 * 1 + 1 * r_2A^2 - 2 r_1A * r_2A ],    k(u) = u g(u).
So every f-dressing field is a reweighted dressing of the EXISTING kernel k:
    Psi^f_P(r) = r_A^2 Psi^k_P(r) + Psi^k_{r_A^2 P}(r) - 2 r_A Psi^k_{r_A P}(r),
and the product kernels needed by sigma^2 / h_Vee are again separable:
    f^2   = k2(u) t^4,  k2 = u^2 g^2,  t^4 = sum_m C(4,m)(-1)^m r_1A^(4-m) r_2A^m   (rank 5)
    f/u   = g(u) t^2                                                            (rank 3).
With g = e^{-gam u}: k = linexp kernel (existing), k2 = u^2 e^{-2 gam u}, f/u = t^2 x Slater kernel.

This script feeds those REWEIGHTED densities/fields into the package's own generic contraction
(geovac.lih_r12ci.hVee.cov_FA_FB, the sigma^2 / h_Vee engine) on the 2-MO reference
Phi0 = |m0^2 m1^2| (Loewdin 1s_A(2.7), 1s_B(1.0), R = 3.015), and checks against brute-force
real-space Monte-Carlo with EXACT i.i.d. sampling of |Phi0|^2 (rejection; no Metropolis bias).

Quantities: Fbar, sigma^2 = Var F (contains the 3-body share-vertex <f12 f13> and the
block-correlated 4-body <f13 f24>), h_Vne = Cov[F, V_ne], h_Vee = Cov[F, V_ee] (contains the
3-body <f12/r13> Coulomb share-vertex, and the 4-body disjoint <f12/r34>).
Positive control: the plain exp(-0.5 u) geminal through the same wrapper vs the banked
Fbar = 1.880743, sigma^2 = 0.13699.

Run:  python debug/r12t2_step1_separable.py [--nmc 2e8]
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
from math import comb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from geovac.lih_r12ci import kernels as KG                        # noqa: E402
from geovac.lih_r12ci.kernels import dens, grid_int, build_kernel_gen, geo_f, a   # noqa: E402
from geovac.lih_r12ci import hVee as HV                            # noqa: E402
from geovac.lih_r12ci.energy import ZA, ZB, N_A, N_B, Z_A, Z_B     # noqa: E402
from geovac.lih_r12ci.fourbody import CENTER_A, CENTER_B  # noqa: E402
from geovac.lih_r12ci.basis import X                               # noqa: E402

AOP = ['aa', 'ab', 'bb']
NPHI = 28
rA_f = KG.rA.reshape(-1); rB_f = KG.rB.reshape(-1)


# --------------------------------------------------------------------------- #
# separable geminal -> package kernel object {W (3x3 AO-pair), Psi (dressing fields)}
#   f(1,2) = K(u) * sum_c coef * wa(r1) * wb(r2)   (symmetric overall)
# --------------------------------------------------------------------------- #
def sep_kernel_object(Kmat, terms):
    Psi = {}
    for k in AOP:
        field = np.zeros_like(dens[k])
        for coef, wa, wb in terms:
            field += coef * wa * (a ** 3 * (Kmat @ (geo_f * wb * dens[k])))
        Psi[k] = field
    W = np.array([[grid_int(dens[u] * Psi[v]) for v in AOP] for u in AOP])
    return dict(W=W, Psi=Psi)


def t2n_terms(n):
    """t^(2n) = sum_m C(2n,m) (-1)^m r1A^(2n-m) r2A^m."""
    return [(comb(2 * n, m) * (-1) ** m, rA_f ** (2 * n - m), rA_f ** m) for m in range(2 * n + 1)]


ONE = [(1.0, np.ones_like(rA_f), np.ones_like(rA_f))]


def analytic(gam, kind):
    t0 = time.time()
    if kind == 'exp':           # control: f = e^{-gam u}
        Kf = build_kernel_gen(lambda d: np.exp(-gam * d), nphi=NPHI)
        Kf2 = build_kernel_gen(lambda d: np.exp(-2 * gam * d), nphi=NPHI)
        A = sep_kernel_object(Kf, ONE); C2 = sep_kernel_object(Kf2, ONE)
        CY = HV.make_kernel_Y(gam)      # package Yukawa (closed iso + MC abab): grid diag singular
    else:                       # f = u t^2 e^{-gam u}
        Kk = build_kernel_gen(lambda d: d * np.exp(-gam * d), nphi=NPHI)            # linexp kernel
        Kk2 = build_kernel_gen(lambda d: d * d * np.exp(-2 * gam * d), nphi=NPHI)
        Kg = build_kernel_gen(lambda d: np.exp(-gam * d), nphi=NPHI)
        A = sep_kernel_object(Kk, t2n_terms(1))
        C2 = sep_kernel_object(Kk2, t2n_terms(2))
        CY = sep_kernel_object(Kg, t2n_terms(1))        # f/u = t^2 e^{-gam u}
    sig2, info = HV.cov_FA_FB(A, A, C2)
    Fbar = info['FbarA']
    Kc = HV.make_kernel_coul()
    hVee, infoV = HV.cov_FA_FB(A, Kc, CY)
    E2 = infoV['FbarB']
    # h_Vne (kernels.__main__ formula, generic in the dressing fields)
    P00, P11, P01, rho = HV.P00, HV.P11, HV.P01, HV.rho_g
    vne = -Z_A / rA_f - Z_B / rB_f
    psi = lambda c: HV._psi(A['Psi'], c)
    Psi00, Psi11, Psi01, Psirho = psi(HV.c00), psi(HV.c11), psi(HV.c01), psi(HV.crho)
    a1, b1 = HV._a1b1(A['W'])
    Vne_exp = 2.0 * grid_int(rho * vne); vbar = 0.5 * grid_int(rho * vne)
    Pd = {'00': P00, '11': P11, '01': P01}; Psd = {'00': Psi00, '11': Psi11, '01': Psi01}
    comps = [('00', '11', 1.0), ('11', '00', 1.0), ('01', '01', -2.0)]
    Ta = 0.5 * sum(c * grid_int(vne * Pd[g] * Psd[h]) for g, h, c in comps)
    Tb = vbar * a1
    Tc = 0.25 * grid_int(vne * rho * Psirho)
    Td = 0.5 * sum(c * grid_int(vne * Pd[g]) * grid_int(Pd[h] * Psirho / 2.0) for g, h, c in comps)
    hVne = 4.0 * (Ta + Tb + 2 * Tc + 2 * Td) - Fbar * Vne_exp
    return dict(Fbar=Fbar, sigma2=sig2, hVne=hVne, hVee=hVee, E2=E2, Vne=Vne_exp,
                parts_sigma2={k: info[k] for k in ('same', 'dj_intra', 'share_delta', 'share_gsh', 'dj')},
                parts_hVee={k: infoV[k] for k in ('same', 'dj_intra', 'share_delta', 'share_gsh', 'dj')},
                seconds=time.time() - t0)


# --------------------------------------------------------------------------- #
# brute-force MC: exact i.i.d. sampling of |Phi0|^2 = D(1,2) D(3,4)
#   D(r1,r2) ∝ (A1 B2 - B1 A2)^2 (span{m0,m1} = span{1s_A,1s_B});
#   proposal q = 1/2 [A1^2 B2^2 + B1^2 A2^2], accept prob (A1B2-B1A2)^2 / (2(A1^2B2^2+B1^2A2^2)).
# --------------------------------------------------------------------------- #
def _A(r): return N_A * np.exp(-ZA * np.linalg.norm(r - CENTER_A, axis=-1))
def _B(r): return N_B * np.exp(-ZB * np.linalg.norm(r - CENTER_B, axis=-1))


def sample_1s(n, Z, center, rng):
    r = rng.gamma(3.0, 1.0 / (2 * Z), size=n)
    ct = rng.uniform(-1, 1, size=n); phi = rng.uniform(0, 2 * np.pi, size=n)
    st = np.sqrt(np.maximum(1 - ct * ct, 0.0))
    return np.stack([r * st * np.cos(phi), r * st * np.sin(phi), r * ct], axis=1) + center


def sample_pairs(n, rng):
    out1, out2, got = [], [], 0
    while got < n:
        m = int(1.15 * (n - got) / 0.45) + 1000
        swap = rng.random(m) < 0.5
        x = sample_1s(m, ZA, CENTER_A, rng); y = sample_1s(m, ZB, CENTER_B, rng)
        r1 = np.where(swap[:, None], y, x); r2 = np.where(swap[:, None], x, y)
        A1, B1, A2, B2 = _A(r1), _B(r1), _A(r2), _B(r2)
        num = (A1 * B2 - B1 * A2) ** 2; den = 2 * (A1 ** 2 * B2 ** 2 + B1 ** 2 * A2 ** 2)
        acc = rng.random(m) < num / den
        out1.append(r1[acc]); out2.append(r2[acc]); got += acc.sum()
    return np.concatenate(out1)[:n], np.concatenate(out2)[:n]


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def mc(gams, n_total, batch=1_000_000, seed=20261005):
    """returns per-(gam,kind) MC estimates with batch-means errors."""
    rng = np.random.default_rng(seed)
    acc = {key: [] for key in gams}
    nb = int(n_total // batch)
    t0 = time.time()
    for b in range(nb):
        u1, u2 = sample_pairs(batch, rng); d1, d2 = sample_pairs(batch, rng)
        pos = np.stack([u1, u2, d1, d2], axis=1)                  # (M,4,3)
        rA = np.linalg.norm(pos - CENTER_A, axis=-1); rB = np.linalg.norm(pos - CENTER_B, axis=-1)
        Vne = (-Z_A / rA - Z_B / rB).sum(axis=1)
        d = {ij: np.linalg.norm(pos[:, ij[0]] - pos[:, ij[1]], axis=-1) for ij in PAIRS}
        Vee = sum(1.0 / d[ij] for ij in PAIRS)
        for (gam, kind) in gams:
            if kind == 'exp':
                F = sum(np.exp(-gam * d[ij]) for ij in PAIRS)
            else:
                F = sum(d[ij] * (rA[:, ij[0]] - rA[:, ij[1]]) ** 2 * np.exp(-gam * d[ij]) for ij in PAIRS)
            Fm = F.mean()
            acc[(gam, kind)].append([Fm, F.var(), np.mean((F - Fm) * (Vne - Vne.mean())),
                                     np.mean((F - Fm) * (Vee - Vee.mean())), Vne.mean(), Vee.mean()])
        if b % 10 == 0:
            print(f"   MC batch {b+1}/{nb}  ({time.time()-t0:.0f}s)", flush=True)
    res = {}
    for key, rows in acc.items():
        R_ = np.array(rows); m = R_.mean(0); e = R_.std(0, ddof=1) / np.sqrt(len(R_))
        res[key] = dict(Fbar=(m[0], e[0]), sigma2=(m[1], e[1]), hVne=(m[2], e[2]),
                        hVee=(m[3], e[3]), Vne=(m[4], e[4]), E2=(m[5], e[5]))
    return res, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmc', type=float, default=1e8)
    ap.add_argument('--exact-neumann', action='store_true')
    ap.add_argument('--nphi', type=int, default=28)
    ap.add_argument('--tag', default='')
    args = ap.parse_args()
    KG.USE_EXACT_NEUMANN = bool(args.exact_neumann)
    NPHI = args.nphi
    CASES = [(0.5, 'exp'), (1.5, 'ut2'), (0.5, 'ut2')]
    out = {'nmc': args.nmc, 'exact_neumann': KG.USE_EXACT_NEUMANN, 'nphi': NPHI}
    ana = {}
    for gam, kind in CASES:
        r = analytic(gam, kind); ana[(gam, kind)] = r
        print(f"[analytic] {kind} gam={gam}: Fbar={r['Fbar']:.6f} sigma2={r['sigma2']:.6f} "
              f"hVne={r['hVne']:+.6f} hVee={r['hVee']:+.6f}  (E2={r['E2']:.5f}, Vne={r['Vne']:.5f}) "
              f"{r['seconds']:.0f}s", flush=True)
    res, tmc = mc(CASES, args.nmc)
    print(f"\nMC {args.nmc:.1e} configs, {tmc:.0f}s")
    print(f"{'case':14s} {'qty':7s} {'separable':>12s} {'MC':>12s} {'MC err':>9s} {'rel dev':>9s} {'n sigma':>8s}")
    table = []
    for gam, kind in CASES:
        for q in ('Fbar', 'sigma2', 'hVne', 'hVee', 'E2', 'Vne'):
            av = ana[(gam, kind)][q]; mv, me = res[(gam, kind)][q]
            rel = abs(av - mv) / abs(mv); ns = abs(av - mv) / me
            table.append(dict(gam=gam, kind=kind, qty=q, analytic=av, mc=mv, mc_err=me, rel=rel, nsig=ns))
            print(f"{kind+' g='+str(gam):14s} {q:7s} {av:12.6f} {mv:12.6f} {me:9.1e} {rel:9.1e} {ns:8.1f}")
    out['table'] = table
    out['analytic_parts'] = {f"{k[1]}_{k[0]}": {kk: vv for kk, vv in v.items() if kk.startswith('parts')}
                             for k, v in ana.items()}
    tag = ('_exactN' if KG.USE_EXACT_NEUMANN else '') + f'_nphi{NPHI}' + args.tag
    with open(os.path.join(ROOT, 'debug', 'data', f'r12t2_step1_separable{tag}.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)
    print("wrote debug/data/r12t2_step1_separable%s.json" % tag)
