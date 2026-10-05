r"""Step 3 diagnostic: core/valence pair-split geminals vs the orbital basis (VMC linear method).

Extends debug/ci_relax_step2_lm.py (same sampler, same guide density, same estimators).

Pair classes by a REAL-SPACE partition.  w(r) = 1/(1+(r_Li/r0)^n) is a smooth radial switch
around the Li nucleus (z = -R/2).  For a pair term f(r_ij) the three class weights are
    a_cc = w_i w_j,   a_cv = w_i + w_j - 2 w_i w_j,   a_vv = (1-w_i)(1-w_j),
with a_cc + a_cv + a_vv = 1 identically, so G_cc + G_cv + G_vv = F (same gamma) exactly.
Each weight is affine in w_i and in w_j separately, so
    grad_i (a f) = a_i' grad w_i f + a grad_i f
    lap_i  (a f) = a_i' lap w_i f + 2 a_i' grad w_i . grad_i f + a lap_i f.

Basis columns (one sampled S,H serves every ansatz -> exact correlated sampling):
    [0]        Phi_I                         N   (determinants of the truncated reference)
    [1]        F_{1.5} Phi_I                 N   (step-2 ansatz B block)
    [2..]      dressed blocks D_k Phi_I      N each, D_k in DRESS (class geminals at r0_d)
    [tail]     G^{r0}_{c,gamma} Psi_CI       one column per (r0, class, gamma)

mode 'lm'  : writes debug/data/pairgem_<tag>.npz (block S, H) + the column labels.
mode 'vmc' : plain VMC of Psi = v . X (v from an npz), standard local energy.
--ref 90det : the step-2 90-det pi reference; --ref m17 : NO-truncated M=17 cached wavefunction
             (debug/data/lih_vmc_wf_m17.pkl, natural orbitals of the full CI, |c|>trunc).
"""
from __future__ import annotations

import os
import sys
import time
import argparse
import pickle

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")

import geovac.lih_vmc as V
from geovac.lih_r12ci import pi_m16 as PM
import ci_relax_step2_lm as S2

GUIDE_NPZ = os.path.join(DATA, "lih_pi_m16_energy_exp1.5.npz")
GAMMAS = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
RADII = (0.6, 1.0, 1.5)
NSW = 6                      # switch steepness exponent
CLASSES = ("cc", "cv", "vv")
R0_D = 1.0
DRESS = (("cc", 3.0), ("cc", 5.0), ("ncc", 3.0))   # dressed per-det blocks at R0_D


# --------------------------------------------------------------------------------------------
def switch(r, R, r0, n=NSW):
    """w, grad w (Nw,Ne,3), lap w for w = 1/(1+(rho/r0)^n), rho = |r - Li|."""
    zc = np.array([0.0, 0.0, -R / 2.0])
    dv = r - zc
    rho = np.maximum(np.linalg.norm(dv, axis=-1), 1e-12)
    x = rho / r0; u = x ** n
    w = 1.0 / (1.0 + u)
    wp_over_r = -n * x ** (n - 2) / (r0 * r0 * (1.0 + u) ** 2)       # w'/rho (smooth at 0)
    wpp = (-n * (n - 1) * x ** (n - 2) / (1.0 + u) ** 2
           + 2.0 * n * n * x ** (2 * n - 2) / (1.0 + u) ** 3) / (r0 * r0)
    gw = wp_over_r[..., None] * dv
    lw = wpp + 2.0 * wp_over_r
    return w, gw, lw


def class_geminals(r, R, gammas, radii, need_derivs=True):
    """Returns dict (r0, cls, gamma) -> (G (Nw,), gradG (Nw,Ne,3), sum_i lap_i G (Nw,)).
    cls in cc/cv/vv; also 'ncc' = cv+vv and 'all' (= plain F, r0-independent check)."""
    Nw, Ne, _ = r.shape
    g = np.asarray(gammas)
    nG = len(g)
    # pair geometry, all gammas at once
    pairs = [(i, j) for i in range(Ne) for j in range(i + 1, Ne)]
    E = {}; FP = {}; LF = {}; U = {}
    for (i, j) in pairs:
        dv = r[:, i] - r[:, j]
        d = np.maximum(np.linalg.norm(dv, axis=-1), 1e-12)
        e = np.exp(-g[:, None] * d[None])                 # (nG,Nw)
        E[i, j] = e
        if need_derivs:
            FP[i, j] = -g[:, None] * e                      # f'
            LF[i, j] = g[:, None] ** 2 * e - 2.0 * g[:, None] * e / d[None]   # lap_i f (= lap_j f)
            U[i, j] = dv / d[:, None]
    out = {}
    for r0 in radii:
        w, gw, lw = switch(r, R, r0)
        acc = {c: [np.zeros((nG, Nw)), np.zeros((nG, Nw, Ne, 3)), np.zeros((nG, Nw))]
               for c in CLASSES}
        for (i, j) in pairs:
            wi, wj = w[:, i], w[:, j]
            A = {"cc": (wi * wj, wj, wi),
                 "cv": (wi + wj - 2 * wi * wj, 1 - 2 * wj, 1 - 2 * wi),
                 "vv": ((1 - wi) * (1 - wj), -(1 - wj), -(1 - wi))}
            e = E[i, j]
            for c in CLASSES:
                a, ai, aj = A[c]
                acc[c][0] += a[None] * e
                if need_derivs:
                    fp = FP[i, j]; u = U[i, j]
                    gif = fp[..., None] * u[None]           # grad_i f (nG,Nw,3); grad_j f = -gif
                    acc[c][1][:, :, i] += (ai[None, :, None] * gw[None, :, i] * e[..., None]
                                           + a[None, :, None] * gif)
                    acc[c][1][:, :, j] += (aj[None, :, None] * gw[None, :, j] * e[..., None]
                                           - a[None, :, None] * gif)
                    gwi_gif = np.einsum('wd,gwd->gw', gw[:, i], gif)
                    gwj_gjf = -np.einsum('wd,gwd->gw', gw[:, j], gif)
                    acc[c][2] += (ai[None] * lw[None, :, i] * e + 2 * ai[None] * gwi_gif
                                  + aj[None] * lw[None, :, j] * e + 2 * aj[None] * gwj_gjf
                                  + 2.0 * a[None] * LF[i, j])
        for c in CLASSES:
            for k, gam in enumerate(g):
                out[(r0, c, float(gam))] = tuple(x[k] for x in acc[c])
        for k, gam in enumerate(g):
            out[(r0, "ncc", float(gam))] = tuple(acc["cv"][m][k] + acc["vv"][m][k] for m in range(3))
    return out


# --------------------------------------------------------------------------------------------
def load_ref(name, trunc=1e-3):
    if name == "90det":
        z = np.load(S2.REFS["90det"], allow_pickle=True)
        wf = PM.build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
        return wf, float(z['E_trunc'])
    if name in ("m17", "m16pkl"):
        tag = "m17" if name == "m17" else "m16"
        with open(os.path.join(DATA, f"lih_vmc_wf_{tag}.pkl"), "rb") as f:
            d = pickle.load(f)
        wf = V._wf_from_dict(d) if isinstance(d, dict) else d
        return no_truncate(wf, trunc), float("nan")
    raise ValueError(name)


def no_truncate(wf, trunc):
    """Rotate the full-CI wavefunction to its natural orbitals (exact rotation of the CI tensor
    through 2x2 minors) and keep |c| > trunc, renormalized."""
    Mk = wf.Mk; ap = wf.apairs; bp = wf.bpairs
    na = len(ap)
    def antisym(vec):
        A = np.zeros((Mk, Mk))
        for k, (p, q) in enumerate(ap):
            A[p, q] = vec[k]; A[q, p] = -vec[k]
        return A
    CS = wf.CS
    ga = np.zeros((Mk, Mk)); gb = np.zeros((Mk, Mk))
    for b in range(CS.shape[1]):
        A = antisym(CS[:, b]); ga += A @ A.T
    for a in range(CS.shape[0]):
        B = antisym(CS[a, :]); gb += B @ B.T
    D = ga + gb
    occ, U = np.linalg.eigh(D)
    o = np.argsort(occ)[::-1]; occ = occ[o]; U = U[:, o]
    # chi_r = sum_k chi'_k U[r,k]  ->  |chi_p chi_q| = sum_{k<l} (U_pk U_ql - U_pl U_qk) |chi'_k chi'_l|
    W = np.zeros((na, na))
    for i, (p, q) in enumerate(ap):
        for j, (k, l) in enumerate(ap):
            W[i, j] = U[p, k] * U[q, l] - U[p, l] * U[q, k]
    CSn = W.T @ CS @ W
    print(f"  NO occupations: {np.round(occ, 5)}")
    print(f"  NO rotation: |CS|={np.linalg.norm(CS):.6f} |CSn|={np.linalg.norm(CSn):.6f}")
    keep = np.abs(CSn) > trunc
    nk = float(np.sum(CSn[keep] ** 2))
    CSt = np.where(keep, CSn, 0.0) / np.sqrt(nk)
    print(f"  trunc {trunc}: {int(keep.sum())} dets, kept norm {nk:.6f}")
    return V.LiHWavefunction(wf.prims, wf.T @ U, None, ap, bp, CSt, wf.R, wf.Z_A, wf.Z_B, wf.Vnn,
                             wf.E_fci, Mk, wf.na, wf.nb)


# --------------------------------------------------------------------------------------------
class Layout:
    def __init__(self, N, dress, radii, gammas):
        self.N = N
        self.labels = ["Phi"] * N + ["F1.5Phi"] * N
        for (c, gm) in dress:
            self.labels += [f"D_{c}{gm}_{R0_D}"] * N
        self.psi_cols = [(r0, c, gm) for r0 in radii for c in CLASSES for gm in gammas]
        self.labels += [f"P_{c}{gm}_{r0}" for (r0, c, gm) in self.psi_cols]
        self.K = len(self.labels)


def build_X(bas, r, lay, cvec, dress, radii, gammas, need_H=True):
    N = bas.N
    if need_H:
        Phi, gPhi, lPhi = bas.phi_full(r)
        Vp = V.potential(bas.wf, r)
        HPhi = Vp[None] * Phi - 0.5 * lPhi
    else:
        Phi = bas.phi_vals(r)
    F, gF, lF = S2.geminal(r, 1.5, need_derivs=need_H)
    gem = class_geminals(r, bas.wf.R, gammas, radii, need_derivs=need_H)
    Xs = [Phi, F[None] * Phi]; HXs = []
    if need_H:
        HXs = [HPhi, F[None] * HPhi - 0.5 * lF[None] * Phi - np.einsum('Wed,NWed->NW', gF, gPhi)]
    for (c, gm) in dress:
        G, gG, lG = gem[(R0_D, c, gm)]
        Xs.append(G[None] * Phi)
        if need_H:
            HXs.append(G[None] * HPhi - 0.5 * lG[None] * Phi - np.einsum('Wed,NWed->NW', gG, gPhi))
    psi = cvec @ Phi
    Xp = []; HXp = []
    if need_H:
        gpsi = np.einsum('N,NWed->Wed', cvec, gPhi); hpsi = cvec @ HPhi
    for key in lay.psi_cols:
        G, gG, lG = gem[key]
        Xp.append(G * psi)
        if need_H:
            HXp.append(G * hpsi - 0.5 * lG * psi - np.einsum('Wed,Wed->W', gG, gpsi))
    X = np.concatenate(Xs + ([np.array(Xp)] if Xp else []))
    if not need_H:
        return X, None
    HX = np.concatenate(HXs + ([np.array(HXp)] if HXp else []))
    return X, HX


def run_lm(bas, args, lay, c_g, Fbar_g, dress, radii, gammas):
    N = bas.N; eps = args.eps_mass / N; cvec = bas.cvec
    def logp(rr):
        Phi = bas.phi_vals(rr)
        psi = cvec @ Phi
        F, _, _ = S2.geminal(rr, 1.5, need_derivs=False)
        g = 1.0 + c_g * (F - Fbar_g)
        P = np.abs(g * psi) ** 2 + eps * np.sum(np.abs(Phi) ** 2, axis=0)
        return np.log(P + 1e-300)
    K = lay.K
    blocks_S, blocks_H = [], []
    accS = np.zeros((K, K), complex); accH = np.zeros_like(accS); cnt = [0]
    def measure(r):
        X, HX = build_X(bas, r, lay, cvec, dress, radii, gammas)
        Phi = X[:N]; psi = cvec @ Phi
        Fv, _, _ = S2.geminal(r, 1.5, need_derivs=False)
        g = 1.0 + c_g * (Fv - Fbar_g)
        P = np.abs(g * psi) ** 2 + eps * np.sum(np.abs(Phi) ** 2, axis=0)
        Xc = X.conj() * (1.0 / P)[None]
        accS[...] += Xc @ X.T / r.shape[0]; accH[...] += Xc @ HX.T / r.shape[0]
        cnt[0] += 1
        if cnt[0] == args.blk:
            blocks_S.append(accS / cnt[0]); blocks_H.append(accH / cnt[0])
            accS[...] = 0; accH[...] = 0; cnt[0] = 0
    acc, sd = S2.metropolis(bas, logp, args.nwalk, args.nsweep, args.nburn, args.seed, args.thin,
                            measure)
    return np.array(blocks_S), np.array(blocks_H), acc, sd


def run_vmc(bas, args, lay, v, dress, radii, gammas):
    cvec = bas.cvec
    nz = np.abs(v) > 0
    def logp(rr):
        X, _ = build_X(bas, rr, lay, cvec, dress, radii, gammas, need_H=False)
        return np.log(np.abs(v[nz] @ X[nz]) ** 2 + 1e-300)
    est = []
    def measure(r):
        X, HX = build_X(bas, r, lay, cvec, dress, radii, gammas)
        est.append(((v[nz] @ HX[nz]) / (v[nz] @ X[nz])).real.mean())
    acc, sd = S2.metropolis(bas, logp, args.nwalk, args.nsweep, args.nburn, args.seed, args.thin,
                            measure)
    return np.array(est), acc, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='lm', choices=['lm', 'vmc', 'check'])
    ap.add_argument('--ref', default='90det')
    ap.add_argument('--trunc', type=float, default=1e-3)
    ap.add_argument('--lite', action='store_true', help='m17 probe: only Phi, F1.5 Phi columns')
    ap.add_argument('--gammas', default=None, help='comma list overriding GAMMAS (Psi columns)')
    ap.add_argument('--radii', default=None, help='comma list overriding RADII')
    ap.add_argument('--nodress', action='store_true')
    ap.add_argument('--nwalk', type=int, default=1000)
    ap.add_argument('--nsweep', type=int, default=9000)
    ap.add_argument('--nburn', type=int, default=4000)
    ap.add_argument('--thin', type=int, default=5)
    ap.add_argument('--blk', type=int, default=40)
    ap.add_argument('--eps_mass', type=float, default=0.05)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--vec', default=None)
    ap.add_argument('--vkey', default='v')
    ap.add_argument('--tag', default=None)
    args = ap.parse_args()

    wf, E_trunc = load_ref(args.ref, args.trunc)
    bas = S2.Basis(wf)
    if args.lite:
        dress, radii, gammas = (), (), ()
    else:
        dress, radii, gammas = DRESS, RADII, GAMMAS
        if args.gammas: gammas = tuple(float(x) for x in args.gammas.split(','))
        if args.radii: radii = tuple(float(x) for x in args.radii.split(','))
        if args.nodress: dress = ()
    lay = Layout(bas.N, dress, radii, gammas)
    zg = np.load(GUIDE_NPZ, allow_pickle=True)
    c_g = float(zg['exp_c']); Fbar_g = float(zg['exp_Fbar'])
    print(f"[{args.mode}] ref={args.ref} N={bas.N} K={lay.K} seed={args.seed} nwalk={args.nwalk} "
          f"nsweep={args.nsweep} nburn={args.nburn}", flush=True)
    t0 = time.time()
    if args.mode == 'check':
        check(bas, lay, dress, radii, gammas)
        return
    if args.mode == 'lm':
        tag = args.tag or f"{args.ref}_lm_s{args.seed}"
        Sb, Hb, acc, sd = run_lm(bas, args, lay, c_g, Fbar_g, dress, radii, gammas)
        np.savez(os.path.join(DATA, f"pairgem_{tag}.npz"), S=Sb, H=Hb, cvec=bas.cvec, N=bas.N,
                 labels=np.array(lay.labels), psi_cols=np.array(lay.psi_cols, dtype=object),
                 dress=np.array(dress, dtype=object), acc=acc, step=sd, E_trunc=E_trunc,
                 nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn, seed=args.seed)
        print(f"[lm] done {time.time()-t0:.0f}s blocks={len(Sb)} acc={acc:.2f}", flush=True)
    else:
        tag = args.tag or f"{args.ref}_vmc_{args.vkey}_s{args.seed}"
        v = np.load(args.vec, allow_pickle=True)[args.vkey]
        est, acc, sd = run_vmc(bas, args, lay, v, dress, radii, gammas)
        nb = 20; m = len(est) // nb
        bl = est[:m * nb].reshape(nb, m).mean(1)
        E = bl.mean(); err = bl.std(ddof=1) / np.sqrt(nb)
        np.savez(os.path.join(DATA, f"pairgem_vmc_{tag}.npz"), est=est, E=E, err=err, acc=acc,
                 seed=args.seed)
        print(f"[vmc] {args.vkey}: E = {E:.6f} +/- {err:.6f}  ({time.time()-t0:.0f}s acc={acc:.2f})",
              flush=True)


def check(bas, lay, dress, radii, gammas):
    """Instrument checks: class sum == F; analytic grad/lap of class geminals vs finite diff;
    HX of a Psi column vs finite-difference Laplacian."""
    rng = np.random.default_rng(7)
    r = V._init_walkers(bas.wf, 6, rng)
    R = bas.wf.R
    gem = class_geminals(r, R, gammas, radii)
    for r0 in radii:
        for gm in gammas:
            s = sum(gem[(r0, c, gm)][0] for c in CLASSES)
            F, gF, lF = S2.geminal(r, gm)
            sg = sum(gem[(r0, c, gm)][1] for c in CLASSES); sl = sum(gem[(r0, c, gm)][2] for c in CLASSES)
            assert np.allclose(s, F, rtol=1e-12, atol=1e-14), "class sum"
            assert np.allclose(sg, gF, atol=1e-11) and np.allclose(sl, lF, atol=1e-9), "class sum derivs"
    print("  class partition sums to F (value, grad, lap): OK")
    h = 1e-4; worst_g = worst_l = 0.0
    for key in [(1.0, "cc", 3.0), (0.6, "cv", 5.0), (1.5, "vv", 1.0), (1.0, "ncc", 3.0)]:
        G, gG, lG = gem[key]
        lap_fd = np.zeros_like(G)
        for e in range(r.shape[1]):
            for dd in range(3):
                rp = r.copy(); rp[:, e, dd] += h; rm = r.copy(); rm[:, e, dd] -= h
                Gp = class_geminals(rp, R, gammas, radii, False)[key][0]
                Gm = class_geminals(rm, R, gammas, radii, False)[key][0]
                worst_g = max(worst_g, np.max(np.abs((Gp - Gm) / (2 * h) - gG[:, e, dd])))
                lap_fd += (Gp - 2 * G + Gm) / h ** 2
        worst_l = max(worst_l, np.max(np.abs(lap_fd - lG) / (np.abs(lG) + 1e-3)))
    print(f"  class geminal FD: max|grad err| = {worst_g:.2e}, max rel lap err = {worst_l:.2e}")
    # HX vs FD kinetic for one Psi column and one dressed column
    X, HX = build_X(bas, r, lay, bas.cvec, dress, radii, gammas)
    Vp = V.potential(bas.wf, r)
    for col in [lay.K - 1, lay.K - 30, 2 * bas.N + 5]:
        lap = np.zeros(r.shape[0], complex)
        for e in range(r.shape[1]):
            for dd in range(3):
                rp = r.copy(); rp[:, e, dd] += h; rm = r.copy(); rm[:, e, dd] -= h
                Xp, _ = build_X(bas, rp, lay, bas.cvec, dress, radii, gammas, need_H=False)
                Xm, _ = build_X(bas, rm, lay, bas.cvec, dress, radii, gammas, need_H=False)
                lap += (Xp[col] - 2 * X[col] + Xm[col]) / h ** 2
        Hfd = Vp * X[col] - 0.5 * lap
        rel = np.max(np.abs(Hfd - HX[col])) / np.max(np.abs(HX[col]))
        print(f"  H X[{col}] ({lay.labels[col]}) vs FD: rel err {rel:.2e}")


if __name__ == "__main__":
    main()
