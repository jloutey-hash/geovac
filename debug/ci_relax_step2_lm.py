r"""Step 2 diagnostic: relax the CI coefficients alongside the geminal (VMC linear method).

Basis (2N functions, N = kept determinants of the truncated M=16 pi reference):
    X_I       = Phi_I                    (I = 0..N-1)
    X_{N+I}   = F Phi_I,   F = sum_{i<j} exp(-gamma r_ij)   (gamma = 1.5)
Phi_I = dA[a_I] dB[b_I] (no sign; Psi_CI = sum_I CS[a_I,b_I] Phi_I).
Ansatz subspaces (all projections of ONE set of sampled matrices -> correlated sampling):
    (C) span{Psi_CI, F Psi_CI}                  2 params   (frozen CI; the banked 2x2)
    (A) span{Phi_I, F Psi_CI}                   N+1 params (CI free; geminal on frozen CI)
    (B) span{Phi_I, F Phi_I}                    2N params  (geminal-dressed copy per det)
(F - Fbar) Psi_CI and F Psi_CI span the same space once Psi_CI is in the span.

mode 'lm' : Metropolis on P = |g Psi_CI|^2 + eps sum_I |Phi_I|^2 (g = 1 + c(F-Fbar), the frozen
            marriage multiplier; the eps term bounds every ratio X/sqrt(P) so S,H have finite
            variance).  Accumulates block averages of S_kl = <X_k* X_l / P>, H_kl = <X_k* (H X_l)/P>
            (non-symmetric, standard-form local action).  Writes debug/data/cirelax_<tag>.npz.
mode 'vmc': plain VMC of Psi = v . X with coefficient vector v (from an npz), sampling |Psi|^2,
            standard local energy (H Psi)/Psi.  Writes debug/data/cirelax_vmc_<tag>.npz.
--nsweep counts TOTAL sweeps including --nburn.
"""
from __future__ import annotations

import os
import sys
import time
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

import geovac.lih_vmc as V
from geovac.lih_r12ci import pi_m16 as PM

REFS = {
    "9det": os.path.join(DATA, "lih_pi_m16_energy_trunc1e-2.npz"),
    "90det": os.path.join(DATA, "lih_pi_m16_energy_exp1.5.npz"),
}
GUIDE_NPZ = os.path.join(DATA, "lih_pi_m16_energy_exp1.5.npz")   # c, Fbar of the gamma=1.5 marriage
GAMMA = 1.5


class Basis:
    def __init__(self, wf):
        self.wf = wf
        ia, ib = np.nonzero(wf.CS)
        self.ia = ia; self.ib = ib; self.N = len(ia)
        self.cvec = wf.CS[ia, ib].astype(complex)            # Psi_CI in the Phi basis
        self.ua = sorted(set(ia.tolist())); self.ub = sorted(set(ib.tolist()))
        self.mapa = {a: k for k, a in enumerate(self.ua)}
        self.mapb = {b: k for k, b in enumerate(self.ub)}
        self.ka = np.array([self.mapa[a] for a in ia]); self.kb = np.array([self.mapb[b] for b in ib])

    # ---- MOs -------------------------------------------------------------------------------
    def _mo_vals(self, r):
        Phi = V.eval_values_batch(self.wf.prims, r, self.wf.R)
        return np.einsum('ap,aWe->pWe', self.wf.T, Phi)

    def _mo_vgl(self, r):
        Nw, Ne, _ = r.shape; M = len(self.wf.prims)
        Phi = np.empty((M, Nw, Ne), complex); gPhi = np.empty((M, Nw, Ne, 3), complex)
        lPhi = np.empty((M, Nw, Ne), complex)
        for a, p in enumerate(self.wf.prims):
            v, g, l = V.orbital_vgl(p, r, self.wf.R)
            Phi[a] = v; gPhi[a] = g; lPhi[a] = l
        T = self.wf.T
        return (np.einsum('ap,aWe->pWe', T, Phi), np.einsum('ap,aWed->pWed', T, gPhi),
                np.einsum('ap,aWe->pWe', T, lPhi))

    # ---- values only (Metropolis) ----------------------------------------------------------
    def phi_vals(self, r):
        MO = self._mo_vals(r)
        dA = np.array([MO[p, :, 0] * MO[q, :, 1] - MO[q, :, 0] * MO[p, :, 1]
                       for (p, q) in (self.wf.apairs[a] for a in self.ua)])
        dB = np.array([MO[p, :, 2] * MO[q, :, 3] - MO[q, :, 2] * MO[p, :, 3]
                       for (p, q) in (self.wf.bpairs[b] for b in self.ub)])
        return dA[self.ka] * dB[self.kb]                      # (N, Nw)

    # ---- full (measurement) ----------------------------------------------------------------
    def phi_full(self, r):
        """Phi (N,Nw), gradPhi (N,Nw,Ne,3), sum_i lap_i Phi (N,Nw)."""
        MO, gMO, lMO = self._mo_vgl(r)
        dA, gdA, ldA = V._minors_full(MO, gMO, lMO, [self.wf.apairs[a] for a in self.ua], [0, 1])
        dB, gdB, ldB = V._minors_full(MO, gMO, lMO, [self.wf.bpairs[b] for b in self.ub], [2, 3])
        A = dA[self.ka]; B = dB[self.kb]
        gA = gdA[:, self.ka]; gB = gdB[:, self.kb]            # (2,N,Nw,3)
        lA = ldA[:, self.ka].sum(0); lB = ldB[:, self.kb].sum(0)
        Nw = r.shape[0]
        G = np.empty((self.N, Nw, 4, 3), complex)
        G[:, :, 0] = gA[0] * B[:, :, None]; G[:, :, 1] = gA[1] * B[:, :, None]
        G[:, :, 2] = A[:, :, None] * gB[0]; G[:, :, 3] = A[:, :, None] * gB[1]
        L = lA * B + A * lB
        return A * B, G, L


def geminal(r, gamma=GAMMA, need_derivs=True):
    """F, grad_i F (Nw,Ne,3), sum_i lap_i F for f = exp(-gamma d)."""
    Nw, Ne, _ = r.shape
    F = np.zeros(Nw); G = np.zeros((Nw, Ne, 3)); L = np.zeros(Nw)
    for i in range(Ne):
        for j in range(i + 1, Ne):
            dv = r[:, i] - r[:, j]
            d = np.maximum(np.linalg.norm(dv, axis=-1), 1e-12)
            e = np.exp(-gamma * d)
            F += e
            if need_derivs:
                fp = -gamma * e
                u = dv / d[:, None]
                G[:, i] += fp[:, None] * u; G[:, j] -= fp[:, None] * u
                L += 2.0 * (gamma * gamma * e - 2.0 * gamma * e / d)
    return F, G, L


def basis_X_HX(bas, r):
    """X (2N,Nw), HX (2N,Nw) complex."""
    Phi, gPhi, lPhi = bas.phi_full(r)
    Vp = V.potential(bas.wf, r)
    F, gF, lF = geminal(r)
    HPhi = Vp[None] * Phi - 0.5 * lPhi
    cross = np.einsum('Wed,NWed->NW', gF, gPhi)
    HFPhi = F[None] * HPhi - 0.5 * lF[None] * Phi - cross
    X = np.concatenate([Phi, F[None] * Phi]); HX = np.concatenate([HPhi, HFPhi])
    return X, HX


def metropolis(bas, logp, nwalk, nsweep, nburn, seed, thin, measure, step=0.45):
    rng = np.random.default_rng(seed)
    wf = bas.wf; Ne = wf.na + wf.nb
    r = V._init_walkers(wf, nwalk, rng)
    lp = logp(r)
    acc = ntry = 0; sd = step; t0 = time.time()
    for sweep in range(nsweep):
        for i in range(Ne):
            rp = r.copy(); rp[:, i] += rng.normal(0, sd, size=(nwalk, 3))
            lpp = logp(rp)
            a = np.log(rng.uniform(size=nwalk)) < (lpp - lp)
            r[a] = rp[a]; lp[a] = lpp[a]
            acc += int(a.sum()); ntry += nwalk
        if sweep < nburn and sweep % 50 == 49:
            ar = acc / max(ntry, 1)
            if ar > 0.6: sd *= 1.1
            elif ar < 0.4: sd /= 1.1
            acc = ntry = 0
        if sweep >= nburn and (sweep - nburn) % thin == 0:
            measure(r)
        if sweep % 1000 == 0:
            print(f"   sweep {sweep}/{nsweep} step={sd:.3f} t={time.time()-t0:.0f}s", flush=True)
    return acc / max(ntry, 1), sd


def run_lm(bas, args, c_g, Fbar_g):
    N = bas.N; eps = args.eps_mass / N
    def logp(rr):
        Phi = bas.phi_vals(rr)
        psi = bas.cvec @ Phi
        F, _, _ = geminal(rr, need_derivs=False)
        g = 1.0 + c_g * (F - Fbar_g)
        P = np.abs(g * psi) ** 2 + eps * np.sum(np.abs(Phi) ** 2, axis=0)
        return np.log(P + 1e-300)
    blocks_S, blocks_H = [], []
    accS = np.zeros((2 * N, 2 * N), complex); accH = np.zeros_like(accS); cnt = [0]
    def measure(r):
        X, HX = basis_X_HX(bas, r)
        Phi = X[:N]; psi = bas.cvec @ Phi
        Fv, _, _ = geminal(r, need_derivs=False)
        g = 1.0 + c_g * (Fv - Fbar_g)
        P = np.abs(g * psi) ** 2 + eps * np.sum(np.abs(Phi) ** 2, axis=0)
        w = 1.0 / P
        Xc = X.conj() * w[None]
        accS[...] += Xc @ X.T / r.shape[0]; accH[...] += Xc @ HX.T / r.shape[0]
        cnt[0] += 1
        if cnt[0] == args.blk:
            blocks_S.append(accS / cnt[0]); blocks_H.append(accH / cnt[0])
            accS[...] = 0; accH[...] = 0; cnt[0] = 0
    acc, sd = metropolis(bas, logp, args.nwalk, args.nsweep, args.nburn, args.seed, args.thin, measure)
    return np.array(blocks_S), np.array(blocks_H), acc, sd


def run_vmc(bas, args, v):
    N = bas.N
    def logp(rr):
        Phi = bas.phi_vals(rr)
        F, _, _ = geminal(rr, need_derivs=False)
        psi = v[:N] @ Phi + v[N:] @ (F[None] * Phi)
        return np.log(np.abs(psi) ** 2 + 1e-300)
    est = []
    def measure(r):
        X, HX = basis_X_HX(bas, r)
        EL = ((v @ HX) / (v @ X)).real
        est.append(EL.mean())
    acc, sd = metropolis(bas, logp, args.nwalk, args.nsweep, args.nburn, args.seed, args.thin, measure)
    return np.array(est), acc, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='lm', choices=['lm', 'vmc'])
    ap.add_argument('--ref', default='9det', choices=list(REFS))
    ap.add_argument('--nwalk', type=int, default=1000)
    ap.add_argument('--nsweep', type=int, default=6000)
    ap.add_argument('--nburn', type=int, default=3000)
    ap.add_argument('--thin', type=int, default=5)
    ap.add_argument('--blk', type=int, default=10, help='measurements per stored block (lm)')
    ap.add_argument('--eps_mass', type=float, default=0.05)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--vec', default=None, help='npz with key v (vmc mode)')
    ap.add_argument('--vkey', default='v')
    ap.add_argument('--tag', default=None)
    args = ap.parse_args()

    z = np.load(REFS[args.ref], allow_pickle=True)
    wf = PM.build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
    bas = Basis(wf)
    zg = np.load(GUIDE_NPZ, allow_pickle=True)
    c_g = float(zg['exp_c']); Fbar_g = float(zg['exp_Fbar'])
    t0 = time.time()
    print(f"[{args.mode}] ref={args.ref} N={bas.N} seed={args.seed} nwalk={args.nwalk} "
          f"nsweep={args.nsweep} nburn={args.nburn} E_trunc={float(z['E_trunc']):.6f}", flush=True)
    if args.mode == 'lm':
        tag = args.tag or f"{args.ref}_lm_s{args.seed}"
        Sb, Hb, acc, sd = run_lm(bas, args, c_g, Fbar_g)
        np.savez(os.path.join(DATA, f"cirelax_{tag}.npz"), S=Sb, H=Hb, cvec=bas.cvec, N=bas.N,
                 ia=bas.ia, ib=bas.ib, acc=acc, step=sd, E_trunc=float(z['E_trunc']),
                 nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn, seed=args.seed)
        print(f"[lm] done {time.time()-t0:.0f}s blocks={len(Sb)} acc={acc:.2f}", flush=True)
    else:
        tag = args.tag or f"{args.ref}_vmc_{args.vkey}_s{args.seed}"
        v = np.load(args.vec, allow_pickle=True)[args.vkey]
        est, acc, sd = run_vmc(bas, args, v)
        nb = 20; m = len(est) // nb
        bl = est[:m * nb].reshape(nb, m).mean(1)
        E = bl.mean(); err = bl.std(ddof=1) / np.sqrt(nb)
        np.savez(os.path.join(DATA, f"cirelax_vmc_{tag}.npz"), est=est, E=E, err=err, acc=acc,
                 nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn, seed=args.seed)
        print(f"[vmc] {args.vkey}: E = {E:.6f} +/- {err:.6f}  ({time.time()-t0:.0f}s acc={acc:.2f})",
              flush=True)


if __name__ == "__main__":
    main()
