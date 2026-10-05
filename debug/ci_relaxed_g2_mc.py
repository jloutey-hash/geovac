r"""Gate G2 for geovac/lih_r12ci/ci_relaxed.py: the transition columns vs direct real-space MC.

Real-space orbitals from geovac.lih_vmc (analytic primitives -> NOs, complex pi), the same basis as
the step-2 LM diagnostic: chi_I = dA[a_I] dB[b_I] / 2, psi = sum_I d_I chi_I.  Metropolis on
    P = |psi|^2 + (1/N) sum_I |chi_I|^2           (normalisation Z = 2 exactly: orthonormal NOs)
and the unbiased estimators  <chi_I|O|psi> = Z < chi_I^* (O psi) / P >  for
    O psi in { psi, F psi, T(F psi), V F psi (V = V_ne + V_ee, no V_NN), H psi }
with T(F psi) = F T psi - 1/2 (sum_i lap_i F) psi - sum_i grad_i F . grad_i psi  (analytic Laplacians).
F = sum_{i<j} exp(-1.5 r_ij).  Writes debug/data/ci_relaxed_g2_<ref>_s<seed>.npz (block means).

Usage: python debug/ci_relaxed_g2_mc.py 1e-2 --seed 1 --nwalk 1000 --nsweep 4000 --nburn 1000
"""
from __future__ import annotations

import argparse
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")
GAMMA = 1.5


def geminal(r, gamma=GAMMA):
    Nw, Ne, _ = r.shape
    F = np.zeros(Nw); G = np.zeros((Nw, Ne, 3)); L = np.zeros(Nw)
    for i in range(Ne):
        for j in range(i + 1, Ne):
            dv = r[:, i] - r[:, j]
            d = np.maximum(np.linalg.norm(dv, axis=-1), 1e-12)
            e = np.exp(-gamma * d)
            F += e
            fp = -gamma * e
            u = dv / d[:, None]
            G[:, i] += fp[:, None] * u; G[:, j] -= fp[:, None] * u
            L += 2.0 * (gamma * gamma * e - 2.0 * gamma * e / d)
    return F, G, L


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trunc", type=float)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nwalk", type=int, default=1000)
    ap.add_argument("--nsweep", type=int, default=4000, help="total sweeps incl. burn-in")
    ap.add_argument("--nburn", type=int, default=1000)
    ap.add_argument("--thin", type=int, default=5)
    ap.add_argument("--blk", type=int, default=10)
    a = ap.parse_args()
    import geovac.lih_vmc as V
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.ref_tensor as P0
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=a.trunc,
                                        verbose=False)
    wf = M.build_wf_from_no(info['ct'], info['T_no'], int(info['Mk']), float(info['E_trunc']))
    CS = wf.CS
    ia, ib = np.nonzero(CS); d = CS[ia, ib]; N = len(d)
    ua = sorted(set(ia.tolist())); ub = sorted(set(ib.tolist()))
    ka = np.array([ua.index(x) for x in ia]); kb = np.array([ub.index(x) for x in ib])
    apl = [wf.apairs[x] for x in ua]; bpl = [wf.bpairs[x] for x in ub]
    T = wf.T; prims = wf.prims; R = wf.R

    def mo_vals(r):
        Phi = V.eval_values_batch(prims, r, R)
        return np.einsum('ap,aWe->pWe', T, Phi)

    def chi_vals(r):
        MO = mo_vals(r)
        dA = np.array([MO[p, :, 0] * MO[q, :, 1] - MO[q, :, 0] * MO[p, :, 1] for (p, q) in apl])
        dB = np.array([MO[p, :, 2] * MO[q, :, 3] - MO[q, :, 2] * MO[p, :, 3] for (p, q) in bpl])
        return 0.5 * dA[ka] * dB[kb]

    def chi_full(r):
        Mn = len(prims); Nw = r.shape[0]
        Ph = np.empty((Mn, Nw, 4), complex); gP = np.empty((Mn, Nw, 4, 3), complex)
        lP = np.empty((Mn, Nw, 4), complex)
        for k, p in enumerate(prims):
            v, g, l = V.orbital_vgl(p, r, R)
            Ph[k] = v; gP[k] = g; lP[k] = l
        MO = np.einsum('ap,aWe->pWe', T, Ph); gMO = np.einsum('ap,aWed->pWed', T, gP)
        lMO = np.einsum('ap,aWe->pWe', T, lP)
        dA, gdA, ldA = V._minors_full(MO, gMO, lMO, apl, [0, 1])
        dB, gdB, ldB = V._minors_full(MO, gMO, lMO, bpl, [2, 3])
        A = dA[ka]; B = dB[kb]
        G = np.empty((N, Nw, 4, 3), complex)
        G[:, :, 0] = gdA[0][ka] * B[:, :, None]; G[:, :, 1] = gdA[1][ka] * B[:, :, None]
        G[:, :, 2] = A[:, :, None] * gdB[0][kb]; G[:, :, 3] = A[:, :, None] * gdB[1][kb]
        L = ldA[:, ka].sum(0) * B + A * ldB[:, kb].sum(0)
        return 0.5 * A * B, 0.5 * G, 0.5 * L

    Z = 2.0
    eps = 1.0 / N

    def logp(r):
        X = chi_vals(r)
        psi = d @ X
        return np.log(np.abs(psi) ** 2 + eps * np.sum(np.abs(X) ** 2, axis=0) + 1e-300)

    names = ['S0', 'SF', 'KF', 'VF', 'Hpsi', 'Tpsi']
    acc_blk = {k: np.zeros(N) for k in names}
    blocks = {k: [] for k in names}
    cnt = [0]

    def measure(r):
        X, gX, lX = chi_full(r)
        psi = d @ X; gpsi = np.einsum('I,IWed->Wed', d, gX); lpsi = d @ lX
        P = np.abs(psi) ** 2 + eps * np.sum(np.abs(X) ** 2, axis=0)
        F, gF, lF = geminal(r)
        Vp = V.potential(wf, r) - wf.Vnn
        Tpsi = -0.5 * lpsi
        TF = F * Tpsi - 0.5 * lF * psi - np.einsum('Wed,Wed->W', gF, gpsi)
        Xc = X.conj() * (Z / P)[None]
        vals = dict(S0=psi, SF=F * psi, KF=TF, VF=Vp * F * psi, Hpsi=Tpsi + Vp * psi, Tpsi=Tpsi)
        for k in names:
            acc_blk[k] += (Xc @ vals[k]).real / r.shape[0]
        cnt[0] += 1
        if cnt[0] == a.blk:
            for k in names:
                blocks[k].append(acc_blk[k] / cnt[0]); acc_blk[k][:] = 0
            cnt[0] = 0

    rng = np.random.default_rng(a.seed)
    r = V._init_walkers(wf, a.nwalk, rng)
    lp = logp(r); sd = 0.45; acc = ntry = 0; t0 = time.time()
    for sweep in range(a.nsweep):
        for i in range(4):
            rp = r.copy(); rp[:, i] += rng.normal(0, sd, size=(a.nwalk, 3))
            lpp = logp(rp)
            ok = np.log(rng.uniform(size=a.nwalk)) < (lpp - lp)
            r[ok] = rp[ok]; lp[ok] = lpp[ok]
            acc += int(ok.sum()); ntry += a.nwalk
        if sweep < a.nburn and sweep % 50 == 49:
            ar = acc / max(ntry, 1)
            sd = sd * 1.1 if ar > 0.6 else (sd / 1.1 if ar < 0.4 else sd)
            acc = ntry = 0
        if sweep >= a.nburn and (sweep - a.nburn) % a.thin == 0:
            measure(r)
        if sweep % 1000 == 0:
            print(f"  sweep {sweep}/{a.nsweep} step={sd:.3f} t={time.time()-t0:.0f}s", flush=True)
    tag = f"{a.trunc:.0e}"
    out = os.path.join(DATA, f"ci_relaxed_g2_{tag}_s{a.seed}.npz")
    np.savez(out, d=d, ia=ia, ib=ib, **{k: np.array(v) for k, v in blocks.items()},
             nwalk=a.nwalk, nsweep=a.nsweep, nburn=a.nburn, seed=a.seed, acc=acc / max(ntry, 1))
    print(f"done {time.time()-t0:.0f}s blocks={len(blocks['SF'])} -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
