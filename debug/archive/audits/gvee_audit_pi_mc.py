r"""g_Vee defect audit, step D: lean Metropolis MC of every <f_p f_q / r_r> on the pi trunc-1e-2
reference |Psi_CI|^2 (NO guide, no geminal; single-electron moves updating only the moved electron's
orbitals).  Independent of the analytic path (no grid / kernel / Neumann).  Per sweep it stores the
walker means of f_i, v_k, f_i v_k, f_i f_j v_k for each gamma, so any linear combination (e.g. the
centered g_Vee at a fixed Fbar, or a topology group) gets a block-jackknife error later.

    python debug/gvee_audit_pi_mc.py --gam 0.25 0.5 --nwalk 4000 --nsweep 6000 --nburn 1500 --seed 1
"""
from __future__ import annotations

import os
import sys
import time
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp')
    ap.add_argument('--gam', type=float, nargs='+', default=[0.25, 0.5])
    ap.add_argument('--nwalk', type=int, default=4000)
    ap.add_argument('--nsweep', type=int, default=6000)
    ap.add_argument('--nburn', type=int, default=1500)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    import lih_multigeminal_vmc as MG
    import lih_vmc
    from lih_pi_m16_vmc import _build_wf_from_no
    z = np.load(MG.REF_NPZ, allow_pickle=True)
    wf = _build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
    R = wf.R
    # restrict to the pairs / NOs that carry weight
    CS = np.asarray(wf.CS)
    ia = np.nonzero(np.any(CS != 0, axis=1))[0]; ib = np.nonzero(np.any(CS != 0, axis=0))[0]
    ap_ = [wf.apairs[k] for k in ia]; bp_ = [wf.bpairs[k] for k in ib]
    CSr = CS[np.ix_(ia, ib)]
    used = sorted(set(x for pq in ap_ + bp_ for x in pq))
    Tu = np.asarray(wf.T)[:, used]
    loc = {g: i for i, g in enumerate(used)}
    ap_l = [(loc[p], loc[q]) for p, q in ap_]; bp_l = [(loc[p], loc[q]) for p, q in bp_]
    print(f"  dets {np.count_nonzero(CSr)}  alpha pairs {len(ap_l)}  beta pairs {len(bp_l)}  NOs {used}",
          flush=True)
    assert np.allclose(CSr, CSr.T) if CSr.shape[0] == CSr.shape[1] else True

    def mo(re):                                     # (Nw,3) -> (nNO, Nw)
        Phi = lih_vmc.eval_values_batch(wf.prims, re, R)       # (Nprim, Nw)
        return np.einsum('ap,aW->pW', Tu, Phi)

    def dets(M1, M2, pl):
        return np.array([M1[p] * M2[q] - M1[q] * M2[p] for p, q in pl])   # (npair, Nw)

    rng = np.random.default_rng(5000 + a.seed)
    Nw = a.nwalk
    r = lih_vmc._init_walkers(wf, Nw, rng)
    MO = [mo(r[:, e]) for e in range(4)]
    dA = dets(MO[0], MO[1], ap_l); dB = dets(MO[2], MO[3], bp_l)
    psi = np.einsum('aW,ab,bW->W', dA, CSr, dB)
    lp = 2 * np.log(np.abs(psi) + 1e-300)
    sd = 0.45; acc = ntry = 0
    ng = len(a.gam)
    recs = []
    t0 = time.time()
    for sweep in range(a.nsweep):
        for e in range(4):
            rp = r[:, e] + rng.normal(0, sd, size=(Nw, 3))
            MOe = mo(rp)
            if e < 2:
                M1, M2 = (MOe, MO[1]) if e == 0 else (MO[0], MOe)
                dAp = dets(M1, M2, ap_l); psip = np.einsum('aW,ab,bW->W', dAp, CSr, dB)
            else:
                M1, M2 = (MOe, MO[3]) if e == 2 else (MO[2], MOe)
                dBp = dets(M1, M2, bp_l); psip = np.einsum('aW,ab,bW->W', dA, CSr, dBp)
            lpp = 2 * np.log(np.abs(psip) + 1e-300)
            A = np.log(rng.uniform(size=Nw)) < (lpp - lp)
            r[A, e] = rp[A]; MO[e][:, A] = MOe[:, A]; lp[A] = lpp[A]
            if e < 2:
                dA[:, A] = dAp[:, A]
            else:
                dB[:, A] = dBp[:, A]
            acc += int(A.sum()); ntry += Nw
        if sweep < a.nburn and sweep % 50 == 49:
            ar = acc / ntry
            sd = sd * 1.1 if ar > 0.6 else (sd / 1.1 if ar < 0.4 else sd)
            acc = ntry = 0
        if sweep >= a.nburn:
            d = np.stack([np.linalg.norm(r[:, i - 1] - r[:, j - 1], axis=1) for i, j in PAIRS6], axis=1)
            v = 1.0 / d
            rec = [v.mean(0)]
            for g in a.gam:
                f = d * np.exp(-g * d) if a.geminal == 'linexp' else np.exp(-g * d)
                rec += [f.mean(0), (np.einsum('ni,nk->ik', f, v) / Nw).ravel(),
                        (np.einsum('ni,nj,nk->ijk', f, f, v, optimize=True) / Nw).ravel()]
            recs.append(np.concatenate(rec))
        if sweep % 500 == 0:
            print(f"   sweep {sweep}/{a.nsweep} sd={sd:.3f} t={time.time()-t0:.0f}s", flush=True)
    X = np.array(recs)
    tag = f"{a.geminal}_g{'-'.join(f'{g:g}' for g in a.gam)}_s{a.seed}"
    np.savez(os.path.join(DATA, f"gvee_audit_pi_mc_{tag}.npz"), X=X, gam=np.array(a.gam), nwalk=Nw,
             acc=acc / max(ntry, 1))
    print(f"done {time.time()-t0:.0f}s nsamp={len(X)}", flush=True)


if __name__ == "__main__":
    main()
