r"""Job 2 (c-offset) v2: per-COMPONENT VMC, basis {Psi, F_lx0.25 Psi, F_lx0.5 Psi}, with BOTH kinetic estimators (non-symmetric -1/2 lap and symmetric 1/2 |grad|^2); derived from lih_c_offset_vmc.py (same guide). on the trunc-1e-2 pi reference.

Guide |Psi_CI (1 + c (F - Fbar))|^2 (c, Fbar = the analytic values, as lih_multigeminal_vmc).
Per sample (walker-averaged, weight w = 1/g^2) stores, in the basis phi0 = Psi_CI, phi1 = F Psi_CI:
    S (2,2) = <w A_a A_b>
    H^X (2,2) = <w A_a  H_X(A_b Psi)/Psi>   for X in {T, ne, ee}  (V_NN is a constant, added later)
  H_T(A Psi)/Psi = A T_L0 - 1/2 lap A - grad A . grad ln Psi,   T_L0 = -1/2 lap Psi / Psi.
Then g_X = [H^X_11 - Fbar (H^X_01 + H^X_10) + Fbar^2 H^X_00]/S00  and
h_X(col) = (H^X_10 - Fbar H^X_00)/S00, h_X(row) = (H^X_01 - Fbar H^X_00)/S00.
--nsweep is TOTAL sweeps including --nburn.  Writes debug/data/cvmc_s<seed>.npz.
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
sys.path.insert(0, HERE)
import lih_multigeminal_vmc as MG          # namespace stub + _f + REF_NPZ
import lih_vmc
from lih_pi_m16_vmc import _build_wf_from_no

DATA = MG.DATA
KIND, GAM = "linexp", 0.5
BASIS = [("linexp", 0.25), ("linexp", 0.5)]


def F_terms(r, kind=KIND, gam=GAM):
    Nw, Ne, _ = r.shape
    A = np.zeros(Nw); G = np.zeros((Nw, Ne, 3)); L = np.zeros(Nw)
    for i in range(Ne):
        for j in range(i + 1, Ne):
            dvec = r[:, i] - r[:, j]
            d = np.maximum(np.linalg.norm(dvec, axis=-1), 1e-12)
            u = dvec / d[:, None]
            f, fp, fl = MG._f(kind, gam, d)
            A += f; G[:, i] += fp[:, None] * u; G[:, j] -= fp[:, None] * u; L += 2.0 * fl
    return A, G, L


def split_potential(wf, r):
    R = wf.R; Nw, Ne, _ = r.shape
    zc = np.array([0.0, 0.0, -R / 2.0]); zh = np.array([0.0, 0.0, R / 2.0])
    rA = np.linalg.norm(r - zc, axis=-1); rB = np.linalg.norm(r - zh, axis=-1)
    Vne = -np.sum(wf.Z_A / rA + wf.Z_B / rB, axis=1)
    Vee = np.zeros(Nw)
    for i in range(Ne):
        for j in range(i + 1, Ne):
            Vee += 1.0 / np.linalg.norm(r[:, i] - r[:, j], axis=-1)
    return Vne, Vee


def sample(wf, r, cg, Fb):
    psi, gpsi, lpsi = lih_vmc.psi_ci_full(wf, r)
    glnpsi = gpsi / psi[:, None, None]
    TL0 = (-0.5 * np.sum(lpsi, axis=1) / psi)
    Vne, Vee = split_potential(wf, r)
    Nw = r.shape[0]
    As, Gs, HTns = [np.ones(Nw)], [np.zeros((Nw,) + r.shape[1:])], [TL0.real]
    for kind, gam in BASIS:
        F, gF, lF = F_terms(r, kind, gam)
        As.append(F); Gs.append(gF)
        HTns.append((F * TL0 - 0.5 * lF - np.einsum('Wed,Wed->W', gF, glnpsi)).real)
    A = np.array(As); HT = np.array(HTns)
    # symmetric kinetic: 1/2 sum_i Re[conj(grad A_a + A_a grad lnPsi) . (grad A_b + A_b grad lnPsi)]
    V = np.array([Gs[k] + As[k][:, None, None] * glnpsi for k in range(len(As))])   # complex (nb,Nw,Ne,3)
    Tsym = 0.5 * np.einsum('aWed,bWed->abW', V.conj(), V).real
    Fg = F_terms(r)[0]
    w = 1.0 / (1.0 + cg * (Fg - Fb)) ** 2
    S = np.einsum('aW,bW,W->ab', A, A, w) / Nw
    out = [S, np.einsum('aW,bW,W->ab', A, HT, w) / Nw, np.einsum('abW,W->ab', Tsym, w) / Nw,
           np.einsum('aW,bW,W->ab', A, A * Vne[None, :], w) / Nw,
           np.einsum('aW,bW,W->ab', A, A * Vee[None, :], w) / Nw]
    return np.array(out)                                        # (5,3,3): S, T_ns, T_sym, Vne, Vee


def run(wf, cg, Fb, nwalk, nsweep, nburn, seed, thin, step=0.45):
    rng = np.random.default_rng(seed)
    Ne = wf.na + wf.nb
    r = lih_vmc._init_walkers(wf, nwalk, rng)

    def logp(rr):
        psi = lih_vmc.psi_ci_val(wf, rr)
        F = F_terms(rr)[0]
        g = 1.0 + cg * (F - Fb)
        return 2.0 * np.log(np.abs(psi) + 1e-300) + 2.0 * np.log(np.abs(g) + 1e-300)

    lp = logp(r); acc = ntry = 0; sd = step; out = []; t0 = time.time()
    for sweep in range(nsweep):
        for i in range(Ne):
            rp = r.copy(); rp[:, i] += rng.normal(0, sd, size=(nwalk, 3))
            lpp = logp(rp)
            a = np.log(rng.uniform(size=nwalk)) < (lpp - lp)
            r[a] = rp[a]; lp[a] = lpp[a]; acc += int(a.sum()); ntry += nwalk
        if sweep < nburn and sweep % 50 == 49:
            ar = acc / max(ntry, 1)
            if ar > 0.6: sd *= 1.1
            elif ar < 0.4: sd /= 1.1
            acc = ntry = 0
        if sweep >= nburn and (sweep - nburn) % thin == 0:
            out.append(sample(wf, r, cg, Fb))
        if sweep % 1000 == 0:
            print(f"   sweep {sweep}/{nsweep} step={sd:.3f} t={time.time()-t0:.0f}s", flush=True)
    return np.array(out), acc / max(ntry, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nwalk', type=int, default=1000)
    ap.add_argument('--nsweep', type=int, default=9000, help='TOTAL incl. burn-in')
    ap.add_argument('--nburn', type=int, default=4000)
    ap.add_argument('--thin', type=int, default=5)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    assert a.nsweep > a.nburn
    z = np.load(MG.REF_NPZ, allow_pickle=True)
    cg = float(z['linexp_c']); Fb = float(z['linexp_Fbar'])
    wf = _build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
    t0 = time.time()
    X, acc = run(wf, cg, Fb, a.nwalk, a.nsweep, a.nburn, a.seed, a.thin)
    np.savez(os.path.join(DATA, f"cvmc2_s{a.seed}.npz"), X=X, seed=a.seed, nwalk=a.nwalk,
             nsweep=a.nsweep, nburn=a.nburn, acc=acc, c_guide=cg, Fbar_guide=Fb, Vnn=wf.Vnn)
    print(f"[cvmc] seed {a.seed} done {time.time()-t0:.0f}s acc={acc:.2f} nsamp={len(X)}", flush=True)


if __name__ == "__main__":
    main()
