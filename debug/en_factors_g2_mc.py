r"""Gate G2 for geovac/lih_r12ci/en_factors.py: the new e-n matrix elements vs direct real-space MC.

Real-space orbitals from geovac.lih_vmc (as debug/ci_relaxed_g2_mc.py): chi_I = dA dB / 2,
psi = sum_I d_I chi_I.  Metropolis on P = |psi|^2 + (1/N) sum_I |chi_I|^2 (Z = 2 exactly) and the
unbiased estimators <A|B> = Z < A^* B / P >.  Columns X_k = A_k psi for
    A in [1, F, W1, W2, Fcc, Fc1]
(F = sum exp(-1.5 r_ij); W1 = sum_i w_i; W2 = sum_{i<j} w_i w_j; Fcc = sum w_i w_j f_ij;
Fc1 = sum (w_i + w_j) f_ij), w = the parameter-free core switch (core_switch.CoreSwitch, core-NO
density).  H(A psi) = A H psi - 1/2 (sum_i lap_i A) psi - sum_i grad_i A . grad_i psi, all analytic.
Electronic H (V_NN removed).  Saved block means:
    tS[I,k] = <chi_I|X_k>,  tH[I,k] = <chi_I|H X_k>,  S[k,l] = <X_k|X_l>,  H[k,l] = <X_k|H X_l>,
    and the V_ee-only parts tV[I,k] = <chi_I|V_ee X_k>, V[k,l] = <X_k|V_ee X_l> (low variance), the V_ne part
    Vn[k,l], and the kinetic part in GRADIENT form Tg[k,l] = 1/2 sum_i <grad_i X_k|grad_i X_l> (low variance).
Usage: python debug/en_factors_g2_mc.py 1e-2 --seed 1 [--nwalk 1000 --nsweep 6000 --nburn 1500]
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
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")

import geovac.lih_vmc as V  # noqa: E402
import ci_relax_step2_lm as S2  # noqa: E402
import en_factors_stage0_lm as E0  # noqa: E402
from geovac.lih_r12ci import pi_m16 as PM  # noqa: E402
from geovac.lih_r12ci.core_switch import CoreSwitch  # noqa: E402

REFS = {1e-2: S2.REFS["9det"], 1e-3: S2.REFS["90det"]}
COLS = ["1", "F", "W1", "W2", "Fcc", "Fc1"]
GAMMA = 1.5


def en_ops(w, gw, lw):
    """W1, W2 with per-electron gradients and summed Laplacians.  w (Nw,Ne), gw (Nw,Ne,3), lw (Nw,Ne)."""
    W1 = w.sum(1)
    W2 = 0.5 * (W1 ** 2 - (w ** 2).sum(1))
    rest = W1[:, None] - w                                     # sum_{j != i} w_j
    return [(W1, gw.copy(), lw.sum(1)), (W2, gw * rest[..., None], (lw * rest).sum(1))]


def weighted_geminals(r, w, gw, lw, gamma=GAMMA):
    """(G, grad G, sum lap G) for Fcc = sum w_i w_j f_ij and Fc1 = sum (w_i + w_j) f_ij."""
    Nw, Ne, _ = r.shape
    out = []
    for kind in ("cc", "c1"):
        G = np.zeros(Nw); gG = np.zeros((Nw, Ne, 3)); lG = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                dv = r[:, i] - r[:, j]
                d = np.maximum(np.linalg.norm(dv, axis=-1), 1e-12)
                f = np.exp(-gamma * d); u = dv / d[:, None]
                gif = (-gamma * f)[:, None] * u                     # grad_i f; grad_j f = -gif
                lf = gamma * gamma * f - 2.0 * gamma * f / d        # lap_i f = lap_j f
                if kind == "cc":
                    a = w[:, i] * w[:, j]; dai = w[:, j]; daj = w[:, i]
                else:
                    a = w[:, i] + w[:, j]; dai = np.ones(Nw); daj = np.ones(Nw)
                G += a * f
                gG[:, i] += dai[:, None] * gw[:, i] * f[:, None] + a[:, None] * gif
                gG[:, j] += daj[:, None] * gw[:, j] * f[:, None] - a[:, None] * gif
                lG += (dai * lw[:, i] + daj * lw[:, j]) * f \
                    + 2.0 * dai * np.einsum('wd,wd->w', gw[:, i], gif) \
                    - 2.0 * daj * np.einsum('wd,wd->w', gw[:, j], gif) + 2.0 * a * lf
        out.append((G, gG, lG))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trunc", type=float)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nwalk", type=int, default=1000)
    ap.add_argument("--nsweep", type=int, default=6000)
    ap.add_argument("--nburn", type=int, default=1500)
    ap.add_argument("--thin", type=int, default=5)
    ap.add_argument("--blk", type=int, default=20)
    ap.add_argument("--check", action="store_true", help="finite-difference check of the columns")
    ap.add_argument("--vmc-npz", default=None, help="analytic npz holding the eigenvector (VMC mode)")
    ap.add_argument("--vmc-key", default="vec_A_W12")
    a = ap.parse_args()
    z = np.load(REFS[a.trunc], allow_pickle=True)
    wf = PM.build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
    bas = S2.Basis(wf)
    T_no, occ, m_no, R = E0.switch_reference()
    sw = CoreSwitch.from_reference(wf.prims, T_no, occ, m_no, R)
    N = bas.N; cvec = bas.cvec.real
    Z = 2.0
    eps = 1.0 / N
    K = len(COLS)

    def columns(r, need_H=True):
        if need_H:
            Phi, gPhi, lPhi = bas.phi_full(r)
        else:
            Phi = bas.phi_vals(r)
        chi = Phi / 2.0
        psi = cvec @ chi
        F, gF, lF = S2.geminal(r, GAMMA, need_derivs=need_H)
        if need_H:
            w, gw, lw = sw.vgl(r)
        else:
            w = sw.value(r); gw = np.zeros(r.shape); lw = np.zeros(w.shape)
        mults = [(np.ones(r.shape[0]), np.zeros(r.shape), np.zeros(r.shape[0])), (F, gF, lF)]
        mults += en_ops(w, gw, lw) + weighted_geminals(r, w, gw, lw)
        X = np.array([m[0] * psi for m in mults])
        if not need_H:
            return chi, X, None, None, None, None
        gchi = gPhi / 2.0; lchi = lPhi / 2.0
        Vp = V.potential(wf, r) - wf.Vnn
        gpsi = np.einsum('N,NWed->Wed', cvec, gchi); lpsi = cvec @ lchi
        Hpsi = Vp * psi - 0.5 * lpsi
        HX = np.array([A * Hpsi - 0.5 * lA * psi - np.einsum('Wed,Wed->W', gA, gpsi) for (A, gA, lA) in mults])
        Vee = np.zeros(r.shape[0])
        for i in range(r.shape[1]):
            for j in range(i + 1, r.shape[1]):
                Vee += 1.0 / np.linalg.norm(r[:, i] - r[:, j], axis=-1)
        GX = np.array([gA * psi[:, None, None] + A[:, None, None] * gpsi for (A, gA, _lA) in mults])
        Vne = Vp - Vee
        return chi, X, HX, Vee[None] * X, GX, Vne[None] * X

    if a.check:
        rng = np.random.default_rng(7)
        r = V._init_walkers(wf, 5, rng)
        chi, X, HX, _VX, _GX, _VnX = columns(r)
        Vp = V.potential(wf, r) - wf.Vnn
        h = 1e-4; lap = np.zeros_like(X)
        for e in range(r.shape[1]):
            for dd in range(3):
                rp = r.copy(); rp[:, e, dd] += h; rm = r.copy(); rm[:, e, dd] -= h
                lap += (columns(rp, False)[1] - 2 * X + columns(rm, False)[1]) / h ** 2
        Hfd = Vp[None] * X - 0.5 * lap
        for k in range(K):
            print(f"  H X[{COLS[k]}] vs FD: rel {np.max(np.abs(Hfd[k] - HX[k])) / np.max(np.abs(HX[k])):.2e}")
        return 0

    if a.vmc_npz:
        # plain VMC of the FIXED analytic vector c in [chi_I, xi, W1 psi, W2 psi, (Fcc psi, Fc1 psi)]
        import json
        an = np.load(a.vmc_npz, allow_pickle=True)
        c = np.asarray(an[a.vmc_key], float)
        gold = json.load(open(os.path.join(ROOT, "tests", "data", "lih_golden", "golden.json"), encoding="utf-8"))
        items = gold.get("items", gold)
        Fbar = float(items["pi_marriage_exp1.5_1e-3" if N > 20 else "pi_marriage_exp1.5_1e-2"]["values"]
                     ["pieces"]["value"]["Fbar"]["repr"])
        nx = len(c) - N - 1

        def psi_val(chi, X):
            return (c[:N] @ chi + c[N] * (X[1] - Fbar * X[0]) + c[N + 1:N + 1 + nx] @ X[2:2 + nx])

        def logpv(rr):
            chi, X, _h, _v, _g, _n = columns(rr, need_H=False)
            return np.log(np.abs(psi_val(chi, X)) ** 2 + 1e-300)
        est = []

        def meas(r):
            Phi, gPhi, lPhi = bas.phi_full(r)
            chi = Phi / 2.0; lchi = lPhi / 2.0
            Vp = V.potential(wf, r) - wf.Vnn
            Hchi = Vp[None] * chi - 0.5 * lchi
            _chi, X, HX, _VX, _GX, _VnX = columns(r)
            num = c[:N] @ Hchi + c[N] * (HX[1] - Fbar * HX[0]) + c[N + 1:N + 1 + nx] @ HX[2:2 + nx]
            est.append(float((num / psi_val(chi, X)).real.mean()) + wf.Vnn)
        t0 = time.time()
        accr, sd = S2.metropolis(bas, logpv, a.nwalk, a.nsweep, a.nburn, a.seed, a.thin, meas)
        est = np.array(est); nb = 20; m = len(est) // nb
        bl = est[:m * nb].reshape(nb, m).mean(1)
        E = bl.mean(); err = bl.std(ddof=1) / np.sqrt(nb)
        tag = os.path.splitext(os.path.basename(a.vmc_npz))[0]
        np.savez(os.path.join(DATA, f"en_factors_vmc_{tag}_{a.vmc_key}_s{a.seed}.npz"), est=est, E=E, err=err)
        print(f"[vmc] {a.vmc_key}: E = {E:.6f} +/- {err:.6f}  ({time.time()-t0:.0f}s acc={accr:.2f})", flush=True)
        return 0

    def logp(rr):
        Phi = bas.phi_vals(rr) / 2.0
        psi = cvec @ Phi
        return np.log(np.abs(psi) ** 2 + eps * np.sum(np.abs(Phi) ** 2, axis=0) + 1e-300)

    acc = dict(tS=np.zeros((N, K)), tH=np.zeros((N, K)), S=np.zeros((K, K)), H=np.zeros((K, K)),
               tV=np.zeros((N, K)), V=np.zeros((K, K)), Tg=np.zeros((K, K)), Vn=np.zeros((K, K)),
               c0=np.zeros(N), PP=0.0)
    blocks = {k: [] for k in acc}
    cnt = [0]

    def measure(r):
        chi, X, HX, VX, GX, VnX = columns(r)
        P = np.abs(cvec @ chi) ** 2 + eps * np.sum(np.abs(chi) ** 2, axis=0)
        q = Z / P
        nw = r.shape[0]
        cc = chi.conj() * q[None]
        Xc = X.conj() * q[None]
        acc["tS"] += (cc @ X.T).real / nw
        acc["tH"] += (cc @ HX.T).real / nw
        acc["S"] += (Xc @ X.T).real / nw
        acc["H"] += (Xc @ HX.T).real / nw
        acc["tV"] += (cc @ VX.T).real / nw
        acc["V"] += (Xc @ VX.T).real / nw
        acc["Vn"] += (Xc @ VnX.T).real / nw
        acc["Tg"] += 0.5 * np.einsum('kwed,lwed,w->kl', GX.conj(), GX, q).real / nw     # gradient form
        acc["c0"] += (cc @ X[0]).real / nw
        acc["PP"] += np.mean(q * np.abs(X[0]) ** 2)
        cnt[0] += 1
        if cnt[0] == a.blk:
            for k in acc:
                blocks[k].append(acc[k] / cnt[0])
                acc[k] = acc[k] * 0
            cnt[0] = 0

    t0 = time.time()
    print(f"[g2] trunc={a.trunc} N={N} seed={a.seed} core NO {sw.no_index} phi0={sw.phi0:+.6f}", flush=True)
    accr, sd = S2.metropolis(bas, logp, a.nwalk, a.nsweep, a.nburn, a.seed, a.thin, measure)
    out = {k: np.array(v) for k, v in blocks.items()}
    np.savez(os.path.join(DATA, f"en_factors_g2_{a.trunc:.0e}_s{a.seed}.npz"), cvec=cvec,
             cols=np.array(COLS), **out)
    print(f"[g2] done {time.time()-t0:.0f}s blocks={len(blocks['c0'])} acc={accr:.2f}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
