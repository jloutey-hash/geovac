"""Stage B: is the one-electron separation constant a useful (approximate) label for H2?

Uses the production Paper 12 re-conditioned prolate CI (geovac.prolate_recondition, direct
engine, Laguerre x Legendre) UNCHANGED for S, H1, V_ee.  Adds, per azimuthal block mu, the
one-electron operators in the SAME orthogonal basis (built from the module's own exact mpf
1D blocks):

    W1 = r2 (x) a0 - r0 (x) a2
    Hm1 = (2/R^2) [ Fx (x) a0 + r0 (x) Gy ],  Fx = Kr + mu^2 r0' - 2R r1,  Gy = Ka + mu^2 a0'
    N1 = r2 (x) Gy + Fx (x) a2                    (Lambda for Z_A = Z_B = 1)

(Z_- = 0: the eta^1 block is not needed.)  The joint (Hm1, N1) eigenvectors are the H2+
molecular orbitals in this basis, each carrying an exact Lambda.  The CI ground state is
re-expressed on MO-pair configurations (orthonormal), giving:
  * ||[H, Lambda_1]|| (all of it from V_ee: the one-body part commutes to machine precision);
  * <Lambda_1 + Lambda_2>, Var, the weight distribution over Lambda-configurations;
  * Lambda-truncated CI vs energy-truncated CI at equal configuration count.

Usage: python debug/prolate_separation_constant_stageB.py j_max l_max mu_max alpha
"""
from __future__ import annotations

import json
import os
import sys
import time
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from geovac import neumann_vee_general_m as ngm  # noqa: E402
from geovac import prolate_recondition as pr  # noqa: E402

R = pr.R_DEFAULT
BASIS = "laguerre_legendre"


def one_e_ops(j_max: int, l_max: int, mu: int, alpha: float, dps: int = 40
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with mp.workdps(dps):
        n_mom = 6 * max(j_max, l_max) + 6 * (mu + 2) + 20
        A = ngm._mono_moments(2.0 * alpha, n_mom)
        r0, r1, r2, Kr, r0p, r2p = pr._radial_blocks(j_max, mu, alpha, A, BASIS)
        a0, a2, Ka, a0p, a2p = pr._angular_blocks(l_max, mu, BASIS)
        f = pr._to_f64
        r0, r1, r2, Kr, a0, a2, Ka = map(f, (r0, r1, r2, Kr, a0, a2, Ka))
        if mu > 0:
            Kr = Kr + mu * mu * f(r0p)
            Ka = Ka + mu * mu * f(a0p)
    Fx = Kr - 2.0 * R * r1            # Z_+ = 2
    Gy = Ka                           # Z_- = 0
    W1 = np.kron(r2, a0) - np.kron(r0, a2)
    Hm1 = (2.0 / R ** 2) * (np.kron(Fx, a0) + np.kron(r0, Gy))
    N1 = np.kron(r2, Gy) + np.kron(Fx, a2)
    return W1, Hm1, N1


def mo_basis(W1: np.ndarray, Hm1: np.ndarray, N1: np.ndarray, l_max: int
             ) -> Dict[str, np.ndarray]:
    """Joint eigenbasis (W1-orthonormal).  Diagonalize H; read Lambda; report commutator."""
    d = np.sqrt(np.diag(W1))
    Wn, Hn, Nn = (M / np.outer(d, d) for M in (W1, Hm1, N1))
    E, V = eigh(0.5 * (Hn + Hn.T), 0.5 * (Wn + Wn.T))
    lamM = V.T @ Nn @ V
    lam = np.diag(lamM).copy()
    off = float(np.abs(lamM - np.diag(lam)).max() / np.abs(lam).max())
    # commutator in the W-orthonormal frame
    w, U = np.linalg.eigh(Wn)
    X = U / np.sqrt(w)
    Hh, Nh = X.T @ Hn @ X, X.T @ Nn @ X
    comm = float(np.linalg.norm(Hh @ Nh - Nh @ Hh, 2) /
                 (np.linalg.norm(Hh, 2) * np.linalg.norm(Nh, 2)))
    L1 = l_max + 1
    ang = np.arange(W1.shape[0]) % L1
    # eta-parity of each MO = parity of its Legendre support (Z_- = 0 keeps them pure)
    par = []
    imp = []
    for k in range(V.shape[1]):
        v = V[:, k]
        pe = np.sum(v[ang % 2 == 0] ** 2); po = np.sum(v[ang % 2 == 1] ** 2)
        par.append(0 if pe > po else 1)
        imp.append(min(pe, po) / (pe + po))
    return dict(imp=float(max(imp)), E=E, lam=lam, V=V / d[:, None], off=off, comm=comm, par=np.array(par),
                cond=float(w.max() / w.min()))


def main() -> None:
    j_max, l_max, mu_max = (int(a) for a in sys.argv[1:4])
    alpha = float(sys.argv[4])
    tag = f"{j_max}{l_max}{mu_max}_a{alpha}"
    t0 = time.time()
    idx = pr._product_index(j_max, l_max, mu_max)
    fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    Nr = (j_max + 1) ** 2
    Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
    Nmu = mu_max + 1
    l_neu = 2 * l_max + 4 * mu_max + 10
    with mp.workdps(pr.DEFAULT_DPS):
        Tr, Ta = pr._transforms_per_mu(BASIS, j_max, l_max, mu_max, alpha)
        S_o, H1_o = pr.build_one_body_direct(j_max, l_max, mu_max, alpha, R, BASIS)
        V = pr.vee_mp(fns, alpha, R, l_neu, False)
        V_o = pr._to_f64(pr._factored_cob(V, Nmu, Nr, Na, Tr, Ta))
    H_o = H1_o + V_o + (1.0 / R) * S_o
    print(f"built N={len(idx)} in {time.time() - t0:.0f}s", flush=True)
    E_ref = pr._normalized_solve(S_o, H_o)[0]

    L1 = l_max + 1
    idx_arr = np.array(idx)
    h6 = (R / 2) ** 6
    # MO-pair transform T: product-basis index -> MO-pair configuration
    cfg: List[Tuple[int, int, int]] = []   # (mu, p, q)
    cols = []
    mos = {}
    for mu in range(Nmu):
        W1, Hm1, N1 = one_e_ops(j_max, l_max, mu, alpha)
        mo = mo_basis(W1, Hm1, N1, l_max)
        mos[mu] = mo
        cc = 4 * np.pi ** 2 if mu == 0 else 2 * np.pi ** 2
        # consistency: product-basis S block must equal h6*cc*W1 (x) W1 on this mu
        rows = np.where(idx_arr[:, 4] == mu)[0]
        se1 = idx_arr[rows, 0] * L1 + idx_arr[rows, 1]
        se2 = idx_arr[rows, 2] * L1 + idx_arr[rows, 3]
        Sblk = h6 * cc * W1[np.ix_(se1, se1)] * W1[np.ix_(se2, se2)]
        dS = float(np.abs(Sblk - S_o[np.ix_(rows, rows)]).max() / np.abs(S_o[np.ix_(rows, rows)]).max())
        Vm = mo["V"]
        nse = Vm.shape[0]
        for p in range(nse):
            for q in range(nse):
                if (mo["par"][p] + mo["par"][q]) % 2:
                    continue              # gerade pair space only
                col = np.zeros(len(idx))
                col[rows] = Vm[se1, p] * Vm[se2, q] / np.sqrt(h6 * cc)
                cols.append(col)
                cfg.append((mu, p, q))
        print(f"  mu={mu}: MO comm={mo['comm']:.1e} offdiag={mo['off']:.1e} "
              f"cond(W1)={mo['cond']:.1e} parity-impurity={mo['imp']:.0e} S-block check {dS:.1e}; "
              f"MO eps[:4]={np.round(mo['E'][:4], 5)} lam[:4]={np.round(mo['lam'][:4], 4)}",
              flush=True)
    T = np.array(cols).T
    S_mo = T.T @ S_o @ T
    H_mo = T.T @ H_o @ T
    ortho_dev = float(np.abs(S_mo - np.eye(len(cfg))).max())
    H_mo = 0.5 * (H_mo + H_mo.T)
    Eall, Call = np.linalg.eigh(H_mo)
    E0, c = float(Eall[0]), Call[:, 0]
    print(f"S_MO orthonormality dev {ortho_dev:.1e}; E0(MO basis)={E0:.8f} vs "
          f"E_ref(recondition solve)={E_ref:.8f}", flush=True)

    lam_cfg = np.array([mos[mu]["lam"][p] + mos[mu]["lam"][q] for (mu, p, q) in cfg])
    eps_cfg = np.array([mos[mu]["E"][p] + mos[mu]["E"][q] for (mu, p, q) in cfg])
    lam1 = np.array([mos[mu]["lam"][p] for (mu, p, q) in cfg])
    w = c ** 2
    mean = float(w @ lam_cfg)
    var = float(w @ (lam_cfg - mean) ** 2)
    # commutator [H, Lambda_1] in the orthonormal MO-pair basis
    C1 = H_mo * (lam1[None, :] - lam1[:, None])
    Hv = H_mo.copy()
    rel_comm = float(np.linalg.norm(C1, 2) / (np.linalg.norm(H_mo, 2) * np.abs(lam1).max()))
    gs_comm = float(np.linalg.norm(C1 @ c))
    # one-body part alone commutes: build H without V_ee
    H1_mo = T.T @ (H1_o + (1.0 / R) * S_o) @ T
    C1b = H1_mo * (lam1[None, :] - lam1[:, None])
    rel_comm_1b = float(np.linalg.norm(C1b, 2) / (np.linalg.norm(H1_mo, 2) * np.abs(lam1).max()))

    # weights by configuration (unordered pair, mu)
    agg: Dict[Tuple, float] = {}
    for (mu, p, q), wi in zip(cfg, w):
        key = (mu,) + tuple(sorted((p, q)))
        agg[key] = agg.get(key, 0.0) + wi
    top = sorted(agg.items(), key=lambda kv: -kv[1])[:10]
    print(f"<L1+L2> = {mean:.5f}  Var = {var:.5f}  sd = {np.sqrt(var):.5f}")
    print(f"rel ||[H,L1]|| = {rel_comm:.3e}   one-body only = {rel_comm_1b:.1e}   "
          f"||[H,L1] psi|| = {gs_comm:.4e}")
    tops = []
    for (key, wt) in top:
        mu, p, q = key
        lp, lq = mos[mu]["lam"][p], mos[mu]["lam"][q]
        ep, eq = mos[mu]["E"][p], mos[mu]["E"][q]
        tops.append(dict(mu=mu, p=p, q=q, w=wt, lam_p=lp, lam_q=lq, eps_p=ep, eps_q=eq))
        print(f"  cfg mu={mu} ({p},{q}) w={wt:.6f}  Lambda=({lp:.4f},{lq:.4f})  eps=({ep:.4f},{eq:.4f})")

    # Lambda-sum distribution in bins
    order = np.argsort(lam_cfg)
    # truncation: keep configurations with smallest key, compare Lambda vs energy key
    ncfg = len(cfg)
    trunc = []
    for frac in (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05):
        nk = max(1, int(round(frac * ncfg)))
        out = dict(frac_kept=frac, n_kept=nk)
        for name, keyv in (("Lambda", lam_cfg), ("energy", eps_cfg),
                           ("Lambda_max", np.array([max(mos[mu]["lam"][p], mos[mu]["lam"][q])
                                                     for (mu, p, q) in cfg]))):
            sel = np.argsort(keyv, kind="stable")[:nk]
            e = float(np.linalg.eigvalsh(H_mo[np.ix_(sel, sel)])[0])
            out[name] = (e - E0) * 1000.0
        rng = np.random.default_rng(0)
        rr = []
        for _ in range(10):   # control: random configuration subsets always containing the HF-like (0,0)
            sel = np.concatenate([[int(np.argmin(eps_cfg))],
                                  rng.choice(np.delete(np.arange(ncfg), int(np.argmin(eps_cfg))),
                                             nk - 1, replace=False)]) if nk > 1 else np.array([int(np.argmin(eps_cfg))])
            rr.append((float(np.linalg.eigvalsh(H_mo[np.ix_(sel, sel)])[0]) - E0) * 1000.0)
        out["random_mean"] = float(np.mean(rr))
        trunc.append(out)
        print(f"  keep {frac:4.2f} ({nk:5d}/{ncfg}): dE Lambda-sum {out['Lambda']:8.3f} mHa | "
              f"energy-sum {out['energy']:8.3f} | Lambda-max {out['Lambda_max']:8.3f} | "
              f"random {out['random_mean']:8.3f}")
    # ---- orbital-level truncation: Lambda-ordered MOs vs energy-ordered MOs vs natural orbitals
    cfg_arr = np.array(cfg)
    pos = {tuple(x): i for i, x in enumerate(cfg)}
    Tno = np.zeros((ncfg, ncfg))
    occ = {}
    for mu in range(Nmu):
        nse = mos[mu]["V"].shape[0]
        Cm = np.zeros((nse, nse))
        for (m_, p, q), ci in zip(cfg, c):
            if m_ == mu:
                Cm[p, q] = ci
        gam = Cm @ Cm.T + Cm.T @ Cm
        nocc, U = np.linalg.eigh(gam)
        o = np.argsort(-nocc)
        nocc, U = nocc[o], U[:, o]
        occ[mu] = nocc
        rows = [i for i, x in enumerate(cfg) if x[0] == mu]
        for i in rows:
            _, p, q = cfg[i]
            for j in rows:
                _, a, b = cfg[j]
                Tno[i, j] = U[p, a] * U[q, b]
    H_no = Tno.T @ H_mo @ Tno
    orbs = [(mu, p) for mu in range(Nmu) for p in range(mos[mu]["V"].shape[0])]
    keys = {"Lambda": [mos[mu]["lam"][p] for (mu, p) in orbs],
            "energy": [mos[mu]["E"][p] for (mu, p) in orbs],
            "NO": [-occ[mu][p] for (mu, p) in orbs]}
    orb_trunc = []
    for n_orb in (2, 4, 6, 8, 10, 15, 20, 30, 45):
        row = dict(n_orb=n_orb)
        for name, kv in keys.items():
            keep = {orbs[i] for i in np.argsort(kv, kind="stable")[:n_orb]}
            sel = [i for i, (mu, p, q) in enumerate(cfg) if (mu, p) in keep and (mu, q) in keep]
            Hsel = (H_no if name == "NO" else H_mo)[np.ix_(sel, sel)]
            row[name] = ((float(np.linalg.eigvalsh(Hsel)[0]) - E0) * 1000.0) if sel else None
            row["ncfg_" + name] = len(sel)
        orb_trunc.append(row)
        print(f"  n_orb={n_orb:3d}: dE Lambda {row['Lambda']:9.3f} mHa ({row['ncfg_Lambda']} cfg) | "
              f"energy {row['energy']:9.3f} ({row['ncfg_energy']}) | NO {row['NO']:9.3f} ({row['ncfg_NO']})")
    res = dict(orb_trunc=orb_trunc, basis=[j_max, l_max, mu_max], alpha=alpha, N=len(idx), n_cfg=ncfg,
               E0=E0, E_ref=E_ref, ortho_dev=ortho_dev, mean=mean, var=var,
               rel_comm=rel_comm, rel_comm_one_body=rel_comm_1b, gs_comm=gs_comm,
               mo={str(mu): dict(E=mos[mu]["E"][:12].tolist(), lam=mos[mu]["lam"][:12].tolist(),
                                 par=mos[mu]["par"][:12].tolist(), comm=mos[mu]["comm"])
                   for mu in mos},
               top=tops, trunc=trunc, seconds=time.time() - t0)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", f"prolate_separation_constant_stageB_{tag}.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
