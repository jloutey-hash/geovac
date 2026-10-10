"""Paper 60 [OPEN] item: does a basis-size-bounded preconditioner exist for the two-centre
Shibuya-Wulfman (SW) metric on the full sigma basis (n <= N, l < n, m = 0 per centre)?

Diagnostic driver (2026-10-09).  Imports the tracked builder geovac/sturmian_sigma_basis.py.
Output: debug/data/p60_sigma_preconditioner.json.  Memo: debug/p60_sigma_preconditioner_memo.md.

Sections
  A  spectra of the cross block C and of C.Pi (Pi = diag((-1)^l)); involution test (C Pi)^2 -> I
  B  soft-subspace characterisation (l / n / centre weights)
  C  preconditioner levers: centre-block-diagonal (lemma), l-block Jacobi, n-shell Jacobi,
     per-l tri(1,2,1) in n, s-sector tri on l=0 only, fixed-rank ideal deflation
  D  canonical orthogonalisation: H2+ R=2 Galerkin energy vs threshold (accuracy cost)
  E  quadrature floor (160 vs 240 points)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geovac.sturmian_sigma_basis as sb  # noqa: E402

H2P_R2 = -1.102634214
KRS = (1.0, 3.0, 8.0)
NS = (2, 3, 4, 5, 6, 7)


def tridiag(n):
    return 2 * np.eye(n) + np.eye(n, k=1) + np.eye(n, k=-1)


def inv_sqrt(M):
    ev, U = np.linalg.eigh(0.5 * (M + M.T))
    return U @ np.diag(ev ** -0.5) @ U.T


def pcond(A, P):
    X = inv_sqrt(P)
    ev = np.linalg.eigvalsh(X @ A @ X)
    return float(ev[-1] / ev[0])


def cond(A):
    ev = np.linalg.eigvalsh(0.5 * (A + A.T))
    return float(ev[-1] / ev[0])


def blocks(S, labels):
    m = len(labels)
    C = S[:m, m:]
    Pi = np.diag([(-1.0) ** l for (_, l) in labels])
    return C, Pi


def l_tri_precond(labels):
    """block-Toeplitz tri(1,2,1) in n, one block per l (labels ordered by n then l)."""
    m = len(labels)
    P = np.zeros((m, m))
    for l in sorted({l for _, l in labels}):
        idx = [i for i, (n, ll) in enumerate(labels) if ll == l]
        P[np.ix_(idx, idx)] = tridiag(len(idx))
    return P


def l0_tri_precond(labels):
    m = len(labels)
    P = np.eye(m)
    idx = [i for i, (n, ll) in enumerate(labels) if ll == 0]
    P[np.ix_(idx, idx)] = tridiag(len(idx))
    return P


def block_jacobi(A, groups):
    P = np.zeros_like(A)
    for g in groups:
        P[np.ix_(g, g)] = A[np.ix_(g, g)]
    return P


def main():
    t0 = time.time()
    out = {"meta": {"date": "2026-10-09", "builder": "geovac/sturmian_sigma_basis.py",
                    "quadrature": f"{sb.NX}x{sb.NE}"}}

    # ---------------- E: quadrature floor
    floor = []
    for N in (5, 7):
        for kR in KRS:
            S1, _, lab = sb.build(N, 1.0, kR)
            m = len(lab)
            intra = float(np.max(np.abs(S1[:m, :m] - np.eye(m))))
            floor.append({"N": N, "kR": kR, "max_intra_dev": intra})
    # 240-point comparison at N=7 (module globals swapped then restored)
    saved = (sb.NX, sb.NE, sb._TX, sb._WX, sb._EE, sb._WE)
    S160 = {kR: sb.build(7, 1.0, kR)[0] for kR in KRS}
    sb.NX, sb.NE = 176, 320   # Gauss-Laguerre overflows exp(t) beyond ~178 nodes
    sb._TX, sb._WX = np.polynomial.laguerre.laggauss(176)
    sb._EE, sb._WE = np.polynomial.legendre.leggauss(320)
    for kR in KRS:
        S240 = sb.build(7, 1.0, kR)[0]
        e1 = np.linalg.eigvalsh(S160[kR])
        e2 = np.linalg.eigvalsh(S240)
        floor.append({"N": 7, "kR": kR, "max_entry_diff_vs_176x320": float(np.abs(S160[kR] - S240).max()),
                      "lam_min_160": float(e1[0]), "lam_min_176x320": float(e2[0]),
                      "intra_dev_176x320": float(np.max(np.abs(S240[:28, :28] - np.eye(28))))})
    sb.NX, sb.NE, sb._TX, sb._WX, sb._EE, sb._WE = saved
    out["quadrature_floor"] = floor

    # ---------------- A/C: spectra + levers
    spec, lever, soft = [], [], []
    for kR in KRS:
        for N in NS:
            S, _, lab = sb.build(N, 1.0, kR)
            m = len(lab)
            C, Pi = blocks(S, lab)
            sig = np.linalg.svd(C, compute_uv=False)
            CP = C @ Pi
            asym = float(np.abs(CP - CP.T).max())
            muCP = np.linalg.eigvalsh(0.5 * (CP + CP.T))
            evS = np.linalg.eigvalsh(S)
            Sg, Su = sb.gerade_ungerade(S, lab)
            evg, evu = np.linalg.eigvalsh(Sg), np.linalg.eigvalsh(Su)
            # involution test on the N=2 (first 3) labels: ((C Pi)^2)_low -> I
            CP2 = CP @ CP
            inv_dev = float(np.abs(CP2[:3, :3] - np.eye(3)).max())
            spec.append({
                "kR": kR, "N": N, "fns_per_centre": m, "cond_S": float(evS[-1] / evS[0]),
                "lam_min_S": float(evS[0]), "cond_g": float(evg[-1] / evg[0]),
                "cond_u": float(evu[-1] / evu[0]), "lam_min_g": float(evg[0]),
                "lam_min_u": float(evu[0]),
                "n_S_below_1e-2": int((evS < 1e-2).sum()), "n_S_above_2m1e-2": int((evS > 2 - 1e-2).sum()),
                "n_sigma_gt_0.99": int((sig > 0.99).sum()), "n_sigma_gt_0.9": int((sig > 0.9).sum()),
                "frac_sigma_gt_0.9": float((sig > 0.9).mean()),
                "frac_|mu(CPi)|_gt_0.9": float((np.abs(muCP) > 0.9).mean()),
                "n_mu_lt_-0.99": int((muCP < -0.99).sum()), "n_mu_gt_0.99": int((muCP > 0.99).sum()),
                "CPi_asym": asym, "involution_dev_low3": inv_dev,
            })

            # ---- levers
            row = {"kR": kR, "N": N, "fns_per_centre": m,
                   "raw_S": float(evS[-1] / evS[0]), "raw_g": float(evg[-1] / evg[0]),
                   "raw_u": float(evu[-1] / evu[0])}
            # (1) centre-block-diagonal preconditioners on full S (lemma: no gain possible)
            PA = l_tri_precond(lab)
            row["S_blockdiag_ltri"] = pcond(S, np.block([[PA, 0 * PA], [0 * PA, PA]]))
            rng = np.random.default_rng(60 + N)
            Ga = rng.normal(size=(m, m)); Gb = rng.normal(size=(m, m))
            Pr = np.block([[Ga @ Ga.T + m * np.eye(m), np.zeros((m, m))],
                           [np.zeros((m, m)), Gb @ Gb.T + m * np.eye(m)]])
            row["S_blockdiag_random"] = pcond(S, Pr)
            # (2) on gerade / ungerade blocks
            lgroups = [[i for i, (n, ll) in enumerate(lab) if ll == l] for l in range(N)]
            ngroups = [[i for i, (n, ll) in enumerate(lab) if n == nn] for nn in range(1, N + 1)]
            for name, A in (("g", Sg), ("u", Su)):
                row[f"{name}_lJacobi"] = pcond(A, block_jacobi(A, lgroups))
                row[f"{name}_nJacobi"] = pcond(A, block_jacobi(A, ngroups))
                row[f"{name}_ltri"] = pcond(A, l_tri_precond(lab))
                row[f"{name}_l0tri"] = pcond(A, l0_tri_precond(lab))
                row[f"{name}_cond_worst_lblock"] = max(cond(A[np.ix_(g, g)]) for g in lgroups)
                ev = np.linalg.eigvalsh(A)
                for r in (1, 3, 5):
                    row[f"{name}_deflate_r{r}"] = float(ev[-1] / ev[r]) if r < len(ev) else None
                row[f"{name}_rank_to_cond100"] = int((ev < ev[-1] / 100).sum())
                row[f"{name}_rank_to_cond1e4"] = int((ev < ev[-1] / 1e4).sum())
            lever.append(row)

            # ---- B: soft-vector characterisation (gerade block, eigenvalue < 1e-2)
            if kR == 3.0 and N in (4, 5, 6, 7):
                ev, V = np.linalg.eigh(Sg)
                sel = ev < 1e-2
                W = (V[:, sel] ** 2)
                lw = {l: float(W[[i for i, (n, ll) in enumerate(lab) if ll == l]].sum() / max(sel.sum(), 1))
                      for l in range(N)}
                nw = {n: float(W[[i for i, (nn, ll) in enumerate(lab) if nn == n]].sum() / max(sel.sum(), 1))
                      for n in range(1, N + 1)}
                lfrac = {l: sum(1 for (_, ll) in lab if ll == l) / m for l in range(N)}
                purity = [float(max(W[[i for i, (n, ll) in enumerate(lab) if ll == l], j].sum()
                                    for l in range(N))) for j in range(W.shape[1])]
                # stiff (near-2) side for the symmetric Slepian picture
                soft.append({"N": N, "n_soft": int(sel.sum()), "n_stiff_gt_1.98": int((ev > 1.98).sum()),
                             "l_weight_soft": lw, "l_weight_uniform": lfrac, "n_weight_soft": nw,
                             "max_l_purity_per_vector": purity,
                             "soft_eigs": [float(x) for x in ev[sel]]})
            print(f"kR={kR} N={N} m={m} condS={evS[-1]/evS[0]:.3e} t={time.time()-t0:.0f}s", flush=True)
    out["spectra"], out["levers"], out["soft_vectors"] = spec, lever, soft

    # ---------------- D: canonical orthogonalisation accuracy cost, H2+ R=2
    def top(W, S, tol):
        ev, V = np.linalg.eigh(S)
        keep = ev > tol * ev[-1]
        X = V[:, keep] / np.sqrt(ev[keep])
        return float(np.linalg.eigvalsh(X.T @ (0.5 * (W + W.T)) @ X)[-1]), int((~keep).sum())

    canon = []
    for N in (3, 4, 5, 6, 7):
        for tol in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
            try:
                k = brentq(lambda kk: top(*sb.build(N, kk, 2.0)[1::-1], tol)[0] - kk,
                           0.9, 2.5, xtol=1e-12)
            except ValueError:
                canon.append({"N": N, "tol": tol, "err": None})
                continue
            S, W, lab = sb.build(N, k, 2.0)
            _, ndrop = top(W, S, tol)
            ev = np.linalg.eigvalsh(S)
            canon.append({"N": N, "fns_total": 2 * len(lab), "tol": tol, "k": k,
                          "E": -0.5 * k * k, "err": -0.5 * k * k - H2P_R2, "n_dropped": ndrop,
                          "cond_kept": float(ev[-1] / ev[ndrop])})
            print(f"canon N={N} tol={tol:g} err={-0.5*k*k-H2P_R2:.3e} drop={ndrop}", flush=True)
    out["canonical_orth_h2plus_R2"] = canon
    out["meta"]["wall_s"] = time.time() - t0
    Path("debug/data").mkdir(exist_ok=True)
    with open("debug/data/p60_sigma_preconditioner.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
