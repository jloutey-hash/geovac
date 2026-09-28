"""Stage C: the heteronuclear (Z_A=3, Z_B=1) separation constant on the Route C LiH orbitals.

Engines UNCHANGED (imported): debug/prolate_allelectron_c4.py (orbitals, one_body_general,
build_eri_tensor_m) and the FCI Slater-Condon machinery of debug/prolate_allelectron_fci.py.
Set: lih_orbitals_c4(R, extended=True, n_pi=1) at R = 3.015 (the logged C4 point,
E_tot = -8.00794, 10 orbitals, 2025 dets).

New here: the one-electron Lambda matrix over that orbital set.  Each orbital is a product
f(xi) g(eta) e^{i m phi}, so with the engine's own moment helpers

    Lambda_ij = Ns pref_S [ r2 * Gy + Fx * a2 ],
    Fx = Kxi + mu^2 X(P, mu-1) - R Z_+ X(P+1, mu),  Gy = Kang + mu^2 E(0, mu-1) + R Z_- E(1, mu),
    r2 = X(P+2, mu), a2 = E(2, mu)

and the same pieces re-assemble h1 = (2/R^2) pref_S Ns [Fx a0 + r0 Gy] (gate: must equal the
engine's h1).  The LiH set is NOT a tensor-product basis, so [h1, Lambda] no longer vanishes
exactly (Atkinson needs the product structure); its size is measured.

Diagnostics (as stage B): MOs = eigenvectors of the bare one-electron (Li^3+ + H^+) h1 in the
canonical-orthogonalized basis, each tagged with <Lambda>; FCI vector in that MO basis;
determinant weights by Sum Lambda; Lambda- vs energy- vs random- vs natural-orbital truncation.
Usage: python debug/prolate_separation_constant_stageC.py [build|analyse]
"""
from __future__ import annotations

import os
import sys
import time
from itertools import combinations
from typing import List

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import prolate_allelectron_c4 as c4  # noqa: E402
from prolate_allelectron_fci import _dets, _matel  # noqa: E402
from geovac import neumann_vee_general_m as ngm  # noqa: E402
from geovac import prolate_recondition as pr  # noqa: E402

R = 3.015
ZA, ZB = 3.0, 1.0
NPZ = os.path.join(HERE, "data", "prolate_separation_constant_stageC_ints.npz")


def lambda_matrix(orbs, R_, Z_A, Z_B):
    _pm, _shift, _pder = ngm._polymul, ngm._shift, c4._pder
    with mp.workdps(60):
        R_ = mp.mpf(R_)
        pref_S = 2 * mp.pi * (R_ / 2) ** 3
        Zp, Zm = mp.mpf(Z_A + Z_B), mp.mpf(Z_A - Z_B)
        n = len(orbs)
        L = np.zeros((n, n))
        Hc = np.zeros((n, n))
        for i in range(n):
            oi = orbs[i]
            for j in range(i, n):
                oj = orbs[j]
                if oi.mu != oj.mu or oi.msign != oj.msign:
                    continue
                mu = oi.mu
                gi, gj = list(oi.eta_poly), list(oj.eta_poly)
                A = ngm._mono_moments(oi.alpha + oj.alpha, oi.xi_power + oj.xi_power + 2 * mu + 8)
                P = oi.xi_power + oj.xi_power
                E = _pm(gi, gj)

                def X(nn, w):
                    return pr._mom_xi(_shift(pr._xi2m1(w), nn), A)

                def Em(dQ, w):
                    return pr._mom_eta(_pm(_shift(E, dQ), pr._meta2(w)))

                if mu == 0:
                    nxi_i, nxi_j = pr._mx_poly(oi.xi_power, oi.alpha), pr._mx_poly(oj.xi_power, oj.alpha)
                    ny_i, ny_j, wk = _pder(gi), _pder(gj), 1
                else:
                    nxi_i = pr._nx_poly(oi.xi_power, mu, oi.alpha)
                    nxi_j = pr._nx_poly(oj.xi_power, mu, oj.alpha)
                    ny_i, ny_j, wk = c4._ny_poly(gi, mu), c4._ny_poly(gj, mu), mu - 1
                Kxi = pr._mom_xi(_pm(_pm(nxi_i, nxi_j), pr._xi2m1(wk)), A)
                Kang = pr._mom_eta(_pm(_pm(ny_i, ny_j), pr._meta2(wk)))
                az_x = mu * mu * X(P, mu - 1) if mu > 0 else 0
                az_y = mu * mu * Em(0, mu - 1) if mu > 0 else 0
                Fx = Kxi + az_x - R_ * Zp * X(P + 1, mu)
                Gy = Kang + az_y + R_ * Zm * Em(1, mu)
                Ns = oi.norm * oj.norm * pref_S
                lam = Ns * (X(P + 2, mu) * Gy + Fx * Em(2, mu))
                h = Ns * (2 / R_ ** 2) * (Fx * Em(0, mu) + X(P, mu) * Gy)
                L[i, j] = L[j, i] = float(lam)
                Hc[i, j] = Hc[j, i] = float(h)
    return L, Hc


def build() -> None:
    t0 = time.time()
    orbs = c4.lih_orbitals_c4(R, extended=True, n_pi=1)
    S, h1 = c4.one_body_general(orbs, R, ZA, ZB)
    eri = c4.build_eri_tensor_m(orbs, R, verbose=True)
    Lm, Hc = lambda_matrix(orbs, R, ZA, ZB)
    np.savez(NPZ, S=S, h1=h1, eri=eri, L=Lm, Hc=Hc,
             mu=np.array([o.mu for o in orbs]), msign=np.array([o.msign for o in orbs]))
    print(f"built in {time.time() - t0:.0f}s -> {NPZ}", flush=True)


def fci_vec(h1, eri, M, nelec=4):
    na = nb = nelec // 2
    dets = _dets(M, na, nb)
    nd = len(dets)
    H = np.zeros((nd, nd))
    for a in range(nd):
        for b in range(a, nd):
            H[a, b] = H[b, a] = _matel(h1, eri, dets[a], dets[b])
    w, V = np.linalg.eigh(H)
    return w[0], V[:, 0], dets, H


def analyse() -> None:
    d = np.load(NPZ)
    S, h1, eri, Lm, Hc = d["S"], d["h1"], d["eri"], d["L"], d["Hc"]
    msign = d["msign"]
    gate = float(np.abs(Hc - h1).max() / np.abs(h1).max())
    print(f"gate: re-assembled h1 vs engine h1: rel max dev {gate:.1e}")
    w, U = np.linalg.eigh(S)
    X = U / np.sqrt(w)                  # all kept (cond ~6.5e4, as C4 with tol 1e-10)
    ho, Lo = X.T @ h1 @ X, X.T @ Lm @ X
    comm = float(np.linalg.norm(ho @ Lo - Lo @ ho, 2) / (np.linalg.norm(ho, 2) * np.linalg.norm(Lo, 2)))
    print(f"one-body rel ||[h1, Lambda]|| on the LiH set (non-product basis) = {comm:.2e}")
    eps, C = np.linalg.eigh(ho)
    lamM = C.T @ Lo @ C
    lam = np.diag(lamM)
    leak = np.sqrt(np.maximum(np.sum(lamM ** 2, 0) - lam ** 2, 0))
    # m of each MO
    ms = [int(round(np.sum((X @ C[:, k]) ** 2 * msign) / np.sum((X @ C[:, k]) ** 2))) for k in range(len(eps))]
    for k in range(len(eps)):
        print(f"  MO {k}: eps={eps[k]:+.5f}  <Lambda>={lam[k]:+9.3f}  leak||(L-<L>)phi||={leak[k]:.3f}  m={ms[k]:+d}")
    M = len(eps)
    h1m = C.T @ ho @ C
    Xc = X @ C
    erim = np.einsum("ap,bq,cr,ds,abcd->pqrs", Xc, Xc, Xc, Xc, eri, optimize=True)
    E0, v, dets, H = fci_vec(h1m, erim, M)
    Vnn = ZA * ZB / R
    print(f"FCI E_tot = {E0 + Vnn:.5f} (C4 log -8.00794), ndet={len(dets)}")
    wts = v ** 2
    # per determinant: sum of MO Lambda / eps over occupied spin-orbitals (so index = 2*spatial+spin)
    lam_d = np.array([sum(lam[i // 2] for i in D) for D in dets])
    eps_d = np.array([sum(eps[i // 2] for i in D) for D in dets])
    mean = float(wts @ lam_d)
    var = float(wts @ (lam_d - mean) ** 2)
    o = np.argsort(-wts)
    print(f"<Sum Lambda> = {mean:.3f}  sd = {np.sqrt(var):.3f}")
    for k in o[:8]:
        print(f"  det {dets[k]} w={wts[k]:.5f} SumLambda={lam_d[k]:.2f} SumEps={eps_d[k]:.4f}")
    nd = len(dets)
    rng = np.random.default_rng(0)
    i_hf = int(np.argmin(eps_d))
    for frac in (0.5, 0.3, 0.2, 0.1, 0.05):
        nk = max(1, int(round(frac * nd)))
        row = []
        for key in (lam_d, eps_d):
            sel = np.argsort(key, kind="stable")[:nk]
            row.append((np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0] - E0) * 1000)
        rr = []
        for _ in range(5):
            sel = np.concatenate([[i_hf], rng.choice(np.delete(np.arange(nd), i_hf), nk - 1, replace=False)])
            rr.append((np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0] - E0) * 1000)
        print(f"  keep {frac:4.2f} ({nk:4d}/{nd}): dE SumLambda {row[0]:8.3f} mHa | SumEps {row[1]:8.3f} | random {np.mean(rr):8.3f}")
    # orbital-level: keep n spatial orbitals ordered by Lambda / eps / NO occupation
    gam = np.zeros((M, M))
    idx = {D: k for k, D in enumerate(dets)}
    for k, D in enumerate(dets):       # diagonal part of the 1-RDM is enough to rank NOs? no: full
        pass
    # full 1-RDM via single excitations
    for a, Da in enumerate(dets):
        if abs(v[a]) < 1e-12:
            continue
        for so in Da:
            for t in range(2 * M):
                if t % 2 != so % 2 or (t != so and t in Da):
                    continue
                Db = tuple(sorted([x for x in Da if x != so] + [t]))
                b = idx.get(Db)
                if b is None:
                    continue
                lst = [x for x in Da if x != so]
                sgn = (-1) ** (Da.index(so) + Db.index(t))
                gam[t // 2, so // 2] += sgn * v[a] * v[b]
    gam = 0.5 * (gam + gam.T)
    nocc, Uno = np.linalg.eigh(gam)
    print(f"  1-RDM trace = {np.trace(gam):.6f} (4 expected); NO occupations {np.round(nocc[::-1], 5)}")
    ono = np.argsort(-nocc)
    Uno = Uno[:, ono]
    h_no = Uno.T @ h1m @ Uno
    e_no = np.einsum("ap,bq,cr,ds,abcd->pqrs", Uno, Uno, Uno, Uno, erim, optimize=True)
    for n_orb in (4, 5, 6, 7, 8):
        out = []
        for order, hh, ee in ((np.argsort(lam, kind="stable"), h1m, erim),
                              (np.arange(M), h1m, erim),
                              (np.arange(M), h_no, e_no)):
            sel = np.sort(order[:n_orb])
            Es, _ = __import__("prolate_allelectron_fci").fci_energy(hh[np.ix_(sel, sel)], ee[np.ix_(sel, sel, sel, sel)], n_orb, 4)
            out.append((Es - E0) * 1000)
        print(f"  n_orb={n_orb}: dE Lambda-ordered {out[0]:8.3f} mHa | eps-ordered {out[1]:8.3f} | NO {out[2]:8.3f}")


if __name__ == "__main__":
    (build if (len(sys.argv) > 1 and sys.argv[1] == "build") else analyse)()
