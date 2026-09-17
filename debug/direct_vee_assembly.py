"""V_ee X+V assembly (sigma sector, mu=0), float64, in the orthogonal basis.

Reuses neumann_vee_general_m's validated Xtab (mpf, incl. the IBP corr) and the
eta moments Ytab, re-bases them per (l,m,s) to the PRODUCT-Laguerre / product-
Legendre index (small blocks, mpf), then assembles V_orth in vectorized float64.
Validates against prolate_recondition's V (mu=0) re-based (= ground truth).

mu=0 has only m=0 (phi_cec(0,0,m) != 0 only for m=0), so this is the m=0 Neumann
path -- the cleanest case to validate the X+V mechanism before the mu>0 extension.
"""
from __future__ import annotations
import time
import numpy as np
import mpmath as mp
from geovac import prolate_recondition as pr
from geovac import neumann_vee_general_m as ngm

R = pr.R_DEFAULT


def build_vorth_mu0(j_max, l_max, alpha, l_neumann=None):
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, 0)
        N = len(idx)
        c = mp.mpf(2.0 * alpha)
        s = 0                                            # mu_i=mu_j=0, m=0 -> s=0
        p_max = 2 * j_max + 2                            # room for the +2 jac shift
        q_max = 2 * l_max + 2
        if l_neumann is None:
            l_neumann = q_max                            # eta selection caps here
        l_hi = min(l_neumann, q_max)                     # l > q_max+2s-m=q_max -> zero
        # --- ngm X table (mpf), m=0,s=0 ---
        Xtab = pr._build_Xtab_mp([(0, 0)], l_hi, p_max, alpha,
                                 {(0, 0): l_hi})          # X[(l,0,0)][P1][P2]
        # --- eta moments Ytab[(l,Q)] (mpf) ---
        Ytab = {}
        for l in range(l_hi + 1):
            pmpoly = list(ngm._RP_poly(l, 0))
            yp = [mp.mpf(1)]                              # (1-eta^2)^0
            for Qq in range(q_max + 3):
                if l > Qq or (Qq + l) % 2 != 0:
                    Ytab[(l, Qq)] = mp.mpf(0)
                else:
                    Ytab[(l, Qq)] = pr._mom_eta(pr._pm(pr._shift(yp, Qq), pmpoly))
        # --- product-poly coeffs (mpf): radial L_a.L_a', angular P_b.P_b' ---
        n_r, n_a = j_max + 1, l_max + 1
        Lc = [pr.laguerre_coeffs(a, alpha, n_r) for a in range(n_r)]
        Pc = [pr.legendre_coeffs(b, n_a) for b in range(n_a)]
        PPrad = {}
        for a in range(n_r):
            for a2 in range(n_r):
                PPrad[(a, a2)] = ngm._polymul(Lc[a], Lc[a2])    # up to deg 2*j_max
        PPang = {}
        for b in range(n_a):
            for b2 in range(n_a):
                PPang[(b, b2)] = ngm._polymul(Pc[b], Pc[b2])

        # --- re-base X and Y to orthogonal PAIR index, per shift, downcast float64 ---
        radpairs = [(a, a2) for a in range(n_r) for a2 in range(n_r)]
        angpairs = [(b, b2) for b in range(n_a) for b2 in range(n_a)]
        # X_orth[l][(dP1,dP2)] : dict (radpair, radpair) -> float
        Xo = {}
        for l in range(l_hi + 1):
            X = Xtab.get((l, 0, 0))
            if X is None:
                continue
            for dP1 in (0, 2):
                for dP2 in (0, 2):
                    M = np.zeros((len(radpairs), len(radpairs)))
                    for i1, (a, a2) in enumerate(radpairs):
                        w1 = PPrad[(a, a2)]
                        for i2, (cc, cc2) in enumerate(radpairs):
                            w2 = PPrad[(cc, cc2)]
                            tot = mp.mpf(0)
                            for P1 in range(len(w1)):
                                if w1[P1] == 0 or P1 + dP1 > p_max:
                                    continue
                                for P2 in range(len(w2)):
                                    if w2[P2] == 0 or P2 + dP2 > p_max:
                                        continue
                                    tot += w1[P1] * w2[P2] * X[P1 + dP1][P2 + dP2]
                            M[i1, i2] = float(tot)
                    Xo[(l, dP1, dP2)] = M
        # Y_orth[l][dQ] : angpair -> float
        Yo = {}
        for l in range(l_hi + 1):
            for dQ in (0, 2):
                v = np.zeros(len(angpairs))
                for ip, (b, b2) in enumerate(angpairs):
                    w = PPang[(b, b2)]
                    tot = mp.mpf(0)
                    for Q in range(len(w)):
                        if w[Q] != 0:
                            tot += w[Q] * Ytab.get((l, Q + dQ), mp.mpf(0))
                    v[ip] = float(tot)
                Yo[(l, dQ)] = v

    # --- vectorized float64 assembly ---
    h6 = float((mp.mpf(R) / 2) ** 6)
    pref = float((2 / mp.mpf(R)) * h6)
    fphi = float(pr._phi_cec(0, 0, 0))                   # m=0
    idx_arr = np.array(idx)
    rp_index = {p: i for i, p in enumerate(radpairs)}
    ap_index = {p: i for i, p in enumerate(angpairs)}
    # per basis function: radial pair elec1 (j), elec2 (k); angular pair elec1 (l), elec2 (m)
    ja, la, ka, ma = idx_arr[:, 0], idx_arr[:, 1], idx_arr[:, 2], idx_arr[:, 3]
    jac = [(+1, 2, 0, 2, 0), (-1, 2, 0, 0, 2), (-1, 0, 2, 2, 0), (+1, 0, 2, 0, 2)]
    V = np.zeros((N, N))
    for l in range(l_hi + 1):
        if (l, 0, 0) not in Xtab:
            continue
        npre = float(pr._neumann_prefactor(l, 0))
        coef = pref * fphi * npre
        for sgn, dP1, dQ1, dP2, dQ2 in jac:
            Xm = Xo[(l, dP1, dP2)]
            Y1 = Yo[(l, dQ1)]
            Y2 = Yo[(l, dQ2)]
            # radial-pair indices for elec1 (a_i,a_j) and elec2 (c_i,c_j) for all i,j
            rp1 = (ja[:, None] * (j_max + 1) + ja[None, :])   # (a_i, a_j) flat
            rp2 = (ka[:, None] * (j_max + 1) + ka[None, :])
            ap1 = (la[:, None] * (l_max + 1) + la[None, :])
            ap2 = (ma[:, None] * (l_max + 1) + ma[None, :])
            V += sgn * coef * Xm[rp1, rp2] * Y1[ap1] * Y2[ap2]
    return V, idx


def ground_truth_vorth_mu0(j_max, l_max, alpha, l_neumann=None):
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, 0)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        if l_neumann is None:
            l_neumann = 2 * l_max + 10
        V = pr.vee_mp(fns, alpha, R, l_neumann)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu("laguerre_legendre", j_max, l_max, 0, alpha)
        V_o = pr._factored_cob(V, 1, Nr, Na, Tr, Ta)
        N = len(fns)
        return np.array([[float(V_o[i, j]) for j in range(N)] for i in range(N)])


if __name__ == "__main__":
    alpha = 1.0
    for (j_max, l_max) in [(2, 2), (3, 3)]:
        t0 = time.time()
        Vd, idx = build_vorth_mu0(j_max, l_max, alpha)
        td = time.time() - t0
        t1 = time.time()
        Vg = ground_truth_vorth_mu0(j_max, l_max, alpha)
        tg = time.time() - t1
        rel = np.max(np.abs(Vd - Vg) / np.maximum(np.abs(Vg), 1e-10))
        print(f"({j_max},{l_max}) mu=0 N={Vd.shape[0]:3d}  relV={rel:.2e}  "
              f"[direct {td:.2f}s vs mpf {tg:.1f}s]", flush=True)
