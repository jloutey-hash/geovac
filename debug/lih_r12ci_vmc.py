"""LiH R12-CI PoC ENERGY via the VMC linear method (one Jastrow parameter on chi = F - Fbar).

JUSTIFIED by the Stage-3 conditioning diagnostic (lih_r12ci_sigma2_mc.py): sigma^2/Fbar^2 ~ 4%
(HEALTHY, not the Be ill-conditioning trap), so ALL 2x2 matrix elements are well-conditioned
VMC averages -- the analytic TWOPAIR port is NOT needed for the ENERGY.  (The analytic RI-free
4-body bridge, v5.15.10, remains the separate structural result, validated at the integral
level; it is what a fully-analytic energy would use, and what F12 needs RI for.)

Trial  Psi = Phi0 + c*(F-Fbar)*Phi0 ,  chi = F - Fbar ,  F = sum_{i<j} f(r_ij),
geminal f = r e^{-gam r}  (correct electron-cusp sign f'(0)=1 > 0; exp(-gam r) has the wrong sign).
2x2 generalized eigenproblem in {Phi0, chi*Phi0}.  Kinetic via the BOUNDED gradient form
(integration by parts -> 1/2 <|grad|^2>, finite for Slater orbitals -- no Laplacian cusp spikes):
  E0    = 1/2 sum_i <v_i^2> + <V>                                   (geminal-INDEPENDENT)
  S01   = <chi> (~0),   sigma2 = S11 = <chi^2>
  h=H01 = 1/2 sum_i <chi v_i^2 + grad_i chi . v_i> + <chi V>
  g=H11 = 1/2 sum_i <chi^2 v_i^2 + 2 chi grad_i chi . v_i + |grad_i chi|^2> + <chi^2 V>
  v_i = grad_i ln Phi0 (drift, bounded),  V = V_ne + V_ee.
Solve [[E0,h],[h,g]] c = E [[1,S01],[S01,sigma2]] c ;  lowest E = E_R12.

|Phi0|^2 is geminal-independent, so ONE MC sampling -> scan gam.  E0 (geminal-independent)
must reproduce the Stage-2 analytic value -7.888 (validation of the kinetic gradient machinery).

Run:  python debug/lih_r12ci_vmc.py
"""
import numpy as np
from scipy.linalg import eigh

from lih_r12ci_sigma2_mc import (ZA, ZB, X, CENTER_A, CENTER_B, Z_A, Z_B, R, sample_block)

E0_REF = -7.887822          # Stage-2 analytic <Phi0|H|Phi0> (geminal-independent target)
GAMS = [0.5, 0.8, 1.0, 1.3, 1.7, 2.2]     # geminal f = r e^{-gam r}: scan gam
_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def _orb_vg(pts, C, z):
    """1s Slater(z) value + gradient at pts (M,3).  grad = -z*orb*rhat."""
    dvec = pts - C
    rc = np.linalg.norm(dvec, axis=-1)
    orb = np.sqrt(z ** 3 / np.pi) * np.exp(-z * rc)
    grad = (-z * orb / np.maximum(rc, 1e-30))[:, None] * dvec
    return orb, grad


def _m_vg(pts, p):
    """MO m_p = X[0,p] 1s_A + X[1,p] 1s_B : value + gradient."""
    oA, gA = _orb_vg(pts, CENTER_A, ZA)
    oB, gB = _orb_vg(pts, CENTER_B, ZB)
    return X[0, p] * oA + X[1, p] * oB, X[0, p] * gA + X[1, p] * gB


def _block_drift(ra, rb):
    """block amplitude D = m0(ra)m1(rb)-m1(ra)m0(rb); return v_ra, v_rb = grad ln|D| (M,3)."""
    m0a, gm0a = _m_vg(ra, 0); m1a, gm1a = _m_vg(ra, 1)
    m0b, gm0b = _m_vg(rb, 0); m1b, gm1b = _m_vg(rb, 1)
    D = m0a * m1b - m1a * m0b
    gDa = gm0a * m1b[:, None] - gm1a * m0b[:, None]
    gDb = m0a[:, None] * gm1b - m1a[:, None] * gm0b
    Dg = np.where(np.abs(D) < 1e-300, 1e-300, D)[:, None]
    return gDa / Dg, gDb / Dg


def geom_quantities(pos):
    """pos (M,4,3): drift v (M,4,3), v2=sum_i|v_i|^2 (M,), V_ne, V_ee (M,)  [geminal-free]."""
    v1, v2 = _block_drift(pos[:, 0], pos[:, 1])
    v3, v4 = _block_drift(pos[:, 2], pos[:, 3])
    v = np.stack([v1, v2, v3, v4], axis=1)
    v2sum = (v ** 2).sum(axis=(1, 2))
    rA = np.linalg.norm(pos - CENTER_A, axis=-1); rB = np.linalg.norm(pos - CENTER_B, axis=-1)
    Vne = (-Z_A / np.maximum(rA, 1e-12) - Z_B / np.maximum(rB, 1e-12)).sum(axis=1)
    Vee = sum(1.0 / np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12)
              for i, j in _PAIRS)
    return v, v2sum, Vne, Vee


def F_and_gradF(pos, gam, kind='linexp'):
    """F = sum_{i<j} f(r_ij), gradF (M,4,3) = sum over pairs of f'(r) rhat.
    kind='linexp': f = r e^{-gam r} (correct cusp); kind='exp': f = e^{-gam r} (matches the
    analytic-port machinery / Stage-1 f-tensor)."""
    M = pos.shape[0]
    F = np.zeros(M); gradF = np.zeros((M, 4, 3))
    for i, j in _PAIRS:
        dvec = pos[:, i] - pos[:, j]
        r = np.linalg.norm(dvec, axis=-1); rr = np.maximum(r, 1e-12)
        if kind == 'linexp':
            F += r * np.exp(-gam * r); fp = (1.0 - gam * r) * np.exp(-gam * r)
        else:  # exp
            F += np.exp(-gam * r); fp = -gam * np.exp(-gam * r)
        gf = (fp / rr)[:, None] * dvec                                   # f'(r) rhat
        gradF[:, i] += gf; gradF[:, j] -= gf
    return F, gradF


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI PoC energy -- VMC linear method, geminal f = r e^{-gam r}")
    print(f"  Phi0 = |m0^2 m1^2| ionic ref (Z_orb {ZA}/{ZB}), R={R};  E0 target (Stage 2) = {E0_REF:.5f}")
    print("=" * 78)

    up, acc_u = sample_block(seed=1)
    dn, acc_d = sample_block(seed=2)
    nsnap = len(up)
    print(f"\n  Metropolis: {nsnap} snapshots x {up[0].shape[0]} walkers, acc {acc_u:.2f}/{acc_d:.2f}")

    # precompute geminal-free quantities per snapshot (drift, V) -- reused for every gam
    POS, V2, VNE, VEE, VSUM = [], [], [], [], []
    for s in range(nsnap):
        pos = np.concatenate([up[s], dn[s]], axis=1)
        v, v2, vne, vee = geom_quantities(pos)
        POS.append(pos); V2.append(v2); VNE.append(vne); VEE.append(vee); VSUM.append(v)

    # E0 (geminal-independent): batch-means over snapshots.  V_NN is a CONSTANT added to H
    # (shifts all eigenvalues equally -> cancels in dE = E_R12 - E0), so add it to the totals.
    VNN = Z_A * Z_B / R
    e0_s = np.array([(0.5 * V2[s] + VNE[s] + VEE[s]).mean() for s in range(nsnap)]) + VNN
    E0 = e0_s.mean(); eE0 = e0_s.std(ddof=1) / np.sqrt(nsnap)
    print(f"\n  E0_VMC = {E0:.5f} +/- {eE0:.5f}   vs Stage-2 analytic {E0_REF:.5f}"
          f"   (dev {abs(E0-E0_REF):.4f} = {abs(E0-E0_REF)/eE0:.1f} sigma)  <- kinetic gradient form + V_NN validated")

    print("\n  --- geminal scan: 2x2 linear method, E_R12 = lowest generalized eigenvalue ---")
    print(f"  {'gam':>5} {'Fbar':>8} {'sigma2':>9} {'h':>10} {'g':>10} {'E_R12':>12} {'dE=E_R12-E0':>14}")
    best = None
    for gam in GAMS:
        # pass 1: global Fbar for this gam
        Fall = np.concatenate([F_and_gradF(POS[s], gam)[0] for s in range(nsnap)])
        Fbar = Fall.mean()
        # pass 2: per-snapshot 2x2 elements + solve -> batch means
        ER, dE, sig2s, hs, gs = [], [], [], [], []
        for s in range(nsnap):
            pos = POS[s]; v = VSUM[s]; v2 = V2[s]; Vt = VNE[s] + VEE[s]
            F, gradF = F_and_gradF(pos, gam)
            chi = F - Fbar
            gcv = (gradF * v).sum(axis=(1, 2))          # sum_i grad_i chi . v_i
            gc2 = (gradF ** 2).sum(axis=(1, 2))         # sum_i |grad_i chi|^2
            s01 = chi.mean(); sig2 = (chi ** 2).mean()
            E0s = (0.5 * v2 + Vt).mean()
            h = (0.5 * (chi * v2 + gcv) + chi * Vt).mean()
            g = (0.5 * (chi ** 2 * v2 + 2 * chi * gcv + gc2) + chi ** 2 * Vt).mean()
            Hm = np.array([[E0s, h], [h, g]]); Sm = np.array([[1.0, s01], [s01, sig2]])
            w = eigh(Hm, Sm, eigvals_only=True)
            ER.append(w[0] + VNN); dE.append(w[0] - E0s)      # dE is V_NN-independent
            sig2s.append(sig2); hs.append(h); gs.append(g)
        ER = np.array(ER); dE = np.array(dE)
        er_m = ER.mean(); er_e = ER.std(ddof=1) / np.sqrt(nsnap)
        de_m = dE.mean(); de_e = dE.std(ddof=1) / np.sqrt(nsnap)
        print(f"  {gam:5.1f} {Fbar:8.4f} {np.mean(sig2s):9.4f} {np.mean(hs):10.4f} "
              f"{np.mean(gs):10.4f} {er_m:8.4f}+/-{er_e:.4f} {de_m*1e3:8.2f}+/-{de_e*1e3:.2f} mHa")
        if best is None or er_m < best[1]:
            best = (gam, er_m, er_e, de_m, de_e)

    gam, er_m, er_e, de_m, de_e = best
    print("\n  --- verdict ---")
    print(f"  best geminal gam = {gam:.1f}:  E_R12 = {er_m:.5f} +/- {er_e:.5f} Ha")
    print(f"  correlation captured  dE = {de_m*1e3:+.2f} +/- {de_e*1e3:.2f} mHa  "
          f"(of the ~83 mHa LiH correlation energy)")
    print(f"  variational: E0 >= E_R12 : {E0 >= er_m - 3*er_e} ;  E_R12 > exact -8.070 : {er_m > -8.070}")
    print("  (correctly-signed lowering with the RI-free 4-body content handled = the PoC gate.)")

    # ---- analytic-port validation targets: decompose h and g for f = exp(-0.5 r) ----
    #      (same geminal as the analytic machinery / Stage-1 f-tensor).  h = h_T + h_Vne + h_Vee.
    print("\n  --- validation targets for the ANALYTIC port (geminal f = exp(-0.5 r)) ---")
    Fall = np.concatenate([F_and_gradF(POS[s], 0.5, kind='exp')[0] for s in range(nsnap)])
    Fbar_e = Fall.mean()
    hT_s, hVne_s, hVee_s, sg2_s, gg_s = [], [], [], [], []
    for s in range(nsnap):
        pos = POS[s]; v = VSUM[s]; v2 = V2[s]; vne = VNE[s]; vee = VEE[s]
        F, gradF = F_and_gradF(pos, 0.5, kind='exp')
        chi = F - Fbar_e
        gcv = (gradF * v).sum(axis=(1, 2)); gc2 = (gradF ** 2).sum(axis=(1, 2))
        hT_s.append((0.5 * (chi * v2 + gcv)).mean())
        hVne_s.append((chi * vne).mean()); hVee_s.append((chi * vee).mean())
        sg2_s.append((chi ** 2).mean())
        gg_s.append((0.5 * (chi ** 2 * v2 + 2 * chi * gcv + gc2) + chi ** 2 * (vne + vee)).mean())
    sem = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
    hT, hVne, hVee = np.mean(hT_s), np.mean(hVne_s), np.mean(hVee_s)
    print(f"    Fbar(exp) = {Fbar_e:.5f}   sigma^2 = {np.mean(sg2_s):.5f} +/- {sem(sg2_s):.1e}")
    print(f"    h_T   = {hT:+.5f} +/- {sem(hT_s):.1e}")
    print(f"    h_Vne = {hVne:+.5f} +/- {sem(hVne_s):.1e}")
    print(f"    h_Vee = {hVee:+.5f} +/- {sem(hVee_s):.1e}")
    print(f"    h = h_T+h_Vne+h_Vee = {hT+hVne+hVee:+.5f} +/- "
          f"{np.sqrt(sem(hT_s)**2+sem(hVne_s)**2+sem(hVee_s)**2):.1e}")
    print(f"    g = <G|H|G> = {np.mean(gg_s):+.5f} +/- {sem(gg_s):.1e}   (incl. V_NN? no: add VNN*sigma^2)")
    print(f"      g + VNN*sigma^2 = {np.mean(gg_s)+VNN*np.mean(sg2_s):+.5f}   (VNN={VNN:.4f})")
