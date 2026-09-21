"""LiH R12-CI STAGE 3, part 1 -- sigma^2 = <Phi0|(F-Fbar)^2|Phi0> ground truth (VMC).

sigma^2 is the metric element S_11 of the 2x2 (S=[[1,0],[0,sigma2]]); it is exactly the
VARIANCE of F = sum_{i<j} f(r_ij) under the determinant density |Phi0|^2 -- a POSITIVE
quantity, so a Monte-Carlo estimate has NO cancellation (unlike forming <F^2>-Fbar^2 from
separately-noisy pieces, the Be ill-conditioning trap).  This is the well-conditioned ground
truth the analytic prolate TWOPAIR port (next increment) must reproduce, and it settles the
geminal choice empirically.

|Phi0|^2 for the closed-shell |m0^2 m1^2| factorizes into two INDEPENDENT 2-electron block
densities  D_p(r_a,r_b) = (m0(r_a)m1(r_b) - m1(r_a)m0(r_b))^2  (up block {1,2}, down block
{3,4}).  So one 6-D Metropolis sampler of D_p, run as two independent chains, gives the full
4-electron configuration.  m0,m1 = Loewdin(1s_A, 1s_B), the Stage-2 ionic reference.

Cross-checks (independent, whole-determinant, of the earlier stages):
  <F>_MC       vs analytic Fbar   = 1.880743  (Stage 1)
  <V_ee>_MC    vs analytic E2     = 3.614436  (Stage 2 two-electron energy)
Deliverable:
  sigma^2 = Var[F], with the conditioning ratio sigma^2/Fbar^2, for the current geminal
  f=exp(-0.5 r) AND the Be-recommended short-range f=r e^{-g r} (same samples; |Phi0|^2 is
  geminal-independent) -> the geminal recommendation for the analytic build.

Run from root:  python debug/lih_r12ci_sigma2_mc.py
"""
import numpy as np

from lih_r12ci_energy import ZA, ZB, N_A, N_B, S_AB, Z_A, Z_B
from lih_r12_4body_integral import R, a, CENTER_A, CENTER_B, f_gem   # f_gem = exp(-0.5 r)

rng = np.random.default_rng(20260921)

FBAR_REF = 1.880743      # Stage 1 analytic  <Phi0|F|Phi0>            (f = exp(-0.5 r))
VEE_REF = 3.614436       # Stage 2 analytic  E2 = <Phi0|V_ee|Phi0>
GAM_SR = 1.0             # short-range test geminal  f_sr = r e^{-GAM_SR r}

# --- Loewdin AO->MO (2x2), MOs m0,m1 = symmetric-orthogonalized 1s_A, 1s_B --------------- #
_sv, _U = np.linalg.eigh(np.array([[1.0, S_AB], [S_AB, 1.0]]))
X = _U @ np.diag(1.0 / np.sqrt(_sv)) @ _U.T          # AO_mu -> MO_p : m_p = sum_mu X[mu,p] AO_mu


def _oA(r):
    return N_A * np.exp(-ZA * np.linalg.norm(r - CENTER_A, axis=-1))


def _oB(r):
    return N_B * np.exp(-ZB * np.linalg.norm(r - CENTER_B, axis=-1))


def m0(r):
    return X[0, 0] * _oA(r) + X[1, 0] * _oB(r)


def m1(r):
    return X[0, 1] * _oA(r) + X[1, 1] * _oB(r)


def block_density(r1, r2):
    """D_p(r1,r2) = (m0(r1)m1(r2) - m1(r1)m0(r2))^2  (>=0), r1,r2 shape (...,3)."""
    return (m0(r1) * m1(r2) - m1(r1) * m0(r2)) ** 2


def sample_block(seed, nw=30000, burn=6000, nsnap=20, thin=60,
                 step_small=0.22, step_big=1.7, p_big=0.35):
    """Single-electron-update Metropolis on D_p with a MIXTURE proposal: mostly small steps
    (resolve the tight Li core, zeta=2.7) + occasional large steps (~R, hop the A<->B basins
    so the walker is ergodic across both nuclei -- the fix for the drift-form kinetic bias).
    Returns nsnap snapshots, each (nw,2,3), and the mean acceptance."""
    rg = np.random.default_rng(seed)
    r = np.empty((nw, 2, 3))
    # seed HALF the walkers with e1 near A / e2 near B and half swapped, to populate both basins
    half = nw // 2
    r[:half, 0] = CENTER_A + 0.4 * rg.standard_normal((half, 3))
    r[:half, 1] = CENTER_B + 0.9 * rg.standard_normal((half, 3))
    r[half:, 0] = CENTER_B + 0.9 * rg.standard_normal((nw - half, 3))
    r[half:, 1] = CENTER_A + 0.4 * rg.standard_normal((nw - half, 3))
    cur = block_density(r[:, 0], r[:, 1])
    snaps = []; accs = []
    total = burn + nsnap * thin
    for it in range(total):
        a = 0.0
        for k in range(2):                                  # single-electron updates
            big = rg.random(nw) < p_big
            sd = np.where(big, step_big, step_small)[:, None]
            prop = r[:, k] + sd * rg.standard_normal((nw, 3))
            new = block_density(prop, r[:, 1]) if k == 0 else block_density(r[:, 0], prop)
            m = rg.random(nw) < (new / np.maximum(cur, 1e-300))
            r[m, k] = prop[m]; cur[m] = new[m]
            a += m.mean()
        accs.append(0.5 * a)
        if it >= burn and (it - burn) % thin == 0:
            snaps.append(r.copy())
    return snaps, float(np.mean(accs[burn:]))


_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def _config_estimators(pos):
    """pos (M,4,3): up e1,e2 = 0,1 ; down e3,e4 = 2,3.  Returns F_exp, F_sr, V_ee, V_ne (M,)."""
    d = {ij: np.linalg.norm(pos[:, ij[0]] - pos[:, ij[1]], axis=-1) for ij in _PAIRS}
    F_exp = sum(f_gem(d[ij]) for ij in _PAIRS)                       # f = exp(-0.5 r)
    F_sr = sum(d[ij] * np.exp(-GAM_SR * d[ij]) for ij in _PAIRS)     # f = r e^{-GAM_SR r}
    V_ee = sum(1.0 / np.maximum(d[ij], 1e-12) for ij in _PAIRS)
    rA = np.linalg.norm(pos - CENTER_A, axis=-1)                     # (M,4)
    rB = np.linalg.norm(pos - CENTER_B, axis=-1)
    V_ne = (-Z_A / np.maximum(rA, 1e-12) - Z_B / np.maximum(rB, 1e-12)).sum(axis=1)
    return F_exp, F_sr, V_ee, V_ne


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI Stage 3 (part 1): sigma^2 = Var_{|Phi0|^2}[F]  (well-conditioned VMC)")
    print(f"  Phi0 = |m0^2 m1^2|, m=Loewdin(1s_A(Z={ZA}), 1s_B(Z={ZB})), <a|b>={S_AB:.6f}, R={R}")
    print("=" * 78)

    up, acc_u = sample_block(seed=1)
    dn, acc_d = sample_block(seed=2)
    nsnap = len(up)
    print(f"\n  Metropolis: {nsnap} snapshots x {up[0].shape[0]} walkers/block, "
          f"acceptance up/down = {acc_u:.2f}/{acc_d:.2f}")

    # per-snapshot estimator means (snapshots ~ independent -> batch-means error)
    Fe_s, Fs_s, Vee_s, Vne_s = [], [], [], []      # per-snapshot MEANS
    v2e_s, v2s_s = [], []                           # per-snapshot VARIANCES of F
    allFe, allFs = [], []
    for s in range(nsnap):
        pos = np.concatenate([up[s], dn[s]], axis=1)               # (nw,4,3)
        Fe, Fs, Vee, Vne = _config_estimators(pos)
        Fe_s.append(Fe.mean()); Fs_s.append(Fs.mean())
        Vee_s.append(Vee.mean()); Vne_s.append(Vne.mean())
        v2e_s.append(Fe.var(ddof=1)); v2s_s.append(Fs.var(ddof=1))
        allFe.append(Fe); allFs.append(Fs)
    Fe_s = np.array(Fe_s); Fs_s = np.array(Fs_s)
    Vee_s = np.array(Vee_s); Vne_s = np.array(Vne_s)
    v2e_s = np.array(v2e_s); v2s_s = np.array(v2s_s)
    sem = lambda x: x.std(ddof=1) / np.sqrt(len(x))

    Fbar = Fe_s.mean();  eF = sem(Fe_s)
    Vee = Vee_s.mean();  eVee = sem(Vee_s)
    Vne = Vne_s.mean();  eVne = sem(Vne_s)
    sig2 = np.concatenate(allFe).var(ddof=1); esig2 = sem(v2e_s)
    Fbar_sr = Fs_s.mean(); eFsr = sem(Fs_s)
    sig2_sr = np.concatenate(allFs).var(ddof=1); esig2sr = sem(v2s_s)

    print("\n--- independent whole-determinant cross-checks of Stages 1 & 2 ---")
    print(f"  <F>_MC   = {Fbar:.6f} +/- {eF:.1e}   vs Stage-1 analytic Fbar = {FBAR_REF:.6f}"
          f"   (dev {abs(Fbar-FBAR_REF):.1e} = {abs(Fbar-FBAR_REF)/eF:.1f} sigma)")
    print(f"  <V_ee>_MC= {Vee:.6f} +/- {eVee:.1e}   vs Stage-2 analytic E2   = {VEE_REF:.6f}"
          f"   (dev {abs(Vee-VEE_REF):.1e} = {abs(Vee-VEE_REF)/eVee:.1f} sigma)")
    print(f"  <V_ne>_MC= {Vne:.6f} +/- {eVne:.1e}   (info; V_ne part of E1)")

    print("\n--- DELIVERABLE: sigma^2 = Var[F] (metric element S_11) ---")
    print(f"  current geminal  f = exp(-0.5 r):  Fbar = {Fbar:.5f},  "
          f"sigma^2 = {sig2:.6f} +/- {esig2:.1e}")
    print(f"     conditioning ratio sigma^2 / Fbar^2 = {sig2/Fbar**2:.4f}   "
          f"(catastrophic if << 1e-3; healthy if O(1e-1))")
    print(f"  short-range      f = r e^(-{GAM_SR:.1f} r): Fbar_sr = {Fbar_sr:.5f}, "
          f"sigma^2_sr = {sig2_sr:.6f} +/- {esig2sr:.1e}")
    print(f"     conditioning ratio sigma^2_sr / Fbar_sr^2 = {sig2_sr/Fbar_sr**2:.4f}")

    print("\n--- verdict ---")
    ok_stage = (abs(Fbar - FBAR_REF) < 5 * eF) and (abs(Vee - VEE_REF) < 5 * eVee)
    print(f"  Stages 1&2 cross-validated by independent VMC: {ok_stage}")
    better = "short-range r e^{-g r}" if (sig2_sr / Fbar_sr ** 2) > (sig2 / Fbar ** 2) else "exp(-0.5 r)"
    print(f"  better-conditioned geminal (larger variance/mean^2 -> less cancellation): {better}")
    print("  NEXT: analytic prolate TWOPAIR port must reproduce this sigma^2; then h, then g.")
