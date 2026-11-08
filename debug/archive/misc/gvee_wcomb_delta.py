r"""Re-assemble an EXISTING sparse-marriage npz with the corrected triangle weight, without a rerun.

The 2026-09-26 g_Vee defect (debug/sprint_gvee_defect_audit_memo.md) is one weight in the general-m
triangle: combos (A>=1, B=0) got pi/2 instead of pi.  Only <F^2 V_ee>'s triangle class changes, so
    F2Vee_new = F2Vee_old + Delta,  Delta = SUM_triangle-triples [the (A>=1,B=0) combos at weight pi/2],
and g_Vee_new = g_Vee_old + Delta (Vee, FVee, Fbar, every other piece untouched).  This script replays
ONLY the triangle triples of <F^2 V_ee> (lih_gvee_sparse replay, same reference build) and evaluates
the triangle entries with the DELTA weight (new - old), skipping every combo whose delta is zero (so
only the non-symmetric (A>=1,B=0) randomized SVDs run).  It then re-solves the 2x2 (M.assemble_c).

    python debug/gvee_wcomb_delta.py --npz debug/data/lih_exp_geminal_exp_g1.5_..._1e-02.npz
    (geminal, gamma, trunc and grid are read from the npz)
Prints old/new g_Vee, E_R12, c; writes <npz>_wcombfix.npz next to it (the original is not touched).
"""
from __future__ import annotations

import os
import sys
import time
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def delta_w(A: int, B: int) -> float:
    new = 2 * np.pi * (1.0 if A == 0 else 0.5) * (1.0 if B == 0 else 0.5)
    old = 2 * np.pi if (A == 0 and B == 0) else 0.5 * (2 * np.pi if A == 0 else np.pi)
    return new - old


def tri_entries_delta(M, S, ctx, cm, il, ir):
    """lih_gvee_sparse.tri_entries with the delta weight; zero-delta combos skipped."""
    refm, gem = ctx.refm, ctx.gem
    mu = refm.rho[cm]; mu_mid = int(refm.mA[cm])
    rho = refm.rho; mArows = refm.mA
    gmu = S.P2.GEO * mu
    out = np.zeros(len(il))
    mL = mArows[il]; mR = mArows[ir]
    for A in range(M.MMAX + 1):
        for B in range(M.MMAX + 1):
            dw = delta_w(A, B)
            if dw == 0.0:
                continue
            for sa in ((0,) if A == 0 else (+1, -1)):
                for sb in ((0,) if B == 0 else (+1, -1)):
                    if sa * A + sb * B != mu_mid:
                        continue
                    sig, U, W = M._rand_svd_middle(gem.Fm[A], gem.Fm[B], gmu, rank=30)
                    if sig.size == 0:
                        continue
                    for muL in sorted(set(int(x) for x in mL)):
                        c = abs(muL + sa * A)
                        psel = (mL == muL) & (mR == -(muL + mu_mid))
                        if not psel.any():
                            continue
                        uL, linv = np.unique(il[psel], return_inverse=True)
                        uR, rinv = np.unique(ir[psel], return_inverse=True)
                        nR = len(uR)
                        rhoR = rho[uR]; rhoL = rho[uL]
                        blk = np.zeros((len(uL), nR))
                        rk = sig.size
                        rstep = max(1, ctx.chunk_rows // max(nR, 1))
                        for r0 in range(0, rk, rstep):
                            r1 = min(r0 + rstep, rk)
                            Bst = (rhoR[None, :, :] * W[:, r0:r1].T[:, None, :]).reshape((r1 - r0) * nR, S.NG)
                            Pr = M.coul_mode_stack_any(Bst, int(c)).reshape(r1 - r0, nR, S.NG)
                            for r in range(r0, r1):
                                Lr = rhoL * (U[:, r] * S._GW2)[None, :]
                                blk += dw * sig[r] * (Lr @ Pr[r - r0].T)
                        out[np.nonzero(psel)[0]] += blk[linv, rinv]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npz', required=True)
    ap.add_argument('--mem-cap-gb', type=float, default=3.0)
    a = ap.parse_args()
    z = np.load(a.npz, allow_pickle=True)
    gam = float(z['gam']); geminal = str(z['geminal']); trunc = float(z['trunc'])
    grid = (int(z['NXI']), int(z['NETA']), int(z['NPHI']))
    zg = lambda k, d: int(z[k]) if k in z.files else d
    assert zg('tri_mmax', 4) == 4 and zg('rank', 30) == 30, "delta tool assumes default MMAX/rank"
    import lih_exp_geminal_run as RUN
    KG = RUN.install_kernels(gam, None if grid == (72, 44, 28) else grid)
    import lih_marriage_phase2 as P2
    import lih_pi_m16 as M
    import lih_gvee_sparse as S
    if zg('nphi_gen', 28) != 28:
        P2.build_kernel_gen.__defaults__ = (zg('nphi_gen', 28),)
    if zg('nphi_mode', 48) != 48:
        P2.build_kernel_m_gen.__defaults__ = (zg('nphi_mode', 48),)
        RUN_TR = __import__('geovac.lih_r12ci.triangle', fromlist=['x'])
        RUN_TR.build_kernel_m.__defaults__ = (zg('nphi_mode', 48),)
    KG.USE_EXACT_NEUMANN = True
    t0 = time.time()
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=trunc, verbose=False)
    assert abs(info['E_trunc'] - float(z['E_trunc'])) < 1e-9
    gem = M.GeminalM(geminal)
    ctx = S.Ctx(gem, refm, a.mem_cap_gb)
    col = S.TriCollector()
    ntri = 0
    for p in PAIRS6:
        for q in PAIRS6:
            for r in PAIRS6:
                rem_f, kk, gk = S._triple_args(p, q, r, gem)
                if kk != 'coul' or len(rem_f) != 2:
                    continue
                verts = set(rem_f[0]) | set(rem_f[1])
                mid = set(rem_f[0]) & set(rem_f[1])
                if len(mid) != 1 or (verts - mid) != set(r):
                    continue                                  # not a triangle triple
                n0 = len(col.recs)
                val = S.vee_reduce_sparse(rem_f, r, kk, gk, ctx, col)
                assert val == 0.0 and len(col.recs) > n0
                ntri += 1
    im = np.concatenate([x[0] for x in col.recs]); il = np.concatenate([x[1] for x in col.recs])
    ir = np.concatenate([x[2] for x in col.recs]); w = np.concatenate([x[3] for x in col.recs])
    n = ctx.nact
    code = (im * n + il) * n + ir
    uc, inv = np.unique(code, return_inverse=True)
    wsum = np.bincount(inv, weights=w)
    ucm = uc // (n * n); uil = (uc // n) % n; uir = uc % n
    delta = 0.0
    for cm in col.order:
        s = ucm == cm
        if int(refm.mA[cm]) == 0:
            continue                                          # mu_mid = 0 has no (A>=1,B=0) combo
        vals = tri_entries_delta(M, S, ctx, cm, uil[s], uir[s])
        delta += float((wsum[s] * vals).sum())
    p = {k: float(z[k]) for k in z.files if z[k].ndim == 0 and np.issubdtype(z[k].dtype, np.number)}
    E_old, c_old = float(z['E_R12']), float(z['c'])
    p_new = dict(p); p_new['F2Vee'] = p['F2Vee'] + delta; p_new['g_Vee'] = p['g_Vee'] + delta
    E_new, c_new, h_new = M.assemble_c(p_new, p_new['E0_grid'])
    E_chk, c_chk, _ = M.assemble_c(p, p['E0_grid'])
    print(f"{geminal} gamma={gam} trunc={trunc:.0e}: {ntri} triangle triples, Delta F2Vee = {delta:+.8f}  "
          f"[{time.time()-t0:.0f}s]")
    print(f"  re-assembly check (old pieces): E {E_chk:.9f} vs banked {E_old:.9f}")
    print(f"  g_Vee  {p['g_Vee']:+.8f} -> {p_new['g_Vee']:+.8f}")
    print(f"  E_R12  {E_old:.9f} -> {E_new:.9f}   (dE {1e3*(E_new-E_old):+.4f} mHa)")
    print(f"  c      {c_old:+.6f} -> {c_new:+.6f}")
    out = dict(z); out.update(F2Vee=p_new['F2Vee'], g_Vee=p_new['g_Vee'], E_R12=E_new, c=c_new, h=h_new,
                              wcomb_delta=delta, E_R12_old=E_old, c_old=c_old, g_Vee_old=p['g_Vee'])
    np.savez(a.npz[:-4] + "_wcombfix.npz", **out)


if __name__ == "__main__":
    main()
