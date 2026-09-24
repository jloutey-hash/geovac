r"""Phase 3 (final) of the LiH "marriage" build (2026-09-23; plan debug/lih_marriage_build_plan.md,
memo debug/sprint_lih_marriage_memo.md).  Two milestones:

M1 -- g_Vee productionization + the 74-determinant marriage energy E_R12.
    The Phase-2 machinery (debug/lih_marriage_phase2.py) is G1-validated on the 2-MO reference;
    its ONLY wall-clock blocker on the 74-det (55 active pair) reference was g_Vee's per-triple
    recomputation of the triangle mode-m Coulomb solves.  Phase 3 caches + batches those
    (get_tri_matrix / coul_mode_stack in phase2; NOT a math change).  This driver re-checks the
    optimization is EXACT on the 2-MO reference (<0.001 mHa vs the frozen G1 value) BEFORE running
    the 74-det case, then computes E_R12 = assemble(all_pieces(gem, ref_74det), E0_grid) on the
    exact-Neumann path, the correlation lowering vs E_trunc = -8.0229, and the 2x2 mixing c.

M2 -- G2 VMC cross-check.  Psi_T = (1 + c (F - Fbar)) Psi_CI is real-space evaluable
    (lih_vmc.LinearGeminalJastrow / vmc_linear, STANDARD local energy, analytic Laplacian).
    The reference Psi_CI is the SAME truncated 74-det natural-orbital vector (from the Phase-0b
    artifact).  Gate: E_VMC(c) == E_R12 within +-1 mHa, and E_VMC(c=0) == E_trunc (gate-6 sanity).

Run from root:
    python debug/lih_marriage_phase3.py --stage all > debug/data/lih_marriage_phase3.log 2>&1
    echo $? > debug/data/lih_marriage_phase3.exit
Options: --stage {m1,m2,all}  --gems linexp,exp  --vmc-nwalk N --vmc-nsweep N --vmc-nburn N
         --vmc-seeds s1,s2  --probe (time one 74-det g_Vee + one VMC batch and stop)
Writes debug/data/lih_marriage_phase3.npz  (pieces, E_R12, c, VMC results, verdict).

NOTE ON THREADS: this driver spawns NO Monte-Carlo worker processes (the VMC is a single-process
vectorised Metropolis), so it does NOT pin *_NUM_THREADS -- multi-threaded BLAS is safe and speeds
the g_Vee matmuls.  (phase2 pins them at its import; that is inherited but harmless here.)
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from itertools import combinations
from typing import Dict, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

REF_EXACT = os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz")
OUT = os.path.join(DATA, "lih_marriage_phase3.npz")

# frozen G1 (2-MO) exact-path references (debug/sprint_lih_marriage_memo.md Phase 2)
FROZEN_2MO = {"linexp": -7.917198, "exp": -7.943536}
EXACT_BELOW = -8.0705      # variational floor (frozen falsifier)
ORB_CEIL = -8.032          # orbital ceiling to beat


def _t(t0: float) -> str:
    return f"[{time.time() - t0:7.1f}s]"


# =========================================================================== #
# M2 support: reconstruct the truncated 74-det natural-orbital wavefunction as a
# lih_vmc.LiHWavefunction (restricted to the active alpha/beta pairs for speed).
# =========================================================================== #
def build_truncated_wf():
    """LiHWavefunction for the SAME truncated 74-det NO reference the analytic 2x2 uses.
    E_fci is set to E_trunc so the c=0 VMC gate-6 target is E_trunc = -8.0229."""
    import lih_vmc as V
    z = np.load(REF_EXACT, allow_pickle=True)
    prims = list(z['prims'])
    T_no = np.asarray(z['T_no'])                        # (M, Mk)
    CS_t = np.asarray(z['CS_t'])                        # (n_ap, n_bp) folded, sign-included
    apairs_all = [tuple(int(x) for x in pr) for pr in z['apairs']]   # combinations(range(Mk),2)
    act_a = np.asarray(z['act_a']).astype(int)
    act_b = np.asarray(z['act_b']).astype(int)
    R = float(z['R']); Vnn = float(z['Vnn']); Z_A = float(z['Z_A']); Z_B = float(z['Z_B'])
    E_trunc = float(z['E_trunc'])
    Mk = int(np.asarray(z['occ_no']).size)
    # restrict to the active pairs (all nonzero CS entries kept; psi identical, far fewer minors)
    apairs = [apairs_all[i] for i in act_a]
    bpairs = [apairs_all[j] for j in act_b]
    CS = CS_t[np.ix_(act_a, act_b)]
    na = nb = 2
    wf = V.LiHWavefunction(prims, T_no, None, apairs, bpairs, CS, R, Z_A, Z_B, Vnn,
                           E_trunc - Vnn, Mk, na, nb)
    return wf, E_trunc, dict(nact_a=len(act_a), nact_b=len(act_b), n_nz=int(np.count_nonzero(CS)),
                             Mk=Mk, M=len(prims))


# =========================================================================== #
# 2x2 mixing coefficient c (lower eigenvector's (F-Fbar)Phi0-to-Phi0 ratio)
# =========================================================================== #
def mixing_c(pieces: dict, E0_tot: float, V_NN: float) -> Tuple[float, float, dict]:
    sig2 = pieces['sigma2']
    h = pieces['h_T'] + pieces['h_Vne'] + pieces['h_Vee']
    g_elec = pieces['g_T'] + pieces['g_Vne'] + pieces['g_Vee']
    aa = E0_tot - V_NN                                  # electronic reference energy
    cc = g_elec / sig2
    bb = h / np.sqrt(sig2)
    lam_minus = 0.5 * (aa + cc - np.sqrt((aa - cc) ** 2 + 4 * bb ** 2))
    E_R12 = lam_minus + V_NN
    ratio = (lam_minus - aa) / bb                       # x1/x0 (normalised-basis eigenvector)
    c_mix = ratio / np.sqrt(sig2)                       # (F-Fbar)Phi0 : Phi0 coefficient ratio
    return E_R12, c_mix, dict(aa=aa, bb=bb, cc=cc, sig2=sig2, h=h, g_elec=g_elec,
                              lam_minus=lam_minus, ratio=ratio)


def _print_pieces(tag: str, p: dict, E_R12: float, c_mix: float, E0_tot: float, V_NN: float):
    print(f"\n  --- {tag} pieces ---")
    for k in ('Fbar', 'sigma2', 'h_T', 'h_Vne', 'h_Vee', 'g_T', 'g_Vne', 'g_Vee',
              'Vne_exp', 'Vee', 'T_kin', 'E0_grid'):
        print(f"    {k:10s} = {p[k]:+.6f}")
    print(f"    E0_tot(used) = {E0_tot:+.6f}   V_NN = {V_NN:+.6f}")
    print(f"    >>> E_R12    = {E_R12:.6f}   c(mix) = {c_mix:+.6f}")


# =========================================================================== #
def run_m1(gems, T0, probe=False):
    import geovac.lih_r12ci.kernels as KG
    import lih_marriage_phase2 as P2
    from geovac.lih_r12ci import gVee as GV
    from geovac.lih_r12ci.triangle import _MODE_CACHE

    V_NN = P2.V_NN
    print("=" * 100)
    print("M1: g_Vee PRODUCTIONIZATION (batched/cached triangle) + the 74-det marriage energy E_R12")
    print(f"  V_NN={V_NN:.6f}  R={P2.R_ENG}  GAM={P2.GAM}  frozen 2-MO (exact): {FROZEN_2MO}")
    print("=" * 100, flush=True)

    KG.USE_EXACT_NEUMANN = True
    GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); P2._TRI_MAT_CACHE.clear()

    # ---------------- M1a: exactness re-check on the 2-MO reference ----------------
    print("\n[M1a] EXACTNESS re-check of the batched/cached g_Vee on the 2-MO reference:")
    ref2 = P2.build_ref_2mo()
    recheck = {}
    for gname in gems:
        gem = P2.Geminal(gname)
        p = P2.all_pieces(gem, ref2, T0, verbose=False)
        E_R12_fixed, _ = P2.assemble(p, -7.887822)      # assembly convention E0 for direct compare
        dfroz = (E_R12_fixed - FROZEN_2MO[gname]) * 1e3
        print(f"    {gname:7s}: g_Vee={p['g_Vee']:+.6f}  E_R12(fixed-E0)={E_R12_fixed:.6f}  "
              f"vs frozen {FROZEN_2MO[gname]:.6f}  d={dfroz:+.4f} mHa  "
              f"{'PASS' if abs(dfroz) < 1e-3 else 'FAIL (>0.001 mHa)'}", flush=True)
        recheck[gname] = dict(g_Vee=float(p['g_Vee']), E_R12=float(E_R12_fixed), d_mHa=float(dfroz))
    exact_ok = all(abs(recheck[g]['d_mHa']) < 1e-3 for g in gems)
    print(f"  [M1a] exactness: {'PASS' if exact_ok else 'FAIL'} (all geminals <0.001 mHa vs frozen)")

    # ---------------- M1b: the 74-det energy ----------------
    print("\n[M1b] 74-det marriage energy (exact-Neumann path):")
    GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); P2._TRI_MAT_CACHE.clear()
    refA = P2.build_ref_artifact()
    m1 = {}
    for gi, gname in enumerate(gems):
        print(f"\n  === geminal = {gname} ===", flush=True)
        t_g = time.time()
        gem = P2.Geminal(gname)
        # clear only the flag-independent triangle/mode caches between geminals (densities differ)
        _MODE_CACHE.clear(); P2._TRI_MAT_CACHE.clear()
        p = P2.all_pieces(gem, refA, T0, verbose=True)
        E0_grid = p['E0_grid']
        E_R12, c_mix, info = mixing_c(p, E0_grid, V_NN)
        _print_pieces(f"74-det {gname}", p, E_R12, c_mix, E0_grid, V_NN)
        dcorr = (E_R12 - E0_grid) * 1e3
        above_floor = E_R12 > EXACT_BELOW
        beats_ceil = E_R12 < ORB_CEIL
        in_window = -8.055 <= E_R12 <= -8.045
        print(f"    E_trunc(E0_grid) = {E0_grid:.6f} ;  correlation lowering = {dcorr:+.4f} mHa")
        print(f"    variational bound E_R12 > {EXACT_BELOW}: {'OK' if above_floor else 'VIOLATED (BUG!)'}")
        print(f"    beats orbital ceiling {ORB_CEIL}: {beats_ceil} ;  in -8.045..-8.055 window: {in_window}")
        print(f"    geminal wall {_t(t_g)}", flush=True)
        if not above_floor:
            print(f"    *** STOP: E_R12 = {E_R12:.6f} < {EXACT_BELOW} -- variational bound violated; "
                  f"likely a bug.  Reporting and halting this geminal.")
        m1[gname] = dict(pieces={k: (float(v) if isinstance(v, (int, float, np.floating)) else None)
                                 for k, v in p.items()},
                         E_R12=float(E_R12), c_mix=float(c_mix), E0_grid=float(E0_grid),
                         Fbar=float(p['Fbar']), dcorr_mHa=float(dcorr),
                         above_floor=bool(above_floor), beats_ceil=bool(beats_ceil),
                         in_window=bool(in_window), info={k: float(v) for k, v in info.items()})
        if probe:
            print("  [--probe] stopping after the first geminal's g_Vee.")
            break
    return dict(recheck=recheck, exact_ok=bool(exact_ok), m1=m1, V_NN=float(V_NN))


# =========================================================================== #
def run_m2(m1: dict, gems, T0, nwalk, nsweep, nburn, seeds, probe=False):
    import lih_vmc as V
    print("\n" + "=" * 100)
    print("M2: G2 VMC CROSS-CHECK  Psi_T = (1 + c (F - Fbar)) Psi_CI  (standard local energy)")
    print("=" * 100, flush=True)

    wf, E_trunc, meta = build_truncated_wf()
    print(f"  truncated wf: M={meta['M']} Mk={meta['Mk']} active pairs a/b={meta['nact_a']}/{meta['nact_b']} "
          f"(CS nnz={meta['n_nz']})  E_trunc={E_trunc:.6f}  R={wf.R}", flush=True)

    # gate-6 sanity: VMC of Psi_CI (c=0) reproduces E_trunc.  Coarse by design (a +-4 sigma bar,
    # not the G2 +-1 mHa bar) -> 1 seed, reduced sweeps, so the wall-clock budget goes to the G2 runs.
    print("\n  [gate-6 sanity] VMC(Psi_CI, c=0) vs E_trunc:")
    lin0 = V.LinearGeminalJastrow('linexp', 0.5, 0.0, 0.0, wf.na, wf.nb)
    g6_nsweep = 400 if probe else min(nsweep, 2000)
    g6_nburn = 100 if probe else min(nburn, 500)
    E0s = []
    for s in seeds[:1]:
        E0, e0, a0 = V.vmc_linear(wf, lin0, nwalk=nwalk, nsweep=g6_nsweep,
                                  nburn=g6_nburn, seed=s, verbose=True)
        E0s.append((E0, e0))
    Eg6 = np.mean([x[0] for x in E0s]); eg6 = np.mean([x[1] for x in E0s])
    dev6 = (Eg6 - E_trunc) / max(eg6, 1e-9)
    print(f"  gate-6: E_VMC(c=0) = {Eg6:.5f} +/- {eg6:.5f}  vs E_trunc {E_trunc:.5f}  "
          f"dev = {dev6:+.1f} sigma  {'PASS' if abs(dev6) < 4 else 'FAIL'}", flush=True)

    if probe:
        print("  [--probe] stopping after the gate-6 batch.")
        return dict(gate6=dict(E=float(Eg6), err=float(eg6), dev=float(dev6), E_trunc=float(E_trunc)))

    # G2 at the analytic c for each geminal
    m2 = {}
    for gname in gems:
        if gname not in m1:
            continue
        c_mix = m1[gname]['c_mix']; Fbar = m1[gname]['Fbar']; E_R12 = m1[gname]['E_R12']
        print(f"\n  --- G2 for geminal={gname}: c={c_mix:+.6f} Fbar={Fbar:.6f} target E_R12={E_R12:.6f} ---",
              flush=True)
        lin = V.LinearGeminalJastrow(gname, 0.5, c_mix, Fbar, wf.na, wf.nb)
        Es, errs = [], []
        for s in seeds:
            E, err, acc = V.vmc_linear(wf, lin, nwalk=nwalk, nsweep=nsweep, nburn=nburn,
                                       seed=s, verbose=True)
            Es.append(E); errs.append(err)
        Es = np.array(Es); errs = np.array(errs)
        Emean = Es.mean()
        Eerr = (Es.std(ddof=1) / np.sqrt(len(Es))) if len(Es) > 1 else errs[0]
        diff = abs(Emean - E_R12) * 1e3
        combo_sig = np.sqrt(Eerr ** 2)
        g2_pass = diff <= 1.0
        print(f"  G2 [{gname}]: E_VMC = {Emean:.5f} +/- {Eerr:.5f}   E_R12 = {E_R12:.6f}   "
              f"|E_VMC - E_R12| = {diff:.3f} mHa   {'PASS' if g2_pass else 'FAIL'} (bar 1.0 mHa)",
              flush=True)
        m2[gname] = dict(E_VMC=float(Emean), E_VMC_err=float(Eerr), E_R12=float(E_R12),
                         diff_mHa=float(diff), g2_pass=bool(g2_pass), c=float(c_mix))
    return dict(gate6=dict(E=float(Eg6), err=float(eg6), dev=float(dev6), E_trunc=float(E_trunc)),
                m2=m2)


# =========================================================================== #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', default='all', choices=['m1', 'm2', 'all'])
    ap.add_argument('--gems', default='linexp,exp')
    ap.add_argument('--vmc-nwalk', type=int, default=3000)
    ap.add_argument('--vmc-nsweep', type=int, default=6000)
    ap.add_argument('--vmc-nburn', type=int, default=1200)
    ap.add_argument('--vmc-seeds', default='10,11')
    ap.add_argument('--probe', action='store_true')
    args = ap.parse_args()
    gems = [g.strip() for g in args.gems.split(',') if g.strip()]
    seeds = [int(s) for s in args.vmc_seeds.split(',') if s.strip()]

    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    print("#" * 100)
    print(f"# LiH MARRIAGE -- PHASE 3  (stage={args.stage}, gems={gems}, seeds={seeds})")
    print("#" * 100, flush=True)

    save = {}
    m1res = None
    if args.stage in ('m1', 'all'):
        m1res = run_m1(gems, T0, probe=args.probe)
        save['exact_ok'] = m1res['exact_ok']
        save['V_NN'] = m1res['V_NN']
        for g, r in m1res['m1'].items():
            save[f"m1_{g}_E_R12"] = r['E_R12']
            save[f"m1_{g}_c_mix"] = r['c_mix']
            save[f"m1_{g}_E0_grid"] = r['E0_grid']
            save[f"m1_{g}_Fbar"] = r['Fbar']
            save[f"m1_{g}_dcorr_mHa"] = r['dcorr_mHa']
            for k, v in r['pieces'].items():
                if v is not None:
                    save[f"m1_{g}_piece_{k}"] = v
        # persist m1 for a separate --stage m2 run
        np.savez(os.path.join(DATA, "lih_marriage_phase3_m1.npz"),
                 m1=np.array(m1res, dtype=object))

    if args.stage == 'm2' and m1res is None:
        z = np.load(os.path.join(DATA, "lih_marriage_phase3_m1.npz"), allow_pickle=True)
        m1res = z['m1'].item()
        save['exact_ok'] = m1res.get('exact_ok'); save['V_NN'] = m1res.get('V_NN')
        for g, r in m1res['m1'].items():        # carry M1 forward so the merged npz keeps both
            save[f"m1_{g}_E_R12"] = r['E_R12']; save[f"m1_{g}_c_mix"] = r['c_mix']
            save[f"m1_{g}_E0_grid"] = r['E0_grid']; save[f"m1_{g}_Fbar"] = r['Fbar']
            save[f"m1_{g}_dcorr_mHa"] = r['dcorr_mHa']

    if args.stage in ('m2', 'all'):
        m1 = m1res['m1'] if m1res else {}
        m2res = run_m2(m1, gems, T0, args.vmc_nwalk, args.vmc_nsweep, args.vmc_nburn,
                       seeds, probe=args.probe)
        save['gate6_E'] = m2res['gate6']['E']; save['gate6_dev'] = m2res['gate6']['dev']
        save['E_trunc'] = m2res['gate6']['E_trunc']
        for g, r in m2res.get('m2', {}).items():
            save[f"m2_{g}_E_VMC"] = r['E_VMC']; save[f"m2_{g}_E_VMC_err"] = r['E_VMC_err']
            save[f"m2_{g}_diff_mHa"] = r['diff_mHa']; save[f"m2_{g}_g2_pass"] = r['g2_pass']

    # ---------------- verdict ----------------
    if args.stage in ('all', 'm2') and not args.probe and m1res is not None:
        print("\n" + "#" * 100)
        print("# DECISION GATE")
        print("#" * 100)
        prim = 'linexp' if 'linexp' in m1res['m1'] else (gems[0] if gems else None)
        if prim and prim in m1res['m1']:
            r1 = m1res['m1'][prim]
            m2r = None
            for g, rr in save.items():
                pass
            E_R12 = r1['E_R12']
            g2 = save.get(f"m2_{prim}_g2_pass", None)
            beats = r1['beats_ceil']; above = r1['above_floor']
            if not above:
                verdict = "STOP (variational bound violated -- bug)"
            elif g2 and beats:
                verdict = "GO"
            elif g2 and not beats:
                verdict = "BORDERLINE (G2 agrees but E_R12 >= -8.032; widen the det truncation)"
            elif g2 is None:
                verdict = "M1 only (no G2 in this run)"
            else:
                verdict = "G2 FAIL -- investigate"
            print(f"  primary geminal = {prim}")
            print(f"  E_R12 = {E_R12:.6f}  (E_trunc {r1['E0_grid']:.6f}, corr {r1['dcorr_mHa']:+.3f} mHa)")
            print(f"  beats -8.032: {beats} ; above -8.0705: {above} ; G2 pass: {g2}")
            print(f"  VERDICT: {verdict}")
            save['verdict'] = verdict

    os.makedirs(DATA, exist_ok=True)
    np.savez(OUT, **save)
    print(f"\nsaved {OUT}")
    print(f"TOTAL wall {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
