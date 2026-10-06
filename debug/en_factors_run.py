r"""Driver: analytic e-n extension of the CI-relaxed LiH eigenproblem (geovac/lih_r12ci/en_factors.py).

Phase 1: {chi_I} U {(F-Fbar) psi} U {W1 psi, W2 psi}, W_k from the parameter-free core switch.
Frozen 2x2 pieces from the golden oracle (as debug/ci_relaxed_run.py).

Usage:  python debug/en_factors_run.py 1e-2 [--switch dens|frac] [--wone] [--save debug/data/en_factors_1e-2.npz]
  --wone : identity mode, w == 1 (lap w == 0): W1 psi = 4 psi, W2 psi = 6 psi -> every new column must
           collapse onto the reference column (S_I = 4 d, H_I = 4 H_CI d, ...).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trunc", type=float)
    ap.add_argument("--switch", default="dens", choices=["dens", "frac"])
    ap.add_argument("--wone", action="store_true")
    ap.add_argument("--phase2", action="store_true")
    ap.add_argument("--save", default=None)
    a = ap.parse_args()
    from geovac.lih_r12ci import configure_geminal
    configure_geminal(1.5)
    import numpy as np
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.ci_relaxed as CR
    import geovac.lih_r12ci.en_factors as EN
    from ci_relaxed_run import golden_pieces
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=a.trunc,
                                        verbose=False)
    gem = M.GeminalM('exp')
    gold = golden_pieces(a.trunc)
    pieces = gold["pieces"]
    print(f"trunc={a.trunc} N_det={len(info['kept'])}  [{time.time()-T0:.1f}s]", flush=True)
    cr = CR.ci_relaxed_energy(gem, refm, info, pieces, verbose=True)
    print(f"  ci_relaxed E_A = {cr['E_A']:.10f}   [{time.time()-T0:.1f}s]", flush=True)
    if a.wone:
        w, lw, meta = EN.switch_on_grid(info, a.switch)
        w = np.ones_like(w); lw = np.zeros_like(lw)
        en = EN.en_columns(gem, info, float(pieces['Fbar']), w, lw, verbose=True)
        d = en['d']; N = en['N']
        Hd = cr['H_CI'] @ d; Sd = cr['S_CI_grid'] @ d
        print("w == 1 identity checks (should be ~0):")
        for k, (mult) in enumerate((4.0, 6.0)):
            print(f"  W{k+1}: |S_I - {mult} S d| = {np.max(np.abs(en['S_I'][:, k] - mult * Sd)):.2e}"
                  f"  |H_I - {mult} H d| = {np.max(np.abs(en['H_I'][:, k] - mult * Hd)):.2e}"
                  f"  |LA_I| = {np.max(np.abs(en['LA_I'][:, k])):.2e}")
            print(f"       S_F = {en['S_F'][k]:+.2e}   H_F - {mult} d.t = {en['H_F'][k] - mult * float(d @ cr['t']):+.2e}")
        E0 = float(d @ Hd)
        print(f"  S_ab - [16,24;24,36] = {np.max(np.abs(en['S_ab'] - np.array([[16, 24], [24, 36]]))):.2e}")
        print(f"  H_ab - [16,24;24,36] E0 = {np.max(np.abs(en['H_ab'] - np.array([[16, 24], [24, 36]]) * E0)):.2e}")
        print(f"wall {time.time()-T0:.1f}s")
        return 0
    r = EN.en_energy(gem, refm, info, pieces, switch=a.switch, verbose=True, cr=cr)
    print(f"\nswitch: {r['switch']}")
    print("identity diagnostics:")
    for k, v in r['identity'].items():
        print(f"  {k:22s} {v:+.3e}")
    print("energies (total, Ha):")
    for k, v in r['E'].items():
        print(f"  {k:8s} {v:.10f}   vs A: {(v - r['E']['A'])*1e3:+.4f} mHa")
    en = r['en']
    print("parts_F:", {k: {kk: round(vv, 8) for kk, vv in p.items()} for k, p in en['parts_F'].items()})
    print("parts_ab:", {str(k): {kk: round(vv, 8) for kk, vv in p.items()} for k, p in en['parts_ab'].items()})
    print(f"wall {time.time()-T0:.1f}s  en columns {en['wall_s']:.1f}s")
    try:
        import psutil
        mi = psutil.Process().memory_info()
        print(f"peak RSS {getattr(mi, 'peak_wset', mi.rss)/1e9:.2f} GB")
    except Exception:
        pass
    if a.save:
        np.savez(a.save, S=r['S'], H=r['H'], E=json.dumps(r['E']), identity=json.dumps(r['identity']),
                 switch=json.dumps(r['switch']), S_I=en['S_I'], H_I=en['H_I'], T_I=en['T_I'],
                 Vne_I=en['Vne_I'], Vee_I=en['Vee_I'], S_F=en['S_F'], H_F=en['H_F'], S_ab=en['S_ab'],
                 H_ab=en['H_ab'], d=en['d'], w=r['w'], lapw=r['lapw'],
                 parts_F=json.dumps({str(k): v for k, v in en['parts_F'].items()}),
                 parts_ab=json.dumps({str(k): v for k, v in en['parts_ab'].items()}),
                 E_A_cr=cr['E_A'], vec_A_W12=r['vecs']['A_W12'][1])
    if a.phase2:
        r2 = EN.en2_energy(gem, refm, info, pieces, verbose=True, ph1=r)
        print("\nPHASE 2 energies (total, Ha):")
        for k, v in r2['E'].items():
            print(f"  {k:8s} {v:.10f}   vs A: {(v - r2['E']['A'])*1e3:+.4f} mHa")
        e2 = r2['en2']
        print(f"  variants {e2['variants']}")
        print(f"wall {time.time()-T0:.1f}s  en2 columns {e2['wall_s']:.1f}s")
        try:
            import psutil
            mi = psutil.Process().memory_info()
            print(f"peak RSS {getattr(mi, 'peak_wset', mi.rss)/1e9:.2f} GB")
        except Exception:
            pass
        if a.save:
            base = a.save[:-4] if a.save.endswith('.npz') else a.save
            np.savez(base + "_phase2.npz", S=r2['S'], H=r2['H'], E=json.dumps(r2['E']),
                     S_I=e2['S_I'], H_I=e2['H_I'], HPG=e2['HPG'], SPG=e2['SPG'], d=en['d'],
                     blk=json.dumps({f"{k[0]}|{k[1]}": v for k, v in e2['blk'].items()}),
                     parts_I=json.dumps({str(k): {kk: vv.tolist() for kk, vv in v.items()}
                                         for k, v in e2['parts_I'].items()}),
                     vec_A_W12_G=r2['vecs']['A_W12_G'][1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
