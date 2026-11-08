r"""g_Vee defect audit, step C: per-triple analytic <f_p f_q / r_r> on the pi trunc reference via the
O(nnz) sparse path (debug/lih_gvee_sparse), triangles resolved per triple; then the SYMMETRY-ORBIT
test: |Psi|^2 is invariant under the electron relabelings (12), (34) and (13)(24) (singlet), so every
triple in an orbit must give the SAME value although each is reduced along a different path
(different kept pair / vertex order / dressing modes).  A path-dependent bug shows up as orbit spread.

    python debug/gvee_audit_pi_triples.py --geminal linexp --gam 0.25 --trunc 1e-2
Writes debug/data/gvee_audit_pi_<gem>_g<gam>_<trunc>.npz
"""
from __future__ import annotations

import os
import sys
import time
import argparse
import itertools

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def perm_pair(perm, p):
    return PAIRS6.index(tuple(sorted((perm[p[0]], perm[p[1]]))))


def orbits():
    gens = [{1: 2, 2: 1, 3: 3, 4: 4}, {1: 1, 2: 2, 3: 4, 4: 3}, {1: 3, 2: 4, 3: 1, 4: 2}]
    seen = {}
    orb = []
    for t in itertools.product(range(6), repeat=3):
        if t in seen:
            continue
        stack = [t]; members = {t}
        while stack:
            x = stack.pop()
            cands = [(x[1], x[0], x[2])]
            for g in gens:
                cands.append(tuple(perm_pair(g, PAIRS6[i]) for i in x))
            for y in cands:
                if y not in members:
                    members.add(y); stack.append(y)
        for m in members:
            seen[m] = len(orb)
        orb.append(sorted(members))
    return orb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp')
    ap.add_argument('--gam', type=float, default=0.25)
    ap.add_argument('--trunc', type=float, default=1e-2)
    ap.add_argument('--mem-cap-gb', type=float, default=4.0)
    ap.add_argument('--ref', default='pi', choices=['pi', 'sigma74'])
    a = ap.parse_args()
    import lih_exp_geminal_run as RUN
    KG = RUN.install_kernels(a.gam, None)
    import lih_marriage_phase2 as P2
    import lih_pi_m16 as M
    import lih_gvee_sparse as S
    KG.USE_EXACT_NEUMANN = True
    t0 = time.time()
    if a.ref == 'pi':
        ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
        refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=a.trunc,
                                            verbose=False)
    else:
        refm = M._wrap_p2_ref(P2.build_ref_artifact())
    gem = M.GeminalM(a.geminal)
    Fbar = M.fbar(gem, refm)
    ctx = S.Ctx(gem, refm, a.mem_cap_gb)
    col = S.TriCollector()
    Vee_r = np.array([S.vee_reduce_sparse([], r, 'coul', None, ctx, col) for r in PAIRS6])
    FV = np.zeros((6, 6))
    for i, p in enumerate(PAIRS6):
        for k, r in enumerate(PAIRS6):
            if p == r:
                FV[i, k] = S.vee_reduce_sparse([], r, 'yuk', gem.g, ctx, col)
            else:
                FV[i, k] = S.vee_reduce_sparse([tuple(p)], r, 'coul', None, ctx, col)
    assert not col.recs
    F2V = np.zeros((6, 6, 6))
    tid_of_rec = []
    for i, p in enumerate(PAIRS6):
        for j, q in enumerate(PAIRS6):
            for k, r in enumerate(PAIRS6):
                n0 = len(col.recs)
                rem_f, kk, gk = S._triple_args(p, q, r, gem)
                F2V[i, j, k] = S.vee_reduce_sparse(rem_f, r, kk, gk, ctx, col)
                tid_of_rec += [(i * 6 + j) * 6 + k] * (len(col.recs) - n0)
    print(f"  non-triangle done [{time.time()-t0:.0f}s]", flush=True)
    # triangles, per triple
    im = np.concatenate([x[0] for x in col.recs]); il = np.concatenate([x[1] for x in col.recs])
    ir = np.concatenate([x[2] for x in col.recs]); w = np.concatenate([x[3] for x in col.recs])
    tid = np.concatenate([np.full(len(x[0]), t) for x, t in zip(col.recs, tid_of_rec)])
    n = ctx.nact
    code = (im * n + il) * n + ir
    uc, inv = np.unique(code, return_inverse=True)
    vals = np.zeros(len(uc))
    ucm = uc // (n * n); uil = (uc // n) % n; uir = uc % n
    for cm in col.order:
        s = ucm == cm
        vals[s] = S.tri_entries(ctx, cm, uil[s], uir[s])
    tri_t = np.bincount(tid, weights=w * vals[inv], minlength=216).reshape(6, 6, 6)
    F2V = F2V + tri_t
    Vee = Vee_r.sum(); FVee = FV.sum(); F2Vee = F2V.sum()
    g = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    print(f"{a.geminal} g={a.gam} {a.ref}: Fbar={Fbar:.10f} Vee={Vee:.10f} FVee={FVee:.10f} F2Vee={F2Vee:.10f} "
          f"g_Vee={g:.10f}  [{time.time()-t0:.0f}s]", flush=True)
    # orbit test
    print("\n  symmetry-orbit spread (max-min within each orbit of equivalent triples):")
    worst = []
    for ob in orbits():
        v = np.array([F2V[t] for t in ob])
        worst.append((v.max() - v.min(), ob, v))
    worst.sort(key=lambda x: -x[0])
    for spread, ob, v in worst[:15]:
        lab = " ".join(f"{PAIRS6[t[0]]}{PAIRS6[t[1]]}/{PAIRS6[t[2]]}" for t in ob[:2])
        print(f"   spread {spread:.3e}  n={len(ob):2d}  mean {v.mean():+.8f}  e.g. {lab}")
    print("  FV orbits:")
    for i in range(6):
        pass
    tag = f"{a.ref}_{a.geminal}_g{a.gam:g}_{a.trunc:.0e}"
    np.savez(os.path.join(DATA, f"gvee_audit_pi_{tag}.npz"), Fbar=Fbar, Vee_r=Vee_r, FV=FV, F2V=F2V,
             tri=tri_t, g_Vee=g)


if __name__ == "__main__":
    main()
