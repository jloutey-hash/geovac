r"""pi/M=16 marriage at a chosen geminal range gamma (and optionally a finer kernel grid).

Runs the O(nnz) sparse analytic marriage (geovac.lih_r12ci.gvee_sparse.all_pieces_sparse) exactly as
debug/lih_pi_m16_energy_sparse.py does, at geminal range gamma set through the real API
geovac.lih_r12ci.configure_geminal (2026-10-04; it replaced the kernels.py GAM source override,
bit-identically against tests/data/lih_golden).  Only the optional kernel GRID (NXI, NETA, NPHI =
72, 44, 28) and the triangle cutoffs still go through geovac.lih_r12ci.source_override (a source exec
before the consumers import).  After import the driver ASSERTS that every namespace sees the new value
(P2.GAM, M.GAM, GV.GAM, hT/gT/triangle GAM, GV._KS keys).

    python debug/lih_exp_geminal_run.py --geminal exp --gam 2.0 --trunc 1e-2 [--grid 96 60 40]

Writes debug/data/lih_exp_geminal_<gem>_g<gam>_<NXI>x<NETA>x<NPHI>_<trunc>.npz.
"""
from __future__ import annotations

import os
import re
import sys
import time
import argparse
import importlib.util

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")


# install_kernels / install_triangle moved to geovac/lih_r12ci/source_override.py (LiH
# productionization step 2c, 2026-10-04); re-exported here so XR.install_kernels keeps working
# (install_kernels(gam, None) is now just geovac.lih_r12ci.configure_geminal(gam)).
from geovac.lih_r12ci.source_override import install_kernels, install_triangle  # noqa: E402,F401


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='exp', choices=['linexp', 'exp'])
    ap.add_argument('--gam', type=float, required=True)
    ap.add_argument('--trunc', type=float, default=1e-2)
    ap.add_argument('--grid', type=int, nargs=3, default=None, help='NXI NETA NPHI')
    ap.add_argument('--mem-cap-gb', type=float, default=6.0)
    ap.add_argument('--nphi-gen', type=int, default=None,
                    help='override the hard-coded nphi=28 of P2.build_kernel_gen (linexp kernels etc.)')
    ap.add_argument('--nphi-mode', type=int, default=None,
                    help='override nphi=48 of the mode/smooth kernels (P2.build_kernel_m_gen, '
                         'triangle.build_kernel_m, gVee.build_kernel_smooth -> _KS rebuilt)')
    ap.add_argument('--tri-mmax', type=int, default=None, help='triangle azimuthal mode cutoff (4)')
    ap.add_argument('--tri-lmax', type=int, default=None,
                    help='Neumann l cutoff (34): triangle LMAX, M._LMAX_TRI, P2.neumann_stack')
    ap.add_argument('--rank', type=int, default=None, help='triangle low-rank cutoff (30)')
    args = ap.parse_args()

    if args.grid is None:
        from geovac.lih_r12ci import configure_geminal
        configure_geminal(args.gam)
        import geovac.lih_r12ci.kernels as KG
    else:
        KG = install_kernels(args.gam, args.grid)      # grid source override, then configure_geminal
    if args.tri_mmax is not None or args.tri_lmax is not None:
        install_triangle(args.tri_mmax, args.tri_lmax)
    import geovac.lih_r12ci.enumerators as P2
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.gvee_sparse as S
    from geovac.lih_r12ci import gVee as GV, hT as HT, gT as GT, triangle as TR, hVee as HV

    # --- propagation asserts (every namespace that bound GAM / the grid at import) ---
    g = float(args.gam)
    for nm, v in [("KG", KG.GAM), ("P2", P2.GAM), ("M", M.GAM), ("GV", GV.GAM), ("HT", HT.GAM),
                  ("GT", GT.GAM), ("TR", TR.GAM), ("HV", HV.GAM)]:
        assert v == g, f"{nm}.GAM = {v} != {g}"
    assert g in GV._KS and 2 * g in GV._KS, f"GV._KS keys {list(GV._KS)}"
    assert M.NXI == KG.NXI and M.NETA == KG.NETA and M.NG == KG.NG == KG.NXI * KG.NETA
    assert P2.NG == KG.NG and TR.NG == KG.NG
    if args.nphi_gen is not None:
        P2.build_kernel_gen.__defaults__ = (args.nphi_gen,)
    if args.nphi_mode is not None:
        P2.build_kernel_m_gen.__defaults__ = (args.nphi_mode,)
        TR.build_kernel_m.__defaults__ = (args.nphi_mode,)
        assert M._bkm is TR.build_kernel_m
        GV._KS = {g: GV.build_kernel_smooth(g, nphi=args.nphi_mode),
                  2 * g: GV.build_kernel_smooth(2 * g, nphi=args.nphi_mode)}
    assert M.MMAX == P2.MMAX == TR.MMAX
    if args.tri_lmax is not None:
        assert TR.LMAX == args.tri_lmax
        M._LMAX_TRI = args.tri_lmax
        P2.neumann_stack.__defaults__ = (args.tri_lmax,)
    if args.rank is not None:
        _km, _rs = TR.Kmu_modes, M._rand_svd_middle
        R_ = args.rank

        def _kmu(Fm, mu_flat, m, key, rank=30, **kw):
            return _km(Fm, mu_flat, m, key, rank=R_, **kw)

        def _rsvd(FA, FB, gmu, rank=30, **kw):
            return _rs(FA, FB, gmu, rank=R_, **kw)
        M.Kmu_modes = _kmu; P2.Kmu_modes = _kmu; M._rand_svd_middle = _rsvd
    trunc_tag = f"_m{M.MMAX}_l{M._LMAX_TRI}_r{args.rank or 30}"
    nphi_gen = P2.build_kernel_gen.__defaults__[0]
    nphi_mode = TR.build_kernel_m.__defaults__[0]

    tag = (f"{args.geminal}_g{g:g}_{KG.NXI}x{KG.NETA}x{KG.NPHI}_p{nphi_gen}-{nphi_mode}"
           f"{trunc_tag if trunc_tag != '_m4_l34_r30' else ''}_{args.trunc:.0e}")
    out_path = os.path.join(DATA, f"lih_exp_geminal_{tag}.npz")
    if os.path.exists(out_path):
        out_path = out_path[:-4] + time.strftime("_%Y%m%d_%H%M%S") + ".npz"

    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    print("=" * 96)
    print(f"pi/M=16 marriage  geminal={args.geminal}  gamma={g}  grid {KG.NXI}x{KG.NETA} nphi {KG.NPHI}"
          f"  nphi gen/mode {nphi_gen}/{nphi_mode}  trunc={args.trunc:.0e}  (propagation asserts PASSED)")
    print("=" * 96, flush=True)
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=args.trunc)
    print(f"  reference: E_trunc={info['E_trunc']:.6f}  nact={refm.nact}  T nnz={len(refm.T_val)}",
          flush=True)
    gem = M.GeminalM(args.geminal)
    assert gem.g == g
    p = S.all_pieces_sparse(gem, refm, T0, mem_cap_gb=args.mem_cap_gb)
    E_R12, c, h = M.assemble_c(p, p['E0_grid'])
    corr = (E_R12 - info['E_trunc']) * 1e3
    print(f"\n  >>> E_R12 ({args.geminal}, gamma={g}) = {E_R12:.6f} Ha  (corr {corr:+.3f} mHa vs E_trunc "
          f"{info['E_trunc']:.6f}; c = {c:+.6f}; Fbar = {p['Fbar']:.6f}; E0_grid = {p['E0_grid']:.6f})",
          flush=True)
    out = dict(gam=g, geminal=args.geminal, trunc=args.trunc, NXI=KG.NXI, NETA=KG.NETA, NPHI=KG.NPHI,
               nphi_gen=nphi_gen, nphi_mode=nphi_mode, tri_mmax=M.MMAX, tri_lmax=M._LMAX_TRI,
               rank=args.rank or 30,
               E_trunc=info['E_trunc'], nact=refm.nact, E_R12=E_R12, c=c, h=h)
    for k, v in p.items():
        if isinstance(v, (int, float, np.floating)):
            out[k] = float(v)
    out['wall_s'] = time.time() - T0
    out['peak_rss_gb'] = S.peak_rss_gb()
    np.savez(out_path, **out)
    for k in sorted(out):
        if isinstance(out[k], float):
            print(f"    {k:>12s} = {out[k]:+.9f}")
    print(f"  saved {out_path}\n  peak RSS {out['peak_rss_gb']:.2f} GB  wall {out['wall_s']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
