r"""DIAGNOSTIC (measure only) -- can g_Vee in the LiH "marriage" explicit-r12 engine be made much
cheaper by algebra BEFORE anyone builds the brute O(nnz) rewrite?  (2026-09-25)

g_Vee = <F^2 V_ee> - 2 Fbar <F V_ee> + Fbar^2 <V_ee>, evaluated term-by-term (216 + 36 + 6 edge
multisets over 4 electrons) by the production reducers (lih_pi_m16.vee_reduce_m for the pi reference,
lih_marriage_phase2.vee_reduce for the sigma 74-det reference).  Nothing in geovac/ or the two
production drivers is edited -- timing/memory instrumentation is by monkeypatching module attributes
IN THIS PROCESS ONLY.

Subcommands (run ONE heavy job at a time; each appends to debug/data/lih_gvee_diag.log):
  d0 --ref pi|sigma   D0 cost/memory profile by term topology + D1 cheap pair-space route
  d2                  D2 derivative identities + gamma-series / exact-power-route convergence
  nnz                 T-sparsity probe: distinct index triples vs nact^3 (O(nnz) memory estimate)
  merge               combine the part files into debug/data/lih_gvee_diag.npz

Run from repo root (Git Bash), e.g.
  python debug/lih_gvee_diag.py d0 --ref pi >> debug/data/lih_gvee_diag.log 2>&1; echo $? > debug/data/lih_gvee_diag_d0pi.exit
"""
from __future__ import annotations

import os
import sys
import time
import argparse
import tracemalloc
from collections import defaultdict
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

import geovac.lih_r12ci.kernels as KG
import lih_marriage_phase2 as P2
import lih_pi_m16 as M

PAIRS6 = P2.PAIRS6
PART = lambda name: os.path.join(DATA, f"lih_gvee_diag_{name}.npz")


def _t(t0: float) -> str:
    return f"[{time.time() - t0:7.1f}s]"


# =========================================================================== #
# exclusive-time profiler (monkeypatch wrappers)
# =========================================================================== #
class Prof:
    def __init__(self) -> None:
        self.excl: Dict[str, float] = defaultdict(float)
        self.stack: List[float] = []
        self.max_field_bytes = 0

    def wrap(self, owner, name: str, label: str) -> None:
        orig = getattr(owner, name)
        prof = self

        def w(*a, **k):
            t0 = time.perf_counter(); prof.stack.append(0.0)
            try:
                out = orig(*a, **k)
            finally:
                el = time.perf_counter() - t0; child = prof.stack.pop()
                prof.excl[label] += el - child
                if prof.stack:
                    prof.stack[-1] += el
            if name == 'cf_mul':
                prof.max_field_bytes = max(prof.max_field_bytes, int(out[0].nbytes))
            return out
        setattr(owner, name, w)

    def snap(self) -> Dict[str, float]:
        return dict(self.excl)


def install_prof() -> Prof:
    pr = Prof()
    pr.wrap(P2, 'dress_mat', 'kernel_matvec(m0 NGxNG)')
    pr.wrap(M.GeminalM, 'dress_mode', 'kernel_matvec(mode-m NGxNG)')
    pr.wrap(P2, 'neumann_stack', 'neumann(m0)')
    pr.wrap(P2, 'coul_mode_stack', 'neumann(mode-m)')
    pr.wrap(M, 'coul_mode_stack_any', 'neumann(mode-m,any)')
    pr.wrap(M, 'Kmu_modes', 'tri: Kmu randomized eigh (NG^2 x rank)')
    pr.wrap(P2, 'Kmu_modes', 'tri: Kmu randomized eigh (NG^2 x rank)')
    pr.wrap(M, '_rand_svd_middle', 'tri: Kmu randomized SVD (A!=B)')
    pr.wrap(M, 'get_tri_matrix_m', 'tri: matrix assembly')
    pr.wrap(P2, 'get_tri_matrix', 'tri: matrix assembly')
    pr.wrap(P2, 'cf_mul', 'compound-field product (cf_mul)')
    return pr


# =========================================================================== #
# term enumeration + topology classification
# =========================================================================== #
def classify(fedges: Sequence[Tuple[int, int]], r: Tuple[int, int]) -> str:
    """Topology of the edge multiset {f-edges} + {coul r} on 4 electrons.
    PAIR-prod / PAIR-disc : every connected component spans ONE electron pair (pair-space W contraction)
    PROD3  : an f sits on r and another f hangs off r (3-body path, Yukawa kept kernel)
    TREE3  : no f on r, connected 3-vertex tree (share-vertex)
    TREE4  : 4-vertex tree that is not the bridge (star / chain through a kept vertex)
    BRIDGE : f_ab . coul_bc . f_cd  (disjoint f's joined by the Coulomb)
    TRI    : 3-cycle f_ik f_jk / r_ij  (the RI-free triangle)"""
    R = frozenset(r); fs = [frozenset(x) for x in fedges]
    on_r = sum(1 for f in fs if f == R)
    rem = [f for f in fs if f != R]
    # connected components of all distinct edges
    edges = set(fs) | {R}
    comp: List[set] = []
    for ed in edges:
        hit = [c for c in comp if c & ed]
        new = set(ed)
        for c in hit:
            new |= c; comp.remove(c)
        comp.append(new)
    if all(len(c) == 2 for c in comp):
        return 'PAIR-prod' if on_r else 'PAIR-disc'
    if len(rem) == 1:
        return 'PROD3' if on_r else 'TREE3'
    a_, b_ = rem
    if a_ == b_:                                   # f^2 on a pair adjacent to r
        return 'TREE3'
    verts = a_ | b_ | R
    if len(verts) == 3:
        return 'TRI' if (a_ != R and b_ != R and len({a_, b_, R}) == 3) else 'TREE3'
    if not (a_ & b_) and (a_ & R) and (b_ & R):
        return 'BRIDGE'
    return 'TREE4'


COARSE = {'PAIR-disc': '(a) separable', 'TREE3': '(a) separable', 'TREE4': '(a) separable',
          'PAIR-prod': '(b) same-pair product', 'PROD3': '(b) same-pair product',
          'TRI': '(c) triangle', 'BRIDGE': '(d) bridge'}


def enumerate_terms() -> List[Tuple[str, tuple, Tuple[int, int]]]:
    out = [('V', (), r) for r in PAIRS6]
    out += [('FV', (p,), r) for p in PAIRS6 for r in PAIRS6]
    out += [('F2', (p, q), r) for p in PAIRS6 for q in PAIRS6 for r in PAIRS6]
    return out


def term_args(fedges, r, g: float):
    rset = frozenset(r); fs = [frozenset(x) for x in fedges]
    n_on_r = sum(1 for ed in fs if ed == rset)
    rem_f = [tuple(sorted(ed)) for ed in fs if ed != rset]
    kk = 'coul' if n_on_r == 0 else 'yuk'
    gam_k = None if n_on_r == 0 else (g if n_on_r == 1 else 2 * g)
    return rem_f, kk, gam_k


# =========================================================================== #
# references
# =========================================================================== #
def load_ref(which: str, trunc: float = 1e-2):
    if which == 'pi':
        ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
        refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi',
                                            trunc=trunc, verbose=False)
        z = np.load(os.path.join(DATA, "lih_pi_m16_energy_trunc1e-2.npz"))
        pieces = {k[len('linexp_'):]: float(z[k]) for k in z.files if k.startswith('linexp_')
                  and z[k].ndim == 0}
        return refm, info, pieces, None
    p2ref = P2.build_ref_artifact()
    refm = M._wrap_p2_ref(p2ref)
    z = np.load(os.path.join(DATA, "lih_marriage_phase3_m1.npz"), allow_pickle=True)
    pieces = dict(z['m1'].item()['m1']['linexp']['pieces'])
    return refm, None, pieces, p2ref


def e_sensitivity(pieces: dict) -> Tuple[float, float]:
    """E_R12 and dE_R12/dg_Vee (central difference) from the banked pieces."""
    E0 = pieces['E0_grid']
    E, _, _ = M.assemble_c(pieces, E0)
    d = 1e-6
    pp = dict(pieces); pp['g_Vee'] = pieces['g_Vee'] + d
    pm = dict(pieces); pm['g_Vee'] = pieces['g_Vee'] - d
    return E, (M.assemble_c(pp, E0)[0] - M.assemble_c(pm, E0)[0]) / (2 * d)


# =========================================================================== #
# D1 cheap route: pair-space W matrices + O(nnz) T contraction
# =========================================================================== #
def pair_W(gem: "M.GeminalM", refm: "M.RefM", nf: int, has_c: bool) -> np.ndarray:
    """W[A,B] = INT rho_A(1) rho_B(2) K(r12) with K = f^nf * (1/r12)^has_c, selection mA+mB=0.
    Built with the SAME primitives the reducers use (dress_rows_m / _kept_pot), O(nact) solves."""
    rho, mA, g = refm.rho, refm.mA, gem.g
    if not has_c:
        psi = M.dress_rows_m(gem, 'f' if nf == 1 else ('f', 'f'), rho, mA)
    else:
        kk = 'coul' if nf == 0 else 'yuk'
        gk = None if nf == 0 else (g if nf == 1 else 2 * g)
        psi = np.empty_like(rho)
        am = np.abs(mA)
        for mode in sorted(set(am.tolist())):
            rows = am == mode
            psi[rows] = M._kept_pot(gem, rho[rows], kk, gk, int(mode))
    W = (rho * P2.GI_W[None, :]) @ psi.T
    W *= (mA[:, None] + mA[None, :]) == 0
    return W


def pair_value(fedges, r, Wc: Dict[Tuple[int, bool], np.ndarray], refm) -> float:
    idx = {1: refm.T_iA, 2: refm.T_iB, 3: refm.T_iC, 4: refm.T_iD}
    cnt: Dict[frozenset, List[int]] = defaultdict(lambda: [0, 0])
    for f in fedges:
        cnt[frozenset(f)][0] += 1
    cnt[frozenset(r)][1] += 1
    w = refm.T_val.copy()
    used = set()
    for pr, (nf, nc) in cnt.items():
        e1, e2 = sorted(pr)
        w = w * Wc[(nf, bool(nc))][idx[e1], idx[e2]]
        used |= set(pr)
    for e in (1, 2, 3, 4):
        if e not in used:
            w = w * refm.d[idx[e]]
    return float(w.sum())


# =========================================================================== #
# D0 / D1 runner
# =========================================================================== #
def run_d0(which: str) -> None:
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    print("=" * 100)
    print(f"D0/D1 -- g_Vee cost profile by topology, reference = {which}  (linexp, EXACT path)")
    print("=" * 100, flush=True)
    refm, info, pieces, p2ref = load_ref(which)
    gem = M.GeminalM('linexp')
    Fbar = pieces['Fbar']
    E_R12, dEdg = e_sensitivity(pieces)
    print(f"  ref: nact={refm.nact}  T nnz={len(refm.T_val)}  m!=0 pairs={int(np.sum(refm.mA != 0))}  "
          f"Fbar(banked)={Fbar:.6f}  E_R12(banked 2x2)={E_R12:.6f}  dE/dg_Vee={dEdg:+.5f}  {_t(T0)}",
          flush=True)
    pr = install_prof()
    reducer = (lambda rem_f, r, kk, gk: P2.vee_reduce(rem_f, r, kk, gk, gem, p2ref)) if which == 'sigma' \
        else (lambda rem_f, r, kk, gk: M.vee_reduce_m(rem_f, r, kk, gk, gem, refm))
    terms = enumerate_terms()
    rec = []
    tracemalloc.start()
    base_mem = tracemalloc.get_traced_memory()[0]
    for n, (kind, fed, r) in enumerate(terms):
        rem_f, kk, gk = term_args(fed, r, gem.g)
        cls = classify(fed, r)
        s0 = pr.snap(); pr.max_field_bytes = 0
        tracemalloc.reset_peak()
        m0 = tracemalloc.get_traced_memory()[0]
        t0 = time.perf_counter()
        val = reducer(rem_f, r, kk, gk)
        dt = time.perf_counter() - t0
        peak = tracemalloc.get_traced_memory()[1] - m0
        s1 = pr.snap()
        dprof = {k: s1.get(k, 0.0) - s0.get(k, 0.0) for k in s1}
        wg = 1.0 if kind == 'F2' else (-2 * Fbar if kind == 'FV' else Fbar ** 2)
        rec.append(dict(kind=kind, fed=fed, r=r, cls=cls, val=val, g=wg * val, t=dt, peak=peak,
                        field=pr.max_field_bytes, prof=dprof))
        if dt > 5 or n % 40 == 0:
            print(f"    [{n:3d}] {kind:2s} f={fed} r={r} {cls:9s} val={val:+.8f} t={dt:7.1f}s "
                  f"peak={peak/1e9:6.2f}GB field={pr.max_field_bytes/1e9:6.2f}GB {_t(T0)}", flush=True)
    tracemalloc.stop()
    Vee = sum(x['val'] for x in rec if x['kind'] == 'V')
    FVee = sum(x['val'] for x in rec if x['kind'] == 'FV')
    F2Vee = sum(x['val'] for x in rec if x['kind'] == 'F2')
    g = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    ttot = sum(x['t'] for x in rec)
    print(f"\n  CONSISTENCY: Vee={Vee:+.8f} (banked {pieces['Vee']:+.8f})  FVee={FVee:+.8f} "
          f"(banked {pieces['FVee']:+.8f})  F2Vee={F2Vee:+.8f} (banked {pieces['F2Vee']:+.8f})")
    print(f"  g_Vee recomputed = {g:+.8f}  banked = {pieces['g_Vee']:+.8f}  d = {(g-pieces['g_Vee']):+.2e} Ha"
          f"   total reducer wall = {ttot:.1f}s", flush=True)

    # ---- D0 table ----
    print("\n  D0 -- by topology class (fine):")
    print(f"  {'class':10s} {'coarse':22s} {'nterm':>5s} {'wall_s':>8s} {'%wall':>6s} {'peakGB':>7s} "
          f"{'fieldGB':>7s} {'F2Vee':>11s} {'FVee':>11s} {'Vee':>9s} {'g contrib':>11s} {'mHa(2x2)':>9s}")
    classes = ['PAIR-disc', 'PAIR-prod', 'PROD3', 'TREE3', 'TREE4', 'BRIDGE', 'TRI']
    table = {}
    for c in classes:
        xs = [x for x in rec if x['cls'] == c]
        if not xs:
            continue
        row = dict(n=len(xs), t=sum(x['t'] for x in xs), peak=max(x['peak'] for x in xs),
                   field=max(x['field'] for x in xs),
                   F2=sum(x['val'] for x in xs if x['kind'] == 'F2'),
                   FV=sum(x['val'] for x in xs if x['kind'] == 'FV'),
                   V=sum(x['val'] for x in xs if x['kind'] == 'V'),
                   g=sum(x['g'] for x in xs))
        row['mHa'] = row['g'] * dEdg * 1e3
        table[c] = row
        print(f"  {c:10s} {COARSE[c]:22s} {row['n']:5d} {row['t']:8.1f} {100*row['t']/ttot:6.1f} "
              f"{row['peak']/1e9:7.2f} {row['field']/1e9:7.2f} {row['F2']:+11.5f} {row['FV']:+11.5f} "
              f"{row['V']:+9.5f} {row['g']:+11.6f} {row['mHa']:+9.3f}")
    print("\n  D0 -- by coarse class:")
    coarse_tab = {}
    for cc in sorted(set(COARSE.values())):
        xs = [x for x in rec if COARSE[x['cls']] == cc]
        if not xs:
            continue
        row = dict(t=sum(x['t'] for x in xs), peak=max(x['peak'] for x in xs),
                   field=max(x['field'] for x in xs), g=sum(x['g'] for x in xs))
        coarse_tab[cc] = row
        print(f"    {cc:24s} wall {row['t']:8.1f}s ({100*row['t']/ttot:5.1f}%)  peak {row['peak']/1e9:6.2f} GB"
              f"  max field {row['field']/1e9:6.2f} GB  g {row['g']:+.6f} Ha = {row['g']*dEdg*1e3:+.3f} mHa")
    # dense-field owners
    big = sorted(rec, key=lambda x: -x['field'])[:12]
    print("\n  largest compound fields (the dense nact^k x NG owners):")
    for x in big:
        print(f"    {x['kind']:2s} f={x['fed']} r={x['r']} {x['cls']:9s} field={x['field']/1e9:6.2f} GB "
              f"peak={x['peak']/1e9:6.2f} GB t={x['t']:.1f}s")
    # primitive profile by class
    print("\n  exclusive primitive time by class (s):")
    labs = sorted({k for x in rec for k in x['prof']})
    for c in classes:
        xs = [x for x in rec if x['cls'] == c]
        if not xs:
            continue
        agg = {k: sum(x['prof'].get(k, 0.0) for x in xs) for k in labs}
        s = "  ".join(f"{k}={v:.1f}" for k, v in agg.items() if v > 0.05)
        print(f"    {c:10s} {s}")
    aggall = {k: sum(x['prof'].get(k, 0.0) for x in rec) for k in labs}
    print("    ALL        " + "  ".join(f"{k}={v:.1f}" for k, v in aggall.items() if v > 0.05))

    # ---- D1 ----
    print("\n  D1 -- cancellation:")
    ratio = abs(F2Vee) / abs(g)
    print(f"    |F2Vee|/|g_Vee| = {ratio:.1f}   |2 Fbar FVee|/|g| = {abs(2*Fbar*FVee)/abs(g):.1f}   "
          f"|Fbar^2 Vee|/|g| = {abs(Fbar**2*Vee)/abs(g):.1f}")
    for c, row in table.items():
        gross = row['F2'] + abs(2 * Fbar * row['FV']) + abs(Fbar ** 2 * row['V'])
        print(f"    {c:10s} gross |parts| {abs(row['F2'])+abs(2*Fbar*row['FV'])+abs(Fbar**2*row['V']):10.5f}"
              f"  net g {row['g']:+.6f}  (net/gross {abs(row['g'])/max(gross,1e-30):.3f})")
    # cheap pair-space route on PAIR-* terms
    t1 = time.perf_counter()
    Wc = {}
    for key in [(1, False), (2, False), (0, True), (1, True), (2, True)]:
        Wc[key] = pair_W(gem, refm, key[0], key[1])
    tW = time.perf_counter() - t1
    for key, Wm in Wc.items():
        print(f"      W{key}: max|W|={np.max(np.abs(Wm)):.3e}  max|W-W^T|={np.max(np.abs(Wm - Wm.T)):.2e}")
    maxdev = 0.0; gcheap = 0.0; gprod = 0.0; t2 = time.perf_counter()
    for x in rec:
        if x['cls'].startswith('PAIR'):
            v = pair_value(x['fed'], x['r'], Wc, refm)
            wg = x['g'] / x['val'] if x['val'] != 0 else (1.0 if x['kind'] == 'F2' else
                                                          (-2 * Fbar if x['kind'] == 'FV' else Fbar ** 2))
            maxdev = max(maxdev, abs(v - x['val']))
            gcheap += wg * v; gprod += x['g']
    tC = time.perf_counter() - t2
    tpair = sum(x['t'] for x in rec if x['cls'].startswith('PAIR'))
    print(f"    cheap pair-space route (5 W matrices, {refm.nact}^2 each + O(nnz) contraction):")
    print(f"      W build {tW:.1f}s + contraction {tC:.2f}s  vs production on PAIR-* terms {tpair:.1f}s "
          f"({100*tpair/ttot:.1f}% of g_Vee wall)")
    print(f"      max |cheap - production| per term = {maxdev:.2e} Ha ;  g contrib cheap {gcheap:+.10f} vs "
          f"prod {gprod:+.10f} (d {gcheap-gprod:+.2e} Ha)")
    tri_t = sum(x['t'] for x in rec if x['cls'] == 'TRI')
    tri_g = sum(x['g'] for x in rec if x['cls'] == 'TRI')
    junc_t = sum(x['t'] for x in rec if x['cls'] in ('PROD3', 'TREE3', 'TREE4', 'BRIDGE'))
    print(f"    wall share: PAIR {100*tpair/ttot:.1f}%  |  junction/tree/bridge (O(nnz) w/o triangle) "
          f"{100*junc_t/ttot:.1f}%  |  TRI {100*tri_t/ttot:.1f}%")
    print(f"    TRI contribution to g = {tri_g:+.6f} Ha = {tri_g*dEdg*1e3:+.3f} mHa in E_R12", flush=True)

    np.savez(PART(f"d0{which}"), which=which, nact=refm.nact, nnz=len(refm.T_val), Fbar=Fbar,
             dEdg=dEdg, E_R12=E_R12, g=g, g_banked=pieces['g_Vee'], Vee=Vee, FVee=FVee, F2Vee=F2Vee,
             ttot=ttot, rec=np.array(rec, dtype=object), table=np.array(table, dtype=object),
             coarse=np.array(coarse_tab, dtype=object), cheap_maxdev=maxdev, cheap_tW=tW, cheap_tC=tC,
             cheap_g=gcheap, prod_pair_g=gprod, prof_all=np.array(aggall, dtype=object))
    print(f"  saved {PART('d0'+which)}  {_t(T0)}", flush=True)


# =========================================================================== #
# D2 -- derivative identities, gamma series, exact-power routes
# =========================================================================== #
def _kernel_rows(func: Callable, rows: np.ndarray, nphi: int = 28) -> np.ndarray:
    """rows-subset of the axially averaged kernel (same quadrature as build_kernel_gen)."""
    from numpy.polynomial.legendre import leggauss
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = P2.rho_cyl_f; z = P2.zc_f
    rc2 = rc[rows, None] * rc[None, :]
    base = rc[rows, None] ** 2 + rc[None, :] ** 2 + (z[rows, None] - z[None, :]) ** 2
    K = np.zeros((len(rows), rc.size))
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 0.0))
        K += 2.0 * w * func(d)
    return K


def _bil(L: np.ndarray, K: np.ndarray, R: np.ndarray) -> np.ndarray:
    return np.sum((L @ K) * R, axis=1)


def grid_power_moments(vecsL: np.ndarray, vecsR: np.ndarray, nmax: int, gam: float,
                       nphi: int = 28) -> Dict[str, np.ndarray]:
    """For each pair k: W_n[k] = INT INT L_k R_k r12^n (n=0..nmax), plus the reference kernels
    W_lin = r e^{-g r}, W_exp = e^{-g r}, and the |.|-bounds with e^{+g r}.  Full NG x NG, per phi."""
    from numpy.polynomial.legendre import leggauss
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = P2.rho_cyl_f; z = P2.zc_f
    rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    pref = 2 * np.pi * P2.A3 * P2.A3                       # grid_int(L * a^3 K @ geo R)
    gL = vecsL * P2.GEO[None, :]; gR = vecsR * P2.GEO[None, :]
    aL = np.abs(gL); aR = np.abs(gR)
    nk = vecsL.shape[0]
    Wn = np.zeros((nmax + 1, nk)); Wabs = np.zeros((nmax + 1, nk))
    Wlin = np.zeros(nk); Wexp = np.zeros(nk); Blin = np.zeros(nk); Bexp = np.zeros(nk)
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 0.0))
        ww = 2.0 * w * pref
        e = np.exp(-gam * d)
        Wexp += ww * _bil(gL, e, gR)
        Wlin += ww * _bil(gL, d * e, gR)
        ep = np.exp(gam * d)
        Bexp += ww * _bil(aL, ep, aR)
        Blin += ww * _bil(aL, d * ep, aR)
        Dn = np.ones_like(d)
        for n in range(nmax + 1):
            if n:
                Dn *= d
            Wn[n] += ww * _bil(gL, Dn, gR)
            Wabs[n] += ww * _bil(aL, Dn, aR)
    return dict(Wn=Wn, Wabs=Wabs, Wlin=Wlin, Wexp=Wexp, Blin=Blin, Bexp=Bexp)


def series_table(Wn: np.ndarray, Wabs: np.ndarray, ref: float, gam: float, shift: int,
                 use_mp: bool = False) -> Dict[str, object]:
    """partial sums S_N = sum_{k<=N} (-g)^k/k! W_{k+shift}; terms needed for 1e-6 rel (and staying)."""
    import mpmath as mp
    mp.mp.dps = 40
    nmax = Wn.shape[0] - 1 - shift
    S = mp.mpf(0) if use_mp else 0.0
    parts, terms, absterms = [], [], []
    fact = 1.0
    for k in range(nmax + 1):
        if k:
            fact *= k
        tk = (mp.mpf(-gam) ** k / mp.factorial(k) * mp.mpf(Wn[k + shift])) if use_mp else \
            ((-gam) ** k / fact * Wn[k + shift])
        S = S + tk
        parts.append(float(S)); terms.append(float(tk))
        absterms.append(gam ** k / fact * Wabs[k + shift])
    rel = np.abs(np.array(parts) - ref) / abs(ref)
    ok = np.nonzero(rel <= 1e-6)[0]
    need = None
    for k in ok:
        if np.all(rel[k:] <= 1e-6):
            need = int(k) + 1; break
    return dict(parts=np.array(parts), terms=np.array(terms), absterms=np.array(absterms), rel=rel,
                need=need, maxterm=float(np.max(np.abs(terms))), final_rel=float(rel[-1]))


def even_power_moment_route(rhoA: np.ndarray, rhoB: np.ndarray, smax: int, use_mp: bool
                            ) -> Tuple[np.ndarray, np.ndarray]:
    """W_{2s} by the SEPARABLE polynomial route (the grid analog of the H2 engine's even r12^{2s}):
    r12^2 = r1^2 + r2^2 - 2 z1 z2 - 2 rho1 rho2 cos(dphi); multinomial expansion; one-electron moments
    mu[p,q] = INT rho rhocyl^p z^q.  Returns (W_2s, sum|terms|) for s=0..smax."""
    from math import comb, factorial
    import mpmath as mp
    mp.mp.dps = 40
    P = 2 * smax
    rc = P2.rho_cyl_f; z = P2.zc_f; w = P2.GI_W
    if use_mp:
        wA = [mp.mpf(float(x)) * mp.mpf(float(y)) for x, y in zip(rhoA * 0 + w, rhoA)]
        wB = [mp.mpf(float(x)) * mp.mpf(float(y)) for x, y in zip(rhoB * 0 + w, rhoB)]
        rcm = [mp.mpf(float(x)) for x in rc]; zm = [mp.mpf(float(x)) for x in z]
        muA = [[mp.mpf(0)] * (P + 1) for _ in range(P + 1)]
        muB = [[mp.mpf(0)] * (P + 1) for _ in range(P + 1)]
        for i in range(len(rcm)):
            rp = mp.mpf(1)
            for p in range(P + 1):
                zq = mp.mpf(1)
                for q in range(P + 1 - p):             # need p+q <= 2s only
                    muA[p][q] += wA[i] * rp * zq; muB[p][q] += wB[i] * rp * zq
                    zq *= zm[i]
                rp *= rcm[i]
    else:
        muA = np.zeros((P + 1, P + 1)); muB = np.zeros((P + 1, P + 1))
        rp = np.ones_like(rc)
        for p in range(P + 1):
            zq = np.ones_like(z)
            for q in range(P + 1 - p):
                muA[p, q] = np.sum(w * rhoA * rp * zq); muB[p, q] = np.sum(w * rhoB * rp * zq)
                zq = zq * z
            rp = rp * rc

    def M1(mu, i, k, l):                               # INT rho (rc^2+z^2)^i z^k rc^l
        acc = 0
        for a_ in range(i + 1):
            acc = acc + comb(i, a_) * mu[2 * a_ + l][2 * (i - a_) + k]
        return acc
    Ws, Aa = [], []
    for s in range(smax + 1):
        tot = mp.mpf(0) if use_mp else 0.0; ab = 0.0
        for i in range(s + 1):
            for j in range(s + 1 - i):
                for k in range(s + 1 - i - j):
                    l = s - i - j - k
                    if l % 2:
                        continue
                    cl = comb(l, l // 2) / 2 ** l
                    coef = factorial(s) // (factorial(i) * factorial(j) * factorial(k) * factorial(l))
                    coef = coef * (-2) ** (k + l) * cl
                    t = coef * M1(muA, i, k, l) * M1(muB, j, k, l)
                    tot = tot + t; ab += abs(float(t))
        Ws.append(float(tot)); Aa.append(ab)
    return np.array(Ws), np.array(Aa)


def odd_power_neumann_route(rhoA: np.ndarray, rhoB: np.ndarray, n: int) -> Tuple[float, float, int]:
    """W_n (n odd) = INT INT rhoA rhoB r12^{n+1}/r12 : multinomial of r12^{n+1} (separable monomials x
    cos^l dphi) -> cos^l = sum_m beta_{l,m} cos(m dphi) -> mode-m prolate-Neumann Coulomb E_m (exact,
    the grid analog of the H2 engine's A K0 - B K1).  Returns (W_n, sum|terms|, n_coulomb_solves)."""
    from math import comb, factorial
    s1 = (n + 1) // 2
    rc = P2.rho_cyl_f; z = P2.zc_f
    GW2 = M._GW2_flat
    tot = 0.0; ab = 0.0; nsolve = 0
    for i in range(s1 + 1):
        for j in range(s1 + 1 - i):
            for k in range(s1 + 1 - i - j):
                l = s1 - i - j - k
                coef = factorial(s1) // (factorial(i) * factorial(j) * factorial(k) * factorial(l))
                coef = coef * (-2.0) ** (k + l)
                fA = rhoA * (rc ** 2 + z ** 2) ** i * z ** k * rc ** l
                fB = rhoB * (rc ** 2 + z ** 2) ** j * z ** k * rc ** l
                for m in range(l % 2, l + 1, 2):
                    beta = comb(l, (l - m) // 2) / 2 ** l * (2 if m > 0 else 1)
                    pot = M.coul_mode_stack_any(fB[None, :], m)[0]
                    Em = float(np.sum(GW2 * fA * pot))
                    nsolve += 1
                    t = coef * beta * Em
                    tot += t; ab += abs(t)
    return tot, ab, nsolve


def closed_1s_control(zeta: float, gam: float, nmax: int, dps: int = 30) -> Dict[str, object]:
    """Exact (mpmath) W_n = <r12^n> and W_lin = <r12 e^{-g r12}> for two same-centre 1s(zeta) densities.
    Angular average: <r12^n>_Omega = [(r1+r2)^{n+2} - |r1-r2|^{n+2}] / (2 r1 r2 (n+2))."""
    import mpmath as mp
    mp.mp.dps = dps
    z = mp.mpf(zeta); g = mp.mpf(gam)
    rad = lambda r: 4 * z ** 3 * r ** 2 * mp.e ** (-2 * z * r)     # 4 pi r^2 |1s|^2

    def avg_pow(r1, r2, n):
        return ((r1 + r2) ** (n + 2) - abs(r1 - r2) ** (n + 2)) / (2 * r1 * r2 * (n + 2))

    def F(s):                                             # antiderivative of s^2 e^{-g s}
        return -mp.e ** (-g * s) * (s ** 2 / g + 2 * s / g ** 2 + 2 / g ** 3)

    def avg_lin(r1, r2):                                  # <r12 e^{-g r12}>_Omega = INT s^2 e^{-gs} ds/(2 r1 r2)
        return (F(r1 + r2) - F(abs(r1 - r2))) / (2 * r1 * r2)

    # r2 = r1 t (t in [0,1], x2 for symmetry): the r1 integral of the power moment is a Gamma function,
    #   W_n = 2 INT_0^1 8 z^6 t [(1+t)^{n+2}-(1-t)^{n+2}] (n+5)! / ((n+2) (2z(1+t))^{n+6}) dt
    def Wn_exact(n):
        return 2 * mp.quad(lambda t: 8 * z ** 6 * t * ((1 + t) ** (n + 2) - (1 - t) ** (n + 2))
                           * mp.factorial(n + 5) / ((n + 2) * (2 * z * (1 + t)) ** (n + 6)), [0, 1])
    Wn = [Wn_exact(n) for n in range(nmax + 1)]
    Wlin = 2 * mp.quad(lambda t, r1: r1 * rad(r1) * rad(r1 * t)
                       * avg_lin(r1, r1 * t), [0, 1], [0, mp.inf])
    return dict(Wn=np.array([float(x) for x in Wn]), Wlin=float(Wlin), Wn_mp=Wn, Wlin_mp=Wlin)


def run_d2() -> None:
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    import mpmath as mp
    print("=" * 100)
    print("D2 -- one-kernel closed-form route: gamma-derivative identities + power series in r12")
    print("=" * 100, flush=True)
    g = KG.GAM
    # ---- (i) derivative identities on the grid kernels (row subset) ----
    rng = np.random.default_rng(7)
    rows = np.sort(rng.choice(P2.NG, 150, replace=False))
    h = 1e-3
    S = lambda gg: (lambda d: (1.0 - np.exp(-gg * d)) / np.maximum(d, 1e-300) * (d > 0)
                    + gg * (d == 0))                      # smooth Yukawa complement, limit g at d=0
    Sp, S0, Sm = (_kernel_rows(S(g + h), rows), _kernel_rows(S(g), rows), _kernel_rows(S(g - h), rows))
    Kexp = _kernel_rows(lambda d: np.exp(-g * d), rows)
    Klin = _kernel_rows(lambda d: d * np.exp(-g * d), rows)
    # production kernel for reference: GeminalM linexp Kf rows vs Klin (same quadrature)
    gemK = P2.build_kernel_gen(lambda d: d * np.exp(-g * d))[rows]
    d1 = (Sp - Sm) / (2 * h); d2 = (Sp - 2 * S0 + Sm) / h ** 2
    e1 = np.max(np.abs(d1 - Kexp)) / np.max(np.abs(Kexp))
    e2 = np.max(np.abs(-d2 - Klin)) / np.max(np.abs(Klin))
    e3 = np.max(np.abs(gemK - Klin)) / np.max(np.abs(Klin))
    print(f"  (i) identities on grid kernels (150 rows, FD h={h}): "
          f"dS/dg == K[e^-gr] rel {e1:.1e};  -d2S/dg2 == K[r e^-gr] rel {e2:.1e}  (O(h^2) FD floor);  "
          f"production linexp Kf == K[r e^-gr] rel {e3:.1e}")
    print("      (Y_g = Coulomb - S_g, so e^-gr = -dY/dg and r e^-gr = d2Y/dg2 hold; the derivative acts on the")
    print("       SMOOTH part S_g, which is itself an NG x NG grid kernel -- no Yukawa algebra is gained.)",
          flush=True)
    del Sp, S0, Sm, Kexp, Klin, gemK, d1, d2

    # ---- (ii) densities: real signed pi-reference pair + isotropic 1s control ----
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=1e-2,
                                        verbose=False)
    nval = info['nval']; m_no = info['m_no']
    rhoA = (nval[0] * nval[1]).real; rhoB = (nval[1] * nval[1]).real
    print(f"  (ii) A = rho_(NO0,NO1) (m={m_no[1]-m_no[0]}, INT={P2.gint(rhoA):+.2e}, signed), "
          f"B = rho_(NO1,NO1) (INT={P2.gint(rhoB):.6f});  control: |1s(zeta=1)|^2 on H x itself", flush=True)
    zeta = 1.0
    rB = P2.rB_f
    rho1s = (zeta ** 3 / np.pi) * np.exp(-2 * zeta * rB)
    NMAX = 60
    t1 = time.time()
    gm = grid_power_moments(np.stack([rhoA, rho1s]), np.stack([rhoB, rho1s]), NMAX, g)
    print(f"      grid W_n (n=0..{NMAX}) built in {time.time()-t1:.0f}s", flush=True)
    t1 = time.time()
    cc = closed_1s_control(zeta, g, 40, dps=20)
    print(f"      closed-form 1s control (mpmath) built in {time.time()-t1:.0f}s", flush=True)
    print(f"      grid vs closed 1s: W_lin grid {gm['Wlin'][1]:.12f}  closed {cc['Wlin']:.12f}  rel "
          f"{abs(gm['Wlin'][1]-cc['Wlin'])/cc['Wlin']:.1e};  W_10 rel "
          f"{abs(gm['Wn'][10,1]-cc['Wn'][10])/cc['Wn'][10]:.1e};  W_40 rel "
          f"{abs(gm['Wn'][40,1]-cc['Wn'][40])/cc['Wn'][40]:.1e}")

    res = {}
    for k, nm in enumerate(['signed A.B', '1s control']):
        for geminal, shift, refval, bound in (('linexp r e^-gr', 1, gm['Wlin'][k], gm['Blin'][k]),
                                             ('exp e^-gr', 0, gm['Wexp'][k], gm['Bexp'][k])):
            st = series_table(gm['Wn'][:, k], gm['Wabs'][:, k], refval, g, shift)
            stm = series_table(gm['Wn'][:, k], gm['Wabs'][:, k], refval, g, shift, use_mp=True)
            fl = 2.2e-16 * bound / abs(refval)
            res[(nm, geminal)] = dict(st=st, stm=stm, ref=refval, floor=fl)
            print(f"\n  series  {nm:11s} {geminal:15s}: grid W = {refval:+.12e}; float64 roundoff floor "
                  f"eps*INT|rhoA rhoB| e^(+g r)/|W| = {fl:.1e};  terms for 1e-6 rel: float64 {st['need']}, "
                  f"mpmath {stm['need']};  max|term|/|W| = {st['maxterm']/abs(refval):.2e}; "
                  f"final rel (N={len(st['parts'])}) {st['final_rel']:.1e} / {stm['final_rel']:.1e}")
            print(f"    {'N':>3s} {'partial sum':>20s} {'rel err':>9s} {'term':>11s} {'|term| bound':>12s}")
            for N in (0, 1, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, len(st['parts']) - 1):
                if N < len(st['parts']):
                    print(f"    {N+1:3d} {st['parts'][N]:+20.12e} {st['rel'][N]:9.1e} {st['terms'][N]:+11.3e} "
                          f"{st['absterms'][N]:12.3e}")
    # control: series built from the EXACT closed-form moments (no grid) -> intrinsic convergence
    stc = series_table(np.array(cc['Wn']), np.array(cc['Wn']), cc['Wlin'], g, 1)
    print(f"\n  1s control, series of EXACT closed-form moments vs exact W_lin: terms for 1e-6 = {stc['need']}, "
          f"final rel (40 terms) {stc['final_rel']:.1e}", flush=True)

    # ---- (iii) exact-power routes (the H2-engine analog on grid densities) ----
    print("\n  (iii) EXACT algebraic routes for r12^n on grid densities (grid analog of the H2 r12 engine):")
    smax = 15
    t1 = time.time(); Wev64, Aev = even_power_moment_route(rhoA, rhoB, smax, use_mp=False); te64 = time.time() - t1
    t1 = time.time(); Wevmp, _ = even_power_moment_route(rhoA, rhoB, smax, use_mp=True); temp = time.time() - t1
    print(f"    EVEN r12^2s via separable one-electron moments (float64 {te64:.1f}s, mpmath-40 {temp:.1f}s):")
    print(f"    {'n':>3s} {'grid kernel W_n':>18s} {'moment f64':>18s} {'rel f64':>8s} {'moment mp40':>18s} "
          f"{'rel mp':>8s} {'sum|t|/|W|':>10s}")
    ev_rows = []
    for s in range(0, smax + 1):
        n = 2 * s; wg = gm['Wn'][n, 0]
        r64 = abs(Wev64[s] - wg) / abs(wg); rmp = abs(Wevmp[s] - wg) / abs(wg)
        ev_rows.append((n, wg, Wev64[s], r64, Wevmp[s], rmp, Aev[s] / abs(wg)))
        if s in (0, 1, 2, 3, 5, 7, 10, 12, 15):
            print(f"    {n:3d} {wg:+18.10e} {Wev64[s]:+18.10e} {r64:8.1e} {Wevmp[s]:+18.10e} {rmp:8.1e} "
                  f"{Aev[s]/abs(wg):10.1e}")
    print("    ODD r12^n via multinomial x mode-m prolate-Neumann Coulomb (float64):")
    od_rows = []
    for n in (1, 3, 5, 7, 9, 11):
        t1 = time.time()
        Wo, Ao, ns = odd_power_neumann_route(rhoA, rhoB, n)
        wg = gm['Wn'][n, 0]
        od_rows.append((n, wg, Wo, abs(Wo - wg) / abs(wg), Ao / abs(wg), ns, time.time() - t1))
        print(f"    n={n:2d}: grid {wg:+.10e}  Neumann {Wo:+.10e}  rel {abs(Wo-wg)/abs(wg):.1e}  "
              f"sum|t|/|W| {Ao/abs(wg):.1e}  ({ns} mode-Coulomb solves, {time.time()-t1:.1f}s)", flush=True)
    np.savez(PART("d2"), ident=np.array([e1, e2, e3]), NMAX=NMAX, Wn=gm['Wn'], Wabs=gm['Wabs'],
             Wlin=gm['Wlin'], Wexp=gm['Wexp'], Blin=gm['Blin'], Bexp=gm['Bexp'],
             closed_Wn=cc['Wn'], closed_Wlin=cc['Wlin'],
             series=np.array({f"{a}|{b}": dict(need64=v['st']['need'], needmp=v['stm']['need'],
                                                rel=v['st']['rel'], relmp=v['stm']['rel'],
                                                terms=v['st']['terms'], maxterm=v['st']['maxterm'],
                                                floor=v['floor'], ref=v['ref'])
                              for (a, b), v in res.items()}, dtype=object),
             even_rows=np.array(ev_rows), odd_rows=np.array(od_rows), t_even64=te64, t_evenmp=temp)
    print(f"\n  saved {PART('d2')}  {_t(T0)}", flush=True)


# =========================================================================== #
# nnz probe: what an O(nnz) compound field would cost vs the dense nact^3 x NG
# =========================================================================== #
def run_nnz() -> None:
    T0 = time.time()
    print("=" * 100)
    print("NNZ probe -- distinct index combos in T vs dense nact^k (O(nnz) field-memory estimate)")
    print("=" * 100, flush=True)
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    out = []
    print(f"  {'trunc':>7s} {'ndet':>5s} {'nact':>5s} {'nnz':>8s} {'max uniq3':>10s} {'nact^3':>9s} "
          f"{'dense GB':>9s} {'O(nnz) GB':>9s} {'E_trunc':>10s}")
    for tr in (1e-2, 5e-3, 3e-3, 2e-3, 1e-3):
        refm, info = M.build_ref_from_prims(ref, "probe", cfg_tag='pi', trunc=tr, verbose=False)
        idx = [refm.T_iA, refm.T_iB, refm.T_iC, refm.T_iD]
        u3 = max(np.unique(np.stack([idx[a], idx[b], idx[c]], axis=1), axis=0).shape[0]
                 for a, b, c in combinations(range(4), 3))
        n = refm.nact; nnz = len(refm.T_val)
        dense = 1.4 * n ** 3 * P2.NG * 8 / 1e9
        onnz = 3.0 * u3 * P2.NG * 8 / 1e9                 # field + dressed copy + product temp
        out.append((tr, len(info['kept']), n, nnz, u3, dense, onnz, info['E_trunc']))
        print(f"  {tr:7.0e} {len(info['kept']):5d} {n:5d} {nnz:8d} {u3:10d} {n**3:9d} {dense:9.1f} "
              f"{onnz:9.2f} {info['E_trunc']:10.6f}  {_t(T0)}", flush=True)
    np.savez(PART("nnz"), rows=np.array(out))


def run_merge() -> None:
    out = {}
    for nm in ("d0pi", "d0sigma", "d2", "nnz"):
        p = PART(nm)
        if os.path.exists(p):
            z = np.load(p, allow_pickle=True)
            for k in z.files:
                out[f"{nm}__{k}"] = z[k]
    np.savez(os.path.join(DATA, "lih_gvee_diag.npz"), **out)
    print(f"merged {len(out)} arrays -> debug/data/lih_gvee_diag.npz")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['d0', 'd2', 'nnz', 'merge'])
    ap.add_argument('--ref', default='pi', choices=['pi', 'sigma'])
    a = ap.parse_args()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    {'d0': lambda: run_d0(a.ref), 'd2': run_d2, 'nnz': run_nnz, 'merge': run_merge}[a.cmd]()
