r"""O(nnz) rewrite of g_Vee for the LiH many-determinant explicit-r12 "marriage" engine.

Plan: debug/lih_pi_m16_build_plan.md (STATUS 2026-09-25); diagnostic: debug/lih_gvee_diag.py.

WHAT IS SPARSIFIED
------------------
The dense reducer `lih_pi_m16.vee_reduce_m` materialises every compound field (a kept electron dressed
by an f-chain / f-star, e.g. f12 f13 / r34) over ALL nact^k index combos x NG grid points (~nact^3 NG,
the TREE4 memory wall), and builds the full (nact x nact) triangle matrix per middle vertex (mode-m
Neumann solves over every rank x nact row).  But every piece is finally contracted against T's
nonzeros, so only the index tuples that OCCUR in T are ever read.  This module:

 1. replays the dense reduction order SYMBOLICALLY (same vertex order, same dressing keys, same
    integration masks) into a small expression tree of field nodes;
 2. evaluates each field / integrated factor ONLY at the unique projections of T's nonzero index
    tuples onto that field's axes -- row-wise, because every dressing is a row-wise linear map (the
    dressing mode is |total circulation of the row|, a per-row property);
 3. streams the kept-pair Coulomb/Yukawa over R-tuple chunks bounded by a memory cap, dotting only
    the (L,R) tuple pairs present in T (no ndL x ndR Emat);
 4. defers every triangle term, then evaluates the triangle ONLY at the (cm, L, R) combos present in
    T -- the mode-m Neumann solves run on the needed R rows only.  The per-middle-vertex randomized
    eigh / SVD calls are replayed in EXACTLY the dense order (first-encounter cm order, every combo),
    so the random streams are identical and sparse == dense to roundoff, not to ~1e-7.

Compound dressed rows are memoised across the 216 triples by (structure signature, tuple code) in a
memory-bounded cache.  Nothing in lih_pi_m16 / lih_marriage_phase2 / geovac is modified -- the dense
path stays the oracle.

API:  all_pieces_sparse(gem, refm, T0, mem_cap_gb=6.0)  -> same dict as lih_pi_m16.all_pieces_m
      hg_Vee_sparse(gem, refm, Fbar, mem_cap_gb=6.0)    -> (h_Vee, g_Vee, Vee, FVee, F2Vee)
Gates:  python debug/lih_gvee_sparse.py r1   (pi trunc 1e-2 vs banked dense)
        python debug/lih_gvee_sparse.py r2   (sigma 74-det vs banked -8.037236)
"""
from __future__ import annotations

import os
import sys
import time
from typing import Dict, List, Optional, Tuple

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
import geovac.lih_r12ci.triangle as TRI
import lih_marriage_phase2 as P2
import lih_pi_m16 as M

NG = P2.NG
GI_W = P2.GI_W
_GW2 = M._GW2_flat
PAIRS6 = P2.PAIRS6

# plant hooks for the fire test (all False in production)
_DROP_NESTED = False          # drop kept-pair terms whose field is a nested (depth>=2) chain (TREE4)


# =========================================================================== #
# memory helpers
# =========================================================================== #
def peak_rss_gb() -> float:
    try:
        import psutil
        mi = psutil.Process().memory_info()
        return float(getattr(mi, 'peak_wset', mi.rss)) / 1e9
    except Exception:
        return float('nan')


def rss_gb() -> float:
    try:
        import psutil
        return float(psutil.Process().memory_info().rss) / 1e9
    except Exception:
        return float('nan')


# =========================================================================== #
# symbolic field nodes
# =========================================================================== #
class FNode:
    """field of vertex v: rho_v(r) * prod_children Dress_keys(child field)(r).
    axes = (v,) + concatenated child axes (the dense cf_mul order)."""

    def __init__(self, v: int, children: Optional[List[Tuple["FNode", object]]] = None):
        self.v = v
        self.children = list(children or [])
        ax = [v]
        for ch, _k in self.children:
            ax.extend(ch.axes)
        self.axes: Tuple[int, ...] = tuple(ax)
        self.sig = "r" + "".join("(" + M._kind_of(k) + ":" + ch.sig + ")" for ch, k in self.children)

    def depth(self) -> int:
        return 0 if not self.children else 1 + max(ch.depth() for ch, _ in self.children)


class Ctx:
    """per-(gem, refm) evaluation context: memo of dressed rows + memory cap."""

    def __init__(self, gem: "M.GeminalM", refm: "M.RefM", mem_cap_gb: float):
        self.gem = gem
        self.refm = refm
        self.nact = refm.nact
        self.cap = float(mem_cap_gb) * 1e9
        self.chunk_rows = max(64, int(self.cap / (NG * 8 * 10)))
        self.cache_budget = self.cap / 3.0
        self.full: Dict[str, np.ndarray] = {}                  # kind -> dressed rho (nact, NG)
        self.dc: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}  # sig -> (sorted codes, rows)
        self.stats = dict(dress_rows=0, pot_rows=0, tri_solve_rows=0, tri_cm=0)

    def cache_bytes(self) -> int:
        return sum(r.nbytes for (_c, r) in self.dc.values())

    # ---- tuple helpers ---- #
    def codes(self, tup: np.ndarray) -> np.ndarray:
        k = tup.shape[1]
        if k == 1:
            return tup[:, 0].astype(np.int64)
        return np.ravel_multi_index(tuple(tup.T.astype(np.int64)), (self.nact,) * k).astype(np.int64)

    def uncode(self, codes: np.ndarray, k: int) -> np.ndarray:
        if k == 1:
            return codes[:, None].astype(np.int64)
        return np.stack(np.unravel_index(codes, (self.nact,) * k), axis=1).astype(np.int64)

    def circ(self, tup: np.ndarray) -> np.ndarray:
        return self.refm.mA[tup].sum(axis=1).astype(np.int64)

    # ---- field evaluation ---- #
    def eval_field(self, node: FNode, tup: np.ndarray) -> np.ndarray:
        """rows (n, NG) of the field of `node` at index tuples tup (n, len(node.axes))."""
        out = self.refm.rho[tup[:, 0]].copy()
        col = 1
        for ch, keys in node.children:
            kc = len(ch.axes)
            sub = tup[:, col:col + kc]
            col += kc
            out *= self.dressed(ch, keys, sub)
        return out

    def dressed(self, ch: FNode, keys, sub: np.ndarray) -> np.ndarray:
        kind = M._kind_of(keys)
        if not ch.children:                                   # single-axis: full nact table, cached
            if kind not in self.full:
                self.full[kind] = M.dress_rows_m(self.gem, keys, self.refm.rho, self.refm.mA)
                self.stats['dress_rows'] += self.nact
            return self.full[kind][sub[:, 0]]
        sig = kind + ":" + ch.sig

        def compute(mt: np.ndarray) -> np.ndarray:
            blocks = []
            for s in range(0, len(mt), self.chunk_rows):
                t = mt[s:s + self.chunk_rows]
                rows = self.eval_field(ch, t)
                blocks.append(M.dress_rows_m(self.gem, keys, rows, self.circ(t)))
                self.stats['dress_rows'] += len(t)
            return np.concatenate(blocks, axis=0) if len(blocks) > 1 else blocks[0]
        return self.memo(sig, sub, compute)

    def pot(self, node: FNode, kernel_kind: str, gam_k, tup: np.ndarray) -> np.ndarray:
        """kept-pair potential rows (mode |circ| of each row) at tuples, memoised by structure."""
        sig = f"pot:{kernel_kind}:{gam_k}:" + node.sig

        def compute(mt: np.ndarray) -> np.ndarray:
            out = np.empty((len(mt), NG))
            for s in range(0, len(mt), self.chunk_rows):
                t = mt[s:s + self.chunk_rows]
                Rrows = self.eval_field(node, t)
                mR = np.abs(self.circ(t))
                for mode in sorted(set(mR.tolist())):
                    rows = mR == mode
                    out[s:s + len(t)][rows] = M._kept_pot(self.gem, Rrows[rows], kernel_kind, gam_k,
                                                          int(mode))
                self.stats['pot_rows'] += len(t)
            return out
        return self.memo(sig, tup, compute)

    def memo(self, sig: str, sub: np.ndarray, compute) -> np.ndarray:
        """rows for index tuples `sub`, from the bounded (sig, tuple-code) cache; misses via compute."""
        k = sub.shape[1]
        cod = self.codes(sub)
        ucod, inv = np.unique(cod, return_inverse=True)
        cached = self.dc.get(sig)
        if cached is not None:
            ccod, crows = cached
            pos = np.searchsorted(ccod, ucod)
            pos_c = np.minimum(pos, len(ccod) - 1)
            have = ccod[pos_c] == ucod
        else:
            ccod = np.zeros(0, np.int64); crows = np.zeros((0, NG)); have = np.zeros(len(ucod), bool)
            pos_c = np.zeros(len(ucod), np.int64)
        miss = ucod[~have]
        newrows = None
        if len(miss):
            newrows = compute(self.uncode(miss, k))
        # assemble the requested unique rows
        urows = np.empty((len(ucod), NG))
        if have.any():
            urows[have] = crows[pos_c[have]]
        if newrows is not None:
            urows[~have] = newrows
        # merge into the cache if the budget allows (else keep only this signature's new rows)
        if newrows is not None:
            add_bytes = newrows.nbytes
            if self.cache_bytes() + add_bytes > self.cache_budget:
                self.dc.clear()
                ccod = np.zeros(0, np.int64); crows = np.zeros((0, NG))
            if newrows.nbytes <= self.cache_budget:
                allc = np.concatenate([ccod, miss]); allr = np.concatenate([crows, newrows], axis=0)
                o = np.argsort(allc, kind='stable')
                self.dc[sig] = (allc[o], allr[o])
        return urows[inv]

    # ---- integrated factor at tuples ---- #
    def integrate(self, node: FNode, tup: np.ndarray) -> np.ndarray:
        vals = np.empty(len(tup))
        for s in range(0, len(tup), self.chunk_rows):
            t = tup[s:s + self.chunk_rows]
            vals[s:s + len(t)] = self.eval_field(node, t) @ GI_W
        return vals * (self.circ(tup) == 0)


# =========================================================================== #
# symbolic replay of vee_reduce_m (identical vertex order / keys)
# =========================================================================== #
def _replay(rem_f, coul):
    e, gg = coul
    kept = {e, gg}
    field = {x: FNode(x) for x in (1, 2, 3, 4)}
    edge_list = [set(ed) for ed in rem_f]
    alive = {1, 2, 3, 4}
    factors: List[FNode] = []
    tri: Optional[FNode] = None

    def dnb(u):
        s = set()
        for ed in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    while True:
        nk = [u for u in alive if u not in kept]
        if not nk:
            break
        u = next((c for c in sorted(nk) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(nk) if len(dnb(c)) == 1), None)
        if u is None:
            u = next((c for c in sorted(nk) if dnb(c) == kept), None)
            if u is None:
                raise RuntimeError(f"vee: no reducible vertex; rem_f={rem_f} coul={coul}")
            if field[u].children:
                raise RuntimeError("triangle vertex carries dressings")
            tri = field[u]
            edge_list = [ed for ed in edge_list if u not in ed]; alive.discard(u)
            continue
        inc = [ed for ed in edge_list if u in ed]
        nb = dnb(u)
        if not nb:
            factors.append(field[u])
        else:
            w = next(iter(nb)); mult = len(inc)
            keys = 'f' if mult == 1 else ('f', 'f')
            field[w] = FNode(field[w].v, field[w].children + [(field[u], keys)])
        alive.discard(u)
        edge_list = [ed for ed in edge_list if u not in ed]
    return factors, field[e], field[gg], tri


class TriCollector:
    def __init__(self) -> None:
        self.recs: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
        self.order: List[int] = []                            # first-encounter cm order (dense rng order)
        self._seen: set = set()

    def add(self, im, il, ir, w) -> None:
        self.recs.append((im, il, ir, w))
        for cm in np.unique(im):                              # dense: `for cm in np.unique(im)`
            c = int(cm)
            if c not in self._seen:
                self._seen.add(c); self.order.append(c)


def vee_reduce_sparse(rem_f, coul, kernel_kind, gam_k, ctx: Ctx, tri_col: TriCollector) -> float:
    """O(nnz) twin of lih_pi_m16.vee_reduce_m.  Triangle terms are deferred into tri_col."""
    refm = ctx.refm; gem = ctx.gem
    factors, Lnode, Rnode, tri = _replay(rem_f, coul)
    idxmap = {1: refm.T_iA, 2: refm.T_iB, 3: refm.T_iC, 4: refm.T_iD}

    def tuples(node: FNode, sel=None) -> np.ndarray:
        cols = [idxmap[a] if sel is None else idxmap[a][sel] for a in node.axes]
        return np.stack(cols, axis=1).astype(np.int64)

    wnnz = refm.T_val.copy()
    for F in factors:
        tup = tuples(F)
        cod = ctx.codes(tup)
        ucod, inv = np.unique(cod, return_inverse=True)
        vals = ctx.integrate(F, ctx.uncode(ucod, tup.shape[1]))
        wnnz = wnnz * vals[inv]
    keepnz = np.abs(wnnz) > 0
    if not np.any(keepnz):
        return 0.0
    wnnz = wnnz[keepnz]

    if tri is not None:
        assert len(Lnode.axes) == 1 and len(Rnode.axes) == 1
        im = idxmap[tri.axes[0]][keepnz]; il = idxmap[Lnode.axes[0]][keepnz]
        ir = idxmap[Rnode.axes[0]][keepnz]
        tri_col.add(im.astype(np.int64), il.astype(np.int64), ir.astype(np.int64), wnnz)
        return 0.0

    if _DROP_NESTED and (Lnode.depth() >= 2 or Rnode.depth() >= 2):
        return 0.0

    TL = tuples(Lnode, keepnz); TR = tuples(Rnode, keepnz)
    sel = (ctx.circ(TL) + ctx.circ(TR)) == 0
    if not sel.any():
        return 0.0
    cL = ctx.codes(TL[sel]); cR = ctx.codes(TR[sel]); ww = wnnz[sel]
    # group identical (L,R) tuple pairs
    kL = TL.shape[1]; kR = TR.shape[1]
    uR, invR = np.unique(cR, return_inverse=True)
    uL, invL = np.unique(cL, return_inverse=True)
    pair = invL.astype(np.int64) * len(uR) + invR
    up, pinv = np.unique(pair, return_inverse=True)
    wp = np.bincount(pinv, weights=ww)
    pl = up // len(uR); pr = up % len(uR)
    tupR = ctx.uncode(uR, kR); tupL = ctx.uncode(uL, kL)
    modeR = np.abs(ctx.circ(tupR))
    total = 0.0
    cr = ctx.chunk_rows
    for s in range(0, len(uR), cr):
        rsl = slice(s, min(s + cr, len(uR)))
        pot = ctx.pot(Rnode, kernel_kind, gam_k, tupR[rsl])
        psel = (pr >= rsl.start) & (pr < rsl.stop)
        if not psel.any():
            continue
        pls = pl[psel]; prs = pr[psel] - rsl.start; wps = wp[psel]
        ul2, linv = np.unique(pls, return_inverse=True)
        Lrows = ctx.eval_field(Lnode, tupL[ul2]) * GI_W[None, :]
        for t0 in range(0, len(pls), cr):
            t1 = min(t0 + cr, len(pls))
            dots = np.einsum('ij,ij->i', Lrows[linv[t0:t1]], pot[prs[t0:t1]])
            total += float((wps[t0:t1] * dots).sum())
        del Lrows, pot
    return float(total)


# =========================================================================== #
# triangle: evaluate only the needed (cm, L, R) entries, replaying the dense rng order
# =========================================================================== #
def tri_entries(ctx: Ctx, cm: int, il: np.ndarray, ir: np.ndarray) -> np.ndarray:
    """Triangle matrix entries M_cm[il, ir] (== lih_pi_m16.get_tri_matrix_m(...)[il, ir]).
    Replays every Kmu_modes / _rand_svd_middle call of the dense builder IN ORDER (identical random
    streams); only the per-(muL) Coulomb solves are restricted to the needed rows."""
    gem, refm = ctx.gem, ctx.refm
    mu = refm.rho[cm]; mu_mid = int(refm.mA[cm])
    rho = refm.rho; mArows = refm.mA
    gmu = P2.GEO * mu
    key = f"{refm.tag}_trim_mu{int(cm)}_{int(KG.USE_EXACT_NEUMANN)}"
    out = np.zeros(len(il))
    combos = []
    MMAX = M.MMAX
    for A in range(MMAX + 1):
        for B in range(MMAX + 1):
            for sa in ((0,) if A == 0 else (+1, -1)):
                for sb in ((0,) if B == 0 else (+1, -1)):
                    if sa * A + sb * B == mu_mid:
                        combos.append((A, sa, B, sb))

    wcomb = M.tri_wcomb                                       # shared weight (A>=1,B=0 fix 2026-09-26)

    mL = mArows[il]; mR = mArows[ir]
    used_keys = []
    for (A, sa, B, sb) in combos:
        if A == B:
            kk = f"{key}_A{A}"
            lam, V = M.Kmu_modes(gem.Fm[A], mu, A, kk)
            used_keys.append((kk, A, 30))
            d_A = 2 * np.pi if A == 0 else np.pi
            sig = lam / d_A; U = V; W = V
        else:
            sig, U, W = M._rand_svd_middle(gem.Fm[A], gem.Fm[B], gmu, rank=30)
        if sig.size == 0:
            continue
        w = wcomb(A, sa, B, sb)
        rk = sig.size
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
            # batch the Coulomb solves over (rank x needed-R) rows, bounded by the chunk size
            rstep = max(1, ctx.chunk_rows // max(nR, 1))
            for r0 in range(0, rk, rstep):
                r1 = min(r0 + rstep, rk)
                Bst = (rhoR[None, :, :] * W[:, r0:r1].T[:, None, :]).reshape((r1 - r0) * nR, NG)
                Pr = M.coul_mode_stack_any(Bst, int(c)).reshape(r1 - r0, nR, NG)
                ctx.stats['tri_solve_rows'] += (r1 - r0) * nR
                for r in range(r0, r1):
                    Lr = rhoL * (U[:, r] * _GW2)[None, :]
                    blk += w * sig[r] * (Lr @ Pr[r - r0].T)
                del Pr, Bst
            out[np.nonzero(psel)[0]] += blk[linv, rinv]
    for ck in used_keys:                                      # free this cm's mode cache (memory)
        TRI._MODE_CACHE.pop(ck, None)
    ctx.stats['tri_cm'] += 1
    return out


def eval_triangles(ctx: Ctx, col: TriCollector) -> float:
    if not col.recs:
        return 0.0
    im = np.concatenate([r[0] for r in col.recs]); il = np.concatenate([r[1] for r in col.recs])
    ir = np.concatenate([r[2] for r in col.recs]); w = np.concatenate([r[3] for r in col.recs])
    n = ctx.nact
    code = (im * n + il) * n + ir
    uc, inv = np.unique(code, return_inverse=True)
    wsum = np.bincount(inv, weights=w)
    ucm = uc // (n * n); uil = (uc // n) % n; uir = uc % n
    total = 0.0
    for cm in col.order:
        s = ucm == cm
        vals = tri_entries(ctx, cm, uil[s], uir[s])
        total += float((wsum[s] * vals).sum())
    return total


# =========================================================================== #
# h_Vee / g_Vee and the full pieces dict
# =========================================================================== #
def _triple_args(p, q, r, gem):
    rset = frozenset(r); fedges = [frozenset(p), frozenset(q)]
    n_on_r = sum(1 for ed in fedges if ed == rset)
    rem_f = [tuple(ed) for ed in fedges if ed != rset]
    kk = 'coul' if n_on_r == 0 else 'yuk'
    gam_k = None if n_on_r == 0 else (gem.g if n_on_r == 1 else 2 * gem.g)
    return rem_f, kk, gam_k


def hg_Vee_sparse(gem: "M.GeminalM", refm: "M.RefM", Fbar: float, mem_cap_gb: float = 6.0,
                  verbose: bool = False, T0: Optional[float] = None):
    ctx = Ctx(gem, refm, mem_cap_gb)
    col = TriCollector()
    t0 = time.time() if T0 is None else T0
    Vee = sum(vee_reduce_sparse([], r, 'coul', None, ctx, col) for r in PAIRS6)
    FVee = 0.0
    for p in PAIRS6:
        for r in PAIRS6:
            if frozenset(p) == frozenset(r):
                FVee += vee_reduce_sparse([], r, 'yuk', gem.g, ctx, col)
            else:
                FVee += vee_reduce_sparse([tuple(p)], r, 'coul', None, ctx, col)
    if col.recs:
        raise RuntimeError("unexpected triangle in <V_ee>/<F V_ee>")
    F2_nt = 0.0
    for p in PAIRS6:
        for q in PAIRS6:
            for r in PAIRS6:
                rem_f, kk, gk = _triple_args(p, q, r, gem)
                F2_nt += vee_reduce_sparse(rem_f, r, kk, gk, ctx, col)
    if verbose:
        print(f"      [sparse g_Vee] non-triangle done: F2(non-tri)={F2_nt:+.10f}  "
              f"tri cm={len(col.order)}  cache={ctx.cache_bytes()/1e9:.2f} GB  "
              f"rss={rss_gb():.2f} GB  [{time.time()-t0:7.1f}s]", flush=True)
    F2_tri = eval_triangles(ctx, col)
    F2Vee = F2_nt + F2_tri
    if verbose:
        print(f"      [sparse g_Vee] triangle F2={F2_tri:+.10f}  stats={ctx.stats}  "
              f"[{time.time()-t0:7.1f}s]", flush=True)
    h_Vee = FVee - Fbar * Vee
    g_Vee = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    return h_Vee, g_Vee, Vee, FVee, F2Vee


def all_pieces_sparse(gem: "M.GeminalM", refm: "M.RefM", T0: Optional[float] = None,
                      mem_cap_gb: float = 6.0, verbose: bool = True) -> dict:
    """Same pieces dict as lih_pi_m16.all_pieces_m (cheap pieces via the identical code; g_Vee O(nnz))."""
    T0 = time.time() if T0 is None else T0
    _t = M._t
    p = {}
    Fbar = M.fbar(gem, refm); p['Fbar'] = Fbar
    sig2 = M.sigma2(gem, refm, Fbar); p['sigma2'] = sig2
    if verbose:
        print(f"    Fbar={Fbar:.6f}  sigma2={sig2:.6f}  {_t(T0)}", flush=True)
    h_Vne, g_Vne, Vne_exp = M.hg_Vne(gem, refm, Fbar)
    p['h_Vne'] = h_Vne; p['g_Vne'] = g_Vne; p['Vne_exp'] = Vne_exp
    if verbose:
        print(f"    h_Vne={h_Vne:+.6f}  g_Vne={g_Vne:+.6f}  <V_ne>={Vne_exp:+.6f}  {_t(T0)}", flush=True)
    h_T, g_T, tinfo = M.hg_T(gem, refm, Fbar, sig2)
    p['h_T'] = h_T; p['g_T'] = g_T; p.update({f"T_{k}": v for k, v in tinfo.items()})
    if verbose:
        print(f"    h_T={h_T:+.6f}  g_T={g_T:+.6f}  {_t(T0)}", flush=True)
    h_Vee, g_Vee, Vee, FVee, F2Vee = hg_Vee_sparse(gem, refm, Fbar, mem_cap_gb, verbose, T0)
    p['h_Vee'] = h_Vee; p['g_Vee'] = g_Vee; p['Vee'] = Vee; p['FVee'] = FVee; p['F2Vee'] = F2Vee
    if verbose:
        print(f"    h_Vee={h_Vee:+.6f}  g_Vee={g_Vee:+.6f}  <V_ee>={Vee:+.6f}  {_t(T0)}", flush=True)
    T_kin = 0.5 * p['T_S2']
    p['T_kin'] = T_kin
    p['E0_grid'] = T_kin + Vne_exp + Vee + M.V_NN
    if verbose:
        print(f"    <T>={T_kin:+.6f}  E0_grid={p['E0_grid']:+.6f}  {_t(T0)}", flush=True)
    return p


# =========================================================================== #
# gates R1 / R2
# =========================================================================== #
_PIECES = ['Fbar', 'sigma2', 'h_Vne', 'g_Vne', 'h_T', 'g_T', 'Vee', 'FVee', 'F2Vee', 'h_Vee', 'g_Vee',
           'E0_grid']


def _table(p: dict, ref: dict) -> float:
    worst = 0.0
    print(f"  {'piece':8s} {'dense (banked)':>22s} {'sparse':>22s} {'diff':>10s}")
    for k in _PIECES:
        if k not in ref:
            continue
        d = p[k] - ref[k]
        if k in ('Vee', 'FVee', 'F2Vee', 'h_Vee', 'g_Vee'):
            worst = max(worst, abs(d))
        print(f"  {k:8s} {ref[k]:+22.14f} {p[k]:+22.14f} {d:+10.2e}")
    return worst


def run_r1(mem_cap_gb: float) -> bool:
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    print("=" * 96); print("GATE R1 -- pi trunc 1e-2: sparse g_Vee vs banked dense"); print("=" * 96, flush=True)
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=1e-2,
                                        verbose=False)
    print(f"  nact={refm.nact} nnz={len(refm.T_val)} E_trunc={info['E_trunc']:.6f}  {M._t(T0)}", flush=True)
    z = np.load(os.path.join(DATA, "lih_pi_m16_energy_trunc1e-2.npz"))
    bank = {k[len('linexp_'):]: float(z[k]) for k in z.files if k.startswith('linexp_') and z[k].ndim == 0}
    gem = M.GeminalM('linexp')
    tg = time.time()
    p = all_pieces_sparse(gem, refm, T0, mem_cap_gb=mem_cap_gb)
    E, c, _h = M.assemble_c(p, p['E0_grid'])
    worst = _table(p, bank)
    dE = E - bank['E_R12']
    print(f"\n  E_R12 sparse = {E:.9f}   banked dense = {bank['E_R12']:.9f}   dE = {dE:+.2e} Ha   c={c:+.6f}")
    print(f"  worst Vee-piece diff = {worst:.2e}   wall(pieces) = {time.time()-tg:.1f} s   "
          f"peak RSS = {peak_rss_gb():.2f} GB", flush=True)
    ok = abs(dE) < 1e-9
    print(f"  GATE R1: {'PASS' if ok else ('PASS(<1e-6)' if abs(dE) < 1e-6 else 'FAIL')}", flush=True)
    return abs(dE) < 1e-6


def run_r2(mem_cap_gb: float) -> bool:
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    print("=" * 96); print("GATE R2 -- sigma 74-det: sparse g_Vee vs banked -8.037236"); print("=" * 96, flush=True)
    p2ref = P2.build_ref_artifact()
    refm = M._wrap_p2_ref(p2ref)
    z = np.load(os.path.join(DATA, "lih_marriage_phase3_m1.npz"), allow_pickle=True)
    bank = dict(z['m1'].item()['m1']['linexp']['pieces'])
    gem = M.GeminalM('linexp')
    tg = time.time()
    p = all_pieces_sparse(gem, refm, T0, mem_cap_gb=mem_cap_gb)
    E, c, _h = M.assemble_c(p, p['E0_grid'])
    worst = _table(p, bank)
    dE = E - (-8.037236)
    print(f"\n  E_R12 sparse = {E:.9f}   banked = -8.037236   dE = {dE:+.2e} Ha (banked is 6 dp)")
    Eb, _, _ = M.assemble_c(bank, bank['E0_grid'])
    print(f"  E_R12 from banked pieces = {Eb:.9f}  sparse - that = {E - Eb:+.2e} Ha")
    print(f"  worst Vee-piece diff = {worst:.2e}   wall(pieces) = {time.time()-tg:.1f} s   "
          f"peak RSS = {peak_rss_gb():.2f} GB", flush=True)
    ok = abs(E - Eb) < 1e-6
    print(f"  GATE R2: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('gate', choices=['r1', 'r2'])
    ap.add_argument('--mem-cap-gb', type=float, default=6.0)
    a = ap.parse_args()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    ok = {'r1': run_r1, 'r2': run_r2}[a.gate](a.mem_cap_gb)
    sys.exit(0 if ok else 1)
