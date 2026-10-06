r"""Electron-nucleus (core-region) extension of the CI-relaxed RI-free LiH eigenproblem.

Backs (Paper 19, LiH explicit-r12 marriage): ansatz (A) of ci_relaxed.py -- {chi_I} U {(F - Fbar) psi}
-- extended by r12-FREE core-switched electron-nucleus columns (phase 1)

    W1 psi = sum_i w(r_i) psi,      W2 psi = sum_{i<j} w(r_i) w(r_j) psi,

and (phase 2, separate entry point) core-weighted geminal columns.  w(r) is the parameter-free core
switch of core_switch.py (the reference's own core natural-orbital density, w = 1 at the Li nucleus;
PI ruling 2026-10-05).  Diagnostic origin: debug/sprint_pair_geminal_step3_memo.md (VMC linear
method: with the geminal present the r12-free core columns are worth ~4-5 mHa).

REDUCTION.  Every new matrix element is an expectation / transition value of products of ONE-ELECTRON
grid weights (w, lap w, V_ne) with at most one pair kernel -- i.e. reweighted pair densities in the
existing pi_m16 / ci_relaxed reductions; no new kernel and no new reduction:
  * one-body weights multiply the vertex's own ordered-pair base:  base 'rho*w', 'grad*w*w', ...
    (ENRef.base; w is sigma, so the circulation labels mA are unchanged);
  * <V_ee> pieces go through gvee_sparse.vee_reduce_sparse on an EXTENDED pair space whose rows are
    the weight variants rho * w^k (VariantRef) -- one stacked coefficient tensor per operator.
Kinetic (u = bra, v = ket, B a real multiplicative operator; per-electron sums implied):
  transition (ci_relaxed's formula with F -> W):
      <u|T|W v> = 1/2 INT W grad u*.grad v - 1/4 INT (lap W) u* v + 1/4 INT W (u* lap v - v lap u*)
  same reference both sides (psi real up to a phase -> no current term):
      <A psi|T|B psi> = 1/2 INT A B |grad psi|^2 - 1/4 INT (A lap B + B lap A) |psi|^2
  with -1/4 lap F = -1/2 g^2 F + g Y  (exp geminal; Y = f/r, the pi_m16 Yukawa kernel).
The ci_relaxed (N+1) block is taken verbatim from ci_relaxed.ci_relaxed_energy (its numerics unchanged);
dropping the new columns reproduces it identically (gate G1).

API:
    switch_on_grid(info) -> (w, lapw, meta)
    en_columns(gem, info, Fbar, w, lapw) -> dict of the new S/H blocks
    en_energy(gem, refm, info, pieces) -> dict (E_A, E_A_en, sub-ansatz energies, S, H, identity checks)
Geminal kind: 'exp' only (asserted).
"""
from __future__ import annotations

import time
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg as sl

import geovac.lih_r12ci.enumerators as P2
import geovac.lih_r12ci.pi_m16 as M
import geovac.lih_r12ci.gvee_sparse as GS
import geovac.lih_r12ci.ci_relaxed as CR
from geovac.lih_r12ci.core_switch import CoreSwitch, grid_switch, grid_fraction_switch, select_core_no
from geovac.lih_r12ci.transition_tensor import (transition_stack, split_owners, determinants,
                                                bra_det_ket_psi)

PAIRS6 = P2.PAIRS6
GI_W = P2.GI_W
VNE = P2.VNE
V_NN = M.V_NN

from geovac.lih_r12ci.en_algebra import (Op, ONE, ELEC, W1_op, W2_op, op_add, op_mul, lap_op,  # noqa: E402,F401
                                          base_name, TOKCIRC, GOp, _edge, g_from, g_geminal, g_mul,
                                          g_lap, _route, op_mul_pairs, _single)


# =========================================================================== #
# reference objects
# =========================================================================== #
class ENRef(CR.FullPairRef):
    """FullPairRef + one-electron-weighted bases '<kind>*<tok>*...' (kind in rho/grad/rhov/lapA)."""

    def __init__(self, info: dict, tag: str, mults: Dict[str, np.ndarray]):
        super().__init__(info, tag)
        self.mults = dict(mults)
        self._xb: Dict[str, np.ndarray] = {}

    def base(self, name: str) -> np.ndarray:
        if '*' not in name:
            return super().base(name)
        if name not in self._xb:
            kind, *toks = name.split('*')
            b = super().base(kind).copy()
            for t in toks:
                b = b * self.mults[t][None, :]
            self._xb[name] = b
        return self._xb[name]


class VariantRef:
    """Extended ordered-pair space for gvee_sparse: row v*npair + A = rho_A * variant_v (grid weight)."""

    def __init__(self, base_ref: CR.FullPairRef, variants: Sequence[Tuple[str, np.ndarray]], tag: str):
        self.npair = base_ref.rho.shape[0]
        self.vnames = [n for n, _ in variants]
        self.vindex = {n: k for k, (n, _) in enumerate(variants)}
        self.rho = np.concatenate([base_ref.rho * arr[None, :] for _, arr in variants], axis=0)
        self.mA = np.tile(base_ref.mA, len(variants))
        self.nact = self.rho.shape[0]
        self.tag = tag
        self.T_iA = self.T_iB = self.T_iC = self.T_iD = np.zeros(0, int)
        self.T_val = np.zeros(0)

    def view(self, T_local) -> "VariantRef":
        v = object.__new__(VariantRef)
        v.__dict__.update(self.__dict__)
        v.T_iA, v.T_iB, v.T_iC, v.T_iD, v.T_val = T_local
        return v


def _toks_variant(toks: Tuple[str, ...]) -> str:
    return "1" if not toks else "*".join(toks)


def stacked_operator(op: Op, stack, vref: VariantRef):
    """Stacked COO (iA..iD, val, owner) on the variant space: sum over the operator's terms of the
    base stack with each electron's index shifted to its weight variant and val scaled by the coef."""
    iA, iB, iC, iD, val, ow = stack
    parts = []
    for t, c in op.items():
        sh = [vref.vindex[_toks_variant(t[e])] * vref.npair for e in range(4)]
        parts.append((iA + sh[0], iB + sh[1], iC + sh[2], iD + sh[3], c * val, ow))
    return tuple(np.concatenate([p[k] for p in parts]) for k in range(6))


def vee_operator(gem, vref: VariantRef, op: Op, stack, n_owner: int, with_F: bool,
                 mem_cap_gb: float = 6.0, ctx: Optional["GS.Ctx"] = None) -> np.ndarray:
    """Per-owner <op V_ee> (with_F False) or <op F V_ee> (with_F True; exp geminal, Yukawa when the f
    edge and the Coulomb pair coincide -- the ci_relaxed._vee_owner rule)."""
    st = stacked_operator(op, stack, vref)
    parts = split_owners(st, n_owner)
    if ctx is None:
        ctx = GS.Ctx(gem, vref, mem_cap_gb)
    col = GS.TriCollector()
    out = np.zeros(n_owner)
    for o, T in enumerate(parts):
        if len(T[4]) == 0:
            continue
        ctx.refm = vref.view(T)
        if not with_F:
            out[o] = sum(GS.vee_reduce_sparse([], r, 'coul', None, ctx, col) for r in PAIRS6)
        else:
            tot = 0.0
            for p in PAIRS6:
                for r in PAIRS6:
                    if frozenset(p) == frozenset(r):
                        tot += GS.vee_reduce_sparse([], r, 'yuk', gem.g, ctx, col)
                    else:
                        tot += GS.vee_reduce_sparse([tuple(p)], r, 'coul', None, ctx, col)
            out[o] = tot
    if col.recs:
        raise RuntimeError("unexpected triangle in a one-edge e-n V_ee element")
    ctx.refm = vref
    return out


# =========================================================================== #
# diagram evaluator over two stacks (transition chi_I|.|psi and reference psi|.|psi)
# =========================================================================== #
def diagram_sparse(edges, base_names: Dict[int, str], gem, ref: "ENRef", stack, n_owner: int,
                   chunk: int = 4096) -> np.ndarray:
    """Owner-resolved value of the diagram (== ci_relaxed.contract_owners(ci_relaxed.diagram_factors(...)))
    evaluated SPARSELY: each integrated vertex (own base x leaf dressings, L_z-masked) only at the index
    tuples that occur in the stack (no dense nact^k factor tensors).  Same reduction order / dressings."""
    bn = {e: base_names.get(e, 'rho') for e in (1, 2, 3, 4)}
    dressings: Dict[int, List[Tuple[np.ndarray, int]]] = {e: [] for e in (1, 2, 3, 4)}
    edge_list = [[set((i, j)), key] for (i, j, key) in edges]
    alive = {1, 2, 3, 4}
    vertices = []

    def dnb(u):
        s = set()
        for (ed, _k) in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    while alive:
        u = next((c for c in sorted(alive) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(alive) if len(dnb(c)) == 1 and not dressings[c]), None)
        if u is None:
            raise RuntimeError(f"unexpected topology: {edges}")
        inc = [(ed, k) for (ed, k) in edge_list if u in ed]
        nb = dnb(u)
        if not nb:
            vertices.append((u, ref.base(bn[u]), list(dressings[u])))
        else:
            if dressings[u]:
                raise RuntimeError(f"leaf {u} carries dressings (chain>2): {edges}")
            w_ = next(iter(nb))
            keys = [k for (ed, k) in inc]
            keys = keys[0] if len(keys) == 1 else tuple(keys)
            dressings[w_].append((M._refm_dressed(ref, gem, keys, bn[u]), u))
        alive.discard(u)
        edge_list = [[ed, k] for (ed, k) in edge_list if u not in ed]
    iA, iB, iC, iD, val, ow = stack
    idxmap = {1: iA, 2: iB, 3: iC, 4: iD}
    res = val.copy()
    npair = ref.rho.shape[0]
    for (u, own, dl) in vertices:
        axes = [u] + [e for (_p, e) in dl]
        tup = np.stack([idxmap[e] for e in axes], axis=1).astype(np.int64)
        code = np.ravel_multi_index(tuple(tup.T), (npair,) * len(axes))
        uc, inv = np.unique(code, return_inverse=True)
        ut = np.stack(np.unravel_index(uc, (npair,) * len(axes)), axis=1)
        vals = np.empty(len(uc))
        for s in range(0, len(uc), chunk):
            t = ut[s:s + chunk]
            rows = own[t[:, 0]] * GI_W[None, :]
            for k, (P, _e) in enumerate(dl):
                rows = rows * P[t[:, k + 1]]
            vals[s:s + len(t)] = rows.sum(axis=1)
        circ = ref.mA[ut].sum(axis=1)
        vals = vals * (circ == 0)
        res = res * vals[inv.reshape(-1)]
    return np.bincount(ow, weights=res, minlength=n_owner)


class _Eval:
    def __init__(self, gem, ref: ENRef, stI, N: int, stP):
        self.gem = gem; self.ref = ref; self.stI = stI; self.N = N; self.stP = stP
        self._res: Dict[tuple, Tuple[np.ndarray, float]] = {}

    def value(self, edges, bn: Dict[int, str], which: str = 'both') -> Tuple[np.ndarray, float]:
        """(per-owner transition values, reference value); `which` in {'I','P','both'} limits the work
        (the skipped one is returned as NaN)."""
        key = (tuple(edges), tuple(sorted((e, b) for e, b in bn.items() if b != 'rho')))
        vI, vP = self._res.get(key, (None, None))
        if vI is None and which in ('I', 'both'):
            vI = diagram_sparse(list(edges), dict(bn), self.gem, self.ref, self.stI, self.N)
        if vP is None and which in ('P', 'both'):
            vP = float(diagram_sparse(list(edges), dict(bn), self.gem, self.ref, self.stP, 1)[0])
        self._res[key] = (vI, vP)
        return (vI if vI is not None else np.full(self.N, np.nan),
                vP if vP is not None else float('nan'))

    def op(self, A: Op, kind: Optional[str] = None, edges=()) -> Tuple[np.ndarray, float]:
        """< A * (sum_i kind_i) * prod(edges) >  (kind None: no one-body kind; else summed over i)."""
        accI = np.zeros(self.N); accP = 0.0
        elecs = (None,) if kind is None else ELEC
        for t, c in A.items():
            for ei in elecs:
                bn = {}
                for e in ELEC:
                    k = kind if (ei is not None and e == ei) else 'rho'
                    bn[e] = base_name(k, t[e - 1])
                vI, vP = self.value(edges, bn)
                accI += c * vI; accP += c * vP
        return accI, accP

    def opF(self, A: Op, kind: Optional[str] = None, key: str = 'f') -> Tuple[np.ndarray, float]:
        accI = np.zeros(self.N); accP = 0.0
        for (i, j) in PAIRS6:
            vI, vP = self.op(A, kind, edges=((i, j, key),))
            accI += vI; accP += vP
        return accI, accP


def bra_psi_ket_psi(ia, ib, d):
    N = len(d)
    bra = [(ia[I], ib[I], d[I]) for I in range(N) for J in range(N)]
    ket = [(ia[J], ib[J], d[J]) for I in range(N) for J in range(N)]
    return bra, ket, [0] * (N * N)


# =========================================================================== #
# the switch on the grid
# =========================================================================== #
def switch_on_grid(info: dict, kind: str = 'dens') -> Tuple[np.ndarray, np.ndarray, dict]:
    """w, lap w on the kernels grid.  kind 'dens' (the ruled switch) or 'frac' (variant)."""
    from geovac.lih_r12ci.pi_m16 import build_primitives, CFG, PRIM_CACHE_PI
    prims = build_primitives(CFG, PRIM_CACHE_PI, verbose=False)['prims']
    core = select_core_no(info['occ'], info['m_no'], prims, info['T_no'], info['R'])
    sw = CoreSwitch(prims, info['T_no'][:, core], info['R'])
    if kind == 'dens':
        w, lw = grid_switch(info['nval'][core], info['ngrad'][core], info['nlap'][core], sw.phi0)
    elif kind == 'frac':
        w, lw = grid_fraction_switch(info['nval'], info['ngrad'], info['nlap'], info['occ'], core)
    else:
        raise ValueError(kind)
    return w, lw, dict(core_no=core, occ_core=float(info['occ'][core]), phi0=sw.phi0, kind=kind)


# =========================================================================== #
# phase 1: the W1 / W2 columns
# =========================================================================== #
def en_columns(gem, info: dict, Fbar: float, w: np.ndarray, lapw: np.ndarray, mem_cap_gb: float = 6.0,
               verbose: bool = False) -> dict:
    """New S / H blocks for the columns omega_a = W_a psi (a = 1, 2).  Electronic H (no V_NN)."""
    assert gem.name == 'exp', "en_factors supports the exp geminal only"
    t0 = time.time()
    ref = ENRef(info, "M=16 pi full ordered pairs + e-n weights", dict(w=w, lapw=lapw))
    ia, ib, d, apairs = determinants(info['ct'], int(info['Mk']))
    N = len(d); Mk = ref.Mk
    stI = transition_stack(*bra_det_ket_psi(ia, ib, d), apairs, Mk)
    stP = transition_stack(*bra_psi_ket_psi(ia, ib, d), apairs, Mk)
    E = _Eval(gem, ref, stI, N, stP)
    g = gem.g
    W = [W1_op(), W2_op()]
    LW = [lap_op(x) for x in W]
    na = len(W)
    out: dict = dict(N=N, d=d)
    # ---- transition columns <chi_I| . |psi> ----
    S_I = np.zeros((N, na)); T_I = np.zeros((N, na)); Vne_I = np.zeros((N, na)); LA_I = np.zeros((N, na))
    S_P = np.zeros(na)
    for a in range(na):
        S_I[:, a], S_P[a] = E.op(W[a])
        gI, _ = E.op(W[a], 'grad')
        lI, _ = E.op(LW[a])
        aI, _ = E.op(W[a], 'lapA')
        T_I[:, a] = 0.5 * gI - 0.25 * lI + 0.25 * aI
        LA_I[:, a] = aI
        Vne_I[:, a], _ = E.op(W[a], 'rhov')
    if verbose:
        print(f"    [en] transition one-body columns  [{time.time()-t0:6.1f}s]", flush=True)
    vref = VariantRef(ref, [("1", np.ones_like(w)), ("w", w), ("w*w", w * w)], ref.tag + " variants")
    ctx = GS.Ctx(gem, vref, mem_cap_gb)
    Vee_I = np.zeros((N, na))
    for a in range(na):
        Vee_I[:, a] = vee_operator(gem, vref, W[a], stI, N, False, ctx=ctx)
    if verbose:
        print(f"    [en] transition V_ee columns  [{time.time()-t0:6.1f}s]", flush=True)
    H_I = T_I + Vne_I + Vee_I
    # ---- reference-side blocks ----
    S_F = np.zeros(na); H_F = np.zeros(na)
    S_ab = np.zeros((na, na)); H_ab = np.zeros((na, na))
    parts_F = {}
    for a in range(na):
        _, FW = E.opF(W[a])
        S_F[a] = FW - Fbar * S_P[a]
        _, FWg = E.opF(W[a], 'grad')
        _, YW = E.opF(W[a], key='Y')
        _, FLW = E.opF(LW[a])
        _, FWv = E.opF(W[a], 'rhov')
        FWee = vee_operator(gem, vref, W[a], stP, 1, True, ctx=ctx)[0]
        kin = 0.5 * FWg + (-0.5 * g ** 2 * FW + g * YW) - 0.25 * FLW
        HFW = kin + FWv + FWee                                   # <F psi|H|W psi>
        HPW = float(d @ H_I[:, a])                               # <psi|H|W psi>
        H_F[a] = HFW - Fbar * HPW
        parts_F[a] = dict(FW=FW, FWg=FWg, YW=YW, FLW=FLW, FWv=FWv, FWee=FWee, kin=kin, HPW=HPW)
    if verbose:
        print(f"    [en] F-W blocks  [{time.time()-t0:6.1f}s]", flush=True)
    parts_ab = {}
    for a in range(na):
        for b in range(a, na):
            AB = op_mul(W[a], W[b])
            _, s = E.op(AB)
            _, gg = E.op(AB, 'grad')
            _, l1 = E.op(op_add(op_mul(W[a], LW[b]), op_mul(W[b], LW[a])))
            _, v = E.op(AB, 'rhov')
            vee = vee_operator(gem, vref, AB, stP, 1, False, ctx=ctx)[0]
            S_ab[a, b] = S_ab[b, a] = s
            H_ab[a, b] = H_ab[b, a] = 0.5 * gg - 0.25 * l1 + v + vee
            parts_ab[(a, b)] = dict(S=s, grad=gg, lap=l1, Vne=v, Vee=vee)
    if verbose:
        print(f"    [en] W-W blocks  [{time.time()-t0:6.1f}s]", flush=True)
    out.update(S_I=S_I, H_I=H_I, T_I=T_I, Vne_I=Vne_I, Vee_I=Vee_I, LA_I=LA_I, S_P=S_P,
               S_F=S_F, H_F=H_F, S_ab=S_ab, H_ab=H_ab, parts_F=parts_F, parts_ab=parts_ab,
               wall_s=time.time() - t0, E=E, ref=ref)
    return out


def lowest_root(H: np.ndarray, S: np.ndarray, tol: float = 1e-10) -> Tuple[float, np.ndarray]:
    """Lowest generalized eigenpair by canonical orthogonalisation of S (drops S eigenvalues
    < tol * max)."""
    Hs = 0.5 * (H + H.T); Ss = 0.5 * (S + S.T)
    s, U = np.linalg.eigh(Ss)
    keep = s > tol * s.max()
    X = U[:, keep] / np.sqrt(s[keep])
    e, V = np.linalg.eigh(X.T @ Hs @ X)
    return float(e[0]), X @ V[:, 0]


def assemble(cr: dict, en: dict) -> Tuple[np.ndarray, np.ndarray]:
    """(N+1+2) electronic S, H: ci_relaxed's (N+1) block + the W1, W2 columns (order chi, xi, W1, W2)."""
    N = cr['N']; K = N + 3
    S = np.zeros((K, K)); H = np.zeros((K, K))
    S[:N + 1, :N + 1] = cr['S']; H[:N + 1, :N + 1] = cr['H']
    S[:N, N + 1:] = en['S_I']; S[N + 1:, :N] = en['S_I'].T
    H[:N, N + 1:] = en['H_I']; H[N + 1:, :N] = en['H_I'].T
    S[N, N + 1:] = S[N + 1:, N] = en['S_F']
    H[N, N + 1:] = H[N + 1:, N] = en['H_F']
    S[N + 1:, N + 1:] = en['S_ab']; H[N + 1:, N + 1:] = en['H_ab']
    return S, H


def en_energy(gem, refm, info: dict, pieces: dict, switch: str = 'dens', mem_cap_gb: float = 6.0,
              verbose: bool = False, cr: Optional[dict] = None) -> dict:
    """Phase 1: ci_relaxed (A) + {W1 psi, W2 psi}.  Returns total energies (Ha) of the sub-ansatzes
    and the identity diagnostics (G1 and the internal consistency checks)."""
    if cr is None:
        cr = CR.ci_relaxed_energy(gem, refm, info, pieces, mem_cap_gb=mem_cap_gb, verbose=verbose)
    w, lw, meta = switch_on_grid(info, switch)
    en = en_columns(gem, info, float(pieces['Fbar']), w, lw, mem_cap_gb=mem_cap_gb, verbose=verbose)
    S, H = assemble(cr, en)
    N = cr['N']
    idx = dict(A=list(range(N + 1)), A_W1=list(range(N + 2)), A_W2=list(range(N + 1)) + [N + 2],
               A_W12=list(range(N + 3)), CI_W12=list(range(N)) + [N + 1, N + 2])
    E = {}
    vecs = {}
    for k, ii in idx.items():
        e, v = lowest_root(H[np.ix_(ii, ii)], S[np.ix_(ii, ii)])
        E[k] = e + V_NN; vecs[k] = (ii, v)
    d = en['d']
    ident = dict(
        G1_A_vs_ci_relaxed=E['A'] - cr['E_A'],
        S_P_vs_dS_I=float(d @ en['S_I'][:, 0]) - en['S_P'][0],
        LA_ref_W1=float(d @ en['LA_I'][:, 0]), LA_ref_W2=float(d @ en['LA_I'][:, 1]),
        HPW_vs_kinsym_W1=float(d @ en['H_I'][:, 0]) - _sym_HPW(en, 0),
        HPW_vs_kinsym_W2=float(d @ en['H_I'][:, 1]) - _sym_HPW(en, 1),
    )
    return dict(E=E, vecs=vecs, S=S, H=H, cr=cr, en=en, switch=meta, identity=ident, w=w, lapw=lw)


def _sym_HPW(en: dict, a: int) -> float:
    """<psi|H|W_a psi> from the bra=ket symmetric kinetic formula (A = 1): cross-check of the
    transition column contracted with d (the current term must vanish on the reference)."""
    E: _Eval = en['E']
    Wa = [W1_op(), W2_op()][a]
    _, gg = E.op(Wa, 'grad')
    _, l1 = E.op(lap_op(Wa))
    _, v = E.op(Wa, 'rhov')
    return 0.5 * gg - 0.25 * l1 + v + float(en['d'] @ en['Vee_I'][:, a])


# =========================================================================== #
# =====================  PHASE 2: core-weighted geminals  =================== #
# =========================================================================== #
# Columns Fcc psi = sum_{i<j} w_i w_j f_ij psi and Fc1 psi = sum_{i<j} (w_i + w_j) f_ij psi.  With F
# already in the basis, {Fcc, Fc1} spans the step-3 class geminals {F_cc, F_cv, F_vv} (F_cv = Fc1 - 2 Fcc,
# F_vv = F - Fc1 + Fcc).
#
# Generalised operator: {(toks4, edges): coef}, edges = sorted tuple of (i, j, key), key in
# {'f', 'Y', 'coul'} (f = e^{-g r}, Y = e^{-g r}/r).  Laplacian of a weighted pair term a(w) f_kl:
#   lap(a f) = (lap a) f + 2 sum_k (da/dw_k) grad w_k . grad_k f + a (lap_k + lap_l) f,
#   grad_k f = -g Y (r_k - r_l),   grad w(r_k).(r_k - r_l) = u(r_k) - dzw(r_k) z_l
#                                                      - 1/2 drw(r_k) rho_l (e^{i(phi_k-phi_l)} + c.c.),
# u = rho d_rho w + z d_z w.  The e^{+-i phi} factors are one-electron weights carrying circulation
# +-1 (tokens drwP/drwM, rhoP/rhoM): they shift the pair's L_z label, which the extended-pair-space
# route (gvee_sparse + VariantRef2 with mA offsets) handles exactly; no new kernel is needed.
class VariantRef2:
    """Extended ordered-pair space: row v*npair + A = rho_A * variant_v, circulation mA_A + shift_v."""

    def __init__(self, base_ref: CR.FullPairRef, arrays: Dict[str, np.ndarray], toks_list, tag: str):
        self.npair = base_ref.rho.shape[0]
        self.vindex = {}
        rows = []; mAs = []
        for tt in toks_list:
            arr = np.ones(base_ref.rho.shape[1]); sh = 0
            for x in tt:
                arr = arr * arrays[x]; sh += TOKCIRC[x]
            self.vindex[tt] = len(rows)
            rows.append(base_ref.rho * arr[None, :]); mAs.append(base_ref.mA + sh)
        self.rho = np.concatenate(rows, axis=0)
        self.mA = np.concatenate(mAs).astype(int)
        self.nact = self.rho.shape[0]
        self.tag = tag
        self.T_iA = self.T_iB = self.T_iC = self.T_iD = np.zeros(0, int)
        self.T_val = np.zeros(0)

    def view(self, T_local) -> "VariantRef2":
        v = object.__new__(VariantRef2)
        v.__dict__.update(self.__dict__)
        v.T_iA, v.T_iB, v.T_iC, v.T_iD, v.T_val = T_local
        return v


def _gs_struct(edges, g: float):
    """(rem_f, kept, kernel_kind, gam) for vee_reduce_sparse from an edge multiset with ONE non-f edge."""
    kept = [e for e in edges if e[2] != 'f']
    assert len(kept) == 1, edges
    (ki, kj, kk) = kept[0]
    fs = [(i, j) for (i, j, k) in edges if k == 'f']
    n_on = sum(1 for (i, j) in fs if (i, j) == (ki, kj))
    rem = [(i, j) for (i, j) in fs if (i, j) != (ki, kj)]
    if kk == 'coul':
        kind = 'coul' if n_on == 0 else 'yuk'
        gam = None if n_on == 0 else n_on * g
        nexp = n_on
    else:                                                       # 'Y' = Yukawa(g)
        kind = 'yuk'; gam = (1 + n_on) * g
        nexp = 1 + n_on
    assert nexp <= 2, "kernel exponent > 2g not supported for m != 0 modes"
    return tuple(rem), (ki, kj), kind, gam


class Phase2Eval:
    """Evaluate GOp expectations over the transition stack (owners I) and the reference stack."""

    def __init__(self, gem, ref: ENRef, E1: _Eval, arrays: Dict[str, np.ndarray], mem_cap_gb: float = 6.0):
        self.gem = gem; self.ref = ref; self.E1 = E1; self.arrays = arrays
        self.cap = mem_cap_gb
        self.vref = None; self.ctx = None

    def prepare(self, ops: Sequence[GOp]) -> None:
        toks = {()}
        for op in ops:
            for term in op:
                if _route(term) == 'G':
                    for tt in term[0]:
                        toks.add(tuple(tt))
        self.toks = sorted(toks, key=lambda x: (len(x), x))
        self.vref = VariantRef2(self.ref, self.arrays, self.toks, self.ref.tag + " phase2 variants")
        self.ctx = GS.Ctx(self.gem, self.vref, self.cap)

    def _gs(self, terms: List[Tuple[tuple, float]], stack, n_owner: int) -> np.ndarray:
        """GS-route terms, grouped by edge structure; per-owner vee_reduce_sparse (+ triangles)."""
        iA, iB, iC, iD, val, ow = stack
        groups: Dict[tuple, list] = {}
        for (toks, edges), c in terms:
            groups.setdefault(_gs_struct(edges, self.gem.g), []).append((toks, c))
        out = np.zeros(n_owner)
        npair = self.vref.npair
        for (rem, kept, kind, gam), lst in groups.items():
            parts = []
            for toks, c in lst:
                sh = [self.vref.vindex[tuple(toks[e])] * npair for e in range(4)]
                parts.append((iA + sh[0], iB + sh[1], iC + sh[2], iD + sh[3], c * val, ow))
            st = tuple(np.concatenate([p[k] for p in parts]) for k in range(6))
            for o, T in enumerate(split_owners(st, n_owner)):
                if len(T[4]) == 0:
                    continue
                self.ctx.refm = self.vref.view(T)
                col = GS.TriCollector()
                v = GS.vee_reduce_sparse(list(rem), kept, kind, gam, self.ctx, col)
                if col.recs:
                    v += GS.eval_triangles(self.ctx, col)
                out[o] += v
        self.ctx.refm = self.vref
        return out

    def ev(self, op: GOp, kind: Optional[str] = None, which: str = 'P') -> np.ndarray:
        """<op * sum_i kind_i> per owner ('I': transition, N owners) or reference ('P': length 1)."""
        E = self.E1
        st, n = (E.stI, E.N) if which == 'I' else (E.stP, 1)
        acc = np.zeros(n)
        gterms = []
        for term, c in op.items():
            toks, edges = term
            if _route(term) == 'G':
                assert kind is None, "one-body kinds are not supported on the GS route"
                gterms.append((term, c))
                continue
            elecs = (None,) if kind is None else ELEC
            for ei in elecs:
                bn = {e: base_name(kind if (ei is not None and e == ei) else 'rho', toks[e - 1])
                      for e in ELEC}
                vI, vP = E.value(list(edges), bn, which=which)
                acc += c * (vI if which == 'I' else np.array([vP]))
        if gterms:
            acc += self._gs(gterms, st, n)
        return acc

    def vee_batch(self, ops: Dict[str, GOp]) -> Dict[str, float]:
        """Reference-side <op V_ee> for several operators, the triangle class evaluated ONCE per middle
        vertex for the union of all operators' (left, right) index pairs (tri_entries is T-independent)."""
        iA, iB, iC, iD, val, ow = self.E1.stP
        npair = self.vref.npair
        out: Dict[str, float] = {}
        recs: Dict[str, GS.TriCollector] = {}
        for name, op in ops.items():
            tot: GOp = {}
            for (toks, edges), c in op.items():
                for (i, j) in PAIRS6:
                    k = (toks, tuple(sorted(edges + (_edge(i, j, 'coul'),))))
                    tot[k] = tot.get(k, 0.0) + c
            groups: Dict[tuple, list] = {}
            for (toks, edges), c in tot.items():
                groups.setdefault(_gs_struct(edges, self.gem.g), []).append((toks, c))
            col = GS.TriCollector(); acc = 0.0
            for (rem, kept, kind, gam), lst in groups.items():
                parts = []
                for toks, c in lst:
                    sh = [self.vref.vindex[tuple(toks[e])] * npair for e in range(4)]
                    parts.append((iA + sh[0], iB + sh[1], iC + sh[2], iD + sh[3], c * val))
                T = tuple(np.concatenate([p[k] for p in parts]) for k in range(5))
                self.ctx.refm = self.vref.view(T)
                acc += GS.vee_reduce_sparse(list(rem), kept, kind, gam, self.ctx, col)
            out[name] = acc; recs[name] = col
        # union triangle evaluation
        n = self.vref.nact
        order: List[int] = []; seen = set()
        per = {}
        for name, col in recs.items():
            if not col.recs:
                continue
            im = np.concatenate([r[0] for r in col.recs]); il = np.concatenate([r[1] for r in col.recs])
            ir = np.concatenate([r[2] for r in col.recs]); w_ = np.concatenate([r[3] for r in col.recs])
            per[name] = (im, il, ir, w_)
            for cm in col.order:
                if cm not in seen:
                    seen.add(cm); order.append(cm)
        if per:
            allc = np.unique(np.concatenate([(p[0] * n + p[1]) * n + p[2] for p in per.values()]))
            ucm = allc // (n * n); uil = (allc // n) % n; uir = allc % n
            vals = np.zeros(len(allc))
            for cm in order:
                s = ucm == cm
                vals[s] = GS.tri_entries(self.ctx, cm, uil[s], uir[s])
            for name, (im, il, ir, w_) in per.items():
                pos = np.searchsorted(allc, (im * n + il) * n + ir)
                out[name] += float(np.sum(w_ * vals[pos]))
        self.ctx.refm = self.vref
        return out

    def vee(self, op: GOp, which: str = 'P') -> np.ndarray:
        """<op V_ee>: a Coulomb edge on every pair, through the GS route."""
        tot: GOp = {}
        for (toks, edges), c in op.items():
            for (i, j) in PAIRS6:
                k = (toks, tuple(sorted(edges + (_edge(i, j, 'coul'),))))
                tot[k] = tot.get(k, 0.0) + c
        E = self.E1
        st, n = (E.stI, E.N) if which == 'I' else (E.stP, 1)
        return self._gs(list(tot.items()), st, n)


def grid_geometry(info: dict, core: int, phi0: float) -> Dict[str, np.ndarray]:
    """One-electron grid arrays for the phase-2 tokens (pts = (rho, 0, z) at phi_azi = 0)."""
    pts = info['grid']['pts']
    rho = pts[:, 0]; z = pts[:, 2]
    ph = np.asarray(info['nval'][core]).real
    gph = np.asarray(info['ngrad'][core]).real
    gw = 2.0 * ph[:, None] * gph / phi0 ** 2
    return dict(u=rho * gw[:, 0] + z * gw[:, 2], dzw=gw[:, 2], z=z, drwP=gw[:, 0], drwM=gw[:, 0],
                rhoP=rho, rhoM=rho)


def en2_columns(gem, info: dict, Fbar: float, w: np.ndarray, lapw: np.ndarray, core: int, phi0: float,
                mem_cap_gb: float = 6.0, verbose: bool = False, en1: Optional[dict] = None) -> dict:
    """Phase-2 blocks for the columns G psi, G in {Fcc, Fc1} (order after W1, W2)."""
    assert gem.name == 'exp'
    t0 = time.time()
    if en1 is None:
        en1 = en_columns(gem, info, Fbar, w, lapw, mem_cap_gb=mem_cap_gb, verbose=verbose)
    E1: _Eval = en1['E']; ref: ENRef = en1['ref']; d = en1['d']; N = en1['N']
    g = gem.g
    arrays = dict(w=w, lapw=lapw, **grid_geometry(info, core, phi0))
    P2E = Phase2Eval(gem, ref, E1, arrays, mem_cap_gb)
    Fop = g_geminal('F'); G = [g_geminal('cc'), g_geminal('c1')]
    Wop = [g_from(W1_op()), g_from(W2_op())]
    lapF = g_lap(Fop, g); lapG = [g_lap(x, g) for x in G]; lapW = [g_lap(x, g) for x in Wop]
    others = [('F', Fop, lapF), ('W1', Wop[0], lapW[0]), ('W2', Wop[1], lapW[1]),
              ('Fcc', G[0], lapG[0]), ('Fc1', G[1], lapG[1])]
    allops = list(lapG)
    for (_n, A, LA) in others:
        for b in range(2):
            allops += [g_mul(A, lapG[b]), g_mul(G[b], LA)]
    for b in range(2):                                          # V_ee operators (any coul pair)
        for (_n, A, _LA) in others + [('1', g_from(ONE), {})]:
            AB = g_mul(A, G[b])
            allops.append({(t, tuple(sorted(e + (_edge(1, 2, 'coul'),)))): c for (t, e), c in AB.items()})
    P2E.prepare(allops)
    if verbose:
        print(f"    [en2] variants: {P2E.toks}  nact={P2E.vref.nact}  [{time.time()-t0:6.1f}s]", flush=True)
    nb = 2
    S_I = np.zeros((N, nb)); H_I = np.zeros((N, nb)); parts_I = {}
    for b in range(nb):
        S_I[:, b] = P2E.ev(G[b], which='I')
        gg = P2E.ev(G[b], 'grad', 'I'); lG = P2E.ev(lapG[b], which='I'); la = P2E.ev(G[b], 'lapA', 'I')
        vne = P2E.ev(G[b], 'rhov', 'I'); vee = P2E.vee(G[b], 'I')
        T = 0.5 * gg - 0.25 * lG + 0.25 * la
        H_I[:, b] = T + vne + vee
        parts_I[b] = dict(T=T, Vne=vne, Vee=vee, LA=la, grad=gg, lap=lG)
        if verbose:
            print(f"    [en2] transition column G{b}  [{time.time()-t0:6.1f}s]", flush=True)
    HPG = d @ H_I                                              # <psi|H|G psi>
    SPG = d @ S_I
    pairs_blk = [(na_, A, LA, b) for (na_, A, LA) in others for b in range(nb)
                 if not (na_ in ('Fcc', 'Fc1') and ['Fcc', 'Fc1'].index(na_) > b)]
    veeb = P2E.vee_batch({f"{na_}|{b}": g_mul(A, G[b]) for (na_, A, _LA, b) in pairs_blk})
    if verbose:
        print(f"    [en2] reference V_ee blocks (batched triangles)  [{time.time()-t0:6.1f}s]", flush=True)
    blk = {}
    for (na_, A, LA, b) in pairs_blk:
        AB = g_mul(A, G[b])
        s = float(P2E.ev(AB)[0])
        gg = float(P2E.ev(AB, 'grad')[0])
        l1 = float(P2E.ev(g_mul(A, lapG[b]))[0]); l2 = float(P2E.ev(g_mul(G[b], LA))[0])
        v = float(P2E.ev(AB, 'rhov')[0]); vee = float(veeb[f"{na_}|{b}"])
        h = 0.5 * gg - 0.25 * (l1 + l2) + v + vee
        blk[(na_, b)] = dict(S=s, H=h, grad=gg, lapAB=l1, lapBA=l2, Vne=v, Vee=vee)
        if verbose:
            print(f"    [en2] <{na_}|G{b}>: S={s:+.8f} H={h:+.8f}  [{time.time()-t0:6.1f}s]", flush=True)
    return dict(S_I=S_I, H_I=H_I, parts_I=parts_I, HPG=HPG, SPG=SPG, blk=blk, Fbar=Fbar,
                wall_s=time.time() - t0, variants=P2E.toks)


def assemble2(cr: dict, en: dict, en2: dict) -> Tuple[np.ndarray, np.ndarray]:
    """(N+5) electronic S, H in the order [chi_I, xi, W1, W2, Fcc, Fc1]."""
    S1, H1 = assemble(cr, en)
    N = cr['N']; K = N + 5; Fbar = en2['Fbar']
    S = np.zeros((K, K)); H = np.zeros((K, K))
    S[:N + 3, :N + 3] = S1; H[:N + 3, :N + 3] = H1
    S[:N, N + 3:] = en2['S_I']; S[N + 3:, :N] = en2['S_I'].T
    H[:N, N + 3:] = en2['H_I']; H[N + 3:, :N] = en2['H_I'].T
    blk = en2['blk']
    for b in range(2):
        c = N + 3 + b
        S[N, c] = S[c, N] = blk[('F', b)]['S'] - Fbar * en2['SPG'][b]
        H[N, c] = H[c, N] = blk[('F', b)]['H'] - Fbar * en2['HPG'][b]
        for a, nm in enumerate(('W1', 'W2')):
            S[N + 1 + a, c] = S[c, N + 1 + a] = blk[(nm, b)]['S']
            H[N + 1 + a, c] = H[c, N + 1 + a] = blk[(nm, b)]['H']
    for (nm, b), v in blk.items():
        if nm in ('Fcc', 'Fc1'):
            a = ['Fcc', 'Fc1'].index(nm)
            S[N + 3 + a, N + 3 + b] = S[N + 3 + b, N + 3 + a] = v['S']
            H[N + 3 + a, N + 3 + b] = H[N + 3 + b, N + 3 + a] = v['H']
    return S, H


def en2_energy(gem, refm, info: dict, pieces: dict, mem_cap_gb: float = 6.0, verbose: bool = False,
               cr: Optional[dict] = None, ph1: Optional[dict] = None) -> dict:
    """Phase 2 (density switch): (A) + {W1, W2} + {Fcc, Fc1}.  Sub-ansatz energies (total, Ha)."""
    if ph1 is None:
        ph1 = en_energy(gem, refm, info, pieces, switch='dens', mem_cap_gb=mem_cap_gb, verbose=verbose, cr=cr)
    cr = ph1['cr']; en = ph1['en']; meta = ph1['switch']
    assert meta['kind'] == 'dens', "phase-2 grid geometry is built for the density switch"
    en2 = en2_columns(gem, info, float(pieces['Fbar']), ph1['w'], ph1['lapw'], meta['core_no'], meta['phi0'],
                      mem_cap_gb=mem_cap_gb, verbose=verbose, en1=en)
    S, H = assemble2(cr, en, en2)
    N = cr['N']
    idx = dict(A=list(range(N + 1)), A_W12=list(range(N + 3)), A_G=list(range(N + 1)) + [N + 3, N + 4],
               A_W12_G=list(range(N + 5)))
    E = {}; vecs = {}
    for k, ii in idx.items():
        e, v = lowest_root(H[np.ix_(ii, ii)], S[np.ix_(ii, ii)])
        E[k] = e + V_NN; vecs[k] = (ii, v)
    return dict(E=E, vecs=vecs, S=S, H=H, en2=en2, ph1=ph1)
