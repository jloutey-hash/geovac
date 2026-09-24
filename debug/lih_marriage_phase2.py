r"""Phase 2 of the LiH "marriage" build (2026-09-23; plan debug/lih_marriage_build_plan.md, memo
debug/sprint_lih_marriage_memo.md) -- the SEVEN enumerator generalisations over the pair-index
coefficient tensor T, and GATE G1 (regression on the 2-orbital reference).

The 2-MO assembly (geovac/lih_r12ci/assembly.py) is welded to |Phi0|^2 = D_p(1,2) D_p(3,4)/4 for the
single closed-shell reference |phi_0^2 phi_1^2|; every explicit-r12 energy piece is a sum over the AO
pair-density components (DP_COMPS / KD_COMPS / COMPS).  Audit B: the reference enters ONLY through the
4-index tensor T[A,B,C,D] with  P(1,2,3,4) = |Psi|^2/4 = sum_{ABCD} T[A,B,C,D] rho_A(1) rho_B(2)
rho_C(3) rho_D(4)  (electrons 1,2 alpha ; 3,4 beta).  Every primitive takes ONE separable product of
one-electron pair densities as an argument, and is multilinear in the four densities, so

    <O> = sum_{ABCD} T[A,B,C,D] . contract(O ; rho_A, rho_B, rho_C, rho_D)

with NO other reference dependence.  This module rewrites each enumerator as such a T-contraction and
runs GATE G1: the 2-MO reference fed through the generalised enumerators must reproduce the banked
per-piece analytic references (assembly.ANALYTIC_REF) and the assembled E_R12.

The seven generalisations (see assembly.py line refs in the plan):
  (1) CI vector -> T + spectator traces        -- DONE in Phase 0 (loaded, not rebuilt)
  (2) Fbar, sigma^2                             -- the 0/1/2-edge no-operator enumerator
  (3) h_Vne, g_Vne                              -- V_ne as a one-body multiplier on one electron
  (4) h_T, g_T                                  -- drift^2 (gradient pair-densities) + IBP PartB/gT23
  (5) g_Vee                                     -- the kept-Coulomb-pair reducer (+ RI-free triangle)
  (6) Cov[F,Y], Cov[F,E]                        -- the mixed 1-edge x kernel enumerator
  (7) E0                                        -- grid-consistent <T>+<V_ne>+<V_ee>+V_NN

Structural note (audit B): every diagram is reduced by message-passing to pair-space kernel matrices
W[A,B] = INT rho_A K rho_B and 3-index junctions J3[A,B,C] = INT rho_A Psi^K_B Psi^K'_C, then contracted
with T over its 38,786 nonzeros -- O(nnz) per piece, never the dense 55^4 tensor.

Run from root (Git Bash):
    python debug/lih_marriage_phase2.py > debug/data/lih_marriage_phase2.log 2>&1; echo $? > debug/data/lih_marriage_phase2.exit
Options: --preview  (also run the 74-det tensor through the enumerators, UNGATED)  --no-exp / --no-linexp
Writes debug/data/lih_marriage_phase2.npz  (the generalised pieces + G1 comparison table + preview).
Nothing in geovac/lih_r12ci/ or debug/lih_vmc.py is modified.
"""
from __future__ import annotations

import os
import sys
import time
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# Pin BLAS to one thread BEFORE numpy import (Phase-1 hit oversubscription crashes).
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
REF_EXACT = os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz")
OUT = os.path.join(DATA, "lih_marriage_phase2.npz")

# electron blocks: 1,2 alpha ; 3,4 beta.  The 6 interaction pairs:
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def _t(t0: float) -> str:
    return f"[{time.time() - t0:7.1f}s]"


# =========================================================================== #
# geovac primitives  (import once; toggle USE_EXACT_NEUMANN per path)
# =========================================================================== #
import geovac.lih_r12ci.kernels as KG
from geovac.lih_r12ci.kernels import (a as A_FOCUS, GAM, geo_f, grid_int, rho_cyl_f, zc_f,
                                       NXI, NETA, NG, build_kernel, Xg as _Xg, Eg as _Eg,
                                       ETA as _ETA1D, WETA as _WETA1D)
from geovac.lih_r12ci.hVee import neumann_potential, P00, P01, P11
from geovac.lih_r12ci.hT import Gpq, cvec, d_aa, d_ab, d_bb, rA_f, rB_f
from geovac.lih_r12ci.energy import R as R_ENG, Z_A, Z_B
from geovac.lih_r12ci import gVee as GV
from geovac.lih_r12ci.triangle import (build_kernel_m, triangle_raw, MMAX, _MODE_CACHE,
                                        Kmu_modes, coul_mode_potential, GW2 as _GW2,
                                        Peta as _PETA, JAC2 as _JAC2T, _NORM as _TRINORM,
                                        LMAX as _TRILMAX)
from scipy.special import eval_legendre as _eval_legendre

_JAC2 = (_Xg ** 2 - _Eg ** 2)                    # (NXI, NETA)  d3r Jacobian factor

A3 = A_FOCUS ** 3
GEO = geo_f
GI_W = 2.0 * np.pi * A3 * GEO           # grid_int(field) == GI_W @ field
VNE = (-Z_A / rA_f - Z_B / rB_f)        # one-body V_ne on the grid
V_NN = Z_A * Z_B / R_ENG


def gint(field: np.ndarray) -> float:
    return float(GI_W @ field)


# ---- geminal kernel builders (copied verbatim from assembly.py) ---- #
def build_kernel_gen(func: Callable, nphi: int = 28) -> np.ndarray:
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f; rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG, NG))
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 0.0))
        K += 2.0 * w * func(d)
    return K


def build_kernel_m_gen(func: Callable, m: int, nphi: int = 48) -> np.ndarray:
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f; rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG, NG))
    for ph, w in zip(phi, wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * np.cos(ph), 0.0))
        K += 2.0 * w * func(d) * np.cos(m * ph)
    return K / (2 * np.pi) if m == 0 else K / np.pi


def dress_mat(K: np.ndarray, stack: np.ndarray) -> np.ndarray:
    """a^3 (K @ (geo * h)) for a stack (npair, NG) -> (npair, NG)  (K symmetric)."""
    return A3 * ((GEO[None, :] * stack) @ K)


def neumann_stack(stack: np.ndarray, LMAX: int = 34) -> np.ndarray:
    """Coulomb (prolate-Neumann) potential of each row.  On the EXACT path this is BATCHED over the
    whole stack through the l-loop (op.radial is a matmul A[0,l] @ g, so it vectorises exactly);
    on the legacy path it falls back to the per-row neumann_potential (bit-identical to G1)."""
    stack = np.ascontiguousarray(stack)
    M = stack.shape[0]
    if not KG.USE_EXACT_NEUMANN:
        out = np.empty_like(stack)
        for i in range(M):
            out[i] = neumann_potential(stack[i].reshape(NXI, NETA)).reshape(-1)
        return out
    op = KG.exact_neumann(LMAX, 0)
    pref = (2.0 / R_ENG) * (2.0 * np.pi) * A3
    D = stack.reshape(M, NXI, NETA)
    W = _JAC2[None, :, :] * D                                  # (M, NXI, NETA)
    V = np.zeros((M, NXI, NETA))
    for l in range(LMAX + 1):
        Pl = _eval_legendre(l, _ETA1D)                         # (NETA,)
        g_l = (W * Pl[None, None, :]) @ _WETA1D                 # (M, NXI)
        radial = g_l @ op.A[0, l].T                            # (M, NXI)  == op.radial per row
        V += (2 * l + 1) * radial[:, :, None] * Pl[None, None, :]
    return (pref * V).reshape(M, NG)


def psi_yuk_stack(stack: np.ndarray, gam: float) -> np.ndarray:
    """Yukawa(gam) potential e^{-gam r}/r of each row = psi_coul - smooth (batched); == GV.psi_yuk."""
    return neumann_stack(stack) - A3 * ((GEO[None, :] * stack) @ GV._KS[gam])


# =========================================================================== #
# PHASE 3 g_Vee productionization: batched mode-m Coulomb + cached triangle matrices
# =========================================================================== #
# The Phase-2 triangle path recomputes coul_mode_potential once per (cm, m, rank, R-index)
# INSIDE each of the ~24-48 triangle triples of <F^2 V_ee>.  The mode-builds (Kmu_modes) are
# already shared across triples via triangle._MODE_CACHE (keyed by the middle vertex cm + mode).
# The Coulomb solves were NOT shared and were the wall-clock blocker on the 55-active-pair 74-det
# reference.  Fix (exact, NOT a math change -- verified <0.001 mHa vs the per-row path on the 2-MO
# reference before any 74-det run): (i) BATCH coul_mode_potential over all (rank x R-density) rows
# of one (cm, m) block into a single vectorised l-loop; (ii) CACHE the full (nact x nact) triangle
# matrix  M_cm[a,b] = sum_m sum_rr lam_rr <rho_a v_rr, C_m(rho_b v_rr)>_grid2  per middle vertex cm,
# so every triangle triple through cm just reads M_cm[il, ir].
_JAC2T_ = _JAC2T                                       # (NXI, NETA)
_WETA_STACK = _WETA1D                                  # 1D eta weights (kernels.WETA)


def coul_mode_stack(Bstack: np.ndarray, m: int) -> np.ndarray:
    """Batched EXACT mode-m prolate-Neumann Coulomb potential of each row of Bstack (nB, NG).

    Row-for-row identical to triangle.coul_mode_potential(B2d, m, exact=True); only the l-loop is
    vectorised over the nB inputs (one matmul g @ A[m,l].T per l instead of nB separate solves)."""
    if not KG.USE_EXACT_NEUMANN:
        # legacy fallback (used only if ever run off the exact path) -- per-row, bit-identical
        return np.stack([coul_mode_potential(Bstack[i].reshape(NXI, NETA), m).reshape(-1)
                         for i in range(Bstack.shape[0])])
    op = KG.exact_neumann(_TRILMAX, m)
    Wm = (2 * np.pi) ** 2 if m == 0 else 2 * np.pi ** 2
    pref = (2.0 - (1.0 if m == 0 else 0.0)) * (2.0 / R_ENG) * A_FOCUS ** 6 * Wm
    nB = Bstack.shape[0]
    B = Bstack.reshape(nB, NXI, NETA)
    WB = _JAC2T_[None, :, :] * B                        # (nB, NXI, NETA)
    V = np.zeros((nB, NXI, NETA))
    sign = (-1.0) ** m
    for l in range(m, _TRILMAX + 1):
        Plm_e = _PETA[m, l]                             # (NETA,)
        gB = (WB * Plm_e[None, None, :]) @ _WETA_STACK  # (nB, NXI)
        radial = gB @ op.A[m, l].T                      # (nB, NXI)  == op.radial(gB_row, l, m) per row
        V += (sign * (2 * l + 1) * _TRINORM[(m, l)]) * (radial[:, :, None] * Plm_e[None, None, :])
    return (pref * V).reshape(nB, NG)


_TRI_MAT_CACHE: Dict[tuple, np.ndarray] = {}


def get_tri_matrix(gem: "Geminal", ref: "Ref", cm: int) -> np.ndarray:
    """Cached (nact, nact) triangle matrix for middle-vertex density rho_cm:
        M[a, b] = sum_m sum_rr lam_rr  E_m(rho_a v_rr, rho_b v_rr)
                = sum_m sum_rr lam_rr  (rho_a v_rr) . GW2 . C_m(rho_b v_rr) .
    Reproduces triangle_raw's mode expansion exactly (batched over the rank x nact solves)."""
    ck = (ref.tag, gem.name, int(KG.USE_EXACT_NEUMANN), int(cm))
    if ck in _TRI_MAT_CACHE:
        return _TRI_MAT_CACHE[ck]
    nact = ref.nact
    rho = ref.rho                                       # (nact, NG)
    mu = rho[cm]                                         # (NG,)
    gw2f = _GW2.reshape(-1)                              # (NG,) grid2 weight
    key = f"{gem.name}_mu{int(cm)}_{int(KG.USE_EXACT_NEUMANN)}"
    M = np.zeros((nact, nact))
    for m in range(MMAX + 1):
        lam, Vmodes = Kmu_modes(gem.Fm[m], mu, m, key)  # (rank,), (NG, rank)
        rk = lam.size
        if rk == 0:
            continue
        # batch the coulomb over all (rank x nact) rows of rho_b * v_rr in one vectorised l-loop
        Bstack = (rho[None, :, :] * Vmodes.T[:, None, :]).reshape(rk * nact, NG)   # (rk*nact, NG)
        Pm = coul_mode_stack(Bstack, m).reshape(rk, nact, NG)                       # (rk, nact, NG)
        for rr in range(rk):
            Lm = rho * (Vmodes[:, rr] * gw2f)[None, :]                              # (nact, NG)
            M += lam[rr] * (Lm @ Pm[rr].T)                                          # (nact, nact)
    _TRI_MAT_CACHE[ck] = M
    return M


# =========================================================================== #
# geminal kernel bundle
# =========================================================================== #
class Geminal:
    """Kernels + dressing functions for one geminal (exp or linexp).  dress(key, stack)->stack.

    Kernel keys and their dressings (Psi^K of a density stack):
      'f'  the geminal f              'f2' f^2         'E' exp(-g r) (linexp only)
      'Y'  Yukawa(g) e^{-gr}/r        product keys 'fY','fE','ff'->'f2'
    Yukawa dressings use GV.psi_yuk (Neumann-based; reads USE_EXACT_NEUMANN).  On the EXACT path
    psi_yuk == the closed-radial yukawa_pot_iso to ~1e-13, so the enumerator matches the
    exact-kernel references; on the LEGACY path it carries the documented ~0.5% psi_yuk-vs-closed
    bias in the aa/bb self-terms (gT.py:92-102) -- reported, not a generalisation bug.
    """

    def __init__(self, name: str):
        self.name = name
        g = GAM
        if name == "exp":
            self.Kf = build_kernel(g); self.Kf2 = build_kernel(2 * g)
            self.Fm = [build_kernel_m(g, m) for m in range(MMAX + 1)]
        elif name == "linexp":
            f_lin = lambda d: d * np.exp(-g * d)
            self.Kf = build_kernel_gen(f_lin)
            self.Kf2 = build_kernel_gen(lambda d: (d * np.exp(-g * d)) ** 2)
            self.Fm = [build_kernel_m_gen(f_lin, m) for m in range(MMAX + 1)]
        else:
            raise ValueError(name)
        self.Kexp_g = build_kernel_gen(lambda d: np.exp(-g * d))       # e^{-g d}
        self.Kexp_2g = build_kernel_gen(lambda d: np.exp(-2 * g * d))  # e^{-2g d}
        self.Klin_2g = build_kernel_gen(lambda d: d * np.exp(-2 * g * d))  # d e^{-2g d}
        self.g = g

    # -- leaf dressing of a base stack by a single kernel key -- #
    def dress(self, key: str, stack: np.ndarray) -> np.ndarray:
        g = self.g
        if key == 'f':
            return dress_mat(self.Kf, stack)
        if key == 'f2':
            return dress_mat(self.Kf2, stack)
        if key == 'E':
            return dress_mat(self.Kexp_g, stack)
        if key == 'Y':
            return psi_yuk_stack(stack, g)
        raise ValueError(key)

    # -- product dressing (two edges on the same electron pair, different kernels) -- #
    def dress_prod(self, keys: Tuple[str, str], stack: np.ndarray) -> np.ndarray:
        ks = tuple(sorted(keys))
        g = self.g
        if ks == ('f', 'f'):
            return dress_mat(self.Kf2, stack)
        if ks == ('E', 'f'):        # f * E  (linexp): d e^{-gr} * e^{-gr} = d e^{-2gr}
            return dress_mat(self.Klin_2g, stack)
        if ks == ('Y', 'f'):        # f * Y
            if self.name == "exp":  # e^{-gr} * e^{-gr}/r = e^{-2gr}/r = Yukawa(2g)
                return psi_yuk_stack(stack, 2 * g)
            else:                    # d e^{-gr} * e^{-gr}/r = e^{-2gr}
                return dress_mat(self.Kexp_2g, stack)
        raise ValueError(ks)

    # -- kept-pair Coulomb potential of Rt (Coulomb / Yukawa; geminal-specific yuk) -- #
    def coul_pot(self, stack: np.ndarray) -> np.ndarray:
        return neumann_stack(stack)

    def yuk_pot(self, stack: np.ndarray, gam_k: float) -> np.ndarray:
        g = self.g
        if self.name == "exp":       # f/r -> Yukawa(gam_k)
            return psi_yuk_stack(stack, gam_k)
        # linexp: f/r -> e^{-g r} (Kexp_g);  f^2/r -> d e^{-2g d} (Klin_2g)
        K = self.Kexp_g if abs(gam_k - g) < 1e-9 else self.Klin_2g
        return dress_mat(K, stack)


# =========================================================================== #
# reference bundle: active pair list, base density stacks, T over its nonzeros
# =========================================================================== #
class Ref:
    def __init__(self, act: np.ndarray, rho: np.ndarray, grad: np.ndarray,
                 d: np.ndarray, T_local, tag: str):
        self.act = act                      # global pair indices (for record)
        self.nact = len(act)
        self.rho = rho                      # (nact, NG)
        self.grad = grad                    # (nact, NG)
        self.rhov = rho * VNE[None, :]      # (nact, NG)  V_ne-weighted
        self.d = d                          # (nact,)  spectator integral of each pair density
        self.T_iA, self.T_iB, self.T_iC, self.T_iD, self.T_val = T_local   # LOCAL nnz indices
        self.tag = tag
        # cache dressed fields per (geminal name, kernel key, base name)
        self._cache: Dict[tuple, np.ndarray] = {}

    def base(self, name: str) -> np.ndarray:
        return {'rho': self.rho, 'grad': self.grad, 'rhov': self.rhov}[name]

    def dressed(self, gem: Geminal, key: str, bname: str) -> np.ndarray:
        ck = (gem.name, key, bname)
        if ck not in self._cache:
            self._cache[ck] = gem.dress(key, self.base(bname))
        return self._cache[ck]

    def dressed_prod(self, gem: Geminal, keys, bname: str) -> np.ndarray:
        ck = (gem.name, ('P',) + tuple(sorted(keys)), bname)
        if ck not in self._cache:
            self._cache[ck] = gem.dress_prod(keys, self.base(bname))
        return self._cache[ck]

    # ---- contract a set of already-computed factors with T over its nonzeros ---- #
    def contract(self, factors: List[Tuple[Tuple[int, ...], np.ndarray]]) -> float:
        idxmap = {1: self.T_iA, 2: self.T_iB, 3: self.T_iC, 4: self.T_iD}
        res = self.T_val.copy()
        for axes, tens in factors:
            sel = tuple(idxmap[e] for e in axes)
            res = res * tens[sel]
        return float(res.sum())


def _grid_einsum(ownw: np.ndarray, dress_list: List[Tuple[np.ndarray, int]]):
    """INT geo * own * prod(dressings)  -> (own_axis, dress_axes...) tensor.
    ownw = own * GI_W (weight folded in), shape (nown, NG).  dress_list = [(Psi (nd,NG), elec_id)]."""
    letters = "abcdefgh"
    subs_in = ["A" + "z"]  # own uses index 'A' for its pair axis, 'z' for grid... use explicit
    # build einsum: own 'Az', each dressing 'Bz','Cz',... , output 'ABC...'
    op_in = ["Az"]
    out = "A"
    arrs = [ownw]
    for k, (psi, _e) in enumerate(dress_list):
        lab = letters[k]
        op_in.append(lab + "z")
        out += lab
        arrs.append(psi)
    sub = ",".join(op_in) + "->" + out
    return np.einsum(sub, *arrs, optimize=True)


# =========================================================================== #
# generalised all-integrated f-diagram contractor (sigma^2, V_ne, drift, covariances)
# =========================================================================== #
def contract_diagram(edges: List[Tuple[int, int, str]], base_names: Dict[int, str],
                     gem: Geminal, ref: Ref) -> float:
    """<prod_edges kernel . one-body ops>  under |Psi|^2 (normalised).  Returns the T-contracted value.
    edges: (i, j, kernel_key).  base_names: electron -> {'rho','grad','rhov'} (default 'rho')."""
    bn = {e: base_names.get(e, 'rho') for e in (1, 2, 3, 4)}
    # per-electron state: own base stack (nact,NG) + accumulated dressings [(Psi,src_elec)]
    own = {e: ref.base(bn[e]) for e in (1, 2, 3, 4)}
    dressings: Dict[int, List[Tuple[np.ndarray, int]]] = {e: [] for e in (1, 2, 3, 4)}
    edge_list = [[set((i, j)), key] for (i, j, key) in edges]
    alive = {1, 2, 3, 4}
    factors: List[Tuple[Tuple[int, ...], np.ndarray]] = []

    def dnb(u):
        s = set()
        for (ed, _k) in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    while alive:
        # 1) integrate any electron with no live edges (spectator or a ready junction);
        # 2) else dress a PURE leaf (one neighbour, no accumulated dressings).  This ordering
        #    keeps a junction (deg-2 vertex) unreduced until BOTH its leaves have dressed onto it.
        u = next((c for c in sorted(alive) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(alive) if len(dnb(c)) == 1 and not dressings[c]), None)
        if u is None:
            raise RuntimeError(f"unexpected topology (no reducible vertex): {edges}")
        inc = [(ed, k) for (ed, k) in edge_list if u in ed]
        nb = dnb(u)
        if not nb:                                   # spectator / junction -> integrate
            ownw = own[u] * GI_W[None, :]
            tens = _grid_einsum(ownw, dressings[u])
            axes = (u,) + tuple(e for (_p, e) in dressings[u])
            factors.append((axes, tens))
        else:
            if dressings[u]:
                raise RuntimeError(f"leaf {u} carries dressings (chain>2): {edges}")
            w = next(iter(nb))
            keys = [k for (ed, k) in inc]
            if len(keys) == 1:
                Psi = ref.dressed(gem, keys[0], bn[u])
            elif len(keys) == 2:
                Psi = ref.dressed_prod(gem, keys, bn[u])
            else:
                raise RuntimeError(f"electron {u} has {len(keys)} edges to one neighbour")
            dressings[w].append((Psi, u))
        alive.discard(u)
        edge_list = [[ed, k] for (ed, k) in edge_list if u not in ed]

    return ref.contract(factors)


# =========================================================================== #
# no-operator + one-body-operator moments  (Fbar, sigma^2, V_ne, drift)
# =========================================================================== #
def moment_F(k: int, gem: Geminal, ref: Ref, op_elec: Optional[int] = None,
             op_base: str = 'rho', fkey: str = 'f') -> float:
    """<F^k O> where O is the identity (op_elec None) or a one-body op carried on op_elec (op_base).
    F = sum over the 6 pairs of f(r_ij).  k in {0,1,2}."""
    def bn():
        return {} if op_elec is None else {op_elec: op_base}
    if k == 0:
        return contract_diagram([], bn(), gem, ref)
    if k == 1:
        return sum(contract_diagram([(i, j, fkey)], bn(), gem, ref) for (i, j) in PAIRS6)
    tot = 0.0
    for (i, j) in PAIRS6:
        for (p, q) in PAIRS6:
            tot += contract_diagram([(i, j, fkey), (p, q, fkey)], bn(), gem, ref)
    return tot


def fbar(gem: Geminal, ref: Ref) -> float:
    return moment_F(1, gem, ref)


def sigma2(gem: Geminal, ref: Ref, Fbar: float) -> float:
    return moment_F(2, gem, ref) - Fbar ** 2


def hg_Vne(gem: Geminal, ref: Ref, Fbar: float):
    """h_Vne = <F V_ne> - Fbar <V_ne> ; g_Vne = <F^2 V_ne> - 2 Fbar <F V_ne> + Fbar^2 <V_ne>."""
    F0 = sum(moment_F(0, gem, ref, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    F1 = sum(moment_F(1, gem, ref, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    F2 = sum(moment_F(2, gem, ref, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    h_Vne = F1 - Fbar * F0
    g_Vne = F2 - 2 * Fbar * F1 + Fbar ** 2 * F0
    return h_Vne, g_Vne, F0


def drift_moments(gem: Geminal, ref: Ref):
    """S2 = <sum_i v_i^2>, FS = <F sum v^2>, F2S = <F^2 sum v^2>  (gradient on each electron)."""
    S2 = sum(moment_F(0, gem, ref, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    FS = sum(moment_F(1, gem, ref, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    F2S = sum(moment_F(2, gem, ref, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    return S2, FS, F2S


# ---- mixed 1-edge x kernel covariance (Cov[F,Y], Cov[F,E]) ---- #
def kbar(gem: Geminal, ref: Ref, kkey: str) -> float:
    """<sum_ij K(r_ij)> for a single kernel key (Yukawa 'Y' or exp 'E')."""
    return sum(contract_diagram([(i, j, kkey)], {}, gem, ref) for (i, j) in PAIRS6)


def cov_FK(gem: Geminal, ref: Ref, kkey: str, Fbar: float, Kbar: float) -> float:
    """Cov[F, K_sum] = <(sum f)(sum K)> - Fbar Kbar,  K in {Y, E}."""
    tot = 0.0
    for (i, j) in PAIRS6:
        for (p, q) in PAIRS6:
            tot += contract_diagram([(i, j, 'f'), (p, q, kkey)], {}, gem, ref)
    return tot - Fbar * Kbar


def hg_T(gem: Geminal, ref: Ref, Fbar: float, sig2: float):
    """h_T = PartA + PartB ; g_T = gT1 + gT23.  PartA/gT1 = drift; PartB/gT23 geminal-specific."""
    S2, FS, F2S = drift_moments(gem, ref)
    PartA = 0.5 * (FS - Fbar * S2)
    gT1 = 0.5 * (F2S - 2 * Fbar * FS + Fbar ** 2 * S2)
    g = gem.g
    Ybar = kbar(gem, ref, 'Y')
    CovFY = cov_FK(gem, ref, 'Y', Fbar, Ybar)
    if gem.name == "exp":
        PartB = -0.5 * g ** 2 * Fbar + g * Ybar
        gT23 = -g ** 2 * sig2 + 2 * g * CovFY
        Ebar = CovFE = None
    else:
        Ebar = kbar(gem, ref, 'E')
        CovFE = cov_FK(gem, ref, 'E', Fbar, Ebar)
        PartB = -0.5 * g ** 2 * Fbar + 2 * g * Ebar - Ybar
        gT23 = -g ** 2 * sig2 + 4 * g * CovFE - 2 * CovFY
    h_T = PartA + PartB
    g_T = gT1 + gT23
    return h_T, g_T, dict(PartA=PartA, PartB=PartB, gT1=gT1, gT23=gT23, S2=S2, FS=FS, F2S=F2S,
                          Ybar=Ybar, Ebar=Ebar, CovFY=CovFY, CovFE=CovFE)


# =========================================================================== #
# compound-field helpers  (data with trailing grid axis + a tuple of electron-id pair axes)
# =========================================================================== #
def cf_mul(A, B):
    """Elementwise over grid, outer over the distinct pair axes."""
    dA, axA = A; dB, axB = B
    NGl = dA.shape[-1]
    newA = dA.reshape(dA.shape[:-1] + (1,) * len(axB) + (NGl,))
    newB = dB.reshape((1,) * len(axA) + dB.shape)
    return (newA * newB, axA + axB)


def cf_dress(A, key: str, gem: Geminal):
    dA, axA = A
    flat = dA.reshape(-1, dA.shape[-1])
    return (gem.dress(key, flat).reshape(dA.shape), axA)


def cf_gint(A):
    """INT compound field d3r -> tensor over its pair axes."""
    dA, axA = A
    return (np.tensordot(dA, GI_W, axes=([-1], [0])), axA)


# =========================================================================== #
# V_ee reducer (kept Coulomb pair + RI-free triangle), generalised over T
# =========================================================================== #
def vee_reduce(rem_f: List[Tuple[int, int]], coul: Tuple[int, int], kernel_kind: str,
               gam_k: Optional[float], gem: Geminal, ref: Ref) -> float:
    """<(prod rem_f f-edges) . Coulomb(coul)>  under |Psi|^2, generalised over T.
    coul = (e, gg) kept pair (Coulomb/Yukawa on it).  Triangle when a non-kept vertex f-bonds BOTH."""
    e, gg = coul
    kept = {e, gg}
    # compound field per electron: (data, axes); starts as its own base density (axis = its id)
    field = {x: (ref.base('rho'), (x,)) for x in (1, 2, 3, 4)}
    edge_list = [set(ed) for ed in rem_f]
    alive = {1, 2, 3, 4}
    factors: List[Tuple[Tuple[int, ...], np.ndarray]] = []
    tri = None                                          # (data (nact,NG), axis) of the middle vertex

    def dnb(u):
        s = set()
        for ed in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    def ndressed(u):                                    # True if u has accumulated pair axes beyond its own
        return len(field[u][1]) > 1

    while True:
        nk = [u for u in alive if u not in kept]
        if not nk:
            break
        # 1) integrate a ready non-kept vertex; 2) dress any non-kept leaf; 3) else triangle.
        u = next((c for c in sorted(nk) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(nk) if len(dnb(c)) == 1), None)
        if u is None:                                   # triangle: vertex whose nbrs are BOTH kept
            u = next((c for c in sorted(nk) if dnb(c) == kept), None)
            if u is None:
                raise RuntimeError(f"vee: no reducible non-kept vertex; rem_f={rem_f} coul={coul}")
            if ndressed(u):
                raise RuntimeError("triangle vertex carries dressings")
            tri = field[u]
            edge_list = [ed for ed in edge_list if u not in ed]
            alive.discard(u)
            continue
        inc = [ed for ed in edge_list if u in ed]
        nb = dnb(u)
        if not nb:                                       # integrate (spectator or ready junction)
            data, axes = cf_gint(field[u])
            factors.append((axes, data))
        else:
            w = next(iter(nb)); mult = len(inc)
            key = 'f' if mult == 1 else 'f2'
            Psi = cf_dress(field[u], key, gem)           # dress the (possibly compound) field of u
            field[w] = cf_mul(field[w], Psi)
        alive.discard(u)
        edge_list = [ed for ed in edge_list if u not in ed]

    L = field[e]; Rt = field[gg]
    L_axes = L[1]; R_axes = Rt[1]
    Lflat = L[0].reshape(-1, NG); Rflat = Rt[0].reshape(-1, NG)

    if tri is not None:
        # SPARSE triangle: the kept pair carries no f-dressings (both its f-edges meet the middle
        # vertex), so L, R, mu are single-axis density stacks.  Evaluate the triangle ONLY at the
        # distinct (mu,L,R) index combos surviving in T's nonzeros (dense 55^3 is infeasible).
        # PHASE 3: read the cached, batched per-middle-vertex triangle matrix M_cm (get_tri_matrix)
        # instead of re-running the mode/rank/coul loop per triple -- exact, shared across triples.
        mu_data, mu_axis = tri
        assert len(L_axes) == 1 and len(R_axes) == 1, "triangle with kept-pair dressings unsupported"
        idxmap = {1: ref.T_iA, 2: ref.T_iB, 3: ref.T_iC, 4: ref.T_iD}
        w = ref.T_val.copy()
        for axes, tens in factors:                        # fold in the spectator scalar factors
            w = w * tens[tuple(idxmap[e] for e in axes)]
        im = idxmap[mu_axis[0]]; il = idxmap[L_axes[0]]; ir = idxmap[R_axes[0]]
        keep = np.abs(w) > 0
        im, il, ir, w = im[keep], il[keep], ir[keep], w[keep]
        if len(w) == 0:
            return 0.0
        total = 0.0
        for cm in np.unique(im):
            s = (im == cm)
            M = get_tri_matrix(gem, ref, int(cm))         # (nact, nact), cached + batched
            total += float((w[s] * M[il[s], ir[s]]).sum())
        return float(total)

    # SPARSE Coulomb/Yukawa: the kept-pair fields L, Rt can be compound (a kept electron dressed by
    # up to two f-edges -> up to nact^3 combos).  Building the full potential stack is infeasible at
    # 55 pairs, so -- as with the triangle -- fold the spectator factors into the nnz weight and
    # evaluate the potential ONLY at the distinct R-combos that survive in T's nonzeros.
    idxmap = {1: ref.T_iA, 2: ref.T_iB, 3: ref.T_iC, 4: ref.T_iD}
    w = ref.T_val.copy()
    for axes, tens in factors:
        w = w * tens[tuple(idxmap[e] for e in axes)]
    keep = np.abs(w) > 0
    if not np.any(keep):
        return 0.0
    w = w[keep]
    L_shape = L[0].shape[:-1]; R_shape = Rt[0].shape[:-1]
    Lidx = np.ravel_multi_index([idxmap[e][keep] for e in L_axes], L_shape) if len(L_axes) > 1 \
        else idxmap[L_axes[0]][keep]
    Ridx = np.ravel_multi_index([idxmap[e][keep] for e in R_axes], R_shape) if len(R_axes) > 1 \
        else idxmap[R_axes[0]][keep]
    dL, invL = np.unique(Lidx, return_inverse=True)
    dR, invR = np.unique(Ridx, return_inverse=True)
    RsubFlat = Rflat[dR]                                 # (ndR, NG)
    potR = neumann_stack(RsubFlat) if kernel_kind == 'coul' else gem.yuk_pot(RsubFlat, gam_k)
    Lsub = Lflat[dL] * GI_W[None, :]                     # (ndL, NG)
    Emat = Lsub @ potR.T                                 # (ndL, ndR) = INT L . pot(Rt)
    return float((w * Emat[invL, invR]).sum())


def _vee_triple(p, q, r, gem, ref):
    """<f_p f_q coul_r>  generalised."""
    rset = frozenset(r); fedges = [frozenset(p), frozenset(q)]
    n_on_r = sum(1 for ed in fedges if ed == rset)
    rem_f = [tuple(ed) for ed in fedges if ed != rset]
    kk = 'coul' if n_on_r == 0 else 'yuk'
    gam_k = None if n_on_r == 0 else (gem.g if n_on_r == 1 else 2 * gem.g)
    return vee_reduce(rem_f, r, kk, gam_k, gem, ref)


def _vee_pair_fc(p, r, gem, ref):
    """<f_p coul_r>  generalised."""
    rset = frozenset(r); pf = frozenset(p)
    if pf == rset:
        return vee_reduce([], r, 'yuk', gem.g, gem, ref)
    return vee_reduce([tuple(p)], r, 'coul', None, gem, ref)


def hg_Vee(gem: Geminal, ref: Ref, Fbar: float):
    """h_Vee = <F V_ee> - Fbar <V_ee> ; g_Vee = <F^2 V_ee> - 2 Fbar <F V_ee> + Fbar^2 <V_ee>."""
    Vee = sum(vee_reduce([], r, 'coul', None, gem, ref) for r in PAIRS6)
    FVee = sum(_vee_pair_fc(p, r, gem, ref) for p in PAIRS6 for r in PAIRS6)
    F2Vee = sum(_vee_triple(p, q, r, gem, ref) for p in PAIRS6 for q in PAIRS6 for r in PAIRS6)
    h_Vee = FVee - Fbar * Vee
    g_Vee = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    return h_Vee, g_Vee, Vee, FVee, F2Vee


# =========================================================================== #
# reference builders
# =========================================================================== #
def build_ref_2mo() -> Ref:
    """The 2-orbital closed-shell reference |phi_0^2 phi_1^2| on the general machinery.
    pairs = [(0,0),(0,1),(1,1)] ; rho = [P00,P01,P11] ; grad = [G00,G01,G11] ; d = [1,0,1]."""
    from lih_marriage_phase0 import fold_cs, pair_index, build_Eexp, build_T
    K = 2
    pairs, idx = pair_index(K)                      # [(0,0),(0,1),(1,1)]
    npair = len(pairs)
    apairs = list(combinations(range(K), 2))        # [(0,1)]
    n_ap = len(apairs)
    # single closed-shell determinant: c = 1 on the only det
    c = np.array([1.0])
    CS, _, _ = fold_cs(c, K)
    Eexp = build_Eexp(K, apairs, pairs, idx)
    Tfull, act_a, act_b = build_T(CS, Eexp, n_ap, npair)
    T4 = Tfull.reshape(npair, npair, npair, npair)
    # cross-check against the direct (1/4) DPmat (x) DPmat construction
    i00, i01, i11 = idx[(0, 0)], idx[(0, 1)], idx[(1, 1)]
    DP = np.zeros((npair, npair))
    DP[i00, i11] = 1.0; DP[i11, i00] = 1.0; DP[i01, i01] = -2.0
    T4_direct = 0.25 * np.einsum('AB,CD->ABCD', DP, DP)
    dev = np.max(np.abs(T4 - T4_direct))
    print(f"  [2-MO] build_T vs (1/4)DPmat(x)DPmat: max|dev| = {dev:.2e}")
    # grid densities (MO pair densities on the kernels.py grid)
    rho = np.stack([P00, P01, P11])                 # order matches pairs
    grad = np.stack([Gpq(0, 0), Gpq(0, 1), Gpq(1, 1)])
    d = np.array([1.0, 0.0, 1.0])
    # nnz of T4 -> local indices (already local: npair=3)
    iA, iB, iC, iD = np.nonzero(T4)
    val = T4[iA, iB, iC, iD]
    T_local = (iA, iB, iC, iD, val)
    print(f"  [2-MO] npair={npair}, T nnz={len(val)}, INT P = "
          f"{np.einsum('ABCD,A,B,C,D', T4, d, d, d, d):.12f}")
    return Ref(np.array([0, 1, 2]), rho, grad, d, T_local, "2-MO |phi0^2 phi1^2|")


def build_ref_artifact() -> Ref:
    """The 74-det CI natural-orbital reference from the Phase-0b artifact (active pairs only)."""
    z = np.load(REF_EXACT, allow_pickle=True)
    npair = z['pairs'].shape[0]
    row = z['T_row']; col = z['T_col']; data = z['T_data']
    A = (row // npair).astype(int); B = (row % npair).astype(int)
    C = (col // npair).astype(int); D = (col % npair).astype(int)
    act = np.array(sorted(set(A.tolist()) | set(B.tolist()) | set(C.tolist()) | set(D.tolist())))
    loc = {g: i for i, g in enumerate(act)}
    iA = np.array([loc[x] for x in A]); iB = np.array([loc[x] for x in B])
    iC = np.array([loc[x] for x in C]); iD = np.array([loc[x] for x in D])
    rho = np.asarray(z['rho_no'])[act]                              # (nact, NG)
    ng = np.asarray(z['no_grad'])                                   # (Mk, NG, 3)
    pairs = [tuple(int(x) for x in pr) for pr in z['pairs']]
    grad = np.empty((len(act), NG))
    for k, gpair in enumerate(act):
        p, q = pairs[gpair]
        grad[k] = np.einsum('gd,gd->g', ng[p], ng[q])
    d = np.asarray(z['delta_pair'])[act]
    T_local = (iA, iB, iC, iD, data)
    print(f"  [74-det] npair active={len(act)} of {npair}, T nnz={len(data)}, INT P = "
          f"{float((data * d[iA] * d[iB] * d[iC] * d[iD]).sum()):.9f}")
    return Ref(act, rho, grad, d, T_local, "74-det CI natural-orbital")


# =========================================================================== #
# assemble the 2x2 -> E_R12
# =========================================================================== #
def assemble(pieces: dict, E0_tot: float) -> Tuple[float, float]:
    sig2 = pieces['sigma2']
    h = pieces['h_T'] + pieces['h_Vne'] + pieces['h_Vee']
    g_elec = pieces['g_T'] + pieces['g_Vne'] + pieces['g_Vee']
    E0_elec = E0_tot - V_NN
    aa = E0_elec; cc = g_elec / sig2; bb = h / np.sqrt(sig2)
    E_R12 = 0.5 * (aa + cc - np.sqrt((aa - cc) ** 2 + 4 * bb ** 2)) + V_NN
    return E_R12, h


def all_pieces(gem: Geminal, ref: Ref, T0: float, verbose: bool = True) -> dict:
    p = {}
    Fbar = fbar(gem, ref); p['Fbar'] = Fbar
    sig2 = sigma2(gem, ref, Fbar); p['sigma2'] = sig2
    if verbose:
        print(f"    Fbar={Fbar:.6f}  sigma2={sig2:.6f}  {_t(T0)}", flush=True)
    h_Vne, g_Vne, Vne_exp = hg_Vne(gem, ref, Fbar)
    p['h_Vne'] = h_Vne; p['g_Vne'] = g_Vne; p['Vne_exp'] = Vne_exp
    if verbose:
        print(f"    h_Vne={h_Vne:+.6f}  g_Vne={g_Vne:+.6f}  <V_ne>={Vne_exp:+.6f}  {_t(T0)}", flush=True)
    h_T, g_T, tinfo = hg_T(gem, ref, Fbar, sig2)
    p['h_T'] = h_T; p['g_T'] = g_T; p.update({f"T_{k}": v for k, v in tinfo.items()})
    if verbose:
        print(f"    h_T={h_T:+.6f} (PartA {tinfo['PartA']:+.5f} PartB {tinfo['PartB']:+.5f})  "
              f"g_T={g_T:+.6f} (gT1 {tinfo['gT1']:+.5f} gT23 {tinfo['gT23']:+.5f})  {_t(T0)}", flush=True)
    h_Vee, g_Vee, Vee, FVee, F2Vee = hg_Vee(gem, ref, Fbar)
    p['h_Vee'] = h_Vee; p['g_Vee'] = g_Vee; p['Vee'] = Vee; p['FVee'] = FVee; p['F2Vee'] = F2Vee
    if verbose:
        print(f"    h_Vee={h_Vee:+.6f}  g_Vee={g_Vee:+.6f}  <V_ee>={Vee:+.6f}  {_t(T0)}", flush=True)
    # E0 grid-consistent: <T> = S2/2, <V_ne>, <V_ee>
    T_kin = 0.5 * p['T_S2']
    E0_grid = T_kin + Vne_exp + Vee + V_NN
    p['T_kin'] = T_kin; p['E0_grid'] = E0_grid
    if verbose:
        print(f"    <T>={T_kin:+.6f}  E0_grid = <T>+<V_ne>+<V_ee>+V_NN = {E0_grid:+.6f}  {_t(T0)}",
              flush=True)
    return p


# =========================================================================== #
def run_gate(paths: Sequence[str], gems: Sequence[str], preview: bool):
    from geovac.lih_r12ci.assembly import ANALYTIC_REF
    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    print("=" * 100)
    print("LiH MARRIAGE -- PHASE 2: seven enumerator generalisations over T + GATE G1 (2-MO regression)")
    print(f"  R={R_ENG}  V_NN={V_NN:.6f}  grid {NXI}x{NETA} (NG={NG})  GAM={GAM}")
    print("=" * 100, flush=True)

    EXACT_REF = {'linexp': -7.917199, 'exp': -7.943605}   # Phase 0b G-e exact-kernel references

    results = {}
    ref2 = None
    for path in paths:
        KG.USE_EXACT_NEUMANN = (path == 'exact')
        # clear per-path Neumann caches AND the Ref dressing cache (Yukawa dressings are
        # flag-dependent; the Ref persists across paths, so its cache MUST be invalidated).
        GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); _TRI_MAT_CACHE.clear()
        if ref2 is not None:
            ref2._cache.clear()
        print(f"\n{'#' * 96}\n# KERNEL PATH: {path}  (USE_EXACT_NEUMANN={KG.USE_EXACT_NEUMANN})\n{'#' * 96}", flush=True)
        for gname in gems:
            print(f"\n--- geminal = {gname} | path = {path} ---", flush=True)
            # LIVE oracle: assembly.energy() on THIS kernel path (path-consistent per-piece truth).
            from geovac.lih_r12ci.assembly import energy as assembly_energy
            oracle = assembly_energy(gname)
            op = oracle.pieces
            print(f"  [oracle] assembly.energy({gname}) on path={path}: E_R12={oracle.E_R12:.6f} "
                  f"sigma2={oracle.sigma2:.6f} {_t(T0)}", flush=True)
            gem = Geminal(gname)
            if ref2 is None:
                ref2 = build_ref_2mo()
            p = all_pieces(gem, ref2, T0)
            ref = ANALYTIC_REF[gname]
            E_R12, hval = assemble(p, -7.887822)     # fixed-E0 (assembly convention) for direct compare
            E_R12_gridE0, _ = assemble(p, p['E0_grid'])
            target_E = EXACT_REF[gname] if path == 'exact' else ref['E_R12']
            # compare each piece to the LIVE oracle (primary) and ANALYTIC_REF (legacy hardcoded)
            rows = [
                ('sigma2', p['sigma2'], op['sigma2'], ref['sigma2']),
                ('h_Vne', p['h_Vne'], op['h_Vne'], ref['h_Vne']),
                ('h_T', p['h_T'], op['h_T'], ref['h_T']),
                ('h_Vee', p['h_Vee'], op['h_Vee'], ref['h_Vee']),
                ('g_Vne', p['g_Vne'], op['g_Vne'], ref['g_Vne']),
                ('g_T', p['g_T'], op['g_T'], ref['g_T']),
                ('g_Vee', p['g_Vee'], op['g_Vee'], ref['g_Vee']),
            ]
            print(f"\n  G1 piece table (path={path}, geminal={gname}) -- generalized vs LIVE oracle:")
            print(f"    {'piece':8s} {'generalized':>13s} {'oracle(live)':>13s} {'dOracle(mHa)':>13s} "
                  f"{'ANALYTIC_REF':>13s} {'dRef(mHa)':>11s}")
            worst = 0.0; worst_name = ''
            for nm, got, olive, rf in rows:
                do = (got - olive) * 1e3
                dr = (got - rf) * 1e3
                if abs(do) > worst:
                    worst = abs(do); worst_name = nm
                print(f"    {nm:8s} {got:13.6f} {olive:13.6f} {do:+13.4f} {rf:13.6f} {dr:+11.4f}")
            print(f"    {'E0_grid':8s} {p['E0_grid']:13.6f} (vs hardcoded -7.887822: "
                  f"{(p['E0_grid'] + 7.887822) * 1e3:+.4f} mHa)")
            dev_o = abs(E_R12 - oracle.E_R12) * 1e3
            dev_t = abs(E_R12 - target_E) * 1e3
            print(f"    E_R12 (fixed-E0) = {E_R12:.6f}  vs oracle {oracle.E_R12:.6f} (d {dev_o:+.4f} mHa)"
                  f"  vs path-target {target_E:.6f} (d {dev_t:+.4f} mHa)")
            print(f"    E_R12 (grid-E0)  = {E_R12_gridE0:.6f}")
            results[(path, gname)] = dict(pieces=p, E_R12=E_R12, E_R12_gridE0=E_R12_gridE0,
                                          oracle_ER12=oracle.E_R12, oracle_pieces=dict(op),
                                          worst=worst, worst_name=worst_name)
            print(f"  worst piece dev vs LIVE oracle: {worst:.4f} mHa ({worst_name})", flush=True)

    # ---- preview: 74-det tensor through the enumerators on the exact path (UNGATED) ----
    preview_res = {}
    if preview:
        KG.USE_EXACT_NEUMANN = True
        GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); _TRI_MAT_CACHE.clear()
        print(f"\n{'#' * 96}\n# 74-DET PREVIEW (exact path, UNGATED -- awaits Phase 3 G2 VMC cross-check)\n{'#' * 96}", flush=True)
        refA = build_ref_artifact()
        for gname in gems:
            print(f"\n--- PREVIEW geminal = {gname} ---", flush=True)
            gem = Geminal(gname)
            p = all_pieces(gem, refA, T0)
            E_R12, _ = assemble(p, p['E0_grid'])
            inwin = -8.055 <= E_R12 <= -8.045
            print(f"  >>> PREVIEW E_R12 ({gname}) = {E_R12:.6f} Ha  "
                  f"[{'in' if inwin else 'OUTSIDE'} the -8.045..-8.055 window]", flush=True)
            preview_res[gname] = dict(pieces=p, E_R12=E_R12, in_window=inwin)

    # ---- save ----
    save = dict(V_NN=V_NN, R=R_ENG, EXACT_REF=str(EXACT_REF))
    for (path, gname), r in results.items():
        pre = f"{path}_{gname}_"
        save[pre + 'E_R12'] = r['E_R12']
        save[pre + 'E_R12_gridE0'] = r['E_R12_gridE0']
        save[pre + 'oracle_ER12'] = r['oracle_ER12']
        save[pre + 'worst_mHa'] = r['worst']
        for k, v in r['pieces'].items():
            if isinstance(v, (int, float)) or (v is None):
                save[pre + k] = np.nan if v is None else v
        for k, v in r['oracle_pieces'].items():
            if isinstance(v, (int, float)) or (v is None):
                save[pre + 'oracle_' + k] = np.nan if v is None else v
    for gname, r in preview_res.items():
        save[f"preview_{gname}_E_R12"] = r['E_R12']
        for k, v in r['pieces'].items():
            if isinstance(v, (int, float)) or (v is None):
                save[f"preview_{gname}_" + k] = np.nan if v is None else v
    os.makedirs(DATA, exist_ok=True)
    np.savez(OUT, **save)
    print(f"\nsaved {OUT}")
    print(f"wall time {time.time() - T0:.0f} s")
    return results, preview_res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--preview', action='store_true')
    ap.add_argument('--no-exp', action='store_true')
    ap.add_argument('--no-linexp', action='store_true')
    ap.add_argument('--paths', default='legacy,exact')
    args = ap.parse_args()
    gems = [g for g in ('exp', 'linexp') if not getattr(args, f'no_{g}')]
    paths = [p.strip() for p in args.paths.split(',') if p.strip()]
    run_gate(paths, gems, args.preview)
