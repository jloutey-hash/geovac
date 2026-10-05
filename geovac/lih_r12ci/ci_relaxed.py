r"""CI coefficients relaxed alongside the geminal -- analytic ansatz (A), RI-free and deterministic.

Backs (Paper 19): the LiH explicit-r12 "marriage" with the CI coefficients of the truncated M=16 pi
natural-orbital reference RELAXED together with the geminal (ansatz (A) of the CI-relaxation step-2
diagnostic, debug/sprint_ci_relax_step2_memo.md, whose VMC linear-method numbers it reproduces
analytically).  The frozen-coefficient marriage (pi_m16 / gvee_sparse, the {Psi_CI, (F-Fbar) Psi_CI}
2x2; -8.048141 nine-det, -8.049045 ninety-det, e^{-1.5 r}) is the 2-dimensional restriction of the
problem solved here.

PROBLEM.  Normalised determinants chi_I = dA[a_I] dB[b_I] / 2 (unsigned grouped minors, exactly the
LM basis), reference psi = sum_I d_I chi_I with d_I = CS[a_I, b_I] (sum d^2 = 1), xi = (F - Fbar) psi.
The (N+1)-dimensional generalized eigenproblem in {chi_I} U {xi}:

    S = [[ 1_N , s ], [ s^T, sigma^2 ]],      H = [[ H_CI, t ], [ t^T, g ]]   (electronic; + V_NN)
    s_I = <chi_I|F|psi> - Fbar d_I,           t_I = <chi_I|H|xi>,
    H_CI[I,J] = <chi_I|H|chi_J>  on the SAME grid quadrature as the 2x2 (E0_grid),
    g, sigma^2 = the frozen 2x2's pieces (the expensive O(nnz) g_Vee is NOT recomputed).

Every new element is a TRANSITION contraction: the pair-index coefficient tensor of pi_m16 built with
two different coefficient matrices (bra = one determinant, ket = psi or another determinant).  It is
NOT bra/ket symmetric, so ordered pair indices (x = bra orbital, y = ket orbital) are kept as such.
All one-electron factor tensors (dressings, junctions, kept-pair potentials) are T-independent and are
built once (on the full set of Mk^2 ordered pairs) and contracted against every transition tensor.

KINETIC TRANSITION (the genuinely new piece).  With u = chi_I, v = psi and sums over electrons,
    <u|T|F v> = 1/2 INT F grad u* . grad v  -  1/4 INT (lap F) u* v  +  1/4 INT F (u* lap v - v lap u*).
The first two terms are the existing h_T formula (drift 'grad' base + the geminal PartB kernel) with
the transition tensor; the third is the antisymmetric current term, IBP'd to a SCALAR ordered-pair
base  LA_(x,y) = phi_x* lap phi_y - (lap phi_x)* phi_y  (orbital Laplacians on the grid).  It vanishes
when contracted with the reference (bra = ket, LA antisymmetric), so sum_I d_I t_I == h identically and
the frozen restriction reproduces the 2x2 to roundoff (gate G1).

API:
    ci_relaxed_energy(gem, refm, info, pieces) -> dict  (E_A, E_frozen, s, t, H_CI, vector, ...)
    columns(gem, info) -> transition columns only
Side effects: none beyond pi_m16 / gvee_sparse (reads KG.GAM at GeminalM construction; set the
geminal range with geovac.lih_r12ci.configure_geminal BEFORE building the GeminalM).
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg as sl

import geovac.lih_r12ci.enumerators as P2
import geovac.lih_r12ci.pi_m16 as M
import geovac.lih_r12ci.gvee_sparse as GS
from geovac.lih_r12ci.transition_tensor import (transition_stack, split_owners, determinants,
                                                bra_det_ket_psi, bra_det_ket_det)

NG = M.NG
GI_W = P2.GI_W
PAIRS6 = P2.PAIRS6
V_NN = M.V_NN


# =========================================================================== #
# full ordered-pair reference (bases rho, grad, rhov + the new antisymmetric Laplacian base lapA)
# =========================================================================== #
class FullPairRef(M.RefM):
    """RefM over ALL Mk^2 ordered pairs (x,y), index x*Mk+y; T is supplied per contraction."""

    def __init__(self, info: dict, tag: str):
        Mk = int(info['Mk'])
        nval = info['nval']; ngrad = info['ngrad']; nlap = info['nlap']; m_no = info['m_no']
        pairs = [(x, y) for x in range(Mk) for y in range(Mk)]
        rho = np.empty((len(pairs), NG)); grad = np.empty((len(pairs), NG))
        lapA = np.empty((len(pairs), NG))
        for k, (x, y) in enumerate(pairs):
            rho[k] = (nval[x] * nval[y]).real                       # == pi_m16.build_ref_from_prims
            grad[k] = np.einsum('gd,gd->g', np.conj(ngrad[x]), ngrad[y]).real
            lapA[k] = (nval[x] * nlap[y]).real - (nlap[x] * nval[y]).real
        d = np.array([1.0 if x == y else 0.0 for (x, y) in pairs])
        mA = np.array([m_no[y] - m_no[x] for (x, y) in pairs], dtype=int)
        empty = (np.zeros(0, int),) * 4 + (np.zeros(0),)
        super().__init__(np.arange(len(pairs)), rho, grad, d, mA, empty, tag)
        self.lapA = lapA
        self.Mk = Mk
        self.pairs = pairs

    def base(self, name: str) -> np.ndarray:
        if name == 'lapA':
            return self.lapA
        return super().base(name)

    def view(self, T_local) -> "FullPairRef":
        """Shallow copy carrying a different coefficient tensor; shares fields AND the dressing cache."""
        v = object.__new__(FullPairRef)
        v.__dict__.update(self.__dict__)
        v.T_iA, v.T_iB, v.T_iC, v.T_iD, v.T_val = T_local
        return v


# =========================================================================== #
# T-independent diagram factors (the loop of pi_m16.contract_diagram_m, returning the factors)
# =========================================================================== #
def diagram_factors(edges, base_names: Dict[int, str], gem: "M.GeminalM", ref: FullPairRef):
    """Factor list of pi_m16.contract_diagram_m(edges, base_names, gem, ref) -- same reduction order,
    same mode-aware dressings, same L_z masks -- without contracting against a tensor."""
    bn = {e: base_names.get(e, 'rho') for e in (1, 2, 3, 4)}
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
        u = next((c for c in sorted(alive) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(alive) if len(dnb(c)) == 1 and not dressings[c]), None)
        if u is None:
            raise RuntimeError(f"unexpected topology: {edges}")
        inc = [(ed, k) for (ed, k) in edge_list if u in ed]
        nb = dnb(u)
        if not nb:
            ownw = own[u] * GI_W[None, :]
            tens = P2._grid_einsum(ownw, dressings[u])
            axes = (u,) + tuple(e for (_p, e) in dressings[u])
            tens = tens * M._circ_mask(ref.mA, len(axes))
            factors.append((axes, tens))
        else:
            if dressings[u]:
                raise RuntimeError(f"leaf {u} carries dressings (chain>2): {edges}")
            w = next(iter(nb))
            keys = [k for (ed, k) in inc]
            keys = keys[0] if len(keys) == 1 else tuple(keys)
            Psi = M._refm_dressed(ref, gem, keys, bn[u])
            dressings[w].append((Psi, u))
        alive.discard(u)
        edge_list = [[ed, k] for (ed, k) in edge_list if u not in ed]
    return factors


def contract_owners(factors, stack, n_owner: int) -> np.ndarray:
    """sum over each owner's nonzeros of val * prod factor[...]  (owner-resolved _refm_contract)."""
    iA, iB, iC, iD, val, ow = stack
    idxmap = {1: iA, 2: iB, 3: iC, 4: iD}
    res = val.copy()
    for axes, tens in factors:
        res = res * tens[tuple(idxmap[e] for e in axes)]
    return np.bincount(ow, weights=res, minlength=n_owner)


class _Diagrams:
    """factor cache keyed by the diagram, so each diagram is reduced once for all owners."""

    def __init__(self, gem, ref: FullPairRef):
        self.gem = gem; self.ref = ref
        self._f: Dict[tuple, list] = {}

    def value(self, edges, base_names, stack, n_owner) -> np.ndarray:
        key = (tuple(edges), tuple(sorted(base_names.items())))
        if key not in self._f:
            self._f[key] = diagram_factors(edges, base_names, self.gem, self.ref)
        return contract_owners(self._f[key], stack, n_owner)

    def moment(self, k: int, stack, n_owner, op_elec: Optional[int] = None, op_base: str = 'rho',
               fkey: str = 'f') -> np.ndarray:
        """owner-resolved pi_m16.moment_F(k in {0,1}, ...)."""
        bn = {} if op_elec is None else {op_elec: op_base}
        if k == 0:
            return self.value([], bn, stack, n_owner)
        if k == 1:
            return sum(self.value([(i, j, fkey)], bn, stack, n_owner) for (i, j) in PAIRS6)
        raise ValueError(k)

    def onebody(self, k: int, base: str, stack, n_owner) -> np.ndarray:
        return sum(self.moment(k, stack, n_owner, op_elec=i, op_base=base) for i in (1, 2, 3, 4))


def _vee_owner(gem, ref: FullPairRef, parts: List[tuple], one_edge: bool, mem_cap_gb: float = 6.0,
               ctx: Optional["GS.Ctx"] = None) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Per-owner <V_ee> and (one_edge) <F V_ee>, through gvee_sparse.vee_reduce_sparse (the O(nnz)
    twin of pi_m16.vee_reduce_m) with one shared memo context; no triangle can occur here."""
    if ctx is None:
        ctx = GS.Ctx(gem, ref, mem_cap_gb)
    col = GS.TriCollector()
    Vee = np.zeros(len(parts)); FVee = np.zeros(len(parts)) if one_edge else None
    for o, T in enumerate(parts):
        ctx.refm = ref.view(T)
        Vee[o] = sum(GS.vee_reduce_sparse([], r, 'coul', None, ctx, col) for r in PAIRS6)
        if one_edge:
            tot = 0.0
            for p in PAIRS6:
                for r in PAIRS6:
                    if frozenset(p) == frozenset(r):
                        tot += GS.vee_reduce_sparse([], r, 'yuk', gem.g, ctx, col)
                    else:
                        tot += GS.vee_reduce_sparse([tuple(p)], r, 'coul', None, ctx, col)
            FVee[o] = tot
    if col.recs:
        raise RuntimeError("unexpected triangle in a one-edge V_ee transition")
    ctx.refm = ref
    return Vee, FVee


# =========================================================================== #
# the transition columns and the grid CI matrix
# =========================================================================== #
def _partB(gem, SF: np.ndarray, Y: np.ndarray, E: Optional[np.ndarray]) -> np.ndarray:
    """-1/4 <lap F> with the transition tensor (pi_m16.hg_T's PartB, Fbar/Ybar/Ebar -> per owner)."""
    g = gem.g
    if gem.name == 'exp':
        return -0.5 * g ** 2 * SF + g * Y
    return -0.5 * g ** 2 * SF + 2 * g * E - Y


def columns(gem: "M.GeminalM", info: dict, Fbar: float, mem_cap_gb: float = 6.0,
            verbose: bool = False, ref: Optional[FullPairRef] = None) -> dict:
    """Transition columns <chi_I|F|psi>, <chi_I|H(F-Fbar)|psi> (by part) and the grid CI matrix."""
    t0 = time.time()
    ref = ref if ref is not None else FullPairRef(info, "M=16 pi full ordered pairs")
    ia, ib, d, apairs = determinants(info['ct'], int(info['Mk']))
    N = len(d); Mk = ref.Mk
    D = _Diagrams(gem, ref)
    stI = transition_stack(*bra_det_ket_psi(ia, ib, d), apairs, Mk)      # bra chi_I, ket psi: owner I
    SF = D.moment(1, stI, N)                                       # <chi_I|F|psi>
    S0 = D.moment(0, stI, N)                                       # <chi_I|psi> (grid)
    Y = sum(D.value([(i, j, 'Y')], {}, stI, N) for (i, j) in PAIRS6)
    E = sum(D.value([(i, j, 'E')], {}, stI, N) for (i, j) in PAIRS6) if gem.name == 'linexp' else None
    S2 = D.onebody(0, 'grad', stI, N); FS = D.onebody(1, 'grad', stI, N)
    V0 = D.onebody(0, 'rhov', stI, N); V1 = D.onebody(1, 'rhov', stI, N)
    LA1 = D.onebody(1, 'lapA', stI, N); LA0 = D.onebody(0, 'lapA', stI, N)
    if verbose:
        print(f"    [ci_relaxed] one-body / F columns done  [{time.time()-t0:6.1f}s]", flush=True)
    ctx = GS.Ctx(gem, ref, mem_cap_gb)
    Vee, FVee = _vee_owner(gem, ref, split_owners(stI, N), True, ctx=ctx)
    if verbose:
        print(f"    [ci_relaxed] V_ee columns done  [{time.time()-t0:6.1f}s]", flush=True)
    PartB = _partB(gem, SF, Y, E)
    KF = 0.5 * FS + PartB + 0.25 * LA1                            # <chi_I| T F |psi>
    tT = KF - Fbar * 0.5 * S2                                     # kinetic part of t_I
    tVne = V1 - Fbar * V0
    tVee = FVee - Fbar * Vee
    # ---- grid CI matrix: owners (I,J) ----
    stIJ = transition_stack(*bra_det_ket_det(ia, ib), apairs, Mk)        # bra chi_I, ket chi_J
    HT = 0.5 * D.onebody(0, 'grad', stIJ, N * N)
    HV = D.onebody(0, 'rhov', stIJ, N * N)
    SIJ = D.moment(0, stIJ, N * N)
    HVee, _ = _vee_owner(gem, ref, split_owners(stIJ, N * N), False, ctx=ctx)
    H_CI = (HT + HV + HVee).reshape(N, N)
    if verbose:
        print(f"    [ci_relaxed] grid CI matrix done (N={N})  [{time.time()-t0:6.1f}s]", flush=True)
    return dict(N=N, ia=ia, ib=ib, d=d, SF=SF, S0=S0, Y=Y, E=E, S2=S2, FS=FS, V0=V0, V1=V1,
                LA1=LA1, LA0=LA0, Vee=Vee, FVee=FVee, PartB=PartB, KF=KF,
                tT=tT, tVne=tVne, tVee=tVee, t=tT + tVne + tVee, s=SF - Fbar * d,
                H_CI=H_CI, H_CI_T=HT.reshape(N, N), H_CI_Vne=HV.reshape(N, N),
                H_CI_Vee=HVee.reshape(N, N), S_CI_grid=SIJ.reshape(N, N), wall_s=time.time() - t0)


def lowest_root(H: np.ndarray, S: np.ndarray) -> Tuple[float, np.ndarray]:
    """Lowest generalized eigenpair of the real symmetric pencil (H, S), S positive definite."""
    w, V = sl.eigh(0.5 * (H + H.T), 0.5 * (S + S.T))
    return float(w[0]), V[:, 0]


def ci_relaxed_energy(gem: "M.GeminalM", refm: "M.RefM", info: dict, pieces: dict,
                      mem_cap_gb: float = 6.0, verbose: bool = False) -> dict:
    """Ansatz (A): the {chi_I} U {(F - Fbar) psi} generalized eigenproblem.  `pieces` = the frozen 2x2
    pieces of this reference/geminal (gvee_sparse.all_pieces_sparse output, or the golden oracle):
    uses Fbar, sigma2, g_T+g_Vne+g_Vee, and (identity checks only) h_*, E0_grid.  `refm` is only used
    for the bookkeeping check that the full-pair tensor reproduces it.  Returns energies (total, Ha)
    and the identity diagnostics of gate G1."""
    Fbar = float(pieces['Fbar']); sig2 = float(pieces['sigma2'])
    g = float(pieces['g_T'] + pieces['g_Vne'] + pieces['g_Vee'])
    h_ref = float(pieces['h_T'] + pieces['h_Vne'] + pieces['h_Vee'])
    col = columns(gem, info, Fbar, mem_cap_gb=mem_cap_gb, verbose=verbose)
    N = col['N']; d = col['d']; s = col['s']; t = col['t']; H_CI = col['H_CI']
    # ---- frozen restriction (the existing 2x2, from the new code) ----
    aa = float(d @ H_CI @ d); h_new = float(d @ t)
    E_frozen, c_frozen, _ = M.assemble_c(dict(sigma2=sig2, h_T=h_new, h_Vne=0.0, h_Vee=0.0,
                                              g_T=g, g_Vne=0.0, g_Vee=0.0), aa + V_NN)
    E_frozen_ref, _, _ = M.assemble_c(pieces, pieces['E0_grid'])
    # ---- (A): N+1 ----
    S = np.eye(N + 1); S[:N, N] = s; S[N, :N] = s; S[N, N] = sig2
    H = np.zeros((N + 1, N + 1)); H[:N, :N] = H_CI; H[:N, N] = t; H[N, :N] = t; H[N, N] = g
    eA, vA = lowest_root(H, S)
    E_A = eA + V_NN
    eC, vC = lowest_root(H_CI, np.eye(N))
    if float(vA[:N] @ d) < 0:                                     # phase: align with the reference
        vA = -vA
    dA = vA[:N] / np.linalg.norm(vA[:N])
    ident = dict(
        aa_vs_E0grid=aa + V_NN - float(pieces['E0_grid']),
        h_vs_pieces=h_new - h_ref,
        hT_vs=float(d @ col['tT']) - float(pieces['h_T']),
        hVne_vs=float(d @ col['tVne']) - float(pieces['h_Vne']),
        hVee_vs=float(d @ col['tVee']) - float(pieces['h_Vee']),
        Fbar_vs=float(d @ col['SF']) - Fbar,
        s_dot_d=float(d @ s),
        S2_vs=float(d @ col['S2']) - float(pieces['T_S2']),
        FS_vs=float(d @ col['FS']) - float(pieces['T_FS']),
        Ybar_vs=float(d @ col['Y']) - float(pieces['T_Ybar']),
        LA_current_ref=float(d @ col['LA1']),
        Vee_vs=float(d @ col['Vee']) - float(pieces['Vee']),
        FVee_vs=float(d @ col['FVee']) - float(pieces['FVee']),
        E_frozen_vs_2x2=E_frozen - E_frozen_ref,
    )
    out = dict(E_A=E_A, E_frozen=E_frozen, E_frozen_2x2=E_frozen_ref, c_frozen=c_frozen,
               E_CI_relaxed_only=eC + V_NN, E0_grid=aa + V_NN, vec=vA, d_relaxed=dA, d_ref=d,
               overlap=float(abs(dA @ d)), c_geminal=float(vA[N] / (vA[:N] @ d)),
               S=S, H=H, identity=ident, N=N)
    out.update({k: col[k] for k in ('s', 't', 'tT', 'tVne', 'tVee', 'SF', 'KF', 'H_CI', 'S_CI_grid',
                                    'V1', 'FVee', 'S0', 'wall_s')})
    return out
